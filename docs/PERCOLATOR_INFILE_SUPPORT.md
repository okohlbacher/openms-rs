# PercolatorInfile support

`FORMAT/PercolatorInfile.h` and `FORMAT/PercolatorInfile.cpp` at `bc9cc12`,
ported to `src/format/percolator_infile.rs` and tested by
`tests/percolator_infile.rs`. Fixtures, hashes and source anchors are in
`tests/data/percolator_infile_provenance.json`.

A `.pin` file is one tab-separated header line naming the columns, then one
line per peptide-spectrum match. The contract Percolator parses is: three
mandatory leading columns `SpecId`, `Label`, `ScanNr`; then the per-PSM feature
columns; then `Peptide` and `Proteins`, in that order, last. The module has no
Cargo feature gate, because it needs nothing beyond the crate's identification
records and its `TextFile`/`CsvFile` ports.

## API mapping

Every public and protected member of the header.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `PercolatorInfile()` (implicit) | `PinOptions::default()`, `ReadOptions::default()` | The class is a bag of static functions; its constructor and destructor are only checked for existence. There is no object here. |
| `~PercolatorInfile()` (implicit) | not ported: nothing to release. |
| `using PinFeatureMetaValueMap = std::map<std::pair<size_t, size_t>, std::set<std::string>>` | `StampReport::added_meta_values` (`BTreeMap<(usize, usize), BTreeSet<String>>`) | Same shape; folded into the report rather than an out-parameter. |
| `static void store(const std::string&, const PeptideIdentificationList&, const StringList&, const std::string&, int, int)` | `store` | Returns the `PreparedPin` instead of discarding it, so the caller sees what the source only logged. `enz`, `min_charge` and `max_charge` move into `PinOptions`. |
| `static PeptideIdentificationList load(const std::string&, bool, const std::string&, const StringList&, StringList&, std::string, double, bool)` | `load`, `read` | The `filenames` out-parameter becomes `PinDocument::filenames`. `higher_score_better`, `score_name`, `extra_scores`, `decoy_prefix`, `threshold` and `SageAnnotation` become `ReadOptions` fields (`spectrum_q_threshold`, `sage_annotation`). `read` takes a stream and refuses `sage_annotation`, which needs a path. |
| `static std::string getScanIdentifier(const PeptideIdentification&, size_t)` | `scan_identifier` | |
| `static StringList getStandardFeatureSet(int, int)` | `standard_feature_set` | Returns `Result`, because the charge span is now bounded. |
| `static std::set<std::pair<size_t, size_t>> stampPinFeaturesOnHits(PeptideIdentificationList&, const std::string&, int, int)` | `stamp_pin_features` | The skipped set is `StampReport::skipped`. |
| `static std::set<std::pair<size_t, size_t>> stampPinFeaturesOnHits(..., PinFeatureMetaValueMap&)` | `stamp_pin_features` | The two overloads collapse into one function: the added-keys record is always produced, since the source's four-argument form merely discards it. |
| `static TextFile preparePin_(const PeptideIdentificationList&, const StringList&, const std::string&, int, int)` | `prepare_pin` | Protected in the source, public here: it is the only way to get the `.pin` text without writing a file, which a tool that streams its output needs. Returns `PreparedPin`, carrying the lines plus the drop counts the source only logged. |
| `static bool isEnz_(const char&, const char&, const std::string&)` | `is_enzymatic` | Public here: the enzyme table is a specification a caller may want to query, and the port's own test needs it. |
| `static Size countEnzymatic_(const std::string&, const std::string&)` | `count_enzymatic` | Public for the same reason. |
| (none) | `MANDATORY_COLUMNS`, `MASS_COLUMNS`, `ENZYME_COLUMNS`, `TRAILING_COLUMNS` | The column contract as named constants; the source spells the same strings inline in `getStandardFeatureSet` and in every adapter. |
| (none) | `ScanNumberPattern` | `SpectrumNativeIDParser::getRegExFromNativeID` plus `extractScanNumber`, reduced to the literal prefix each pattern looks for. The source builds a `boost::regex`; this crate has no regular-expression dependency and may not add one, and every pattern the source can return has the shape `<prefix>=(\d+)` or a bare `(\d+)`. |
| (none) | `bracket_sequence` | `AASequence::toBracketString(false, true)`, the `Peptide` column's middle part. `AASequence.h` is not this package's header and the crate's `AASequence` has no bracket rendering, so the one form the `.pin` writer needs is reimplemented here. |
| (none) | `write` | `prepare_pin` plus a stream write. |
| (none) | `PinOptions`, `ReadOptions`, `PreparedPin`, `StampReport`, `PinDocument` | Grouped arguments, returned diagnostics and explicit ceilings. |
| (none) | `MAX_CHARGE_COLUMNS`, `MAX_HITS`, `MAX_ROWS`, `MAX_COLUMNS`, `MAX_BYTES` | The default ceilings. |

## Preserved source conventions

- **The column contract.** `standard_feature_set` returns exactly
  `SpecId, Label, ScanNr, ExpMass, CalcMass, mass, peplen`, then `charge<c>` for
  every `c` in `min_charge..=max_charge`, then `enzN, enzC, enzInt, dm, absdm`.
  A `max_charge` below `min_charge` yields no charge column at all, because the
  source's loop body never runs.
- **`bool` features read as `1`/`0`.** `stampPinFeaturesOnHits` stamps
  `charge<i>`, `enzN` and `enzC` with a C++ `bool`, which `DataValue` has no
  alternative for and therefore promotes to an `int`. The port stamps `i32`
  values so the column text is identical.
- **Float features use `StringUtils::toStr(double)`.** That is what
  `DataValue::toString()` calls, and the crate already reproduces it as
  `crate::param::value::format_float(v, true)`.
- **`getScanIdentifier`'s fallback chain.** Spectrum reference (MS-GF+), then
  `scan=` plus the `spectrum_id` meta value (X!Tandem, one-based), then
  `index=` plus the index. All space, tab, CR and LF characters are then
  removed, as `StringUtils::removeWhitespaces` does.
- **The scan-number pattern comes from the first identification only.**
  `stampPinFeaturesOnHits` derives it from `peptide_ids[0]` and applies it to
  every identification, so a list mixing `scan=` and `index=` identifiers
  yields no scan number for the minority shape. Preserved, and called out at
  the item.
- **An unextractable scan number is stamped as `ScanNr = -1`.** The source
  calls `extractScanNumber(..., no_error = true)`
  (`PercolatorInfile.cpp:420`, `SpectrumNativeIDParser.cpp:80`), which returns
  `-1` instead of throwing, and stamps that value. So the minority shape above
  gets `-1` and its row is still written. Preserved and tested
  (`an_unmatched_scan_identifier_stamps_the_sources_minus_one`). Until this
  fix the port returned `Error::MissingInformation` there, which refused input
  the source writes a row for, while this document claimed the source's
  behaviour was preserved.
- **The last match wins**, and a digit run that does not fit `i32` yields no
  scan number — the source catches its own `ConversionError` and falls through
  to the same `-1`.
- **An existing `CalcMass` is reused**, not recomputed.
- **Isotope-error correction.** An `IsotopeError` meta value (the legacy MS-GF+
  adapter's spelling before OpenMS 2.6) or the current `isotope_error` shifts
  `ExpMass` down by `isotope_error * C13C12_MASSDIFF_U / charge`. The source
  reads the value through `toString()` and then `toFloat`, so a string-typed
  meta value is accepted; so is it here.
- **Skipping rules.** A hit with no `PeptideEvidence`, and a hit whose
  target/decoy status is unknown, are left untouched and reported, with the
  source's two warnings about incomplete `PeptideIndexing`.
- **Only newly introduced keys are recorded** as added, so a caller that stamps
  temporarily can remove its own values without touching input metadata of the
  same name.
- **The flanking residues come from the first evidence only.**
- **`enzN`/`enzC` are computed before the terminus markers are remapped.** The
  source computes both from the raw `[`/`]` markers at `PercolatorInfile.cpp:524`
  and only maps them to Percolator's `-` at line 534, for the `Peptide` column.
  A protein-terminal PSM is therefore reported as non-enzymatic. Preserved; see
  the C++ finding below.
- **`Proteins` joins the accessions with tabs**
  (`PercolatorInfile.cpp:549`), so a hit with several evidences renders one
  field per accession past the declared column count. That is what
  Percolator's trailing protein list is, and it is why the writer's separator
  guard accepts a tab in a trailing `Proteins` column — see the native
  difference below. Until this fix the guard rejected the value stamping had
  just produced and aborted the whole file, which made the writer unusable on
  shared-peptide data; `a_two_accession_psm_writes_percolators_trailing_protein_list`
  now covers it.
- **A short hit is dropped, not padded.** `preparePin_` writes a row only when
  every declared feature has a meta value, and otherwise counts the hit and
  records which names were missing.
- **An empty input yields a header-only file**, with the source's "Creating
  empty percolator input" warning.
- **`load`'s row grouping.** Consecutive rows sharing a `SpecId` become one
  identification. `RT` is the `retentiontime` column times 60, because search
  engines typically write minutes; the spectrum reference is the `ScanNr`
  column, which may be an integer or a whole vendor identifier; the original
  `SpecId` is kept as `PinSpecId`.
- **Charge column spellings.** `charge<c>` and Sage's `z=<c>`, plus Sage's
  `z=other` for matches outside the searched range.
- **`m/z` is `ExpMass / |charge| + PROTON_MASS_U`**, and is left unset when no
  charge could be determined.
- **Decoy reannotation.** With a `decoy_prefix`, the accessions decide:
  mixed → `target+decoy`, all prefixed → `decoy`, none → `target`. Without one,
  the `Label` column is trusted (`1` → target, anything else → decoy).
- **Percolator's terminal bracket spelling.** `]-` becomes `].` and `-[`
  becomes `.[` before the sequence is parsed.
- **The Sage `inf` workaround.** A `ln(-poisson)` value of `inf` is replaced by
  `3.5`.
- **Extra scores are strings.** The source stores the raw column text, not a
  number, and iterates the found names in lexicographic order (a `std::set`),
  which matters for a comparable idXML; `BTreeSet` gives the same order.
- **`DeltaMass`** is `ExpMass - CalcMass`, stamped on every hit.
- **The enzyme table**, verbatim, including the `-` terminus shortcut and the
  `else` branch that makes every bond enzymatic for an enzyme name the source
  does not know — which is how the adapters pass `no_enzyme`.

## Native differences

- **A `.pin` file without a `FileName` column is refused.** The source reads
  that column optionally but then calls `map_filename_to_idx.at(raw_file_name)`
  unconditionally (`PercolatorInfile.cpp:273`), so a file without it looks up
  the never-inserted default name `UNKNOWN` and throws `std::out_of_range`,
  which can be caught by its caller. This port returns
  `Error::MissingInformation` naming the column.
- **Required columns are named, not dereferenced.** `SpecId`, `ScanNr`,
  `Label`, `Peptide`, `Proteins`, `retentiontime`, `ExpMass`, `CalcMass`,
  `FileName` and the configured score column are all looked up once, up front,
  and a missing one is `Error::MissingInformation`. The source reads each with
  `std::unordered_map::at` inside the row loop, which throws
  `std::out_of_range` for an absent key.
- **A duplicate column name is refused.** The source builds an
  `unordered_map`, so a duplicate silently resolves to the last occurrence.
- **The first row always opens an identification.** The source compares the
  row's `SpecId` against a string that starts empty, so a first row with an
  empty `SpecId` calls `pids.back()` on an empty vector.
- **One-hot charge columns are scanned in ascending charge order.** The source
  iterates an `unordered_map` and breaks at the first column equal to `"1"`, so
  a malformed row with two columns set picks an unspecified charge; ascending
  order makes the same row deterministic.
- **A `rank` column below one is refused.** The source computes `rank - 1` into
  an `int`, giving `-1` for zero and a signed overflow for `i32::MIN`; the
  crate's rank is a `u32`, and the decrement is `checked_sub` before it is
  narrowed, so the whole `i32` range the column accepts is an
  `Error::Parse` naming the line rather than a debug panic or a release wrap to
  `i32::MAX`. Covered by
  `an_extreme_rank_column_is_an_error_not_an_overflow`, through both `read` and
  `load`.
- **A non-finite numeric field is refused**, on reading and on stamping. The
  source's `PeptideIdentification` defaults both `mz_` and `rt_` to NaN and
  writes whatever it holds, so an identification with neither would give
  `ExpMass`, `mass`, `dm` and `retentiontime` columns reading `nan`; this port
  returns `Error::MissingInformation` naming the coordinate.
- **A feature name or value containing a CR or LF is refused**, and so is a tab
  in every column but a trailing `Proteins`. The source writes all of them
  through unescaped: a line break ends the row early and a tab shifts every
  column after it. The trailing `Proteins` column is the one place a tab is
  meaningful, because Percolator reads that column as a tab-separated protein
  list; `Proteins` declared anywhere but last is therefore still refused when
  the hit has more than one accession.
- **The Sage sibling paths are derived by stripping suffixes.** The source
  computes them with `StringUtils::substr(pin_file, 0, pin_file.size() - 3)`
  and `pin_file.size() - 15` — unchecked byte-offset arithmetic on the path
  string. A path shorter than `results.sage.pin` wraps the unsigned
  subtraction; the substring count is clamped, producing the wrong sibling
  name. A positive cutoff can split a UTF-8 character if the suffix is absent. This port
  requires the literal suffixes (`strip_suffix`) and returns
  `Error::InvalidValue` otherwise, so neither a short nor a non-ASCII path can
  misbehave. `tests/percolator_infile.rs` exercises both `a.pin` and
  `日本語.pin`.
- **The charge span is bounded.** `standard_feature_set` returns
  `Error::InvalidRange` when `min_charge..=max_charge` would exceed
  `MAX_CHARGE_COLUMNS` (1024) columns. The source builds the list unbounded,
  and the range is caller input.
- **Stamping is atomic.** The computation runs on a temporary that is committed
  only on success, so a refused call leaves the identifications unchanged. The
  source mutates in place and leaves partial state behind on any throw.
- **Warnings are returned, not logged.** `StampReport::warnings`,
  `PreparedPin::warnings` and `PinDocument::warnings` carry what the source
  wrote to `OPENMS_LOG_WARN`.
- **`bracket_sequence` refuses what it cannot render.** An annotation with no
  known delta mass, and a modified `X` residue whose absolute internal mass the
  crate cannot resolve, are `Error::Unsupported`. The source substitutes the
  absolute internal mass for `X` and writes it unsigned; that branch is
  reproduced where the mass is resolvable.
- **`count_enzymatic` counts characters, not bytes.** The source indexes the
  peptide string per byte. For the unmodified sequences it is given the two
  agree; for a caller-supplied non-ASCII string, iterating characters cannot
  split a codepoint.
- **No OpenMP gap.** Neither the header nor its implementation carries a
  `#pragma omp`.

## Checked boundaries and evidence

| Ceiling | Default | What it bounds |
|---|---|---|
| `PinOptions::max_hits` | 5,000,000 | total peptide hits one stamping or writing call may process, counted before anything is touched |
| `PinOptions::max_bytes` | 256 MiB | cumulative owned payload of stamping |
| `ReadOptions::max_rows` | 5,000,000 | data rows of one `.pin` or Sage sibling |
| `ReadOptions::max_columns` | 10,000 | header columns, and the accession list of one row |
| `ReadOptions::max_bytes` | 256 MiB | cumulative owned payload of reading |
| `MAX_CHARGE_COLUMNS` | 1,024 | one-hot charge columns a feature set may declare |

A zero ceiling is `Error::InvalidValue`. `prepare_pin` also bounds the declared
column count. Before CSV materialization, `max_bytes` also caps each input file
(and each Sage sibling), subject to the stricter `CsvFile` hard ceilings.
These per-file staging bounds are separate from the parsed-payload budget.
`max_rows` and `max_columns` are checked after bounded CSV staging; skipped
comments do not consume the data-row limit. The CSV line and field ceilings
also apply. Invalid zero limits are
rejected before input is consumed.

No string is byte-sliced on file-derived data: the Sage sibling paths use
`strip_suffix`, the peptide spelling fix uses `str::replace`, the scan-number
extraction uses `str::find` plus `get`, and every field is parsed whole.
`scan_identifier` uses `String::retain`, which is codepoint-aware.

**Evidence: tier 3, source review.** The expectations in
`tests/percolator_infile.rs` are transcribed class-test literals: the exact
ordered feature sets for charges 2–4 and 3–3; the assembled header beginning
with `SpecId` and ending with `Peptide`, `Proteins`; the written row's
`SpecId` `scan=529`, `ScanNr` `529`, `Label` `1`, `peplen` `7`, `charge2` `1`,
`charge3` `0`, `enzN` `1`, `enzC` `1`, `Peptide` `K.SAMPLER.S` and `Proteins`
`PROT1`; and, for the reader, 9 identifications, 2 filenames, spectrum
references `30381` and `spectrum=2041`, and the `DECOY_` reannotation of the
eighth entry. Transcribed literals detect transcription drift but cannot
falsify a misread algorithm. No C++ was built or executed and no C++ output was
retained, so nothing here is tier 1 or 2. The per-enzyme specificity table, the
scan-number pattern table, the isotope-error shift, the bracket rendering, the
Sage annotation path, the non-ASCII inputs, the resource ceilings, the
multi-accession trailing protein list, the separator refusals, the `ScanNr`
`-1` sentinel and the out-of-range `rank` column are independently derived
(tier 4) from the implementation, since the class test reaches none of them:
its one stored PSM has a single accession, a matching `scan=` identifier and no
`rank` column at all, which is why the writer's own guard could reject its own
output and the `rank` decrement could overflow without any test noticing.

### Section accounting

All five `START_SECTION`s of `PercolatorInfile_test.cpp` are ported.

| Section | Assertion macros | Rust test | One reproduced value |
|---|---|---|---|
| `PercolatorInfile()` | 1 | `the_default_options_are_the_constructed_state` | `MANDATORY_COLUMNS == ["SpecId", "Label", "ScanNr"]` |
| `~PercolatorInfile()` | 0 | `the_default_options_are_the_constructed_state` | nothing is asserted in the source section; it only deletes the pointer |
| `load(pin_file, higher_score_better, score_name, decoy_prefix)` | 5 | `loading_a_sage_pin_file_groups_rows_and_reannotates_decoys` | `pids[6].getSpectrumReference() == "spectrum=2041"` |
| `static StringList getStandardFeatureSet(int, int)` | 10 (two inside `for` loops over the expected lists) | `the_standard_feature_set_is_the_exact_ordered_column_contract` | the charge 2–4 set is `SpecId, Label, ScanNr, ExpMass, CalcMass, mass, peplen, charge2, charge3, charge4, enzN, enzC, enzInt, dm, absdm` |
| `static void store(...)` | 31 | `storing_writes_the_column_contract_and_the_computed_features` | the written row's `Peptide` is `K.SAMPLER.S` |

### Known gaps

- **A `.pin` file this module writes for a multi-accession PSM cannot be read
  back by this module.** The writer emits Percolator's tab-separated trailing
  protein list, so such a row holds more fields than the header; `read` and
  `load` require a rectangular table, because the source does
  (`PercolatorInfile.cpp:239`) — and the source's own `load` splits `Proteins`
  on `;`, which is Sage's spelling, not the `\t` its `store` writes. So
  `store` → `load` is not a round trip in C++ either. This port reproduces
  both halves rather than inventing a reader the source does not have; the
  parse error now names the trailing-protein-list case, and
  `a_surplus_field_is_a_parse_error_that_names_the_protein_list` pins it.
  Recorded as OPENMS-PERCIN-007. Relaxing the column-count check would let a
  genuinely misaligned row through silently, so it is left to the integrator
  together with the upstream fix.

### C++ findings

Recorded for `OpenMS_CPP_ISSUES.md`; the integrator owns that file. Suggested
IDs are noted so the Rust test comments can cite them.

- **OPENMS-PERCIN-001 — `load` throws on a `.pin` file with no `FileName`
  column.** `PercolatorInfile.cpp:245` fills `map_filename_to_idx` only inside
  `if (file_name_column_index >= 0)`, and line 273 then calls
  `map_filename_to_idx.at(raw_file_name)` unconditionally with the initial
  default name `"UNKNOWN"`. `std::map::at` throws `std::out_of_range`, which this function does not
  translate to a file-format error. A caller may catch it; process termination
  is not inevitable. Confirmed in the shared log as CPP-161.
  Proposed fix: insert `"UNKNOWN"` on first use, or set the merge index only
  when the column exists. Rust handling: `Error::MissingInformation` naming the
  column.
- **OPENMS-PERCIN-002 — an unchecked path subtraction builds the Sage sibling
  paths.** `PercolatorInfile.cpp:95` and `:100` compute
  `pin_file.size() - 3` and `pin_file.size() - std::string("results.sage.pin").length()`
  without checking the path length. A path shorter than fifteen bytes wraps the unsigned requested
  substring length. `StringUtils::substr` delegates to `std::string::substr`,
  which clamps that count and keeps the whole input, producing an incorrect
  sibling name rather than an out-of-bounds read or a subtraction exception.
  When the expected suffix is absent, a positive byte cutoff can also split
  a UTF-8 character. Proposed fix: strip the literal suffixes and report a mismatch.
  Rust handling: `strip_suffix` plus `Error::InvalidValue`; tested with `a.pin`
  and `日本語.pin`.
- **OPENMS-PERCIN-003 — protein-terminal PSMs are reported as non-enzymatic.**
  `PercolatorInfile.cpp:524` computes `enzN` and `enzC` from
  `PeptideEvidence`'s `[` and `]` terminus markers, and only line 534 maps them
  to the `-` that `isEnz_` recognises — and that remapping is used solely for
  the `Peptide` column. So a peptide at a protein terminus gets `enzN = 0`
  where Percolator's own convention would give `1`, feeding Percolator a
  systematically wrong feature for exactly the peptides at protein ends.
  Proposed fix: remap the markers before computing `enzN`/`enzC`. Rust
  handling: reproduced, tested, and documented at the item.
- **OPENMS-PERCIN-004 — `pids.back()` on an empty vector for an empty first
  `SpecId`.** `PercolatorInfile.cpp:268` compares the row's `SpecId` against a
  `spec_id` that starts as the empty string, so a first data row whose `SpecId`
  field is empty takes the `else` path and line 273 dereferences the back of an
  empty `PeptideIdentificationList`. Proposed fix: track whether an
  identification is open. Rust handling: the first row always opens one.
- **OPENMS-PERCIN-005 — `retentiontime` and the other `.at`-read columns are
  undocumented requirements.** `PercolatorInfile.cpp:257`, `:267`, `:274`,
  `:285`, `:286`, `:288` and `:370` all use `std::unordered_map::at` on the
  column index map, so a `.pin` file missing any of them throws
  `std::out_of_range` instead of the documented parse error. A caller may catch it. The header's `@throws` documents only
  `Exception::ParseError` for a wrong column count. Proposed fix: check the
  header up front and throw `ParseError`. Rust handling: all of them are
  checked up front, as `Error::MissingInformation`.
- **OPENMS-PERCIN-006 — the winning one-hot charge column is unspecified.**
  `PercolatorInfile.cpp:293` iterates `col_name_to_charge`, an
  `unordered_map`, and breaks at the first column equal to `"1"`. A row with
  more than one one-hot column set therefore yields a charge that depends on
  the hash order, so the same file can read differently across builds.
  Proposed fix: iterate in ascending charge order. Rust handling: ascending
  charge order.
- **OPENMS-PERCIN-007 — `store` and `load` disagree about the `Proteins`
  column.** `PercolatorInfile.cpp:549` writes the accessions of one PSM joined
  with `\t`, which is Percolator's trailing protein list and therefore makes a
  data row wider than the header. `load` then rejects exactly that row:
  `:239` throws `Exception::ParseError` when the field count differs from the
  header, and `:290` splits the `Proteins` field on `;` rather than `\t`. So
  OpenMS cannot read back the `.pin` files it writes for shared peptides — the
  normal case — and a Sage `.pin` (semicolon-separated, rectangular) is the
  only shape `load` accepts. Proposed fix: treat the final `Proteins` column as
  a variable-length tab-separated list on reading, and accept `;` inside it for
  Sage. Rust handling: the writer emits the source's tab-joined list, the
  reader keeps the source's rectangular requirement and now names this case in
  the parse error, and the gap is recorded above.
