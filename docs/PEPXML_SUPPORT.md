# Native pepXML interchange

`format::pepxml` reads and writes Trans-Proteomic Pipeline pepXML search results
with the crate's existing identification types. It ports
`src/openms/include/OpenMS/FORMAT/PepXMLFile.h` and
`src/openms/source/FORMAT/PepXMLFile.cpp` at OpenMS4-core revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The module is gated on the existing
**`idxml`** feature, which is enabled by default and already provides the
`quick-xml` dependency and the identification stack this format maps onto. It
calls no C++ code and starts no threads; `PepXMLFile.cpp` contains no
`#pragma omp`, so there is no parallelism gap to record.

```rust
use openms::format::pepxml;

# fn example() -> openms::Result<()> {
let mut options = pepxml::ReadOptions::default();
options.experiment_name = "PepXMLFile_test".into();
let document = pepxml::load_with_options("run.pepxml", &options)?;
for warning in &document.warnings {
    eprintln!("pepxml: {warning}");
}
pepxml::store("copy.pepxml", &document)?;
# Ok(())
# }
```

[`PepXmlDocument`](../src/format/pepxml.rs) carries `protein_identifications`,
`peptide_identifications`, `warnings` and `suppressed_warnings`. The two
identification collections are linked by
`ProteinIdentification::identifier`, exactly as the source's two output
containers are. `warnings` is the native channel for what the source writes to
its log stream through `XMLHandler::error` and `::warning`; the source loses
those on a non-interactive run, and nothing in this port writes to a log.

## API mapping

Every public member declared in `PepXMLFile.h`:

| C++ member | Rust counterpart |
|---|---|
| `PepXMLFile()` | No object to construct. `ReadOptions::default()` / `WriteOptions::default()` carry the constructor's state (the cached hydrogen element becomes an `EmpiricalFormula::parse("H")` per read; `keep_native_name_`, `analysis_summary_`, `search_score_summary_` and the empty preferred-modification lists become option fields and reader state) |
| `~PepXMLFile()` | Not ported: the source destructor is `= default` and the class owns no external resource. `PepXmlDocument` and the reader state are ordinary owned values |
| `void load(const std::string&, std::vector<ProteinIdentification>&, PeptideIdentificationList&, const std::string& experiment_name, const SpectrumMetaDataLookup& lookup)` | `load_with_options`, `load_with_registry`, `read_with_options`, `read_with_registry`; `experiment_name` is `ReadOptions::experiment_name` and `lookup` is `ReadOptions::lookup`. The out-parameters become the returned `PepXmlDocument` |
| `void load(const std::string&, std::vector<ProteinIdentification>&, PeptideIdentificationList&, const std::string& experiment_name = "")` | `load`, `read`; the defaulted name is the default empty `ReadOptions::experiment_name`, and an absent lookup is the default empty `SpectrumIndex` |
| `void store(const std::string&, std::vector<ProteinIdentification>&, PeptideIdentificationList&, const std::string& mz_file = "", const std::string& mz_name = "", bool peptideprophet_analyzed = false, double rt_tolerance = 0.01)` | `store`, `store_with_options`, `write`, `write_with_options`. `mz_file`, `mz_name` and `peptideprophet_analyzed` are `WriteOptions` fields; `rt_tolerance` moves to `SpectrumIndex::set_rt_tolerance`, because this port does not build the lookup itself. The source's non-const `std::vector&` arguments are read-only here |
| `void keepNativeSpectrumName(bool keep)` | `ReadOptions::keep_native_spectrum_name` and `WriteOptions::keep_native_spectrum_name`. Its `@note` is preserved: a `pepxml_spectrum_name` meta value carrying the original TPP-format spectrum name is added to each `PeptideIdentification` |
| `void setPreferredFixedModifications(const std::vector<const ResidueModification*>&)` | `ReadOptions::preferred_fixed_modifications`, a `Vec<Arc<ResidueModification>>` obtained from `ModificationsDB::get_modification_handle` or `find_handles` |
| `void setPreferredVariableModifications(const std::vector<const ResidueModification*>&)` | `ReadOptions::preferred_variable_modifications` |
| `void setParseUnknownScores(bool)` | `ReadOptions::parse_unknown_scores` |
| `protected void onStartElement(const char16_t*, const Internal::XMLAttributes&)` | Not public API. Internal `ReaderState::start`, dispatching in the same order over the same element names |
| `protected void onEndElement(const char16_t*)` | Not public API. Internal `ReaderState::end` |
| Inherited from `public Internal::XMLFile`: `getVersion()`, `isValid(filename, os)` | Not ported. `getVersion` returns the handler's `"1.12"` and `isValid` validates against the bundled `pepXML_v114.xsd` through the Xerces schema validator. This port neither carries a schema version string nor validates against a schema; XML-schema validation lives in `format::mzml_schema` behind the `mzml-schema` feature and has no pepXML entry point |

The private nested class and the private helpers, because they define the
behaviour the public surface exposes:

| C++ private member | Rust counterpart |
|---|---|
| `struct AminoAcidModification` (private, so not public API) | Private `HeaderModification`, with the same fields: `aminoacid_`/`amino_acid`, `massdiff_`/`mass_diff`, `mass_`/`mass`, `is_variable_`/`variable`, `description_`/`description`, `terminus_`/`terminus`, `is_protein_terminus_`/`protein_terminus`, `term_spec_`/`term` (`NUMBER_OF_TERM_SPECIFICITY` becomes `None`), `errors_`/`warnings`, `registered_mod_`/`resolved` |
| `AminoAcidModification(aminoacid, massdiff, mass, variable, description, terminus, protein_terminus, preferred_fixed, preferred_var, tolerance)` | `HeaderModification::new` with the same argument order and the same `Exception::MissingInformation` |
| `AminoAcidModification() = delete`, copy ctor, `operator=`, virtual destructor | Not needed: the Rust type has no default, derives `Clone` and needs no virtual destructor |
| `lookupModInPreferredMods_(preferred, aminoacid, massdiff, description, term_spec, tolerance)` | `lookup_preferred`, with the same two passes (full identifier, then residue/terminus/mass window) |
| `toUnimodLikeString()` | Not ported: the whole pinned tree contains only its declaration and its definition, no call. Its output — a signed full-precision mass and then, in parentheses, an optional `Protein ` prefix, the upper-cased terminus plus `-term`, and the upper-cased amino acid, so `+57.021464 (C)` or `-18.010565 (Protein N-term E)` — is neither a UniMod identifier nor the `getFullId()` that every consumer of this class uses |
| `getDescription()` | `ResolvedModification::full_id` — the source returns `registered_mod_->getFullId()`, not `description_` |
| `isVariable()`, `getMassDiff()`, `getMass()`, `getTerminus()`, `getAminoAcid()`, `getErrors()`, `getRegisteredMod()` | `HeaderModification` fields, read directly |
| `readRTMZCharge_(attributes)` | `ReaderState::read_rt_mz_charge` |
| `lookupAddFromHeader_(modification_mass, modification_position, header_mods)` | `lookup_from_header` |
| `makeScanMap_()` | Not ported: declared in the header, never defined and never called. `scan_map_` is constructed and cleared but never written or read, and `use_precursor_data_` is never used at all |
| `mod_tol_`, `xtandem_artificial_mod_tol_` | Public `MODIFICATION_TOLERANCE` (0.002) and `XTANDEM_ARTIFICIAL_MODIFICATION_TOLERANCE` (0.0005) |
| `getIsotopeErrorsFromIntSetting_` | Not ported: commented out in the source |
| parser state (`proteins_`, `peptides_`, `lookup_`, `exp_name_`, `search_engine_`, `native_spectrum_name_`, `experiment_label_`, `swath_assay_`, `status_`, `analysis_summary_`, `search_score_summary_`, `search_summary_`, `wrong_experiment_`, `seen_experiment_`, `checked_base_name_`, `has_decoys_`, `decoy_prefix_`, `current_base_name_`, `current_ms_run_path_`, `current_proteins_`, `params_`, `enzyme_`, `current_peptide_`, `current_analysis_result_`, `peptide_hit_`, `current_sequence_`, `rt_`, `mz_`, `scannr_`, `charge_`, `search_id_`, `prot_id_`, `date_`, `hydrogen_`, `hydrogen_mass_`, `current_modifications_`, `fixed_modifications_`, `variable_modifications_`, `preferred_*_modifications_`) | Private `ReaderState` and `RunState` fields of the same names. `enzyme_cuttingsite_` is declared and never used and has no counterpart |

`SpectrumMetaDataLookup` and its `SpectrumLookup` base are separate unported
headers. `SpectrumIndex` re-creates only what pepXML needs and is deliberately
**not** named `SpectrumLookup`:

| C++ | Rust |
|---|---|
| `SpectrumMetaDataLookup::readSpectra(spectra)` | `SpectrumIndex::read_spectra(&[MSSpectrum])`, or `SpectrumIndex::from_entries` when the caller already has the metadata |
| `SpectrumLookup::empty()` | `SpectrumIndex::is_empty`, plus `len` |
| `SpectrumLookup::findByNativeID(id)` | `SpectrumIndex::find_by_native_id` — `Option`, not `Exception::ElementNotFound` |
| `SpectrumLookup::findByScanNumber(n)` | `SpectrumIndex::find_by_scan_number` |
| `SpectrumLookup::findByRT(rt)` | `SpectrumIndex::find_by_rt`, with `rt_tolerance`/`set_rt_tolerance` |
| `SpectrumMetaDataLookup::getSpectrumMetaData(index, meta)` | `SpectrumIndex::get(index) -> Option<&SpectrumMetaData>` |
| `SpectrumLookup::isNativeID(id)` | `SpectrumIndex::is_native_id` |
| `SpectrumLookup::extractScanNumber` with `default_scan_regexp` | Private `extract_scan_number` |
| `findByIndex`, `addReferenceFormat`, `findByReference`, `setScanRegExp_`, the `SpectrumMetaData` precursor fields | Not ported |

## Preserved source conventions

- **m/z is recomputed, never read.** `readRTMZCharge_` computes
  `(precursor_neutral_mass + hydrogen_mass * charge) / charge`. The hydrogen mass
  is monoisotopic or average according to `search_summary/@precursor_mass_type`,
  and **average is assumed when `search_summary` is absent**. That is why
  `538.605` and `585.3166250319` are what the upstream test expects.
- **`end_scan` is ignored.** A merged spectrum query is not supported; only the
  start scan is parsed.
- **Ranks are converted.** pepXML `hit_rank` is 1-based and OpenMS ranks are
  0-based, in both directions.
- **`massdiff` becomes an isotope error.** `round(massdiff / C13C12_MASSDIFF_U)`
  is stored as the `isotope_error` meta value on every hit.
- **Scores are engine-specific.** `expect`, `mvh`, `xcorr` and `fval` set the
  hit score and the identification's score type and direction; `hyperscore` and
  `nextscore` become meta values; `xcorr` is ignored as a score for MyriMatch,
  which reports one of its own. PSI-MS accessions are attached exactly as the
  source attaches them, including its choice of `MS:1001155` (SEQUEST:xcorr) as
  the generic xcorr term and `MS:1001330` (X!Tandem:expect) for MSFragger.
  `analysis_result/peptideprophet_result` overwrites the search score unless
  InterProphet already did, and `interprophet_result` always overwrites.
- **Modification masses, not names, are the reference.** `mod_aminoacid_mass`
  and `mod_nterm_mass`/`mod_cterm_mass` are matched against the header
  declarations' absolute `mass` within `MODIFICATION_TOLERANCE` (0.002 Da), fixed
  declarations before variable ones, and a header declaration additionally has
  to list the annotated residue in its `aminoacid` string. Only then is the
  registry consulted, least-specific specificity first.
- **A header declaration is resolved in four steps**: the caller's preferred
  list (by full identifier, then by residue/terminus/mass), the `description`
  attribute, a registry mass search (`ANYWHERE` first only when the declaration
  left the terminus open), and finally an anonymous mass annotation. A
  declaration with `massdiff="0"` that nothing explains is dropped entirely.
- **`protein_terminus` is read case-sensitively.** The schema allows only `""`,
  `n` or `c`, but many tools write `Y`/`N`; lower-case `y` and `c` and exact `n`
  mark a protein terminus, exact `N` clears it.
- **`mass == massdiff` is repaired.** Such a declaration is wrong; the absolute
  mass is recomputed from `massdiff` plus the internal-to-terminus mass, or plus
  the residue's internal monoisotopic mass when no terminus is given.
- **Artefact terminal modifications are dropped.** A `terminal_modification`
  whose `|massdiff|` is below `XTANDEM_ARTIFICIAL_MODIFICATION_TOLERANCE`
  (0.0005) is discarded, because some X!Tandem versions annotate a spurious
  near-zero fixed terminal modification that interferes with the real ones.
- **Explicit annotations beat implicit fixed ones.** Modifications annotated on
  a `search_hit` are applied first; the header's fixed declarations then fill
  every residue and terminus that is still unannotated.
- **Both `search_id` numbering schemes work.** pepXML numbers searches either
  per `msms_run_summary` (TPP) or sequentially across them (ProteomeDiscoverer),
  so `search_id` may exceed the number of runs recorded for the current run
  summary; the last recorded run is then the right one.
- **The run date advances one second per `search_summary`.** idXML identifies a
  run by search engine plus date, and two runs must differ, so the source's
  work-around is reproduced verbatim.
- **The MS run path prefers the spectra file.** `msms_run_summary/@base_name`
  plus `@raw_data` (a leading dot is added only when missing) becomes the run's
  primary MS run path; `search_summary/@base_name` is used only when the run
  summary had none, because it can name the identification result file.
- **The Mascot work-around is preserved.** A `msms_run_summary` with an empty
  `base_name` defers the experiment check to `search_summary`, which rolls the
  run back when it does not match.
- **Duplicate protein hits are removed per run**, keeping the first occurrence
  of each accession.
- **`stricttrypsin` is MSFragger's name for `Trypsin/P`** in both
  `sample_enzyme` and `enzymatic_search_constraint`.
- **`min_number_termini` is the specificity.** `EnzymaticDigestion::Specificity`
  is numbered so an engine can report the number of required termini: 0 none,
  1 semi, 2 full.
- **Comet's `peptide_mass_units`** is 0 amu, 1 mmu, 2 ppm; only 2 selects a
  relative precursor tolerance. `fragment_bin_tol` is halved.
- **The writer's constants are the source's.** `date="2007-12-05T17:49:46"`,
  the PeptideProphet `analysis_timestamp` at `17:49:52`, `num_tot_proteins="1"`,
  `num_matched_ions="0"`, `tot_num_ions="0"`, `massdiff="0.0"`,
  `num_missed_cleavages="0"`, `is_rejected="0"` and
  `protein_descr="Protein No. 1"` are all written unconditionally, and
  `num_tol_term` is 2 only when the leading evidence's preceding residue is R or
  K and the enzyme is Trypsin.
- **The writer emits one `spectrum_query` per hit**, not per identification, so
  a multi-hit identification becomes several single-hit queries sharing a scan
  number. The `spectrum` attribute is `base_name.scan.scan.charge`, which
  PeptideProphet requires or the TPP reports a parsing error, and every `.` in
  `base_name` is replaced by `_` so the charge can be read back.
- **Numbers are formatted two ways.** Values the source passes through
  `precisionWrapper` use the full-precision printer (fixed 15 fractional digits
  inside `[1e-2, 1e4)`, shortest round-trip scientific outside it); values
  streamed directly use the stream's precision of 15, which is C's `%.15g`. Both
  are reproduced, the first by `crate::param::value::format_float` and the
  second by a local `general`.
- **Only the first `ProteinIdentification`'s search parameters are written**, and
  the search engine must be the same for all of them. The source warns; this
  port has no log stream in the writer and states the restriction here.
- **`XTandem` is written as `X! Tandem` and `Mascot` as `MASCOT`**; every other
  engine name passes through.

## Native differences

Each of these is also stated at the item in `src/format/pepxml.rs`.

- **A uniquely resolved modification is applied, not dropped.** In
  `mod_aminoacid_mass` the source appends the registry match only inside its
  `mods.size() > 1` branch (`PepXMLFile.cpp:1745`), so a modification that
  resolves to *exactly one* registry entry and was not declared in the pepXML
  header is silently lost. This port applies it. Reported as PEPXML-005.
- **Nothing is registered globally.** The source registers an unexplained mass
  as a `MASS_ONLY` record in the process-wide `ModificationsDB`, which mutates
  global state and grows for the life of the process. This port produces an
  anonymous `MassTag` with the *identical* identifier — `M[+1.0]`, `.n[+2.5]`,
  `.c[+3.4]`, spelled by the same signed full-precision rule — and never mutates
  a registry, so reading the same document twice yields identical results.
- **Numeric annotations may be named.** A resolved anonymous mass is attached
  through the crate's numeric-bracket convention, which resolves a known
  modification within the *written* precision (`AASequence` documents this for
  every numeric annotation). The source keeps its `MASS_ONLY` record. The mass
  is preserved either way.
- **A modification that cannot sit where it was annotated is kept as a mass.**
  The source assigns the resolved `ResidueModification*` directly and never
  re-checks it, so a terminal mass search that returns, say,
  `Ammonia-loss (N-term C)` for an E-terminal peptide annotates the wrong
  chemistry. This port prefers a candidate whose origin is compatible with the
  peptide's terminal residue and, if the named attachment is still invalid,
  retains the mass difference as an anonymous annotation and records a
  diagnostic.
- **One-based indices of zero are errors.** `hit_rank="0"`, `position="0"`, a
  `position` past the peptide, `search_id="0"` and an empty run all underflow an
  unsigned subtraction in the source and then index out of bounds
  (`PepXMLFile.cpp:1381`, `:1695`). Reported as PEPXML-001 and PEPXML-002.
- **`assumed_charge="0"` is an error.** The source divides by the charge
  unconditionally and produces an infinite m/z.
- **`end_scan` is actually compared.** The source reads `endscan` from the
  `start_scan` attribute (`PepXMLFile.cpp:935`), so its
  "endscan not equal to startscan" error can never fire. This port compares the
  two attributes as intended and warns. Reported as PEPXML-006.
- **A corrupted `date` is length-checked.** The source tests `date[4]`,
  `date[7]` and `date[10]` with no bounds check (`PepXMLFile.cpp:2004`);
  a shorter value reads past the string. Reported as PEPXML-007.
- **An unknown enzyme name does not abort the load.** `sample_enzyme` guards its
  `ProteaseDB` lookup with `hasEnzyme`, but `search_summary` does not
  (`PepXMLFile.cpp:1947`), so the source throws `Exception::ElementNotFound` and
  cannot read such a file at all. This port leaves `unknown_enzyme`. Reported as
  PEPXML-003.
- **Enzyme names are matched case-insensitively.** `DigestionEnzymeDB` indexes
  every enzyme under its name, its lower-cased name and each synonym, so
  `name="trypsin"` resolves to `Trypsin`; the Rust `ProteaseDB` matches exactly,
  so the lower-cased key is applied inside this module.
- **A fixed protein-C-terminal modification is applied to the C terminus.** The
  source's C-terminal branch tests `C_TERM || PROTEIN_N_TERM`
  (`PepXMLFile.cpp:2134`), a copy/paste slip that sends a fixed
  `PROTEIN_C_TERM` declaration into the internal-residue branch instead, where
  it is applied to every matching residue. Reported as PEPXML-008.
- **The C-terminal diagnostic names the C-terminal mass.** The source prints
  `mod_nterm_mass`, a different and possibly uninitialised local, and says
  "N-terminal" (`PepXMLFile.cpp:1642`). Reported as PEPXML-009.
- **The cleavage expression is split safely.** The writer splits the enzyme's
  expression on `)` and the source then reads `sub_regex[1]` without a bounds
  check (`PepXMLFile.cpp:427`); `StringUtils::split` returns one element when
  the separator is absent, which a user-defined enzyme's expression — the bare
  cut residues — always is. Reported as PEPXML-010.
- **`digestion_regex` stays empty for a named enzyme.** idXML refuses a nonempty
  custom expression and the enzyme name determines the rule, so the reader
  stores only the name and the writer reads the expression back from
  `ProteaseDB`. Only a user-defined enzyme stores an expression.
- **`search_engine_version` is retained.** The source reads it into a local,
  uses it only to detect MSFragger and drops it with a `TODO`; this port also
  stores it on `ProteinIdentification::search_engine_version`.
- **Flanking residues are typed.** pepXML spells a protein terminus `-`, which
  becomes `FlankingResidue::NTerminus`/`CTerminus` on read and `-` again on
  write. The source stores the character itself, so an OpenMS marker character
  reaching `store` is written verbatim there and as `-` here.
- **`>` is not escaped and `&`, `<` and `"` are.** The source escapes nothing,
  which produces a malformed document for an accession containing `&`; `>` is
  legal inside an attribute value and the retained outputs contain
  `description="Glu->pyro-Glu (N-term E)"`, so it is left alone.
- **Store is atomic.** `store_with_options` serializes into a sibling file and
  publishes only on success, so a failure never leaves a truncated pepXML. The
  source writes directly into an `ofstream`.
- **The writer never reads spectra.** The source `store` loads `mz_file` to
  build its retention-time lookup; here `mz_file` supplies only `base_name` and
  `raw_data`, and `WriteOptions::lookup` supplies the metadata. This keeps the
  module independent of the `mzml` feature.
- **Diagnostics are returned, and are bounded in bytes.** Everything the
  source passes to `XMLHandler::error`/`::warning` is collected in
  `PepXmlDocument::warnings`, bounded by `ReadOptions::max_warnings` **and** by
  `ReadOptions::max_warning_bytes`, with the overflow counted in
  `suppressed_warnings`. None of them is fatal, exactly as in the source. A
  message that quotes the peptide quotes at most its first 32 residues and then
  its byte length, where the source logs the whole sequence: several of these
  messages are emitted once per annotation, so quoting a file-derived peptide of
  unbounded length in full let 114,941 bytes of input retain 64,091,000 bytes of
  diagnostics. Building and retaining a message is charged against the meter.
- **One `search_hit`'s annotations are planned and then installed in one pass.**
  The source stores a resolved `ResidueModification*` per slot, so its loop over
  `mod_aminoacid_mass` entries can consult and mutate the peptide as it goes;
  every `AASequence` setter here clones the peptide and recomputes its formula
  and its mass from every residue, which made one hit cost its peptide length
  per annotation — measured at 121.2 s for a 946,785 byte document with one
  17,400-residue peptide carrying 17,400 annotations. The annotations are now
  decided first, against the free slots rather than against intermediate
  peptides, and then installed together: the planned mass annotations are
  respelled as bracket annotations on the peptide and parsed in one pass, which
  resolves each of them through the same function the setter would have called;
  the planned named modifications are resolved with the same lookup the setter
  performs and installed in one call; and the at most two terminal annotations
  keep their own setters. The same document now loads in 0.18 s. A `peptide`
  attribute that carries an annotation of its own cannot be respelled and keeps
  one setter per modification, which is charged its real cost, so a peptide long
  enough for the product of sites and length to matter is refused rather than
  run.
- **Ambiguous residue codes have no monoisotopic mass.** The internal residue
  mass a `mass == massdiff` repair and a `mod_aminoacid_mass` difference need is
  computed from the crate's chemistry, which has no monoisotopic composition for
  B/Z/X where the source supplies an averaged placeholder. The repair then keeps
  the declared mass and records a diagnostic; the source dereferences a null
  `Residue*`.
- **`makeScanMap_`, `scan_map_`, `use_precursor_data_` and
  `enzyme_cuttingsite_` have no counterpart**, being dead members of the source
  class. Reported as PEPXML-011.

## Checked boundaries and evidence

Every ceiling is checked before the corresponding allocation, following
`src/identification/run_mapping.rs`, so exceeding one leaves no partially built
document and the input is untouched:

| `ReadOptions` | Default | Bounds |
|---|---|---|
| `max_input_bytes` | 512 MiB | Decoded input, and the work/allocation meter that all other accounting draws from |
| `max_elements` | 20,000,000 | XML elements |
| `max_depth` | 64 | XML nesting |
| `max_identifications` | 2,000,000 | `search_result` records |
| `max_hits` | 10,000 | `search_hit` records per `search_result` |
| `max_modifications` | 100,000 | Header declarations per run, and annotations per hit |
| `max_protein_hits` | 5,000,000 | `search_hit` plus `alternative_protein` references per run |
| `max_warnings` | 1,000 | Retained diagnostics; the rest are counted |
| `max_warning_bytes` | 1 MiB | Total length of the retained diagnostics |

`WriteOptions::max_output_bytes` (512 MiB) and `max_identifications` bound
serialization. Zero-valued or over-large limits are rejected before any input is
read.

The work and allocation budget is derived from `max_input_bytes`, which is a
ceiling and not a measurement: with the default ceiling a small document would
be entitled to gigabytes of work. Once the input is decoded the budget is
lowered to what its length entitles it to — 2,048 work units and 512 allocation
units per decoded byte, never below the 16 MiB floor the ceiling-derived budget
already has, and never raised. Total work is therefore linear in the input
length, which is what bounds the per-`search_hit` annotation work described
above. Reading the upstream fixtures spends 13 to 22 work units per byte and 7
to 15 allocation units per byte; a document that is almost entirely one long
peptide spends an order of magnitude more, which is the headroom those two
factors leave.

No index, slice or arithmetic on file-derived data is unchecked. Positions and
ranks are `checked_sub`; an annotated position addresses the peptide through the
`Vec<char>` walked once per `search_hit`, never a byte offset and never a fresh
`chars().nth` walk per annotation; every number must parse and be finite; the
isotope-error and nominal-mass conversions are range-checked before the
`f64`-to-integer cast. No string this module did not construct is ever
byte-sliced — the bracket annotations a plan is respelled into are built from
the module's own numeric spellings, and a spelling that is not plain decimal is
left to the setter instead — and `tests/pepxml.rs` reads a pepXML whose
accessions are Japanese, stores a document whose base name is Japanese and loads
a Japanese path.

**Evidence.** Tier 1 on the writer and tier 3 on the reader; see
`tests/data/pepxml_provenance.json` for the hashes, the 28 source anchors and
the full statement.

- **Tier 1, executed differential.** The three upstream `store` sections compare
  against three retained C++ outputs — `PepXMLFile_test_out.pepxml`,
  `PepXMLFile_test_out_1.pepxml` and `PepXMLFile_test_out_mzML.pepxml`, all
  copied into `tests/data/` and hashed. `tests/pepxml.rs` reproduces each of
  those invocations and compares against the retained bytes with the sections'
  own `FuzzyStringComparator` policy: a number pair passes when the absolute
  difference is within 1e-7 **or** the ratio is within 1 + 1e-7, the `or` read
  from `FuzzyStringComparator.cpp:616-684` and itself unit-tested. Byte equality
  is not the contract and could not be — the two retained files disagree with
  each other on `2269.1954242316001` against `2269.195424231600555` for the same
  peptide, so they were produced by two different OpenMS builds.
- **Tier 3, source review.** Every load expectation is a literal transcribed
  from `PepXMLFile_test.cpp`.
- **Tier 4, independently derived.** The resource ceilings, the one-based index
  and non-ASCII boundary tests, the malformed-document rejections, the
  error-variant choices and the diagnostics channel; no upstream section reaches
  them.

No C++ was built or executed for this package.

## Ported test sections

All ten `START_SECTION`s of `PepXMLFile_test.cpp` are ported; none is merely
mapped. Macro counts are assertion macros inside the section, comments excluded.

| # | Section | Macros | Rust test in `tests/pepxml.rs` |
|---|---|---:|---|
| 1 | `PepXMLFile()` | 1 | `default_construction_and_drop` |
| 2 | `~PepXMLFile()` | 0 | `default_construction_and_drop` (the `drop` and the load-twice equality) |
| 3 | `load(..., const SpectrumMetaDataLookup& lookup)` | 4 | `load_with_a_spectrum_lookup` |
| 4 | `load(..., experiment_name = "")` | 62 | `load_two_search_runs_of_one_experiment`, plus `load_without_an_experiment_name_accepts_every_run` |
| 5 | `[EXTRA] load(..., experiment_name = "")` | 37 | `load_extended_fixture_with_native_names_and_analysis_results` |
| 6 | `store(...)` | 2 | `store_reproduces_the_retained_peptideprophet_output`, `store_reproduces_the_retained_raw_output` |
| 7 | `[EXTRA] store(...)` | 37 | `store_and_reload_the_extended_fixture`, plus `native_spectrum_names_are_only_kept_on_request` for its second block |
| 8 | `store(..., mz_file = "PepXMLFile_test.mzML", ...)` | 1 | `store_with_spectra_metadata_reproduces_the_retained_output` |
| 9 | `keepNativeSpectrumName(bool)` (`NOT_TESTABLE`) | 0 | `native_spectrum_names_are_only_kept_on_request` |
| 10 | `[EXTRA] checking pepxml transformation to reusable identifications` | 7 | `pepxml_search_parameters_survive_an_idxml_round_trip` |

Sections 4, 5, 7 and 10 are above the five-macro threshold and are ported
assertion by assertion. `VALIDATE_TMP_FILES` at the end of the upstream file
validates the temporary outputs against the pepXML schema; this port carries no
schema validator for pepXML and does not reproduce it.

The remaining tests in `tests/pepxml.rs` cover the native boundaries no upstream
section reaches: the resource ceilings, one-based indices of zero, malformed
documents, non-ASCII input, the corrupted-date repair, unknown scores, the
decoy-prefix annotation, the Mascot base-name roll-back, an unknown enzyme, a
custom cleavage expression, the modification diagnostics, the preferred
modification lists, the `SpectrumIndex` queries, the retention-time fallback,
Percolator's missing posterior error probability, the atomic store, and the four
bounds on the annotation pass: that a densely annotated 946,804 byte peptide
loads within a budget only linear in the input, that the diagnostics are bounded
in bytes as well as in number, that a `peptide` attribute carrying its own
annotation still annotates through the per-setter path, and that the work budget
follows the decoded input rather than the declared ceiling.

## Not covered

- No schema validation, and no `pepXML_v114.xsd`/`pepXML_v117.xsd` handling.
- `PepXMLFileMascot.h`, a separate header with its own class test, is unported.
- `analysis_summary` subtrees are skipped, as in the source.
- `spectrum_query/@index`, `search_hit/@num_tot_proteins`,
  `@num_matched_ions`, `@tot_num_ions`, `@calc_neutral_pep_mass`,
  `@num_missed_cleavages`, `@is_rejected` and `@protein_descr` are not read; the
  source ignores them too, and the writer emits its fixed literals.
- `<parameter>` elements outside `search_score_summary` and `search_summary` are
  ignored, and inside `search_summary` only `fragment_bin_tol`,
  `peptide_mass_tolerance`, `peptide_mass_units`, `decoy_search` and
  `decoy_prefix` are interpreted — exactly the source's set.
- `ProteinIdentification` score type, direction and significance threshold are
  never set in either direction; a pepXML carries no protein-level score.
