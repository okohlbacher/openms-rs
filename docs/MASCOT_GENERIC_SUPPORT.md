# Mascot generic format (MGF) support

Port of `FORMAT/MascotGenericFile.h` and `FORMAT/MascotGenericFile.cpp` at
source revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| | |
|---|---|
| Rust | `src/format/mascot_generic.rs` |
| Tests | `tests/mascot_generic.rs` |
| Provenance | `tests/data/mascot_generic_provenance.json` |
| Fixtures | `tests/data/MascotInfile_test.mascot_in`, `tests/data/MascotGenericFile_GNPS.mgf` |
| Feature | none — the module is unconditional, like `crate::format::mgf` |

MGF is the peak-list format Mascot consumes: a sequence of `BEGIN IONS` …
`END IONS` blocks, each with `KEY=value` header lines followed by
whitespace-separated m/z–intensity pairs. Matrix Science documents the format at
`data_file_help.html#GEN`. `MascotGenericFile` is both the reader and the writer
of the *submission* file, which is why it is a `DefaultParamHandler`: its
parameters are the Mascot search settings that become the file's parameter
header.

## Relation to `crate::format::mgf`

The crate already had `src/format/mgf.rs`, a native MGF interchange adapter that
rejects repeated fields, ambiguous charges and stray content, and does not write
a search header. It is a *stricter* reader by design, and it is unchanged.

`crate::format::mascot_generic` is the source-fidelity port: it reproduces what
`MascotGenericFile` accepts and rejects, including the surprising cases below,
and it writes the Mascot parameter header. Use it when parity with OpenMS
matters; use `mgf` when strict interchange matters.

## API mapping — `MascotGenericFile.h`

Every member the header declares, public and protected.

| C++ member | Rust counterpart |
|---|---|
| `class MascotGenericFile` | `MascotGenericFile` |
| base `ProgressLogger` | **not ported**: progress reporting is not threaded through this port. `src/concept/progress_logger.rs` exists; wiring it into file I/O is a separate decision, as for `ImzMLHandler` |
| base `DefaultParamHandler` | `crate::param::DefaultParamHandler`, held as a private field and surfaced through `parameters`, `defaults` and `set_parameters` |
| `MascotGenericFile()` | `MascotGenericFile::new` (uses the shared `ModificationsDB`) and `MascotGenericFile::with_modifications` (caller-owned registry). Fallible, because the parameter tree is resource-checked |
| `~MascotGenericFile()` | `Drop`, implicit |
| `void updateMembers_()` | folded into `set_parameters` / `with_modifications`, which rebuild the specificity-group map atomically with the new parameters; the result is readable through `special_modification_groups` |
| `void store(const std::string& filename, const PeakMap&, bool compact)` | `MascotGenericFile::store`, plus `store_with_options` for explicit ceilings |
| `void store(std::ostream&, const std::string& filename, const PeakMap&, bool compact)` | `MascotGenericFile::store_to` |
| `template <typename MapType> void load(const std::string&, MapType&)` | free `load` / `load_with_options`, plus `MascotGenericFile::load` for a one-to-one API map. The template parameter disappears: this crate has one `MSExperiment` |
| `std::pair<std::string,std::string> getHTTPPeakListEnclosure(const std::string&) const` | `MascotGenericFile::http_peak_list_enclosure` |
| `void writeSpectrum(std::ostream&, const PeakSpectrum&, const std::string& filename, const std::string& native_id_type_accession)` | `MascotGenericFile::write_spectrum`, returning `SpectrumOutcome` instead of writing to `cerr`/`cout` |
| protected `bool store_compact_` | private field, readable through `MascotGenericFile::store_compact` |
| protected `std::map<std::string,std::string> mod_group_map_` | private `BTreeMap`, readable through `MascotGenericFile::special_modification_groups` |
| protected `void writeParameterHeader_(const std::string&, std::ostream&)` | private `write_parameter_header` |
| protected `void writeModifications_(const std::vector<std::string>&, std::ostream&, bool)` | private `write_modifications` |
| protected `void writeHeader_(std::ostream&)` | `MascotGenericFile::write_header_to`, public because `internal:content = header_only` makes it the whole output |
| protected `void writeMSExperiment_(std::ostream&, const std::string&, const PeakMap&)` | private `write_experiment`, which fills `WriteReport` |
| protected `template <typename SpectrumType> bool getNextSpectrum_(std::ifstream&, SpectrumType&, Size& line_number, const Size& spectrum_number)` | `MascotGenericReader::next_block`, driven by `Iterator::next`. The out-parameter line number is the reader's own counter; the boolean return becomes `Option` |
| `@htmlinclude OpenMS_MascotGenericFile.parameters` | this document's parameter table |

Rust-only surface, all native: `MascotGenericReader`, `ReadOptions`,
`CarryOver`, `WriteOptions` (alias of the shared `Limits`), `WriteReport`,
`SpectrumOutcome`, `ExperimentCollector`, the free `read`, `read_with_options`,
`read_into`, `read_into_with_options`, `consume`, `consume_with_options`,
`write`, `store`, and the constants `MAX_WRITTEN_PEAKS`, `MAX_SEQ_ENTRIES`,
`MAX_MODIFICATIONS`, `TITLE_KEY`, `SEQ_KEY`, `SCAN_ID_KEY`,
`GNPS_SPECTRUM_ID_KEY`, `UNKNOWN_NATIVE_ID_TYPE`.

### Parameters

The constructor's `defaults_`, reproduced verbatim including the `internal:`
section the source hides from TOPP users.

| Key | Default | Restriction |
|---|---|---|
| `database` | `MSDB` | |
| `search_type` | `MIS` (advanced) | `MIS`, `SQ`, `PMF` |
| `enzyme` | `Trypsin` | |
| `instrument` | `Default` | |
| `missed_cleavages` | `1` | ≥ 0 |
| `precursor_mass_tolerance` | `3.0` | ≥ 0.0 |
| `precursor_error_units` | `Da` | `%`, `ppm`, `mmu`, `Da` |
| `fragment_mass_tolerance` | `0.3` | ≥ 0.0 |
| `fragment_error_units` | `Da` | `mmu`, `Da` |
| `charges` | `1,2,3` | |
| `taxonomy` | `All entries` | |
| `fixed_modifications` | empty list | every UniMod-backed modification identifier |
| `variable_modifications` | empty list | every UniMod-backed modification identifier |
| `special_modifications` | `Cation:Na (DE),Deamidated (NQ),Oxidation (HW),Phospho (ST),Sulfo (ST)` (advanced) | |
| `mass_type` | `monoisotopic` | `monoisotopic`, `average` |
| `number_of_hits` | `0` (AUTO) | ≥ 0 |
| `skip_spectrum_charges` | `false` | `true`, `false` |
| `decoy` | `false` | `true`, `false` |
| `search_title` | `OpenMS_search` (advanced) | |
| `username` | `OpenMS` (advanced) | |
| `email` | empty | |
| `internal:format` | `Mascot generic` (advanced) | `Mascot generic`, `mzData (.XML)`, `mzML (.mzML)` |
| `internal:boundary` | `GZWgAaYKjHFeUaLOLEIOMq` (advanced) | |
| `internal:HTTP_format` | `false` (advanced) | `true`, `false` |
| `internal:content` | `all` (advanced) | `all`, `peaklist_only`, `header_only` |

The two modification lists are restricted to
`ModificationsDB::getAllSearchModifications`: the full IDs of records carrying
a UniMod record ID, sorted case-insensitively with the shorter of two otherwise
equal names first.

## Preserved source conventions

Every one of these is reproduced, and each is covered by a named test in
`tests/mascot_generic.rs`.

**Reading**

1. *Only blocks are read.* The outer loop looks for `BEGIN IONS` and skips every
   other line, so a parameter header written by `store` — including its global
   `CHARGE=1,2,3` — is not read back
   (`lines_outside_a_block_are_ignored_including_the_parameter_header`).
2. *The peak list starts at the first ASCII-digit-leading line.* Everything
   after it must be a peak or blank; a header line there is a parse error
   (`a_header_line_after_the_peak_list_is_a_parse_error`). A line starting with
   `i` (as in `inf`) is therefore not a peak line at all.
3. *A block with no peak line does not end.* `END IONS` matches no header prefix
   and is ignored, so the block runs on into the next one and the two merge into
   a single spectrum (`a_block_with_no_peak_line_merges_into_the_next`).
4. *Truncation is asymmetric.* A block cut off before its first peak line is
   dropped silently; one cut off after it is a parse error
   (`truncation_before_the_peak_list_is_silent_and_after_it_is_an_error`).
5. *A peak line is whitespace-collapsed, then split on single spaces.* Runs of
   space, tab, CR and LF collapse to one space (the standard allows double
   spaces; tabs are an accepted extension). A line with no whitespace at all is
   the only "does not contain m/z and intensity" error. Only the first two
   fields are converted: a third one — nominally the optional per-peak charge —
   is *ignored without conversion*, which is how the upstream fixture carries
   `#.two.spaces.(allowed.by.the.standard)` comments through the reader, and why
   `100 1 garbage` loads (`load_upstream_infile_fixture`). An earlier revision
   of this list claimed the source parsed and discarded that field; it does not
   touch it.
6. *`PEPMASS=` takes one or two values*; three or more is an error naming the
   count (`pepmass_accepts_one_or_two_fields_and_rejects_three`). Its value is
   tab-substituted but **not** whitespace-collapsed — `simplify` is called only
   on the peak lines — so `PEPMASS=500  10` splits into three fields with an
   empty middle one and is exactly that error
   (`pepmass_is_not_whitespace_collapsed`).
7. *`CHARGE=` strips every `+` and then requires one complete `i32`.* `2+`,
   `+2`, `2` and `-2` parse; a charge list (`1,2,3`) and a trailing `-` (`2-`)
   are conversion errors that propagate out of the load
   (`charge_accepts_a_sign_suffix_but_rejects_lists_and_negatives`).
8. *`TITLE=` has two halves.* A title containing `min` anywhere is treated as a
   Bruker export: the line is split on `,` and the first whitespace token of
   every chunk containing `min` sets the retention time in minutes, so the last
   such chunk wins and **no** `TITLE` meta value is stored at all. If a
   conversion fails, the text between the first and second `=` is stored as
   `TITLE`, truncating a title with a second `=`. `setRT` is called *inside*
   that loop and the single enclosing `catch` does not undo it, so a retention
   time an earlier chunk committed survives a later chunk's failure **and** the
   fallback title is stored beside it: `TITLE=run, 2 min, bad min, 3 min` keeps
   120 s, and the fourth chunk is never visited
   (`a_committed_title_retention_time_survives_a_later_failure`). Otherwise the
   value is the
   text after the first `=` at or after index 4 with `_<native ID>` appended, so
   repeated titles stay distinguishable — unless the stored `TITLE` already
   contains the native ID, in which case a second `TITLE=` line in the same
   block replaces it without the suffix
   (`title_with_minutes_sets_the_retention_time_and_stores_no_title`,
   `title_with_minutes_but_unparsable_falls_back_to_the_first_equals_split`,
   `a_second_title_line_replaces_the_first_without_the_native_id_suffix`).
9. *Key prefixes are matched in the source's order and each assumes `=` right
   after the key.* So `ADDUCT=` and `ION_MODE=` — both present in the upstream
   GNPS fixture — match nothing and are silently dropped, while `IONMODE=`
   is kept (`gnps_library_spectrum`). The value offset is a byte count, and
   `StringUtils::substr` clamps it to the string length, so a line shorter than
   its own key — a bare `NAME` or `MSLEVEL` — yields an *empty* value rather
   than an error (`a_header_line_shorter_than_its_key_yields_an_empty_value`).
10. *Meta key renaming.* `NAME` and `COMPOUND_NAME` both become
    `Metabolite_Name`; `INCHI` becomes `Inchi_String`; `SMILES` becomes
    `SMILES_String`; `SPECTRUMID` becomes `GNPS_Spectrum_ID`; `SCANS` becomes
    `Scan_ID`. `IONMODE`, `SOURCE_INSTRUMENT`, `ORGANISM`, `PI`,
    `DATACOLLECTOR` and `LIBRARYQUALITY` keep their names.
11. *`MSLEVEL=` is read with `std::stoi`*, so a numeric prefix wins and the rest
    is ignored; an unparsable value falls back to MS 2 *and* records
    `MSLEVEL="2"` as a meta value, while an out-of-range value falls back
    silently (`mslevel_falls_back_to_two_when_unparsable`).
12. *`SEQ=` is always a string list*, even for one line, and accumulates in
    order across repeated lines within one query, because the Mascot
    specification makes each `SEQ` an independent sequence filter
    (`many_seq_lines_in_one_query_are_linear`).
13. *Numeric conversion is `std::from_chars` strict.* `toDouble`/`toInt32` skip
    space, tab, CR and LF, consume one leading `+`, and then require a complete
    token: a *second* `+` (`++5`) and an overflowing decimal literal (`1e999`,
    which `from_chars` reports as `result_out_of_range`) are conversion errors,
    not the value Rust's own parser would produce
    (`numeric_conversion_refuses_a_second_plus_and_an_overflowing_literal`).
    `CHARGE=` is the exception: it removes *every* `+` before converting, so a
    second one never reaches the conversion.
14. *Every spectrum is MS 2, centroided, has one precursor and native ID
    `index=<n>`.* MGF is centroided by definition, so the type is asserted
    rather than inferred.

**Writing**

15. *Ten thousand peaks is the ceiling*, with the source's message naming
    profile data as the likely cause
    (`ten_thousand_peaks_are_refused_by_the_writer`).
16. *A precursor m/z of exactly zero skips the spectrum.*
17. *Only MS level 2 is written.* Level 0 warns; every other level is dropped in
    silence (`a_zero_precursor_and_non_ms2_levels_are_skipped_with_a_report`).
18. *A stored `TITLE` is written verbatim*, because it was either parsed from an
    MGF or set to be written to one. Otherwise the title is
    `<m/z>_<RT>_<native ID>_<file name stem>`, the stem stripped of every
    non-alphanumeric character.
19. *Precision.* The default form writes `precisionWrapper`, which is
    `StringUtils::toStr(value, true)`: 15 significant digits for a `double` and
    6 for a `float`, fixed below 1e4 and above 1e-2, scientific outside, always
    with at least one fractional digit — so `1998` becomes `1998.0` and `25.379`
    becomes `25.379000000000001`. The compact form uses five fixed decimals for
    peak m/z, three for peak intensity, and omits zero-intensity peaks.
    The header's `TOL`/`ITOL` use plain `ostream <<`, i.e. `%g` with precision
    6, so `3.0` is written as `3`.
20. *Compact `PEPMASS=`/`RTINSECONDS=` depend on a stream flag.* The source
    streams `fixed` only in the branch that *generates* a `TITLE=` line, and
    again on every compact peak line. A compact spectrum that already carries a
    `TITLE` meta value is therefore written while the stream is still in its
    default float format, so `setprecision(5)`/`setprecision(3)` mean five and
    three *significant* digits: precursor m/z 901.234567 becomes `901.23` and
    retention time 234.5678 becomes `235`. Because the flag lives in the
    `ostream` it stays set, so within one `store` only the spectra before the
    first generated title or written peak line are affected — the second such
    spectrum gets `901.23457`/`234.568`
    (`compact_output_with_a_stored_title_uses_significant_digits`). Recorded as
    `MGF-03` in the C++ issue log: the compact form loses four significant
    digits of precursor m/z for exactly the files that came from an MGF.
21. *`SCANS=`.* With the `UNKNOWN` accession sentinel — which is what an
    experiment with no source file or an empty accession produces — the value is
    the text after the native ID's last `=`. Otherwise it goes through the
    accession table `SpectrumLookup::extractScanNumber` uses, whose failure
    sentinel `-1` is written verbatim
    (`scans_follows_the_native_id_type_accession_table`). The regex token
    iterator collects *every* match and converts only `matches.back()`, with
    `toInt32`: so the last match wins even when its digits overflow 32 bits, an
    earlier convertible match is not a fallback, and the sentinel is written
    with the source's own warning. The WIFF branch collects two subgroups and
    inspects only the final `cycle=`/`experiment=` pair, so an earlier pair with
    an experiment of 1000 or more is never seen — and when the final pair has
    one, `Exception::InvalidValue` is raised and *not* caught (the handler
    catches only `ConversionError`), aborting the whole store
    (`the_last_scan_number_match_wins_before_conversion`).
22. *`FORMAT` stays within the first five lines*, because that is how OpenMS
    recognises its own MGF files when the suffix is not `.mgf`.
23. *Optional header lines.* `COM` only for a non-empty `search_title`,
    `USEREMAIL` only for a non-empty `email`, `DECOY=1` only when `decoy` is
    true, `REPORT=AUTO` when `number_of_hits` is zero.
24. *Modification rewriting.* Each configured modification is looked up in the
    specificity-group map (`Deamidated (N)` → `Deamidated (NQ)`) and the result
    collected into a set, so the output is sorted and a group named twice
    appears once (`duplicate_modification_groups_collapse_to_one_line`).
25. *`skip_spectrum_charges`* suppresses the per-spectrum `CHARGE=` line while
    leaving the header's general `CHARGE` untouched.
26. *Extension check.* `store` to a path requires the `mgf` suffix; a stream
    is written unchecked.

## Native differences

| Difference | Why |
|---|---|
| **Carry-over is off by default.** `getNextSpectrum_` clears only peaks, native ID, `TITLE` and `SEQ`, so retention time, precursor m/z, precursor intensity, precursor charge, MS level and every other meta value leak from one block into the next. `CarryOver::Reset` (the default) starts each block fresh; `CarryOver::Source` reproduces the leak. | The upstream class test asserts only that `SEQ` does not leak, which is the one field the source explicitly resets — the author clearly knew about the others. An inherited charge is not recoverable by a caller and is almost never what the file meant. Both behaviours are tested (`carry_over_reproduces_the_source_bleed_and_reset_prevents_it`). Recorded as `MGF-01` in the C++ issue log. |
| **Non-finite numbers are refused.** The source's `toDouble` accepts `inf` and `nan` and stores them silently. | Every consumer of an `MSSpectrum` rejects a non-finite coordinate, so the failure is moved to the point of parsing. `nan(payload)` is likewise not accepted. |
| **A value offset landing inside a multi-byte character is a parse error.** The source slices raw bytes and produces an ill-formed string. | Rust string slicing at a non-boundary aborts the process, and half a character is not a usable value. Tested with `PEPMASSé=1.0`. An offset *past* the end of the line is not this case: `StringUtils::substr` clamps it, and so does this port (`a_header_line_shorter_than_its_key_yields_an_empty_value`). |
| **A `SEQ=` list is accumulated in the reader**, not round-tripped through the meta value on every line. | The source reads the whole list out of the meta value, appends one entry and writes it back for each `SEQ=` line, which is quadratic: the 100,000-line ceiling costs about five billion string copies, 434 s in release mode on the gate node, for a 1.2 MB block. Accumulating and storing once is 0.04 s and stores the identical list (`many_seq_lines_in_one_query_are_linear`). Recorded as `MGF-04`. |
| **A scan number that overflows 32-bit arithmetic yields the `-1` sentinel and a warning.** The source's `index=` branch computes `toInt32(value) + 1` and its WIFF branch `cycle * 1000 + experiment` in `int`, both of which are signed overflow — undefined behaviour — for a large native ID. | Checked arithmetic cannot wrap, so the port takes the same path as every other extraction failure and says so in the warning. Recorded as `MGF-05` (`the_last_scan_number_match_wins_before_conversion`). |
| **Byte, line, spectrum and peak ceilings.** | Shared with the other text adapters through `Limits`; the source has no ceiling other than the 10,000-peak write guard. Output bytes are counted before the destination file is created, so a refused store leaves no truncated file. |
| **Console output becomes a return value.** `WriteReport` and `SpectrumOutcome` carry the messages the source sends to `cerr`, `cout` and `OPENMS_LOG_WARN`. | The crate has no global streams. |
| **`isdigit(line[0])` is replaced by a byte test.** | Passing a negative `char` to `isdigit` is undefined behaviour; recorded as `MGF-02`. Tested with a whole non-ASCII line. |
| **A non-positive `MSLEVEL` is refused** rather than stored. | The source assigns whatever `std::stoi` returned, and MS level 0 makes the record invalid for every consumer. |
| **Three range filters and a consumer interface.** `ReadOptions::{rt_range, mz_range, intensity_range}` and `consume`. | The source MGF adapter consumes no `PeakFileOptions` and has no consumer interface, only the whole-file `load`. MGF is a streaming format and the rest of the crate offers both, so they are offered here as documented native additions; the semantics mirror `crate::format::dta2d`, where the source does consume the three filters. |
| **`SpectrumLookup::extractScanNumber` is reproduced in subset form.** `METADATA/SpectrumLookup.h` and `METADATA/SpectrumNativeIDParser.h` are separate, unported headers. Only the accession table the MGF writer reaches is implemented, with hand-coded scanning instead of Boost regular expressions. | The writer cannot fill `SCANS=` without it and this package may not add a regular-expression dependency. The table is `MS:1000768/769/771/772/776/1002818` → `scan=`, `MS:1000773/775` → `file=`, `MS:1000774` → `index=` plus one, `MS:1001508` → `scanId=`, `MS:1000777` → `spectrum=`, `MS:1001530` → the last bare number, `MS:1000770` → `cycle * 1000 + experiment`. When `SpectrumLookup` is ported, this local subset should be deleted in favour of it. |
| **No OpenMP.** | The source `writeMSExperiment_` is serial too; the header includes `omp.h` but uses no `#pragma omp`. |

## Checked boundaries and evidence

**Evidence tier 3 (source review).** Every expectation comes from the pinned
`MascotGenericFile_test.cpp` literals or from reading
`MascotGenericFile.cpp`/`.h`; no C++ was built or executed and no C++ output was
retained for MGF, so this is not a tier-1 differential. The upstream suite has
no retained MGF output: its `store` sections compare substrings of an
in-memory stream, and both transcribed fixtures are inputs.

Transcribed literals, all reproduced: one spectrum with nine peaks from
`MascotInfile_test.mascot_in`; the written block
`TITLE=Testtitle_index=0`, `PEPMASS=1998.0`,
`RTINSECONDS=25.379000000000001`, `SCANS=0`, the nine peak lines
`1.0 1.0` … `9.0 81.0`; the generated title
`1998.0_25.379000000000001_index=0_test`; `MODS=Carbamidomethyl (C)`,
`MODS=Phospho (ST)`, `IT_MODS=Deamidated (NQ)`, `IT_MODS=Oxidation (M)`; the
compact block `TITLE=901.23457_234.568_index=250_test`, `PEPMASS=901.23457`,
`RTINSECONDS=234.568`, `SCANS=250`, `890.12346 2345.679`; `Caffeine` under
`Metabolite_Name` with two peaks, precursor 500.0 and charge 1; the GNPS
spectrum's `CCMSLIB00000001547`, `3-Des-Microcystein_LR`, precursor 981.54,
charge 1, 43 peaks and the base peak 599.352783 at intensity 764523.0; and the
`SEQ` single, multiple, programmatic and non-bleeding cases.

All 8 upstream `START_SECTION`s are ported — none is merely mapped. The
destructor section has no assertion upstream; the Rust analogue asserts that the
value is owned, cloneable and droppable without shared state.

**Independently derived (tier 4).** The resource ceilings, the non-ASCII cases,
the multi-byte-offset refusal, the substr clamp, the numeric-conversion
strictness, the `PEPMASS` whitespace rule, the title retention-time commit
order, the compact-precision stream flag, the carry-over comparison, the
block-merge and truncation asymmetry, the `CHARGE` rejection set, the range
filters, the consumer interface, the accession table (including the
last-match-wins, 32-bit and WIFF-pair rules) and the parameter-header line order
are Rust-only checks derived from reading the implementation.

**Second-model review.** An adversarial review by another model raised nine
findings against this port. Eight changed behaviour, each with a named
regression test: the `PEPMASS` whitespace collapsing, the title retention-time
commit order, the compact precision with a stored `TITLE`, the `substr` clamp,
the numeric-conversion strictness, the scan number's 32-bit width and
last-match-wins rule, the WIFF final-pair rule with its uncaught
`InvalidValue`, and the quadratic `SEQ` accumulation. The ninth was a
documentation correction: the third field of a peak line is ignored, not parsed
and discarded. One reported difference was *not* changed — on a libc++ build
the source's `toDouble` accepts `0x10` as 16 through `strtod`, while the
`std::from_chars` path every other platform takes rejects it, as this port
does; reproducing a platform-specific accident is not fidelity.

**Tolerances.** Text output is compared exactly. Precursor m/z and retention
time comparisons on round-trip use 1e-9 absolute, which is the precision the
writer's 15-significant-digit form preserves.

## Deferred

- `ProgressLogger` is not threaded through, so a large file reports no progress.
- `SpectrumLookup` / `SpectrumNativeIDParser` remain unported; see the subset
  note above.
- `MascotRemoteQuery.h` and `MascotInfile.h` are separate headers and are not
  touched here; `MascotGenericFile` replaced `MascotInfile` upstream, and the
  upstream fixture is still named `MascotInfile_test.mascot_in`.
- The writer has no `PeakFileOptions` hook, matching the source.

### Stream failure boundary

`store_with_options` preflights semantic errors and output size before creating
a file. `store_to`, `write_spectrum` and `write_header_to` write directly to the
caller-owned stream: a later validation or I/O error can leave an output prefix,
and `store_to` does not apply the path writer's output-byte ceiling. Callers
requiring atomic publication should stage the stream before publishing it.
