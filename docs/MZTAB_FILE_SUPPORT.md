# MzTab file adapter

Source: `src/openms/include/OpenMS/FORMAT/MzTabFile.h` (208 lines) and
`src/openms/source/FORMAT/MzTabFile.cpp` (3,352 lines) at
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

Rust: `src/format/mztab_file.rs`. Tests: `tests/mztab_file.rs`. Provenance:
`tests/data/mztab_file_provenance.json`.

The data model this adapter reads into and writes out of — every cell type,
every record struct and `MzTab` itself — is `src/format/mztab.rs`, documented in
`docs/MZTAB_SUPPORT.md`. This document covers only the file layer.

Status: **partial**. Both `MzTab`-shaped entry points, `load` and
`store(filename, MzTab)`, are ported with every column family of all seven
sections, and every public and protected member of the header is accounted for
below. The two streaming `store` overloads are not ported: they consume
`MzTab::IDMzTabStream` and `MzTab::CMMzTabStream`, the exporters that turn a
`ConsensusMap` or a list of `ProteinIdentification`/`PeptideIdentification` into
MzTab rows, and neither exporter is in the crate yet. Fifteen further protected
members are declared in the header and defined nowhere in the pinned tree.

---

## API mapping

Every public and protected member the header declares, in header order.

### Public

| C++ | Rust |
|---|---|
| `MzTabFile()` | `MzTabFile::new` / `MzTabFile::default` — all sixteen flags `false` |
| `~MzTabFile()` | not ported: `= default`, and the Rust adapter is sixteen `bool`s with no owned resource, so `std::mem::needs_drop::<MzTabFile>()` is `false` |
| `typedef … MapAccPepType` | not ported: only the fifteen undefined protected statics below use it |
| `void store(const std::string&, const MzTab&) const` | `MzTabFile::store` |
| `void store(const std::string&, const std::vector<ProteinIdentification>&, const PeptideIdentificationList&, bool, bool, bool, const std::string&)` | not ported: builds rows through `MzTab::IDMzTabStream`, which is not in the crate. See "Not ported" below |
| `void store(const std::string&, const ConsensusMap&, const bool, const bool, const bool, const bool, const bool, const bool) const` | not ported: builds rows through `MzTab::CMMzTabStream`, which is not in the crate |
| `void storeProteinReliabilityColumn(bool)` | public field `MzTabFile::store_protein_reliability` |
| `void storePeptideReliabilityColumn(bool)` | public field `MzTabFile::store_peptide_reliability` |
| `void storePSMReliabilityColumn(bool)` | public field `MzTabFile::store_psm_reliability` |
| `void storeSmallMoleculeReliabilityColumn(bool)` | public field `MzTabFile::store_small_molecule_reliability` |
| `void storeProteinUriColumn(bool)` | public field `MzTabFile::store_protein_uri` |
| `void storePeptideUriColumn(bool)` | public field `MzTabFile::store_peptide_uri` |
| `void storePSMUriColumn(bool)` | public field `MzTabFile::store_psm_uri` |
| `void storeSmallMoleculeUriColumn(bool)` | public field `MzTabFile::store_small_molecule_uri` |
| `void storeProteinGoTerms(bool)` | public field `MzTabFile::store_protein_go_terms` |
| `void load(const std::string&, MzTab&)` | `MzTabFile::load`, which returns the document instead of filling an out-parameter |

None of the nine setters validates anything, so the skill's rule makes each a
public field rather than an `x()`/`set_x(v)` pair.

### Protected data members

All sixteen are public fields in Rust, because seven of them have no setter at
all in C++ and are therefore unreachable from outside the class.

| C++ | Rust |
|---|---|
| `store_protein_reliability_` | `store_protein_reliability` |
| `store_peptide_reliability_` | `store_peptide_reliability` |
| `store_psm_reliability_` | `store_psm_reliability` |
| `store_smallmolecule_reliability_` | `store_small_molecule_reliability` |
| `store_protein_uri_` | `store_protein_uri` |
| `store_peptide_uri_` | `store_peptide_uri` |
| `store_psm_uri_` | `store_psm_uri` |
| `store_smallmolecule_uri_` | `store_small_molecule_uri` |
| `store_protein_goterms_` | `store_protein_go_terms` |
| `store_nucleic_acid_reliability_` | `store_nucleic_acid_reliability` (no C++ setter) |
| `store_oligonucleotide_reliability_` | `store_oligonucleotide_reliability` (no C++ setter) |
| `store_osm_reliability_` | `store_osm_reliability` (no C++ setter) |
| `store_nucleic_acid_uri_` | `store_nucleic_acid_uri` (no C++ setter) |
| `store_oligonucleotide_uri_` | `store_oligonucleotide_uri` (no C++ setter) |
| `store_osm_uri_` | `store_osm_uri` (no C++ setter) |
| `store_nucleic_acid_goterms_` | `store_nucleic_acid_go_terms` (no C++ setter) |
| — | `MzTabFile::lossless`, native: all sixteen set, which is what a lossless round trip needs |

### Protected member functions — defined in the `.cpp`

| C++ | Rust |
|---|---|
| `generateMzTabMetaDataSection_(const MzTabMetaData&, StringList&) const` | `MzTabFile::metadata_section_lines`, returning the lines rather than appending to an out-parameter |
| `generateMzTabProteinHeader_(const MzTabProteinSectionRow&, Size, const std::vector<std::string>&, const MzTabMetaData&, size_t&) const` | `MzTabFile::protein_header(&SectionLayout)`; the reference row, the score count, the optional column names and the metadata all fold into `SectionLayout`, and the `n_columns` out-parameter is unnecessary because the layout guarantees the count |
| `generateMzTabSectionRow_(const MzTabProteinSectionRow&, …) const` | `MzTabFile::protein_row` |
| `generateMzTabPeptideHeader_(Size, Size, Size, Size, Size, const std::vector<std::string>&, size_t&) const` | `MzTabFile::peptide_header(&SectionLayout)` |
| `generateMzTabSectionRow_(const MzTabPeptideSectionRow&, …) const` | `MzTabFile::peptide_row` |
| `generateMzTabPSMHeader_(Size, const std::vector<std::string>&, size_t&) const` | `MzTabFile::psm_header(&SectionLayout)` |
| `generateMzTabSectionRow_(const MzTabPSMSectionRow&, …) const` | `MzTabFile::psm_row` |
| `generateMzTabSmallMoleculeHeader_(Size, Size, Size, Size, Size, const std::vector<std::string>&, size_t&) const` | `MzTabFile::small_molecule_header(&SectionLayout)` |
| `generateMzTabSectionRow_(const MzTabSmallMoleculeSectionRow&, …) const` | `MzTabFile::small_molecule_row` |
| `generateMzTabNucleicAcidHeader_(Size, Size, Size, const std::vector<std::string>&, size_t&) const` | `MzTabFile::nucleic_acid_header(&SectionLayout)` |
| `generateMzTabSectionRow_(const MzTabNucleicAcidSectionRow&, …) const` | `MzTabFile::nucleic_acid_row` |
| `generateMzTabOligonucleotideHeader_(Size, Size, Size, const std::vector<std::string>&, size_t&) const` | `MzTabFile::oligonucleotide_header(&SectionLayout)` |
| `generateMzTabSectionRow_(const MzTabOligonucleotideSectionRow&, …) const` | `MzTabFile::oligonucleotide_row` |
| `generateMzTabOSMHeader_(Size, const std::vector<std::string>&, size_t&) const` | `MzTabFile::osm_header(&SectionLayout)` |
| `generateMzTabSectionRow_(const MzTabOSMSectionRow&, …) const` | `MzTabFile::osm_row` |
| `template <typename SectionRow> generateMzTabSection_(…) const` | the per-section loops inside `MzTabFile::document_lines`; its `Exception::Postcondition` becomes the `Error::InvalidValue` that `emit_section` raises when a header and a row disagree on the column count — unreachable by construction, and checked rather than assumed |
| `static void addOptionalColumnsToSectionRow_(const std::vector<std::string>&, const std::vector<MzTabOptionalColumnEntry>&, StringList&)` | `optional_column_cells`, a free function returning the cells |
| `static std::pair<int, int> extractIndexPairsFromBrackets_(const std::string&)` | `extract_index_pairs_from_brackets`, returning `Result<(usize, usize)>` instead of substituting `0` for a missing group |
| file-static `extractBracketIndex(std::string, const std::string&)` in `MzTabFile.cpp:29` | `extract_bracket_index`, returning `Result<usize>` |

### Protected member functions — declared and never defined

These fifteen appear in the header and have no definition anywhere in the pinned
tree; `grep -c` over `MzTabFile.cpp` finds zero occurrences of each, and no
other translation unit defines them. Taking the address of any one, or calling
it, is a link error. None is ported, and none should be: they describe an
identification-export path that `MzTab::IDMzTabStream` replaced. Recorded as
issue **CPP-MZTABFILE-01** below.

| C++ | Rust |
|---|---|
| `static void sortPSM_(PeptideIdentificationList::iterator, PeptideIdentificationList::iterator)` | not ported: declared, never defined |
| `static void keepFirstPSM_(PeptideIdentificationList::iterator, PeptideIdentificationList::iterator)` | not ported: declared, never defined |
| `static void partitionIntoRuns_(…)` | not ported: declared, never defined |
| `static void createProteinToPeptideLinks_(…)` | not ported: declared, never defined |
| `static std::string extractProteinAccession_(const PeptideHit&)` | not ported: declared, never defined |
| `static std::string extractPeptideModifications_(const PeptideHit&)` | not ported: declared, never defined |
| `static std::string mapSearchEngineToCvParam_(const std::string&)` | not ported: declared, never defined |
| `static std::string mapSearchEngineScoreToCvParam_(const std::string&, double, std::string)` | not ported: declared, never defined |
| `static std::string extractNumPeptides_(…)` | not ported: declared, never defined |
| `static std::string extractNumPeptidesDistinct_(…)` | not ported: declared, never defined |
| `static std::string extractNumPeptidesUnambiguous_(…)` | not ported: declared, never defined |
| `static std::map<std::string, Size> extractNumberOfSubSamples_(…)` | not ported: declared, never defined |
| `static void writePeptideHeader_(SVOutStream&, std::map<std::string, Size>)` | not ported: declared, never defined |
| `static void writeProteinHeader_(SVOutStream&, std::map<std::string, Size>)` | not ported: declared, never defined |
| `static void writeProteinData_(…)` | not ported: declared, never defined |

### Private and file scope

| C++ | Rust |
|---|---|
| `friend class MzTabMFile` | not ported: Rust has no friendship. the sibling `mztab_m` module implements `MzTabMFile`; native generators are public without a friendship mechanism |
| `class SVOutStream;` forward declaration | not ported: only the two undefined `write*Header_` statics mention the type |

### Native additions

| Rust | Why |
|---|---|
| `SectionLayout` and its seven `for_*` constructors | The column layout of a section is a value in its own right. The source recomputes it in three mutually inconsistent ways; making it explicit is what removes the four `Exception::Postcondition` throw sites |
| `MzTabFile::load_reader`, `MzTabFile::load_str` | Reading from anything that is `BufRead`, and from a string. The source only reads files |
| `MzTabFile::load_reporting` | Returns the mandatory-column diagnostics the source writes to `std::cout` |
| `MzTabFile::write_to_string`, `MzTabFile::document_lines` | Rendering without a filesystem. The source only writes files |
| `MzTabFile::MAX_LINES`, `MAX_BYTES`, `MAX_COLUMNS`, `MAX_INDEX` | Resource ceilings; the source has none |

---

## Preserved source conventions

- **Line numbering counts every line.** The source increments `line_number` in
  the `for`-statement, so a `continue` still advances it
  (`MzTabFile.cpp:223`). Comment and empty-row positions are therefore indices
  into the whole file, and `MzTab::comment_rows` / `MzTab::empty_rows` keep that
  meaning.
- **A line shorter than three bytes is an empty row, not data.**
  `StringUtils::trim(s).size() < 3` (`MzTabFile.cpp:229`) records the position
  and moves on, so a one- or two-character line is discarded rather than parsed.
- **A data line with fewer than three tab-separated cells is a parse error**
  (`MzTabFile.cpp:249`), with the source's own wording about the tabulator.
- **An unknown three-letter tag is ignored**, and so is an unrecognised `MTD`
  key: the source's if/else chain simply falls through.
- **Optional column names come back in name order, not header order.** The
  source keys its discovery map by name (`std::map<std::string, Size>`,
  `MzTabFile.cpp:109`), so `row.opt_` is filled in lexicographic order whatever
  order the header used. Preserved with a `BTreeMap`, because
  `MzTab::getProteinOptionalColumnNames` then derives the *written* order from
  `opt_`, and changing the read order would change the written order.
- **`addOptionalColumnsToSectionRow_` takes the first match and writes `null`
  for a miss** (`MzTabFile.cpp:2873`), so a duplicated optional column name is
  shadowed by the earlier entry.
- **`MzTabString("null")` is the null cell.** `MzTabString::set` stores nothing
  for the literal text `null` in any case, so the source's "missing optional
  column" cell and a genuinely null one are the same value.
- **Metadata write order**, key for key: `mzTab-version`, `mzTab-mode`,
  `mzTab-type` unconditionally — a null cell writes the text `null` — then
  `title` and `mzTab-ID` only when set, then `description` unconditionally, then
  `sample_processing`, the seven `*_search_engine_score` families,
  `instrument`, `software` (written even when the term is null),
  `false_discovery_rate`, `publication`, `contact`, `uri`, `fixed_mod`,
  `variable_mod`, `quantification_method`, the three `*-quantification_unit`
  keys, `ms_run`, `custom`, `sample`, `assay`, `study_variable`, `cv`, and the
  four `colunit-*` keys.
- **`Complete` mode declares a score column for every `ms_run[n]` of the
  metadata**, whether or not a row carries a value for it
  (`MzTabFile.cpp:3168`). In `Summary` mode only the runs the rows carry are
  declared.
- **Section write order** `PRT`, `PEP`, `PSM`, `SML`, `NUC`, `OLI`, `OSM`, each
  preceded by a blank line and introduced by its header row; an empty section
  contributes nothing.
- **The `PSM` and `OSM` sections have a flat `search_engine_score[n]` family**
  while the other five have the two-dimensional
  `search_engine_score[s]_ms_run[r]`.
- **`MzTabDouble` renders with fifteen fractional digits, trimmed**, so
  `51.9678841193106` is written `51.967884119310597`; and `|value| >= 1e4` or
  `< 1e-2` switches to scientific notation, so `5035500000` is written
  `5.0355e09`. Both come from the data model's `format_float`, which is
  `StringUtils::toStr(double)`.

---

## Native differences

Each of these is a place where the source's reader and its writer disagree, or
where the source cannot read back what it writes. The line numbers are
`MzTabFile.cpp` at the pinned revision. The round trip the task asks for is not
reachable without them.

1. **`search_engine_score[s]_ms_run[r]` column order.** The header writes them
   runs-outer, scores-inner (`:2037`, `:2231`, `:2472`, `:2588`, `:2718`); every
   row writer emits the values scores-outer, runs-inner (`:2119`, `:2334`,
   `:2540`, `:2656`, `:2767`). The counts agree, so the postcondition stays
   silent, and with two or more score types *and* two or more runs the values
   land under the wrong headers. This port uses scores-outer, runs-inner on both
   sides — which agrees with the source on every reference file, all of which
   have exactly one score type. Test:
   `two_score_types_over_two_runs_keep_header_and_values_aligned`.

2. **Header and row column counts are made equal by construction.** The source
   derives its counts from the metadata in the protein section, from the *first*
   row in the peptide and small-molecule sections, and from each row's own maps
   in the row writers; where they disagree it throws
   `Exception::Postcondition` from `generateMzTabSection_` (`:139`) or one of
   the three copies in `store` (`:2958`, `:2987`, `:3086`). Here one
   `SectionLayout` per section is the union of every row's keys with the
   metadata keys the source consults, and each row is rendered against it with a
   `null` cell for a member it lacks.

3. **The small-molecule row writer emits no assay cells** (`:2548`), although
   its header declares one `smallmolecule_abundance_assay[n]` per assay
   (`:2482`), so any document with at least one assay cannot be written at all.
   This port emits them. Test: `small_molecule_assay_cells_are_written`.

4. **The peptide and small-molecule study-variable triples are driven by the
   layout, not by three iterators in lock-step** (`:2367`, `:2555`). The source
   stops at the first exhausted map, so a row with values but no standard
   deviations writes no study-variable columns at all. Test:
   `peptide_study_variable_triple_survives_a_missing_stdev_map`.

5. **`num_osms_ms_run`, `num_oligos_distinct_ms_run` and
   `num_oligos_unique_ms_run` are numbered from one**, from the keys the rows
   carry. The source's nucleic-acid header numbers them from *zero* (`:2603`),
   alone among all bracketed families, while its row writer emits the values
   keyed from one. Test:
   `nucleic_acid_count_columns_are_numbered_from_the_data`.

6. **`generateMzTabNucleicAcidHeader_` is called with its score and best-score
   counts transposed** at `:3257`: the signature is `(search_ms_runs,
   n_best_search_engine_scores, n_search_engine_scores, …)` and the call passes
   `(search_ms_runs, n_search_engine_score, n_best_search_engine_score, …)`.
   This port does not transpose them.

7. **The `PEH` and `PSH` `uri` columns are read into their own row.** The source
   assigns both column indices to `protein_uri_index` (`:1089`, `:1265`), a
   copy-paste, and then never reads the peptide or PSM `uri` cell at all (both
   assignments are commented out, `:1160`, `:1324`). Worse, when a `PEH` or
   `PSH` follows a `PRH` in the same file, the corrupted index makes a later
   `PRT` row read its protein `uri` out of the peptide section's column number.
   Test: `uri_columns_are_read_into_their_own_section`.

8. **`PSM` optional columns are discovered by the `opt_` prefix.** The source
   tests `cells[i] == "opt_"` for exact equality (`:1287`) where the other three
   sections test the prefix, so no `PSM` section's optional column is ever read
   and every `opt_` cell of a PSM row is silently dropped. Test:
   `psm_optional_columns_are_discovered_from_the_header`.

9. **`colunit-*` is written with a tab and read without an index.** The source
   concatenates the key and the value with no separator at all (`:1984`), so its
   own reader sees a two-cell line; and even a correctly separated line makes
   its reader throw, because it hands the literal key `colunit` to
   `StringUtils::toInt32` looking for a bracketed index a `colunit` key never
   carries (`:685`). Any `MTD colunit-…` line therefore aborts
   `MzTabFile::load` with a `ConversionError`. The source also writes
   `colunit-PSM` while its reader matches the lower-case `psm`; this port writes
   the specification's `colunit-psm` and accepts either case. Test:
   `colunit_keys_round_trip`.

10. **`nucleic_acid_search_engine_score[n]`,
    `oligonucleotide_search_engine_score[n]` and `osm_search_engine_score[n]`
    are read.** The source writes all three (`:1608`-`:1624`) and reads none of
    them, and it does not read the `NUH`/`NUC`, `OLH`/`OLI` or `OSH`/`OSM`
    sections at all — `load` handles only `MTD`, `COM`, `PRH`/`PRT`, `PEH`/`PEP`,
    `PSH`/`PSM` and `SMH`/`SML`. This port reads all seven sections and all seven
    score families, which is what makes the extension sections round-trip. Test:
    `every_section_type_and_metadata_key_round_trips`.

11. **A column the header did not declare leaves its field at the default.** The
    source uses `0` as the "column absent" sentinel, which is also the index of
    the section tag, and then reads `cells[0]` into the field: a `PRT` row whose
    header omitted `taxid` parses the text `PRT` as an integer and throws. Test:
    `optional_columns_tolerate_any_header_order_and_unknown_names`.

12. **A row shorter than its header leaves the missing cells unset** and adds a
    diagnostic. The source indexes `cells[index]` unconditionally, which is an
    out-of-range `std::vector::operator[]`. Test:
    `a_row_shorter_than_its_header_is_reported_and_leaves_cells_unset`.

13. **A bracketed index must be a positive integer no greater than
    `MzTabFile::MAX_INDEX`.** The source's `extractBracketIndex` returns a signed
    `Int` that every caller casts to `Size`, so `assay[0]` becomes key `0` and
    `assay[-1]` becomes key `18446744073709551615`, both of which the writer
    then emits as a column name the format does not define; and
    `extractIndexPairsFromBrackets_` substitutes `0` for a bracketed group its
    regex does not match. Tests:
    `bracket_index_extraction_follows_the_source_and_refuses_zero`,
    `index_pair_extraction_needs_two_bracketed_groups`,
    `a_bad_index_in_a_metadata_key_is_a_parse_error`,
    `a_bad_index_in_a_section_header_is_a_parse_error`.

14. **A key with fewer `-` fields than a branch inspects fails to match.** The
    source indexes `meta_key_fields[1]` in thirty-two branches and
    `meta_key_fields[2]` in two more without checking the field count, so a key
    such as `MTD instrument[1] …` — no suffix — reads past the end of the
    vector; and it reads field zero of an empty vector when the key cell itself
    is empty.

15. **A section header replaces the previous one for that section.** The source
    accumulates into the same column maps, so a file with two `PRH` lines keeps
    the first one's column numbers alongside the second's. Test:
    `a_second_section_header_replaces_the_first`.

16. **Blank lines are restored without being doubled.** The source walks the
    generated lines and *inserts* a blank whenever the current position was
    blank in the original file (`FORMAT/MzTabFile.cpp:3328`), even though the
    generated lines already
    carry a blank before every section, so every separator is emitted twice. Its
    own round-trip test cannot see this, because
    `FuzzyStringComparator::readNextLine_` skips blank lines outright
    (`FuzzyStringComparator.cpp:843`). Here a recorded blank position consumes a
    generated blank when there is one. The source also stops as soon as the
    generated lines run out, dropping any comment or blank recorded past that
    point; those are emitted here. Their recorded positions only survive when
    the recorded tail is contiguous with the end of the generated lines —
    nothing fills the gap left by a recorded metadata key that is not written
    back, so a tail behind such a gap keeps its order and its content but moves
    up by the width of the gap. Tests:
    `store_restores_comments_and_blank_lines_in_place`,
    `store_emits_a_comment_and_a_blank_recorded_past_the_generated_lines`.

17. **A cell carrying a tab or a line break is refused.** Either would corrupt
    the row; the source concatenates unconditionally and produces a file it
    cannot read back. Test: `a_failed_store_leaves_the_previous_file_intact`.

18. **The `psm_search_engine_score` count is not capped at one.** The source
    caps it with the comment "we currently only store one search engine score
    per PSM" (`:3200`); that limit belongs to its `ConsensusMap` and
    identification exporters, not to the file format.

19. **`store` renders and validates the whole document before creating
    anything**, then publishes the bytes by renaming a sibling temporary file.
    The source opens the destination with `ios::trunc` first and writes as it
    goes, so a failure halfway leaves a truncated file where a good one stood.
    Test: `a_failed_store_leaves_the_previous_file_intact`.

20. **`load` decompresses a gzip or bzip2 input transparently**, through the
    crate's shared `path_io::open`. The source reads plain text only.

21. **Diagnostics are returned, not printed.** The twenty mandatory-column
    checks of `:832`-`:934` write to `std::cout` once per data row — a
    thousand-row section prints the same line a thousand times — and the caller
    never sees them. `MzTabFile::load_reporting` returns them, once per section,
    capped at 1024. Test:
    `load_reporting_returns_the_mandatory_column_diagnostics`.

22. **`sections_present` is not ported.** The source fills it while parsing
    (and inserts `"PRT"` from the `PEP` branch, `:1130`, another copy-paste) but
    its only consumer, `hasMandatoryMetaDataKeys_`, is commented out at `:1547`.

### One thing that is lost, on purpose

The `reliability`, `uri` and `go_terms` columns are optional in the MzTab
specification, and the sixteen flags decide whether they are written. With a
flag cleared, a row that carries the corresponding cell does not round-trip —
the cell is simply not in the file. `MzTabFile::lossless` sets all sixteen.
Test: `clearing_an_optional_column_flag_drops_its_cells`.

### What the format itself cannot express

A section's header declares one column set for every row. A member of an indexed
family that only some rows carry therefore becomes a `null` cell — and so a null
entry — in the rest: the file cannot distinguish "this row has no
`protein_abundance_assay[2]`" from "this row's `protein_abundance_assay[2]` is
null". A document built in memory is canonicalised in exactly that way by its
first write, and the round trip is a fixed point from then on. The blank
separator line before each section is recorded by the reader for the same
reason. Both are asserted item by item in
`the_canonicalisation_a_write_performs_is_exactly_the_declared_columns`.

One more asymmetry, shared with the source and not fixed here: an indexed
metadata group survives a round trip only if at least one of its keys has a
non-null value, because every optional metadata key is skipped when null. A
document whose `assay[1]` is declared only by `MTD assay[1]-sample_ref null`
loses the key, and with it the `<section>_abundance_assay[1]` column the layout
would have declared. Closing it would need a record of which keys were
*declared* as distinct from which carry a value, which neither the source nor
`MzTabMetaData` has. No reference document is affected: in all five, every
`assay[n]` and `study_variable[n]` key is introduced by a line with a value.

### Parallelism

`MzTabFile.cpp` carries no `#pragma omp`, so there is no OpenMP gap to record
here. Reading and writing are both single-pass and serial in the source and in
the port.

---

## Checked boundaries and evidence

### Resource ceilings

| Constant | Value | Checked before |
|---|---|---|
| `MzTabFile::MAX_LINES` | 4,000,000 | reading or writing another line; also caps the recorded empty-row list |
| `MzTabFile::MAX_BYTES` | 512 MiB | accumulating another line, in both directions |
| `MzTabFile::MAX_COLUMNS` | 200,000 | splitting a line, parsing a header, computing a `SectionLayout`, and splitting a comma-separated reference list — each by counting separators first, so nothing is collected before the ceiling is charged |
| `MzTabFile::MAX_INDEX` | 1,000,000 | inserting an indexed metadata key or column family member |
| `MzTab::MAX_ROWS` | from the data model | appending a row to a section |
| `MzTab::MAX_OPTIONAL_COLUMNS` | from the data model | recording an optional column, in either direction |
| `TextFile::Limits::max_line_bytes` | 16 MiB | reading one line, through `TextFile::get_line_with_limits` |

Every per-cell ceiling — `MAX_CELL_BYTES`, `MAX_CELL_ITEMS` — belongs to the
data model and is enforced by each cell's `read_cell`.

A `SectionLayout` is computed and range-checked before a single row is rendered,
and `store` renders the whole document into memory before the destination is
touched, so a refusal leaves the input document and the previous output
unchanged.

No `unsafe`, no threads, no panicking index or slice on file-derived data: every
cell access is `slice::get`, the three-byte section tag is `str::get(..3)`
(which yields `None` rather than panicking when byte three falls inside a
character), and every arithmetic step on a counter is `checked_add`.

### Evidence

**Tier 3, source review.** `MzTabFile_test.cpp` has five `START_SECTION`s and
none has more than four assertion macros; all five are ported, not mapped:

| Section | Assertion macros | Rust test |
|---|---|---|
| `MzTabFile()` | 1 (`TEST_NOT_EQUAL`) | `default_construction_disables_every_optional_column` |
| `void load(const std::string&, MzTab&)` | 0 | `load_reads_the_silac_reference_document`, `load_reads_every_reference_document`, `load_reads_the_cytidine_small_molecule_document` |
| `void store(const std::string&, MzTab&)` | 1 (`TEST_FILE_SIMILAR`, over five files) | `store_reproduces_every_reference_document`, `store_round_trips_every_reference_document_through_a_file`, `store_restores_comments_and_blank_lines_in_place` |
| `~MzTabFile()` | 0 | `dropping_the_adapter_releases_nothing` |
| `generateMzTabPSMSectionRow_` | 4 (`TEST_EQUAL`) | `psm_row_fills_requested_optional_columns_and_nulls_the_rest` |

The four transcribed literals of the last section are the last four cells of a
PSM row — `null`, `NDYKAPPQPAPGK`, `0.0420992`, `null` — and the Rust test
reproduces the whole row it renders them from, including the
`51.967884119310597` spelling of `51.9678841193106`.

The store section's comparison is reproduced exactly rather than approximated:
`similar_lines` sorts the lines, removes every space, drops the blank ones as
`FuzzyStringComparator::readNextLine_` does, and compares numbers numerically
wherever a numeric literal begins on both sides — which is how `46` matches
`46.0` and `5035500000` matches `5.0355e09`. The numeric comparison here is
*stricter* than upstream's: `TEST_FILE_SIMILAR` never uses
`FuzzyStringComparator`'s constructor defaults, because `TEST::isFileSimilar`
overrides them with `absdiff_max_allowed` = 1e-5 and `ratio_max_allowed` =
1 + 1e-5 (`ClassTest.cpp:591-592`, values at `ClassTest.cpp:35, 38`, and
`MzTabFile_test.cpp` sets no `TOLERANCE_*`), whereas `similar_line` requires
exact `f64` equality. Nothing upstream tolerates is refused here that the
reference files actually contain, but the comparison is not tolerance-faithful.

`store_round_trips_every_reference_document_through_a_file` is the stronger
statement, and it is the one that had to be earned: the five reference documents
are written to disk, read back, and compared as *models*, which covers the
`null`/`NaN`/`Inf` state of every numeric cell, the optional columns of every
row, and the recorded comment and blank-line positions — none of which a sorted,
space-stripped text comparison can see.

**Tier 4, independently derived.** The resource ceilings, the hostile headers
and metadata keys, every error-variant choice, the non-ASCII fixture and output
paths, and every test of a divergence listed above. No C++ was built or executed
and no C++ output was retained, so this is **not** a tier 1 differential.

### Fixtures

The five upstream reference documents are copied unmodified into `tests/data/`,
with their sha256 in the provenance manifest.
`tests/data/MzTabFile_unicode.mzTab` is new: a Japanese `mzTab-ID`, description
and comment, a URL with Japanese path segments, a `colunit-psm` key, an
`opt_global_備考` column, a sequence of two-byte characters, and the three
numeric states. Every value in it has an exact binary representation and its
metadata is already in the writer's order, so
`non_ascii_content_round_trips_byte_for_byte` asserts byte equality — the
tightest check in the package.

---

## OpenMS C++ issues found

Proposed for `OpenMS_CPP_ISSUES.md`; the integrating agent owns that file.

- **CPP-MZTABFILE-01** — fifteen protected static members of `MzTabFile.h`
  (lines 151-201) are declared and defined nowhere in the tree. Dead
  declarations; any use is a link error. Fix: delete them, and the
  `class SVOutStream;` forward declaration that only they need.
- **CPP-MZTABFILE-02** — `colunit-*` cannot round-trip: written without a
  separator (`:1984`), and read by feeding the literal key `colunit` to
  `StringUtils::toInt32` (`:685`), which throws. Any `MTD colunit-…` line aborts
  `load`. Fix: write `"MTD\tcolunit-protein\t" + value`; read the value with no
  index, into a `std::vector<std::string>`; and write `colunit-psm` to match the
  reader and the specification.
- **CPP-MZTABFILE-03** — `search_engine_score[s]_ms_run[r]` header and row
  orders are transposed relative to each other (`:2037` vs `:2119`, and the
  same in four more sections). With ≥2 score types and ≥2 runs the values land
  under the wrong headers, and the column-count postcondition does not notice.
  Fix: use one nesting order in both.
- **CPP-MZTABFILE-04** — the small-molecule row writer never emits
  `smallmolecule_abundance_assay[n]` (`:2548`) although the header declares one
  per assay (`:2482`), so any document with an assay throws
  `Exception::Postcondition`. Fix: emit the cells, as the peptide writer does.
- **CPP-MZTABFILE-05** — the `PEH` and `PSH` `uri` column index is stored in
  `protein_uri_index` (`:1089`, `:1265`). The peptide and PSM `uri` cells are
  never read, and a `PRT` row following a `PEH`/`PSH` section reads its protein
  `uri` from the wrong column. Fix: separate indices, and re-enable the two
  commented-out assignments.
- **CPP-MZTABFILE-06** — `PSH` optional columns are matched with
  `cells[i] == "opt_"` instead of the `opt_` prefix (`:1287`), so every
  optional column of a PSM section is dropped on load. Fix:
  `StringUtils::hasPrefix(cells[i], "opt_")`.
- **CPP-MZTABFILE-07** — `load` never parses the `NUH`/`NUC`, `OLH`/`OLI` or
  `OSH`/`OSM` sections, nor the three `*_search_engine_score` metadata families
  of those sections, all of which `store` writes. Every nucleic-acid,
  oligonucleotide and OSM row is lost by a load/store cycle. Fix: add the
  reader branches.
- **CPP-MZTABFILE-08** — the nucleic-acid header numbers `num_osms_ms_run`,
  `num_oligos_distinct_ms_run` and `num_oligos_unique_ms_run` from zero
  (`:2603`) while the row writer emits values keyed from one. Fix: `i + 1`, as
  every other family does.
- **CPP-MZTABFILE-09** — `generateMzTabNucleicAcidHeader_` is called with its
  second and third arguments swapped (`:3257`). Fix: swap them back.
- **CPP-MZTABFILE-10** — blank lines are duplicated by the comment restoration
  (`:3328`): a recorded blank is inserted even when the generated line at that
  position is already the section separator. Masked by
  `FuzzyStringComparator`, which skips blank lines. Fix: consume the generated
  blank when there is one.
- **CPP-MZTABFILE-11** — out-of-range reads. `meta_key_fields[1]` is indexed in
  thirty-two `MTD` branches and `meta_key_fields[2]` in two more without
  checking the field count; `meta_key_fields[0]` is indexed when
  `StringUtils::split` cleared the vector for an empty key cell; and every
  section row reads `cells[index]` with no bounds check, so a row shorter than
  its header reads past the end. Fix: check the sizes.
- **CPP-MZTABFILE-12** — `extractBracketIndex` returns a signed `Int` that every
  caller casts to `Size`, so `assay[0]` yields key `0` and `assay[-1]` yields
  `SIZE_MAX`, which the writer emits as a column name the format does not
  define. `extractIndexPairsFromBrackets_` likewise substitutes `0` for a
  bracketed group its regex does not match. Fix: reject a non-positive index.
- **CPP-MZTABFILE-13** — the `PEP` branch inserts `"PRT"` into
  `sections_present` (`:1130`). No effect today, because the only consumer is
  commented out at `:1547`. Fix: insert `"PEP"`, or delete the set.
- **CPP-MZTABFILE-14 — withdrawn**: the earlier character-narrowing claim
  overlooked `StringUtils.h`'s numeric string operators. Source and native
  modification positions render as decimal digits; this is not a source bug.
