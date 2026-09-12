# MSstatsFile support

`FORMAT/MSstatsFile.h` and `FORMAT/MSstatsFile.cpp` at `bc9cc12`, ported to
`src/format/msstats.rs` and tested by `tests/msstats.rs`. Fixtures, hashes,
the upstream test definitions and source anchors are in
`tests/data/msstats_provenance.json`.

MSstats consumes a long-format CSV: one row per quantified peptide ion per run,
carrying the protein, the modified peptide sequence, the precursor charge, the
experimental annotation taken from the experimental design, the intensity and a
reference string naming where the intensity came from. Two layouts exist —
label-free (`storeLFQ`) and isobaric/MSstatsTMT (`storeISO`). The module has no
Cargo feature gate; `tests/msstats.rs` is gated on `consensusxml`, which is
what its fixtures need to be read.

MSstats counts a run as one `(spectra file, fraction)` pair, while OpenMS splits
one run into fractions; that translation is `assemble_run_map`.

## API mapping

Every member of the header, public and private, plus the two namespace-level
aliases it declares.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `using IndProtGrp = ProteinIdentification::ProteinGroup` | `crate::identification::ProteinGroup` | The crate already has the type; no alias is added. |
| `using IndProtGrps = std::vector<IndProtGrp>` | `&[ProteinGroup]` | |
| `MSstatsFile() = default` | not ported: the class holds no state, so there is nothing to construct. The two operations are free functions. |
| `~MSstatsFile() = default` | not ported: nothing to release. |
| `void storeLFQ(filename, consensus_map, design, reannotate_filenames, is_isotope_label_type, bioreplicate, condition, retention_time_summarization_method, remove_shared_peptides = true)` | `store_lfq`, `write_lfq`, `prepare_lfq` | The last five arguments become `LfqOptions`. `prepare_lfq` returns the rows without writing; `write_lfq` writes to a stream; `store_lfq` publishes a file only after the bytes are complete. All three return `MSstatsReport`. |
| `void storeISO(filename, consensus_map, design, reannotate_filenames, bioreplicate, condition, mixture, retention_time_summarization_method, remove_shared_peptides = true)` | `store_iso`, `write_iso`, `prepare_iso` | Same shape, with `IsoOptions`. |
| `typedef Peak2D::IntensityType Intensity` | `f32` | |
| `typedef Peak2D::CoordinateType Coordinate` | `f64` | |
| `static const std::string na_string_` | `NA` | Public here: it is part of the written format. |
| `static const char delim_ = ','` | `DELIMITER` | Public, same reason. |
| `static const char accdelim_ = ';'` | `ACCESSION_DELIMITER` | Public, same reason. |
| `static const char quote_ = '"'` | `QUOTE` | Public, same reason. |
| `struct AggregatedConsensusInfo` | private `Aggregated` | |
| `AggregatedConsensusInfo::consensus_feature_filenames` | `Aggregated::filenames` | |
| `AggregatedConsensusInfo::consensus_feature_intensities` | `Aggregated::intensities` | |
| `AggregatedConsensusInfo::consensus_feature_retention_times` | `Aggregated::retention_times` | |
| `AggregatedConsensusInfo::consensus_feature_labels` | `Aggregated::labels` | |
| `AggregatedConsensusInfo::features` | not ported: the source's third copy of the features, which its own `//TODO` at `MSstatsFile.cpp:78` questions ("why do we need this method and store everything three times"). The loop indexes `map.features` directly. |
| `AggregatedConsensusInfo aggregateInfo_(const ConsensusMap&, const std::vector<std::string>&)` | private `aggregate_info` | |
| `static void checkConditionLFQ_(const SampleSection&, const std::string&, const std::string&)` | `check_condition_lfq` | Public here: a tool wants to validate a design before it starts converting. Returns the paired-design warnings instead of logging them. |
| `static void checkConditionISO_(const SampleSection&, const std::string&, const std::string&, const std::string&)` | `check_condition_iso` | Public, same reason. |
| `static void assembleRunMap_(std::map<std::pair<std::string, unsigned>, unsigned>&, const ExperimentalDesign&)` | `assemble_run_map` | Public here: the MSstats-run numbering is a documented translation a caller may need to report. The out-parameter becomes the return value. |
| `static bool isSubsetOf_(const std::vector<std::string>&, const std::vector<std::string>&)` | private `is_subset_of` | |
| `static void warnOnSubsetFiles_(const std::vector<std::string>&, const std::vector<std::string>&)` | the subset warning inside `prepare_common` | Text preserved; returned in `MSstatsReport::warnings` rather than logged. |
| `Peak2D::IntensityType sumIntensity_(const std::set<IntensityType>&) const` | `RetentionTimeSummarization::Sum` in `combine` | |
| `Peak2D::IntensityType meanIntensity_(const std::set<IntensityType>&) const` | `RetentionTimeSummarization::Mean` in `combine` | |
| `class MSstatsLine_` | private `LineKey` plus the rendered prefix string | The class is private in the source too, and its only observable products are its ordering and its `toString()`. Both are reproduced; the type is not. |
| `MSstatsLine_::MSstatsLine_(has_fraction, accession, sequence, precursor_charge, fragment_ion, frag_charge, isotope_label_type, condition, bioreplicate, run, fraction)` | the `prefix` built in `prepare_lfq` | |
| `MSstatsLine_::accession()`, `sequence()`, `precursor_charge()`, `run()` | fields of `LineKey` | The four accessors exist only to build the sanity-check triple, which is dead (see below). |
| `MSstatsLine_::toString()` | the `prefix` format string in `prepare_lfq` | |
| `operator<(const MSstatsLine_&, const MSstatsLine_&)` | `LineKey`'s derived `Ord` over `[accession, run, condition, bioreplicate, precursor_charge, sequence]` | Same fields, same order, same text comparison of the charge. |
| `MSstatsLine_`'s eleven private fields | the prefix's eleven rendered columns | |
| `class MSstatsTMTLine_` | private `LineKey` plus the rendered prefix string in `prepare_iso` | |
| `MSstatsTMTLine_::MSstatsTMTLine_(accession, sequence, precursor_charge, channel, condition, bioreplicate, run, mixture, techrepmixture, fraction)` | the `prefix` built in `prepare_iso` | |
| `MSstatsTMTLine_::accession()`, `sequence()`, `precursor_charge()`, `run()` | fields of `LineKey` | |
| `MSstatsTMTLine_::toString()` | the `prefix` format string in `prepare_iso` | |
| `operator<(const MSstatsTMTLine_&, const MSstatsTMTLine_&)` | `LineKey`'s derived `Ord` over `[accession, run, condition, bioreplicate, mixture, precursor_charge, sequence, channel]` | |
| `MSstatsTMTLine_`'s ten private fields | the prefix's ten rendered columns | |
| `template <class LineType> void constructFile_(retention_time_summarization_method, rt_summarization_manual, TextFile&, const std::set<std::string>&, LineType&) const` | private `Rows::render` | The template parameter is not needed: both layouts reduce to an ordered key plus a rendered prefix. |
| `static std::unordered_map<std::string, const IndProtGrp*> getAccessionToGroupMap_(const IndProtGrps&)` | private `accession_to_group` | Returns a group index rather than a pointer, so the map borrows nothing mutable. |
| `bool isQuantifyable_(const std::set<std::string>&, const std::unordered_map<std::string, const IndProtGrp*>&) const` | private `is_quantifyable` | |
| (none) | `RetentionTimeSummarization` | The free-form method string becomes an enum with `from_name` and `name`. |
| (none) | `MSstatsReport` | The rows, the warnings the source logged, the shared-peptide drop count and the run-to-fraction-group map the source printed to standard output. |
| (none) | `LfqOptions`, `IsoOptions` | Grouped arguments plus explicit ceilings. |
| (none) | `MAX_FEATURES`, `MAX_LINES`, `MAX_BYTES` | The default ceilings. |
| (none) | private `spectra_paths`, `prepare_common`, `lookup`, `Limits`, `Rows`, `Sample`, `LineSamples` | The shared prologue both layouts run, factored once. |

## Preserved source conventions

- **The two headers, verbatim.** Label-free:
  `ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run,Intensity,Reference`,
  with `RetentionTime,` prefixed when the summarization is manual and
  `Fraction,` inserted before `Intensity` when the design is fractionated.
  Isobaric:
  `RetentionTime,ProteinName,PeptideSequence,Charge,Channel,Condition,BioReplicate,Run,Mixture,TechRepMixture,Fraction,Intensity,Reference`,
  always with the leading `RetentionTime` because the isobaric layout is always
  manual.
- **`FragmentIon` is `NA` and `ProductCharge` is `0`**, because neither is used
  for DDA data.
- **`IsotopeLabelType`** is `L` for endogenous peptides and `H` when
  `is_isotope_label_type` is set; the source carries a `@todo` doubting DDA
  label-free is ever `H`, and that doubt is repeated at the field.
- **An empty accession set becomes `NA`**, which the source notes should not
  matter because such peptides are unquantifiable anyway.
- **Decoy hits are skipped.**
- **The sequence is the modified `toString()` form**, so
  `.(Carbamidomethyl)AETAAQDVQQK` and not `AETAAQDVQQK`.
- **Shared peptides are dropped by default**, with the source's
  "Use -remove_shared_peptides false to keep them" warning and the drop count.
- **`isQuantifyable_`'s rules.** No accession → not quantifiable; exactly one →
  quantifiable; several → only when every accession resolves to the same
  indistinguishable group, with an accession outside every group assumed to be
  a singleton and therefore disqualifying.
- **The label-free precursor charge is the hit's charge unchanged**, while the
  isobaric one is `max(charge, 0)` — MSstats user manual 3.7.3 documents an
  unknown precursor charge as `0`. The asymmetry is the source's.
- **The channel is `channel_id + 1`.** `aggregateInfo_` pushes the raw
  `channel_id` meta value of the column header, or `1` when it has none, and
  `storeISO` adds one. A map whose headers carry no channel id therefore yields
  channel `2`, and `storeLFQ` uses the raw value as the design label, so a
  label-free map only resolves when that value is `1`.
- **`TechRepMixture` is `<mixture>_<fraction group>`** and `Run` is
  `<TechRepMixture>_<fraction>`.
- **The isobaric reference is `<run basename>_<native id>`**, with the native
  id taken from the `spectrum_reference` meta value or the literal
  `NONATIVEID`. The label-free reference is the run basename alone.
- **The reference is quoted with `"`, unescaped.**
- **`FileFilter` gaps.** The run-path vector is sized by the largest column
  index and the run paths are consumed in column order, so a column past the
  last path keeps an empty name.
- **Subset designs warn and proceed.** The active column basenames must be a
  subset of the design's basenames; a strict subset only warns, naming the
  missing files.
- **`storeLFQ` refuses a multi-label design** ("Too many labels for a
  label-free quantitation experiments").
- **An empty protein identification list is refused**; more than one run, and a
  first run without inference data, only warn.
- **A recurring biological replicate warns**, because it legitimately encodes a
  paired design (upstream #9864 B7), and the warning text is preserved.
- **`storeISO` forces manual summarization**, with the source's "Reverting to
  'manual'" warning, because MSstatsTMT aggregates itself.
- **Manual summarization writes every stored sample**, including the duplicate
  retention times the deduplication loop only warned about. Preserved exactly:
  the warning loop and the writing loop read different collections.
- **Aggregation is over *distinct* intensities.** The source collects them into
  a `std::set`, so `sum` ignores repetitions and `mean` divides by the number
  of distinct values. The arithmetic is single precision, and the accumulation
  order is ascending intensity — all three reproduced.
- **The aggregated reference is the *lowest* sample.** The source takes
  `get<2>(*line.second.begin())`, i.e. the first element of the set ordered by
  `(intensity, retention time, reference)`, not the first sample inserted.
- **Ordering.** The output is sorted by the modified sequence (a `std::set` of
  strings, so byte-lexicographic), then by the line key, whose numeric fields
  are compared as text: `Run` `10` sorts before `Run` `2`.
- **A line that compares equal to one already present keeps the first
  rendering**, as the source's `std::map::operator[]` insertion does.

## Native differences

- **A missing `(file, label)` design row is refused.** The source indexes
  `path_label_to_sample`, `path_label_to_fraction`,
  `path_label_to_fractiongroup` and `run_map` with `std::map::operator[]`,
  which inserts a zero for a key the design does not hold; the sample index `0`
  is then read as sample-section row `0`, so the row is silently annotated with
  another sample's condition and biological replicate. This port returns
  `Error::MissingInformation` naming the file and the label.
- **An unknown summarization method is refused.** The source's accumulator is
  initialised to `0` and no branch assigns it unless the name is `max`, `min`,
  `mean` or `sum`, so any other name silently writes the intensity `0` for
  every row. `RetentionTimeSummarization::from_name` returns
  `Error::InvalidValue`.
- **Warnings and the run mapping are returned, not printed.** The source writes
  its warnings to `OPENMS_LOG_WARN` and the MSstats-run-to-fraction-group
  mapping to `std::cout`. Both are fields of `MSstatsReport`. The isobaric
  mapping is empty, reproducing the source's own dead code.
- **The triple sanity check is not ported.** `constructFile_` builds a
  `set<tuple<sequence, precursor_charge, run>>` whose comment says it is a
  "sanity check that the triples ... only appears once", inserts into it and
  never reads it. There is nothing to reproduce.
- **`AggregatedConsensusInfo::features` is not ported.** It is the third copy
  of the features the source's own `//TODO` questions; the loop indexes
  `map.features`.
- **`getAccessionToGroupMap_` returns indices, not pointers.** Group identity
  is compared by index, which is equivalent and borrows nothing.
- **An empty accession is excluded from the joined accession string**, because
  the crate's `PeptideHit::protein_accessions` filters empty accessions out of
  the set. `extractProteinAccessionsSet` does not; a hit carrying an evidence
  with an empty accession would therefore produce a leading `;` in C++.
- **The intensity text is the *current* `StringUtils::toStr(float)`.** The
  retained reference outputs were written before
  `Internal::NumericFormatting::appendNumeric` moved to shortest-round-trip
  scientific notation above 1e4, so their `Intensity` fields are in a fixed
  six-decimal form neither today's C++ nor this port emits. The numbers are the
  same; see the tolerance note below.
- **Ceilings exist.** The source has none.
- **No OpenMP gap.** Neither the header nor its implementation carries a
  `#pragma omp`.

## Checked boundaries and evidence

| Ceiling | Default | What it bounds |
|---|---|---|
| `max_features` | 5,000,000 | consensus features processed, and the largest column index accepted |
| `max_lines` | 20,000,000 | output rows |
| `max_bytes` | 512 MiB | cumulative owned payload |

A zero ceiling is `Error::InvalidValue`, checked before the consensus map is
touched. Nothing is byte-sliced: basenames come from
`crate::system::file::basename`, which splits on `/` and `\`, and every numeric
field is produced by the crate's formatters rather than parsed.

**Evidence: tier 1, executed differential from retained C++ output.** The class
test cannot serve as an oracle — both of its sections have an empty body
carrying only the comment "tested via MSstatsConverter tool", so there is no
assertion and no literal to transcribe. Instead all three retained
`MSstatsConverter` outputs are reproduced row for row and field for field:

| Upstream test | Retained output | Rows | Rust test |
|---|---|---|---|
| `TOPP_MSstatsConverter_1` (`CMakeLists.txt:1642`, diffed at `:1643`) | `MSstatsConverter_1_out.csv` | 2029 | `lfq_reproduces_the_retained_cpp_output` |
| `TOPP_MSstatsConverter_2` (`:1646`, diffed at `:1647`) | `MSstatsConverter_2_out.csv` | 471 | `iso_reproduces_the_retained_cpp_output` |
| `TOPP_MSstatsConverter_3` (`:1650`, diffed at `:1651`) | `MSstatsConverter_3_out.csv` | 471 | `iso_reproduces_the_retained_cpp_output_for_a_second_design` |

The two isobaric cases share one input and differ only in the design — one
fraction group with two fractions versus two fraction groups with one fraction
each — and their retained outputs differ in the `Run`, `TechRepMixture` and
`Fraction` columns. Both pass, so the comparison cannot be satisfied by a
producer that ignores the design.

**Comparison policy.** Byte equality is not the upstream contract: the upstream
`${DIFF}` is `FuzzyDiff` configured (`FuzzyDiff.ini`) with a 1% relative or
0.01 absolute tolerance, and the retained files prove some tolerance is
necessary because their single-precision intensity text predates the current
formatter. `tests/msstats.rs` compares field by field: every field exactly
unless both sides parse as a number, in which case the relative difference must
stay below 1e-6 with a 1e-9 absolute floor. That is far tighter than upstream,
and 1e-6 is the decimal resolution of the `f32` intensities. Row count, field
count per row, and every non-numeric field — protein, sequence, charge,
channel, condition, biological replicate, run, mixture, technical replicate,
fraction and the quoted reference — are exact.

The unknown-summarization refusal, the distinct-intensity aggregation, the
design-column checks, the paired-design warning, the multi-label refusal and
the resource ceilings are independently derived (tier 4).

### Section accounting

| Section | Assertion macros | Status | Rust evidence |
|---|---|---|---|
| `storeLFQ(...)` | 0 | ported | The section's body is the comment "tested via MSstatsConverter tool" and asserts nothing, so there is no literal to reproduce. `lfq_reproduces_the_retained_cpp_output` reproduces all 2029 rows of the retained C++ output of the tool the section defers to, which is strictly stronger. |
| `storeISO(...)` | 0 | ported | As above; `iso_reproduces_the_retained_cpp_output` and `..._for_a_second_design` reproduce all 471 rows of each of the two retained outputs. |

### Known gaps

- The consensusXML read budget shipped in `src/format/consensusxml.rs` cannot
  load either fixture at its defaults: the isobaric map exhausts `max_work`
  while parsing its modified peptide sequences (each `AASequence` parse charges
  the modification registry's size) and the label-free map exhausts
  `max_payload_bytes`. `tests/msstats.rs` raises both, which are
  caller-supplied ceilings. Raising the shipped defaults is a change to a
  module this package does not own, and is flagged for the integrator.
- `TOPP_MSstatsConverter_subset` (`CMakeLists.txt:1655`) is not ported as a
  differential: its input is 5.4 MB and it adds only the `warnOnSubsetFiles_`
  path, which `tests/msstats.rs` already exercises through the produced
  warning. Adding it is a fixture-size decision.

### C++ findings

Recorded for `OpenMS_CPP_ISSUES.md`; the integrator owns that file.

- **OPENMS-MSSTATS-001 — a design gap silently annotates rows with the wrong
  sample.** `MSstatsFile.cpp:440` (and the isobaric counterpart at `:701`)
  indexes `path_label_to_sample`, `path_label_to_fraction`,
  `path_label_to_fractiongroup` and `run_map` with `std::map::operator[]`. A
  `(file, label)` pair the design does not declare is inserted with the value
  `0`, and `SampleSection::getFactorValue(0, factor)`
  (`ExperimentalDesign.cpp:971`) then returns row `0`'s condition and
  biological replicate. The output looks complete and is wrong. Proposed fix:
  use `at` inside a membership check and throw
  `Exception::MissingInformation`. Rust handling: `Error::MissingInformation`
  naming the file and the label.
- **OPENMS-MSSTATS-002 — an unrecognised summarization method writes zero
  intensities.** `MSstatsFile.cpp:178` initialises
  `MSstatsFile::Intensity intensity(0)` and assigns it only in the `max`,
  `min`, `mean` and `sum` branches. Any other non-`manual` name produces a
  complete CSV in which every intensity is `0`. The TOPP tool restricts the
  parameter, but the library API does not. Proposed fix: throw for an
  unrecognised name. Rust handling: `RetentionTimeSummarization::from_name`
  returns `Error::InvalidValue`.
- **OPENMS-MSSTATS-003 — `storeISO`'s run mapping is never populated.**
  `MSstatsFile.cpp:534` declares `msstats_run_to_openms_fractiongroup`, the
  isobaric loop never writes to it, and the informational loop at `:744` that
  prints "MSstats run N corresponds to OpenMS TechRepMixture M" therefore
  always prints nothing. Proposed fix: fill it, or delete it and its printout.
  Rust handling: the map is reported empty for the isobaric layout, and
  `tests/msstats.rs` asserts that.
- **OPENMS-MSSTATS-004 — the triple sanity check is dead.**
  `MSstatsFile.cpp:138` declares
  `set<tuple<std::string, std::string, std::string>> peptideseq_precursor_charge_run`
  described as a "sanity check that the triples (peptide_sequence,
  precursor_charge, run) only appears once", `:162` inserts into it, and
  nothing ever reads it or its size. The check the comment promises is not
  performed. Proposed fix: perform it or remove it. Rust handling: not ported.
- **OPENMS-MSSTATS-005 — `mean` and `sum` silently ignore repeated
  intensities.** `MSstatsFile.cpp:146` collects the intensities of one peptide
  ion in one run into a `std::set<Intensity>`, so two features with the same
  intensity at different retention times contribute once, and
  `meanIntensity_` divides by the number of *distinct* values rather than the
  number of samples. For exact float duplicates this changes the reported
  quantity. Proposed fix: collect into a vector. Rust handling: reproduced
  exactly, documented at the item, and unit-tested.
- **OPENMS-MSSTATS-006 — the label-free channel default cannot match a design
  label of anything but 1.** `MSstatsFile.cpp:106` and `:115` push the raw
  `channel_id` when the column header has one and `1u` when it does not, and
  `storeLFQ` uses that value directly as the design `Label`. A label-free map
  whose headers happen to carry `channel_id = 0` therefore looks up label `0`,
  which no design declares, and falls into OPENMS-MSSTATS-001. The source's own
  comment at `:113` notes the value "could be missing due to other reasons".
  Proposed fix: derive the label from the design's label count. Rust handling:
  reproduced, and the missing lookup is an explicit error rather than a silent
  zero.
