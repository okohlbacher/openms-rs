# Open items: every one decided

This is the answer to requirement **DEBT-03** of milestone M1: each of the 116
deduplicated open items, the 18 decisions owed and the stale records of the
open-item inventory, re-checked against the tree and given exactly one verdict.

- **Inventory checked:** `.planning/OPEN-ITEMS.md` (compiled 2026-09-21
  against `c3ac959`; planning notes, not tracked). Its numbering (`OI-nn`,
  `ST-nn`, `DEC-n`, `B1`-`B7`) is kept so the two can be read side by side.
- **Tree checked:** `main` at `f0990c0` (2026-09-27), 130 commits after the
  inventory. Where a wave-4 lane that has not merged yet is working on an item,
  the lane is named; its branch is evidence of intent, not of closure.
- **What "checked" means:** every cited site was opened again in the current
  tree. A line number in this document is the current one, not the inventory's.

## Verdicts

| Verdict | Meaning |
|---|---|
| **resolved** | Done in the tree; the commit is cited. Any document that still says otherwise is a stale record (section 2). |
| **scheduled** | Real work, owned by an M1 requirement: phase letter and requirement ID from `.planning/REQUIREMENTS.md` (A = DEBT/EVID, B = BLOCK/SHARED/CLI, C = CORE-01..05, D = CORE-06..12, E = CORE-13..18, F = STRUCT, G = UPSTREAM). |
| **deferred** | Deliberately not done in M1, with the reason. Nothing is lost: the reason says what would reopen it. |
| **decision owed** | Needs a ruling nobody below the lead may take. Each has a question, the options and what each costs, in section 5. |

An item whose question an existing ruling already answers is **scheduled**, not
**decision owed**, and the ruling is named. The rulings used that way are
Series-B **D1** (refuse exactly where the source is out of bounds, races, ends
the process or never terminates), **D16** (an unspecified but deterministic
`std::sort` permutation is in scope to reproduce), **D18** (D16 applied to
FileInfo `-c`), Series-A **D10** (source-compatibility load options; library
defaults stay strict), the user's **2026-09-15 platform decision** (the port
matches the Linux x86_64 Release build on every target), the user's **M1
definition of 2026-09-22** (every one of the 786 core headers closes) and the
project rule that **parallel lanes' duplicated logic is collapsed at
integration**.

Two D-series exist and collide: Series A is `docs/EARLY_TOPP_WORK_PACKAGES.md`
"Decisions taken" (D1-D10, process), Series B is the lead's list in
`docs/VALIDATION.md` "Lead decisions D1-D18". Every bare "D*n*" below is
Series B.

### Counts

| | resolved | scheduled | deferred | decision owed | total |
|---|---:|---:|---:|---:|---:|
| Open items OI-01..OI-116 | 10 | 60 | 23 | 23 | 116 |
| Lead decisions DEC-1..DEC-18 | 2 | 1 | 0 | 15 | 18 |
| New decisions found here (DEC-19..DEC-21) | | | | 3 | 3 |

Ten further items are **resolved in part**; each is counted under the verdict
of what remains, and the part already done is cited. The 23 items marked
**decision owed** collapse onto 13 decisions (for example OI-10 onto DEC-4,
six tool-level items onto DEC-21); each decision is stated once, in section 5,
and the item row points to it.

---

## 1. The 116 items

### 1.1 Shared math, sorting and NaN

| # | Item | Verdict | Evidence in the tree today | Next |
|---|---|---|---|---|
| OI-01 | `compute_rank` / `rank_correlation_coefficient` refuse a NaN on two grounds | **scheduled** A · DEBT-04, lane A3a | Still refused: `src/math/statistic_functions.rs:1378`, `:1441`. Ground (a), "the ranking sorts", is the argument D16 overturned; ground (b), the relative tie test, is ordinary defined `f64` arithmetic. | A3a ports both legs and measures the tie branch in its `../oracle/A3a-source-order` driver (mode `rank`). |
| OI-02 | The `_sorted` entry points return `UnsortedData` for a NaN | **decision owed** → DEC-19 | `statistic_functions.rs:511-516`, `:670-675`, `:728`. The C++ does not verify a caller's "already sorted"; it reads positions. | |
| OI-03 | `covariance`'s `mulsd` operand order is assumed, not measured | **scheduled** E · CORE-18 (MATH) | `statistic_functions.rs:1128-1147` still reads "not measured … No oracle row measures it". Lane A2b's EVID-02 re-audit moves `MATH/StatisticFunctions.h` from `native_equivalent` to `partial`; this is one of its named gaps. | One Release-build call on ibminode06 with two NaN payloads, read back as bits. |
| OI-04 | `median_of_sorted`'s `addsd` operand order, same hazard | **scheduled** E · CORE-18 (MATH) | `statistic_functions.rs:524-528`. | The same oracle row as OI-03. |
| OI-05 | Two Pearsons coexist; `pearson_correlation_coefficient` / `matthews_correlation_coefficient` substitute `f64::NAN` | **scheduled** E · CORE-18 (MATH) | `statistic_functions.rs:62-71` still "a known gap"; the bit-faithful one is `analysis::feature_finder_picked::scoring::source_pearson`. | Converge on `source_pearson` behind its own oracle row; the duplication rule decides the direction. |
| OI-06 | `consensus.rs`'s `it_aad` NaN was host-shaped | **resolved** `6aed3c3` | Routed through `math::x86_64::add`. `docs/VALIDATION.md` was corrected in `3942b51`; the work-packages bullet in `b6dee21` (this lane). | — |
| OI-07 | `range_utils.rs` sorts with `f64::total_cmp`, without D17's guard | **scheduled** A · DEBT-04, lane A3a | `src/kernel/range_utils.rs:721`. Same site as stale refusal B1. | A3a: `source_sort_by` + the ported `lower_bound`. |
| OI-08 | `feature_overlap_filter.rs` emulates `std::stable_sort` with a library `sort_by` | **scheduled** A · DEBT-04, lane A3a | `src/processing/feature_overlap_filter.rs:1105`. A `sort_by` with a comparator that is not a total order may also panic on Rust >= 1.81. | A3a: `source_stable_sort_permutation`. |

### 1.2 FileInfo and the mzML reader

| # | Item | Verdict | Evidence in the tree today | Next |
|---|---|---|---|---|
| OI-09 | Package A8: the rest of `FileInfo.h` | **scheduled** C · CORE-01 — resolved in part | Done: the mzXML, mzData, MGF and MS2 peak branches (`39fdeb6`, `eac9ac6`). Open: `-v` (`src/format/file_info/report.rs:44`, `Error::Unsupported`) although XSD validation now exists behind `xml-schema` (`3897138`); the pepXML, mzTab, trafoXML and PQP branches (`report.rs:434`); sqMass, XMass and MSP (`report.rs:437-443`). | `-v` can call `format::xml_schema`; the sqMass branch waits for lane F4's reader, MSP for its own port. |
| OI-10 | `-c` cannot report a repeated auxiliary array name or an MS-level-0 mass spectrum | **decision owed** → DEC-4 | `src/format/mzml.rs:1069` ("duplicate auxiliary array name"), `src/kernel.rs:811`. | |
| OI-11 | Three mzML reader gaps behind `#[ignore]` | **scheduled** C · CORE-01 (`FORMAT/HANDLERS/MzMLHandler.h`) | Ignores at `tests/file_info.rs:168`, `:325`, `:351`; the tripwire comment at `:1234`. Series-A D10 gives the shape: a source-compatibility read option, strict default. | |
| OI-12 | A3 request 5: selected-ion drift time onto the MS2 spectrum | **scheduled** C · CORE-01 (`MzMLHandler.h`) | Decision taken 2026-09-15 (follow the executed source); `tests/file_info.rs:494` still `#[ignore]`d; no lane opened. | |
| OI-13 | `file_info::Options` defaults strict on the dangling `defaultDataProcessingRef` | **scheduled** C · CORE-01 (`FORMAT/FileInfo.h`) | `src/format/file_info/model.rs:601` defaults `false`, which is the correct library default under Series-A D10; the tool passes `true`. The ignored case `tests/file_info.rs:210` runs the library at its default. | Run that case with the option set, as the tool does; do not flip the default. |
| OI-14 | Native differences 11 and 12 of the `-i` index check | **deferred** | `docs/FILE_INFO_CHECKS_SUPPORT.md:225`, `:248`. Reproducing either changes the decoder every indexed reader shares, and 12 would break random access; the C++ behaviour is a defect (`CPP-337`) and no OpenMS writer produces such a file. | Reopen if a reference writer is found that produces one. `CPP-337` goes to G · UPSTREAM-01. |
| OI-15 | `CPP-337`'s singleton consequence is a probe, not a pinned case | **deferred** | `docs/VALIDATION.md:1357`. Cost stated there: a fixture, two expectations, a test and five recounts. | Fold into the next re-execution of `../oracle/a6-fileinfo`. |
| OI-16 | `../oracle/a6-fileinfo/manifest.json` keeps the narrower `CPP-337` wording | **deferred** | `docs/VALIDATION.md:1173`. Correcting it re-executes 59 Release cases. | The next a6 re-execution. |
| OI-17 | `../oracle/a7-fileinfo` driver keeps an over-generalising sentence | **deferred** | `docs/VALIDATION.md:933`: "the next re-emission's to fix". | The next a7 re-emission. |
| OI-18 | A forced type the file contradicts exits 3/8 in C++, 6/11 here | **scheduled** B · BLOCK-01 (`FileHandler.h`) | `docs/TOPP_FILE_INFO_SUPPORT.md:136-144`. The `FileHandler` branch B0 is integrating does not close it (its ranked remainder is in `docs/FILE_HANDLER_SUPPORT.md` on that branch). | |
| OI-19 | `-log`, `-instance`, the timing line and `writeDebug_` | **resolved** `e64ade6`, `fb0ae8e`, `287c89f` | Tool messages reach the `-log` file (`e64ade6`); the timing line is printed (`fb0ae8e`); tool bodies write their debug lines (`287c89f`). `-instance` reproduces the source (`docs/TOPP_CLI_SUPPORT.md`). | — |

### 1.3 Identification XML, formats and readers

| # | Item | Verdict | Evidence in the tree today | Next |
|---|---|---|---|---|
| OI-20 | More than 14 modified peptide hits cannot be loaded (`parse_with_budget`) | **scheduled** C · CORE-01, lane F3a | `src/chemistry/sequence.rs:1147` unchanged on `main`; F3a replaces the per-sequence charge (`parse_metered`) and adds the upstream-fixture differential. | |
| OI-21 | consensusXML read budget refuses both MSstats fixtures at its defaults | **scheduled** C · CORE-01 (`ConsensusXMLFile.h`) | `docs/MSSTATS_SUPPORT.md:258-265`; `tests/msstats.rs:40-41` raise both ceilings. The `max_work` half has OI-20's cause (the modification-registry charge per parse), so F3a's fix should move it; the `max_payload_bytes` half needs a size-derived default. | Re-measure after F3a; then size the payload default from the input, as the mzML reader's are. |
| OI-22 | A `.pin` file written for a multi-accession PSM cannot be read back | **deferred** | `docs/PERCOLATOR_INFILE_SUPPORT.md:275`, `:335` (`OPENMS-PERCIN-007`). The source's own `store` then `load` is not a round trip either; the port is faithful to it, and relaxing the column check would admit misaligned rows. | An upstream fix, via G · UPSTREAM-01. |
| OI-23 | Three mzML-reader strictnesses block a combined imzML load | **scheduled** B · BLOCK-06 (`MzMLFile.h`) | `src/format/mzml.rs:2010` (ISO-8859-1 "requires transcoding"), `:2270` (`version` must start `1.1.`), and the forward `softwareRef`. | Series-A D10 shape: source-compatibility options, then `ImzMLFile` (C · CORE-01) passes them. |
| OI-24 | `ImzMLFile::preflight` does not see inline peaks | **scheduled** C · CORE-01 (`ImzMLFile.h`, `partial`) | `src/format/imzml_file.rs:1536`. | |
| OI-25 | `OpenMS/IMAGING/` has no Rust owner | **decision owed** → DEC-8 | The ledger has moved: `MSImagingExperiment.h` is now `complete`, `MSImagingGeometry.h`, `MSImagingRegion.h` and `IonImage.h` `partial`, each "owned by a later IMAGING package"; the types still live in `src/kernel/on_disc_imzml_experiment.rs:25-30`. The work itself is E · CORE-18 (IMAGING 3). | |
| OI-26 | `ImzMLFile::store`'s geometry overload is not ported | **scheduled** E · CORE-18 (IMAGING) + C · CORE-01 | `docs/IMZML_WRITER_SUPPORT.md:475-481`. Reachable once DEC-8 places `MSImagingGeometry`. | |
| OI-27 | `QcMLFile::collectQCData` (926 lines, ≥ 7 defects) | **scheduled** C · CORE-01 (`QcMLFile.h`, `partial`), metrics with D · CORE-10 | The ledger scope still says "explicitly unported". The per-defect question the inventory raised is answered by D1 and D16: reproduce every deterministic in-bounds behaviour, refuse exactly the UB. Bundling `qc-cv.obo` / `qc-cv-legacy.obo` is a data dependency, recorded like a crate. | |
| OI-28 | The `MzTab.h` export surface | **scheduled** C · CORE-01, lane F2 | `docs/MZTAB_SUPPORT.md`; F2 is on MzTab, MzTab-M, mzQC and OMS. | |
| OI-29 | `OMSFile.h` and `MzTabFile.h` | **scheduled** C · CORE-01, lane F2 | `OMSFile.h` `unmapped`, `MzTabFile.h` `partial`. `MzTabMFile.h` is `native_equivalent` and is re-audited under A · EVID-02. | |
| OI-30 | No mzXML schema embedded for `isValid` | **resolved** `3897138`, `4d59f17` | `src/format/xml_schema.rs:57`, `:105`, `:188` bundle `mzXML_idx_3.1.xsd`; `docs/MZXML_SUPPORT.md:73` maps `MzXMLFile::is_valid`. | — |
| OI-31 | `MzDataFile::is_semantically_valid` is `Error::Unsupported` | **scheduled** C · CORE-01 (`MzDataFile.h`, `partial`) | `src/format/mzdata.rs:672`; needs `psi-mzdata.obo` bundled and the `semantic-validation` machinery wired. | |
| OI-32 | mzData cannot carry integer/string data arrays | **scheduled** C · CORE-01 (`MzDataFile.h`) | `src/format/mzdata.rs:45` (`discard_unrepresentable`). The format has no element for them; the port's refusal-unless-discard matches the source's upstream `@improvement`. | Close as a documented native difference when the header closes. |
| OI-33 | `b_spline`, `polynomial` and `akima` transformation models | **scheduled** C · CORE-05 (`TransformationModelBSpline.h` `unmapped`, `TransformationModelInterpolated.h` `partial`) | `src/format/transformation_xml.rs:175`, `:352` return `Error::Unsupported`. | |
| OI-34 | `SpectrumLookup`, `SpectrumMetaDataLookup`, `SpectrumNativeIDParser`; a caller's `scan_regex` | **scheduled** D · CORE-08 (METADATA) | All three are now `partial` with consumer-local subsets. The stated reason for refusing `scan_regex` was stale: `crate::concept::boost_regex` landed in `854d996`, the day after `docs/MASCOT_XML_SUPPORT.md` was written (sentence corrected, `4754da8`). | |
| OI-35 | `IndexedMzMLFileLoader` and the `Interfaces::` pointer overloads | **scheduled** C · CORE-01 (`IndexedMzMLFileLoader.h`, `partial`) + E · CORE-18 (`INTERFACES/*`) | Ledger scope: every operation exists elsewhere. The pointer overloads wait for the `Interfaces` layer. | |
| OI-36 | `ReactionMonitoringTransition`, `MRMTransitionGroupPicker` | **scheduled** A · DEBT-02, then C · CORE-02 | `ReactionMonitoringTransition.h` is `evidence_requires_review` (DEBT-02 decides it); `MRMTransitionGroupPicker.h` `unmapped` in OpenSWATH. `src/kernel/mrm.rs:28`, `:103`. | |
| OI-37 | `-write_ctd`, `-write_cwl`, `-write_json` refused on every tool | **resolved** `1a773e4`, `c36d087` | `-write_ctd` writes the CTD (`src/cli.rs:1163-1167`); the four CWL/JSON writers end with exit 8 and the Release build's own "TDL support is not available" message, as the Release build does without TDL (`src/cli.rs:1169-1174`). | — |
| OI-38 | `get_type` fails on a directory where C++ reports an unknown type | **scheduled** B · BLOCK-01 (`FileHandler.h`) | `src/cli.rs:2008-2011`: the tool layer already treats the failure as an undetermined format; the library still fails. | |

### 1.4 Kernel, processing and feature finding

| # | Item | Verdict | Evidence in the tree today | Next |
|---|---|---|---|---|
| OI-39 | `MSSpectrum` has no stored ion-mobility format or peak type | **scheduled** E · CORE-18 (`KERNEL/MSSpectrum.h`, `partial`) | `src/kernel.rs:176`, `:291` carry `drift_time_unit` only. | |
| OI-40 | Three deferred `MSExperiment` members | **scheduled** B · BLOCK-03 | Ported on the `MSExperiment` branch B0 is integrating (`set_2d_data_mass_traces`, `first_product_spectrum_index`, `Display`, with Release oracle cases). Not on `main` yet. | B0 merges; the verifier re-checks. |
| OI-41 | The rest of `IMDataConverter.h` | **scheduled** E · CORE-18 (IONMOBILITY) | `src/kernel/im_data_converter.rs` has none of `reshapeIMFrameTo*`, `splitExperimentByIonMobility`, `set/getIMUnit`. Waits for `IonMobilityBinning`, `FeatureFinderMetabo`, `FeatureFinderMultiplex`, `SpectraMerger`, `IMDataArrayUtils`. | |
| OI-42 | The cross-voltage FAIMS merge has no C++ oracle | **decision owed** → DEC-5 | Unchanged. | |
| OI-43 | `Math::RandomShuffler` duplicates the `pub(crate)` `DecoyRandom` | **scheduled** B · BLOCK-06 (`MathFunctions.h`) | `src/chemistry/decoy_random.rs:33`. The duplication rule answers the inventory's question: promote `DecoyRandom` into `math` (a sink, so no cycle) and have both consumers use it. | |
| OI-44 | `TraceFitterParams`' `Param` round trip is `i32`-bounded | **scheduled** B · SHARED-01, lane B6 | `tests/trace_fitter.rs:381-389` pins the disagreement. `crate::param`'s restriction check converts to `i32`; the `DefaultParamHandler` trait lane owns the `Param` contract. | |
| OI-45 | LM solver fidelity beyond the fixtures, macOS arm64 | **deferred** | The user's 2026-09-15 platform decision makes Linux x86_64 Release the reference on every target; there the solver is bit-exact on all 141 traced fits (`docs/TRACE_FITTER_SUPPORT.md`, *Known gap*). `tests/gauss_trace_fitter.rs:1593` is report-only by design. | Nothing, unless the platform decision is revisited. |
| OI-46 | Is bitwise equality with Eigen required? | **decision owed** → DEC-6 | `docs/EGH_TRACE_FITTER_SUPPORT.md:369` "not yet recorded". | |
| OI-47 | `atan` fallback off the reference node | **deferred** | `src/analysis/feature_finder_picked/glibc_libm.rs`: glibc's `__atan_fma` is LGPL with no licence-clean upstream; the measured departure is recorded. | A licence-clean source for `__atan_fma`. |
| OI-48 | `2^32 - 5` score arrays never run uncapped | **deferred** | `docs/VALIDATION.md`: ~928 GiB, 93 % of the reference node. | A node with the memory. |
| OI-49 | `allow_duplicate_positions` / `allow_unsorted_positions` refused | **scheduled** E · CORE-15 (`PeakPickerIterative.h`, `partial`) | `src/processing/iterative.rs:105-107`, `:193-194`. D16 and D1 answer the inventory's question: a `std::map<double,double>` overwrite is deterministic and in bounds, so it is ported, not refused. | Port it after A3a's `iterative.rs` change (B3) lands. |
| OI-50 | Chromatogram and iterative pickers call the strict S/N `estimate` | **scheduled** E · CORE-15 + C · CORE-02 | `docs/SIGNAL_TO_NOISE_SUPPORT.md:589-600`; `PeakPickerChromatogram.h` (OpenSWATH) and `PeakPickerIterative.h` are both `partial`. The doc names the exact call. | |
| OI-51 | Signal-to-noise signed-overflow sites stay refused | **deferred** | D1: each site is UB and the refusal is as wide as the UB; emulating needs > 2^31 points and 64-90 GB per run. | Reopen only if a real input approaches 2^31 points. |
| OI-52 | B6 (c): the overall seed score against platform `powf` | **resolved** — see DEC-7 | `60d050b`. | — |
| OI-53 | Weighted intensity spread for high-resolution binned sums | **deferred** | `src/kernel/gap_closures.rs:97-99`: an upstream `@todo` with no upstream implementation, so there is nothing to port. `BinnedSumAgreeingIntensities.h` is `complete`. | An upstream implementation. |
| OI-54 | `FFT` is not bit-identical to evergreen | **deferred** | `src/math/fft.rs:57-62`: pinned against a naive DFT (an independent invariant). `KernelDensityEstimation.h` is `complete`; its upstream expectations pass. | An oracle row through `KernelDensityEstimation` on the Release build, if a KDE consumer ever needs bit-identity. |
| OI-55 | MD5 is ported inline, not taken from a crate | **scheduled** C · CORE-01 (`ImzMLFile.h`, `partial`) | `src/format/imzml_writer.rs:812-817`. Crates-before-code applies; `md-5` is small and pure Rust, so it may be a plain dependency with a `docs/THIRD_PARTY_CRATE_DECISIONS.md` entry. | |

### 1.5 Ledger, tooling, CI and infrastructure

| # | Item | Verdict | Evidence in the tree today | Next |
|---|---|---|---|---|
| OI-56 | The decision register stops at D17 | **scheduled** A · DEBT-03 (the lead's file) — resolved in part | `9466368` retitled both lists "Lead decisions D1-D18" (`docs/VALIDATION.md:485`, `:1865`). Three links still point at the old anchor `#lead-decisions-d1-d17`: one fixed here (`docs/EARLY_TOPP_WORK_PACKAGES.md:1533`, `4754da8`); `docs/VALIDATION.md:204` and `:1869` are the lead's. | Change both to `#lead-decisions-d1-d18`. |
| OI-57 | Provenance manifests not registered for hash verification | **scheduled** A · EVID-03 — resolved in part | The five "never added" (`consensus_display`, `im_types`, `on_disc_experiment` in `a80a370`; `imzml_handler`, `on_disc_imzml` in `051f475`) were registered before the inventory was compiled; their support documents were the stale part (fixed, `4754da8`). Still outside `current_sdk_reference_manifests`: the six with a foreign `target_verification` key (`feature_finder_picked`, `lm_eigen_path_differential`, `mzml_reader_scale`, `lm_budget_differential`, `mzml_header_leniency`, `isotopes_source_precision`). | Rename their key, then register (`SOURCE_PROVENANCE.json` is the lead's). |
| OI-58 | Ledger, coverage and CI wiring for the three imzML packages | **resolved** `1b6c30a`, `051f475` | CI lines in `1b6c30a` (now generated by `tools/ci_matrix.py`, `5c4a44a`); `docs/doc-coverage.json` and the ledger entries in `051f475`. The three support documents were corrected (`4754da8`). | — |
| OI-59 | Ledger says `setSkipXMLChecks` is not ported | **resolved** `a80a370` | `docs/core-sdk-reviewed-apis.json` now reads "setSkipXMLChecks IS ported". | — |
| OI-60 | Citation-checker item 7, `Pins.packaged()`, the paraphrased A6 defect | **scheduled** A · EVID-04, lane A4 — resolved in part | `cf899a0` measured item 7 (66 bare ranges down to 3). Remaining: 3 bare ranges and ~35 unreachable names (A4's count), `Pins.packaged()` without its own test, the paraphrase class. | |
| OI-61 | `is_x86_feature_detected!` / `cfg!(target_feature)` fold at compile time | **deferred** | A standing hazard, documented in `docs/VALIDATION.md` and the porting skill; no lint exists. | A repository check that refuses both outside an allow-list would make it enforceable; nobody has asked for one. |
| OI-62 | `Record::finish` has two code-generation states | **deferred** | A measurement rule for performance comparisons, not work. | — |
| OI-63 | `check_schema_feature_graph.py` runs only in CI | **decision owed** → DEC-9 | `openms-ci-sweep.sh` sweeps only the `cargo` lines of the `portable-feature-graph` job, so the checker itself still never runs outside CI. | |
| OI-64 | Gate and node hazards | **deferred** — resolved in part outside the repository | Done: a second gate host (`OPENMS_GATE_HOST` takes kim, dax, spock, ibminode05, ibminode06) and a sweep of every CI line (`openms-ci-sweep.sh`). The rest (`--delete-excluded` per slot, the `/usr/local/bin/cc` shadow, linker via `CARGO_TARGET_*_LINKER`) are standing operating rules. | — |
| OI-65 | Test binaries are not guarded against a non-FMA processor | **deferred** | Accepted by the lead; residual risk in `docs/FMA_BUILD_FLAG.md`. | — |
| OI-66 | SIMD base64 decoder | **decision owed** → DEC-3 | `docs/THIRD_PARTY_CRATE_DECISIONS.md:46` still "Pending; measured, not adopted". | |
| OI-67 | `mimalloc` measurement withdrawn | **deferred** | Only matters if DEC-14 is reopened; then it needs an override build or a real `#[global_allocator]`. | — |
| OI-68 | Stale `TOPP_DTAExtractor_{1,2,3}_output.dta` could not be checked | **resolved** by this lane (`OPEN_ITEMS.md`) | The test data is now readable at `0cb15f2`: `topp/DTAExtractor_1_output.dta` and `_2` hold `120 100`, `_3` holds `110 100` / `120 200`, where the pinned tool writes `120.0 100`. The reference files are textually stale and pass upstream only through FuzzyDiff's numeric comparison. | A `CPP-` entry (test-data hygiene, cosmetic) for the lead's log; `docs/VALIDATION.md:2404` ("could not be checked") is the lead's to update. |
| OI-69 | The cycle gate counts a `crate::` path in a comment | **decision owed** → DEC-20 | `tools/check_module_cycles.py` still reads raw text; `055b4e3` pinned what it counts but not comments. | |
| OI-70 | The `format` ↔ `system` cycle forces a duplicated CPUID read | **scheduled** F · STRUCT-01 | `docs/module-cycles.json` still lists `format`/`system` as mutually dependent; the duplicate is at `glibc_libm.rs:825-852`. | Collapse onto `system::cpu_features` once the pair is cut. |
| OI-71 | A resource-refusal exit code | **decision owed** → DEC-17 | Unchanged. | |
| OI-72 | `render_list_parameters` lives in the wrong module | **scheduled** C · CORE-01 (`MzMLHandler.h`, the writer) | `src/cli/tools/peak_picker_hi_res.rs:537`. The source's mzML writer renders a list `DataValue` as text; following the source puts the fix in the writer, which answers "where". | |
| OI-73 | Five private float formatters; `param/value.rs` breaks exact ties upward | **scheduled** B · BLOCK-06 (`StringUtils.h`) | `src/format/file_info/text_format.rs:216-223` ("Follow-up"): `ostream_g` (`mascot_generic.rs:1022`), `pepxml.rs`, `mzxml.rs`, `posterior_error_probability.rs`, `param/value.rs:521`, `:560`. The duplication rule decides "consolidate". The tie defect is a live wrong digit. | |
| OI-74 | `TraceFitter` signature changes need the integrator | **deferred** | A coordination rule, not work; no change is pending and `TraceFitter.h`, `GaussTraceFitter.h`, `EGHTraceFitter.h` are `complete`. | — |
| OI-75 | Four carried-open framework questions on `TOPPBase.h` | **scheduled** B · BLOCK-05 | Never enumerated (`docs/VALIDATION.md:2311`). Superseded by BLOCK-05's ranked remainder (`-write_ini` float spelling, JSON syntax-error text, parse-error wording, `MissingInformation` exit 7 vs 8, Windows width). | |
| OI-76 | `SummaryLimits::default().max_work` refuses a file the reader reads | **scheduled** B · BLOCK-03 | Size-derived defaults on the `MSExperiment` branch B0 is integrating ("Size-derived default ceilings for the linear experiment summaries"). | |
| OI-77 | The chunk merge's stability is untested at a chunk boundary | **scheduled** E · CORE-18 (`KERNEL/MSSpectrum.h`) | `docs/VALIDATION.md:3474-3476`. | One fixture with equal m/z spanning a chunk boundary. |
| OI-78 | 2 of 28 adversarial picker cases not bit-identical (zero intensity total) | **scheduled** E · CORE-15 (`PeakPickerHiRes.h`, `partial`) | `docs/VALIDATION.md:2781` "open note". | |
| OI-79 | WindowMower unresolved at n = 3; `CPP-338` one observed crash | **deferred** | A benchmark-precision question, not a fidelity one; `CPP-338` goes to G · UPSTREAM-01. | More repetitions when a benchmark wave runs. |
| OI-80 | HPC-only `#[ignore]`d tests | **deferred** | Environmental by design: they read `/ceph/ibmi/abi/oliver/bench/openms4/…`. 27 `#[ignore]` attributes today (`mzml_reader_scale` 8, `file_info` 5, `mzml_writer_scale` 3, `topp_threads` 3, `topp_baseline_filter_edges` 3, one each in five more files). | — |
| OI-81 | `ProgressLogger` not threaded; `LogStream` not wired to the kernel | **scheduled** C · CORE-01, E · CORE-18, B · BLOCK-02 — resolved in part | Done: Mascot generic (`ba02d5a`), mzXML (`25b6342`), imzML writer (`3fc42af`). Open: `src/format/imzml_handler.rs` and `src/kernel/on_disc_imzml_experiment.rs`; `LogStream` to the kernel is BLOCK-02. | |
| OI-82 | Three private re-implementations | **scheduled** D · CORE-11 + B · BLOCK-03 — resolved in part | `DataProcessingUtils.h` is `complete`, with `has_isobaric_analyzer` private by ledger decision until a second consumer (`src/kernel/map_operations.rs:1494`). Open: the `UniqueIdIndexer.h` subset (`partial`, CONCEPT) and the `MSExperiment` path helper (`:1505-1509`). | |
| OI-83 | `KERNEL/DimMapper.h` deferred to TOPPView | **scheduled** E · CORE-18 | The user's M1 definition counts every core header, and CORE-18's "KERNEL 4" includes `DimMapper.h` (`unmapped`). The exclusion predates it; `docs/IM_TYPES_SUPPORT.md` now says so (`4754da8`). | |

### 1.6 Lane items carried forward

| # | Item | Verdict | Evidence in the tree today | Next |
|---|---|---|---|---|
| OI-84 | `ConvexHull2D::add_points` merges with `f64::min`/`max` | **scheduled** D · CORE-06 (`ConvexHull2D.h`, `partial`) | `src/kernel/geometry.rs:264-265` and the batch path from `:284`. | The source's keep-first rule for `-0.0` vs `+0.0`. |
| OI-85 | DTAExtractor's reversed `-rt` range | **decision owed** → DEC-21 — resolved in part | BaselineFilter's centroided-input warning, the pair's other half, closed in `aaf32e8`. The range: `docs/TOPP_CLI_SUPPORT.md:883-891`, "no test covers it yet". | |
| OI-86 | The generator cannot express "reviewed, not ported" | **resolved** `94715ce` | Reviewed ledger entries now carry each header's real port: `IsotopeCluster.h` `partial`, `ProcessingStep.h` `complete`, `KERNEL/MassTrace.h` `native_equivalent` (re-audited under A · EVID-02, lane A2b). None is left at `evidence_requires_review` by name match. | — |
| OI-87 | Clippy lints under 1.85, invisible to CI | **scheduled** A · DEBT-03 follow-up (the lead's CI file) | `.github/workflows/rust.yml` still runs clippy on `stable` only (`:146`, `:182`). Measured here, `cargo +1.85.0 clippy --locked --all-features --all-targets` on dax, 2026-09-27: **12** distinct warnings, not 24 — 7 `format_collect`, 3 `comparison_chain`, 2 `excessive_precision`; 3 in the library (`src/analysis/feature_finder_picked/extension.rs:584`, `fitting.rs:315`, `src/concept/math_functions.rs:510`), 2 in `b_spline.rs`'s unit tests (`src/processing/spline/b_spline.rs:1031`, `:1087`), 7 in `tests/`. None of the inventory's `multiple_testing.rs` or `comparison.rs` sites remains. | Fix the 12, then either a 1.85 clippy line or a stated reason not to have one. |
| OI-88 | A3 reader minors | **scheduled** B · BLOCK-01 + B · BLOCK-06 | Mobility `" -35"` still a parse error; `get_type` still reads at most 64 KiB (`docs/MZML_MOBILITY_SUPPORT.md:151`); the nine `getType` literals of `FileHandler_test.cpp` are unported. | The closing lanes re-check each sub-point; the list is in `docs/EARLY_TOPP_WORK_PACKAGES.md` (A3 minors). |
| OI-89 | C3 / FuzzyDiff minors | **decision owed** → DEC-21 | Still divergent on the reference platform: the Linux Release build rejects `1e-400`, the port accepts it (`docs/FUZZY_STRING_COMPARATOR_SUPPORT.md:168-177`). FuzzyDiff is the ninth validated workflow; no core header owns it. | |
| OI-90 | CLI-1 requests 7 and 8 | **scheduled** B · CLI-01 + B · BLOCK-06 | `update_diagnostics` still at `src/cli.rs:1629`; `to_int32` twice (`src/cli/context.rs:849`, `src/system/update_check.rs:219`). | |
| OI-91 | `-in /dev/null` exits 2 with a false message | **decision owed** → DEC-13 | `docs/THIRD_PARTY_CRATE_DECISIONS.md:47` still "Pending; candidate only". | |
| OI-92 | Doc-accuracy debts | **scheduled** A · DEBT-01 (lanes A5a/A5b) — resolved in part | Fixed here: the B3b reference in `docs/EGH_TRACE_FITTER_SUPPORT.md` (`4754da8`). Still stale, in code and so outside this lane: `src/math/multiple_testing.rs:1121` ("2.5 million", measured 5 million). The rest are listed in `docs/EARLY_TOPP_WORK_PACKAGES.md` per lane. | Each partial header's re-read corrects its own; the rest go with their owners. |
| OI-93 | P1 picker minors, incl. the 256 MiB acquisition-copy ledger | **decision owed** → DEC-17 | `src/processing.rs:73` (`256 * 1024 * 1024`). The other sub-points are E · CORE-15 (`PeakPickerHiRes.h`). | |
| OI-94 | The mzML reader refuses duplicate record ids (`CPP-261`) | **scheduled** B · BLOCK-06 (`MzMLFile.h`) | `src/format/mzml.rs:2431`, `:3903`, `:3925`. Series-A D10 answers the "D10 candidate" question: a source-compatibility read option, strict default. | |
| OI-95 | P2 header minors | **scheduled** C · CORE-01 (`MzMLHandler.h`, `IndexedMzMLHandler.h`) | Unchanged; list in `docs/EARLY_TOPP_WORK_PACKAGES.md` (P2 minors). | |
| OI-96 | CLI-2 minors | **decision owed** → DEC-21 | Tool-level divergences in validated tools (directory `-in` 8 vs 3, NaN `struc_elem_length` 6 vs 0). The directory case is also EVID-01's known divergence. | |
| OI-97 | A5 minors | **decision owed** → DEC-21 | `detect_type` still duplicated (`src/cli/tools/file_info.rs:257`, `src/format/file_info/report.rs:390`); the duplicate itself is collapsed by the duplication rule wherever DEC-21 lands. | |
| OI-98 | C5 minors | **decision owed** → DEC-21 | Tool-level (`FeatureFinderCentroided`). | |
| OI-99 | B7 minors, incl. a panic in `extend_mass_traces` | **decision owed** → DEC-21 | **Live panic on public input**: `src/analysis/feature_finder_picked/extension.rs:364-367` indexes `pattern.peak[p]` then `pattern.spectrum[p]`; a public `IsotopePattern` with `spectrum` shorter than `peak` panics. `FeatureFinderAlgorithmPicked.h` is `complete`, so no requirement revisits it. | Fix first whatever DEC-21 decides: "no panics on untrusted input" is a standing constraint. |
| OI-100 | B6/B8/B9/B11 minors | **decision owed** → DEC-21 | Unchanged; all in `complete` FFAP/FAIMS headers. | |
| OI-101 | B10-FFC-ACCEPT's tight comparison and the LM budget sweep | **decision owed** → DEC-16 | Unchanged. | |
| OI-102 | Neither implementation finishes the 43,745-spectrum run | **decision owed** → DEC-16 | Unchanged. | |
| OI-103 | `map_collect` builds a pool unconditionally; `Threads::from_cli` on a negative count | **scheduled** D · CORE-11 (CONCEPT) | `src/concept/parallel.rs:99`, `:128`. Coupled to DEC-12. | |
| OI-104 | Two full scans on the picker's path | **deferred** | Performance, not fidelity; the inventory's own record says a fusion must be measured first (one cost 1.0 s of wall). | A benchmark wave. |
| OI-105 | `MSExperiment::max_intensity()` missing | **scheduled** B · BLOCK-03 | Only `src/kernel/ranges.rs:761` and `mass_trace.rs:426` define one; the `MSExperiment` branch does not add it. | Port `getMaxIntensity`'s combined semantics or cross-reference `ranges`. |
| OI-106 | The window mower's input-point ceiling refuses real spectra | **resolved** `2f173d9` | `src/processing/window_mower.rs:62-75`: size-derived (`max_work + work_per_point * n`). | — |
| OI-107 | Where the acquisition-copy ledger meters | **decision owed** → DEC-17 | `src/processing.rs:52` meters per experiment, `baseline.rs:219` per spectrum; no `AcquisitionCopies::for_records`. | |
| OI-108 | Wave-3a tools never opted into `-threads` | **decision owed** → DEC-21 — resolved in part | `PeakPickerHiRes` in `3b943e4`; FuzzyDiff too. Still missing: `FileInfo`, `FeatureFinderCentroided`; `tests/topp_threads.rs:248` lists six tools. | |
| OI-109 | Five formatters print `nan` for a sign-bit NaN | **scheduled** B · BLOCK-06 (`StringUtils.h`) | `src/format/mascot_generic.rs:1023-1024` and the four others. The 2026-09-20 spelling rule (glibc's `-nan`) already decides the answer; consolidating onto `text_format` (OI-73) fixes all five. | |
| OI-110 | Benchmark-harness debts | **deferred** | Performance tooling; M1 is fidelity. `Cargo.toml` has no `debug = 1` release profile. | The next benchmark wave. |
| OI-111 | Two manifests cannot be registered as written | **scheduled** A · EVID-03 | `tests/data/topp_feature_finder_centroided_provenance.json:388` cites an absolute `/ceph` path; `baseline_filter_edges` and `mzml_mobility` still lack `external_reference_note` and `origin_key`. | |
| OI-112 | `.reference/` exists only in the main checkout | **deferred** | Documented; every worktree lane meets it. | — |
| OI-113 | C1 oracle follow-ups and the macOS policy note | **scheduled** A · EVID-01 | The early-bundle oracle is one of the three EVID-01 re-runs on the Linux Release build, which retires the macOS policy note and is the natural place to bind the hashes and the crash exit-code class. | |
| OI-114 | Cluster infrastructure never closed | **deferred** | Outside the repository (`perf_event_paranoid = 4`, `/usr/local/bin/cc`, the idle runner on ibminode05). The admin wishlist is its home. | — |
| OI-115 | The early-bundle deferral list; the S1 review "pending authorization" | **scheduled** B · BLOCK-01, B · SHARED-02, C · CORE-01 (lane F4 for SQLite) | M1's roadmap supersedes the early bundle's order-only list. The review awaiting authorization is DEC-18 (c). The S0 "auxiliary-array recovery contract" candidate goes with F4's sqMass work. | |
| OI-116 | Package C6 never reported merged | **decision owed** → DEC-16 | What C6 lacked has landed; what remains is B10's acceptance statement, which DEC-16 governs. | |

---

## 2. Stale records

### 2.1 The inventory's twelve (ST-01..ST-12)

| # | Now |
|---|---|
| ST-01..ST-04 | Fixed in `3942b51`. |
| ST-05 | `docs/VALIDATION.md:2461` is a gate record of what held at `0510382`; it is correct as history. Annotating it "NaN half reversed by D18" is the lead's call. |
| ST-06 | Still stale, in code: `src/analysis/feature_finder_picked/helper_structs.rs:10-11` says the three trace fitters are "not ported yet". Left for a code lane. |
| ST-07 | Still stale, in code: `src/kernel/map_operations.rs:1313` says `ConsensusFeature` "has not yet gained" a `Display`; `src/kernel/consensus_display.rs` has it. Left for a code lane. |
| ST-08 | Historical wave record; superseded as the inventory says. |
| ST-09 | Already fixed in both documents. |
| ST-10 | Historical; the log now holds 354 entries (lane A4). DEBT-05 regenerates the index. |
| ST-11, ST-12 | Historical wave records. B6 (c), which ST-12 thought survived, is also closed (DEC-7), and both of its "still open" sentences now say so (`4754da8`). |

### 2.2 Stale records this pass found and fixed (one sentence each)

| Document | What it said | Commit |
|---|---|---|
| `docs/CONSENSUS_DISPLAY_SUPPORT.md`, `docs/IM_TYPES_SUPPORT.md`, `docs/ON_DISC_EXPERIMENT_SUPPORT.md`, `docs/ON_DISC_IMZML_SUPPORT.md` | Manifest not registered for hash verification | `4754da8` |
| `docs/ON_DISC_IMZML_SUPPORT.md`, `docs/IMZML_WRITER_SUPPORT.md`, `docs/IMZML_FILE_SUPPORT.md` | CI wiring and doc coverage "the integrator's step" | `4754da8` |
| `docs/ON_DISC_EXPERIMENT_SUPPORT.md` | The `setSkipXMLChecks` scope string is wrong | `4754da8` |
| `docs/EGH_TRACE_FITTER_SUPPORT.md` | Solver root cause "under investigation in lane B3b" | `4754da8` |
| `docs/EARLY_TOPP_WORK_PACKAGES.md` (twice) | B6 (c) still open | `4754da8` |
| `docs/EARLY_TOPP_WORK_PACKAGES.md` | Link to `#lead-decisions-d1-d17` | `4754da8` |
| `docs/EARLY_TOPP_WORK_PACKAGES.md` | `it_aad` still host-shaped | `b6dee21` |
| `docs/MASCOT_XML_SUPPORT.md` | "This crate has no regular-expression engine" | `4754da8` |
| `docs/IM_TYPES_SUPPORT.md` | `DimMapper.h` is the KERNEL domain's deliberate exclusion | `4754da8` |

### 2.3 Stale records in files this lane may not edit

| File | Record | Correct statement |
|---|---|---|
| `docs/VALIDATION.md:204`, `:1869` | `#lead-decisions-d1-d17` | `#lead-decisions-d1-d18` (`9466368` renamed the heading) |
| `docs/VALIDATION.md:2404-2409` | DTAExtractor references "could not be checked" | Checked (OI-68): stale `120` against the tool's `120.0` |
| `OpenMS_CPP_ISSUES.md:4822` (`CPP-271`) | "whether the port follows the source is open for the lead" | Decided 2026-09-15, followed since `1c27421` (DEC-2) |
| `src/analysis/feature_finder_picked/helper_structs.rs:10-11`, `src/kernel/map_operations.rs:1313` | ST-06, ST-07 | See 2.1 |
| `src/math/multiple_testing.rs:1121` | "2.5 million epsilons" | 5 million (OI-92) |

---

## 3. Refusals

### 3.1 The seven stale refusals (inventory section 3B) — A · DEBT-04, lane A3a

All seven are **scheduled** under DEBT-04 and are lane A3a's brief: B1
`range_utils.rs` (with OI-07), B2 `spectrum_mobility.rs` (`is_sorted_by_im`
and `sort_by_ion_mobility`), B3 the `std::stable_sort` leg in `iterative.rs`,
B4 `compute_rank` (OI-01), B5 `faims_helper.rs` (A3a models the source's
`std::set` insertion instead of porting `_Rb_tree`), B6 `debug.rs`'s abort map
and B7 `chromatogram_tools.rs`. A3a also takes OI-08.

### 3.2 Adjacent

- **`src/kernel/spectrum_mobility.rs`'s length-mismatch refusal** (the
  "borderline" case): the source's answer is defined, so this is a strictness
  question, not a D16 one. It is folded into **DEC-4**.
- **`src/processing/spline/smoothing.rs`** needs a one-line note that its probe
  ran on libc++, whose sorting networks `source_sort.rs` does not reproduce.
  **scheduled** E · CORE-15.
- **~17 order-substituting divergences** that refuse nothing
  (`src/comparison.rs`, `src/format/qcml.rs`, `src/chemistry/ims_alphabet.rs`,
  `src/analysis/feature_finding_metabo.rs`, and the documents listed from
  `docs/CONVERSION_HELPER_SUPPORT.md`): D16 puts each in scope. Those in
  `partial` headers go with their header's requirement; those in `complete`
  headers are **DEC-21**.

## 4. The C++ defect log

- **Index covering 58 of 354 entries; header naming one revision:** A · DEBT-05,
  lane A4. A4 found that `82ce5b3` *is* used, by 76 early links, and is an
  ancestor of `bc9cc12`, so the header is wrong only in claiming it for every
  entry.
- **`CPP-059`, `CPP-060` without a status line; two undeclared evidence
  tiers:** A · DEBT-05.
- **Triage, the 48 executed memory-safety findings, `CPP-341` first:**
  G · UPSTREAM-01..04. The port files nothing.

---

## 5. Decisions

### 5.1 Resolved or scheduled since the inventory

- **DEC-2 — `mass_trace:min_spectra = 1`: resolved.** The lead decided on
  2026-09-15 to follow the source; `1c27421` ported it and
  `tests/feature_finder_picked_seeds.rs:820`
  (`min_spectra_one_follows_the_source_and_finds_no_seed`) pins 0 seeds, 0
  features, exit 0 (`docs/FEATURE_FINDER_PICKED_SUPPORT.md:1090`). The inventory
  read `OpenMS_CPP_ISSUES.md:4822`, which is stale.
- **DEC-7 — B6 (c), the overall seed score: resolved.** `60d050b` ported glibc's
  `__powf_fma` from Arm optimized-routines as `glibc_powf.rs`; it equals the
  executed library on every probed pair, all `2^32` bases with exponent `1/3`
  included (`CPP-272`'s *Rust handling*, wave-5 decision D5). Nothing was left
  to rule on.
- **DEC-10 — reporting the C++ defects upstream: scheduled** G · UPSTREAM-01..04.
  The user's M1 definition settles it: the port prepares the reports, `CPP-341`
  first, and the user decides what is filed.

### 5.2 Owed, carried from the inventory

Each is one question. Options are the ones on the table; the cost is the one
the records support.

**DEC-1 — mzML whole-document writer: reproduce the C++'s dangling
`dataProcessingRef`s (`CPP-172`), or keep valid output?**
A: reproduce by default — every `mzml::write` caller then emits invalid mzML.
B: keep content deduplication — a permanent, pinned divergence from
`MzMLFile::store` (`tests/mzml_source_file_round_trip.rs:547`).
C: a third, opt-in write policy beside `SourceDangling` — one more policy to
document and test. *C matches Series-A D10's shape (source behaviour opt-in,
strict default).*

**DEC-3 — adopt `base64-simd` for the mzML decode path?**
A: adopt, behind a feature — 2.03x on the hot path (≈ −0.9 s on 2.3 GB);
internal `unsafe` in a dependency, and the crate must be in the nodes'
`--offline` cache. B: decline — the headroom stays. Now also closes
`Base64.h`, which `1712a03` made `partial` for exactly this pair.

**DEC-4 — strictness where the source's answer is defined: may the reader and
kernel refuse input the C++ loads?** Covers FileInfo `-c`'s two classes (a
repeated auxiliary array name, `src/format/mzml.rs:1069`; an MS-level-0 mass
spectrum, `src/kernel.rs:811`) and `spectrum_mobility.rs`'s length-mismatch
refusal. A: a source-compatibility read option that `-c` sets (Series-A D10) —
one more `ReadOptions` switch, strict default kept. B: relax the kernel
invariant for MS level 0 — every other consumer loses it. C: leave it — `-c`
permanently cannot report two of the corruption classes it exists for. *D18's
reasoning points at A.*

**DEC-5 — is a patched C++ build acceptable evidence** for the cross-voltage
FAIMS merge (`CPP-278` stops the pinned tool first)? A: yes — needs provenance
rules for a patched build and an amendment to "no C++ in the repository" saying
how one is pinned. B: no — the merge stays the lane's only below-tier-1
evidence. C: wait for an upstream fix — indefinite; nothing is filed yet.

**DEC-6 — trace fitters: is the residual off-glibc gap accepted?** The user's
2026-09-15 decision already fixes the reference (Linux x86_64 Release, every
target), and on it all 141 traced fits are bit-identical. What is still
unrecorded is the acceptance criterion elsewhere: 1e-14 relative with the
mandated `libm` crate (10,380 of 11,023 bitwise). A: accept 1e-14 off glibc and
record it — the criterion stops being re-litigated. B: require bitwise
everywhere — a correctly-rounded `exp` port, and `atan` stays blocked by LGPL
(OI-47).

**DEC-8 — where do the IMAGING types live?** E · CORE-18 counts three IMAGING
headers; four types sit inside `src/kernel/on_disc_imzml_experiment.rs` and
`src/format/imzml_file.rs`. A: move them to their own module and port the class
tests (`IonImage_test.cpp`, `MSImagingGeometry_test.cpp`, …) — touches two
landed modules. B: keep them where they are and port the class tests there —
cheapest, but `OnDiscImzMLExperiment` then owns types it does not declare.
Option "leave them ad hoc" is gone: M1 closes every core header.

**DEC-9 — give the gates a way to run a repository checker on a node?**
A: extend `openms-ci-sweep.sh` to run a job's non-`cargo` lines — shared
tooling every lane uses. B: leave `check_schema_feature_graph.py` CI-only — a
pre-push sweep cannot catch its class of failure.

**DEC-11 — read as the source does in every tool?** Only `PeakPickerHiRes`
passes `mzml::ReadOptions::source()`; the other tools read 32-bit minute time
arrays at full `f64` (`CPP-306`), and for `MzMLSplitter` that lands in written
output. A: a tool-framework default (library defaults stay strict) — one change
moves seven tools' output at once. B: per-tool opt-in — seven chances to
forget. C: leave it — a measured written-output divergence.

**DEC-12 — `ToolContext::in_thread_pool` at one worker:** A: a one-worker fast
path — a stray `par_iter` on the global pool becomes unbounded. B: let a tool
that scopes its own pool skip the outer one — unenforced. C: leave it — 0.68 s
on every large single-threaded run. Coupled to OI-103.

**DEC-13 — admit `rustix` for `access(R_OK)`?** A: yes, as a dependency —
`unsafe` inside a dependency, cache it for `--offline`. B: no — `-in /dev/null`
keeps exit 2 and a false message where C++ exits 4
(`docs/THIRD_PARTY_CRATE_DECISIONS.md:47`).

**DEC-14 — the allocator.** A: set glibc's three `MALLOC_*` knobs in the
launcher, all or none — glibc-only, −0.79 s. B: nothing — the lane's own
recommendation. Forbidden either way: `libc`, or relaxing
`#![forbid(unsafe_code)]`.

**DEC-15 — a scan-level non-FAIMS ion mobility of `-1`** loads but the writer
refuses it. A: refuse on read — rejects a file C++ reads. B: read it as unset —
drops a value silently. C: document it as a native difference of
`docs/MZML_MOBILITY_SUPPORT.md` — the round trip stays impossible, with a name.

**DEC-16 — budget node time for full-size FeatureFinderCentroided runs?**
Covers OI-101, OI-102, OI-116. A: about an hour per cell on a reserved node,
both sides, plus the never-run `FFC_max_iterations_*` sweep. B: accept the
4,000-spectrum subset and say so — the tool's acceptance stays qualified.

**DEC-17 — may a resource ceiling refuse data the reference build completes?**
One ruling settles OI-71, OI-93, OI-107 and A4's SRM conversion ceiling (about
430,000 SRM points), and OI-76 is already being fixed the size-derived way.
A: no — every default ceiling is size-derived, as the mzML reader's and the
window mower's (`2f173d9`) are; a ceiling then refuses only a pathological
ratio, and OI-71's exit code means "pathological input". B: yes, with each
bound documented and tested — named, permanent refusals of real data.

**DEC-18 — process rules.** (a) Does an integrator request naming an
integrator-owned file get its own checklist entry? (b) Must a new lane-owned
`docs/*.md` be linked from a shared document in the same pass? Nothing measures
Markdown reachability today. (c) The S1 SQLite review that
`docs/PORTING_WAVES.md` records as "pending authorization": authorize or
withdraw. Costs of not deciding: stale integrator paragraphs, unreachable
support documents, and a review that neither happens nor is cancelled.

### 5.3 Owed, new

**DEC-19 — may the public `_sorted` statistics refuse input the source
computes?** `median_sorted`, `quantile1st_sorted`, `quantile3rd_sorted` and
`quantile` return `UnsortedData` for a NaN or an unsorted slice
(`src/math/statistic_functions.rs:511-516`, `:670-675`, `:728`). The C++ does
not check a caller's "already sorted"; it reads positions, which is in bounds
and deterministic, so D1's letter says reproduce. A: follow the source — the
functions return what C++ returns for any input, including a garbage answer for
a false claim. B: keep the check for an unsorted slice, drop it for a NaN
(which "sorted" cannot describe) — the port stays stricter than the source, with
a named native difference. C: keep today's refusal — the last place where
permutation indeterminacy alone justifies a refusal, which D16 was taken to end.

**DEC-20 — does a rustdoc intra-doc link count as a module edge?**
`tools/check_module_cycles.py` counts every `crate::` path in the text,
comments included (OI-69). A doc link *is* a dependency once the workspace is
split (STRUCT-02): rustdoc cannot resolve it across crates without one. A plain
`//` comment is not. A: count doc links, ignore plain comments — a small
tokenizer and new gate tests. B: ignore all comments — cycles through doc links
reappear at the split. C: keep counting everything — the gate keeps failing on
prose.

**DEC-21 — work that no M1 requirement owns.** M1 is header-scoped, so two kinds
of known defect have no requirement: tool-level divergences in the nine
validated workflows (OI-85, OI-89, OI-96, OI-97, OI-98, OI-108) and defects in
headers already `complete` (OI-99's panic, OI-100, the order-substituting
divergences of section 3.2). A: add a Phase A requirement (say DEBT-06) and run
them as one lane — M1 grows by a lane, and the validated count stays honest.
B: defer to M2's TOOLS-01 — cheapest, but a validated workflow carries a known
divergence until then, and OI-99 is a live panic on public input. C: un-validate
a workflow while its divergence stands — honest, and the count drops.
*Whatever is chosen, OI-99 is a standing-constraint violation and should not
wait for it.*

---

## 6. By requirement

| Requirement | Items |
|---|---|
| A · DEBT-01 | OI-92 |
| A · DEBT-02 | OI-36 |
| A · DEBT-03 (this document; CI and lead-file follow-ups) | OI-56, OI-87 |
| A · DEBT-04 (lane A3a) | OI-01, OI-07, OI-08, B1-B7 |
| A · DEBT-05 (lane A4) | section 4 |
| A · EVID-01 | OI-113 |
| A · EVID-03 | OI-57, OI-111 |
| A · EVID-04 (lane A4) | OI-60 |
| B · BLOCK-01 | OI-18, OI-38, OI-88, OI-115 |
| B · BLOCK-02 | OI-81 (LogStream) |
| B · BLOCK-03 | OI-40, OI-76, OI-82, OI-105 |
| B · BLOCK-05 | OI-75 |
| B · BLOCK-06 | OI-23, OI-43, OI-73, OI-88, OI-90, OI-94, OI-109 |
| B · SHARED-01 (lane B6) | OI-44 |
| B · SHARED-02 | OI-115 (ZIP) |
| B · CLI-01 | OI-90 |
| C · CORE-01 | OI-09, OI-11, OI-12, OI-13, OI-20 (F3a), OI-21, OI-24, OI-26, OI-27, OI-28 (F2), OI-29 (F2), OI-31, OI-32, OI-35, OI-55, OI-72, OI-81, OI-95, OI-115 (F4) |
| C · CORE-02 | OI-36, OI-50 |
| C · CORE-05 | OI-33 |
| D · CORE-06 | OI-84 |
| D · CORE-08 | OI-34 |
| D · CORE-10 | OI-27 (metrics) |
| D · CORE-11 | OI-82, OI-103 |
| E · CORE-15 | OI-49, OI-50, OI-78, section 3.2 (smoothing note) |
| E · CORE-18 | OI-03, OI-04, OI-05, OI-26, OI-35, OI-39, OI-41, OI-77, OI-81, OI-83 |
| F · STRUCT-01 | OI-70 |
| G · UPSTREAM-01..04 | DEC-10; `CPP-337` (OI-14), `OPENMS-PERCIN-007` (OI-22), `CPP-338` (OI-79) |
