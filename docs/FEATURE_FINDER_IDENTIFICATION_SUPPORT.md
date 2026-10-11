# FeatureFinderIdentification and FeatureFinderMetaboIdent

Ports, at core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`, the two targeted MS1 feature
finders of `FEATUREFINDER/`. Both drive the OpenSWATH chromatogram extraction, the transition
group picker and the MRM feature scorer (`MRMFeatureFinderScoring`) and then annotate, select and
quantify the picked candidates.

| Header | Rust |
| --- | --- |
| `FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h` (+ `.cpp`) | [`src/analysis/featurefinder/feature_finder_identification.rs`](../src/analysis/featurefinder/feature_finder_identification.rs) |
| `FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h` (+ `.cpp`) | [`src/analysis/featurefinder/feature_finder_metabo_ident.rs`](../src/analysis/featurefinder/feature_finder_metabo_ident.rs) |
| (the steps the two share in the source's two bodies) | [`src/analysis/featurefinder/ms1_targeted.rs`](../src/analysis/featurefinder/ms1_targeted.rs) (not public API) |

Tests: [`tests/feature_finder_identification.rs`](../tests/feature_finder_identification.rs) (class
test sections, parameters, behaviour on small synthetic MS1 data; no feature),
[`tests/feature_finder_identification_topp.rs`](../tests/feature_finder_identification_topp.rs)
(the retained upstream TOPP outputs; `mzml` + `featurexml` + `idxml`) and
[`tests/feature_finder_identification_oracle.rs`](../tests/feature_finder_identification_oracle.rs)
(executed Release build, 29 cases; same features). Support code
[`tests/support/ffid_inputs.rs`](../tests/support/ffid_inputs.rs). Fixtures:
[`tests/data/feature_finder_identification/`](../tests/data/feature_finder_identification/).
Manifest:
[`tests/data/feature_finder_identification_provenance.json`](../tests/data/feature_finder_identification_provenance.json).

The source has no OpenMP pragma in either file, so both are serial here. No crate was added.

One addition outside the two files: `AASequence::formula_for_source`
([`src/chemistry/sequence.rs`](../src/chemistry/sequence.rs)) and the one line of
[`openswath_scoring.rs`](../src/analysis/openswath/openswath_scoring.rs) that uses it
(see *Preserved source conventions*, "an anonymous modification has no atoms"), and
`MRMFeatureFinderScoring::set_shared_unique_id_generator`
([`src/analysis/openswath/mrm_feature_finder_scoring.rs`](../src/analysis/openswath/mrm_feature_finder_scoring.rs)).

## API mapping

### `FeatureFinderIdentificationAlgorithm`

| Source | Rust |
| --- | --- |
| `class ... : public DefaultParamHandler`; default constructor installing the defaults | `FeatureFinderIdentificationAlgorithm::new()`, `defaults()`, `from_param`, `from_param_with_warnings`, `set_parameters`, `parameters()` |
| `run(PeptideIdentificationList peptides, proteins, FeatureMap& features, seeds, spectra_file)` | `run(Vec<PeptideIdentification>, &[ProteinIdentification], &mut FeatureMap, &FeatureMap, &str) -> Result<()>` (staged: `features` is unchanged on error) |
| `runOnCandidates(FeatureMap&)` | `run_on_candidates(&mut FeatureMap) -> Result<()>` |
| `getMSData()` const and non-const | `ms_data()`, `ms_data_mut()` |
| `setMSData(const PeakMap&)`, `setMSData(PeakMap&&)` | `set_ms_data(MSExperiment)` (both overloads: keeps only MS1) |
| `getChromatograms()` const and non-const | `chromatograms()` (empty after a run, as the source) |
| `getProgressLogger()` const and non-const | `progress_logger()`, `progress_logger_mut()` |
| `getLibrary()` const and non-const | `library()` (the accumulated `output_library_`; empty for multi-FAIMS data) |
| protected `n_peptides_`, `statistics_` (reached by the class test through a subclass) | `peptide_count()`, `set_peptide_count()`, `statistics(&FeatureMap) -> FeatureFinderIdentificationStatistics` |
| protected `rt_window_`, `min_peak_width_`, `global_im_stats_` | `rt_window()`, `min_peak_width()`, `global_im_stats()` |
| the log lines (`OPENMS_LOG_INFO`, `_WARN`, `_ERROR`) | `report() -> &FeatureFinderIdentificationReport` (numbers and `messages`) |
| process-wide `UniqueIdGenerator` | `set_unique_id_generator(UniqueIdGenerator)`; shared with the FAIMS group instances |
| `typedef MassTrace(s)`, `RTMap`, `ChargeMap`, `PeptideMap`, `PeptideRefRTMap` | private type aliases over indices into the run's identification list |
| `struct IMStats` | `IMStats` (`median`, `min`, `max`, all `-1` by default) |
| `struct RTRegion` | private `RtRegion` |
| `FeatureFilterQuality`, `feature_filter_quality_` | not ported: declared and never used in the source |
| `FeatureFilterPeptides`, `PeptideCompare`, `FeatureCompare` | closures `filter_features`, `peptide_compare`, `feature_compare` |
| `updateMembers_`, parameter-derived members | `Settings::from_param` (staged by `set_parameters_with`) |
| `generateTransitions_`, `addPeptideRT_`, `getRTRegions_`, `getRTRegionIMStats_`, `calculateGlobalIMStats_`, `createAssayLibrary_` | `generate_transitions`, `add_peptide_rt`, `rt_regions`, `region_im_stats`, `calculate_global_im_stats`, `create_assay_library` |
| `annotateFeaturesFinalizeAssay_`, `annotateFeatures_`, `ensureConvexHulls_`, `postProcess_`, `filterFeatures_` | `finalize_assay`, `annotate_features`, `ms1_targeted::ensure_convex_hulls`, `post_process`, `filter_features` |
| `initializeFeatureFinder_`, `calculateRTWindow_`, `calculateRTBounds_`, `removeSeedPseudoIDs_`, `isSeedPseudoHit_` | `initialize_feature_finder`, `calculate_rt_window`, `calculate_rt_bounds`, `remove_seed_pseudo_ids`, `is_seed_pseudo_hit` |
| `addPeptideToMap_`, `addSeeds_`, `addOffsetPeptides_`, `runSingleGroup_`, `chunk_` | `add_peptide_to_map`, `add_seeds`, `add_offset_peptides`, `run_single_group`, `chunk_ranges` |
| `static getIMUnitFromExperiment` | `first_im_unit` |

### `FeatureFinderAlgorithmMetaboIdent`

| Source | Rust |
| --- | --- |
| `class ... : public DefaultParamHandler` | `FeatureFinderAlgorithmMetaboIdent::new()`, `defaults()`, `from_param`, `from_param_with_warnings`, `set_parameters`, `parameters()` |
| `struct FeatureFinderMetaboIdentCompound` (constructor, nine getters) | `FeatureFinderMetaboIdentCompound::new(...)`, `name()`, `formula()`, `mass()`, `charges()`, `rts()`, `rt_ranges()`, `isotope_distribution()`, `ion_mobilities()`, `adduct()` |
| `run(const std::vector<...Compound>&, FeatureMap&, spectra_file)` | `run(&[Compound], &mut FeatureMap, &str) -> Result<()>` (staged) |
| `getMSData()`, `setMSData(const&)`, `setMSData(&&)` | `ms_data()`, `ms_data_mut()`, `set_ms_data` (keeps only MS1) |
| `getChromatograms()` | `chromatograms()` (kept after the run, unlike the other class) |
| `getLibrary()`, `getTransformations()`, `getNShared()` | `library()`, `transformations()`, `n_shared()` |
| `static prettyPrintCompound(const Compound&)` | `FeatureFinderAlgorithmMetaboIdent::pretty_print_compound(&Compound) -> Result<String>` |
| `struct MassTraceBounds`, `FeatureBoundsMap`, `MassTrace(s)`, `FeatureGroup` | not ported: declared and never used in the source |
| `FeatureFilterQuality feature_filter_` (marker `FFMetId_remove`) | the retained-test in `select_features_from_candidates` |
| `FeatureCompare` | `feature_compare` |
| `extractTransformations_`, `addTargetToLibrary_`, `addTargetAnnotations_`, `addTargetRT_`, `calculateMZ_`, `generateTransitions_`, `annotateFeatures_`, `ensureConvexHulls_`, `selectFeaturesFromCandidates_`, `runSingleGroup_` | `extract_transformations`, `add_target_to_library`, `add_target_annotations`, `add_target_rt`, `calculate_mz`, `generate_transitions`, `annotate_features`, `ms1_targeted::ensure_convex_hulls`, `select_features_from_candidates`, `run_single_group` |
| `updateMembers_`, `iso_gen_` (`CoarseIsotopePatternGenerator`) | `Settings::from_param`; the generator is built per target with `n_isotopes` bins |
| `debug_level_` | read in the source and used nowhere in this class; not kept |

## Preserved source conventions

* **Targets are indices.** The source maps sequence, charge and RT to a *pointer* to the
  identification, stored in the caller's list. The port keeps the list in the run and maps to
  indices, so seeds and offset peptides appended later do not invalidate anything. The map keys
  are ordered as `AASequence::operator<` orders them (length first, then the terminals, then
  residue letters and modifications), not as the native `Ord` does, because the order is the
  order of the batches, of the library and of the unique IDs.
* **Batches.** `chunk_`: `total / batch_size` batches, the last one extended to the end; one
  batch for fewer, also for none. The chromatograms of a batch are cleared after the
  features are picked, so `chromatograms()` is always empty after a run (the TOPP tool's
  `-chrom_out` therefore writes an empty file; candidate C++ defect, below).
* **Per-batch state.** The seed counter of the library builder restarts for every batch; the
  `IMStats` and isotope probabilities accumulate over the run.
* **One run per instance.** `rt_window_` and `min_peak_width_` are updated in place by a run
  (a relative minimum width below one is multiplied by the peak width each time a run starts), the
  unassigned decoys and the IM statistics are not reset. A second run on one instance sees these.
* **Seeds and offset peptides** are identifications of the pseudo sequence `XXX[10000]` with a
  `SeedFeatureID`; the library gives them an empty sequence and a group label of their own. The
  pseudo hits are removed from the output; features whose identifications are all offset peptides
  get `OffsetPeptide`, and offset peptides are exempt from the seed-apex filter.
* **IDs are mapped by RT only.** Identifications between `leftWidth` and `rightWidth`, widened by
  `detect:mapping_tolerance` (relative below one), belong to a feature; the feature with most
  identifications is "positive", ties broken by intensity; other features are "ambiguous" or
  "negative". `n_total_ids` counts the identifications of the region.
* **An anonymous modification has no atoms.** `YIC[47.021464]DK` carries a mass the source adds to
  its database as an "unknown modification"; its `getFormula` leaves it out and only the weight
  carries the shift. The native `AASequence::formula` refuses such a sequence, so the isotope
  pattern here (and the MS1 isotope score of OpenSWATH) uses `formula_for_source`.
* **Isotope arithmetic is `float`.** The distributions are computed, trimmed (`trim_left_source`,
  `trim_right`) and renormalised in the source's binary32 arithmetic
  (`ProbabilityPrecision::SourceSingle`); library intensities and `isotope_probability` are those
  `float` values widened. `float(counter)` and `float(mass)`/`float(rt)` in the ids are kept
  (`StringUtils::toStr(float)`).
* **MetaboIdent keeps stale target state.** The target object is reused across charges and
  RTs: `expected_im` and `adduct` persist when a later RT has none; a failed adduct parse leaves
  the charge on the formula for the following charges.
* **MetaboIdent's selection.** The candidates of an assay are consecutive; the one whose range
  contains the expected RT wins, else the nearest, ties by intensity; the loser is marked
  `FFMetId_remove` and removed. `rt_deviation` is `expected - found`. Overlapping features are
  then resolved with `FeatureOverlapFilter::filter` at trace level, the comparator comparing
  `(|rt_deviation|, truncated intensity)` and the overlap callback annotating `label`,
  `alt_PeptideRef` and `overlap_removed` (the latter only when the truncated intensities are
  equal, as written).
* **An empty candidate map throws.** `FeatureOverlapFilter::filter` reads the m/z range of the
  map, so a run without any candidate ends in `Exception::InvalidRange` ("Empty or uninitalized
  range object"), executed on the Release build (`mi_no_candidates`). The port returns
  `Error::InvalidRange` and leaves the features unchanged.

## Native differences

| Source | Here | Why |
| --- | --- | --- |
| log streams | `report()` | no module of this crate logs |
| `run` modifies `features` as it goes | staged copy, committed on success | an error must not leave a half-updated map |
| `min_psm_cutoff` is a string parameter read with `double(param.getValue(...))`, which reads the `double` member of the value union: the address of the string | the number the text spells; a non-number is an error | undefined behaviour in the source (candidate defect) |
| `peptide.sort()` on NaN scores | `Error::InvalidValue` | the native sort validates |
| identification without hits in `addOffsetPeptides_` | skipped | the source reads `getHits()[0]` unchecked |
| hit without charge | `Error::InvalidValue` | the m/z would be infinite |
| `candidates_out` as FeatureParquet | `Error::Unsupported` before anything runs | FeatureParquet writing is a backlog item (dependency: the Arrow/Parquet FORMAT port) |
| `debug >= 668` dumps `debug.traml` | the dump is not written, a message says so | the TraML writer (`FORMAT/TraMLFile.h`) is not ported; the run is unaffected |
| `MSExperiment::updateRanges` before `getMinMobility` | ranges are computed on demand | no range state in this port |
| `setModificationByDiffMonoMass` (`CalcMass` meta value) | the best modification within 0.002 Da of the registry, else an anonymous mass tag | the unknown-modification database entry is not mutated |
| unique IDs from the process-wide generator | one `Arc<Mutex<UniqueIdGenerator>>` shared by the algorithm, its finder and the FAIMS groups | the stream is the source's when the generator is seeded as `-test` |
| `std::sort` of features | stable sort | ties of (`PeptideRef`, RT) are not reproduced by an unstable sort anyway |
| `unassignedIDs_` moved from, the husks stay | drained | the husks are empty identifications |
| `ProgressLogger::NONE` forced on `feat_finder_` | no progress of the finder | the finder has no progress output here |

## Checked boundaries

* `run` refuses without spectra, with more than `MAX_TARGETS` (2^26) identifications and seeds
  (offset peptides counted), with `candidates_out` that cannot be written, before any state changes.
* MetaboIdent refuses more than `MAX_LIBRARY_ENTRIES` (2^24) library entries before building them.
* A target that violates the table contract is skipped with a message and counted
  (`targets_skipped`); a charge of 0 skips that charge only.
* The report keeps at most 10,000 lines; the rest is counted.

## Evidence

| Evidence | Tier | Where |
| --- | --- | --- |
| `FeatureFinderIdentificationAlgorithm_test.cpp`, all three sections | 3 | `feature_finder_identification.rs` (`class_test_*`) |
| `TOPP_FeatureFinderIdentification_1`, `_3`, `_5` (both outputs), `_summary_count`; `TOPP_FeatureFinderMetaboIdent_1`; reprocessing of the retained candidates | 1 (upstream retained outputs, `FuzzyDiff` 1.01 / 0.01) | `feature_finder_identification_topp.rs` |
| 29 cases of the Release build (below), values bit for bit | 1 | `feature_finder_identification_oracle.rs` |
| parameters, RT window rule, relative minimum width, decoys, cutoff direction, seed overlap rule, offset peptides, determinism, table contract, adduct defaults | 4 | `feature_finder_identification.rs` |

The oracle driver is `../oracle/FF3-ffid/driver.cpp`, run on ibminode06 against the Linux x86-64
Release install `openms4-release-bc9cc12-c19e494-174b576`, twice in fresh processes with identical
standard output (the written `result.featureXML` differs between the runs in the identifiers
of subordinate features, which are drawn from an unseeded thread-local stream and are not
compared). The identification run identifiers of the loaded idXML (`<engine>_<date>_<random>`) are
not compared either. Cases:

| Case | What it adds to the TOPP tests |
| --- | --- |
| `ffid_none`, `ffid_symmetric`, `ffid_asymmetric`, `ffid_batch10_candidates` | the TOPP cases in full precision; the asymmetric (EGH) model; candidates file |
| `ffid_mapping_tolerance`, `ffid_pmin`, `ffid_ppm` | relative mapping tolerance with a set RT window; `isotope_pmin`; a ppm window with three isotopes |
| `ffid_seeds`, `ffid_seeds_offset`, `ffid_seeds_apex_off` | seeds (overlapping and new); offset peptides; the seed-apex filter off |
| `ffid_decoys`, `ffid_decoys_quantified`, `ffid_psm_cutoff_string`, `ffid_calcmass` | decoys; the (defective) string cutoff; the `CalcMass` shift with a known and an unknown mass |
| `ffid_faims`, `ffid_faims_ids`, `ffid_im`, `ffid_im_off` | FAIMS split, merge and per-voltage identifications; per-peak ion mobility with and without the window |
| `mi_tsv`, `mi_none_candidates`, `mi_asymmetric`, `mi_pmin`, `mi_rt_window_default` | the TOPP table; no model with candidates; asymmetric; `isotope_pmin`; the default RT window |
| `mi_variants`, `mi_variants_none` | charges 1, 2 and -1; adducts (`M+Na;1+`, `M+K;1+`, an unparsable one); mass without formula; isotope distribution given; RT ranges; overlapping targets; contract violations |
| `mi_faims`, `mi_im`, `mi_im_off` | FAIMS; ion mobility per target with the window; the window off |
| `mi_no_candidates` | the empty-range exception |

`ffid_psm_cutoff_string` documents the defect: with `min_psm_cutoff` set to a string, the Release
build treats the cutoff as a number near zero and sets nearly every identification aside. It is not
compared (`tests` assert the intended behaviour at tier 4); the record file keeps the outcome.

## Candidate C++ defects

1. `min_psm_cutoff`: `double(ParamValue)` on a string value reads the union's `double` member (the
   address of the string). Undefined behaviour; a cutoff other than `none` acts as about `+0`.
2. `FeatureFinderIdentificationAlgorithm::run` clears `chrom_data_` after every batch, so
   `getChromatograms()` and the tool's `-chrom_out` are always empty (the member is documented as
   "accumulated chromatograms").
3. `FeatureFinderIdentification` on a second `run`: `unassignedIDs_` is not cleared and the moved-from
   elements are re-added; `min_peak_width_` is multiplied again.
4. `FeatureFinderAlgorithmMetaboIdent::run` throws `InvalidRange` ("Empty or uninitalized range
   object. Did you forget to call updateRanges()?") when no candidate feature was found, from
   `FeatureOverlapFilter::filter` on an empty map.
5. `addTargetToLibrary_` keeps `expected_im`/`adduct` meta values and the formula charge across
   iterations of the reused target object (listed above).
6. `addTargetAnnotations_`: `library_.getCompounds().size() - found_refs.size()` is an unsigned
   subtraction that wraps when more references than compounds were found (duplicate targets).
7. `annotateFeatures_` ends the feature loop with `ref_rt_map[peptide_ref]` for a possibly empty
   `peptide_ref`, creating an empty map entry (harmless).

## Findings in other lanes' code (not fixed here)

* The native idXML reader refuses the attribute `prot` that the source ignores (two hits of
  `FeatureFinderIdentification_1_input.idXML`).
* The idXML and featureXML readers spend a flat 50 million work units, and every annotated sequence
  is charged a scan of the whole registry (about 7,000 modifications). A 23 KB document with 21
  `(Carbamidomethyl)` hits exhausts it ("peptide parent parsing resource limit exceeded"), so the
  retained input of `TOPP_FeatureFinderIdentification_1` cannot be read with the global registry. The
  tests read it against a one-modification registry (`tests/support/ffid_inputs.rs`).
* `AASequence::formula_for` refuses an anonymous mass modification (the source omits it): the MS1
  isotope score of `OpenSwathScoring` failed on the two such peptides of that input. Fixed through
  `formula_for_source`.
* The native featureXML writer requires `ProteinIdentification::date_time`, where the source writes
  the zero date of an unset `DateTime`; `FeatureFinderMetaboIdent` sets `0000-00-00T00:00:00` on
  the run it creates.
