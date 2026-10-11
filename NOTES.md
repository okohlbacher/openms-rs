# FF3 notes (lane FeatureFinderIdentification + MetaboIdent)

Branch worktree-wf_97065d91-48e-8. Prereq branches merged (8fd08226): FF1 models, O3b MRM scoring, extraction.
No earlier FF3 incarnation existed (no NOTES.md on any prerequisite).

## Sources (pinned checkout, read here)
/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12/
- src/openms/{include/OpenMS,source}/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.{h,cpp}
- ... /FeatureFinderAlgorithmMetaboIdent.{h,cpp}
- class test (only FFId has one, 3 tiny sections): src/tests/class_tests/openms/source/FeatureFinderIdentificationAlgorithm_test.cpp
- TOPP: `git -C ../OpenMS4-tests/packages/topp show 174b576:src/FeatureFinderIdentification.cpp` (and MetaboIdent; TSV parse in readTargets_)
- Retained outputs (tier 1): test-data topp/CMakeLists.txt:434-459 (FFId_1 -model:type none, _3 symmetric, _5 batch_size 10 + candidates_out)
  and :4267-4270 (MetaboIdent_1). Inputs/outputs in ../OpenMS4-tests/packages/test-data/topp/FeatureFinder{Identification,MetaboIdent}_*.
  FFId_4 is disabled upstream (##). Only whitelist: "spectra_data" "featureMap" (+ "feature id" for _5, "id=" for MetaboIdent); TOPP -test sets unique-id seed.

## Rust dependencies already present (all verified)
- openswath::mrm_feature_finder_scoring::MRMFeatureFinderScoring (parameters/set_parameters/set_ms1_map/pick_experiment(chroms,&mut fm,&TargetedExperiment,&TransformationDescription,&MSExperiment)); `ids: Arc<Mutex<UniqueIdGenerator>>` private -> I add set_shared_unique_id_generator
- openswath::chromatogram_extractor::{ChromatogramExtractor, prepare_coordinates, return_chromatogram}, data_access_helper::convert_targeted_exp, SimpleOpenMSSpectraFactory::spectrum_access_openms_ptr
- featurefinder::elution_model_fitter::ElutionModelFitter (defaults/to_param/set_parameters/fit_elution_models(&mut FeatureMap)->Result<Report>)
- kernel::im_data_converter::ImDataConverter::split_by_faims_cv(&mut MSExperiment), kernel::faims_helper::FaimsHelper::filter_peptides_by_faims_cv(peps, cv, 0.01)
- processing::feature_overlap_filter::FeatureOverlapFilter::{filter, merge_faims_features}
- chemistry::isotopes (use ProbabilityPrecision::SourceSingle, trim_left_source, renormalize_with), AdductInfo::parse + mz, EmpiricalFormula::parse
- FeatureMap::set_primary_ms_run_path_from_experiment; TransformationDescription::set_data_points; FileHandler::store_feature_map (feature featurexml)
- tests/support/{decoded_compare,fuzzy_string_comparator,release_runs}.rs

## Design
- peptide_map: BTreeMap<SourceSeqKey (C++ AASequence::operator< : size, n-term, letters+mods, c-term), BTreeMap<i32, Vec<(f64 rt, usize idx)>>>; idx -> local Vec<PeptideIdentification>.
- ensureUniqueId on a FeatureMap only draws the MAP's own id (UniqueIdInterface), not features.
- log lines -> run report struct (crate does not log).

## Status
- [x] prerequisites merged, C++ read, Rust deps surveyed
- [ ] feature_finder_identification.rs
- [ ] feature_finder_metabo_ident.rs
- [ ] tests (class test sections, oracle replay of TOPP cases, own Release oracle on ibminode06)
- [ ] support doc, provenance, ledger fragment, lane report, CI slices (via report only; ci_matrix is lead's)

## Next step
1. Add `set_shared_unique_id_generator` to mrm_feature_finder_scoring.rs.
2. Write src/analysis/featurefinder/ffid_common.rs (shared helpers), feature_finder_identification.rs.
