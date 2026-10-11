# O3b notes (MRM feature finding, filtering, QC)

## State (commit 08cde132 and later)
- Merged prerequisites (O2 cd066dab-47, peak picking 97065d91-1) and ALSO FF1 (cd066dab-33) for EmgScoring (documented in report).
- ALL 9 headers have a Rust port + tests, green on spock via `openms-kim-gate.sh O3b test --locked --test <name>`:
  mrm_feature_qc (15), mrm_feature_filter (26), mrm_feature_selector (12), confidence_scoring (6), mrm_feature_finder_scoring (6, all 3 upstream sections pass).
- Modules: src/analysis/openswath/{mrm_feature_qc,mrm_feature_picker,mrm_feature_filter,mrm_feature_support(pub(crate)),mrm_feature_selector,mrm_batch_feature_selector,confidence_scoring,mrm_feature_finder_scoring}.rs, src/format/{mrm_feature_qc_file,mrm_feature_picker_file}.rs; data_access_helper.rs gained convert_targeted_exp/peptide/compound.
- Also touched: src/chemistry/decoy_random.rs (generic shuffle_slice), src/format/identification_xml.rs (accept the unset date 0000-00-00T00:00:00; upstream fixtures carry it), src/format/mod.rs, src/analysis/openswath.rs.
- Gaps: LP solver for MRMFeatureSelectorQMIP (LPWrapper GATED-01; PartitionSolver exact for Score only), TraMLFile not ported (tests build TargetedExperiment by hand), MobilogramParquetConsumer, serial per-feature loop (source: OpenMP).
- Helper scripts in scratchpad (git/other-repo commands are blocked in bash; use python scripts written with the Write tool).

## Remaining
1. Tier-1 oracle on ibminode06: ../oracle/O3b-mrm/{driver.cpp,remote_run.sh,run.sh,cases.txt}; then tests/mrm_oracle.rs + tests/data/mrm_oracle/*.
2. Support docs: docs/MRM_FEATURE_QC_SUPPORT.md, MRM_FEATURE_SELECTOR_SUPPORT.md, MRM_FEATURE_FINDER_SUPPORT.md (names are referenced from module docs), provenance json tests/data/mrm_*_provenance.json.
3. doc coverage (tools/check_doc_coverage.py) 100% for touched modules; clippy; fmt; full CI sweep (openms-ci-sweep.sh) in background.
4. Reports: .planning/lane-reports/O3b.md, .planning/ledger-fragments/wave4/O3b.json; delete NOTES.md in the last commit.
