# F1c notes (QPX lane)

Merged prerequisite worktree-wf_cd066dab-0a7-5 (F1a arrow core, fast-forward). No earlier incarnation branch.

C++ sources at pin copied to scratchpad (read via `git -C .../OpenMS4-tests/packages/core show bc9cc12:<path> > file`,
one plain git command per call; compound git forms are refused by the worktree guard):
- QPXValueValidation.h/.cpp (136/722), test 384
- QPXIdentity.h/.cpp (208/453), test 244
- QPXFile.h/.cpp (146?/2158), test 1889
- QPXCollectionExport.h/.cpp (146/237), test 538

Existing Rust: src/format/arrow_io_helpers.rs has QPX naming helpers, qpx_file_metadata, and a
`validate_qpx_values_for_write` stub that refuses psm/feature/pg with Error::Unsupported (F1a gap).
Schemas QPXPSMSchema/QPXFeatureSchema/QPXPgSchema in src/format/arrow_schema_registry.rs.

Order: 1 QPXValueValidation (+ wire into arrow_io_helpers gate), 2 QPXIdentity, 3 QPXFile, 4 QPXCollectionExport.

## Findings
- QPXValueValidation: needs only schema registry + qpx_is_canonical_intensity_label (exists). Tables are
  RecordBatch (single chunk) so combinedColumn is trivial. Replace the Unsupported stub in
  arrow_io_helpers::validate_qpx_values_for_write with the real validator (need C++ ArrowIOHelpers.cpp
  validateQPXValuesForWrite message text).
- QPXIdentity: pure; class test values come from the qpx Python reference (independent, tier upstream).
  BLAKE2b-64 -> crate `blake2` (=0.10.6, RustCrypto, Blake2bVar::new(8)); digest 0.10 already in lock.
  Plan: optional dep enabled by `arrow` feature. Float repr: Rust `{:e}` on f64 gives shortest digits.
- QPXFile deps: IDScoreSwitcherAlgorithm (partial port in src/analysis/scores.rs: find_peptide_score,
  ScoreType), Scores::isKnownScoreType (check), IdentifierMSRunMapper (src/identification/run_mapping.rs),
  ProForma from_aa_sequence (src/chemistry/proforma_conversion.rs), PSMSchema/QPXPSMSchema (registry).
  format<->analysis already mutually dependent, so format->analysis::scores costs no new cycle.
- QPXCollectionExport: Transaction is self-contained (port fully). requireExportable needs
  ConsensusMapArrowExport::{requireUnambiguousIdentities,requireResolvableIdRuns},
  ProteinGroupArrowExport_impl (QuantificationUnits, sampleAbundances, ...), ExperimentalDesign -> gap.
- Worktree guard: run git commands alone (no chaining, no 'git' inside heredocs in the same call).

## Done
- QPXIdentity + QPXValueValidation committed (tests pass on dax: qpx_identity 7, qpx_value_validation 17,
  arrow_parquet 23, lib 368). Cargo.lock fetched from dax (ssh ... cat /scratch/kohlbach/agents/F1c/Cargo.lock).
- IdentifierMSRunMapper::validate_merge_index now uses the C++ messages + MissingInformation (tests/run_mapping.rs updated).

## QPXFile design decisions
- rank: Rust PeptideHit.rank field <-> C++ "rank" metavalue (only when != 0).
- Warnings returned (QpxOutput{value,warnings}); streaming writes per-batch row groups via
  ArrowWriter::into_serialized_writer over Vec<u8> drained into path_io::write_plain.
- operator double on a string DataValue is UB-ish garbage in C++ -> Rust refuses (defect candidate).
- isKnownScoreType = ScoreType::from_normalized_name; higher_better = ScoreType::parse(key).ok().

## Oracle (done)
- ../oracle/F1c-qpx (driver.cpp, node.sh, run.sh; `bash run.sh push` then `pull`), ibminode06 Release,
  run twice identical records. Outputs copied to tests/data/qpx/ (records, import_records, 5 parquet).
- Findings: C++ import resolves mass deltas only via process-global ModificationsDB mutated by
  AASequence::fromString (fresh process = Rust behaviour); AASequence mz differs 1-2 ULP (sequence port);
  ProForma parser error text differs (parser port); mass-tag mod name is "" in C++.

## Next step
1. Get tests/qpx_file.rs green (calculated_mz tolerance 1e-15 declared; unordered metadata lists sorted).
2. docs/QPX_SUPPORT.md, tests/data/qpx_provenance.json, THIRD_PARTY blake2 row, ARROW_PARQUET_SUPPORT update.
3. MSRV + clippy + doc + ci sweep; report + ledger fragment; delete NOTES.md in final commit.
