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

## Next step
1. Port QPXIdentity (src/format/qpx_identity.rs) + blake2 dep; commit.
2. Port QPXValueValidation (src/format/qpx_value_validation.rs), wire into arrow_io_helpers; commit.
3. QPXFile (src/format/qpx_file.rs). 4. QPXCollectionExport::Transaction.
