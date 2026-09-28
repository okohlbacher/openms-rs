# Lane I4 notes (IDMapper + metabolomics ID)

## Status
- GREEN + committed: IDMapper (21 tests), IonIdentityMolecularNetworking (8), SiriusFragmentAnnotation (6);
  clippy all-features/all-targets clean, MSRV 1.85 green, doc green, no-default check green, doctests green.
- MetaboliteSpectralMatching: code + tests written (partial: merge_spectra=true needs SpectraMerger); gating t7.
- docs/ID_MAPPER_SUPPORT.md, tests/data/*_provenance.json (4) written; generators in scratchpad/I4/{provenance,ledger}.py.
- Ledger fragment written: .planning/ledger-fragments/wave4/I4.json (4 entries).
- Scratchpad is SHARED with other lanes: use scratchpad/I4/ only.
- Worktree guard: one plain git command per Bash call; no shell variables as paths in git commands.
- Module cycles: analysis must NOT import crate::system; analysis->data_structures edge recorded.
- Doc coverage: new modules 100%; module roots need `///` on each `pub mod`; do NOT run --write (lead's file).
- core_sdk_coverage.py reports docs/core-sdk-coverage.json stale -> lead regenerates.
- Blocked headers: SiriusMSConverter + SiriusExportAlgorithm (FeatureMapping.h unmapped), AccurateMassSearchEngine
  (MzTab/MzTabM export + IdentificationData + HMDB mapping files), FIAMSDataProcessor/Scheduler (AccurateMassSearchEngine,
  SpectrumAddition, SignalToNoiseEstimatorMedianRapid).

## Next step
1. Confirm t7 (MSM) green; clippy + MSRV for tests/metabolite_spectral_matching.rs.
2. Write the lane report (.planning/lane-reports/I4.md) with slices and pass counts; delete NOTES.md in final commit.
