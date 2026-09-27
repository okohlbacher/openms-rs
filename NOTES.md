# Lane A5b NOTES (restart checkpoint)

Scope list: `/Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/a5-lists.json` key `A5b` (37 headers):
DATASTRUCTURES/Matrix.h, IMAGING/IonImage.h, IMAGING/MSImagingGeometry.h, IMAGING/MSImagingRegion.h,
IONMOBILITY/IMDataArrayUtils.h, IONMOBILITY/IMDataConverter.h, KERNEL/MSSpectrum.h, KERNEL/OnDiscImzMLExperiment.h,
MATH/STATISTICS/GumbelDistributionFitter.h, MATH/STATISTICS/MultipleTesting.h, MATH/STATISTICS/PosteriorErrorProbabilityModel.h,
METADATA/ChromatogramSettings.h, METADATA/DataArrays.h, METADATA/ID/IdentificationData.h, METADATA/ID/IdentificationDataConverter.h,
METADATA/MetaInfo.h, METADATA/MetaInfoDescription.h, METADATA/MetaInfoInterface.h, METADATA/MetaInfoRegistry.h,
METADATA/PeptideHit.h, METADATA/Precursor.h, METADATA/Product.h, METADATA/ProteinHit.h, METADATA/Software.h,
METADATA/SpectrumLookup.h, METADATA/SpectrumMetaDataLookup.h, METADATA/SpectrumNativeIDParser.h, METADATA/SpectrumSettings.h,
METADATA/USI.h, ML/SVM/SimpleSVM.h, SYSTEM/BuildInfo.h, SYSTEM/JavaInfo.h, SYSTEM/PythonInfo.h, SYSTEM/RWrapper.h,
SYSTEM/SIMDe.h, SYSTEM/StopWatch.h, SYSTEM/SysInfo.h

Outputs: report -> main repo .planning/lane-reports/A5b.md ; ledger -> main repo .planning/ledger-fragments/wave4/A5b.json

## Findings so far
- C++ sources dumped to scratchpad cpp/ via scratchpad/dump_cpp.sh (git -C with a LITERAL path works; shell variables
  in git -C, or python heredocs mentioning git, are refused by the worktree guard).
- Ledger statuses allowed: complete | partial | native_equivalent (tools/core_sdk_coverage.py:234).
- JavaInfo/PythonInfo/RWrapper: every member mapped in docs/*_SUPPORT.md API tables; remaining 'partial' reasons are API
  shape differences (verbose param, FileContext) -> mark COMPLETE.
- Gumbel: git grep at bc9cc12 confirms GumbelDistributionFitter::fitWeighted (h:81) and GumbelDistributionFitResult::eval
  (h:51) have NO definition (only GumbelMaxLikelihoodFitter::fitWeighted exists, cpp:74; its test :175 calls the
  MaxLikelihood one). -> record 'complete' with the impossibility evidence (no 'impossible' status exists).
- MultipleTesting: only gap = Pi0Smoother seam (math cannot import processing). Close by a library adapter in processing
  (new file + pub mod line) implementing Pi0Smoother over processing::spline::smoothing::BSplineSmoothingSpline; switch
  tests/multiple_testing.rs's local struct to the library one.
- SysInfo/StopWatch/BuildInfo: partial for platform bindings (macOS/Windows); leave partial, name crates (rustix getrusage,
  memory-stats, sysinfo). cpu-time 1.0.0 already a target dep.
- fancy-regex =0.19.2 IS a dependency now -> SpectrumNativeIDParser/SpectrumLookup scopes saying 'no regex engine' are STALE.

## Next step
1. Dump ledger entries (docs/core-sdk-coverage.json) for the 37 headers into scratchpad.
2. Triage: cheap confirmations first (MultipleTesting, SysInfo, PythonInfo, RWrapper, JavaInfo, BuildInfo, StopWatch, SIMDe), Gumbel impossible.
3. Then USI.h, MetaInfoRegistry.h ports.
