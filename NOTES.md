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
(none yet)

## Next step
1. Dump ledger entries (docs/core-sdk-coverage.json) for the 37 headers into scratchpad.
2. Triage: cheap confirmations first (MultipleTesting, SysInfo, PythonInfo, RWrapper, JavaInfo, BuildInfo, StopWatch, SIMDe), Gumbel impossible.
3. Then USI.h, MetaInfoRegistry.h ports.
