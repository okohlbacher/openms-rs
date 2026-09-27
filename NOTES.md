# Lane F5 NOTES (restart-safety; delete in final commit)

Branch: worktree-wf_ee239125-5be-15, base f0990c0. C++ pin: `git -C /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core show bc9cc12:<path>`

## Missing (all `unmapped` in docs/core-sdk-coverage.json except ImzMLFile = partial)
- FORMAT/XMassFile (+ HANDLERS/AcqusHandler, HANDLERS/FidHandler): Bruker XMass fid/acqus. Pure Rust, no dep. Class test XMassFile_test.cpp + data src/tests/class_tests/openms/data/XMassFile_test/{acqu,acqus,fid,...} + XMassFile_test_data.txt
- FORMAT/HDF5Connector: needs HDF5 crate (hdf5 / hdf5-metno = C lib) -> optional feature. Test HDF5_test.cpp
- FORMAT/BrukerTimsFile (TOPP FileConverter, PeakPickerIM), BrukerTimsImagingFile, RationalScan2ImConverter, HANDLERS/PASEFHillCentroider: opentims (C++). TDF = SQLite (rusqlite already in tree?) + tdf_bin (zstd-compressed frames). timsrust crate exists (pure Rust). Tests BrukerTimsFile_test.cpp, BrukerTimsImagingFile_test.cpp, PASEFHillCentroider_test.cpp
- FORMAT/ThermoRawFile (+ HANDLERS/ThermoRawFileMetadata not in list): openms-thermo-bridge (.NET). External reader process / behind feature. Test ThermoRawFile_test.cpp
- FORMAT/ImzMLFile partial: gap = mzML-level metadata (src/format/mzml.rs strictnesses)

## Findings (learned, keep)
- Incarnation 2: lane sources extracted fresh to scratchpad/F5/src/... (session scratchpad c1467273...; shared with other lanes - only use F5/ subdir).
- C++ lane sources extracted to scratchpad/cpp/ via ONE literal `git -C <core> archive -o <tar> bc9cc12 <paths>` (sandbox refuses git in loops/pipes/heredocs; also no `cd .../source` paths - use absolute paths).
- opentims v1.2.0b4 source (exact tag OpenMS fetches): /Users/kohlbach/Claude/OpenMS-mzpeak-build/build/_deps/opentims-src/src/opentims++/ (opentims.cpp, tof2mz_converter.cpp, scan2inv_ion_mobility_converter.cpp).
- Small TDF fixture (1 MB): "/Users/kohlbach/Claude/mzPeak OpenMS/hupo-mzpeak/test/data/diaPASEF.d" (origin/license to check). Upstream integration data = DDA/DIA_HeLa_50ng_5_6min.d.zip (archive.openms.de, whitelisted HTTP).
- timsrust 0.4.1 unsuitable: rusqlite 0.31 (libsqlite3-sys `links` clash with our 0.40.2), zstd C, own linear calibration. Plan: port opentims frame decode on rusqlite 0.40 + zstd decoder crate.
- thermorawfilereader 0.7/0.8 crate (mobiusklein; .NET host) in local registry = candidate optional dep for ThermoRawFile.
- Cargo: `sqlite = [dep:rusqlite]` (0.40.2 bundled). No zstd/hdf5 dep yet.
- RationalScan2ImConverter::getCalibration fallback = unordered_map::begin() -> hash-order dependent (C++ defect candidate).

## Next step
1. Read the C++ headers/cpps/tests (sizes), decide per format.
2. Port XMassFile + AcqusHandler + FidHandler first (pure Rust, class test data available).
3. RationalScan2ImConverter / PASEFHillCentroider if pure algorithm.
4. BrukerTims via rusqlite+zstd (or timsrust) behind feature `bruker-tdf`.
5. HDF5Connector behind feature `hdf5`; ThermoRawFile behind external process.
