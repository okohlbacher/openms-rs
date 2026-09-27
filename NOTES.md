# Lane F5 NOTES (restart-safety; delete in final commit)

Branch: worktree-wf_ee239125-5be-15, base f0990c0. C++ pin: `git -C /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core show bc9cc12:<path>`

## Missing (all `unmapped` in docs/core-sdk-coverage.json except ImzMLFile = partial)
- FORMAT/XMassFile (+ HANDLERS/AcqusHandler, HANDLERS/FidHandler): Bruker XMass fid/acqus. Pure Rust, no dep. Class test XMassFile_test.cpp + data src/tests/class_tests/openms/data/XMassFile_test/{acqu,acqus,fid,...} + XMassFile_test_data.txt
- FORMAT/HDF5Connector: needs HDF5 crate (hdf5 / hdf5-metno = C lib) -> optional feature. Test HDF5_test.cpp
- FORMAT/BrukerTimsFile (TOPP FileConverter, PeakPickerIM), BrukerTimsImagingFile, RationalScan2ImConverter, HANDLERS/PASEFHillCentroider: opentims (C++). TDF = SQLite (rusqlite already in tree?) + tdf_bin (zstd-compressed frames). timsrust crate exists (pure Rust). Tests BrukerTimsFile_test.cpp, BrukerTimsImagingFile_test.cpp, PASEFHillCentroider_test.cpp
- FORMAT/ThermoRawFile (+ HANDLERS/ThermoRawFileMetadata not in list): openms-thermo-bridge (.NET). External reader process / behind feature. Test ThermoRawFile_test.cpp
- FORMAT/ImzMLFile partial: gap = mzML-level metadata (src/format/mzml.rs strictnesses)

## Next step
1. Read the C++ headers/cpps/tests (sizes), decide per format.
2. Port XMassFile + AcqusHandler + FidHandler first (pure Rust, class test data available).
3. RationalScan2ImConverter / PASEFHillCentroider if pure algorithm.
4. BrukerTims via rusqlite+zstd (or timsrust) behind feature `bruker-tdf`.
5. HDF5Connector behind feature `hdf5`; ThermoRawFile behind external process.
