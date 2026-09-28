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

## Incarnation 3 (started)
- HEAD had WIP drafts of xmass_file.rs, acqus_handler.rs, fid_handler.rs (not yet wired into src/format/mod.rs, not gated).

## Incarnation 4 (started)
- C++ sources already in scratchpad/F5/src/... (list with find, no shell vars in commands: sandbox refuses).
- Drafts xmass/acqus/fid in HEAD~1 still unwired; incarnation 3 left no further notes -> it stalled while exploring.

## Incarnation 4 progress
- XMassFile/AcqusHandler/FidHandler: tests/xmass_file.rs 8/8 PASS on spock (--no-default-features). Remaining: docs/XMASS_FILE_SUPPORT.md, provenance json, unit tests gate.
- PASEFHillCentroider: src/format/pasef_hill_centroider.rs + tests/pasef_hill_centroider.rs (all 15 class-test sections). Gate running.
- RationalScan2ImConverter: pure formula + tryCreateRationalConverter (SQLite -> `sqlite` feature). Class test sections live in BrukerTimsFile_test.cpp lines 129-293 (forward 1.070429, scan0 1.63883, roundtrip, per-frame, description, singularity).
- Gate cmd: separate `--test X` and `--lib -- filter` runs (a filter applies to every binary).

## Incarnation 5 progress
- Gate (spock, --no-default-features --features sqlite): rational_scan2im_converter 10/10, xmass_file 8/8, pasef_hill_centroider 16/16. Logs: session scratchpad F5/gate2.log.
- C++ sources re-extracted to session scratchpad F5/src/openms/... (tar via one git archive).
- zstd decision: ruzstd 0.8.2 (pure Rust, no declared MSRV; 0.8.3+ need 1.87) default-features=false features=[std]. zstd 0.13.3 (C) in local registry as alternative. Cargo.lock must be regenerated on spock and copied back via ssh -J sshgw.
- opentims frame format: opentims.cpp:89-240 (tdf_bin at TimsId: u32 packet size, u32 num_scans, zstd from +8).
- Shell sandbox: no python heredocs mentioning git, no compound cd && git; use Edit tool + plain git commands.

- Oracle ../oracle/F5-xmass (ibminode06 Release, gcc 14.4.0, run twice identical): XMassFile peaks + PASEF 7 synthetic cases. Rust tests bit-exact: xmass 9/9, pasef 18/18 (gate3.log). Release install has NO opentims (no RationalScan2Im / BrukerTimsFile symbols).
- TDF fixture diaPASEF.d (hupo-mzpeak, no license: do NOT commit) = external check only. Committed test must synthesize a TDF (ruzstd encoder compress_to_vec).
- ruzstd 0.8.2 manifest: resolver 3 (ok >= 1.84), edition 2018, std feature, no deps w/o hash.

- DONE: src/format/tims_data_handle.rs + bruker_tims_file.rs (feature bruker-tdf = sqlite + ruzstd 0.8.2 [std,hash] -> twox-hash 2.1.4). Cargo.lock regenerated on spock (fetch) and copied back. tests/bruker_tims_file.rs 14/14 (gate6.log).
- External cross-check (not in repo): diaPASEF.d (sha256 tdf 23809f8a..., tdf_bin 9692d7ee...) decoded by timsrust 0.4.1 vs port: 291531 peaks, 0 mismatches (frame, scan, tof, intensity=trunc(raw*100/AccTime+0.5)); /scratch/kohlbach/f5-tdfcheck on spock; ts.tsv sha 1b764e23..., om.tsv f1c69b50...
- C++ defect candidates: default Config warns "ms1_centroid_mz_ppm and ms1_centroid_im_pct must both be > 0" on every DIA/DDA load (default 10/0); sortSpectra uses std::sort -> DIA equal-RT window spectra order unspecified; opentims short decompression reuses stale buffer; RationalScan2Im getCalibration fallback hash order.

## Next step
1. (done) RationalScan2ImConverter.
2. Docs: docs/XMASS_FILE_SUPPORT.md, docs/BRUKER_TIMS_SUPPORT.md; provenance jsons.
3. TDF reader core (rusqlite + zstd decoder crate) behind feature `bruker-tdf`: decide ruzstd vs zstd.
4. HDF5Connector / ThermoRawFile / BrukerTimsImagingFile / ImzML: decide + record.
5. Lane report + ledger fragment.

## Incarnation 5 (started)
- Branch now worktree-wf_cd066dab-0a7-10 (fast-forwarded onto wf_ee239125-5be-15 @ 06e47f9).
- HEAD has src/format/rational_scan2im_converter.rs (570 lines, lead-saved WIP, ungated).

## Incarnation 7 (started)
- Incarnation 6 left uncommitted: src/format/bruker_tims_imaging_file.rs draft (unwired, ungated), tests/support/tdf_fixture.rs (shared fixture), THIRD_PARTY rows for opentims/HDF5/Thermo. Committed as 93c4d197.
- Next: finish+wire+gate BrukerTimsImagingFile; docs/BRUKER_TIMS_SUPPORT.md; decide HDF5/Thermo/ImzML (record partial); report + ledger.
- Incarnation 7: BrukerTimsImagingFile wired (cfg bruker-tdf+mzml), tests/bruker_tims_imaging_file.rs 9/9 + bruker_tims_file 14/14 on spock (scratch F5/gate7.log). Parsers made StringUtils-exact.
- HDF5: no hdf5 crate in local or spock registry, no libhdf5 on spock -> backlog (named dep hdf5-metno, feature hdf5). Thermo: thermorawfilereader 0.8 default feature downloads nethost at build; deps (dotnetrawfilereader-sys, netcorehost, flatbuffers) not cached; RawFileReader license -> backlog.
- Next: docs/BRUKER_TIMS_SUPPORT.md, provenance, report+ledger. ImzML mzml.rs gap optional.
