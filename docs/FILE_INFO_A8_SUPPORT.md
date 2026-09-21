# FileInfo: the mzXML, mzData, MGF and MS2 peak files (A8, the cheap half)

Package **A8-FILEINFO**, its first half. A4 ran the peak-file branch of
`OpenMS::FileInfo` for DTA, DTA2D and mzML and refused seven more peak types
before loading them. Four of the seven already had readers in this crate, so
closing them is wiring, not porting:

| Type | Source load step (`FORMAT/FileHandler.cpp`) | Rust reader |
| --- | --- | --- |
| mzXML | `:886-893`, `MzXMLFile` with `f.getOptions() = options_` | `src/format/mzxml.rs`, `ReadOptions::default()` |
| mzData | `:895-902`, `MzDataFile` with `f.getOptions() = options_` | `src/format/mzdata.rs`, default `PeakFileOptions` and `ReadLimits` |
| MGF | `:923-929`, `MascotGenericFile::load` | `src/format/mascot_generic.rs`, `CarryOver::Source` and `source_ms_level` |
| MS2 | `:931-937`, `MS2File::load` | `src/format/ms2.rs` |

The report written from the loaded experiment is A4's peak-file branch
(`FORMAT/FileInfo.cpp:1532-1965`, with the `-m`, `-p` and `-s` arms at
`:2005-2081`, `:2115-2127` and `:2384-2440`) and A6's `-i`, `-d` and `-c`; none
of it changes. What changes is `src/format/file_info/report.rs`, which routes
the four types to that branch, and `src/format/file_info/peaks.rs`, which loads
them.

**Still refused**, with `FileInfo peak-file branch for <type> input is not
ported`: sqMass, XMass (`fid`) and MSP, which have no reader here that fills an
`MSExperiment`, and Thermo RAW and Bruker TDF, as before. The rest of A8 — `-v`
and the pepXML, mzTab, trafoXML and PQP branches — is the other half.

Evidence and hashes: `tests/data/file_info_a8_provenance.json`.
Tests: `tests/file_info_a8.rs` (19), and three older files that change with the
scope: `tests/topp_file_info.rs` reproduces TOPP_FileInfo_4, _5 and _6 through
FuzzyDiff against the retained upstream outputs instead of listing them as not
ported; `tests/file_info.rs` drops the four types from its refusal table and
keeps a tripwire that fails if one of them goes back to refusing.
Oracle: `../oracle/a8-fileinfo`, the A7 machinery, **70** cases of the Release
C++ FileInfo of `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`
on ibminode06, plus **14** cases of `driver/fileinfo_driver.cpp`, which calls
`OpenMS::FileInfo::run` against the same install; all 84 run twice and
reproduced. **61** reports are compared byte for byte, text and TSV, with only
the input path normalised.

---

## 1. API mapping

Every source member A8 reaches, and its counterpart.

| Source | Rust | Note |
| --- | --- | --- |
| `FileInfo::run`, the `else // peaks` arm for `MZXML`, `MZDATA`, `MGF`, `MS2` (`FORMAT/FileInfo.cpp:1532-1534`) | `report::FileInfo::run` → `peaks::report` | the four types now map to `Branch::Peaks` |
| `FileHandler::loadExperiment(in, exp, {in_type}, …)` for the four types (`FORMAT/FileHandler.cpp:849-937`) | `peaks::load_source_reader` (private) | not `FileHandler::load_experiment_with_options`: that loader has no mzXML or mzData case and reads MGF with the strict native `mgf` adapter, not `MascotGenericFile` |
| `getType(filename)` and the allowed-type check (`:856`, `:858-864`) | `FileHandler::get_type`, then `Error::InvalidValue` for a type other than the forced or detected one | native difference 1 |
| `MzXMLFile::load` with default `PeakFileOptions` | `mzxml::load_with_options(path, &ReadOptions::default())` | |
| `MzDataFile::load` with default `PeakFileOptions` | `mzdata::load_with_options(path, &PeakFileOptions::default(), &ReadLimits::default())` | |
| `MascotGenericFile::load` (`FORMAT/MascotGenericFile.h:74-104`) | `mascot_generic::load_with_options` with `carry_over: CarryOver::Source` and `source_ms_level: true` | section 2 |
| `spectrum.setMSLevel(std::stoi(tmp))` (`FORMAT/MascotGenericFile.h:325-331`) | **new** `mascot_generic::ReadOptions::source_ms_level` | the one reader change of this package, section 2.2 |
| `MS2File::load` (`FORMAT/MS2File.h`, header-inline) | `ms2::load` | |
| `static_cast<Int>(l)` for the result's per-level keys (`FORMAT/FileInfo.cpp:1645-1657`) | `peaks::to_int`, now modular | section 2.3 |

The FileInfo tool (`OpenMS4-topp/src/FileInfo.cpp`, topp `174b576`) is not
changed; its `-in` and `-in_type` lists (`OpenMS4-topp/src/FileInfo.cpp:83-87`)
already contained `mzData`, `mzXML` and `mgf` and never `ms2`, and its refusal
of `-i` on anything but mzML (`OpenMS4-topp/src/FileInfo.cpp:118-121`) was
already ported. Both are compared with the Release build below.

## 2. Preserved source conventions

### 2.1 The source's readers, with the source's options

FileInfo hands `loadExperiment` a `FileHandler` whose `PeakFileOptions` are the
defaults, and the handler copies them into `MzXMLFile` and `MzDataFile` only.
`MascotGenericFile` and `MS2File` receive none (`FORMAT/FileHandler.cpp:923-937`),
so no filter reaches them here either.

The MGF reader is the one with a choice to make. `MascotGenericFile::load`
declares one spectrum outside its loop, and `getNextSpectrum_` clears only its
peaks, native ID, `TITLE` and `SEQ` (`FORMAT/MascotGenericFile.h:89-104`,
`:141-157`), so a block that omits `CHARGE=`, `PEPMASS=`, `RTINSECONDS=` or
`MSLEVEL=` keeps the previous block's value — and the report counts it.
`mascot_generic`'s native default starts every block fresh (`CarryOver::Reset`,
`docs/MASCOT_GENERIC_SUPPORT.md`); FileInfo passes `CarryOver::Source`.
`a8_mgf_carry.mgf` is the measurement: four blocks whose charges the Release
build counts as 2, 2, 3, 3 and whose MS levels as 2, 2, 1, 1, where the native
default would count charges 2, 0, 3, 0.

### 2.2 `MSLEVEL=-1` is the MS level 4294967295

`std::stoi` returns `-1` and `setMSLevel` stores it in a `UInt`, so the Release
build prints `MS levels: 4294967295`, a `MS Level 4294967295 Ranges:` block and
`number of MS4294967295 spectra` in the TSV (`g_mslevel_negative_all`). The MGF
reader refuses a non-positive level by default, because such a record is invalid
for every consumer. The new `ReadOptions::source_ms_level` stores the value as
the source does — `int` to `UInt` is modular since C++20 — and FileInfo sets it.
The change to the reader is one field of `ReadOptions` and one guarded match
arm in the `MSLEVEL` branch of `read_header_line`, both in
`src/format/mascot_generic.rs`; the default is unchanged.

### 2.3 The result's `Int` keys wrap as `static_cast<Int>` does

The structured `PeakInfo` keys its per-level maps by `Int`, which the source
fills with `static_cast<Int>(level)` (`FORMAT/FileInfo.cpp:1645-1657`). A4 refused
a level above `i32::MAX` there; the conversion is defined as modular since
C++20, so `peaks::to_int` now keeps the bits and MS level 4294967295 is the key
`-1`, as in the source (`mgf_negative_ms_level_wraps_as_in_the_source`).

### 2.4 Readers' own source behaviour, visible through the report

Each of these is the reader's, measured here through the report the Release
build writes, and needed no change:

- **mzXML.** `msLevel="0"` is read as MS1 with a warning
  (`FORMAT/HANDLERS/MzXMLHandler.cpp:247-252`); an `activationMethod` short name
  is looked up and an unknown one dropped; every precursor of a scan counts its
  activation methods, the first alone its charge; a scan's data processing is
  the document's, so `-p` reports the file-level `<dataProcessing>` and a
  `centroided="1"` there makes an unknown spectrum type read as centroid.
- **mzData.** A second `ChargeState` resets the charge to 0
  (`FORMAT/HANDLERS/MzDataHandler.cpp:1194-1204`); an activation method outside
  the handler's fixed vocabulary is read as its first entry, CID; an invalid
  `spectrumType` leaves the type unknown, which `PSI:1000127` peak picking in
  the data processing then turns into centroid; `TimeInMinutes` is converted;
  `.mzDat` is no mzData extension and is found by content.
- **MGF.** A `TITLE` holding `min` sets the retention time from its comma
  field, and one that does not convert keeps the previous time; a block without
  a peak line runs on into the next and the two become one spectrum; tabs,
  runs of spaces and a third, charge column on a peak line are accepted; MGF is
  centroided by definition (`FORMAT/MascotGenericFile.h:95`).
- **MS2.** A peak line before the first `S` line is dropped, an `S` line with no
  peaks is an empty spectrum, `Z` and `I` lines are skipped, so the charge is
  always 0 and the retention time the default `-1`.

### 2.5 What the tool refuses before the class runs

The tool's `-in` accepts a fixed list of formats without `ms2`, so every MS2
case the oracle ran through the tool exits 6 with `has invalid format 'ms2'`,
in the Release build and here alike
(`an_ms2_input_is_refused_by_the_input_format_check`). The library branch is
reached only through the class, which is why the oracle has a driver: its MS2
reports are compared in `ms2_reports_match_the_release_build_class`. Likewise
`-i` on mzXML: the tool refuses it (`x4_i`), the class runs the indexed-mzML
check on the file and reports that it has no index (`lib_x4_i`), and the port's
class does the same.

## 3. Native differences

1. **A forced type the loader's detection contradicts is `Error::InvalidValue`,
   exit 6.** The source throws `ParseError` (`FORMAT/FileHandler.cpp:858-864`)
   and the tool exits 3 (`g_forced_on_mzxml`, `lib_m_forced_mgf`). The port
   answers as it already does for every other peak type, which
   `docs/TOPP_FILE_INFO_SUPPORT.md` records as native difference 3; giving the
   four new types the source's exit code alone would make the code depend on
   which type was forced. Neither build writes a report.
2. **`CHARGE=2-` is a parse error, exit 3.** The source's `StringUtils::toInt32`
   throws `ConversionError`, which leaves `load` unconverted and which the tool
   reports as an unexpected internal error, exit 8 (`g_charge_minus`). Both
   refuse the file and write no report.
3. **Error messages are this crate's.** Where both builds refuse a malformed MGF
   or MS2 file with a parse error (`g_truncated`, `g_one_value`, `lib_m_bad_s`,
   `lib_m_bad_peak`), the exit code agrees and the message text does not; the
   reports agree in being absent.

## 4. What this port refuses where the Release build reports

Two inputs load and are then refused by the kernel's range computation, which
validates each spectrum in full (`MSSpectrum::range_manager`,
`src/kernel/ranges.rs:1270-1272`, calling `MSSpectrum::validate`) where the
source's `MSSpectrum::updateRanges` reads only the peaks and the ion mobility
(`KERNEL/MSSpectrum.cpp:573-621`). Neither is an out-of-bounds read or a crash;
both reports the Release build writes are deterministic and are kept as
`tests/data/file_info_a8/expected/d_mzdata1*.{txt,tsv}` and
`g_mslevel_zero_all.{txt,tsv}`, so the day the kernel relaxes, the tests that
now assert the refusal turn into byte comparisons with no new oracle run.

1. **An mzData scan window that begins after it ends.** `MzDataFile_1.mzData`,
   an upstream class-test fixture, has a spectrum with `mzRangeStart="110"` and
   no `mzRangeStop`. `MzDataHandler` keeps the window `[110, 0]` because one
   bound is non-zero (`FORMAT/HANDLERS/MzDataHandler.cpp:392-395`), and so does
   the port's reader, deliberately. `ScanWindow::validate` refuses it
   (`src/metadata/acquisition.rs:351`): `Error::InvalidValue("scan window begin
   exceeds end")`, three oracle cases (`d_mzdata1`, `_all`, `_dc`).
2. **An MGF spectrum of MS level 0.** `MSLEVEL=0` loads as level 0, as in the
   source, and `MSSpectrum::validate` refuses an MS level of 0 outside the three
   optical scan modes (`src/kernel.rs:809-824`): `Error::InvalidValue("spectrum
   MS level must be positive")`, one case (`g_mslevel_zero_all`). A6 recorded
   the same kernel boundary for mzML.

Both belong to the kernel, not to this branch or the readers: the minimal
change is for `MSSpectrum::range_manager` to validate what it reads — the RT,
the peaks and the mobility values — rather than the whole spectrum. That file
is outside this package, so the change is left to its owner.

3. **A gzip-compressed MGF.** `FileHandler::getType` strips `.gz` and answers
   MGF in both builds; `MascotGenericFile::load` then reads with a plain
   `std::ifstream` (`FORMAT/MascotGenericFile.h:81-84`), so the source sees the
   compressed bytes as text, finds no `BEGIN IONS` line and reports an empty
   map with exit 0 (`g_gzipped_all`, kept as `expected/g_gzipped_all.*`). The
   MGF reader does not decompress either, and its line reader refuses bytes
   that are not UTF-8 text: `Error::Io` of kind `InvalidData`, which the tool
   reports as an unexpected error, exit 8. This one is a choice, not a kernel
   boundary: reproducing the source would take a byte-line reader in the
   shared text input of `src/format/ms2.rs`, to report an empty map for a file
   that was never read. The XML readers decompress in both builds, and the
   gzip-compressed `a8_mzxml_edges` is compared byte for byte (`x_edges_gz`).

## 5. Checked boundaries and evidence

`../oracle/a8-fileinfo/manifest.json` records, per case, argv, exit code, stdout,
stderr and every file the case left, with sha256, for two runs; both runs agree
on every case. The fixtures are 12 test-data inputs at `0cb15f2` (`topp/`,
`spectra_spectrast.mzXML` from `topp/THIRDPARTY/`), 14 core class-test inputs
at `bc9cc12` (`src/tests/class_tests/openms/data/`) and 21 generated by
`scripts/make_a8_fixtures.py`, which rewrites them byte for byte.

| Format | Tool cases | Driver cases | Compared byte for byte | Refused here, reported there |
| --- | --- | --- | --- | --- |
| mzXML | 21 | 2 | 21 (+ `x4_bare` on stdout) | — |
| mzData | 18 | 1 | 16 | 3 (section 4.1) |
| MGF | 22 | 1 | 17 | 2 (sections 4.2, 4.3) |
| MS2 | 9, all refused by `-in` in both | 10 | 7 | — |

Every kind of refusal the Release build answers with is asserted here too,
with the difference, if any, in section 3: the malformed MGF and MS2 files and
the forced types in `tests/file_info_a8.rs`, the `ms2` extension (`m_test`, and
the eight other `m_*` tool cases that take the same path) and `-i` on mzXML
(`x4_i`) in `tests/topp_file_info.rs`. Upstream registrations: TOPP_FileInfo_4 (`topp/CMakeLists.txt:890-892`,
`-m`), _5 (`:893-895`, `-in_type mzData -m -s` on `FileInfo_5_input.mzDat`) and
_6 (`:896-898`, `-d -s`) pass FuzzyDiff with the whitelist `File name` against
their retained outputs, which the Release build's `-out` also equals apart from
the file name line; none of the three passed before, all three refused.
