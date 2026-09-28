# Bruker timsTOF (TDF) support

Bruker timsTOF `.d` directories: `analysis.tdf` (SQLite metadata) plus
`analysis.tdf_bin` (zstd-compressed frames). Source headers at core pin
`bc9cc12`:

| Source header | Rust module | Feature | Status |
|---|---|---|---|
| `FORMAT/BrukerTimsFile.h` / `.cpp` | `src/format/bruker_tims_file.rs` | `bruker-tdf` | partial |
| `FORMAT/BrukerTimsImagingFile.h` / `.cpp` | `src/format/bruker_tims_imaging_file.rs` | `bruker-tdf` + `mzml` | partial (inherits the `BrukerTimsFile` gaps) |
| `FORMAT/RationalScan2ImConverter.h` / `.cpp` | `src/format/rational_scan2im_converter.rs` | none; `try_create_rational_converter` needs `sqlite` | complete |
| `FORMAT/HANDLERS/PASEFHillCentroider.h` / `.cpp` | `src/format/pasef_hill_centroider.rs` | none | complete |
| opentims v1.2.0b4 (third-party C++, fetched by the source build) | `src/format/tims_data_handle.rs` | `bruker-tdf` | the parts `BrukerTimsFile` calls |

## Dependencies

OpenMS links the C++ library opentims for the frame codec, `TimsDataHandle`
and the linear converters, and optionally loads Bruker's closed-source
`libtimsdata` at run time. No crate replaces opentims (`timsrust` 0.4.1 pins
`rusqlite` 0.31, whose `libsqlite3-sys` cannot link beside this crate's 0.40,
and brings its own calibration and the C `zstd`), so the opentims layer is
ported onto crates:

- `rusqlite =0.40.2` (already a dependency, feature `sqlite`) for
  `analysis.tdf`;
- `ruzstd =0.8.1` (pure Rust, `std` + `hash`; `hash` pulls `twox-hash
  =2.1.4`) for the frames. 0.8.1 because 0.8.2 calls `is_multiple_of`
  (Rust 1.87) without declaring a `rust-version`; the two decoders are the
  same code. See `docs/THIRD_PARTY_CRATE_DECISIONS.md`.

The feature `bruker-tdf = ["sqlite", "dep:ruzstd"]` is off by default.
The Bruker SDK is not loaded: loading a shared library needs `unsafe`, which
`#![forbid(unsafe_code)]` excludes; see the gaps below.

## The TDF frame format

As opentims `TimsFrame::decompress` / `save_to_buffs` (vendored
`opentims.cpp:89-240`): `Frames.TimsId` is the byte offset of the packet in
`analysis.tdf_bin`; a packet is a little-endian `u32` size (header
included), a `u32` scan count, then one zstd frame. Decompressed, it is
`num_scans + 2 * num_peaks` 32-bit words stored as four byte planes. Word
`s + 1` is twice the peak count of scan `s`; the last scan takes the rest.
Each peak is a TOF delta (running sum from `u32::MAX`) and an intensity.
Intensities are scaled by `100 / Frames.AccumulationTime` and rounded half
up (`intensity * correction + 0.5`, truncated), as the vendored opentims
does.

m/z: `OpenSourceTof2MzConverter`, `mz = (sqrt(mz_min) + (sqrt(mz_max) -
sqrt(mz_min)) / tof_max * tof)^2` from `GlobalMetadata`. 1/K0: the
per-frame rational `TimsCalibration` model (`RationalScan2ImConverter`) when
the table is usable, else opentims `OpenSourceScan2ImConverter`, `1/K0 =
im_max + (im_min - im_max) / scan_max * scan`.

## API mapping

### `BrukerTimsFile` (a `ProgressLogger`)

| Source member | Rust | Notes |
|---|---|---|
| `Config` (all 30 fields, defaults) | `bruker_tims_file::Config`, `Default` | every field, same default; plus native `limits: TdfLimits` |
| `Config::CentroidAlgo`, `ExportMode`, `TimsCalibrationStrategy`, `PressureCompensation` | `CentroidAlgo`, `ExportMode`, `TimsCalibrationStrategy`, `PressureCompensation` | |
| `load(path, exp)` / `load(path, exp, config)` | `BrukerTimsFile::load(path, &Config) -> Result<MSExperiment>` | `Config::default()` for the one-argument form |
| `transform(path, consumer)` / `transform(path, consumer, config)` | `BrukerTimsFile::transform(path, &mut dyn MSDataConsumer, &Config)` | |
| `loadFrames_` | private `load_frames` | `FRAME` export, MS1 `Off` / `HillBased`, isotopic prefilter, hill-bounds arrays |
| `loadDIA_` | private `load_dia` | raw per-frame path only (no aggregation) |
| `loadDDA_` | none | **gap**: `Error::Unsupported` |
| `isDIA_` | `is_dia(&Connection)` | public free function |
| `loadExperimentalSettings_` | private, inside `load` / `transform` | source file, instrument IM unit |
| `resolveEffectiveConfig` (anonymous namespace) | `resolve_effective_config` | frame / RT range resolution |
| `DIAStreamingMetadata`, `readDIAMetadata`, `loadDIAStreaming` | none | **gap**: typed on OpenSWATH `SwathMap` / `FullSwathFileConsumer` |
| `ProgressLogger::setLogType` / `getLogType` | `set_log_type` / `log_type` | |
| — | `DiaWindow`, `read_dia_windows`, `read_frame_to_window_group`, `IM_ARRAY_NAME` | helpers of `loadDIA_` made public for the tests |

Calibration strategy (the source's `openTimsDataHandle`): `Rational` needs a
usable `TimsCalibration` table, else refused; `Linear` always uses the
`GlobalMetadata` range; `Auto` tries the SDK (never available here, as when
the source's SDK load fails), then rational, then linear; `BrukerSdk` is
refused.

### `BrukerTimsImagingFile` (a `ProgressLogger`)

| Source member | Rust | Notes |
|---|---|---|
| `MaldiPixel` | `MaldiPixel` | |
| `Config` (`inner_config`, `strict_imaging_only`) | `bruker_tims_imaging_file::Config` | |
| `load(path, exp)` / `load(path, exp, config)` | `load(path, &Config) -> Result<ImagingExperiment>` | `export_mode` forced to `FRAME`, `load_ms1` to true |
| `isImagingDataset` | `is_imaging_dataset` | `false` on any I/O error |
| `readMaldiFrameInfo` | `read_maldi_frame_info` | |
| `readGlobalMetadataValue` | `read_global_metadata_value` | |
| private `computeBoundingBox_`, `resolveTdfPath_` | private helpers | |

### `RationalScan2ImConverter` and opentims `Scan2InvIonMobilityConverter`

| Source member | Rust | Notes |
|---|---|---|
| base class `Scan2InvIonMobilityConverter` | trait `Scan2InvIonMobilityConverter` | returns a new `Vec` instead of writing a caller's array |
| `Coefficients` (`c0`..`c9`) | `Coefficients` | |
| constructor `(calibrations, frame_to_cal)` | `RationalScan2ImConverter::new` -> `Result` | an empty map is refused (the source dereferences `begin()` of it) |
| `convert` (double and `uint32_t` scans) | `convert`, `convert_indices` | |
| `inverse_convert` | `inverse_convert` | |
| `description` | `description` | |
| private `calibrations_`, `frame_to_cal_` | `calibrations()`, `frame_to_cal()` | `BTreeMap` for `unordered_map` |
| private `getCalibration` | `calibration` | fallback: smallest Id (see divergences) |
| static `applyFormula`, `invertFormula` | `apply_formula`, `invert_formula` | |
| free `tryCreateRationalConverter` | `try_create_rational_converter` (`sqlite`) | |
| opentims `OpenSourceScan2ImConverter` | `LinearScan2ImConverter` | |

### `PASEFHillCentroider`

| Source member | Rust |
|---|---|
| `Params` (six fields, defaults 8.0 / 1.3 / 1 / 0.0 / 0.0 / 0) | `Params`, `Default` |
| `Centroid` (eight fields) | `Centroid` |
| `centroidFrame(mz, intensity, im, n, params, scan_indices = nullptr)` | `centroid_frame(mz, intensity, im, params, Option<&[u32]>)` |

## Native differences and ceilings

- `TdfLimits`: at most 2^25 peaks per frame and 2^24 frames, checked before
  allocation; a packet that runs past `analysis.tdf_bin`, a size field
  smaller than its header, a decompressed length other than the header's
  word count, a zstd content checksum mismatch and scan counts that
  disagree with the payload are all refused (`Error`), never read past.
- `try_create_rational_converter` refuses a `Frames.Id` above 2^24 (the
  source allocates a vector sized by the largest Id).
- Sorting of spectra is stable; the source's `std::sort` leaves the order of
  equal-RT DIA window spectra unspecified.

## Gaps (each refused with `Error::Unsupported`, never silently different)

1. DDA-PASEF export (`loadDDA_`: `SPECTRUM` mode, and `AUTO` on DDA data).
2. Frame aggregation (`ms1_n_neighbors > 0` outside `FRAME`,
   `dia_ms2_n_neighbors > 0`): needs the source's `FrameAggregator`.
3. `Greedy2D` centroiding, which `dia_ms2_centroid` also selects: needs the
   source's `FrameCentroider`.
4. The Bruker SDK calibration (`BrukerSdk`, the SDK step of `Auto`): a
   closed-source shared library. `pressure_compensation` acts only through
   the SDK; otherwise it warns and is ignored, as in the source.
5. `readDIAMetadata` / `loadDIAStreaming`: typed on OpenSWATH classes.
6. TOPP wiring: `FileConverter` and `PeakPickerIM` do not yet read `.d`.

Of the 37 sections of `BrukerTimsFile_test.cpp`, `load` (missing path) and
the five `RationalScan2ImConverter` sections are ported; `[FileHandler]
BRUKER_TDF detection` belongs to `FileHandler` (`tests/file_types.rs`). The
other 30 need `DDA_HeLa_50ng_5_6min.d` / `DIA_HeLa_50ng_5_6min.d`
(archive.openms.de) and are not replicated; most exercise gaps 1-5.

## Evidence

| Test | Tier | Source |
|---|---|---|
| `tests/rational_scan2im_converter.rs` (10) | 3 | `BrukerTimsFile_test.cpp:129-293`: forward 1.070429, scan 0 1.63883, round trip, per-frame, description, singularities; the linear converter on binary-exact values |
| `tests/pasef_hill_centroider.rs` (18) | 1 + 3 | all 15 sections of `PASEFHillCentroider_test.cpp`; 7 synthetic cases bit-identical to the executed Release build (`../oracle/F5-xmass`, ibminode06, gcc 14.4.0, run twice, identical), outputs in `tests/data/pasef_hill_centroider_oracle/` |
| `tests/bruker_tims_file.rs` (14) | 4 | synthetic TDF (`tests/support/tdf_fixture.rs`) encoded independently of the decoder, expected values from the opentims formulas; class-test sections `load` (missing path) and FRAME / range behaviour |
| `tests/bruker_tims_imaging_file.rs` (9) | 3 + 4 | the 6 sections of `BrukerTimsImagingFile_test.cpp` (the class test's `createFakeD`) plus pixel binding on synthetic frames |
| unit tests in the four modules | 4 | formulas evaluated outside Rust (Python IEEE doubles) or exact in binary |

The Release install (`openms4-release-bc9cc12-c19e494-174b576`) is built
without opentims, so `BrukerTimsFile` and `RationalScan2ImConverter` have no
executable C++ oracle there.

**External cross-check (not in the repository).** The 5-frame diaPASEF run
`hupo-mzpeak/test/data/diaPASEF.d` (`analysis.tdf` sha256
`23809f8a6677117b7a745fd3c577106152734336edf2473bde20f10c43f4aa12`,
`analysis.tdf_bin` `9692d7eee2dd2b8cefd07fdb50a945094ba08343fa06122260f2a5d286efbb12`;
no license statement, so not committed) was decoded on spock by the
independent `timsrust` 0.4.1 decoder and by this port (first with `ruzstd`
0.8.2, then again with the pinned 0.8.1, byte-identical output): 291,531
peaks, zero mismatches in
frame, scan, TOF and intensity (`trunc(raw * 100 / AccumulationTime + 0.5)`).
Outputs `ts.tsv` (`1b764e23fc1b1931541b3fb0e4cbc94bdfb271fd47df4d666974214bb4de717d`)
and `om.tsv` / `om081.tsv` (both `f1c69b50eb741f75e26fd6eaf54f81f58ab8632136dd23b1ae48397fe0038a09`)
in `spock:/scratch/kohlbach/f5-tdfcheck/tdfcheck/`.

## C++ defect candidates (unnumbered)

- The default `Config` warns "ms1_centroid_mz_ppm and ms1_centroid_im_pct
  must both be > 0" on every DIA/DDA load (defaults 10 / 0), although the
  default algorithm is `OFF`.
- `sortSpectra` uses `std::sort`: equal-RT DIA window spectra come out in an
  unspecified order.
- The vendored opentims reuses a stale buffer when zstd decompresses to fewer
  bytes than expected, instead of refusing the frame.
- `RationalScan2ImConverter::getCalibration` falls back to
  `unordered_map::begin()`, whose element depends on the hash layout when
  there is more than one calibration segment.
