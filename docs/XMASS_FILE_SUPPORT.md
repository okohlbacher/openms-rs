# XMassFile, AcqusHandler and FidHandler support

Bruker "XMass Analysis (fid)" spectra, the native format of the Bruker Flex
series (MALDI-TOF). Source headers at core pin `bc9cc12`:

| Source header | Rust module |
|---|---|
| `FORMAT/XMassFile.h` (header-only; `XMassFile.cpp` holds the constructor and destructor) | `src/format/xmass_file.rs` |
| `FORMAT/HANDLERS/AcqusHandler.h` / `.cpp` | `src/format/acqus_handler.rs` |
| `FORMAT/HANDLERS/FidHandler.h` / `.cpp` | `src/format/fid_handler.rs` |

No third-party library is involved: the source reads both files with the
standard library only. No feature gate.

## API mapping

### `XMassFile` (a `ProgressLogger`)

| Source member | Rust | Notes |
|---|---|---|
| `XMassFile()` | `XMassFile::new`, `Default` | |
| `~XMassFile()` | drop | |
| `load(filename, MSSpectrum&)` | `XMassFile::load` -> `Result<MSSpectrum>` | returns a new spectrum; default `LoadOptions` |
| (same) | `XMassFile::load_with_options` | native options: truncated-`fid` handling and a point ceiling |
| `importExperimentalSettings(filename, PeakMap&)` | `XMassFile::import_experimental_settings(path, &mut MSExperiment)` | `exp` unchanged on error |
| `store(filename, const MSSpectrum&)` | `XMassFile::store` | always `Error::Unsupported`, the source's `NotImplemented` |
| `ProgressLogger::setLogType` / `getLogType` | `set_log_type` / `log_type` | no member reports progress, as in the source |
| — | `LoadOptions`, `LoadOptions::source`, `TruncatedFid`, `DEFAULT_MAX_POINTS` | native |

### `Internal::AcqusHandler`

| Source member | Rust | Notes |
|---|---|---|
| `AcqusHandler(const std::string& filename)` | `AcqusHandler::open(path)` | |
| (stream constructor, none in source) | `AcqusHandler::from_reader` | native, for in-memory text |
| `~AcqusHandler()` | drop | |
| `getPosition(Size index) const` | `position(u64) -> f64` | same operations, same order, `f64` |
| `getParam(const std::string&)` | `param(&str) -> String` | computed `mzMin` / `mzMax` as in the source |
| `getSize() const` | `size() -> u64` | `$TD` |
| private `params_`, `dw_`, `delay_`, `ml1_`, `ml2_`, `ml3_`, `td_` | private fields; `params()` exposes the map read-only | |
| — | `MAX_ACQUS_BYTES` | native input ceiling |

### `Internal::FidHandler` (an `std::ifstream`)

| Source member | Rust | Notes |
|---|---|---|
| `FidHandler(const std::string& filename)` | `FidHandler::open(path)` | |
| (any stream) | `FidHandler::new(reader)` | native, generic over `Read` |
| `~FidHandler()` | drop | |
| `getIndex() const` | `index() -> u64` | |
| `getIntensity()` | `intensity() -> Result<u32>` | refuses a short read |
| (same) | `intensity_source() -> Result<u32>` | the source's short-read value |
| inherited `std::ifstream` interface | not exposed | only the two declared members are used by `XMassFile` |

## Preserved source conventions

- The `acqus` path is the `fid` path without its last three bytes plus
  `acqus`, whatever those bytes are.
- `acqus` line grammar exactly as `AcqusHandler.cpp:45-61` (`##` prefix,
  exactly one `=`, value trimmed, last repeat wins).
- The m/z calibration of Titulaer *et al.* (BMC Bioinformatics 7:403, 2006),
  operation for operation; NaN on a negative discriminant, unchecked.
- `$TD - 1` for `mzMax` wraps in unsigned arithmetic when `$TD` is 0.
- Intensities are the `fid` 32-bit integers converted to `f32`.
- Every metadata assignment of `load` and `importExperimentalSettings`,
  including the `LD+` / `LD-` polarity rule, the direct-inlet rule that leaves
  the ionization method unknown, and the `<...>` stripping of `$ID_raw`,
  `$ACQMETH`, `$TgIDS` and `$AQ_DATE`.

## Native differences

- A `fid` shorter than `$TD` points is refused by default
  (`TruncatedFid::Refuse`). The source reads past the end and pads with a
  partial value and zeros; `LoadOptions::source()` reproduces that.
- `$TD` above `LoadOptions::max_points` (default 2^24) is refused before
  allocation; a negative `$TD`, which the source's cast wraps to about 2^64,
  always is.
- `acqus` larger than `MAX_ACQUS_BYTES` (16 MiB) is refused before reading.
- A non-UTF-8 `acqus` line is decoded as ISO-8859-1 (reversible); the source
  has no notion of encoding.
- `load` returns a new spectrum instead of filling a caller's; an error
  leaves nothing half-filled.
- The spectrum comment `no comment` is not reproduced: the kernel
  `MSSpectrum` has no comment field (`SpectrumSettings::comment`, a kernel
  gap outside this group). The Release oracle prints it; the Rust test does
  not assert it.
- `store` documents `FileNotWritable` in the source but throws
  `NotImplemented`; the code is followed.

## Checked boundaries and evidence

`tests/xmass_file.rs` (9 tests, no feature gate):

| Test | Tier | Source |
|---|---|---|
| `load_matches_the_release_build_bit_for_bit` | 1 | every peak's m/z and intensity bits (FNV-1a digest plus ten sampled bit patterns) and the source-file size against the executed Linux x86_64 Release build (`../oracle/F5-xmass`, ibminode06, run twice, identical) |
| `load_matches_the_class_test_expectation_file` | 3 | `XMassFile_test.cpp` and `XMassFile_test_data.txt`, 80478 pairs, `TOLERANCE_ABSOLUTE(0.001)` |
| `load_sets_the_source_metadata`, `import_experimental_settings_reads_acqus` | 3 | source assignments in `XMassFile.h`; the Release oracle printed the same values |
| `missing_file_is_not_found`, `store_is_not_implemented`, `handlers_read_the_class_test_files` | 3 | class test sections |
| `truncated_fid_is_refused_by_default_and_padded_on_request`, `oversized_td_is_refused_before_allocation` | 4 | native bounds |

The oracle driver, case generator and outputs live outside the repository in
`../oracle/F5-xmass/` (`driver.cpp`, `remote_run.sh`, `results/`); their
sha256 are recorded in `tests/data/xmass_file_provenance.json`.
