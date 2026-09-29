# DTA, MS2 and DTA2D peak lists

`format::ms2` and `format::dta2d` provide buffered `read`, `read_with_options`,
atomic `read_into`/`read_into_with_options`, `load`/`load_with_options`, and
`write`/`write_with_options`/`store`/`store_with_options`. The path wrappers delegate
to the same adapters and, as the source's `load` does, record the document's
identity in `settings.document`: `loaded_file_path` is the path as given (a
relative path prefixed with the current directory) and `loaded_file_type` the
type recognised from the file's content (`DTA2DFile.h:97-98`, `MS2File.h:66-67`;
`DTA2DFile_test.cpp:72-73` and `MS2File_test.cpp:49-50` are
`path_loads_record_the_loaded_file_path_and_type`). A stream read records
neither, since a stream has no name.

The pinned [MS2File header](https://github.com/okohlbacher/OpenMS4-core/blob/6bfc0e4711105f4eda2fea86812a83af7c7e791f/src/openms/include/OpenMS/FORMAT/MS2File.h)
is an **input-only** adapter. Its Rust writer is an explicit native extension.
[DTA2DFile](https://github.com/okohlbacher/OpenMS4-core/blob/6bfc0e4711105f4eda2fea86812a83af7c7e791f/src/openms/include/OpenMS/FORMAT/DTA2DFile.h)
has source load, store and TIC-store operations, all represented here. No native
dependency or external converter is involved.

## API mapping

| Source member | Rust | Notes |
|---|---|---|
| `MS2File()`, `~MS2File()` | none | the class has no state; its operations are free functions in `format::ms2` |
| `template <typename MapType> void MS2File::load(filename, exp)` | `ms2::load`, `ms2::load_with_options` | `MapType` is `MSExperiment`, see below; stream forms `read`, `read_with_options`, `read_into`, `read_into_with_options` |
| `MS2File`'s `ProgressLogger` base | none owed | its only `startProgress` is commented out (`MS2File.h:52`); the Release build makes no call (tier 1, `tests/progress_format_readers.rs`) |
| (none) | `ms2::write`, `write_with_options`, `store`, `store_with_options` | native extension; the source has no MS2 writer |
| `DTA2DFile()`, `~DTA2DFile()` | none | the only state is `options_`, passed as an argument |
| `PeakFileOptions& getOptions()`, `const PeakFileOptions& getOptions() const` | `dta2d::ReadOptions` (read), `dta2d::WriteOptions` (write) | the adapter consumes only `rt_range`, `mz_range`, `intensity_range` of `PeakFileOptions`; `ReadOptions::from_peak_options` takes a whole source `PeakFileOptions` and keeps those three, the rest having no DTA2D behaviour; the source's store reads no option |
| `template <typename MapType> void DTA2DFile::load(filename, map)` | `dta2d::load`, `load_with_options`, `load_with_progress` | stream forms `read`, `read_with_options`, `read_into`, `read_into_with_options` |
| `template <typename MapType> void DTA2DFile::store(filename, map) const` | `dta2d::store`, `store_with_options`, `store_with_progress` | stream forms `write`, `write_with_options`; `WriteOptions::source()` is the source's behaviour |
| `template <typename MapType> void DTA2DFile::storeTIC(filename, map) const` | `dta2d::store_tic`, `store_tic_with_options`, `store_tic_with_progress` | stream forms `write_tic`, `write_tic_with_options` |
| `DTA2DFile`'s `ProgressLogger` base | the three `*_with_progress` functions | see Progress |

**`MapType`.** All four templates are documented as taking "a MSExperiment or
[a type with] the same interface". Every instantiation in the pinned source is
`PeakMap`, i.e. `MSExperiment` (`FileHandler.cpp:877-883`, `:931-937` and
`:1159-1165`, and both class tests), and this crate has one experiment type,
so the port fixes `MapType` to `MSExperiment` rather than inventing a trait
with one implementor.

## DTA

`format::dta` ports `FORMAT/DTAFile.h`, the single-spectrum DTA file: its
`load` template is `dta::load` (path) and `dta::read` (stream), its `store`
template `dta::store`/`store_with_options` (path) and `write`,
`write_with_options`, `write_with_convention` (stream). The constructor and
destructor have no counterpart, since the class's only member,
`default_ms_level_`, is set to 2 by the constructor (`DTAFile.cpp:16-19`) and
changed by nothing; it is the constant `dta::DEFAULT_MS_LEVEL`.
`SpectrumType` is `MSSpectrum`, the only instantiation in the pinned source
(`FileHandler.cpp:796`, `:832`, `:873`, `:1155`).

- **Reading.** The first line is the MH+ mass and the charge; a charge of 0
  keeps the mass as the precursor m/z, any other charge `z` gives
  `(MH+ - PROTON_MASS_U) / z + PROTON_MASS_U` (`DTAFile.h:118-126`). Every
  line is trimmed of spaces, tabs, CR and LF, split on a tab if it has one and
  on single spaces otherwise, and must give exactly two fields, whose
  surrounding whitespace the number conversion skips (`StringUtils::toDouble`,
  `toInt32`). Two adjacent spaces therefore make a line of three fields and a
  parse error, as in the source. Blank lines after the first are skipped.
  `dta::load` names the spectrum after the file's base name
  (`DTAFile.h:164`); a stream read leaves the name empty. Nonfinite values are
  a checked parse error, where the source would store them.
- **Writing.** The numeric text, the two mass conventions and the lossless
  default with `WriteOptions::source()` for the source's silent discard are
  described in the `format::dta` module documentation. `store_with_options`
  produces and checks the text before it creates the file, so a refused
  spectrum leaves no file; the source creates the file first. The source's
  `std::cerr` warning when it drops a second precursor has no counterpart; the
  default refuses that spectrum instead.
- **Evidence.** `DTAFile_test.cpp`'s load and store sections are
  `dta_class_test_load` and `dta_class_test_store` (tier 3, with the class
  tests' `TEST_REAL_SIMILAR` tolerances); the line splitting is
  `dta_lines_split_as_the_source_splits_them` (tier 4 from `DTAFile.h:87-160`
  and `StringUtils.cpp`). The writer's text is tier 1 through DTAExtractor
  (`tests/topp_dta_extractor.rs`, `docs/TOPP_CLI_SUPPORT.md`).

## MS2 conventions

The source reader trims ASCII space/tab/CR/LF and collapses that whitespace
between fields. A line starting `S` must contain four fields; the fourth is the
precursor m/z. The scan-number fields are not interpreted. Each scan gets MS
level 2, one precursor with unknown charge/zero intensity, RT -1 and an
`index=N` native ID. Empty scans survive, and peak order is unchanged.

Lines starting `H`, `I`, `Z` or `D` are ignored regardless of their remaining
text, exactly as in the pinned adapter. In particular, **Z charge and I RTime
are not loaded**. Valid peak lines before the first S record are parsed and
counted against limits but discarded; a file with no S records has no spectra.
Other nonempty lines require exactly two finite numeric fields.

The native writer emits S records and peaks only. It requires MS2, unset RT and
exactly one precursor without charge, intensity or acquisition metadata. It
rejects populated names, spectrum type, annotations, identifications and
experiment metadata/chromatograms before producing output. Empty and ordinary
`index=digits` IDs are accepted as source-generated bookkeeping and regenerated;
other IDs are rejected. These restrictions make the stored scientific fields
readable by the source-compatible reader without silent loss.

## DTA2D conventions and filters

Data rows have exactly three fields: RT in seconds, m/z and intensity by default.
If a line contains a tab, the source uses tab delimiters; otherwise it uses single
spaces. Repeated delimiters create empty fields and are not silently collapsed.
Trailing whitespace is trimmed. Every `#` line is a header, not a free comment.
The first three header fields select any order of `SEC`/`MIN`, `MZ`, `INT`,
case-insensitively. Deprecated `RT`, `RETENTION_TIME`, `MASS-TO-CHARGE`, `IT` and
`INTENSITY` remain accepted. Additional header fields are ignored. A short or
invalid header returns a checked parse error.

Minutes multiply RT by 60. The source's minutes flag is sticky: after a MIN
header, a later SEC header does not reset it. A new spectrum begins only when
`abs(row_rt - spectrum_anchor_rt) > 0.0001`; the anchor is the first RT in the
group, not the preceding row. The initial anchor is the source's -1 sentinel.
Rows within its tolerance retain RT -1 and an empty native ID until a real RT
change. Later IDs count all raw groups, including those removed by filtering.
Spectra are MS1, unsorted input order is preserved, and empty filtered groups
are omitted.

`dta2d::ReadOptions` exposes `rt_range`, `mz_range`, `intensity_range` as optional
Rust ranges, plus `limits`. All three ranges include their start and exclude
their end, following OpenMS DRange. Intensity filtering uses the parsed f32
value widened to f64. Empty ranges select nothing; reversed/nonfinite endpoints
are errors. These are the only PeakFileOptions fields consumed by the C++
adapter. Generic MS-level, metadata-only, compression and precision flags have
no DTA2D behavior and are not presented as working filters.

Full storage emits `#SEC\tMZ\tINT` and one row per peak. Every number is spelled
as the source spells it: `DTA2DFile::store` writes
`precisionWrapper(rt)`, `precisionWrapper(mz)` and `precisionWrapper(intensity)`,
and `PrecisionWrapper`'s `operator<<` is `StringUtils::toStr(value, true)`
(`CONCEPT/PrecisionWrapper.h`). For a double that is 15 *fraction* digits in
`[0.01, 10000)` and the shortest round-trip scientific text elsewhere
(`4711.100000000000364`, `230.02000000000001`); for the float intensity 6
fraction digits in that range and the shortest float scientific text elsewhere
(`318.519989`, `4.721889e04`). The port reuses
`format::file_info::text_format::{to_str, to_str_f32}`, which FileInfo already
proved against the Release build. The stored file of `DTA2DFile_test_1.dta2d`,
and its TIC, are byte-identical to the Release build's
(`dta2d_store_writes_the_release_build_bytes`, tier 1).

That text is not always exact: a double below about 8, or a float in
`[0.01, 1)` with more than 6 significant fraction digits, reads back as a
different value. **By default the writers refuse such a number** with
`Error::Unsupported`, before any byte is written, so a stored file always reads
back as the experiment it came from. They also require nonempty MS1 spectra
without precursors and the same metadata restrictions as MS2, and reject
adjacent spectra within the grouping tolerance or a first RT that would be
changed by the -1 sentinel. Empty experiments are valid and produce the header.

`dta2d::WriteOptions::source()` sets `discard_unrepresentable` and writes what
the source writes: every spectrum's peaks under its retention time whatever its
MS level, precursors or metadata, an empty spectrum contributing no row, and
every number in the source text even where it rounds. Nonfinite values and the
output ceilings are refused in both modes; the source would write `NaN` or
`inf` text.

## TIC projection

`write_tic`, `write_tic_with_options`, `store_tic` and `store_tic_with_options`
explicitly project MS1 spectra into RT/zero-m/z/intensity rows. They use f32
intensity accumulation in peak order, include zero-intensity points for empty
MS1 spectra, and preserve scan order and duplicate RTs. They do not resample or
use stored chromatograms. Higher MS levels and unrelated metadata are
unconsumed, as in `calculateTIC()` with source defaults. Duplicate-RT rows can
be grouped if the resulting file is read as DTA2D spectra. The retention time and
the summed intensity are spelled as in full storage, the m/z as the literal `0`,
and the same number-text check applies unless `discard_unrepresentable` is set.

## Progress

`DTA2DFile` derives from `ProgressLogger`. `dta2d::load_with_progress`,
`store_with_progress` and `store_tic_with_progress` make the calls of its
`load`, `store` and `storeTIC` on a caller's logger:
`startProgress(0, 0, "loading DTA2D file")` before the file opens and
`setProgress(0)` as each spectrum begins; `startProgress(0, spectra, "storing
DTA2D file")` before the output is created and `setProgress(i)` per spectrum,
none for the TIC; `endProgress()` at the end (`DTA2DFile.h:72-248`, `:258-287`,
`:297-320`). A failure after the start leaves the section open, as the source's
exception does. They are the Release build's calls, call for call (tier 1,
`tests/progress_format_readers.rs`; see `docs/PROGRESS_LOGGER_SUPPORT.md`). The
other entry points run the same code and report nothing. `MS2File` owes no
call: its only one is commented out (`MS2File.h:52`), and the Release build
makes none.

## Validation and limits

Finite signed coordinates/intensities are supported. Intensity tokens parse
directly to f32, avoiding an extra f64 rounding step. Nonfinite values, nonzero
values that underflow to zero, minute-conversion overflow and nonfinite TIC sums
are checked errors. Parsing returns owned data only on success; replacement
leaves the destination unchanged on parse, stream or limit errors.

Both adapters expose the same `Limits` type. Defaults and hard ceilings are
256 MiB input, 1 MiB per line, 20 million lines, 100,000 raw spectra/groups,
10 million raw peaks and 512 MiB output. Limits can be lowered. Ignored headers,
pre-scan peaks and filtered data still count. A bounded line reader uses at most
one extra byte to detect an exceeded input/line limit; vectors have checked
allocation. Peak/scan limits bound retained storage independently of input size.

Writers validate the entire represented payload and count formatted output
before writing any bytes or creating/truncating a path. All writers flush and
report flush errors. An actual destination I/O failure can still leave partial
output; this is not a transactional filesystem writer.

[Tests](../tests/text_peak_lists.rs) cover four exact source fixtures, source
filter results, scientific rounding, malformed inputs, atomic replacement,
limits, metadata-loss guards, the source text and its refusal, the source mode,
the recorded document identity and flush failures. The
[provenance manifest](../tests/data/text_peak_lists_provenance.json) records
14 source hashes, six fixture hashes (two of them Release build output) and
derived versus literal expectations.

## Evidence

- **Tier 1, retained Release output.** `dta2d_store_release.dta2d` and
  `dta2d_store_tic_release.dta2d` are what the OpenMS4 Release build
  (`openms4-release-bc9cc12-c19e494-174b576`, ibminode06, conda-forge g++ 14.4.0)
  wrote from `DTA2DFile_test_1.dta2d` in `../oracle/progress-format-readers`,
  three repetitions and both logger modes byte-identical. The port writes them
  byte for byte; `tests/progress_format_readers.rs` now also compares the
  output sizes it used to exempt. The progress calls of all three operations
  are tier 1 as described under Progress.
- **Tier 3, class tests.** `DTA2DFile_test.cpp` and `MS2File_test.cpp`
  literals: spectrum counts, native IDs, values, filters and the recorded
  document identity.
- **Tier 4.** The limits, the refusal of inexact text and the source mode's
  handling of unrepresentable spectra are independently derived.
The first fixture TIC is exactly f32 `141649.671875`; the source class test's
short printed expectation is `141650`.
