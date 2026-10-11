# ChromeleonFile support

Port of `FORMAT/ChromeleonFile.h` / `ChromeleonFile.cpp` (core `bc9cc12`),
Chromeleon HPLC text exports read into one chromatogram. Rust:
`src/format/chromeleon.rs` (header expressions through the crate's
`BoostRegex` facade; lines through `src/format/line_reader.rs`). Tests:
`tests/chromeleon_file.rs` (11; the mzML round trips need `mzml`).
Provenance: `tests/data/chromeleon_file_provenance.json`.

## API mapping

| Source member | Rust |
|---|---|
| `ChromeleonFile()`, `~ChromeleonFile()` | `ChromeleonFile::new()` / `Default`, drop |
| `load(filename, MSExperiment&) const` | `load(path) -> MSExperiment`, `load_with_limits`, `read(reader)`, `read_with_limits` |
| `removeCommasAndParseDouble(std::string&) const` | `remove_commas_and_parse_double(&mut String) -> Result<f64>` |

## Preserved source conventions

- Header keys and their targets are the source's (`Injection` -> `mzml_id`,
  `Channel` -> `acq_method_name`, `Processing Method` -> instrument software
  name, `Instrument Method` -> instrument name, `Injection Date/Time`,
  `Detector`, `Signal Quantity/Unit/Info` -> meta values), matched in the
  source's order over the whole line, later lines overwriting; a value with a
  form feed does not match (`.` under `no_mod_s`), as in the source.
- `Raw Data:` or `Chromatogram Data:` ends the header; the caption after it
  is skipped; three tab-separated fields give (time, intensity) with the step
  unread; the time stays in **minutes**, as the source stores it.
- Commas are removed from the caller's string even when conversion fails; the
  conversion is `StringUtils::toDouble`.
- The experiment always receives exactly one chromatogram.

## Native differences

- A NaN, infinite or `f32`-overflowing point is `Error::Parse`; the source
  stores it and later consumers refuse the map.
- A non-UTF-8 header value is `Error::Parse`; the source stores the bytes.
- Resource ceilings: `ReadLimits` (line length, input size, 100 million
  points).
- The source's `FileNotFound`/`FileNotReadable`/`IOException` are
  `Error::Io` with the matching kind.
- Writing the loaded experiment as mzML needs the instrument name cleared: the
  crate's mzML writer refuses a non-CV instrument name, where the C++ writer
  writes the generic `MS:1000031` and drops the name
  (`MzMLHandler.cpp:3961-3969`). The class test's round trip does so.

## Checked boundaries and evidence

- Tier 3: every `START_SECTION` of `ChromeleonFile_test.cpp` with its literals
  over the byte-identical upstream fixtures
  (`20171013_HMP_C61_ISO_P1_GA1_UV_VIS_2.txt`, `ChromeleonFile_commas.txt`,
  `ChromeleonFile_new_header.txt`).
- Tier 4: whole-line matching of markers and keys, a file without data
  marker, the `Chromatogram Data:` marker, malformed data lines, the point
  limit and a missing file. The non-finite refusal is not exercised by a
  test.
