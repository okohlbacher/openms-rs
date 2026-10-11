# Spectrum lookup

`metadata::spectrum_lookup` covers `METADATA/SpectrumLookup.h` and
`METADATA/SpectrumLookup.cpp` at SDK `bc9cc12`. It builds on
[`SPECTRUM_NATIVE_ID_PARSER_SUPPORT.md`](SPECTRUM_NATIVE_ID_PARSER_SUPPORT.md) and
runs its regular expressions through `concept::boost_regex::BoostRegex`.

Tests: [`tests/spectrum_lookup.rs`](../tests/spectrum_lookup.rs), every
`START_SECTION` of `SpectrumLookup_test.cpp` (12; the two reference-format
sections share a fixture) plus derivations from the source.

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `static const std::string& default_scan_regexp` | `DEFAULT_SCAN_REGEXP` (`=(?<SCAN>\d+)$`) | |
| `std::vector<boost::regex> reference_formats` | `pub reference_formats: Vec<BoostRegex>` | |
| `double rt_tolerance` | `pub rt_tolerance: f64`, default `DEFAULT_RT_TOLERANCE` (0.01) | |
| `SpectrumLookup()`, `~SpectrumLookup()` | `new`, `Default`, `Drop` | |
| private copy constructor/assignment (not implemented) | `Clone` | Harmless for a value that owns its indexes. |
| `bool empty() const` | `is_empty` | |
| `template <class C> void readSpectra(const C&, const std::string& scan_regexp = default_scan_regexp)` | `read_spectra(&mut self, &[MSSpectrum], &str) -> Result<Vec<String>>` | Returns the per-spectrum warnings the source logs. |
| `Size findByRT(double) const` | `find_by_rt(f64) -> Result<usize>` | |
| `Size findByNativeID(const std::string&) const` | `find_by_native_id(&str) -> Result<usize>` | |
| `Size findByIndex(Size, bool count_from_one = false) const` | `find_by_index(usize, bool) -> Result<usize>` | |
| `Size findByScanNumber(Size) const` | `find_by_scan_number(usize) -> Result<usize>` | |
| `Size findByReference(const std::string&) const` | `find_by_reference(&str) -> Result<usize>` | |
| `void addReferenceFormat(const std::string&)` | `add_reference_format(&str) -> Result<()>` | |
| deprecated `static extractScanNumber(native_id, boost::regex, no_error)` | deprecated `extract_scan_number_with` / `try_extract_scan_number_with` | Forward to `spectrum_native_id_parser`. |
| deprecated `static extractScanNumber(native_id, accession)` | deprecated `extract_scan_number` | Forwards. |
| deprecated `static getRegExFromNativeID`, `static isNativeID` | deprecated `regex_from_native_id`, `is_native_id` | Forward. |
| protected `regexp_names_` | `REGEXP_NAMES` (`INDEX0 INDEX1 SCAN ID RT`) | |
| protected `addEntry_`, `findByRegExpMatch_`, `setScanRegExp_`, the maps | private | |

`Exception::ElementNotFound` is `Error::InvalidRange` (the convention
`format::mascot_xml` set for the same class), `IllegalArgument` and
`ParseError` are `Error::InvalidValue`, `MissingInformation` is
`Error::MissingInformation`.

## Preserved source conventions

- **Later spectra win** for a duplicate retention time, native ID or scan
  number, as the source's map assignments overwrite.
- **Retention-time look-up** takes the nearest time at or below the query only
  when it is strictly closer than the nearest above; a tie goes up. A distance
  equal to the tolerance matches.
- **Scan numbers** are extracted with the expression given (the default takes
  the digits after the last `=` at the end of the ID); every failure gets the
  source's warning text; `-1` is never stored; an empty expression disables
  scan look-up.
- **A reference format needs one of the five group names** as `?<NAME>`; the
  first format that matches decides, even when it carries nothing usable, and
  its first non-empty group in the order `INDEX0`, `INDEX1`, `SCAN`, `ID`, `RT`
  is used. Numbers convert as `StringUtils::toInt32`/`toDouble`
  (`concept::string_conversion`), so a conversion failure is an error, not a
  miss.
- **Index 0 counted from one is not found**, which the source reaches through
  unsigned wrap-around.

## Native differences

- **`read_spectra` is atomic.** The source clears its indexes and sets the new
  count before it validates the scan expression, so a refused expression leaves
  it half-reset; here nothing changes on error.
- **Non-finite retention times are refused.** The source keys a
  `std::map<double, Size>` with them, whose ordering a NaN breaks. A NaN query
  is simply not found, as in the source.
- **Negative numbers** from a caller-supplied expression: the source converts a
  negative index or scan number to a huge `Size`; here a negative index or scan
  in a reference is not found, and a negative scan number other than `-1` is
  stored under its signed value (still warned about, as the source does), where
  only a wrapped `Size` query could reach it in the source.
- **`MAX_SPECTRA`** (5·10^7) bounds one `read_spectra`, and every search has
  `BoostRegex`'s work budget; exceeding either is an error.

## Other implementations in the crate

`format::pepxml::SpectrumIndex` and `format::mascot_xml::SpectrumTitleLookup`
each carry the subset their reader needs (written before this module existed).
Collapsing them onto `SpectrumLookup` is a follow-up for the owners of those
formats.

## Consumers

NovorAdapter and PercolatorAdapter at the TOPP pin. Neither is ported yet.

## Evidence

Tier 3: every class-test literal transcribed unchanged; the rules above that go
beyond the class test are derived from `SpectrumLookup.cpp` and marked at
`source_rules_beyond_the_class_test`. No C++ was executed.
