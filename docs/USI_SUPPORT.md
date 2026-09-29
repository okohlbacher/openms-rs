# Universal Spectrum Identifiers (USI)

`identification::usi` ports `src/openms/include/OpenMS/METADATA/USI.h` and
`src/openms/source/METADATA/USI.cpp` at revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. It was ported as the return type
of `PeptideIdentification::buildUSI` (see
[PEPTIDE_IDENTIFICATION_SUPPORT.md](PEPTIDE_IDENTIFICATION_SUPPORT.md)); its
only direct TOPP consumer is `TextExporter`.

A USI is `mzspec:<collection>:<ms_run>:<index_type>:<index>[:interpretation]`
(PSI-MS `MS:1003063`, <https://www.psidev.info/usi>).

Rust file: `src/identification/usi.rs`. Tests: `tests/peptide_identification.rs`
(seven USI tests) and one private unit test. The module lives beside
`PeptideIdentification` rather than in `metadata` because it reuses the
identification module's `StringUtils::toInt32` port; `identification` already
depends on `metadata`, and the reverse edge would close a cycle.

## API mapping

| C++ member | Rust | Notes |
|---|---|---|
| `enum class IndexType { SCAN, INDEX, NATIVEID }` | `UsiIndexType { Scan, Index, NativeId }` | Declaration order kept, so the derived `Ord` is the source enum order `operator<` uses. |
| `USI()` | `Usi::default()` | Empty, invalid, index type `scan`. |
| `USI(collection, ms_run, index_type, index, interpretation = "")` | `Usi::new(..., interpretation)` | Pass `""` for the default. Nothing validated, as in the source. |
| `explicit USI(const std::string&)` (throws `ParseError`) | `Usi::parse` / `FromStr` | `Error::Parse` with the source text `Invalid USI format in: <text>`. |
| copy/move constructors and assignments, destructor | `Clone`, moves, `Drop` | |
| `operator==`, `operator!=` | derived `PartialEq`/`Eq` | All five components. |
| `operator<` | derived `Ord`/`PartialOrd` | Collection, MS run, index type, index, interpretation; strings compared bytewise, so index `"10" < "2"`. |
| `isValid()` | `is_valid()` | Non-empty collection, MS run and index. |
| `static isValidUSI(s)` | `Usi::is_valid_usi` | |
| `static tryParse(s)` | `Usi::try_parse` | `Option<Usi>`. |
| `getCollection`/`setCollection` | field `collection` | |
| `getMSRun`/`setMSRun` | field `ms_run` | |
| `getIndexType`/`setIndexType` | field `index_type` | |
| `getIndex`/`setIndex` | field `index` | |
| `getInterpretation`/`setInterpretation` | field `interpretation` | |
| `hasInterpretation()` | `has_interpretation()` | |
| `toString()` | `Display` (`to_string()`) | `""` when invalid. |
| `fromString(s)` | `from_string(&mut self, s) -> bool` | Source partial state kept (below). |
| `static createFromScanNumber(dataset, file, int, interpretation = "")` | `Usi::from_scan_number` | Negative numbers are written with their sign, as in the source. |
| `static createFromNativeID(dataset, file, native_id)` | `Usi::from_native_id` | |
| `static extractScanNumberFromNativeID(id)` | `Usi::extract_scan_number_from_native_id` | `Option<i32>`. |
| `static indexTypeToString(t)` | `UsiIndexType::as_str` / `Display` | |
| `static indexTypeFromString(s)` (throws `InvalidValue`) | `UsiIndexType::parse` / `FromStr` | ASCII case-insensitive; `Error::InvalidValue` with the source message. |
| `static extractBasename(path)` | `Usi::extract_basename` | Borrows from the input. |
| `static getCVAccession()`, `getCVName()` | `Usi::CV_ACCESSION`, `Usi::CV_NAME` | |
| protected `USI_PREFIX`, `CV_ACCESSION`, `CV_NAME` | `Usi::PREFIX`, `Usi::CV_ACCESSION`, `Usi::CV_NAME` | Public constants. |
| `operator<<(ostream&, const USI&)` | `Display` | Same text as `toString`. |
| (none) | `Usi::LOCAL_COLLECTION` | Native: the `"local"` default argument of the `buildUSI` overloads. |

Every public member of the header is mapped; the status is `complete`.

## Preserved source conventions

1. `fromString` resets all components, requires the case-sensitive `mzspec:`
   prefix and splits the rest on every `:`. Fewer than four parts fail with
   everything reset. The fifth and later parts are joined back with `:`, so an
   interpretation such as `EM[UNIMOD:35]K/2` survives and a trailing `:`
   leaves it empty (`mzspec:a:b:scan:1:` is valid). The return value is
   `isValid()`, so an empty collection, MS run or index fails.
2. **Partial state.** The source assigns components as it parses. With an
   unknown index type the collection and MS run are already set while the index
   stays empty; the port leaves the same observable state. `parse` and
   `try_parse` never expose it.
3. `extractScanNumberFromNativeID` first applies the pattern
   `SpectrumNativeIDParser::getRegExFromNativeID` chooses when the identifier
   starts with one of the ten `isNativeID` prefixes, takes the digits of the
   **last** match and converts them with `toInt32`; otherwise, or when that
   fails, it converts the whole text with `toInt32` (blanks, one optional `+`,
   an optionally negative `i32`). Consequences the oracle confirms:
   `sample=1 period=1 cycle=42 experiment=3` gives 3 (the bare `(\d+)` pattern,
   not the WIFF cycle*1000+experiment rule), `index=7` gives 7 (not index+1),
   `scan=12 scan=34` gives 34, `frame=5` and `scan=-3` give nothing,
   ` 12345 ` and `+12` give 12345 and 12, `-5` gives -5, and any digit run that
   overflows `i32` gives nothing.
4. `extractBasename` removes `file://`, then one more `/` when the third
   character is `:`, and returns the text after the last `/` or `\`.
5. `toString` of an invalid USI is the empty string.

## Native differences

* The patterns are all a literal prefix followed by `(?<GROUP>\d+)`, so they are
  matched directly (leftmost, non-overlapping, greedy ASCII digits) rather than
  compiled; `usi_scan_prefix_matching_agrees_with_the_boost_patterns` compiles
  the source patterns with `concept::boost_regex::BoostRegex` and checks the
  two agree on every oracle input.
* `SpectrumNativeIDParser` itself is not ported as a type; the two helpers USI
  needs are private here, beside the partial copies in
  `format::{pepxml, mascot_generic, percolator_infile}` (see that header's
  ledger entry).
* Getters and setters are public fields.

## Evidence

Tier 1 (executed differential): `../oracle/B8b-peptide-id/driver.cpp`, case
`usi_class`, compiled with libOpenMS's Release flags against the Linux x86_64
Release install `openms4-release-bc9cc12-c19e494-174b576` and run twice on
ibminode06 with byte-identical results. It covers 20 parse inputs (including
partial state and error text), nine index-type names, 15 basenames, 32 native
IDs, the three factory functions, the CV constants, a 7x7 `operator<` matrix
and `operator<<`. Cases `usi_build` and `usi_mapping` exercise USI through
`buildUSI`. Fixture: `tests/data/peptide_identification_oracle.tsv`; hashes in
`tests/data/peptide_identification_provenance.json`.

Tier 3: `usi_class_test_literals` transcribes the values of `USI_test.cpp`.
All 36 `START_SECTION`s of that file are covered by the tier-1 cases or the
literals (the copy/move/assignment sections by `Clone`).
