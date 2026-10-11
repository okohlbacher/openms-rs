# Universal Spectrum Identifier

`metadata::usi` covers `METADATA/USI.h` and `METADATA/USI.cpp` at SDK `bc9cc12`.
It builds on [`SPECTRUM_NATIVE_ID_PARSER_SUPPORT.md`](SPECTRUM_NATIVE_ID_PARSER_SUPPORT.md).

Tests: [`tests/usi.rs`](../tests/usi.rs), every `START_SECTION` of `USI_test.cpp`
(36 sections, the ten accessor sections gathered in one test).

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `enum class IndexType { SCAN, INDEX, NATIVEID }` | `IndexType::{Scan, Index, NativeId}` | Declaration order kept; `Default` is `Scan`. |
| `USI()` | `Usi::default()` | Every string empty, `Scan`. |
| `USI(collection, ms_run, index_type, index, interpretation = "")` | `Usi::new(...)` | Rust has no default arguments; pass `""`. |
| `explicit USI(const std::string&)` | `Usi::parse(&str) -> Result<Usi>`, `FromStr` | `Exception::ParseError` is `Error::Parse` with the source text `Invalid USI format in: <text>` (executed, below). |
| copy/move constructors and assignments, `~USI()` | `Clone`, moves, `Drop` | |
| `operator==`, `operator!=` | derived `PartialEq`/`Eq` | All five fields. |
| `operator<` | derived `PartialOrd`/`Ord` | Field order collection, MS run, index type, index, interpretation — the source's comparison chain; strings compare byte-wise in both. |
| `isValid()` | `Usi::is_valid()` | |
| `static isValidUSI(const std::string&)` | `Usi::is_valid_usi(&str)` | |
| `static std::optional<USI> tryParse(const std::string&)` | `Usi::try_parse(&str) -> Option<Usi>` | |
| `getCollection`/`setCollection`, `getMSRun`/`setMSRun`, `getIndexType`/`setIndexType`, `getIndex`/`setIndex`, `getInterpretation`/`setInterpretation` | public fields `collection`, `ms_run`, `index_type`, `index`, `interpretation` | The setters validate nothing, so a field is the same contract. |
| `hasInterpretation()` | `Usi::has_interpretation()` | |
| `toString()` | `Usi::to_usi_string()`, and `Display` | Empty for an invalid USI. |
| `fromString(const std::string&) -> bool` | `Usi::set_from_string(&mut self, &str) -> bool` | Renamed because a `from_` method taking `self` breaks the Rust naming convention. |
| `static createFromScanNumber(dataset_id, filename, scan_number, interpretation = "")` | `Usi::from_scan_number(&str, &str, i32, &str)` | |
| `static createFromNativeID(dataset_id, filename, native_id)` | `Usi::from_native_id(&str, &str, &str) -> Result<Usi>` | `Result` only for the regular-expression budget. |
| `static std::optional<int> extractScanNumberFromNativeID(const std::string&)` | `Usi::extract_scan_number_from_native_id(&str) -> Result<Option<i32>>` | As above. |
| `static indexTypeToString(IndexType)` | `IndexType::as_str()`, `Display` | The unreachable `default: "scan"` arm has no counterpart. |
| `static indexTypeFromString(const std::string&)` | `IndexType::parse(&str)`, `FromStr` | `Exception::InvalidValue` is `Error::InvalidValue` with the source text (`the value '<text>' was used but is not valid; Invalid USI index type`). |
| `static extractBasename(const std::string&)` | `Usi::extract_basename(&str) -> String` | |
| `static getCVAccession()`, `static getCVName()` | `CV_ACCESSION`, `CV_NAME` | |
| (none) | `Usi::LOCAL_COLLECTION` | Native: the `"local"` default argument of the `PeptideIdentification::buildUSI` overloads. |
| protected `USI_PREFIX`, `CV_ACCESSION`, `CV_NAME` | `USI_PREFIX`, `CV_ACCESSION`, `CV_NAME` | Public constants. |
| `operator<<(std::ostream&, const USI&)` | `Display` | Writes `to_usi_string()`. |

## Preserved source conventions

- **Parsing** requires the `mzspec:` prefix, splits the rest on every `:` keeping
  empty pieces, needs at least four pieces, reads the index type
  case-insensitively and joins any further pieces back with `:` as the
  interpretation, so `EM[UNIMOD:35]K/2` survives. The result must be valid.
- **A rejected parse can leave the value partly filled**, as in the source:
  fields are reset, then assigned as they are read, and an unknown index type
  stops after collection and MS run. The boolean is the answer, not the state.
- **A trailing `:` is an empty interpretation**, so `to_usi_string` drops it.
- **The index type is written back canonically** (`scan`, `index`, `nativeId`)
  whatever case it was read in.
- **`extract_scan_number_from_native_id`** first applies, for a native ID, the
  expression `regex_from_native_id` picks and takes its last match — so
  `index=100` is `100`, not the accession overload's `101` — and otherwise
  converts the whole string with `toInt32`. A WIFF ID has no prefix-specific
  expression and yields its last number, the experiment; that is the source's
  behaviour and is pinned by a derived test.
- **`extract_basename`** strips `file://`, then the `/` before a Windows drive
  letter, then everything up to the later of the last `/` and `\`.

## Native differences

- **`Result` on the native-ID paths.** The regular expressions run through
  `BoostRegex`, whose searches have a work budget the source does not; a search
  over the budget is `Error::InvalidValue` rather than a silently different
  answer. The fixed expressions used here need no backtracking and do not reach
  it in practice.

## Consumers

`PeptideIdentification::build_usi` and `build_usi_with_mapping`
(`identification::peptide_identification`, see
[PEPTIDE_IDENTIFICATION_SUPPORT.md](PEPTIDE_IDENTIFICATION_SUPPORT.md)) return
this type. `TextExporter` is the one TOPP consumer at the pin. It is not ported
yet; this module gives it the complete type.

## Evidence

Tier 1 (executed differential): `../oracle/B8b-peptide-id/driver.cpp`, case
`usi_class`, compiled with libOpenMS's Release flags against the Linux x86_64
Release install `openms4-release-bc9cc12-c19e494-174b576` and run twice on
ibminode06 with byte-identical results. It covers 20 parse inputs (including
partial state and error text), nine index-type names, 15 basenames, 32 native
IDs, the three factory functions, the CV constants, a 7x7 `operator<` matrix
and `operator<<`; `tests/peptide_identification.rs` compares every row against
this module. Cases `usi_build` and `usi_mapping` exercise USI through
`buildUSI`. Fixture: `tests/data/peptide_identification_oracle.tsv`; hashes in
`tests/data/peptide_identification_provenance.json`.

Tier 3: every class-test literal of `USI_test.cpp` is transcribed unchanged in
`tests/usi.rs` (36 sections, the ten accessor sections gathered in one test)
and, again, in `usi_class_test_literals`. The partial-fill, trailing-colon,
ordering and WIFF cases of `tests/usi.rs` are derived from `USI.cpp` and say so
at their tests; the executed rows above cover the same behaviour.

Duplicate collapsed at integration: group B had ported the class a second time
as `identification::usi` (literal prefix matching instead of the regular
expressions). That module was removed; `identification` uses this one, and the
executed rows run against it.
