# Spectrum native-ID parser

`metadata::spectrum_native_id_parser` covers `METADATA/SpectrumNativeIDParser.h`
and `METADATA/SpectrumNativeIDParser.cpp` at SDK `bc9cc12`. The class holds four
static functions and no state, so it becomes a module of free functions.

Tests: [`tests/spectrum_native_id_parser.rs`](../tests/spectrum_native_id_parser.rs),
one test per `START_SECTION` of `SpectrumNativeIDParser_test.cpp` (all six), plus
unit tests of the `toInt32` port and of every fixed expression compiling.

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class SpectrumNativeIDParser` | the module | |
| `static bool isNativeID(const std::string& id)` | `is_native_id(&str) -> bool` | The ten prefixes are public as `NATIVE_ID_PREFIXES`, `frame=` included. |
| `static std::string getRegExFromNativeID(const std::string& native_id)` | `regex_from_native_id(&str) -> &'static str` | Returns the source's Boost text unchanged, `(?<GROUP>\d+)` included. |
| `static Int extractScanNumber(const std::string&, const boost::regex&, bool no_error = false)` with `no_error = false` | `extract_scan_number_with(&str, &BoostRegex) -> Result<i32>` | `Exception::ParseError` becomes `Error::InvalidValue`. |
| the same with `no_error = true` | `try_extract_scan_number_with(&str, &BoostRegex) -> Result<Option<i32>>` | The `-1` sentinel becomes `None`. |
| `static Int extractScanNumber(const std::string&, const std::string& native_id_type_accession)` | `extract_scan_number(&str, &str) -> Result<ScanNumberExtraction>` | `ScanNumberExtraction { scan: Option<i32>, warning: Option<String> }`; `source_value()` gives the source's `Int`. |

Native additions: `NATIVE_ID_PREFIXES`, `SCAN_ACCESSIONS`, `FILE_ACCESSIONS`,
and the named accessions `WIFF_ACCESSION`, `INDEX_ACCESSION`, `SCAN_ID_ACCESSION`,
`SPECTRUM_ACCESSION`, `NUMERIC_ACCESSION` — the source's literals, named.

## Preserved source conventions

- **Regular expressions are Boost's**, executed by `concept::boost_regex::BoostRegex`,
  whose `tokens` is `boost::sregex_token_iterator`. "The last matching sub-group"
  is therefore exactly the source's: every non-overlapping match contributes its
  group 1 (the header asks for a group named `SCAN`; the code reads group 1
  whatever its name), and the last token is converted. `\d` and `\s` are ASCII,
  as in Boost's C locale.
- **The WIFF expression keeps its two groups both named `GROUP`**, which Boost
  allows; the token iterator addresses them as 1 and 2, and the last match's
  pair gives `cycle * 1000 + experiment`. The experiment is converted and range
  checked before the cycle is converted, as in the source.
- **`index=` under MS:1000774 returns `index + 1`**, the pepXML convention the
  source notes; `regex_from_native_id` plus the regex overload does not add one,
  which `USI::extractScanNumberFromNativeID` relies on.
- **Conversion is `StringUtils::toInt32`**: surrounding space, tab, newline and
  carriage return are skipped, one `+` may precede the number, the rest must be
  a complete `i32`. A digit run beyond `i32` is the source's `ConversionError`
  branch.
- **The three warning texts** the source writes to `OPENMS_LOG_WARN` (unknown
  accession, fewer tokens than groups, conversion failure) are returned verbatim
  in `ScanNumberExtraction::warning`; the crate has no global log stream.
- **An experiment number of 1000 or more is an error**, `Error::InvalidValue`,
  where the source throws `Exception::InvalidValue` out of its `ConversionError`
  handler.

## Native differences

- **Overflow is checked.** `cycle * 1000 + experiment` and `index + 1` are `int`
  arithmetic in the source, so a cycle above 2147483 or an index of `INT_MAX` is
  signed overflow, which C++ leaves undefined. Here the arithmetic is checked
  and returns `scan: None` with a warning naming the overflow. The same rule as
  `format::mascot_generic`'s MGF writer.
- **`None` for the `-1` sentinel.** A caller-supplied expression that captures a
  literal `-1` reads as failure in the source and as `Some(-1)` here.
- **Regular-expression work is bounded** by `BoostRegex`'s budget; a search that
  exceeds it is `Error::InvalidValue`. The source has no such limit. The fixed
  expressions need no backtracking and never reach it.

## Other implementations in the crate

Three consumers carried their own subsets before this module existed, and still
do: `format::mascot_generic::extract_scan_number` (the accession table, by hand),
`format::pepxml::SpectrumIndex::is_native_id` (which additionally accepts
`frame=`, now also the source's behaviour), and
`format::percolator_infile::ScanNumberPattern::for_native_id` (the prefix only,
written when the crate had no regular-expression engine). Collapsing them onto
this module is a follow-up for the owners of those formats.

## Evidence

Tier 3: all six class-test sections are transcribed with unchanged literals.
Boundaries beyond the class test (warning texts, the experiment limit, overflow,
too few WIFF tokens) are derived from the source branches and say so at the test.
No C++ was executed for this module.
