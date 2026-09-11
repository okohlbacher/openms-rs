# DateTime

`openms::data_structures::DateTime` ports the complete class-specific value/operation surface of OpenMS4-core `DateTime.h/.cpp` at revision `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. The module is also public as `data_structures::datetime`. It reuses the existing chrono clock dependency. No experiment aggregate or mzML timestamp adapter is changed.

## State and public operations

The value is fixed-size, privately storing signed year/month/day/hour/minute/second/millisecond fields and an independent validity flag. `Default`, `Copy`, `Clone`, assignment and `Eq` map source construction/copy/move/equality. `clear` resets every field to zero and validity to false. `is_null` means the flag is false. It is not a test for calendar completeness: setting only a time marks a zero-date value valid, as in the source.

| Source operation | Native operation |
| --- | --- |
| `set(string)` | `set(&str)`; `parse`/`FromStr` construct through the same automatic parser |
| numeric `set(month, day, year, hour, minute, second)` | `set_components`, using six u32 arguments in the same order |
| string/numeric `setDate`, `setTime` | `set_date`, `set_date_components`, `set_time`, `set_time_components` |
| component-output `get`, `getDate`, `getTime` | `components`, `date_components`, `time_components`, plus `millisecond` |
| string `get`, `getDate`, `getTime` | `get`, `date_string`, `time_string`; Display uses `get` |
| default ISO `toString` | `iso_string`; explicit formats use `format` |
| `fromString(input, format)` | `from_format`; pass `DATETIME_FORMATS[0]` for the source default ISO format |
| `isValid`, `isNull`, `clear` | `is_valid`, `is_null`, `clear` |
| `addSecs(int)` | checked `add_seconds(i32) -> Result<&mut Self>` |
| `now`, `nowUTC` | `now`, `now_utc` |
| `operator<` | `source_less`; full Eq/Ord cannot consistently model the source's distinct comparison keys |
| `std::hash` | Rust Hash of the same ISO-millisecond rendered identity; no C++ numeric digest guarantee |

Component getters return signed i32 tuples in source month/day/year and hour/minute/second order. This makes pre-year-one arithmetic readable; C++ casts negative years into unsigned outputs. Numeric setters reject values outside the valid positive i32-year/calendar/time domain. They reset milliseconds only for the full six-component setter. Partial setters preserve milliseconds and untouched fields.

Equality includes every field and validity. `source_less` compares chronological fields and milliseconds while ignoring validity. An invalid default and a valid midnight-only value can therefore be unequal but neither less than the other. Invalid values hash as an empty string even if arithmetic changed their private fields; source hash has the same allowed collision. Rust hashes remain hasher-dependent.

## Parsing and rendering

`DATETIME_FORMATS` lists the seven exact source spellings: `yyyy-MM-ddThh:mm:ss`, `yyyy-MM-ddThh:mm:ss.zzz`, `yyyy-MM-dd hh:mm:ss`, `yyyy-MM-dd+hh:mm`, `yyyy-MM-ddThh:mm:ssZ`, `yyyy-MM-dd`, `hh:mm:ss`. No locale pattern language or timezone state is added. Invalid values format as empty strings even for unknown format text. Valid values reject unknown formats. The `get` family instead returns source zero-date/time fallback strings when invalid.

Automatic parsing preserves source branch priority: German dot/no-T, slash, ISO/legacy hyphen branches, then two legacy weekday/month-name fallbacks. ISO text containing a plus is cut at its first plus; offset text is ignored, without clock conversion. The legacy date-plus-time form uses plus as a separator. Date-only automatic parsing requires its source final-Z condition; explicit `from_format` can parse a plain date. The source's successful-assignment count means a missing final Z literal, trailing text and some odd prefixes are accepted in specific branches. An arbitrary three-byte weekday token is permitted, while English month abbreviations are case-sensitive.

The borrowed scanner reproduces source `%d`, `%3s`, literal and whitespace rules: optional single signs, ASCII C whitespace (including vertical tab), zero-or-more format whitespace, successful prefix parsing and NUL termination of the scanned string. Branch selection still sees the complete supplied string, including bytes after NUL, because source StringUtils uses std::string. Unicode digits/whitespace are not added as scanf alternatives. Integer conversion overflow rejects instead of reproducing C undefined conversion behavior.

Fractions count consecutive ASCII digits immediately after the first dot. Fewer than three digits scale the parsed integer, then signed remainder 1000 is applied. Thus `.4` becomes 400, `.46` becomes 460 and `.12345` becomes 345. Source scanf also permits a sign/whitespace before the integer; the raw digit count can then be zero and change normalization. Checked scaling overflow rejects instead of executing source signed-overflow undefined behavior.

`set` clears the old value before a syntax/calendar error, matching source error publication. Partial setters and the numeric full setter leave the old value intact on failure. `from_format` returns an invalid default value for unknown formats, invalid calendars, syntax/numeric/scaling failure; an input-size violation remains an operation error. This keeps invalid values distinct from failing bounded operations.

## Arithmetic and host behavior

`add_seconds` uses deterministic proleptic Gregorian UTC arithmetic, normalizes zero/partial date fields, and retains the independent validity flag and milliseconds. It never uses local DST. All i64 intermediates fit throughout the private i32-year domain; a resulting year outside i32 returns an error before mutation. For example, adding zero to default fields produces private date -1-11-30 while the value stays invalid; its fallback string remains all zeros.

The actual source delegates arithmetic to unchecked `timegm`/`gmtime` (Windows uses different CRT functions). The macOS reference run returned time_t -1 for input years below 1900 and the source then continued from 1969-12-31 23:59:59. This is a platform-dependent source failure, not a portable calendar rule. Native arithmetic deliberately corrects it. All 35 affected reference rows are retained and checked against independent Python calendar month-stepping expectations; they are not skipped or called C++ matches. Calendar years outside host libc's useful range remain supported within the checked native domain.

`now` reads local wall-clock time and `now_utc` reads UTC through the existing chrono clock support. Both set milliseconds to zero and validity true. Other operations do not mutate or depend on process locale/timezone. Clock values depend on the actual host time; they are not fixture outputs.

## Bounds and evidence

Supplied date/time input is capped at `MAX_DATETIME_INPUT_BYTES` (1 MiB), including NUL suffix bytes. This preflight occurs before scanning or clearing state, a documented native exception to source clear-first publication. Parsing uses fixed arrays and borrowed slices, with a fixed number of linear passes; no allocation proportional to input is performed. Formatting produces less than 64 bytes for every reachable state. Calendar arithmetic is constant-size. Ordinary Copy/Hash/format operations need no heap-size configuration.

[The tests](../tests/datetime.rs) cover the complete source public operation families, 19 published class-test string rows, valid/invalid/partial state, formatting, C-string/scanner corners, full numeric range, fractional overflow, Gregorian reversibility, clock precision and hashing. The [raw C++ fixture](../tests/data/datetime_cpp_probe.tsv) contains 301 executed rows: **266 exact matches and 35 corrected host-dependent rows**, with seven output formats and component/fallback strings. [Independent corrected rows](../tests/data/datetime_calendar_corrections.tsv) use Python `calendar.monthrange` with month-by-month stepping, not the native March-era day-number formula.

The source DateTime.cpp, DateTime.h and HashUtils.h were compiled unmodified in a small UBSan-enabled probe. Four narrow include adapters supply export macros, UInt, exceptions and the exact used StringUtils has/prefix operations; integer-to-text exception diagnostics use std::to_string. This validates the DateTime scientific body, not the complete SDK dependency stack, full source exception messages, clocks or numeric hash ABI. [Probe instructions and code](../tools/datetime_probe/README.md) and [provenance](../tests/data/datetime_provenance.json) record compiler, source, adapter, binary and fixture hashes. The first differential failure that exposed macOS behavior is retained in the integration validation record. Root independently reviewed the complete native implementation and test policy with no remaining findings.

Two confirmed upstream defects are documented with reproducible inputs and actual source diagnostics in [the C++ issue evidence](../tests/data/datetime_cpp_issues.json): unchecked macOS timegm failure and signed overflow in all three fractional-scaling branches. The sanitizer-abort cases are separate from the 301 ordinary probe rows.
