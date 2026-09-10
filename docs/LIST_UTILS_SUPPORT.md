# Primitive and string list helpers

`data_structures::list` and `data_structures::string_list` cover the public operation families of pinned OpenMS `ListUtils` and `StringListUtils`. They are functions on native slices/iterators and strings, with no collection wrapper or added dependency.

## API mapping

| Source operation | Native operation |
| --- | --- |
| `create<T>(text, splitter)` | `list::create::<T>(text, splitter: u8)` |
| `create<T>(vector<string>)` | `list::create_from_strings::<T>(&values)` |
| `toStringList` | `list::to_string_list(iterable)` |
| Vector/container `concatenate` | `list::concatenate(iterable, glue)` |
| Generic `contains` | `list::contains(slice, &element)` |
| Double `contains` | `contains_f64` (source default `1e-5`) or `contains_approx` |
| String `contains` with case option | `contains_string` and `CaseSensitivity` |
| `getIndex` | `get_index`, returning the first `Option<usize>` |
| Mutable/const container prefix/suffix searches | `string_list::search_prefix` / `search_suffix` |
| Mutable/const iterator-range searches | `search_prefix_in` / `search_suffix_in` with checked half-open ranges |
| `toUpper`, `toLower` | `to_upper`, `to_lower` on mutable string slices |

Pass `b','` and `""` explicitly for the source default delimiter and glue. Search functions return positions in the supplied slice; range variants return positions in the original full slice. `None` replaces iterator-end and `-1` sentinels, avoiding the C++ narrowing of large indices to `Int`.

## Literal source conventions

`create` splits on a single byte, ignores quote structure, preserves empty trailing fields, and returns an empty vector for empty input. String conversion never trims, including the vector overload whose generic header comment says otherwise. Numeric conversion trims exactly space, tab, LF and CR and requires the entire remaining token. The four conversions actually defined by the source are `String`, `i32`, `f32` and `f64`; `ListParse` permits caller-defined extensions. There is no implicit conversion of arbitrary Rust types.

Integer parsing preserves the source's one-leading-plus removal: `+-1` is accepted, `++1` is not. Floating parsing supports decimal/scientific notation and explicit case-insensitive NaN/infinity spellings. The source's separate unsigned NaN branch accepts arbitrary parenthesized payload text; signed NaN payloads use ASCII letters, digits and underscore. Overflows and residual text fail. Binary32 values parse directly to `f32`, avoiding double rounding. Representable subnormals and true signed-zero inputs are accepted. A nonzero decimal mantissa that underflows to zero is rejected, matching the source standard `from_chars` path; the libc++ fallback accepts that case, so this native policy deliberately selects the stricter portable boundary. Hexadecimal floats accidentally accepted by that fallback are rejected, consistent with the source `from_chars` general-format path. Cross-platform C++ library parsing identity outside these documented rules is not claimed.

`ListFormat` covers strings, primitive integers, ASCII `char`, `bool`, `f32`, `f64`, `ParamValue` and the native `MetaValue` counterpart of DataValue. Units are not appended. Floats use the shared source StringUtils formatter: ordinary fixed notation has 15 fractional digits for `f64` or 6 for `f32`, with trailing zero removal and a mandatory fractional digit; nonzero magnitude outside `[0.01, 10000)` uses shortest scientific notation. Exponents omit `+` and contain at least two digits. NaN/inf use `NaN`, `inf`, `-inf`. Iterator order is retained. C++ `long double`, arbitrary raw bytes and locale-specific non-ASCII character conversions have no direct portable Rust equivalent; a caller may supply its own `ListFormat`.

Floating containment uses the literal strict test `abs(value - target) < tolerance`. It does not normalize a negative, zero, NaN or infinite tolerance; equal infinities do not match. Exact `get_index` never applies this approximate comparison. Case-insensitive containment and case conversion use deterministic ASCII/C-locale rules while retaining other Unicode characters unchanged. Prefix/suffix searches remain case-sensitive. Optional trimming applies to both query and each candidate, never mutating either; a trimmed-empty query matches the first element.

## Checked boundaries and evidence

Allocating helpers cap one million items and a conservative 64 MiB of retained storage (two vector-slot and text-payload allowances). They precharge known input/output counts and copies, check overflow/allocation failure, and return fresh results atomically. Formatting custom trait implementations remains the caller's responsibility; the enclosing helper bounds returned content. Nonallocating borrowed searches and in-place ASCII case conversion remain ordinary linear operations on caller-owned slices. Invalid ranges, numeric conversions and UTF-8 byte splits return errors. A source byte delimiter that would split a Unicode character is not silently replaced with a Unicode delimiter.

[Primitive tests](../tests/list_utils.rs) retain all distinct ListUtils class-test literals and add parser, source formatter, direct-binary32, strict tolerance, ownership and limit cases. [String tests](../tests/string_list_utils.rs) cover the source prefix/suffix positions, range exclusions, trim/case literals and checked boundaries, reusing the exact [TextFile source fixture](../tests/data/text_file_source.txt) whose eleven lines equal the StringListUtils test's inline list. The [manifest](../tests/data/list_utils_provenance.json) pins headers, implementations, tests and that reused fixture. Derived edge expectations come from the source operations and explicit IEEE/UTF-8 arithmetic; no C++ execution or Rust-generated goldens are claimed.
