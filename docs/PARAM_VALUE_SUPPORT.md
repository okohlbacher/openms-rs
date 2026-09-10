# Parameter values

`openms::param::ParamValue` implements the seven storage alternatives in the pinned Core SDK's `ParamValue`: `Empty`, `String(String)`, `Integer(i64)`, `Float(f64)`, `StringList(Vec<String>)`, `IntegerList(Vec<i32>)`, and `FloatList(Vec<f64>)`. It is separate from identification metadata and accepts NaN, positive infinity, and negative infinity. Stored strings may contain Unicode, embedded NUL, commas, or line breaks. A format adapter can impose additional restrictions at its own boundary.

```rust
use openms::param::ParamValue;

fn main() -> openms::Result<()> {
    let threshold = ParamValue::from(47.11);
    assert_eq!(threshold.to_text(true)?, "47.109999999999999");
    assert_eq!(threshold.to_text(false)?, "47.11");
    assert!(ParamValue::from("true").to_bool()?);
    Ok(())
}
```

## Construction, conversion, and ownership

`Default` and `ParamValue::EMPTY` are empty values. An empty string or empty list is a different value. `value_type()` reports the source type discriminant through `ParamValueType`; `is_empty()` tests only the empty alternative.

Owned variants, `From`, Rust assignment, moves, and `Clone` replace the source constructors and assignment operators. Signed integer types and `u8`/`u16`/`u32` fit the signed storage. `TryFrom<u64>`, `TryFrom<usize>`, and `TryFrom<isize>` check the storage range. A source `float` is promoted to binary64; it is not reparsed from decimal. Rust has no portable C++ `long double` counterpart: callers must explicitly supply the binary64 value used by this source class.

`as_str`, `as_string_list`, `as_integer_list`, and `as_float_list` return borrowed values of exactly the requested storage type. `to_string_vector`, `to_int_vector`, `to_double_vector`, and `TryFrom<&ParamValue>` provide checked owned conversions. A string such as `"1,2"` does not implicitly become a list. A floating-point value does not implicitly become an integer.

The `to_i*`/`to_u*` methods cover all native 8/16/32/64-bit and pointer-sized integer targets. Narrowing and negative-to-unsigned conversions fail when out of range. These checked errors replace source narrowing that can truncate or depend on the C++ host's integer widths. `to_f64` accepts integer or floating-point storage; `to_f32` preserves explicit NaN/infinity and ordinary rounding, but rejects finite overflow. Integer-to-f32 conversion is direct, avoiding a second rounding through f64. Finite underflow to f32 zero remains allowed.

`to_bool` accepts only the exact strings `"true"` and `"false"`. `to_char` returns `Result<Option<&str>>`: empty becomes `None`, a string becomes a complete borrowed string, and other types fail. It does not expose a raw pointer or truncate embedded NUL. Source floating-point casts accidentally read an inactive union field for nonnumeric strings/lists; the native API rejects those casts.

## Formatting and comparison

`to_text(full_precision)` reproduces source `toString`. Empty becomes `""` despite the source header's contradictory exception comment. String lists are unquoted and unescaped (`[a, b]`); this display is not a reversible list encoding. Integer text is decimal. Floating-point output preserves negative zero and uses `NaN`, `inf`, and `-inf` for special values.

For nonzero magnitudes below 0.01 or at least 10,000, full output uses the shortest scientific representation that round-trips; low precision uses three fractional digits. Other magnitudes use 15 fractional digits at full precision or three at low precision. Trailing zeros are removed while retaining one digit after the decimal point. Scientific exponents have at least two digits and omit the positive sign. These are the current implementation's thresholds, superseding an older header comment. Rust's standard floating-point formatter performs the decimal conversion; source literal and boundary tests pin the resulting text.

`to_stream_text()` separately provides the source output operator's default classic-locale behavior: six significant digits, with ordinary `5.0` displayed as `5`. The source operator delegates to the caller's stream, despite a comment promising full precision. Arbitrary mutable C++ stream locale/flags are not modeled; callers needing other presentation can format borrowed scalar/list values using Rust's formatting facilities.

Equality compares storage types and complete values. Consequently integer 5 differs from floating 5.0; NaN is unequal to itself. There is no `Eq` or `Ord` implementation. `source_less` and `source_greater` exactly preserve the source comparisons: differing types are neither less nor greater, and lists compare **length only**. Rust `PartialOrd` preserves those less/greater results while returning `None` for equal-length unequal lists, avoiding a false `Equal` result inconsistent with equality.

`Hash` normalizes signed zero and retains every other floating-point bit. `source_hash64()` exposes the source's FNV-1a and hash-combine recurrence for a 64-bit little-endian host; it is not a portable source `size_t` serialization. Equal values have equal hashes. NaN storage prevents treating arbitrary `ParamValue` as an `Eq` key in Rust hash collections. Explicit type/value access remains available for applications choosing their own key policy.

## Bounds and evidence

Checked cloning, owned conversions, formatting, and `source_hash64()` preflight values before their allocation/traversal. Limits are one million list elements, 64 MiB logical value/operation allocation, and 50 million charged work units. Formatting includes conservative temporary numeric-buffer accounting, so a large list can reach the work bound before the element cap. Shared parameter-tree formatting uses the same cumulative work and allocation budget across values. Borrowed accessors, ordinary enum construction, `Clone`, equality, comparison traits, and `Hash` retain normal Rust caller-owned behavior; they do not secretly impose fallible resource checks.

The eight integration tests in [param_value.rs](../tests/param_value.rs) and two private budget/formatter tests cover source constructor/conversion and formatting literals, all storage alternatives, source comparison quirks, IEEE special values, native narrowing errors, ownership, and resource failures. Eight exact hash constants were independently calculated from explicit little-endian input bytes and the source recurrence. [param_value_provenance.json](../tests/data/param_value_provenance.json) records the pinned source hashes, literal assertion locations, derived inputs, and methodology. No C++ build or execution was used.
