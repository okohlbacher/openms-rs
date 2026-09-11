# Peak value types: native header equivalents

`Peak1D` and `ChromatogramPeak` provide native equivalents for the complete public
operation groups in Core SDK `54a232fe2cae9c590d5c997fa49d20e7769860fb`
`KERNEL/Peak1D.h` and `KERNEL/ChromatogramPeak.h`. This describes value-type
behavior, not C++ ABI, numeric stream settings or hash-digest compatibility.

| Source operation | Native representation |
| --- | --- |
| One-dimensional position, double coordinate, float intensity | Public f64 `mz`/`rt` and f32 `intensity` fields |
| Default and position/intensity constructors | `Default` and `const new(position, intensity)` |
| Copy/move/assignment/destruction | `Copy`, `Clone`, ownership/assignment and ordinary drop |
| Position/intensity getters, setters and mutable references | Direct field access and borrowing |
| Equality/inequality | `PartialEq`, comparing both fields with ordinary float equality |
| Intensity, position and m/z/RT comparator overloads | Scalar `<` comparisons on the corresponding public field, including mixed peak/coordinate comparisons |
| Stream output | `Display`: `POS: <position> INT: <intensity>`, without a newline |
| `std::hash` specialization | Rust `std::hash::Hash` over position and intensity, normalizing signed zero |

`Display` uses Rust float formatting. Without a precision option, components use
Rust's default round-trippable display, rather than a C++ stream's default six
significant digits. An explicit precision applies Rust fixed-decimal formatting
to both numbers: `format!("{:.2}", Peak1D::new(123.456, 7.25))` gives
`POS: 123.46 INT: 7.25`. Source C++ locale, stream flags and width state are not
emulated. Other numeric formats remain available directly through the fields.
NaN/infinity and the sign of zero follow Rust formatting; constructing or
formatting a standalone peak does not validate its scientific suitability.

Hashing normalizes each zero-valued field to positive zero, because -0 and +0
compare equal. Other floating-point bit patterns are passed to the caller's
Rust `Hasher` in position-then-intensity order. No stable numeric digest,
cross-process value, or cross-language equivalence to the source FNV/golden-ratio
algorithm is promised. Hashes are not collision-free scientific identifiers.

Neither record implements `Eq` or a whole-record total order. NaN makes ordinary
float equality non-reflexive, and adding `Eq` would assert a false contract.
Consequently these records are not direct `HashMap`/`HashSet` keys merely because
they implement `Hash`; callers needing keys must choose a validated or explicit
bitwise representation with their intended equality. The source hash-container
test's finite-value use is represented here by equality/hash assertions, not by
pretending all IEEE floating values satisfy Rust `Eq`.

[`kernel_value_traits.rs`](../tests/kernel_value_traits.rs) covers the source
hash-test values, exact output labels, native precision, signed-zero equality and
hash normalization, individual changed fields, NaN partial equality, and the
scalar comparator mappings. Existing [`kernel.rs` tests](../tests/kernel.rs)
exercise defaults, copies, aligned sorting and selection using the same values.
[`kernel_value_traits_provenance.json`](../tests/data/kernel_value_traits_provenance.json)
pins the inspected headers, source stream implementations and class tests. No
C++ execution or reference digest golden was used.

Product now has a reviewed [native hash implementation](METADATA_HASH_SUPPORT.md).
Its full inherited source-model equivalence remains limited by the native
CVTerm/DataValue unit representation and metadata registry identity.
