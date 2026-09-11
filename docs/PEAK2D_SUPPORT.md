# Two-dimensional peak values

`kernel::Peak2D`, `MobilityPeak2D` and `RichPeak2D` follow their direct public
value surfaces in Core SDK `54a232fe2cae9c590d5c997fa49d20e7769860fb`.
[`peak2d.rs`](../src/kernel/peak2d.rs) uses owned native data; no C++ hierarchy or
generic DPeak/DPosition algebra is introduced.

`Peak2D` and `MobilityPeak2D` store public `position: [f64; 2]` and
`intensity: f32`. The array is the safe mutable two-coordinate position: callers
can read, replace or mutate either component directly. `new(first, mz,
intensity)` and `from_position([first, mz], intensity)` are available. `rt()` /
`set_rt()` or `mobility()` / `set_mobility()`, plus `mz()` / `set_mz()`, map the
source scalar accessors. Default values are zero. Ordinary Rust copy/move/
assignment implement the source value ownership operations. Construction and
field mutation retain arbitrary IEEE values; these are unvalidated raw records.

Both types expose `MZ = 1`, `DIMENSION = 2` and respectively `RT = 0` or `IM = 0`.
The four generic `short_dimension_name`, `full_dimension_name`,
`short_dimension_unit` and `full_dimension_unit` methods take a usize and return
`Result<&'static str>`. Their eight named RT/IM/MZ counterparts return strings
directly. Invalid indices return an error rather than source out-of-bounds
array access. Every source spelling and capitalization is preserved:

| Type/dimension | Short name | Full name | Short unit | Full unit |
| --- | --- | --- | --- | --- |
| Peak2D RT | RT | retention time | sec | Seconds |
| MobilityPeak2D IM | IM | ion mobility | ? | ? |
| Either MZ | MZ | mass-to-charge | Th | Thomson |

The mobility value does not contain an RT, unit field or spectrum attachment.
The question-mark unit is the source literal; it is not replaced with assumed
milliseconds. Spectrum mobility transport, mobility area selection, generic
range algebra and RichMobilityPeak2D remain separate work.

All four overload forms of the source intensity, coordinate and position
comparators map to ordinary native scalar/array comparisons, using getters or
public fields as needed. Array comparison is lexicographic in dimension order,
including partial ordering for NaNs. No arbitrary whole-peak ordering or `Eq`
is added. `PartialEq` compares every stored numeric field exactly. `Hash`
normalizes both signs of zero in every component; equal values therefore hash
equally. No stable digest, C++ digest or cross-version hash promise is made.

Display preserves the labels and spacing: `RT: 1 MZ: 2 INT: 3` and
`IM: 1 MZ: 2 INT: 3`. Numeric rendering uses Rust formatting, and a requested
formatter precision applies to all three numbers. C++ locale/stream precision
state is not emulated.

## Rich values and inherited state

`RichPeak2D { peak: Peak2D, metadata: MetaInfo, unique_id: u64 }` owns its metadata
and implements `Deref`/`DerefMut` to its plain point. Rich constructors and
`From<Peak2D>` / `From<&Peak2D>` start with empty metadata and ID zero. Clone
copies metadata independently. Equality and native Hash include the point,
metadata and unique ID. Metadata follows the existing native exact typed-value,
unit and string-key contract; C++ metadata registry numeric IDs or its floating
epsilon comparison are not newly reproduced here.

`replace_from_peak(plain)` implements plain-value assignment, clearing the
current metadata and ID. It returns the **entire previous RichPeak2D** by move;
no unrelated metadata is cloned or destroyed inside this constant-time method.
Dropping or reusing that returned ownership is the caller's operation. Rust's
borrow rules prevent passing an aliased base reference while mutably replacing
the same object; ordinary Rust self-value moves preserve the whole value.
Rich Display delegates to the point and does not print metadata or ID.

The two-accessor `HasUniqueId` implementation reuses the separate native
UniqueIdInterface/UniqueIdGenerator group. Import that trait to use validity,
clear, suffix-text parsing, ID-only swap and caller-owned generator assignment/
ensure operations. It introduces no second stored ID or duplicate generator.
This document covers RichPeak2D's direct value/conversion surface and that trait
integration, not complete parity with every inherited C++ metadata or global
singleton service.

Six focused [value tests](../tests/peak2d.rs) pin source constructor/mutation
values, every dimension string and named wrapper, comparator ordering literals,
Display labels, rich metadata ownership/assignment, ID state, units and signed
zero hashing. The [provenance manifest](../tests/data/peak2d_provenance.json)
records exact source hashes. No C++ reference program was built or executed.
