# Native equivalents for aliases and standard containers

This review covers every declaration in seven small public headers from Core SDK
`6bfc0e4711105f4eda2fea86812a83af7c7e791f`; the `StandardTypes.h` and `DPeak.h`
sections were re-reviewed member by member at
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` (see
[kernel_gap_closures_provenance.json](../tests/data/kernel_gap_closures_provenance.json)). A `native_equivalent` ledger entry
means that Rust's existing types and standard library provide the header's data
and container operations.

## What the tier means, decided 2026-09-21

The tier drifted: it was defined here for seven alias-and-container headers, and
grew to carry **90** ledger entries, **89 of which have Rust files** — although
the definition means Rust's own types already do the job. Because
`tools/core_sdk_coverage.py` filters `native_equivalent` out of open work
*identically to* `complete`, that drift inflated the apparent completion of the
port. The lead's decision, in one sentence:

> **`native_equivalent` means the Rust language and its standard library
> discharge the header, with no port owed and no third-party crate required.**

Consequences, and why this line and not another:

- **A crate is not the standard library.** Where a crate discharges a header, a
  port exists and its implementation happens to be a dependency; it takes
  `complete` or `partial` on member coverage, and the crate is recorded in
  [THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md) under the
  crates-first policy. The distinction that matters to a planner is whether
  work is owed, and a dependency is a decision and a maintenance surface in a
  way `Vec<T>` replacing a container alias is not. `native_equivalent` bypasses
  the recorded-decision requirement, which is reason enough not to stretch it.
- **A header with any unported member is not `native_equivalent`**, whatever
  discharges the rest; it is `partial` with the member named.
- Applied on the day of the decision: `FORMAT/Base64.h` → `partial` (the
  `base64` crate, and the private SIMD encoder/decoder pair has no port),
  `FORMAT/ZlibCompression.h` → `complete` (the `flate2` crate, all four members
  covered), `CONCEPT/Macros.h` → `partial` (it flagged itself as outside this
  document's seven and concedes `OPENMS_THREAD_CRITICAL` is unported).
  `DATASTRUCTURES/MapUtilities.h` stays `native_equivalent`: a CRTP mixin of
  four pure traversal templates, no state, no crate.

**The 90 entries that predate this decision have not been re-audited against
it.** That audit is open work, not a claim this document makes; until it is
done, the `native_equivalent` count should be read as an upper bound on how
much is genuinely owed nothing. It does not promise C++ symbol names, binary layout,
allocator behavior, template inheritance, invalid iterator behavior, or identical
observability of moved-from objects. It also does not certify the complete APIs
of the domain classes stored in those containers.

| Header | Review result |
| --- | --- |
| `KERNEL/StandardTypes.h` | Native equivalent: three aliases |
| `KERNEL/DPeak.h` | Native equivalent: `Peak1D` / `Peak2D` are the two concrete `DPeak<D>::Type` results |
| `METADATA/PeptideIdentificationList.h` | Native equivalent: `Vec<PeptideIdentification>` |
| `DATASTRUCTURES/ExposedVector.h` | Native equivalent: `Vec<T>`, slices, iterators and `std::mem` |
| `DATASTRUCTURES/TypeAliases.h` | Native equivalent: three vector aliases |
| `CONCEPT/Types.h` | Partial: primitive/ASCII mappings exist; the entire precision/time surface is not implemented |
| `DATASTRUCTURES/ConstRefVector.h` | Partial: borrowed views exist; custom identity/capacity/comparison behavior is not represented |

No compatibility wrapper or production code was added for this review. The
existing Rust tests cited below exercise these representations in real library
operations. This is not a claim that every C++ class test was copied, nor a
separate validation of the Rust standard library.

## Standard spectrum names

`StandardTypes.h` (sha256 `434a332d5ff5388768e66b463b6e32814d6f44990d8fa1b998041259c3b4047a`)
declares four forward declarations and three typedefs. Every declaration:

| Source declaration | Rust counterpart | Notes |
| --- | --- | --- |
| `class MSSpectrum;` | `openms::MSSpectrum` | forward declaration; adds no member or runtime behavior |
| `class MSChromatogram;` | `openms::MSChromatogram` | forward declaration |
| `class Mobilogram;` | `openms::Mobilogram` | forward declaration only; the header defines **no** `Mobilogram` alias |
| `class MSExperiment;` | `openms::MSExperiment` | forward declaration |
| `typedef MSSpectrum PeakSpectrum;` | `openms::MSSpectrum` | not ported as a second Rust name: one name per type. Every consumer that takes a source `PeakSpectrum` (for example `comparison::BinnedSpectrum::new`) takes `&MSSpectrum` |
| `typedef MSExperiment PeakMap;` | `openms::MSExperiment` | same |
| `typedef MSChromatogram Chromatogram;` | `openms::MSChromatogram` | same |

Behavioural differences: none. A C++ typedef is the same type, and so is the
Rust mapping; the only difference is that Rust code spells the `MS*` name. This
review neither claims `Mobilogram` support through this header nor marks
`MSSpectrum` / `MSExperiment` / `MSChromatogram` themselves feature-complete.

The Rust types are exported by [lib.rs](../src/lib.rs) and implemented in
[kernel.rs](../src/kernel.rs). `StandardTypes_test.cpp` (sha256
`5daffef668fe54e6bce15348745502908830e253567005b4390fe7fd2d522c73`) has one
`NOT_TESTABLE` section plus the `GOOD_TYPEDEF` macro, which expands to a
construct section and a delete section and is invoked for `PeakSpectrum` and
`PeakMap` twice each (nine section invocations, five distinct sections);
`Chromatogram` is never instantiated by the source test. Section mapping into
[tests/kernel_aliases.rs](../tests/kernel_aliases.rs), tier 3 (source review):

| Source section | Rust test |
| --- | --- |
| `StandardTypes` (`NOT_TESTABLE`) | mapped: no assertion exists; `peak_spectrum_typedef_is_ms_spectrum` exercises the typedef target and asserts the source default `rt == -1.0`, `ms_level == 1` |
| `PeakSpectrum()` / `~PeakSpectrum()` (x2) | ported: `peak_spectrum_typedef_is_ms_spectrum` |
| `PeakMap()` / `~PeakMap()` (x2) | ported: `peak_map_typedef_is_ms_experiment` |
| (none) | extra: `chromatogram_typedef_is_ms_chromatogram` |

Self-audit (`StandardTypes.h`): 4 ported, 1 mapped-with-evidence, 0
mapped-without-evidence, 0 unaccounted.

## DPeak metafunction

`KERNEL/DPeak.h` (sha256 `82b6a53f382fdb9111a17ae1162ef0f33e90cfe34eee37e35bfa5e45c956e60f`)
declares a compile-time selector `DPeak<dimensions>::Type` with exactly two
specialisations; `DPeak.cpp` (sha256
`6b182e99286f8647d4493040c43eb0288a4036a1bdad10024b4d69b4e0477504`) only
instantiates four namespace-scope default objects. Rust has no use for the
metafunction: C++ code that is generic over dimensionality is written against
the two concrete types.

| Source declaration | Rust counterpart |
| --- | --- |
| `template <UInt dimensions> struct DPeak {}` (primary template, no members) | not ported: an unspecialised `DPeak<N>` has no `Type` and is a compile error in C++; nothing to represent |
| `DPeak<1>::Type` = `Peak1D` | `openms::Peak1D` ([kernel.rs](../src/kernel.rs)) |
| `DPeak<2>::Type` = `Peak2D` | `openms::kernel::Peak2D` ([kernel/peak2d.rs](../src/kernel/peak2d.rs)) |
| `DPeak.cpp` globals `default_dpeak_1`, `default_dpeak_1_type`, `default_dpeak_2`, `default_dpeak_2_type` | not ported: translation-unit anchors with no API; `Peak1D::default()` and `Peak2D::default()` are the observable values, all zero |

Behavioural differences: none at runtime. Representation differs as documented
for the concrete types: `Peak1D::PositionType` is `DPosition<1>` in C++ and
`Peak1D` stores `mz: f64` directly; `Peak2D` stores `position: [f64; 2]`
([PEAK2D_SUPPORT.md](PEAK2D_SUPPORT.md)).

`DPeak_test.cpp` (sha256
`32f56ebcc47e75bb040d9b4848122104c41eaf40577c75e968893492c7dd2185`) has four
sections, all construct/delete. Mapping into
[tests/kernel_aliases.rs](../tests/kernel_aliases.rs), tier 3:

| Source section | Rust test |
| --- | --- |
| `DPeak()` / `~DPeak()` (`DPeak<1>::Type`) | ported: `dpeak_1_type_is_peak1d` (`Peak1D::default() == Peak1D::new(0.0, 0.0)`) |
| `[EXTRA]DPeak()` / `[EXTRA]~DPeak()` (`DPeak<2>::Type`) | ported: `dpeak_2_type_is_peak2d` (`position == [0.0, 0.0]`, `DIMENSION == 2`) |
| (`DPeak.cpp` globals) | extra: `dpeak_cpp_globals_have_only_the_two_concrete_defaults` |

Self-audit (`DPeak.h`): 4 ported, 0 mapped-with-evidence, 0
mapped-without-evidence, 0 unaccounted.

## Peptide-identification collections

`Vec<openms::identification::PeptideIdentification>` is already the collection
used by spectrum/feature metadata, ID filtering, protein inference and format
transport. Algorithms that only need a view accept `[PeptideIdentification]`
slices. The element type is `Clone + Default + PartialEq`; its internal domain
API is a separate review from the list container.

Every added declaration in `PeptideIdentificationList.h` is covered:

| Source operation | Native operation |
| --- | --- |
| `EXPOSED_VECTOR_INTERFACE(PeptideIdentification)` and inherited constructors | The `Vec<T>` operations listed below |
| Constructor from `const std::vector&` | `vec.clone()` or `slice.to_vec()` |
| Constructor from `std::vector&&` | Ownership transfer of a `Vec`; `std::mem::take(&mut vec)` when the source binding must remain usable |
| Initializer-list constructor | `vec![a, b]` or `Vec::from([a, b])` |
| Copy assignment from vector | `destination.clone_from(&source)` or assignment of `source.clone()` |
| Move assignment from vector | Ownership assignment; `std::mem::take` where needed |
| Initializer-list assignment | `destination = vec![a, b]` |

Source distinction: the constructor taking `std::vector&&` uses move iterators,
so it moves elements without shrinking that source vector. The inherited
PeptideIdentificationList move constructor instead transfers the underlying
vector; its source tests inspect the moved-from vector's zero capacity. Rust
ownership transfer makes the moved-from binding unavailable, or `mem::take`
replaces it with an empty vector. C++ moved-from-object inspection is an explicit
representation boundary, not an additional identification operation.

Existing evidence includes [identification tests](../tests/identification.rs),
[filter tests](../tests/id_filter.rs), [protein-inference tests](../tests/protein_inference.rs)
and the [idXML workflow](../tests/identification_pipeline.rs): these construct,
clone, compare, iterate, reorder and filter vectors while retaining the actual
identification fields. The source list class test adds only standard container
construction/copy/move/assignment, size, clear/resize and equality assertions.

## Every ExposedVector operation

The template stores exactly one `std::vector<T>` and delegates its operations to
that member. `LessThanComparable` and the interface macro are compile-time
plumbing; Rust uses ordinary trait bounds and the vector/slice types.

| Source declarations | Native equivalent |
| --- | --- |
| `VecMember`, `value_type` | `Vec<T>`, `T` |
| `iterator`, `const_iterator` | Mutable/immutable slice iterators, or an owning `IntoIter<T>` |
| `reverse_iterator`, `const_reverse_iterator` | `.iter_mut().rev()` / `.iter().rev()` |
| `size_type`, `difference_type` | `usize`, `isize` |
| `pointer`, `reference`, `const_reference` | Typed pointers where required; normally `&mut T` / `&T` |
| Default constructor | `Vec::new()` / `Default` |
| Count constructor | `resize_with(n, T::default)` on an empty vector |
| Count-and-value constructor | `vec![value; n]`, requiring `Clone` |
| Iterator-range constructor | Collect the owned iterator, or clone borrowed elements |
| Copy construction/assignment | `clone` / `clone_from` |
| Move construction/assignment | Rust ownership transfer; `mem::take` or `mem::replace` for an accessible old binding |
| Destruction | Automatic drop |
| Mutable/const `begin`, `end`, `cbegin`, `cend` | `.iter_mut()` / `.iter()` and exhaustion; slices/index ranges when random positions are needed |
| Mutable/const `rbegin`, `rend`, `crbegin`, `crend` | Reversed double-ended iterators |
| `size`, `empty` | `len`, `is_empty` |
| `resize` | `resize_with(n, T::default)` or `resize(n, value)` |
| `reserve` | `try_reserve` / `reserve` with the additional count needed to reach the requested total |
| `capacity`, `shrink_to_fit` | `capacity`, `shrink_to_fit` |
| `max_size` | Native allocation ceiling: at most `isize::MAX` bytes for nonzero-sized `T`; actual success checked with `try_reserve` |
| Mutable/const `operator[]`, `at` | Indexing or `get`/`get_mut` |
| Mutable/const `front`, `back` | `first`/`first_mut`, `last`/`last_mut` |
| Copy/move `push_back` | `push(value.clone())` / `push(value)` |
| `emplace_back` | Construct `T`, `push`, then `last_mut` if the returned reference is needed |
| `pop_back` | `pop` |
| Single/range `erase` | `remove` / `drain`; retain the index or obtain the following element with `get` |
| Single copy/move `insert` | `insert(index, value.clone())` / `insert(index, value)` |
| Iterator-range/count/initializer-list `insert` | `splice(index..index, values)` with owned/cloned values or repeated clones |
| `emplace` | Construct `T`, insert at the index and borrow that element if needed |
| `clear` | `clear` |
| `swap` | `std::mem::swap` |
| Iterator/count/initializer-list `assign` | Assign a collected vector, or clear and extend from the corresponding values |
| Mutable/const `getData` | The vector itself, or mutable/immutable slices |
| `==`, `!=` | Elementwise vector equality when `T: PartialEq` |
| Conditional `<`, `<=`, `>`, `>=` | Lexicographic vector comparison when `T: PartialOrd`; `Ord` for total ordering |

Source vector `reserve` takes a target total capacity; Rust takes an additional
count relative to length. For target `n`, use `n.saturating_sub(vec.len())`.
Neither C++ nor Rust guarantees identical capacity growth or that shrink-to-fit
releases all spare space. C++ allocator-specific `max_size` values are not
reproduced: Rust's byte limit, checked allocation and zero-sized-type rules are
used instead. C++ has no zero-sized object analogue to Rust's special Vec case.
Out-of-range C++ indexing and empty front/back/pop are not reproduced as undefined
behavior; Rust offers checked `Option` access or indexing panics. Rust borrowing
also prevents using an iterator/reference across an invalidating mutation.

The native library already uses these standard operations in
[kernel.rs](../src/kernel.rs), [ID filtering](../src/analysis/id_filter.rs),
[protein inference](../src/analysis/protein_inference.rs) and their tests. This
classification does not certify a full generic C++ subclass/ABI emulation.

## Basic list aliases

`TypeAliases.h` has exactly three declarations:

| Source | Native type |
| --- | --- |
| `IntList = std::vector<Int>` | `Vec<i32>` on the SDK's supported 32-bit-C-int platforms |
| `DoubleList = std::vector<double>` | `Vec<f64>` |
| `StringList = std::vector<std::string>` | `Vec<String>` for text; `Vec<Vec<u8>>` when arbitrary byte strings must be retained |

A C++ string is not necessarily UTF-8; conversion of arbitrary bytes to `String`
is not implicitly lossless. These are container aliases, not a claim that all
StringUtils/ListUtils behavior exists. The finite-number checks in typed
`MetaValue` are additional domain rules, not restrictions of `Vec<f64>` itself.
In particular, source IntList is not the native metadata `Vec<i64>` alternative.

`Vec<i32>` appears in precursor possible-charge states, and `Vec<f64>` /
`Vec<String>` in the typed metadata alternatives. [Metadata tests](../tests/metadata.rs)
and [precursor workflows](../tests/precursor_workflow.rs) exercise these stored
lists, including order and duplicates.

## Types.h remains partial

All primitive aliases and ASCII declarations have straightforward native
representations, but this header also contains behavior that cannot be marked
complete from the existing port:

| Source declaration | Mapping / boundary |
| --- | --- |
| `Int32`, `UInt32`, `Int64`, `UInt64` | `i32`, `u32`, `i64`, `u64` |
| `Int`, `UInt` | `i32`, `u32` on supported SDK targets; `std::ffi::c_int` / `c_uint` when C ABI width is specifically required |
| `Byte`, `UID` | `u8`, `u64`; this does not implement UniqueIdGenerator/UniqueIdInterface |
| `Size`, `SignedSize` | `usize`, `isize` |
| `Time = time_t` | No portable time_t alias or whole-header time ABI bridge; calendar CompletionTime is a different domain type |
| `ASCII__BACKSPACE`, `ASCII__BELL` | `b'\x08'`, `b'\x07'` |
| `ASCII__CARRIAGE_RETURN`, `ASCII__HORIZONTAL_TAB`, `ASCII__TAB` | `b'\r'`, `b'\t'`, `b'\t'` |
| `ASCII__NEWLINE`, `ASCII__RETURN`, `ASCII__SPACE`, `ASCII__VERTICAL_TAB` | `b'\n'`, `b'\n'`, `b' '`, `b'\x0b'` |
| `ASCII__COLON`, `ASCII__COMMA`, `ASCII__EXCLAMATION_MARK` | `b':'`, `b','`, `b'!'` |
| `ASCII__POINT`, `ASCII__QUESTION_MARK`, `ASCII__SEMICOLON` | `b'.'`, `b'?'`, `b';'` |
| `writtenDigits<float>`, `writtenDigits<double>` | Rust `f32::DIGITS == 6`, `f64::DIGITS == 15` represent these values |
| `writtenDigits<int>`, `writtenDigits<unsigned int>` | Source decimal precision is 9 on supported 32-bit-int targets |
| `writtenDigits<long int>`, `writtenDigits<unsigned long int>` | Source deliberately returns the **int/unsigned-int** precision, not long's precision; the whole source helper is absent |
| `writtenDigits<long double>` | Non-Windows uses the platform long-double precision; Windows forces double precision. Stable Rust 1.85 has no matching cross-platform long-double type |
| Generic `writtenDigits<T>` declaration/default argument and fallback | Source accepts an unused/default-constructed value and falls back to 6; no equivalent generic port helper is supplied |

Existing Rust shortest-roundtrip formatting must not be labeled identical to the
source's significant-decimal-digit policy. The primitive types are widely used
and tested, but that evidence does not complete the extra helper/platform API.

## ConstRefVector is not just Vec of references

The complete header was inspected, including both custom iterator classes,
constructors/assignment, vector operations, comparisons, swaps and all four sort
paths (ascending/descending intensity, position, custom comparator).
`Vec<&T>` with slice sorting captures its primary borrowed-view use case, but does
not reproduce these additional source rules:

- Equality checks base-container pointer identity, length, dynamic element type,
  then element inequality. Two equal-valued views from different source
  containers compare unequal.
- `<` and `>` compare **only lengths**; `<=` / `>=` combine length ordering with
  the custom equality. Equal-sized distinct views can be incomparable.
- Capacity is `max(size, capacity_)` with a separately maintained cached value.
  Swapping transfers the pointer vector but leaves cached capacity and source
  container identity in place. This can change later equality/capacity results.
- Count construction and growing resize can create null reference slots; iterator
  positions narrow to unsigned int. Const/mutable iterators expose custom
  position arithmetic, equality and swap semantics.

Those are concrete unimplemented behaviors, not safe grounds to classify the
whole header as a standard-library equivalent. No unsafe/null-reference wrapper
was added. The source class test is not evidence that these behaviors are present
in the native library; its overall header remains partial.
