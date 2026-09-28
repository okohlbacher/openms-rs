# MetaInfoDescription and the data arrays

`metadata::meta_info_description` covers `METADATA/MetaInfoDescription.h` and
`METADATA/MetaInfoDescription.cpp`; together with `kernel::DataArray<T>` it covers
`METADATA/DataArrays.h` (`DataArrays.cpp` is an empty translation unit), all at SDK `bc9cc12`.
Earlier work on the arrays' descriptions, their preservation through processing and their XML
transport is in [`MOBILOGRAM_SUPPORT.md`](MOBILOGRAM_SUPPORT.md).

Tests: [`tests/data_arrays.rs`](../tests/data_arrays.rs), every `START_SECTION` of
`DataArrays_test.cpp` and `MetaInfoDescription_test.cpp` with unchanged literals, and tests
derived from `MetaInfoDescription.cpp:29-78` for each step of the ordering and from
`DataArrays.h` for the tie-break on the values.

## API mapping

| C++ | Rust | Notes |
|---|---|---|
| `class MetaInfoDescription : public MetaInfoInterface` | `MetaInfoDescription { name, metadata, data_processing }` | Public fields: `getName`/`setName`, the `MetaInfoInterface` base, `getDataProcessing` (both) / `setDataProcessing`. |
| constructors, assignments, destructor | `Default`, `Clone`, moves, `Drop` | |
| `operator==`, `operator!=` | `PartialEq` | Data processing compared by pointee (`Helpers::cmpPtrSafe`). |
| `operator<`, `<=`, `>`, `>=` | `source_lt`, `source_le`, `source_gt`, `source_ge` `(&rhs, &MetaInfoRegistry) -> Result<bool>` | Not `PartialOrd`: the order is not consistent with `==` (see below), and it needs the registry. |
| `DataArrays::FloatDataArray`, `IntegerDataArray`, `StringDataArray` | `DataArray<f32>`, `DataArray<i32>`, `DataArray<String>` | The `MetaInfoDescription` base is the array's `name`, `metadata` and `data_processing`; the `std::vector` base is `data`. |
| their `operator==`/`!=` | derived `PartialEq` | Description and values; `NaN != NaN` as in the source's class test. |
| their `operator<`, `<=`, `>`, `>=` | `DataArray::source_lt`, `source_le`, `source_gt`, `source_ge` | Descriptions first (when unequal), then `std::lexicographical_compare` on the values, reproduced element by element so that NaN is neither less nor greater. |
| — | `DataArray::description_view()`, `DataArray::description()`, `DescriptionView`, `description_lt`, `data_value_text` | The borrowed and owned description of an array, the ordering over views, and `DataValue::toString()` as the ordering uses it. |

## The ordering

`MetaInfoDescription::operator<` (`MetaInfoDescription.cpp:29-78`) is reproduced step by step:
fewer data-processing entries first; then, **only when the meta values differ**, a description
with values before one without, the name, the keys as registry indices, and each value's
`DataValue::toString()` text; otherwise the processing-action counts entry by entry; otherwise
`false`. `<=` is `< || ==`, `>` is `!(<=)`, `>=` is `!(<)`, exactly as the source composes them.

Consequences kept from the source:

- Two descriptions that differ only in their names are unequal, yet neither is less, and each is
  greater than the other. `DataArrays_test.cpp`'s mixed-metadata section passes on exactly this:
  `A < B` is false and `A > B` is true. **Source defect candidate**: not a strict weak ordering,
  so a `std::set` or `std::sort` over such arrays is undefined.
- Values compare by their text, so an integer 10 sorts before 9.
- The keys are ordered by registry index (registration order), not by name, which is why the
  comparison takes the caller's `MetaInfoRegistry`; a name missing from it is an
  `Error::InvalidValue` when the order needs it.

A null data-processing pointer (the source's class test resizes the vector to one null entry)
cannot exist in an `Arc`; the two null branches of the ordering have no counterpart, and the
class-test sections use a default `DataProcessing` record instead.

## Evidence

Tier 3: both class tests transcribed with unchanged literals; the ordering steps are derived
from `MetaInfoDescription.cpp` and cite its lines. No C++ was executed.
