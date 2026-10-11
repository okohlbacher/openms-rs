# CONCEPT basics: CommonEnums, Types, HashUtils

`CONCEPT/CommonEnums.h`, `CONCEPT/Types.h` and `CONCEPT/HashUtils.h` at core
`bc9cc12` are ported in `src/concept/common_enums.rs`, `src/concept/types.rs`
and `src/concept/hash_utils.rs`. None has a `.cpp` with behaviour beyond the
header and none has a class test; the test is `tests/concept_basics.rs`.
`concept` depends on no other module of the crate, so every module can use
these items.

## API mapping

### `CommonEnums.h`

| Source | Rust |
|---|---|
| `enum class DIM_UNIT` (`RT`, `MZ`, `INT`, `IM_MS`, `IM_VSSC`, `FAIMS_CV`, `IM_CCS`) | `concept::common_enums::DimUnit`, source discriminants 0 to 6 |
| `DIM_UNIT::SIZE_OF_DIM_UNITS` | `DimUnit::SIZE_OF_DIM_UNITS`, the length of `DimUnit::ALL` |
| `DIM_NAMES[]`, `DIM_NAMES_SHORT[]` | `DIM_NAMES`, `DIM_NAMES_SHORT`; `DimUnit::name`, `DimUnit::short_name` |
| `enum class MZ_UNITS` (`DA`, `PPM`), `SIZE_OF_MZ_UNITS` | `MzUnits`, `MzUnits::SIZE_OF_MZ_UNITS` |
| `MZ_UNIT_NAMES[]` | `MZ_UNIT_NAMES`, `MzUnits::name`, `FromStr` |
| `enum class HydrophobicityScaleMethod` | `chemistry::hydrophobicity::HydrophobicityScale` (all seven, discriminants 0 to 6) |
| `(int)X` casts | `X as usize`, `index()`, `from_index()` |

### `Types.h`

| Source | Rust |
|---|---|
| `Int32`, `UInt32`, `Int64`, `UInt64`, `Time`, `UInt`, `Int`, `Byte`, `UID`, `Size`, `SignedSize` | the type aliases of `concept::types` (`Time` is `i64`, `Uid` spells `UID`) |
| `enum ASCII` and its 15 `ASCII__*` enumerators | the constants of `concept::types::ascii` |
| `writtenDigits<T>()`, the generic fallback and the `float`, `double`, `int`, `unsigned int`, `long int`, `unsigned long int`, `long double` specialisations | `written_digits::<T>()` over `WrittenDigits`; `LONG_DOUBLE_WRITTEN_DIGITS` |

### `HashUtils.h`

| Source | Rust |
|---|---|
| `fnv1a_hash_bytes`, `fnv1a_hash_string` | `fnv1a_hash_bytes(&[u8])`, `fnv1a_hash_string(&str)` (`const fn`) |
| `hash_combine` | `hash_combine(&mut u64, u64)` |
| `hash_int<T>` | `hash_int` over `HashInt` (`bool`, `char`, every integer width) |
| `hash_char` | `hash_char(u8)` |
| `hash_float<T>` | `hash_float` over `HashFloat` (`f32`, `f64`), with the signed-zero normalisation |
| the FNV constants and the golden ratio | `FNV_OFFSET_BASIS`, `FNV_PRIME`, `GOLDEN_RATIO` |

## Preserved source conventions

- Name tables, discriminants and the order of every enum are the header's.
- `writtenDigits` returns `numeric_limits<T>::digits10`, including the
  source's `long int`/`unsigned long int` specialisations that return the
  `int` values.
- `hash_float` hashes `+0.0` for `-0.0`, as the header's `if (value == T(0))`
  does; NaN payloads hash as their bytes.

## Native differences

- The `SIZE_OF_*` enumerators are associated constants, not values of the
  enum, so no unit value indexes past the name tables.
- `HydrophobicityScaleMethod` lives with its table in `chemistry`, because
  `concept` may not depend on `chemistry`.
- Hashes are 64-bit and use little-endian bytes on every target; the source's
  `size_t` is 32-bit on 32-bit hosts and its bytes are native-endian. On the
  64-bit little-endian hosts the source is built for the values are equal.
- `long double` has no Rust type; `LONG_DOUBLE_WRITTEN_DIGITS` records the
  x86_64 Linux value (18).

## Checked boundaries and evidence

Tier 3: the header literals (name tables, discriminants, ASCII codes,
`digits10`). Tier 4: the published FNV-1a 64-bit test vectors and FNV-1a and
Boost `hash_combine` recomputed in Python over explicit little-endian bytes.
Tier 1 indirectly: every `source_hash64` in the crate built on these
functions (`EmpiricalFormula`, `Element`, `IsotopeDistribution`,
`Ribonucleotide`) matches the executed Release build
(`tests/chemistry_members_oracle.rs`). All functions are constant work or
linear in their input bytes.
