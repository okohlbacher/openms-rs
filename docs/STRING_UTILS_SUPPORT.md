# StringUtils support

`DATASTRUCTURES/StringUtils.h` and its numeric core `CONCEPT/Detail/NumericFormatting.h`
at core `bc9cc12`, ported by lane B8a (Phase 3 wave 4). The member-by-member API
table is the rustdoc of `openms::data_structures::string_utils`; this document
records where the port lives, what it deliberately does differently and the
evidence.

## Where it lives

| Rust | Source | Role |
|---|---|---|
| `src/concept/numeric_formatting.rs` | `CONCEPT/Detail/NumericFormatting.h` | `appendNumeric` for `f64` and `f32`, both precisions; C `%.*g` stream text (`ostream_g`) |
| `src/data_structures/string_utils.rs` | `DATASTRUCTURES/StringUtils.h`, `DATASTRUCTURES/StringUtils.cpp` | every public member: parsing, `toStr`/`appendToStr`, `number`, `numberLength`, `random`, predicates, accessors, mutators, quoting, splitting, `concatenate`, the `StringConversions` forwarding namespace, the global `std::string` + number operators |
| `src/metadata/value_text.rs` | `toStr(const DataValue&, bool)` | `MetaValue::to_text`, beside the value type |
| `src/param/value.rs` | `toStr(const ParamValue&, bool)` | `ParamValue::to_text` |
| `src/format/file_info/text_format.rs` | FileInfo's use of the same routines | re-exports the moved items, so its paths stay valid |

Module placement follows the module graph: `concept` and `param` are sinks,
`data_structures` depends on `param` only and `metadata` on `data_structures`,
so the numeric core sits in `concept` (where `NumericFormatting.h` lives in the
source) and every module can reach it; `python3 tools/check_module_cycles.py`
passes.

## Native differences

1. **Text model.** A C++ `std::string` is bytes, a Rust `String` is UTF-8.
   Every function keeps the source's byte semantics on ASCII text, which is
   all the source's callers feed it. A `char` argument that is not ASCII
   (a quote, fill character or splitter, which is one byte of the text in the
   source) is refused, and so is a byte count that would cut a multi-byte
   character (`prefix`, `suffix`, `substr`, `chop`) and `reverse` of non-ASCII
   text, because the source's result there is not UTF-8.
2. **Ceilings.** Functions whose result can outgrow their input (`fill_left`,
   `fill_right`, `substitute`, `quote`, `concatenate`, `random`) refuse a
   result above `MAX_TEXT_BYTES` before allocating and reserve fallibly;
   `split` refuses more than `MAX_SPLIT_PARTS` parts. Every refusal leaves the
   input unchanged. The source has no ceiling.
3. **`number`.** `number` refuses what the source's 64-byte `snprintf`
   buffer would silently cut (CPP-253) and a digit count that becomes a
   negative `printf` precision (CPP-254); `number_truncated` reproduces the
   source bytes.
4. **`numberLength`.** A non-finite value is refused (the source loops forever
   on an infinity and converts a NaN to `int`, which is undefined), and so is a
   truncated magnitude that does not fit an `int` (undefined in the source; x86
   gives `INT_MIN`). The oracle's 19 such records are classified, not skipped.
5. **`random`.** `random_with` takes the uniform source from the caller;
   `random` seeds a Mersenne Twister (`rand_mt`) from the clock and a call counter. The source uses `rand()`, whose
   sequence is the C library's, so no expected text exists for it; the tests
   check length and alphabet.
6. **`long double`.** `toStr(long double, bool)`, `appendToStr(long double)`,
   `appendToStrLowP(long double)` and the `long double` operators have no Rust
   type to port to. The `f64` functions cover every value a Rust caller holds.
7. **`toStr(DataValue)` with a non-finite value.** The public `MetaValue`
   constructors refuse NaN and the infinities (`MetaValue::validate`), so the
   nine oracle records of non-finite `DataValue`s cannot be built through the
   public API; the same text is covered by the `toStr(double)` records, and the
   crate-private `MetaValue::source_float` path the XML readers use writes it.

## Source behaviour kept

- Parsing is `std::from_chars` after the source's own blank and `+` handling:
  `+-1` is `-1`, `++1` and `-+1` fail, `nan(...)` payload forms are read, an
  underflow to zero is a `ConversionError` while a subnormal result parses (the
  Release build's libstdc++ `from_chars`).
- Error texts are the source's `ConversionError` messages, verbatim.
- `toStr(double)`/`toStr(float)` are `appendNumeric`: fixed notation in
  `[1e-2, 1e4)` with trailing zeros trimmed to one fraction digit, shortest
  round-trip scientific text with `std::to_chars`' ties to even otherwise, an
  integer mantissa keeping its `.0` (`1.0e06`), `NaN` for either sign.

## Evidence

- **Tier 1:** `../oracle/B8a-stringutils` (driver `driver.cpp`) links the Linux
  x86_64 Release install `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`
  on ibminode06 and records 5,068 calls over the whole header; two runs were
  byte-identical. `tests/data/string_utils_oracle.tsv` is its `records.tsv`,
  and `tests/string_utils.rs::oracle_records_match` recomputes every record and
  pins how each is classified (compared, non-UTF-8 input, non-UTF-8 result,
  undefined in the source, non-finite `DataValue`). The earlier
  `../oracle/integ-F-stringutils` (246 parsing rows, same install and host)
  uses inputs that are all among these records.
- **Tier 3:** `tests/string_utils.rs` ports the literals of
  `StringUtils_test.cpp` (`class_test_*`) and `String_test.cpp`
  (`string_test_*`) at core `bc9cc12`, except the `std::string` constructor
  sections and the `long double` sections, which have no port counterpart.
- The FileInfo scope (`tests/file_info_text_format.rs`, tier 1 against
  `../oracle/file-info-text-format` and `../oracle/text-format`) is unchanged
  and still passes through the re-exports.
