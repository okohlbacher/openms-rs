# Peak-file option values

`format::PeakFileOptions` provides the complete option state exposed by
OpenMS4-core `54a232fe2cae9c590d5c997fa49d20e7769860fb`'s `PeakFileOptions`.
The type, `NumpressCompression`, `NumpressConfig` and constants are available in
`format::peak_options` without optional features or new dependencies.

**This is a value API. Existing format readers and writers do not consume this
new type yet.** Setting a filter, precision, Numpress mode or compatibility flag
neither executes that operation nor silently alters existing `mzml::ReadOptions`
or `WriteOptions`. Resource-limit options remain separate. Unsupported executing
behavior must be rejected explicitly when adapters are wired in a later batch.

## Scalar state and defaults

Public Rust fields replace the source's mechanical getter/setter pairs. `Default`
and `new` have identical state; `Clone` makes an independent copy. `PartialEq` is
a native convenience using ordinary floating-point equality, so NaN-bearing
values need not compare equal to themselves.

| Native field | Source getter/setter stem | Default |
| --- | --- | --- |
| `metadata_only` | `MetadataOnly` | false |
| `force_mq_compatibility` | `ForceMQCompatability` | false |
| `force_tpp_compatibility` | `ForceTPPCompatability` | false |
| `write_supplemental_data` | `WriteSupplementalData` | true |
| `mz_32_bit` | `Mz32Bit` | false |
| `intensity_32_bit` | `Intensity32Bit` | true |
| `zlib_compression` | `Compression` | false |
| `always_append_data` | `AlwaysAppendData` | false |
| `skip_xml_checks` | `SkipXMLChecks` | false |
| `sort_spectra_by_mz` | `SortSpectraByMZ` | true |
| `sort_chromatograms_by_rt` | `SortChromatogramsByRT` | true |
| `fill_data` | `FillData` | true |
| `write_index` | `WriteIndex` | true |
| `max_data_pool_size` | `MaxDataPoolSize` | 100 |
| `precursor_mz_selected_ion` | `PrecursorMZSelectedIon` | true |
| `skip_chromatograms` | `SkipChromatograms` | false |

The MQ switch is source mzXML-specific; TPP and precursor selection are
mzML-specific; supplemental writing is mzData-specific. The m/z precision flag
also controls chromatogram retention-time precision. Pool size is a processing
batch size, not a security/resource limit: the value API accepts zero and the
complete native `usize` range, without allocating a pool.

## Ranges and MS levels

Methods `set_rt_range`, `rt_range`, `has_rt_range` and the corresponding
`mz_range`, `intensity_range`, `precursor_mz_range` methods preserve all four
source states. Endpoints use existing `kernel::NumericRange` as a small value
container. Source filtering uses **half-open** DRange membership rather than the
inclusive interpretation used for kernel extrema. This module performs no range
membership tests or scientific filtering.

`EMPTY_PEAK_FILE_RANGE` is exactly `{ min: f64::MAX, max: f64::MIN }`, using finite
extrema. The DRange header's statement about zero initialization is stale; its
DIntervalBase constructor uses these extrema. All range getters initially return
that sentinel and all four `has_*` predicates initially return false.

Setting RT to the exact sentinel disables RT filtering. Equal endpoints,
infinities, NaNs and all other endpoint values enable it. Setting m/z, intensity
or precursor-m/z **always enables that option**, even with the empty sentinel;
these source setters do not provide a clearing operation. A fresh/default options
value resets the full state. `has_filters` is precisely:

```text
has_rt_range OR has_ms_levels OR has_precursor_mz_range
```

It intentionally excludes peak m/z/intensity filtering, metadata-only loading,
skipping chromatograms and every writing choice. It is not a general predicate
for whether an adapter's output can change.

Setters copy supplied endpoints without validation or normalization, as the
source option setters copy an already-constructed DRange. The source DRange
constructor taking two points sorts inverted coordinates; constructing a raw
NumericRange does not. Native callers emulating that constructor must order their
endpoints first. The explicit empty sentinel must stay inverted. Storing raw
inverted non-sentinel bounds is an additional native representable state; it is
not evidence for equivalent source point-constructor behavior.

`set_ms_levels(&[i32])`, `add_ms_level`, `clear_ms_levels`, `ms_levels`,
`has_ms_levels` and `contains_ms_level` mirror the source vector operations.
They preserve signed values, zero, duplicates and caller order. Empty membership
returns false; an adapter must test `has_ms_levels` before treating membership
as a filter. `set_ms_levels` copies its input and fails atomically when the native
`MAX_PEAK_FILE_MS_LEVELS` bound of 1,000,000 elements is exceeded. Adding beyond
that bound also leaves state intact. Membership remains a bounded linear search;
no hidden set changes ordering or duplicates.

## Required Numpress configuration values

`NumpressCompression::{None, Linear, Pic, Slof}` has the source numeric identities
0 through 3. `ALL`, `name()` and `FromStr` provide the exact case-sensitive names
`none`, `linear`, `pic`, `slof`. Whitespace is not trimmed. The source enum's count
sentinel is not a usable mode. `NumpressConfig::set_compression` changes only the
mode and leaves the configuration unchanged on an unknown string.

| NumpressConfig field | Source field | Default |
| --- | --- | --- |
| `fixed_point` | `numpressFixedPoint` | 0.0 |
| `error_tolerance` | `numpressErrorTolerance` | 0.0001 |
| `compression` | `np_compression` | None |
| `estimate_fixed_point` | `estimate_fixed_point` | true |
| `linear_fp_mass_acc` | `linear_fp_mass_acc` | -1.0 |

These are stored values, so negative and nonfinite numeric values remain
representable exactly as in the source configuration. No numerical algorithm
runs here. Future encoders must enforce algorithm-specific requirements when
consuming them. The source error-tolerance zero disables its encode/check step;
linear mass-accuracy -1 means no requested accuracy. Merely storing those values
here does not claim an encoder or validator implementation.

Independent configurations are exposed through
`set_numpress_configuration_mass_time`, `set_numpress_configuration_intensity`,
`set_numpress_configuration_float_data_array` and matching getters without `set_`.
All getters return copies. The mass/time setter retains PIC/SLOF configurations
and returns `Some(NUMPRESS_MASS_TIME_WARNING)` containing the exact source warning
text; other modes return none. This replaces unconditional source stderr output
with explicit caller-controlled routing. Intensity and float-array setters return
no warning. None of the MSNumpressCoder encode/decode operations is implemented
or implied by this configuration support.

## Reference coverage and remaining integration

[Nine tests](../tests/peak_file_options.rs) cover every header default, the source
class-test ranges `[2,4]`, `[3,5]`, `[400,1200]`, levels `[1,3,5]`, pool size 250,
and the three literal `hasFilters` examples. Separate source-derived tests check
sentinel asymmetry, zero-width ranges, signed/order/duplicate state, complete
Numpress copies and exact warning text. NaN-bit preservation and allocation
bounds are native checks. Numpress configuration numeric defaults come from its
header; its class test only checks construction, and is not presented as a
numerical output oracle. No C++ was built or executed.

[Provenance](../tests/data/peak_file_options_provenance.json) records nine source
hashes, literal values and source line anchors, with source evidence distinct
from derived/native checks. The existing [mzML support](MZML_SUPPORT.md) and
[DTA2D/MS2 support](TEXT_PEAK_LIST_SUPPORT.md) remain authoritative about current
adapter behavior; this state API does not expand those format promises.

The next adapter batch should pass scientific options alongside each adapter's
existing resource limits. It must apply RT/MS-level/precursor filters at the
correct source record stage, m/z/intensity filters with aligned auxiliary-array
selection, and sorting with annotation alignment. Metadata-only, fill-data,
consumer append, and skip-chromatogram behavior need separate streaming/counting
contracts. Precursor selection must retain the distinction between selected-ion
m/z and isolation target. Encoding options require actual selected precision,
zlib, index and compatibility implementations; unavailable Numpress and format
choices must error before output. `skip_xml_checks` must never disable native
memory/work limits or unsafe-input guards. Repeated MS-level membership work must be charged across all records. Batch size
must not become an unbounded allocation request. No such wiring is part of this staged change.
