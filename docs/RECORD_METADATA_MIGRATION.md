# Spectrum and chromatogram metadata migration

`MSSpectrum.metadata` and `MSChromatogram.metadata` are now `MetaInfo`, the same typed map already used by experiment settings, features and acquisition descriptions. No duplicate legacy map exists. The source C++ records likewise inherit typed MetaInfo storage.

Old string values retain their type. Convert an existing String map with `metadata::meta_from_strings`, which copies strings without interpreting their contents, or move each pair with `value.into()`:

```rust
use openms::{MSSpectrum, metadata::MetaValue};
let mut spectrum = MSSpectrum::default();
spectrum.metadata.insert("sample".into(), "42".into());
assert_eq!(spectrum.metadata["sample"].as_str()?, "42");
spectrum.metadata.insert("count".into(), 42_i64.into());
spectrum.metadata.insert("time".into(), MetaValue::try_from(2.5)?);
spectrum.metadata.insert("grid".into(), MetaValue::try_from(vec![1., 3.])?);
# Ok::<(), openms::Error>(())
```

Replace direct string indexing/borrowing with `as_str()?` when a key is known to contain String; use `as_i64`, `as_f64`, `as_string_list`, `as_integer_list` or `as_float_list` for typed values. Float constructors reject nonfinite values. MetaValue also supports Empty and units. Copying a record owns its metadata independently; references behind processing-history Arc handles remain shared as before.

Existing caller-created strings are never automatically parsed, even for numeric-looking or list-looking values. New SpectrumAnnotator output now uses Float `fragment_mass_tolerance` and Integer `fragment_mass_tolerance_ppm` (0/1), matching the source. New alignment output uses Float `original_RT` and FloatList `original_rt`; a preexisting original-value key of any type is retained. Feature and peptide metadata already used MetaInfo and retain their existing behavior.

The [mzML typed transport group](MZML_TYPED_TRANSPORT_SUPPORT.md) writes scalar String/Integer/Float values and their units. Generic lists and Empty values are explicitly rejected before output because source generic XML writing stringifies them. Three spectrum noise keys have a dedicated FloatList binary route with independent lengths and f64 precision. mzML primary array selectors remain ordinary unit-free Strings in the same record map.

MGF still preserves its original unit-free string fields. Numeric, list, Empty or unit-bearing MGF values now fail preflight instead of being stringified. DTA rejects any nonempty spectrum metadata because that format cannot preserve it. Existing MS2/DTA2D restrictions and XML transport guards remain. Converting a value to text before insertion is an explicit application choice with the same limits as before.

Copying transformations now include complete typed record metadata in their existing shared acquisition-copy allowance: 50 million work units and a conservative 256 MiB storage estimate. Nested and batched clones use the same allowance. Theoretical append and EMG copy paths likewise preflight the new owned payload. Metadata insertion/replacement accounts for sparse tree slots, key comparison and replaced list/string payload before mutation. In-place operations which retain or ignore metadata do not gain extra copy/update accounting beyond ordinary record validation. The bounds are logical estimates rather than measured physical allocator usage.

The existing scientific number formats, peak order, aligned arrays and native metadata ownership remain unchanged outside these stated producer/transport differences. The update is an API change for typed access, not a conversion of historic stored text into inferred values.
