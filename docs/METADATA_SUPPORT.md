# Typed metadata and acquisition settings

`openms::metadata` provides native owned values, controlled-vocabulary terms and acquisition records derived from OpenMS4-core revision `7c029e8`. These are functional data models with validation and source-derived merge behavior. They do not load an ontology or claim coverage of the complete OpenMS metadata/identification object graph.

## Values and string interoperability

`MetaInfo` is `BTreeMap<String, MetaValue>`. Keys are owned strings; there is no process-wide numeric name registry. Missing keys use ordinary `Option` lookup. A present empty value is distinct from a missing key, an empty string and an empty list.

`MetaValueData` represents the seven DataValue alternatives: `Empty`, `String`, `Integer(i64)`, `Float(f64)`, `StringList`, `IntegerList` and `FloatList`. `MetaValue` owns a private, validated alternative and an optional `Unit`. Integers/strings and their lists support `From`; floating scalars/lists use `TryFrom` to reject NaN and infinities. `MetaValue::new` provides the enum-based constructor. Unsigned uint64 conversion is checked against int64 range. No mutable reference to stored floating data is exposed.

`as_str`, `as_i64` and list accessors require the matching type. `as_f64` accepts float or integer, with the usual loss of precision for integers beyond exactly representable float64 values. Strings are not implicitly parsed as numbers. `to_bool` follows the source convention: only the exact strings `true` and `false` are accepted; integer 0/1 and alternative spellings are errors.

`Display` provides lenient stringification for logging and source-style spectrum references. Empty values format as an empty string; lists use `[first, second]`. Scalar floats use Rust's round-trippable representation, not C++'s historical 17-digit/scientific formatting thresholds. List formatting is human-readable and does not escape commas or brackets inside strings, so it is not a lossless serialization format. Units are not appended.

`Unit` stores a validated nonempty accession, name and vocabulary reference. `Unit::from_ontology` maps source UO/MS numeric IDs into explicit accessions; arbitrary vocabularies use `Unit::new`. Missing units are `None`, replacing DataValue's numeric -1 sentinel and `OTHER` tag. Unit fields are immutable after construction. No unit conversion, vocabulary membership validation or canonical-name lookup is inferred.

Equality compares exact type, exact values and unit identity. This deliberately differs from source DataValue scalar-double equality (`abs(a-b) < 1e-6`), which is not transitive and differs from the source's exact list equality. No approximate hash or cross-type ordering is exposed. Finite -0 and +0 compare equal as ordinary Rust floats do.

Spectra, chromatograms and experiments still use string metadata through explicit bridges. Feature, consensus, map and column records now use `MetaInfo` directly:

| Function | Behavior |
| --- | --- |
| `meta_from_strings` | Copies every value as a string; does not guess numeric types or units |
| `meta_to_strings` | Accepts only strings without units; errors rather than losing type/unit information |
| `meta_to_strings_lossy` | Explicitly discards types and units using Display |
| `validate_meta` | Checks every typed value; available to nested identification and acquisition models |
| `merge_meta` | Supports overwrite, keep-existing and reject-conflicts policies; validates before mutation |

Overwrite matches `MetaInfoInterface::addMetaValues` / `MetaInfo::operator+=` at the pinned revision. Reject-conflicts permits equal shared values and rejects unequal shared values transactionally.

## Controlled vocabulary terms

`CVTerm` stores accession, name, vocabulary reference and a `MetaValue`; its value's optional unit is also the term's unit. This unifies the C++ CVTerm and DataValue unit mechanisms instead of allowing two contradictory unit fields. An empty value with a unit is supported: `has_value` remains false while `has_unit` is true. A present empty string counts as a value, matching the source distinction.

`CVTermList` owns an accession-indexed map of ordered term vectors and ordinary `metadata`. Its map is private and immutable through readers, preventing disagreement between a map key and a contained accession.

- `add` and `add_terms` append, retaining duplicates. The latter maps the additive behavior of C++ `setCVTerms`, whose name can suggest replacement.
- `replace` replaces one accession with one term; `replace_accession` replaces its whole vector; `replace_all` replaces the complete CV map and preserves ordinary metadata.
- `consume` appends all source term vectors without deduplication and leaves ordinary metadata unchanged, matching C++ `consumeCVTerms`.
- An explicitly stored empty vector still makes `contains(accession)` true. `is_empty` tests CV keys only, independently of ordinary metadata.

Batch mutations validate all terms before modification. Empty/whitespace accessions and mismatched accession-map keys return errors; these correct unchecked source cases. Terms otherwise remain opaque: unknown accessions are preserved, and no CV spelling validation is claimed.

## Acquisition model coverage

| Native type | Source fields and behavior |
| --- | --- |
| `Precursor` / `PrecursorInfo` | Kernel precursor owns m/z, charge/intensity, activation-method set, energy, isolation offsets/optional distinct target, drift value/unit/offsets, possible charge states, CV terms and parent spectrum reference; the compatibility wrapper owns only that record |
| `Product` | Target m/z, isolation offsets and CV terms |
| `ScanWindow` | Inclusive begin/end m/z and typed metadata |
| `InstrumentSettings` | Scan mode, zoom flag, polarity, scan windows and metadata |
| `SourceFile` | Filename/path, size in MB, file type, checksum/type, native ID type/accession and CV terms |
| `Acquisition` / `AcquisitionInfo` | Acquisition identifiers, ordered acquisitions, combination method and metadata |
| `Software` / `DataProcessing` | Software name/version/CV terms, processing-action set, optional completion time and metadata |
| `SpectrumSettings` | Type, native ID/comment, instrument/acquisition/source information, precursor/product/processing lists, mobility format/state and metadata |
| `ChromatogramSettings` | Type, native ID/comment, instrument/acquisition/source information, precursor/product, processing list and metadata |

Enums preserve source ordering and exact names for 19 activation methods, 15 scan modes, 22 processing actions, 10 chromatogram types, polarity, checksum types and ion-mobility units/formats/states. Activation parsing accepts either the source full name or short name, case-sensitively. For example, `HCD` means beam-type CID; `HCID` is the separate high-energy CID enum. No invalid `SIZE_OF_*` sentinel can be constructed. `ALL` exposes the supported entries. Defaults follow the source: unknown scan/polarity/spectrum type, no zoom, zero energy/offsets, empty vectors, unknown checksum, unknown mobility format/state, and **mass chromatogram**, not unknown chromatogram.

`kernel::Precursor` now stores acquisition data directly in both spectra and chromatograms. `PrecursorInfo { peak }` remains usable in standalone settings and delegates field access to the same owned record through `Deref`/`DerefMut`; `From` converts in either direction. There is no duplicated acquisition state. `Precursor` is `Clone`, no longer `Copy`, and its constructor is no longer const. Old complete three-field literals require `..Precursor::default()`; old `PrecursorInfo` construction with separate acquisition fields should populate a `Precursor` and convert it.

`isolation_target_mz: None` means the selected m/z; a different target can be retained explicitly. The source-style `isolation_window()` and purity algorithms center their calculations on selected m/z, independently of this interchange field. `spectrum_reference` retains a native scan ID. `MSExperiment::precursor_spectrum_index` searches earlier spectra at one lower MS level for this reference, then falls back to the most recent earlier spectrum at that level. It returns `None` when no parent exists and an error for an invalid input index; it does not infer future parents.

Precursor drift absence uses `None` instead of the C++ -1 sentinel. Finite signed drift values are retained, including negative FAIMS compensation voltages. No implicit mobility conversion or inference is performed. Possible charge states retain input order and duplicates; zero means unknown. The ordinary precursor charge remains independent of that list.

`Precursor::uncharged_mass` (also available through an instance of `PrecursorInfo`) intentionally preserves `Precursor::getUnchargedMass`: unknown charge assumes 2 and the calculation uses the signed charge, `mz * charge - charge * proton_mass`. Thus negative charge yields a signed result. The separate feature decharge calculation uses absolute charge and should be used when that neutral-mass convention is intended. All arithmetic is checked for finite results.

Isolation and drift-window offsets and activation energy must be finite and nonnegative when validated. m/z and intensity follow the existing kernel's finite signed-value policy. Isolation range calculation detects overflow; scan windows require finite ordered bounds and include both endpoints. This is stronger validation than the source's unrestricted setters. Source-file size must be finite and nonnegative. A nonempty declared SHA-1/MD5 digest must have the matching hexadecimal length; an absent digest remains allowed. Validation neither opens a source file nor calculates its digest.

`CompletionTime` is a validated Gregorian wall-clock timestamp, parsed from `YYYY-MM-DD HH:MM:SS` or the equivalent `T` separator and displayed in the source-style space form. It supports years 1–9999 and leap years, rejects leap-second 60 and does not infer a time zone. Absence is `None`.

Whether the day exists is decided by `chrono::NaiveDate::from_ymd_opt`, which uses the proleptic Gregorian calendar. The strict 19-byte parse, the year-zero refusal and the hour, minute and second bounds stay in this crate. `chrono::NaiveDateTime::parse_from_str` with `%Y-%m-%d %H:%M:%S` is not used because it accepts `+1234-01-01 00:00:00`, the unpadded `2024-1-01 00:00:00` and second 60. The chrono call replaced a hand-written month-length table on 2026-09-13 without changing which timestamps are accepted. The crate survey found no difference in any of the 4,620,000 combinations of years 0–9999, months 0–13 and days 0–32. `completion_time_days_follow_the_proleptic_gregorian_calendar` checks this against an independent leap-year rule. It tries every month and day in 16 years that cover each rule, plus the February and 30/31-day boundaries of every year. It passed before and after the change. The check is integer arithmetic, so results are the same on every machine. This is a bounded replacement for the optional timestamp field, not a port of general OpenMS DateTime arithmetic, locale parsing or timezone handling.

## Settings merge and native attachments

`SpectrumSettings::unify` follows the source's asymmetric merge: incoming ordinary metadata overwrites; comments concatenate directly without a delimiter; precursor/product/processing lists append; unequal spectrum types become unknown. The receiving native ID, instrument, acquisition, source and mobility settings remain unchanged. The complete incoming structure is validated before any assignment, so errors are atomic. Records are owned values, replacing shared-pointer storage for processing records.

Equality includes mobility format and peak state, correcting their omission in the pinned SpectrumSettings equality implementation. `SpectrumSettings` also has a native `peptide_identifications` attachment, which validates and appends during unify. The source stores this attachment on `MSSpectrum` rather than on SpectrumSettings itself; the Rust extension allows standalone settings to carry the same records without claiming source-class layout parity.

Standalone settings do not silently synchronize duplicate scalar fields with an existing kernel spectrum/chromatogram, and they are not implicitly flattened into string metadata. Format support remains governed by each format's documented contract. The [mzML adapter](MZML_SUPPORT.md) preserves its documented precursor subset directly through the kernel record. DTA/MGF writers reject nondefault precursor acquisition fields. Other standalone settings do not implicitly expand file-format coverage.

## Verification and provenance

`tests/metadata.rs` covers source empty/string/list/bool values, exact equality and numeric boundaries, units and string bridges, metadata merge policies, CV append/replace/consume and empty-key behavior, enum names and defaults, precursor mass/isolation/drift conventions, scan/checksum/date validation, source settings-unify expectations, nested errors and the native identification attachment.

```sh
cargo test --offline --no-default-features --test metadata
cargo clippy --offline --no-default-features --lib --test metadata -- -D warnings
```

No C++ build was performed. Tests were derived from pinned implementations and class-test expectations rather than a live C++ differential run. Selected implementation SHA-256 values, under source `src/openms/source/`:

| Source | SHA-256 |
| --- | --- |
| `DATASTRUCTURES/DataValue.cpp` | `ab04af951f41dc9b767f73287557d3c5a3f7b7edcacc6611fefffe3f02a2d89b` |
| `METADATA/MetaInfoInterface.cpp` | `4784d93bf3e62f4774f926520f1f9c7b0ab3598e686091c5c83ecf1569368699` |
| `METADATA/MetaInfo.cpp` | `bafa2699024ef964a166234e0706bc0b29d5314abd45ff516046dd568ba42e18` |
| `METADATA/CVTerm.cpp` | `63b99060f543aa90c6ac4f60a25900711b974d38fb01ff1d7728dbf73dfd4649` |
| `METADATA/CVTermList.cpp` | `c064f61d73b4ae0792156a6730d599198630a64b82043686dd889c4f06b67fba` |
| `METADATA/Precursor.cpp` | `009d54ab86dd4dd85d11f61a227bc54c9aef9f061dd97f3aee78568b12350ec1` |
| `METADATA/InstrumentSettings.cpp` | `f412905ada7edc511bb89f974be4a694163e3a247db483999a9cf28757c64632` |
| `METADATA/DataProcessing.cpp` | `283a93d7d36c6bf76da711d19c373daad164d601c8b3810503a5b99cb14da4b4` |
| `METADATA/SpectrumSettings.cpp` | `1c0f999341b221548a751e51adbde3025db84d5a7e18bf25e020ce43a6ec2f7b` |
| `METADATA/ChromatogramSettings.cpp` | `2af3b20decd7a06a4546b60056f4a86f7ea42ef81df73534c646743080db849f` |
