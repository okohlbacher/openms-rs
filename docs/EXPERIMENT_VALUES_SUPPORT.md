# Sample and instrument values

`openms::metadata` provides the complete stored values and class-specific public operation groups for OpenMS4-core `Sample`, `IonSource`, `MassAnalyzer`, `IonDetector` and `Instrument`, audited at revision `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. They are independent owned values. This batch does not attach them to `MSExperiment`, implement `ExperimentalSettings`, or extend mzML header transport.

## Public mapping

| Native value | Complete owned state |
| --- | --- |
| `Sample` | `name`, `organism`, `number`, `comment`, `state`, `mass`, `volume`, `concentration`, ordered recursive `subsamples`, `metadata` |
| `IonSource` | `inlet_type`, `ionization_method`, reused `Polarity`, signed `order`, `metadata` |
| `MassAnalyzer` | `analyzer_type`, `resolution_method`, `resolution_type`, `scan_direction`, `scan_law`, `reflectron_state`, `resolution`, `accuracy`, `scan_rate`, `scan_time`, `tof_total_path_length`, `isolation_width`, `final_ms_exponent`, `magnetic_field_strength`, signed `order`, `metadata` |
| `IonDetector` | `detector_type`, `acquisition_mode`, `resolution`, `adc_sampling_frequency`, signed `order`, `metadata` |
| `Instrument` | `name`, `vendor`, `model`, `customizations`, ordered `ion_sources`, `mass_analyzers`, `ion_detectors`, owned `Software`, `ion_optics`, `metadata` |

Public fields and ordinary vector/map operations replace source scalar setters and const/mutable getters. `Default`, `Clone`, Rust moves and `PartialEq` provide default construction, independent copying/assignment and equality/inequality. Every source-owned field participates in equality, including metadata, nested software, ordered subsamples and component order. Component `order` is independent of its position in a vector; setting it does not sort the instrument.

Defaults are empty strings/collections/metadata/software, `Unknown` enums, positive floating zero and signed integer zero. Sample units are grams, milliliters and grams per liter; analyzer documentation retains source units and definitions. Scalars are stored unchanged. Negative values, infinities, signed zeros and NaNs are permitted, matching source assignment and direct `double ==` behavior. A NaN in a scalar field makes that record unequal to its own clone; both zero signs compare equal. Integer component orders and the final MS exponent retain the full `i32` domain.

`MetaInfo` and `Software` reuse the existing native representations. Their established exact typed metadata comparison, unit identity and validation policies apply; this does not reintroduce the source DataValue scalar epsilon comparison. See [metadata support](METADATA_SUPPORT.md).

## Enum names

Every enum has `ALL` in source order, `name()`, `Display`, and exact case-sensitive `FromStr`. `ALL.iter().map(|value| value.name())` replaces allocating source name-list helpers. Unknown strings return `Error::InvalidValue`; parsing does not trim or accept symbolic enum aliases. Invalid source `SIZE_OF_*` sentinels are excluded from the Rust types. These value types do not promise C++ enum ABI layout.

| Family | Names |
| --- | ---: |
| `SampleState` | 7 |
| `InletType` | 20 |
| `IonizationMethod` | 52 |
| reused `Polarity` | 3 |
| `AnalyzerType` | 15 |
| `ResolutionMethod` | 4 |
| `ResolutionType` | 3 |
| `ScanDirection` | 3 |
| `ScanLaw` | 4 |
| `ReflectronState` | 4 |
| `DetectorType` | 22 |
| `DetectorAcquisitionMode` | 5 |
| `IonOpticsType` | 12 |

The 154 names preserve exact source spellings, including `Membrane sparator`, `Collsion induced decomposition`, `Collsiona activated decomposition` and `Linar`. Most families use `Unknown`; source polarity names are lowercase `unknown`, `positive`, `negative`. `DetectorAcquisitionMode` is separate from mzML's scan materialization options. The native optics `FromStr` is a convenience beyond the source instrument's name-list helper.

## Hashing and ownership

The three component types implement Rust `Hash` with the fields selected by their source specializations. `IonSource` includes metadata; `MassAnalyzer` and `IonDetector` omit metadata even though it remains part of equality. Equal positive/negative zeros hash identically. Hashing does not validate or change a stored scalar. Numeric Rust hashes are hasher-dependent and are not C++ digest values or a persistence format. No `Eq`/total ordering is asserted for these complete records, which may contain IEEE partial values. Sample and Instrument have no new hashing operation.

Public construction, field assignment, `Clone`, vector mutation and recursive Sample equality/drop retain normal Rust allocation/stack costs. This value batch adds no arbitrary physical limits or hidden budget ledger. A future untrusted parser or bounded aggregate-copy operation must meter its consumed payload and subsample depth before copying. There is no reference sharing between cloned owned sample/instrument records; ordinary contained metadata and software clone independently.

## Verification and source boundaries

[The 154-row TSV](../tests/data/experiment_value_names.tsv) is a literal extraction from source `NamesOf*` tables paired with header enum symbols. It stores source paths and declaration line numbers, and includes all entries before each SIZE sentinel. Tests compare every entry, its position and each complete native `ALL` length; expected names are not generated from Rust.

[Seven focused tests](../tests/experiment_values.rs) also use published class-test values for sample numbers and strings, component values `47.11`–`47.17`, detector sampling frequency `47.21`, order `45`, software `sn`, and instrument component ordering. Independent tests cover each field's equality participation, deep ownership, unsorted signed component orders, scalar bit retention, NaNs, infinities and zero hashing. Hash checks distinguish field inclusion from a portable digest oracle.

The [provenance manifest](../tests/data/experiment_values_provenance.json) records all five source headers, implementations and class tests, the reused metadata dependencies and the literal fixture hash. No C++ build or runtime differential execution was performed for this batch. Rust 1.98/all-features and Rust 1.85/no-default-features targeted tests and strict scoped lint are recorded in the isolated integration handoff.

Complete experiment aggregation, general DateTime, DocumentIdentifier integration, controlled-vocabulary registry resolution and header/consumer transport remain separate operation groups. No format support is inferred by adding these standalone values.
