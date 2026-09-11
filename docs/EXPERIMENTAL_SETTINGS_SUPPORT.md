# ExperimentalSettings and experiment ownership

The native [aggregate](../src/metadata/experimental_settings.rs) represents every
field and class-specific operation of `METADATA/ExperimentalSettings` at reduced
SDK `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. Public owned fields replace source
getters/setters, ordinary Clone and moves replace copy/assignment, and Display
emits the exact two diagnostic marker lines. This is an owned value group;
complete mzML header transport and streaming consumers are separate work.

`ExperimentalSettings` owns the document identifier/provenance, Sample and nested
subsamples, source files, ordered contacts, default instrument, additional
instrument configurations keyed by ID, HPLC/Gradient, complete DateTime, comment,
fraction identifier and typed MetaInfo. Source defaults are preserved, including
HPLC temperature 21 and invalid zero-field DateTime. Equality uses all these
values; inherited DocumentIdentifier equality compares only the identifier and
ignores loaded path/type. Clone and moves retain that provenance nevertheless.
Floating fields retain their bits and source IEEE PartialEq behavior. Resource
validation adds no physical positivity, finite-value, calendar completeness,
component-order or Gradient-shape restrictions to this owned aggregate.

## Single metadata owner and migration

`MSExperiment` now has `settings: ExperimentalSettings`. The old independently
owned `metadata: BTreeMap<String,String>` field is removed. The authoritative
run metadata is `experiment.settings.metadata`, whose values are `MetaValue`.
There is no second mutable map or precedence/synchronization rule. Existing
callers must update experiment metadata accesses and explicit struct literals:

```rust
use openms::{MSExperiment, metadata::ExperimentalSettings};
let mut experiment = MSExperiment {
    settings: ExperimentalSettings::default(),
    ..Default::default()
};
experiment.settings.metadata.insert("sample".into(), "sample A".into());
assert_eq!(experiment.settings.metadata["sample"].as_str()?, "sample A");
# Ok::<(), openms::Error>(())
```

Spectrum/chromatogram scalar metadata stays at its existing location; this
migration changes the experiment owner only. Callers converting a preexisting
string map can use `metadata::meta_from_strings` explicitly without guessing
numeric types or units.

`clear(false)` clears spectra/chromatograms and preserves all settings.
`clear(true)` resets the complete aggregate. Whole experiment Clone and standard
swap retain every field. Assigning `experiment.settings` replaces only settings.
Checked 2D imports replace the entire experiment and return previous ownership,
including original metadata buffers; they do not traverse old settings. Sorting,
range queries, in-place alignment and ChromatogramTools operate on contents and
leave settings untouched.

Spectrum-only filter batches continue to stage/replace just spectra. Gaussian,
Savitzky–Golay, HiRes and iterative whole-experiment outputs clone all settings;
their existing shared `AcquisitionCopies` ledger now charges the aggregate once
for each actual whole-experiment clone. Nested picker calls keep the same
counters. Their numerical budgets are unchanged.

## Checked resource operations

`validate`, `validate_with_limits`, `checked_clone`, and
`checked_clone_with_limits` traverse the full aggregate before allowing a
checked clone. Defaults are 1,000,000 structural records, sample depth 64,
50,000,000 work units and 256 MiB of cumulative payload/scratch allowance.
`ExperimentalSettingsLimits` can adjust those limits, with a hard depth ceiling
of 64 before recursive derived Clone. Metadata/list entries also consume work
and bytes. Source files, all component metadata, Software/CV payload, units,
string/numeric lists and the Gradient's complete stored percentage rows are
included, even when its axes were cleared and storage is stale.

The sample walk is iterative; scratch grows with depth, not sibling count.
Record/vector/map bounds are charged before their payloads are inspected.
Sparse BTree storage and cloned string/vector buffers use the existing shared
metadata meter. Callers copying several aggregates can share counters through
the private helper; processing does so. Metadata lists and units retain owned
identity through cloning. Ordinary public field construction, Clone, equality,
clear and Drop keep ordinary Rust costs, including caller-created recursive
sample trees; checked operations reject excess depth before recursive cloning.

## Current transport boundary

mzML run userParams now use the existing typed scalar codec, preserving string,
integer, float, and supported MS/UO unit values. Spectrum/chromatogram userParams
retain their existing textual contract. Numeric input follows the existing
source XSD scalar parsing rules; unrecognized XSD kinds remain text. Empty and
list MetaValues return Unsupported before writing until the complete source
header metadata codec is implemented. Reserved record-name keys remain errors.

Both ordinary and Numpress mzML writers reject nondefault unrepresented header
settings before output. The predicate inspects only fixed scalar fields and
container lengths. It distinguishes stored negative zero, partial DateTime,
stale Gradient rows and nonempty default-valued lists. It does not traverse
caller-controlled trees to compare them with Default. Persistent document
identifier is included in the loss guard; loaded path/type are local provenance
and are excluded. A normal loaded experiment can therefore be written again.
Run scalar metadata validation/rendering has a conservative shared work/byte
precharge, independent of binary-array encoding budgets.

Successful `mzml::load*` and FileHandler path loads populate document loaded
path and content-detected type using the shared DocumentIdentifier operations.
Stream readers leave path/type unset. Relative paths become absolute by the
existing lexical filesystem policy; type recognition reads a bounded 64 KiB
plain/decompressed prefix separately from parser limits. Document identifier
remains independent. Publication into a destination remains atomic on errors.

DTA/MS2/DTA2D and MGF experiment writers reject unrepresented persistent settings
and run metadata before emitting bytes. This closes the former MGF run-metadata
loss path. DTA2D `write_tic`/`store_tic` remain explicitly lossy signal-only
projections: unconsumed settings are preserved in the caller and ignored by the
projection. Full header read/write mappings, arbitrary list/Empty mzML metadata,
additional instrument references and consumer setup remain explicit gaps; this
aggregate does not imply their implementation.
