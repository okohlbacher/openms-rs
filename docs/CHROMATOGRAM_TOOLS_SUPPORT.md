# Chromatogram conversion and attached acquisition settings

The two `kernel::ChromatogramTools` operations implement the complete concrete
MSExperiment conversion surface from Core SDK
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. See
[the implementation](../src/kernel/chromatogram_tools.rs),
[tests](../tests/chromatogram_tools.rs), and
[source hashes](../tests/data/chromatogram_tools_provenance.json).
No C++ program was built or executed for this increment.

## Attached fields

Both `MSSpectrum` and `MSChromatogram` now own `instrument_settings:
InstrumentSettings`, `acquisition_info: AcquisitionInfo`, `source_file:
SourceFile`, and `data_processing: Vec<Arc<DataProcessing>>`. Spectra also have
`products: Vec<Product>`; chromatograms have `chromatogram_type:
ChromatogramType`, defaulting to `Mass`. These reuse existing native metadata
records, with no second copy of RT, MS level, precursor, native ID or spectrum
type. The independent SpectrumSettings/ChromatogramSettings records remain
standalone; no automatic synchronization or full settings bridge is introduced.
Comments and scan mobility attachment remain separate gaps.

**Public struct literals must use `..Default::default()` or supply the new
fields.** Existing constructors populate source defaults. Clone copies owned
acquisition data and shares the processing handles; equality compares all stored
fields, including processing contents through Arc, rather than pointer identity.
Ordinary Rust swap moves the entire record. Existing `clear(false)` removes peaks
and all parallel arrays while retaining acquisition data; `clear(true)` resets
the complete native record. No new clear implementation is introduced.

`has_acquisition_settings()` is an O(1) unsupported-transport predicate: it checks
scalar fields and collection lengths without visiting referenced metadata. A
SourceFile size of negative zero is treated as stored information here, although
source and native numerical equality equate it with positive zero. The mzML
writer rejects any of the new nondefault fields before output until their
transport is implemented; see [XML loss guards](MZML_ACQUISITION_GUARDS.md).

## Public conversion methods

`ChromatogramTools { limits: ChromatogramConversionLimits }` is defaultable.

- `convert_chromatograms_to_spectra(&mut MSExperiment)` appends spectra and
  removes all input chromatograms.
- `convert_spectra_to_chromatograms(&mut MSExperiment, remove_spectra,
  force_conversion)` appends grouped chromatograms; both flags default to false
  in source and are explicit native arguments.

Both return `Result<ChromatogramConversionReport>`. The report provides added
spectrum/chromatogram counts, the count of selected but unconvertible SRM scans,
and **unchanged removed records by ownership** in `removed_spectra` and
`removed_chromatograms`. The caller chooses when to release or reuse their full
annotations. This avoids cloning or recursively dropping unrelated record graphs
inside the checked conversion. Ignoring/dropping the result has ordinary Rust
ownership cost.

### Chromatograms to spectra

Each input point produces one MS2 spectrum, in chromatogram order followed by
point order. Its RT/intensity come from the point, and its single m/z is the
chromatogram Product m/z. One complete Precursor and Product, InstrumentSettings,
AcquisitionInfo and SourceFile are copied from the chromatogram. SRM/SIM
chromatogram types override scan mode to their corresponding spectrum modes;
other types retain the copied mode.

Source does not copy native ID, name, general metadata, processing history or
parallel arrays into these new spectra. Those remain available on the removed
originals returned in the report. Existing spectra and experiment metadata are
retained. Empty chromatograms produce no scans but are still removed; when no
point is added, the existing spectrum buffer is unchanged.

### Spectra to chromatograms

Without forcing, only SRM scan mode is selected; MS level is not consulted. A
selected spectrum with exactly one precursor and at least one peak contributes
one point per peak to the exact `(precursor m/z, peak m/z)` SRM group. Multi-peak
pseudo-spectra are split logically, avoiding the source's unnecessary full
spectrum copies. Duplicate peak keys contribute repeated points.

When forced, all scan modes are selected. A one-precursor nonempty scan still
uses the SRM branch. Every other selected scan contributes XIC points grouped by
peak m/z, including scans with **more than one precursor**. The source comment
mentions missing precursors, but its actual fallback branch is broader.

New XICs are appended first, ordered by m/z; new SRM groups follow in ascending
precursor then product m/z order. Existing chromatograms stay as a prefix. Within
each group points retain encounter order, including unsorted RTs and exact
duplicates. No tolerance or RT sort is applied. Signed-zero keys compare equal;
the first key/payload retains its sign.

An SRM chromatogram copies its first contributing scan's complete first Precursor,
InstrumentSettings, AcquisitionInfo and SourceFile. Its Product is newly created
from the grouping m/z, ignoring any spectrum Products. It has SRM type and native
ID `chromatogram=` plus the first scan's native ID. XICs copy the same acquisition
fields from their first scan; their new precursor m/z is the group key, with
metadata `description = "XIC @ " + StringUtils-compatible float text`. Their
Product, native ID and chromatogram type remain source defaults (Mass, not an
invented selected-ion type). Generated processing history and general annotations
remain empty in both branches.

`remove_spectra=true` removes **every original SRM-mode scan**, including skipped
or empty ones, and retains forced non-SRM scans even when they were converted.
Removed scans are returned unchanged. Unselected/unconsumed annotations are not
validated or copied by the converter. Skipped-count reporting replaces source
warning side effects; ordinary unselected spectra are not counted as warnings.

## Checked boundaries and limits

Only consumed coordinates and intensities must be finite. Negative finite values
are allowed. Nonfinite map keys, RTs or emitted intensities produce a checked
atomic error; source ordered maps cannot safely order NaN keys. Copied metadata
is retained as stored, without speculative validation of unused values.

Defaults allow one million input records, one million final destination records
including an existing prefix, ten million generated points, 50 million weighted
work visits, and 256 MiB cumulative logical new allocation. Limits are public and
configurable per conversion. Metadata names, values, units, CV terms and nested
acquisition records are charged before every deep copy. Repeated source metadata
copies share the operation's counters. Map nodes, comparison work, point-vector
growth, destination slots and removed-record slots are precharged. No source
input controls recursion. Unchanged prefix records move without deep cloning.
No checked error mutates the experiment; only after staging succeeds are buffers
moved into their final owners. These are logical payload bounds, not a promise
about allocator bookkeeping or process resident memory.

Existing full-record processing transformations retain the new fields; see
[processing acquisition support](PROCESSING_ACQUISITION_SUPPORT.md). They use
a separate acquisition-copy allowance (50 million work visits / 256 MiB), shared
across actual copies within each operation, without changing their numerical
algorithm budgets. Processing Arc copies charge handle slots only. General
record validation separately meters shared processing contents before reading
them, under a per-record 50 million / 256 MiB validation allowance. This does not
replace earlier numerical or unrelated metadata resource contracts. Ordinary
public Clone, equality, swap, clear and caller drop have ordinary Rust costs.

TSG append and EMG fit receive the same bounded acquisition-copy preflight.
Concrete 2D import resets new output settings and returns the complete old map;
plain export and TIC/XIC summaries retain their source omission/default behavior.
Existing range and selection operations modify point data while retaining the
attached settings. No scientific algorithm or sample precision changes.

## Verification

The 13 direct conversion/ownership tests and 74 adjacent kernel, 2D, summary,
aggregation, theoretical-spectrum and EMG tests pass on Rust 1.98 with all
features and Rust 1.85 without default features. Scoped strict library/test
Clippy passes on both compilers, and changed-file rustfmt checks pass. The XML
guards and processing copy paths have separate focused suites documented in
their linked support pages. Independent source and resource reviews closed
without outstanding findings. All 23 manifest hashes were checked against the
pinned source checkout.
