# Acquisition settings through processing

Native spectrum and chromatogram transformations retain the attached acquisition
fields described in [ChromatogramTools support](CHROMATOGRAM_TOOLS_SUPPORT.md).
This extends the existing processing adapters; it introduces no new public
processing interface or numerical algorithm.

Both record types own InstrumentSettings, AcquisitionInfo and SourceFile, and
share DataProcessing records through `Vec<Arc<DataProcessing>>`. Spectrum outputs
also retain Products. Chromatogram outputs retain their Product and
ChromatogramType. A copied output can change its owned strings or acquisition
entries without changing its input. Processing handles retain their allocation
identity, so cloning a spectrum does not duplicate a referenced processing graph.
These rules also apply to empty records and to unselected records retained by an
experiment operation.

## Operations and ownership

All existing processing full-record copies are covered, including:

- The default `SpectrumFilter::filter_experiment`, which copies spectra and
  retains the original chromatograms.
- Gaussian and Savitzky–Golay experiment filtering, which copy the experiment and
  smooth both spectra and chromatograms.
- WindowMower's owned result and experiment operation.
- Deisotoper and AveragineDeisotoper, including their empty-input shortcuts and
  experiment operations.
- PeakPickerHiRes spectrum, chromatogram and experiment outputs.
- PeakPickerIterative spectrum and experiment outputs.
- PeakPickerChromatogram's returned smoothed and picked chromatograms.

Intensity-only and selection operations continue to change their existing
scientific fields in place. No reset of acquisition settings or new processing
history entry is invented. Existing choices about spectrum representation,
retained/omitted auxiliary arrays, selected MS levels, centroids, smoothing and
integration are unchanged. Iterative peak picking still constructs a temporary
peak-only seed spectrum; it derives the final output settings from the original
input, as before.

The source correspondence is settings ownership, not a new promise of full C++
object or ABI identity. For example, the pinned
[PeakPickerHiRes implementation](https://github.com/okohlbacher/OpenMS4-core/blob/54a232fe2cae9c590d5c997fa49d20e7769860fb/src/openms/source/PROCESSING/CENTROIDING/PeakPickerHiRes.cpp)
copies spectrum metadata and the complete ChromatogramSettings portion when
preparing picked outputs. The native flat fields represent the supported settings
from those source records. Existing metadata-model and processing differences
remain documented in the individual algorithm support pages.

## Copy accounting and failure behavior

One private acquisition-copy ledger follows each top-level bundled operation.
Its fixed ceilings are **50 million weighted work visits** and **256 MiB of
cumulative logical copied payload**. These are additional bounds for the newly
attached acquisition fields. They do not change the meaning or defaults of an
algorithm's numerical `max_work`, peak-count, isotope or sampling limits, and do
not replace earlier limits for other metadata and arrays. They are not exact
allocator or process-resident-memory ceilings.

The preflight accounts for owned names, values, units, CV terms, scan windows,
acquisition entries and spectrum Products before each full-record clone. It
charges DataProcessing vector handle slots without traversing or copying the
shared processing contents. Record validation has a separate bounded traversal
of those contents because validation reads them. Public Clone, equality and
caller destruction retain ordinary Rust behavior outside a processing call.

Copy multiplicity follows the implementation, including temporary copies. An
ordinary direct Deisotoper output costs one acquisition copy. A nonempty
AveragineDeisotoper call costs two: its selected working spectrum and final
output; its empty shortcut costs one. An experiment's initial copy and each
subsequent picked or deisotoped replacement all consume the same ledger.
PeakPickerChromatogram charges both its smoothed copy and the nested HiRes
output. Iterative picking shares the ledger with its HiRes seed call and final
output; the internal seed has default acquisition settings, but its fixed record
fields are still charged. Limits do not reset for every record or nested picker.

Copy preflight occurs before the corresponding clone, sometimes after numerical
work has already succeeded. Any checked failure preserves the caller's spectrum,
chromatogram or complete experiment. Local staged results are dropped rather
than published. A batch failure can therefore happen after an earlier staged
record was successfully processed, without changing the original input.

The default trait method can account for its own initial spectrum-vector copy.
A caller-defined `SpectrumFilter` implementation remains responsible for work
and extra copies inside its callback; the private ledger is not an extensible
public execution framework. Bundled filters that produce further full-record
copies override or delegate through the shared internal path.

## Verification

[Four focused public tests](../tests/processing_acquisition.rs) check direct and
empty outputs, all affected bundled experiment filters, nested chromatogram
picking, independent owned-field edits, shared processing identity, and rollback
when a later spectrum fails. Three private tests in
[processing.rs](../src/processing.rs) use small deterministic budgets to check
cumulative copies across records and nested calls, byte exhaustion before
cloning, and charging shared handles without deep-copying their payload.
The tests exercise the native ownership and resource extension; they are not
executed C++ differential tests. Existing source-derived numerical tests remain
in the individual processing suites.

See [signal processing](SIGNAL_PROCESSING.md),
[WindowMower](WINDOW_MOWER_SUPPORT.md), [deisotoping](DEISOTOPING_SUPPORT.md),
[averagine deisotoping](AVERAGINE_DEISOTOPING_SUPPORT.md),
[HiRes picking](PEAK_PICKING_SUPPORT.md),
[iterative picking](ITERATIVE_PICKING_SUPPORT.md), and
[chromatogram picking](CHROMATOGRAM_PICKING_SUPPORT.md) for unchanged scientific
behavior and source fixture provenance.
