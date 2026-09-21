# MorphologicalFilter support

Source: `src/openms/include/OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h` at
core `bc9cc12` (header-only; `MorphologicalFilter.cpp` is an empty translation
unit), and its consumer `OpenMS4-topp/src/BaselineFilter.cpp` at topp `174b576`.
Rust: `src/processing/baseline.rs`, used by `src/cli/tools/baseline_filter.rs`.
Evidence and fixtures: `tests/data/baseline_filter_edges_provenance.json`.

## API mapping

| C++ member | Rust | Note |
| --- | --- | --- |
| `MorphologicalFilter()` | `MorphologicalFilter::default`, `::new` | `new` validates the element |
| `~MorphologicalFilter()` | — | no resources held |
| `MorphologicalFilter(const MorphologicalFilter&)` (declared, never defined) | `Clone`, `Copy` | the class is stateless between calls |
| `setParameters` / `param_` (`struc_elem_length`, `struc_elem_unit`, `method`) | `structuring_element: StructuringElement`, `method: MorphologicalMethod` | typed fields instead of a `DefaultParamHandler` |
| `filterRange(input_begin, input_end, output_begin)` | `MorphologicalFilter::filter_range` | returns a `Vec<f32>` instead of writing an output range |
| `filter(MSSpectrum&)` | `SpectrumFilter::filter_spectrum` | |
| `filterExperiment(PeakMap&)` | `SpectrumFilter::filter_experiment` | |
| `applyErosion_` / `applyDilation_` (van Herk) | `source_extrema` (private) | monotonic deque, same results (see below) |
| `applyErosionSimple_` / `applyDilationSimple_` | `MorphologicalMethod::ErosionSimple` / `DilationSimple`, and `clipped_extrema` (private) | |
| `struct_size_in_datapoints_`, its `= 0` reset in `filterRange` | `MorphologicalFilter::effective_window` | a computed value instead of a cached member |
| `Internal::IntensityIteratorWrapper`, `intensityIteratorWrapper` | — | slices replace the iterator adapter |
| `ProgressLogger` base, as `filterExperiment` uses it (`startProgress`/`setProgress`/`endProgress`) | `MorphologicalFilter::filter_experiment_with_progress` | the caller passes the logger for the call; the same calls, values (`0` to `n - 1`) and command output as the Release build, replayed in `tests/progress_consumers.rs`; the trait's `filter_experiment` reports nothing |
| — | `BASELINE_PROGRESS_LABEL` | the source label `filtering baseline` |
| — | `MorphologicalFilter::filter_chromatogram`, `StructuringElement::Seconds` | native extension with no source counterpart |

## Preserved source conventions

- **Element conversion.** `Thomson` uses the global average spacing, exactly as
  source `filter`: `ceil(width * (N - 1) / (last_mz - first_mz))`, then an even
  count is rounded up to odd. `DataPoints` is used as a count, and a fractional
  parameter truncates, as the source's `(UInt)(double)` cast does.
- **Short records.** A spectrum with fewer than two peaks keeps its intensities
  and is still marked `Profile`, so a singleton top-hat *spectrum* is unchanged
  while a singleton top-hat *range* is zero. Both follow the source.
- **Signed bottom-hat.** `bothat` is input minus closing, generally
  non-positive, as OpenMS defines it (not the usual closing minus input).
- **The single-sample element.** Source `applyErosion_` and `applyDilation_`
  take van Herk's block method when `size > struc_size && size > 5`. With
  `struc_size == 1` that method writes every output sample **but the last**:
  the lower margin loop starts past the output, the middle blocks cover
  `0 ..= size - 2`, and the higher margin loop is empty. The last sample keeps
  whatever the output buffer held — zero for a spectrum, because
  `filter(MSSpectrum&)` allocates `std::vector<Peak1D::IntensityType>
  output(spectrum.size())` for each call. The results are therefore
  - `erosion`, `dilation`, `opening`, `closing`: the input, last sample zeroed;
  - `tophat`, `bothat`: zero everywhere, last sample kept;
  - `gradient`: zero everywhere, last sample the negated stale buffer value;
  - `erosion_simple`, `dilation_simple`, `identity`: the input unchanged.

  This is what the benchmark hit: `BaselineFilter -struc_elem_length 1` on
  `sub_profile_uk222_first600.mzML` leaves the last peak of 189 of the 198
  centroided MS2 spectra, whose peaks are tens of Thomson apart, and the port
  zeroed it. Spectra of at most five peaks take the direct-window method and
  are zeroed completely, which is the other nine.
- **The shared buffer.** Source `filterRange` keeps a function-local `static`
  buffer, shared by every call in the process and only ever grown with
  zero-filled elements. `opening`, `closing`, `gradient`, `tophat` and `bothat`
  write their intermediate result into it. Only `gradient` reads a stale
  sample: with a one-sample element on more than five samples its last output
  is `0 - buffer[n - 1]`, a value left by an earlier spectrum of the same run.
  `filter_experiment` carries one buffer through the spectra in order, which
  reproduces the `BaselineFilter` tool exactly (executed: the
  `history_gradient_th003` tool run and the `hist_*` experiment groups).
  Spectra with fewer than two peaks never reach `filterRange` and leave the
  buffer untouched. A single `filter_spectrum`, `filter_range` or
  `filter_chromatogram` call starts from an empty buffer, as the first source
  call in a process does; a C++ program that calls the filter repeatedly sees
  the whole process history instead, and that is not reproduced.
- **Arithmetic.** Subtractions are `f32`, as in the source; intensities are
  `f32` and coordinates `f64` throughout.

## Native differences

| Behaviour | Source | Here |
| --- | --- | --- |
| Unsorted positions | documented precondition (`@note The data must be sorted according to ascending m/z!`), unchecked; the filter runs over the sample order | `Error::UnsortedData`, record unchanged. The executed C++ results for the unsorted shapes are still asserted, through `filter_range`, which has no positions |
| Element of zero samples (`DataPoints 0`, a `struc_elem_length` below 1, `Thomson 0`) | `filter` turns the count into 1 | `Error::InvalidValue`. The tests assert that `DataPoints(1)` then reproduces the C++ result exactly |
| Negative `struc_elem_length` | converting a negative double to `UInt` is undefined behaviour in C++, and a negative `Int` element makes `applyErosion_` loop out of bounds | refused by `validate_options` |
| Element of 2^31 samples or more | passed to `applyErosion_` as a negative `Int`, undefined behaviour | `usize` throughout; such an element simply covers the whole signal |
| Zero or non-finite m/z span for a `Thomson` element | divides by zero and converts the result to `UInt` | `Error::InvalidValue` |
| Even element count in `filterRange` | used as given, so van Herk's blocks produce asymmetric windows (measured: 32,628 differing samples in the sweep) | rounded up to odd, as `filter(MSSpectrum&)` does. Documented on `filter_range`; the tests assert the even count equals the executed result for the next odd count |
| Non-finite intensity, overflowing subtraction | computed, storing an infinity or a NaN | `Error::InvalidValue`, record unchanged |
| `filterRange` with a method outside the valid strings | silently writes nothing (`@exception Exception::IllegalArgument` is documented but not thrown) | unrepresentable: `MorphologicalMethod` is an enum |
| van Herk's three comparisons per sample | prefix and suffix blocks plus a `static` scratch buffer of the element length | a monotonic deque of indices, also linear, with the same results |
| Progress logging | `ProgressLogger` base, type `NONE` by default | a caller-owned logger passed to `filter_experiment_with_progress`; the metadata-copy preflight runs before the section, so an experiment refused there prints nothing, and a failure inside the section still ends it (the source cannot fail there) |

`BaselineFilter`'s `method` parameter maps each of the ten source names onto its
own `MorphologicalMethod` in `src/cli/tools/baseline_filter.rs`, `erosion_simple`
and `dilation_simple` included: those two source methods are **not** the same
operation as `erosion` and `dilation` once the element is one sample wide, so
the four names cannot share two variants. Executed on the Release tool over
`tests/data/baseline_filter_edges_edges.mzML` with `-struc_elem_unit DataPoints
-struc_elem_length 1`, `erosion_simple` differs from `erosion` at the last peak
of 11 of the 20 spectra and nowhere else (likewise `dilation_simple`); both runs
are recorded in `baseline_filter_edges_tool_expected.tsv` and asserted by
`simple_method_names_match_the_release_tool`.

## Checked boundaries and evidence

Tier 1, executed against the OpenMS4 Release build (core `bc9cc12`) on
ibminode06; the driver, its build and the run environment are in the provenance
file, and the driver source is
`../oracle/baseline-filter-edges/drivers/morph_edges.cpp`.

- **Edge shapes**, `tests/baseline_filter_edges.rs`: 1,698 executed cases over
  29 shapes and all ten methods — an empty spectrum, one peak, two peaks,
  centroided spectra of 3 to 20 widely spaced peaks, the benchmark's own
  reproducer spectrum, uniform profile spectra, an element wider than the
  spectrum, an element exactly as long as the spectrum, one sample below it,
  plateaus and ties, and unsorted input in three shapes. Elements are given in
  Thomson and in DataPoints, including fractional and zero lengths. Intensities
  are compared bitwise.
- **Buffer history**, same file: twelve `filterExperiment` groups (one per
  method, plus a reversed and a `DataPoints 1` ordering) reproduce the source's
  stale-buffer values sample for sample.
- **Sweep**, same file: 51,780 executed cases over 3,341,420 samples (signal
  lengths 0 to 48 against every element length up to `2n + 3`, and lengths 97,
  128, 255, 256, 1000 and 4097 against fifteen element lengths) recorded every
  sample where the Release build departs from a clipped-window filter. There
  are 528 such samples after the even-length `filterRange` rows are set aside,
  all of them the last sample of a one-sample-element case, and the port
  reproduces exactly that set: no sample more, no sample fewer.
- **Tool**, `tests/topp_baseline_filter_edges.rs`: eight Release
  `BaselineFilter` runs over the same mzML inputs — the reproducer with the
  benchmark INI, the edge shapes with the benchmark INI, with the tool defaults
  and with `-method erosion`, `-method erosion_simple` and `-method
  dilation_simple` at `-struc_elem_unit DataPoints -struc_elem_length 1`, and
  the history file with `gradient` and with `tophat`. The first six were
  re-executed unchanged when the two `*_simple` runs were added.
- **Scale**, same file, `#[ignore]` because they read `/ceph` on the IBMI
  nodes: the benchmark's 600-spectrum UK222 slice (30 MB) and the whole 2.3 GB
  UK222 run (40,856 spectra) against the Release tool's output, every intensity
  bitwise.
- **Upstream class test**: `tests/baseline.rs` keeps the complete 40-row
  `MorphologicalFilter_test_1.txt` table and the class test's direct-window
  comparison across signal lengths and element lengths. The upstream test's
  element lengths start at 3, which is why the single-sample behaviour was
  never visible there.

Bounded work: the filters allocate one output vector per record and one shared
buffer per experiment, both the length of the record; the element length never
drives an allocation, so `DataPoints(usize::MAX)` is answered with the clipped
window over the whole signal. `filter_experiment` builds the filtered spectra
in a copy and commits them only on success, and meters that copy with the
processing acquisition-copy budget **per spectrum**, where the trait default
meters one ledger for the whole experiment. A single ledger makes the ceiling
shrink as the run grows, which rejects real data: the full 2.3 GB UK222 run
(40,856 spectra) fails the shared ledger with "data array description resource
limit exceeded" although every one of its spectra is far inside the per-record
allowance, and the metadata being copied is already held in memory. The same
shared ledger is used by `SpectrumFilter`'s default `filter_experiment` in
`src/processing.rs`, so the other spectrum filters still have this limit; that
is an integrator request, not a change made here.
