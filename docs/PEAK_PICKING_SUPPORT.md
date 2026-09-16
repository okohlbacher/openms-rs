# High-resolution peak picking

`processing::peak_picking` ports `PeakPickerHiRes` and hosts the median
signal-to-noise estimator it uses, whose API mapping, conventions, native
differences and evidence are in [SIGNAL_TO_NOISE_SUPPORT.md](SIGNAL_TO_NOISE_SUPPORT.md).
The reference is OpenMS4-core commit
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` (re-pinned from `7c029e8`; the diff
between the two is empty for every file below). No C++ library is called or
built by the crate.

| Source | Rust |
|---|---|
| `PROCESSING/CENTROIDING/PeakPickerHiRes.h`, `PeakPickerHiRes.cpp` | `src/processing/peak_picking.rs` |
| `PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h` (a header template; its `.cpp` only defines a global default object) | `src/processing/peak_picking/noise.rs` (see `SIGNAL_TO_NOISE_SUPPORT.md`) |
| `PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h`, `SignalToNoiseEstimator.cpp` (the base and `estimateNoiseFromRandomScans`) | `src/processing/noise_estimation.rs` (see `SIGNAL_TO_NOISE_SUPPORT.md`) |
| `MATH/MISC/CubicSpline2d`, `MATH/MISC/SplineBisection.h` | `src/processing/spline/` (read-only here; see `CUBIC_SPLINE2D_SUPPORT.md` and `SPLINE_BISECTION_SUPPORT.md`) |

```rust
use openms::processing::peak_picking::{FwhmUnit, PeakPickerHiRes, PickingCompatibility};
use openms::MSSpectrum;

fn centroid(profile: &MSSpectrum) -> openms::Result<MSSpectrum> {
    let picker = PeakPickerHiRes {
        signal_to_noise: 1.0,
        report_fwhm: Some(FwhmUnit::Absolute),
        ..Default::default()
    };
    let result = picker.pick_spectrum(profile)?;
    // result.boundaries follows the centroid order.
    // result.omitted_arrays identifies profile annotations without an aggregation rule.
    Ok(result.spectrum)
}

// A TOPP tool reproducing the C++ output maps the `algorithm` subsection and
// adopts the source's acceptance of degenerate input.
fn from_tool_parameters(algorithm: &openms::param::Param) -> openms::Result<PeakPickerHiRes> {
    let mut picker = PeakPickerHiRes::from_param(algorithm)?;
    picker.compatibility = PickingCompatibility::source();
    Ok(picker)
}
```

`pick_spectrum`, `pick_chromatogram` and `pick_experiment` return newly owned
results with boundaries and the names of omitted profile arrays.
`SpectrumFilter::filter_spectrum`, `filter_experiment` and `filter_chromatogram`
replace the input only after the whole operation succeeds; they discard the
omission report. `pick_experiment_in_place` centroids a whole run without
building a second experiment, for a caller that does not need the profile data
afterwards; it gives up that atomicity in exchange.

## API mapping

### `PeakPickerHiRes`

| Source member | Rust | Notes |
|---|---|---|
| `PeakPickerHiRes()` | `PeakPickerHiRes::default()` | The typed defaults equal `getDefaults()`; checked against the stored C++ defaults. |
| `~PeakPickerHiRes()` | `Drop` (implicit) | |
| `struct PeakBoundary { mz_min, mz_max }` | `PeakBoundary { min, max }` | Retention times for chromatograms, as in the source. |
| `pick(const MSSpectrum&, MSSpectrum&)` | `pick_spectrum` | Returns boundaries too. |
| `pick(const MSSpectrum&, MSSpectrum&, boundaries, check_spacings = true)` | `pick_spectrum_with_spacing` | |
| `pick(const MSChromatogram&, MSChromatogram&)` | `pick_chromatogram` | |
| `pick(const MSChromatogram&, MSChromatogram&, boundaries, check_spacings = false)` | `pick_chromatogram_with_spacing` | |
| `pick(const Mobilogram&, Mobilogram&[, boundaries, check_spacings])` | not ported | No TOPP tool on the bundle path calls it. |
| `pickExperiment(const PeakMap&, PeakMap&, check_spectrum_type = true)` | `pick_experiment` | `check_spectrum_type` is a field. |
| `pickExperiment(input, output, boundaries_spec, boundaries_chrom, check_spectrum_type)` | `pick_experiment` | `spectrum_boundaries` holds one `Option` per input spectrum; the source appends boundaries for picked spectra only. |
| `pickExperiment(OnDiscMSExperiment&, PeakMap&, check_spectrum_type)` | not ported | Low-memory processing is package P4. That overload sorts each spectrum and consults the stored type only (`getType()` without data), unlike the in-memory overload. |
| — | `pick_experiment_in_place` | Native streaming form of `pick_experiment`: same records, same rules, bit-identical centroids, written back over the input so the profile samples are released per record. Returns `PickedExperimentReport` and is not atomic. |
| `pick_` (protected template) | private `pick_signal` over `Peak1D` and `ChromatogramPeak` | |
| `updateMembers_` | `PeakPickerHiRes::from_param` | |
| members `signal_to_noise_`, `spacing_difference_gap_`, `spacing_difference_`, `missing_`, `ms_levels_`, `report_FWHM_`, `report_FWHM_as_ppm_`, `allow_missing_flank_` | public fields `signal_to_noise`, `spacing_difference_gap`, `spacing_difference`, `missing`, `ms_levels`, `report_fwhm` (with `inactive_fwhm_unit`), `allow_missing_flank` | Zero spacings stay zero in the fields; the picker treats them as the source's infinity. |
| `DefaultParamHandler::getDefaults` | `PeakPickerHiRes::defaults` | |
| `DefaultParamHandler::setParameters` | `PeakPickerHiRes::from_param`, `from_param_with_warnings` | Type and restriction violations are `Error::InvalidValue`; unknown names are warnings. |
| `DefaultParamHandler::getParameters` | `PeakPickerHiRes::to_param` | |
| `ProgressLogger` base | not ported | No progress output. |
| — | `PickingCompatibility`, `ion_mobility_array`, `max_points`, `max_work`, `max_metadata_per_record`, `omitted_arrays` | Native. |
| — | `PEAK_PICKER_HI_RES_NAME`, `CENTROIDED_INPUT_MESSAGE` | The handler name and the source exception text, for tools. |

### `SignalToNoiseEstimatorMedian` and its base

[SIGNAL_TO_NOISE_SUPPORT.md](SIGNAL_TO_NOISE_SUPPORT.md#api-mapping) maps all 13
public declarations of `SignalToNoiseEstimatorMedian.h`, its protected members,
the inherited `DefaultParamHandler` and `ProgressLogger` members it uses, and
every member of `SignalToNoiseEstimator.h` and `SignalToNoiseEstimator.cpp`.
The picker holds its estimator as `PeakPickerHiRes::noise_estimator` and passes
its own `compatibility`, whose `noise` field (`NoiseCompatibility`) selects the
estimator's source behaviours.

## Parameter contract

`PeakPickerHiRes::defaults()` builds, in declaration order: `signal_to_noise`
(double `0.0`, min `0.0`), `spacing_difference_gap` (double `4.0`, min `0.0`,
advanced), `spacing_difference` (double `1.5`, min `0.0`, advanced), `missing`
(int `1`, min `0`, advanced), `ms_levels` (empty int list, min `1`),
`report_FWHM` (string `false` restricted to `true,false`), `report_FWHM_unit`
(string `relative` restricted to `relative,absolute`), `allow_missing_flank`
(string `false` restricted to `true,false`, advanced) and the subsection
`SignalToNoise:` with `max_intensity` (int `-1`, min `-1`, advanced),
`auto_max_stdev_factor` (double `3.0`, `0.0..999.0`, advanced),
`auto_max_percentile` (int `95`, `0..100`, advanced), `auto_mode` (int `0`,
`-1..1`, advanced), `win_len` (double `200.0`, min `1.0`), `bin_count` (int `30`,
min `3`), `min_required_elements` (int `10`, min `1`), `noise_for_empty_window`
(double `1e20`, advanced) and `write_log_messages` (string `true` restricted to
`true,false`). Descriptions are the source strings, including the double space
after "too small as well." in `max_intensity`.

Evidence: the tree equals the product SDK's `ParamXMLFile::store` of
`PeakPickerHiRes().getDefaults()` read back (`tests/data/peak_picking/defaults.ini`),
the estimator's own defaults equal `noise_defaults.ini`, and the `algorithm`
node of the retained `WRITE_INI_OUT.ini` (`TOPPWRITEINI_OVERWRITE`) equals the
defaults with `signal_to_noise = 1.0`, the value that test's update keeps from
`WRITE_INI_IN.ini`. The writer shows `report_FWHM` and `allow_missing_flank` as
`bool` and `write_log_messages` as `string`: the source writes `bool` only for a
`true,false` string whose value is `false`.

`from_param` follows `DefaultParamHandler::setParameters`: missing values take
their defaults, a wrong value type or a violated restriction is an error, an
unknown name is a warning (`from_param_with_warnings`). `report_FWHM_unit`
other than `absolute` is ppm, as `updateMembers_` tests `!= "absolute"`.
`SignalToNoise:auto_mode` `0` is the standard-deviation range, `1` the
percentile range and every other value (only `-1` passes the restriction)
the manual range. `to_param(from_param(p))` equals `p` completed with the
defaults, including values the selected modes ignore.

## Preserved source conventions

- **Peak cores.** Sample `i` with `2 <= i < n - 2` is a core when both
  neighbours have magnitude at least `f64::EPSILON`, it is a strict maximum and
  the three samples reach the signal-to-noise threshold. Fewer than five
  samples give an empty record, but the ion mobility and FWHM arrays are still
  created, empty, as the source creates them before its size check.
- **Spacing.** `min_spacing` is the smaller apex spacing (a ternary `<`, as the
  source writes it); a neighbour is present when its spacing is below
  `spacing_difference * min_spacing`; both must be, or one with
  `allow_missing_flank`. A zero `spacing_difference` or `spacing_difference_gap`
  is infinity; both infinite disable spacing checks. Without spacing checks
  `min_spacing` stays zero and every spacing test passes, as for chromatograms.
- **Satellites.** A core flanked by more intense samples at `i - 2` and `i + 2`
  is skipped together with the next sample.
- **Extension.** Each side proceeds while the next intensity does not exceed the
  intensity stored at the current end of the support map, the spacing to that end
  is below `spacing_difference_gap * min_spacing`, the previous sample was not
  zero and at most `missing` samples failed. A failing sample is still added
  while the missing count allows it, and it always moves the boundary. Picking
  resumes after the last sample the right extension visited.
- **Support map.** The source keeps the support in a `std::map<double, double>`.
  The port keeps it in two reusable vectors with the same observable semantics:
  keys in position order, `map[key] = value` overwriting the value of an equal
  key, the extension comparing against the first and last stored entries, and
  the spline built from the sorted unique keys. At least three keys are needed.
- **Maximum.** `spline_bisection(spline, left, right, 1e-6)` with the neighbours,
  or the core where a neighbour is missing; the literal `(l + r) / 2` port. For
  finite, non-subnormal brackets this equals the former native
  `left / 2 + right / 2`, so the change is textual fidelity and moved no
  oracle value.
- **FWHM.** Half height is `max / 2` with tolerance `0.01 * (max / 2)`. A support
  end above half height is the crossing. The left bisection midpoint is
  `left / 2 + center / 2` and the right one `(right + center) / 2`, as the source
  writes each side; the loop stops on `!(|value - half| > tolerance)`. The width
  in ppm is `width / position * 1e6`.
- **Ion mobility.** `MSSpectrum::contains_im_data` and `im_data` pick the array:
  the first float array whose name is a PSI-MS child of `MS:1002893` or carries
  one of the vendor prefixes. The weight of each added sample is the float32
  product `mobility * intensity`, widened and summed in insertion order (core,
  left neighbour, right neighbour, left extension, right extension), divided by
  the sum of the stored intensities in key order, and narrowed to float32.
  With duplicate positions the weights keep every added sample while the total
  only the surviving values, as in the source. Three oracle cases distinguish
  this from double-precision products.
- **Narrowing.** Positions stay double; intensities, FWHM and mobility are
  narrowed to float32.
- **Output metadata.** A picked spectrum keeps every field of the input except
  its peaks and arrays and is marked centroided (`copySpectrumMeta`); a picked
  chromatogram keeps settings, metadata and name.
- **Experiments.** Automatic mode copies a spectrum whose `get_type(true)` is
  centroid (stored type, then a `PeakPicking` processing record, then the
  `PeakTypeEstimator` heuristic) and picks the rest. Manual mode copies unlisted
  MS levels and, with `check_spectrum_type`, refuses a listed centroid spectrum
  with `CENTROIDED_INPUT_MESSAGE` before picking it. Every chromatogram is picked
  without spacing checks.
- **Noise estimation.** The picker estimates each record it picks, only when
  `signal_to_noise > 0`, and compares the ratio of each apex against that
  threshold. The estimator's conventions (range, windows, bins and the CPP-257
  conversion, the median walk and interpolation, the percentages, warnings and
  the negative-range early return) are in
  [SIGNAL_TO_NOISE_SUPPORT.md](SIGNAL_TO_NOISE_SUPPORT.md#preserved-source-conventions).
- **Serial.** `PeakPickerHiRes.cpp` and `SignalToNoiseEstimatorMedian.h` have no
  OpenMP; the port is serial and its result does not depend on the thread count.

## Native differences

1. **Input refusals by default, source behaviour on request.** The native default
   refuses duplicate positions, decreasing positions, negative intensities, a
   non-positive spline maximum, a ppm width at a non-positive position and a
   second ion mobility array. Each `PickingCompatibility` flag selects the
   source behaviour for one of them and `PickingCompatibility::source()` for
   all; the oracle's `source_*` cases pin both sides. The source documents
   sorted input as a precondition without checking it.
2. **Ion mobility output description.** The native default copies the input
   array's metadata and processing handles onto the weighted mobility array
   (`MOBILOGRAM_SUPPORT.md` records this as a retention correction); the source,
   and `PickingCompatibility::source_mobility_arrays`, set only the name.
3. **Explicit mobility array.** `ion_mobility_array` selects a float array by
   exact name, a native extension.
4. **Non-termination.** The source's FWHM bisection never terminates once its
   midpoint equals a bracket end without meeting the tolerance, which a
   non-positive maximum always reaches. The port detects the fixed point and
   returns `Error::InvalidValue`; a defensive ceiling of
   `MAX_BISECTION_STEPS` (4096) halvings is never reached by a search that
   terminates in the source.
5. **Percentile range.** `auto_mode = 1` is computed exactly on the input
   domain where the source is defined (a non-empty record whose intensities
   lie roughly in `[m, m + 1)` for a minimum `m > 100 / 101`), and estimation
   returns `Error::Unsupported` everywhere else, where the source reads or
   writes out of bounds (CPP-256; the derivation is in
   `SIGNAL_TO_NOISE_SUPPORT.md`). A picker with `signal_to_noise = 0` never
   estimates, as in the source. Real spectra are outside that domain, the
   orbitrap class-test spectrum included.
6. **Estimator options and values.** As in the source, estimator options are
   only checked when estimation runs; a manual range with `max_intensity <= 0`
   is `Error::InvalidValue` there, the source's `Exception::InvalidValue`. The
   values the native profile refuses and `PickingCompatibility::source()`
   accepts (an infinite or NaN `win_len`, any `noise_for_empty_window`, a NaN
   `auto_max_stdev_factor`, non-finite results) are listed in
   `SIGNAL_TO_NOISE_SUPPORT.md`, native differences 2 to 5.
7. **Finite values.** Non-finite positions, intensities or mobility values are
   refused; an output that overflows float32 (the source stores infinity) or a
   non-finite mobility quotient is an error.
8. **Logging.** The source logs picked/total spectra per MS level; the port
   returns the boundaries instead. The estimator's three warnings, including
   the ungated negative-range one, are returned in `NoiseEstimates::log`, and
   `write_log_messages` gates the other two as in the source.
9. **Resource limits.** One million points per record, ten million core and
   extension visits, one million histogram bins and fifty million histogram
   updates, all configurable and checked before use. The source has none.
   The experiment entry points additionally charge every record's acquisition
   metadata to the shared `AcquisitionCopies` ledger. That ledger's fixed part
   is a ceiling on the
   *number* of records rather than on any one record: an ordinary
   vendor-converted run spends about 8 KiB of it per spectrum, so the fixed part
   alone stops at roughly 34 000 spectra — reached by the 2.3 GB Q Exactive run
   of the benchmark set, whose 40 856 spectra it refused outright.
   `max_metadata_per_record` (default 64 KiB) is added to the fixed part once
   per input record, so the budget follows the input. The allowance is pooled,
   as the shared ledger is, so one record may spend another's share; a
   single-record experiment whose metadata dwarfs it is still refused. Zero pins
   the fixed part, which is the behaviour before the field existed. Like
   `max_points` and `max_work` the field is Rust-API-only: it is not in
   `defaults()`, so `to_param` does not emit it and `from_param` cannot set it,
   and a caller driving the picker from a TOPP `.ini` always gets the default.
   That matters more here than for the other two, because this is the documented
   way out of a refusal that real data can provoke; whether the tool should
   expose it is a question for the `PeakPickerHiRes` tool lane.
10. **Boundaries per spectrum.** `spectrum_boundaries` has an entry for every
    input spectrum (`None` for a copied one).
11. **Experiment copy.** `pick_experiment` builds its output record by record,
    as the source does, so it never holds a second copy of the profile data. It
    does own a copy of every record it does not pick, which is what returning an
    owned experiment from a borrowed one means; `pick_experiment_in_place` is
    the streaming entry point that avoids even that, at the cost of atomicity
    (see the benchmark notes). A picked record is built from the input's
    metadata field by field, which is source `copySpectrumMeta`
    (`SpectrumHelper.cpp:14-25`); the port used to clone the whole input record
    and overwrite its samples one line later, copying every profile sample and
    annotation array of the record only to drop them. The fields are listed
    exhaustively, with no `..` rest, so a member added to `MSSpectrum` or
    `MSChromatogram` is a compile error there rather than a silently dropped
    one, and `the_picked_record_carries_every_metadata_field_the_input_had`
    pins the equivalence against the recipe it replaced.
12. **`estimate_spectrum_type`.** The public helper keeps the picker's strict
    input contract; `MSSpectrum::get_type(true)`, which `pick_experiment` uses,
    classifies any finite data as the source does.
13. **Chromatogram times come from the reader, not from here.** Nothing in the
    picker treats a chromatogram's retention times differently from a
    spectrum's m/z. The one place a picked chromatogram can differ from the
    source's while every centroid of every spectrum agrees is the reader's
    conversion of a 32-bit `time array` in minutes, which the source narrows to
    `f32` and this port keeps in `f64` unless
    `mzml::ReadOptions::source_time_array_precision` is set; see
    [the mzML reader's section](MZML_SUPPORT.md#the-minute-conversion-of-a-32-bit-time-array)
    for the mechanism and the measurement. The amplification is the picker's:
    the apex comes from a cubic spline through the support points and a
    bisection of its first derivative, so a relative change of up to 6.2e-3 in
    the point spacings moved the picked apexes of the benchmark TIC by up to
    3.19e-3 s and their intensities by up to 1.75e-3 relative.
14. **A parallel spectrum loop, where the source has none.**
    `PeakPickerHiRes.cpp` carries no `#pragma omp` at all: its spectrum loop
    (`:504`), chromatogram loop (`:548`) and on-disc loop (`:584`) are plain
    `for` statements and the header describes consecutive scans. (What OpenMS
    does parallelise for this workload is the mzML reader,
    `MzMLHandler.cpp:206`.) Picking one record reads nothing but that record, so
    `pick_experiment_with_threads` and `pick_experiment_in_place_with_threads`
    run the numerical half of the spectrum loop on a worker pool behind the
    `parallel` feature. `pick_experiment` and `pick_experiment_in_place` are
    those at `Threads::serial()`, so no existing caller changes behaviour.

    The determinism contract of `src/concept/parallel.rs` holds by
    construction, not by tolerance. The arithmetic is per record and never
    crosses records, so nothing is re-associated; results come back from an
    indexed parallel iterator in input order; the pooled acquisition-metadata
    ledger is charged **serially in input order** after the parallel pass,
    because a ledger charged from several workers would exhaust at a
    schedule-dependent record; and the first error is the first in input order,
    because the parallel pass collects a `Result` per record and the serial pass
    resolves them by index — `Result`'s own `FromParallelIterator` explicitly
    does not promise which of several errors it returns. The in-place form keeps
    its documented failure detail for the same reason: a batch is prepared in
    parallel and committed in input order, so exactly the records before the
    first failing one are replaced, where a parallel write-back would leave a
    schedule-dependent subset centroided.

    The parallel pass runs in bounded batches, at most
    `PARALLEL_BATCH_RECORDS` (4096) records or `PARALLEL_BATCH_POINTS`
    (16,000,000) input samples, whichever comes first, and always at least one
    record. That is a native memory bound with no counterpart in the serial
    source loop: it is what keeps the extra memory of a parallel pick at one
    batch of centroids rather than at the whole run, and independent of the
    worker count. Only the spectrum loop is parallel; chromatograms are picked
    serially after it, because the runs this is measured on carry a few thousand
    chromatogram points against hundreds of millions of profile samples.

    The workers are opened **once per call**, not once per batch (`BatchWorkers`
    in `src/processing/peak_picking.rs`). A caller already running inside a
    rayon pool — a TOPP tool is — uses that pool and builds nothing; a library
    caller outside one builds a single pool of the requested width for the whole
    call, where the first version of this code built one per batch and so paid
    the thread starts a dozen times over on an instrument-scale run. At one
    worker, or for an experiment of fewer than two spectra, no pool is built and
    the loop runs on the calling thread.
15. **One validation pass in the experiment path.** `start_experiment` validates
    the whole experiment through `MSExperiment::validate`, which visits every
    spectrum and every one of its peaks. Selected records used to be validated a
    second time on the way into `pick_`, over the same, unchanged samples; that
    pass is gone. Nothing mutates between the two — `pick_experiment` never
    mutates the input, `pick_experiment_in_place` only replaces records it has
    already picked — so an invalid record is refused with the same error, one
    pass over the samples earlier, and
    `an_invalid_record_is_refused_by_the_experiment_entry_points` pins both the
    picked and the merely copied case. The single-record entry points still
    validate what they are given.

## Class-test accounting

`PeakPickerHiRes_test.cpp` (16 sections):

| # | Section | Rust |
|---|---|---|
| 1 | `PeakPickerHiRes()` | `parameters_round_trip_through_the_typed_members` (default equals `from_param(defaults)`) |
| 2 | `~PeakPickerHiRes()` | not testable; `Drop` |
| 3 | `pick(input, output)`: dummy spectrum, dummy ion mobility with `Ion Mobility` and `raw inverse reduced ion mobility array`, orbitrap S/N 1 first spectrum | `upstream_dummy_peak_and_weighted_mobility`, `upstream_real_data_centroids_match_source_fixtures`; oracle `class_dummy`, `class_dummy_im`, `class_dummy_im_raw`, `class_orbitrap_sn1_spectrum0` |
| 4 | `pick(input, output, boundaries)`: orbitrap S/N 1, boundary literals 25 and 26 | `class_test_boundary_literals`, `source_boundaries_and_missing_flank_mobility`; oracle `class_orbitrap_sn1_spectrum0` |
| 5, 6 | `[EXTRA] pickExperiment` overloads | `NOT_TESTABLE` in the source; covered by 7, 9, 12 |
| 7 | `pickExperiment` orbitrap S/N 1 against `PeakPickerHiRes_orbitrap_sn1_out.mzML` | `class_test_pick_experiment_sections_match_the_retained_outputs`; oracle `class_orbitrap_sn1_experiment` |
| 8 | `[EXTRA] pick` orbitrap S/N 4 | same test (first spectrum); oracle `class_orbitrap_sn4_spectrum0` |
| 9 | `[EXTRA] pickExperiment` orbitrap S/N 4 | same test; oracle `class_orbitrap_sn4_experiment` |
| 10 | `[EXTRA] pick` FTMS S/N 1 | same test (and the whole experiment); oracle `class_ftms_sn1_spectrum0`, `extra_ftms_sn1_experiment` |
| 11 | `[EXTRA] pick` FTMS S/N 4 | same test; oracle `class_ftms_sn4_spectrum0` |
| 12 | `[EXTRA] pickExperiment` FTMS S/N 4 | same test; oracle `class_ftms_sn4_experiment` |
| 13 | `[EXTRA]` spectrum level selection, `ms_levels` 2, 1 and 1,2 | `class_test_spectrum_level_selection`; oracle `class_selection_ms2`, `class_selection_ms1`, `class_selection_ms1_ms2` |
| 14 | boundaries on simulation data: 167 peaks, literals 146, 148, 158, 159 | `class_test_boundary_literals`; oracle `class_simulation_boundaries` |
| 15 | boundaries on orbitrap data: 82 peaks, literals 14, 37, 54, 55 and the shared boundary | `class_test_boundary_literals`; oracle `class_orbitrap_boundaries` |
| 16 | `[EXTRA] allow_missing_flank`: symmetric, missing left, missing right, ion mobility | `source_boundaries_and_missing_flank_mobility`; oracle `class_flank_*` |

`SignalToNoiseEstimatorMedian_test.cpp` (5 sections) and
`SignalToNoiseEstimator_test.cpp` (6 sections) are accounted for in
`SIGNAL_TO_NOISE_SUPPORT.md`; the `[EXTRA] init` section is
`class_test_noise_estimator_section` here (through the parameters) and
`upstream_noise_values_match_source_fixture`, and oracle `class_noise_init` bit
for bit.

Two class-test inputs need care. `PeakPickerHiRes_spectrum_selection.mzML` holds
three MS2 spectra (`scan=5537`, `5541`, `5544`) with one decreasing m/z step
each: `MzMLFile::load` sorts them (`PeakFileOptions` sorts by default,
`MzMLHandler.cpp:218`), and the tests load with the same default through
`mzml::read_with_load_options`; loading unsorted changes centroid 63 of
`scan=5537`. `PeakPickerHiRes_ftms.mzML` and its two outputs repeat the id
`spectrum=1`, which the native mzML reader refuses; the committed copies rename
the second id and change no other byte.

The core data directory also holds `PeakPickerHiRes_{orbitrap,ftms}_sn0_out.mzML`
and four `*_ppmax.mzML` files that no class test references. The current C++
gives 82/112/89 and 314/319 centroids at `signal_to_noise 0`, not the stored
679/860/640 and 9359/9384; they are not used.

## Evidence

- **Tier 1, executed.** `tests/peak_picking_experiment.rs` compares 93 cases of
  `tests/data/peak_picking/cases.tsv` with the unmodified product SDK
  (`../oracle/peak-picker-hires/`), which is the **macOS arm64 Debug** build of
  core `4fdec46` (its diff to `bc9cc12` is empty for these sources), run twice
  byte-identically: 8,795 centroids with their positions, intensities,
  boundaries, FWHM and ion mobility values, the copied-record equality flags and
  the four exceptions, all as IEEE-754 bit patterns; the comparison itself runs
  on every platform the crate's tests run on. Every non-`source_` case runs in
  both the native default and `PickingCompatibility::source()`; every `source_`
  case must be refused by the default. Inputs: the class-test mzML and DTA
  files, the TOPP workflow 1, 2 and 6 inputs with their parameter values, and
  synthetic records for flanks, spacing, FWHM units, satellites, ion mobility
  names and float32 product rounding, histories, marked centroids, negative,
  duplicate and unsorted data.

  The noise estimator's share of this oracle is smaller than the case count
  suggests: 11 cases carry estimator ratios and percentages directly (the 12
  `noise` cases less `extra_noise_manual_invalid`, which throws), and their
  percentages are only ever `0` or `100`; 19 picker cases run the estimator
  indirectly (`signal_to_noise > 0`), among them 6 chromatogram records (5 in
  `extra_topp2_parameters`, 1 in `extra_topp2_chromatogram0_check`), where
  only the picked peaks are compared; the other 62 cases never estimate,
  because `signal_to_noise` defaults to `0`. The estimator's direct Linux
  x86-64 Release evidence is in `SIGNAL_TO_NOISE_SUPPORT.md`.
- **Tier 1, the same cases on Linux x86-64 Release.** The unchanged P1 driver
  (sha256 `2d06db2f…`) re-run against `openms4-release-bc9cc12-c19e494-174b576`
  on `ibminode06` (`../oracle/sne-completion/p1/`, twice, byte-identically)
  printed the same 93 cases byte for byte (sha256 `9eb8f249…`) and the same two
  parameter files; its standard error differs only by six mzML loader lines
  that the Debug build alone prints. The fixture is therefore also the Release
  build's output. Every case also runs with the clamp-first bin conversion and
  gives the same bits, and
  `no_executed_case_bins_a_quotient_outside_the_int_range` shows that no input
  record of an estimating case reaches CPP-257 at all.
- **Tier 1, executed, chromatogram time units.**
  `tests/data/peak_picking/chromatogram_time_oracle.tsv` is the output of the
  prebuilt C++ **Release** build `openms4-release-bc9cc12-c19e494-174b576` over
  the seven cases of `tests/data/peak_picking/chromatogram_time/` (driver,
  cases and hashes in `../oracle/picked-chromatogram/`, run on `dax`). Each case
  holds the same 40 real time and intensity values of the benchmark TIC
  chromatogram and differs only in the time array's encoding (32-bit or 64-bit)
  and unit (minute or second), plus the empty, single-point and unsorted edge
  shapes. Both the values `MzMLFile::load` stores and the chromatogram
  `pickExperiment` produces are compared as IEEE-754 bit patterns. The same TSV
  pins **both** modes of `mzml::ReadOptions::source_time_array_precision`: the
  `min32` rows are what the source stores for a converted 32-bit array, and the
  `min64` rows are the C++'s own answer for the full-precision seconds, which is
  what the library default must reproduce from `min32`. On the whole 2.3 GB
  input the two tools' picked TIC chromatograms are bit-identical with the
  switch set; the measurement is in
  [the mzML reader's section](MZML_SUPPORT.md#the-minute-conversion-of-a-32-bit-time-array).
- **Tier 1, parameters.** `defaults.ini` and `noise_defaults.ini` from the
  product SDK, and the retained `WRITE_INI_OUT.ini`.
- **Tier 1, retained class-test outputs.** The orbitrap and FTMS S/N 1 and S/N 4
  outputs under `TEST_REAL_SIMILAR` (`ClassTest.cpp` defaults: absolute
  difference `1e-5` or ratio `1 + 1e-5`), the simulation and orbitrap boundary
  literals, and the historical noise output.
- **Oracle driver note.** libOpenMS is compiled with `-ffp-contract=off`; the
  driver uses the same flag, because the header-only estimator it instantiates
  otherwise fuses some `a * b + c` into FMA on arm64 and moved six noise cases
  by one unit in the last place. Which expression fused was not recorded; the
  candidates are the range `mean + sqrt(variance) * factor` and the
  interpolation `lower_edge + ratio * bin_size`. A Release C++ build for
  benchmarking needs the same contraction policy to be comparable.
- **Tier 4.** Native refusals, each compatibility flag in isolation, resource
  limits and atomicity (`tests/peak_picking.rs`). The acquisition-copy ledger
  and the streaming entry point are covered by two files.
  `src/processing/peak_picking.rs` holds the synthetic ledger probe:
  `the_acquisition_ledger_follows_the_record_count` reads the ledger
  `acquisition_ledger` opens and requires it to equal the shared fixed part plus
  `max_metadata_per_record` per record in both work and bytes, to be strictly
  increasing in the record count, and to collapse to exactly the fixed part when
  the allowance is zero; `an_overflowing_metadata_allowance_is_refused` pins the
  two checked multiplications. It probes the arithmetic rather than building an
  experiment that exhausts the ledger because the meter charges a record within
  about a factor of two of what that record costs in memory, so crossing the
  256 MiB fixed part end to end costs the test process on the order of a hundred
  mebibytes, and four times the crossing point costs of the order of a gibibyte.
  That the ledger now admits a whole run is evidence, not a test: the real
  40 856-spectrum file is measured in the benchmark notes below.
  `tests/peak_picking_experiment.rs` covers the end-to-end behaviour that is
  cheap: `an_overflowing_metadata_allowance_is_refused_rather_than_wrapping`
  drives both entry points through the overflow refusal and requires the
  streaming one to leave its input untouched,
  `a_record_whose_metadata_outweighs_the_input_is_still_refused` keeps the
  adversarial single record refused -- the one test here that does reach the
  fixed part, and it builds its record by moving rather than cloning it, which
  is what keeps it to 95 MiB and 0.31 s (`/usr/bin/time -v` on the debug binary
  on `dax`, `--test-threads 1`) -- and
  `in_place_picking_matches_the_borrowing_entry_point` requires
  `pick_experiment_in_place` to produce the same experiment and the same four
  report vectors as `pick_experiment` over six class-test inputs in three MS
  level modes. `examples/peak_picking_scale.rs` is the harness behind the
  benchmark notes; it reports stage wall times, peak RSS and a bitwise digest of
  every centroid, and its `--ledger-probe` mode measures the ledger ceiling for
  a given run's metadata. Its body is behind the `mzml` feature, because
  `cargo test` builds every example and the crate must still compile with no
  default features.

  What the ledger tests cost, measured with `/usr/bin/time -v` on the debug
  binaries on `dax`: the synthetic probe 3 MiB and under 10 ms, the overflow
  test 3 MiB and under 10 ms, and the whole `peak_picking_experiment` target
  113 MiB and 0.86 s, almost all of it the single-record refusal. An earlier
  form of the probe drove the property end to end -- it searched for the record
  count the fixed part refuses and then built an experiment four times that
  size -- and cost 835 MiB and 2.68 s for the target, against 36 MiB and 0.95 s
  for the target before this work. That is the cost the factor-of-two metering
  implies, and the reason the property is probed rather than provoked.

`tests/data/peak_picking_provenance.json` records the pinned sources, every
fixture hash with its origin, the oracle artifacts outside the repository and
the earlier first-spectrum TSV fixtures, which are unchanged.

## Benchmark notes

The pick loop allocates nothing per peak besides what the source allocates:
the support map lives in two vectors reused for every peak of a record, and
positions and intensities are read in place. Remaining per-peak and
per-record costs, reported and not changed here:

- `CubicSpline2d::with_max_points` (read-only module) allocates **eight** vectors
  per peak — `h`, `mu`, `z`, `b`, `c`, `d`, `x.to_vec()` and `y[..n].to_vec()` —
  and checks every intermediate for finiteness; the source allocates its map
  nodes and five spline vectors. On the benchmark run that is 22,776,198 peaks
  × 8 allocate/free pairs, and `alloc::alloc` is 271,144,136 of the
  590,918,793 instructions the construction costs, i.e. 46% of it. This is the
  largest cost the picker still pays that the source does not, and it is not in
  this module: `src/processing/spline/cubic.rs` is shared with the retention-time
  transformations. A `fit_into(&mut self, x, y, max_points)` that reuses the
  spline's own vectors plus a small scratch, with the recurrence unchanged term
  by term, would be bit-identical and remove all eight; integrator request.
- `pick_experiment` builds its output record by record, as source
  `pickExperiment` does; it no longer clones the whole input experiment first,
  and since `perf/peak-picker` it no longer clones each picked record either
  (native difference 11).
- `MSSpectrum::get_type(true)` copies positions and intensities of each
  unknown-type spectrum to estimate its type. The source's `pickExperiment`
  makes the same query, so this is not a port-only cost; it is a second full
  pass over the samples of every record whose type is not stored.
- The noise estimator scans the input once more for its own checks. It is off
  on the benchmark INI (`signal_to_noise` 0 disables estimation entirely) and
  so is not on that hot path at all.
- `validate_points` still scans a record inside `pick_`. Under
  `PickingCompatibility::source()`, which every tool path uses, that is one pass
  (the order and duplicate scans are skipped); under the native default it is
  three, which cannot be fused without changing which refusal wins.

Measured on `ibminode06` against the 2.3 GB, 40 856-spectrum Q Exactive run
`profile_hr_qe_silac_uk222/UK222.mzML`, with `PickingCompatibility::source()`
and the C++ Release build `openms4-release-bc9cc12-c19e494-174b576` as the
reference. Load and pick only: when these runs were made the Rust mzML writer
had its own fixed budget and refused an output this size, so no Rust end-to-end
number exists here. The writer has since been given per-record budgets
(`fix/mzml-writer-scale-parity`, merged into `integrate/wave2`), so the
load + pick + write number that would face the C++ tool's 3883 MiB and 28.4 s
is now measurable and has not been measured;
`examples/peak_picking_scale.rs --in-place --out` is the harness for it.

| | peak RSS | wall (one run) |
| --- | --- | --- |
| C++ `PeakPickerHiRes`, load + pick + write | 3883 MiB | 28.4 s |
| C++ `FileInfo`, load only | 3201 MiB | 12.1 s |
| Rust before, 40 856 spectra | — | refused by the ledger |
| Rust before, first 30 000 spectra | 4268 MiB | 29.8 s |
| Rust after, first 30 000 spectra | 3409 MiB | 28.2 s |
| Rust after, 40 856 spectra, `pick_experiment` | 4309 MiB | 35.1 s |
| Rust after, 40 856 spectra, `pick_experiment_in_place` | 3410 MiB | 34.4 s |

On the 30 000-spectrum subset, which is the largest the old ledger admits, peak
RSS falls by 20% and the picked output is bit for bit what it was. After the
change the peak of the whole run is the peak of loading it: picking adds nothing
to the high-water mark under `pick_experiment_in_place`.

**How far these numbers carry.** Every peak RSS above was reproduced
independently on the same node, to better than 0.2%; those figures are the
result. The wall column is single-run and is not: a second, independent set of
runs gave C++ `PeakPickerHiRes` 31.4 s against 28.4 s here, C++ `FileInfo`
12.9 s against 12.1 s, and on the 30 000-spectrum subset a total that moved the
other way, because the load stage alone differed by 2.9 s between the two
builds' runs. No claim of a whole-run wall improvement is made, and none should
be read out of the table. The wall effect that does reproduce is in the pick
stage with the load stage subtracted: 11.8 s to 10.0 s, about -15%, which is the
whole-experiment clone no longer being made. The ordering of `pick_experiment`
and `pick_experiment_in_place` on the full run also flipped between the two sets
of runs, so the two are wall-indistinguishable here and differ only in memory.

### The parallel spectrum loop, end to end (`perf/peak-picker`)

The end-to-end numbers now exist and live with the tool that produced them:
`TOPP_PEAK_PICKER_HI_RES_SUPPORT.md`, *Instrument scale across thread counts*.
Same node, same 2.3 GB input, same C++ Release build, medians of three pinned
runs with `-test`. In short:

* **Bit-identical.** The `PeakPickerHiRes` tool wrote the same 535,613,726
  bytes (sha256 `bb13eecf…`) at 1, 8 and 32 workers, with the `parallel` feature
  off, and from the `fabd4b9` build that predates this branch. The determinism
  contract holds on real data at instrument scale, and the single-thread changes
  here moved no byte.
* **1.66x at 32 workers** on the whole tool (37.32 → 24.40 → 22.53 s), each run
  pinned to as many CPUs as it was given workers. The ceiling is the serial
  mzML read, about 45% of the run.
* **The C++ tool's own `-threads`, without the CPU-count confound.** Its picker
  has no OpenMP — `PeakPickerHiRes.cpp` at `bc9cc12` carries no `#pragma omp` at
  all, verified by grep, and its three loops are plain `for` statements — but
  that does not make the tool flat in the setting. Pinned so that every run has
  the same eight CPUs and only `-threads` changes, medians of three, the C++
  tool goes **27.39 s at one worker to 25.45 s at eight**: 1.94 s, 7.1%, with no
  extra CPU to spin on. Which region that comes from is not separated by this
  benchmark and is not claimed here; the one OpenMP region the run passes
  through is the mzML loader (`MzMLHandler.cpp:205-206`). On those same eight
  CPUs the port's parallel picking is worth more than that whole response:
  36.76 → 23.94 s.
* **Peak RSS flat in the worker count**: 3,527,720 KiB at one worker,
  3,519,504 KiB at 32, and 4,414,500 KiB before this branch. The parallel path
  holds one bounded batch beyond a serial pick, not one per worker.
* **The two single-thread changes** — the metadata-only record construction
  (native difference 11) and the removed second validation pass (15) — are worth
  **-1.84% instructions, -5.05% data references and -9.27% L1 data misses** over
  a 682-spectrum slice in the `parallel`-off build, and **-1.71%, -4.85% and
  -9.18%** in the configuration an ordinary run uses, the shipped default at
  `-threads 1`. Callgrind, because wall time on that node could not resolve them
  under its foreign load. The default keeps the saving only because the pool is
  scoped to the picking call: with the pool around the whole tool body the same
  build gave it back in full (-0.02% instructions). Fewer data references than
  instructions is what removing copies rather than computation looks like.

The library's own harness for the same shape is
`examples/peak_picking_scale.rs --in-place --out`.
