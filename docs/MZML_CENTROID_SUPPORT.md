# Spectrum type queries and mzML centroid inspection

Source pin: `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. This increment implements the complete `MSSpectrum::getType(bool)` and `MzMLFile::getCentroidInfo`/`SpecInfo` operations. It does not claim completion of the remaining `MSSpectrum` or `MzMLFile` headers. The existing public peak-picker estimator retains its stricter input contract.

## Native API and precedence

`MSSpectrum::get_type(query_data)` and `get_type_with_limits(query_data, SpectrumTypeQueryLimits)` borrow the spectrum and do not cache their answer. These methods and their limits are available without format features. As in `MSSpectrum.cpp:142–165`:

1. A stored Profile or Centroid type wins immediately, even over PeakPicking history.
2. Otherwise, any visited record-level `DataProcessing` action set containing PeakPicking produces Centroid.
3. Otherwise, `query_data=false` produces Unknown; `true` invokes the source shoulder estimator.

No unrelated metadata, auxiliary arrays, acquisition settings, processing software or timestamps are validated or copied. History handles retain their `Arc` identity. Fewer than five peaks produce Unknown before peak values are consumed, including a short record containing nonfinite values. With at least five peaks, only finite m/z and intensity values are required. Finite negative intensities, duplicate coordinates and unsorted coordinates retain source behavior; the getter does not sort.

The extracted scalar helper preserves the existing native implementation of `PeakTypeEstimator.h:43–156`. It considers up to five strict first maxima, stopping after more than half the total intensity is explained. Shoulder intensity ratios are strictly greater than 0.1, and the m/z window is strictly within one Th; final Profile evidence is computed in f32 and must be strictly greater than 0.75. Thus three Profile candidates out of four are Centroid, while four out of five are Profile. All-zero or nonpositive signals with at least five points take the source zero-evidence comparison and return Centroid. Finite extreme m/z values retain ordinary floating-point rounding of `mz ± 1`.

`processing::estimate_spectrum_type` still performs its original whole-record validation, nonnegative-intensity checks and strictly increasing-coordinate checks before calling that same helper. No picker input domain or numeric expression was widened by the extraction. The generic source template's other point types and its inactive historical alternative algorithm are not new public native APIs.

## File inspection

With the `mzml` feature, `format::mzml::centroid_info(path)` uses the source default quota of ten recognized spectra. `centroid_info_with_options(path, quota, &LoadOptions, &ReadOptions, CentroidInfoLimits)` returns an owned `BTreeMap<u32, SpecInfo>`. `SpecInfo` has the source's three zero-default counters: `count_centroided`, `count_profile` and `count_unknown`. `CentroidInfoLimits` is an alias for the same query-limit type.

The positive quota is global across all MS levels after scientific filtering. Unknown records are counted but do not consume it. The last recognized record is counted before the successful soft stop. Empty input, no delivered spectra, or `metadata_only=true` returns an empty map. The actual stored u32 MS level is retained.

Inspection opens the input once through the existing nonretaining consumer, suppressing both setup callbacks. A local copy forces `fill_data=true`; all other load options remain effective, including filters, sorting, metadata-only and chromatogram skipping. Plain/gzip/bzip2 input uses the existing content-sniffed path reader. Pool decoding remains observable: an invalid later record in the same pool may fail before the first callback even when the requested quota is one. A soft stop leaves the remaining pool and unread XML/compression tail uninspected. The operation does not promise to decode only the requested number of records.

The caller's options, their vector allocation and all external input models remain untouched on success or error. A partial result map is never returned after genuine failure. The existing reader's strict XML/scientific domain and complete-pool predecode policy remain in force; the standalone getter's broader finite scalar domain does not weaken mzML decoding.

## Checked limits and original-source corrections

Defaults are one million points per actual estimate, 50 million cumulative classification work units, and 256 MiB cumulative classification/result bytes. Actual estimation precharges both f64 workspaces and conservative validation/copy/sum/five-scan work before allocation. History charges only visited handles and bounded action lookups. Explicit stored types require no history/peak budget. The file operation shares these counters across all classifications, option-vector copying and result-map bookkeeping/storage, including sparse BTree roots. They are conservative logical allowances, not allocator- or RSS-exact measurements. Parser, header, binary, selection and consumer-administrative limits remain their existing separate operation allowances.

Consumed nonfinite scalars, arithmetic size/count overflow and resource exhaustion produce checked errors. No new native dependency or numerical backend is introduced. [OpenMS_CPP_ISSUES.md](../OpenMS_CPP_ISSUES.md) records the two existing defects corrected here:

- **CPP-018:** source file inspection restores its temporary fill-data option only on success. The native operation uses local options, preserving caller state on every path.
- **CPP-048:** source accepts a zero unsigned quota, whose decrement wraps on the first recognized spectrum. Native zero quota is `InvalidValue` before even requesting the path or opening input.

No additional original C++ defect was found in this increment.

## Evidence

`tests/spectrum_type.rs` retains all type-query assertions from `MSSpectrum_test.cpp:1382–1418`. The three unchanged original DTA fixtures preserve the Profile/Profile/Centroid expectations from `PeakTypeEstimator_test.cpp:40–56`, plus its four-point Unknown assertion. Their bytes match the pinned checkout exactly.

Independent hand-derived examples exercise the strict shoulder/window boundaries, first-maximum ties, the 3/4 versus 4/5 evidence threshold, finite signed/duplicate/unsorted shapes, nonfinite short-circuiting and exact resource boundaries. `tests/mzml_centroid.rs` supplies explicit record inputs and independently stated expected counts/order; the native writer constructs transport documents only, never expected scientific answers. It covers default quota ten, multiple levels, Unknown interleaving, history/metadata precedence, filtering, metadata-only, forced population, caller state on four failure classes, pre-path rejection, pool predecode, unread tails, compressed inputs and cumulative classification/map/history limits.

The 16 new public tests are accompanied by the existing picker, consumer and count suites. Exact source and fixture hashes are recorded in `tests/data/mzml_centroid_provenance.json`. No C++ build or Rust-generated numerical oracle is claimed.
