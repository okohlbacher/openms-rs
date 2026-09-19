# Iterative peak-picking reference review

This review covers `PeakPickerIterative` at OpenMS4-core revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. The implementation is in [PeakPickerIterative.h](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h); its `.cpp` contains no algorithm. The [source class test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/PeakPickerIterative_test.cpp) checks construction and destruction, with numerical sections marked TODO. There are no upstream iterative numerical goldens to claim. No C++ was built or executed.

The reference fixtures come from an independent Python translation of the header's arithmetic and control flow, exercised on explicit raw traces and seed lists. They are derived algorithm references, with separate manually interpretable boundary cases. [iterative_provenance.json](../tests/data/iterative_provenance.json) records the seven source hashes, four fixture hashes and the provenance of two reused real input spectra. The complete core port remains in progress.

## Derived numerical cases

The four TSV files separate [106 raw samples](../tests/data/iterative_traces.tsv), [16 seeds](../tests/data/iterative_seeds.tsv), [fourteen configurations](../tests/data/iterative_options.tsv) and [expected surviving candidates](../tests/data/iterative_expected.tsv). Intensities and reported centroids include binary32 bit patterns. Expected records preserve the original seed ordinal, initial and final raw centers, inclusive raw boundary indices, unrounded binary64 centroid and integrated intensity.

The source first calls `PeakPickerHiRes`, passing the iterative signal-to-noise and spacing settings. Refinement uses its separate median-noise estimator, with a default 20-unit window rather than the seed estimator's default 200-unit window. Explicit seeds are injected only into a private refinement helper for the independent tests; no test-only public API was introduced.

The cases pin these details:

- A seed exactly at 103 maps to the raw point at 104, because association uses the first coordinate strictly greater than the seed. With extension disabled by strict spacing equality, its three-point core has sum 15 and centroid 103.4. It is not assigned to the nearest raw sample.
- Two seeds between the same raw samples advance at most one seed per raw-loop iteration. They receive successive raw indices. The fixture also distinguishes original seed priority from the resulting integrated intensity.
- The central raw sample and its immediate neighbors always enter the core. Further points require strict spacing and either strictly falling intensity or strict inclusion in the configured half width. Equal widths and equal intensities therefore stop extension unless the other condition permits it.
- Minimum spacing stays fixed at the smaller of the two original core gaps for each recentering pass. A smaller outer gap does not tighten subsequent extension spacing.
- A seven-point symmetric peak has integrated intensity 20 and centroid 103. The reported intensity is the inclusive sample sum, not a trapezoidal area.
- The retained right-search condition is `i-m > 0` despite accessing `i+m`. In the asymmetric case, an initial center at raw index 2 moves only to index 3, although the centroid is closest to index 6. This finite source behavior is preserved, with an additional bounds check.
- A symmetric peak around `100000003.25` recenters at the correct raw sample using the unrounded centroid. The stored binary32 candidate rounds to `100000000`, which lies below the full-precision left boundary `100000000.25`. Output and suppression use that rounded candidate; raw recentering uses the unrounded value.
- Sparse seven-point noise windows use the source default noise value `1e20`. An extension intensity of four passes a threshold of `4/1e20`; the next larger representable threshold rejects it. Signal-to-noise equality is accepted.
- Internal width validation accepts equality and rejects spacing above the configured width. Invalidated candidates are omitted.
- In the overlapping-priority case, the seed with original intensity 100 and integrated intensity 43 suppresses a seed with original intensity one and integrated intensity 57. Sorting by the new integrated intensity would change the scientific result.

Source sorting does not specify the order of equal-priority seeds. A separate test asserts the native policy of retaining original seed order. That expected tie order is an intentional deterministic policy, not an upstream numerical assertion.

## Public behavior and checked boundaries

The [private review tests](../src/processing/iterative/review_tests.rs) compare all fourteen explicit-seed cases with the derived records and check deterministic ties. The four [public reference tests](../tests/iterative_picking_reference.rs) exercise the full HiRes-seeded picker, metadata and array handling, experiment selection, failures and real input spectra.

Selected spectra become centroid spectra. The three generated float arrays are `IntegratedIntensity`, `leftWidth` and `rightWidth`, aligned with final sorted peaks. Width arrays contain absolute boundary coordinates rounded to binary32. The typed boundaries retain binary64 coordinates, and typed regions retain original raw indices. Input float, integer and string annotations have no defined aggregation rule and are reported as omitted; ordinary spectrum metadata remains intact.

Experiment tests verify `ms1_only`, clearing generated float arrays only for selected spectra when `clear_meta_data` is enabled, preservation of unselected spectra, and preservation of chromatograms. The direct spectrum API retains generated arrays even when that experiment option is enabled. The native experiment result preserves chromatograms; the C++ experiment overload only constructs its spectrum output.

The review found no unbounded refinement or suppression loop: association, recentering, sample integration, recenter searches, sorting charges and candidate comparisons consume the configured work budget. The budget is per processing stage, not a single aggregate for all seed/noise/refinement stages or the whole experiment. Noise histogram setup is charged before construction, and the noise estimators retain their own allocation and visit bounds. Input and seed counts are checked before their main allocations. Excessive iteration products use checked arithmetic.

Atomic failure checks exercise invalid later spectra, point and work limits, iteration-count overflow and nonfinite configuration. The mutable spectrum and experiment filter interfaces leave inputs unchanged on failure. Native safety checks also reject negative profile m/z, duplicate raw coordinates, invalid seed neighbor access, nonpositive integrated intensity, nonfinite numerical results and values that cannot be represented in the output float arrays; negative intensities are rejected by the native profile only, since `PickingCompatibility::allow_negative_intensities` accepts them. Rejecting negative m/z avoids treating a legitimate negative left boundary as the source's invalid-candidate sentinel.

## Real spectra and limits of the evidence

The existing first-spectrum Orbitrap and FTMS input fixtures are reused directly, with their original source and fixture hashes linked in the provenance record. Their HiRes output fixtures are not iterative-picker references. Instead, the public tests independently reconstruct each iterative output intensity and rounded centroid from its reported inclusive raw region, verify region and peak alignment, and validate the resulting spectrum. The manually derived symmetric peak supplies a complete end-to-end numerical case.

The source retains several unusual choices, including strict-next-sample seed association, original-seed suppression priority, binary32 candidate m/z and asymmetric recenter search. These choices are covered explicitly rather than replaced with plausible alternative algorithms. The native tie policy, checked failure cases, short-input metadata preservation and typed annotations are documented boundaries. Broader scientific equivalence on additional instruments remains unverified; the absence of source numerical tests is not evidence of full parity.

The focused suite consists of two private tests and four public tests. Current-toolchain and Rust 1.85 checks cover the supported feature configurations used in the repository validation matrix. See [ITERATIVE_PICKING_SUPPORT.md](ITERATIVE_PICKING_SUPPORT.md) for the public API and supported behavior.
