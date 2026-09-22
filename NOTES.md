# Lane NOTES: close KERNEL/MSExperiment.h (delete before final commit)

Scratch copies of pin sources: $SCRATCH/MSExperiment.{h,cpp}
(core bc9cc12: src/openms/include/OpenMS/KERNEL/MSExperiment.h, src/openms/source/KERNEL/MSExperiment.cpp)
Upstream test: core bc9cc12 src/tests/class_tests/openms/source/MSExperiment_test.cpp

## Missing (from ledger scope + docs/EXPERIMENT_MOBILITY_SUPPORT.md "Deferrals")
1. set2DData<add_mass_traces=true> -> src/kernel/experiment_2d.rs (h:265-289, ContainerAdd_<.,true> h:1368-1392)
   - meta num_of_masstraces; per i: masstrace_intensity_<i> (missing -> Precondition), mz + C13C12_MASSDIFF_U/charge*i, charge 0->1
   - constant: src/concept/constants.rs; meta names in src/analysis/feature_finding_metabo.rs
   - test section 11
2. getFirstProductSpectrum(ConstIterator)/(int) -> src/kernel.rs beside precursor_spectrum_index (cpp:997-1060)
   - end/last -> end; stop at lower level; level diff==1 & first precursor meta "spectrum_ref" == parent native id; multiple precursors -> LOG_WARN
   - test sections 42, 43
3. operator<< -> src/kernel.rs Display (cpp:1365-1388): "-- MSEXPERIMENT BEGIN --\n", ExperimentalSettings <<, each spectrum <<, each chrom <<, "-- MSEXPERIMENT END --\n"
   - test section 54; need Release run for exact text
4. "missing-value NaNs": check gap vs divergence (store_metadata_names path, h:1355-1364 NaN when meta missing)
5. Ledger also notes SummaryLimits max_work fixed 50M (open, larger than one tool) - report, not necessarily this lane

## Next step
1. Read MSExperiment_test.cpp sections 11, 42, 43, 54 at pin
