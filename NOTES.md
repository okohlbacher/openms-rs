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

## Findings (session 2)
- Scratch: /private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad
  (use literal paths; sandbox refuses $VAR in sed/git and python heredocs mentioning git; use Edit tool + plain git cmds)
- Item 4 NaNs: NOT a gap. experiment_2d.rs import() writes f32::NAN for a missing key. Represented.
- Sec 11 (test.cpp:308-395): fm f1(int 7.5,rt 8.5,mz 9.5), f2(17.5,18.5,19.5; NUM_OF_MASSTRACES=2, masstrace_intensity_0=11f, _1=12f, charge 2)
  -> get2DData(FeatureMap) == [f1, f2{int 11}, f2{int 12, mz + C13C12/2}], same charge? (f2_x copies f2 which has charge 0! fm2 copied BEFORE setCharge)
- Sec 42/43 (test.cpp:1211-1340): levels 1,2,1,2,2; native ids scan=1..5; precursor spectrum_ref scan=1/scan=3/scan=3.
  expect 0->1, 1->end, 2->3, 3->end, 4->end, end->end; all MS1 -> all end.
  Rust Precursor has typed `spectrum_reference: Option<String>` (kernel.rs:182) == C++ meta "spectrum_ref".
  One Rust fn first_product_spectrum_index(usize)->Option<usize> covers both overloads (len == end).
  >1 precursors -> OPENMS_LOG_WARN "Spectrum at index N has multiple precursors. Only the first precursor will be considered." (cpp:1029-1031)
- Sec 54 (test.cpp:1641-1667): substrings MSEXPERIMENT BEGIN, MSSPECTRUM BEGIN, MSCHROMATOGRAM BEGIN, 47.11, 10.77.
  MSChromatogram Display exists (src/kernel/chromatogram_merge.rs:585); ExperimentalSettings Display exists
  (src/metadata/experimental_settings.rs:200, only BEGIN/END lines -- verify vs C++);
  MSSpectrum Display MISSING (docs/SPECTRUM_MOBILITY_SUPPORT.md:125) -> add in src/kernel.

## Design decisions (session 2)
- Mass traces: only caller topp FileConverter.cpp:590 (`exp.set2DData<true>(fm)`); Rust cli has no counterpart (not mine).
  C++ h:1374-1391: `Size mts = getMetaValue` -> DataValue::operator unsigned long (INT only, >=0, else ConversionError; DataValue.cpp:554-567);
  intensity via operator float (EMPTY->error, INT->float, else data_.dou_ (union garbage for string: record as divergence), DataValue.cpp:480-491).
  Spectrum created on RT change BEFORE addData_ -> a feature with num_of_masstraces=0 still opens an (empty) spectrum. Keep that!
  Error order per item: sort check, then meta. Missing name -> Error::InvalidValue("Meta value 'masstrace_intensity_i' expected but not found in container.").
  Implement as separate import_mass_traces in experiment_2d.rs (groups (rt,count) + flat points), public set_2d_data_mass_traces(_with_limits)(&[Feature]).
- getFirstProductSpectrum: warn via `let _ = crate::concept::log_stream::log_message(LogLevel::Warn, format_args!(..))` (precedent src/format/mzml_header/read.rs:102).
  Test capture pattern: tests/mzml_header_leniency.rs:571-596 (Capture sink on thread-local Warn log).
- operator<<: Peak1D/ChromatogramPeak Display use Rust number formatting (ledger: recorded divergence, native_equivalent); C++ uses %g prec 6.
  %g helpers: crate::format::file_info::text_format::ostream_g(v,6) (kernel->format edge already recorded in docs/module-cycles.json, so no new edge).
  Plan: Display layout exact; decide whether `{:#}` gives byte-exact %g. Release oracle for both.

## Next step
1. Read C++ MSSpectrum operator<< + SpectrumSettings operator<< + ExperimentalSettings operator<< at pin; chromatogram Display style.
2. first_product_spectrum_index in src/kernel.rs + tests.
3. MSSpectrum Display + MSExperiment Display; Release oracle on ibminode06 for exact text.
4. set_2d_data mass traces from FeatureMap.
