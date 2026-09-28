# FF4 lane notes (Biosaur2)

## Scope
- Header: OpenMS/FEATUREFINDER/Biosaur2Algorithm.h (+ source/FEATUREFINDER/Biosaur2Algorithm.cpp) @ core bc9cc12
- Class test: src/tests/class_tests/openms/source/Biosaur2Algorithm_test.cpp @ core bc9cc12
- TOPP consumer: topp 174b576 src/FeatureFinderLFQ.cpp (222 lines, wraps Biosaur2Algorithm)
- TOPP test: test-data 0cb15f2 topp/CMakeLists.txt:429-431 TOPP_FeatureFinderLFQ_1
  (-in FeatureFinderCentroided_1_input.mzML -mini 20 -htol 100 -itol 100 -minlh 2 -cmin 2 -cmax 2;
   expected topp/FeatureFinderLFQ_1_output.featureXML, FuzzyDiff whitelist "id=")
- Also used by ProteomicsLFQ -Seeding:algorithm biosaur2 (not this lane)
- Rust target: src/analysis/featurefinder/biosaur2.rs (append `pub mod` in src/analysis/featurefinder.rs)

## Git access to pins
Worktree guard refuses `git -C`; use `cd <pkg> && git show <sha>:<path>` as separate simple commands.
Pins root: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages

## Next step (updated session 1)
1. DONE: port drafted in src/analysis/featurefinder/biosaur2{.rs,/hills.rs,/isotopes.rs,/preprocess.rs,/output.rs,/averagine_table.rs}
2. DONE: tests/biosaur2.rs + tests/data/biosaur2 (oracle fixtures staged by ../oracle/FF4-biosaur2/stage_fixtures.py)
3. NOW: get tests/biosaur2.rs green on dax (gate slot FF4, logs in scratchpad gate_*.log)
4. THEN: TOPP FeatureFinderLFQ tool (src/cli/tools/feature_finder_lfq.rs + src/bin/FeatureFinderLFQ.rs + Cargo [[bin]]),
   tier-1 vs retained FeatureFinderLFQ_1_output.featureXML and ../oracle/FF4-biosaur2/results/topp/*
5. THEN: docs/BIOSAUR2_SUPPORT.md, tests/data/biosaur2_provenance.json, doc coverage, clippy, MSRV, ci sweep, report+ledger.

## Findings (session 1)
- Sources extracted to scratchpad; header 722 lines, cpp 3373, test 446, TOPP 222.
- Rust building blocks that exist: crate::math::source_sort::source_sort_by (libstdc++ std::sort replica, use for
  the 3 std::sort calls: peaks by desc intensity in linkScanToHills_, hills by mz_weighted_mean, candidates by
  (n_isotopes desc, cos desc)); kernel::faims_helper::FaimsHelper::get_compensation_voltages;
  kernel::im_data_converter::ImDataConverter::split_by_faims_cv (returns FaimsSplit{groups,dropped_chromatograms,..});
  processing::feature_overlap_filter::FeatureOverlapFilter::merge_faims_features; metadata::ImTypes::determine_im_format;
  MSSpectrum::im_data()/contains_im_data(); MSExperiment::sort_spectra.
- Rust MSSpectrum has NO stored IM peak type (only SpectrumSettings has) -> IM_CENTROIDED skip logic is a gap to document.
- Peak1D.intensity is f32 (as C++ float).
- TOPP template: src/cli/tools/feature_finder_centroided.rs + src/bin/FeatureFinderCentroided.rs; [[bin]] in Cargo.toml.
- boost::math binomial pdf is used for the averagine table (n_C trials, p=0.0107, k=0..9).
- OpenMP: FAIMS groups parallel for; simd reduction in cosineCorrelation1D_ (check Release bits via oracle).

## Oracle (session 1)
- ../oracle/FF4-biosaur2/: gen_cases.py (47 run cases incl. class-test inputs, synthetic, mzML; 378 helper probes),
  driver.cpp (#define private public; bit-pattern records), remote_run.sh (ibminode06, Release, runs twice, also
  Release FeatureFinderLFQ TOPP + -write_ini), run.sh (stage+start detached), fetch.sh (when results/DONE).
- Remote dir ibminode06:/scratch/kohlbach/FF4-biosaur2. Started; fetch with `bash ../oracle/FF4-biosaur2/fetch.sh`.
- Experiment text format (.exp) documented in gen_cases.py docstring; Rust test must implement same reader.
