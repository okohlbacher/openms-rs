# F1b notes (Arrow exports and XIC/XIM/XIPM)

Prereq merged: worktree-wf_cd066dab-0a7-5 (F1a) fast-forwarded, HEAD 8bb62562. No earlier F1b incarnation.

C++ at pin bc9cc12 readable in /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.reference/openms4-core-bc9cc12
(sandbox refuses complex git; use that checkout, cite with git show one file per call).

Sizes (cpp lines): ConsensusMapArrowExport 1733, MSExperimentArrowExport 1336, ProteinGroupArrowExport 1448 (+418 impl.h),
XICParquetFile 798, XIMParquetFile 851, XIPMParquetFile 348, XIPMParquetConsumer 771,
MSChromatogramParquetConsumer 748, MobilogramParquetConsumer 704.
Class tests: src/tests/class_tests/openms/source/<Name>_test.cpp; data XICParquetFile_{1,2}_input.xic,
XIMParquetFile_23_input.xim, XIPMParquetFile_1_input.xipm, MSChromatogramParquetConsumer_1_{input.TraML,output.chrom.mzML}.

F1a infra: src/format/{arrow_io_helpers,arrow_schema_registry,parquet_file}.rs, feature `arrow`.

## Findings
- Merged lane OA (worktree-wf_cd066dab-0a7-27) too: LightTargetedExperiment (src/analysis/openswath/transition_experiment.rs),
  openswath_helper::compute_transition_group_id. REPORT this extra merge.
- No TraML reader in Rust -> MSChromatogramParquetConsumer sections 1,2,4 (TraML input) not runnable; 3,5 constructed are.
- XIPMParquetConsumer needs PeakMapExtractor::ExtractedPeakMap (ANALYSIS/OPENSWATH/PeakMapExtractor.h) - unported.
- Private C++ helper src/openms/source/FORMAT/DATAACCESS/ParquetReaderHelpers.{h,cpp} (1262 lines) is shared by
  XIC/XIM/XIPM readers: port as crate-private src/format/parquet_reader_helpers.rs.
  Arrow compute filter == manual filter for monotone exprs (null->false), left-to-right no precedence.
- XIM reader = XIC reader + MOBILOGRAM_TYPE/FEATURE_ID/FEATURE_RT; pruning opts {MOBILITY,INTENSITY,FEATURE_RT},true,true;
  nest_transitions requires precursor discriminator. XIPM: no filter string; Subset schema validation; target_rt midpoint.
- Codec ids 0..7 (0 raw,1 zlib,2/3/4 np lin/slof/pic,5/6/7 +zlib); XIPM has no 7 ("Unknown compression code").
- Reuse numpress_coder::{zlib_encode,zlib_decode,decode_raw,encode_raw,Work} (pub(crate)); sqMass adapter
  (mzml_sqlite_handler.rs store_array) guards linear <3 values (fixed point 0 -> RT lost) -> lossless code 1.
  C++ defect candidate: XIC consumer writes 1-2 point chromatograms with fixed point 0 => RT decode to 0.
- Fixtures copied to tests/data/openswath_parquet/ (sha256 recorded via shasum). pyarrow available locally.

## Next step
1. Write src/format/parquet_reader_helpers.rs (crate-private) + xic_parquet_file.rs, tests/openswath_parquet.rs.
2. xim_parquet_file.rs, xipm_parquet_file.rs.
3. ms_chromatogram_parquet_consumer.rs, mobilogram_parquet_consumer.rs (XIPM consumer partial: needs ExtractedPeakMap).
4. ArrowExports last (likely untouched).
