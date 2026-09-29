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

## Next step
1. Read XICParquetFile/XIMParquetFile/XIPMParquetFile headers + tests; port files first (consumers depend on them).
2. Then consumers, then the three ArrowExports.
