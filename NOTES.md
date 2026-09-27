# Lane F4 notes (delete in the final commit)

Branch: worktree-wf_ee239125-5be-14 (second incarnation; merged worktree-wf_bcdd783d-f3f-14 ff), base f0990c0.

Scratchpad (session-shared): /private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/F4/
- cpp/<Name>.{h,cpp,test} : dumped C++ at core bc9cc12
- data/ : class-test data (they live in CORE src/tests/class_tests/openms/data at bc9cc12, not test-data)
- g.py : `python3 g.py <pkg> show|ls|grep|fetch|sha ...` wrapper for pinned git show
- gate.sh <log> <cargo args> : runs the dax gate from this worktree, log in scratchpad/F4/<log>.log
  (run with run_in_background; the worktree guard refuses complex bash with git/python+git)
- Edits: write python edit scripts to scratchpad and run them; do not mix heredoc python with git/cargo.

## Done (committed)
- EDTAFile -> src/format/edta.rs, tests/edta_file.rs (16 pass)
- ChromeleonFile -> src/format/chromeleon.rs, tests/chromeleon_file.rs (11 pass)
- shared: src/format/string_utils.rs (moved from cli/context.rs), src/format/line_reader.rs
- GNPSQuantificationFile, GNPSMetaValueFile -> gnps_quantification.rs, gnps_meta_value.rs, tests/gnps_tables.rs (8 pass)
- IBSpectraFile -> src/format/ibspectra.rs, tests/ibspectra_file.rs (7 pass)
- MSPGenericFile -> src/format/msp_generic.rs, tests/msp_generic_file.rs (9 pass); stod/stoul in string_utils; LineReader::lf_only
- CachedMzMLHandler + CachedMzML -> cached_mzml_handler.rs, cached_mzml.rs, tests/cached_mzml.rs (14 pass)
- SqMassFile -> sqmass_file.rs, tests/sqmass_file.rs (7 pass, --features sqmass); FIX numpress_coder zlib state release
  (gap: convertToXICParquet needs Arrow/Parquet + MSChromatogramParquetConsumer + LightTargetedExperiment)
- Gate commands: filter args like `format::` apply to ALL test targets -> run integration tests without filters.

## Findings for the report
- Rust mzML writer refuses a non-CV instrument name; C++ writes generic MS:1000031 and drops it
  (MzMLHandler.cpp:3961-3969). Chromeleon->mzML needs name cleared.
- GNPSExport_1_out_quant.txt: 52 fields on 44 lines spelled non-exactly (old formatter); same doubles.
- IBSpectra: hits[0] unchecked in C++ (UB on empty hits) -> refused here.

## Next step
1. MSPGenericFile, MSPFile (text readers; BoostRegex facade).
2. CachedMzMLHandler + CachedMzML (TICCalculator).
3. SqMassFile (wrap mzml_sqlite_handler).
4. GNPSMGFFile, FLASHDeconv files, SwathFile/SwathFileConsumer.
5. Support docs, provenance manifests, report F4.md, ledger wave4/F4.json.
