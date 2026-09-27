# F1a NOTES (lane scratch, delete in final commit)

Lane: F1a — Arrow/Parquet core. Branch: worktree-wf_ee239125-5be-10 (fast-forwarded onto the earlier
incarnation worktree-wf_bcdd783d-f3f-10 at 27d502a; no merge commit was needed).

## Done
- 52953b5: Cargo.toml `arrow` feature, arrow-array/schema/buffer/select + parquet =59.3.0 (MSRV verified on
  dax with +1.85.0; 60.0.0 needs 1.88), Cargo.lock fetched back from dax (144 -> 180 packages),
  THIRD_PARTY_CRATE_DECISIONS row (fidelity/determinism sentences promise tests in tests/arrow_parquet.rs
  and docs/ARROW_PARQUET_SUPPORT.md -- keep that promise or edit the row).

## Facts found (expensive)
- C++ dumps in scratchpad F1a-w4/ (headers at src/openms/include/OpenMS/FORMAT/, sources at
  src/openms/source/FORMAT/). Use a script file run with `bash script.sh` for git show in other repos
  (the worktree guard refuses inline `git show` in other dirs).
- Gate: Cargo.lock is NOT synced back; fetch it with
  `ssh -J kohlbach@sshgw.cs.uni-tuebingen.de kohlbach@dax cat /scratch/kohlbach/agents/F1a/Cargo.lock`.
- ParquetWriteConfig lives in MSExperimentArrowExport.h (not my header); row_group_size is documented
  as BYTES but passed as ROWS to WriteTable; C++ caps chunk_size at max_row_group_length (1Mi rows).
- Retained C++-written fixtures (tier 1 for READ): test-data 0cb15f2 topp/FileConverter_9_output.featureparquet/,
  FileConverter_14_output.consensusparquet/, IDMerger_idparquet_input{1,2}.idparquet/ (5/4 files each).
- PSMArrowIO + the PSM tables of FeatureMap/ConsensusMapArrowIO depend on QPXFile::exportToArrow /
  importFromArrow (not in scope) and ModificationDefinitionIO::encodeByRun/collect (check Rust).
- writeTableToParquet validates QPX tables via QPXValueValidation (not ported) -> refuse with Unsupported.
- DataValue::toString for doubles == crate::param::value::format_float(v, true).
- Rust: MetaInfo = BTreeMap<String, MetaValue>; MetaValue::source_float for non-finite; stem_name in
  src/system/file.rs; UniqueIdGenerator get_uuid in src/concept/unique_id.rs; DateTime::now_utc.
- Table type decision: arrow_array::RecordBatch (one chunk per column).

## Next step
1. src/format/parquet_filter.rs (ungated, complete) + tests.
2. src/format/arrow_schema_registry.rs (complete) + tests; ParquetFile (complete).
3. ArrowIOHelpers (partial: QPX validation, isobaric labels depend on unported pieces).
4. FeatureMapArrowIO / ConsensusMapArrowIO features tables; read C++ fixtures.
5. docs/ARROW_PARQUET_SUPPORT.md, tests/data/arrow_parquet_provenance.json, report + ledger fragment.

## Incarnation 3 (19:22) found uncommitted WIP from incarnation 2, committed as-is
- src/format/parquet_filter.rs + tests/parquet_filter.rs (18:47-18:49, not yet gated)
- src/format/arrow_schema_registry.rs (19:10; generated body from scratchpad F1a-w4/gen_registry.py ->
  registry_gen.rs + hand-written registry_head.rs), not in mod.rs yet, not gated
- src/format/parquet_file.rs (19:16, possibly unfinished), not in mod.rs yet
