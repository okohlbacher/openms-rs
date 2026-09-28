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
- ORACLE ALREADY RAN (incarnation 2): ../oracle/F1a-arrow/{driver.cpp,node.sh,run.sh,results/}.
  ibminode06, Linux x86_64 Release bc9cc12-c19e494-174b576, g++ 14.4 -O2 -DNDEBUG, Arrow C++ 23.0.1,
  run twice, results/differ.txt EMPTY (identical records + parquet files). results/records.tsv
  (889 lines, TSV kind<TAB>case<TAB>value) covers: all 22 schemas field-by-field (Field::ToString), the
  nested *Type() helpers, DataType::Equals cases, validate() messages (strict/subset), ArrowIOHelpers
  getters (getStringList/getInt64/getDouble/getBool/getString), jsonEscape, readMetaValues typing (82
  cases), qpxRunFileName/ScanFormat/ScanComponents/IntensityLabel(s)/IsCanonical, isobaric channel lists,
  qpxFileMetadata keys, writeTableToParquet codecs/levels/row groups (+ files/*.parquet written by C++),
  ParquetFile writeTable/readTable/rowCount errors. results/files/*.parquet = tier-1 C++-written fixtures.
  -> Rust tests should transcribe records.tsv lines (cite file + sha256 from results/hashes.txt).

## Incarnation 4 (19:39) committed incarnation 3's uncommitted edits as 8ddbef2 (ungated)
- mod.rs now declares arrow_schema_registry + parquet_file (cfg arrow). parquet_file.rs covers every
  ParquetFile.h member. First dax compile gate running (scratchpad F1a-w4/gates/g1.log).
- g1/g2 on dax: lib + tests/parquet_filter compile with --features arrow; parquet_filter 19/19 pass.
  ParquetFilter checked member-by-member -> complete.
- Oracle part-1 caveats: (a) "exists=" suffixes in writeTableToParquet/concatenate records were evaluated
  BEFORE the write (g++ evaluates operator+ args right-to-left) -> meaningless, do not transcribe;
  (b) qpxPsmRunFileName records are all "fallback" because QPXRunFileNameKeys was built before the names
  were registered (MetaInfoRegistry index -1). Part 2 (../oracle/F1a-arrow2, launched 19:58 on ibminode06,
  pull with `bash ../oracle/F1a-arrow2/run.sh pull`) registers first + records a nested list round trip.
- Reuse: format::percolator_infile::ScanNumberPattern == getRegExFromNativeID+extractScanNumber;
  data_structures::list::create::<i32/f64/String> == ListUtils::create; system::file::stem_name;
  param::value::format_float(v,true) == DataValue::toString double; kernel::features::ColumnHeader.
  MetaValue::source_float_list cfg needs `feature = "arrow"` added (src/metadata/value.rs).
- readMetaValues int = std::stoll prefix semantics; double = glibc strtod (hex floats, ERANGE on
  inexact subnormal/zero and overflow). Write local parsers in arrow_io_helpers.rs.
- Isobaric channel names: transcribe from records.tsv "isobaric" lines (tier 1); method names from
  IsobaricQuantitationMethod.cpp METHOD_REGISTRY.
- Next: write src/format/arrow_io_helpers.rs, then tests/arrow_parquet.rs reading
  tests/data/arrow_parquet/cpp_release_records.tsv (copy of oracle records.tsv) + files/*.parquet.

## Incarnation 5 (worktree wf_cd066dab-0a7-5, branch worktree-wf_cd066dab-0a7-5) started
- fast-forwarded onto worktree-wf_ee239125-5be-10 (410faf0). Continue from incarnation 4 Next.
- Scratch: USE scratchpad/F1a/ (scratchpad is SHARED with other lanes; gates/ is not mine).
- Guard: `git -C <core> show bc9cc12:<path> > <scratch file>` works as ONE plain command (no scripts,
  no variables). Scratchpad now: /private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/*.cpp
- Oracle part 2 finished (DONE, ibminode06, same build): qpxPsmRunFileName precedence + nested list
  (parquet path l.list.element; timestamp[s] reads back ms; large_string restored).
- Plan inc.5: arrow_io_helpers.rs with ParquetWriteConfig (from MSExperimentArrowExport.h),
  KeyValueMetadata (ordered, duplicates), strtod/stoll emulation for readMetaValues, isobaric table
  transcribed from records.tsv, global UUID generator in a Mutex; then tests/arrow_parquet.rs.
- inc5 status: src/format/arrow_io_helpers.rs written; tests/arrow_parquet.rs (22 tests, record-driven)
  17/22 passed at t3; fixes for seconds timestamps (coerce on write, ms hint on read), one empty
  row group for empty tables, metadata kept on read -> t4 running. Test data in tests/data/arrow_parquet/.
- Remaining NEXT: (a) green t4 + clippy + MSRV; (b) src/format/protein_identification_arrow_io.rs
  (C++ in scratchpad/F1a/ProteinIdentificationArrowIO.cpp; Rust ProteinIdentification in
  src/identification.rs: date_time Option<String>, Tolerance enum, digestion_enzyme String,
  ModificationsDB caller-owned: modification_definitions::register_search_parameters);
  (c) FeatureMapArrowIO features table (+ TOPP fixture FileConverter_9_output.featureparquet),
  (d) docs/ARROW_PARQUET_SUPPORT.md, provenance json, THIRD_PARTY row check, report, ledger.
