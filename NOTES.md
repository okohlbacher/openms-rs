# F7a lane notes (restart-safety scratch; delete in final commit)

Lane: XML handlers + MS data consumers. Worktree branch: worktree-wf_ee239125-5be-17, base f0990c0.

## Coverage status at base (docs/core-sdk-coverage.json)
- DATAACCESS: MSDataAggregatingConsumer, MSDataCachedConsumer (TOPP FileConverter), MSDataChainingConsumer,
  MSDataSqlConsumer, MSDataStoringConsumer, NoopMSDataConsumer = unmapped;
  MSDataTransformingConsumer = evidence_requires_review (NovorAdapter); MSDataWritingConsumer = partial
  (src/format/ms_data_writing_consumer.rs; gaps: per-record semantic validation + progress).
- HANDLERS: ConsensusXMLHandler, FeatureXMLHandler, XMLHandler, MzML*, MzData, MzXML = partial (see review scopes).
- FORMAT: XMLFile partial, IndexedMzMLFileLoader partial (missing type + options_ member + OnDiscPeakMap typedef).

## Tooling (worktree guard quirks)
- Pinned C++ flattened at $SP/F7a/cpp/<Name>.{h,cpp,test} (SP = session scratchpad). Re-dump with
  `python3 $SP/F7a/dump.py FORMAT/DATAACCESS/X ...` (uses symlink $SP/F7a/cs -> .reference/openms4-core-bc9cc12/src).
- The guard refuses: a `source/` path with pipes/&&; git mixed with python/heredoc; git -C only literal path, one git per call.
- Other lanes: CachedMzMLHandler not in F7a scope; SpectrumAddition is OpenSWATH (O lanes). Rust trait =
  src/interfaces.rs MSDataConsumer (4 ops, ControlFlow). copy_spectrum_meta in src/kernel/spectrum_helper.rs:504.
  mzml_sqlite_handler.rs exists (for MSDataSqlConsumer).

## Plan (dependency order)
1. Consumer trait family: find existing IMSDataConsumer equivalent (grep ms_data_writing_consumer.rs).
2. Port Noop, Storing, Chaining, Transforming, Aggregating, Cached, Sql consumers (small headers) -> src/format/dataaccess/ or alongside.
3. IndexedMzMLFileLoader type (join load_experiment + writing consumer).
4. Re-read handler partial scopes; close small gaps.

## Design decisions (made; do not re-derive)
- Files (src/format/): noop_ms_data_consumer.rs, ms_data_storing_consumer.rs, ms_data_chaining_consumer.rs,
  ms_data_transforming_consumer.rs, ms_data_aggregating_consumer.rs (all ungated), ms_data_sql_consumer.rs
  (sqmass), indexed_mzml_file_loader.rs (mzml). Tests: tests/ms_data_consumers.rs (mzml gate, MzMLFile_1 at
  tests/data/mzml_validator/), tests/ms_data_sql_consumer.rs (sqmass), tests/indexed_mzml_file_loader.rs
  (mzml; tests/data/indexed_mzml/IndexedmzMLFile_1.mzML). Doc: docs/MS_DATA_CONSUMERS_SUPPORT.md,
  docs/INDEXED_MZML_FILE_LOADER_SUPPORT.md.
- Chaining: Vec<&'a mut dyn MSDataConsumer>; first Err/Break stops the chain for that record.
- Transforming: Option<Box<dyn FnMut(..)->Result<()> + 'a>> x3, set_*/clear_*.
- Aggregating: generic merge closure F: FnMut(&[MSSpectrum])->Result<MSSpectrum> (SpectrumAddition is an
  O-lane header; format cannot depend on analysis) -> PARTIAL, named dep ANALYSIS/OPENSWATH/SpectrumAddition.h.
  Source does NOT forward setExpectedSize/setExperimentalSettings. Explicit finish() flush.
- Storing: set_experimental_settings replaces exp.settings only; reserve is a capped hint (no abort).
- SQL consumer: handler Reject policy refuses rich metadata before RUN_EXTRA exists; plan: under full_meta
  write batches with Source policy after own aux-array/pepID checks. Run-level info written at finish.
  No class test exists -> tier 4 round-trip invariants.
- Cached consumer: needs CachedMzMLHandler writeSpectrum_/writeChromatogram_ (unported, not F7a) -> leave.
- IndexedMzMLFileLoader store(OnDisc): writer flags beyond zlib -> Unsupported (streaming consumer only has
  zlib); store() sets options.write_index = true persistently (source does). path_io::write_plain is
  idxml-gated; widen cfg to mzml.

## Next step
1. Write the five ungated consumers + tests/ms_data_consumers.rs; append pub mod lines; commit.
2. IndexedMzMLFileLoader; 3. SQL consumer; 4. docs + provenance; 5. gate on spock; 6. report + ledger.

## Incarnation 3 (branch worktree-wf_cd066dab-0a7-12), started after merging wf_ee239125-5be-17 (ff to 6fe0ff0)
- The lead's WIP commit 6fe0ff0 already holds: 5 ungated consumers, sql consumer, indexed_mzml_file_loader.rs,
  tests/ms_data_consumers.rs (824 lines), mod.rs pub mods, path_io cfg widen. NOT yet gated on spock.
- Next step: (a) review WIP files compile-sanity by reading; (b) gate on spock (--features mzml / sqmass);
  (c) tests for loader + sql consumer; (d) docs; (e) report + ledger.
- gate2 (spock, --features mzml,sqmass --test ms_data_consumers): 28/29 pass; the failure was a test input on the
  strict float boundary (d = c + RT_TOLERANCE rounds to < tol) -> input moved to 2*tol. Re-gate.
- FINDING: MSDataWritingConsumer progress gap is vacuous: ctor binds MzMLHandler's `const ProgressLogger& logger_`
  to a temporary `ProgressLogger()` (MSDataWritingConsumer.cpp:20, MzMLHandler.h:468 -> dangling ref, C++ defect
  candidate), and the streaming path (writeHeader_/writeSpectrum_/writeChromatogram_) never calls logger_ (only
  load at MzMLHandler.cpp:966-1497 and writeTo at :4763-4838 do). validator_ gap = record-level writeUserParam_
  CV promotion (validateCV_, MzMLHandler.cpp:3576-3669); Rust record writer mzml.rs write_scalar_metadata writes
  plain userParams only -> MzMLHandler writer gap, keep WritingConsumer partial naming it.
- Reference C++ readable at literal path .../OpenMS4-R/.reference/openms4-core-bc9cc12/src (guard refuses $VARS in sed).
