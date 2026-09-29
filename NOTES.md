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
- DONE (commits c7b4c68..42e7609): loader tests (12, incl. XSD under mzml-schema), SQL consumer tests (11), docs
  MS_DATA_CONSUMERS_SUPPORT.md + INDEXED_MZML_FILE_LOADER_SUPPORT.md, manifest tests/data/ms_data_consumers_provenance.json,
  writing-consumer doc/manifest amended (progress vacuous; validator = record writer gap).
- Next: sweep F7a-sweep running (log scratchpad/F7a/sweep); then measure slices: mzml(ms_data_consumers,
  indexed_mzml_file_loader), sqmass(ms_data_sql_consumer), mzml-schema(indexed_mzml_file_loader); then report+ledger;
  delete NOTES.md in final commit.
- IN PROGRESS: writing consumer now takes PeakFileOptions (set_peak_file_options; Output::set_streamed_index in
  mzml_write_options.rs), loader forwards options -> no more Unsupported refusal. Tests: loader byte-identity x9,
  writing consumer byte-identity x6 + frozen-after-header. gate7 running. Pre-change slices all green
  (mzml 12+29, sqmass 11, schema 12, msrv 12+29+11, nofeat lib 316).

## Incarnation 4 (same branch, 19:57)
- HEAD cef7ab1 committed clean; slices (17:50-17:54, after cef7ab1) all green: mzml 13+29+38, sqmass 11, schema 13,
  msrv 13+29+11+38, nofeat lib 316. Ledger fragment written (9 entries). Report NOT written. sweep2 died on spock timeout.
- sweep4 started 19:58 (scratchpad/F7a/sweep4.out). Next: report F7a.md; consider record-level CV promotion
  (validateCV_) to close the writing-consumer validator gap; delete NOTES.md in final commit.

## Incarnation 5 (branch worktree-wf_a10fdd71-2d2-4, 2026-09-29), merged wf_cd066dab-0a7-12 (ff to cf60145)
- sweep3/sweep4 produced only the header line (spock timeout). Report F7a.md still NOT written.
- Next: (1) retry sweep (sweep5) in background; (2) write report F7a.md from slice logs; (3) decide on validateCV_
  record-level CV promotion (optional); (4) delete NOTES.md in final commit.
- FINDING (inc 5): the "record-level metadata written as userParam only" claim (doc row validator_, native-diff bullet,
  ledger MzMLHandler + WritingConsumer scopes) is WRONG: spectrum/chromatogram metadata goes through
  mzml_header/write.rs:443 x.metadata("spectrum"/"chromatogram") -> same promotion as header owners, limited to
  reversible promotions (xml.rs:265-267; mapping.rs:7 SPECTRUM_CV; chromatogram never). Plain write_scalar_metadata
  only for scanList/scan/scanWindow (mzml_settings.rs:392,431,460), precursor (mzml_precursor.rs:546), product
  (mzml.rs:3309). C++ paths: 5374 spec, 5389 scanList, 5443 scan, 5462/5491 window, 4736 precursor, 4755 product, 5930 chrom.
  Fix docs + ledger wording. B8a already did ZIP input (don't duplicate). F4 did CachedMzMLHandler (7bfad445, not in
  base) -> MSDataCachedConsumer unblocked next wave. sweep5 running (task b62uh81ma).
- inc5: record promotion is REVERSIBLE-ONLY (spectrum SPECTRUM_CV, precursor ACTIVATION_METADATA, header owners);
  non-reversible source promotions (scan->spectrum meta, chromatogram title->name, scanList/scanWindow/isolationWindow
  dropped as Unhandled cvParam) deliberately not done. SUSPECTED PORT DIVERGENCE: precursor "percent collision energy"
  -> Rust promotes to cvParam MS:1000138 (obo name "normalized collision energy"), C++ checkAndGetTermByName fails ->
  userParam. Oracle ../oracle/F7a-record-cv (FileConverter -process_lowmemory + plain, ibminode06, run.sh pull).
  Added test record_metadata_is_promoted_only_where_the_reader_routes_it_back (tests/ms_data_writing_consumer.rs), ungated yet.
