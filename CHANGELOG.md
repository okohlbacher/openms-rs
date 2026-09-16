# Changelog

## Unreleased

- Integrated the wave-4 performance and correctness work (2026-09-16): the
  parallel peak picker, the rewritten mzML reader, the spline scratch buffers,
  the removal of a dead validation loop, and five tool fixes, recorded in
  `docs/VALIDATION.md` and measured in `docs/BENCHMARKS.md`.
  - `PeakPickerHiRes` parallelises its spectrum loop behind the `parallel`
    feature. `-threads` is wired through a pool scoped to the picking call
    rather than to the tool body, each worker owns its spline fitter and
    scratch, and the output is **byte-identical at every worker count** and to
    the serial result. 1.66x at 32 workers in the lane's own measurement and
    2.01x end to end on the benchmark input, where the C++ picker gains 1.06x;
    peak RSS falls 886 MB. `tests/topp_threads.rs` now covers the picker as its
    sixth tool, including that it starts no worker at `-threads 1`.
  - mzML reader: rewritten for **-47.2 %** of the instructions its load path
    executed — buffer reuse across records, `decode_slice` instead of
    `decode_vec`'s zero fill, and a base64 accumulation that no longer pushes
    one `char` at a time. No decoded value changes. The port's load of a 1.2 GB
    mzML is now faster than the C++ reader's.
  - `CubicSpline2dFitter`: the cubic spline takes caller-owned scratch instead
    of allocating eight heap vectors per constructed spline. The recurrence is
    unchanged term by term and in its original evaluation order, so the
    coefficients are bit-identical by construction.
  - mzML reader: the per-record `spectrum.validate()` is **removed as dead**,
    not weakened. Every decoded float is already refused at decode time if it is
    nonfinite, so both halves of that per-peak loop were unreachable over all
    197,765,338 points of the benchmark input. The argument is written into the
    rustdoc at the call site and pinned by a 36-probe, mutation-checked test.
  - `MapNormalizer` normalised against the spectrum maximum where the source
    uses the combined maximum including chromatograms — 13.24x on 7,302 of
    87,492 intensity arrays of a 1.2 GB run. Fixed; all 87,492 now agree
    bitwise. The empty-range refusal is restored and was executed against the
    C++, which refuses the same input (`CPP-308`).
  - `MzMLSplitter` and `SpectraFilterWindowMower` process full-size input. The
    writer refused precursor references that do not resolve inside a split part
    — which the source emits and mzML 1.1 forbids (`CPP-311`) — and the window
    mower applied its 1,000,000-point cap to the whole map instead of per
    spectrum (88.4 M peaks, 88x over).
  - featureXML: read ceilings are size-derived
    (`src/format/featurexml_scaling.rs`) instead of a fixed ~12.5 MB formed by
    three limits combined with `min()`, and features are streamed rather than
    retained. The 59.6 MiB and 2.06 GiB benchmark maps load and round-trip, at
    150 MiB and 5.39 GiB peak RSS, where the port previously could not read the
    featureXML its own `FeatureFinderCentroided` writes.
  - `DTAExtractor`/`DTAFile`: the writer reproduces the source's two 15-digit
    numeric rules — 15 fraction digits for the peak m/z through
    `precisionWrapper`, 15 significant digits for the precursor mass and the
    intensity through the stream (`CPP-310`) — so the 36,443-file,
    836,505,793-byte extraction of a 1.2 GB run is byte-identical to the C++
    tool's. It previously wrote 21.2 % less. `FORMAT/DTAFile.h` enters the
    ledger as `partial`, with the load/store proton-mass asymmetry named
    (`CPP-309`).
  - `docs/BENCHMARKS.md` is replaced by the wave-4 run: all eight tools on
    full-size instrument data at 1 and 32 threads, both tables, sampled thread
    behaviour, per-tool equivalence judged separately for data and metadata, and
    the comparison against wave 3. The reviewer's corrections are carried, not
    the runner's prose: the output-byte deficit is 92-93 % XML indentation the
    port does not write; the C++ 32-thread cells' surplus threads are an idle
    library pool, not a doubled compute budget; the SHA-1 cost shares are an
    arithmetic projection, not a measured ablation; and the SpectraFilterWindowMower
    ratios are quoted as ranges.
  - Six new upstream findings, each checked against the pinned source first:
    `CPP-308` (MapNormalizer's unguarded division), `CPP-309` and `CPP-310`
    (the two `DTAFile` defects), `CPP-311` (MzMLSplitter's unresolvable
    `spectrumRef`), `CPP-312` (FeatureFinderAlgorithmPicked dividing a
    zero-width RT range) and `CPP-313` (FeatureXMLHandler's 1e5 reservation cap).

- Integrated early TOPP bundle wave 3 (2026-09-15): the trace fitters, the
  wave-3a tools, the FeatureFinderAlgorithmPicked feature stage, the Boost.Regex
  facade and seven fix lanes, recorded in `docs/VALIDATION.md`.
  - Three new TOPP tools: `PeakPickerHiRes` (in-memory centroiding, the
    `algorithm` subsection, the registered parameter failures and `-write_ini`;
    `-processOption lowmemory` refused until P4), `FileInfo` (the peak-file and
    featureXML reports with `-m`, `-p` and `-s`) and `FeatureFinderCentroided`
    (features end to end).
  - `FeatureFinderAlgorithmPicked` is complete as an algorithm: mass-trace
    extension, the isotope fit, the Gauss and EGH trace fits, the quality
    checks, cropping, abort bookkeeping and overlap resolution (B7). Two lead
    decisions applied: `mass_trace:min_spectra = 1` now follows the source
    (CPP-271) and a changed isotope abundance computes the intended two-isotope
    override (CPP-247).
  - `TraceFitter`, `GaussTraceFitter` and `EGHTraceFitter` (B4, B5), with the
    shared start-value steps, the `Param` mapping and a residual-work ceiling.
  - The Levenberg-Marquardt transcription now reproduces Eigen's own reduction
    kernels, so all 141 traced trace fits are bit-identical to Eigen 5.0.1 as
    the Linux x86_64 Release build compiles it, in every evaluation argument,
    the final parameters, the status, `nfev` and `njev` (B3b). Matching that
    build everywhere is the user's decision of 2026-09-15; three NaN-only
    transcription defects were corrected in the same pass.
  - `src/concept/boost_regex.rs`: one Boost.Regex-compatible facade over
    `fancy-regex`. Every expression it compiles gives Boost's answer, or
    construction refuses it with `Error::Unsupported`. 6,531,674 compared cases
    over 134,075 patterns with 0 mismatches, and no OpenMS-derived expression
    is refused.
  - mzML reader: ceilings are derived from the document's own size instead of
    being fixed, and `run/@startTimeStamp="-infinity"` is read as the source
    reads it, so instrument-sized files load through the tools.
  - mzML reader, `ReadOptions::source_time_array_precision`: a time array given
    in minutes at 32-bit precision is converted and narrowed back to `f32`,
    which is what the source does (`CPP-306`) and what the full-precision
    library default does not. Set by `ReadOptions::source`, so every tool path
    that reproduces source loading has it. This closed the last open
    instrument-scale difference: the picked TIC chromatogram of the 2.3 GB
    benchmark run is now bit-identical, and the defect was in the reader, not
    in `PeakPickerHiRes::pick_chromatogram`.
  - mzML writer: per-record budgets, and `indexedmzML` with real record offsets
    and a real SHA-1 `fileChecksum` by default, as `MzMLFile::store` does. The
    invented zero precursor intensity is gone.
  - `MorphologicalFilter` reproduces the source's single-sample-element end
    behaviour, and `BaselineFilter -method erosion_simple` / `dilation_simple`
    now select the simple variants, which differ there (CPP-303, CPP-304). The
    f32 subtraction change is faithfulness and readability only: it is
    bit-identical to what it replaced, and no expected value moved.
  - `-threads` reaches the five wave-2 tool bodies through a worker pool, a
    non-positive count means every available processor as the executed C++
    does, and `OMP_NUM_THREADS`/`RAYON_NUM_THREADS` are ignored. The three new
    wave-3a tools are serial and do not open a pool yet.
    `docs/TOPP_CLI_SUPPORT.md` now states that contract instead of the
    superseded one and links the new `docs/TOPP_THREADS_SUPPORT.md`.
  - The picker's acquisition-copy ledger is derived from the input, and
    `pick_experiment_in_place` avoids the whole-experiment clone.
  - [BUILD] `levenberg-marquardt` and `nalgebra` moved to `[dev-dependencies]`:
    the crate was measured and not adopted, and only the ignored candidate gate
    uses it.
  - New `docs/BENCHMARKS.md`: the C++ Release reference build, the harness, and
    the first instrument-scale comparison — `PeakPickerHiRes` on a 2.3 GB Q
    Exactive run, 22,776,198 centroids bit-identical in m/z and intensity, the
    port 38.7 s and 4,309 MiB against 26.8-28.9 s and 3,884 MiB, and the picked
    TIC chromatogram bit-identical as well once the reader reproduces `CPP-306`.
  - Logged CPP-289 to CPP-307, registered 64 oracle artifacts and four more
    reference manifests, gave the three new tools feature-sliced CI lines along
    with eleven other targets, recorded the acyclic `cli -> analysis` edge, and
    updated the ledger: `TraceFitter.h`, `GaussTraceFitter.h`,
    `EGHTraceFitter.h` and `MorphologicalFilter.h` are complete, and
    `PeakPickerHiRes`, `FileInfo` and `FeatureFinderCentroided` are validated
    TOPP workflows.

- Integrated early TOPP bundle wave 2 (2026-09-14): eight verifier-approved
  branches on `integrate/wave2`, recorded in `docs/VALIDATION.md`.
  - PeakPickerHiRes and SignalToNoiseEstimatorMedian parameter contract
    (`defaults`, `from_param`, `to_param`) and fidelity fixes: `get_type(true)`
    in `pick_experiment`, spline bisection, the source FWHM midpoints, float32
    ion-mobility products, and `PickingCompatibility` for the source's
    acceptance of degenerate input (P1).
  - mzML reader: `ReadOptions::source_dangling_references` reads dangling
    `softwareRef` and data-processing references as MzMLHandler does, and a
    whitespace-padded placeholder `softwareRef` no longer panics (P2).
  - TOPPBase lifecycle part 2: `-write_ini` files equal the C++ tools', with the
    ISO-8859-1 declaration through `paramxml::WriteOptions::source`; input
    formats through `FileHandler::get_type`; usage on stderr with the product
    version line; BaselineFilter, MapNormalizer, SpectraFilterWindowMower and
    MzMLSplitter attach their DataProcessing records (decision D4); an `-ini`
    that exists but cannot be opened exits 8 (CLI-2).
  - FileInfo library preview: DTA, DTA2D, mzML (SRM spectra converted to
    chromatograms) and featureXML with `-m`, `-p` and `-s` in text and TSV (A4).
  - FeatureFinderAlgorithmPicked front half: parameters, input checks,
    intensity, trace and isotope-pattern scores, pattern precalculation and seed
    selection; `run` returns `Unsupported` after seed selection (B6).
  - `IMDataConverter::split_by_faims_cv` (B8) and FeatureOverlapFilter with its
    quadtree, source mode only (B9).
  - The Levenberg-Marquardt crate failed its gate (B3): the Eigen transcription
    stays, with a C2 evaluation-budget differential over 29,004 budgets.
  - Logged CPP-256 to CPP-288, registered the packages' manifests and oracle
    artifacts, added their test targets to the minimum-Rust CI job, recorded
    five acyclic module edges and the quadtree's MIT notice, and updated the
    ledger: FeatureOverlapFilter.h is complete; FeatureFinderAlgorithmPicked.h,
    FileInfo.h, PeakPickerHiRes.h, both signal-to-noise estimator headers and
    IMDataConverter.h are partial.

- Integrated early TOPP bundle wave 1 and crate wave 1 (2026-09-14): nine
  verifier-approved branches, recorded in `docs/VALIDATION.md`.
  - PeakTypeEstimator's public API and FAIMSHelper (A1); mzML spectrum and scan
    mobility, the representation reset, unit-bearing ion-mobility arrays and
    FileHandler type detection with option-taking loaders (A3); the
    FeatureFinderAlgorithmPicked helper structures (B1); opt-in source-precision
    isotope patterns, the source `trimLeft` and bounding-box predicates (B2).
  - TOPPBase lifecycle closure part 1 with TOPPBase exit codes for every tool
    (CLI-1, decision D3): a bare invocation exits 6, an unreadable `-ini` 2 and a
    directory `-ini` 3.
  - FuzzyStringComparator, FuzzyDiff and a decoded featureXML/mzML comparator as
    shared test support (C3).
  - MultipleTesting's probit quantile now uses statrs `erfc_inv` in Boost's
    statement order; SpectrumCheapDPCorr divides by the run-time `sqrt(2*pi)`
    Boost uses instead of the `root_two_pi` literal (cross score +2 ulp).
  - [REFACTOR] DecoyGenerator and UniqueIdGenerator take MT19937-64 from
    `rand_mt` 6.0.3; output unchanged.
  - [REFACTOR] Reference, URL-host and calendar-day checks take quick-xml,
    `http::Uri` and chrono helpers. Signed character references are now refused
    in the shared ID and map XML reader, and empty-host and unparseable URLs are
    refused by the network preflight.
  - The hand-rolled digamma stays, because the `special` crate measured further
    from the executed libOpenMS; the unused `special` dependency is removed from
    `Cargo.toml` and `Cargo.lock`.
  - Logged 23 C++ issue candidates as CPP-230 to CPP-252, registered the new
    provenance manifests and oracle artifacts, added the packages' test targets
    to CI, and updated the ledger: PeakTypeEstimator.h is complete;
    FAIMSHelper.h, the helper structures, DBoundingBox.h and both coarse isotope
    headers are partial.
  - Wave-1 fix rounds and the final shared-file pass:
    - FileInfo's numeric text (A2): `StringUtils::number`, full-precision `toStr`
      with the shortest round-trip scientific text and `std::to_chars` ties to
      even, integer and string `toStr`, vector `operator<<` and stream `%g`
      output. Apple libc `%g` ties, the NaN sign and the INT_MAX-byte `%f` band
      are documented platform differences.
    - FAIMS targets (A1): an infinite target is filtered as the executed C++
      filters it, keeping only unannotated identifications even at an infinite
      tolerance; a NaN target stays refused. The FAIMSHelper_test reader section
      runs again, and `FAIMSHelper.h` is complete.
    - `-ini` exit codes (CLI-1): a `-ini` that is neither a regular file nor a
      directory skips the readability precheck, which reported such a file
      unreadable without opening it; the load opens it instead. `/dev/null`
      exits 3 and a FIFO the user cannot open exits 2, before a run and with
      `-write_ini`, as in the C++ tool.
    - Source-precision isotopes (B2): `ProbabilityPrecision::SourceSingle`
      reproduces the SDK bit for bit only in runs that iterate elements in the
      port's order, the majority order for natural elements, and never for a
      formula containing iridium, which the SDK builds from rhenium's tables.
    - Logged CPP-253 to CPP-255 and extended CPP-245; registered the A1, CLI-1
      and A2 oracle artifacts; `StringUtils.h` and `ListUtilsIO.h` are partial and
      `FileInfo.h` is unmapped again. Tests and doctests that wrote under fixed
      temporary names now use their own directories.

- Resumed the integrated Claude checkpoint with sqMass handler hardening: invalid
  Numpress quantization fails atomically, writes recheck database limits before
  commit, snapshot/selection allocations are charged, and schema-name shadowing
  cannot authorize views. SWATH endpoints retain documented full binary precision.
- Added nine precursor activation metadata routes and selected-ion intensity-unit
  transport, with spectrum/chromatogram regressions and explicit source differences.
- Reconciled retained worktrees, source evidence and C++ defect reports. The next
  storage stage remains SqMassFile, streaming consumers and spectrum access;
  this checkpoint does not complete the SDK.
- Fixed intermittent `FailedToStart` results in `tests/system_process.rs`. Tests
  running in parallel raced on `ETXTBSY`: a child forked by one test held the
  write handle of a script another test was about to execute. It reproduced
  under Rust 1.85 and 1.96 alike (18 of 25 parallel runs). The binary's tests now
  run serially, and 0 of 25 parallel runs fail under either toolchain.
- Ran the requested Fable review of the resume checkpoint and fixed its
  confirmed findings. Metadata-only mzML `<activation>` blocks now put every
  cvParam before the fallback userParam, as the schema requires; a database
  storing more than the 512 MiB allocation budget of array data now opens, with
  its stored size bounded by the 2 GiB database limit instead; the main sqMass
  handler's full-precision RT bounds are documented and pinned by a test. Also
  corrected CPP-219's CV accessions and a wrong feature name in the ImzMLFile
  docs. The review record is in `tests/data/sqlite_s1_resume_validation/`.
- Fixed `File::writable` reporting a file unwritable when it could be created
  within a few bytes of the platform path limit. The probe's shortest name was
  the unique-name counter, whose width grows with use: a six-digit counter
  missed five path lengths below Linux's 4095-byte limit, which made
  `writable_agrees_with_the_operating_system_at_every_path_depth` fail
  intermittently in CI. One-byte probe names now close the band, and the test
  sweeps every length one byte apart with a six-digit counter.

- Added the optional SQLite connector: three open modes, checked table/column queries and row counts, SQL batches and binary bindings. Names are treated as literal identifiers, statements use owned cleanup, and operational errors retain their SQLite cause. SQL transaction control remains with the caller; sqMass, OSW and OMS adapters are subsequent work.

- Integrated the fifteen-module FORMAT wave: mzTab records and file adapters,
  mzTab-M, streaming mzML consumers, separated-value output, mzXML, mzData,
  pepXML, qcML, MSstats, Percolator input, transformation XML, Mascot and
  mzIdentML. Public API coverage and remaining gaps are recorded in
  `docs/FORMAT_WAVE_SUPPORT.md`; this does not complete the Core SDK.
- Added bounded pepXML annotation processing and qcML child handling, and
  tightened malformed-document and encoding checks in legacy XML readers.
  Review corrections and regression evidence are recorded with the wave.

- Fixed (imzML reader): a spectrum whose m/z and intensity arrays are not both external
  now decodes — the external side from the `.ibd`, the other from its inline base64, as
  `ImzMLInterceptConsumer` does — instead of reporting a length mismatch;
  `ImzMLHandler::mz_array` / `intensity_array` now apply that same `IMS:1000101` rule, so
  `OnDiscImzMLExperiment::extract_ion_image` and `spectrum` can no longer disagree about
  the same pixel; and `ImzMLReadLimits::max_xml_bytes` now caps the reader's input, so it
  bounds the XML parser's peak buffer and not only its cumulative progress.
- Fixed: the imzML writer now writes every float cvParam with the source's 15-digit
  NumericFormatting rule (so an unset retention time is `-1.0` and a pixel size of 1e5 is
  `1.0e05`), reads its six vocabulary meta keys with the source's lenient
  `DataValue::toString()` instead of refusing a non-string value, and skips a misaligned
  auxiliary data array under every `PeakFileOptions` rather than failing whenever a sort
  or trimming filter happens to run.
- Fixed: the imzML class-test suite's float oracle omitted the opposite-sign branch of
  `ClassTest::isRealSimilar`, so it accepted any sign error whose magnitudes matched and
  all 38 assertions resting on it were weaker than they read; the oracle is now ported
  branch for branch, guarded against the two defects in upstream's own version, and pinned
  by 13 tests.
- imzML loads now charge inline peak arrays against `max_loaded_peaks`. The preflight only
  sees the index, whose lengths come from the external-array CV params, so a spectrum
  storing its peaks inline reached the caller uncounted once the reader learned to decode
  inline arrays; the only other bound was the 512 MiB XML ceiling.
- Fixed: mzML reading no longer rejects a list whose declared `count` attribute disagrees
  with the actual number of child elements. Upstream's `MzMLHandler` reads the attribute
  only for progress reporting and capacity hints and never compares it, and real files —
  including OpenMS's own class-test input `MzMLFile_1.mzML`, which declares two binary
  data arrays and carries four — carry wrong counts. The attribute is still required and
  must be numeric, every declared count that bounds a resource still rejects hostile
  values before allocation, and writing still emits the true count.
- Ported the 22 `MSExperiment_test.cpp` sections that the kernel WP7 package had mapped
  without naming a test or an asserted value — including copy and move assignment, which
  had no Rust evidence at all — rewrote the class-test accounting so every remaining
  mapped section cites a test function and one concrete value (81 sections: 56 ported, 21
  mapped, 4 unaccounted), corrected the claim that a reversed ion-mobility range "silently
  selects nothing" in the source (it throws `Exception::InvalidRange`), re-verified every
  `MSExperiment.cpp` line anchor, limited the serial rasterizer's "same image" claim to
  `RasterAggregation::Max`, and recorded the source's unreachable negative-bin skip
  branch.
- Fixed: support docs no longer disclaim ConsensusFeature's Display and Ratio::description
  or SpectrumSettings' four ion-mobility accessors as unported; the on-disc failure-path
  tests use the unmodified upstream MzMLFile_1.mzML again instead of a substitute, and the
  documented duplicate-native-identifier rule is now covered by a test that actually
  contains a duplicate.
- Kernel WP10: `FeatureMap.h` and `ConsensusMap.h` audited member by member for the first
  time and completed, and `ConversionHelper.h` ported. `kernel::map_operations` adds
  `AnnotationStatistics` with the source stream layout, `merged`/`append` for
  `operator+`/`operator+=` (caller-owned `UniqueIdGenerator`, returning the redraw count),
  `swap`/`swap_features_only`, `find_protein_identification`, the primary MS run path trio
  on both containers, the container-level `applyMemberFunction` walk,
  `annotation_statistics`, `unassigned_id_matches`, `SplitMeta`/`split`, `append_rows`,
  `append_columns`, the checked `set_experiment_type`,
  `sort_peptide_identifications_by_map_index`, `ConsensusMap::with_size` and `Display` for
  both maps; `kernel::conversion_helper` adds the three `MapConversion` overloads as named
  functions returning the new container. Identification data attaches by reference, as for
  the feature types. Source quirks preserved and documented: `swap` leaves the meta values
  behind, `appendRows` pairs column headers positionally, `appendColumns` shifts by the
  header count, `setPrimaryMSRunPath` writes by position through a default-inserting map,
  `split` default-inserts for an unmatched `map_index` and tests `COPY_FIRST` against the
  handle index. Two documented divergences: a map cleared with `clear(false)` compares
  equal to a default map (no range cache), and the peak-map conversion clamps against the
  points actually collected. All 74 class-test sections (FeatureMap 32, ConsensusMap 39,
  ConversionHelper 3) are ported; nine further C++ defects recorded.
- Kernel wave 2.
- Ion mobility on `MSSpectrum`: drift-time accessors over the `-1` sentinel, ion-mobility
  data-array detection and retrieval (`contains_im_data`, `im_data`, `maybe_im_data`),
  `sort_by_ion_mobility`, `is_sorted_by_im`, the chunked `sort_by_position_presorted` with
  `Chunk`/`Chunks`, and `rasterize_im_frame` with bounded pixel and item ceilings
  (`src/kernel/spectrum_mobility.rs`).
- Kernel `MSChromatogram.h` and `Mobilogram.h` closed. `kernel::chromatogram_merge` ports
  `mergePeaks` and its `setSumSimilarUnion` helper, exposing the `round(rt * 1000.0)`
  merge bucket as `merge_rt_key`; adds `rt_begin_in`/`rt_end_in` for the eight subrange
  search overloads, `sort_by` for the predicate sort, `source_equal`,
  `chromatogram_mz_less` and `Display`. Unsorted input is checked instead of undefined, a
  non-finite summed intensity is an error, and `MergedDataArrays` turns the source's
  silently misaligned annotation arrays into an explicit refuse/drop/replicate choice. All
  43 `MSChromatogram_test` and 48 `Mobilogram_test` sections are ported;
  `docs/MOBILOGRAM_SUPPORT.md` gains the complete member table and the inherited
  `RangeManager` mapping onto the on-demand `range_manager()`. Six further C++ defects
  recorded.
- Kernel WP8: the identification surface of `BaseFeature.h`, `Feature.h` and
  `ConsensusFeature.h` — annotation state with the source names, checked peptide-
  identification sorting, the `map_index` copy, primary IDs, observation-match sets,
  atomic reference translation, the `applyMemberFunction` subordinate traversal and
  consensus ratios. Identification data attaches by reference, so the graph or its
  `ReferenceTranslator` is a parameter instead of a `FeatureMap` member, a documented
  divergence from `FeatureMap.h:294` that keeps `Clone` and `PartialEq` on both map types.
  All 92 class-test sections of the three headers are ported; eight further C++ defects
  recorded.
- `kernel::mrm` ports `MRMFeature.h` and `MRMTransitionGroup.h`: the SRM/MRM peak group
  with its per-transition and precursor feature lists, the OpenSWATH score records, and
  the transition group with its three parallel key maps, both subset operations and the
  consistency checks. The group is generic over a native `Transition` trait because
  `ReactionMonitoringTransition` is not ported, and `is_internally_consistent` now
  actually returns `false`, where the source's release build returns `true`
  unconditionally because its three checks are `OPENMS_PRECONDITION`s.
- `format::indexed_mzml_handler` ports `IndexedMzMLHandler.h` and the record-decoding half
  of `MzMLSpectrumDecoder.h`: random access to one spectrum or chromatogram at its index
  offset. The record's byte range is trimmed at its own closing tag and wrapped in the
  file's cached header, so the existing mzML reader decodes it and the record keeps RT, MS
  level, precursors and auxiliary arrays, where the source decodes only `binaryDataArray`
  payloads and the `id`. `PeakFileOptions` filtering runs through the existing load path.
  A malformed index is rejected against explicit ceilings and the file length before any
  allocation; the source hands `endidx - startidx` straight to `new char[]` and never
  inspects the read result. Six further C++ defects recorded.
- Adopt `cargo-nextest` as the development test runner. Profiling showed 82% of the
  sweep was test execution, not compilation: `cargo test` runs each of the 225 test
  binaries in turn, nextest runs all tests in one pool. Full sweep 42 s on a 384-core
  node, test phase 48 s to 9 s there and 56 s to 18 s locally. CI keeps `cargo test`.
- Kernel wave 1. `kernel::ranges` ports `RangeManager.h`, `SpectrumRangeManager.h` and
  `ChromatogramRangeManager.h` as pure value types; `range_manager()` accessors on spectra,
  chromatograms and mobilograms and three experiment roles replace the source's cached
  `updateRanges()`, and the combined role now includes chromatogram retention time,
  intensity and product m/z, which the previous `ranges()` omitted.
- `kernel::range_utils` ports all fifteen `RangeUtils.h` predicates with the source reverse
  flag, plus `retain_spectra` and `retain_peaks_where`.
- `kernel::spectrum_helper` ports the `SpectrumHelper.h` free functions over spectra and
  chromatograms; source metadata loss is available only behind an explicit option.
- Member-by-member review of `DPeak.h`, `StandardTypes.h`, `RichPeak2D.h`, `FeatureHandle.h`
  and `BinnedSpectrum.h`, with `FeatureHandle::from_peak`, `Display`, `Hash` and `HasUniqueId`
  and the `BinnedSpectrum::DEFAULT_BIN_*` constants added.
- New `DPosition`, `DIntervalBase` and `DRange` with the exact finite-extrema empty sentinel.
- Kernel headers closed or native-equivalent: 8 of 34 before this wave, 18 after.
- Fifteen further C++ defects recorded as CPP-061 to CPP-075.
- Restore the Rust 1.85 minimum: five let-chain sites (stabilised in 1.88, silently
  accepted by newer compilers under edition 2024) are rewritten as nested `if`s.
  CI `minimum-rust` was red on `5688775`.
- Cut `target/` from 27 GB to 4 GB with `[profile.dev] debug = "line-tables-only"`;
  record the first build and test timing baseline (16 s build, 176 s suite).
- Kernel scaffold: `MSSpectrum::{drift_time, drift_time_unit}`,
  `MSExperiment::sql_run_id`, `BaseFeature::{primary_id, id_matches}`,
  `ConsensusFeature::ratios` and `Ratio`, plus `Error::{InvalidRange,
  MissingInformation}`. Every existing struct literal already used
  `..Default::default()`, so no caller changed.

## 0.2.0 — 2026-09-11

First release with executable TOPP tools.

### Core SDK

- Target advanced to `bc9cc12`. Registered public headers: 786
  (12 complete, 62 native-equivalent,
  21 partial, 167 evidence-requires-review,
  524 unmapped).
- Typed record metadata, mzML primary-array roles and an optional libxml2-backed
  `mzml-schema` validator.
- The complete experimental design and its tab-separated reader.
- The whole kernel is documented from the C++: all sixteen `src/kernel*` modules
  at 100% rustdoc coverage, crate-wide 47.5%, enforced by a per-module ratchet.

### TOPP

- A native `TOPPBase`: registration, defaults/INI/command-line resolution,
  `-write_ini`, usage text and all fifteen source exit codes.
- **Five executable tools**, each reproducing its upstream test against retained
  C++ output: `BaselineFilter`, `DTAExtractor`, `MapNormalizer`, `MzMLSplitter`
  and `SpectraFilterWindowMower`. `validated_topp_workflows` is 5, the first
  nonzero value the ledger has reported.
- Algorithm subsections, porting `getSubsectionDefaults_`.

### Evidence

- `docs/DIFFERENTIAL_VALIDATION.md` defines four evidence tiers and a
  canonical, tolerance-based comparison policy. Byte equality is not the
  contract for mzML, because the upstream suite itself uses FuzzyDiff.
- No C++ is committed to this repository; probe sources live outside it and are
  recorded by hash, enforced by `tools/check_core_sdk.py`.
- The top-level differential-testing claim was corrected: 2310 executed cases
  against extracted OpenMS translation units, and no algorithm with SDK
  dependency closure compared against running C++.

### Known limits

- 524 of 786 headers are unmapped; 141 of 146 TOPP tools are not ported.
- The port is serial; the source parallelises with OpenMP in 36 files.
- Vendor formats, HDF5 and ONNX are out of scope.

## 0.1.0 — Ongoing native Rust SDK port

- Add algorithm subsections to the TOPP framework (`Tool::subsection_defaults`, porting `getSubsectionDefaults_`) and two more tools: `MapNormalizer` and `SpectraFilterWindowMower`, both reproducing their upstream tests against retained C++ output.
- Backport the OpenMS documentation across the whole kernel: all sixteen `src/kernel*` modules reach 100% rustdoc coverage, crate-wide 42% to 48%.
- Add a rustdoc coverage ratchet (`tools/check_doc_coverage.py`, CI-enforced) and backport the full MassTrace documentation: `src/kernel/mass_trace.rs` goes from 9% to 100% documented. Tool bodies move into `src/cli/tools/` so the shipped binary and its differential test share one definition.
- Advance the Core SDK target to `bc9cc12`. Only the Parquet reader changed, which the port does not implement; registered public headers stay at 786 and no pinned reference bytes moved. Adds `tools/core_sdk_retarget.py`, which refuses to run when a build-registration input changes.
- Add the TOPP command-line framework (the native `OpenMS4-cli` TOPPBase) and the first TOPP tool, DTAExtractor. Its three upstream tests are reproduced byte-for-byte against the retained C++ outputs, making it the first validated TOPP workflow. Header list `count` attributes are now advisory on mzML reading, and `dta::WriteOptions::source()` selects the source writer conventions.
- Add the complete experimental design: both public classes, all five path/label mappings, both sample-grouping rules, the consensus/feature/identification constructors, column-header annotation and the tab-separated reader in both source table layouts. Ragged rows and negative indices are rejected rather than read out of bounds or wrapped.
- Transport ordered mzML spectrum acquisitions with explicit read normalization, source combination/CV metadata, bounded header references and schema-valid parameter order.
- Add native MassTraceDetection with all source settings, mobility-aware extension, both termination criteria, area input and cumulative bounded atomic results.
- Add the complete ProForma annotation AST and both text writers; compare 160 numerical/formatting cases with compiled exact-source writer extraction. Parser and scientific backends remain open.
- Add protein-run inference metadata, settings export, singleton groups and metadata-only copying with preserved target result ownership.

- Add the complete built-in monosaccharide database, preserving 24 records, 12 aliases and exact source mass/formula fields without a runtime dependency.
- Transport mzML scan modes, polarity, zoom, scan windows, spectrum Product lists and chromatogram types, including source file-content summaries and bounded validation.
- Add the complete MassTrace container and centroid/area/FWHM operations, preserving source cache and numerical conventions with bounded native failures.
- Expose every source numeric constant and metadata key, with all 129 values checked against an executed unchanged C++ header and existing chemistry mass paths preserved.

- Update the SDK target to `82ce5b3`, verifying the unchanged numeric formatter after its upstream relocation and retaining original fixture pins.
- Add source-compatible spectrum/chromatogram conversion, attached acquisition settings and shared processing handles. Existing constructors retain defaults; struct literals require the new fields or defaults.
- Preserve newly attached fields through processing under cumulative copy bounds; reject unsupported mzML settings before output.

- Add exact overlapping unmodified peptide-to-protein sequence coverage with source percentage arithmetic and bounded native work.

- Add all six mzML Numpress read transports and configurable writing with bounded preparation, source precision repairs and ordinary fallback.
- Add owned unique ID/UUID generation and a reusable ID value interface using the existing source-compatible random engine.
- Add plain/mobility/rich 2D peaks and checked unfiltered experiment import/export with source metadata conversion and RT grouping.
- Add public IMS integer/real decomposer APIs with source table/order/endpoint semantics, constrained queries and checked nontermination/resource boundaries.

- Add native IMS isotope distributions, elements, alphabets and replaceable text parsers with explicit source arithmetic and bounded operations.
- Add checked peak indices, borrowed scalar area traversal and filtered bulk exports preserving source append and RT grouping behavior.
- Add the independent Numpress wrapper feature for base64/zlib transport, source estimation/accuracy checks and explicit rejection/fallback reports.
- Add all four raw Numpress codecs and fixed-point helpers, with 295 executed C++ reference cases and retained original component license notices.
- Add checked mobility peak/mobilogram values, searches, stable aligned sorting, selection and summaries.
- Preserve generic array metadata and shared processing descriptions; reject lossy XML output. Existing DataArray struct literals require defaults for the new fields.
- Complete the standalone IMSWeights scaling/GCD/rounding/parent-mass utility with checked native boundaries.
- Add direct mzML file loading with scientific options, compressed input and atomic replacement/output.

- Execute source mzML loading filters and aligned sorting before native float conversion; preserve 26 canonical auxiliary binary-array roles.

- Port MassDecompositionAlgorithm with all source settings, private residue-table solver, literal count tests and independent integer-composition checks.

- Add checked source TIC binning, combined/chromatogram ranges, aligned sorting and experiment summary operations.

- Add bounded idXML path loading, atomic plain-file publication and identification dispatch with source extension/allowlist rules.

- Add complete native peak-file option state, metadata/Product hash traits, and source-compatible MassDecomposition count records with checked arithmetic.

- Update the SDK target to `54a232f`, retaining historical fixture pins and explicit changed-source compatibility reviews.
- Complete native FASTA parsing/file/stream/seek/progress lifecycle, including source modified sequences and PEFF prologue handling.
- Add bounded experiment aggregation and XIC extraction with all four source reducers and product m/z metadata.
- Preserve chromatogram product isolation/scalar metadata through mzML; add bounded indexed-mzML footer/offset parsing and index detection.
- Validate with Rust 1.98 and 1.85 and correct the newer compiler test-formatting check.

- Add owned logging with source routing, prefixes, duplicate caching and notifications, and replaceable progress backends with process CPU timing.
- Complete public modification collection for feature/consensus maps and nested subordinates; charge even empty search-name lookups.
- Add gzip/bzip2 INI loading with content detection, decoded limits and atomic parameter updates, preserving source plain output.

- Add featureXML and consensusXML adapters with typed feature metadata, processing history, assigned/unassigned identifications, checked custom chemistry, source filters, and protein-group quantity ownership guards.
- Add reusable gzip/bzip2 transport using Rust backends, atomic file output and feature/consensus FileHandler dispatch.
- Add portable modification definition records and owned-registry registration shared across all three identification XML formats.
- Add native filesystem helpers, explicit runtime/data/configuration paths and owned temporary resources, with documented platform conventions.
- Native API migration: feature/map/column metadata now uses typed MetaInfo; map records gain processing history and loaded-file path/type.

- Add native parameter values, hierarchical parameters, defaults/restrictions/update/copy/merge, forward traces and both command-line parsers, with source quirks and transactional errors documented.
- Add default-parameter lifecycle handling, pure typed-state callbacks, source warning policies and atomic leaf-key metadata export.
- Add OpenMS parameter XML/INI read/write, legacy fixtures, UTF-8/Latin-1/UTF-16 input, full-precision special floats and checked metadata preservation.
- Add complete source TextFile/CsvFile helpers and ListUtils/StringListUtils operations, using native streams, slices and standard containers.
- Document four reviewed standard-container/alias equivalents separately from scientific class coverage.
- Add an exhaustive completion ledger for 786 registered public SDK headers and direct dependencies from 146 TOPP sources, with reviewed coverage distinguished from names and source-reference evidence.
- Port all 73 file-type descriptors and lexical filename helpers; add native experiment dispatch, bounded format sniffing, gzip transport and atomic path output.
- Add source MS2/DTA2D readers, DTA2D filters/storage/TIC, and a checked native MS2 writer with original source fixtures.
- Add mzML referenceable parameter groups, forward header references, shared inline validation and bounded definition/reference expansion.
- Add graph referential cleanup with all five source switches and conditional predicate filtering; surviving IDs remain stable and removed IDs cannot be reused.
- Add parent/match grouping, legacy sequence/evidence converters and charged peptide fragment mass queries, with focused source-derived tests.

- Update the target to reduced Core SDK 4.0.0 at `6bfc0e4`; record retained/removed source scope and 220 unchanged historical reference paths, expose the exact target in Rust, and verify SDK/RNA resource consistency in CI.
- Add the native identification sequence/provenance graph: typed graph-owned references, input/software/search/score history, parent/peptide/oligo registration, translated merge/copy, and inclusive sequence coverage with checked atomic updates.
- Add both RNase graph digestion operations, preserving custom RNA chemistry, parent flanks and processing history with cumulative parsing/digestion limits. Add independent source references, cross-parent budget/rollback regressions and a provenance-aware RNA example.
- Add typed peptide fragment formula dispatch plus identification observations, compounds, adducts, molecule references and observation matches. Registration, source merge semantics, score and annotation history, best-match queries, graph-owned translation and atomic resource checks are covered by focused source-derived tests.

- Add all fourteen RNA enzymes and native digestion, fixed/variable modified-RNA enumeration, and both annotated RNA spectrum-generation operations, with source-derived processing references and digestion-to-mzML workflows.
- Add RNA nucleoside records, the complete pinned 378-record registry, bounded TSV and optional MODOMICS JSON providers, and owned nucleic-acid sequences with source terminal/linkage/slicing and fragment mass conventions.
- Add independent RNA source formulas, masses, slices and complete registry projections, negative-ion isotope/mzML workflows, data provenance and a regeneration tool. MODOMICS redistribution terms remain separately tracked for release.
- Add the complete native Tagger API for measured spectra and m/z arrays: source mass/charge/tolerance traversal, fixed/variable modifications, I/L alternatives, deterministic two-stage registry lookup and bounded atomic append.
- Add all six original tag-count and 120 membership assertions, independent edge/mass-table references, digestion/indexing/mzML workflows and an extraction example.
- Add native AdductInfo parsing, electron-aware mass conversion and mono/average shifts, signed formula compatibility, source integer/whitespace grammar and independent ion-composition/mzML workflows.

- Add native DecoyGenerator with all source reversal/shuffle operations, seeded MT19937-64 and Boost interval mapping, preserved cache/cleavage quirks, transactional state and shared work/output limits.
- Add all thirteen literal decoy reference cases, independently reproducible integer RNG checks, FASTA/indexing/FDR/idXML workflows and a FASTA decoy example; preserve the helper's Boost license notices.

- Added all three SpectrumAnnotator operations and four IonNaming helpers: matched spectrum arrays, hit peak annotations, matching statistics and bounded charge display/parsing. Preserved ppm duplicate conventions, source no-op flags and f32 errors, with atomic updates and checked undefined statistics.
- Shared actual isotope/alignment work and precharged fragment/loss work across candidate peptides; coarse convolution also shares its allowance across theoretical envelopes. Anonymous mass-tag spellings now use shared immutable ownership during peptide slicing without changing their text or equality.
- Added original annotation fixtures, independent error/ratio references, modified-digestion and XML workflows, cumulative-work regressions and an annotation example.

- Added charge and isoelectric point with four pKa scales; seven hydrophobicity scales, GRAVY, moving profiles and moments; ten AAindex scales and source residue indicators; and gas-phase basicity. Preserved parent-residue and terminal-annotation conventions with checked inputs and work limits.
- Gas basicity preserves ordinary source arithmetic, uses a stable energy-domain retry for overflow and corrects empty-sequence high-temperature cancellation. Added complete table fixtures, independent numerical references, digestion/idXML workflows and a peptide-property example.

- Added an owning fallible fine-isotope iterator, absolute/relative threshold selection, raw/log probability access, checked peak conversion and custom binary64 isotope populations for both streaming and materialized patterns. Iterator formula input follows the raw wrapper's charge-ignoring convention; the high-level generator retains its natural-H convention.
- Shared the existing bounded isotope search with incremental iteration, preserving deferred expansion, source rounding and cumulative theoretical-spectrum work limits. Added an enrichment/streaming example.

- Added native fine isotope configurations with absolute/relative thresholds and total-probability coverage, preserving source abundance/output rounding, fixed labels and natural-H charge handling without a C++ dependency.
- Added fine theoretical fragment, neutral-loss and precursor envelopes with shared work limits and atomic append. Isotope loss intensities now retain the full floating-point product until final f32 storage, including the existing coarse path.
- Native API: `TheoreticalIsotopeModel::Fine { unexplained_probability }` selects fine envelopes; this enum now implements `PartialEq` without `Eq` because the probability is floating point. IsoSpec layer order remains outside the implemented surface.

- Added scalar precursor purity, SPS fragment matching, fuzzy scan purity, RT interpolation and experiment scoring with source numerical conventions and shared work limits.
- Connected precursor acquisition fields directly to spectrum/chromatogram records, including isolation offsets, activation, mobility and spectrum references. Parent lookup honors earlier referenced scans before acquisition-order fallback.
- Expanded the mzML subset to preserve supported precursor acquisition fields and reject unrepresentable native fields before writing. DTA/MGF exports reject richer precursor metadata rather than discarding it.
- Native API migration: `Precursor` is now `Clone` rather than `Copy`; `Precursor::new` is no longer const. Complete struct literals need `..Precursor::default()`. `PrecursorInfo` owns only `peak: Precursor`; field access delegates to that record through `Deref`/`DerefMut`, with explicit conversions in both directions.

- Added internal b/a fragments, abundant immonium peaks, activation-method presets and the compact mass-only ladder helper, preserving source charge, terminal, loss, ordering and floating-point conventions with bounded atomic output.

- Added shared ownership for modification records, caller-owned registry parsing/setters/generation, and exact custom-registry idXML interchange without leaked storage.
- Added bounded native OBO loading, PSI-MOD alias semantics, pinned XLMOD monolinks and a separate crosslink lookup database with preserved specificity and mass conventions.
- Added anonymous modification definitions and inference, retaining spelling, source compatibility behavior, and full-residue/terminal absolute-mass anchors.
- Added absolute-formula replacement for changed residues, including restoration of unknown B/Z/X chemistry; unchanged-residue records keep their original formula and mass.
- Preserved distinct same-name custom chemistry in protein modification observations, sequence duplicate filtering, peptide identity keys and conflict resolution.
- Native API migration: registry UniMod IDs are optional, registry entries and known sequence annotations use `Arc<ResidueModification>`, and `to_unimod_string()` returns `Result<String>` with non-UniMod mass fallback. `to_accession_string()` retains vocabulary accessions; peptide identity keys now own complete sequence values.

- Added fixed and variable peptide-modification generation with bounded combinatorial output, source terminal/ordering conventions, preserved existing annotations and atomic append.
- Added modification definition sets, compatibility checks, inference from peptide identifications and delta/absolute mass matching, including explicit stored absolute masses on owned modification records.
- Added a digest-variant enumeration example and a workflow verifying chemical formulas, independent fragment shifts, inferred search definitions and idXML round trips. The identification example now uses the fixed-modification generator.
- Preserved custom formula-free modification masses through peptide generation and monoisotopic fragments, with explicit unavailable composition and source absolute-mass precedence.
- idXML writing now rejects typed modification placements that cannot round-trip through its sequence syntax before touching the destination.

- Added native iterative peak picking with HiRes seeds, original-index regions, source refinement/suppression conventions and aligned integrated-intensity/width annotations.
- Added sliding and jumping window filters with source window/duplicate rules, stable ties and atomic spectrum/experiment updates.
- Added iterative mean noise estimation with three clipping passes, source histogram arithmetic, sparse-window diagnostics and checked historical percentile behavior.
- Added an end-to-end profile noise → iterative centroiding → window filtering example and exact mzML annotation round-trip tests. Independent reviews cover fourteen refinement cases, 480 window-selection cases and 192 mean-noise configurations.

- Added native EMG fitting with source training selection, analytic gradients, iRprop+ optimization, best-fit diagnostics and bounded truncated-side reconstruction.
- Added optional EMG preprocessing for peak integration, baseline estimation and shape metrics, with typed parameters and explicit finite bounds.
- Added a cropped-peak fitting example, independent area/centroid identities and a fitted-trace mzML workflow. The scientific core now uses the pure Rust `libm` dependency for the complementary error function.

- Added native legacy/corrected chromatogram peak picking with independent seed/boundary noise settings, source overlap handling, exact original indices and aligned peak annotations.
- Added sampled peak integration, all source integration/baseline choices and shape metrics for spectra and chromatograms, preserving source floating-point and Simpson conventions.
- Added a chromatogram integration example and a workflow covering exact boundaries, raw intensity sums, time-weighted areas and mzML annotation interchange.
- Added named mzML float/integer/ASCII string arrays, preserved empty strings and annotation placeholders, and bounded cumulative decoded array storage with validation before output.

- Added source Poisson/KL deisotoping with threshold/top-N preprocessing, longest-cluster selection, original-index membership, source shared-isotope behavior, optional disjoint clusters and checked aligned annotations.
- Corrected the threshold filter default to source 0.05; preserved direct Poisson recurrence rounding with a stable overflow fallback.
- Added explicit isotope sharing to the simple deisotoper and retained exact source precursor arithmetic in both methods.

- Expanded `AASequence` to unresolved B/Z/X and owned numeric mass tags with source precision-dependent registry lookup; formula and mass APIs now return `Result`.
- Preserved numeric annotations through digestion, indexing, idXML, modification mapping and conflict/filter operations; monoisotopic fragments support known mass without invented composition.

- Added peptide-to-protein indexing with ambiguity/mismatch rules, enzyme-aware evidence, I/L handling, decoy inference and atomic run reconstruction.
- Added basic protein score aggregation, representative counts, indistinguishable groups and greedy resolution for vector, single-run and consensus-map inputs.
- Added feature/spectrum identification conflict resolution, file-origin partitions and a FASTA-to-inference/FDR/idXML workflow test.
- Added score categories and atomic score switching, HyperScore/Morpheus fragment scoring, and native peptide/protein identification filters.
- Added legacy/Basic target-decoy FDR, peptide/protein q-values, picked proteins/groups, posterior estimates and ROC with pinned-source fixtures.
- Added optional bounded idXML 1.5 read/write with typed metadata, run/evidence references, independent schema validation and a complete identification-processing round trip.
- Added typed metadata, CV terms and validated acquisition/settings records.
- Added peptide/protein identification records, evidence coverage, observed protein modifications and spectrum/feature attachments; peak writers reject unsupported identification loss.
- Added retention-time regression, interpolation and robust LOWESS models and atomic experiment/feature/consensus transformations with source-derived reference tests.
- Expanded native digestion to all 33 pinned enzymes with full/semi/nonspecific enumeration and checked validity/count operations.
- Added a synthetic peptide-identification workflow example with checked fragment matches, evidence and protein coverage.
- Added native spectra, chromatograms, experiments and basic precursor metadata.
- Added feature/consensus containers, checked maps and IDs, scan-envelope geometry, containment and consensus/decharge summaries.
- Added formula/element chemistry, modified peptide masses, b/y fragments and tryptic digestion.
- Added an embedded UniMod/custom modification registry with preserved data/source licenses.
- Added coarse isotope patterns, convolution, enrichment, averagine and conditional fragment estimates.
- Added configurable theoretical peptide spectra with terminal ion series, neutral losses, precursor peaks, coarse isotope envelopes and aligned annotations.
- Added spectrum alignment, sparse binning, and common spectrum similarity scores.
- Added Gaussian/Savitzky–Golay smoothing and morphological baseline correction.
- Added HiRes spectrum/chromatogram centroiding, natural cubic splines, histogram-median noise estimation, FWHM and peak boundaries.
- Added simple C13-spacing deisotoping with optional charge conversion, intensity summation and annotations.
- Added normalization, filtering, scaling and linear resampling.
- Added DTA, FASTA and MGF interchange plus an optional documented mzML subset.
- Added source inventory, coverage notes, fixtures, examples, tests and CI configuration.
- Added a modified-peptide mass, isotope and theoretical-spectrum example.
- Recorded deliberate differences and remaining C++ library scope.

This version is an initial API and may change. It has not been published.
