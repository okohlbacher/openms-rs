# Early TOPP bundle: work packages and decisions

Working plan as of 2026-09-13. It breaks the [early TOPP build plan](EARLY_TOPP_BUILD_PLAN.md)
(FileInfo, FeatureFinderCentroided, PeakPickerHiRes) into packages that run in
parallel git worktrees. Five read-only tracers followed each tool through the
pinned C++ and this port, and a synthesizer turned their results into packages.
The synthesizer checked the most load-bearing claims against the sources itself.
Statuses below change as packages merge; the ledger and
[VALIDATION](VALIDATION.md) remain the record of what has landed.

## Pins and oracles

- **Pins.** Recorded with hashes in
  [`topp_early_bundle_provenance.json`](../tests/data/topp_early_bundle_provenance.json):
  cli `c19e494`, topp `174b576`, test-data `0cb15f2`, core `bc9cc12`. The package
  checkouts have moved past these pins, so packages read sources with
  `git show <pin>:<path>`.
- **Tool-level oracle (C1).** `../oracle/topp-early-bundle/` runs every registration and
  regression against the built product SDK. It records exit code, streams,
  output hashes and the C++ FuzzyDiff verdict, and executes everything twice to
  prove reproduction.
- **Class-level oracle (C2).** `../oracle/featurefinder-picked/` holds CMake drivers
  linked against the product SDK. They cover the helper structures, the trace
  fitters, the Levenberg-Marquardt evaluation budget, isotope patterns, picked
  scoring, and a FAIMS closure.
- **Oracle caveats.** The product SDK is a Debug build of core `4fdec46`, 23
  commits before the pin; the traced code paths are unchanged between the two.
  - Outputs are labelled `oracle-generated (tier 1 executed differential)`, or
    `adapted` when a driver was modified.
  - Debug-only precondition exits are flagged and never become Rust expectations.
  - Bitwise comparison holds only on macOS arm64, the oracle's own platform.
    Elsewhere the tolerance is 1e-9 relative.

## Decisions taken

Most of these follow from rules the project already has, so they were taken as
defaults. **D3 and D5 change or define observable tool behaviour; say so if you want
them decided differently.**

| | Question | Choice | Basis |
|---|---|---|---|
| D1 | The module ratchet rejected every new cross-module edge | Reject only edges that close a cycle (`b6e59d8`) | A cycle-free edge cannot block the workspace split. The plan's own `processing -> concept` edge would close a cycle (concept reaches processing through chemistry, kernel and analysis) and stays forbidden. |
| D2 | Levenberg-Marquardt backend for the trace fitters | Crate behind package B3's gate, with Eigen `maxfev` emulated exactly; if the gate fails, the transcription stays for FEATUREFINDER | Crates-first policy; for analytic Jacobians the evaluation accounting already matches Eigen |
| D3 | TOPPBase strictness and exit codes | Adopt the C++ behaviour for every tool, including the five existing ones | Fidelity. Strict INI and command-line validation. Phase-aware exit codes: initialisation 6, run-phase parse 3, empty input 4, other exceptions 8, bare invocation 6. One SpectraFilterWindowMower test and three missing-parameter assertions flip. |
| D4 | Four existing tools omit the DataProcessing entry C++ adds | Add it in wave 2 and re-validate without regenerating expectations | Fidelity |
| D5 | FAIMS in FeatureFinderCentroided | **Open until wave 5.** The preview refuses FAIMS input explicitly | C++ exits 8 on every FAIMS input, and its merge removes every feature. Reproducing that or shipping a corrected closure is a scientific choice. |
| D6 | Comparing XML outputs | Decoded content with FuzzyDiff's per-number rule and id exclusion, plus exact structure | [DIFFERENTIAL_VALIDATION](DIFFERENTIAL_VALIDATION.md); line-layout parity (package C7) only if reversed |
| D7 | Oracle identity | Accept the product SDK with the labels above; correct the documents that say TOPP binaries do not run (that was `topp-build/bin`) | Built C++ is an accepted development-time oracle |
| D8 | Source pins | Keep them | Refreshing is a separate checkpoint |
| D9 | GPT's uncommitted baseline | Committed in `2e85db5` | |
| D10 | Native readers stricter than the source on tool paths | Explicit source-compatibility load options on tool paths; strict library defaults | Fidelity without loosening the library |

## Ownership

- **Integrator only.** The integrator alone edits `Cargo.toml`, `src/lib.rs`,
  `src/error.rs`, module roots and registration lines, the ledger and coverage
  JSON, `SOURCE_PROVENANCE.json`, the shared provenance files, the C++ issue log,
  README, CHANGELOG, VALIDATION and the CI workflow.
- **Packages.** Each package edits only its own files. Modules are registered as
  documented stubs before a package starts, and new crates arrive as integrator
  commits.
- **Merges.** Every package is reviewed by an independent verifier before it
  merges; merges follow the dependencies below.

## Packages

Letters mark lanes:
- A: FileInfo.
- B: FeatureFinderCentroided's numerical chain (the critical path).
- C: oracles, comparison support and wrappers.
- P: PeakPickerHiRes.

| Wave | Package | Depends on |
|---|---|---|
| 1 | C1 tool-level oracle | none |
| 1 | C2 class-level C++ drivers | none |
| 1 | C3 FuzzyStringComparator port and decoded-XML comparator (test support) | none |
| 1 | CLI-1 TOPPBase lifecycle: spec suite (C4) and closure part 1 | none |
| 1 | B1 FeatureFinderAlgorithmPicked helper structures | module stubs |
| 1 | B2 source-precision isotope patterns, source trimLeft, bounding-box operations | none |
| 1 | B3 Levenberg-Marquardt crate adapter with exact `maxfev` emulation | dependency commit (the digamma lane was dropped) |
| 1 | A1 PeakTypeEstimator API and FAIMSHelper | module stubs |
| 1 | A2 C++ stream and StringUtils numeric formatting | module stubs |
| 1 | A3 mzML spectrum and scan mobility, FAIMS CV, type reset, unit-bearing IM arrays; FileHandler type detection and load options | baseline |
| 2 | B4 TraceFitter trait and GaussTraceFitter | B1, C2 |
| 2 | B5 EGHTraceFitter | B1, C2, B4's trait |
| 2 | B6 FeatureFinderAlgorithmPicked front half: parameters, validation, scoring, seeds | B1, B2, C2 |
| 2 | A4 FileInfo library: model, peak-file and featureXML branches, `-m/-p/-s`, text and TSV | A1, A2, A3, C1, C3 |
| 2 | CLI part 2: `-write_ini` parity, format checks, DataProcessing retrofit (D4) | CLI-1, A3 |
| 3 | B7 FeatureFinderAlgorithmPicked back half: isotope fit, extension, fitting, quality, parallel seed loop, overlap | B4, B5, B6, C2, C3 |
| 3 | A5 FileInfo tool, binary and preview workflows | A4, C1, C3, CLI parts 1 and 2 |
| 3 | C5 FeatureFinderCentroided wrapper: load, error branches, FAIMS refusal, annotations | B6, A1, A3, C1, C3, CLI parts 1 and 2 |
| 3 | P1 PeakPickerHiRes and SignalToNoiseEstimatorMedian parameter contract and fidelity fixes | none |
| 3 | P2 source-compatible mzML header references (D10) | baseline |
| 3 | B8 IMDataConverter splitByFAIMSCV | A1, A3 |
| 3 | B9 FeatureOverlapFilter | C2 |
| 4 | B10 FeatureFinderCentroided acceptance: FFC_1, seeds, asymmetric, debug, threads, LM budget | B7, C5, B3 |
| 4 | A6 FileInfo `-i`, `-d`, `-c` | A5, A3, C1 |
| 4 | P3 PeakPickerHiRes tool, parameter failures, `-write_ini` | P1, P2, C1, C3, CLI part 2 |
| 5 | B11 FeatureFinderCentroided FAIMS closure | B10, B8, B9, C2, D5 |
| 5 | P4 PeakPickerHiRes low-memory mode | P3 |
| 5 | C6 end-to-end chain and eight-executable release bundle | A6, B10, P3 |
| 6 | A7 FileInfo consensusXML, idXML/mzid and FASTA branches | A6 |
| 6 | A8 FileInfo schema validation, mzXML/mzData, trafoXML | A7 |
| 6 | C7 featureXML writer layout parity | B10, and only if D6 is reversed |

**Preview criteria.**
- FileInfo: A5 and A6.
- FeatureFinderCentroided on non-FAIMS centroided mzML: B10 with C5's error branches.
- PeakPickerHiRes in memory: P3.
- The bundle: C6.

The three tools and FileInfo.h stay `partial` in the ledger until waves 5 and 6
close them.

## Wave 1 status

Status on 2026-09-14. `integrate/wave1` merges ten verifier-approved package
branches onto `3e171b2`: nine, then the shared-file pass, then four `--no-ff`
merges of approved fix rounds, of which A2's is that package's first merge. The
shared files (ledger, provenance, CI, licences, C++ issue log and
documentation) follow on `integrate/wave1-shared`. [VALIDATION](VALIDATION.md)
records each package's evidence.

**Merged** (branch commit, then merge commit):

| Package | Branch | Merge |
|---|---|---|
| crate/quantile: normal quantile from `statrs` | `29bfbc4` | `823c8e9` |
| crate/mt64: MT19937-64 from `rand_mt` | `1d0e9c8` | `f913cf5` |
| crate/existing-crates: reference, URL-host and calendar-day checks | `be2bb15` | `db19767` |
| C3 FuzzyStringComparator and decoded comparison | `30d82db` | `72650a5` |
| A3 mzML mobility and FileHandler loading | `ff44438` | `c8b0141` |
| B1 FeatureFinderAlgorithmPicked helper structures | `1bd8686` | `09575c3` |
| A1 PeakTypeEstimator and FAIMSHelper | `cd7b83e` | `3d88650` |
| B2 source-precision isotopes and bounding boxes | `799f610` | `3982561` |
| CLI-1 TOPPBase lifecycle part 1 | `5ae7826` | `f6bdd99` |
| A2 FileInfo numeric text formatting, after two fix rounds | `6cec951` | `7f0a8e0` |

**Fix rounds merged** (fix commit, then merge commit):

| Package | What changed | Fix | Merge |
|---|---|---|---|
| B2 | the iridium exception in the `SourceSingle` bit-identity wording; documentation only | `42f142f` | `05ca267` |
| A1 | infinite and NaN FAIMS targets matched to the executed C++; the `estimate_type_with_limits` parity qualified | `f3c29cb` | `cab8976` |
| CLI-1 | a `-ini` that is neither a regular file nor a directory skips the readability precheck; `/dev/null` exits 3 | `86d2733` | `9612872` |

**Dropped:** crate/digamma (`d94274c`). The hand-rolled digamma stays, because
`special` measured further from the executed libOpenMS than the port's series
from the class-test start; `special` is removed. See
[THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md). B3 no longer
waits for a digamma lane.

**C1 done** (tool-level oracle). Its follow-ups:
- bind the `case.json` and `run.json` hashes in the manifest guard;
- update the moved `_work` paths;
- record crashes as an exit-code class: PeakPickerHiRes `auto_mode 1` ends with
  SIGBUS 138 or SIGSEGV 139, varying between attempts;
- state the scope of each normalisation rule.

**C2 done** (class-level drivers). C1 and C2 manifests are registered in
`SOURCE_PROVENANCE.json` at W4.4.

**D7 recorded as accepted.** The product SDK is the development-time tier-1
oracle; [DIFFERENTIAL_VALIDATION](DIFFERENTIAL_VALIDATION.md) no longer implies
that no TOPP binary runs, and records the driver build recipe.

**Open follow-ups:**
- **B1, kernel hull.** `ConvexHull2D::add_points` in `src/kernel/geometry.rs`
  merges one retention time's m/z range with `f64::min` and `f64::max`. For
  `-0.0` against `+0.0` the result differs from the source's keep-first rule
  (`ConvexHull2D.cpp:165-171`, `DBoundingBox.h:102-103`) and between machines.
  Replace them with strict comparisons and add a signed-zero hull test. The
  helper-structures ledger row stays `partial` until then.
- **CLI-1, DTAExtractor.** A reversed `-rt` range such as `70:50` must be
  ordered before filtering, as `DRange<1>` does (`DTAExtractor.cpp:144`,
  `DIntervalBase.h:85-90`); the m/z bounds stay unswapped. The C++ product SDK
  writes `DTA_RT60.0.dta` and the port writes nothing. Add a test with the fix.
- **Clippy 1.85.** `cargo +1.85.0 clippy --all-features --lib --tests -- -D warnings`
  fails with 24 older-clippy lints (`excessive_precision`, `comparison_chain`)
  in `src/math/multiple_testing.rs`, `src/processing/spline/b_spline.rs`,
  `src/comparison.rs` and `src/concept/math_functions.rs`, all predating wave 1.
  CI runs clippy only on stable.
- **A1, FAIMSHelper ledger promotion: resolved.** The `#[ignore]` on
  `get_compensation_voltages_section_through_the_mzml_reader` is removed, and
  the test passes on stable and 1.85.0 after the A3 hand-off. The infinite and
  NaN targets follow the executed C++ (`f3c29cb`), and `FAIMSHelper.h` is
  `complete`.
- **B1, IsotopeCluster.h remap (deferred).** Its `evidence_requires_review`
  status comes only from the name match with `pub struct IsotopeCluster` in
  `src/processing/deisotoping.rs`; the header is not ported. A review entry can
  set only `complete`, `partial` or `native_equivalent`, so recording it as not
  ported needs a generator change (a reviewed not-ported state) or a rename of
  the Deisotoper struct.
- **B1, MassTrace.h name collision (deferred).** `KERNEL/MassTrace.h` lists
  `src/analysis/feature_finder_picked/helper_structs.rs` as a candidate file
  because both declare `pub struct MassTrace`. That is not coverage; the header
  keeps `native_equivalent` through its own review. It needs the same generator
  change.
- **FileInfo.h evidence: resolved.** A3's citation went in the audit repair.
  A1's two anchors (`FileInfo.cpp` lines 1673 and 1599 in
  `tests/data/faims_helper_provenance.json` and
  `tests/data/peak_type_estimator_provenance.json`) and A2's `sources` entry in
  `tests/data/file_info_text_format_provenance.json` are now `context_sources`
  spelled `FORMAT/FileInfo.cpp`, with the pinned sha256.
  `core_sdk_coverage.py --write` changed only `FileInfo.h`, back to `unmapped`
  (`unmapped` +1, `evidence_requires_review` -1). Nothing of FileInfo is
  ported; A4 owns it.
- **Fixture ownership.** `tests/data/peak_type_estimator/` and
  `tests/data/faims_helper/` are A1-owned; A3, A6 and B8 reuse
  `tests/data/faims_helper/IM_FAIMS_test.mzML` instead of copying it.
  `tests/data/feature_finder_picked_helper_structs_oracle.tsv` is accepted as
  B1-owned test data.
- **Final-merge items: resolved.** `SOURCE_PROVENANCE.json` registers A1's
  (`../oracle/pte-faims-helper`), CLI-1's (`../oracle/topp-cli-lifecycle`) and
  A2's (`../oracle/file-info-text-format`, `../oracle/text-format`) oracle
  artifacts with recomputed sha256 values. The CLI-1 group manifest now lives
  at `tests/data/topp_cli_lifecycle_provenance.json`; the move changes no ledger
  row, because it cites no `src/openms/` path. The `FileInfo.cpp` anchors are
  respelled as described under FileInfo.h evidence.

**Verifier notes carried forward.** Findings of the wave-1 verifiers that no
follow-up in flight covers, each with its file and owner.
- **Open for the lead: selected-ion drift time (A3 request 5).** C++ copies a
  selected ion's drift time onto the spectrum (`MzMLFile_test.cpp:418-421`, the
  a3-format-io oracle). On an MS2 spectrum with a scan-level FAIMS voltage it
  replaces that voltage (`MzMLHandler.cpp:1871-1875`). The reader in
  `src/format/mzml.rs` does neither, and `tests/precursor_workflow.rs:157` and
  `:188` assert the unpropagated round trip, so adopting the source flips them.
  Either way, add the MS2 case to native difference 1 in
  `docs/MZML_MOBILITY_SUPPORT.md`.
- **A3, `src/format/mzml.rs`:**
  - A scan-level non-FAIMS mobility value of `-1` loads, but the writer
    preflight (line 2714) refuses the spectrum it produced; C++ writes nothing
    for it (`MzMLHandler.cpp:5417`). Refuse it on read, read it as unset, or
    document it in native difference 6.
  - A mobility value with leading whitespace (`" -35"`) is a parse error; C++
    loads it.
- **A3, `src/format/file_handler.rs`:**
  - `get_type` returns an I/O error for a directory without a recognised
    extension; C++ returns `UNKNOWN`.
  - `load_experiment_with_options` opens the file before `check_allowed`, so a
    missing file with a known but disallowed extension gives an I/O error
    instead of the disallowed-type error (`FileHandler.cpp:855-863`).
- **A3, `docs/MZML_MOBILITY_SUPPORT.md`:** native difference 7 should say that a
  file whose first five lines exceed 64 KiB can be reported `Unknown`.
- **A3, `tests/file_handler_type_detection.rs`:** the nine `getType` literals of
  `FileHandler_test.cpp:135-167` are not ported, and the support doc has no
  START_SECTION table for `FileHandler_test` and `MzMLFile_test`.
- **A3, `tests/mzml_mobility.rs` and `tests/file_handler_type_detection.rs`:**
  the duplicated `hex()` helper returns 0 for a subnormal `%a` value, because
  its single `2^(exponent - 4 * digits)` scale underflows. No current oracle
  value is subnormal.
- **C3, `tests/support/decoded_compare.rs`:** `Walk::debug_text` compares
  strings inside nested metadata (identifications, settings, source files, CV
  terms) under the FuzzyDiff number and whitespace rules, while
  `docs/FUZZY_STRING_COMPARATOR_SUPPORT.md:213` and the provenance say strings
  must agree.
- **C3, `tests/support/fuzzy_string_comparator.rs`:**
  - Only `tok_hex_vs_decimal` is documented as a token divergence. Hex floats
    are a whole class (`0X1P-2`); a signed or `+`-prefixed `nan(...)` with
    characters outside `[A-Za-z0-9_]` is consumed to a different length; and
    accepting underflow to zero is not the `std::from_chars` contract, which
    libstdc++ rejects.
  - An INI value that fails conversion exits 6; TOPPBase stores 0 and applies
    only the range check.
  - `from_ini` and `load_ini` accept content that is not ParamXML and return the
    defaults; C++ FuzzyDiff exits 3.
  - The `fuzzy_diff` rustdoc promises the C++ exit code, but two directory
    inputs exit 10 where C++ exits 0.
- **crate/quantile, `src/math/multiple_testing.rs` and
  `docs/MULTIPLE_TESTING_SUPPORT.md`:**
  - F2: the `TOLERANCE` doc (line 1121) says the replaced Acklam code missed by
    up to 2.5 million epsilons at `f64::MIN_POSITIVE`; the largest miss was
    about 5 million, at `p = 1 - 1e-13`.
  - F3: the Accuracy rustdoc (line 748), "Machines" and the `determinism` text
    in `tests/data/math_kde_provenance.json` say machines differ by one or two
    units in the last place; up to 3 were measured. The replaced code was not
    identical across machines either, and only `ln`, not `sqrt`, comes from the
    platform.
  - F4: "Limits" says the port returns `+inf` at `eps <= 2^-54`; that is the
    quantile, while `lfdr` returns `Error::InvalidValue`. Boost's `domain_error`
    for NaN, `p < 0` and `p > 1` is not mapped to the port's results.
- **crate/mt64, `docs/DECOY_REFERENCE_REVIEW.md`, `src/chemistry/decoy_random.rs`
  and `tests/data/decoy_provenance.json` (F2):**
  - "lengths 0 to 1,000" (review line 54) was 12 fixed lengths: 0, 1, 2, 3, 5,
    11, 13, 20, 31, 64, 200 and 1000.
  - The equality argument (rustdoc line 31, provenance `normalization`) covers
    only seeded states compared with seeded states.
  - "Engine in the state of `boost::mt19937_64(seed)`" (rustdoc line 38) should
    say output-equivalent.
- **crate/existing-crates:**
  - F1, `docs/IDXML_SUPPORT.md:67`: attribute values are not resolved by
    `resolve_char_ref`/`resolve_xml_entity` but by `escape::unescape`
    (`identification_xml.rs:474`). It shares the `CharRef` parser, but an
    undeclared entity there is `Error::Parse`, not `Error::Unsupported`.
  - F2, `docs/NETWORK_GET_REQUEST_SUPPORT.md:174`: `http://:[::1]:80/` cannot
    come from the fuzz generator, whose alphabet has no `8` or `0`; use a
    measured example such as `HTTPS://:[]:.`.
- **CLI-1, `src/param` and `src/cli.rs` (request 7):** align the
  `ParamUpdateReport` messages and `ParamEntry::validation_error` with
  `Param.cpp:56-166` and `1216-1374`, then delete `update_diagnostics`
  (`src/cli.rs:891`), which rebuilds the source wording.
- **CLI-1, `src/cli/context.rs` and `src/system/update_check.rs` (request 8):**
  `to_int32` is ported twice (`update_check.rs:219`, `context.rs:461`, with
  `to_double` at 515); one shared `StringUtils::toInt32`/`toDouble` port removes
  the duplicate.
- **Integration, `tests/ms_data_writing_consumer.rs` (the 1.85 failure A3's
  verifier saw): fixed.** Concurrent runs raced on a fixed `/tmp` directory;
  neither the toolchain nor the writer was at fault. The test now uses its own
  `TempDir`; the reproduction is in [VALIDATION](VALIDATION.md).
- **Integration lane, the other fixed temporary names: fixed.** The same
  pattern, which two overlapping runs can break the same way, was in
  `tests/mascot_generic.rs:293` and `:1537`, `tests/sv_out_stream.rs:434`,
  `tests/mztab_m.rs:958`, `:2425` and `:2531`, and the doctests at
  `src/format/imzml_file.rs:580` and `src/format/imzml_writer.rs:890`. Each now
  uses its own `TempDir::new_in(std::env::temp_dir(), false)`, removed on drop;
  no assertion changed. A re-grep of `src`, `tests`, `examples` and `benches`
  finds no other written path under a fixed name: every other `temp_dir()` use
  carries the process id or goes through `TempDir`, and
  `src/system/java_info.rs:214` only probes a path that must not exist.
- **CLI-1, `-ini` open failures other than NotFound and PermissionDenied
  (fix-round-2 verdict minor 1).** A `-ini` that is neither a regular file nor a
  directory and whose open fails otherwise exits 8 (`UNKNOWN_ERROR`) in the
  executed C++: a Unix socket (`EOPNOTSUPP`) and `/dev/tty` without a
  controlling terminal (`ENXIO`), both on the run path and with `-write_ini`,
  because `TextFile.cpp:44-47` throws `IOException` when `File::readable` holds
  and `TOPPBase.cpp:495-500` (cli `c19e494`) maps it to `UNKNOWN_ERROR`. The
  port exits 3, and `docs/TOPP_CLI_SUPPORT.md:143-144` and the `load_ini`
  rustdoc (`src/cli.rs:719-721`) say exit 3. Either map such an open failure to
  exit 8 with `Error: Unexpected internal error (IO error for file '<path>')`,
  keeping later read failures at 3, and add socket and tty cases to
  `ini_read_failures.sh`; or narrow both texts and record the executed exit 8
  as a documented divergence. Owner: CLI-1.
- **CLI-1, single-writer FIFO (fix-round-2 verdict minor 2).**
  `an_ini_fifo_is_read_once_a_writer_opens_it` (docstring
  `tests/topp_cli_lifecycle.rs:1206-1209`), `docs/TOPP_CLI_SUPPORT.md:105` and
  `:147-149`, and the rustdoc at `src/cli.rs:727-729` present an openable FIFO as
  a case without an oracle. It is a divergence: the C++ tool opens the INI
  twice, once for the compression peek (`XMLFile.cpp:141-147`) and again for
  xerces (`:166`), so a writer that opens the FIFO once leaves the C++ tool
  blocked (killed by the alarm, exit 142, nothing written) while the port exits
  0. State that as a deliberate tier-4 divergence, and optionally record the
  executed single-writer result as an oracle observation. Owner: CLI-1.
- **SYSTEM and CLI-1, `-in /dev/null` (fix-round-2 verdict minor 3).**
  `cli::input_file_readable` (`src/cli/context.rs:378`) asks `file::readable`,
  which answers `false` for a device or FIFO without opening it, so
  `-in /dev/null` exits 2 with a false "not readable for the current user"
  message; the executed C++ exits 4 (`INPUT_FILE_EMPTY`, "Error: File empty").
  The fix needs a readability query that never opens the file, as
  `access(R_OK)` in `File.cpp:506-514`; `unsafe` is forbidden and no `libc` or
  `rustix` dependency exists, so the crate decision comes first (Pending row in
  [THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md)). Then add the
  oracle cases `-in /dev/null` (exit 4) and a mode-000 FIFO as `-in` (exit 2).
  Owners: the SYSTEM owner for `src/system/file.rs`, CLI-1 for the cases.
- **A2, `src/format/file_info/text_format.rs` Follow-up (verdict minor).** The
  section says every caller of `param/value.rs` `format_float` and
  `format_float32` inherits the tie defect, then lists the callers
  incompletely: add the `XIC @ ...` label of `kernel::chromatogram_tools`
  (`src/kernel/chromatogram_tools.rs:268`) and `identification::protein_run`
  (`src/identification/protein_run.rs:199`), or write "including". Owner: A2.
  The consolidation itself (the private copies in `format/mascot_generic.rs`,
  `format/pepxml.rs`, `format/mzxml.rs`, `math/posterior_error_probability.rs`
  and `param/value.rs`) stays with those modules' owners.
- **A4, FileInfo documentation (A2 request 6).** `docs/FILE_INFO_SUPPORT.md`
  should link the `text_format.rs` module documentation, which is A2's support
  document, and `tests/data/file_info_text_format_provenance.json`. Optional,
  owner A2: move the embedded oracle tables (`ORACLE_DRIVER_TSV`, 69,970 bytes;
  `SWEEP_ORACLE_TSV`, 32,818 bytes) and the retained reports out of the roughly
  170 KB `tests/file_info_text_format.rs` into `tests/data` with
  `include_str!`.
- **A4, A5 and C1, using the text formatting (A2 request 7).**
  - On the TOPP tool path use `fixed_truncated`, not `fixed`: `fixed` refuses
    what `StringUtils::number` cuts (CPP-253, CPP-254), so FileInfo would be
    stricter than the source.
  - Track stream precision explicitly: default 6, then `WRITTEN_DIGITS_F32` or
    `WRITTEN_DIGITS_F64` per statistics block, persisting as it does in
    `FileInfo.cpp` (lines 2224-2254, 2330-2362 and 2404).
  - Do not copy the six-space `intensity:` padding of the retained
    `FileInfo_3` and `FileInfo_7` reports; `FileInfo.cpp:140` and `:189` at
    `bc9cc12` write one space.
- **C1 and A4, macOS report comparisons (A2 platform note).** Default-stream
  `%g` text differs between Apple libc and glibc on one class of exact ties: an
  integer-valued double below `1e15` whose round-down neighbour ends in `0`.
  Apple libc keeps the trailing zeros (a total ion current of `1463805` at
  precision 6 prints `1.46380e+06`), glibc, the C standard and the port strip
  them (`1.4638e+06`). C1's macOS comparisons must not count these, or the
  `nan`/`-nan` sign difference, as port defects.

**Still in flight:** the Boost.Regex facade (`crate/regex-facade`), in review;
the lead merges it separately.
