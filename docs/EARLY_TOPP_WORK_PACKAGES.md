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
| D5 | FAIMS in FeatureFinderCentroided | **Settled in wave 6 (2026-09-18): the corrected closure.** Package B11 replaced the refusal; the tool splits by compensation voltage, runs once per voltage, annotates `FAIMS_CV` and merges across voltages with `FaimsMergeFidelity::Corrected` | C++ exits 8 on every FAIMS input, and past that failure its merge removes every feature. Reproducing that or shipping a corrected closure was the scientific choice; the corrected one was taken, with the source merge kept and still tested beside it. `IMDataConverter` stays `partial` for the members D5 left out. |
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
| 1, carried into 2 | B3 Levenberg-Marquardt crate adapter with exact `maxfev` emulation (gate failed; D2 fallback merged in wave 2) | dependency commit (the digamma lane was dropped) |
| 1 | A1 PeakTypeEstimator API and FAIMSHelper | module stubs |
| 1 | A2 C++ stream and StringUtils numeric formatting | module stubs |
| 1 | A3 mzML spectrum and scan mobility, FAIMS CV, type reset, unit-bearing IM arrays; FileHandler type detection and load options | baseline |
| 2 | B4 TraceFitter trait and GaussTraceFitter | B1, C2 |
| 2 | B5 EGHTraceFitter | B1, C2, B4's trait |
| 2 | B6 FeatureFinderAlgorithmPicked front half: parameters, validation, scoring, seeds | B1, B2, C2 |
| 2 | A4 FileInfo library: model, peak-file and featureXML branches, `-m/-p/-s`, text and TSV | A1, A2, A3, C1, C3 |
| 2 | CLI-2 TOPPBase lifecycle part 2: `-write_ini` parity, format checks, DataProcessing retrofit (D4) | CLI-1, A3 |
| 2 | P1 PeakPickerHiRes and SignalToNoiseEstimatorMedian parameter contract and fidelity fixes (pulled forward from wave 3) | none |
| 2 | P2 source-compatible mzML header references (D10) (pulled forward from wave 3) | baseline |
| 2 | B8 IMDataConverter splitByFAIMSCV (pulled forward from wave 3) | A1, A3 |
| 2 | B9 FeatureOverlapFilter (pulled forward from wave 3) | C2 |
| 3 | B7 FeatureFinderAlgorithmPicked back half: isotope fit, extension, fitting, quality, parallel seed loop, overlap | B4, B5, B6, C2, C3 |
| 3 | A5 FileInfo tool, binary and preview workflows | A4, C1, C3, CLI parts 1 and 2 |
| 3 | C5 FeatureFinderCentroided wrapper: load, error branches, FAIMS refusal, annotations | B6, A1, A3, C1, C3, CLI parts 1 and 2 |
| 4 | B10 FeatureFinderCentroided acceptance: FFC_1, seeds, asymmetric, debug, threads, LM budget | B7, C5, B3 |
| 4 | A6 FileInfo `-i`, `-d`, `-c` | A5, A3, C1 |
| 4 | P3 PeakPickerHiRes tool, parameter failures, `-write_ini` | P1, P2, C1, C3, CLI part 2 |
| 5 | B11 FeatureFinderCentroided FAIMS closure: the split, one run and seed filter per voltage, the `FAIMS_CV` annotation and the corrected cross-voltage merge | B10, B8, B9, C2, D5 |
| 5 | P4 PeakPickerHiRes low-memory mode | P3 |
| 5 | C6 end-to-end chain and eight-executable release bundle | A6, B10, P3 |
| 6 | ~~A7 FileInfo consensusXML, idXML/mzid and FASTA branches~~ **landed**, wave 8 (`port/a7-fileinfo` `7801a05`, merge `b2a173c`) | A6 |
| 6 | A8 FileInfo schema validation, mzXML/mzData, trafoXML | A7 |
| 6 | C7 featureXML writer layout parity | B10, and only if D6 is reversed |

**Preview criteria.**
- FileInfo: A5, A6 and A7 (all done; the three flags and the consensusXML, identification and FASTA branches are ported).
- FeatureFinderCentroided on non-FAIMS centroided mzML: B10 with C5's error branches.
- PeakPickerHiRes in memory: P3.
- The bundle: C6.

The three tools stay `partial` in the ledger until waves 5 and 6 close them;
FileInfo.h is `partial` since A4 and stays so until A8 lands; A6 closed `-i`, `-d` and `-c` and A7 the consensusXML, identification and FASTA branches, which leaves `-v` the only refused flag and pepXML, mzTab, trafoXML and PQP the only refused branches.

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
  (fix-round-2 verdict minor 1). Closed by CLI-2 (`f886d90`): exit 8 with the
  source message, socket and tty oracle cases added.** A `-ini` that is neither a regular file nor a
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
- **CLI-1, single-writer FIFO (fix-round-2 verdict minor 2). Closed by CLI-2
  (`f886d90`): documented as a deliberate difference with the executed
  observation (CPP-265).**
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

## Wave 2

Plan as of 2026-09-14, on `main` after `1d80eed`. Ten packages run in parallel.
The integrator's scaffold commit registers every new module they need, so no
package edits a module root or a registration line.

The user's priority of 2026-09-14 is a Release benchmark of the Rust tools
against C++: PeakPickerHiRes first, then FeatureFinderCentroided. P1, P2, B8 and
B9 are therefore pulled forward from wave 3.

### Packages and dependencies

Every dependency below is merged on `main` or done outside the repository: B1
`09575c3`, B2 `3982561`, A1 `3d88650`, A2 `7f0a8e0`, A3 `c8b0141`, C3 `72650a5`,
CLI-1 `f6bdd99` (fix round `9612872`), C1 and C2 under `../oracle/`.

| Package | Scope | Depends on | Merge order |
|---|---|---|---|
| B3-LM | `levenberg-marquardt` adapter behind the unchanged `minimize`, exact Eigen `maxfev` emulation, budget differential | the pinned crates (`levenberg-marquardt =0.14.0`, `nalgebra =0.33.3`), in `Cargo.toml`; moved to `[dev-dependencies]` in wave 3 | first in the fitter lane if its gate passes; if not, D2's fallback, reported to the integrator. Outcome: the gate failed on parameter fidelity, D2's fallback applied, and the candidate survives only as an `#[ignore]`d measurement |
| B4-GAUSS | TraceFitter driver, defaults and `Param` mapping; GaussTraceFitter | B1, C2, the scaffold's trait | before B5. Outcome: B3 merged first (`14284fb`) and changed no backend, so no rerun was needed there; the reruns that mattered came at B3b's merge, where every B4 bound was re-measured (one rose from 1e-3 to 1e-2, the rest fell) |
| B5-EGH | EGHTraceFitter | B1, C2, the scaffold's trait, B4's `optimize` | after B4, rebased on it |
| B6-FFAP-SEEDS | FeatureFinderAlgorithmPicked parameters, `run()` validation, scoring, pattern precalculation, seed selection | B1, B2, C2 | independent; `algorithm.rs` passes to B7 in wave 3 |
| A4-FILEINFO-CORE | FileInfo library: model, peak-file and featureXML branches, `-m/-p/-s`, text and TSV | A1, A2, A3, C1, C3 | independent; before A5 |
| CLI-2 | TOPPBase lifecycle part 2: input and output format checks, `-write_ini` parity, usage on stderr, DataProcessing retrofit (D4) | CLI-1, A3 | independent; before A5, C5 and P3 |
| P1-PICKER-LIB | PeakPickerHiRes and SignalToNoiseEstimatorMedian parameter contract and fidelity fixes | none | before P3 |
| P2-MZML-LENIENCY | Source-compatible dangling `softwareRef` and `defaultDataProcessingRef` (D10) | baseline | before P3; the integrator lands the FileHandler call site |
| B8-IMSPLIT | `IMDataConverter::splitByFAIMSCV` only (D5) | A1, A3 | independent; before B11 |
| B9-OVERLAP | FeatureOverlapFilter with its quadtree, source mode only (D5) | C2 | independent; before B11 |
| B11-FAIMS | The FeatureFinderCentroided FAIMS closure; `FaimsMergeFidelity` in `FeatureOverlapFilter` | B8, B9, C5/B10 | closes D5 for this tool; the split and the filter are unchanged apart from the one added mode |

### Decisions in force

- **D1.** The module ratchet rejects only a new edge that closes a cycle. Each
  package reports its new `crate::X` edges; the integrator records them at merge.
- **D2, updated.** `levenberg-marquardt =0.14.0` and `nalgebra =0.33.3` are pinned
  in `Cargo.toml`. Only B3 edits `src/math/fitters/levenberg_marquardt.rs`.
- **D3.** TOPP tools adopt TOPPBase's strictness and its phase-aware exit codes.
- **D4.** CLI-2 retrofits the DataProcessing entry into the existing tools.
- **D5, updated.** The FAIMS closure of FeatureFinderCentroided is deferred. B8
  ports only the library `splitByFAIMSCV`, and B9 only the source mode.
- **D6.** XML outputs are compared by decoded content (FuzzyDiff's per-number
  rule, id exclusion, exact structure), never by line layout.
- **D7.** The product SDK (Debug, core `4fdec46`) is the development-time oracle.
  Its outputs are `tier 1 executed differential`; Debug-only preconditions are
  flagged `debug_only`.
- **D8.** Pins cli `c19e494`, topp `174b576`, test-data `0cb15f2`, core `bc9cc12`.
  Package sources are read with `git show <pin>:<path>`.
- **D10.** Tool paths get explicit source-compatibility load options; library
  defaults stay strict.
- **Benchmark.** P1, P2, B8 and B9, and the B-lane code the benchmark later runs,
  avoid needless per-peak allocations and quadratic passes on hot paths.
  Fidelity comes first, and a parallel result stays bit-identical to the serial
  one.

### The scaffold

New files, each with a `//!` doc naming the package that fills it. All but
`trace_fitter.rs` export nothing yet; that one holds the contract described
below.

- `src/analysis/feature_finder_picked/`: `trace_fitter.rs` (with the trait
  below), `gauss_trace_fitter.rs`, `egh_trace_fitter.rs`, `algorithm.rs`,
  `scoring.rs` and `seeds.rs`, registered in `src/analysis/feature_finder_picked.rs`.
- `src/format/file_info/`: `model.rs`, `report.rs`, `peaks.rs` and `features.rs`,
  registered in `src/format/file_info.rs`.
- `src/kernel/im_data_converter.rs`, registered in `src/kernel.rs`.
- `src/processing/feature_overlap_filter.rs`, registered in `src/processing.rs`,
  and `src/processing/feature_overlap_filter/quadtree.rs`, registered in that
  file.

Registrations are unconditional. A package whose whole module needs a feature
starts its file with an inner `#![cfg(feature = "...")]`, for example A4's
`features.rs` on `featurexml`.

Not registered, because no wave-2 package fills them: the tool modules
`file_info`, `feature_finder_centroided` and `peak_picker_hi_res` (for A5, C5 and
P3), and B7's `extension`, `fitting` and `resolution`. They are registered when
those packages start. P1 needs no tool module: `processing -> param` closes no
cycle, so its `Param` defaults stay in `processing`.

**Rustdoc links in module docs.** Each registration carries an outer `///`,
because `tools/check_doc_coverage.py` counts `pub mod` lines as items. rustdoc
then resolves the intra-doc links of that module's `//!` docs from the parent
module. The stubs contain no links. A package filling a `//!` doc writes its
links crate-absolute (`crate::...`), as `helper_structs.rs` does, or not at all.
Links in item-level `///` docs are unaffected.

**The TraceFitter contract.** `trace_fitter.rs` declares `pub trait
TraceFitter` and `pub struct TraceFitterParams`, transcribed from `TraceFitter.h`
at core `bc9cc12`. They hold signatures and documentation only: the trait has no
provided method, and the record has no `Default`. B4 adds the rest of the file:
the `optimize` driver, the defaults, the `Param` mapping and the shared helpers.
B5 builds against the declarations and does not edit the file. **A change to a
trait signature or to the record's fields needs the integrator.** The contract
differs from the draft in the bundle plan's B4 notes as follows:

- `TraceFitterParams::max_iteration` is `i64`, not `usize`. The source member is
  `SignedSize`, the parameter has no minimum, and a value of zero or less fails
  every fit, because Eigen refuses `maxfev <= 0`.
- `compute_theoretical` returns `Result<f64>`. An out-of-range `k` gives
  `Error::InvalidValue` where the source indexes without a check. It is a
  required method, and implementors keep the source formula.
- `parameters()` and `set_parameters()` are added. They are the typed form of
  the inherited `getParameters` and `setParameters`, which the algorithm calls
  on the fitter `chooseTraceFitter_` returns. With them it can use a
  `dyn TraceFitter`, since the trait is object safe.
- `fit` borrows the traces immutably, and `area`, both span checks and
  `gnuplot_formula` take `&self`. The source declares them non-const, but
  neither subclass writes.
- The span checks are documented as the subclasses implement them: `true`
  reports the violation that `checkFeatureQuality_` rejects. The header's
  prose says the opposite.
- `optimize(x, values, residual, jacobian, &TraceFitterParams) -> Result<()>` is
  B4's free function. Its signature is fixed in the module documentation. Until
  B4 merges, B5 may call `minimize` through a private stand-in, and must delete
  it before its own merge.

`tools/core_sdk_coverage.py` matches declarations to headers by name, so `pub
trait TraceFitter` moves `TraceFitter.h` from `unmapped` to
`evidence_requires_review`, with `trace_fitter.rs` as its candidate file; the
ledger is regenerated in the scaffold commit (unmapped 429 to 428). Only the
interface is declared. B4 earns a reviewed status when it merges.

**Module edges.** The scaffold adds none: `check_module_cycles.py` reports 58
edges and 13 mutual pairs, as at `1d80eed`. The table lists the edges packages
are expected to add. Each was tested against the recorded graph with the
checker's `reaches`, and none closes a cycle.

| Edge | Package | Why |
|---|---|---|
| `format -> math` | A4 | `SummaryStatistics` |
| `analysis -> math` | B4, B6 | `minimize` (B5 reaches it through B4's `optimize`), `pearson_correlation_coefficient` |
| `analysis -> param` | B4, B6 | `Param` defaults and mapping |
| `processing -> param` | P1 | `PeakPickerHiRes` and `SignalToNoiseEstimatorMedian` defaults |
| `processing -> concept` | B9 | the `FAIMS_CV` meta-value key |

`math`, `param` and `concept` name no other top-level module, so these edges
cannot close a cycle together either. **No wave-2 package may add an edge out of
`math`, `param` or `concept`.** At `1d80eed` these edges would close a cycle and
are refused: `processing -> format`, `processing -> analysis`,
`processing -> system`, `kernel -> processing` and `analysis -> system`. The D1
row above says `processing -> concept` closes a cycle. That held when D1 was
taken, through `concept -> chemistry`; `a18c1f1` removed that edge, and at
`1d80eed` `concept` reaches nothing. The checker re-evaluates the edge when B9
merges.

### Ownership

Each file has exactly one owning package. A path ending in `/` covers the
directory.

| Package | Files it owns |
|---|---|
| B3-LM | `src/math/fitters/levenberg_marquardt.rs`, `tests/math_distribution_fitters.rs`, `tests/lm_budget_differential.rs`, `tests/data/lm_budget_differential/` (added at integration), `tests/data/lm_budget_differential_provenance.json` (added at integration), `docs/DISTRIBUTION_FITTERS_SUPPORT.md` |
| B4-GAUSS | `src/analysis/feature_finder_picked/trace_fitter.rs` (the declarations are the contract), `src/analysis/feature_finder_picked/gauss_trace_fitter.rs`, `tests/trace_fitter.rs`, `tests/gauss_trace_fitter.rs`, `tests/data/gauss_trace_fitter/`, `tests/data/gauss_trace_fitter_provenance.json`, `docs/TRACE_FITTER_SUPPORT.md` |
| B5-EGH | `src/analysis/feature_finder_picked/egh_trace_fitter.rs`, `tests/egh_trace_fitter.rs`, `tests/data/egh_trace_fitter/`, `tests/data/egh_trace_fitter_provenance.json`, `docs/EGH_TRACE_FITTER_SUPPORT.md` |
| B6-FFAP-SEEDS | `src/analysis/feature_finder_picked/algorithm.rs`, `src/analysis/feature_finder_picked/scoring.rs`, `src/analysis/feature_finder_picked/seeds.rs`, `tests/feature_finder_picked_seeds.rs`, `tests/data/feature_finder_picked/`, `tests/data/feature_finder_picked_provenance.json`, `docs/FEATURE_FINDER_PICKED_SUPPORT.md` |
| A4-FILEINFO-CORE | `src/format/file_info/model.rs`, `src/format/file_info/report.rs`, `src/format/file_info/peaks.rs`, `src/format/file_info/features.rs`, `tests/file_info.rs`, `tests/data/file_info/`, `tests/data/file_info_provenance.json`, `docs/FILE_INFO_SUPPORT.md` |
| CLI-2 | `src/cli.rs` (except its `mod` lines), `src/cli/context.rs`, `src/cli/spec.rs`, `src/cli/parameter.rs`, `src/cli/usage.rs`, `src/cli/processing.rs`, `src/cli/tools/baseline_filter.rs`, `src/cli/tools/dta_extractor.rs`, `src/cli/tools/map_normalizer.rs`, `src/cli/tools/mzml_splitter.rs`, `src/cli/tools/spectra_filter_window_mower.rs`, `src/format/paramxml.rs` (writer options only), `tests/paramxml.rs` (writer-option tests only), `docs/PARAMXML_SUPPORT.md` (writer-option section only), `tests/topp_cli_lifecycle.rs`, `tests/data/topp_cli_lifecycle/`, `tests/data/topp_cli_lifecycle_provenance.json`, `tests/topp_baseline_filter.rs`, `tests/topp_dta_extractor.rs`, `tests/topp_map_normalizer.rs`, `tests/topp_mzml_splitter.rs`, `tests/topp_spectra_filter_window_mower.rs`, `docs/TOPP_CLI_SUPPORT.md` |
| P1-PICKER-LIB | `src/processing/peak_picking.rs`, `src/processing/peak_picking/noise.rs`, `src/processing/iterative.rs`, `src/processing/chromatogram.rs`, `tests/peak_picking.rs`, `tests/peak_picking_experiment.rs`, `tests/data/peak_picking/`, `tests/data/peak_picking_provenance.json`, `docs/PEAK_PICKING_SUPPORT.md`; only where a fidelity fix forces it: `tests/iterative_picking.rs`, `tests/iterative_picking_reference.rs`, `tests/iterative_workflow.rs`, `tests/chromatogram_picking.rs`, `tests/chromatogram_processing_reference.rs`, `tests/chromatogram_workflow.rs`, `tests/data/iterative_provenance.json`, `tests/data/chromatogram_processing_provenance.json`, `docs/ITERATIVE_PICKING_SUPPORT.md`, `docs/ITERATIVE_REFERENCE_REVIEW.md`, `docs/CHROMATOGRAM_PICKING_SUPPORT.md` |
| P2-MZML-LENIENCY | `src/format/mzml_header/read.rs`, `src/format/mzml_header.rs` (forwarding the option only), `src/format/mzml.rs` (`ReadOptions` and the header-registry call only), `src/format/mzml_counts.rs` (the header-registry call only), `tests/mzml_header_leniency.rs`, `tests/data/mzml_header_leniency/`, `tests/data/mzml_header_leniency_provenance.json`, `docs/MZML_HEADER_SUPPORT.md` |
| B8-IMSPLIT | `src/kernel/im_data_converter.rs`, `tests/im_data_converter.rs`, `docs/IM_DATA_CONVERTER_SUPPORT.md`, `tests/data/im_data_converter/` (added at integration), `tests/data/im_data_converter_provenance.json` |
| B9-OVERLAP | `src/processing/feature_overlap_filter.rs` (except the `quadtree` registration), `src/processing/feature_overlap_filter/quadtree.rs`, `tests/feature_overlap_filter.rs`, `docs/FEATURE_OVERLAP_FILTER_SUPPORT.md`, `tests/data/feature_overlap_filter/` (added at integration), `tests/data/feature_overlap_filter_provenance.json` |

**Boundaries.**

- **Trace fitters (B3, B4, B5).**
  - B3 keeps the public surface of `levenberg_marquardt.rs` unchanged:
    `minimize`, `LmParameters`, `LmStatus` and `DenseMatrix`. B4 codes against
    it while B3 swaps the backend.
  - B3's `tests/lm_budget_differential.rs` writes out the Gauss and EGH
    residuals and Jacobians itself, from the class-test functors, rather than
    importing B4's or B5's modules.
  - B3 hands its `tests/data/distribution_fitters_provenance.json` delta and its
    `docs/THIRD_PARTY_CRATE_DECISIONS.md` status to the integrator as text.
  - `src/math/fitters/{mod,gauss,gamma,gumbel,gumbel_max_likelihood}.rs` stay
    unedited. `minimize` keeps its signature, so they need no change.
  - B4 may extend the rustdoc of the trait and the record, and add impls, but
    changes no signature.
  - B5 reads `trace_fitter.rs` and never edits it. EGH's `tau` and `sigma` are
    inherent methods.
- **B6.**
  - `helper_structs.rs` (B1), `src/chemistry/isotopes.rs` and
    `src/kernel/geometry.rs` (B2) are read-only. The open signed-zero hull
    follow-up in `geometry.rs` is not B6's.
  - `run()` past seed selection returns `Error::Unsupported` until B7.
- **A4.**
  - `src/format/file_info.rs` and its re-exports are the integrator's; A4 asks
    for the re-exports at merge.
  - `text_format.rs` is A2's and read-only; `checks.rs` is A6's, and A4 does not
    create it.
  - `FileHandler` and the mzML and featureXML readers are read-only.
  - The A2 notes carried forward above apply: use `fixed_truncated`, track the
    stream precision, write one space after `intensity:`, and do not count the
    macOS `%g` ties as defects.
- **CLI-2.**
  - The `mod` and `pub mod` lines of `src/cli.rs`, all of `src/cli/tools.rs`,
    `src/bin/` and `Cargo.toml` stay the integrator's. A new `cli` module is
    registered on request.
  - The five tools' shared manifest `tests/data/topp_cli_provenance.json` stays
    the integrator's; CLI-2 hands its delta over.
  - `FileHandler::get_type` is read-only. CLI-2 handles its error at the call
    site; that a directory gives an I/O error where C++ gives `UNKNOWN` stays an
    open A3 follow-up.
  - In `src/format/paramxml.rs`, CLI-2 adds writer options only; the reader and
    the existing defaults stay unchanged.
  - The CLI-1 notes carried forward that fall inside these files (the `-ini`
    open failure that C++ exits 8 on, the single-writer FIFO) may be closed
    here. Request 7 (`src/param`), request 8 (`src/system/update_check.rs`) and
    `-in /dev/null` (`src/system/file.rs`) are outside CLI-2 and stay open.
  - The retrofit meets the native difference "Non-finite values in processing
    records" in `docs/TOPP_CLI_SUPPORT.md`: it records such values in a
    representable form or keeps that documented difference.
- **P1.**
  - P1 keeps the picker API used outside its files: the public fields of
    `PeakPickerHiRes`, `pick_spectrum`, `pick_experiment`,
    `pick_spectrum_with_acquisition` (the acquisition tests in
    `src/processing.rs`), `CubicSpline2d` (`src/analysis/transformations.rs`) and
    `estimate_spectrum_type` (A1's `src/format/peak_type_estimator.rs`,
    `tests/spectrum_type.rs`, `tests/peak_type_estimator.rs`).
  - The struct literals in `tests/processing_acquisition.rs` and
    `tests/data_array_descriptions.rs` use `..Default::default()`, so new fields
    break nothing.
  - The conditional files are edited only when a fidelity fix changes their
    results, with the root cause stated; retained fixtures such as
    `tests/data/peak_picking_*.tsv` are never regenerated.
  - `src/processing/spline/` is read-only.
- **P2.**
  - The strict reader stays the default.
  - The option is added to `mzml::ReadOptions` and passed through to the header
    registry. P2 changes nothing else in `mzml.rs`, and does not take up A3's
    open `mzml.rs` follow-ups.
  - The one-line call site in `src/format/file_handler.rs` is the integrator's,
    at merge.
  - `tests/mzml_header.rs` stays green and unedited.
- **B8.**
  - `src/kernel/faims_helper.rs` and `tests/data/faims_helper/` (A1) are
    read-only; B8 reuses `IM_FAIMS_test.mzML` in place.
  - The missing `updateRanges` in the C++ groups is documented as a native
    difference, not emulated.
  - `kernel -> processing` closes a cycle and is refused.
- **B9.**
  - The `pub mod quadtree;` line and its doc in `feature_overlap_filter.rs` are
    the integrator's.
  - A qualifying quadtree crate is a dependency request, not a `Cargo.toml`
    edit.
  - Source mode only (D5).

**Held by no wave-2 package.** A change to one of these goes to the integrator,
who assigns it:

- `src/format/file_handler.rs`, and `src/format/mzml.rs` outside P2's part
- `src/kernel/faims_helper.rs`, `src/format/peak_type_estimator.rs` and
  `src/format/file_info/text_format.rs`
- `src/analysis/feature_finder_picked/helper_structs.rs`,
  `src/chemistry/isotopes.rs` and `src/kernel/geometry.rs`
- `src/math/statistic_functions.rs`, the other `src/math/fitters/` files,
  `src/processing/spline/`, `src/param/` and `src/system/`
- `tests/support/`, `tests/processing_acquisition.rs`,
  `tests/data_array_descriptions.rs`, `tests/experimental_settings.rs`,
  `tests/workflows.rs`, `tests/spectrum_type.rs`, `tests/peak_type_estimator.rs`,
  `tests/mean_noise.rs`, `tests/mzml_header.rs`,
  `tests/file_handler_type_detection.rs` and `tests/mzml_mobility.rs`

**Integrator only.** The list under Ownership above applies. It also covers:

- every registration line of this scaffold, including the `quadtree` line and
  the `mod` lines of `src/cli.rs`
- `src/cli/tools.rs` and `src/bin/`
- `docs/module-cycles.json` and `docs/doc-coverage.json` (packages run the
  checkers and report; the integrator records)
- `docs/core-sdk-coverage.json` and `docs/core-sdk-reviewed-apis.json`
- `tests/data/topp_cli_provenance.json`,
  `tests/data/distribution_fitters_provenance.json` and
  `tests/data/topp_early_bundle_provenance.json`

### Gate hosts

Run cargo gates through `~/.local/bin/openms-kim-gate.sh <slot> <cargo args>`.
`/scratch` is node-local, so each slot stays on one host.

| Role | Host | Slot |
|---|---|---|
| Implementers | dax (`OPENMS_GATE_HOST=dax`) | `w2-<package>`, for example `w2-b4-gauss` |
| Reviewers | spock (`OPENMS_GATE_HOST=spock`) | `w2-<package>-review` |
| Re-reviews of fix rounds | kim (`OPENMS_GATE_HOST=kim`) | `w2-<package>-rereview` |
| Integrator | kim | `w2-scaffold`, `w2-integrate` |

- Run gates one at a time per slot. Exit 255 is a dropped SSH connection;
  rerun the gate.
- A package's set: `cargo fmt --all -- --check` locally.
  `+1.85.0 check --locked --all-features --all-targets`.
  `clippy --locked --all-features --all-targets -- -D warnings`. Its own tests on
  stable and on `+1.85.0`, under the feature line CI will use.
  `doc --locked --all-features --no-deps` (the script sets `-D warnings`).
  The `tools/` checkers locally.

## Wave 2 status

Status on 2026-09-14. `integrate/wave2` (`1c14d60`) merges eight
verifier-approved package branches with `--no-ff` onto the scaffold `42c21e7`.
A4 and B9 merged after one fix round each. The shared files (module graph,
ledger, provenance, CI, crate register, licences, C++ issue log and
documentation) follow on `integrate/wave2-shared`.
[VALIDATION](VALIDATION.md) records each package's evidence, the verifiers'
reruns, the lead's kim results and this pass's gates.

**Merged** (branch commit, then merge commit):

| Package | Branch | Merge | Outcome |
|---|---|---|---|
| B3-LM | `f604b4a` | `14284fb` | partial: the crate gate failed and D2's fallback keeps the Eigen transcription; the budget differential against C2 is committed |
| P1-PICKER-LIB | `1566a5f` | `6386e10` | done |
| P2-MZML-LENIENCY | `15a5f09` | `fabaea4` | done; FileHandler does not pass the option yet |
| CLI-2 | `f886d90` | `b4f4456` | done (W2.2); closes the CLI-1 carried-forward `-ini` open failure (exit 8) and single-writer FIFO notes |
| A4-FILEINFO-CORE | `bb28715` (fix round on `61ad528`) | `93bd214` | done; acceptance 4 holds on derived inputs only, for mzML reader gaps outside A4 |
| B6-FFAP-SEEDS | `80bbdf1` | `a532c31` | done; three decisions open (below) |
| B8-IMSPLIT | `c527c5f` | `09ce3d4` | done |
| B9-OVERLAP | `f1c7785` (fix round on `413ec7a`) | `1c14d60` | done |

**Not merged:**

- **B4-GAUSS** is in fix round 3. W2.4's `gauss_trace_fitter` and `trace_fitter`
  CI lines wait for it, and the ledger rows of `TraceFitter.h` and
  `GaussTraceFitter.h` stay unreviewed.
- **B5-EGH** merges after B4, rebased on it (the TraceFitter contract); its CI
  line and `EGHTraceFitter.h` review wait too.
- **Lane B3b** is root-causing the transcription's own first-step divergence
  from Eigen, under investigation. B3's verifier notes below go there.
- **Wave 3a (P3, A5, C5)** is being implemented on a scaffold commit above
  `integrate/wave2`; the lead merges their shared files later.
- **The Boost.Regex facade** (`crate/regex-facade`) stays with the lead.

B4's budget-boundary tests need no rerun for B3: the backend did not change.

**Ownership amendments** recorded in the table above: B3-LM owns
`tests/data/lm_budget_differential/` and
`tests/data/lm_budget_differential_provenance.json`; B8-IMSPLIT owns
`tests/data/im_data_converter/`; B9-OVERLAP owns
`tests/data/feature_overlap_filter/`. No other package touched these paths.
CLI-2's rustdoc on the paramxml reader items and one general encoding line in
`docs/PARAMXML_SUPPORT.md` exceed its "writer options only" boundary; they are
doc-only and accepted.

**Decisions open for the lead:**

- **B6 (a)** `mass_trace:min_spectra = 1`: the port refuses it; the executed
  C++ is defined (NaN trace scores, 0 seeds, exit 0; CPP-271). Following the
  source is a one-line change in `SeedStage::compute`.
- **B6 (b)** a changed `abundance_12C` or `abundance_14N`: refused by default,
  `AbundanceOverride::Intended` builds the intended distribution; the source's
  stray `(0, 1)` peak cannot be reproduced (CPP-247). The same question sits in
  front of B10.
- **B6 (c)** the overall score: the port's correctly rounded, platform-independent
  power against the platform `powf` (CPP-272). Either way, restate acceptance
  criterion 4 ("one binary32 step from the Apple `powf` oracle" or a
  macOS-only bitwise contract).
- **P2 request 1 and D10.** P2 asks to pass
  `ReadOptions { source_dangling_references: true }` unconditionally in
  `FileHandler::load_experiment_with_options` (`src/format/file_handler.rs:297`).
  Its verifier notes that FileHandler is a library API and D10 keeps library
  defaults strict, so a tool-side load option is the alternative. Not applied
  in this pass; P3's `TOPP_PeakPickerHiRes_5` and A4's `c1_empty_mzml_mps` depend
  on it.
- **B3 request 2.** `levenberg-marquardt` and `nalgebra` are used only by the
  ignored candidate in `tests/lm_budget_differential.rs`: move both pins to
  `[dev-dependencies]` and reword their `Cargo.toml` comment, or remove them with
  the candidate. Not applied; lane B3b may still need them.
- **Ledger calls made here, open to reversal.** `GaussTraceFitter.h` and
  `EGHTraceFitter.h` are `evidence_requires_review` only through B3's manifest,
  which cites their functors as `sources` (B3 verifier, major finding). The
  pass keeps that, because the test transcribes those functors; the wave-1
  alternative, respelling the citations as `context_sources` as done for
  `FileInfo.cpp`, would edit a manifest lane B3b may still change.
  `FeatureOverlapFilter.h` is `complete` (every member, tier 1), where B9's
  first report proposed `partial` until B11. `SignalToNoiseEstimatorMedian.h`
  stays `partial` for its `AUTOMAXBYPERCENT` refusal.

**Carried forward** (owner in bold at the end of each item):

- **The six new ignored tests.** Each names a documented gap; none hides an
  unexplained failure. `levenberg_marquardt_crate_candidate_gate_report` (a
  measurement, **B3b**). `c1_file_info_9_mzml_mps` and
  `a4_file_info_9_default_flags` (repeated userParam, processing on primary
  arrays, a 64-bit float charge array) and `a4_indexed_file_info_12_all_flags`
  (the float charge array): **mzML reader owner, D10**. `c1_empty_mzml_mps`:
  the FileHandler call site above and a load-options field on
  `file_info::model::Options`, **integrator and A5**. `a4_mzml_file_1_all_flags`:
  the selected-ion drift time (A3 request 5), **the lead**. The tripwire
  `reader_gaps_behind_the_ignored_cases_are_still_present` fails when a gap
  closes.
- **B3, `src/math/fitters/levenberg_marquardt.rs` and
  `docs/DISTRIBUTION_FITTERS_SUPPORT.md`:**
  - the `minimize` rustdoc (lines 829-834) and the tier-1 paragraph say the four
    degenerate fits match at every max_fev and that the termination order is
    checked; they are checked at budget 500 only, and the fixture reaches
    statuses 1 to 5 only;
  - section 8 names `1.49012e-8` as the crate's default tolerance, which is the
    `minpack-compat` default; without it the default is `30 * f64::EPSILON`;
  - `tests/lm_budget_differential.rs`: `Run` and `run_distribution` count
    njev on numerical Jacobians, which Eigen does not; and
    `the_written_out_functors_reproduce_the_c2_start_evaluations` asserts the
    tolerance, not the bit identity the docs claim (measured bit-identical on
    Linux) and does not count Jacobian comparisons;
  - performance: `minimize` allocates short-lived vectors in every outer and
    `lmpar` iteration (Jacobian column copies for `blue_norm`, the QR copy,
    Householder vectors, projected residuals, step and candidate vectors, the
    `lmpar`/`qrsolv` workspaces). Hoisting them without changing the arithmetic
    order is the likely way to close the 1.26x (dax) to 1.28x (spock) gap to the
    rejected crate on FeatureFinderCentroided's hot path;
    `tests/lm_budget_differential.rs` and the class tests guard the bits;
  - `minimize` has no internal allocation ceiling, so B4's TraceFitter must run
    `preflight_points` on peak counts; C2's `trace_fitters_degenerate_values_lt_inputs_*`
    records are B4's boundary;
  - the **evaluation-budget boundaries themselves are B10's**, not B3's and not B4's: B4 pins the
    boundary of every class-test fit and every FeatureFinderCentroided_1 seed, but the sweep C1 defines
    (`FFC_max_iterations_symmetric_01..60` and `_asymmetric_01..10`, where
    `-algorithm:fit:max_iterations` 40 differs from 500 and 50 equals it, and 6 differs and 8 equals in
    the asymmetric case) has never been run on either side;
  - the x tolerance margin (6.36e-10 against 1e-9) was measured on Linux only;
    run the cross-platform `workflow_dispatch` job before relying on it on
    macOS or Windows.
  **B3b** (performance: **owner to assign**; cross-platform run: **the lead**).
- **P1, `src/processing/peak_picking.rs`, `noise.rs` and their tests:**
  - `PickingCompatibility::source()` still refuses non-finite narrowed outputs:
    a zero intensity total gives the C++ mobility `-inf`, and the port returns
    "intensity overflow"; document it at `source()` and give it its own message;
  - no test pins the percentage formula (`count * 100 / n` against
    `count * (100 / n)`); use 3 of 7 windows;
  - `agrees()` in `tests/peak_picking_experiment.rs` compares only the error
    class for `InvalidValue`; require the source fragment;
  - benchmark: `pick_experiment` clones the whole input and
    `pick_spectrum_with_acquisition` clones each picked spectrum again; one fixed
    `AcquisitionCopies` ledger (256 MiB) charges every spectrum twice and refuses
    about 400k spectra that C++ picks; records are validated twice and
    `get_type(true)` copies each unknown-type spectrum;
  - P3 hand-off: load with `FileHandler::load_experiment_with_options` and default
    `PeakFileOptions`, build with `PeakPickerHiRes::from_param` and
    `PickingCompatibility::source()`, map `CENTROIDED_INPUT_MESSAGE` to exit 8;
  - a Release C++ benchmark must use libOpenMS's `-ffp-contract=off` for any
    header-only template built outside it, and may bin differently from arm64
    for quotients beyond `INT_MAX` (CPP-257).
  **P1** (hand-off: **P3**; benchmark build: **benchmark lane**).
- **Stale docs outside P1** (P1 request 4): `docs/CUBIC_SPLINE2D_SUPPORT.md:35`,
  `docs/SPLINE_BISECTION_SUPPORT.md:30-35` and the `CubicSpline2d::peak_maximum`
  rustdoc in `src/processing/spline/cubic.rs` still say peak picking uses
  `peak_maximum` (it uses `spline_bisection`); `docs/MOBILOGRAM_SUPPORT.md:168-172`
  says the weighted mobility arithmetic is unchanged (it now uses float32
  products in source order). **Spline and mobilogram owners (held by no wave-2
  package)**.
- **mzML reader, duplicate record ids** (P1 request 6): the native reader refuses
  the repeated `spectrum=1` of the FTMS class-test files (CPP-261) that
  `MzMLHandler` loads; P1 uses id-renamed copies. A D10 candidate. **mzML
  reader owner**.
- **P2, `src/format/mzml_header/read.rs` and `docs/MZML_HEADER_SUPPORT.md`:**
  - whitespace-padded IDs are normalised and resolve under the option, where the
    source treats them as dangling; add a native difference (resolving against
    the raw attribute would touch `parameter_id` in `mzml.rs`, integrator-held);
  - one `transform` call warns twice (setup and data pass each build a registry);
    share the de-duplication or document it, and correct the field rustdoc, the
    support doc and the P3 note;
  - a dangling lookup is charged twice against the header work allowance; drop
    the second `id.len() * 64` charge;
  - an empty optional `dataProcessingRef` is treated as absent by the source
    (`XMLHandler.h:563-571`); the "malformed IDs" native difference and the
    `ReadOptions` rustdoc should say so;
  - untested: `IndexedMzMLHandler::open_with_limits` with the flag probably warns
    on every fetched record (**P4**);
  - the option leaves `sampleRef`, `defaultInstrumentConfigurationRef` and scan
    `instrumentConfigurationRef` strict, where C++ is equally lenient; a
    follow-up if a tool input needs them. The `sourceFileRef` variants were on
    this list and came off it in wave 8 under decision D14, because the port's
    own writer produces them; with the option on, the Rust tool prints warnings
    where C++ is silent, or warns once per ID where C++ warns once per
    occurrence (**P3**).
  **P2**.
- **CLI-2, `docs/TOPP_CLI_SUPPORT.md`, `src/cli/usage.rs` and
  `tests/topp_cli_lifecycle.rs`:**
  - the D4 comparisons normalise the port's software name (`MS:1000799` with the
    name against C++'s `MS:1002146` "TOPP SpectraFilterWindowMower"),
    processingMethod order and completion-time seconds into the C++ writer's
    form; move these from preserved conventions to native differences, and give
    `src/format/mzml_header/write.rs` a writer-side follow-up (**mzML writer
    owner**);
  - a NaN BaselineFilter `struc_elem_length` exits 0 in C++ (recording NaN) and 6
    in the port; the "non-finite values" note covers only inf;
  - a directory `-in` exits 8 in the port and 3 in C++ (also named `d.mzML`);
    record it next to the unknown-content case (exit 6 against 3), which is a
    FileHandler follow-up (**A3 FileHandler owner**);
  - `an_undetermined_input_format_only_warns` asserts only "not 0 and not 1";
    pin `IllegalParameters` with a comment;
  - usage text wraps valid strings in quotes without `StringUtils::quote`'s
    escaping of `\` and `'` (latent);
  - still open from CLI-1: DTAExtractor's reversed `-rt` range and BaselineFilter's
    centroided-input warning.
  **CLI-2**.
- **A4, `src/format/file_info/`:**
  - the SRM conversion uses `ChromatogramTools::default()` limits; the 5e7 work
    ceiling refuses about 430k SRM points (4,300 spectra of 100 product m/z) that
    C++ reports; size the limits from the reader's, or document the effective
    bound in native difference 7 and `run()` `# Errors`, with a boundary test;
  - re-exports in `src/format/file_info.rs` (`FileInfo`, `Options`,
    `FileInfoResult` and the model types), which A5 imports, are not applied in
    this pass (**the lead, with A5's scaffold**);
  - a native load-options field on `Options`, strict by default, once the
    FileHandler call site takes read options (**A5 or the integrator**);
  - `ProcessingStep.h` now lists `src/format/file_info/model.rs` as a candidate
    through the `pub struct ProcessingStep` name match, the same generator
    limitation as `MassTrace.h` (**deferred with it**).
  - `FileInfo::run` returns an unknown-type result for a directory whose name
    gives no type, as C++ does, so the FileInfo tool need not handle
    `FileHandler::get_type`'s I/O error at its call site; that A3 follow-up stays
    open for other callers (**A5**);
  - A2's wave-1 requests 6 and 7 are closed: the support document links
    `text_format.rs` and its manifest, and the reports use `fixed_truncated`,
    explicit stream precision and one space after `intensity:`.
  **A4**.
- **B6, `src/analysis/feature_finder_picked/` and
  `docs/FEATURE_FINDER_PICKED_SUPPORT.md`:**
  - one isotope window failing in the generator (first at about 273,750 Da of
    `max_mz * charge_high`) fails the whole run, where the SDK continues with NaN
    weights; add a native difference or keep NaN windows for the affected peaks;
  - a zero-intensity peak takes the Release NaN intensity score where the Debug
    SDK throws a postcondition (`.cpp:1892`); flag it `debug_only`, and correct
    line 198: the tool filters negative intensities only (CPP-237);
  - `overall_score_is_the_float_cube_root_of_the_product` repeats the
    implementation's expression; assert an independent literal;
  - B7 hand-off: reuse `tests/data/feature_finder_picked/FeatureFinderAlgorithmPicked.mzML`
    and `.ini`; build step 3.3 on `SeedStage`; interleave the "Found N feature
    candidates" lines per charge; seeds with an exact `f32` intensity tie are
    ordered stably, where libc++ `std::sort` gives a different order that feature
    suppression depends on; assign ownership for any B7 edit of `seeds.rs` or
    `scoring.rs`;
  - B6's notes for this document: the plan's "undefined behaviour" premise for
    `min_spectra = 1` was refuted by execution; the 12C = 90% part of B2's pending
    pattern comparison cannot be reproduced by design; the tolerance for the
    overall score off macOS is one binary32 step, with the oracle's `powf`
    misrounding listed in `overall_rounding.tsv`.
  **B6 and B7**.
- **B8, `src/kernel/im_data_converter.rs` and
  `docs/IM_DATA_CONVERTER_SUPPORT.md`:**
  - the settings ceiling covers one copy while the split makes one per group,
    and `messages.push` is infallible; charge all copies against one budget or
    narrow the "allocation failure is refused" wording;
  - a group keyed by an infinite voltage gets `InvalidValue` from the range
    manager; qualify "ranges are always correct" (**also a B11 note**);
  - `FaimsSplit` derives `Default` with empty `groups` although the field doc
    says "never empty";
  - the derived C2 fixtures depend on A3's `FeatureFinderCentroided_1_input.mzML`
    staying byte-identical (a length and FNV-1a check fails loudly);
  - B11 hand-off: `FaimsSplit::has_faims()`, the groups with
    `FaimsGroupKey::volts()` and `messages` to relay to the tool log; seed
    filtering and `mergeFAIMSFeatures` stay with B11 and B9 (**B11**);
  - optional: a `kernel.rs` re-export of the module, now used by path
    (**integrator, if wanted**).
  **B8**.
- **B9, `src/processing/feature_overlap_filter.rs`:**
  - **Closed.** B9 request 4 asked for `git diff 4fdec46 bc9cc12 -- src/openms/extern/Quadtree`, which
    it could not run from an isolated worktree. The wave-3 pass ran it against
    `.reference/openms4-core-bc9cc12`: the diff is **empty** — all five files
    (`CMakeLists.txt`, `LICENSE`, `include/{Box,Quadtree,Vector2}.h`) are byte-identical between the
    oracle's core revision and the pin. The vendored quadtree's identity no longer rests on hashes plus
    replica equality alone;
  - `merge_faims_features` reproduces the uid-0 wipe, so FeatureFinderCentroided
    output without unique IDs loses every FAIMS feature unless IDs are assigned
    before the merge or D5 selects a corrected mode (a new opt-in API)
    (**B11, D5**);
  - a feature holding a non-empty and an empty hull is refused ("extent does not
    fit f32"); a Release C++ build probably completes; not executed;
  - hull-mode boxes are recomputed per `get_box` call where C++ caches the
    multi-hull box (constant factor, off the benchmark path);
  - `feature_overlap_filter_oracle.tsv` is 1.4 MB, already losslessly compacted,
    the largest fixture after the two controlled-vocabulary tables.
  **B9**.

## Wave 3 status

Status on 2026-09-15. `integrate/wave2` is at `864b295`, 82 commits above
`origin/main` `1d80eed` and unpushed. On top of the wave-2 tip `1c14d60` it now
carries the two trace fitters, the wave-3a scaffold and its three tools, the
FeatureFinderAlgorithmPicked feature stage, the Boost.Regex facade, seven fix
lanes and the Levenberg-Marquardt rewrite. The shared files (CI, ledger,
provenance, C++ issue log, crate register, benchmarks and documentation) follow
on `integrate/wave3-shared`. [VALIDATION](VALIDATION.md) records each package's
evidence, its verifier's reruns and this pass's gates;
[BENCHMARKS](BENCHMARKS.md) is new.

**Merged** (branch commit, then merge commit where it differs):

| Package | Branch | Merge | Outcome |
|---|---|---|---|
| B4-GAUSS (fix round 3) | `2fc61e5` | — | done; `TraceFitter.h` and `GaussTraceFitter.h` now have real ledger entries at `complete` |
| B5-EGH | `7d0c975` | — | done, rebased on B4; `EGHTraceFitter.h` `complete` |
| wave-3a scaffold | `8f0bb3e` | — | integrator-owned: the three tool registrations, their `[[bin]]` entries and `FileHandler::load_experiment_with_read_options` |
| B7-FFAP-FEATURES | `ba913da` | — | done; `run()` produces features end to end. It also applied the lead's B6 (a) and (b) decisions in `seeds.rs` and its tests, with the reports disclosing the out-of-scope edits |
| P3-PICKER-TOOL | final commit | `4c2806d` | done for the preview: in-memory mode; `-processOption lowmemory` is P4 |
| A5-FILEINFO-TOOL | `a4eb586` | — | done for stage 1; `-i`, `-d` and `-c` need A6 |
| P4-PICKER-LOWMEM | `4293aab` | `4293aab` | done: `-processOption lowmemory`, the source's `PPHiResMzMLConsumer` through `MzMLFile::transform`; the mode's divergences from the in-memory mode reproduced, including the source's dangling header references and the partial document a failing run leaves; the index it writes ported into `MSDataWritingConsumer`; `CPP-172` promoted to executed and widened, `CPP-339` and `CPP-340` found |
| A6-FILEINFO | `0510382` | `0510382` | done: `-i`, `-d` and `-c`, tier 1 on 59 executed Release-build cases of which 38 compare both reports byte for byte; `CPP-335`, `CPP-336` and `CPP-337` found |
| C5-FFC-WRAPPER | `ba67aa9` | — | superseded in part by `fix/ffc-integration` |
| crate/regex-facade | `2a29bd0` | — | done after six review rounds; the crate register row is closed |
| fix/mzml-reader-scale | `f03ef85` | — | done; size-derived reader allowances and the source's timestamp leniency |
| fix/mzml-writer-scale-parity | `522ce8f` | — | done; per-record writer budgets and `indexedmzML` by default |
| fix/tool-threads | `a442694` | — | done; `-threads` reaches a real pool and a non-positive count means every processor, as the executed C++ does |
| fix/baseline-filter-last-point | `4a597cf` | — | done; the one-sample-element end behaviour, and `erosion_simple`/`dilation_simple` now select the simple variants |
| fix/picker-scale | `6b55773` | — | done; input-derived acquisition ledger and `pick_experiment_in_place` |
| bundle/B3b-LM-FIDELITY | `dc56a9f` | `f7c9157` | done; the solver matches the Linux x86_64 Release Eigen bit for bit on all 141 traced fits, and the user took the platform decision |
| fix/ffc-integration | `4d53a7e` | `4091665` | done; the six FeatureFinderCentroided tool tests that B7 turned red are re-derived against the Release build |
| fix/picked-chromatogram | `e269586` | `b1700de` | done; the picked TIC chromatogram divergence was in the **mzML reader**, not in the picker, and the instrument-scale run is now bit-identical on both data and picked chromatogram |

**Nothing from wave 3 is in flight.** `fix/picked-chromatogram` was the last
lane out; the lead merged it as `b1700de` (with the tool-document correction
`864b295`) while this shared-file pass was under audit, and the pass was rebased
onto that tip. Its result closes the only open instrument-scale finding: the
picked TIC chromatogram of the 2.3 GB PeakPickerHiRes run had differed on 8,173
of 8,174 retention times and 7,891 of 8,174 intensities while all 22,776,198
spectrum centroids were bit-identical, and the cause was not
`PeakPickerHiRes::pick_chromatogram` at all. `MzMLHandlerHelper` applies the
minute multiplier of an `MS:1000595` time array in place, and for a 32-bit
array the element is a `float&`, so the `double` product is narrowed back to
`f32` (`CPP-306`); this input's TIC time array is exactly that case and the
picker's spline apex amplifies it.
`mzml::ReadOptions::source_time_array_precision`, set by `ReadOptions::source`,
reproduces the narrowing, and the library default keeps the precision.
[BENCHMARKS](BENCHMARKS.md) §3.7 now records the closed row, and the wave-4
run confirms it at full size on both thread counts.

**What B10-FFC-ACCEPT still has to close** (from the `fix/ffc-integration`
report):

1. **Tight comparison against the C1 oracle outputs themselves.** This branch
   pins FFC_1 against the in-repo retained expectation at 1e-9 relative on
   rt/`score_fit`/`score_correlation`. B10 clause 1 asks for 1e-9 against the
   C1 oracle output file, which lives outside the repository; it was measured
   (5.3e-7 worst, entirely the Debug writer's print precision; 2.2e-10 on
   `score_fit`) but not encoded as a test. B10 owns how that file enters the
   repository or the comparison harness.
2. **Seeds, asymmetric and debug 5 against the C1 outputs.** Measured here
   against the Release build (asymmetric exact to about 1e-15), but the in-repo
   assertions for those three modes are counts plus stdout, not a numeric
   comparison against a retained file, because no retained expectation exists.
   B10 clause 3 (24 seeds / 8 features, the `EGH_*` meta values, 1,054 hull
   points) is satisfied by measurement.
3. **Threads.** Clause 2 is satisfied and asserted in-repo
   (`the_output_is_byte_identical_at_every_thread_count`: `-threads` 1/2/4/8/0
   byte-identical, ids included, on a 384-core node and on macOS). What remains
   is pinning it against the C1 `FFC_1_threads_0/1/2/4/8` outputs, if B10 wants
   that.
4. **The LM budget boundary.** Untouched. `-algorithm:fit:max_iterations` 40
   differs from 500 and 50 equals it (symmetric); 6 differs and 8 equals
   (asymmetric), following C1's `FFC_max_iterations_symmetric_01..60` and
   `_asymmetric_01..10` sweep. None of those cases has been run on either side.
5. **The ledger.** Clause 5's wording is now applied:
   `FeatureFinderCentroided` is a validated TOPP workflow, and the tool stays
   `partial` for exactly one reason, the D5 FAIMS refusal (B11).
6. **New for B10: the zero-width-RT divergence.** The port refuses a zero-width
   retention-time range (`CPP-274`); the C++ Release build computes non-finite
   bin bounds and writes an empty feature map. The divergence is documented and
   pinned by an `#[ignore]`d test. Deciding it is an FFAP decision that
   surfaces through this tool.

**Decisions of 2026-09-15 applied in this pass:**

- `mass_trace:min_spectra = 1` follows the source (B6 (a), `CPP-271`) and a
  changed isotope abundance computes the intended override (B6 (b),
  `CPP-247`) — both landed in B7's branch, and the ledger and support documents
  now say so. B6 (c), the overall score, is still open.
- The dangling-reference leniency is a tool-side option (P2 request 1, D10):
  the library default stays strict, `FileHandler` does not enable it, and A5's
  FileInfo tool passes it on its own load path.
- `levenberg-marquardt` and `nalgebra` moved to `[dev-dependencies]` (B3
  request 2).
- `FeatureOverlapFilter.h` stays `complete`; `SignalToNoiseEstimatorMedian.h`
  stays `partial` while `AUTOMAXBYPERCENT` is refused. The wave-2 note that
  `GaussTraceFitter.h` and `EGHTraceFitter.h` should return to `unmapped` is
  superseded: B4 and B5 merged, so both have real review entries now.
- The numerics platform policy: match the Linux x86_64 Release build
  everywhere. Taken by the user after lane B3b reported, and recorded in
  [DISTRIBUTION_FITTERS_SUPPORT](DISTRIBUTION_FITTERS_SUPPORT.md) §1.

**Still open for the lead:**

- **Promote `FeatureFinderAlgorithmPicked.h` to `complete`?** Every public and
  protected member is ported except `writeFeatureDebugInfo_` and
  `abort_reasons_`, which are reachable only through `write_debug`, which the
  port refuses because the source throws there. This pass kept it `partial`,
  because decision 5 of 2026-09-15 keeps `SignalToNoiseEstimatorMedian.h`
  partial for the same shape of refusal. One rule should cover both.
- **A3 request 5** (the selected-ion drift time onto MS2) was decided in favour
  of the executed source, but no lane has been opened;
  `tests/file_info.rs::a4_mzml_file_1_all_flags` stays `#[ignore]`d until one
  is. **The lead.**
- **B6 (c)**, the overall seed score against the platform `powf` (`CPP-272`).
- **A resource-refusal exit code.** A ceiling refusal surfaces as exit 8
  ("Unexpected internal error") or, on the reader, exit 3
  (INPUT_FILE_CORRUPT). Neither reads as "out of resources", and one
  framework-level mapping would serve every tool. **CLI-1/CLI-2.**

**Carried forward, with owners in bold:**

- **Three manifests cannot be registered for source verification.**
  `tests/data/{feature_finder_picked,lm_eigen_path_differential,mzml_reader_scale}_provenance.json`
  (and, from wave 2, `lm_budget_differential`, `mzml_header_leniency` and
  `isotopes_source_precision`) use the key `target_verification` with a
  different meaning from the one `tools/check_core_sdk.py` reads (it expects
  `{revision, source_revisions, source_changes}` and finds a gate list or a
  free-form record), so adding them to `current_sdk_reference_manifests` would
  fail the checker. `baseline_filter_edges` and `mzml_mobility` additionally
  have `external_reference_artifacts` without an `external_reference_note`, and
  `baseline_filter_edges` without `origin_key` on its entries; the
  `file_info_text_format` fixture list names test-file anchors rather than
  paths. Their oracle artifacts are registered at `SOURCE_PROVENANCE.json`
  level with recomputed hashes, so nothing is unverifiable — but their `sources`
  are not re-hashed against the pinned SDK checkout. Rename the key (for
  example to `verification_runs`) and fill the two missing fields.
  **Each manifest's owning lane.**
- **A5's manifest was corrected by this pass**, which is an integrator edit of
  a lane file: `evidence_tier` was the bare number `1`, which the coverage
  generator cannot read as a tier (and crashed on), so it became the
  conventional descriptive string; and its "No case reaches a Debug-only
  precondition" sentence became "no compared *value* depends on one", because
  eight recorded stderr streams do carry the Debug-only "Update ranges was
  called but ranges were already up-to-date" assertion message. **A5/A6 to
  sign off.**
- **P3's six open verifier minors**, unchanged by the two follow-up lanes: a
  directory as `-in` exits 8 against C++'s 3 (a `FileHandler` follow-up); a
  non-mzML XML file under an `.mzML` name exits 3 against C++'s 11; the
  workflow-1 and C++-INI inputs make the C++ reader warn "Ill formed absolute
  or relative sourceFile path" where this port is silent, and that is not in
  the native-difference list; the C1 case `PPHR_cli_signal_to_noise_2_auto_levels`
  has no Rust-side test although it passes; and `render_list_parameters` is
  tested only for the empty list. `render_list_parameters` itself should move
  into `cli::processing::processing_info` or the mzML writer, because any tool
  with a list parameter is in the same position today. **P3/P4, and CLI-1/CLI-2
  for the shared rendering.**
- **A5's five open verifier minors:** the C++ loader warnings are not ported,
  so the tool's error stream is a subset of the C++ one on ordinary successful
  runs and two tests assert `err.is_empty()`; `c1_invalid_in_type_exits_6` uses
  `contains` and hides that the port emits the two diagnostic lines in the
  reverse order (a CLI-framework difference); `detect_type` is duplicated
  verbatim between `src/cli/tools/file_info.rs` and
  `src/format/file_info/report.rs`; `-in /dev/null` exits 2 where C++ exits 4
  (the `rustix::fs::access` candidate in the crate register); and `-write_ini
  <existing directory>` exits 8 where C++ exits 5, while the tool's own `-out`
  path agrees with C++. **A6, and CLI-1/CLI-2 for the last three.**
- **C5's five open verifier minors:** three C++-verbatim diagnostics are
  asserted against the port's own constants instead of literals; the tool
  citations (`To cite FeatureFinderCentroided:` plus Sturm 2010 and Weisser
  2013) are the first real `--help` divergence, which makes
  `docs/TOPP_CLI_SUPPORT.md:63` and its byte-for-byte help claim stale; the
  framework validates every registered input format before the tool body where
  the source checks `-seeds` late, so exit codes differ on an input that fails
  both; the wrapper writes the apex warning to stdout where the source writes
  it to stderr; and the module is the only tool module that widens the crate's
  public API, part of it with native wording B11 will change. **C5/B10, and
  CLI-1/CLI-2 for the help text.**
- **B7's four open verifier minors:** `extend_mass_traces` does not validate
  `pattern.spectrum` against `pattern.peak`, so an inconsistently built public
  `IsotopePattern` panics instead of returning an error;
  `preflight_seed_loop`'s per-seed estimate ignores the peaks inside the
  isotope window, so the ceiling under-refuses a dense adversarial input (real
  data is an order of magnitude inside it);
  `cropping_follows_the_source_position_rules` asserts `<= 1` where the source
  determines exactly 1; and `RunOutput::log` mixes the source's `LOG_INFO` and
  `LOG_WARN` streams (the blank lines were fixed by `fix/ffc-integration`).
  **B7/B10.**
- **`fix/picker-scale`'s two open minors:** `examples/peak_picking_scale.rs`
  `--ledger-probe` bounds its doubling search by a record count rather than by
  memory, so it peaked at 56.7 GiB against the 3,410 MiB of the run it
  measures; and `docs/PEAK_PICKING_SUPPORT.md`'s wall column was measured
  before the wave-2 merge and no longer describes the committed tree (every
  peak-RSS cell does reproduce to within 1 MiB). **P1/P4.**
- **`fix/tool-threads`' two library requests:** `Threads::from_cli` still maps a
  negative count to one worker and its doc still claims the source treats
  nonsensical values as single-threaded, which executed C++ disproves; and
  `map_collect`/`sum_in_order` build a fresh rayon pool per call, so calling
  them from inside `ToolContext::in_thread_pool` nests a second pool of the same
  size. No caller does that yet; the first wave-3 tool that parallelises will.
  **`src/concept/parallel.rs`'s owner.**
- **The acquisition-copy ledger's home.** `MorphologicalFilter::filter_experiment`
  meters a fresh `AcquisitionCopies` per spectrum while the trait default in
  `src/processing.rs:37` meters one for the whole experiment. Per-record is what
  real data needs (the 40,856-spectrum UK222 run fails with the shared ledger,
  and 200,000 three-peak spectra pass in 353 ms with the per-record one), and
  the per-record ledger still refuses a 60 MiB `source_file.name`. Apply the
  decision to the trait default instead of leaving one filter as the exception,
  and give `AcquisitionCopies` a named constructor
  (`AcquisitionCopies::for_records`) plus `#[derive(Debug)]` while doing it.
  **`src/processing.rs`'s owner.**
- **`SpectraFilterWindowMower`'s input-point ceiling**
  (`src/processing/window_mower.rs:69`) refuses 5,000 real spectra that the C++
  completes in 5.93 s. It is the last blocker for that tool's benchmark row.
  **Window-mower owner.**
- **Three stale "under investigation in lane B3b" references** remain in other
  lanes' files and are now wrong:
  `docs/EGH_TRACE_FITTER_SUPPORT.md:35`,
  `src/analysis/feature_finder_picked/gauss_trace_fitter.rs:102`,
  `src/analysis/feature_finder_picked/egh_trace_fitter.rs:85`, plus
  `src/analysis/feature_finder_picked/trace_fitter.rs:494` and softer phrasings
  in `tests/feature_finder_picked.rs:166`,
  `docs/EGH_TRACE_FITTER_SUPPORT.md:422` and this document. All should point at
  [DISTRIBUTION_FITTERS_SUPPORT](DISTRIBUTION_FITTERS_SUPPORT.md) §1 when their
  lanes are next opened. **B4/B5/B7.**
- **The benchmark harness needs three fixes before its next run**: the peak-RSS
  measurement (caveat 1 of [BENCHMARKS](BENCHMARKS.md)), data-and-metadata
  verdicts instead of stopping at the first difference, and one INI per tool
  for both implementations in `plans/full.json`, which today runs
  PeakPickerHiRes, FeatureFinderCentroided and FileInfo on each side's own
  defaults. **Benchmark lane.**
- **ibminode06's `/usr/local/bin/cc`** is a 2023 admin shell script that reports
  Ceph quotas and shadows the real compiler, so every cargo build on that node
  fails with a misleading "build-script-build (never executed)". All three of
  `CC=/usr/bin/gcc`, `CXX=/usr/bin/g++` and `-C linker=/usr/bin/gcc` are needed;
  the linker flag alone is not enough, because cc-rs then misdetects the script
  as MSVC. Ask the admins to rename it or add the three to
  `/scratch/kohlbach/openms-rs-env.sh`. Also: one agent overwrote that env file
  and the previous content is unknown, and ibminode05 carries an idle GitHub
  Actions runner that would disturb timing if a job landed on it.
  **Infrastructure.**
- **Only `PeakPickerHiRes` passes `mzml::ReadOptions::source()`.** `grep -rn
  "ReadOptions::source()" src/` matches `src/cli/tools/peak_picker_hi_res.rs`
  and nothing else, so `BaselineFilter`, `MzMLSplitter`, `MapNormalizer`,
  `SpectraFilterWindowMower`, `FileInfo`, `DTAExtractor` and
  `FeatureFinderCentroided` still read a 32-bit minute time array at full `f64`
  precision and will differ from their C++ counterparts by up to 2.44e-4 s per
  chromatogram point (`CPP-306`); for an I/O tool such as `MzMLSplitter` that
  difference lands in the written output. Raised by `fix/picked-chromatogram`
  and confirmed by its verifier; it is a cross-lane decision, not a defect of
  that commit, and the same question covers the two older switches
  (`source_dangling_references`, `source_invalid_timestamps`). The scaffold
  entry point exists: `FileHandler::load_experiment_with_read_options`. **Each
  tool's lane, with CLI-1/CLI-2 for a framework-level default.**
- **`fix/picked-chromatogram`'s two open verifier minors:** the `sec64` case of
  `tests/data/peak_picking/chromatogram_time/` carries the `f32`-narrowed
  seconds rather than the full-precision seconds, so it duplicates `sec32`'s
  physical values and is not the independent control the case-set sentence in
  `docs/PEAK_PICKING_SUPPORT.md` describes (the other six cases and every
  assertion are unaffected); and the commit message cites `MzMLFile.cpp:170`
  for the writer precision on the `store` path, where `store` actually reaches
  it through `XMLFile::save_` (`XMLFile.cpp:362` and `:369`) and `:170` is
  `storeBuffer` — `CPP-307` in this pass cites all three. **P1/P4.**
- **The three wave-3a tools never opted into the `-threads` pool.**
  `fix/tool-threads` request 6 asked each of `PeakPickerHiRes`, `FileInfo` and
  `FeatureFinderCentroided` to wire `Tool::run` as
  `ctx.in_thread_pool(|| Self::run_in_pool(ctx))?` and to pass
  `ctx.thread_policy()` to its parallel algorithms. Only the second half
  happened: `grep -rn in_thread_pool src/cli/tools/` matches the five wave-2
  tools and none of the three, `PeakPickerHiRes` and `FeatureFinderCentroided`
  take the policy into their algorithms without wrapping their bodies, and
  `FileInfo` uses neither. `tests/topp_threads.rs`'s `TOOLS` list is still the
  five, so nothing tests the three. Nothing is computed wrongly — all three are
  serial and byte-identical at every thread count — but `-threads` sizes no
  pool for them, and the first tool that parallelises will need the wiring
  anyway. Found while fixing the CLI document below, not by a lane.
  **P3/P4, A6 and B10, one line each.**
- **`docs/TOPP_CLI_SUPPORT.md`'s `-threads` paragraph was stale for a whole
  pass** and is fixed here: `fix/tool-threads` asked for it by name in round 1
  (request 3, `wf_126d6fd5-e02`), the lane landed the new contract and its own
  `docs/TOPP_THREADS_SUPPORT.md`, and the first shared-file pass recorded
  neither, leaving the shared CLI document describing the superseded mapping and
  the 169-line support document reachable from no shared file. Two process
  points follow, both **the lead's**: an integrator request naming an
  integrator-owned file needs a checklist entry of its own, and a new
  lane-owned `docs/*.md` needs a link from a shared document in the same pass —
  `tools/check_doc_coverage.py` measures rustdoc coverage, not whether a
  Markdown file is reachable.

## Wave 4 status

Status on 2026-09-16. `main` is at `9a392fe`, **26 commits above
`origin/main` `fabd4b9` and unpushed**. Unlike waves 1 to 3 this window adds
almost no new surface: it is performance on real data, and correctness on the
inputs the port previously refused. The shared files (CI, ledger, provenance,
C++ issue log, benchmarks and documentation) follow on
`integrate/wave4-shared`. [VALIDATION](VALIDATION.md) records each lane's
evidence, its verifier's reruns and this pass's gates;
[BENCHMARKS](BENCHMARKS.md) is rewritten around the wave-4 run.

**Merged** (branch commit, then merge commit where it differs):

| Lane | Branch | Merge | Outcome |
|---|---|---|---|
| `perf/peak-picker` (3 review rounds) | `6a278f6` | `2491afc` | done; the picking loop is parallel behind `parallel`, 1.66x at 32 workers and byte-identical at every worker count, peak RSS −886 MB, `-threads` wired through the pool, and the picker is the sixth tool in `tests/topp_threads.rs` |
| `perf/mzml-reader` | `5f762aa` | `821e783` | done; −47.2 % of the load path's instructions, no decoded value changed. The port's load of a 1.2 GB mzML is now **faster than the C++ reader's** |
| `perf/spline-scratch` | `6b43790` | `e766311` | done; `CubicSpline2dFitter` replaces eight heap vectors per constructed spline, bit-identical by construction |
| `perf/validation` (3 rounds) | `58dbff1` | `4389b96` | done; the reader's per-peak validation loop removed as **dead**, with the argument in the rustdoc and a 36-probe mutation-checked test. It also produced this window's measurement hazard, below |
| `fix/map-normalizer` (2 rounds) | `663d78b` | `c2ecade` | done; normalised against the spectrum maximum where the source uses the combined one (13.24x on 7,302 of 87,492 arrays), and the empty-range refusal restored and **executed** against the C++ |
| `fix/tool-limits` | `2f173d9` | `5856c63` | done; `MzMLSplitter` and `SpectraFilterWindowMower` process full-size input |
| `fix/featurexml-scale` | `82d13af` | `23b6519` | done; the 12.5 MB decode ceiling replaced by size-derived ceilings, both benchmark featureXML maps load |
| `fix/dta-precision` | `943a5de` | `db00dbf` | done; `DTAExtractor` output is byte-identical to the C++ tool's over 36,443 files and 836,505,793 bytes |
| `fix/ffc-integration`, `fix/picked-chromatogram` | — | wave 3 | recorded in wave 3; both are confirmed at full size by the wave-4 benchmark |

**Nothing from wave 4 is in flight.** The last lane out was `perf/validation`;
its third round restated every figure against the current merge base after the
reviewer showed that one arm had not been rebuilt in the same batch.

### What remains

- **B10-FFC-ACCEPT.** Unchanged and still the largest open item; its five
  clauses are listed under wave 3. The wave-4 benchmark adds one fact to
  clause 5: `FeatureFinderCentroided` **cannot** be measured at full size —
  neither implementation finishes the 43,745-spectrum run, the C++ was killed at
  a 600 s cap and an independent Rust re-run at 700 s was killed with no output.
  Any acceptance statement about that tool is a statement about a 4,000-spectrum
  subset until someone budgets an hour per cell. **B10.**
- **A6.** Unchanged: `-i`, `-d` and `-c` for `FileInfo`, and the
  `file_info::Options` strict default that keeps `c1_empty_mzml_mps` ignored
  although the source-compatible option now exists and A5's tool passes it.
  **A6.**
- **The FAIMS closure.** Unchanged: `FeatureFinderCentroided` is still
  `partial` for exactly one reason, the D5 FAIMS refusal, and
  `tests/topp_feature_finder_centroided.rs` still carries the zero-width RT
  divergence as its one ignored test — now with the C++ side logged as executed
  in `CPP-312`. **B10/C5.**
- **The bundle.** Eight tools are built and agree with the C++ Release build on
  the data of every one — seven of them measured on full-size instrument data,
  `FeatureFinderCentroided` on the documented 4,000-spectrum subset, because
  neither implementation finishes the full 43,745-spectrum run. What the bundle still
  lacks is not a tool: `-processOption lowmemory` (P4) and the three A6 FileInfo
  flags are ported as of wave 7, which leaves the acceptance statement B10 owns.

### Carried forward, with owners

New in this window, or restated because this window changed them:

- **THE COMPILER CODE-GENERATION QUANTUM — a measurement hazard for every
  future comparison.** `Record::finish` in `src/format/mzml.rs` has exactly two
  code-generation states on the PeakPickerHiRes workload, 88,619,093 and
  68,006,957 instructions of self cost, a quantum of **20,612,136**, and it
  flips between them on source perturbations that have nothing to do with it:
  the two `main` commits `e766311` and `8889ece` sit in different states and the
  only source file that differs between them is
  `src/cli/tools/map_normalizer.rs`, which the picker never calls.
  `MSSpectrum::validate`'s self cost is identical in both, so this is not an
  inlining transfer between those two symbols. Three rules follow, and they cost
  this lane a whole measurement cycle each: **(1)** a per-function instruction
  diff is not attribution — the same lane's `check_sorted` fusion showed −20.6 M
  in `Record::finish` and −5.3 M in `prepare_spectrum`, neither of which calls
  it, so only an ablation built on each arm attributes a saving; **(2)** every
  arm must be rebuilt in the same batch — reusing a binary from an earlier
  session moved a headline by 1.08 %, which the reviewer demonstrated by
  rebuilding one arm from `git archive`; **(3)** a layout control is mandatory
  and wall clock resolves nothing small here — a control with two unrelated
  functions swapped in source order executes the *identical* instruction count,
  yet a plain rebuild of an identical tree is worth a couple of tenths of a
  second either way and the two `main` commits sit 0.100 s apart in the
  **opposite** direction to their instruction counts. Nothing under about 1 s is
  resolvable on a loaded node, about 0.2 s on a quiet one. **Every future perf
  lane, and the benchmark lane.**
- **`ToolContext::in_thread_pool` at one worker.** Its doc still says a pool is
  built even for one worker so the body runs on a pool thread, which is what
  bounds a stray `par_iter` — and that costs **0.68 s** on a gigabyte-scale body
  through glibc's per-thread arenas. `PeakPickerHiRes` calls it
  (`peak_picker_hi_res.rs:249`) but scopes it to the picking call rather than
  the tool body, and skips it entirely at one worker. Decide once for the
  framework: a documented one-worker fast path (accepting that a stray
  global-pool `par_iter` would then be unbounded), or a note that a tool scoping
  its pool to its parallel region may skip it. Raised in all three picker
  rounds. **CLI-1/CLI-2 owner.**
- **`concept::parallel::map_collect` builds a `ThreadPoolBuilder`
  unconditionally**, with no ambient-pool check and no reuse across calls, which
  is why the picker could not call it and carries its own `BatchWorkers` in
  `src/processing/peak_picking.rs`. Moving `BatchWorkers` into the shared helper
  would give `src/comparison.rs` and
  `src/analysis/feature_finder_picked/algorithm.rs` the same behaviour.
  **Concept/parallel owner.**
- **`SummaryLimits::default().max_work` is a fixed 50,000,000** while
  `combined_ranges`, `chromatogram_ranges` and `calculate_tic_binned` charge one
  unit per spectrum, per peak, per chromatogram and per chromatogram point, so
  the 1.2 GB benchmark run charges **88,521,983** and `combined_ranges()`
  refuses a file the port's own mzML reader reads without complaint. Same class
  as the window-mower cap and the featureXML ceiling this window fixed.
  `MapNormalizer::range_limits` is a local workaround and should be deleted once
  the shared default is size-derived in the `src/format/mzml_scaling.rs` style.
  **`src/kernel/experiment_summary.rs` owner.**
- **Two remaining full scans on the picker's own path**, both measured and
  neither this window's to take: the experiment-level `input.validate()` in
  `start_experiment` (`peak_picking.rs:1631`) at 47,230,109 instructions, 1.24 %
  of the program and 16.0 Ir per peak — the largest single scan left — and
  `validate_points` (`:140`), still three passes of which
  `PickingCompatibility::source()` leaves one. The first cannot simply be
  deleted, because `MSExperiment` has public fields and a library caller can
  hand the picker anything; the technique that worked for the reader applies —
  a `validate_given_finite_peaks`-style entry point for the tool path, **with
  the argument written down** the way the reader's was. The second, if fused,
  must be **measured**: the same lane fused exactly this pattern in `kernel.rs`
  and it cost 1.0 s of wall despite executing 0.82 % fewer instructions.
  **Picker owner.**
- **`MSExperiment::max_intensity()`, or a doc cross-reference on
  `MSExperiment::ranges`.** `ranges(ms_level)` is spectra-only by design and the
  source's `updateRanges()`/`getMaxIntensity()` is not; nothing at the
  definition said so loudly enough to stop one tool making the substitution, and
  that substitution was a 13x wrong answer on real data. `grep -rn "\.ranges("
  src/cli/ src/bin/` matches one file today, so the window is small.
  **`src/kernel.rs` owner.**
- **A SIMD base64 decoder is the whole of what is left in the reader**, and it
  is a dependency decision, not a code one: `base64` 0.22's `GeneralPurpose`
  decoder is 330,022,649 Ir on a 49.7 MB slice against OpenMS's own
  `Base64::stringSimdDecoder_` at 162,210,512 — **2.03x**, and 16.5 % of
  everything the reader now executes, worth roughly −0.9 s on a 2.3 GB input.
  `base64-simd` is the obvious candidate; it is `no_std`-capable but uses
  internal `unsafe`, which `#![forbid(unsafe_code)]` permits in a dependency.
  Note it must also be in the nodes' `~/.cargo` registry cache, because the
  benchmark build is `--locked --offline`. **Crate-decision owner.**
- **The allocator question is smaller than it was, and the recommendation is
  "nothing in the crate".** On the current tree the untuned binary is 9.8 s
  faster, takes 38 % fewer minor faults and holds 940 MB less peak RSS than when
  it was first measured — purely from the merged reader and spline work — and
  its fault count is now level with the C++ reference while its peak RSS is
  487 MB below it. The glibc knobs are down from −1.78 s to −0.79 s. If they are
  wanted, set `MALLOC_MMAP_THRESHOLD_=4194304`,
  `MALLOC_TRIM_THRESHOLD_=67108864` and `MALLOC_TOP_PAD_=8388608` in whatever
  launches the tools, **all three or none**, and note they are glibc-only. Do
  **not** add `libc` and do **not** relax `#![forbid(unsafe_code)]` for
  `mallopt`. Separately, the earlier profiling lane's "`mimalloc` `LD_PRELOAD`
  gave −1.18 s" is **withdrawn**: that library defines none of
  `malloc`/`free`/`calloc`/`realloc`/`posix_memalign`, so under `LD_PRELOAD` it
  initialised and intercepted nothing, and three runs matched untuned glibc to
  72 minor faults in 1.85 million. Any future `mimalloc` decision needs an
  override build or a real `#[global_allocator]`. **Crate-decision owner.**
- **`tests/feature_finder_picked.rs` does not compile under
  `--no-default-features --features mzml,paramxml`**: it is gated
  `#![cfg(all(feature = "mzml", feature = "paramxml"))]` yet imports
  `openms::format::featurexml`, so that feature set fails the whole suite with
  `rc=101`. Pre-existing, hit by two lanes and by a reviewer in this window; it
  is why several lanes' no-default gates name targets explicitly instead of
  building everything. **B7/B10.**
- **The bit-identity gate needs no new definition, and the earlier request to
  change it is withdrawn.** The recipe is `-test`, which the port already
  implements for exactly this purpose (`src/cli/processing.rs:25-39`): it
  suppresses the wall-clock completion time, and with it the SHA-1 over the
  file. `bb13eecf…` reproduces from the merge base, the branch tip and a
  rebuilt tree, five runs. The replacement invariant an earlier pass proposed
  was itself path-dependent (its masked whole-file digest includes
  `parameter: out`), i.e. strictly worse than the gate it would have replaced.
  Recorded so it is not proposed again. **No owner needed.**
- **A `debug = 1` release profile for the benchmark harness.** The release
  builds carry no debug info, so callgrind gives function-level but not
  line-level attribution; that cost `perf/validation` a build-and-ablate cycle.
  **Benchmark lane.**
- **Infrastructure, re-confirmed by four lanes this window and still open:**
  `kernel.perf_event_paranoid = 4` cluster-wide, so no unprivileged `perf` and
  no hardware counters anywhere; `/usr/local/bin/cc` on ibminode06 is the Ceph
  quota script that shadows the compiler (see the wave-3 entry below — every
  lane that touched that node needed the same PATH shim); and `ibminode06` has
  no `~/.ssh/config` entry, unlike kim and dax. One correction to the wave-3
  wording: `ibminode06` is **not** unreachable — it answers SSH from ibminode05
  and rejects the key, and from the workstation the bare name does not resolve.
  **Infrastructure.**

## Wave 5 status

Status on 2026-09-17. `main` is at `59e0e1c`, pushed and green; this wave is
collected on `integrate/wave5`, which merges `port/ffap-complete` and
`port/signal-to-noise` (the third lane, `port/progress-logger-release-range`, is
already inside the first). The two branches share **no file**, so both merges
were conflict-free; everything integrator-owned was left to this pass.
[VALIDATION](VALIDATION.md) records the lanes, the six fix rounds and their
verdicts, the gates and the ignored-test inventory; [BENCHMARKS](BENCHMARKS.md)
§4 is the wave-5 measurement.

This is a completion wave, not a new-surface wave: it moves three headers from
`partial` to `complete` and closes the port's last documented divergences from
the Linux x86_64 Release build on the feature-finder path.

**Merged:**

| Lane | Branch | Merge | Outcome |
|---|---|---|---|
| `port/ffap-complete` (merge of `port/ffap-instrumentation` `c79e66f`, `port/ffap-semantics` `511d29e` and `port/progress-logger-release-range` `f89d5d4`; six combined fix rounds and a minors pass, each adversarially verified by two independent lenses) | `a11fc26` | `087ebd6` | done; `FEATUREFINDER/FeatureFinderAlgorithmPicked.h` **partial -> complete** |
| `port/signal-to-noise` (2 review rounds and a follow-up) | `be70a98` | `1a98104` | done; `SignalToNoiseEstimatorMedian.h` and `SignalToNoiseEstimator.h` both **partial -> complete** |
| `port/progress-logger-release-range` | `f89d5d4` | inside `087ebd6` | done; the Debug-only `OPENMS_PRECONDITION(begin <= end)` is no longer a refusal, and `CONCEPT/ProgressLogger.h` stays `native_equivalent` |

Review state across the ledger: complete 60 -> **63**, partial 62 -> **59**,
`native_equivalent` 90 unchanged, 786 registered public headers unchanged.

### Lead decisions of this wave (D1-D13)

These governed all six fix rounds and are the reason the three promotions are
defensible. They are also recorded in [VALIDATION](VALIDATION.md).

- **D1, the undefined-behaviour rule.** Reproduce the Linux x86_64 Release
  outcome when it is measured, repeatable, explained by the executed
  instructions and **in bounds**; refuse out-of-bounds reads or writes, data
  races, process termination and loops that never end. Round 4 applies it to the
  non-finite widths and meta values the source stores (`MetaValue::source_float`).
- **D2.** libstdc++'s binary-search probes on NaN keys.
- **D3.** Every `std::sort` as libstdc++'s introsort, and every
  `std::stable_sort` as libstdc++ 14.4.0's, NaN and equal keys included.
- **D4.** glibc's `-nan` in every ported `%g`/ostream formatter of FFAP and the
  trace fitters; the five formatters outside them are reported, not changed
  (carried forward below).
- **D5.** The reference build's glibc `powf`, ported licence-clean from Arm
  optimized-routines. It was reproduced exactly, so no fallback was needed.
- **D6.** `std::length_error`'s text above `vector::max_size()` in step 2.5, and
  the native ceiling below it.
- **D7.** `ChargedIndexSet` equality by index sets.
- **D8.** `RejectedParameters::Shown` stays the default, the heap-address bound
  is scoped to the reference platform, and macOS-only fit departures are platform
  notes. Since D10 only the EGH area's `atan` and the sign of a NaN the solver
  creates remain.
- **D9.** FAIMS is out of scope (a FeatureFinderCentroided tool decision, B11).
- **D10.** Every libm transcendental on the FFAP path follows the reference
  host's glibc: `exp` and `log` ported from Arm optimized-routines with the FMA
  fusion of `__ieee754_exp_fma` and `__ieee754_log_fma` (exact against 670
  million probed inputs on Linux x86_64 and macOS arm64); `atan` (`__atan_fma`,
  the IBM Accurate Mathematical Library, LGPL only) takes the fallback — the
  host's on x86_64 Linux with glibc, the `libm` crate's elsewhere, with a
  measured-maximum bound.
- **D11, the one accepted exception to D1.** The multi-thread race on `aborts_`,
  `abort_reasons_` and `log_` gives the **single-thread result**, because the
  determinism contract requires parallel output to equal serial output and a
  refusal would block essentially every parallel run.
- **D12.** Every `UInt`/`int` wrap that leads to an out-of-bounds write is
  refused exactly where it wraps, not through a raisable ceiling (the score-array
  count `3 + 2 * charge_count`); every wrap that stays in bounds is reproduced
  (step 1's `startProgress` range). Round 6 applies the same rule to the record
  of where the executed process ends: that line is a crate constant, not a
  `Limits` field. The round-6 minors then moved it from 1 GiB of arrays — read
  off runs under the oracle harness's 16 GB address-space cap — to the bytes of
  1,000,000,003 arrays, the largest count measured to die on the **uncapped**
  reference node.
- **D13** (after round 4). The crate-private `MetaValue::source_float` is
  accepted (the public API still refuses non-finite values); the featureXML
  writer's refusal of non-finite values and the CLI's "Unable to read file" text
  for a write failure are a separate task, and TOPP native difference 16 stays as
  recorded; the longer `feature_finder_picked` test time is accepted, with no
  assertion dropped; the charge counts `-2` and `-3` stay refused
  unconditionally; `atan` off x86_64 with glibc keeps the `libm` crate; and the
  step-3.3.5 termination was searched for once more and found.

### What remains

- **B10-FFC-ACCEPT.** Still the largest open item, and this wave narrows it
  rather than closing it. Its clause about `FeatureFinderAlgorithmPicked` being
  `partial` is now answered: the header is `complete`, the zero-width-range and
  short-input divergences are gone (the port follows the Release build), and
  `write_debug` is ported byte for byte. What is left is the **tool** level:
  TOPP native differences 14 and 16 (see below).
- **The `-C target-feature=+fma` build-flag question is SETTLED: on by default
  since `port/fma-default` (2026-09-18).**
  Measured on dax: the completed `FeatureFinderAlgorithmPicked` costs 21 % at one
  thread on a default x86-64 baseline build, the flag removes it and more
  (0.715x, 1.055x the C++ Release build at one thread and 0.894x at 32) with
  bitwise-identical output, and the control — main with the flag — gains only
  2.3 %. The flag also enables AVX and SSE4.2 crate-wide and raises the minimum
  CPU to Haswell/Piledriver and later. **Nothing was changed**: no `RUSTFLAGS`,
  no `[profile.release]` override, no `.cargo/config.toml`, on any branch or in
  CI. The alternatives are the flag, narrower `#[target_feature]` dispatch on the
  three ported replica functions (not measured), or accepting the 21 %.
  The user chose the flag. `.cargo/config.toml` sets it for
  `cfg(target_arch = "x86_64")` and a startup CPUID check refuses a processor
  without FMA with exit 12 and the rebuild command; output is unchanged
  ([FMA_BUILD_FLAG](FMA_BUILD_FLAG.md)).
- **TOPP native difference 16.** The featureXML writer refuses the non-finite
  feature values the C++ build writes as `inf`, so the port's
  `FeatureFinderCentroided` exits 3 without an output file where the C++ exits 0.
  Split off by D13 as `task_d0b10659`, together with the CLI reporting a **write**
  failure as "Unable to read file (parse error on line 0: …)" with exit 3.
  **featureXML and CLI owners.**
- **TOPP native difference 14.** At the port's own isotope-window ceiling the
  tool exits 8 with its own message where the Release build reports
  `std::bad_alloc` with exit 12. Memory-dependent, and deliberate under D6.
  **B10.**
- **The signal-to-noise signed-overflow refusals stay refusals.**
  `SignalToNoiseEstimator.h:123` and `SignalToNoiseEstimatorMedian.h:216`, `:228`,
  `:324` and `:365` each need more than `2^31` points per spectrum and 64-90 GB
  per evidence run, so emulating them is a stated cost/benefit decision, not an
  oversight. Reopening needs another evidence run. The beyond-`INT_MAX` full-path
  evidence lives outside CI in `../oracle/sne-fix/harness`. **Lead.**
- **The two consumer call sites of the noise estimator** —
  `PeakPickerChromatogram` (`src/processing/chromatogram.rs:181`) and
  `PeakPickerIterative` (`src/processing/iterative.rs:298`) — still call the
  strict `estimate` instead of the source-mode
  `estimate_peaks(..., &PickingCompatibility::source(), None)`.
  `PeakPickerIterative` has no source mode at all yet. Split off as
  `task_7dcbc815`. **Centroiding owner.**
- **`2^32 - 5` score arrays was not run uncapped.** On the measured model (about
  232 bytes per array, because the pattern loop names and `assign`s every
  in-bounds array first) it needs about 928 GiB, 93 % of the shared reference
  node's memory, so it was deliberately not run. That is why a recording ceiling
  exists at all. **No owner; documented.**

### Carried forward, with owners

- **CLOSED this wave:** the wave-4 item "`tests/feature_finder_picked.rs` does
  not compile under `--no-default-features --features mzml,paramxml`". The file
  is now gated `#![cfg(all(feature = "mzml", feature = "paramxml", feature =
  "featurexml"))]`, which matches what it imports, and the feature slices build
  clean. `tests/topp_threads.rs`'s unused picking constants on macOS and
  `kernel::validate_given_finite_peaks` under `--no-default-features` were fixed
  in the same round.
- **Five ported `%g`/ostream formatters still print `nan` for a NaN whose sign
  bit is set, where glibc prints `-nan`** (lead decision D4 applies only to FFAP
  and the trace fitters, which are fixed). They are
  `src/format/file_info/text_format.rs:615-623`,
  `src/math/posterior_error_probability.rs:1093-1096`,
  `src/format/mascot_generic.rs:908-911`, `src/format/pepxml.rs:3135-3137` and
  `src/format/sv_out_stream.rs:177-179`, `:208-210` and `:378`. The last one
  prints `NaN` and `nan`; check both against the source's spelling.
  `src/param/value.rs:612-619` and `src/processing/peak_picking/noise.rs`'s
  `stream_double` are already glibc-conformant. **Format owners.**
- **`tests/data_array_xml.rs` warns in two feature slices.**
  `cargo check --no-default-features --features mzml,paramxml --all-targets`
  reports `struct NoOutput` (`:16`) never constructed and `fn unsupported`
  (`:41`) never used; adding `featurexml` leaves the second. The file is gated
  `cfg(any(mzml, consensusxml, idxml, featurexml))` and compiles helpers only the
  idxml/consensusxml tests use. Not this wave's file and not changed here.
  **Owner of that test.**
- **CLOSED after the integration pass:
  `tools/check_core_sdk.py`'s "still present in repo" guard compared base names,
  not paths.** `check_external_artifacts` asserted
  `not (ROOT / Path(path).name).exists()`, so an oracle artifact called
  `Cargo.toml` or `LICENSE` was rejected because the crate has files of those
  names, which is why the integrator could not add
  `tests/data/signal_to_noise_provenance.json` to
  `current_sdk_reference_manifests` (the signal-to-noise lane's request R5) and
  reported the limitation instead of weakening the checker. The lead fixed the
  guard in `1c4f174`: it now compares the artifact's own path with the leading
  `../` dropped, and for a C++ artifact additionally refuses that file name
  anywhere in the crate, which is what the guard was written for in `82f2924`
  when the probes were moved out. Both halves are mutation-checked (a planted
  mirrored copy and a planted `.cpp` are each caught, and removing them restores
  exit 0). R5 is applied: the manifest is registered, and the added source
  references go 1,720 -> 1,728. `../oracle/ffap-complete-fix1/port-harness/Cargo.toml`
  and `../oracle/ffap-complete-fix3/upstream/LICENSE` are unblocked too.
- **`tests/data/topp_feature_finder_centroided_provenance.json` cites an absolute
  path.** It records the Release `FeatureFinderCentroided` binary as
  `/ceph/ibmi/abi/oliver/opt/…/bin/FeatureFinderCentroided`, which is not an
  `../oracle/` path, so that manifest can never enter
  `current_sdk_reference_manifests` as written. **TOPP FFC owner.**
- **`OpenMS_CPP_ISSUES.md`'s header and summary table have fallen behind the
  log.** The header says findings "were checked against OpenMS4-core
  `82ce5b3…`", while every entry from the CPP-2xx range on is checked against
  `bc9cc12`; and the summary table at the top stops at CPP-058, so 276 of the 334
  entries have no index row. Both are pre-existing and neither was changed here,
  because fixing them touches every entry. **C++ issue-log owner.**
- **Infrastructure, unchanged from wave 4:** `kernel.perf_event_paranoid = 4`
  cluster-wide; `/usr/local/bin/cc` on ibminode06 shadows the compiler; and
  `ibminode06` has no `~/.ssh/config` entry. Wave 5's oracle work ran on
  ibminode06 with the same PATH shim. **Infrastructure.**
- **`.reference/` is gitignored, so it exists only in the main checkout.**
  `tools/check_core_sdk.py --source .reference/openms4-core-bc9cc12` therefore
  fails from a worktree under `.claude/worktrees/`; pass the absolute path in the
  main checkout instead. Likewise `../oracle/` resolves against the **main**
  repository root, not a worktree's. **No owner; documented so the next
  integrator does not lose a cycle to it.**

## Wave 6 status

Status on 2026-09-18. `main` is at `e1c3115`, pushed and green; this wave is
collected on `integrate/wave6`, which merges `port/b11-faims` (`7921409`) and
`port/fma-default` (`86b0337`). The two branches share **no file** — 13 files
against 10, `comm -12` of their name-only diffs empty, re-checked at those two
heads — so both merges were conflict-free, and the merged tree changes exactly
23 files, which is 13 + 10. Everything integrator-owned was left to this pass.
[VALIDATION](VALIDATION.md) records the lanes, their verdicts, the guard's
evidence and its documented limits; [BENCHMARKS](BENCHMARKS.md) §4.5 records the
build-flag decision now that it is taken.

**B11 — the FeatureFinderCentroided FAIMS closure (`port/b11-faims`).** The
refusal of decision D5 is gone: the tool splits by compensation voltage (B8),
runs the algorithm once per voltage on that voltage's seeds, annotates each
feature with `FAIMS_CV` and merges across voltages (B9). Three defects lie on
that path and each is named at its item in
[TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT](TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md),
native difference 1: `CPP-278` cannot arise, because native ranges are computed
on demand; `CPP-282` is answered by assigning unique ids before the merge;
`CPP-283` by `FaimsMergeFidelity::Corrected`, the library option added beside
the faithful one, in the `AbundanceOverride` pattern. The corrected path has no
whole-tool C++ oracle — the C++ exits 8 on every FAIMS input — so each voltage
group was written as its own single-voltage mzML and run through the C++ Release
build, and the port's per-group features must equal that run (D6); the merge is
pinned against the derivation in *What the merge is meant to do* and
hand-derived numbers. 19 executed runs in `../oracle/b11-faims`, with a
`manifest.json`. The tool is no longer `partial` for the FAIMS reason.

One thing changed outside the FAIMS path: the tool's `OPENMS_LOG_WARN` lines now
go to **stderr**, where the executed C++ writes them (`faims_partial_cv`: the
skip warning and its `occurred 56 times` line are on the C++ stderr). This
closes the C5/B10 minor "the wrapper writes the apex warning to stdout where the
source writes it to stderr". No test asserted the old destination.

**The FMA build flag (`port/fma-default`).** The wave-5 open question is taken:
x86_64 builds set `-C target-feature=+fma` from the tracked `.cargo/config.toml`,
and `cli::run` refuses a processor without FMA with exit 12 and the exact rebuild
command rather than letting it die on `SIGILL`. No result moves;
[FMA_BUILD_FLAG](FMA_BUILD_FLAG.md) §7 measures that and says what would
invalidate it. `raw-cpuid =11.6.0` is the one new crate
([THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md)).

### What remains

Everything listed under *Wave 5 status → What remains* that this wave did not
touch still stands, in particular **B10-FFC-ACCEPT**, TOPP native differences 14
and 16, and the signal-to-noise refusals. New or changed here:

- **The cross-voltage merge has no C++ oracle and cannot get one from the pinned
  build.** Every per-voltage group is pinned against an executed C++ Release run,
  but the merge across voltages is pinned against a derivation and hand-derived
  numbers, because the C++ tool exits 8 before reaching it (`CPP-278`). The
  composite three-fix statement in `OpenMS_CPP_ISSUES.md` now says plainly that
  **no patched C++ build was ever executed**. Closing this needs a patched build,
  which is a decision about running modified C++, not a porting task. **Lead.**
- **The `is_x86_feature_detected!` trap is now live for every lane.** With the
  flag on by default, that macro and `cfg!(target_feature = …)` answer a question
  about the *build*, not the processor, for `fma`, `avx`, `sse3`, `ssse3`,
  `sse4.1` and `sse4.2` on `x86_64-unknown-linux-gnu`, and `-O` then deletes the
  guard written around them with no diagnostic. It had already, silently, folded
  `atan_is_reference()` in `glibc_libm.rs` to a constant `true`. Recorded as a
  standing hazard in [VALIDATION](VALIDATION.md) and in the porting skill.
  **Every lane owner.**
- **Test binaries are not guarded.** They do not go through `cli::run`, so a test
  binary on a processor without FMA dies on `SIGILL` with no message where a tool
  binary would print the requirement and exit 12. Accepted and documented; it is
  also the CI failure mode if a runner ever lacks FMA. **Documented, no owner.**
- **Pre-AVX processors get best effort, measured rather than guaranteed.** A
  separate baseline-built launcher is out of scope. **Lead.**
- **The benchmark nodes' shared `~/.cargo` cache needs `raw-cpuid 11.6.0`**
  before `build/build_rust_orig.sh` can keep building `--offline`. Being
  populated by the lead. **Lead.**
- **The `format` ↔ `system` module cycle forces a four-line duplicated CPUID
  read** in `src/analysis/feature_finder_picked/glibc_libm.rs`, because `analysis`
  may not name `crate::system` without closing a cycle that
  `tools/check_module_cycles.py` refuses. When that cycle is unpicked the
  duplication collapses into a single call to `cpu_features::cpu_provides_fma`.
  **Workspace-split owner.**
- **`tools/check_module_cycles.py` counts a module path written in a comment as a
  dependency**, because it scans raw file text with
  `re.findall(r"\bcrate::(\w+)")`. That is enough to fail the gate on
  documentation alone, and it did, twice, in the FMA lane. The heuristic is cheap
  and conservative, so this is a note, not a bug report; stripping `//` lines
  before the scan would be a two-line change if it ever bites a docstring nobody
  can reword. **Tooling owner.**

## Wave 7 status

Status on 2026-09-19. `main` is at `36c26a0`, pushed and green at 5,317 tests;
this wave is collected on `integrate/wave7`, which merges
`bench/wave6-refresh` (`3096196`), `port/a6-fileinfo` (`0510382`) and
`port/p4-lowmemory` (`4293aab`). All three merged without a conflict; the merged
tree changes 119 files against `main` before this pass's own records, and no
integrator-owned file was touched by any lane.
[VALIDATION](VALIDATION.md) records the three lanes, their verdicts, the two
majors applied in the merged branch and their re-measurement, the corrected
wave-4 memory basis, the lead decisions of this wave, the gates and the
ignored-test inventory; [BENCHMARKS](BENCHMARKS.md) §3 is the wave-6 refresh.

**A6 — the FileInfo `-i`, `-d` and `-c` checks (`port/a6-fileinfo`).** The three
sections of `FileInfo::report_` that inspect a file rather than summarise it are
ported, with a Release-build executed differential: 59 cases, 38 of them
compared on both reports byte for byte, run twice and reproduced. `-i` ends the
report where the source returns, so the tool reproduces `TOPP_FileInfo_11`'s
non-zero exit, and `TOPP_FileInfo_19` is reproduced verbatim. Two places where
the source is undefined are refused rather than guessed. Three C++ defects were
found and are filed as `CPP-335`, `CPP-336` and `CPP-337`.

**P4 — the low-memory picker (`port/p4-lowmemory`).** `-processOption lowmemory`
is ported through the port's own `MSDataWritingConsumer`, with the mode's
divergences from the in-memory mode reproduced where the source is defined and
recorded where it is not. 81.7 MiB resident against the in-memory mode's 3.29
GiB on a 2.3 GB run, with a byte-identical mzML body. `CPP-172` is promoted to
executed and widened; `CPP-339` and `CPP-340` are new.

**The benchmark refresh (`bench/wave6-refresh`).** One Markdown file, no
measurement changed in its closing round, nothing run on ibminode05. Its major
was a false memory claim, now replaced by a paragraph that prints wave 4 beside
wave 6 so a reader can check it.

### What this wave leaves for the next one

- **A7 and A8** are what `FileInfo.h` is still `partial` for. A6 closed `-i`,
  `-d` and `-c`; `-v` is the only refused flag left, and the consensusXML,
  identification, FASTA, mzXML, mzData and trafoXML branches are still open.
- **`-c` cannot report two classes of corruption**, because this port's reader
  and kernel refuse them before `-c` sees them
  (`src/format/mzml.rs:1038`, `src/kernel.rs:810-821`). The C++ loads both.
  **mzML reader / kernel owners.**
- ~~**The mzML reader's dangling-`sourceFileRef` strictness.**~~ **Closed in
  wave 8 as decision D14** (`fix/reader-round-trip`): `Registry::source` has the
  two-policy treatment, the reversal is recorded, and the round trip is executed
  both ways against the Release build.
- ~~**`SourceDangling` reproduces the source by rendered content where the
  source compares by pointer.**~~ **Closed in wave 8**, and not the way this
  bullet expected: no identifier had to be carried through the reader, because
  `Registry::processing` already shares one `Arc` per `dataProcessingRef`
  exactly as `processing_[ref]` shares one `shared_ptr`, so element-wise
  `Arc::ptr_eq` *is* the source's comparison. What remains is the whole-document
  writer, which still deduplicates by content on both sides; see wave 8 in
  `VALIDATION.md`. **mzML writer owner.**
- ~~**A mechanical citation check would pay for itself.**~~ **Built in wave 8**
  (`tools/citation-checker`): `tools/check_source_citations.py` resolves 3,344
  citations against the pins and reads them back, and it caught a defect A7
  shipped and three review rounds missed. Read its limits before reading a green
  run for more than it says — the A6 instance that motivated it is still **not**
  caught, because all eight places that write it paraphrase the source, and a
  paraphrase cannot be read back. The lever for the rest is a convention: quote
  the source verbatim in the span beside the citation. **Tooling owner.**
- **A robust suite count belongs next to the gate brief, not in each lane.** The
  ssh capture on the gate path drops lines, and both damage kinds matter: an
  eaten `... ok` line leaves the count low and the `test result:` line right,
  while an eaten `test result:` line leaves the window covering two binaries, so
  the count is right and the surviving result line is low by a whole binary.
  Walking the log and taking the larger of each pair repairs both; detecting the
  disagreement alone repairs neither. **Gate-brief owner.**
- **`/usr/local/bin/cc` on ibminode06 is an admin ceph-quota shell script**, not
  a compiler, and it shadows `/usr/bin/cc` on `PATH`. `rustc` links a build
  script through it, writes no binary and still exits 0, so cargo fails later
  with "could not execute process .../build-script-build (never executed)". Name
  the linker with `CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=/usr/bin/gcc`,
  **not** with `RUSTFLAGS`, which replaces the repository `.cargo/config.toml`
  entry and silently drops `-C target-feature=+fma`. **Gate-brief owner.**
- **A second gate host is worth assigning when lanes close in parallel.** One
  lane's first full-suite attempt exited 255 from an ssh drop mid-run while dax
  carried another lane's test gate at load 11.4. **Gate-brief owner.**

## Wave 8 status

Status on 2026-09-19. `main` is at `2ba9c1d`, pushed and green at 5,401 tests;
this wave is collected on `integrate/wave8`, which merges
`fix/reader-round-trip` (`42064c6`), `fix/featurexml-nonfinite` (`511fafa`),
`fix/picker-noise-consumers` (`8d3e052`), `port/a7-fileinfo` (`7801a05`) and
`tools/citation-checker` (`f9c692f`). Three of the five touch `src/format/`, so
disjointness was checked rather than assumed: the five name-only diffs against
`main` cover 237 paths and **no path appears in more than one lane**. All five
merged without a conflict, and the merged tree changes 248 files, which is those
237 plus the 11 this pass's own records touch. No lane touched an
integrator-owned file. [VALIDATION](VALIDATION.md) records the five lanes, their
verdicts, the A7 major applied in the merged branch and its re-measurement, the
two new lead decisions, the gates and the ignored-test inventory.

**A7 — the FileInfo consensusXML, identification and FASTA branches
(`port/a7-fileinfo`, merge `b2a173c`).** The three content branches A4 left open
are ported with their `-m`, `-p` and `-s` arms, against a Release-build executed
differential: 75 cases, 53 of them compared on both reports, run twice and
reproduced. libstdc++'s `std::hash<std::string>` is implemented because the
FASTA duplicate detection is sensitive to it. Seven C++ defects are filed,
`CPP-341` to `CPP-347`.

**The reader's round trip (`fix/reader-round-trip`, merge `bfd54a3`).** Decision
D14, and the pointer rule wave 7 left open, closed together — the second without
the reader change that bullet anticipated, because the identity was already in
the model.

**The featureXML writer (`fix/featurexml-nonfinite`, merge `008a09a`).** The
source's `inf`, `-inf` and `NaN` are written and read back, a failed store is a
write failure, and TOPP native difference 16 is closed.

**The picker consumers (`fix/picker-noise-consumers`, merge `c9d5b75`).**
Decision D15. `CPP-348` is new.

**The citation checker (`tools/citation-checker`, merge `ca960d1`).** A
repository checker that reads C++ citations back against the pins, plus a fix to
`check_core_sdk.py`, which had been answering "yes, C++ is in the tree" in the
one worktree that has the pinned checkouts.

### What this wave leaves for the next one

- **A8 is the last package `FileInfo.h` is `partial` for**, and it owes exactly
  four things now: `-v`, the pepXML, mzTab, trafoXML and PQP branches, the
  mzXML, mzData, MGF, MS2, sqMass, XMass and MSP peak files, and the schema
  validation A8 was scoped for. A6 closed `-i`, `-d` and `-c`; A7 closed
  consensusXML, idXML/mzIdentML and FASTA.
- **The x86_64-emulation promotion, carried forward with its three parts.**
  This is the lead's route-(b) decision of this round made concrete: the port
  spells every NaN the FileInfo text layer prints `nan` where glibc spells a
  sign-bit NaN `-nan` (native difference 5), and spelling the sign honestly
  requires the *value* to stop depending on the host first. Three parts, all
  prerequisites, in this order:
  1. promote the x86_64 emulation out of
     `analysis::feature_finder_picked::scoring` into shared math, where it can
     be reused;
  2. give `src/math/statistic_functions.rs` an x86_64-faithful
     `variance_with_mean` built on it, so the NaN a variance produces carries a
     host-independent sign;
  3. re-capture A2's oracle row against the **Linux Release build** rather than
     the macOS SDK — `../oracle/file-info-text-format/results/driver.tsv:92` is
     `D fff8000000000000 nan nan nan nan nan NaN`, measured with Apple libc, and
     `tests/file_info_text_format.rs:145` asserts it verbatim.

  Only after all three can `text_format`'s `nonfinite` rule change; until then
  changing it would contradict an executed oracle row and make every frozen
  expectation architecture-dependent. **Owner: shared math
  (`src/math/statistic_functions.rs`), with A2 as second party.**
- **A libstdc++-faithful `sort_ascending`, the other half of the same file.**
  `SummaryStatistics::new` reproduces the two NaN shapes whose set of outputs
  has one member and refuses a NaN next to a number — a **deferral**, not a D1
  refusal, because nothing is out of bounds and the values are stable. The same
  `std::sort` equivalence without a NaN is *accepted* and its divergence pinned
  (native difference 6, signed zeros). Closing either means porting libstdc++'s
  permutation into `sort_ascending`, which every `SummaryStatistics` caller
  consumes, and pinning it against a build that is free to change it. Evidence
  is already on disk: oracle cases `c_nan_one_s`, `c_nan_two_s`,
  `c_nan_then_finite_s`, `c_finite_then_nan_s` and `c_zero_swapped_s`, all
  frozen. **What the lead has to decide is whether reproducing an unspecified
  `std::sort` permutation is in scope at all**; refusing is what the port does
  today for the NaN half, and the signed-zero half shows that refusing
  everything would mean refusing inputs the Release build handles in bounds.
  **Owner: shared math.** Run this as one wave with the promotion above: both
  land in the same file, both are consumed by landed ports, and one wave can
  re-capture A2's oracle row once instead of twice.
- **The whole-document mzML writer still deduplicates by content**, where
  `MzMLFile::store` dangles. Reproducing it would emit references mzML forbids
  by default from every write path, with nothing to opt into, so it is the
  lead's call rather than a fidelity fix. **mzML writer owner.**
- **The identification-XML reader's modified-hit budget.** More than 14 modified
  peptide hits in one idXML is refused, because `AASequence::parse_with_budget`
  charges a per-modified-sequence preflight against one document-wide
  50-million budget. Two A7 oracle cases have no differential because of it, and
  the oracle already records what the C++ prints for both, so closing it needs
  no reference-build run. **Identification-XML reader owner.**
- **The citation checker's item 7** — a bare range under a bare file name — is
  the largest unchecked population left and the largest false-positive surface
  in the tool. It needs its own measurement pass over the unresolved ranges
  before a line of it is written. **Tooling owner.**
- **One gate-slot rule, learned twice this wave.** `openms-kim-gate.sh` rsyncs
  with `--delete-excluded`, which deletes the remote `target/`, so two batteries
  sharing one slot produce failures that look like code failures — "extern
  location for approx does not exist", a build-script panic on a missing
  `OUT_DIR`. One slot, one battery at a time; check `ps -ax | grep
  openms-kim-gate` before launching. **Gate-brief owner.**
