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

The three tools stay `partial` in the ledger until waves 5 and 6 close them;
FileInfo.h is `partial` since A4 and stays so until A6-A8 land.

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
| B3-LM | `levenberg-marquardt` adapter behind the unchanged `minimize`, exact Eigen `maxfev` emulation, budget differential | the pinned crates (`levenberg-marquardt =0.14.0`, `nalgebra =0.33.3`), already in `Cargo.toml` | first in the fitter lane if its gate passes; if not, D2's fallback, reported to the integrator |
| B4-GAUSS | TraceFitter driver, defaults and `Param` mapping; GaussTraceFitter | B1, C2, the scaffold's trait | before B5. If B3 merges later, B4's budget-boundary tests run again at B3's merge |
| B5-EGH | EGHTraceFitter | B1, C2, the scaffold's trait, B4's `optimize` | after B4, rebased on it |
| B6-FFAP-SEEDS | FeatureFinderAlgorithmPicked parameters, `run()` validation, scoring, pattern precalculation, seed selection | B1, B2, C2 | independent; `algorithm.rs` passes to B7 in wave 3 |
| A4-FILEINFO-CORE | FileInfo library: model, peak-file and featureXML branches, `-m/-p/-s`, text and TSV | A1, A2, A3, C1, C3 | independent; before A5 |
| CLI-2 | TOPPBase lifecycle part 2: input and output format checks, `-write_ini` parity, usage on stderr, DataProcessing retrofit (D4) | CLI-1, A3 | independent; before A5, C5 and P3 |
| P1-PICKER-LIB | PeakPickerHiRes and SignalToNoiseEstimatorMedian parameter contract and fidelity fixes | none | before P3 |
| P2-MZML-LENIENCY | Source-compatible dangling `softwareRef` and `defaultDataProcessingRef` (D10) | baseline | before P3; the integrator lands the FileHandler call site |
| B8-IMSPLIT | `IMDataConverter::splitByFAIMSCV` only (D5) | A1, A3 | independent; before B11 |
| B9-OVERLAP | FeatureOverlapFilter with its quadtree, source mode only (D5) | C2 | independent; before B11 |

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
  - the option leaves `sampleRef`, `defaultInstrumentConfigurationRef`, the
    `sourceFileRef` variants and scan `instrumentConfigurationRef` strict, where
    C++ is equally lenient; a follow-up if a tool input needs them; with the option
    on, the Rust tool prints warnings where C++ is silent (**P3**).
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
