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
| 1 | B3 Levenberg-Marquardt crate adapter with exact `maxfev` emulation | dependency commit; the digamma refactor lane |
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

- **Approved by their verifiers:** C2, C3 and A3.
- **A1:** approved; its commit is being checked against the reviewed hashes.
- **Fixes under way:**
  - CLI-1: the bare-invocation exit code and an unreadable `-ini`.
  - B1: provenance references that would shift unported headers in the ledger.
  - B2: the bit-for-bit claim is narrowed, because the C++ SDK iterates elements
    in a pointer-keyed map whose order changes between runs.
  - A2: exact decimal ties must round half to even, as `std::to_chars` does.
- **C1:** a full uninterrupted reproduction run is in progress.
- **B3:** starts after the digamma refactor lane merges.
