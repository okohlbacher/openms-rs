# Differential validation against executed C++

This document defines how ported behavior is compared against the C++ SDK, what
counts as evidence, and what may enter this repository. It replaces the implicit
convention under which most operation groups were validated by reading the C++
rather than running it.

## Why this exists

Through the first 32 commits, evidence for nearly every operation group was
transcribed class-test literals, independently derived expectations and Rust-only
invariant tests. Those detect transcription drift. They cannot detect a misread
algorithm: a systematic error in, say, an optimizer's update order or a smoothing
window would be reproduced identically in the port, in the fixture and in the
support document, and every test would pass.

Executed comparison against C++ existed only for six isolated leaf probes
(2,310 cases against extracted OpenMS translation units, plus 488 against
external LIBSVM). No algorithm with SDK dependency closure — no smoother,
picker, aligner, inference algorithm or feature finder — had ever been compared
against running C++. `SOURCE_PROVENANCE.json` now states that precisely.

## The oracle lives outside this repository

**No C++ is committed here.** The prebuilt C++ tree and every oracle driver live
outside the crate, under `../oracle/`:

| Path | Contents |
|---|---|
| `../oracle/datetime/` | The relocated DateTime probe sources and stub headers |
| `../oracle/drivers/<group>.cpp` | A small `main` linking the prebuilt `core-build` library, reading a JSON case file and printing canonical full-precision results |
| `../oracle/run.py` | Builds and runs a driver, emitting a fixture into `tests/data/` |

Committing oracle drivers would re-import the C++ build system that this port
exists to escape, and would make the crate's evidence depend on a toolchain the
crate deliberately does not require.

What enters this repository is only: the generated fixture, its SHA-256, the
driver's SHA-256, the compiler identification, and the C++ source revision — all
recorded in the group's `tests/data/*_provenance.json`. The evidence chain is the
hash, not the C++.

## Evidence tiers

Every provenance manifest declares which tier its expectations come from. The
existing four-tier vocabulary from `AGENTS.md` is retained and sharpened:

1. **Executed differential** — a driver linking the built C++ produced the values,
   with the driver and fixture hashed. The only tier that can falsify a misread
   algorithm.
2. **Executed probe** — an extracted or adapted translation unit was compiled and
   run. Weaker than tier 1 because the compiled unit is an adaptation: stub
   headers or a capture-only error adapter may stand in for SDK infrastructure.
3. **Source review** — expectations transcribed from C++ class-test literals or
   read off the implementation.
4. **Independently derived / Rust-only** — expectations computed from the
   specification, or native invariants and resource bounds with no C++ analogue.

A group may not claim tier 1 or 2 without a hashed artifact in its manifest.

## Comparison policy

**Byte equality is not the criterion.** Compressed mzML will never match
byte-for-byte, and requiring it would produce false failures that mask real ones.
`docs/REPOSITORY_ANALYSIS.md` already states this. Compare canonical forms:

| Quantity | Comparison |
|---|---|
| Peak and array values | Decoded to `f64`, compared elementwise within a declared per-quantity tolerance |
| m/z, RT, intensity | Relative tolerance, declared per group; absolute floor for values near zero |
| Counts, indices, charges, MS levels | Exact |
| Metadata maps, CV terms | Compared as sorted key/value sets, not in document order |
| Text formats (DTA, MS2, TSV) | Exact after newline normalisation |
| Element order within a spectrum | Exact — ordering is part of the contract |

Each group's support document states its tolerances and justifies anything looser
than exact. A tolerance is a claim about floating-point accumulation order, not a
budget for disagreement: if a tolerance has to be widened to make a test pass,
that is a finding, not a fix.

## Current state of the C++ reference

Verified on this machine at the time of writing:

- `core-build/lib/libOpenMS.dylib` links, so oracle drivers are buildable now.
- 710 compiled class-test binaries run; `DateTime_test`, `MSSpectrum_test` and
  `GaussFilter_test` pass.
- The 130 TOPP executables in `topp-build/bin` did not run when this section was
  first written. They rejected the first data row of a well-formed 131-row tool
  manifest (`share/openms4/tools/topp.tools.tsv`, 4 tab-separated fields, no
  duplicates) with "Invalid or duplicate tool package manifest entry", for every
  prefix tried. That observation is limited to `topp-build/bin`, and it remains a
  candidate entry for `OpenMS_CPP_ISSUES.md` once diagnosed.
- The product SDK (`../product-sdk`, a Debug build of core `4fdec46`) runs. Its
  TOPP binaries execute, and class-level drivers link its `libOpenMS` and
  `libOpenMSTestFramework.a`. Early-TOPP wave 1 executed FuzzyDiff, FileInfo,
  FeatureFinderCentroided, SpectraFilterWindowMower and DTAExtractor from it.
  Decision D7 in [the work packages](EARLY_TOPP_WORK_PACKAGES.md) accepts it as a
  development-time tier-1 oracle: outputs are labelled
  `oracle-generated (tier 1 executed differential)`, Debug-only precondition
  exits never become Rust expectations, and bitwise comparison holds only on
  macOS arm64.

Oracle drivers built against the product SDK on this host:

- Compile with `-ffp-contract=off`, as OpenMS's `cmake/compiler_flags.cmake`
  does, whenever a driver instantiates OpenMS or Eigen templates or replicates
  library arithmetic. AppleClang's default FMA contraction otherwise changes last
  bits, and the replica differs from `libOpenMS`.
- Configure with `Boost_DIR=/opt/homebrew/Cellar/boost/1.90.0_1/lib/cmake/Boost-1.90.0`
  (the SDK requires Boost 1.90.0 exactly; the Homebrew `opt` link is 1.92),
  `Arrow_DIR` and `Parquet_DIR` from `Cellar/apache-arrow/25.0.0_1` (25.0.0
  exactly), `CMAKE_FIND_FRAMEWORK=NEVER` with the Homebrew CURL paths (to avoid a
  stray `/Library/Frameworks/libcurl.framework`), and
  `OpenMP_ROOT=/opt/homebrew/opt/libomp`.

An earlier diagnosis attributed the TOPP launch failure to a stale Homebrew
abseil dylib. That was correct when observed and has since resolved through a
`re2` relink; the manifest rejection above is a separate and current issue.

## Scope

This policy governs new operation groups from this point, and the retroactive
validation of already-ported groups in descending order of consumer fan-out. It
does not retroactively invalidate existing work: source-review evidence remains
recorded as what it is, and is upgraded to tier 1 as drivers are written.
