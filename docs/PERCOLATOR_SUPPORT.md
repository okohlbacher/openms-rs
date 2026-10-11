# Percolator support

Three headers of `ANALYSIS/ID` at `bc9cc12`, with the vendored Percolator
engine the last one wraps:

| Header | Rust | Tests |
|---|---|---|
| `PercolatorTypes.h` (header only) | `src/analysis/id/percolator_types.rs` | `tests/percolator.rs` |
| `PercolatorFeatureSetHelper.h` / `.cpp` | `src/analysis/id/percolator_feature_set_helper.rs` | `tests/percolator_feature_set_helper.rs` |
| `Percolator.h` / `Percolator.cpp` / `PercolatorImpl.h` | `src/analysis/id/percolator.rs` | `tests/percolator.rs`, `tests/percolator_oracle.rs` |
| `thirdparty/percolator/*` (vendored, upstream `eb157f7`) | `src/analysis/id/percolator/{engine,pep,svm,training}.rs` | the above, plus unit tests in `pep.rs` and `engine.rs` |

Hashes, fixtures and the external oracle artefacts are in
`tests/data/percolator_provenance.json`. None of the modules needs a Cargo
feature; `tests/percolator_feature_set_helper.rs` needs `idxml` for its
fixtures.

## Third-party code

`Percolator.h` wraps a vendored subset of [Percolator](https://github.com/percolator/percolator)
(Copyright 2006-2012 Lukas Käll and the Percolator authors, Apache License
2.0) whose SVM solver is liblinear's TRON (Copyright (c) 2007-2015 The
LIBLINEAR Project, BSD 3-Clause), as integrated by J. Halloran. No crate
implements Percolator, its semi-supervised cross-validation, or the
I-spline/PAVA PEP estimators; a general SVM crate would not reproduce TRON's
iterates. The vendored algorithm is therefore translated, file by file, into
the four private submodules, each of which carries the upstream notice in its
header. The four BLAS-1 kernels and the dense products the vendored code runs
through Eigen are reproduced with the crate's existing Eigen kernel
transcriptions (`eigen_dot`, `eigen_squared_norm`, `eigen_gemv_row`,
`eigen_sum` in `src/math/fitters/levenberg_marquardt.rs`, made `pub(crate)`
for this) plus a transcription of the serial Eigen 5 GEMM for `X^T X`
(`pep.rs::eigen_xtx`). No dependency was added.

## API mapping

### `PercolatorTypes.h`

| C++ | Rust | Notes |
|---|---|---|
| `struct RescoreInput` | `RescoreInput` | All nine fields, same names; `std::vector<bool>`/`int`/`double` become `Vec<bool>`/`Vec<i32>`/`Vec<f64>`, `StringList` `Vec<String>`. |
| (none) | `RescoreInput::{n_rows, n_features, validate}` | Native: the documented row/column contract as a check (see "Native differences"). |
| `struct RescoreOutput` | `RescoreOutput` | `scores`, `q_values`, `peps`. |
| `struct PercolatorModel` | `PercolatorModel` | `format_version` (default 1), `feature_names`, `weights`, `normalizer_type`, `seed` (default 0) as `Default`. |
| (none) | `PercolatorModel::{validate, bias}` | Native. |

### `PercolatorFeatureSetHelper.h`

| C++ | Rust | Notes |
|---|---|---|
| `static void mergeMULTISEProteinIds(vector<ProteinIdentification>&, vector<ProteinIdentification>&)` | `merge_multise_protein_ids`, `merge_multise_protein_ids_at` | The `_at` form takes the creation time the source reads from `DateTime::now()`. |
| `static void addMSGFFeatures(PeptideIdentificationList&, StringList&)` | `add_msgf_features` | Infallible, as in the source. |
| `static void addXTANDEMFeatures(...)` | `add_xtandem_features` | `Result`. |
| `static void addCOMETFeatures(...)` | `add_comet_features` | `Result`. |
| `static void addANDESFeatures(...)` | `add_andes_features` | Takes `&[PeptideIdentification]`: the source never mutates the ids. |
| `static void addMASCOTFeatures(...)` | `add_mascot_features`, `add_mascot_features_source` | See the defect below. |
| `static void checkExtraFeatures(const vector<PeptideHit>&, StringList&)` | `check_extra_features` | Returns the removed names the source only logs. |
| `static void addMSFRAGGERFeatures(StringList&)` | `add_msfragger_features` | |
| protected `static double rescaleFragmentFeature_(double, int)` | `rescale_fragment_feature` | Public here; the source never calls it. |
| protected `static void assignDeltaScore_(vector<PeptideHit>&, const string&, const string&)` | `assign_delta_score` | Public here. |
| protected `struct lq_ProteinHit` | `accession_order` | The comparator as a function. |

### `Percolator.h`

| C++ | Rust | Notes |
|---|---|---|
| `Percolator()` | `Percolator::new`, `Default`, `Percolator::with_parameters` | Parameters live in the public `PercolatorParameters` (one field per source parameter). |
| `~Percolator()` | (drop) | The source resets process-wide statics here; the port has none. |
| deleted copy constructor / assignment | `Clone` | Allowed because no state is shared (see "Native differences"). |
| `void rescore(vector<PeptideIdentification>&, const StringList& = {})` | `rescore_peptide_ids` | Atomic. |
| `RescoreOutput rescore(const RescoreInput&)` | `rescore` | |
| `const vector<vector<double>>& getSvmWeights() const` | `svm_weights` | |
| `double getPi0() const` | `pi0` | |
| `PercolatorModel train(const RescoreInput&)` | `train` | |
| `RescoreOutput score(const RescoreInput&, const PercolatorModel&)` | `score` | |
| `static void fillPINCompatibleFields(const vector<PeptideIdentification>&, bool, RescoreInput&)` | `Percolator::fill_pin_compatible_fields` | Uses `format::percolator_infile::{scan_identifier, ScanNumberPattern}` for the `.pin` conventions. |
| `static void saveModel(const PercolatorModel&, const string&)` | `save_model`, `write_model` | Free functions (the source's are static). |
| `static PercolatorModel loadModel(const string&)` | `load_model`, `read_model` | |
| protected `void updateMembers_()` | `Percolator::from_param`, `from_param_with_warnings` | |
| `DefaultParamHandler` `getDefaults` / `getParameters` / `setParameters` | `Percolator::defaults`, `to_param`, `from_param` | The shape of `src/processing/peak_picking.rs`. All 17 parameters, with the source's descriptions, restrictions and valid strings, in declaration order. |
| private `struct Impl` | private fields of `Percolator` | |

The vendored classes are internal (`OpenMS::Internal::Percolator`) and are not
public API; their mapping is by file: `Scores`/`ScoreHolder`/`PseudoRandom`/
`PosteriorEstimator` (q-values, p-values, pi0)/`InProcessSampler` →
`engine.rs`; `InferPEP`/`PAVARegressor`/`ISplineTRRRegressor` → `pep.rs`;
`ssl`/`tron` → `svm.rs`; the normalizers, `SanityCheck` and `CrossValidation`
→ `training.rs`. Only the paths the wrapper reaches are ported:
`quickValidation`, per-iteration test reporting, RESET decoy subsampling
(`decoyFractionTraining < 1`, pseudo-targets), weight-file input, `pin.xml`
default weights, the `.pin` reader, protein inference, the XML/pepXML writers,
`weedOutRedundant` (peptide level) and `weedOutRedundantMixMax` are
unreachable from `Percolator.h` and absent.

## Preserved source conventions

- **Everything is bit-identical to the serial C++ Release build** on the 13
  oracle cases: SVM weights, pi0, scores, q-values and PEPs (see the evidence
  section). That required keeping, among others:
  - every `std::sort` through `math::source_sort::source_sort_by` (libstdc++
    introsort order for ties), every `std::unique` with its un-erased tail
    (`generatePositiveTrainingSet` with `train_best_positive` sorts and
    de-duplicates the fold's own score list without erasing, so later
    iterations see the duplicated tail; the port keeps that);
  - `calcScore` summing the features from the last to the first, but
    `scoreAndAddPSM` from the first to the last;
  - the process-wide LCG (`seed * 279470273 mod 4294967291`), reseeded at the
    start of `rescore` and `train` but **not** of `score`;
  - `OrderScanHash` with the 32-bit shift before widening;
  - `getInitDirection` resetting all but the bias entry of the direction, and
    `postIterationProcessing` rescoring every test fold with fold 0's weights
    after a failed direction check;
  - `normalizeScores` folding the rescaling into the fold weights before
    `getAvgWeights` averages them;
  - Eigen's reduction orders in TRON and the I-spline solver, and its serial
    GEMM for `X^T X` (panel depth `kc` from a 32 KiB L1, the even/odd split
    accumulators of the one-packet micro kernel, coefficient-based products
    below 20 combined dimensions, column blocks of 16 in the column-major
    GEMV from 128 columns on).
- The wrapper always sets the vendored `NO_TERMINATE`, so the vendored
  `MyException` paths are warnings: no targets or decoys, too good separation
  (pi0 = 1), no initial direction (first feature), median decoy above the FDR
  score (shift by one instead of rescaling), too few PEP bins (skip). The port
  follows them silently. The only vendored errors that remain are an
  unresolvable `initial_direction` and an out-of-range `std::vector::at`,
  which both become `Error::InvalidValue`.
- `pep_method = logistic_regression` (the default) is, despite its name, the
  I-spline fit of the target q-values (`calcPep(false, false, false)`);
  `isotonic` is PAVA.
- Rows lost to target-decoy competition are returned with score 0,
  q-value 1 and PEP 1.
- `rescore_peptide_ids` discovers features by sorted meta-value name, skips
  the documented blocklist, and uses identification index `i` as the
  cross-validation key of all its hits.
- The model file is the source's byte for byte (17 significant digits, the em
  dash in the second comment line); the reader is as strict as the source's,
  including `stoi`/`stod` prefix parsing.
- `PercolatorFeatureSetHelper`: feature registration order and names, the
  X!Tandem ion probe on the first hit of the first identification only, Comet
  `lnNumSP` falling back to the *unlogged* `MS:1002255`, `set_union` keeping
  the target run's hit for equal accessions, and `SE:<engine>` settings
  recorded once per engine.

## Native differences

- **State.** The source's process-wide statics (generator, feature count,
  normalizer, sanity-check direction) become fields of one `Percolator`, so
  instances are independent and `Send`. `score` continues from the generator
  state this instance's previous call left (seed 1 for a fresh instance); the
  source continues from whatever any Percolator in the process left. The two
  agree for one instance used alone, which is how the class test and the
  oracle use it.
- **Threads.** Results do not depend on `num_threads` (cross-validation
  candidates are trained in parallel and collected in order). The source's
  I-spline PEPs *do* depend on the OpenMP thread count: Eigen parallelises the
  `X^T X` product and then uses a different depth blocking (`kc` 320 instead
  of 304 for 600 targets). The port follows the serial C++
  (`OMP_NUM_THREADS=1`); `../oracle/I2-percolator/results/omp3` records the
  difference (7 of 13 cases, PEPs only, up to 1.6% relative).
- **Host.** Eigen derives the GEMM depth blocking from the L1 size it reads
  with CPUID; the port fixes 32 KiB, the reference host's value. On a host
  with another L1 size the source's PEPs for more than 504 targets differ from
  both.
- **Input checks.** `RescoreInput::validate` refuses non-finite features or
  masses and a non-empty `feature_names` of the wrong length; `rescore`,
  `train` and `score` also refuse zero feature columns and more than
  `MAX_INPUT_VALUES` values or `MAX_ROWS` rows. The source checks only the
  shape and lets NaN through the sorts.
- **Unported option.** `pep_method = nonparametric`
  (`PosteriorEstimator::estimatePEP`: a roughness-penalised logistic spline
  over `BaseSpline`, `PackedMatrix`, `PackedVector`) returns
  `Error::Unsupported`. It is the one remaining gap of `Percolator.h`.
- **Termination guards.** TRON's inner CG (at most `20 n + 1000` steps) and
  outer loop (100,000 steps) and the fold dealer (ten million draws) return
  `Error::InvalidValue` where the source would loop forever; only NaN in the
  data can reach them, and the input check excludes it.
- **Atomicity.** Every mutating entry point works on a copy: an error leaves
  identifications, runs and feature lists unchanged (the source throws
  halfway through its loops).
- **Model files.** `save_model` refuses what `loadModel` cannot read back: a
  feature name with a tab, a line break or a leading `#`, and a normalizer
  other than `stdv`/`uni`/`none`.
- **`PercolatorFeatureSetHelper`:** see the defect below for Mascot; empty
  X!Tandem input, an identification without hits and an empty sequence are
  errors where the source has undefined behaviour or divides by zero; a
  logarithm of zero or a division by zero in the Comet formulas stores the
  non-finite value, as the source, via `MetaValue::source_float`.

## C++ defect candidates

1. `PercolatorFeatureSetHelper::addMASCOTFeatures` computes
   `MASCOT:delta_score` and `MASCOT:hasMod` on a copy of each
   identification's hits (`std::vector<PeptideHit> hits = it->getHits();`)
   and discards it: only the sort survives, while the three features are
   registered. `add_mascot_features` writes them; `add_mascot_features_source`
   reproduces the source.
2. `Percolator` I-spline PEPs depend on `OMP_NUM_THREADS` through Eigen's
   parallel GEMM (above): the same input and seed give different PEPs on one
   and on three threads.
3. `Percolator.h` documents `rescore`/`train`/`score` as throwing on "too few
   decoys, no discriminative feature"; with `NO_TERMINATE` set by the wrapper
   itself they never do.

## Checked boundaries and evidence

- **Tier 1 (executed differential).** `tests/percolator_oracle.rs` replays 13
  cases produced by `../oracle/I2-percolator/driver.cpp` linked against the
  Linux x86_64 Release `libOpenMS.so` on ibminode06 (two serial runs,
  identical bytes): default, one thread, isotonic, `uni` with the Cpos/Cneg
  grid, `none`, nested cross-validation, three candidates per spectrum with
  and without TDC, reservoir subsampling, a negated initial direction,
  `use_pi0 = false`, the TRON fixture, and train-on-A/score-on-B. Contract:
  bit-identical weights, pi0, scores, q-values and PEPs, and byte-identical
  model files.
- **Tier 2 (executed probe).** `../oracle/I2-ispline-probe` compiles the
  vendored `ISplineTRRRegressor.h` with the `OpenMS_Percolator` target's
  Release flags (it reproduces libOpenMS's PEPs bit for bit) and runs it
  serially on ten q-value vectors, including 1, 5 and 17 targets (Eigen's
  coefficient-based product), 600 and 1500 (several depth panels) and 20,000
  (145 design columns, the blocked GEMV). `pep.rs` asserts every PEP bit for
  bit and Eigen's blocking sizes for seven shapes.
- **Tier 3 (class-test literals and fixtures).** Every section of
  `Percolator_test.cpp` is ported in `tests/percolator.rs`, with the source's
  glibc `rand()` data reproduced (`GlibcRand`, checked against glibc's first
  three outputs); the TRON literals at the source's 1e-8 relative tolerance.
  `Percolator_subprocess_parity_test.cpp` needs the external percolator binary
  and is not ported. `PercolatorFeatureSetHelper_test.cpp`: every section,
  with the MS-GF+, X!Tandem and Comet outputs compared value by value (not
  only by key, as the class test does) against the retained
  `*.topperc_check.idXML`.
- **Tier 4 (Rust-only).** Atomicity on error, parameter round trip, model file
  layout, `fill_pin_compatible_fields`, thread-count independence, the
  libstdc++ `unique` tail, the LCG recurrence.
