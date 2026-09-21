# What CI runs, what it costs, and the gate that keeps it from running less

Phase 2e, 2026-09-21, measured at `ab9eb55`. Three questions: can caching
make the Rust workflow cheaper, which integration tests never run under a
reduced feature set, and can the hand-kept `--test` lists be generated without
dropping anything. The answers are no, [154 on the minimum-Rust job](#2-the-reduced-feature-audit),
and [yes, proved](#3-the-generated-jobs-and-the-proof-that-they-drop-nothing).
What does make CI cheaper is [running its long lines side by side](#4-where-the-time-goes-and-what-shortens-it).

## 1. Caching: the dependency cache already works, and the workspace cannot be cached

The claim to test was that `Swatinem/rust-cache@v2` evicts the workspace's
own artifacts, so the ~228k-line crate and its 332 test binaries rebuild on
every run. It does, by design - and there is nothing to gain by stopping it.

**The dependency cache is already a complete hit.** In run 35565954724
(`1712a03`, all jobs green) the 67 cargo steps of `test`, `minimum-rust` and
`quality` compiled one crate between them, `openms`, once per feature set; no
dependency was rebuilt in any of them. Restoring took 5-9 s per job. The repository holds 12
cache entries of 126-379 MiB, 2.6 GiB of the 10 GiB quota, so nothing is being
evicted.

**The workspace's artifacts are invalid on nine pushes in ten.** Of the 63
pushes to `main` between 2026-09-10 and 2026-09-21, 56 changed `src/`,
`Cargo.toml`, `Cargo.lock` or `.cargo/`. The library is one crate and
rust-cache sets `CARGO_INCREMENTAL=0`, so each of those rebuilds it from
scratch - and every integration test depends on it, so every test binary
rebuilds with it. Only the other 7 (5 touching docs or tools, 2 only
`tests/`) could have reused anything.

**Keeping them anyway would cost more than it saves**, for three reasons read
out of the source of rust-cache v2.9.2 - the commit `@v2` resolved to in that
run, `6323deb` - rather than its README:

- *Its key has no source component.* It hashes the toolchain, the
  compiler-related environment and the parsed `Cargo.toml`/`Cargo.lock`, and
  on an exact hit it does not save at all ("Cache up-to-date", the last line of
  the `test` job above). With `cache-workspace-crates: true` the cached
  `openms` artifacts would be those of the first green run after each manifest
  change, never refreshed; every later run would download them and rebuild.
- *It would not keep the test binaries.* For a workspace member it keeps only
  the lib-like targets (`lib`, `rlib`, `dylib`, `cdylib`, `staticlib`,
  `proc-macro`) and deletes every other `deps/` entry, so all 332 integration
  test binaries rebuild regardless.
- *It would not fit.* A full sweep's `target/` on kim is 74 GB, 54 GB of it in
  `debug/deps`, with a 150-200 MB `libopenms` rlib per feature set; the whole
  repository quota is 10 GiB.

The other inputs named in the brief do not help either: `cache-all-crates`
keeps registry `.crate` files for crates outside the dependency graph (CI
tooling), not workspace artifacts; `shared-key` would put `test` and `quality`,
both on stable, under one key, where whichever saves first wins and the other
loses its artifacts; `save-if: main` changes nothing measurable (8 of 72 runs
were not on `main`, and there is no eviction to prevent); and
`cache-on-failure: true` would be harmful, because a cache saved from a failed
run holds only the dependencies of the steps that ran, and with exact-hit
semantics it would then stay the cache for that manifest until the manifest
changes.

**The one cache change made** follows from section 4: the legs of one job share
`GITHUB_JOB`, so each gets `key: ${{ matrix.leg }}`. Without it they would race
to save one entry and the loser's dependencies would never be cached. The first
run after this change starts cold, once per leg.

## 2. The reduced-feature audit

A test that uses a feature-gated path without carrying the gate's `#[cfg]`
compiles and passes under `--all-features` and breaks only in a build that has
the test's own gate on and the other feature off. That is how `f2bf2ed`
(`topp_file_info` under `mzml paramxml featurexml`, no `consensusxml`) and
`b1b52f7` (`statistic_functions` under no features at all) reached CI. So the
question is not whether a target is *named* on a reduced-feature line, but
whether it is ever compiled **non-empty** under one: a target whose crate-level
`#![cfg(...)]` fails on every reduced line it is on compiles to a binary with no
tests, and a missing `#[cfg]` inside it cannot show.

`python3 tools/ci_coverage.py --report` computes this from the workflow and the
test sources. At `ab9eb55`, for the 332 integration-test targets:

| job | non-empty under a reduced set | only as an empty binary | on no reduced line |
|---|---:|---:|---:|
| `test` (stable) | 285 | 47 | 0 |
| `minimum-rust` (1.85) | 178 | 0 | **154** |

The brief's "149 of 332 never exercised under any reduced feature set" is
332 minus the 183 names in `--test` lists. It misses that `test`'s bare
`--no-default-features` and `--features idxml` lines select every target, and
that five listed names appear only in `test`. The figures that matter are
above.

**The 154 on `minimum-rust`'s no line**: 17 have a crate-level feature gate, 19
have no gate but feature-conditional items or `cfg!` branches, 118 name no
feature at all.

- Crate-gated: `identification_pipeline`, `idxml`, `idxml_definitions`,
  `protein_workflow` (`idxml`); `mzml_auxiliary_review`,
  `mzml_auxiliary_writer_review`, `mzml_param_groups`, `mzml_review`,
  `on_disc_experiment` (`mzml`); `topp_baseline_filter`, `topp_dta_extractor`,
  `topp_map_normalizer`, `topp_mzml_splitter`,
  `topp_spectra_filter_window_mower` (`mzml paramxml`); `network`,
  `network_get_request`, `update_check` (`network`).
- Item-level gates: `adduct_workflow`, `chromatogram_workflow`,
  `emg_workflow`, `fine_isotope_stream_workflow`, `fine_isotope_workflow`,
  `iterative_workflow`, `precursor_workflow`, `rna_processing_workflow`,
  `rna_workflow`, `tagger_workflow` (`mzml`); `custom_modification_workflow`,
  `decoy_workflow`, `modification_generation_workflow`,
  `modified_peptides_reference`, `peptide_properties_workflow`,
  `sequence_workflow` (`idxml`); `spectrum_annotation_workflow` (`idxml`,
  `mzml`); `ribonucleotides`, `rna_reference` (`rna-json`).

**Eight targets were never compiled non-empty under any reduced set, in either
job** - the ones in which the defect of `f2bf2ed` could sit unseen:
`mzml_auxiliary_review`, `mzml_auxiliary_writer_review`, `mzml_param_groups`,
`mzml_review` and `on_disc_experiment` (gate `mzml`), and `network`,
`network_get_request` and `update_check` (gate `network`, a feature outside
`default` that only `--all-features` turned on).

**The 47 empty on `test`** are all crate-gated: 36 on `mzml`, the 3 on
`network`, `feature_finder_picked` and `feature_finder_picked_instrumentation`
(`mzml paramxml featurexml`), `feature_finder_picked_seeds` (`mzml paramxml`),
`msstats` (`consensusxml`), `numpress_coder` (`numpress`), `proforma_json`
(`proforma-json`), `qcml` (`paramxml`) and `transformation_xml` (`consensusxml`
or `featurexml`).
Stable runs them only under `--all-features`; after this change each of them
runs non-empty on `minimum-rust` (section 3).

**Whether the gaps hid anything**, measured on kim: every target under each
smallest feature set its gate admits, 17 lines on `+1.85.0` and the same 17 on
stable 1.96.0, as a generated job would run them. See section 3.4 for the
result.

## 3. The generated jobs, and the proof that they drop nothing

### 3.1 What is generated

`tools/ci_matrix.py` writes the `test` and `minimum-rust` jobs of
`.github/workflows/rust.yml` between two marker comments. Every cargo line
stays written out in the file, because `~/.local/bin/openms-ci-sweep.sh` reads
the lines it runs out of it; the rest of the file is edited by hand as before.

**`minimum-rust` derives its slices.** Every integration test runs there under
each smallest feature set its crate-level gate admits - `{}` for a test with no
gate, which is one bare `cargo test --locked --no-default-features`, and each
alternative of an `any(...)`. That is the build in which a missing `#[cfg]`
shows, and it is the build both incidents broke in. What a gate cannot say is
written down in `JOBS`, and nothing else is: the generator refuses a listed name
the derivation would add by itself, so each of the 28 names it lists, in 16
slices, is a choice somebody made -

- a test run under a larger set than its gate's smallest: `numpress_coder`
  (gate `numpress`) under `mzml`, `mzml_precursor_activation` (gate `mzml`)
  under `sqmass`, `file_info_checks` (gate `mzml`) and
  `feature_finder_picked_seeds` (gate `mzml paramxml`) under
  `mzml paramxml featurexml`, `cv_mapping_file` and `semantic_validator` under
  the validation features that contain their gate;
- a test with no gate run under the features its items, or the library paths
  it reaches, depend on: `peak_type_estimator`, `faims_helper`,
  `im_data_converter` and `record_metadata` under `mzml`, `file_handler` under
  `featurexml` and under `consensusxml`, `file_info_a7` under
  `consensusxml idxml`, `signal_to_noise` under `mzml paramxml`, and others;
- the dependency feature `rusqlite/extra_check` with `sqlite` and `sqmass`, and
  the serial/parallel split `mzml paramxml parallel`.

**`test` is written out whole**, line for line as it was: it is the stable
job, its slices are the stable counterpart of chosen `minimum-rust` ones, and
deriving there too would run every reduced slice twice for no class of defect
the minimum-Rust job does not already catch.

Lines with the same feature set are one line. `minimum-rust` had five
`--features mzml` lines and thirteen with no features; it has one of each.

### 3.2 The proof

`tools/ci_coverage.py` reads a workflow the way the sweep does - every `run:`
line of every step, leading `VAR=value` words stripped - and expands each cargo
line into the units it covers: (job, when it runs, runner, toolchain,
environment, command, feature set, target), with a line that names no target
covering every target cargo selects for it. `--superset OLD NEW` then compares:

```
$ git show ab9eb55:.github/workflows/rust.yml > /tmp/rust-ab9eb55.yml
$ python3 tools/ci_coverage.py --superset /tmp/rust-ab9eb55.yml .github/workflows/rust.yml
before: {"cargo lines": 83, "units": 3624, "(job, features, test) pairs": 2572, "(features, test) pairs": 1089}
after:  {"cargo lines": 64, "units": 3877, "(job, features, test) pairs": 2821, "(features, test) pairs": 1101}
superset: every one of the 3624 units before is covered after (253 units added)
```

"Cargo lines" counts each line once per runner it runs on, so the tag-only
`cross-platform` and `schema-minimum-platforms` jobs count twice. The lines the
pre-push sweep runs by default go from 67 to 48; the file's distinct cargo
lines from 75 to 56. Nothing that ran before stopped running:

| | before | after |
|---|---:|---:|
| `test` (job, features, test) pairs | 1037 | 1037 |
| `minimum-rust` (job, features, test) pairs | 529 | 778 |
| distinct (features, test) pairs, all jobs | 1089 | 1101 |

The 253 added units are all on `minimum-rust`: 228 integration tests plus the
library's unit tests, its binaries, examples and doctests under no features
(the bare line), and 21 tests under their smallest admitted set, 9 of which
`test` already ran under the same features on stable. The 12 new
(features, test) pairs are `network`, `network_get_request` and `update_check`
under `network`; `mzml_auxiliary_review`, `mzml_auxiliary_writer_review`,
`mzml_param_groups`, `mzml_review`, `on_disc_experiment` and `file_info_checks`
under `mzml`; `data_array_xml` under `consensusxml` and under `featurexml`; and
`feature_finder_picked_seeds` under `mzml paramxml`. After it, all 332 targets
run non-empty under a reduced set on `minimum-rust`, and no target is left that
no reduced line of either job compiles non-empty.

### 3.3 The gate that keeps it that way

A proof for one change says nothing about the next. The `quality` job now runs:

- `ci_matrix.py --check` - the generated jobs are what `JOBS` and the test
  gates produce; a new test file, or a changed gate, fails here until
  `ci_matrix.py --write` is run;
- `ci_coverage.py --check` - the workflow covers exactly the units recorded in
  `tools/ci_coverage_baseline.json`. Fewer is a drop; more is a unit a later
  change could drop unseen, so it is recorded (`ci_matrix.py --write` records
  additions itself). `--write` refuses to record a drop without
  `--allow-drop`, so running less is always a reviewed diff of that file;
- the same check refuses a cargo invocation the sweep cannot read (`cd x &&
  cargo test`), a job that runs cargo on every push but is not one of the four
  the sweep sweeps, and anything it does not model - a step `if:` other than a
  matrix leg, `continue-on-error`, a matrix `include`, a `cargo test` filter -
  rather than guess;
- `test_ci_matrix.py` and `test_ci_coverage.py`, which read the committed
  workflow back through the same parser and check that every derived
  (feature set, test) pair really is in `minimum-rust`, that every line runs in
  exactly one leg, and that the drops above are seen and the non-drops - a line
  moved to another leg, two lines folded into one - are not.

The coverage check needs PyYAML, as the sweep does; the `quality` job installs
`python3-yaml` rather than rely on the runner image having it.

**Adding a test** is now: add the file, run `python3 tools/ci_matrix.py
--write`, commit what it changes (the workflow and the baseline). The test lands
in `minimum-rust` under its gate's smallest feature set by itself. Only a slice
its gate cannot express - an item-level gate, a dependency feature, a TOPP tool
that should also run on stable - goes into `JOBS` by hand; the generated block
of the workflow is never edited directly, and `ci_matrix.py --check` fails if it
is.

### 3.4 Whether the added coverage is green today

Yes. On kim, at `ab9eb55`'s sources, every integration test under each smallest
feature set its gate admits - 17 lines, one per set, run on `+1.85.0` and again
on stable 1.96.0 - exited 0 on all 34, with 5368 tests passed and none failed
on each toolchain:

| smallest set | binaries | passed | `+1.85.0` | 1.96.0 |
|---|---:|---:|---:|---:|
| `{}` (bare: every target, lib, doctests) | 333 | 3673 | 153 s | 125 s |
| `mzml` | 40 | 850 | 112 s | 105 s |
| `mzml paramxml` | 11 | 201 | 107 s | 90 s |
| `featurexml mzml paramxml` | 4 | 140 | 327 s | 225 s |
| `idxml` | 9 | 151 | 29 s | 29 s |
| `paramxml` | 2 | 78 | 29 s | 26 s |
| `consensusxml` | 4 | 43 | 29 s | 24 s |
| `sqlite` | 2 | 39 | 24 s | 22 s |
| `featurexml` | 3 | 37 | 27 s | 24 s |
| `network` | 3 | 37 | 33 s | 22 s |
| `sqmass` | 1 | 33 | 27 s | 24 s |
| `mzml-schema` | 2 | 19 | 37 s | 30 s |
| `mzml-validation` | 1 | 16 | 27 s | 24 s |
| `semantic-validation` | 1 | 15 | 22 s | 20 s |
| `numpress` | 1 | 14 | 21 s | 20 s |
| `proforma-json` | 1 | 12 | 25 s | 25 s |
| `cv-mapping` | 1 | 10 | 21 s | 20 s |

So the gaps hid no defect today, including in the eight targets no reduced line
had ever compiled non-empty; what the derivation adds is that the next one
cannot hide. The generated `minimum-rust` job itself, line for line with
`cargo +1.85.0`, and the full pre-push sweep of the new workflow are recorded
in section 6.

## 4. Where the time goes, and what shortens it

From the same run, split at each step's `Finished` line:

| job | wall | compiling | running tests |
|---|---:|---:|---:|
| `test` | 31.6 min | 12.1 min | 19.0 min |
| `minimum-rust` | 50.9 min | 16.1 min | 34.0 min |
| `quality` | 2.0 min | 1.4 min | - |

Test execution, not compilation, is two thirds of it, and one binary dominates:
`feature_finder_picked` ran 5.3 min on stable under `--all-features`, 8.1 min
on 1.85 under `--all-features` and 8.6 min on 1.85 under
`mzml paramxml featurexml` - 22 runner-minutes, and 16.7 of the 51 on the
critical path. Making it cheaper is a change to the test or the test profile,
outside this lane.

What the workflow can do is stop running `minimum-rust`'s 43 lines one after
another. Each generated job is a matrix of legs that run in parallel, the lines
assigned by where their time goes:

| leg | lines | estimate from the run above |
|---|---|---:|
| `test` / `all-features` | `--all-features --all-targets`, `--doc` | 15.4 min |
| `test` / `slices` | the 18 reduced lines | 16.9 min |
| `minimum-rust` / `all-features` | `--all-features --all-targets`, `--doc` | 22.6 min |
| `minimum-rust` / `mzml-paramxml` | the three `mzml paramxml` slices | 15.9 min |
| `minimum-rust` / `slices` | the other 19 reduced lines, bare line included | about 17 min |

Setup - checkout, toolchain, cache restore, `libxml2` - is about 0.6 min a leg.
The workflow's wall time should fall from 51 min to about 23, bounded by
`minimum-rust`'s `--all-features` leg; runner time rises by about 6 min, three
legs' setup and the 253 added units. The legs use `fail-fast: false`, so a
failure in one no longer hides what the others would have found. The sweep
ignores `if:`, so it still runs every line of every leg, one after another.

These are estimates from one run's step times; the bare no-feature line is the
only line with no measured 1.85 counterpart on GitHub, and is estimated at
5.5 min from its stable counterpart scaled by the 1.85/stable ratios of the
`--all-features` line.

## 5. Two gates, checked before relying on them

- **The pre-push sweep does not test the minimum Rust.** It runs every line of
  `minimum-rust` with kim's default toolchain, rustc 1.96.0, because the line
  itself names no toolchain and the job's `dtolnay/rust-toolchain@1.85.0` step
  is not a `run:` line. So the sweep has never run a line of that job on
  1.85.0; the gate list's `cargo +1.85.0 check --all-features --all-targets`
  neither runs tests nor compiles any reduced feature set, and lanes have run
  chosen lines on `+1.85.0` by hand (docs/VALIDATION.md records several). The
  sweep lives outside this repository and is not changed here; section 3.4 ran
  every derived line on `+1.85.0` by hand.
- **The sweep skips any job not in its list** (`quality`,
  `portable-feature-graph`, `test`, `minimum-rust` by default). Splitting
  `minimum-rust` into new jobs would have silently dropped them from every
  pre-push sweep; that is why the legs are a matrix inside the existing jobs,
  and why `ci_coverage.py --check` now refuses a push job with cargo lines the
  sweep does not sweep.

## 6. Gates run on this change

- `python3 tools/ci_matrix.py --check`, `python3 tools/ci_coverage.py --check`,
  `tools/test_ci_matrix.py` (13 tests) and `tools/test_ci_coverage.py` (42
  tests): pass, on Python 3.14 and on the system Python 3.9. Each of twelve
  mutations of `ci_coverage.py` - legs ignored, the job condition or the runner
  dropped from a unit, a hidden cargo call or a test filter accepted, `all()`
  read as `any()`, implied features ignored, a job's `env` ignored, and others -
  fails at least one of them.
- `actionlint` 1.7.12 on the new workflow: no findings, as on the old one.
- `check_core_sdk.py`, `check_module_cycles.py` (60 edges, 8 mutually-dependent
  pairs), `check_source_citations.py` (3523 citations; 31 unreachable and 3 bare
  ranges, both unchanged) and `check_doc_coverage.py --report`: unchanged, since
  nothing under `src/` changed.
