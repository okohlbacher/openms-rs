# Validation of the ongoing Rust port

## SQLite S0 checkpoint (2026-09-13)

The [public SQLite connector](SQLITE_CONNECTOR_SUPPORT.md) was validated on
IBMI `kim`, using node-local `/scratch` sources and build targets, 32 compile
jobs, Rust 1.96.0 and minimum Rust 1.85.0. SQLite remains optional.

| Check | Final executed result |
|---|---|
| All features, current Rust (`nextest` and all-target `cargo test`) | 3,854 passed |
| All features/all targets, Rust 1.85 | 3,854 passed |
| No default features, current Rust | 2,640 passed |
| SQLite-only, both compilers, with and without `rusqlite/extra_check` | 19 passed in each of four selections |
| Doctests, each compiler | 26 passed |
| Formatting, full strict Clippy and rustdoc | passed |
| SQLite-only strict Clippy with `extra_check` | passed in the focused worker run |
| Provenance, coverage and feature-boundary gates | passed |

All eight upstream connector class-test sections and 17 assertion macros are
mapped to named native tests. The 19 native tests also cover binary bindings,
literal identifiers, external non-UTF-8 schema names, lock errors and recovery,
transaction cleanup, late row errors, ignored binding tails and dependency
feature unification. The feature-boundary check on Rust 1.85 verifies that
default/no-default builds exclude SQLite and libxml, and SQLite-only excludes
libxml while selecting the bundled SQLite dependency.

The initial full run passed before review refinements. Claude Fable 5.1 then
reviewed the connector and re-reviewed the changes. Byte comparison removed
an unnecessary UTF-8 failure, a stronger destructor test proved release of the
write lock, and bound execution was made independent of the dependency's
optional precheck. An overstatement in the first review was corrected against
the pinned dependency source. Both reports and their dispositions are
[retained](../tests/data/sqlite_connector_validation/reviews.json).

A separate adapted C++ probe reproduced CPP-184, CPP-185 and the CPP-189
exception-documentation mismatch using the exact pinned connector implementation
and declarations, substitute support headers and host SQLite 3.45.1. Its
sources, binary and logs are retained outside this repository and hashed in
[the provenance manifest](../tests/data/sqlite_connector_provenance.json).
This is isolated adapted execution, not a full SDK build, original exception ABI
test or execution of the upstream class-test binary.

[The validation record](../tests/data/sqlite_connector_validation/summary.json)
contains final runtime hashes, commands, outcomes and retained remote log paths.
The final full run has 19 successful checks; the initial run has 17. Strong local
verification checked 2,007 distinct pinned source/registration/reference files.
No new macOS or Windows execution is claimed. The next storage stages remain
unported; their source findings and fixture inspection are recorded separately
in the [S1 plan](SQLITE_STORAGE_PLAN.md), without claiming native fixes.


## FORMAT integration checkpoint (2026-09-12)

The fifteen-module [FORMAT wave](FORMAT_WAVE_SUPPORT.md) was built and tested
on IBMI `kim`, with sources and target directories on node-local NVMe `/scratch`.
Rust 1.96.0 and minimum Rust 1.85.0 were used with 32 compile jobs. Independent
workers used separate scratch directories; final checks used frozen snapshots.

| Check | Executed result |
|---|---|
| All features, current Rust (`nextest` and `cargo test --all-targets`) | 3,835 passed |
| No default features, current Rust | 2,640 passed |
| All features/all targets, Rust 1.85 | 3,835 passed |
| Doctests, each compiler | 26 passed |
| Six minimum-feature FORMAT selections on Rust 1.85 | all passed |
| Full strict Clippy, rustdoc and formatting | passed after fixes |
| Six provenance/coverage/feature-graph checks | passed |
| Ten scientific-data/reference regeneration checks | passed |

The first combined run found a module-order formatting change, three test-style
Clippy diagnostics and a DTA documentation link that needed qualification.
Those corrections were followed by green checks, another full current-Rust test
run, the affected minimum-feature tests and an MSRV all-target build check.
Runtime source and Rust-test hashes match the final tested snapshot; later edits
record evidence and documentation. Strong local source verification checked
2,004 distinct pinned source/registration/reference files.

Five completed Claude Fable 5.1 review reports and their dispositions are
[retained with the validation artifacts](../tests/data/format_wave_validation/reviews.json).
Review found a Percolator line-limit regression and helped check order-sensitive
pepXML chemistry. Unsupported review claims were corrected against source and
regressions; a model's approval is not an SDK-wide correctness proof.

[The machine-readable record](../tests/data/format_wave_validation/summary.json)
contains commands, outcomes, snapshot hashes and remote log locations. Full logs
remain under `/ceph/ibmi/abi/oliver/openms-rs/results/format-final-20260912-214524`
and `format-final-20260912-215205`. Initial failed quality checks are retained,
not erased. This checkpoint adds no executed full-SDK C++ differential and makes
no new macOS/Windows validation claim. Public API gaps remain in the coverage
ledger and per-format support documents.

## Historical throughput measurement before FORMAT (2026-09-12)

The full sweep on the remote host was profiled after wave 1 rather than tuned by
assumption. Splitting `cargo test --all-features --all-targets` showed a rebuild
after touching `src/lib.rs` costs 10 s while running the already-built tests
costs 47 s: **82% of the wall clock was test execution**, not compiling.

`cargo test` runs each test binary in turn. With 225 integration binaries the
per-binary serialisation dominates, and no amount of build parallelism touches
it. `cargo-nextest` runs every test from every binary in one work-stealing pool:

| Runner | kim (384 c) | Mac (16 c) |
| --- | --- | --- |
| `cargo test --all-features --all-targets` | 48 s | 56 s |
| `cargo nextest run --all-features` | **9 s** | **18 s** |
| `cargo test --no-default-features` | 38 s | — |
| `cargo nextest run --no-default-features` | **8 s** | — |

At that earlier, smaller snapshot, the complete sweep — build, both test
selections, doctests, clippy, the MSRV 1.85 gate, rustdoc, fmt and the six
Python gates — took **42 s** on kim. This is not a timing claim for the later
FORMAT integration or for a cold build.

Two things did not help and are recorded so they are not retried:

- **More build jobs.** The default thread count already reaches 7.7 s; forcing
  384 gives 7.7 s. Compilation is 10 s incrementally, so raising
  `CARGO_BUILD_JOBS` past the current 96 changes nothing measurable.
- **A RAM-backed `target/`.** `/dev/shm` offers 1.2 TB, but a cold build there
  took 29 s against 22 s on node-local NVMe. The node's 2.2 TB of RAM already
  page-caches `/scratch`, so tmpfs only adds a copy.

This corrects an earlier judgement in this project's own plan, which stated that
nextest should not be installed because "full-suite time is compile/link of 225
test binaries, not test execution". That was asserted without measurement and is
wrong by a factor of four on both hosts.

Doctests are not a nextest feature and keep `cargo test --locked --all-features
--doc`. CI continues to use `cargo test`, so the runner change affects local and
remote development loops only; every test still runs in both. nextest runs each
test in its own process, which is stricter than `cargo test`'s shared-process
threads, and all 2,506 tests pass under it.

## imzML family: the five-header format, and what four audits found (2026-09-12)

The user withdrew an earlier deferral of `KERNEL/OnDiscImzMLExperiment.h` and then
set a standing rule that **all file formats stay in core**. The imzML family was
ported as a staged DAG on that basis: the reader core, then the writer and the
kernel-level facade in parallel, then `ImzMLFile` carrying the family's whole
class-test suite. It is not five headers but nine — the facade pulls in
`MSImagingGeometry`, `MSImagingRegion` and `IonImage`, and a tenth,
`IonImageExtraction.h`, is outside the registered union and recorded separately.

| Package | Headers | Class-test sections | Rust |
| --- | --- | --- | --- |
| Reader core | ImzMLHandlerHelper, ImzMLHandler | own tests; family suite is stage 3's | `src/format/imzml_handler.rs` |
| Writer | ImzMLWriter | round-trip against the reader | `src/format/imzml_writer.rs` |
| On-disc facade | OnDiscImzMLExperiment + 3 imaging headers | 17 ported, 22 mapped | `src/kernel/on_disc_imzml_experiment.rs` |
| File adapter | ImzMLFile | **50 ported, 0 mapped, 0 unaccounted** | `src/format/imzml_file.rs` |

All four audits returned zero blockers and `bounds_enforced: true`, which was the
property that mattered: an imzML dataset is two files, and every array offset and
length comes out of the XML and indexes into the companion `.ibd`.

**One finding should have been a blocker and two auditors found it independently.**
`infer_ibd_path` tested its suffix by byte-slicing `text[text.len() - 6..]`, which
panics whenever the sixth-from-last byte is a UTF-8 continuation byte. Every public
entry point in the family calls it first — `load`, `load_experiment`,
`load_into_consumer`, `load_spectra_index`, `store` and the facade — so
`dir/日本語.txt` aborted the process before a single validation ran. Reproduced,
then fixed with a char-boundary-safe `str::get`. The same rewrite closed a second
divergence: the suffix is now replaced by truncation, matching the source's
`p.substr(0, p.size() - 6) + ".ibd"`, where `PathBuf::set_extension` had treated a
name that is entirely `.imzML` as an extensionless hidden file and appended.

**Eight further majors were fixed rather than carried.** The reader returned `Err`
for a spectrum with exactly one external peak array, where ImzMLHandler.cpp:198-234
fills the non-external side from the inline peaks and succeeds; `extract_ion_image`
read arrays that `spectrum()` refused for the same pixel, because one path tested
the externality flags and the other did not; and the writer's float cvParam text
used Rust's `Display` where the source uses `std::to_chars` with `chars_format::fixed`
precision 15 inside [1e-2, 1e4) and shortest-round-trip scientific outside it
(`NumericFormatting.h:26-135`), its six vocabulary meta keys errored where
`DataValue::toString()` is deliberately lenient, and a misaligned auxiliary array was
skipped under default options but fatal under any sort or trimming filter.

**The test oracle was weaker than it claimed, and the re-run cleared the port.**
The ported `ClassTest::isRealSimilar` omitted the opposite-sign branch of
ClassTest.cpp:439-451, so `close(-1.0, 1.0)` returned true where C++ returns false,
and it backed 38 assertions. Ported faithfully, **all 38 still pass** — the weak
oracle was not masking a defect. The same package checked the three sibling suites
that carry a similar helper (`binned_spectrum`, `feature_handle`, `rich_peak2d`) and
established that their `|1.0 - ratio| <= 1e-5` form rejects negative ratios, so they
never had the defect; that was independently re-derived here before accepting it.

**One bound was closed by the integrator because no package owned the file.**
`ImzMLFile::preflight` sums `mz_length` from the index, but an inline array carries
its length in the XML, so inline peaks reached the caller uncounted once the reader
began decoding them — bounded only by the 512 MiB XML ceiling, roughly 96 million
`f32` peaks. A `PeakBudget` now charges every decoded spectrum against
`max_loaded_peaks` in both load loops, failing as soon as the ceiling is crossed.

Nineteen further C++ defects are recorded as CPP-124 to CPP-142.

Gates on the Linux node at the integration commit: nextest all-features **3,193
passed**, no-default-features 2,376 passed, doctests 16 passed, clippy `-D warnings`
clean, `cargo +1.85.0 check --all-features --all-targets` clean, rustdoc
`-D warnings` clean, `cargo fmt --check` clean, all six Python gates green.

A correction to the preceding checkpoint: the sweep it cites for kernel wave 3 was
killed by an ssh disconnect after the all-features run, so its no-default-features,
clippy, MSRV, rustdoc and Python gate results were never obtained. Those gates are
verified here, on a tree that contains that work.

## Kernel wave 3 completion and its audit fixes (2026-09-12)

The three packages a session limit had killed were re-run from the integrated
base and merged: MSExperiment/AreaIterator residuals, the MSSpectrum ion-mobility
quartet with ConsensusFeature's last two members, and the OnDiscMSExperiment
facade. Kernel headers closed or native-equivalent rise from 27 of 34 to 30 of 34,
and `IMTypes.h` closed with a first entry for `SpectrumSettings.h` — the ion-mobility
quartet turned out to be declared in METADATA, not KERNEL.

Three audits returned, none with a blocker, and one reported
`section_audit_honest: false`. Every finding was fixed in a follow-up wave rather
than merged as-is. Three are worth recording because each was a claim the code did
not support.

**The test-mapping rule was circumvented.** The MSExperiment package reported "63
sections mapped with a cited asserted value". Its support doc accounted for all 54
mapped `MSExperiment_test.cpp` sections with a bare list of nine test *files* — no
function, no value. The auditor counted assertion macros per section and found 22
above the five-macro threshold that mandates porting, one of them with 64 macros;
sections 5 and 6 (copy and move assignment, asserting `getMinMZ 5.0`, `getMaxMZ
10.0` and a moved-from size of 0) had no Rust evidence anywhere in the repo. All 22
are now ported. The honest accounting is 56 ported, 21 mapped with a named function
and value, and **4 unaccounted** — `set2DData<add_mass_traces=true>`, both
`getFirstProductSpectrum` overloads and `operator<<`, whose members are unported.
`MSExperiment.h` stays `partial` for exactly those four.

**A rustdoc claim misdescribed the C++ it cited.** `area_iteration.rs` said a
reversed low/high ion-mobility pair "silently selects nothing" upstream.
`AreaIterator.h:277` builds `RangeMobility{low_im_, high_im_}`, and
`RangeBase(min,max)` (RangeManager.h:48-52) *throws* `InvalidRange` when `min > max`.
The port's `Err` agreed with the source by accident, not by the stated reasoning.
Several line anchors had drifted and an OpenMP note claimed the serial rasterizer
"computes the same image" — true only for `Max`, since the parallel branch merges
per-thread f32 buffers and so differs for `Sum`.

**An audit found a defect in this port's own mzML reader.** A fixture substitution
in the OnDisc package was documented as forced by an unavailable upstream file. The
file is committed in-tree, byte-identical to the pinned copy; what blocked it was
this reader rejecting `binaryDataArrayList count="2"` with four children. Checking
every list handler in `MzMLHandler.cpp` showed upstream **never** compares a declared
count against the actual number of children: `binaryDataArrayList` feeds only
`bin_data_.reserve` (:1015), `selectedIonList`'s count only warns when above one
(:1371), and `precursorList`, `productList`, `scanWindowList` and
`referenceableParamGroupList` have no list handler at all. This port hard-errored at
four sites. All are advisory on reading now, in both the reader and the `loadSize`
counting path, while every declared count remains a **resource ceiling** enforced
before allocation and writing still emits the true count. The earlier header-list fix
(`src/format/mzml_header.rs:107-116`) had addressed only one instance of this defect
class; this is the general case, and it had blocked real data twice.

The last item cost a deliberate reversal. The count fix initially kept
`referenceableParamGroupList` strict so `read` and `read_size` would agree, and pinned
that with a test. Relaxing only the reader would have left `read` accepting a document
`read_size` rejects, so both paths were relaxed together and the test rewritten to
assert that the two readers still agree — the property the strict check had existed to
protect.

Gates at the integration commit, on the 384-core node: build 21 s, nextest
all-features 2,917 passed, no-default-features 2,362 passed, doctests green, clippy
`-D warnings` clean, `cargo +1.85.0 check --all-features --all-targets` clean, rustdoc
`-D warnings` clean, `cargo fmt --check` clean, all six Python gates green. Four further
C++ defects recorded as CPP-120 to CPP-123.

## Kernel wave 3, partial: the map containers (2026-09-12)

Wave 3 launched four work packages. **Three were killed mid-run by an account
session limit** (WP7 MSExperiment/AreaIterator residuals, WP11 residual closure,
WP12b the OnDisc facade); they produced no commits and are queued for a clean
re-run. WP10 completed and is integrated here. A concurrently launched imzML wave
died the same way before its first stage committed.

| Work package | Headers | Class-test sections | Rust |
| --- | --- | --- | --- |
| Map containers | FeatureMap, ConsensusMap, ConversionHelper | 74 ported | `src/kernel/map_operations.rs`, `src/kernel/conversion_helper.rs` |

FeatureMap.h and ConsensusMap.h had never been reviewed member by member — both
sat at `evidence_requires_review` with no entry in the reviewed-API ledger at all.
All 74 upstream sections (32 + 39 + 3) are ported with transcribed literals, tier 3
evidence under [the differential validation policy](DIFFERENTIAL_VALIDATION.md);
none is merely mapped. Nine further C++ defects are recorded as CPP-111 to CPP-119.

**The audit found three documentation defects, all fixed before merge.** Each was
a claimed equivalence the code does not have: `isMapConsistent` was documented as
a plain mapping to `validate_consistency()`, which is in fact strictly stricter —
the source checks only duplicate column descriptions and unregistered handle map
indices, while the Rust additionally rejects a bad `experiment_type`, non-finite
coordinates, duplicate unique IDs and invalid attached records, so a map the
source calls consistent can return `Err`. Both `updateRanges()` rows omitted that
`ranges()` is fallible where the source cannot fail, because it opens with
`self.validate()?`. And the five `FeatureMap` sorts are stable here but use
`std::sort` upstream, which is not — a divergence the document was meticulous
about elsewhere and silent about here. The range *content* was verified faithful
in both maps.

**One finding was a defect in the integrator's own instructions.** The package
added four lines to `src/kernel.rs` where the rule allowed two. It registered two
modules, and the rule had assumed one module per package; the wording is corrected
for the re-run rather than charged against the package.

Gates, independently re-run by the auditor at the package commit: MSRV
`cargo +1.85.0 check --all-features --all-targets` clean, nextest all-features
2,825 passed, no-default-features 2,309 passed, clippy `-D warnings` clean,
rustdoc `-D warnings` clean, and 100% rustdoc on both new modules.

## Kernel wave 2: ion mobility, chromatogram merging, feature identification, MRM and indexed mzML (2026-09-12)

Five work packages ported in parallel git worktrees, each followed by an
**independent adversarial audit in its own worktree** that re-ran every gate
itself rather than trusting the port's report. Kernel headers closed or
native-equivalent rise from 18 of 34 to 24 of 34.

| Work package | Headers | Class-test sections | Rust |
| --- | --- | --- | --- |
| Ion mobility | MSSpectrum (mobility surface) | 68 ported, 3 mapped | `src/kernel/spectrum_mobility.rs` |
| Chromatogram merging | MSChromatogram, Mobilogram | 91 ported | `src/kernel/chromatogram_merge.rs` |
| Feature identification | BaseFeature, Feature, ConsensusFeature | 92 ported | `src/kernel/feature_identification.rs` |
| MRM | MRMFeature, MRMTransitionGroup | 44 ported | `src/kernel/mrm.rs` |
| Indexed mzML | FORMAT/HANDLERS/IndexedMzMLHandler | 16 ported | `src/format/indexed_mzml_handler.rs` |

311 upstream class-test sections were ported and 3 mapped with cited evidence,
0 unaccounted. All transcribed literals, which is tier 3 evidence under
[the differential validation policy](DIFFERENTIAL_VALIDATION.md). No C++ was
executed. Thirty-five further C++ defects were recorded as CPP-076 to CPP-110,
including an `MRMFeature` lookup by unknown key that default-inserts into the map
and returns the first feature, a `BaseFeature::sortPeptideIdentifications`
comparator that is not a strict weak ordering and mutates its own arguments, and
an `IndexedMzMLHandler::openFile` that accumulates index state instead of
replacing it.

**No audit returned a blocker.** Every auditor independently re-ran
`cargo +1.85.0 check --locked --all-features --all-targets`, the full nextest
sweep, clippy, rustdoc and the doc-coverage report, enumerated each owned
header's public members by hand against the support doc's API table, and counted
`START_SECTION`s itself. Four verdicts were `accept_with_fixes` and one `accept`.

**Both major findings were claim defects, not code defects, and both were
fixed before merge.** The ion-mobility package proposed a ledger entry for
`MSSpectrum.h` that would have *replaced* the accumulated record of three earlier
waves — `documentation` and `scope` pointed only at the new support document.
Applied verbatim it would have erased the InstrumentSettings, AcquisitionInfo,
`getType(bool)` precedence and record-metadata review history. The integrator
merged instead: `rust` and `tests` unioned, the prior scope preserved and the new
wave's sentence appended naming its document. The indexed-mzML package documented
`setSkipXMLChecks` as "not ported" in three table rows and its provenance
manifest; the auditor traced `options_mut().skip_xml_checks` through
`mzml::read_with_load_options` to the Base64 whitespace strip at
`src/format/mzml.rs:2386`, which is the entire effect of the flag in the source
as well. All four claims were corrected to record the member as ported.

Three minor findings against the ion-mobility package are recorded rather than
fixed: an unreachable `DriftTimeUnit::None` branch documented as reachable, a
non-finite rejection whose stated rationale covers NaN but not infinities, and
the bottom-up chunk merge whose stability argument is untested because both
presorted fixtures use all-distinct m/z. The last is a real coverage gap — no
case has equal m/z spanning a chunk boundary.

Full sweep on the 384-core node at the merge commit: build 32 s, nextest
all-features 2,763 passed, doctests 8 passed, nextest no-default-features 2,247
passed, clippy `-D warnings` clean, `cargo +1.85.0 check --all-features
--all-targets` clean, rustdoc `-D warnings` clean, `cargo fmt --check` clean, and
all six Python gates green after the integrator regenerated the three generated
files.

## Kernel wave 1: ranges, predicates, helpers, gap-0 review and geometry (2026-09-12)

Five work packages ported in parallel git worktrees, then rebased onto wave 0,
merged and verified together. Kernel headers closed or native-equivalent rise
from 8 of 34 to 18 of 34.

| Work package | Headers | Class-test sections | Rust |
| --- | --- | --- | --- |
| Ranges | RangeManager, SpectrumRangeManager, ChromatogramRangeManager | 48 ported | `src/kernel/ranges.rs` |
| Predicates | RangeUtils | 33 ported | `src/kernel/range_utils.rs` |
| Helpers | SpectrumHelper | 9 ported | `src/kernel/spectrum_helper.rs` |
| Gap-0 review | DPeak, StandardTypes, RichPeak2D, FeatureHandle, BinnedSpectrum | 51 ported | `src/kernel/gap_closures.rs` and existing modules |
| Geometry | DPosition, DIntervalBase, DRange | 104 ported | `src/data_structures/{dposition,dinterval,drange}.rs` |

245 upstream class-test sections were ported with transcribed literals, which is
tier 3 evidence under [the differential validation policy](DIFFERENTIAL_VALIDATION.md).
No C++ was executed. Fifteen further C++ defects were recorded as CPP-061 to
CPP-075, including a `makePeakPositionUnique` swap that discards the whole
spectrum record while warning only about data arrays, and a `DRange::united` of
two empty ranges that returns the universal range.

**Ranges are computed on demand, deliberately.** The source caches ranges in a
mutable member refreshed by `updateRanges()`. Peak vectors are public here, so a
cache cannot be invalidated soundly; the algebra is ported as pure values and the
inherited container surface becomes `range_manager()` accessors. All 21 TOPP
`updateRanges` call sites were audited before choosing this: every one is a plain
update-then-read, and the mutable `getRange()` has no mutating caller in core, so
no ported tool's numbers change. The combined experiment role now folds in
chromatogram retention time, intensity and **product m/z**, which the previous
`MSExperiment::ranges()` omitted entirely; a ported `FileInfo` would have printed
a narrower m/z range.

**Two defects came from stale worktree bases.** Three of the five agents branched
from `a463e3e`, 26 commits behind, so their own green gates were green against a
tree without the typed record metadata, the drift-time fields or the 1.85 fixes.
Rebasing exposed both: two SpectrumHelper test assertions used the pre-migration
metadata API, and `copySpectrumMeta` documented a drift-time deferral that wave 0
had already made obsolete. The implementation was correct by construction; the
rustdoc, support document and API table were not, and nothing tested it. Each
branch was rebased onto wave 0 and re-verified before merging.

**The no-default-features CI line was broken by wave 0 and is now fixed.** Gating
`pub mod cli` on `paramxml` left the five `tests/topp_*.rs` files using
`openms::cli` without a gate of their own, so `cargo test --locked
--no-default-features` (rust.yml line 22) failed to compile. Verified after wave 0
were `--all-features` and the Python gates, not that line. Each file now carries
`#![cfg(all(feature = "mzml", feature = "paramxml"))]`.

Recorded checks on the integrated tree (16-core Apple Silicon, cargo 1.96):

| Check | Result |
| --- | --- |
| `cargo build --locked --all-features --all-targets` | clean |
| `cargo test --locked --no-default-features` | clean, after the gate fix |
| `cargo fmt --all -- --check` | clean |
| `python3 tools/check_core_sdk.py` | 1,060 added source references agree |
| `python3 tools/core_sdk_coverage.py --write` | complete 16, native-equivalent 70, partial 20, evidence-requires-review 159, unmapped 521 |
| `python3 tools/check_doc_coverage.py --write` | 1,764 of 3,401 public items = 51.9%, all new modules at 100% |

**A remote build host caught a defect the local integration missed.** The full
sweep was also run on an IBMI HPC node (`kim`, 384 cores). It reported two
unresolved intra-doc links in `spectrum_helper` that the local run had not been
repeated after merging. The cause was the wave-1 merge resolution itself: adding
an outer `///` doc comment on `pub mod spectrum_helper` in `src/kernel.rs` makes
rustdoc resolve that module's inner `//!` links in the parent `kernel` module,
so `[`PeakContainer`]` reported "no item named `PeakContainer` in module
`kernel`" although the trait exists in the module. One work package had warned
of this mechanism in its own notes. Both links now use full crate paths, and the
module records why.

### Remote build host

`kim` was provisioned as a second verification host: both toolchains in shared
CephFS home (install once, visible on every node), source tree and `target/` on
node-local NVMe `/scratch`, and libxml2 2.15.4 with pkg-config in a shared
micromamba environment under `/ceph/ibmi/abi/oliver/envs/rustbuild`. Three
environment facts had to be discovered and are recorded in
`/scratch/kohlbach/openms-rs-env.sh`: `pkg-config` must be on `PATH` for the
`libxml` build script, the node ships `libclang.so.1` without its resource
headers so bindgen needs GCC 13's include directory, and the default
`ulimit -n` of 1024 starves parallel `rustc`.

| Check | kim (384 c, 96 jobs) | Mac (16 c) |
| --- | --- | --- |
| `build --all-features --all-targets` | 22 s | 16 s |
| `test --all-features --all-targets` | **48 s**, 2,506 passed | **176 s**, 2,235 passed |
| `test --no-default-features` | 31 s, 2,028 passed | — |
| `clippy --all-features --all-targets` | 20 s | — |
| `cargo +1.85.0 check --all-features --all-targets` | 16 s | clean |
| six Python gates | all pass | all pass |

The test suite runs 3.7 times faster; the cold build does not, because the
toolchain is read over CephFS. `target/` reaches 23 GB there against 4 GB
locally, which node-local `/scratch` absorbs.


## Kernel wave 0: MSRV, build baseline and scaffold (2026-09-12)

Preparation for the parallel kernel port. Three findings are recorded because
each corrects a claim made earlier in this project.

**The minimum Rust version was broken.** `cargo +1.85.0 check --locked
--all-features --all-targets` failed on release `5688775` with `E0658` at five
let-chain sites: `src/cli.rs:128,345`, `src/metadata/experimental_design.rs:163,721`
and `src/format/experimental_design_file.rs:178`. Let-chains stabilised in Rust
1.88; edition 2024 accepts the syntax, so rustc 1.96 never reported it and the
crate's own `rust-version = "1.85"` was not enforced locally. CI `minimum-rust`
was red. All five are rewritten as nested `if`s; the 35 tests covering those
sites pass unchanged. An adversarial review found three of the five; the 1.85
compiler found the other two, and is now the gate — not a grep.

**Build baseline, the first recorded.** Before any change `target/` held 27 GB
(53.7 GiB of files) for 225 test binaries with full DWARF. With
`[profile.dev] debug = "line-tables-only"` (test semantics unchanged):

| Check | Result |
| --- | --- |
| `cargo test --locked --all-features --all-targets --no-run --timings` | 16 s, clean tree |
| `cargo test --locked --all-features --all-targets` | 2,231 passed, 176 s |
| `cargo test --locked --all-features --doc` | 4 passed (excluded by `--all-targets`; 2,235 total) |
| `target/` after the full build | 4.0 GB |
| `cargo +1.85.0 check --locked --all-features --all-targets` | clean |
| `cargo +1.85.0 check --locked --no-default-features` | clean |
| `cargo fmt --all -- --check` | clean |
| `check_core_sdk`, `core_sdk_coverage`, `test_core_sdk_coverage`, `check_schema_feature_graph`, `check_doc_coverage` | pass |

Machine: 16-core Apple Silicon, cargo 1.96. `cargo-timing.html` is kept outside
the repository.

**The struct-literal risk was phantom.** The plan feared ~400 struct literals
would break when fields were added. A regex over `MSSpectrum {` counted return
types, `impl` blocks and closure bodies; the compiler reports **zero** `E0063`
missing-field errors after adding `MSSpectrum::{drift_time, drift_time_unit}`,
`MSExperiment::sql_run_id`, `BaseFeature::{primary_id, id_matches}` and
`ConsensusFeature::ratios`. Every real literal already used
`..Default::default()`, the crate's existing convention. An automated fixer
built on the same regex was reverted in full rather than patched.

Scaffold additions: the fields above, `Ratio` (ports `ConsensusFeature::Ratio`:
`ratio_value`, `denominator_ref`, `numerator_ref`), and `Error::{InvalidRange,
MissingInformation}` mapped to `ILLEGAL_PARAMETERS` and `MISSING_PARAMETERS`.
`src/kernel.rs`, `src/kernel/features.rs` and `src/error.rs` remain at 100%
rustdoc coverage. Identification data is attached to maps **by reference**, a
deliberate divergence from `FeatureMap.h:294`, because the graph is not `Clone`
and embedding it would strip `Clone`/`PartialEq` from both map types.

## Indexed mzML writing and binary normalization (2026-09-11)

[Recorded checks](mzml-output-validation.json) cover complete represented
[source writer options](MZML_WRITE_OPTIONS_SUPPORT.md) and the
[binary whitespace option](MZML_NORMALIZATION_SUPPORT.md). The group adds 24 tests;
MzMLFile remains partial while typed/noise/detector transport and XSD validation
are implemented separately. Spectrum mobility/IMPeakType and additional precursor
activation metadata/unit routes also remain.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,139 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/minimal mzML-validation selection | 373 passed |
| Strict Clippy on both compilers; release library, Rustdoc, Rustfmt | Passed |
| Source audit, extraction integrity and two independent fixture regenerations | Passed |
| Completion checks and optional checksum dependency boundary | Passed |

All 15 integrated checks passed on their first run. Selected counts overlap.
Twenty-six final extraction files remain exact; Cargo.toml retains only the
previous semantic-feature registration in addition to the frozen writer change.
The final reader's normalization patch starts from the exact reviewed writer
and earlier centroid/semantic registrations. All 28 source hashes and nine
fixture/generator/reused-resource hashes were independently verified.

The writer reuses one prepared header/binary payload across two precharged markup
passes. Independent Python checks verify the actual SHA-1 prefix and every
record/index ID and byte offset, including partial writes and UTF-8. The unchanged
indexed XSD passes actual validation; this is separate from index integrity.
Original Numpress literals and independent precision/whitespace input fixtures
are retained. CPP-049/050 are corrected, and both normalization settings retain
checked malformed-Base64 rejection for CPP-055. No C++ method execution or full
C++ SDK build is claimed.

The [published ProForma checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34588912206)
passed CI. The ledger remains at 72 complete/native-equivalent headers and 714 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,945 distinct current files, 220 historical references, 21 graph references
and 952 added references. The source issue log contains 56 entries, including the newly source-reviewed
CPP-056 synthetic-scan mobility omission. Full SDK
completion remains outstanding.

## mzML inspection, isolation loading and semantic validation (2026-09-11)

[Recorded checks](mzml-operations-validation.json) cover
[spectrum type/centroid inspection](MZML_CENTROID_SUPPORT.md),
[isolation-target loading](MZML_ISOLATION_SUPPORT.md) and the complete native
[MzMLValidator specialization](MZML_VALIDATOR_SUPPORT.md). This checkpoint adds
47 tests and closes one public header; MzMLFile itself remains partial.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,115 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/minimal mzML-validation selection | 349 passed |
| Spectrum type queries without format features, Rust 1.85 | Seven passed |
| Strict Clippy on both compilers; release library, Rustdoc, Rustfmt | Passed |
| Source audit, source projection, extraction integrity and completion checks | Passed |

All 15 integrated checks passed on their first run. Selected test counts overlap.
Of 34 integrated extraction files, 33 remain exact; the shared reader merges
only the frozen isolation implementation and centroid/semantic registrations.
All 47 source hashes and 16 fixture/resource/tool/projection hashes were checked,
including the unchanged DTA files, original validator inputs, mapping and reused
ontologies. Native tests verify cumulative limits, source event order, caller
option preservation and fresh validator state across documents. No C++ method
execution or full SDK build is claimed.

The [published XLMS checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34586771858)
passed CI. The issue log now contains 55 entries. New findings concern mixed
chromatogram primary roles, ineffective indexed-schema ID references and unchecked
source Base64 alphabet bytes. The [index-schema probe](mzml-index-schema-probe.json)
records actual acceptance of a dangling reference by the unchanged XSD; it does
not claim byte-offset or checksum validity.

The ledger records 72 complete/native-equivalent headers and 714 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,942 distinct current files, 220 historical references, 21 graph references
and 924 added references. Indexed writing, complete writer options, typed/noise
transport, runtime XSD validation, whitespace-normalization options and broader
SDK work remain outstanding.

## ProForma spectrum generation (2026-09-11)

[Recorded checks](proforma-spectra-validation.json) cover all six
[ProForma spectrum operations](PROFORMA_SPECTRA_SUPPORT.md), completing the pinned
ProForma public header's native operation groups. Source compatibility and
checked boundaries remain explicit; this is not full notation-standard certification.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,068 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/no-default library and adjacent chemistry suites | 274 passed |
| Strict Clippy on both compilers; release library, Rustdoc, Rustfmt | Passed |
| Source audit, ProForma fixture regeneration and completion checks | Passed |

All 12 integrated checks passed on their first run. Selected counts overlap;
this group adds 21 tests. Twelve frozen extraction files remain exact; the
thirteenth has only a root module-documentation correction. All 30 source hashes
and both fixture/generator hashes were independently verified. The fixture
retains ten full source sections and 18 literal assertions, including six
generation calls; it does not invent numerical spectrum oracles. Separate tests
exercise real-backend composition, repeated resolution/warning order, shared
resource limits and atomic registry publication. CPP-038/047 finite source
behavior is retained beside independent chemical/position expectations.
No C++ spectrum execution or full SDK build is claimed.

The [published general-validator checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34585986090)
passed all CI jobs. The source issue log contains 52 entries, including a
cross-document mzML parameter-group state leak and four-line XSD schema-selection
defect. The latter has independent pinned-schema validation evidence in
[mzml-schema-selection-probe.json](mzml-schema-selection-probe.json). The ledger now
records 71 complete/native-equivalent headers and 715 requiring implementation
or review, with zero certified TOPP workflows. Source verification covers 1,933
distinct current files, 220 historical references, 21 graph references and 877
added references. Remaining SDK work is tracked in the completion ledger.

## Streaming mzML consumers (2026-09-11)

[Recorded checks](consumer-validation.json) cover the complete unconditional
MSDataConsumer interface, both source [transform operation groups](MZML_CONSUMER_SUPPORT.md)
and disabled scientific data population with `fill_data=false`.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,047 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/mzML-only library and adjacent selection | 236 passed |
| Rust 1.85/unconditional interface without format features | One passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, consumer fixture regeneration and completion checks | Passed |

All 14 integrated checks passed on their first run. Selected counts overlap;
the full suite adds 24 tests. All 15 frozen extraction files remain exact.
Independent regeneration reproduces both fixture files after five documented
source child-count repairs, retaining original scientific data. Both transform
modes reproduce the source test's four spectra, 40 peaks and TIC 350. Source
control-flow review and independent tests cover setup order, separate pools,
mutation/stop/error boundaries, descriptor validation and atomic destination
publication. No C++ consumer execution or full SDK build is claimed.

The [published header checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34584903921)
passed cross-platform, minimum-Rust and quality CI. The C++ issue log contains 50
entries, including the indexed writer's placeholder checksum and empty dummy index. The ledger records
70 complete/native-equivalent headers, 716 requiring implementation or review,
and zero certified TOPP workflows. Source verification covers 1,933 distinct
current files, 220 historical references, 21 graph references and 847 added
references. Remaining mzML source options/validation, centroid inspection,
indexed output, ProForma wrappers and broader SDK work are tracked separately.

## Crosslink spectrum generation (2026-09-11)

[Recorded checks](xlms-validation.json) cover the complete class-specific
[XLMS generator](THEORETICAL_XLMS_SUPPORT.md), all 25 source options, and the
ProteinProteinCrossLink record/reaction enum. Other OPXL records remain separate.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,023 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/no-default library and adjacent selection | 176 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, XLMS fixture extraction and completion checks | Passed |

All 12 integrated checks passed on their first run. Selected counts overlap;
the full suite adds 25 tests. All ten frozen extraction files remain exact.
The independent extractor reproduces 52 literal source masses and 113 allowed
annotation strings. Native tests additionally cover source branches, cumulative
limits, sequence identity and atomic aligned appending. Review corrected linked
loss endpoint guards, empty-alpha pair behavior and portable charge-span handling
before integration. No C++ spectrum execution or full SDK build is claimed.

The generator explicitly retains the finite source suffix-loss and precursor
isotope defects (CPP-042/043), with separately calculated chemical expectations.
The C++ issue log contains 48 entries. The ledger records 69 complete or
native-equivalent headers, 717 requiring implementation or review and zero
certified TOPP workflows. Source verification covers 1,932 distinct current
files, 220 historical references, 21 graph references and 835 added references.
ProForma wrappers, mzML consumers and broader XLMS analysis remain ongoing work.

## General semantic validation (2026-09-11)

[Recorded checks](semantic-validation.json) cover the complete class-specific
[SemanticValidator group](SEMANTIC_VALIDATOR_SUPPORT.md) and the shared XML reader
used by CV mappings. Enable the optional `semantic-validation` feature.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,998 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/semantic feature and adjacent selection | 163 passed |
| Rust 1.85/CV mapping feature and library selection | 127 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, fixture projections and completion checks | Passed |

All 16 integrated checks passed on their first run. Selected counts overlap;
the full suite adds 19 tests. All 14 reviewed extraction files remain
byte-identical after integration. The historical valid fixture produces no
diagnostics; the corrupt fixture reproduces all five error and four warning
messages in source order using its original 738-term vocabulary. Four raw
fixtures retain their source bytes. No C++ execution or full SDK build is claimed.

The native validator corrects descendant-unit lookup, failed-parse state leakage
and history-dependent missing-path lookup (CPP-039/040/044). It retains the source
date conversion behavior (CPP-046); complete XSD conformance and derived format
validators are separate work. The C++ issue log contains 47 entries, including
the newly recorded ProForma flattened-range position defect.

The ledger records 68 complete/native-equivalent headers and 718 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,930 distinct current files, 220 historical references, 21 graph
references and 824 added references. Streaming consumers and the separately
reviewed XLMS/ProForma spectrum work are outside this checkpoint.

## mzML headers and peptide-evidence keys (2026-09-11)

[Recorded checks](header-evidence-validation.json) cover the complete
[source-supported header/reference group](MZML_HEADER_SUPPORT.md), metadata-only
reading and [PeptideEvidence value/key operations](PEPTIDE_EVIDENCE_SUPPORT.md).
DataProcessing timestamps now use DateTime across mzML, FeatureXML and ConsensusXML.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,979 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/header and adjacent selection | 321 passed |
| Rust 1.85/native selection | 171 passed |
| Final Rust 1.85/path and precursor workflows | 14 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, reference regeneration and completion checks | Passed |

Selected totals overlap. The full suite adds 26 tests. Two old rejection tests
needed capability updates: metadata-only paths and scalar precursor metadata
now succeed. Their replacement assertions verify the new results and retain
unsupported-value/error checks. An initial positive path test used a historical
projection lacking a required processing reference; the final test uses a complete
writer-produced document. All failures and successful reruns are retained in the
record, without claiming first-run success.

All 35 reviewed header files remain byte-identical after integration. Six
fixtures and the formatted mapping table regenerate from the pinned C++ source;
116 reader/writer instrument pairs agree. Both rich-header writer outputs pass
actual independent XSD validation. No new C++ execution is claimed.

The [published parent checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34582838457)
passed all five cross-platform/minimum/quality jobs. The C++ issue log now has
46 entries. The completion ledger records 67 complete/native-equivalent headers,
719 requiring implementation or review, and zero certified TOPP workflows.
Source verification covers 1,926 distinct current files, 220 historical references,
21 graph references and 799 added references. Consumers, ProForma spectra and
separately staged semantic validation remain outside this checkpoint.

## ProForma conversion and CV mappings (2026-09-11)

[Recorded checks](conversion-mapping-validation.json) cover complete ProForma
AASequence conversion in both directions and the five-class CV mapping group,
adding 30 tests against SDK `82ce5b3`.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,953 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Rust 1.85/CV mapping and library selection | 126 passed |
| Rust 1.85/native and adjacent selection | 260 passed |
| Rust 1.85/JSON and mass selection | 142 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, independent projections and completion ledger | Passed |

The [previous published checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34580543352)
passed all Linux, macOS, Windows, minimum-Rust and quality jobs. This is distinct
from CI for the present increment.

All 21 integrated checks passed on their first run. Selected totals overlap.
Every one of 1,078 executed C++ formatter cases is exercised through native
reverse conversion, writing and parsing. Only the exact source formatter helper
was compiled; this is not an executed C++ conversion or spectrum oracle. Mapping
tests compare all 683 full records projected independently from six unchanged
source XML files. Complete operation and evidence review is recorded in the
[conversion](PROFORMA_CONVERSION_SUPPORT.md) and [mapping](CV_MAPPING_SUPPORT.md)
support documents.

The [C++ issue log](../OpenMS_CPP_ISSUES.md) contains 38 entries. New source-reviewed
findings document charge loss during formula combination and duplicate linker
mass in the still-unported ProForma XLMS spectrum wrapper. Proposed C++ fixes,
executed evidence and native compatibility behavior remain separate.

The ledger records 66 complete or native-equivalent headers and 720 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,925 distinct current files, 220 historical references, 21 graph references
and 776 added references. Full mzML headers/consumers, ProForma spectra and general
semantic validation remain ongoing groups; staged work is not certified here.

## ProForma mass and controlled vocabularies (2026-09-11)

[Recorded checks](vocabulary-mass-validation.json) cover the complete ProForma
mass operation group and ControlledVocabulary at SDK `82ce5b3`, adding 37 tests.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,923 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Rust 1.85/native unit and adjacent selection | 240 passed |
| Rust 1.85/ProForma JSON and mass selection | 138 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, independent projections and completion ledger | Passed |

All 17 integrated checks passed. Selected totals overlap. Root review covered
the complete new operation groups, including shared registry transactions,
source mass pass order, iterative graph behavior and cumulative allocation
limits. The vocabulary tests compare every field of all 9,254 final terms and
all 16,852 name aliases with an independent source-loop projection. Five raw
providers retain their original bytes and separate data-license notices.
No executed C++ mass or vocabulary differential result is claimed.

The [previous checkpoint's CI](https://github.com/okohlbacher/openms-rs/actions/runs/34578513534)
passed every Linux, macOS, Windows, minimum-Rust and quality job, confirming the
Windows resource checkout correction. That result is distinct from CI for
this later change.

The [C++ issue log](../OpenMS_CPP_ISSUES.md) now contains 36 stable entries,
including ten newly source-reviewed defects found during subsequent header,
conversion and mapping work. Executed evidence, source deductions, proposed
upstream fixes and native handling remain distinguished.

The completion ledger records 61 complete or native-equivalent headers and
725 still requiring implementation or review, with zero certified TOPP workflows.
Source verification covers 1,915 distinct current files, 220 historical references,
21 graph references and 734 added references. Full mzML headers/consumers,
ProForma sequence conversion/spectra and CV mapping/semantic validation remain
separate ongoing groups.

## Experiment settings, DateTime and ProForma resolution (2026-09-11)

[Recorded checks](settings-resolution-validation.json) cover the integrated
settings ownership migration, DateTime operations and ProForma resolver at
SDK `82ce5b3`, with 42 additional tests.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,886 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Rust 1.85/native unit and adjacent selection | 400 passed |
| Rust 1.85/JSON and RNA selection | 155 passed |
| Rust 1.85/mzML and settings selection | 103 passed |
| Final private-test repeats, current/minimum | 110/102 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, independent projections and completion ledger | Passed |

Selected totals overlap the full suite. The record retains an initial minimum
lint failure and its successful retry after correcting a test initializer.
The original-file metabolite example still produces 81 features; its output
passes the original schema and 1,458 scalar/typed-metadata comparisons with the
source expected file. Default hull omission remains covered by direct tests.

DateTime retains 301 executed unmodified-source probe rows: 266 native matches
and 35 deliberately corrected early-year calendar results. Independent Python
month stepping regenerates all 35 corrections exactly. Three separate UBSan
executions expose signed fractional-second overflow. The
[C++ issue log](../OpenMS_CPP_ISSUES.md) records these and source-reviewed defects,
with an independent sorted-weight reproduction of the IMS nonprogressing witness.
Proposed upstream fixes are distinct from the native compatibility policy.

The [previous checkpoint's CI](https://github.com/okohlbacher/openms-rs/actions/runs/34576479236)
failed on Windows because Git converted the original SVM resources to CRLF.
Explicit byte-preservation rules fix that checkout issue; all four resources
now remain identical under Windows-style conversion. Cross-platform CI for
this checkpoint is a separate check, not inferred from local success.

The ledger records 60 complete or native-equivalent headers and 726 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,910 distinct current files, 220 historical references, 21 graph references
and 702 added references. Full header transport, streaming consumers and the
remaining ProForma scientific operations are ongoing work.

## Metabolite feature finding and experiment values (2026-09-11)

[Recorded checks](metabo-values-validation.json) cover the integrated feature
finder, DocumentIdentifier, five sample/instrument value types and the mzML
ASCII-subset Latin-1 fix against SDK `82ce5b3`.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,844 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Native unit/adjacent tests, Rust 1.85/no defaults | 306 passed |
| Rust 1.85/JSON and RNA selection | 134 passed |
| Rust 1.85/mzML count and adjacent selection | 70 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed; warnings denied |
| Source audit, model projections, ledger and regression checks | Passed |

The batch adds 41 tests; selected totals overlap the full suite. The feature
finder reproduces the source 83/81/80 counts and all 81 expected scientific
records, including typed metadata and compressed hulls. Both fixed classifiers
match 488 separately executed LIBSVM reference cases; all 7,755 model constants
are checked by bits. Four independent peptide scores cover source f32 rounding
and the declared native f64 precision boundary.

The file example reads the unchanged original mzML and writes 81 features.
Its featureXML passes the original schema; all scalar and typed metadata fields
match the source expected output, including 1,053 numerical comparisons. This
workflow exposed an ASCII-only Latin-1 declaration compatibility gap, now fixed
with explicit non-ASCII rejection. It does not certify complete TOPP behavior.

The ledger now has 58 complete or native-equivalent headers and 728 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,908 distinct current files, 220 historical references, 21 graph references
and 660 added references. The [previous checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34574946683)
passed every GitHub job. Full experiment/header transport, generic SVM support
and ProForma scientific backends remain separate work.


## Experiment metadata, feature hypotheses, ProForma JSON and mzML counts (2026-09-11)

All four groups are integrated against SDK `82ce5b3`. [Recorded checks](annotation-counts-validation.json)
retain exact commands, outcomes, log hashes, source scope and review fixes.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,803 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Native unit/adjacent tests, Rust 1.85/no defaults | 217 passed |
| Rust 1.85/JSON and RNA selection | 122 passed |
| Rust 1.85/mzML count and adjacent selection | 55 passed |
| Strict Clippy | All targets on current/all-features and minimum/no-defaults; minimum JSON/mzML targets also passed |
| Release library, Rustdoc and Rustfmt | Passed; documentation warnings denied |
| Source audit, projection, ledger and regression checks | Passed |

The batch adds 57 tests and one lifetime doctest. Selected totals overlap the full
suite. Source metadata quirks, ordered borrowing, all reachable JSON schema branches,
and the five source count pairs are covered, with independent resource and boundary
regressions. Reviews corrected retained allocation accounting, JSON exponent handling
and sparse map allowances, as well as XML declaration/attribute validation and tiny
input-chunk handling. The original scientific fixtures remain unchanged.

Source verification covers 1,894 distinct current files, 220 historical references,
21 graph references and 601 added references. The ledger records 51 complete or
native-equivalent headers, with 735 requiring implementation or review and zero
certified TOPP workflows. FeatureFindingMetabo orchestration, ProForma scientific
backends and full experiment/header/consumer support remain separate work. No new
C++ execution is claimed for this batch. The [elution and mapping checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34572998432)
passed all GitHub CI jobs.


## ProForma text parsing and structured errors (2026-09-11)

The complete source single-chain and ion grammars, all error codes and diagnostic
operations are integrated against SDK 82ce5b3.
[Recorded checks](proforma-parser-validation.json) retain commands, outcomes,
log hashes, staged tests and the C++ extraction boundary.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,746 tests passed |
| Doctests | Three passed |
| Native unit and adjacent tests, Rust 1.85/no defaults | 195 passed |
| Strict Clippy, all targets | Both compiler/feature configurations passed |
| Release library, Rustdoc and Rustfmt | Passed; documentation warnings denied |
| Source audit, projection, ledger and regression checks | Passed |

The 16 new tests include all 198 upstream grammar cases, 62 source component calls,
structured diagnostics and cumulative resource failures. A separately compiled
[exact-source probe](../tests/data/proforma_parser_probe_provenance.json) runs 238
inputs through both C++ grammars. All 476 comparisons match: 382 accepted cases
match both text modes, and 94 errors match code, byte position and original message.
The probe includes unchanged AST/tokenizer/parser/writer and prefix-helper blocks;
a capture-only exception adapter replaces SDK exception infrastructure. It does
not execute source exception formatting or a full C++ SDK. Frozen outputs add
no C++ dependency to native builds or CI.

Source verification covers 1,886 distinct current files, 220 historical references,
21 graph references and 556 added references. ProForma remains a partial SDK header:
JSON, resolution/conversion, mass/mz and spectrum methods are outstanding in this
checkpoint. The ledger remains at 48 complete/native-equivalent headers and 738
requiring implementation or review, with zero certified TOPP workflows.
The [earlier detection/acquisition checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34571964157)
passed every CI job, including Linux, macOS, Windows and minimum Rust.


## Elution-peak detection and identification run mapping (2026-09-11)

Both operation groups are integrated against SDK 82ce5b3.
[Recorded checks](elution-mapping-validation.json) include commands, outcomes,
log hashes, review scope and the resolved integration-formatting mismatch.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,730 tests passed |
| Doctests | Three passed |
| Native unit and adjacent tests, Rust 1.85/no defaults | 179 passed |
| Strict Clippy, all targets | Both compiler/feature configurations passed |
| Release library, Rustdoc and Rustfmt | Passed; documentation warnings denied |
| Source audit, projection, ledger and regression checks | Passed |

The 24 new tests cover all six elution options, source trace splitting and
smoothing, extrema and width/noise calculations, as well as complete run/path
mapping and merged-file selection. Elution also has 49 distinct passing staged
checks per compiler, overlapping the full suite. Source numerical fixtures and
independent small-case oracles remain unchanged. Both scientific results and
input trace updates roll back on operation failure; external progress output
has a separate documented scope.

Run mapping preserves the source's deliberate duplicate-error state: complete
forward mappings remain available while reverse mappings stop before the first
collision. Resource failures preserve the prior mapping. Source and native
implementations were independently reviewed. No C++ execution is claimed for
these two groups.

Source verification covers 1,886 distinct current files, 220 historical references,
21 graph references and 548 added references. The ledger now records 48 complete
or native-equivalent headers and 738 requiring implementation or review, with
zero certified TOPP workflows.


## Detection, ProForma writing, protein runs and mzML acquisitions (2026-09-11)

All four additions are integrated against SDK 82ce5b3. [Recorded checks](detection-acquisition-validation.json)
include exact commands, outcomes, log hashes, staged validation and review fixes.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,706 tests passed |
| Doctests | Three passed |
| Native unit and adjacent tests, Rust 1.85/no defaults | 155 passed |
| Staged mzML acquisition and adjacent tests | 101 current / 100 minimum passed |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Staged mzML-only strict Clippy, Rust 1.85 | Passed |
| Release library, Rustdoc, Rustfmt | Passed; documentation warnings denied |
| Source audit, projection, ledger and regression checks | Passed |

This group adds 48 tests; selected and staged totals overlap the full suite.
The [ProForma extraction probe](../tests/data/proforma_writer_probe_provenance.json)
actually compiled the exact C++ annotation declarations and complete writer,
with no scientific substitutions or backend linkage. All 160 output cases match
native output, covering both modes, twenty precise float bit patterns and four
formatting/chain scenarios. This is not a full C++ SDK build. The other three
increments use source fixtures, branch analysis and independent native oracles.

Mass-trace detection retains source growth/termination and metadata rules, with
explicit atomic failure and reusable-state corrections. Protein-run helpers
preserve target result ownership; native lexical metadata ordering is distinguished
from C++ registry ordering. The mzML tests include real independent XSD checks
for acquisition metadata combined with zoom and scan windows in both writers.
Zoom CV placement repairs a source ordering bug; all supported parameters precede
scan windows. Historical fixture bytes and pins are retained.

Source verification covers 1,883 distinct current files, 220 historical references,
21 graph references and 533 added references. The ledger records 46 complete or
native-equivalent headers and 740 requiring implementation or review; no TOPP
workflow is certified. The [preceding published checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34570080887)
passed all CI jobs, including Linux, macOS, Windows and minimum Rust.


## Mass traces, constants, monosaccharides and mzML settings (2026-09-11)

The four additions are integrated against SDK 82ce5b3.
[Recorded checks](trace-settings-validation.json) include commands, outcomes,
log hashes, staged scope and the executed C++ constants probe.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,658 tests passed |
| Doctests | Three passed |
| Native unit and adjacent tests, Rust 1.85/no defaults | 153 passed |
| Staged mzML settings and adjacent tests | 91 passed on Rust 1.98; 90 on Rust 1.85/mzML only |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Staged mzML-only strict Clippy, Rust 1.85 | Passed |
| Release library, Rustdoc, Rustfmt | Passed; documentation warnings denied |
| Source audit, projection and ledger/regression checks | Passed |

This group adds 29 tests. Selected and staged totals overlap the full suite.
All 38 numeric declarations and 91 metadata strings match an actually compiled,
unchanged C++ Constants header. Its unused configuration include needs only an
empty shim, with no scientific declarations replaced. This is not a full SDK
build. The other three ports do not claim C++ execution.

MassTrace tests preserve source numerical/cache quirks and use independently
computed expected values where broad upstream comparison constants are stale.
Monosaccharides retain every source field, literal mass and synonym precedence.
The mzML projection preserves scientific literals and explicitly repairs its
original Product count mismatch; ordinary and Numpress writer outputs pass
independent schema tests. Acquisition fields outside this expanded representation
remain guarded before output.

The [preceding acquisition checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34568959479)
passed every CI job. The ledger now records 45 complete/native-equivalent headers,
741 requiring implementation or review, and no certified TOPP workflow.


## Acquisition, chromatogram conversion and SDK refresh (2026-09-11)

The acquisition fields, [ChromatogramTools](CHROMATOGRAM_TOOLS_SUPPORT.md),
[processing propagation](PROCESSING_ACQUISITION_SUPPORT.md) and
[mzML write guards](MZML_ACQUISITION_GUARDS.md) are integrated against
[SDK 82ce5b3](CORE_SDK_82CE5B3_REVIEW.md). [Recorded checks](acquisition-validation.json)
include exact commands, outcomes, log hashes and separate staged scope.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,629 tests passed |
| Doctests | Three passed |
| Native acquisition and adjacent tests, Rust 1.85/no defaults | 167 passed |
| Staged mzML guards and adjacent tests | 79 passed on Rust 1.98; 78 on Rust 1.85/mzML only |
| Final signed-zero guard regression | Seven passed on each Rust version |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Staged mzML-only strict Clippy, Rust 1.85 | Passed |
| Release library, Rustdoc, Rustfmt | Passed; documentation warnings denied |
| Source audit and ledger/regression checks | Passed; 1,872 distinct current files verified |

This group adds 27 tests. The full count also includes the preceding seven-test
SequenceCoverage addition; selected totals overlap. Review covered exact grouping
and encounter order, atomic conversion failures, preservation of owned metadata
and shared processing identity, cumulative nested copy budgets, and early XML
rejection for fields that cannot yet be represented. Existing numerical
processing operations are unchanged.

Both preceding published checkpoints passed every CI job:
[5a90169](https://github.com/okohlbacher/openms-rs/actions/runs/34567597493) and
[4d2b648](https://github.com/okohlbacher/openms-rs/actions/runs/34567824647).
These are separate from the local checks for this new group. No new C++ execution
or full SDK build is claimed. The ledger records 42 complete/native-equivalent
headers, 744 requiring implementation or review, and no certified TOPP workflow.


## Standalone sequence coverage (2026-09-11)

[SequenceCoverage](SEQUENCE_COVERAGE_SUPPORT.md) is integrated. Seven new tests
cover the complete source operation, independent positional enumeration and
bounded failures. [Recorded checks](sequence-coverage-validation.json) include
82 unit/selected tests on Rust 1.98/all features and 78 on Rust 1.85/no defaults,
strict scoped Clippy on both, and the release library build. Formatting, source
and ledger checks also pass; 1,860 distinct current source files were verified.

This standalone addition leaves prior operations unchanged. Its focused totals
overlap existing tests and supplement the preceding full 1,595-test/three-doctest
checkpoint. No new C++ execution or full combined suite is claimed here. The
ledger now has 41 complete/native-equivalent headers; 745 require implementation
or review, and no TOPP workflow is certified.


## Unique IDs, 2D conversion, IMS solvers and mzML Numpress (2026-09-11)

The four operation groups are integrated. [Recorded checks](transport-values-validation.json)
include exact commands, outcomes, log hashes and separate staged validation scope.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,595 tests passed |
| Doctests | Three passed |
| Native/decoy tests, Rust 1.85/no defaults | 108 passed |
| mzML Numpress and adjacent staged tests, both Rust versions | 69 passed per configuration |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Scoped mzML-only strict Clippy, Rust 1.85 | Passed |
| Release library, Rustdoc, Rustfmt | Passed; documentation warnings denied |
| Source audit and ledger/regression checks | Passed; 1,859 distinct current files verified |

This group adds 48 tests. Selected totals overlap the full suite. Independent
schema tests validate both Numpress writer variants, and 36 unchanged upstream
binary payloads decode to 342 independently checked points. No new C++ execution
or full SDK build is claimed; the earlier raw-codec probes remain separate evidence.

Reviews covered ID word consumption and native-endian UUID layout, 2D grouping
and direct integer metadata conversion, source decomposition table/order quirks,
and mzML codec/type/fallback behavior. A source zero-witness loop now returns an
error; binary validation work is bounded before mzML scalar traversal. Ordinary
fallback retains its existing numeric precision and optional zlib setting.

The ledger records 40 complete/native-equivalent headers and 746 requiring
implementation or review. Full SDK parity and TOPP readiness remain open;
no TOPP workflow is certified.


## IMS foundations, peak traversal and Numpress (2026-09-11)

The isotope/element/alphabet, area/peak-export/index, raw Numpress and configurable
base64/zlib wrapper APIs are integrated. [Recorded checks](foundation-validation.json)
include commands, outcomes, log hashes and the scope of C++ reference execution.

| Check | Result |
| --- | --- |
| New/adjacent native suites, Rust 1.98/all features | 167 tests passed |
| Corresponding native suites, Rust 1.85/no defaults | 164 tests passed |
| Wrapper/raw/options/mzML, Rust 1.98/all features | 121 tests passed |
| Wrapper/raw/options, Rust 1.85/Numpress only | 104 tests passed |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Numpress-only strict Clippy, Rust 1.85 | Passed |
| Doctests, Rustdoc, release library, Rustfmt | Passed; three doctests and documentation warnings denied |
| Source audit and review/ledger regressions | Passed; 1,850 distinct current files verified |

This additive group introduces 81 integration tests. Selected totals overlap;
they supplement the preceding full 1,466-test suite rather than establish an
additive total. The [preceding commit](https://github.com/okohlbacher/openms-rs/actions/runs/34565942433)
passed every Linux, macOS, Windows, minimum-Rust and quality CI job. A missing
cached dependency interrupted the first documentation attempt; an isolated build
completed it and all remaining checks without production changes.

Raw Numpress tests compare 295 cases against an actually compiled, unmodified
pinned C++ implementation. Of those, 287 execute its decoder; eight empty Safe
cases avoid undefined source decoding. Encoded bytes and fixed-point helpers are
exact; SLOF decoding permits a documented host-math tolerance. Wrapper transport
fixtures are independently derived Python projections, not C++ wrapper runs.
The full C++ SDK has not been built or differentially validated.

Source reviews covered isotope convolution order, portable parser boundaries,
borrowed mutable traversal, source RT grouping, and compression rejection/fallback
semantics. The ledger records 33 complete/native-equivalent headers, with 753
requiring implementation or review. mzML Numpress wiring and full TOPP readiness
remain open; no TOPP workflow is certified.

## Mobility, array descriptions, weights and mzML paths (2026-09-11)

Mobilogram operations, IMSWeights, generic array descriptions and mzML file APIs
are integrated. [Recorded checks](mobility-validation.json) include commands,
results and log hashes.

| Check | Result |
| --- | --- |
| Rust 1.98, all features and all targets | 1,466 tests passed; examples compiled |
| Rust 1.85, selected native and adjacent processing suites | 137 tests passed |
| Rust 1.85, selected XML and path suites | 72 tests passed |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Selected XML Clippy, Rust 1.85 | Passed |
| Documentation examples, Rustdoc, release library, Rustfmt | Passed; three doctests and documentation warnings denied |
| Source audit and review/ledger regressions | Passed; 1,842 distinct current files verified |

Focused counts overlap the full suite. Independent source review covered mobility
sorting/search and annotation alignment, weight quantization/GCD quirks, description
preservation and XML rejection before publication, and compressed path I/O. Array
struct literals now need the description fields or `..Default::default()`; existing
constructors retain their use. Full array-description XML transport remains open.

The preceding [scientific-operations commit](https://github.com/okohlbacher/openms-rs/actions/runs/34564499936)
passed all CI jobs on Linux, macOS, Windows, minimum Rust and quality checks.
The ledger records 25 complete/native-equivalent headers and 761 still requiring
implementation or review. No TOPP workflow is certified.

## Scientific loading, decomposition and experiment operations (2026-09-11)

The native mass-decomposition solver, mzML scientific filtering and canonical
array types, experiment summaries, and idXML filesystem/dispatch operations are
integrated. [Recorded checks](scientific-operations-validation.json) include the
commands, outcomes and log hashes.

| Check | Result |
| --- | --- |
| Rust 1.98, all features and all targets | 1,432 tests passed; examples compiled |
| Rust 1.85, mzML-only plus solver/summary/value/unit suites | 147 tests passed |
| Rust 1.85, idXML-only paths/definitions plus summaries | 39 tests passed |
| Strict Clippy, all targets | Passed on Rust 1.98 with all features and Rust 1.85 without defaults |
| Focused mzML/solver/summary Clippy, Rust 1.85 | Passed |
| Documentation examples, Rustdoc, optimized library, Rustfmt | Passed; three doctests and documentation warnings denied |
| Source/fixture audit | 1,840 distinct current source/registration/reference files verified |
| Source-review and completion-ledger regressions | Three source-review tests and two ledger tests passed |

Selected test counts overlap; they are not added to the full-suite count.
Independent source review covered decomposition order and finite residue-table
limitations, mzML raw-precision filtering and aligned annotations, summary
floating-point accumulation, and idXML plain-output/extension behavior. Native
regressions cover bounded work, corrupt input and atomic failure. A preexisting
Rust 1.85 test-expression lint was corrected without changing fixture bytes.

The preceding published [extraction increment](https://github.com/okohlbacher/openms-rs/actions/runs/34561162608)
and [value-API increment](https://github.com/okohlbacher/openms-rs/actions/runs/34561620258)
passed every CI job: Linux, macOS, Windows, minimum Rust and quality checks.

The completion ledger records 23 complete/native-equivalent headers and 763
requiring implementation or review. The solver's scientific API is represented;
standalone IMS utilities are separately tracked. Full mzML metadata, codecs,
consumer/transform behavior and many other SDK APIs remain. No C++ runtime
comparison or certified TOPP workflow is claimed.


## Peak-file, metadata and composition values (2026-09-11)

PeakFileOptions, equality-compatible metadata/Product hashing, and the complete
MassDecomposition count-container API are integrated. [Recorded checks](value-apis-validation.json)
include commands, outcomes and log hashes.

| Check | Result |
| --- | --- |
| Combined new/adjacent suites, Rust 1.98, all features | 58 tests passed |
| Same applicable suites, Rust 1.85 without default features | 49 tests passed |
| Strict Clippy, all targets | Passed on Rust 1.98 with all features and Rust 1.85 without defaults |
| Doctests, optimized library, Rustfmt | Passed; three doctests |
| Source audit and completion ledger | Passed |

These selected runs supplement the preceding full 1,356-test suite; their counts
overlap and are not an additive full-suite total. Independent reviews covered
source option activation/defaults and the distinct `+` versus `+=` cached-maximum
semantics in MassDecomposition. Native hash tests record every equality-significant
field and signed-zero normalization without assuming cross-language digest values.

At this checkpoint PeakFileOptions and MassDecomposition were value APIs. The
following scientific-operations increment adds loading execution and the solver. Product's own hash
operation is now present, while inherited CVTerm/DataValue independent-unit states
and numeric registry semantics remain under review. The inventory records 22
reviewed complete/native-equivalent headers and 764 still requiring work or review.
No C++ differential execution or certified TOPP workflow is claimed.

## SDK update and extraction operations (2026-09-11)

The port now targets Core SDK 4.0.0 at `54a232f`. This increment completes the
reviewed native FASTA lifecycle and indexed mzML offset decoder, adds experiment
aggregation and XIC extraction with mzML Product interchange, and implements
peak display/hash traits. The original provenance is retained for carried-forward
fixtures; changed source files receive explicit review records.
[Recorded results](extraction-validation.json) contain the commands and log hashes.

| Check | Result |
| --- | --- |
| Rust 1.98, all features and all targets | 1,356 tests passed; examples compiled |
| Documentation examples | Three passed |
| Strict Clippy, all targets | Passed on Rust 1.98 with all features and Rust 1.85 without defaults |
| Rustdoc, optimized library, Rustfmt | Passed; documentation warnings denied |
| Indexed mzML alone, Rust 1.85 | Nine tests passed |
| Source/fixture audit | 1,834 distinct source/registration/reference files verified |
| Source-review and coverage regressions | Three source-review tests and two ledger tests passed |

All added components also passed their targeted current/minimum compiler tests
before integration. These overlapping runs are not added to the full-suite count.
Independent review closed XML lexical validation, duplicate-attribute work limits,
FASTA byte-wrapping and mzML Product metadata issues. Existing EMG regression tests
cover the private fitted-trace storage adjustment required by the larger Product
representation.

The preceding published commit `8fb47cb` passed Linux, macOS, Windows and Rust 1.85
CI, including 1,310 tests in the Linux all-target suite. Its quality job found a
new Rust 1.98 test-expression lint; the corrected expression and entire combined
crate now pass strict local Rust 1.98 checks.

The completion inventory still records 766 headers requiring implementation or
review. Twenty headers have a reviewed complete implementation or native equivalent;
this does not establish full dependency or TOPP workflow parity. No C++ executable
was built or run and no TOPP workflow is certified.

## Logging, progress and remaining map helpers (2026-09-10)

Owned logging, replaceable progress reporting, public feature/consensus
modification collection and compressed INI loading are integrated. Real local
timestamps and process CPU timing use pinned safe Rust adapters. The optimized
all-feature library builds. [Recorded results](runtime-foundations-validation.json)
include each command, outcome and log hash.

| Check | Result |
| --- | --- |
| Combined runtime, XML, definition, filesystem and unit suites, current Rust | 206 tests passed |
| Selected runtime, definition, filesystem and unit suites, Rust 1.85 without defaults | 128 tests passed |
| INI feature alone plus filesystem integration, Rust 1.85 | 32 tests passed |
| Final unit and CSV regression checks | 73 tests passed |
| Strict Clippy, all targets | Passed on current Rust with all features and Rust 1.85 without defaults |
| Doctests, Rustdoc, optimized library, Rustfmt | Passed; three doctests and documentation warnings denied |
| Source/fixture audit | 1,827 distinct current source/registration/reference files; four unchanged logging fixtures |
| Windows-style Git checkout filters | All four tested source fixtures retain exact bytes |
| Completion ledger and regressions | All 786 headers accounted for; two regression tests passed |

These selected runs overlap and supplement the complete 1,276-test suite in the
preceding increment. They are not a new complete all-target run. Reviews closed
an empty-name work-accounting gap and a logging flush-failure path that allowed
later output. Logging and progress also passed their isolated staged tests on
both compilers before integration.

CI for the preceding commit passed macOS and Rust 1.85. Windows exposed CSV
fixture checkout conversion, and the quality job exposed a Clippy test-expression
diagnostic. This increment protects all scientific fixture bytes and uses a
portable mutable byte array; local regression checks passed. Updated CI results
are tracked separately from these local results.

The completion gate still fails for 770 headers requiring implementation or
explicit review. Public modification collection and the native progress API now
have complete reviewed mappings. ZIP input, some platform logging/filesystem
behavior, full format/API coverage and executed C++/TOPP parity remain open.


## Map interchange and runtime resources (2026-09-10)

The target remains Core SDK 4.0.0 at `6bfc0e4`. Native featureXML and consensusXML,
portable modification definitions, typed feature metadata, filesystem helpers,
and shared gzip/bzip2 transport are integrated. The complete current-Rust suite
passes **1,276 tests**. [Machine-readable results](map-interchange-validation.json)
record the commands, outcomes and log hashes.

| Check | Result |
| --- | --- |
| Current Rust, all features and all targets | 1,276 tests passed; examples compiled |
| Rust 1.85, consensusXML alone plus unit/dispatch suites | 86 tests passed |
| Rust 1.85, featureXML alone plus unit/dispatch suites | 85 tests passed |
| Rust 1.85, idXML alone plus unit/custom-definition suites | 100 tests passed |
| Rust 1.85, no default features, selected unit/filesystem/definition suites | 95 tests passed |
| Strict Clippy, all targets | Passed on current Rust with all features and Rust 1.85 without default features |
| Final strict idXML-only Clippy, Rust 1.85 | Passed after narrowing a map-only helper's feature gate |
| Doctests, Rustdoc, optimized library, Rustfmt | Passed; three doctests, documentation warnings denied |
| Source/fixture audit | 1,825 distinct current source/registration/reference files; 11 batch fixtures verified |
| Coverage inventory and its regressions | All 786 headers accounted for; two regression tests passed |

Focused test counts overlap and must not be added to the full-suite count. Both
map dialects were tested against original fixtures; consensusXML output also
passed the original XSD through the locally available `xmllint`. Metadata,
protein references, portable chemistry, compressed paths and atomic failures
have independent regression coverage. Reviews caught and closed reference
expansion and validation-order resource-limit gaps.

The completion review now records ZIP input as missing from both map XML
loaders, and compressed input as missing from the INI loader. The INI entry was
therefore corrected from complete to partial. Map-based public modification
collection helpers, platform-specific filesystem behavior, inherited XML schema
validation, and broad SDK parity remain open. The completion gate still fails
for 772 headers requiring implementation or explicit review; it is not a count
of wholly absent Rust classes. No C++ executable was built or run, and no TOPP
workflow is yet certified as port-ready.


## Historical SDK foundations increment (2026-09-10)

The target remains Core SDK 4.0.0 at `6bfc0e4`. This increment adds identification
cleanup, reusable mzML parameter groups, MS2/DTA2D, file dispatch, text/list
utilities and the parameter/INI lifecycle used by TOPP tools. The final optimized
all-feature library builds successfully. Results are recorded in
[sdk-completion-validation.json](sdk-completion-validation.json).

| Check | Result |
| --- | --- |
| Full all-feature/all-target suite before the additive configuration layer | 1,128 tests passed |
| Final unit, parameter, INI, text and list suites, current Rust | 137 tests passed |
| Same selected suites, Rust 1.85 without default features | 122 tests passed |
| INI feature alone, Rust 1.85 | 14 tests passed |
| Final file dispatch and unit regression checks | 67 tests passed |
| Selected file, graph and unit suites, Rust 1.85 without default features | 99 tests passed |
| Strict Clippy, all targets | Passed with all features on current Rust and no default features on Rust 1.85 |
| Doctests, Rustdoc, optimized library, Rustfmt | Passed; three doctests, documentation warnings denied |
| Pinned source/fixture audit | 1,813 distinct source, registration and reference files verified |
| Completion inventory and its two regression tests | Passed; all 786 registered SDK headers accounted for |

These overlapping test counts are separate runs, not an additive total. The
full 1,128-test run preceded the additive parameter/text layer; focused tests
cover that layer and final whole-crate lint, documentation and release checks
cover the combined implementation. An initial inventory check found a stale
generated ledger during development; regeneration and the final check passed.

The completion gate deliberately still fails: 772 headers need implementation
or further review. Fourteen headers have a reviewed complete implementation or
native equivalent. An unmapped header can have Rust functionality that still
needs explicit review; this is not a count of wholly missing classes. The
[completion ledger](CORE_SDK_COMPLETION.md) is an inventory, not a completion
percentage. Source fixtures, schema validation and native library workflows do
not establish executed C++ differential parity or certify a complete TOPP tool.

## Historical reduced-SDK graph increment (2026-09-10)

The port now targets Core SDK 4.0.0 revision `6bfc0e4`. The optimized all-feature
library builds successfully. The complete all-feature/all-target suite passed
1,034 tests on Rust 1.96 and Rust 1.85 before this additive observation/match
layer; the no-default current run passed 905 tests before the layer was added.
The additive graph/formula checks now pass 54 private tests, 12 graph-operation
tests, 12 graph-record tests, 9 source-reference tests, 8 observation-match
tests, and 8 peptide-formula tests on the current compiler. The corresponding
focused suites pass on Rust 1.85 with no default features. Strict focused Clippy,
Rustfmt, source-target inventory checks and the 442-hash/438-link audit pass.

This increment adds observation and compound records, graph adducts, typed
molecules, observation matches, best-match queries, translated ownership and
typed peptide fragment formulas. At that stage, graph groups, cleanup/persistence
and the legacy converter were outside the port. Later increments add groups,
cleanup and the sequence/evidence conversion bridge; graph persistence and full
conversion remain outstanding.

The machine-readable results for this increment are in
[graph-validation-results.json](graph-validation-results.json).

Validated locally on **2026-09-10, macOS ARM64**. The complete all-feature suite
passes **993 tests on Rust 1.96.0: 947 integration tests across 110 suites and
46 unit tests**. The declared minimum Rust **1.85.0** passes all 129 selected
RNA and unit tests with all features. Both compilers pass 124 selected RNA and
unit tests without default features, plus three all-feature doctests each.
Full-library parity remains in progress.

The [machine-readable results](validation-results.json) retain the earlier
complete 938-test runs on both compilers and identify the current checks with
`rna_processing_` names. The new complete 993-test run covers every implemented
chemistry, processing, kernel and format suite on Rust 1.96. The minimum-compiler
and no-default-feature runs in this increment are focused checks, not complete
993-test runs.

## Historical RNA processing checks

| Check | Result |
| --- | --- |
| All features and all targets, Rust 1.96 | 993 tests passed; all fifteen examples compiled |
| All features, Rust 1.85, all RNA suites and unit tests | 129 tests passed |
| No default features, both compilers, all RNA suites and unit tests | 124 tests passed on each compiler |
| Doctests, all features | Three passed on each compiler |
| Clippy, all features/targets, Rust 1.96, warnings denied | Passed for the complete crate, tests and fifteen examples |
| Clippy, no default features/all targets, Rust 1.85, warnings denied | Passed |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| RNA processing example | Executed; RNase_T1 digestion, uridine variants and annotated negative b/y fragments |
| RNA registry regeneration | Original JSON/custom TSV reproduce all 378 embedded records and 375 distinct codes |
| RNA enzyme regeneration | Original XML reproduces all fourteen embedded enzymes and 84 field values |
| Package file list | 428 files; no package build or publication; RNA data release terms tracked separately |
| Current provenance and documentation | 421 source/fixture hashes and 413 local documentation links verified; RNA processing source lines, registry fields and ion bits checked alongside all retained scientific fixtures |

This increment adds fifty integration tests in five suites: twelve RNase,
fourteen modification-generation, twelve spectrum-generation, nine independent
source-reference and three workflow tests. Five new private tests verify the
shared formula work and allocation allowance, including the spectrum generator's
formula-only precursor path. The selected matrix includes all earlier RNA
record/provider, sequence, source-reference and workflow tests to check the
shared sequence helper refactor on both compilers.

Independent references cover fourteen enzyme records, 6,048 pattern/code cases,
twelve original digestion cases with 38 products, all four modification counts
(7, 6, 27 and 432), and all 126 fragment masses actually compared in the source.
The workflows connect modified cleavage, positional products and variable
modifications to annotated spectra and both mzML compression modes. Tests also
cover retained metadata, source charge/sulfur conventions, resource exhaustion
and atomic failures. No expected scientific value was changed to fit an
implementation result. See [RNase support](RNASE_SUPPORT.md),
[RNA modification generation](RNA_MODIFICATION_SUPPORT.md),
[RNA spectra](RNA_SPECTRUM_SUPPORT.md) and the
[independent review](RNA_PROCESSING_REFERENCE_REVIEW.md).

The checks below are retained historical evidence, with their original scope.

## Retained complete RNA foundation checks

| Check | Result |
| --- | --- |
| All features and all targets, Rust 1.96 and 1.85 | 938 tests passed on each compiler: 897 integration tests across 105 suites and 41 unit tests |
| No default features and all targets, Rust 1.96 | 863 tests passed |
| No default features, Rust 1.85, original RNA suites and unit tests | 70 tests passed |
| JSON-only feature, Rust 1.96, original RNA suites and unit tests | 73 tests passed |
| Doctests, all features | Three passed on each compiler |
| Clippy, all features/targets, warnings denied | Passed for the complete crate, tests and fourteen examples |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| RNA example | Executed; charged formula/mass/m/z, sulfur-aware suffix and coarse isotope probabilities |
| RNA registry regeneration | All 378 embedded records and 375 distinct codes reproduced |
| Package file list | 403 files at that stage; no package build or publication |
| Provenance and documentation at that stage | 359 source/fixture hashes and 382 local documentation links verified |

The foundation added thirteen record/provider tests, nine sequence tests,
eight independent source-reference tests and three workflows. These cover all
fifteen original formula assertions, eighteen mono/five average mass assertions,
ten positive slice examples and every registry entry. The JSON dependency
exposed two test-only empty-array type-inference ambiguities, resolved with
equivalent `is_empty()` assertions. A narrow precision-lint allowance preserves
the original carbon-13 mass literal in one workflow; the affected workflows
were rerun on both compilers. See [RNA support](RNA_SUPPORT.md) and the
[foundation reference review](RNA_REFERENCE_REVIEW.md).

## Earlier Tagger-specific checks

| Check | Result |
| --- | --- |
| All features, Rust 1.96 and 1.85 | 58 tests passed on each compiler: all 41 unit tests plus 17 Tagger integration tests |
| No default features, Rust 1.96 and 1.85 | 57 tests passed on each compiler; the mzML workflow excluded |
| Clippy, all features/targets, warnings denied | Passed for the final crate, tests and thirteen examples |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| Tag extraction example | Executed; five measured peaks produce the exact sorted tags EP, EPT, PE, PEP, PEPT and PT |
| Package file list | 380 files, including Tagger, its example, source fixtures and documentation; no package build or publication |
| Provenance and documentation at that stage | 293 source/fixture hashes and 353 local documentation links verified; six Tagger count rows and 120 membership rows verified at their original source lines alongside retained scientific fixtures |

The independent suite reproduces all six original Tagger counts and all 120
membership assertions using the source's theoretical-spectrum settings. The
source input sizes, 357 and 180 peaks, are checked first. Derived tests isolate
strict mass-window boundaries, ties, exact collisions, modified I/L behavior,
signed/unsorted coordinates and zero/inverted settings. Twelve private mass-table
tests cover the 19 base residues, two-stage resolution, deterministic provider
order, empty short IDs and source free-residue arithmetic.

Three private traversal/sort tests include a 9,999-residue path without recursion
and failures after partial traversal that preserve old string allocations. The
workflows connect measured tags to exact target/decoy substring matching,
modified digestion fragments and both mzML compression modes. See
[Tagger support](TAGGER_SUPPORT.md) and its [independent review](TAGGER_REFERENCE_REVIEW.md).

## Retained decoy and adduct checks before Tagger

| Check | Result |
| --- | --- |
| All features, Rust 1.96 and 1.85 | 58 tests passed on each compiler: all 26 unit tests plus 32 new integration tests |
| No default features, Rust 1.96 and 1.85 | 56 tests passed on each compiler; both XML workflows excluded |
| idXML without mzML, Rust 1.96 | 57 tests passed; decoy idXML workflow included |
| Clippy, all features/targets, warnings denied | Passed for the final crate, tests and twelve examples |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| Decoy FASTA example | Executed; two targets and two deterministic variants per target, six records total |
| Package file list | 369 files, including full Boost license notices and the reproducible RNG oracle; no package build or publication |
| Provenance and documentation at that stage | 277 source/fixture hashes and 336 local documentation links verified; new source/Boost hashes and RNG oracle checked alongside all retained scientific fixtures |

Decoy reference tests replay all thirteen unique original sequence cases,
including their shared random/cache history. Nine private RNG tests cover four
complete cycles for three seeds, a published standard value, rejection draws,
shuffle order and state handling. Three additional private tests cover identity
and transactional failures. The packaged integer oracle regenerates independent
checkpoints, checksums and permutations without running C++.

Adduct tests preserve the original scalar examples, then distinguish their
rounded tolerance from exact atomic/electron arithmetic. They check parser
edge cases, signed formulas, finite negative masses and native limits. The
workflows compare sodium-adduct isotope masses with complete ion composition,
preserve charge/metadata through mzML, and connect decoy FASTA to indexing,
synthetic FDR values and idXML. See [decoy support](DECOY_GENERATION_SUPPORT.md),
[adduct support](ADDUCT_SUPPORT.md) and their linked reference reviews.

## Retained full baseline before decoy and adduct additions

| Check | Result |
| --- | --- |
| All features and all targets, Rust 1.96 | 815 integration tests and 14 unit tests passed; eleven examples compiled |
| Documentation examples/API encapsulation | 3 doctests passed |
| No default features | 747 integration tests, 14 unit tests and 3 doctests passed; mzML/idXML excluded |
| All features and all targets, Rust 1.85 | Same 815 integration tests and 14 unit tests passed |
| Rust 1.85 documentation tests | 3 passed |
| idXML without mzML | 772 integration tests, 14 unit tests and 3 doctests passed; no base64/zlib feature dependency |
| Clippy, all features/targets, warnings denied | Passed |
| Rustfmt check | Passed |
| Rustdoc, warnings denied | Passed |
| Spectrum processing example | 121 input peaks → 14 retained peaks; normalized TIC 1.000000 |
| FASTA digestion example | Executed; positional peptide/mass/m/z TSV produced |
| Modified peptide analysis example | Executed; masses, five-bin isotope envelope and annotated b/y fragments produced; invalid input independently checked |
| Synthetic identification example | Executed; selected AC(Carbamidomethyl)DMK with seven fragment matches and 26.32% protein coverage; precursor/evidence/metadata checked |
| Chromatogram integration example | Executed; two synthetic peaks, original sample boundaries, raw intensity sums, time-weighted areas and baseline estimates produced |
| EMG fitting example | Executed; 36 cropped samples → 80 fitted samples, 700 iterations; sampled area 2,121.476414 → 2,501.460284 versus known complete synthetic area 2,506.628275 |
| Profile processing example | Executed; 401 profile samples → four iterative centroids → two retained peaks with aligned integration and width annotations |
| Modified peptide enumeration example | Executed; two digest products → four alkylated/oxidized variants plus one unchanged product; masses and doubly charged m/z values printed |
| Isotope streaming/enrichment example | Executed; first five natural glucose configurations and two enriched-carbon configurations above absolute threshold 0.1, with raw/log probabilities |
| Peptide-property example | Executed; modified and unmodified peptides produce finite charge, pI, GRAVY and gas basicity at both 500 K and 100 K |
| Spectrum annotation example | Executed; eleven source IFSQVGK measurements annotated with source b/y labels, exact widened intensity 12.100000262260437, longest y-series six and finite matching statistics |
| mzML input processing example | Executed using the independent mixed-precision fixture |
| Independent XML XSD checks | Populated/empty mzML, native idXML, auxiliary-array and populated-precursor mzML writer output validated with xmllint and their pinned schemas |
| Generated modification data | All 3,035 specificity rows regenerate exactly; original XML matches pinned archive byte for byte |
| Generated enzyme data | All 33 source records and 602,316 independent Python regex contexts regenerate exactly; bundled XML is byte-identical to the pinned source |
| Package file list | 351 files; spectrum annotation/ion naming, shared-work regressions, peptide-property tables/numerical fixtures, original data, schemas, isotope streaming/custom inputs and fine isotope/precursor purity integration, theoretical extensions, XLMOD/OBO and prior scientific fixtures, examples and licenses included; local reference snapshot and target directory excluded |
| Source/fixture provenance and documentation | 252 source/fixture hashes and 307 local documentation links verified; 41 spectrum-annotation rows and source/derived bits, 633 peptide-property rows and their binary64 values, isotope streaming/custom inputs, fine isotope/precursor purity and theoretical-spectrum sources, all 148 XLMOD mass literals/IEEE values, byte-identical ontology, prior scientific fixtures and libm archive/lock checksum checked; earlier chromatogram/schema checks retained |

A constructor-level source review during validation corrected `precursor_in_ms2`
to integer 0/1 and the spectrum ppm flag to the string projection `"0"`/`"1"`.
The 21 affected annotation/reference/workflow tests were rerun on both compilers
with all features; the optional-feature matrix, lint, documentation and examples
completed after the correction.

The unchanged modification/enzyme generators and data retain their successful
checks from earlier on the same date; they were not rerun for these API edits.

The package-list check does not publish or build a package. Cargo reported only
optional missing documentation/homepage/repository metadata for this local port.

## Test coverage

| Suite | Tests | What it establishes |
| --- | ---: | --- |
| RNase records, registry and digestion | 12 | Fourteen source enzymes, complete record identity, replacement/alias semantics, modified-code cleavage, end gains, product order/coordinates and atomic limits |
| RNA modification generation | 14 | Fixed/variable selection, terminal-site and maximum-one source quirks, alternative identity/order, exact preflight counts/bytes and atomic failures |
| RNA spectrum generation | 12 | Nine ion series, single/multiple charge branches, source mass/sulfur/precursor conventions, aligned annotations and shared resource limits |
| Independent RNA processing references | 9 | All fourteen enzymes, 6,048 predicates, 38 source digestion products, modification counts/identities and 126 compared source ion literals |
| RNA processing workflows | 3 | Modified cleavage, positional digests and variant mass shifts, annotated spectra and both mzML compression modes |
| RNA formula and spectrum budget internals | 5 | Shared work/allocation depletion, preserved mass bits and input state, empty-sequence behavior and cumulative formula-only precursors |
| RNA records and providers | 13 | Complete fields, identity, ordered duplicate registry, ambiguity resolution, JSON/TSV differences, diagnostics and checked limits |
| RNA sequences | 9 | Source parser/display, custom ownership, terminal and sulfur slicing, all fragment variants, natural-H/electron mass convention and atomic failures |
| Independent RNA references | 8 | Fifteen original formulas, 23 mono/average mass references, ten slices, 378-record registry projection and source branch/identity cases |
| RNA workflows | 3 | Independent charged formulas, custom carbon-13 registry lifetime, coarse/fine isotopes and both mzML compression modes |
| Sequence tags | 8 | Defaults/setter, exact append ownership/order, finite source option edges, ignored spectrum fields, registry lifetime and atomic resource failures |
| Independent Tagger references | 6 | Six exact source counts, 120 membership assertions, source peak counts/rounded traces, strict bounds, nearest ties, collisions, I/L, signed/unsorted and no-op behavior |
| Tagger mass-table internals | 12 | Nineteen free-residue formulas, fixed/variable order, exact keys, two-pass provider resolution, terminal/wildcard/anonymous records, formula/mass precedence and construction limits |
| Tagger traversal/sort internals | 3 | A 9,999-edge heap traversal, failure after emitted paths, string-byte sort work and allocation-preserving atomic errors, strict lookup/ties |
| Tagger workflows | 3 | Exact substring target/decoy matching, modified digest fragment gaps and both mzML compression modes |
| Decoy generation | 7 | All registered enzymes, reversal/shuffle source rules, composition where applicable, modified/empty inputs, cache/resource boundaries and atomic failures |
| Independent decoy references | 7 | All thirteen literal source cases, ordered RNG/cache history, reseeding, cross-context reuse, zero attempts, short products and unspecific order |
| Decoy RNG internals | 9 | Three seeds over four cycles, standard 10000th word, rejection mapping, shuffle order, clone/reseed and no-draw boundaries |
| Decoy identity/state internals | 3 | Forward/reverse maximum identity, strict ties and late cache/work/allocation rollback |
| Decoy workflows | 3 | FASTA, target/decoy evidence and synthetic FDR, composition through digestion and idXML transport |
| Molecular adducts | 8 | Complete parser/components/getters/equality, electron and n-mer conversions, mono/average shifts, compatibility and checked limits |
| Independent adduct references | 5 | Literal source cases, tighter atomic constants, isotope/charge/whitespace grammar, signed containment and source numerical boundaries |
| Adduct workflows | 2 | Complete sodium-ion composition, unchanged isotope probabilities, charge-dependent spacing and both mzML compression modes |
| Chemistry | 12 | Source-derived mass/formula and digestion values, residue/isotope tables, formula grammar/algebra limits, b/y mass conservation |
| Charge and isoelectric point | 10 | Four source pKa scales, terminal overrides/suppression, parent PTMs, U/O/ambiguous residues, saturated finite pH, exact endpoint/midpoint behavior and resource/convergence guards |
| pI internal work accounting | 1 | Endpoint and midpoint evaluations consume one precharged budget |
| Amino-acid indices and gas basicity | 7 | Source accessions/indicators, ordinary split arithmetic, independent low-temperature values, tied maxima, empty/high-temperature identity, formula-free annotations and input limits |
| Hydrophobicity profiles | 7 | Seven scales, literal GRAVY/window/moment values, sliding order, all annotation types, window/angle semantics and preallocated work limits |
| Independent property references | 7 | 503 source constants including 42 sentinels, isolated pKas, 30 source scalar assertions, 420 singleton/pair GB expressions, 100 Decimal-derived extreme-temperature cases and independent moments |
| Peptide-property workflows | 4 | Digestion retains only original terminal caps, mass-only chemistry needs no formula, source generator attachment distinctions and typed property metadata through idXML |
| Ion naming | 9 | Every source charge/ordinal case, first-line and caret priority, field/overflow fallbacks, about 18,000 round trips, Unicode text, output limits and allocation-free parsing |
| Spectrum annotation | 10 | Three source operations, array replacement, ppm last-match/duplicate branches, no-op flags, final sorting, safe small-list statistics, finite source special cases, scoped validation and atomic failures |
| Independent annotation references | 7 | Original measured bits and 17 literal assertions, independently derived fragment errors, exact binary top-N padding/sample variance/ratios, source label grammars and precursor absolute tolerance in ppm mode |
| Annotation workflows | 4 | Modified digestion against independent fragment formulas, distinct chemistry with identical sequence text, peak annotations/statistics through idXML and aligned arrays/acquisition through both mzML compression modes |
| Shared generation/alignment work | 3 | Actual alignment initialization/cells and shared residue/loss/fine/coarse allowances survive successive calls; standalone calls start fresh budgets |
| Annotation grammar internals | 1 | Complete source label regex semantics, including commas, losses, ordinal boundaries and terminal-series distinctions |
| Sequence chemistry and numeric tags | 14 | Source integer/decimal registry lookup, absolute/internal/H/OH masses, unresolved residue representation, formula/mass availability, owned tags and shared immutable spelling through independently owned slices, stable attachment, slicing, atomic setters and numeric parser limits |
| Sequence identification operations | 8 | B/Z/X non-mass filters, known-mass precursor filtering, exact modification IDs, owned protein observations, distinct resolver keys and atomic failures |
| Sequence chemistry workflow | 8 | Ambiguous digestion/indexing, numeric-tag idXML round trips, independent fragment/loss shifts, formula-dependent append rejection and extreme-mass cancellation regressions |
| Modifications | 11 | Full pinned registry load, specificity/name/mass lookup, neutral losses, isotope labels, terminal notation, atomic setters, modified digestion/fragments |
| Modified peptide generation | 12 | Fixed/variable source goldens, stable alternative order, terminal overwrite/duplicates, existing anonymous annotations, bounded combinations, atomic errors and custom formula-free mass/absolute/no-op rules |
| Independent modification generation | 9 | Seven fixed and nineteen variable source cases, full weighted-site order, exact typed terminal/fragment-mass distinctions, append limits, empty idXML rejection and independent custom mass branches |
| Modification definitions | 13 | Owned definition identity, fixed/variable set semantics, compatibility, count non-enforcement, absolute/delta matching, inference, conflicting full-ID chemistry and atomic limits |
| Independent definition reference | 5 | Seven literal compatibility cases, partition/merged precedence, negative-delta tolerance endpoints, stored/fallback absolute masses and named/anonymous all-hit inference |
| Modified peptide workflow | 5 | Digestion → variants → formulas and independent fragment shifts → inferred search definitions → idXML; mass-only fragments, isotope rejection and serialization preflight |
| Owned modification registry | 7 | Caller record validation, optional vocabulary IDs, shared handles, bounded atomic OBO appends, formula fields and complete Eq/Ord contracts including signed zero |
| Independent OBO registry reference | 10 | All 148 XLMOD records in exact source order with identities, specificities, synonyms and mass bits; PSI aliases and absent targets, stanza/EOF handling, literal empty versus zero formula and alias-work limits |
| Crosslink lookup | 6 | Separate 56-record registry, source DSS/BS3/EDC mass/site goldens, reactive-side union, terminal conversion, search eligibility and bounded caller-owned extensions |
| Anonymous modification definitions | 6 | Owned exact-spelling annotations, all-hit inference, source empty-short-ID compatibility, full-residue/H/OH mass anchors, unresolved deltas and atomic errors |
| Caller-owned chemistry workflow | 7 | Registry lifetime release and sharing, numeric/name resolution, absolute-formula X rescue, no-change/delta/terminal precedence, non-UniMod mass export, full-ID-only records and exact custom-registry idXML |
| Custom chemical identity | 6 | Same text with distinct formulas/vocabularies remains distinct in protein observations, owned peptide keys, sequence duplicate filtering, rank/spectrum conflicts; atomic late errors |
| Internal theoretical fragments | 7 | Source interval and ten-residue boundaries, ordered annotations, numeric residue and terminal chemistry, loss replacement, isotope independence, intensities, atomic append and custom declaration/storage limits |
| Compact mass-only spectra | 7 | Full-length six-series ladders, source terminal/float operation order, observed masses, independent suffix accumulation, sorted no-ops and combined/precision limits with late-error atomicity |
| Activation presets and immonium ions | 6 | Source mass goldens, all activation enum cases, inferred precursor metadata, unmodified residue eligibility and L-only branch, charge/intensity/isotope independence and aligned atomic append |
| Independent theoretical extension review | 9 | Analytical internal interval/mass/loss oracles, exact immonium constants, source rounded CID/partial ECD tables and independently rounded f32 compact-helper values |
| Scalar purity and SPS matching | 9 | Source isolation totals/residuals, inclusive doubled tolerances, charge/nearest conventions, f32 fragment windows, input immutability and shared work/unused-annotation limits |
| Fuzzy and interpolated purity | 11 | Neutron-spacing source successor lookup, f32 sums/division, strict and half-weight boundaries, ratios above one, zero-window/empty-parent behavior, RT extrapolation/fallbacks and nonfinite-field checks |
| Independent precursor purity review | 13 | Exact decoded source peaks and seven scalar/map goldens, SPS float32 bounds and counts, independent fuzzy arithmetic and parent reference/acquisition-order rules |
| Precursor metadata and interchange workflow | 10 | Shared native precursor ownership, referenced parents, scalar batches, original complete mzML fixture, all activation/mobility quantities, malformed XML/loss guards and populated acquisition XSD validation |
| Native fine isotope enumeration | 8 | Small exhaustive configurations, fixed labels/natural zero abundances, charge adducts, source rounding/coverage, deterministic thresholds and checked atom/output/frontier resource limits |
| Independent fine isotope references | 7 | All 44 source counts, 14 fructose mass/probability rows, 6 bromine configurations, subnormal/full-support tails and the 19,615-state insulin f32 coverage boundary |
| Independent fine theoretical spectra | 7 | Literal 10/5/50/12 source counts and exact mass/intensity tables, independent CHNO oracles, neutral H/charge division, terminal/loss rules and single/shared-budget atomic failures |
| Fine isotope workflow | 3 | Modified peptide formulas and charge adducts, aligned spectrum selection, unknown-composition errors and both compressed/uncompressed mzML round trips with precursor acquisition |
| Fine isotope internals | 2 | Deterministic equal-probability state identity and cumulative work allowance across calls |
| Fine isotope stream and custom inputs | 8 | Materializer equivalence, original-mode thresholds, charge-ignoring raw adapter, owned lifetimes, fused errors, checked conversion, zero-count validation and large-support prefixes |
| Independent stream references | 6 | Original 14-row/2548-state fructose and 10000-state insulin ordered cases, threshold counts, independently reconstructed natural tables, custom f64 precision and direct multinomial products |
| Fine isotope stream boundary review | 10 | Extreme finite weights, log underflow, late mass overflow, getter-consistent cutoff equality and adjacent values, cumulative custom dimensions, deferred expansion and 103041-state full streaming |
| Custom isotope spectrum workflow | 3 | Owned raw-coverage prefixes and continuation, distinct equal-mass configurations through materialization, aligned spectrum selection and both mzML compression modes |
| Digestion registry and specificity | 12 | All 33 enzyme predicates vs independent regex hashes, full/semi/nonspecific order and small-sequence oracle, validity/missed cleavages, modifications and resource limits |
| Typed metadata and acquisition | 17 | Scalar/list/unit invariants, explicit string bridges, CV merging, source settings defaults/unify, drift/isolation conventions, date/checksum validation and atomic errors |
| Identification records | 18 | Stable scores, typed evidence/metadata, charge ranges, grouping, interval coverage vs independent oracle, observed modifications, kernel attachments and file loss rejection |
| Score categories and switching | 9 | All 29 names/six categories, source lookup/backup conventions, restoration, heterogeneous records, reserved labels and atomic map/slice errors |
| Peptide-spectrum scoring | 8 | Source HyperScore/Morpheus goldens, overload-specific precision/charge/boundary rules, ion annotations/ordinals, error means and checked failures |
| Identification filters | 16 | Source cutoff/top-N/dense-rank/tie behavior, modifications, exact/sequence duplicates, signed charge/precursor errors, run-aware references/groups and map atomicity |
| FDR and q-values | 13 | OMSSA1534 and XTandem source thresholds, picked goldens, actual legacy/Basic formulas, score ties/directions, peptide representatives, posterior estimates/ROC, reserved labels and atomic limits |
| idXML | 16 | Exact pinned fixtures, native metadata/run/evidence round trips, independent XSD validation, malformed XML/reserved encodings, byte/element/list guards and preflight/flush errors |
| Identification processing workflow | 1 | Independent XML → score switching/top-hit selection → hand-computable target/decoy q-values → threshold/reference cleanup → coverage → exact native idXML round trip |
| Peptide-to-protein indexing | 17 | 121 matcher/17 decoy/15 enzyme oracle cases, ambiguity/mismatch/I-L rules, original-engine recovery, per-run evidence/protein reconstruction, metadata placement and preflight/atomic limits |
| Basic protein inference | 21 | Source merger scores/groups, best/product/mean, charge/modification representatives, greedy graph ties/negative scores, score restoration, vector/single/consensus conventions and atomic limits |
| Identification conflicts | 11 | Source rank aggregation and spectrum reports, best/matching selection, modified-sequence/charge keys, intensity winners, metadata/subordinate retention and atomic map errors |
| Identification origin partitions | 7 | Three origin formats, source ordering/path reduction/hit union, basename/key collisions, run-specific protein values, skipped records, group omission and checked errors |
| Combined protein workflow | 1 | FASTA → modified/shared peptide evidence → basic inference → distinct protein/PSM q-values → global coverage → native idXML round trip → file-origin partitions |
| Retention-time models | 12 | Source weighted/unweighted fits, interpolation/extrapolation goldens, LOWESS sine/cars/original tables, inverses, diagnostics/windows and resource/atomic errors |
| Container RT transformations | 7 | Source experiment/feature/consensus/peptide values, original metadata, nested hulls/subordinates/handles, optional spectrum IDs, ordering and late-error atomicity |
| Isotope distributions | 13 | Coarse glucose/heavy/bromine goldens, convolution, enrichment, peptide/RNA/DNA averagine, exact Poisson recurrence rounding and overflow fallback, conditional fragment probabilities, limits |
| Theoretical spectra | 12 | All ion-series mass tables, losses, precursors, isotope conventions, modified peptides, aligned atomic append and resource limits |
| Kernel | 19 | Source search/TIC values, ties/ranges, metadata alignment, mutation checks and malformed inputs |
| Feature/consensus geometry | 18 | Scan-envelope hulls, source containment quirks, consensus means/charges/decharge, map identities/ranges, subordinates and checked mutations |
| Processing | 14 | Original 121-peak fixture, normalization, source 0.05 threshold default and adjacent-f32 cutoff, top-N/rank behavior, atomic errors and resampling conservation/boundaries |
| Smoothing | 10 | Gaussian goldens/trapezoids/ppm behavior; Savitzky–Golay goldens, asymmetric edges and polynomial reproduction |
| Baseline morphology | 8 | All ten source operations, boundary windows, signed bottom-hat, width conversion and aligned atomic updates |
| Peak picking | 12 | Orbitrap/FTMS centroid goldens, full historical noise values, natural spline values/derivatives, FWHM, mobility, boundaries and experiment selection |
| Iterative peak picking | 8 | HiRes seed/refinement noise separation, source centroid and integration conventions, width/spacing rules, aligned annotations, selected experiments and atomic limits |
| Internal iterative conventions | 2 | Fourteen independent refinement configurations and deterministic ties; strict-next-sample association, asymmetric recenter search, original seed priority and rounded centroid storage |
| Independent iterative reference | 4 | Complete synthetic result, original Orbitrap/FTMS raw-region sums and centroid rounding, omitted annotations, experiment behavior and checked failures |
| Iterative profile workflow | 2 | Noise estimation → centroiding → local peak filtering, independent raw sums, preserved metadata and exact mzML width/intensity array round trip |
| Window filtering | 10 | Source 56/30 counts, exact triangle indices, 480 independent cases, sliding/jumping boundaries and order, duplicate membership, stable ties and atomic resource failures |
| Iterative mean noise | 7 | All 2,526 historical source outputs, hand-computed fixed-denominator clipping, strict window edge, histogram ceiling, signed/empty inputs, checked legacy percentile and exact work limits |
| Independent mean-noise reference | 3 | Fresh histogram rescans for 192 configurations, signed global statistics, effective bin ceiling and legacy-percentile f32 arithmetic |
| Chromatogram picking | 9 | Legacy/corrected methods, Gauss/SG smoothing, independent seed/boundary noise, raw sums, exact index regions, aligned arrays, overlap and checked failures |
| Internal chromatogram conventions | 2 | Source right-side closest-sample ties/end sentinel and sequential overlap midpoint assignments |
| Peak integration | 12 | Inclusive sampled sum/trapezoid/Simpson, all baseline choices, signed/empty/singleton inputs, shape metrics, duplicates, bounds and numeric/resource failures |
| Independent chromatogram processing | 12 | All 146 source trace samples and literal integration/background/shape/picker goldens, nonuniform/even Simpson behavior, f32 arithmetic, -1 sentinel exclusion and source noise/overlap conventions |
| Chromatogram workflow | 2 | Raw picking → exact full-precision boundaries → intensity sum/time-weighted area/baseline/shape → mzML peak-annotation round trip; f32 boundary rounding cannot change integration |
| EMG fitting | 8 | Source cutoff fit, three model branches, metadata and omitted arrays, zero bounds, best iterations, extrapolation, shared evaluation limits and numerical failures |
| Internal EMG conventions | 3 | All four analytic gradients against central differences, source training collection order and iRprop sign/zero/rollback rules |
| Independent EMG reference | 8 | Seven complete source fits and four full raw loss goldens, true relative parameter tolerances, ordered training losses, scalar tails/subnormals, branch boundaries and inclusive budgets |
| EMG integration | 6 | Both container types, all integration/baseline methods, unchanged supplied shape inputs, cropped and expanded spans, literal zero bounds and preserved inputs on errors |
| EMG workflow | 2 | Independent continuous area/centroid identities and cropped fit → mzML → sampled integration without unaligned parameter annotations |
| Deisotoping | 8 | Unknown/unequal precursor-charge source regressions, charge priority, ladders, intensities, annotations, disjoint membership and error atomicity |
| Poisson/KL deisotoping | 11 | Source defaults and threshold/top-N preprocessing, original indices and aligned arrays, both sharing policies on both algorithms, low-mass noise, bounded work and transactional errors |
| Complete deisotoping fixture | 3 | All 5,407 source input peaks retained; all 104 expected output peaks bit-exact, independent seed/charge mapping, shared membership, exact 103-peak disjoint subset and retained metadata/arrays |
| Independent deisotoping review | 8 | Adjacent-f32 KL thresholds for sizes 2–7, mixed-precision accumulation, longest/highest-charge selection, selected counts/sums, inclusive ppm endpoints, nearest ties and exact precursor arithmetic on both methods |
| Spectrum comparison | 19 | Source alignment/scoring/binning goldens; 1,000 compact-DP/reference-map comparisons; directed ppm and weighted rounding boundaries |
| Text formats | 12 | DTA/FASTA/MGF parsing, streaming/round trips, invalid input, writer preflight/flush errors and malformed-input corpus |
| mzML | 14 | Independent/upstream fixtures, codecs/precision/units, metadata subset, limits, malformed XML/arrays and independent XSD checks |
| Independent format/processing review | 5 | Regressions for nonfinite filtering, empty resampling, full-range charges and metadata key loss |
| Independent mzML review | 3 | Conflicting scientific CV fields, reserved names and forbidden XML characters |
| Independent mzML auxiliary reader | 10 | Both float/signed integer widths, exact integer range, empty ASCII elements/placeholders, malformed/unsupported metadata, variable string decompression and cumulative byte/element/array limits |
| Independent mzML auxiliary writer | 4 | Validation before first output, exact source binary encodings, zero-length distinctions, type/name order and metadata, all four compression/empty-record XSD cases |
| Synthetic identification workflow | 1 | Embedded FASTA → fixed modification → digestion → literal observed fragment matching → peptide evidence/annotations → protein coverage/modifications; precursor consistency and retained metadata |
| Workflows | 4 | FASTA → digestion/fragments → MGF; DTA → processing → MGF; modified peptide → isotope envelopes → deisotoping → alignment; profile → centroid → consensus |

Spectrum annotation and ion naming received separate source, numerical and
integration reviews. The four new integration suites contain 30 tests, with one
additional sequence-ownership regression and four private grammar/work tests.
Original source fixtures contain thirteen measured peaks and seventeen scalar
or label assertions. Eleven independently derived masses/errors distinguish
rounded source references, including its zero-MSE assertion, from actual
nonzero error statistics. Exact binary examples independently check padding,
sample deviation, quartile indices and current ratios. See the
[annotation support](SPECTRUM_ANNOTATION_SUPPORT.md),
[reference review](SPECTRUM_ANNOTATION_REFERENCE_REVIEW.md) and
[fixture provenance](../tests/data/spectrum_annotation_provenance.json).

Native annotations preserve source last-match array behavior, duplicate
matched-only records, parameter no-ops, sort order and raw absolute precursor
tolerance in ppm mode. Small-list quartiles and undefined enabled statistics
have documented safe policies. Workflows verify modified fragments from
independent formulas, complete custom-chemistry identity, and XML transport.

Shared theoretical and alignment counters now persist across candidate hits;
coarse convolution also retains one allowance across theoretical envelopes.
Focused exhaustion tests verify that calls cannot restart these counters.
Anonymous modification strings are shared immutably, with independent slices
and setters retaining their owned lifetime, value and ordering semantics.

Peptide physicochemical utilities received independent source-table and numerical
reviews. Five integration suites add 35 tests, and a private pI test verifies its
shared evaluation budget. All 200 public AAindex values, 182 hydrophobicity cells
and 59 pKa constants are checked. Private gas-basicity tables feed independent
source expressions for every canonical singleton and all 400 ordered pairs.
Thirty original scalar assertions retain their source spelling and tolerance
context. A separate 100-case grid uses 90-digit Decimal arithmetic, including
low temperatures, tied maxima, subnormal products and the largest finite inputs.
See the [property support](PEPTIDE_PROPERTIES_SUPPORT.md),
[reference review](PEPTIDE_PROPERTIES_REFERENCE_REVIEW.md) and
[fixture provenance](../tests/data/peptide_properties_provenance.json).

Ordinary gas basicity preserves source evaluation order. Stable overflow retries,
the underflowed-product limit and the empty-sequence identity are documented
native numerical differences. Terminal annotation tests preserve source pI
suppression without inventing PTM-specific pKas; hydrophobicity and gas basicity
continue to use parent residues. Digestion and identification interchange tests
exercise those rules together. The current complete Rust 1.96 run includes all
110 integration suites and 46 private unit tests. Three doctests pass on both
compilers, alongside the focused RNA processing checks detailed above.

Owning isotope streaming and custom binary64 populations received independent
source and numerical reviews. The source ordered fructose support (2,548 states),
fourteen literal configurations and insulin's 10,000-state prefix pass alongside
threshold and custom natural-table materialization checks. Custom weights retain
extra binary64 precision; hand-derived products preserve separate equal-mass
configurations. Raw formula charge is ignored while the high-level generator's
hydrogen-adduct rules and all prior fine-spectrum references remain unchanged.

Boundary tests distinguish finite logarithms from underflowed raw probabilities,
check late overflow and fused errors, and validate every custom cell even for
zero-count populations. A 10,000-category case succeeds when coverage stops at
the first configuration, then fails with a checked work error when its raw stream
requests more. Raising its threshold avoids expanding an ineligible suffix.
A separate stream drains 103,041 configurations although materialization rejects
that full support. Threshold tests retain a returned raw probability exactly and
reject the next larger binary64 value; materialized logarithmic cutoff behavior
is preserved. Distinct equal-mass peaks and their configuration indices survive
selection and both mzML compression modes.

The materializer boundary assertion checks the analytic logarithmic condition
rather than assuming a particular platform's `ln`/`exp` rounding direction. Both
compilers pass it alongside the raw-value equality assertions in the full matrix.

The existing search is shared without a new dependency. Iterator cutoffs use
returned raw values, with logarithmic fallbacks for relative underflow or an
unrepresentable mode. Source layers, performance hints and untrimmed layer-based
membership remain outside the current surface. Streaming retains a frontier and
visited set, so atom, work, state and memory limits bound its lifetime.

Native fine-isotope enumeration received independent source, mathematical and
theoretical-spectrum integration reviews. All 44 selected class-test count
assertions, the 14-row high-precision fructose table and six bromine configurations
are checked separately from older loose/header illustrations. The 19,615-state
insulin boundary verifies accumulation of stored f32 probabilities; full-support
tests retain configurations whose stored probabilities underflow to zero. A small
exhaustive combinatorial oracle checks masses and probabilities independently of
the native heap algorithm, including separate fixed labels and natural H adducts.

The theoretical integration reproduces the source's 10-, 5-, 50- and 12-peak
fine-spectrum cases and checks retained terminal chemistry, source loss formulas,
isotope annotation alignment and natural-H mass/charge conventions. Loss intensity
products now remain f64 until final storage, avoiding an intermediate f32 rounding
in both fine and coarse envelopes; existing coarse regressions still pass. A
105-prefix test fails specifically at the shared fine-work limit while each
single envelope is valid and total possible output is only 105 peaks. Both this
failure and late atom-limit errors leave appended observations unchanged. Modified
fine spectra preserve their annotations and precursor acquisition data through
both mzML compression modes.

This validates materialized native configurations with explicit deterministic
ties, inclusive threshold equality and checked resource limits. It does not
claim identical IsoSpec layering or untrimmed layer-dependent selections, or
an executed C++ differential comparison. No dependency was added.

Precursor purity received independent numerical, source and acquisition-format
reviews. Exact binary64 coordinates and binary32 intensities from the original
mzML support seven source scalar/map cases; eight SPS cases and independent fuzzy
arithmetic are also checked. Rounded upstream literals are kept separate from
independently calculated results. The full original fixture reproduces the same
scores through the native mzML reader after an explicit in-memory encoding
declaration change; the bundled source bytes remain unchanged. All 19 activation
methods, four mobility quantities, distinct isolation targets and spectrum
references are covered, including populated writer output checked by xmllint.

Review found full spectrum validation repeatedly traversed unused annotations
outside the purity work budget. Purity now validates only consumed numerical
fields; a 10,000-child/10,000-placeholder regression verifies this bound. Scoped
precursor validation remains explicit before fuzzy empty-parent shortcuts.
Malformed container nesting and loss of CV-list ordinary metadata on export are
also rejected. The writer rejects a negative effective isolation target before
output, including selected-m/z fallback; spectrum/chromatogram regression cases
prevent producing a document that the reader would reject. Native source-undefined iterator/arithmetic cases return checked
errors or documented end-candidate fallbacks; finite source overcount and RT
extrapolation are retained. No C++ binary was executed for these comparisons.

Theoretical-spectrum extensions received independent source and integration
reviews. Analytical references cover internal fragment masses, counts, source
start/end conventions and first-residue loss omission. Literal source immonium
masses and rounded activation tables are checked separately; the compact helper
has an independent bit-exact f32 analytical oracle, rather than captured C++
output. Coarse-envelope independence, custom mass tags, retained termini,
annotation alignment and error atomicity are covered across the focused suites.

Review found repeated custom loss declarations could consume unbounded work
before deduplication and distinct losses could fill templates before peak
preflight. Declaration visits now receive a shared preflight estimate, and
unique loss-template entries are charged before formula cloning. Regression
tests verify both limits for ordinary and internal generation while preserving
existing output. Fine isotope support was added subsequently and is covered above.

Owned modification records, OBO loading and crosslink lookup received independent
source and consumer reviews. All 148 projected XLMOD records match exact source
accession/site order, names, synonym sets and f64 mass bits. Synthetic PSI-style
cases verify all-specificity UniMod aliases, absent targets and absolute formula
precedence. The historical PSI-MOD snapshot remains unbundled; its source hash
is verified separately from packaged fixtures.

Review corrected source empty-sequence export, literal empty versus explicit
zero-formula behavior and repeated alias expansion work accounting. A further
consumer review found same-name custom chemistry being merged in protein
observations, sequence duplicate filtering and conflict keys. Complete chemical
value ordering now preserves those distinctions and agrees with equality,
including signed zero. Peptide identity keys retain owned sequence values.
Source-intentional textual keys in other identification algorithms remain.
Caller-owned registry tests verify record release after the last owner is dropped,
absolute-formula restoration of unknown residues, checked vocabulary/mass export,
and exact idXML chemistry validation before output.

Modified-peptide generation and definition sets received independent source and
integration reviews. Literal cases include 26, 71 and 199 variants, supplemented
by a complete small weighted-site oracle. Source maximum-one/general terminal
placement and formula-derived versus declared mass distinctions are preserved.
Review found and corrected inference merging of chemically different records
sharing a full ID; conflicting records now fail atomically in either order.
Custom formula-free records retain declared/absolute monoisotopic mass through
generation and fragment calculations, with explicit unavailable composition.
The source no-change rule and formula precedence are checked separately.
Normal generated variants and search definitions round-trip through idXML;
unrepresentable typed terminal states are rejected before any writer call.

Iterative peak picking received an independent source and numerical review.
The source class test has no numerical assertions; fourteen derived refinement
configurations explicitly cover source control flow and mixed precision. Two
real input profiles independently verify centroid rounding and raw-region sums;
they are not treated as iterative output goldens. Window selection matches
literal upstream counts and 480 independent small cases. Iterative mean noise
matches all 2,526 historical outputs at the source tolerance and fresh histogram
rescans over 192 configurations. The combined profile workflow verifies that
filtering preserves the picker’s aligned integration and width arrays through
mzML. The implemented source quirks and native checked-error policies are
recorded in the corresponding support documents.

EMG fitting received an independent numerical review. All twelve branch
expressions across four gradients match the source arithmetic, with separate
central-difference checks. Seven source fits retain their literal minute/second
coordinates, f32 container intensities and true relative parameter checks;
four source loss assertions distinguish full raw-f64 loss from training loss.
Derived training-order fixtures and erfc tail/subnormal tests independently
check operation ordering and special-function behavior. The strict numerical
error policy is documented, including failures after a previous finite best.
An additional mathematical workflow verifies total EMG area and centroid, then
checks cropped reconstruction and exact fitted-trace mzML interchange.

Independent reviews inspected the numerical modules and their source conventions.
Chromatogram integration and picking preserve literal source reference values,
mixed f32/f64 operation order, sampled boundaries and shape conventions. Reviews
caught an uncharged Gaussian coefficient-table allocation and verified the
corrected work limit. The combined workflow additionally preserves exact original
sample indices when f32 boundary metadata cannot represent their coordinates.

Auxiliary mzML arrays received independent reader and writer reviews. Tests
verify binary bytes, signed integer limits, ASCII terminators and empty elements,
validation before any write/flush, and cumulative decoded-storage bounds. A
review regression rejects units attached to encoding/compression terms, preventing
unrepresented metadata loss. The provenance audit distinguished the mzML schema
fixture's CRLF bytes from the immutable source's LF bytes and verified identical
normalized schema content; both hashes are now recorded.

Poisson/KL deisotoping received separate numerical and data-flow reviews. The
complete real-spectrum fixture revealed two source clusters sharing a heavy
isotope; the native default preserves both, while an explicit disjoint policy
omits only the independently identified overlapping seed. Source binary arrays,
decimal fixture values and all 104 seed/charge mappings were independently
verified. Review also corrected precursor operation order and low-mass noise
handling, and checked mixed f32/f64 KL acceptance boundaries.

Sequence review checks the distinction between known mass and known formula,
annotation ownership, position-preserving round trips, and source numeric lookup
precision. Regressions cover two large-mass cancellation errors: independently
accumulated suffix fragments preserve a small alanine fragment beside a large
tag, and equal/opposite terminal deltas preserve the peptide residue mass.

Protein inference, peptide indexing, conflict resolution, origin splitting and
the combined protein workflow received independent scientific reviews. The
matcher oracle includes 121 raw cases, 17 decoy-inference cases and 15 enzyme
boundaries. Review regressions address undefined-negative greedy scores,
original search-engine recovery, search-parameter metadata placement and
pre-clone sequence limits.

Identification scoring/filtering/FDR and idXML received separate peer reviews.
Regressions reject target/decoy label replacement, aliased analysis-result
indices, malformed declarations/processing instructions and encoded group lists
over their configured limit. Other reviews checked digestion context, typed
metadata, evidence coordinates, transformation inverses and atomic container
edits. The synthetic identification precursor was corrected against an
independent mass calculation.

All 84 elements and 283 isotope data points were compared against the source
declarations. All 2,341 inventory entries were verified against immutable Git
objects at the pinned revision; one CRLF checkout conversion is recorded
separately. Checks establish the documented surface; they do not count fully
ported C++ classes.

## Limits of this evidence

The C++ reference was inspected, **not built or executed**. These are
source-derived reference tests, independently reconstructed numerical
expectations, fixture/schema checks and Rust invariants. They are not runtime
differential tests against OpenMS or proof of full workflow equivalence.

Linux, macOS, Windows and Rust 1.85 CI have passed for the preceding published
commit, as recorded above; newer local changes have separate validation records.
Linux CI installs xmllint; on other hosts the schema tests explicitly report when
that optional executable is unavailable. It was available and executed locally.
Rust 1.85 compatibility was also exercised locally.

There are no runtime performance/memory benchmarks, external application
interoperability certification or full PSI controlled-vocabulary validation.
The capability documents describe current exclusions and deliberate source
corrections. No C++ source, contrib tree or vendored dependency was modified.
