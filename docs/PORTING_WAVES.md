# Resumed porting waves

Prioritize working FileInfo and FeatureFinderCentroided executables, followed
by PeakPickerHiRes, using the existing core and five TOPP integrations. Use core revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4` until an
explicit source-refresh checkpoint updates the inventory and manifests.
The [completion ledger](CORE_SDK_COMPLETION.md), not a module-name list, remains
the record of outstanding public operations.

## Completed checkpoint: FORMAT integration

[The fifteen-module wave](FORMAT_WAVE_SUPPORT.md) combines the preserved
`fmt-pkg-*`/`fx-*` branches with source, documentation and regression review.
It includes minimum-feature CI slices and explicit partial-header assessments.
The [validation record](VALIDATION.md) identifies executed checks. Known gaps
remain work for subsequent waves, including mzTab exporters, format-specific
validators, qcML metrics and missing FileHandler dispatch.

## Recovered Claude checkpoints (2026-09-13)

The resume began at clean `45c424a`, with all 42 retained Claude worktree heads
already merged. `068a1f5` preserved interrupted SQLite S1 work; `3b8ca94` and
`4e8934b` integrated numerical, comparison and system waves A/B. `d396f42` adopted
rustfft with the scalar planner; `45c424a` records the crate-first dependency
policy. Default Rayon parallelism and the serial/parallel bitwise contract remain.
These commits are preserved; their reported test counts are historical author
reports. The resume's whole-tree remote run is summarised in the
[resume assessment](RESUME_ASSESSMENT_2026-09-13.md); its four
`tests/system_process.rs` failures were a test-harness `ETXTBSY` race, since
fixed. The final resume snapshot's full gates are recorded in
[VALIDATION](VALIDATION.md).

## Current priority: early TOPP bundle

The [early TOPP plan](EARLY_TOPP_BUILD_PLAN.md) supersedes the previous
storage-first sequence. Build FileInfo for mzML/featureXML while porting the
FeatureFinderAlgorithmPicked helper and trace-fitting chain in parallel.
FeatureFinderCentroided is the principal scientific target; PeakPickerHiRes
follows to enable profile mzML -> centroided mzML -> featureXML. Keep the five
existing tools building throughout.

Preserve and validate the current SQLite checkpoint. Its process-test failures
are diagnosed and fixed: an `ETXTBSY` race between parallel tests. Neither the complete SDK nor SQLite S2/OSW/OMS is a gate for
the early tools. Build real supported paths early, explicitly label preview
limitations, and expand toward full tool contracts using the source test matrix.
The detailed plan defines dependencies, ownership and acceptance cases.

## Queued storage waves: SQLite family

The optional `sqlite` feature and pinned `rusqlite` dependency are already in
Cargo. Keeping SQLite optional is a dependency boundary, not a scope exclusion.

S0 now has a [native connector](SQLITE_CONNECTOR_SUPPORT.md), all eight source
class-test sections mapped, nineteen native regression tests and an adapted
C++ defect probe. The [validation record](VALIDATION.md) records the executed
integration checks. S1 code is integrated. The current closeout adds precursor metadata transport,
write-limit rollback and schema-guard regressions, with fresh remote tests and
source review. The requested Fable review remains pending authorization; no
review approval or full database-file adapter closure is inferred from S1.

| Stage | Work | Dependency and verification boundary |
|---|---|---|
| S0 | `SqliteConnector` public operations | Three open modes, table/column queries, row counts, statement and blob binding; port all eight source class-test sections. Raw helpers in `SqliteConnector_impl.h` are private implementation, not another public-header closure. |
| S1 | `MzMLSqliteHandler` and `MzMLSqliteSwathHandler` | Pin schemas, compression and identifier rules; check malformed/truncated blobs and read-only behavior against retained sqMass inputs. |
| S1 parallel | OSW records and `OSWFile` | Map `OSWData`, export and inference value/config dependencies before claiming the file API; parameterized SQL and schema checks must preserve identifiers and relationships. |
| S2 | `SqMassFile`, consumer and `SpectrumAccessSqMass` | Reuse S1 storage, existing experiment/consumer interfaces and Numpress; test random access, full load, transforms and metadata retention. |
| S2 parallel | `OMSFileStore`, `OMSFileLoad`, `OMSFile` | Map identification-graph, FeatureMap and ConsensusMap persistence; test foreign-key relationships, version handling and round trips against retained source databases. |
| S3 | Integrate and review the whole family | Wire supported FileHandler operations, run all features plus SQLite-only/minimum combinations, compare retained source results and update documentation before committing. |

The [S1 implementation plan](SQLITE_STORAGE_PLAN.md) pins the schema, fixtures,
API mapping and regression matrix. The initial findings are recorded as
CPP-190 through CPP-202; the auxiliary-array recovery contract remains an
unconfirmed candidate. The adapted SWATH probe reproduces CPP-196/197/202;
CPP-203 adds unchecked query-step errors, reviewed from source and tested in
Rust. The shared issue log distinguishes native handling from proposed C++
fixes; no upstream fixes are claimed.

The [S2 consumer and spectrum-access plan](SQLITE_CONSUMER_PLAN.md) inventories
the queued public APIs, Parquet and OpenSwath prerequisites, and streaming tests.
Its [source-review manifest](../tests/data/sqlite_s2_review.json) records
CPP-204 through CPP-217 with explicit defect/candidate labels. No S2 runtime
implementation or execution is claimed by that preparation.

Before each stage, inspect the pinned header, implementation, registrations and
class tests. A dependency discovered during implementation becomes an explicit
prerequisite; do not silently omit its methods. Distinct agents may own leaf
modules in parallel, while the integrator owns shared registration, provenance,
coverage and the C++ issue log.

## Following waves: full core, CLI and TOPP closure

Finish open FORMAT contracts and remaining scientific/core dependencies in the
ledger, prioritizing dependencies used by multiple TOPP tools. All vendor,
container, database and binary formats stay in scope. Keep public docs explicit
when an adapter cannot preserve a field or validate a schema.

For CLI work, pin `openms4-cli` separately, map argument/parameter registration,
configuration, diagnostics and exit behavior, then reuse that layer across TOPP
ports. For each TOPP tool, pin `openms4-topp` and its test-data revision, port the
actual upstream workflow, and compare outputs at the upstream tolerances.
The current five validated workflows do not imply closure of the other tools.

Run compilation and test batches on `kim` in separate `/scratch` workspaces.
Use current Rust and MSRV 1.85, strict clippy/rustdoc, doctests, feature checks,
and provenance/coverage gates. Request bounded Claude Fable 5.1 reviews with
concrete source context, reproduce actionable findings, and retain dispositions;
a model review alone is not correctness evidence. Commit and push tested wave
checkpoints without presenting the unfinished SDK as complete.
