# SQLite consumers and spectrum access: S2 plan

This source-reviewed plan follows the S1 sqMass handler. It targets
`OpenMS4-core` revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4` and does not
claim S2 implementation or executed C++ evidence. The shared
[issue log](../OpenMS_CPP_ISSUES.md) records confirmed findings separately from
unconfirmed contract candidates.

## Scope and implementation order

First finish S1 validation. Then implement `SqMassFile` loading, storing and
streaming alongside `MSDataSqlConsumer`, reusing the existing
`interfaces::MSDataConsumer` trait and the handler's storage and metadata code.
Complete the OpenSwath data/access interface before presenting
`SpectrumAccessSqMass` as a replacement usable by TOPP consumers.

`SqMassFile` also exposes XIC Parquet conversion. That method needs
`MSChromatogramParquetConsumer` and `LightTargetedExperiment` metadata support.
Loading, storing and streaming alone leave the public header partial. An
unsupported placeholder does not complete a method implemented by C++.

## SqMassFile

Map the default constructor/destructor, `MapType` alias, `SqMassConfig`,
`setConfig`, `load`, `store`, `transform`, and `convertToXICParquet`.
Configuration defaults are full metadata, lossless arrays and linear mass
accuracy **−1**. Preserve the documented accuracy sentinel: when Numpress is
selected, nonpositive target accuracy uses ordinary fixed-point estimation.
Replacing it with 0.0001 silently changes the configuration.

`store` takes the experiment's SQL run ID and replaces the destination in C++.
The native wrapper should stage complete output and replace the destination
only after success. Document this behavior and the chosen metadata-loss policy.

`transform` sends expected counts and experimental settings before spectrum
callbacks, followed by chromatogram callbacks. Both source skip flags are
explicitly ignored. Read bounded batches without collecting all primary arrays.
Avoid the source's empty final batch at zero or an exact multiple of 500
(CPP-200). Use actual record IDs so external databases with gaps work. Specify
callback cancellation, mutation and failure behavior and the consistency
boundary if separate batches use separate transactions.

XIC conversion defaults the source filename to the input path, passes run ID and
transition metadata to the Parquet consumer, streams, then finalizes. Missing
transition metadata must produce a checked error. Completion requires tests of
real Parquet output and metadata, not merely the delegation call.

## MSDataSqlConsumer

Map construction with buffer size, full metadata, compression and mass accuracy;
`flush`, `addRun`, `setRunId`; and all four consumer callbacks. The constructor
creates the database. Source consumption copies a record into a buffer, clears
its caller-visible primary and auxiliary arrays, retains metadata when enabled,
and flushes both buffers when either reaches its threshold.

Use an owned handler, two buffers and retained descriptive metadata. Validate
configuration before file replacement or allocation. Bound buffer records,
primary values, metadata and total retained bytes. Flush may continue accepting
data afterward; define whether its two record kinds commit together or
separately, and retain failed buffers predictably.

Provide explicit `finish() -> Result` for confirmed persistence. Rust destruction
cannot report I/O failure, so it must not be the only completion mechanism.
Preserve received experimental settings deliberately: the source callback is a
no-op. Define run-switch behavior for already buffered records. Source `addRun`
writes an empty snapshot and suppresses the final accumulated snapshot, so run
registration and metadata finalization need separate states. Multiple RUN rows
also conflict with the handler's single-run reader; do not advertise a multi-run
round trip without resolving that contract.

## SpectrumAccessSqMass and OpenSwath dependencies

Map the three constructors (all records, absolute subset, and parent-relative
subset), copy/light clone, destruction, individual spectrum and metadata access,
bulk spectra, RT selection, both counts, and both chromatogram accessors. The
source deliberately throws `NotImplemented` for chromatogram payload and native
ID access; checked native errors can preserve that declared scope.

Keep visible positions separate from SQL IDs. Empty source selections mean all
records or inheritance, rather than an empty view. Preserve reversed selections
and repeated positions, or explicitly document checked rejection. S1 returns
sorted unique IDs, so bulk view access needs an ID-to-record mapping followed by
view-order reconstruction. Validate every index and populate metadata indices.
Clones should retain read configuration and view IDs while opening independent
connections. Define the source's special zero-delta RT behavior explicitly.

The dependency is **OpenSwath** `ISpectrumAccess`, not the unrelated
`OpenMS::Interfaces` reader/writer interfaces. Minimal shared closure includes
binary arrays, spectrum/chromatogram arrays and metadata, the access trait,
light cloning, drift filtering and nearest multiple-spectrum helpers. Missing
mobility and mismatched array lengths require checked errors. Define counts and
ties for multiple-spectrum selection. Other accessors and scoring algorithms
are not required to complete this bounded dependency group.

## Test and evidence responsibilities

| Source class | Sections | Active assertion macros |
| --- | ---: | ---: |
| SqMassFile | 6 | 87 within sections, plus 8 helper macros |
| SpectrumAccessSqMass | 7 | 39 |
| MSDataSqlConsumer | No dedicated class test found | No translated assertion claim |

Exclude commented assertions; static counts are not runtime counts. Reuse the
retained sqMass and paired mzML fixtures. Check all primary values with source
tolerances, precursor/product literals, settings and run IDs. Add independent
tests for configuration sentinels, batching at 0/499/500/501, callback order and
cancellation, finishing failures, run transitions, reversed/nested/duplicate
views, ID gaps, cloning, metadata indices and unsupported mobility. Neither
section mapping nor native execution establishes a full C++ differential.
