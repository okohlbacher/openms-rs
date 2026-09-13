# Early TOPP builds: FileInfo and FeatureFinderCentroided

This plan supersedes the storage-first execution order in the September 13
resume assessment. The objective is an early, useful Rust tool bundle, led by
**FileInfo** and **FeatureFinderCentroided**, while the complete SDK remains the
long-term scope. Completing SQLite S2, every FileHandler format, all CLI base
classes or entire scientific domains is not a prerequisite for these builds.

## Deliverables and order

Keep the five existing tools building: BaselineFilter, DTAExtractor,
MapNormalizer, SpectraFilterWindowMower and MzMLSplitter. Add FileInfo first,
FeatureFinderCentroided as the main scientific target, and PeakPickerHiRes as
the next complementary tool. This yields an eight-tool initial bundle and the
useful chain:

```text
profile mzML -> PeakPickerHiRes -> centroided mzML
                                      |
                                      v
                            FeatureFinderCentroided -> featureXML
                                      |                    |
                                      +------ FileInfo ----+
```

PeakPickerHiRes is not required when the input is already centroided. Its port
must not delay FeatureFinderCentroided. FileInfo should inspect both the input
and resulting feature map, making the earliest bundle useful for real work.

An early build means a working, tested command on real input with a documented
capability table, not a binary that only prints help. Unsupported formats or
options must fail explicitly. A preview with limited input modes remains
partial in the ledger; it is not advertised as a complete upstream tool.

## Source basis and current gaps

The dependency trace uses Core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`,
`../OpenMS4-tests/packages/topp/src/{FileInfo,FeatureFinderCentroided,PeakPickerHiRes}.cpp`,
and the test registrations in
`../OpenMS4-tests/packages/test-data/topp/CMakeLists.txt`. CLI, TOPP and test data
are independently pinned in `tests/data/topp_cli_provenance.json`. Before coding,
hash the selected source/tests against those pins and record any deliberate
source refresh as a separate checkpoint.

| Dependency group | Current evidence | Work needed for early tools |
|---|---|---|
| TOPP registration, INI, exit codes | Existing framework and five tools | Reuse; close the actual flags, subsection defaults, processing annotation and output behavior these tools call. Wire requested thread policy. |
| mzML, Feature/FeatureMap, featureXML | Native readers/models/writers; featureXML header partial | Verify required metadata, hulls, subordinates, processing, primary-run paths and seeds. Complete only missing operations on these paths first. |
| FileInfo library | `FORMAT/FileInfo.h` unmapped | Port structured results/options, computation and text/TSV rendering; keep executable thin like C++. |
| FileInfo statistics and peak types | Numerical helpers and spectrum-type code exist | Verify SummaryStatistics conventions, empty ranges, type estimation and formatting against the source reporter. |
| Picked feature algorithm | `FeatureFinderAlgorithmPicked.h` unmapped | Port the actual centroided-peptide algorithm; existing FeatureFindingMetabo is a different workflow. |
| Picked helper types and fitters | HelperStructs, TraceFitter, GaussTraceFitter and EGHTraceFitter unmapped | Implement the helper contracts and both trace models before claiming full algorithm support. Audit existing LM/isotope/math routines for reuse. |
| FAIMS and overlap | IMTypes complete; IMDataConverter requires review; FAIMSHelper and FeatureOverlapFilter unmapped | Implement CV grouping, seed filtering and cross-CV merging used by the wrapper. Existing mobility containers alone do not supply these operations. |
| PeakPickerHiRes | Native processing implementation; header evidence requires review | Audit full experiment/metadata behavior and add the TOPP wrapper plus upstream workflow comparisons. |

Do not replace the picked algorithm with mass-trace detection/FeatureFindingMetabo,
or assume GaussFitter implements GaussTraceFitter: the objectives, data and
initialization contracts differ. Reuse compatible internals after comparing the
actual call paths. Follow the crate-first policy for external numerical libraries;
any new dependency is integrated centrally and measured against source results.

## Milestone 0: preserve the baseline and start focused builds

Preserve the current S1 changes and their validation evidence. The four SYSTEM
process-test failures first seen under Rust 1.85 are resolved: they were an
`ETXTBSY` race between parallel tests in `tests/system_process.rs`, reproduced
under both 1.85 and 1.96, and fixed by running that binary's tests serially (see
the update in the resume assessment). The full gates still have to pass on the
final checkpoint before it is published.
Do not expand this work into SQLite S2 or unrelated network repairs before the
target tools can run.

Inventory only the shared CLI methods the two wrappers execute. In particular,
FeatureFinderCentroided needs algorithm subsection defaults, seed input, `-force`,
`-test`, processing annotations and debug-dependent output. FileInfo needs
optional text/TSV outputs, quiet progress, forced type and mode-dependent errors.
Reuse existing facilities rather than completing all descriptor writers,
ToolHandler discovery or unrelated search-engine base classes first.

As soon as the real implementations link, add Cargo binary registrations and
focused CI builds. The initial mzML/featureXML bundle should use the narrow
`mzml,paramxml,featurexml` feature set; add optional formats/validation features
only as their tested modes land. Confirm the final feature closure rather than
assuming this initial set is sufficient for every mode.

## Milestone 1: useful FileInfo preview

Port the reusable FileInfo result/options/report layer and the thin CLI together.
First support mzML and featureXML with type detection/forced type, counts and
ranges, metadata (`-m`), processing (`-p`), summary statistics (`-s`), stdout,
`-out` and `-out_tsv`. The featureXML report must work on the eventual feature
finder output. Preserve the distinction between absent and zero-valued ranges.

Reuse existing DTA and DTA2D readers for small, fast reference cases. Add detailed
peak listings (`-d`), corrupt-data checks (`-c`) and indexed-mzML checks (`-i`) as
bounded follow-ons; existing index/validation code should be audited and reused.
Schema/semantic validation (`-v`) is a separate capability: never turn a missing
validator into a successful validation report.

First acceptance cases are upstream FileInfo_1 (DTA), _2 (DTA2D), _3
(featureXML with metadata/statistics/processing) and _9 (mzML with those flags),
plus TSV, empty-map, unknown-type and output-error regressions. FileInfo_11/_12
exercise index behavior; recover their expected exit statuses from the test
registration. Pin source formatting and upstream FuzzyDiff tolerances/whitelists;
do not widen tolerances or ignore arbitrary metadata differences.

The inspected test registration has 19 FileInfo invocations (1–7 and 9–20),
including legacy XML, consensusXML, idXML, mzIdentML validation, transformation
XML, FASTA and FAIMS. Extend these after the first useful build. Full FileInfo
also retains sqMass, PQP and the other source-advertised formats in its backlog;
those dependencies do not gate the mzML/featureXML preview. Keep `FileInfo.h`
partial until every required branch and public result/helper contract is covered.

## Milestone 2: FeatureFinderCentroided scientific dependency chain

Begin this lane immediately alongside FileInfo; it is the longer critical path.

1. Port picked-helper seeds, isotope patterns and mass-trace collections, with
   source class tests. Verify the required isotope-distribution/averagine,
   IsotopeCluster, hull and FeatureMap operations already present in Rust.
2. Port TraceFitter, GaussTraceFitter and EGHTraceFitter. Compare residuals,
   Jacobians, initialization, stopping, fit failures and integrated areas using
   executed C++ probes before trusting whole-feature agreement. Default symmetric
   fitting uses Gaussian; `feature:rt_shape=asymmetric` uses EGH.
3. Port FeatureFinderAlgorithmPicked: preprocessing/scoring, seed selection,
   isotope matching, trace extension, fitting/cropping, quality checks and feature
   overlap/selection. Preserve seed ordering, tie handling, charge hypotheses and
   numerical operation order where the source contract needs it.
4. Connect the actual executable pipeline: load MS1 only, filter nonpositive
   intensities, reject empty MS1 input, reject per-peak mobility, enforce the
   profile-input/`-force` rule, load featureXML seeds, forward algorithm parameters,
   and write featureXML with IDs, primary-run path and QUANTITATION processing.
   Preserve hull bounding-box/subordinate cleanup under the source debug rules.

The first working preview may explicitly support non-FAIMS centroided mzML.
It must reject unsupported FAIMS input rather than pool compensation voltages
silently. Ship that preview when it passes the retained workflow, while finishing
FAIMS grouping, CV-matched/unannotated seed handling, annotation and optional
cross-CV merge (`faims_merge_features`). Conditional vendor RAW input remains a
separate platform capability, not a reason to hold back ordinary mzML builds.

Acceptance starts with TOPP_FeatureFinderCentroided_1: the exact retained input,
INI and featureXML expectation registered at CMakeLists.txt:426–428. Its output
comparison permits `id=` differences; reproduce the source numerical comparison
policy explicitly. Also compare feature count, RT, m/z, intensity, charge,
quality, hulls, subordinates and processing metadata. This single source workflow
does not test the entire tool: add seeds, asymmetric fitting, empty input,
profile rejection/force, rejected per-peak mobility and FAIMS cases. Compare
thread counts using the project's serial/parallel bitwise contract.

Acceptance for a complete supported configuration requires these wrapper modes,
not just the happy-path fixture. Preserve missing modes as explicit partial
status until verified. Both library-level and tool-level comparisons are needed.

## Milestone 3: PeakPickerHiRes and an end-to-end bundle

Audit the existing picker rather than reimplementing it. Add its CLI defaults,
MS-level selection, metadata preservation and output annotations; use the six
registered upstream workflows, including low-memory modes, as the closure list.
An initial regular-memory mzML mode may land earlier, with the other modes
explicitly partial. This third new executable can progress when a lane is free;
it must not take the numerical implementation lane away from the requested
feature finder.

Run the real chain on retained and representative data:
FileInfo(profile) -> PeakPickerHiRes -> FileInfo(centroided) ->
FeatureFinderCentroided -> FileInfo(featureXML).
Retain reports, parameter files, feature outputs and resource measurements. Build
all eight executables in release mode and rerun the existing five workflows.
This is the first bundle milestone, not a complete SDK release.

## Parallel ownership and validation

Use three bounded lanes: (A) FileInfo/reporting and required format dispatch,
(B) picked-helper/trace-fit/algorithm work, (C) C++ reference execution
and fixtures. After common helper interfaces stabilize,
Gauss/EGH implementations can be split; start PeakPickerHiRes when capacity frees.
The integrator owns Cargo/module registration, shared CLI/kernel changes,
provenance, coverage and the C++ issue log. Avoid two agents editing the same
shared file or creating duplicate scientific models.

Build/test remotely on a freshly checked IBMI node, with scratch sources/targets
and durable Ceph evidence. Run focused current/MSRV tests as each slice lands,
then full gates before claiming a release-ready bundle. Scope independent model
reviews to these actual tool paths; retain findings and dispositions. Existing
Fable authorization/authentication constraints still apply, but do not prevent
local porting, native tests or C++ comparisons. C++ drivers stay in `../oracle/`
with hashed manifests, and source bugs/candidates retain stable log IDs.

## Work moved behind the early bundle

Broad SQLite S2/S3, OSW/OMS, full FileHandler format coverage, unrelated CLI
base classes/descriptor exporters and domain-wide audits follow the early tool
milestones. Pull a specific operation forward only when the traced tool path
needs it (for example SqMassFile for FileInfo's sqMass mode). FileConverter and
FileFilter remain useful later candidates, but their broad option/format surface
must not delay the named targets. All formats and remaining SDK/TOPP APIs stay
in scope; only the order changes.
