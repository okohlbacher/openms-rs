# Upstream Core SDK delta: 6bfc0e4 → 54a232f

Read-only audit, 2026-09-10. Base: `6bfc0e4711105f4eda2fea86812a83af7c7e791f`; new target: `54a232fe2cae9c590d5c997fa49d20e7769860fb`, repository `okohlbacher/OpenMS4-core`.

**Conclusion:** no changed scientific constants, datasets, identification-graph records, or current Rust f32/f64 arithmetic require a production update. The update fixes several C++ error/lifetime/portability paths, improves long-double formatting outside the native Rust value surface, and simplifies upstream builds/delivery. Keep the current Rust implementation, audit/update source-target compatibility metadata separately, and use the new failure behavior when porting the still-missing components. Do not relabel historical fixture hashes.

## Evidence and scope

The GitHub compare API reports **10 commits ahead, 0 behind**, and **158 changed paths: 12 added, 77 modified, 69 removed**. This is below the API's 300-file comparison ceiling. Exact changed paths and patches are available in the [upstream comparison](https://github.com/okohlbacher/OpenMS4-core/compare/6bfc0e4711105f4eda2fea86812a83af7c7e791f...54a232fe2cae9c590d5c997fa49d20e7769860fb).

Changes within scientific code comprise **five existing public headers**, one generated config template, and **22 implementation paths** (21 modified plus the new private `SYSTEM/ProcessWait.h`). **No registered public header is added or removed.** No public-header/source registration manifests change; the current 786-header scope remains applicable. `src/tests/class_tests/openms/executables.cmake` removes only the obsolete `Boost_dependent_tests` linkage list, not its listed test executables.

There are **24 modified class-test source files**, no changed class-test data files, and no changed `share/OpenMS` runtime resources. OpenSWATH algorithm code is unchanged; only `src/openswathalgo/CMakeLists.txt` changes. No identification graph, kernel, NASequence, residue registry, isotope data, XML schema, or vendored implementation changes appear.

The major behavior-fix commit is `56b181fe33933828be11ac95be1b702a48c55247`; the extended-long-double fix is `1b38a997784e957374a47402ddf42d22e058ac8e`. The final `54a232f` commit adds the MSVC `/utf-8` compiler flag. Other commits primarily change build, CI, packaging, or test path handling. Upstream commit messages report their own validation; this audit did not execute or independently certify those tests.

## Actual API and behavior changes

Paths in the next two tables are relative to `src/openms`.

| Public header/config | Change and Rust consequence |
| --- | --- |
| `include/OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h` | `AA::fromIndex` loses `noexcept`, allowing its existing invalid-index precondition to throw instead of terminating. Domain and return values are unchanged. Current native peptide indexing uses checked Rust representations; retain checked errors if exposing this low-level API later. |
| `include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h` | Adds the correct `override` marker to `getMetaValue`; no signature or runtime behavior change. |
| `include/OpenMS/FEATUREFINDER/MassTraceDetection.h` | Marks unused `mass_error_da_` field `[[maybe_unused]]`, explicitly retaining layout; ppm tolerance behavior unchanged. |
| `include/OpenMS/FORMAT/HANDLERS/UniProtXMLHandler.h` | Marks two unused state booleans `[[maybe_unused]]`, retaining layout; capture-based name selection unchanged. |
| `include/OpenMS/MATH/MISC/BSplineSmoothingSpline.h` | Marks unused degree field `k_` `[[maybe_unused]]`, retaining layout; cubic fitting unchanged. |
| `include/OpenMS/config.h.in` | Removes unused `ENABLE_UPDATE_CHECK` macro. No Rust scientific feature or value changes. |

| Implementation path | Actual delta and current-port impact |
| --- | --- |
| `source/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.cpp` | Holds the protein description in an owned `std::string` during regex matching, eliminating a dangling `c_str()` pointer. No numerical formula change. This normalization class is currently unmapped; its future Rust port should use a borrowed/owned live description. |
| `source/ANALYSIS/OPENSWATH/OpenSwathResultsExporter.cpp` | Checks every Arrow append/null append and array finalization through `ParquetFile` helpers instead of discarding status/using unchecked values. Column schema/order and normal data values remain unchanged. This exporter is currently unmapped; future implementation must propagate failures, including empty-output schema retention. |
| `source/FORMAT/ZipArchiveFile.cpp` | Uses current UTF-8 path helper; uses `zip_file_replace`; captures error text before freeing archive handles and uses `zip_discard` on failed modifications. Fixes use-after-close diagnostics and prevents error-path commit. This archive-writing class is currently unmapped. The separate `ZipIfstream`/`ZipInputStream` single-entry read contract is unchanged, so the retained ZIP-input brief remains applicable. |
| `source/SYSTEM/ProcessWait.h` (new), `source/SYSTEM/JavaInfo.cpp`, `source/SYSTEM/PythonInfo.cpp` | Replaces deprecated/unreliable Boost timed waits with a monotonic deadline and 10 ms nonblocking status polling, then collects exit status. Existing 30-second probe deadlines and caller termination behavior remain. JavaInfo/PythonInfo are unmapped; use safe native timed process handling when implementing them. |
| `source/DATASTRUCTURES/StringUtils.cpp` | On libc++ with extended `long double`, uses classic-locale streams to preserve values beyond double precision/range. Removes two redundant branches with identical outcomes. **float/double thresholds, precision and scientific/fixed choices do not change.** Current `param/value.rs` and list formatting are f32/f64 and need no change. Portable long-double support remains an explicitly documented existing boundary, not newly implemented Rust functionality. |
| `source/IONMOBILITY/IMTypesExperiment.cpp` | Uses `std::cmp_not_equal` for signed/unsigned MS-level comparison. Negative requested levels no longer accidentally equal a wrapped unsigned level. This specific class is not implemented; preserve typed/checked comparison in future work. |
| `source/CHEMISTRY/IsoelectricPoint.cpp` | Adds explicit empty residue-override maps to three pKa aggregate initializers. Previously omitted aggregate fields were already value-initialized empty. All numerical constants and computations are identical; native `chemistry/isoelectric_point.rs` needs no change. |
| `source/CHEMISTRY/ProForma.cpp` | Removes a counter that is never read. Parsing behavior unchanged; the native general ProForma gap remains. |
| `source/COMPARISON/PeakAlignment.cpp` | Changes eight matrix-loop counters from platform-dependent `long int` to `Size`. No recurrence/order change for representable bounded dimensions. Current native comparisons use bounded Rust indices; no arithmetic update required. The coverage ledger still does not claim the whole PeakAlignment header. |
| `source/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.cpp`, `source/ANALYSIS/QUANTITATION/ItraqConstants.cpp` | Matrix-loop counters become `Size`; constants and operations unchanged. Both classes remain unmapped. |
| `source/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.cpp`, `source/ANALYSIS/OPENSWATH/TransitionParquetFile.cpp` | Replace deprecated `u8path` calls with existing `OpenMS::to_path`; first file also drops unused overload/helpers. No normal output schema change. These exporters remain unmapped. |
| `source/ANALYSIS/OPENSWATH/OpenSwathPercolatorScoring.cpp` | Simplifies a self-assignment branch to an equivalent condition; score-name normalization unchanged. |
| `source/ANALYSIS/ID/FragmentIndex.cpp`, `source/ANALYSIS/ID/Percolator.cpp`, `source/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.cpp` | Remove unused private functions only; no public feature removal or runtime change. |
| `source/FORMAT/HANDLERS/ImzMLHandlerHelper.cpp` | Compiles byte-swap helpers only on big-endian hosts; decoding behavior unchanged. |
| `source/FORMAT/HANDLERS/UniProtXMLHandler.cpp` | Removes unused private `isEmpty` helper. |
| `source/FORMAT/ParamCWLFile.cpp` | Compiles private replacement helper only when TDL is enabled; output behavior unchanged. |

## New or strengthened source reference cases

* `AhoCorasickAmbiguous_test.cpp`: indices 0 and last valid produce A/V; the first invalid index throws a precondition error.
* `ConsensusMapNormalizerAlgorithmMedian_test.cpp`: heap-allocated protein description containing 4096 `x` characters matches the intended description/accession regex and rejects mismatches.
* `OpenSwathResultsExporter_test.cpp`: checks populated/null optional columns, 64-bit alignment ID `1234567890123`, bool fields, and an empty Parquet result retaining all 59 schema columns.
* `ZipArchiveFile_test.cpp`: invalid UTF-8 entry-name addition returns the libzip diagnostic safely.
* `String_test.cpp`: extended `9007199254740993.0L`, long-double maximum and minimum survive textual round trips. These are **not f64 goldens**.
* New `tests/runtime/ProcessWait_test.cpp`: success/failure exit statuses and a 20 ms timeout returning control while the child remains alive.
* `ChromatogramProcessor_test.cpp` and `DefaultChromHandler_test.cpp`: replace vacuous unsigned-size checks with specific empty/mapped counts, native IDs/intensity, and meaningful invalid-parameter exceptions. Production algorithms did not change.

Other class-test deltas are null-pointer spelling, explicit initialization, logical `&&` instead of boolean `&`, literal escape cleanup, avoiding needless string temporaries, and removal of unused local helpers/data. AASequence's four edits only replace `0` with `nullptr`. MzMLFile's edit changes the local accession loop variable to `const char*`. No existing scientific fixture values were changed.

## Build/delivery scope changes

* Root CMake minimum actually rises from 3.21 to 3.24 (AGENTS already recommended 3.24); C++23 remains. Presets are reduced to Core Debug/Release, using native CMake unity batching.
* Core removes unused Boost date_time/iostreams link dependencies; regex and Boost headers remain. Removes stale architecture-probing and installer/build wrappers. New dependency-default module prefers ordinary libraries over Apple frameworks unless explicitly configured otherwise; installed config includes it.
* Tests share `build/bin` beside Core DLLs; Windows relies on its active dependency environment instead of copying all dependency DLLs. `/utf-8` fixes MSVC literal encoding. This is distinct from Rust's Git checkout CRLF fixture issue.
* Adds five-platform Release CI (Linux x64/arm64, macOS x64/arm64, Windows x64), conda-forge dependency environment, tested SDK archives/release validation, and per-platform concurrency without canceling active builds. These are C++ SDK delivery changes; no Rust dependency addition is required by them.
* Removes 69 obsolete orchestration/packaging assets, predominantly `cmake/MacOSX`, `cmake/Windows`, `cmake/knime`, Jenkins/Gitpod and old build/dependency scripts. **No scientific class or runtime data is removed.** `source-provenance.json` now explicitly says its imported-file list is historical, not a live inventory.
* Latest AGENTS adds `core-release` delivery and warns that Windows conda Release dependencies must not be linked into an ordinary MSVC Debug build. Existing no-third-party-edits/build-authorization rules remain.

## Required Rust follow-up

1. Update the SDK target/version and compatibility report deliberately, preserving every historical source/fixture manifest's original revision/hash. Changed source files such as IsoelectricPoint.cpp, StringUtils.cpp, AhoCorasickAmbiguous.h and several class tests must be recorded as **changed but reviewed for the implemented native behavior**, not falsely byte-identical.
2. No immediate Rust scientific production patch is justified by this delta. Rerun the existing affected peptide-properties, indexing, ParamValue/list/text, and comparison references as part of the normal Rust validation after a pin update; no new C++ build is needed to establish this source comparison.
3. Use the strengthened tests/error semantics above when implementing currently unmapped normalization, Arrow/Parquet exporters, ZIP archive editing, process probes, and ion-mobility APIs. Do not upgrade their completion status merely because upstream code changed.
4. Keep earlier ZIP-input and current runtime-batch conclusions: `ZipIfstream`, SYSTEM/File, LogStream, ProgressLogger, DefaultParamHandler, ParamXML, FeatureXML, ConsensusXML and ModificationDefinitionIO implementations/headers are unchanged in this upstream range. Their existing native differences and outstanding gaps remain.

The reference checkout, Rust files, dependencies, ledger and fixtures were not modified. No C++ or Rust compilation or tests ran during this audit.

## Applied source-target update

The Rust target inventory now records `54a232f`. Historical fixture revisions and hashes remain unchanged. Each earlier current-SDK manifest records a `target_verification` with the reviewed revision, original revision and any changed source hashes. Seven of the 220 historical reference paths changed and carry explicit original/target hashes plus behavior reviews; the remaining 213 are byte-identical. The mzML class-test source changed only its accession loop variable type and has a corresponding manifest review. `tools/check_core_sdk.py` checks this mapping and optionally verifies every target hash against the clean independent checkout `.reference/openms4-core-54a232f`. These checks establish traceability, not executed C++ parity.
