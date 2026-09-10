# Historical OpenMS4 core: Rust port analysis

**Current target:** the [reduced SDK update](CORE_SDK_UPDATE.md) supersedes the revision and scope counts below. This original analysis and its inventory remain pinned to `7c029e8` as historical evidence.

This repository starts a **native Rust scientific library**, with spectra, chemistry, and common processing as the first implementation scope. The C++ repository is much larger than this starting scope: its first-party library headers and implementations contain **484,760 physical lines**, and its class-test sources contain another **235,900 lines**. A complete replacement requires incremental scientific validation, format compatibility work, and decisions about native dependencies. This document is an inventory and engineering plan, not a claim that all C++ functionality has been ported.

## Source identity and scope

The analyzed source is [okohlbacher/OpenMS4-core at commit `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`](https://github.com/okohlbacher/OpenMS4-core/tree/7c029e8cdba6abab503708ecdd56f6ab55e38ce4), inspected from an existing clean checkout on 2026-09-10. Its recorded upstream extraction point is `ca32296038839459d8c9b075b759e285913d6294`, and its core SDK version is 4.0.0. These revisions identify different things: the latter is upstream provenance; the former identifies the actual extracted package being ported.

The extraction builds OpenMS and OpenSwathAlgo, optional test support, and scientific class tests. It deliberately excludes GUI products, TOPP executables, Python bindings, and generated documentation from the build. `ConsoleUtils` is the one remaining `APPLICATIONS` class. The source describes itself as an experimental extraction candidate, with C++ binary acceptance still outstanding. No C++ configuration, compilation, or C++ test execution was performed for this analysis. See the pinned [extraction notes](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/CORE_IMPLEMENTATION_NOTES.md) and [source target list](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/CMakeLists.txt).

## Measured inventory

[source-inventory.json](source-inventory.json) records every counted path, category, domain, byte count, line counts, and SHA-256. The file set comes from `git ls-files` at the pinned revision. Public `.h` files, private `.h` files under the implementation tree, `.cpp` implementations, and class-test `.cpp` sources are counted separately. The three generated `.h.in` templates add 235 lines and are excluded from the runtime line total below.

“Lines” means physical lines, including comments and blanks. The JSON additionally records nonblank lines; neither metric removes comments or purports to measure executable statements. File counts are not counts of classes, test cases, or enabled build targets. Conditional implementations and auxiliary class-test programs remain in the inventory.

All 2,341 recorded source/template/test hashes were verified against immutable Git objects at the pinned revision. Canonical byte counts and hashes use Git blob content; the one CRLF checkout conversion (`FORMAT/MSNUMPRESS/MSNumpress.cpp`) is recorded separately. The clean-worktree observation refers to the initial inventory, so later work in the source checkout does not change this pinned analysis.

| Domain | Public headers | Internal headers | Implementations | Runtime lines | Test sources | Test lines |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ANALYSIS | 214 | 11 | 212 | 145,313 | 190 | 65,936 |
| APPLICATIONS | 1 | 0 | 1 | 297 | 1 | 163 |
| CHEMISTRY | 66 | 1 | 59 | 38,931 | 58 | 22,829 |
| COMPARISON | 13 | 0 | 13 | 2,899 | 13 | 1,523 |
| CONCEPT | 21 | 0 | 16 | 7,127 | 14 | 2,610 |
| DATASTRUCTURES | 44 | 0 | 36 | 23,415 | 38 | 15,204 |
| DEPENDENCY_TESTS | 0 | 0 | 0 | 0 | 3 | 723 |
| FEATUREFINDER | 46 | 0 | 43 | 24,826 | 42 | 6,697 |
| FORMAT | 186 | 1 | 184 | 143,998 | 153 | 55,135 |
| IMAGING | 5 | 0 | 4 | 1,283 | 4 | 874 |
| INTERFACES | 3 | 0 | 0 | 425 | 0 | 0 |
| INTERFACES_IMPL | 0 | 0 | 1 | 131 | 0 | 0 |
| IONMOBILITY | 4 | 0 | 5 | 1,177 | 4 | 686 |
| KERNEL | 34 | 0 | 32 | 20,323 | 34 | 18,823 |
| MATH | 17 | 0 | 16 | 9,012 | 17 | 5,168 |
| METADATA | 66 | 1 | 48 | 24,336 | 49 | 16,445 |
| ML | 32 | 1 | 31 | 10,587 | 27 | 5,631 |
| OPENSWATHALGO | 13 | 0 | 9 | 2,540 | 7 | 1,092 |
| PROCESSING | 31 | 0 | 31 | 17,475 | 31 | 9,255 |
| QC | 19 | 0 | 19 | 5,700 | 18 | 3,449 |
| SYSTEM | 14 | 0 | 14 | 4,965 | 13 | 1,935 |
| TEST_SUPPORT | 0 | 0 | 0 | 0 | 4 | 1,722 |
| **Total** | **829** | **15** | **774** | **484,760** | **720** | **235,900** |

Most test domains are assigned from a uniquely matching public-header basename. Nonmatching cases use inspected includes, with explicit reviewed overrides for helper programs, dependency tests, and classes whose filenames differ. Assignment reasons are recorded in the JSON. Seven sources test the separate OpenSwathAlgo domain, although only five reside in its own suite directory. Domain counts express source ownership rather than coverage percentages.

The tracked tree has 5,823 paths. The inventory excludes 2,213 files under `src/openms/extern/` and 84 under `src/openms/thirdparty/`; those are vendored dependencies, not a Rust translation backlog. There are also 579 files under class-test data directories and 185 under `share/OpenMS/`. Test-framework library code, support headers, build scripts, and the separate Python SDK-contract test file are outside the table. All source and vendored dependency trees were left unchanged.

`ANALYSIS` and `FORMAT` account for about 60% of the measured runtime lines. The difficult part of the port is therefore not only peak storage or mass arithmetic: it includes identification, quantification, OpenSWATH workflows, specialized formats, and interoperability with scientific infrastructure.

## Actual C++ build and dependencies

The pinned root configuration requires CMake 3.21 and sets the SDK to 4.0.0. Contrary to the generic repository agent notes' C++20 description, the actual target configuration exports **C++23** through `cxx_std_23`. The source files are authoritative for the analyzed revision. See [root configuration](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/CMakeLists.txt), [compiler flags](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/cmake/compiler_flags.cmake#L72), and [library target helper](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/cmake/add_library_macros.cmake#L150).

| Dependency area | Observed C++ configuration | Rust direction |
| --- | --- | --- |
| Containers, strings, errors | STL, OpenMS helpers, Boost; Boost 1.81 minimum | Rust collections, owned strings, slices, iterators, `Result`, typed configuration |
| Linear algebra | Eigen 3.4 minimum, public target | Choose Rust numerical backends per algorithm after parity benchmarks; avoid a mandatory backend for basic spectra |
| XML and compressed files | XercesC, zlib, bzip2, libzip | Optional Rust format modules; parse XML as events and decode binary arrays explicitly |
| Columnar data | Arrow and Parquet >=23 required; dataset component probed | Optional native Rust Arrow/Parquet support; preserve on-disk schemas and metadata |
| Database and transport | SQLiteCpp, libcurl | Separate adapters from the in-memory scientific core |
| Optimization and learning | LIBSVM >=2.91; COIN-OR, GLPK, or HiGHS; bundled Percolator | Explicit optional solver/inference interfaces; do not substitute numerical solvers without validating outcomes |
| Scientific helpers | IsoSpec, eol-bspline, Evergreen, GTE, Quadtree, SIMDe | Port only owned algorithms where justified; isolate any retained native dependency behind a small optional boundary |
| Parallelism | OpenMP enabled by default | Serial reference behavior first; optional Rust parallel iterators after deterministic tests |
| Optional platform features | HDF5 and ONNX off; opentims on; Thermo bridge generally on, default off on Linux aarch64 | Feature-gated adapters with documented platform and runtime requirements |

The dependency list is derived from the pinned [discovery module](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/cmake/cmake_findExternalLibs.cmake) and [OpenMS link dependencies](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/CMakeLists.txt). The LP solver selection tries COIN-OR, then GLPK, then installed or fetched HiGHS. Qt discovery is conditional on GUI, which the extracted root disables. CMake can also fetch vendor bridges and models, so configuration itself is not equivalent to a harmless inventory command.

The Rust default should remain small: scientific containers, chemistry, and common transformations do not need the complete C++ dependency closure. For XML, `quick-xml` provides streaming reader/writer APIs; this is suitable infrastructure for mzML parsing but does not implement its scientific semantics. Native Rust Arrow provides columnar types; its version numbers should not be equated with Arrow C++ version numbers. Rayon offers optional parallel iteration. These are dependency choices for the Rust architecture, not claims of OpenMS compatibility. Current primary documentation: [quick-xml](https://docs.rs/quick-xml/latest/quick_xml/), [Arrow Rust](https://docs.rs/arrow/latest/arrow/), [Rayon](https://docs.rs/rayon/latest/rayon/). Freeze tested dependency versions in the project lockfile, with a documented Rust toolchain requirement.

## Native Rust architecture

Use one coherent crate initially, organized by domain, and split crates only when dependency boundaries or compile times warrant it. The supported surface should be enumerated in the README and executable examples; an entire C++ class name is too coarse a unit for a compatibility claim.

| Rust module | Responsibility and boundary |
| --- | --- |
| `kernel` | Peaks, spectra, chromatograms, experiments, ranges, sorted search, aligned auxiliary arrays |
| `metadata` | Precursors, scan properties, typed values, and extensible annotations; composition replaces inheritance |
| `chemistry` | Elements, formulas, residues, peptides, modifications, digestion, isotope and fragment calculations |
| `processing` | Stateless or explicitly configured transformations over spectra and experiments |
| `format` | Readers/writers over `BufRead` and `Write`, including bounded optional mzML/idXML subsets; stream ownership and unsupported data are explicit |
| `identification` | Owned peptide/protein search-result records, evidence coverage and modification mapping; search/inference algorithms are explicit operations in `analysis` |
| `analysis` | Retention-time models/application, identification score switching, HyperScore/Morpheus, peptide/protein filtering and target/decoy FDR, peptide indexing, basic protein inference/grouping, conflict resolution and origin partitions; landmark discovery, database search and probabilistic inference remain separate work |
| Later `math`, `qc` modules | Scientific algorithms added with behavior fixtures and their own explicit configuration |
| Optional adapters | Arrow/Parquet, vendor readers, databases, inference runtimes, solver backends |

Important design commitments:

- Preserve scientific units explicitly: m/z and retention time use `f64`; the source's `Peak1D` intensity uses `float`, so `f32` is the natural compatibility baseline. Accumulation may require `f64` and tolerance-based comparison. Keep retention time in seconds at the in-memory boundary.
- Store owned values and expose borrowed slices. Sorting, filtering, and selecting peaks must apply the same permutation to every associated data array. Invalid array lengths should be errors rather than partial mutations.
- Expose sorted-search preconditions clearly. Either validate order when searching or represent sorted views; silently assuming order produces plausible but wrong scientific results.
- Replace exception hierarchies and debug-only preconditions with documented fallible operations. Unsupported formats, unknown residue/modification syntax, and unimplemented compression methods must return explicit errors.
- Keep algorithms independent of filesystem and global environment state. Share immutable chemistry tables safely and record their provenance; do not reconstruct C++ global registries or data discovery paths by accident.
- Prefer typed parameter structs and enums for new Rust APIs. Preserve a mapping from C++ parameter names/defaults where comparison or imported workflows require it.
- Treat Rust API ergonomics and scientific result compatibility separately. Native methods need not mimic C++ inheritance, mutable output parameters, exported symbols, or ABI. Interoperability can be added later without making the core a C++ wrapper.

## Behavior parity and available reference cases

The initial oracle is the **pinned C++ test source and committed fixtures**. A Rust test copied from a known C++ assertion provides source-derived reference coverage. It does not prove that the pinned C++ checkout compiles or that the implementations have been run side by side. Record that distinction in release notes.

Useful reference cases available without building C++:

| Area | Existing source-derived expectation | Reference |
| --- | --- | --- |
| Threshold filtering | `Transformers_tests.dta` has 121 peaks; threshold 1 retains 121, threshold 10 retains 14 | [ThresholdMower test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/ThresholdMower_test.cpp) |
| Normalization | The same fixture has maximum intensity 46; maximum normalization yields 1 and TIC normalization sums to approximately 1 | [Normalizer test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/Normalizer_test.cpp) |
| Top-N filtering | A 100-peak triangle retains descending-intensity indices `50,51,49,52,48,53,47,54,46,55`; integer and string arrays follow these indices | [NLargest test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/NLargest_test.cpp) |
| Formula mass and charge | Formula order is immaterial; empty formula mass is zero; source charge semantics add/subtract proton mass | [EmpiricalFormula test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/EmpiricalFormula_test.cpp) |
| Trypsin digestion | `ACKDE` and `ACRDE` split after K/R; `ACKPDE` and `ACRPDE` exercise the proline exception; missed-cleavage tests specify order | [ProteaseDigestion test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/ProteaseDigestion_test.cpp) |
| DTA loading | Fixture assertions include precursor m/z approximately 582.40666 and the first peak `(139.42, 318.52)` | [DTAFile test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/DTAFile_test.cpp) |
| FASTA handling | Tests cover descriptions, multiline sequences, whitespace, PEFF headers, empty lines, and malformed symbols | [FASTAFile test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/FASTAFile_test.cpp) |

Each adopted fixture should retain its original license, pinned source path, SHA-256, expected results, and numerical tolerance. Prefer full peak arrays and metadata assertions over only output lengths. Add boundary cases that test behavior, including empty inputs, equal intensities, duplicate positions, malformed input, unsupported syntax, nonfinite numbers, and inconsistent data arrays. For chemistry, test mass/formula consistency and digestion cleavage boundaries independently of parser internals.

Subsequent differential validation should use an explicitly authorized, independently built C++ reference executable at the pinned revision. Feed the same compact cases into both implementations and compare a canonical, versioned representation of results. Exact comparison is appropriate for identities, sequence order, integer counts, and flags; floating-point outputs need documented absolute/relative tolerances and a deliberate nonfinite-value policy. Do not require byte-identical compressed files: the C++ dependency module already notes zlib versus zlib-ng differences. Compare decoded values, metadata, and format semantics instead.

## Portability and scientific risks

C++ APIs contain templates, iterator-heavy algorithms, inheritance, shared ownership, and configuration cached after mutation. Rust requires explicit ownership and borrowing decisions; simple textual translation will not preserve lifetime or mutation semantics. In particular, metadata alignment and references between identification objects need stable identifiers or explicit ownership rather than borrowed pointers stored across mutations.

Numerical differences can arise from `float` storage, accumulation order, fused operations, solver tolerances, random generators, and parallel reductions. Deterministic serial reference tests should precede parallel optimization. A faster result is useful only after the scientific equivalence criteria are stated and met.

File readers carry substantial semantics. mzML requires controlled vocabulary interpretation, unit conversion, reference groups, binary precision, compression, spectrum/chromatogram metadata, and potentially indexing. MGF and FASTA also have dialects and edge cases. An initial reader for common records must list its supported subset; merely parsing an example or writing a self-readable file does not establish full standards compliance or lossless OpenMS round trips.

Chemistry is data as well as code. Element masses, isotope abundances, residue definitions, modification databases, and enzyme rules must have pinned provenance. Supporting common unmodified peptides or a handful of modifications is a useful milestone but leaves a substantial part of `AASequence`, `ProForma`, nucleic-acid chemistry, and isotope modeling outstanding.

Vendor raw readers, HDF5, ONNX, and numerical solvers can reintroduce native build and distribution constraints. Evaluate them separately from the portable Rust core, including Linux aarch64 and macOS arm64. Keep optional native adapters visibly distinct from native Rust implementations in the capability matrix.

## Staged roadmap and acceptance gates

1. **Foundation and a working scientific slice.** Provide the native kernel, chemistry foundations, common processing, text I/O, executable examples, and source-derived reference tests. Gate: the supported surface is documented, tests run in Rust, and unsupported behavior is explicit. This is the intended first delivery; it does not complete the whole-library port.
2. **Interchange and chemistry depth.** Extend mzML with a documented coverage matrix, richer metadata, indexed/streaming access where needed, complete modification handling, and broader digestion/fragment/isotope semantics. Gate: independent input fixtures, error-path tests, validated output interoperability, and data provenance.
3. **Numerical processing and feature workflows.** Port smoothing, baseline correction, centroiding, alignment, feature finding, and related quantitative operations in dependency order. Gate: per-algorithm golden outputs, tolerances, memory profiles, and representative performance baselines.
4. **Identification, quantification, and OpenSWATH.** Port inference, scoring, targeted workflows, database/columnar interchange, and solver-dependent algorithms. Gate: workflow-level comparisons and scientifically meaningful metrics, not only unit-test counts.
5. **Release acceptance across platforms.** Validate Linux, macOS, and Windows; test optional feature combinations; establish published API stability, fixture provenance, dependency/license notices, and explicit compatibility claims. Gate: a reviewed method-level coverage matrix and differential validation against the authorized C++ reference.

The work should be measured by validated capabilities and workflow coverage. The source inventory gives a reproducible denominator for planning; it is not a meaningful percentage-complete calculation for a native Rust redesign.
