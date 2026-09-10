# Current reduced Core SDK target

The Rust port now targets OpenMS Core SDK **4.0.0** at
[`6bfc0e4711105f4eda2fea86812a83af7c7e791f`](https://github.com/okohlbacher/OpenMS4-core/tree/6bfc0e4711105f4eda2fea86812a83af7c7e791f),
fetched from the repository's default `codex/package-split` branch on 2026-09-10.
The original reference archive at `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`
is retained for the existing scientific fixtures. The package's upstream
extraction revision `ca32296038839459d8c9b075b759e285913d6294` identifies a
separate historical input, not the current SDK package.

The [inventory and comparison](core-sdk-update.json) records every current
scientific source path, hash, size and build-registration role, all changed or
removed paths, and **220 unchanged historical reference paths**. It distinguishes
current target identity from the origins of existing source-derived tests.
`openms::CORE_SDK_VERSION` and `openms::CORE_SDK_REVISION` expose this target in Rust;
they do not assert a C++ ABI, a Rust artifact hash or complete API parity.

## What the reduction changes

The current SDK removes product-specific backends and retains a broad scientific
library. The comparable physical inventory is:

| Source category | Earlier reference | Current SDK | Removed |
| --- | ---: | ---: | ---: |
| Core include-directory `.h` files | 816 | 794 | 22 |
| Core `.cpp` implementations | 765 | 746 | 19 |
| Core source-side private headers | 15 | 15 | 0 |
| OpenSwathAlgo include headers | 13 | 13 | 0 |
| OpenSwathAlgo `.cpp` files | 9 | 9 | 0 |
| Core class-test `.cpp` files | 715 | 703 | 12 |

Including private headers and OpenSwathAlgo, runtime physical lines decrease
from **484,760 to 468,372**. Of 1,618 old/current runtime paths, 1,569 are
unchanged, 41 removed and eight changed. File and line counts are not method
coverage percentages.

Current reachable build registrations contain 774 Core and twelve OpenSwathAlgo
public headers, 746 Core and eight OpenSwathAlgo implementations, plus twenty
private headers. This is a static union of optional branches, not the result of
configuring a portable native profile. Generated export/configuration headers
are additional. Physical include-directory files are therefore not an exact
public API count. Six ONNX/PeptDeep headers have a declaration list that is not
reached from current `includes.cmake`, although optional implementations are
attached separately; the inventory records that distinction.

These classes have left Core and are excluded from its remaining Rust backlog:

- `ProSEAlgorithm`.
- `NuXLAnnotateAndLocate`, `NuXLAnnotatedHit`, `NuXLConstants`, `NuXLDeisotoper`,
  `NuXLFDR`, `NuXLFragmentAdductDefinition`, `NuXLFragmentAnnotationHelper`,
  `NuXLFragmentIonGenerator`, `NuXLModificationsGenerator`, `NuXLParameterParsing`
  and `NuXLPresets`.
- `FLASHDeconvAlgorithm`, `MassFeatureTrace`, `Qvalue`, `SpectralDeconvolution`
  and `TopDownIsobaricQuantification`.
- `DDAWorkflowCommons`, `MascotRemoteQuery`, `ParquetTableComparator`,
  `CometNativeIDRemapper` and `DBSuitability`.

None had a corresponding implementation in this Rust port, so updating scope
does not require deleting working Rust APIs. Core still owns `NuXLReport`,
`NuXLMarkerIonExtractor`, `PeakGroup`, `PeakGroupScoring`, `DeconvolvedSpectrum`,
FLASH records/file writers, Comet modification records and reusable readers and
writers. The current [SDK ownership notes](https://github.com/okohlbacher/OpenMS4-core/blob/6bfc0e4711105f4eda2fea86812a83af7c7e791f/README.md#tool-backend-extraction)
are authoritative for this split.

## Changes relevant to the port

All chemistry, kernel, processing, comparison, mathematical and data-structure
scientific sources used by the existing Rust implementation are unchanged.
The historical reference comparison also verifies the existing test inputs,
source assertions, CV/mapping files and XML schemas. No mass-table regeneration,
fixture replacement or numerical tolerance change is justified by this update.

The retained changes are:

| Source | Consequence |
| --- | --- |
| `PeakGroup` and `PeakGroupScoring` | Isotope-cosine helpers move from the extracted FLASH backend into retained scoring. Fractional isotope residuals now use floating-point absolute values; the new source regression must guide their eventual Rust implementation. |
| `VersionInfo` | Full package revision, source dirtiness and native build identity are exposed. The Rust target constants identify the reference SDK separately from the Rust package/build. |
| `SYSTEM/File` | Installed data discovery and error/retry behavior change. Embedded Rust chemistry resources do not need C++ runtime path discovery. The broader filesystem API is still unported. |
| `ParamCTDFile` | An explicit `StringUtils` include repairs dependency closure; scientific serialization behavior is unchanged. |

The identification graph remains registered and byte-identical. Its
[sequence and provenance layer](IDENTIFICATION_GRAPH_SUPPORT.md) now supports
stable owned references, registration and score histories, translated merge/copy,
parent coverage and atomic RNase integration. Observations, compounds, adducts and
observation matches are also implemented with typed molecule dispatch and source
merge semantics. Match groups, parent groups, full cleanup/persistence and the
legacy converter remain separate requirements.

Other retained gaps include ProForma/MzPAF, mass decomposition and XLMS chemistry,
additional formats and Arrow/Parquet, native readers, feature finding/grouping,
calibration, further alignment and numerical models, targeted/quantitative and
OpenSWATH workflows, ML/QC and remaining search/inference/scoring methods. Existing
partial formats/models retain the limits in [Porting status](PORTING_STATUS.md).
The historical [repository analysis](REPOSITORY_ANALYSIS.md) is retained as the
original baseline; this document supersedes its target revision and scope counts.

## Reproduce the source-target checks

The packaged check verifies target constants, inventory totals, scope changes
and existing-manifest links without a C++ checkout:

```bash
python3 tools/check_core_sdk.py
```

With a clean independent checkout of the pinned SDK, it also verifies the exact
revision, full scientific file set and 1,788 distinct source/registration/reference
hashes:

```bash
python3 tools/check_core_sdk.py --source .reference/openms4-core-current
```

The old archive is not a Git checkout. The check explicitly verifies the supplied
Git root to avoid accidentally querying an enclosing repository. Source inspection
and these hash checks do not establish a native C++ build or differential numerical
parity. Rust build and test results are recorded separately in [Validation](VALIDATION.md).
