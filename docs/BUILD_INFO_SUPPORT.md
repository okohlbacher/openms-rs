# Build and runtime platform identity

`system::build_info` covers `SYSTEM/BuildInfo.h` and `SYSTEM/BuildInfo.cpp` at
SDK `bc9cc12`. Its one direct TOPP consumer is `OpenMSInfo`.

The header answers two questions and keeps them apart: what platform is
*running* (`OpenMSOSInfo`), and what the compiler *emitted* (`getBinaryArchitecture`,
`getActiveSIMDExtensions`, `OpenMSBuildInfo`). Much of the second half has no
meaningful Rust counterpart, and this document says which and why.

## API mapping

Every public member of the header and its `.cpp`, with its Rust counterpart.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `enum class OpenMS_OS` | `build_info::OperatingSystem` | `OS_UNKNOWN` → `Unknown`, `OS_MACOS` → `MacOs`, `OS_WINDOWS` → `Windows`, `OS_LINUX` → `Linux`. |
| `OpenMS_OS::SIZE_OF_OPENMS_OS` | — | Not ported: a sentinel that exists only to size the parallel name array. |
| `OpenMS_OSNames[]` | `OperatingSystem::as_str` | The array and the index into it become one method. |
| `enum class OpenMS_Architecture` | `build_info::Architecture` | `ARCH_UNKNOWN` → `Unknown`, `ARCH_32BIT` → `Bits32`, `ARCH_64BIT` → `Bits64`. |
| `OpenMS_Architecture::SIZE_OF_OPENMS_ARCHITECTURE` | — | Not ported: array-sizing sentinel. |
| `OpenMS_ArchNames[]` | `Architecture::as_str` | |
| `std::string getOSVersionString_()` | `build_info::os_version_string() -> String` | Public here; Linux only, see *Native differences*. |
| `class OpenMSOSInfo` | `build_info::OsInfo` | |
| `OpenMSOSInfo::OpenMSOSInfo()` | `OsInfo::default` | Everything unknown. |
| `~OpenMSOSInfo()` | — | Implicit on both sides. |
| `std::string getOSAsString() const` | `OsInfo::os` + `OperatingSystem::as_str` | Public field plus label method. |
| `std::string getArchAsString() const` | `OsInfo::arch` + `Architecture::as_str` | |
| `std::string getOSVersionAsString() const` | `OsInfo::os_version` | Public `String` field. |
| `static std::string getBinaryArchitecture()` | `build_info::binary_architecture() -> Architecture` | From `size_of::<usize>()`, as the source uses `sizeof(size_t)`. |
| `static std::string getActiveSIMDExtensions()` | `build_info::active_simd_extensions() -> String` | Rust target features replace the SIMDe probes; same order, same labels. |
| `static OpenMSOSInfo getOSInfo()` | `OsInfo::detect() -> OsInfo` | |
| `OpenMSOSInfo::os_`, `os_version_`, `arch_` (private) | `OsInfo::os`, `os_version`, `arch` (public) | See *Native differences*. |
| `struct OpenMSBuildInfo` | the module itself | A stateless namespace of static methods becomes free functions. |
| `static bool isOpenMPEnabled()` | `build_info::openmp_enabled() -> bool` | Always `false`. |
| `static std::string getBuildType()` | `build_info::build_type() -> &'static str` | `"Debug"` / `"Release"` from `debug_assertions`. |
| `static Size getOpenMPMaxNumThreads()` | `build_info::openmp_max_num_threads() -> usize` | Always `1`. |
| `static void setOpenMPNumThreads(Int)` | — | **Not ported**: a no-op in the source's own non-OpenMP build, and there is no thread pool here to size. |

Native additions: `build_info::UNKNOWN`.

## Preserved source conventions

**The name tables are exact**: `"unknown"`, `"MacOS"`, `"Windows"`, `"Linux"`,
and `"unknown"`, `"32 bit"`, `"64 bit"` — including the capitalisation of
`MacOS` and the space in `32 bit`.

**A default instance is entirely unknown**, including the version string, which
is the literal `"unknown"` rather than an empty string. The class test asserts
all three.

**The OS probe order is the source's**: Windows, then macOS, then any other
Unix, which is therefore reported as `Linux`. The source's own `TODO` asks
whether FreeBSD should be distinguished; it is not, and neither is it here.

**The architecture of a probe is never unknown.** `getOSInfo` maps four-byte
pointers to 32-bit and *everything else* to 64-bit, so unlike
`getBinaryArchitecture` — which has a genuine `unknown` case — the probed field
cannot report it. That asymmetry is reproduced rather than tidied away.

**`VERSION_ID` parsing is the source's**: a prefix match on `VERSION_ID=`, then
one surrounding pair of double quotes removed, and only when the value is at
least two characters long, so a lone `"` is returned unchanged.

**The SIMD list keeps the source's order and labels** — `neon`, `SSE`, `SSE2`,
`SSE3`, `SSE4.1`, `SSE4.2`, `AVX`, `AVX2`, `FMA` — joined with `", "`, and is
empty when none is active.

## Native differences

**Fields instead of accessors.** The source keeps `os_`, `os_version_` and
`arch_` private behind three `*AsString` accessors. `OsInfo` exposes them,
because the enums carry their own labels and hiding the values would only stop a
caller from matching on them. The label methods remain, so the accessors have an
exact counterpart.

**The OS version is Linux-only.**

| Platform | Source | Port |
|---|---|---|
| Linux / other Unix | `VERSION_ID` from `/etc/os-release` | same |
| macOS | `sysctlbyname("kern.osproductversion")` | `"unknown"` |
| Windows | `RtlGetVersion` from `ntdll.dll` | `"unknown"` |

Both missing probes need a platform binding this crate does not take. The value
is a display string and nothing in the port branches on it, so an honest
`"unknown"` is preferable to a guess — and it is exactly what the source itself
returns when its own probe fails. The class test only requires a non-empty
string, which the fallback satisfies. Adding them is a deferral, not a gap in
the mapping.

**SIMD reporting is re-based, not re-implemented.** The source reads the SIMDe
`SIMDE_ARCH_*` macros, which record what was defined when the C++ was compiled.
The port reads the corresponding Rust target features, which the compiler
enables from the target definition and from `-C target-feature` / `-C
target-cpu`. Both answer the same question — what the compiler was allowed to
emit — and both are fixed for a given binary, which is the only property the
class test asserts. Neither is a claim that the code *uses* those instructions;
this port has no hand-written SIMD at all.

**OpenMP is absent, and says so.** The port is deliberately serial: the source's
36 files with `#pragma omp` have no threaded counterpart here. So
`openmp_enabled()` is `false` and `openmp_max_num_threads()` is `1` — which are
not approximations but *exactly* what the source returns when `_OPENMP` is not
defined, so this is the one part of the header the port reproduces perfectly by
having less. `OMP_NUM_THREADS` has nothing to act on, and `setOpenMPNumThreads`
is dropped rather than shipped as a no-op nobody can observe.

**`build_type` is the nearest observable fact.** The source returns CMake's
`OPENMS_BUILD_TYPE` string verbatim. Cargo bakes no profile name into a library,
so the port reports `"Debug"` when debug assertions are on and `"Release"`
otherwise. That separates Cargo's stock `dev` and `release` profiles, but it
follows the `debug-assertions` *setting*: a custom profile that turns assertions
on in an optimised build reports `"Debug"`.

**`/etc/os-release` reads are capped** at 64 KiB; the source reads the file line
by line with no bound.

## Checked boundaries and evidence

| Boundary | Value | Source behaviour |
|---|---|---|
| `/etc/os-release` read | ≤ 64 KiB | unbounded |

All eight class-test sections are mapped in `tests/build_info.rs`. The three
`"unknown"` defaults and the seven name-table labels are transcribed (tier 3).
The binary-architecture and OS-probe sections are derived exactly as the source
test derives its own expectations — from `sizeof(size_t)` and from the accepted
set `{Windows, MacOS, Linux}` — rather than pinned to whichever host runs them
(tier 4). The SIMD section asserts the determinism the source test asserts, plus
that every emitted label is one of the nine the source can emit and that the
join produces no empty element.

`OpenMS/build_config.h` also carries the library version, which is not part of
this header's surface; the crate's counterpart is `openms::CORE_SDK_VERSION` in
`src/lib.rs`.

No C++ execution is claimed. Source hashes, line anchors and the class-test
review are in [the provenance record](../tests/data/build_info_provenance.json).
