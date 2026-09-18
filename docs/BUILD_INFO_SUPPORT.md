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
| Linux | `VERSION_ID` from `/etc/os-release` | same |
| other Unix (`__unix__`) | `VERSION_ID` from `/etc/os-release` | `"unknown"` |
| macOS | `sysctlbyname("kern.osproductversion")` | `"unknown"` |
| Windows | `RtlGetVersion` from `ntdll.dll` | `"unknown"` |

`/etc/os-release` is read under `cfg(target_os = "linux")`, not under
`cfg(unix)`, so the second row is narrower than the source's `__unix__` branch:
a FreeBSD build that publishes the file still reports `"unknown"` here. The OS
*name* is unaffected — the port maps every non-macOS Unix to `Linux`, exactly as
the source does — so only the version string differs, and only on a platform
this crate does not currently target.

The two missing platform probes need a binding this crate does not take. The
value is a display string and nothing in the port branches on it, so an honest
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

**OpenMP is absent, and says so.** No OpenMP runtime is linked here: the 36
files under `src/openms/source` that carry `#pragma omp` have no counterpart,
and nothing in the crate is threaded by OpenMP. So
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

## The required-CPU-feature check (native, no source counterpart)

`.cargo/config.toml` builds every **x86** target of this repository with
`-C target-feature=+fma`. The measurement behind that decision is section 4 of
[BENCHMARKS](BENCHMARKS.md): the completed `FeatureFinderAlgorithmPicked` costs
21 % at one thread on a baseline x86-64 build, because every `f64::mul_add` of
the ported glibc `powf`, `exp` and `log` becomes two indirect calls into
`compiler_builtins`' `fma` stub instead of one instruction; the flag removes
that cost and more (0.715 at one thread, 1.055 against the C++ build), and the
output is **bitwise identical** with and without it at both measured thread
counts. AArch64 is untouched: `fma` is an x86 target-feature name, and AArch64
has fused multiply-add in its base instruction set.

rustc implies `avx`, `sse3`, `ssse3`, `sse4.1` and `sse4.2` from `fma`, so the
whole crate may be compiled with VEX encoding, and such a binary needs Intel
Haswell (2013) or AMD Piledriver (2012) or later. Three items answer that:

| Item | What it does |
|---|---|
| `RequiredCpuFeatures` | `None` (not x86, or built without the flag), `FmaPresent`, `FmaMissing` |
| `required_cpu_features()` | `cfg!(target_feature = "fma")` for the build, `std::arch::is_x86_feature_detected!("fma")` for the processor. Both are safe; the crate keeps `#![forbid(unsafe_code)]` |
| `check_required_cpu_features(err)` | `true` and nothing printed when the binary can run; otherwise `FMA_MISSING_MESSAGE` on `err` and `false`. `cli::run_with` calls it first and returns `ExitCode::InternalError` (12) on `false`, so every TOPP executable of this repository refuses before it reads its arguments |

`FMA_MISSING_MESSAGE` is two lines: what is missing and which processors have
it, then how to build a binary that runs here
(`RUSTFLAGS='' cargo build --release`, which replaces the section of
`.cargo/config.toml` wholesale).

**It is a courtesy, not a guarantee, and the documentation says so.** The whole
crate is compiled with those instructions, so a processor without them can fault
on one before the check is reached; the check is the first thing each tool does,
but it cannot be the first instruction of the process. What it guarantees is
that a machine which gets that far is told the cause and the remedy instead of
being left with `SIGILL`. A narrower alternative — `#[target_feature(enable =
"fma")]` on the three ported replica functions with runtime dispatch — needs
`unsafe fn` at this MSRV and was not taken; BENCHMARKS section 4.4 records it as
unmeasured.

`active_simd_extensions` and `required_cpu_features` read the same build flag
from the two sides, and `tests/build_info.rs` asserts that they agree: an x86
binary lists `FMA` exactly when it requires it.

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
