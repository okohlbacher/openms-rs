// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Build and runtime platform identity of the Core SDK `SYSTEM/BuildInfo.h`.
//!
//! The header answers two different questions and keeps them apart.
//! [`OsInfo`](crate::system::build_info::OsInfo) describes the *running*
//! platform — operating system, its version, and the pointer width the binary
//! was built for — while
//! [`binary_architecture`](crate::system::build_info::binary_architecture),
//! [`active_simd_extensions`](crate::system::build_info::active_simd_extensions)
//! and [`build_type`](crate::system::build_info::build_type) describe what the
//! compiler emitted.
//!
//! The source's OpenMP accessors survive as honest constants: this port is
//! serial by design, so
//! [`openmp_enabled`](crate::system::build_info::openmp_enabled) is always
//! `false` and
//! [`openmp_max_num_threads`](crate::system::build_info::openmp_max_num_threads)
//! is always `1`, which is exactly what the source's own non-OpenMP build
//! reports. See `docs/BUILD_INFO_SUPPORT.md`.

/// The label the source uses for anything it could not determine.
pub const UNKNOWN: &str = "unknown";

/// Largest `/etc/os-release` this module will read, in bytes.
#[cfg(target_os = "linux")]
const MAX_OS_RELEASE_BYTES: u64 = 64 * 1024;

/// The operating systems the source distinguishes.
///
/// The source's `SIZE_OF_OPENMS_OS` sentinel exists only to size its parallel
/// name array and has no counterpart: [`OperatingSystem::as_str`] replaces the
/// array and the index into it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OperatingSystem {
    /// The running platform was not recognised.
    #[default]
    Unknown,
    /// Apple macOS.
    MacOs,
    /// Microsoft Windows.
    Windows,
    /// Linux, and — as in the source — any other Unix that is not macOS.
    Linux,
}

impl OperatingSystem {
    /// The source's `OpenMS_OSNames` label for this value.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unknown => UNKNOWN,
            Self::MacOs => "MacOS",
            Self::Windows => "Windows",
            Self::Linux => "Linux",
        }
    }
}

/// The pointer widths the source distinguishes.
///
/// As with [`OperatingSystem`], the source's `SIZE_OF_OPENMS_ARCHITECTURE`
/// sentinel sizes a name array and has no counterpart here.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Architecture {
    /// The pointer width was not recognised.
    #[default]
    Unknown,
    /// Four-byte pointers.
    Bits32,
    /// Eight-byte pointers.
    Bits64,
}

impl Architecture {
    /// The source's `OpenMS_ArchNames` label for this value.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unknown => UNKNOWN,
            Self::Bits32 => "32 bit",
            Self::Bits64 => "64 bit",
        }
    }
}

/// A snapshot of the running operating system.
///
/// [`Default`] gives the source's freshly constructed instance: everything
/// unknown. [`detect`](Self::detect) fills it in, and is the counterpart of the
/// source's static `getOSInfo`.
///
/// The fields are public where the source keeps them private behind
/// `get*AsString` accessors; the enums carry their own labels through
/// [`OperatingSystem::as_str`] and [`Architecture::as_str`], so the accessors
/// would only hide the values a caller may well want to match on.
///
/// ```
/// use openms::system::build_info::{Architecture, OperatingSystem, OsInfo};
///
/// let fresh = OsInfo::default();
/// assert_eq!(fresh.os, OperatingSystem::Unknown);
/// assert_eq!(fresh.arch, Architecture::Unknown);
/// assert_eq!(fresh.os_version, "unknown");
///
/// let running = OsInfo::detect();
/// assert_eq!(running.arch, openms::system::build_info::binary_architecture());
/// ```
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct OsInfo {
    /// The running operating system.
    pub os: OperatingSystem,
    /// The operating system version, or [`UNKNOWN`] when it could not be read.
    pub os_version: String,
    /// The pointer width, derived from the size of a pointer in this binary.
    pub arch: Architecture,
}

impl Default for OsInfo {
    /// Everything unknown, as the source's default constructor leaves it.
    fn default() -> Self {
        Self {
            os: OperatingSystem::Unknown,
            os_version: UNKNOWN.to_string(),
            arch: Architecture::Unknown,
        }
    }
}

impl OsInfo {
    /// Probe the running platform.
    ///
    /// The operating system comes from the compile-time target, as the source's
    /// preprocessor branches do, so a binary run under emulation reports what it
    /// was built for. Any Unix that is not macOS is reported as `Linux`, which
    /// is the source's own behaviour and the subject of its `TODO`. The version
    /// comes from [`os_version_string`], and the architecture from the pointer
    /// width: the source treats four-byte pointers as 32-bit and *everything
    /// else* as 64-bit, which is reproduced here rather than corrected, so this
    /// field is never `Unknown` while [`binary_architecture`] can be.
    pub fn detect() -> Self {
        let os = if cfg!(target_os = "windows") {
            OperatingSystem::Windows
        } else if cfg!(target_os = "macos") {
            OperatingSystem::MacOs
        } else if cfg!(unix) {
            OperatingSystem::Linux
        } else {
            OperatingSystem::Unknown
        };
        let arch = if size_of::<*const ()>() == 4 {
            Architecture::Bits32
        } else {
            Architecture::Bits64
        };
        Self {
            os,
            os_version: os_version_string(),
            arch,
        }
    }
}

/// The pointer width this binary was compiled for, from the size of [`usize`].
///
/// Four bytes is 32-bit, eight is 64-bit and anything else is
/// [`Architecture::Unknown`]. This is the build-time view; [`OsInfo::detect`]
/// reports the same quantity from a pointer's size and never says unknown.
pub fn binary_architecture() -> Architecture {
    match size_of::<usize>() {
        4 => Architecture::Bits32,
        8 => Architecture::Bits64,
        _ => Architecture::Unknown,
    }
}

/// The SIMD instruction sets compiled into this binary, comma separated.
///
/// The source lists the SIMDe build-time probes that were defined, in the order
/// `neon, SSE, SSE2, SSE3, SSE4.1, SSE4.2, AVX, AVX2, FMA`, and returns an empty
/// string when none is. This port reads the corresponding Rust target features,
/// which the compiler enables from the target definition and from
/// `-C target-feature` / `-C target-cpu`, and keeps the source's order and
/// labels. The result is fixed for a given binary, so two calls always agree.
///
/// This is not a claim that the port *uses* those instruction sets: the port has
/// no hand-written SIMD, and the value reports only what the compiler was
/// allowed to emit.
pub fn active_simd_extensions() -> String {
    let mut active: Vec<&'static str> = Vec::new();
    if cfg!(target_feature = "neon") {
        active.push("neon");
    }
    if cfg!(target_feature = "sse") {
        active.push("SSE");
    }
    if cfg!(target_feature = "sse2") {
        active.push("SSE2");
    }
    if cfg!(target_feature = "sse3") {
        active.push("SSE3");
    }
    if cfg!(target_feature = "sse4.1") {
        active.push("SSE4.1");
    }
    if cfg!(target_feature = "sse4.2") {
        active.push("SSE4.2");
    }
    if cfg!(target_feature = "avx") {
        active.push("AVX");
    }
    if cfg!(target_feature = "avx2") {
        active.push("AVX2");
    }
    if cfg!(target_feature = "fma") {
        active.push("FMA");
    }
    active.join(", ")
}

/// The running operating system's version, or [`UNKNOWN`].
///
/// Linux reads `VERSION_ID` from `/etc/os-release` and strips a surrounding
/// pair of double quotes, exactly as the source does; a distribution that omits
/// the key — several rolling releases do — yields [`UNKNOWN`].
///
/// The source additionally queries `RtlGetVersion` from `ntdll.dll` on Windows
/// and `sysctlbyname("kern.osproductversion")` on macOS. Neither is reachable
/// without a platform binding this crate does not take, so both report
/// [`UNKNOWN`] here rather than a guess; the value is a display string, and no
/// behaviour in the port depends on it.
pub fn os_version_string() -> String {
    read_version_id().unwrap_or_else(|| UNKNOWN.to_string())
}

/// Whether the library was built with OpenMP support.
///
/// Always `false`. The port is deliberately serial — the source's 36 files with
/// `#pragma omp` have no threaded counterpart here — so this reports what the
/// source's own non-OpenMP build reports.
pub fn openmp_enabled() -> bool {
    false
}

/// The maximum number of OpenMP threads.
///
/// Always `1`, which is what the source returns when `_OPENMP` is not defined.
/// The `OMP_NUM_THREADS` environment variable has no effect here, and the
/// source's `setOpenMPNumThreads` — a no-op in a non-OpenMP build — is not
/// ported.
pub fn openmp_max_num_threads() -> usize {
    1
}

/// The build configuration this library was compiled with.
///
/// `"Debug"` when debug assertions are on and `"Release"` otherwise. The source
/// returns CMake's `OPENMS_BUILD_TYPE` string verbatim; Cargo bakes no profile
/// name into a library, so this reports the nearest observable fact, which
/// separates Cargo's stock `dev` and `release` profiles but follows the
/// `debug-assertions` setting rather than the profile's name.
pub fn build_type() -> &'static str {
    if cfg!(debug_assertions) {
        "Debug"
    } else {
        "Release"
    }
}

#[cfg(target_os = "linux")]
fn read_version_id() -> Option<String> {
    use std::io::Read;
    let mut text = String::new();
    std::fs::File::open("/etc/os-release")
        .ok()?
        .take(MAX_OS_RELEASE_BYTES)
        .read_to_string(&mut text)
        .ok()?;
    let value = text
        .lines()
        .find_map(|line| line.strip_prefix("VERSION_ID="))?;
    Some(unquoted(value).to_string())
}

#[cfg(not(target_os = "linux"))]
fn read_version_id() -> Option<String> {
    None
}

/// Strip one surrounding pair of double quotes, as the source does.
#[cfg(target_os = "linux")]
fn unquoted(value: &str) -> &str {
    match value.strip_prefix('"').and_then(|v| v.strip_suffix('"')) {
        Some(inner) => inner,
        None => value,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_follow_the_source_name_tables() {
        assert_eq!(OperatingSystem::Unknown.as_str(), "unknown");
        assert_eq!(OperatingSystem::MacOs.as_str(), "MacOS");
        assert_eq!(OperatingSystem::Windows.as_str(), "Windows");
        assert_eq!(OperatingSystem::Linux.as_str(), "Linux");
        assert_eq!(Architecture::Unknown.as_str(), "unknown");
        assert_eq!(Architecture::Bits32.as_str(), "32 bit");
        assert_eq!(Architecture::Bits64.as_str(), "64 bit");
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn version_id_loses_only_a_surrounding_quote_pair() {
        assert_eq!(unquoted("\"22.04\""), "22.04");
        assert_eq!(unquoted("22.04"), "22.04");
        assert_eq!(unquoted("\"22.04"), "\"22.04");
        assert_eq!(unquoted("\""), "\"");
    }
}
