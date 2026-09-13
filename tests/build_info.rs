// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `SYSTEM/BuildInfo.h`: running-platform identity and build configuration.
//!
//! The class test's eight sections are covered here. Its expectations are
//! themselves platform-derived — `sizeof(size_t)`, a set of accepted OS names,
//! a deterministic SIMD string — so the same derivations are reproduced rather
//! than pinned to whichever host happens to run them.

use openms::system::build_info::{
    Architecture, OperatingSystem, OsInfo, UNKNOWN, active_simd_extensions, binary_architecture,
    build_type, openmp_enabled, openmp_max_num_threads, os_version_string,
};

/// Class-test section `OpenMSOSInfo()`.
///
/// A freshly constructed instance reports everything as `unknown`.
#[test]
fn a_default_instance_reports_everything_unknown() {
    let info = OsInfo::default();
    assert_eq!(info.os, OperatingSystem::Unknown);
    assert_eq!(info.arch, Architecture::Unknown);
    assert_eq!(info.os_version, "unknown");
    assert_eq!(info.os.as_str(), "unknown");
    assert_eq!(info.arch.as_str(), "unknown");
    assert_eq!(UNKNOWN, "unknown");
}

/// Class-test section `~OpenMSOSInfo()`.
///
/// The source's destructor is implicit; the port's value owns only its version
/// string, so a clone outlives the original's scope with the same contents.
#[test]
fn an_instance_owns_only_its_version_string() {
    let escaped = {
        let info = OsInfo::detect();
        let copy = info.clone();
        drop(info);
        copy
    };
    assert_eq!(escaped, OsInfo::detect());
}

/// Class-test section `std::string getOSAsString() const`.
#[test]
fn the_operating_system_label_follows_the_source_name_table() {
    assert_eq!(OsInfo::default().os.as_str(), "unknown");
    assert_eq!(OperatingSystem::MacOs.as_str(), "MacOS");
    assert_eq!(OperatingSystem::Windows.as_str(), "Windows");
    assert_eq!(OperatingSystem::Linux.as_str(), "Linux");
}

/// Class-test section `std::string getArchAsString() const`.
#[test]
fn the_architecture_label_follows_the_source_name_table() {
    assert_eq!(OsInfo::default().arch.as_str(), "unknown");
    assert_eq!(Architecture::Bits32.as_str(), "32 bit");
    assert_eq!(Architecture::Bits64.as_str(), "64 bit");
}

/// Class-test section `std::string getOSVersionAsString() const`.
#[test]
fn the_version_of_a_default_instance_is_unknown_and_a_probe_is_never_empty() {
    assert_eq!(OsInfo::default().os_version, "unknown");
    // The source test only requires a non-empty string; the fallback satisfies it.
    assert!(!os_version_string().is_empty());
    assert_eq!(os_version_string(), OsInfo::detect().os_version);
}

/// Class-test section `static std::string getBinaryArchitecture()`.
///
/// Derived from `sizeof(size_t)`, exactly as the source test derives its own
/// expectation, so it holds on any host.
#[test]
fn the_binary_architecture_comes_from_the_size_of_usize() {
    let expected = match size_of::<usize>() {
        4 => Architecture::Bits32,
        8 => Architecture::Bits64,
        _ => Architecture::Unknown,
    };
    assert_eq!(binary_architecture(), expected);
    assert_ne!(binary_architecture(), Architecture::Unknown);
}

/// Class-test section `static OpenMSOSInfo getOSInfo()`.
#[test]
fn a_probe_recognises_the_running_platform() {
    let info = OsInfo::detect();
    // The architecture is derived from the pointer width, so it must agree.
    assert_eq!(info.arch, binary_architecture());
    // The running OS must be recognised on any supported build platform.
    assert!(["Windows", "MacOS", "Linux"].contains(&info.os.as_str()));
    assert_ne!(info.os, OperatingSystem::Unknown);
    // A version string is always reported, at worst the `unknown` fallback.
    assert!(!info.os_version.is_empty());
}

/// Class-test section `static std::string getActiveSIMDExtensions()`.
///
/// Build-dependent content, but a stable string for a given binary.
#[test]
fn the_simd_extension_list_is_deterministic_and_well_formed() {
    let first = active_simd_extensions();
    assert_eq!(first, active_simd_extensions());
    if !first.is_empty() {
        // The source joins with ", " and never emits an empty element.
        assert!(!first.starts_with(", "));
        assert!(!first.ends_with(", "));
        for name in first.split(", ") {
            assert!(!name.is_empty());
            assert!(
                [
                    "neon", "SSE", "SSE2", "SSE3", "SSE4.1", "SSE4.2", "AVX", "AVX2", "FMA"
                ]
                .contains(&name),
                "unexpected SIMD label {name}"
            );
        }
    }
    // Every 64-bit x86 target has SSE2 in its baseline, and every aarch64 has neon.
    #[cfg(target_arch = "x86_64")]
    assert!(first.contains("SSE2"));
    #[cfg(target_arch = "aarch64")]
    assert!(first.contains("neon"));
}

/// `OpenMSBuildInfo`: the port is serial, so its OpenMP answers are constants.
#[test]
fn openmp_is_absent_and_reports_a_single_thread() {
    assert!(!openmp_enabled());
    // The source's non-OpenMP build gives the same two answers. There is no
    // thread pool to size, so `OMP_NUM_THREADS` has nothing to act on and the
    // source's `setOpenMPNumThreads` has no counterpart to call.
    assert_eq!(openmp_max_num_threads(), 1);
}

/// `OpenMSBuildInfo::getBuildType`, mapped onto Cargo's debug-assertion setting.
#[test]
fn the_build_type_is_one_of_the_two_cmake_names_the_source_reports() {
    let expected = if cfg!(debug_assertions) {
        "Debug"
    } else {
        "Release"
    };
    assert_eq!(build_type(), expected);
    assert!(["Debug", "Release"].contains(&build_type()));
}
