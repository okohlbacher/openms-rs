// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The processor features this binary was compiled to require, checked once
//! before the first instruction that needs them.
//!
//! **This module has no counterpart in the C++ SDK.** CMake builds OpenMS for
//! the baseline architecture with no `-march`, so a C++ tool binary runs on any
//! processor its operating system runs on and never has to ask. This port does
//! not: `.cargo/config.toml` builds x86_64 with `-C target-feature=+fma`,
//! because the ported glibc `exp`, `log` and `powf` are built out of
//! [`f64::mul_add`] and a baseline x86-64 target has to call out of line for
//! every one of them (docs/BENCHMARKS.md section 4, docs/FMA_BUILD_FLAG.md).
//!
//! The price of that flag is a minimum processor. Without a guard, running such
//! a binary on an older CPU kills it with `SIGILL` and no explanation at
//! whatever arithmetic happens to come first. [`unsupported_cpu`] turns that
//! into a diagnosis: it reads CPUID, and the caller prints the message and
//! exits with a documented status. [`crate::cli::run`] calls it as its very
//! first statement — before it even reads the command line — and
//! [`crate::cli::run_with`] calls it again for a caller that drives a tool in
//! process, so every ported TOPP executable is covered before it touches a
//! floating-point number.
//!
//! **`std::arch::is_x86_feature_detected!` cannot be used here.** It is
//! documented to answer `true` without asking the processor whenever the
//! feature is already enabled at compile time, which is precisely this case: a
//! `+fma` build folds `is_x86_feature_detected!("fma")` to a constant and the
//! whole guard disappears. That was measured, not assumed — under `rustc -O`
//! the probe function is eliminated outright in a `+fma` build and keeps its
//! CPUID cache lookup in a baseline one (`docs/FMA_BUILD_FLAG.md`, section *The
//! standard library cannot answer this*). The question is therefore put to
//! `raw-cpuid`, which always executes `cpuid`.
//!
//! What that does and does not guarantee is in `docs/FMA_BUILD_FLAG.md`:
//! `fma` implies `avx`, so the whole crate is compiled with VEX encoding, and a
//! processor old enough to lack AVX as well can fault on an instruction that is
//! not a fused multiply-add. The guard is exact for the processors that have
//! AVX but not FMA, and measured — not guaranteed by construction — for the
//! older ones.
//!
//! The check costs one relaxed atomic load after the first call — the standard
//! library caches the CPUID result — and is compiled out entirely on a build
//! that does not set the flag, because [`build_requires_fma`] is then a
//! compile-time `false`.

/// What a binary that requires FMA prints on a processor that has none.
///
/// The text names the requirement, the processor generations that satisfy it
/// and the exact command that rebuilds without it. It is a constant rather than
/// a formatted string so that nothing on the path to it can allocate or do
/// arithmetic.
///
/// The command works because `RUSTFLAGS` *replaces* the `rustflags` of
/// `.cargo/config.toml` rather than adding to them: one variable is the whole
/// opt-out. `docs/FMA_BUILD_FLAG.md` records the verification of that command.
pub const FMA_UNSUPPORTED_MESSAGE: &str = "\
This build of OpenMS requires a processor with the x86-64 FMA3 instruction set \
(Intel Haswell, AMD Piledriver or newer). This processor does not have it, so \
the executable would die on an illegal instruction as soon as it did any \
arithmetic, and it stops here instead.

It was built with `-C target-feature=+fma` from `.cargo/config.toml`. To build \
it for this processor, run this in the OpenMS4-R checkout:

    RUSTFLAGS=\"-C target-feature=-fma\" cargo build --release --locked

`RUSTFLAGS` replaces the flags in `.cargo/config.toml` rather than adding to \
them, so that one command is the whole change. The results are bitwise \
identical either way; only the speed of the ported glibc `exp`, `log` and \
`powf` differs.";

/// Whether this build may emit fused-multiply-add instructions the processor
/// has to provide.
///
/// True only on x86 or x86_64 compiled with the `fma` target feature, which is
/// what `.cargo/config.toml` sets for x86_64 and what
/// `RUSTFLAGS="-C target-feature=-fma"` removes. Every other target either has
/// fused multiply-add in its baseline — aarch64 does, and rejects the `fma`
/// feature name — or never has the compiler emit it, so there is nothing to
/// check and this is a compile-time `false`.
///
/// This is a property of *this crate's* compilation, not of the process: a
/// dependency compiled with different flags would not change it. All of the
/// port's fused arithmetic is in this crate, and `rustflags` apply to the whole
/// build, so the two agree.
///
/// ```
/// # use openms::system::cpu_features::build_requires_fma;
/// // A statement about the build, so it is the same answer every time.
/// assert_eq!(build_requires_fma(), build_requires_fma());
/// ```
pub fn build_requires_fma() -> bool {
    cfg!(all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_feature = "fma"
    ))
}

/// Whether the processor running this code provides those instructions.
///
/// On x86 and x86_64 this executes `cpuid` through `raw-cpuid` and reads bit 12
/// of `ECX` in leaf 1, the architectural FMA bit. It is not
/// [`std::arch::is_x86_feature_detected!`], which would be a compile-time
/// `true` on exactly the builds that need the answer; see the module
/// documentation.
///
/// A processor whose CPUID does not even report leaf 1 — nothing since the 486
/// — is taken to provide it. The guard exists to turn a crash into an
/// explanation, never to invent a refusal on a machine it cannot diagnose.
///
/// On every other architecture the question does not arise, because nothing
/// there requires FMA at build time (see [`build_requires_fma`]), and this
/// reports `true` rather than a guess about a processor it has no instruction
/// to interrogate.
///
/// ```
/// # use openms::system::cpu_features::{build_requires_fma, cpu_provides_fma};
/// // Whatever this host is, a binary that runs at all satisfies its own build.
/// assert!(!build_requires_fma() || cpu_provides_fma());
/// ```
pub fn cpu_provides_fma() -> bool {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        match raw_cpuid::CpuId::new().get_feature_info() {
            Some(info) => info.has_fma(),
            None => true,
        }
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        true
    }
}

/// The message to print when this binary cannot run on this processor, or
/// `None` when it can.
///
/// `Some` exactly when [`build_requires_fma`] and not [`cpu_provides_fma`].
/// Callers print the message to standard error and exit; they must not
/// continue, because the very next floating-point operation may be the illegal
/// instruction. [`crate::cli::run_with`] does this with
/// [`ExitCode::InternalError`](crate::cli::ExitCode::InternalError), the status
/// `docs/TOPP_CLI_SUPPORT.md` records for a failure of the build rather than of
/// the input.
///
/// ```
/// # use openms::system::cpu_features::unsupported_cpu;
/// // This process is running, so its own processor satisfies its own build.
/// assert!(unsupported_cpu().is_none());
/// ```
pub fn unsupported_cpu() -> Option<&'static str> {
    message_for(build_requires_fma(), cpu_provides_fma())
}

/// The decision itself, with both answers supplied.
///
/// Split out so the four combinations can be tested on any host: neither
/// [`build_requires_fma`] nor [`cpu_provides_fma`] can be made to say anything
/// other than what this machine and this build say.
fn message_for(build_requires_fma: bool, cpu_provides_fma: bool) -> Option<&'static str> {
    if build_requires_fma && !cpu_provides_fma {
        Some(FMA_UNSUPPORTED_MESSAGE)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The whole truth table, with the "does the CPU have FMA" answer injected.
    #[test]
    fn only_a_requiring_build_on_a_processor_without_fma_is_refused() {
        assert_eq!(message_for(true, false), Some(FMA_UNSUPPORTED_MESSAGE));
        assert_eq!(message_for(true, true), None);
        assert_eq!(message_for(false, false), None);
        assert_eq!(message_for(false, true), None);
    }

    /// The message is useless without the command that fixes it, so pin both
    /// the requirement it names and the command verbatim.
    #[test]
    fn the_message_names_the_requirement_and_the_exact_rebuild_command() {
        assert!(FMA_UNSUPPORTED_MESSAGE.contains("FMA3"));
        assert!(FMA_UNSUPPORTED_MESSAGE.contains("Intel Haswell, AMD Piledriver or newer"));
        assert!(
            FMA_UNSUPPORTED_MESSAGE
                .contains("RUSTFLAGS=\"-C target-feature=-fma\" cargo build --release --locked")
        );
        // The opt-out only works because RUSTFLAGS replaces rather than adds.
        assert!(FMA_UNSUPPORTED_MESSAGE.contains("replaces the flags"));
        // One line per idea; a wrapped paragraph is unreadable on a terminal.
        assert!(FMA_UNSUPPORTED_MESSAGE.lines().count() >= 5);
        assert!(!FMA_UNSUPPORTED_MESSAGE.ends_with('\n'));
    }

    /// A binary that is running has already executed whatever its build emits.
    #[test]
    fn this_process_satisfies_its_own_build() {
        assert_eq!(unsupported_cpu(), None);
        assert!(!build_requires_fma() || cpu_provides_fma());
    }

    /// x86 is the only architecture the requirement can apply to.
    #[test]
    fn no_target_but_x86_requires_anything() {
        if !cfg!(any(target_arch = "x86", target_arch = "x86_64")) {
            assert!(!build_requires_fma());
            assert!(cpu_provides_fma());
        }
    }

    /// The detector must be the processor's answer and not the build's. On
    /// Linux the kernel publishes the same CPUID bit as a `flags` word, so the
    /// two can be compared; they must agree on any build, with the flag or
    /// without it.
    #[cfg(all(target_os = "linux", any(target_arch = "x86", target_arch = "x86_64")))]
    #[test]
    fn the_detector_agrees_with_what_the_kernel_reports() {
        let Ok(text) = std::fs::read_to_string("/proc/cpuinfo") else {
            return; // No procfs: nothing to compare against, and nothing wrong.
        };
        let Some(flags) = text
            .lines()
            .find_map(|line| line.strip_prefix("flags").and_then(|l| l.split_once(':')))
        else {
            return;
        };
        let kernel_says = flags.1.split_whitespace().any(|flag| flag == "fma");
        assert_eq!(
            cpu_provides_fma(),
            kernel_says,
            "CPUID and /proc/cpuinfo disagree about FMA"
        );
    }
}
