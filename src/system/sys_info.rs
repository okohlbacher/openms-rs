// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Process and system memory reporting of the Core SDK `SYSTEM/SysInfo.h`.
//!
//! The source implements every reading three times, once per operating system,
//! and its own documentation warns that outside Windows the numbers "might be
//! very unreliable". This port is Linux-first: on Linux it reads `/proc`, which
//! needs no libc binding, and on every other platform each reading returns
//! `None` — an explicit "this platform does not report it" — rather than a
//! fabricated zero. [`bytes_to_human_readable`](crate::system::sys_info::bytes_to_human_readable)
//! and [`process_id`](crate::system::sys_info::process_id) are portable and
//! work everywhere.
//!
//! Every memory figure is in kibibytes (the source calls them "KB"), as the
//! source's `size_t` out-parameters are. See `docs/SYS_INFO_SUPPORT.md` for the
//! platform matrix and the `/proc` field derivations.

use crate::{Error, Result};

/// Largest `/proc` record this module will read, in bytes.
#[cfg(target_os = "linux")]
const MAX_PROC_BYTES: u64 = 256 * 1024;

/// Largest event label accepted by [`MemUsage::delta`], in bytes.
///
/// The source imposes no bound on its `event` argument; this port refuses an
/// oversized label before building the report, so a caller-controlled string
/// cannot drive an unbounded allocation.
pub const MAX_EVENT_LABEL_BYTES: usize = 1024 * 1024;

/// Binary unit suffixes, in the source's order, from bytes up to pebibytes.
const UNITS: [&str; 6] = ["byte", "KiB", "MiB", "GiB", "TiB", "PiB"];

/// Convert a byte count to a human readable unit, for example `45.34 MiB`.
///
/// The value is divided by 1024 until it is below 1024, then printed with four
/// significant digits — `2 byte`, `2 KiB`, `1000 byte`, `45.34 MiB` — followed
/// by the unit. The unit name is never pluralised, so one byte reads `1 byte`
/// and so does any other count in that range.
///
/// A count of 2^60 or more exhausts the source's six-entry unit table. The
/// source then returns its literal apology, `Congrats. That's a lot of bytes: `
/// followed by the count, and so does this; the string is reproduced rather
/// than improved because it is the source's observable output.
///
/// ```
/// use openms::system::sys_info::bytes_to_human_readable;
/// assert_eq!(bytes_to_human_readable(2), "2 byte");
/// assert_eq!(bytes_to_human_readable(2048), "2 KiB");
/// assert_eq!(bytes_to_human_readable(2048 << 10), "2 MiB");
/// ```
pub fn bytes_to_human_readable(bytes: u64) -> String {
    let mut value = bytes as f64;
    for unit in UNITS {
        if value < 1024.0 {
            return format!("{} {unit}", four_significant_digits(value));
        }
        value /= 1024.0;
    }
    // Reached from 2^60 bytes upwards, where six units are not enough.
    format!("Congrats. That's a lot of bytes: {bytes}")
}

/// Current memory consumption of this process, in kibibytes.
///
/// `None` means the platform does not report it here, not that the process uses
/// no memory. Linux reads `VmRSS` from `/proc/self/status`; the source computes
/// the same quantity as `statm.resident * sysconf(_SC_PAGESIZE) / 1024`, which
/// the kernel derives from the very counter `VmRSS` publishes, already in
/// kibibytes and without a libc call. On Windows the source reports
/// `WorkingSetSize` and on macOS the Mach `resident_size`; neither is reachable
/// from this crate without a platform binding, so both report `None`.
///
/// The source's own note applies unchanged: outside Windows this figure can be
/// unreliable, and how promptly it reflects a release depends on the kernel and
/// the allocator.
pub fn process_memory_consumption() -> Option<u64> {
    status_field("VmRSS:")
}

/// Peak memory consumption of this process, in kibibytes.
///
/// `None` means the platform does not report it here. Linux reads `VmHWM` —
/// the peak resident set size — from `/proc/self/status`, which is the same
/// counter the source obtains as `getrusage(...).ru_maxrss`, already in
/// kibibytes. Other platforms report `None`.
///
/// The source's `MemUsage` note that peak usage "is only supported on
/// WindowsOS" understates what its own Linux and macOS branches do: both call
/// `getrusage` and both succeed. It is the peak *working set* of Windows that
/// has no exact counterpart elsewhere.
pub fn process_peak_memory_consumption() -> Option<u64> {
    status_field("VmHWM:")
}

/// Physical system memory currently available, in kibibytes.
///
/// `None` means the platform does not report it here. Linux prefers
/// `MemAvailable` from `/proc/meminfo`, so reclaimable page cache counts as
/// available. That is not a substitution: it is the source's own *primary*
/// path, the same file and the same key, which returns on the first match
/// (`SysInfo.cpp:107-133`). Only the fallback differs. Where the running kernel
/// is too old to publish `MemAvailable`, this port falls back to `MemFree` in
/// the same file, while the source falls back to
/// `sysconf(_SC_AVPHYS_PAGES) * sysconf(_SC_PAGESIZE) / 1024`. Both fallbacks
/// answer the narrower question — physically free pages, without the
/// reclaimable cache `MemAvailable` adds — but no claim is made here that they
/// agree to the byte, because `_SC_AVPHYS_PAGES` is answered by whichever libc
/// the source was linked against.
pub fn free_system_memory() -> Option<u64> {
    meminfo_field("MemAvailable:").or_else(|| meminfo_field("MemFree:"))
}

/// The process identifier of the current process.
///
/// The source returns `Int64` because `getpid` and `_getpid` have different
/// signed types; every platform identifier fits in [`u32`], which is what
/// [`std::process::id`] returns, so no sign or width question arises.
pub fn process_id() -> u32 {
    std::process::id()
}

/// Absolute or delta memory usage between two recorded time points.
///
/// Working-set and peak memory are sampled at a "before" and an "after" point.
/// [`delta`](Self::delta) reports the change between them and
/// [`usage`](Self::usage) only the absolute "after" value; both record the
/// second point first if it is still missing. Construction records the first
/// point, so the common use is to build the value, do the work, and print.
///
/// The source uses `0` for "not recorded yet", which makes a genuine reading of
/// zero indistinguishable from a missing one and re-samples every time it is
/// printed. This port uses [`Option`], so "not recorded" and "the platform does
/// not report this" are distinct, and neither is a number.
///
/// ```
/// use openms::system::sys_info::MemUsage;
/// let mut usage = MemUsage::new();
/// let report = usage.delta("loading")?;
/// assert!(report.starts_with("Memory usage (loading): "));
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MemUsage {
    /// Working set at the first time point, in kibibytes.
    pub mem_before: Option<u64>,
    /// Peak working set at the first time point, in kibibytes.
    pub mem_before_peak: Option<u64>,
    /// Working set at the second time point, in kibibytes.
    pub mem_after: Option<u64>,
    /// Peak working set at the second time point, in kibibytes.
    pub mem_after_peak: Option<u64>,
}

impl Default for MemUsage {
    /// Same as [`MemUsage::new`]: constructing records the first time point.
    fn default() -> Self {
        Self::new()
    }
}

impl MemUsage {
    /// Record the first time point, as the source's constructor does.
    pub fn new() -> Self {
        let mut usage = Self {
            mem_before: None,
            mem_before_peak: None,
            mem_after: None,
            mem_after_peak: None,
        };
        usage.before();
        usage
    }

    /// Forget all four readings; [`before`](Self::before) must be called again.
    pub fn reset(&mut self) {
        *self = Self {
            mem_before: None,
            mem_before_peak: None,
            mem_after: None,
            mem_after_peak: None,
        };
    }

    /// Record the first time point.
    pub fn before(&mut self) {
        self.mem_before = process_memory_consumption();
        self.mem_before_peak = process_peak_memory_consumption();
    }

    /// Record the second time point.
    pub fn after(&mut self) {
        self.mem_after = process_memory_consumption();
        self.mem_after_peak = process_peak_memory_consumption();
    }

    /// Report the change in memory usage between the two time points.
    ///
    /// The second point is recorded first if it is still missing. Peak usage is
    /// appended only when the platform reported a non-zero peak, which is the
    /// source's own test for support.
    ///
    /// Differences are printed in mebibytes, truncated towards zero. Where the
    /// working set shrank this port prints a leading `-`; the source builds that
    /// sign into a local string and then overwrites the string instead of
    /// appending to it, so every source report reads as an increase. That is a
    /// defect, not a convention, and is not reproduced.
    ///
    /// A reading the platform does not provide prints as `n/a` rather than as a
    /// number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `event` exceeds
    /// [`MAX_EVENT_LABEL_BYTES`]. The source applies no bound.
    pub fn delta(&mut self, event: &str) -> Result<String> {
        if event.len() > MAX_EVENT_LABEL_BYTES {
            return Err(Error::InvalidValue(
                "memory usage event label exceeds the checked length".into(),
            ));
        }
        if self.mem_after.is_none() {
            self.after();
        }
        let mut report = format!(
            "Memory usage ({event}): {} (working set delta)",
            difference_string(self.mem_before, self.mem_after)
        );
        if self.mem_after_peak.is_some_and(|peak| peak > 0) {
            report.push_str(&format!(
                ", {} (peak working set delta)",
                difference_string(self.mem_before_peak, self.mem_after_peak)
            ));
        }
        Ok(report)
    }

    /// Report the absolute memory usage at the second time point.
    ///
    /// The second point is recorded first if it is still missing. Peak usage is
    /// appended only when the platform reported a non-zero peak.
    pub fn usage(&mut self) -> String {
        if self.mem_after.is_none() {
            self.after();
        }
        let mut report = format!(
            "Memory usage: {} (working set)",
            difference_string(Some(0), self.mem_after)
        );
        if self.mem_after_peak.is_some_and(|peak| peak > 0) {
            report.push_str(&format!(
                ", {} (peak working set)",
                difference_string(Some(0), self.mem_after_peak)
            ));
        }
        report
    }
}

/// The source's `diff_str_`, in mebibytes, with the lost sign restored.
fn difference_string(before: Option<u64>, after: Option<u64>) -> String {
    match (before, after) {
        (Some(before), Some(after)) => {
            let difference = i128::from(after) - i128::from(before);
            // The source takes the absolute value before dividing; for truncation
            // towards zero that is the same magnitude either way.
            let megabytes = difference.abs() / 1024;
            if difference < 0 {
                format!("-{megabytes} MB")
            } else {
                format!("{megabytes} MB")
            }
        }
        _ => "n/a".to_string(),
    }
}

/// Four significant digits with trailing zeros removed, as `std::setprecision(4)`.
///
/// The only caller divides by 1024 until the value is below 1024 and stops at
/// the first quotient, so the argument is either exactly `0.0` or lies in
/// `[1.0, 1024)`. Over that domain the general format never selects scientific
/// notation and the significant-digit count is fixed by the width of the
/// integer part, which is what this computes. Outside it — for a value strictly
/// between zero and one — the significant digits would start at the first
/// non-zero decimal instead, and this would print one digit too few; that case
/// is unreachable from [`bytes_to_human_readable`].
fn four_significant_digits(value: f64) -> String {
    let whole = value.trunc() as u64;
    let integer_digits = match whole {
        0..=9 => 1,
        10..=99 => 2,
        100..=999 => 3,
        _ => 4,
    };
    let decimals = 4usize.saturating_sub(integer_digits);
    let text = format!("{value:.decimals$}");
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        text
    }
}

#[cfg(target_os = "linux")]
fn read_proc(path: &str) -> Option<String> {
    use std::io::Read;
    let mut text = String::new();
    std::fs::File::open(path)
        .ok()?
        .take(MAX_PROC_BYTES)
        .read_to_string(&mut text)
        .ok()?;
    Some(text)
}

/// Parse `<key> <number> kB` from one `/proc` line.
///
/// The unit is required, which the source does not check; every key this module
/// asks for publishes kibibytes, so a line without the suffix is not the record
/// that was wanted.
#[cfg(target_os = "linux")]
fn kibibytes(line: &str, key: &str) -> Option<u64> {
    let mut fields = line.strip_prefix(key)?.split_ascii_whitespace();
    let value: u64 = fields.next()?.parse().ok()?;
    (fields.next() == Some("kB")).then_some(value)
}

#[cfg(target_os = "linux")]
fn status_field(key: &str) -> Option<u64> {
    read_proc("/proc/self/status")?
        .lines()
        .find_map(|line| kibibytes(line, key))
}

#[cfg(not(target_os = "linux"))]
fn status_field(_key: &str) -> Option<u64> {
    None
}

#[cfg(target_os = "linux")]
fn meminfo_field(key: &str) -> Option<u64> {
    read_proc("/proc/meminfo")?
        .lines()
        .find_map(|line| kibibytes(line, key))
}

#[cfg(not(target_os = "linux"))]
fn meminfo_field(_key: &str) -> Option<u64> {
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn four_significant_digits_matches_the_general_format() {
        assert_eq!(four_significant_digits(0.0), "0");
        assert_eq!(four_significant_digits(2.0), "2");
        assert_eq!(four_significant_digits(45.34), "45.34");
        assert_eq!(four_significant_digits(1000.0), "1000");
        assert_eq!(four_significant_digits(1023.0), "1023");
        assert_eq!(four_significant_digits(1.5), "1.5");
    }

    #[test]
    fn difference_string_keeps_the_sign_the_source_discards() {
        assert_eq!(difference_string(Some(0), Some(2048)), "2 MB");
        assert_eq!(difference_string(Some(2048), Some(0)), "-2 MB");
        assert_eq!(difference_string(Some(0), Some(1023)), "0 MB");
        assert_eq!(difference_string(None, Some(1)), "n/a");
        assert_eq!(difference_string(Some(1), None), "n/a");
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn proc_lines_parse_with_their_unit() {
        assert_eq!(kibibytes("VmRSS:\t   12345 kB", "VmRSS:"), Some(12345));
        assert_eq!(
            kibibytes("MemAvailable:   99 kB", "MemAvailable:"),
            Some(99)
        );
        assert_eq!(kibibytes("HugePages_Total:  0", "HugePages_Total:"), None);
        assert_eq!(kibibytes("VmHWM:\t   1 kB", "VmRSS:"), None);
    }
}
