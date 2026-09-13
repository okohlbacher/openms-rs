// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Lexical path helpers of the Core SDK `SYSTEM/PathUtils.h`.
//!
//! The source header carries exactly two entry points and no `.cpp`: an inline
//! `PathUtils::basename` that deliberately duplicates `File::basename` so that
//! callers need not include `File.h`, and `to_path`, which converts a UTF-8
//! `std::string` into a `std::filesystem::path` without going through the
//! Windows ANSI code page.
//!
//! Both are purely lexical: neither touches the filesystem, normalises `.`/`..`
//! segments, resolves symbolic links, or collapses repeated separators.
//! See `docs/PATH_UTILS_SUPPORT.md` for the member-by-member mapping.

use crate::{Error, Result};
use std::path::PathBuf;

/// Largest path accepted by [`to_path`], in bytes.
///
/// The source applies no bound; this port refuses oversized input before
/// allocating, matching `system::file`'s ceiling so the two agree.
pub const MAX_PATH_BYTES: usize = super::file::MAX_PATH_BYTES;

/// Return the part of `file` after the last `/` or `\`, on every platform.
///
/// A path without a separator is returned whole, and a path that ends in a
/// separator — including a path that is nothing but a separator — yields the
/// empty string. The source obtains that by `substr(find_last_of("\\/") + 1)`
/// and relies on the unsigned wraparound of `npos + 1` to zero, a quirk its own
/// comment calls out; this port reaches the same three results without
/// arithmetic on a sentinel.
///
/// This delegates to [`file::basename`](crate::system::file::basename): the
/// source's `PathUtils::basename` is a verbatim copy of `File::basename` kept
/// only so that the header has no dependency on `File.h`, and duplicating the
/// body here would let the two drift apart.
///
/// ```
/// use openms::system::path_utils::basename;
/// assert_eq!(basename("/data/run1.mzML"), "run1.mzML");
/// assert_eq!(basename("run1.mzML"), "run1.mzML");
/// assert_eq!(basename("C:\\data\\run1.mzML"), "run1.mzML");
/// assert_eq!(basename("/data/"), "");
/// assert_eq!(basename("/"), "");
/// assert_eq!(basename(""), "");
/// ```
pub fn basename(file: &str) -> &str {
    super::file::basename(file)
}

/// Convert a UTF-8 string into a [`PathBuf`], preserving its bytes exactly.
///
/// The source exists to work around a C++ problem that does not arise here:
/// `std::filesystem::path(std::string)` interprets its argument in the active
/// Windows code page rather than as UTF-8, so the source constructs from
/// `std::u8string` and falls back to the code-page constructor when the bytes
/// are not valid UTF-8 (a filename taken from `argv` under an ANSI locale).
/// A Rust `&str` is UTF-8 by construction and [`PathBuf`] stores it as UTF-8 on
/// Unix and as WTF-8 on Windows, so neither the re-encoding nor its fallback
/// has a counterpart: the conversion is lossless on every platform and the
/// source's `std::system_error` branch is unreachable.
///
/// The returned path round-trips: for any accepted `s`, the path's string form
/// is byte-identical to `s`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `s` exceeds [`MAX_PATH_BYTES`] or
/// contains a NUL byte. The source checks neither; a NUL cannot reach any
/// platform filesystem call intact, so rejecting it here turns a failure at
/// every later use into one checked error at the conversion.
///
/// ```
/// use openms::system::path_utils::to_path;
/// let path = to_path("日本語.mzML")?;
/// assert_eq!(path.to_str(), Some("日本語.mzML"));
/// # Ok::<(), openms::Error>(())
/// ```
pub fn to_path(s: &str) -> Result<PathBuf> {
    if s.len() > MAX_PATH_BYTES {
        return Err(Error::InvalidValue("filesystem path limit exceeded".into()));
    }
    if s.as_bytes().contains(&0) {
        return Err(Error::InvalidValue("filesystem path contains NUL".into()));
    }
    Ok(PathBuf::from(s))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basename_handles_trailing_and_bare_separators() {
        assert_eq!(basename("a/b/c.txt"), "c.txt");
        assert_eq!(basename("a\\b\\c.txt"), "c.txt");
        assert_eq!(basename("a/b\\c.txt"), "c.txt");
        assert_eq!(basename("c.txt"), "c.txt");
        assert_eq!(basename("dir/"), "");
        assert_eq!(basename("dir\\"), "");
        assert_eq!(basename("/"), "");
        assert_eq!(basename("\\"), "");
        assert_eq!(basename(""), "");
    }

    #[test]
    fn to_path_refuses_nul_and_oversized_input() {
        assert!(to_path("a\0b").is_err());
        assert!(to_path(&"a".repeat(MAX_PATH_BYTES + 1)).is_err());
        assert!(to_path(&"a".repeat(MAX_PATH_BYTES)).is_ok());
    }
}
