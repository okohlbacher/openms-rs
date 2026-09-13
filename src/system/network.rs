// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! One-shot URL download of the Core SDK `SYSTEM/Network.h`.
//!
//! The source is a single static method, `Network::downloadFile`, wrapping
//! `NetworkGetRequest` with a ten-minute timeout, plus a file-static
//! `saveFileName_` that turns a URL into a never-overwriting destination name.
//! Both are ported:
//! [`download_file`](crate::system::network::download_file) and
//! [`save_file_name`](crate::system::network::save_file_name). The naming rules
//! are observable through the file that appears on disk, so the helper is public
//! here where the source keeps it in an anonymous namespace.
//!
//! The module is behind the non-default `network` feature, with the transport
//! injected as a
//! [`HttpTransport`](crate::system::network_get_request::HttpTransport) so that
//! the naming and write paths can be tested against a recorded response instead
//! of the internet. See `docs/NETWORK_SUPPORT.md`.

use crate::system::network_get_request::{HttpTransport, NetworkGetRequest, UreqTransport};
use crate::{Error, Result};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Whole-request timeout the source gives the download, in seconds.
///
/// The source's comment calls it a ten-minute timeout and passes `600` to
/// `NetworkGetRequest::setTimeout`.
pub const DOWNLOAD_TIMEOUT_SECONDS: i32 = 600;

/// Fallback destination name when the URL yields no basename.
pub const FALLBACK_BASENAME: &str = "download";

/// Largest URL [`download_file`] and [`save_file_name`] will look at, in bytes.
///
/// The source bounds neither. This ceiling is checked before anything is copied.
pub const MAX_URL_BYTES: usize = 64 * 1024;

/// Numeric suffixes tried before [`save_file_name`] gives up.
///
/// The source's `while (fs::exists(...)) ++i;` has no upper bound, so a
/// directory already holding every candidate makes it spin, probing the
/// filesystem forever. This port stops and reports instead.
pub const MAX_NAME_SUFFIX: u32 = 10_000;

/// Destination file name for `url` inside `dest_folder`, never overwriting.
///
/// The source's `saveFileName_`, rule for rule:
///
/// 1. anything from the first `?` on is dropped, then anything from the first
///    `#` on — in that order, so a `#` inside a query string is already gone;
/// 2. the basename of what remains is taken;
/// 3. an empty basename becomes [`FALLBACK_BASENAME`];
/// 4. if `dest_folder` holds no such name, that is the answer;
/// 5. otherwise `.0`, `.1`, `.2`, … are appended in order until one is unused.
///
/// Percent escapes are not decoded, by either side, so `%20` survives into the
/// file name.
///
/// The result is a bare file name, never a path: join it onto `dest_folder`
/// yourself, as [`download_file`] does.
///
/// ```
/// use openms::system::network::save_file_name;
/// use std::path::Path;
///
/// let empty = Path::new("");
/// assert_eq!(save_file_name("https://host/a/b/c.tsv", empty)?, "c.tsv");
/// assert_eq!(save_file_name("https://host/c.tsv?v=2#top", empty)?, "c.tsv");
/// assert_eq!(save_file_name("https://host/", empty)?, "download");
/// # Ok::<(), openms::Error>(())
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `url` exceeds [`MAX_URL_BYTES`], when
/// the derived basename is `.`, `..`, or contains a path separator or a NUL
/// byte, or when every suffix up to [`MAX_NAME_SUFFIX`] is taken. The source
/// checks none of these: `.` and `..` reach `std::ofstream` and fail there with
/// an `IOException`, a backslash is a directory separator on Windows and an
/// ordinary character elsewhere, and the suffix search does not terminate.
///
/// # Notes
///
/// The basename is taken by splitting on `/` only, on every platform, because
/// `/` is the path separator of a URL (RFC 3986) whatever the host filesystem
/// does. The source calls `std::filesystem::path::filename`, which additionally
/// treats `\` as a separator on Windows, so `http://h/a\b.txt` names `b.txt`
/// there and `a\b.txt` on Linux. This port names `a\b.txt` everywhere — and
/// then refuses it, per the separator rule above.
///
/// A URL with no `/` after the scheme, such as `https://host`, has `host` as its
/// last `/`-separated component, so that is the file name the source derives and
/// the one derived here.
pub fn save_file_name(url: &str, dest_folder: &Path) -> Result<String> {
    if url.len() > MAX_URL_BYTES {
        return Err(Error::InvalidValue("URL length limit exceeded".into()));
    }
    let path_part = match url.find('?') {
        Some(at) => &url[..at],
        None => url,
    };
    let path_part = match path_part.find('#') {
        Some(at) => &path_part[..at],
        None => path_part,
    };
    let basename = match path_part.rfind('/') {
        Some(at) => &path_part[at + 1..],
        None => path_part,
    };
    let basename = if basename.is_empty() {
        FALLBACK_BASENAME
    } else {
        basename
    };
    if basename == "." || basename == ".." {
        return Err(Error::InvalidValue(
            "URL yields a relative directory name".into(),
        ));
    }
    if basename.contains('/') || basename.contains('\\') || basename.contains('\0') {
        return Err(Error::InvalidValue(
            "URL yields a file name with a path separator".into(),
        ));
    }
    let folder = folder_or_current(dest_folder);
    if !crate::system::file::exists(folder.join(basename)) {
        return Ok(basename.to_owned());
    }
    for suffix in 0..=MAX_NAME_SUFFIX {
        let candidate = format!("{basename}.{suffix}");
        if !crate::system::file::exists(folder.join(&candidate)) {
            return Ok(candidate);
        }
    }
    Err(Error::InvalidValue(
        "no unused download file name below the suffix limit".into(),
    ))
}

/// Download `url` into `download_folder` over the network, returning the path written.
///
/// Equivalent to [`download_file_with`] using
/// [`UreqTransport`], which
/// is the closest counterpart of the source's hard-wired libcurl. This is the
/// only function in the group that opens a socket, and nothing in the crate's
/// own tests calls it.
///
/// # Errors
///
/// As [`download_file_with`].
pub fn download_file(url: &str, download_folder: impl AsRef<Path>) -> Result<PathBuf> {
    download_file_with(&UreqTransport::new(), url, download_folder)
}

/// Download `url` into `download_folder` through `transport`, returning the path written.
///
/// The transfer is synchronous and gets the source's
/// [`DOWNLOAD_TIMEOUT_SECONDS`] deadline. The destination name comes from
/// [`save_file_name`], so an existing file is never overwritten. An empty
/// `download_folder` means the current directory, as the source's `""` → `"./"`.
///
/// Returns the path actually written. The source returns `void` and emits two
/// `OPENMS_LOG_INFO` lines, one of which is that path; returning it is both more
/// useful and the only option here, since `system` may not reach into the
/// logging module.
///
/// # Errors
///
/// * [`Error::Io`] when the request fails — a transport error or an HTTP status
///   of 400 or above — carrying the source's message,
///   `Download of '<url>' failed!. Error: <error>`.
/// * [`Error::Io`] when the destination cannot be opened or written.
/// * [`Error::InvalidValue`] from [`save_file_name`].
///
/// # Notes
///
/// The source leaves a partially written file behind when the write fails, and
/// says so in its documentation. This port removes it, so a failed download
/// leaves the directory as it found it; a removal that itself fails is not
/// allowed to mask the write error.
pub fn download_file_with(
    transport: &dyn HttpTransport,
    url: &str,
    download_folder: impl AsRef<Path>,
) -> Result<PathBuf> {
    let folder = download_folder.as_ref();
    // Fail on an unusable URL before opening a socket, not after.
    if url.len() > MAX_URL_BYTES {
        return Err(Error::InvalidValue("URL length limit exceeded".into()));
    }

    let mut request = NetworkGetRequest::new();
    request.set_url(url);
    request.set_timeout(DOWNLOAD_TIMEOUT_SECONDS);
    request.run(transport);
    if request.has_error() {
        return Err(Error::Io(std::io::Error::other(format!(
            "Download of '{url}' failed!. Error: {}",
            request.error_string()
        ))));
    }

    let filename = folder_or_current(folder).join(save_file_name(url, folder)?);
    let mut file = fs::File::create(&filename).map_err(|error| {
        Error::Io(std::io::Error::new(
            error.kind(),
            format!(
                "Failed to open output file: {}: {error}",
                filename.display()
            ),
        ))
    })?;
    if let Err(error) = file
        .write_all(request.response_binary())
        .and_then(|()| file.flush())
    {
        drop(file);
        let _ = fs::remove_file(&filename);
        return Err(Error::Io(std::io::Error::new(
            error.kind(),
            format!(
                "Failed to write downloaded data to: {}: {error}",
                filename.display()
            ),
        )));
    }
    Ok(filename)
}

/// The source's `download_folder.empty() ? "./" : download_folder`.
///
/// The source then builds `folder + "/" + name` textually, so `"./"` yields
/// `".//name"`; `Path::join` yields `"./name"`, which names the same file.
fn folder_or_current(folder: &Path) -> &Path {
    if folder.as_os_str().is_empty() {
        Path::new(".")
    } else {
        folder
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn query_and_fragment_are_stripped_before_the_basename() {
        let empty = Path::new("");
        assert_eq!(
            save_file_name("https://h/a/b/c.tsv", empty).unwrap(),
            "c.tsv"
        );
        assert_eq!(
            save_file_name("https://h/c.tsv?a=1&b=2", empty).unwrap(),
            "c.tsv"
        );
        assert_eq!(
            save_file_name("https://h/c.tsv#frag", empty).unwrap(),
            "c.tsv"
        );
        // '?' is cut first, so a '#' inside the query never matters.
        assert_eq!(
            save_file_name("https://h/c.tsv?q=#x#y", empty).unwrap(),
            "c.tsv"
        );
    }

    #[test]
    fn an_empty_basename_falls_back_to_download() {
        let empty = Path::new("");
        assert_eq!(save_file_name("https://h/", empty).unwrap(), "download");
        assert_eq!(save_file_name("https://h/?q=1", empty).unwrap(), "download");
        assert_eq!(save_file_name("", empty).unwrap(), "download");
    }

    #[test]
    fn a_host_only_url_names_the_host() {
        // std::filesystem::path("https://host").filename() is "host": the port
        // reproduces that rather than treating the authority specially.
        assert_eq!(
            save_file_name("https://host", Path::new("")).unwrap(),
            "host"
        );
    }

    #[test]
    fn relative_directory_names_and_separators_are_refused() {
        let empty = Path::new("");
        assert!(save_file_name("https://h/a/.", empty).is_err());
        assert!(save_file_name("https://h/a/..", empty).is_err());
        assert!(save_file_name("https://h/a\\b.txt", empty).is_err());
    }

    #[test]
    fn an_oversized_url_is_refused_before_any_work() {
        let url = format!("https://h/{}", "a".repeat(MAX_URL_BYTES));
        assert!(save_file_name(&url, Path::new("")).is_err());
    }
}
