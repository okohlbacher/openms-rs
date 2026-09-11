// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use crate::{Error, Result, format::FileType};
use std::{io::Read, path::Path};

/// Source document identity and independently stored file provenance.
/// Equality compares only `identifier`, as in the source. Clone and swap retain
/// all three fields. Public field assignments have ordinary Rust costs.
#[derive(Clone, Debug)]
pub struct DocumentIdentifier {
    pub identifier: String,
    pub loaded_file_path: String,
    pub loaded_file_type: FileType,
}

impl Default for DocumentIdentifier {
    fn default() -> Self {
        Self {
            identifier: String::new(),
            loaded_file_path: String::new(),
            loaded_file_type: FileType::Unknown,
        }
    }
}
impl PartialEq for DocumentIdentifier {
    fn eq(&self, other: &Self) -> bool {
        self.identifier == other.identifier
    }
}
impl Eq for DocumentIdentifier {}

impl DocumentIdentifier {
    pub fn new() -> Self {
        Self::default()
    }

    /// Preserve absolute spellings exactly; prefix relative paths with the
    /// current directory without canonicalizing, checking existence or resolving
    /// symlinks. A leading `/` is absolute on every platform, matching source.
    /// Empty input records the current directory. Failure preserves all fields.
    pub fn set_loaded_file_path(&mut self, file_name: &str) -> Result<()> {
        if file_name.len() > crate::system::file::MAX_PATH_BYTES || file_name.contains('\0') {
            return Err(invalid("invalid or oversized document path"));
        }
        let path = Path::new(file_name);
        let value = if path.is_absolute() || file_name.starts_with('/') {
            copy(file_name)?
        } else {
            let absolute = if file_name.is_empty() {
                std::env::current_dir()?
            } else {
                crate::system::file::absolute_path(path)?
            };
            let value = absolute
                .to_str()
                .ok_or_else(|| invalid("absolute document path is not UTF-8"))?;
            if value.len() > crate::system::file::MAX_PATH_BYTES {
                return Err(invalid("absolute document path exceeds limit"));
            }
            // Source filesystem::generic_string uses slash separators for
            // relative-to-absolute conversions, but preserves absolute input.
            let value = copy(value)?;
            #[cfg(windows)]
            let value = value.replace('\\', "/");
            value
        };
        self.loaded_file_path = value;
        Ok(())
    }

    /// Inspect content, independently of the extension and stored path. Reads
    /// at most 64 KiB of plain/decompressed input using the shared detector.
    /// gzip/bzip2 need `file-compression`; ZIP remains unsupported. Recognition
    /// is heuristic and does not validate the document. Failure retains type.
    pub fn set_loaded_file_type(&mut self, file_name: impl AsRef<Path>) -> Result<()> {
        let path = file_name.as_ref();
        if path.as_os_str().len() > crate::system::file::MAX_PATH_BYTES {
            return Err(invalid("document type input path exceeds limit"));
        }
        let mut input = crate::format::path_io::open(path)?.take(65_536);
        let mut preview = Vec::new();
        preview
            .try_reserve_exact(65_536)
            .map_err(|_| invalid("document type preview allocation failure"))?;
        input.read_to_end(&mut preview)?;
        self.loaded_file_type = crate::format::file_handler::type_by_content(&preview);
        Ok(())
    }
}

fn copy(value: &str) -> Result<String> {
    let mut result = String::new();
    result
        .try_reserve_exact(value.len())
        .map_err(|_| invalid("document path allocation failure"))?;
    result.push_str(value);
    Ok(result)
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
