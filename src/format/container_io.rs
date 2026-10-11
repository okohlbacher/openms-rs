// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! File-system operations for containers that are a directory or an archive
//! of one.
//!
//! The OpenSWATH Parquet containers (`.oswpq`: `TransitionParquetFile`,
//! `OpenSwathOSWParquetWriter`, `OpenSwathOSWParquetReader`) stage their
//! tables in a temporary directory, create, copy, rename and remove
//! directories, and hand archives to
//! [`ZipArchiveFile`](crate::format::zip_archive_file::ZipArchiveFile),
//! whose extraction API already owns a
//! [`TempDir`](crate::system::file::TempDir). This module re-exports
//! the `File.h` operations they use (ported in
//! [`crate::system::file`]), so that those formats reach the file system
//! through the format layer, as their source does through `FORMAT`
//! headers, instead of adding an `analysis` to `system` module edge (which
//! would close a module cycle, see `tools/check_module_cycles.py`).
//!
//! Every item is the `system::file` item itself; nothing is wrapped.

pub use crate::system::file::{
    CopyOptions, TempDir, absolute_path, copy_dir_recursively, exists, is_directory, make_dir,
    remove, remove_dir_recursively, rename,
};
