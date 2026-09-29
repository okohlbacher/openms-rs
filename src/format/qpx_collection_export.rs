// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! All-or-nothing writing of a QPX collection (`-out_qpx`).
//!
//! Ports `QPXCollectionExport` (`FORMAT/QPXCollectionExport.h` and
//! `QPXCollectionExport.cpp`) in part: its [`Transaction`] scope guard. See
//! `docs/QPX_SUPPORT.md`.
//!
//! A QPX collection is three files written in a fixed order --
//! `quantms.feature.parquet`, then `quantms.psm.parquet`, then
//! `quantms.pg.parquet`. Each view refuses input it cannot represent, but it
//! does so before opening *its own* file, so a refusal in the second or third
//! view leaves the earlier files on disk: a partial collection that looks
//! like a successful export. The source handles that from two sides:
//!
//! * `requireExportable(cmap, design)` runs every deterministic refusal of
//!   all three views before the first write. **Not ported**: it needs the
//!   consensus-map and protein-group exporters (`ConsensusMapArrowExport`'s
//!   `requireUnambiguousIdentities` and `requireResolvableIdRuns`,
//!   `ProteinGroupArrowExport`'s quantification-unit checks) and
//!   `ExperimentalDesign`, none of which exist here yet.
//! * [`Transaction`] removes whatever this export created when it does not
//!   run to completion -- a refusal only discovered while building a row, or
//!   an I/O failure.
//!
//! The source marks the class experimental ("may change in future
//! versions").

use crate::concept::log_stream::{LogLevel, with_global_log};
use std::io::Write;
use std::path::{Path, PathBuf};

/// The three files of a QPX collection, in the order the views write them.
pub const QPX_COLLECTION_FILES: [&str; 3] = [
    "quantms.feature.parquet",
    "quantms.psm.parquet",
    "quantms.pg.parquet",
];

/// Scope guard that leaves no partial QPX collection behind.
///
/// Ports `QPXCollectionExport::Transaction`. Construct it before the first
/// write and call [`commit`](Self::commit) once the last one has succeeded.
/// On any other exit -- a view returning an error, a panic unwinding, an
/// early return -- dropping the guard removes the collection files *this
/// export created*, so the export is all-or-nothing. Enumerating the failure
/// sites instead would silently miss the next one.
///
/// Files that already existed when the guard was constructed are never
/// removed: exporting into a directory that still holds an earlier
/// collection is the normal way to re-run a pipeline, and a failed re-run
/// must not take the previous result with it. Such a run does leave the
/// directory mixed -- the views that ran overwrote their files before the
/// failure -- so the rollback names the surviving files in its warning rather
/// than pretending the directory is still coherent.
///
/// ```no_run
/// use openms::format::qpx_collection_export::Transaction;
/// # fn export_views(_: &std::path::Path) -> openms::Result<()> { Ok(()) }
/// let out_qpx = std::path::Path::new("out_qpx");
/// let qpx = Transaction::new(out_qpx);
/// export_views(out_qpx)?; // feature, psm and pg views, in that order
/// qpx.commit();
/// # Ok::<(), openms::Error>(())
/// ```
///
/// The source's guard is neither copyable nor movable; this one is not
/// `Clone`, and moving it is harmless in Rust (the destructor runs once).
#[derive(Debug)]
pub struct Transaction {
    directory: PathBuf,
    /// Collection files that were already on disk when the guard started;
    /// never removed.
    preexisting: Vec<PathBuf>,
    committed: bool,
}

/// One line the rollback logs, at the level the source logs it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RollbackMessage {
    /// `Info` for a removed file, `Error` for one that could not be removed,
    /// `Warn` for the files that predate the guard.
    pub level: LogLevel,
    /// The message, as the source writes it.
    pub text: String,
}

impl Transaction {
    /// Start guarding a QPX collection directory.
    ///
    /// Ports `Transaction(directory)`: records which of
    /// [`QPX_COLLECTION_FILES`] already exist in `directory`, so that only
    /// the ones this export creates can be removed again. The source joins
    /// the names to the directory as strings (`directory + "/quantms..."`);
    /// this joins paths, so an empty `directory` means the current directory
    /// rather than the file system root.
    pub fn new(directory: impl AsRef<Path>) -> Self {
        let directory = directory.as_ref().to_path_buf();
        let preexisting = QPX_COLLECTION_FILES
            .iter()
            .map(|name| directory.join(name))
            .filter(|path| path.exists())
            .collect();
        Self {
            directory,
            preexisting,
            committed: false,
        }
    }

    /// The guarded directory.
    pub fn directory(&self) -> &Path {
        &self.directory
    }

    /// The collection files that predate the guard, which a rollback keeps.
    pub fn preexisting(&self) -> &[PathBuf] {
        &self.preexisting
    }

    /// Keep the collection: the export completed. Ports `commit()`.
    pub fn commit(mut self) {
        self.committed = true;
    }

    /// Roll back now and return what the source would log: remove the
    /// collection files this export created, and warn about the ones that
    /// predate the guard.
    ///
    /// Native: the source only rolls back in its destructor, which logs these
    /// lines; dropping an uncommitted guard here does the same, writing them
    /// to the global log. This explicit form hands them to the caller
    /// instead. It never fails; a file that cannot be removed is reported.
    pub fn roll_back(mut self) -> Vec<RollbackMessage> {
        self.committed = true;
        self.remove_created()
    }

    fn remove_created(&self) -> Vec<RollbackMessage> {
        let mut messages = Vec::new();
        let mut kept = Vec::new();
        for name in QPX_COLLECTION_FILES {
            let path = self.directory.join(name);
            if !path.exists() {
                continue;
            }
            if self.preexisting.contains(&path) {
                kept.push(path.display().to_string());
                continue;
            }
            if std::fs::remove_file(&path).is_ok() {
                messages.push(RollbackMessage {
                    level: LogLevel::Info,
                    text: format!(
                        "QPXCollectionExport: removed {} -- the QPX collection was not written in \
                         full.",
                        path.display()
                    ),
                });
            } else {
                messages.push(RollbackMessage {
                    level: LogLevel::Error,
                    text: format!(
                        "QPXCollectionExport: failed to remove {}, which is part of an incomplete \
                         QPX collection. Delete the directory before using it.",
                        path.display()
                    ),
                });
            }
        }
        if !kept.is_empty() {
            messages.push(RollbackMessage {
                level: LogLevel::Warn,
                text: format!(
                    "QPXCollectionExport: the QPX collection was not written in full, and {} \
                     already held a collection before this export. These files predate the \
                     guard and were left in place, but a view that ran before the failure may \
                     have overwritten them, so the directory may now mix two exports: {}",
                    self.directory.display(),
                    kept.join(", ")
                ),
            });
        }
        messages
    }
}

impl Drop for Transaction {
    /// Remove the collection files this export created, unless the guard
    /// was committed, and log what was done. Never panics: it runs during
    /// unwinding when a writer panicked, and a failure to log is ignored.
    fn drop(&mut self) {
        if self.committed {
            return;
        }
        for message in self.remove_created() {
            let _ = with_global_log(message.level, |log| {
                log.write_all(message.text.as_bytes())?;
                log.write_all(b"\n")
            });
        }
    }
}
