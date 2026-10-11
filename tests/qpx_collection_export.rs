// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `QPXCollectionExport::Transaction` (`FORMAT/QPXCollectionExport.h`).
//!
//! The section `[EXTRA] Transaction removes a partially written collection
//! unless committed` of the upstream class test `QPXCollectionExport_test.cpp`
//! (core `bc9cc12`), step by step; the rollback messages are the source's
//! (`QPXCollectionExport.cpp`, source review). The class test's
//! `requireExportable` sections need the consensus-map and protein-group
//! exporters, which are not ported.
#![cfg(feature = "arrow")]

use openms::concept::log_stream::LogLevel;
use openms::format::qpx_collection_export::{QPX_COLLECTION_FILES, Transaction};
use std::path::{Path, PathBuf};

fn directory(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("qpx_collection_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn files(dir: &Path) -> Vec<PathBuf> {
    QPX_COLLECTION_FILES
        .iter()
        .map(|name| dir.join(name))
        .collect()
}

fn write_stubs(files: &[PathBuf], n: usize) {
    for file in &files[..n] {
        std::fs::write(file, "stub").unwrap();
    }
}

#[test]
fn collection_file_names_and_order() {
    assert_eq!(
        QPX_COLLECTION_FILES,
        [
            "quantms.feature.parquet",
            "quantms.psm.parquet",
            "quantms.pg.parquet"
        ]
    );
}

/// `[EXTRA] Transaction removes a partially written collection unless committed`.
#[test]
fn transaction_removes_a_partial_collection_unless_committed() {
    let dir = directory("class");
    let files = files(&dir);

    // Abandoned after the first view.
    {
        let _txn = Transaction::new(&dir);
        write_stubs(&files, 1);
        assert!(files[0].exists());
    }
    assert!(!files[0].exists());

    // Abandoned after two of three views.
    {
        let _txn = Transaction::new(&dir);
        write_stubs(&files, 2);
        assert!(files[0].exists() && files[1].exists());
    }
    assert!(!files[0].exists() && !files[1].exists());

    // Committed: the complete collection is kept.
    {
        let txn = Transaction::new(&dir);
        write_stubs(&files, 3);
        txn.commit();
    }
    assert!(files.iter().all(|f| f.exists()));

    // A failed export does not take a previous, complete collection with it.
    {
        let _txn = Transaction::new(&dir);
    }
    assert!(files.iter().all(|f| f.exists()));

    // Mixed: only the newly created file is rolled back.
    std::fs::remove_file(&files[2]).unwrap();
    {
        let txn = Transaction::new(&dir);
        assert_eq!(txn.preexisting(), &files[..2]);
        std::fs::write(&files[2], "stub").unwrap();
        assert!(files[2].exists());
    }
    assert!(files[0].exists() && files[1].exists());
    assert!(!files[2].exists());
    std::fs::remove_file(&files[0]).unwrap();
    std::fs::remove_file(&files[1]).unwrap();

    // The guard cleans up while a panic unwinds, and does not panic itself.
    let unwound = std::panic::catch_unwind(|| {
        let _txn = Transaction::new(&dir);
        write_stubs(&files, 3);
        panic!("boom");
    });
    assert!(unwound.is_err());
    assert!(files.iter().all(|f| !f.exists()));

    // A directory that never received a file is left alone.
    {
        let _txn = Transaction::new(&dir);
    }
    assert!(dir.exists());
    std::fs::remove_dir_all(&dir).unwrap();
}

/// The explicit rollback returns the lines the source's destructor logs.
#[test]
fn rollback_messages_are_the_sources() {
    let dir = directory("messages");
    let files = files(&dir);
    write_stubs(&files, 1);
    let txn = Transaction::new(&dir);
    assert_eq!(txn.directory(), dir.as_path());
    std::fs::write(&files[1], "stub").unwrap();
    let messages = txn.roll_back();
    assert_eq!(messages.len(), 2);
    assert_eq!(messages[0].level, LogLevel::Info);
    assert_eq!(
        messages[0].text,
        format!(
            "QPXCollectionExport: removed {} -- the QPX collection was not written in full.",
            files[1].display()
        )
    );
    assert_eq!(messages[1].level, LogLevel::Warn);
    assert_eq!(
        messages[1].text,
        format!(
            "QPXCollectionExport: the QPX collection was not written in full, and {} already \
             held a collection before this export. These files predate the guard and were left \
             in place, but a view that ran before the failure may have overwritten them, so the \
             directory may now mix two exports: {}",
            dir.display(),
            files[0].display()
        )
    );
    assert!(files[0].exists());
    assert!(!files[1].exists());

    // Nothing to do, nothing to say.
    let messages = Transaction::new(directory("quiet")).roll_back();
    assert!(messages.is_empty());
    std::fs::remove_dir_all(&dir).unwrap();
}
