// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A failed `FASTAFile` load or store must not block every other loader.
//!
//! Progress nesting is process-wide, as the source's `static int
//! recursion_depth_` is (`CONCEPT/ProgressLogger.h:105`), and the port caps it
//! at [`MAX_PROGRESS_DEPTH`]. Every other reader reaches its logger through a
//! [`ProgressReporter`](openms::concept::progress_logger::ProgressReporter),
//! which, when a call ends, *abandons* the sections the call started and left
//! open. An abandoned section still counts for the logger that abandoned it —
//! so that logger stays indented as the source's would — but not for any
//! other: another logger dispatches at every *open* section plus only its own
//! abandoned ones. `ProgressNesting::depth` keeps counting abandoned sections
//! for as long as the logger that abandoned them lives, and drops them when it
//! goes.
//!
//! `FASTAFile::load` and `FASTAFile::store` called their own `ProgressLogger`
//! directly, bypassing the reporter. Their failed sections were therefore never
//! abandoned: they stayed **open**, and an open section counts for every logger.
//! After [`MAX_PROGRESS_DEPTH`] failures through one `FASTAFile`, every
//! progress-reporting load in the process — through any other logger, of any
//! format — was refused with "progress nesting limit exceeded".
//!
//! Found by adversarial verification of the Phase 3 progress repair, which had
//! fixed the same defect on every other reader.
//!
//! This is its own test binary on purpose: it drives the process-wide nesting
//! to its limit, which would disturb any test sharing the process, and it is a
//! single sequential test because tests in one binary run on parallel threads.

use openms::concept::progress_logger::{MAX_PROGRESS_DEPTH, ProgressLogger, ProgressNesting};
use openms::format::{dta2d, fasta::FASTAFile};
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// A progress-reporting load through a fresh logger, of a different format:
/// exactly the call the leak refused.
fn fresh_dta2d_load_succeeds() -> bool {
    dta2d::load_with_progress(
        data("progress_format_readers/DTA2DFile_test_1.dta2d"),
        &Default::default(),
        &mut ProgressLogger::new(),
    )
    .is_ok_and(|experiment| !experiment.spectra.is_empty())
}

#[test]
fn failed_fasta_loads_and_stores_leave_every_other_loader_alone() {
    let before = ProgressNesting::global().depth();
    let entries = FASTAFile::new()
        .load(data("fasta_source.fasta"))
        .expect("the FASTA fixture loads");
    assert!(!entries.is_empty());
    {
        let mut loader = FASTAFile::new();
        for _ in 0..MAX_PROGRESS_DEPTH {
            assert!(loader.load(data("no_such.fasta")).is_err());
        }
        assert!(
            fresh_dta2d_load_succeeds(),
            "failed FASTA loads left sections open that block every other logger"
        );

        // The store path had the same defect. A destination that cannot be
        // created fails after its section has started.
        let mut storer = FASTAFile::new();
        let unwritable = data("no_such_directory/out.fasta");
        for _ in 0..MAX_PROGRESS_DEPTH {
            assert!(storer.store(&unwritable, &entries).is_err());
        }
        assert!(
            fresh_dta2d_load_succeeds(),
            "failed FASTA stores left sections open that block every other logger"
        );
    }
    // Once the objects that failed are gone, nothing they abandoned remains.
    assert_eq!(ProgressNesting::global().depth(), before);
}
