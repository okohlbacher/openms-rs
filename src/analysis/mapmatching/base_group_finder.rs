// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The interface of all element group finding algorithms.
//!
//! Native equivalent of `ANALYSIS/MAPMATCHING/BaseGroupFinder.h`. See
//! `docs/FEATURE_GROUPING_SUPPORT.md`.
//!
//! A group finder takes one or several consensus maps and finds corresponding
//! features across the maps (or within one map), returning one consensus map
//! of the consensus features found. The element indices of the result's
//! handles are the container access indices of the input maps, and their map
//! indices are the indices in the input map vector.
//!
//! The source class is an abstract `DefaultParamHandler` and `ProgressLogger`
//! whose copy and assignment are suppressed. Here it is a trait; each
//! implementation carries its own parameters in the crate's
//! `defaults`/`from_param`/`to_param` shape, and progress reporting, which the
//! source's implementations do not use, is not ported.

use std::collections::BTreeSet;

use crate::error::{Error, Result};
use crate::kernel::features::ConsensusMap;

/// An element group finding algorithm (source `BaseGroupFinder`).
pub trait BaseGroupFinder {
    /// Run the algorithm on `input`, writing the found consensus features into
    /// `result` (source `run`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the input is not valid (source
    /// `Exception::IllegalArgument`).
    fn run(&mut self, input: &[ConsensusMap], result: &mut ConsensusMap) -> Result<()>;
}

/// Check that the column headers of all `maps` use disjoint map identifiers
/// (source `checkIds_`).
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a map identifier occurs in the column
/// headers of two maps (source `Exception::IllegalArgument`, "file ids have to
/// be unique").
pub fn check_ids(maps: &[ConsensusMap]) -> Result<()> {
    let mut used = BTreeSet::new();
    for map in maps {
        for id in map.column_headers.keys() {
            if !used.insert(*id) {
                return Err(Error::InvalidValue("file ids have to be unique".into()));
            }
        }
    }
    Ok(())
}
