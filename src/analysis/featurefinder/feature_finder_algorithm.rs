// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The feature finder algorithm interface and its shared definitions.
//!
//! Ports `FEATUREFINDER/FeatureFinderAlgorithm.h` (header only: the fitting
//! [`Summary`] and the abstract [`FeatureFinderAlgorithm`]) and
//! `FEATUREFINDER/FeatureFinderDefs.h` (header only). See
//! `docs/FEATUREFINDER_MODELS_SUPPORT.md`.
//!
//! `FeatureFinderDefs.h` declares, token for token, the struct that
//! `FeatureFinderAlgorithmPicked.h` also declares; that declaration is already
//! ported as [`crate::analysis::feature_finder_picked::defs`], and the two
//! source headers cannot even be included together (a redefinition). The
//! module [`defs`] re-exports that port under this header's name rather than
//! duplicating it.

use crate::Result;
use crate::error::Error;
use crate::kernel::MSExperiment;
use crate::kernel::features::FeatureMap;
use crate::param::Param;
use std::collections::BTreeMap;

/// Source `FeatureFinderDefs.h`: `IndexPair`, `IndexSet`, `ChargedIndexSet`,
/// the `Flag` enum (`UNUSED`, `USED`) and the `NoSuccessor` exception.
pub mod defs {
    pub use crate::analysis::feature_finder_picked::defs::{
        ChargedIndexSet, Flag, IndexPair, IndexSet, NoSuccessor,
    };
}

/// Summary of fitting results (source struct `Summary`).
#[derive(Clone, Debug, PartialEq)]
pub struct Summary {
    /// Count of exceptions, by name (source `std::map<std::string, UInt>`).
    pub exception: BTreeMap<String, u32>,
    /// Number of fits without an exception, initially 0.
    pub no_exceptions: u32,
    /// Count of the m/z models used, by name.
    pub mz_model: BTreeMap<String, u32>,
    /// Count of the m/z standard deviations used, in ascending order of the
    /// deviation (source `std::map<float, UInt>`); see
    /// [`Summary::count_mz_stdev`].
    pub mz_stdev: Vec<(f32, u32)>,
    /// Count of the charges used, indexed by charge.
    pub charge: Vec<u32>,
    /// Mean correlation, initially 0 (box plot of the correlations).
    pub corr_mean: f64,
    /// Maximum correlation, initially 0.
    pub corr_max: f64,
    /// Minimum correlation, initially 1.
    pub corr_min: f64,
}

impl Default for Summary {
    /// The source's initial values: no counts, `corr_mean` 0, `corr_max` 0,
    /// `corr_min` 1.
    fn default() -> Self {
        Self {
            exception: BTreeMap::new(),
            no_exceptions: 0,
            mz_model: BTreeMap::new(),
            mz_stdev: Vec::new(),
            charge: Vec::new(),
            corr_mean: 0.0,
            corr_max: 0.0,
            corr_min: 1.0,
        }
    }
}

impl Summary {
    /// Count one use of the m/z standard deviation `stdev`, as the source's
    /// `++mz_stdev[stdev]` does: an entry is found by `std::map`'s
    /// equivalence (neither is less than the other), else inserted in order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for NaN, which is not a strict weak
    /// ordering key (undefined behaviour in the source's map).
    pub fn count_mz_stdev(&mut self, stdev: f32) -> Result<()> {
        if stdev.is_nan() {
            return Err(Error::InvalidValue(
                "a NaN m/z standard deviation cannot be a map key".into(),
            ));
        }
        let position = self.mz_stdev.partition_point(|(key, _)| *key < stdev);
        match self.mz_stdev.get_mut(position) {
            Some((key, count)) if !(stdev < *key) => *count = count.saturating_add(1),
            _ => self.mz_stdev.insert(position, (stdev, 1)),
        }
        Ok(())
    }
}

/// The message the source's default `setSeeds` throws with.
pub const SEEDS_NOT_SUPPORTED: &str =
    "The used feature detection algorithm does not support user-specified seed lists!";

/// A feature finder algorithm (source abstract class `FeatureFinderAlgorithm`,
/// a `DefaultParamHandler` and `ProgressLogger`).
///
/// The source stores pointers to the input map and output features with
/// `setData` and then calls `run()`; here [`FeatureFinderAlgorithm::run`]
/// takes both. The source's copy constructor and assignment operator are
/// private and unimplemented, so implementors need not be `Clone`.
pub trait FeatureFinderAlgorithm {
    /// The default parameters (source `getDefaultParameters`, which returns
    /// `defaults_`; an algorithm with sub-algorithms includes theirs).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    fn default_parameters(&self) -> Result<Param>;

    /// Provide a seed list (source `setSeeds`). The default accepts only an
    /// empty one.
    ///
    /// # Errors
    ///
    /// The default returns [`Error::InvalidValue`] with
    /// [`SEEDS_NOT_SUPPORTED`] for a non-empty `seeds` (the source's
    /// `Exception::IllegalArgument`).
    fn set_seeds(&mut self, seeds: &FeatureMap) -> Result<()> {
        if seeds.features.is_empty() {
            Ok(())
        } else {
            Err(Error::InvalidValue(SEEDS_NOT_SUPPORTED.into()))
        }
    }

    /// Run the algorithm on `map`, writing into `features` (source `setData`
    /// followed by `run`).
    ///
    /// # Errors
    ///
    /// As the implementing algorithm documents.
    fn run(&mut self, map: &MSExperiment, features: &mut FeatureMap) -> Result<()>;
}
