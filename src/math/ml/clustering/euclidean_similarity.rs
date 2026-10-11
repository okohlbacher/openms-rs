// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A similarity for points in the plane derived from their Euclidean distance.
//!
//! Native equivalent of `ML/CLUSTERING/EuclideanSimilarity.h`. See
//! `docs/CLUSTERING_SUPPORT.md`.

use crate::{Error, Result};

/// Compare functor for 2D points (source `EuclideanSimilarity`).
///
/// Each point is a pair of `f32` coordinates `(x, y)`. The similarity is
/// `1 - distance / scale`, so it lies in `[0, 1]` for every pair whose
/// distance does not exceed the scale. The scale defaults to 1.
///
/// The source's copy constructor, assignment and destructor map to `Clone`,
/// `Copy` and ordinary drop.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EuclideanSimilarity {
    /// Divisor of the distance (source `scale_`).
    scale: f32,
}

impl Default for EuclideanSimilarity {
    /// Scale 1 (source default constructor).
    fn default() -> Self {
        Self { scale: 1.0 }
    }
}

impl EuclideanSimilarity {
    /// A functor with scale 1 (source default constructor).
    pub const fn new() -> Self {
        Self { scale: 1.0 }
    }

    /// Similarity of the points `a` and `b` (source two-argument
    /// `operator()`): `1 - sqrtf(dx*dx + dy*dy) / scale`, all in `f32`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the scale is 0 and the scaling therefore
    /// inapplicable (source `Exception::DivisionByZero`).
    pub fn similarity(&self, a: (f32, f32), b: (f32, f32)) -> Result<f32> {
        if self.scale == 0.0 {
            return Err(Error::InvalidValue(
                "EuclideanSimilarity: division by zero, the scale is 0".into(),
            ));
        }
        let dx = a.0 - b.0;
        let dy = a.1 - b.1;
        Ok(1.0 - (dx * dx + dy * dy).sqrt() / self.scale)
    }

    /// Self-similarity of `c` (source one-argument `operator()`), computed as
    /// [`similarity(c, c)`](Self::similarity).
    ///
    /// The source documents that this "will yield 0"; the distance of a point
    /// to itself is 0, so the similarity it computes is 1 (NaN for a NaN
    /// coordinate), and this port returns the same.
    ///
    /// # Errors
    ///
    /// As [`similarity`](Self::similarity).
    pub fn self_similarity(&self, c: (f32, f32)) -> Result<f32> {
        self.similarity(c, c)
    }

    /// Sets the scale (source `setScale`).
    ///
    /// Choose it so that the greatest distance in the chosen set scales to 1,
    /// that is `x` = the greatest possible distance in the set, so that
    /// similarities are correctly calculated from distances. A scale of 0 is
    /// accepted here and refused by [`similarity`](Self::similarity), as in
    /// the source.
    pub fn set_scale(&mut self, x: f32) {
        self.scale = x;
    }
}
