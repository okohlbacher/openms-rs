// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A brute-force stand-in for the source's `KDTreeFeatureMaps`, for the
//! `MapAlignmentAlgorithmKD` tests: the same queries answered by a linear
//! scan instead of a kd-tree. `MapAlignmentAlgorithmKD` does not depend on
//! the order of a neighbourhood, so the answers are those of the source's
//! tree.

use openms::Result;
use openms::analysis::mapmatching::alignment_kd::KdFeatureMaps;
use openms::analysis::transformations::LowessModel;

/// One pooled feature.
#[derive(Clone, Copy, Debug)]
pub struct PoolFeature {
    pub map: usize,
    pub rt: f64,
    pub mz: f64,
    pub intensity: f32,
    pub charge: i32,
}

/// Features of several maps with a cached, transformable retention time.
#[derive(Clone, Debug, Default)]
pub struct Pool {
    pub features: Vec<PoolFeature>,
    pub rt: Vec<f64>,
}

impl Pool {
    pub fn new(features: Vec<PoolFeature>) -> Self {
        let rt = features.iter().map(|f| f.rt).collect();
        Self { features, rt }
    }
}

/// `Math::getTolWindow` at the pin (`MATH/MathFunctions.h`).
pub fn tol_window(value: f64, tol: f64, ppm: bool) -> (f64, f64) {
    if ppm {
        (value - value * tol * 1e-6, value / (1.0 - tol * 1e-6))
    } else {
        (value - tol, value + tol)
    }
}

impl KdFeatureMaps for Pool {
    fn len(&self) -> usize {
        self.features.len()
    }
    fn rt(&self, index: usize) -> f64 {
        self.rt[index]
    }
    fn charge(&self, index: usize) -> i32 {
        self.features[index].charge
    }
    fn map_index(&self, index: usize) -> usize {
        self.features[index].map
    }
    fn neighborhood(
        &self,
        index: usize,
        rt_tol: f64,
        mz_tol: f64,
        mz_ppm: bool,
        include_features_from_same_map: bool,
        max_pairwise_log_fc: f64,
    ) -> Result<Vec<usize>> {
        let (rt_low, rt_high) = tol_window(self.rt[index], rt_tol, false);
        let (mz_low, mz_high) = tol_window(self.features[index].mz, mz_tol, mz_ppm);
        let own_map = self.features[index].map;
        let own_intensity = f64::from(self.features[index].intensity);
        let mut result = Vec::new();
        for (j, feature) in self.features.iter().enumerate() {
            let rt = self.rt[j];
            // The kd-tree region is closed: `low <= value <= high`.
            if rt < rt_low || rt_high < rt || feature.mz < mz_low || mz_high < feature.mz {
                continue;
            }
            if !include_features_from_same_map && feature.map == own_map {
                continue;
            }
            if max_pairwise_log_fc >= 0.0 {
                let fold = (f64::from(feature.intensity) / own_intensity).log10().abs();
                // NaN and +inf fold changes are dropped, as in the source.
                let compatible = fold <= max_pairwise_log_fc;
                if !compatible {
                    continue;
                }
            }
            result.push(j);
        }
        Ok(result)
    }
    fn apply_transformations(&mut self, transformations: &[LowessModel]) -> Result<()> {
        for (i, feature) in self.features.iter().enumerate() {
            self.rt[i] = transformations[feature.map].apply(feature.rt)?;
        }
        Ok(())
    }
}
