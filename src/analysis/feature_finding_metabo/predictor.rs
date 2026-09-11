// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// Fixed-model RBF prediction follows LIBSVM (Copyright 2000-2023 Chang and Lin).
// Its retained notice is resources/metabolite_isotope_models/LIBSVM_COPYRIGHT.txt.

use crate::{Error, Result};

#[derive(Clone, Copy, Debug)]
pub(super) enum Model {
    Noise2,
    Noise5,
}
struct SupportVector {
    coefficient: f64,
    features: [f64; 4],
}
struct ModelData {
    gamma: f64,
    rho: f64,
    centers: [f64; 4],
    scales: [f64; 4],
    support: &'static [SupportVector],
}
static NOISE2: ModelData = include!("../../../resources/metabolite_isotope_models/noise2.rs");
static NOISE5: ModelData = include!("../../../resources/metabolite_isotope_models/noise5.rs");

/// Predict source label 2 using the two exact immutable FFM model tables.
/// The caller supplies capped mass and three isotope/mono ratios (missing=0).
/// All support-vector work is precharged; successful prediction allocates no heap
/// storage, so the caller's byte ledger is unchanged. Raw inputs must be finite.
pub(super) fn predict(
    model: Model,
    raw_features: [f64; 4],
    remaining_work: &mut usize,
    _remaining_bytes: &mut usize,
) -> Result<bool> {
    let model = match model {
        Model::Noise2 => &NOISE2,
        Model::Noise5 => &NOISE5,
    };
    Ok(decision(model, raw_features, remaining_work)? > 0.0)
}
fn decision(model: &ModelData, raw: [f64; 4], remaining: &mut usize) -> Result<f64> {
    let cost = model
        .support
        .len()
        .checked_mul(32)
        .and_then(|n| n.checked_add(16))
        .ok_or_else(resource)?;
    *remaining = remaining.checked_sub(cost).ok_or_else(resource)?;
    if raw.iter().any(|x| !x.is_finite()) {
        return Err(Error::InvalidValue(
            "FFM predictor features must be finite".into(),
        ));
    }
    let mut features = [0.0; 4];
    for i in 0..4 {
        features[i] = (raw[i] - model.centers[i]) / model.scales[i];
    }
    let mut sum = 0.0;
    for vector in model.support {
        let mut squared_distance = 0.0;
        for (x, y) in features.iter().zip(vector.features) {
            let difference = x - y;
            squared_distance += difference * difference;
        }
        // Keep source operations separate, in file order. Finite raw values
        // can overflow scaling/distance to infinity: exp(-infinity)=0 is the
        // defined source RBF limit, and still yields a finite decision.
        let kernel = (-model.gamma * squared_distance).exp();
        sum += vector.coefficient * kernel;
    }
    Ok(sum - model.rho)
}
fn resource() -> Error {
    Error::InvalidValue("FFM predictor shared work limit exceeded".into())
}

#[cfg(test)]
#[path = "predictor_tests.rs"]
mod tests;
