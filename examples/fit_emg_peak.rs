// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reconstruct a cropped synthetic EMG peak and compare integrated areas.

use openms::Result;
use openms::analysis::emg::{EmgGradientDescent, EmgParameters};
use openms::analysis::peak_integrator::{IntegrationMethod, PeakIntegrator};
use openms::kernel::{ChromatogramPeak, MSChromatogram};
use std::io::{self, Write};

fn main() -> Result<()> {
    let parameters = EmgParameters {
        h: 1000.0,
        mu: 100.0,
        sigma: 1.0,
        tau: 2.0,
    };
    let positions: Vec<_> = (0..=35).map(|i| 97.0 + f64::from(i) * 0.2).collect();
    let curve = EmgGradientDescent {
        compute_additional_points: false,
        ..Default::default()
    }
    .apply_parameters(&positions, parameters)?;
    let input = MSChromatogram::from_peaks(
        curve
            .positions
            .iter()
            .zip(&curve.intensities)
            .map(|(&rt, &y)| ChromatogramPeak::new(rt, y as f32))
            .collect(),
    );
    let fitter = EmgGradientDescent::default();
    let fitted = fitter.fit_chromatogram(&input, None, None)?;
    let measured = PeakIntegrator {
        integration_method: IntegrationMethod::Trapezoid,
        ..Default::default()
    }
    .integrate_chromatogram(&input, 97.0, 104.0)?;
    let reconstructed = PeakIntegrator {
        integration_method: IntegrationMethod::Trapezoid,
        emg: Some(fitter),
        ..Default::default()
    }
    .integrate_chromatogram(&input, 97.0, 104.0)?;
    let estimate = &fitted.estimate;
    let p = estimate.parameters;
    let mut out = io::BufWriter::new(io::stdout().lock());
    writeln!(
        out,
        "input_points\tfitted_points\titerations\tbest_iteration\tconverged\ttraining_mse"
    )?;
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{:.6}",
        input.len(),
        fitted.chromatogram.len(),
        estimate.iterations,
        estimate.best_iteration,
        estimate.converged,
        estimate.loss
    )?;
    writeln!(out, "h\tmu_seconds\tsigma_seconds\ttau_seconds")?;
    writeln!(out, "{:.6}\t{:.6}\t{:.6}\t{:.6}", p.h, p.mu, p.sigma, p.tau)?;
    writeln!(
        out,
        "cropped_area\treconstructed_sampled_area\tcomplete_synthetic_area"
    )?;
    // The normalized exponential convolution preserves the Gaussian's total area.
    writeln!(
        out,
        "{:.6}\t{:.6}\t{:.6}",
        measured.area,
        reconstructed.area,
        parameters.h * parameters.sigma * (2.0 * std::f64::consts::PI).sqrt()
    )?;
    out.flush()?;
    Ok(())
}
