// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Pick and integrate a synthetic trace, or every chromatogram in an input mzML.
//! The output distinguishes raw intensity sum from time-weighted area.

use openms::analysis::peak_integrator::{IntegrationMethod, PeakIntegrator};
use openms::kernel::{ChromatogramPeak, MSChromatogram};
use openms::processing::chromatogram::PeakPickerChromatogram;
use openms::{Error, Result};
use std::io::{self, Write};

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 1 {
        return Err(Error::InvalidValue(
            "usage: integrate_chromatogram [input.mzML]".into(),
        ));
    }
    let chromatograms = if let Some(path) = args.first() {
        read_chromatograms(path)?
    } else {
        let mut trace = MSChromatogram::from_peaks(
            (0..=300)
                .map(|i| {
                    let rt = f64::from(i);
                    let intensity = 10.0
                        + 1000.0 * (-0.5 * ((rt - 100.0) / 9.0).powi(2)).exp()
                        + 600.0 * (-0.5 * ((rt - 200.0) / 12.0).powi(2)).exp();
                    ChromatogramPeak::new(rt, intensity as f32)
                })
                .collect(),
        );
        trace.native_id = "synthetic_transition".into();
        vec![trace]
    };
    if chromatograms.is_empty() {
        return Err(Error::InvalidValue(
            "input contains no chromatograms".into(),
        ));
    }
    let picker = PeakPickerChromatogram::default();
    let sum_integrator = PeakIntegrator::default();
    let area_integrator = PeakIntegrator {
        integration_method: IntegrationMethod::Trapezoid,
        ..Default::default()
    };
    let mut out = io::BufWriter::new(io::stdout().lock());
    writeln!(
        out,
        "chromatogram\tpeak\tapex_seconds\tleft_seconds\tright_seconds\traw_intensity_sum\tarea_intensity_seconds\tbaseline_area\tnet_area"
    )?;
    for (trace_index, trace) in chromatograms.iter().enumerate() {
        let result = picker.pick_chromatogram(trace)?;
        for (peak_index, region) in result.regions.iter().enumerate() {
            // Use original indices, not the rounded f32 leftWidth/rightWidth arrays.
            let left = trace.peaks[region.left_index].rt;
            let right = trace.peaks[region.right_index].rt;
            let sum = sum_integrator.integrate_chromatogram(trace, left, right)?;
            let area = area_integrator.integrate_chromatogram(trace, left, right)?;
            let background = area_integrator.estimate_background_chromatogram(
                trace,
                left,
                right,
                sum.apex_pos,
            )?;
            writeln!(
                out,
                "{trace_index}\t{peak_index}\t{:.6}\t{left:.6}\t{right:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                sum.apex_pos,
                sum.area,
                area.area,
                background.area,
                area.area - background.area
            )?;
        }
    }
    out.flush()?;
    Ok(())
}

#[cfg(feature = "mzml")]
fn read_chromatograms(path: &str) -> Result<Vec<MSChromatogram>> {
    let input = io::BufReader::new(std::fs::File::open(path)?);
    Ok(openms::format::mzml::read(input)?.chromatograms)
}

#[cfg(not(feature = "mzml"))]
fn read_chromatograms(_: &str) -> Result<Vec<MSChromatogram>> {
    Err(Error::Unsupported(
        "reading mzML requires the mzml feature".into(),
    ))
}
