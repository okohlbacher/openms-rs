// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Run with no arguments for the bundled sample, or supply input.mgf/input.mzML/
//! input.dta and optionally an output.mgf path. Output contains filtered spectra.

use openms::format::{dta, mgf};
use openms::processing::{
    NLargest, NormalizationMethod, Normalizer, SpectrumFilter, ThresholdMower,
};
use openms::{Error, MSExperiment, Result};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() > 2 {
        return Err(Error::InvalidValue(
            "usage: process_spectra [input.mgf|input.mzML|input.dta] [output.mgf]".into(),
        ));
    }
    let mut experiment = if let Some(path) = args.first() {
        let reader = BufReader::new(File::open(path)?);
        match Path::new(path)
            .extension()
            .and_then(|e| e.to_str())
            .map(str::to_ascii_lowercase)
            .as_deref()
        {
            Some("mgf") => mgf::read(reader)?,
            Some("dta") => MSExperiment {
                spectra: vec![dta::read(reader)?],
                ..Default::default()
            },
            #[cfg(feature = "mzml")]
            Some("mzml") => openms::format::mzml::read(reader)?,
            _ => {
                return Err(Error::Unsupported(
                    "input must be .mgf, .dta, or .mzML (with mzml feature)".into(),
                ));
            }
        }
    } else {
        MSExperiment {
            spectra: vec![dta::read(
                &include_bytes!("../tests/data/Transformers_tests.dta")[..],
            )?],
            ..Default::default()
        }
    };
    let initial: usize = experiment.spectra.iter().map(|s| s.len()).sum();
    ThresholdMower { threshold: 10.0 }.filter_experiment(&mut experiment)?;
    NLargest { n: 100 }.filter_experiment(&mut experiment)?;
    Normalizer {
        method: NormalizationMethod::ToTic,
    }
    .filter_experiment(&mut experiment)?;
    for spectrum in &mut experiment.spectra {
        spectrum.sort_by_position()?;
    }
    let retained: usize = experiment.spectra.iter().map(|s| s.len()).sum();
    println!(
        "{} spectra: {initial} input peaks, {retained} retained peaks",
        experiment.spectra.len()
    );
    for (i, spectrum) in experiment.spectra.iter().enumerate().take(10) {
        println!(
            "spectrum {i}: MS{}, {} peaks, normalized TIC {:.6}",
            spectrum.ms_level,
            spectrum.len(),
            spectrum.calculate_tic()
        );
    }
    if let Some(path) = args.get(1) {
        if Path::new(path)
            .extension()
            .and_then(|e| e.to_str())
            .is_none_or(|e| !e.eq_ignore_ascii_case("mgf"))
        {
            return Err(Error::InvalidValue("output path must end with .mgf".into()));
        }
        // Serialize first: a format validation error must not truncate a file.
        let mut encoded = Vec::new();
        mgf::write(&mut encoded, &experiment)?;
        let mut output = BufWriter::new(File::create(path)?);
        output.write_all(&encoded)?;
        output.flush()?;
        println!("Wrote {path}");
    }
    Ok(())
}
