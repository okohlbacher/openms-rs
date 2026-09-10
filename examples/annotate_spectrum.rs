// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Annotate the pinned SpectrumAnnotator class-test peptide and measured peaks.

use openms::chemistry::{AASequence, SpectrumAnnotator, TheoreticalSpectrumGenerator, ion_naming};
use openms::comparison::{SpectrumAlignment, Tolerance};
use openms::identification::{PeptideHit, PeptideIdentification};
use openms::{MSSpectrum, Peak1D, Result};
use std::io::{BufWriter, Write};

fn main() -> Result<()> {
    // Original measured decimal positions and binary32 intensities. These are
    // experimental-side inputs, not a generated spectrum used as its own target.
    let mut spectrum = MSSpectrum::from_peaks(
        [
            147.113, 204.135, 303.203, 431.262, 518.294, 665.362, 261.16, 348.192, 476.251,
            575.319, 632.341,
        ]
        .into_iter()
        .map(|mz| Peak1D::new(mz, 1.1))
        .collect(),
    );
    let mut identification = PeptideIdentification {
        hits: vec![PeptideHit::new(0.0, 1, 2, AASequence::parse("IFSQVGK")?)?],
        ..Default::default()
    };
    let generator = TheoreticalSpectrumGenerator {
        add_metainfo: true,
        ..Default::default()
    };
    let alignment = SpectrumAlignment {
        tolerance: Tolerance::Absolute(0.1),
        ..Default::default()
    };
    let annotator = SpectrumAnnotator::default();
    annotator.add_ion_match_statistics(
        &mut identification,
        &mut spectrum,
        &generator,
        &alignment,
    )?;
    let hit = &mut identification.hits[0];
    annotator.add_peak_annotations(hit, &spectrum, &generator, &alignment, true)?;
    let mut output = BufWriter::new(std::io::stdout().lock());
    writeln!(output, "mz\tintensity\tion\tcharge")?;
    for peak in &hit.peak_annotations {
        writeln!(
            output,
            "{:.6}\t{:.6}\t{}\t{}",
            peak.mz,
            peak.intensity,
            ion_naming::with_charge(&peak.annotation, peak.charge)?,
            peak.charge
        )?;
    }
    for key in [
        "matched_ion_number",
        "matched_intensity",
        "max_series_type",
        "max_series_size",
        "median_fragment_error",
        "topN_meanfragmenterror",
        "NTermIonCurrentRatio",
        "CTermIonCurrentRatio",
    ] {
        writeln!(output, "# {key}\t{}", hit.metadata[key])?;
    }
    output.flush()?;
    Ok(())
}
