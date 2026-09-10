// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Named modifications -> peptide masses -> coarse isotopes -> annotated b/y ions.
//! Run without arguments for a modified peptide, or supply one sequence argument.

use openms::chemistry::{
    AASequence, CoarseIsotopePatternGenerator, CoarseMassMode, TheoreticalSpectrumGenerator,
};
use openms::{Error, Result};
use std::io::{BufWriter, Write};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() > 1 {
        return Err(Error::InvalidValue(
            "usage: peptide_analysis ['modified peptide sequence']".into(),
        ));
    }
    let sequence = args
        .first()
        .map(String::as_str)
        .unwrap_or("(Acetyl)AC(Carbamidomethyl)M(Oxidation)K");
    let peptide = AASequence::parse(sequence)?;
    if peptide.is_empty() {
        return Err(Error::InvalidValue("peptide must contain a residue".into()));
    }
    let formula = peptide.formula()?;
    let isotopes =
        CoarseIsotopePatternGenerator::new(Some(5), CoarseMassMode::Approximate)?.run(&formula)?;
    let spectrum = TheoreticalSpectrumGenerator {
        add_metainfo: true,
        ..Default::default()
    }
    .generate(&peptide, 1, 2, Some(3))?;
    // Annotations are generated together with peaks and remain aligned after sorting.
    let names = &spectrum
        .string_data_arrays
        .iter()
        .find(|a| a.name == "IonNames")
        .ok_or_else(|| Error::InvalidValue("missing generated ion names".into()))?
        .data;
    let charges = &spectrum
        .integer_data_arrays
        .iter()
        .find(|a| a.name == "Charges")
        .ok_or_else(|| Error::InvalidValue("missing generated charges".into()))?
        .data;

    let mut output = BufWriter::new(std::io::stdout().lock());
    writeln!(output, "peptide\t{peptide}")?;
    writeln!(output, "neutral_formula\t{formula}")?;
    writeln!(output, "neutral_mono_mass\t{:.6}", peptide.mono_mass()?)?;
    writeln!(output, "mz_2plus\t{:.6}", peptide.mz(2)?)?;
    writeln!(
        output,
        "\n# Neutral coarse envelope; probabilities normalized over retained bins"
    )?;
    writeln!(output, "isotope_bin\tneutral_mass\tprobability")?;
    for (index, peak) in isotopes.peaks().iter().enumerate() {
        writeln!(output, "{index}\t{:.6}\t{:.8}", peak.mass, peak.probability)?;
    }
    writeln!(output, "\n# Monoisotopic b/y fragments, charges 1 and 2")?;
    writeln!(output, "ion\tcharge\tmz\tintensity")?;
    for ((peak, name), charge) in spectrum.peaks.iter().zip(names).zip(charges) {
        writeln!(
            output,
            "{name}\t{charge}\t{:.6}\t{:.3}",
            peak.mz, peak.intensity
        )?;
    }
    output.flush()?;
    Ok(())
}
