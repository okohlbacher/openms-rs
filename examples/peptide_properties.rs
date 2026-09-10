// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Inspect source-model peptide properties, including modified parent residues.

use openms::chemistry::{AAIndex, AASequence, HydrophobicityProfile, IsoelectricPoint};
use openms::{Error, Result};
use std::io::{BufWriter, Write};

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 1 {
        return Err(Error::InvalidValue(
            "usage: peptide_properties ['peptide sequence']".into(),
        ));
    }
    let sequences: Vec<&str> = if let Some(sequence) = args.first() {
        vec![sequence.as_str()]
    } else {
        vec![
            "(Acetyl)AC(Carbamidomethyl)M(Oxidation)K",
            "PEPTIDER",
            "GTVVTGR",
        ]
    };
    let calculator = IsoelectricPoint::default();
    let mut output = BufWriter::new(std::io::stdout().lock());
    writeln!(
        output,
        "peptide\tcharge_pH7\tpI_Lehninger\tGRAVY\tGB_500K\tGB_100K"
    )?;
    for sequence in sequences {
        let peptide = AASequence::parse(sequence)?;
        // Charge/pI use terminal annotation presence. These source property
        // models otherwise retain parent-residue values despite modifications.
        writeln!(
            output,
            "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            peptide,
            calculator.compute_charge(&peptide, 7.0)?,
            calculator.compute_pi(&peptide)?,
            HydrophobicityProfile::compute_gravy(&peptide)?,
            AAIndex::calculate_gb(&peptide, 500.0)?,
            AAIndex::calculate_gb(&peptide, 100.0)?,
        )?;
    }
    output.flush()?;
    Ok(())
}
