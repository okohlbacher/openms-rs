// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Digest a small protein and enumerate fixed/variable peptide modifications.

use openms::Result;
use openms::chemistry::{AASequence, ModifiedPeptideGenerator, ProteaseDigestion};
use std::io::{BufWriter, Write};

fn main() -> Result<()> {
    let protein = AASequence::parse("ACMMKAGHIK")?;
    let generator = ModifiedPeptideGenerator::default();
    let fixed = ModifiedPeptideGenerator::get_modifications(&["Carbamidomethyl (C)"])?;
    let variable = ModifiedPeptideGenerator::get_modifications(&["Oxidation (M)"])?;
    let mut output = BufWriter::new(std::io::stdout().lock());
    writeln!(
        output,
        "start\tend_exclusive\tpeptide\tneutral_mass\tmz_charge_2"
    )?;
    for product in ProteaseDigestion::default().digest(&protein)? {
        let mut peptide = product.sequence;
        generator.apply_fixed_modifications(&fixed, &mut peptide)?;
        for variant in generator.variable_modifications(&variable, &peptide, 2, true)? {
            writeln!(
                output,
                "{}\t{}\t{}\t{:.6}\t{:.6}",
                product.start,
                product.end,
                variant,
                variant.mono_mass()?,
                variant.mz(2)?
            )?;
        }
    }
    output.flush()?;
    Ok(())
}
