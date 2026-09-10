// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Digest RNA, enumerate modified products and generate annotated negative ions.
use openms::chemistry::{
    ModifiedNASequenceGenerator, NASequence, NucleicAcidSpectrumGenerator, RNaseDigestion,
    RibonucleotideDB,
};

fn main() -> openms::Result<()> {
    let input: NASequence = "AGUACG".parse()?;
    let database = RibonucleotideDB::global();
    let alternatives = [database.get("s4U")?, database.get("m3U")?];
    let modifications = ModifiedNASequenceGenerator::default();
    let spectrum_generator = NucleicAcidSpectrumGenerator {
        add_metainfo: true,
        add_first_prefix_ion: true,
        ..Default::default()
    };
    println!("RNA: {input}; RNase_T1, no missed cleavages");
    for product in RNaseDigestion::default().digest_with_positions(&input)? {
        for variant in
            modifications.variable_modifications(&alternatives, &product.sequence, 1, true)?
        {
            let spectrum = spectrum_generator.generate(&variant, -1, -1)?;
            println!(
                "{}..{}\t{}\t{} annotated singly negative b/y peaks",
                product.start,
                product.end,
                variant,
                spectrum.len()
            );
            for (peak, name) in spectrum
                .peaks
                .iter()
                .zip(&spectrum.string_data_arrays[0].data)
            {
                println!("  {name}\t{:.6}\t{:.1}", peak.mz, peak.intensity);
            }
        }
    }
    Ok(())
}
