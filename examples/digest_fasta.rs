// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Streaming FASTA -> tryptic peptides -> doubly protonated m/z.
use openms::chemistry::{AASequence, ProteaseDigestion};
use openms::format::fasta::FastaReader;
use openms::{Error, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() > 1 {
        return Err(Error::InvalidValue(
            "usage: digest_fasta [input.fasta]".into(),
        ));
    }
    let reader: Box<dyn BufRead> = if let Some(path) = args.first() {
        Box::new(BufReader::new(File::open(path)?))
    } else {
        Box::new(b">example Demonstration protein\nDFPIANGERACDEK\n".as_slice())
    };
    let digestion = ProteaseDigestion {
        missed_cleavages: 1,
        min_length: 5,
        ..Default::default()
    };
    println!("protein\tstart\tend\tpeptide\tneutral_mass\tmz_2plus");
    for record in FastaReader::new(reader) {
        let record = record?;
        let sequence = AASequence::parse(&record.sequence)?;
        for peptide in digestion.digest(&sequence)? {
            println!(
                "{}\t{}\t{}\t{}\t{:.6}\t{:.6}",
                record.identifier,
                peptide.start,
                peptide.end,
                peptide.sequence,
                peptide.sequence.mono_mass()?,
                peptide.sequence.mz(2)?
            );
        }
    }
    Ok(())
}
