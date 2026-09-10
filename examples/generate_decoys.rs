// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{AASequence, DecoyGenerator, Protease};
use openms::format::fasta::{self, FASTAEntry};

fn main() -> openms::Result<()> {
    let generator = DecoyGenerator::with_seed(4711);
    let enzyme = Protease::from_name("Trypsin/P")?;
    let mut database = Vec::new();
    for (accession, sequence) in [("P1", "TESTPEPTIDE"), ("P2", "TESTRPEPTRIDE")] {
        database.push(FASTAEntry {
            identifier: accession.into(),
            description: "target".into(),
            sequence: sequence.into(),
        });
        for (variant, decoy) in generator
            .shuffle(&AASequence::parse(sequence)?, enzyme, 2)?
            .into_iter()
            .enumerate()
        {
            database.push(FASTAEntry {
                identifier: format!("DECOY_{accession}_{variant}"),
                description: format!("{} shuffled variant {variant}", enzyme.name()),
                sequence: decoy.as_str().to_owned(),
            });
        }
    }
    fasta::write(std::io::stdout().lock(), &database)
}
