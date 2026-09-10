// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Inspect a natural-isotope prefix and an enriched carbon population.

use openms::Result;
use openms::chemistry::{EmpiricalFormula, FineIsotopeIterator, element};
use std::io::{BufWriter, Write};

fn main() -> Result<()> {
    let mut output = BufWriter::new(std::io::stdout().lock());
    writeln!(
        output,
        "population\tneutral_mass\tprobability\tlog_probability"
    )?;
    let glucose = EmpiricalFormula::parse("C6H12O6")?;
    // Taking a prefix does not require materializing the complete support.
    for item in FineIsotopeIterator::from_formula(&glucose)?.take(5) {
        let config = item?;
        writeln!(
            output,
            "natural_glucose\t{:.9}\t{:.12}\t{:.9}",
            config.mass, config.probability, config.log_probability
        )?;
    }

    let masses: Vec<_> = element("C")
        .unwrap()
        .isotopes()
        .iter()
        .map(|i| i.mass)
        .collect();
    // Two carbon atoms, each with a 75% carbon-13 abundance. Custom weights
    // retain f64 precision; the iterator does not normalize the supplied row.
    let enriched = FineIsotopeIterator::from_isotopes(&[2], &[masses], &[vec![0.25, 0.75]])?
        .with_absolute_threshold(0.1)?;
    for item in enriched {
        let config = item?;
        writeln!(
            output,
            "enriched_carbon\t{:.9}\t{:.12}\t{:.9}",
            config.mass, config.probability, config.log_probability
        )?;
    }
    output.flush()?;
    Ok(())
}
