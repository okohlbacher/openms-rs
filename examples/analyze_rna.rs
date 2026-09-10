// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::chemistry::{CoarseIsotopePatternGenerator, NAFragmentType, NASequence};

fn main() -> openms::Result<()> {
    let sequence = NASequence::parse("pA[C*]Gp")?;
    let charge = -2;
    let formula = sequence.formula(NAFragmentType::Full, charge)?;
    let mass = sequence.mono_mass(NAFragmentType::Full, charge)?;
    println!("Sequence: {sequence}");
    println!("Ion formula at charge {charge}: {formula}");
    println!("Ion mass: {mass:.6} Da; m/z: {:.6}", mass / 2.0);
    println!("Suffix after sulfur linkage: {}", sequence.suffix(1)?);
    println!("First five coarse isotope probabilities:");
    let distribution = CoarseIsotopePatternGenerator::default().run(&formula)?;
    for (index, peak) in distribution.peaks().iter().take(5).enumerate() {
        println!("{index}\t{:.8}", peak.probability);
    }
    Ok(())
}
