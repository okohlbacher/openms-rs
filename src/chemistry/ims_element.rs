// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Named IMS isotope distribution, retaining the source's independent sequence label.

use super::ims_isotope_distribution::{MAX_TEXT, Work, finite, invalid};
use super::{IMSIsotopeDistribution, IMSIsotopeOptions};
use crate::Result;

const MAX_LABEL: usize = 1024 * 1024;

/// An owned IMS element. Names and sequences are arbitrary bounded UTF-8 text.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IMSElement {
    name: String,
    sequence: String,
    isotopes: IMSIsotopeDistribution,
}

impl IMSElement {
    /// Source IMS constant, intentionally distinct from newer chemistry constants.
    pub const ELECTRON_MASS_IN_U: f64 = 0.00054858;

    /// A name and an empty nominal-mass distribution; sequence initially equals name.
    pub fn new(name: &str, nominal_mass: u32) -> Result<Self> {
        Self::from_distribution(name, IMSIsotopeDistribution::new(nominal_mass))
    }

    pub fn from_mass(name: &str, mass: f64) -> Result<Self> {
        Self::from_distribution(name, IMSIsotopeDistribution::from_mass(mass)?)
    }

    pub fn from_distribution(name: &str, isotopes: IMSIsotopeDistribution) -> Result<Self> {
        check_label(name)?;
        Ok(Self {
            name: name.to_owned(),
            sequence: name.to_owned(),
            isotopes,
        })
    }

    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn sequence(&self) -> &str {
        &self.sequence
    }
    pub fn nominal_mass(&self) -> u32 {
        self.isotopes.nominal_mass()
    }
    pub fn isotope_distribution(&self) -> &IMSIsotopeDistribution {
        &self.isotopes
    }

    /// Changes only the name, preserving the independent sequence label.
    pub fn set_name(&mut self, name: &str) -> Result<()> {
        check_label(name)?;
        self.name = name.to_owned();
        Ok(())
    }

    pub fn set_sequence(&mut self, sequence: &str) -> Result<()> {
        check_label(sequence)?;
        self.sequence = sequence.to_owned();
        Ok(())
    }

    pub fn set_isotope_distribution(&mut self, isotopes: IMSIsotopeDistribution) {
        self.isotopes = isotopes;
    }
    pub fn mass(&self, index: usize) -> Result<f64> {
        self.isotopes.mass(index)
    }
    pub fn average_mass(&self) -> Result<f64> {
        self.isotopes.average_mass()
    }

    /// Mass of isotope zero minus the signed number of missing electrons times 0.00054858.
    /// The source default of one electron is supplied explicitly by callers.
    pub fn ion_mass(&self, electrons: i32) -> Result<f64> {
        finite(
            self.mass(0)? - f64::from(electrons) * Self::ELECTRON_MASS_IN_U,
            "IMS ion mass",
        )
    }

    /// Source labels and six-significant-digit isotope lines, with a final blank line.
    pub fn to_text(&self, options: IMSIsotopeOptions) -> Result<String> {
        let mut work = Work::default();
        let labels = self.name.len()
            + self.sequence.len()
            + "name:\t\nsequence:\t\nisotope distribution:\n\n".len();
        work.consume(labels)?;
        work.allocate(labels)?;
        let distribution = self.isotopes.text_with_work(options, &mut work)?;
        let bytes = labels
            .checked_add(distribution.len())
            .ok_or_else(|| invalid("IMS element text overflow"))?;
        if bytes > MAX_TEXT {
            return Err(invalid("IMS element text exceeds 8 MiB bound"));
        }
        // The distribution remains live while its text is copied into the result.
        work.consume(distribution.len())?;
        work.allocate(distribution.len())?;
        let mut text = String::new();
        text.try_reserve_exact(bytes)
            .map_err(|_| invalid("cannot allocate IMS element text"))?;
        text.push_str("name:\t");
        text.push_str(&self.name);
        text.push_str("\nsequence:\t");
        text.push_str(&self.sequence);
        text.push_str("\nisotope distribution:\n");
        text.push_str(&distribution);
        text.push('\n');
        Ok(text)
    }
}

fn check_label(label: &str) -> Result<()> {
    if label.len() > MAX_LABEL {
        return Err(invalid("IMS element label exceeds 1 MiB"));
    }
    Ok(())
}
