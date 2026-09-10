// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Activation presets and the seven fixed abundant immonium ions from
//! TheoreticalSpectrumGenerator.cpp at OpenMS4-core revision 7c029e8.

use super::{
    AASequence, PeakRow, TheoreticalIonSeries as Ion, TheoreticalSpectrumGenerator, push_peak,
};
use crate::kernel::MSSpectrum;
use crate::metadata::ActivationMethod;
use crate::{Error, Result};
use std::collections::BTreeSet;

impl TheoreticalSpectrumGenerator {
    /// Generate with the source activation preset and otherwise default settings.
    ///
    /// CID uses b/y; HCID and HCD add a. ECD and ETD use c/z./z', while
    /// ETciD and EThcD include all eight series. Other methods are unsupported.
    /// Zero charge is treated as two. Charges up to two produce singly charged
    /// fragments; higher charges add doubly charged fragments. As in the source
    /// factory, the precursor metadata charge is inferred as two or three from
    /// that fragment range, rather than copied from the supplied charge.
    pub fn generate_for_activation(
        method: ActivationMethod,
        peptide: &AASequence,
        precursor_charge: u16,
    ) -> Result<MSSpectrum> {
        let ion_series = match method {
            ActivationMethod::Cid => vec![Ion::B, Ion::Y],
            ActivationMethod::Hcid | ActivationMethod::Hcd => vec![Ion::A, Ion::B, Ion::Y],
            ActivationMethod::Ecd | ActivationMethod::Etd => {
                vec![Ion::C, Ion::ZPlusOne, Ion::ZPlusTwo]
            }
            ActivationMethod::Etcid | ActivationMethod::Ethcd => vec![
                Ion::A,
                Ion::B,
                Ion::C,
                Ion::X,
                Ion::Y,
                Ion::Z,
                Ion::ZPlusOne,
                Ion::ZPlusTwo,
            ],
            _ => {
                return Err(Error::Unsupported(format!(
                    "theoretical spectrum activation method {} is not supported",
                    method.short_name()
                )));
            }
        };
        // Replacing zero by two leaves this source charge-range branch unchanged.
        let max_charge = if precursor_charge <= 2 { 1 } else { 2 };
        Self {
            ion_series,
            ..Self::default()
        }
        .generate(peptide, 1, max_charge, None)
    }

    pub(super) fn emit_immonium(
        &self,
        rows: &mut Vec<PeakRow>,
        peptide: &AASequence,
    ) -> Result<()> {
        let mut unmodified = BTreeSet::new();
        for (index, residue) in peptide.as_str().bytes().enumerate() {
            if peptide.residue_modification(index)?.is_none() {
                unmodified.insert(residue);
            }
        }
        // The C++ has(Residue) predicate requires an unmodified residue value.
        // It checks L only despite the shared isoleucine/leucine annotation.
        for (residue, mz, name) in [
            (b'P', 70.0656, "iP+"),
            (b'C', 76.0221, "iC+"),
            (b'L', 86.09698, "iL/I+"),
            (b'H', 110.0718, "iH+"),
            (b'F', 120.0813, "iF+"),
            (b'Y', 136.0762, "iY+"),
            (b'W', 159.0922, "iW+"),
        ] {
            if unmodified.contains(&residue) {
                push_peak(rows, mz, 1.0, name, 1)?;
            }
        }
        Ok(())
    }
}
