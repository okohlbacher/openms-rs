// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{MSChromatogram, MSSpectrum, Precursor, data_array::Meter};
use crate::Result;
use crate::metadata::{AcquisitionInfo, InstrumentSettings, Product, SourceFile};

pub(super) fn instrument(m: &mut Meter<'_>, v: &InstrumentSettings) -> Result<()> {
    m.slots::<InstrumentSettings>(1)?;
    m.meta(&v.metadata)?;
    m.slots::<crate::metadata::ScanWindow>(v.scan_windows.len())?;
    for window in &v.scan_windows {
        m.meta(&window.metadata)?;
    }
    Ok(())
}
pub(super) fn acquisition(m: &mut Meter<'_>, v: &AcquisitionInfo) -> Result<()> {
    m.slots::<AcquisitionInfo>(1)?;
    m.text(&v.method_of_combination)?;
    m.meta(&v.metadata)?;
    m.slots::<crate::metadata::Acquisition>(v.acquisitions.len())?;
    for a in &v.acquisitions {
        m.text(&a.identifier)?;
        m.meta(&a.metadata)?;
    }
    Ok(())
}
pub(crate) fn source(m: &mut Meter<'_>, v: &SourceFile) -> Result<()> {
    m.slots::<SourceFile>(1)?;
    for text in [
        &v.name,
        &v.path,
        &v.file_type,
        &v.checksum,
        &v.native_id_type,
        &v.native_id_type_accession,
    ] {
        m.text(text)?;
    }
    m.cv(&v.cv_terms)
}
pub(super) fn precursor(m: &mut Meter<'_>, v: &Precursor) -> Result<()> {
    m.slots::<Precursor>(1)?;
    m.tree::<crate::metadata::ActivationMethod>(v.activation_methods.len())?;
    m.slots::<i32>(v.possible_charge_states.len())?;
    if let Some(text) = &v.spectrum_reference {
        m.text(text)?;
    }
    m.cv(&v.cv_terms)
}
pub(super) fn product(m: &mut Meter<'_>, v: &Product) -> Result<()> {
    m.slots::<Product>(1)?;
    m.cv(&v.cv_terms)
}
fn nondefault(i: &InstrumentSettings, a: &AcquisitionInfo, s: &SourceFile) -> bool {
    i.scan_mode != crate::metadata::ScanMode::Unknown
        || i.zoom_scan
        || i.polarity != crate::metadata::Polarity::Unknown
        || !i.scan_windows.is_empty()
        || !i.metadata.is_empty()
        || !a.acquisitions.is_empty()
        || !a.method_of_combination.is_empty()
        || !a.metadata.is_empty()
        || !s.name.is_empty()
        || !s.path.is_empty()
        || s.size_mb.to_bits() != 0
        || !s.file_type.is_empty()
        || !s.checksum.is_empty()
        || s.checksum_type != crate::metadata::ChecksumType::Unknown
        || !s.native_id_type.is_empty()
        || !s.native_id_type_accession.is_empty()
        || !s.cv_terms.terms().is_empty()
        || !s.cv_terms.metadata.is_empty()
}
macro_rules! settings {
    ($ty:ty) => {
        impl $ty {
            /// Charge newly attached owned acquisition data before a full-record
            /// clone. Processing Arc payloads are shared; only handle slots copy.
            pub(crate) fn acquisition_with_budget(
                &self,
                work: &mut usize,
                bytes: &mut usize,
            ) -> Result<()> {
                let mut m = Meter { work, bytes };
                instrument(&mut m, &self.instrument_settings)?;
                acquisition(&mut m, &self.acquisition_info)?;
                source(&mut m, &self.source_file)?;
                m.slots::<std::sync::Arc<crate::metadata::DataProcessing>>(
                    self.data_processing.len(),
                )?;
                self.products_with_budget(&mut m)
            }
            pub(crate) fn validate_acquisition_settings(&self) -> Result<()> {
                // Validation visits shared processing payload; unlike cloning,
                // this pass must meter the referenced records themselves.
                let (mut work, mut bytes) = (50_000_000, 256 * 1024 * 1024);
                self.acquisition_with_budget(&mut work, &mut bytes)?;
                let mut meter = Meter {
                    work: &mut work,
                    bytes: &mut bytes,
                };
                for p in &self.data_processing {
                    meter.slots::<crate::metadata::DataProcessing>(1)?;
                    meter.text(&p.software.name)?;
                    meter.text(&p.software.version)?;
                    meter.cv(&p.software.cv_terms)?;
                    meter.tree::<crate::metadata::ProcessingAction>(p.actions.len())?;
                    meter.meta(&p.metadata)?;
                }
                self.instrument_settings.validate()?;
                self.acquisition_info.validate()?;
                self.source_file.validate()?;
                self.validate_products()?;
                for p in &self.data_processing {
                    p.validate()?;
                }
                Ok(())
            }
            /// True if a legacy record transport would omit attached settings.
            /// This check examines only scalar fields and container lengths.
            pub fn has_acquisition_settings(&self) -> bool {
                nondefault(
                    &self.instrument_settings,
                    &self.acquisition_info,
                    &self.source_file,
                ) || !self.data_processing.is_empty()
                    || self.has_products_or_type()
            }
        }
    };
}
settings!(MSSpectrum);
settings!(MSChromatogram);
impl MSSpectrum {
    fn products_with_budget(&self, m: &mut Meter<'_>) -> Result<()> {
        m.slots::<Product>(self.products.len())?;
        for p in &self.products {
            product(m, p)?;
        }
        Ok(())
    }
    fn validate_products(&self) -> Result<()> {
        for p in &self.products {
            p.validate()?;
        }
        Ok(())
    }
    fn has_products_or_type(&self) -> bool {
        !self.products.is_empty()
    }
}
impl MSChromatogram {
    fn products_with_budget(&self, _: &mut Meter<'_>) -> Result<()> {
        Ok(())
    }
    fn validate_products(&self) -> Result<()> {
        Ok(())
    }
    fn has_products_or_type(&self) -> bool {
        self.chromatogram_type != crate::metadata::ChromatogramType::Mass
    }
}
