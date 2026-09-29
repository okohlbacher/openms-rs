// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Conversions between OpenMS kernel records and the OpenSWATH light records.
//!
//! Ports the spectrum and chromatogram half of
//! `ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h` (class
//! `OpenSwathDataAccessHelper`, whose static members become free functions).
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).
//!
//! Not ported yet, because they read the heavyweight
//! `ANALYSIS/TARGETED/TargetedExperiment` and `TargetedExperimentHelper`, which
//! the crate does not have: `convertTargetedExp`, both `convertTargetedCompound`
//! overloads and `convertPeptideToAASequence`.
//!
//! # Lengths
//!
//! The source walks the m/z (time) array and reads the intensity array at the
//! same index without checking its length. The port checks, and returns
//! [`Error::InvalidValue`](crate::Error::InvalidValue) before modifying the destination.

use super::data_structures::{Chromatogram, Spectrum};
use crate::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};
use crate::{Error, Result};

fn check_pair(position: usize, intensity: usize, what: &str) -> Result<()> {
    if intensity != position {
        return Err(Error::InvalidValue(format!(
            "the {what} array has {position} entries but the intensity array {intensity}"
        )));
    }
    Ok(())
}

/// Convert an OpenSWATH spectrum into the peaks of an OpenMS spectrum.
///
/// The peaks of `spectrum` are replaced by one peak per m/z value; its
/// metadata is kept. As the kernel's `clear(false)` does, the data arrays of
/// `spectrum` are cleared too, since they would no longer align with the new
/// peaks; further arrays of `sptr` (such as ion mobility) are not carried over,
/// as in the source. Intensities are narrowed to `f32`, the kernel's intensity
/// type, as in the source.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the intensity array is not as long as the m/z
/// array; `spectrum` is then unchanged.
pub fn convert_to_openms_spectrum(sptr: &Spectrum, spectrum: &mut MSSpectrum) -> Result<()> {
    let mz = &sptr.mz_array().data;
    let intensity = &sptr.intensity_array().data;
    check_pair(mz.len(), intensity.len(), "m/z")?;
    spectrum.clear(false);
    spectrum.peaks.reserve(mz.len());
    spectrum.peaks.extend(
        mz.iter()
            .zip(intensity)
            .map(|(&m, &i)| Peak1D::new(m, i as f32)),
    );
    Ok(())
}

/// Convert an OpenMS spectrum into an OpenSWATH spectrum holding only its m/z
/// and intensity arrays. Source `convertToSpectrumPtr`.
///
/// Data arrays are **not** carried over; the spectrum access view uses
/// [`spectrum_with_data_arrays`](super::spectrum_access_openms::spectrum_with_data_arrays),
/// which does carry them.
pub fn convert_to_spectrum(spectrum: &MSSpectrum) -> Spectrum {
    Spectrum::from_arrays(
        spectrum.peaks.iter().map(|p| p.mz).collect(),
        spectrum
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect(),
    )
}

/// Convert an OpenMS chromatogram into an OpenSWATH chromatogram holding only
/// its time and intensity arrays. Source `convertToChromatogramPtr`.
pub fn convert_to_chromatogram(chromatogram: &MSChromatogram) -> Chromatogram {
    Chromatogram::from_arrays(
        chromatogram.peaks.iter().map(|p| p.rt).collect(),
        chromatogram
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect(),
    )
}

/// Convert an OpenSWATH chromatogram into the peaks of an OpenMS chromatogram.
///
/// The rules of [`convert_to_openms_spectrum`] apply: peaks and data arrays are
/// replaced, metadata kept, intensities narrowed to `f32`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the intensity array is not as long as the time
/// array; `chromatogram` is then unchanged.
pub fn convert_to_openms_chromatogram(
    cptr: &Chromatogram,
    chromatogram: &mut MSChromatogram,
) -> Result<()> {
    convert_to_openms_chromatogram_filter(chromatogram, cptr, f64::NEG_INFINITY, f64::INFINITY)
        .map(|_| ())
}

/// Convert the points of an OpenSWATH chromatogram with `rt_min <= RT <=
/// rt_max` into the peaks of an OpenMS chromatogram.
///
/// Rules as [`convert_to_openms_chromatogram`]; points outside the window are
/// skipped. A NaN RT is never outside the window by the source's test
/// (`rt < rt_min || rt > rt_max` is false), so it is kept, as in the source.
/// Returns the number of points kept.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the intensity array is not as long as the time
/// array; `chromatogram` is then unchanged.
pub fn convert_to_openms_chromatogram_filter(
    chromatogram: &mut MSChromatogram,
    cptr: &Chromatogram,
    rt_min: f64,
    rt_max: f64,
) -> Result<usize> {
    let rt = &cptr.time_array().data;
    let intensity = &cptr.intensity_array().data;
    check_pair(rt.len(), intensity.len(), "time")?;
    chromatogram.clear(false);
    chromatogram.peaks.reserve(rt.len());
    for (&t, &i) in rt.iter().zip(intensity) {
        if t < rt_min || t > rt_max {
            continue;
        }
        chromatogram.peaks.push(ChromatogramPeak::new(t, i as f32));
    }
    Ok(chromatogram.peaks.len())
}
