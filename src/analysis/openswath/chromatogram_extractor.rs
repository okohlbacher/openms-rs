// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Chromatogram extraction for targeted experiments: coordinates from a
//! transition list, extraction, and conversion back to annotated
//! chromatograms.
//!
//! Ports `ANALYSIS/OPENSWATH/ChromatogramExtractor.h` with
//! `ChromatogramExtractor.cpp`. The class extracts chromatographic traces
//! (intensity against retention time) at given m/z, RT and optionally ion
//! mobility coordinates, for MS1 (precursor) and MS2 (fragment, SWATH/DIA)
//! extraction. [`prepare_coordinates`] derives the coordinates from a
//! [`LightTargetedExperiment`]; callers holding a heavyweight
//! `TargetedExperiment` convert it first, as the source's
//! `OpenSwathDataAccessHelper::convertTargetedExp` does.
//! [`ChromatogramExtractor::extract_chromatograms`] runs
//! [`ChromatogramExtractorAlgorithm`], and [`return_chromatogram`] turns the
//! light results into annotated [`MSChromatogram`]s. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXTRACTION_SUPPORT.md).
//!
//! Not ported: the `return_chromatogram` instantiation over the heavyweight
//! `ANALYSIS/TARGETED/TargetedExperiment` (not yet ported), and the private
//! members `prepareSpectra_`, `outsideExtractionWindow_`, `getFilterNr_`,
//! `populatePeptideRTMap_` and `PeptideRTMap_`, which no public member of the
//! pinned source calls.

use super::chromatogram_extractor_algorithm::ChromatogramExtractorAlgorithm;
pub use super::chromatogram_extractor_algorithm::ExtractionCoordinates;
use super::data_access_helper::convert_to_openms_chromatogram;
use super::data_structures::Chromatogram;
use super::openswath_helper::{compute_precursor_id, compute_transition_group_id};
use super::spectrum_access::SpectrumAccess;
use super::transition_experiment::{LightCompound, LightTargetedExperiment, LightTransition};
use crate::constants::C13C12_MASSDIFF_U;
use crate::kernel::{MSChromatogram, MSSpectrum, Precursor};
use crate::metadata::{ChromatogramType, MetaValue, Product};
use crate::{Error, Result};
use std::collections::HashMap;
use std::sync::Arc;

/// Largest number of coordinates [`prepare_coordinates`] may produce.
///
/// Native ceiling; the source has none. It is checked before anything is
/// allocated.
pub const MAX_COORDINATES: usize = 1 << 28;

/// The source's `ChromatogramExtractor`, a thin front end of
/// [`ChromatogramExtractorAlgorithm`] with no state of its own.
///
/// The source class is also a `ProgressLogger`, but none of its members
/// report progress; use
/// [`ChromatogramExtractorAlgorithm::extract_chromatograms_with_progress`]
/// for that.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ChromatogramExtractor;

impl ChromatogramExtractor {
    /// A new extractor. Source default constructor.
    pub fn new() -> Self {
        Self
    }

    /// Extract chromatograms at the m/z and RT of each coordinate, without ion
    /// mobility.
    ///
    /// Source `extractChromatograms` (six arguments), which runs
    /// [`ChromatogramExtractorAlgorithm::extract_chromatograms`] with an ion
    /// mobility window of -1. `output` must be as long as
    /// `extraction_coordinates` (use [`prepare_coordinates`]).
    /// `mz_extraction_window` is the full window in Th, or ppm when `ppm`
    /// (50 ppm extracts 25 ppm on either side); `filter` is `"tophat"`
    /// (`"bartlett"` is declared but not implemented).
    ///
    /// # Errors
    ///
    /// As [`ChromatogramExtractorAlgorithm::extract_chromatograms`].
    #[allow(clippy::too_many_arguments)]
    pub fn extract_chromatograms(
        &self,
        input: &dyn SpectrumAccess,
        output: &mut [Chromatogram],
        extraction_coordinates: &[ExtractionCoordinates],
        mz_extraction_window: f64,
        ppm: bool,
        filter: &str,
    ) -> Result<()> {
        ChromatogramExtractorAlgorithm::new().extract_chromatograms(
            input,
            output,
            extraction_coordinates,
            mz_extraction_window,
            ppm,
            -1.0,
            filter,
        )
    }

    /// Extract chromatograms at the m/z, ion mobility and RT of each
    /// coordinate.
    ///
    /// Source `extractChromatograms` (seven arguments), a direct call of
    /// [`ChromatogramExtractorAlgorithm::extract_chromatograms`];
    /// `im_extraction_window` is the full ion-mobility window. The source
    /// notes that the algorithm class should be used directly whenever
    /// possible.
    ///
    /// # Errors
    ///
    /// As [`ChromatogramExtractorAlgorithm::extract_chromatograms`].
    #[allow(clippy::too_many_arguments)]
    pub fn extract_chromatograms_im(
        &self,
        input: &dyn SpectrumAccess,
        output: &mut [Chromatogram],
        extraction_coordinates: &[ExtractionCoordinates],
        mz_extraction_window: f64,
        ppm: bool,
        im_extraction_window: f64,
        filter: &str,
    ) -> Result<()> {
        ChromatogramExtractorAlgorithm::new().extract_chromatograms(
            input,
            output,
            extraction_coordinates,
            mz_extraction_window,
            ppm,
            im_extraction_window,
            filter,
        )
    }
}

/// Apply the source's RT-window rule to `coord` for compound `pep`.
fn apply_rt_window(
    coord: &mut ExtractionCoordinates,
    pep: &LightCompound,
    rt_extraction_window: f64,
) -> Result<()> {
    if rt_extraction_window >= 0.0 {
        if pep.rt.is_nan() {
            return Err(Error::InvalidValue(format!(
                "Error: Peptide {} does not have retention time information which is necessary to perform an RT-limited extraction",
                pep.id
            )));
        }
        coord.rt_start = pep.rt - rt_extraction_window / 2.0;
        coord.rt_end = pep.rt + rt_extraction_window / 2.0;
    } else if rt_extraction_window.is_nan() {
        if pep.rt_start.is_nan() || pep.rt_end.is_nan() {
            return Err(Error::InvalidValue(format!(
                "Error: rt_extraction_window=NaN requires the compound '{}' to carry a retention-time range (rt_start, rt_end), but none is set.",
                pep.id
            )));
        }
        coord.rt_start = pep.rt_start;
        coord.rt_end = pep.rt_end;
    }
    Ok(())
}

/// Build the extraction coordinates (and matching empty chromatograms) of a
/// targeted experiment, sorted by m/z.
///
/// Source static `prepare_coordinates`. For MS2 (`ms1 == false`) there is one
/// coordinate per transition: the product m/z, the precursor m/z as
/// `mz_precursor`, and the transition's native id. For MS1 there is one per
/// compound: the precursor m/z of the compound's first transition and the id
/// [`compute_precursor_id`]`(compound id, 0)`, plus `ms1_isotopes` further
/// coordinates at `+k * C13C12_MASSDIFF_U` with ids `..._Precursor_i<k>`. A
/// compound without transitions still gets one coordinate, with m/z 0, the
/// whole RT range and ion mobility 0, and no isotopes (the source logs a
/// warning for it).
///
/// The RT range of each coordinate follows `rt_extraction_window`:
///
/// - non-negative: the full window, centred on the compound's `rt`;
/// - negative: the whole RT range (`rt_start = 0`, `rt_end = -1`);
/// - NaN: the compound's own `rt_start` and `rt_end`, which must be set (as
///   `convertTargetedExp` sets them from a compound with exactly two RTs).
///
/// Each coordinate's `ion_mobility` is the compound's drift time.
///
/// One empty chromatogram per coordinate is appended to
/// `output_chromatograms`, and the coordinates to `coordinates`; the whole of
/// `coordinates`, including any entries it held before, is then stably
/// sorted by m/z, as in the source. A non-positive `ms1_isotopes` adds no
/// isotopes.
///
/// # Errors
///
/// All errors leave both vectors unchanged; the source may have appended some
/// entries before it throws.
///
/// - [`Error::InvalidValue`] when an RT-limited extraction meets a compound
///   without RT, or a NaN window meets a compound without an RT range (the
///   source's `Exception::IllegalArgument`).
/// - [`Error::MissingInformation`] in MS2 mode when a transition references
///   a compound the experiment does not contain; the source dereferences a
///   null pointer there.
/// - [`Error::InvalidRange`] when more than [`MAX_COORDINATES`] coordinates
///   would be produced.
pub fn prepare_coordinates(
    output_chromatograms: &mut Vec<Chromatogram>,
    coordinates: &mut Vec<ExtractionCoordinates>,
    transition_exp_used: &LightTargetedExperiment,
    rt_extraction_window: f64,
    ms1: bool,
    ms1_isotopes: i32,
) -> Result<()> {
    let transitions = transition_exp_used.transitions();
    let compounds = transition_exp_used.compounds();
    let mut pep2tr: HashMap<&str, Vec<&LightTransition>> = HashMap::new();
    for transition in transitions {
        pep2tr
            .entry(transition.peptide_ref())
            .or_default()
            .push(transition);
    }
    let tr2pep = transition_exp_used.compound_reference_map();

    let itersize = if ms1 {
        compounds.len()
    } else {
        transitions.len()
    };
    let isotopes = if ms1 && ms1_isotopes > 0 {
        usize::try_from(ms1_isotopes).unwrap_or(usize::MAX)
    } else {
        0
    };
    let expected = isotopes
        .checked_add(1)
        .and_then(|per| itersize.checked_mul(per))
        .filter(|&n| n <= MAX_COORDINATES)
        .ok_or_else(|| {
            Error::InvalidRange(format!(
                "more than {MAX_COORDINATES} extraction coordinates requested"
            ))
        })?;

    let mut new_coords: Vec<ExtractionCoordinates> = Vec::with_capacity(expected);
    for i in 0..itersize {
        let mut coord = ExtractionCoordinates {
            rt_start: 0.0,
            rt_end: -1.0,
            ..ExtractionCoordinates::default()
        };
        let pep: &LightCompound = if ms1 {
            let pep = &compounds[i];
            coord.id = compute_precursor_id(&pep.id, 0);
            match pep2tr.get(pep.id.as_str()) {
                // The precursor m/z lives on the transitions only.
                Some(trs) => coord.mz = trs[0].precursor_mz(),
                None => {
                    new_coords.push(coord);
                    continue;
                }
            }
            pep
        } else {
            let transition = &transitions[i];
            coord.mz = transition.product_mz();
            coord.mz_precursor = transition.precursor_mz();
            coord.id = transition.native_id().to_string();
            tr2pep
                .get(transition.peptide_ref())
                .copied()
                .ok_or_else(|| {
                    Error::MissingInformation(format!(
                        "transition '{}' references compound '{}', which the experiment does not contain",
                        transition.native_id(),
                        transition.peptide_ref()
                    ))
                })?
        };
        apply_rt_window(&mut coord, pep, rt_extraction_window)?;
        coord.ion_mobility = pep.drift_time();
        let (pep_id, base_mz) = (pep.id.as_str(), coord.mz);
        let template = coord.clone();
        new_coords.push(coord);
        for k in 1..=isotopes {
            // `k` is at most `i32::MAX`, as `ms1_isotopes` is.
            let k_i32 = i32::try_from(k).unwrap_or(i32::MAX);
            let mut coord_new = template.clone();
            coord_new.id = compute_precursor_id(pep_id, k_i32);
            coord_new.mz = base_mz + f64::from(k_i32) * C13C12_MASSDIFF_U;
            new_coords.push(coord_new);
        }
    }

    output_chromatograms.extend(std::iter::repeat_with(Chromatogram::new).take(new_coords.len()));
    coordinates.extend(new_coords);
    coordinates.sort_by(ExtractionCoordinates::cmp_by_mz);
    Ok(())
}

/// The id (peptide sequence, else compound name) and charge of compound
/// `id`, source `extract_id_<LightTargetedExperiment>`.
fn extract_id(exp: &LightTargetedExperiment, id: &str) -> Result<(String, i32)> {
    let comp = exp.compound_by_ref(id).ok_or_else(|| {
        Error::MissingInformation(format!(
            "compound '{id}' is not part of the targeted experiment"
        ))
    })?;
    // An empty sequence falls through to the compound name, which may itself
    // be empty (iRT calibration peptides); consumers expect an empty
    // peptide_sequence then.
    let name = if comp.sequence.is_empty() {
        comp.compound_name.clone()
    } else {
        comp.sequence.clone()
    };
    Ok((name, comp.charge))
}

/// Convert light chromatograms into annotated [`MSChromatogram`]s and append
/// them to `output_chromatograms`.
///
/// Source template `return_chromatogram<LightTargetedExperiment>`. For
/// chromatogram `i` and coordinate `coordinates[i]` it sets
///
/// 1. the peaks, from the light time and intensity arrays;
/// 2. the native id, from the coordinate id;
/// 3. the precursor: in MS1 mode (`ms1`) the coordinate m/z, and a base-peak
///    chromatogram type; in MS2 mode the precursor m/z of the transition with
///    the coordinate's id, the isolation window offsets of the first
///    precursor of `settings` (when it has one), a product at the transition's
///    product m/z, and a selected-reaction-monitoring type;
/// 4. the precursor charge and a `peptide_sequence` meta value, from the
///    compound: in MS1 mode the transition group id of the coordinate
///    ([`compute_transition_group_id`]; skipped when that is empty), in MS2
///    mode the transition's peptide reference, else its compound reference;
///    the value is the compound's sequence, or its name when the sequence is
///    empty;
/// 5. when the coordinate's `ion_mobility` is non-negative and
///    `im_extraction_width` is positive, the precursor drift time and a
///    drift-time window of `im_extraction_width / 2` on either side;
/// 6. the instrument settings, acquisition info and source file of
///    `settings`, and its data-processing entries, each marked with the meta
///    value `performed_on_spectra = "true"`.
///
/// `settings` is the kernel spectrum whose settings the source passes as a
/// `SpectrumSettings` (for example the first spectrum of the SWATH map). The
/// source marks the caller's shared data-processing objects themselves; the
/// port marks copies and leaves `settings` untouched.
///
/// # Errors
///
/// Nothing is appended on error.
///
/// - [`Error::InvalidValue`] when there are fewer coordinates than
///   chromatograms (the source reads past the end), or a light chromatogram's
///   arrays differ in length.
/// - [`Error::MissingInformation`] when an MS2 coordinate id names no
///   transition, or a compound reference names no compound; the source
///   dereferences a null pointer in both cases.
pub fn return_chromatogram(
    chromatograms: &[Chromatogram],
    coordinates: &[ExtractionCoordinates],
    transition_exp_used: &LightTargetedExperiment,
    settings: &MSSpectrum,
    output_chromatograms: &mut Vec<MSChromatogram>,
    ms1: bool,
    im_extraction_width: f64,
) -> Result<()> {
    if coordinates.len() < chromatograms.len() {
        return Err(Error::InvalidValue(format!(
            "{} chromatograms but only {} extraction coordinates",
            chromatograms.len(),
            coordinates.len()
        )));
    }
    let mut trans_map: HashMap<&str, &LightTransition> = HashMap::new();
    for transition in transition_exp_used.transitions() {
        trans_map.insert(transition.native_id(), transition);
    }
    let data_processing: Vec<Arc<crate::metadata::DataProcessing>> = settings
        .data_processing
        .iter()
        .map(|dp| {
            let mut marked = (**dp).clone();
            marked
                .metadata
                .insert("performed_on_spectra".into(), MetaValue::from("true"));
            Arc::new(marked)
        })
        .collect();

    let mut converted = Vec::with_capacity(chromatograms.len());
    for (chromptr, coord) in chromatograms.iter().zip(coordinates) {
        let mut chrom = MSChromatogram::default();
        convert_to_openms_chromatogram(chromptr, &mut chrom)?;
        chrom.native_id = coord.id.clone();

        let mut prec = Precursor::default();
        if ms1 {
            prec.mz = coord.mz;
            chrom.chromatogram_type = ChromatogramType::BasePeak;
            let transition_group_id = compute_transition_group_id(&coord.id);
            if !transition_group_id.is_empty() {
                let (name, charge) = extract_id(transition_exp_used, &transition_group_id)?;
                prec.charge = charge;
                prec.cv_terms
                    .metadata
                    .insert("peptide_sequence".into(), MetaValue::from(name));
            }
        } else {
            let transition = trans_map.get(coord.id.as_str()).copied().ok_or_else(|| {
                Error::MissingInformation(format!(
                    "no transition with native id '{}' in the targeted experiment",
                    coord.id
                ))
            })?;
            prec.mz = transition.precursor_mz();
            if let Some(first) = settings.precursors.first() {
                prec.isolation_window_lower_offset = first.isolation_window_lower_offset;
                prec.isolation_window_upper_offset = first.isolation_window_upper_offset;
            }
            chrom.product = Product {
                mz: transition.product_mz(),
                ..Product::default()
            };
            chrom.chromatogram_type = ChromatogramType::SelectedReactionMonitoring;
            let reference = if transition.peptide_ref().is_empty() {
                transition.compound_ref()
            } else {
                transition.peptide_ref()
            };
            let (name, charge) = extract_id(transition_exp_used, reference)?;
            prec.charge = charge;
            prec.cv_terms
                .metadata
                .insert("peptide_sequence".into(), MetaValue::from(name));
        }

        if coord.ion_mobility >= 0.0 && im_extraction_width > 0.0 {
            prec.drift_time = Some(coord.ion_mobility);
            prec.drift_window_lower_offset = im_extraction_width / 2.0;
            prec.drift_window_upper_offset = im_extraction_width / 2.0;
        }
        chrom.precursor = prec;
        chrom.instrument_settings = settings.instrument_settings.clone();
        chrom.acquisition_info = settings.acquisition_info.clone();
        chrom.source_file = settings.source_file.clone();
        chrom
            .data_processing
            .extend(data_processing.iter().cloned());
        converted.push(chrom);
    }
    output_chromatograms.extend(converted);
    Ok(())
}
