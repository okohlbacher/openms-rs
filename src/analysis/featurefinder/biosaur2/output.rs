// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Biosaur2 output: the `FeatureMap` conversion (`convertToFeatureMap_`) and
//! the two TSV exports (`writeTSV`, `writeHills`).

use super::{
    Biosaur2Algorithm, Biosaur2LogLevel, Biosaur2Message, ConvexHullMode, Hill, PeptideFeature, g6,
    isotopes,
};
use crate::concept::constants::user_param;
use crate::concept::unique_id::{HasUniqueId, UniqueIdGenerator};
use crate::identification::ProteinIdentification;
use crate::kernel::{ConvexHull2D, Feature, FeatureMap, Point2D};
use crate::metadata::{ImTypes, MetaValue};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::Path;

/// Convert peptide features to `FeatureMap` entries, the source
/// `convertToFeatureMap_`.
///
/// Each feature takes the monoisotopic m/z, the apex retention time, the apex
/// intensity (narrowed to `f32`), the charge and, as overall quality, the
/// isotope count; the meta values `mass_calib`, `n_isotopes`, `n_scans`,
/// `intensity_sum`, `FAIMS_CV` (unless the drift time is unset) and, for real
/// ion-mobility data, `IM`, `IM_min` and `IM_max`. The convex hulls follow
/// [`ConvexHullMode`]: one hull per contributing hill (in hill-index order) or
/// one bounding box around all their points; a feature whose hills cannot be
/// resolved gets a two-point hull along its retention-time span. Each feature,
/// then the map, draws a unique id from `generator`, and the map carries one
/// empty protein identification. The isotope assignments are checked with
/// [`isotopes::debug_check_isotope_consistency`] first.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when a hull point or a meta value is not
/// finite, which finite input cannot produce.
pub fn convert_to_feature_map(
    algorithm: &Biosaur2Algorithm,
    features: &[PeptideFeature],
    hills: &[Hill],
    generator: &mut UniqueIdGenerator,
    messages: &mut Vec<Biosaur2Message>,
) -> Result<FeatureMap> {
    let mut feature_map = FeatureMap::new();
    let use_mass_trace_hulls = algorithm.convex_hulls == ConvexHullMode::MassTraces;
    let lookup: BTreeMap<usize, &Hill> = hills.iter().map(|hill| (hill.hill_idx, hill)).collect();

    for f in features {
        if let Some(mono) = lookup.get(&f.mono_hill_idx) {
            for iso in &f.isotopes {
                if let Some(iso_hill) = lookup.get(&iso.hill_idx) {
                    isotopes::debug_check_isotope_consistency(
                        algorithm,
                        "convertToFeatureMap_",
                        mono.mz_weighted_mean,
                        mono.rt_apex,
                        mono.hill_idx,
                        f.charge,
                        algorithm.itol,
                        iso_hill,
                        iso.isotope_number,
                        messages,
                    );
                }
            }
        }

        let mut feature = Feature::new(f.rt_apex, f.mz, f.intensity_apex as f32);
        feature.charge = f.charge;
        feature.quality = f.n_isotopes as f32;

        let mut pattern_hill_ids = BTreeSet::new();
        pattern_hill_ids.insert(f.mono_hill_idx);
        pattern_hill_ids.extend(f.isotopes.iter().map(|iso| iso.hill_idx));

        // (min rt, min mz, max rt, max mz)
        let mut feature_box: Option<(f64, f64, f64, f64)> = None;
        let mut im_min = f64::MAX;
        let mut im_max = f64::MIN;
        let mut have_im_values = false;
        for hill_id in &pattern_hill_ids {
            let Some(hill) = lookup.get(hill_id) else {
                continue;
            };
            if hill.rt_values.is_empty() || hill.mz_values.is_empty() {
                continue;
            }
            for &im in &hill.ion_mobilities {
                if im < 0.0 {
                    continue;
                }
                have_im_values = true;
                im_min = if im < im_min { im } else { im_min };
                im_max = if im_max < im { im } else { im_max };
            }
            let n_pts = hill.rt_values.len().min(hill.mz_values.len());
            let points: Vec<Point2D> = hill.rt_values[..n_pts]
                .iter()
                .zip(&hill.mz_values[..n_pts])
                .map(|(&rt, &mz)| Point2D::new(rt, mz))
                .collect();
            if use_mass_trace_hulls {
                feature
                    .convex_hulls
                    .push(ConvexHull2D::from_points(&points)?);
            } else {
                for point in &points {
                    feature_box = Some(match feature_box {
                        None => (point.rt, point.mz, point.rt, point.mz),
                        Some((min_rt, min_mz, max_rt, max_mz)) => (
                            min_rt.min(point.rt),
                            min_mz.min(point.mz),
                            max_rt.max(point.rt),
                            max_mz.max(point.mz),
                        ),
                    });
                }
            }
        }
        if let (false, Some((min_rt, min_mz, max_rt, max_mz))) = (use_mass_trace_hulls, feature_box)
        {
            feature.convex_hulls.push(ConvexHull2D::from_points(&[
                Point2D::new(min_rt, min_mz),
                Point2D::new(max_rt, min_mz),
                Point2D::new(min_rt, max_mz),
                Point2D::new(max_rt, max_mz),
            ])?);
        }
        if feature.convex_hulls.is_empty() {
            feature.convex_hulls.push(ConvexHull2D::from_points(&[
                Point2D::new(f.rt_start, f.mz),
                Point2D::new(f.rt_end, f.mz),
            ])?);
        }

        let metadata = &mut feature.metadata;
        metadata.insert("mass_calib".into(), MetaValue::try_from(f.mass_calib)?);
        metadata.insert(
            "n_isotopes".into(),
            MetaValue::try_from(f.n_isotopes as u64)?,
        );
        metadata.insert("n_scans".into(), MetaValue::try_from(f.n_scans as u64)?);
        metadata.insert(
            "intensity_sum".into(),
            MetaValue::try_from(f.intensity_sum)?,
        );
        if f.drift_time != ImTypes::DRIFTTIME_NOT_SET {
            metadata.insert("FAIMS_CV".into(), MetaValue::try_from(f.drift_time)?);
        }
        if f.ion_mobility >= 0.0 && have_im_values {
            metadata.insert(user_param::IM.into(), MetaValue::try_from(f.ion_mobility)?);
            metadata.insert("IM_min".into(), MetaValue::try_from(im_min)?);
            metadata.insert("IM_max".into(), MetaValue::try_from(im_max)?);
        }
        feature.unique_id.ensure_unique_id(generator);
        feature_map.features.push(feature);
    }
    finish_feature_map(&mut feature_map, generator);
    Ok(feature_map)
}

/// The source's closing `applyMemberFunction(&UniqueIdInterface::ensureUniqueId)`,
/// `ensureUniqueId()` and `getProteinIdentifications().resize(1)`: the map,
/// then every feature without one, draws a unique id, and the map is left with
/// exactly one protein identification.
pub(crate) fn finish_feature_map(feature_map: &mut FeatureMap, generator: &mut UniqueIdGenerator) {
    feature_map.unique_id.ensure_unique_id(generator);
    let mut pending: Vec<&mut Feature> = feature_map.features.iter_mut().collect();
    while let Some(feature) = pending.pop() {
        feature.unique_id.ensure_unique_id(generator);
        pending.extend(feature.subordinates.iter_mut().rev());
    }
    feature_map
        .protein_identifications
        .resize_with(1, ProteinIdentification::default);
}

/// The Biosaur2 feature TSV, as [`write_tsv`] writes it.
///
/// A header `massCalib rtApex intensityApex intensitySum charge nIsotopes nScans
/// mz rtStart rtEnd FAIMS IM` (tab-separated) and one line per feature, the
/// doubles as `std::ostream` prints them at precision 6. The format matches
/// the output of the reference Python implementation.
pub fn tsv_text(features: &[PeptideFeature]) -> String {
    let mut out = String::from(
        "massCalib\trtApex\tintensityApex\tintensitySum\tcharge\tnIsotopes\tnScans\tmz\trtStart\trtEnd\tFAIMS\tIM\n",
    );
    for f in features {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            g6(f.mass_calib),
            g6(f.rt_apex),
            g6(f.intensity_apex),
            g6(f.intensity_sum),
            f.charge,
            f.n_isotopes,
            f.n_scans,
            g6(f.mz),
            g6(f.rt_start),
            g6(f.rt_end),
            g6(f.drift_time),
            g6(f.ion_mobility)
        ));
    }
    out
}

/// The hill diagnostics TSV, as [`write_hills`] writes it.
///
/// A header `hill_idx mz rtStart rtEnd rtApex intensityApex intensitySum nScans
/// IM_min IM_max` and one line per hill, the ion-mobility range over the
/// hill's non-negative values (-1 and -1 when it has none).
pub fn hills_text(hills: &[Hill]) -> String {
    let mut out = String::from(
        "hill_idx\tmz\trtStart\trtEnd\trtApex\tintensityApex\tintensitySum\tnScans\tIM_min\tIM_max\n",
    );
    for hill in hills {
        let mut im_min = f64::MAX;
        let mut im_max = f64::MIN;
        let mut have_im = false;
        for &im in &hill.ion_mobilities {
            if im < 0.0 {
                continue;
            }
            have_im = true;
            im_min = if im < im_min { im } else { im_min };
            im_max = if im_max < im { im } else { im_max };
        }
        if !have_im {
            im_min = -1.0;
            im_max = -1.0;
        }
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            hill.hill_idx,
            g6(hill.mz_weighted_mean),
            g6(hill.rt_start),
            g6(hill.rt_end),
            g6(hill.rt_apex),
            g6(hill.intensity_apex),
            g6(hill.intensity_sum),
            hill.length,
            g6(im_min),
            g6(im_max)
        ));
    }
    out
}

fn write_file(path: &Path, text: &str) -> Result<()> {
    let mut file = std::fs::File::create(path).map_err(Error::Io)?;
    file.write_all(text.as_bytes()).map_err(Error::Io)?;
    file.flush().map_err(Error::Io)
}

/// Export peptide features to a Biosaur2-compatible TSV file, the source
/// `writeTSV`: *the TSV file format matches the output of the reference Python
/// implementation, allowing for easy comparison and downstream processing
/// with Biosaur2-compatible tools* ([`tsv_text`]).
///
/// Returns the info line the source logs, `Wrote N features to TSV file:
/// <filename>`.
///
/// # Errors
///
/// Returns [`Error::Io`] when the file cannot be created or written (the
/// source's `Exception::UnableToCreateFile`).
pub fn write_tsv(features: &[PeptideFeature], filename: &Path) -> Result<Biosaur2Message> {
    write_file(filename, &tsv_text(features))?;
    Ok(Biosaur2Message {
        level: Biosaur2LogLevel::Info,
        text: format!(
            "Wrote {} features to TSV file: {}",
            features.len(),
            filename.display()
        ),
    })
}

/// Export the detected hills as TSV for diagnostic purposes, the source
/// `writeHills` ([`hills_text`]).
///
/// Returns the info line the source logs, `Wrote N hills to: <filename>`.
///
/// # Errors
///
/// As [`write_tsv`].
pub fn write_hills(hills: &[Hill], filename: &Path) -> Result<Biosaur2Message> {
    write_file(filename, &hills_text(hills))?;
    Ok(Biosaur2Message {
        level: Biosaur2LogLevel::Info,
        text: format!("Wrote {} hills to: {}", hills.len(), filename.display()),
    })
}
