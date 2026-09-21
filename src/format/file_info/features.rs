// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The FileInfo summary of featureXML feature maps
//! (`FORMAT/FileInfo.cpp:1077-1145`, `:1978-1984`, `:2097-2100`, `:2201-2256`).
//!
//! The featureXML branch of the report. The map is loaded through
//! [`crate::format::FileHandler::load_feature_map_with_options`] with convex
//! hulls and subordinate features off, as the source sets
//! `setLoadConvexHull(false)` and `setLoadSubordinates(false)` before
//! `loadFeatures`. The branch then writes, in the source order:
//!
//! 1. the feature count;
//! 2. one range block over retention time, m/z and intensity, the source
//!    `FeatureMap::updateRanges`: every feature position and intensity, then
//!    the bounding box of every feature's hulls (none are loaded here). A
//!    feature map has no mobility dimension;
//! 3. the total ion current, the feature intensities summed as `double` in file
//!    order and printed at the stream precision;
//! 4. the charge distribution;
//! 5. the number of features per peptide-identification count, and the
//!    assigned and unassigned identification counts.
//!
//! `-m` adds the document identifier; `-p` the map's data processing; `-s` five
//! statistics blocks at precision [`crate::format::file_info::text_format::WRITTEN_DIGITS_F32`]:
//! intensity, FWHM in RT (the feature width), overall quality, RT quality and
//! m/z quality, each promoted from `float`.
//!
//! The structured [`crate::format::file_info::model::FeatureInfo`] and ranges
//! are filled alongside, as the source fills its `Result`.
//!
//! See `docs/FILE_INFO_SUPPORT.md` for the evidence and the native
//! differences.

#![cfg(feature = "featurexml")]

use super::model::{FeatureInfo, FileInfoResult, Options, Range, RangeSet, Ranges};
use super::report::{
    ReportStream, statistics_buffer, summarize, write_charge_distribution, write_meta_title,
    write_processing, write_processing_title, write_ranges_text, write_ranges_tsv,
    write_statistics_title, write_summary_text,
};
use super::text_format::WRITTEN_DIGITS_F32;
use crate::format::FileHandler;
use crate::format::featurexml::FeatureFileOptions;
use crate::kernel::ranges::RangeBase;
use crate::kernel::{Feature, FeatureMap};
use crate::math::statistic_functions::SummaryStatistics;
use crate::math::x86_64;
use crate::metadata::DataProcessing;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

/// Load the feature map and write the featureXML branch, `-m`, `-p` and `-s`.
pub(crate) fn report(
    path: &Path,
    options: &Options,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    result: &mut FileInfoResult,
) -> Result<()> {
    let load = FeatureFileOptions {
        load_convex_hulls: false,
        load_subordinates: false,
        ..FeatureFileOptions::default()
    };
    let map = FileHandler::load_feature_map_with_options(path, &load)?;
    let ranges = feature_map_ranges(&map)?;

    os.text("Number of features: ")
        .value(map.len())
        .text("\n\n");
    os_tsv
        .text("general: number of features\t")
        .value(map.len())
        .text("\n");
    os.text("Ranges:\n");
    write_ranges_text(os, &ranges, false);
    write_ranges_tsv(os_tsv, "general: ranges: ", &ranges);

    let mut charges: BTreeMap<i32, u64> = BTreeMap::new();
    let mut ids_per_feature: BTreeMap<u64, u64> = BTreeMap::new();
    let mut tic = 0.0_f64;
    let mut assigned = 0_u64;
    for feature in &map.features {
        *charges.entry(feature.charge).or_insert(0) += 1;
        // `FORMAT/FileInfo.cpp:1098,1103`: `double tic = 0.0; tic += feat[i].getIntensity()`,
        // where `getIntensity()` is a `float` — a `cvtss2sd` and an `addsd`.
        // Both can produce a NaN, and plain Rust would leave its bits to the
        // host, which the report's spelling now reads.
        //
        // Not reachable through `FileInfo` today: featureXML does read `inf`,
        // `-inf` and `NaN` for an intensity (`format::featurexml`), but
        // `FeatureMap::ranges` folds intensity into its ranges and refuses a
        // non-finite one, so such a document is declined before this loop runs
        // (`tests/featurexml.rs::a_nonfinite_map_reads_and_its_ranges_are_a_checked_error`,
        // which also pins what the Release build prints for one: `nan`, its
        // sign bit clear, because the file's first NaN intensity propagates
        // quieted). Routed through the emulation regardless, so that the value
        // is the Release build's if the kernel's finite invariant is ever
        // relaxed — that is a kernel change with its own evidence, and this
        // line should not be the thing that then has to be found.
        tic = x86_64::add(tic, x86_64::widen(feature.intensity));
        let ids = to_count(feature.peptide_identifications.len())?;
        *ids_per_feature.entry(ids).or_insert(0) += 1;
        assigned = assigned.checked_add(ids).ok_or_else(count_overflow)?;
    }
    let unassigned = to_count(map.unassigned_peptide_identifications.len())?;

    os.text("Total ion current in features: ")
        .double(tic)
        .text("\n");
    os_tsv
        .text("general: total ion current in features\t")
        .double(tic)
        .text("\n");
    os.text("\n");

    write_charge_distribution(os, os_tsv, "Charge", &charges);

    os.text("Distribution of peptide identifications (IDs) per feature:\n");
    for (ids, count) in &ids_per_feature {
        os.text("  ")
            .value(ids)
            .text(" IDs: ")
            .value(count)
            .text("\n");
        os_tsv
            .text("general: distribution of peptide identifications (IDs) per feature: IDs: ")
            .value(ids)
            .text("\t")
            .value(count)
            .text("\n");
    }
    os.text("\nAssigned peptide identifications: ")
        .value(assigned)
        .text("\n");
    os_tsv
        .text("general: assigned peptide identifications\t")
        .value(assigned)
        .text("\n");
    os.text("Unassigned peptide identifications: ")
        .value(unassigned)
        .text("\n");
    os_tsv
        .text("general: unassigned peptide identifications\t")
        .value(unassigned)
        .text("\n");

    result.feature = Some(FeatureInfo {
        is_consensus: false,
        num_features: to_count(map.len())?,
        tic,
        charges,
        ids_per_element: ids_per_feature,
        assigned_ids: assigned,
        unassigned_ids: unassigned,
        ..FeatureInfo::default()
    });
    result.ranges = Ranges {
        combined: ranges,
        is_experiment: false,
        ..Ranges::default()
    };

    if options.meta {
        write_meta_title(os);
        os.text("Document ID: ").text(&map.identifier).text("\n\n");
        os_tsv
            .text("meta: document ID\t")
            .text(&map.identifier)
            .text("\n");
    }
    if options.processing {
        write_processing_title(os);
        let processing: Vec<&DataProcessing> = map.data_processing.iter().collect();
        write_processing(os, os_tsv, &processing, result);
    }
    if options.statistics {
        write_statistics(&map, os, os_tsv)?;
    }
    Ok(())
}

/// The ranges of `FeatureMap::updateRanges` (`FeatureMap.cpp:260-290`): every
/// feature's RT, m/z and intensity, then the RT and m/z bounds of every
/// feature's hull bounding box, each extended with the source's keep-first rule
/// on equal endpoints ([`RangeBase::extend_value`]).
///
/// [`FeatureMap::ranges`] computes the same bounds but first validates the
/// whole map, which the source does not, and merges with `f64::min` and
/// `f64::max`, whose result for `-0.0` against `0.0` is unspecified. The hull
/// bounding box of one feature is still the kernel's
/// [`Feature::hull_bounding_box`], with that caveat; FileInfo loads no hulls, so
/// the box is always absent here.
fn feature_map_ranges(map: &FeatureMap) -> Result<RangeSet> {
    let mut rt = RangeBase::new();
    let mut mz = RangeBase::new();
    let mut intensity = RangeBase::new();
    for feature in &map.features {
        rt.extend_value(feature.rt)?;
        mz.extend_value(feature.mz)?;
        intensity.extend_value(f64::from(feature.intensity))?;
    }
    for feature in &map.features {
        if let Some(bounds) = feature.hull_bounding_box() {
            rt.extend_value(bounds.rt_range().min)?;
            rt.extend_value(bounds.rt_range().max)?;
            mz.extend_value(bounds.mz_range().min)?;
            mz.extend_value(bounds.mz_range().max)?;
        }
    }
    Ok(RangeSet {
        rt: to_range(&rt)?,
        mz: to_range(&mz)?,
        mobility: None,
        intensity: to_range(&intensity)?,
        has_mobility: false,
    })
}

fn to_range(range: &RangeBase) -> Result<Option<Range>> {
    if range.is_empty() {
        return Ok(None);
    }
    Ok(Some(Range {
        min: range.min()?,
        max: range.max()?,
    }))
}

fn to_count(count: usize) -> Result<u64> {
    u64::try_from(count).map_err(|_| count_overflow())
}

fn count_overflow() -> Error {
    Error::InvalidValue("FileInfo feature count overflows 64 bits".into())
}

/// One featureXML statistics block: text title, TSV title and the `float`
/// feature field it summarises.
type StatisticsBlock = (&'static str, &'static str, fn(&Feature) -> f32);

/// The five `-s` blocks: intensities, FWHM in RT, overall quality, RT quality
/// and m/z quality, all at `writtenDigits<float>()` precision.
fn write_statistics(
    map: &FeatureMap,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
) -> Result<()> {
    write_statistics_title(os);
    let blocks: [StatisticsBlock; 5] = [
        ("Intensities", "intensities", |feature| feature.intensity),
        (
            "Feature FWHM in RT dimension",
            "feature FWHM in RT dimension",
            |feature| feature.width,
        ),
        ("Overall qualities", "overall qualities", |feature| {
            feature.quality
        }),
        (
            "Qualities in retention time dimension",
            "qualities in retention time dimension",
            |feature| feature.quality_rt,
        ),
        (
            "Qualities in mass-to-charge dimension",
            "qualities in mass-to-charge dimension",
            |feature| feature.quality_mz,
        ),
    ];
    for (title, tsv_title, field) in blocks {
        let mut values = statistics_buffer(map.len())?;
        values.extend(map.features.iter().map(|feature| f64::from(field(feature))));
        let stats = summarize(&mut values)?;
        os.set_precision(WRITTEN_DIGITS_F32);
        os.text(title).text(":\n");
        write_summary_text(os, &stats);
        os.text("\n");
        os_tsv.set_precision(WRITTEN_DIGITS_F32);
        write_summary_tsv(os_tsv, &stats, tsv_title);
    }
    Ok(())
}

/// `writeSummaryStatisticsMachineReadable_`: eight `statistics: <title>: ...`
/// TSV lines with the doubles at the stream's current precision.
fn write_summary_tsv(os_tsv: &mut ReportStream, stats: &SummaryStatistics, title: &str) {
    os_tsv
        .text("statistics: ")
        .text(title)
        .text(": num. of values\t")
        .value(stats.count)
        .text("\n");
    for (label, value) in [
        ("mean", stats.mean),
        ("minimum", stats.min),
        ("lower quartile", stats.lowerq),
        ("median", stats.median),
        ("upper quartile", stats.upperq),
        ("maximum", stats.max),
        ("variance", stats.variance),
    ] {
        os_tsv
            .text("statistics: ")
            .text(title)
            .text(": ")
            .text(label)
            .text("\t")
            .double(value)
            .text("\n");
    }
}
