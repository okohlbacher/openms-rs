// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The FileInfo summary of consensusXML consensus maps
//! (`FORMAT/FileInfo.cpp:1146-1311`, `:1985-1989`, `:2101-2104`, `:2257-2372`).
//!
//! The consensusXML branch of the report. The map is loaded through
//! [`FileHandler::load_consensus_map`](crate::format::FileHandler::load_consensus_map), and the branch then writes, in the
//! source order:
//!
//! 1. the consensus-feature size histogram, in *descending* size, each row
//!    carrying the number of consensus features of that size, the sub-features
//!    they account for, and the same two counts restricted to the consensus
//!    features that carry at least one peptide identification. The size column
//!    is right-aligned in `largest_size / 10 + 1` characters — the source's own
//!    expression, which is not a digit count;
//! 2. one row per number of maps a peptide (sequence and charge together) was
//!    observed in, again descending;
//! 3. the two total lines;
//! 4. one range block over retention time, m/z and intensity. A consensus map
//!    has no mobility dimension. The histogram, the peptide rows, the totals
//!    and the ranges are all skipped when the map holds no consensus feature,
//!    which is reported by two lines of its own instead;
//! 5. the column headers, when there are any;
//! 6. the assigned and unassigned peptide identification counts, which are
//!    written whether or not the map is empty.
//!
//! `-m` adds the document identifier, with no TSV twin, unlike the featureXML
//! branch; `-p` the map's data processing; `-s` eleven statistics blocks, none of
//! which is written to the TSV.
//!
//! The structured [`FeatureInfo`](crate::format::file_info::model::FeatureInfo) and ranges are filled alongside, as the
//! source fills its `Result`.
//!
//! # What this port refuses
//!
//! `FileInfo.cpp:1176-1183` sizes one occurrence vector from the number of
//! column headers and then indexes it with each sub-feature's map index. A map
//! index at or beyond the header count — a map with no `<mapList>` at all
//! among them — is an out-of-bounds `std::vector::operator[]`, which the
//! reference build answers with a segmentation fault or with whatever occupies
//! the adjacent heap. Lead decision D1 refuses exactly there, so
//! [`report`] returns [`Error::InvalidValue`](crate::Error::InvalidValue) before any of the report is
//! written. The same file without a peptide identification on the offending
//! consensus feature never reaches the indexing and is reported normally, as
//! it is by the reference build.
//!
//! See `docs/FILE_INFO_A7_SUPPORT.md` for the evidence and the native
//! differences.

#![cfg(feature = "consensusxml")]

use super::model::{FeatureInfo, FileInfoResult, MapColumn, Options, Range, RangeSet, Ranges};
use super::report::{
    ReportStream, check_statistics_values, statistics_buffer, summarize, write_meta_title,
    write_processing, write_processing_title, write_ranges_text, write_ranges_tsv,
    write_statistics_title, write_summary_text,
};
use super::text_format::{WRITTEN_DIGITS_F32, WRITTEN_DIGITS_F64};
use crate::format::FileHandler;
use crate::format::file_types::FileType;
use crate::kernel::ranges::RangeBase;
use crate::kernel::{ConsensusFeature, ConsensusMap};
use crate::math::x86_64;
use crate::metadata::DataProcessing;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

/// Load the consensus map and write the consensusXML branch, `-m`, `-p` and `-s`.
pub(crate) fn report(
    path: &Path,
    options: &Options,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    result: &mut FileInfoResult,
) -> Result<()> {
    let map = FileHandler::load_consensus_map(path, &[FileType::ConsensusXml])?;
    let ranges = consensus_map_ranges(&map)?;
    let columns = map.column_headers.len();

    // FileInfo.cpp:1153-1186. `size_with_id` is read through operator[] in the
    // print loop, so an absent size reads as zero; a BTreeMap lookup with a
    // zero default is the same reading.
    let mut size_histogram: BTreeMap<u64, u64> = BTreeMap::new();
    let mut size_with_id: BTreeMap<u64, u64> = BTreeMap::new();
    let mut assigned_ids = 0_u64;
    // map<pair<string, UInt>, vector<int>>: the charge is stored as an int and
    // the key converts it to UInt, which for a negative charge is the
    // two's-complement value. The map's order never reaches the report, because
    // the aggregation below only sums into other maps.
    let mut occurrences: BTreeMap<(String, u32), Vec<u64>> = BTreeMap::new();

    for feature in &map.features {
        let size = count(feature.handles().len())?;
        *size_histogram.entry(size).or_insert(0) += 1;
        let identifications = &feature.peptide_identifications;
        assigned_ids = assigned_ids
            .checked_add(count(identifications.len())?)
            .ok_or_else(|| overflow("FileInfo assigned identification count overflows 64 bits"))?;
        let Some(first) = identifications.first() else {
            continue;
        };
        *size_with_id.entry(size).or_insert(0) += 1;
        let Some(hit) = first.hits.first() else {
            continue;
        };
        reject_out_of_range_map_index(feature, columns)?;
        let key = (hit.sequence.to_string(), unsigned_charge(hit.charge));
        let row = occurrences.entry(key).or_insert_with(|| vec![0; columns]);
        for handle in feature.handles() {
            let index = usize::try_from(handle.map_index)
                .map_err(|_| overflow("consensus map index overflows this platform's usize"))?;
            let slot = row
                .get_mut(index)
                .ok_or_else(|| overflow("consensus map index is outside the column headers"))?;
            *slot += 1;
        }
    }

    // FileInfo.cpp:1193-1209: the same sequence and charge seen in n maps
    // counts once, and contributes every sub-feature it was seen in.
    let mut aggregated_consensus: BTreeMap<u64, u64> = BTreeMap::new();
    let mut aggregated_features: BTreeMap<u64, u64> = BTreeMap::new();
    for row in occurrences.values() {
        let maps = row.iter().filter(|&&seen| seen != 0).count();
        let features: u64 = row.iter().sum();
        let maps = count(maps)?;
        *aggregated_consensus.entry(maps).or_insert(0) += 1;
        *aggregated_features.entry(maps).or_insert(0) += features;
    }

    if size_histogram.is_empty() {
        os.text("\nNumber of consensus features: 0\n");
        os.text("No consensus features found, map is empty!\n\n");
    } else {
        let largest = *size_histogram
            .keys()
            .next_back()
            .expect("the histogram is not empty");
        // FileInfo.cpp:1221: the source's own expression, not a digit count.
        let field_width = usize::try_from(largest / 10 + 1)
            .map_err(|_| overflow("consensus size column width overflows this platform's usize"))?;
        os.text("\nNumber of consensus features:\n");

        let mut number_features = 0_u64;
        let mut consensus_with_id = 0_u64;
        let mut features_with_id = 0_u64;
        for (&size, &features_of_size) in size_histogram.iter().rev() {
            let sub_features = size
                .checked_mul(features_of_size)
                .ok_or_else(|| overflow("consensus sub-feature count overflows 64 bits"))?;
            let with_id = size_with_id.get(&size).copied().unwrap_or(0);
            let sub_features_with_id = size
                .checked_mul(with_id)
                .ok_or_else(|| overflow("consensus sub-feature count overflows 64 bits"))?;
            number_features = number_features
                .checked_add(sub_features)
                .ok_or_else(|| overflow("consensus sub-feature count overflows 64 bits"))?;
            features_with_id = features_with_id
                .checked_add(sub_features_with_id)
                .ok_or_else(|| overflow("consensus sub-feature count overflows 64 bits"))?;
            consensus_with_id = consensus_with_id
                .checked_add(with_id)
                .ok_or_else(|| overflow("consensus feature count overflows 64 bits"))?;
            os.text("  of size ")
                .text(&right_aligned(size, field_width))
                .text(": ")
                .value(features_of_size)
                .text("\t (features: ")
                .value(sub_features)
                .text(" )\t with at least one ID: ")
                .value(with_id)
                .text("\t (features: ")
                .value(sub_features_with_id)
                .text(" )\n");
        }

        for ((&maps, &consensus_count), &feature_count) in aggregated_consensus
            .iter()
            .rev()
            .zip(aggregated_features.values().rev())
        {
            os.text("  peptides (with different mod. and charge) observed in ")
                .text(&right_aligned(maps, field_width))
                .text(" maps: ")
                .value(consensus_count)
                .text("\t (features: ")
                .value(feature_count)
                .text(" )\n");
        }

        os.text("  total consensus features:    ")
            .value(map.features.len())
            .text("  with at least one ID: ")
            .value(consensus_with_id)
            .text("\n  total features:              ")
            .value(number_features)
            .text("  with at least one ID: ")
            .text(&" ".repeat(field_width))
            .value(features_with_id)
            .text("\n");

        os.text("Ranges:\n");
        write_ranges_text(os, &ranges, false);
        write_ranges_tsv(os_tsv, "general: ranges: ", &ranges);
    }

    if !map.column_headers.is_empty() {
        os.text("File descriptions:\n");
        for (identifier, header) in &map.column_headers {
            os.text("  ")
                .text(&header.filename)
                .text(":\n    identifier: ")
                .value(identifier)
                .text("\n    label:      ")
                .text(&header.label)
                .text("\n    size:       ")
                .value(header.size)
                .text("\n");
        }
        os.text("\n");
    }

    let unassigned_ids = count(map.unassigned_peptide_identifications.len())?;
    os.text("Assigned peptide identifications: ")
        .value(assigned_ids)
        .text("\n");
    os_tsv
        .text("general: assigned peptide identifications\t")
        .value(assigned_ids)
        .text("\n");
    os.text("Unassigned peptide identifications: ")
        .value(unassigned_ids)
        .text("\n");
    os_tsv
        .text("general: unassigned peptide identifications\t")
        .value(unassigned_ids)
        .text("\n");

    result.feature = Some(FeatureInfo {
        is_consensus: true,
        num_features: count(map.features.len())?,
        size_distribution: size_histogram,
        assigned_ids,
        unassigned_ids,
        map_columns: map
            .column_headers
            .iter()
            .map(|(identifier, header)| {
                Ok(MapColumn {
                    filename: header.filename.clone(),
                    identifier: identifier.to_string(),
                    label: header.label.clone(),
                    size: count(header.size)?,
                })
            })
            .collect::<Result<Vec<_>>>()?,
        ..FeatureInfo::default()
    });
    result.ranges = Ranges {
        combined: ranges,
        is_experiment: false,
        ..Ranges::default()
    };

    if options.meta {
        // FileInfo.cpp:1985-1989: unlike the featureXML arm, no TSV line.
        write_meta_title(os);
        os.text("Document ID: ").text(&map.identifier).text("\n\n");
    }
    if options.processing {
        write_processing_title(os);
        let processing: Vec<&DataProcessing> = map.data_processing.iter().collect();
        write_processing(os, os_tsv, &processing, result);
    }
    if options.statistics {
        write_statistics(&map, os)?;
    }
    Ok(())
}

/// Refuse a sub-feature whose map index is at or beyond the number of column
/// headers, which `FileInfo.cpp:1183` would use to index a vector sized from
/// that number.
///
/// # Errors
///
/// [`Error::InvalidValue`](crate::Error::InvalidValue) naming the index and the header count.
fn reject_out_of_range_map_index(feature: &ConsensusFeature, columns: usize) -> Result<()> {
    let columns_available = u64::try_from(columns)
        .map_err(|_| overflow("consensus column-header count overflows 64 bits"))?;
    for handle in feature.handles() {
        if handle.map_index >= columns_available {
            return Err(Error::InvalidValue(format!(
                "FileInfo consensusXML branch: a consensus feature with a peptide \
                 identification holds a sub-feature of map {} while the map has {columns} column \
                 header(s); the source indexes an occurrence vector of that length with the map \
                 index (FileInfo.cpp:1176-1183), which is out of bounds",
                handle.map_index
            )));
        }
    }
    Ok(())
}

/// `ConsensusMap::updateRanges` (`ConsensusMap.cpp`): for every consensus
/// feature, its own position and intensity and then those of each of its
/// sub-features, in that interleaved order, each extended with the source's
/// keep-first rule on equal endpoints ([`RangeBase::extend_value`]).
///
/// A consensus map carries no mobility dimension.
fn consensus_map_ranges(map: &ConsensusMap) -> Result<RangeSet> {
    let mut rt = RangeBase::new();
    let mut mz = RangeBase::new();
    let mut intensity = RangeBase::new();
    for feature in &map.features {
        rt.extend_value(feature.rt)?;
        mz.extend_value(feature.mz)?;
        intensity.extend_value(f64::from(feature.intensity))?;
        for handle in feature.handles() {
            rt.extend_value(handle.rt)?;
            mz.extend_value(handle.mz)?;
            intensity.extend_value(f64::from(handle.intensity))?;
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

/// The samples the `-s` block summarises (`FileInfo.cpp:2259-2328`).
///
/// `qualities` and `widths` reproduce a source defect: both are declared as
/// `vector<double> qualities(size)` — `size` zero-initialised values — and then
/// appended to, so each ends up holding twice as many values as there are
/// consensus features, the first half of them `0.0`. `intensities` is declared
/// empty and merely reserved, so it is the length one would expect. The
/// upstream reference output `FileInfo_7_output.txt` records the difference:
/// five consensus features, `num. of values: 5` for the intensities and
/// `num. of values: 10` for the qualities. `widths` is collected and never
/// printed, so its defect is not observable; it is not collected here.
struct Samples {
    intensities: Vec<f64>,
    qualities: Vec<f64>,
    rt_delta_by_elems: Vec<f64>,
    rt_aad_by_elems: Vec<f64>,
    rt_aad_by_cfs: Vec<f64>,
    mz_delta_by_elems: Vec<f64>,
    mz_aad_by_elems: Vec<f64>,
    mz_aad_by_cfs: Vec<f64>,
    it_delta_by_elems: Vec<f64>,
    it_aad_by_elems: Vec<f64>,
    it_aad_by_cfs: Vec<f64>,
}

fn collect(map: &ConsensusMap) -> Result<Samples> {
    let size = map.features.len();
    let elements: usize = map
        .features
        .iter()
        .try_fold(0_usize, |total, feature| {
            total.checked_add(feature.handles().len())
        })
        .ok_or_else(|| overflow("consensus sub-feature count overflows this platform's usize"))?;
    check_statistics_values(elements)?;
    let mut samples = Samples {
        intensities: statistics_buffer(size)?,
        // The pre-sized half the source zero-initialises, then the values.
        qualities: {
            let mut values = statistics_buffer(size.checked_mul(2).ok_or_else(|| {
                overflow("consensus quality sample overflows this platform's usize")
            })?)?;
            values.resize(size, 0.0);
            values
        },
        rt_delta_by_elems: statistics_buffer(elements)?,
        rt_aad_by_elems: statistics_buffer(elements)?,
        rt_aad_by_cfs: statistics_buffer(size)?,
        mz_delta_by_elems: statistics_buffer(elements)?,
        mz_aad_by_elems: statistics_buffer(elements)?,
        mz_aad_by_cfs: statistics_buffer(size)?,
        it_delta_by_elems: statistics_buffer(elements)?,
        it_aad_by_elems: statistics_buffer(elements)?,
        it_aad_by_cfs: statistics_buffer(size)?,
    };
    for feature in &map.features {
        let mut rt_aad = 0.0_f64;
        let mut mz_aad = 0.0_f64;
        let mut it_aad = 0.0_f64;
        samples.intensities.push(f64::from(feature.intensity));
        samples.qualities.push(f64::from(feature.quality));
        let centre_intensity = f64::from(feature.intensity);
        let denominator = if centre_intensity > 0.0 {
            centre_intensity
        } else {
            1.0
        };
        for handle in feature.handles() {
            let rt_diff = handle.rt - feature.rt;
            samples.rt_delta_by_elems.push(rt_diff);
            let rt_diff = if rt_diff < 0.0 { -rt_diff } else { rt_diff };
            samples.rt_aad_by_elems.push(rt_diff);
            rt_aad += rt_diff;

            let mz_diff = handle.mz - feature.mz;
            samples.mz_delta_by_elems.push(mz_diff);
            let mz_diff = if mz_diff < 0.0 { -mz_diff } else { mz_diff };
            samples.mz_aad_by_elems.push(mz_diff);
            mz_aad += mz_diff;

            // FileInfo.cpp:2310-2315. A sub-feature of intensity 0 under a
            // centroid of positive intensity gives 1 / 0 = +inf here, and the
            // summary of a sample holding it has an infinite mean and a NaN
            // variance; one of intensity -0.0 gives 1 / -0.0 = -inf, so the
            // `it_aad` below can be (-inf) + (+inf) = NaN and the NaN reaches
            // the per-consensus-feature sample itself. All of that is in bounds
            // and reproducible, so D1 keeps it, and it is pinned by
            // `consensus_zero_intensity_sub_feature_makes_the_variance_a_nan`
            // and `consensus_nan_in_the_statistics_sample`. Since D16 there is
            // no shape here that `SummaryStatistics` declines: it sorts with
            // `std::sort(begin, end)` itself and reads its order statistics out
            // of whatever that leaves. Section 5.2 of
            // docs/FILE_INFO_A7_SUPPORT.md records the measurement.
            let it_ratio = f64::from(handle.intensity) / denominator;
            samples.it_delta_by_elems.push(it_ratio);
            let it_ratio = if it_ratio < 1.0 {
                1.0 / it_ratio
            } else {
                it_ratio
            };
            samples.it_aad_by_elems.push(it_ratio);
            // The one place in this module where plain Rust arithmetic can
            // *generate* a NaN, so the one place whose NaN bits would otherwise
            // be the host's: `addsd` answers an invalid (-inf) + (+inf) with
            // SSE2's default NaN `0xfff8000000000000`, whose sign bit is set,
            // while an AArch64 host produces the positive default NaN. The text
            // layer spells those two differently (`-nan` against `nan`), so the
            // spelling is only honest once the value is.
            //
            // Nothing else here needs it. `map.validate()` in the consensusXML
            // reader refuses a non-finite rt, m/z, intensity or width
            // (`ConsensusFeature::validate` -> `validate_values`,
            // `src/kernel/features.rs:954-964`), so `rt_diff` and `mz_diff` are
            // differences of finite values and `it_ratio` a quotient of one by a
            // non-zero one; the `rt_aad`/`mz_aad` sums accumulate values already
            // made non-negative, which cannot cancel to a NaN; and the three
            // `/= handles` divisors are finite counts of at least one. An
            // overflow to an infinity is IEEE-determined and so host-independent
            // either way.
            it_aad = x86_64::add(it_aad, it_ratio);
        }
        if !feature.handles().is_empty() {
            let handles = as_double(count(feature.handles().len())?);
            rt_aad /= handles;
            mz_aad /= handles;
            it_aad /= handles;
        }
        samples.rt_aad_by_cfs.push(rt_aad);
        samples.mz_aad_by_cfs.push(mz_aad);
        samples.it_aad_by_cfs.push(it_aad);
    }
    Ok(samples)
}

/// The eleven `-s` blocks, in the source order and at the source precisions:
/// `writtenDigits<float>()` for the two consensus-feature blocks and the three
/// intensity-ratio blocks, `writtenDigits<double>()` for the six retention-time
/// and mass-to-charge blocks. None of them is written to the TSV report.
fn write_statistics(map: &ConsensusMap, os: &mut ReportStream) -> Result<()> {
    write_statistics_title(os);
    let mut samples = collect(map)?;
    let blocks: [(u32, &str, &mut Vec<f64>); 11] = [
        (
            WRITTEN_DIGITS_F32,
            "Intensities of consensus features:",
            &mut samples.intensities,
        ),
        (
            WRITTEN_DIGITS_F32,
            "Qualities of consensus features:",
            &mut samples.qualities,
        ),
        (
            WRITTEN_DIGITS_F64,
            "Retention time differences (\"element - center\", weight 1 per element):",
            &mut samples.rt_delta_by_elems,
        ),
        (
            WRITTEN_DIGITS_F64,
            "Absolute retention time differences (\"|element - center|\", weight 1 per element):",
            &mut samples.rt_aad_by_elems,
        ),
        (
            WRITTEN_DIGITS_F64,
            "Average absolute differences of retention time within consensus features \
             (\"|element - center|\", weight 1 per consensus features):",
            &mut samples.rt_aad_by_cfs,
        ),
        (
            WRITTEN_DIGITS_F64,
            "Mass-to-charge differences (\"element - center\", weight 1 per element):",
            &mut samples.mz_delta_by_elems,
        ),
        (
            WRITTEN_DIGITS_F64,
            "Absolute differences of mass-to-charge (\"|element - center|\", weight 1 per element):",
            &mut samples.mz_aad_by_elems,
        ),
        (
            WRITTEN_DIGITS_F64,
            "Average absolute differences of mass-to-charge within consensus features \
             (\"|element - center|\", weight 1 per consensus features):",
            &mut samples.mz_aad_by_cfs,
        ),
        (
            WRITTEN_DIGITS_F32,
            "Intensity ratios (\"element / center\", weight 1 per element):",
            &mut samples.it_delta_by_elems,
        ),
        (
            WRITTEN_DIGITS_F32,
            "Relative intensity error (\"max{(element / center), (center / element)}\", \
             weight 1 per element):",
            &mut samples.it_aad_by_elems,
        ),
        (
            WRITTEN_DIGITS_F32,
            "Average relative intensity error within consensus features \
             (\"max{(element / center), (center / element)}\", weight 1 per consensus features):",
            &mut samples.it_aad_by_cfs,
        ),
    ];
    for (precision, title, values) in blocks {
        let stats = summarize(values)?;
        os.set_precision(precision);
        os.text(title).text("\n");
        write_summary_text(os, &stats);
        os.text("\n");
    }
    Ok(())
}

/// `setw(width)` on a `Size`: the decimal text right-aligned in `width`
/// characters, and left as it is when it is already longer.
fn right_aligned(value: u64, width: usize) -> String {
    let text = value.to_string();
    if text.len() >= width {
        return text;
    }
    let mut padded = " ".repeat(width - text.len());
    padded.push_str(&text);
    padded
}

#[expect(
    clippy::cast_precision_loss,
    reason = "the source converts its Size counts to double in exactly this place"
)]
fn as_double(value: u64) -> f64 {
    value as f64
}

/// The `int` to `UInt` conversion the source's `make_pair(s, z)` performs when
/// it builds a `pair<std::string, UInt>` key from an `int` charge: the
/// two's-complement value, which is what the C++ standard specifies for a
/// conversion to an unsigned type.
#[expect(
    clippy::cast_sign_loss,
    reason = "the source converts the int charge to UInt in exactly this place"
)]
fn unsigned_charge(charge: i32) -> u32 {
    charge as u32
}

fn count(value: usize) -> Result<u64> {
    u64::try_from(value).map_err(|_| overflow("consensus count overflows 64 bits"))
}

fn overflow(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
