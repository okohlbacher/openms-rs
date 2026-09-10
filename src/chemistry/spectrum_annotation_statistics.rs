// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source statistics and label grammars, scoped to staged annotation arrays.

use super::*;
use crate::kernel::{Precursor, nearest};

pub(super) fn compute(
    options: &SpectrumAnnotator,
    spectrum: &mut MSSpectrum,
    hit: &PeptideHit,
    precursors: &[Precursor],
    alignment: &SpectrumAlignment,
    work: &mut AnnotationWork,
) -> Result<MetaInfo> {
    if hit.sequence.is_empty() {
        return Err(invalid("annotation statistics require a nonempty peptide"));
    }
    if options.precursor_statistics && precursors.len() > MAX_ANNOTATION_PRECURSORS {
        return Err(invalid("annotation precursor limit exceeded"));
    }
    // At most twenty fixed-size metadata nodes/keys, plus separately charged labels.
    work.allocate::<u8>(20 * 256)?;
    sort_annotated(spectrum, true, work)?;
    let n = spectrum.len();
    work.consume(n)?;
    work.allocate::<&str>(n)?;
    work.allocate::<f64>(n)?;
    let mut ions = Vec::with_capacity(n);
    let mut errors = Vec::with_capacity(n);
    let mut total = 0.0;
    let mut matched = 0.0;
    let mut nterm = 0.0;
    let mut cterm = 0.0;
    let series_len = hit.sequence.len() - 1;
    if options.max_series {
        work.allocate::<bool>(6 * series_len)?;
        work.consume(6 * series_len)?;
    }
    let mut series: [Vec<bool>; 6] = std::array::from_fn(|_| {
        if options.max_series {
            vec![false; series_len]
        } else {
            Vec::new()
        }
    });
    for (index, peak) in spectrum.peaks.iter().enumerate() {
        total += f64::from(peak.intensity);
        let label = &spectrum.string_data_arrays[0].data[index];
        if label.is_empty() {
            continue;
        }
        work.consume(label.len())?;
        errors.push(f64::from(spectrum.float_data_arrays[0].data[index]));
        matched += f64::from(peak.intensity);
        if options.terminal_series_match_ratio {
            match terminal_series(label) {
                Some(true) => nterm += f64::from(peak.intensity),
                Some(false) => cterm += f64::from(peak.intensity),
                None => {}
            }
        }
        if options.max_series {
            if let Some((kind, ordinal)) = series_position(label) {
                if ordinal == 0 || ordinal > series_len {
                    // Source skips list inclusion AFTER recording error/intensity.
                    continue;
                }
                series[kind][ordinal - 1] = true;
            }
        }
        ions.push(label.as_str());
    }
    let mut updates = MetaInfo::new();
    if options.basic_statistics {
        let bytes = ions
            .iter()
            .try_fold(0_usize, |total, name| total.checked_add(name.len()))
            .and_then(|bytes| bytes.checked_add(ions.len().saturating_sub(1)))
            .ok_or_else(|| invalid("matched-ion metadata length overflow"))?;
        work.copy_label(bytes)?;
        updates.insert("matched_ions".into(), ions.join(",").into());
        set_float(&mut updates, "matched_intensity", matched)?;
        updates.insert("matched_ion_number".into(), (ions.len() as i64).into());
        updates.insert("peak_number".into(), (n as i64).into());
        set_float(&mut updates, "sum_intensity", total)?;
    }
    if options.terminal_series_match_ratio {
        set_float(&mut updates, "NTermIonCurrentRatio", nterm / matched)?;
        set_float(&mut updates, "CTermIonCurrentRatio", cterm / matched)?;
    }
    if options.top_n_fragment_errors != 0 {
        error_statistics(&mut updates, &errors, options.top_n_fragment_errors, work)?;
    }
    if options.max_series {
        let mut longest = 0;
        let mut kind = "";
        for (name, positions) in ["a", "b", "c", "x", "y", "z"].into_iter().zip(series) {
            let mut stretch = 0;
            for present in positions {
                stretch = if present { stretch + 1 } else { 0 };
                if stretch > longest {
                    longest = stretch;
                    kind = name;
                }
            }
        }
        updates.insert("max_series_type".into(), kind.into());
        updates.insert("max_series_size".into(), (longest as i64).into());
    }
    if options.sn_statistics {
        let sn = if n == ions.len() {
            0.0_f32
        } else {
            ((matched / ions.len() as f64) / ((total - matched) / (n - ions.len()) as f64)) as f32
        };
        set_float(&mut updates, "sn_by_matched_intensity", f64::from(sn))?;
        let median = if n % 2 == 0 {
            (spectrum.peaks[n / 2 - 1].intensity + spectrum.peaks[n / 2].intensity) / 2.0_f32
        } else {
            spectrum.peaks[n / 2].intensity
        };
        let (mut signal, mut noise) = (0.0_f32, 0.0_f32);
        let (mut signal_count, mut noise_count) = (0_usize, 0_usize);
        work.consume(n)?;
        for peak in &spectrum.peaks {
            if peak.intensity <= median {
                noise += peak.intensity;
                noise_count += 1;
            } else {
                signal += peak.intensity;
                signal_count += 1;
            }
        }
        let sn = if signal_count == 0 || noise_count == 0 {
            0.0_f32
        } else {
            (signal / signal_count as f32) / (noise / noise_count as f32)
        };
        set_float(&mut updates, "sn_by_median_intensity", f64::from(sn))?;
    }
    if options.precursor_statistics {
        let mut found = false;
        if !precursors.is_empty() {
            // Repeating this sort for every precursor has no additional effect.
            sort_annotated(spectrum, false, work)?;
        }
        let tolerance = raw_tolerance(alignment);
        work.consume(
            precursors
                .len()
                .checked_mul(usize::BITS as usize + 2)
                .ok_or_else(|| invalid("precursor search work overflow"))?,
        )?;
        for precursor in precursors {
            let mz = finite(precursor.mz, "annotation precursor m/z is nonfinite")?;
            if let Some(index) = nearest(&spectrum.peaks, mz, |p| p.mz) {
                let inside = |i: usize| {
                    spectrum
                        .peaks
                        .get(i)
                        .is_some_and(|p| p.mz >= mz - tolerance && p.mz <= mz + tolerance)
                };
                if inside(index) {
                    found = true;
                } else {
                    let other = if spectrum.peaks[index].mz < mz {
                        index.checked_add(1)
                    } else {
                        index.checked_sub(1)
                    };
                    if other.is_some_and(inside) {
                        found = true;
                    }
                }
            }
        }
        updates.insert("precursor_in_ms2".into(), i32::from(found).into());
    }
    Ok(updates)
}

fn error_statistics(
    updates: &mut MetaInfo,
    errors: &[f64],
    top_n: usize,
    work: &mut AnnotationWork,
) -> Result<()> {
    if errors.is_empty() {
        for key in [
            "median_fragment_error",
            "IQR_fragment_error",
            "topN_meanfragmenterror",
            "topN_MSEfragmenterror",
            "topN_stddevfragmenterror",
        ] {
            updates.insert(key.into(), 0_i32.into());
        }
        return Ok(());
    }
    if top_n == 1 {
        return Err(invalid(
            "sample fragment-error deviation requires top N at least two",
        ));
    }
    work.allocate::<f64>(errors.len())?;
    work.sort::<f64>(errors.len())?;
    let mut sorted = errors.to_vec();
    sorted.sort_by(f64::total_cmp);
    let n = sorted.len();
    let mid = n / 2;
    let median = if n % 2 == 0 {
        (sorted[mid] + sorted[mid - 1]) / 2.0
    } else {
        sorted[mid]
    };
    set_float(updates, "median_fragment_error", median)?;
    // Safe extension of source intended order statistics for lists of1..3.
    set_float(
        updates,
        "IQR_fragment_error",
        sorted[n / 4 + mid] - sorted[n / 4],
    )?;
    work.allocate::<f64>(top_n)?;
    work.consume(
        top_n
            .checked_mul(3)
            .ok_or_else(|| invalid("top-N statistics work overflow"))?,
    )?;
    let mut top = vec![0.0; top_n];
    for (slot, error) in top.iter_mut().zip(errors.iter().rev()) {
        *slot = *error;
    }
    let mean = top.iter().sum::<f64>() / top_n as f64;
    let mut variance = 0.0;
    let mut square_sum = 0.0;
    for error in top {
        let difference = error - mean;
        variance += difference * difference;
        square_sum += error * error;
    }
    set_float(updates, "topN_meanfragmenterror", mean)?;
    set_float(updates, "topN_MSEfragmenterror", square_sum / top_n as f64)?;
    set_float(
        updates,
        "topN_stddevfragmenterror",
        (variance / (top_n - 1) as f64).sqrt(),
    )?;
    Ok(())
}
fn set_float(updates: &mut MetaInfo, key: &str, value: f64) -> Result<()> {
    if !value.is_finite() {
        return Err(invalid(&format!(
            "nonfinite enabled annotation statistic {key}"
        )));
    }
    updates.insert(key.into(), MetaValue::try_from(value)?);
    Ok(())
}

/// Source complete-name [a,b,c]/[x,y,z] digit+ plus* grammars. Comma belongs to
/// BOTH character classes; the N-terminal branch has priority as in the source.
fn terminal_series(label: &str) -> Option<bool> {
    let bytes = label.as_bytes();
    let first = *bytes.first()?;
    if !b"abcxyz,".contains(&first) {
        return None;
    }
    let end = bytes[1..]
        .iter()
        .position(|b| !b.is_ascii_digit())
        .map_or(bytes.len(), |i| i + 1);
    if end == 1 || !bytes[end..].iter().all(|&b| b == b'+') {
        return None;
    }
    Some(b"abc,".contains(&first))
}
/// Source ordinal regex permits a sign/comma run, ASCII word tail, then pluses.
/// It deliberately does not parse internal-fragment colon or mzPAF labels.
fn series_position(label: &str) -> Option<(usize, usize)> {
    let bytes = label.as_bytes();
    let kind = b"abcxyz".iter().position(|b| Some(b) == bytes.first())?;
    let mut end = 1;
    let mut ordinal = 0_usize;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        ordinal = ordinal
            .checked_mul(10)?
            .checked_add(usize::from(bytes[end] - b'0'))?;
        end += 1;
    }
    if end == 1 {
        return None;
    }
    while end < bytes.len() && b"+,-".contains(&bytes[end]) {
        end += 1;
    }
    while end < bytes.len() && (bytes[end].is_ascii_alphanumeric() || bytes[end] == b'_') {
        end += 1;
    }
    while end < bytes.len() && bytes[end] == b'+' {
        end += 1;
    }
    (end == bytes.len()).then_some((kind, ordinal))
}

fn sort_annotated(
    spectrum: &mut MSSpectrum,
    intensity: bool,
    work: &mut AnnotationWork,
) -> Result<()> {
    let n = spectrum.len();
    work.sort::<usize>(n)?;
    work.allocate::<usize>(n)?;
    work.allocate::<usize>(n)?;
    work.consume(
        n.checked_mul(6)
            .ok_or_else(|| invalid("annotation permutation work overflow"))?,
    )?;
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        if intensity {
            spectrum.peaks[a]
                .intensity
                .partial_cmp(&spectrum.peaks[b].intensity)
                .unwrap()
        } else {
            spectrum.peaks[a]
                .mz
                .partial_cmp(&spectrum.peaks[b].mz)
                .unwrap()
        }
    });
    let mut destination = vec![0; n];
    for (new, old) in order.into_iter().enumerate() {
        destination[old] = new;
    }
    // In-place cycles move owned labels rather than cloning them through selection.
    for index in 0..n {
        while destination[index] != index {
            let next = destination[index];
            spectrum.peaks.swap(index, next);
            spectrum.string_data_arrays[0].data.swap(index, next);
            spectrum.float_data_arrays[0].data.swap(index, next);
            spectrum.integer_data_arrays[0].data.swap(index, next);
            destination.swap(index, next);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_label_grammars_are_distinct_from_charge_name_parsing() {
        assert_eq!(terminal_series(",3++"), Some(true));
        assert_eq!(terminal_series("y3++"), Some(false));
        for label in [
            "b3-H2O++",
            "b3+9",
            "b3,foo+",
            "y3^2",
            "b2:AB++",
            "b3-[H2O]+",
        ] {
            assert_eq!(terminal_series(label), None, "{label}");
        }
        for label in ["b3-H2O++", "b3+9", "b3,foo+"] {
            assert_eq!(series_position(label), Some((1, 3)), "{label}");
        }
        for label in ["y3^2", "b2:AB++", "b3-[H2O]+", ",3++"] {
            assert_eq!(series_position(label), None, "{label}");
        }
        assert_eq!(series_position("b0+"), Some((1, 0)));
    }
}
