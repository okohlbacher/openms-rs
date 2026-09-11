// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source spectrum/identification annotations with atomic, scoped updates.
//! Matching replaces every measured data array; unrelated acquisition metadata,
//! attached identifications and existing identification metadata are not cloned.

use super::theoretical::{TheoreticalGenerationWork, TheoreticalSpectrumGenerator};
use crate::comparison::{AlignmentWork, SpectrumAlignment, Tolerance};
use crate::identification::{PeakAnnotation, PeptideHit, PeptideIdentification};
use crate::kernel::{DataArray, MSSpectrum, Peak1D};
use crate::metadata::{MetaInfo, MetaValue};
use crate::{Error, Result};
use std::mem::size_of;

#[path = "spectrum_annotation_statistics.rs"]
mod statistics;

pub const MAX_ANNOTATION_PEAKS: usize = 1_000_000;
pub const MAX_ANNOTATION_HITS: usize = 1_000;
pub const MAX_ANNOTATION_TOP_N: usize = 1_000_000;
pub const MAX_ANNOTATION_PRECURSORS: usize = 100_000;
/// Cumulative visits/comparisons for the annotation/statistics stages.
pub const MAX_ANNOTATION_WORK: usize = 50_000_000;
/// Cumulative copied label bytes, including matched-ion metadata strings.
pub const MAX_ANNOTATION_LABEL_BYTES: usize = 64 * 1024 * 1024;
/// Conservative cumulative annotation allocation payload, separate from the
/// generator's and aligner's own bounded working storage.
pub const MAX_ANNOTATION_BYTES: usize = 128 * 1024 * 1024;

/// Options corresponding to all eight SpectrumAnnotator source parameters.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpectrumAnnotator {
    pub basic_statistics: bool,
    /// Stored source option with no effect: basic_statistics controls matched_ions.
    pub list_of_ions_matched: bool,
    pub max_series: bool,
    pub sn_statistics: bool,
    pub precursor_statistics: bool,
    /// Zero disables ALL five error statistics; short matched lists are padded
    /// with zeros to exactly this length. A nonempty list with N=1 errors.
    pub top_n_fragment_errors: usize,
    /// Stored source option with no effect: top_n_fragment_errors controls this block.
    pub fragment_error_statistics: bool,
    pub terminal_series_match_ratio: bool,
}
impl Default for SpectrumAnnotator {
    fn default() -> Self {
        Self {
            basic_statistics: true,
            list_of_ions_matched: true,
            max_series: true,
            sn_statistics: true,
            precursor_statistics: true,
            top_n_fragment_errors: 7,
            fragment_error_statistics: true,
            terminal_series_match_ratio: true,
        }
    }
}

impl SpectrumAnnotator {
    /// Sort measured peaks by m/z and replace ALL data arrays with IonNames,
    /// Charges and f32 absolute-Dalton IonMatchError. Reused ppm targets retain
    /// the LAST aligned theoretical label. Unmatched entries are empty/zero.
    /// Required theoretical arrays may be absent only when there are no matches.
    pub fn annotate_matches(
        &self,
        spectrum: &mut MSSpectrum,
        hit: &PeptideHit,
        generator: &TheoreticalSpectrumGenerator,
        alignment: &SpectrumAlignment,
    ) -> Result<()> {
        let mut work = AnnotationWork::default();
        let matched = match_spectrum(&spectrum.peaks, hit, generator, alignment, &mut work)?;
        let staged = matched.into_annotations(&mut work)?;
        let updates = tolerance_metadata(alignment)?;
        for (key, _) in &updates {
            crate::kernel::data_array::Meter {
                work: &mut work.remaining,
                bytes: &mut work.bytes,
            }
            .meta_update(&spectrum.metadata, key)?;
        }
        spectrum.array_descriptions_with_budget(&mut work.remaining, &mut work.bytes)?;
        commit_spectrum(spectrum, staged, updates);
        Ok(())
    }

    /// Replace hit peak annotations, leaving measured input and all other hit
    /// fields intact. Matched-only output keeps one entry per alignment (including
    /// repeated measured targets); all-peaks output uses the last match per target
    /// and is ordered by measured m/z. Missing theoretical labels/charges use
    /// empty/zero values, as in the source's named-array lookup.
    pub fn add_peak_annotations(
        &self,
        hit: &mut PeptideHit,
        spectrum: &MSSpectrum,
        generator: &TheoreticalSpectrumGenerator,
        alignment: &SpectrumAlignment,
        include_unmatched_peaks: bool,
    ) -> Result<()> {
        let mut work = AnnotationWork::default();
        let matched = match_spectrum(&spectrum.peaks, hit, generator, alignment, &mut work)?;
        let names = matched
            .theoretical
            .string_data_arrays
            .iter()
            .find(|a| a.name == "IonNames");
        let charges = matched
            .theoretical
            .integer_data_arrays
            .iter()
            .find(|a| a.name == "Charges");
        let count = if include_unmatched_peaks {
            matched.measured.len()
        } else {
            matched.pairs.len()
        };
        work.allocate::<PeakAnnotation>(count)?;
        work.consume(count)?;
        let last_matches = if include_unmatched_peaks {
            work.allocate::<Option<usize>>(matched.measured.len())?;
            work.consume(matched.measured.len() + matched.pairs.len())?;
            let mut last = vec![None; matched.measured.len()];
            for &(theoretical, observed) in &matched.pairs {
                last[observed] = Some(theoretical);
            }
            Some(last)
        } else {
            None
        };
        let mut annotations = Vec::with_capacity(count);
        let mut add = |observed: usize, theoretical: Option<usize>| -> Result<()> {
            let peak = matched.measured.peaks[observed];
            let label = theoretical
                .and_then(|i| names.and_then(|a| a.data.get(i)))
                .map_or("", String::as_str);
            work.copy_label(label.len())?;
            annotations.push(PeakAnnotation {
                mz: peak.mz,
                intensity: f64::from(peak.intensity),
                charge: theoretical
                    .and_then(|i| charges.and_then(|a| a.data.get(i)))
                    .copied()
                    .unwrap_or(0),
                annotation: label.to_owned(),
            });
            Ok(())
        };
        if let Some(last_match) = last_matches {
            // The distinct mapping is deliberate: matched-only preserves duplicates.
            for (observed, theoretical) in last_match.into_iter().enumerate() {
                add(observed, theoretical)?;
            }
        } else {
            for &(theoretical, observed) in &matched.pairs {
                add(observed, Some(theoretical))?;
            }
        }
        hit.peak_annotations = annotations;
        Ok(())
    }

    /// Add all enabled source statistics and leave the final hit's annotation
    /// arrays on the measured spectrum. Empty spectra or identifications with no
    /// hits are immediate no-ops. Only affected fields are staged; any later error
    /// leaves BOTH inputs unchanged. Nonfinite enabled statistics error rather
    /// than writing NaN/Inf metadata. Small-list quartiles use the source's intended
    /// sorted indices safely, replacing invalid C++ nth_element subranges.
    pub fn add_ion_match_statistics(
        &self,
        identification: &mut PeptideIdentification,
        spectrum: &mut MSSpectrum,
        generator: &TheoreticalSpectrumGenerator,
        alignment: &SpectrumAlignment,
    ) -> Result<()> {
        if spectrum.is_empty() || identification.hits.is_empty() {
            return Ok(());
        }
        if identification.hits.len() > MAX_ANNOTATION_HITS {
            return Err(invalid("annotation hit limit exceeded"));
        }
        if self.top_n_fragment_errors > MAX_ANNOTATION_TOP_N {
            return Err(invalid("annotation top-N limit exceeded"));
        }
        let mut work = AnnotationWork::default();
        work.allocate::<MetaInfo>(identification.hits.len())?;
        let mut hit_updates = Vec::with_capacity(identification.hits.len());
        let mut last: Option<MSSpectrum> = None;
        for hit in &identification.hits {
            let peaks = last
                .as_ref()
                .map_or(spectrum.peaks.as_slice(), |s| s.peaks.as_slice());
            let matched = match_spectrum(peaks, hit, generator, alignment, &mut work)?;
            let mut staged = matched.into_annotations(&mut work)?;
            let updates = statistics::compute(
                self,
                &mut staged,
                hit,
                &spectrum.precursors,
                alignment,
                &mut work,
            )?;
            hit_updates.push(updates);
            last = Some(staged);
        }
        let updates = tolerance_metadata(alignment)?;
        for (key, _) in &updates {
            crate::kernel::data_array::Meter {
                work: &mut work.remaining,
                bytes: &mut work.bytes,
            }
            .meta_update(&spectrum.metadata, key)?;
        }
        let tolerance = MetaValue::try_from(raw_tolerance(alignment))?;
        spectrum.array_descriptions_with_budget(&mut work.remaining, &mut work.bytes)?;
        // No fallible scientific operation remains once mutations start.
        for (hit, updates) in identification.hits.iter_mut().zip(hit_updates) {
            hit.metadata.extend(updates);
        }
        identification
            .metadata
            .insert("fragment_match_tolerance".into(), tolerance);
        commit_spectrum(
            spectrum,
            last.expect("nonempty hits produce a stage"),
            updates,
        );
        Ok(())
    }
}

fn commit_spectrum(
    spectrum: &mut MSSpectrum,
    staged: MSSpectrum,
    updates: [(String, MetaValue); 2],
) {
    spectrum.peaks = staged.peaks;
    spectrum.float_data_arrays = staged.float_data_arrays;
    spectrum.integer_data_arrays = staged.integer_data_arrays;
    spectrum.string_data_arrays = staged.string_data_arrays;
    spectrum.metadata.extend(updates);
}
fn raw_tolerance(alignment: &SpectrumAlignment) -> f64 {
    match alignment.tolerance {
        Tolerance::Absolute(value) | Tolerance::Ppm(value) => value,
    }
}
fn tolerance_metadata(alignment: &SpectrumAlignment) -> Result<[(String, MetaValue); 2]> {
    let tolerance = finite(
        raw_tolerance(alignment),
        "annotation tolerance is nonfinite",
    )?;
    Ok([
        (
            "fragment_mass_tolerance".into(),
            MetaValue::try_from(tolerance)?,
        ),
        (
            "fragment_mass_tolerance_ppm".into(),
            i64::from(matches!(alignment.tolerance, Tolerance::Ppm(_))).into(),
        ),
    ])
}

struct MatchedSpectrum {
    measured: MSSpectrum,
    theoretical: MSSpectrum,
    pairs: Vec<(usize, usize)>,
}
impl MatchedSpectrum {
    fn into_annotations(self, work: &mut AnnotationWork) -> Result<MSSpectrum> {
        let mut measured = self.measured;
        let names = self.theoretical.string_data_arrays.first();
        let charges = self.theoretical.integer_data_arrays.first();
        if !self.pairs.is_empty() && (names.is_none() || charges.is_none()) {
            return Err(invalid(
                "matched theoretical spectrum lacks annotation arrays",
            ));
        }
        let n = measured.len();
        work.allocate::<f32>(n)?;
        work.allocate::<i32>(n)?;
        work.allocate::<String>(n)?;
        work.consume(n + self.pairs.len())?;
        let mut errors = vec![0.0; n];
        let mut labels = vec![String::new(); n];
        let mut charge_values = vec![0; n];
        for (theoretical, observed) in self.pairs {
            let label = names
                .and_then(|a| a.data.get(theoretical))
                .ok_or_else(|| invalid("theoretical ion-name array is not aligned"))?;
            let charge = charges
                .and_then(|a| a.data.get(theoretical))
                .ok_or_else(|| invalid("theoretical charge array is not aligned"))?;
            let error =
                (measured.peaks[observed].mz - self.theoretical.peaks[theoretical].mz).abs() as f32;
            finite(f64::from(error), "annotation mass error exceeds f32 range")?;
            work.copy_label(label.len())?;
            labels[observed] = label.clone();
            errors[observed] = error;
            charge_values[observed] = *charge;
        }
        measured
            .float_data_arrays
            .push(DataArray::new("IonMatchError", errors));
        measured
            .integer_data_arrays
            .push(DataArray::new("Charges", charge_values));
        measured
            .string_data_arrays
            .push(DataArray::new("IonNames", labels));
        Ok(measured)
    }
}
fn match_spectrum(
    peaks: &[Peak1D],
    hit: &PeptideHit,
    generator: &TheoreticalSpectrumGenerator,
    alignment: &SpectrumAlignment,
    work: &mut AnnotationWork,
) -> Result<MatchedSpectrum> {
    if peaks.len() > MAX_ANNOTATION_PEAKS {
        return Err(invalid("annotation peak limit exceeded"));
    }
    if hit.charge <= 0 {
        return Err(invalid("annotation requires positive peptide charge"));
    }
    // Scoped numeric validation and the subsequent sortedness scan both visit peaks.
    work.consume(peaks.len() * 2)?;
    for peak in peaks {
        if !peak.mz.is_finite() || peak.mz < 0.0 || !peak.intensity.is_finite() {
            return Err(invalid(
                "annotation peaks require finite nonnegative m/z and finite intensity",
            ));
        }
    }
    let theoretical = generator.generate_with_work(
        &hit.sequence,
        1,
        hit.charge.min(2) as u8,
        None,
        &mut work.generation,
    )?;
    work.allocate::<Peak1D>(peaks.len())?;
    let mut measured = MSSpectrum {
        peaks: peaks.to_vec(),
        ..Default::default()
    };
    if !measured.is_sorted() {
        work.sort::<Peak1D>(measured.len())?;
        measured
            .peaks
            .sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
    }
    // Every source alignment uses each theoretical peak at most once.
    work.allocate::<(usize, usize)>(theoretical.len())?;
    let pairs = alignment.align_with_work(&theoretical, &measured, &mut work.alignment)?;
    Ok(MatchedSpectrum {
        measured,
        theoretical,
        pairs,
    })
}

struct AnnotationWork {
    generation: TheoreticalGenerationWork,
    alignment: AlignmentWork,
    remaining: usize,
    bytes: usize,
    labels: usize,
}
impl Default for AnnotationWork {
    fn default() -> Self {
        Self {
            generation: TheoreticalGenerationWork::default(),
            alignment: AlignmentWork::default(),
            remaining: MAX_ANNOTATION_WORK,
            bytes: MAX_ANNOTATION_BYTES,
            labels: MAX_ANNOTATION_LABEL_BYTES,
        }
    }
}
impl AnnotationWork {
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| invalid("annotation work limit exceeded"))?;
        Ok(())
    }
    fn allocate<T>(&mut self, count: usize) -> Result<()> {
        let bytes = count
            .checked_mul(size_of::<T>())
            .ok_or_else(|| invalid("annotation allocation size overflow"))?;
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("annotation allocation limit exceeded"))?;
        Ok(())
    }
    fn copy_label(&mut self, bytes: usize) -> Result<()> {
        self.consume(bytes)?;
        self.labels = self
            .labels
            .checked_sub(bytes)
            .ok_or_else(|| invalid("annotation label limit exceeded"))?;
        self.allocate::<u8>(bytes)
    }
    fn sort<T>(&mut self, count: usize) -> Result<()> {
        let factor = if count < 2 {
            1
        } else {
            usize::BITS as usize - (count - 1).leading_zeros() as usize + 1
        };
        self.consume(
            count
                .checked_mul(factor)
                .ok_or_else(|| invalid("annotation sort work overflow"))?,
        )?;
        // Stable sort scratch space is bounded by one element per input item.
        self.allocate::<T>(count)
    }
}
fn finite(value: f64, message: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(message))
    }
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
