// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Comparators, hashes and stream output of peptide and protein hits.
//!
//! Completes the ports of `src/openms/include/OpenMS/METADATA/PeptideHit.h`
//! and `src/openms/include/OpenMS/METADATA/ProteinHit.h` (core SDK `bc9cc12`)
//! for the members that are not record fields: the comparator functors, the
//! hash functors and `std::hash` specialisations, `operator<<`, the fragment
//! annotation writer and `ProteinHit::operator=(const MetaInfoInterface&)`.
//! The records themselves are [`PeptideHit`], [`ProteinHit`] and
//! [`PeakAnnotation`]; see `docs/IDENTIFICATION_SUPPORT.md`.
//!
//! # Hashing
//!
//! Every `Hash` here is consistent with the record's derived equality: equal
//! values hash equally, floating values with both zero signs identified as
//! `==` identifies them. Hash *values* are not the source's FNV-1a bits, which
//! no caller can observe except through a C++ hash table's iteration order.
//! Where the source hashes a field through a type this crate does not hash
//! (an `AASequence`, a `ResidueModification`), the sequence's canonical text or
//! the modification's position stands in: equal records still hash equally.
//! The records hold floating values and are therefore not `Eq`, so a hash-keyed
//! container of whole hits needs a key type; [`SequenceChargeKey`] is the one
//! the source itself provides.

use std::fmt;
use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher, Hash, Hasher};

use super::{AnalysisResult, PeakAnnotation, PeptideHit, PeptideIdentification, ProteinHit};
use crate::Result;
use crate::metadata::MetaInfo;
use crate::param::value::format_float;

/// The source's forward-looking name for a peptide hit, `SpectrumMatch`
/// (`PeptideHit.h`: `using SpectrumMatch = PeptideHit;`).
pub type SpectrumMatch = PeptideHit;

/// `hash_float`: both zero signs hash alike, as `==` identifies them.
fn hash_f64<H: Hasher>(value: f64, state: &mut H) {
    (if value == 0.0 { 0 } else { value.to_bits() }).hash(state);
}

impl Hash for PeakAnnotation {
    /// The source's `std::hash<PeptideHit::PeakAnnotation>` covers the same
    /// four fields.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.annotation.hash(state);
        self.charge.hash(state);
        hash_f64(self.mz, state);
        hash_f64(self.intensity, state);
    }
}

impl Hash for AnalysisResult {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.score_type.hash(state);
        self.higher_is_better.hash(state);
        hash_f64(self.main_score, state);
        self.sub_scores.len().hash(state);
        for (name, score) in &self.sub_scores {
            name.hash(state);
            hash_f64(*score, state);
        }
    }
}

impl PeakAnnotation {
    /// The source's `PeakAnnotation::operator<`: order by m/z, then charge,
    /// then annotation text, then intensity.
    ///
    /// Floating values compare with `f64::total_cmp`, which agrees with the
    /// source's `<` on every finite value except that `-0.0` sorts before
    /// `0.0`; [`PeakAnnotation::validate`] refuses the non-finite values the
    /// source's comparison leaves undefined.
    pub fn cmp_source_order(&self, other: &Self) -> std::cmp::Ordering {
        self.mz
            .total_cmp(&other.mz)
            .then(self.charge.cmp(&other.charge))
            .then_with(|| self.annotation.cmp(&other.annotation))
            .then(self.intensity.total_cmp(&other.intensity))
    }

    /// Fragment annotations as the `fragment_annotation` text idXML stores,
    /// the source's `PeakAnnotation::writePeakAnnotationsString_`.
    ///
    /// The annotations are stably sorted by [`PeakAnnotation::cmp_source_order`]
    /// and written as `mz,intensity,charge,"annotation"` joined by `|`; the
    /// numbers use `StringUtils::toStr(double)` (fixed 15 fractional digits,
    /// trailing zeros removed, scientific outside `[1e-2, 1e4)`) and the
    /// annotation is quoted with backslash escapes as `StringUtils::quote`
    /// does. An empty list gives an empty string. The source appends to an
    /// out-parameter; this returns the text to append.
    ///
    /// # Errors
    ///
    /// [`crate::Error::InvalidValue`] when an annotation fails
    /// [`PeakAnnotation::validate`] (a non-finite m/z or intensity), where the
    /// source's sort is undefined.
    pub fn write_peak_annotations_string(annotations: &[PeakAnnotation]) -> Result<String> {
        for annotation in annotations {
            annotation.validate()?;
        }
        let mut sorted: Vec<&PeakAnnotation> = annotations.iter().collect();
        sorted.sort_by(|a, b| a.cmp_source_order(b));
        let mut out = String::new();
        for (i, annotation) in sorted.iter().enumerate() {
            if i > 0 {
                out.push('|');
            }
            out.push_str(&format_float(annotation.mz, true));
            out.push(',');
            out.push_str(&format_float(annotation.intensity, true));
            out.push(',');
            out.push_str(&annotation.charge.to_string());
            out.push(',');
            out.push('"');
            out.push_str(
                &annotation
                    .annotation
                    .replace('\\', "\\\\")
                    .replace('"', "\\\""),
            );
            out.push('"');
        }
        Ok(out)
    }
}

impl Hash for PeptideHit {
    /// The source's `std::hash<PeptideHit>` covers the metadata (rank
    /// included, which it stores as a meta value), sequence, score, charge,
    /// evidences and peak annotations. Here rank is a field and is hashed as
    /// one, and the analysis results are hashed too, because this crate's
    /// equality compares them. The sequence is hashed through its canonical
    /// text, modifications included.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.sequence.to_string().hash(state);
        hash_f64(self.score, state);
        self.rank.hash(state);
        self.charge.hash(state);
        self.evidences.hash(state);
        self.peak_annotations.hash(state);
        self.analysis_results.hash(state);
        self.metadata.hash(state);
    }
}

impl PeptideHit {
    /// The source's `PeptideHit::ScoreMore` predicate: `a.score > b.score`.
    pub fn score_more(a: &Self, b: &Self) -> bool {
        a.score > b.score
    }

    /// The source's `PeptideHit::ScoreLess` predicate: `a.score < b.score`.
    pub fn score_less(a: &Self, b: &Self) -> bool {
        a.score < b.score
    }

    /// The source's `PeptideHit::RankLess` predicate: `a.rank < b.rank`.
    pub fn rank_less(a: &Self, b: &Self) -> bool {
        a.rank < b.rank
    }

    /// The source's `PeptideHit::SequenceLessComparator`: whether `a`'s
    /// sequence text, modifications included, sorts before `b`'s.
    ///
    /// The source declares the functor's `operator()` without an access
    /// specifier in a `class`, so it is private and no caller can use it; this
    /// is what it was written to do.
    pub fn sequence_less(a: &Self, b: &Self) -> bool {
        a.sequence.to_string() < b.sequence.to_string()
    }
}

impl fmt::Display for PeptideHit {
    /// The source's `operator<<`: `peptide hit with sequence '<seq>', charge
    /// <z>, score <score>`, the score written by `StringUtils::toStr(double)`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "peptide hit with sequence '{}', charge {}, score {}",
            self.sequence,
            self.charge,
            format_float(self.score, true)
        )
    }
}

/// A peptide hit keyed by sequence and charge only, the source's
/// `PeptideHit::SequenceChargeHash` and `SequenceChargeEqual` pair.
///
/// Hash and equality both look at the modified sequence and the charge and
/// nothing else, which is the "identity" of a peptide hit for most practical
/// purposes; they are deliberately not consistent with the hit's own equality,
/// which also compares score, rank, evidences, annotations and metadata. The
/// source's `std::unordered_set<PeptideHit, SequenceChargeHash,
/// SequenceChargeEqual>` is `HashSet<SequenceChargeKey>` here. Equality is
/// [`PeptideHit::same_sequence_and_charge`].
#[derive(Clone, Copy, Debug)]
pub struct SequenceChargeKey<'a>(pub &'a PeptideHit);

impl PartialEq for SequenceChargeKey<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.0.same_sequence_and_charge(other.0)
    }
}

impl Eq for SequenceChargeKey<'_> {}

impl Hash for SequenceChargeKey<'_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.sequence.to_string().hash(state);
        self.0.charge.hash(state);
    }
}

impl Hash for ProteinHit {
    /// The source's `std::hash<ProteinHit>` hashes score, rank, accession,
    /// sequence, coverage and each modification's position and full ID. This
    /// hashes the same fields; a modification contributes its position only,
    /// which keeps the hash consistent with equality, and the metadata, which
    /// equality compares, is hashed too.
    fn hash<H: Hasher>(&self, state: &mut H) {
        hash_f64(self.score, state);
        self.rank.hash(state);
        self.accession.hash(state);
        self.sequence.hash(state);
        self.coverage.is_some().hash(state);
        if let Some(coverage) = self.coverage {
            hash_f64(coverage, state);
        }
        self.modifications.len().hash(state);
        for modification in &self.modifications {
            modification.position.hash(state);
        }
        self.metadata.hash(state);
    }
}

impl ProteinHit {
    /// The source's `ProteinHit::ScoreMore` predicate: `(a.score, a.accession)
    /// > (b.score, b.accession)`, so equal scores fall back to the accession.
    pub fn score_more(a: &Self, b: &Self) -> bool {
        (a.score, a.accession.as_str()) > (b.score, b.accession.as_str())
    }

    /// The source's `ProteinHit::ScoreLess` predicate: `(a.score, a.accession)
    /// < (b.score, b.accession)`.
    pub fn score_less(a: &Self, b: &Self) -> bool {
        (a.score, a.accession.as_str()) < (b.score, b.accession.as_str())
    }

    /// Hash of the accession alone, the source's `ProteinHitAccessionHash`;
    /// `ProteinHitPtrAccessionHash`, the same through a pointer, is this
    /// function on a reference.
    ///
    /// The hasher is deterministic (`DefaultHasher` with its fixed keys), so
    /// the value is stable within one build, as the source's `std::hash` is.
    pub fn accession_hash(&self) -> u64 {
        BuildHasherDefault::<DefaultHasher>::default().hash_one(self.accession.as_str())
    }

    /// Replace the metadata by a copy of `source`, leaving every other field
    /// alone: the source's `ProteinHit::operator=(const MetaInfoInterface&)`.
    pub fn assign_meta_info(&mut self, source: &MetaInfo) {
        self.metadata.clone_from(source);
    }
}

impl fmt::Display for ProteinHit {
    /// The source's `operator<<`: `protein hit with accession '<acc>', score
    /// <score>`, the score written by `StringUtils::toStr(double)`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "protein hit with accession '{}', score {}",
            self.accession,
            format_float(self.score, true)
        )
    }
}

impl Hash for PeptideIdentification {
    /// Consistent with the derived equality: every field, floating values with
    /// both zero signs alike. Lets a `SpectrumSettings`, which carries the
    /// identifications as a native attachment, hash consistently with its
    /// equality.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.identifier.hash(state);
        self.hits.hash(state);
        self.score_type.hash(state);
        self.higher_score_better.hash(state);
        hash_f64(self.significance_threshold, state);
        for coordinate in [self.rt, self.mz] {
            coordinate.is_some().hash(state);
            if let Some(value) = coordinate {
                hash_f64(value, state);
            }
        }
        self.metadata.hash(state);
    }
}
