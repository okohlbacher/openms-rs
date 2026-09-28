// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `ResidueModification::SourceClassification` of `CHEMISTRY/ResidueModification.h`.
//!
//! A [`ResidueModification`] keeps its classification as the text it was
//! loaded with (Unimod's `classification` attribute, or the compiled table's
//! column). The source converts that text to this enum when it is set and
//! converts back when it is read, so text it does not know reads back empty.
//!
//! # API mapping
//!
//! | Source (`ResidueModification.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `enum SourceClassification` (`ARTIFACT` 0 ... `UNKNOWN` 15) | [`SourceClassification`], [`SourceClassification::value`] |
//! | `NUMBER_OF_SOURCE_CLASSIFICATIONS` | [`SourceClassification::COUNT`] |
//! | `setSourceClassification(const std::string&)` (the text-to-enum rule) | [`SourceClassification::from_text`] |
//! | `getSourceClassification()` | [`ResidueModification::source_classification`] |
//! | `getSourceClassificationName(SourceClassification)` | [`SourceClassification::name`]; for a record, [`ResidueModification::source_classification_name`] |

use super::ResidueModification;

/// Source `ResidueModification::SourceClassification`: the biological or
/// technical origin of a modification.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SourceClassification {
    /// `ARTIFACT` (0): a technical or chemical artefact.
    Artifact,
    /// `HYPOTHETICAL` (1).
    Hypothetical,
    /// `NATURAL` (2).
    Natural,
    /// `POSTTRANSLATIONAL` (3).
    PostTranslational,
    /// `MULTIPLE` (4): a combination of multiple modifications.
    Multiple,
    /// `CHEMICAL_DERIVATIVE` (5).
    ChemicalDerivative,
    /// `ISOTOPIC_LABEL` (6): a chemical label (artificial).
    IsotopicLabel,
    /// `PRETRANSLATIONAL` (7).
    PreTranslational,
    /// `OTHER_GLYCOSYLATION` (8): neither N- nor O-linked.
    OtherGlycosylation,
    /// `NLINKED_GLYCOSYLATION` (9).
    NLinkedGlycosylation,
    /// `AA_SUBSTITUTION` (10): an amino-acid substitution that presents like a
    /// modification.
    AaSubstitution,
    /// `OTHER` (11).
    Other,
    /// `NONSTANDARD_RESIDUE` (12): a non-standard amino acid.
    NonstandardResidue,
    /// `COTRANSLATIONAL` (13).
    CoTranslational,
    /// `OLINKED_GLYCOSYLATION` (14).
    OLinkedGlycosylation,
    /// `UNKNOWN` (15), also for any text the source does not recognise.
    Unknown,
}

impl SourceClassification {
    /// Source `NUMBER_OF_SOURCE_CLASSIFICATIONS`.
    pub const COUNT: usize = 16;

    /// Every value in source order.
    pub const ALL: [Self; Self::COUNT] = [
        Self::Artifact,
        Self::Hypothetical,
        Self::Natural,
        Self::PostTranslational,
        Self::Multiple,
        Self::ChemicalDerivative,
        Self::IsotopicLabel,
        Self::PreTranslational,
        Self::OtherGlycosylation,
        Self::NLinkedGlycosylation,
        Self::AaSubstitution,
        Self::Other,
        Self::NonstandardResidue,
        Self::CoTranslational,
        Self::OLinkedGlycosylation,
        Self::Unknown,
    ];

    /// The source enum value.
    pub fn value(self) -> usize {
        self as usize
    }

    /// Source `setSourceClassification(const std::string&)`: the text
    /// compared ASCII-case-insensitively with the fifteen known names, where
    /// `artifact` and `artefact` are both [`Self::Artifact`]; anything else is
    /// [`Self::Unknown`].
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::chemistry::modification_classification::SourceClassification;
    ///
    /// assert_eq!(SourceClassification::from_text("Artefact"), SourceClassification::Artifact);
    /// assert_eq!(SourceClassification::from_text("chemical DERIVATIVE"), SourceClassification::ChemicalDerivative);
    /// assert_eq!(SourceClassification::from_text("Synth. pep. protect. gp."), SourceClassification::Unknown);
    /// ```
    pub fn from_text(text: &str) -> Self {
        let lower = text.to_ascii_lowercase();
        match lower.as_str() {
            "artifact" | "artefact" => Self::Artifact,
            "natural" => Self::Natural,
            "hypothetical" => Self::Hypothetical,
            "post-translational" => Self::PostTranslational,
            "multiple" => Self::Multiple,
            "chemical derivative" => Self::ChemicalDerivative,
            "isotopic label" => Self::IsotopicLabel,
            "pre-translational" => Self::PreTranslational,
            "other glycosylation" => Self::OtherGlycosylation,
            "n-linked glycosylation" => Self::NLinkedGlycosylation,
            "aa substitution" => Self::AaSubstitution,
            "other" => Self::Other,
            "non-standard residue" => Self::NonstandardResidue,
            "co-translational" => Self::CoTranslational,
            "o-linked glycosylation" => Self::OLinkedGlycosylation,
            _ => Self::Unknown,
        }
    }

    /// Source `getSourceClassificationName`: `Artefact` (British spelling),
    /// `Post-translational`, ..., and the empty string for [`Self::Unknown`].
    pub fn name(self) -> &'static str {
        match self {
            Self::Artifact => "Artefact",
            Self::Hypothetical => "Hypothetical",
            Self::Natural => "Natural",
            Self::PostTranslational => "Post-translational",
            Self::Multiple => "Multiple",
            Self::ChemicalDerivative => "Chemical derivative",
            Self::IsotopicLabel => "Isotopic label",
            Self::PreTranslational => "Pre-translational",
            Self::OtherGlycosylation => "Other glycosylation",
            Self::NLinkedGlycosylation => "N-linked glycosylation",
            Self::AaSubstitution => "AA substitution",
            Self::Other => "Other",
            Self::NonstandardResidue => "Non-standard residue",
            Self::CoTranslational => "Co-translational",
            Self::OLinkedGlycosylation => "O-linked glycosylation",
            Self::Unknown => "",
        }
    }
}

impl ResidueModification {
    /// The classification as the source enum: source
    /// `getSourceClassification()`, from the stored text by
    /// [`SourceClassification::from_text`].
    pub fn source_classification(&self) -> SourceClassification {
        SourceClassification::from_text(self.classification())
    }

    /// The classification name the source reports: source
    /// `getSourceClassificationName()`. Text the source does not recognise,
    /// such as Unimod's `Synth. pep. protect. gp.`, reads back empty, while
    /// [`ResidueModification::classification`] keeps it.
    pub fn source_classification_name(&self) -> &'static str {
        self.source_classification().name()
    }
}
