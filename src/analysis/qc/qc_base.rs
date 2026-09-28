// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The base of every quality-control metric (source `QC/QCBase.h`).
//!
//! A metric declares the input data its computation needs as a [`Status`], a
//! [`FlagSet`] of [`Requires`] values; a caller collects the inputs it has into
//! another `Status` and asks [`QCBase::is_runnable`]. The source's abstract
//! class is the [`QCBase`] trait; its static members are the associated
//! functions and constants of this module. See `docs/QC_SUPPORT.md`.

use crate::data_structures::flag_set::{FlagEnum, FlagSet};
use crate::kernel::{ConsensusMap, FeatureMap, MSExperiment};
use crate::{Error, Result};
use std::collections::BTreeMap;

/// One kind of input data a metric can require, encoded as a bit.
///
/// Source `QCBase::Requires`, an `enum class` over `UInt64` whose values are
/// the bit positions `0..=6`. `Nothing` is position 0 and therefore a real
/// bit: a metric requiring `Nothing` is runnable only when the caller's
/// status also holds `Nothing`, which the empty status does not. The source
/// behaves the same; no metric in either code base requires `Nothing`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum Requires {
    /// Default, does not require anything.
    Nothing = 0,
    /// An mzML file is required.
    RawMzml = 1,
    /// Features with FDR-filtered peptide identifications.
    PostFdrFeat = 2,
    /// Features with unfiltered peptide identifications.
    PreFdrFeat = 3,
    /// A contaminant database.
    Contaminants = 4,
    /// transformationXML files for RT alignment.
    TrafoAlign = 5,
    /// An idXML file with protein identifications.
    Id = 6,
}

impl Requires {
    /// Source `Requires::SIZE_OF_REQUIRES`: the number of values.
    pub const SIZE_OF_REQUIRES: usize = 7;
    /// Every value in bit order.
    pub const ALL: [Requires; Requires::SIZE_OF_REQUIRES] = [
        Requires::Nothing,
        Requires::RawMzml,
        Requires::PostFdrFeat,
        Requires::PreFdrFeat,
        Requires::Contaminants,
        Requires::TrafoAlign,
        Requires::Id,
    ];
    /// The input file each value stands for, source
    /// `QCBase::names_of_requires[value]`.
    pub fn name(self) -> &'static str {
        NAMES_OF_REQUIRES[self as usize]
    }
}

impl FlagEnum for Requires {
    fn bit(self) -> u32 {
        u32::from(self as u8)
    }
}

/// Source `QCBase::names_of_requires`: one entry per [`Requires`] value, in
/// bit order. The source's class test guards this length because an earlier
/// version omitted `id.idXML` and `isRunnable` read past the array.
pub const NAMES_OF_REQUIRES: [&str; Requires::SIZE_OF_REQUIRES] = [
    "fail",
    "raw.mzML",
    "postFDR.featureXML",
    "preFDR.featureXML",
    "contaminants.fasta",
    "trafoAlign.trafoXML",
    "id.idXML",
];

/// The unit of a mass tolerance, source `QCBase::ToleranceUnit`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum ToleranceUnit {
    /// Chosen automatically.
    Auto = 0,
    /// Parts per million.
    Ppm = 1,
    /// Dalton.
    Da = 2,
}

impl ToleranceUnit {
    /// Source `ToleranceUnit::SIZE_OF_TOLERANCEUNIT`: the number of values.
    pub const SIZE_OF_TOLERANCE_UNIT: usize = 3;
    /// The unit's name, source `QCBase::names_of_toleranceUnit[value]`.
    pub fn name(self) -> &'static str {
        NAMES_OF_TOLERANCE_UNIT[self as usize]
    }
}

/// Source `QCBase::names_of_toleranceUnit`, one entry per [`ToleranceUnit`].
pub const NAMES_OF_TOLERANCE_UNIT: [&str; ToleranceUnit::SIZE_OF_TOLERANCE_UNIT] =
    ["auto", "ppm", "da"];

/// The input data a metric requires or a caller provides, source
/// `QCBase::Status = FlagSet<Requires>`.
pub type Status = FlagSet<Requires>;

/// Map from a spectrum's native ID to its index in an experiment.
///
/// Source `QCBase::SpectraMap`. When several spectra share a native ID the
/// last one wins, as the source's repeated `operator[]` assignment does.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SpectraMap {
    nativeid_to_index: BTreeMap<String, u64>,
}

impl SpectraMap {
    /// An empty map, source `SpectraMap()`.
    pub fn new() -> Self {
        Self::default()
    }
    /// A map indexing `exp` at once, source `explicit SpectraMap(const MSExperiment&)`.
    pub fn from_experiment(exp: &MSExperiment) -> Self {
        let mut map = Self::new();
        map.calculate_map(exp);
        map
    }
    /// Index `exp`, replacing whatever the map held, source `calculateMap`.
    pub fn calculate_map(&mut self, exp: &MSExperiment) {
        self.nativeid_to_index.clear();
        for (index, spectrum) in exp.spectra.iter().enumerate() {
            self.nativeid_to_index
                .insert(spectrum.native_id.clone(), index as u64);
        }
    }
    /// The index of the spectrum with native ID `identifier`, source `at`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] with the source's message when no
    /// spectrum has that ID; the source throws `Exception::ElementNotFound`.
    pub fn at(&self, identifier: &str) -> Result<u64> {
        self.nativeid_to_index
            .get(identifier)
            .copied()
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "No spectrum with identifier '{identifier}' in MSExperiment!"
                ))
            })
    }
    /// Remove every entry, source `clear`.
    pub fn clear(&mut self) {
        self.nativeid_to_index.clear();
    }
    /// Whether the map is empty, source `empty`.
    pub fn is_empty(&self) -> bool {
        self.nativeid_to_index.is_empty()
    }
    /// The number of distinct native IDs, source `size`.
    pub fn len(&self) -> usize {
        self.nativeid_to_index.len()
    }
}

/// A quality-control metric: its name and the input data it requires.
///
/// Source abstract class `QCBase`. Its two pure virtual members are the
/// required methods; `isRunnable` is provided.
pub trait QCBase {
    /// The name of the metric, source `getName`.
    fn name(&self) -> &str;
    /// The input data the metric's computation requires, source
    /// `requirements`.
    fn requirements(&self) -> Status;

    /// Whether `status` provides every input the metric requires, source
    /// `isRunnable`.
    ///
    /// The source also writes one `OPENMS_LOG_WARN` line per missing input;
    /// [`QCBase::runnable_warnings`] returns those lines instead of logging
    /// them.
    fn is_runnable(&self, status: &Status) -> bool {
        status.is_super_set_of(&self.requirements())
    }

    /// The warnings the source's `isRunnable` logs for `status`: none when the
    /// metric is runnable, otherwise one line, in bit order, per required input
    /// that `status` lacks, in the source's words.
    fn runnable_warnings(&self, status: &Status) -> Vec<String> {
        if self.is_runnable(status) {
            return Vec::new();
        }
        let required = self.requirements();
        Requires::ALL
            .iter()
            .filter(|&&input| {
                required.is_super_set_of_flag(input) && !status.is_super_set_of_flag(input)
            })
            .map(|input| {
                format!(
                    "Note: Metric '{}' cannot run because input data '{}' is missing!",
                    self.name(),
                    input.name()
                )
            })
            .collect()
    }
}

/// Whether the IsobaricAnalyzer TOPP tool created `cm`, source static
/// `QCBase::isLabeledExperiment`, which asks
/// `DataProcessingUtils::hasIsobaricAnalyzer`: some data-processing record's
/// software is named exactly `IsobaricAnalyzer`.
pub fn is_labeled_experiment(cm: &ConsensusMap) -> bool {
    crate::kernel::map_operations::has_isobaric_analyzer(&cm.data_processing)
}

/// Whether `fmap` holds a peptide identification, either unassigned or on a
/// top-level feature, source template `QCBase::hasPepID` for a feature map.
/// Subordinate features are not inspected, as in the source.
pub fn has_pep_id_feature_map(fmap: &FeatureMap) -> bool {
    !fmap.unassigned_peptide_identifications.is_empty()
        || fmap
            .features
            .iter()
            .any(|f| !f.peptide_identifications.is_empty())
}

/// Whether `cmap` holds a peptide identification, either unassigned or on a
/// consensus feature, source template `QCBase::hasPepID` for a consensus map.
pub fn has_pep_id_consensus_map(cmap: &ConsensusMap) -> bool {
    !cmap.unassigned_peptide_identifications.is_empty()
        || cmap
            .features
            .iter()
            .any(|f| !f.peptide_identifications.is_empty())
}
