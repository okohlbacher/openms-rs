// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! mzQC quality-control files (source `FORMAT/MzQCFile.h`).
//!
//! The source's one member, `MzQCFile::store`, computes a fixed set of
//! quality metrics from an experiment, an optional feature map and optional
//! identifications, and writes them with their controlled-vocabulary names
//! as an mzQC 1.0.0 JSON document. The metrics come from the QC classes in
//! [`crate::analysis::qc`]; the names from the PSI-MS vocabulary and the
//! pinned PSI QC vocabulary `resources/cv/qc-cv.obo`. The TOPP tool
//! QCCalculator writes its `mzQC` output through this class.
//!
//! This module is compiled with the `mzml` feature, which brings the SHA-1
//! implementation the input-file checksum needs; an mzQC document always
//! describes an mzML input. The API mapping, preserved source behaviour and
//! evidence are in `docs/MZQC_FILE_SUPPORT.md`.

use super::controlled_vocabulary::ControlledVocabulary;
use super::file_handler::FileHandler;
use super::nlohmann_json::Json;
use crate::analysis::qc::{
    FeatureSummary, IdentificationSummary, QCBase, Requires, SpectrumCount, Status, Tic, TicResult,
};
use crate::data_structures::DateTime;
use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::kernel::{FeatureMap, MSExperiment};
use crate::system::file::{absolute_path, basename};
use crate::{Error, Result};
use std::path::Path;
use std::sync::OnceLock;

/// The pinned PSI QC vocabulary, `share/OpenMS/CV/qc-cv.obo` at the core SDK
/// pin (data-version 0.1.2); see `resources/cv/NOTICE.md`.
const QC_CV: &[u8] = include_bytes!("../../resources/cv/qc-cv.obo");

/// The document inputs besides the identifications, source `store`'s
/// parameters `input_file` to `feature_map`.
#[derive(Clone, Copy, Debug)]
pub struct MzQcInput<'a> {
    /// The mzML input file's name, as given: written as the input file's
    /// location (made absolute) and name, and hashed. Empty when there is no
    /// input file, which also withholds the metrics that need raw data.
    pub input_file: &'a str,
    /// The experiment to extract QC data from. The source requires
    /// `sortSpectra()` and `updateRanges()` beforehand; ranges are computed
    /// here on demand, so only the sort remains the caller's.
    pub experiment: &'a MSExperiment,
    /// Name of the person creating the file; omitted when empty.
    pub contact_name: &'a str,
    /// Contact address (mail, e-mail or phone) of that person; omitted when
    /// empty.
    pub contact_address: &'a str,
    /// Description of and comments about the file's contents; omitted when
    /// empty.
    pub description: &'a str,
    /// A unique and informative label for the run.
    pub label: &'a str,
    /// The features of a featureXML file; an empty map withholds the feature
    /// metrics.
    pub feature_map: &'a FeatureMap,
}

/// What [`MzQcFile::store`] reports besides the file.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MzQcReport {
    /// The lines the source writes to its log and to standard output, in
    /// order, without line ends: a metric that cannot run for lack of an input
    /// (`OPENMS_LOG_WARN`), the missed-cleavage warnings of the identification
    /// summary, and an accession missing from the vocabularies (`cout`).
    pub log: Vec<String>,
}

/// File adapter for mzQC files, source `MzQCFile`.
///
/// The source class is stateless and only stores; there is no mzQC reader in
/// the source either.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MzQcFile;

impl MzQcFile {
    /// A new adapter, source `MzQCFile()`.
    pub fn new() -> Self {
        Self
    }

    /// Store QC data of `input`, `prot_ids` and `pep_ids` as an mzQC file in
    /// JSON format at `output_file`, source `store`.
    ///
    /// The document is [`MzQcFile::text`] with the current local time as its
    /// creation date. It is written as plain UTF-8 whatever the file name's
    /// extension, as the source's `std::ofstream` writes it.
    ///
    /// The identification summary annotates the first hit of each peptide
    /// identification with the meta value `missed_cleavages`, which is why
    /// `pep_ids` is mutable; the source takes both lists by non-const
    /// reference for the same reason.
    ///
    /// The source opens (and so truncates) the output before computing
    /// anything, and a failing metric leaves an empty file behind; this port
    /// computes the whole document first and publishes it atomically, so a
    /// failure leaves an existing file unchanged.
    ///
    /// # Errors
    ///
    /// As [`MzQcFile::text`], and [`Error::Io`] when the file cannot be
    /// written (the source: `Exception::UnableToCreateFile`).
    pub fn store(
        &self,
        output_file: impl AsRef<Path>,
        input: &MzQcInput<'_>,
        prot_ids: &[ProteinIdentification],
        pep_ids: &mut [PeptideIdentification],
    ) -> Result<MzQcReport> {
        let (text, report) = self.text(input, prot_ids, pep_ids, &DateTime::now())?;
        super::path_io::write_plain(output_file.as_ref(), |writer| {
            writer.write_all(text.as_bytes())?;
            Ok(())
        })?;
        Ok(report)
    }

    /// The mzQC document [`MzQcFile::store`] writes, with `creation_date` as
    /// its creation date, and the report.
    ///
    /// The metrics a document holds depend on the inputs present, as the
    /// source's `QCBase::Status` decides: the raw-data metrics need a non-empty
    /// `input_file`, the feature metrics a non-empty feature map, and the
    /// identification metrics non-empty `prot_ids` *and* `pep_ids`. In order:
    ///
    /// | Accession | Metric | Needs |
    /// |---|---|---|
    /// | QC:4000059, QC:4000060 | number of MS1 and MS2 spectra | input file |
    /// | QC:4000135 | number of chromatograms | - |
    /// | QC:4000053 | RT duration: max - min spectrum RT, truncated to an integer | - |
    /// | QC:4000138 | m/z acquisition range: min and max peak m/z, truncated | - |
    /// | QC:4000067, QC:4000077 | TIC chromatogram over all MS levels and its area | input file, spectra |
    /// | QC:4000069, QC:4000172, QC:4000173 | MS1 TIC chromatogram, 10x jumps, 10x falls | input file, MS1 spectra |
    /// | QC:4000070 | MS2 TIC chromatogram | input file, MS2 spectra |
    /// | QC:4000257, QC:4000262 | detected compounds, RT mean shift when not 0 | feature map |
    /// | QC:4000186, QC:4000187, QC:4000214, QC:4000209, QC:4000185, QC:4000204 | PSMs, identified peptides and FDR threshold, mean peptide length, mean missed cleavages, identified proteins and FDR threshold, mean protein score | identifications |
    ///
    /// A chromatogram's value is an object with the relative intensities and
    /// retention times of [`TicResult`]. RT duration and m/z range are written
    /// whatever the inputs; without spectra (or without peaks, for the m/z
    /// range) the source converts an infinite or `DBL_MAX` bound to an
    /// unsigned integer, which is undefined, while this port writes 0 for the
    /// duration and saturates the bounds, giving `[4294967295, 0]`.
    ///
    /// The run's metadata holds the label, one input file (absolute location,
    /// base name, the format accession `MS:10000584` - sic, one digit too many,
    /// as the source writes it - with the name `mzML format`, and three file
    /// properties: the experiment's date and time as the completion time, the
    /// SHA-1 of the file ([`FileHandler::compute_file_hash`], empty when it
    /// cannot be read) and the instrument's name as the instrument model), and
    /// the analysis software QCCalculator with the accession `MS:1009001` and
    /// the core SDK version ([`crate::CORE_SDK_VERSION`]; the source writes its
    /// library version). The document's controlled vocabularies are always the
    /// PSI QC ontology 0.1.2 and PSI-MS 4.1.155.
    ///
    /// Text follows `nlohmann::ordered_json::dump(2)`: members in insertion
    /// order, two-space indentation, doubles in the library's shortest form
    /// (a `float` widened to `double` first), and no final newline.
    ///
    /// # Errors
    ///
    /// Returns the metric's error when a feature's `rt_deviation` is not
    /// numeric ([`FeatureSummary::compute`]) or the identifications name no
    /// usable enzyme ([`IdentificationSummary::compute`]), both of which the
    /// source throws; [`MSExperiment::ranges`]'s error for a spectrum that
    /// fails validation, which the source does not check; the vocabulary's
    /// error if the embedded vocabularies fail to load; and
    /// [`absolute_path`]'s error when the current directory cannot be resolved
    /// for the input file's location.
    pub fn text(
        &self,
        input: &MzQcInput<'_>,
        prot_ids: &[ProteinIdentification],
        pep_ids: &mut [PeptideIdentification],
        creation_date: &DateTime,
    ) -> Result<(String, MzQcReport)> {
        let cv = vocabulary()?;
        let exp = input.experiment;
        let mut report = MzQcReport::default();

        let mut status = Status::new();
        if !input.input_file.is_empty() {
            status |= Requires::RawMzml;
        }
        if !input.feature_map.features.is_empty() {
            status |= Requires::PreFdrFeat;
        }
        if !prot_ids.is_empty() && !pep_ids.is_empty() {
            status |= Requires::Id;
        }

        let mut metrics = Metrics {
            cv,
            list: Vec::new(),
            log: &mut report.log,
        };

        let spectrum_count = SpectrumCount::new();
        if runnable(&spectrum_count, &status, metrics.log) {
            let counts = spectrum_count.compute(exp);
            let count =
                |level: usize| Json::Unsigned(counts.get(&level).copied().unwrap_or(0).into());
            // Number of MS1 spectra, number of MS2 spectra.
            metrics.add("QC:4000059", count(1));
            metrics.add("QC:4000060", count(2));
        }
        // Number of chromatograms.
        metrics.add("QC:4000135", Json::Unsigned(exp.chromatograms.len() as u64));
        let ranges = exp.ranges(0)?;
        // Run time (RT duration).
        let duration = ranges.rt.map_or(0, |rt| (rt.max - rt.min) as u32);
        metrics.add("QC:4000053", Json::Unsigned(duration.into()));
        // MZ acquisition range.
        let (min_mz, max_mz) = ranges
            .mz
            .map_or((u32::MAX, 0), |mz| (mz.min as u32, mz.max as u32));
        metrics.add(
            "QC:4000138",
            Json::Array(vec![
                Json::Unsigned(min_mz.into()),
                Json::Unsigned(max_mz.into()),
            ]),
        );

        let tic = Tic::new();
        if runnable(&tic, &status, metrics.log) {
            // Complete TIC (all MS levels) with its area.
            let result = tic.compute(exp, 0.0, 0)?;
            if !result.intensities.is_empty() {
                metrics.add("QC:4000067", chromatogram(&result));
                metrics.add("QC:4000077", Json::Unsigned(result.area));
            }
            // MS1 TIC with the signal jump and fall counts.
            let result = tic.compute(exp, 0.0, 1)?;
            if !result.intensities.is_empty() {
                metrics.add("QC:4000069", chromatogram(&result));
                metrics.add("QC:4000172", Json::Unsigned(result.jump.into()));
                metrics.add("QC:4000173", Json::Unsigned(result.fall.into()));
            }
            // MS2 TIC.
            let result = tic.compute(exp, 0.0, 2)?;
            if !result.intensities.is_empty() {
                metrics.add("QC:4000070", chromatogram(&result));
            }
        }

        // Metabolomics: detected compounds from the featureXML file.
        let feature_summary = FeatureSummary::new();
        if runnable(&feature_summary, &status, metrics.log) {
            let result = feature_summary.compute(input.feature_map)?;
            metrics.add("QC:4000257", Json::Unsigned(result.feature_count.into()));
            if result.rt_shift_mean != 0.0 {
                metrics.add("QC:4000262", Json::Float(result.rt_shift_mean.into()));
            }
        }

        // Peptides and proteins from the idXML file.
        let identification_summary = IdentificationSummary::new();
        if runnable(&identification_summary, &status, metrics.log) {
            let result = identification_summary.compute(prot_ids, pep_ids)?;
            let unique = |id: crate::analysis::qc::UniqueId| {
                Json::Array(vec![
                    Json::Unsigned(id.count.into()),
                    Json::Float(id.fdr_threshold.into()),
                ])
            };
            metrics.add(
                "QC:4000186",
                Json::Unsigned(result.peptide_spectrum_matches.into()),
            );
            metrics.add("QC:4000187", unique(result.unique_peptides));
            metrics.add("QC:4000214", Json::Float(result.peptide_length_mean));
            metrics.add(
                "QC:4000209",
                Json::Float(result.missed_cleavages_mean.into()),
            );
            metrics.add("QC:4000185", unique(result.unique_proteins));
            metrics.add("QC:4000204", Json::Float(result.protein_hit_scores_mean));
        }
        let quality_metrics = metrics.into_json();

        let mut mzqc = vec![
            (
                "creationDate".to_owned(),
                Json::string(creation_date.iso_string()),
            ),
            ("version".to_owned(), Json::string("1.0.0")),
        ];
        for (name, value) in [
            ("contactName", input.contact_name),
            ("contactAddress", input.contact_address),
            ("description", input.description),
        ] {
            if !value.is_empty() {
                mzqc.push((name.to_owned(), Json::string(value)));
            }
        }
        let location = absolute_path(input.input_file)?;
        let settings = &exp.settings;
        let completion_time = format!(
            "{}T{}",
            settings.date_time.date_string(),
            settings.date_time.time_string()
        );
        let input_file = Json::object([
            ("location", Json::string(generic_string(&location))),
            ("name", Json::string(basename(input.input_file))),
            (
                "fileFormat",
                Json::object([
                    ("accession", Json::string("MS:10000584")),
                    ("name", Json::string("mzML format")),
                ]),
            ),
            (
                "fileProperties",
                Json::Array(vec![
                    property("MS:1000747", "completion time", completion_time),
                    property(
                        "MS:1000569",
                        "SHA-1",
                        FileHandler::compute_file_hash(input.input_file),
                    ),
                    property(
                        "MS:1000031",
                        "instrument model",
                        settings.instrument.name.clone(),
                    ),
                ]),
            ),
        ]);
        let software = Json::object([
            ("accession", Json::string("MS:1009001")),
            ("name", Json::string("QCCalculator")),
            ("version", Json::string(crate::CORE_SDK_VERSION)),
            ("uri", Json::string("https://www.openms.de")),
        ]);
        let run = Json::object([
            (
                "metadata",
                Json::object([
                    ("label", Json::string(input.label)),
                    ("inputFiles", Json::Array(vec![input_file])),
                    ("analysisSoftware", Json::Array(vec![software])),
                ]),
            ),
            ("qualityMetrics", quality_metrics),
        ]);
        mzqc.push(("runQualities".to_owned(), Json::Array(vec![run])));
        mzqc.push((
            "controlledVocabularies".to_owned(),
            Json::Array(vec![
                Json::object([
                    (
                        "name",
                        Json::string("Proteomics Standards Initiative Quality Control Ontology"),
                    ),
                    (
                        "uri",
                        Json::string(
                            "https://raw.githubusercontent.com/HUPO-PSI/mzQC/master/cv/qc-cv.obo",
                        ),
                    ),
                    ("version", Json::string("0.1.2")),
                ]),
                Json::object([
                    (
                        "name",
                        Json::string("Proteomics Standards Initiative Mass Spectrometry Ontology"),
                    ),
                    (
                        "uri",
                        Json::string("http://purl.obolibrary.org/obo/ms/psi-ms.obo"),
                    ),
                    ("version", Json::string("4.1.155")),
                ]),
            ]),
        ));
        let document = Json::object([("mzQC", Json::Object(mzqc))]);
        Ok((document.dump(2), report))
    }
}

/// The vocabularies the source loads for the metric names:
/// `loadFromOBO("PSI-MS", psi-ms.obo)` and then `loadFromOBO("QC", qc-cv.obo)`
/// into one vocabulary, built once from the embedded copies.
fn vocabulary() -> Result<&'static ControlledVocabulary> {
    static CV: OnceLock<std::result::Result<ControlledVocabulary, String>> = OnceLock::new();
    let built = CV.get_or_init(|| {
        let mut cv = ControlledVocabulary::psi_ms_only()
            .map_err(|e| e.to_string())?
            .clone();
        cv.load_obo_reader("QC", QC_CV).map_err(|e| e.to_string())?;
        Ok(cv)
    });
    match built {
        Ok(cv) => Ok(cv),
        Err(message) => Err(Error::InvalidValue(message.clone())),
    }
}

/// Source `isRunnable`, with its warnings appended to `log`.
fn runnable(metric: &dyn QCBase, status: &Status, log: &mut Vec<String>) -> bool {
    log.extend(metric.runnable_warnings(status));
    metric.is_runnable(status)
}

/// The quality metrics collected so far, source `quality_metrics` and its
/// `addMetric` lambda.
struct Metrics<'a> {
    cv: &'a ControlledVocabulary,
    list: Vec<Json>,
    log: &'a mut Vec<String>,
}

impl Metrics<'_> {
    /// Add one metric with its vocabulary name, or log that the accession is
    /// not in the vocabulary and skip it, as the source does.
    fn add(&mut self, accession: &str, value: Json) {
        match self.cv.get_term(accession) {
            Ok(term) => self.list.push(Json::object([
                ("accession", Json::string(accession)),
                ("name", Json::string(term.name.clone())),
                ("value", value),
            ])),
            Err(_) => self.log.push(format!("{accession} not found in CV.")),
        }
    }

    /// The `qualityMetrics` value: the list, or `null` when no metric was
    /// added, since the source's `json quality_metrics = {}` is a null value
    /// until its first `push_back`.
    fn into_json(self) -> Json {
        if self.list.is_empty() {
            Json::Null
        } else {
            Json::Array(self.list)
        }
    }
}

/// A TIC chromatogram metric value.
fn chromatogram(result: &TicResult) -> Json {
    let floats =
        |values: &[f32]| Json::Array(values.iter().map(|&v| Json::Float(v.into())).collect());
    Json::object([
        ("Relative intensity", floats(&result.relative_intensities)),
        ("Retention time", floats(&result.retention_times)),
    ])
}

/// One input-file property.
fn property(accession: &str, name: &str, value: String) -> Json {
    Json::object([
        ("accession", Json::string(accession)),
        ("name", Json::string(name)),
        ("value", Json::String(value)),
    ])
}

/// A path as `std::filesystem::path::generic_string` writes it: `/` as the
/// separator on every platform.
fn generic_string(path: &Path) -> String {
    let text = path.to_string_lossy();
    if cfg!(windows) {
        text.replace('\\', "/")
    } else {
        text.into_owned()
    }
}
