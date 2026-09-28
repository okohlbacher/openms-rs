// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Writer for the consensus-spectrum MGF file consumed by GNPS Feature-Based
//! and Ion-Identity Molecular Networking, ported from `FORMAT/GNPSMGFFile.h`.
//!
//! For each consensus feature of a [`ConsensusMap`], the GNPS workflow expects
//! one representative MS2 spectrum. This writer reads the spectrum references
//! stored on each consensus element, resolves them to the original mzML files,
//! optionally merges several MS2 scans per feature into one binned consensus
//! spectrum, and writes the result as an MGF file ready for upload to GNPS.
//! The `GNPSExport` TOPP tool uses it. See `docs/GNPS_MGF_FILE_SUPPORT.md`.

use crate::comparison::{
    BinConfig, BinUnit, BinnedSpectralContrastAngle, BinnedSpectrum, BinnedSpectrumCompareFunctor,
};
use crate::format::file_info::text_format::to_str;
use crate::kernel::on_disc_experiment::OnDiscMSExperiment;
use crate::kernel::{ConsensusFeature, ConsensusMap, MSSpectrum, Peak1D};
use crate::metadata::MetaValueData;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Name of the source `DefaultParamHandler`.
const HANDLER_NAME: &str = "GNPSMGFFile";
/// Default cosine-similarity threshold for merging MS2 scans into one
/// representative spectrum (source `DEF_COSINE_SIMILARITY`).
pub const DEF_COSINE_SIMILARITY: f64 = 0.9;
/// Default m/z bin width, the high-resolution default of `BinnedSpectrum`
/// (source `DEF_MERGE_BIN_SIZE`, the `float` 0.02 widened to `double`).
pub const DEF_MERGE_BIN_SIZE: f64 = 0.02_f32 as f64;
/// Default number of most intense peptides per consensus element for
/// `merged_spectra` output (source `DEF_PEPT_CUTOFF`).
pub const DEF_PEPT_CUTOFF: i32 = 5;
/// Declared by the source but wired to no parameter; kept for source
/// compatibility (source `DEF_MSMAP_CACHE`).
pub const DEF_MSMAP_CACHE: i32 = 50;
/// Bin offset of `BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES`.
const BIN_OFFSET_HIRES: f32 = 0.0;

/// What each spectrum block holds (parameter `output_type`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OutputType {
    /// The most intense MS2 spectrum of the feature (`most_intense`).
    #[default]
    MostIntense,
    /// The MS2 spectra similar to the most intense one, merged
    /// (`merged_spectra`).
    MergedSpectra,
}

impl OutputType {
    /// The parameter spelling.
    pub fn name(self) -> &'static str {
        match self {
            Self::MostIntense => "most_intense",
            Self::MergedSpectra => "merged_spectra",
        }
    }
}

/// Writer for the consensus-spectrum MGF of GNPS molecular networking (source
/// `GNPSMGFFile`, a `DefaultParamHandler` and `ProgressLogger`).
///
/// The parameters are `output_type` (`most_intense` or `merged_spectra`),
/// `peptide_cutoff` (how many of the most intense peptide annotations per
/// consensus element `merged_spectra` considers, -1 for all), `ms2_bin_size`
/// (bin width in Da for merging fragment ions) and
/// `merged_spectra:cos_similarity` (the cosine threshold for merging). The
/// source's progress reporting has no counterpart here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GNPSMGFFile {
    /// Parameter `output_type`.
    pub output_type: OutputType,
    /// Parameter `peptide_cutoff`; -1 considers every annotation.
    pub peptide_cutoff: i32,
    /// Parameter `ms2_bin_size`, in Da.
    pub ms2_bin_size: f64,
    /// Parameter `merged_spectra:cos_similarity`.
    pub cos_similarity: f64,
}

impl Default for GNPSMGFFile {
    fn default() -> Self {
        Self {
            output_type: OutputType::MostIntense,
            peptide_cutoff: DEF_PEPT_CUTOFF,
            ms2_bin_size: DEF_MERGE_BIN_SIZE,
            cos_similarity: DEF_COSINE_SIMILARITY,
        }
    }
}

fn float_entry(name: &str, value: f64, description: &str) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value: ParamValue::Float(value),
        min_float: 0.0,
        ..ParamEntry::default()
    }
}

impl GNPSMGFFile {
    /// A writer with the default parameters (source default constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// The default parameters (installed by the source constructor).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current parameters in the layout of [`Self::defaults`] (source
    /// `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut peptide_cutoff = ParamEntry {
            name: "peptide_cutoff".into(),
            description: "Number of most intense peptides to consider per consensus element; '-1' to consider all identifications.".into(),
            value: ParamValue::Integer(i64::from(self.peptide_cutoff)),
            ..ParamEntry::default()
        };
        peptide_cutoff.min_int = -1;
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                ParamEntry {
                    name: "output_type".into(),
                    description: "specificity of mgf output information".into(),
                    value: ParamValue::String(self.output_type.name().into()),
                    valid_strings: vec!["merged_spectra".into(), "most_intense".into()],
                    ..ParamEntry::default()
                },
                peptide_cutoff,
                float_entry(
                    "ms2_bin_size",
                    self.ms2_bin_size,
                    "Bin size (Da) for fragment ions when merging ms2 scans.",
                ),
            ],
            nodes: vec![ParamNode {
                name: "merged_spectra".into(),
                description:
                    "Options for exporting mgf file with merged spectra per consensusElement".into(),
                entries: vec![float_entry(
                    "cos_similarity",
                    self.cos_similarity,
                    "Cosine similarity threshold for merged_spectra output.",
                )],
                nodes: Vec::new(),
            }],
        })
    }

    /// A writer from parameters, as source `setParameters` followed by the
    /// parameter reads of `store`. Missing parameters take their defaults;
    /// unknown ones are ignored, as the source only warns (see
    /// [`Self::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction (the source's `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`Self::from_param`], also returning the unknown-parameter warnings
    /// the source logs.
    ///
    /// # Errors
    ///
    /// As [`Self::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(HANDLER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            let output_type = match param.value("output_type")?.as_str()? {
                "merged_spectra" => OutputType::MergedSpectra,
                _ => OutputType::MostIntense,
            };
            let peptide_cutoff = match param.value("peptide_cutoff")? {
                ParamValue::Integer(value) => i32::try_from(*value)
                    .map_err(|_| Error::InvalidValue("peptide_cutoff out of range".into()))?,
                _ => {
                    return Err(Error::InvalidValue(
                        "peptide_cutoff must be an integer".into(),
                    ));
                }
            };
            let float = |key: &str| -> Result<f64> {
                match param.value(key)? {
                    ParamValue::Float(value) => Ok(*value),
                    ParamValue::Integer(value) => Ok(*value as f64),
                    _ => Err(Error::InvalidValue(format!("{key} must be a number"))),
                }
            };
            Ok(Self {
                output_type,
                peptide_cutoff,
                ms2_bin_size: float("ms2_bin_size")?,
                cos_similarity: float("merged_spectra:cos_similarity")?,
            })
        })
    }

    /// Write the consensus-spectrum MGF for GNPS molecular networking (source
    /// `store`).
    ///
    /// The consensusXML at `consensus_file_path` must reference MS2 spectra in
    /// `mzml_file_paths`, in the order of its map indices (typically the mzML
    /// files that fed the feature finder), and those must be indexed mzML. The
    /// layout is that of [`Self::write`].
    ///
    /// # Errors
    ///
    /// The errors of the consensusXML reader and of [`Self::write`], and
    /// [`Error::Io`] when `out` cannot be created or written; the source
    /// writes to an unchecked stream.
    pub fn store(
        &self,
        consensus_file_path: impl AsRef<Path>,
        mzml_file_paths: &[impl AsRef<Path>],
        out: impl AsRef<Path>,
    ) -> Result<()> {
        let consensus_map = crate::format::consensusxml::load(consensus_file_path)?;
        let mut writer = BufWriter::new(File::create(out.as_ref())?);
        self.write(&consensus_map, mzml_file_paths, &mut writer)?;
        writer.flush()?;
        Ok(())
    }

    /// Write the MGF of `consensus_map` to `writer`, reading the MS2 spectra
    /// from `mzml_file_paths` on demand.
    ///
    /// For each consensus feature, in map order:
    ///
    /// 1. The charge is the largest of the feature's charge and its handles'.
    /// 2. The handles are ordered by intensity, most intense first, and each
    ///    handle's map contributes the first peptide identification of the
    ///    feature whose `map_index` meta value names that map and which has a
    ///    `spectrum_index` meta value. `most_intense` keeps the first such
    ///    annotation; `merged_spectra` keeps `peptide_cutoff` of them (all for
    ///    -1).
    /// 3. The best spectrum is that of the first annotation. A best spectrum
    ///    without peaks skips the feature (some Bruker files have such MS2
    ///    spectra).
    /// 4. The block header is `BEGIN IONS`, `OUTPUT=<output_type>`,
    ///    `SCANS=<feature position + 1>`, `FEATURE_ID=e_<unique id>`,
    ///    `MSLEVEL=2`, `CHARGE=<|charge|, 1 for 0><+ or ->`,
    ///    `PEPMASS=<feature m/z>`, `FILE_INDEX=<best spectrum index>` and
    ///    `RTINSECONDS=<best spectrum RT>`, numbers in the source's
    ///    `StringUtils::toStr` form.
    /// 5. With `merged_spectra`, every kept annotation's spectrum (the best
    ///    one included, so it is taken twice) whose binned cosine with the best
    ///    spectrum reaches `cos_similarity` is added, with bins of
    ///    `ms2_bin_size` Da, spread 1 and offset 0.
    /// 6. All peaks are sorted by m/z and binned: a bin closes when a peak
    ///    lies more than `ms2_bin_size` above the bin's first m/z, and becomes
    ///    one peak at the mean m/z and mean intensity; a closed bin whose
    ///    intensity sum is not positive is dropped (the last bin never is).
    ///    Each peak is written as `<mz>\t<intensity>` with four decimals,
    ///    followed by `END IONS` and an empty line.
    ///
    /// The source sorts handles and peaks with `std::sort`, which leaves the
    /// order of equal intensities, and of equal m/z values, unspecified; this
    /// port keeps their input order.
    ///
    /// # Errors
    ///
    /// - [`Error::MissingInformation`] when a feature has no usable peptide
    ///   annotation (or `peptide_cutoff` is 0), where the source reads past
    ///   an empty vector.
    /// - [`Error::InvalidValue`] when an annotation's `map_index` names no
    ///   entry of `mzml_file_paths` or its `spectrum_index` is negative or not
    ///   below the file's spectrum count (the source indexes without checks),
    ///   or an annotation meta value is not an integer (the source's
    ///   `DataValue` conversion throws).
    /// - [`Error::Parse`] when an mzML file is not indexed, and the errors of
    ///   the on-disc reader and of the binned comparison.
    /// - [`Error::Io`] when writing fails.
    pub fn write(
        &self,
        consensus_map: &ConsensusMap,
        mzml_file_paths: &[impl AsRef<Path>],
        mut writer: impl Write,
    ) -> Result<()> {
        let merged = self.output_type == OutputType::MergedSpectra;
        let pept_cutoff = if merged { self.peptide_cutoff } else { 1 };
        let mut files: BTreeMap<usize, OnDiscMSExperiment> = BTreeMap::new();
        let contrast = BinnedSpectralContrastAngle::new()?;
        let bins = BinConfig {
            size: self.ms2_bin_size as f32,
            unit: BinUnit::Absolute,
            spread: 1,
            offset: BIN_OFFSET_HIRES,
            ..BinConfig::default()
        };

        for (cons_i, feature) in consensus_map.features.iter().enumerate() {
            // Determine the feature's charge as the maximum feature handle charge.
            let charge = feature
                .handles()
                .iter()
                .map(|handle| handle.charge)
                .fold(feature.base.charge, i32::max);

            // Most intense peptide identifications (by precursor intensity).
            let mut pepts = element_peptides(feature)?;
            // Discard poorer precursor spectra.
            if pept_cutoff != -1 {
                pepts.truncate(usize::try_from(pept_cutoff).unwrap_or(0));
            }
            let Some(&(best_mapi, best_speci)) = pepts.first() else {
                return Err(Error::MissingInformation(format!(
                    "consensus feature {cons_i} has no MS2 peptide annotation with map_index and spectrum_index"
                )));
            };

            // Open the on-disc experiments of every kept annotation.
            for &(map_index, _) in &pepts {
                if !files.contains_key(&map_index) {
                    let path = mzml_file_paths.get(map_index).ok_or_else(|| {
                        Error::InvalidValue(format!(
                            "map_index {map_index} names no input mzML file ({} given)",
                            mzml_file_paths.len()
                        ))
                    })?;
                    files.insert(map_index, OnDiscMSExperiment::open(path)?);
                }
            }

            let best_spec = spectrum(&mut files, best_mapi, best_speci)?;
            if best_spec.peaks.is_empty() {
                continue;
            }

            write!(
                writer,
                "BEGIN IONS\nOUTPUT={}\nSCANS={}\nFEATURE_ID=e_{}\nMSLEVEL=2\nCHARGE={}{}\nPEPMASS={}\nFILE_INDEX={}\nRTINSECONDS={}\n",
                self.output_type.name(),
                cons_i + 1,
                feature.base.unique_id,
                if charge == 0 {
                    1
                } else {
                    charge.unsigned_abs()
                },
                if charge >= 0 { "+" } else { "-" },
                to_str(feature.base.mz),
                best_speci,
                to_str(best_spec.rt),
            )?;

            let mut spectra = vec![best_spec];
            if merged {
                let binned_best = BinnedSpectrum::new(&spectra[0], bins)?;
                for &(map_index, spec_index) in &pepts {
                    let test_spec = spectrum(&mut files, map_index, spec_index)?;
                    let binned = BinnedSpectrum::new(&test_spec, bins)?;
                    if contrast.score(&binned_best, &binned)? >= self.cos_similarity {
                        spectra.push(test_spec);
                    }
                }
            }

            for peak in flatten_and_bin(&spectra, self.ms2_bin_size) {
                writeln!(writer, "{:.4}\t{:.4}", peak.mz, f64::from(peak.intensity))?;
            }
            writer.write_all(b"END IONS\n\n")?;
        }
        Ok(())
    }
}

/// An integer meta value of a peptide identification; the source's `(int)`
/// conversion of a `DataValue` throws `Exception::ConversionError` for any
/// other type.
fn int_meta(value: &crate::metadata::MetaValue, key: &str) -> Result<i64> {
    match value.data() {
        MetaValueData::Integer(v) => Ok(*v),
        _ => Err(Error::InvalidValue(format!(
            "peptide identification meta value '{key}' is not an integer"
        ))),
    }
}

/// The source's `sortElementMapsByIntensity_` followed by
/// `getElementPeptideIdentificationsByElementIntensity_`: `(map index,
/// spectrum index)` per handle, most intense handle first.
fn element_peptides(feature: &ConsensusFeature) -> Result<Vec<(usize, usize)>> {
    let mut element_maps: Vec<(u64, f32)> = feature
        .handles()
        .iter()
        .map(|handle| (handle.map_index, handle.intensity))
        .collect();
    element_maps.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut pepts = Vec::new();
    for (element_map, _) in element_maps {
        for identification in &feature.base.peptide_identifications {
            let (Some(spectrum_index), Some(map_index)) = (
                identification.metadata.get("spectrum_index"),
                identification.metadata.get("map_index"),
            ) else {
                continue;
            };
            let map_index = int_meta(map_index, "map_index")?;
            if u64::try_from(map_index).ok() != Some(element_map) {
                continue;
            }
            let spec_index = int_meta(spectrum_index, "spectrum_index")?;
            let spec_index = usize::try_from(spec_index).map_err(|_| {
                Error::InvalidValue(format!("negative spectrum_index {spec_index}"))
            })?;
            pepts.push((map_index as usize, spec_index));
            break;
        }
    }
    Ok(pepts)
}

fn spectrum(
    files: &mut BTreeMap<usize, OnDiscMSExperiment>,
    map_index: usize,
    spec_index: usize,
) -> Result<MSSpectrum> {
    let file = files
        .get_mut(&map_index)
        .ok_or_else(|| Error::InvalidValue(format!("map_index {map_index} is not open")))?;
    if spec_index >= file.spectrum_count() {
        return Err(Error::InvalidValue(format!(
            "spectrum_index {spec_index} is not below the {} spectra of map {map_index}",
            file.spectrum_count()
        )));
    }
    file.spectrum(spec_index)
}

/// The source's `flattenAndBinSpectra_` and `binPeaks_`.
fn flatten_and_bin(spectra: &[MSSpectrum], bin_width: f64) -> Vec<Peak1D> {
    let mut flat: Vec<Peak1D> = spectra
        .iter()
        .flat_map(|s| s.peaks.iter().copied())
        .collect();
    flat.sort_by(|a, b| a.mz.total_cmp(&b.mz));
    let mut binned = Vec::new();
    let Some(first) = flat.first() else {
        return binned;
    };
    let mut last_mz = first.mz;
    let mut sum_mz = 0.0_f64;
    let mut sum_intensity = 0.0_f64;
    let mut count = 0usize;
    for peak in &flat {
        if count > 0 && peak.mz - last_mz > bin_width {
            if sum_intensity > 0.0 {
                binned.push(Peak1D::new(
                    sum_mz / count as f64,
                    (sum_intensity / count as f64) as f32,
                ));
            }
            last_mz = peak.mz;
            sum_mz = 0.0;
            sum_intensity = 0.0;
            count = 0;
        }
        sum_mz += peak.mz;
        sum_intensity += f64::from(peak.intensity);
        count += 1;
    }
    if count > 0 {
        binned.push(Peak1D::new(
            sum_mz / count as f64,
            (sum_intensity / count as f64) as f32,
        ));
    }
    binned
}
