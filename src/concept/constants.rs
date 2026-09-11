// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Source SDK mathematical/physical values and metadata keys.
//!
//! Historical values, expressions and spellings are retained for compatibility;
//! these are not substituted with newer physical reference values. See
//! `docs/CONSTANTS_SUPPORT.md` for units and the mutable epsilon mapping.

// Original source precision and spelling make the public mapping auditable.
#![allow(
    clippy::approx_constant,
    clippy::excessive_precision,
    non_upper_case_globals
)]

use std::sync::atomic::{AtomicU64, Ordering};

pub const PI: f64 = 3.14159265358979323846;
pub const E: f64 = 2.718281828459045235;
pub const ELEMENTARY_CHARGE: f64 = 1.60217738E-19;
pub const e0: f64 = ELEMENTARY_CHARGE;
pub const ELECTRON_MASS: f64 = 9.1093897E-31;
pub const ELECTRON_MASS_U: f64 = 1.0 / 1822.8885020477;
pub const PROTON_MASS: f64 = 1.6726230E-27;
pub const PROTON_MASS_U: f64 = 1.0072764667710;
pub const C13C12_MASSDIFF_U: f64 = 1.0033548378;
pub const ISOTOPE_MASSDIFF_55K_U: f64 = 1.002371;
pub const NEUTRON_MASS: f64 = 1.6749286E-27;
pub const NEUTRON_MASS_U: f64 = 1.00866491566;
pub const AVOGADRO: f64 = 6.0221367E+23;
pub const NA: f64 = AVOGADRO;
pub const MOL: f64 = AVOGADRO;
pub const BOLTZMANN: f64 = 1.380657E-23;
pub const k: f64 = BOLTZMANN;
pub const PLANCK: f64 = 6.6260754E-34;
pub const h: f64 = PLANCK;
pub const GAS_CONSTANT: f64 = NA * k;
pub const R: f64 = GAS_CONSTANT;
pub const FARADAY: f64 = NA * e0;
pub const F: f64 = FARADAY;
pub const BOHR_RADIUS: f64 = 5.29177249E-11;
pub const a0: f64 = BOHR_RADIUS;
pub const VACUUM_PERMITTIVITY: f64 = 8.85419E-12;
pub const VACUUM_PERMEABILITY: f64 = 4.0 * PI * 1E-7;
pub const SPEED_OF_LIGHT: f64 = 2.99792458E+8;
pub const c: f64 = SPEED_OF_LIGHT;
pub const GRAVITATIONAL_CONSTANT: f64 = 6.67259E-11;
pub const FINE_STRUCTURE_CONSTANT: f64 = 7.29735E-3;
pub const DEG_PER_RAD: f64 = 57.2957795130823209;
pub const RAD_PER_DEG: f64 = 0.0174532925199432957;
pub const MM_PER_INCH: f64 = 25.4;
pub const M_PER_FOOT: f64 = 3.048;
pub const JOULE_PER_CAL: f64 = 4.184;
pub const CAL_PER_JOULE: f64 = 1.0 / 4.184;

/// Initial value of the source's writable comparison threshold.
pub const DEFAULT_EPSILON: f64 = 1e-6;
static EPSILON_BITS: AtomicU64 = AtomicU64::new(DEFAULT_EPSILON.to_bits());

/// Read the shared source-compatible comparison threshold.
/// This value is separate from machine precision (`f64::EPSILON`) and from
/// explicitly configured numerical tolerances in individual operations.
pub fn epsilon() -> f64 {
    f64::from_bits(EPSILON_BITS.load(Ordering::Relaxed))
}

/// Replace the shared threshold, preserving every f64 bit pattern.
/// Atomic storage avoids the source writable global's data race; concurrent
/// callers observe a complete old or new value, without locking or validation.
pub fn set_epsilon(value: f64) {
    EPSILON_BITS.store(value.to_bits(), Ordering::Relaxed);
}

/// Exact metadata keys from `OpenMS::Constants::UserParam`.
pub mod user_param {
    pub const IM: &str = "IM";
    pub const FAIMS_CV: &str = "FAIMS_CV";
    pub const ION_MOBILITY: &str = "Ion Mobility";
    pub const INVERSE_REDUCED_ION_MOBILITY: &str = "inverse reduced ion mobility";
    pub const MEAN_INVERSE_REDUCED_ION_MOBILITY_ARRAY: &str =
        "mean inverse reduced ion mobility array";
    pub const ION_MOBILITY_CENTROID: &str = "Ion Mobility Centroid";
    pub const FWHM_IM: &str = "IM Peak FWHM";
    pub const FWHM_IM_AVG: &str = "FWHM_im_avg";
    pub const FWHM_MZ_ppm: &str = "FWHM_ppm";
    pub const FWHM_MZ_AVG: &str = "FWHM_mz_avg";
    pub const SD: &str = "SD";
    pub const SD_ppm: &str = "SD_ppm";
    pub const IonNames: &str = "IonNames";
    pub const CONCAT_PEPTIDE: &str = "concatenated_peptides";
    pub const LOCALIZED_MODIFICATIONS_USERPARAM: &str = "localized_modifications";
    pub const MERGED_CHROMATOGRAM_MZS: &str = "merged_chromatogram_mzs";
    pub const PRECURSOR_ERROR_PPM_USERPARAM: &str = "precursor_mz_error_ppm";
    pub const FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM: &str = "fragment_mz_error_median_ppm";
    pub const FRAGMENT_ERROR_PPM_USERPARAM: &str = "fragment_mass_error_ppm";
    pub const FRAGMENT_ERROR_DA_USERPARAM: &str = "fragment_mass_error_da";
    pub const FRAGMENT_ANNOTATION_USERPARAM: &str = "fragment_annotation";
    pub const PSM_EXPLAINED_ION_CURRENT_USERPARAM: &str = "PSM_explained_ion_current";
    pub const MATCHED_PREFIX_IONS_FRACTION: &str = "matched_prefix_ions_fraction";
    pub const MATCHED_SUFFIX_IONS_FRACTION: &str = "matched_suffix_ions_fraction";
    pub const MATCHED_PREFIX_IONS: &str = "matched_prefix_ions";
    pub const MATCHED_SUFFIX_IONS: &str = "matched_suffix_ions";
    pub const NUM_MATCHED_PEAKS: &str = "num_matched_peaks";
    pub const LONGEST_PEPTIDE_ION_SEQUENCE: &str = "longest_peptide_ion_sequence";
    pub const MATCHED_ION_CURRENT: &str = "matched_ion_current";
    pub const SPECTRUM_REFERENCE: &str = "spectrum_reference";
    pub const ID_MERGE_INDEX: &str = "id_merge_index";
    pub const TARGET_DECOY: &str = "target_decoy";
    pub const DELTA_SCORE: &str = "delta_score";
    pub const ISOTOPE_ERROR: &str = "isotope_error";
    pub const HYPERSCORE_ZSCORE: &str = "hyperscore_zscore";
    pub const LN_NUM_CANDIDATES: &str = "ln_num_candidates";
    pub const MATCHED_ION_CURRENT_FRACTION: &str = "matched_ion_current_fraction";
    pub const COMPLEMENTARY_IONS_FRACTION: &str = "complementary_ions_fraction";
    pub const PEPTIDE_Q_VALUE: &str = "peptide q-value";
    pub const OPENPEPXL_SCORE: &str = "OpenPepXL:score";
    pub const OPENPEPXL_BETA_SEQUENCE: &str = "sequence_beta";
    pub const OPENPEPXL_BETA_ACCESSIONS: &str = "accessions_beta";
    pub const OPENPEPXL_XL_POS1: &str = "xl_pos1";
    pub const OPENPEPXL_XL_POS2: &str = "xl_pos2";
    pub const OPENPEPXL_XL_POS1_PROT: &str = "xl_pos1_protein";
    pub const OPENPEPXL_XL_POS2_PROT: &str = "xl_pos2_protein";
    pub const OPENPEPXL_XL_TYPE: &str = "xl_type";
    pub const OPENPEPXL_XL_RANK: &str = "xl_rank";
    pub const OPENPEPXL_XL_MOD: &str = "xl_mod";
    pub const OPENPEPXL_XL_MASS: &str = "xl_mass";
    pub const OPENPEPXL_XL_TERM_SPEC_ALPHA: &str = "xl_term_spec_alpha";
    pub const OPENPEPXL_XL_TERM_SPEC_BETA: &str = "xl_term_spec_beta";
    pub const OPENPEPXL_HEAVY_SPEC_RT: &str = "spec_heavy_RT";
    pub const OPENPEPXL_HEAVY_SPEC_MZ: &str = "spec_heavy_MZ";
    pub const OPENPEPXL_HEAVY_SPEC_REF: &str = "spectrum_reference_heavy";
    pub const OPENPEPXL_TARGET_DECOY_ALPHA: &str = "xl_target_decoy_alpha";
    pub const OPENPEPXL_TARGET_DECOY_BETA: &str = "xl_target_decoy_beta";
    pub const OPENPEPXL_BETA_PEPEV_PRE: &str = "BetaPepEv:pre";
    pub const OPENPEPXL_BETA_PEPEV_POST: &str = "BetaPepEv:post";
    pub const OPENPEPXL_BETA_PEPEV_START: &str = "BetaPepEv:start";
    pub const OPENPEPXL_BETA_PEPEV_END: &str = "BetaPepEv:end";
    pub const SIRIUS_MZ: &str = "mz";
    pub const SIRIUS_EXACTMASS: &str = "exact_mass";
    pub const SIRIUS_EXPLANATION: &str = "explanation";
    pub const SIRIUS_SCORE: &str = "score";
    pub const SIRIUS_PEAKMZ: &str = "peak_mz";
    pub const SIRIUS_ANNOTATED_SUMFORMULA: &str = "annotated_sumformula";
    pub const SIRIUS_ANNOTATED_ADDUCT: &str = "annotated_adduct";
    pub const SIRIUS_DECOY: &str = "decoy";
    pub const SIRIUS_FEATURE_ID: &str = "feat_id";
    pub const XFDR_FDR: &str = "XFDR:FDR";
    pub const IIMN_BEST_ION: &str = "best ion";
    pub const IIMN_ADDUCT_PARTNERS: &str = "partners";
    pub const IIMN_ROW_ID: &str = "row ID";
    pub const IIMN_ANNOTATION_NETWORK_NUMBER: &str = "annotation network number";
    pub const ADDUCT_GROUP: &str = "Group";
    pub const IIMN_LINKED_GROUPS: &str = "LinkedGroups";
    pub const DC_CHARGE_ADDUCTS: &str = "dc_charge_adducts";
    pub const NUM_OF_MASSTRACES: &str = "num_of_masstraces";
    pub const NUM_OF_DATAPOINTS: &str = "num_of_datapoints";
    pub const MSM_METABOLITE_NAME: &str = "Metabolite_Name";
    pub const MSM_INCHI_STRING: &str = "Inchi_String";
    pub const MSM_SMILES_STRING: &str = "SMILES_String";
    pub const MSM_PRECURSOR_ADDUCT: &str = "Precursor_Ion";
    pub const MSM_SUM_FORMULA: &str = "Sum_Formula";
    pub const MSM_CCS: &str = "CCS";
    pub const BASE_NAME: &str = "base_name";
    pub const SIGNIFICANCE_THRESHOLD: &str = "significance_threshold";
    pub const RANK: &str = "rank";
    pub const NUM_PEAKS: &str = "num_peaks";
    pub const MODIFICATION_DEFINITIONS: &str = "modification_definitions";
}
