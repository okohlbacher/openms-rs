// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::constants::{self, user_param};
use std::collections::BTreeMap;

// Expected bits and strings were emitted by the actual, unchanged C++ header.
fn reference(kind: &str) -> BTreeMap<&'static str, &'static str> {
    include_str!("data/constants_reference.tsv")
        .lines()
        .filter_map(|line| {
            let mut parts = line.splitn(3, '\t');
            let group = parts.next()?;
            let name = parts.next()?;
            let value = parts.next()?;
            (group == kind).then_some((name, value))
        })
        .collect()
}

#[test]
fn all_numeric_source_bits_aliases_and_mutable_epsilon() {
    let actual = [
        ("PI", constants::PI),
        ("E", constants::E),
        ("EPSILON", constants::epsilon()),
        ("ELEMENTARY_CHARGE", constants::ELEMENTARY_CHARGE),
        ("e0", constants::e0),
        ("ELECTRON_MASS", constants::ELECTRON_MASS),
        ("ELECTRON_MASS_U", constants::ELECTRON_MASS_U),
        ("PROTON_MASS", constants::PROTON_MASS),
        ("PROTON_MASS_U", constants::PROTON_MASS_U),
        ("C13C12_MASSDIFF_U", constants::C13C12_MASSDIFF_U),
        ("ISOTOPE_MASSDIFF_55K_U", constants::ISOTOPE_MASSDIFF_55K_U),
        ("NEUTRON_MASS", constants::NEUTRON_MASS),
        ("NEUTRON_MASS_U", constants::NEUTRON_MASS_U),
        ("AVOGADRO", constants::AVOGADRO),
        ("NA", constants::NA),
        ("MOL", constants::MOL),
        ("BOLTZMANN", constants::BOLTZMANN),
        ("k", constants::k),
        ("PLANCK", constants::PLANCK),
        ("h", constants::h),
        ("GAS_CONSTANT", constants::GAS_CONSTANT),
        ("R", constants::R),
        ("FARADAY", constants::FARADAY),
        ("F", constants::F),
        ("BOHR_RADIUS", constants::BOHR_RADIUS),
        ("a0", constants::a0),
        ("VACUUM_PERMITTIVITY", constants::VACUUM_PERMITTIVITY),
        ("VACUUM_PERMEABILITY", constants::VACUUM_PERMEABILITY),
        ("SPEED_OF_LIGHT", constants::SPEED_OF_LIGHT),
        ("c", constants::c),
        ("GRAVITATIONAL_CONSTANT", constants::GRAVITATIONAL_CONSTANT),
        (
            "FINE_STRUCTURE_CONSTANT",
            constants::FINE_STRUCTURE_CONSTANT,
        ),
        ("DEG_PER_RAD", constants::DEG_PER_RAD),
        ("RAD_PER_DEG", constants::RAD_PER_DEG),
        ("MM_PER_INCH", constants::MM_PER_INCH),
        ("M_PER_FOOT", constants::M_PER_FOOT),
        ("JOULE_PER_CAL", constants::JOULE_PER_CAL),
        ("CAL_PER_JOULE", constants::CAL_PER_JOULE),
    ];
    let expected = reference("number");
    assert_eq!(actual.len(), expected.len());
    for (name, value) in actual {
        assert_eq!(
            value.to_bits(),
            u64::from_str_radix(expected[name], 16).unwrap(),
            "{name}"
        );
    }
    assert_eq!(
        openms::chemistry::PROTON_MASS_U.to_bits(),
        constants::PROTON_MASS_U.to_bits()
    );
    assert_eq!(
        openms::chemistry::ELECTRON_MASS_U.to_bits(),
        constants::ELECTRON_MASS_U.to_bits()
    );
    assert_eq!(
        openms::chemistry::C13C12_MASSDIFF_U.to_bits(),
        constants::C13C12_MASSDIFF_U.to_bits()
    );
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            constants::set_epsilon(constants::DEFAULT_EPSILON);
        }
    }
    let _reset = Reset;
    for bits in [
        0,
        (-0.0f64).to_bits(),
        1e-9f64.to_bits(),
        (-2.0f64).to_bits(),
        f64::INFINITY.to_bits(),
        0x7ff8_1234_5678_9abc,
    ] {
        constants::set_epsilon(f64::from_bits(bits));
        assert_eq!(constants::epsilon().to_bits(), bits);
    }
}

#[test]
fn all_user_parameter_spellings_match_the_source() {
    let actual = [
        ("IM", user_param::IM),
        ("FAIMS_CV", user_param::FAIMS_CV),
        ("ION_MOBILITY", user_param::ION_MOBILITY),
        (
            "INVERSE_REDUCED_ION_MOBILITY",
            user_param::INVERSE_REDUCED_ION_MOBILITY,
        ),
        (
            "MEAN_INVERSE_REDUCED_ION_MOBILITY_ARRAY",
            user_param::MEAN_INVERSE_REDUCED_ION_MOBILITY_ARRAY,
        ),
        ("ION_MOBILITY_CENTROID", user_param::ION_MOBILITY_CENTROID),
        ("FWHM_IM", user_param::FWHM_IM),
        ("FWHM_IM_AVG", user_param::FWHM_IM_AVG),
        ("FWHM_MZ_ppm", user_param::FWHM_MZ_ppm),
        ("FWHM_MZ_AVG", user_param::FWHM_MZ_AVG),
        ("SD", user_param::SD),
        ("SD_ppm", user_param::SD_ppm),
        ("IonNames", user_param::IonNames),
        ("CONCAT_PEPTIDE", user_param::CONCAT_PEPTIDE),
        (
            "LOCALIZED_MODIFICATIONS_USERPARAM",
            user_param::LOCALIZED_MODIFICATIONS_USERPARAM,
        ),
        (
            "MERGED_CHROMATOGRAM_MZS",
            user_param::MERGED_CHROMATOGRAM_MZS,
        ),
        (
            "PRECURSOR_ERROR_PPM_USERPARAM",
            user_param::PRECURSOR_ERROR_PPM_USERPARAM,
        ),
        (
            "FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM",
            user_param::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM,
        ),
        (
            "FRAGMENT_ERROR_PPM_USERPARAM",
            user_param::FRAGMENT_ERROR_PPM_USERPARAM,
        ),
        (
            "FRAGMENT_ERROR_DA_USERPARAM",
            user_param::FRAGMENT_ERROR_DA_USERPARAM,
        ),
        (
            "FRAGMENT_ANNOTATION_USERPARAM",
            user_param::FRAGMENT_ANNOTATION_USERPARAM,
        ),
        (
            "PSM_EXPLAINED_ION_CURRENT_USERPARAM",
            user_param::PSM_EXPLAINED_ION_CURRENT_USERPARAM,
        ),
        (
            "MATCHED_PREFIX_IONS_FRACTION",
            user_param::MATCHED_PREFIX_IONS_FRACTION,
        ),
        (
            "MATCHED_SUFFIX_IONS_FRACTION",
            user_param::MATCHED_SUFFIX_IONS_FRACTION,
        ),
        ("MATCHED_PREFIX_IONS", user_param::MATCHED_PREFIX_IONS),
        ("MATCHED_SUFFIX_IONS", user_param::MATCHED_SUFFIX_IONS),
        ("NUM_MATCHED_PEAKS", user_param::NUM_MATCHED_PEAKS),
        (
            "LONGEST_PEPTIDE_ION_SEQUENCE",
            user_param::LONGEST_PEPTIDE_ION_SEQUENCE,
        ),
        ("MATCHED_ION_CURRENT", user_param::MATCHED_ION_CURRENT),
        ("SPECTRUM_REFERENCE", user_param::SPECTRUM_REFERENCE),
        ("ID_MERGE_INDEX", user_param::ID_MERGE_INDEX),
        ("TARGET_DECOY", user_param::TARGET_DECOY),
        ("DELTA_SCORE", user_param::DELTA_SCORE),
        ("ISOTOPE_ERROR", user_param::ISOTOPE_ERROR),
        ("HYPERSCORE_ZSCORE", user_param::HYPERSCORE_ZSCORE),
        ("LN_NUM_CANDIDATES", user_param::LN_NUM_CANDIDATES),
        (
            "MATCHED_ION_CURRENT_FRACTION",
            user_param::MATCHED_ION_CURRENT_FRACTION,
        ),
        (
            "COMPLEMENTARY_IONS_FRACTION",
            user_param::COMPLEMENTARY_IONS_FRACTION,
        ),
        ("PEPTIDE_Q_VALUE", user_param::PEPTIDE_Q_VALUE),
        ("OPENPEPXL_SCORE", user_param::OPENPEPXL_SCORE),
        (
            "OPENPEPXL_BETA_SEQUENCE",
            user_param::OPENPEPXL_BETA_SEQUENCE,
        ),
        (
            "OPENPEPXL_BETA_ACCESSIONS",
            user_param::OPENPEPXL_BETA_ACCESSIONS,
        ),
        ("OPENPEPXL_XL_POS1", user_param::OPENPEPXL_XL_POS1),
        ("OPENPEPXL_XL_POS2", user_param::OPENPEPXL_XL_POS2),
        ("OPENPEPXL_XL_POS1_PROT", user_param::OPENPEPXL_XL_POS1_PROT),
        ("OPENPEPXL_XL_POS2_PROT", user_param::OPENPEPXL_XL_POS2_PROT),
        ("OPENPEPXL_XL_TYPE", user_param::OPENPEPXL_XL_TYPE),
        ("OPENPEPXL_XL_RANK", user_param::OPENPEPXL_XL_RANK),
        ("OPENPEPXL_XL_MOD", user_param::OPENPEPXL_XL_MOD),
        ("OPENPEPXL_XL_MASS", user_param::OPENPEPXL_XL_MASS),
        (
            "OPENPEPXL_XL_TERM_SPEC_ALPHA",
            user_param::OPENPEPXL_XL_TERM_SPEC_ALPHA,
        ),
        (
            "OPENPEPXL_XL_TERM_SPEC_BETA",
            user_param::OPENPEPXL_XL_TERM_SPEC_BETA,
        ),
        (
            "OPENPEPXL_HEAVY_SPEC_RT",
            user_param::OPENPEPXL_HEAVY_SPEC_RT,
        ),
        (
            "OPENPEPXL_HEAVY_SPEC_MZ",
            user_param::OPENPEPXL_HEAVY_SPEC_MZ,
        ),
        (
            "OPENPEPXL_HEAVY_SPEC_REF",
            user_param::OPENPEPXL_HEAVY_SPEC_REF,
        ),
        (
            "OPENPEPXL_TARGET_DECOY_ALPHA",
            user_param::OPENPEPXL_TARGET_DECOY_ALPHA,
        ),
        (
            "OPENPEPXL_TARGET_DECOY_BETA",
            user_param::OPENPEPXL_TARGET_DECOY_BETA,
        ),
        (
            "OPENPEPXL_BETA_PEPEV_PRE",
            user_param::OPENPEPXL_BETA_PEPEV_PRE,
        ),
        (
            "OPENPEPXL_BETA_PEPEV_POST",
            user_param::OPENPEPXL_BETA_PEPEV_POST,
        ),
        (
            "OPENPEPXL_BETA_PEPEV_START",
            user_param::OPENPEPXL_BETA_PEPEV_START,
        ),
        (
            "OPENPEPXL_BETA_PEPEV_END",
            user_param::OPENPEPXL_BETA_PEPEV_END,
        ),
        ("SIRIUS_MZ", user_param::SIRIUS_MZ),
        ("SIRIUS_EXACTMASS", user_param::SIRIUS_EXACTMASS),
        ("SIRIUS_EXPLANATION", user_param::SIRIUS_EXPLANATION),
        ("SIRIUS_SCORE", user_param::SIRIUS_SCORE),
        ("SIRIUS_PEAKMZ", user_param::SIRIUS_PEAKMZ),
        (
            "SIRIUS_ANNOTATED_SUMFORMULA",
            user_param::SIRIUS_ANNOTATED_SUMFORMULA,
        ),
        (
            "SIRIUS_ANNOTATED_ADDUCT",
            user_param::SIRIUS_ANNOTATED_ADDUCT,
        ),
        ("SIRIUS_DECOY", user_param::SIRIUS_DECOY),
        ("SIRIUS_FEATURE_ID", user_param::SIRIUS_FEATURE_ID),
        ("XFDR_FDR", user_param::XFDR_FDR),
        ("IIMN_BEST_ION", user_param::IIMN_BEST_ION),
        ("IIMN_ADDUCT_PARTNERS", user_param::IIMN_ADDUCT_PARTNERS),
        ("IIMN_ROW_ID", user_param::IIMN_ROW_ID),
        (
            "IIMN_ANNOTATION_NETWORK_NUMBER",
            user_param::IIMN_ANNOTATION_NETWORK_NUMBER,
        ),
        ("ADDUCT_GROUP", user_param::ADDUCT_GROUP),
        ("IIMN_LINKED_GROUPS", user_param::IIMN_LINKED_GROUPS),
        ("DC_CHARGE_ADDUCTS", user_param::DC_CHARGE_ADDUCTS),
        ("NUM_OF_MASSTRACES", user_param::NUM_OF_MASSTRACES),
        ("NUM_OF_DATAPOINTS", user_param::NUM_OF_DATAPOINTS),
        ("MSM_METABOLITE_NAME", user_param::MSM_METABOLITE_NAME),
        ("MSM_INCHI_STRING", user_param::MSM_INCHI_STRING),
        ("MSM_SMILES_STRING", user_param::MSM_SMILES_STRING),
        ("MSM_PRECURSOR_ADDUCT", user_param::MSM_PRECURSOR_ADDUCT),
        ("MSM_SUM_FORMULA", user_param::MSM_SUM_FORMULA),
        ("MSM_CCS", user_param::MSM_CCS),
        ("BASE_NAME", user_param::BASE_NAME),
        ("SIGNIFICANCE_THRESHOLD", user_param::SIGNIFICANCE_THRESHOLD),
        ("RANK", user_param::RANK),
        ("NUM_PEAKS", user_param::NUM_PEAKS),
        (
            "MODIFICATION_DEFINITIONS",
            user_param::MODIFICATION_DEFINITIONS,
        ),
    ];
    let expected = reference("text");
    assert_eq!(actual.len(), expected.len());
    for (name, value) in actual {
        assert_eq!(value, expected[name], "{name}");
    }
}
