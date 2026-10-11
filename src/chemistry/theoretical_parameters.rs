// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `DefaultParamHandler` contract of `TheoreticalSpectrumGenerator`: the
//! thirty parameters its constructor declares
//! (`TheoreticalSpectrumGenerator.cpp:25-112`) and the members its
//! `updateMembers_` derives from them (`:1263-1298`).

use super::{
    TheoreticalIonSeries, TheoreticalIsotopeModel, TheoreticalIsotopeParameters,
    TheoreticalSpectrumGenerator,
};
use crate::param::handler::{int_parameter, set_flag};
use crate::param::{DefaultParamHandling, Param, ParamValue};
use crate::{Error, Result};

/// The ion-series flags in declaration order, with their descriptions.
const SERIES: [(&str, TheoreticalIonSeries, &str); 8] = [
    (
        "add_y_ions",
        TheoreticalIonSeries::Y,
        "Add peaks of y-ions to the spectrum",
    ),
    (
        "add_b_ions",
        TheoreticalIonSeries::B,
        "Add peaks of b-ions to the spectrum",
    ),
    (
        "add_a_ions",
        TheoreticalIonSeries::A,
        "Add peaks of a-ions to the spectrum",
    ),
    (
        "add_c_ions",
        TheoreticalIonSeries::C,
        "Add peaks of c-ions to the spectrum",
    ),
    (
        "add_x_ions",
        TheoreticalIonSeries::X,
        "Add peaks of  x-ions to the spectrum",
    ),
    (
        "add_z_ions",
        TheoreticalIonSeries::Z,
        "Add peaks of z-ions to the spectrum (sometimes observed in CID and for some AAs in ExD due to H abstraction)",
    ),
    (
        "add_zp1_ions",
        TheoreticalIonSeries::ZPlusOne,
        "Add peaks of z+1-radical cations (also [z+H]*^{+} or simply z*) to the spectrum (often observed in ExD)",
    ),
    (
        "add_zp2_ions",
        TheoreticalIonSeries::ZPlusTwo,
        "Add peaks of z+2-radical cations (also [z+2H]*^{2+} or simply z') to the spectrum (often observed in ExD esp. with higher precursor charges >3 and smaller z-ions.)",
    ),
];

/// The order in which `update_members` lists
/// the selected series: the source's generation order.
const GENERATION_ORDER: [TheoreticalIonSeries; 8] = [
    TheoreticalIonSeries::B,
    TheoreticalIonSeries::Y,
    TheoreticalIonSeries::A,
    TheoreticalIonSeries::C,
    TheoreticalIonSeries::X,
    TheoreticalIonSeries::Z,
    TheoreticalIonSeries::ZPlusOne,
    TheoreticalIonSeries::ZPlusTwo,
];

/// The intensity parameters after the series flags, in declaration order:
/// name, description, maximum (the minimum is always 0).
const INTENSITIES: [(&str, &str, Option<f64>); 10] = [
    ("y_intensity", "Intensity of the y-ions", None),
    ("b_intensity", "Intensity of the b-ions", None),
    ("a_intensity", "Intensity of the a-ions", None),
    ("c_intensity", "Intensity of the c-ions", None),
    ("x_intensity", "Intensity of the x-ions", None),
    ("z_intensity", "Intensity of the z-ions", None),
    (
        "relative_loss_intensity",
        "Intensity of loss ions, in relation to the intact ion intensity",
        Some(1.0),
    ),
    (
        "precursor_intensity",
        "Intensity of the precursor peak",
        None,
    ),
    (
        "precursor_H2O_intensity",
        "Intensity of the H2O loss peak of the precursor",
        None,
    ),
    (
        "precursor_NH3_intensity",
        "Intensity of the NH3 loss peak of the precursor",
        None,
    ),
];

impl TheoreticalSpectrumGenerator {
    /// The value of each intensity parameter of [`INTENSITIES`], in order.
    fn intensity_values(&self) -> [f64; 10] {
        let i = self.intensities;
        [
            i.y,
            i.b,
            i.a,
            i.c,
            i.x,
            i.z,
            self.relative_loss_intensity,
            i.precursor,
            i.precursor_water_loss,
            i.precursor_ammonia_loss,
        ]
    }

    /// The boolean parameters other than the series flags, in declaration
    /// order, with their values.
    fn flags(&self) -> [(&'static str, bool, &'static str); 8] {
        [
            (
                "add_metainfo",
                self.add_metainfo,
                "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++",
            ),
            (
                "add_losses",
                self.add_losses,
                "Adds common losses to those ion expect to have them, only water and ammonia loss is considered",
            ),
            (
                "add_term_losses",
                self.add_terminal_losses,
                "Adds common N- and C-term losses (only if add_losses=true and isotope_model=none), only water and ammonia loss is considered.",
            ),
            (
                "add_internal_fragments",
                self.add_internal_fragments,
                "Add b and a type fragments from internal cleavage events.",
            ),
            (
                "sort_by_position",
                self.sort_by_position,
                "Sort output by position",
            ),
            (
                "add_precursor_peaks",
                self.add_precursor_peaks,
                "Adds peaks of the unfragmented precursor ion to the spectrum",
            ),
            (
                "add_all_precursor_charges",
                self.add_all_precursor_charges,
                "Adds precursor peaks with all charges in the given range",
            ),
            (
                "add_abundant_immonium_ions",
                self.add_abundant_immonium_ions,
                "Add most abundant immonium ions (for Proline, Cystein, Iso/Leucine, Histidin, Phenylalanin, Tyrosine, Tryptophan)",
            ),
        ]
    }

    /// The parameter tree of this generator: the constructor's declarations
    /// with this generator's values.
    fn parameter_tree(&self) -> Result<Param> {
        let mut param = Param::new();
        let (model, max_isotope, probability) = match self.isotope_model {
            TheoreticalIsotopeModel::None => (
                "none",
                ParamValue::Integer(i64::from(self.isotope_parameters.max_isotope)),
                self.isotope_parameters.max_isotope_probability,
            ),
            TheoreticalIsotopeModel::Coarse { max_peaks } => (
                "coarse",
                int_parameter(max_peaks, "max_isotope")?,
                self.isotope_parameters.max_isotope_probability,
            ),
            TheoreticalIsotopeModel::Fine {
                unexplained_probability,
            } => (
                "fine",
                ParamValue::Integer(i64::from(self.isotope_parameters.max_isotope)),
                unexplained_probability,
            ),
        };
        param.set_value(
            "isotope_model",
            ParamValue::String(model.into()),
            "Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add accurate isotopic peaks. Note that adding isotopic peaks is very slow.",
            &[],
        )?;
        param.set_valid_strings(
            "isotope_model",
            &["none".into(), "coarse".into(), "fine".into()],
        )?;
        param.set_value(
            "max_isotope",
            max_isotope,
            "Defines the maximal isotopic peak which is added if 'isotope_model' is 'coarse'",
            &[],
        )?;
        param.set_value(
            "max_isotope_probability",
            ParamValue::Float(probability),
            "Defines the maximal isotopic probability to cover if 'isotope_model' is 'fine'",
            &[],
        )?;
        let flags = self.flags();
        // add_metainfo .. sort_by_position, then the three precursor/immonium
        // flags, then add_first_prefix_ion, then the series.
        for (key, value, description) in flags {
            set_flag(&mut param, key, value, description, &[])?;
        }
        set_flag(
            &mut param,
            "add_first_prefix_ion",
            self.add_first_prefix_ion,
            "If set to true e.g. b1 ions are added",
            &[],
        )?;
        for (key, series, description) in SERIES {
            set_flag(
                &mut param,
                key,
                self.ion_series.contains(&series),
                description,
                &[],
            )?;
        }
        for ((key, description, max), value) in INTENSITIES.into_iter().zip(self.intensity_values())
        {
            param.set_value(key, ParamValue::Float(value), description, &[])?;
            param.set_min_float(key, 0.0)?;
            if let Some(max) = max {
                param.set_max_float(key, max)?;
            }
        }
        Ok(param)
    }
}

/// Source `TheoreticalSpectrumGenerator`'s `DefaultParamHandler` base.
///
/// | Parameter | Default | Field |
/// |---|---|---|
/// | `isotope_model` (`none`, `coarse`, `fine`) | `none` | [`isotope_model`](TheoreticalSpectrumGenerator::isotope_model) |
/// | `max_isotope` | 2 | `Coarse { max_peaks }`, and [`isotope_parameters`](TheoreticalSpectrumGenerator::isotope_parameters) |
/// | `max_isotope_probability` | 0.05 | `Fine { unexplained_probability }`, and [`isotope_parameters`](TheoreticalSpectrumGenerator::isotope_parameters) |
/// | `add_metainfo` | `false` | [`add_metainfo`](TheoreticalSpectrumGenerator::add_metainfo) |
/// | `add_losses` | `false` | [`add_losses`](TheoreticalSpectrumGenerator::add_losses) |
/// | `add_term_losses` | `false` | [`add_terminal_losses`](TheoreticalSpectrumGenerator::add_terminal_losses) |
/// | `add_internal_fragments` | `false` | [`add_internal_fragments`](TheoreticalSpectrumGenerator::add_internal_fragments) |
/// | `sort_by_position` | `true` | [`sort_by_position`](TheoreticalSpectrumGenerator::sort_by_position) |
/// | `add_precursor_peaks` | `false` | [`add_precursor_peaks`](TheoreticalSpectrumGenerator::add_precursor_peaks) |
/// | `add_all_precursor_charges` | `false` | [`add_all_precursor_charges`](TheoreticalSpectrumGenerator::add_all_precursor_charges) |
/// | `add_abundant_immonium_ions` | `false` | [`add_abundant_immonium_ions`](TheoreticalSpectrumGenerator::add_abundant_immonium_ions) |
/// | `add_first_prefix_ion` | `false` | [`add_first_prefix_ion`](TheoreticalSpectrumGenerator::add_first_prefix_ion) |
/// | `add_y_ions`, `add_b_ions` | `true` | [`ion_series`](TheoreticalSpectrumGenerator::ion_series) |
/// | `add_a_ions`, `add_c_ions`, `add_x_ions`, `add_z_ions`, `add_zp1_ions`, `add_zp2_ions` | `false` | [`ion_series`](TheoreticalSpectrumGenerator::ion_series) |
/// | `y_intensity` ... `z_intensity` (min 0) | 1.0 | [`intensities`](TheoreticalSpectrumGenerator::intensities) |
/// | `relative_loss_intensity` (0..1) | 0.1 | [`relative_loss_intensity`](TheoreticalSpectrumGenerator::relative_loss_intensity) |
/// | `precursor_intensity`, `precursor_H2O_intensity`, `precursor_NH3_intensity` (min 0) | 1.0 | [`intensities`](TheoreticalSpectrumGenerator::intensities) |
impl DefaultParamHandling for TheoreticalSpectrumGenerator {
    const NAME: &'static str = "TheoreticalSpectrumGenerator";

    fn defaults() -> Result<Param> {
        Self::default().parameter_tree()
    }

    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a coarse model's peak count exceeds the
    /// source's 32-bit `int`.
    fn to_param(&self) -> Result<Param> {
        self.parameter_tree()
    }

    /// Source `updateMembers_`.
    ///
    /// The selected series are listed in the source's generation order. Values
    /// the generator refuses (a zero coarse peak count, terminal losses without
    /// `add_losses` or with isotopes) are accepted here, as by the source, and
    /// refused by generation.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a negative `max_isotope` with
    /// `isotope_model` = `coarse`: the source converts it to a huge `Size`
    /// peak count, which this port cannot hold.
    fn update_members(param: &Param) -> Result<Self> {
        let flag = |key: &str| -> Result<bool> { param.value(key)?.to_bool() };
        let float = |key: &str| -> Result<f64> { param.value(key)?.to_f64() };
        let isotope_parameters = TheoreticalIsotopeParameters {
            max_isotope: param.value("max_isotope")?.to_i32()?,
            max_isotope_probability: float("max_isotope_probability")?,
        };
        let isotope_model = match param.value("isotope_model")?.as_str()? {
            "coarse" => TheoreticalIsotopeModel::Coarse {
                max_peaks: usize::try_from(isotope_parameters.max_isotope).map_err(|_| {
                    Error::InvalidValue(
                        "TheoreticalSpectrumGenerator: a negative 'max_isotope' has no coarse isotope peak count"
                            .into(),
                    )
                })?,
            },
            "fine" => TheoreticalIsotopeModel::Fine {
                unexplained_probability: isotope_parameters.max_isotope_probability,
            },
            // 'none' (the restriction admits nothing else): the source resets
            // isotope_model_ to 0 so that no earlier setting lingers.
            _ => TheoreticalIsotopeModel::None,
        };
        let mut ion_series = Vec::new();
        for series in GENERATION_ORDER {
            let key = SERIES
                .iter()
                .find(|(_, candidate, _)| *candidate == series)
                .map(|(key, _, _)| *key)
                .ok_or_else(|| Error::InvalidValue("unknown ion series".into()))?;
            if flag(key)? {
                ion_series.push(series);
            }
        }
        let mut generator = Self {
            ion_series,
            isotope_model,
            isotope_parameters,
            add_first_prefix_ion: flag("add_first_prefix_ion")?,
            add_losses: flag("add_losses")?,
            add_terminal_losses: flag("add_term_losses")?,
            relative_loss_intensity: float("relative_loss_intensity")?,
            add_precursor_peaks: flag("add_precursor_peaks")?,
            add_all_precursor_charges: flag("add_all_precursor_charges")?,
            add_internal_fragments: flag("add_internal_fragments")?,
            add_abundant_immonium_ions: flag("add_abundant_immonium_ions")?,
            add_metainfo: flag("add_metainfo")?,
            sort_by_position: flag("sort_by_position")?,
            ..Self::default()
        };
        let i = &mut generator.intensities;
        i.y = float("y_intensity")?;
        i.b = float("b_intensity")?;
        i.a = float("a_intensity")?;
        i.c = float("c_intensity")?;
        i.x = float("x_intensity")?;
        i.z = float("z_intensity")?;
        i.precursor = float("precursor_intensity")?;
        i.precursor_water_loss = float("precursor_H2O_intensity")?;
        i.precursor_ammonia_loss = float("precursor_NH3_intensity")?;
        Ok(generator)
    }
}
