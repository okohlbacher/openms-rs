// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The `DefaultParamHandler` contract of `NucleicAcidSpectrumGenerator`: the
//! twenty-three parameters its constructor declares
//! (`NucleicAcidSpectrumGenerator.cpp:18-74`) and the members its
//! `updateMembers_` reads from them (`:485-510`).

use super::NucleicAcidSpectrumGenerator;
use crate::Result;
use crate::param::handler::set_flag;
use crate::param::{DefaultParamHandling, Param, ParamValue};

/// A boolean member of the generator.
type FlagField = fn(&mut NucleicAcidSpectrumGenerator) -> &mut bool;
/// An intensity member of the generator.
type IntensityField = fn(&mut NucleicAcidSpectrumGenerator) -> &mut f64;

/// The thirteen flags in declaration order: key, description, member.
const FLAGS: [(&str, &str, FlagField); 13] = [
    (
        "add_metainfo",
        "Adds the type of peaks as meta information to the peaks, e.g. c1, y2, a3-B",
        |g| &mut g.add_metainfo,
    ),
    (
        "add_precursor_peaks",
        "Adds peaks of the unfragmented precursor ion to the spectrum",
        |g| &mut g.add_precursor_peaks,
    ),
    (
        "add_all_precursor_charges",
        "Adds precursor peaks with all charges in the given range",
        |g| &mut g.add_all_precursor_charges,
    ),
    (
        "add_first_prefix_ion",
        "If set to true a1, b1, ..., z1 ions are added",
        |g| &mut g.add_first_prefix_ion,
    ),
    (
        "add_a_ions",
        "Add peaks of a-ions to the spectrum",
        |g| &mut g.add_a_ions,
    ),
    (
        "add_b_ions",
        "Add peaks of b-ions to the spectrum",
        |g| &mut g.add_b_ions,
    ),
    (
        "add_c_ions",
        "Add peaks of c-ions to the spectrum",
        |g| &mut g.add_c_ions,
    ),
    (
        "add_d_ions",
        "Add peaks of d-ions to the spectrum",
        |g| &mut g.add_d_ions,
    ),
    (
        "add_w_ions",
        "Add peaks of w-ions to the spectrum",
        |g| &mut g.add_w_ions,
    ),
    (
        "add_x_ions",
        "Add peaks of  x-ions to the spectrum",
        |g| &mut g.add_x_ions,
    ),
    (
        "add_y_ions",
        "Add peaks of y-ions to the spectrum",
        |g| &mut g.add_y_ions,
    ),
    (
        "add_z_ions",
        "Add peaks of z-ions to the spectrum",
        |g| &mut g.add_z_ions,
    ),
    (
        "add_a-B_ions",
        "Add peaks of a-B-ions to the spectrum",
        |g| &mut g.add_a_minus_b_ions,
    ),
];

/// The ten intensities in declaration order, none restricted: key,
/// description, member.
const INTENSITIES: [(&str, &str, IntensityField); 10] = [
    ("a_intensity", "Intensity of the a-ions", |g| {
        &mut g.a_intensity
    }),
    ("b_intensity", "Intensity of the b-ions", |g| {
        &mut g.b_intensity
    }),
    ("c_intensity", "Intensity of the c-ions", |g| {
        &mut g.c_intensity
    }),
    ("d_intensity", "Intensity of the d-ions", |g| {
        &mut g.d_intensity
    }),
    ("w_intensity", "Intensity of the w-ions", |g| {
        &mut g.w_intensity
    }),
    ("x_intensity", "Intensity of the x-ions", |g| {
        &mut g.x_intensity
    }),
    ("y_intensity", "Intensity of the y-ions", |g| {
        &mut g.y_intensity
    }),
    ("z_intensity", "Intensity of the z-ions", |g| {
        &mut g.z_intensity
    }),
    ("a-B_intensity", "Intensity of the a-B-ions", |g| {
        &mut g.a_minus_b_intensity
    }),
    (
        "precursor_intensity",
        "Intensity of the precursor peak",
        |g| &mut g.precursor_intensity,
    ),
];

impl NucleicAcidSpectrumGenerator {
    /// The constructor's declarations carrying this generator's values.
    fn parameter_tree(&self) -> Result<Param> {
        // The accessors take `&mut`; read through a scratch copy.
        let mut values = self.clone();
        let mut param = Param::new();
        for (key, description, field) in FLAGS {
            set_flag(&mut param, key, *field(&mut values), description, &[])?;
        }
        for (key, description, field) in INTENSITIES {
            param.set_value(
                key,
                ParamValue::Float(*field(&mut values)),
                description,
                &[],
            )?;
        }
        Ok(param)
    }
}

/// Source `NucleicAcidSpectrumGenerator`'s `DefaultParamHandler` base.
///
/// Every field is a parameter; the key is the field name, except
/// `add_a-B_ions` ([`add_a_minus_b_ions`](NucleicAcidSpectrumGenerator::add_a_minus_b_ions))
/// and `a-B_intensity`
/// ([`a_minus_b_intensity`](NucleicAcidSpectrumGenerator::a_minus_b_intensity)).
/// Only `add_b_ions` and `add_y_ions` default to `true`; every intensity
/// defaults to 1.0 and has no restriction.
impl DefaultParamHandling for NucleicAcidSpectrumGenerator {
    const NAME: &'static str = "NucleicAcidSpectrumGenerator";

    fn defaults() -> Result<Param> {
        Self::default().parameter_tree()
    }

    fn to_param(&self) -> Result<Param> {
        self.parameter_tree()
    }

    /// Source `updateMembers_`: every member read back from its parameter.
    fn update_members(param: &Param) -> Result<Self> {
        let mut generator = Self::default();
        for (key, _, field) in FLAGS {
            *field(&mut generator) = param.value(key)?.to_bool()?;
        }
        for (key, _, field) in INTENSITIES {
            *field(&mut generator) = param.value(key)?.to_f64()?;
        }
        Ok(generator)
    }
}
