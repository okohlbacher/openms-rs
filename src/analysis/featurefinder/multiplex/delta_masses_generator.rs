// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The complete list of mass shifts that isotopic labelling can cause.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexDeltaMassesGenerator.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use std::collections::BTreeMap;

use super::delta_masses::{DeltaMass, LabelSet, MultiplexDeltaMasses};
use crate::chemistry::AASequence;
use crate::format::file_info::text_format::ostream_g;
use crate::math::source_sort::source_sort_by;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};

/// Name of the source's parameter handler (`DefaultParamHandler("labels")`).
pub const LABELS_PARAM_HANDLER_NAME: &str = "labels";

/// Ceiling on the number of mass shift patterns one generator may build,
/// checked from `max_nr_labelled_aas` before any pattern is generated. The
/// source has no bound and loops `(max + 2)^2` times for SILAC.
pub const MAX_DELTA_MASS_PATTERNS: usize = 100_000;

/// The labels the source accepts outside numeric labelling, as one
/// space-separated string (source `MultiplexDeltaMassesGenerator.cpp:309`).
/// A label is accepted when it is a **substring** of this string; see
/// [`MultiplexDeltaMassesGenerator::new`].
pub const KNOWN_LABELS: &str = "Arg6 Arg10 Lys4 Lys6 Lys8 Leu3 Dimethyl0 Dimethyl4 Dimethyl6 Dimethyl8 ICPL0 ICPL4 ICPL6 ICPL10 no_label";

/// Complete label information (source `MultiplexDeltaMassesGenerator::Label`).
#[derive(Clone, Debug, PartialEq)]
pub struct Label {
    /// Short name as used in the `labels` parameter, e.g. `Arg10` (source
    /// `short_name`).
    pub short_name: String,
    /// UniMod name as it appears in peptide sequences, e.g.
    /// `Label:13C(6)15N(4)` (source `long_name`).
    pub long_name: String,
    /// Description, also the parameter description (source `description`).
    pub description: String,
    /// Mass shift in Da (source `delta_mass`).
    pub delta_mass: f64,
}

impl Label {
    /// A label (source constructor).
    pub fn new(
        short_name: impl Into<String>,
        long_name: impl Into<String>,
        description: impl Into<String>,
        delta_mass: f64,
    ) -> Self {
        Self {
            short_name: short_name.into(),
            long_name: long_name.into(),
            description: description.into(),
            delta_mass,
        }
    }
}

/// The master list of all labels, in source order (source private
/// `fillLabelMasterList_`, `MultiplexDeltaMassesGenerator.cpp`).
pub fn label_master_list() -> Vec<Label> {
    vec![
        Label::new(
            "Arg6",
            "Label:13C(6)",
            "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188",
            6.0201290268,
        ),
        Label::new(
            "Arg10",
            "Label:13C(6)15N(4)",
            "Label:13C(6)15N(4)  |  C(-6) 13C(6) N(-4) 15N(4)  |  unimod #267",
            10.008268600,
        ),
        Label::new(
            "Lys4",
            "Label:2H(4)",
            "Label:2H(4)  |  H(-4) 2H(4)  |  unimod #481",
            4.0251069836,
        ),
        Label::new(
            "Lys6",
            "Label:13C(6)",
            "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188",
            6.0201290268,
        ),
        Label::new(
            "Lys8",
            "Label:13C(6)15N(2)",
            "Label:13C(6)15N(2)  |  C(-6) 13C(6) N(-2) 15N(2)  |  unimod #259",
            8.0141988132,
        ),
        Label::new(
            "Leu3",
            "Label:2H(3)",
            "Label:2H(3)  |  H(-3) 2H(3)  |  unimod #262",
            3.018830,
        ),
        Label::new(
            "Dimethyl0",
            "Dimethyl",
            "Dimethyl  |  H(4) C(2)  |  unimod #36",
            28.031300,
        ),
        Label::new(
            "Dimethyl4",
            "Dimethyl:2H(4)",
            "Dimethyl:2H(4)  |  2H(4) C(2)  |  unimod #199",
            32.056407,
        ),
        Label::new(
            "Dimethyl6",
            "Dimethyl:2H(4)13C(2)",
            "Dimethyl:2H(4)13C(2)  |  2H(4) 13C(2)  |  unimod #510",
            34.063117,
        ),
        Label::new(
            "Dimethyl8",
            "Dimethyl:2H(6)13C(2)",
            "Dimethyl:2H(6)13C(2)  |  H(-2) 2H(6) 13C(2)  |  unimod #330",
            36.075670,
        ),
        Label::new(
            "ICPL0",
            "ICPL",
            "ICPL  |  H(3) C(6) N O  |  unimod #365",
            105.021464,
        ),
        Label::new(
            "ICPL4",
            "ICPL:2H(4)",
            "ICPL:2H(4)  |  H(-1) 2H(4) C(6) N O  |  unimod #687",
            109.046571,
        ),
        Label::new(
            "ICPL6",
            "ICPL:13C(6)",
            "ICPL:13C(6)  |  H(3) 13C(6) N O  |  unimod #364",
            111.041593,
        ),
        Label::new(
            "ICPL10",
            "ICPL:13C(6)2H(4)",
            "ICPL:13C(6)2H(4)  |  H(-1) 2H(4) 13C(6) N O  |  unimod #866",
            115.066700,
        ),
    ]
}

/// The label mass shifts, the state of the source's default-constructed
/// generator as a `DefaultParamHandler("labels")`.
///
/// The source registers one float parameter per master-list label, named by
/// its short name, valued by its mass shift, described by its description and
/// restricted to `0:`. The source generator never reads these parameters
/// itself; `FeatureFinderMultiplexAlgorithm` embeds them as its `labels:`
/// section and passes the values to [`MultiplexDeltaMassesGenerator::new`].
#[derive(Clone, Debug, PartialEq)]
pub struct LabelMassShifts {
    /// Short label name to mass shift in Da, one entry per master-list label.
    pub masses: BTreeMap<String, f64>,
}

impl Default for LabelMassShifts {
    fn default() -> Self {
        Self {
            masses: label_master_list()
                .into_iter()
                .map(|l| (l.short_name, l.delta_mass))
                .collect(),
        }
    }
}

impl LabelMassShifts {
    /// The source defaults (source `MultiplexDeltaMassesGenerator()` followed
    /// by `getDefaults`), in master-list order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters carrying the current masses, in the layout of
    /// [`Self::defaults`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a mass is not in the master list,
    /// which has no parameter to carry it.
    pub fn to_param(&self) -> Result<Param> {
        for name in self.masses.keys() {
            if !label_master_list().iter().any(|l| &l.short_name == name) {
                return Err(Error::InvalidValue(format!(
                    "label '{name}' has no source parameter"
                )));
            }
        }
        let entries = label_master_list()
            .into_iter()
            .map(|label| {
                let value = self
                    .masses
                    .get(&label.short_name)
                    .copied()
                    .unwrap_or(label.delta_mass);
                ParamEntry {
                    name: label.short_name,
                    description: label.description,
                    value: ParamValue::Float(value),
                    min_float: 0.0,
                    ..ParamEntry::default()
                }
            })
            .collect();
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries,
            nodes: Vec::new(),
        })
    }

    /// The masses from `parameters`, missing ones at their defaults (source
    /// `setParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a value of the wrong type or below
    /// zero (the source's `Exception::InvalidParameter`).
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
        let mut handler = DefaultParamHandler::new(LABELS_PARAM_HANDLER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            let mut masses = BTreeMap::new();
            for label in label_master_list() {
                let value = param.value(&label.short_name)?.to_f64()?;
                masses.insert(label.short_name, value);
            }
            Ok(Self { masses })
        })
    }
}

/// The kind of labelling a `labels` string describes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Labelling {
    /// SILAC with Arg and/or Lys labels.
    Silac,
    /// Leucine labelling.
    Leu,
    /// Dimethyl labelling.
    Dimethyl,
    /// ICPL labelling.
    Icpl,
    /// Numeric mass shifts given directly as labels.
    Numeric,
    /// No labelling: singlet detection.
    None,
}

/// Generates the complete list of all possible mass shifts due to isotopic
/// labelling (source `MultiplexDeltaMassesGenerator`).
///
/// Isotopic labelling shifts peptide masses. In a Lys8/Arg10 SILAC sample,
/// peptides with one Arg show a mass shift of 10 Da between light and heavy
/// partners. This type constructs every mass shift pattern that the labels and
/// the permitted number of labelled amino acids can produce.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiplexDeltaMassesGenerator {
    /// The labels string (source `labels_`).
    labels: String,
    /// Flat list of all labels other than `no_label` (source `labels_list_`).
    labels_list: Vec<String>,
    /// The labels of each sample (source `samples_labels_`).
    samples_labels: Vec<Vec<String>>,
    /// Maximum number of labelled amino acids per peptide, minus one (source
    /// `max_nr_labelled_aas_`).
    max_nr_labelled_aas: i32,
    /// All mass shift patterns (source `delta_masses_list_`).
    delta_masses_list: Vec<MultiplexDeltaMasses>,
    /// Label to mass shift (source `label_delta_mass_`).
    label_delta_mass: BTreeMap<String, f64>,
    /// Short to long label (source `label_short_long_`).
    label_short_long: BTreeMap<String, String>,
    /// Long to short label (source `label_long_short_`).
    label_long_short: BTreeMap<String, String>,
    /// The detected kind of labelling.
    labelling: Labelling,
}

/// Source `boost::split(out, s, boost::is_any_of(separators))` with the
/// default `token_compress_off`: empty tokens between adjacent separators and
/// at either end are kept.
fn split_any(text: &str, separators: &str) -> Vec<String> {
    text.split(|c: char| separators.contains(c))
        .map(str::to_owned)
        .collect()
}

/// Outcome of `std::stod` on one label.
enum Stod {
    Value(f64),
    /// `std::invalid_argument`: no conversion could be performed.
    Invalid,
    /// `std::out_of_range`: `strtod` set `ERANGE`.
    OutOfRange,
}

/// `std::stod(text)`: `strtod` on the longest valid prefix after leading C
/// whitespace, throwing when nothing converts or the result is out of range.
///
/// Decimal numbers, `inf`/`infinity` and `nan` (optionally `nan(...)`) are
/// recognised case-insensitively, as glibc does. **Hexadecimal floats are not**:
/// `0x1p3` converts to `0` here (the prefix `0`), where glibc reads 8. Labels
/// are split at spaces, commas, colons and semicolons, so this is the only
/// `strtod` input class the port does not reproduce. Underflow to zero or a
/// subnormal counts as out of range, as glibc's `ERANGE` does.
fn stod(text: &str) -> Stod {
    let bytes = text.as_bytes();
    let mut start = 0;
    while start < bytes.len() && matches!(bytes[start], b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
    {
        start += 1;
    }
    let rest = &text[start..];
    let body = rest.strip_prefix(['+', '-']).unwrap_or(rest);
    let sign_len = rest.len() - body.len();
    let negative = rest.starts_with('-');
    let lower = body.to_ascii_lowercase();
    if lower.starts_with("inf") {
        return Stod::Value(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    if lower.starts_with("nan") {
        return Stod::Value(if negative { -f64::NAN } else { f64::NAN });
    }
    let b = body.as_bytes();
    let mut end = 0;
    let mut mantissa_digits = 0;
    let mut nonzero_digit = false;
    while end < b.len() && b[end].is_ascii_digit() {
        nonzero_digit |= b[end] != b'0';
        end += 1;
        mantissa_digits += 1;
    }
    if end < b.len() && b[end] == b'.' {
        end += 1;
        while end < b.len() && b[end].is_ascii_digit() {
            nonzero_digit |= b[end] != b'0';
            end += 1;
            mantissa_digits += 1;
        }
    }
    if mantissa_digits == 0 {
        return Stod::Invalid;
    }
    if end < b.len() && (b[end] == b'e' || b[end] == b'E') {
        let mut e = end + 1;
        if e < b.len() && (b[e] == b'+' || b[e] == b'-') {
            e += 1;
        }
        let digits_start = e;
        while e < b.len() && b[e].is_ascii_digit() {
            e += 1;
        }
        if e > digits_start {
            end = e;
        }
    }
    let literal = &rest[..sign_len + end];
    match literal.parse::<f64>() {
        Ok(value) if value.is_infinite() => Stod::OutOfRange,
        Ok(value) if nonzero_digit && (value == 0.0 || value.is_subnormal()) => Stod::OutOfRange,
        Ok(value) => Stod::Value(value),
        Err(_) => Stod::Invalid,
    }
}

impl MultiplexDeltaMassesGenerator {
    /// Parse `labels` and generate every mass shift pattern (source
    /// `MultiplexDeltaMassesGenerator(std::string labels, int
    /// max_nr_labelled_aas, std::map<std::string,double> label_mass_shift)`).
    ///
    /// `labels` describes the labels of each sample: `[...]` (or `(...)` or
    /// `{...}`) holds the labels of one sample. `[][Lys8,Arg10]` is a standard
    /// SILAC experiment: the light sample is unlabelled, and lysines and
    /// arginines of the heavy sample carry the labels `Lys8` and `Arg10`.
    /// Within a sample, labels are separated by `,`, `;`, `:` or a space, and
    /// empty labels between adjacent separators are kept, as `boost::split`
    /// keeps them. An empty sample is `no_label`; an empty string or no sample
    /// at all is one unlabelled sample (singlet detection).
    /// `max_nr_labelled_aas` is the maximum number of labelled amino acids per
    /// peptide minus one (the number of missed cleavages for tryptic SILAC), and
    /// `label_mass_shift` maps labels (e.g. `Lys8`) to their mass shifts (e.g.
    /// `8.0141988132`). A label missing from it contributes a zero shift, as the
    /// source's `operator[]` inserts zero.
    ///
    /// The labelling is SILAC when the string contains `Arg` or `Lys`, Leu for
    /// `Leu`, dimethyl for `Dimethyl` and ICPL for `ICPL`; exactly one of them
    /// may apply. Otherwise, when every label parses as a number (by
    /// `std::stod`, which accepts a numeric prefix: `3.5abc` is 3.5), the
    /// numbers are the mass shifts; a label that does not parse logs the
    /// source's warning "Unrecognized non-numeric label found. Assuming
    /// label-free." and then fails as below. Outside numeric labelling every
    /// label must be a **substring** of [`KNOWN_LABELS`], so `Arg`, `6 A` and
    /// the empty label are accepted, as in the source.
    ///
    /// Patterns: for SILAC, every combination of up to `max + 1` labelled
    /// arginines and lysines (their sum at most `max + 1`); for Leu, the second
    /// sample's first label `mc + 1` times against the first sample's first
    /// label, for `mc` in `0..=max`; for dimethyl, ICPL and numeric labelling,
    /// each sample's first label `mc + 1` times against the first sample's;
    /// without labelling, the single zero shift. The list is then sorted by
    /// [`MultiplexDeltaMasses::source_less`] with
    /// [`source_sort_by`](crate::math::source_sort::source_sort_by), which is
    /// the Release build's `std::sort`, so equivalent patterns end in the
    /// source's order.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] "Unknown labelling. Neither SILAC, Leu,
    ///   Dimethyl nor ICPL." when no or several labelling kinds apply (source
    ///   `Exception::IllegalArgument`), and "The label ... is unknown." for a
    ///   label outside [`KNOWN_LABELS`] (source `Exception::InvalidParameter`).
    /// - [`Error::InvalidValue`] for a negative `max_nr_labelled_aas`, which the
    ///   source casts to `unsigned` and loops over for billions of iterations.
    /// - [`Error::InvalidValue`] for Leu labelling with a single sample, where
    ///   the source reads a second sample that does not exist (undefined
    ///   behaviour).
    /// - [`Error::InvalidRange`] when more than [`MAX_DELTA_MASS_PATTERNS`]
    ///   patterns would be generated.
    pub fn new(
        labels: &str,
        max_nr_labelled_aas: i32,
        label_mass_shift: BTreeMap<String, f64>,
    ) -> Result<Self> {
        let max = u32::try_from(max_nr_labelled_aas)
            .map_err(|_| Error::InvalidValue("max_nr_labelled_aas must not be negative".into()))?;
        let masters = label_master_list();
        let mut label_short_long = BTreeMap::new();
        let mut label_long_short = BTreeMap::new();
        for label in &masters {
            label_short_long
                .entry(label.short_name.clone())
                .or_insert_with(|| label.long_name.clone());
            label_long_short
                .entry(label.long_name.clone())
                .or_insert_with(|| label.short_name.clone());
        }
        let mut label_delta_mass = label_mass_shift;

        // split the labels string
        let prepared = labels
            .replace("[]", "no_label")
            .replace("()", "no_label")
            .replace("{}", "no_label");
        let mut samples_labels: Vec<Vec<String>> = Vec::new();
        for sample in split_any(&prepared, "[](){}") {
            if sample.is_empty() {
                continue;
            }
            if sample == "no_label" {
                samples_labels.push(vec!["no_label".into()]);
            } else {
                samples_labels.push(split_any(&sample, ",;: "));
            }
        }
        if samples_labels.is_empty() {
            samples_labels.push(vec!["no_label".into()]);
        }

        // What kind of labelling do we have?
        let no_label = samples_labels.len() == 1
            && samples_labels[0].len() == 1
            && samples_labels[0][0] == "no_label";
        let labelling_silac = labels.contains("Arg") || labels.contains("Lys");
        let labelling_leu = labels.contains("Leu");
        let labelling_dimethyl = labels.contains("Dimethyl");
        let labelling_icpl = labels.contains("ICPL");
        let mut labelling_numeric = false;
        if !(no_label || labelling_silac || labelling_leu || labelling_dimethyl || labelling_icpl) {
            let mut all_numeric = true;
            for sample in &samples_labels {
                for label in sample {
                    match stod(label) {
                        Stod::Value(mass_shift) => {
                            label_delta_mass.entry(label.clone()).or_insert(mass_shift);
                            label_short_long
                                .entry(label.clone())
                                .or_insert_with(|| label.clone());
                            label_long_short
                                .entry(label.clone())
                                .or_insert_with(|| label.clone());
                        }
                        Stod::Invalid | Stod::OutOfRange => {
                            // Source: OPENMS_LOG_WARN "Unrecognized non-numeric
                            // label found. Assuming label-free." and leave this
                            // sample's loop.
                            all_numeric = false;
                            break;
                        }
                    }
                }
            }
            labelling_numeric = all_numeric;
        }
        let labelling_none =
            labels.is_empty() || labels == "[]" || labels == "()" || labels == "{}";
        let flags = [
            labelling_silac,
            labelling_leu,
            labelling_dimethyl,
            labelling_icpl,
            labelling_numeric,
            labelling_none,
        ];
        if flags.iter().filter(|&&f| f).count() != 1 {
            return Err(Error::InvalidValue(
                "Unknown labelling. Neither SILAC, Leu, Dimethyl nor ICPL.".into(),
            ));
        }
        let labelling = match flags.iter().position(|&f| f) {
            Some(0) => Labelling::Silac,
            Some(1) => Labelling::Leu,
            Some(2) => Labelling::Dimethyl,
            Some(3) => Labelling::Icpl,
            Some(4) => Labelling::Numeric,
            _ => Labelling::None,
        };

        // Check that the labels are known, unless the labelling is numeric.
        if labelling != Labelling::Numeric {
            for label in samples_labels.iter().flatten() {
                if !KNOWN_LABELS.contains(label.as_str()) {
                    return Err(Error::InvalidValue(format!(
                        "The label {label} is unknown."
                    )));
                }
            }
        }

        let max = max as usize;
        let pattern_count = match labelling {
            Labelling::Silac => (max + 2).checked_mul(max + 3).map(|n| n / 2),
            Labelling::None => Some(1),
            _ => max.checked_add(1),
        };
        if pattern_count.is_none_or(|n| n > MAX_DELTA_MASS_PATTERNS) {
            return Err(Error::InvalidRange(format!(
                "more than {MAX_DELTA_MASS_PATTERNS} mass shift patterns"
            )));
        }
        let mass =
            |map: &BTreeMap<String, f64>, label: &str| map.get(label).copied().unwrap_or(0.0);

        let mut delta_masses_list = Vec::new();
        match labelling {
            Labelling::Silac => {
                // The first sample is assumed unlabelled, even when its "[]" is
                // missing from the labels string.
                for arg_per_peptide in 0..=(max + 1) {
                    for lys_per_peptide in 0..=(max + 1) {
                        if arg_per_peptide + lys_per_peptide > max + 1 {
                            continue;
                        }
                        let mut pattern = MultiplexDeltaMasses::new();
                        pattern
                            .delta_masses
                            .push(DeltaMass::with_label(0.0, "no_label"));
                        for sample in &samples_labels {
                            let mut mass_shift = 0.0;
                            let mut label_set = LabelSet::new();
                            // An amino acid counted in the peptide but without
                            // a label in this sample must not produce a shift.
                            let mut go_ahead_lys = false;
                            let mut go_ahead_arg = false;
                            for label in sample {
                                let arg6 = label.contains("Arg6");
                                let arg10 = label.contains("Arg10");
                                let lys4 = label.contains("Lys4");
                                let lys6 = label.contains("Lys6");
                                let lys8 = label.contains("Lys8");
                                for (present, name, count) in [
                                    (arg6, "Arg6", arg_per_peptide),
                                    (arg10, "Arg10", arg_per_peptide),
                                    (lys4, "Lys4", lys_per_peptide),
                                    (lys6, "Lys6", lys_per_peptide),
                                    (lys8, "Lys8", lys_per_peptide),
                                ] {
                                    if present {
                                        for _ in 0..count {
                                            label_set.insert(name);
                                        }
                                    }
                                }
                                let b = |x: bool| if x { 1.0 } else { 0.0 };
                                let arg = b(arg6) * mass(&label_delta_mass, "Arg6")
                                    + b(arg10) * mass(&label_delta_mass, "Arg10");
                                let lys = b(lys4) * mass(&label_delta_mass, "Lys4")
                                    + b(lys6) * mass(&label_delta_mass, "Lys6")
                                    + b(lys8) * mass(&label_delta_mass, "Lys8");
                                mass_shift = mass_shift
                                    + arg_per_peptide as f64 * arg
                                    + lys_per_peptide as f64 * lys;
                                go_ahead_arg =
                                    go_ahead_arg || !(arg_per_peptide != 0 && !arg6 && !arg10);
                                go_ahead_lys = go_ahead_lys
                                    || !(lys_per_peptide != 0 && !lys4 && !lys6 && !lys8);
                            }
                            if go_ahead_arg && go_ahead_lys && mass_shift != 0.0 {
                                pattern
                                    .delta_masses
                                    .push(DeltaMass::new(mass_shift, label_set));
                            }
                        }
                        if pattern.delta_masses.len() > 1 {
                            delta_masses_list.push(pattern);
                        }
                    }
                }
            }
            Labelling::Leu => {
                // Each sample is labelled once; only the first label counts.
                if samples_labels.len() < 2 {
                    return Err(Error::InvalidValue(
                        "Leu labelling needs a second sample (the source reads past the sample list)"
                            .into(),
                    ));
                }
                let heavy = samples_labels[1][0].clone();
                let light = samples_labels[0][0].clone();
                for mc in 0..=max {
                    let mut pattern = MultiplexDeltaMasses::new();
                    pattern
                        .delta_masses
                        .push(DeltaMass::with_label(0.0, "no_label"));
                    let mass_shift = (mc + 1) as f64
                        * (mass(&label_delta_mass, &heavy) - mass(&label_delta_mass, &light));
                    let label_set: LabelSet = std::iter::repeat_n(heavy.clone(), mc + 1).collect();
                    pattern
                        .delta_masses
                        .push(DeltaMass::new(mass_shift, label_set));
                    delta_masses_list.push(pattern);
                }
            }
            Labelling::Dimethyl | Labelling::Icpl | Labelling::Numeric => {
                // Each sample is labelled once; only the first label counts.
                let first = samples_labels[0][0].clone();
                for mc in 0..=max {
                    let mut pattern = MultiplexDeltaMasses::new();
                    for sample in &samples_labels {
                        let mass_shift = (mc + 1) as f64
                            * (mass(&label_delta_mass, &sample[0])
                                - mass(&label_delta_mass, &first));
                        let label_set: LabelSet =
                            std::iter::repeat_n(sample[0].clone(), mc + 1).collect();
                        pattern
                            .delta_masses
                            .push(DeltaMass::new(mass_shift, label_set));
                    }
                    delta_masses_list.push(pattern);
                }
            }
            Labelling::None => {
                let mut pattern = MultiplexDeltaMasses::new();
                pattern
                    .delta_masses
                    .push(DeltaMass::with_label(0.0, "no_label"));
                delta_masses_list.push(pattern);
            }
        }

        // from few miscleavages (simple explanations) to many
        source_sort_by(&mut delta_masses_list, |a, b| a.source_less(b))?;

        let labels_list = samples_labels
            .iter()
            .flatten()
            .filter(|l| l.as_str() != "no_label")
            .cloned()
            .collect();

        Ok(Self {
            labels: labels.to_owned(),
            labels_list,
            samples_labels,
            max_nr_labelled_aas,
            delta_masses_list,
            label_delta_mass,
            label_short_long,
            label_long_short,
            labelling,
        })
    }

    /// Add the patterns that arise when one or several peptides of the
    /// multiplet are absent (source `generateKnockoutDeltaMasses`): for a
    /// triplet experiment the doublets and singlets, for a quadruplet the
    /// triplets, doublets and singlets.
    ///
    /// The multiplicity `n` is the length of the first pattern. For doublets a
    /// single singlet with the label `any_label_set` is added. For triplets,
    /// each original pattern adds the doublets (0,1), (1,2) and (0,2), then one
    /// `any_label_set` singlet follows. For quadruplets, each adds the triplets
    /// (1,2,3), (0,2,3) and (0,1,2) — the source comments out (0,1,3) as
    /// "previously forgotten", and so does this port — and the doublets (0,1),
    /// (0,2), (0,3), (1,2), (1,3), (2,3), then one singlet. The list is sorted
    /// again as in [`Self::new`]. A first pattern with no mass shift adds
    /// nothing and only sorts, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an empty list (source
    /// `Exception::InvalidSize`), for singlets ("Knock-outs for singlet
    /// detection not relevant.") and for more than four samples ("Knock-outs
    /// for multiplex experiments with more than 4 samples not supported."),
    /// both source `Exception::IllegalArgument`; and when an original pattern
    /// is shorter than the first, where the source reads past its end. On
    /// error the list is unchanged.
    pub fn generate_knockout_delta_masses(&mut self) -> Result<()> {
        let Some(first) = self.delta_masses_list.first() else {
            return Err(Error::InvalidValue(
                "delta_masses_list_ must not be empty".into(),
            ));
        };
        let n = first.delta_masses.len();
        let combinations: &[&[usize]] = match n {
            0 => &[],
            1 => {
                return Err(Error::InvalidValue(
                    "Knock-outs for singlet detection not relevant.".into(),
                ));
            }
            2 => &[],
            3 => &[&[0, 1], &[1, 2], &[0, 2]],
            4 => &[
                &[1, 2, 3],
                &[0, 2, 3],
                &[0, 1, 2],
                &[0, 1],
                &[0, 2],
                &[0, 3],
                &[1, 2],
                &[1, 3],
                &[2, 3],
            ],
            _ => {
                return Err(Error::InvalidValue(
                    "Knock-outs for multiplex experiments with more than 4 samples not supported."
                        .into(),
                ));
            }
        };
        let mut list = self.delta_masses_list.clone();
        if n >= 2 {
            for original in &self.delta_masses_list {
                for combination in combinations {
                    let mut pattern = MultiplexDeltaMasses::new();
                    for &k in *combination {
                        let delta = original.delta_masses.get(k).ok_or_else(|| {
                            Error::InvalidValue(
                                "a mass shift pattern is shorter than the first (the source reads past its end)"
                                    .into(),
                            )
                        })?;
                        pattern.delta_masses.push(delta.clone());
                    }
                    list.push(pattern);
                }
            }
            let mut singlet = MultiplexDeltaMasses::new();
            singlet
                .delta_masses
                .push(DeltaMass::with_label(0.0, "any_label_set"));
            list.push(singlet);
        }
        source_sort_by(&mut list, |a, b| a.source_less(b))?;
        self.delta_masses_list = list;
        Ok(())
    }

    /// The labels of each sample, one line per sample, as source
    /// `printSamplesLabelsList` writes them: a leading newline, then
    /// `sample <i>:    ` and every label followed by four spaces.
    pub fn samples_labels_text(&self) -> String {
        let mut text = String::from("\n");
        for (i, sample) in self.samples_labels.iter().enumerate() {
            text.push_str(&format!("sample {}:    ", i + 1));
            for label in sample {
                text.push_str(label);
                text.push_str("    ");
            }
            text.push('\n');
        }
        text
    }

    /// Write [`Self::samples_labels_text`] to `stream` (source
    /// `printSamplesLabelsList(std::ostream&)`).
    ///
    /// # Errors
    ///
    /// The stream's write error.
    pub fn print_samples_labels_list(&self, stream: &mut impl std::io::Write) -> Result<()> {
        stream.write_all(self.samples_labels_text().as_bytes())?;
        Ok(())
    }

    /// All mass shift patterns, one line per pattern, as source
    /// `printDeltaMassesList` writes them: a leading newline, then
    /// `mass shift <i>:    ` and each shift as `<mass> (<labels>)    `, the mass
    /// formatted as `std::ostream << double` does (six significant digits) and
    /// the labels comma-separated; a trailing empty line.
    pub fn delta_masses_text(&self) -> String {
        let mut text = String::from("\n");
        for (i, pattern) in self.delta_masses_list.iter().enumerate() {
            text.push_str(&format!("mass shift {}:    ", i + 1));
            for delta in &pattern.delta_masses {
                text.push_str(&ostream_g(delta.delta_mass, 6));
                text.push_str(" (");
                text.push_str(&delta.label_set.as_slice().join(","));
                text.push_str(")    ");
            }
            text.push('\n');
        }
        text.push('\n');
        text
    }

    /// Write [`Self::delta_masses_text`] to `stream` (source
    /// `printDeltaMassesList(std::ostream&)`).
    ///
    /// # Errors
    ///
    /// The stream's write error.
    pub fn print_delta_masses_list(&self, stream: &mut impl std::io::Write) -> Result<()> {
        stream.write_all(self.delta_masses_text().as_bytes())?;
        Ok(())
    }

    /// The list of mass shift patterns (source `getDeltaMassesList`, both
    /// overloads).
    pub fn delta_masses_list(&self) -> &[MultiplexDeltaMasses] {
        &self.delta_masses_list
    }

    /// The labels of each sample (source `getSamplesLabelsList`, both
    /// overloads). In a standard SILAC experiment: `[["no_label"], ["Lys8",
    /// "Arg10"]]`.
    pub fn samples_labels_list(&self) -> &[Vec<String>] {
        &self.samples_labels
    }

    /// The short label for a long (UniMod) label such as `Label:13C(6)15N(4)`
    /// (source `getLabelShort`).
    ///
    /// `Label:13C(6)` is the long name of both `Arg6` and `Lys6`; the source
    /// map keeps the first insertion, so it is always `Arg6`. An unknown label
    /// is `None`, where the source's `operator[]` inserts and returns an empty
    /// string.
    pub fn label_short(&self, label: &str) -> Option<&str> {
        self.label_long_short.get(label).map(String::as_str)
    }

    /// The long (UniMod) label for a short label such as `Arg10` (source
    /// `getLabelLong`). An unknown label is `None`, where the source returns an
    /// empty string.
    pub fn label_long(&self, label: &str) -> Option<&str> {
        self.label_short_long.get(label).map(String::as_str)
    }

    /// The label set of `sequence` (source `extractLabelSet`), from its string
    /// form; see [`Self::extract_label_set_from_str`].
    pub fn extract_label_set(&self, sequence: &AASequence) -> LabelSet {
        self.extract_label_set_from_str(&sequence.to_string())
    }

    /// The label set of a sequence string: for every label of the samples, the
    /// number of occurrences of `(<long label>)`, removed from the string as
    /// they are counted, in label order. `VLSEEEIDDNFK(Label:13C(6)15N(2))AQR(Label:13C(6)15N(4))`
    /// holds `Lys8` and `Arg10`. A set with no label found is `{no_label}`,
    /// whether or not `no_label` is among the samples' labels.
    ///
    /// Because occurrences are removed, a label counted earlier hides later
    /// labels sharing its text, and a label without a long name searches for
    /// `()`, as the source's empty-string lookup does.
    pub fn extract_label_set_from_str(&self, sequence: &str) -> LabelSet {
        let mut s = sequence.to_owned();
        let mut label_set = LabelSet::new();
        for label in &self.labels_list {
            let pattern = format!("({})", self.label_long(label).unwrap_or(""));
            if s.contains(&pattern) {
                let before = s.len();
                s = s.replace(&pattern, "");
                let multiple = (before - s.len()) / pattern.len();
                for _ in 0..multiple {
                    label_set.insert(label.clone());
                }
            }
        }
        if label_set.is_empty() {
            label_set.insert("no_label");
        }
        label_set
    }

    /// The labels string this generator was built from (source `labels_`).
    pub fn labels(&self) -> &str {
        &self.labels
    }

    /// The flat list of all labels other than `no_label` (source
    /// `labels_list_`).
    pub fn labels_list(&self) -> &[String] {
        &self.labels_list
    }

    /// Maximum number of labelled amino acids per peptide, minus one (source
    /// `max_nr_labelled_aas_`).
    pub fn max_nr_labelled_aas(&self) -> i32 {
        self.max_nr_labelled_aas
    }

    /// Label to mass shift, including the numeric labels the constructor added
    /// (source `label_delta_mass_`).
    pub fn label_delta_mass(&self) -> &BTreeMap<String, f64> {
        &self.label_delta_mass
    }

    /// The detected kind of labelling.
    pub fn labelling(&self) -> Labelling {
        self.labelling
    }
}
