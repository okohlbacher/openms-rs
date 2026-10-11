// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Selection of MRM features by linear programming.
//!
//! Ports `ANALYSIS/OPENSWATH/MRMFeatureSelector.h` and its `.cpp`: the
//! abstract selector [`MRMFeatureSelector`], its two implementations
//! [`MRMFeatureSelectorQMIP`] and [`MRMFeatureSelectorScore`], the
//! [`SelectorParameters`] they share and the [`LambdaScore`] weightings.
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_FEATURE_SELECTOR_SUPPORT.md).
//!
//! # Algorithm
//!
//! 1. The features are grouped by component group (or component) name and
//!    sorted by assay retention time.
//! 2. The sorted list is cut into overlapping segments of
//!    [`SelectorParameters::segment_window_length`] entries, advancing by
//!    [`SelectorParameters::segment_step_length`].
//! 3. Each segment becomes a linear program that is handed to a solver; the
//!    columns above [`SelectorParameters::optimal_threshold`] are the selected
//!    features.
//! 4. The selections of all segments are merged.
//!
//! # The linear-programming solver
//!
//! The source builds its problems in `LPWrapper` (`DATASTRUCTURES/LPWrapper.h`)
//! and solves them with GLPK or COIN-OR, which this crate does not provide:
//! the solver choice for `LPWrapper` is a recorded backlog item. The problem
//! construction is ported in full against [`LinearProgram`], a minimal model
//! of the problem-building half of `LPWrapper`, and the solve is a call to an
//! [`LpSolver`]. [`PartitionSolver`] solves exactly the problems of
//! [`MRMFeatureSelectorScore`]; [`MRMFeatureSelectorQMIP`] needs a
//! general solver supplied by the caller.
//!
//! # Example
//!
//! ```
//! use openms::analysis::openswath::mrm_feature_selector::{
//!     LambdaScore, MRMFeatureSelector, MRMFeatureSelectorScore, PartitionSolver,
//!     SelectorParameters,
//! };
//! use openms::kernel::{Feature, FeatureMap};
//!
//! fn candidate(id: u64, sn: f64) -> Feature {
//!     let mut f = Feature::default();
//!     f.unique_id = id;
//!     f.metadata.insert("PeptideRef".into(), "peptide".into());
//!     f.metadata.insert("assay_rt".into(), 1.0_f64.try_into().unwrap());
//!     f.metadata.insert("sn_ratio".into(), sn.try_into().unwrap());
//!     let mut sub = Feature::default();
//!     sub.metadata.insert("native_id".into(), "transition".into());
//!     f.subordinates.push(sub);
//!     f
//! }
//!
//! let features = FeatureMap::from_features(vec![candidate(1, 20.0), candidate(2, 5.0)]);
//! let mut parameters = SelectorParameters::default();
//! parameters.segment_window_length = -1;
//! parameters.segment_step_length = -1;
//! parameters.score_weights.insert("sn_ratio".into(), LambdaScore::Linear);
//! let selector = MRMFeatureSelectorScore::new(PartitionSolver);
//! let selected = selector.select_mrm_feature(&features, &parameters).unwrap();
//! // the score is minimised: the smaller sn_ratio wins
//! assert_eq!(selected.features.len(), 1);
//! assert_eq!(selected.features[0].unique_id, 2);
//! ```
//!
//! [`MRMFeatureSelector`]: crate::analysis::openswath::mrm_feature_selector::MRMFeatureSelector
//! [`MRMFeatureSelectorQMIP`]: crate::analysis::openswath::mrm_feature_selector::MRMFeatureSelectorQMIP
//! [`MRMFeatureSelectorScore`]: crate::analysis::openswath::mrm_feature_selector::MRMFeatureSelectorScore
//! [`SelectorParameters`]: crate::analysis::openswath::mrm_feature_selector::SelectorParameters
//! [`LambdaScore`]: crate::analysis::openswath::mrm_feature_selector::LambdaScore
//! [`SelectorParameters::segment_window_length`]: crate::analysis::openswath::mrm_feature_selector::SelectorParameters::segment_window_length
//! [`SelectorParameters::segment_step_length`]: crate::analysis::openswath::mrm_feature_selector::SelectorParameters::segment_step_length
//! [`SelectorParameters::optimal_threshold`]: crate::analysis::openswath::mrm_feature_selector::SelectorParameters::optimal_threshold
//! [`LinearProgram`]: crate::analysis::openswath::mrm_feature_selector::LinearProgram
//! [`LpSolver`]: crate::analysis::openswath::mrm_feature_selector::LpSolver
//! [`PartitionSolver`]: crate::analysis::openswath::mrm_feature_selector::PartitionSolver

use crate::analysis::openswath::mrm_feature_support::{meta_double, meta_string};
use crate::kernel::{Feature, FeatureMap};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// Whether the selector's variables are integers or continuous.
///
/// Source `MRMFeatureSelector::VariableType`; the values are the source's
/// `1` and `2`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum VariableType {
    /// Integer variables, an exact selection.
    Integer = 1,
    /// Continuous variables, the relaxed problem (the source default).
    #[default]
    Continuous = 2,
}

/// The transformation applied to a score before it enters the objective.
///
/// Source `MRMFeatureSelector::LambdaScore`; the values are the source's
/// `1` to `5`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum LambdaScore {
    /// The score itself.
    Linear = 1,
    /// `1 / score`.
    Inverse = 2,
    /// `ln(score)`.
    Log = 3,
    /// `1 / ln(score)`.
    InverseLog = 4,
    /// `1 / log10(score)`.
    InverseLog10 = 5,
}

/// The parameters of a selector call (source
/// `MRMFeatureSelector::SelectorParameters`).
#[derive(Clone, Debug, PartialEq)]
pub struct SelectorParameters {
    /// Nearest neighbour threshold: the number of components or component
    /// groups to the left and right to include in the optimisation problem.
    pub nn_threshold: i32,
    /// Weight compounds with a nearer retention time greater than compounds
    /// with a further one.
    pub locality_weight: bool,
    /// Use component groups instead of components for the retention time
    /// optimisation.
    pub select_transition_group: bool,
    /// Number of components or component groups in a segment. `-1` together
    /// with a step of `-1` means a single segment of everything.
    pub segment_window_length: i32,
    /// Number of components or component groups the window advances by.
    pub segment_step_length: i32,
    /// Integer or continuous variables.
    pub variable_type: VariableType,
    /// Value at or above which a variable counts as selected, `0 < x < 1`.
    pub optimal_threshold: f64,
    /// Metavalue name to the transformation applied to it.
    pub score_weights: BTreeMap<String, LambdaScore>,
}

impl Default for SelectorParameters {
    fn default() -> Self {
        Self {
            nn_threshold: 4,
            locality_weight: false,
            select_transition_group: true,
            segment_window_length: 8,
            segment_step_length: 4,
            variable_type: VariableType::Continuous,
            optimal_threshold: 0.5,
            score_weights: BTreeMap::new(),
        }
    }
}

/// Whether a column is bounded, source `LPWrapper::VariableType` values used
/// by the selectors (`UNBOUNDED`, `DOUBLE_BOUNDED`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundType {
    /// No bounds.
    Unbounded,
    /// Only a lower bound.
    LowerBoundOnly,
    /// Only an upper bound.
    UpperBoundOnly,
    /// Both bounds.
    DoubleBounded,
    /// Both bounds, equal.
    Fixed,
}

/// Source `LPWrapper::VariableType`: whether a column is continuous or integer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ColumnType {
    /// A real-valued column.
    Continuous,
    /// An integer-valued column.
    Integer,
}

/// A column (variable) of a [`LinearProgram`].
#[derive(Clone, Debug, PartialEq)]
pub struct Column {
    /// The column name.
    pub name: String,
    /// The bound kind.
    pub bound_type: BoundType,
    /// Lower bound; meaningful for the bound kinds that have one.
    pub lower: f64,
    /// Upper bound.
    pub upper: f64,
    /// Integer or continuous.
    pub column_type: ColumnType,
    /// The objective coefficient.
    pub objective: f64,
}

/// A row (constraint) of a [`LinearProgram`].
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    /// The row name.
    pub name: String,
    /// Column indices of the nonzero coefficients.
    pub indices: Vec<usize>,
    /// The coefficients, one per index.
    pub values: Vec<f64>,
    /// Lower bound.
    pub lower: f64,
    /// Upper bound.
    pub upper: f64,
    /// Which bounds apply.
    pub bound_type: BoundType,
}

/// The objective sense (source `LPWrapper::Sense`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sense {
    /// Minimise the objective.
    Min,
    /// Maximise the objective.
    Max,
}

/// A linear program as the selectors build it: the problem-building half of
/// source `LPWrapper`.
#[derive(Clone, Debug, PartialEq)]
pub struct LinearProgram {
    /// The objective sense.
    pub sense: Sense,
    /// The columns, in creation order.
    pub columns: Vec<Column>,
    /// The rows, in creation order.
    pub rows: Vec<Row>,
}

impl LinearProgram {
    /// An empty program, minimising.
    pub fn new() -> Self {
        Self {
            sense: Sense::Min,
            columns: Vec::new(),
            rows: Vec::new(),
        }
    }

    /// Index of the column with `name`, the first one if the name repeats.
    pub fn column_index(&self, name: &str) -> Option<usize> {
        self.columns.iter().position(|c| c.name == name)
    }
}

impl Default for LinearProgram {
    fn default() -> Self {
        Self::new()
    }
}

/// How a solve ended (source `LPWrapper::SolverStatus`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolverStatus {
    /// An optimal solution was found.
    Optimal,
    /// A feasible, not proven optimal, solution was found.
    Feasible,
    /// The problem is infeasible.
    Infeasible,
    /// The problem is unbounded.
    Unbounded,
    /// The solver failed or its outcome is undefined.
    Undefined,
}

/// The outcome of [`LpSolver::solve`]: a status and one value per column.
#[derive(Clone, Debug, PartialEq)]
pub struct LpSolution {
    /// How the solve ended.
    pub status: SolverStatus,
    /// The value of each column, indexed like [`LinearProgram::columns`].
    pub column_values: Vec<f64>,
}

/// A linear (mixed integer) programming solver.
pub trait LpSolver {
    /// Solve `problem`.
    ///
    /// # Errors
    ///
    /// An implementation may report a failure to run as an [`Error`]; a
    /// problem that has no solution is a [`SolverStatus`], not an error.
    fn solve(&self, problem: &LinearProgram) -> Result<LpSolution>;
}

/// An exact solver for programs whose rows are disjoint, fixed to one and
/// whose columns are bounded to `[0, 1]`: choose, in each row, the column with
/// the smallest objective.
///
/// These are the programs of [`MRMFeatureSelectorScore`]. For such a program
/// the continuous and the integer optimum coincide. When several columns of a
/// row have the same smallest objective the **last** is chosen: the source
/// leaves that choice to the solver's pivoting, and the Release build's GLPK
/// was measured to take the last of two and of three tied columns on the
/// upstream data in six of the seven parameter sets of `tests/mrm_oracle.rs`;
/// in the seventh (continuous, overlapping segments) it took different tied
/// columns in the two segments, which a one-row-at-a-time solver cannot
/// reproduce. That is an observation about one solver, not a guarantee.
/// A program of any other shape is refused with [`Error::Unsupported`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PartitionSolver;

impl LpSolver for PartitionSolver {
    fn solve(&self, problem: &LinearProgram) -> Result<LpSolution> {
        if problem.sense != Sense::Min {
            return Err(Error::Unsupported(
                "PartitionSolver only minimises".to_owned(),
            ));
        }
        let mut owner: Vec<Option<usize>> = vec![None; problem.columns.len()];
        for (r, row) in problem.rows.iter().enumerate() {
            let fixed_to_one = row.bound_type == BoundType::Fixed
                && row.lower == 1.0
                && row.upper == 1.0
                && row.values.iter().all(|&v| v == 1.0);
            if !fixed_to_one {
                return Err(Error::Unsupported(
                    "PartitionSolver solves only rows of the form sum(x) = 1".to_owned(),
                ));
            }
            for &c in &row.indices {
                match owner.get_mut(c) {
                    Some(slot @ None) => *slot = Some(r),
                    _ => {
                        return Err(Error::Unsupported(
                            "PartitionSolver needs every column in exactly one row".to_owned(),
                        ));
                    }
                }
            }
        }
        for column in &problem.columns {
            if column.bound_type != BoundType::DoubleBounded
                || column.lower != 0.0
                || column.upper != 1.0
            {
                return Err(Error::Unsupported(
                    "PartitionSolver needs every column bounded to [0, 1]".to_owned(),
                ));
            }
        }
        let mut values = vec![0.0; problem.columns.len()];
        for row in &problem.rows {
            let mut best: Option<usize> = None;
            for &c in &row.indices {
                let better = match best {
                    None => true,
                    Some(b) => problem.columns[c].objective <= problem.columns[b].objective,
                };
                if better {
                    best = Some(c);
                }
            }
            match best {
                Some(c) => values[c] = 1.0,
                None => {
                    return Ok(LpSolution {
                        status: SolverStatus::Infeasible,
                        column_values: values,
                    });
                }
            }
        }
        Ok(LpSolution {
            status: SolverStatus::Optimal,
            column_values: values,
        })
    }
}

/// A selector of MRM features (source abstract class `MRMFeatureSelector`).
///
/// An implementation provides [`optimize`](Self::optimize); the segmenting and
/// merging of [`select_mrm_feature`](Self::select_mrm_feature) is shared.
pub trait MRMFeatureSelector {
    /// Set up the linear program of one segment, solve it and return the names
    /// of the variables it selects.
    ///
    /// `time_to_name` holds (retention time, component or component group
    /// name) pairs; `feature_name_map` maps a name to its candidate features.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a name of `time_to_name` is missing from
    /// `feature_name_map` or the solver finds no feasible selection (the
    /// source's `Exception::FailedAPICall`).
    fn optimize(
        &self,
        time_to_name: &[(f64, String)],
        feature_name_map: &BTreeMap<String, Vec<Feature>>,
        parameters: &SelectorParameters,
    ) -> Result<Vec<String>>;

    /// Select features: sort them by retention time, split them into segments
    /// of `segment_window_length` entries advancing by `segment_step_length`,
    /// [`optimize`](Self::optimize) each segment and keep the subordinates
    /// whose name was selected.
    ///
    /// Features without a selected subordinate are dropped; the others are
    /// returned with only their selected subordinates. When both lengths are
    /// `-1` there is one segment holding everything. An empty input gives an
    /// empty output.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `segment_step_length` is `0` (an integer
    /// division by zero in the source), when `assay_rt` of a feature is not a
    /// number, and the errors of [`optimize`](Self::optimize).
    fn select_mrm_feature(
        &self,
        features: &FeatureMap,
        parameters: &SelectorParameters,
    ) -> Result<FeatureMap> {
        let mut selected_filtered = FeatureMap::default();
        if features.is_empty() {
            return Ok(selected_filtered);
        }
        let (mut time_to_name, feature_name_map) =
            construct_targ_trans_list(features, parameters.select_transition_group)?;
        time_to_name.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

        let mut window_length = i64::from(parameters.segment_window_length);
        let mut step_length = i64::from(parameters.segment_step_length);
        if window_length == -1 && step_length == -1 {
            window_length = time_to_name.len() as i64;
            step_length = window_length;
        }
        // The source computes with `Size` (unsigned 64-bit): a negative `Int`
        // converts to a huge value. Reproduce that, and refuse the division by
        // zero it would crash on.
        let step = step_length as u64;
        if step == 0 {
            return Err(Error::InvalidValue(
                "segment_step_length is 0; the source divides by it".to_owned(),
            ));
        }
        let size = time_to_name.len() as u64;
        let mut n_segments = size / step;
        if size % step != 0 {
            n_segments += 1;
        }
        let window = window_length as u64;
        let mut result_names: Vec<String> = Vec::new();
        for i in 0..n_segments {
            let start = step.wrapping_mul(i);
            let end = start.wrapping_add(window).min(size);
            if start > end {
                return Err(Error::InvalidValue(
                    "segment window and step lengths select an empty or inverted slice".to_owned(),
                ));
            }
            let slice = &time_to_name[start as usize..end as usize];
            result_names.extend(self.optimize(slice, &feature_name_map, parameters)?);
        }
        let result_names: BTreeSet<String> = result_names.into_iter().collect();
        for feature in &features.features {
            let mut subordinates_filtered = Vec::new();
            for subordinate in &feature.subordinates {
                let feature_name = if parameters.select_transition_group {
                    format!(
                        "{}_{}",
                        remove_spaces(&meta_string(&feature.metadata, "PeptideRef")),
                        feature.unique_id
                    )
                } else {
                    format!(
                        "{}_{}",
                        remove_spaces(&meta_string(&subordinate.metadata, "native_id")),
                        feature.unique_id
                    )
                };
                if result_names.contains(&feature_name) {
                    subordinates_filtered.push(subordinate.clone());
                }
            }
            if !subordinates_filtered.is_empty() {
                let mut feature_filtered = feature.clone();
                feature_filtered.subordinates = subordinates_filtered;
                selected_filtered.features.push(feature_filtered);
            }
        }
        Ok(selected_filtered)
    }
}

/// Remove the space characters of `text` (source `removeSpaces_`).
pub fn remove_spaces(text: &str) -> String {
    text.chars().filter(|&c| c != ' ').collect()
}

/// What [`construct_targ_trans_list`] returns: the `(assay retention time,
/// name)` pairs in feature order and, per name, the candidate features.
pub type TargTransList = (Vec<(f64, String)>, BTreeMap<String, Vec<Feature>>);

/// Build the retention time / name list and the name to candidate map (source
/// `constructTargTransList_`).
///
/// A feature is filed under its `PeptideRef` (spaces removed), once in
/// `time_to_name` with its `assay_rt`. With `select_transition_group` false its
/// subordinates are also filed under their `native_id`; the source tests
/// `names.contains(component_name)` where it plainly means the negation, so a
/// component name is added to `time_to_name` only when it has been seen before
/// and then again for every further occurrence. That is reproduced.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `assay_rt` is missing or not a number.
pub fn construct_targ_trans_list(
    features: &FeatureMap,
    select_transition_group: bool,
) -> Result<TargTransList> {
    let mut time_to_name: Vec<(f64, String)> = Vec::new();
    let mut feature_name_map: BTreeMap<String, Vec<Feature>> = BTreeMap::new();
    let mut names: BTreeSet<String> = BTreeSet::new();
    for feature in &features.features {
        let component_group_name = remove_spaces(&meta_string(&feature.metadata, "PeptideRef"));
        let assay_retention_time = meta_double(&feature.metadata, "assay_rt")?;
        if !names.contains(&component_group_name) {
            time_to_name.push((assay_retention_time, component_group_name.clone()));
            names.insert(component_group_name.clone());
        }
        feature_name_map
            .entry(component_group_name)
            .or_default()
            .push(feature.clone());
        if select_transition_group {
            continue;
        }
        for subordinate in &feature.subordinates {
            let component_name = remove_spaces(&meta_string(&subordinate.metadata, "native_id"));
            if names.contains(&component_name) {
                time_to_name.push((assay_retention_time, component_name.clone()));
                names.insert(component_name.clone());
            }
            feature_name_map
                .entry(component_name)
                .or_default()
                .push(subordinate.clone());
        }
    }
    Ok((time_to_name, feature_name_map))
}

/// Transform `score` with `lambda_score` (source `weightScore_`).
///
/// The source throws `IllegalArgument` for a value outside the enumeration;
/// every [`LambdaScore`] is handled here, so no error exists.
pub fn weight_score(score: f64, lambda_score: LambdaScore) -> f64 {
    match lambda_score {
        LambdaScore::Linear => score,
        LambdaScore::Inverse => 1.0 / score,
        LambdaScore::Log => score.ln(),
        LambdaScore::InverseLog => 1.0 / score.ln(),
        LambdaScore::InverseLog10 => 1.0 / score.log10(),
    }
}

/// The score of a feature under the given weights (source `computeScore_`).
///
/// The product of the transformed metavalues, taken in the key order of
/// `score_weights`. A metavalue that is missing is skipped (the source logs a
/// warning), and a transformed value that is not positive, NaN or infinite is
/// left out of the product. An empty weight map scores `1.0`.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a present metavalue is not numeric.
pub fn compute_score(
    feature: &Feature,
    score_weights: &BTreeMap<String, LambdaScore>,
) -> Result<f64> {
    let mut score_1 = 1.0;
    for (metavalue_name, &lambda_score) in score_weights {
        if !feature.metadata.contains_key(metavalue_name) {
            continue;
        }
        let value = weight_score(
            meta_double(&feature.metadata, metavalue_name)?,
            lambda_score,
        );
        if value > 0.0 && !value.is_nan() && !value.is_infinite() {
            score_1 *= value;
        }
    }
    Ok(score_1)
}

/// Add a variable to a program (source `addVariable_`); returns its index.
///
/// `bounded` gives `[0, 1]` with both bounds, otherwise the column is
/// unbounded (the `0, 1` bounds the source still passes are not applied).
pub fn add_variable(
    problem: &mut LinearProgram,
    name: &str,
    bounded: bool,
    obj: f64,
    variable_type: VariableType,
) -> usize {
    problem.columns.push(Column {
        name: name.to_owned(),
        bound_type: if bounded {
            BoundType::DoubleBounded
        } else {
            BoundType::Unbounded
        },
        lower: 0.0,
        upper: 1.0,
        column_type: match variable_type {
            VariableType::Integer => ColumnType::Integer,
            VariableType::Continuous => ColumnType::Continuous,
        },
        objective: obj,
    });
    problem.columns.len() - 1
}

/// Add a constraint to a program (source `addConstraint_`).
pub fn add_constraint(
    problem: &mut LinearProgram,
    indices: &[usize],
    values: &[f64],
    name: &str,
    lb: f64,
    ub: f64,
    bound_type: BoundType,
) {
    problem.rows.push(Row {
        name: name.to_owned(),
        indices: indices.to_vec(),
        values: values.to_vec(),
        lower: lb,
        upper: ub,
        bound_type,
    });
}

fn check_solution(solution: &LpSolution) -> Result<()> {
    match solution.status {
        SolverStatus::Optimal | SolverStatus::Feasible => Ok(()),
        status => Err(Error::InvalidValue(format!(
            "The LP solver did not find a feasible feature selection (solver status: {status:?})."
        ))),
    }
}

fn named_rows<'a>(
    feature_name_map: &'a BTreeMap<String, Vec<Feature>>,
    name: &str,
) -> Result<&'a Vec<Feature>> {
    feature_name_map.get(name).ok_or_else(|| {
        Error::InvalidValue(format!(
            "no features are filed under '{name}' (the source's std::map::at throws)"
        ))
    })
}

/// Selection by a score-weighted linear program (source
/// `MRMFeatureSelectorScore`): one variable per candidate feature, scored by
/// [`compute_score`], and a row per name forcing exactly one candidate.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MRMFeatureSelectorScore<S> {
    /// The solver of the segment programs.
    pub solver: S,
}

impl<S: LpSolver> MRMFeatureSelectorScore<S> {
    /// A selector solving with `solver`.
    pub fn new(solver: S) -> Self {
        Self { solver }
    }

    /// The program [`optimize`](MRMFeatureSelector::optimize) solves for one
    /// segment, with its variable names in column order.
    ///
    /// # Errors
    ///
    /// As [`optimize`](MRMFeatureSelector::optimize), without the solver's.
    pub fn build_problem(
        &self,
        time_to_name: &[(f64, String)],
        feature_name_map: &BTreeMap<String, Vec<Feature>>,
        parameters: &SelectorParameters,
    ) -> Result<LinearProgram> {
        let mut variables: BTreeSet<String> = BTreeSet::new();
        let mut problem = LinearProgram::new();
        for (_, name) in time_to_name {
            let mut constraints = Vec::new();
            for feature in named_rows(feature_name_map, name)? {
                let name1 = format!("{name}_{}", feature.unique_id);
                if !variables.contains(&name1) {
                    let score = compute_score(feature, &parameters.score_weights)?;
                    let col_idx =
                        add_variable(&mut problem, &name1, true, score, parameters.variable_type);
                    constraints.push(col_idx);
                    variables.insert(name1);
                }
            }
            let constraints_values = vec![1.0; constraints.len()];
            add_constraint(
                &mut problem,
                &constraints,
                &constraints_values,
                &format!("{name}_constraint"),
                1.0,
                1.0,
                BoundType::Fixed,
            );
        }
        Ok(problem)
    }
}

impl<S: LpSolver> MRMFeatureSelector for MRMFeatureSelectorScore<S> {
    /// Build the program, solve it and return the columns whose value is at
    /// least `optimal_threshold`.
    ///
    /// The retention times, `nn_threshold` and `locality_weight` are not
    /// consulted.
    fn optimize(
        &self,
        time_to_name: &[(f64, String)],
        feature_name_map: &BTreeMap<String, Vec<Feature>>,
        parameters: &SelectorParameters,
    ) -> Result<Vec<String>> {
        let problem = self.build_problem(time_to_name, feature_name_map, parameters)?;
        let solution = self.solver.solve(&problem)?;
        check_solution(&solution)?;
        let mut result = Vec::new();
        for (c, column) in problem.columns.iter().enumerate() {
            let value = solution.column_values.get(c).copied().ok_or_else(|| {
                Error::InvalidValue("the solver returned fewer values than columns".to_owned())
            })?;
            if value >= parameters.optimal_threshold {
                result.push(column.name.clone());
            }
        }
        Ok(result)
    }
}

/// Selection by a quadratic programming formulation over relative retention
/// time (source `MRMFeatureSelectorQMIP`).
///
/// For each name it considers the neighbours within `nn_threshold` and adds,
/// per pair of candidates, a continuous variable with its absolute value
/// linearisation and five constraints scoring the deviation from the expected
/// retention time difference. Scores are raised to `1 / n` for `n > 1` score
/// weights. It needs a general mixed integer solver; none ships with this
/// crate.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MRMFeatureSelectorQMIP<S> {
    /// The solver of the segment programs.
    pub solver: S,
}

impl<S: LpSolver> MRMFeatureSelectorQMIP<S> {
    /// A selector solving with `solver`.
    pub fn new(solver: S) -> Self {
        Self { solver }
    }

    /// The program [`optimize`](MRMFeatureSelector::optimize) solves for one
    /// segment, and the set of names of the candidate variables (the other
    /// columns are the pair and absolute value variables).
    ///
    /// # Errors
    ///
    /// As [`optimize`](MRMFeatureSelector::optimize), without the solver's.
    pub fn build_problem(
        &self,
        time_to_name: &[(f64, String)],
        feature_name_map: &BTreeMap<String, Vec<Feature>>,
        parameters: &SelectorParameters,
    ) -> Result<(LinearProgram, BTreeSet<String>)> {
        let mut variables: BTreeSet<String> = BTreeSet::new();
        let mut problem = LinearProgram::new();
        let n_names = time_to_name.len();
        let nn = i64::from(parameters.nn_threshold);
        for cnt1 in 0..n_names {
            let c1 = cnt1 as i64;
            // Size start_iter = max(cnt1 - nn, 0);
            let start_iter = (c1 - nn).max(0) as u64;
            // Size stop_iter = min(Size(cnt1 + nn + 1), size) -- a negative
            // sum converts to a huge `Size`.
            let stop_iter = ((c1 + nn + 1) as u64).min(n_names as u64);
            let mut constraints: Vec<usize> = Vec::new();
            let feature_row1 = named_rows(feature_name_map, &time_to_name[cnt1].1)?;
            for (i, feature1) in feature_row1.iter().enumerate() {
                let name1 = format!("{}_{}", time_to_name[cnt1].1, feature1.unique_id);
                if !variables.contains(&name1) {
                    constraints.push(add_variable(
                        &mut problem,
                        &name1,
                        true,
                        0.0,
                        parameters.variable_type,
                    ));
                    variables.insert(name1.clone());
                } else {
                    constraints.push(first_index(&problem, &name1)?);
                }
                let mut score_1 = compute_score(feature1, &parameters.score_weights)?;
                let n_score_weights = parameters.score_weights.len();
                if n_score_weights > 1 {
                    score_1 = score_1.powf(1.0 / n_score_weights as f64);
                }
                let index1 = first_index(&problem, &name1)?;

                let mut cnt2 = start_iter;
                while cnt2 < stop_iter {
                    let c2 = cnt2 as usize;
                    if cnt1 == c2 {
                        cnt2 += 1;
                        continue;
                    }
                    let feature_row2 = named_rows(feature_name_map, &time_to_name[c2].1)?;
                    let locality_weight = if parameters.locality_weight {
                        let distance = ((start_iter + cnt2) as i64 - c1).abs();
                        1.0 / ((nn - distance + 1) as f64)
                    } else {
                        1.0
                    };
                    let tr_delta_expected = time_to_name[cnt1].0 - time_to_name[c2].0;
                    for (j, feature2) in feature_row2.iter().enumerate() {
                        let name2 = format!("{}_{}", time_to_name[c2].1, feature2.unique_id);
                        if !variables.contains(&name2) {
                            add_variable(&mut problem, &name2, true, 0.0, parameters.variable_type);
                            variables.insert(name2.clone());
                        }
                        let var_qp_name = format!(
                            "{}_{}-{}_{}",
                            time_to_name[cnt1].1, i, time_to_name[c2].1, j
                        );
                        let index_var_qp = add_variable(
                            &mut problem,
                            &var_qp_name,
                            true,
                            0.0,
                            VariableType::Continuous,
                        );
                        let index_var_abs = add_variable(
                            &mut problem,
                            &format!("{var_qp_name}-ABS"),
                            false,
                            1.0,
                            VariableType::Continuous,
                        );
                        let index2 = first_index(&problem, &name2)?;
                        let mut score_2 = compute_score(feature2, &parameters.score_weights)?;
                        if n_score_weights > 1 {
                            score_2 = score_2.powf(1.0 / n_score_weights as f64);
                        }
                        let tr_delta = feature1.rt - feature2.rt;
                        let score =
                            locality_weight * score_1 * score_2 * (tr_delta - tr_delta_expected);
                        add_constraint(
                            &mut problem,
                            &[index1, index_var_qp],
                            &[1.0, -1.0],
                            &format!("{var_qp_name}-QP1"),
                            0.0,
                            1.0,
                            BoundType::LowerBoundOnly,
                        );
                        add_constraint(
                            &mut problem,
                            &[index2, index_var_qp],
                            &[1.0, -1.0],
                            &format!("{var_qp_name}-QP2"),
                            0.0,
                            1.0,
                            BoundType::LowerBoundOnly,
                        );
                        add_constraint(
                            &mut problem,
                            &[index1, index2, index_var_qp],
                            &[1.0, 1.0, -1.0],
                            &format!("{var_qp_name}-QP3"),
                            0.0,
                            1.0,
                            BoundType::UpperBoundOnly,
                        );
                        let indices_abs = [index_var_abs, index_var_qp];
                        add_constraint(
                            &mut problem,
                            &indices_abs,
                            &[-1.0, score],
                            &format!("{var_qp_name}-obj+"),
                            -1.0,
                            0.0,
                            BoundType::UpperBoundOnly,
                        );
                        add_constraint(
                            &mut problem,
                            &indices_abs,
                            &[-1.0, -score],
                            &format!("{var_qp_name}-obj-"),
                            -1.0,
                            0.0,
                            BoundType::UpperBoundOnly,
                        );
                    }
                    cnt2 += 1;
                }
            }
            let constraints_values = vec![1.0; constraints.len()];
            add_constraint(
                &mut problem,
                &constraints,
                &constraints_values,
                &format!("{}_constraint", time_to_name[cnt1].1),
                1.0,
                1.0,
                BoundType::Fixed,
            );
        }
        Ok((problem, variables))
    }
}

fn first_index(problem: &LinearProgram, name: &str) -> Result<usize> {
    problem
        .column_index(name)
        .ok_or_else(|| Error::InvalidValue(format!("column '{name}' was not created")))
}

impl<S: LpSolver> MRMFeatureSelector for MRMFeatureSelectorQMIP<S> {
    /// Build the program, solve it and return the candidate columns whose
    /// value is **greater than** `optimal_threshold` (the Score selector uses
    /// "at least").
    fn optimize(
        &self,
        time_to_name: &[(f64, String)],
        feature_name_map: &BTreeMap<String, Vec<Feature>>,
        parameters: &SelectorParameters,
    ) -> Result<Vec<String>> {
        let (problem, variables) =
            self.build_problem(time_to_name, feature_name_map, parameters)?;
        let solution = self.solver.solve(&problem)?;
        check_solution(&solution)?;
        let mut result = Vec::new();
        for (c, column) in problem.columns.iter().enumerate() {
            let value = solution.column_values.get(c).copied().ok_or_else(|| {
                Error::InvalidValue("the solver returned fewer values than columns".to_owned())
            })?;
            if value > parameters.optimal_threshold && variables.contains(&column.name) {
                result.push(column.name.clone());
            }
        }
        Ok(result)
    }
}
