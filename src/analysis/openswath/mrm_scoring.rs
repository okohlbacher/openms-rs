// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Scores for peaks found in SRM/MRM, from `ANALYSIS/OPENSWATH/MRMScoring.h`
//! (class `OpenSwath::MRMScoring`).
//!
//! The scores use parameters of the peaks from the individual transitions,
//! following Reiter L, Rinner O, Picotti P, Huettenhain R, Beck M, Brusniak MY,
//! Hengartner MO, Aebersold R. *mProphet: automated data processing and
//! statistical validation for large-scale SRM experiments.* Nat Methods. 2011
//! May;8(5):430-5. Implemented:
//!
//! - *xcorr_coelution*: cross-correlation of the transitions;
//! - *xcorr_shape*: whether the maximal cross-correlation coincides with the
//!   maximal intensity;
//! - *library_rmsd* and *library_correlation*: normalised RMSD and correlation
//!   of the measured to the expected intensities;
//! - *rt_score*: deviation from the expected retention time;
//! - signal-to-noise and ranked mutual-information scores.
//!
//! [`MRMScoring`](crate::analysis::openswath::mrm_scoring::MRMScoring) holds the precomputed cross-correlation and
//! mutual-information matrices; the `initialize_*` methods build them and the
//! `calc_*` methods reduce them to scores. Cross-correlations are computed for
//! lags up to 10 points (the source's `XCORR_MAX_DELAY`), or the trace length if
//! shorter.
//!
//! # Native differences
//!
//! - The matrices are [`ScoreMatrix`](crate::analysis::openswath::mrm_scoring::ScoreMatrix) values, a column-major stand-in for the
//!   source's `OpenMS::Matrix` (from `DATASTRUCTURES/Matrix.h`, which this port
//!   does not provide yet). Every initialisation starts from a zeroed matrix;
//!   the source resizes in place, which leaves earlier values in entries it does
//!   not write (only lower triangles no score reads).
//! - Whole-matrix sums add in column-major order, left to right; the source's
//!   Eigen `sum()` adds in vectorised blocks, so results agree to rounding.
//! - The source's `OPENSWATH_PRECONDITION`s on matrix shape are debug-only
//!   `assert`s; the release build computes on regardless, and so does the port
//!   (an empty matrix reduces to NaN or 0 as in the source's arithmetic).
//! - The source keeps per-thread scratch pools; the port allocates per call.
//!
//! See `docs/OPENSWATH_SCORING_SUPPORT.md`.

use super::itransition::{MRMFeatureAccess, SignalToNoisePtr};
use super::scoring::{
    XCorrArrayType, compute_rank_vector, normalize_sum, normalized_cross_correlation_post,
    normalized_manhattan_dist, ranked_mutual_information, root_mean_square_deviation,
    spectral_angle, standardize_data, xcorr_array_get_max_peak,
};
use super::stats_helpers::{MeanAndStddev, dotprod_scoring, manhattan_scoring};
use super::transition_experiment::{LightCompound, LightTransition};
use crate::math::statistic_functions::pearson_correlation_coefficient;
use crate::{Error, Result};

/// Maximum lag, in data points, of the cross-correlations (source
/// `XCORR_MAX_DELAY`); lags beyond it are physically meaningless for coelution.
pub const XCORR_MAX_DELAY: i32 = 10;

/// A dense matrix stored column by column, standing in for the source's
/// `OpenMS::Matrix` in the scoring matrices of [`MRMScoring`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ScoreMatrix<T> {
    rows: usize,
    cols: usize,
    data: Vec<T>,
}

impl<T: Clone + Default> ScoreMatrix<T> {
    /// A `rows x cols` matrix of default values (zero for numbers).
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![T::default(); rows * cols],
        }
    }
}

impl<T> ScoreMatrix<T> {
    /// Number of rows.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Number of entries, `rows * cols`.
    pub fn size(&self) -> usize {
        self.data.len()
    }

    /// The entry at `(row, col)`, or `None` outside the matrix.
    pub fn get(&self, row: usize, col: usize) -> Option<&T> {
        (row < self.rows && col < self.cols).then(|| &self.data[col * self.rows + row])
    }

    /// The entries in column-major order, as the source's `data()`.
    pub fn data(&self) -> &[T] {
        &self.data
    }

    fn set(&mut self, row: usize, col: usize, value: T) {
        self.data[col * self.rows + row] = value;
    }

    fn at(&self, row: usize, col: usize) -> &T {
        &self.data[col * self.rows + row]
    }
}

/// A cross-correlation matrix (source `XCorrMatrixType`).
pub type XCorrMatrixType = ScoreMatrix<XCorrArrayType>;

/// One cross-correlation matrix with its per-entry apex lag and apex value.
#[derive(Clone, Debug, Default, PartialEq)]
struct XCorrSet {
    xcorr: XCorrMatrixType,
    max_peak: ScoreMatrix<i32>,
    max_peak_sec: ScoreMatrix<f64>,
}

impl XCorrSet {
    fn new(rows: usize, cols: usize) -> Self {
        Self {
            xcorr: ScoreMatrix::new(rows, cols),
            max_peak: ScoreMatrix::new(rows, cols),
            max_peak_sec: ScoreMatrix::new(rows, cols),
        }
    }

    fn fill(&mut self, i: usize, j: usize, a: &[f64], b: &[f64]) -> Result<()> {
        let delay = XCORR_MAX_DELAY.min(i32::try_from(a.len()).unwrap_or(i32::MAX));
        let array = normalized_cross_correlation_post(a, b, delay, 1)?;
        let apex = array.data[xcorr_array_get_max_peak(&array)?];
        self.max_peak.set(i, j, apex.0.abs());
        self.max_peak_sec.set(i, j, apex.1);
        self.xcorr.set(i, j, array);
        Ok(())
    }
}

/// SRM/MRM peak-group scores (source `OpenSwath::MRMScoring`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MRMScoring {
    xcorr: XCorrSet,
    xcorr_contrast: XCorrSet,
    xcorr_precursor: XCorrSet,
    xcorr_precursor_contrast: XCorrSet,
    xcorr_precursor_combined: XCorrSet,
    mi: ScoreMatrix<f64>,
    mi_contrast: ScoreMatrix<f64>,
    mi_precursor: ScoreMatrix<f64>,
    mi_precursor_contrast: ScoreMatrix<f64>,
    mi_precursor_combined: ScoreMatrix<f64>,
}

fn missing(kind: &str, id: &str) -> Error {
    Error::MissingInformation(format!("no {kind} feature '{id}'"))
}

/// The intensity traces of the transition features `ids` (source
/// `fillIntensityFromFeature`).
fn feature_intensities(feature: &dyn MRMFeatureAccess, ids: &[String]) -> Result<Vec<Vec<f64>>> {
    ids.iter()
        .map(|id| {
            feature
                .feature(id)
                .ok_or_else(|| missing("transition", id))?
                .intensity_values()
        })
        .collect()
}

/// The intensity traces of the precursor features `ids` (source
/// `fillIntensityFromPrecursorFeature`).
fn precursor_intensities(feature: &dyn MRMFeatureAccess, ids: &[String]) -> Result<Vec<Vec<f64>>> {
    ids.iter()
        .map(|id| {
            feature
                .precursor_feature(id)
                .ok_or_else(|| missing("precursor", id))?
                .intensity_values()
        })
        .collect()
}

fn standardized(mut data: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    for trace in &mut data {
        standardize_data(trace);
    }
    data
}

/// Upper-triangle (with diagonal) cross-correlation of `data` against itself.
fn triangle(data: &[Vec<f64>]) -> Result<XCorrSet> {
    let mut set = XCorrSet::new(data.len(), data.len());
    for i in 0..data.len() {
        for j in i..data.len() {
            set.fill(i, j, &data[i], &data[j])?;
        }
    }
    Ok(set)
}

/// Full cross-correlation of every row of `a` against every row of `b`.
fn rectangle(a: &[Vec<f64>], b: &[Vec<f64>]) -> Result<XCorrSet> {
    let mut set = XCorrSet::new(a.len(), b.len());
    for (i, row) in a.iter().enumerate() {
        for (j, col) in b.iter().enumerate() {
            set.fill(i, j, row, col)?;
        }
    }
    Ok(set)
}

/// `w^T M_sym w` for an upper-triangular `M` (source `weightedTriangularSum`).
fn weighted_triangular_sum<T: Copy + Into<f64>>(mat: &ScoreMatrix<T>, w: &[f64]) -> Result<f64> {
    if mat.rows != mat.cols || mat.rows != w.len() {
        return Err(Error::InvalidValue(
            "matrix must be square and match the weight vector size".into(),
        ));
    }
    let mut sum = 0.0;
    for i in 0..w.len() {
        let wi = w[i];
        sum += wi * wi * (*mat.at(i, i)).into();
        for (j, &wj) in w.iter().enumerate().skip(i + 1) {
            sum += 2.0 * wi * wj * (*mat.at(i, j)).into();
        }
    }
    Ok(sum)
}

/// Mean plus sample standard deviation of the upper triangle with diagonal
/// (source `meanStdevUpperTriangle`).
fn mean_stdev_upper_triangle(mat: &ScoreMatrix<i32>) -> f64 {
    let mut msc = MeanAndStddev::new();
    for i in 0..mat.rows {
        for j in i..mat.cols {
            msc.push(f64::from(*mat.at(i, j)));
        }
    }
    msc.mean() + msc.sample_stddev()
}

/// Mean plus sample standard deviation of all entries, in storage order
/// (source `meanStdevAll`).
fn mean_stdev_all(mat: &ScoreMatrix<i32>) -> f64 {
    let msc: MeanAndStddev = mat.data.iter().map(|&v| f64::from(v)).collect();
    msc.mean() + msc.sample_stddev()
}

/// Mean of the upper triangle with diagonal; 0 for an empty matrix (source
/// `meanUpperTriangle`).
fn mean_upper_triangle(mat: &ScoreMatrix<f64>) -> f64 {
    let mut sum = 0.0;
    let mut n = 0usize;
    for i in 0..mat.rows {
        for j in i..mat.cols {
            sum += *mat.at(i, j);
            n += 1;
        }
    }
    if n > 0 { sum / n as f64 } else { 0.0 }
}

fn matrix_sum(mat: &ScoreMatrix<f64>) -> f64 {
    mat.data.iter().fold(0.0, |acc, v| acc + v)
}

fn row_means(mat: &ScoreMatrix<f64>) -> Vec<f64> {
    (0..mat.rows)
        .map(|i| {
            let mut sum = 0.0;
            for j in 0..mat.cols {
                sum += *mat.at(i, j);
            }
            sum / mat.cols as f64
        })
        .collect()
}

/// Number of entries in the upper triangle with diagonal of an `n x n` matrix,
/// in the source's integer arithmetic `n*n/2 + (n+1)/2`.
fn triangular_count(n: usize) -> f64 {
    (n * n / 2 + n.div_ceil(2)) as f64
}

fn ranks(data: &[Vec<f64>]) -> Result<(Vec<Vec<u32>>, Vec<u32>)> {
    let mut rank_vec = Vec::new();
    let max = compute_rank_vector(data, &mut rank_vec)?;
    Ok((rank_vec, max))
}

impl MRMScoring {
    /// An empty scorer: every matrix is `0 x 0`.
    pub fn new() -> Self {
        Self::default()
    }

    /// The cross-correlation matrix (source `getXCorrMatrix`).
    pub fn xcorr_matrix(&self) -> &XCorrMatrixType {
        &self.xcorr.xcorr
    }

    /// The contrast cross-correlation matrix (source `getXCorrContrastMatrix`).
    pub fn xcorr_contrast_matrix(&self) -> &XCorrMatrixType {
        &self.xcorr_contrast.xcorr
    }

    /// The precursor cross-correlation matrix. The source has the member but
    /// no accessor; this one is native.
    pub fn xcorr_precursor_matrix(&self) -> &XCorrMatrixType {
        &self.xcorr_precursor.xcorr
    }

    /// The precursor-contrast cross-correlation matrix (source
    /// `getXCorrPrecursorContrastMatrix`).
    pub fn xcorr_precursor_contrast_matrix(&self) -> &XCorrMatrixType {
        &self.xcorr_precursor_contrast.xcorr
    }

    /// The precursor-combined cross-correlation matrix (source
    /// `getXCorrPrecursorCombinedMatrix`, documented there as the "precursor
    /// contrast matrix", a copy-paste slip).
    pub fn xcorr_precursor_combined_matrix(&self) -> &XCorrMatrixType {
        &self.xcorr_precursor_combined.xcorr
    }

    /// The mutual-information matrix (source `getMIMatrix`); only the upper
    /// triangle with diagonal is written, the rest is 0.
    pub fn mi_matrix(&self) -> &ScoreMatrix<f64> {
        &self.mi
    }

    /// The contrast mutual-information matrix (source `getMIContrastMatrix`).
    pub fn mi_contrast_matrix(&self) -> &ScoreMatrix<f64> {
        &self.mi_contrast
    }

    /// The precursor mutual-information matrix. The source has the member but
    /// no accessor; this one is native.
    pub fn mi_precursor_matrix(&self) -> &ScoreMatrix<f64> {
        &self.mi_precursor
    }

    /// The precursor-contrast mutual-information matrix (source
    /// `getMIPrecursorContrastMatrix`).
    pub fn mi_precursor_contrast_matrix(&self) -> &ScoreMatrix<f64> {
        &self.mi_precursor_contrast
    }

    /// The precursor-combined mutual-information matrix, filled symmetrically
    /// (source `getMIPrecursorCombinedMatrix`).
    pub fn mi_precursor_combined_matrix(&self) -> &ScoreMatrix<f64> {
        &self.mi_precursor_combined
    }

    /// Build the cross-correlation matrix from raw traces (source
    /// `initializeXCorrMatrix(const std::vector<std::vector<double>>&)`).
    ///
    /// Each trace is standardised (a copy; `data` is unchanged) and every pair
    /// `i <= j` cross-correlated; the lower triangle stays empty.
    ///
    /// # Errors
    ///
    /// The cross-correlation errors: traces of different or zero length.
    pub fn initialize_xcorr_matrix_from_data(&mut self, data: &[Vec<f64>]) -> Result<()> {
        self.xcorr = triangle(&standardized(data.to_vec()))?;
        Ok(())
    }

    /// Build the cross-correlation matrix of the transitions `native_ids` of
    /// `feature` (source `initializeXCorrMatrix(IMRMFeature*, ...)`).
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] for an unknown native ID (the source
    /// dereferences a null pointer or throws), the feature's trace errors, and
    /// the cross-correlation errors.
    pub fn initialize_xcorr_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        native_ids: &[String],
    ) -> Result<()> {
        self.xcorr = triangle(&standardized(feature_intensities(feature, native_ids)?))?;
        Ok(())
    }

    /// Build the cross-correlation matrix of the transitions `set1` (for
    /// example identification transitions) against `set2` (for example
    /// detection transitions) (source `initializeXCorrContrastMatrix`).
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::initialize_xcorr_matrix`].
    pub fn initialize_xcorr_contrast_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        set1: &[String],
        set2: &[String],
    ) -> Result<()> {
        let a = standardized(feature_intensities(feature, set1)?);
        let b = standardized(feature_intensities(feature, set2)?);
        self.xcorr_contrast = rectangle(&a, &b)?;
        Ok(())
    }

    /// Build the cross-correlation matrix of the precursor traces
    /// `precursor_ids` (source `initializeXCorrPrecursorMatrix`).
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::initialize_xcorr_matrix`], for precursor features.
    pub fn initialize_xcorr_precursor_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        precursor_ids: &[String],
    ) -> Result<()> {
        self.xcorr_precursor = triangle(&standardized(precursor_intensities(
            feature,
            precursor_ids,
        )?))?;
        Ok(())
    }

    /// Build the cross-correlation matrix of the precursor isotope traces
    /// against the transitions (source
    /// `initializeXCorrPrecursorContrastMatrix(IMRMFeature*, ...)`).
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::initialize_xcorr_matrix`].
    pub fn initialize_xcorr_precursor_contrast_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        precursor_ids: &[String],
        native_ids: &[String],
    ) -> Result<()> {
        let a = standardized(precursor_intensities(feature, precursor_ids)?);
        let b = standardized(feature_intensities(feature, native_ids)?);
        self.xcorr_precursor_contrast = rectangle(&a, &b)?;
        Ok(())
    }

    /// Build the precursor-contrast cross-correlation matrix from raw traces
    /// (source `initializeXCorrPrecursorContrastMatrix(const
    /// std::vector<std::vector<double>>&, ...)`); the inputs are unchanged.
    ///
    /// # Errors
    ///
    /// The cross-correlation errors.
    pub fn initialize_xcorr_precursor_contrast_matrix_from_data(
        &mut self,
        data_precursor: &[Vec<f64>],
        data_fragments: &[Vec<f64>],
    ) -> Result<()> {
        let a = standardized(data_precursor.to_vec());
        let b = standardized(data_fragments.to_vec());
        self.xcorr_precursor_contrast = rectangle(&a, &b)?;
        Ok(())
    }

    /// Build the cross-correlation matrix of the precursor traces followed by
    /// the transitions, as one combined upper triangle (source
    /// `initializeXCorrPrecursorCombinedMatrix`).
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::initialize_xcorr_matrix`].
    pub fn initialize_xcorr_precursor_combined_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        precursor_ids: &[String],
        native_ids: &[String],
    ) -> Result<()> {
        let mut combined = standardized(precursor_intensities(feature, precursor_ids)?);
        combined.extend(standardized(feature_intensities(feature, native_ids)?));
        self.xcorr_precursor_combined = triangle(&combined)?;
        Ok(())
    }

    /// The cross-correlation coelution score: mean plus sample standard
    /// deviation of the apex lags over the upper triangle. A distance; 0 means
    /// perfect coelution (source `calcXcorrCoelutionScore`).
    pub fn calc_xcorr_coelution_score(&self) -> f64 {
        mean_stdev_upper_triangle(&self.xcorr.max_peak)
    }

    /// The weighted cross-correlation coelution score, `w^T L w` over the
    /// symmetric apex-lag matrix `L`; non-perfect coelution of low-intensity
    /// transitions weighs less (source `calcXcorrCoelutionWeightedScore`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the weights do not match the matrix size
    /// (a debug-only precondition in the source, which otherwise reads out of
    /// bounds).
    pub fn calc_xcorr_coelution_weighted_score(
        &self,
        normalized_library_intensity: &[f64],
    ) -> Result<f64> {
        weighted_triangular_sum(&self.xcorr.max_peak, normalized_library_intensity)
    }

    /// Per row of the contrast matrix, the mean apex lag (source
    /// `calcSeparateXcorrContrastCoelutionScore`).
    pub fn calc_separate_xcorr_contrast_coelution_score(&self) -> Vec<f64> {
        let m = &self.xcorr_contrast.max_peak;
        (0..m.rows)
            .map(|i| {
                let mut sum = 0.0;
                for j in 0..m.cols {
                    sum += f64::from(*m.at(i, j));
                }
                sum / m.cols as f64
            })
            .collect()
    }

    /// The precursor cross-correlation coelution score (source
    /// `calcXcorrPrecursorCoelutionScore`).
    pub fn calc_xcorr_precursor_coelution_score(&self) -> f64 {
        mean_stdev_upper_triangle(&self.xcorr_precursor.max_peak)
    }

    /// The precursor-contrast coelution score against the transitions: mean
    /// plus sample standard deviation of all apex lags; 0 means perfect
    /// coelution (source `calcXcorrPrecursorContrastCoelutionScore`).
    pub fn calc_xcorr_precursor_contrast_coelution_score(&self) -> f64 {
        mean_stdev_all(&self.xcorr_precursor_contrast.max_peak)
    }

    /// The precursor-contrast coelution score against the summed transitions;
    /// computed exactly as
    /// [`MRMScoring::calc_xcorr_precursor_contrast_coelution_score`], the
    /// source differing only in its debug assertion (source
    /// `calcXcorrPrecursorContrastSumFragCoelutionScore`).
    pub fn calc_xcorr_precursor_contrast_sum_frag_coelution_score(&self) -> f64 {
        mean_stdev_all(&self.xcorr_precursor_contrast.max_peak)
    }

    /// The combined precursor and transition coelution score (source
    /// `calcXcorrPrecursorCombinedCoelutionScore`).
    pub fn calc_xcorr_precursor_combined_coelution_score(&self) -> f64 {
        mean_stdev_upper_triangle(&self.xcorr_precursor_combined.max_peak)
    }

    /// The cross-correlation shape score: the mean apex correlation over the
    /// upper triangle; 1 is perfect correlation, 0 none (source
    /// `calcXcorrShapeScore`).
    pub fn calc_xcorr_shape_score(&self) -> f64 {
        mean_upper_triangle(&self.xcorr.max_peak_sec)
    }

    /// The weighted cross-correlation shape score, `w^T C w` over the
    /// symmetric apex-correlation matrix (source
    /// `calcXcorrShapeWeightedScore`).
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::calc_xcorr_coelution_weighted_score`].
    pub fn calc_xcorr_shape_weighted_score(
        &self,
        normalized_library_intensity: &[f64],
    ) -> Result<f64> {
        weighted_triangular_sum(&self.xcorr.max_peak_sec, normalized_library_intensity)
    }

    /// The contrast shape score: the **sum** of all apex correlations of the
    /// contrast matrix, not their mean (source `calcXcorrContrastShapeScore`).
    pub fn calc_xcorr_contrast_shape_score(&self) -> f64 {
        matrix_sum(&self.xcorr_contrast.max_peak_sec)
    }

    /// Per row of the contrast matrix, the mean apex correlation (source
    /// `calcSeparateXcorrContrastShapeScore`).
    pub fn calc_separate_xcorr_contrast_shape_score(&self) -> Vec<f64> {
        row_means(&self.xcorr_contrast.max_peak_sec)
    }

    /// The precursor shape score (source `calcXcorrPrecursorShapeScore`).
    pub fn calc_xcorr_precursor_shape_score(&self) -> f64 {
        mean_upper_triangle(&self.xcorr_precursor.max_peak_sec)
    }

    /// The precursor-contrast shape score against the summed transitions, the
    /// mean of all apex correlations (source
    /// `calcXcorrPrecursorContrastSumFragShapeScore`).
    pub fn calc_xcorr_precursor_contrast_sum_frag_shape_score(&self) -> f64 {
        let m = &self.xcorr_precursor_contrast.max_peak_sec;
        matrix_sum(m) / m.size() as f64
    }

    /// The precursor-contrast shape score against the transitions, the mean of
    /// all apex correlations (source `calcXcorrPrecursorContrastShapeScore`).
    pub fn calc_xcorr_precursor_contrast_shape_score(&self) -> f64 {
        let m = &self.xcorr_precursor_contrast.max_peak_sec;
        matrix_sum(m) / m.size() as f64
    }

    /// The combined precursor and transition shape score (source
    /// `calcXcorrPrecursorCombinedShapeScore`).
    pub fn calc_xcorr_precursor_combined_shape_score(&self) -> f64 {
        mean_upper_triangle(&self.xcorr_precursor_combined.max_peak_sec)
    }

    /// Library scores of a peak group against the transitions' library
    /// intensities (source static `calcLibraryScore`).
    ///
    /// Negative library intensities count as 0; the experimental intensity of
    /// each transition is its feature's apex intensity. Computed, in the
    /// source's order: the Manhattan and dot-product scores (square-rooted
    /// intensities), the spectral angle (NaN becomes 0), then on sum-normalised
    /// vectors the normalised Manhattan distance, the RMSD and the Pearson
    /// correlation (NaN becomes -1).
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] for a transition without feature, and
    /// [`Error::InvalidValue`] or [`Error::InvalidRange`] for no transitions,
    /// where the source's release build divides by zero.
    pub fn calc_library_score(
        feature: &dyn MRMFeatureAccess,
        transitions: &[LightTransition],
    ) -> Result<LibraryScores> {
        let mut library = Vec::with_capacity(transitions.len());
        let mut experimental = Vec::with_capacity(transitions.len());
        for transition in transitions {
            let id = transition.native_id();
            let f = feature
                .feature(id)
                .ok_or_else(|| missing("transition", id))?;
            experimental.push(f64::from(f.intensity()));
            let intensity = transition.library_intensity;
            library.push(if intensity < 0.0 { 0.0 } else { intensity });
        }
        let manhattan = manhattan_scoring(experimental.clone(), library.clone())?;
        let dotprod = dotprod_scoring(experimental.clone(), library.clone())?;
        let mut angle = spectral_angle(&experimental, &library)?;
        if angle.is_nan() {
            angle = 0.0;
        }
        normalize_sum(&mut experimental);
        normalize_sum(&mut library);
        let norm_manhattan = normalized_manhattan_dist(&mut experimental, &mut library)?;
        let rmsd = root_mean_square_deviation(&experimental, &library)?;
        let mut correlation = pearson_correlation_coefficient(&experimental, &library)?;
        if correlation.is_nan() {
            correlation = -1.0;
        }
        Ok(LibraryScores {
            correlation,
            norm_manhattan,
            manhattan,
            dotprod,
            spectral_angle: angle,
            rmsd,
        })
    }

    /// The retention-time score, `|normalized_experimental_rt - peptide.rt|`,
    /// or 0 when the expected RT is `<= -1000` (source static `calcRTScore`).
    ///
    /// A NaN expected RT (the port's "unset") compares false and yields NaN,
    /// as in the source.
    pub fn calc_rt_score(peptide: &LightCompound, normalized_experimental_rt: f64) -> f64 {
        if peptide.rt <= -1000.0 {
            return 0.0;
        }
        (normalized_experimental_rt - peptide.rt).abs()
    }

    /// The mean signal-to-noise ratio of the estimators at the feature RT; 0
    /// without estimators (source static `calcSNScore`).
    pub fn calc_sn_score(feature: &dyn MRMFeatureAccess, estimators: &[SignalToNoisePtr]) -> f64 {
        if estimators.is_empty() {
            return 0.0;
        }
        let rt = feature.rt();
        let mut sum = 0.0;
        for estimator in estimators {
            sum += estimator.value_at_rt(rt);
        }
        sum / estimators.len() as f64
    }

    /// The natural logarithm of each estimator's S/N at the feature RT, with
    /// values below 1 set to 0 so the logarithm is safe (source static
    /// `calcSeparateSNScore`).
    pub fn calc_separate_sn_score(
        feature: &dyn MRMFeatureAccess,
        estimators: &[SignalToNoisePtr],
    ) -> Vec<f64> {
        let rt = feature.rt();
        estimators
            .iter()
            .map(|estimator| {
                let sn = estimator.value_at_rt(rt);
                if sn < 1.0 { 0.0 } else { sn.ln() }
            })
            .collect()
    }

    /// Build the ranked mutual-information matrix of the transitions (source
    /// `initializeMIMatrix`); the lower triangle is 0.
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::initialize_xcorr_matrix`], and the mutual-information
    /// errors (traces of different or zero length).
    pub fn initialize_mi_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        native_ids: &[String],
    ) -> Result<()> {
        let (rank_vec, max) = ranks(&feature_intensities(feature, native_ids)?)?;
        self.mi = mi_triangle(&rank_vec, &max, false)?;
        Ok(())
    }

    /// Build the mutual-information matrix of `set1` against `set2` (source
    /// `initializeMIContrastMatrix`).
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::initialize_mi_matrix`].
    pub fn initialize_mi_contrast_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        set1: &[String],
        set2: &[String],
    ) -> Result<()> {
        let (r1, m1) = ranks(&feature_intensities(feature, set1)?)?;
        let (r2, m2) = ranks(&feature_intensities(feature, set2)?)?;
        self.mi_contrast = mi_rectangle(&r1, &m1, &r2, &m2)?;
        Ok(())
    }

    /// Build the mutual-information matrix of the precursor traces (source
    /// `initializeMIPrecursorMatrix`); the lower triangle is 0.
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::initialize_mi_matrix`], for precursor features.
    pub fn initialize_mi_precursor_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        precursor_ids: &[String],
    ) -> Result<()> {
        let (rank_vec, max) = ranks(&precursor_intensities(feature, precursor_ids)?)?;
        self.mi_precursor = mi_triangle(&rank_vec, &max, false)?;
        Ok(())
    }

    /// Build the mutual-information matrix of the precursor traces against the
    /// transitions (source `initializeMIPrecursorContrastMatrix`).
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::initialize_mi_matrix`].
    pub fn initialize_mi_precursor_contrast_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        precursor_ids: &[String],
        native_ids: &[String],
    ) -> Result<()> {
        let (r1, m1) = ranks(&precursor_intensities(feature, precursor_ids)?)?;
        let (r2, m2) = ranks(&feature_intensities(feature, native_ids)?)?;
        self.mi_precursor_contrast = mi_rectangle(&r1, &m1, &r2, &m2)?;
        Ok(())
    }

    /// Build the symmetric mutual-information matrix of the precursor traces
    /// followed by the transitions (source
    /// `initializeMIPrecursorCombinedMatrix`).
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::initialize_mi_matrix`].
    pub fn initialize_mi_precursor_combined_matrix(
        &mut self,
        feature: &dyn MRMFeatureAccess,
        precursor_ids: &[String],
        native_ids: &[String],
    ) -> Result<()> {
        let mut rank_vec = Vec::new();
        let mut max = compute_rank_vector(
            &precursor_intensities(feature, precursor_ids)?,
            &mut rank_vec,
        )?;
        max.extend(compute_rank_vector(
            &feature_intensities(feature, native_ids)?,
            &mut rank_vec,
        )?);
        self.mi_precursor_combined = mi_triangle(&rank_vec, &max, true)?;
        Ok(())
    }

    /// The mean mutual information over the upper triangle with diagonal
    /// (source `calcMIScore`). The whole matrix is summed (its lower triangle
    /// is 0) and divided by `n*n/2 + (n+1)/2` in integer arithmetic; an empty
    /// matrix gives NaN.
    pub fn calc_mi_score(&self) -> f64 {
        matrix_sum(&self.mi) / triangular_count(self.mi.rows)
    }

    /// The weighted mutual-information score, `w^T M w` over the symmetric MI
    /// matrix (source `calcMIWeightedScore`).
    ///
    /// # Errors
    ///
    /// As [`MRMScoring::calc_xcorr_coelution_weighted_score`].
    pub fn calc_mi_weighted_score(&self, normalized_library_intensity: &[f64]) -> Result<f64> {
        weighted_triangular_sum(&self.mi, normalized_library_intensity)
    }

    /// The mean precursor mutual information (source `calcMIPrecursorScore`),
    /// computed as [`MRMScoring::calc_mi_score`].
    pub fn calc_mi_precursor_score(&self) -> f64 {
        matrix_sum(&self.mi_precursor) / triangular_count(self.mi_precursor.rows)
    }

    /// The mean of all precursor-contrast mutual informations (source
    /// `calcMIPrecursorContrastScore`).
    pub fn calc_mi_precursor_contrast_score(&self) -> f64 {
        matrix_sum(&self.mi_precursor_contrast) / self.mi_precursor_contrast.size() as f64
    }

    /// The mean of all entries of the symmetric combined matrix (source
    /// `calcMIPrecursorCombinedScore`).
    pub fn calc_mi_precursor_combined_score(&self) -> f64 {
        matrix_sum(&self.mi_precursor_combined) / self.mi_precursor_combined.size() as f64
    }

    /// Per row of the contrast MI matrix, the mean (source
    /// `calcSeparateMIContrastScore`).
    pub fn calc_separate_mi_contrast_score(&self) -> Vec<f64> {
        row_means(&self.mi_contrast)
    }
}

fn mi_triangle(rank_vec: &[Vec<u32>], max: &[u32], symmetric: bool) -> Result<ScoreMatrix<f64>> {
    let n = rank_vec.len();
    let mut m = ScoreMatrix::new(n, n);
    for i in 0..n {
        for j in i..n {
            let value = ranked_mutual_information(&rank_vec[i], &rank_vec[j], max[i], max[j])?;
            m.set(i, j, value);
            if symmetric && i != j {
                m.set(j, i, value);
            }
        }
    }
    Ok(m)
}

fn mi_rectangle(
    r1: &[Vec<u32>],
    m1: &[u32],
    r2: &[Vec<u32>],
    m2: &[u32],
) -> Result<ScoreMatrix<f64>> {
    let mut m = ScoreMatrix::new(r1.len(), r2.len());
    for i in 0..r1.len() {
        for j in 0..r2.len() {
            m.set(
                i,
                j,
                ranked_mutual_information(&r1[i], &r2[j], m1[i], m2[j])?,
            );
        }
    }
    Ok(m)
}

/// The library scores of [`MRMScoring::calc_library_score`], named as the
/// source's out-parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LibraryScores {
    /// Pearson correlation of the sum-normalised intensities; -1 when
    /// undefined.
    pub correlation: f64,
    /// Normalised Manhattan distance of the sum-normalised intensities (the
    /// `library_rmsd` score of OpenSWATH).
    pub norm_manhattan: f64,
    /// Manhattan score of the square-rooted, sum-normalised intensities.
    pub manhattan: f64,
    /// Dot-product score of the square-rooted, unit-normalised intensities.
    pub dotprod: f64,
    /// Spectral angle in radians; 0 when undefined.
    pub spectral_angle: f64,
    /// Root mean square deviation of the sum-normalised intensities.
    pub rmsd: f64,
}
