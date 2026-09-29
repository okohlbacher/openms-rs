// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The analysis half of `ANALYSIS/OPENSWATH/OpenSwathScores.h`: the
//! `OpenSwath_Scores_Usage` switchboard and the four LDA scoring methods of
//! `OpenSwath_Scores`.
//!
//! The score records themselves (`OpenSwath_Scores` and `OpenSwath_Ind_Scores`)
//! are ported in [`crate::kernel::mrm`], because `MRMFeature` stores one by
//! value; this module adds the methods from `OpenSwathScores.cpp` to
//! [`OpenSwathScores`] and ports [`OpenSwathScoresUsage`].
//!
//! See `docs/OPENSWATH_SCORING_SUPPORT.md`.

use crate::kernel::mrm::OpenSwathScores;

/// Which scores the OpenSWATH algorithm should compute; turns individual
/// scores on or off (source `OpenSwath_Scores_Usage`).
///
/// Every switch defaults to `true` except
/// [`use_peak_shape_metrics`](Self::use_peak_shape_metrics). The source names
/// most members with a trailing underscore (`use_coelution_score_`) and three
/// without (`use_im_scores`, `use_ms1_correlation`, ...); the port drops the
/// underscores.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct OpenSwathScoresUsage {
    /// Cross-correlation coelution score (`use_coelution_score_`).
    pub use_coelution_score: bool,
    /// Cross-correlation shape score (`use_shape_score_`).
    pub use_shape_score: bool,
    /// Retention-time score (`use_rt_score_`).
    pub use_rt_score: bool,
    /// Library scores (`use_library_score_`).
    pub use_library_score: bool,
    /// Elution-model fit score (`use_elution_model_score_`).
    pub use_elution_model_score: bool,
    /// Intensity score (`use_intensity_score_`).
    pub use_intensity_score: bool,
    /// Total XIC score (`use_total_xic_score_`).
    pub use_total_xic_score: bool,
    /// Total mutual-information score (`use_total_mi_score_`).
    pub use_total_mi_score: bool,
    /// Number-of-peaks score (`use_nr_peaks_score_`).
    pub use_nr_peaks_score: bool,
    /// Signal-to-noise score (`use_sn_score_`).
    pub use_sn_score: bool,
    /// Mutual-information score (`use_mi_score_`).
    pub use_mi_score: bool,
    /// Full-spectrum DIA scores (`use_dia_scores_`).
    pub use_dia_scores: bool,
    /// Ion-mobility scores (`use_im_scores`).
    pub use_im_scores: bool,
    /// MS1 cross-correlation scores (`use_ms1_correlation`).
    pub use_ms1_correlation: bool,
    /// MS1 full-scan scores (`use_ms1_fullscan`).
    pub use_ms1_fullscan: bool,
    /// MS1 mutual-information scores (`use_ms1_mi`).
    pub use_ms1_mi: bool,
    /// Unique-ion-signature (IPF) scores (`use_uis_scores`).
    pub use_uis_scores: bool,
    /// Ion-series scores (`use_ionseries_scores`).
    pub use_ionseries_scores: bool,
    /// MS2 isotope scores (`use_ms2_isotope_scores`).
    pub use_ms2_isotope_scores: bool,
    /// Per-transition peak-shape metrics (`use_peak_shape_metrics`); off by
    /// default.
    pub use_peak_shape_metrics: bool,
}

impl Default for OpenSwathScoresUsage {
    fn default() -> Self {
        Self {
            use_coelution_score: true,
            use_shape_score: true,
            use_rt_score: true,
            use_library_score: true,
            use_elution_model_score: true,
            use_intensity_score: true,
            use_total_xic_score: true,
            use_total_mi_score: true,
            use_nr_peaks_score: true,
            use_sn_score: true,
            use_mi_score: true,
            use_dia_scores: true,
            use_im_scores: true,
            use_ms1_correlation: true,
            use_ms1_fullscan: true,
            use_ms1_mi: true,
            use_uis_scores: true,
            use_ionseries_scores: true,
            use_ms2_isotope_scores: true,
            use_peak_shape_metrics: false,
        }
    }
}

impl OpenSwathScores {
    /// A quick LDA score from six non-DIA scores (source
    /// `get_quick_lda_score`); better when **more negative**.
    ///
    /// The weights come from manual evaluation of 80 chromatograms, an average
    /// model over 100 twofold cross-validated runs (0.85 TPR / 0.17 FDR; true
    /// mean 4.2, false mean -0.075, both with SD 1.055). The source notes that
    /// a cut-off of -0.5 removes about 30% of the peaks, 0 about 50%, 0.5 about
    /// 70%, 1.0 about 85%, 1.5 about 93% and 2.0 about 97%.
    ///
    /// The source is a `const` member that ignores `this`; this is an
    /// associated function.
    pub fn get_quick_lda_score(
        library_corr: f64,
        library_norm_manhattan: f64,
        norm_rt_score: f64,
        xcorr_coelution_score: f64,
        xcorr_shape_score: f64,
        log_sn_score: f64,
    ) -> f64 {
        library_corr * -0.5319046
            + library_norm_manhattan * 2.1643962
            + norm_rt_score * 8.0353047
            + xcorr_coelution_score * 0.1458914
            + xcorr_shape_score * -1.6901925
            + log_sn_score * -0.8002824
    }

    /// A quick LDA average model using only non-DIA scores: `library_corr`,
    /// `library_norm_manhattan`, `norm_rt_score`, `xcorr_coelution_score`,
    /// `xcorr_shape_score`, `log_sn_score` and `elution_model_fit_score`
    /// (0.91 TPR / 0.20 FDR); better when more negative (source
    /// `calculate_lda_prescore`).
    ///
    /// The source takes the scores as an argument and ignores `this`; every
    /// source caller passes the object itself, so this scores `self`.
    pub fn calculate_lda_prescore(&self) -> f64 {
        self.library_corr * -0.34664267
            + self.library_norm_manhattan * 2.98700722
            + self.norm_rt_score * 7.05496384
            + self.xcorr_coelution_score * 0.09445371
            + self.xcorr_shape_score * -5.71823862
            + self.log_sn_score * -0.72989582
            + self.elution_model_fit_score * 1.88443209
    }

    /// A manually derived model for single-transition peak groups using only
    /// `norm_rt_score`, `log_sn_score` and `elution_model_fit_score`; better
    /// when more negative (source `calculate_lda_single_transition`; see
    /// [`OpenSwathScores::calculate_lda_prescore`] for the receiver).
    pub fn calculate_lda_single_transition(&self) -> f64 {
        self.norm_rt_score * 7.05496384
            + self.log_sn_score * -0.72989582
            + self.elution_model_fit_score * -1.08443209
    }

    /// A full LDA average model using DIA and non-DIA scores (0.76 TPR / 0.20
    /// FDR, without elution model); better when more negative (source
    /// `calculate_swath_lda_prescore`; see
    /// [`OpenSwathScores::calculate_lda_prescore`] for the receiver).
    pub fn calculate_swath_lda_prescore(&self) -> f64 {
        self.library_corr * -0.19011762
            + self.library_norm_manhattan * 2.47298914
            + self.norm_rt_score * 5.63906731
            + self.isotope_correlation * -0.62640133
            + self.isotope_overlap * 0.36006925
            + self.massdev_score * 0.08814003
            + self.xcorr_coelution_score * 0.13978311
            + self.xcorr_shape_score * -1.16475032
            + self.yseries_score * -0.19267813
            + self.log_sn_score * -0.61712054
    }
}
