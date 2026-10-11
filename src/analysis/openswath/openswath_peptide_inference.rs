// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSWATH peptide-level context inference.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathPeptideInference.h` and its translation unit: a
//! level check in front of
//! [`LevelContextInference::infer`](crate::analysis::openswath::level_context_inference::LevelContextInference::infer).
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_INFERENCE_SUPPORT.md).

use crate::analysis::openswath::level_context_inference::LevelContextInference;
use crate::analysis::openswath::openswath_inference_config::LevelContextInferenceConfig;
use crate::analysis::openswath::openswath_inference_data::{
    InferenceLevel, LevelContextInputRow, LevelContextResultRow,
};
use crate::{Error, Result};

/// OpenSWATH peptide-level context inference. Stateless.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OpenSwathPeptideInference;

impl OpenSwathPeptideInference {
    /// Infer peptide-level statistics for the provided rows.
    ///
    /// Delegates to
    /// [`LevelContextInference::infer`](crate::analysis::openswath::level_context_inference::LevelContextInference::infer);
    /// see there for the rows, the result and the statistics.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] (the source's
    /// `Exception::Precondition`) when `config.level` is not
    /// [`InferenceLevel::Peptide`], checked before anything else, even for an
    /// empty input; otherwise every error of `LevelContextInference::infer`.
    pub fn infer(
        &self,
        input: &[LevelContextInputRow],
        config: &LevelContextInferenceConfig,
    ) -> Result<Vec<LevelContextResultRow>> {
        if config.level != InferenceLevel::Peptide {
            return Err(Error::InvalidValue(
                "OpenSwathPeptideInference requires config.level = InferenceLevel::Peptide."
                    .to_string(),
            ));
        }
        LevelContextInference::infer(input, config)
    }
}
