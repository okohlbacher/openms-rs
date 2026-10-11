// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Write scored OpenSWATH feature and transition tables to Parquet.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathParquetExporter.h` and its translation
//! unit: one row per [`OpenSwathFeatureScoreRow`] or
//! [`OpenSwathTransitionScoreRow`], in order, written with
//! [`write_table_to_parquet`] and its default [`ParquetWriteConfig`] (the
//! source's `ArrowIOHelpers::writeTableToParquet`). The column layout is
//! fixed except for the dynamically discovered `FEATURE_MS1_<name>`,
//! `FEATURE_MS2_<name>` and `FEATURE_TRANSITION_<name>` columns, which the
//! table names. An empty text is written as null, an absent optional value
//! as null, and a value vector shorter than its column list leaves the
//! remaining columns null. Integer charges, ordinals and ranks are widened to
//! `int64`. An empty table is written with its full schema.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_PARQUET_SUPPORT.md).
//!
//! # Differences from the source
//!
//! Every failure to write is reported as [`Error::Io`](crate::Error::Io), carrying the
//! underlying reason; the source maps the `false` of `writeTableToParquet` to
//! `Exception::FileNotWritable` and loses it.
//!
//! [`OpenSwathFeatureScoreRow`]: crate::analysis::openswath::openswath_export_data::OpenSwathFeatureScoreRow
//! [`OpenSwathTransitionScoreRow`]: crate::analysis::openswath::openswath_export_data::OpenSwathTransitionScoreRow
//! [`write_table_to_parquet`]: crate::format::arrow_io_helpers::write_table_to_parquet
//! [`ParquetWriteConfig`]: crate::format::arrow_io_helpers::ParquetWriteConfig

use crate::analysis::openswath::openswath_export_data::{
    OpenSwathFeatureScoreRow, OpenSwathFeatureScoreTable, OpenSwathTransitionScoreRow,
    OpenSwathTransitionScoreTable,
};
use crate::format::arrow_io_helpers::{ParquetWriteConfig, write_table_to_parquet};
use crate::{Error, Result};
use arrow_array::{ArrayRef, BooleanArray, Float64Array, Int64Array, RecordBatch, StringArray};
use arrow_schema::{DataType, Field, Schema};
use std::path::Path;
use std::sync::Arc;

/// Writer of the OpenSWATH Parquet score exports. Stateless; the source's
/// class has only static members.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OpenSwathParquetExporter;

/// Columns of a table under construction, in order.
#[derive(Default)]
struct Columns {
    fields: Vec<Field>,
    arrays: Vec<ArrayRef>,
}

impl Columns {
    fn push(&mut self, name: &str, data_type: DataType, array: ArrayRef) {
        self.fields.push(Field::new(name, data_type, true));
        self.arrays.push(array);
    }

    fn i64<T>(&mut self, name: &str, rows: &[T], value: impl Fn(&T) -> Option<i64>) {
        let array = Int64Array::from(rows.iter().map(value).collect::<Vec<_>>());
        self.push(name, DataType::Int64, Arc::new(array));
    }

    fn f64<T>(&mut self, name: &str, rows: &[T], value: impl Fn(&T) -> Option<f64>) {
        let array = Float64Array::from(rows.iter().map(value).collect::<Vec<_>>());
        self.push(name, DataType::Float64, Arc::new(array));
    }

    fn bool<T>(&mut self, name: &str, rows: &[T], value: impl Fn(&T) -> Option<bool>) {
        let array = BooleanArray::from(rows.iter().map(value).collect::<Vec<_>>());
        self.push(name, DataType::Boolean, Arc::new(array));
    }

    /// The source's `appendString_`: an empty text is null.
    fn text<T>(&mut self, name: &str, rows: &[T], value: impl Fn(&T) -> &str) {
        let array = StringArray::from(
            rows.iter()
                .map(|row| Some(value(row)).filter(|s| !s.is_empty()))
                .collect::<Vec<_>>(),
        );
        self.push(name, DataType::Utf8, Arc::new(array));
    }

    fn write(self, filename: &Path) -> Result<()> {
        let table = RecordBatch::try_new(Arc::new(Schema::new(self.fields)), self.arrays)
            .map_err(|e| Error::InvalidValue(format!("Failed to build export table: {e}")))?;
        write_table_to_parquet(&table, filename, &ParquetWriteConfig::default()).map_err(|e| {
            Error::Io(std::io::Error::other(format!(
                "the file '{}' could not be created: {e}",
                filename.display()
            )))
        })
    }
}

/// A feature-level score column: (name, value).
type FeatureScore = (&'static str, fn(&OpenSwathFeatureScoreRow) -> Option<f64>);

/// The feature-level score columns after the MS1 and MS2 ranks, in the
/// source's order.
const PEPTIDE_PROTEIN_GENE_SCORES: [FeatureScore; 36] = [
    ("SCORE_PEPTIDE_GLOBAL_SCORE", |r| {
        r.score_peptide_global_score
    }),
    ("SCORE_PEPTIDE_GLOBAL_PVALUE", |r| {
        r.score_peptide_global_pvalue
    }),
    ("SCORE_PEPTIDE_GLOBAL_QVALUE", |r| {
        r.score_peptide_global_qvalue
    }),
    ("SCORE_PEPTIDE_GLOBAL_PEP", |r| r.score_peptide_global_pep),
    ("SCORE_PEPTIDE_EXPERIMENT_WIDE_SCORE", |r| {
        r.score_peptide_experiment_wide_score
    }),
    ("SCORE_PEPTIDE_EXPERIMENT_WIDE_PVALUE", |r| {
        r.score_peptide_experiment_wide_pvalue
    }),
    ("SCORE_PEPTIDE_EXPERIMENT_WIDE_QVALUE", |r| {
        r.score_peptide_experiment_wide_qvalue
    }),
    ("SCORE_PEPTIDE_EXPERIMENT_WIDE_PEP", |r| {
        r.score_peptide_experiment_wide_pep
    }),
    ("SCORE_PEPTIDE_RUN_SPECIFIC_SCORE", |r| {
        r.score_peptide_run_specific_score
    }),
    ("SCORE_PEPTIDE_RUN_SPECIFIC_PVALUE", |r| {
        r.score_peptide_run_specific_pvalue
    }),
    ("SCORE_PEPTIDE_RUN_SPECIFIC_QVALUE", |r| {
        r.score_peptide_run_specific_qvalue
    }),
    ("SCORE_PEPTIDE_RUN_SPECIFIC_PEP", |r| {
        r.score_peptide_run_specific_pep
    }),
    ("SCORE_PROTEIN_GLOBAL_SCORE", |r| {
        r.score_protein_global_score
    }),
    ("SCORE_PROTEIN_GLOBAL_PVALUE", |r| {
        r.score_protein_global_pvalue
    }),
    ("SCORE_PROTEIN_GLOBAL_QVALUE", |r| {
        r.score_protein_global_qvalue
    }),
    ("SCORE_PROTEIN_GLOBAL_PEP", |r| r.score_protein_global_pep),
    ("SCORE_PROTEIN_EXPERIMENT_WIDE_SCORE", |r| {
        r.score_protein_experiment_wide_score
    }),
    ("SCORE_PROTEIN_EXPERIMENT_WIDE_PVALUE", |r| {
        r.score_protein_experiment_wide_pvalue
    }),
    ("SCORE_PROTEIN_EXPERIMENT_WIDE_QVALUE", |r| {
        r.score_protein_experiment_wide_qvalue
    }),
    ("SCORE_PROTEIN_EXPERIMENT_WIDE_PEP", |r| {
        r.score_protein_experiment_wide_pep
    }),
    ("SCORE_PROTEIN_RUN_SPECIFIC_SCORE", |r| {
        r.score_protein_run_specific_score
    }),
    ("SCORE_PROTEIN_RUN_SPECIFIC_PVALUE", |r| {
        r.score_protein_run_specific_pvalue
    }),
    ("SCORE_PROTEIN_RUN_SPECIFIC_QVALUE", |r| {
        r.score_protein_run_specific_qvalue
    }),
    ("SCORE_PROTEIN_RUN_SPECIFIC_PEP", |r| {
        r.score_protein_run_specific_pep
    }),
    ("SCORE_GENE_GLOBAL_SCORE", |r| r.score_gene_global_score),
    ("SCORE_GENE_GLOBAL_PVALUE", |r| r.score_gene_global_pvalue),
    ("SCORE_GENE_GLOBAL_QVALUE", |r| r.score_gene_global_qvalue),
    ("SCORE_GENE_GLOBAL_PEP", |r| r.score_gene_global_pep),
    ("SCORE_GENE_EXPERIMENT_WIDE_SCORE", |r| {
        r.score_gene_experiment_wide_score
    }),
    ("SCORE_GENE_EXPERIMENT_WIDE_PVALUE", |r| {
        r.score_gene_experiment_wide_pvalue
    }),
    ("SCORE_GENE_EXPERIMENT_WIDE_QVALUE", |r| {
        r.score_gene_experiment_wide_qvalue
    }),
    ("SCORE_GENE_EXPERIMENT_WIDE_PEP", |r| {
        r.score_gene_experiment_wide_pep
    }),
    ("SCORE_GENE_RUN_SPECIFIC_SCORE", |r| {
        r.score_gene_run_specific_score
    }),
    ("SCORE_GENE_RUN_SPECIFIC_PVALUE", |r| {
        r.score_gene_run_specific_pvalue
    }),
    ("SCORE_GENE_RUN_SPECIFIC_QVALUE", |r| {
        r.score_gene_run_specific_qvalue
    }),
    ("SCORE_GENE_RUN_SPECIFIC_PEP", |r| {
        r.score_gene_run_specific_pep
    }),
];

impl OpenSwathParquetExporter {
    /// Write the feature-level score rows to `filename`.
    ///
    /// Ports `writeFeatureScores`. The columns, in order: `PROTEIN_ID`,
    /// `PEPTIDE_ID`, `IPF_PEPTIDE_ID`, `PRECURSOR_ID`, `PROTEIN_ACCESSION`,
    /// `UNMODIFIED_SEQUENCE`, `MODIFIED_SEQUENCE`, `PRECURSOR_TRAML_ID`,
    /// `PRECURSOR_GROUP_LABEL`, `PRECURSOR_MZ`, `PRECURSOR_CHARGE`,
    /// `PRECURSOR_LIBRARY_INTENSITY`, `PRECURSOR_LIBRARY_RT`,
    /// `PRECURSOR_LIBRARY_DRIFT_TIME`, `GENE_ID`, `GENE_NAME`, `GENE_DECOY`,
    /// `PROTEIN_DECOY`, `PEPTIDE_DECOY`, `PRECURSOR_DECOY`, `RUN_ID`,
    /// `FILENAME`, `FEATURE_ID`, `EXP_RT`, `EXP_IM`, `NORM_RT`, `DELTA_RT`,
    /// `LEFT_WIDTH`, `RIGHT_WIDTH`, `IM_leftWidth`, `IM_rightWidth`, one
    /// `FEATURE_MS1_<name>` per MS1 column name and one `FEATURE_MS2_<name>`
    /// per MS2 column name, `SCORE_MS1_SCORE`, `SCORE_MS1_RANK`,
    /// `SCORE_MS1_PVALUE`, `SCORE_MS1_QVALUE`, `SCORE_MS1_PEP`,
    /// `SCORE_MS2_SCORE`, `SCORE_MS2_PEAK_GROUP_RANK`, `SCORE_MS2_PVALUE`,
    /// `SCORE_MS2_QVALUE`, `SCORE_MS2_PEP`,
    /// `SCORE_IPF_PRECURSOR_PEAKGROUP_PEP`, `SCORE_IPF_PEP`,
    /// `SCORE_IPF_QVALUE`, the score, p-value, q-value and PEP of the
    /// peptide, protein and gene levels in the global, experiment-wide and
    /// run-specific contexts, `alignment_group_id`,
    /// `alignment_reference_feature_id`, `alignment_reference_rt`,
    /// `alignment_pep`, `alignment_qvalue` and `from_alignment`.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be written (the source's
    /// `Exception::FileNotWritable`).
    pub fn write_feature_scores(
        filename: impl AsRef<Path>,
        table: &OpenSwathFeatureScoreTable,
    ) -> Result<()> {
        let rows: &[OpenSwathFeatureScoreRow] = &table.rows;
        let mut c = Columns::default();
        c.i64("PROTEIN_ID", rows, |r| Some(r.protein_id));
        c.i64("PEPTIDE_ID", rows, |r| Some(r.peptide_id));
        c.i64("IPF_PEPTIDE_ID", rows, |r| r.ipf_peptide_id);
        c.i64("PRECURSOR_ID", rows, |r| Some(r.precursor_id));
        c.text("PROTEIN_ACCESSION", rows, |r| &r.protein_accession);
        c.text("UNMODIFIED_SEQUENCE", rows, |r| &r.unmodified_sequence);
        c.text("MODIFIED_SEQUENCE", rows, |r| &r.modified_sequence);
        c.text("PRECURSOR_TRAML_ID", rows, |r| &r.precursor_traml_id);
        c.text("PRECURSOR_GROUP_LABEL", rows, |r| &r.precursor_group_label);
        c.f64("PRECURSOR_MZ", rows, |r| Some(r.precursor_mz));
        c.i64("PRECURSOR_CHARGE", rows, |r| {
            Some(i64::from(r.precursor_charge))
        });
        c.f64("PRECURSOR_LIBRARY_INTENSITY", rows, |r| {
            r.precursor_library_intensity
        });
        c.f64("PRECURSOR_LIBRARY_RT", rows, |r| r.precursor_library_rt);
        c.f64("PRECURSOR_LIBRARY_DRIFT_TIME", rows, |r| {
            r.precursor_library_drift_time
        });
        c.i64("GENE_ID", rows, |r| r.gene_id);
        c.text("GENE_NAME", rows, |r| &r.gene_name);
        c.bool("GENE_DECOY", rows, |r| r.gene_decoy);
        c.bool("PROTEIN_DECOY", rows, |r| Some(r.protein_decoy));
        c.bool("PEPTIDE_DECOY", rows, |r| Some(r.peptide_decoy));
        c.bool("PRECURSOR_DECOY", rows, |r| Some(r.precursor_decoy));
        c.i64("RUN_ID", rows, |r| Some(r.run_id));
        c.text("FILENAME", rows, |r| &r.filename);
        c.i64("FEATURE_ID", rows, |r| Some(r.feature_id));
        c.f64("EXP_RT", rows, |r| Some(r.exp_rt));
        c.f64("EXP_IM", rows, |r| r.exp_im);
        c.f64("NORM_RT", rows, |r| Some(r.norm_rt));
        c.f64("DELTA_RT", rows, |r| Some(r.delta_rt));
        c.f64("LEFT_WIDTH", rows, |r| Some(r.left_width));
        c.f64("RIGHT_WIDTH", rows, |r| Some(r.right_width));
        c.f64("IM_leftWidth", rows, |r| r.im_left_width);
        c.f64("IM_rightWidth", rows, |r| r.im_right_width);
        for (i, name) in table.feature_ms1_column_names.iter().enumerate() {
            c.f64(&format!("FEATURE_MS1_{name}"), rows, |r| {
                r.feature_ms1_values.get(i).copied()
            });
        }
        for (i, name) in table.feature_ms2_column_names.iter().enumerate() {
            c.f64(&format!("FEATURE_MS2_{name}"), rows, |r| {
                r.feature_ms2_values.get(i).copied()
            });
        }
        c.f64("SCORE_MS1_SCORE", rows, |r| r.score_ms1_score);
        c.i64("SCORE_MS1_RANK", rows, |r| r.score_ms1_rank.map(i64::from));
        c.f64("SCORE_MS1_PVALUE", rows, |r| r.score_ms1_pvalue);
        c.f64("SCORE_MS1_QVALUE", rows, |r| r.score_ms1_qvalue);
        c.f64("SCORE_MS1_PEP", rows, |r| r.score_ms1_pep);
        c.f64("SCORE_MS2_SCORE", rows, |r| r.score_ms2_score);
        c.i64("SCORE_MS2_PEAK_GROUP_RANK", rows, |r| {
            r.score_ms2_peak_group_rank.map(i64::from)
        });
        c.f64("SCORE_MS2_PVALUE", rows, |r| r.score_ms2_pvalue);
        c.f64("SCORE_MS2_QVALUE", rows, |r| r.score_ms2_qvalue);
        c.f64("SCORE_MS2_PEP", rows, |r| r.score_ms2_pep);
        c.f64("SCORE_IPF_PRECURSOR_PEAKGROUP_PEP", rows, |r| {
            r.score_ipf_precursor_peakgroup_pep
        });
        c.f64("SCORE_IPF_PEP", rows, |r| r.score_ipf_pep);
        c.f64("SCORE_IPF_QVALUE", rows, |r| r.score_ipf_qvalue);
        for (name, value) in PEPTIDE_PROTEIN_GENE_SCORES {
            c.f64(name, rows, value);
        }
        c.i64("alignment_group_id", rows, |r| r.alignment_group_id);
        c.i64("alignment_reference_feature_id", rows, |r| {
            r.alignment_reference_feature_id
        });
        c.f64("alignment_reference_rt", rows, |r| r.alignment_reference_rt);
        c.f64("alignment_pep", rows, |r| r.alignment_pep);
        c.f64("alignment_qvalue", rows, |r| r.alignment_qvalue);
        c.bool("from_alignment", rows, |r| Some(r.from_alignment));
        c.write(filename.as_ref())
    }

    /// Write the transition-level score rows to `filename`.
    ///
    /// Ports `writeTransitionScores`. The columns, in order: `RUN_ID`,
    /// `IPF_PEPTIDE_ID`, `PRECURSOR_ID`, `TRANSITION_ID`,
    /// `TRANSITION_TRAML_ID`, `PRODUCT_MZ`, `TRANSITION_CHARGE`,
    /// `TRANSITION_TYPE`, `TRANSITION_ORDINAL`, `ANNOTATION`,
    /// `TRANSITION_DETECTING`, `TRANSITION_LIBRARY_INTENSITY`,
    /// `TRANSITION_DECOY`, `FEATURE_ID`, one `FEATURE_TRANSITION_<name>` per
    /// column name, `SCORE_TRANSITION_SCORE`, `SCORE_TRANSITION_RANK`,
    /// `SCORE_TRANSITION_PVALUE`, `SCORE_TRANSITION_QVALUE` and
    /// `SCORE_TRANSITION_PEP`.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be written (the source's
    /// `Exception::FileNotWritable`).
    pub fn write_transition_scores(
        filename: impl AsRef<Path>,
        table: &OpenSwathTransitionScoreTable,
    ) -> Result<()> {
        let rows: &[OpenSwathTransitionScoreRow] = &table.rows;
        let mut c = Columns::default();
        c.i64("RUN_ID", rows, |r| r.run_id);
        c.i64("IPF_PEPTIDE_ID", rows, |r| r.ipf_peptide_id);
        c.i64("PRECURSOR_ID", rows, |r| Some(r.precursor_id));
        c.i64("TRANSITION_ID", rows, |r| Some(r.transition_id));
        c.text("TRANSITION_TRAML_ID", rows, |r| &r.transition_traml_id);
        c.f64("PRODUCT_MZ", rows, |r| Some(r.product_mz));
        c.i64("TRANSITION_CHARGE", rows, |r| {
            Some(i64::from(r.transition_charge))
        });
        c.text("TRANSITION_TYPE", rows, |r| &r.transition_type);
        c.i64("TRANSITION_ORDINAL", rows, |r| {
            Some(i64::from(r.transition_ordinal))
        });
        c.text("ANNOTATION", rows, |r| &r.annotation);
        c.bool("TRANSITION_DETECTING", rows, |r| {
            Some(r.transition_detecting)
        });
        c.f64("TRANSITION_LIBRARY_INTENSITY", rows, |r| {
            r.transition_library_intensity
        });
        c.bool("TRANSITION_DECOY", rows, |r| Some(r.transition_decoy));
        c.i64("FEATURE_ID", rows, |r| r.feature_id);
        for (i, name) in table.feature_transition_column_names.iter().enumerate() {
            c.f64(&format!("FEATURE_TRANSITION_{name}"), rows, |r| {
                r.feature_transition_values.get(i).copied().flatten()
            });
        }
        c.f64("SCORE_TRANSITION_SCORE", rows, |r| r.score_transition_score);
        c.i64("SCORE_TRANSITION_RANK", rows, |r| {
            r.score_transition_rank.map(i64::from)
        });
        c.f64("SCORE_TRANSITION_PVALUE", rows, |r| {
            r.score_transition_pvalue
        });
        c.f64("SCORE_TRANSITION_QVALUE", rows, |r| {
            r.score_transition_qvalue
        });
        c.f64("SCORE_TRANSITION_PEP", rows, |r| r.score_transition_pep);
        c.write(filename.as_ref())
    }
}
