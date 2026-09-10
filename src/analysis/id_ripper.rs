// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Partition merged identifications by their annotated input file.
//! This is an in-memory port of IDRipper; no files are opened or written.

use crate::identification::{PeptideIdentification, ProteinIdentification};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OriginAnnotationFormat {
    FileOrigin,
    MapIndex,
    IdMergeIndex,
}
impl OriginAnnotationFormat {
    pub fn key(self) -> &'static str {
        match self {
            Self::FileOrigin => "file_origin",
            Self::MapIndex => "map_index",
            Self::IdMergeIndex => "id_merge_index",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RipFileIdentifier {
    /// None combines runs. Some preserves the input protein-run index.
    pub identification_run_index: Option<usize>,
    pub file_origin_index: usize,
    pub output_basename: String,
    pub origin_fullname: String,
}
#[derive(Clone, Debug, PartialEq)]
pub struct RippedFile {
    pub identifier: RipFileIdentifier,
    pub protein_identifications: Vec<ProteinIdentification>,
    pub peptide_identifications: Vec<PeptideIdentification>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct RippingResult {
    pub origin_format: OriginAnnotationFormat,
    /// Sorted by run index and then file-origin index, like the source map.
    pub files: Vec<RippedFile>,
    /// Empty IDs and IDs without protein accessions, omitted by source IDRipper.
    /// Retained here with original metadata so callers can account for them.
    pub skipped_peptide_identifications: Vec<PeptideIdentification>,
    /// Groups omitted across output run copies. Input groups remain available.
    pub groups_not_copied: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IDRipper {
    /// Allow repeated basename strings when numeric keys distinguish outputs.
    pub numeric_filenames: bool,
    pub split_identification_runs: bool,
}
fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

/// Require one common, unambiguous source annotation on every ID, including
/// empty IDs. The source `unknown` annotation is explicitly unsupported.
pub fn detect_origin_annotation_format(
    ids: &[PeptideIdentification],
) -> Result<OriginAnnotationFormat> {
    let mut mode = None;
    for id in ids {
        id.validate()?;
        if id.metadata.contains_key("unknown") {
            return Err(bad("unknown origin annotation format"));
        }
        let formats: Vec<_> = [
            OriginAnnotationFormat::FileOrigin,
            OriginAnnotationFormat::MapIndex,
            OriginAnnotationFormat::IdMergeIndex,
        ]
        .into_iter()
        .filter(|f| id.metadata.contains_key(f.key()))
        .collect();
        if formats.len() != 1 || mode.is_some_and(|mode| mode != formats[0]) {
            return Err(bad(
                "identifications require one consistent origin annotation format",
            ));
        }
        mode = Some(formats[0]);
    }
    mode.ok_or_else(|| bad("cannot detect origin annotation on an empty identification list"))
}

impl IDRipper {
    /// Split owned output copies without changing either input. Unknown run or
    /// protein references and output-name collisions return an error, with no
    /// partial result. Selected protein hits are resolved within their own run.
    pub fn rip(
        &self,
        proteins: &[ProteinIdentification],
        peptides: &[PeptideIdentification],
    ) -> Result<RippingResult> {
        let format = detect_origin_annotation_format(peptides)?;
        let mut run_indices = BTreeMap::new();
        let mut hits_by_run = Vec::with_capacity(proteins.len());
        for (i, run) in proteins.iter().enumerate() {
            run.validate()?;
            if run_indices.insert(run.identifier.as_str(), i).is_some() {
                return Err(bad("identification run IDs must be unique"));
            }
            let mut hits = BTreeMap::new();
            for hit in &run.hits {
                if hits.insert(hit.accession.as_str(), hit).is_some() {
                    return Err(bad("protein accessions must be unique within a run"));
                }
            }
            hits_by_run.push(hits);
        }
        let mut origins = BTreeMap::new();
        if format == OriginAnnotationFormat::FileOrigin {
            for id in peptides {
                let origin = id.metadata[format.key()].as_str()?;
                let index = origins.len();
                origins.entry(origin).or_insert(index);
            }
        }
        let mut files = BTreeMap::<(Option<usize>, usize), RippedFile>::new();
        let mut basenames = BTreeMap::new();
        let mut names_by_key = BTreeMap::new();
        let mut skipped = Vec::new();
        let mut groups_not_copied = 0;
        for id in peptides {
            let run_index = *run_indices
                .get(id.identifier.as_str())
                .ok_or_else(|| bad("unknown peptide identification run"))?;
            let run = &proteins[run_index];
            let (origin_index, origin) = match format {
                OriginAnnotationFormat::FileOrigin => {
                    let origin = id.metadata[format.key()].as_str()?;
                    (origins[origin], origin)
                }
                _ => {
                    // The source converts DataValue to text, then to signed int32.
                    let index = id.metadata[format.key()]
                        .to_string()
                        .trim()
                        .parse::<i32>()
                        .ok()
                        .and_then(|index| usize::try_from(index).ok())
                        .ok_or_else(|| {
                            bad("origin index must be a nonnegative signed 32-bit integer")
                        })?;
                    let origin = run.primary_ms_run_paths.get(index).ok_or_else(|| {
                        bad("origin index has no corresponding spectra_data entry")
                    })?;
                    (index, origin.as_str())
                }
            };
            if origin.is_empty() {
                return Err(bad("file origin must not be empty"));
            }
            let basename = Path::new(origin)
                .file_stem()
                .and_then(|v| v.to_str())
                .unwrap_or("")
                .to_owned();
            if !self.numeric_filenames && basename.is_empty() {
                return Err(bad("file origin has no usable basename"));
            }
            let key = (
                self.split_identification_runs.then_some(run_index),
                origin_index,
            );
            if let Some(previous) = names_by_key.insert(key, origin) {
                if previous != origin {
                    return Err(bad(
                        "one numeric origin refers to different files; split identification runs",
                    ));
                }
            }
            if !self.numeric_filenames {
                if let Some(previous) = basenames.insert(basename.clone(), key) {
                    if previous != key {
                        return Err(bad(
                            "output basenames are not unique; use numeric filenames",
                        ));
                    }
                }
            }
            let accessions: BTreeSet<_> = id
                .hits
                .iter()
                .flat_map(|hit| hit.protein_accessions())
                .collect();
            if id.hits.is_empty() || accessions.is_empty() {
                skipped.push(id.clone());
                continue;
            }
            let mut selected = Vec::with_capacity(accessions.len());
            for accession in accessions {
                selected.push(
                    (*hits_by_run[run_index].get(accession).ok_or_else(|| {
                        bad("peptide evidence references a missing protein in its run")
                    })?)
                    .clone(),
                );
            }
            let file = files.entry(key).or_insert_with(|| RippedFile {
                identifier: RipFileIdentifier {
                    identification_run_index: key.0,
                    file_origin_index: key.1,
                    output_basename: basename,
                    origin_fullname: origin.to_owned(),
                },
                protein_identifications: Vec::new(),
                peptide_identifications: Vec::new(),
            });
            if let Some(output_run) = file
                .protein_identifications
                .iter_mut()
                .find(|p| p.identifier == run.identifier)
            {
                let mut known: BTreeSet<_> = output_run
                    .hits
                    .iter()
                    .map(|h| h.accession.clone())
                    .collect();
                for hit in selected {
                    if known.insert(hit.accession.clone()) {
                        output_run.hits.push(hit);
                    }
                }
            } else {
                let mut output_run = run.clone();
                output_run.hits = selected;
                output_run.metadata.remove(format.key());
                groups_not_copied +=
                    output_run.protein_groups.len() + output_run.indistinguishable_groups.len();
                output_run.protein_groups.clear();
                output_run.indistinguishable_groups.clear();
                if format != OriginAnnotationFormat::FileOrigin {
                    output_run.primary_ms_run_paths = vec![origin.to_owned()];
                }
                file.protein_identifications.push(output_run);
            }
            let mut output_id = id.clone();
            output_id.metadata.remove(format.key());
            file.peptide_identifications.push(output_id);
        }
        Ok(RippingResult {
            origin_format: format,
            files: files.into_values().collect(),
            skipped_peptide_identifications: skipped,
            groups_not_copied,
        })
    }
}
