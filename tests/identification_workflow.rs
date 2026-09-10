// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

// Exercise the executable's workflow itself so the example and its test agree.
#[path = "../examples/identify_peptides.rs"]
mod example;

use openms::identification::FlankingResidue;
use std::collections::BTreeSet;

#[test]
fn modified_digest_identifies_known_peptide_and_retains_evidence_and_metadata() {
    let result = example::identify_demo().unwrap();
    assert_eq!(result.candidates.hits.len(), 4);
    let candidate = &result.candidates.hits[0];
    assert_eq!(candidate.sequence.to_string(), "AC(Carbamidomethyl)DMK");
    assert!(candidate.score > result.candidates.hits[1].score);
    // Seven unit-intensity matches against nine observed and seven theoretical
    // peaks: this uses OpenMS's alignment score formula, not a probability.
    assert!((candidate.score - 7.0 / 63.0_f64.sqrt()).abs() < 1e-12);

    let spectrum = &result.spectrum;
    assert_eq!(spectrum.metadata["sample"], "synthetic digest");
    assert_eq!(spectrum.peptide_identifications.len(), 1);
    let id = &spectrum.peptide_identifications[0];
    assert_eq!(id.hits.len(), 1);
    assert_eq!(id.identifier, result.proteins.identifier);
    assert_eq!(id.spectrum_reference(), spectrum.native_id);
    assert_eq!(id.rt, Some(120.0));
    assert_eq!(id.mz, Some(spectrum.precursors[0].mz));
    assert_eq!(id.metadata["sample"].as_str().unwrap(), "synthetic digest");
    assert_eq!(
        id.metadata["selection"].as_str().unwrap(),
        "highest fragment score in synthetic example"
    );
    let hit = &id.hits[0];
    assert_eq!(hit.sequence, candidate.sequence);
    assert_eq!((hit.rank, hit.charge), (1, 2));
    assert!((hit.sequence.mz(hit.charge).unwrap() - spectrum.precursors[0].mz).abs() < 0.01);
    assert_eq!(hit.metadata["fragment_charge"].as_i64().unwrap(), 1);
    assert_eq!(hit.evidences.len(), 1);
    let evidence = &hit.evidences[0];
    assert_eq!(evidence.protein_accession, "P_DEMO");
    assert_eq!(evidence.positions().unwrap(), 9..=13);
    assert_eq!(evidence.aa_before, FlankingResidue::Residue('R'));
    assert_eq!(evidence.aa_after, FlankingResidue::Residue('A'));
    let ion_names: BTreeSet<_> = hit
        .peak_annotations
        .iter()
        .map(|annotation| annotation.annotation.as_str())
        .collect();
    assert_eq!(
        ion_names,
        BTreeSet::from(["b2+", "b3+", "b4+", "y1+", "y2+", "y3+", "y4+"])
    );
    assert!(hit.peak_annotations.iter().all(|annotation| {
        annotation.charge == 1
            && spectrum.peaks.iter().any(|peak| {
                peak.mz == annotation.mz && f64::from(peak.intensity) == annotation.intensity
            })
    }));

    let protein = result.proteins.find_hit("P_DEMO").unwrap();
    assert_eq!(protein.description(), "Synthetic example protein");
    assert!((protein.coverage.unwrap() - 100.0 * 5.0 / 19.0).abs() < 1e-12);
    assert_eq!(protein.modifications.len(), 1);
    assert_eq!(protein.modifications[0].position, 10);
    assert_eq!(
        protein.modifications[0].modification.full_id(),
        "Carbamidomethyl (C)"
    );
    assert_eq!(
        result.proteins.find_hit("P_OTHER").unwrap().coverage,
        Some(0.0)
    );
    assert_eq!(
        result.proteins.search_parameters.fixed_modifications,
        ["Carbamidomethyl (C)"]
    );
}
