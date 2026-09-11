// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! FASTA -> modified tryptic peptides -> fragment matching -> identification/coverage.
//! This fixed synthetic example ranks every candidate by fragment alignment; it
//! does not implement a search engine, precursor filtering or confidence/FDR estimation.

use openms::chemistry::{
    AASequence, ModifiedPeptideGenerator, ProteaseDigestion, TheoreticalSpectrumGenerator,
};
use openms::comparison::{SpectrumAlignment, SpectrumAlignmentScore, Tolerance};
use openms::format::fasta::FastaReader;
use openms::identification::{
    EnzymeTermSpecificity, FlankingResidue, PeakAnnotation, PeptideEvidence, PeptideHit,
    PeptideIdentification, ProteinHit, ProteinIdentification, SearchParameters,
};
use openms::{Error, MSSpectrum, Peak1D, Precursor, Result};
use std::collections::BTreeSet;
use std::io::Cursor;

const FASTA: &str = ">P_DEMO Synthetic example protein\nMPEPTIDERACDMKAGHIK\n\
                     >P_OTHER Synthetic distractor protein\nVVVVK\n";

/// The ranked candidates, annotated observed spectrum and protein records.
pub struct DemoIdentification {
    pub candidates: PeptideIdentification,
    pub spectrum: MSSpectrum,
    pub proteins: ProteinIdentification,
}

/// Run the deterministic example; retained records are also used by its workflow test.
pub fn identify_demo() -> Result<DemoIdentification> {
    // Seven singly charged b/y peaks of AC(Carbamidomethyl)DMK, rounded and
    // shifted by 0.001--0.003 Da, plus two unrelated peaks. This literal fixture
    // is independent of the candidate spectra generated below. Intensities are 1.
    let mut spectrum = MSSpectrum {
        peaks: [
            147.114804, 180.010000, 232.073039, 278.156289, 347.100982, 393.182232, 410.120000,
            478.139466, 553.213880,
        ]
        .into_iter()
        .map(|mz| Peak1D::new(mz, 1.0))
        .collect(),
        rt: 120.0,
        ms_level: 2,
        native_id: "scan=1".into(),
        precursors: vec![Precursor {
            mz: 312.627635,
            intensity: 1.0,
            charge: 2,
            ..Precursor::default()
        }],
        metadata: [("sample".into(), "synthetic digest".into())].into(),
        ..Default::default()
    };
    let digest = ProteaseDigestion::default();
    let modifications = ModifiedPeptideGenerator::default();
    let fixed_modifications =
        ModifiedPeptideGenerator::get_modifications(&["Carbamidomethyl (C)"])?;
    let generator = TheoreticalSpectrumGenerator {
        add_metainfo: true,
        ..Default::default()
    };
    let scorer = SpectrumAlignmentScore {
        alignment: SpectrumAlignment {
            tolerance: Tolerance::Absolute(0.02),
            ..Default::default()
        },
        ..Default::default()
    };
    let mut candidates = PeptideIdentification {
        identifier: "synthetic-run".into(),
        score_type: "SpectrumAlignmentScore".into(),
        rt: Some(spectrum.rt),
        mz: Some(spectrum.precursors[0].mz),
        // The bridge is explicit: legacy spectrum strings become typed strings.
        metadata: spectrum.metadata.clone(),
        ..Default::default()
    };
    candidates.set_spectrum_reference(spectrum.native_id.clone());
    let mut proteins = ProteinIdentification {
        identifier: candidates.identifier.clone(),
        search_engine: "identify_peptides demonstration".into(),
        search_parameters: SearchParameters {
            database: "embedded synthetic FASTA".into(),
            digestion_enzyme: "Trypsin".into(),
            enzyme_specificity: EnzymeTermSpecificity::Full,
            fixed_modifications: vec!["Carbamidomethyl (C)".into()],
            charges: "2".into(),
            fragment_tolerance: scorer.alignment.tolerance,
            ..Default::default()
        },
        ..Default::default()
    };

    for entry in FastaReader::new(Cursor::new(FASTA)) {
        let entry = entry?;
        let mut protein = AASequence::parse(&entry.sequence)?;
        // Fixed cysteine alkylation is applied before digestion. Subsequence
        // extraction carries each modification into the appropriate peptide.
        modifications.apply_fixed_modifications(&fixed_modifications, &mut protein)?;
        let mut protein_hit = ProteinHit::new(0.0, 0, &entry.identifier, &entry.sequence)?;
        protein_hit.set_description(&entry.description);
        proteins.hits.push(protein_hit);

        for peptide in digest.digest(&protein)? {
            let theoretical = generator.generate(&peptide.sequence, 1, 1, Some(2))?;
            let score = scorer.score(&theoretical, &spectrum)?;
            let matches = scorer.alignment.align(&theoretical, &spectrum)?;
            let names = &theoretical
                .string_data_arrays
                .iter()
                .find(|array| array.name == "IonNames")
                .ok_or_else(|| Error::InvalidValue("missing generated ion names".into()))?
                .data;
            let charges = &theoretical
                .integer_data_arrays
                .iter()
                .find(|array| array.name == "Charges")
                .ok_or_else(|| Error::InvalidValue("missing generated ion charges".into()))?
                .data;
            // Digestion ends are exclusive; PeptideEvidence ends are inclusive.
            let mut evidence =
                PeptideEvidence::new(&entry.identifier, peptide.start..=peptide.end - 1)?;
            evidence.aa_before = if peptide.start == 0 {
                FlankingResidue::NTerminus
            } else {
                FlankingResidue::Residue(entry.sequence.as_bytes()[peptide.start - 1] as char)
            };
            evidence.aa_after = if peptide.end == entry.sequence.len() {
                FlankingResidue::CTerminus
            } else {
                FlankingResidue::Residue(entry.sequence.as_bytes()[peptide.end] as char)
            };
            let mut hit = PeptideHit::new(score, 0, 2, peptide.sequence)?;
            hit.evidences.push(evidence);
            // Each pair indexes the theoretical spectrum first, observed second.
            for (theoretical_index, observed_index) in matches {
                let observed = spectrum.peaks[observed_index];
                hit.peak_annotations.push(PeakAnnotation {
                    mz: observed.mz,
                    intensity: f64::from(observed.intensity),
                    charge: charges[theoretical_index],
                    annotation: names[theoretical_index].clone(),
                });
            }
            hit.metadata.insert("fragment_charge".into(), 1_i64.into());
            candidates.hits.push(hit);
        }
    }
    candidates.sort()?;
    let mut selected = candidates.clone();
    selected.hits.truncate(1);
    let selected_hit = selected
        .hits
        .first_mut()
        .ok_or_else(|| Error::InvalidValue("example produced no peptide candidates".into()))?;
    selected_hit.rank = 1;
    selected.metadata.insert(
        "selection".into(),
        "highest fragment score in synthetic example".into(),
    );
    spectrum.peptide_identifications.push(selected);

    // Coverage consumes every supplied hit: only the selected identification is
    // passed here. Candidate alternatives must not inflate covered residues.
    proteins.compute_coverage(&spectrum.peptide_identifications)?;
    proteins.compute_modifications(&spectrum.peptide_identifications, &BTreeSet::new())?;
    spectrum.validate()?;
    proteins.validate()?;
    Ok(DemoIdentification {
        candidates,
        spectrum,
        proteins,
    })
}

#[cfg(not(test))]
fn main() -> Result<()> {
    use std::io::{BufWriter, Write};

    let result = identify_demo()?;
    let mut output = BufWriter::new(std::io::stdout().lock());
    writeln!(
        output,
        "Synthetic identification example; scores do not estimate confidence or FDR."
    )?;
    writeln!(output, "\npeptide\tfragment_score\tmatched_ions")?;
    for hit in &result.candidates.hits {
        writeln!(
            output,
            "{}\t{:.6}\t{}",
            hit.sequence,
            hit.score,
            hit.peak_annotations.len()
        )?;
    }
    let id = &result.spectrum.peptide_identifications[0];
    let hit = &id.hits[0];
    writeln!(
        output,
        "\nselected\t{}\tspectrum={}\tcharge={}",
        hit.sequence,
        id.spectrum_reference(),
        hit.charge
    )?;
    for protein in &result.proteins.hits {
        writeln!(
            output,
            "protein\t{}\tcoverage={:.2}%\t{}",
            protein.accession,
            protein.coverage.unwrap_or(0.0),
            protein.description()
        )?;
    }
    output.flush()?;
    Ok(())
}
