// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#![cfg(feature = "idxml")]

use openms::analysis::false_discovery_rate::FalseDiscoveryRate;
use openms::analysis::id_filter::{
    filter_peptides_by_score, keep_n_best_peptide_hits, remove_dangling_protein_references,
    remove_decoy_protein_hits, remove_empty_peptide_identifications, remove_unreferenced_proteins,
    update_protein_groups,
};
use openms::analysis::scores::ScoreSwitcher;
use openms::format::idxml;
use openms::identification::{FlankingResidue, TargetDecoyType};
use std::io::Cursor;

// Independent source-format literal, not emitted by the Rust writer. Group
// and evidence encodings follow the pinned IdXMLFile fixtures. The deliberate
// tiny score population tests the source D/T calculation, not empirical confidence.
const INPUT: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<IdXML version="1.5" id="counting-example">
  <SearchParameters id="SP_0" db="synthetic.fasta" db_version="1" taxonomy="" mass_type="monoisotopic" charges="2" enzyme="Trypsin" missed_cleavages="0" precursor_peak_tolerance="0.01" peak_mass_tolerance="0.02">
    <FixedModification name="Carbamidomethyl (C)"/>
    <UserParam type="string" name="database_note" value="independent fixture"/>
  </SearchParameters>
  <IdentificationRun date="2026-09-10T12:00:00" search_engine="synthetic" search_engine_version="1" search_parameters_ref="SP_0">
    <ProteinIdentification score_type="Mascot" higher_score_better="true" significance_threshold="0">
      <ProteinHit id="PH_A" accession="A" score="10" sequence="PEPTIDEAACDMK">
        <UserParam type="string" name="target_decoy" value="target"/>
        <UserParam type="string" name="description" value="retained protein"/>
      </ProteinHit>
      <ProteinHit id="PH_B" accession="B" score="7" sequence="VVVVK">
        <UserParam type="string" name="target_decoy" value="target"/>
      </ProteinHit>
      <ProteinHit id="PH_D" accession="DECOY_A" score="8" sequence="PEPTIDEAACDMK">
        <UserParam type="string" name="target_decoy" value="decoy"/>
      </ProteinHit>
      <ProteinHit id="PH_U" accession="UNREFERENCED" score="0" sequence="AAAAK">
        <UserParam type="string" name="target_decoy" value="target"/>
      </ProteinHit>
      <UserParam type="string" name="operator_note" value="counted &amp; retained"/>
      <UserParam type="string" name="protein_group_0" value="0.8,PH_A,PH_B"/>
      <UserParam type="string" name="protein_group_1" value="0.2,PH_D,PH_U"/>
    </ProteinIdentification>
    <PeptideIdentification score_type="Mascot" higher_score_better="true" significance_threshold="0" RT="10" spectrum_reference="scan=1">
      <PeptideHit score="100" sequence="VVVVK" charge="2" protein_refs="PH_B" start="0" end="4" aa_before="[" aa_after="]">
        <UserParam type="float" name="hyperscore" value="1"/>
        <UserParam type="string" name="target_decoy" value="target"/>
      </PeptideHit>
      <PeptideHit score="1" sequence="PEPTIDEA" charge="2" protein_refs="PH_A PH_D" start="0 0" end="7 7" aa_before="[ [" aa_after="A A">
        <UserParam type="float" name="hyperscore" value="10"/>
        <UserParam type="string" name="target_decoy" value="target+decoy"/>
        <UserParam type="string" name="annotation_note" value="shared evidence"/>
      </PeptideHit>
      <UserParam type="string" name="sample" value="synthetic digest"/>
    </PeptideIdentification>
    <PeptideIdentification score_type="Mascot" higher_score_better="true" significance_threshold="0" RT="20" MZ="312.627635" spectrum_reference="scan=2">
      <PeptideHit score="2" sequence="AC(Carbamidomethyl)DMK" charge="2" protein_refs="PH_A" start="8" end="12" aa_before="A" aa_after="]">
        <UserParam type="float" name="hyperscore" value="9"/>
        <UserParam type="string" name="target_decoy" value="target"/>
        <UserParam type="intList" name="supporting_scans" value="[2, 12]"/>
      </PeptideHit>
    </PeptideIdentification>
    <PeptideIdentification score_type="Mascot" higher_score_better="true" significance_threshold="0" spectrum_reference="scan=3">
      <PeptideHit score="3" sequence="AAA" charge="2" protein_refs="PH_D">
        <UserParam type="float" name="hyperscore" value="8"/>
        <UserParam type="string" name="target_decoy" value="decoy"/>
      </PeptideHit>
    </PeptideIdentification>
    <PeptideIdentification score_type="Mascot" higher_score_better="true" significance_threshold="0" spectrum_reference="scan=4">
      <PeptideHit score="4" sequence="VVVVK" charge="2" protein_refs="PH_B" start="0" end="4">
        <UserParam type="float" name="hyperscore" value="7"/>
        <UserParam type="string" name="target_decoy" value="target"/>
      </PeptideHit>
    </PeptideIdentification>
    <PeptideIdentification score_type="Mascot" higher_score_better="true" significance_threshold="0" spectrum_reference="scan=5">
      <PeptideHit score="5" sequence="CCC" charge="2" protein_refs="PH_D">
        <UserParam type="float" name="hyperscore" value="6"/>
        <UserParam type="string" name="target_decoy" value="decoy"/>
      </PeptideHit>
    </PeptideIdentification>
  </IdentificationRun>
</IdXML>"#;

#[test]
fn source_counting_pipeline_preserves_selected_identifications_across_idxml_roundtrip() {
    let mut doc = idxml::read(Cursor::new(INPUT.as_bytes())).unwrap();
    assert_eq!(doc.peptide_identifications.len(), 5);
    assert_eq!(doc.protein_identifications[0].hits.len(), 4);
    assert_eq!(doc.protein_identifications[0].protein_groups.len(), 2);
    let run_id = doc.protein_identifications[0].identifier.clone();
    assert!(
        doc.peptide_identifications
            .iter()
            .all(|id| id.identifier == run_id)
    );

    // Select after switching: scan1's high Mascot-scored alternative must lose
    // to the independently supplied hyperscore of the true fixture candidate.
    assert_eq!(
        ScoreSwitcher::new("hyperscore", true)
            .switch_peptides(&mut doc.peptide_identifications)
            .unwrap(),
        6
    );
    keep_n_best_peptide_hits(&mut doc.peptide_identifications, 1).unwrap();
    assert_eq!(
        doc.peptide_identifications[0].hits[0].sequence.as_str(),
        "PEPTIDEA"
    );
    assert_eq!(
        doc.peptide_identifications[0].hits[0].metadata["Mascot"]
            .as_f64()
            .unwrap(),
        1.0
    );

    // The source legacy cumulative D/T counts at target thresholds 10,9,7
    // are 0/1, 0/2, 1/3. Their q-values are therefore 0,0,1/3. target+decoy
    // counts as target. Default application removes decoy hits but keeps IDs.
    FalseDiscoveryRate::default()
        .apply_peptides(&mut doc.peptide_identifications, false)
        .unwrap();
    assert_eq!(doc.peptide_identifications[0].hits[0].score, 0.0);
    assert_eq!(doc.peptide_identifications[1].hits[0].score, 0.0);
    assert_eq!(doc.peptide_identifications[3].hits[0].score, 1.0 / 3.0);
    assert!(doc.peptide_identifications[2].hits.is_empty());
    assert!(doc.peptide_identifications[4].hits.is_empty());
    assert!(
        doc.peptide_identifications
            .iter()
            .all(|id| id.score_type == "q-value" && !id.higher_score_better)
    );
    filter_peptides_by_score(&mut doc.peptide_identifications, 0.1).unwrap();
    remove_empty_peptide_identifications(&mut doc.peptide_identifications).unwrap();
    assert_eq!(doc.peptide_identifications.len(), 2);

    remove_decoy_protein_hits(&mut doc.protein_identifications).unwrap();
    remove_dangling_protein_references(
        &mut doc.peptide_identifications,
        &doc.protein_identifications,
        true,
    )
    .unwrap();
    remove_unreferenced_proteins(
        &mut doc.protein_identifications,
        &doc.peptide_identifications,
    )
    .unwrap();
    let protein = &mut doc.protein_identifications[0];
    assert!(!update_protein_groups(&mut protein.protein_groups, &protein.hits).unwrap());
    protein
        .compute_coverage(&doc.peptide_identifications)
        .unwrap();
    assert_eq!(protein.hits.len(), 1);
    assert_eq!(protein.hits[0].accession, "A");
    assert_eq!(protein.hits[0].coverage, Some(100.0));
    assert_eq!(protein.protein_groups.len(), 1);
    assert_eq!(protein.protein_groups[0].accessions, ["A"]);
    assert_eq!(protein.protein_groups[0].probability, 0.8);

    let mut xml = Vec::new();
    idxml::write(&mut xml, &doc).unwrap();
    let reread = idxml::read(Cursor::new(xml)).unwrap();
    assert_eq!(reread, doc); // Transport PH/SP IDs can change; all native records must survive.
    assert_eq!(reread.document_id, "counting-example");
    let protein = &reread.protein_identifications[0];
    assert_eq!(protein.identifier, run_id);
    assert_eq!(
        protein.metadata["operator_note"].as_str().unwrap(),
        "counted & retained"
    );
    assert_eq!(
        protein.search_parameters.metadata["database_note"]
            .as_str()
            .unwrap(),
        "independent fixture"
    );
    assert_eq!(
        protein.search_parameters.fixed_modifications,
        ["Carbamidomethyl (C)"]
    );
    let peptides = &reread.peptide_identifications;
    assert_eq!(
        peptides
            .iter()
            .map(|id| id.spectrum_reference())
            .collect::<Vec<_>>(),
        ["scan=1", "scan=2"]
    );
    assert_eq!(
        peptides[0].metadata["sample"].as_str().unwrap(),
        "synthetic digest"
    );
    assert_eq!(
        peptides[0].hits[0].metadata["annotation_note"]
            .as_str()
            .unwrap(),
        "shared evidence"
    );
    // Cleanup does not recompute target/decoy annotations. Mixed status remains
    // as recorded even though its decoy evidence has now been removed.
    assert_eq!(
        peptides[0].hits[0].target_decoy_type().unwrap(),
        TargetDecoyType::TargetAndDecoy
    );
    for (id, original_main, original_hyper) in [(&peptides[0], 1.0, 10.0), (&peptides[1], 2.0, 9.0)]
    {
        assert_eq!(id.identifier, run_id);
        assert_eq!(id.score_type, "q-value");
        assert!(!id.higher_score_better);
        let hit = &id.hits[0];
        assert_eq!(hit.score, 0.0);
        assert_eq!(hit.metadata["Mascot"].as_f64().unwrap(), original_main);
        assert_eq!(
            hit.metadata["hyperscore_score"].as_f64().unwrap(),
            original_hyper
        );
        assert_eq!(hit.evidences.len(), 1);
        assert_eq!(hit.evidences[0].protein_accession, "A");
    }
    assert_eq!(peptides[0].hits[0].evidences[0].positions().unwrap(), 0..=7);
    let modified = &peptides[1].hits[0];
    assert_eq!(modified.sequence.to_string(), "AC(Carbamidomethyl)DMK");
    assert_eq!(modified.evidences[0].positions().unwrap(), 8..=12);
    assert_eq!(modified.evidences[0].aa_after, FlankingResidue::CTerminus);
    assert_eq!(
        modified.metadata["supporting_scans"]
            .as_integer_list()
            .unwrap(),
        [2, 12]
    );
    assert_eq!(peptides[1].rt, Some(20.0));
    assert_eq!(peptides[1].mz, Some(312.627635));
}
