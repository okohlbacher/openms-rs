// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Retain parent positions and processing provenance when digesting RNA.
use openms::chemistry::RNaseDigestion;
use openms::identification::graph::{
    IdentificationData, MoleculeType, ParentSequence, ProcessingSoftware, ProcessingStep,
};

fn main() -> openms::Result<()> {
    let mut graph = IdentificationData::new()?;
    let mut parent = ParentSequence::new("rna_1");
    parent.molecule_type = MoleculeType::RNA;
    parent.sequence = "pAUGUCGCAG".into();
    let parent_id = graph.register_parent_sequence(parent)?;
    let software = graph.register_processing_software(ProcessingSoftware::new("RNase_T1", "1"))?;
    let step = graph.register_processing_step(ProcessingStep::new(software), None)?;
    graph.set_current_processing_step(step)?;
    RNaseDigestion::default().digest_identification_data(&mut graph)?;
    graph.calculate_coverages(true)?;
    for (_, oligo) in graph.oligos() {
        for matched in &oligo.parent_matches[&parent_id] {
            println!(
                "{}\t{}..={}\t{} | {}\t{} processing step(s)",
                oligo.sequence,
                matched.start_pos.expect("known digestion start"),
                matched.end_pos.expect("known digestion end"),
                matched.left_neighbor,
                matched.right_neighbor,
                oligo.result.steps_and_scores.len()
            );
        }
    }
    println!(
        "Parent coverage: {:.0}%",
        graph.parent(parent_id)?.coverage * 100.0
    );
    Ok(())
}
