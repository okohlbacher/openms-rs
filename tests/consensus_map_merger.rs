// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

//! Every section of `ConsensusMapMergerAlgorithm_test.cpp` at core `bc9cc12`
//! (on `BSA.consensusXML`: six label-free columns, one identification run
//! each), plus invariants derived from the source `id_merge_index` handling.

use openms::analysis::id::consensus_map_merger::ConsensusMapMergerAlgorithm;
use openms::param::ParamValue;

#[test]
fn parameters_match_the_source_declarations() {
    let d = ConsensusMapMergerAlgorithm::defaults().unwrap();
    assert_eq!(
        d.value("annotate_origin").unwrap(),
        &ParamValue::String("true".into())
    );
    let off = ConsensusMapMergerAlgorithm {
        annotate_origin: false,
        ..ConsensusMapMergerAlgorithm::default()
    };
    assert_eq!(
        ConsensusMapMergerAlgorithm::from_param(&off.to_param().unwrap()).unwrap(),
        off
    );
}

#[cfg(feature = "consensusxml")]
mod bsa {
    use super::*;
    use openms::format::consensusxml;
    use openms::kernel::ConsensusMap;
    use openms::metadata::{ExperimentalDesign, SampleSection};
    use std::collections::BTreeMap;
    use std::path::Path;

    const DIR: &str = "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/";

    /// The default parse-work ceiling (50,000,000) is exhausted by this
    /// 136 KB file: each of its 42 modified sequences is charged about 1.8
    /// million units for modification lookups. The ceiling is raised here.
    fn load() -> ConsensusMap {
        let options = consensusxml::ReadOptions {
            max_work: 1_000_000_000,
            ..consensusxml::ReadOptions::default()
        };
        consensusxml::load_with_options(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/BSA.consensusXML"),
            &options,
        )
        .unwrap()
    }

    fn merge_indices(map: &ConsensusMap) -> Vec<(String, i64)> {
        map.features
            .iter()
            .flat_map(|f| f.peptide_identifications.iter())
            .chain(&map.unassigned_peptide_identifications)
            .map(|p| {
                (
                    p.identifier.clone(),
                    p.metadata
                        .get("id_merge_index")
                        .map_or(-1, |v| v.as_i64().unwrap()),
                )
            })
            .collect()
    }

    #[test]
    fn merge_all_id_runs() {
        let mut cmap = load();
        let before = cmap.protein_identifications.clone();
        ConsensusMapMergerAlgorithm::default()
            .merge_all_id_runs(&mut cmap)
            .unwrap();
        assert_eq!(cmap.protein_identifications.len(), 1);
        // Derived: the merged run lists every old run's files in order, and a
        // peptide's index is its old index plus its run's offset.
        let merged = &cmap.protein_identifications[0];
        assert_eq!(merged.identifier, "merged");
        let paths: Vec<_> = before
            .iter()
            .flat_map(|r| r.primary_ms_run_paths.iter().cloned())
            .collect();
        assert_eq!(merged.primary_ms_run_paths, paths);
        assert!(
            merge_indices(&cmap)
                .iter()
                .all(|(run, index)| run == "merged" && (0..6).contains(index))
        );
    }

    #[test]
    fn merge_proteins_across_fractions_and_replicates_no_design() {
        let mut cmap = load();
        let design = ExperimentalDesign::from_consensus_map(&cmap).unwrap();
        ConsensusMapMergerAlgorithm::default()
            .merge_proteins_across_fractions_and_replicates(&mut cmap, &design)
            .unwrap();
        // Without a sample-level design, runs count as replicates or fractions
        // and all are merged.
        assert_eq!(cmap.protein_identifications.len(), 1);
        assert_eq!(
            cmap.protein_identifications[0].primary_ms_run_paths.len(),
            6
        );
    }

    fn three_conditions() -> (ConsensusMap, ExperimentalDesign) {
        let cmap = load();
        let mut design = ExperimentalDesign::from_consensus_map(&cmap).unwrap();
        let text = |v: &[&str]| v.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>();
        let section = SampleSection::from_table(
            vec![text(&["1", "C1"]), text(&["2", "C2"]), text(&["3", "C3"])],
            BTreeMap::from([("1".into(), 0), ("2".into(), 1), ("3".into(), 2)]),
            BTreeMap::from([("Sample".into(), 0), ("Condition".into(), 1)]),
        )
        .unwrap();
        design.set_sample_section(section);
        (cmap, design)
    }

    #[test]
    fn merge_proteins_across_fractions_and_replicates() {
        let (mut cmap, design) = three_conditions();
        ConsensusMapMergerAlgorithm::default()
            .merge_proteins_across_fractions_and_replicates(&mut cmap, &design)
            .unwrap();
        assert_eq!(cmap.protein_identifications.len(), 3);
        for (run, sample) in cmap.protein_identifications.iter().zip(1..) {
            assert_eq!(
                run.primary_ms_run_paths,
                vec![
                    format!("{DIR}BSA{sample}_F1.mzML"),
                    format!("{DIR}BSA{sample}_F2.mzML")
                ]
            );
        }
    }

    #[test]
    fn id_merge_index_offsets_native_and_source() {
        // Derived from the source mergeProteinIDRuns: the second fraction's
        // run sits at offset 1 in its condition run. The source looks the
        // offset up under the new run identifier and always finds 0.
        let (mut native, design) = three_conditions();
        let old: Vec<_> = native
            .features
            .iter()
            .flat_map(|f| f.peptide_identifications.iter())
            .chain(&native.unassigned_peptide_identifications)
            .map(|p| p.identifier.clone())
            .collect();
        let second_fraction_runs: Vec<_> = native
            .protein_identifications
            .iter()
            .filter(|r| {
                r.primary_ms_run_paths
                    .iter()
                    .any(|p| p.ends_with("_F2.mzML"))
            })
            .map(|r| r.identifier.clone())
            .collect();
        let mut source = native.clone();
        ConsensusMapMergerAlgorithm::default()
            .merge_proteins_across_fractions_and_replicates(&mut native, &design)
            .unwrap();
        ConsensusMapMergerAlgorithm {
            source_merge_index_offsets: true,
            ..ConsensusMapMergerAlgorithm::default()
        }
        .merge_proteins_across_fractions_and_replicates(&mut source, &design)
        .unwrap();
        let native = merge_indices(&native);
        let source = merge_indices(&source);
        assert!(!old.is_empty());
        let mut seen_second = false;
        for ((run, (_, n)), (_, s)) in old.iter().zip(&native).zip(&source) {
            assert_eq!(*s, 0);
            if second_fraction_runs.contains(run) {
                seen_second = true;
                assert_eq!(*n, 1);
            } else {
                assert_eq!(*n, 0);
            }
        }
        assert!(seen_second);
    }
}
