# Lane M1 notes (map alignment) -- delete in the final commit

## Scratch
C++ pins extracted (git archive) under
`/private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/M1/`
with symlinks `mct` (class tests), `mcpp` (C++ implementations), `mhdr` (headers), `mdata` (MapAlign* class-test data).
THE SCRATCHPAD ROOT IS SHARED WITH OTHER SESSIONS: only ever create things inside `M1/`.
The worktree guard refuses commands where a path ends in the C++-implementation dir name or that use heredocs:
use the symlinks, and write files with the Write tool.
Gate logs: `M1/gate*.log`.

## Done (committed)
1. 48b83c7 transformations: B-spline model, Akima, LOWESS auto-span (+ libstdc++ sort of anchors), Param surface
   (defaults/from_param/to_param per model, ModelConfig::from_param, TransformationDescription::fit_model_with_param,
   model_parameters, model_type(s), print_summary, weights), TransformationModelDefaults; TrafoXML b_spline+akima.
   tests/transformation_models.rs (17 tests) green on dax; tests/transformation_xml.rs b_spline assertion updated.
2. a22dad0 mapmatching: superimposer.rs (BaseSuperimposer trait + histogram/tophat/nth_element), pose_clustering_shift.rs,
   pose_clustering_affine.rs, alignment_evaluation.rs. tests/pose_clustering.rs (5), tests/map_alignment_evaluation.rs (5,
   needs feature consensusxml). docs/module-cycles.json re-recorded (new edge analysis->data_structures).

## Gaps so far
- Superimposer `dump_buckets`/`dump_pairs` debug outputs: refused (Unsupported) -> superimposers partial.
- InterpolatedModel::fit_pairs reproduces the source zero-padding defect (CPP defect candidate).
- MapAlignmentAlgorithmPoseClustering blocked on StablePairFinder (lane M2).

## Next step
1. MapAlignmentAlgorithmIdentification (src/analysis/mapmatching/alignment_identification.rs): &mut self state
   machine (reference_, better_, min_run_occur_), inputs FeatureMap/ConsensusMap/Vec<PeptideIdentification>;
   IdentificationData overload = gap unless identification::graph has pickScoreType/getBestMatchPerObservation.
   Class test mct/MapAlignmentAlgorithmIdentification_test.cpp, data mdata/*.idXML.
2. Then docs: docs/TRANSFORMATIONS_SUPPORT.md update, docs/MAP_ALIGNMENT_SUPPORT.md new, provenance json(s),
   THIRD_PARTY_CRATE_DECISIONS entry (Akima), doc coverage ratchet, clippy/doc/MSRV gates, report + ledger fragment.
