# Lane M2 notes — feature grouping and normalisation

Prereq merged: worktree-wf_cd066dab-0a7-34 (FF2a: ML clustering, GridFeature, QTCluster) — fast-forward.

## Reading C++
- `git -C .../OpenMS4-tests/packages/core show bc9cc12:<path>` works only as a plain literal command
  (the worktree guard refuses variables/pipes around git). Working tree of core is at 291c5ba, NOT the pin;
  FeatureGroupingAlgorithmUnlabeled.cpp and KERNEL/ConsensusMap differ from the pin — always use `show bc9cc12:`.

## Existing Rust
- kernel types: src/kernel/features.rs (BaseFeature, Feature, FeatureHandle, ConsensusFeature, FeatureMap, ConsensusMap)
- src/analysis/mapmatching/{grid_feature.rs, qt_cluster.rs}; src/math/ml/clustering/* (GridBasedClustering, HashGrid)
- Param pattern: src/processing/peak_picking.rs:1853-1985

## Scope order (dependency first)
1. ConsensusMapNormalizerAlgorithm{Median,Quantile,Threshold} (self-contained, class tests)
2. FeatureDistance -> BaseGroupFinder -> StablePairFinder, QTClusterFinder, LabeledPairFinder
3. FeatureGroupingAlgorithm (+Unlabeled, QT, Labeled, KD, WNet), WNetMatcher, FeatureMapping, PipEchoAlgorithm

## Next step
1. Read headers/tests for the normalizers, port them.
