# FF2a lane notes (ML clustering)

Worktree branch: worktree-wf_cd066dab-0a7-34 (off main f0990c0).
C++ pin: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core @ bc9cc12
(use `git -C <abs path> show bc9cc12:<path>` as its own plain command; compound git commands get refused).

## Scope (all missing in Rust at start; src/math/ml.rs is an empty module root)
ML/CLUSTERING: AverageLinkage, ClusterAnalyzer, ClusterFunctor, ClusterHierarchical, ClusteringGrid,
CompleteLinkage, EuclideanSimilarity, GridBasedCluster, GridBasedClustering, HashGrid, SingleLinkage
DATASTRUCTURES: GridFeature, KDTree, QTCluster (+ BinaryTreeNode is a dependency of the linkage family)

Class tests exist for every header except KDTree.h (header-only third-party libkdtree++ wrapper).

## Dependency order
1. BinaryTreeNode (DATASTRUCTURES), ClusterFunctor, EuclideanSimilarity, HashGrid, ClusteringGrid, GridBasedCluster
2. SingleLinkage/CompleteLinkage/AverageLinkage (ClusterFunctor + BinaryTreeNode)
3. ClusterAnalyzer (BinaryTreeNode), ClusterHierarchical
4. GridBasedClustering (HashGrid/ClusteringGrid/GridBasedCluster)
5. GridFeature/QTCluster need Feature/ConsensusFeature/PeptideIdentification -> not math; place in data_structures or analysis

## Next step
1. Read C++ headers + tests; port step 1 into src/math/ml/clustering/...
