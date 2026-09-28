# FF2a lane notes (ML clustering)

Worktree branch: worktree-wf_cd066dab-0a7-34 (off main f0990c0).
C++ pin: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core @ bc9cc12
Read with ONE plain command: `git -C <abs core path> show bc9cc12:<p1> bc9cc12:<p2>` (no loops/pipes
into git, no `archive`; the worktree guard refuses them). Class tests at
src/tests/class_tests/openms/source/<Name>_test.cpp. KDTree.h has no class test.

## Findings (sources read)
- Deps not in scope but needed + unported: DATASTRUCTURES/DistanceMatrix.h, DATASTRUCTURES/BinaryTreeNode.h
  -> port into src/data_structures/ (data_structures -> param only; math -> data_structures closes no cycle).
- Module graph: math has NO out-edges today; math -> data_structures/concept/param close no cycle.
  data_structures -> kernel WOULD close a cycle (kernel->metadata->data_structures). So GridFeature
  (needs BaseFeature + AASequence + PeptideIdentification) and QTCluster go to src/analysis/mapmatching/.
- `crate::math::source_sort::source_sort_by` emulates libstdc++ std::sort (needed: SingleLinkage sorts
  pre-tree by distance with std::sort, unstable; ties observable).
- `crate::math::libstdcxx::lower_bound` for ClusteringGrid::getIndex (std::lower_bound with less_equal).
- ClusterHierarchical: generic `cluster` fits math; binned-PeakSpectrum overload needs BinnedSpectrum
  (src/kernel.rs) + BinnedSpectrumCompareFunctor (src/comparison.rs?) -> cannot live in math.
- HashGrid uses boost::unordered iteration order (not reproducible) -> BTreeMap, documented divergence.
- GridBasedClustering quirks: multiset<MinimumDistance> ordered by distance only (equal -> insertion
  order after existing); extendClustersY uses std::map<GridBasedCluster,int> keyed by centre y (equal y
  collapse -> erase(end()) UB), stale-list merge chain, only grid y-cell 1 scanned (y==ymin skipped).
- EuclideanSimilarity doc says self similarity "will yield 0" but code yields 1 -> doc defect candidate.
- Complete/AverageLinkage fill tail with dummy nodes distance -1 up to capacity n-1.

## Plan / status
1. [ ] DistanceMatrix + BinaryTreeNode (data_structures)
2. [ ] math/ml/clustering: ClusterFunctor, Single/Complete/AverageLinkage, ClusterAnalyzer,
       ClusterHierarchical, EuclideanSimilarity
3. [ ] ClusteringGrid, GridBasedCluster, GridBasedClustering (+ DBoundingBox need: check drange.rs)
4. [ ] HashGrid
5. [ ] GridFeature, QTCluster (analysis/mapmatching)
6. [ ] KDTree (crate search: kiddo/kdtree) - likely untouched/partial

## Next step
1. Read class tests for step 1-2, write Rust, gate on dax with slot FF2a.
