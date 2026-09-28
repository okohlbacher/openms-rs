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
1. [x] DistanceMatrix + BinaryTreeNode (data_structures) - committed 4dfa4de8
2. [x] math/ml/clustering: ClusterFunctor, Single/Complete/AverageLinkage, ClusterAnalyzer,
       ClusterHierarchical (generic in math; binned overload = src/analysis/spectra_clustering.rs),
       EuclideanSimilarity. tests/ml_clustering.rs (no features). Gate: dax test run pending/passing?
   - check_module_cycles.py counts `crate::x` inside doc comments too! avoid crate:: paths to format/analysis in math docs.
   - recorded new acyclic edges with --write (docs/module-cycles.json) -> mention in report.
   - newick_tree takes Option<&dyn Fn(f32)->String>; test passes format::file_info::text_format::to_str_f32.
   - Grid: extendClustersY only used by class test (commented out in MultiplexClustering.cpp:165);
     MultiplexDistance = sqrt(dx^2 + s^2 dy^2) (MultiplexClustering.cpp:184-187), not ported in Rust.
3. [x] ClusteringGrid, GridBasedCluster, GridBasedClustering (DBoundingBox = DIntervalBase + local enlarge)
4. [x] HashGrid (BTreeMap cells, insertion-order cell content; boost hash order NOT reproduced -> partial)
5. [x] GridFeature, QTCluster in src/analysis/mapmatching/ (SourceOrderedSequence emulates AASequence operator< length-first)
6. [ ] KDTree: untouched (2208-line libkdtree++; consumer KDTreeFeatureMaps; needs nth_element-faithful optimize) -> report as backlog
Tests: tests/ml_clustering.rs (21), tests/ml_grid_clustering.rs (9), tests/mapmatching_qt_cluster.rs (4); all pass on dax,
  also under +1.85.0 --no-default-features. clippy --all-features clean.
check_doc_coverage: passes but floor file docs/doc-coverage.json is LEAD's -> do not --write; report.
Oracle: ../oracle/FF2a-clustering (driver.cpp, run.sh, remote_run.sh) running on ibminode06; results/oracle.tsv
  -> copy to tests/data/ml_clustering_oracle.tsv + write tests/ml_clustering_oracle.rs replaying LCG inputs.

## Next step
1. fetch oracle results, write tests/ml_clustering_oracle.rs, compare.
2. docs/CLUSTERING_SUPPORT.md, tests/data/ml_clustering_provenance.json, ledger fragment, lane report.
