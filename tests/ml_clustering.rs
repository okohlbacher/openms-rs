// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Section-by-section ports of the hierarchical-clustering class tests at
//! source revision bc9cc12514c768385ce121d6ca4bb710fe1983c4
//! (`src/tests/class_tests/openms/source/`): `DistanceMatrix_test.cpp`,
//! `BinaryTreeNode_test.cpp`, `ClusterFunctor_test.cpp`,
//! `SingleLinkage_test.cpp`, `CompleteLinkage_test.cpp`,
//! `AverageLinkage_test.cpp`, `ClusterAnalyzer_test.cpp`,
//! `ClusterHierarchical_test.cpp` and `EuclideanSimilarity_test.cpp`.
//!
//! Evidence tier 3: every expected value is a transcribed class-test literal
//! (or, where marked, an invariant derived from the source's documented
//! contract). No expected value was derived from this port's output.

use openms::comparison::BinnedSharedPeakCount;
use openms::data_structures::binary_tree_node::BinaryTreeNode;
use openms::data_structures::distance_matrix::DistanceMatrix;
use openms::format::dta;
use openms::format::file_info::text_format::to_str_f32;
use openms::math::ml::clustering::average_linkage::AverageLinkage;
use openms::math::ml::clustering::cluster_analyzer::{ClusterAnalyzer, compare_binary_tree_node};
use openms::math::ml::clustering::cluster_functor::{ClusterFunctor, DEFAULT_CLUSTER_THRESHOLD};
use openms::math::ml::clustering::cluster_hierarchical::ClusterHierarchical;
use openms::math::ml::clustering::complete_linkage::CompleteLinkage;
use openms::math::ml::clustering::euclidean_similarity::EuclideanSimilarity;
use openms::math::ml::clustering::single_linkage::SingleLinkage;
use openms::{Error, MSSpectrum, Peak1D};

/// `TEST_REAL_SIMILAR` under `TOLERANCE_ABSOLUTE(abs)` with the default
/// relative tolerance 1e-5: passes when either bound holds.
fn real_similar(actual: f64, expected: f64, abs: f64) {
    let absdiff = (actual - expected).abs();
    let ratio = if actual.abs() > expected.abs() {
        actual / expected
    } else {
        expected / actual
    };
    assert!(
        absdiff <= abs || (1.0 - ratio).abs() <= 1e-5,
        "{actual} !~ {expected} (absdiff {absdiff})"
    );
}

/// The 6 x 6 matrix every linkage and analyzer class test builds with
/// `DistanceMatrix<float> matrix(6,666)` and fifteen `setValue` calls.
fn class_test_matrix(five_zero: f32) -> DistanceMatrix<f32> {
    let mut matrix = DistanceMatrix::with_value(6, 666.0_f32).unwrap();
    let cells = [
        (1, 0, 0.5_f32),
        (2, 0, 0.8),
        (2, 1, 0.3),
        (3, 0, 0.6),
        (3, 1, 0.8),
        (3, 2, 0.8),
        (4, 0, 0.8),
        (4, 1, 0.8),
        (4, 2, 0.8),
        (4, 3, 0.4),
        (5, 0, five_zero),
        (5, 1, 0.8),
        (5, 2, 0.8),
        (5, 3, 0.8),
        (5, 4, 0.8),
    ];
    for (i, j, value) in cells {
        matrix.set(i, j, value).unwrap();
    }
    matrix
}

/// The tree `SingleLinkage_test.cpp` expects, which `ClusterAnalyzer_test.cpp`
/// analyzes.
fn single_linkage_tree() -> Vec<BinaryTreeNode> {
    vec![
        BinaryTreeNode::new(1, 2, 0.3),
        BinaryTreeNode::new(3, 4, 0.4),
        BinaryTreeNode::new(0, 1, 0.5),
        BinaryTreeNode::new(0, 3, 0.6),
        BinaryTreeNode::new(0, 5, 0.7),
    ]
}

fn assert_tree(result: &[BinaryTreeNode], expected: &[BinaryTreeNode]) {
    assert_eq!(result.len(), expected.len(), "{result:?}");
    for (r, e) in result.iter().zip(expected) {
        assert_eq!(r.left_child, e.left_child, "{result:?}");
        assert_eq!(r.right_child, e.right_child, "{result:?}");
        real_similar(f64::from(r.distance), f64::from(e.distance), 0.0001);
    }
}

// ---------------------------------------------------------------- DistanceMatrix

#[test]
fn distance_matrix_class_test() {
    // DistanceMatrix()
    let empty = DistanceMatrix::<f64>::new();
    assert_eq!(empty.dimension_size(), 0);

    // DistanceMatrix(SizeType dimensionsize, Value value=Value())
    let mut dm = DistanceMatrix::with_value(8, 1.0_f64).unwrap();
    assert_eq!(dm.dimension_size(), 8);
    assert_eq!(dm.get(6, 7).unwrap(), 1.0);

    // DistanceMatrix(const DistanceMatrix &source)
    let mut dm2 = dm.clone();
    assert_eq!(dm2.dimension_size(), 8);
    assert_eq!(dm2.get(2, 3).unwrap(), 1.0);

    // void resize(SizeType dimensionsize, Value value=Value())
    dm2.resize(5, 0.0).unwrap();
    assert_eq!(dm2.dimension_size(), 5);
    // SizeType dimensionsize() const
    assert_eq!(dm2.dimension_size(), 5);

    // void setValue(SizeType i, SizeType j, ValueType value)
    for (i, j, v) in [
        (0, 1, 10.0),
        (0, 2, 9.0),
        (0, 3, 8.0),
        (0, 4, 7.0),
        (1, 2, 6.0),
        (1, 3, 5.0),
        (1, 4, 4.0),
        (2, 3, 3.0),
        (2, 4, 2.0),
        (3, 4, 0.5),
    ] {
        dm.set(i, j, v).unwrap();
    }
    assert_eq!(dm.get(0, 1).unwrap(), 10.0);
    let (a, b) = dm.min_element_coordinates().unwrap();
    assert_eq!(dm.get(a, b).unwrap(), 0.5);
    dm.set(3, 4, 1.0).unwrap();
    let (a, b) = dm.min_element_coordinates().unwrap();
    assert_eq!(dm.get(a, b).unwrap(), 1.0);
    // Regression (issue #9488, DATA-24): a smaller value in a cell sharing one
    // index with the current minimum cell must update the cached minimum.
    dm.set(2, 4, 0.1).unwrap();
    let (a, b) = dm.min_element_coordinates().unwrap();
    real_similar(dm.get(a, b).unwrap(), 0.1, 1e-5);
    real_similar(dm.get(2, 4).unwrap(), 0.1, 1e-5);
    dm.set(2, 4, 2.0).unwrap();
    let (a, b) = dm.min_element_coordinates().unwrap();
    assert_eq!(dm.get(a, b).unwrap(), 1.0);
    assert_eq!(dm.get(2, 4).unwrap(), 2.0);

    // getValue, const and non-const
    let expected = [
        (0, 1, 10.0),
        (0, 2, 9.0),
        (0, 3, 8.0),
        (0, 4, 7.0),
        (1, 2, 6.0),
        (1, 3, 5.0),
        (1, 4, 4.0),
        (2, 3, 3.0),
        (2, 4, 2.0),
        (3, 4, 1.0),
    ];
    for (i, j, v) in expected {
        assert_eq!(dm.get(i, j).unwrap(), v);
    }

    // void clear()
    dm2.clear();
    assert_eq!(dm2.dimension_size(), 0);

    // void setValueQuick(SizeType i, SizeType j, ValueType value)
    for (k, (i, j)) in [
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (1, 2),
        (1, 3),
        (1, 4),
        (2, 3),
        (2, 4),
        (3, 4),
    ]
    .into_iter()
    .enumerate()
    {
        dm.set_quick(i, j, (k + 1) as f64).unwrap();
    }
    for (k, (i, j)) in [
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (1, 2),
        (1, 3),
        (1, 4),
        (2, 3),
        (2, 4),
        (3, 4),
    ]
    .into_iter()
    .enumerate()
    {
        assert_eq!(dm.get(i, j).unwrap(), (k + 1) as f64);
        // operator(), const and non-const, is `get`, and reads the mirror.
        assert_eq!(dm.get(j, i).unwrap(), dm.get(i, j).unwrap());
    }

    // void reduce(SizeType j)
    dm.reduce(2).unwrap();
    assert_eq!(dm.get(0, 1).unwrap(), 1.0);
    assert_eq!(dm.get(0, 2).unwrap(), 3.0);
    assert_eq!(dm.get(0, 3).unwrap(), 4.0);
    assert_eq!(dm.get(1, 2).unwrap(), 6.0);
    assert_eq!(dm.get(1, 3).unwrap(), 7.0);
    assert_eq!(dm.get(2, 3).unwrap(), 10.0);
    assert_eq!(dm.dimension_size(), 7);

    // std::pair<SizeType,SizeType> getMinElementCoordinates() const
    dm.update_min_element().unwrap();
    assert_eq!(dm.min_element_coordinates().unwrap(), (1, 0));

    // void updateMinElement()
    dm.set_quick(2, 3, 0.5).unwrap();
    dm.update_min_element().unwrap();
    assert_eq!(dm.min_element_coordinates().unwrap(), (3, 2));

    // bool operator==(DistanceMatrix< ValueType > const &rhs) const
    let dm3 = dm.clone();
    assert!(dm == dm3);
    let dm_small = DistanceMatrix::with_value(3, 1.0_f64).unwrap();
    assert!(dm != dm_small);
    assert!(dm_small != dm);
}

#[test]
fn distance_matrix_boundaries() {
    // Source `Exception::OutOfRange` on every accessor.
    let mut dm = DistanceMatrix::with_value(3, 1.0_f32).unwrap();
    assert!(matches!(dm.get(3, 0), Err(Error::InvalidRange(_))));
    assert!(matches!(dm.set(0, 3, 0.0), Err(Error::InvalidRange(_))));
    assert!(matches!(
        dm.set_quick(3, 3, 0.0),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(dm.reduce(3), Err(Error::InvalidRange(_))));
    // The diagonal reads zero and ignores writes.
    dm.set(1, 1, 5.0).unwrap();
    assert_eq!(dm.get(1, 1).unwrap(), 0.0);
    let mut empty = DistanceMatrix::<f32>::new();
    assert!(matches!(
        empty.update_min_element(),
        Err(Error::InvalidRange(_))
    ));
    assert!(matches!(
        empty.min_element_coordinates(),
        Err(Error::InvalidRange(_))
    ));
    // The allocation ceiling is checked before anything is allocated, and a
    // refused resize leaves the matrix unchanged.
    assert!(matches!(
        DistanceMatrix::with_value(1 << 20, 0.0_f32),
        Err(Error::InvalidValue(_))
    ));
    assert!(dm.resize(1 << 20, 0.0).is_err());
    assert_eq!(dm.dimension_size(), 3);
    // A set after reduce, whose tracked minimum row no longer exists, is
    // defined here (the source reads freed memory).
    let mut dm = DistanceMatrix::with_value(3, 1.0_f32).unwrap();
    dm.set(2, 1, 0.1).unwrap();
    dm.reduce(2).unwrap();
    dm.set(1, 0, 0.7).unwrap();
    assert_eq!(dm.min_element_coordinates().unwrap(), (1, 0));
    assert_eq!(dm.get(1, 0).unwrap(), 0.7);
}

#[test]
fn distance_matrix_display_matches_source_stream_layout() {
    // Derived from the source `operator<<`: precision 6, width 8 consumed by
    // the first value only, a tab after every value, a newline per row.
    let mut dm = DistanceMatrix::with_value(2, 0.5_f32).unwrap();
    dm.set(1, 0, 0.25).unwrap();
    assert_eq!(dm.to_string(), "       0\t0.25\t\n0.25\t0\t\n");
}

// ---------------------------------------------------------------- BinaryTreeNode

#[test]
fn binary_tree_node_class_test() {
    // The source class test's sections are all `// TODO`; this checks the
    // constructor stores its arguments as given.
    let node = BinaryTreeNode::new(4, 2, 0.25);
    assert_eq!(
        (node.left_child, node.right_child, node.distance),
        (4, 2, 0.25)
    );
    let copy = node;
    assert_eq!(copy, node);
}

// ---------------------------------------------------------------- linkage functors

#[test]
fn single_linkage_class_test() {
    let mut matrix = class_test_matrix(0.7);
    let result = SingleLinkage::new()
        .cluster(&mut matrix, DEFAULT_CLUSTER_THRESHOLD)
        .unwrap();
    assert_tree(&result, &single_linkage_tree());
}

#[test]
fn complete_linkage_class_test() {
    let mut matrix = class_test_matrix(0.7);
    let mut matrix2 = matrix.clone();
    let result = CompleteLinkage::new().cluster(&mut matrix, 1.0).unwrap();
    let mut tree = vec![
        BinaryTreeNode::new(1, 2, 0.3),
        BinaryTreeNode::new(3, 4, 0.4),
        BinaryTreeNode::new(0, 5, 0.7),
        BinaryTreeNode::new(0, 1, 0.8),
        BinaryTreeNode::new(0, 3, 0.8),
    ];
    assert_tree(&result, &tree);

    tree.truncate(2);
    tree.push(BinaryTreeNode::new(0, 1, -1.0));
    tree.push(BinaryTreeNode::new(0, 3, -1.0));
    tree.push(BinaryTreeNode::new(0, 5, -1.0));
    let result = CompleteLinkage::new().cluster(&mut matrix2, 0.7).unwrap();
    assert_tree(&result, &tree);
}

#[test]
fn average_linkage_class_test() {
    // `0.7000001f`: the class test's "minimal adjustment for gcc 4 with -o2".
    let mut matrix = class_test_matrix(0.700_000_1);
    let mut matrix2 = matrix.clone();
    let result = AverageLinkage::new().cluster(&mut matrix, 1.0).unwrap();
    let mut tree = vec![
        BinaryTreeNode::new(1, 2, 0.3),
        BinaryTreeNode::new(3, 4, 0.4),
        BinaryTreeNode::new(0, 1, 0.65),
        BinaryTreeNode::new(0, 3, 0.766_667),
        BinaryTreeNode::new(0, 5, 0.78),
    ];
    assert_tree(&result, &tree);

    tree.truncate(3);
    tree.push(BinaryTreeNode::new(0, 3, -1.0));
    tree.push(BinaryTreeNode::new(0, 5, -1.0));
    let result = AverageLinkage::new().cluster(&mut matrix2, 0.7).unwrap();
    assert_tree(&result, &tree);
}

#[test]
fn cluster_functors_refuse_insufficient_input() {
    // Source `ClusterFunctor::InsufficientInput` for fewer than two elements.
    let functors: [&dyn ClusterFunctor; 3] = [
        &SingleLinkage::new(),
        &CompleteLinkage::new(),
        &AverageLinkage::new(),
    ];
    for functor in functors {
        for n in [0, 1] {
            let mut matrix = DistanceMatrix::with_value(n, 0.5_f32).unwrap();
            let error = functor.cluster(&mut matrix, 1.0).unwrap_err();
            assert!(matches!(&error, Error::InvalidValue(m) if m.contains("InsufficentInput")));
        }
    }
    // Source `Exception::NotImplemented` for a single-linkage threshold below 1.
    let mut matrix = class_test_matrix(0.7);
    assert!(matches!(
        SingleLinkage::new().cluster(&mut matrix, 0.5),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn every_linkage_tree_has_n_minus_one_nodes() {
    // Invariant of the source's documented contract: a clustering of n
    // elements has n - 1 steps, placeholders included, whatever the threshold.
    for threshold in [0.0_f32, 0.35, 0.55, 1.0] {
        let mut matrix = class_test_matrix(0.7);
        assert_eq!(
            CompleteLinkage::new()
                .cluster(&mut matrix, threshold)
                .unwrap()
                .len(),
            5
        );
        let mut matrix = class_test_matrix(0.7);
        assert_eq!(
            AverageLinkage::new()
                .cluster(&mut matrix, threshold)
                .unwrap()
                .len(),
            5
        );
    }
}

// ---------------------------------------------------------------- ClusterAnalyzer

#[test]
fn cluster_analyzer_average_silhouette_width() {
    let matrix = class_test_matrix(0.7);
    let asw = [0.170_833_f32, 0.309_722, 0.306_412, 0.125_744, 0.0];
    let result = ClusterAnalyzer::new()
        .average_silhouette_width(&single_linkage_tree(), &matrix)
        .unwrap();
    assert_eq!(result.len(), asw.len());
    for (r, e) in result.iter().zip(asw) {
        real_similar(f64::from(*r), f64::from(e), 0.001);
    }
}

#[test]
fn cluster_analyzer_dunn_indices() {
    let matrix = class_test_matrix(0.7);
    let di = [0.4_f32 / 0.3, 0.5 / 0.4, 0.6 / 0.8, 0.7 / 0.8, 0.0];
    for single in [false, true] {
        let result = ClusterAnalyzer::new()
            .dunn_indices(&single_linkage_tree(), &matrix, single)
            .unwrap();
        assert_eq!(result.len(), di.len());
        for (r, e) in result.iter().zip(di) {
            real_similar(f64::from(*r), f64::from(e), 0.001);
        }
    }
}

#[test]
fn cluster_analyzer_cohesion() {
    let matrix = class_test_matrix(0.7);
    let analyzer = ClusterAnalyzer::new();
    let clusters = vec![vec![0, 1, 2], vec![3, 4], vec![5]];
    let result = analyzer.cohesion(&clusters, &matrix).unwrap();
    for (r, e) in result.iter().zip([0.533_f32, 0.4, 0.7]) {
        real_similar(f64::from(*r), f64::from(e), 0.001);
    }
    assert_eq!(result.len(), 3);
    let clusters = vec![vec![0, 1, 2, 3], vec![4], vec![5]];
    let result = analyzer.cohesion(&clusters, &matrix).unwrap();
    assert_eq!(result.len(), 3);
    for (r, e) in result.iter().zip([0.633_f32, 0.7, 0.7]) {
        real_similar(f64::from(*r), f64::from(e), 0.001);
    }
}

#[test]
fn cluster_analyzer_average_population_aberration() {
    let result = ClusterAnalyzer::new()
        .average_population_aberration(3, &single_linkage_tree())
        .unwrap();
    // TEST_REAL_SIMILAR at the default absolute tolerance 1e-5.
    real_similar(f64::from(result), 2.0 / 3.0, 1e-5);
}

fn trunk() -> Vec<BinaryTreeNode> {
    vec![
        BinaryTreeNode::new(4, 6, 0.1),
        BinaryTreeNode::new(2, 3, 0.11),
        BinaryTreeNode::new(5, 8, 0.111),
        BinaryTreeNode::new(4, 7, 0.2),
        BinaryTreeNode::new(2, 9, 0.22),
        BinaryTreeNode::new(1, 10, 0.222),
        BinaryTreeNode::new(2, 11, 0.3),
        BinaryTreeNode::new(1, 5, 0.33),
        BinaryTreeNode::new(1, 12, 0.333),
        BinaryTreeNode::new(0, 1, -1.0),
        BinaryTreeNode::new(0, 2, -1.0),
        BinaryTreeNode::new(0, 4, -1.0),
    ]
}

#[test]
fn cluster_analyzer_cut_clusters() {
    let analyzer = ClusterAnalyzer::new();
    let clusters = analyzer.cut(3, &single_linkage_tree()).unwrap();
    assert_eq!(clusters, vec![vec![0, 1, 2], vec![3, 4], vec![5]]);

    let clusters = analyzer.cut(4, &trunk()).unwrap();
    assert_eq!(
        clusters,
        vec![
            vec![0],
            vec![1, 5, 8, 10, 12],
            vec![2, 3, 9, 11],
            vec![4, 6, 7]
        ]
    );
}

#[test]
fn cluster_analyzer_cut_subtrees() {
    let ts = ClusterAnalyzer::new().cut_subtrees(4, &trunk()).unwrap();
    let c_ts = vec![
        vec![],
        vec![
            BinaryTreeNode::new(5, 8, 0.111),
            BinaryTreeNode::new(1, 10, 0.222),
            BinaryTreeNode::new(1, 5, 0.33),
            BinaryTreeNode::new(1, 12, 0.333),
        ],
        vec![
            BinaryTreeNode::new(2, 3, 0.11),
            BinaryTreeNode::new(2, 9, 0.22),
            BinaryTreeNode::new(2, 11, 0.3),
        ],
        vec![
            BinaryTreeNode::new(4, 6, 0.1),
            BinaryTreeNode::new(4, 7, 0.2),
        ],
    ];
    // TEST_EQUAL on every field, the distance included.
    assert_eq!(ts, c_ts);
}

#[test]
fn cluster_analyzer_newick_tree() {
    let analyzer = ClusterAnalyzer::new();
    let tree = single_linkage_tree();
    assert_eq!(
        analyzer.newick_tree(&tree, None).unwrap(),
        "( ( ( 0 , ( 1 , 2 ) ) , ( 3 , 4 ) ) , 5 )"
    );
    assert_eq!(
        analyzer.newick_tree(&tree, Some(&to_str_f32)).unwrap(),
        "( ( ( 0:0.5 , ( 1:0.3 , 2:0.3 ):0.5 ):0.6 , ( 3:0.4 , 4:0.4 ):0.6 ):0.7 , 5:0.7 )"
    );
}

#[test]
fn cluster_analyzer_refuses_invalid_partitions() {
    let analyzer = ClusterAnalyzer::new();
    let tree = single_linkage_tree();
    let matrix = class_test_matrix(0.7);
    // Source `Exception::InvalidParameter` in each case.
    assert!(analyzer.average_silhouette_width(&[], &matrix).is_err());
    assert!(analyzer.dunn_indices(&[], &matrix, false).is_err());
    assert!(analyzer.cut(0, &tree).is_err());
    assert!(analyzer.cut(7, &tree).is_err());
    assert_eq!(analyzer.cut(6, &tree).unwrap().len(), 6);
    assert!(analyzer.cut_subtrees(0, &tree).is_err());
    assert!(analyzer.cut_subtrees(6, &tree).is_err());
    assert!(analyzer.average_population_aberration(0, &tree).is_err());
    assert!(analyzer.average_population_aberration(6, &tree).is_err());
    assert!(analyzer.cohesion(&[], &matrix).is_err());
    assert!(analyzer.cohesion(&vec![vec![0]; 7], &matrix).is_err());
    // Out-of-range element indices are errors, not out-of-bounds reads.
    let bad = vec![
        BinaryTreeNode::new(0, 9, 0.1),
        BinaryTreeNode::new(0, 1, 0.2),
    ];
    assert!(analyzer.dunn_indices(&bad, &matrix, false).is_err());
    assert!(analyzer.average_population_aberration(1, &bad).is_err());
    assert!(analyzer.newick_tree(&[], None).is_err());
    // compareBinaryTreeNode is `x.distance < y.distance`.
    assert!(compare_binary_tree_node(&tree[0], &tree[1]));
    assert!(!compare_binary_tree_node(&tree[1], &tree[0]));
}

// ---------------------------------------------------------------- ClusterHierarchical

/// The class test's `LowlevelComparator`: similarity `1 - d` for the class-test
/// matrix, 0 on the diagonal and for any unlisted pair.
fn lowlevel_comparator(first: &usize, second: &usize) -> openms::Result<f64> {
    let (x, y) = ((*first).min(*second), (*first).max(*second));
    Ok(match (x, y) {
        (0, 1) => 1.0 - 0.5,
        (0, 2) => 1.0 - 0.8,
        (0, 3) => 1.0 - 0.6,
        (0, 4) => 1.0 - 0.8,
        (0, 5) => 1.0 - 0.7,
        (1, 2) => 1.0 - 0.3,
        (1, 3) | (1, 4) | (1, 5) => 1.0 - 0.8,
        (2, 3) | (2, 4) | (2, 5) => 1.0 - 0.8,
        (3, 4) => 1.0 - 0.4,
        (3, 5) => 1.0 - 0.8,
        (4, 5) => 1.0 - 0.8,
        (0..=4, _) => 0.0,
        _ => 666.0,
    })
}

#[test]
fn cluster_hierarchical_threshold() {
    let mut ch = ClusterHierarchical::new();
    ch.set_threshold(66.6);
    let copy = ch;
    assert_eq!(copy.threshold(), 66.6);
    let mut ch = ClusterHierarchical::default();
    assert_eq!(ch.threshold(), 1.0);
    ch.set_threshold(0.666);
    assert_eq!(ch.threshold(), 0.666);
}

#[test]
fn cluster_hierarchical_generic_cluster() {
    let d: Vec<usize> = (0..6).collect();
    let mut matrix = DistanceMatrix::new();
    let result = ClusterHierarchical::new()
        .cluster(&d, lowlevel_comparator, &SingleLinkage::new(), &mut matrix)
        .unwrap();
    assert_tree(&result, &single_linkage_tree());
    assert_eq!(matrix.dimension_size(), 6);
}

#[test]
fn cluster_hierarchical_binned_spectra() {
    use openms::analysis::spectra_clustering::cluster_binned_spectra;

    let s1: MSSpectrum =
        dta::read(include_bytes!("data/comparison_dfpianger.dta").as_slice()).unwrap();
    let mut s2 = s1.clone();
    let mut s3 = s1.clone();
    s2.peaks.pop();
    s3.peaks.pop();
    let peak = Peak1D {
        mz: 666.66,
        intensity: 999.99,
    };
    s2.peaks.push(peak);
    s2.sort_by_position().unwrap();
    s3.peaks.push(peak);
    s3.sort_by_position().unwrap();

    let data = vec![s1, s2, s3];
    let tree = vec![
        BinaryTreeNode::new(1, 2, 0.0),
        BinaryTreeNode::new(0, 1, 0.008_498_58),
    ];
    let mut matrix = DistanceMatrix::new();
    let result = cluster_binned_spectra(
        &ClusterHierarchical::new(),
        &data,
        &BinnedSharedPeakCount::new().unwrap(),
        1.5,
        2,
        0.4, // BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES
        &SingleLinkage::new(),
        &mut matrix,
    )
    .unwrap();
    assert_tree(&result, &tree);
}

// ---------------------------------------------------------------- EuclideanSimilarity

#[test]
fn euclidean_similarity_class_test() {
    let es = EuclideanSimilarity::new();
    real_similar(
        f64::from(es.similarity((2.0, 2.0), (4.0, 4.0)).unwrap()),
        1.0 - 8.0_f64.sqrt(),
        0.0001,
    );
    real_similar(
        f64::from(es.similarity((9.0, 0.1), (2.8, 2.0)).unwrap()),
        1.0 - 42.05_f64.sqrt(),
        0.0001,
    );
    real_similar(
        f64::from(es.similarity((12.0, 0.0), (2.0, 0.0)).unwrap()),
        1.0 - 100.0_f64.sqrt(),
        0.0001,
    );

    // float operator()(const std::pair<float, float>& c) const
    real_similar(
        f64::from(es.self_similarity((9.0, 0.1)).unwrap()),
        1.0,
        0.0001,
    );
    real_similar(
        f64::from(es.self_similarity((2.8, 2.0)).unwrap()),
        1.0,
        0.0001,
    );

    // void setScale(float x)
    let mut es = EuclideanSimilarity::default();
    es.set_scale(10.0);
    real_similar(
        f64::from(es.similarity((2.0, 2.0), (4.0, 4.0)).unwrap()),
        1.0 - 8.0_f64.sqrt() / 10.0,
        0.0001,
    );
    real_similar(
        f64::from(es.similarity((9.0, 0.1), (2.8, 2.0)).unwrap()),
        1.0 - 42.05_f64.sqrt() / 10.0,
        0.0001,
    );
    real_similar(
        f64::from(es.similarity((12.0, 0.0), (2.0, 0.0)).unwrap()),
        1.0 - 100.0_f64.sqrt() / 10.0,
        0.0001,
    );
    es.set_scale(233.28_f32.sqrt());
    real_similar(
        f64::from(es.similarity((0.1, 0.1), (10.9, 10.9)).unwrap()),
        0.0,
        0.0001,
    );

    // Source `Exception::DivisionByZero` for scale 0.
    es.set_scale(0.0);
    assert!(es.similarity((0.0, 0.0), (1.0, 1.0)).is_err());
}
