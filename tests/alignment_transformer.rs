// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::analysis::alignment_transformer::MapAlignmentTransformer;
use openms::analysis::transformations::*;
use openms::identification::PeptideIdentification;
use openms::kernel::features::{
    BaseFeature, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap,
};
use openms::kernel::geometry::{ConvexHull2D, Point2D};
use openms::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D,
};

fn linear(slope: f64, intercept: f64) -> TransformationDescription {
    let mut description = TransformationDescription::default();
    description
        .fit_model(ModelConfig::Linear(LinearOptions {
            coefficients: Some(LinearCoefficients { slope, intercept }),
            ..Default::default()
        }))
        .unwrap();
    description
}
fn peptide(rt: Option<f64>) -> PeptideIdentification {
    PeptideIdentification {
        rt,
        ..Default::default()
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} != {b}");
}

#[test]
fn upstream_experiment_retention_times_and_original_value_preservation() {
    // MapAlignmentTransformer_test.cpp anchors (0,1), (1,3).
    let transformation = linear(2., 1.);
    let times = [11.1, 11.5, 12.2, 12.5];
    let expected = [23.2, 24., 25.4, 26.];
    let mut experiment = MSExperiment {
        spectra: times
            .iter()
            .enumerate()
            .map(|(i, &rt)| MSSpectrum {
                rt,
                ms_level: 1 + (i as u32 % 2),
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    };
    MapAlignmentTransformer::default()
        .transform_experiment(&mut experiment, &transformation)
        .unwrap();
    for (spectrum, expected) in experiment.spectra.iter().zip(expected) {
        close(spectrum.rt, expected);
        assert!(!spectrum.metadata.contains_key("original_RT"));
    }
    let transformer = MapAlignmentTransformer {
        store_original_rt: true,
        ..Default::default()
    };
    transformer
        .transform_experiment(&mut experiment, &transformation)
        .unwrap();
    transformer
        .transform_experiment(&mut experiment, &transformation)
        .unwrap();
    for (spectrum, expected) in experiment.spectra.iter().zip(expected) {
        close(spectrum.metadata["original_RT"].as_f64().unwrap(), expected);
    }
    close(experiment.ranges(0).unwrap().rt.unwrap().min, 95.8);
}

#[test]
fn upstream_feature_consensus_and_peptide_values() {
    let transformation = linear(2., 1.);
    let times = [11.1, 11.5, 12.2, 12.5];
    let expected = [23.2, 24., 25.4, 26.];
    let mut features =
        FeatureMap::from_features(times.iter().map(|&rt| Feature::new(rt, 100., 1.)).collect());
    let mut consensus = ConsensusMap::from_features(
        times
            .iter()
            .map(|&rt| ConsensusFeature::from(BaseFeature::new(rt, 100., 1.)))
            .collect(),
    );
    let mut ids: Vec<_> = times.iter().map(|&rt| peptide(Some(rt))).collect();
    let transformer = MapAlignmentTransformer::default();
    transformer
        .transform_feature_map(&mut features, &transformation)
        .unwrap();
    transformer
        .transform_consensus_map(&mut consensus, &transformation)
        .unwrap();
    transformer
        .transform_peptide_identifications(&mut ids, &transformation)
        .unwrap();
    for i in 0..4 {
        close(features.features[i].rt, expected[i]);
        close(consensus.features[i].rt, expected[i]);
        close(ids[i].rt.unwrap(), expected[i]);
    }
    let transformer = MapAlignmentTransformer {
        store_original_rt: true,
        ..Default::default()
    };
    for _ in 0..2 {
        transformer
            .transform_feature_map(&mut features, &transformation)
            .unwrap();
        transformer
            .transform_consensus_map(&mut consensus, &transformation)
            .unwrap();
        transformer
            .transform_peptide_identifications(&mut ids, &transformation)
            .unwrap();
    }
    for i in 0..4 {
        close(
            features.features[i].metadata["original_RT"]
                .as_f64()
                .unwrap(),
            expected[i],
        );
        close(
            consensus.features[i].metadata["original_RT"]
                .as_f64()
                .unwrap(),
            expected[i],
        );
        close(
            ids[i].metadata["original_RT"].as_f64().unwrap(),
            expected[i],
        );
    }
}

#[test]
fn chromatogram_arrays_order_ranges_and_optional_spectrum_identifications() {
    let mut spectrum = MSSpectrum::from_peaks(vec![Peak1D::new(100., 2.)]);
    spectrum.rt = 2.;
    spectrum.peptide_identifications.push(peptide(Some(3.)));
    let mut chromatogram = MSChromatogram::from_peaks(vec![
        ChromatogramPeak::new(1., 4.),
        ChromatogramPeak::new(3., 2.),
    ]);
    chromatogram.string_data_arrays.push(DataArray::new(
        "sample",
        vec!["first".into(), "second".into()],
    ));
    let mut experiment = MSExperiment {
        spectra: vec![spectrum],
        chromatograms: vec![chromatogram],
        ..Default::default()
    };
    let transformer = MapAlignmentTransformer {
        store_original_rt: true,
        ..Default::default()
    };
    transformer
        .transform_experiment(&mut experiment, &linear(-2., 10.))
        .unwrap();
    assert_eq!(
        experiment.spectra[0].peptide_identifications[0].rt,
        Some(3.)
    );
    assert_eq!(
        experiment.chromatograms[0]
            .peaks
            .iter()
            .map(|p| p.rt)
            .collect::<Vec<_>>(),
        [8., 4.]
    );
    assert!(!experiment.chromatograms[0].is_sorted());
    assert_eq!(
        experiment.chromatograms[0].metadata["original_rt"]
            .as_float_list()
            .unwrap(),
        &[1., 3.]
    );
    assert_eq!(
        experiment.chromatograms[0].string_data_arrays[0].data,
        ["first", "second"]
    );
    assert_eq!(
        experiment.chromatograms[0]
            .ranges()
            .unwrap()
            .rt
            .unwrap()
            .min,
        4.
    );
    let transformer = MapAlignmentTransformer {
        transform_spectrum_identifications: true,
        ..transformer
    };
    transformer
        .transform_experiment(&mut experiment, &linear(2., 1.))
        .unwrap();
    assert_eq!(
        experiment.spectra[0].peptide_identifications[0].rt,
        Some(7.)
    );
    assert_eq!(
        experiment.spectra[0].peptide_identifications[0].metadata["original_RT"]
            .as_f64()
            .unwrap(),
        3.
    );
    assert_eq!(
        experiment.chromatograms[0].metadata["original_rt"]
            .as_float_list()
            .unwrap(),
        &[1., 3.]
    );
    assert_eq!(experiment.spectra[0].peaks, [Peak1D::new(100., 2.)]);
}

#[test]
fn feature_subordinates_hulls_assigned_and_unassigned_ids() {
    let mut feature = Feature::new(2., 500., 100.);
    feature.peptide_identifications = vec![peptide(Some(2.5)), peptide(None)];
    feature.subordinates.push(Feature::new(1., 499., 50.));
    feature.subordinates[0]
        .peptide_identifications
        .push(peptide(Some(1.5)));
    feature.convex_hulls.push(
        ConvexHull2D::from_points(&[
            Point2D::new(1., 499.),
            Point2D::new(1., 501.),
            Point2D::new(3., 499.),
            Point2D::new(3., 501.),
        ])
        .unwrap(),
    );
    let before = feature.convex_hulls[0].hull_points();
    let mut map = FeatureMap::from_features(vec![feature]);
    map.unassigned_peptide_identifications
        .push(peptide(Some(4.)));
    MapAlignmentTransformer {
        store_original_rt: true,
        ..Default::default()
    }
    .transform_feature_map(&mut map, &linear(2., 1.))
    .unwrap();
    let feature = &map.features[0];
    assert_eq!(feature.rt, 5.);
    assert_eq!(feature.subordinates[0].rt, 3.);
    assert_eq!(feature.peptide_identifications[0].rt, Some(6.));
    assert_eq!(
        feature.subordinates[0].peptide_identifications[0].rt,
        Some(4.)
    );
    assert_eq!(map.unassigned_peptide_identifications[0].rt, Some(9.));
    assert!(feature.peptide_identifications[1].metadata.is_empty());
    assert_eq!(feature.convex_hulls[0].scan_count(), 0);
    for (actual, old) in feature.convex_hulls[0].hull_points().iter().zip(before) {
        assert_eq!(actual.rt, 2. * old.rt + 1.);
        assert_eq!(actual.mz, old.mz);
    }
    assert_eq!(map.ranges().unwrap().rt.unwrap().min, 3.);
    assert_eq!(map.ranges().unwrap().rt.unwrap().max, 7.);
    assert_eq!(feature.mz, 500.);
    assert_eq!(feature.intensity, 100.);
}

#[test]
fn consensus_handles_transform_without_recomputing_centroid_or_identities() {
    let mut feature = ConsensusFeature::from(BaseFeature::new(5., 100., 10.));
    feature.peptide_identifications.push(peptide(Some(4.)));
    feature
        .insert(FeatureHandle {
            map_index: 7,
            unique_id: 42,
            rt: 1.,
            mz: 99.,
            intensity: 2.,
            ..Default::default()
        })
        .unwrap();
    feature
        .insert(FeatureHandle {
            map_index: 2,
            unique_id: 84,
            rt: 3.,
            mz: 101.,
            intensity: 4.,
            ..Default::default()
        })
        .unwrap();
    let mut map = ConsensusMap::from_features(vec![feature]);
    map.unassigned_peptide_identifications
        .push(peptide(Some(6.)));
    MapAlignmentTransformer {
        store_original_rt: true,
        ..Default::default()
    }
    .transform_consensus_map(&mut map, &linear(2., 1.))
    .unwrap();
    let feature = &map.features[0];
    assert_eq!(feature.rt, 11.);
    assert_eq!(
        feature
            .handles()
            .iter()
            .map(|h| (h.map_index, h.unique_id, h.rt))
            .collect::<Vec<_>>(),
        [(2, 84, 7.), (7, 42, 3.)]
    );
    assert_eq!(feature.peptide_identifications[0].rt, Some(9.));
    assert_eq!(map.unassigned_peptide_identifications[0].rt, Some(13.));
    assert_eq!(feature.handle_ranges().rt.unwrap().min, 3.);
}

fn logarithm() -> TransformationDescription {
    let mut transformation = TransformationDescription::default();
    transformation
        .fit_model(ModelConfig::Linear(LinearOptions {
            x_weight: CoordinateWeight {
                function: WeightFunction::Log,
                ..Default::default()
            },
            coefficients: Some(LinearCoefficients {
                slope: 1.,
                intercept: 0.,
            }),
            ..Default::default()
        }))
        .unwrap();
    transformation
}
#[test]
fn every_container_is_unchanged_after_late_transform_failure() {
    let transform = logarithm();
    let transformer = MapAlignmentTransformer {
        store_original_rt: true,
        ..Default::default()
    };
    let mut experiment = MSExperiment {
        spectra: vec![MSSpectrum {
            rt: 2.,
            ..Default::default()
        }],
        chromatograms: vec![MSChromatogram::from_peaks(vec![
            ChromatogramPeak::new(2., 1.),
            ChromatogramPeak::new(0., 1.),
        ])],
        ..Default::default()
    };
    let before = experiment.clone();
    assert!(
        transformer
            .transform_experiment(&mut experiment, &transform)
            .is_err()
    );
    assert_eq!(experiment, before);
    let mut feature = Feature::new(2., 100., 1.);
    feature.convex_hulls.push(
        ConvexHull2D::from_points(&[Point2D::new(0., 100.), Point2D::new(1., 101.)]).unwrap(),
    );
    let mut features = FeatureMap::from_features(vec![Feature::new(2., 100., 1.), feature]);
    let before = features.clone();
    assert!(
        transformer
            .transform_feature_map(&mut features, &transform)
            .is_err()
    );
    assert_eq!(features, before);
    let mut consensus =
        ConsensusMap::from_features(vec![ConsensusFeature::from(BaseFeature::new(2., 100., 1.))]);
    consensus
        .unassigned_peptide_identifications
        .push(peptide(Some(0.)));
    let before = consensus.clone();
    assert!(
        transformer
            .transform_consensus_map(&mut consensus, &transform)
            .is_err()
    );
    assert_eq!(consensus, before);
    let mut ids = vec![peptide(Some(2.)), peptide(Some(0.))];
    let before = ids.clone();
    assert!(
        transformer
            .transform_peptide_identifications(&mut ids, &transform)
            .is_err()
    );
    assert_eq!(ids, before);
}

#[test]
fn checked_limits_fail_without_changes_and_existing_original_keys_are_retained() {
    let mut feature = Feature::new(2., 100., 1.);
    feature.subordinates.push(Feature::new(1., 100., 1.));
    let before = feature.clone();
    for transformer in [
        MapAlignmentTransformer {
            max_subordinate_depth: 0,
            ..Default::default()
        },
        MapAlignmentTransformer {
            max_records: 1,
            ..Default::default()
        },
        MapAlignmentTransformer {
            max_rt_values: 1,
            ..Default::default()
        },
        MapAlignmentTransformer {
            max_rt_values: 0,
            ..Default::default()
        },
    ] {
        assert!(
            transformer
                .transform_feature(&mut feature, &linear(2., 1.))
                .is_err()
        );
        assert_eq!(feature, before);
    }
    let mut ids = vec![peptide(None), peptide(Some(2.))];
    ids[1]
        .metadata
        .insert("original_RT".into(), "already stored".into());
    MapAlignmentTransformer {
        store_original_rt: true,
        ..Default::default()
    }
    .transform_peptide_identifications(&mut ids, &linear(2., 1.))
    .unwrap();
    assert_eq!(
        ids[1].metadata["original_RT"].as_str().unwrap(),
        "already stored"
    );
    assert!(ids[0].metadata.is_empty());
    let mut spectrum = MSSpectrum {
        rt: 1.,
        ..Default::default()
    };
    spectrum
        .float_data_arrays
        .push(DataArray::new("bad", vec![1.]));
    let mut experiment = MSExperiment {
        spectra: vec![spectrum],
        ..Default::default()
    };
    let before = experiment.clone();
    assert!(
        MapAlignmentTransformer::default()
            .transform_experiment(&mut experiment, &linear(2., 1.))
            .is_err()
    );
    assert_eq!(experiment, before);
}
