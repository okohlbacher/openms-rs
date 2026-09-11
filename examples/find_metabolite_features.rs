// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Find metabolite features in centroided mzML using the source default chain.
//! This SDK example is not a complete TOPP command-line replacement.

#[cfg(all(feature = "mzml", feature = "featurexml"))]
fn main() -> openms::Result<()> {
    use openms::{
        analysis::{
            elution_peak_detection::ElutionPeakDetection,
            feature_finding_metabo::FeatureFindingMetabo, mass_trace_detection::MassTraceDetection,
        },
        concept::unique_id::UniqueIdGenerator,
        format::{featurexml, mzml},
    };
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 2 {
        return Err(openms::Error::InvalidValue(
            "usage: find_metabolite_features input.mzML output.featureXML".into(),
        ));
    }
    let experiment = mzml::load(&args[0])?;
    let mut traces = MassTraceDetection::new().run(&experiment, 0)?;
    let mut split = ElutionPeakDetection::new().detect_peaks_many(&mut traces)?;
    let output = FeatureFindingMetabo::new()?.run(&mut split, &mut UniqueIdGenerator::new())?;
    featurexml::store(&args[1], &output.features)?;
    println!(
        "Wrote {} metabolite features",
        output.features.features.len()
    );
    Ok(())
}

#[cfg(not(all(feature = "mzml", feature = "featurexml")))]
fn main() -> openms::Result<()> {
    Err(openms::Error::Unsupported(
        "find_metabolite_features requires the mzml and featurexml features".into(),
    ))
}
