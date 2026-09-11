// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::{
    MSChromatogram, MSSpectrum, Result, interfaces::MSDataConsumer, metadata::ExperimentalSettings,
};
use std::ops::ControlFlow;
struct Consumer(usize);
impl MSDataConsumer for Consumer {
    fn set_expected_size(&mut self, s: usize, c: usize) -> Result<()> {
        self.0 += s + c;
        Ok(())
    }
    fn set_experimental_settings(&mut self, _: &ExperimentalSettings) -> Result<()> {
        self.0 += 1;
        Ok(())
    }
    fn consume_spectrum(&mut self, s: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        s.name = "s".into();
        Ok(ControlFlow::Continue(()))
    }
    fn consume_chromatogram(&mut self, c: &mut MSChromatogram) -> Result<ControlFlow<()>> {
        c.name = "c".into();
        Ok(ControlFlow::Break(()))
    }
}
#[test]
fn four_required_callbacks_are_object_safe_without_any_format_feature() {
    let mut consumer = Consumer(0);
    let object: &mut dyn MSDataConsumer = &mut consumer;
    object.set_expected_size(2, 3).unwrap();
    object
        .set_experimental_settings(&ExperimentalSettings::default())
        .unwrap();
    let mut s = MSSpectrum::default();
    let mut c = MSChromatogram::default();
    assert_eq!(
        object.consume_spectrum(&mut s).unwrap(),
        ControlFlow::Continue(())
    );
    assert_eq!(
        object.consume_chromatogram(&mut c).unwrap(),
        ControlFlow::Break(())
    );
    assert_eq!((s.name.as_str(), c.name.as_str()), ("s", "c"));
    assert_eq!(consumer.0, 6);
}
