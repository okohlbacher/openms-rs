// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Processing annotations for TOPP tool outputs, the native form of the source
//! `TOPPBase::getProcessingInfo_` and `TOPPBase::addDataProcessing_`
//! (`TOPPBase.cpp:527-605`).
//!
//! A tool builds one [`DataProcessing`] record with
//! [`ToolContext::processing_info`](super::ToolContext::processing_info) and
//! attaches it to its output with
//! [`ToolContext::add_data_processing`](super::ToolContext::add_data_processing).
//! See `docs/TOPP_CLI_SUPPORT.md` for the supported source subset.

use crate::Result;
use crate::data_structures::DateTime;
use crate::kernel::{ConsensusMap, FeatureMap, MSExperiment};
use crate::metadata::{DataProcessing, MetaValue, ProcessingAction};
use crate::param::{Param, ParamValue};
use crate::system::file;
use std::sync::Arc;

/// Software version recorded under `-test`, as `getProcessingInfo_` sets it.
pub const TEST_MODE_VERSION: &str = "version_string";
/// Completion time recorded under `-test`, as `getProcessingInfo_` sets it.
pub const TEST_MODE_COMPLETION_TIME: &str = "1999-12-31 23:59:59";
/// The only metadata key recorded under `-test`.
pub const TEST_MODE_PARAMETER_KEY: &str = "parameter: mode";
/// The value of [`TEST_MODE_PARAMETER_KEY`].
pub const TEST_MODE_PARAMETER_VALUE: &str = "test_mode";

/// Build the processing record of one tool run.
///
/// Source `getProcessingInfo_`. The software is named after the tool. Under
/// `-test` the version is [`TEST_MODE_VERSION`], the completion time
/// [`TEST_MODE_COMPLETION_TIME`] and the only metadata entry
/// `parameter: mode` = `test_mode`, so test outputs do not depend on the clock
/// or on the machine's parameter values. Otherwise the version is the tool's
/// product version, the completion time is now in local time, and every
/// resolved parameter is recorded as `parameter: <name>` with its typed value.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`](crate::Error::InvalidValue) when a
/// floating-point parameter cannot be represented as metadata, and propagates
/// parameter-tree traversal failures.
pub(crate) fn processing_info(
    tool_name: &str,
    version: &str,
    test_mode: bool,
    param: &Param,
    actions: &[ProcessingAction],
) -> Result<DataProcessing> {
    let mut processing = DataProcessing {
        actions: actions.iter().cloned().collect(),
        ..DataProcessing::default()
    };
    processing.software.name = tool_name.to_owned();
    if test_mode {
        processing.software.version = TEST_MODE_VERSION.to_owned();
        processing.completion_time = Some(DateTime::parse(TEST_MODE_COMPLETION_TIME)?);
        processing.metadata.insert(
            TEST_MODE_PARAMETER_KEY.to_owned(),
            MetaValue::from(TEST_MODE_PARAMETER_VALUE),
        );
    } else {
        processing.software.version = version.to_owned();
        processing.completion_time = Some(DateTime::now());
        for item in param.iter()? {
            processing.metadata.insert(
                format!("parameter: {}", item.key),
                meta_value(&item.entry.value)?,
            );
        }
    }
    Ok(processing)
}

/// The typed metadata form of a parameter value; source `ParamValue` converts
/// to `DataValue` implicitly.
fn meta_value(value: &ParamValue) -> Result<MetaValue> {
    Ok(match value {
        ParamValue::Empty => MetaValue::default(),
        ParamValue::String(text) => MetaValue::from(text.clone()),
        ParamValue::Integer(number) => MetaValue::from(*number),
        ParamValue::Float(number) => MetaValue::try_from(*number)?,
        ParamValue::StringList(values) => MetaValue::from(values.clone()),
        ParamValue::IntegerList(values) => {
            MetaValue::from(values.iter().map(|v| i64::from(*v)).collect::<Vec<i64>>())
        }
        ParamValue::FloatList(values) => MetaValue::try_from(values.clone())?,
    })
}

/// An output a TOPP tool attaches its processing record to.
///
/// Source `addDataProcessing_`, which is overloaded for the three map types.
/// The work is one clone or one shared handle per record already held in
/// memory, so it needs no separate bound.
pub trait AddDataProcessing {
    /// Attach `processing`; `test_mode` is the tool's `-test` flag.
    fn add_data_processing(&mut self, processing: &DataProcessing, test_mode: bool);
}

impl AddDataProcessing for FeatureMap {
    /// Appends the record to the map's own processing list.
    fn add_data_processing(&mut self, processing: &DataProcessing, _test_mode: bool) {
        self.data_processing.push(processing.clone());
    }
}

impl AddDataProcessing for ConsensusMap {
    /// Appends the record to the map's processing list. Under `-test` every
    /// column header's file name is reduced to its base name, as the source
    /// removes absolute map paths from test outputs.
    fn add_data_processing(&mut self, processing: &DataProcessing, test_mode: bool) {
        self.data_processing.push(processing.clone());
        if test_mode {
            for header in self.column_headers.values_mut() {
                header.filename = file::basename(&header.filename).to_owned();
            }
        }
    }
}

impl AddDataProcessing for MSExperiment {
    /// Appends one shared handle to every spectrum and every chromatogram, as
    /// the source pushes one `shared_ptr` everywhere. Experiment-level settings
    /// are left alone.
    fn add_data_processing(&mut self, processing: &DataProcessing, _test_mode: bool) {
        let shared = Arc::new(processing.clone());
        for spectrum in &mut self.spectra {
            spectrum.data_processing.push(Arc::clone(&shared));
        }
        for chromatogram in &mut self.chromatograms {
            chromatogram.data_processing.push(Arc::clone(&shared));
        }
    }
}
