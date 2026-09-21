// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Native functions for primitive and string list operations, dates and
//! `D`-dimensional positions and ranges.

/// Native ListUtils operations and conversions.
pub mod list;
/// String list prefix/suffix searches and C-locale case conversion.
pub mod string_list;

/// Naive source DateTime state, parsing and Gregorian arithmetic.
pub mod datetime;
pub use datetime::DateTime;

/// Owned CV mapping records and ordered container operations.
pub mod cv_mapping;
pub use cv_mapping::{
    CVMappingRule, CVMappingTerm, CVMappings, CVReference, CombinationsLogic, RequirementLevel,
};

/// `D`-dimensional coordinates (source `DATASTRUCTURES/DPosition.h`).
pub mod dposition;
pub use dposition::{DPosition, DPosition1, DPosition2};

/// Closed `D`-dimensional intervals (source `DATASTRUCTURES/DIntervalBase.h`).
pub mod dinterval;
pub use dinterval::{DIntervalBase, DIntervalBase1, DIntervalBase2};

/// Half-open `D`-dimensional ranges (source `DATASTRUCTURES/DRange.h`).
pub mod drange;
pub use drange::{DRange, DRange1, DRange2, DRangeIntersection};

/// Tool descriptions for the TOPP tool registry (source
/// `DATASTRUCTURES/ToolDescription.h`).
pub mod tool_description;
pub use tool_description::{
    FileMapping, MappingParam, ToolDescription, ToolDescriptionInternal, ToolExternalDetails,
};

/// Tool metadata for the tool-description writers (source
/// `DATASTRUCTURES/ToolInfo.h`).
pub mod tool_info;
pub use tool_info::ToolInfo;
