// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The last two members of `KERNEL/ConsensusFeature.h`: the stream operator and
//! the ratio description.
//!
//! Closes the two gaps the feature-identification package left open, because
//! each needed something it did not own — a
//! [`Display`](std::fmt::Display) implementation for
//! [`ConsensusFeature`](crate::kernel::features::ConsensusFeature) reproducing
//! `operator<<(std::ostream&, const ConsensusFeature&)`
//! (`ConsensusFeature.cpp:391-418`), and checked operations over
//! [`Ratio::description`](crate::kernel::features::Ratio::description), the
//! source `Ratio::description_`. See `docs/CONSENSUS_DISPLAY_SUPPORT.md`.

use super::features::{ConsensusFeature, Ratio};
use crate::{Error, Result};
use std::fmt;

/// Source `operator<<(std::ostream&, const ConsensusFeature&)`: a banner, the
/// position, intensity and quality, one indented block per grouped handle in
/// handle order, the meta values in key order, then a closing banner.
///
/// The layout is reproduced exactly, including the two trailing spaces the
/// source leaves behind: `"Grouped features: "` and `"Meta information: "` end
/// in a space before the newline because the source streams `"Grouped features:
/// \n"` with the space inside the literal, and the closing banner line likewise
/// ends `"----------------- "`. The `"Intensity "` and `"Quality "` labels carry
/// no colon while every other label does; that too is the source.
///
/// Three differences from the C++ stream, all observable and none of them a
/// choice about content:
///
/// * **Number formatting.** The source wraps intensity, quality and every handle
///   coordinate in `precisionWrapper`, which sets the stream precision to
///   `writtenDigits<T>` — 6 significant digits for the `f32` intensity and
///   quality, 15 for the `f64` coordinates — while the position goes through
///   `DPosition`'s own operator, also `precisionWrapper`. This writes Rust's
///   shortest round-trip form instead, as every other kernel `Display` in the
///   port does, so `1.5` prints as `1.5` and not `1.50000000000000`.
/// * **Meta value order.** The source `getKeys` yields names in
///   `MetaInfoRegistry` index order, i.e. the order in which the process first
///   registered each name. The port's `MetaInfo` is a `BTreeMap`, so names come
///   out sorted. The set of printed pairs is the same.
/// * **Ratios are not printed**, because the source does not print them either.
///
/// `ConsensusMap`'s own `Display` (`src/kernel/map_operations.rs`) writes the
/// same per-feature block through a private helper, because that module could
/// not call an implementation that did not exist yet. `tests/consensus_display.rs`
/// asserts the two agree character for character; collapsing the helper into
/// `{feature}` is left to whichever package next owns `map_operations.rs`.
impl fmt::Display for ConsensusFeature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "---------- CONSENSUS ELEMENT BEGIN -----------------")?;
        writeln!(f, "Position: {} {}", self.base.rt, self.base.mz)?;
        writeln!(f, "Intensity {}", self.base.intensity)?;
        writeln!(f, "Quality {}", self.base.quality)?;
        writeln!(f, "Grouped features: ")?;
        for handle in self.handles() {
            writeln!(f, " - Map index: {}", handle.map_index)?;
            writeln!(f, "   Feature id: {}", handle.unique_id)?;
            writeln!(f, "   RT: {}", handle.rt)?;
            writeln!(f, "   m/z: {}", handle.mz)?;
            writeln!(f, "   Intensity: {}", handle.intensity)?;
        }
        writeln!(f, "Meta information: ")?;
        for (key, value) in &self.base.metadata {
            writeln!(f, "   {key}: {value}")?;
        }
        writeln!(f, "---------- CONSENSUS ELEMENT END ----------------- ")
    }
}

/// Checked operations over the source `ConsensusFeature::Ratio::description_`.
///
/// The source member is a bare public `std::vector<std::string>` that nothing in
/// the SDK writes; it exists so that a quantification algorithm can record how a
/// ratio was formed. The Rust counterpart is the equally public
/// [`Ratio::description`](crate::kernel::features::Ratio::description) field —
/// direct assignment stays available and unchecked, exactly as in C++ — and
/// these operations are the checked path over it.
impl Ratio {
    /// Description lines one ratio may hold.
    ///
    /// Native ceiling; the source vector is unbounded. A description is
    /// free-text provenance for a single ratio, so a caller that reaches this
    /// has a bug rather than a large dataset.
    pub const MAX_DESCRIPTION_LINES: usize = 65_536;

    /// Total bytes of description one ratio may hold, counting the UTF-8 length
    /// of every line and nothing else.
    ///
    /// Native ceiling, for the reason given on
    /// [`Ratio::MAX_DESCRIPTION_LINES`].
    pub const MAX_DESCRIPTION_BYTES: usize = 1 << 20;

    /// The description lines, in insertion order.
    ///
    /// The source has no accessor at all: `description_` is read through the
    /// public member. This is the read half of that, and
    /// [`Ratio::description`](crate::kernel::features::Ratio::description)
    /// remains available for mutation without checks.
    pub fn description(&self) -> &[String] {
        &self.description
    }

    /// Append one description line.
    ///
    /// Duplicate lines are appended, as a `push_back` on the source vector
    /// would be; the lines are provenance, not a set. An empty line is a valid
    /// line, because the source places no constraint on the strings.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the line would take the ratio past
    /// [`Ratio::MAX_DESCRIPTION_LINES`] lines or
    /// [`Ratio::MAX_DESCRIPTION_BYTES`] bytes. Both ceilings are checked before
    /// anything is appended, so a rejected call leaves the description
    /// unchanged.
    ///
    /// The byte preflight sums the stored lines on every call, because a
    /// `Ratio` keeps no cached total; the sum is bounded by
    /// [`Ratio::MAX_DESCRIPTION_LINES`] additions and allocates nothing.
    pub fn add_description(&mut self, line: &str) -> Result<()> {
        if self.description.len() >= Self::MAX_DESCRIPTION_LINES {
            return Err(description_limit("lines", Self::MAX_DESCRIPTION_LINES));
        }
        let bytes = description_bytes(self.description.iter().map(String::as_str))?
            .checked_add(line.len())
            .ok_or_else(|| description_limit("bytes", Self::MAX_DESCRIPTION_BYTES))?;
        if bytes > Self::MAX_DESCRIPTION_BYTES {
            return Err(description_limit("bytes", Self::MAX_DESCRIPTION_BYTES));
        }
        self.description.push(line.to_string());
        Ok(())
    }

    /// Replace every description line.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `lines` exceeds
    /// [`Ratio::MAX_DESCRIPTION_LINES`] or
    /// [`Ratio::MAX_DESCRIPTION_BYTES`]. The whole input is measured before the
    /// stored description is touched, so a rejected call leaves the previous
    /// lines in place.
    pub fn set_description(&mut self, lines: Vec<String>) -> Result<()> {
        check_description(&lines)?;
        self.description = lines;
        Ok(())
    }

    /// Check that the stored description is within both ceilings.
    ///
    /// [`Ratio::validate`](crate::kernel::features::Ratio::validate) checks the
    /// ratio value only and does not call this, so a description assigned
    /// straight to the public field is not measured by
    /// [`ConsensusFeature::add_ratio`](crate::kernel::features::ConsensusFeature::add_ratio)
    /// or
    /// [`ConsensusFeature::set_ratios`](crate::kernel::features::ConsensusFeature::set_ratios)
    /// either. That mirrors the source, where nothing validates a `Ratio` at
    /// all; call this explicitly when the description came from untrusted input.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the stored description exceeds
    /// [`Ratio::MAX_DESCRIPTION_LINES`] or [`Ratio::MAX_DESCRIPTION_BYTES`].
    pub fn validate_description(&self) -> Result<()> {
        check_description(&self.description)
    }
}

/// The error both ceilings report.
fn description_limit(unit: &str, ceiling: usize) -> Error {
    Error::InvalidValue(format!(
        "consensus feature ratio description exceeds {ceiling} {unit}"
    ))
}

/// Total UTF-8 length of the given lines, refusing an overflowing sum.
fn description_bytes<'a>(lines: impl Iterator<Item = &'a str>) -> Result<usize> {
    let mut total = 0usize;
    for line in lines {
        total = total
            .checked_add(line.len())
            .ok_or_else(|| description_limit("bytes", Ratio::MAX_DESCRIPTION_BYTES))?;
    }
    Ok(total)
}

/// Preflight a candidate description against both ceilings.
fn check_description(lines: &[String]) -> Result<()> {
    if lines.len() > Ratio::MAX_DESCRIPTION_LINES {
        return Err(description_limit("lines", Ratio::MAX_DESCRIPTION_LINES));
    }
    if description_bytes(lines.iter().map(String::as_str))? > Ratio::MAX_DESCRIPTION_BYTES {
        return Err(description_limit("bytes", Ratio::MAX_DESCRIPTION_BYTES));
    }
    Ok(())
}
