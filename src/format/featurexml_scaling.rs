// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Size-derived resource ceilings for the featureXML reader and writer.
//!
//! Source `FeatureXMLHandler` is a Xerces SAX handler with no resource
//! ceilings: it keeps pushing features into the map until the document ends.
//! The port bounds every cumulative quantity instead, but a *fixed* ceiling
//! refuses real files once they are large enough. Before this module the
//! reader's effective input ceiling was 12,500,000 bytes — `max_xml_bytes`
//! (64 MiB) reduced by `max_payload_bytes / 8` (32 MiB) and by
//! `max_work / 4` (12.5 MB) — so a 59.6 MiB benchmark featureXML failed with
//! `identification XML byte limit exceeded` and a 2.06 GiB one could not be
//! opened at all, while the C++ `FileInfo` summarised both.
//!
//! An [`Allowance`] instead grows with the input the reader has actually
//! consumed, so its ceiling is `floor + units * (consumed / per_bytes)`.
//! Work and storage stay linear in the input, a document cannot amplify a few
//! bytes into unbounded work, and a document of any realistic size fits. The
//! writer uses the same shape against the size of the map it is asked to write.
//! See `docs/FEATUREXML_SCALE_SUPPORT.md` for the measured ratios the defaults
//! are derived from.
//!
//! This mirrors `src/format/mzml_scaling.rs`, the established pattern, and
//! repeats its [`Allowance`] rather than importing it: that module is compiled
//! only with the `mzml` feature, which `featurexml` does not imply.

/// A cumulative allowance that grows with the size of the work item.
///
/// The ceiling after `consumed` units of input is
/// `floor + units * (consumed / per_bytes)`, saturating at `usize::MAX`.
/// Growth is credited from input that has already been measured, never in
/// advance, so a small document is held to roughly its floor however large the
/// quantities it declares. `per_bytes` lets a rate below one unit per byte be
/// expressed exactly: an XML element needs a start tag, so element ceilings
/// grow once per group of bytes rather than once per byte.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Allowance {
    /// Units available before the first input byte.
    pub floor: usize,
    /// Units credited for each consumed group of `per_bytes` input units.
    pub units: usize,
    /// The consumed-input group that credits `units`. Zero is read as one.
    pub per_bytes: usize,
}
impl Allowance {
    /// An allowance of `floor` units plus `per_byte` units per consumed unit.
    pub const fn new(floor: usize, per_byte: usize) -> Self {
        Self {
            floor,
            units: per_byte,
            per_bytes: 1,
        }
    }
    /// An allowance of `floor` units plus one unit per consumed `per_bytes`
    /// input units. `per_bytes` of zero is read as one.
    pub const fn every(floor: usize, per_bytes: usize) -> Self {
        Self {
            floor,
            units: 1,
            per_bytes: if per_bytes == 0 { 1 } else { per_bytes },
        }
    }
    /// A fixed allowance that does not grow with the input.
    pub const fn fixed(floor: usize) -> Self {
        Self {
            floor,
            units: 0,
            per_bytes: 1,
        }
    }
    /// The ceiling after `consumed` input units, saturating at `usize::MAX`.
    #[must_use]
    pub fn after(self, consumed: usize) -> usize {
        let groups = consumed / self.per_bytes.max(1);
        self.floor.saturating_add(self.units.saturating_mul(groups))
    }
    /// The ceiling after `consumed` input units, never above `cap`.
    #[must_use]
    pub fn capped(self, consumed: usize, cap: usize) -> usize {
        self.after(consumed).min(cap)
    }
}

/// Size-derived allowances for every cumulative quantity the featureXML reader
/// charges, applied against the decoded size of the document in bytes.
///
/// Each allowance applies together with the matching absolute ceiling in
/// [`super::Limits`]: a charge fails when it exceeds either one. The absolute
/// ceilings default to "unbounded", so the defaults here decide; a caller that
/// sets an explicit absolute ceiling keeps it exactly.
///
/// The floors are the fixed ceilings this reader used before, so a small
/// document is bounded as tightly as it was. The per-byte growth is at least
/// seven times the largest ratio measured on the benchmark featureXML inputs
/// (a 2.06 GiB `MassTraceExtractor` map with 826,019 features, a 59.6 MiB
/// `FeatureFinderCentroided` map with 42,789 features and a 188 KiB sanity
/// map), whose charges are recorded in `docs/FEATUREXML_SCALE_SUPPORT.md`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InputScaling {
    /// XML elements in the document. The densest benchmark input opens one
    /// element per 42.5 bytes; default 1,000,000 plus one per 8 bytes, 5.3
    /// times that density. The floor is the former fixed ceiling.
    pub records: Allowance,
    /// Entries in one metadata or evidence list. Lists are a property of a
    /// single value rather than of the document, so this grows slowly:
    /// default 1,000,000, the former fixed ceiling, plus one per 1,024 bytes.
    pub list_items: Allowance,
    /// Cumulative parser and conversion work. Measured at most 8.08 units per
    /// decoded byte; default 50,000,000, the former fixed ceiling, plus 64 per
    /// byte, 7.9 times that.
    pub work: Allowance,
    /// Cumulative payload estimate: the conservative storage charge the reader
    /// makes before every copy, which is a multiple of the storage actually
    /// held and not a byte count. Streaming a featureXML charges at most 8.2
    /// per decoded byte, 16.9 on a file too small for the ratio to settle;
    /// default 256 MiB, the former fixed ceiling, plus 32 per byte. That is
    /// four times the streaming charge and below the roughly 34 per byte a
    /// whole retained tree would cost, so a document whose bulk cannot be
    /// streamed — one feature that is the entire file, say — is refused rather
    /// than held.
    pub payload_bytes: Allowance,
}
impl Default for InputScaling {
    fn default() -> Self {
        Self {
            records: Allowance::every(1_000_000, 8),
            list_items: Allowance::every(1_000_000, 1_024),
            work: Allowance::new(50_000_000, 64),
            payload_bytes: Allowance::new(256 * 1024 * 1024, 32),
        }
    }
}
impl InputScaling {
    /// The same floors with no growth, which reproduces the reader's former
    /// fixed ceilings.
    #[must_use]
    pub fn fixed(self) -> Self {
        Self {
            records: Allowance::fixed(self.records.floor),
            list_items: Allowance::fixed(self.list_items.floor),
            work: Allowance::fixed(self.work.floor),
            payload_bytes: Allowance::fixed(self.payload_bytes.floor),
        }
    }
}

/// Size-derived allowances for the featureXML writer, applied against the
/// counted size of the map being written.
///
/// The counted size is one unit per feature, one per convex-hull point, one
/// per identification and one per metadata entry, summed over subordinates —
/// the quantities the writer charges for. As for the reader, each allowance
/// applies together with the matching absolute ceiling in [`super::Limits`],
/// whose defaults are unbounded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OutputScaling {
    /// XML elements rendered. A hull point and a feature are one element each;
    /// default 1,000,000 plus 4 per counted unit.
    pub records: Allowance,
    /// Entries in one metadata list; default 1,000,000 plus one per 64 units.
    pub list_items: Allowance,
    /// Cumulative rendering work; default 50,000,000 plus 1,024 per unit.
    pub work: Allowance,
    /// Cumulative payload estimate. The writer charges up to 2,048 bytes per
    /// hull point and about 5,000 per feature before rendering; default
    /// 256 MiB plus 16,384 per unit.
    pub payload_bytes: Allowance,
}
impl Default for OutputScaling {
    fn default() -> Self {
        Self {
            records: Allowance::new(1_000_000, 4),
            list_items: Allowance::every(1_000_000, 64),
            work: Allowance::new(50_000_000, 1_024),
            payload_bytes: Allowance::new(256 * 1024 * 1024, 16_384),
        }
    }
}
impl OutputScaling {
    /// The same floors with no growth, which reproduces the writer's former
    /// fixed ceilings.
    #[must_use]
    pub fn fixed(self) -> Self {
        Self {
            records: Allowance::fixed(self.records.floor),
            list_items: Allowance::fixed(self.list_items.floor),
            work: Allowance::fixed(self.work.floor),
            payload_bytes: Allowance::fixed(self.payload_bytes.floor),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_ceiling_grows_with_the_consumed_input_and_saturates() {
        let allowance = Allowance::new(10, 2);
        assert_eq!(allowance.after(0), 10);
        assert_eq!(allowance.after(3), 16);
        assert_eq!(Allowance::new(1, usize::MAX).after(usize::MAX), usize::MAX);
        assert_eq!(Allowance::new(usize::MAX, 1).after(1), usize::MAX);
    }

    #[test]
    fn a_rate_below_one_unit_per_byte_is_exact() {
        let allowance = Allowance::every(4, 512);
        assert_eq!(allowance.after(0), 4);
        assert_eq!(allowance.after(511), 4);
        assert_eq!(allowance.after(512), 5);
        assert_eq!(allowance.after(4_096), 12);
        assert_eq!(Allowance::every(0, 0).after(7), 7); // zero reads as one
    }

    #[test]
    fn a_fixed_allowance_and_a_cap_hold_the_floor() {
        assert_eq!(Allowance::fixed(7).after(u32::MAX as usize), 7);
        assert_eq!(Allowance::new(10, 2).capped(1_000, 42), 42);
        assert_eq!(Allowance::new(10, 2).capped(1, 42), 12);
        let fixed = InputScaling::default().fixed();
        assert_eq!(
            fixed.work.after(1 << 40),
            InputScaling::default().work.floor
        );
        assert_eq!(fixed.records.after(1 << 40), 1_000_000);
        let fixed = OutputScaling::default().fixed();
        assert_eq!(fixed.payload_bytes.after(1 << 40), 256 * 1024 * 1024);
        assert_eq!(fixed.list_items.after(1 << 40), 1_000_000);
    }
}
