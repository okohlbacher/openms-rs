// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Size-derived resource ceilings for the mzML reader.
//!
//! Source `MzMLHandler` has no resource ceilings: it reserves whatever a
//! document declares and keeps reading. The port bounds every cumulative
//! allowance, but a fixed ceiling rejects real files once they are large
//! enough; every benchmark input from 0.5 to 2.3 GB failed that way. An
//! [`Allowance`] instead grows with the XML bytes the reader has actually
//! consumed, so its ceiling is `floor + rate * consumed`. Work and storage
//! then stay linear in the input, a document cannot amplify a few bytes into
//! unbounded work, and a document of any realistic size fits. See
//! `docs/MZML_READER_SCALE_SUPPORT.md` for the measured ratios the defaults are
//! derived from.

/// A cumulative reader allowance that grows with the consumed input.
///
/// The ceiling after `consumed` XML bytes is
/// `floor + units * (consumed / per_bytes)`, saturating at `usize::MAX`.
/// Growth is credited as bytes are consumed, never in advance, so a small
/// document is held to roughly its floor however large the quantities it
/// declares. `per_bytes` lets a rate below one unit per byte be expressed
/// exactly: records, binary arrays and parameter groups each need their own
/// start tag, so their ceilings grow once per group of bytes rather than once
/// per byte.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Allowance {
    /// Units available before the first input byte.
    pub floor: usize,
    /// Units credited for each consumed group of `per_bytes` input bytes.
    pub units: usize,
    /// The consumed-byte group that credits `units`. Zero is read as one.
    pub per_bytes: usize,
}
impl Allowance {
    /// An allowance of `floor` units plus `per_byte` units per consumed byte.
    pub const fn new(floor: usize, per_byte: usize) -> Self {
        Self {
            floor,
            units: per_byte,
            per_bytes: 1,
        }
    }
    /// An allowance of `floor` units plus one unit per consumed `per_bytes`
    /// input bytes. `per_bytes` of zero is read as one.
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
    /// The ceiling after `consumed` XML bytes, saturating at `usize::MAX`.
    pub fn after(self, consumed: u64) -> usize {
        let consumed = usize::try_from(consumed).unwrap_or(usize::MAX);
        let groups = consumed / self.per_bytes.max(1);
        self.floor.saturating_add(self.units.saturating_mul(groups))
    }
}

/// Size-derived allowances for every cumulative quantity the mzML reader and
/// counter charge.
///
/// Each allowance applies together with the matching absolute ceiling in
/// [`super::ReadOptions`] or [`super::LoadOptions`]: a charge fails when it
/// exceeds either one. The absolute ceilings default to "unbounded", so the
/// defaults here decide; a caller that sets an explicit absolute ceiling keeps
/// it exactly.
///
/// The floors are the fixed ceilings this reader used before, so a small
/// document is bounded as tightly as it was. The per-byte growth is at least
/// eight times the largest ratio measured on the benchmark inputs (a 2.3 GB
/// uncompressed Q Exactive run, a 1.2 GB zlib-compressed LTQ Orbitrap Velos
/// run, a 1.5 GB LTQ Orbitrap XL run, a 0.5 GB centroided run and two
/// sanity files), and also covers encodings those inputs do not use, such as
/// Numpress and highly compressible zlib arrays.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InputScaling {
    /// Declared spectrum and chromatogram points, charged when a record opens.
    /// Measured at most 0.087 per byte; default 10,000,000 plus 8 per byte.
    pub peaks: Allowance,
    /// Decoded bytes across all binary arrays. Measured at most 1.19 per byte
    /// (zlib); default 512 MiB plus 64 per byte.
    pub array_bytes: Allowance,
    /// Elements across all binary arrays. Measured at most 0.17 per byte;
    /// default 20,000,000 plus 16 per byte.
    pub array_elements: Allowance,
    /// Spectra plus chromatograms, charged when a record opens. Each one needs
    /// its own start tag, so the rate is one record per consumed byte group.
    /// The densest benchmark input writes one record per 3,951 bytes; default
    /// 1,000,000 plus one per 512 bytes, 7.7 times that density. The floor is
    /// the former fixed ceiling and already exceeds the largest benchmark input
    /// (40,857 records) 24-fold.
    pub records: Allowance,
    /// Binary arrays, including empty placeholders. The densest benchmark input
    /// writes one per 3,559 bytes; default 1,000,000 plus one per 256 bytes,
    /// 13.9 times that density and twice the record rate, since a record
    /// carries at least two arrays.
    pub arrays: Allowance,
    /// Referenceable parameter group definitions. Benchmark inputs define at
    /// most two in a gigabyte; default 100,000, the former fixed ceiling, plus
    /// one per 4,096 bytes.
    pub param_groups: Allowance,
    /// Parameters, parameter-group expansions and acquisition descriptors.
    /// Measured at most 0.007 per byte; default 10,000,000 plus 4 per byte.
    pub params: Allowance,
    /// The parameter storage estimate, which charges transient attribute maps
    /// as well as retained values. Measured at most 19.4 per byte, on a small
    /// file of sparse MS2 spectra; default 512 MiB plus 256 per byte.
    pub param_bytes: Allowance,
    /// Header, reference-resolution and record-metadata work. Measured at
    /// most 0.11 per byte; default 50,000,000 plus 64 per byte.
    pub metadata_work: Allowance,
    /// Header, reference-resolution and record-metadata storage. Measured at
    /// most 0.11 per byte; default 256 MiB plus 64 per byte.
    pub metadata_bytes: Allowance,
    /// Scientific selection work: range and MS-level checks, sorting and
    /// aligned permutation. Default 500,000,000 plus 512 per byte.
    pub selection_work: Allowance,
    /// Scientific selection index storage, held per record. Default 256 MiB
    /// plus 64 per byte.
    pub selection_bytes: Allowance,
    /// Counting work of `read_size`: one unit per XML event plus reference
    /// expansion. Default 50,000,000 plus 16 per byte.
    pub count_work: Allowance,
}
impl Default for InputScaling {
    fn default() -> Self {
        Self {
            peaks: Allowance::new(10_000_000, 8),
            array_bytes: Allowance::new(512 * 1024 * 1024, 64),
            array_elements: Allowance::new(20_000_000, 16),
            records: Allowance::every(1_000_000, 512),
            arrays: Allowance::every(1_000_000, 256),
            param_groups: Allowance::every(100_000, 4_096),
            params: Allowance::new(10_000_000, 4),
            param_bytes: Allowance::new(512 * 1024 * 1024, 256),
            metadata_work: Allowance::new(super::header::MAX_WORK, 64),
            metadata_bytes: Allowance::new(super::header::MAX_BYTES, 64),
            selection_work: Allowance::new(500_000_000, 512),
            selection_bytes: Allowance::new(256 * 1024 * 1024, 64),
            count_work: Allowance::new(super::counts::MAX_COUNT_WORK, 16),
        }
    }
}
impl InputScaling {
    /// The same floors with no growth: every allowance is fixed at its floor,
    /// which reproduces the reader's former fixed ceilings.
    pub fn fixed(self) -> Self {
        let fixed = |a: Allowance| Allowance::fixed(a.floor);
        Self {
            peaks: fixed(self.peaks),
            array_bytes: fixed(self.array_bytes),
            array_elements: fixed(self.array_elements),
            records: fixed(self.records),
            arrays: fixed(self.arrays),
            param_groups: fixed(self.param_groups),
            params: fixed(self.params),
            param_bytes: fixed(self.param_bytes),
            metadata_work: fixed(self.metadata_work),
            metadata_bytes: fixed(self.metadata_bytes),
            selection_work: fixed(self.selection_work),
            selection_bytes: fixed(self.selection_bytes),
            count_work: fixed(self.count_work),
        }
    }
}

/// One allowance reconciled against a plain `usize` counter that the reader
/// spends from directly.
///
/// The counter always holds the smaller of the room left under the absolute
/// ceiling and the room left under the size-derived one. Spending and refunds
/// between two reconciliations are read back from the counter.
struct Track {
    /// The absolute ceiling the counter held when it was attached.
    cap: usize,
    /// The size-derived allowance.
    allowance: Allowance,
    /// Everything charged so far, less everything refunded.
    spent: usize,
    /// The counter value written at the previous reconciliation.
    last: usize,
}
impl Track {
    /// Attach to `counter`, which holds the room left under its absolute
    /// ceiling, and clamp it to the floor.
    fn attach(counter: &mut usize, allowance: Allowance) -> Self {
        let cap = *counter;
        *counter = cap.min(allowance.floor);
        Self {
            cap,
            allowance,
            spent: 0,
            last: *counter,
        }
    }
    /// Read back what was spent or refunded since the previous call and write
    /// the room left after `consumed` input bytes to `counter`.
    ///
    /// The room is recomputed from the totals rather than credited in steps, so
    /// a rate below one unit per byte is exact however the input is chunked.
    fn sync(&mut self, counter: &mut usize, consumed: u64) {
        if *counter <= self.last {
            self.spent = self.spent.saturating_add(self.last - *counter);
        } else {
            self.spent = self.spent.saturating_sub(*counter - self.last);
        }
        *counter = self
            .cap
            .saturating_sub(self.spent)
            .min(self.allowance.after(consumed).saturating_sub(self.spent));
        self.last = *counter;
    }
}

/// Size-derived reconciliation for a fixed set of counters.
///
/// Attach once with each counter holding the room left under its absolute
/// ceiling, then call [`Ledger::sync`] with the input position after every XML
/// event. A counter that the reader decrements with `checked_sub` then fails
/// exactly when the charge exceeds either ceiling.
pub(super) struct Ledger<const N: usize> {
    /// One track per counter, in attachment order.
    tracks: [Track; N],
    /// Input bytes already credited.
    consumed: u64,
}
impl<const N: usize> Ledger<N> {
    /// Attach every `(counter, allowance)` pair, clamping each counter to its
    /// allowance's floor.
    pub(super) fn attach(counters: [(&mut usize, Allowance); N]) -> Self {
        Self {
            tracks: counters.map(|(counter, allowance)| Track::attach(counter, allowance)),
            consumed: 0,
        }
    }
    /// Reconcile the counters, in attachment order, with the input position
    /// `position` (total XML bytes consumed so far).
    pub(super) fn sync(&mut self, position: u64, counters: [&mut usize; N]) {
        self.consumed = self.consumed.max(position);
        for (track, counter) in self.tracks.iter_mut().zip(counters) {
            track.sync(counter, self.consumed);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_counter_holds_the_smaller_of_both_ceilings_and_grows_with_input() {
        let mut counter = 100;
        let mut ledger = Ledger::attach([(&mut counter, Allowance::new(10, 2))]);
        assert_eq!(counter, 10);
        counter -= 4;
        ledger.sync(3, [&mut counter]);
        assert_eq!(counter, 12); // 10 + 2 * 3 - 4
        counter -= 12;
        ledger.sync(3, [&mut counter]);
        assert_eq!(counter, 0);
        ledger.sync(1_000, [&mut counter]);
        // The absolute ceiling 100 minus the 16 spent wins over the scaled room.
        assert_eq!(counter, 84);
    }

    #[test]
    fn refunds_are_returned_to_both_ceilings() {
        let mut counter = 50;
        let mut ledger = Ledger::attach([(&mut counter, Allowance::fixed(20))]);
        counter -= 15;
        ledger.sync(0, [&mut counter]);
        assert_eq!(counter, 5);
        counter += 10;
        ledger.sync(0, [&mut counter]);
        assert_eq!(counter, 15);
    }

    #[test]
    fn a_rate_below_one_unit_per_byte_is_exact_however_the_input_is_chunked() {
        // One unit per 512 bytes, floor 4: 4 units at 0 bytes, 5 at 512.
        let allowance = Allowance::every(4, 512);
        assert_eq!(allowance.after(0), 4);
        assert_eq!(allowance.after(511), 4);
        assert_eq!(allowance.after(512), 5);
        assert_eq!(allowance.after(4_096), 12);
        assert_eq!(Allowance::every(0, 0).after(7), 7); // zero reads as one
        // Reaching 4,096 bytes in 8 steps or in one gives the same room, where
        // per-step crediting would have truncated every step to zero.
        for step in [1, 7, 512, 4_096] {
            let mut counter = usize::MAX;
            let mut ledger = Ledger::attach([(&mut counter, allowance)]);
            let mut position = 0;
            while position < 4_096 {
                position = (position + step).min(4_096);
                ledger.sync(position, [&mut counter]);
            }
            assert_eq!(counter, 12, "step {step}");
        }
    }

    #[test]
    fn saturating_arithmetic_never_panics() {
        let mut counter = usize::MAX;
        let mut ledger = Ledger::attach([(&mut counter, Allowance::new(usize::MAX, usize::MAX))]);
        ledger.sync(u64::MAX, [&mut counter]);
        assert_eq!(counter, usize::MAX);
        assert_eq!(Allowance::new(1, usize::MAX).after(u64::MAX), usize::MAX);
        assert_eq!(Allowance::new(7, 3).after(5), 22);
    }
}
