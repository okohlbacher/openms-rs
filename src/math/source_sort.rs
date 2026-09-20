// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The orders in which the C++ Release build's `std::sort` and
//! `std::stable_sort` leave the elements they sort, for every sort the picked
//! feature finder (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`) reaches.
//!
//! `FeatureFinderAlgorithmPicked::run` and `run_` sort with `std::sort` six
//! times:
//!
//! - the spectra by retention time (`MSExperiment::sortSpectra`,
//!   `MSExperiment.cpp:793`) and the chromatograms by product m/z
//!   (`MSExperiment::sortChromatograms`, `:813`) of an unsorted input
//!   (`FeatureFinderAlgorithmPicked.cpp:1084-1085`);
//! - the user seeds by m/z (`FeatureMap::sortByMZ`, `:190`);
//! - the intensities of each step-1 cell (`std::sort` of a `std::vector<double>`,
//!   `:270`);
//! - the seeds by descending intensity (`:548`);
//! - the feature map by m/z (`FeatureMap::sortByMZ`, `:866`) and by descending
//!   intensity (`FeatureMap::sortByIntensity(true)`, `:991`).
//!
//! `std::sort` is not stable, and the standard leaves the order of equal
//! elements unspecified. That order is observable everywhere above: equal
//! retention times decide which scans are neighbours, equal product m/z the
//! chromatogram order, a `-0.0` next to a `+0.0` the sign of a stored quantile,
//! and equal features which one step 4 keeps (a reused instance or a caller's
//! non-empty map produces such ties routinely: running twice on one input
//! doubles every feature).
//!
//! The sorts of the spectra and chromatograms are followed by
//! `MSSpectrum::sortByPosition` and `MSChromatogram::sortByPosition`
//! (`MSSpectrum.cpp:444-460`, `MSChromatogram.cpp:165-183`), which sort each
//! unsorted spectrum's peaks by m/z and each chromatogram's peaks by retention
//! time with `std::stable_sort`: of the peaks when the spectrum has no data
//! array, and of an index vector (`MSSpectrum::sort`, `MSSpectrum.h:365-376`)
//! otherwise. Both compare the same keys in the same order, so they leave the
//! same permutation.
//!
//! The port follows the Linux x86_64 Release build, whose algorithms are those
//! of the conda-forge GCC 14.4.0 `libstdc++` headers it was compiled with
//! (`bits/stl_algo.h`, sha256 `0598c5b1...`; `bits/stl_heap.h`, `f18f83b2...`;
//! `bits/stl_tempbuf.h`, `f92e0ecf...`; `bits/stl_algobase.h`, `0ec2358c...`;
//! `bits/predefined_ops.h`, `25478342...`):
//!
//! - `std::sort` is an introsort: a quicksort that moves the median of the
//!   first, middle and last element to the front and partitions around it,
//!   recursing into the upper part, down to ranges of at most 16 elements or
//!   down to a recursion budget of `2 * floor(log2(n))`, where a heapsort takes
//!   over; and a final insertion sort over the whole range, whose first 16
//!   elements are inserted with a lower-bound check and the rest without one
//!   ([`source_sort_permutation`](crate::math::source_sort::source_sort_permutation));
//! - `std::stable_sort` asks `std::get_temporary_buffer` for `(n + 1) / 2`
//!   elements, which halves its request after every failed
//!   `operator new(nothrow)` until one succeeds or nothing is left. With the
//!   whole buffer it sorts both halves by an insertion sort over chunks of 7
//!   and buffered merges of doubling width, and merges the halves with
//!   `__merge_adaptive`; with a smaller buffer it recurses
//!   (`__stable_sort_adaptive_resize`, whose merges rotate through the buffer
//!   when a side does not fit); without one it merges in place
//!   (`__inplace_stable_sort`, insertion sort below 15 elements,
//!   `__merge_without_buffer` with `std::rotate`)
//!   ([`source_stable_sort_permutation`](crate::math::source_sort::source_stable_sort_permutation)).
//!
//! Both are reproduced comparison by comparison and move by move on a
//! permutation, so equal elements land where the executed C++ puts them. The
//! oracle drivers compare them with the library's own sorts:
//! `../oracle/ffap-instr-completion/drivers/ffap_instr_driver.cpp` (mode
//! `sort`, `FeatureMap::sortByMZ` and `sortByIntensity(true)` on 130 inputs)
//! and `../oracle/ffap-complete-fix1/drivers/sort_probe.cpp` (both
//! `sortByPosition`, with and without a data array, under a full, a partial and
//! no temporary buffer; `sortSpectra`, `sortChromatograms` and `sortByMZ`;
//! 2,272 inputs with ties, signed zeros, infinities and NaN keys of four bit
//! patterns). Where the port allocates,
//! [`TemporaryBuffer::Allocate`](crate::math::source_sort::TemporaryBuffer::Allocate)
//! halves its request on failure as `get_temporary_buffer` does; which
//! requests fail depends on the memory of the running process and is not
//! reproducible.
//!
//! # Undefined behaviour
//!
//! Both sorts require a strict weak ordering; `<` on floating-point keys that
//! include NaN is not one. The libstdc++ algorithms still run on such keys.
//! Every comparison and move of `std::stable_sort` stays inside the range
//! whatever the comparisons return (its merges check both ends, its binary
//! searches and rotations are bounded, and its insertion sort never passes the
//! first element, which the element it inserts was just found not to be less
//! than), so the port reproduces it on any keys. The introsort's partition
//! and final insertion loops have no bound checks: they rely on the ordering to
//! stop. The port follows them while they read inside the vector, and
//! returns [`Error::InvalidValue`] at exactly the step where the C++ would read
//! outside it, which is undefined behaviour with no reproducible result.
//!
//! That guard is never reached by a deterministic *asymmetric* comparison
//! (`less(a, b)` excludes `less(b, a)`), which every sort of the picked
//! feature finder uses: `<` on `f32` or `f64` keys is asymmetric with NaN keys
//! too. For such a comparison:
//!
//! - **Partition.** Whichever of the three elements `__move_median_to_first`
//!   picks as the pivot `p`, one of the other two is not less than `p` and
//!   stays in the range, so the upward scan stops inside it; every swap puts
//!   an element that is not less than `p` at the upper end, where it stops
//!   the later scans; and the downward scan stops at `p` itself, which is
//!   not greater than itself. Each partition therefore stays inside its range
//!   and its cut lies inside it.
//! - **Final insertion.** An element at position `i >= 16` is either in the
//!   upper part of one of the top-level partitions, whose elements are all not
//!   less than that partition's pivot `q`, while `q` stays below that part;
//!   or it was placed by the heapsort that takes over when the recursion budget
//!   of the top level runs out. That heapsort pops a root only while the heap
//!   still holds an element the root is not less than: the sibling the sift
//!   chose the root over, or the element the sift had moved to the top and the
//!   root was then pushed above, and neither leaves the heap before the root
//!   does. Either way an element the
//!   inserted one is not less than lies to its left, and the walk stops there.
//!
//! The executed and generated inputs of the oracles bear this out: none
//! reached the guard (`docs/FEATURE_FINDER_PICKED_SUPPORT.md`). A comparator
//! that is not asymmetric, such as `<=`, does reach it: the upward scan can
//! then leave its range, and the source's signed iterator distance ends the
//! recursion on such a range, which the port follows while it reads inside the
//! vector. No other input is refused, and no comparator makes the port panic.
//!
//! [`Error::InvalidValue`]: crate::Error::InvalidValue

use crate::math::libstdcxx;
use crate::{Error, Result};

/// Ranges up to this size are left to the insertion sort: libstdc++
/// `_S_threshold`.
const THRESHOLD: usize = 16;

/// The permutation the C++ Release build's `std::sort` applies to `len`
/// elements under the strict comparison `less(a, b)` on original positions.
///
/// Element `k` of the result is the original position of the element that ends
/// at position `k`. `less` is called with original positions, in the order the
/// libstdc++ introsort makes its comparisons; it must be deterministic, or
/// the port may stop a walk the C++ would continue.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the algorithm would read outside the
/// vector. Only a comparison that is not asymmetric can cause that; `<` on
/// floating-point keys, NaN included, never does (see the module
/// documentation).
pub fn source_sort_permutation(
    len: usize,
    less: impl FnMut(usize, usize) -> bool,
) -> Result<Vec<usize>> {
    let mut order = Vec::new();
    order
        .try_reserve_exact(len)
        .map_err(|_| Error::InvalidValue("cannot allocate the sort permutation".into()))?;
    order.extend(0..len);
    let mut sorter = Introsort {
        order,
        less,
        #[cfg(test)]
        heap_fallbacks: 0,
    };
    sorter.sort()?;
    Ok(sorter.order)
}

/// Sort `items` as the C++ Release build's `std::sort` does under `less`.
///
/// `items` is left unchanged when an error is returned. Every reachable error
/// is raised by [`source_sort_permutation`], before a single element has been
/// moved.
///
/// # Errors
///
/// As [`source_sort_permutation`]; and — unreachably, but refused rather than
/// left to panic or to spin — [`Error::InvalidValue`] if the permutation the
/// sort returned is not one, in which case `items` may hold a partially
/// applied permutation. `source_sort_permutation` starts from `0..len` and
/// only ever swaps within it, so that second error cannot be provoked from
/// here.
pub fn source_sort_by<T>(items: &mut [T], mut less: impl FnMut(&T, &T) -> bool) -> Result<()> {
    let mut order = source_sort_permutation(items.len(), |a, b| less(&items[a], &items[b]))?;
    apply_permutation(items, &mut order)
}

/// Sort `items` as `std::sort(items.rbegin(), items.rend())` does under
/// `less`: the reversed sequence is sorted ascending, which leaves the
/// sequence itself descending.
///
/// `items` is left unchanged when an error is returned, on the same terms as
/// [`source_sort_by`].
///
/// # Errors
///
/// As [`source_sort_by`].
pub fn source_sort_reversed_by<T>(
    items: &mut [T],
    mut less: impl FnMut(&T, &T) -> bool,
) -> Result<()> {
    let len = items.len();
    // Position `k` of the reversed view is element `len - 1 - k`.
    let mut order =
        source_sort_permutation(len, |a, b| less(&items[len - 1 - a], &items[len - 1 - b]))?;
    // The reversed view ends as `order`, so the sequence itself is its
    // reverse, read back through the same index map. Rewriting `order` in
    // place rather than collecting the composition keeps this function's
    // allocation to the one permutation `source_sort_permutation` reserved
    // fallibly.
    order.reverse();
    for position in order.iter_mut() {
        *position = len - 1 - *position;
    }
    apply_permutation(items, &mut order)
}

/// Reorder `items` so that position `k` holds the element that was at
/// `order[k]`.
///
/// Applied **in place**, cycle by cycle, with [`slice::swap`] and no auxiliary
/// buffer: `order` is the caller's own vector and is used as the scratch that
/// marks progress, every entry it has visited being set to its own index. The
/// buffered predecessor staged a `Vec<Option<T>>` — 16 bytes per `f64` — with
/// an infallible `collect`, so at `FileInfo::MAX_STATISTICS_VALUES` an
/// allocation failure **aborted the process**, which is the outcome lead
/// decision D1 refuses to reproduce. There is now nothing to fail.
///
/// `order` is left as the identity on success and is otherwise unspecified.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `order` is not a permutation of
/// `0..items.len()`, which is unreachable from either caller above:
/// [`source_sort_permutation`] starts from `0..len` and only ever swaps within
/// it. It is checked rather than assumed because the two ways a walk can leave
/// the cycle structure of a permutation are an index outside the slice, which
/// would make `swap` panic, and a second visit to a position already marked,
/// which would make the walk spin forever. Neither is acceptable on an input,
/// so both end the call instead. `items` may then hold a partially applied
/// permutation — the elements are all still there, in some order.
///
/// Completing without that error *proves* `order` was a permutation: each walk
/// marks every position it visits and refuses on a revisit, so the walks are
/// disjoint cycles, and every position no walk visited was already a fixed
/// point.
fn apply_permutation<T>(items: &mut [T], order: &mut [usize]) -> Result<()> {
    let len = items.len();
    if order.len() != len {
        return Err(Error::InvalidValue(
            "the sort permutation does not cover the range".into(),
        ));
    }
    for start in 0..len {
        if order[start] == start {
            // A fixed point, or a position an earlier walk already filled.
            continue;
        }
        // The cycle through `start` is `start -> order[start] -> ...`, and the
        // value position `k` must end with is the one at `order[k]`, so the
        // values rotate one step along it. Walking it with `swap` carries the
        // value that belongs at the *last* position of the cycle ahead of the
        // walk, and leaves it there when the cycle closes.
        let mut hole = start;
        loop {
            let source = order[hole];
            if source >= len {
                return Err(Error::InvalidValue(
                    "the sort permutation leaves the range".into(),
                ));
            }
            if source == hole && hole != start {
                return Err(Error::InvalidValue(
                    "the sort permutation visits a position twice".into(),
                ));
            }
            order[hole] = hole;
            if source == start {
                // The cycle is closed: `items[hole]` is already the element
                // that started at `start`, which is where it belongs.
                break;
            }
            items.swap(hole, source);
            hole = source;
        }
    }
    Ok(())
}

/// How [`source_stable_sort_permutation`] obtains `std::stable_sort`'s
/// temporary buffer.
pub enum TemporaryBuffer<'a> {
    /// Allocate it, halving the element count after every failed allocation
    /// as `std::get_temporary_buffer` does. What the port runs.
    Allocate,
    /// Grant a request of `count` elements exactly when `grant(count)` is
    /// true: the model of a C++ allocator that refuses some requests, with
    /// which the executed fallbacks are replayed. `grant` sees every request
    /// in order.
    Model(&'a mut dyn FnMut(usize) -> bool),
}

/// The permutation the C++ Release build's `std::stable_sort` applies to
/// `len` elements under the strict comparison `less(a, b)` on original
/// positions.
///
/// Element `k` of the result is the original position of the element that ends
/// at position `k`. `less` is called with original positions, in the order the
/// libstdc++ algorithm makes its comparisons; it must be deterministic. The
/// algorithm depends on how many of the `(len + 1) / 2` buffer elements it
/// asked for it obtained (module documentation), which `buffer` decides.
///
/// Every input is accepted: the algorithm stays inside the range whatever
/// `less` returns.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] only when the permutation itself cannot be
/// allocated.
pub fn source_stable_sort_permutation(
    len: usize,
    less: impl FnMut(usize, usize) -> bool,
    buffer: TemporaryBuffer<'_>,
) -> Result<Vec<usize>> {
    let mut order = Vec::new();
    order
        .try_reserve_exact(len)
        .map_err(|_| Error::InvalidValue("cannot allocate the sort permutation".into()))?;
    order.extend(0..len);
    if len == 0 {
        return Ok(order);
    }
    let requested = len.div_ceil(2);
    let buffer = temporary_buffer(requested, buffer);
    let granted = buffer.len();
    let mut sorter = StableSort {
        order,
        buffer,
        less,
    };
    if granted == requested {
        sorter.adaptive(0, granted, len);
    } else if granted == 0 {
        sorter.in_place(0, len);
    } else {
        sorter.adaptive_resize(0, len, granted);
    }
    Ok(sorter.order)
}

/// `std::get_temporary_buffer`: the first count in the sequence `requested`,
/// `(requested + 1) / 2`, and so on down to 1, that can be allocated, as a
/// buffer of that many slots, or an empty one. The request never exceeds
/// `PTRDIFF_MAX / sizeof(T)` here, since it is half of a vector's length.
fn temporary_buffer(requested: usize, policy: TemporaryBuffer<'_>) -> Vec<usize> {
    let mut count = requested;
    match policy {
        TemporaryBuffer::Allocate => {
            while count > 0 {
                let mut buffer = Vec::new();
                if buffer.try_reserve_exact(count).is_ok() {
                    buffer.resize(count, 0);
                    return buffer;
                }
                count = if count == 1 { 0 } else { count.div_ceil(2) };
            }
            Vec::new()
        }
        TemporaryBuffer::Model(grant) => {
            while count > 0 {
                if grant(count) {
                    return vec![0; count];
                }
                count = if count == 1 { 0 } else { count.div_ceil(2) };
            }
            Vec::new()
        }
    }
}

/// libstdc++'s `std::stable_sort` on a permutation. Positions index `order`
/// or `buffer`; `less` compares the original positions stored there.
struct StableSort<F> {
    order: Vec<usize>,
    buffer: Vec<usize>,
    less: F,
}

/// `std::__insertion_sort` (the one `std::sort` uses as well) on `order`.
///
/// An element less than the first moves to the front; every other one is
/// inserted from the right by `std::__unguarded_linear_insert`, whose walk
/// stops at the first element at the latest: that element is unchanged and
/// the inserted one was just found not to be less than it, so the bound below
/// never ends a walk the C++ would continue.
fn insertion_sort(order: &mut [usize], less: &mut impl FnMut(usize, usize) -> bool) {
    for i in 1..order.len() {
        if less(order[i], order[0]) {
            order[..=i].rotate_right(1);
        } else {
            let value = order[i];
            let mut hole = i;
            while hole > 0 && less(value, order[hole - 1]) {
                order[hole] = order[hole - 1];
                hole -= 1;
            }
            order[hole] = value;
        }
    }
}

/// `std::__move_merge`: merge `source[a..a_end]` and `source[b..b_end]` into
/// `target` from `out`, the second run's element first only when it is less.
/// Returns the position after the last element written.
#[allow(clippy::too_many_arguments)]
fn move_merge(
    source: &[usize],
    (mut a, a_end): (usize, usize),
    (mut b, b_end): (usize, usize),
    target: &mut [usize],
    mut out: usize,
    less: &mut impl FnMut(usize, usize) -> bool,
) -> usize {
    while a != a_end && b != b_end {
        if less(source[b], source[a]) {
            target[out] = source[b];
            b += 1;
        } else {
            target[out] = source[a];
            a += 1;
        }
        out += 1;
    }
    for &value in source[a..a_end].iter().chain(&source[b..b_end]) {
        target[out] = value;
        out += 1;
    }
    out
}

/// `std::__merge_sort_loop`: merge the runs of `step` elements of `source`
/// pairwise into `target`.
fn merge_sort_loop(
    source: &[usize],
    target: &mut [usize],
    step: usize,
    less: &mut impl FnMut(usize, usize) -> bool,
) {
    let last = source.len();
    let two_steps = 2 * step;
    let mut first = 0;
    let mut out = 0;
    while last - first >= two_steps {
        out = move_merge(
            source,
            (first, first + step),
            (first + step, first + two_steps),
            target,
            out,
            less,
        );
        first += two_steps;
    }
    let step = step.min(last - first);
    move_merge(
        source,
        (first, first + step),
        (first + step, last),
        target,
        out,
        less,
    );
}

/// `std::_S_chunk_size`.
const CHUNK: usize = 7;

impl<F: FnMut(usize, usize) -> bool> StableSort<F> {
    /// `std::__chunk_insertion_sort` with chunks of [`CHUNK`].
    fn chunk_insertion_sort(&mut self, mut first: usize, last: usize) {
        while last - first >= CHUNK {
            insertion_sort(&mut self.order[first..first + CHUNK], &mut self.less);
            first += CHUNK;
        }
        insertion_sort(&mut self.order[first..last], &mut self.less);
    }

    /// `std::__merge_sort_with_buffer`: `order[first..last]` fits the buffer.
    fn merge_sort_with_buffer(&mut self, first: usize, last: usize) {
        let len = last - first;
        let mut step = CHUNK;
        self.chunk_insertion_sort(first, last);
        let Self {
            order,
            buffer,
            less,
        } = self;
        while step < len {
            merge_sort_loop(&order[first..last], &mut buffer[..len], step, less);
            step *= 2;
            merge_sort_loop(&buffer[..len], &mut order[first..last], step, less);
            step *= 2;
        }
    }

    /// `std::__stable_sort_adaptive`: both halves fit the buffer.
    fn adaptive(&mut self, first: usize, middle: usize, last: usize) {
        self.merge_sort_with_buffer(first, middle);
        self.merge_sort_with_buffer(middle, last);
        self.merge_adaptive(first, middle, last, middle - first, last - middle);
    }

    /// `std::__stable_sort_adaptive_resize` with a buffer of `buffer_size`.
    fn adaptive_resize(&mut self, first: usize, last: usize, buffer_size: usize) {
        let len = (last - first).div_ceil(2);
        let middle = first + len;
        if len > buffer_size {
            self.adaptive_resize(first, middle, buffer_size);
            self.adaptive_resize(middle, last, buffer_size);
            self.merge_adaptive_resize(
                first,
                middle,
                last,
                middle - first,
                last - middle,
                buffer_size,
            );
        } else {
            self.adaptive(first, middle, last);
        }
    }

    /// `std::__merge_adaptive`: the shorter run, which fits the buffer, is
    /// copied there and merged back, forwards when it is the first run.
    fn merge_adaptive(
        &mut self,
        first: usize,
        middle: usize,
        last: usize,
        len1: usize,
        len2: usize,
    ) {
        let Self {
            order,
            buffer,
            less,
        } = self;
        if len1 <= len2 {
            buffer[..len1].copy_from_slice(&order[first..middle]);
            // `std::__move_merge_adaptive`.
            let (mut a, mut b, mut out) = (0, middle, first);
            while a != len1 && b != last {
                if less(order[b], buffer[a]) {
                    order[out] = order[b];
                    b += 1;
                } else {
                    order[out] = buffer[a];
                    a += 1;
                }
                out += 1;
            }
            order[out..out + (len1 - a)].copy_from_slice(&buffer[a..len1]);
        } else {
            buffer[..len2].copy_from_slice(&order[middle..last]);
            // `std::__move_merge_adaptive_backward`.
            if first == middle {
                order[last - len2..last].copy_from_slice(&buffer[..len2]);
                return;
            }
            if len2 == 0 {
                return;
            }
            let (mut a, mut b, mut out) = (middle - 1, len2 - 1, last);
            loop {
                if less(buffer[b], order[a]) {
                    out -= 1;
                    order[out] = order[a];
                    if a == first {
                        order[out - (b + 1)..out].copy_from_slice(&buffer[..=b]);
                        return;
                    }
                    a -= 1;
                } else {
                    out -= 1;
                    order[out] = buffer[b];
                    if b == 0 {
                        return;
                    }
                    b -= 1;
                }
            }
        }
    }

    /// The cut of `std::__merge_adaptive_resize` and
    /// `std::__merge_without_buffer`: halve the longer run and find where its
    /// middle element belongs in the other with `std::__lower_bound` (after the
    /// equal elements of the first run) or `std::__upper_bound`. Returns
    /// `(first_cut, second_cut, len11, len22)`.
    fn cut(
        &mut self,
        (first, middle, last): (usize, usize, usize),
        len1: usize,
        len2: usize,
    ) -> (usize, usize, usize, usize) {
        let Self { order, less, .. } = self;
        if len1 > len2 {
            let len11 = len1 / 2;
            let first_cut = first + len11;
            let value = order[first_cut];
            let len22 = libstdcxx::lower_bound(&order[middle..last], |&e| less(e, value));
            (first_cut, middle + len22, len11, len22)
        } else {
            let len22 = len2 / 2;
            let second_cut = middle + len22;
            let value = order[second_cut];
            let len11 = libstdcxx::upper_bound(&order[first..middle], |&e| less(value, e));
            (first + len11, second_cut, len11, len22)
        }
    }

    /// `std::__merge_adaptive_resize` with a buffer of `buffer_size`.
    fn merge_adaptive_resize(
        &mut self,
        first: usize,
        middle: usize,
        last: usize,
        len1: usize,
        len2: usize,
        buffer_size: usize,
    ) {
        if len1 <= buffer_size || len2 <= buffer_size {
            self.merge_adaptive(first, middle, last, len1, len2);
            return;
        }
        let (first_cut, second_cut, len11, len22) = self.cut((first, middle, last), len1, len2);
        let new_middle = self.rotate_adaptive(
            first_cut,
            middle,
            second_cut,
            len1 - len11,
            len22,
            buffer_size,
        );
        self.merge_adaptive_resize(first, first_cut, new_middle, len11, len22, buffer_size);
        self.merge_adaptive_resize(
            new_middle,
            second_cut,
            last,
            len1 - len11,
            len2 - len22,
            buffer_size,
        );
    }

    /// `std::__rotate_adaptive`: rotate `[first, middle)` behind
    /// `[middle, last)`, through the buffer when the shorter (second first)
    /// side fits it, else with `std::rotate`. Returns the new middle.
    fn rotate_adaptive(
        &mut self,
        first: usize,
        middle: usize,
        last: usize,
        len1: usize,
        len2: usize,
        buffer_size: usize,
    ) -> usize {
        let Self { order, buffer, .. } = self;
        if len1 > len2 && len2 <= buffer_size {
            if len2 == 0 {
                return first;
            }
            buffer[..len2].copy_from_slice(&order[middle..last]);
            order.copy_within(first..middle, last - len1);
            order[first..first + len2].copy_from_slice(&buffer[..len2]);
            first + len2
        } else if len1 <= buffer_size {
            if len1 == 0 {
                return last;
            }
            buffer[..len1].copy_from_slice(&order[first..middle]);
            order.copy_within(middle..last, first);
            order[last - len1..last].copy_from_slice(&buffer[..len1]);
            last - len1
        } else {
            order[first..last].rotate_left(middle - first);
            first + (last - middle)
        }
    }

    /// `std::__inplace_stable_sort`.
    fn in_place(&mut self, first: usize, last: usize) {
        if last - first < 15 {
            insertion_sort(&mut self.order[first..last], &mut self.less);
            return;
        }
        let middle = first + (last - first) / 2;
        self.in_place(first, middle);
        self.in_place(middle, last);
        self.merge_without_buffer(first, middle, last, middle - first, last - middle);
    }

    /// `std::__merge_without_buffer`.
    fn merge_without_buffer(
        &mut self,
        first: usize,
        middle: usize,
        last: usize,
        len1: usize,
        len2: usize,
    ) {
        if len1 == 0 || len2 == 0 {
            return;
        }
        if len1 + len2 == 2 {
            if (self.less)(self.order[middle], self.order[first]) {
                self.order.swap(first, middle);
            }
            return;
        }
        let (first_cut, second_cut, len11, len22) = self.cut((first, middle, last), len1, len2);
        // `std::rotate(first_cut, middle, second_cut)`.
        self.order[first_cut..second_cut].rotate_left(middle - first_cut);
        let new_middle = first_cut + (second_cut - middle);
        self.merge_without_buffer(first, first_cut, new_middle, len11, len22);
        self.merge_without_buffer(new_middle, second_cut, last, len1 - len11, len2 - len22);
    }
}

/// The refusal of a read outside the vector, which the C++ performs without
/// a check when the comparison is not asymmetric (module documentation).
fn unordered_read(where_: &str) -> Error {
    Error::InvalidValue(format!(
        "std::sort: the comparison is not a strict weak ordering, and the C++ \
         introsort reads {where_} here, which is undefined behaviour"
    ))
}

/// `std::__lg`: the position of the highest set bit, for `n >= 1`.
fn floor_log2(n: usize) -> usize {
    (usize::BITS - 1 - n.leading_zeros()) as usize
}

/// The libstdc++ introsort on a permutation. Every position below is an index
/// into `order`; `less` compares the original positions stored there.
struct Introsort<F> {
    order: Vec<usize>,
    less: F,
    /// How often the depth budget ran out and the heapsort took over.
    #[cfg(test)]
    heap_fallbacks: usize,
}

impl<F: FnMut(usize, usize) -> bool> Introsort<F> {
    /// `__comp(a, b)` on two positions.
    fn lt(&mut self, a: usize, b: usize) -> bool {
        let (x, y) = (self.order[a], self.order[b]);
        (self.less)(x, y)
    }

    /// `std::__sort`.
    fn sort(&mut self) -> Result<()> {
        let len = self.order.len();
        if len == 0 {
            return Ok(());
        }
        self.introsort_loop(0, len, floor_log2(len) * 2)?;
        self.final_insertion_sort(0, len)
    }

    /// `std::__introsort_loop`. The recursion is at most `2 * log2(len)` deep.
    ///
    /// The source compares the signed iterator distance `last - first` with
    /// the threshold. A partition whose cut lies beyond `last`, which only a
    /// comparator that is not asymmetric can produce, makes that distance
    /// negative in the recursive call, which then returns at once; the loop
    /// here does the same instead of wrapping.
    fn introsort_loop(
        &mut self,
        first: usize,
        mut last: usize,
        mut depth_limit: usize,
    ) -> Result<()> {
        while last > first && last - first > THRESHOLD {
            if depth_limit == 0 {
                #[cfg(test)]
                {
                    self.heap_fallbacks += 1;
                }
                // `__partial_sort(first, last, last)`: `__heap_select` with an
                // empty tail is `__make_heap`, followed by `__sort_heap`.
                self.make_heap(first, last);
                self.sort_heap(first, last);
                return Ok(());
            }
            depth_limit -= 1;
            let cut = self.unguarded_partition_pivot(first, last)?;
            self.introsort_loop(cut, last, depth_limit)?;
            last = cut;
        }
        Ok(())
    }

    /// `std::__unguarded_partition_pivot`.
    fn unguarded_partition_pivot(&mut self, first: usize, last: usize) -> Result<usize> {
        let mid = first + (last - first) / 2;
        self.move_median_to_first(first, first + 1, mid, last - 1);
        self.unguarded_partition(first + 1, last, first)
    }

    /// `std::__move_median_to_first`.
    fn move_median_to_first(&mut self, result: usize, a: usize, b: usize, c: usize) {
        let median = if self.lt(a, b) {
            if self.lt(b, c) {
                b
            } else if self.lt(a, c) {
                c
            } else {
                a
            }
        } else if self.lt(a, c) {
            a
        } else if self.lt(b, c) {
            c
        } else {
            b
        };
        self.order.swap(result, median);
    }

    /// `std::__unguarded_partition` around the element at `pivot`.
    ///
    /// Neither scan has a bound in the source. For a deterministic asymmetric
    /// comparison, `<` on floating-point keys included, both stay inside
    /// `[pivot, last)`: the median selection leaves an element that stops the
    /// upward scan, every swap moves such an element to the upper end, and the
    /// downward scan stops at the pivot itself, which is not less than itself.
    /// Another comparison can move the upward scan past `last`; the port
    /// follows it while it reads inside the vector, as the C++ does, and
    /// reports a scan that would leave the vector as the undefined read it is.
    fn unguarded_partition(
        &mut self,
        mut first: usize,
        mut last: usize,
        pivot: usize,
    ) -> Result<usize> {
        let len = self.order.len();
        loop {
            while self.lt(first, pivot) {
                first += 1;
                if first >= len {
                    return Err(unordered_read("past the last element"));
                }
            }
            last = last
                .checked_sub(1)
                .ok_or_else(|| unordered_read("before the first element"))?;
            while self.lt(pivot, last) {
                last = last
                    .checked_sub(1)
                    .ok_or_else(|| unordered_read("before the first element"))?;
            }
            if first >= last {
                return Ok(first);
            }
            self.order.swap(first, last);
            first += 1;
        }
    }

    /// `std::__final_insertion_sort`.
    fn final_insertion_sort(&mut self, first: usize, last: usize) -> Result<()> {
        if last - first > THRESHOLD {
            self.insertion_sort(first, first + THRESHOLD);
            for i in first + THRESHOLD..last {
                self.unguarded_linear_insert(i)?;
            }
            Ok(())
        } else {
            self.insertion_sort(first, last);
            Ok(())
        }
    }

    /// `std::__insertion_sort` on `[first, last)`.
    fn insertion_sort(&mut self, first: usize, last: usize) {
        insertion_sort(&mut self.order[first..last], &mut self.less);
    }

    /// `std::__unguarded_linear_insert`: move the element at `last` left while
    /// it is less than its left neighbour, with no lower bound.
    fn unguarded_linear_insert(&mut self, last: usize) -> Result<()> {
        let value = self.order[last];
        let mut hole = last;
        loop {
            let next = hole
                .checked_sub(1)
                .ok_or_else(|| unordered_read("before the first element"))?;
            if !(self.less)(value, self.order[next]) {
                break;
            }
            self.order[hole] = self.order[next];
            hole = next;
        }
        self.order[hole] = value;
        Ok(())
    }

    /// `std::__make_heap` over `[first, last)`.
    fn make_heap(&mut self, first: usize, last: usize) {
        let len = last - first;
        if len < 2 {
            return;
        }
        let mut parent = (len - 2) / 2;
        loop {
            let value = self.order[first + parent];
            self.adjust_heap(first, parent, len, value);
            if parent == 0 {
                return;
            }
            parent -= 1;
        }
    }

    /// `std::__sort_heap` over `[first, last)`.
    fn sort_heap(&mut self, first: usize, mut last: usize) {
        while last - first > 1 {
            last -= 1;
            // `__pop_heap(first, last, last)`.
            let value = self.order[last];
            self.order[last] = self.order[first];
            self.adjust_heap(first, 0, last - first, value);
        }
    }

    /// `std::__adjust_heap`: sift the hole at `hole` down to a leaf, then push
    /// `value` up from there.
    fn adjust_heap(&mut self, first: usize, hole: usize, len: usize, value: usize) {
        let top = hole;
        let mut hole = hole;
        let mut second = hole;
        while second < (len - 1) / 2 {
            second = 2 * (second + 1);
            if self.lt(first + second, first + (second - 1)) {
                second -= 1;
            }
            self.order[first + hole] = self.order[first + second];
            hole = second;
        }
        if len % 2 == 0 && second == (len - 2) / 2 {
            second = 2 * (second + 1);
            self.order[first + hole] = self.order[first + (second - 1)];
            hole = second - 1;
        }
        // `__push_heap(first, hole, top, value)` with the iterator-value
        // comparison. The signed `(hole - 1) / 2` truncates to 0 at hole 0,
        // where the loop condition fails anyway.
        let parent_of = |h: usize| h.saturating_sub(1) / 2;
        let mut parent = parent_of(hole);
        while hole > top && (self.less)(self.order[first + parent], value) {
            self.order[first + hole] = self.order[first + parent];
            hole = parent;
            parent = parent_of(hole);
        }
        self.order[first + hole] = value;
    }

    #[cfg(test)]
    fn heap_sort_only(&mut self) {
        let len = self.order.len();
        self.make_heap(0, len);
        self.sort_heap(0, len);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn keys_sorted(keys: &[f64]) -> Vec<usize> {
        source_sort_permutation(keys.len(), |a, b| keys[a] < keys[b]).unwrap()
    }

    #[test]
    fn a_total_order_is_sorted() {
        let keys: Vec<f64> = (0..200).map(|i| ((i * 7919) % 211) as f64).collect();
        let order = keys_sorted(&keys);
        for pair in order.windows(2) {
            assert!(keys[pair[0]] <= keys[pair[1]]);
        }
        let mut seen = order.clone();
        seen.sort_unstable();
        assert_eq!(seen, (0..200).collect::<Vec<_>>());
    }

    #[test]
    fn the_heap_fallback_sorts() {
        let keys: Vec<f64> = (0..100).map(|i| ((i * 37) % 17) as f64).collect();
        let mut sorter = Introsort {
            order: (0..keys.len()).collect(),
            less: |a: usize, b: usize| keys[a] < keys[b],
            heap_fallbacks: 0,
        };
        sorter.heap_sort_only();
        for pair in sorter.order.windows(2) {
            assert!(keys[pair[0]] <= keys[pair[1]]);
        }
    }

    /// McIlroy's adversary ("A Killer Adversary for Quicksort", 1999):
    /// values are fixed lazily so that every pivot is extreme. Returns keys
    /// that make the introsort take the same path again.
    pub(super) fn adversarial_keys(len: usize) -> (Vec<f64>, usize) {
        use std::cell::RefCell;
        const GAS: usize = usize::MAX;
        let values = RefCell::new(vec![GAS; len]);
        let solid = RefCell::new(0usize);
        let candidate = RefCell::new(0usize);
        let less = |x: usize, y: usize| {
            let mut values = values.borrow_mut();
            if values[x] == GAS && values[y] == GAS {
                let frozen = if x == *candidate.borrow() { x } else { y };
                let mut solid = solid.borrow_mut();
                values[frozen] = *solid;
                *solid += 1;
            }
            if values[x] == GAS {
                *candidate.borrow_mut() = x;
            } else if values[y] == GAS {
                *candidate.borrow_mut() = y;
            }
            values[x] < values[y]
        };
        let mut sorter = Introsort {
            order: (0..len).collect(),
            less,
            heap_fallbacks: 0,
        };
        sorter.sort().unwrap();
        let fallbacks = sorter.heap_fallbacks;
        drop(sorter);
        let values = values.into_inner();
        let rest = solid.into_inner();
        let keys = values
            .iter()
            .map(|&value| {
                if value == GAS {
                    rest as f64
                } else {
                    value as f64
                }
            })
            .collect();
        (keys, fallbacks)
    }

    #[test]
    fn an_adversarial_input_reaches_the_heap_fallback() {
        let file = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/feature_finder_picked_instrumentation/sort_keys.txt"
        ))
        .unwrap();
        let mut recorded = file.lines().skip(1);
        for len in [100, 500, 2000] {
            let (keys, fallbacks) = adversarial_keys(len);
            assert!(fallbacks > 0, "{len}");
            // The concrete keys drive the same path.
            let mut sorter = Introsort {
                order: (0..len).collect(),
                less: |a: usize, b: usize| keys[a] < keys[b],
                heap_fallbacks: 0,
            };
            sorter.sort().unwrap();
            assert_eq!(sorter.heap_fallbacks, fallbacks);
            for pair in sorter.order.windows(2) {
                assert!(keys[pair[0]] <= keys[pair[1]]);
            }
            // Lines 2 to 4 of the executed driver's key file are these keys.
            let text: Vec<String> = keys.iter().map(|k| format!("{k}")).collect();
            assert_eq!(recorded.next(), Some(text.join(" ").as_str()), "{len}");
        }
    }

    /// A deterministic generator for the comparator checks below.
    struct XorShift(u64);

    impl XorShift {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// Asymmetric comparisons never reach the out-of-bounds guard: `<` on
    /// keys with NaN of either sign, few distinct values and ties, and random
    /// asymmetric tournaments (module documentation). Lengths up to 300 cover
    /// the unguarded final insertion; the adversarial keys with NaN mixed in
    /// cover the top-level heap fallback.
    #[test]
    fn asymmetric_comparisons_never_reach_the_guard() {
        let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
        let pool = [0.0, -0.0, 1.0, 2.0, 3.0, f64::NAN, -f64::NAN, f64::INFINITY];
        for round in 0..4000 {
            let len = 1 + rng.below(300) as usize;
            let spread = 1 + rng.below(pool.len() as u64) as usize;
            let keys: Vec<f64> = (0..len)
                .map(|_| {
                    if rng.below(4) == 0 {
                        pool[rng.below(pool.len() as u64) as usize]
                    } else {
                        pool[rng.below(spread as u64) as usize] + rng.below(3) as f64
                    }
                })
                .collect();
            let order = source_sort_permutation(len, |a, b| keys[a] < keys[b])
                .unwrap_or_else(|error| panic!("round {round}: {error}"));
            let mut seen = order.clone();
            seen.sort_unstable();
            assert_eq!(seen, (0..len).collect::<Vec<_>>());
        }
        for round in 0..1000 {
            let len = 1 + rng.below(200) as usize;
            // A random tournament: exactly one of `a < b` and `b < a`, or
            // neither, for every pair; irreflexive.
            let relation: Vec<u8> = (0..len * len).map(|_| rng.below(3) as u8).collect();
            let less = |a: usize, b: usize| {
                if a == b {
                    return false;
                }
                let (low, high) = (a.min(b), a.max(b));
                match relation[low * len + high] {
                    0 => a == low,
                    1 => a == high,
                    _ => false,
                }
            };
            source_sort_permutation(len, less)
                .unwrap_or_else(|error| panic!("tournament {round}: {error}"));
        }
        for len in [100, 500, 2000] {
            let (mut keys, fallbacks) = adversarial_keys(len);
            assert!(fallbacks > 0);
            for (index, key) in keys.iter_mut().enumerate() {
                if index % 7 == 3 {
                    *key = f64::NAN;
                }
            }
            source_sort_permutation(len, |a, b| keys[a] < keys[b]).unwrap();
        }
    }

    /// Comparisons that are not asymmetric may make the C++ read outside the
    /// vector, which the port refuses, but never make the port panic: a
    /// partition whose cut lies beyond its range ends the recursion there, as
    /// the source's signed iterator distance does (verifier report of the
    /// combined fix round, `<=` with few distinct values, a symmetric relation
    /// and a nondeterministic one).
    #[test]
    fn other_comparisons_are_refused_or_sorted_without_a_panic() {
        let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
        let mut refused = 0;
        for _ in 0..3000 {
            let len = 1 + rng.below(120) as usize;
            let keys: Vec<u64> = (0..len).map(|_| rng.below(4)).collect();
            if source_sort_permutation(len, |a, b| keys[a] <= keys[b]).is_err() {
                refused += 1;
            }
            let relation: Vec<bool> = (0..len * len).map(|_| rng.below(2) == 0).collect();
            let symmetric = |a: usize, b: usize| relation[a.min(b) * len + a.max(b)];
            let _ = source_sort_permutation(len, symmetric);
            let mut coin = XorShift(rng.next() | 1);
            let _ = source_sort_permutation(len, |_, _| coin.below(2) == 0);
            let mut coin = XorShift(rng.next() | 1);
            let _ = source_stable_sort_permutation(
                len,
                |_, _| coin.below(2) == 0,
                TemporaryBuffer::Allocate,
            )
            .unwrap();
        }
        assert!(refused > 0, "the `<=` generator reaches the guard");
    }

    #[test]
    fn floor_log2_matches_bit_width_minus_one() {
        assert_eq!(floor_log2(1), 0);
        assert_eq!(floor_log2(2), 1);
        assert_eq!(floor_log2(3), 1);
        assert_eq!(floor_log2(16), 4);
        assert_eq!(floor_log2(17), 4);
    }

    #[test]
    fn the_reversed_sort_is_descending() {
        let mut items = vec![3.0, 1.0, 2.0, 5.0, 4.0];
        source_sort_reversed_by(&mut items, |a, b| a < b).unwrap();
        assert_eq!(items, vec![5.0, 4.0, 3.0, 2.0, 1.0]);
    }

    /// The buffered `apply_permutation` this module replaced: `order[k]` is
    /// read out of a `Vec<Option<T>>` staged from `items`. Kept here, in test
    /// code only, as the reference the in-place walk is differentially
    /// compared against — it is the behaviour that was tier-1 validated, and
    /// the only way to prove the replacement is a replacement.
    fn apply_permutation_buffered<T>(items: &mut Vec<T>, order: &[usize]) {
        let mut slots: Vec<Option<T>> = items.drain(..).map(Some).collect();
        for &position in order {
            if let Some(item) = slots.get_mut(position).and_then(Option::take) {
                items.push(item);
            }
        }
    }

    /// A deterministic 64-bit xorshift, so the property test is the same
    /// sequence on every host and in every run.
    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    // Native, finding F2: `apply_permutation` is applied in place, following
    // cycles with `slice::swap` and no auxiliary buffer, where it used to
    // stage a `Vec<Option<T>>` of 16 bytes per `f64` with an infallible
    // `collect`. This asserts the two agree element for element over random
    // permutations of every length that matters, including the empty one, the
    // identity, the full reversal and a rotation — which between them cover a
    // single long cycle, `n` fixed points and everything in between.
    #[test]
    fn the_in_place_permutation_equals_the_buffered_one() {
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut lengths: Vec<usize> = (0..40).collect();
        lengths.extend([64, 127, 128, 255, 1000]);
        for &len in &lengths {
            // Four shapes per length, then random ones.
            let mut orders: Vec<Vec<usize>> = vec![
                (0..len).collect(),
                (0..len).rev().collect(),
                (0..len).map(|k| (k + 1) % len.max(1)).collect(),
                (0..len).map(|k| (k + len / 2) % len.max(1)).collect(),
            ];
            for _ in 0..8 {
                // Fisher-Yates from the same stream.
                let mut order: Vec<usize> = (0..len).collect();
                for i in (1..len).rev() {
                    let j = (xorshift(&mut state) % (i as u64 + 1)) as usize;
                    order.swap(i, j);
                }
                orders.push(order);
            }
            for order in orders {
                let items: Vec<usize> = (0..len).map(|k| k * 7 + 1).collect();

                let mut buffered = items.clone();
                apply_permutation_buffered(&mut buffered, &order);

                let mut in_place = items.clone();
                let mut scratch = order.clone();
                apply_permutation(&mut in_place, &mut scratch).unwrap();

                assert_eq!(in_place, buffered, "len {len}, order {order:?}");
                // The permutation applied is the one asked for.
                for (k, &source) in order.iter().enumerate() {
                    assert_eq!(in_place[k], items[source]);
                }
                // `order` is left as the identity.
                assert!(scratch.iter().enumerate().all(|(k, &v)| k == v));
            }
        }
    }

    // Native, finding F2: the walk refuses rather than panicking on an index
    // outside the slice or spinning forever on a repeated one. Neither is
    // reachable from `source_sort_permutation`, which starts from `0..len` and
    // only swaps within it, but this module may not panic or hang on any
    // input, so both are checked.
    #[test]
    fn a_non_permutation_is_refused_rather_than_panicking_or_spinning() {
        // An index outside the slice.
        let mut items = vec![1.0, 2.0, 3.0];
        let mut order = vec![0, 3, 2];
        assert!(apply_permutation(&mut items, &mut order).is_err());

        // A repeated index: position 1 is visited twice.
        let mut items = vec![1.0, 2.0];
        let mut order = vec![1, 1];
        assert!(apply_permutation(&mut items, &mut order).is_err());

        // A longer one whose walk leaves the cycle structure.
        let mut items = vec![1.0, 2.0, 3.0, 4.0];
        let mut order = vec![1, 2, 0, 0];
        assert!(apply_permutation(&mut items, &mut order).is_err());

        // A length mismatch.
        let mut items = vec![1.0, 2.0, 3.0];
        let mut order = vec![0, 1];
        assert!(apply_permutation(&mut items, &mut order).is_err());
    }
}
