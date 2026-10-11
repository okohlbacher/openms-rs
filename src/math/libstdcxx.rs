// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The libstdc++ algorithms the source calls on floating-point keys, as the
//! conda-forge GCC 14.4.0 headers the Linux x86_64 Release build was compiled
//! with implement them (`bits/stl_algobase.h`, `bits/stl_algo.h`).
//!
//! The standard requires the keys of a binary search to be partitioned and
//! the keys of a sort to be strictly weakly ordered; a NaN key breaks both
//! whenever the other keys are not all equivalent, and the standard then
//! leaves the result undefined. These functions give the result the library
//! code computes anyway, comparison by comparison: neither algorithm reads
//! outside its range whatever the comparisons return.
//!
//! These are not specific to any one ported header: `std::lower_bound` and
//! `std::upper_bound` decide *positions* in the source, and a position the
//! Release build computes is part of the result every port of it has to
//! reproduce. The module therefore lives under `math` beside
//! [`crate::math::x86_64`] and [`crate::math::source_sort`], and was promoted
//! out of `analysis::feature_finder_picked::scoring` when shared math needed
//! it; not one line of the algorithms changed in the move.
//!
//! **Which citations here are resolvable.** `bits/stl_algo.h` is retained at
//! `../oracle/sne-completion/libstdcxx/stl_algo.h` (sha256 `0598c5b1…`), so
//! [`upper_bound`]'s line range can be checked against it.
//! `bits/stl_algobase.h` is **not** retained: its sha256 is pinned
//! (`0ec2358c…`, `tests/data/feature_finder_picked_provenance.json`), which is
//! what makes a toolchain change detectable, but no copy of the file is kept,
//! so [`lower_bound`]'s line range cannot be resolved from any artefact in this
//! repository and is given below as the header's own numbering rather than as a
//! checked citation. When the header is retained — at
//! `../oracle/sne-completion/libstdcxx/stl_algobase.h`, beside the two that
//! already are — that line range becomes checkable and this paragraph goes.

/// `std::__lower_bound` (`bits/stl_algobase.h:1491-1514`, the header's own
/// numbering; the header is not retained, see the module documentation): the
/// first position whose element does not satisfy `less_than_value`, found by
/// halving.
pub(crate) fn lower_bound<T>(items: &[T], mut less_than_value: impl FnMut(&T) -> bool) -> usize {
    let mut first = 0usize;
    let mut len = items.len();
    while len > 0 {
        let half = len >> 1;
        // `first + len <= items.len()` holds throughout, so `middle` is
        // in range.
        let middle = first + half;
        if less_than_value(&items[middle]) {
            first = middle + 1;
            len = len - half - 1;
        } else {
            len = half;
        }
    }
    first
}

/// `std::__upper_bound` (`bits/stl_algo.h:1980-2003`, checkable against the
/// retained header): the first position
/// whose element satisfies `value_less_than`, found by halving.
pub(crate) fn upper_bound<T>(items: &[T], mut value_less_than: impl FnMut(&T) -> bool) -> usize {
    let mut first = 0usize;
    let mut len = items.len();
    while len > 0 {
        let half = len >> 1;
        let middle = first + half;
        if value_less_than(&items[middle]) {
            len = half;
        } else {
            first = middle + 1;
            len = len - half - 1;
        }
    }
    first
}

/// `std::is_sorted` (`std::is_sorted_until`): no element is `less` than
/// its predecessor.
pub(crate) fn is_sorted_by<T>(items: &[T], mut less: impl FnMut(&T, &T) -> bool) -> bool {
    items.windows(2).all(|pair| !less(&pair[1], &pair[0]))
}

/// `std::nth_element(first, first + nth, last, less)` as the Release
/// toolchain's libstdc++ (GCC 14.4.0, `bits/stl_algo.h` and `bits/stl_heap.h`)
/// implements it, on any `Copy` element: introselect with a depth limit of
/// `2 * floor(log2(n))`, median-of-three pivots, the `__heap_select` fallback
/// and a final insertion sort, every comparison and move in the library's
/// order, so that the permutation it leaves is the one the source's later
/// steps (tie-breaks, kd-tree construction, quantile reads) see.
///
/// An `nth` of `items.len()` or more returns at once, as the library does for
/// `nth == last`. Every scan carries an explicit bound that a strict weak
/// order never reaches; with a comparator that is not one (a NaN key under
/// `<`) the algorithm still reads and writes inside the range only, and the
/// result is the one the library's code computes comparison by comparison.
///
/// This is the one transcription: the `f32` form in
/// `processing::noise_estimation` (checked bit for bit against the Release
/// toolchain on 390 vectors, `tests/data/signal_to_noise/oracle.tsv`), the
/// pose-clustering superimposers' form and the kd-tree's form were three
/// copies of it before integration.
pub(crate) fn nth_element_by<T: Copy>(items: &mut [T], nth: usize, less: impl Fn(&T, &T) -> bool) {
    let n = items.len();
    if n == 0 || nth >= n {
        return;
    }
    // `std::__lg(last - first) * 2`; n > 0.
    let lg = (usize::BITS - 1 - n.leading_zeros()) as usize;
    let (mut first, mut last, mut depth) = (0usize, n, lg * 2);
    while last - first > 3 {
        if depth == 0 {
            heap_select(items, first, nth + 1, last, &less);
            // Place the nth largest element in its final position.
            items.swap(first, nth);
            return;
        }
        depth -= 1;
        // `__unguarded_partition_pivot`.
        let mid = first + (last - first) / 2;
        move_median_to_first(items, first, first + 1, mid, last - 1, &less);
        let cut = unguarded_partition(items, first + 1, last, first, &less);
        if cut <= nth {
            first = cut;
        } else {
            last = cut;
        }
    }
    insertion_sort(items, first, last, &less);
}

/// `__move_median_to_first(result, a, b, c)`: swaps the median of `v[a]`,
/// `v[b]` and `v[c]` under `less` into `result`.
fn move_median_to_first<T: Copy>(
    v: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    less: &impl Fn(&T, &T) -> bool,
) {
    let median = if less(&v[a], &v[b]) {
        if less(&v[b], &v[c]) {
            b
        } else if less(&v[a], &v[c]) {
            c
        } else {
            a
        }
    } else if less(&v[a], &v[c]) {
        a
    } else if less(&v[b], &v[c]) {
        c
    } else {
        b
    };
    v.swap(result, median);
}

/// `__unguarded_partition(first, last, pivot)`.
fn unguarded_partition<T: Copy>(
    v: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: usize,
    less: &impl Fn(&T, &T) -> bool,
) -> usize {
    let end = last;
    loop {
        while first < end && less(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while last > pivot && less(&v[pivot], &v[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

/// `__insertion_sort(first, last)` with `__unguarded_linear_insert`.
fn insertion_sort<T: Copy>(v: &mut [T], first: usize, last: usize, less: &impl Fn(&T, &T) -> bool) {
    if first == last {
        return;
    }
    for i in first + 1..last {
        let value = v[i];
        if less(&value, &v[first]) {
            // `_GLIBCXX_MOVE_BACKWARD3(first, i, i + 1)`.
            v.copy_within(first..i, first + 1);
            v[first] = value;
        } else {
            // `__unguarded_linear_insert`, stopped at `first` by the bound.
            let mut hole = i;
            let mut next = i - 1;
            while less(&value, &v[next]) {
                v[hole] = v[next];
                hole = next;
                if next == first {
                    break;
                }
                next -= 1;
            }
            v[hole] = value;
        }
    }
}

/// `__heap_select(first, middle, last)`.
fn heap_select<T: Copy>(
    v: &mut [T],
    first: usize,
    middle: usize,
    last: usize,
    less: &impl Fn(&T, &T) -> bool,
) {
    make_heap(v, first, middle, less);
    for i in middle..last {
        if less(&v[i], &v[first]) {
            // `__pop_heap(first, middle, i)`.
            let value = v[i];
            v[i] = v[first];
            adjust_heap(v, first, 0, middle - first, value, less);
        }
    }
}

/// `__make_heap(first, last)`.
fn make_heap<T: Copy>(v: &mut [T], first: usize, last: usize, less: &impl Fn(&T, &T) -> bool) {
    if last - first < 2 {
        return;
    }
    let len = last - first;
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[first + parent];
        adjust_heap(v, first, parent, len, value, less);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

/// `__adjust_heap(first, hole, len, value)` followed by its `__push_heap`;
/// `(hole - 1) / 2` truncates toward zero in the signed source, so a hole of
/// `0` has parent `0`.
fn adjust_heap<T: Copy>(
    v: &mut [T],
    first: usize,
    mut hole: usize,
    len: usize,
    value: T,
    less: &impl Fn(&T, &T) -> bool,
) {
    let top = hole;
    let mut second = hole;
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if less(&v[first + second], &v[first + second - 1]) {
            second -= 1;
        }
        v[first + hole] = v[first + second];
        hole = second;
    }
    if len & 1 == 0 && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[first + hole] = v[first + second - 1];
        hole = second - 1;
    }
    let mut parent = hole.saturating_sub(1) / 2;
    while hole > top && less(&v[first + parent], &value) {
        v[first + hole] = v[first + parent];
        hole = parent;
        parent = hole.saturating_sub(1) / 2;
    }
    v[first + hole] = value;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tier 4: `nth_element` leaves the `nth` element in its sorted
    /// position with no greater element before and no smaller one after, for
    /// every position of several inputs, including the heap-select path
    /// (reached by the adversarial orders). The permutation itself is
    /// checked against the Release toolchain by
    /// `processing::noise_estimation` (390 vectors, 30 of them through
    /// `__heap_select`).
    #[test]
    fn nth_element_partitions() {
        let inputs: Vec<Vec<i32>> = vec![
            vec![5, 1, 4, 2, 3],
            vec![1, 1, 1, 1, 1, 1, 1],
            (0..50).rev().collect(),
            (0..50).map(|i| (i * 37) % 23).collect(),
            vec![3, 1, 2],
        ];
        for input in inputs {
            let mut sorted = input.clone();
            sorted.sort();
            for nth in 0..input.len() {
                let mut v = input.clone();
                nth_element_by(&mut v, nth, |a: &i32, b: &i32| a < b);
                assert_eq!(v[nth], sorted[nth]);
                assert!(v[..nth].iter().all(|x| *x <= v[nth]));
                assert!(v[nth + 1..].iter().all(|x| *x >= v[nth]));
            }
        }
    }

    #[test]
    fn heap_select_keeps_smallest() {
        let mut v: Vec<i32> = vec![9, 3, 7, 1, 8, 2, 6, 4, 5, 0];
        let n = v.len();
        heap_select(&mut v, 0, 4, n, &|a: &i32, b: &i32| a < b);
        let mut head = v[..4].to_vec();
        head.sort();
        assert_eq!(head, vec![0, 1, 2, 3]);
        // A max-heap on the first four: the root is the largest of them.
        assert_eq!(v[0], 3);
    }

    /// A comparator that is not a strict weak order (NaN keys under `<`)
    /// still terminates inside the range.
    #[test]
    fn nan_keys_stay_inside_the_range() {
        let mut v = vec![f32::NAN, 1.0, f32::NAN, 0.5, 2.0, f32::NAN, 3.0, 0.1];
        let n = v.len();
        for nth in 0..=n {
            let mut w = v.clone();
            nth_element_by(&mut w, nth, |a: &f32, b: &f32| a < b);
            assert_eq!(w.len(), n);
        }
        v.clear();
        nth_element_by(&mut v, 0, |a: &f32, b: &f32| a < b);
    }
}
