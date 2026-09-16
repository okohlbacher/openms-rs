// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The order in which the C++ Release build's `std::sort` leaves equal
//! elements, for the three sorts of the picked feature finder
//! (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`) whose ties are observable.
//!
//! `FeatureFinderAlgorithmPicked::run_` sorts with `std::sort` three times: the
//! seeds by descending intensity (`FeatureFinderAlgorithmPicked.cpp:548`), and
//! the feature map by m/z (`FeatureMap::sortByMZ`, `:866`) and by descending
//! intensity (`FeatureMap::sortByIntensity(true)`, `:991`). `std::sort` is not
//! stable, and the standard leaves the order of equal elements unspecified.
//! Which of two equal features comes first decides which one step 4 keeps, and
//! a reused instance or a caller's non-empty map produces such ties routinely:
//! running twice on one input doubles every feature.
//!
//! The port follows the Linux x86_64 Release build, whose `std::sort` is the
//! introsort of the conda-forge GCC 14.4.0 `libstdc++` headers it was compiled
//! with (`bits/stl_algo.h`, sha256 `0598c5b1...`, and `bits/stl_heap.h`, sha256
//! `f18f83b2...`): a quicksort that moves the median of the first, middle and
//! last element to the front and partitions around it, recursing into the
//! upper part, down to ranges of at most 16 elements or down to a recursion
//! budget of `2 * floor(log2(n))`, where a heapsort takes over; and a final
//! insertion sort over the whole range, whose first 16 elements are inserted
//! with a lower-bound check and the rest without one. This module reproduces
//! that algorithm comparison by comparison and move by move on a permutation,
//! so equal elements land where the executed C++ puts them. The oracle driver
//! `../oracle/ffap-instr-completion/drivers/ffap_instr_driver.cpp` (mode `sort`)
//! compares it with `FeatureMap::sortByMZ` and `sortByIntensity(true)` of the
//! Release `libOpenMS.so` on 130 generated inputs.
//!
//! # Undefined behaviour
//!
//! `std::sort` requires a strict weak ordering; `<` on floating-point keys that
//! include NaN is not one. The libstdc++ algorithm still runs on such keys, and
//! its partition and insertion loops have no bound checks: they rely on the
//! ordering to stop. This port follows the algorithm for such keys too, since
//! every comparison and move stays defined while it reads inside the vector,
//! and returns [`Error::InvalidValue`] at exactly the step where the C++ would
//! read outside it, which is undefined behaviour with no reproducible result.
//! No other input is refused.
//!
//! [`Error::InvalidValue`]: crate::Error::InvalidValue

use crate::{Error, Result};

/// Ranges up to this size are left to the insertion sort: libstdc++
/// `_S_threshold`.
const THRESHOLD: usize = 16;

/// The permutation the C++ Release build's `std::sort` applies to `len`
/// elements under the strict comparison `less(a, b)` on original positions.
///
/// Element `k` of the result is the original position of the element that ends
/// at position `k`. `less` is called with original positions, in the order the
/// libstdc++ introsort makes its comparisons; it must be deterministic.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the algorithm would read outside the
/// vector, which only an ordering that is not a strict weak ordering can cause
/// (see the module documentation).
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
/// `items` is left unchanged when an error is returned.
///
/// # Errors
///
/// As [`source_sort_permutation`].
pub fn source_sort_by<T>(items: &mut Vec<T>, mut less: impl FnMut(&T, &T) -> bool) -> Result<()> {
    let order = source_sort_permutation(items.len(), |a, b| less(&items[a], &items[b]))?;
    apply_permutation(items, &order);
    Ok(())
}

/// Sort `items` as `std::sort(items.rbegin(), items.rend())` does under
/// `less`: the reversed sequence is sorted ascending, which leaves the
/// sequence itself descending.
///
/// `items` is left unchanged when an error is returned.
///
/// # Errors
///
/// As [`source_sort_permutation`].
pub fn source_sort_reversed_by<T>(
    items: &mut Vec<T>,
    mut less: impl FnMut(&T, &T) -> bool,
) -> Result<()> {
    let len = items.len();
    // Position `k` of the reversed view is element `len - 1 - k`.
    let order =
        source_sort_permutation(len, |a, b| less(&items[len - 1 - a], &items[len - 1 - b]))?;
    // The reversed view ends as `order`, so the sequence is its reverse.
    let original: Vec<usize> = order.iter().rev().map(|&k| len - 1 - k).collect();
    apply_permutation(items, &original);
    Ok(())
}

/// Reorder `items` so that position `k` holds the element that was at
/// `order[k]`; `order` is a permutation of `0..items.len()`.
fn apply_permutation<T>(items: &mut Vec<T>, order: &[usize]) {
    let mut slots: Vec<Option<T>> = items.drain(..).map(Some).collect();
    for &position in order {
        if let Some(item) = slots.get_mut(position).and_then(Option::take) {
            items.push(item);
        }
    }
}

/// The refusal of a read outside the vector, which the C++ performs without
/// a check when the keys are not strictly weakly ordered.
fn unordered_read(where_: &str) -> Error {
    Error::InvalidValue(format!(
        "std::sort: the keys are not strictly weakly ordered (a NaN key), and the C++ \
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
    fn introsort_loop(
        &mut self,
        first: usize,
        mut last: usize,
        mut depth_limit: usize,
    ) -> Result<()> {
        while last - first > THRESHOLD {
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
    /// Neither scan has a bound in the source. For any deterministic `<` on
    /// floating-point keys both stay inside the vector: the median selection
    /// leaves an element that stops the upward scan, every swap moves such an
    /// element to the upper end, and the downward scan stops at the pivot
    /// itself, which is not less than itself. A scan that would leave the
    /// vector anyway is reported as the undefined read it is.
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
            self.insertion_sort(first, first + THRESHOLD)?;
            for i in first + THRESHOLD..last {
                self.unguarded_linear_insert(i)?;
            }
            Ok(())
        } else {
            self.insertion_sort(first, last)
        }
    }

    /// `std::__insertion_sort`: an element less than the first moves to the
    /// front; every other one is inserted from the right.
    fn insertion_sort(&mut self, first: usize, last: usize) -> Result<()> {
        if first == last {
            return Ok(());
        }
        for i in first + 1..last {
            if self.lt(i, first) {
                self.order[first..=i].rotate_right(1);
            } else {
                // The walk cannot pass `first`: the element is not less than it.
                self.unguarded_linear_insert(i)?;
            }
        }
        Ok(())
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
}
