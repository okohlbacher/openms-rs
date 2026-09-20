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
