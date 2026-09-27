// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The superimposer interface and the histogram machinery the pose-clustering
//! superimposers share.
//!
//! Ports `ANALYSIS/MAPMATCHING/BaseSuperimposer.h`: the abstract base class
//! becomes the [`BaseSuperimposer`](crate::analysis::mapmatching::superimposer::BaseSuperimposer) trait. The helpers below are the parts of
//! `Math::LinearInterpolation` (the voting histograms), of
//! `MorphologicalFilter::filterRange` (the top-hat baseline removal on
//! `double` values) and of the two superimposers' shared statistics that both
//! `PoseClusteringShiftSuperimposer.cpp` and
//! `PoseClusteringAffineSuperimposer.cpp` repeat. See
//! `docs/MAP_ALIGNMENT_SUPPORT.md`.

use crate::analysis::transformations::TransformationDescription;
use crate::concept::progress_logger::ProgressReporter;
use crate::kernel::ConsensusMap;
use crate::math::basic_statistics::BasicStatistics;
use crate::param::{Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use std::collections::VecDeque;

/// A superimposer: estimates the transformation that maps the elements of a
/// scene map as near as possible onto those of a model map (source
/// `BaseSuperimposer`, a `DefaultParamHandler` named `BaseSuperimposer` and a
/// `ProgressLogger`).
///
/// The source's `run` fills a `TransformationDescription` out-parameter; here
/// it is returned, together with the lines the source writes to its warning
/// log, in a [`Superimposition`].
pub trait BaseSuperimposer {
    /// Estimate the transformation from `scene` onto `model`, reporting
    /// progress as the source's `ProgressLogger` base does (source `run`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the input maps are invalid (the
    /// source's `Exception::IllegalArgument`), and as the implementation
    /// documents.
    fn superimpose(
        &self,
        model: &ConsensusMap,
        scene: &ConsensusMap,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Superimposition>;

    /// [`BaseSuperimposer::superimpose`] without progress, returning only the
    /// transformation.
    ///
    /// # Errors
    ///
    /// As [`BaseSuperimposer::superimpose`].
    fn run(&self, model: &ConsensusMap, scene: &ConsensusMap) -> Result<TransformationDescription> {
        Ok(self
            .superimpose(model, scene, &mut ProgressReporter::silent())?
            .transformation)
    }
}

/// What a superimposer returns.
#[derive(Clone, Debug)]
pub struct Superimposition {
    /// The estimated transformation from the scene onto the model.
    pub transformation: TransformationDescription,
    /// The lines the source writes to `OPENMS_LOG_WARN`, without line ends,
    /// in order.
    pub warnings: Vec<String>,
}

/// Number of samples of the top-hat structuring element that removes the
/// histogram baseline (the sources' `struc_elem_length_datapoints`).
pub(super) const TOPHAT_LENGTH: usize = 21;
/// Slope divisor of the noise-level line through the sorted histogram (the
/// sources' `scaling_histogram_crossing_slope`).
pub(super) const CROSSING_SLOPE: f64 = 3.0;
/// Multiplier of the standard deviation in the outlier cut-off (the sources'
/// `scaling_cutoff_stdev_multiplier`).
pub(super) const CUTOFF_STDEV_MULTIPLIER: f64 = 1.5;
/// Iterations of the mean and standard-deviation cut-off (the sources'
/// `loops_mean_stdev_cutoff`).
pub(super) const CUTOFF_LOOPS: usize = 3;
/// A window of elements with similar m/z gets unit weight; this is subtracted
/// from the weight of each element in it, so that windows with ten or more
/// elements are skipped (the sources' `winlength_factor_baseline`).
pub(super) const WINDOW_LENGTH_BASELINE: f64 = 0.1;

/// A voting histogram over an affinely mapped key axis: the part of
/// `Math::LinearInterpolation<double, double>` the superimposers use.
#[derive(Clone, Debug)]
pub(super) struct Histogram {
    pub(super) data: Vec<f64>,
    scale: f64,
    offset: f64,
}

impl Histogram {
    /// `data.resize(buckets)` and `setMapping(scale, inside, outside)`: key
    /// `outside` maps onto bucket `inside`, consecutive buckets are `scale`
    /// apart.
    pub(super) fn new(buckets: usize, scale: f64, inside: f64, outside: f64) -> Self {
        Self {
            data: vec![0.0; buckets],
            scale,
            offset: outside - scale * inside,
        }
    }

    /// `key2index`: the fractional bucket of a key; 0 for a zero scale.
    pub(super) fn key_to_index(&self, key: f64) -> f64 {
        if self.scale != 0.0 {
            (key - self.offset) / self.scale
        } else {
            0.0
        }
    }

    /// `index2key`: the key of a (fractional) bucket.
    pub(super) fn index_to_key(&self, index: f64) -> f64 {
        index * self.scale + self.offset
    }

    /// `addValue`: split `value` between the two buckets around `key` in
    /// proportion to its fractional position; a key more than one bucket
    /// outside the range is dropped, one within a bucket of either end
    /// contributes only to that end.
    ///
    /// The source converts the integer part with `DiffType(left_key)`, which
    /// is undefined for keys beyond the `ptrdiff_t` range; the saturating
    /// conversion here puts such keys outside the range, where they are
    /// dropped.
    pub(super) fn add_value(&mut self, key: f64, value: f64) {
        let position = self.key_to_index(key);
        // std::modf: integral part toward zero, fraction with the same sign.
        let whole = position.trunc();
        let fraction = position - whole;
        let left = whole as i64;
        if position < 0.0 {
            if left == 0 {
                self.data[0] += (1.0 + fraction) * value;
            }
            return;
        }
        let back = self.data.len() as i64 - 1;
        if left >= back {
            if left == back {
                self.data[left as usize] += (1.0 - fraction) * value;
            }
            return;
        }
        self.data[left as usize + 1] += fraction * value;
        self.data[left as usize] += (1.0 - fraction) * value;
    }

    /// Remove the baseline with a top-hat filter of [`TOPHAT_LENGTH`] samples
    /// (source `MorphologicalFilter::filterRange` with `tophat` and
    /// `DataPoints`).
    pub(super) fn tophat(&mut self) {
        self.data = tophat(&self.data, TOPHAT_LENGTH);
    }

    /// Estimate the noise level of the filtered histogram and set the buckets
    /// below it to zero; returns the cut-off.
    ///
    /// The histogram's values, sorted in descending order, are walked from the
    /// second while they lie on or above the line through the first bucket's
    /// value with slope `(last - first) / n / 3`; the cut-off is the last
    /// value on or above it. A zero slope gives a cut-off of zero.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when every sorted value lies on or above
    /// the line: the source's loop then reads past the end of its buffer,
    /// which is undefined behaviour.
    pub(super) fn cut_noise(&mut self) -> Result<f64> {
        let n = self.data.len();
        if n == 0 {
            return Ok(0.0);
        }
        let mut sorted = self.data.clone();
        sorted.sort_by(|a, b| b.total_cmp(a));
        let intercept = self.data[0];
        let slope = (self.data[n - 1] - self.data[0]) / n as f64 / CROSSING_SLOPE;
        let cutoff = if slope == 0.0 {
            0.0
        } else {
            let mut index = 1;
            loop {
                if index >= n {
                    return Err(Error::InvalidValue(
                        "pose clustering histogram has no noise level (the source reads past its buffer)"
                            .into(),
                    ));
                }
                if sorted[index] >= intercept + slope * index as f64 {
                    index += 1;
                } else {
                    break;
                }
            }
            sorted[index - 1]
        };
        for value in &mut self.data {
            if *value < cutoff {
                *value = 0.0;
            }
        }
        Ok(cutoff)
    }

    /// The iterative mean and standard-deviation cut-off: [`CUTOFF_LOOPS`]
    /// times, the weighted mean and standard deviation of the bucket indices
    /// within the current range are computed, and the range narrowed to
    /// `[floor(mean - 1.5 sd - low_pad), ceil(mean + 1.5 sd + high_pad))`,
    /// clipped to the histogram. Returns the key and the standard deviation
    /// in key units of the last mean.
    ///
    /// The shift superimposer and the affine superimposer's scaling and low
    /// end use pads of 0 and 1; the affine superimposer's high end uses 1 and
    /// 2.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the statistics are not finite.
    pub(super) fn centroid(&self, low_pad: f64, high_pad: f64) -> Result<(f64, f64)> {
        let mut statistics = BasicStatistics::new();
        let size = self.data.len();
        let (mut begin, mut end) = (0usize, size);
        let mut result = (0.0, 0.0);
        for _ in 0..CUTOFF_LOOPS {
            statistics.update(&self.data[begin..end])?;
            let mean = statistics.mean() + begin as f64;
            let stdev = statistics.variance().sqrt();
            if !mean.is_finite() || !stdev.is_finite() {
                return Err(Error::InvalidValue(
                    "pose clustering histogram statistics are not finite".into(),
                ));
            }
            let low = (mean - CUTOFF_STDEV_MULTIPLIER * stdev - low_pad)
                .max(0.0)
                .floor();
            let high = (mean + CUTOFF_STDEV_MULTIPLIER * stdev + high_pad)
                .min(size as f64)
                .ceil();
            begin = low as usize;
            end = (high as usize).max(begin).min(size);
            result = (self.index_to_key(mean), stdev * self.scale);
        }
        Ok(result)
    }
}

/// Top-hat filter of a `double` range with a structuring element of `length`
/// samples: the input minus its opening (an erosion followed by a dilation),
/// each over windows clipped at the ends.
///
/// The source's `MorphologicalFilter::filterRange` is a template that the
/// superimposers instantiate on `std::vector<double>` iterators, so the
/// histogram is filtered in double precision; the crate's
/// [`MorphologicalFilter`](crate::processing::baseline::MorphologicalFilter)
/// filters `f32` intensities and cannot be used without rounding the votes.
/// Its van Herk erosion and dilation compute the clipped-window extrema, which
/// is what this computes with monotonic deques.
pub(super) fn tophat(values: &[f64], length: usize) -> Vec<f64> {
    let eroded = clipped_extrema(values, length, false);
    let opened = clipped_extrema(&eroded, length, true);
    values.iter().zip(opened).map(|(v, o)| v - o).collect()
}

/// Minima (or maxima) over the windows of `length / 2` samples either side of
/// each sample, clipped at the ends.
fn clipped_extrema(values: &[f64], length: usize, maximum: bool) -> Vec<f64> {
    let radius = length / 2;
    let mut result = Vec::with_capacity(values.len());
    let mut queue: VecDeque<usize> = VecDeque::new();
    let mut next = 0;
    for center in 0..values.len() {
        let start = center.saturating_sub(radius);
        let end = center
            .saturating_add(radius)
            .saturating_add(1)
            .min(values.len());
        while next < end {
            while queue.back().is_some_and(|&i| {
                if maximum {
                    values[i] < values[next]
                } else {
                    values[i] > values[next]
                }
            }) {
                queue.pop_back();
            }
            queue.push_back(next);
            next += 1;
        }
        while queue.front().is_some_and(|&i| i < start) {
            queue.pop_front();
        }
        result.push(values[queue.front().copied().unwrap_or(center)]);
    }
    result
}

/// `std::nth_element(first, first + nth, last, comp)` as the Release
/// toolchain's libstdc++ (GCC 14.4.0, `bits/stl_algo.h` and `bits/stl_heap.h`)
/// implements it, on any `Copy` element: introselect with median-of-three
/// pivots and a heap-select fallback, every comparison and move in the
/// library's order.
///
/// This is the generic form of the `f32`-only port in
/// `processing::noise_estimation::libstdcxx`, which is checked bit for bit
/// against the Release toolchain on 390 vectors; the two should be merged.
pub(super) fn nth_element_by<T: Copy>(v: &mut [T], nth: usize, lt: &impl Fn(&T, &T) -> bool) {
    let n = v.len();
    if n == 0 || nth >= n {
        return;
    }
    let lg = (usize::BITS - 1 - n.leading_zeros()) as usize;
    let (mut first, mut last, mut depth) = (0usize, n, lg * 2);
    while last - first > 3 {
        if depth == 0 {
            heap_select(v, first, nth + 1, last, lt);
            v.swap(first, nth);
            return;
        }
        depth -= 1;
        let mid = first + (last - first) / 2;
        move_median_to_first(v, first, first + 1, mid, last - 1, lt);
        let cut = unguarded_partition(v, first + 1, last, first, lt);
        if cut <= nth {
            first = cut;
        } else {
            last = cut;
        }
    }
    insertion_sort(v, first, last, lt);
}

fn move_median_to_first<T: Copy>(
    v: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    lt: &impl Fn(&T, &T) -> bool,
) {
    let median = if lt(&v[a], &v[b]) {
        if lt(&v[b], &v[c]) {
            b
        } else if lt(&v[a], &v[c]) {
            c
        } else {
            a
        }
    } else if lt(&v[a], &v[c]) {
        a
    } else if lt(&v[b], &v[c]) {
        c
    } else {
        b
    };
    v.swap(result, median);
}

fn unguarded_partition<T: Copy>(
    v: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: usize,
    lt: &impl Fn(&T, &T) -> bool,
) -> usize {
    let end = last;
    loop {
        while first < end && lt(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while last > pivot && lt(&v[pivot], &v[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

fn insertion_sort<T: Copy>(v: &mut [T], first: usize, last: usize, lt: &impl Fn(&T, &T) -> bool) {
    if first == last {
        return;
    }
    for i in first + 1..last {
        if lt(&v[i], &v[first]) {
            let value = v[i];
            v.copy_within(first..i, first + 1);
            v[first] = value;
        } else {
            let value = v[i];
            let mut hole = i;
            let mut next = i - 1;
            while lt(&value, &v[next]) {
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

fn heap_select<T: Copy>(
    v: &mut [T],
    first: usize,
    middle: usize,
    last: usize,
    lt: &impl Fn(&T, &T) -> bool,
) {
    make_heap(v, first, middle, lt);
    for i in middle..last {
        if lt(&v[i], &v[first]) {
            let value = v[i];
            v[i] = v[first];
            adjust_heap(v, first, 0, middle - first, value, lt);
        }
    }
}

fn make_heap<T: Copy>(v: &mut [T], first: usize, last: usize, lt: &impl Fn(&T, &T) -> bool) {
    if last - first < 2 {
        return;
    }
    let len = last - first;
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[first + parent];
        adjust_heap(v, first, parent, len, value, lt);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

fn adjust_heap<T: Copy>(
    v: &mut [T],
    first: usize,
    mut hole: usize,
    len: usize,
    value: T,
    lt: &impl Fn(&T, &T) -> bool,
) {
    let top = hole;
    let mut second = hole;
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if lt(&v[first + second], &v[first + second - 1]) {
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
    // __push_heap
    let mut parent = hole.saturating_sub(1) / 2;
    while hole > top && lt(&v[first + parent], &value) {
        v[first + hole] = v[first + parent];
        hole = parent;
        parent = hole.saturating_sub(1) / 2;
    }
    v[first + hole] = value;
}

/// A parameter entry with the `advanced` tag when asked.
pub(super) fn entry(
    name: &str,
    value: ParamValue,
    description: &str,
    advanced: bool,
) -> ParamEntry {
    let mut e = ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    };
    if advanced {
        e.tags.insert("advanced".into());
    }
    e
}

/// A flat parameter tree.
pub(super) fn tree(entries: Vec<ParamEntry>) -> Result<Param> {
    Param::from_root(ParamNode {
        name: String::new(),
        description: String::new(),
        entries,
        nodes: Vec::new(),
    })
}

/// A float value, converting an integer as `(double)param.getValue` does.
pub(super) fn float(params: &Param, key: &str) -> Result<f64> {
    params
        .value(key)?
        .to_f64()
        .map_err(|_| Error::InvalidValue(format!("parameter '{key}' is not a number")))
}

/// An integer value, as `(Int)param.getValue` does.
pub(super) fn int(params: &Param, key: &str) -> Result<i64> {
    params
        .value(key)?
        .to_i64()
        .map_err(|_| Error::InvalidValue(format!("parameter '{key}' is not an integer")))
}

/// A string value, as `param.getValue(key).toString()` of a string does.
pub(super) fn string(params: &Param, key: &str) -> Result<String> {
    Ok(params
        .value(key)?
        .as_str()
        .map_err(|_| Error::InvalidValue(format!("parameter '{key}' is not a string")))?
        .to_owned())
}

/// The number of elements a `num_used_points` value keeps: `-1` (any negative
/// value, which the source casts to the largest `Size`) keeps all.
pub(super) fn used_points(num_used_points: i64) -> usize {
    usize::try_from(num_used_points).unwrap_or(usize::MAX)
}

/// The intensity-similarity weight of one element pair: the smaller
/// intensity over the larger, after scaling the scene intensity by the ratio
/// of the total intensities.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when it is not finite: two zero
/// intensities make it `0/0`, and the source would carry the NaN into the
/// histogram.
pub(super) fn similarity(model: f32, scene: f32, ratio: f64) -> Result<f64> {
    let int_i = f64::from(model);
    let int_k = f64::from(scene) * ratio;
    let s = if int_i < int_k {
        int_i / int_k
    } else {
        int_k / int_i
    };
    if s.is_finite() {
        Ok(s)
    } else {
        Err(Error::InvalidValue(
            "pose clustering intensity similarity is not finite (zero intensities?)".into(),
        ))
    }
}
