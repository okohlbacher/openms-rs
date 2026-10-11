// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The Wasserstein-network aligner that `WNetMatcher` and
//! `FeatureGroupingAlgorithmWNet` call: `WNetAligner<2>` of the `wnetalign`
//! library on top of `WassersteinNetwork` of the `wnet` library.
//!
//! Third-party dependency of the source, not an OpenMS header. The pinned
//! releases are `wnetalign` v0.9.8 (`cff6a19a6b540247d57044e06b1852afe24346a0`),
//! `wnet` v0.9.11 (`18a15250adb7ed478ef40d26d736bcd873265c74`) and the LEMON
//! copy of `pylmcf` v0.9.8 (`d2c9c52bd67d7198ae17b389d77884357260f114`), which
//! OpenMS fetches in `cmake/cmake_findExternalLibs.cmake`
//! (`WITH_WNETALIGN`). The solver is in
//! [`lemon_network_simplex`](crate::analysis::mapmatching::lemon_network_simplex);
//! `docs/FEATURE_GROUPING_SUPPORT.md` records what was read and why no crate
//! stands in for it.
//!
//! # The model
//!
//! The empirical spectrum `A` is one side and each theoretical spectrum the
//! other. All quantities are scaled to integers first: with
//! `s = sqrt(2^60 / (max(sum A, sum B) * trash_cost))` every position and
//! every intensity is multiplied by `s`, intensities are then truncated to
//! `int64`, and so are the distances computed from the scaled positions.
//! Empirical peak `i` and theoretical peak `j` are joined by an arc of cost
//! `trunc(distance)` when that distance is at most `trunc(max_distance * s)`
//! (and the first coordinates lie within that window, which the source's
//! sorted-window iteration tests in `double`). The arcs split into connected
//! components; a peak without any arc belongs to none and takes no part.
//!
//! Per component, `set_point` builds a flow network of a source, a sink, one
//! node per peak and
//!
//! - source -> empirical peak, capacity its intensity, cost 0;
//! - empirical -> theoretical peak, capacity
//!   `min(theoretical intensity * point, empirical intensity)`, cost the
//!   scaled distance;
//! - theoretical peak -> sink, capacity `trunc(intensity * point)`, cost 0;
//! - source -> sink (the "simple trash"), capacity the total flow, cost
//!   `trunc(trash_cost * s)`;
//!
//! with supply `total = max(sum of empirical capacities, sum of theoretical
//! capacities)` at the source and `-total` at the sink, and solves it by
//! network simplex. The matched pairs of a theoretical spectrum are the
//! arcs carrying flow, and the *consensus* selects among them greedily,
//! heaviest flow first, one pair per peak on either side.
//!
//! # Differences from the source
//!
//! - Only two-dimensional points (the `DIM = 2` instance both OpenMS classes
//!   use).
//! - The optional caller-supplied `scale_factor` constructor argument is not
//!   offered; the source never passes it.
//! - The one-dimensional chain factory (`create_1d`) and the derivative and
//!   residual-distance members are not ported: no OpenMS code reaches them.
//! - The source converts floating-point values to `int64_t` without a range
//!   check, which is undefined behaviour for NaN, infinite and out-of-range
//!   values (`max_distance * s` too large, a non-finite position, a distance
//!   of NaN as `m/z` scaled by an infinite factor produces). Those inputs are
//!   [`Error::InvalidValue`] here.
//! - Resource ceilings the source lacks: [`MAX_PEAKS`] and
//!   [`MAX_MATCHING_ARCS`], checked before the network is allocated.

use crate::error::{Error, Result};
use crate::math::source_sort::source_sort_permutation;

use super::lemon_network_simplex::{Problem, solve};

/// Largest number of peaks (empirical and theoretical together) one aligner
/// accepts. The source's ceiling is the range of `int`; a smaller one keeps
/// the arrays of the simplex addressable on a 32-bit index.
pub const MAX_PEAKS: usize = 1 << 25;

/// Largest number of empirical-to-theoretical arcs one aligner builds. The
/// source has none and allocates them all.
pub const MAX_MATCHING_ARCS: usize = 1 << 24;

fn bad(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

/// The metric positions are compared with (source `DistanceMetric` of `wnet`,
/// `WNetMatcher::DistanceMetric`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DistanceMetric {
    /// Sum of the coordinate differences.
    L1,
    /// Euclidean distance.
    L2,
    /// Largest coordinate difference.
    Linf,
}

/// A set of weighted two-dimensional points (source `Spectrum<2>`).
#[derive(Clone, Debug, PartialEq)]
pub struct Spectrum {
    positions: Vec<[f64; 2]>,
    intensities: Vec<f64>,
}

impl Spectrum {
    /// A spectrum of the given points and weights.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the two vectors differ in length
    /// (the source throws `std::invalid_argument` when the network is built).
    pub fn new(positions: Vec<[f64; 2]>, intensities: Vec<f64>) -> Result<Self> {
        if positions.len() != intensities.len() {
            return Err(bad("Positions and intensities must have the same size"));
        }
        Ok(Self {
            positions,
            intensities,
        })
    }

    /// The number of points.
    pub fn len(&self) -> usize {
        self.intensities.len()
    }

    /// Whether the spectrum has no points.
    pub fn is_empty(&self) -> bool {
        self.intensities.is_empty()
    }

    /// The positions.
    pub fn positions(&self) -> &[[f64; 2]] {
        &self.positions
    }

    /// The intensities.
    pub fn intensities(&self) -> &[f64] {
        &self.intensities
    }

    /// The sum of the intensities, left to right from `0.0`
    /// (`std::accumulate`).
    pub fn sum_intensities(&self) -> f64 {
        let mut sum = 0.0;
        for value in &self.intensities {
            sum += *value;
        }
        sum
    }

    /// The spectrum with every intensity divided by their sum
    /// (`normalized()`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the sum is exactly zero (the
    /// source throws `std::runtime_error`).
    pub fn normalized(&self) -> Result<Self> {
        let total = self.sum_intensities();
        if total == 0.0 {
            return Err(bad(
                "Cannot normalize a spectrum with total intensity of 0.",
            ));
        }
        Ok(Self {
            positions: self.positions.clone(),
            intensities: self.intensities.iter().map(|value| value / total).collect(),
        })
    }
}

/// What a flow-network node stands for.
#[derive(Clone, Copy, Debug)]
enum NodeKind {
    Source,
    Sink,
    Empirical {
        peak: usize,
        intensity: i64,
    },
    Theoretical {
        spectrum: usize,
        peak: usize,
        intensity: i64,
    },
}

/// What a flow-network arc stands for.
#[derive(Clone, Copy, Debug)]
enum ArcKind {
    SrcToEmpirical,
    TheoreticalToSink,
    Matching { cost: i64 },
    Trash { cost: i64 },
}

#[derive(Clone, Copy, Debug)]
struct NetworkArc {
    start: usize,
    end: usize,
    kind: ArcKind,
}

/// One connected component of the matching graph with its flow network
/// (`WassersteinNetworkSubgraph`).
struct Subgraph {
    nodes: Vec<NodeKind>,
    /// Arcs by ascending `(start, end)`, which is the `StaticDigraph` arc id
    /// order of the source's `build()`.
    arcs: Vec<NetworkArc>,
    empirical_intensity: i64,
    flows: Option<Vec<i64>>,
    total_cost: i64,
}

/// The aligner of one empirical spectrum against theoretical spectra
/// (`WNetAligner<2>`).
pub struct WNetAligner {
    scale_factor: f64,
    no_theoretical: usize,
    subgraphs: Vec<Subgraph>,
}

/// `static_cast<int64_t>(value)`, refused where the conversion is undefined.
fn to_i64(value: f64, what: &str) -> Result<i64> {
    // 2^63 is exactly representable; the conversion is defined on
    // (-2^63 - 1, 2^63).
    if !value.is_finite()
        || !(-9_223_372_036_854_775_808.0..9_223_372_036_854_775_808.0).contains(&value)
    {
        return Err(bad(format!(
            "{what} ({value}) does not fit a 64-bit integer after scaling"
        )));
    }
    Ok(value as i64)
}

/// `compute_scale_factor`.
fn compute_scale_factor(
    empirical: &Spectrum,
    theoretical: &[&Spectrum],
    trash_cost: f64,
) -> Result<f64> {
    if trash_cost < 0.0 {
        return Err(bad(format!(
            "trash_cost must be non-negative, got {trash_cost}"
        )));
    }
    let almost_maxint = (1u64 << 60) as f64;
    let empirical_sum = empirical.sum_intensities();
    let mut theoretical_sum = 0.0;
    for spectrum in theoretical {
        theoretical_sum += spectrum.sum_intensities();
    }
    // std::max(a, b) is (a < b) ? b : a.
    let max_sum = if empirical_sum < theoretical_sum {
        theoretical_sum
    } else {
        empirical_sum
    };
    if max_sum <= 0.0 {
        return Err(bad(format!(
            "max intensity sum must be positive (got {max_sum}). Are all spectra empty?"
        )));
    }
    let product = max_sum * trash_cost;
    if product.is_infinite() {
        return Err(bad(format!(
            "max_sum * trash_cost overflows double ({max_sum} * {trash_cost})"
        )));
    }
    let factor = (almost_maxint / product).sqrt();
    if !factor.is_finite() || factor <= 0.0 {
        return Err(bad(format!(
            "the scale factor of intensity sum {max_sum} and trash cost {trash_cost} is not finite and positive"
        )));
    }
    Ok(factor)
}

/// A distance of two points under `metric`, in `double` as the source
/// computes it.
fn distance(metric: DistanceMetric, a: &[f64; 2], b: &[f64; 2]) -> f64 {
    let d0 = a[0] - b[0];
    let d1 = a[1] - b[1];
    match metric {
        DistanceMetric::L1 => d0.abs() + d1.abs(),
        DistanceMetric::L2 => ((d0 * d0) + (d1 * d1)).sqrt(),
        DistanceMetric::Linf => {
            // std::max({x, y}) keeps x unless x < y.
            let x = d0.abs();
            let y = d1.abs();
            if x < y { y } else { x }
        }
    }
}

/// A union-find over the peaks, for the connected components.
struct Components {
    parent: Vec<u32>,
}

impl Components {
    fn new(len: usize) -> Self {
        Self {
            parent: (0..len as u32).collect(),
        }
    }
    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let grand = self.parent[self.parent[x as usize] as usize];
            self.parent[x as usize] = grand;
            x = grand;
        }
        x
    }
    fn union(&mut self, a: u32, b: u32) {
        let a = self.find(a);
        let b = self.find(b);
        if a != b {
            let (low, high) = if a < b { (a, b) } else { (b, a) };
            self.parent[high as usize] = low;
        }
    }
}

impl WNetAligner {
    /// Build the network of `empirical` against `theoretical` (source
    /// constructor): scale, find the pairs closer than `max_distance`, split
    /// into components and add the trash arc. Nothing is solved until
    /// [`Self::set_point`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a negative `trash_cost`, a
    /// non-positive intensity sum, an overflowing scale factor, no or an empty
    /// theoretical spectrum, an empty empirical spectrum, mismatched vector
    /// lengths, a value that does not fit `int64` after scaling, or a network
    /// beyond [`MAX_PEAKS`] or [`MAX_MATCHING_ARCS`].
    pub fn new(
        empirical: &Spectrum,
        theoretical: &[&Spectrum],
        metric: DistanceMetric,
        max_distance: f64,
        trash_cost: f64,
    ) -> Result<Self> {
        let scale_factor = compute_scale_factor(empirical, theoretical, trash_cost)?;
        if theoretical.is_empty() {
            return Err(bad("Need at least one theoretical spectrum"));
        }
        if empirical.is_empty() {
            return Err(bad("Empirical spectrum is empty"));
        }
        for (index, spectrum) in theoretical.iter().enumerate() {
            if spectrum.is_empty() {
                return Err(bad(format!(
                    "Theoretical distribution at index {index} is empty."
                )));
            }
        }
        let mut total_peaks = empirical.len();
        for spectrum in theoretical {
            total_peaks = total_peaks.saturating_add(spectrum.len());
        }
        if total_peaks > MAX_PEAKS {
            return Err(bad(format!("more than {MAX_PEAKS} peaks")));
        }

        // Scaled integer intensities and scaled positions.
        let scaled_intensities = |spectrum: &Spectrum| -> Result<Vec<i64>> {
            spectrum
                .intensities
                .iter()
                .map(|value| to_i64(value * scale_factor, "an intensity"))
                .collect()
        };
        let scaled_positions = |spectrum: &Spectrum| -> Result<Vec<[f64; 2]>> {
            let mut out = Vec::with_capacity(spectrum.len());
            for position in &spectrum.positions {
                let scaled = [position[0] * scale_factor, position[1] * scale_factor];
                if !scaled[0].is_finite() || !scaled[1].is_finite() {
                    return Err(bad("a position is not finite after scaling"));
                }
                out.push(scaled);
            }
            Ok(out)
        };
        let empirical_intensities = scaled_intensities(empirical)?;
        let empirical_positions = scaled_positions(empirical)?;
        let mut theoretical_intensities = Vec::with_capacity(theoretical.len());
        let mut theoretical_positions = Vec::with_capacity(theoretical.len());
        for spectrum in theoretical {
            theoretical_intensities.push(scaled_intensities(spectrum)?);
            theoretical_positions.push(scaled_positions(spectrum)?);
        }

        // network.create(..., static_cast<int64_t>(max_distance * scale)).
        let max_dist = to_i64(max_distance * scale_factor, "the maximum distance")?;
        let window = max_dist as f64;

        // The pairs within the window and at most max_dist apart.
        let empirical_count = empirical.len();
        let mut theoretical_offsets = Vec::with_capacity(theoretical.len());
        let mut offset = 0usize;
        for spectrum in theoretical {
            theoretical_offsets.push(offset);
            offset += spectrum.len();
        }
        let theoretical_count = offset;
        // (empirical peak, global theoretical peak, scaled distance)
        let mut pairs: Vec<(u32, u32, i64)> = Vec::new();
        for (spectrum_index, positions) in theoretical_positions.iter().enumerate() {
            let mut sorted: Vec<u32> = (0..positions.len() as u32).collect();
            sorted.sort_by(|a, b| {
                positions[*a as usize][0]
                    .total_cmp(&positions[*b as usize][0])
                    .then(a.cmp(b))
            });
            let sorted_first: Vec<f64> = sorted
                .iter()
                .map(|index| positions[*index as usize][0])
                .collect();
            for (e, e_position) in empirical_positions.iter().enumerate() {
                let low = e_position[0] - window;
                let high = e_position[0] + window;
                let start = sorted_first.partition_point(|value| *value < low);
                for k in start..sorted.len() {
                    if sorted_first[k] > high {
                        break;
                    }
                    let t = sorted[k] as usize;
                    let d = distance(metric, e_position, &positions[t]);
                    let scaled = to_i64(d, "a distance")?;
                    if scaled <= max_dist {
                        if pairs.len() >= MAX_MATCHING_ARCS {
                            return Err(bad(format!(
                                "more than {MAX_MATCHING_ARCS} candidate pairs"
                            )));
                        }
                        pairs.push((
                            e as u32,
                            (theoretical_offsets[spectrum_index] + t) as u32,
                            scaled,
                        ));
                    }
                }
            }
        }
        // The arcs of a subgraph are sorted by (start, end) in build(); one
        // global sort by (empirical, theoretical) gives every subgraph's order.
        pairs.sort_unstable_by_key(|pair| (pair.0, pair.1));

        // Connected components over the peaks that have a pair; the source
        // numbers a component by its lowest node id.
        let mut components = Components::new(empirical_count + theoretical_count);
        let mut has_pair = vec![false; empirical_count + theoretical_count];
        for &(e, t, _) in &pairs {
            let t_node = empirical_count as u32 + t;
            components.union(e, t_node);
            has_pair[e as usize] = true;
            has_pair[t_node as usize] = true;
        }
        let mut component_of_root: std::collections::HashMap<u32, usize> =
            std::collections::HashMap::new();
        let mut component_nodes: Vec<Vec<u32>> = Vec::new();
        for node in 0..(empirical_count + theoretical_count) as u32 {
            if !has_pair[node as usize] {
                continue;
            }
            let root = components.find(node);
            let index = *component_of_root.entry(root).or_insert_with(|| {
                component_nodes.push(Vec::new());
                component_nodes.len() - 1
            });
            component_nodes[index].push(node);
        }
        let mut component_pairs: Vec<Vec<(u32, u32, i64)>> =
            vec![Vec::new(); component_nodes.len()];
        for pair in &pairs {
            let root = components.find(pair.0);
            component_pairs[component_of_root[&root]].push(*pair);
        }
        drop(pairs);

        // The theoretical spectrum a global theoretical peak belongs to.
        let spectrum_of = |global: usize| -> (usize, usize) {
            let index = theoretical_offsets
                .partition_point(|start| *start <= global)
                .saturating_sub(1);
            (index, global - theoretical_offsets[index])
        };

        let trash_scaled = to_i64(trash_cost * scale_factor, "the trash cost")?;
        let mut subgraphs = Vec::with_capacity(component_nodes.len());
        let mut local_of = vec![0u32; empirical_count + theoretical_count];
        for (nodes, component_arcs) in component_nodes.iter().zip(&component_pairs) {
            let mut kinds: Vec<NodeKind> = Vec::with_capacity(nodes.len() + 2);
            kinds.push(NodeKind::Source);
            kinds.push(NodeKind::Sink);
            let mut empirical_locals: Vec<usize> = Vec::new();
            let mut theoretical_locals: Vec<usize> = Vec::new();
            for &node in nodes {
                let local = kinds.len();
                local_of[node as usize] = local as u32;
                if (node as usize) < empirical_count {
                    kinds.push(NodeKind::Empirical {
                        peak: node as usize,
                        intensity: empirical_intensities[node as usize],
                    });
                    empirical_locals.push(local);
                } else {
                    let global = node as usize - empirical_count;
                    let (spectrum, peak) = spectrum_of(global);
                    kinds.push(NodeKind::Theoretical {
                        spectrum,
                        peak,
                        intensity: theoretical_intensities[spectrum][peak],
                    });
                    theoretical_locals.push(local);
                }
            }

            // Arcs in (start, end) order: the trash arc (0, 1), the source's
            // arcs to the empirical peaks, the matching arcs by empirical and
            // theoretical peak, and the theoretical peaks' arcs to the sink.
            let mut arcs: Vec<NetworkArc> =
                Vec::with_capacity(1 + nodes.len() + component_arcs.len());
            arcs.push(NetworkArc {
                start: 0,
                end: 1,
                kind: ArcKind::Trash { cost: trash_scaled },
            });
            let mut empirical_intensity: i64 = 0;
            for &local in &empirical_locals {
                if let NodeKind::Empirical { intensity, .. } = kinds[local] {
                    empirical_intensity = empirical_intensity.wrapping_add(intensity);
                }
                arcs.push(NetworkArc {
                    start: 0,
                    end: local,
                    kind: ArcKind::SrcToEmpirical,
                });
            }
            for &(e, t, cost) in component_arcs {
                arcs.push(NetworkArc {
                    start: local_of[e as usize] as usize,
                    end: local_of[empirical_count + t as usize] as usize,
                    kind: ArcKind::Matching { cost },
                });
            }
            for &local in &theoretical_locals {
                arcs.push(NetworkArc {
                    start: local,
                    end: 1,
                    kind: ArcKind::TheoreticalToSink,
                });
            }
            subgraphs.push(Subgraph {
                nodes: kinds,
                arcs,
                empirical_intensity,
                flows: None,
                total_cost: 0,
            });
        }

        Ok(Self {
            scale_factor,
            no_theoretical: theoretical.len(),
            subgraphs,
        })
    }

    /// The factor positions, intensities and costs were multiplied by.
    pub fn scale_factor(&self) -> f64 {
        self.scale_factor
    }

    /// The number of theoretical spectra.
    pub fn no_theoretical_spectra(&self) -> usize {
        self.no_theoretical
    }

    /// The number of connected components (source `no_subgraphs()`).
    pub fn no_subgraphs(&self) -> usize {
        self.subgraphs.len()
    }

    /// Fix the proportion of every theoretical spectrum and solve every
    /// component (source `set_point`; the network is rebuilt per call, which
    /// is what the source's `set_point` does with its solver).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `point` does not have one entry
    /// per theoretical spectrum, an entry is not finite or makes a capacity
    /// not fit `int64`, or the solver reports the problem infeasible or
    /// unbounded (not reachable from this network).
    pub fn set_point(&mut self, point: &[f64]) -> Result<()> {
        if point.len() != self.no_theoretical {
            return Err(bad(format!(
                "Point dimension: {} does not match number of target distributions: {}",
                point.len(),
                self.no_theoretical
            )));
        }
        for subgraph in &mut self.subgraphs {
            let node_num = subgraph.nodes.len();
            let arc_num = subgraph.arcs.len();
            let mut upper = vec![0i64; arc_num];
            let mut cost = vec![0i64; arc_num];
            let mut theoretical_intensity: i64 = 0;
            let mut trash_arc = 0usize;
            for (index, arc) in subgraph.arcs.iter().enumerate() {
                match arc.kind {
                    ArcKind::Matching { cost: arc_cost } => {
                        let NodeKind::Theoretical {
                            spectrum,
                            intensity: theoretical,
                            ..
                        } = subgraph.nodes[arc.end]
                        else {
                            return Err(bad("a matching arc does not end at a theoretical peak"));
                        };
                        let NodeKind::Empirical {
                            intensity: empirical,
                            ..
                        } = subgraph.nodes[arc.start]
                        else {
                            return Err(bad("a matching arc does not start at an empirical peak"));
                        };
                        // std::min<double>(a, b) is (b < a) ? b : a.
                        let a = theoretical as f64 * point[spectrum];
                        let b = empirical as f64;
                        let capped = if b < a { b } else { a };
                        upper[index] = to_i64(capped, "a capacity")?;
                        cost[index] = arc_cost;
                    }
                    ArcKind::SrcToEmpirical => {
                        let NodeKind::Empirical { intensity, .. } = subgraph.nodes[arc.end] else {
                            return Err(bad("a source arc does not end at an empirical peak"));
                        };
                        upper[index] = intensity;
                    }
                    ArcKind::TheoreticalToSink => {
                        let NodeKind::Theoretical {
                            spectrum,
                            intensity,
                            ..
                        } = subgraph.nodes[arc.start]
                        else {
                            return Err(bad("a sink arc does not start at a theoretical peak"));
                        };
                        let capacity = to_i64(intensity as f64 * point[spectrum], "a capacity")?;
                        upper[index] = capacity;
                        theoretical_intensity = theoretical_intensity.wrapping_add(capacity);
                    }
                    ArcKind::Trash { cost: arc_cost } => {
                        trash_arc = index;
                        cost[index] = arc_cost;
                    }
                }
            }
            // std::max<int64_t>(empirical, theoretical).
            let total_flow = if subgraph.empirical_intensity < theoretical_intensity {
                theoretical_intensity
            } else {
                subgraph.empirical_intensity
            };
            upper[trash_arc] = total_flow;
            let mut supply = vec![0i64; node_num];
            supply[0] = total_flow;
            supply[1] = total_flow.wrapping_neg();
            let sources: Vec<usize> = subgraph.arcs.iter().map(|arc| arc.start).collect();
            let targets: Vec<usize> = subgraph.arcs.iter().map(|arc| arc.end).collect();
            let flows = solve(&Problem {
                node_num,
                arc_source: &sources,
                arc_target: &targets,
                upper: &upper,
                cost: &cost,
                supply: &supply,
            })?;
            let mut total: i64 = 0;
            for (flow, arc_cost) in flows.iter().zip(&cost) {
                total = total.wrapping_add(flow.wrapping_mul(*arc_cost));
            }
            subgraph.total_cost = total;
            subgraph.flows = Some(flows);
        }
        Ok(())
    }

    /// The optimal total transport cost of the last [`Self::set_point`],
    /// divided by the squared scale factor (source `total_cost()`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] before `set_point`.
    pub fn total_cost(&self) -> Result<f64> {
        let mut total: i64 = 0;
        for subgraph in &self.subgraphs {
            if subgraph.flows.is_none() {
                return Err(bad(
                    "You must call build() and set_point() before calling total_cost().",
                ));
            }
            total = total.wrapping_add(subgraph.total_cost);
        }
        Ok(total as f64 / (self.scale_factor * self.scale_factor))
    }

    /// The matched pairs of one theoretical spectrum: the empirical and
    /// theoretical peak indices and the flow of every matching arc that
    /// carries flow, component by component and arc by arc (source
    /// `flows_for_target`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] before `set_point`.
    pub fn flows_for_target(&self, target: usize) -> Result<(Vec<usize>, Vec<usize>, Vec<i64>)> {
        let mut empirical = Vec::new();
        let mut theoretical = Vec::new();
        let mut flows = Vec::new();
        for subgraph in &self.subgraphs {
            let Some(solved) = &subgraph.flows else {
                return Err(bad("flows requested before set_point"));
            };
            for (arc, flow) in subgraph.arcs.iter().zip(solved) {
                if *flow == 0 {
                    continue;
                }
                if let ArcKind::Matching { .. } = arc.kind {
                    if let (
                        NodeKind::Empirical { peak: e, .. },
                        NodeKind::Theoretical {
                            spectrum, peak: t, ..
                        },
                    ) = (subgraph.nodes[arc.start], subgraph.nodes[arc.end])
                    {
                        if spectrum == target {
                            empirical.push(e);
                            theoretical.push(t);
                            flows.push(*flow);
                        }
                    }
                }
            }
        }
        Ok((empirical, theoretical, flows))
    }

    /// The greedy one-to-one consensus of one theoretical spectrum: its
    /// matching arcs sorted by descending flow (the C++ Release build's
    /// `std::sort`, so equal flows end in its order), then taken while
    /// neither peak is used (source `consensus_for_target`).
    ///
    /// # Errors
    ///
    /// As [`Self::flows_for_target`].
    pub fn consensus_for_target(&self, target: usize) -> Result<(Vec<usize>, Vec<usize>)> {
        let (empirical, theoretical, flows) = self.flows_for_target(target)?;
        let order = source_sort_permutation(flows.len(), |a, b| flows[a] > flows[b])?;
        let mut used_empirical = std::collections::HashSet::new();
        let mut used_theoretical = std::collections::HashSet::new();
        let mut consensus_empirical = Vec::new();
        let mut consensus_theoretical = Vec::new();
        for index in order {
            let e = empirical[index];
            let t = theoretical[index];
            if !used_empirical.contains(&e) && !used_theoretical.contains(&t) {
                used_empirical.insert(e);
                used_theoretical.insert(t);
                consensus_empirical.push(e);
                consensus_theoretical.push(t);
            }
        }
        Ok((consensus_empirical, consensus_theoretical))
    }
}
