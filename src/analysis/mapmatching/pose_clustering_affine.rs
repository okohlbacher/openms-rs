// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Affine superimposer by pose clustering.
//!
//! Ports `ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h` and its
//! `.cpp`, including the file-local helpers `initializeHashTables`,
//! `affineTransformationHashing`, `scalingEstimate`, `shiftEstimate` and
//! `computeIntensityRatio`. See `docs/MAP_ALIGNMENT_SUPPORT.md`.

use super::superimposer::{
    BaseSuperimposer, Dump, Histogram, Superimposition, WINDOW_LENGTH_BASELINE, entry, float, g,
    int, nth_element_by, similarity, string, tree, used_points,
};
use std::sync::atomic::{AtomicI32, Ordering};

/// The source's function-local `static Int dump_buckets_serial` of `run`,
/// incremented on every run past the empty-map check and naming the debug
/// dumps.
static AFFINE_DUMP_SERIAL: AtomicI32 = AtomicI32::new(0);
use crate::analysis::transformations::TransformationDescription;
use crate::concept::fuzzy_string_comparator::format_g;
use crate::concept::progress_logger::ProgressReporter;
use crate::kernel::{ConsensusMap, Peak2D};
use crate::math::source_sort::source_sort_by;
use crate::param::{DefaultParamHandler, Param, ParamValue};
use crate::{Error, Result};

/// Parameter handler name of the source class.
pub const POSE_CLUSTERING_AFFINE_SUPERIMPOSER_NAME: &str = "PoseClusteringAffineSuperimposer";

/// A superimposer that finds an affine retention-time transformation by pose
/// clustering (source `PoseClusteringAffineSuperimposer`).
///
/// It works on two maps and computes the affine transformation that maps the
/// elements of the scene as near as possible onto those of the model. Every
/// pair of model elements and pair of scene elements of similar m/z defines a
/// potential "pose" of the scene relative to the model; the scaling, and the
/// images of the start and end of the retention-time range, are hashed into
/// histograms weighted by intensity similarity, and the clusters of those
/// histograms give the transformation, hence the name pose clustering. See
/// also [`PoseClusteringShiftSuperimposer`](super::pose_clustering_shift::PoseClusteringShiftSuperimposer).
///
/// The estimation runs in two rounds: the first hashes only the scaling and
/// estimates a likely range for it; the second hashes the images of the
/// range ends for the quadruples whose scaling falls into that range.
#[derive(Clone, Debug, PartialEq)]
pub struct PoseClusteringAffineSuperimposer {
    /// Maximum m/z deviation of corresponding elements in different maps
    /// (source `mz_pair_max_distance`, default 0.5, minimum 0). This applies
    /// to the pairs considered in hashing.
    pub mz_pair_max_distance: f64,
    /// Minimum separation of the pairs considered within each map, as a
    /// fraction of the total elution time interval (source
    /// `rt_pair_distance_fraction`, default 0.1, from 0 to 1, advanced).
    pub rt_pair_distance_fraction: f64,
    /// Maximum number of elements considered in each map, selected by
    /// intensity (source `num_used_points`, default 2000, minimum -1). Reduces
    /// the running time and disregards weak signals; -1 uses all points.
    pub num_used_points: i64,
    /// Size of the histogram buckets of the (logarithmic) scaling (source
    /// `scaling_bucket_size`, default 0.005, minimum 0). A good choice is a
    /// bit smaller than the error expected from repeated runs.
    pub scaling_bucket_size: f64,
    /// Size of the histogram buckets of the shifts at the ends of the
    /// retention-time interval (source `shift_bucket_size`, default 3.0,
    /// minimum 0). A good choice is about the time between consecutive MS
    /// scans.
    pub shift_bucket_size: f64,
    /// Maximal shift considered during histogramming, in seconds, in both
    /// directions (source `max_shift`, default 1000, minimum 0, advanced).
    pub max_shift: f64,
    /// Maximal scaling considered during histogramming; the minimal scaling
    /// is its reciprocal (source `max_scaling`, default 2.0, minimum 1,
    /// advanced).
    pub max_scaling: f64,
    /// `[DEBUG]` base file name for dumps of the histograms (source
    /// `dump_buckets`, advanced). When non-empty, each run writes
    /// `<dump_buckets>_scale_<serial>`, `_low_<serial>` and `_high_<serial>`:
    /// each histogram before and after the top-hat filter and after the
    /// noise cut-off, and each cut-off round.
    pub dump_buckets: String,
    /// `[DEBUG]` base file name for dumps of the hashed quadruples (source
    /// `dump_pairs`, advanced; large). When non-empty, each run writes
    /// `<dump_pairs>_phase_two_<serial>`, one line per quadruple of the second
    /// round; the first round opens the same file and so leaves only its
    /// header, as in the source.
    ///
    /// The serial numbers the runs of this class in this process, starting at
    /// 1: the source's function-local `static` counter, which is a
    /// process-wide counter here too.
    pub dump_pairs: String,
    /// Maximum number of buckets of each histogram (native limit).
    pub max_buckets: usize,
    /// Maximum number of element quadruples examined in one hashing round
    /// (native limit; the work grows with the fourth power of the number of
    /// elements within the m/z windows). Counted as the quadruples are
    /// visited, so an input over the limit fails early.
    pub max_quadruples: u64,
}

impl Default for PoseClusteringAffineSuperimposer {
    fn default() -> Self {
        Self {
            mz_pair_max_distance: 0.5,
            rt_pair_distance_fraction: 0.1,
            num_used_points: 2000,
            scaling_bucket_size: 0.005,
            shift_bucket_size: 3.0,
            max_shift: 1000.0,
            max_scaling: 2.0,
            dump_buckets: String::new(),
            dump_pairs: String::new(),
            max_buckets: 10_000_000,
            max_quadruples: 20_000_000_000,
        }
    }
}

impl PoseClusteringAffineSuperimposer {
    /// The source parameter defaults, as the constructor sets them:
    /// `mz_pair_max_distance`, `rt_pair_distance_fraction`,
    /// `num_used_points`, `scaling_bucket_size`, `shift_bucket_size`,
    /// `max_shift`, `max_scaling`, `dump_buckets` and `dump_pairs`, with their
    /// descriptions, restrictions and `advanced` tags, in declaration order.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this superimposer in the layout of
    /// [`PoseClusteringAffineSuperimposer::defaults`] (source
    /// `getParameters`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn to_param(&self) -> Result<Param> {
        let mut mz = entry(
            "mz_pair_max_distance",
            ParamValue::Float(self.mz_pair_max_distance),
            "Maximum of m/z deviation of corresponding elements in different maps.  This condition applies to the pairs considered in hashing.",
            false,
        );
        mz.min_float = 0.0;
        let mut fraction = entry(
            "rt_pair_distance_fraction",
            ParamValue::Float(self.rt_pair_distance_fraction),
            "Within each of the two maps, the pairs considered for pose clustering must be separated by at least this fraction of the total elution time interval (i.e., max - min).  ",
            true,
        );
        fraction.min_float = 0.0;
        fraction.max_float = 1.0;
        let mut points = entry(
            "num_used_points",
            ParamValue::Integer(self.num_used_points),
            "Maximum number of elements considered in each map (selected by intensity).  Use this to reduce the running time and to disregard weak signals during alignment.  For using all points, set this to -1.",
            false,
        );
        points.min_int = -1;
        let mut scaling = entry(
            "scaling_bucket_size",
            ParamValue::Float(self.scaling_bucket_size),
            "The scaling of the retention time interval is being hashed into buckets of this size during pose clustering.  A good choice for this would be a bit smaller than the error you would expect from repeated runs.",
            false,
        );
        scaling.min_float = 0.0;
        let mut shift_bucket = entry(
            "shift_bucket_size",
            ParamValue::Float(self.shift_bucket_size),
            "The shift at the lower (respectively, higher) end of the retention time interval is being hashed into buckets of this size during pose clustering.  A good choice for this would be about the time between consecutive MS scans.",
            false,
        );
        shift_bucket.min_float = 0.0;
        let mut shift = entry(
            "max_shift",
            ParamValue::Float(self.max_shift),
            "Maximal shift which is considered during histogramming (in seconds).  This applies for both directions.",
            true,
        );
        shift.min_float = 0.0;
        let mut max_scaling = entry(
            "max_scaling",
            ParamValue::Float(self.max_scaling),
            "Maximal scaling which is considered during histogramming.  The minimal scaling is the reciprocal of this.",
            true,
        );
        max_scaling.min_float = 1.0;
        tree(vec![
            mz,
            fraction,
            points,
            scaling,
            shift_bucket,
            shift,
            max_scaling,
            entry(
                "dump_buckets",
                ParamValue::String(self.dump_buckets.clone()),
                "[DEBUG] If non-empty, base filename where hash table buckets will be dumped to.  A serial number for each invocation will be appended automatically.",
                true,
            ),
            entry(
                "dump_pairs",
                ParamValue::String(self.dump_pairs.clone()),
                "[DEBUG] If non-empty, base filename where the individual hashed pairs will be dumped to (large!).  A serial number for each invocation will be appended automatically.",
                true,
            ),
        ])
    }

    /// Build a superimposer from parameters, as source `setParameters`.
    /// Missing parameters take their defaults; unknown ones are ignored.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its restriction (the source's `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`PoseClusteringAffineSuperimposer::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`PoseClusteringAffineSuperimposer::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(POSE_CLUSTERING_AFFINE_SUPERIMPOSER_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |p| {
            Ok(Self {
                mz_pair_max_distance: float(p, "mz_pair_max_distance")?,
                rt_pair_distance_fraction: float(p, "rt_pair_distance_fraction")?,
                num_used_points: int(p, "num_used_points")?,
                scaling_bucket_size: float(p, "scaling_bucket_size")?,
                shift_bucket_size: float(p, "shift_bucket_size")?,
                max_shift: float(p, "max_shift")?,
                max_scaling: float(p, "max_scaling")?,
                dump_buckets: string(p, "dump_buckets")?,
                dump_pairs: string(p, "dump_pairs")?,
                ..Self::default()
            })
        })
    }

    fn validate(&self) -> Result<()> {
        let positive = |value: f64| value.is_finite() && value > 0.0;
        if self.mz_pair_max_distance < 0.0
            || !self.mz_pair_max_distance.is_finite()
            || !self.rt_pair_distance_fraction.is_finite()
            || !positive(self.scaling_bucket_size)
            || !positive(self.shift_bucket_size)
            || !self.max_shift.is_finite()
            || !positive(self.max_scaling)
        {
            // Zero bucket sizes make the source convert an infinite bucket
            // count to Int, which is undefined behaviour.
            return Err(Error::InvalidValue(
                "invalid pose clustering distances, bucket sizes or maxima".into(),
            ));
        }
        Ok(())
    }

    /// Estimate the affine transformation from `scene` onto `model` given as
    /// plain points (source `run(const std::vector<Peak2D>&, const
    /// std::vector<Peak2D>&, TransformationDescription&)`).
    ///
    /// The steps, as in the source:
    ///
    /// 1. Keep the `num_used_points` most intense points of each map
    ///    (`std::nth_element` on the reversed vector, as the Release build's
    ///    libstdc++ does it), then sort each by m/z (`std::sort`).
    /// 2. Take `rt_low` and `rt_high` as the means of the two maps' smallest
    ///    and largest retention times, over *all* points; pairs within a map
    ///    must be `rt_pair_distance_fraction * (rt_high - rt_low)` apart.
    ///    When the maps' crude scaling or shift exceeds `max_scaling` or
    ///    `max_shift` (with a factor 1.2 of slack), a warning is returned.
    /// 3. Size the histograms: `ceil(ln(max_scaling) / scaling_bucket_size) + 1`
    ///    buckets either side of scaling 1 (in log space), and
    ///    `4 + 2 * ceil(max_shift * max_scaling / shift_bucket_size)` either
    ///    side of `rt_low` and `rt_high`.
    /// 4. Hash every quadruple `(i, j)` of model and `(k, l)` of scene points
    ///    whose partners lie within `mz_pair_max_distance`, whose RT
    ///    differences are at least the minimum distance and have the same
    ///    sign: the scaling is `(rt_j - rt_i) / (rt_l - rt_k)`, weighted by
    ///    both pairs' intensity similarity (scene intensities scaled by the
    ///    ratio of the total intensities) and by the inverse number of
    ///    elements in each m/z window, less 0.1; windows of ten or more
    ///    elements are skipped. The first round hashes `ln(scaling)`; a top-hat
    ///    filter, a noise cut-off and three rounds of a mean and
    ///    standard-deviation cut-off give the likely scaling range
    ///    `exp(mean -/+ sd)`. The second round hashes, for quadruples in that
    ///    range, the images of `rt_low` and `rt_high`, whose centroids are
    ///    estimated the same way (with a wider range at the high end).
    /// 5. The transformation is the line through `(rt_low, image_low)` and
    ///    `(rt_high, image_high)`, fitted as `linear` from explicit slope and
    ///    intercept without data.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an empty map (the source's
    /// `Exception::IllegalArgument`), a non-finite slope or intercept (the
    /// source's `Exception::InvalidValue`, "Superimposer could not compute an
    /// initial transformation!"), an invalid parameter value, a histogram of
    /// more than [`max_buckets`](Self::max_buckets) or more quadruples than
    /// [`max_quadruples`](Self::max_quadruples), a non-finite total intensity
    /// ratio or pair similarity, or a histogram whose noise level the source
    /// cannot determine without reading past its buffer; [`Error::Io`] when a
    /// debug dump cannot be written (the source's stream fails silently).
    /// Progress errors are returned as well.
    pub fn superimpose_points(
        &self,
        model: &[Peak2D],
        scene: &[Peak2D],
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Superimposition> {
        if model.is_empty() || scene.is_empty() {
            return Err(Error::InvalidValue(
                "One of the input maps is empty! This is not allowed!".into(),
            ));
        }
        self.validate()?;
        for point in model.iter().chain(scene) {
            if !point.position[0].is_finite()
                || !point.position[1].is_finite()
                || !point.intensity.is_finite()
            {
                return Err(Error::InvalidValue(
                    "pose clustering input coordinates must be finite".into(),
                ));
            }
        }
        progress.section(0, 100, "affine pose clustering", |progress| {
            self.estimate(model, scene, progress)
        })
    }

    fn estimate(
        &self,
        map_model: &[Peak2D],
        map_scene: &[Peak2D],
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Superimposition> {
        let mut actual = 0i64;
        for _ in 0..3 {
            actual += 1;
            progress.set(actual)?;
        }
        // Step 1: the most abundant points, sorted by m/z.
        let limit = used_points(self.num_used_points);
        let intensity_less = |a: &Peak2D, b: &Peak2D| a.intensity < b.intensity;
        let truncate = |map: &[Peak2D]| -> Vec<Peak2D> {
            let mut points = map.to_vec();
            if points.len() > limit {
                // std::nth_element(rbegin, rbegin + (size - limit), rend, IntensityLess)
                let nth = points.len() - limit;
                points.reverse();
                nth_element_by(&mut points, nth, &intensity_less);
                points.reverse();
                points.truncate(limit);
            }
            points
        };
        let mut model = truncate(map_model);
        actual += 1;
        progress.set(actual)?;
        let mut scene = truncate(map_scene);
        actual += 1;
        progress.set(actual)?;
        source_sort_by(&mut model, |a, b| a.position[1] < b.position[1])?;
        source_sort_by(&mut scene, |a, b| a.position[1] < b.position[1])?;
        actual = 10;
        progress.set(actual)?;

        // Preprocessing: the RT ends over all points (Peak2D::RTLess).
        let min_rt = |map: &[Peak2D]| {
            map.iter()
                .map(|p| p.position[0])
                .fold(f64::INFINITY, f64::min)
        };
        let max_rt = |map: &[Peak2D]| {
            map.iter()
                .map(|p| p.position[0])
                .fold(f64::NEG_INFINITY, f64::max)
        };
        let (model_min, scene_min) = (min_rt(map_model), min_rt(map_scene));
        let (model_max, scene_max) = (max_rt(map_model), max_rt(map_scene));
        let rt_low = (model_min + scene_min) / 2.0;
        let rt_high = (model_max + scene_max) / 2.0;

        // Sanity check of a crude estimate of the shift and slope.
        let mut warnings = Vec::new();
        {
            let shift = (model_min - scene_min).abs();
            let slope = (model_max - model_min) / (scene_max - scene_min);
            if self.max_scaling < slope * 1.2 || 1.0 / self.max_scaling > slope / 1.2 {
                warnings.push(format!(
                    "WARNING: your map likely has a scaling around {} but your parameters only allow for a maximal scaling of {}",
                    format_g(slope),
                    format_g(self.max_scaling)
                ));
                warnings.push("It is strongly advised to adjust your max_scaling factor".into());
            }
            if self.max_shift < shift * 1.2 {
                warnings.push(format!(
                    "WARNING: your map likely has a shift around {} but your parameters only allow for a maximal shift of {}",
                    format_g(shift),
                    format_g(self.max_shift)
                ));
                warnings.push("It is strongly advised to adjust your max_shift factor".into());
            }
        }
        let rt_pair_min_distance = self.rt_pair_distance_fraction * (rt_high - rt_low);

        // Step 2: the histograms (initializeHashTables).
        let scaling_half = (self.max_scaling.ln() / self.scaling_bucket_size).ceil() + 1.0;
        let rt_half =
            4.0 + 2.0 * (self.max_shift * self.max_scaling / self.shift_bucket_size).ceil();
        let cap = (self.max_buckets / 2) as f64;
        // Both counts are finite: the validated sizes are positive and finite.
        if scaling_half > cap || rt_half > cap || scaling_half < 0.0 || !rt_half.is_finite() {
            return Err(Error::InvalidValue(
                "pose clustering histogram exceeds its bucket limit".into(),
            ));
        }
        let (scaling_half, rt_half) = (scaling_half as usize, rt_half as usize);
        let mut hashes = Hashes {
            scaling_1: Histogram::new(
                2 * scaling_half + 1,
                self.scaling_bucket_size,
                scaling_half as f64,
                0.0,
            ),
            scaling_2: Histogram::new(
                2 * scaling_half + 1,
                self.scaling_bucket_size,
                scaling_half as f64,
                0.0,
            ),
            rt_low: Histogram::new(
                1 + 2 * rt_half,
                self.shift_bucket_size,
                rt_half as f64,
                rt_low,
            ),
            rt_high: Histogram::new(
                1 + 2 * rt_half,
                self.shift_bucket_size,
                rt_half as f64,
                rt_high,
            ),
        };
        actual += 1;
        progress.set(actual)?;

        // Step 3: the ratio of the total intensities (computeIntensityRatio).
        let total = |map: &[Peak2D]| {
            let mut sum = 0.0;
            for point in map {
                sum += f64::from(point.intensity);
            }
            sum
        };
        let ratio = total(&model) / total(&scene);
        progress.set(20)?;
        // The source's function-local `static Int dump_buckets_serial`.
        let serial = AFFINE_DUMP_SERIAL.fetch_add(1, Ordering::Relaxed) + 1;
        if !ratio.is_finite() {
            return Err(Error::InvalidValue(
                "the total intensity ratio of the maps is not finite".into(),
            ));
        }
        let open = |base: &str, infix: &str| -> Result<Dump> {
            if base.is_empty() {
                Ok(Dump::off())
            } else {
                Dump::create(&format!("{base}{infix}{serial}"))
            }
        };

        // Step 4: two rounds of hashing.
        let context = Hashing {
            model: &model,
            scene: &scene,
            rt_pair_min_distance,
            mz_pair_max_distance: self.mz_pair_max_distance,
            ratio,
            rt_low,
            rt_high,
            max_quadruples: self.max_quadruples,
        };
        // Both rounds open (and so truncate) the same pairs file; only the
        // second writes pairs.
        let mut pairs = open(&self.dump_pairs, "_phase_two_")?;
        pairs.write("# i j k l \n")?;
        context.hash(&mut hashes, None, &mut pairs)?;
        pairs.close()?;
        progress.set(30)?;

        // scalingEstimate
        let mut scale_dump = open(&self.dump_buckets, "_scale_")?;
        if scale_dump.is_open() {
            scale_dump.write(&format!(
                "# rt scale hash table buckets dump ( scale, height ) : {}_scale_{serial}\n# unfiltered hash data\n",
                self.dump_buckets
            ))?;
        }
        let histogram = &mut hashes.scaling_1;
        histogram.dump_rows(&mut scale_dump, 0, true)?;
        scale_dump.write("\n")?;
        histogram.tophat();
        scale_dump.write("# tophat filtered hash data\n")?;
        histogram.dump_rows(&mut scale_dump, 1, true)?;
        scale_dump.write("\n")?;
        let cutoff = histogram.cut_noise()?;
        scale_dump.write(&format!("# after freq_cutoff, which is: {}\n", g(cutoff)))?;
        histogram.dump_rows(&mut scale_dump, 2, true)?;
        scale_dump.write("\n")?;
        let steps = histogram.centroid_steps(0.0, 1.0)?;
        for (round, step) in steps.iter().enumerate() {
            let (mean, sd) = (step.mean, step.stdev);
            scale_dump.write(&format!(
                "# loop: {round}  mean: {} [{}]  stdev: {} [{}]  (mean-stdev): {} [{}]  (mean+stdev): {} [{}]  data_range_begin: {}  data_range_end: {}\n",
                g(mean),
                g(mean.exp()),
                g(sd),
                g(mean.exp()),
                g(mean - sd),
                g((mean - sd).exp()),
                g(mean + sd),
                g((mean + sd).exp()),
                step.begin,
                step.end
            ))?;
        }
        scale_dump.write("# EOF\n")?;
        scale_dump.close()?;
        let (log_mean, log_stdev) = steps.last().map_or((0.0, 0.0), |s| (s.mean, s.stdev));
        let scale_low = (log_mean - log_stdev).exp();
        let scale_high = (log_mean + log_stdev).exp();
        progress.set(40)?;

        let mut pairs = open(&self.dump_pairs, "_phase_two_")?;
        pairs.write("# i j k l \n")?;
        context.hash(&mut hashes, Some((scale_low, scale_high)), &mut pairs)?;
        pairs.close()?;
        progress.set(50)?;

        // shiftEstimate
        let mut ends = [
            (
                &mut hashes.rt_low,
                open(&self.dump_buckets, "_low_")?,
                "low",
                0.0,
                1.0,
            ),
            (
                &mut hashes.rt_high,
                open(&self.dump_buckets, "_high_")?,
                "high",
                1.0,
                2.0,
            ),
        ];
        for (histogram, dump, end, _, _) in ends.iter_mut() {
            if dump.is_open() {
                dump.write(&format!(
                    "# rt {end} hash table buckets dump ( scale, height ) : {}_{end}_{serial}\n# unfiltered hash data\n",
                    self.dump_buckets
                ))?;
            }
            histogram.dump_rows(dump, 0, true)?;
            dump.write("\n")?;
        }
        for (histogram, dump, _, _, _) in ends.iter_mut() {
            histogram.tophat();
            dump.write("# tophat filtered hash data\n")?;
            histogram.dump_rows(dump, 1, true)?;
            dump.write("\n")?;
        }
        for (histogram, dump, _, _, _) in ends.iter_mut() {
            let cutoff = histogram.cut_noise()?;
            dump.write(&format!("# after freq_cutoff, which is: {}\n", g(cutoff)))?;
            histogram.dump_rows(dump, 2, true)?;
            dump.write("\n")?;
        }
        let mut images = [0.0; 2];
        for (index, (histogram, dump, _, low_pad, high_pad)) in ends.iter_mut().enumerate() {
            let steps = histogram.centroid_steps(*low_pad, *high_pad)?;
            for (round, step) in steps.iter().enumerate() {
                dump.write(&format!(
                    "# loop: {round}  mean: {}  stdev: {}  (mean-stdev): {}  (mean+stdev): {}  data_range_begin: {}  data_range_end: {}\n",
                    g(step.mean),
                    g(step.stdev),
                    g(step.mean - step.stdev),
                    g(step.mean + step.stdev),
                    step.begin,
                    step.end
                ))?;
            }
            images[index] = steps.last().map_or(0.0, |s| s.mean);
        }
        for (_, dump, _, _, _) in ends.iter_mut() {
            dump.write("# EOF\n")?;
            dump.close()?;
        }
        let [low_image, high_image] = images;
        progress.set(80)?;

        // Step 5: the transformation.
        actual = 50;
        actual += 1;
        progress.set(actual)?;
        let slope = (high_image - low_image) / (rt_high - rt_low);
        let intercept = low_image - rt_low * slope;
        if !slope.is_finite() || !intercept.is_finite() {
            return Err(Error::InvalidValue(
                "Superimposer could not compute an initial transformation!You can try to increase 'max_num_peaks_considered' to solve this.".into(),
            ));
        }
        let mut params = Param::new();
        params.set_value("slope", ParamValue::Float(slope), "", &[])?;
        params.set_value("intercept", ParamValue::Float(intercept), "", &[])?;
        let mut transformation = TransformationDescription::default();
        transformation.fit_model_with_param("linear", &params)?;
        actual += 1;
        progress.set(actual)?;
        Ok(Superimposition {
            transformation,
            warnings,
        })
    }
}

impl BaseSuperimposer for PoseClusteringAffineSuperimposer {
    /// Estimate the affine transformation from `scene` onto `model` (source
    /// `run(const ConsensusMap&, const ConsensusMap&, TransformationDescription&)`):
    /// each consensus feature becomes a point of its RT, m/z and intensity,
    /// and [`PoseClusteringAffineSuperimposer::superimpose_points`] runs on
    /// them. The source's `updateRanges` precondition is not needed, because
    /// the ranges are taken from the points.
    ///
    /// # Errors
    ///
    /// As [`PoseClusteringAffineSuperimposer::superimpose_points`].
    fn superimpose(
        &self,
        model: &ConsensusMap,
        scene: &ConsensusMap,
        progress: &mut ProgressReporter<'_>,
    ) -> Result<Superimposition> {
        let points = |map: &ConsensusMap| -> Vec<Peak2D> {
            map.features
                .iter()
                .map(|f| Peak2D::new(f.rt, f.mz, f.intensity))
                .collect()
        };
        self.superimpose_points(&points(model), &points(scene), progress)
    }
}

/// The four histograms of the affine estimation.
struct Hashes {
    scaling_1: Histogram,
    scaling_2: Histogram,
    rt_low: Histogram,
    rt_high: Histogram,
}

/// The inputs of `affineTransformationHashing` that do not change between
/// the rounds.
struct Hashing<'a> {
    model: &'a [Peak2D],
    scene: &'a [Peak2D],
    rt_pair_min_distance: f64,
    mz_pair_max_distance: f64,
    ratio: f64,
    rt_low: f64,
    rt_high: f64,
    max_quadruples: u64,
}

impl Hashing<'_> {
    /// One round of `affineTransformationHashing`: `bounds` is `None` for the
    /// first round (hash `ln(scaling)` only) and the likely scaling range for
    /// the second (hash `ln(scaling)` and the images of both RT ends for
    /// quadruples within it).
    fn hash(
        &self,
        hashes: &mut Hashes,
        bounds: Option<(f64, f64)>,
        pairs: &mut Dump,
    ) -> Result<()> {
        let (model, scene) = (self.model, self.scene);
        let (n, m) = (model.len(), scene.len());
        let d = self.mz_pair_max_distance;
        let rt = |p: &Peak2D| p.position[0];
        let mz = |p: &Peak2D| p.position[1];
        let mut visited: u64 = 0;
        let (mut i_low, mut i_high, mut k_low, mut k_high) = (0usize, 0usize, 0usize, 0usize);
        for i in 0..n - 1 {
            let mz_i = mz(&model[i]);
            while i_low < n && mz(&model[i_low]) < mz_i - d {
                i_low += 1;
            }
            while i_high < n && mz(&model[i_high]) <= mz_i + d {
                i_high += 1;
            }
            let i_weight = 1.0 / (i_high - i_low) as f64 - WINDOW_LENGTH_BASELINE;
            if i_weight <= 0.0 {
                continue;
            }
            while k_low < m && mz(&scene[k_low]) < mz_i - d {
                k_low += 1;
            }
            while k_high < m && mz(&scene[k_high]) <= mz_i + d {
                k_high += 1;
            }
            for k in k_low..k_high {
                let k_weight = 1.0 / (k_high - k_low) as f64 - WINDOW_LENGTH_BASELINE;
                if k_weight <= 0.0 {
                    continue;
                }
                let mut similarity_ik =
                    similarity(model[i].intensity, scene[k].intensity, self.ratio)?;
                similarity_ik *= i_weight;
                similarity_ik *= k_weight;
                // The j window starts from i's lower bound, and both j and l
                // windows are centred on model[i] and model[j] as in the
                // source (the j window deliberately on model[i]).
                let (mut j_low, mut j_high, mut l_low, mut l_high) = (i_low, i_low, k_low, k_high);
                for j in i + 1..n {
                    let diff_model = rt(&model[j]) - rt(&model[i]);
                    if diff_model.abs() < self.rt_pair_min_distance {
                        continue;
                    }
                    while j_low < n && mz(&model[j_low]) < mz_i - d {
                        j_low += 1;
                    }
                    while j_high < n && mz(&model[j_high]) <= mz_i + d {
                        j_high += 1;
                    }
                    let j_weight = 1.0 / (j_high - j_low) as f64 - WINDOW_LENGTH_BASELINE;
                    if j_weight <= 0.0 {
                        continue;
                    }
                    let mz_j = mz(&model[j]);
                    while l_low < m && mz(&scene[l_low]) < mz_j - d {
                        l_low += 1;
                    }
                    while l_high < m && mz(&scene[l_high]) <= mz_j + d {
                        l_high += 1;
                    }
                    for l in l_low..l_high {
                        visited += 1;
                        if visited > self.max_quadruples {
                            return Err(Error::InvalidValue(
                                "pose clustering exceeds its quadruple limit".into(),
                            ));
                        }
                        let l_weight = 1.0 / (l_high - l_low) as f64 - WINDOW_LENGTH_BASELINE;
                        if l_weight <= 0.0 {
                            continue;
                        }
                        let diff_scene = rt(&scene[l]) - rt(&scene[k]);
                        if diff_scene.abs() < self.rt_pair_min_distance
                            || ((diff_model > 0.0) != (diff_scene > 0.0))
                        {
                            continue;
                        }
                        let scaling = diff_model / diff_scene;
                        let shift = rt(&model[i]) - rt(&scene[k]) * scaling;
                        let mut similarity_jl =
                            similarity(model[j].intensity, scene[l].intensity, self.ratio)?;
                        similarity_jl *= j_weight;
                        similarity_jl *= l_weight;
                        let weight = similarity_ik * similarity_jl;
                        match bounds {
                            None => hashes.scaling_1.add_value(scaling.ln(), weight),
                            Some((low, high)) => {
                                if scaling >= low && scaling <= high {
                                    hashes.scaling_2.add_value(scaling.ln(), weight);
                                    hashes
                                        .rt_low
                                        .add_value(shift + self.rt_low * scaling, weight);
                                    hashes
                                        .rt_high
                                        .add_value(shift + self.rt_high * scaling, weight);
                                    if pairs.is_open() {
                                        pairs.write(&format!(
                                            "{i} {} {} {j} {} {} {k} {} {} {l} {} {} {} \n",
                                            g(rt(&model[i])),
                                            g(mz(&model[i])),
                                            g(rt(&model[j])),
                                            g(mz(&model[j])),
                                            g(rt(&scene[k])),
                                            g(mz(&scene[k])),
                                            g(rt(&scene[l])),
                                            g(mz(&scene[l])),
                                            g(weight)
                                        ))?;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
