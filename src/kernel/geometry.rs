// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Checked RT/m/z geometry and the OpenMS scan-envelope hull representation.

use super::NumericRange;
use crate::{Error, Result};

/// A position in retention time (seconds) and mass-to-charge ratio.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Point2D {
    pub rt: f64,
    pub mz: f64,
}

impl Point2D {
    /// A point at the given retention time and m/z.
    pub const fn new(rt: f64, mz: f64) -> Self {
        Self { rt, mz }
    }

    /// Check that both coordinates are finite.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when either coordinate is not finite.
    pub fn validate(self) -> Result<()> {
        if !self.rt.is_finite() || !self.mz.is_finite() {
            return Err(Error::InvalidValue(
                "point coordinates must be finite".into(),
            ));
        }
        Ok(())
    }
}

/// A nonempty, inclusive rectangular bound; use `Option` for an empty bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundingBox2D {
    rt: NumericRange,
    mz: NumericRange,
}

impl BoundingBox2D {
    /// A bounding box spanning `min` to `max`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a coordinate is not finite or a
    /// minimum exceeds its maximum.
    pub fn new(min: Point2D, max: Point2D) -> Result<Self> {
        min.validate()?;
        max.validate()?;
        if min.rt > max.rt || min.mz > max.mz {
            return Err(Error::InvalidValue(
                "bounding box minima exceed maxima".into(),
            ));
        }
        Ok(Self {
            rt: NumericRange {
                min: min.rt,
                max: max.rt,
            },
            mz: NumericRange {
                min: min.mz,
                max: max.mz,
            },
        })
    }

    /// The retention time extent.
    pub const fn rt_range(self) -> NumericRange {
        self.rt
    }
    /// The mass-to-charge extent.
    pub const fn mz_range(self) -> NumericRange {
        self.mz
    }
    /// The lower corner.
    pub const fn min(self) -> Point2D {
        Point2D::new(self.rt.min, self.mz.min)
    }
    /// The upper corner.
    pub const fn max(self) -> Point2D {
        Point2D::new(self.rt.max, self.mz.max)
    }

    /// Whether `point` lies inside the box, borders included.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a coordinate is not finite.
    pub fn encloses(self, point: Point2D) -> Result<bool> {
        point.validate()?;
        Ok(self.rt.min <= point.rt
            && point.rt <= self.rt.max
            && self.mz.min <= point.mz
            && point.mz <= self.mz.max)
    }

    /// Smallest rectangle containing both boxes, including any gap between them.
    pub fn union(self, other: Self) -> Self {
        Self {
            rt: NumericRange {
                min: self.rt.min.min(other.rt.min),
                max: self.rt.max.max(other.rt.max),
            },
            mz: NumericRange {
                min: self.mz.min.min(other.mz.min),
                max: self.mz.max.max(other.mz.max),
            },
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Scan {
    rt: f64,
    min_mz: f64,
    max_mz: f64,
}

/// OpenMS's hull: one m/z interval per RT scan, linearly joined between scans.
///
/// This representation can be non-convex. [`Self::from_points`] builds scan
/// intervals, while [`Self::set_hull_points`] retains an ordered outline only;
/// as in C++, outline-only hulls cannot answer containment queries. Operations
/// never maintain an observable cache, so reading the outline does not change
/// equality or leave stale points after compression.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ConvexHull2D {
    scans: Vec<Scan>,
    outline: Vec<Point2D>,
}

impl ConvexHull2D {
    /// An empty hull with no scans.
    pub fn new() -> Self {
        Self::default()
    }

    /// Groups equal RT values and retains their inclusive minimum/maximum m/z.
    pub fn from_points(points: &[Point2D]) -> Result<Self> {
        let mut hull = Self::new();
        hull.add_points(points)?;
        Ok(hull)
    }

    /// Whether the hull holds no points.
    pub fn is_empty(&self) -> bool {
        self.scans.is_empty() && self.outline.is_empty()
    }
    /// The number of scans whose envelopes the hull records.
    pub fn scan_count(&self) -> usize {
        self.scans.len()
    }
    /// Allocation-free upper bound before exporting or cloning hull points.
    #[cfg(feature = "featurexml")]
    pub(crate) fn point_count_bound(&self) -> usize {
        if self.scans.is_empty() {
            self.outline.len()
        } else {
            self.scans.len().saturating_mul(2)
        }
    }
    /// Remove every point and scan.
    pub fn clear(&mut self) {
        self.scans.clear();
        self.outline.clear();
    }

    /// Replaces the hull with an ordered outline, without inventing scan data.
    pub fn set_hull_points(&mut self, points: &[Point2D]) -> Result<()> {
        for point in points {
            point.validate()?;
        }
        self.outline = points.to_vec();
        self.scans.clear();
        Ok(())
    }

    /// Lower m/z values in ascending RT, then upper values in descending RT.
    pub fn hull_points(&self) -> Vec<Point2D> {
        if self.scans.is_empty() {
            return self.outline.clone();
        }
        let mut points: Vec<_> = self
            .scans
            .iter()
            .map(|s| Point2D::new(s.rt, s.min_mz))
            .collect();
        for (i, scan) in self.scans.iter().enumerate().rev() {
            if (i == 0 || i == self.scans.len() - 1) && scan.min_mz == scan.max_mz {
                continue;
            }
            points.push(Point2D::new(scan.rt, scan.max_mz));
        }
        points
    }

    /// Extends a scan interval; returns false if this point was already covered
    /// at exactly the same RT. Outline-only hulls must be explicitly cleared
    /// first, preventing the silent outline loss of C++ `addPoint`.
    pub fn add_point(&mut self, point: Point2D) -> Result<bool> {
        point.validate()?;
        self.require_scan_mode()?;
        match self
            .scans
            .binary_search_by(|scan| scan.rt.partial_cmp(&point.rt).unwrap())
        {
            Ok(i) => {
                let scan = &mut self.scans[i];
                let changed = point.mz < scan.min_mz || point.mz > scan.max_mz;
                scan.min_mz = scan.min_mz.min(point.mz);
                scan.max_mz = scan.max_mz.max(point.mz);
                Ok(changed)
            }
            Err(i) => {
                self.scans.insert(
                    i,
                    Scan {
                        rt: point.rt,
                        min_mz: point.mz,
                        max_mz: point.mz,
                    },
                );
                Ok(true)
            }
        }
    }

    /// Validates the complete input before changing the hull; batch sorting
    /// avoids quadratic insertion cost for an unsorted point cloud.
    pub fn add_points(&mut self, points: &[Point2D]) -> Result<()> {
        for point in points {
            point.validate()?;
        }
        if points.is_empty() {
            return Ok(());
        }
        self.require_scan_mode()?;
        let mut scans = self.scans.clone();
        scans.extend(points.iter().map(|p| Scan {
            rt: p.rt,
            min_mz: p.mz,
            max_mz: p.mz,
        }));
        scans.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap());
        scans.dedup_by(|later, earlier| {
            if later.rt == earlier.rt {
                earlier.min_mz = earlier.min_mz.min(later.min_mz);
                earlier.max_mz = earlier.max_mz.max(later.max_mz);
                true
            } else {
                false
            }
        });
        self.scans = scans;
        Ok(())
    }

    /// The enclosing bounding box, or `None` when the hull is empty.
    pub fn bounding_box(&self) -> Option<BoundingBox2D> {
        let mut bounds = None;
        let mut add = |point| {
            let box_ =
                BoundingBox2D::new(point, point).expect("private hull coordinates are finite");
            bounds = Some(bounds.map_or(box_, |old: BoundingBox2D| old.union(box_)));
        };
        if self.scans.is_empty() {
            for &point in &self.outline {
                add(point);
            }
        } else {
            for scan in &self.scans {
                add(Point2D::new(scan.rt, scan.min_mz));
                add(Point2D::new(scan.rt, scan.max_mz));
            }
        }
        bounds
    }

    /// Inclusive source containment. At an exact interior scan, a point outside
    /// that scan is also checked against the interpolation of its two neighbors,
    /// preserving an unusual but observable C++ behavior.
    pub fn encloses(&self, point: Point2D) -> Result<bool> {
        point.validate()?;
        self.require_scan_mode()?;
        let index = self.scans.partition_point(|scan| scan.rt < point.rt);
        let exact = self.scans.get(index).filter(|scan| scan.rt == point.rt);
        if let Some(scan) = exact {
            if scan.min_mz <= point.mz && point.mz <= scan.max_mz {
                return Ok(true);
            }
        }
        let upper = index + usize::from(exact.is_some());
        if index == 0 || upper == self.scans.len() {
            return Ok(false);
        }
        let lower = self.scans[index - 1];
        let higher = self.scans[upper];
        let fraction = (point.rt - lower.rt) / (higher.rt - lower.rt);
        let min = lower.min_mz + fraction * (higher.min_mz - lower.min_mz);
        let max = lower.max_mz + fraction * (higher.max_mz - lower.max_mz);
        if !fraction.is_finite()
            || !min.is_finite()
            || !max.is_finite()
            || !(higher.rt - lower.rt).is_finite()
        {
            return Err(Error::InvalidValue("hull interpolation overflow".into()));
        }
        Ok(min <= point.mz && point.mz <= max)
    }

    /// Removes an interior scan only when its m/z span equals both neighbors.
    pub fn compress(&mut self) -> usize {
        if self.scans.len() < 3 {
            return 0;
        }
        let before = self.scans.len();
        let same_span = |a: Scan, b: Scan| a.min_mz == b.min_mz && a.max_mz == b.max_mz;
        let mut compressed = vec![self.scans[0]];
        for triple in self.scans.windows(3) {
            if !same_span(triple[0], triple[1]) || !same_span(triple[1], triple[2]) {
                compressed.push(triple[1]);
            }
        }
        compressed.push(self.scans[before - 1]);
        self.scans = compressed;
        before - self.scans.len()
    }

    /// Replaces the hull with its rectangle. An empty hull remains empty.
    pub fn expand_to_bounding_box(&mut self) {
        if let Some(bounds) = self.bounding_box() {
            *self = Self::from_bounding_box(bounds);
        }
    }

    /// A hull covering exactly the given bounding box.
    pub fn from_bounding_box(bounds: BoundingBox2D) -> Self {
        let mut scans = vec![Scan {
            rt: bounds.rt.min,
            min_mz: bounds.mz.min,
            max_mz: bounds.mz.max,
        }];
        if bounds.rt.max != bounds.rt.min {
            scans.push(Scan {
                rt: bounds.rt.max,
                min_mz: bounds.mz.min,
                max_mz: bounds.mz.max,
            });
        }
        Self {
            scans,
            outline: Vec::new(),
        }
    }

    fn require_scan_mode(&self) -> Result<()> {
        if !self.outline.is_empty() {
            Err(Error::Unsupported(
                "outline-only hull has no scan intervals; clear or expand it first".into(),
            ))
        } else {
            Ok(())
        }
    }
}
