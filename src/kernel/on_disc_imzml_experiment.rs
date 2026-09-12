// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

#![cfg(feature = "mzml")]

//! Random-access, on-disc reader for imzML mass-spectrometry-imaging datasets.
//!
//! Ports `KERNEL/OnDiscImzMLExperiment.h`. See `docs/ON_DISC_IMZML_SUPPORT.md`.
//!
//! This is the kernel-level façade over the imzML reader in
//! [`imzml_handler`](crate::format::imzml_handler): it opens a dataset once,
//! builds the pixel grid eagerly, and then answers per-pixel questions by
//! seeking into the companion `.ibd`. Nothing but the XML index is held in
//! memory, which is what makes a multi-gigabyte image browsable.
//!
//! An imzML dataset is two files. The `.imzML` is mzML XML carrying metadata
//! and, per spectrum, the byte offset, element count and encoded length of that
//! spectrum's m/z and intensity arrays inside the `.ibd`; the `.ibd` opens with
//! a 16-byte UUID that must equal the XML's `IMS:1000080`. Both storage modes
//! are supported: *continuous*, where every pixel names the same shared m/z
//! offset, and *processed*, where each pixel has its own. The parsing of all
//! that belongs to [`imzml_handler`](crate::format::imzml_handler) and is not
//! repeated here; what this module adds is the grid, the pixel lookup and the
//! ion-image sweep.
//!
//! Three types the source reaches for live in `OpenMS/IMAGING/`, a directory no
//! Rust package owns yet:
//! [`ImagingGeometry`](crate::kernel::on_disc_imzml_experiment::ImagingGeometry)
//! (source `MSImagingGeometry`),
//! [`ImagingRegion`](crate::kernel::on_disc_imzml_experiment::ImagingRegion)
//! (source `MSImagingRegion`) and
//! [`IonImage`](crate::kernel::on_disc_imzml_experiment::IonImage). They are
//! reproduced here because `getGeometry()` and `extractIonImage()` are part of
//! the façade's public surface and cannot be ported without them; they are
//! *not* a claim on those headers, whose own class tests were not used as an
//! oracle. See the support document for what that means for a later IMAGING
//! package.
//!
//! Every offset and length the sweep follows came out of a file, and the pixel
//! count is a product of two file-supplied dimensions, so nothing is allocated
//! before it is checked: the grid is refused above
//! [`MAX_IMAGE_PIXELS`](crate::kernel::on_disc_imzml_experiment::MAX_IMAGE_PIXELS),
//! and every `.ibd` range is preflighted against
//! [`ImzMLReadLimits`](crate::format::imzml_handler::ImzMLReadLimits) and the
//! actual file length by the layer below.
//!
//! Nothing here starts a thread. Every decode takes `&mut self`, so the
//! compiler enforces the exclusive access to the `.ibd` handle that the
//! source's `@note` about thread safety only documents.

use crate::format::imzml_handler::{
    DecodedSpectrum, ImzMLHandler, ImzMLIndex, ImzMLMeta, ImzMLReadLimits, ImzMLSpectrumIndex,
    UuidStatus, infer_ibd_path,
};
use crate::kernel::MSSpectrum;
use crate::kernel::ranges::RangeBase;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// Largest pixel grid this module will declare, allocate or walk.
///
/// 4096 x 4096 cells. The source has no such ceiling: `MSImagingGeometry`
/// stores whatever `setDimensions` is handed, and `IonImage::resize` then
/// allocates `width * height` doubles plus a parallel mask from it. Both
/// dimensions come from the file's `IMS:1000042` / `IMS:1000043` (raised to the
/// largest observed coordinate), so a hostile or corrupt header sizes that
/// allocation. This value is comfortably above any real acquisition — the
/// reader's own `max_spectra` ceiling is 5,000,000 pixels — and matches
/// [`MSSpectrum::MAX_RASTER_PIXELS`](crate::kernel::MSSpectrum::MAX_RASTER_PIXELS),
/// the ceiling the ion-mobility rasterizer applies to the same product-of-two-
/// file-values shape.
pub const MAX_IMAGE_PIXELS: usize = 16_777_216;

/// Largest number of regions one geometry may carry.
///
/// Native ceiling with no source counterpart. Each
/// [`ImagingGeometry::add_region`] tests the newcomer against every region
/// already present, and a masked footprint's test walks the intersection of the
/// two bounding boxes, so an unbounded region list would make insertion
/// quadratic in caller-supplied work.
pub const MAX_REGIONS: usize = 4_096;

/// How many offending spectra a [`GeometryReport`] list names before it stops
/// collecting and only counts.
///
/// Source `buildImagingGeometry` uses exactly this cap (`const Size max_listed
/// = 20`) for its duplicate-pixel warning, with the comment that a broken
/// converter can put every spectrum on the same pixel.
pub const MAX_LISTED_PROBLEMS: usize = 20;

// ---------------------------------------------------------------------------
// Imaging region (source IMAGING/MSImagingRegion.h)
// ---------------------------------------------------------------------------

/// How a region's footprint is represented.
///
/// Source `MSImagingRegion::Shape`.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum RegionShape {
    /// Grid-aligned bounding box; every cell inside the box belongs to the
    /// region.
    Rectangle,
    /// Per-pixel bitmask inside the bounding box.
    Mask,
}

/// A spatial region within an imaging dataset, in global pixel coordinates.
///
/// Source `MSImagingRegion`, a pure-geometry footprint that knows nothing about
/// acquired pixels or spectra and is reusable as a bare annotation. Coordinates
/// are zero-based and bounding boxes are inclusive on both ends. Build one with
/// [`rectangle`](Self::rectangle) or [`from_mask`](Self::from_mask); the source
/// likewise has no public constructor and only those two factories.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ImagingRegion {
    id: usize,
    shape: RegionShape,
    name: String,
    mask: Vec<bool>,
    min_x: u32,
    min_y: u32,
    max_x: u32,
    max_y: u32,
}

impl ImagingRegion {
    /// A rectangular region spanning the inclusive bounding box
    /// `[min_x, max_x] x [min_y, max_y]`.
    ///
    /// # Arguments
    ///
    /// * `id` — region identifier; [`ImagingGeometry::NO_REGION`] is rejected
    ///   by [`ImagingGeometry::add_region`], not here, exactly as in the source.
    /// * `name` — human-readable region name.
    /// * `min_x`, `min_y` — leftmost column and top row, inclusive, zero-based.
    /// * `max_x`, `max_y` — rightmost column and bottom row, inclusive.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidRange`] when `max_x < min_x` or `max_y < min_y`, where
    /// the source raises `Exception::InvalidValue`; the dedicated range variant
    /// says more about the same condition.
    pub fn rectangle(
        id: usize,
        name: impl Into<String>,
        min_x: u32,
        min_y: u32,
        max_x: u32,
        max_y: u32,
    ) -> Result<Self> {
        if max_x < min_x || max_y < min_y {
            return Err(Error::InvalidRange(format!(
                "imaging region coordinate maximum is below the minimum: min=({min_x},{min_y}) max=({max_x},{max_y})"
            )));
        }
        Ok(Self {
            id,
            shape: RegionShape::Rectangle,
            name: name.into(),
            mask: Vec::new(),
            min_x,
            min_y,
            max_x,
            max_y,
        })
    }

    /// A masked region from a row-major bitmask over a bounding box.
    ///
    /// The global bounding box is
    /// `[origin_x, origin_x + width - 1] x [origin_y, origin_y + height - 1]`
    /// and `mask` is stored bounding-box-local, `true` meaning inside.
    ///
    /// # Arguments
    ///
    /// * `id` — region identifier.
    /// * `name` — human-readable region name.
    /// * `origin_x`, `origin_y` — leftmost column and top row of the bounding
    ///   box, zero-based.
    /// * `width`, `height` — bounding-box extent in pixels; both must be
    ///   positive.
    /// * `mask` — row-major bitmask of exactly `width * height` entries, at
    ///   least one of them `true`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `width` or `height` is zero, when
    /// `mask.len()` is not `width * height`, when the mask is all-false, or
    /// when `width * height` exceeds [`MAX_IMAGE_PIXELS`]. The source raises
    /// `Exception::InvalidValue` for the first three and has no ceiling for the
    /// fourth. [`Error::InvalidRange`] when the far corner
    /// `origin + extent - 1` leaves `u32`; the source computes `origin_x +
    /// width - 1` in `UInt` arithmetic, which wraps silently and produces a
    /// bounding box that starts above where it ends.
    pub fn from_mask(
        id: usize,
        name: impl Into<String>,
        origin_x: u32,
        origin_y: u32,
        width: u32,
        height: u32,
        mask: Vec<bool>,
    ) -> Result<Self> {
        if width == 0 || height == 0 {
            return Err(Error::InvalidValue(format!(
                "imaging region extent must be positive: (w,h)=({width},{height})"
            )));
        }
        let cells = u64::from(width) * u64::from(height);
        let cells = usize::try_from(cells)
            .ok()
            .filter(|&cells| cells <= MAX_IMAGE_PIXELS)
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "imaging region covers {cells} cells, above the ceiling of {MAX_IMAGE_PIXELS}"
                ))
            })?;
        if mask.len() != cells {
            return Err(Error::InvalidValue(format!(
                "imaging region mask holds {} entries, not the {cells} its {width}x{height} box needs",
                mask.len()
            )));
        }
        if !mask.iter().any(|&inside| inside) {
            return Err(Error::InvalidValue(
                "imaging region mask selects no pixel at all".into(),
            ));
        }
        let max_x = far_corner(origin_x, width, "column")?;
        let max_y = far_corner(origin_y, height, "row")?;
        Ok(Self {
            id,
            shape: RegionShape::Mask,
            name: name.into(),
            mask,
            min_x: origin_x,
            min_y: origin_y,
            max_x,
            max_y,
        })
    }

    /// The region identifier, as source `getId()`.
    pub fn id(&self) -> usize {
        self.id
    }

    /// The region name, as source `getName()`.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Whether the footprint is a rectangle or a bitmask, as source
    /// `getShape()`.
    pub fn shape(&self) -> RegionShape {
        self.shape
    }

    /// Leftmost column of the global bounding box, inclusive.
    pub fn min_x(&self) -> u32 {
        self.min_x
    }

    /// Top row of the global bounding box, inclusive.
    pub fn min_y(&self) -> u32 {
        self.min_y
    }

    /// Rightmost column of the global bounding box, inclusive.
    pub fn max_x(&self) -> u32 {
        self.max_x
    }

    /// Bottom row of the global bounding box, inclusive.
    pub fn max_y(&self) -> u32 {
        self.max_y
    }

    /// Bounding-box width, `max_x - min_x + 1`.
    pub fn bbox_width(&self) -> u32 {
        self.max_x - self.min_x + 1
    }

    /// Bounding-box height, `max_y - min_y + 1`.
    pub fn bbox_height(&self) -> u32 {
        self.max_y - self.min_y + 1
    }

    /// The bounding-box-local row-major bitmask, empty for a
    /// [`RegionShape::Rectangle`].
    ///
    /// Branch on [`shape`](Self::shape) before reading it, as the source's own
    /// documentation instructs.
    pub fn mask(&self) -> &[bool] {
        &self.mask
    }

    /// Number of pixels inside the region: the bounding-box area for a
    /// rectangle, the set-bit count for a mask.
    ///
    /// The source returns `Size`; `usize` here cannot overflow because both
    /// factors are `u32` and the product is formed in `u64`.
    pub fn area(&self) -> usize {
        match self.shape {
            RegionShape::Rectangle => {
                let area = u64::from(self.bbox_width()) * u64::from(self.bbox_height());
                usize::try_from(area).unwrap_or(usize::MAX)
            }
            RegionShape::Mask => self.mask.iter().filter(|&&inside| inside).count(),
        }
    }

    /// Whether the global coordinate `(x, y)` lies inside the footprint: the
    /// bounding box for a rectangle, a set bit for a mask.
    pub fn contains(&self, x: u32, y: u32) -> bool {
        if x < self.min_x || x > self.max_x || y < self.min_y || y > self.max_y {
            return false;
        }
        if self.shape == RegionShape::Rectangle {
            return true;
        }
        let row = u64::from(y - self.min_y);
        let index = row * u64::from(self.bbox_width()) + u64::from(x - self.min_x);
        usize::try_from(index)
            .ok()
            .and_then(|index| self.mask.get(index))
            .copied()
            .unwrap_or(false)
    }

    /// Whether this footprint geometrically overlaps `other`, i.e. whether some
    /// global coordinate lies inside both.
    ///
    /// Symmetric and independent of any acquired pixels, as the source's
    /// documentation states. Two rectangles answer from their bounding boxes
    /// alone; any pair involving a mask walks the intersection of the two
    /// boxes, which is bounded by the smaller box and therefore by
    /// [`MAX_IMAGE_PIXELS`], since a mask's box is exactly its mask.
    pub fn intersects(&self, other: &Self) -> bool {
        if self.max_x < other.min_x
            || self.max_y < other.min_y
            || self.min_x > other.max_x
            || self.min_y > other.max_y
        {
            return false;
        }
        if self.shape == RegionShape::Rectangle && other.shape == RegionShape::Rectangle {
            return true;
        }
        let lo_x = self.min_x.max(other.min_x);
        let hi_x = self.max_x.min(other.max_x);
        let lo_y = self.min_y.max(other.min_y);
        let hi_y = self.max_y.min(other.max_y);
        for y in lo_y..=hi_y {
            for x in lo_x..=hi_x {
                if self.contains(x, y) && other.contains(x, y) {
                    return true;
                }
            }
        }
        false
    }
}

// ---------------------------------------------------------------------------
// Imaging geometry (source IMAGING/MSImagingGeometry.h)
// ---------------------------------------------------------------------------

/// One pixel of the imaging grid, linked to one spectrum of the dataset.
///
/// Source `MSImagingGeometry::Pixel`.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct ImagingPixel {
    /// Column index, zero-based.
    pub x: u32,
    /// Row index, zero-based.
    pub y: u32,
    /// Position of the bound spectrum in the dataset's index.
    pub spectrum_index: usize,
}

/// Pixel-grid metadata plus the `(x, y) -> spectrum index` lookup.
///
/// Source `MSImagingGeometry`. Coordinates are **zero-based**: imzML files are
/// one-based and the loader normalises them, so `IMS:1000050` value 1 becomes
/// column 0. Three-dimensional imaging is deliberately not modelled — the
/// source's own note says serial sections belong in a collection of
/// experiments, one per section — so only the `z == 1` plane of a dataset
/// reaches this grid.
///
/// The default geometry is empty with 1.0 x 1.0 micrometer pixels, which is the
/// source's member-initialiser state and what [`clear`](Self::clear) restores.
#[derive(Clone, Debug, PartialEq)]
pub struct ImagingGeometry {
    width: u32,
    height: u32,
    pixel_size_x: f64,
    pixel_size_y: f64,
    pixel_size_unit: String,
    pixels: Vec<ImagingPixel>,
    lookup: BTreeMap<u64, usize>,
    regions: Vec<ImagingRegion>,
    region_positions: BTreeMap<usize, usize>,
}

impl Default for ImagingGeometry {
    fn default() -> Self {
        Self {
            width: 0,
            height: 0,
            pixel_size_x: 1.0,
            pixel_size_y: 1.0,
            pixel_size_unit: "micrometer".into(),
            pixels: Vec::new(),
            lookup: BTreeMap::new(),
            regions: Vec::new(),
            region_positions: BTreeMap::new(),
        }
    }
}

impl ImagingGeometry {
    /// The identifier a region may not use, as source
    /// `MSImagingGeometry::NO_REGION`.
    ///
    /// The source returns this sentinel from `regionOf()` for a coordinate that
    /// belongs to no region; [`region_of`](Self::region_of) returns `None`
    /// instead, so the value survives here only as the identifier
    /// [`add_region`](Self::add_region) refuses.
    pub const NO_REGION: usize = usize::MAX;

    /// An empty geometry, as the source's default constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the image dimensions.
    ///
    /// # Arguments
    ///
    /// * `width` — number of columns.
    /// * `height` — number of rows.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `width * height` exceeds
    /// [`MAX_IMAGE_PIXELS`]. The source accepts any pair and defers the cost to
    /// `IonImage::resize`, which allocates the product; because both values
    /// come from the file, that allocation is refused here instead. A dataset
    /// declaring a grid above the ceiling therefore fails
    /// [`OnDiscImzMLExperiment::open`] rather than opening and failing at the
    /// first extraction.
    pub fn set_dimensions(&mut self, width: u32, height: u32) -> Result<()> {
        grid_cells(width, height)?;
        self.width = width;
        self.height = height;
        Ok(())
    }

    /// Number of columns, as source `getWidth()`.
    pub fn width(&self) -> u32 {
        self.width
    }

    /// Number of rows, as source `getHeight()`.
    pub fn height(&self) -> u32 {
        self.height
    }

    /// Record the physical pixel size and its unit.
    ///
    /// # Arguments
    ///
    /// * `x`, `y` — pixel extent along each axis, in `unit`.
    /// * `unit` — length unit; the source defaults it to `"micrometer"` and
    ///   every caller inside OpenMS passes that, since imzML's `IMS:1000046` /
    ///   `IMS:1000047` are defined in micrometres. Rust has no default
    ///   arguments, so pass it explicitly.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when either extent is not finite. The source
    /// stores any `double`, including a NaN that would then propagate into
    /// every physical-distance computation a caller performs.
    pub fn set_pixel_size(&mut self, x: f64, y: f64, unit: impl Into<String>) -> Result<()> {
        finite(x, "pixel size along x")?;
        finite(y, "pixel size along y")?;
        self.pixel_size_x = x;
        self.pixel_size_y = y;
        self.pixel_size_unit = unit.into();
        Ok(())
    }

    /// Physical pixel extent along x, as source `getPixelSizeX()`.
    pub fn pixel_size_x(&self) -> f64 {
        self.pixel_size_x
    }

    /// Physical pixel extent along y, as source `getPixelSizeY()`.
    pub fn pixel_size_y(&self) -> f64 {
        self.pixel_size_y
    }

    /// Unit of the pixel size, as source `getPixelSizeUnit()`.
    pub fn pixel_size_unit(&self) -> &str {
        &self.pixel_size_unit
    }

    /// Bind the spectrum at `spectrum_index` to the pixel `(x, y)`.
    ///
    /// # Arguments
    ///
    /// * `x`, `y` — zero-based column and row.
    /// * `spectrum_index` — position of the spectrum in the dataset's index.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] on a duplicate coordinate, when dimensions have
    /// been set and `(x, y)` is outside `[0, width) x [0, height)`, or when the
    /// grid already holds [`MAX_IMAGE_PIXELS`] pixels. The source throws
    /// `Exception::InvalidValue` for the first two and has no third condition.
    ///
    /// Note that the bounds test is skipped entirely while either dimension is
    /// zero, which is the source's behaviour and is what lets
    /// [`build_imaging_geometry`] insert pixels first and derive the extent from
    /// them afterwards.
    pub fn add_pixel(&mut self, x: u32, y: u32, spectrum_index: usize) -> Result<()> {
        if self.pixels.len() >= MAX_IMAGE_PIXELS {
            return Err(Error::InvalidValue(format!(
                "imaging geometry already holds the maximum of {MAX_IMAGE_PIXELS} pixels"
            )));
        }
        if self.width > 0 && self.height > 0 && (x >= self.width || y >= self.height) {
            return Err(Error::InvalidValue(format!(
                "pixel coordinate ({x},{y}) is outside the configured {}x{} geometry",
                self.width, self.height
            )));
        }
        let key = pack_key(x, y);
        if self.lookup.contains_key(&key) {
            return Err(Error::InvalidValue(format!(
                "duplicate pixel coordinate ({x},{y})"
            )));
        }
        self.lookup.insert(key, self.pixels.len());
        self.pixels.push(ImagingPixel {
            x,
            y,
            spectrum_index,
        });
        Ok(())
    }

    /// Whether a pixel was inserted at `(x, y)`, as source `hasPixel()`.
    pub fn has_pixel(&self, x: u32, y: u32) -> bool {
        self.lookup.contains_key(&pack_key(x, y))
    }

    /// The spectrum index recorded for the pixel `(x, y)`, or `None` when no
    /// pixel was inserted there.
    ///
    /// Source `getSpectrumIndex` throws `Exception::ElementNotFound` for an
    /// absent coordinate. An absent pixel is the normal state of most of a
    /// tissue image, so it is reported as `None` rather than as an error;
    /// [`has_pixel`](Self::has_pixel) answers the same question without the
    /// lookup result.
    pub fn spectrum_index(&self, x: u32, y: u32) -> Option<usize> {
        self.lookup
            .get(&pack_key(x, y))
            .and_then(|&position| self.pixels.get(position))
            .map(|pixel| pixel.spectrum_index)
    }

    /// Every bound pixel, in insertion order, as source `getPixels()`.
    pub fn pixels(&self) -> &[ImagingPixel] {
        &self.pixels
    }

    /// How many pixels carry a bound spectrum, as source
    /// `getNumberOfPixels()`.
    pub fn number_of_pixels(&self) -> usize {
        self.pixels.len()
    }

    /// Whether the grid holds no pixel at all.
    ///
    /// Native convenience; the source expresses it as
    /// `getNumberOfPixels() == 0`.
    pub fn is_empty(&self) -> bool {
        self.pixels.is_empty()
    }

    /// Reset dimensions, pixel size, pixels, lookup and regions to the default
    /// state, as source `clear()`.
    pub fn clear(&mut self) {
        *self = Self::default();
    }

    /// Add a region as a decoupled overlay.
    ///
    /// Membership of acquired pixels is derived on demand through
    /// [`ImagingRegion::contains`]; no per-pixel region state is stored, which
    /// is the source's design and the reason a region may be added long after
    /// the pixels.
    ///
    /// # Arguments
    ///
    /// * `region` — footprint to add, moved in. The source copies it.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the identifier equals [`Self::NO_REGION`],
    /// when it is already present, when the footprint geometrically overlaps an
    /// existing region, or when the geometry already holds [`MAX_REGIONS`]
    /// regions. The source throws `Exception::InvalidValue` for the first three
    /// and has no fourth condition.
    pub fn add_region(&mut self, region: ImagingRegion) -> Result<()> {
        if region.id == Self::NO_REGION {
            return Err(Error::InvalidValue(
                "the NO_REGION sentinel is not a valid region identifier".into(),
            ));
        }
        if self.regions.len() >= MAX_REGIONS {
            return Err(Error::InvalidValue(format!(
                "imaging geometry already holds the maximum of {MAX_REGIONS} regions"
            )));
        }
        if self.region_positions.contains_key(&region.id) {
            return Err(Error::InvalidValue(format!(
                "duplicate imaging region identifier {}",
                region.id
            )));
        }
        if self
            .regions
            .iter()
            .any(|existing| existing.intersects(&region))
        {
            return Err(Error::InvalidValue(format!(
                "imaging regions must be disjoint; region {} overlaps an existing one",
                region.id
            )));
        }
        self.region_positions.insert(region.id, self.regions.len());
        self.regions.push(region);
        Ok(())
    }

    /// Remove the region with identifier `id`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no region carries that identifier, where
    /// the source throws `Exception::ElementNotFound`.
    pub fn remove_region(&mut self, id: usize) -> Result<()> {
        let position = self
            .region_positions
            .get(&id)
            .copied()
            .ok_or_else(|| unknown_region(id))?;
        self.regions.remove(position);
        self.region_positions.clear();
        for (position, region) in self.regions.iter().enumerate() {
            self.region_positions.insert(region.id, position);
        }
        Ok(())
    }

    /// Remove every region, leaving the acquired pixels untouched, as source
    /// `clearRegions()`.
    pub fn clear_regions(&mut self) {
        self.regions.clear();
        self.region_positions.clear();
    }

    /// Every region, in insertion order, as source `getRegions()`.
    pub fn regions(&self) -> &[ImagingRegion] {
        &self.regions
    }

    /// Whether a region with identifier `id` exists.
    ///
    /// Native convenience, so that a caller can ask without handling an error;
    /// the source has no such member.
    pub fn has_region(&self, id: usize) -> bool {
        self.region_positions.contains_key(&id)
    }

    /// The region with identifier `id`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when no region carries that identifier, where
    /// the source throws `Exception::ElementNotFound`. This one keeps the error
    /// rather than returning an `Option`, because the region-scoped ion-image
    /// extraction has to distinguish an unknown region from an empty one and
    /// the upstream suite asserts that it does.
    pub fn region(&self, id: usize) -> Result<&ImagingRegion> {
        self.region_positions
            .get(&id)
            .and_then(|&position| self.regions.get(position))
            .ok_or_else(|| unknown_region(id))
    }

    /// How many regions the overlay holds, as source `getNumberOfRegions()`.
    pub fn number_of_regions(&self) -> usize {
        self.regions.len()
    }

    /// The identifier of the region owning the *acquired* pixel `(x, y)`, or
    /// `None` when `(x, y)` carries no spectrum or belongs to no region.
    ///
    /// Source `regionOf` returns the [`Self::NO_REGION`] sentinel for both of
    /// those cases; `None` says the same thing without a magic value. The
    /// source's "acquired" precondition is kept: a coordinate inside a region's
    /// footprint but with no bound spectrum answers `None`. The first matching
    /// region in insertion order wins, which regions being disjoint makes
    /// unambiguous.
    pub fn region_of(&self, x: u32, y: u32) -> Option<usize> {
        if !self.has_pixel(x, y) {
            return None;
        }
        self.regions
            .iter()
            .find(|region| region.contains(x, y))
            .map(ImagingRegion::id)
    }

    /// Positions in [`pixels`](Self::pixels) of the acquired pixels belonging
    /// to region `id`.
    ///
    /// # Errors
    ///
    /// As [`region`](Self::region).
    pub fn region_pixels(&self, id: usize) -> Result<Vec<usize>> {
        let region = self.region(id)?;
        Ok(self
            .pixels
            .iter()
            .enumerate()
            .filter(|(_, pixel)| region.contains(pixel.x, pixel.y))
            .map(|(position, _)| position)
            .collect())
    }

    /// Spectrum indices of the acquired pixels belonging to region `id`.
    ///
    /// # Errors
    ///
    /// As [`region`](Self::region).
    pub fn region_spectrum_indices(&self, id: usize) -> Result<Vec<usize>> {
        let region = self.region(id)?;
        Ok(self
            .pixels
            .iter()
            .filter(|pixel| region.contains(pixel.x, pixel.y))
            .map(|pixel| pixel.spectrum_index)
            .collect())
    }
}

// ---------------------------------------------------------------------------
// Ion image (source IMAGING/IonImage.h)
// ---------------------------------------------------------------------------

/// A dense `width x height` grid of ion intensities with a per-pixel mask.
///
/// Source `IonImage`. Storage is row-major, `index = y * width + x`. Every cell
/// starts masked out; [`set_intensity`](Self::set_intensity) marks one present.
/// The m/z window the image was extracted from travels with the data for
/// traceability. Three-dimensional imaging is deliberately not modelled, as in
/// the source: a serial-section experiment is a collection of independent
/// images, one per section.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IonImage {
    width: u32,
    height: u32,
    intensities: Vec<f64>,
    mask: Vec<bool>,
    mz_range: RangeBase,
}

impl IonImage {
    /// A zero-initialised `width x height` image with every pixel invalid.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `width * height` exceeds
    /// [`MAX_IMAGE_PIXELS`] or the allocation fails. The source's constructor
    /// allocates the product unconditionally.
    pub fn new(width: u32, height: u32) -> Result<Self> {
        let mut image = Self::default();
        image.resize(width, height)?;
        Ok(image)
    }

    /// Resize and zero-initialise; every pixel becomes invalid again.
    ///
    /// # Errors
    ///
    /// As [`new`](Self::new). The image is left unchanged when the new size is
    /// refused, so a rejected resize cannot discard the current contents.
    pub fn resize(&mut self, width: u32, height: u32) -> Result<()> {
        let cells = grid_cells(width, height)?;
        let mut intensities = Vec::new();
        intensities
            .try_reserve_exact(cells)
            .map_err(|_| Error::InvalidValue(format!("cannot allocate {cells} image pixels")))?;
        intensities.resize(cells, 0.0);
        let mut mask = Vec::new();
        mask.try_reserve_exact(cells)
            .map_err(|_| Error::InvalidValue(format!("cannot allocate {cells} image mask bits")))?;
        mask.resize(cells, false);
        self.width = width;
        self.height = height;
        self.intensities = intensities;
        self.mask = mask;
        Ok(())
    }

    /// Number of columns, as source `getWidth()`.
    pub fn width(&self) -> u32 {
        self.width
    }

    /// Number of rows, as source `getHeight()`.
    pub fn height(&self) -> u32 {
        self.height
    }

    /// Whether the cell at `(x, y)` has been written.
    ///
    /// `false` for an out-of-bounds coordinate and for a cell never written, as
    /// in the source: this accessor does not report a range error.
    pub fn has_pixel(&self, x: u32, y: u32) -> bool {
        self.linear_index(x, y)
            .ok()
            .and_then(|index| self.mask.get(index))
            .copied()
            .unwrap_or(false)
    }

    /// The intensity stored at `(x, y)`, `0.0` when the cell was never written.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an out-of-bounds coordinate, where the
    /// source throws `Exception::IndexOverflow`.
    pub fn intensity(&self, x: u32, y: u32) -> Result<f64> {
        let index = self.linear_index(x, y)?;
        Ok(self.intensities[index])
    }

    /// Store `intensity` at `(x, y)` and mark the cell valid.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an out-of-bounds coordinate, where the
    /// source throws `Exception::IndexOverflow`, and also when `intensity` is
    /// not finite. The source stores any `double`; refusing a NaN here means a
    /// `.ibd` holding one cannot silently produce a pixel that compares unequal
    /// to itself.
    pub fn set_intensity(&mut self, x: u32, y: u32, intensity: f64) -> Result<()> {
        finite(intensity, "image pixel intensity")?;
        let index = self.linear_index(x, y)?;
        self.intensities[index] = intensity;
        self.mask[index] = true;
        Ok(())
    }

    /// Record the m/z window the image was extracted from, as source
    /// `setMzRange()`.
    pub fn set_mz_range(&mut self, range: RangeBase) {
        self.mz_range = range;
    }

    /// The m/z window the image was extracted from; empty when never set.
    ///
    /// Source `getMzRange()` returns a `RangeMZ`, the m/z specialisation of
    /// `RangeBase`; this crate's [`RangeBase`] is that base type and carries
    /// the same empty-by-default state.
    pub fn mz_range(&self) -> &RangeBase {
        &self.mz_range
    }

    /// The raw row-major intensity buffer of `width * height` entries, as
    /// source `getData()`.
    pub fn data(&self) -> &[f64] {
        &self.intensities
    }

    /// The parallel pixel mask, indexed exactly as [`data`](Self::data), as
    /// source `getMask()`.
    pub fn mask(&self) -> &[bool] {
        &self.mask
    }

    /// Row-major offset of `(x, y)`, checked against the image extent.
    fn linear_index(&self, x: u32, y: u32) -> Result<usize> {
        if x >= self.width || y >= self.height {
            return Err(Error::InvalidValue(format!(
                "image coordinate ({x},{y}) is outside the {}x{} image",
                self.width, self.height
            )));
        }
        let index = u64::from(y) * u64::from(self.width) + u64::from(x);
        usize::try_from(index)
            .map_err(|_| Error::InvalidValue("image cell offset exceeds usize".into()))
    }
}

// ---------------------------------------------------------------------------
// Geometry construction
// ---------------------------------------------------------------------------

/// What [`build_imaging_geometry`] had to leave out of the grid.
///
/// Source `ImzMLFile::buildImagingGeometry` writes each of these to
/// `OPENMS_LOG_WARN` and loads the dataset anyway; the header's `open()`
/// documentation lists them as warnings that do not stop the load. Returning
/// them keeps that policy while letting a caller that wants strict conformance
/// enforce it — the same choice
/// [`UuidStatus`] makes one layer
/// down. Each list names at most [`MAX_LISTED_PROBLEMS`] spectra, as the
/// source's warning does, while the counts cover every occurrence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct GeometryReport {
    /// Spectra dropped because their coordinates are not both `>= 1`. imzML
    /// coordinates are one-based, so a zero is non-conformant.
    pub non_positive_coordinates: Vec<usize>,
    /// How many spectra were dropped for that reason.
    pub non_positive_count: usize,
    /// Spectra dropped because their pixel lies outside the declared
    /// `IMS:1000042` x `IMS:1000043` grid.
    pub out_of_grid: Vec<usize>,
    /// How many spectra were dropped for that reason.
    pub out_of_grid_count: usize,
    /// Spectra dropped because an earlier spectrum already claimed their pixel.
    /// Only the first spectrum per pixel enters the grid; every spectrum stays
    /// reachable by index.
    pub duplicate_pixels: Vec<usize>,
    /// How many spectra were dropped for that reason.
    pub duplicate_count: usize,
    /// How many spectra were skipped because their `IMS:1000052` is not 1. The
    /// grid is two-dimensional, so only the first plane is mapped. The source
    /// skips these silently, without even a warning.
    pub other_plane_count: usize,
}

impl GeometryReport {
    /// Whether every spectrum of the dataset reached the grid.
    ///
    /// Native convenience: the source has no report to summarise.
    pub fn is_clean(&self) -> bool {
        self.non_positive_count == 0
            && self.out_of_grid_count == 0
            && self.duplicate_count == 0
            && self.other_plane_count == 0
    }
}

/// Build the two-dimensional imaging grid from a parsed imzML index.
///
/// Source `ImzMLFile::buildImagingGeometry(const std::vector<ImzMLSpectrumIndex>&,
/// const ImzMLMeta&, MSImagingGeometry&)`, which its own comment calls the
/// source-of-truth builder shared by the in-memory and on-disc paths so the two
/// cannot diverge. Coordinates come straight from the index, not from
/// `imzml:x` / `imzml:y` meta values.
///
/// The declared `IMS:1000042` / `IMS:1000043` extent is applied first, so that
/// a pixel outside it can be recognised; spectra are then walked in document
/// order, one-based coordinates are converted to zero-based, and the extent is
/// finally raised to cover the pixels that were inserted when the file declared
/// none. A positive `IMS:1000046` / `IMS:1000047` pixel size is copied in, in
/// micrometres.
///
/// # Arguments
///
/// * `index` — per-spectrum entries in document order, as
///   [`ImzMLHandler::index`](crate::format::imzml_handler::ImzMLHandler::index)
///   returns them.
/// * `meta` — dataset metadata, for the declared grid and pixel size.
///
/// # Returns
///
/// The grid and a [`GeometryReport`] naming every spectrum it could not place.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the declared or derived grid exceeds
/// [`MAX_IMAGE_PIXELS`], or when the index holds more placeable spectra than
/// that. The source has no ceiling and no failure mode here at all.
pub fn build_imaging_geometry(
    index: &[ImzMLSpectrumIndex],
    meta: &ImzMLMeta,
) -> Result<(ImagingGeometry, GeometryReport)> {
    let mut geometry = ImagingGeometry::new();
    let mut report = GeometryReport::default();

    let mut width = meta.max_count_x;
    let mut height = meta.max_count_y;
    if width > 0 && height > 0 {
        geometry.set_dimensions(width, height)?;
    }

    let mut max_x = 0_u32;
    let mut max_y = 0_u32;
    for (position, entry) in index.iter().enumerate() {
        if entry.z != 1 {
            // The grid is 2-D: only the first plane is mapped. The source skips
            // these without a warning; the count is reported instead.
            report.other_plane_count += 1;
            continue;
        }
        if entry.x < 1 || entry.y < 1 {
            report.non_positive_count += 1;
            push_capped(&mut report.non_positive_coordinates, position);
            continue;
        }
        let x = entry.x - 1;
        let y = entry.y - 1;
        if width > 0 && height > 0 && (x >= width || y >= height) {
            // The geometry's own rule is an error; the loader softens it to a
            // warning so one stray pixel cannot abort a whole dataset. The
            // spectrum stays reachable by index.
            report.out_of_grid_count += 1;
            push_capped(&mut report.out_of_grid, position);
            continue;
        }
        max_x = max_x.max(x);
        max_y = max_y.max(y);
        if geometry.has_pixel(x, y) {
            report.duplicate_count += 1;
            push_capped(&mut report.duplicate_pixels, position);
            continue;
        }
        geometry.add_pixel(x, y, position)?;
    }

    if width == 0 && (max_x > 0 || geometry.number_of_pixels() > 0) {
        width = max_x
            .checked_add(1)
            .ok_or_else(|| Error::InvalidValue("derived image width overflows u32".into()))?;
    }
    if height == 0 && (max_y > 0 || geometry.number_of_pixels() > 0) {
        height = max_y
            .checked_add(1)
            .ok_or_else(|| Error::InvalidValue("derived image height overflows u32".into()))?;
    }
    if width > 0 && height > 0 && (geometry.width() != width || geometry.height() != height) {
        geometry.set_dimensions(width, height)?;
    }

    if meta.pixel_size_x > 0.0 && meta.pixel_size_y > 0.0 {
        geometry.set_pixel_size(meta.pixel_size_x, meta.pixel_size_y, "micrometer")?;
    }

    Ok((geometry, report))
}

// ---------------------------------------------------------------------------
// The façade
// ---------------------------------------------------------------------------

/// Whether the `.ibd` is open, and where the index lives while it is not.
#[derive(Debug)]
enum Backing {
    /// The dataset is open: the handler owns the index, the metadata and the
    /// `.ibd` handle.
    Open(Box<ImzMLHandler>),
    /// No `.ibd` handle: either nothing was opened yet or `close()` released
    /// it. The index and metadata parsed at `open()` are retained, as the
    /// source's `close()` retains `index_` and `meta_`.
    Closed(Box<ImzMLIndex>),
}

impl Default for Backing {
    fn default() -> Self {
        Self::Closed(Box::default())
    }
}

/// Random-access, on-disc reader for an imzML imaging dataset.
///
/// Source `OnDiscImzMLExperiment`. The counterpart of `OnDiscMSExperiment` for
/// indexed mzML, built for imzML's two-file layout: [`open`](Self::open) parses
/// the XML index and builds the pixel grid, and no peak array is touched until
/// [`spectrum`](Self::spectrum) or an extraction asks for one.
///
/// # Ownership and copying
///
/// The source deletes the copy constructor and copy assignment and defines move
/// construction and move assignment, which leave the moved-from object holding
/// a fresh empty `Impl` so that it stays usable rather than null. Rust moves by
/// default and a moved-from binding is not usable at all, so there is nothing to
/// reset and no move constructor to write; the source's "fresh `Impl`" step has
/// no counterpart here. This type is deliberately **not** [`Clone`]: it owns an
/// open file handle, and duplicating one would hand two readers the same
/// independent seek position. Clone the pieces that are values —
/// [`meta`](Self::meta), [`geometry`](Self::geometry) — when you need a copy.
///
/// # Thread safety
///
/// The source's `@note` says the class is not thread-safe and that concurrent
/// readers need one instance each. Here every decoding method takes `&mut self`,
/// so that is a compile-time property rather than a documented hazard.
///
/// # See also
///
/// The header's `@see` list is `ImzMLFile`, `ImzMLMeta`, `ImzMLSpectrumIndex`
/// and `OnDiscMSExperiment`. Two of them exist here:
/// [`ImzMLMeta`] and
/// [`ImzMLSpectrumIndex`].
/// `FORMAT/ImzMLFile.h`, the in-memory loader and writer front end, is not
/// ported; [`ImzMLHandler`] is the
/// layer this façade actually uses. `KERNEL/OnDiscMSExperiment.h`, the indexed
/// mzML analogue named as the design precedent, is not ported either; the
/// closest thing the crate has is
/// [`indexed_mzml_handler`](crate::format::indexed_mzml_handler), which reads
/// one record at its index offset.
///
/// # Examples
///
/// The workflow the header's `@code` block sketches, on the upstream continuous
/// fixture: nine pixels on a 3 x 3 grid sharing one m/z array.
///
/// ```
/// use openms::kernel::on_disc_imzml_experiment::OnDiscImzMLExperiment;
///
/// let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
///     .join("tests/data/ImzMLFile_1_Example_Continuous.imzML");
///
/// let mut experiment = OnDiscImzMLExperiment::new();
/// experiment.open(path)?; // parses the XML index only, no peak I/O
///
/// assert!(experiment.is_open());
/// assert_eq!(experiment.len(), 9);
/// assert_eq!((experiment.grid_width(), experiment.grid_height()), (3, 3));
/// assert_eq!(experiment.meta().imaging_mode.unwrap().as_str(), "continuous");
///
/// // The grid is built during open(): zero-based, 2-D, one pixel per spectrum.
/// assert_eq!(experiment.geometry().number_of_pixels(), 9);
/// assert_eq!(experiment.geometry().spectrum_index(0, 0), Some(0));
/// assert!(experiment.geometry_report().is_clean());
///
/// // One seek and read per array, on demand.
/// let by_index = experiment.spectrum(0)?;
/// let by_coord = experiment.spectrum_at_pixel(1, 1)?; // imzML coords are 1-based
/// assert_eq!(by_index.peaks.len(), 8399);
/// assert_eq!(by_coord.peaks.len(), by_index.peaks.len());
///
/// // The index entry itself costs no .ibd read.
/// assert_eq!(experiment.index(0)?.mz_offset, 16);
///
/// let image = experiment.extract_ion_image(100.0, 1000.0)?;
/// assert_eq!((image.width(), image.height()), (3, 3));
/// assert!(image.has_pixel(0, 0));
///
/// experiment.close();
/// assert!(!experiment.is_open());
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Debug, Default)]
pub struct OnDiscImzMLExperiment {
    state: Backing,
    geometry: ImagingGeometry,
    report: GeometryReport,
    uuid: Option<UuidStatus>,
    imzml_path: PathBuf,
    ibd_path: PathBuf,
}

impl OnDiscImzMLExperiment {
    /// An experiment with nothing opened: no spectra, an empty grid and default
    /// metadata.
    ///
    /// Source default constructor, which allocates an empty `Impl`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Open an imzML dataset, deriving the `.ibd` path from `imzml_path`.
    ///
    /// The derivation is the source's `ImzMLFile::inferIbdPath_`: a
    /// case-insensitive `.imzML` suffix is replaced by `.ibd`, and any other
    /// name simply gains `.ibd`. The source's `open()` spells this out a second
    /// time, in a branch that cannot be reached — see the note on
    /// [`open_with_limits`](Self::open_with_limits).
    ///
    /// # Errors
    ///
    /// As [`open_with_limits`](Self::open_with_limits).
    pub fn open(&mut self, imzml_path: impl AsRef<Path>) -> Result<()> {
        let imzml_path = imzml_path.as_ref();
        let ibd_path = infer_ibd_path(imzml_path);
        self.open_with_limits(imzml_path, ibd_path, ImzMLReadLimits::default())
    }

    /// Open an imzML dataset with an explicit `.ibd` path.
    ///
    /// Source `open(imzml_path, ibd_path)` with a non-empty second argument.
    /// The override is threaded through the index load *and* the UUID check, so
    /// both target the file that will actually be read rather than an inferred
    /// sibling that may be missing or stale; the upstream suite guards exactly
    /// that case.
    ///
    /// # Errors
    ///
    /// As [`open_with_limits`](Self::open_with_limits).
    pub fn open_with_ibd(
        &mut self,
        imzml_path: impl AsRef<Path>,
        ibd_path: impl AsRef<Path>,
    ) -> Result<()> {
        self.open_with_limits(imzml_path, ibd_path, ImzMLReadLimits::default())
    }

    /// Open an imzML dataset with an explicit `.ibd` path and explicit resource
    /// ceilings.
    ///
    /// Parses the XML index, verifies the `.ibd` UUID header, builds the
    /// two-dimensional pixel grid and leaves the `.ibd` open. No peak data is
    /// read: arrays are decoded lazily, one seek and read per
    /// [`spectrum`](Self::spectrum) call. The grid is built eagerly here — an
    /// in-memory pass over the already-parsed index, negligible next to the XML
    /// parse — so a problem with the coordinate grid surfaces at `open()`
    /// rather than on the first pixel query, which is the source's stated
    /// reason for doing it here.
    ///
    /// A UUID mismatch between the `.ibd` header and the XML's `IMS:1000080`,
    /// an out-of-grid pixel, a coordinate below 1 and a duplicate pixel are all
    /// non-fatal: the dataset still loads, and of duplicates only the first
    /// spectrum per pixel enters the grid. The source writes each to its warning
    /// log; here they are retained, in [`uuid_status`](Self::uuid_status) and
    /// [`geometry_report`](Self::geometry_report).
    ///
    /// # Arguments
    ///
    /// * `imzml_path` — path to the `.imzML` file.
    /// * `ibd_path` — path to the companion `.ibd`. The source defaults this to
    ///   the empty string and then derives it; [`open`](Self::open) is that
    ///   overload.
    /// * `limits` — ceilings the layer below preflights every `.ibd` range
    ///   against.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when either file cannot be opened or read, where the source
    /// throws `Exception::FileNotFound`; [`Error::Parse`] when the `.imzML` is
    /// malformed, when an IMS value does not parse, or when a spectrum declares
    /// no pixel coordinate, where the source throws `Exception::ParseError`;
    /// [`Error::InvalidValue`] when the file exceeds one of `limits` or declares
    /// a grid above [`MAX_IMAGE_PIXELS`].
    ///
    /// On failure this experiment is left exactly as it was. The source instead
    /// replaces its `Impl` on the first line of `open()`, so a failed open there
    /// discards whatever dataset was previously loaded; building into a
    /// temporary and committing at the end is this crate's atomicity rule.
    ///
    /// Note also that the source re-derives the `.ibd` path after the index
    /// load, from `meta_.ibd_file_path`, and falls back to its own copy of
    /// `inferIbdPath_` when that is empty. It never is: `loadImpl_` assigns
    /// `meta_.ibd_file_path` unconditionally from the override or from
    /// `inferIbdPath_`. That fallback is unreachable duplicated logic and has no
    /// counterpart here.
    pub fn open_with_limits(
        &mut self,
        imzml_path: impl AsRef<Path>,
        ibd_path: impl AsRef<Path>,
        limits: ImzMLReadLimits,
    ) -> Result<()> {
        let mut handler = ImzMLHandler::open_with_limits(imzml_path, ibd_path, limits)?;
        let uuid = handler.uuid_status()?;
        let (geometry, report) = build_imaging_geometry(handler.index(), handler.meta())?;
        self.imzml_path = handler.imzml_path().to_path_buf();
        self.ibd_path = handler.ibd_path().to_path_buf();
        self.geometry = geometry;
        self.report = report;
        self.uuid = Some(uuid);
        self.state = Backing::Open(Box::new(handler));
        Ok(())
    }

    /// Close the companion `.ibd` and release the on-disc resources.
    ///
    /// Afterwards [`is_open`](Self::is_open) is `false` and every decoding call
    /// fails until [`open`](Self::open) runs again. The parsed index and
    /// metadata survive, so [`len`](Self::len), [`index`](Self::index) and
    /// [`meta`](Self::meta) keep answering; the grid is cleared, as the source's
    /// `close()` clears `geometry_` and retains `index_` and `meta_`.
    ///
    /// Source `close()` is `noexcept`; this cannot fail either.
    ///
    /// Retaining the index means copying it out of the reader that owned it, so
    /// this call is O(number of spectra) where the source's is O(1). The copy
    /// buys the source's post-close contract without keeping the file handle
    /// alive.
    pub fn close(&mut self) {
        self.state = match std::mem::take(&mut self.state) {
            Backing::Open(handler) => Backing::Closed(Box::new(handler.parsed().clone())),
            closed => closed,
        };
        self.geometry.clear();
    }

    /// Whether the companion `.ibd` is open, as source `isOpen()`.
    pub fn is_open(&self) -> bool {
        matches!(self.state, Backing::Open(_))
    }

    /// The number of indexed spectra.
    ///
    /// Source `getNrSpectra()`, and `size()`, which is defined inline as a
    /// synonym for it; one method covers both.
    pub fn len(&self) -> usize {
        self.parsed().len()
    }

    /// Whether the dataset indexed no spectra.
    ///
    /// Native convenience; the source expresses it as `getNrSpectra() == 0`.
    pub fn is_empty(&self) -> bool {
        self.parsed().is_empty()
    }

    /// The `.imzML` path this experiment parsed.
    ///
    /// Native accessor; the source keeps only the `.ibd` path in its `Impl` and
    /// exposes neither.
    pub fn imzml_path(&self) -> &Path {
        &self.imzml_path
    }

    /// The `.ibd` path this experiment opened.
    ///
    /// Source `Impl::ibd_path_`, reachable there only as
    /// `getImzMLMeta().ibd_file_path`, which holds the same value.
    pub fn ibd_path(&self) -> &Path {
        &self.ibd_path
    }

    /// The index entry for spectrum `index`: pixel coordinates and `.ibd` byte
    /// offsets, with no `.ibd` read.
    ///
    /// # Arguments
    ///
    /// * `index` — zero-based spectrum position.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `index` is not below [`len`](Self::len),
    /// where the source throws `Exception::IndexOverflow`.
    pub fn index(&self, index: usize) -> Result<&ImzMLSpectrumIndex> {
        self.parsed()
            .get(index)
            .ok_or_else(|| out_of_range(index, self.len()))
    }

    /// Imaging metadata parsed during `open()`, with no `.ibd` read, as source
    /// `getImzMLMeta()`.
    pub fn meta(&self) -> &ImzMLMeta {
        &self.parsed().meta
    }

    /// The declared image width, as source `gridWidth()`, a shorthand for
    /// `getImzMLMeta().max_count_x`.
    ///
    /// This is `IMS:1000042` as the reader recorded it, raised to the largest
    /// observed x, and it is not necessarily
    /// [`geometry().width()`](ImagingGeometry::width) — the grid derives its own
    /// extent when the file declares none.
    pub fn grid_width(&self) -> u32 {
        self.meta().max_count_x
    }

    /// The declared image height, as source `gridHeight()`, a shorthand for
    /// `getImzMLMeta().max_count_y`.
    pub fn grid_height(&self) -> u32 {
        self.meta().max_count_y
    }

    /// The two-dimensional imaging grid: pixel extent plus the
    /// `(x, y) -> spectrum index` map.
    ///
    /// Source `getGeometry() const`, which returns the shared
    /// `MSImagingGeometry` so that on-disc and in-memory access expose pixel
    /// coordinates the same way. Coordinates are **zero-based**: imzML's
    /// one-based values are normalised, and only the `z == 1` plane is
    /// represented. Built during `open()` from the parsed index with no `.ibd`
    /// read, so this is an O(1) borrow and any coordinate problem already
    /// surfaced in [`geometry_report`](Self::geometry_report).
    pub fn geometry(&self) -> &ImagingGeometry {
        &self.geometry
    }

    /// Mutable access to the imaging grid.
    ///
    /// Source `getGeometry()`, the non-const overload, whose documented purpose
    /// is annotating a loaded dataset with user-defined regions:
    ///
    /// ```
    /// use openms::kernel::on_disc_imzml_experiment::{
    ///     ImagingRegion, OnDiscImzMLExperiment,
    /// };
    ///
    /// let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
    ///     .join("tests/data/ImzMLFile_1_Example_Continuous.imzML");
    /// let mut experiment = OnDiscImzMLExperiment::new();
    /// experiment.open(path)?;
    ///
    /// let roi = ImagingRegion::rectangle(1, "roi", 0, 0, 1, 1)?;
    /// experiment.geometry_mut().add_region(roi)?;
    /// assert_eq!(experiment.geometry().region_pixels(1)?.len(), 4);
    /// # Ok::<(), openms::Error>(())
    /// ```
    ///
    /// Rust cannot overload on the receiver's mutability, so the two source
    /// overloads become two names.
    pub fn geometry_mut(&mut self) -> &mut ImagingGeometry {
        &mut self.geometry
    }

    /// What the grid construction had to leave out, retained from `open()`.
    ///
    /// The source emits these as warnings and keeps nothing; see
    /// [`GeometryReport`]. Unchanged by [`close`](Self::close), which clears the
    /// grid itself.
    pub fn geometry_report(&self) -> &GeometryReport {
        &self.report
    }

    /// The verdict of comparing the `.ibd` UUID header with the XML's
    /// `IMS:1000080`, or `None` before the first successful `open()`.
    ///
    /// Source `ImzMLFile::verifyIbdUuid_`, which every read path runs after the
    /// parse and which warns for each non-matching outcome and loads anyway.
    /// `open()` here rejects none of them either; the verdict is retained so
    /// that a caller can.
    pub fn uuid_status(&self) -> Option<&UuidStatus> {
        self.uuid.as_ref()
    }

    /// Decode spectrum `index` from the `.ibd`.
    ///
    /// The returned spectrum carries the peaks decoded from the `.ibd`, the
    /// pixel coordinates as the `imzml:x`, `imzml:y` and `imzml:z` meta values,
    /// and one float data array per decoded auxiliary external array — the same
    /// contract as the in-memory loader, so a viewer needs no second code path.
    /// mzML scan metadata — retention time, MS level, instrument — is not
    /// loaded in on-disc mode, in the source either.
    ///
    /// A peak array that declares no `IMS:1000101` is decoded from its own
    /// inline base64 instead, as source `ImzMLInterceptConsumer` fills the
    /// non-external side from the peaks `MzMLHandler` decoded. Conformant
    /// imzML 1.1.0 always stores both arrays externally, so that only arises
    /// for a non-conformant dataset;
    /// [`extract_ion_image`](Self::extract_ion_image) resolves each array by
    /// the same rule, so both paths answer alike for such a file.
    ///
    /// Peaks are sorted by m/z before returning, after the auxiliary arrays are
    /// attached, so those arrays stay aligned with peak order. That ordering is
    /// the source's: `decodeSpectrum` attaches first and calls
    /// `sortByPosition()` last, with a comment saying why.
    ///
    /// # Arguments
    ///
    /// * `index` — zero-based spectrum position.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `index` is not below [`len`](Self::len),
    /// where the source throws `Exception::IndexOverflow`; [`Error::Io`] with
    /// [`std::io::ErrorKind::NotFound`] when the `.ibd` is not open, where the
    /// source throws `Exception::FileNotFound`; [`Error::Unsupported`] when the
    /// m/z, intensity or any auxiliary array is compressed — anything other
    /// than `MS:1000576` — and [`Error::Parse`] when the decoded m/z and
    /// intensity arrays have different lengths, both of which the source
    /// reports as `Exception::ParseError`.
    ///
    /// The source additionally raises `ParseError` with "failed to decode
    /// spectrum peaks from .ibd" when a declared non-zero array decodes to
    /// nothing. That is unreachable here: the layer below preflights each range
    /// against the `.ibd` length and reports a short read as an error, so a
    /// non-zero count cannot silently produce an empty array.
    ///
    /// One further source step has no counterpart: `decodeSpectrum` marks a
    /// spectrum carrying per-peak ion mobility as `IMPeakType::IM_PROFILE`.
    /// [`MSSpectrum`] has no such field — see `docs/SPECTRUM_MOBILITY_SUPPORT.md`
    /// — so the flag is not set and a caller cannot distinguish profile from
    /// centroided mobility here.
    pub fn spectrum(&mut self, index: usize) -> Result<MSSpectrum> {
        self.decoded_spectrum(index).map(|decoded| decoded.spectrum)
    }

    /// Decode spectrum `index` and also report what the decode left out.
    ///
    /// Exactly [`spectrum`](Self::spectrum), returning the reader's
    /// [`DecodedSpectrum`] rather than only its spectrum, so that the auxiliary
    /// arrays the source warns about and drops — unnamed, zero-length, of the
    /// wrong length, of an unsupported type, or inline while the peaks are
    /// external — are visible instead of logged.
    ///
    /// # Errors
    ///
    /// As [`spectrum`](Self::spectrum).
    pub fn decoded_spectrum(&mut self, index: usize) -> Result<DecodedSpectrum> {
        if index >= self.len() {
            return Err(out_of_range(index, self.len()));
        }
        let ibd_path = self.ibd_path.clone();
        let handler = self.handler_mut(&ibd_path)?;
        let mut decoded = handler.spectrum(index)?;
        decoded.spectrum.sort_by_position()?;
        Ok(decoded)
    }

    /// Decode the spectrum at the imzML pixel coordinate `(x, y, z)`.
    ///
    /// Coordinates are imzML-native **one-based**, matching
    /// [`ImzMLSpectrumIndex::x`](crate::format::imzml_handler::ImzMLSpectrumIndex::x).
    /// The lookup goes through the shared grid, which is zero-based and
    /// two-dimensional, so `(x, y)` is converted and only the `z == 1` plane is
    /// addressable — consistent with the in-memory loader, which likewise builds
    /// a two-dimensional grid. The grid is built once during `open()`, so this
    /// is an O(1) lookup plus one peak decode.
    ///
    /// # Arguments
    ///
    /// * `x` — pixel column, one-based.
    /// * `y` — pixel row, one-based.
    /// * `z` — depth slice, one-based; only `1` is supported. The source
    ///   defaults it to 1 and Rust has no default arguments, so
    ///   [`spectrum_at_pixel`](Self::spectrum_at_pixel) is that overload.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] with [`std::io::ErrorKind::NotFound`] when the `.ibd` is
    /// not open — checked before the coordinate, as in the source;
    /// [`Error::InvalidValue`] when no spectrum exists at `(x, y, z)`, where the
    /// source throws `Exception::ElementNotFound`, which is also the answer for
    /// any `z` other than 1 and for a coordinate below 1; otherwise as
    /// [`spectrum`](Self::spectrum).
    ///
    /// A spectrum on another plane, or at a duplicated or out-of-grid pixel, is
    /// unreachable by coordinate and reachable by index. The lower-level
    /// [`ImzMLHandler::spectrum_at_coord`](crate::format::imzml_handler::ImzMLHandler::spectrum_at_coord)
    /// scans the index instead of the grid and so addresses every `z` in the
    /// file; use it when the grid's two-dimensionality is the problem.
    pub fn spectrum_at_coord(&mut self, x: u32, y: u32, z: u32) -> Result<MSSpectrum> {
        if !self.is_open() {
            return Err(not_open(&self.ibd_path));
        }
        let position = if z == 1 && x >= 1 && y >= 1 {
            self.geometry.spectrum_index(x - 1, y - 1)
        } else {
            None
        };
        let position = position.ok_or_else(|| {
            Error::InvalidValue(format!(
                "no imzML spectrum is mapped to the imaging pixel ({x},{y},{z})"
            ))
        })?;
        self.spectrum(position)
    }

    /// Decode the spectrum at the imzML pixel coordinate `(x, y)` of the first
    /// plane.
    ///
    /// The source's `getSpectrumAtCoord(x, y)` with `z` left at its default of
    /// 1.
    ///
    /// # Errors
    ///
    /// As [`spectrum_at_coord`](Self::spectrum_at_coord).
    pub fn spectrum_at_pixel(&mut self, x: u32, y: u32) -> Result<MSSpectrum> {
        self.spectrum_at_coord(x, y, 1)
    }

    /// Extract an ion image over the whole dataset by summing peak intensities
    /// inside `[mz - dm, mz + dm]`, with `dm = mz * tolerance_ppm * 1e-6`.
    ///
    /// The on-disc counterpart of the in-memory `MSImagingExperiment::extractIonImage`:
    /// it walks the shared grid and decodes each pixel's peaks from the `.ibd`
    /// on demand, so the full dataset never has to be held in memory, which is
    /// what makes single-mass ion images of large datasets viewable. The window
    /// is inclusive at both ends, because the source sums from `MZBegin(lo)` —
    /// the first peak at or above `lo` — to `MZEnd(hi)`, the first peak strictly
    /// above `hi`.
    ///
    /// Only m/z and intensity are decoded. Auxiliary arrays do not contribute to
    /// the sum, so a compressed auxiliary array is *not* reported here, unlike
    /// [`spectrum`](Self::spectrum) — the source makes the same distinction and
    /// documents it.
    ///
    /// Pixels absent from the grid stay invalid in the returned image. A pixel
    /// with a spectrum but no peak in the window is marked valid with intensity
    /// 0. The image's m/z range is set to `[mz - dm, mz + dm]`.
    ///
    /// # Arguments
    ///
    /// * `mz` — m/z centre of the extraction window; must be finite and `>= 0`.
    /// * `tolerance_ppm` — half-window width in ppm; must be finite and `>= 0`.
    ///
    /// # Returns
    ///
    /// An image with the same dimensions as
    /// [`geometry()`](Self::geometry) — not necessarily
    /// [`grid_width`](Self::grid_width) x [`grid_height`](Self::grid_height),
    /// which are the declared extent rather than the grid's.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `mz` or `tolerance_ppm` is negative or not
    /// finite, when a pixel references a spectrum index at or above
    /// [`len`](Self::len), or when the grid is larger than
    /// [`MAX_IMAGE_PIXELS`]; the source throws `Exception::InvalidValue` for the
    /// first two. [`Error::Io`] with [`std::io::ErrorKind::NotFound`] when the
    /// `.ibd` is not open, where the source throws `Exception::FileNotFound`.
    /// Otherwise as [`spectrum`](Self::spectrum), minus the auxiliary-array
    /// conditions.
    ///
    /// A grid carrying pixels while both its dimensions are zero yields
    /// `Error::InvalidValue` on the first pixel, because a 0 x 0 image has
    /// nowhere to put it; the source reaches the same dead end as
    /// `Exception::IndexOverflow` from `IonImage::linearIndex_`. Only a
    /// hand-built grid can be in that state — [`build_imaging_geometry`] always
    /// derives an extent that covers the pixels it inserted.
    ///
    /// The source walks the pixels with a `#pragma`-free serial loop, and so
    /// does this; neither parallelises the sweep.
    pub fn extract_ion_image(&mut self, mz: f64, tolerance_ppm: f64) -> Result<IonImage> {
        self.extract(mz, tolerance_ppm, None)
    }

    /// Extract an ion image from the acquired pixels of one region only.
    ///
    /// The same summation and window semantics as
    /// [`extract_ion_image`](Self::extract_ion_image), limited to the pixels
    /// belonging to `region_id` and decoded lazily from the `.ibd`. Pixels
    /// outside the region stay invalid, so the image keeps the grid's full
    /// dimensions.
    ///
    /// # Arguments
    ///
    /// * `mz` — m/z centre of the extraction window; finite and `>= 0`.
    /// * `tolerance_ppm` — half-window width in ppm; finite and `>= 0`.
    /// * `region_id` — the region to extract, as added through
    ///   [`geometry_mut`](Self::geometry_mut).
    ///
    /// # Errors
    ///
    /// As [`extract_ion_image`](Self::extract_ion_image), plus
    /// [`Error::InvalidValue`] when `region_id` is unknown, where the source
    /// throws `Exception::ElementNotFound` out of `getRegionPixels`.
    pub fn extract_ion_image_in_region(
        &mut self,
        mz: f64,
        tolerance_ppm: f64,
        region_id: usize,
    ) -> Result<IonImage> {
        let positions = self.geometry.region_pixels(region_id)?;
        self.extract(mz, tolerance_ppm, Some(&positions))
    }

    /// The shared extraction kernel: source `Internal::extractIonImage`, whose
    /// only caller-specific part is how a spectrum is obtained.
    ///
    /// `positions` is `None` for the whole grid, matching the source's
    /// `std::iota` over every pixel without materialising the index vector.
    fn extract(
        &mut self,
        mz: f64,
        tolerance_ppm: f64,
        positions: Option<&[usize]>,
    ) -> Result<IonImage> {
        if !self.is_open() {
            return Err(not_open(&self.ibd_path));
        }
        if !mz.is_finite() || !tolerance_ppm.is_finite() || mz < 0.0 || tolerance_ppm < 0.0 {
            return Err(Error::InvalidValue(format!(
                "ion image mz and tolerance_ppm must be finite and non-negative: mz={mz}, tolerance_ppm={tolerance_ppm}"
            )));
        }
        let dm = mz * tolerance_ppm * 1e-6;
        let lo = mz - dm;
        let hi = mz + dm;

        let mut image = IonImage::new(self.geometry.width(), self.geometry.height())?;
        image.set_mz_range(RangeBase::from_min_max(lo, hi)?);

        let spectra = self.len();
        let count = positions.map_or_else(|| self.geometry.number_of_pixels(), <[usize]>::len);
        for step in 0..count {
            let position = positions.map_or(step, |positions| positions[step]);
            let pixel = *self.geometry.pixels().get(position).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "imaging pixel {position} is not below the {} pixels of the grid",
                    self.geometry.number_of_pixels()
                ))
            })?;
            if pixel.spectrum_index >= spectra {
                return Err(Error::InvalidValue(format!(
                    "imaging pixel ({},{}) references the missing spectrum {}",
                    pixel.x, pixel.y, pixel.spectrum_index
                )));
            }
            let spectrum = self.peaks_only(pixel.spectrum_index)?;
            let begin = spectrum.mz_begin(lo)?;
            let end = spectrum.mz_end(hi)?;
            let sum: f64 = spectrum.peaks[begin..end]
                .iter()
                .map(|peak| f64::from(peak.intensity))
                .sum();
            image.set_intensity(pixel.x, pixel.y, sum)?;
        }
        Ok(image)
    }

    /// Decode only the m/z and intensity arrays of spectrum `index`, sorted.
    ///
    /// Source `Impl::decodePeaks`, which exists because ion-image extraction
    /// sums inside an m/z window and therefore needs neither the auxiliary
    /// arrays nor the pixel meta values that `decodeSpectrum` has to produce.
    /// Skipping the auxiliary reads is also why a compressed auxiliary array
    /// does not fail an extraction.
    ///
    /// The two peak arrays come from
    /// [`ImzMLHandler::mz_array`](crate::format::imzml_handler::ImzMLHandler::mz_array)
    /// and its intensity counterpart, which resolve an array to the `.ibd` or
    /// to its inline base64 by the same `IMS:1000101` rule
    /// [`spectrum`](Self::spectrum) uses. That shared rule is what makes an
    /// ion image and a decoded spectrum report the same peaks for the same
    /// pixel: an earlier revision had those accessors read the `.ibd`
    /// unconditionally, so a non-conformant file with an inline peak array
    /// that still carried an `IMS:1000102` offset produced an ion image from
    /// bytes `spectrum` never returned.
    fn peaks_only(&mut self, index: usize) -> Result<MSSpectrum> {
        let (x, y, z) = {
            let entry = self.index(index)?;
            (entry.x, entry.y, entry.z)
        };
        let ibd_path = self.ibd_path.clone();
        let handler = self.handler_mut(&ibd_path)?;
        let mz = handler.mz_array(index)?;
        let intensity = handler.intensity_array(index)?;
        if mz.len() != intensity.len() {
            return Err(Error::Parse {
                line: 0,
                message: format!(
                    "m/z and intensity array length mismatch at pixel ({x},{y},{z}): mz={} intensity={}",
                    mz.len(),
                    intensity.len()
                ),
            });
        }
        let mut spectrum = MSSpectrum {
            peaks: mz
                .iter()
                .zip(&intensity)
                .map(|(&mz, &intensity)| crate::kernel::Peak1D::new(mz, intensity))
                .collect(),
            ..MSSpectrum::default()
        };
        spectrum.sort_by_position()?;
        Ok(spectrum)
    }

    /// The parsed index, whether or not the `.ibd` is open.
    fn parsed(&self) -> &ImzMLIndex {
        match &self.state {
            Backing::Open(handler) => handler.parsed(),
            Backing::Closed(index) => index,
        }
    }

    /// The open reader, or the source's `FileNotFound` equivalent.
    fn handler_mut(&mut self, ibd_path: &Path) -> Result<&mut ImzMLHandler> {
        match &mut self.state {
            Backing::Open(handler) => Ok(handler),
            Backing::Closed(_) => Err(not_open(ibd_path)),
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Row-major grid key, as source `MSImagingGeometry::packKey_`.
fn pack_key(x: u32, y: u32) -> u64 {
    (u64::from(y) << 32) | u64::from(x)
}

/// `width * height` as a `usize`, refused above [`MAX_IMAGE_PIXELS`].
fn grid_cells(width: u32, height: u32) -> Result<usize> {
    let cells = u64::from(width) * u64::from(height);
    usize::try_from(cells)
        .ok()
        .filter(|&cells| cells <= MAX_IMAGE_PIXELS)
        .ok_or_else(|| {
            Error::InvalidValue(format!(
                "an image of {width}x{height} needs {cells} pixels, above the ceiling of {MAX_IMAGE_PIXELS}"
            ))
        })
}

/// The inclusive far corner `origin + extent - 1`, refused when it leaves `u32`.
fn far_corner(origin: u32, extent: u32, what: &str) -> Result<u32> {
    origin
        .checked_add(extent - 1)
        .ok_or_else(|| Error::InvalidRange(format!("imaging region {what} extent overflows u32")))
}

/// Reject a non-finite scalar, as the crate does on every numeric entry point.
fn finite(value: f64, what: &str) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::InvalidValue(format!("{what} must be finite")))
    }
}

/// Append to a problem list only while it is below [`MAX_LISTED_PROBLEMS`].
fn push_capped(list: &mut Vec<usize>, position: usize) {
    if list.len() < MAX_LISTED_PROBLEMS {
        list.push(position);
    }
}

/// Source `Exception::IndexOverflow` on a spectrum index.
fn out_of_range(index: usize, len: usize) -> Error {
    Error::InvalidValue(format!(
        "imzML spectrum {index} is not below the {len} spectra of the dataset"
    ))
}

/// Source `Exception::ElementNotFound` on a region identifier.
fn unknown_region(id: usize) -> Error {
    Error::InvalidValue(format!("no imaging region carries the identifier {id}"))
}

/// Source `Exception::FileNotFound` on the companion `.ibd`.
fn not_open(ibd_path: &Path) -> Error {
    Error::Io(std::io::Error::new(
        std::io::ErrorKind::NotFound,
        format!(
            "the companion .ibd '{}' is not open; call open() first",
            ibd_path.display()
        ),
    ))
}
