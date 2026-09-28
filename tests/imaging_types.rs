// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Section-by-section ports of the three imaging value-type class tests at
//! core SDK `bc9cc12`, with unchanged literals:
//!
//! * `src/tests/class_tests/openms/source/MSImagingRegion_test.cpp`
//! * `src/tests/class_tests/openms/source/MSImagingGeometry_test.cpp`
//! * `src/tests/class_tests/openms/source/IonImage_test.cpp`
//!
//! The Rust types live in `openms::kernel::on_disc_imzml_experiment`
//! (`ImagingRegion`, `ImagingGeometry`, `IonImage`); the API mapping is in
//! `docs/ON_DISC_IMZML_SUPPORT.md`. Where the source throws, the Rust call
//! returns `Err` (or `None` for the two lookups that return an `Option`), as
//! each section notes. `NOT_TESTABLE` sections and destructor sections have no
//! body to port and are listed where they occur.

use openms::kernel::on_disc_imzml_experiment::{
    ImagingGeometry, ImagingRegion, IonImage, RegionShape,
};
use openms::kernel::ranges::RangeBase;

fn similar(actual: f64, expected: f64) {
    // TEST_REAL_SIMILAR's default relative tolerance is 1e-5; every value
    // here is stored and read back unchanged, so exact equality holds.
    assert_eq!(actual, expected);
}

// ===========================================================================
// MSImagingRegion_test.cpp
// ===========================================================================

// START_SECTION((MSImagingRegion())) / ((~MSImagingRegion())): construction and
// Drop; covered by the defaults section.

// START_SECTION((MSImagingRegion() defaults))
#[test]
fn region_defaults() {
    let r = ImagingRegion::default();
    assert_eq!(r.id(), 0);
    assert_eq!(r.bbox_width(), 1);
    assert_eq!(r.bbox_height(), 1);
    assert_eq!(r.min_x(), 0);
    assert_eq!(r.min_y(), 0);
    assert_eq!(r.max_x(), 0);
    assert_eq!(r.max_y(), 0);
    assert_eq!(r.area(), 1);
    assert_eq!(r.mask().len(), 0);
    assert_eq!(r.shape(), RegionShape::Rectangle);
}

// START_SECTION((MSImagingRegion() rectangle))
#[test]
fn region_rectangle() {
    let rectangle = ImagingRegion::rectangle(7, "rect", 2, 3, 5, 6).unwrap();
    assert_eq!(rectangle.id(), 7);
    assert_eq!(rectangle.name(), "rect");
    assert_eq!(rectangle.bbox_width(), 4);
    assert_eq!(rectangle.bbox_height(), 4);
    assert_eq!(rectangle.min_x(), 2);
    assert_eq!(rectangle.min_y(), 3);
    assert_eq!(rectangle.max_x(), 5);
    assert_eq!(rectangle.max_y(), 6);
    assert_eq!(rectangle.area(), 16);
    assert_eq!(rectangle.mask().len(), 0);
    assert_eq!(rectangle.shape(), RegionShape::Rectangle);
    assert!(rectangle.contains(2, 3));
    assert!(!rectangle.contains(1, 4));

    let single_cell = ImagingRegion::rectangle(8, "single_cell", 1, 1, 1, 1).unwrap();
    assert_eq!(single_cell.area(), 1);
    assert!(single_cell.contains(1, 1));
    assert!(!single_cell.contains(1, 2));
    assert!(!single_cell.contains(2, 1));

    // TEST_EXCEPTION(Exception::InvalidValue, ...)
    assert!(ImagingRegion::rectangle(9, "exc", 10, 1, 1, 1).is_err());
}

// START_SECTION((MSImagingRegion() fromMask))
#[test]
fn region_from_mask() {
    let mask = vec![true, false, false, true];
    let from_mask = ImagingRegion::from_mask(0, "mask", 10, 20, 2, 2, mask.clone()).unwrap();
    assert_eq!(from_mask.mask(), mask.as_slice());
    assert_eq!(from_mask.shape(), RegionShape::Mask);
    assert!(!from_mask.contains(11, 20));
    assert!(from_mask.contains(11, 21));
    assert!(!from_mask.contains(1, 2));

    // Four TEST_EXCEPTION(Exception::InvalidValue, ...) cases.
    assert!(ImagingRegion::from_mask(0, "mask", 10, 20, 0, 2, mask.clone()).is_err());
    assert!(ImagingRegion::from_mask(0, "mask", 10, 20, 2, 0, mask.clone()).is_err());
    assert!(ImagingRegion::from_mask(0, "mask", 10, 20, 2, 2, vec![false; 4]).is_err());
    assert!(
        ImagingRegion::from_mask(
            0,
            "mask",
            10,
            20,
            2,
            2,
            vec![false, false, false, false, true]
        )
        .is_err()
    );
}

// START_SECTION((bool intersects(const MSImagingRegion& other) const))
#[test]
fn region_intersects() {
    let a = ImagingRegion::rectangle(1, "a", 1, 1, 1, 1).unwrap();
    let b = ImagingRegion::rectangle(2, "b", 1, 1, 5, 5).unwrap();
    assert!(a.intersects(&b));
    assert!(b.intersects(&a));
    let c = ImagingRegion::rectangle(3, "c", 10, 10, 11, 11).unwrap();
    assert!(!c.intersects(&a));
    let m1 =
        ImagingRegion::from_mask(4, "mask", 1, 1, 2, 2, vec![true, false, false, true]).unwrap();
    assert!(m1.intersects(&a));
    let m2 =
        ImagingRegion::from_mask(5, "mask2", 1, 1, 2, 2, vec![false, true, true, false]).unwrap();
    assert!(!m2.intersects(&m1));
    let m3 =
        ImagingRegion::from_mask(6, "mask3", 1, 1, 2, 2, vec![true, true, true, false]).unwrap();
    assert!(m3.intersects(&m2));
}

// START_SECTION(area() & contains()): NOT_TESTABLE, covered above.

// ===========================================================================
// MSImagingGeometry_test.cpp
// ===========================================================================

// START_SECTION((MSImagingGeometry()))
#[test]
fn geometry_default_constructor() {
    let g = ImagingGeometry::new();
    assert_eq!(g.width(), 0);
    assert_eq!(g.height(), 0);
    assert_eq!(g.number_of_pixels(), 0);
    similar(g.pixel_size_x(), 1.0);
    similar(g.pixel_size_y(), 1.0);
    assert_eq!(g.pixel_size_unit(), "micrometer");
}

// START_SECTION((~MSImagingGeometry())): Drop.

// START_SECTION((void setDimensions(UInt width, UInt height)))
#[test]
fn geometry_set_dimensions() {
    let mut g = ImagingGeometry::new();
    g.set_dimensions(10, 8).unwrap();
    assert_eq!(g.width(), 10);
    assert_eq!(g.height(), 8);
}

// START_SECTION((void setPixelSize(double x, double y, const std::string& unit)))
#[test]
fn geometry_set_pixel_size() {
    let mut g = ImagingGeometry::new();
    g.set_pixel_size(25.0, 25.0, "micrometer").unwrap();
    similar(g.pixel_size_x(), 25.0);
    similar(g.pixel_size_y(), 25.0);
    assert_eq!(g.pixel_size_unit(), "micrometer");
    g.set_pixel_size(50.0, 75.0, "nanometer").unwrap();
    similar(g.pixel_size_x(), 50.0);
    similar(g.pixel_size_y(), 75.0);
    assert_eq!(g.pixel_size_unit(), "nanometer");
}

// START_SECTION((void addPixel(UInt x, UInt y, Size spectrum_index)))
#[test]
fn geometry_add_pixel() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(0, 0, 42).unwrap();
    g.add_pixel(1, 0, 7).unwrap();
    assert_eq!(g.number_of_pixels(), 2);
    assert!(g.has_pixel(0, 0));
    assert!(g.has_pixel(1, 0));
    assert!(!g.has_pixel(2, 0));
    assert_eq!(g.spectrum_index(0, 0), Some(42));
    assert_eq!(g.spectrum_index(1, 0), Some(7));
    // duplicate insertion: TEST_EXCEPTION(Exception::InvalidValue, ...)
    assert!(g.add_pixel(0, 0, 99).is_err());

    let mut g2 = ImagingGeometry::new();
    g2.set_dimensions(3, 2).unwrap();
    g2.add_pixel(2, 1, 5).unwrap();
    assert!(g2.add_pixel(3, 0, 6).is_err());
    assert!(g2.add_pixel(0, 2, 7).is_err());
}

// START_SECTION((bool hasPixel(UInt x, UInt y) const))
#[test]
fn geometry_has_pixel() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(5, 9, 100).unwrap();
    assert!(g.has_pixel(5, 9));
    assert!(!g.has_pixel(0, 0));
}

// START_SECTION((Size getSpectrumIndex(UInt x, UInt y) const))
#[test]
fn geometry_spectrum_index() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(2, 3, 17).unwrap();
    assert_eq!(g.spectrum_index(2, 3), Some(17));
    // TEST_EXCEPTION(Exception::ElementNotFound, ...): None here.
    assert_eq!(g.spectrum_index(0, 0), None);
}

// START_SECTION((const std::vector<Pixel>& getPixels() const))
#[test]
fn geometry_pixels_in_insertion_order() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(2, 0, 10).unwrap();
    g.add_pixel(0, 1, 20).unwrap();
    g.add_pixel(1, 1, 30).unwrap();
    let pixels = g.pixels();
    assert_eq!(pixels.len(), 3);
    assert_eq!(
        (pixels[0].x, pixels[0].y, pixels[0].spectrum_index),
        (2, 0, 10)
    );
    assert_eq!(
        (pixels[1].x, pixels[1].y, pixels[1].spectrum_index),
        (0, 1, 20)
    );
    assert_eq!(
        (pixels[2].x, pixels[2].y, pixels[2].spectrum_index),
        (1, 1, 30)
    );
}

// START_SECTION((Size getNumberOfPixels() const)): NOT_TESTABLE.

// START_SECTION((void clear()))
#[test]
fn geometry_clear() {
    let mut g = ImagingGeometry::new();
    g.set_dimensions(4, 4).unwrap();
    g.set_pixel_size(10.0, 10.0, "micrometer").unwrap();
    g.add_pixel(0, 0, 1).unwrap();
    g.add_pixel(1, 1, 2).unwrap();
    g.clear();
    assert_eq!(g.width(), 0);
    assert_eq!(g.height(), 0);
    assert_eq!(g.number_of_pixels(), 0);
    assert!(!g.has_pixel(0, 0));
    similar(g.pixel_size_x(), 1.0);
    assert_eq!(g.pixel_size_unit(), "micrometer");
}

// getWidth/getHeight/getPixelSizeX/getPixelSizeY/getPixelSizeUnit: NOT_TESTABLE.

// START_SECTION((addRegion()))
#[test]
fn geometry_add_region() {
    let mut g = ImagingGeometry::new();
    g.set_dimensions(6, 6).unwrap();
    g.add_pixel(0, 0, 1).unwrap();
    g.add_pixel(1, 0, 2).unwrap();
    g.add_pixel(5, 5, 3).unwrap();
    g.add_region(ImagingRegion::rectangle(1, "name", 0, 0, 1, 0).unwrap())
        .unwrap();
    assert_eq!(g.number_of_regions(), 1);
    let spec_idx = g.region_spectrum_indices(1).unwrap();
    assert_eq!(spec_idx, [1, 2]);
    assert!(
        g.add_region(ImagingRegion::rectangle(2, "collision", 0, 0, 10, 10).unwrap())
            .is_err()
    );
    assert_eq!(g.number_of_regions(), 1);
    assert!(
        g.add_region(ImagingRegion::rectangle(1, "dup", 5, 5, 6, 6).unwrap())
            .is_err()
    );
    assert!(
        g.add_region(
            ImagingRegion::rectangle(ImagingGeometry::NO_REGION, "bad id", 5, 5, 6, 6).unwrap()
        )
        .is_err()
    );
    g.add_region(ImagingRegion::rectangle(2, "far away", 100, 100, 100, 100).unwrap())
        .unwrap();
    assert_eq!(g.number_of_regions(), 2);
    assert_eq!(g.region_pixels(2).unwrap().len(), 0);
}

// START_SECTION((regionOf())): NO_REGION is `None` here.
#[test]
fn geometry_region_of() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(0, 0, 0).unwrap();
    g.add_pixel(5, 5, 1).unwrap();
    g.add_region(ImagingRegion::rectangle(1, "name", 0, 0, 1, 0).unwrap())
        .unwrap();
    assert_eq!(g.region_of(0, 0), Some(1));
    assert_eq!(g.region_of(5, 5), None);
    assert_eq!(g.region_of(10, 10), None);
}

// START_SECTION((getRegionPixels()))
#[test]
fn geometry_region_pixels() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(0, 0, 1).unwrap();
    g.add_pixel(1, 0, 2).unwrap();
    g.add_region(ImagingRegion::rectangle(1, "name", 0, 0, 1, 0).unwrap())
        .unwrap();
    let pixel_idx = g.region_pixels(1).unwrap();
    assert_eq!(pixel_idx[0], 0);
    assert_eq!(pixel_idx[1], 1);
}

// START_SECTION((getRegionSpectrumIndices()))
#[test]
fn geometry_region_spectrum_indices() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(0, 0, 1).unwrap();
    g.add_pixel(1, 0, 2).unwrap();
    g.add_region(ImagingRegion::rectangle(1, "name", 0, 0, 1, 0).unwrap())
        .unwrap();
    let spec_idx = g.region_spectrum_indices(1).unwrap();
    assert_eq!(spec_idx.len(), 2);
    assert_eq!(spec_idx[0], 1);
    assert_eq!(spec_idx[1], 2);
}

// START_SECTION((region overlay : Mask region))
#[test]
fn geometry_mask_region_overlay() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(10, 20, 1).unwrap();
    g.add_pixel(11, 20, 2).unwrap();
    g.add_pixel(11, 21, 3).unwrap();
    g.add_region(
        ImagingRegion::from_mask(1, "m", 10, 20, 2, 2, vec![true, false, false, true]).unwrap(),
    )
    .unwrap();
    assert_eq!(g.region_of(10, 20), Some(1));
    assert_eq!(g.region_of(11, 20), None);
    assert_eq!(g.region_of(11, 21), Some(1));
    assert_eq!(g.region_spectrum_indices(1).unwrap().len(), 2);
}

// START_SECTION((removeRegion()))
#[test]
fn geometry_remove_region() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(0, 0, 1).unwrap();
    g.add_pixel(5, 5, 2).unwrap();
    g.add_region(ImagingRegion::rectangle(1, "A", 0, 0, 0, 0).unwrap())
        .unwrap();
    g.add_region(ImagingRegion::rectangle(2, "B", 5, 5, 5, 5).unwrap())
        .unwrap();
    g.remove_region(1).unwrap();
    assert_eq!(g.number_of_regions(), 1);
    assert_eq!(g.region_of(0, 0), None);
    assert_eq!(g.region_of(5, 5), Some(2));
    assert_eq!(g.number_of_pixels(), 2);
}

// START_SECTION((clearRegions()))
#[test]
fn geometry_clear_regions() {
    let mut g = ImagingGeometry::new();
    g.add_pixel(0, 0, 1).unwrap();
    g.add_region(ImagingRegion::rectangle(1, "A", 0, 0, 0, 0).unwrap())
        .unwrap();
    g.clear_regions();
    assert_eq!(g.number_of_regions(), 0);
    assert_eq!(g.region_of(0, 0), None);
    assert_eq!(g.number_of_pixels(), 1);
}

// START_SECTION((getRegion()))
#[test]
fn geometry_get_region() {
    let mut g = ImagingGeometry::new();
    g.add_region(ImagingRegion::rectangle(1, "A", 0, 0, 1, 1).unwrap())
        .unwrap();
    assert_eq!(g.region(1).unwrap().id(), 1);
    // Three TEST_EXCEPTION(Exception::ElementNotFound, ...) cases.
    assert!(g.region(99).is_err());
    assert!(g.region_pixels(99).is_err());
    assert!(g.remove_region(99).is_err());
}

// START_SECTION((membership staleness test))
#[test]
fn geometry_membership_is_live() {
    let mut g = ImagingGeometry::new();
    g.add_region(ImagingRegion::rectangle(1, "A", 0, 0, 0, 0).unwrap())
        .unwrap();
    g.add_pixel(0, 0, 1).unwrap();
    assert_eq!(g.region_of(0, 0), Some(1));
    assert_eq!(g.region_spectrum_indices(1).unwrap().len(), 1);
}

// START_SECTION((addRegion rejects geometrically overlapping footprints regardless of pixels))
#[test]
fn geometry_rejects_overlapping_footprints() {
    let mut g = ImagingGeometry::new();
    g.add_region(ImagingRegion::rectangle(1, "A", 0, 0, 5, 5).unwrap())
        .unwrap();
    assert!(
        g.add_region(ImagingRegion::rectangle(2, "B", 3, 3, 8, 8).unwrap())
            .is_err()
    );
    assert_eq!(g.number_of_regions(), 1);
    g.add_region(ImagingRegion::rectangle(3, "C", 100, 100, 110, 110).unwrap())
        .unwrap();
    assert_eq!(g.number_of_regions(), 2);
}

// ===========================================================================
// IonImage_test.cpp
// ===========================================================================

// START_SECTION((IonImage()))
#[test]
fn ion_image_default_constructor() {
    let img = IonImage::default();
    assert_eq!(img.width(), 0);
    assert_eq!(img.height(), 0);
    assert!(!img.has_pixel(0, 0));
}

// START_SECTION((~IonImage())): Drop.

// START_SECTION((IonImage(UInt width, UInt height)))
#[test]
fn ion_image_sized_constructor() {
    let img = IonImage::new(3, 2).unwrap();
    assert_eq!(img.width(), 3);
    assert_eq!(img.height(), 2);
    assert_eq!(img.data().len(), 6);
    assert_eq!(img.mask().len(), 6);
    for y in 0..2 {
        for x in 0..3 {
            assert!(!img.has_pixel(x, y));
            similar(img.intensity(x, y).unwrap(), 0.0);
        }
    }
}

// START_SECTION((void resize(UInt width, UInt height)))
#[test]
fn ion_image_resize() {
    let mut img = IonImage::default();
    img.resize(2, 2).unwrap();
    img.set_intensity(0, 0, 1.0).unwrap();
    img.resize(4, 3).unwrap();
    assert_eq!(img.width(), 4);
    assert_eq!(img.height(), 3);
    assert_eq!(img.data().len(), 12);
    assert!(!img.has_pixel(0, 0));
    similar(img.intensity(0, 0).unwrap(), 0.0);
}

// START_SECTION((void setIntensity(UInt x, UInt y, double intensity)))
#[test]
fn ion_image_set_intensity() {
    let mut img = IonImage::new(3, 2).unwrap();
    img.set_intensity(1, 0, 7.5).unwrap();
    assert!(img.has_pixel(1, 0));
    similar(img.intensity(1, 0).unwrap(), 7.5);
    assert!(!img.has_pixel(2, 1));
    similar(img.intensity(2, 1).unwrap(), 0.0);
    // index = y * W + x = 1 * 3 + 2 = 5
    img.set_intensity(2, 1, 42.0).unwrap();
    similar(img.data()[5], 42.0);
    assert!(img.mask()[5]);
}

// START_SECTION((void setMzRange(const RangeMZ& range)))
#[test]
fn ion_image_set_mz_range() {
    let mut img = IonImage::default();
    img.set_mz_range(RangeBase::from_min_max(499.9, 500.1).unwrap());
    similar(img.mz_range().min().unwrap(), 499.9);
    similar(img.mz_range().max().unwrap(), 500.1);
}

// START_SECTION((const RangeMZ& getMzRange() const)): NOT_TESTABLE.

// START_SECTION((double getIntensity(UInt x, UInt y) const))
#[test]
fn ion_image_intensity_out_of_bounds() {
    let img = IonImage::new(3, 2).unwrap();
    // TEST_EXCEPTION(Exception::IndexOverflow, ...)
    assert!(img.intensity(5, 0).is_err());
    assert!(img.intensity(0, 5).is_err());
}

// START_SECTION((bool hasPixel(UInt x, UInt y) const))
#[test]
fn ion_image_has_pixel_out_of_bounds() {
    let img = IonImage::new(3, 2).unwrap();
    assert!(!img.has_pixel(5, 0));
    assert!(!img.has_pixel(0, 5));
}

// getWidth/getHeight/getData/getMask: NOT_TESTABLE.
