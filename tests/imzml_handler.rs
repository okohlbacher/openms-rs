// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Coverage for `FORMAT/HANDLERS/ImzMLHandlerHelper.h` and
//! `FORMAT/HANDLERS/ImzMLHandler.h`.
//!
//! Neither header has a class test of its own upstream; the imzML family is
//! tested through `ImzMLFile_test.cpp` and `ImzMLFile_all_modes_test.cpp`, which
//! exercise `ImzMLFile` and `OnDiscImzMLExperiment`. The literals below that
//! come from that suite or from the two unmodified upstream fixtures are
//! transcribed source review (tier 3): the 3x3 grid, 100 µm pixels, the 300 µm
//! extents, `continuous`/`processed`, the two UUIDs, the MD5 and SHA-1 strings,
//! `float32` for both arrays, the `negative` polarity with `top down` /
//! `horizontal` / `left-right` acquisition geometry, and the first m/z of pixel
//! (1,1) in each mode — 100.0 continuous and 100.083336 processed.
//!
//! Two checks are independent of that suite. The decoded intensity sum of
//! processed pixel (1,1) is compared with the `MS:1000285` total ion current
//! that the same fixture's XML declares, which ties the `.ibd` decode to a
//! number written by the instrument's own exporter; and the declared
//! `IMS:1000091` SHA-1 of the processed `.ibd` is recomputed here, which no
//! OpenMS read path does. The synthetic documents, the resource ceilings and
//! the hostile offsets are independently derived (tier 4), because no upstream
//! fixture reaches them.

#![cfg(feature = "mzml")]

use openms::Error;
use openms::format::imzml_handler::{
    AuxSkipReason, ChecksumStatus, IBD_UUID_BYTES, ImagingMode, ImzMLBinaryIO, ImzMLDataType,
    ImzMLHandler, ImzMLIndex, ImzMLReadLimits, UuidStatus, infer_ibd_path, read_index,
    read_index_with_limits, uuid_bytes, write_float32_array, write_float64_array,
    write_mz_as_float32, write_mz_as_float64,
};
use openms::system::file::TempDir;
use std::io::BufReader;
use std::path::PathBuf;

const CONTINUOUS: &str = "ImzMLFile_1_Example_Continuous.imzML";
const CONTINUOUS_IBD: &str = "ImzMLFile_1_Example_Continuous.ibd";
const PROCESSED: &str = "ImzMLFile_2_Example_Processed.imzML";
const PROCESSED_IBD: &str = "ImzMLFile_2_Example_Processed.ibd";
const CONTINUOUS_IBD_LEN: u64 = 335_976;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn open(name: &str) -> ImzMLHandler {
    ImzMLHandler::open(data(name)).expect("upstream imzML fixture opens")
}

fn index_of(name: &str) -> ImzMLIndex {
    let file = BufReader::new(std::fs::File::open(data(name)).unwrap());
    read_index(file).expect("upstream imzML fixture indexes")
}

/// A minimal well-formed `.imzML` around `body`, in the shape of the upstream
/// fixtures: the IMS terms this parser reads and nothing else.
fn document(body: &str) -> String {
    format!(
        concat!(
            "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
            "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">",
            "{}",
            "</mzML>\n"
        ),
        body
    )
}

/// One synthetic `<spectrum>` with pixel coordinates and an external m/z and
/// intensity pair addressing `.ibd` bytes. The defaults address the upstream
/// continuous `.ibd`.
struct Spec<'a> {
    id: &'a str,
    x: u32,
    y: u32,
    z: Option<u32>,
    mz_offset: u64,
    int_offset: u64,
    mz_length: u64,
    int_length: u64,
    type_accession: &'a str,
    extra: &'a str,
}

impl Default for Spec<'_> {
    fn default() -> Self {
        Self {
            id: "s=1",
            x: 1,
            y: 1,
            z: None,
            mz_offset: 16,
            int_offset: 33_612,
            mz_length: 1,
            int_length: 1,
            type_accession: "MS:1000521",
            extra: "",
        }
    }
}

impl Spec<'_> {
    fn with_length(length: u64) -> Self {
        Self {
            mz_length: length,
            int_length: length,
            ..Self::default()
        }
    }

    fn xml(&self) -> String {
        let z = match self.z {
            Some(z) => {
                format!("<cvParam accession=\"IMS:1000052\" name=\"position z\" value=\"{z}\"/>")
            }
            None => String::new(),
        };
        format!(
            concat!(
                "<spectrum id=\"{id}\" index=\"0\" defaultArrayLength=\"0\">",
                "<scanList count=\"1\"><scan>",
                "<cvParam accession=\"IMS:1000050\" name=\"position x\" value=\"{x}\"/>",
                "<cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"{y}\"/>",
                "{z}",
                "</scan></scanList>",
                "<binaryDataArrayList count=\"2\">",
                "{mz}",
                "{int}",
                "{extra}",
                "</binaryDataArrayList>",
                "</spectrum>"
            ),
            id = self.id,
            x = self.x,
            y = self.y,
            z = z,
            mz = peak_array(
                "MS:1000514",
                self.type_accession,
                self.mz_offset,
                self.mz_length
            ),
            int = peak_array(
                "MS:1000515",
                self.type_accession,
                self.int_offset,
                self.int_length
            ),
            extra = self.extra
        )
    }

    fn document(&self) -> String {
        document(&format!(
            "<run><spectrumList count=\"1\">{}</spectrumList></run>",
            self.xml()
        ))
    }
}

fn peak_array(role: &str, type_accession: &str, offset: u64, length: u64) -> String {
    format!(
        concat!(
            "<binaryDataArray encodedLength=\"0\">",
            "<cvParam accession=\"{role}\" name=\"peak array\"/>",
            "<cvParam accession=\"{type_accession}\" name=\"binary type\"/>",
            "<cvParam accession=\"MS:1000576\" name=\"no compression\"/>",
            "<cvParam accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>",
            "<cvParam accession=\"IMS:1000102\" name=\"external offset\" value=\"{offset}\"/>",
            "<cvParam accession=\"IMS:1000103\" name=\"external array length\" value=\"{length}\"/>",
            "<binary/>",
            "</binaryDataArray>"
        ),
        role = role,
        type_accession = type_accession,
        offset = offset,
        length = length
    )
}

/// Write `xml` next to the unmodified continuous `.ibd` and open both.
fn open_against_continuous_ibd(xml: &str) -> (TempDir, ImzMLHandler) {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("synthetic.imzML");
    std::fs::write(&path, xml).unwrap();
    let handler = ImzMLHandler::open_with_ibd(&path, data(CONTINUOUS_IBD))
        .expect("synthetic document opens against the upstream .ibd");
    (dir, handler)
}

// ---------------------------------------------------------------------------
// Dataset metadata, both upstream fixtures
// ---------------------------------------------------------------------------

#[test]
fn continuous_fixture_metadata_matches_the_upstream_literals() {
    let handler = open(CONTINUOUS);
    let meta = handler.meta();
    assert_eq!(meta.imaging_mode, Some(ImagingMode::Continuous));
    assert_eq!(meta.imaging_mode.unwrap().as_str(), "continuous");
    assert_eq!(
        (meta.max_count_x, meta.max_count_y, meta.max_count_z),
        (3, 3, 1)
    );
    assert_eq!(meta.pixel_size_x, 100.0);
    assert_eq!(meta.pixel_size_y, 100.0);
    assert_eq!(meta.max_dim_x, 300.0);
    assert_eq!(meta.max_dim_y, 300.0);
    assert_eq!(meta.uuid, "12345678-1234-1234-1234-123456789012");
    assert_eq!(meta.ibd_md5, "4b5dd9fa84fafc955cfdd301f9ed55d7");
    // The upstream suite asserts that this fixture declares no SHA-1.
    assert!(meta.ibd_sha1.is_empty());
    assert_eq!(meta.mz_data_type, ImzMLDataType::Float32);
    assert_eq!(meta.int_data_type, ImzMLDataType::Float32);
    // The continuous fixture declares no polarity and no acquisition geometry.
    assert!(meta.polarity.is_empty());
    assert!(meta.scan_pattern.is_empty());
    assert!(meta.scan_direction.is_empty());
    assert!(meta.line_scan_direction.is_empty());
    assert_eq!(meta.ibd_file_path, data(CONTINUOUS_IBD));
}

#[test]
fn processed_fixture_metadata_matches_the_upstream_literals() {
    let handler = open(PROCESSED);
    let meta = handler.meta();
    assert_eq!(meta.imaging_mode, Some(ImagingMode::Processed));
    assert_eq!(meta.imaging_mode.unwrap().as_str(), "processed");
    assert_eq!(
        (meta.max_count_x, meta.max_count_y, meta.max_count_z),
        (3, 3, 1)
    );
    assert_eq!(meta.pixel_size_x, 100.0);
    assert_eq!(meta.pixel_size_y, 100.0);
    assert_eq!(meta.max_dim_x, 300.0);
    assert_eq!(meta.max_dim_y, 300.0);
    assert_eq!(meta.uuid, "9d501bdc53444916b7e97e795b02c856");
    assert_eq!(meta.ibd_sha1, "7e8fdb93053915d3edb51b70aa0619ac209964df");
    assert!(meta.ibd_md5.is_empty());
    assert_eq!(meta.scan_pattern, "top down");
    assert_eq!(meta.scan_direction, "horizontal");
    assert_eq!(meta.line_scan_direction, "left-right");
    // MS:1000129 reaches the dataset block from a referenceableParamGroup
    // applied at spectrum level, where neither the scan nor the array branch
    // claims it.
    assert_eq!(meta.polarity, "negative");
    assert_eq!(meta.mz_data_type, ImzMLDataType::Float32);
    assert_eq!(meta.int_data_type, ImzMLDataType::Float32);
}

/// The processed fixture declares `encoding="ISO-8859-1"` and carries a
/// non-UTF-8 byte in an `MS:1000590` contact affiliation. Skipping the value of
/// every param no imzML rule can act on is what lets it index.
#[test]
fn a_non_utf8_declared_encoding_does_not_stop_the_index() {
    let raw = std::fs::read(data(PROCESSED)).unwrap();
    assert!(std::str::from_utf8(&raw).is_err());
    assert_eq!(index_of(PROCESSED).len(), 9);
}

/// A non-UTF-8 byte inside a value this parser must read is an error, not a
/// silent replacement: the source relies on Xerces to transcode from the
/// declared encoding and this port implements no transcoder.
#[test]
fn a_non_utf8_value_on_a_read_param_is_rejected() {
    let mut xml = document(&format!(
        "<fileDescription><fileContent>\
         <cvParam accession=\"IMS:1000080\" name=\"uuid\" value=\"PLACEHOLDER\"/>\
         </fileContent></fileDescription><run><spectrumList count=\"1\">{}</spectrumList></run>",
        Spec::default().xml()
    ))
    .into_bytes();
    let at = xml
        .windows(11)
        .position(|window| window == b"PLACEHOLDER")
        .unwrap();
    xml[at] = 0xfc;
    let error = read_index(BufReader::new(xml.as_slice())).unwrap_err();
    assert!(matches!(error, Error::Parse { .. }), "{error:?}");
}

// ---------------------------------------------------------------------------
// The index itself
// ---------------------------------------------------------------------------

/// Continuous mode stores one m/z array for the whole image, so every entry
/// names the same `IMS:1000102`, while intensity offsets advance by the stored
/// array length.
#[test]
fn continuous_mode_shares_one_mz_array_across_the_image() {
    let handler = open(CONTINUOUS);
    assert_eq!(handler.len(), 9);
    assert_eq!(handler.index()[0].mz_offset, 16);
    assert!(handler.index().iter().all(|entry| entry.mz_offset == 16));
    assert!(handler.index().iter().all(|entry| entry.mz_length == 8399));
    assert!(
        handler
            .index()
            .iter()
            .all(|entry| entry.mz_encoded_bytes == 33_596)
    );
    let intensity: Vec<u64> = handler.index().iter().map(|e| e.int_offset).collect();
    assert_eq!(intensity[0], 33_612);
    for pair in intensity.windows(2) {
        assert_eq!(pair[1] - pair[0], 33_596);
    }
}

/// Processed mode gives every pixel its own m/z array, so no two entries share
/// an offset.
#[test]
fn processed_mode_gives_every_pixel_its_own_mz_array() {
    let handler = open(PROCESSED);
    assert_eq!(handler.len(), 9);
    let mut offsets: Vec<u64> = handler.index().iter().map(|e| e.mz_offset).collect();
    let count = offsets.len();
    offsets.sort_unstable();
    offsets.dedup();
    assert_eq!(offsets.len(), count);
    assert_eq!(offsets[0], 16);
    // m/z and intensity alternate, each 33,596 stored bytes long.
    for entry in handler.index() {
        assert_eq!(entry.int_offset, entry.mz_offset + 33_596);
        assert_eq!(entry.mz_length, 8399);
        assert_eq!(entry.int_length, 8399);
    }
}

/// Both fixtures raster a 3x3 grid row by row, in 1-based imzML coordinates.
#[test]
fn index_entries_carry_document_order_and_pixel_coordinates() {
    for name in [CONTINUOUS, PROCESSED] {
        let handler = open(name);
        for (position, entry) in handler.index().iter().enumerate() {
            assert_eq!(entry.index as usize, position, "{name}");
            assert_eq!(entry.x, (position as u32 % 3) + 1, "{name}");
            assert_eq!(entry.y, (position as u32 / 3) + 1, "{name}");
            // The processed fixture declares no IMS:1000052 at all; z defaults
            // to 1, as the source's per-spectrum reset does.
            assert_eq!(entry.z, 1, "{name}");
            assert!(entry.mz_external && entry.int_external, "{name}");
            assert!(!entry.mz_compressed && !entry.int_compressed, "{name}");
            assert!(entry.aux.is_empty(), "{name}");
            assert_eq!(entry.unnamed_aux, 0, "{name}");
            assert!(entry.inline_aux_names.is_empty(), "{name}");
        }
    }
}

#[test]
fn native_ids_come_from_the_spectrum_id_attribute() {
    let continuous = open(CONTINUOUS);
    assert_eq!(continuous.index()[0].native_id, "spectrum=1");
    assert_eq!(continuous.index()[8].native_id, "spectrum=9");
    let processed = open(PROCESSED);
    assert_eq!(processed.index()[0].native_id, "Scan=1");
    assert_eq!(processed.index()[8].native_id, "Scan=9");
}

/// `IMS:1000042`/`IMS:1000043` are a declaration; the source raises the
/// bounding box to the largest coordinate it actually sees, and `max_count_z`
/// has no CV term at all.
#[test]
fn observed_coordinates_raise_the_declared_bounding_box() {
    let body = format!(
        "<scanSettingsList count=\"1\"><scanSettings id=\"s1\">\
         <cvParam accession=\"IMS:1000042\" name=\"max count of pixels x\" value=\"1\"/>\
         <cvParam accession=\"IMS:1000043\" name=\"max count of pixels y\" value=\"1\"/>\
         </scanSettings></scanSettingsList><run><spectrumList count=\"2\">{}{}</spectrumList></run>",
        Spec::default().xml(),
        // A second pixel beyond the declared grid, with an explicit z slice.
        Spec {
            id: "s=2",
            x: 7,
            y: 4,
            z: Some(5),
            ..Spec::default()
        }
        .xml()
    );
    let index = read_index(BufReader::new(document(&body).as_bytes())).unwrap();
    assert_eq!(index.meta.max_count_x, 7);
    assert_eq!(index.meta.max_count_y, 4);
    assert_eq!(index.meta.max_count_z, 5);
    assert_eq!(index.spectra[1].z, 5);
}

/// A referenceable parameter group is replayed in the context of the reference,
/// which is how both fixtures give their arrays a type and an offset base.
#[test]
fn referenceable_parameter_groups_are_replayed_at_the_reference() {
    let body = concat!(
        "<referenceableParamGroupList count=\"1\">",
        "<referenceableParamGroup id=\"mzArray\">",
        "<cvParam accession=\"MS:1000514\" name=\"m/z array\"/>",
        "<cvParam accession=\"MS:1000523\" name=\"64-bit float\"/>",
        "<cvParam accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>",
        "</referenceableParamGroup>",
        "</referenceableParamGroupList>",
        "<run><spectrumList count=\"1\">",
        "<spectrum id=\"s=1\" index=\"0\" defaultArrayLength=\"0\">",
        "<scanList count=\"1\"><scan>",
        "<cvParam accession=\"IMS:1000050\" name=\"position x\" value=\"2\"/>",
        "<cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"3\"/>",
        "</scan></scanList>",
        "<binaryDataArrayList count=\"1\">",
        "<binaryDataArray encodedLength=\"0\">",
        "<referenceableParamGroupRef ref=\"mzArray\"/>",
        "<cvParam accession=\"IMS:1000102\" name=\"external offset\" value=\"24\"/>",
        "<cvParam accession=\"IMS:1000103\" name=\"external array length\" value=\"3\"/>",
        "<binary/>",
        "</binaryDataArray>",
        // An unknown reference is ignored, as the source's failed map lookup is.
        "<binaryDataArray encodedLength=\"0\">",
        "<referenceableParamGroupRef ref=\"neverDeclared\"/>",
        "<cvParam accession=\"MS:1000515\" name=\"intensity array\"/>",
        "<cvParam accession=\"MS:1000521\" name=\"32-bit float\"/>",
        "<cvParam accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>",
        "<cvParam accession=\"IMS:1000102\" name=\"external offset\" value=\"40\"/>",
        "<cvParam accession=\"IMS:1000103\" name=\"external array length\" value=\"3\"/>",
        "<binary/>",
        "</binaryDataArray>",
        "</binaryDataArrayList></spectrum></spectrumList></run>"
    );
    let index = read_index(BufReader::new(document(body).as_bytes())).unwrap();
    let entry = &index.spectra[0];
    assert_eq!((entry.x, entry.y, entry.z), (2, 3, 1));
    assert_eq!(entry.mz_type, ImzMLDataType::Float64);
    assert_eq!(entry.mz_offset, 24);
    assert!(entry.mz_external);
    assert_eq!(entry.int_type, ImzMLDataType::Float32);
    assert!(entry.int_external);
}

/// Source `onEndElement("spectrum")` raises `Exception::ParseError` when a
/// spectrum declares neither `IMS:1000050` nor `IMS:1000051`, which is the
/// upstream `load rejects spectrum missing pixel coordinates` section.
#[test]
fn a_spectrum_without_pixel_coordinates_is_rejected() {
    let body = concat!(
        "<run><spectrumList count=\"1\">",
        "<spectrum id=\"s=1\" index=\"0\" defaultArrayLength=\"0\">",
        "<scanList count=\"1\"><scan/></scanList>",
        "<binaryDataArrayList count=\"0\"/>",
        "</spectrum></spectrumList></run>"
    );
    let error = read_index(BufReader::new(document(body).as_bytes())).unwrap_err();
    match error {
        Error::Parse { message, .. } => assert!(message.contains("pixel coordinate"), "{message}"),
        other => panic!("{other:?}"),
    }

    // y alone is not enough either.
    let body = body.replace(
        "<scan/>",
        "<scan><cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"1\"/></scan>",
    );
    assert!(matches!(
        read_index(BufReader::new(document(&body).as_bytes())),
        Err(Error::Parse { .. })
    ));
}

/// A coordinate param is only a coordinate inside a `<scan>` inside a
/// `<spectrum>`, which is the source's first dispatch block.
#[test]
fn a_position_param_outside_a_scan_is_not_a_coordinate() {
    let body = concat!(
        "<run><spectrumList count=\"1\">",
        "<spectrum id=\"s=1\" index=\"0\" defaultArrayLength=\"0\">",
        "<cvParam accession=\"IMS:1000050\" name=\"position x\" value=\"4\"/>",
        "<cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"4\"/>",
        "<binaryDataArrayList count=\"0\"/>",
        "</spectrum></spectrumList></run>"
    );
    assert!(matches!(
        read_index(BufReader::new(document(body).as_bytes())),
        Err(Error::Parse { .. })
    ));
}

#[test]
fn malformed_ims_values_are_rejected_with_the_source_reason() {
    for (value, reason) in [
        ("", "empty value"),
        ("-1", "negative value not allowed"),
        ("1.5", "not a valid unsigned integer"),
        ("12abc", "not a valid unsigned integer"),
        ("4294967296", "out of range for uint32"),
    ] {
        let body = format!(
            "<run><spectrumList count=\"1\">\
             <spectrum id=\"s=1\" index=\"0\" defaultArrayLength=\"0\">\
             <scanList count=\"1\"><scan>\
             <cvParam accession=\"IMS:1000050\" name=\"position x\" value=\"{value}\"/>\
             <cvParam accession=\"IMS:1000051\" name=\"position y\" value=\"1\"/>\
             </scan></scanList><binaryDataArrayList count=\"0\"/>\
             </spectrum></spectrumList></run>"
        );
        match read_index(BufReader::new(document(&body).as_bytes())).unwrap_err() {
            Error::Parse { message, .. } => {
                assert!(message.contains(reason), "{value:?} gave {message}");
                assert!(message.contains("IMS:1000050"), "{message}");
            }
            other => panic!("{value:?} gave {other:?}"),
        }
    }
}

/// The source's `std::stod` accepts `inf` and `nan`; a non-finite pixel size
/// cannot be used by any consumer, so this port refuses it.
#[test]
fn a_non_finite_pixel_size_is_rejected() {
    for value in ["nan", "inf", "-inf", "", "100.0.0"] {
        let body = format!(
            "<scanSettingsList count=\"1\"><scanSettings id=\"s1\">\
             <cvParam accession=\"IMS:1000046\" name=\"pixel size x\" value=\"{value}\"/>\
             </scanSettings></scanSettingsList>"
        );
        assert!(
            matches!(
                read_index(BufReader::new(document(&body).as_bytes())),
                Err(Error::Parse { .. })
            ),
            "{value:?}"
        );
    }
    let body = "<scanSettingsList count=\"1\"><scanSettings id=\"s1\">\
                <cvParam accession=\"IMS:1000046\" name=\"pixel size x\" value=\"-25.5\"/>\
                </scanSettings></scanSettingsList>";
    let index = read_index(BufReader::new(document(body).as_bytes())).unwrap();
    assert_eq!(index.meta.pixel_size_x, -25.5);
}

// ---------------------------------------------------------------------------
// Decoding one pixel out of the .ibd
// ---------------------------------------------------------------------------

/// Upstream asserts `exp[0][0].getMZ()` is 100.0 for the continuous fixture.
#[test]
fn continuous_pixel_one_one_decodes_the_shared_axis() {
    let mut handler = open(CONTINUOUS);
    let decoded = handler.spectrum(0).unwrap();
    assert!(!decoded.inline_peaks);
    assert!(decoded.skipped_aux.is_empty());
    let peaks = &decoded.spectrum.peaks;
    assert_eq!(peaks.len(), 8399);
    assert!((peaks[0].mz - 100.0).abs() < 1e-6, "{}", peaks[0].mz);
    assert!((peaks[8398].mz - 800.0).abs() < 1e-6, "{}", peaks[8398].mz);
    // float32 storage widened to f64, exactly as the source widens it: every
    // value survives a round trip through f32.
    assert!(
        peaks
            .iter()
            .all(|peak| f64::from(peak.mz as f32) == peak.mz)
    );
    assert!((peaks[1].mz - 100.083_35).abs() < 1e-5, "{}", peaks[1].mz);
    assert!(
        (peaks[0].intensity - 0.003_479_33).abs() < 1e-8,
        "{}",
        peaks[0].intensity
    );
    assert_eq!(decoded.spectrum.native_id, "spectrum=1");
}

/// Upstream asserts `exp[0][0].getMZ()` is 100.083336 for the processed
/// fixture, and the fixture's own `MS:1000285` declares the total ion current
/// of that pixel; summing the decoded intensities must reproduce it.
#[test]
fn processed_pixel_one_one_decodes_its_own_axis() {
    let mut handler = open(PROCESSED);
    let decoded = handler.spectrum(0).unwrap();
    let peaks = &decoded.spectrum.peaks;
    assert_eq!(peaks.len(), 8399);
    assert!((peaks[0].mz - 100.083_336).abs() < 1e-5, "{}", peaks[0].mz);
    let total: f64 = peaks.iter().map(|peak| f64::from(peak.intensity)).sum();
    let declared = 121.850_390_398_684_7_f64;
    assert!((total - declared).abs() < 1e-9, "{total} vs {declared}");
}

#[test]
fn decoded_spectra_carry_their_pixel_coordinates_as_meta_values() {
    use openms::metadata::MetaValueData;
    let mut handler = open(CONTINUOUS);
    let decoded = handler.spectrum(4).unwrap();
    for (key, expected) in [("imzml:x", 2), ("imzml:y", 2), ("imzml:z", 1)] {
        let value = decoded.spectrum.metadata.get(key).expect(key);
        assert_eq!(*value.data(), MetaValueData::Integer(expected), "{key}");
    }
}

/// Every pixel of the continuous fixture reads the same m/z axis, and the
/// per-array accessors agree with the assembled spectrum.
#[test]
fn array_accessors_agree_with_the_assembled_spectrum() {
    let mut handler = open(CONTINUOUS);
    let first = handler.mz_array(0).unwrap();
    let last = handler.mz_array(8).unwrap();
    assert_eq!(first, last);
    let intensity = handler.intensity_array(0).unwrap();
    let decoded = handler.spectrum(0).unwrap();
    assert_eq!(decoded.spectrum.peaks.len(), first.len());
    for (peak, (&mz, &value)) in decoded
        .spectrum
        .peaks
        .iter()
        .zip(first.iter().zip(&intensity))
    {
        assert_eq!(peak.mz, mz);
        assert_eq!(peak.intensity, value);
    }
}

// ---------------------------------------------------------------------------
// Coordinate lookup
// ---------------------------------------------------------------------------

#[test]
fn every_pixel_of_the_grid_is_addressable_by_coordinate() {
    let mut handler = open(PROCESSED);
    for y in 1..=3 {
        for x in 1..=3 {
            let position = handler.index_at_coord(x, y, 1).expect("pixel is present");
            assert_eq!(handler.index()[position].x, x);
            assert_eq!(handler.index()[position].y, y);
            let by_coord = handler.spectrum_at_coord(x, y, 1).unwrap();
            let by_index = handler.spectrum(position).unwrap();
            assert_eq!(by_coord.spectrum, by_index.spectrum);
        }
    }
    assert_eq!(handler.index_at_coord(4, 1, 1), None);
    assert_eq!(handler.index_at_coord(1, 1, 2), None);
    let error = handler.spectrum_at_coord(9, 9, 1).unwrap_err();
    assert!(matches!(error, Error::InvalidValue(_)), "{error:?}");
}

/// The source's geometry builder gives a duplicated pixel to the first
/// spectrum in document order, which the upstream suite asserts.
#[test]
fn the_first_spectrum_wins_a_duplicated_pixel() {
    let body = format!(
        "<run><spectrumList count=\"2\">{}{}</spectrumList></run>",
        Spec::default().xml(),
        Spec {
            id: "s=2",
            ..Spec::default()
        }
        .xml()
    );
    let index = read_index(BufReader::new(document(&body).as_bytes())).unwrap();
    assert_eq!(index.len(), 2);
    assert_eq!(index.index_at_coord(1, 1, 1), Some(0));
}

#[test]
fn an_out_of_range_spectrum_index_is_an_error() {
    let mut handler = open(CONTINUOUS);
    assert!(handler.entry(9).is_err());
    assert!(matches!(handler.spectrum(9), Err(Error::InvalidValue(_))));
    assert!(matches!(handler.mz_array(100), Err(Error::InvalidValue(_))));
    assert!(matches!(
        handler.intensity_array(100),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        handler.aux_array(0, 0),
        Err(Error::InvalidValue(_))
    ));
}

// ---------------------------------------------------------------------------
// UUID and checksums
// ---------------------------------------------------------------------------

#[test]
fn both_fixtures_ibd_headers_match_their_declared_uuid() {
    for name in [CONTINUOUS, PROCESSED] {
        let mut handler = open(name);
        assert_eq!(handler.uuid_status().unwrap(), UuidStatus::Match, "{name}");
    }
}

/// The source recomputes neither declared checksum. SHA-1 is verifiable here;
/// MD5 is parsed only.
#[test]
fn the_declared_ibd_sha1_is_reproducible_and_md5_is_only_parsed() {
    let mut processed = open(PROCESSED);
    assert_eq!(processed.verify_ibd_sha1().unwrap(), ChecksumStatus::Match);
    assert_eq!(
        processed.ibd().sha1_hex().unwrap(),
        "7e8fdb93053915d3edb51b70aa0619ac209964df"
    );

    let mut continuous = open(CONTINUOUS);
    assert_eq!(
        continuous.verify_ibd_sha1().unwrap(),
        ChecksumStatus::NotDeclared
    );
    // Parsed, never verified: nothing here computes MD5.
    assert_eq!(
        continuous.meta().ibd_md5,
        "4b5dd9fa84fafc955cfdd301f9ed55d7"
    );
}

/// A declared digest that does not describe the file is reported, and the
/// comparison ignores hex case.
#[test]
fn a_declared_sha1_that_does_not_match_is_reported() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let ibd = dir.path().join("wrong.ibd");
    std::fs::write(&ibd, [0u8; 32]).unwrap();
    let one_peak = Spec {
        mz_offset: 0,
        int_offset: 4,
        ..Spec::default()
    };
    let declared = "7E8FDB93053915D3EDB51B70AA0619AC209964DF";
    let xml = dir.path().join("wrong.imzML");
    std::fs::write(
        &xml,
        document(&format!(
            "<fileDescription><fileContent>\
             <cvParam accession=\"IMS:1000091\" name=\"ibd SHA-1\" value=\"{declared}\"/>\
             </fileContent></fileDescription>\
             <run><spectrumList count=\"1\">{}</spectrumList></run>",
            one_peak.xml()
        )),
    )
    .unwrap();
    let mut handler = ImzMLHandler::open_with_ibd(&xml, &ibd).unwrap();
    let real = handler.ibd().sha1_hex().unwrap();
    match handler.verify_ibd_sha1().unwrap() {
        ChecksumStatus::Mismatch {
            found,
            declared: reported,
        } => {
            assert_eq!(found, real);
            assert_eq!(reported, declared);
        }
        other => panic!("{other:?}"),
    }

    // The same digest in the other case matches.
    std::fs::write(
        &xml,
        document(&format!(
            "<fileDescription><fileContent>\
             <cvParam accession=\"IMS:1000091\" name=\"ibd SHA-1\" value=\"{}\"/>\
             </fileContent></fileDescription>\
             <run><spectrumList count=\"1\">{}</spectrumList></run>",
            real.to_uppercase(),
            one_peak.xml()
        )),
    )
    .unwrap();
    let mut handler = ImzMLHandler::open_with_ibd(&xml, &ibd).unwrap();
    assert_eq!(handler.verify_ibd_sha1().unwrap(), ChecksumStatus::Match);
}

#[test]
fn a_checksum_over_the_configured_limit_is_refused() {
    let limits = ImzMLReadLimits {
        max_checksum_bytes: 16,
        ..ImzMLReadLimits::default()
    };
    let mut ibd = ImzMLBinaryIO::open_with_limits(data(PROCESSED_IBD), limits).unwrap();
    assert!(matches!(ibd.sha1_hex(), Err(Error::InvalidValue(_))));
    assert_eq!(ibd.limits().max_checksum_bytes, 16);
}

#[test]
fn uuid_strings_parse_with_dashes_braces_or_neither() {
    let plain = uuid_bytes("12345678123412341234123456789012").unwrap();
    assert_eq!(
        uuid_bytes("12345678-1234-1234-1234-123456789012").unwrap(),
        plain
    );
    assert_eq!(
        uuid_bytes("{12345678-1234-1234-1234-123456789012}").unwrap(),
        plain
    );
    assert_eq!(
        uuid_bytes("9D501BDC53444916B7E97E795B02C856").unwrap()[0],
        0x9d
    );
    assert_eq!(plain.len(), IBD_UUID_BYTES);
    // Too short, too long and non-hex are all "not a UUID".
    assert_eq!(uuid_bytes(""), None);
    assert_eq!(uuid_bytes("1234"), None);
    assert_eq!(uuid_bytes("123456781234123412341234567890123"), None);
    assert_eq!(uuid_bytes("1234567812341234123412345678901z"), None);
}

#[test]
fn a_mismatched_or_missing_uuid_is_reported_not_rejected() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let one_peak = Spec {
        mz_offset: 0,
        int_offset: 4,
        ..Spec::default()
    };

    // A .ibd whose header is 16 bytes of a different UUID.
    let mismatched = dir.path().join("mismatch.ibd");
    std::fs::write(&mismatched, [0xabu8; 64]).unwrap();
    let xml = dir.path().join("mismatch.imzML");
    std::fs::write(
        &xml,
        document(&format!(
            "<fileDescription><fileContent>\
             <cvParam accession=\"IMS:1000080\" name=\"uuid\" \
             value=\"12345678-1234-1234-1234-123456789012\"/>\
             </fileContent></fileDescription>\
             <run><spectrumList count=\"1\">{}</spectrumList></run>",
            one_peak.xml()
        )),
    )
    .unwrap();
    let mut handler = ImzMLHandler::open_with_ibd(&xml, &mismatched).unwrap();
    match handler.uuid_status().unwrap() {
        UuidStatus::Mismatch { found, declared } => {
            assert_eq!(found, "abababababababababababababababab");
            assert_eq!(declared, "12345678123412341234123456789012");
        }
        other => panic!("{other:?}"),
    }
    // The dataset still opens and decodes, as the source intends.
    assert_eq!(handler.spectrum(0).unwrap().spectrum.peaks.len(), 1);

    // A .ibd shorter than the header.
    let short = dir.path().join("short.ibd");
    std::fs::write(&short, [0u8; 8]).unwrap();
    let mut handler = ImzMLHandler::open_with_ibd(&xml, &short).unwrap();
    assert_eq!(handler.uuid_status().unwrap(), UuidStatus::IbdTooShort);
    assert_eq!(handler.ibd().uuid().unwrap(), None);

    // No IMS:1000080 at all, and an empty one, are both "not declared": the
    // source keeps the previous value for an empty value.
    for uuid_param in [
        "",
        "<cvParam accession=\"IMS:1000080\" name=\"uuid\" value=\"\"/>",
    ] {
        let bare = dir.path().join("bare.imzML");
        std::fs::write(
            &bare,
            document(&format!(
                "<fileDescription><fileContent>{uuid_param}</fileContent></fileDescription>\
                 <run><spectrumList count=\"1\">{}</spectrumList></run>",
                one_peak.xml()
            )),
        )
        .unwrap();
        let mut handler = ImzMLHandler::open_with_ibd(&bare, &mismatched).unwrap();
        assert_eq!(handler.uuid_status().unwrap(), UuidStatus::NotDeclared);
        assert!(handler.meta().uuid.is_empty());
    }
}

// ---------------------------------------------------------------------------
// Bounded work: every offset and length is preflighted
// ---------------------------------------------------------------------------

/// An offset and length that leave the `.ibd` are refused before any
/// allocation. The source resizes the output first and discovers the truncation
/// when `fread` comes up short.
#[test]
fn a_range_past_the_end_of_the_ibd_is_refused() {
    let (_dir, mut handler) = open_against_continuous_ibd(
        &Spec {
            mz_offset: CONTINUOUS_IBD_LEN - 8,
            mz_length: 8399,
            int_length: 8399,
            ..Spec::default()
        }
        .document(),
    );
    match handler.spectrum(0).unwrap_err() {
        Error::Parse { message, .. } => {
            assert!(message.contains("extends past"), "{message}");
            assert!(
                message.contains(&CONTINUOUS_IBD_LEN.to_string()),
                "{message}"
            );
        }
        other => panic!("{other:?}"),
    }
}

/// A count above the ceiling never reaches an allocation. 100,000,001 float32
/// values would be 400 MB in the source, which allocates them and then fails.
#[test]
fn an_element_count_above_the_ceiling_is_refused_before_allocating() {
    let (_dir, mut handler) =
        open_against_continuous_ibd(&Spec::with_length(100_000_001).document());
    match handler.spectrum(0).unwrap_err() {
        Error::InvalidValue(message) => {
            assert!(message.contains("100000001"), "{message}");
            assert!(message.contains("100000000"), "{message}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_byte_length_above_the_configured_limit_is_refused() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("input.imzML");
    std::fs::write(&path, Spec::with_length(8399).document()).unwrap();
    let limits = ImzMLReadLimits {
        max_array_bytes: 1024,
        ..ImzMLReadLimits::default()
    };
    let mut handler = ImzMLHandler::open_with_limits(&path, data(CONTINUOUS_IBD), limits).unwrap();
    match handler.spectrum(0).unwrap_err() {
        Error::InvalidValue(message) => assert!(message.contains("byte limit"), "{message}"),
        other => panic!("{other:?}"),
    }
}

/// An offset near `u64::MAX` overflows the end of the range. Checked
/// arithmetic turns that into an error instead of a wrapped, plausible range.
#[test]
fn an_offset_that_overflows_the_range_is_refused() {
    let (_dir, mut handler) = open_against_continuous_ibd(
        &Spec {
            mz_offset: u64::MAX,
            int_offset: u64::MAX,
            mz_length: 4,
            int_length: 4,
            ..Spec::default()
        }
        .document(),
    );
    match handler.spectrum(0).unwrap_err() {
        Error::InvalidValue(message) => assert!(message.contains("overflows"), "{message}"),
        other => panic!("{other:?}"),
    }
}

/// A zero-length array never touches the file, so even a nonsense offset on one
/// is not an error. Source `readMzArray` returns before its seek.
#[test]
fn a_zero_length_array_never_reads_and_never_fails() {
    let (_dir, mut handler) = open_against_continuous_ibd(
        &Spec {
            mz_offset: u64::MAX,
            int_offset: u64::MAX,
            mz_length: 0,
            int_length: 0,
            ..Spec::default()
        }
        .document(),
    );
    let decoded = handler.spectrum(0).unwrap();
    assert!(decoded.spectrum.peaks.is_empty());
    assert!(handler.mz_array(0).unwrap().is_empty());
    assert!(handler.intensity_array(0).unwrap().is_empty());
}

/// An array with no `MS:1000521`/`523`/`519`/`522` is `Unknown`, which is not
/// an indexing error and is a decoding one.
#[test]
fn an_array_without_a_binary_data_type_cannot_be_decoded() {
    // An accession that is a real CV term but not a binary data type.
    let (_dir, mut handler) = open_against_continuous_ibd(
        &Spec {
            type_accession: "MS:1000127",
            mz_length: 4,
            int_length: 4,
            ..Spec::default()
        }
        .document(),
    );
    assert_eq!(handler.index()[0].mz_type, ImzMLDataType::Unknown);
    match handler.spectrum(0).unwrap_err() {
        Error::Unsupported(message) => assert!(message.contains("data type"), "{message}"),
        other => panic!("{other:?}"),
    }
}

/// A truncated `.ibd` fails the preflight rather than the read, because the
/// range is checked against the file's real length.
#[test]
fn a_truncated_ibd_fails_the_preflight() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let ibd = dir.path().join("truncated.ibd");
    std::fs::write(&ibd, [0u8; 20]).unwrap();
    let xml = dir.path().join("truncated.imzML");
    std::fs::write(&xml, Spec::with_length(8399).document()).unwrap();
    let mut handler = ImzMLHandler::open_with_ibd(&xml, &ibd).unwrap();
    assert_eq!(handler.ibd().len(), 20);
    assert!(!handler.ibd().is_empty());
    assert!(matches!(handler.spectrum(0), Err(Error::Parse { .. })));
}

/// A pair of external arrays whose declared lengths disagree cannot become a
/// peak list; the source names the pixel in the same message.
#[test]
fn mismatched_mz_and_intensity_lengths_name_the_pixel() {
    let (_dir, mut handler) = open_against_continuous_ibd(
        &Spec {
            x: 2,
            y: 3,
            mz_length: 6,
            int_length: 4,
            ..Spec::default()
        }
        .document(),
    );
    match handler.spectrum(0).unwrap_err() {
        Error::Parse { message, .. } => {
            assert!(message.contains("(2,3,1)"), "{message}");
            assert!(message.contains("mz=6"), "{message}");
            assert!(message.contains("intensity=4"), "{message}");
        }
        other => panic!("{other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Compression is indexed and refused at decode time
// ---------------------------------------------------------------------------

/// Upstream `load and OnDisc reject zlib-compressed external m/z and intensity`
/// asserts both index flags and the decode failure.
#[test]
fn compressed_external_peak_arrays_are_indexed_and_then_refused() {
    let compressed = Spec::with_length(4).document().replace(
        "<cvParam accession=\"MS:1000576\" name=\"no compression\"/>",
        "<cvParam accession=\"MS:1000574\" name=\"zlib compression\"/>",
    );
    let (_dir, mut handler) = open_against_continuous_ibd(&compressed);
    assert!(handler.index()[0].mz_compressed);
    assert!(handler.index()[0].int_compressed);
    for error in [
        handler.spectrum(0).unwrap_err(),
        handler.mz_array(0).unwrap_err(),
        handler.intensity_array(0).unwrap_err(),
    ] {
        match error {
            Error::Unsupported(message) => {
                assert!(message.contains("MS:1000576"), "{message}");
                assert!(
                    message.contains("Re-export without compression"),
                    "{message}"
                );
            }
            other => panic!("{other:?}"),
        }
    }
}

/// `MS:1000576` clears the flag, and a numpress term sets it, because the rule
/// is "any child of `MS:1000572` other than `MS:1000576`" rather than a list.
#[test]
fn any_non_uncompressed_child_of_the_compression_term_counts() {
    for (accession, expected) in [
        ("MS:1000576", false),
        ("MS:1000574", true),
        ("MS:1002312", true), // MS-Numpress linear prediction
        ("MS:1002746", true), // MS-Numpress linear prediction with zlib
    ] {
        let xml = Spec::with_length(4).document().replace(
            "accession=\"MS:1000576\" name=\"no compression\"",
            &format!("accession=\"{accession}\" name=\"compression\""),
        );
        let index = read_index(BufReader::new(xml.as_bytes())).unwrap();
        assert_eq!(index.spectra[0].mz_compressed, expected, "{accession}");
        assert_eq!(index.spectra[0].int_compressed, expected, "{accession}");
    }
}

// ---------------------------------------------------------------------------
// Auxiliary arrays
// ---------------------------------------------------------------------------

const EXTERNAL: &str = "<cvParam accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>";

fn aux_array(accession: &str, length: u64, extra: &str) -> String {
    format!(
        concat!(
            "<binaryDataArray encodedLength=\"0\">",
            "<cvParam accession=\"{accession}\" name=\"aux\" value=\"free text\" ",
            "unitAccession=\"MS:1002814\"/>",
            "<cvParam accession=\"MS:1000521\" name=\"32-bit float\"/>",
            "<cvParam accession=\"MS:1000576\" name=\"no compression\"/>",
            "{extra}",
            "<cvParam accession=\"IMS:1000102\" name=\"external offset\" value=\"16\"/>",
            "<cvParam accession=\"IMS:1000103\" name=\"external array length\" value=\"{length}\"/>",
            "<binary/>",
            "</binaryDataArray>"
        ),
        accession = accession,
        length = length,
        extra = extra
    )
}

fn with_aux(aux: &str) -> String {
    Spec {
        mz_length: 4,
        int_length: 4,
        extra: aux,
        ..Spec::default()
    }
    .document()
}

/// An ion-mobility array is named from its own CV term, not from the XML `name`
/// attribute, and keeps the `unitAccession` the array-identity param carried.
#[test]
fn an_external_ion_mobility_array_decodes_into_a_named_float_array() {
    let (_dir, mut handler) =
        open_against_continuous_ibd(&with_aux(&aux_array("MS:1003006", 4, EXTERNAL)));
    let entry = &handler.index()[0];
    assert_eq!(entry.aux.len(), 1);
    assert_eq!(entry.aux[0].name, "mean inverse reduced ion mobility array");
    assert_eq!(entry.aux[0].accession, "MS:1003006");
    assert_eq!(entry.aux[0].unit_accession, "MS:1002814");
    assert_eq!(entry.aux[0].data_type, ImzMLDataType::Float32);
    assert_eq!(entry.aux[0].length, 4);

    let raw = handler.aux_array(0, 0).unwrap();
    let decoded = handler.spectrum(0).unwrap();
    assert!(decoded.skipped_aux.is_empty());
    assert_eq!(decoded.spectrum.float_data_arrays.len(), 1);
    let array = &decoded.spectrum.float_data_arrays[0];
    assert_eq!(array.name, "mean inverse reduced ion mobility array");
    assert_eq!(array.data, raw);
    assert!(array.metadata.contains_key("unit_accession"));
    // The array shares the m/z offset in this document, so its values are the
    // first four m/z of the shared axis narrowed to f32.
    assert_eq!(array.data[0], 100.0_f32);
}

/// `MS:1000786` names an array from the param's `value`, as the source does so
/// that a free-text array keeps the name a writer gave it.
#[test]
fn a_non_standard_array_is_named_from_its_value() {
    let (_dir, handler) =
        open_against_continuous_ibd(&with_aux(&aux_array("MS:1000786", 4, EXTERNAL)));
    assert_eq!(handler.index()[0].aux[0].name, "free text");
    assert_eq!(handler.index()[0].aux[0].accession, "MS:1000786");
}

/// Each of the source's warn-and-skip conditions drops the array and keeps the
/// spectrum, and is reported rather than only logged.
#[test]
fn skippable_auxiliary_arrays_are_reported_and_the_spectrum_survives() {
    let cases: [(&str, u64, &str, AuxSkipReason); 4] = [
        ("MS:1003006", 0, EXTERNAL, AuxSkipReason::ZeroLength),
        (
            "MS:1003006",
            9,
            EXTERNAL,
            AuxSkipReason::LengthMismatch {
                length: 9,
                peaks: 4,
            },
        ),
        // No external flag: the payload would be inline base64 while the peaks
        // come from the .ibd.
        ("MS:1003006", 4, "", AuxSkipReason::Inline),
        // A CV term that is not a binary data array at all leaves the array
        // unnamed, so the index drops it and only the count survives.
        ("MS:1000127", 4, EXTERNAL, AuxSkipReason::Unnamed),
    ];
    for (accession, length, external, reason) in cases {
        let (_dir, mut handler) =
            open_against_continuous_ibd(&with_aux(&aux_array(accession, length, external)));
        let decoded = handler.spectrum(0).unwrap();
        assert_eq!(decoded.spectrum.peaks.len(), 4, "{accession} {length}");
        assert!(
            decoded.spectrum.float_data_arrays.is_empty(),
            "{accession} {length}"
        );
        assert_eq!(decoded.skipped_aux.len(), 1, "{accession} {length}");
        assert_eq!(decoded.skipped_aux[0].reason, reason);
        // The source keeps only named arrays in the index entry; the count of
        // the nameless ones is what makes that drop visible.
        let entry = handler.entry(0).unwrap();
        if reason == AuxSkipReason::Unnamed {
            assert!(entry.aux.is_empty());
            assert_eq!(entry.unnamed_aux, 1);
        } else {
            assert_eq!(entry.unnamed_aux, 0);
        }
    }
}

/// An auxiliary array with no supported binary data type is skipped rather than
/// aborting the load, which is the source's explicit choice.
#[test]
fn an_auxiliary_array_without_a_data_type_is_skipped() {
    let aux = aux_array("MS:1003006", 4, EXTERNAL).replace(
        "<cvParam accession=\"MS:1000521\" name=\"32-bit float\"/>",
        "",
    );
    let (_dir, mut handler) = open_against_continuous_ibd(&with_aux(&aux));
    let decoded = handler.spectrum(0).unwrap();
    assert_eq!(
        decoded.skipped_aux[0].reason,
        AuxSkipReason::UnknownDataType
    );
    assert_eq!(decoded.spectrum.peaks.len(), 4);
}

/// A non-external auxiliary array is recorded by name so the drop can be
/// reported, which is the source's `inline_aux_names` list.
#[test]
fn inline_auxiliary_arrays_are_named_in_the_index() {
    let (_dir, handler) = open_against_continuous_ibd(&with_aux(&aux_array("MS:1003006", 4, "")));
    let entry = &handler.index()[0];
    assert!(entry.aux.is_empty());
    assert_eq!(
        entry.inline_aux_names,
        vec!["mean inverse reduced ion mobility array".to_string()]
    );
}

/// A compressed auxiliary array is the one condition the source raises rather
/// than skips, and its message names the `IMS:1000104` length.
#[test]
fn a_compressed_auxiliary_array_is_an_error() {
    let aux = aux_array("MS:1003006", 0, EXTERNAL)
        .replace(
            "<cvParam accession=\"MS:1000576\" name=\"no compression\"/>",
            "<cvParam accession=\"MS:1000574\" name=\"zlib compression\"/>",
        )
        .replace(
            "<binary/>",
            "<cvParam accession=\"IMS:1000104\" name=\"external encoded length\" \
             value=\"64\"/><binary/>",
        );
    let (_dir, mut handler) = open_against_continuous_ibd(&with_aux(&aux));
    assert!(handler.index()[0].aux[0].compressed);
    assert_eq!(handler.index()[0].aux[0].encoded_bytes, 64);
    assert_eq!(handler.index()[0].aux[0].length, 0);
    match handler.spectrum(0).unwrap_err() {
        Error::Unsupported(message) => {
            assert!(message.contains("mean inverse reduced"), "{message}");
            assert!(message.contains("encoded length=64"), "{message}");
        }
        other => panic!("{other:?}"),
    }
    assert!(matches!(
        handler.aux_array(0, 0),
        Err(Error::Unsupported(_))
    ));
}

// ---------------------------------------------------------------------------
// Resource ceilings on the index scan
// ---------------------------------------------------------------------------

type LimitAdjust = fn(&mut ImzMLReadLimits);

#[test]
fn the_index_scan_honours_its_ceilings() {
    let xml = document(&format!(
        "<run><spectrumList count=\"2\">{}{}</spectrumList></run>",
        Spec {
            extra: &aux_array("MS:1003006", 1, EXTERNAL),
            ..Spec::default()
        }
        .xml(),
        Spec {
            id: "s=2",
            x: 2,
            ..Spec::default()
        }
        .xml()
    ));
    assert_eq!(read_index(BufReader::new(xml.as_bytes())).unwrap().len(), 2);

    let checks: [(LimitAdjust, &str); 5] = [
        (|l| l.max_spectra = 1, "spectra exceed"),
        (|l| l.max_xml_bytes = 64, "byte limit"),
        (|l| l.max_text_bytes = 4, "text exceeds"),
        (|l| l.max_cv_lookups = 0, "vocabulary lookups"),
        (|l| l.max_aux_arrays = 0, "auxiliary arrays"),
    ];
    for (adjust, needle) in checks {
        let mut limits = ImzMLReadLimits::default();
        adjust(&mut limits);
        let error =
            read_index_with_limits(BufReader::new(xml.as_bytes()), &limits).expect_err(needle);
        match error {
            Error::InvalidValue(message) => assert!(message.contains(needle), "{message}"),
            other => panic!("{needle}: {other:?}"),
        }
    }
}

#[test]
fn the_total_auxiliary_array_ceiling_is_enforced() {
    let xml = with_aux(&aux_array("MS:1003006", 1, EXTERNAL));
    let limits = ImzMLReadLimits {
        max_total_aux_arrays: 0,
        ..ImzMLReadLimits::default()
    };
    match read_index_with_limits(BufReader::new(xml.as_bytes()), &limits).unwrap_err() {
        Error::InvalidValue(message) => assert!(message.contains("auxiliary arrays"), "{message}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn group_and_parameter_ceilings_are_enforced() {
    let mut groups = String::new();
    for id in 0..4 {
        groups.push_str(&format!(
            "<referenceableParamGroup id=\"g{id}\">\
             <cvParam accession=\"MS:1000521\" name=\"32-bit float\"/>\
             </referenceableParamGroup>"
        ));
    }
    let xml = document(&format!(
        "<referenceableParamGroupList count=\"4\">{groups}</referenceableParamGroupList>"
    ));
    let limits = ImzMLReadLimits {
        max_param_groups: 2,
        ..ImzMLReadLimits::default()
    };
    match read_index_with_limits(BufReader::new(xml.as_bytes()), &limits).unwrap_err() {
        Error::InvalidValue(message) => assert!(message.contains("groups exceed"), "{message}"),
        other => panic!("{other:?}"),
    }

    let limits = ImzMLReadLimits {
        max_group_params: 2,
        ..ImzMLReadLimits::default()
    };
    match read_index_with_limits(BufReader::new(xml.as_bytes()), &limits).unwrap_err() {
        Error::InvalidValue(message) => {
            assert!(message.contains("parameters exceed"), "{message}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn malformed_xml_is_a_parse_error() {
    for text in [
        "<mzML><spectrum></mzML>",
        "<mzML><spectrum id=\"a\"</spectrum></mzML>",
        "not xml at all <<<",
    ] {
        assert!(
            matches!(
                read_index(BufReader::new(text.as_bytes())),
                Err(Error::Parse { .. })
            ),
            "{text}"
        );
    }
}

// ---------------------------------------------------------------------------
// ImzMLBinaryIO: types, writers and round trips
// ---------------------------------------------------------------------------

#[test]
fn data_type_widths_and_names_match_the_source_spellings() {
    assert_eq!(ImzMLDataType::Float32.width(), Some(4));
    assert_eq!(ImzMLDataType::Float64.width(), Some(8));
    assert_eq!(ImzMLDataType::Int32.width(), Some(4));
    assert_eq!(ImzMLDataType::Int64.width(), Some(8));
    assert_eq!(ImzMLDataType::Unknown.width(), None);
    assert_eq!(ImzMLDataType::Float32.name(), "float32");
    assert_eq!(ImzMLDataType::Float64.name(), "float64");
    assert_eq!(ImzMLDataType::Int32.name(), "int32");
    assert_eq!(ImzMLDataType::Int64.name(), "int64");
    assert_eq!(ImzMLDataType::Unknown.name(), "unknown");
    assert_eq!(ImzMLDataType::default(), ImzMLDataType::Unknown);
    for (accession, expected) in [
        ("MS:1000521", ImzMLDataType::Float32),
        ("MS:1000523", ImzMLDataType::Float64),
        ("MS:1000519", ImzMLDataType::Int32),
        ("MS:1000522", ImzMLDataType::Int64),
    ] {
        assert_eq!(ImzMLDataType::from_accession(accession), Some(expected));
    }
    assert_eq!(ImzMLDataType::from_accession("MS:1000514"), None);
}

/// The four writers emit little-endian payloads that the four readers recover,
/// including the source's narrowing of m/z to float32.
#[test]
fn the_binary_writers_round_trip_through_the_readers() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("round-trip.ibd");
    let limits = ImzMLReadLimits::default();
    let mz = [100.0_f64, 200.5, 300.25];
    let intensity = [1.5_f32, 2.5, 3.5];
    {
        let mut file = std::fs::File::create(&path).unwrap();
        write_mz_as_float64(&mut file, &mz, &limits).unwrap();
        write_mz_as_float32(&mut file, &mz, &limits).unwrap();
        write_float32_array(&mut file, &intensity, &limits).unwrap();
        write_float64_array(&mut file, &[9.5_f64], &limits).unwrap();
    }
    let mut ibd = ImzMLBinaryIO::open(&path).unwrap();
    assert_eq!(ibd.len(), 8 * 3 + 4 * 3 + 4 * 3 + 8);
    assert_eq!(ibd.path(), path);
    assert_eq!(
        ibd.read_mz_array(0, 3, ImzMLDataType::Float64).unwrap(),
        mz.to_vec()
    );
    assert_eq!(
        ibd.read_mz_array(24, 3, ImzMLDataType::Float32).unwrap(),
        mz.to_vec()
    );
    assert_eq!(
        ibd.read_intensity_array(36, 3, ImzMLDataType::Float32)
            .unwrap(),
        intensity.to_vec()
    );
    assert_eq!(
        ibd.read_aux_array(48, 1, ImzMLDataType::Float64, "aux")
            .unwrap(),
        vec![9.5_f32]
    );
    // An empty name falls back to the source's generic label; the read itself
    // is unchanged.
    assert_eq!(
        ibd.read_aux_array(48, 1, ImzMLDataType::Float64, "")
            .unwrap(),
        vec![9.5_f32]
    );
    // Little-endian on every host.
    let raw = std::fs::read(&path).unwrap();
    assert_eq!(&raw[..8], &100.0_f64.to_le_bytes());
}

#[test]
fn integer_arrays_widen_and_narrow_as_the_source_casts_them() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("integers.ibd");
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(-7_i32).to_le_bytes());
    bytes.extend_from_slice(&11_i32.to_le_bytes());
    bytes.extend_from_slice(&(-9_i64).to_le_bytes());
    bytes.extend_from_slice(&13_i64.to_le_bytes());
    std::fs::write(&path, &bytes).unwrap();
    let mut ibd = ImzMLBinaryIO::open(&path).unwrap();
    assert_eq!(
        ibd.read_mz_array(0, 2, ImzMLDataType::Int32).unwrap(),
        vec![-7.0, 11.0]
    );
    assert_eq!(
        ibd.read_mz_array(8, 2, ImzMLDataType::Int64).unwrap(),
        vec![-9.0, 13.0]
    );
    assert_eq!(
        ibd.read_intensity_array(0, 2, ImzMLDataType::Int32)
            .unwrap(),
        vec![-7.0_f32, 11.0]
    );
    assert_eq!(
        ibd.read_intensity_array(8, 2, ImzMLDataType::Int64)
            .unwrap(),
        vec![-9.0_f32, 13.0]
    );
    // An unknown type is not decodable in either reader.
    assert!(matches!(
        ibd.read_mz_array(0, 2, ImzMLDataType::Unknown),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        ibd.read_intensity_array(0, 2, ImzMLDataType::Unknown),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn the_writers_refuse_a_count_above_the_ceiling() {
    let limits = ImzMLReadLimits {
        max_array_elements: 2,
        ..ImzMLReadLimits::default()
    };
    let mut sink = Vec::new();
    assert!(write_float32_array(&mut sink, &[1.0, 2.0], &limits).is_ok());
    assert!(matches!(
        write_float32_array(&mut sink, &[1.0, 2.0, 3.0], &limits),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        write_mz_as_float32(&mut sink, &[1.0, 2.0, 3.0], &limits),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        write_float64_array(&mut sink, &[1.0, 2.0, 3.0], &limits),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        write_mz_as_float64(&mut sink, &[1.0, 2.0, 3.0], &limits),
        Err(Error::InvalidValue(_))
    ));
    // An empty array writes nothing, as the source's early return does.
    let mut empty = Vec::new();
    write_float32_array(&mut empty, &[], &limits).unwrap();
    write_mz_as_float64(&mut empty, &[], &limits).unwrap();
    assert!(empty.is_empty());
}

#[test]
fn an_empty_ibd_is_empty_and_has_no_uuid() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let path = dir.path().join("empty.ibd");
    std::fs::write(&path, []).unwrap();
    let mut ibd = ImzMLBinaryIO::open(&path).unwrap();
    assert!(ibd.is_empty());
    assert_eq!(ibd.len(), 0);
    assert_eq!(ibd.uuid().unwrap(), None);
    // SHA-1 of the empty input, the published constant.
    assert_eq!(
        ibd.sha1_hex().unwrap(),
        "da39a3ee5e6b4b0d3255bfef95601890afd80709"
    );
}

// ---------------------------------------------------------------------------
// Path inference and handler plumbing
// ---------------------------------------------------------------------------

#[test]
fn the_ibd_sibling_is_inferred_case_insensitively() {
    assert_eq!(
        infer_ibd_path("/data/tissue.imzML"),
        PathBuf::from("/data/tissue.ibd")
    );
    assert_eq!(
        infer_ibd_path("/data/tissue.IMZML"),
        PathBuf::from("/data/tissue.ibd")
    );
    assert_eq!(
        infer_ibd_path("/data/tissue.imzml"),
        PathBuf::from("/data/tissue.ibd")
    );
    // Any other name simply gains the extension, as the source does.
    assert_eq!(
        infer_ibd_path("/data/tissue"),
        PathBuf::from("/data/tissue.ibd")
    );
    assert_eq!(
        infer_ibd_path("/data/tissue.mzML"),
        PathBuf::from("/data/tissue.mzML.ibd")
    );
    assert_eq!(infer_ibd_path("short"), PathBuf::from("short.ibd"));
}

#[test]
fn an_explicit_ibd_override_is_what_gets_read() {
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    // The inferred sibling does not exist, so only the override can work.
    let xml = dir.path().join("no-sibling.imzML");
    std::fs::copy(data(PROCESSED), &xml).unwrap();
    assert!(ImzMLHandler::open(&xml).is_err());
    let mut handler = ImzMLHandler::open_with_ibd(&xml, data(PROCESSED_IBD)).unwrap();
    assert_eq!(handler.len(), 9);
    assert_eq!(handler.imzml_path(), xml);
    assert_eq!(handler.ibd_path(), data(PROCESSED_IBD));
    assert_eq!(handler.meta().ibd_file_path, data(PROCESSED_IBD));
    assert_eq!(handler.uuid_status().unwrap(), UuidStatus::Match);
    assert_eq!(handler.spectrum(0).unwrap().spectrum.peaks.len(), 8399);
}

#[test]
fn accessors_expose_the_parsed_index_and_the_configured_limits() {
    let handler = open(CONTINUOUS);
    assert!(!handler.is_empty());
    assert_eq!(handler.index().len(), handler.len());
    assert_eq!(handler.parsed().spectra.len(), 9);
    assert_eq!(handler.parsed().meta, *handler.meta());
    assert_eq!(handler.parsed().get(0), Some(handler.entry(0).unwrap()));
    assert_eq!(handler.parsed().get(9), None);
    assert_eq!(
        handler.limits().max_array_elements,
        ImzMLReadLimits::default().max_array_elements
    );

    let empty = read_index(BufReader::new(document("<run/>").as_bytes())).unwrap();
    assert!(empty.is_empty());
    assert_eq!(empty.len(), 0);
    assert_eq!(empty.index_at_coord(1, 1, 1), None);
    // Nothing was observed, so the bounding box stays at its zero default.
    assert_eq!(empty.meta.max_count_z, 0);
}

#[test]
fn a_missing_file_is_an_io_error() {
    let missing = data("ImzMLFile_does_not_exist.imzML");
    assert!(matches!(ImzMLHandler::open(&missing), Err(Error::Io(_))));
    assert!(matches!(ImzMLBinaryIO::open(&missing), Err(Error::Io(_))));
}
