// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded mzML 1.1 subset I/O. See `docs/MZML_SUPPORT.md` for metadata limitations.
//! Numeric peak arrays are little-endian f32/f64, optionally zlib compressed.
//! This is an event parser, but the returned experiment is held in memory.

use crate::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Precursor,
    SpectrumType,
};
use crate::{Error, Result};
#[path = "mzml_precursor.rs"]
mod precursor_metadata;
use base64::{Engine, engine::general_purpose::STANDARD};
use flate2::{Compression, Decompress, FlushDecompress, Status, write::ZlibEncoder};
use quick_xml::{
    NsReader,
    events::{BytesStart, Event},
    name::ResolveResult,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::{BufRead, Write},
};

const NS: &[u8] = b"http://psi.hupo.org/ms/mzml";
const NAME_KEY: &str = "openms-rust:name";

/// Resource limits apply before allocation from declared lengths and during decoding.
#[derive(Clone, Copy, Debug)]
pub struct ReadOptions {
    /// Maximum XML input bytes, including metadata and encoded arrays.
    pub max_xml_bytes: u64,
    /// Maximum compressed or decoded bytes for each binary array.
    pub max_array_bytes: usize,
    /// Maximum total spectrum and chromatogram peaks in the returned experiment.
    pub max_total_peaks: usize,
    /// Maximum decoded bytes across all primary and auxiliary arrays.
    pub max_total_array_bytes: usize,
    /// Maximum elements across all arrays, including empty string elements.
    pub max_total_array_elements: usize,
    /// Maximum binary arrays, including empty placeholders.
    pub max_total_arrays: usize,
    /// Maximum total spectra plus chromatograms.
    pub max_records: usize,
}
impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            max_xml_bytes: 512 * 1024 * 1024,
            max_array_bytes: 64 * 1024 * 1024,
            max_total_peaks: 10_000_000,
            max_records: 1_000_000,
            max_total_array_bytes: 512 * 1024 * 1024,
            max_total_array_elements: 20_000_000,
            max_total_arrays: 1_000_000,
        }
    }
}

/// Writer stores f64 coordinates and f32 intensities; compression changes only encoding.
#[derive(Clone, Copy, Debug, Default)]
pub struct WriteOptions {
    pub zlib_compression: bool,
}

fn invalid(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}
fn number<T: std::str::FromStr>(value: &str, label: &str) -> Result<T> {
    value
        .parse()
        .map_err(|_| invalid(format!("invalid {label}: {value}")))
}
fn finite(value: &str, label: &str) -> Result<f64> {
    let value: f64 = number(value, label)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(format!("nonfinite {label}")))
    }
}
fn intensity(value: f64) -> Result<f32> {
    let result = value as f32;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(invalid("intensity overflows f32"))
    }
}
fn attributes(
    element: &BytesStart<'_>,
    decoder: quick_xml::encoding::Decoder,
) -> Result<BTreeMap<String, String>> {
    element
        .attributes()
        .map(|attribute| {
            let attribute = attribute.map_err(|e| invalid(e.to_string()))?;
            let key = std::str::from_utf8(attribute.key.as_ref())
                .map_err(|e| invalid(e.to_string()))?
                .to_owned();
            let value = attribute
                .decode_and_unescape_value(decoder)
                .map_err(|e| invalid(e.to_string()))?
                .into_owned();
            xml_string(&value)?;
            Ok((key, value))
        })
        .collect()
}
fn required<'a>(attrs: &'a BTreeMap<String, String>, key: &str) -> Result<&'a str> {
    attrs
        .get(key)
        .map(String::as_str)
        .ok_or_else(|| invalid(format!("missing {key} attribute")))
}
fn seconds(unit: Option<&String>) -> Result<f64> {
    match unit.map(String::as_str) {
        Some("UO:0000010") => Ok(1.0),
        Some("UO:0000031") => Ok(60.0),
        other => Err(Error::Unsupported(format!(
            "time unit {other:?}; explicit seconds/minutes required"
        ))),
    }
}

#[derive(Clone, PartialEq, Eq)]
enum Kind {
    Mz,
    Time,
    Intensity,
    Auxiliary(String),
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum Encoding {
    Float32,
    Float64,
    Int32,
    Int64,
    Ascii,
}
impl Encoding {
    fn width(self) -> Option<usize> {
        match self {
            Self::Float32 | Self::Int32 => Some(4),
            Self::Float64 | Self::Int64 => Some(8),
            Self::Ascii => None,
        }
    }
}
enum Values {
    Floats(Vec<f64>),
    Integers(Vec<i32>),
    Strings(Vec<String>),
}
#[derive(Default)]
struct Binary {
    kind: Option<Kind>,
    encoding: Option<Encoding>,
    compressed: Option<bool>,
    encoded: String,
    encoded_length: usize,
    array_length: Option<usize>,
    time_scale: f64,
    has_binary: bool,
}
impl Binary {
    fn cv(&mut self, attrs: &BTreeMap<String, String>) -> Result<()> {
        let accession = required(attrs, "accession")?;
        let kind = match accession {
            "MS:1000514" => Some(Kind::Mz),
            "MS:1000515" => Some(Kind::Intensity),
            "MS:1000595" => {
                self.time_scale = seconds(attrs.get("unitAccession"))?;
                Some(Kind::Time)
            }
            "MS:1000786" => {
                let name = required(attrs, "value")?;
                if name.is_empty() {
                    return Err(invalid("auxiliary array name is empty"));
                }
                if ["unitAccession", "unitCvRef", "unitName"]
                    .iter()
                    .any(|k| attrs.contains_key(*k))
                {
                    return Err(Error::Unsupported(
                        "units on auxiliary arrays are not represented".into(),
                    ));
                }
                Some(Kind::Auxiliary(name.to_owned()))
            }
            _ => None,
        };
        if let Some(kind) = kind {
            if self.kind.replace(kind).is_some() {
                return Err(invalid("multiple binary array types"));
            }
        } else {
            if ["unitAccession", "unitCvRef", "unitName"]
                .iter()
                .any(|key| attrs.contains_key(*key))
            {
                return Err(Error::Unsupported(
                    "units on binary encoding/compression terms are not represented".into(),
                ));
            }
            let encoding = match accession {
                "MS:1000521" => Some(Encoding::Float32),
                "MS:1000523" => Some(Encoding::Float64),
                "MS:1000519" => Some(Encoding::Int32),
                "MS:1000522" => Some(Encoding::Int64),
                "MS:1001479" => Some(Encoding::Ascii),
                _ => None,
            };
            if let Some(encoding) = encoding {
                if self.encoding.replace(encoding).is_some() {
                    return Err(invalid("multiple binary precisions"));
                }
            } else if matches!(accession, "MS:1000574" | "MS:1000576") {
                if self.compressed.replace(accession == "MS:1000574").is_some() {
                    return Err(invalid("multiple binary compression terms"));
                }
            } else {
                return Err(Error::Unsupported(format!("binary array CV {accession}")));
            }
        }
        Ok(())
    }
    fn decode(
        self,
        default_count: usize,
        options: &ReadOptions,
        remaining_bytes: &mut usize,
        remaining_elements: &mut usize,
    ) -> Result<(Kind, Values)> {
        if !self.has_binary {
            return Err(invalid("missing binary element"));
        }
        let kind = self
            .kind
            .ok_or_else(|| invalid("missing binary array type"))?;
        let encoding = self
            .encoding
            .ok_or_else(|| invalid("missing binary precision"))?;
        let compressed = self
            .compressed
            .ok_or_else(|| invalid("missing binary compression term"))?;
        let auxiliary = matches!(kind, Kind::Auxiliary(_));
        if !auxiliary && !matches!(encoding, Encoding::Float32 | Encoding::Float64) {
            return Err(Error::Unsupported(
                "primary peak arrays must use floating-point encoding".into(),
            ));
        }
        let count = self.array_length.unwrap_or(default_count);
        if count != default_count && !(auxiliary && count == 0) {
            return Err(invalid(
                "nonempty arrayLength differs from defaultArrayLength",
            ));
        }
        *remaining_elements = remaining_elements
            .checked_sub(count)
            .ok_or_else(|| invalid("total binary element limit exceeded"))?;
        let byte_limit = options.max_array_bytes.min(*remaining_bytes);
        let expected = encoding
            .width()
            .map(|width| {
                count
                    .checked_mul(width)
                    .filter(|&n| n <= byte_limit)
                    .ok_or_else(|| invalid("binary array exceeds configured byte limit"))
            })
            .transpose()?;
        // Every string has at least its NUL terminator. Check before allocating.
        if expected.is_none() && count > byte_limit {
            return Err(invalid("string array exceeds configured byte limit"));
        }
        if self.encoded.len() != self.encoded_length {
            return Err(invalid(
                "encodedLength does not match base64 character count",
            ));
        }
        let bytes = STANDARD
            .decode(self.encoded.as_bytes())
            .map_err(|e| invalid(format!("invalid base64: {e}")))?;
        if bytes.len() > options.max_array_bytes {
            return Err(invalid("encoded array exceeds configured byte limit"));
        }
        let limit = expected.unwrap_or(byte_limit);
        let decoded = if compressed {
            // Bounded chunks also handle string arrays, whose byte lengths are
            // variable. Reject concatenated/truncated streams and expansion past
            // either the declared numeric length or cumulative storage budget.
            let mut decoder = Decompress::new(true);
            let mut output = Vec::new();
            let mut chunk = [0u8; 8192];
            loop {
                let before_in = decoder.total_in();
                let before_out = decoder.total_out();
                let capacity = limit
                    .saturating_sub(output.len())
                    .saturating_add(1)
                    .min(chunk.len());
                let status = decoder
                    .decompress(
                        &bytes[before_in as usize..],
                        &mut chunk[..capacity],
                        FlushDecompress::None,
                    )
                    .map_err(|e| invalid(format!("invalid zlib array: {e}")))?;
                let written = (decoder.total_out() - before_out) as usize;
                if output.len().checked_add(written).is_none_or(|n| n > limit) {
                    return Err(invalid("decoded binary array exceeds byte limit"));
                }
                output.extend_from_slice(&chunk[..written]);
                if status == Status::StreamEnd {
                    if decoder.total_in() != bytes.len() as u64 {
                        return Err(invalid("trailing compressed data"));
                    }
                    break;
                }
                if decoder.total_in() == before_in && written == 0 {
                    return Err(invalid("incomplete zlib stream"));
                }
            }
            output
        } else {
            bytes
        };
        if decoded.len() > limit || expected.is_some_and(|n| decoded.len() != n) {
            return Err(invalid(
                "binary byte length differs from declared length or exceeds limit",
            ));
        }
        *remaining_bytes = remaining_bytes
            .checked_sub(decoded.len())
            .ok_or_else(|| invalid("total binary byte limit exceeded"))?;
        let values = match encoding {
            Encoding::Ascii => {
                if !decoded.is_ascii() || (!decoded.is_empty() && decoded.last() != Some(&0)) {
                    return Err(invalid(
                        "string array must contain NUL-terminated ASCII strings",
                    ));
                }
                if decoded.iter().filter(|&&b| b == 0).count() != count {
                    return Err(invalid(
                        "string array element count differs from declared length",
                    ));
                }
                let strings = if decoded.is_empty() {
                    Vec::new()
                } else {
                    decoded[..decoded.len() - 1]
                        .split(|&b| b == 0)
                        .map(|value| String::from_utf8(value.to_vec()).expect("validated ASCII"))
                        .collect()
                };
                Values::Strings(strings)
            }
            Encoding::Int32 | Encoding::Int64 => {
                let width = encoding.width().unwrap();
                let mut values = Vec::with_capacity(count);
                for chunk in decoded.chunks_exact(width) {
                    let value = if width == 4 {
                        i64::from(i32::from_le_bytes(chunk.try_into().unwrap()))
                    } else {
                        i64::from_le_bytes(chunk.try_into().unwrap())
                    };
                    values.push(
                        i32::try_from(value)
                            .map_err(|_| invalid("integer array value overflows native i32"))?,
                    );
                }
                Values::Integers(values)
            }
            Encoding::Float32 | Encoding::Float64 => {
                let width = encoding.width().unwrap();
                let scale = if kind == Kind::Time {
                    self.time_scale
                } else {
                    1.0
                };
                let mut values = Vec::with_capacity(count);
                for chunk in decoded.chunks_exact(width) {
                    let value = if width == 4 {
                        f64::from(f32::from_le_bytes(chunk.try_into().unwrap()))
                    } else {
                        f64::from_le_bytes(chunk.try_into().unwrap())
                    } * scale;
                    if !value.is_finite() {
                        return Err(invalid("nonfinite binary value"));
                    }
                    values.push(value);
                }
                Values::Floats(values)
            }
        };
        Ok((kind, values))
    }
}

struct Record {
    spectrum: Option<MSSpectrum>,
    chromatogram: Option<MSChromatogram>,
    count: usize,
    array_names: BTreeSet<String>,
    coordinates: Option<Vec<f64>>,
    intensities: Option<Vec<f64>>,
    precursor: Option<Precursor>,
    selected_ion: bool,
    selected_mz: bool,
    rt_seen: bool,
    precursors_seen: usize,
    name_seen: bool,
    seen_fields: BTreeSet<&'static str>,
    precursor_fields: BTreeSet<&'static str>,
    precursor_list_seen: bool,
    selected_list_seen: bool,
    scan_list_seen: bool,
}
impl Record {
    fn metadata(&mut self) -> &mut BTreeMap<String, String> {
        if let Some(s) = &mut self.spectrum {
            &mut s.metadata
        } else {
            &mut self.chromatogram.as_mut().unwrap().metadata
        }
    }
    fn cv(&mut self, parent: &str, attrs: &BTreeMap<String, String>) -> Result<()> {
        let accession = required(attrs, "accession")?;
        let value = attrs.get("value").map(String::as_str).unwrap_or("");
        let field = match (parent, accession) {
            ("spectrum", "MS:1000511") => Some("ms_level"),
            ("spectrum", "MS:1000127" | "MS:1000128") => Some("spectrum_type"),
            ("selectedIon", "MS:1000744" | "MS:1000040") => Some("precursor_mz"),
            ("selectedIon", "MS:1000041") => Some("precursor_charge"),
            ("selectedIon", "MS:1000042") => Some("precursor_intensity"),
            ("isolationWindow", "MS:1000827") if self.precursor.is_some() => {
                Some("isolation_target")
            }
            _ if self.precursor.is_some() => precursor_metadata::field(parent, accession),
            _ => None,
        };
        if let Some(field) = field {
            let seen = if matches!(parent, "selectedIon" | "isolationWindow" | "activation") {
                &mut self.precursor_fields
            } else {
                &mut self.seen_fields
            };
            if !seen.insert(field) {
                return Err(invalid(format!("duplicate scientific field {field}")));
            }
        }
        if let Some(p) = self.precursor.as_mut() {
            if precursor_metadata::read_cv(p, parent, accession, value, attrs, self.selected_mz)? {
                return Ok(());
            }
        }
        match (parent, accession) {
            ("spectrum", "MS:1000511") => {
                self.spectrum.as_mut().unwrap().ms_level = number(value, "MS level")?
            }
            ("spectrum", "MS:1000127") => {
                self.spectrum.as_mut().unwrap().spectrum_type = SpectrumType::Centroid
            }
            ("spectrum", "MS:1000128") => {
                self.spectrum.as_mut().unwrap().spectrum_type = SpectrumType::Profile
            }
            ("scan", "MS:1000016") => {
                if self.rt_seen {
                    return Err(Error::Unsupported(
                        "multiple scan start times in one spectrum".into(),
                    ));
                }
                let rt = finite(value, "scan start time")? * seconds(attrs.get("unitAccession"))?;
                if !rt.is_finite() {
                    return Err(invalid("scan start time overflow"));
                }
                self.spectrum
                    .as_mut()
                    .ok_or_else(|| invalid("scan outside spectrum"))?
                    .rt = rt;
                self.rt_seen = true;
            }
            ("selectedIon", "MS:1000744" | "MS:1000040") => {
                self.precursor
                    .as_mut()
                    .ok_or_else(|| invalid("selected ion outside precursor"))?
                    .mz = finite(value, "precursor m/z")?;
                self.selected_mz = true;
            }
            ("selectedIon", "MS:1000041") => {
                self.precursor
                    .as_mut()
                    .ok_or_else(|| invalid("selected ion outside precursor"))?
                    .charge = number(value, "precursor charge")?
            }
            ("selectedIon", "MS:1000042") => {
                self.precursor
                    .as_mut()
                    .ok_or_else(|| invalid("selected ion outside precursor"))?
                    .intensity = intensity(finite(value, "precursor intensity")?)?
            }
            _ => {} // Acquisition CVs outside the supported model are intentionally not retained.
        }
        Ok(())
    }
    fn finish(mut self, experiment: &mut MSExperiment) -> Result<()> {
        let (positions, intensities) = match (self.coordinates.take(), self.intensities.take()) {
            (Some(p), Some(i)) => (p, i),
            (None, None) if self.count == 0 => (Vec::new(), Vec::new()),
            _ => return Err(invalid("missing coordinate or intensity array")),
        };
        if let Some(mut spectrum) = self.spectrum {
            spectrum.peaks = positions
                .into_iter()
                .zip(intensities)
                .map(|(mz, i)| Ok(Peak1D::new(mz, intensity(i)?)))
                .collect::<Result<_>>()?;
            spectrum.validate()?;
            experiment.spectra.push(spectrum);
        } else {
            let mut chromatogram = self.chromatogram.unwrap();
            chromatogram.peaks = positions
                .into_iter()
                .zip(intensities)
                .map(|(rt, i)| Ok(ChromatogramPeak::new(rt, intensity(i)?)))
                .collect::<Result<_>>()?;
            chromatogram.validate()?;
            experiment.chromatograms.push(chromatogram);
        }
        Ok(())
    }
}

/// Read a complete experiment with default resource limits.
pub fn read(reader: impl BufRead) -> Result<MSExperiment> {
    read_with_options(reader, &ReadOptions::default())
}

/// Read the documented mzML subset, rejecting unsupported binary encodings.
/// Errors have line zero when the XML parser cannot supply a line number.
pub fn read_with_options(reader: impl BufRead, options: &ReadOptions) -> Result<MSExperiment> {
    let limit = options
        .max_xml_bytes
        .checked_add(1)
        .ok_or_else(|| Error::InvalidValue("XML byte limit must be below u64::MAX".into()))?;
    let mut reader = NsReader::from_reader(reader.take(limit));
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().enable_all_checks(true);
    let mut buffer = Vec::new();
    let mut stack = Vec::<String>::new();
    let mut experiment = MSExperiment::new();
    let mut record: Option<Record> = None;
    let mut binary: Option<Binary> = None;
    let mut seen_root = false;
    let mut seen_mzml = false;
    let mut seen_run = false;
    let mut total_peaks = 0usize;
    let mut records = 0usize;
    let mut spectrum_count = None;
    let mut chromatogram_count = None;
    let mut array_count = None;
    let mut arrays_seen = 0usize;
    let mut total_arrays = 0usize;
    let mut remaining_array_bytes = options.max_total_array_bytes;
    let mut remaining_array_elements = options.max_total_array_elements;
    let mut ids = BTreeSet::new();
    let mut counted_lists: Vec<(usize, &str, usize, usize)> = Vec::new();
    let mut seen_declaration = false;
    let mut ascii_only = false;
    loop {
        let decoder = reader.decoder();
        let (namespace, event) = reader
            .read_resolved_event_into(&mut buffer)
            .map_err(|e| invalid(e.to_string()))?;
        let namespace_ok = matches!(namespace, ResolveResult::Bound(ns) if ns.as_ref() == NS);
        if ascii_only && !event.is_ascii() {
            return Err(invalid("non-ASCII bytes in US-ASCII XML"));
        }
        if reader.buffer_position() > options.max_xml_bytes {
            return Err(invalid("XML exceeds configured byte limit"));
        }
        match event {
            Event::Start(element) => {
                if !namespace_ok {
                    return Err(invalid("element is outside the mzML namespace"));
                }
                let tag = std::str::from_utf8(element.local_name().as_ref())
                    .map_err(|e| invalid(e.to_string()))?
                    .to_owned();
                let attrs = attributes(&element, decoder)?;
                let parent = stack.last().map(String::as_str).unwrap_or("");
                if stack.is_empty() {
                    if seen_root || !matches!(tag.as_str(), "mzML" | "indexedmzML") {
                        return Err(invalid("expected a single mzML or indexedmzML root"));
                    }
                    seen_root = true;
                }
                if stack.len() >= 128 {
                    return Err(invalid("XML nesting exceeds 128 levels"));
                }
                if let Some((depth, child, _, actual)) = counted_lists.last_mut() {
                    if *depth == stack.len() && *child == tag {
                        *actual += 1;
                    }
                }
                match tag.as_str() {
                    "indexedmzML" if !parent.is_empty() => {
                        return Err(invalid("nested indexedmzML wrapper"));
                    }
                    "isolationWindow" | "activation" => {
                        if parent == "product" && tag == "isolationWindow" {
                            if record.as_ref().is_some_and(|r| r.precursor.is_some()) {
                                return Err(invalid("product nested inside precursor"));
                            }
                        } else {
                            if parent != "precursor" {
                                return Err(invalid("misplaced precursor acquisition element"));
                            }
                            let r = record
                                .as_mut()
                                .filter(|r| r.precursor.is_some())
                                .ok_or_else(|| invalid("acquisition element outside precursor"))?;
                            let key = if tag == "activation" {
                                "activation_element"
                            } else {
                                "isolation_element"
                            };
                            if !r.precursor_fields.insert(key) {
                                return Err(invalid("duplicate precursor acquisition element"));
                            }
                        }
                    }
                    "precursorList" | "selectedIonList" | "scanList" => {
                        let expected_parent = match tag.as_str() {
                            "precursorList" | "scanList" => "spectrum",
                            _ => "precursor",
                        };
                        if parent != expected_parent {
                            return Err(invalid("misplaced precursor/scan list"));
                        }
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("list outside spectrum"))?;
                        let seen = match tag.as_str() {
                            "precursorList" => &mut r.precursor_list_seen,
                            "selectedIonList" => &mut r.selected_list_seen,
                            _ => &mut r.scan_list_seen,
                        };
                        if *seen {
                            return Err(invalid("duplicate precursor/scan list"));
                        }
                        *seen = true;
                        let expected: usize = number(required(&attrs, "count")?, "list count")?;
                        let child = match tag.as_str() {
                            "precursorList" => "precursor",
                            "selectedIonList" => "selectedIon",
                            _ => "scan",
                        };
                        counted_lists.push((stack.len() + 1, child, expected, 0));
                    }
                    "mzML" => {
                        if seen_mzml || !matches!(parent, "" | "indexedmzML") {
                            return Err(invalid("unexpected mzML element"));
                        }
                        if !required(&attrs, "version")?.starts_with("1.1.") {
                            return Err(Error::Unsupported("only mzML 1.1 is supported".into()));
                        }
                        seen_mzml = true;
                    }
                    "run" => {
                        if parent != "mzML" || seen_run {
                            return Err(invalid("expected exactly one mzML run"));
                        }
                        seen_run = true;
                    }
                    "referenceableParamGroup" | "referenceableParamGroupRef" => {
                        return Err(Error::Unsupported(
                            "mzML referenceable parameter groups".into(),
                        ));
                    }
                    "spectrumList" | "chromatogramList" => {
                        if parent != "run" {
                            return Err(invalid("record list outside run"));
                        }
                        let slot = if tag == "spectrumList" {
                            &mut spectrum_count
                        } else {
                            &mut chromatogram_count
                        };
                        if slot.is_some() {
                            return Err(invalid("duplicate record list"));
                        }
                        *slot = Some(number::<usize>(required(&attrs, "count")?, "record count")?);
                    }
                    "spectrum" | "chromatogram" => {
                        let expected_parent = if tag == "spectrum" {
                            "spectrumList"
                        } else {
                            "chromatogramList"
                        };
                        if parent != expected_parent || record.is_some() {
                            return Err(invalid("misplaced spectrum/chromatogram"));
                        }
                        records = records
                            .checked_add(1)
                            .filter(|&n| n <= options.max_records)
                            .ok_or_else(|| invalid("record count exceeds configured limit"))?;
                        let count: usize = number(
                            required(&attrs, "defaultArrayLength")?,
                            "defaultArrayLength",
                        )?;
                        total_peaks = total_peaks
                            .checked_add(count)
                            .filter(|&n| n <= options.max_total_peaks)
                            .ok_or_else(|| invalid("peak count exceeds configured limit"))?;
                        let id = required(&attrs, "id")?.to_owned();
                        if id.is_empty() || !ids.insert((tag.clone(), id.clone())) {
                            return Err(invalid("empty or duplicate record id"));
                        }
                        record = Some(Record {
                            spectrum: (tag == "spectrum").then(|| MSSpectrum {
                                native_id: id.clone(),
                                ..Default::default()
                            }),
                            chromatogram: (tag == "chromatogram").then(|| MSChromatogram {
                                native_id: id,
                                ..Default::default()
                            }),
                            count,
                            array_names: BTreeSet::new(),
                            coordinates: None,
                            intensities: None,
                            precursor: None,
                            selected_ion: false,
                            selected_mz: false,
                            rt_seen: false,
                            precursors_seen: 0,
                            name_seen: false,
                            seen_fields: BTreeSet::new(),
                            precursor_fields: BTreeSet::new(),
                            precursor_list_seen: false,
                            selected_list_seen: false,
                            scan_list_seen: false,
                        });
                        array_count = None;
                        arrays_seen = 0;
                    }
                    "binaryDataArrayList" => {
                        if !matches!(parent, "spectrum" | "chromatogram")
                            || record.is_none()
                            || array_count.is_some()
                        {
                            return Err(invalid("misplaced/duplicate binary array list"));
                        }
                        array_count = Some(number::<usize>(
                            required(&attrs, "count")?,
                            "binary array count",
                        )?);
                    }
                    "binaryDataArray" => {
                        if parent != "binaryDataArrayList" || binary.is_some() {
                            return Err(invalid("misplaced binary array"));
                        }
                        let encoded_length: usize =
                            number(required(&attrs, "encodedLength")?, "encodedLength")?;
                        let max_encoded = options
                            .max_array_bytes
                            .saturating_add(2)
                            .saturating_div(3)
                            .saturating_mul(4);
                        if encoded_length > max_encoded {
                            return Err(invalid("encodedLength exceeds configured byte limit"));
                        }
                        binary = Some(Binary {
                            encoded_length,
                            array_length: attrs
                                .get("arrayLength")
                                .map(|n| number(n, "arrayLength"))
                                .transpose()?,
                            ..Default::default()
                        });
                        arrays_seen += 1;
                        total_arrays = total_arrays
                            .checked_add(1)
                            .filter(|&n| n <= options.max_total_arrays)
                            .ok_or_else(|| invalid("total binary array count limit exceeded"))?;
                    }
                    "binary" => {
                        if parent != "binaryDataArray" {
                            return Err(invalid("misplaced binary element"));
                        }
                        let b = binary
                            .as_mut()
                            .ok_or_else(|| invalid("binary element outside array"))?;
                        if b.has_binary {
                            return Err(invalid("duplicate binary element"));
                        }
                        b.has_binary = true;
                    }
                    "precursor" => {
                        if !matches!(parent, "precursorList" | "chromatogram") {
                            return Err(invalid("misplaced precursor"));
                        }
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("precursor outside record"))?;
                        r.precursors_seen += 1;
                        if r.chromatogram.is_some() && r.precursors_seen > 1 {
                            return Err(invalid("multiple chromatogram precursors"));
                        }
                        if r.precursor.replace(Precursor::default()).is_some() {
                            return Err(invalid("nested precursor"));
                        }
                        r.precursor.as_mut().unwrap().spectrum_reference =
                            attrs.get("spectrumRef").cloned();
                        r.selected_ion = false;
                        r.selected_mz = false;
                        r.selected_list_seen = false;
                        r.precursor_fields.clear();
                    }
                    "selectedIon" => {
                        if parent != "selectedIonList" {
                            return Err(invalid("misplaced selected ion"));
                        }
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("selected ion outside record"))?;
                        if r.precursor.is_none() {
                            return Err(invalid("selected ion outside precursor"));
                        }
                        if r.selected_ion {
                            return Err(Error::Unsupported(
                                "multiple selected ions in one precursor".into(),
                            ));
                        }
                        r.selected_ion = true;
                    }
                    "cvParam" => {
                        if parent == "binaryDataArray" {
                            binary
                                .as_mut()
                                .ok_or_else(|| invalid("CV outside binary array"))?
                                .cv(&attrs)?;
                        } else if let Some(r) = &mut record {
                            r.cv(parent, &attrs)?;
                        }
                    }
                    "userParam" if parent == "binaryDataArray" => {
                        return Err(Error::Unsupported(
                            "metadata on binary arrays is not represented".into(),
                        ));
                    }
                    "userParam" if matches!(parent, "run" | "spectrum" | "chromatogram") => {
                        let name = required(&attrs, "name")?.to_owned();
                        let value = attrs.get("value").cloned().unwrap_or_default();
                        if name == NAME_KEY && parent == "run" {
                            return Err(invalid("reserved record name userParam at run level"));
                        }
                        if name == NAME_KEY && parent != "run" {
                            let r = record
                                .as_mut()
                                .ok_or_else(|| invalid("name outside record"))?;
                            if r.name_seen {
                                return Err(invalid("duplicate record name userParam"));
                            }
                            r.name_seen = true;
                            if let Some(s) = &mut r.spectrum {
                                s.name = value;
                            } else {
                                r.chromatogram.as_mut().unwrap().name = value;
                            }
                        } else {
                            let metadata = if parent == "run" {
                                &mut experiment.metadata
                            } else {
                                record
                                    .as_mut()
                                    .ok_or_else(|| invalid("metadata outside record"))?
                                    .metadata()
                            };
                            if metadata.insert(name.clone(), value).is_some() {
                                return Err(Error::Unsupported(format!(
                                    "duplicate userParam name {name}"
                                )));
                            }
                        }
                    }
                    _ => {}
                }
                stack.push(tag);
            }
            Event::End(_) => {
                if counted_lists
                    .last()
                    .is_some_and(|(depth, _, _, _)| *depth == stack.len())
                {
                    let (_, child, expected, actual) = counted_lists.pop().unwrap();
                    if expected != actual {
                        return Err(invalid(format!("{child} list count mismatch")));
                    }
                }
                let tag = stack
                    .pop()
                    .ok_or_else(|| invalid("unmatched closing tag"))?;
                match tag.as_str() {
                    "binaryDataArray" => {
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("array outside record"))?;
                        let (kind, values) = binary
                            .take()
                            .ok_or_else(|| invalid("missing binary array state"))?
                            .decode(
                                r.count,
                                options,
                                &mut remaining_array_bytes,
                                &mut remaining_array_elements,
                            )?;
                        if let Kind::Auxiliary(name) = kind {
                            if !r.array_names.insert(name.clone()) {
                                return Err(invalid("duplicate auxiliary array name"));
                            }
                            let (floats, integers, strings) =
                                if let Some(spectrum) = &mut r.spectrum {
                                    (
                                        &mut spectrum.float_data_arrays,
                                        &mut spectrum.integer_data_arrays,
                                        &mut spectrum.string_data_arrays,
                                    )
                                } else {
                                    let chromatogram = r.chromatogram.as_mut().unwrap();
                                    (
                                        &mut chromatogram.float_data_arrays,
                                        &mut chromatogram.integer_data_arrays,
                                        &mut chromatogram.string_data_arrays,
                                    )
                                };
                            match values {
                                Values::Floats(values) => floats.push(DataArray::new(
                                    name,
                                    values.into_iter().map(intensity).collect::<Result<_>>()?,
                                )),
                                Values::Integers(values) => {
                                    integers.push(DataArray::new(name, values))
                                }
                                Values::Strings(values) => {
                                    strings.push(DataArray::new(name, values))
                                }
                            }
                        } else {
                            let Values::Floats(values) = values else {
                                return Err(invalid("non-floating primary binary array"));
                            };
                            let slot = match kind {
                                Kind::Intensity => &mut r.intensities,
                                Kind::Mz if r.spectrum.is_some() => &mut r.coordinates,
                                Kind::Time if r.chromatogram.is_some() => &mut r.coordinates,
                                _ => {
                                    return Err(invalid(
                                        "coordinate array has wrong type for record",
                                    ));
                                }
                            };
                            if slot.replace(values).is_some() {
                                return Err(invalid("duplicate coordinate/intensity array"));
                            }
                        }
                    }
                    "binaryDataArrayList" => {
                        if array_count != Some(arrays_seen) {
                            return Err(invalid("binary array count mismatch"));
                        }
                    }
                    "precursor" => {
                        let r = record
                            .as_mut()
                            .ok_or_else(|| invalid("precursor outside record"))?;
                        let mut p = r
                            .precursor
                            .take()
                            .ok_or_else(|| invalid("missing precursor state"))?;
                        if p.isolation_target_mz == Some(p.mz) {
                            p.isolation_target_mz = None;
                        }
                        p.validate()?;
                        if let Some(s) = &mut r.spectrum {
                            s.precursors.push(p);
                        } else {
                            r.chromatogram.as_mut().unwrap().precursor = p;
                        }
                    }
                    "spectrum" | "chromatogram" => record
                        .take()
                        .ok_or_else(|| invalid("missing record state"))?
                        .finish(&mut experiment)?,
                    _ => {}
                }
            }
            Event::Text(text) => {
                let text = text.decode().map_err(|e| invalid(e.to_string()))?;
                xml_string(&text)?;
                if stack.last().is_some_and(|tag| tag == "binary") {
                    let b = binary
                        .as_mut()
                        .ok_or_else(|| invalid("text outside binary array"))?;
                    for byte in text.bytes().filter(|b| !b.is_ascii_whitespace()) {
                        if b.encoded.len() >= b.encoded_length {
                            return Err(invalid("binary text exceeds encodedLength"));
                        }
                        if !byte.is_ascii() {
                            return Err(invalid("non-ASCII base64 text"));
                        }
                        b.encoded.push(char::from(byte));
                    }
                } else if stack.is_empty() && !text.trim().is_empty() {
                    return Err(invalid("text outside XML root"));
                }
            }
            Event::Decl(declaration) => {
                if seen_root || seen_declaration {
                    return Err(invalid("misplaced or duplicate XML declaration"));
                }
                seen_declaration = true;
                if declaration
                    .version()
                    .map_err(|e| invalid(e.to_string()))?
                    .as_ref()
                    != b"1.0"
                {
                    return Err(Error::Unsupported("only XML 1.0 is supported".into()));
                }
                if let Some(encoding) = declaration.encoding() {
                    let encoding = encoding.map_err(|e| invalid(e.to_string()))?;
                    ascii_only = encoding.eq_ignore_ascii_case(b"US-ASCII");
                    if !encoding.eq_ignore_ascii_case(b"UTF-8")
                        && !encoding.eq_ignore_ascii_case(b"US-ASCII")
                    {
                        return Err(Error::Unsupported(
                            "only UTF-8/US-ASCII mzML XML encodings are supported".into(),
                        ));
                    }
                }
            }
            Event::Comment(text) => {
                xml_string(&text.decode().map_err(|e| invalid(e.to_string()))?)?;
            }
            Event::PI(text) => {
                xml_string(
                    std::str::from_utf8(text.as_ref()).map_err(|e| invalid(e.to_string()))?,
                )?;
            }
            Event::DocType(_) => {
                return Err(Error::Unsupported("XML DTDs are not supported".into()));
            }
            Event::GeneralRef(_) | Event::CData(_) => {
                return Err(Error::Unsupported(
                    "XML entity references in text and CDATA are not supported".into(),
                ));
            }
            Event::Eof => break,
            _ => {}
        }
        buffer.clear();
    }
    if !seen_mzml || !seen_run || !stack.is_empty() || record.is_some() {
        return Err(invalid("incomplete mzML document"));
    }
    if spectrum_count.is_some_and(|n| n != experiment.spectra.len())
        || chromatogram_count.is_some_and(|n| n != experiment.chromatograms.len())
    {
        return Err(invalid("declared record count mismatch"));
    }
    Ok(experiment)
}

fn escape(value: &str) -> String {
    quick_xml::escape::escape(value)
        .replace('\n', "&#10;")
        .replace('\r', "&#13;")
        .replace('\t', "&#9;")
}
fn xml_string(value: &str) -> Result<()> {
    if value.chars().all(|c| {
        matches!(c, '\t' | '\n' | '\r')
            || ('\u{20}'..='\u{d7ff}').contains(&c)
            || ('\u{e000}'..='\u{fffd}').contains(&c)
            || c >= '\u{10000}'
    }) {
        Ok(())
    } else {
        Err(Error::InvalidValue(
            "text contains invalid XML 1.0 characters".into(),
        ))
    }
}
fn user_params(w: &mut impl Write, metadata: &BTreeMap<String, String>, name: &str) -> Result<()> {
    if !name.is_empty() {
        writeln!(
            w,
            "<userParam name=\"{NAME_KEY}\" value=\"{}\" type=\"xsd:string\"/>",
            escape(name)
        )?;
    }
    for (key, value) in metadata {
        writeln!(
            w,
            "<userParam name=\"{}\" value=\"{}\" type=\"xsd:string\"/>",
            escape(key),
            escape(value)
        )?;
    }
    Ok(())
}
fn cv(w: &mut impl Write, accession: &str, name: &str, value: &str, unit: &str) -> Result<()> {
    writeln!(
        w,
        "<cvParam cvRef=\"MS\" accession=\"{accession}\" name=\"{name}\" value=\"{}\"{unit}/>",
        escape(value)
    )?;
    Ok(())
}
const SECOND: &str = " unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"";
fn write_precursor(w: &mut impl Write, precursor: &Precursor) -> Result<()> {
    precursor_metadata::write_start(w, precursor)?;
    cv(
        w,
        "MS:1000744",
        "selected ion m/z",
        &precursor.mz.to_string(),
        " unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"",
    )?;
    if precursor.charge != 0 {
        cv(
            w,
            "MS:1000041",
            "charge state",
            &precursor.charge.to_string(),
            "",
        )?;
    }
    cv(
        w,
        "MS:1000042",
        "peak intensity",
        &precursor.intensity.to_string(),
        "",
    )?;
    precursor_metadata::write_end(w, precursor)
}
fn write_array(
    w: &mut impl Write,
    bytes: Vec<u8>,
    kind: Kind,
    encoding: Encoding,
    array_length: Option<usize>,
    options: &WriteOptions,
) -> Result<()> {
    let bytes = if options.zlib_compression {
        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(&bytes)?;
        encoder.finish()?
    } else {
        bytes
    };
    let encoded = STANDARD.encode(bytes);
    write!(w, "<binaryDataArray encodedLength=\"{}\"", encoded.len())?;
    if let Some(length) = array_length {
        write!(w, " arrayLength=\"{length}\"")?;
    }
    writeln!(w, ">")?;
    let (accession, name) = match encoding {
        Encoding::Float32 => ("MS:1000521", "32-bit float"),
        Encoding::Float64 => ("MS:1000523", "64-bit float"),
        Encoding::Int32 => ("MS:1000519", "32-bit integer"),
        Encoding::Int64 => ("MS:1000522", "64-bit integer"),
        Encoding::Ascii => ("MS:1001479", "null-terminated ASCII string"),
    };
    cv(w, accession, name, "", "")?;
    if options.zlib_compression {
        cv(w, "MS:1000574", "zlib compression", "", "")?;
    } else {
        cv(w, "MS:1000576", "no compression", "", "")?;
    }
    match kind {
        Kind::Mz => cv(
            w,
            "MS:1000514",
            "m/z array",
            "",
            " unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"",
        )?,
        Kind::Time => cv(w, "MS:1000595", "time array", "", SECOND)?,
        Kind::Intensity => cv(w, "MS:1000515", "intensity array", "", "")?,
        Kind::Auxiliary(name) => cv(w, "MS:1000786", "non-standard data array", &name, "")?,
    }
    writeln!(w, "<binary>{encoded}</binary></binaryDataArray>")?;
    Ok(())
}

fn write_auxiliary_arrays(
    w: &mut impl Write,
    floats: &[DataArray<f32>],
    integers: &[DataArray<i32>],
    strings: &[DataArray<String>],
    options: &WriteOptions,
) -> Result<()> {
    for array in floats {
        write_array(
            w,
            array.data.iter().flat_map(|v| v.to_le_bytes()).collect(),
            Kind::Auxiliary(array.name.clone()),
            Encoding::Float32,
            Some(array.data.len()),
            options,
        )?;
    }
    for array in integers {
        // OpenMS writes its integer annotations in signed 64-bit representation.
        write_array(
            w,
            array
                .data
                .iter()
                .flat_map(|&v| i64::from(v).to_le_bytes())
                .collect(),
            Kind::Auxiliary(array.name.clone()),
            Encoding::Int64,
            Some(array.data.len()),
            options,
        )?;
    }
    for array in strings {
        write_array(
            w,
            array
                .data
                .iter()
                .flat_map(|v| v.bytes().chain(std::iter::once(0)))
                .collect(),
            Kind::Auxiliary(array.name.clone()),
            Encoding::Ascii,
            Some(array.data.len()),
            options,
        )?;
    }
    Ok(())
}

fn check_auxiliary_arrays(
    floats: &[DataArray<f32>],
    integers: &[DataArray<i32>],
    strings: &[DataArray<String>],
) -> Result<()> {
    let mut names = BTreeSet::new();
    for name in floats
        .iter()
        .map(|a| &a.name)
        .chain(integers.iter().map(|a| &a.name))
        .chain(strings.iter().map(|a| &a.name))
    {
        xml_string(name)?;
        if name.is_empty() || !names.insert(name) {
            return Err(Error::InvalidValue(
                "empty or duplicate auxiliary array name".into(),
            ));
        }
    }
    if floats.iter().flat_map(|a| &a.data).any(|v| !v.is_finite()) {
        return Err(Error::InvalidValue(
            "nonfinite auxiliary float value".into(),
        ));
    }
    for value in strings.iter().flat_map(|a| &a.data) {
        if !value.is_ascii() || value.as_bytes().contains(&0) {
            return Err(Error::Unsupported(
                "mzML string arrays require ASCII without embedded NUL".into(),
            ));
        }
    }
    Ok(())
}

/// Write plain mzML 1.1 XML with uncompressed binary arrays.
pub fn write(writer: impl Write, experiment: &MSExperiment) -> Result<()> {
    write_with_options(writer, experiment, &WriteOptions::default())
}

/// Write the supported data model after preflight validation.
/// Named float, integer and ASCII string arrays are preserved. Unsupported
/// metadata or unrepresentable array values are rejected before output.
pub fn write_with_options(
    mut w: impl Write,
    experiment: &MSExperiment,
    options: &WriteOptions,
) -> Result<()> {
    experiment.validate()?;
    let check_metadata = |metadata: &BTreeMap<String, String>, name: &str| -> Result<()> {
        xml_string(name)?;
        for (key, value) in metadata {
            xml_string(key)?;
            xml_string(value)?;
            if key == NAME_KEY {
                return Err(Error::InvalidValue(format!(
                    "reserved metadata name {NAME_KEY}"
                )));
            }
        }
        Ok(())
    };
    check_metadata(&experiment.metadata, "")?;
    let mut spectrum_ids = BTreeSet::new();
    let mut chromatogram_ids = BTreeSet::new();
    for (i, s) in experiment.spectra.iter().enumerate() {
        for precursor in &s.precursors {
            precursor_metadata::validate_write(precursor)?;
        }
        if !s.peptide_identifications.is_empty() {
            return Err(Error::Unsupported(
                "mzML writer cannot store peptide identifications".into(),
            ));
        }
        check_metadata(&s.metadata, &s.name)?;
        let id = if s.native_id.is_empty() {
            format!("index={i}")
        } else {
            s.native_id.clone()
        };
        xml_string(&id)?;
        if id.split(' ').any(|part| {
            part.split_once('=').is_none_or(|(key, value)| {
                key.is_empty()
                    || value.is_empty()
                    || key.chars().any(char::is_whitespace)
                    || value.chars().any(char::is_whitespace)
            })
        }) {
            return Err(Error::InvalidValue(
                "spectrum native_id must have mzML key=value form".into(),
            ));
        }
        if !spectrum_ids.insert(id) {
            return Err(Error::InvalidValue("duplicate spectrum native_id".into()));
        }
        check_auxiliary_arrays(
            &s.float_data_arrays,
            &s.integer_data_arrays,
            &s.string_data_arrays,
        )?;
    }
    for (i, c) in experiment.chromatograms.iter().enumerate() {
        precursor_metadata::validate_write(&c.precursor)?;
        check_metadata(&c.metadata, &c.name)?;
        let id = if c.native_id.is_empty() {
            format!("chromatogram={i}")
        } else {
            c.native_id.clone()
        };
        xml_string(&id)?;
        if !chromatogram_ids.insert(id) {
            return Err(Error::InvalidValue(
                "duplicate chromatogram native_id".into(),
            ));
        }
        check_auxiliary_arrays(
            &c.float_data_arrays,
            &c.integer_data_arrays,
            &c.string_data_arrays,
        )?;
    }
    for precursor in experiment
        .spectra
        .iter()
        .flat_map(|s| &s.precursors)
        .chain(experiment.chromatograms.iter().map(|c| &c.precursor))
    {
        if precursor
            .spectrum_reference
            .as_ref()
            .is_some_and(|reference| !spectrum_ids.contains(reference))
        {
            return Err(invalid(
                "precursor spectrum reference does not name an output spectrum",
            ));
        }
    }
    writeln!(
        w,
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" version=\"1.1.0\">"
    )?;
    writeln!(
        w,
        "<cvList count=\"2\"><cv id=\"MS\" fullName=\"PSI-MS\" URI=\"https://purl.obolibrary.org/obo/ms.obo\"/><cv id=\"UO\" fullName=\"Unit Ontology\" URI=\"https://purl.obolibrary.org/obo/uo.obo\"/></cvList>\n<fileDescription><fileContent>"
    )?;
    cv(&mut w, "MS:1000294", "mass spectrum", "", "")?;
    writeln!(
        w,
        "</fileContent></fileDescription>\n<softwareList count=\"1\"><software id=\"openms_rust\" version=\"{}\">",
        env!("CARGO_PKG_VERSION")
    )?;
    cv(
        &mut w,
        "MS:1000799",
        "custom unreleased software tool",
        "OpenMS Rust",
        "",
    )?;
    writeln!(
        w,
        "</software></softwareList>\n<instrumentConfigurationList count=\"1\"><instrumentConfiguration id=\"unknown_instrument\"><userParam name=\"original instrument information unavailable\"/></instrumentConfiguration></instrumentConfigurationList>\n<dataProcessingList count=\"1\"><dataProcessing id=\"conversion\"><processingMethod order=\"0\" softwareRef=\"openms_rust\">"
    )?;
    cv(&mut w, "MS:1000544", "Conversion to mzML", "", "")?;
    writeln!(
        w,
        "</processingMethod></dataProcessing></dataProcessingList>\n<run id=\"run\" defaultInstrumentConfigurationRef=\"unknown_instrument\">"
    )?;
    user_params(&mut w, &experiment.metadata, "")?;
    if !experiment.spectra.is_empty() {
        writeln!(
            w,
            "<spectrumList count=\"{}\" defaultDataProcessingRef=\"conversion\">",
            experiment.spectra.len()
        )?;
        for (i, spectrum) in experiment.spectra.iter().enumerate() {
            let id = if spectrum.native_id.is_empty() {
                format!("index={i}")
            } else {
                spectrum.native_id.clone()
            };
            writeln!(
                w,
                "<spectrum id=\"{}\" index=\"{i}\" defaultArrayLength=\"{}\">",
                escape(&id),
                spectrum.len()
            )?;
            cv(
                &mut w,
                "MS:1000511",
                "ms level",
                &spectrum.ms_level.to_string(),
                "",
            )?;
            match spectrum.spectrum_type {
                SpectrumType::Centroid => cv(&mut w, "MS:1000127", "centroid spectrum", "", "")?,
                SpectrumType::Profile => cv(&mut w, "MS:1000128", "profile spectrum", "", "")?,
                SpectrumType::Unknown => {}
            }
            user_params(&mut w, &spectrum.metadata, &spectrum.name)?;
            if spectrum.rt != -1.0 {
                writeln!(w, "<scanList count=\"1\">")?;
                cv(&mut w, "MS:1000795", "no combination", "", "")?;
                writeln!(w, "<scan>")?;
                cv(
                    &mut w,
                    "MS:1000016",
                    "scan start time",
                    &spectrum.rt.to_string(),
                    SECOND,
                )?;
                writeln!(w, "</scan></scanList>")?;
            }
            if !spectrum.precursors.is_empty() {
                writeln!(w, "<precursorList count=\"{}\">", spectrum.precursors.len())?;
                for precursor in &spectrum.precursors {
                    write_precursor(&mut w, precursor)?;
                }
                writeln!(w, "</precursorList>")?;
            }
            writeln!(
                w,
                "<binaryDataArrayList count=\"{}\">",
                2 + spectrum.float_data_arrays.len()
                    + spectrum.integer_data_arrays.len()
                    + spectrum.string_data_arrays.len()
            )?;
            write_array(
                &mut w,
                spectrum
                    .peaks
                    .iter()
                    .flat_map(|p| p.mz.to_le_bytes())
                    .collect(),
                Kind::Mz,
                Encoding::Float64,
                None,
                options,
            )?;
            write_array(
                &mut w,
                spectrum
                    .peaks
                    .iter()
                    .flat_map(|p| p.intensity.to_le_bytes())
                    .collect(),
                Kind::Intensity,
                Encoding::Float32,
                None,
                options,
            )?;
            write_auxiliary_arrays(
                &mut w,
                &spectrum.float_data_arrays,
                &spectrum.integer_data_arrays,
                &spectrum.string_data_arrays,
                options,
            )?;
            writeln!(w, "</binaryDataArrayList></spectrum>")?;
        }
        writeln!(w, "</spectrumList>")?;
    }
    if !experiment.chromatograms.is_empty() {
        writeln!(
            w,
            "<chromatogramList count=\"{}\" defaultDataProcessingRef=\"conversion\">",
            experiment.chromatograms.len()
        )?;
        for (i, chromatogram) in experiment.chromatograms.iter().enumerate() {
            let id = if chromatogram.native_id.is_empty() {
                format!("chromatogram={i}")
            } else {
                chromatogram.native_id.clone()
            };
            writeln!(
                w,
                "<chromatogram id=\"{}\" index=\"{i}\" defaultArrayLength=\"{}\">",
                escape(&id),
                chromatogram.len()
            )?;
            user_params(&mut w, &chromatogram.metadata, &chromatogram.name)?;
            if chromatogram.precursor != Precursor::default() {
                write_precursor(&mut w, &chromatogram.precursor)?;
            }
            writeln!(
                w,
                "<binaryDataArrayList count=\"{}\">",
                2 + chromatogram.float_data_arrays.len()
                    + chromatogram.integer_data_arrays.len()
                    + chromatogram.string_data_arrays.len()
            )?;
            write_array(
                &mut w,
                chromatogram
                    .peaks
                    .iter()
                    .flat_map(|p| p.rt.to_le_bytes())
                    .collect(),
                Kind::Time,
                Encoding::Float64,
                None,
                options,
            )?;
            write_array(
                &mut w,
                chromatogram
                    .peaks
                    .iter()
                    .flat_map(|p| p.intensity.to_le_bytes())
                    .collect(),
                Kind::Intensity,
                Encoding::Float32,
                None,
                options,
            )?;
            write_auxiliary_arrays(
                &mut w,
                &chromatogram.float_data_arrays,
                &chromatogram.integer_data_arrays,
                &chromatogram.string_data_arrays,
                options,
            )?;
            writeln!(w, "</binaryDataArrayList></chromatogram>")?;
        }
        writeln!(w, "</chromatogramList>")?;
    }
    writeln!(w, "</run></mzML>")?;
    w.flush()?;
    Ok(())
}
