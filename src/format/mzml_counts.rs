// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source loadSize counting without constructing peak arrays or an experiment.

use super::{
    NS, Parameter, ParameterBudget, ReadOptions, attributes, invalid, parameter_id, required,
    xml_string,
};
use crate::data_structures::list::ListParse;
use crate::format::peak_options::PeakFileOptions;
use crate::{Error, Result};
use quick_xml::{NsReader, events::Event, name::ResolveResult};
use std::{
    collections::BTreeMap,
    io::{self, BufRead, Read},
    path::Path,
};

/// Source list counts, or source event-based filtered record counts.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MzMLCounts {
    pub spectra: usize,
    pub chromatograms: usize,
}

/// Maximum buffered markup, comment, CDATA or entity-reference event.
pub const MAX_COUNT_EVENT_BYTES: usize = 1024 * 1024;
/// Maximum XML nesting, including ignored record payload.
pub const MAX_COUNT_XML_DEPTH: usize = 256;
/// Maximum semantic descriptor, expanded-parameter and MS-level membership work.
pub const MAX_COUNT_WORK: usize = 50_000_000;

/// Read source-default declared counts; stop once both list counts are known.
pub fn read_size(reader: impl BufRead) -> Result<MzMLCounts> {
    read_size_with_options(reader, &PeakFileOptions::default(), &ReadOptions::default())
}
/// Count using source RT/MS-level/precursor event ordering. Binary data is never
/// decoded. XML, record and parameter limits apply; array limits and acquisition
/// materialization are inapplicable. See MZML_COUNTS_SUPPORT.md for source quirks.
pub fn read_size_with_options(
    input: impl BufRead,
    scientific: &PeakFileOptions,
    limits: &ReadOptions,
) -> Result<MzMLCounts> {
    parse(input, scientific, limits)
}
/// Count a plain/gzip/bzip2 file selected by content magic, with source defaults.
pub fn load_size(path: impl AsRef<Path>) -> Result<MzMLCounts> {
    load_size_with_options(path, &PeakFileOptions::default(), &ReadOptions::default())
}
/// Path form of read_size_with_options. A successful early stop does not drain
/// compressed trailers. ZIP containers remain unsupported by the shared reader.
pub fn load_size_with_options(
    path: impl AsRef<Path>,
    scientific: &PeakFileOptions,
    limits: &ReadOptions,
) -> Result<MzMLCounts> {
    read_size_with_options(
        crate::format::path_io::open(path.as_ref())?,
        scientific,
        limits,
    )
}

// quick-xml intentionally leaves these XML lexical constraints to callers.
// Scan borrowed event bytes before the shared helper allocates attribute maps.
fn checked_attributes(
    element: &quick_xml::events::BytesStart<'_>,
    decoder: quick_xml::encoding::Decoder,
    budget: &mut ParameterBudget,
) -> Result<BTreeMap<String, String>> {
    let mut quote = None;
    let mut after_quote = false;
    for &byte in element.attributes_raw() {
        if let Some(q) = quote {
            if byte == q {
                quote = None;
                after_quote = true;
            }
        } else {
            if after_quote && !matches!(byte, b' ' | b'\t' | b'\r' | b'\n') {
                return Err(invalid("missing whitespace between XML attributes"));
            }
            after_quote = false;
            if matches!(byte, b'\'' | b'"') {
                quote = Some(byte);
            }
        }
    }
    for attribute in element.attributes().with_checks(false) {
        let attribute = attribute.map_err(|e| invalid(e.to_string()))?;
        let name =
            std::str::from_utf8(attribute.key.as_ref()).map_err(|e| invalid(e.to_string()))?;
        let mut parts = name.split(':');
        parameter_id(parts.next().unwrap_or(""))?;
        if let Some(local) = parts.next() {
            parameter_id(local)?;
        }
        if parts.next().is_some() || attribute.value.contains(&b'<') {
            return Err(invalid("invalid XML attribute name or raw value"));
        }
    }
    attributes(element, decoder, Some(budget))
}

fn spend(work: &mut usize, amount: usize) -> Result<()> {
    *work = work
        .checked_sub(amount)
        .ok_or_else(|| invalid("mzML count work limit exceeded"))?;
    Ok(())
}
fn integer(value: &str) -> Result<i32> {
    i32::from_list_item(value)
}
fn number(value: &str) -> Result<f64> {
    let n = f64::from_list_item(value)?;
    if n.is_finite() {
        Ok(n)
    } else {
        Err(invalid("nonfinite mzML count value"))
    }
}
fn increment(value: &mut usize, maximum: usize) -> Result<()> {
    *value = value
        .checked_add(1)
        .filter(|v| *v <= maximum)
        .ok_or_else(|| invalid("mzML count exceeds record limit"))?;
    Ok(())
}

// A per-event quota prevents quick-xml allocating an entire oversized event.
// Direct text consumption uses the same total byte counter. The original
// BufRead may prefetch, but only bytes actually consumed are charged here.
struct Input<R> {
    inner: R,
    bytes: u64,
    event: u64,
    consumed: u64,
}
impl<R: BufRead> BufRead for Input<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        let bytes = self.inner.fill_buf()?;
        if bytes.is_empty() {
            return Ok(bytes);
        }
        let n = self.bytes.min(self.event).min(bytes.len() as u64) as usize;
        if n == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "mzML count XML/event byte limit exceeded",
            ));
        }
        Ok(&bytes[..n])
    }
    fn consume(&mut self, n: usize) {
        self.inner.consume(n);
        self.bytes -= n as u64;
        self.event -= n as u64;
        self.consumed += n as u64;
    }
}
impl<R: BufRead> Read for Input<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        let n = {
            let bytes = self.fill_buf()?;
            let n = bytes.len().min(out.len());
            out[..n].copy_from_slice(&bytes[..n]);
            n
        };
        self.consume(n);
        Ok(n)
    }
}

#[derive(Clone, Copy)]
enum Encoding {
    Utf8,
    Ascii,
    Latin1Ascii,
}
impl Encoding {
    fn check(self, ascii: bool) -> Result<()> {
        if ascii || matches!(self, Self::Utf8) {
            return Ok(());
        }
        if matches!(self, Self::Latin1Ascii) {
            Err(Error::Unsupported(
                "non-ASCII ISO-8859-1 counting requires transcoding".into(),
            ))
        } else {
            Err(invalid("non-ASCII bytes in US-ASCII XML"))
        }
    }
}

// Constant-size UTF-8 decoder also checks XML Char, forbidden literal ']]>',
// and document-level whitespace. It spans arbitrary BufRead chunk boundaries.
#[derive(Default)]
struct TextCheck {
    continuation: u8,
    scalar: u32,
    minimum: u32,
    brackets: u8,
}
impl TextCheck {
    fn byte(&mut self, b: u8, whitespace: bool, encoding: Encoding) -> Result<()> {
        encoding.check(b.is_ascii())?;
        let scalar = if self.continuation != 0 {
            if b & 0xc0 != 0x80 {
                return Err(invalid("invalid UTF-8 XML text"));
            }
            self.scalar = (self.scalar << 6) | u32::from(b & 0x3f);
            self.continuation -= 1;
            if self.continuation != 0 {
                return Ok(());
            }
            if self.scalar < self.minimum {
                return Err(invalid("overlong UTF-8 XML text"));
            }
            self.scalar
        } else {
            match b {
                0..=0x7f => u32::from(b),
                0xc2..=0xdf => {
                    self.continuation = 1;
                    self.scalar = u32::from(b & 0x1f);
                    self.minimum = 0x80;
                    return Ok(());
                }
                0xe0..=0xef => {
                    self.continuation = 2;
                    self.scalar = u32::from(b & 0x0f);
                    self.minimum = 0x800;
                    return Ok(());
                }
                0xf0..=0xf4 => {
                    self.continuation = 3;
                    self.scalar = u32::from(b & 0x07);
                    self.minimum = 0x10000;
                    return Ok(());
                }
                _ => return Err(invalid("invalid UTF-8 XML text")),
            }
        };
        let c = char::from_u32(scalar).ok_or_else(|| invalid("invalid UTF-8 XML scalar"))?;
        valid_char(c, whitespace)?;
        if c == '>' && self.brackets == 2 {
            return Err(invalid("forbidden ]]> in XML text"));
        }
        self.brackets = if c == ']' {
            (self.brackets + 1).min(2)
        } else {
            0
        };
        Ok(())
    }
    fn finish(&self) -> Result<()> {
        if self.continuation == 0 {
            Ok(())
        } else {
            Err(invalid("truncated UTF-8 XML text"))
        }
    }
}
fn valid_char(c: char, whitespace: bool) -> Result<()> {
    if whitespace && !matches!(c, ' ' | '\t' | '\r' | '\n') {
        return Err(invalid(
            "non-whitespace text outside XML root or inside parameter",
        ));
    }
    if matches!(c, '\t' | '\n' | '\r')
        || ('\u{20}'..='\u{d7ff}').contains(&c)
        || ('\u{e000}'..='\u{fffd}').contains(&c)
        || c >= '\u{10000}'
    {
        Ok(())
    } else {
        Err(invalid("invalid XML 1.0 character"))
    }
}
fn discard_text<R: BufRead>(
    input: &mut Input<R>,
    whitespace: bool,
    encoding: Encoding,
) -> Result<()> {
    input.event = u64::MAX;
    let mut check = TextCheck::default();
    loop {
        let (n, stop) = {
            let bytes = input.fill_buf()?;
            let bytes = &bytes[..bytes.len().min(8192)];
            let n = bytes
                .iter()
                .position(|b| matches!(b, b'<' | b'&'))
                .unwrap_or(bytes.len());
            for &b in &bytes[..n] {
                check.byte(b, whitespace, encoding)?;
            }
            (n, n != bytes.len() || bytes.is_empty())
        };
        input.consume(n);
        if stop {
            return check.finish();
        }
    }
}

#[derive(Default)]
struct Record {
    chromatogram: bool,
    skip: bool,
    target: f64,
    selected_ions: usize,
}
struct State<'a> {
    scientific: &'a PeakFileOptions,
    limits: &'a ReadOptions,
    raw: bool,
    metadata_only: bool,
    declared: [i32; 2],
    counts: MzMLCounts,
    records: usize,
    skip_depth: Option<usize>,
    record: Option<Record>,
    groups: BTreeMap<String, Vec<Parameter>>,
    group: Option<(String, Vec<Parameter>)>,
    group_list: Option<(usize, usize)>,
    group_list_seen: bool,
    pending: Vec<String>,
    work: usize,
    budget: ParameterBudget,
    seen_mzml: bool,
    seen_run: bool,
    seen_lists: [bool; 2],
}
impl<'a> State<'a> {
    fn new(scientific: &'a PeakFileOptions, limits: &'a ReadOptions) -> Self {
        Self {
            scientific,
            limits,
            raw: !scientific.has_filters(),
            metadata_only: scientific.metadata_only,
            declared: [-1, -1],
            counts: MzMLCounts::default(),
            records: 0,
            skip_depth: None,
            record: None,
            groups: BTreeMap::new(),
            group: None,
            group_list: None,
            group_list_seen: false,
            pending: Vec::new(),
            work: MAX_COUNT_WORK,
            budget: ParameterBudget {
                remaining: limits.max_total_params,
                bytes: limits.max_param_bytes,
            },
            seen_mzml: false,
            seen_run: false,
            seen_lists: [false; 2],
        }
    }
    fn result(&self) -> MzMLCounts {
        if self.raw {
            MzMLCounts {
                spectra: self.declared[0].max(0) as usize,
                chromatograms: if self.scientific.skip_chromatograms {
                    0
                } else {
                    self.declared[1].max(0) as usize
                },
            }
        } else {
            self.counts
        }
    }
    fn skipped(&self) -> bool {
        self.skip_depth.is_some() || self.record.as_ref().is_some_and(|r| r.skip)
    }
    fn start(
        &mut self,
        tag: &str,
        attrs: BTreeMap<String, String>,
        stack: &[String],
    ) -> Result<bool> {
        let parent = stack.last().map(String::as_str).unwrap_or("");
        let grandparent = stack.iter().rev().nth(1).map(String::as_str).unwrap_or("");
        // Physical records remain bounded even when their body is never consumed.
        if matches!(tag, "spectrum" | "chromatogram") {
            increment(&mut self.records, self.limits.max_records)?;
            if (tag == "spectrum" && parent != "spectrumList")
                || (tag == "chromatogram" && parent != "chromatogramList")
            {
                return Err(invalid("misplaced mzML count record"));
            }
        }
        if self.skipped() {
            return Ok(false);
        }
        match tag {
            "indexedmzML" if !stack.is_empty() => return Err(invalid("misplaced indexedmzML")),
            "mzML" => {
                if self.seen_mzml
                    || !(stack.is_empty() || (stack.len() == 1 && parent == "indexedmzML"))
                {
                    return Err(invalid("misplaced or duplicate mzML root"));
                }
                if required(&attrs, "version")? != "1.1.0" {
                    return Err(Error::Unsupported("counting supports mzML 1.1.0".into()));
                }
                self.seen_mzml = true;
            }
            "run" => {
                if parent != "mzML" || self.seen_run {
                    return Err(invalid("misplaced or duplicate mzML run"));
                }
                self.seen_run = true;
                for id in &self.pending {
                    spend(&mut self.work, id.len() + 1)?;
                    if !self.groups.contains_key(id) {
                        return Err(invalid("unresolved header parameter reference"));
                    }
                }
                self.pending.clear();
            }
            "spectrumList" | "chromatogramList" => {
                if parent != "run" {
                    return Err(invalid("record list outside mzML run"));
                }
                let i = usize::from(tag == "chromatogramList");
                if self.seen_lists[i] {
                    return Err(invalid("duplicate mzML record list"));
                }
                self.seen_lists[i] = true;
                required(&attrs, "defaultDataProcessingRef")?;
                if self.metadata_only {
                    return Ok(true);
                }
                self.declared[i] = integer(required(&attrs, "count")?)?;
                let declared = self
                    .declared
                    .iter()
                    .map(|n| (*n).max(0) as usize)
                    .sum::<usize>();
                if declared > self.limits.max_records {
                    return Err(invalid("declared mzML counts exceed record limit"));
                }
                if self.raw {
                    if self.declared[1 - i] != -1 {
                        return Ok(true);
                    }
                    self.skip_depth = Some(stack.len() + 1);
                } else if i == 1 && self.scientific.skip_chromatograms {
                    self.skip_depth = Some(stack.len() + 1);
                }
            }
            "spectrum" | "chromatogram" => {
                required(&attrs, "id")?;
                integer(required(&attrs, "defaultArrayLength")?)?;
                let chromatogram = tag == "chromatogram";
                if chromatogram {
                    increment(&mut self.counts.chromatograms, self.limits.max_records)?;
                }
                self.record = Some(Record {
                    chromatogram,
                    skip: chromatogram,
                    ..Record::default()
                });
            }
            "referenceableParamGroupList" => {
                if parent != "mzML" || self.group_list_seen || self.seen_run {
                    return Err(invalid("misplaced or duplicate parameter group list"));
                }
                let count = integer(required(&attrs, "count")?)?;
                let count = usize::try_from(count)
                    .map_err(|_| invalid("negative parameter group count"))?;
                if count > self.limits.max_param_groups {
                    return Err(invalid("parameter group limit exceeded"));
                }
                self.group_list_seen = true;
                self.group_list = Some((count, 0));
            }
            "referenceableParamGroup" => {
                if parent != "referenceableParamGroupList" || self.group.is_some() {
                    return Err(invalid("misplaced parameter group"));
                }
                let (_, actual) = self
                    .group_list
                    .as_mut()
                    .ok_or_else(|| invalid("missing parameter group list"))?;
                increment(actual, self.limits.max_param_groups)?;
                let id = parameter_id(required(&attrs, "id")?)?;
                if self.groups.contains_key(id) {
                    return Err(invalid("duplicate parameter group ID"));
                }
                self.group = Some((id.to_owned(), Vec::new()));
            }
            "cvParam" | "userParam" => {
                if matches!(
                    parent,
                    "cvParam" | "userParam" | "referenceableParamGroupRef"
                ) {
                    return Err(invalid("nested XML parameter"));
                }
                if tag == "cvParam" {
                    required(&attrs, "accession")?;
                    required(&attrs, "name")?;
                }
                if parent == "referenceableParamGroup" {
                    self.group
                        .as_mut()
                        .ok_or_else(|| invalid("missing parameter group"))?
                        .1
                        .push(Parameter {
                            tag: if tag == "cvParam" {
                                "cvParam"
                            } else {
                                "userParam"
                            },
                            attrs,
                        });
                } else if tag == "cvParam" {
                    apply_cv(
                        &attrs,
                        parent,
                        grandparent,
                        self.scientific,
                        &mut self.record,
                        &mut self.counts,
                        self.limits.max_records,
                        &mut self.work,
                    )?;
                }
            }
            "referenceableParamGroupRef" => {
                if !super::parameter_context(parent) {
                    return Err(invalid("parameter reference outside supported XML context"));
                }
                let id = parameter_id(required(&attrs, "ref")?)?;
                if !self.seen_run {
                    self.pending.push(id.to_owned());
                } else {
                    let parameters = self
                        .groups
                        .get(id)
                        .ok_or_else(|| invalid("unknown parameter group reference"))?;
                    spend(&mut self.work, parameters.len() + id.len())?;
                    // Source deliberately does not check the skip flag between
                    // terms in one reference expansion. Multiple RTs can count.
                    for parameter in parameters {
                        self.budget.charge(&parameter.attrs)?;
                        if parameter.tag == "cvParam" {
                            apply_cv(
                                &parameter.attrs,
                                parent,
                                grandparent,
                                self.scientific,
                                &mut self.record,
                                &mut self.counts,
                                self.limits.max_records,
                                &mut self.work,
                            )?;
                        }
                    }
                }
            }
            "precursor" if parent == "precursorList" => {
                if let Some(r) = &mut self.record {
                    r.target = 0.0;
                    r.selected_ions = 0;
                }
            }
            "selectedIon" if parent == "selectedIonList" => {
                if let Some(r) = &mut self.record {
                    r.selected_ions = r
                        .selected_ions
                        .checked_add(1)
                        .ok_or_else(|| invalid("selected ion count overflow"))?;
                }
            }
            _ => {}
        }
        if self
            .counts
            .spectra
            .checked_add(self.counts.chromatograms)
            .is_none_or(|n| n > self.limits.max_records)
        {
            return Err(invalid("mzML counts exceed record limit"));
        }
        Ok(false)
    }
    fn end(&mut self, tag: &str, depth: usize) -> Result<()> {
        if self.skip_depth == Some(depth) {
            self.skip_depth = None;
        }
        if matches!(tag, "spectrum" | "chromatogram") {
            self.record = None;
        }
        if self.skipped() {
            return Ok(());
        }
        if tag == "referenceableParamGroup" {
            let (id, parameters) = self
                .group
                .take()
                .ok_or_else(|| invalid("missing parameter group"))?;
            self.groups.insert(id, parameters);
        } else if tag == "referenceableParamGroupList" {
            let (declared, actual) = self
                .group_list
                .take()
                .ok_or_else(|| invalid("missing group list"))?;
            if declared != actual {
                return Err(invalid("parameter group count mismatch"));
            }
        }
        Ok(())
    }
}

// Keep borrowed group storage disjoint from mutable count state.
#[allow(clippy::too_many_arguments)]
fn apply_cv(
    attrs: &BTreeMap<String, String>,
    parent: &str,
    grandparent: &str,
    scientific: &PeakFileOptions,
    record: &mut Option<Record>,
    counts: &mut MzMLCounts,
    max_records: usize,
    work: &mut usize,
) -> Result<()> {
    spend(work, 1)?;
    let Some(record) = record.as_mut() else {
        return Ok(());
    };
    if record.chromatogram {
        return Ok(());
    }
    let accession = required(attrs, "accession")?;
    let value = attrs.get("value").map(String::as_str).unwrap_or("");
    match (parent, accession) {
        ("spectrum", "MS:1000511") => {
            let level = integer(value)?;
            spend(work, scientific.ms_levels().len())?;
            if scientific.has_ms_levels() && !scientific.contains_ms_level(level) {
                record.skip = true;
            }
        }
        ("scan", "MS:1000016") => {
            let mut rt = number(value)?;
            if attrs
                .get("unitAccession")
                .is_some_and(|u| u == "UO:0000031")
            {
                rt *= 60.0;
            }
            if !rt.is_finite() {
                return Err(invalid("scan time overflows seconds"));
            }
            record.skip = true;
            if !scientific.has_rt_range() || super::load::contains(scientific.rt_range(), rt) {
                increment(&mut counts.spectra, max_records)?;
            }
        }
        ("isolationWindow", "MS:1000827") if grandparent == "precursor" => {
            record.target = number(value)?;
            if !scientific.precursor_mz_selected_ion
                && scientific.has_precursor_mz_range()
                && !super::load::contains(scientific.precursor_mz_range(), record.target)
            {
                record.skip = true;
            }
        }
        ("selectedIon", "MS:1000744") if record.selected_ions <= 1 => {
            let mz = number(value)?;
            if mz != record.target && scientific.precursor_mz_selected_ion {
                record.target = mz;
                if scientific.has_precursor_mz_range()
                    && !super::load::contains(scientific.precursor_mz_range(), mz)
                {
                    record.skip = true;
                }
            }
        }
        _ => {}
    }
    Ok(())
}

// The transform setup pass uses the same counting lexer/state and the ordinary
// header decoder. It does not construct records or decode the skipped arrays.
#[derive(Default)]
struct Setup {
    experiment: crate::MSExperiment,
    draft: super::header::Draft,
    registry: super::header::Registry,
    work: super::header::Work,
}
impl Setup {
    fn start(
        &mut self,
        tag: &str,
        attrs: BTreeMap<String, String>,
        stack: &[String],
        state: &mut State<'_>,
    ) -> Result<Option<BTreeMap<String, String>>> {
        if state.skipped() {
            return Ok(Some(attrs));
        }
        let parent = stack.last().map_or("", String::as_str);
        if self.draft.captures(tag, parent) {
            if state.seen_run {
                return Err(invalid("header element after run start"));
            }
            self.draft.start(tag, attrs, &mut self.work)?;
            return Ok(None);
        }
        match tag {
            "mzML" => {
                if let Some(accession) = attrs.get("accession") {
                    self.experiment.settings.document.identifier = self.work.copy(accession)?;
                }
                if let Some(id) = attrs.get("id") {
                    self.work
                        .meter()
                        .tree::<(String, crate::metadata::MetaValue)>(1)?;
                    self.experiment
                        .settings
                        .metadata
                        .insert("mzml_id".into(), self.work.copy(id)?.into());
                }
            }
            "run" => {
                if parent != "mzML" || state.seen_run {
                    return Err(invalid("misplaced or duplicate mzML run"));
                }
                self.registry = std::mem::take(&mut self.draft).finish(
                    &attrs,
                    &state.groups,
                    &mut state.budget,
                    &mut self.work,
                    &mut self.experiment.settings,
                )?;
            }
            "spectrumList" | "chromatogramList" => {
                self.registry.processing(
                    required(&attrs, "defaultDataProcessingRef")?,
                    &mut self.work,
                )?;
            }
            "cvParam" | "userParam" if parent == "run" => {
                super::apply_parameter(
                    tag,
                    &attrs,
                    parent,
                    &mut None,
                    &mut None,
                    &mut self.experiment,
                    None,
                )?;
                return Ok(None);
            }
            "referenceableParamGroupRef" if parent == "run" => {
                let id = parameter_id(required(&attrs, "ref")?)?;
                spend(&mut state.work, id.len() + 1)?;
                let parameters = state
                    .groups
                    .get(id)
                    .ok_or_else(|| invalid("unknown parameter group reference"))?;
                super::apply_group(
                    parameters,
                    parent,
                    &mut state.budget,
                    &mut None,
                    &mut None,
                    &mut self.experiment,
                    None,
                )?;
                return Ok(None);
            }
            "fileDescription"
            | "sourceFileList"
            | "sourceFile"
            | "contact"
            | "fileContent"
            | "sampleList"
            | "sample"
            | "softwareList"
            | "software"
            | "instrumentConfigurationList"
            | "instrumentConfiguration"
            | "componentList"
            | "source"
            | "analyzer"
            | "detector"
            | "softwareRef"
            | "dataProcessingList"
            | "dataProcessing"
            | "processingMethod" => {
                return Err(invalid("misplaced mzML header element"));
            }
            _ => {}
        }
        Ok(Some(attrs))
    }
}
pub(super) fn setup(
    input: impl BufRead,
    scientific: &PeakFileOptions,
    limits: &ReadOptions,
    metadata_only: bool,
) -> Result<(crate::metadata::ExperimentalSettings, MzMLCounts)> {
    let mut setup = Setup::default();
    let counts = parse_impl(input, scientific, limits, Some(&mut setup), metadata_only)?;
    Ok((setup.experiment.settings, counts))
}

fn parse(
    input: impl BufRead,
    scientific: &PeakFileOptions,
    limits: &ReadOptions,
) -> Result<MzMLCounts> {
    parse_impl(input, scientific, limits, None, scientific.metadata_only)
}

fn parse_impl(
    input: impl BufRead,
    scientific: &PeakFileOptions,
    limits: &ReadOptions,
    mut setup: Option<&mut Setup>,
    metadata_only: bool,
) -> Result<MzMLCounts> {
    let mut state = State::new(scientific, limits);
    if setup.is_some() {
        state.raw = true;
    }
    state.metadata_only = metadata_only;
    // quick-xml only strips a BOM present in one fill_buf result. Replay a
    // bounded prefix so a BOM split across even one-byte readers is recognized.
    let mut input = input;
    let mut prefix = [0u8; 3];
    let mut length = 0;
    let maximum = limits.max_xml_bytes.min(3) as usize;
    while length < maximum {
        match input.read(&mut prefix[length..maximum]) {
            Ok(0) => break,
            Ok(n) => length += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
    }
    let start = if length == 3 && prefix == [0xef, 0xbb, 0xbf] {
        3
    } else {
        0
    };
    let replay = io::Cursor::new(&prefix[start..length]).chain(input);
    let mut reader = NsReader::from_reader(Input {
        inner: replay,
        bytes: limits.max_xml_bytes - start as u64,
        event: u64::MAX,
        consumed: start as u64,
    });
    reader.config_mut().enable_all_checks(true);
    // Empty events are handled locally. This keeps the next parser state at
    // text, so bounded direct text reads never bypass a pending synthetic end.
    reader.config_mut().expand_empty_elements = false;
    let mut buffer = Vec::new();
    let mut stack = Vec::<String>::new();
    let mut seen_root = false;
    let mut declaration = false;
    let mut encoding = Encoding::Utf8;
    let mut can_discard = true;
    loop {
        let whitespace = stack.is_empty()
            || (setup.is_some() && !state.seen_run)
            || (!state.skipped()
                && stack.last().is_some_and(|s| {
                    matches!(
                        s.as_str(),
                        "referenceableParamGroupList"
                            | "referenceableParamGroup"
                            | "referenceableParamGroupRef"
                            | "cvParam"
                            | "userParam"
                    )
                }));
        if can_discard {
            discard_text(reader.get_mut(), whitespace, encoding)?;
        }
        // Reserve a conservative temporary allowance before quick-xml can copy
        // bytes or namespace slots; return unused allowance after the event.
        let quota = (state.budget.bytes.saturating_sub(1024) / 16).min(MAX_COUNT_EVENT_BYTES);
        if quota == 0 {
            return Err(invalid("mzML count event storage limit exceeded"));
        }
        let allowance = quota * 16 + 1024;
        state.budget.spend(allowance)?;
        reader.get_mut().event = quota as u64;
        let before = reader.get_mut().consumed;
        let decoder = reader.decoder();
        let (namespace, event) = reader
            .read_resolved_event_into(&mut buffer)
            .map_err(|e| invalid(e.to_string()))?;
        let namespace_ok = matches!(namespace,ResolveResult::Bound(ns) if ns.as_ref()==NS);
        let used = (reader.get_mut().consumed - before) as usize;
        state.budget.bytes += allowance - (used * 16 + 1024);
        spend(&mut state.work, 1)?;
        encoding.check(event.is_ascii())?;
        can_discard = !matches!(event, Event::Text(_));
        let empty = matches!(&event, Event::Empty(_));
        match event {
            Event::Start(element) | Event::Empty(element) => {
                if !namespace_ok {
                    return Err(invalid("element outside mzML namespace"));
                }
                if stack.len() >= MAX_COUNT_XML_DEPTH {
                    return Err(invalid("mzML count XML depth limit exceeded"));
                }
                let local = element.local_name();
                let tag =
                    std::str::from_utf8(local.as_ref()).map_err(|e| invalid(e.to_string()))?;
                parameter_id(tag)?;
                if stack.is_empty() {
                    if seen_root || !matches!(tag, "mzML" | "indexedmzML") {
                        return Err(invalid("invalid or duplicate mzML document root"));
                    }
                    seen_root = true;
                }
                let attrs = checked_attributes(&element, decoder, &mut state.budget)?;
                let attrs = if let Some(setup) = setup.as_deref_mut() {
                    setup.start(tag, attrs, &stack, &mut state)?
                } else {
                    Some(attrs)
                };
                if let Some(attrs) = attrs {
                    if state.start(tag, attrs, &stack)? {
                        return Ok(state.result());
                    }
                }
                stack.push(tag.to_owned());
                if empty {
                    if let Some(setup) = setup.as_deref_mut() {
                        setup.draft.end(tag)?;
                    }
                    state.end(tag, stack.len())?;
                    stack.pop();
                }
            }
            Event::End(element) => {
                if !namespace_ok {
                    return Err(invalid("element outside mzML namespace"));
                }
                let tag = std::str::from_utf8(element.local_name().as_ref())
                    .map_err(|e| invalid(e.to_string()))?
                    .to_owned();
                if let Some(setup) = setup.as_deref_mut() {
                    setup.draft.end(&tag)?;
                }
                state.end(&tag, stack.len())?;
                if stack.pop().as_deref() != Some(tag.as_str()) {
                    return Err(invalid("mismatched XML end element"));
                }
            }
            Event::Text(text) => {
                let mut check = TextCheck::default();
                for &b in text.as_ref() {
                    check.byte(b, whitespace, encoding)?;
                }
                check.finish()?;
            }
            Event::GeneralRef(reference) => {
                if stack.is_empty() {
                    return Err(invalid("XML reference outside root"));
                }
                let c = if let Some(c) = reference
                    .resolve_char_ref()
                    .map_err(|e| invalid(e.to_string()))?
                {
                    c
                } else {
                    match reference.as_ref() {
                        b"lt" => '<',
                        b"gt" => '>',
                        b"amp" => '&',
                        b"apos" => '\'',
                        b"quot" => '"',
                        _ => return Err(invalid("unknown XML text entity")),
                    }
                };
                valid_char(c, whitespace)?;
            }
            Event::CData(text) => {
                if stack.is_empty() {
                    return Err(invalid("CDATA outside XML root"));
                }
                let text = text.decode().map_err(|e| invalid(e.to_string()))?;
                for c in text.chars() {
                    valid_char(c, whitespace)?;
                }
            }
            Event::Decl(d) => {
                if seen_root || declaration || before != start as u64 {
                    return Err(invalid("misplaced XML declaration"));
                }
                declaration = true;
                let raw = std::str::from_utf8(d.as_ref()).map_err(|e| invalid(e.to_string()))?;
                let decl = quick_xml::events::BytesStart::from_content(raw, 3);
                let attrs = checked_attributes(&decl, decoder, &mut state.budget)?;
                if attrs
                    .keys()
                    .any(|k| !matches!(k.as_str(), "version" | "encoding" | "standalone"))
                    || attrs
                        .get("standalone")
                        .is_some_and(|v| !matches!(v.as_str(), "yes" | "no"))
                {
                    return Err(invalid("invalid XML declaration attribute"));
                }
                // XMLDecl has an ordered grammar, unlike ordinary attributes.
                let mut previous = 0;
                for attribute in decl.attributes().with_checks(false) {
                    let attribute = attribute.map_err(|e| invalid(e.to_string()))?;
                    let position = match attribute.key.as_ref() {
                        b"version" => 0,
                        b"encoding" => 1,
                        _ => 2,
                    };
                    if position < previous {
                        return Err(invalid("XML declaration attributes out of order"));
                    }
                    previous = position;
                }
                if d.version().map_err(|e| invalid(e.to_string()))?.as_ref() != b"1.0" {
                    return Err(Error::Unsupported("only XML 1.0 is supported".into()));
                }
                if let Some(label) = d.encoding() {
                    let label = label.map_err(|e| invalid(e.to_string()))?;
                    // The pinned literal source fixture uses this declaration
                    // with entirely ASCII bytes. Non-ASCII Latin-1 is checked,
                    // never interpreted as UTF-8 or silently transcoded.
                    encoding = if label.eq_ignore_ascii_case(b"US-ASCII") {
                        Encoding::Ascii
                    } else if label.eq_ignore_ascii_case(b"ISO-8859-1") {
                        Encoding::Latin1Ascii
                    } else if label.eq_ignore_ascii_case(b"UTF-8") {
                        Encoding::Utf8
                    } else {
                        return Err(Error::Unsupported("only UTF-8 or ASCII-compatible US-ASCII/ISO-8859-1 counting is supported".into()));
                    };
                }
            }
            Event::Comment(text) => {
                xml_string(&text.decode().map_err(|e| invalid(e.to_string()))?)?;
            }
            Event::PI(text) => {
                let target =
                    std::str::from_utf8(text.target()).map_err(|e| invalid(e.to_string()))?;
                if target.eq_ignore_ascii_case("xml") {
                    return Err(invalid("reserved XML PI target"));
                }
                parameter_id(&target.replace(':', "_"))?;
                xml_string(
                    std::str::from_utf8(text.as_ref()).map_err(|e| invalid(e.to_string()))?,
                )?;
            }
            Event::DocType(_) => {
                return Err(Error::Unsupported("XML DTDs are not supported".into()));
            }
            Event::Eof => break,
        }
        buffer.clear();
    }
    if !seen_root || !state.seen_mzml || !state.seen_run || !stack.is_empty() {
        return Err(invalid("incomplete mzML count document"));
    }
    Ok(state.result())
}
