// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Source writer flags, bounded two-pass markup and indexed UTF-8 output.
use super::*;
use crate::format::{numpress_coder::NumpressCoderLimits, peak_options::PeakFileOptions};
use crate::kernel::data_array::Meter;
use sha1::{Digest, Sha1};
use std::{borrow::Cow, io, path::Path};

/// Binary preparation and cumulative markup/index allowances are independent of
/// the existing fixed header allowance. Both markup passes are precharged.
#[derive(Clone, Copy, Debug)]
pub struct PeakWriteLimits {
    pub binary: NumpressCoderLimits,
    pub max_xml_bytes: u64,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for PeakWriteLimits {
    fn default() -> Self {
        Self {
            binary: NumpressCoderLimits::default(),
            max_xml_bytes: 512 * 1024 * 1024,
            max_work: 2_000_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}
/// Encoding outcomes and the exact uncompressed output length.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PeakWriteReport {
    pub binary: NumpressWriteReport,
    pub indexed: bool,
    pub xml_bytes: u64,
}
pub fn write_with_peak_options(
    writer: impl Write,
    experiment: &MSExperiment,
    options: &PeakFileOptions,
) -> Result<PeakWriteReport> {
    write_with_peak_options_and_limits(writer, experiment, options, &PeakWriteLimits::default())
}
/// Read-only flags do not filter, sort or copy input on store. Indexed empty
/// output is rejected before writing (CPP-050); SHA-1 covers actual UTF-8 bytes
/// through the opening checksum tag (CPP-049). Stream I/O is not transactional.
pub fn write_with_peak_options_and_limits(
    writer: impl Write,
    experiment: &MSExperiment,
    options: &PeakFileOptions,
    limits: &PeakWriteLimits,
) -> Result<PeakWriteReport> {
    let prepared = prepare(experiment, options, limits)?;
    emit(writer, experiment, &prepared)
}
pub fn store_with_peak_options(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &PeakFileOptions,
) -> Result<PeakWriteReport> {
    store_with_peak_options_and_limits(path, experiment, options, &PeakWriteLimits::default())
}
/// Prepares before opening any temporary output; the shared path transport
/// publishes atomically after compression/flush. Offsets address decoded XML.
pub fn store_with_peak_options_and_limits(
    path: impl AsRef<Path>,
    experiment: &MSExperiment,
    options: &PeakFileOptions,
    limits: &PeakWriteLimits,
) -> Result<PeakWriteReport> {
    let prepared = prepare(experiment, options, limits)?;
    let mut report = None;
    crate::format::path_io::write(path.as_ref(), |writer| {
        report = Some(emit(writer, experiment, &prepared)?);
        Ok(())
    })?;
    report.ok_or_else(resource)
}
fn resource() -> Error {
    Error::InvalidValue("mzML writer execution resource limit exceeded".into())
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(resource)
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(resource)
}
struct Work {
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn meter(&mut self) -> Meter<'_> {
        Meter {
            work: &mut self.remaining,
            bytes: &mut self.bytes,
        }
    }
    fn charge(&mut self, n: usize, bytes: usize) -> Result<()> {
        self.meter().charge(n, bytes)
    }
    fn vector(&mut self, n: usize) -> Result<Vec<u64>> {
        self.meter().slots::<u64>(n)?;
        let mut v = Vec::new();
        v.try_reserve_exact(n).map_err(|_| resource())?;
        Ok(v)
    }
    fn precursor(&mut self, p: &Precursor) -> Result<()> {
        let mut m = self.meter();
        m.slots::<Precursor>(1)?;
        m.tree::<crate::metadata::ActivationMethod>(p.activation_methods.len())?;
        m.slots::<i32>(p.possible_charge_states.len())?;
        if let Some(id) = &p.spectrum_reference {
            m.text(id)?;
        }
        m.cv(&p.cv_terms)
    }
    fn record(
        &mut self,
        id: &str,
        name: &str,
        meta: &BTreeMap<String, String>,
        floats: &[DataArray<f32>],
        integers: &[DataArray<i32>],
        strings: &[DataArray<String>],
    ) -> Result<()> {
        let mut m = self.meter();
        m.text(id)?;
        m.text(name)?;
        m.tree::<(String, String)>(meta.len())?;
        for (key, value) in meta {
            m.text(key)?;
            m.text(value)?;
        }
        m.slots::<String>(add(floats.len(), add(integers.len(), strings.len())?)?)?;
        for name in floats
            .iter()
            .map(|v| &v.name)
            .chain(integers.iter().map(|v| &v.name))
            .chain(strings.iter().map(|v| &v.name))
        {
            m.text(name)?;
        }
        Ok(())
    }
    fn emission_inputs(&mut self, experiment: &MSExperiment) -> Result<()> {
        let initial_work = self.remaining;
        let initial_bytes = self.bytes;
        self.meter().slots::<MSSpectrum>(experiment.spectra.len())?;
        self.meter()
            .slots::<MSChromatogram>(experiment.chromatograms.len())?;
        for s in &experiment.spectra {
            self.record(
                &s.native_id,
                &s.name,
                &s.metadata,
                &s.float_data_arrays,
                &s.integer_data_arrays,
                &s.string_data_arrays,
            )?;
            s.acquisition_with_budget(&mut self.remaining, &mut self.bytes)?;
            self.meter().slots::<Precursor>(s.precursors.len())?;
            for p in &s.precursors {
                self.precursor(p)?;
            }
        }
        for c in &experiment.chromatograms {
            self.record(
                &c.native_id,
                &c.name,
                &c.metadata,
                &c.float_data_arrays,
                &c.integer_data_arrays,
                &c.string_data_arrays,
            )?;
            c.acquisition_with_budget(&mut self.remaining, &mut self.bytes)?;
            self.precursor(&c.precursor)?;
            self.meter().slots::<Product>(1)?;
            self.meter().cv(&c.product.cv_terms)?;
        }
        let work = initial_work - self.remaining;
        let bytes = initial_bytes - self.bytes;
        // Before either pass, cover repeated metadata walks, escaped strings,
        // default IDs and decimal scalar scratch. Metered struct/tree storage
        // gives conservative headroom for fixed scalars (including f64 text).
        self.charge(add(mul(work, 15)?, mul(bytes, 4)?)?, mul(bytes, 15)?)?;
        let records = add(experiment.spectra.len(), experiment.chromatograms.len())?;
        self.charge(mul(records, 4096)?, mul(records, 4096)?)
    }
}
struct Layout {
    indexed: bool,
    spectra: Vec<u64>,
    chromatograms: Vec<u64>,
    index_list: u64,
    xml_bytes: u64,
}
struct Prepared {
    payload: numpress_transport::Prepared,
    binary: WriteOptions,
    tpp: bool,
    layout: Layout,
}
fn prepare(
    experiment: &MSExperiment,
    options: &PeakFileOptions,
    limits: &PeakWriteLimits,
) -> Result<Prepared> {
    if limits.max_xml_bytes > u64::MAX / 8 {
        return Err(resource());
    }
    if options.write_index && experiment.spectra.is_empty() && experiment.chromatograms.is_empty() {
        return Err(Error::InvalidValue(
            "indexed mzML requires a spectrum or chromatogram".into(),
        ));
    }
    let mut work = Work {
        remaining: limits.max_work,
        bytes: limits.max_bytes,
    };
    work.emission_inputs(experiment)?;
    let mut layout = Layout {
        indexed: options.write_index,
        spectra: work.vector(if options.write_index {
            experiment.spectra.len()
        } else {
            0
        })?,
        chromatograms: work.vector(if options.write_index {
            experiment.chromatograms.len()
        } else {
            0
        })?,
        index_list: 0,
        xml_bytes: 0,
    };
    // Only the source writer's fixed-size flags/configurations are consumed.
    // Read filters and the caller's potentially large MS-level Vec are untouched.
    let binary_options = NumpressWriteOptions {
        binary: WriteOptions {
            zlib_compression: options.zlib_compression,
        },
        mass_time: options.numpress_configuration_mass_time(),
        intensity: options.numpress_configuration_intensity(),
        float_data_array: options.numpress_configuration_float_data_array(),
        limits: limits.binary,
    };
    let mass_np = binary_options.mass_time.compression != NumpressCompression::None;
    let precision = |is32| {
        if is32 && !mass_np {
            Encoding::Float32
        } else {
            Encoding::Float64
        }
    };
    let payload = numpress_transport::prepare(
        experiment,
        &binary_options,
        precision(options.mz_32_bit),
        precision(options.intensity_32_bit),
    )?;
    let mut sink = Output {
        writer: io::sink(),
        position: 0,
        hash: None,
        mode: Mode::Measure {
            work: &mut work,
            layout: &mut layout,
            max_xml_bytes: limits.max_xml_bytes,
        },
    };
    write_document(
        &mut sink,
        experiment,
        &binary_options.binary,
        &mut Some(payload.arrays.iter()),
        &payload.header,
        options.force_tpp_compatibility,
    )?;
    let length = sink.position;
    layout.xml_bytes = length;
    Ok(Prepared {
        payload,
        binary: binary_options.binary,
        tpp: options.force_tpp_compatibility,
        layout,
    })
}
fn emit(writer: impl Write, experiment: &MSExperiment, p: &Prepared) -> Result<PeakWriteReport> {
    let mut output = Output {
        writer,
        position: 0,
        hash: p.layout.indexed.then(Sha1::new),
        mode: Mode::Emit(&p.layout),
    };
    write_document(
        &mut output,
        experiment,
        &p.binary,
        &mut Some(p.payload.arrays.iter()),
        &p.payload.header,
        p.tpp,
    )?;
    if output.position != p.layout.xml_bytes {
        return Err(invalid("prepared mzML output length changed"));
    }
    Ok(PeakWriteReport {
        binary: p.payload.report,
        indexed: p.layout.indexed,
        xml_bytes: output.position,
    })
}
pub(super) fn native_id(id: &str, index: usize, chrom: bool) -> Cow<'_, str> {
    if id.is_empty() {
        Cow::Owned(if chrom {
            format!("chromatogram={index}")
        } else {
            format!("index={index}")
        })
    } else {
        Cow::Borrowed(id)
    }
}
enum Mode<'a> {
    Legacy,
    Measure {
        work: &'a mut Work,
        layout: &'a mut Layout,
        max_xml_bytes: u64,
    },
    Emit(&'a Layout),
}
pub(super) struct Output<'a, W> {
    writer: W,
    position: u64,
    hash: Option<Sha1>,
    mode: Mode<'a>,
}
impl<W: Write> Output<'_, W> {
    pub(super) fn legacy(writer: W) -> Self {
        Self {
            writer,
            position: 0,
            hash: None,
            mode: Mode::Legacy,
        }
    }
    fn indexed(&self) -> bool {
        match &self.mode {
            Mode::Legacy => false,
            Mode::Measure { layout, .. } => layout.indexed,
            Mode::Emit(layout) => layout.indexed,
        }
    }
    pub(super) fn header(&mut self, prefix: &str) -> Result<()> {
        self.write_all(b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")?;
        if self.indexed() {
            self.write_all(b"<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd\">\n")?;
        }
        self.write_all(prefix.as_bytes())?;
        Ok(())
    }
    pub(super) fn record(&mut self, chrom: bool) -> Result<()> {
        if let Mode::Measure { layout, .. } = &mut self.mode {
            if layout.indexed {
                let offsets = if chrom {
                    &mut layout.chromatograms
                } else {
                    &mut layout.spectra
                };
                if offsets.len() == offsets.capacity() {
                    return Err(resource());
                }
                offsets.push(self.position);
            }
        }
        Ok(())
    }
    fn offset(&self, chrom: bool, i: usize) -> Result<u64> {
        let layout = match &self.mode {
            Mode::Measure { layout, .. } => &**layout,
            Mode::Emit(layout) => layout,
            Mode::Legacy => return Err(resource()),
        };
        let offsets = if chrom {
            &layout.chromatograms
        } else {
            &layout.spectra
        };
        offsets.get(i).copied().ok_or_else(resource)
    }
    pub(super) fn footer(&mut self, experiment: &MSExperiment) -> Result<()> {
        if !self.indexed() {
            return Ok(());
        }
        let index_list = self.position;
        if let Mode::Measure { layout, .. } = &mut self.mode {
            layout.index_list = index_list;
        }
        if let Mode::Emit(layout) = &self.mode {
            if index_list != layout.index_list {
                return Err(invalid("prepared index offset changed"));
            }
        }
        let count = usize::from(!experiment.spectra.is_empty())
            + usize::from(!experiment.chromatograms.is_empty());
        writeln!(self, "<indexList count=\"{count}\">")?;
        if !experiment.spectra.is_empty() {
            writeln!(self, "<index name=\"spectrum\">")?;
            for (i, s) in experiment.spectra.iter().enumerate() {
                let offset = self.offset(false, i)?;
                writeln!(
                    self,
                    "<offset idRef=\"{}\">{offset}</offset>",
                    escape(&native_id(&s.native_id, i, false))
                )?;
            }
            writeln!(self, "</index>")?;
        }
        if !experiment.chromatograms.is_empty() {
            writeln!(self, "<index name=\"chromatogram\">")?;
            for (i, c) in experiment.chromatograms.iter().enumerate() {
                let offset = self.offset(true, i)?;
                writeln!(
                    self,
                    "<offset idRef=\"{}\">{offset}</offset>",
                    escape(&native_id(&c.native_id, i, true))
                )?;
            }
            writeln!(self, "</index>")?;
        }
        writeln!(
            self,
            "</indexList>\n<indexListOffset>{index_list}</indexListOffset>"
        )?;
        self.write_all(b"<fileChecksum>")?;
        if let Some(hash) = self.hash.take() {
            for byte in hash.finalize() {
                write!(self, "{byte:02x}")?;
            }
        } else {
            self.write_all(b"0000000000000000000000000000000000000000")?;
        }
        self.write_all(b"</fileChecksum>\n</indexedmzML>\n")?;
        Ok(())
    }
}
impl<W: Write> Write for Output<'_, W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let error = || io::Error::other("mzML output resource limit exceeded");
        let full = self
            .position
            .checked_add(bytes.len() as u64)
            .ok_or_else(error)?;
        if let Mode::Measure {
            work,
            layout,
            max_xml_bytes,
        } = &mut self.mode
        {
            if full > *max_xml_bytes {
                return Err(error());
            }
            // Both markup traversals plus the incremental checksum byte work
            // are paid on the first pass, before any external output exists.
            let factor = if layout.indexed { 3 } else { 2 };
            work.charge(bytes.len().checked_mul(factor).ok_or_else(error)?, 0)
                .map_err(|_| error())?;
        }
        let n = self.writer.write(bytes)?;
        let written = bytes
            .get(..n)
            .ok_or_else(|| io::Error::other("invalid writer byte count"))?;
        self.position = self.position.checked_add(n as u64).ok_or_else(error)?;
        if let Some(hash) = &mut self.hash {
            hash.update(written);
        }
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}
