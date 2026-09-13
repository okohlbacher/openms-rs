// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Numpress configuration, rejection diagnostics and base64/zlib transport.
//! No mzML reader/writer is enabled or changed by these operations.

use super::numpress::{self as raw, NumpressLimits};
pub use super::peak_options::{NumpressCompression, NumpressConfig};
use crate::{Error, Result};
use base64::{Engine, engine::general_purpose::STANDARD};
use flate2::{Compress, Compression, Decompress, FlushCompress, FlushDecompress, Status};

/// Failed error check at the highest failing input index, as in the source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NumpressAccuracyFailure {
    pub index: usize,
    pub original: f64,
    pub decoded: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub enum NumpressRejection {
    /// The raw encoder, factor estimation or verification decoder rejected data.
    Codec(String),
    Accuracy(NumpressAccuracyFailure),
}
#[derive(Clone, Debug, PartialEq)]
pub enum NumpressEncodeStatus {
    Encoded,
    EmptyInput,
    Disabled,
    Rejected(NumpressRejection),
}
/// Encoding rejections are explicit, without source stderr side effects.
/// `output` is empty unless status is Encoded; rejection never means ordinary
/// binary encoding was substituted. `fixed_point` is absent for PIC/skip paths.
#[derive(Clone, Debug, PartialEq)]
pub struct NumpressEncodeReport<T> {
    pub output: T,
    pub status: NumpressEncodeStatus,
    pub fixed_point: Option<f64>,
    pub used_maximal_fixed_point_fallback: bool,
}
impl<T> NumpressEncodeReport<T> {
    pub fn is_encoded(&self) -> bool {
        self.status == NumpressEncodeStatus::Encoded
    }
}

/// Raw limits plus cumulative wrapper allocations and base64 text length.
/// `raw.max_work` is shared across promotion, estimation, raw codec/verification
/// passes, validation, base64 and optional zlib processing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NumpressCoderLimits {
    pub raw: NumpressLimits,
    pub max_text_bytes: usize,
    pub max_total_bytes: usize,
}
impl Default for NumpressCoderLimits {
    fn default() -> Self {
        Self {
            raw: NumpressLimits::default(),
            max_text_bytes: 192 * 1024 * 1024,
            max_total_bytes: 512 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MSNumpressCoder {
    pub limits: NumpressCoderLimits,
}

impl MSNumpressCoder {
    pub fn encode_raw(
        &self,
        input: &[f64],
        config: &NumpressConfig,
    ) -> Result<NumpressEncodeReport<Vec<u8>>> {
        encode_raw(input, config, &mut Work::new(self.limits))
    }
    /// Source raw-output semantics: skipped/rejected encoding preserves output.
    /// Native resource errors also preserve it.
    pub fn encode_raw_into(
        &self,
        input: &[f64],
        output: &mut Vec<u8>,
        config: &NumpressConfig,
    ) -> Result<NumpressEncodeStatus> {
        let report = self.encode_raw(input, config)?;
        if report.is_encoded() {
            *output = report.output;
        }
        Ok(report.status)
    }
    pub fn encode(
        &self,
        input: &[f64],
        zlib: bool,
        config: &NumpressConfig,
    ) -> Result<NumpressEncodeReport<String>> {
        encode_text(input, zlib, config, &mut Work::new(self.limits))
    }
    /// Source float overload promotes each f32 exactly to f64 before any codec.
    pub fn encode_f32(
        &self,
        input: &[f32],
        zlib: bool,
        config: &NumpressConfig,
    ) -> Result<NumpressEncodeReport<String>> {
        let mut work = Work::new(self.limits);
        work.value_count(input.len())?;
        work.spend(input.len())?;
        let mut promoted = work.vector(input.len())?;
        promoted.extend(input.iter().map(|&value| f64::from(value)));
        encode_text(&promoted, zlib, config, &mut work)
    }
    /// Successful skip/rejection clears the destination as source encodeNP does.
    /// Native resource/transport errors preserve it rather than clearing early.
    pub fn encode_into(
        &self,
        input: &[f64],
        output: &mut String,
        zlib: bool,
        config: &NumpressConfig,
    ) -> Result<NumpressEncodeStatus> {
        let report = self.encode(input, zlib, config)?;
        *output = report.output;
        Ok(report.status)
    }
    pub fn decode_raw(&self, input: &[u8], config: &NumpressConfig) -> Result<Vec<f64>> {
        decode_raw(input, config.compression, &mut Work::new(self.limits))
    }
    pub fn decode(&self, input: &str, zlib: bool, config: &NumpressConfig) -> Result<Vec<f64>> {
        decode_text(input, zlib, config, &mut Work::new(self.limits))
    }

    /// Native atomic output replacement; source decode clears output first.
    pub fn decode_into(
        &self,
        input: &str,
        output: &mut Vec<f64>,
        zlib: bool,
        config: &NumpressConfig,
    ) -> Result<()> {
        *output = self.decode(input, zlib, config)?;
        Ok(())
    }
}

pub(crate) fn decode_text(
    input: &str,
    zlib: bool,
    config: &NumpressConfig,
    work: &mut Work,
) -> Result<Vec<f64>> {
    work.text(input.len())?;
    // Source Base64::decodeSingleString returns without decoding short text,
    // even for zlib=true or nonalphabet bytes. The local raw result is empty.
    if input.len() < 4 {
        return Ok(Vec::new());
    }
    if input.len() % 4 != 0 {
        return Err(invalid("invalid base64 length"));
    }
    let padding = if input.ends_with("==") {
        2
    } else {
        usize::from(input.ends_with('='))
    };
    let count = mul(input.len() / 4, 3)?
        .checked_sub(padding)
        .ok_or_else(|| invalid("invalid base64 padding"))?;
    work.binary(count)?;
    work.spend(input.len())?;
    let mut binary = work.vector(count)?;
    binary.resize(count, 0);
    let length = STANDARD
        .decode_slice(input.as_bytes(), &mut binary)
        .map_err(|_| invalid("invalid base64 alphabet, length or padding"))?;
    binary.truncate(length);
    let raw = if zlib {
        zlib_decode(&binary, work)?
    } else {
        binary
    };
    decode_raw(&raw, config.compression, work)
}

fn report(status: NumpressEncodeStatus) -> NumpressEncodeReport<Vec<u8>> {
    NumpressEncodeReport {
        output: Vec::new(),
        status,
        fixed_point: None,
        used_maximal_fixed_point_fallback: false,
    }
}
fn rejected(
    mut report: NumpressEncodeReport<Vec<u8>>,
    error: Error,
) -> NumpressEncodeReport<Vec<u8>> {
    report.status = NumpressEncodeStatus::Rejected(NumpressRejection::Codec(error.to_string()));
    report
}
pub(crate) fn encode_raw(
    input: &[f64],
    config: &NumpressConfig,
    work: &mut Work,
) -> Result<NumpressEncodeReport<Vec<u8>>> {
    if input.is_empty() {
        return Ok(report(NumpressEncodeStatus::EmptyInput));
    }
    if config.compression == NumpressCompression::None {
        return Ok(report(NumpressEncodeStatus::Disabled));
    }
    work.value_count(input.len())?;
    let mut report = report(NumpressEncodeStatus::Encoded);
    let mode = config.compression;
    let mut factor = config.fixed_point;
    if mode != NumpressCompression::Pic {
        if config.estimate_fixed_point {
            let limits = work.raw_call(input.len(), 16, 0)?;
            let estimate = match mode {
                NumpressCompression::Linear if config.linear_fp_mass_acc > 0.0 => {
                    raw::optimal_linear_fixed_point_mass_with_limits(
                        input,
                        config.linear_fp_mass_acc,
                        &limits,
                    )
                }
                NumpressCompression::Linear => {
                    raw::optimal_linear_fixed_point_with_limits(input, &limits)
                }
                NumpressCompression::Slof => {
                    raw::optimal_slof_fixed_point_with_limits(input, &limits)
                }
                _ => unreachable!(),
            };
            factor = match estimate {
                Ok(value) => value,
                Err(error) => return Ok(rejected(report, error)),
            };
            if mode == NumpressCompression::Linear
                && config.linear_fp_mass_acc > 0.0
                && factor < 0.0
            {
                report.used_maximal_fixed_point_fallback = true;
                let limits = work.raw_call(input.len(), 16, 0)?;
                factor = match raw::optimal_linear_fixed_point_with_limits(input, &limits) {
                    Ok(value) => value,
                    Err(error) => return Ok(rejected(report, error)),
                };
            }
        }
        report.fixed_point = Some(factor);
    }
    let capacity = match mode {
        NumpressCompression::Linear => add(mul(input.len(), 5)?, 8)?,
        NumpressCompression::Pic => mul(input.len(), 5)?,
        NumpressCompression::Slof => add(mul(input.len(), 2)?, 8)?,
        _ => unreachable!(),
    };
    work.binary(capacity)?;
    let limits = work.raw_call(input.len(), 64, capacity)?;
    let encoded = match mode {
        NumpressCompression::Linear => raw::encode_linear_with_limits(input, factor, &limits),
        NumpressCompression::Pic => raw::encode_pic_with_limits(input, &limits),
        NumpressCompression::Slof => raw::encode_slof_with_limits(input, factor, &limits),
        _ => unreachable!(),
    };
    let encoded = match encoded {
        Ok(bytes) => bytes,
        Err(error) => return Ok(rejected(report, error)),
    };
    if config.error_tolerance > 0.0 {
        let mut limits = work.raw_call(encoded.len(), 64, mul(input.len(), 8)?)?;
        // Valid encoding produces exactly this many values. Checking it bounds
        // verification allocation without charging arbitrary output heuristics.
        limits.max_values = input.len();
        let decoded = raw_decode(&encoded, mode, &limits);
        let decoded = match decoded {
            Ok(values) => values,
            Err(error) => return Ok(rejected(report, error)),
        };
        if decoded.len() != input.len() {
            return Ok(rejected(report, invalid("verification count mismatch")));
        }
        work.spend(mul(input.len(), 12)?)?;
        if let Some(failure) = accuracy_failure(input, &decoded, config) {
            report.status = NumpressEncodeStatus::Rejected(NumpressRejection::Accuracy(failure));
            return Ok(report);
        }
    }
    report.output = encoded;
    Ok(report)
}
fn accuracy_failure(
    input: &[f64],
    decoded: &[f64],
    config: &NumpressConfig,
) -> Option<NumpressAccuracyFailure> {
    for index in (0..input.len()).rev() {
        let original = input[index];
        let actual = decoded[index];
        let failed = if config.compression == NumpressCompression::Pic {
            !actual.is_finite() || (original - actual).abs() >= 1.0
        } else if !original.is_finite() || !actual.is_finite() {
            true
        } else if original == 0.0 {
            actual.abs() > config.error_tolerance
        } else if actual == 0.0 {
            original.abs() > config.error_tolerance
        } else {
            (1.0 - original / actual).abs() > config.error_tolerance
        };
        if failed {
            return Some(NumpressAccuracyFailure {
                index,
                original,
                decoded: actual,
            });
        }
    }
    None
}
pub(crate) fn encode_text(
    input: &[f64],
    zlib: bool,
    config: &NumpressConfig,
    work: &mut Work,
) -> Result<NumpressEncodeReport<String>> {
    let raw = encode_raw(input, config, work)?;
    let output = if raw.is_encoded() {
        encode_binary(raw.output, zlib, work)?
    } else {
        String::new()
    };
    Ok(NumpressEncodeReport {
        output,
        status: raw.status,
        fixed_point: raw.fixed_point,
        used_maximal_fixed_point_fallback: raw.used_maximal_fixed_point_fallback,
    })
}
pub(crate) fn encode_binary(input: Vec<u8>, zlib: bool, work: &mut Work) -> Result<String> {
    work.binary(input.len())?;
    let binary = if zlib {
        zlib_encode(&input, work)?
    } else {
        input
    };
    let length = mul(binary.len().div_ceil(3), 4)?;
    work.text(length)?;
    work.spend(add(binary.len(), length)?)?;
    let mut output = work.vector(length)?;
    output.resize(length, 0);
    let written = STANDARD
        .encode_slice(&binary, &mut output)
        .map_err(|_| invalid("base64 output bound"))?;
    output.truncate(written);
    String::from_utf8(output).map_err(|_| invalid("base64 emitted invalid UTF-8"))
}

fn raw_decode(data: &[u8], mode: NumpressCompression, limits: &NumpressLimits) -> Result<Vec<f64>> {
    match mode {
        NumpressCompression::None => Ok(Vec::new()),
        NumpressCompression::Linear => raw::decode_linear_with_limits(data, limits),
        NumpressCompression::Pic => raw::decode_pic_with_limits(data, limits),
        NumpressCompression::Slof => raw::decode_slof_with_limits(data, limits),
    }
}
pub(crate) fn decode_raw(
    input: &[u8],
    mode: NumpressCompression,
    work: &mut Work,
) -> Result<Vec<f64>> {
    if input.is_empty() || mode == NumpressCompression::None {
        return Ok(Vec::new());
    }
    work.binary(input.len())?;
    let upper = match mode {
        NumpressCompression::Linear | NumpressCompression::Pic => mul(input.len(), 2)?,
        NumpressCompression::Slof => input.len() / 2,
        _ => 0,
    }
    .min(work.limits.raw.max_values);
    let limits = work.raw_call(input.len(), 64, mul(upper, 8)?)?;
    raw_decode(input, mode, &limits)
}

const CHUNK: usize = 16 * 1024;
const ZLIB_STATE_BYTES: usize = 1024 * 1024;
pub(crate) fn zlib_encode(input: &[u8], work: &mut Work) -> Result<Vec<u8>> {
    work.allocate(ZLIB_STATE_BYTES)?;
    let mut encoder = Compress::new(Compression::default(), true);
    let mut result = Vec::new();
    let mut offset = 0usize;
    let mut chunk = [0u8; CHUNK];
    loop {
        let end = offset.saturating_add(CHUNK).min(input.len());
        work.spend(add(CHUNK, end - offset)?)?;
        let before_in = encoder.total_in();
        let before_out = encoder.total_out();
        let flush = if end == input.len() {
            FlushCompress::Finish
        } else {
            FlushCompress::None
        };
        let status = encoder
            .compress(&input[offset..end], &mut chunk, flush)
            .map_err(|_| invalid("zlib compression failed"))?;
        let consumed = (encoder.total_in() - before_in) as usize;
        let produced = (encoder.total_out() - before_out) as usize;
        offset += consumed;
        work.append(&mut result, &chunk[..produced])?;
        if status == Status::StreamEnd {
            if offset != input.len() {
                return Err(invalid("unfinished zlib input"));
            }
            return Ok(result);
        }
        if consumed == 0 && produced == 0 {
            return Err(invalid("stalled zlib compressor"));
        }
    }
}
pub(crate) fn zlib_decode(input: &[u8], work: &mut Work) -> Result<Vec<u8>> {
    work.allocate(ZLIB_STATE_BYTES)?;
    let mut decoder = Decompress::new(true);
    let mut result = Vec::new();
    let mut offset = 0usize;
    let mut chunk = [0u8; CHUNK];
    loop {
        let end = offset.saturating_add(CHUNK).min(input.len());
        work.spend(add(CHUNK, end - offset)?)?;
        let before_in = decoder.total_in();
        let before_out = decoder.total_out();
        let status = decoder
            .decompress(&input[offset..end], &mut chunk, FlushDecompress::None)
            .map_err(|_| invalid("invalid zlib stream"))?;
        let consumed = (decoder.total_in() - before_in) as usize;
        let produced = (decoder.total_out() - before_out) as usize;
        offset += consumed;
        work.append(&mut result, &chunk[..produced])?;
        if status == Status::StreamEnd {
            if offset != input.len() {
                return Err(invalid("trailing or concatenated zlib data"));
            }
            return Ok(result);
        }
        if consumed == 0 && produced == 0 {
            return Err(invalid("truncated or stalled zlib stream"));
        }
    }
}

pub(crate) struct Work {
    pub(crate) limits: NumpressCoderLimits,
    work: usize,
    bytes: usize,
}
impl Work {
    #[cfg(feature = "sqmass")]
    pub(crate) fn remaining_bytes(&self) -> usize {
        self.bytes
    }

    #[cfg(feature = "sqmass")]
    pub(crate) fn remaining_work(&self) -> usize {
        self.work
    }

    pub(crate) fn new(limits: NumpressCoderLimits) -> Self {
        Self {
            work: limits.raw.max_work,
            bytes: limits.max_total_bytes,
            limits,
        }
    }
    pub(crate) fn spend(&mut self, n: usize) -> Result<()> {
        self.work = self
            .work
            .checked_sub(n)
            .ok_or_else(|| invalid("cumulative work limit exceeded"))?;
        Ok(())
    }
    pub(crate) fn allocate(&mut self, n: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(n)
            .ok_or_else(|| invalid("cumulative allocation limit exceeded"))?;
        Ok(())
    }
    pub(crate) fn value_count(&self, n: usize) -> Result<()> {
        if n > self.limits.raw.max_values {
            return Err(invalid("value count limit exceeded"));
        }
        Ok(())
    }
    pub(crate) fn binary(&self, n: usize) -> Result<()> {
        if n > self.limits.raw.max_encoded_bytes {
            return Err(invalid("binary byte limit exceeded"));
        }
        Ok(())
    }
    pub(crate) fn text(&self, n: usize) -> Result<()> {
        if n > self.limits.max_text_bytes {
            return Err(invalid("base64 text limit exceeded"));
        }
        Ok(())
    }
    fn raw_call(&mut self, n: usize, factor: usize, bytes: usize) -> Result<NumpressLimits> {
        let cost = mul(n, factor)?;
        self.spend(cost)?;
        self.allocate(bytes)?;
        Ok(NumpressLimits {
            max_work: cost,
            ..self.limits.raw
        })
    }
    pub(crate) fn vector<T>(&mut self, n: usize) -> Result<Vec<T>> {
        self.allocate(mul(n, std::mem::size_of::<T>())?)?;
        let mut value = Vec::new();
        value
            .try_reserve_exact(n)
            .map_err(|_| invalid("allocation failed"))?;
        Ok(value)
    }
    fn append(&mut self, output: &mut Vec<u8>, bytes: &[u8]) -> Result<()> {
        let required = add(output.len(), bytes.len())?;
        self.binary(required)?;
        if required > output.capacity() {
            let capacity = required.max(
                output
                    .capacity()
                    .saturating_mul(2)
                    .min(self.limits.raw.max_encoded_bytes),
            );
            self.allocate(capacity)?;
            self.spend(output.len())?;
            output
                .try_reserve_exact(capacity - output.len())
                .map_err(|_| invalid("allocation failed"))?;
        }
        output.extend_from_slice(bytes);
        Ok(())
    }
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(|| invalid("size overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(|| invalid("size overflow"))
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(format!("MSNumpressCoder: {message}"))
}
