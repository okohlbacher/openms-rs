/*
        Native Rust MSNumpress codecs
        johan.teleman@immun.lth.se

        This distribution goes under the BSD 3-clause license. If you prefer to use Apache
        version 2.0, that is also available at https://github.com/fickludd/ms-numpress
        Copyright (c) 2013, Johan Teleman
        All rights reserved.

        Redistribution and use in source and binary forms, with or without modification,
        are permitted provided that the following conditions are met:

*         Redistributions of source code must retain the above copyright notice, this list
        of conditions and the following disclaimer.
*        Redistributions in binary form must reproduce the above copyright notice, this
        list of conditions and the following disclaimer in the documentation and/or other
        materials provided with the distribution.
*        Neither the name of the Lund University nor the names of its contributors may be
        used to endorse or promote products derived from this software without specific
        prior written permission.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
        EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
        OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
        SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
        OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
        HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
        OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// $Maintainer: OpenMS Rust contributors $

//! Raw MSNumpress linear, positive-integer, short-logged-float and Safe codecs.
//! This module does not apply base64, zlib or mzML CV/transport conventions.

use crate::{Error, Result};

/// Independent input/output and work bounds for one raw operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NumpressLimits {
    pub max_values: usize,
    pub max_encoded_bytes: usize,
    pub max_work: usize,
}
impl Default for NumpressLimits {
    fn default() -> Self {
        Self {
            max_values: 10_000_000,
            max_encoded_bytes: 128 * 1024 * 1024,
            max_work: 500_000_000,
        }
    }
}
impl NumpressLimits {
    fn values(self, count: usize) -> Result<()> {
        if count > self.max_values {
            return Err(invalid("value count limit exceeded"));
        }
        Ok(())
    }
    fn bytes(self, count: usize) -> Result<()> {
        if count > self.max_encoded_bytes {
            return Err(invalid("encoded byte limit exceeded"));
        }
        Ok(())
    }
    fn work(self, count: usize, factor: usize) -> Result<()> {
        if multiply(count, factor)? > self.max_work {
            return Err(invalid("work limit exceeded"));
        }
        Ok(())
    }
    fn encoding(self, count: usize, bytes: usize) -> Result<Vec<u8>> {
        self.values(count)?;
        self.bytes(bytes)?;
        self.work(count, 64)?;
        vector(bytes)
    }
    fn decoding(self, data: &[u8]) -> Result<()> {
        self.bytes(data.len())?;
        // Includes both the validation/counting pass and the materialization pass.
        self.work(data.len(), 64)
    }
}

/// Source maximal linear fixed point; empty input returns zero.
pub fn optimal_linear_fixed_point(data: &[f64]) -> Result<f64> {
    optimal_linear_fixed_point_with_limits(data, &NumpressLimits::default())
}
pub fn optimal_linear_fixed_point_with_limits(
    data: &[f64],
    limits: &NumpressLimits,
) -> Result<f64> {
    limits.values(data.len())?;
    limits.work(data.len(), 16)?;
    for &value in data {
        finite(value)?;
    }
    let Some(&first) = data.first() else {
        return Ok(0.0);
    };
    let mut maximum = first;
    if data.len() > 1 {
        maximum = maximum.max(data[1]);
    }
    for window in data.windows(3) {
        let extrapolated = window[1] + (window[1] - window[0]);
        maximum = maximum.max(((window[2] - extrapolated).abs() + 1.0).ceil());
    }
    finite((f64::from(i32::MAX) / maximum).floor())
}

/// Source accuracy helper: fewer than three values returns zero without
/// consuming their values or accuracy; impossible requested accuracy returns -1.
pub fn optimal_linear_fixed_point_mass(data: &[f64], accuracy: f64) -> Result<f64> {
    optimal_linear_fixed_point_mass_with_limits(data, accuracy, &NumpressLimits::default())
}
pub fn optimal_linear_fixed_point_mass_with_limits(
    data: &[f64],
    accuracy: f64,
    limits: &NumpressLimits,
) -> Result<f64> {
    limits.values(data.len())?;
    if data.len() < 3 {
        return Ok(0.0);
    }
    finite(accuracy)?;
    let requested = 0.5 / accuracy;
    let maximum = optimal_linear_fixed_point_with_limits(data, limits)?;
    if requested > maximum {
        Ok(-1.0)
    } else {
        finite(requested)
    }
}

/// Source SLOF helper uses 65534 and a minimum maximum-log value of one.
/// As in the source, nonpositive log-domain inputs do not raise that maximum;
/// this helper does not guarantee those values can subsequently be encoded.
pub fn optimal_slof_fixed_point(data: &[f64]) -> Result<f64> {
    optimal_slof_fixed_point_with_limits(data, &NumpressLimits::default())
}
pub fn optimal_slof_fixed_point_with_limits(data: &[f64], limits: &NumpressLimits) -> Result<f64> {
    limits.values(data.len())?;
    limits.work(data.len(), 16)?;
    if data.is_empty() {
        return Ok(0.0);
    }
    let mut maximum = 1.0_f64;
    for &value in data {
        finite(value)?;
        let logged = (value + 1.0).ln();
        if maximum < logged {
            maximum = logged;
        }
    }
    finite((65534.0 / maximum).floor())
}

pub fn encode_linear(data: &[f64], fixed_point: f64) -> Result<Vec<u8>> {
    encode_linear_with_limits(data, fixed_point, &NumpressLimits::default())
}
/// The first two quantized i64 values serialize their low 32 bits, as in source.
/// A caller-selected factor can therefore lose high bits; use the optimal helper
/// for ordinary positive m/z data. Later signed i32 residuals are checked.
pub fn encode_linear_with_limits(
    data: &[f64],
    fixed_point: f64,
    limits: &NumpressLimits,
) -> Result<Vec<u8>> {
    let capacity = add(multiply(data.len(), 5)?, 8)?;
    let mut output = limits.encoding(data.len(), capacity)?;
    output.extend_from_slice(&fixed_point.to_be_bytes());
    if data.is_empty() {
        return Ok(output);
    }
    finite(fixed_point)?;
    let mut previous = quantize(data[0], fixed_point)?;
    output.extend_from_slice(&(previous as u32).to_le_bytes());
    if data.len() == 1 {
        return Ok(output);
    }
    let mut current = quantize(data[1], fixed_point)?;
    output.extend_from_slice(&(current as u32).to_le_bytes());
    let mut nibbles = NibbleWriter {
        output: &mut output,
        pending: None,
    };
    for &value in &data[2..] {
        let next = quantize(value, fixed_point)?;
        let predicted = predict(previous, current)?;
        let delta = next.checked_sub(predicted).ok_or_else(overflow)?;
        let delta = i32::try_from(delta).map_err(|_| invalid("linear residual exceeds i32"))?;
        nibbles.integer(delta as u32);
        previous = current;
        current = next;
    }
    nibbles.finish();
    Ok(output)
}

pub fn decode_linear(data: &[u8]) -> Result<Vec<f64>> {
    decode_linear_with_limits(data, &NumpressLimits::default())
}
pub fn decode_linear_with_limits(data: &[u8], limits: &NumpressLimits) -> Result<Vec<f64>> {
    decode(data, limits, walk_linear)
}

pub fn encode_pic(data: &[f64]) -> Result<Vec<u8>> {
    encode_pic_with_limits(data, &NumpressLimits::default())
}
/// Preserves the compiled source guard: -0.5 <= input and input+0.5 <= INT_MAX.
/// The source header's advertised unsigned maximum is wider than that guard.
pub fn encode_pic_with_limits(data: &[f64], limits: &NumpressLimits) -> Result<Vec<u8>> {
    let mut output = limits.encoding(data.len(), multiply(data.len(), 5)?)?;
    let mut nibbles = NibbleWriter {
        output: &mut output,
        pending: None,
    };
    for &value in data {
        finite(value)?;
        if value < -0.5 || value + 0.5 > f64::from(i32::MAX) {
            return Err(invalid("PIC value exceeds source INT_MAX rounding bounds"));
        }
        nibbles.integer((value + 0.5) as u32);
    }
    nibbles.finish();
    Ok(output)
}
pub fn decode_pic(data: &[u8]) -> Result<Vec<f64>> {
    decode_pic_with_limits(data, &NumpressLimits::default())
}
/// Decoding retains all u32 patterns, including values above the encoder's cap.
pub fn decode_pic_with_limits(data: &[u8], limits: &NumpressLimits) -> Result<Vec<f64>> {
    decode(data, limits, walk_pic)
}

pub fn encode_slof(data: &[f64], fixed_point: f64) -> Result<Vec<u8>> {
    encode_slof_with_limits(data, fixed_point, &NumpressLimits::default())
}
pub fn encode_slof_with_limits(
    data: &[f64],
    fixed_point: f64,
    limits: &NumpressLimits,
) -> Result<Vec<u8>> {
    let mut output = limits.encoding(data.len(), add(multiply(data.len(), 2)?, 8)?)?;
    output.extend_from_slice(&fixed_point.to_be_bytes());
    if !data.is_empty() {
        finite(fixed_point)?;
    }
    for &value in data {
        finite(value)?;
        // Literal source log(1+x), not ln_1p; changes near zero affect bytes.
        let scaled = finite((value + 1.0).ln() * fixed_point)?;
        if scaled > 65535.0 {
            return Err(invalid("SLOF value exceeds USHRT_MAX"));
        }
        let rounded = scaled + 0.5;
        // C++ floating-to-unsigned conversion is defined when truncation fits.
        if !(rounded > -1.0 && rounded < 65536.0) {
            return Err(invalid("SLOF rounded value is outside u16"));
        }
        output.extend_from_slice(&(rounded as u16).to_le_bytes());
    }
    Ok(output)
}
pub fn decode_slof(data: &[u8]) -> Result<Vec<f64>> {
    decode_slof_with_limits(data, &NumpressLimits::default())
}
pub fn decode_slof_with_limits(data: &[u8], limits: &NumpressLimits) -> Result<Vec<f64>> {
    decode(data, limits, walk_slof)
}

pub fn encode_safe(data: &[f64]) -> Result<Vec<u8>> {
    encode_safe_with_limits(data, &NumpressLimits::default())
}
/// Source double-residual codec. "Safe" is the source name; arbitrary f64
/// sequences need not roundtrip bit-exactly after prediction and subtraction.
pub fn encode_safe_with_limits(data: &[f64], limits: &NumpressLimits) -> Result<Vec<u8>> {
    let mut output = limits.encoding(data.len(), multiply(data.len(), 8)?)?;
    for (index, &value) in data.iter().enumerate() {
        finite(value)?;
        let encoded = if index < 2 {
            value
        } else {
            finite(value - (data[index - 1] + (data[index - 1] - data[index - 2])))?
        };
        output.extend_from_slice(&encoded.to_be_bytes());
    }
    Ok(output)
}
pub fn decode_safe(data: &[u8]) -> Result<Vec<f64>> {
    decode_safe_with_limits(data, &NumpressLimits::default())
}
pub fn decode_safe_with_limits(data: &[u8], limits: &NumpressLimits) -> Result<Vec<f64>> {
    decode(data, limits, walk_safe)
}

// A validation/counting pass bounds the exact decoded output before allocation.
// The second pass never invokes a user callback, and errors expose no partial Vec.
type Walk = fn(&[u8], &mut dyn FnMut(f64) -> Result<()>) -> Result<()>;
fn decode(data: &[u8], limits: &NumpressLimits, walk: Walk) -> Result<Vec<f64>> {
    limits.decoding(data)?;
    let mut count = 0;
    walk(data, &mut |_| {
        count = add(count, 1)?;
        limits.values(count)
    })?;
    let mut output = vector(count)?;
    walk(data, &mut |value| {
        output.push(value);
        Ok(())
    })?;
    Ok(output)
}
fn walk_linear(data: &[u8], emit: &mut dyn FnMut(f64) -> Result<()>) -> Result<()> {
    if data.len() < 8 {
        return Err(truncated());
    }
    if data.len() == 8 {
        return Ok(());
    } // unused fixed-point bits
    let fixed = fixed_point(data)?;
    if data.len() < 12 {
        return Err(truncated());
    }
    let mut previous = i64::from(u32::from_le_bytes(data[8..12].try_into().unwrap()));
    emit(finite(previous as f64 / fixed)?)?;
    if data.len() == 12 {
        return Ok(());
    }
    if data.len() < 16 {
        return Err(truncated());
    }
    let mut current = i64::from(u32::from_le_bytes(data[12..16].try_into().unwrap()));
    emit(finite(current as f64 / fixed)?)?;
    let mut nibbles = NibbleReader::new(&data[16..]);
    while let Some(bits) = nibbles.integer()? {
        let next = predict(previous, current)?
            .checked_add(i64::from(bits as i32))
            .ok_or_else(overflow)?;
        emit(finite(next as f64 / fixed)?)?;
        previous = current;
        current = next;
    }
    Ok(())
}
fn walk_pic(data: &[u8], emit: &mut dyn FnMut(f64) -> Result<()>) -> Result<()> {
    let mut nibbles = NibbleReader::new(data);
    while let Some(bits) = nibbles.integer()? {
        emit(f64::from(bits))?;
    }
    Ok(())
}
fn walk_slof(data: &[u8], emit: &mut dyn FnMut(f64) -> Result<()>) -> Result<()> {
    if data.len() < 8 || data.len() % 2 != 0 {
        return Err(truncated());
    }
    if data.len() == 8 {
        return Ok(());
    }
    let fixed = fixed_point(data)?;
    for pair in data[8..].chunks_exact(2) {
        let value = u16::from_le_bytes(pair.try_into().unwrap());
        emit(finite((f64::from(value) / fixed).exp() - 1.0)?)?;
    }
    Ok(())
}
fn walk_safe(data: &[u8], emit: &mut dyn FnMut(f64) -> Result<()>) -> Result<()> {
    if data.len() % 8 != 0 {
        return Err(truncated());
    }
    let mut previous = 0.0;
    let mut current = 0.0;
    for (index, bytes) in data.chunks_exact(8).enumerate() {
        let encoded = finite(f64::from_be_bytes(bytes.try_into().unwrap()))?;
        let next = if index < 2 {
            encoded
        } else {
            finite((current + (current - previous)) + encoded)?
        };
        emit(next)?;
        previous = current;
        current = next;
    }
    Ok(())
}

fn predict(previous: i64, current: i64) -> Result<i64> {
    current
        .checked_sub(previous)
        .and_then(|delta| current.checked_add(delta))
        .ok_or_else(overflow)
}
fn quantize(value: f64, fixed: f64) -> Result<i64> {
    finite(value)?;
    let rounded = finite(value * fixed + 0.5)?;
    if !(-9223372036854775808.0..9223372036854775808.0).contains(&rounded) {
        return Err(overflow());
    }
    Ok(rounded as i64)
}
fn fixed_point(data: &[u8]) -> Result<f64> {
    // Callers check header length and selected results. SLOF with negative-zero
    // fixed point and a positive word has the finite source result -1.
    finite(f64::from_be_bytes(data[..8].try_into().unwrap()))
}
struct NibbleWriter<'a> {
    output: &'a mut Vec<u8>,
    pending: Option<u8>,
}
impl NibbleWriter<'_> {
    fn push(&mut self, nibble: u8) {
        if let Some(high) = self.pending.take() {
            self.output.push((high << 4) | (nibble & 15));
        } else {
            self.pending = Some(nibble & 15);
        }
    }
    fn integer(&mut self, value: u32) {
        let zeros = value.leading_zeros() / 4;
        let ones = (value.leading_ones() / 4).min(7);
        let leading = if zeros != 0 {
            self.push(zeros as u8);
            zeros
        } else if ones != 0 {
            self.push((ones + 8) as u8);
            ones
        } else {
            self.push(0);
            0
        };
        for index in 0..8 - leading {
            self.push((value >> (4 * index)) as u8);
        }
    }
    fn finish(self) {
        if let Some(high) = self.pending {
            self.output.push(high << 4);
        }
    }
}
struct NibbleReader<'a> {
    data: &'a [u8],
    index: usize,
}
impl<'a> NibbleReader<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, index: 0 }
    }
    fn read(&mut self) -> Result<u8> {
        let byte = *self.data.get(self.index / 2).ok_or_else(truncated)?;
        let value = if self.index % 2 == 0 {
            byte >> 4
        } else {
            byte & 15
        };
        self.index += 1;
        Ok(value)
    }
    fn integer(&mut self) -> Result<Option<u32>> {
        if self.index / 2 == self.data.len() {
            return Ok(None);
        }
        if self.index % 2 == 1
            && self.index / 2 + 1 == self.data.len()
            && self.data[self.index / 2] & 15 == 0
        {
            return Ok(None);
        }
        let head = self.read()?;
        let leading = if head <= 8 { head } else { head - 8 };
        let mut value = if head > 8 {
            u32::MAX << (4 * (8 - leading))
        } else {
            0
        };
        for index in 0..8 - leading {
            value |= u32::from(self.read()?) << (4 * index);
        }
        Ok(Some(value))
    }
}
fn vector<T>(capacity: usize) -> Result<Vec<T>> {
    multiply(capacity, std::mem::size_of::<T>())?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| invalid("allocation failed"))?;
    Ok(result)
}
fn multiply(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(|| invalid("size overflow"))
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(|| invalid("size overflow"))
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite numerical value"))
    }
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(format!("MSNumpress: {message}"))
}
fn overflow() -> Error {
    invalid("linear i64 conversion or prediction overflow")
}
fn truncated() -> Error {
    invalid("truncated or malformed byte stream")
}
