// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Binary transport only; ordinary mzML options and path APIs are unchanged.
use super::*;
use coder::{NumpressCoderLimits, NumpressEncodeStatus};

/// Three source Numpress configurations, separate from ordinary writer options.
/// Limits cover all binary arrays, codec attempts, verification, fallback and
/// retained encoded text across the operation, rather than resetting per array.
#[derive(Clone, Debug, Default)]
pub struct NumpressWriteOptions {
    pub binary: WriteOptions,
    pub mass_time: NumpressConfig,
    pub intensity: NumpressConfig,
    pub float_data_array: NumpressConfig,
    pub limits: NumpressCoderLimits,
}
/// Array counts describe the final encoding, including checked ordinary fallback.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NumpressWriteReport {
    pub encoded_arrays: usize,
    /// Requested Numpress but emitted ordinary data: empty input, rejected codec
    /// or accuracy check, or a canonical role that cannot declare float64.
    pub fallback_arrays: usize,
    /// Numpress was not requested (including all integer/string arrays).
    pub ordinary_arrays: usize,
}

/// Preflight and encode every array before emitting XML. Codec/accuracy rejection
/// falls back to ordinary binary as in the source; transport/resource errors
/// return an error before output. I/O errors may leave partial output, as usual.
pub fn write_with_numpress(
    writer: impl Write,
    experiment: &MSExperiment,
    options: &NumpressWriteOptions,
) -> Result<NumpressWriteReport> {
    experiment_header_guard(experiment)?;
    let mut work = coder::Work::new(options.limits);
    work.spend(add(
        experiment.spectra.len(),
        experiment.chromatograms.len(),
    )?)?;
    let mut count = 0usize;
    for (floats, integers, strings) in experiment
        .spectra
        .iter()
        .map(|s| {
            (
                &s.float_data_arrays,
                &s.integer_data_arrays,
                &s.string_data_arrays,
            )
        })
        .chain(experiment.chromatograms.iter().map(|c| {
            (
                &c.float_data_arrays,
                &c.integer_data_arrays,
                &c.string_data_arrays,
            )
        }))
    {
        count = add(
            count,
            add(2, add(floats.len(), add(integers.len(), strings.len())?)?)?,
        )?;
    }
    work.spend(count)?;
    let arrays = work.vector(count)?;
    // Bound lengths and charge validation before the general validator visits
    // any peak or auxiliary scalar. Per-array limits cannot be deferred until
    // encoding: a late invalid value must not force an over-limit scan first.
    for spectrum in &experiment.spectra {
        preflight_values(
            &mut work,
            spectrum.len(),
            &spectrum.float_data_arrays,
            &spectrum.integer_data_arrays,
            &spectrum.string_data_arrays,
        )?;
    }
    for chromatogram in &experiment.chromatograms {
        preflight_values(
            &mut work,
            chromatogram.len(),
            &chromatogram.float_data_arrays,
            &chromatogram.integer_data_arrays,
            &chromatogram.string_data_arrays,
        )?;
    }
    // The existing whole-document preflight checks unsupported metadata and all
    // native values before any output, including description metadata on arrays.
    validate_write(experiment)?;
    let mut prepared = Preparation {
        work,
        options,
        arrays,
        report: NumpressWriteReport::default(),
    };
    for spectrum in &experiment.spectra {
        prepared.floats(
            spectrum.peaks.iter().map(|p| p.mz),
            Encoding::Float64,
            &options.mass_time,
            None,
        )?;
        prepared.floats(
            spectrum.peaks.iter().map(|p| f64::from(p.intensity)),
            Encoding::Float32,
            &options.intensity,
            None,
        )?;
        prepared.auxiliary(
            &spectrum.float_data_arrays,
            &spectrum.integer_data_arrays,
            &spectrum.string_data_arrays,
        )?;
    }
    for chromatogram in &experiment.chromatograms {
        prepared.floats(
            chromatogram.peaks.iter().map(|p| p.rt),
            Encoding::Float64,
            &options.mass_time,
            None,
        )?;
        prepared.floats(
            chromatogram.peaks.iter().map(|p| f64::from(p.intensity)),
            Encoding::Float32,
            &options.intensity,
            None,
        )?;
        prepared.auxiliary(
            &chromatogram.float_data_arrays,
            &chromatogram.integer_data_arrays,
            &chromatogram.string_data_arrays,
        )?;
    }
    write_impl(
        writer,
        experiment,
        &options.binary,
        &mut Some(prepared.arrays.into_iter()),
    )?;
    Ok(prepared.report)
}

fn preflight_values(
    work: &mut coder::Work,
    count: usize,
    floats: &[DataArray<f32>],
    integers: &[DataArray<i32>],
    strings: &[DataArray<String>],
) -> Result<()> {
    work.value_count(count)?;
    work.spend(mul(count, 2)?)?;
    let arrays = add(floats.len(), add(integers.len(), strings.len())?)?;
    work.spend(mul(arrays, 8)?)?;
    for array in floats {
        work.value_count(array.data.len())?;
        work.spend(array.data.len())?;
    }
    for array in integers {
        work.value_count(array.data.len())?;
        work.spend(array.data.len())?;
    }
    for array in strings {
        work.value_count(array.data.len())?;
        work.spend(array.data.len())?;
        for value in &array.data {
            work.spend(mul(value.len(), 2)?)?;
        }
    }
    // Array descriptions are unsupported. Reject with an O(1) check before the
    // general validator can traverse their metadata/processing payload.
    if floats.iter().any(DataArray::has_description_metadata)
        || integers.iter().any(DataArray::has_description_metadata)
        || strings.iter().any(DataArray::has_description_metadata)
    {
        return Err(Error::Unsupported(
            "mzML array description metadata or processing is not represented".into(),
        ));
    }
    // Name validation and ordered duplicate-name comparisons are binary-array
    // work too. The logarithmic factor bounds the existing BTreeSet traversal.
    let factor = mul(8, add(arrays.max(1).ilog2() as usize, 2)?)?;
    for name in floats
        .iter()
        .map(|a| &a.name)
        .chain(integers.iter().map(|a| &a.name))
        .chain(strings.iter().map(|a| &a.name))
    {
        work.spend(mul(name.len(), factor)?)?;
    }
    Ok(())
}

pub(super) struct PreparedArray {
    pub(super) encoded: String,
    pub(super) encoding: Encoding,
    pub(super) mode: NumpressCompression,
}
struct Preparation<'a> {
    work: coder::Work,
    options: &'a NumpressWriteOptions,
    arrays: Vec<PreparedArray>,
    report: NumpressWriteReport,
}
impl Preparation<'_> {
    fn floats(
        &mut self,
        input: impl ExactSizeIterator<Item = f64>,
        encoding: Encoding,
        config: &NumpressConfig,
        name: Option<&str>,
    ) -> Result<()> {
        let count = input.len();
        self.work.value_count(count)?;
        self.work.spend(count)?;
        let mut values = self.work.vector(count)?;
        values.extend(input);
        let requested = config.compression != NumpressCompression::None;
        // A canonical float32-only role cannot be serialized with the source's
        // unconditional float64 Numpress CV. Preserve its identity via fallback.
        let compatible =
            name.is_none_or(|n| check_canonical_encoding(n, Encoding::Float64).is_ok());
        if requested && compatible {
            let encoded = coder::encode_text(
                &values,
                self.options.binary.zlib_compression,
                config,
                &mut self.work,
            )?;
            if encoded.status == NumpressEncodeStatus::Encoded {
                self.arrays.push(PreparedArray {
                    encoded: encoded.output,
                    encoding: Encoding::Float64,
                    mode: config.compression,
                });
                self.report.encoded_arrays += 1;
                return Ok(());
            }
        }
        if requested {
            self.report.fallback_arrays += 1;
        } else {
            self.report.ordinary_arrays += 1;
        }
        self.work.spend(count)?;
        let length = mul(count, encoding.width().expect("floating width"))?;
        self.work.binary(length)?;
        let mut bytes = self.work.vector(length)?;
        for value in values {
            if encoding == Encoding::Float32 {
                bytes.extend_from_slice(&(value as f32).to_le_bytes());
            } else {
                bytes.extend_from_slice(&value.to_le_bytes());
            }
        }
        self.ordinary(bytes, encoding)
    }
    fn ordinary(&mut self, bytes: Vec<u8>, encoding: Encoding) -> Result<()> {
        let encoded =
            coder::encode_binary(bytes, self.options.binary.zlib_compression, &mut self.work)?;
        self.arrays.push(PreparedArray {
            encoded,
            encoding,
            mode: NumpressCompression::None,
        });
        Ok(())
    }
    fn auxiliary(
        &mut self,
        floats: &[DataArray<f32>],
        integers: &[DataArray<i32>],
        strings: &[DataArray<String>],
    ) -> Result<()> {
        for array in floats {
            self.floats(
                array.data.iter().map(|&v| f64::from(v)),
                Encoding::Float32,
                &self.options.float_data_array,
                Some(&array.name),
            )?;
        }
        for array in integers {
            self.work.value_count(array.data.len())?;
            self.work.spend(array.data.len())?;
            let encoding = integer_array_encoding(&array.name);
            let length = mul(array.data.len(), encoding.width().unwrap())?;
            self.work.binary(length)?;
            let mut bytes = self.work.vector(length)?;
            for value in &array.data {
                if encoding == Encoding::Int32 {
                    bytes.extend_from_slice(&value.to_le_bytes());
                } else {
                    bytes.extend_from_slice(&i64::from(*value).to_le_bytes());
                }
            }
            self.report.ordinary_arrays += 1;
            self.ordinary(bytes, encoding)?;
        }
        for array in strings {
            self.work.value_count(array.data.len())?;
            self.work.spend(array.data.len())?;
            let mut length = 0usize;
            for value in &array.data {
                length = add(length, add(value.len(), 1)?)?;
            }
            self.work.spend(length)?;
            self.work.binary(length)?;
            let mut bytes = self.work.vector(length)?;
            for value in &array.data {
                bytes.extend_from_slice(value.as_bytes());
                bytes.push(0);
            }
            self.report.ordinary_arrays += 1;
            self.ordinary(bytes, Encoding::Ascii)?;
        }
        Ok(())
    }
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| invalid("Numpress size overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b)
        .ok_or_else(|| invalid("Numpress size overflow"))
}

pub(super) fn compression(accession: &str) -> Option<(NumpressCompression, bool)> {
    use NumpressCompression::*;
    Some(match accession {
        "MS:1002312" => (Linear, false),
        "MS:1002313" => (Pic, false),
        "MS:1002314" => (Slof, false),
        "MS:1002746" => (Linear, true),
        "MS:1002747" => (Pic, true),
        "MS:1002748" => (Slof, true),
        _ => return Option::None,
    })
}
pub(super) fn compression_term(
    mode: NumpressCompression,
    zlib: bool,
) -> Option<(&'static str, &'static str)> {
    use NumpressCompression::*;
    Some(match (mode, zlib) {
        (Linear, false) => ("MS:1002312", "MS-Numpress linear prediction compression"),
        (Pic, false) => ("MS:1002313", "MS-Numpress positive integer compression"),
        (Slof, false) => ("MS:1002314", "MS-Numpress short logged float compression"),
        (Linear, true) => (
            "MS:1002746",
            "MS-Numpress linear prediction compression followed by zlib compression",
        ),
        (Pic, true) => (
            "MS:1002747",
            "MS-Numpress positive integer compression followed by zlib compression",
        ),
        (Slof, true) => (
            "MS:1002748",
            "MS-Numpress short logged float compression followed by zlib compression",
        ),
        (None, _) => return Option::None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn decode(work: &mut coder::Work) -> Result<(Kind, Values)> {
        let text = "QWR64UAAAADo//8/0P//f1kSgA==";
        Binary {
            kind: Some(Kind::Mz),
            encoding: Some(Encoding::Float64),
            numpress: Some(NumpressCompression::Linear),
            encoded: text.into(),
            encoded_length: text.len(),
            has_binary: true,
            ..Default::default()
        }
        .decode(4, &ReadOptions::default(), &mut 1024, &mut 100, work)
    }
    #[test]
    fn decoder_counters_span_arrays_and_fail_before_resetting() {
        let mut limits = NumpressCoderLimits::default();
        limits.raw.max_work = 2000;
        let mut work = coder::Work::new(limits);
        assert!(decode(&mut work).is_ok());
        assert!(decode(&mut work).is_err());
        let limits = NumpressCoderLimits {
            max_total_bytes: 80,
            ..Default::default()
        };
        let mut work = coder::Work::new(limits);
        assert!(decode(&mut work).is_ok());
        assert!(decode(&mut work).is_err());
    }
}
