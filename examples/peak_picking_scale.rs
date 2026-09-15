// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Scale harness for [`PeakPickerHiRes`]: load a profile mzML run, centroid it,
//! and report the peak resident set size and wall time of every stage.
//!
//! The harness exists to produce the numbers behind the peak-picking scale work:
//! it is the only in-tree caller that drives `pick_experiment` over a
//! multi-gigabyte run, and it is deliberately free of the TOPP tool's parameter
//! handling so that the measurement covers loading, picking and writing only.
//!
//! ```text
//! peak_picking_scale <input.mzML> [--out <output.mzML>] [--in-place]
//!                                 [--ledger-probe] [--limit <n>] [--threads <n>]
//! ```
//!
//! * `--out` writes the centroided run; without it the harness only digests the
//!   picked peaks, so that the writer's own cost stays out of the measurement.
//! * `--in-place` drives [`PeakPickerHiRes::pick_experiment_in_place`] instead of
//!   [`PeakPickerHiRes::pick_experiment`].
//! * `--threads` picks on that many workers, as the TOPP `-threads` parameter
//!   does (`0` means every available core); the default is one. The digest it
//!   reports is the same at every count, which is the determinism contract of
//!   `openms::concept::parallel` measured rather than asserted.
//! * `--ledger-probe` reports the largest spectrum count the acquisition-copy
//!   ledger admits for this run's metadata, by replicating its metadata-only
//!   records.
//! * `--limit` keeps only the first `n` spectra of the run.
//!
//! Resident set sizes come from `/proc/self/status` and are reported only where
//! that file exists; wall times and digests are reported everywhere.
//!
//! The harness reads mzML, so its body is behind the `mzml` feature; without it
//! the example still builds and explains that, which keeps the crate's
//! `--no-default-features` builds — and `cargo test`, which builds every
//! example — compiling.

/// The harness proper, compiled only where the mzML reader exists.
#[cfg(feature = "mzml")]
mod harness {
    use openms::concept::parallel::Threads;
    use openms::format::mzml;
    use openms::kernel::{MSExperiment, MSSpectrum};
    use openms::processing::peak_picking::{PeakPickerHiRes, PickingCompatibility};
    use openms::{Error, Result};
    use std::io::Write;
    use std::time::Instant;

    /// Read options generous enough for a multi-gigabyte profile run. The picker
    /// lane does not own the reader's defaults; the harness states what it needs.
    fn read_options() -> mzml::ReadOptions {
        mzml::ReadOptions {
            max_xml_bytes: 64 * 1024 * 1024 * 1024,
            max_array_bytes: 4 * 1024 * 1024 * 1024,
            max_total_peaks: 20_000_000_000,
            max_records: 100_000_000,
            max_total_array_bytes: 64 * 1024 * 1024 * 1024,
            max_total_array_elements: 20_000_000_000,
            max_total_arrays: 100_000_000,
            max_param_groups: 10_000_000,
            max_total_params: 1_000_000_000,
            max_param_bytes: 8 * 1024 * 1024 * 1024,
            ..mzml::ReadOptions::default()
        }
    }

    /// The picker a TOPP tool reproducing the C++ output configures: the source
    /// behaviours, since vendor profile data carries repeated m/z positions that the
    /// native default refuses.
    fn picker() -> PeakPickerHiRes {
        PeakPickerHiRes {
            compatibility: PickingCompatibility::source(),
            ..PeakPickerHiRes::default()
        }
    }

    /// One `/proc/self/status` size in kibibytes, or `None` where it is unreadable.
    fn status_kib(key: &str) -> Option<u64> {
        let text = std::fs::read_to_string("/proc/self/status").ok()?;
        text.lines().find_map(|line| {
            line.strip_prefix(key)?
                .split_whitespace()
                .next()?
                .parse()
                .ok()
        })
    }

    /// Print one stage line: elapsed wall time, current RSS and peak RSS.
    fn stage(label: &str, started: Instant) -> Result<()> {
        let mebibytes = |key: &str| {
            status_kib(key).map_or_else(|| "n/a".to_string(), |kib| format!("{}", kib / 1024))
        };
        let mut out = std::io::stdout().lock();
        writeln!(
            out,
            "stage\t{label}\twall_s\t{:.3}\trss_mib\t{}\tpeak_rss_mib\t{}",
            started.elapsed().as_secs_f64(),
            mebibytes("VmRSS:"),
            mebibytes("VmHWM:")
        )
        .map_err(|e| Error::InvalidValue(format!("cannot write the stage report: {e}")))
    }

    /// An order-sensitive FNV-1a digest of every picked coordinate and intensity.
    ///
    /// Bit patterns are hashed, not decimal renderings, so the digest changes on any
    /// change to a picked value, including a sign of zero.
    fn digest(experiment: &MSExperiment) -> u64 {
        let mut hash = 0xcbf2_9ce4_8422_2325_u64;
        let mut eat = |value: u64| {
            for byte in value.to_le_bytes() {
                hash ^= u64::from(byte);
                hash = hash.wrapping_mul(0x0100_0000_01b3);
            }
        };
        eat(experiment.spectra.len() as u64);
        for spectrum in &experiment.spectra {
            eat(spectrum.peaks.len() as u64);
            for peak in &spectrum.peaks {
                eat(peak.mz.to_bits());
                eat(u64::from(peak.intensity.to_bits()));
            }
            for array in &spectrum.float_data_arrays {
                eat(array.data.len() as u64);
                for &value in &array.data {
                    eat(u64::from(value.to_bits()));
                }
            }
        }
        eat(experiment.chromatograms.len() as u64);
        for chromatogram in &experiment.chromatograms {
            eat(chromatogram.peaks.len() as u64);
            for peak in &chromatogram.peaks {
                eat(peak.rt.to_bits());
                eat(u64::from(peak.intensity.to_bits()));
            }
        }
        hash
    }

    /// The number of picked centroids, for a size report next to the digest.
    fn centroids(experiment: &MSExperiment) -> usize {
        experiment.spectra.iter().map(MSSpectrum::len).sum()
    }

    /// Report the largest replicated spectrum count the acquisition ledger admits.
    ///
    /// The records keep this run's acquisition metadata and lose their peaks, so the
    /// search isolates the ledger from the picking itself.
    fn ledger_probe(input: &MSExperiment) -> Result<()> {
        let picker = picker();
        let mut bare: Vec<MSSpectrum> = input.spectra.clone();
        for spectrum in &mut bare {
            spectrum.peaks.clear();
            spectrum.float_data_arrays.clear();
            spectrum.integer_data_arrays.clear();
            spectrum.string_data_arrays.clear();
        }
        if bare.is_empty() {
            return Err(Error::InvalidValue(
                "the run has no spectra to probe".into(),
            ));
        }
        let admits = |count: usize| -> bool {
            let experiment = MSExperiment {
                spectra: (0..count).map(|i| bare[i % bare.len()].clone()).collect(),
                chromatograms: Vec::new(),
                settings: input.settings.clone(),
                sql_run_id: input.sql_run_id,
            };
            picker.pick_experiment(&experiment).is_ok()
        };
        let mut high = 1usize;
        while high <= 8_000_000 && admits(high) {
            high *= 2;
        }
        let mut out = std::io::stdout().lock();
        if high > 8_000_000 {
            return writeln!(out, "ledger\tadmits\t>8000000\tspectra")
                .map_err(|e| Error::InvalidValue(format!("cannot write the ledger report: {e}")));
        }
        let (mut low, mut high) = (high / 2, high);
        while low + 1 < high {
            let middle = low + (high - low) / 2;
            if admits(middle) {
                low = middle;
            } else {
                high = middle;
            }
        }
        writeln!(out, "ledger\tadmits\t{low}\tspectra\trefuses\t{high}")
            .map_err(|e| Error::InvalidValue(format!("cannot write the ledger report: {e}")))
    }

    /// Drive the harness from the command line.
    pub fn run() -> Result<()> {
        let started = Instant::now();
        let mut arguments = std::env::args().skip(1);
        let mut input_path = None;
        let mut output_path = None;
        let (mut in_place, mut probe, mut limit) = (false, false, usize::MAX);
        let mut threads = Threads::serial();
        while let Some(argument) = arguments.next() {
            match argument.as_str() {
                "--out" => output_path = arguments.next(),
                "--in-place" => in_place = true,
                "--ledger-probe" => probe = true,
                "--threads" => {
                    threads = arguments
                        .next()
                        .and_then(|value| value.parse().ok())
                        .map(Threads::from_cli)
                        .ok_or_else(|| Error::InvalidValue("--threads needs a count".into()))?;
                }
                "--limit" => {
                    limit = arguments
                        .next()
                        .and_then(|value| value.parse().ok())
                        .ok_or_else(|| Error::InvalidValue("--limit needs a count".into()))?;
                }
                _ => input_path = Some(argument),
            }
        }
        let input_path =
            input_path.ok_or_else(|| Error::InvalidValue("usage: peak_picking_scale <in.mzML> [--out <out.mzML>] [--in-place] [--ledger-probe] [--limit <n>] [--threads <n>]".into()))?;

        let file = std::fs::File::open(&input_path)
            .map_err(|e| Error::InvalidValue(format!("cannot open {input_path}: {e}")))?;
        let mut experiment = mzml::read_with_options(
            std::io::BufReader::with_capacity(1 << 20, file),
            &read_options(),
        )?;
        if experiment.spectra.len() > limit {
            experiment.spectra.truncate(limit);
        }
        stage("loaded", started)?;
        {
            let mut out = std::io::stdout().lock();
            writeln!(
                out,
                "input\tspectra\t{}\tchromatograms\t{}\tpoints\t{}",
                experiment.spectra.len(),
                experiment.chromatograms.len(),
                centroids(&experiment)
            )
            .map_err(|e| Error::InvalidValue(format!("cannot write the input report: {e}")))?;
        }

        if probe {
            ledger_probe(&experiment)?;
            return stage("probed", started);
        }

        let picker = picker();
        let picked = if in_place {
            picker.pick_experiment_in_place_with_threads(&mut experiment, threads)?;
            experiment
        } else {
            let result = picker.pick_experiment_with_threads(&experiment, threads)?;
            drop(experiment);
            result.experiment
        };
        stage("picked", started)?;
        {
            let mut out = std::io::stdout().lock();
            writeln!(
                out,
                "picked\tcentroids\t{}\tdigest\t{:016x}",
                centroids(&picked),
                digest(&picked)
            )
            .map_err(|e| Error::InvalidValue(format!("cannot write the digest: {e}")))?;
        }

        if let Some(path) = output_path {
            let file = std::fs::File::create(&path)
                .map_err(|e| Error::InvalidValue(format!("cannot create {path}: {e}")))?;
            let mut writer = std::io::BufWriter::with_capacity(1 << 20, file);
            mzml::write_with_options(&mut writer, &picked, &mzml::WriteOptions::default())?;
            writer
                .flush()
                .map_err(|e| Error::InvalidValue(format!("cannot flush {path}: {e}")))?;
            stage("written", started)?;
        }
        Ok(())
    }
}

#[cfg(feature = "mzml")]
fn main() -> openms::Result<()> {
    harness::run()
}

#[cfg(not(feature = "mzml"))]
fn main() -> openms::Result<()> {
    Err(openms::Error::Unsupported(
        "peak_picking_scale requires the mzml feature".into(),
    ))
}
