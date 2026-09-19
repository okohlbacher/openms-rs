#![cfg(feature = "mzml")]
// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! A `sourceFileRef` that names no definition, and the round trip that needs it.
//!
//! Wave 7 gave the streaming [`MSDataWritingConsumer`] a
//! `ReferencePolicy::SourceDangling` that reproduces the `sf_sp_<s>` and
//! `dp_sp_<s>` references source `MzMLHandler::writeSpectrum_` writes for a
//! record the header cannot declare. The reader kept refusing a dangling
//! `sourceFileRef`, so on an input with per-record source files this crate
//! wrote a file it would not read back, while the C++ reader read both its own
//! output and this one. Decision D14 closes that: the reader accepts what the
//! Release build accepts under
//! [`ReadOptions::source_dangling_references`], and the refusal stays as the
//! default strict profile.
//!
//! Evidence: `data/mzml_source_file_round_trip/source_file_refs_oracle.tsv` is
//! the stdout of `../oracle/reader-roundtrip/probe_source_refs.cpp` linked
//! against the Release C++ install at the pins (core `bc9cc12`, cli `c19e494`,
//! topp `174b576`) and run on `ibminode06`
//! (`../oracle/reader-roundtrip/logs/probe_06.log`). It reports, per record,
//! the source file the C++ reader ends up with, whether a scan's and a
//! precursor's `source_file_name`/`source_file_path` metadata exist and what
//! they hold, and which records share one `DataProcessing` allocation — the
//! pointer identity `writeSpectrum_` compares.

use openms::concept::log_stream::{LogColor, LogLevel, LogSink, with_thread_local_log};
use openms::format::ms_data_writing_consumer::{PlainMSDataWritingConsumer, ReferencePolicy};
use openms::format::mzml::{self, AcquisitionMode, ReadOptions};
use openms::kernel::MSExperiment;
use openms::metadata::{MetaInfo, SourceFile};
use openms::{Error, Result};
use std::io::{self, Cursor, Write};
use std::sync::{Arc, Mutex};

const ORACLE: &str = include_str!("data/mzml_source_file_round_trip/source_file_refs_oracle.tsv");
const DANGLING: &str =
    include_str!("data/mzml_source_file_round_trip/record_source_file_dangling.mzML");
const REFS: &str = include_str!("data/peak_picking/PeakPickerHiRes_refs_input.mzML");
/// Four spectra and two chromatograms, with five `sourceFile` definitions
/// declared in the header; untouched by this wave.
const CHROM: &str = include_str!("data/mzml_validator/MzMLFile_1.mzML");
/// The opening tag of that file's first chromatogram, verbatim.
const CHROM_TAG: &str = r#"<chromatogram index="0" id="tic native" defaultArrayLength="15" >"#;
const DUPDP: &str = include_str!("data/peak_picking/PeakPickerHiRes_dupdp_input.mzML");

/// The oracle's cases, in the order `probe_06.sh` passes them to the driver.
const CASES: [(&str, &str); 3] = [("dangling_sf", DANGLING), ("refs", REFS), ("dupdp", DUPDP)];

fn source() -> ReadOptions {
    ReadOptions {
        source_dangling_references: true,
        ..Default::default()
    }
}

/// The oracle comparison reads in [`AcquisitionMode::Source`].
///
/// The default `Canonical` mode collapses a lone `scan` that carries no
/// identifier, no metadata and no meaningful combination
/// (`src/format/mzml_acquisition.rs`, `normalize`), which the `refs` fixture's
/// records have and the source keeps. That normalization is native, predates
/// this group and is documented with the mode itself; reading in `Source` mode
/// takes it out of the comparison so every line below is about references.
fn source_scans() -> ReadOptions {
    ReadOptions {
        acquisition_mode: AcquisitionMode::Source,
        ..source()
    }
}

/// Send this test thread's warning route to a discarding sink, so the expected
/// dangling-reference warnings do not flood the test output. Only
/// `each_distinct_dangling_source_file_reference_warns_once_per_read` inspects
/// them.
fn discard_warnings() {
    with_thread_local_log(LogLevel::Warn, |log| {
        log.remove_all_streams()?;
        log.insert(&LogSink::new(io::sink()))
    })
    .unwrap();
}

fn parse_message(result: Result<impl Sized>) -> String {
    match result {
        Err(Error::Parse { message, .. }) => message,
        Err(other) => panic!("expected a parse error, got {other}"),
        Ok(_) => panic!("expected a parse error, got a result"),
    }
}

fn unsupported_message(result: Result<impl Sized>) -> String {
    match result {
        Err(Error::Unsupported(message)) => message,
        Err(other) => panic!("expected an unsupported error, got {other}"),
        Ok(_) => panic!("expected an unsupported error, got a result"),
    }
}

// Rendering in the oracle driver's line format.

fn text(value: &str) -> &str {
    if value.is_empty() { "<empty>" } else { value }
}

/// A metadata entry the way the driver prints it: `has_X yes/no`, the key's
/// short label, and either the value or `<absent>`.
fn entry(label: &str, value: Option<String>) -> String {
    match value {
        Some(v) => format!(
            "yes\t{label}\t{}",
            if v.is_empty() { "<empty>" } else { &v }
        ),
        None => format!("no\t{label}\t<absent>"),
    }
}

fn oracle_lines(label: &str, experiment: &MSExperiment) -> Vec<String> {
    let mut out = Vec::new();
    out.push(format!(
        "{label}\tLOAD\tok\tspectra\t{}\tchromatograms\t{}\texp_source_files\t{}",
        experiment.spectra.len(),
        experiment.chromatograms.len(),
        experiment.settings.source_files.len()
    ));
    for (i, sf) in experiment.settings.source_files.iter().enumerate() {
        out.push(format!(
            "{label}\tEXPSF\t{i}\tname\t{}\tpath\t{}",
            text(&sf.name),
            text(&sf.path)
        ));
    }
    let default = SourceFile::default();
    for spectrum in &experiment.spectra {
        let sf = &spectrum.source_file;
        out.push(format!(
            "{label}\tSPEC\t{}\tsf_name\t{}\tsf_path\t{}\tsf_checksum\t{}\tsf_type\t{}\tsf_nativeid\t{}\tis_default\t{}",
            spectrum.native_id,
            text(&sf.name),
            text(&sf.path),
            text(&sf.checksum),
            text(&sf.file_type),
            text(&sf.native_id_type),
            if *sf == default { "yes" } else { "no" }
        ));
        for (a, acquisition) in spectrum.acquisition_info.acquisitions.iter().enumerate() {
            let get = |key: &str| acquisition.metadata.get(key).map(ToString::to_string);
            out.push(format!(
                "{label}\tSCAN\t{}\t{a}\thas_name\t{}\thas_path\t{}",
                spectrum.native_id,
                entry("name", get("source_file_name")),
                entry("path", get("source_file_path"))
            ));
        }
        for (p, precursor) in spectrum.precursors.iter().enumerate() {
            let get = |key: &str| {
                precursor
                    .cv_terms
                    .metadata
                    .get(key)
                    .map(ToString::to_string)
            };
            out.push(format!(
                "{label}\tPREC\t{}\t{p}\thas_name\t{}\thas_path\t{}",
                spectrum.native_id,
                entry("name", get("source_file_name")),
                entry("path", get("source_file_path"))
            ));
        }
    }
    for chromatogram in &experiment.chromatograms {
        let sf = &chromatogram.source_file;
        out.push(format!(
            "{label}\tCHROM\t{}\tsf_name\t{}\tsf_path\t{}\tis_default\t{}",
            chromatogram.native_id,
            text(&sf.name),
            text(&sf.path),
            if *sf == default { "yes" } else { "no" }
        ));
    }
    // The pointer identity `writeSpectrum_` compares: which records share one
    // `DataProcessing` allocation. `Arc::ptr_eq` here answers what
    // `shared_ptr::operator==` answers there.
    for (s, spectrum) in experiment.spectra.iter().enumerate() {
        let mut line = format!(
            "{label}\tDPPTR\t{}\tsize\t{}",
            spectrum.native_id,
            spectrum.data_processing.len()
        );
        for (j, entry) in spectrum.data_processing.iter().enumerate() {
            let mut shared: i64 = -1;
            for t in 0..s {
                let earlier = &experiment.spectra[t].data_processing;
                if earlier.get(j).is_some_and(|e| Arc::ptr_eq(e, entry)) {
                    shared = t as i64;
                    break;
                }
            }
            line.push_str(&format!("\tslot{j}_same_as_record\t{shared}"));
        }
        let first = &experiment.spectra[0].data_processing;
        let same = spectrum.data_processing.len() == first.len()
            && std::iter::zip(&spectrum.data_processing, first).all(|(a, b)| Arc::ptr_eq(a, b));
        line.push_str(&format!(
            "\tequals_record0\t{}",
            if s == 0 {
                "self"
            } else if same {
                "yes"
            } else {
                "no"
            }
        ));
        out.push(line);
    }
    out
}

/// The port reads each of the three inputs into exactly what the executed
/// Release C++ reads them into, line for line.
#[test]
fn source_option_reproduces_the_executed_cpp_oracle() {
    discard_warnings();
    let mut rust = Vec::new();
    for (label, xml) in CASES {
        let experiment = mzml::read_with_options(Cursor::new(xml), &source_scans()).unwrap();
        rust.extend(oracle_lines(label, &experiment));
    }
    let oracle: Vec<&str> = ORACLE.lines().collect();
    assert_eq!(rust.len(), oracle.len(), "{rust:#?}");
    for (line, (rust, oracle)) in rust.iter().zip(&oracle).enumerate() {
        assert_eq!(rust, oracle, "oracle line {}", line + 1);
    }
}

/// The strict default keeps refusing every shape of dangling `sourceFileRef`:
/// on a spectrum, on a scan and on a precursor.
///
/// This is the native profile the source has no equivalent of, kept because a
/// caller that did not ask for the loss should not receive a record whose
/// source file silently became empty.
#[test]
fn strict_default_refuses_every_dangling_source_file_reference() {
    discard_warnings();
    assert!(!ReadOptions::default().source_dangling_references);
    // The whole fixture: the spectrum reference is the first one reached.
    assert_eq!(
        parse_message(mzml::read(Cursor::new(DANGLING))),
        "unresolved sourceFileRef"
    );
    // Each shape on its own, so no earlier refusal can mask a later one.
    for (what, keep) in [
        ("spectrum", "scan=2"),
        ("scan", "scan=3"),
        ("precursor", "scan=4"),
    ] {
        let only = one_record(DANGLING, keep);
        assert_ne!(only, DANGLING, "{what}");
        assert_eq!(
            parse_message(mzml::read(Cursor::new(&only))),
            "unresolved sourceFileRef",
            "{what}"
        );
        // And the same document reads under the source policy, into the value
        // the executed C++ leaves behind for that shape.
        let experiment = mzml::read_with_options(Cursor::new(&only), &source()).unwrap();
        assert_eq!(experiment.spectra.len(), 1, "{what}");
        let spectrum = &experiment.spectra[0];
        assert_eq!(spectrum.source_file, SourceFile::default(), "{what}");
        let carried: Vec<&MetaInfo> = spectrum
            .acquisition_info
            .acquisitions
            .iter()
            .map(|a| &a.metadata)
            .chain(spectrum.precursors.iter().map(|p| &p.cv_terms.metadata))
            .collect();
        match what {
            // The spectrum's own reference leaves no metadata anywhere: the
            // record simply keeps a default `SourceFile`.
            "spectrum" => assert!(carried.is_empty(), "{carried:?}"),
            // A scan's and a precursor's set both keys, to empty strings.
            _ => {
                assert_eq!(carried.len(), 1, "{what}: {carried:?}");
                assert_eq!(carried[0].len(), 2, "{what}: {carried:?}");
                assert_eq!(carried[0]["source_file_name"].to_string(), "", "{what}");
                assert_eq!(carried[0]["source_file_path"].to_string(), "", "{what}");
            }
        }
    }
}

/// A chromatogram's `sourceFileRef` is refused before any policy is consulted,
/// and refused whether the reference resolves or dangles.
///
/// mzML 1.1 puts `sourceFileRef` on `SpectrumType`, `ScanType` and
/// `PrecursorType` and **not** on `ChromatogramType`
/// (`share/OpenMS/SCHEMAS/mzML_1_10.xsd` of the pinned core `bc9cc12`, whose
/// `ChromatogramType` carries only `id`, `index`, `defaultArrayLength` and
/// `dataProcessingRef`), so a document with one is not mzML at all and there is
/// no reference for a dangling-reference policy to have an opinion about.
///
/// This is a **divergence from the Release reader**, which is why it is pinned
/// here rather than left to the prose that asserts it. `MzMLHandler.cpp:937-941`
/// reads the attribute with an unguarded
/// `chromatogram_.setSourceFile(source_files_[source_file_ref]);`. The
/// spectrum arm at `:896-906` asks
/// `if (source_files_.contains(source_file_ref))` first and warns when the
/// answer is no. The chromatogram arm does neither, so the C++ accepts the
/// document and silently default-constructs an empty `SourceFile` for an id it
/// has never seen. Three places in this wave describe that difference as
/// measured — the [`ReadOptions::source_dangling_references`] rustdoc, the
/// header reader and the deferred list. `mzml_header.rs`'s
/// `chromatogram_source_file_is_rejected_instead_of_emitting_an_invalid_attribute`
/// already executed the branch on the writing side; this test executes it on
/// the reading side, under both policies, which is what those three places
/// describe.
///
/// Neither writer in this crate emits the attribute, so the round trip does
/// not reach this; it is the boundary beside the round trip, not inside it.
#[test]
fn a_chromatograms_source_file_ref_is_refused_under_both_policies() {
    discard_warnings();
    // The fixture itself reads, so the refusal below is about the attribute.
    let clean = mzml::read(Cursor::new(CHROM)).unwrap();
    assert_eq!(clean.spectra.len(), 4);
    assert_eq!(clean.chromatograms.len(), 2);

    for (what, id) in [
        // `sf1` IS declared in this file's `sourceFileList`, so this document
        // has nothing dangling in it at all and is still refused.
        ("a reference the header declares", "sf1"),
        ("a reference nothing declares", "sf_not_declared"),
    ] {
        let doc = CHROM.replacen(
            CHROM_TAG,
            &CHROM_TAG.replace(" >", &format!(" sourceFileRef=\"{id}\">")),
            1,
        );
        assert_ne!(doc, CHROM, "{what}: the chromatogram tag was not rewritten");
        assert!(doc.contains(&format!("sourceFileRef=\"{id}\"")), "{what}");

        for (policy, options) in [
            ("the strict default", ReadOptions::default()),
            ("source_dangling_references", source()),
        ] {
            assert_eq!(
                unsupported_message(mzml::read_with_options(Cursor::new(&doc), &options)),
                "chromatogram sourceFileRef is not permitted by mzML",
                "{what}, under {policy}"
            );
        }
    }
}

/// Drop every `spectrum` element except the one whose id is `keep`, and
/// correct the declared count.
fn one_record(xml: &str, keep: &str) -> String {
    let mut out = String::new();
    let mut rest = xml;
    while let Some(at) = rest.find("<spectrum ") {
        out.push_str(&rest[..at]);
        let body = &rest[at..];
        let tag_end = body.find('>').expect("a start tag") + 1;
        let end = if body[..tag_end].ends_with("/>") {
            tag_end
        } else {
            body.find("</spectrum>").expect("a closing tag") + "</spectrum>".len()
        };
        let record = &body[..end];
        if record.contains(&format!("id=\"{keep}\"")) {
            out.push_str(record);
        }
        rest = &body[end..];
    }
    out.push_str(rest);
    out.replace("<spectrumList count=\"5\"", "<spectrumList count=\"1\"")
}

#[derive(Clone, Default)]
struct Capture(Arc<Mutex<Vec<u8>>>);
impl Write for Capture {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(buffer);
        Ok(buffer.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// The port warns once per distinct dangling ID per read, where the source
/// warns once per *occurrence* and only for a spectrum's reference.
///
/// Measured on the same fixture: the executed C++ reports
/// `Error: unregistered source file reference sf_sp_1.` twice, once for each of
/// the two spectra naming it, and says nothing at all about `sf_sp_2`, which
/// only a scan and a precursor name (`../oracle/reader-roundtrip/logs/probe_06.log`,
/// the `dangling_sf` stderr). The port reports each of the two IDs once. The
/// warning is native either way; it exists so a caller that opted into the loss
/// can still see it, and the read result does not depend on it.
#[test]
fn each_distinct_dangling_source_file_reference_warns_once_per_read() {
    let capture = Capture::default();
    let sink = LogSink::new(capture.clone());
    with_thread_local_log(LogLevel::Warn, |log| {
        log.remove_all_streams()?;
        log.set_color(None);
        log.insert(&sink)
    })
    .unwrap();
    let take = || {
        with_thread_local_log(LogLevel::Warn, |log| log.clear_cache()).unwrap();
        String::from_utf8(std::mem::take(&mut *capture.0.lock().unwrap())).unwrap()
    };

    mzml::read_with_options(Cursor::new(DANGLING), &source()).unwrap();
    assert_eq!(
        take(),
        "Warning: mzML sourceFileRef 'sf_sp_1' names no definition; source-compatible reading uses an empty source file.\n\
         Warning: mzML sourceFileRef 'sf_sp_2' names no definition; source-compatible reading uses an empty source file.\n"
    );
    // A second read starts from a fresh record, and a strict read says nothing
    // because it fails instead.
    mzml::read_with_options(Cursor::new(DANGLING), &source()).unwrap();
    assert_eq!(take().lines().count(), 2);
    assert!(mzml::read(Cursor::new(DANGLING)).is_err());
    assert_eq!(take(), "");

    with_thread_local_log(LogLevel::Warn, |log| {
        log.remove_all_streams()?;
        log.set_color(Some(LogColor::Yellow));
        log.insert(&LogSink::stderr())
    })
    .unwrap();
}

/// What this lane exists for: the crate reads back the low-memory file it
/// writes itself.
///
/// Every record of the `refs` fixture but the first carries a `sourceFileRef`,
/// so `ReferencePolicy::SourceDangling` renumbers four of them into `sf_sp_1`
/// through `sf_sp_4`, which the streamed header does not declare. Before
/// decision D14 the crate's own reader refused that file — the round trip this
/// test performs failed at the second step while the C++ reader completed it.
#[test]
fn the_low_memory_output_of_a_per_record_source_file_input_reads_back() {
    discard_warnings();
    let input = mzml::read_with_options(Cursor::new(REFS), &source()).unwrap();
    assert!(input.spectra.iter().all(|s| !s.source_file.name.is_empty()));

    let mut consumer = PlainMSDataWritingConsumer::plain(Vec::new())
        .with_reference_policy(ReferencePolicy::SourceDangling);
    consumer.set_experimental_settings(&input.settings).unwrap();
    consumer.set_expected_size(input.spectra.len(), 0).unwrap();
    for spectrum in &input.spectra {
        consumer.consume_spectrum(&mut spectrum.clone()).unwrap();
    }
    let written = consumer.finish().unwrap();
    let text = String::from_utf8(written.clone()).unwrap();
    for id in ["sf_sp_1", "sf_sp_2", "sf_sp_3", "sf_sp_4"] {
        assert!(text.contains(&format!("sourceFileRef=\"{id}\"")), "{id}");
        assert!(!text.contains(&format!("<sourceFile id=\"{id}\"")), "{id}");
    }

    // The strict profile still refuses it, which is the whole reason the file
    // is worth a test: the references really do name nothing. Records 1 and 2
    // carry a renumbered `dataProcessingRef` as well, and that one is reached
    // first, so the message is whichever the reader meets first.
    let refused = parse_message(mzml::read(Cursor::new(written.clone())));
    assert!(
        matches!(
            refused.as_str(),
            "unresolved sourceFileRef" | "unresolved dataProcessingRef"
        ),
        "{refused}"
    );
    // With the processing references removed by hand, the source-file one is
    // what is left, and it is refused on its own.
    let only_source_files = text
        .split(" dataProcessingRef=\"dp_sp_")
        .enumerate()
        .map(|(i, part)| {
            if i == 0 {
                part.to_owned()
            } else {
                part[part.find('"').map_or(0, |q| q + 1)..].to_owned()
            }
        })
        .collect::<String>();
    assert_ne!(only_source_files, text);
    assert_eq!(
        parse_message(mzml::read(Cursor::new(&only_source_files))),
        "unresolved sourceFileRef"
    );

    // The source profile reads it, and the peak data survives the round trip
    // unchanged (decision D6: decoded content, not bytes).
    let back = mzml::read_with_options(Cursor::new(written), &source()).unwrap();
    assert_eq!(back.spectra.len(), input.spectra.len());
    for (before, after) in std::iter::zip(&input.spectra, &back.spectra) {
        assert_eq!(before.native_id, after.native_id);
        assert_eq!(before.peaks.len(), after.peaks.len());
        for (a, b) in std::iter::zip(&before.peaks, &after.peaks) {
            assert_eq!(a.mz.to_bits(), b.mz.to_bits());
            assert_eq!(a.intensity.to_bits(), b.intensity.to_bits());
        }
    }
    // The first record's source file is the one the header declares and
    // survives; the four renumbered references named nothing and are dropped,
    // which is the loss `CPP-172` describes and this crate now reproduces on
    // both sides of the round trip rather than on the writing side alone.
    assert_eq!(back.spectra[0].source_file.name, "part_one.mzML");
    assert!(
        back.spectra[1..]
            .iter()
            .all(|s| s.source_file == SourceFile::default())
    );
}

/// The whole-document writer still decides by content, and this pins it.
///
/// `MzMLHandler::writeHeader_` deduplicates histories by **content**
/// (`Helpers::cmpPtrContainer`, `MzMLHandler.cpp:5060-5069`) while
/// `writeSpectrum_` compares them by **pointer** (`:5258`), so the C++
/// whole-document `MzMLFile::store` dangles on an input carrying two
/// textually identical `dataProcessing` entries under different identifiers —
/// the defect `CPP-172` records, in the path that is not the streaming
/// consumer. Measured on `ibminode06`
/// (`../oracle/reader-roundtrip/logs/roundtrip_06.log`, section `dupdp`): the
/// C++ in-memory output of `PeakPickerHiRes_dupdp_input` declares ten
/// referenceable ids against eleven references, **two** of them dangling,
/// `dp_sp_1` and `dp_sp_2`; this crate's declares ten against eight, none
/// dangling, because it deduplicates by content on both sides and records 1
/// and 2 fall back to the list default.
///
/// The streaming consumer's half of this is closed —
/// `ReferencePolicy::SourceDangling` now compares by pointer, and the two
/// low-memory outputs carry the same two dangling identifiers. Closing the
/// whole-document half would make `mzml::write` emit references that mzML
/// forbids on every write path, by default and with nothing to opt into,
/// which is a change to the crate's default output validity rather than a
/// fidelity fix inside one policy. Recorded as a divergence and left for the
/// lead; the decoded content is identical either way and each implementation
/// reads the other's file.
#[test]
fn the_whole_document_writer_still_deduplicates_by_content() {
    discard_warnings();
    let experiment = mzml::read_with_options(Cursor::new(DUPDP), &source()).unwrap();
    // The input really is the pointer-vs-content case: records 0 and 1 hold
    // equal but distinct histories.
    let (first, second) = (
        &experiment.spectra[0].data_processing,
        &experiment.spectra[1].data_processing,
    );
    assert_eq!(first, second);
    assert!(!Arc::ptr_eq(&first[0], &second[0]));

    let mut written = Vec::new();
    mzml::write(&mut written, &experiment).unwrap();
    let text = String::from_utf8(written).unwrap();
    // No record carries a reference the document does not declare: the port's
    // whole-document output stays valid mzML where the source's does not.
    let declared: Vec<&str> = text
        .match_indices("<dataProcessing id=\"")
        .map(|(at, pattern)| {
            let rest = &text[at + pattern.len()..];
            &rest[..rest.find('"').expect("a quoted id")]
        })
        .collect();
    let referenced: Vec<&str> = text
        .match_indices("dataProcessingRef=\"")
        .map(|(at, pattern)| {
            let rest = &text[at + pattern.len()..];
            &rest[..rest.find('"').expect("a quoted id")]
        })
        .collect();
    for id in &referenced {
        assert!(declared.contains(id), "{id} dangles: {declared:?}");
    }
    // The five content-equal histories collapse to one declaration and no
    // record carries a reference at all, where the source declares the same
    // one and still writes `dp_sp_1` on record 1 and `dp_sp_2` on record 2.
    assert_eq!(declared.len(), 1, "{declared:?}");
    assert!(referenced.is_empty(), "{referenced:?}");
    assert!(!text.contains("dp_sp_1"), "{text:.400}");
    assert!(!text.contains("dp_sp_2"), "{text:.400}");
}
