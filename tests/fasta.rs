// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use openms::concept::progress_logger::{
    ProgressBackend, ProgressLogger, ProgressNesting, ProgressTime,
};
use openms::format::fasta::{self, FASTAEntry, FASTAFile, FastaOptions, FastaReader, FastaWriter};
use openms::system::file::TempFile;
use openms::{Error, Result};
use std::io::{self, BufRead, BufReader, Cursor, Read, Seek, SeekFrom, Write};
use std::sync::{
    Arc, Mutex,
    atomic::{AtomicI64, AtomicUsize, Ordering},
};

fn entry(id: &str, desc: &str, seq: &str) -> FASTAEntry {
    FASTAEntry::new(id, desc, seq)
}
fn fixture() -> &'static [u8] {
    include_bytes!("data/fasta_source.fasta")
}

#[test]
fn all_fourteen_source_field_literals_and_modified_chemistry() {
    let records = fasta::read(fixture()).unwrap();
    assert_eq!(records.len(), 5);
    for row in include_str!("data/fasta_source_fields.tsv").lines().skip(1) {
        let fields: Vec<_> = row.splitn(4, '\t').collect();
        let record = &records[fields[0].parse::<usize>().unwrap()];
        let actual = match fields[1] {
            "identifier" => &record.identifier,
            "description" => &record.description,
            "sequence" => &record.sequence,
            _ => unreachable!(),
        };
        assert_eq!(actual, fields[3], "source line {}", fields[2]);
    }
    let modified = openms::chemistry::AASequence::parse(&records[3].sequence).unwrap();
    assert!(modified.n_terminal_modification().is_some());
    assert_eq!(modified.as_str(), records[2].sequence);
    let mut bytes = Vec::new();
    fasta::write(&mut bytes, &records).unwrap();
    assert_eq!(fasta::read(bytes.as_slice()).unwrap(), records);
}

#[test]
fn entry_construction_equality_and_source_comparison_helpers() {
    assert_eq!(FASTAEntry::default(), entry("", "", ""));
    let original = entry("id", "desc", "SEQ");
    let different_sequence = entry("id", "desc", "OTHER");
    let different_header = entry("id2", "desc", "SEQ");
    assert!(original.header_matches(&different_sequence));
    assert!(!original.sequence_matches(&different_sequence));
    assert!(!original.header_matches(&different_header));
    assert!(original.sequence_matches(&different_header));
    assert_eq!(original.clone(), original);
}

#[test]
fn literal_peff_whitespace_asterisks_and_opaque_sequence_rules() {
    assert_eq!(
        fasta::read(b"# PEFF 1.0\n  # comment\n\t>TEST_HEADER Test description\nSEQ\n".as_slice())
            .unwrap(),
        [entry("TEST_HEADER", "Test description", "SEQ")]
    );
    assert_eq!(
        fasta::read(b">Header1\n\nSEQ1\n\n>Header2\n\nSEQ2\n\n".as_slice()).unwrap(),
        [entry("Header1", "", "SEQ1"), entry("Header2", "", "SEQ2")]
    );
    let value = "GDREQLLQRAR*LAEQ*AERYDDMASAMKAVTEL";
    assert_eq!(
        fasta::read(format!(">ID\n{value}\n").as_bytes()).unwrap()[0].sequence,
        value
    );
    let records = fasta::read(
        ">\t id  leading\t tabs\r gone\nA 1\t(ICPL:13C(6))\ré;\u{b}\u{c}\u{a0}\n;kept\n".as_bytes(),
    )
    .unwrap();
    assert_eq!(
        records,
        [entry(
            "id",
            " leading tabs gone",
            "A1(ICPL:13C(6))é;\u{b}\u{c}\u{a0};kept"
        )]
    );
    let mut bytes = Vec::new();
    fasta::write(
        &mut bytes,
        &[entry("id", "", "GDREQLLQRAR LAEQ\tAERYDDMASAMKAVTEL")],
    )
    .unwrap();
    assert_eq!(
        fasta::read(bytes.as_slice()).unwrap()[0].sequence,
        "GDREQLLQRARLAEQAERYDDMASAMKAVTEL"
    );
}

#[test]
fn exact_header_boundary_and_eof_state_distinctions() {
    assert_eq!(
        fasta::read(b">empty\n>next\nABC".as_slice()).unwrap(),
        [entry("empty", "", ">nextABC")]
    );
    assert_eq!(
        fasta::read(b">a\nA\n >b\nB".as_slice()).unwrap(),
        [entry("a", "", "A>bB")]
    );
    assert_eq!(
        fasta::read(b">a\nA\r>b\rB".as_slice()).unwrap(),
        [entry("a", "", "A>bB")]
    );
    assert!(
        fasta::read(b"# trailing comment".as_slice())
            .unwrap()
            .is_empty()
    );
    for text in [
        "",
        " \t\r\n",
        "# comment\n",
        ">a",
        ">a desc",
        ">a\n",
        ">\nA",
        ";not a prologue comment\n>a\nA",
    ] {
        assert!(fasta::read(text.as_bytes()).is_err(), "{text:?}");
    }
    let mut reader = FastaReader::new(Cursor::new(b""));
    assert!(reader.at_end().unwrap());
    assert_eq!(reader.position().unwrap(), None);
    assert!(reader.next_entry().unwrap().is_none());
    assert!(fasta::read(b">id\n\xff\n".as_slice()).is_err());
}

#[test]
fn source_positions_seek_replay_and_native_destination_atomicity() {
    let mut reader = FastaReader::new(Cursor::new(fixture()));
    let first_position = reader.position().unwrap().unwrap();
    let mut records = Vec::new();
    let mut positions = Vec::new();
    for _ in 0..4 {
        records.push(reader.next_entry().unwrap().unwrap());
        positions.push(reader.position().unwrap().unwrap());
    }
    assert!(reader.set_position(positions[0]).unwrap());
    for i in 1..4 {
        assert_eq!(reader.next_entry().unwrap().unwrap(), records[i]);
        assert_eq!(reader.position().unwrap(), Some(positions[i]));
    }
    reader.next_entry().unwrap().unwrap();
    assert_eq!(reader.position().unwrap(), None);
    assert!(reader.at_end().unwrap());
    assert!(reader.set_position(positions[0]).unwrap());
    assert!(!reader.set_position(fixture().len() as u64 + 1).unwrap());
    assert_eq!(reader.position().unwrap(), Some(positions[0]));
    assert_eq!(reader.next_entry().unwrap().unwrap(), records[1]);
    assert!(reader.set_position(fixture().len() as u64).unwrap());
    assert_eq!(reader.position().unwrap(), Some(fixture().len() as u64));
    let mut output = entry("keep", "me", "UNCHANGED");
    assert!(reader.read_next(&mut output).is_err());
    assert_eq!(output, entry("keep", "me", "UNCHANGED"));
    assert!(!reader.read_next(&mut output).unwrap());
    assert!(reader.set_position(first_position).unwrap());
    assert!(reader.read_next(&mut output).unwrap());
    let mut peff = FastaReader::new(Cursor::new(b"#header\n>a\nA"));
    assert_eq!(peff.position().unwrap(), Some(8));
    assert!(peff.set_position(0).unwrap());
    assert!(peff.next_entry().is_err()); // seeking does not repeat PEFF initialization
}

struct SeekFailure(Cursor<Vec<u8>>);
impl Read for SeekFailure {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        self.0.read(bytes)
    }
}
impl BufRead for SeekFailure {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.0.fill_buf()
    }
    fn consume(&mut self, count: usize) {
        self.0.consume(count);
    }
}
impl Seek for SeekFailure {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        if position == SeekFrom::Start(1) {
            Err(io::Error::other("seek failed"))
        } else {
            self.0.seek(position)
        }
    }
}
#[test]
fn failed_seek_is_reported_and_replayed_reads_spend_the_same_limits() {
    let mut reader = FastaReader::new(SeekFailure(Cursor::new(b">a\nA".to_vec())));
    assert!(matches!(reader.set_position(1), Err(Error::Io(_))));
    let options = FastaOptions {
        max_records: 1,
        ..Default::default()
    };
    let mut reader = FastaReader::with_options(Cursor::new(b">a\nA"), options);
    assert!(reader.next_entry().unwrap().is_some());
    reader.set_position(0).unwrap();
    assert!(
        reader
            .next_entry()
            .unwrap_err()
            .to_string()
            .contains("record limit")
    );
}

struct Endless {
    consumed: Arc<AtomicUsize>,
}
impl Read for Endless {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        bytes.fill(b' ');
        self.consumed.fetch_add(bytes.len(), Ordering::Relaxed);
        Ok(bytes.len())
    }
}
impl BufRead for Endless {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        Ok(b" ")
    }
    fn consume(&mut self, count: usize) {
        self.consumed.fetch_add(count, Ordering::Relaxed);
    }
}
#[test]
fn streaming_bounds_apply_before_unbounded_lines_or_records_allocate() {
    let consumed = Arc::new(AtomicUsize::new(0));
    let options = FastaOptions {
        max_input_bytes: 12,
        ..Default::default()
    };
    let mut reader = FastaReader::with_options(
        Endless {
            consumed: consumed.clone(),
        },
        options,
    );
    assert!(
        reader
            .next_entry()
            .unwrap_err()
            .to_string()
            .contains("input byte")
    );
    assert_eq!(consumed.load(Ordering::Relaxed), 12);
    for (input, options, message) in [
        (
            b">id\nABC".as_slice(),
            FastaOptions {
                max_entry_bytes: 4,
                ..Default::default()
            },
            "entry byte",
        ),
        (
            b">id\nABC".as_slice(),
            FastaOptions {
                max_work: 3,
                ..Default::default()
            },
            "work limit",
        ),
        (
            b">id\nABC\n>b\nDEF".as_slice(),
            FastaOptions {
                max_records: 1,
                ..Default::default()
            },
            "record limit",
        ),
    ] {
        assert!(
            fasta::read_with_options(input, options)
                .unwrap_err()
                .to_string()
                .contains(message)
        );
    }
    assert_eq!(
        fasta::read(BufReader::with_capacity(1, fixture()))
            .unwrap()
            .len(),
        5
    );
}

#[test]
fn streamed_writes_exact_bytes_empty_sequences_and_preflight() {
    let mut writer = FastaWriter::new(Vec::new());
    writer.write_entry(&entry("empty", "", "")).unwrap();
    writer
        .write_entry(&entry("a", " d", &"A".repeat(81)))
        .unwrap();
    let bytes = writer.finish().unwrap();
    assert_eq!(
        bytes,
        format!(">empty \n>a  d\n{}\nA\n", "A".repeat(80)).into_bytes()
    );
    let records = [entry("a", "", "A"), entry("b", "bad\nheader", "B")];
    let mut bytes = Vec::new();
    assert!(fasta::write(&mut bytes, &records).is_err());
    assert!(bytes.is_empty());
    for options in [
        FastaOptions {
            max_output_bytes: 5,
            ..Default::default()
        },
        FastaOptions {
            max_entry_bytes: 1,
            ..Default::default()
        },
        FastaOptions {
            max_records: 0,
            ..Default::default()
        },
        FastaOptions {
            max_work: 4,
            ..Default::default()
        },
    ] {
        assert!(fasta::write_with_options(&mut bytes, &records[..1], options).is_err());
        assert!(bytes.is_empty());
    }
    let mut writer = FastaWriter::with_options(
        Vec::new(),
        FastaOptions {
            max_records: 1,
            ..Default::default()
        },
    );
    writer.write_entry(&records[0]).unwrap();
    assert!(writer.write_entry(&records[0]).is_err());
    assert_eq!(writer.finish().unwrap(), b">a \nA\n");
    // Byte wrapping can split an encoded character; the parser joins before UTF-8 validation.
    let unicode = entry("utf8", "", &format!("{}é", "A".repeat(79)));
    let mut bytes = Vec::new();
    fasta::write(&mut bytes, std::slice::from_ref(&unicode)).unwrap();
    assert!(std::str::from_utf8(&bytes).is_err());
    assert_eq!(fasta::read(bytes.as_slice()).unwrap(), [unicode]);
}

#[derive(Clone)]
struct FailingWriter {
    writes: Arc<AtomicUsize>,
    flushes: Arc<AtomicUsize>,
    fail_write: bool,
}
impl Write for FailingWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.writes.fetch_add(1, Ordering::Relaxed);
        if self.fail_write {
            Err(io::Error::other("write"))
        } else {
            Ok(bytes.len())
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        self.flushes.fetch_add(1, Ordering::Relaxed);
        Err(io::Error::other("flush"))
    }
}
#[test]
fn writer_io_errors_are_checked_and_no_destructor_retry_is_added() {
    for fail_write in [false, true] {
        let writes = Arc::new(AtomicUsize::new(0));
        let flushes = Arc::new(AtomicUsize::new(0));
        let mut writer = FastaWriter::new(FailingWriter {
            writes: writes.clone(),
            flushes: flushes.clone(),
            fail_write,
        });
        if fail_write {
            assert!(writer.write_entry(&entry("a", "", "A")).is_err());
            assert!(writer.write_entry(&entry("a", "", "A")).is_err());
            assert_eq!(writes.load(Ordering::Relaxed), 1);
            assert!(writer.finish().is_err());
            assert_eq!(flushes.load(Ordering::Relaxed), 0);
        } else {
            writer.write_entry(&entry("a", "", "A")).unwrap();
            assert!(matches!(writer.finish(), Err(Error::Io(_))));
            assert_eq!(flushes.load(Ordering::Relaxed), 1);
        }
    }
}

#[test]
fn file_lifecycle_simultaneous_sessions_extension_and_atomic_load() {
    let input = TempFile::new().unwrap();
    let output = TempFile::new().unwrap();
    let other = TempFile::new().unwrap();
    std::fs::write(input.path(), fixture()).unwrap();
    let expected = fasta::read(fixture()).unwrap();
    let mut file = FASTAFile::new();
    assert!(file.position().is_err());
    assert!(file.at_end().is_err());
    assert!(file.read_next(&mut FASTAEntry::default()).is_err());
    assert!(file.write_next(&expected[0]).is_err());
    file.write_end().unwrap();
    file.read_start(input.path()).unwrap();
    file.write_start(output.path()).unwrap();
    assert!(file.write_start(other.path()).is_err());
    let mut record = FASTAEntry::default();
    assert!(file.read_next(&mut record).unwrap());
    file.write_next(&record).unwrap();
    let position = file.position().unwrap();
    file.store(other.path(), &expected).unwrap();
    assert_eq!(file.load(other.path()).unwrap(), expected);
    assert_eq!(file.position().unwrap(), position);
    while file.read_next(&mut record).unwrap() {
        file.write_next(&record).unwrap();
    }
    file.write_end().unwrap();
    file.write_end().unwrap();
    assert_eq!(file.load(output.path()).unwrap(), expected);
    std::fs::write(other.path(), b">bad").unwrap();
    let mut destination = expected.clone();
    assert!(file.load_into(other.path(), &mut destination).is_err());
    assert_eq!(destination, expected);
    let disallowed = other.path().with_extension("mzML");
    assert!(file.write_start(&disallowed).is_err());
    assert!(!disallowed.exists());
    let unknown = TempFile::new().unwrap();
    file.store(unknown.path(), &expected).unwrap();
    let bad = entry("bad\nheader", "", "A");
    assert!(file.store(unknown.path(), &[bad]).is_err());
    assert_eq!(file.load(unknown.path()).unwrap(), expected);
}

struct Recorder(Arc<Mutex<Vec<String>>>);
impl ProgressBackend for Recorder {
    fn start_progress(&mut self, begin: i64, end: i64, label: &str, _: usize) -> Result<()> {
        self.0
            .lock()
            .unwrap()
            .push(format!("start:{begin}:{end}:{label}"));
        Ok(())
    }
    fn set_progress(&mut self, value: i64, _: usize) -> Result<()> {
        self.0.lock().unwrap().push(format!("set:{value}"));
        Ok(())
    }
    fn next_progress(&mut self) -> Result<i64> {
        self.0.lock().unwrap().push("next".into());
        Ok(1)
    }
    fn end_progress(&mut self, _: usize, bytes: u64) -> Result<()> {
        self.0.lock().unwrap().push(format!("end:{bytes}"));
        Ok(())
    }
}
#[test]
fn source_progress_wrappers_duplicate_updates_and_aggregate_labels() {
    let events = Arc::new(Mutex::new(Vec::new()));
    let tick = Arc::new(AtomicI64::new(0));
    let mut progress = ProgressLogger::with_clock_and_nesting(
        Arc::new(move || {
            Ok(ProgressTime {
                wall_second: tick.fetch_add(1, Ordering::Relaxed),
                wall_seconds: 0.0,
                cpu_seconds: None,
            })
        }),
        ProgressNesting::default(),
    );
    progress.set_logger(Box::new(Recorder(events.clone())));
    let mut file = FASTAFile::new();
    file.progress = progress;
    let input = TempFile::new().unwrap();
    std::fs::write(input.path(), b">a\nA\n>b\nB").unwrap();
    file.read_start_with_progress(input.path(), "Proteins")
        .unwrap();
    let mut record = FASTAEntry::default();
    assert!(file.read_next_with_progress(&mut record).unwrap());
    assert!(file.read_next_with_progress(&mut record).unwrap());
    assert!(!file.read_next_with_progress(&mut record).unwrap());
    assert_eq!(
        *events.lock().unwrap(),
        [
            "start:0:9:Proteins",
            "set:5",
            "set:5",
            "set:-1",
            "set:-1",
            "end:0"
        ]
    );
    events.lock().unwrap().clear();
    let entries = file.load(input.path()).unwrap();
    assert_eq!(
        *events.lock().unwrap(),
        ["start:0:1:Loading FASTA file", "end:0"]
    );
    events.lock().unwrap().clear();
    file.store(input.path(), &entries).unwrap();
    assert_eq!(
        *events.lock().unwrap(),
        [
            "start:0:2:Writing FASTA file",
            "next",
            "set:1",
            "next",
            "set:1",
            "end:0"
        ]
    );
}
