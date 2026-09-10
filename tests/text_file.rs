// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::format::{
    TextFile,
    text::{Limits, ReadOptions},
};
use std::io::{BufReader, Cursor, Read, Write};

fn read(text: impl AsRef<[u8]>, options: &ReadOptions) -> TextFile {
    TextFile::from_reader(Cursor::new(text.as_ref()), options).unwrap()
}
#[test]
fn all_source_text_fixture_positions_and_counts() {
    let source = include_bytes!("data/text_file_source.txt");
    let raw = read(source, &ReadOptions::default());
    assert_eq!(raw.len(), 11);
    assert_eq!(raw.lines()[0], "first_line");
    assert_eq!(raw.lines()[3], "middle_line");
    assert_eq!(raw.lines()[10], "last_line");
    assert!(
        read(
            include_bytes!("data/text_file_empty_source.txt"),
            &ReadOptions::default()
        )
        .is_empty()
    );
    let options = ReadOptions {
        trim_lines: true,
        ..Default::default()
    };
    let trimmed = read(source, &options);
    assert_eq!(
        trimmed.lines(),
        [
            "first_line",
            "",
            "",
            "middle_line",
            "",
            "space_line",
            "tab_line",
            "back_space_line",
            "back_tab_line",
            "",
            "last_line"
        ]
    );
    for n in [1, 3, 4] {
        let first = read(
            source,
            &ReadOptions {
                first_n: n,
                ..options.clone()
            },
        );
        assert_eq!(first.lines(), &trimmed.lines()[..n as usize]);
    }
    let skip = ReadOptions {
        skip_empty_lines: true,
        ..options
    };
    let retained = read(source, &skip);
    assert_eq!(
        retained.lines(),
        [
            "first_line",
            "middle_line",
            "space_line",
            "tab_line",
            "back_space_line",
            "back_tab_line",
            "last_line"
        ]
    );
    assert_eq!(
        read(source, &ReadOptions { first_n: 4, ..skip }).lines(),
        &retained.lines()[..4]
    );
}
#[test]
fn line_endings_final_line_and_source_getline_are_exact() {
    for capacity in [1, 2, 3, 64] {
        let input = b"a\rb\r\nc\n\rd\r\r\nlast";
        let mut reader = BufReader::with_capacity(capacity, Cursor::new(input));
        let mut line = String::from("old");
        for expected in ["a", "b", "c", "", "d", "", "last"] {
            assert!(TextFile::get_line(&mut reader, &mut line).unwrap());
            assert_eq!(line, expected);
        }
        assert!(!TextFile::get_line(&mut reader, &mut line).unwrap());
        assert!(line.is_empty());
    }
    for (input, expected) in [
        ("", vec![]),
        ("\n", vec![""]),
        ("x\n", vec!["x"]),
        ("\n\n", vec!["", ""]),
        ("x", vec!["x"]),
    ] {
        assert_eq!(read(input, &ReadOptions::default()).lines(), expected);
    }
}
#[test]
fn trimming_comments_empty_lines_and_first_n_use_retained_lines() {
    let input = " # spaced\n#raw\n \t\nA\nB\nC";
    let options = ReadOptions {
        comment_symbol: "#".into(),
        first_n: 2,
        ..Default::default()
    };
    assert_eq!(read(input, &options).lines(), [" # spaced", " \t"]);
    let trim = ReadOptions {
        trim_lines: true,
        ..options.clone()
    };
    assert_eq!(read(input, &trim).lines(), ["", "A"]);
    let skip = ReadOptions {
        skip_empty_lines: true,
        ..trim
    };
    assert_eq!(read(input, &skip).lines(), ["A", "B"]);
    for first_n in [0, -1, -2, i32::MIN] {
        assert_eq!(
            read(
                input,
                &ReadOptions {
                    first_n,
                    ..skip.clone()
                }
            )
            .lines(),
            ["A", "B", "C"]
        );
    }
    // Source trims exactly space/tab/CR/LF, not every Unicode whitespace character.
    assert_eq!(
        read(
            "\u{b}v\u{c}\n\u{a0}x\u{a0}",
            &ReadOptions {
                trim_lines: true,
                ..Default::default()
            }
        )
        .lines(),
        ["\u{b}v\u{c}", "\u{a0}x\u{a0}"]
    );
    assert_eq!(
        read(
            "//comment\n/data\nvalue",
            &ReadOptions {
                comment_symbol: "//".into(),
                ..Default::default()
            }
        )
        .lines(),
        ["/data", "value"]
    );
}
#[test]
fn early_stop_does_not_consume_or_validate_the_remaining_input() {
    let mut reader = Cursor::new(b"a\r\n\xffbad\n");
    let file = TextFile::from_reader(
        &mut reader,
        &ReadOptions {
            first_n: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(file.lines(), ["a"]);
    assert_eq!(reader.position(), 3);
    let mut remaining = Vec::new();
    reader.read_to_end(&mut remaining).unwrap();
    assert_eq!(remaining, b"\xffbad\n");
}
#[test]
fn source_store_suffix_rules_and_mutable_reverse_iteration() {
    let mut file = TextFile::new();
    for line in [
        "line1",
        "line2\n",
        "line3\r\n",
        "bare\r",
        "two\n\n",
        "",
        "a\r\nb",
    ] {
        file.add_line(line).unwrap();
    }
    let mut out = Vec::new();
    file.write(&mut out).unwrap();
    let unix = "line1\nline2\nline3\nbare\r\ntwo\n\n\na\r\nb\n";
    let expected = if cfg!(windows) {
        unix.replace('\n', "\r\n")
    } else {
        unix.into()
    };
    assert_eq!(out, expected.as_bytes());
    let mut fixture = read("one\ntwo", &ReadOptions::default());
    assert_eq!(
        fixture.iter().rev().map(String::as_str).collect::<Vec<_>>(),
        ["two", "one"]
    );
    for line in &mut fixture {
        line.push('!');
    }
    assert_eq!(
        (&fixture)
            .into_iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        ["one!", "two!"]
    );
    fixture.clear();
    assert!(fixture.is_empty());
}
#[test]
fn limits_include_skipped_input_and_errors_preserve_destination() {
    let original = read("old", &ReadOptions::default());
    for limits in [
        Limits {
            max_input_bytes: 2,
            ..Default::default()
        },
        Limits {
            max_line_bytes: 1,
            ..Default::default()
        },
        Limits {
            max_lines: 1,
            ..Default::default()
        },
        Limits {
            max_storage_bytes: 1,
            ..Default::default()
        },
    ] {
        let mut destination = original.clone();
        let options = ReadOptions {
            limits,
            comment_symbol: "#".into(),
            ..Default::default()
        };
        assert!(
            destination
                .load_reader(Cursor::new(b"#comment\nvalue\n"), &options)
                .is_err()
        );
        assert_eq!(destination, original);
    }
    let mut destination = original.clone();
    assert!(
        destination
            .load_reader(Cursor::new(b"valid\n\xff"), &ReadOptions::default())
            .is_err()
    );
    assert_eq!(destination, original);
    let exact = Limits {
        max_input_bytes: 3,
        max_line_bytes: 1,
        max_lines: 1,
        ..Default::default()
    };
    TextFile::from_reader(
        Cursor::new(b"x\r\n"),
        &ReadOptions {
            limits: exact,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        TextFile::from_reader(
            Cursor::new(b"x\r\n"),
            &ReadOptions {
                limits: Limits {
                    max_input_bytes: 2,
                    ..exact
                },
                ..Default::default()
            }
        )
        .is_err()
    );
    let mut buffer = "unchanged".to_owned();
    assert!(
        TextFile::get_line_with_limits(
            &mut Cursor::new(b"ab"),
            &mut buffer,
            &Limits {
                max_line_bytes: 1,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert_eq!(buffer, "unchanged");
}
#[derive(Default)]
struct Sink {
    bytes: usize,
    flushes: usize,
}
impl Write for Sink {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes += bytes.len();
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.flushes += 1;
        Ok(())
    }
}
#[test]
fn appended_or_mutated_storage_and_output_are_preflighted() {
    let mut file = TextFile::with_limits(Limits {
        max_lines: 1,
        max_line_bytes: 10,
        max_output_bytes: 3,
        ..Default::default()
    })
    .unwrap();
    file.add_line("abc").unwrap();
    let old = file.clone();
    assert!(file.add_line("next").is_err());
    assert_eq!(file, old);
    let mut sink = Sink::default();
    assert!(file.write(&mut sink).is_err());
    assert_eq!((sink.bytes, sink.flushes), (0, 0));
    file.iter_mut().next().unwrap().push_str("0123456789");
    assert!(file.write(&mut sink).is_err());
    assert_eq!(sink.bytes, 0);
    let mut small = TextFile::with_limits(Limits {
        max_storage_bytes: 200,
        ..Default::default()
    })
    .unwrap();
    small.add_line("a").unwrap();
    small.iter_mut().next().unwrap().reserve(1000);
    assert!(small.add_line("b").is_err());
}
#[test]
fn path_operations_and_stream_failures_are_checked() {
    let path = std::env::temp_dir().join(format!("openms-text-helper-{}.txt", std::process::id()));
    let mut file = TextFile::new();
    for line in ["line1", "line2\n", "line3\r\n"] {
        file.add_line(line).unwrap();
    }
    file.store(&path).unwrap();
    assert_eq!(
        TextFile::from_path(&path, &ReadOptions::default())
            .unwrap()
            .lines(),
        ["line1", "line2", "line3"]
    );
    let mut loaded = TextFile::new();
    loaded.load(&path, &ReadOptions::default()).unwrap();
    std::fs::remove_file(&path).unwrap();
    let before = loaded.clone();
    assert!(loaded.load(&path, &ReadOptions::default()).is_err());
    assert_eq!(loaded, before);
    struct BadFlush;
    impl Write for BadFlush {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("flush"))
        }
    }
    assert!(loaded.write(BadFlush).is_err());
}
