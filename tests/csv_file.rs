// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use openms::format::{
    CsvFile,
    csv::{Limits, ReadOptions},
};
use std::io::Cursor;
fn from(text: &str, options: &ReadOptions) -> CsvFile {
    CsvFile::from_reader(Cursor::new(text), options).unwrap()
}
#[test]
fn all_source_csv_literal_rows_constructor_load_and_store() {
    let options = ReadOptions {
        separator: b'\t',
        ..Default::default()
    };
    let plain = include_str!("data/csv_file_source_1.csv");
    let quoted = include_str!("data/csv_file_source_2.csv");
    let constructor = from(plain, &options);
    let mut loaded = CsvFile::new();
    loaded.load_reader(Cursor::new(plain), &options).unwrap();
    let mut enclosed = CsvFile::new();
    enclosed
        .load_reader(
            Cursor::new(quoted),
            &ReadOptions {
                item_enclosed: true,
                ..options
            },
        )
        .unwrap();
    for file in [&constructor, &loaded, &enclosed] {
        assert_eq!(file.row_count(), 3);
        for (i, expected) in [["hello", "world"], ["the", "dude"], ["spectral", "search"]]
            .iter()
            .enumerate()
        {
            let (split, row) = file.row(i).unwrap();
            assert!(split);
            assert_eq!(&row, expected);
        }
    }
    let mut out = Vec::new();
    enclosed.write(&mut out).unwrap();
    let expected = if cfg!(windows) {
        quoted.replace('\n', "\r\n")
    } else {
        quoted.into()
    };
    assert_eq!(out, expected.as_bytes());
}
#[test]
fn constructor_and_load_preserve_distinct_comment_and_whitespace_rules() {
    let input = "  #comment,field\n#raw,skip\n  a,b  \n\n c,d\n";
    let options = ReadOptions::default();
    let raw = from(input, &options);
    assert_eq!(raw.row_count(), 4);
    assert_eq!(raw.row(0).unwrap().1, ["  #comment", "field"]);
    assert_eq!(raw.row(1).unwrap().1, ["  a", "b  "]);
    let mut loaded = CsvFile::new();
    loaded.load_reader(Cursor::new(input), &options).unwrap();
    assert_eq!(loaded.row_count(), 3);
    assert_eq!(loaded.row(0).unwrap().1, ["a", "b"]);
    assert_eq!(loaded.row(1).unwrap(), (false, vec![]));
    assert_eq!(loaded.row(2).unwrap().1, ["c", "d"]);
    for first_n in [0, -1, -9] {
        assert_eq!(
            from(input, &ReadOptions { first_n, ..options }).row_count(),
            4
        );
    }
    assert_eq!(
        from(
            input,
            &ReadOptions {
                first_n: 1,
                ..options
            }
        )
        .row(0)
        .unwrap()
        .1,
        ["  #comment", "field"]
    );
    loaded
        .load_reader(
            Cursor::new(input),
            &ReadOptions {
                first_n: 1,
                ..options
            },
        )
        .unwrap();
    assert_eq!(loaded.row_count(), 1);
    assert_eq!(loaded.row(0).unwrap().1, ["a", "b"]);
}
#[test]
fn split_boolean_empty_fields_and_enclosure_are_literal_source_rules() {
    let options = ReadOptions {
        item_enclosed: true,
        ..Default::default()
    };
    let file = from(
        "\n\"single\"\n\"a\",\"b\"\nabc,xyz\na,b\n\"\",\"\"\n,\n",
        &options,
    );
    for (i, split, expected) in [
        (0, false, vec![]),
        (1, false, vec!["\"single\""]),
        (2, true, vec!["a", "b"]),
        (3, true, vec!["b", "y"]),
        (4, true, vec!["", ""]),
        (5, true, vec!["", ""]),
    ] {
        assert_eq!(
            file.row(i).unwrap(),
            (split, expected.into_iter().map(str::to_owned).collect())
        );
    }
    let mut output = vec!["keep".to_owned()];
    assert!(file.get_row(6, &mut output).is_err());
    assert_eq!(output, ["keep"]);
    assert!(file.get_row(usize::MAX, &mut output).is_err());
    assert_eq!(output, ["keep"]);
    assert!(!file.get_row(1, &mut output).unwrap());
    assert_eq!(output, ["\"single\""]);
    let plain = from(",a,,\n", &ReadOptions::default());
    assert_eq!(
        plain.row(0).unwrap(),
        (true, vec!["".into(), "a".into(), "".into(), "".into()])
    );
    let rfc = from("\"a,b\",c", &ReadOptions::default());
    assert_eq!(rfc.row(0).unwrap().1, ["\"a", "b\"", "c"]);
}
#[test]
fn add_row_quotes_without_escaping_and_clear_retains_configuration() {
    let mut plain = CsvFile::new();
    plain.add_row(&["first", "second", "third"]).unwrap();
    plain.add_row(&["4", "5", "6"]).unwrap();
    assert_eq!(plain.row(0).unwrap().1, ["first", "second", "third"]);
    assert_eq!(plain.row(1).unwrap().1, ["4", "5", "6"]);
    let options = ReadOptions {
        separator: b';',
        item_enclosed: true,
        ..Default::default()
    };
    let mut quoted = CsvFile::with_options(&options).unwrap();
    quoted.add_row(&["a\"b", "c;d"]).unwrap();
    let mut bytes = Vec::new();
    quoted.write(&mut bytes).unwrap();
    let expected = if cfg!(windows) {
        "\"a\"b\";\"c;d\"\r\n"
    } else {
        "\"a\"b\";\"c;d\"\n"
    };
    assert_eq!(bytes, expected.as_bytes());
    quoted.clear();
    assert_eq!(quoted.row_count(), 0);
    assert_eq!(quoted.separator(), b';');
    assert!(quoted.item_enclosed());
    assert!(quoted.row(0).is_err());
    quoted.add_row::<&str>(&[]).unwrap();
    assert_eq!(quoted.row(0).unwrap(), (false, vec![]));
}
#[test]
fn ascii_byte_delimiters_and_utf8_boundary_checks() {
    let nul = ReadOptions {
        separator: 0,
        ..Default::default()
    };
    assert_eq!(from("a\0b", &nul).row(0).unwrap().1, ["a", "b"]);
    let mut unicode = CsvFile::with_options(&ReadOptions {
        item_enclosed: true,
        ..Default::default()
    })
    .unwrap();
    unicode.add_row(&["α", "β"]).unwrap();
    assert_eq!(unicode.row(0).unwrap().1, ["α", "β"]);
    // The source removes bytes, not Unicode characters. Two-byte tokens
    // become empty; three-byte tokens leave an invalid continuation byte.
    assert_eq!(
        from(
            "α,β",
            &ReadOptions {
                item_enclosed: true,
                ..Default::default()
            }
        )
        .row(0)
        .unwrap()
        .1,
        ["", ""]
    );
    assert!(
        from(
            "€,$",
            &ReadOptions {
                item_enclosed: true,
                ..Default::default()
            }
        )
        .row(0)
        .is_err()
    );
    let mut byte = from(
        "α",
        &ReadOptions {
            separator: 0xb1,
            ..Default::default()
        },
    );
    assert!(byte.row(0).is_err());
    let before = byte.clone();
    assert!(byte.add_row(&["a", "b"]).is_err());
    assert_eq!(byte, before);
    // A delimiter in an item is not escaped; a newline can therefore create
    // multiple rows on reload, exactly as in the source serializer.
    let mut newline = CsvFile::with_options(&ReadOptions {
        separator: b'\n',
        ..Default::default()
    })
    .unwrap();
    newline.add_row(&["a", "b"]).unwrap();
    assert_eq!(newline.row(0).unwrap().1, ["a", "b"]);
}
#[test]
fn load_and_row_failures_are_atomic_and_limits_cover_existing_appends() {
    let mut file = CsvFile::with_options(&ReadOptions {
        limits: Limits {
            max_lines: 1,
            ..Default::default()
        },
        ..Default::default()
    })
    .unwrap();
    file.add_row(&["a", "b"]).unwrap();
    let before = file.clone();
    assert!(file.add_row(&["c", "d"]).is_err());
    assert_eq!(file, before);
    assert!(
        file.load_reader(
            Cursor::new(b"first\n\xff"),
            &ReadOptions {
                separator: b';',
                item_enclosed: true,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert_eq!(file, before);
    let limits = Limits {
        max_line_bytes: 1,
        ..Default::default()
    };
    assert!(
        CsvFile::from_reader(
            Cursor::new("#too long"),
            &ReadOptions {
                limits,
                ..Default::default()
            }
        )
        .is_err()
    );
    let long = ",".repeat(openms::format::csv::MAX_FIELDS);
    let many = from(&long, &ReadOptions::default());
    assert!(many.row(0).is_err());
    let mut destination = vec!["unchanged".into()];
    assert!(many.get_row(0, &mut destination).is_err());
    assert_eq!(destination, ["unchanged"]);
}
#[test]
fn filesystem_and_writer_preflight_cover_source_store_surface() {
    let path = std::env::temp_dir().join(format!("openms-csv-helper-{}.csv", std::process::id()));
    let mut file = CsvFile::new();
    file.add_row(&["hello", "world"]).unwrap();
    file.store(&path).unwrap();
    let mut loaded = CsvFile::new();
    loaded.load(&path, &ReadOptions::default()).unwrap();
    assert_eq!(loaded.row(0).unwrap().1, ["hello", "world"]);
    assert_eq!(
        CsvFile::from_path(&path, &ReadOptions::default())
            .unwrap()
            .row_count(),
        1
    );
    let limited = from(
        "abc",
        &ReadOptions {
            limits: Limits {
                max_output_bytes: 1,
                ..Default::default()
            },
            ..Default::default()
        },
    );
    assert!(limited.store(&path).is_err());
    assert_eq!(
        CsvFile::from_path(&path, &ReadOptions::default())
            .unwrap()
            .row(0)
            .unwrap()
            .1,
        ["hello", "world"]
    );
    std::fs::remove_file(&path).unwrap();
    assert!(file.load(&path, &ReadOptions::default()).is_err());
}
