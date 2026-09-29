// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `InspectOutfile`: the upstream class-test sections
//! (`InspectOutfile_test.cpp` at core `bc9cc12`, fixtures copied from its
//! `src/tests/class_tests/openms/data/`). Where the class test's own check
//! cannot run (its `getSequences` section tests nothing, because its guard
//! expects one sequence for two records), the expectation is derived
//! independently from the fixture, or taken from the executed C++ Release
//! build (`tests/data/search_engine_oracle/inspect_*.out`), and says so.
#![cfg(feature = "mzml")]

use openms::Error;
use openms::data_structures::DateTime;
use openms::format::inspect_outfile::{InspectOutfile, TRIE_DELIMITER};
use openms::identification::{FlankingResidue, PeptideIdentification, ProteinIdentification};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("f3b_inspect_{tag}_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A retained output of the executed C++ Release build (ibminode06,
/// `openms4-release-bc9cc12-c19e494-174b576`; driver
/// `../oracle/F3b-search-engine/driver.cpp`, outside the repository).
fn oracle(name: &str) -> String {
    std::fs::read_to_string(data(&format!("search_engine_oracle/{name}"))).unwrap()
}

/// The oracle driver's rendering of bytes: printable ASCII as is (with line
/// feed and tab when `keep_whitespace`, as its `dumpFile` does), anything else
/// as `\x` and lower-case hex without padding.
fn oracle_escape(bytes: &[u8], keep_whitespace: bool) -> String {
    let mut out = String::new();
    for &b in bytes {
        if (32..127).contains(&b) || (keep_whitespace && (b == b'\n' || b == b'\t')) {
            out.push(char::from(b));
        } else {
            out.push_str(&format!("\\x{b:x}"));
        }
    }
    out
}

fn parse_message(error: Error) -> String {
    match error {
        Error::Parse { message, .. } => message,
        other => panic!("expected a parse error, got {other:?}"),
    }
}

fn assert_not_found(error: Error, name: &str) {
    match error {
        Error::Io(e) => assert_eq!(
            e.to_string(),
            format!("the file '{name}' could not be found")
        ),
        other => panic!("expected a missing-file error, got {other:?}"),
    }
}

/// `TEST_REAL_SIMILAR`: relative tolerance 1e-5, the class tests' default.
fn assert_similar(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= 1e-5 * scale,
        "{actual} vs {expected}"
    );
}

const HEADER: &str = "#SpectrumFile\tScan#\tAnnotation\tProtein\tCharge\tMQScore\tCutScore\tIntenseBY\tBYPresent\tUnused\tp-value\tDeltaScore\tDeltaScoreOther\tRecordNumber\tDBFilePos\tSpecFilePos";

/// The class test's result lines, after the spectrum file.
const LINES: [&str; 19] = [
    "4\tN.EER.N\tQ9CQV8|1433B_MOUSE\t1\t-6889\t1150\t200\t500\t0\t0.95000\t0\t0\t1\t287\t270451",
    "4\tR.EKIE.K\tP68509|1433F_BOVIN\t1\t-1456\t3388\t300\t667\t0\t0.01000\t-1199\t-1199\t0\t87\t276369",
    "4\tE.KKLE.K\tP68509|1433F_BOVIN\t1\t-257\t3820\t450\t333\t0\t0.00001\t1199\t1199\t0\t78\t276369",
    "25\tD.DAIAE.L\tP68509|1433F_BOVIN\t1\t-6526\t354\t250\t375\t0\t0.97500\t-6269\t-6269\t0\t205\t276369",
    "37\tD.EAIAEL.D\tQ9CQV8|1433B_MOUSE\t1\t-2807\t2464\t350\t300\t0\t0.98536\t0\t0\t1\t446\t387302",
    "49\tD.KFLIK.N\tP68509|1433F_BOVIN\t1\t-5308\t1813\t250\t375\t0\t0.98758\t0\t0\t0\t106\t473207",
    "120\tN.EKKLEKVKA.Y\tP68509|1433F_BOVIN\t2\t-9115\t-1187\t185\t125\t0\t0.96296\t0\t0\t0\t77\t1062539",
    "217\tN.EDRNLL.S\tP68509|1433F_BOVIN\t1\t-11835\t140\t100\t0\t0\t0.97143\t0\t0\t0\t41\t1720635",
    "249\tA.LLDKF.L\tP68509|1433F_BOVIN\t1\t-4503\t2554\t300\t250\t0\t0.99043\t0\t0\t0\t103\t1943362",
    "501\tF.DEAIAELDTLNEE.S\tQ9CQV8|1433B_MOUSE\t2\t-7874\t-1829\t231\t0\t0\t0.95652\t0\t0\t1\t452\t3860094",
    "667\tK.LAEQAERYDDMAA.A\tQ9CQV8|1433B_MOUSE\t2\t-9787\t-1362\t77\t125\t0\t0.96552\t0\t0\t1\t266\t5279013",
    "685\tN.LTLWTSENQGDEGDAG.E\tQ9CQV8|1433B_MOUSE\t2\t-9056\t-2174\t125\t0\t0\t0.96296\t0\t0\t1\t479\t5448607",
    "736\tY.QEAFEIS.K\tQ9CQV8|1433B_MOUSE\t1\t-11507\t478\t0\t250\t0\t0.98789\t17\t17\t1\t399\t6018155",
    "736\tY.KEAFEIS.K\tP68509|1433F_BOVIN\t1\t-11524\t457\t0\t250\t0\t0.97059\t-17\t0\t0\t156\t6018155",
    "758\tS.NEDRNLLSVAYKN.V\tP68509|1433F_BOVIN\t2\t-4841\t-481\t256\t208\t0\t0.98948\t0\t0\t0\t43\t6167475",
    "764\tL.LAKQAFDDAIAELDTLNED.S\tP68509|1433F_BOVIN\t2\t-4474\t-560\t140\t250\t0\t0.99043\t0\t0\t0\t201\t6208454",
    "786\tT.MDKSELV.Q\tQ9CQV8|1433B_MOUSE\t1\t-12849\t-1629\t48\t167\t0\t0.97368\t0\t0\t1\t253\t6399323",
    "1045\tS.VFYYEI.Q\tP68509|1433F_BOVIN\t1\t-15475\t-2579\t0\t100\t0\t0.97500\t0\t0\t0\t184\t9842931",
    "1962\tE.AFEIS.K\tP68509|1433F_BOVIN\t1\t-10496\t-742\t100\t375\t0\t0.96774\t0\t0\t0\t159\t19098337",
];

/// The class test's two result files; they name the spectrum file by its
/// absolute path, as the class test's do.
fn result_files(dir: &Path) -> (PathBuf, PathBuf) {
    let spectrum_file1 = data("InspectOutfile_test_1.mzXML");
    let spectrum_file1 = spectrum_file1.to_str().unwrap();
    let mut first = format!("{HEADER}\n");
    for line in LINES {
        first.push_str(&format!("{spectrum_file1}\t{line}\n"));
    }
    let second = format!(
        "{HEADER}\n{spectrum_file1}\tN.EER.N\tQ9CQV8|1433B_MOUSE\t1\t-6889\t1150\t200\t500\t0\t0.95000\t0\t0\t1\t287\t270451\n"
    );
    let (a, b) = (dir.join("result1.txt"), dir.join("result2.txt"));
    std::fs::write(&a, first).unwrap();
    std::fs::write(&b, second).unwrap();
    (a, b)
}

fn run_date() -> DateTime {
    "2026-09-28 12:00:00".parse().unwrap()
}

#[test]
fn equality_is_stateless() {
    // Through a generic, so the `Default` impl itself is what is exercised.
    fn by_default<T: Default>() -> T {
        T::default()
    }
    assert_eq!(InspectOutfile::new(), by_default::<InspectOutfile>());
}

#[test]
fn load_rejects_bad_arguments_and_files() {
    let file = InspectOutfile::new();
    let mut peptides = Vec::new();
    let mut protein = ProteinIdentification::default();
    for p in [2.0, -1.0] {
        assert!(matches!(
            file.load("", &mut peptides, &mut protein, p, None),
            Err(Error::InvalidValue(_))
        ));
    }
    assert_not_found(
        file.load("a", &mut peptides, &mut protein, 0.01, None)
            .unwrap_err(),
        "a",
    );
    let empty = data("Inspect_empty_file.txt");
    assert_eq!(
        parse_message(
            file.load(&empty, &mut peptides, &mut protein, 0.01, None)
                .unwrap_err()
        ),
        format!("the file '{}' is empty", empty.display())
    );
}

#[test]
fn load() {
    let dir = temp_dir("load");
    let (first, second) = result_files(&dir);
    let file = InspectOutfile::new();
    let mut peptides: Vec<PeptideIdentification> = Vec::new();
    let mut protein = ProteinIdentification::default();

    file.load_at(&first, &mut peptides, &mut protein, 0.001, None, run_date())
        .unwrap();
    assert_eq!(peptides.len(), 1);
    let id = &peptides[0];
    assert_eq!(id.hits.len(), 1);
    assert_eq!(id.score_type, "Inspect");
    assert_similar(id.significance_threshold, 0.001);
    assert_eq!(id.identifier, "InsPecT_2026-09-28");
    let hit = &id.hits[0];
    assert_similar(hit.score, -257.0);
    assert_eq!(hit.sequence.to_string(), "KKLE");
    assert_eq!(hit.evidences.len(), 1);
    assert_eq!(
        hit.evidences[0].aa_before,
        FlankingResidue::from_code('E').unwrap()
    );
    assert_eq!(
        hit.evidences[0].aa_after,
        FlankingResidue::from_code('K').unwrap()
    );
    assert_eq!(hit.charge, 1);
    assert_eq!(hit.evidences[0].protein_accession, "P68509");
    assert_eq!(protein.identifier, "InsPecT_2026-09-28");
    assert_eq!(protein.hits.len(), 1);
    assert_eq!(protein.hits[0].accession, "P68509");

    peptides.clear();
    file.load_at(&first, &mut peptides, &mut protein, 0.01, None, run_date())
        .unwrap();
    assert_eq!(peptides.len(), 1);
    let id = &peptides[0];
    assert_eq!(id.hits.len(), 2);
    assert_eq!(id.score_type, "Inspect");
    assert_similar(id.significance_threshold, 0.01);
    let expected = [(-257.0, "KKLE", 'E', 'K'), (-1456.0, "EKIE", 'R', 'K')];
    for (hit, (score, sequence, before, after)) in id.hits.iter().zip(expected) {
        assert_similar(hit.score, score);
        assert_eq!(hit.sequence.to_string(), sequence);
        assert_eq!(hit.evidences.len(), 1);
        assert_eq!(
            hit.evidences[0].aa_before,
            FlankingResidue::from_code(before).unwrap()
        );
        assert_eq!(
            hit.evidences[0].aa_after,
            FlankingResidue::from_code(after).unwrap()
        );
        assert_eq!(hit.charge, 1);
        assert_eq!(hit.evidences[0].protein_accession, "P68509");
    }
    // Scan 4 of the first spectrum file (see precursor_rt_and_mz).
    assert_similar(id.rt.unwrap(), -1.0);
    assert_similar(id.mz.unwrap(), 123.456);

    peptides.clear();
    protein.hits.clear();
    let corrupted = file
        .load_at(&second, &mut peptides, &mut protein, 0.01, None, run_date())
        .unwrap();
    assert_eq!(corrupted, [2]);
    assert!(peptides.is_empty());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn load_reads_protein_sequences_from_the_trie_database() {
    // Independently derived: records 0 and 1 of Inspect_test.trie are the
    // two '*'-separated sequences of the file (P68509 and Q9CQV8, the
    // accessions of RecordNumber 0 and 1 in the result lines).
    let dir = temp_dir("load_db");
    // The class test's scan-4 lines, all of which pass a threshold of 1.
    let spectrum_file1 = data("InspectOutfile_test_1.mzXML");
    let mut text = format!("{HEADER}\n");
    for line in &LINES[..3] {
        text.push_str(&format!("{}\t{line}\n", spectrum_file1.display()));
    }
    let first = dir.join("result.txt");
    std::fs::write(&first, text).unwrap();
    let trie = std::fs::read(data("Inspect_test.trie")).unwrap();
    let records: Vec<String> = trie
        .split(|&b| b == TRIE_DELIMITER)
        .map(|r| String::from_utf8(r.to_vec()).unwrap())
        .collect();
    let mut peptides = Vec::new();
    let mut protein = ProteinIdentification::default();
    InspectOutfile::new()
        .load_at(
            &first,
            &mut peptides,
            &mut protein,
            1.0,
            Some(&data("Inspect_test.trie")),
            run_date(),
        )
        .unwrap();
    // The first passing line names Q9CQV8 (record 1), the second P68509.
    assert_eq!(protein.hits.len(), 2);
    assert_eq!(protein.hits[0].accession, "Q9CQV8");
    assert_eq!(protein.hits[0].sequence, records[1]);
    assert_eq!(protein.hits[1].accession, "P68509");
    assert_eq!(protein.hits[1].sequence, records[0]);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn generate_trie_db() {
    let file = InspectOutfile::new();
    assert_not_found(
        file.generate_trie_db("a", "", "", false, "").unwrap_err(),
        "a",
    );
    let dir = temp_dir("generate");
    assert!(matches!(
        file.generate_trie_db(
            data("Inspect_test.fasta"),
            dir.join("missing/x.trie"),
            dir.join("x.index"),
            false,
            ""
        ),
        Err(Error::Io(_))
    ));
    let mut rendered = String::new();
    for (fasta, trie, index) in [
        (
            "Inspect_test.fasta",
            "Inspect_test.trie",
            "Inspect_test.index",
        ),
        (
            "Inspect_test2.fasta",
            "Inspect_test2.trie",
            "Inspect_test2.index",
        ),
    ] {
        let (out_trie, out_index) = (dir.join("out.trie"), dir.join("out.index"));
        file.generate_trie_db(data(fasta), &out_trie, &out_index, false, "")
            .unwrap();
        let written_trie = std::fs::read(&out_trie).unwrap();
        let written_index = std::fs::read(&out_index).unwrap();
        // The class test's TEST_FILE_EQUAL compares trimmed lines, so the
        // fixture's final line feed, which `generateTrieDB` never writes,
        // does not count; the index matches byte for byte.
        let mut expected_trie = std::fs::read(data(trie)).unwrap();
        assert_eq!(expected_trie.pop(), Some(b'\n'));
        assert_eq!(written_trie, expected_trie);
        assert_eq!(written_index, std::fs::read(data(index)).unwrap());
        let stem = trie.trim_end_matches(".trie");
        for (part, bytes) in [("trie", &written_trie), ("index", &written_index)] {
            rendered.push_str(&format!(
                "case {stem} {part}\nfile_bytes {}\n{}\n<<EOF\n",
                bytes.len(),
                oracle_escape(bytes, true)
            ));
        }
    }
    // The executed Release build writes exactly these bytes.
    assert_eq!(rendered, oracle("inspect_trie.out"));
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn compress_trie_db() {
    let file = InspectOutfile::new();
    let dir = temp_dir("compress");
    let (trie, index) = (data("Inspect_test.trie"), data("Inspect_test.index"));
    let (out_trie, out_index) = (dir.join("out.trie"), dir.join("out.index"));
    let wanted = [0usize];
    assert_eq!(
        parse_message(
            file.compress_trie_db(&trie, &index, &wanted, &trie, "", false)
                .unwrap_err()
        ),
        format!(
            "{} in: Same filename can not be used for original and second database!",
            trie.display()
        )
    );
    assert_eq!(
        parse_message(
            file.compress_trie_db(&trie, &index, &wanted, "", &index, false)
                .unwrap_err()
        ),
        format!(
            "{} in: Same filename can not be used for original and second database!",
            index.display()
        )
    );
    assert_not_found(
        file.compress_trie_db("a", "", &wanted, &out_trie, &out_index, false)
            .unwrap_err(),
        "a",
    );
    assert_not_found(
        file.compress_trie_db(&trie, "b", &wanted, &out_trie, &out_index, false)
            .unwrap_err(),
        "b",
    );
    assert!(matches!(
        file.compress_trie_db(&trie, &index, &wanted, dir.join("missing/x"), "", true),
        Err(Error::Io(_))
    ));
    let empty = data("Inspect_empty_file.txt");
    assert_eq!(
        parse_message(
            file.compress_trie_db(&trie, &empty, &wanted, &out_trie, &out_index, true)
                .unwrap_err()
        ),
        format!("{} in: index file is too short!", empty.display())
    );
    std::fs::remove_file(&out_trie).unwrap();
    std::fs::remove_file(&out_index).unwrap();

    // The class test: append record 0 of the first database, all of the
    // second and record 1 of the first, then compress records 0 and 3 of
    // that into a new database, which is the first database again.
    file.compress_trie_db(&trie, &index, &[0], &out_trie, &out_index, true)
        .unwrap();
    file.compress_trie_db(
        data("Inspect_test2.trie"),
        data("Inspect_test2.index"),
        &[],
        &out_trie,
        &out_index,
        true,
    )
    .unwrap();
    file.compress_trie_db(&trie, &index, &[1], &out_trie, &out_index, true)
        .unwrap();
    let (final_trie, final_index) = (dir.join("final.trie"), dir.join("final.index"));
    file.compress_trie_db(
        &out_trie,
        &out_index,
        &[0, 3],
        &final_trie,
        &final_index,
        false,
    )
    .unwrap();
    assert_eq!(
        std::fs::read(&final_trie).unwrap(),
        std::fs::read(&trie).unwrap()
    );
    assert_eq!(
        std::fs::read(&final_index).unwrap(),
        std::fs::read(&index).unwrap()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn sequences() {
    let file = InspectOutfile::new();
    let wanted: BTreeMap<usize, usize> = [(0, 0), (1, 1)].into_iter().collect();
    let mut found = Vec::new();
    assert_not_found(file.sequences("a", &wanted, &mut found).unwrap_err(), "a");

    // Independently derived: the records are the '*'-separated parts.
    for name in ["Inspect_test.trie", "Inspect_test2.trie"] {
        let bytes = std::fs::read(data(name)).unwrap();
        let records: Vec<String> = bytes
            .split(|&b| b == TRIE_DELIMITER)
            .map(|r| String::from_utf8(r.to_vec()).unwrap())
            .collect();
        assert_eq!(records.len(), 2);
        let mut found = Vec::new();
        let not_found = file.sequences(data(name), &wanted, &mut found).unwrap();
        assert!(not_found.is_empty());
        assert_eq!(found, records);
    }
    // The class test's two literals: the first sequence of Inspect_test.trie
    // and the second of Inspect_test2.trie. The class test never compares the
    // second (its guard is false), and the pinned Release build returns it
    // with the fixture's final line feed, because `getSequences` reads the
    // last record up to the end of the file
    // (search_engine_oracle/inspect_sequences.out).
    let mut found = Vec::new();
    file.sequences(data("Inspect_test.trie"), &wanted, &mut found)
        .unwrap();
    assert_eq!(
        found[0],
        "GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN"
    );
    let mut found = Vec::new();
    file.sequences(data("Inspect_test2.trie"), &wanted, &mut found)
        .unwrap();
    assert_eq!(
        found[1],
        "SAPPSLLVLYFGKKELRAMKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPP\n"
    );

    // The executed Release build, record by record.
    let mut rendered = String::new();
    for name in ["Inspect_test.trie", "Inspect_test2.trie"] {
        let mut found = Vec::new();
        let not_found = file.sequences(data(name), &wanted, &mut found).unwrap();
        rendered.push_str(&format!(
            "case {name} not_found {} found {}\n",
            not_found.len(),
            found.len()
        ));
        for sequence in &found {
            rendered.push_str(&format!(
                "seq {}\n{}\n",
                sequence.len(),
                oracle_escape(sequence.as_bytes(), false)
            ));
        }
    }
    assert_eq!(rendered, oracle("inspect_sequences.out"));

    // A record beyond the database is not found, and so is every later one.
    let beyond: BTreeMap<usize, usize> = [(0, 0), (5, 1), (9, 2)].into_iter().collect();
    let mut found = Vec::new();
    let not_found = file
        .sequences(data("Inspect_test.trie"), &beyond, &mut found)
        .unwrap();
    assert_eq!(not_found, [5, 9]);
    assert_eq!(found.len(), 3);
    assert!(found[1].is_empty() && found[2].is_empty());
}

#[test]
fn ac_and_ac_type() {
    let file = InspectOutfile::new();
    for (line, accession, kind) in [
        (
            ">sp|P02666|CASB_BOVIN Beta-casein precursor - Bos taurus (Bovine).",
            "P02666",
            "SwissProt",
        ),
        (
            ">tr|Q5EEQ7|Q5EEQ7_BOVIN Beta-casein (Fragment) - Bos taurus (Bovine).",
            "Q5EEQ7",
            "SwissProt",
        ),
        ("gi|110174602|gb|DQ660451.1|", "DQ660451.1", "GenBank"),
        ("gi|1655698|emb|Y07752|VCPHEROPH", "Y07752", "EMBL"),
        ("gi|10038695|dbj|BAB12730|", "BAB12730", "DDBJ"),
        ("gi|9628804|ref|NP_043835|", "NP_043835", "NCBI"),
        ("gi|21362794|sp|P58858.0|", "P58858.0", "SwissProt"),
        ("gi|21362794|tr|P58858.0|", "P58858.0", "SwissProt"),
        ("gi|1619818|gnl|PID|d1013471|", "d1013471", "PID"),
        (
            "Q30DX2 Gamma-gliadin/LMW-glutenin chimera Ch7 (Fragment).",
            "Q30DX2",
            "SwissProt",
        ),
        (">P68509|1433F_BOVIN", "P68509", "SwissProt"),
        (">ACBLA (P68509) F_BOVIN", "P68509", "SwissProt"),
    ] {
        assert_eq!(
            file.ac_and_ac_type(line),
            (accession.to_owned(), kind.to_owned()),
            "{line}"
        );
    }
}

#[test]
fn precursor_rt_and_mz() {
    let file = InspectOutfile::new();
    let spectrum_file1 = data("InspectOutfile_test_1.mzXML")
        .to_str()
        .unwrap()
        .to_owned();
    let spectrum_file2 = data("InspectOutfile_test_2.mzXML")
        .to_str()
        .unwrap()
        .to_owned();
    let mut ids = vec![PeptideIdentification::default(); 2];
    let too_far = vec![(spectrum_file1.clone(), vec![(0, 10)])];
    assert_eq!(
        parse_message(file.precursor_rt_and_mz(&too_far, &mut ids).unwrap_err()),
        format!(
            "{spectrum_file1} in: Not enough scans in file! (4 available, should be at least 10)"
        )
    );
    let files = vec![
        (spectrum_file1, vec![(0, 4)]),
        (spectrum_file2, vec![(1, 4)]),
    ];
    file.precursor_rt_and_mz(&files, &mut ids).unwrap();
    assert_similar(ids[0].rt.unwrap(), -1.0);
    assert_similar(ids[0].mz.unwrap(), 123.456);
    assert_similar(ids[1].rt.unwrap(), 180.0);
    assert_similar(ids[1].mz.unwrap(), 123.456);
}

#[test]
fn labels() {
    let file = InspectOutfile::new();
    assert_not_found(file.labels("a").unwrap_err(), "a");
    let unknown = data("Inspect_test1.fasta");
    assert_eq!(
        parse_message(file.labels(&unknown).unwrap_err()),
        format!(
            "{} in: database has unknown file format (neither trie nor FASTA nor swissprot)",
            unknown.display()
        )
    );
    let labels = file.labels(data("Inspect_test.fasta")).unwrap();
    assert_eq!(labels.ac_label, ">");
    assert_eq!(labels.sequence_start_label, ">");
    assert_eq!(labels.sequence_end_label, ">");
    assert_eq!(labels.comment_label, ";");
    assert_eq!(labels.species_label, ">");
}

#[test]
fn wanted_records() {
    let file = InspectOutfile::new();
    for p in [2.0, -1.0] {
        assert!(matches!(
            file.wanted_records("", p),
            Err(Error::InvalidValue(_))
        ));
    }
    assert_not_found(file.wanted_records("a", 0.01).unwrap_err(), "a");
    let empty = data("Inspect_empty_file.txt");
    assert_eq!(
        parse_message(file.wanted_records(&empty, 0.01).unwrap_err()),
        format!("the file '{}' is empty", empty.display())
    );
    let dir = temp_dir("wanted");
    let (first, _) = result_files(&dir);
    assert_eq!(file.wanted_records(&first, 0.01).unwrap(), [0]);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn experiment() {
    let file = InspectOutfile::new();
    // The executed Release build (search_engine_oracle/inspect_experiment.out):
    // both class-test spectrum files load with 20 spectra, and a FASTA file,
    // recognised by content, is refused by `FileHandler::loadExperiment`.
    let (experiment, kind) = file.experiment(data("Inspect.mzXML")).unwrap();
    assert_eq!((kind.as_str(), experiment.spectra.len()), ("mzXML", 20));
    let (experiment, kind) = file.experiment(data("Inspect.mzData")).unwrap();
    assert_eq!((kind.as_str(), experiment.spectra.len()), ("mzData", 20));
    let fasta = data("Inspect_test.fasta");
    assert_eq!(
        parse_message(file.experiment(&fasta).unwrap_err()),
        format!(
            "type is not supported for loading experiments in: {}",
            fasta.display()
        )
    );
    let expected = oracle("inspect_experiment.out");
    assert!(expected.contains("case Inspect.mzXML by_content mzXML\ntype mzXML spectra 20\n"));
    assert!(expected.contains("case Inspect.mzData by_content mzData\ntype mzData spectra 20\n"));
    assert!(expected.contains(
        "case Inspect_test.fasta by_content fasta\nexception Parse Error message type is not supported for loading experiments in: <DATA>/Inspect_test.fasta\n"
    ));

    // Content that no type matches (derived: the source's own branch).
    let dir = temp_dir("experiment");
    let unknown = dir.join("unknown.txt");
    std::fs::write(&unknown, "no spectra here\n").unwrap();
    assert_eq!(
        parse_message(file.experiment(&unknown).unwrap_err()),
        format!(
            "{} in: Could not determine type of the file. Aborting!",
            unknown.display()
        )
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn search_engine_and_version() {
    let file = InspectOutfile::new();
    let mut protein = ProteinIdentification::default();
    // The class test's C++ literals join their lines without line feeds.
    let output = "InsPecT vesrion 20060907  Interpretation of Peptides with Post-translational Modifications.  Copyright 2006, The Regents of the University of California  [See Docs directory for usage manual and copyright information]Sample command-line:Inspect.exe -i Foo.in -o Foo.txt -e ErrorsFoo.txt";
    assert!(file.search_engine_and_version(output, &mut protein));
    assert_eq!(protein.search_engine, "InsPecT");
    assert_eq!(protein.search_engine_version, "20060907");
    let output = "InsPecT version 20100331  Interpretation of Peptides with Post-translational Modifications.  Copyright 2007,2008,2009 The Regents of the University of California";
    assert!(file.search_engine_and_version(output, &mut protein));
    assert_eq!(protein.search_engine, "InsPecT");
    assert_eq!(protein.search_engine_version, "20100331");
    assert!(!file.search_engine_and_version("no version here", &mut protein));
    assert_eq!(protein.search_engine_version, "unknown");
}

#[test]
fn read_out_header() {
    let file = InspectOutfile::new();
    let header_line = "#SpectrumFile\tScan#\tAnnotation\tProtein\tCharge\tMQScore\tLength\tTotalPRMScore\tMedianPRMScore\tFractionY\tFractionB\tIntensity\tNTT\tp-value\tF-Score\tDeltaScore\tDeltaScoreOther\tRecordNumber\tDBFilePos\t";
    assert_eq!(
        parse_message(
            file.read_out_header("dummy_testfile", header_line)
                .unwrap_err()
        ),
        "dummy_testfile in: at least one of the columns '#SpectrumFile', 'Scan#', 'Annotation', 'Protein', 'Charge', 'MQScore', 'p-value', 'RecordNumber', 'DBFilePos' or 'SpecFilePos' is missing!"
    );
    let header_line = "#SpectrumFile\tScan#\tAnnotation\tProtein\tCharge\tMQScore\tLength\tTotalPRMScore\tMedianPRMScore\tFractionY\tFractionB\tIntensity\tNTT\tp-value\tF-Score\tDeltaScore\tDeltaScoreOther\tRecordNumber\tDBFilePos\tSpecFilePos";
    let header = file.read_out_header("dummy_testfile", header_line).unwrap();
    assert_eq!(header.spectrum_file_column, 0);
    assert_eq!(header.scan_column, 1);
    assert_eq!(header.peptide_column, 2);
    assert_eq!(header.protein_column, 3);
    assert_eq!(header.charge_column, 4);
    assert_eq!(header.mq_score_column, 5);
    assert_eq!(header.p_value_column, 13);
    assert_eq!(header.record_number_column, 17);
    assert_eq!(header.db_file_pos_column, 18);
    assert_eq!(header.spec_file_pos_column, 19);
    assert_eq!(header.number_of_columns, 20);
}
