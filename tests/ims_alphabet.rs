use openms::chemistry::{
    IMSAlphabet, IMSAlphabetParser, IMSAlphabetTextParser, IMSElement, IMSIsotopeDistribution,
    IMSIsotopeOptions, IMSIsotopePeak,
};
use openms::system::file::TempDir;
use std::{
    collections::BTreeMap,
    io::{self, BufRead, BufReader, Cursor, Read},
};

fn basic() -> IMSAlphabet {
    IMSAlphabet::from_elements(vec![
        IMSElement::from_mass("hydrogen", 1.0).unwrap(),
        IMSElement::from_mass("oxygen", 16.0).unwrap(),
        IMSElement::from_mass("nitrogen", 14.0).unwrap(),
    ])
    .unwrap()
}

#[test]
fn source_container_lookup_and_mass_assertions() {
    assert!(IMSAlphabet::new().is_empty());
    let a = basic();
    assert_eq!(a.len(), 3);
    assert_eq!(a.name(0).unwrap(), "hydrogen");
    assert_eq!(a.name(1).unwrap(), "oxygen");
    assert_eq!(a.name(2).unwrap(), "nitrogen");
    for (name, index, mass) in [
        ("hydrogen", 0, 1.0),
        ("oxygen", 1, 16.0),
        ("nitrogen", 2, 14.0),
    ] {
        assert_eq!(a.get(name).unwrap(), a.element(index).unwrap());
        assert_eq!(a.mass_by_name(name).unwrap(), mass);
        assert_eq!(a.mass(index).unwrap(), mass);
        assert!(a.has_name(name).unwrap());
    }
    assert_eq!(a.masses(0).unwrap(), vec![1.0, 16.0, 14.0]);
    assert_eq!(a.average_masses().unwrap(), vec![1.0, 16.0, 14.0]);
    assert!(!a.has_name("oxygen2").unwrap());
    assert!(a.get("nitrogen2").is_err());
    assert!(a.mass_by_name("nitrogen2").is_err());
    assert!(a.element(3).is_err());
    assert_eq!(a.clone(), a);
}

#[test]
fn source_replacement_append_erase_and_clear_assertions() {
    let mut a = basic();
    a.set_element("hydrogen", 2.0, false).unwrap();
    assert_eq!(a.len(), 3);
    assert_eq!(a.mass_by_name("hydrogen").unwrap(), 2.0);
    a.set_element("carbon", 12.0, false).unwrap();
    assert_eq!(a.len(), 3);
    a.set_element("carbon", 12.0, true).unwrap();
    assert_eq!(a.len(), 4);
    assert_eq!(a.mass_by_name("carbon").unwrap(), 12.0);
    a.push_mass("carbon", 13.0).unwrap();
    a.push(IMSElement::from_mass("carbon", 14.0).unwrap())
        .unwrap();
    assert_eq!(a.mass(4).unwrap(), 13.0);
    assert_eq!(a.mass(5).unwrap(), 14.0);
    assert!(a.erase("carbon").unwrap());
    assert_eq!(a.mass_by_name("carbon").unwrap(), 13.0);
    for name in ["hydrogen", "oxygen", "nitrogen"] {
        assert!(a.erase(name).unwrap());
        assert!(!a.erase(name).unwrap());
    }
    a.clear();
    assert!(a.is_empty());
    assert!(!a.has_name("oxygen").unwrap());
    a.push_mass("new", 1.0).unwrap();
    assert_eq!(a.len(), 1);
}

#[test]
fn source_name_and_value_sort_assertions() {
    let mut a = basic();
    a.sort_by_names().unwrap();
    assert_eq!(
        a.elements()
            .iter()
            .map(IMSElement::name)
            .collect::<Vec<_>>(),
        ["hydrogen", "nitrogen", "oxygen"]
    );
    a.push_mass("carbon", 12.0).unwrap();
    a.sort_by_mass().unwrap();
    assert_eq!(
        a.elements()
            .iter()
            .map(IMSElement::name)
            .collect::<Vec<_>>(),
        ["hydrogen", "carbon", "nitrogen", "oxygen"]
    );
}

#[test]
fn source_parser_literals_and_file_replacement() {
    let mut p = IMSAlphabetTextParser::new();
    p.parse(Cursor::new(
        "# a comment which should be ignored\nA\t71.03711\nR\t156.10111\n",
    ))
    .unwrap();
    assert_eq!(p.elements().len(), 2);
    assert_eq!(p.elements()["A"], 71.03711);
    assert_eq!(p.elements()["R"], 156.10111);
    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("alphabet.txt");
    std::fs::write(
        &path,
        "# a comment which should be ignored\nhydrogen\t1.0\noxygen\t16.0\nnitrogen\t14.0\n",
    )
    .unwrap();
    let mut a = basic();
    let before = a.clone();
    assert!(a.load(directory.path().join("missing")).is_err());
    assert_eq!(a, before);
    a.load_with_parser(&path, &mut p).unwrap();
    assert_eq!(a.masses(0).unwrap(), [1.0, 14.0, 16.0]);
    assert_eq!(
        p.elements().keys().map(String::as_str).collect::<Vec<_>>(),
        ["hydrogen", "nitrogen", "oxygen"]
    );
    let mut direct = IMSAlphabet::new();
    direct.load(&path).unwrap();
    assert_eq!(direct, a);
    p.load(&path).unwrap();
    assert_eq!(p.elements().len(), 3);
    std::fs::write(&path, "").unwrap();
    a.load(&path).unwrap();
    assert!(a.is_empty());
}

struct Custom {
    elements: BTreeMap<String, f64>,
}
impl IMSAlphabetParser for Custom {
    fn parse(&mut self, _: &mut dyn BufRead) -> openms::Result<()> {
        self.elements.insert("hydrogen".into(), 1.0);
        self.elements.insert("oxygen".into(), 16.0);
        Ok(())
    }
    fn elements(&self) -> &BTreeMap<String, f64> {
        &self.elements
    }
    fn elements_mut(&mut self) -> &mut BTreeMap<String, f64> {
        &mut self.elements
    }
}
#[test]
fn source_custom_parser_extension_and_checked_mutable_maps() {
    let directory = TempDir::new(false).unwrap();
    let path = directory.path().join("data");
    std::fs::write(
        &path,
        "arbitrary custom format ignored by this source test parser",
    )
    .unwrap();
    let mut parser = Custom {
        elements: BTreeMap::new(),
    };
    assert!(parser.load(&directory.path().join("missing")).is_err());
    assert!(parser.elements.is_empty());
    let mut a = basic();
    a.load_with_parser(&path, &mut parser).unwrap();
    assert_eq!(a.masses(0).unwrap(), [1.0, 16.0]);
    let before = a.clone();
    parser.elements_mut().insert("bad".into(), f64::NAN);
    assert!(a.load_with_parser(&path, &mut parser).is_err());
    assert_eq!(a, before);
    assert!(parser.elements()["bad"].is_nan()); // callback state is not rolled back
    parser.elements_mut().remove("bad");
    parser
        .elements_mut()
        .insert("x".repeat(1024 * 1024 + 1), 1.0);
    assert!(a.load_with_parser(&path, &mut parser).is_err());
    assert_eq!(a, before);
}

#[test]
fn duplicate_name_operations_and_replacement_reset_full_element_state() {
    let mut first = IMSElement::from_mass("same", 8.0).unwrap();
    first.set_sequence("modified sequence").unwrap();
    let mut a =
        IMSAlphabet::from_elements(vec![first, IMSElement::from_mass("same", 9.0).unwrap()])
            .unwrap();
    a.set_element("same", 7.0, true).unwrap();
    assert_eq!(a.mass(0).unwrap(), 7.0);
    assert_eq!(a.element(0).unwrap().sequence(), "same");
    assert_eq!(a.mass(1).unwrap(), 9.0);
    assert!(a.erase("same").unwrap());
    assert_eq!(a.mass_by_name("same").unwrap(), 9.0);
    let before = a.clone();
    a.set_element("absent", f64::NAN, false).unwrap();
    assert_eq!(a, before);
    assert!(a.set_element("same", f64::NAN, false).is_err());
    assert_eq!(a, before);
    assert!(a.push_mass("invalid", f64::INFINITY).is_err());
    assert_eq!(a, before);
}

#[test]
fn stable_ties_and_lexical_parser_order_are_deterministic() {
    let a = IMSAlphabet::read(Cursor::new("z 2\na 2\nm 1\na 99\n")).unwrap();
    assert_eq!(
        a.elements()
            .iter()
            .map(IMSElement::name)
            .collect::<Vec<_>>(),
        ["m", "a", "z"]
    );
    assert_eq!(a.mass_by_name("a").unwrap(), 2.0);
    let mut a = IMSAlphabet::from_elements(vec![
        IMSElement::from_mass("z", -0.0).unwrap(),
        IMSElement::from_mass("a", 0.0).unwrap(),
        IMSElement::from_mass("a", 0.0).unwrap(),
    ])
    .unwrap();
    a.sort_by_mass().unwrap();
    assert_eq!(a.name(0).unwrap(), "z");
    a.sort_by_names().unwrap();
    assert_eq!(a.name(2).unwrap(), "z");
    assert!(a.mass(2).unwrap().is_sign_positive()); // source mass adds nominal +0
}

#[test]
fn portable_decimal_prefix_and_exact_comment_whitespace_policy() {
    let input = "\t # comment\n  \t\nA +.5 units\r\nB -2.5E+2#note\nC 1.25suffix\nD -1e-999\nλ 1e2\nE 3.5.7\n";
    let mut p = IMSAlphabetTextParser::new();
    p.parse(Cursor::new(input)).unwrap();
    assert_eq!(p.elements()["A"], 0.5);
    assert_eq!(p.elements()["B"], -250.0);
    assert_eq!(p.elements()["C"], 1.25);
    assert!(p.elements()["D"].is_sign_negative());
    assert_eq!(p.elements()["D"], 0.0);
    assert_eq!(p.elements()["λ"], 100.0);
    assert_eq!(p.elements()["E"], 3.5);
    p.parse(Cursor::new("\u{b}A 3\n")).unwrap(); // stream whitespace after the precheck
    assert_eq!(p.elements()["A"], 3.0);
    let before = p.clone();
    for line in ["\r\n", "\r#comment\n", "\u{b}#comment\n"] {
        assert!(p.parse(Cursor::new(line)).is_err());
        assert_eq!(p, before);
    }
}

#[test]
fn malformed_and_late_io_errors_preserve_parser_and_alphabet() {
    let mut p = IMSAlphabetTextParser::new();
    p.parse(Cursor::new("original 5\n")).unwrap();
    let before = p.clone();
    for value in [
        "", "1e", "1e+", "--1", ".", "1e999", "nan", "inf", "0x1p2", "+0Xff",
    ] {
        let input = format!("A 1\nbad {value}\n");
        let error = p.parse(Cursor::new(&input)).unwrap_err();
        assert!(matches!(error, openms::Error::Parse { line: 2, .. }));
        assert_eq!(p, before);
    }
    assert!(
        p.parse(Cursor::new([b'A', b' ', b'1', b'\n', 0xff]))
            .is_err()
    );
    assert_eq!(p, before);
    struct Broken {
        done: bool,
    }
    impl Read for Broken {
        fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
            if self.done {
                return Err(io::Error::other("late failure"));
            }
            self.done = true;
            out[..4].copy_from_slice(b"A 1\n");
            Ok(4)
        }
    }
    assert!(p.parse(BufReader::new(Broken { done: false })).is_err());
    assert_eq!(p, before);
}

#[test]
fn isotope_indices_averages_and_mass_sort_preflight() {
    let distribution = IMSIsotopeDistribution::from_peaks(
        vec![
            IMSIsotopePeak {
                mass: 0.125,
                abundance: 0.5,
            },
            IMSIsotopePeak {
                mass: 0.25,
                abundance: 0.5,
            },
        ],
        10,
    )
    .unwrap();
    let a = IMSAlphabet::from_elements(vec![
        IMSElement::from_distribution("two", distribution).unwrap(),
    ])
    .unwrap();
    assert_eq!(a.masses(1).unwrap(), [11.25]);
    assert_eq!(a.average_masses().unwrap(), [10.6875]);
    assert!(a.masses(2).is_err());
    let mut a = IMSAlphabet::from_elements(vec![IMSElement::new("empty", 16).unwrap()]).unwrap();
    a.sort_by_mass().unwrap(); // source std::sort never invokes comparator for one
    a.push_mass("mass", 1.0).unwrap();
    let before = a.clone();
    assert!(a.sort_by_mass().is_err());
    assert_eq!(a, before);
    assert_eq!(a.average_masses().unwrap(), [0.0, 1.0]);
}

#[test]
fn verbose_output_is_exact_and_validated_before_writing() {
    let a = IMSAlphabet::read(Cursor::new("B 2\nA 1\n")).unwrap();
    let options = IMSIsotopeOptions::default();
    let expected = "name:\tA\nsequence:\tA\nisotope distribution:\n\n\nname:\tB\nsequence:\tB\nisotope distribution:\n\n\n";
    assert_eq!(a.to_text(options).unwrap(), expected);
    let mut out = Vec::new();
    a.write(&mut out, options).unwrap();
    assert_eq!(out, expected.as_bytes());
    assert!(IMSAlphabet::read(Cursor::new(&out)).is_err()); // deliberately not a flat codec
    let mut out = vec![7, 8];
    assert!(
        a.write(
            &mut out,
            IMSIsotopeOptions {
                size: usize::MAX,
                abundances_sum_error: 0.0
            }
        )
        .is_err()
    );
    assert_eq!(out, [7, 8]);
}

#[test]
fn storage_input_and_shared_comparison_budgets_are_checked_atomically() {
    assert!(IMSAlphabet::from_elements(vec![IMSElement::default(); 100_001]).is_err());
    let mut p = IMSAlphabetTextParser::new();
    p.parse(Cursor::new("old 1\n")).unwrap();
    let before = p.clone();
    assert!(p.parse(Cursor::new("x".repeat(1024 * 1024 + 1))).is_err());
    assert_eq!(p, before);
    let prefix = "x".repeat(16_000);
    let elements = (0..1000)
        .rev()
        .map(|index| IMSElement::from_mass(&format!("{prefix}:{index:03}"), 1.0).unwrap())
        .collect();
    let mut alphabet = IMSAlphabet::from_elements(elements).unwrap();
    assert!(alphabet.sort_by_names().is_err());
    for index in 0..1000 {
        assert!(
            alphabet
                .name(index)
                .unwrap()
                .ends_with(&format!(":{:03}", 999 - index))
        );
    }
}
