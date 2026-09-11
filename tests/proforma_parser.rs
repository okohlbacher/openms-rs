use openms::chemistry::proforma::*;

fn syntax(error: ParseFailure) -> Box<ParseError> {
    match error {
        ParseFailure::Syntax(error) => error,
        ParseFailure::Resource(error) => panic!("unexpected resource failure: {error}"),
    }
}
fn element(pf: &Peptidoform, index: usize) -> &SequenceElement {
    let SequenceSection::Element(element) = &pf.sequence[index] else {
        panic!("not an ordinary element")
    };
    element
}
fn tag(pf: &Peptidoform, index: usize) -> &ModificationTag {
    &element(pf, index).modifications[0].alternatives[0].0
}
fn needs_ion(input: &str) -> bool {
    // Exact source positive-fixture dispatch heuristic (not a native parser).
    if input.contains("//") {
        return true;
    }
    let mut depth = 0i32;
    let bytes = input.as_bytes();
    for (i, byte) in bytes.iter().copied().enumerate() {
        match byte {
            b'[' => depth += 1,
            b']' => depth -= 1,
            b'+' if depth == 0 && i > 0 => {
                if !matches!(bytes[i - 1], b'[' | b':' | b'|') {
                    return true;
                }
            }
            b'/' if depth == 0
                && bytes
                    .get(i + 1)
                    .is_some_and(|b| b.is_ascii_digit() || matches!(b, b'+' | b'-' | b'[')) =>
            {
                return true;
            }
            _ => (),
        }
    }
    false
}

#[test]
fn all_176_source_positive_and_22_single_chain_negative_cases() {
    let positive = include_str!("data/proforma_parser_positive.txt")
        .lines()
        .filter(|s| !s.is_empty() && !s.starts_with('#'))
        .collect::<Vec<_>>();
    assert_eq!(positive.len(), 176);
    for case in positive {
        if needs_ion(case) {
            PeptidoformIon::parse(case).unwrap_or_else(|e| panic!("source ion {case:?}: {e}"));
        } else {
            Peptidoform::parse(case).unwrap_or_else(|e| panic!("source peptide {case:?}: {e}"));
        }
    }
    let negative = include_str!("data/proforma_parser_negative.txt")
        .lines()
        .filter(|s| !s.is_empty() && !s.starts_with('#'))
        .collect::<Vec<_>>();
    assert_eq!(negative.len(), 22);
    for case in negative {
        assert!(
            Peptidoform::parse(case).is_err(),
            "source negative through parse, not parseIon: {case}"
        );
    }
}

#[test]
fn source_component_grammar_success_and_failure_cases() {
    // Source ProFormaParser_test.cpp:2407–2528, 42+14 successes and3+3 failures.
    let plain = [
        "A",
        "U",
        "X",
        "A[#g1]",
        "Asn",
        "A[Formula:N1H3]",
        "A[Formula:NH3]",
        "A[Formula:C-1O-1]",
        "A[Formula:H1H1H1N1]",
        "A[Formula:[15N]H3]",
        "A[Formula:[15N1]H3]",
        "A[Formula:C12H20O2]",
        "A[Formula:HN-1O2]",
        "A[Formula:[13C2][12C-2]H2N]",
        "A[Formula:[13C2]C-2H2N]",
        "A[Formula:Zn1]",
        "A[Formula:Zn:z+2]",
        "A[Formula:Zn1:z+2]",
        "A[Formula:Na:z+1]",
        "A[+23]",
        "A[-23.0]",
        "A[+23.0]",
        "A[U:+23]",
        "A[Obs:+23]",
        "A[+1]",
        "A[U:+1]",
        "A[+23.092]",
        "A[-23.092]",
        "A[14|Obs:+14|UNIMOD:0034|U:Methyl]",
        "A[Formula:C-1O-1]",
        "A[Glycan:Hex1HexNAc2]",
        "A[Formula:AlH-3:z+1]",
        "A[Formula:H-1:z-1]",
        "<[TMT6plex]@K,N-term>A",
        "<[TMT6plex]@K,N-terM>A",
        "<[TMT6plex]@K,n-terM>A",
        "<[TMT6plex]@K,n-tErM>A",
        "<[TMT6plex]@K,N-term:A>AA",
        "<[TMT6plex]@K,N-term:A,N-term:B>AA",
        "A[Glycan:Hex2HexNAc]",
        "A[Glycan:Hex2HexNAc1]",
        "III",
    ];
    assert_eq!(plain.len(), 42);
    for text in plain {
        Peptidoform::parse(text).unwrap_or_else(|e| panic!("{text}: {e}"));
    }
    let ions = [
        "A/0",
        "A/1",
        "A/2",
        "A/1",
        "A/+1",
        "A/-1",
        "A/2",
        "A/+2",
        "A/-10101",
        "A/[Na:z+1]",
        "A/[Na:z+1^2]",
        "A/[Na:z+1,Zn:z+2,H:z+1]",
        "A/[Na:z+1]",
        "A/[Na:z+1^2]",
    ];
    assert_eq!(ions.len(), 14);
    for text in ions {
        PeptidoformIon::parse(text).unwrap_or_else(|e| panic!("{text}: {e}"));
    }
    for text in ["A[#g1", "A[Formula:[15NH3]", "(KKK"] {
        assert!(Peptidoform::parse(text).is_err());
    }
    for text in ["A/+1i", "A/1/1", "A/[Na^1]"] {
        assert!(PeptidoformIon::parse(text).is_err());
    }
}

#[test]
fn source_cv_named_mass_formula_and_glycan_fields() {
    for (text, db, id) in [
        ("A[UNIMOD:35]", CvDatabase::Unimod, "35"),
        ("A[MOD:00046]", CvDatabase::Mod, "00046"),
        ("A[RESID:AA0037]", CvDatabase::Resid, "AA0037"),
        ("A[GNO:G59626AS]", CvDatabase::Gno, "G59626AS"),
        ("A[XLMOD:02001]", CvDatabase::Xlmod, "02001"),
    ] {
        let pf = Peptidoform::parse(text).unwrap();
        assert_eq!(
            tag(&pf, 0),
            &ModificationTag::CvAccession(CvAccession {
                database: db,
                accession: id.into()
            })
        );
    }
    for (text, hint) in [
        ("A[U:Oxidation]", CvDatabase::Unimod),
        ("A[M:Oxidation]", CvDatabase::Mod),
        ("A[R:Oxidation]", CvDatabase::Resid),
        ("A[X:Oxidation]", CvDatabase::Xlmod),
        ("A[G:Oxidation]", CvDatabase::Gno),
    ] {
        assert_eq!(
            tag(&Peptidoform::parse(text).unwrap(), 0),
            &ModificationTag::NamedMod(NamedMod {
                name: "Oxidation".into(),
                cv_hint: Some(hint)
            })
        );
    }
    for (text, source) in [
        ("A[+15.9949]", MassDeltaSource::None),
        ("A[U:+15.9949]", MassDeltaSource::U),
        ("A[Obs:+15.9949]", MassDeltaSource::Obs),
    ] {
        assert_eq!(
            tag(&Peptidoform::parse(text).unwrap(), 0),
            &ModificationTag::MassDelta(MassDelta {
                source,
                mass: 15.9949,
                original_text: "+15.9949".into()
            })
        );
    }
    let pf = Peptidoform::parse("SEQUEN[Formula:[13C2]C-2H2N:z+2]CE").unwrap();
    assert_eq!(
        tag(&pf, 5),
        &ModificationTag::FormulaTag(FormulaTag {
            formula_string: "[13C2]C-2H2N".into(),
            charge: Some(2)
        })
    );
    let pf = Peptidoform::parse("SEQUEN[Glycan:HexNAc1Hex2]CE").unwrap();
    assert_eq!(
        tag(&pf, 5),
        &ModificationTag::GlycanComposition(GlycanComposition {
            components: vec![
                (GlycanComponent::Name("HexNAc".into()), 1),
                (GlycanComponent::Name("Hex".into()), 2)
            ]
        })
    );
    for name in [
        "Cation:Mg[II]",
        "half cystine",
        "O-phospho-L-serine",
        "Glu->pyro-Glu",
        "Label:13C",
    ] {
        if name.contains("Label:") {
            assert!(Peptidoform::parse(&format!("A[{name}]")).is_err());
            continue;
        }
        assert_eq!(
            tag(&Peptidoform::parse(&format!("A[{name}]")).unwrap(), 0),
            &ModificationTag::NamedMod(NamedMod {
                name: name.into(),
                cv_hint: None
            })
        );
    }
}

#[test]
fn source_terminal_prefix_range_ambiguity_and_alternative_structure() {
    let pf=Peptidoform::parse("<13C><15N><[TMT6plex]@K,N-term>[Phospho]^2?{Glycan:Hex}[Acetyl][Methyl]-A(?DQ[+1])(MK[Oxidation])[+2][Phospho]-[Amidated][Methyl]").unwrap();
    assert_eq!(pf.global_mods.len(), 3);
    assert_eq!(pf.unlocalised_mods[0].occurrence, Some(2));
    assert_eq!(pf.labile_mods.len(), 1);
    assert_eq!(pf.n_term_mods.len(), 2);
    assert_eq!(pf.c_term_mods.len(), 2);
    assert_eq!(pf.sequence.len(), 3);
    let SequenceSection::AmbiguousRegion(region) = &pf.sequence[1] else {
        panic!()
    };
    assert_eq!(
        region
            .elements
            .iter()
            .map(|e| e.amino_acid)
            .collect::<String>(),
        "DQ"
    );
    assert_eq!(region.elements[1].modifications.len(), 1);
    let SequenceSection::ModifiedRange(range) = &pf.sequence[2] else {
        panic!()
    };
    assert_eq!(range.elements[1].modifications.len(), 1);
    assert_eq!(range.modifications.len(), 2);
    let pf = Peptidoform::parse("A[Oxidation|+15.99|Position:N-term,C-term,MK|INFO:note]").unwrap();
    let alternatives = &element(&pf, 0).modifications[0].alternatives;
    assert_eq!(alternatives.len(), 4);
    assert_eq!(
        alternatives[2].0,
        ModificationTag::PositionConstraint(PositionConstraint {
            residues: vec!['M', 'K'],
            n_term: true,
            c_term: true
        })
    );
    assert_eq!(
        alternatives[3].0,
        ModificationTag::InfoTag(InfoTag {
            text: "note".into()
        })
    );
    assert!(element(&pf, 0).modifications[0].resolved_mod.is_none());
}

#[test]
fn source_names_crosslinks_labels_and_charge_context() {
    let ion = PeptidoformIon::parse("(>Trypsin)EMEVEESPEK/2+(>Keratin)ELVISLIVER/3").unwrap();
    assert!(ion.is_chimeric);
    assert_eq!(ion.chains[0].name.as_deref(), Some("Trypsin"));
    assert_eq!(ion.chains[1].name.as_deref(), Some("Keratin"));
    assert_eq!(ion.chains[0].charge, Some(ChargeState::Simple(2)));
    assert_eq!(ion.chains[1].charge, Some(ChargeState::Simple(3)));
    assert!(ion.charge.is_none());
    let ion = PeptidoformIon::parse("A[XLMOD:02001#XL1]//K[#XL1]/4").unwrap();
    assert!(!ion.is_chimeric);
    assert_eq!(ion.charge, Some(ChargeState::Simple(4)));
    assert!(ion.chains.iter().all(|c| c.charge.is_none()));
    for (label, kind) in [
        ("XL1", LabelType::Crosslink),
        ("BRANCH", LabelType::Branch),
        ("xl1", LabelType::Ambiguous),
        ("17", LabelType::Ambiguous),
    ] {
        let pf = Peptidoform::parse(&format!("A[#{label}(0.90)]")).unwrap();
        let (tag, label) = &element(&pf, 0).modifications[0].alternatives[0];
        assert_eq!(tag, &ModificationTag::InfoTag(InfoTag::default()));
        assert_eq!(label.as_ref().unwrap().label_type, kind);
        assert_eq!(label.as_ref().unwrap().score, Some(0.9));
    }
    let ion = PeptidoformIon::parse("A/[Na:z+1^2,Zn:z+2,H:z-1]").unwrap();
    assert_eq!(
        ion.charge,
        Some(ChargeState::Adducts(vec![
            AdductIon {
                formula: "Na".into(),
                charge: 1,
                occurrence: Some(2)
            },
            AdductIon {
                formula: "Zn".into(),
                charge: 2,
                occurrence: None
            },
            AdductIon {
                formula: "H".into(),
                charge: -1,
                occurrence: None
            }
        ]))
    );
}

#[test]
fn independent_source_permissive_and_non_roundtrip_boundaries() {
    let pf = Peptidoform::parse("(>>one)(>two(x))A").unwrap();
    assert_eq!(pf.name.as_deref(), Some("one / two(x)"));
    assert!(Peptidoform::parse("(>)A").unwrap().name.is_none());
    let ion = PeptidoformIon::parse("A/1+B//C/3").unwrap();
    assert!(ion.is_chimeric);
    assert_eq!(ion.charge, Some(ChargeState::Simple(3)));
    assert_eq!(ion.chains[0].charge, Some(ChargeState::Simple(1)));
    assert!(ion.name.is_none());
    let pf = Peptidoform::parse("[x][y]^2?A").unwrap();
    assert_eq!(pf.unlocalised_mods.len(), 2);
    assert_eq!(pf.unlocalised_mods[0].modifications.len(), 1);
    assert_eq!(pf.unlocalised_mods[1].occurrence, Some(2));
    assert_eq!(pf.to_text(WriteMode::Lossless).unwrap(), "[x]?[y]^2?A");
    for text in [
        "(?)",
        "(ABC)",
        "A[foo(bar]",
        "A[Formula:)]",
        "A[UNIMOD:-2.5]",
        "<[x]@>A",
        "<[x]@K,>A",
        "A[#XL1]",
        "A[#XL1]//B[#XL2]",
        "A[Glycan:Unknown-2Hex0]",
        "A[Cation:]",
        "Asn",
    ] {
        if text.contains("//") {
            PeptidoformIon::parse(text).unwrap();
        } else {
            Peptidoform::parse(text).unwrap_or_else(|e| panic!("{text}: {e}"));
        }
    }
    for text in [
        "()",
        "A[]",
        "A[x|]",
        "A[|x]",
        "A[Glycan:]",
        "A[Glycan:Formula:C]",
        "{INFO:x}A",
        "A[Position:]",
        "A[1e3]",
        "A/+1",
        "A[Unimod:1]",
        "A[Formula:]",
        "[x]A",
    ] {
        assert!(Peptidoform::parse(text).is_err(), "{text}");
    }
    let ion = PeptidoformIon::parse("A/[Na:z+1]").unwrap();
    assert_eq!(ion.to_text(WriteMode::Lossless).unwrap(), "A/[Na:z+1]1+");
    assert!(PeptidoformIon::parse(&ion.to_text(WriteMode::Lossless).unwrap()).is_err());
}

#[test]
fn decimal_prefix_sign_overflow_and_portable_subnormal_rules() {
    for (text, expected) in [
        ("A/2.5", 2),
        ("A/--1", 1),
        ("A/+-1", -1),
        ("A/-+1", -1),
        ("A/-2147483648", i32::MIN),
        ("A/2147483647", i32::MAX),
    ] {
        assert_eq!(
            PeptidoformIon::parse(text).unwrap().charge,
            Some(ChargeState::Simple(expected))
        );
    }
    let pf = Peptidoform::parse("[Phospho]^2.5?A[Glycan:Hex-2.5]").unwrap();
    assert_eq!(pf.unlocalised_mods[0].occurrence, Some(2));
    assert_eq!(
        tag(&pf, 0),
        &ModificationTag::GlycanComposition(GlycanComposition {
            components: vec![(GlycanComponent::Name("Hex".into()), -2)]
        })
    );
    for text in [
        "A/2147483648",
        "A/--2147483648",
        "A/.5",
        "A/[H:z2147483648]",
        "A/[H:z+1^2147483648]",
    ] {
        assert!(PeptidoformIon::parse(text).is_err(), "{text}");
    }
    for text in ["A[++1]", "A[-+1]", "A[+.]", "[x]^.5?A", "A[Glycan:Hex.5]"] {
        assert!(Peptidoform::parse(text).is_err(), "{text}");
    }
    let tiny = format!("0.{}5", "0".repeat(323));
    let pf = Peptidoform::parse(&format!("A[+{tiny}]")).unwrap();
    let ModificationTag::MassDelta(mass) = tag(&pf, 0) else {
        panic!()
    };
    assert_eq!(mass.mass.to_bits(), 1);
    assert_eq!(mass.original_text, format!("+{tiny}"));
    let below = format!("A[0.{}1]", "0".repeat(400));
    assert_eq!(
        syntax(Peptidoform::parse(&below).unwrap_err()).code(),
        ErrorCode::InvalidMassValue
    );
    let overflow = format!("A[{}]", "9".repeat(400));
    assert_eq!(
        syntax(Peptidoform::parse(&overflow).unwrap_err()).code(),
        ErrorCode::InvalidMassValue
    );
    let pf = Peptidoform::parse("A[-0.000]").unwrap();
    let ModificationTag::MassDelta(mass) = tag(&pf, 0) else {
        panic!()
    };
    assert_eq!(mass.mass.to_bits(), (-0.0f64).to_bits());
}

#[test]
fn unicode_is_preserved_except_unrepresentable_individual_source_fields() {
    let pf = Peptidoform::parse("(>αβ)A[é]|x");
    assert!(pf.is_err());
    for text in [
        "(>αβ)A[é]",
        "A[INFO:αβ]",
        "A[Formula:é]",
        "A[GNO:é]",
        "A[U: é]",
    ] {
        let pf = Peptidoform::parse(text).unwrap();
        assert_eq!(pf.to_text(WriteMode::Lossless).unwrap(), text);
    }
    let error = syntax(Peptidoform::parse("A[Glycan:é]").unwrap_err());
    assert_eq!(error.code(), ErrorCode::UnexpectedCharacter);
    assert_eq!(error.position(), 9);
    assert!(error.message().contains("invalid UTF-8"));
    for text in ["é", "A[Position:é]", "A[foo]é"] {
        assert_eq!(
            syntax(Peptidoform::parse(text).unwrap_err()).code(),
            ErrorCode::InvalidAminoAcid
        );
    }
}

#[test]
fn actual_source_diagnostic_codes_positions_and_empty_expected_found() {
    for (text, code, position) in [
        ("A[+1", ErrorCode::UnexpectedCharacter, 4),
        ("", ErrorCode::EmptySequence, 0),
        ("A1", ErrorCode::UnexpectedCharacter, 1),
        ("A[]", ErrorCode::UnexpectedCharacter, 2),
        ("A[Formula:]", ErrorCode::InvalidFormula, 10),
        ("A[+.]", ErrorCode::InvalidMassValue, 3),
    ] {
        let error = syntax(Peptidoform::parse(text).unwrap_err());
        assert_eq!((error.code(), error.position()), (code, position), "{text}");
        assert_eq!((error.expected(), error.found()), ("", ""));
    }
    let error = syntax(Peptidoform::parse("A[+1").unwrap_err());
    assert_eq!(error.message(), "Expected ']'");
    assert_eq!(
        error.formatted_message().unwrap(),
        "ProForma parse error at position 4: Unexpected character\nContext: A[+1>>><END OF INPUT><<<"
    );
}

#[test]
fn complete_structured_error_api_and_source_context_literals() {
    for (code, text) in [
        (ErrorCode::UnexpectedCharacter, "Unexpected character"),
        (ErrorCode::UnclosedBracket, "Unclosed bracket"),
        (ErrorCode::UnmatchedBracket, "Unmatched closing bracket"),
        (
            ErrorCode::InvalidCvPrefix,
            "Invalid controlled vocabulary prefix",
        ),
        (ErrorCode::InvalidCvAccession, "Invalid CV accession number"),
        (ErrorCode::InvalidAminoAcid, "Invalid amino acid"),
        (ErrorCode::InvalidMassValue, "Invalid mass value"),
        (ErrorCode::InvalidFormula, "Invalid chemical formula"),
        (ErrorCode::UnknownMonosaccharide, "Unknown monosaccharide"),
        (
            ErrorCode::DanglingCrosslinkLabel,
            "Dangling crosslink label",
        ),
        (ErrorCode::EmptySequence, "Empty sequence"),
        (ErrorCode::InvalidCharge, "Invalid charge state"),
        (
            ErrorCode::InvalidOccurrenceSpecifier,
            "Invalid occurrence specifier",
        ),
        (ErrorCode::UnexpectedEndOfInput, "Unexpected end of input"),
        (ErrorCode::InternalError, "Internal parser error"),
    ] {
        assert_eq!(code.as_str(), text);
    }
    let error = ParseError::new(ErrorCode::UnexpectedCharacter, 100, "ABC", "Test error").unwrap();
    assert_eq!(error.position(), 3);
    assert!(error.formatted_message().unwrap().contains("END OF INPUT"));
    let error = ParseError::new(
        ErrorCode::UnexpectedEndOfInput,
        6,
        "ABCDEF",
        "Unexpected end",
    )
    .unwrap();
    assert_eq!(error.context_before(), b"ABCDEF");
    assert!(error.context_after().is_empty());
    let error = ParseError::new(ErrorCode::EmptySequence, 0, "", "Empty sequence").unwrap();
    assert!(error.context_before().is_empty());
    assert!(error.context_after().is_empty());
    let mut error = ParseError::new(
        ErrorCode::UnclosedBracket,
        3,
        "PEP[",
        "specific custom text",
    )
    .unwrap();
    error.set_expected_found("']'", "end of input").unwrap();
    assert_eq!(error.message(), "specific custom text");
    assert_eq!(
        error.formatted_message().unwrap(),
        "ProForma parse error at position 3: Unclosed bracket\nContext: PEP>>>[<<<\nExpected: ']'\nFound: end of input"
    );
    let saved = error.clone();
    assert!(error.set_expected_found(&"x".repeat(65536), "x").is_err());
    assert_eq!(error, saved);
    assert!(ParseError::new(ErrorCode::InternalError, 0, "", &"x".repeat(65537)).is_err());
    let input = "012345678901234567890123456789012345678901234567890123456789";
    let error = ParseError::new(ErrorCode::UnexpectedCharacter, 25, input, "").unwrap();
    assert_eq!(error.context_before(), &input.as_bytes()[5..25]);
    assert_eq!(error.context_after(), &input.as_bytes()[25..45]);
    assert!(
        error
            .formatted_message()
            .unwrap()
            .contains("Context: ...56789012345678901234>>>5<<<6789012345678901234...")
    );
    let error = ParseError::new(ErrorCode::UnexpectedCharacter, 1, "é", "").unwrap();
    assert_eq!(error.position(), 1);
    assert_eq!(error.context_before(), &[0xc3]);
    assert_eq!(error.context_after(), &[0xa9]);
    assert!(error.formatted_message().unwrap().contains('�'));
}

#[test]
fn input_depth_and_lookahead_failures_are_checked() {
    assert!(matches!(
        Peptidoform::parse(&"A".repeat(4 * 1024 * 1024 + 1)),
        Err(ParseFailure::Resource(_))
    ));
    let name = format!("A[foo{}x{}]", "(".repeat(257), ")".repeat(257));
    assert!(matches!(
        Peptidoform::parse(&name),
        Err(ParseFailure::Resource(_))
    ));
    let name = format!("A[foo{}x{}]", "(".repeat(256), ")".repeat(256));
    assert!(Peptidoform::parse(&name).is_ok());
    let many = format!("{}A", "[x]".repeat(200_000));
    assert!(Peptidoform::parse(&many).is_err());
    // Repeated absent delimiter lookahead remains bounded and never publishes.
    let unfinished = format!("{}A", "[".repeat(257));
    assert!(matches!(
        Peptidoform::parse(&unfinished),
        Err(ParseFailure::Resource(_))
    ));
}
