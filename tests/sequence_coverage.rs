use openms::chemistry::{AASequence, SequenceCoverage};

fn sequence(text: &str) -> AASequence {
    AASequence::parse(text).unwrap()
}
fn coverage(protein: &str, peptides: &[&str]) -> f64 {
    let peptides: Vec<_> = peptides.iter().map(|p| sequence(p)).collect();
    SequenceCoverage::get_coverage(&sequence(protein), &peptides).unwrap()
}

#[test]
fn complete_pinned_source_class_test_examples() {
    // All four numerical/edge sections in SequenceCoverage_test.cpp.
    let value = coverage("ACDEFGHIK", &["ACD", "FGH"]);
    assert!((value - 66.6667).abs() < 0.0001);
    assert_eq!(value.to_bits(), (6.0_f64 * 100.0 / 9.0).to_bits());
    assert_eq!(coverage("ACDEFG", &[]), 0.0);
    assert_eq!(coverage("", &["ACD"]), 0.0);
    assert_eq!(coverage("ABCDE", &["ABC", "BCD"]), 80.0);
}

#[test]
fn all_occurrences_include_overlaps_and_union_duplicate_peptides() {
    assert_eq!(coverage("AAAAA", &["AAA"]), 100.0);
    assert_eq!(coverage("ACGACG", &["AC"]), 4.0 * 100.0 / 6.0);
    assert_eq!(
        coverage("ACGACG", &["AC", "AC", "", "GG"]),
        4.0 * 100.0 / 6.0
    );
    assert_eq!(coverage("ACDEFG", &["CDE", "DEF"]), 4.0 * 100.0 / 6.0);
    assert_eq!(coverage("ACDE", &["ACDE", "C"]), 100.0);
    assert_eq!(coverage("ACDE", &["ACDEA", "G", ""]), 0.0);
}

#[test]
fn annotations_and_unknown_chemistry_do_not_change_symbol_matching() {
    let protein = sequence("(Acetyl)AC(Carbamidomethyl)M(Oxidation)K");
    let before = protein.clone();
    let peptide = sequence("CM");
    assert_eq!(
        SequenceCoverage::get_coverage(&protein, &[peptide]).unwrap(),
        50.0
    );
    assert_eq!(protein, before);
    assert_eq!(
        SequenceCoverage::get_coverage(
            &sequence("ACMK"),
            &[sequence("C(Carbamidomethyl)M(Oxidation)")]
        )
        .unwrap(),
        50.0
    );
    let unknown = sequence("AXBZI");
    assert!(unknown.formula().is_err());
    assert_eq!(
        SequenceCoverage::get_coverage(&unknown, &[sequence("X[123.45]"), sequence("BZ")]).unwrap(),
        60.0
    );
    assert_eq!(coverage("IXXL", &["I"]), 25.0); // I and L are distinct.
    assert_eq!(coverage("ACGT", &["X"]), 0.0); // X is not a wildcard.
    assert_eq!(coverage("AXXT", &["X"]), 50.0);
}

// Different oracle: inspect each protein position and ask whether any complete
// peptide interval covering that position has all equal residues. No coverage
// bitmap or search-then-mark loop is shared with the implementation.
fn positional_oracle(protein: &str, peptides: &[String]) -> f64 {
    if protein.is_empty() {
        return 0.0;
    }
    let covered = (0..protein.len())
        .filter(|&position| {
            peptides.iter().any(|peptide| {
                (0..=position).any(|start| {
                    let end = start + peptide.len();
                    start <= position
                        && position < end
                        && end <= protein.len()
                        && protein.as_bytes()[start..end] == *peptide.as_bytes()
                })
            })
        })
        .count();
    covered as f64 * 100.0 / protein.len() as f64
}
#[test]
fn small_exhaustive_positional_oracle_and_peptide_order_invariance() {
    let alternatives = ["", "A", "G", "AA", "AG", "GA", "GG", "AGA", "GAG"];
    for length in 0..=7 {
        for bits in 0..1usize << length {
            let protein: String = (0..length)
                .map(|i| if bits & (1 << i) == 0 { 'A' } else { 'G' })
                .collect();
            let sequence = sequence(&protein);
            for first in 0..alternatives.len() {
                let second = (first * 3 + length) % alternatives.len();
                let strings: Vec<_> = [
                    alternatives[first],
                    alternatives[second],
                    alternatives[first],
                ]
                .iter()
                .map(|s| (*s).to_owned())
                .collect();
                let mut peptides: Vec<_> = strings
                    .iter()
                    .map(|s| AASequence::parse(s).unwrap())
                    .collect();
                let expected = positional_oracle(&protein, &strings);
                assert_eq!(
                    SequenceCoverage::get_coverage(&sequence, &peptides).unwrap(),
                    expected
                );
                peptides.reverse();
                assert_eq!(
                    SequenceCoverage::get_coverage(&sequence, &peptides).unwrap(),
                    expected
                );
            }
        }
    }
}

#[test]
fn conservative_search_budget_can_reject_even_when_comparisons_exit_early() {
    let protein = sequence(&"A".repeat(10_000));
    let peptide = sequence(&"G".repeat(5_000));
    assert_eq!(
        SequenceCoverage::get_coverage(&protein, std::slice::from_ref(&peptide)).unwrap(),
        0.0
    );
    // Two sets of 5,001 windows * (5,000 comparison bytes +1) exceed50M
    // before the coverage vector is allocated, although the first bytes differ.
    assert!(SequenceCoverage::get_coverage(&protein, &[peptide.clone(), peptide]).is_err());
    assert_eq!(
        SequenceCoverage::get_coverage(&protein, &[sequence("A")]).unwrap(),
        100.0
    );
}
