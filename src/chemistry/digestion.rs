// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Native registry and full/semi/nonspecific protein digestion, based on the
//! pinned OpenMS4-core `7c029e8` enzyme data and digestion algorithms.
//! See `docs/DIGESTION_SUPPORT.md` for ordering, compatibility and limitations.

use super::{AASequence, EmpiricalFormula};
use crate::{Error, Result};
use std::fmt;
use std::ops::Range;
use std::str::FromStr;

mod enzymes;
pub use enzymes::Protease;

pub const MAX_DIGESTION_SEQUENCE_LENGTH: usize = 1_000_000;
pub const MAX_DIGESTED_RESIDUES: usize = 10_000_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CleavageRule {
    Never,
    Boundary {
        previous: &'static [u8],
        previous_two: &'static [u8],
        next: &'static [u8],
        forbidden_next: &'static [u8],
    },
}
impl CleavageRule {
    fn matches(self, sequence: &[u8], position: usize) -> bool {
        let Self::Boundary {
            previous,
            previous_two,
            next,
            forbidden_next,
        } = self
        else {
            return false;
        };
        // Only internal positions are passed here; two-residue lookbehind is
        // explicitly guarded at the first bond.
        (previous.is_empty() || previous.contains(&sequence[position - 1]))
            && (previous_two.is_empty()
                || (position >= 2 && previous_two.contains(&sequence[position - 2])))
            && (next.is_empty() || next.contains(&sequence[position]))
            && !forbidden_next.contains(&sequence[position])
    }
}

/// Immutable pinned enzyme metadata. Empty gain formulas/IDs remain empty when
/// the source does not declare them; peptide digestion itself uses subsequences.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DigestionEnzymeProtein {
    protease: Protease,
    name: &'static str,
    regex: &'static str,
    description: &'static str,
    synonyms: &'static [&'static str],
    n_term_gain: &'static str,
    c_term_gain: &'static str,
    psi_id: &'static str,
    xtandem_id: &'static str,
    comet_id: Option<i32>,
    msgf_id: Option<i32>,
    omssa_id: Option<i32>,
    rule: CleavageRule,
}
impl DigestionEnzymeProtein {
    pub fn protease(&self) -> Protease {
        self.protease
    }
    pub fn name(&self) -> &'static str {
        self.name
    }
    /// Source expression, exposed as metadata; it is not executed at runtime.
    pub fn regex(&self) -> &'static str {
        self.regex
    }
    pub fn description(&self) -> &'static str {
        self.description
    }
    pub fn synonyms(&self) -> &'static [&'static str] {
        self.synonyms
    }
    pub fn n_term_gain(&self) -> EmpiricalFormula {
        EmpiricalFormula::parse(self.n_term_gain).expect("validated pinned enzyme formula")
    }
    pub fn c_term_gain(&self) -> EmpiricalFormula {
        EmpiricalFormula::parse(self.c_term_gain).expect("validated pinned enzyme formula")
    }
    pub fn psi_id(&self) -> &'static str {
        self.psi_id
    }
    pub fn xtandem_id(&self) -> &'static str {
        self.xtandem_id
    }
    pub fn comet_id(&self) -> Option<i32> {
        self.comet_id
    }
    pub fn msgf_id(&self) -> Option<i32> {
        self.msgf_id
    }
    pub fn omssa_id(&self) -> Option<i32> {
        self.omssa_id
    }
}

/// Immutable registry. Names and declared synonyms are matched case-sensitively.
#[derive(Clone, Copy, Debug, Default)]
pub struct ProteaseDB;
impl ProteaseDB {
    pub fn global() -> &'static Self {
        &Self
    }
    /// Stable pinned table order, including named entries with identical rules.
    pub fn enzymes(&self) -> &'static [DigestionEnzymeProtein] {
        &enzymes::ENZYMES
    }
    pub fn has_enzyme(&self, name: &str) -> bool {
        self.get_enzyme(name).is_ok()
    }
    pub fn get_enzyme(&self, name: &str) -> Result<&'static DigestionEnzymeProtein> {
        self.enzymes()
            .iter()
            .find(|e| e.name == name || e.synonyms.contains(&name))
            .ok_or_else(|| invalid(format!("unknown protease {name:?}")))
    }
    pub fn names(&self) -> Vec<&'static str> {
        self.enzymes().iter().map(|e| e.name).collect()
    }
    pub fn has_regex(&self, regex: &str) -> bool {
        self.enzymes().iter().any(|e| e.regex == regex)
    }
    /// Lookup of exact pinned expressions, not an arbitrary regex compiler.
    /// Multiple names can share a rule, so all matching entries are returned.
    pub fn enzymes_by_regex(&self, regex: &str) -> Result<Vec<&'static DigestionEnzymeProtein>> {
        let matches: Vec<_> = self.enzymes().iter().filter(|e| e.regex == regex).collect();
        if matches.is_empty() {
            Err(Error::Unsupported(
                "only cleavage expressions present in the pinned enzyme registry are supported"
                    .into(),
            ))
        } else {
            Ok(matches)
        }
    }
    pub fn xtandem_names(&self) -> Vec<&'static str> {
        self.enzymes()
            .iter()
            .filter(|e| !e.xtandem_id.is_empty())
            .map(|e| e.name)
            .collect()
    }
    pub fn comet_names(&self) -> Vec<&'static str> {
        self.enzymes()
            .iter()
            .filter(|e| e.comet_id.is_some())
            .map(|e| e.name)
            .collect()
    }
    pub fn msgf_names(&self) -> Vec<&'static str> {
        self.enzymes()
            .iter()
            .filter(|e| e.msgf_id.is_some())
            .map(|e| e.name)
            .collect()
    }
    pub fn omssa_names(&self) -> Vec<&'static str> {
        self.enzymes()
            .iter()
            .filter(|e| e.omssa_id.is_some())
            .map(|e| e.name)
            .collect()
    }
}
impl Protease {
    pub fn from_name(name: &str) -> Result<Self> {
        Ok(ProteaseDB::global().get_enzyme(name)?.protease)
    }
    pub fn metadata(self) -> &'static DigestionEnzymeProtein {
        enzymes::ENZYMES
            .iter()
            .find(|e| e.protease == self)
            .expect("every protease has a pinned definition")
    }
    pub fn name(self) -> &'static str {
        self.metadata().name
    }
    /// All internal cleavage positions plus protein boundaries 0 and len, without
    /// duplicate boundaries. Empty input returns `[0]`. Accepts uppercase A-Z.
    pub fn cleavage_sites(self, sequence: &str) -> Result<Vec<usize>> {
        validate_sequence(sequence)?;
        let mut cuts = vec![0];
        let rule = self.metadata().rule;
        for position in 1..sequence.len() {
            if rule.matches(sequence.as_bytes(), position) {
                cuts.push(position);
            }
        }
        if !sequence.is_empty() {
            cuts.push(sequence.len());
        }
        Ok(cuts)
    }
    /// Protein termini are valid sites for every enzyme, including no cleavage.
    pub fn is_cleavage_site(self, sequence: &str, position: usize) -> Result<bool> {
        validate_sequence(sequence)?;
        Ok(position <= sequence.len()
            && (position == 0
                || position == sequence.len()
                || self.metadata().rule.matches(sequence.as_bytes(), position)))
    }
    pub fn count_internal_cleavage_sites(self, sequence: &str) -> Result<usize> {
        Ok(self.cleavage_sites(sequence)?.len().saturating_sub(2))
    }
}
impl FromStr for Protease {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        Self::from_name(s)
    }
}
impl fmt::Display for Protease {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DigestionSpecificity {
    /// Both termini must be protein boundaries or cleavage sites.
    #[default]
    Full,
    /// At least one terminus must be a protein boundary or cleavage site.
    Semi,
    /// Enumerate every substring in the length interval, without enzyme limits.
    None,
}
impl DigestionSpecificity {
    pub fn name(self) -> &'static str {
        match self {
            Self::Full => "full",
            Self::Semi => "semi",
            Self::None => "none",
        }
    }
}
impl FromStr for DigestionSpecificity {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        match s {
            "full" => Ok(Self::Full),
            "semi" => Ok(Self::Semi),
            "none" => Ok(Self::None),
            _ => Err(Error::Unsupported(format!(
                "unsupported digestion specificity {s:?}"
            ))),
        }
    }
}
impl fmt::Display for DigestionSpecificity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// Validation-only allowances. Defaults match ProteaseDigestion::isValidProduct.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ProductValidation {
    pub ignore_missed_cleavages: bool,
    /// Treat start 1 or 2 as start 0 when the protein begins with M. Counts then
    /// include any enzyme sites in the restored prefix, matching the source.
    pub allow_nterm_protein_cleavage: bool,
    /// D|P bonds can satisfy termini; they are not extra missed-cleavage sites.
    pub allow_random_asp_pro_cleavage: bool,
}
impl Default for ProductValidation {
    fn default() -> Self {
        Self {
            ignore_missed_cleavages: true,
            allow_nterm_protein_cleavage: false,
            allow_random_asp_pro_cleavage: false,
        }
    }
}

/// Digestion settings. Defaults preserve the original native trypsin API.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProteaseDigestion {
    pub enzyme: Protease,
    pub missed_cleavages: usize,
    pub min_length: usize,
    pub max_length: Option<usize>,
    pub specificity: DigestionSpecificity,
    /// Maximum accepted products for a ranges/count/digest operation.
    pub max_products: usize,
    /// Bound sequence scan plus candidate range checks, including length rejects.
    pub max_work: usize,
}
impl Default for ProteaseDigestion {
    fn default() -> Self {
        Self {
            enzyme: Protease::Trypsin,
            missed_cleavages: 0,
            min_length: 1,
            max_length: None,
            specificity: DigestionSpecificity::Full,
            max_products: 100_000,
            max_work: 10_000_000,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DigestionProduct {
    pub start: usize,
    pub end: usize,
    /// Actual internal enzyme sites, including in unrestricted products.
    pub missed_cleavages: usize,
}
impl DigestionProduct {
    pub fn range(self) -> Range<usize> {
        self.start..self.end
    }
    pub fn len(self) -> usize {
        self.end.saturating_sub(self.start)
    }
    pub fn is_empty(self) -> bool {
        self.start >= self.end
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DigestedPeptide {
    pub sequence: AASequence,
    /// Zero-based inclusive start in the input protein.
    pub start: usize,
    /// Zero-based exclusive end in the input protein.
    pub end: usize,
    pub missed_cleavages: usize,
}

impl ProteaseDigestion {
    pub fn set_enzyme(&mut self, name: &str) -> Result<()> {
        self.enzyme = Protease::from_name(name)?;
        Ok(())
    }
    fn validate(&self) -> Result<()> {
        if self.min_length == 0
            || self.max_length.is_some_and(|max| max < self.min_length)
            || self.max_products == 0
            || self.max_work == 0
        {
            return Err(invalid(
                "digest requires 1 <= min_length <= max_length and positive resource limits",
            ));
        }
        Ok(())
    }
    /// Fully specific products are ordered by missed cleavages then position;
    /// semi-specific variants follow in the source's alternating-end order.
    /// Unrestricted products are ordered by start then length. Modified residues
    /// and original terminal modifications survive applicable subsequences.
    pub fn digest(&self, protein: &AASequence) -> Result<Vec<DigestedPeptide>> {
        let ranges = self.digest_ranges(protein.as_str())?;
        let mut residue_count = 0_usize;
        for product in &ranges {
            residue_count = residue_count
                .checked_add(product.len())
                .ok_or_else(|| invalid("digested residue count overflows"))?;
            if residue_count > MAX_DIGESTED_RESIDUES {
                return Err(invalid(
                    "digested sequences exceed total residue limit; use digest_ranges or narrower length filters",
                ));
            }
        }
        ranges
            .into_iter()
            .map(|p| {
                Ok(DigestedPeptide {
                    sequence: protein.subsequence(p.range())?,
                    start: p.start,
                    end: p.end,
                    missed_cleavages: p.missed_cleavages,
                })
            })
            .collect()
    }
    /// Uppercase A-Z text accepts ambiguous residue codes without assigning them
    /// masses. No modifications, whitespace or separators are accepted here.
    pub fn digest_ranges(&self, sequence: &str) -> Result<Vec<DigestionProduct>> {
        let mut products = Vec::new();
        self.enumerate(sequence, |p| products.push(p))?;
        Ok(products)
    }
    pub fn digest_unmodified<'a>(&self, sequence: &'a str) -> Result<Vec<&'a str>> {
        Ok(self
            .digest_ranges(sequence)?
            .into_iter()
            .map(|p| &sequence[p.range()])
            .collect())
    }
    /// Counts products using the same specificity and length filters as digest.
    /// Does not allocate peptide strings or ranges, but observes work/product limits.
    pub fn peptide_count(&self, protein: &AASequence) -> Result<usize> {
        self.peptide_count_unmodified(protein.as_str())
    }
    pub fn peptide_count_unmodified(&self, sequence: &str) -> Result<usize> {
        self.enumerate(sequence, |_| {})
    }
    pub fn count_internal_cleavage_sites(&self, sequence: &str) -> Result<usize> {
        self.enzyme.count_internal_cleavage_sites(sequence)
    }
    /// Counts sites strictly within the supplied protein range, in full context.
    /// Empty/out-of-bounds/reversed ranges are invalid.
    pub fn count_missed_cleavages(&self, sequence: &str, range: Range<usize>) -> Result<usize> {
        validate_sequence(sequence)?;
        if range.start >= range.end || range.end > sequence.len() {
            return Err(invalid(
                "missed-cleavage range must be nonempty and within the sequence",
            ));
        }
        Ok(internal_count(
            &self.enzyme.cleavage_sites(sequence)?,
            range.start,
            range.end,
        ))
    }
    /// Check termini and optional missed cleavages; length filters are not part
    /// of validity, as in OpenMS. Invalid/empty ranges return false safely.
    pub fn is_valid_product(
        &self,
        protein: &AASequence,
        range: Range<usize>,
        options: ProductValidation,
    ) -> Result<bool> {
        self.is_valid_product_unmodified(protein.as_str(), range, options)
    }
    pub fn is_valid_product_unmodified(
        &self,
        sequence: &str,
        range: Range<usize>,
        options: ProductValidation,
    ) -> Result<bool> {
        self.validate()?;
        validate_sequence(sequence)?;
        if range.start >= range.end || range.end > sequence.len() {
            return Ok(false);
        }
        if self.enzyme == Protease::UnspecificCleavage {
            return Ok(true);
        }
        let mut start = range.start;
        let end = range.end;
        if options.allow_nterm_protein_cleavage && start <= 2 && sequence.starts_with('M') {
            start = 0;
        }
        let cuts = self.enzyme.cleavage_sites(sequence)?;
        let bytes = sequence.as_bytes();
        let random_site = |p: usize| {
            options.allow_random_asp_pro_cleavage
                && p > 0
                && p < bytes.len()
                && bytes[p - 1] == b'D'
                && bytes[p] == b'P'
        };
        let n_term = cuts.binary_search(&start).is_ok() || random_site(start);
        let c_term = cuts.binary_search(&end).is_ok() || random_site(end);
        let specific = match self.specificity {
            DigestionSpecificity::Full => n_term && c_term,
            DigestionSpecificity::Semi => n_term || c_term,
            DigestionSpecificity::None => true,
        };
        Ok(specific
            && (options.ignore_missed_cleavages
                || internal_count(&cuts, start, end) <= self.missed_cleavages))
    }

    fn enumerate(&self, sequence: &str, mut output: impl FnMut(DigestionProduct)) -> Result<usize> {
        self.validate()?;
        validate_sequence(sequence)?;
        let n = sequence.len();
        if n > self.max_work {
            return Err(invalid("digestion sequence scan exceeds work limit"));
        }
        if n < self.min_length {
            return Ok(0);
        }
        let max_length = self.max_length.unwrap_or(n).min(n);
        let cuts = self.enzyme.cleavage_sites(sequence)?;
        let mut work = n;
        let mut count = 0_usize;
        let mut emit = |start: usize, end: usize, missed_cleavages: usize| -> Result<()> {
            work = work
                .checked_add(1)
                .ok_or_else(|| invalid("digestion work overflows"))?;
            if work > self.max_work {
                return Err(invalid("digestion candidate work limit exceeded"));
            }
            let length = end - start;
            if length >= self.min_length && length <= max_length {
                if count >= self.max_products {
                    return Err(invalid("digestion product limit exceeded"));
                }
                count += 1;
                output(DigestionProduct {
                    start,
                    end,
                    missed_cleavages,
                });
            }
            Ok(())
        };
        if self.enzyme == Protease::UnspecificCleavage
            || self.specificity == DigestionSpecificity::None
        {
            // Mirror EnzymaticDigestion::digestUnmodified's unrestricted path.
            // The missed-cleavage ceiling is ignored for this enumeration.
            for start in 0..=n - self.min_length {
                for length in self.min_length..=max_length.min(n - start) {
                    let end = start + length;
                    emit(start, end, internal_count(&cuts, start, end))?;
                }
            }
            return Ok(count);
        }
        let segments = cuts.len() - 1;
        for missed in 0..=self.missed_cleavages.min(segments - 1) {
            for i in 0..segments - missed {
                emit(cuts[i], cuts[i + missed + 1], missed)?;
            }
        }
        if self.specificity == DigestionSpecificity::Semi {
            // At each offset, extend a nonspecific start forwards to successive
            // sites, then a nonspecific end backwards. Full products are excluded.
            let mut forward = 1;
            let mut backward = cuts.len() - 2;
            for shift in 1..n {
                if shift == cuts[forward] {
                    forward += 1;
                } else {
                    let available = cuts.len() - forward;
                    for missed in 0..=self.missed_cleavages.min(available - 1) {
                        emit(shift, cuts[forward + missed], missed)?;
                    }
                }
                let end = n - shift;
                if end == cuts[backward] {
                    backward = backward.saturating_sub(1);
                } else {
                    for missed in 0..=self.missed_cleavages.min(backward) {
                        emit(cuts[backward - missed], end, missed)?;
                    }
                }
            }
        }
        Ok(count)
    }
}

fn internal_count(cuts: &[usize], start: usize, end: usize) -> usize {
    cuts.partition_point(|&p| p < end) - cuts.partition_point(|&p| p <= start)
}
fn validate_sequence(sequence: &str) -> Result<()> {
    if sequence.len() > MAX_DIGESTION_SEQUENCE_LENGTH {
        return Err(invalid("sequence exceeds digestion length limit"));
    }
    if !sequence.bytes().all(|b| b.is_ascii_uppercase()) {
        return Err(invalid(
            "unmodified digestion sequence must contain only uppercase A-Z residues",
        ));
    }
    Ok(())
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
