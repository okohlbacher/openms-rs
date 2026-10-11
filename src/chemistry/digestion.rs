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
mod protease_db;
pub use enzymes::Protease;
pub use protease_db::ProteaseDB;

/// `EnzymaticDigestion::NoCleavage`, the name of the enzyme that never cleaves.
pub const NO_CLEAVAGE: &str = "no cleavage";
/// `EnzymaticDigestion::UnspecificCleavage`, the name of the enzyme that cuts
/// at every site.
pub const UNSPECIFIC_CLEAVAGE: &str = "unspecific cleavage";

/// Longest protein sequence, in residues, that digestion and cleavage-site
/// scans accept; the source has no ceiling.
pub const MAX_DIGESTION_SEQUENCE_LENGTH: usize = 1_000_000;
/// Ceiling on the summed length of the peptides [`ProteaseDigestion::digest`]
/// materialises in one call; the source has no ceiling.
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
    /// The compiled table entry this definition belongs to.
    pub fn protease(&self) -> Protease {
        self.protease
    }
    /// The enzyme name (`getName`).
    pub fn name(&self) -> &'static str {
        self.name
    }
    /// Source expression, exposed as metadata; it is not executed at runtime.
    pub fn regex(&self) -> &'static str {
        self.regex
    }
    /// Prose description of the cleavage rule (`getRegExDescription`).
    pub fn description(&self) -> &'static str {
        self.description
    }
    /// Alternative names the enzyme is registered under (`getSynonyms`),
    /// in the source set's sorted order.
    pub fn synonyms(&self) -> &'static [&'static str] {
        self.synonyms
    }
    /// Formula gained by the N-terminus of a cleavage product
    /// (`getNTermGain`); empty when the source declares none.
    pub fn n_term_gain(&self) -> EmpiricalFormula {
        EmpiricalFormula::parse(self.n_term_gain).expect("validated pinned enzyme formula")
    }
    /// Formula gained by the C-terminus of a cleavage product
    /// (`getCTermGain`); empty when the source declares none.
    pub fn c_term_gain(&self) -> EmpiricalFormula {
        EmpiricalFormula::parse(self.c_term_gain).expect("validated pinned enzyme formula")
    }
    /// PSI-MS controlled-vocabulary accession (`getPSIID`), empty when none.
    pub fn psi_id(&self) -> &'static str {
        self.psi_id
    }
    /// X! Tandem cleavage notation (`getXTandemID`), empty when none.
    pub fn xtandem_id(&self) -> &'static str {
        self.xtandem_id
    }
    /// Comet enzyme number (`getCometID`); the source's `-1` is `None`.
    pub fn comet_id(&self) -> Option<i32> {
        self.comet_id
    }
    /// MS-GF+ enzyme number (`getMSGFID`), counted from 0; the source's `-1`
    /// is `None`.
    pub fn msgf_id(&self) -> Option<i32> {
        self.msgf_id
    }
    /// OMSSA enzyme number (`getOMSSAID`); the source's `-1` is `None`.
    pub fn omssa_id(&self) -> Option<i32> {
        self.omssa_id
    }
}

impl Protease {
    /// Resolves a name, lower-cased name or synonym through
    /// [`ProteaseDB::get_enzyme`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a name no enzyme is registered under.
    pub fn from_name(name: &str) -> Result<Self> {
        Ok(ProteaseDB::global().get_enzyme(name)?.protease)
    }
    /// The full pinned definition of this enzyme.
    pub fn metadata(self) -> &'static DigestionEnzymeProtein {
        enzymes::ENZYMES
            .iter()
            .find(|e| e.protease == self)
            .expect("every protease has a pinned definition")
    }
    /// The enzyme name, as [`DigestionEnzymeProtein::name`].
    pub fn name(self) -> &'static str {
        self.metadata().name
    }
    /// All internal cleavage positions plus protein boundaries 0 and len, without
    /// duplicate boundaries. Empty input returns `[0]`. Accepts uppercase A-Z.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for text that is not uppercase A-Z or exceeds
    /// [`MAX_DIGESTION_SEQUENCE_LENGTH`].
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
    ///
    /// # Errors
    ///
    /// As [`Protease::cleavage_sites`].
    pub fn is_cleavage_site(self, sequence: &str, position: usize) -> Result<bool> {
        validate_sequence(sequence)?;
        Ok(position <= sequence.len()
            && (position == 0
                || position == sequence.len()
                || self.metadata().rule.matches(sequence.as_bytes(), position)))
    }
    /// Number of cleavage sites strictly inside `sequence`, protein termini
    /// excluded.
    ///
    /// # Errors
    ///
    /// As [`Protease::cleavage_sites`].
    pub fn count_internal_cleavage_sites(self, sequence: &str) -> Result<usize> {
        Ok(self.cleavage_sites(sequence)?.len().saturating_sub(2))
    }
    /// The source's `EnzymaticDigestion::countInternalCleavageSites` as
    /// `filterByMissedCleavages` sees it: the number of tokens its
    /// Boost `sregex_token_iterator` split yields, minus one, cast to `Int`.
    ///
    /// It differs from [`Self::count_internal_cleavage_sites`] in one case: a
    /// zero-width match at position 0 (a lookahead-only rule such as Asp-N or
    /// Lys-N on a peptide that starts with D or K) produces an empty leading
    /// token and counts as one more site. A match at the end adds no trailing
    /// token, and `no cleavage` is not split at all (one token, count 0). The
    /// executed C++ oracle behind `tests/id_filter_oracle.rs` checks this.
    ///
    /// # Errors
    ///
    /// As [`Protease::cleavage_sites`].
    pub(crate) fn source_token_cleavage_count(self, sequence: &str) -> Result<i64> {
        let internal = self.count_internal_cleavage_sites(sequence)?;
        let leading = self.matches_with_zero_width_at_start(sequence);
        Ok(i64::try_from(internal).unwrap_or(i64::MAX) + i64::from(leading))
    }

    /// Whether the enzyme's expression matches with zero width before the
    /// first residue: a lookahead-only rule (Asp-N `(?=[DBX])`, Lys-N
    /// `(?=[KX])`) on a sequence that starts with one of its residues.
    fn matches_with_zero_width_at_start(self, sequence: &str) -> bool {
        match self.metadata().rule {
            CleavageRule::Boundary {
                previous,
                previous_two,
                next,
                forbidden_next,
            } => {
                previous.is_empty()
                    && previous_two.is_empty()
                    && sequence.as_bytes().first().is_some_and(|first| {
                        (next.is_empty() || next.contains(first)) && !forbidden_next.contains(first)
                    })
            }
            CleavageRule::Never => false,
        }
    }

    /// `EnzymaticDigestion::tokenize_(sequence)` at `bc9cc12`
    /// (`EnzymaticDigestion.cpp:92-116`): the start of every token Boost's
    /// `sregex_token_iterator` (submatch -1) splits the sequence into, without
    /// the end position.
    ///
    /// It differs from [`Self::cleavage_sites`] in the source's own ways: a
    /// zero-width match at the first residue (see
    /// [`Self::matches_with_zero_width_at_start`]) yields position 0 twice, a
    /// match at the end adds no token, the empty sequence yields no token,
    /// except for the `no cleavage` entry, whose expression `()` is not split
    /// at all and yields `[0]`.
    ///
    /// # Errors
    ///
    /// As [`Self::cleavage_sites`].
    pub(crate) fn source_tokenize(self, sequence: &str) -> Result<Vec<usize>> {
        let cuts = self.cleavage_sites(sequence)?;
        if self == Self::NoCleavage {
            return Ok(vec![0]);
        }
        if sequence.is_empty() {
            return Ok(Vec::new());
        }
        // `cuts` is `[0, internal sites..., len]` for a nonempty sequence.
        let mut positions = Vec::with_capacity(cuts.len());
        if self.matches_with_zero_width_at_start(sequence) {
            positions.push(0);
        }
        positions.extend_from_slice(&cuts[..cuts.len() - 1]);
        Ok(positions)
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

/// How many termini of a product must be enzymatic cleavage sites
/// (`EnzymaticDigestion::Specificity`).
///
/// [`Full`](Self::Full), [`Semi`](Self::Semi) and [`None`](Self::None) are the
/// values a digestion enumerates with. The source enumeration has three more,
/// which a digestion object can hold and which behave as follows there:
/// `SPEC_UNKNOWN`, the value [`DigestionSpecificity::by_name`] answers for a
/// name it does not know, and `SPEC_NOCTERM` and `SPEC_NONTERM`, "currently not
/// supported in the class". `ProteaseDigestion::digest` refuses all three
/// with `Exception::InvalidValue`; `isValidProduct` and `digestUnmodified`
/// treat each as full specificity, as here.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DigestionSpecificity {
    /// Both termini must be protein boundaries or cleavage sites
    /// (`SPEC_FULL`, 2).
    #[default]
    Full,
    /// At least one terminus must be a protein boundary or cleavage site
    /// (`SPEC_SEMI`, 1).
    Semi,
    /// Enumerate every substring in the length interval, without enzyme limits
    /// (`SPEC_NONE`, 0).
    None,
    /// `SPEC_UNKNOWN` (3).
    Unknown,
    /// `SPEC_NOCTERM` (8): no requirement on the C-terminus, "currently not
    /// supported in the class".
    NoCTerm,
    /// `SPEC_NONTERM` (9): no requirement on the N-terminus, "currently not
    /// supported in the class".
    NoNTerm,
}
impl DigestionSpecificity {
    /// `EnzymaticDigestion::NamesOfSpecificity`, indexed by the source's
    /// enumeration value: the entries 3 to 7 of the source's table are all
    /// `unknown`, and only 3 is an enumeration value.
    pub const NAMES_OF_SPECIFICITY: [&'static str; 10] = [
        "none", "semi", "full", "unknown", "unknown", "unknown", "unknown", "unknown", "no-cterm",
        "no-nterm",
    ];

    /// The source's parameter spelling: `full`, `semi`, `none`, `unknown`,
    /// `no-cterm` or `no-nterm`.
    pub fn name(self) -> &'static str {
        Self::NAMES_OF_SPECIFICITY[usize::from(self.value())]
    }

    /// The source's enumeration value: `SPEC_NONE` 0, `SPEC_SEMI` 1,
    /// `SPEC_FULL` 2, `SPEC_UNKNOWN` 3, `SPEC_NOCTERM` 8, `SPEC_NONTERM` 9.
    /// Some engines report just the number of required termini, which is why
    /// the first three are what they are.
    pub fn value(self) -> u8 {
        match self {
            Self::None => 0,
            Self::Semi => 1,
            Self::Full => 2,
            Self::Unknown => 3,
            Self::NoCTerm => 8,
            Self::NoNTerm => 9,
        }
    }

    /// The member of a source enumeration value, or `None` for a number the
    /// enumeration has no member for (4 to 7, and `SIZE_OF_SPECIFICITY`).
    pub fn from_value(value: u8) -> Option<Self> {
        match value {
            0 => Some(Self::None),
            1 => Some(Self::Semi),
            2 => Some(Self::Full),
            3 => Some(Self::Unknown),
            8 => Some(Self::NoCTerm),
            9 => Some(Self::NoNTerm),
            _ => Option::None,
        }
    }

    /// `EnzymaticDigestion::getSpecificityByName`: the first entry of
    /// [`NAMES_OF_SPECIFICITY`](Self::NAMES_OF_SPECIFICITY) equal to `name`, and
    /// [`Unknown`](Self::Unknown) when there is none. Unlike [`FromStr`] it
    /// never refuses a name.
    pub fn by_name(name: &str) -> Self {
        Self::NAMES_OF_SPECIFICITY
            .iter()
            .position(|candidate| *candidate == name)
            .and_then(|index| u8::try_from(index).ok())
            .and_then(Self::from_value)
            .unwrap_or(Self::Unknown)
    }

    /// Whether `ProteaseDigestion::digest` enumerates with this value: full
    /// and semi.
    fn digest_supports(self) -> bool {
        matches!(self, Self::Full | Self::Semi)
    }
}
impl FromStr for DigestionSpecificity {
    type Err = Error;
    /// The three names a digestion is configured with: `full`, `semi` and
    /// `none`. [`DigestionSpecificity::by_name`] is the source's lookup, which
    /// answers `unknown` for any other name instead of refusing it.
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
    /// Accept a product regardless of its missed cleavages (the source's
    /// `ignore_missed_cleavages`, default `true`).
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
    /// The digesting enzyme (default trypsin).
    pub enzyme: Protease,
    /// Maximum number of missed cleavages in a product (default 0).
    pub missed_cleavages: usize,
    /// Shortest product, in residues (default 1; must be at least 1).
    pub min_length: usize,
    /// Longest product, in residues; `None` is unbounded.
    pub max_length: Option<usize>,
    /// Required terminal specificity (default full).
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

/// One digestion product as a residue range of the input protein.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DigestionProduct {
    /// Zero-based inclusive start.
    pub start: usize,
    /// Zero-based exclusive end.
    pub end: usize,
    /// Actual internal enzyme sites, including in unrestricted products.
    pub missed_cleavages: usize,
}
impl DigestionProduct {
    /// The product as a `start..end` range.
    pub fn range(self) -> Range<usize> {
        self.start..self.end
    }
    /// Product length in residues.
    pub fn len(self) -> usize {
        self.end.saturating_sub(self.start)
    }
    /// Whether the range holds no residue.
    pub fn is_empty(self) -> bool {
        self.start >= self.end
    }
}

/// One digestion product with its peptide sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DigestedPeptide {
    /// The product sequence, with the modifications of the covered residues.
    pub sequence: AASequence,
    /// Zero-based inclusive start in the input protein.
    pub start: usize,
    /// Zero-based exclusive end in the input protein.
    pub end: usize,
    /// Internal enzyme sites inside the product.
    pub missed_cleavages: usize,
}

impl ProteaseDigestion {
    /// Selects the enzyme by name, lower-cased name or synonym (`setEnzyme`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for an unknown name; the setting is unchanged.
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
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for invalid settings, a sequence outside the
    /// accepted alphabet or length, or an exceeded product, work or
    /// [`MAX_DIGESTED_RESIDUES`] ceiling.
    pub fn digest(&self, protein: &AASequence) -> Result<Vec<DigestedPeptide>> {
        if !self.specificity.digest_supports() && self.specificity != DigestionSpecificity::None {
            // ProteaseDigestion.cpp:93-98 refuses what is neither full nor
            // semi; the native digest also accepts none, see the module notes.
            return Err(invalid(format!(
                "Specificity value set on current ProteaseDigestion object is not supported by ProteaseDigestion::digest(). Value: {}",
                self.specificity.value()
            )));
        }
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
    ///
    /// # Errors
    ///
    /// As [`Self::digest`], without the residue ceiling.
    pub fn digest_ranges(&self, sequence: &str) -> Result<Vec<DigestionProduct>> {
        let mut products = Vec::new();
        self.enumerate(sequence, |p| products.push(p))?;
        Ok(products)
    }
    /// Digests unmodified one-letter text into borrowed product slices.
    ///
    /// # Errors
    ///
    /// As [`Self::digest_ranges`].
    pub fn digest_unmodified<'a>(&self, sequence: &'a str) -> Result<Vec<&'a str>> {
        Ok(self
            .digest_ranges(sequence)?
            .into_iter()
            .map(|p| &sequence[p.range()])
            .collect())
    }
    /// Counts products using the same specificity and length filters as digest.
    /// Does not allocate peptide strings or ranges, but observes work/product limits.
    ///
    /// # Errors
    ///
    /// As [`Self::digest_ranges`].
    pub fn peptide_count(&self, protein: &AASequence) -> Result<usize> {
        self.peptide_count_unmodified(protein.as_str())
    }
    /// [`Self::peptide_count`] over unmodified one-letter text.
    ///
    /// # Errors
    ///
    /// As [`Self::digest_ranges`].
    pub fn peptide_count_unmodified(&self, sequence: &str) -> Result<usize> {
        self.enumerate(sequence, |_| {})
    }
    /// Number of internal cleavage sites of the configured enzyme.
    ///
    /// # Errors
    ///
    /// As [`Protease::cleavage_sites`].
    pub fn count_internal_cleavage_sites(&self, sequence: &str) -> Result<usize> {
        self.enzyme.count_internal_cleavage_sites(sequence)
    }
    /// Counts sites strictly within the supplied protein range, in full context.
    /// Empty/out-of-bounds/reversed ranges are invalid.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for such a range or an invalid sequence.
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
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for invalid settings or sequence text.
    pub fn is_valid_product(
        &self,
        protein: &AASequence,
        range: Range<usize>,
        options: ProductValidation,
    ) -> Result<bool> {
        self.is_valid_product_unmodified(protein.as_str(), range, options)
    }
    /// [`Self::is_valid_product`] over unmodified one-letter text.
    ///
    /// # Errors
    ///
    /// As [`Self::is_valid_product`].
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
        // The source treats every value but none and semi as full.
        let specific = match self.specificity {
            DigestionSpecificity::Semi => n_term || c_term,
            DigestionSpecificity::None => true,
            DigestionSpecificity::Full
            | DigestionSpecificity::Unknown
            | DigestionSpecificity::NoCTerm
            | DigestionSpecificity::NoNTerm => n_term && c_term,
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

/// What the source's `ProteaseDigestion::digest` (pair overload) returns: the
/// products as `(start, past-the-end)` pairs and its count of discarded
/// candidates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SourceDigestion {
    /// The products, in the source's output order.
    pub ranges: Vec<(usize, usize)>,
    /// The source's return value `wrong_size`: candidates the length filter
    /// rejected, including the empty first token the source's tokenisation
    /// yields when a lookahead-only expression matches at the first residue.
    pub discarded: usize,
}

/// What the source's `isValidProduct_` answers, with the warnings it writes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SourceValidity {
    /// The verdict.
    pub valid: bool,
    /// The lines the source writes to `OPENMS_LOG_WARN` (without their line
    /// ends): a start or an end beyond the sequence, or an empty fragment or
    /// sequence.
    pub warnings: Vec<String>,
}

/// The source's own tokenisation and loops, selected by calling these methods
/// instead of [`ProteaseDigestion::digest_ranges`] and
/// [`ProteaseDigestion::peptide_count`].
///
/// The native operations (see `docs/DIGESTION_SUPPORT.md`) keep each cleavage
/// boundary once and order and count by the settings; the source's
/// `ProteaseDigestion::digest` and `peptideCount` split the sequence with a
/// Boost token iterator and loop over its positions as they come, with
/// unsigned arithmetic. The methods below port those loops literally, so every
/// result is the Release build's, including its defects on a sequence that a
/// lookahead-only expression (Asp-N, Asp-N_ambic, Lys-N) matches at the first
/// residue, and its ordering of an unspecific digest.
impl ProteaseDigestion {
    /// `ProteaseDigestion::digest(const AASequence&, vector<pair<size_t,size_t>>&,
    /// min_length, max_length)` (`ProteaseDigestion.cpp:87-158`), with
    /// [`Self::min_length`] and [`Self::max_length`] (`None` is the source's
    /// 0, unbounded) as the length filter.
    ///
    /// The tokenisation is `Protease::source_tokenize`'s, so a sequence an
    /// Asp-N, Asp-N_ambic or Lys-N expression matches at its first residue
    /// starts with two token positions 0: the first token is empty and is
    /// counted in [`SourceDigestion::discarded`], the first product is repeated
    /// as a product with a missed cleavage, and the semi-specific scan
    /// (`EnzymaticDigestion::semiSpecificDigestion_`, with the source's
    /// `Size` wrap of a negative length) never advances its forward cursor and
    /// loses products. For `unspecific cleavage` the missed-cleavage ceiling
    /// is the source's `numeric_limits<Size>::max()`, so the products come
    /// ordered by length, then position.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for [`DigestionSpecificity::None`], which the
    /// source's `digest` refuses with `Exception::InvalidValue`, for invalid
    /// settings or a sequence outside the accepted alphabet or length, and for
    /// an exceeded product or work ceiling.
    pub fn digest_ranges_source(&self, sequence: &str) -> Result<SourceDigestion> {
        self.validate()?;
        if !self.specificity.digest_supports() {
            return Err(invalid(format!(
                "Specificity value set on current ProteaseDigestion object is not supported by ProteaseDigestion::digest(). Value: {}",
                self.specificity.value()
            )));
        }
        validate_sequence(sequence)?;
        let size = sequence.len();
        let max_length = self.max_length.map_or(size, |max| max.min(size));
        let min_length = self.min_length;
        let mc = if self.enzyme == Protease::UnspecificCleavage {
            usize::MAX
        } else {
            self.missed_cleavages
        };

        let mut positions = self.enzyme.source_tokenize(sequence)?;
        positions.push(size);
        let count = positions.len();

        let mut ranges: Vec<(usize, usize)> = Vec::new();
        let mut discarded = 0_usize;
        let mut work = size;
        if work > self.max_work {
            return Err(invalid("digestion sequence scan exceeds work limit"));
        }
        // Source `variant` / the inline checks: `Size` arithmetic wraps.
        let mut candidate = |start: usize, end: usize| -> Result<()> {
            work = work
                .checked_add(1)
                .ok_or_else(|| invalid("digestion work overflows"))?;
            if work > self.max_work {
                return Err(invalid("digestion candidate work limit exceeded"));
            }
            let length = end.wrapping_sub(start);
            if length >= min_length && length <= max_length {
                if ranges.len() >= self.max_products {
                    return Err(invalid("digestion product limit exceeded"));
                }
                ranges.push((start, end));
            } else {
                discarded += 1;
            }
            Ok(())
        };

        let mut begin = positions[0];
        for &position in &positions[1..] {
            candidate(begin, position)?;
            begin = position;
        }

        // missed cleavages
        if count > 1 && mc != 0 {
            let mut mcs = 1_usize;
            while mcs <= mc && mcs < count - 1 {
                begin = positions[0];
                for j in 1..count - mcs {
                    candidate(begin, positions[j + mcs])?;
                    begin = positions[j];
                }
                mcs += 1;
            }
        }

        // semi-specific variants
        if self.specificity == DigestionSpecificity::Semi {
            // Source: fewer than two positions throws InvalidValue.
            if count < 2 {
                return Err(invalid(
                    "Too few cleavage positions - at least sequence start and end positions are required.",
                ));
            }
            let mc = i64::try_from(self.missed_cleavages).unwrap_or(i64::MAX);
            let first = positions[0];
            let last = positions[count - 1];
            let length = last - first;
            let mut forward = 1_usize;
            // The reverse iterator, as an offset from the end: position
            // `count - 1 - backward`.
            let mut backward = 1_usize;
            for shift in 1..length {
                if first + shift != positions[forward] {
                    let mut i = 0_i64;
                    while i < mc.saturating_add(1) {
                        let Some(&end) = usize::try_from(i)
                            .ok()
                            .and_then(|i| positions.get(forward + i))
                        else {
                            break;
                        };
                        candidate(first + shift, end)?;
                        if end == last {
                            break;
                        }
                        i += 1;
                    }
                } else {
                    forward += 1;
                }
                if last - shift != positions[count - 1 - backward] {
                    let mut i = 0_i64;
                    while i < mc.saturating_add(1) {
                        let Some(&start) = usize::try_from(i)
                            .ok()
                            .and_then(|i| (backward + i).checked_add(1))
                            .and_then(|back| count.checked_sub(back))
                            .and_then(|index| positions.get(index))
                        else {
                            break;
                        };
                        candidate(start, last - shift)?;
                        if start == first {
                            break;
                        }
                        i += 1;
                    }
                } else {
                    backward += 1;
                }
            }
        }
        Ok(SourceDigestion { ranges, discarded })
    }

    /// `EnzymaticDigestion::digestUnmodified`, pair overload
    /// (`EnzymaticDigestion.cpp:505-564`), as `(start, past-the-end)` pairs;
    /// the source writes `(start, length)`.
    ///
    /// It is [`Self::digest_ranges_source`] but for three things the source
    /// does before it tokenises: a sequence shorter than
    /// [`Self::min_length`] yields nothing and counts nothing, `unspecific
    /// cleavage` and specificity `none` enumerate every substring ordered by
    /// start and then length with nothing discarded, and a specificity that is
    /// neither (`SPEC_UNKNOWN`, `SPEC_NOCTERM`, `SPEC_NONTERM`) is digested as
    /// full instead of being refused.
    ///
    /// # Errors
    ///
    /// As [`Self::digest_ranges_source`], without the refusal of a
    /// specificity.
    pub fn digest_unmodified_source(&self, sequence: &str) -> Result<SourceDigestion> {
        self.validate()?;
        validate_sequence(sequence)?;
        if sequence.len() < self.min_length {
            return Ok(SourceDigestion {
                ranges: Vec::new(),
                discarded: 0,
            });
        }
        if self.enzyme == Protease::UnspecificCleavage
            || self.specificity == DigestionSpecificity::None
        {
            let ranges = self
                .digest_ranges(sequence)?
                .into_iter()
                .map(|product| (product.start, product.end))
                .collect();
            return Ok(SourceDigestion {
                ranges,
                discarded: 0,
            });
        }
        if self.specificity.digest_supports() {
            return self.digest_ranges_source(sequence);
        }
        Self {
            specificity: DigestionSpecificity::Full,
            ..self.clone()
        }
        .digest_ranges_source(sequence)
    }

    /// `ProteaseDigestion::peptideCount` (`ProteaseDigestion.cpp:48-69`): the
    /// number of tokens plus those of every missed-cleavage level up to
    /// [`Self::missed_cleavages`], counted from the source's tokenisation
    /// without the length filter and without the semi-specific products;
    /// `unspecific cleavage` is `n (n + 1) / 2`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a sequence outside the accepted alphabet or
    /// length.
    pub fn peptide_count_source(&self, sequence: &str) -> Result<usize> {
        validate_sequence(sequence)?;
        let size = sequence.len();
        if self.enzyme == Protease::UnspecificCleavage {
            // Bounded by MAX_DIGESTION_SEQUENCE_LENGTH, so no overflow.
            return Ok((size + 1) * size / 2);
        }
        let count = self.enzyme.source_tokenize(sequence)?.len();
        let mut sum = count;
        for i in 1..count {
            if i > self.missed_cleavages {
                break;
            }
            sum += count - i;
        }
        Ok(sum)
    }

    /// `EnzymaticDigestion::countInternalCleavageSites`
    /// (`EnzymaticDigestion.cpp:195-198`): `tokenize_(sequence).size() - 1`
    /// in `Size` arithmetic, which is one above
    /// [`Self::count_internal_cleavage_sites`] for a sequence a lookahead-only
    /// expression matches at its first residue and wraps to `usize::MAX` for
    /// the empty sequence (no token) of every enzyme but `no cleavage`.
    ///
    /// # Errors
    ///
    /// As [`Self::peptide_count_source`].
    pub fn count_internal_cleavage_sites_source(&self, sequence: &str) -> Result<usize> {
        validate_sequence(sequence)?;
        Ok(self.enzyme.source_tokenize(sequence)?.len().wrapping_sub(1))
    }

    /// `EnzymaticDigestion::isValidProduct_` (`EnzymaticDigestion.cpp:211-303`),
    /// the routine behind both `isValidProduct` overloads of
    /// `ProteaseDigestion` and `EnzymaticDigestion::isValidProduct`, with the
    /// source's `int` position and length.
    ///
    /// The native [`Self::is_valid_product`] takes a range and answers `false`
    /// for an empty or out-of-bounds one before anything else; the source
    /// applies the N-terminal methionine allowance first (a position of 1 or
    /// 2 with length 0 becomes `0..position`, which is valid), accepts a
    /// negative length, tokenises only the fragment's neighbourhood
    /// (`tokenize_(sequence, pos - 2, end + 2)`, so a lookbehind cannot see
    /// the residue before the window), and writes three warnings. They are
    /// returned in [`SourceValidity::warnings`] instead of being logged.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for invalid settings, a sequence outside the
    /// accepted alphabet or length, an `int` overflow of the source's
    /// arithmetic, and where the source reads `sequence[-1]` or a residue
    /// before the start for the random D|P allowance (undefined there).
    pub fn is_valid_product_source(
        &self,
        sequence: &str,
        position: i32,
        length: i32,
        options: ProductValidation,
    ) -> Result<SourceValidity> {
        validate_sequence(sequence)?;
        let overflow = || invalid("fragment position and length overflow the source's int");
        let bytes = sequence.as_bytes();
        let size = i32::try_from(bytes.len()).map_err(|_| overflow())?;
        let (mut pos, mut length) = (position, length);
        let mut warnings = Vec::new();
        let refuse = |valid: bool, warnings: Vec<String>| Ok(SourceValidity { valid, warnings });
        // `sequence[0] == 'M'`: a std::string's terminator reads '\0'.
        if options.allow_nterm_protein_cleavage && pos <= 2 && bytes.first() == Some(&b'M') {
            length = length.checked_add(pos).ok_or_else(overflow)?;
            pos = 0;
        }
        if pos >= size {
            warnings.push(format!(
                "Warning: start of fragment ({pos}) is beyond end of sequence '{sequence}'!"
            ));
            return refuse(false, warnings);
        }
        let end = pos.checked_add(length).ok_or_else(overflow)?;
        if end > size {
            warnings.push(format!(
                "Warning: end of fragment ({end}) is beyond end of sequence '{sequence}'!"
            ));
            return refuse(false, warnings);
        }
        if length == 0 || sequence.is_empty() {
            warnings.push("Warning: fragment and sequence must not be empty!".to_owned());
            return refuse(false, warnings);
        }
        if self.enzyme == Protease::UnspecificCleavage {
            return refuse(true, warnings);
        }
        let missed = self.missed_cleavages;
        if self.specificity == DigestionSpecificity::None {
            if options.ignore_missed_cleavages {
                return refuse(true, warnings);
            }
            if end >= 0 && end < pos {
                // A negative length: the source builds an iterator range whose
                // end precedes its start (undefined; the Release build
                // crashes). A negative end is reset to the size instead.
                return Err(invalid(
                    "the source tokenises a window whose end precedes its start here (undefined)",
                ));
            }
            let tokens = self.source_tokenize_window(sequence, pos, end)?;
            return refuse(tokens.len().wrapping_sub(1) <= missed, warnings);
        }
        if self.specificity == DigestionSpecificity::Full
            && self.enzyme == Protease::NoCleavage
            && !options.allow_random_asp_pro_cleavage
        {
            return refuse(pos == 0 && size == end, warnings);
        }
        let positions = self.source_tokenize_window(
            sequence,
            pos.checked_sub(2).ok_or_else(overflow)?,
            end.checked_add(2).ok_or_else(overflow)?,
        )?;
        let at = |index: i32| -> Result<u8> {
            usize::try_from(index)
                .ok()
                .and_then(|index| bytes.get(index))
                .copied()
                .ok_or_else(|| {
                    invalid("the source reads a residue before the start of the sequence here (undefined)")
                })
        };
        let pos_usize = usize::try_from(pos).ok();
        let spec_n = if pos_usize.is_some_and(|p| positions.contains(&p)) {
            true
        } else {
            options.allow_random_asp_pro_cleavage
                && at(pos.checked_sub(1).ok_or_else(overflow)?)? == b'D'
                && at(pos)? == b'P'
        };
        let spec_c =
            if end == size || usize::try_from(end).is_ok_and(|end| positions.contains(&end)) {
                true
            } else {
                options.allow_random_asp_pro_cleavage
                    && at(end.checked_sub(1).ok_or_else(overflow)?)? == b'D'
                    && at(end)? == b'P'
            };
        let specific = (spec_n && spec_c)
            || (self.specificity == DigestionSpecificity::Semi && (spec_n || spec_c));
        if !specific {
            return refuse(false, warnings);
        }
        if options.ignore_missed_cleavages {
            return refuse(true, warnings);
        }
        // countMissedCleavages_: the positions strictly inside the fragment.
        let inside = positions
            .iter()
            .filter(|&&p| i64::try_from(p).is_ok_and(|p| i64::from(pos) < p && p < i64::from(end)))
            .count();
        refuse(inside <= missed, warnings)
    }

    /// `EnzymaticDigestion::tokenize_(sequence, start, end)`
    /// (`EnzymaticDigestion.cpp:92-116`): the tokens of the window
    /// `[max(0, start), end)`, `end` being the size when it is negative or
    /// beyond it. The expression sees only the window, so a lookbehind does
    /// not see the residue before it and a lookahead does not see the one
    /// after it.
    fn source_tokenize_window(&self, sequence: &str, start: i32, end: i32) -> Result<Vec<usize>> {
        let size = sequence.len();
        let start = usize::try_from(start.max(0)).unwrap_or(0).min(size);
        let end = usize::try_from(end)
            .ok()
            .filter(|&end| end <= size)
            .unwrap_or(size);
        let window = sequence.get(start..end.max(start)).unwrap_or("");
        let tokens = self.enzyme.source_tokenize(window)?;
        Ok(tokens.into_iter().map(|offset| start + offset).collect())
    }

    /// `EnzymaticDigestion::filterByMissedCleavages`
    /// (`EnzymaticDigestion.cpp:205-209`): `filter` receives
    /// [`Self::count_internal_cleavage_sites_source`] as the source's `Int`
    /// (`const int mc = countInternalCleavageSites(sequence)`, which wraps the
    /// empty sequence's `Size` result to `-1`) and says whether the sequence
    /// is to be filtered out.
    ///
    /// # Errors
    ///
    /// As [`Self::count_internal_cleavage_sites_source`].
    pub fn filter_by_missed_cleavages(
        &self,
        sequence: &str,
        filter: impl FnOnce(i32) -> bool,
    ) -> Result<bool> {
        let count = self.count_internal_cleavage_sites_source(sequence)?;
        // Source: `Size` to `int`, modulo 2^32 on the Linux x86_64 build.
        Ok(filter(count as u32 as i32))
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
