// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/ID/AhoCorasickAmbiguous.h`: an Aho-Corasick trie over
//! amino-acid needles that also finds matches through ambiguous residues
//! (`B`, `J`, `Z`, `X`) and a bounded number of mismatches.
//!
//! Needles (peptides) are added to an [`ACTrie`], which [`ACTrie::compress_trie`]
//! then lays out in breadth-first order with suffix links. A query (protein)
//! is walked with an [`ACTrieState`]: the *primary* follows the query as an
//! ordinary Aho-Corasick automaton, and at every ambiguous residue or
//! potential mismatch it spawns *scouts* ([`ACScout`]) that follow each
//! alternative residue until they lose the residue that spawned them.
//!
//! The eight public types of the header are here: [`AA`] (with the
//! [`AA_TO_CHAR`] and [`CHAR_TO_AA`] tables), [`Hit`], [`Index`], [`Bitset`],
//! [`ACNode`] (with [`DepthHits`]), [`ACScout`], [`ACTrieState`] and [`ACTrie`],
//! plus the free function [`next_valid_aa`]. The mapping, the source defects
//! and the native bounds are in `docs/AHO_CORASICK_AMBIGUOUS_SUPPORT.md`.
//!
//! [`ACTrie`]: crate::analysis::id::aho_corasick_ambiguous::ACTrie
//! [`ACTrie::compress_trie`]: crate::analysis::id::aho_corasick_ambiguous::ACTrie::compress_trie
//! [`ACTrieState`]: crate::analysis::id::aho_corasick_ambiguous::ACTrieState
//! [`ACScout`]: crate::analysis::id::aho_corasick_ambiguous::ACScout
//! [`AA`]: crate::analysis::id::aho_corasick_ambiguous::AA
//! [`AA_TO_CHAR`]: crate::analysis::id::aho_corasick_ambiguous::AA_TO_CHAR
//! [`CHAR_TO_AA`]: crate::analysis::id::aho_corasick_ambiguous::CHAR_TO_AA
//! [`Hit`]: crate::analysis::id::aho_corasick_ambiguous::Hit
//! [`Index`]: crate::analysis::id::aho_corasick_ambiguous::Index
//! [`Bitset`]: crate::analysis::id::aho_corasick_ambiguous::Bitset
//! [`ACNode`]: crate::analysis::id::aho_corasick_ambiguous::ACNode
//! [`DepthHits`]: crate::analysis::id::aho_corasick_ambiguous::DepthHits
//! [`next_valid_aa`]: crate::analysis::id::aho_corasick_ambiguous::next_valid_aa

use crate::{Error, Result};
use std::collections::VecDeque;
use std::ops::{Shl, ShlAssign, Sub};

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

/// Internal amino-acid code to one-letter symbol (source `AAtoChar`).
///
/// Codes 0 to 21 are the unambiguous residues `A` to `V` in the source's order;
/// the ambiguous `B`, `J`, `Z`, `X` are consecutive at 22 to 25, followed by
/// `$` (26, the "super AA" that models a mismatch) and `?` (27, invalid, which
/// must stay last). `D`/`N`, `I`/`L` and `E`/`Q` are neighbours, so the
/// residues an ambiguous code stands for form one range.
pub const AA_TO_CHAR: [u8; 28] = [
    b'A', b'Y', b'C', b'D', b'N', b'F', b'G', b'H', b'I', b'L', b'K', b'W', b'M', b'O', b'P', b'E',
    b'Q', b'R', b'S', b'T', b'U', b'V', b'B', b'J', b'Z', b'X', b'$', b'?',
];

/// Byte to internal amino-acid code (source `CharToAA`), for all 256 byte
/// values.
///
/// Letters map case-insensitively, `$` maps to 26 and every other byte,
/// including bytes 128 to 255, maps to 27, the invalid code.
pub const CHAR_TO_AA: [u8; 256] = {
    let mut table = [27u8; 256];
    let mut code = 0;
    while code < 26 {
        let letter = AA_TO_CHAR[code];
        table[letter as usize] = code as u8;
        table[letter.to_ascii_lowercase() as usize] = code as u8;
        code += 1;
    }
    table[b'$' as usize] = 26;
    table
};

/// One amino acid in the trie's internal code (source struct `AA`).
///
/// Ordered by code, so the ambiguous codes compare above every unambiguous
/// one; the source defines only `==` and `<=`, which this order agrees with.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AA(u8);

impl Default for AA {
    /// The invalid amino acid `?` (source default constructor).
    fn default() -> Self {
        Self::INVALID
    }
}

impl Sub for AA {
    type Output = AA;
    /// Code difference (source `operator-`), wrapping as the source's
    /// `uint8_t` arithmetic does.
    fn sub(self, rhs: AA) -> AA {
        AA(self.0.wrapping_sub(rhs.0))
    }
}

impl AA {
    /// The invalid amino acid, `?` (code 27).
    pub const INVALID: AA = AA(27);
    /// The super amino acid `$` (code 26), which models a mismatch.
    pub const SUPER: AA = AA(26);

    /// The amino acid for a byte (source `AA(char)`).
    ///
    /// Any letter `A`-`Z` or `a`-`z` gives a valid amino acid; `$` gives the
    /// super amino acid, which should only model mismatches; every other byte
    /// gives the invalid `?`.
    pub const fn from_byte(c: u8) -> Self {
        AA(CHAR_TO_AA[c as usize])
    }

    /// The internal 8-bit code (source `operator()`).
    pub const fn code(self) -> u8 {
        self.0
    }

    /// Whether this is `B`, `J`, `Z`, `X` or `$`.
    pub const fn is_ambiguous(self) -> bool {
        self.0 >= 22
    }

    /// Whether this is a letter or `$`, i.e. not `?`.
    pub const fn is_valid(self) -> bool {
        self.0 != 27
    }

    /// Whether this is a letter (`A`-`Z`, `a`-`z`), i.e. neither `$` nor `?`.
    pub const fn is_valid_for_peptide(self) -> bool {
        self.0 <= 25
    }

    /// The next amino acid in code order (source `operator++`), or `None`
    /// past `?`. The source only asserts in debug builds that it does not step
    /// past `?`, after which its table read is out of bounds.
    pub const fn succ(self) -> Option<AA> {
        if self.0 < 27 {
            Some(AA(self.0 + 1))
        } else {
            None
        }
    }

    /// The one-letter symbol: `A`-`V`, `B`, `J`, `Z`, `X`, `$`, or `?` for the
    /// invalid amino acid.
    pub const fn to_char(self) -> char {
        AA_TO_CHAR[self.0 as usize] as char
    }

    /// The unambiguous amino acid with code `index` (0 to 21, `A` to `V`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `index` is not below
    /// [`AA::unambiguous_aa_count`]; the source checks this as an
    /// `OPENMS_PRECONDITION`, which throws `Exception::Precondition` in debug
    /// builds only.
    pub fn from_index(index: usize) -> Result<AA> {
        if index < Self::unambiguous_aa_count() {
            Ok(AA(index as u8))
        } else {
            Err(bad(
                "AA::from_index: index must be below unambiguous_aa_count()",
            ))
        }
    }

    /// The number of unambiguous amino acids (`A`-`Z` without `B`, `J`, `Z`,
    /// `X`): 22.
    pub const fn unambiguous_aa_count() -> usize {
        22
    }

    /// The residue range an ambiguous amino acid stands for (source
    /// `_ambiguitiesOf` in the `.cpp`): `B` is `D`..`N`, `J` is `I`..`L`, `Z` is
    /// `E`..`Q`, `X` is `A`..`V`, and `$` is `A`..`X` (every residue, including
    /// the ambiguous ones, since a needle may contain an `X`).
    fn ambiguities(self) -> (AA, AA) {
        match self.0 {
            22 => (AA::from_byte(b'D'), AA::from_byte(b'N')),
            23 => (AA::from_byte(b'I'), AA::from_byte(b'L')),
            24 => (AA::from_byte(b'E'), AA::from_byte(b'Q')),
            25 => (AA::from_byte(b'A'), AA::from_byte(b'V')),
            _ => (AA::from_byte(b'A'), AA::from_byte(b'X')),
        }
    }
}

/// A needle found in the query (source struct `Hit`).
///
/// The needle with index `needle_index` (in the order of
/// [`ACTrie::add_needle`] calls) and length `needle_length` was found in the
/// query at `query_pos`. The source's nested `T` is `u32`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Hit {
    /// Index of the needle, in the order it was added.
    pub needle_index: u32,
    /// Length of the needle, in residues.
    pub needle_length: u32,
    /// Start of the hit in the query.
    ///
    /// Computed as the end of the match minus the needle length, so it is off
    /// by the number of skipped invalid characters (such as `*`) that lie
    /// inside the match, as in the source.
    pub query_pos: u32,
}

impl Hit {
    /// A hit (source three-argument constructor).
    pub const fn new(needle_index: u32, needle_length: u32, query_pos: u32) -> Self {
        Self {
            needle_index,
            needle_length,
            query_pos,
        }
    }
}

/// A 32-bit node index that can be invalid (source class `Index`).
///
/// The default is invalid (`u32::MAX`). `Hash` stands for the source's
/// `std::hash<Index>` specialisation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Index(u32);

impl Default for Index {
    fn default() -> Self {
        Index(u32::MAX)
    }
}

impl From<u32> for Index {
    fn from(value: u32) -> Self {
        Index(value)
    }
}

impl Index {
    /// An index with value `value` (source constructor from `T`).
    pub const fn new(value: u32) -> Self {
        Index(value)
    }

    /// Whether this index is invalid and must not be dereferenced.
    pub const fn is_invalid(self) -> bool {
        self.0 == u32::MAX
    }

    /// Whether this index is valid, i.e. an actual index into a vector.
    pub const fn is_valid(self) -> bool {
        self.0 != u32::MAX
    }

    /// The value, which may be the invalid marker (source `operator()` and
    /// const `pos()`).
    pub const fn pos(self) -> u32 {
        self.0
    }

    /// Mutable access to the value (source non-const `pos()`).
    pub fn pos_mut(&mut self) -> &mut u32 {
        &mut self.0
    }

    fn at(self) -> usize {
        self.0 as usize
    }
}

/// A bitset stored in one `u32` (source struct `Bitset`).
///
/// Positions are 0 based. The source shifts `1u << pos` with no check, which is
/// undefined for positions of 32 or more; here [`Bitset::set`] and
/// [`Bitset::clear`] ignore such positions and [`Bitset::test`] reports them
/// unset.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Bitset {
    /// The bits.
    pub bits: u32,
}

impl Bitset {
    /// Set the bit at `pos`.
    pub fn set(&mut self, pos: u32) {
        self.bits |= 1u32.checked_shl(pos).unwrap_or(0);
    }

    /// Clear the bit at `pos`.
    pub fn clear(&mut self, pos: u32) {
        self.bits &= !1u32.checked_shl(pos).unwrap_or(0);
    }

    /// Whether the bit at `pos` is set.
    pub fn test(&self, pos: u32) -> bool {
        self.bits.checked_shr(pos).unwrap_or(0) & 1 == 1
    }

    /// The number of set bits (source `pop_count`).
    pub fn pop_count(&self) -> u32 {
        self.bits.count_ones()
    }
}

impl Shl<u32> for Bitset {
    type Output = Bitset;
    /// Shift left; 32 or more positions clear every bit, as the source
    /// special-cases.
    fn shl(mut self, n: u32) -> Bitset {
        self <<= n;
        self
    }
}

impl ShlAssign<u32> for Bitset {
    fn shl_assign(&mut self, n: u32) {
        self.bits = self.bits.checked_shl(n).unwrap_or(0);
    }
}

/// The depth of a node and its hit flag (source `ACNode::DepthHits`, a
/// one-bit and a seven-bit field).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct DepthHits {
    /// Whether a needle ends here, or at any node on this node's suffix chain.
    pub has_hit: bool,
    /// Depth of the node in the trie, at most [`MAX_NEEDLE_LENGTH`] (seven bits
    /// in the source).
    pub depth: u8,
}

/// The number of children a node can record (source `ACNode::ChildCountType`).
pub type ChildCountType = u8;

/// Longest needle an [`ACTrie`] accepts, in residues.
///
/// The source stores a node's depth in a seven-bit field, so a needle of 128
/// residues or more wraps its depth and reports hits at wrong positions; this
/// port refuses such needles.
pub const MAX_NEEDLE_LENGTH: usize = 127;

/// A node of the trie (source struct `ACNode`).
///
/// After [`ACTrie::compress_trie`] the trie is stored in breadth-first order,
/// so a node's children are the `nr_children` nodes from `first_child` on, in
/// edge-label order, and `children_bitset` records their labels.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ACNode {
    /// Labels of the children (breadth-first trie only).
    pub children_bitset: Bitset,
    /// The suffix link.
    pub suffix: Index,
    /// Index of the first child (breadth-first trie only).
    pub first_child: Index,
    /// Label of the edge from the parent to this node.
    pub edge: AA,
    /// Number of children (breadth-first trie only).
    pub nr_children: ChildCountType,
    /// Depth and hit flag.
    pub depth_and_hits: DepthHits,
}

impl Default for ACNode {
    /// Suffix and first child 0, edge `?`, no children, depth 0, no hit: the
    /// source's member initialisers (`AA edge {0}` is `AA(char 0)`, the invalid
    /// amino acid).
    fn default() -> Self {
        Self {
            children_bitset: Bitset::default(),
            suffix: Index(0),
            first_child: Index(0),
            edge: AA::INVALID,
            nr_children: 0,
            depth_and_hits: DepthHits::default(),
        }
    }
}

impl ACNode {
    /// A node reached by edge `label` at `depth` (source two-argument
    /// constructor).
    pub fn new(label: AA, depth: u8) -> Self {
        Self {
            edge: label,
            depth_and_hits: DepthHits {
                has_hit: false,
                depth,
            },
            ..Self::default()
        }
    }
}

/// A search path spun off the primary to explore an ambiguous residue or a
/// mismatch (source struct `ACScout`).
///
/// The source has no default constructor and holds a raw pointer into the
/// query; here `it_query` is the offset into the query of the
/// [`ACTrieState`] the scout belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ACScout {
    /// Offset of the scout's next query byte.
    pub it_query: usize,
    /// The scout's node in the trie.
    pub tree_pos: Index,
    /// Ambiguous amino acids the scout may still consume.
    pub max_aaa_leftover: u8,
    /// Mismatches the scout may still consume.
    pub max_mm_leftover: u8,
    /// Residues the scout may still lose by following suffix links; reaching
    /// zero retires it, because it has lost the residue that spawned it.
    pub max_prefix_loss_leftover: u8,
}

impl ACScout {
    /// A scout (source constructor).
    pub const fn new(
        query_pos: usize,
        tree_pos: Index,
        max_aa: u8,
        max_mm: u8,
        max_prefix_loss: u8,
    ) -> Self {
        Self {
            it_query: query_pos,
            tree_pos,
            max_aaa_leftover: max_aa,
            max_mm_leftover: max_mm,
            max_prefix_loss_leftover: max_prefix_loss,
        }
    }

    /// Where in the query the scout is (source `textPos(state)`, which
    /// subtracts the query start from the scout's pointer).
    pub const fn text_pos(&self) -> usize {
        self.it_query
    }

    /// The scout's next valid amino acid in `query`, advancing past it; the
    /// invalid amino acid when the query is exhausted (source `nextValidAA()`,
    /// which reads through its stored pointer). `query` is the state's query,
    /// [`ACTrieState::get_query`].
    pub fn next_valid_aa(&mut self, query: &str) -> AA {
        next_valid_aa(query.as_bytes(), &mut self.it_query)
    }
}

/// The first valid amino acid at or after `*pos` in `query`, leaving `*pos`
/// just after it; or the invalid amino acid, leaving `*pos` at the end
/// (source free function `nextValidAA`).
///
/// The source walks a NUL-terminated C string, so a NUL byte ends the query
/// early; this port stops at a NUL byte too.
pub fn next_valid_aa(query: &[u8], pos: &mut usize) -> AA {
    let mut result = AA::INVALID;
    while let Some(&byte) = query.get(*pos) {
        if byte == 0 {
            break;
        }
        result = AA::from_byte(byte);
        *pos += 1;
        if result.is_valid() {
            break;
        }
    }
    result
}

/// The dynamic state of one walk through an [`ACTrie`] (source struct
/// `ACTrieState`).
///
/// The trie is immutable once compressed, so several threads can each walk it
/// with their own state.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ACTrieState {
    /// The hits found by the last call.
    pub hits: Vec<Hit>,
    /// Scouts waiting to be processed.
    pub scouts: VecDeque<ACScout>,
    /// The primary's node in the trie.
    pub tree_pos: Index,
    it_q: usize,
    query: String,
}

impl ACTrieState {
    /// Set the query (haystack) to search, reset the primary to the root, and
    /// clear the hits and scouts (source `setQuery`).
    pub fn set_query(&mut self, haystack: &str) {
        self.hits.clear();
        self.query = haystack.to_owned();
        self.it_q = 0;
        self.tree_pos = Index(0);
        self.scouts.clear();
    }

    /// Where in the query the primary is (source `textPos`).
    pub fn text_pos(&self) -> usize {
        self.it_q
    }

    /// The rest of the query from the primary's position (source
    /// `textPosIt`, which returns a pointer to it).
    pub fn text_pos_it(&self) -> &[u8] {
        &self.query.as_bytes()[self.it_q..]
    }

    /// The current query (source `getQuery`).
    pub fn get_query(&self) -> &str {
        &self.query
    }

    /// The primary's next valid amino acid, advancing past it; the invalid
    /// amino acid when the query is exhausted (source `nextValidAA`).
    pub fn next_valid_aa(&mut self) -> AA {
        next_valid_aa(self.query.as_bytes(), &mut self.it_q)
    }
}

/// An Aho-Corasick trie over amino-acid needles with ambiguity and mismatch
/// support (source class `ACTrie`).
///
/// Build it with [`ACTrie::add_needle`] (or the list variants), call
/// [`ACTrie::compress_trie`] once, then search with [`ACTrie::next_hits`] or
/// [`ACTrie::get_all_hits`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ACTrie {
    trie: Vec<ACNode>,
    needle_count: u32,
    max_aaa: u32,
    max_mm: u32,
    index2needles: Vec<Vec<u32>>,
    index2children_naive: Vec<Vec<Index>>,
    compressed: bool,
}

impl Default for ACTrie {
    fn default() -> Self {
        Self::new(0, 0)
    }
}

impl ACTrie {
    /// A trie holding only its root (source constructor).
    ///
    /// # Arguments
    ///
    /// * `max_aaa` - most ambiguous amino acids (`B`, `J`, `Z`, `X`) allowed in
    ///   a hit (source default 0)
    /// * `max_mm` - most mismatched amino acids allowed in a hit (source
    ///   default 0)
    pub fn new(max_aaa: u32, max_mm: u32) -> Self {
        Self {
            trie: vec![ACNode::default()],
            needle_count: 0,
            max_aaa,
            max_mm,
            index2needles: vec![Vec::new()],
            index2children_naive: vec![Vec::new()],
            compressed: false,
        }
    }

    /// Add a needle. Call [`ACTrie::compress_trie`] after the last needle and
    /// before searching.
    ///
    /// Needles may contain the ambiguous `B`, `J`, `Z`, `X`, which then match
    /// only themselves in the query (or anything, through the query's own
    /// ambiguity and mismatch budget); letters are case-insensitive.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the needle contains a character
    /// that is not a letter (such as `*` or `$`; the source
    /// `Exception::InvalidValue`), when it is empty (the source adds it, and
    /// its hit collection then loops forever at the root), when it is longer
    /// than [`MAX_NEEDLE_LENGTH`] (the source's seven-bit depth wraps), when
    /// the trie has already been compressed (the source would corrupt it), or
    /// when the node or needle count would exceed `u32`.
    ///
    /// The needle is checked before anything is added, so a refused needle
    /// leaves the trie unchanged; the source has already created the nodes of
    /// the needle's valid prefix when it throws, which [`ACTrie::size`] would
    /// show.
    pub fn add_needle(&mut self, needle: &str) -> Result<()> {
        if self.compressed {
            return Err(bad("ACTrie: needles cannot be added after compress_trie()"));
        }
        let bytes = needle.as_bytes();
        if let Some(&invalid) = bytes
            .iter()
            .find(|&&c| !AA::from_byte(c).is_valid_for_peptide())
        {
            return Err(bad(&format!(
                "Invalid amino acid: '{}'",
                String::from_utf8_lossy(&[invalid])
            )));
        }
        if bytes.is_empty() {
            return Err(bad("ACTrie: a needle must not be empty"));
        }
        if bytes.len() > MAX_NEEDLE_LENGTH {
            return Err(bad("ACTrie: needle longer than MAX_NEEDLE_LENGTH"));
        }
        if self.needle_count == u32::MAX
            || self.trie.len().saturating_add(bytes.len()) >= u32::MAX as usize
        {
            return Err(bad("ACTrie: node or needle count exceeds u32"));
        }
        let mut cn = Index(0);
        for &c in bytes {
            cn = self.add_child(cn, AA::from_byte(c));
        }
        self.trie[cn.at()].depth_and_hits.has_hit = true;
        if self.index2needles.len() <= cn.at() {
            self.index2needles.resize_with(cn.at() + 1, Vec::new);
        }
        self.index2needles[cn.at()].push(self.needle_count);
        self.needle_count += 1;
        Ok(())
    }

    /// Add several needles in order (source `addNeedles`).
    ///
    /// # Errors
    ///
    /// As [`ACTrie::add_needle`], for the first needle refused. The needles
    /// before it stay added, as in the source; [`ACTrie::get_needle_count`]
    /// tells which needle was refused.
    pub fn add_needles<S: AsRef<str>>(&mut self, needles: &[S]) -> Result<()> {
        for needle in needles {
            self.add_needle(needle.as_ref())?;
        }
        Ok(())
    }

    /// Add several needles, then compress the trie so no more can be added
    /// (source `addNeedlesAndCompress`).
    ///
    /// # Errors
    ///
    /// As [`ACTrie::add_needles`]; the trie is not compressed on error, as
    /// the source's exception skips the compression too.
    pub fn add_needles_and_compress<S: AsRef<str>>(&mut self, needles: &[S]) -> Result<()> {
        self.add_needles(needles)?;
        self.compress_trie();
        Ok(())
    }

    /// The number of nodes.
    pub fn size(&self) -> usize {
        self.trie.len()
    }

    /// Lay the trie out in breadth-first order, which is more compact and
    /// faster to traverse, and create the suffix links (source
    /// `compressTrie`).
    ///
    /// Call it after adding every needle and before searching. Suffix links
    /// are compressed where a suffix node has no hit and no child the node
    /// lacks, so a search skips it.
    ///
    /// A second call does nothing. The source rebuilds from its discarded
    /// construction lists on a second call, which leaves only the root with
    /// its old child bits, and a later search reads out of bounds. The source
    /// also logs the node count and the number of compressed links at info
    /// level; this port does not log.
    pub fn compress_trie(&mut self) {
        if self.compressed {
            return;
        }
        let n = self.trie.len();
        let mut bfs_tree: Vec<ACNode> = Vec::with_capacity(n);
        let mut bfs_needles: Vec<Vec<u32>> = vec![Vec::new(); n];
        let mut tmp_parents: Vec<Index> = Vec::with_capacity(n);
        let mut queue: VecDeque<Index> = VecDeque::new();
        let mut index2needles = std::mem::take(&mut self.index2needles);
        let mut children_naive = std::mem::take(&mut self.index2children_naive);
        let naive = &self.trie;

        let mut bfs_op = |current: Index,
                          bfs_tree: &mut Vec<ACNode>,
                          tmp_parents: &mut Vec<Index>,
                          queue: &mut VecDeque<Index>| {
            let bfs_index = bfs_tree.len() as u32;
            let mut node = naive[current.at()];
            if let Some(children) = children_naive.get_mut(current.at()) {
                if !children.is_empty() {
                    node.nr_children = children.len() as ChildCountType;
                    // Children in edge-label order, the order of the bitset.
                    children.sort_by_key(|child| naive[child.at()].edge);
                    for &child in children.iter() {
                        node.children_bitset
                            .set(u32::from(naive[child.at()].edge.code()));
                        queue.push_back(child);
                        tmp_parents.push(Index(bfs_index));
                    }
                }
            }
            bfs_tree.push(node);
            if let Some(needles) = index2needles.get_mut(current.at()) {
                bfs_needles[bfs_index as usize] = std::mem::take(needles);
            }
        };

        tmp_parents.push(Index(0)); // the root is its own parent
        bfs_op(Index(0), &mut bfs_tree, &mut tmp_parents, &mut queue);
        bfs_tree[0].first_child = Index(1);
        while let Some(current) = queue.pop_front() {
            let last = bfs_tree[bfs_tree.len() - 1];
            bfs_op(current, &mut bfs_tree, &mut tmp_parents, &mut queue);
            let first_child = last.first_child.pos() + u32::from(last.nr_children);
            let new_len = bfs_tree.len();
            bfs_tree[new_len - 1].first_child = Index(first_child);
        }

        self.trie = bfs_tree;
        self.index2needles = bfs_needles;
        self.trie[0].suffix = Index(0);

        // Suffix links, without path compression. Depth 1 already points to
        // the root.
        let first_deep = 1 + usize::from(self.trie[0].nr_children);
        for (i, &parent) in tmp_parents.iter().enumerate().skip(first_deep) {
            let suffix = self.follow(self.trie[parent.at()].suffix, self.trie[i].edge);
            self.trie[i].suffix = suffix;
            let suffix_hit = self.trie[suffix.at()].depth_and_hits.has_hit;
            self.trie[i].depth_and_hits.has_hit |= suffix_hit;
        }
        // Path compression, where the suffix has no hit (there are no separate
        // output links) and no child this node lacks.
        for i in first_deep..self.trie.len() {
            let suffix = self.trie[i].suffix;
            if self.trie[suffix.at()].depth_and_hits.has_hit {
                continue;
            }
            let children_node = self.trie[i].children_bitset.bits;
            let children_suffix = self.trie[suffix.at()].children_bitset.bits;
            if (children_node | children_suffix) == children_node {
                self.trie[i].suffix = self.trie[suffix.at()].suffix;
            }
        }
        self.compressed = true;
    }

    /// How many needles were added.
    pub fn get_needle_count(&self) -> usize {
        self.needle_count as usize
    }

    /// Set the most ambiguous amino acids allowed in a hit. Not to be called
    /// in the middle of a search, or results mix.
    pub fn set_max_aaa_count(&mut self, max_aaa: u32) {
        self.max_aaa = max_aaa;
    }

    /// The most ambiguous amino acids allowed in a hit.
    pub fn get_max_aaa_count(&self) -> u32 {
        self.max_aaa
    }

    /// Set the most mismatches allowed in a hit. Not to be called in the
    /// middle of a search, or results mix.
    pub fn set_max_mm_count(&mut self, max_mm: u32) {
        self.max_mm = max_mm;
    }

    /// The most mismatches allowed in a hit.
    pub fn get_max_mm_count(&self) -> u32 {
        self.max_mm
    }

    fn check_searchable(&self, state: &ACTrieState) -> Result<()> {
        if !self.compressed {
            return Err(bad(
                "ACTrie: compress_trie() must be called before searching",
            ));
        }
        if state.query.len() > u32::MAX as usize {
            return Err(bad("ACTrie: query longer than u32::MAX"));
        }
        if state.tree_pos.is_invalid() || state.tree_pos.at() >= self.trie.len() {
            return Err(bad(
                "ACTrie: state is not positioned in this trie (call ACTrieState::set_query first)",
            ));
        }
        Ok(())
    }

    /// Resume the search where `state` stopped; when the primary's node (or
    /// its suffixes) holds a hit, clear `state.hits`, fill it and return
    /// `true`; return `false` when the query ends without one (source
    /// `nextHits`).
    ///
    /// Hits found by scouts accumulate until the primary next reports hits,
    /// and are returned with them; scout hits found after the primary's last
    /// hit are left in `state.hits` with a `false` return, exactly as in the
    /// source (so a `while next_hits` loop drops them, as the source's does).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the trie was not compressed (a
    /// debug-only assertion in the source), when the query is longer than
    /// `u32::MAX` (hit positions are `u32`), or when `state` is positioned
    /// outside this trie.
    pub fn next_hits(&self, state: &mut ACTrieState) -> Result<bool> {
        self.check_searchable(state)?;
        state.hits.clear();
        self.next_hits_no_clear(state);
        Ok(!state.hits.is_empty())
    }

    /// Collect every hit from the current position to the end of the query
    /// into `state.hits` (source `getAllHits`).
    ///
    /// # Errors
    ///
    /// As [`ACTrie::next_hits`].
    pub fn get_all_hits(&self, state: &mut ACTrieState) -> Result<()> {
        self.check_searchable(state)?;
        state.hits.clear();
        while self.next_hits_no_clear(state) {}
        Ok(())
    }

    fn next_hits_no_clear(&self, state: &mut ACTrieState) -> bool {
        loop {
            let aa = state.next_valid_aa();
            if !aa.is_valid() {
                return false;
            }
            state.tree_pos = self.step_primary(state.tree_pos, aa, state);
            // Process the scouts now, while the caches are hot; each runs to
            // the end of its life and may queue more.
            while let Some(mut scout) = state.scouts.pop_front() {
                while self.step_scout(&mut scout, state) {}
            }
            let text_pos = state.text_pos();
            if self.add_hits(state.tree_pos, text_pos, &mut state.hits) {
                return true;
            }
        }
    }

    fn add_child(&mut self, from: Index, edge: AA) -> Index {
        if self.index2children_naive.len() <= from.at() {
            self.index2children_naive
                .resize_with(from.at() + 1, Vec::new);
        }
        if let Some(child) = self.find_child_naive(from, edge) {
            return child;
        }
        let child = Index(self.trie.len() as u32);
        let depth = self.trie[from.at()].depth_and_hits.depth + 1;
        self.trie.push(ACNode::new(edge, depth));
        self.index2children_naive[from.at()].push(child);
        child
    }

    fn find_child_naive(&self, parent: Index, label: AA) -> Option<Index> {
        self.index2children_naive[parent.at()]
            .iter()
            .copied()
            .find(|child| self.trie[child.at()].edge == label)
    }

    fn find_child_bfs(&self, parent: Index, label: AA) -> Index {
        let node = &self.trie[parent.at()];
        let bit = u32::from(label.code());
        if !node.children_bitset.test(bit) {
            return Index::default();
        }
        let offset = count_set_bits_up_to(node.children_bitset, bit);
        Index(node.first_child.pos() + offset)
    }

    fn follow(&self, mut i: Index, aa: AA) -> Index {
        loop {
            let child = self.find_child_bfs(i, aa);
            if child.is_valid() {
                return child;
            }
            if i.pos() == 0 {
                return Index(0);
            }
            i = self.trie[i.at()].suffix;
        }
    }

    fn collect_hits(
        &self,
        mut i: Index,
        text_pos: usize,
        hits: &mut Vec<Hit>,
        scout: Option<(&ACScout, i32)>,
    ) -> bool {
        let before = hits.len();
        loop {
            let node = &self.trie[i.at()];
            if !node.depth_and_hits.has_hit {
                break;
            }
            let needle_length = node.depth_and_hits.depth;
            if let Some((scout, current_depth)) = scout {
                // A scout reports only hits that include its first ambiguous
                // or mismatched residue; shorter ones are reported by others.
                if current_depth - i32::from(needle_length)
                    >= i32::from(scout.max_prefix_loss_leftover)
                {
                    break;
                }
            }
            // text_pos never exceeds u32::MAX (checked before the search) and
            // is at least the depth, since every level consumed one residue.
            let text_start = text_pos.saturating_sub(usize::from(needle_length)) as u32;
            for &needle in &self.index2needles[i.at()] {
                hits.push(Hit::new(needle, u32::from(needle_length), text_start));
            }
            // The root never has a hit (empty needles are refused), so this
            // chain ends.
            i = node.suffix;
        }
        hits.len() != before
    }

    fn add_hits(&self, i: Index, text_pos: usize, hits: &mut Vec<Hit>) -> bool {
        self.collect_hits(i, text_pos, hits, None)
    }

    fn follow_scout(&self, scout: &mut ACScout, edge: AA, state: &mut ACTrieState) -> bool {
        let j = self.follow(scout.tree_pos, edge);
        let new_depth = i32::from(self.trie[j.at()].depth_and_hits.depth);
        let old_depth = i32::from(self.trie[scout.tree_pos.at()].depth_and_hits.depth);
        let up_count = old_depth - new_depth + 1;
        if up_count >= i32::from(scout.max_prefix_loss_leftover) {
            return false; // the scout lost its ambiguous or mismatched residue
        }
        // 0 <= up_count < max_prefix_loss_leftover <= 255.
        scout.max_prefix_loss_leftover -= up_count as u8;
        scout.tree_pos = j;
        self.collect_hits(
            j,
            scout.text_pos(),
            &mut state.hits,
            Some((scout, new_depth)),
        );
        true
    }

    fn step_primary(&self, i: Index, edge: AA, state: &mut ACTrieState) -> Index {
        let consider_aaa = self.max_aaa != 0;
        let consider_mm = self.max_mm != 0;
        if edge.is_ambiguous() {
            let (mut from, mut to) = (edge, edge);
            if consider_aaa {
                (from, to) = edge.ambiguities();
                self.create_scouts(i, from, to, state, self.max_aaa - 1, self.max_mm);
            }
            if consider_mm {
                self.create_mm_scouts(i, from, to, edge, state, self.max_aaa, self.max_mm - 1);
            }
        } else if consider_mm {
            self.create_mm_scouts(i, edge, edge, edge, state, self.max_aaa, self.max_mm - 1);
        }
        // The primary continues with the residue, whatever it was.
        self.follow(i, edge)
    }

    fn step_scout(&self, scout: &mut ACScout, state: &mut ACTrieState) -> bool {
        loop {
            let edge = next_valid_aa(state.query.as_bytes(), &mut scout.it_query);
            if !edge.is_valid() {
                return false; // end of the query
            }
            let consider_aaa = scout.max_aaa_leftover > 0;
            let consider_mm = scout.max_mm_leftover > 0;
            if edge.is_ambiguous() {
                let (mut from, mut to) = (edge, edge);
                if consider_aaa {
                    (from, to) = edge.ambiguities();
                    let mut prototype = *scout;
                    prototype.max_aaa_leftover -= 1;
                    self.create_sub_scouts(&prototype, from, to, state);
                }
                if consider_mm {
                    let mut prototype = *scout;
                    prototype.max_mm_leftover -= 1;
                    self.create_mm_sub_scouts(&prototype, from, to, edge, state);
                }
            } else if consider_mm {
                let mut prototype = *scout;
                prototype.max_mm_leftover -= 1;
                self.create_mm_sub_scouts(&prototype, edge, edge, edge, state);
            }
            if !self.follow_scout(scout, edge, state) {
                return false;
            }
        }
    }

    /// Every residue from `A` to `X` except the range `except_from..=except_to`
    /// and `except_edge`, in code order (the source's mismatch loops).
    fn mismatch_residues(except_from: AA, except_to: AA, except_edge: AA) -> Vec<AA> {
        let (from, to) = AA::SUPER.ambiguities();
        let mut residues = Vec::with_capacity(26);
        let mut code = from.code();
        while code <= to.code() {
            let aa = AA(code);
            if aa == except_from {
                code = except_to.code() + 1;
                continue;
            }
            if aa != except_edge {
                residues.push(aa);
            }
            code += 1;
        }
        residues
    }

    #[allow(clippy::too_many_arguments)]
    fn create_mm_scouts(
        &self,
        i: Index,
        except_from: AA,
        except_to: AA,
        except_edge: AA,
        state: &mut ACTrieState,
        aaa_left: u32,
        mm_left: u32,
    ) {
        for aa in Self::mismatch_residues(except_from, except_to, except_edge) {
            self.create_scouts(i, aa, aa, state, aaa_left, mm_left);
        }
    }

    fn create_mm_sub_scouts(
        &self,
        prototype: &ACScout,
        except_from: AA,
        except_to: AA,
        except_edge: AA,
        state: &mut ACTrieState,
    ) {
        for aa in Self::mismatch_residues(except_from, except_to, except_edge) {
            self.create_sub_scouts(prototype, aa, aa, state);
        }
    }

    fn create_scouts(
        &self,
        i: Index,
        from: AA,
        to: AA,
        state: &mut ACTrieState,
        aaa_left: u32,
        mm_left: u32,
    ) {
        for code in from.code()..=to.code() {
            let scout_pos = self.follow(i, AA(code));
            if scout_pos.pos() > 0 {
                let depth = self.trie[scout_pos.at()].depth_and_hits.depth;
                // The primary already points past the residue, so the scout
                // starts there. Budgets above 255 saturate (the source
                // narrows them to uint8_t, wrapping).
                let scout = ACScout::new(
                    state.it_q,
                    scout_pos,
                    u8::try_from(aaa_left).unwrap_or(u8::MAX),
                    u8::try_from(mm_left).unwrap_or(u8::MAX),
                    depth,
                );
                state.scouts.push_back(scout);
                self.add_hits(scout_pos, scout.text_pos(), &mut state.hits);
            }
        }
    }

    fn create_sub_scouts(&self, prototype: &ACScout, from: AA, to: AA, state: &mut ACTrieState) {
        for code in from.code()..=to.code() {
            let mut scout = *prototype;
            if self.follow_scout(&mut scout, AA(code), state) {
                state.scouts.push_back(scout);
            }
        }
    }
}

/// The number of set bits below position `i` (the `.cpp`'s
/// `countSetBitsUpTo`).
fn count_set_bits_up_to(bitset: Bitset, i: u32) -> u32 {
    if i == 0 {
        return 0;
    }
    (bitset << (32 - i)).pop_count()
}
