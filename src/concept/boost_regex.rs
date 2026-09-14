// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Boost.Regex-compatible regular expressions, executed by `fancy-regex`.
//!
//! OpenMS takes its regular expressions from Boost.Regex: enzyme and RNase
//! cleavage rules, native-ID scan extraction, Mascot title formats, mzTab
//! column names, the indexedmzML footer, mzIdentML fragment annotations, the
//! pepXML enzyme summary, decoy affixes and ion-name grammars. Several of those
//! expressions need lookaround or repeated group names, which the `regex` crate
//! does not have, so the engine is `fancy-regex`; expressions without such
//! features still run on its `regex-automata` delegate.
//!
//! A bare `fancy-regex` pattern does not behave like Boost.
//! [`BoostRegex`](crate::concept::boost_regex::BoostRegex)
//! parses the Boost perl syntax itself and hands the engine an equivalent
//! pattern in which every Boost convention is spelled out:
//!
//! - `\d`, `\w`, `\s`, `\h`, `\v`, the POSIX classes and bracket expressions
//!   are ASCII, as Boost's `char` traits in the C locale are. The translator
//!   builds each set the way Boost's `basic_regex_creator` fills its byte map,
//!   including its lower-cased range endpoints under `icase` and its single
//!   complement of all negated classes in one bracket expression, and emits the
//!   resulting bytes, so neither the engine's class tables nor its case folding
//!   is consulted for a set;
//! - case-insensitive literals and backreferences, and the word boundaries
//!   `\b`, `\B`, `\<` and `\>`, use the engine's Unicode tables, which agree
//!   with Boost's C-locale tables on every character a haystack can contain here
//!   (ASCII and the transcoded bytes described below);
//! - `.` matches every byte, including line separators, unless `(?-s)` or
//!   [`RegexOptions::no_mod_s`](crate::concept::boost_regex::RegexOptions::no_mod_s) is in effect, in which case it excludes `\n`,
//!   `\r` and `\f`;
//! - `^` and `$` are line anchors whose separators are `\n`, `\r` and `\f`
//!   (Boost's `is_separator<char>`), and neither matches between the `\r` and
//!   `\n` of a CRLF pair;
//! - a lookbehind must have one fixed width, computed with Boost's own rules,
//!   or construction fails as Boost's does;
//! - the token iterator reproduces `boost::sregex_token_iterator`, including
//!   `match_not_initial_null` after an empty match.
//!
//! # Byte matching
//!
//! Boost matches bytes. The engine is run in its Unicode mode instead, over a
//! transcoded haystack: every byte above `0x7F` becomes its own private-use
//! code point (`U+E080` to `U+E0FF`), and every reported offset is mapped back
//! to a byte offset. A pure-ASCII haystack is used as it is. Pattern bytes
//! outside comments are restricted to ASCII, so a transcoded byte is, for every
//! construct the facade accepts, exactly what the original byte is to Boost: one character that no
//! class, literal or case fold names, a non-word character, and distinct from
//! every other byte for backreferences. The detour exists because
//! `fancy-regex` configures `regex-automata` to refuse empty matches inside a
//! UTF-8 sequence even in its bytes mode, which loses Boost's empty matches
//! between the bytes of a multi-byte character and makes `regex-automata`
//! panic on a lone continuation byte; with a transcoded haystack every byte
//! boundary is a character boundary.
//!
//! # Work bounds
//!
//! A search runs either on an automaton or on `fancy-regex`'s backtracking
//! machine, and each has a bound.
//!
//! An expression that needs no backtracking, and whose automaton has at most
//! [`MAX_AUTOMATON_ATOMS`](crate::concept::boost_regex::MAX_AUTOMATON_ATOMS)
//! atoms (divided by the budget's scale factor for a long haystack), is handed
//! whole to an automaton. It spends no budget, and its work is proportional to
//! the haystack length times the atoms, which grow with repeat bounds.
//!
//! Every other search runs on the backtracking machine, which fails a search
//! that takes more backtracking steps than its budget
//! ([`RegexOptions::backtrack_limit`](crate::concept::boost_regex::RegexOptions::backtrack_limit),
//! shared by all start positions of one search and scaled with the haystack
//! length), or that needs more than the 1,000,000 branches of its fixed stack
//! (`MAX_STACK`, reported by the engine as `RuntimeError::StackOverflow`). Both
//! are [`Error::InvalidValue`](crate::Error::InvalidValue): an error, never a
//! different answer. The engine does not count the work it does between two
//! backtracking steps, so the facade spells these expressions such that the
//! budget pays for that work as well:
//!
//! - the expression and every positive lookahead in it end in an always-true
//!   assertion, so the engine runs their sub-expressions on the backtracking
//!   machine, where every loop iteration pushes a branch, instead of handing a
//!   trailing run to an automaton that could scan to the end of the haystack
//!   from every start position;
//! - every iteration of a repeat whose minimum is above one pushes a branch,
//!   which the engine otherwise does only above the minimum;
//! - a repeat inside an atomic group or a negative lookaround, whose branches
//!   the engine discards without counting them, is refused, and so are a
//!   bounded repeat of something that can match the empty string, a lookbehind
//!   wider than [`MAX_LOOKBEHIND_WIDTH`](crate::concept::boost_regex::MAX_LOOKBEHIND_WIDTH)
//!   and a translation longer than
//!   [`MAX_TRANSLATED_BYTES`](crate::concept::boost_regex::MAX_TRANSLATED_BYTES).
//!
//! Between two backtracking steps the engine then does work proportional to the
//! translated pattern and its lookbehind widths, plus the length of the group a
//! backreference names, so a search's work is bounded by its budget times that.
//! A case-sensitive backreference stops comparing at the first byte that
//! differs, as Boost's does; a case-insensitive one first compares exactly and,
//! when that fails, scans its whole group before it compares without case, so
//! it costs the group's length on every attempt even when the first byte
//! differs. With a lazily growing group that makes one search quadratic within
//! its budget where Boost's is not: `([a-c]*?)(?i)[^a](a\1)` over random `a`
//! and `b` took 4.3 s over 1 MiB and 70 s over 4 MiB, where Boost took 0.16 s
//! over 4 MiB.
//!
//! The engine's positive lookaheads are not
//! atomic, as Boost's are: when what follows one fails, the engine backtracks
//! into its body again, so a body with ambiguous repeats followed by a failing
//! continuation can exhaust the budget on a few dozen bytes where Boost answers
//! at once (`(?=(?:a|aa)*)b` on 24 `a` and a `c`). These limits are not Boost's
//! `error_complexity`, so each engine stops some searches the other answers;
//! `docs/BOOST_REGEX_SUPPORT.md` records the measurements.
//!
//! # Engine rewrites
//!
//! Before they compile an expression, `fancy-regex` and `regex-syntax` rewrite
//! some shapes into forms that match the same strings but not with
//! leftmost-first priority or the same captures: adjacent and nested unbounded
//! repeats (`fancy-regex`'s optimizer) and a prefix that every branch of an
//! alternation shares (`regex-syntax`, in what `fancy-regex` hands to an
//! automaton: a whole expression that needs no backtracking, or the body of an
//! atomic group, which then keeps a different match). The facade spells
//! its translation so that neither rewrite applies, checks the result with the
//! engine's own optimizer, and refuses an expression it cannot spell that way.
//!
//! Syntax the facade does not translate is refused with
//! [`Error::Unsupported`](crate::Error::Unsupported), never passed through with different meaning. The
//! executed differential evidence and the list of refused constructs are in
//! `docs/BOOST_REGEX_SUPPORT.md`.

use crate::error::{Error, Result};
use fancy_regex::{BytesMode, RegexBuilder, RegexInput};
use std::fmt::Write as _;
use std::iter::FusedIterator;
use std::ops::Range;
use std::sync::{Arc, OnceLock};

/// Longest pattern, in bytes, that [`BoostRegex::with_options`] accepts.
pub const MAX_PATTERN_BYTES: usize = 64 * 1024;

/// Deepest group nesting the translator accepts.
///
/// `fancy-regex`'s parser refuses a group nested 64 deep, and the translation
/// nests deeper than the pattern: the programs wrap it in two groups, a
/// case-insensitive letter, `.`, a line anchor or a counted, shielded or
/// rewritten repeat adds levels of its own (see "Work bounds" and "Engine
/// rewrites" in the module documentation). So besides this cap on the pattern's
/// own nesting, the translator counts the levels its spellings could add, as if
/// every repeat were counted and shielded, and refuses a pattern that could reach
/// the engine's limit before the engine sees it.
pub const MAX_GROUP_DEPTH: usize = 48;

/// Largest repeat bound in `{n}`, `{n,}` or `{n,m}` that the translator
/// accepts, and the longest shortest match, in bytes, that a repeated
/// sub-expression may require.
///
/// The second bound exists because `fancy-regex`'s analyzer multiplies the
/// shortest match of a repeated sub-expression by the repeat's minimum without
/// an overflow check: `(?:(?:a{999999999}){999999999}){999999999}` would panic
/// in a build with overflow checks. The translator tracks that product with
/// checked arithmetic and refuses such a pattern with
/// [`Error::Unsupported`] before the engine sees it.
pub const MAX_REPEAT: usize = 999_999_999;

/// Largest automaton, in atoms, that a search over at most
/// [`BACKTRACK_LIMIT_BYTES`] bytes runs on; a longer haystack divides it by the
/// factor that scales its budget (4, 16 or [`MAX_BACKTRACK_SCALE`]).
///
/// The atoms of an expression are its one-byte atoms, assertions and
/// backreferences with every repeat written out: a repeat multiplies the atoms
/// of what it repeats by its maximum, or by its minimum (at least 1) when it
/// is unbounded. The branch that guards an alternation, and the extra branch
/// of a repeat shielded from the engine's optimizer, count as one atom each. An expression that needs no backtracking is searched by
/// `fancy-regex`'s `regex-automata` delegate, which spends no backtracking
/// budget. Its slowest mode, the one it falls back to when its lazy DFA runs
/// out of cache, does work proportional to the haystack length times the atoms:
/// `\w{0,9999}b` took 8 s over 100 KB. An expression with more atoms than the
/// haystack's share of this limit is searched with the counted spelling on the
/// backtracking machine instead, where [`RegexOptions::backtrack_limit`] bounds
/// it, so an automaton search does at most about
/// `MAX_AUTOMATON_ATOMS * BACKTRACK_LIMIT_BYTES` atom steps up to
/// `MAX_BACKTRACK_SCALE * BACKTRACK_LIMIT_BYTES` bytes, and
/// `MAX_AUTOMATON_ATOMS / MAX_BACKTRACK_SCALE` per byte beyond.
///
/// An atom step is not a fixed amount of time. A byte above `0x7F` is three
/// bytes of the transcoded haystack (see the module documentation), and a
/// class that accepts such bytes (`.`, `\W`, `[^b]`) is several automaton
/// states, so a haystack of such bytes takes up to about twice as long per byte
/// at the same atom count. The bound is also per search: a token iterator runs
/// one search per match, so its total work is the bound times the number of
/// matches, and `a.{0,2045}a.{0,2046}c|a` (4,096 atoms) over 4 KiB, which
/// matches 2,073 times, took 44 s where Boost took 1 s.
/// `docs/BOOST_REGEX_SUPPORT.md` records the measured times.
pub const MAX_AUTOMATON_ATOMS: usize = 4096;

/// Widest lookbehind, in bytes, that the translator accepts.
///
/// Every time the engine tries a lookbehind it steps back over the whole width,
/// one character at a time, and matches the body forward again: work that no
/// backtracking step pays for. The cap keeps it small per try.
pub const MAX_LOOKBEHIND_WIDTH: usize = 255;

/// Longest engine pattern, in bytes, that a translation may produce.
///
/// The translation spells out every Boost convention, so it is longer than the
/// pattern: a line anchor becomes about 40 bytes and a case-insensitive letter
/// 9, plus the counting constructs described under "Work bounds" in the module
/// documentation. The engine's memory grows with the translation, and so does
/// the work it can do between two backtracking steps, so a longer translation
/// is refused.
pub const MAX_TRANSLATED_BYTES: usize = 8 * MAX_PATTERN_BYTES;

/// Backtracking steps a search over a haystack of at most
/// [`BACKTRACK_LIMIT_BYTES`] bytes may take before it fails, unless
/// [`RegexOptions::backtrack_limit`] says otherwise.
///
/// Only searches on the backtracking engine spend this budget: those of
/// expressions that need it (lookaround, backreferences, atomic groups, line
/// anchors, word boundaries) and of expressions whose automaton would be larger
/// than [`MAX_AUTOMATON_ATOMS`] allows; the others run on an automaton. One
/// search shares the budget across every start position it tries, so a longer
/// haystack scales it (see [`BACKTRACK_LIMIT_BYTES`]).
pub const DEFAULT_BACKTRACK_LIMIT: usize = 1_000_000;

/// Haystack length, in bytes, that one unscaled backtracking budget covers.
///
/// A search over `n` bytes may take `backtrack_limit * s` steps, where `s` is
/// the smallest of 1, 4, 16 and [`MAX_BACKTRACK_SCALE`] with
/// `n <= s * BACKTRACK_LIMIT_BYTES`, or [`MAX_BACKTRACK_SCALE`] when none is.
pub const BACKTRACK_LIMIT_BYTES: usize = 64 * 1024;

/// Largest factor by which a long haystack multiplies the backtracking budget.
pub const MAX_BACKTRACK_SCALE: usize = 64;

/// The budget factors, one set of engine programs each: `fancy-regex` fixes the
/// backtrack limit when it compiles a program.
const BACKTRACK_SCALES: [usize; 4] = [1, 4, 16, MAX_BACKTRACK_SCALE];

/// Entries of `fancy-regex`'s backtracking stack (its fixed `MAX_STACK`).
const ENGINE_STACK_ENTRIES: usize = 1_000_000;

/// First code point of the private-use block a byte above `0x7F` is transcoded
/// into (`U+E000 + byte`).
const TRANSCODE_BASE: u32 = 0xE000;

/// Construction options, the subset of `boost::regex_constants::syntax_option_type`
/// that OpenMS passes, plus the work bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RegexOptions {
    /// `boost::regex::icase`: ASCII case-insensitive matching from the start of
    /// the pattern. `(?i)` and `(?-i)` can still change it inside the pattern.
    pub icase: bool,
    /// `boost::regex::no_mod_s`: `.` does not match the line separators `\n`,
    /// `\r` and `\f`. `(?s)` can still turn matching of separators back on.
    pub no_mod_s: bool,
    /// Backtracking steps a search over at most [`BACKTRACK_LIMIT_BYTES`] bytes
    /// may take; see [`DEFAULT_BACKTRACK_LIMIT`]. Must be positive.
    ///
    /// A longer haystack multiplies the budget by 4, 16 or
    /// [`MAX_BACKTRACK_SCALE`], the smallest factor `s` with
    /// `len <= s * BACKTRACK_LIMIT_BYTES`, because one search shares the budget
    /// across all its start positions: without scaling, a line-anchored scan of
    /// a few hundred kilobytes runs out of budget where Boost, whose own
    /// complexity limit grows with the input, answers at once. The factor stops
    /// at [`MAX_BACKTRACK_SCALE`]. Exceeding the budget is an error rather than
    /// a different answer.
    pub backtrack_limit: usize,
}

impl Default for RegexOptions {
    fn default() -> Self {
        Self {
            icase: false,
            no_mod_s: false,
            backtrack_limit: DEFAULT_BACKTRACK_LIMIT,
        }
    }
}

/// The groups of one match, as `boost::smatch` holds them.
///
/// Offsets are byte positions in the haystack that was searched.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Captures {
    groups: Vec<Option<Range<usize>>>,
    names: Arc<[(String, usize)]>,
}

impl Captures {
    /// Number of groups including group 0, which is `smatch::size()`: one more
    /// than [`BoostRegex::mark_count`].
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Always `false`: a match has at least group 0.
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    /// The whole match, group 0.
    pub fn range(&self) -> Range<usize> {
        self.groups.first().and_then(Clone::clone).unwrap_or(0..0)
    }

    /// Byte range of group `index`, or `None` when the group did not take part
    /// in the match or does not exist.
    ///
    /// Boost reports a group that did not take part as an unmatched
    /// `sub_match` positioned at the end of the searched range; its text is
    /// empty either way.
    pub fn get(&self, index: usize) -> Option<Range<usize>> {
        self.groups.get(index).cloned().flatten()
    }

    /// Byte range of the named group, following
    /// `match_results::named_subexpression`: when several groups share the
    /// name, the first of them, in pattern order, that took part in the match.
    ///
    /// `None` when no group of that name took part (Boost then returns its null
    /// `sub_match`, positioned at the end of the match, whose text is empty).
    pub fn name(&self, name: &str) -> Option<Range<usize>> {
        self.names
            .iter()
            .filter(|(group_name, _)| group_name == name)
            .find_map(|&(_, index)| self.get(index))
    }
}

/// One token of [`BoostRegex::tokens`], a `boost::ssub_match`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubMatch {
    /// Byte range in the haystack, positioned where Boost positions the
    /// `sub_match` even when it did not match (see [`BoostRegex::tokens`]).
    pub range: Range<usize>,
    /// `sub_match::matched`. A text-between-matches token is matched exactly
    /// when it is non-empty, as in Boost.
    pub matched: bool,
}

impl SubMatch {
    /// The token's text in `haystack`, empty when it did not match.
    ///
    /// A range outside `haystack` (a token applied to a different buffer)
    /// yields an empty slice rather than a panic.
    pub fn as_bytes<'h>(&self, haystack: &'h [u8]) -> &'h [u8] {
        if self.matched {
            haystack.get(self.range.clone()).unwrap_or(&[])
        } else {
            &[]
        }
    }
}

/// A haystack as the engine sees it; see the module documentation.
#[derive(Debug)]
enum Haystack<'h> {
    /// Pure ASCII, already valid UTF-8; offsets are shared.
    Ascii(&'h str),
    /// Bytes above `0x7F` transcoded; `offsets[i]` is where original byte `i`
    /// starts, and the last entry is the transcoded length.
    Transcoded { text: String, offsets: Vec<usize> },
}

impl<'h> Haystack<'h> {
    fn new(bytes: &'h [u8]) -> Self {
        if bytes.is_ascii() {
            if let Ok(text) = std::str::from_utf8(bytes) {
                return Self::Ascii(text);
            }
        }
        let mut text = String::with_capacity(bytes.len().saturating_mul(3));
        let mut offsets = Vec::with_capacity(bytes.len().saturating_add(1));
        for &byte in bytes {
            offsets.push(text.len());
            let character = if byte.is_ascii() {
                char::from(byte)
            } else {
                char::from_u32(TRANSCODE_BASE + u32::from(byte))
                    .unwrap_or(char::REPLACEMENT_CHARACTER)
            };
            text.push(character);
        }
        offsets.push(text.len());
        Self::Transcoded { text, offsets }
    }

    fn text(&self) -> &str {
        match self {
            Self::Ascii(text) => text,
            Self::Transcoded { text, .. } => text,
        }
    }

    /// Length in original bytes.
    fn len(&self) -> usize {
        match self {
            Self::Ascii(text) => text.len(),
            Self::Transcoded { offsets, .. } => offsets.len() - 1,
        }
    }

    fn to_engine(&self, original: usize) -> usize {
        match self {
            Self::Ascii(_) => original,
            Self::Transcoded { text, offsets } => {
                offsets.get(original).copied().unwrap_or(text.len())
            }
        }
    }

    fn to_original(&self, engine: usize) -> usize {
        match self {
            Self::Ascii(_) => engine,
            Self::Transcoded { offsets, .. } => match offsets.binary_search(&engine) {
                Ok(index) | Err(index) => index,
            },
        }
    }
}

/// The two engine spellings of one translation.
///
/// `fancy-regex` counts backtracking steps against its budget, but not the
/// work it does between them: an automaton it runs anchored at a position can
/// scan to the end of the haystack, and a counted repeat pushes no backtracking
/// branch until it reaches its minimum. From every start position of a search
/// that work is repeated, so a short expression could do work quadratic in the
/// haystack without spending the budget. The counted spelling closes both
/// gaps: [`FORCE_BACKTRACKING_ROOT`] at the end of the expression and
/// [`FORCE_BACKTRACKING`] at the end of every positive lookahead make the engine
/// run the sub-expressions before them on its backtracking machine, where every
/// loop iteration pushes a branch, and [`COUNT_ITERATION`] makes every iteration
/// of a counted repeat push one.
#[derive(Clone, Debug)]
struct EngineTexts {
    /// The plain spelling, when the expression needs no backtracking.
    plain: Option<String>,
    /// Automaton size of the expression (see [`Translator::last_atoms`]).
    atoms: usize,
    /// The counted spelling, for every search the plain one is not used for,
    /// and for the token iterator's non-empty retry, which the engine always
    /// runs on its backtracking machine.
    counted: String,
}

impl EngineTexts {
    /// Translate `pattern` twice, plainly and in counted form, and once more
    /// with shielded repeats when the engine's optimizer would rewrite either
    /// (see [`optimizer_rewrites`]).
    fn new(pattern: &str, options: RegexOptions) -> Result<(Self, Translation)> {
        let translate = |shield_repeats: bool| -> Result<(Translation, String)> {
            let mut plain = Translator::new(pattern, options, false);
            plain.shield_repeats = shield_repeats;
            let plain = plain.run()?;
            let mut counted = Translator::new(pattern, options, true);
            counted.shield_repeats = shield_repeats;
            let counted = format!("(?:{}){FORCE_BACKTRACKING_ROOT}", counted.run()?.pattern);
            if plain.pattern.len().max(counted.len()) > MAX_TRANSLATED_BYTES {
                return Err(Error::Unsupported(format!(
                    "regular expression {} uses more than MAX_TRANSLATED_BYTES bytes of engine \
                     pattern at byte {}, which the Boost.Regex facade does not translate",
                    quote(pattern),
                    pattern.len()
                )));
            }
            Ok((plain, counted))
        };
        // The plain spelling reaches the engine only when it needs no
        // backtracking; the counted one always does.
        let rewritten = |(plain, counted): &(Translation, String)| {
            (!needs_backtracking(&plain.pattern) && optimizer_rewrites(&plain.pattern))
                || optimizer_rewrites(counted)
        };
        let mut translations = translate(false)?;
        if rewritten(&translations) {
            translations = translate(true)?;
            if rewritten(&translations) {
                return Err(Error::Unsupported(format!(
                    "regular expression {} uses a shape the engine's optimizer rewrites, which the \
                     Boost.Regex facade does not translate",
                    quote(pattern)
                )));
            }
        }
        let (plain, counted) = translations;
        let texts = Self {
            plain: (!needs_backtracking(&plain.pattern)).then(|| plain.pattern.clone()),
            atoms: plain.atoms,
            counted,
        };
        Ok((texts, plain))
    }

    /// The spelling for `regex_search` and, anchored at both ends,
    /// `regex_match` over a haystack that scales the budget by `scale`: the
    /// plain one, which the engine hands whole to an automaton, when the
    /// expression needs no backtracking and has at most
    /// [`MAX_AUTOMATON_ATOMS`]` / scale` atoms, and the counted one otherwise.
    fn search(&self, scale: usize) -> &str {
        match &self.plain {
            Some(plain) if self.atoms <= MAX_AUTOMATON_ATOMS / scale.max(1) => plain,
            _ => &self.counted,
        }
    }
}

/// Whether `fancy-regex`'s optimizer would rewrite the tree of `translated`
/// other than by moving a trailing lookahead out of group 0.
///
/// The check parses the text as the engine does (`RegexBuilder` passes only
/// the Unicode flag), appends an empty node to the root so that the
/// trailing-lookahead rewrite (`optimize_trailing_lookahead`) cannot apply, and
/// runs the engine's own `optimize`. A text the engine cannot parse is left to
/// the compiler, which reports it. The full-match program `\A(?:X)\z` needs no
/// check of its own: its tree holds the tree of `X` as one child between two
/// assertions, and every rewrite looks only at a repeat and its children or at
/// three neighbouring repeats.
///
/// Before it compiles an expression, `fancy-regex` 0.19.2 rewrites nested and
/// adjacent repeats (`optimize_nested_repeats`,
/// `optimize_ambiguous_concat_repeats`), and several of those rewrites change
/// what the expression matches:
///
/// - `X+Y?X+` and `X+Y*X+` become `X+(?:Y{1}X+)?`, which also matches a single
///   `X` (`\w+b?\w+` matches `a`);
/// - `(X+)+` becomes `(X+)` even when a backreference reads the group, so the
///   group can no longer be a later iteration (`(.+)+\1+` on `abb` matches
///   `1..3` instead of Boost's `0..3`);
/// - in an expression without backreferences, `(X)*` whose body is one
///   unbounded repeat becomes `(X)?`, which captures differently when that
///   repeat is lazy (`(\w+?)*?(?<=b)` anchored on `ab` captures `0..2` instead
///   of Boost's `1..2`).
///
/// Every rewrite needs a repeat, unbounded or optional, whose parent or
/// neighbour is another repeat or a group holding exactly one. Shielded
/// repeats ([`Translator::shield_repeats`]) are alternations instead, which no
/// rewrite matches, without a change in what matches or what is captured; each
/// costs one backtracking branch per pass. The facade translates with shielded repeats
/// only when this check finds a rewrite in the ordinary translation, and refuses
/// the pattern if it still finds one.
fn optimizer_rewrites(translated: &str) -> bool {
    let Ok(mut tree) =
        fancy_regex::Expr::parse_tree_with_flags(translated, fancy_regex::internal::FLAG_UNICODE)
    else {
        return false;
    };
    let padded = fancy_regex::Expr::Concat(vec![tree.expr.clone(), fancy_regex::Expr::Empty]);
    tree.expr = padded.clone();
    fancy_regex::internal::optimize(&mut tree);
    tree.expr != padded
}

/// Whether `fancy-regex` runs a translation, or part of it, on its
/// backtracking machine: the translation contains a lookaround, an atomic
/// group, a backreference or a word boundary. Every one of those is spelled
/// with a token below, and no other part of a translation contains one
/// (literals and set members are `\xHH` or alphanumeric).
fn needs_backtracking(translated: &str) -> bool {
    ["(?=", "(?!", "(?<", "(?>", r"\b", r"\B", r"\<", r"\>"]
        .iter()
        .any(|token| translated.contains(token))
        || translated
            .as_bytes()
            .windows(2)
            .any(|pair| pair[0] == b'\\' && (b'1'..=b'9').contains(&pair[1]))
}

/// The engine programs for one backtracking budget: one for `regex_search`, one
/// anchored at both ends for `regex_match`, and one that refuses empty matches
/// for the token iterator's `match_not_initial_null` retry (`None` when the
/// expression can only match the empty string).
#[derive(Clone, Debug)]
struct Programs {
    search: fancy_regex::Regex,
    full: fancy_regex::Regex,
    not_empty: Option<fancy_regex::Regex>,
}

impl Programs {
    /// Compile `texts` for a haystack that scales the budget by `scale`, with
    /// `backtrack_limit`, and check that every program numbers exactly the
    /// pattern's groups.
    fn compile(
        pattern: &str,
        texts: &EngineTexts,
        scale: usize,
        backtrack_limit: usize,
        mark_count: usize,
    ) -> Result<Self> {
        let spelling = texts.search(scale);
        let search = compile(pattern, spelling, backtrack_limit)?;
        let full = compile(pattern, &format!(r"\A(?:{spelling})\z"), backtrack_limit)?;
        let not_empty = match build(&texts.counted, backtrack_limit, true) {
            Ok(regex) => Some(regex),
            Err(fancy_regex::Error::CompileError(error))
                if matches!(*error, fancy_regex::CompileError::PatternCanNeverMatch) =>
            {
                None
            }
            Err(error) => return Err(engine_error(pattern, &error)),
        };
        let miscounted = [Some(&search), Some(&full), not_empty.as_ref()]
            .into_iter()
            .flatten()
            .find(|regex| regex.captures_len() != mark_count + 1);
        if let Some(regex) = miscounted {
            return Err(Error::Unsupported(format!(
                "regular expression {} compiled to {} groups instead of {}",
                quote(pattern),
                regex.captures_len(),
                mark_count + 1
            )));
        }
        Ok(Self {
            search,
            full,
            not_empty,
        })
    }
}

/// A compiled regular expression with `boost::regex` semantics (perl syntax).
///
/// Construction translates the Boost pattern once and compiles the engine
/// programs for the configured backtracking budget: one for `regex_search`, one
/// anchored at both ends for `regex_match`, and one that refuses empty matches
/// for the token iterator's `match_not_initial_null` retry. A haystack longer
/// than [`BACKTRACK_LIMIT_BYTES`] is searched with the same programs compiled
/// for a scaled budget, the first time such a haystack needs them; they are
/// kept for later searches.
#[derive(Clone, Debug)]
pub struct BoostRegex {
    pattern: String,
    options: RegexOptions,
    texts: EngineTexts,
    programs: Programs,
    /// Programs for the budget scaled by 4, 16 and 64, compiled on first use.
    /// `None` if the engine refused one. The unscaled programs are then used
    /// when they have the same spelling, so the search fails earlier rather
    /// than differently, and the search fails otherwise.
    scaled: [OnceLock<Option<Programs>>; 3],
    mark_count: usize,
    names: Arc<[(String, usize)]>,
}

impl BoostRegex {
    /// Compile `pattern` with default options, as `boost::regex(pattern)` does.
    ///
    /// # Errors
    ///
    /// See [`BoostRegex::with_options`].
    pub fn new(pattern: &str) -> Result<Self> {
        Self::with_options(pattern, RegexOptions::default())
    }

    /// Compile `pattern` with `options`, as
    /// `boost::regex(pattern, boost::regex::perl | flags)` does.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when Boost would throw `regex_error` for the
    ///   pattern: unbalanced parentheses or brackets, nothing to repeat, a
    ///   repeat of a zero-width assertion, a backreference to a group that does
    ///   not exist, a reversed range, an unknown POSIX class, a lookbehind
    ///   without one fixed width, an empty lookahead, and similar. Also when the
    ///   pattern is longer than [`MAX_PATTERN_BYTES`] or the backtrack limit is
    ///   zero.
    /// - [`Error::Unsupported`] for valid Boost syntax that the facade does not
    ///   translate (possessive quantifiers, `\Q...\E`, recursion, conditionals,
    ///   `\K`, `\G`, `\Z`, non-ASCII pattern bytes, and the other constructs
    ///   listed in `docs/BOOST_REGEX_SUPPORT.md`), for nesting deeper than
    ///   [`MAX_GROUP_DEPTH`], and when the translated pattern exceeds an engine
    ///   limit.
    pub fn with_options(pattern: &str, options: RegexOptions) -> Result<Self> {
        if pattern.len() > MAX_PATTERN_BYTES {
            return Err(Error::InvalidValue(format!(
                "regular expression of {} bytes exceeds the {MAX_PATTERN_BYTES}-byte limit",
                pattern.len()
            )));
        }
        if options.backtrack_limit == 0 {
            return Err(Error::InvalidValue(
                "regular expression backtrack limit must be positive".to_string(),
            ));
        }
        let (texts, translation) = EngineTexts::new(pattern, options)?;
        let programs = Programs::compile(
            pattern,
            &texts,
            1,
            options.backtrack_limit,
            translation.mark_count,
        )?;
        Ok(Self {
            pattern: pattern.to_string(),
            options,
            texts,
            programs,
            scaled: std::array::from_fn(|_| OnceLock::new()),
            mark_count: translation.mark_count,
            names: translation.names.into(),
        })
    }

    /// The programs whose backtracking budget and spelling cover a haystack of
    /// `length` bytes, with that budget.
    fn programs(&self, length: usize) -> Result<(&Programs, usize)> {
        let needed = length.saturating_sub(1) / BACKTRACK_LIMIT_BYTES + 1;
        let tier = BACKTRACK_SCALES
            .iter()
            .position(|&scale| scale >= needed)
            .unwrap_or(BACKTRACK_SCALES.len() - 1);
        let unscaled = (&self.programs, self.options.backtrack_limit);
        let (Some(cell), Some(&scale)) = (
            tier.checked_sub(1).and_then(|index| self.scaled.get(index)),
            BACKTRACK_SCALES.get(tier),
        ) else {
            return Ok(unscaled);
        };
        let limit = self.options.backtrack_limit.saturating_mul(scale);
        let scaled = cell.get_or_init(|| {
            Programs::compile(&self.pattern, &self.texts, scale, limit, self.mark_count).ok()
        });
        match scaled {
            Some(programs) => Ok((programs, limit)),
            None if self.texts.search(scale) == self.texts.search(1) => Ok(unscaled),
            None => Err(Error::InvalidValue(format!(
                "regular expression {} has no bounded engine program for a haystack of {length} \
                 bytes",
                quote(&self.pattern)
            ))),
        }
    }

    /// The pattern as given, `basic_regex::str()`.
    pub fn as_str(&self) -> &str {
        &self.pattern
    }

    /// The options the expression was compiled with.
    pub fn options(&self) -> RegexOptions {
        self.options
    }

    /// Number of capturing groups, `basic_regex::mark_count()` (group 0 not
    /// counted).
    pub fn mark_count(&self) -> usize {
        self.mark_count
    }

    /// Named groups as `(name, group index)` in pattern order. A name used by
    /// several groups appears once per group.
    pub fn capture_names(&self) -> impl Iterator<Item = (&str, usize)> + '_ {
        self.names
            .iter()
            .map(|(name, index)| (name.as_str(), *index))
    }

    /// `boost::regex_search(haystack, match, regex)`: the leftmost match, with
    /// Boost's alternation priority.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`], never a different answer, when the search hits
    /// one of the engine's two runtime limits; the message names which:
    ///
    /// - the backtracking budget of [`RegexOptions::backtrack_limit`], scaled
    ///   for long haystacks, where Boost would throw `regex_error` with
    ///   `error_complexity`;
    /// - the 1,000,000 entries of `fancy-regex`'s fixed backtracking stack. A
    ///   greedy repeat followed by a construct that needs backtracking (a line
    ///   anchor, a lookaround) pushes one entry per repetition, so
    ///   `=(?<SCAN>\d+)$` fails on `=` followed by a million digits, where
    ///   Boost's growable stack answers.
    ///
    /// Both limits are counted differently from Boost's, so the inputs that hit
    /// them differ; `docs/BOOST_REGEX_SUPPORT.md` records measured haystack
    /// lengths.
    pub fn search(&self, haystack: &[u8]) -> Result<Option<Captures>> {
        self.search_from(&Haystack::new(haystack), 0, false)
    }

    /// `boost::regex_search(haystack, regex)` without match results.
    ///
    /// # Errors
    ///
    /// As [`BoostRegex::search`].
    pub fn is_search_match(&self, haystack: &[u8]) -> Result<bool> {
        let haystack = Haystack::new(haystack);
        let (programs, limit) = self.programs(haystack.len())?;
        programs
            .search
            .is_match(haystack.text())
            .map_err(|error| self.runtime_error(&error, limit))
    }

    /// `boost::regex_match(haystack, match, regex)`: a match of the whole
    /// haystack, backtracking into alternatives until one spans it.
    ///
    /// # Errors
    ///
    /// As [`BoostRegex::search`].
    pub fn full_match(&self, haystack: &[u8]) -> Result<Option<Captures>> {
        let haystack = Haystack::new(haystack);
        let (programs, limit) = self.programs(haystack.len())?;
        self.captures(
            &programs.full,
            &haystack,
            RegexInput::new(haystack.text()),
            limit,
        )
    }

    /// `boost::regex_match(haystack, regex)` without match results.
    ///
    /// # Errors
    ///
    /// As [`BoostRegex::search`].
    pub fn is_full_match(&self, haystack: &[u8]) -> Result<bool> {
        let haystack = Haystack::new(haystack);
        let (programs, limit) = self.programs(haystack.len())?;
        programs
            .full
            .is_match(haystack.text())
            .map_err(|error| self.runtime_error(&error, limit))
    }

    /// `boost::sregex_token_iterator(begin, end, regex, submatches)`.
    ///
    /// For every match the iterator yields one token per entry of `submatches`,
    /// in order, where an entry means what `match_results::operator[]` and the
    /// token iterator make of it:
    ///
    /// - `-1`: the text between the previous match (or the haystack start) and
    ///   this match;
    /// - `-2`: the text from the end of this match to the end of the haystack;
    /// - `0..=mark_count`: that group, or an unmatched token at the end of the
    ///   haystack when the group did not take part;
    /// - any other value: Boost's null `sub_match`, unmatched and empty at the
    ///   end of the match.
    ///
    /// After the last match, when the first entry is `-1`, the rest of the
    /// haystack follows if it is non-empty. After an empty match the next search
    /// refuses a match that ends where the previous one did
    /// (`match_not_initial_null`), exactly as Boost does, so an alternative that
    /// matches text at that position is still found.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `submatches` is empty (Boost indexes past
    /// the end of its vector). Items are [`Error::InvalidValue`] when a search
    /// hits a runtime limit, as for [`BoostRegex::search`]; each search between
    /// two tokens has its own budget, scaled by the whole haystack's length, and
    /// the iterator ends after an error. The work bounds therefore apply per
    /// search, not per iterator: a whole iteration may do the bound of one
    /// search once for every match (see [`MAX_AUTOMATON_ATOMS`]).
    pub fn tokens<'r, 'h>(
        &'r self,
        haystack: &'h [u8],
        submatches: &[i32],
    ) -> Result<Tokens<'r, 'h>> {
        if submatches.is_empty() {
            return Err(Error::InvalidValue(
                "token iterator needs at least one sub-expression index".to_string(),
            ));
        }
        Ok(Tokens {
            regex: self,
            haystack: Haystack::new(haystack),
            submatches: submatches.to_vec(),
            state: TokenState::Start,
        })
    }

    /// One `regex_search` from byte `start`, with `match_not_initial_null` when
    /// `not_initial_null` is set: a match that ends at `start` is refused and
    /// the search backtracks, so a non-empty match starting at `start` is
    /// still preferred to anything further right.
    fn search_from(
        &self,
        haystack: &Haystack<'_>,
        start: usize,
        not_initial_null: bool,
    ) -> Result<Option<Captures>> {
        let text = haystack.text();
        let (programs, limit) = self.programs(haystack.len())?;
        if !not_initial_null {
            let input = RegexInput::new(text).from_pos(haystack.to_engine(start));
            return self.captures(&programs.search, haystack, input, limit);
        }
        if let Some(not_empty) = &programs.not_empty {
            let input = RegexInput::new(text)
                .from_pos(haystack.to_engine(start))
                .anchored(true);
            if let Some(captures) = self.captures(not_empty, haystack, input, limit)? {
                return Ok(Some(captures));
            }
        }
        if start >= haystack.len() {
            return Ok(None);
        }
        let input = RegexInput::new(text).from_pos(haystack.to_engine(start + 1));
        self.captures(&programs.search, haystack, input, limit)
    }

    fn captures(
        &self,
        regex: &fancy_regex::Regex,
        haystack: &Haystack<'_>,
        input: RegexInput<'_, str>,
        limit: usize,
    ) -> Result<Option<Captures>> {
        match regex.captures_input(input) {
            Ok(Some(found)) => Ok(Some(Captures {
                groups: (0..=self.mark_count)
                    .map(|index| {
                        found.get(index).map(|group| {
                            haystack.to_original(group.start())..haystack.to_original(group.end())
                        })
                    })
                    .collect(),
                names: Arc::clone(&self.names),
            })),
            Ok(None) => Ok(None),
            Err(error) => Err(self.runtime_error(&error, limit)),
        }
    }

    /// The error for a search that stopped with `error` under a budget of
    /// `limit` steps.
    fn runtime_error(&self, error: &fancy_regex::Error, limit: usize) -> Error {
        match error {
            fancy_regex::Error::RuntimeError(fancy_regex::RuntimeError::BacktrackLimitExceeded) => {
                Error::InvalidValue(format!(
                    "regular expression {} exceeded the backtracking limit of {limit} steps",
                    quote(&self.pattern)
                ))
            }
            fancy_regex::Error::RuntimeError(fancy_regex::RuntimeError::StackOverflow) => {
                Error::InvalidValue(format!(
                    "regular expression {} exceeded the engine's backtracking stack of \
                     {ENGINE_STACK_ENTRIES} entries",
                    quote(&self.pattern)
                ))
            }
            other => Error::InvalidValue(format!(
                "regular expression {} failed while matching: {other}",
                quote(&self.pattern)
            )),
        }
    }
}

/// Iterator returned by [`BoostRegex::tokens`].
#[derive(Debug)]
pub struct Tokens<'r, 'h> {
    regex: &'r BoostRegex,
    haystack: Haystack<'h>,
    submatches: Vec<i32>,
    state: TokenState,
}

#[derive(Debug)]
enum TokenState {
    Start,
    Match {
        captures: Captures,
        search_start: usize,
        next: usize,
    },
    Done,
}

impl Tokens<'_, '_> {
    fn token(&self, captures: &Captures, search_start: usize, index: i32) -> SubMatch {
        let whole = captures.range();
        let length = self.haystack.len();
        match index {
            -1 => SubMatch {
                matched: search_start != whole.start,
                range: search_start..whole.start,
            },
            -2 => SubMatch {
                matched: whole.end != length,
                range: whole.end..length,
            },
            _ => match usize::try_from(index)
                .ok()
                .filter(|&group| group <= self.regex.mark_count)
            {
                Some(group) => match captures.get(group) {
                    Some(range) => SubMatch {
                        range,
                        matched: true,
                    },
                    None => SubMatch {
                        range: length..length,
                        matched: false,
                    },
                },
                None => SubMatch {
                    range: whole.end..whole.end,
                    matched: false,
                },
            },
        }
    }

    fn found(&mut self, captures: Captures, search_start: usize) -> SubMatch {
        let token = self.token(&captures, search_start, self.submatches[0]);
        self.state = TokenState::Match {
            captures,
            search_start,
            next: 1,
        };
        token
    }
}

impl Iterator for Tokens<'_, '_> {
    type Item = Result<SubMatch>;

    fn next(&mut self) -> Option<Self::Item> {
        let length = self.haystack.len();
        match std::mem::replace(&mut self.state, TokenState::Done) {
            TokenState::Done => None,
            TokenState::Start => match self.regex.search_from(&self.haystack, 0, false) {
                Err(error) => Some(Err(error)),
                Ok(Some(captures)) => Some(Ok(self.found(captures, 0))),
                Ok(None) => (self.submatches[0] == -1 && length > 0).then_some(Ok(SubMatch {
                    range: 0..length,
                    matched: true,
                })),
            },
            TokenState::Match {
                captures,
                search_start,
                next,
            } => {
                if let Some(&index) = self.submatches.get(next) {
                    let token = self.token(&captures, search_start, index);
                    self.state = TokenState::Match {
                        captures,
                        search_start,
                        next: next + 1,
                    };
                    return Some(Ok(token));
                }
                let whole = captures.range();
                match self
                    .regex
                    .search_from(&self.haystack, whole.end, whole.is_empty())
                {
                    Err(error) => Some(Err(error)),
                    Ok(Some(found)) => Some(Ok(self.found(found, whole.end))),
                    Ok(None) => {
                        (whole.end != length && self.submatches[0] == -1).then_some(Ok(SubMatch {
                            range: whole.end..length,
                            matched: true,
                        }))
                    }
                }
            }
        }
    }
}

impl FusedIterator for Tokens<'_, '_> {}

fn build(
    translated: &str,
    backtrack_limit: usize,
    not_empty: bool,
) -> std::result::Result<fancy_regex::Regex, fancy_regex::Error> {
    let mut builder = RegexBuilder::new(translated);
    builder
        .bytes_mode(BytesMode::Unicode)
        .unicode_mode(true)
        .backtrack_limit(backtrack_limit)
        .find_not_empty(not_empty);
    builder.build()
}

fn compile(pattern: &str, translated: &str, backtrack_limit: usize) -> Result<fancy_regex::Regex> {
    build(translated, backtrack_limit, false).map_err(|error| engine_error(pattern, &error))
}

fn engine_error(pattern: &str, error: &fancy_regex::Error) -> Error {
    Error::Unsupported(format!(
        "regular expression {} is valid Boost syntax but the engine refused its translation: {error}",
        quote(pattern)
    ))
}

/// The pattern for error messages, shortened so a 64 KiB pattern does not
/// become a 64 KiB message.
fn quote(pattern: &str) -> String {
    const SHOWN: usize = 120;
    if pattern.len() <= SHOWN {
        return format!("{pattern:?}");
    }
    let mut end = SHOWN;
    while !pattern.is_char_boundary(end) {
        end -= 1;
    }
    format!("{:?}...", &pattern[..end])
}

// ---------------------------------------------------------------------------
// Translation from Boost perl syntax to an explicit fancy-regex pattern
// ---------------------------------------------------------------------------

/// Start of line: buffer start, after `\n` or `\f`, or after `\r` unless a `\n`
/// follows (Boost `match_start_line`).
const START_OF_LINE: &str = r"(?:\A|(?<=[\x0A\x0C])|(?<=\x0D)(?!\x0A))";
/// End of line: buffer end, before `\r` or `\f`, or before `\n` unless a `\r`
/// precedes it (Boost `match_end_line`).
const END_OF_LINE: &str = r"(?:\z|(?=[\x0D\x0C])|(?<!\x0D)(?=\x0A))";
/// `.` under `(?-s)`: any byte except Boost's line separators.
const DOT_NO_SEPARATOR: &str = r"[^\x0A\x0C\x0D]";
/// A zero-width lookahead that always holds, at the end of every positive
/// lookahead body of the counted spelling.
///
/// `fancy-regex` compiles a concatenation that needs backtracking by handing
/// its longest trailing run of sub-expressions that do not to an automaton
/// (`compile_concat`), which runs anchored and can scan to the end of the
/// haystack. A run that ends in this assertion is empty, so everything before it
/// is compiled for the backtracking machine instead.
const FORCE_BACKTRACKING: &str = "(?=)";
/// The same for the whole expression, which needs the assertion twice.
///
/// Unless it refuses empty matches, `fancy-regex` first rewrites an expression
/// that ends in a positive lookahead, `X(?=Y)`, into an explicit group 0 around
/// `X` followed by `Y` (`optimize_trailing_lookahead`). One `(?=)` would become
/// an empty `Y`, and an `X` that needs no backtracking would then be handed to
/// an automaton whole. With two, the rewrite takes the second, the group keeps
/// the first, and the group needs backtracking. Should a crate update change
/// either rewrite, `bounded_spelling_runs_on_the_backtracking_machine` in the
/// module tests and `forced_backtracking_spends_the_budget` in
/// `tests/boost_regex.rs` fail.
const FORCE_BACKTRACKING_ROOT: &str = "(?=)(?=)";
/// An empty alternative beside one that never matches: every pass through it
/// pushes one backtracking branch, which the budget counts when it is undone.
const COUNT_ITERATION: &str = "(?:|(?!))";

/// Refusal of a repeat allowing more than one iteration of a group that can
/// match the empty string; see [`Translator::repeat_refusal`].
const NULLABLE_GROUP_REPEAT: &str = "a repeat of more than one iteration of a group that can match \
                                     the empty string (Boost ends a repeat after an empty iteration)";
/// Refusal of a backreference inside the group it refers to.
///
/// Boost saves a group's previous span when the group starts, sets only its
/// start, and restores the saved span only when the match backtracks past that
/// start (`perl_matcher::match_startmark`, `match_results::set_first`). A
/// backreference inside the group therefore compares against whatever an
/// abandoned attempt at the group left: the span of an alternative that matched
/// and then failed further on (`(a|\1b)x` matches all of `abx`), or, in a later
/// iteration, the range from the new start to the previous end, which is empty
/// when the iterations are adjacent (`(a|b\1)+` matches all of `ab`) and never
/// matches otherwise. `fancy-regex` treats the group as unset in the first case,
/// and in the second slices the haystack from the new start to the previous end,
/// which panics when the start is the larger (`(?:(a|\1b)c)+` on `acbc`).
const OPEN_GROUP_BACKREFERENCE: &str = "a backreference inside the group it refers to (Boost \
                                        compares the span an abandoned attempt at the group left)";
/// Refusal of a lazy repeat that Boost's leading-repeat optimization restarts
/// from; see [`Translator::quantify`].
const LEADING_LAZY_REPEAT: &str = "a lazy repeat with a finite maximum of a one-byte atom that \
                                   starts the expression (Boost's leading-repeat optimization skips \
                                   start positions)";
/// Refusal of a repeat of a group whose body holds a repeat or an alternation
/// compiled under the other case sensitivity.
///
/// Before matching, Boost gives every repeat and alternation state a map of the
/// bytes that can start each of its two ways on
/// (`basic_regex_creator::create_startmaps`), and the matcher consults it before
/// it enters an iteration, leaves a repeat or tries an alternative. The map of a
/// state is built by walking the states that follow it, and whenever that walk
/// reaches a repeat or alternation whose map is not built yet, which happens when
/// it loops back to the repeat of an enclosing group, it recurses into
/// `create_startmap`, which starts again from the case sensitivity of the state
/// whose map is being built (`bool l_icase = m_icase`) rather than from the one in
/// effect where the walk stands. A literal or set further on is then looked up
/// with the wrong case, and the map drops the bytes that start it:
/// `(?i:b.+)*C` on `bcC` matches `2..3` in Boost and `(?i:bc+)+C` does not
/// match, where trying every start position gives `0..3` for both.
const CASE_SWITCHED_REPEAT: &str = "a repeat of a group holding a repeat or alternation under other \
                                    case sensitivity (Boost builds its start map with the wrong case)";
/// Refusal of `\<` in an expression that switches case sensitivity.
///
/// Boost's start-map walk (see [`CASE_SWITCHED_REPEAT`]) recurses at `\<` and
/// `\>` to collect what follows and then removes the bytes the assertion rules
/// out, and the recursion starts from the case sensitivity of the state whose map
/// is being built, or of the options for the expression's own map. After `\<` a
/// letter is looked up with the wrong case: `(?i)\<A` does not match `A` in
/// Boost. (After `\>` only non-word bytes remain, which have no case.)
const WORD_START_CASE_SWITCH: &str = "\\< in an expression that switches case sensitivity (Boost \
                                      builds its start map with the wrong case)";
/// Refusal of `\<` and `\>` in an expression with a repeated group that holds a
/// repeat or an alternation.
///
/// When Boost's start-map walk reaches `\<` or `\>` it recurses and then removes
/// every word byte (or every non-word byte) from the whole map it is filling. When
/// the walk has looped back to an enclosing repeat, that map already holds the
/// bytes that start another iteration, which come before the assertion, and loses
/// them: `(?:\w\w+){2}\>` does not match `aaaa`, and `(\w+?)+\>` on `aaaa` captures
/// group 1 at `0..4` in Boost instead of `3..4`.
const WORD_BOUNDARY_AFTER_LOOP: &str = "\\< or \\> in an expression with a repeated group holding a \
                                        repeat or alternation (Boost's start map drops the bytes that \
                                        start another iteration)";
/// Refusal of an expression whose start maps Boost might not build; see
/// [`Translator::start_map_depth`].
const START_MAP_RECURSION: &str = "more line anchors, word-boundary assertions, alternations and \
                                   repeated groups than Boost's start-map recursion limit allows \
                                   (Boost throws error_complexity)";
/// Refusal of a lookbehind with more alternations than Boost's backstep
/// calculation stacks.
///
/// `basic_regex_creator::calculate_backstep` pushes every alternation it meets on
/// the way to the end of a lookbehind and gives up, rejecting the expression, when
/// it would push one more than `BOOST_REGEX_MAX_BLOCKS` (1,024) of them.
const LOOKBEHIND_ALTERNATIONS: &str = "a lookbehind with more than 1,024 alternations (Boost's \
                                       backstep calculation gives up)";
/// Refusal of a `\x` escape that `std::istream` reads differently from hex digits.
///
/// Boost reads the digits of `\xHH` (at most two bytes) and `\x{H...}` with
/// `cpp_regex_traits::toi`, which is `std::istream >> std::hex >> intmax_t`: it
/// skips white space and accepts a sign and a `0x` prefix, so `\x{+41}`,
/// `\x{ 41}`, `\x{0x41}` and `\x-0` are valid there, and `\x0x` is not.
const HEX_ESCAPE_STREAM_SYNTAX: &str = "a hexadecimal escape whose digits start with a sign, white \
                                        space or 0x (Boost reads them with std::istream)";

/// Refusal of a pattern whose translation could nest groups as deep as
/// [`ENGINE_MAX_NESTING`]; see [`MAX_GROUP_DEPTH`].
const TRANSLATION_TOO_DEEP: &str = "groups and repeats nested deeper than the engine's parser allows \
                                    once translated";

/// `fancy-regex`'s `MAX_RECURSION`: its parser refuses a group at this depth.
const ENGINE_MAX_NESTING: usize = 64;

/// `BOOST_REGEX_MAX_RECURSION_DEPTH`: the deepest recursion of Boost's
/// `create_startmap` before it throws `error_complexity`.
const BOOST_START_MAP_DEPTH: usize = 100;

/// `BOOST_REGEX_MAX_BLOCKS`: the alternations Boost's `calculate_backstep` stacks.
const BOOST_BACKSTEP_BLOCKS: usize = 1024;

/// A Boost character class, as `cpp_regex_traits<char>` defines it in the C
/// locale.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ClassKind {
    Alnum,
    Alpha,
    Blank,
    Cntrl,
    Digit,
    Graph,
    Horizontal,
    Lower,
    Print,
    Punct,
    Space,
    Upper,
    Vertical,
    Word,
    Xdigit,
}

impl ClassKind {
    /// Boost's `isctype` for an ASCII byte. No byte above `0x7F` belongs to
    /// any class in the C locale.
    fn contains(self, byte: u8) -> bool {
        match self {
            Self::Alnum => byte.is_ascii_alphanumeric(),
            Self::Alpha => byte.is_ascii_alphabetic(),
            // isspace && !is_separator: space, tab and vertical tab.
            Self::Blank => matches!(byte, b'\t' | 0x0B | b' '),
            Self::Cntrl => byte.is_ascii_control(),
            Self::Digit => byte.is_ascii_digit(),
            Self::Graph => byte.is_ascii_graphic(),
            // isspace && !is_separator && c != '\v'.
            Self::Horizontal => matches!(byte, b'\t' | b' '),
            Self::Lower => byte.is_ascii_lowercase(),
            Self::Print => byte.is_ascii_graphic() || byte == b' ',
            Self::Punct => byte.is_ascii_punctuation(),
            Self::Space => matches!(byte, b'\t'..=b'\r' | b' '),
            Self::Upper => byte.is_ascii_uppercase(),
            // is_separator || c == '\v'.
            Self::Vertical => matches!(byte, b'\n'..=b'\r'),
            Self::Word => byte.is_ascii_alphanumeric() || byte == b'_',
            Self::Xdigit => byte.is_ascii_hexdigit(),
        }
    }

    fn posix(name: &[u8]) -> Option<Self> {
        Some(match name {
            b"alnum" => Self::Alnum,
            b"alpha" => Self::Alpha,
            b"blank" => Self::Blank,
            b"cntrl" => Self::Cntrl,
            b"digit" => Self::Digit,
            b"graph" => Self::Graph,
            b"lower" => Self::Lower,
            b"print" => Self::Print,
            b"punct" => Self::Punct,
            b"space" => Self::Space,
            b"upper" => Self::Upper,
            b"word" => Self::Word,
            b"xdigit" => Self::Xdigit,
            _ => return None,
        })
    }

    /// The class of a perl escape letter: `\d \s \w \h \v` and their negations.
    fn escape(letter: u8) -> Option<(Self, bool)> {
        Some(match letter {
            b'd' => (Self::Digit, false),
            b'D' => (Self::Digit, true),
            b's' => (Self::Space, false),
            b'S' => (Self::Space, true),
            b'w' => (Self::Word, false),
            b'W' => (Self::Word, true),
            b'h' => (Self::Horizontal, false),
            b'H' => (Self::Horizontal, true),
            b'v' => (Self::Vertical, false),
            b'V' => (Self::Vertical, true),
            _ => return None,
        })
    }
}

/// One member of a bracket expression, as Boost's `basic_char_set` records it.
#[derive(Clone, Copy, Debug)]
enum ClassItem {
    Byte(u8),
    /// A range with its endpoints as written; Boost translates them when it
    /// builds the set.
    Range(u8, u8),
    /// A class, `true` when negated (`\S`, `[:^alpha:]`).
    Class(ClassKind, bool),
}

/// A set that matches no character of a haystack: it excludes every ASCII byte
/// and every transcoded byte above `0x7F` (`U+E080` to `U+E0FF`).
const NO_BYTE: &str = r"[^\x00-\x7F\x{E080}-\x{E0FF}]";

/// Boost `basic_regex_creator::append_set` for `char` in the C locale: which
/// ASCII bytes a bracket expression or class escape matches, and whether the
/// bytes above `0x7F` do.
///
/// Boost fills a map indexed by byte. A single marks every byte that
/// translates to the same byte as it; a range marks the bytes between its
/// translated endpoints; the classes mark their union; the negated classes mark
/// the complement of their union, taken once for all of them; under `icase` a
/// class union that names `lower` or `upper` is widened to `alpha`. `[^` then
/// inverts the map, and matching looks up each input byte after the same
/// translation, which is ASCII lower-casing under `icase` and the identity
/// otherwise. A byte above `0x7F` has no case and belongs to no class, so it is
/// in the set exactly when a negated class is present, inverted by `[^`.
fn class_set(icase: bool, negated: bool, items: &[ClassItem]) -> ([bool; 128], bool) {
    let translate = |byte: u8| {
        if icase {
            byte.to_ascii_lowercase()
        } else {
            byte
        }
    };
    let union = |kinds: &[ClassKind], byte: u8| {
        let widened = icase
            && kinds
                .iter()
                .any(|kind| matches!(kind, ClassKind::Lower | ClassKind::Upper));
        kinds.iter().any(|kind| kind.contains(byte)) || (widened && byte.is_ascii_alphabetic())
    };
    let mut classes = Vec::new();
    let mut negated_classes = Vec::new();
    let mut map = [false; 128];
    for item in items {
        match *item {
            ClassItem::Byte(single) => {
                for (byte, member) in (0u8..).zip(map.iter_mut()) {
                    *member |= translate(byte) == translate(single);
                }
            }
            ClassItem::Range(first, last) => {
                let range = translate(first)..=translate(last);
                for (byte, member) in (0u8..).zip(map.iter_mut()) {
                    *member |= range.contains(&byte);
                }
            }
            ClassItem::Class(kind, false) => classes.push(kind),
            ClassItem::Class(kind, true) => negated_classes.push(kind),
        }
    }
    let any_negated = !negated_classes.is_empty();
    for (byte, member) in (0u8..).zip(map.iter_mut()) {
        *member |= union(&classes, byte) || (any_negated && !union(&negated_classes, byte));
    }
    let mut accepted = [false; 128];
    for (byte, member) in (0u8..).zip(accepted.iter_mut()) {
        *member = map[usize::from(translate(byte))] != negated;
    }
    (accepted, any_negated != negated)
}

/// The engine spelling of a set from [`class_set`], one character wide. A
/// haystack holds only ASCII and transcoded high bytes, so when the high bytes
/// are in the set, a negated class naming the excluded ASCII bytes matches
/// exactly the set.
fn class_text(accepted: &[bool; 128], high: bool) -> String {
    if !high && !accepted.contains(&true) {
        return NO_BYTE.to_string();
    }
    if high && !accepted.contains(&false) {
        return "(?s:.)".to_string();
    }
    let mut runs: Vec<(u8, u8)> = Vec::new();
    for (byte, &member) in (0u8..).zip(accepted) {
        if member == high {
            continue;
        }
        match runs.last_mut() {
            Some((_, last)) if *last + 1 == byte => *last = byte,
            _ => runs.push((byte, byte)),
        }
    }
    let mut text = String::from(if high { "[^" } else { "[" });
    for (first, last) in runs {
        push_set_byte(&mut text, first);
        if last > first + 1 {
            text.push('-');
        }
        if last > first {
            push_set_byte(&mut text, last);
        }
    }
    text.push(']');
    text
}

/// One byte of a bracket expression: letters and digits as themselves,
/// everything else as `\xHH`.
fn push_set_byte(text: &mut String, byte: u8) {
    if byte.is_ascii_alphanumeric() {
        text.push(char::from(byte));
    } else {
        let _ = write!(text, r"\x{byte:02X}");
    }
}

#[derive(Clone, Copy, Debug)]
struct Flags {
    icase: bool,
    dotall: bool,
    multiline: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GroupKind {
    Root,
    Capture,
    NonCapture,
    Atomic,
    LookAhead,
    NegativeLookAhead,
    LookBehind,
    NegativeLookBehind,
}

impl GroupKind {
    fn is_lookaround(self) -> bool {
        matches!(
            self,
            Self::LookAhead | Self::NegativeLookAhead | Self::LookBehind | Self::NegativeLookBehind
        )
    }

    /// Boost runs these as independent sub-matches whose captures it does not
    /// restore when the surrounding match later backtracks past them.
    fn is_independent(self) -> bool {
        self.is_lookaround() || self == Self::Atomic
    }
}

/// What a closed group looks like to a quantifier that follows it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct GroupShape {
    start: usize,
    kind: GroupKind,
    body: Body,
    /// The group is, or contains, a capturing group.
    captures: bool,
    /// The body can match the empty string.
    nullable: bool,
    /// The body asserts positions (`^`, `$`, `\b`, ...) and consumes nothing.
    only_asserts: bool,
    /// The body holds the leading lazy repeat (see [`Translator::leading_repeat`]).
    holds_leading_repeat: bool,
    /// [`Frame::state_icase`] of the body.
    state_icase: [bool; 2],
    /// Line anchors and word-boundary assertions in the body outside every
    /// repeated group of it (see [`StartMapSites`]).
    free_anchors: usize,
    /// Alternations in the body outside every repeated group of it.
    free_alternations: usize,
}

impl GroupShape {
    /// A lookaround, or a non-capturing group around exactly one: the engine
    /// refuses to repeat it, so the translation tries it at most once instead.
    fn is_zero_width(self) -> bool {
        self.kind.is_lookaround()
            || (self.kind == GroupKind::NonCapture && self.body == Body::LookAround)
    }

    /// The translation hands a quantifier on this group to the engine, rather
    /// than dropping an empty group or trying a zero-width one at most once.
    fn is_repeated_by_engine(self) -> bool {
        !(self.is_zero_width() || (self.kind == GroupKind::NonCapture && self.body == Body::Empty))
    }
}

/// What a quantifier would apply to, following Boost's `m_last_state`.
#[derive(Clone, Copy, Debug)]
enum Last {
    /// Pattern start, `(` or `|`: Boost's "nothing to repeat".
    Nothing,
    /// A zero-width assertion or a repeat: Boost's `error_badrepeat`.
    Fixed,
    /// A one-byte atom: literal, class or `.`.
    Byte,
    /// A backreference; `nullable` unless the group it names closed earlier and
    /// cannot match the empty string (a backreference to a group that did not
    /// take part fails, in Boost's perl mode as in the engine).
    Backref { nullable: bool },
    /// A `(?imsx)` group, which Boost compiles to an empty group;
    /// `case_change` when it switches case sensitivity.
    EmptyGroup { case_change: bool },
    /// A closed group.
    Group(GroupShape),
}

/// How `fancy-regex` sees a non-capturing group's body, which decides whether
/// it accepts a quantifier on it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Body {
    /// Nothing but comments and flag groups: the engine's `Empty`.
    Empty,
    /// Exactly one lookaround: the engine's `LookAround`.
    LookAround,
    /// Anything else.
    Other,
}

/// Shortest match of a sub-expression, in bytes, counted for two purposes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct MinLength {
    /// For the repeat rules: a backreference counts zero, so a group is
    /// nullable whenever a backreference could make it so.
    boost: usize,
    /// As `fancy-regex`'s analyzer computes it before multiplying it by a
    /// repeat's minimum: a backreference counts the shortest match of its group
    /// when that group closed earlier, and zero otherwise.
    engine: usize,
}

impl MinLength {
    const BYTE: Self = Self {
        boost: 1,
        engine: 1,
    };

    fn followed_by(self, other: Self) -> Self {
        Self {
            boost: self.boost.saturating_add(other.boost),
            engine: self.engine.saturating_add(other.engine),
        }
    }

    fn shorter(self, other: Self) -> Self {
        Self {
            boost: self.boost.min(other.boost),
            engine: self.engine.min(other.engine),
        }
    }

    fn without(self, other: Self) -> Self {
        Self {
            boost: self.boost.saturating_sub(other.boost),
            engine: self.engine.saturating_sub(other.engine),
        }
    }

    fn repeated(self, count: usize) -> Self {
        Self {
            boost: self.boost.saturating_mul(count),
            engine: self.engine.saturating_mul(count),
        }
    }
}

/// The states of a sub-expression at which Boost's `create_startmap` recurses,
/// counted for [`Translator::start_map_depth`].
///
/// A site is free when no repeated group of the sub-expression encloses it and
/// looped otherwise; a repeat of the whole sub-expression turns its free sites
/// into looped ones.
#[derive(Clone, Copy, Debug, Default)]
struct StartMapSites {
    /// `$` as a line anchor, `\<` and `\>` outside every repeated group.
    free_anchors: usize,
    /// The same inside a repeated group.
    looped_anchors: usize,
    /// Alternation states (one per `|`) outside every repeated group.
    free_alternations: usize,
    /// The same inside a repeated group.
    looped_alternations: usize,
    /// Quantifiers applied to a group (Boost's `syntax_element_rep` around
    /// several states).
    repeated_groups: usize,
}

impl StartMapSites {
    fn add(&mut self, other: Self) {
        self.free_anchors = self.free_anchors.saturating_add(other.free_anchors);
        self.looped_anchors = self.looped_anchors.saturating_add(other.looped_anchors);
        self.free_alternations = self
            .free_alternations
            .saturating_add(other.free_alternations);
        self.looped_alternations = self
            .looped_alternations
            .saturating_add(other.looped_alternations);
        self.repeated_groups = self.repeated_groups.saturating_add(other.repeated_groups);
    }
}

#[derive(Debug)]
struct Frame {
    kind: GroupKind,
    saved_flags: Flags,
    start: usize,
    /// Number of the capturing group this frame is, or 0.
    group: usize,
    /// Width of the current alternative, `None` once Boost's `calculate_backstep`
    /// would give up.
    width: Option<usize>,
    /// Width of the first finished alternative: `None` before any, `Some(None)`
    /// once alternatives disagree or one has no width.
    alternative_width: Option<Option<usize>>,
    /// Shortest match of the current alternative.
    min_length: MinLength,
    /// Shortest match over the finished alternatives.
    alternative_min_length: Option<MinLength>,
    /// Whether Boost has appended any state inside the group (comments append none).
    has_states: bool,
    alternation: bool,
    items: usize,
    last_item_is_lookaround: bool,
    consumes: bool,
    asserts: bool,
    has_capture: bool,
    /// Automaton size of the current alternative (see [`Translator::last_atoms`]).
    atoms: usize,
    /// Automaton size of the finished alternatives.
    alternative_atoms: usize,
    /// The group was open when the leading lazy repeat was recorded.
    holds_leading_repeat: bool,
    /// For a lookaround, whether the states before it were all transparent to
    /// Boost's `probe_leading_repeat`, which skips a lookaround whole.
    leading_before: bool,
    /// Whether the group holds, at any depth, a Boost repeat or alternation
    /// state compiled case-sensitively (index 0) or case-insensitively (index 1)
    /// (see [`CASE_SWITCHED_REPEAT`]). A repeat state takes the case sensitivity in
    /// effect at its quantifier. The alternation state of the first `|` stands at
    /// the start of the group, before a scoped `(?i:` switches, and that of a later
    /// `|` at the start of the alternative before it, where the case sensitivity is
    /// the one in effect at the previous `|`.
    state_icase: [bool; 2],
    /// Case sensitivity of the alternation state the next `|` inserts.
    next_alternation_icase: bool,
    /// Where Boost's start-map walk recurses inside the group.
    sites: StartMapSites,
    /// Deepest group nesting the translation of the group's contents could have
    /// (see [`Translator::last_nest`]).
    nest: usize,
}

impl Frame {
    fn new(kind: GroupKind, saved_flags: Flags, start: usize) -> Self {
        Self {
            kind,
            saved_flags,
            start,
            group: 0,
            width: Some(0),
            alternative_width: None,
            min_length: MinLength::default(),
            alternative_min_length: None,
            has_states: false,
            alternation: false,
            items: 0,
            last_item_is_lookaround: false,
            consumes: false,
            asserts: false,
            has_capture: false,
            atoms: 0,
            alternative_atoms: 0,
            holds_leading_repeat: false,
            leading_before: false,
            state_icase: [false; 2],
            next_alternation_icase: saved_flags.icase,
            sites: StartMapSites::default(),
            nest: 0,
        }
    }

    fn finish_alternative(&mut self) {
        self.alternative_atoms = self.alternative_atoms.saturating_add(self.atoms);
        self.atoms = 0;
        self.alternative_width = Some(match self.alternative_width {
            None => self.width,
            Some(previous) if previous == self.width => previous,
            Some(_) => None,
        });
        self.width = Some(0);
        self.alternative_min_length = Some(
            self.alternative_min_length
                .map_or(self.min_length, |previous| {
                    previous.shorter(self.min_length)
                }),
        );
        self.min_length = MinLength::default();
    }

    fn body(&self) -> Body {
        if self.alternation {
            Body::Other
        } else if self.items == 0 {
            Body::Empty
        } else if self.items == 1 && self.last_item_is_lookaround {
            Body::LookAround
        } else {
            Body::Other
        }
    }
}

struct Translation {
    pattern: String,
    mark_count: usize,
    names: Vec<(String, usize)>,
    /// Automaton size of the whole expression (see [`Translator::last_atoms`]).
    atoms: usize,
}

struct Translator<'p> {
    source: &'p str,
    bytes: &'p [u8],
    position: usize,
    out: String,
    flags: Flags,
    frames: Vec<Frame>,
    mark_count: usize,
    names: Vec<(String, usize)>,
    max_backreference: usize,
    /// Shortest match of every capturing group by number, recorded when the
    /// group closes (index 0 is unused).
    group_min_lengths: Vec<MinLength>,
    last: Last,
    /// Shortest match of `last`, removed again when a quantifier rescales it.
    last_min_length: MinLength,
    /// Where the output of `last` starts, when it is a byte atom or a
    /// backreference.
    last_start: usize,
    /// Emit the counted spelling (see [`EngineTexts`]).
    counted: bool,
    /// Wrap every repeat that is unbounded or optional (`*`, `+`, `?`, `{n,}`,
    /// greedy or lazy) as `(?:X|NO_BYTE)`, an alternation with a branch that
    /// matches no haystack character, so that the engine's optimizer does not
    /// rewrite it (see [`optimizer_rewrites`]).
    shield_repeats: bool,
    /// Automaton size of `last`, removed again when a quantifier rescales it.
    ///
    /// The atoms of a translation are its one-byte atoms, assertions and
    /// backreferences, the guard branch of every guarded alternation and the extra
    /// branch of every shielded repeat, with every repeat written out: a repeat
    /// multiplies the atoms of what it repeats by its maximum, or by its minimum
    /// (at least 1) when unbounded, as the engine's Thompson construction copies
    /// them. The automaton the engine builds for an expression, and the work its
    /// slowest search mode does per haystack byte, grow with this count.
    last_atoms: usize,
    /// Every state Boost has emitted so far is one its `probe_leading_repeat`
    /// steps over: group starts and ends, flag groups that keep case
    /// sensitivity, `^ $ \b \B \< \> \A \z`, and whole lookarounds. A `|`
    /// (which Boost inserts at the start of its group), a flag group that
    /// changes case sensitivity (a `toggle_case` state), a quantifier (a repeat
    /// state in front of what it repeats) or anything that consumes ends it.
    leading_open: bool,
    /// `leading_open` as it was just before the last one-byte atom.
    last_leading: bool,
    /// Where the leading lazy repeat starts, while it is still leading.
    leading_repeat: Option<usize>,
    /// A flag group switches case sensitivity somewhere: Boost has a
    /// `toggle_case` state.
    case_switch: bool,
    /// Where the first `\<` is.
    word_start: Option<usize>,
    /// Where the first `\<` or `\>` is.
    word_boundary: Option<usize>,
    /// Where the first quantifier on a group holding a repeat or an alternation
    /// is, a group whose start map Boost's walk loops back to.
    looped_group: Option<usize>,
    /// Group nesting the translation of `last` could have, removed again when a
    /// quantifier wraps it: an upper bound over both spellings that counts a
    /// one-byte atom or backreference as one level (`(?i:\x61)`, `(?s:.)`), an
    /// assertion as two (the line anchors), a group as one level above its
    /// contents (two for a positive lookahead, which may end in `(?=)`), and a
    /// quantifier as one level above a repeated lookaround or a zero repeat of a
    /// capturing group (which the translation rewrites), and otherwise as one
    /// level above what it repeats for a counted repeat (whose iteration counter
    /// is two levels deep) and one more for a shielded repeat, whether or not the
    /// spelling applies them.
    last_nest: usize,
    /// The expression holds `\z`, `` \' `` or `$` under `(?-m)`, Boost's
    /// `buffer_end`: a start-map walk that ends there adds no byte, so the map of
    /// an alternation all of whose ways end there stays unbuilt and is walked
    /// again each time another walk reaches it.
    buffer_end: bool,
}

impl<'p> Translator<'p> {
    fn new(source: &'p str, options: RegexOptions, counted: bool) -> Self {
        let flags = Flags {
            icase: options.icase,
            dotall: !options.no_mod_s,
            multiline: true,
        };
        Self {
            source,
            bytes: source.as_bytes(),
            position: 0,
            out: String::with_capacity(source.len().saturating_mul(4).saturating_add(16)),
            flags,
            frames: vec![Frame::new(GroupKind::Root, flags, 0)],
            mark_count: 0,
            names: Vec::new(),
            max_backreference: 0,
            group_min_lengths: vec![MinLength::default()],
            last: Last::Nothing,
            last_min_length: MinLength::default(),
            last_start: 0,
            counted,
            shield_repeats: false,
            last_atoms: 0,
            leading_open: true,
            last_leading: false,
            leading_repeat: None,
            case_switch: false,
            word_start: None,
            word_boundary: None,
            looped_group: None,
            last_nest: 0,
            buffer_end: false,
        }
    }

    fn run(mut self) -> Result<Translation> {
        while let Some(&byte) = self.bytes.get(self.position) {
            match byte {
                b'\\' => self.escape()?,
                b'[' => self.class()?,
                b'(' => self.open_group()?,
                b')' => self.close_group()?,
                b'|' => self.alternation(),
                b'*' | b'+' | b'?' => self.simple_quantifier(byte)?,
                b'{' => self.brace()?,
                b'^' => {
                    self.position += 1;
                    let anchor = if self.flags.multiline {
                        START_OF_LINE
                    } else {
                        r"\A"
                    };
                    self.assertion(anchor);
                }
                b'$' => {
                    self.position += 1;
                    let anchor = if self.flags.multiline {
                        // Boost's `end_line` recurses in the start-map walk.
                        self.add_start_map_anchor();
                        END_OF_LINE
                    } else {
                        self.buffer_end = true;
                        r"\z"
                    };
                    self.assertion(anchor);
                }
                b'.' => {
                    self.position += 1;
                    let dot = if self.flags.dotall {
                        "(?s:.)"
                    } else {
                        DOT_NO_SEPARATOR
                    };
                    self.last_start = self.out.len();
                    self.out.push_str(dot);
                    self.byte_atom();
                }
                0x80..=0xFF => return Err(self.unsupported(self.position, "a non-ASCII byte")),
                _ => {
                    self.position += 1;
                    self.literal(byte);
                }
            }
        }
        if self.frames.len() > 1 {
            return Err(self.syntax(self.bytes.len(), "missing ')'"));
        }
        if self.max_backreference > self.mark_count {
            return Err(self.syntax(
                self.bytes.len(),
                "a backreference names a group the expression does not have",
            ));
        }
        // Boost marks a leading repeat only in an expression without
        // backreferences (`basic_regex_creator::probe_leading_repeat`).
        if let (Some(start), 0) = (self.leading_repeat, self.max_backreference) {
            return Err(self.unsupported(start, LEADING_LAZY_REPEAT));
        }
        if let (Some(at), true) = (self.word_start, self.case_switch) {
            return Err(self.unsupported(at, WORD_START_CASE_SWITCH));
        }
        if let (Some(at), Some(_)) = (self.word_boundary, self.looped_group) {
            return Err(self.unsupported(at, WORD_BOUNDARY_AFTER_LOOP));
        }
        let mut root = self
            .frames
            .pop()
            .unwrap_or_else(|| Frame::new(GroupKind::Root, self.flags, 0));
        if self.start_map_depth(root.sites) > BOOST_START_MAP_DEPTH {
            return Err(self.unsupported(self.bytes.len(), START_MAP_RECURSION));
        }
        // The programs wrap the spelling in two groups: `\A(?:(?:X)(?=)(?=))\z`.
        if root.nest.saturating_add(2) >= ENGINE_MAX_NESTING {
            return Err(self.unsupported(self.bytes.len(), TRANSLATION_TOO_DEEP));
        }
        root.finish_alternative();
        self.guard_alternation(&mut root);
        Ok(Translation {
            pattern: self.out,
            mark_count: self.mark_count,
            names: self.names,
            atoms: root.alternative_atoms,
        })
    }

    fn syntax(&self, offset: usize, reason: &str) -> Error {
        Error::InvalidValue(format!(
            "invalid regular expression {} at byte {offset}: {reason}",
            quote(self.source)
        ))
    }

    fn unsupported(&self, offset: usize, what: &str) -> Error {
        Error::Unsupported(format!(
            "regular expression {} uses {what} at byte {offset}, which the Boost.Regex facade does not translate",
            quote(self.source)
        ))
    }

    fn frame(&mut self) -> &mut Frame {
        // The root frame is never popped, so the stack is never empty.
        let last = self.frames.len() - 1;
        &mut self.frames[last]
    }

    fn add_item(&mut self, lookaround: bool) {
        let frame = self.frame();
        frame.items += 1;
        frame.last_item_is_lookaround = lookaround;
        frame.has_states = true;
    }

    fn add_width(&mut self, width: usize) {
        let frame = self.frame();
        frame.width = frame
            .width
            .and_then(|current| current.checked_add(width))
            .filter(|&total| total <= i32::MAX as usize);
    }

    fn clear_width(&mut self) {
        self.frame().width = None;
    }

    fn add_min_length(&mut self, length: MinLength) {
        let frame = self.frame();
        frame.min_length = frame.min_length.followed_by(length);
        self.last_min_length = length;
    }

    fn add_atoms(&mut self, atoms: usize) {
        let frame = self.frame();
        frame.atoms = frame.atoms.saturating_add(atoms);
        self.last_atoms = atoms;
    }

    fn add_nest(&mut self, nest: usize) {
        let frame = self.frame();
        frame.nest = frame.nest.max(nest);
        self.last_nest = nest;
    }

    /// Count a `$` line anchor, `\<` or `\>`, at which Boost's start-map walk
    /// recurses.
    fn add_start_map_anchor(&mut self) {
        let sites = &mut self.frame().sites;
        sites.free_anchors = sites.free_anchors.saturating_add(1);
    }

    /// An upper bound on the recursion depth of Boost's `create_startmap` over
    /// the expression whose sites are `sites`; Boost throws `error_complexity`
    /// at construction when the depth exceeds [`BOOST_START_MAP_DEPTH`].
    ///
    /// Boost builds the maps of its repeat and alternation states from the last to
    /// the first, then the expression's own, each by a walk over the states that
    /// can follow. The walk recurses, one level deeper than the call it is in, at
    /// `$` (line anchor), `\<` and `\>`; and two levels deeper (the second of two
    /// calls) at a repeat or alternation whose map is not built yet. Such a state
    /// is either earlier in the expression, inside a repeated group the walk loops
    /// back into, or an alternation all of whose ways end at a buffer end, whose
    /// map stays empty. A repeat is recursed into at most once per map
    /// (`set_bad_repeat`). A walk moves forward except where it loops back, and
    /// after it loops back to a repeat it stays inside that repeat's group, whose
    /// end leads to the marked repeat and stops. So a chain of recursions meets
    /// each site at most once, except a site it passed on its way out of the
    /// repeated groups around the state whose map is built, which it can meet once
    /// more after looping back; on that way lie only anchors and alternations whose
    /// maps stay empty, since every other state there comes later and has its map.
    /// With `A` the anchors and `L` the alternations, each free (outside every
    /// repeated group) or looped (inside one), and `R` the repeated groups, the
    /// depth is at most `A_free + 2 A_looped + 2 R + 2 L_looped`, and with a buffer
    /// end `2 (L_free + L_looped)` more, plus 2 as a margin.
    fn start_map_depth(&self, sites: StartMapSites) -> usize {
        let mut depth = sites
            .free_anchors
            .saturating_add(sites.looped_anchors.saturating_mul(2))
            .saturating_add(sites.repeated_groups.saturating_mul(2))
            .saturating_add(sites.looped_alternations.saturating_mul(2))
            .saturating_add(2);
        if self.buffer_end {
            depth = depth
                .saturating_add(sites.free_alternations.saturating_mul(2))
                .saturating_add(sites.looped_alternations.saturating_mul(2));
        }
        depth
    }

    fn literal(&mut self, byte: u8) {
        self.last_start = self.out.len();
        if self.flags.icase && byte.is_ascii_alphabetic() {
            let _ = write!(self.out, r"(?i:\x{byte:02X})");
        } else {
            let _ = write!(self.out, r"\x{byte:02X}");
        }
        self.byte_atom();
    }

    fn byte_atom(&mut self) {
        self.add_item(false);
        self.add_width(1);
        self.add_min_length(MinLength::BYTE);
        self.add_atoms(1);
        self.add_nest(1);
        self.frame().consumes = true;
        self.last = Last::Byte;
        self.last_leading = self.leading_open;
        self.leading_open = false;
    }

    fn assertion(&mut self, text: &str) {
        self.out.push_str(text);
        self.add_item(false);
        self.add_min_length(MinLength::default());
        self.add_atoms(1);
        self.add_nest(2);
        self.frame().asserts = true;
        self.last = Last::Fixed;
    }

    fn emit_class(&mut self, negated: bool, items: &[ClassItem]) {
        self.last_start = self.out.len();
        let (accepted, high) = class_set(self.flags.icase, negated, items);
        self.out.push_str(&class_text(&accepted, high));
        self.byte_atom();
    }

    fn escape(&mut self) -> Result<()> {
        let start = self.position;
        let Some(&letter) = self.bytes.get(start + 1) else {
            return Err(self.syntax(start, "escape sequence terminated prematurely"));
        };
        self.position += 2;
        if let Some((kind, negated)) = ClassKind::escape(letter) {
            self.emit_class(negated, &[ClassItem::Class(kind, false)]);
            return Ok(());
        }
        match letter {
            b'b' => self.assertion(r"\b"),
            b'B' => self.assertion(r"\B"),
            b'<' | b'>' => {
                if letter == b'<' {
                    self.word_start.get_or_insert(start);
                }
                self.word_boundary.get_or_insert(start);
                self.add_start_map_anchor();
                self.assertion(if letter == b'<' { r"\<" } else { r"\>" });
            }
            b'A' | b'`' => self.assertion(r"\A"),
            b'z' | b'\'' => {
                self.buffer_end = true;
                self.assertion(r"\z");
            }
            b'a' => self.literal(0x07),
            b'e' => self.literal(0x1B),
            b'f' => self.literal(0x0C),
            b'n' => self.literal(b'\n'),
            b'r' => self.literal(b'\r'),
            b't' => self.literal(b'\t'),
            b'x' => {
                let byte = self.hex_escape(start)?;
                self.literal(byte);
            }
            b'1'..=b'9' => {
                if self
                    .bytes
                    .get(self.position)
                    .is_some_and(u8::is_ascii_digit)
                {
                    return Err(self.unsupported(start, "a backreference with more than one digit"));
                }
                let group = usize::from(letter - b'0');
                if self
                    .frames
                    .iter()
                    .any(|frame| frame.kind == GroupKind::Capture && frame.group == group)
                {
                    return Err(self.unsupported(start, OPEN_GROUP_BACKREFERENCE));
                }
                self.max_backreference = self.max_backreference.max(group);
                self.last_start = self.out.len();
                if self.flags.icase {
                    let _ = write!(self.out, r"(?i:\{group})");
                } else {
                    let _ = write!(self.out, r"\{group}");
                }
                // A group that comes later has no recorded length: zero, as in
                // the engine's analyzer. (One that is still open was refused
                // above.)
                let target = self
                    .group_min_lengths
                    .get(group)
                    .copied()
                    .unwrap_or_default();
                self.add_item(false);
                self.clear_width();
                self.add_min_length(MinLength {
                    boost: 0,
                    engine: target.engine,
                });
                self.add_atoms(1);
                self.add_nest(1);
                self.leading_open = false;
                self.frame().consumes = true;
                self.last = Last::Backref {
                    nullable: target.boost == 0,
                };
            }
            b'Z' => {
                return Err(self.unsupported(
                    start,
                    r"\Z (Boost never tries it at a form feed when it starts the expression)",
                ));
            }
            b'Q' | b'E' => return Err(self.unsupported(start, r"\Q...\E quoting")),
            b'K' => return Err(self.unsupported(start, r"\K")),
            b'G' => return Err(self.unsupported(start, r"\G")),
            b'0' => return Err(self.unsupported(start, "an octal escape")),
            b'c' => return Err(self.unsupported(start, "a control-character escape")),
            b'g' | b'k' => return Err(self.unsupported(start, "a named or relative backreference")),
            b'p' | b'P' => return Err(self.unsupported(start, "a character property escape")),
            b'N' | b'R' | b'X' | b'C' => {
                return Err(
                    self.unsupported(start, "a named-character, line-ending or grapheme escape")
                );
            }
            byte if byte.is_ascii_alphanumeric() => {
                return Err(
                    self.unsupported(start, "an escape letter without a translated meaning")
                );
            }
            0x80..=0xFF => return Err(self.unsupported(start, "a non-ASCII byte")),
            byte => self.literal(byte),
        }
        Ok(())
    }

    /// `\xHH` (one or two hex digits) or `\x{H...}`, after the `x`.
    ///
    /// Digits that start with white space, a sign or `0x` are refused
    /// ([`HEX_ESCAPE_STREAM_SYNTAX`]); every other `\x` escape reads as Boost's
    /// `std::istream` reads it.
    fn hex_escape(&mut self, start: usize) -> Result<u8> {
        let Some(&next) = self.bytes.get(self.position) else {
            return Err(self.syntax(start, "hexadecimal escape sequence terminated prematurely"));
        };
        let digits_start = self.position + usize::from(next == b'{');
        let stream_syntax = match self.bytes.get(digits_start) {
            Some(b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r' | b'+' | b'-') => true,
            Some(b'0') => matches!(self.bytes.get(digits_start + 1), Some(b'x' | b'X')),
            _ => false,
        };
        if stream_syntax {
            return Err(self.unsupported(start, HEX_ESCAPE_STREAM_SYNTAX));
        }
        let value = if next == b'{' {
            self.position += 1;
            let digits = self.bytes[self.position..]
                .iter()
                .take_while(|byte| byte.is_ascii_hexdigit())
                .count();
            if digits == 0 || self.bytes.get(self.position + digits) != Some(&b'}') {
                return Err(self.syntax(start, "invalid hexadecimal escape sequence"));
            }
            let text = &self.source[self.position..self.position + digits];
            self.position += digits + 1;
            u32::from_str_radix(text, 16).unwrap_or(u32::MAX)
        } else {
            let digits = self.bytes[self.position..]
                .iter()
                .take(2)
                .take_while(|byte| byte.is_ascii_hexdigit())
                .count();
            if digits == 0 {
                return Err(self.syntax(start, "escape sequence did not encode a valid character"));
            }
            let text = &self.source[self.position..self.position + digits];
            self.position += digits;
            u32::from_str_radix(text, 16).unwrap_or(u32::MAX)
        };
        match u8::try_from(value) {
            Ok(byte) if byte.is_ascii() => Ok(byte),
            _ => Err(self.unsupported(start, "a hexadecimal escape above 0x7F")),
        }
    }

    fn class(&mut self) -> Result<()> {
        let start = self.position;
        self.position += 1;
        let mut negated = false;
        if self.bytes.get(self.position) == Some(&b'^') {
            negated = true;
            self.position += 1;
        }
        let mut items: Vec<ClassItem> = Vec::new();
        loop {
            let Some(&byte) = self.bytes.get(self.position) else {
                return Err(self.syntax(start, "unterminated character class"));
            };
            if byte == b']' && !items.is_empty() {
                self.position += 1;
                break;
            }
            if byte == b'[' {
                match self.bytes.get(self.position + 1) {
                    Some(b':') => {
                        items.push(self.posix_class(start)?);
                        continue;
                    }
                    Some(b'.') => {
                        return Err(self.unsupported(self.position, "a collating element"));
                    }
                    Some(b'=') => {
                        return Err(self.unsupported(self.position, "an equivalence class"));
                    }
                    _ => {}
                }
            }
            if byte == b'\\' {
                let letter = self.bytes.get(self.position + 1).copied();
                if let Some((kind, negated_class)) = letter.and_then(ClassKind::escape) {
                    if kind == ClassKind::Vertical {
                        // `[\v]` is a vertical tab in Boost; `[\V]` is not translated.
                        if negated_class {
                            return Err(
                                self.unsupported(self.position, r"\V inside a character class")
                            );
                        }
                    } else {
                        self.position += 2;
                        items.push(ClassItem::Class(kind, negated_class));
                        continue;
                    }
                }
            }
            self.set_literal(start, &mut items)?;
        }
        self.emit_class(negated, &items);
        Ok(())
    }

    /// `[:name:]` or `[:^name:]` inside a bracket expression, at its `[`.
    fn posix_class(&mut self, class_start: usize) -> Result<ClassItem> {
        let start = self.position;
        let body = start + 2;
        let Some(length) = self.bytes[body..].windows(2).position(|pair| pair == b":]") else {
            return Err(self.unsupported(start, "a '[:' that does not close with ':]'"));
        };
        let mut name = &self.bytes[body..body + length];
        let negated = name.first() == Some(&b'^');
        if negated {
            name = &name[1..];
        }
        let Some(kind) = ClassKind::posix(name) else {
            return if !name.is_empty() && name.iter().all(u8::is_ascii_alphabetic) {
                Err(self.syntax(class_start, "unknown character class name"))
            } else {
                Err(self.unsupported(start, "a character class name that is not a POSIX class"))
            };
        };
        self.position = body + length + 2;
        Ok(ClassItem::Class(kind, negated))
    }

    /// Boost `parse_set_literal`: one literal or a range, at the current byte.
    fn set_literal(&mut self, class_start: usize, items: &mut Vec<ClassItem>) -> Result<()> {
        let first = self.next_set_literal(class_start, items.is_empty())?;
        if self.position >= self.bytes.len() {
            return Err(self.syntax(class_start, "unterminated character class"));
        }
        if self.bytes[self.position] == b'-' {
            self.position += 1;
            let Some(&after) = self.bytes.get(self.position) else {
                return Err(self.syntax(class_start, "unterminated character class"));
            };
            if after != b']' {
                let last = self.next_set_literal(class_start, items.is_empty())?;
                // Boost orders the endpoints after translating them, which under
                // icase lower-cases both: `[a-Z]` is valid there, `[Z-a]` is not.
                let icase = self.flags.icase;
                let translate = |byte: u8| {
                    if icase {
                        byte.to_ascii_lowercase()
                    } else {
                        byte
                    }
                };
                if translate(first) > translate(last) {
                    return Err(self.syntax(class_start, "invalid range in a character class"));
                }
                items.push(ClassItem::Range(first, last));
                if self.bytes.get(self.position) == Some(&b'-') {
                    match self.bytes.get(self.position + 1) {
                        None => {
                            return Err(self.syntax(class_start, "unterminated character class"));
                        }
                        Some(b']') => {}
                        Some(_) => {
                            return Err(
                                self.syntax(class_start, "invalid range in a character class")
                            );
                        }
                    }
                }
                return Ok(());
            }
            self.position -= 1;
        }
        items.push(ClassItem::Byte(first));
        Ok(())
    }

    /// Boost `get_next_set_literal`.
    fn next_set_literal(&mut self, class_start: usize, set_empty: bool) -> Result<u8> {
        let start = self.position;
        let Some(&byte) = self.bytes.get(start) else {
            return Err(self.syntax(class_start, "unterminated character class"));
        };
        match byte {
            b'-' => {
                if !set_empty && self.bytes.get(start + 1) != Some(&b']') {
                    return Err(self.syntax(class_start, "invalid range in a character class"));
                }
                self.position += 1;
                Ok(b'-')
            }
            b'\\' => {
                let Some(&letter) = self.bytes.get(start + 1) else {
                    return Err(self.syntax(class_start, "unterminated character class"));
                };
                self.position += 2;
                match letter {
                    b'a' => Ok(0x07),
                    b'b' => Ok(0x08),
                    b'e' => Ok(0x1B),
                    b'f' => Ok(0x0C),
                    b'n' => Ok(b'\n'),
                    b'r' => Ok(b'\r'),
                    b't' => Ok(b'\t'),
                    b'v' => Ok(0x0B),
                    b'x' => self.hex_escape(start),
                    0x80..=0xFF => Err(self.unsupported(start, "a non-ASCII byte")),
                    letter if letter.is_ascii_alphanumeric() => Err(self.unsupported(
                        start,
                        "an escape inside a character class without a translated meaning",
                    )),
                    other => Ok(other),
                }
            }
            0x80..=0xFF => Err(self.unsupported(start, "a non-ASCII byte")),
            other => {
                self.position += 1;
                Ok(other)
            }
        }
    }

    fn open_group(&mut self) -> Result<()> {
        let start = self.position;
        if self.frames.len() > MAX_GROUP_DEPTH {
            return Err(self.unsupported(start, "group nesting deeper than MAX_GROUP_DEPTH"));
        }
        self.position += 1;
        match self.bytes.get(self.position) {
            None => Err(self.syntax(start, "unmatched '('")),
            Some(b'?') => self.extension(start),
            Some(b'*') => Err(self.unsupported(start, "a backtracking control verb")),
            Some(_) => self.open_capture(start),
        }
    }

    fn open_capture(&mut self, start: usize) -> Result<()> {
        if self.frames.iter().any(|frame| frame.kind.is_independent()) {
            return Err(self.unsupported(
                start,
                "a capturing group inside a lookaround or atomic group (Boost keeps its \
                 capture when the surrounding match backtracks)",
            ));
        }
        self.mark_count += 1;
        self.group_min_lengths.push(MinLength::default());
        self.push_group(GroupKind::Capture, "(", self.flags);
        self.frame().group = self.mark_count;
        Ok(())
    }

    fn push_group(&mut self, kind: GroupKind, text: &str, saved_flags: Flags) {
        let start = self.out.len();
        self.out.push_str(text);
        self.frame().has_states = true;
        let mut frame = Frame::new(kind, saved_flags, start);
        // Boost's `probe_leading_repeat` enters groups but skips a lookaround
        // whole, so nothing inside one is leading.
        if kind.is_lookaround() {
            frame.leading_before = self.leading_open;
            self.leading_open = false;
        }
        self.frames.push(frame);
        self.last = Last::Nothing;
    }

    fn extension(&mut self, start: usize) -> Result<()> {
        self.position += 1;
        let Some(&kind) = self.bytes.get(self.position) else {
            return Err(self.syntax(start, "incomplete '(?' extension"));
        };
        match kind {
            b'#' => {
                // A comment appends no state and ends at the first ')' (or the pattern end).
                while let Some(&byte) = self.bytes.get(self.position) {
                    self.position += 1;
                    if byte == b')' {
                        break;
                    }
                }
                Ok(())
            }
            b':' => {
                self.position += 1;
                self.push_group(GroupKind::NonCapture, "(?:", self.flags);
                Ok(())
            }
            b'=' => {
                self.position += 1;
                self.push_group(GroupKind::LookAhead, "(?=", self.flags);
                Ok(())
            }
            b'!' => {
                self.position += 1;
                self.push_group(GroupKind::NegativeLookAhead, "(?!", self.flags);
                Ok(())
            }
            b'>' => {
                self.position += 1;
                self.push_group(GroupKind::Atomic, "(?>", self.flags);
                Ok(())
            }
            b'<' => match self.bytes.get(self.position + 1) {
                Some(b'=') => {
                    self.position += 2;
                    self.push_group(GroupKind::LookBehind, "(?<=", self.flags);
                    Ok(())
                }
                Some(b'!') => {
                    self.position += 2;
                    self.push_group(GroupKind::NegativeLookBehind, "(?<!", self.flags);
                    Ok(())
                }
                _ => self.named_group(start, b'>'),
            },
            b'\'' => self.named_group(start, b'\''),
            b'P' => Err(self.syntax(start, "Python-style '(?P' groups are not Boost syntax")),
            b'|' => Err(self.unsupported(start, "a branch-reset group")),
            b'(' => Err(self.unsupported(start, "a conditional expression")),
            b'R' | b'&' | b'+' | b'0'..=b'9' => {
                Err(self.unsupported(start, "a recursive sub-expression"))
            }
            b'-' if self
                .bytes
                .get(self.position + 1)
                .is_some_and(u8::is_ascii_digit) =>
            {
                Err(self.unsupported(start, "a recursive sub-expression"))
            }
            b')' => Err(self.syntax(start, "empty '(?)' extension")),
            _ => self.flag_group(start),
        }
    }

    fn named_group(&mut self, start: usize, delimiter: u8) -> Result<()> {
        self.position += 1;
        let name_start = self.position;
        loop {
            match self.bytes.get(self.position) {
                None => return Err(self.syntax(start, "unterminated named capture")),
                Some(&byte) if byte == delimiter => break,
                Some(&byte) if byte.is_ascii_alphanumeric() || byte == b'_' => self.position += 1,
                Some(_) => {
                    return Err(self.unsupported(start, "a group name outside [A-Za-z0-9_]"));
                }
            }
        }
        let name = self.source[name_start..self.position].to_string();
        self.position += 1;
        self.open_capture(start)?;
        self.names.push((name, self.mark_count));
        Ok(())
    }

    fn flag_group(&mut self, start: usize) -> Result<()> {
        let mut flags = self.flags;
        let mut enable = true;
        loop {
            match self.bytes.get(self.position) {
                None => return Err(self.syntax(start, "unterminated '(?' modifier group")),
                Some(b'i') => flags.icase = enable,
                Some(b'm') => flags.multiline = enable,
                Some(b's') => flags.dotall = enable,
                Some(b'x') => return Err(self.unsupported(start, "the x (extended) modifier")),
                Some(b'-') if enable => enable = false,
                Some(_) => break,
            }
            self.position += 1;
        }
        match self.bytes.get(self.position) {
            Some(b')') => {
                self.position += 1;
                let case_change = flags.icase != self.flags.icase;
                self.case_switch |= case_change;
                self.flags = flags;
                self.frame().has_states = true;
                // Boost appends a `toggle_case` state for a change of case
                // sensitivity, which its `probe_leading_repeat` does not step
                // over.
                if case_change {
                    self.leading_open = false;
                }
                self.last = Last::EmptyGroup { case_change };
                self.last_min_length = MinLength::default();
                self.last_atoms = 0;
                self.last_nest = 0;
                Ok(())
            }
            Some(b':') => {
                self.position += 1;
                let saved = self.flags;
                self.flags = flags;
                self.push_group(GroupKind::NonCapture, "(?:", saved);
                if flags.icase != saved.icase {
                    self.case_switch = true;
                    self.leading_open = false;
                }
                Ok(())
            }
            _ => Err(self.syntax(start, "unknown '(?' extension or modifier")),
        }
    }

    fn close_group(&mut self) -> Result<()> {
        let start = self.position;
        if self.frames.len() == 1 {
            return Err(self.syntax(start, "unmatched ')'"));
        }
        self.position += 1;
        let Some(mut frame) = self.frames.pop() else {
            return Err(self.syntax(start, "unmatched ')'"));
        };
        frame.finish_alternative();
        let width = frame.alternative_width.flatten();
        let min_length = frame.alternative_min_length.unwrap_or_default();
        if let (GroupKind::Capture, Some(recorded)) =
            (frame.kind, self.group_min_lengths.get_mut(frame.group))
        {
            *recorded = min_length;
        }
        if matches!(frame.kind, GroupKind::LookAhead | GroupKind::Atomic) && !frame.has_states {
            return Err(self.syntax(start, "invalid or empty zero-width assertion"));
        }
        if matches!(
            frame.kind,
            GroupKind::LookBehind | GroupKind::NegativeLookBehind
        ) && width.is_none()
        {
            return Err(self.syntax(
                start,
                "invalid lookbehind assertion (no single fixed width)",
            ));
        }
        if matches!(
            frame.kind,
            GroupKind::LookBehind | GroupKind::NegativeLookBehind
        ) && width.is_some_and(|width| width > MAX_LOOKBEHIND_WIDTH)
        {
            return Err(
                self.unsupported(start, "a lookbehind wider than MAX_LOOKBEHIND_WIDTH bytes")
            );
        }
        if matches!(
            frame.kind,
            GroupKind::LookBehind | GroupKind::NegativeLookBehind
        ) && frame
            .sites
            .free_alternations
            .saturating_add(frame.sites.looped_alternations)
            > BOOST_BACKSTEP_BLOCKS
        {
            return Err(self.unsupported(start, LOOKBEHIND_ALTERNATIONS));
        }
        if self.counted && frame.kind == GroupKind::LookAhead {
            self.out.push_str(FORCE_BACKTRACKING);
        }
        self.guard_alternation(&mut frame);
        self.out.push(')');
        self.flags = frame.saved_flags;
        if frame.kind.is_lookaround() {
            self.leading_open = frame.leading_before;
        }
        let body = frame.body();
        let captures = frame.has_capture || frame.kind == GroupKind::Capture;
        if frame.kind.is_lookaround() {
            self.add_item(true);
            self.add_min_length(MinLength::default());
            self.add_atoms(frame.alternative_atoms);
        } else {
            match (frame.kind, body) {
                (GroupKind::NonCapture, Body::Empty) => self.frame().has_states = true,
                (GroupKind::NonCapture, Body::LookAround) => self.add_item(true),
                _ => self.add_item(false),
            }
            match width {
                Some(width) => self.add_width(width),
                None => self.clear_width(),
            }
            self.add_min_length(min_length);
            self.add_atoms(frame.alternative_atoms);
            let parent = self.frame();
            parent.consumes |= frame.consumes;
            parent.asserts |= frame.asserts;
            parent.has_capture |= captures;
        }
        self.add_nest(
            frame
                .nest
                .max(usize::from(frame.kind == GroupKind::LookAhead))
                .saturating_add(1),
        );
        let parent = self.frame();
        parent.state_icase[0] |= frame.state_icase[0];
        parent.state_icase[1] |= frame.state_icase[1];
        parent.sites.add(frame.sites);
        self.last = Last::Group(GroupShape {
            start: frame.start,
            kind: frame.kind,
            body,
            captures,
            nullable: min_length.boost == 0,
            only_asserts: frame.asserts && !frame.consumes,
            holds_leading_repeat: frame.holds_leading_repeat,
            state_icase: frame.state_icase,
            free_anchors: frame.sites.free_anchors,
            free_alternations: frame.sites.free_alternations,
        });
        Ok(())
    }

    /// End the alternation of `frame`, which has just been closed, with a
    /// branch that matches no haystack character ([`NO_BYTE`]) wherever the
    /// engine may hand it to an automaton and its priority matters: in the
    /// plain spelling outside lookarounds, and in both spellings inside an
    /// atomic group, unless the alternation is part of a lookbehind's width.
    ///
    /// `regex-automata` compiles what `fancy-regex` delegates to it from a
    /// `regex-syntax` HIR, and `Hir::alternation` (`lift_common_prefix`,
    /// `regex-syntax` 0.8.11) factors a prefix that all branches share out of
    /// the alternation: `XA|XB` becomes `X(?:A|B)`. That keeps the language but
    /// not leftmost-first priority when `X` can match in more than one way,
    /// because the original tries every way of matching `X` with `A` before it
    /// tries `B`: `[ab]+b|[ab]+c` on `abbc` matches `0..4` after the rewrite and
    /// `0..3` in Boost, and `\w+?a|\w+?()` on `aba` matches `0..1` instead of
    /// `0..3`. The rewrite applies only when every branch is a concatenation,
    /// which the extra branch is not.
    ///
    /// `fancy-regex` delegates three kinds of sub-expression that need no
    /// backtracking (`Compiler::visit`, `compile_concat`): the whole plain
    /// spelling; the body of a lookaround, or its trailing run; and the body of
    /// an atomic group, or a branch or trailing run of it, which it compiles as
    /// if nothing followed. The counted spelling ends in an assertion that needs
    /// backtracking, so outside those bodies it delegates only runs of fixed
    /// size, whose branches all end at the same position. A lookaround only asks
    /// whether its body matches, which the rewrite does not change. An atomic
    /// group keeps the first match of its body and discards the others, so the
    /// rewrite changes which one it keeps: `(?>[ab]?b|[ab]?c)` on `bc` matches
    /// `0..2` rewritten and `0..1` in Boost, also inside a lookahead. An
    /// alternation whose nearest enclosing lookaround is a lookbehind is part of
    /// that lookbehind's width: all its branches have one width, which the extra
    /// branch would break, and every way of matching a body of one width ends at
    /// the same position, so there the branch is not added.
    fn guard_alternation(&mut self, frame: &mut Frame) {
        let enclosing =
            || std::iter::once(frame.kind).chain(self.frames.iter().rev().map(|f| f.kind));
        let nearest_lookaround = enclosing().find(|kind| kind.is_lookaround());
        let in_lookbehind = matches!(
            nearest_lookaround,
            Some(GroupKind::LookBehind | GroupKind::NegativeLookBehind)
        );
        let in_atomic = enclosing().any(|kind| kind == GroupKind::Atomic);
        let on_automaton = !self.counted && nearest_lookaround.is_none();
        if frame.alternation && !in_lookbehind && (in_atomic || on_automaton) {
            self.out.push('|');
            self.out.push_str(NO_BYTE);
            frame.alternative_atoms = frame.alternative_atoms.saturating_add(1);
        }
    }

    fn alternation(&mut self) {
        self.position += 1;
        self.out.push('|');
        // Boost inserts the alternation state at the start of the group (or
        // of the expression), in front of everything the group holds.
        self.leading_open = false;
        let frame = self.frame();
        if frame.holds_leading_repeat {
            self.leading_repeat = None;
        }
        let icase = self.flags.icase;
        let frame = self.frame();
        frame.state_icase[usize::from(frame.next_alternation_icase)] = true;
        frame.next_alternation_icase = icase;
        frame.sites.free_alternations = frame.sites.free_alternations.saturating_add(1);
        frame.finish_alternative();
        frame.alternation = true;
        frame.has_states = true;
        self.last = Last::Nothing;
        self.last_min_length = MinLength::default();
        self.last_atoms = 0;
        self.last_nest = 0;
    }

    fn simple_quantifier(&mut self, symbol: u8) -> Result<()> {
        let start = self.position;
        self.position += 1;
        let (min, max) = match symbol {
            b'*' => (0, None),
            b'+' => (1, None),
            _ => (0, Some(1)),
        };
        self.quantify(start, min, max)
    }

    /// Boost `parse_repeat_range`: a `{` that does not form a valid bound is a
    /// literal.
    fn brace(&mut self) -> Result<()> {
        let start = self.position;
        let skip_space = |bytes: &[u8], mut at: usize| {
            while bytes
                .get(at)
                .is_some_and(|byte| matches!(byte, b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r'))
            {
                at += 1;
            }
            at
        };
        let mut cursor = skip_space(self.bytes, start + 1);
        // A minimum that is missing, negative or outside `intmax_t` leaves the
        // `{` literal.
        let Some((min, after)) = self.bound(cursor).filter(|&(value, _)| value >= 0) else {
            return self.literal_brace();
        };
        cursor = skip_space(self.bytes, after);
        let max = match self.bytes.get(cursor) {
            None => return self.literal_brace(),
            Some(b',') => {
                cursor = skip_space(self.bytes, cursor + 1);
                if cursor >= self.bytes.len() {
                    return self.literal_brace();
                }
                // A missing maximum is unbounded. So is a negative one, which
                // Boost reads and steps over before discarding it; one outside
                // `intmax_t` is not stepped over, so the `}` check below fails.
                match self.bound(cursor) {
                    None => None,
                    Some((value, after)) => {
                        cursor = after;
                        (value >= 0).then_some(value)
                    }
                }
            }
            Some(_) => Some(min),
        };
        cursor = skip_space(self.bytes, cursor);
        if self.bytes.get(cursor) != Some(&b'}') {
            return self.literal_brace();
        }
        if max.is_some_and(|max| min > max) {
            return Err(self.syntax(start, "invalid repeat bounds: minimum above maximum"));
        }
        let bounded = |value: i128| {
            usize::try_from(value)
                .ok()
                .filter(|&value| value <= MAX_REPEAT)
        };
        let (Some(min), Some(max)) = (
            bounded(min),
            max.map_or(Some(None), |max| bounded(max).map(Some)),
        ) else {
            return Err(self.unsupported(start, "a repeat bound above MAX_REPEAT"));
        };
        self.position = cursor + 1;
        self.quantify(start, min, max)
    }

    /// Boost's `cpp_regex_traits::toi`, which reads a bound with
    /// `std::istream >> intmax_t`: an optional sign, then decimal digits, with
    /// leading zeros allowed. `None`, Boost's -1 without consuming anything,
    /// when no digit follows or the value does not fit `intmax_t`; otherwise
    /// the value and the position after its last digit.
    fn bound(&self, at: usize) -> Option<(i128, usize)> {
        let (negative, digits_start) = match self.bytes.get(at) {
            Some(b'-') => (true, at + 1),
            Some(b'+') => (false, at + 1),
            _ => (false, at),
        };
        let digits = self
            .bytes
            .get(digits_start..)?
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count();
        if digits == 0 {
            return None;
        }
        let mut value = 0i128;
        for &digit in self.bytes.get(digits_start..digits_start + digits)? {
            value = value * 10 + i128::from(digit - b'0');
            if value > i128::from(i64::MAX) + 1 {
                return None;
            }
        }
        let value = if negative { -value } else { value };
        (i128::from(i64::MIN)..=i128::from(i64::MAX))
            .contains(&value)
            .then_some((value, digits_start + digits))
    }

    fn literal_brace(&mut self) -> Result<()> {
        self.position += 1;
        self.literal(b'{');
        Ok(())
    }

    fn quantify(&mut self, start: usize, min: usize, max: Option<usize>) -> Result<()> {
        let mut lazy = false;
        if self.bytes.get(self.position) == Some(&b'?') {
            lazy = true;
            self.position += 1;
        }
        if self.bytes.get(self.position) == Some(&b'+') {
            return Err(self.unsupported(start, "a possessive quantifier"));
        }
        if let Some(error) = self.repeat_refusal(start, min, max) {
            return Err(error);
        }
        let mut suffix = match (min, max) {
            (0, None) => "*".to_string(),
            (1, None) => "+".to_string(),
            (0, Some(1)) => "?".to_string(),
            (min, None) => format!("{{{min},}}"),
            (min, Some(max)) if min == max => format!("{{{min}}}"),
            (min, Some(max)) => format!("{{{min},{max}}}"),
        };
        if lazy {
            suffix.push('?');
        }
        let count = self.counted && min >= 2;
        // The shapes `fancy-regex`'s optimizer rewrites are built from repeats
        // that are unbounded or optional (see `optimizer_rewrites`).
        let shield = self.shield_repeats && (max.is_none() || (min, max) == (0, Some(1)));
        let mut shielded = false;
        match self.last {
            // Refused by `repeat_refusal`.
            Last::Nothing | Last::Fixed | Last::EmptyGroup { case_change: true } => {}
            Last::Byte => {
                self.push_repeat(self.last_start, &suffix, count, shield);
                shielded = shield;
                let frame = self.frame();
                frame.width = match max {
                    Some(max) if max == min => frame
                        .width
                        .and_then(|width| width.checked_sub(1))
                        .and_then(|width| width.checked_add(min))
                        .filter(|&total| total <= i32::MAX as usize),
                    _ => None,
                };
            }
            Last::Backref { .. } => {
                self.push_repeat(self.last_start, &suffix, count, shield);
                shielded = shield;
                self.clear_width();
            }
            Last::EmptyGroup { case_change: false } => self.clear_width(),
            Last::Group(shape) => {
                self.clear_width();
                if shape.kind == GroupKind::NonCapture && shape.body == Body::Empty {
                    // The engine has no quantifiable empty group; an empty group
                    // repeated any number of times is still empty.
                } else if shape.is_zero_width() {
                    // The engine refuses to repeat a lookaround. Boost tries a
                    // zero-width body at most once, so the repeat becomes an
                    // optional group (preferring the body when greedy), the body
                    // itself when at least one repetition is required, and a
                    // never-taken alternative when no repetition is allowed.
                    let group = self.out.split_off(shape.start);
                    if max == Some(0) {
                        let _ = write!(self.out, "(?:{group}|(?!)){{0}}");
                        self.frame().last_item_is_lookaround = false;
                    } else if min == 0 {
                        let _ = if lazy {
                            write!(self.out, "(?:|{group})")
                        } else {
                            write!(self.out, "(?:{group}|)")
                        };
                        self.frame().last_item_is_lookaround = false;
                    } else {
                        self.out.push_str(&group);
                    }
                } else if max == Some(0) && shape.captures {
                    // regex-syntax simplifies a zero repeat to nothing and loses
                    // its groups; a never-taken alternative keeps the numbering.
                    let group = self.out.split_off(shape.start);
                    let _ = write!(self.out, "(?:{group}|(?!)){{0}}");
                } else {
                    self.push_repeat(shape.start, &suffix, count, shield);
                    shielded = shield;
                }
            }
        }
        // Boost's `probe_leading_repeat` marks a repeat of a one-byte atom as
        // leading when only transparent states precede it, and a lazy one with
        // room for two more iterations than its minimum then restarts a failed
        // search behind start positions that could still match.
        let leading_lazy = matches!(self.last, Last::Byte)
            && self.last_leading
            && lazy
            && max.is_some_and(|max| max >= min.saturating_add(2));
        if leading_lazy {
            self.leading_repeat = Some(start);
            for frame in &mut self.frames {
                frame.holds_leading_repeat = true;
            }
        } else if let Last::Group(GroupShape {
            holds_leading_repeat: true,
            ..
        }) = self.last
        {
            // Boost inserts the repeat state in front of the group.
            self.leading_repeat = None;
        }
        // Boost puts a repeat state in front of what it repeats. A repeat of a
        // group turns the start-map sites of the group into looped ones.
        let icase = self.flags.icase;
        let group = match self.last {
            Last::Group(shape) => Some(shape),
            _ => None,
        };
        if group.is_some_and(|shape| shape.state_icase.contains(&true)) {
            self.looped_group.get_or_insert(start);
        }
        let repeats_group = matches!(self.last, Last::Group(_) | Last::EmptyGroup { .. });
        let (frame_last, last_nest) = (self.last, self.last_nest);
        let frame = self.frame();
        frame.state_icase[usize::from(icase)] = true;
        if let Some(shape) = group {
            let sites = &mut frame.sites;
            sites.free_anchors = sites.free_anchors.saturating_sub(shape.free_anchors);
            sites.looped_anchors = sites.looped_anchors.saturating_add(shape.free_anchors);
            sites.free_alternations = sites
                .free_alternations
                .saturating_sub(shape.free_alternations);
            sites.looped_alternations = sites
                .looped_alternations
                .saturating_add(shape.free_alternations);
        }
        if repeats_group {
            frame.sites.repeated_groups = frame.sites.repeated_groups.saturating_add(1);
        }
        let wrapped = match frame_last {
            Last::Group(shape) if shape.is_zero_width() || (max == Some(0) && shape.captures) => {
                Some(last_nest.saturating_add(1))
            }
            Last::Group(shape) if !shape.is_repeated_by_engine() => None,
            Last::Byte | Last::Backref { .. } | Last::Group(_) => {
                let mut nest = last_nest;
                if min >= 2 {
                    nest = nest.max(2).saturating_add(1);
                }
                if max.is_none() || (min, max) == (0, Some(1)) {
                    nest = nest.saturating_add(1);
                }
                Some(nest)
            }
            Last::Nothing | Last::Fixed | Last::EmptyGroup { .. } => None,
        };
        if let Some(nest) = wrapped {
            frame.nest = frame.nest.max(nest);
        }
        self.leading_open = false;
        let contribution = self.last_min_length;
        let atoms = self.last_atoms;
        let copies = max.unwrap_or(min.max(1));
        let frame = self.frame();
        frame.min_length = frame
            .min_length
            .without(contribution)
            .followed_by(contribution.repeated(min));
        frame.atoms = frame
            .atoms
            .saturating_sub(atoms)
            .saturating_add(atoms.saturating_mul(copies))
            .saturating_add(usize::from(shielded));
        self.last = Last::Fixed;
        self.last_min_length = MinLength::default();
        self.last_atoms = 0;
        self.last_nest = 0;
        Ok(())
    }

    /// Append `suffix`, a quantifier, to the atom whose output starts at
    /// `start`. With `count`, the atom is first wrapped together with
    /// [`COUNT_ITERATION`]: the engine pushes no backtracking branch for the
    /// iterations below a repeat's minimum, so without it those iterations,
    /// repeated from every start position, would never spend the budget. With
    /// `shield`, the repeat is then wrapped as [`Translator::shield_repeats`]
    /// describes.
    fn push_repeat(&mut self, start: usize, suffix: &str, count: bool, shield: bool) {
        if count {
            let atom = self.out.split_off(start);
            let _ = write!(self.out, "(?:{atom}{COUNT_ITERATION}){suffix}");
        } else {
            self.out.push_str(suffix);
        }
        if shield {
            let repeat = self.out.split_off(start);
            let _ = write!(self.out, "(?:{repeat}|{NO_BYTE})");
        }
    }

    /// Why a repeat of `last` from `min` to `max` is refused, if it is:
    /// Boost's own syntax errors first, then the repeats the facade does not
    /// hand to the engine.
    ///
    /// Boost ends a repeat as soon as one iteration matches the empty string,
    /// even below the minimum: it takes the exit with the iteration accepted
    /// (`repeater_count::check_null_repeat`, `perl_matcher::match_rep`). The
    /// engine behaves differently either way. A bounded repeat iterates up to
    /// its bound, so `(?:\b|a){2}b` matches `ab` there and not in Boost, and
    /// `(?:a{0}\b){999999999}` runs a billion empty iterations. An unbounded
    /// repeat (`RepeatEpsilon`) fails the empty iteration instead and backtracks
    /// into the alternatives of its body, so `(?:b?|a)*` matches all of `ba`
    /// there and only `b` in Boost. Every repeat that allows more than one
    /// iteration of a group that can match empty is therefore refused. A
    /// backreference to a group that is not open where the backreference stands
    /// matches the same text in every iteration of its own repeat, so only its
    /// bounded repeats are refused, for the work; a backreference inside the
    /// group it refers to is refused when it is parsed
    /// ([`OPEN_GROUP_BACKREFERENCE`]). The shortest match of the
    /// repeat, as the engine's analyzer multiplies it out, must also stay within
    /// [`MAX_REPEAT`].
    fn repeat_refusal(&self, start: usize, min: usize, max: Option<usize>) -> Option<Error> {
        let counted = min > 1 || max.is_some_and(|max| max > 1);
        let many = max.is_none_or(|max| max > 1);
        // The engine discards the backtracking branches pushed inside an atomic
        // group or a negative lookaround once its body has matched, without
        // counting them, so a loop there could run to the end of the haystack
        // from every start position without spending the budget. (A loop in a
        // lookahead nested in a negative lookbehind scans forward too.)
        let engine_loop = match self.last {
            Last::Byte | Last::Backref { .. } => true,
            Last::Group(shape) => shape.is_repeated_by_engine(),
            Last::Nothing | Last::Fixed | Last::EmptyGroup { .. } => false,
        };
        let discarded = engine_loop
            && max.is_none_or(|max| max > 1)
            && self.frames.iter().any(|frame| {
                matches!(
                    frame.kind,
                    GroupKind::Atomic
                        | GroupKind::NegativeLookAhead
                        | GroupKind::NegativeLookBehind
                )
            });
        let error = match self.last {
            Last::Nothing => self.syntax(start, "nothing to repeat"),
            Last::Fixed => self.syntax(
                start,
                "a repeat cannot apply to a zero-width assertion or another repeat",
            ),
            Last::EmptyGroup { case_change: true } => self.unsupported(
                start,
                "a repeat of a modifier group that switches case sensitivity (Boost undoes \
                 the switch when the empty repetition is abandoned)",
            ),
            _ if discarded => self.unsupported(
                start,
                "a repeat inside an atomic group or a negative lookaround (the engine discards \
                 its work there without counting it)",
            ),
            Last::Group(shape)
                if shape.captures && shape.nullable && max.is_none_or(|max| max > 1) =>
            {
                self.unsupported(
                    start,
                    "a repeat of a group that can match the empty string and captures (Boost \
                     records a final empty iteration)",
                )
            }
            Last::Group(shape) if shape.only_asserts => {
                self.unsupported(start, "a repeat of a group that only asserts a position")
            }
            Last::Group(shape) if many && shape.nullable && shape.is_repeated_by_engine() => {
                self.unsupported(start, NULLABLE_GROUP_REPEAT)
            }
            Last::Backref { nullable: true } if counted => self.unsupported(
                start,
                "a counted repeat of a backreference that can match the empty string (Boost \
                 ends a repeat after an empty iteration)",
            ),
            Last::Group(shape) if shape.state_icase[usize::from(!self.flags.icase)] => {
                self.unsupported(start, CASE_SWITCHED_REPEAT)
            }
            _ if self
                .last_min_length
                .engine
                .checked_mul(min)
                .is_none_or(|length| length > MAX_REPEAT) =>
            {
                self.unsupported(
                    start,
                    "a repeat whose shortest match is longer than MAX_REPEAT bytes",
                )
            }
            _ => return None,
        };
        Some(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translates_line_anchors_and_dot() {
        let translation = Translator::new("^a.$", RegexOptions::default(), false)
            .run()
            .unwrap();
        assert_eq!(
            translation.pattern,
            format!(r"{START_OF_LINE}\x61(?s:.){END_OF_LINE}")
        );
    }

    #[test]
    fn counted_spelling_counts_iterations_and_forces_backtracking() {
        let translate = |pattern: &str, counted: bool| {
            Translator::new(pattern, RegexOptions::default(), counted)
                .run()
                .unwrap()
                .pattern
        };
        assert_eq!(translate(r"(?=a)b{2}", false), r"(?=\x61)\x62{2}");
        assert_eq!(
            translate(r"(?=a)b{2}", true),
            r"(?=\x61(?=))(?:\x62(?:|(?!))){2}"
        );
        // A minimum of one or less pushes a branch per iteration already.
        assert_eq!(translate(r"(?:ab){0,3}c+", true), r"(?:\x61\x62){0,3}\x63+");
        assert_eq!(translate(r"(a)\1{3,}", true), r"(\x61)(?:\1(?:|(?!))){3,}");
        let (texts, _) = EngineTexts::new(r"(?<=[KR])(?!P)", RegexOptions::default()).unwrap();
        for scale in BACKTRACK_SCALES {
            assert_eq!(texts.search(scale), texts.counted);
        }
        assert!(texts.counted.ends_with(FORCE_BACKTRACKING_ROOT));
        let (texts, _) = EngineTexts::new(r"scan=(\d+)", RegexOptions::default()).unwrap();
        for scale in BACKTRACK_SCALES {
            assert_eq!(texts.search(scale), r"\x73\x63\x61\x6E\x3D([0-9]+)");
        }
    }

    /// The atoms of a translation, and the haystack length from which an
    /// expression that needs no backtracking leaves the automaton.
    #[test]
    fn counts_automaton_atoms() {
        let atoms = |pattern: &str| {
            Translator::new(pattern, RegexOptions::default(), false)
                .run()
                .unwrap()
                .atoms
        };
        assert_eq!(atoms(""), 0);
        assert_eq!(atoms("abc"), 3);
        // An alternation adds the branch that guards it (see `guard_alternation`).
        assert_eq!(atoms("a|bc"), 4);
        assert_eq!(atoms(r"\w{0,9999}b"), 10_000);
        assert_eq!(atoms("(?:ab){3,}c*"), 7);
        assert_eq!(atoms("(?:a{2}|b){5}"), 20);
        assert_eq!(atoms("(?:a(?#c)){4}(?i)"), 4);
        assert_eq!(atoms("(?:(?:a{999}){999}){999}"), 997_002_999);
        let texts = |pattern: &str| {
            EngineTexts::new(pattern, RegexOptions::default())
                .unwrap()
                .0
        };
        let plain =
            |texts: &EngineTexts, scale: usize| texts.plain.as_deref() == Some(texts.search(scale));
        // MAX_AUTOMATON_ATOMS atoms stay on the automaton for the first tier only.
        let widest = texts(&format!(r"\w{{0,{}}}b", MAX_AUTOMATON_ATOMS - 1));
        assert_eq!(widest.atoms, MAX_AUTOMATON_ATOMS);
        assert!(plain(&widest, 1));
        assert!(!plain(&widest, 4));
        let wider = texts(&format!(r"\w{{0,{}}}b", MAX_AUTOMATON_ATOMS));
        assert!(!plain(&wider, 1));
        // The largest expression of the pinned OpenMS sources that needs no
        // backtracking has 46 atoms and stays on the automaton at every length.
        let pinned = texts("TransitionGroupPicker:PeakPickerChromatogram:(.+)");
        assert_eq!(pinned.atoms, 46);
        assert!(plain(&pinned, MAX_BACKTRACK_SCALE));
    }

    /// The counted spelling of an expression that needs no backtracking must
    /// run on the backtracking machine, where every loop iteration spends
    /// budget: a crate update that hands it to an automaton again (see
    /// [`FORCE_BACKTRACKING_ROOT`]) fails here.
    #[test]
    fn bounded_spelling_runs_on_the_backtracking_machine() {
        let (texts, _) = EngineTexts::new("a*b", RegexOptions::default()).unwrap();
        let haystack = "a".repeat(1000);
        let plain = build(texts.plain.as_deref().unwrap(), 1000, false).unwrap();
        assert!(matches!(plain.find(&haystack), Ok(None)));
        let counted = build(&texts.counted, 1000, false).unwrap();
        assert!(matches!(
            counted.find(&haystack),
            Err(fancy_regex::Error::RuntimeError(
                fancy_regex::RuntimeError::BacktrackLimitExceeded
            ))
        ));
        // With one trailing `(?=)` the engine's rewrite hands the whole
        // expression to an automaton, which is why the spelling has two.
        let single = build(
            &format!("(?:{}){FORCE_BACKTRACKING}", texts.plain.unwrap()),
            1000,
            false,
        )
        .unwrap();
        assert!(matches!(single.find(&haystack), Ok(None)));
    }

    /// A backreference inside the group it refers to is refused as soon as it
    /// is read; one to a closed or later group is translated.
    #[test]
    fn refuses_backreferences_inside_their_group() {
        let translate =
            |pattern: &str| Translator::new(pattern, RegexOptions::default(), false).run();
        for pattern in [
            r"(a|\1b)x",
            r"((a)\1)",
            r"(?<n>\1)",
            r"(a(?:b|(?=\1)))",
            r"(a)(b\2)",
        ] {
            match translate(pattern) {
                Err(Error::Unsupported(message)) => {
                    assert!(
                        message.contains(OPEN_GROUP_BACKREFERENCE),
                        "{pattern:?}: {message}"
                    )
                }
                other => panic!(
                    "{pattern:?}: {:?}",
                    other.map(|translation| translation.pattern)
                ),
            }
        }
        for pattern in [r"(a)\1", r"((a)\2)", r"\1(a)", r"(?:(a)|b\1)+", r"(a)(b\1)"] {
            assert!(translate(pattern).is_ok(), "{pattern:?}");
        }
    }

    /// `fancy-regex`'s optimizer rewrites `X+Y?X+`, which matches differently
    /// afterwards; the facade detects it and shields the repeats, and leaves
    /// every expression the optimizer does not touch as it was.
    #[test]
    fn shields_repeats_the_optimizer_would_rewrite() {
        assert!(optimizer_rewrites(r"\x61+\x62?\x61+"));
        assert!(optimizer_rewrites(r"(\x61+)+\1"));
        assert!(!optimizer_rewrites(r"\x61+\x62\x61+"));
        let rewritten = build(r"\x61+\x62?\x61+", 1000, false).unwrap();
        assert!(matches!(rewritten.find("a"), Ok(Some(_))));
        let (texts, _) = EngineTexts::new("a+b?a+", RegexOptions::default()).unwrap();
        let shielded = texts.plain.as_deref().unwrap();
        assert_eq!(
            shielded,
            format!(r"(?:\x61+|{NO_BYTE})(?:\x62?|{NO_BYTE})(?:\x61+|{NO_BYTE})").as_str()
        );
        assert!(!optimizer_rewrites(shielded));
        assert!(!optimizer_rewrites(&texts.counted));
        assert!(matches!(
            build(shielded, 1000, false).unwrap().find("a"),
            Ok(None)
        ));
        let (texts, _) = EngineTexts::new(r"scan=(\d+)", RegexOptions::default()).unwrap();
        assert!(!texts.counted.contains(NO_BYTE));
    }

    /// `regex-syntax` factors a common prefix out of an alternation whose
    /// branches are all concatenations; the plain spelling ends every
    /// alternation outside lookarounds in a branch that matches nothing, and the
    /// counted spelling, which runs on the backtracking machine, does so only
    /// inside atomic groups (see `guards_alternations_inside_atomic_groups`).
    #[test]
    fn guards_alternations_on_the_automaton() {
        let translate = |pattern: &str, counted: bool| {
            Translator::new(pattern, RegexOptions::default(), counted)
                .run()
                .unwrap()
                .pattern
        };
        let plain = translate("[ab]+b|[ab]+c", false);
        assert_eq!(plain, format!(r"[ab]+\x62|[ab]+\x63|{NO_BYTE}"));
        assert_eq!(translate("[ab]+b|[ab]+c", true), r"[ab]+\x62|[ab]+\x63");
        assert_eq!(
            translate("(?<=a|b)(?=c|d)", false),
            r"(?<=\x61|\x62)(?=\x63|\x64)"
        );
        assert_eq!(
            translate("(x|y)z", false),
            format!(r"(\x78|\x79|{NO_BYTE})\x7A")
        );
        let unguarded = build(r"[ab]+\x62|[ab]+\x63", 1000, false).unwrap();
        let guarded = build(&plain, 1000, false).unwrap();
        let range =
            |regex: &fancy_regex::Regex| regex.find("abbc").unwrap().map(|found| found.range());
        assert_eq!(range(&unguarded), Some(0..4));
        assert_eq!(range(&guarded), Some(0..3));
    }

    /// `fancy-regex` hands the body of an atomic group that needs no
    /// backtracking to `regex-automata` whole, where the prefix factoring
    /// changes which match the group keeps: both spellings guard every
    /// alternation inside an atomic group, also inside a lookahead, but not one
    /// that is part of a lookbehind's width.
    #[test]
    fn guards_alternations_inside_atomic_groups() {
        let translate = |pattern: &str, counted: bool| {
            Translator::new(pattern, RegexOptions::default(), counted)
                .run()
                .unwrap()
                .pattern
        };
        for counted in [false, true] {
            assert_eq!(
                translate("(?>[ab]?b|[ab]?c)", counted),
                format!(r"(?>[ab]?\x62|[ab]?\x63|{NO_BYTE})")
            );
            assert_eq!(
                translate("(?>(?:a|ab)c)", counted),
                format!(r"(?>(?:\x61|\x61\x62|{NO_BYTE})\x63)")
            );
        }
        assert!(translate("(?=(?>a|ab)c)", true).contains(NO_BYTE));
        assert!(translate("(?<=(?=(?>a|ab)c).)", true).contains(NO_BYTE));
        assert_eq!(
            translate("(?<=(?>ab|a[bc]))", true),
            r"(?<=(?>\x61\x62|\x61[bc]))"
        );
        // On the engine, the unguarded body keeps `bc`, which Boost does not
        // reach: its first branch matches `b` at 0.
        let unguarded = build(r"(?>[ab]?\x62|[ab]?\x63)(?=)(?=)", 1000, false).unwrap();
        let guarded = build(
            &format!(r"(?>[ab]?\x62|[ab]?\x63|{NO_BYTE})(?=)(?=)"),
            1000,
            false,
        )
        .unwrap();
        let range =
            |regex: &fancy_regex::Regex| regex.find("bc").unwrap().map(|found| found.range());
        assert_eq!(range(&unguarded), Some(0..2));
        assert_eq!(range(&guarded), Some(0..1));
    }

    /// The category of a translation refusal, or `None` when it translates.
    fn refusal(pattern: &str, options: RegexOptions) -> Option<String> {
        match Translator::new(pattern, options, true).run() {
            Ok(_) => None,
            Err(Error::Unsupported(message)) => Some(message),
            Err(other) => panic!("{pattern:?}: {other}"),
        }
    }

    /// Boost's start maps: the case sensitivity of repeat and alternation
    /// states inside a repeated group, `\<` after a case switch, `\<` and `\>`
    /// beside a repeated group with an inner repeat or alternation.
    #[test]
    fn refuses_start_maps_boost_builds_wrongly() {
        let plain = RegexOptions::default();
        let icase = RegexOptions {
            icase: true,
            ..RegexOptions::default()
        };
        let refused = |pattern: &str, options: RegexOptions, category: &str| {
            let message =
                refusal(pattern, options).unwrap_or_else(|| panic!("{pattern:?} translated"));
            assert!(message.contains(category), "{pattern:?}: {message}");
        };
        // A repeat state takes the case sensitivity at its quantifier; the first
        // alternation state stands before a scoped switch, a later one where the
        // previous `|` was.
        for pattern in [
            "(?i:b.+)*C",
            "(?:y|(?i:b.+))+C",
            "(?:(?i)x|y+)*",
            "(?:a(?i)|b|c)+",
            "(?i)(?:(?-i:a|b+)c)+",
        ] {
            refused(pattern, plain, CASE_SWITCHED_REPEAT);
        }
        refused("(?-i:b+)+", icase, CASE_SWITCHED_REPEAT);
        for pattern in [
            "(?i:a|b)*",
            "(?:(?i)a|b)*",
            "(?:(?i:a|b)c+)+",
            "(?:b.+(?i))*C",
            "(?i)(?:b.+)*C",
            "(?:a(?i)|b)+",
            "(?i:b.+)C",
            "(?i)(?:(?-i:a|b)c)+",
        ] {
            assert_eq!(refusal(pattern, plain), None, "{pattern:?}");
        }
        refused(r"a(?i)\<b", plain, WORD_START_CASE_SWITCH);
        refused(r"\<a(?-i)", icase, WORD_START_CASE_SWITCH);
        assert_eq!(refusal(r"(?i)\<a", icase), None);
        assert_eq!(refusal(r"(?i)a\>", plain), None);
        refused(r"(?:\w\w+){2}\>", plain, WORD_BOUNDARY_AFTER_LOOP);
        refused(r"\<(?:a|b)+", plain, WORD_BOUNDARY_AFTER_LOOP);
        assert_eq!(refusal(r"(?:\w\w){2}\>", plain), None);
        assert_eq!(refusal(r"(?:a|b)\>(?:ab)+", plain), None);
    }

    /// The bound on Boost's start-map recursion, at its limit of 100 levels:
    /// anchors count once outside repeated groups and twice inside, as do
    /// alternations inside them and the repeated groups, alternations outside
    /// only with a buffer end, and a margin of 2.
    #[test]
    fn bounds_start_map_recursion() {
        let plain = RegexOptions::default();
        let within_limit = |pattern: &str| match refusal(pattern, plain) {
            None => true,
            Some(message) if message.contains(START_MAP_RECURSION) => false,
            Some(message) => panic!("{pattern:?}: {message}"),
        };
        let alternatives = |count: usize| vec!["b"; count].join("|");
        // (pattern at the bound of 100, the same one site further)
        for (at_limit, beyond) in [
            (
                format!("{}a", "$".repeat(98)),
                format!("{}a", "$".repeat(99)),
            ),
            (
                format!("(?:a{})+", "$".repeat(48)),
                format!("(?:a{})+", "$".repeat(49)),
            ),
            (
                format!("(?:(?:{})x)+", alternatives(49)),
                format!("(?:(?:{})x)+", alternatives(50)),
            ),
            (
                format!("(?:{})\\z", alternatives(50)),
                format!("(?:{})\\z", alternatives(51)),
            ),
            ("(?:a+)+".repeat(49), "(?:a+)+".repeat(50)),
        ] {
            assert!(within_limit(&at_limit), "{at_limit:?}");
            assert!(!within_limit(&beyond), "{beyond:?}");
        }
        // Alternations outside repeated groups count only with a buffer end.
        assert!(within_limit(&format!("(?:{})$", alternatives(200))));
    }

    /// The translator's bound on the group nesting of its spellings refuses a
    /// pattern before the engine's parser would, and not earlier than one level
    /// before: wrapping a pattern in one more group adds one level to both, and
    /// with one wrapper less than the first refused, every program compiles.
    #[test]
    fn bounds_translated_nesting() {
        let options = RegexOptions::default();
        let wrap = |pattern: &str, levels: usize| {
            format!("{}{pattern}{}", "(?:".repeat(levels), ")".repeat(levels))
        };
        // `levels` groups around `inner`, each followed by `quantifier`.
        let nest = |inner: &str, quantifier: &str, levels: usize| {
            format!(
                "{}{inner}{}",
                "(?:".repeat(levels),
                format!("){quantifier}").repeat(levels)
            )
        };
        for pattern in [
            nest("a{2,}", "{2,}", 8),
            nest("a+$", "+$", 20),
            nest("(?=(?:(?=a)b){2,})c", "{2,}", 16),
            nest("(?:(?=a)){0,1}b", "+", 22),
            nest("(?:ab){2,}c?", "+", 16),
        ] {
            let pattern = pattern.as_str();
            let first_refused = (0..MAX_GROUP_DEPTH)
                .find(|&levels| refusal(&wrap(pattern, levels), options).is_some())
                .unwrap_or_else(|| panic!("{pattern:?} is never refused"));
            let message = refusal(&wrap(pattern, first_refused), options).unwrap();
            assert!(
                message.contains(TRANSLATION_TOO_DEEP),
                "{pattern:?}: {message}"
            );
            assert!(first_refused > 0, "{pattern:?}");
            let widest = wrap(pattern, first_refused - 1);
            assert!(BoostRegex::new(&widest).is_ok(), "{widest:?}");
        }
        let deep =
            |levels: usize| format!("{}a{{2,}}{}", "(?:".repeat(levels), "){2,}".repeat(levels));
        assert_eq!(refusal(&deep(19), options), None);
        assert!(refusal(&deep(20), options).is_some_and(|m| m.contains(TRANSLATION_TOO_DEEP)));
    }

    /// `\x` digits that `std::istream` reads with white space, a sign or `0x`.
    #[test]
    fn refuses_hex_escapes_read_as_stream_syntax() {
        for pattern in [
            r"\x{+41}",
            r"\x{ 41}",
            r"\x{0x41}",
            r"\x+4",
            r"[\x-0]",
            r"\x0X",
        ] {
            assert!(
                refusal(pattern, RegexOptions::default())
                    .is_some_and(|message| message.contains(HEX_ESCAPE_STREAM_SYNTAX)),
                "{pattern:?}"
            );
        }
        for pattern in [r"\x41", r"\x{041}", r"\x4x", r"[\x00-\x{7f}]", r"\x0"] {
            assert_eq!(
                refusal(pattern, RegexOptions::default()),
                None,
                "{pattern:?}"
            );
        }
    }

    /// Boost's `probe_leading_repeat` walk, as the translator mirrors it.
    #[test]
    fn finds_leading_lazy_repeats() {
        let leading =
            |pattern: &str| match Translator::new(pattern, RegexOptions::default(), false).run() {
                Ok(_) => false,
                Err(Error::Unsupported(message)) => {
                    assert!(
                        message.contains(LEADING_LAZY_REPEAT),
                        "{pattern:?}: {message}"
                    );
                    true
                }
                Err(error) => panic!("{pattern:?}: {error}"),
            };
        for pattern in [
            r"a{1,3}?\b",
            r"\ba{0,2}?",
            "(?=x)(?<!y)^(?:(a{2,4}?))$",
            "(?#c)()(?s)[ab]{1,3}?(?:x|y)",
        ] {
            assert!(leading(pattern), "{pattern:?}");
        }
        for pattern in [
            r"a{1,2}?\b",
            r"a{1,}?\b",
            r"a{1,3}\b",
            r"x|a{1,3}?\b",
            r"(?:a{1,3}?|x)",
            "(?i)a{1,3}?",
            "(?:a{1,3}?)+",
            "(?=x)?a{1,3}?",
            r"a{1,3}?\b(a)\1",
            "(?=a{1,3}?)",
        ] {
            assert!(!leading(pattern), "{pattern:?}");
        }
    }

    #[test]
    fn detects_translations_that_need_backtracking() {
        let needs = |pattern: &str| {
            needs_backtracking(
                &Translator::new(pattern, RegexOptions::default(), false)
                    .run()
                    .unwrap()
                    .pattern,
            )
        };
        for easy in [r"a{2}[\\1]", r"\Aa\z", r"(?-m:^a$)", r"(a|b)*c", "[^a]"] {
            assert!(!needs(easy), "{easy:?}");
        }
        for hard in [
            r"^a", r"a$", r"\ba", r"a\B", r"\<a\>", r"(a)\1", "(?=a)", "(?!a)", "(?<=a)", "(?>a)",
        ] {
            assert!(needs(hard), "{hard:?}");
        }
    }

    #[test]
    fn counts_and_names_groups() {
        let translation = Translator::new(
            r"cycle=(?<GROUP>\d+)\s+experiment=(?<GROUP>\d+)",
            RegexOptions::default(),
            false,
        )
        .run()
        .unwrap();
        assert_eq!(translation.mark_count, 2);
        assert_eq!(
            translation.names,
            vec![("GROUP".to_string(), 1), ("GROUP".to_string(), 2)]
        );
    }

    #[test]
    fn transcodes_high_bytes_to_single_characters() {
        let haystack = Haystack::new(b"a\xc3\xa9\xff");
        assert_eq!(haystack.len(), 4);
        assert_eq!(haystack.text().chars().count(), 4);
        for original in 0..=4 {
            assert_eq!(haystack.to_original(haystack.to_engine(original)), original);
        }
        assert!(matches!(Haystack::new(b"ascii"), Haystack::Ascii("ascii")));
    }
}
