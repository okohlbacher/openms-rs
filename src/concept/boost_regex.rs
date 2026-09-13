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
//! - `\d`, `\w`, `\s`, `\h`, `\v`, the POSIX classes and case folding are
//!   ASCII, as Boost's `char` traits in the C locale are; they are emitted as
//!   explicit byte ranges, so the engine's own class tables are never consulted;
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
//! to a byte offset. A pure-ASCII haystack is used as it is. Patterns are
//! restricted to ASCII, so a transcoded byte is, for every construct the facade
//! accepts, exactly what the original byte is to Boost: one character that no
//! class, literal or case fold names, a non-word character, and distinct from
//! every other byte for backreferences. The detour exists because
//! `fancy-regex` configures `regex-automata` to refuse empty matches inside a
//! UTF-8 sequence even in its bytes mode, which loses Boost's empty matches
//! between the bytes of a multi-byte character and makes `regex-automata`
//! panic on a lone continuation byte; with a transcoded haystack every byte
//! boundary is a character boundary.
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
use std::sync::Arc;

/// Longest pattern, in bytes, that [`BoostRegex::with_options`] accepts.
pub const MAX_PATTERN_BYTES: usize = 64 * 1024;

/// Deepest group nesting the translator accepts.
///
/// The translation adds up to four levels of its own, and `fancy-regex` stops
/// at 64, so a deeper pattern is refused before it reaches the engine.
pub const MAX_GROUP_DEPTH: usize = 48;

/// Largest repeat bound in `{n}`, `{n,}` or `{n,m}` that the translator reads.
pub const MAX_REPEAT: usize = 999_999_999;

/// Backtracking steps one search may take before it fails, unless
/// [`RegexOptions::backtrack_limit`] says otherwise.
///
/// Only expressions that need the backtracking engine (lookaround,
/// backreferences, atomic groups, line anchors, word boundaries) consume this
/// budget; the others run on an automaton in time linear in the input.
pub const DEFAULT_BACKTRACK_LIMIT: usize = 1_000_000;

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
    /// Backtracking steps one search may take; see [`DEFAULT_BACKTRACK_LIMIT`].
    /// Exceeding it is an error rather than a different answer. Must be
    /// positive.
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

/// A compiled regular expression with `boost::regex` semantics (perl syntax).
///
/// Construction translates the Boost pattern once and compiles three engine
/// programs: one for `regex_search`, one anchored at both ends for
/// `regex_match`, and one that refuses empty matches for the token iterator's
/// `match_not_initial_null` retry.
#[derive(Clone, Debug)]
pub struct BoostRegex {
    pattern: String,
    options: RegexOptions,
    search: fancy_regex::Regex,
    full: fancy_regex::Regex,
    not_empty: Option<fancy_regex::Regex>,
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
        let translation = Translator::new(pattern, options).run()?;
        let search = compile(pattern, &translation.pattern, options, false)?;
        let full = compile(
            pattern,
            &format!(r"\A(?:{})\z", translation.pattern),
            options,
            false,
        )?;
        let not_empty = match build(&translation.pattern, options, true) {
            Ok(regex) => Some(regex),
            Err(fancy_regex::Error::CompileError(error))
                if matches!(*error, fancy_regex::CompileError::PatternCanNeverMatch) =>
            {
                None
            }
            Err(error) => return Err(engine_error(pattern, &error)),
        };
        if search.captures_len() != translation.mark_count + 1
            || full.captures_len() != translation.mark_count + 1
        {
            return Err(Error::Unsupported(format!(
                "regular expression {} compiled to {} groups instead of {}",
                quote(pattern),
                search.captures_len(),
                translation.mark_count + 1
            )));
        }
        Ok(Self {
            pattern: pattern.to_string(),
            options,
            search,
            full,
            not_empty,
            mark_count: translation.mark_count,
            names: translation.names.into(),
        })
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
    /// [`Error::InvalidValue`] when the search exceeds the backtrack limit,
    /// where Boost would throw `regex_error` with `error_complexity`. The two
    /// limits are counted differently, so the inputs that hit them differ.
    pub fn search(&self, haystack: &[u8]) -> Result<Option<Captures>> {
        self.search_from(&Haystack::new(haystack), 0, false)
    }

    /// `boost::regex_search(haystack, regex)` without match results.
    ///
    /// # Errors
    ///
    /// As [`BoostRegex::search`].
    pub fn is_search_match(&self, haystack: &[u8]) -> Result<bool> {
        self.search
            .is_match(Haystack::new(haystack).text())
            .map_err(|error| self.runtime_error(&error))
    }

    /// `boost::regex_match(haystack, match, regex)`: a match of the whole
    /// haystack, backtracking into alternatives until one spans it.
    ///
    /// # Errors
    ///
    /// As [`BoostRegex::search`].
    pub fn full_match(&self, haystack: &[u8]) -> Result<Option<Captures>> {
        let haystack = Haystack::new(haystack);
        self.captures(&self.full, &haystack, RegexInput::new(haystack.text()))
    }

    /// `boost::regex_match(haystack, regex)` without match results.
    ///
    /// # Errors
    ///
    /// As [`BoostRegex::search`].
    pub fn is_full_match(&self, haystack: &[u8]) -> Result<bool> {
        self.full
            .is_match(Haystack::new(haystack).text())
            .map_err(|error| self.runtime_error(&error))
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
    /// exceeds the backtrack limit; the iterator ends after an error.
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
        if !not_initial_null {
            let input = RegexInput::new(text).from_pos(haystack.to_engine(start));
            return self.captures(&self.search, haystack, input);
        }
        if let Some(not_empty) = &self.not_empty {
            let input = RegexInput::new(text)
                .from_pos(haystack.to_engine(start))
                .anchored(true);
            if let Some(captures) = self.captures(not_empty, haystack, input)? {
                return Ok(Some(captures));
            }
        }
        if start >= haystack.len() {
            return Ok(None);
        }
        let input = RegexInput::new(text).from_pos(haystack.to_engine(start + 1));
        self.captures(&self.search, haystack, input)
    }

    fn captures(
        &self,
        regex: &fancy_regex::Regex,
        haystack: &Haystack<'_>,
        input: RegexInput<'_, str>,
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
            Err(error) => Err(self.runtime_error(&error)),
        }
    }

    fn runtime_error(&self, error: &fancy_regex::Error) -> Error {
        match error {
            fancy_regex::Error::RuntimeError(fancy_regex::RuntimeError::BacktrackLimitExceeded) => {
                Error::InvalidValue(format!(
                    "regular expression {} exceeded the backtracking limit of {} steps",
                    quote(&self.pattern),
                    self.options.backtrack_limit
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
    options: RegexOptions,
    not_empty: bool,
) -> std::result::Result<fancy_regex::Regex, fancy_regex::Error> {
    let mut builder = RegexBuilder::new(translated);
    builder
        .bytes_mode(BytesMode::Unicode)
        .unicode_mode(true)
        .backtrack_limit(options.backtrack_limit)
        .find_not_empty(not_empty);
    builder.build()
}

fn compile(
    pattern: &str,
    translated: &str,
    options: RegexOptions,
    not_empty: bool,
) -> Result<fancy_regex::Regex> {
    build(translated, options, not_empty).map_err(|error| engine_error(pattern, &error))
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

/// ASCII members of a Boost character class, spelled for a `regex-syntax`
/// bracket expression.
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
    fn members(self) -> &'static str {
        match self {
            Self::Alnum => "0-9A-Za-z",
            Self::Alpha => "A-Za-z",
            // isspace && !is_separator: space, tab and vertical tab.
            Self::Blank => r"\x09\x0B\x20",
            Self::Cntrl => r"\x00-\x1F\x7F",
            Self::Digit => "0-9",
            Self::Graph => r"\x21-\x7E",
            // isspace && !is_separator && c != '\v'.
            Self::Horizontal => r"\x09\x20",
            Self::Lower => "a-z",
            Self::Print => r"\x20-\x7E",
            Self::Punct => r"\x21-\x2F\x3A-\x40\x5B-\x60\x7B-\x7E",
            Self::Space => r"\x09-\x0D\x20",
            Self::Upper => "A-Z",
            // is_separator || c == '\v'.
            Self::Vertical => r"\x0A-\x0D",
            Self::Word => "0-9A-Za-z_",
            Self::Xdigit => "0-9A-Fa-f",
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

#[derive(Clone, Copy, Debug)]
enum ClassItem {
    Byte(u8),
    Range(u8, u8),
    Class(ClassKind, bool),
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
    /// A backreference.
    Backref,
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

#[derive(Debug)]
struct Frame {
    kind: GroupKind,
    saved_flags: Flags,
    start: usize,
    /// Width of the current alternative, `None` once Boost's `calculate_backstep`
    /// would give up.
    width: Option<usize>,
    /// Width of the first finished alternative: `None` before any, `Some(None)`
    /// once alternatives disagree or one has no width.
    alternative_width: Option<Option<usize>>,
    /// Shortest match of the current alternative.
    min_length: usize,
    /// Shortest match over the finished alternatives.
    alternative_min_length: Option<usize>,
    /// Whether Boost has appended any state inside the group (comments append none).
    has_states: bool,
    alternation: bool,
    items: usize,
    last_item_is_lookaround: bool,
    consumes: bool,
    asserts: bool,
    has_capture: bool,
}

impl Frame {
    fn new(kind: GroupKind, saved_flags: Flags, start: usize) -> Self {
        Self {
            kind,
            saved_flags,
            start,
            width: Some(0),
            alternative_width: None,
            min_length: 0,
            alternative_min_length: None,
            has_states: false,
            alternation: false,
            items: 0,
            last_item_is_lookaround: false,
            consumes: false,
            asserts: false,
            has_capture: false,
        }
    }

    fn finish_alternative(&mut self) {
        self.alternative_width = Some(match self.alternative_width {
            None => self.width,
            Some(previous) if previous == self.width => previous,
            Some(_) => None,
        });
        self.width = Some(0);
        self.alternative_min_length = Some(
            self.alternative_min_length
                .map_or(self.min_length, |previous| previous.min(self.min_length)),
        );
        self.min_length = 0;
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
    last: Last,
    /// Shortest match of `last`, removed again when a quantifier rescales it.
    last_min_length: usize,
}

impl<'p> Translator<'p> {
    fn new(source: &'p str, options: RegexOptions) -> Self {
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
            last: Last::Nothing,
            last_min_length: 0,
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
                        END_OF_LINE
                    } else {
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
        Ok(Translation {
            pattern: self.out,
            mark_count: self.mark_count,
            names: self.names,
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

    fn add_min_length(&mut self, length: usize) {
        let frame = self.frame();
        frame.min_length = frame.min_length.saturating_add(length);
        self.last_min_length = length;
    }

    fn literal(&mut self, byte: u8) {
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
        self.add_min_length(1);
        self.frame().consumes = true;
        self.last = Last::Byte;
    }

    fn assertion(&mut self, text: &str) {
        self.out.push_str(text);
        self.add_item(false);
        self.add_min_length(0);
        self.frame().asserts = true;
        self.last = Last::Fixed;
    }

    fn emit_class(&mut self, negated: bool, items: &[ClassItem]) {
        let mut class = String::from(if negated { "[^" } else { "[" });
        for item in items {
            match *item {
                ClassItem::Byte(byte) => {
                    let _ = write!(class, r"\x{byte:02X}");
                }
                ClassItem::Range(first, last) => {
                    let _ = write!(class, r"\x{first:02X}-\x{last:02X}");
                }
                ClassItem::Class(kind, false) => class.push_str(kind.members()),
                ClassItem::Class(kind, true) => {
                    let _ = write!(class, "[^{}]", kind.members());
                }
            }
        }
        class.push(']');
        if self.flags.icase {
            let _ = write!(self.out, "(?i:{class})");
        } else {
            self.out.push_str(&class);
        }
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
            b'<' => self.assertion(r"\<"),
            b'>' => self.assertion(r"\>"),
            b'A' | b'`' => self.assertion(r"\A"),
            b'z' | b'\'' => self.assertion(r"\z"),
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
                self.max_backreference = self.max_backreference.max(group);
                if self.flags.icase {
                    let _ = write!(self.out, r"(?i:\{group})");
                } else {
                    let _ = write!(self.out, r"\{group}");
                }
                self.add_item(false);
                self.clear_width();
                self.add_min_length(0);
                self.frame().consumes = true;
                self.last = Last::Backref;
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
    fn hex_escape(&mut self, start: usize) -> Result<u8> {
        let Some(&next) = self.bytes.get(self.position) else {
            return Err(self.syntax(start, "hexadecimal escape sequence terminated prematurely"));
        };
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
                if first > last {
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
        self.push_group(GroupKind::Capture, "(", self.flags);
        Ok(())
    }

    fn push_group(&mut self, kind: GroupKind, text: &str, saved_flags: Flags) {
        let start = self.out.len();
        self.out.push_str(text);
        self.frame().has_states = true;
        self.frames.push(Frame::new(kind, saved_flags, start));
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
                self.flags = flags;
                self.frame().has_states = true;
                self.last = Last::EmptyGroup { case_change };
                self.last_min_length = 0;
                Ok(())
            }
            Some(b':') => {
                self.position += 1;
                let saved = self.flags;
                self.flags = flags;
                self.push_group(GroupKind::NonCapture, "(?:", saved);
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
        let min_length = frame.alternative_min_length.unwrap_or(0);
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
        self.out.push(')');
        self.flags = frame.saved_flags;
        let body = frame.body();
        let captures = frame.has_capture || frame.kind == GroupKind::Capture;
        if frame.kind.is_lookaround() {
            self.add_item(true);
            self.add_min_length(0);
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
            let parent = self.frame();
            parent.consumes |= frame.consumes;
            parent.asserts |= frame.asserts;
            parent.has_capture |= captures;
        }
        self.last = Last::Group(GroupShape {
            start: frame.start,
            kind: frame.kind,
            body,
            captures,
            nullable: min_length == 0,
            only_asserts: frame.asserts && !frame.consumes,
        });
        Ok(())
    }

    fn alternation(&mut self) {
        self.position += 1;
        self.out.push('|');
        let frame = self.frame();
        frame.finish_alternative();
        frame.alternation = true;
        frame.has_states = true;
        self.last = Last::Nothing;
        self.last_min_length = 0;
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
        let Some(min) = self.decimal(start, &mut cursor)? else {
            return self.literal_brace();
        };
        cursor = skip_space(self.bytes, cursor);
        let max = match self.bytes.get(cursor) {
            None => return self.literal_brace(),
            Some(b',') => {
                cursor = skip_space(self.bytes, cursor + 1);
                if cursor >= self.bytes.len() {
                    return self.literal_brace();
                }
                self.decimal(start, &mut cursor)?
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
        self.position = cursor + 1;
        self.quantify(start, min, max)
    }

    fn decimal(&self, start: usize, cursor: &mut usize) -> Result<Option<usize>> {
        let digits = self.bytes[*cursor..]
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count();
        if digits == 0 {
            return Ok(None);
        }
        if digits > 9 {
            return Err(self.unsupported(start, "a repeat bound above MAX_REPEAT"));
        }
        let value = self.source[*cursor..*cursor + digits]
            .parse::<usize>()
            .unwrap_or(usize::MAX);
        *cursor += digits;
        Ok(Some(value))
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
        match self.last {
            Last::Nothing => return Err(self.syntax(start, "nothing to repeat")),
            Last::Fixed => {
                return Err(self.syntax(
                    start,
                    "a repeat cannot apply to a zero-width assertion or another repeat",
                ));
            }
            Last::EmptyGroup { case_change: true } => {
                return Err(self.unsupported(
                    start,
                    "a repeat of a modifier group that switches case sensitivity (Boost undoes \
                     the switch when the empty repetition is abandoned)",
                ));
            }
            Last::Group(shape)
                if shape.captures && shape.nullable && max.is_none_or(|max| max > 1) =>
            {
                return Err(self.unsupported(
                    start,
                    "a repeat of a group that can match the empty string and captures (Boost \
                     records a final empty iteration)",
                ));
            }
            Last::Group(shape) if shape.only_asserts => {
                return Err(
                    self.unsupported(start, "a repeat of a group that only asserts a position")
                );
            }
            Last::Byte => {
                self.out.push_str(&suffix);
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
            Last::Backref => {
                self.out.push_str(&suffix);
                self.clear_width();
            }
            Last::EmptyGroup { case_change: false } => self.clear_width(),
            Last::Group(shape) => {
                self.clear_width();
                let zero_width = shape.kind.is_lookaround()
                    || (shape.kind == GroupKind::NonCapture && shape.body == Body::LookAround);
                if shape.kind == GroupKind::NonCapture && shape.body == Body::Empty {
                    // The engine has no quantifiable empty group; an empty group
                    // repeated any number of times is still empty.
                } else if zero_width {
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
                    self.out.push_str(&suffix);
                }
            }
        }
        let contribution = self.last_min_length;
        let frame = self.frame();
        frame.min_length = frame
            .min_length
            .saturating_sub(contribution)
            .saturating_add(contribution.saturating_mul(min));
        self.last = Last::Fixed;
        self.last_min_length = 0;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translates_line_anchors_and_dot() {
        let translation = Translator::new("^a.$", RegexOptions::default())
            .run()
            .unwrap();
        assert_eq!(
            translation.pattern,
            format!(r"{START_OF_LINE}\x61(?s:.){END_OF_LINE}")
        );
    }

    #[test]
    fn counts_and_names_groups() {
        let translation = Translator::new(
            r"cycle=(?<GROUP>\d+)\s+experiment=(?<GROUP>\d+)",
            RegexOptions::default(),
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
