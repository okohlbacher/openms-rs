// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::*;
use std::{fmt, mem::size_of};

/// Maximum input bytes for either text grammar.
pub const MAX_PROFORMA_PARSE_INPUT_BYTES: usize = 4 * 1024 * 1024;
/// Maximum appended AST collection items, shared across every chain.
pub const MAX_PROFORMA_PARSE_NODES: usize = 1_000_000;
/// Maximum charged scanning, lookahead, copying and traversal work.
pub const MAX_PROFORMA_PARSE_WORK: usize = 50_000_000;
/// Maximum cumulative requested AST and temporary buffer capacity in bytes.
pub const MAX_PROFORMA_PARSE_BYTES: usize = 256 * 1024 * 1024;
/// Maximum bracket-lookahead or nested-name depth.
pub const MAX_PROFORMA_PARSE_DEPTH: usize = 256;
/// Maximum combined diagnostic message, expected and found bytes.
pub const MAX_PROFORMA_DIAGNOSTIC_BYTES: usize = 64 * 1024;
type Parsed<T> = std::result::Result<T, ParseFailure>;

/// Complete source parse-error vocabulary; several codes are not emitted by
/// source parsing itself and remain available for constructed diagnostics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ErrorCode {
    UnexpectedCharacter,
    UnclosedBracket,
    UnmatchedBracket,
    InvalidCvPrefix,
    InvalidCvAccession,
    InvalidAminoAcid,
    InvalidMassValue,
    InvalidFormula,
    UnknownMonosaccharide,
    DanglingCrosslinkLabel,
    EmptySequence,
    InvalidCharge,
    InvalidOccurrenceSpecifier,
    UnexpectedEndOfInput,
    InternalError,
}
impl ErrorCode {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::UnexpectedCharacter => "Unexpected character",
            Self::UnclosedBracket => "Unclosed bracket",
            Self::UnmatchedBracket => "Unmatched closing bracket",
            Self::InvalidCvPrefix => "Invalid controlled vocabulary prefix",
            Self::InvalidCvAccession => "Invalid CV accession number",
            Self::InvalidAminoAcid => "Invalid amino acid",
            Self::InvalidMassValue => "Invalid mass value",
            Self::InvalidFormula => "Invalid chemical formula",
            Self::UnknownMonosaccharide => "Unknown monosaccharide",
            Self::DanglingCrosslinkLabel => "Dangling crosslink label",
            Self::EmptySequence => "Empty sequence",
            Self::InvalidCharge => "Invalid charge state",
            Self::InvalidOccurrenceSpecifier => "Invalid occurrence specifier",
            Self::UnexpectedEndOfInput => "Unexpected end of input",
            Self::InternalError => "Internal parser error",
        }
    }
}

/// Owned source diagnostic. Context accessors expose exact byte windows;
/// human-readable formatting replaces only invalid split UTF-8 with U+FFFD.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParseError {
    code: ErrorCode,
    position: usize,
    before: Vec<u8>,
    after: Vec<u8>,
    message: String,
    expected: String,
    found: String,
}
impl ParseError {
    pub fn new(code: ErrorCode, position: usize, input: &str, message: &str) -> Result<Self> {
        if message.len() > MAX_PROFORMA_DIAGNOSTIC_BYTES {
            return Err(invalid("ProForma diagnostic size limit exceeded"));
        }
        let position = position.min(input.len());
        let before = copy_bytes(&input.as_bytes()[position.saturating_sub(20)..position])?;
        let after =
            copy_bytes(&input.as_bytes()[position..position.saturating_add(20).min(input.len())])?;
        Ok(Self {
            code,
            position,
            before,
            after,
            message: copy_string(message)?,
            expected: String::new(),
            found: String::new(),
        })
    }
    pub fn code(&self) -> ErrorCode {
        self.code
    }
    pub fn position(&self) -> usize {
        self.position
    }
    pub fn context_before(&self) -> &[u8] {
        &self.before
    }
    pub fn context_after(&self) -> &[u8] {
        &self.after
    }
    pub fn message(&self) -> &str {
        &self.message
    }
    pub fn expected(&self) -> &str {
        &self.expected
    }
    pub fn found(&self) -> &str {
        &self.found
    }
    pub fn set_expected_found(&mut self, expected: &str, found: &str) -> Result<()> {
        let length = self
            .message
            .len()
            .checked_add(expected.len())
            .and_then(|n| n.checked_add(found.len()));
        if length.is_none_or(|n| n > MAX_PROFORMA_DIAGNOSTIC_BYTES) {
            return Err(invalid("ProForma diagnostic size limit exceeded"));
        }
        let expected = copy_string(expected)?;
        let found = copy_string(found)?;
        self.expected = expected;
        self.found = found;
        Ok(())
    }
    pub fn formatted_message(&self) -> Result<String> {
        use fmt::Write;
        let mut output = String::new();
        output
            .try_reserve_exact(self.expected.len() + self.found.len() + 512)
            .map_err(|_| invalid("cannot allocate ProForma diagnostic"))?;
        write!(output, "{self}").map_err(|_| invalid("cannot format ProForma diagnostic"))?;
        Ok(output)
    }
}
impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "ProForma parse error at position {}: {}\nContext: ",
            self.position,
            self.code.as_str()
        )?;
        if self.position > self.before.len() {
            f.write_str("...")?;
        }
        f.write_str(&String::from_utf8_lossy(&self.before))?;
        if let Some((first, rest)) = self.after.split_first() {
            write!(
                f,
                ">>>{}<<<{}",
                String::from_utf8_lossy(std::slice::from_ref(first)),
                String::from_utf8_lossy(rest)
            )?;
        } else {
            f.write_str(">>><END OF INPUT><<<")?;
        }
        if self.after.len() >= 20 {
            f.write_str("...")?;
        }
        if !self.expected.is_empty() {
            write!(f, "\nExpected: {}", self.expected)?;
        }
        if !self.found.is_empty() {
            write!(f, "\nFound: {}", self.found)?;
        }
        Ok(())
    }
}
impl std::error::Error for ParseError {}

/// Syntax diagnostics remain distinct from checked work/allocation failures.
#[derive(Debug)]
pub enum ParseFailure {
    Syntax(Box<ParseError>),
    Resource(Error),
}
impl fmt::Display for ParseFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Syntax(error) => error.fmt(f),
            Self::Resource(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for ParseFailure {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Syntax(error) => error.as_ref(),
            Self::Resource(error) => error,
        })
    }
}
impl From<Error> for ParseFailure {
    fn from(error: Error) -> Self {
        Self::Resource(error)
    }
}
fn copy_bytes(value: &[u8]) -> Result<Vec<u8>> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(value.len())
        .map_err(|_| invalid("cannot allocate ProForma diagnostic"))?;
    output.extend_from_slice(value);
    Ok(output)
}
fn copy_string(value: &str) -> Result<String> {
    String::from_utf8(copy_bytes(value.as_bytes())?)
        .map_err(|_| invalid("invalid ProForma diagnostic UTF-8"))
}

impl Peptidoform {
    /// Parse the complete source single-chain grammar, without resolution.
    pub fn parse(input: &str) -> Parsed<Self> {
        let mut parser = Parser::new(input)?;
        let result = parser.chain()?;
        parser.require_end("Unexpected characters after peptidoform")?;
        Ok(result)
    }
}
impl PeptidoformIon {
    /// Parse the complete source ion grammar, including mixed chain separators.
    /// Labels are retained without implicit cross-link partner validation.
    pub fn parse(input: &str) -> Parsed<Self> {
        Parser::new(input)?.ion()
    }
}

struct Work {
    remaining: usize,
    bytes: usize,
    nodes: usize,
}
impl Work {
    fn new() -> Self {
        Self {
            remaining: MAX_PROFORMA_PARSE_WORK,
            bytes: MAX_PROFORMA_PARSE_BYTES,
            nodes: MAX_PROFORMA_PARSE_NODES,
        }
    }
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| invalid("ProForma parse work limit exceeded"))?;
        Ok(())
    }
    fn allocate(&mut self, count: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(count)
            .ok_or_else(|| invalid("ProForma parse allocation limit exceeded"))?;
        Ok(())
    }
    fn reserve<T>(&mut self, values: &mut Vec<T>, wanted: usize) -> Result<()> {
        if wanted > values.capacity() {
            self.consume(values.len())?;
            let capacity = wanted.max(values.capacity().saturating_mul(2)).max(4);
            self.allocate(
                capacity
                    .checked_mul(size_of::<T>())
                    .ok_or_else(|| invalid("ProForma parse size overflow"))?,
            )?;
            values
                .try_reserve_exact(capacity - values.len())
                .map_err(|_| invalid("cannot allocate ProForma AST"))?;
        }
        Ok(())
    }
    fn push<T>(&mut self, values: &mut Vec<T>, value: T) -> Result<()> {
        self.nodes = self
            .nodes
            .checked_sub(1)
            .ok_or_else(|| invalid("ProForma parse node limit exceeded"))?;
        self.consume(1)?;
        self.reserve(
            values,
            values
                .len()
                .checked_add(1)
                .ok_or_else(|| invalid("ProForma parse size overflow"))?,
        )?;
        values.push(value);
        Ok(())
    }
    fn append(&mut self, values: &mut Vec<u8>, bytes: &[u8]) -> Result<()> {
        self.consume(bytes.len())?;
        self.reserve(
            values,
            values
                .len()
                .checked_add(bytes.len())
                .ok_or_else(|| invalid("ProForma parse size overflow"))?,
        )?;
        values.extend_from_slice(bytes);
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Lb,
    Rb,
    Lp,
    Rp,
    Lc,
    Rc,
    La,
    Ra,
    Plus,
    Minus,
    Slash,
    Pipe,
    Hash,
    Colon,
    Comma,
    Caret,
    Question,
    At,
    Number,
    Id,
    End,
}
#[derive(Clone, Copy)]
struct Token {
    kind: Kind,
    start: usize,
    end: usize,
}
#[derive(Clone)]
struct Lexer {
    position: usize,
    peek: Option<Token>,
}
impl Lexer {
    fn new(position: usize) -> Self {
        Self {
            position,
            peek: None,
        }
    }
    fn peek(&mut self, input: &[u8], work: &mut Work) -> Result<Token> {
        work.consume(1)?;
        if self.peek.is_none() {
            self.peek = Some(self.scan(input, work)?);
        }
        Ok(self.peek.unwrap())
    }
    fn next(&mut self, input: &[u8], work: &mut Work) -> Result<Token> {
        work.consume(1)?;
        if let Some(token) = self.peek.take() {
            Ok(token)
        } else {
            self.scan(input, work)
        }
    }
    fn bump(&mut self, work: &mut Work) -> Result<()> {
        work.consume(2)?;
        self.position += 1;
        Ok(())
    }
    fn scan(&mut self, input: &[u8], work: &mut Work) -> Result<Token> {
        let start = self.position;
        let Some(&c) = input.get(start) else {
            return Ok(Token {
                kind: Kind::End,
                start,
                end: start,
            });
        };
        let punctuation = match c {
            b'[' => Some(Kind::Lb),
            b']' => Some(Kind::Rb),
            b'(' => Some(Kind::Lp),
            b')' => Some(Kind::Rp),
            b'{' => Some(Kind::Lc),
            b'}' => Some(Kind::Rc),
            b'<' => Some(Kind::La),
            b'>' => Some(Kind::Ra),
            b'/' => Some(Kind::Slash),
            b'|' => Some(Kind::Pipe),
            b'#' => Some(Kind::Hash),
            b':' => Some(Kind::Colon),
            b',' => Some(Kind::Comma),
            b'^' => Some(Kind::Caret),
            b'?' => Some(Kind::Question),
            b'@' => Some(Kind::At),
            _ => None,
        };
        let kind = if let Some(kind) = punctuation {
            self.bump(work)?;
            kind
        } else if (c == b'+' || c == b'-')
            && !input
                .get(start + 1)
                .is_some_and(|c| c.is_ascii_digit() || *c == b'.')
        {
            self.bump(work)?;
            if c == b'+' { Kind::Plus } else { Kind::Minus }
        } else if c.is_ascii_digit()
            || c == b'+'
            || c == b'-'
            || (c == b'.' && input.get(start + 1).is_some_and(u8::is_ascii_digit))
        {
            if c == b'+' || c == b'-' {
                self.bump(work)?;
            }
            while input.get(self.position).is_some_and(u8::is_ascii_digit) {
                self.bump(work)?;
            }
            if input.get(self.position) == Some(&b'.')
                && input.get(self.position + 1).is_some_and(u8::is_ascii_digit)
            {
                self.bump(work)?;
                while input.get(self.position).is_some_and(u8::is_ascii_digit) {
                    self.bump(work)?;
                }
            }
            Kind::Number
        } else if c.is_ascii_alphabetic() {
            while input
                .get(self.position)
                .is_some_and(u8::is_ascii_alphabetic)
            {
                self.bump(work)?;
            }
            Kind::Id
        } else {
            self.bump(work)?;
            Kind::Id
        };
        Ok(Token {
            kind,
            start,
            end: self.position,
        })
    }
}

struct Parser<'a> {
    input: &'a str,
    lexer: Lexer,
    current: Option<Token>,
    work: Work,
}
impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Parsed<Self> {
        if input.len() > MAX_PROFORMA_PARSE_INPUT_BYTES {
            return Err(invalid("ProForma parse input limit exceeded").into());
        }
        let mut work = Work::new();
        work.consume(input.len())?;
        Ok(Self {
            input,
            lexer: Lexer::new(0),
            current: None,
            work,
        })
    }
    fn current(&mut self) -> Parsed<Token> {
        self.work.consume(1)?;
        if self.current.is_none() {
            self.current = Some(self.lexer.next(self.input.as_bytes(), &mut self.work)?);
        }
        Ok(self.current.unwrap())
    }
    fn advance(&mut self) -> Parsed<Token> {
        let token = self.current()?;
        self.current = None;
        Ok(token)
    }
    fn check(&mut self, kind: Kind) -> Parsed<bool> {
        Ok(self.current()?.kind == kind)
    }
    fn take(&mut self, kind: Kind) -> Parsed<bool> {
        if self.check(kind)? {
            self.advance()?;
            Ok(true)
        } else {
            Ok(false)
        }
    }
    fn expect(&mut self, kind: Kind, description: &'static str) -> Parsed<Token> {
        let token = self.current()?;
        if token.kind != kind {
            return Err(self.syntax(
                ErrorCode::UnexpectedCharacter,
                token.start,
                &format!("Expected {description}"),
            ));
        }
        self.advance()
    }
    fn require_end(&mut self, message: &'static str) -> Parsed<()> {
        if !self.check(Kind::End)? {
            return self.fail(ErrorCode::UnexpectedCharacter, message);
        }
        Ok(())
    }
    fn syntax(&mut self, code: ErrorCode, position: usize, message: &str) -> ParseFailure {
        let measured = self
            .work
            .consume(message.len() + 256)
            .and_then(|_| self.work.allocate(message.len() + 512));
        match measured.and_then(|_| ParseError::new(code, position, self.input, message)) {
            Ok(error) => ParseFailure::Syntax(Box::new(error)),
            Err(error) => ParseFailure::Resource(error),
        }
    }
    fn fail<T>(&mut self, code: ErrorCode, message: &'static str) -> Parsed<T> {
        let position = self.current()?.start;
        Err(self.syntax(code, position, message))
    }
    fn text(&self, token: Token) -> &'a [u8] {
        &self.input.as_bytes()[token.start..token.end]
    }
    fn is_text(&self, token: Token, text: &[u8]) -> bool {
        self.text(token) == text
    }
    fn finish(&mut self, bytes: Vec<u8>, position: usize) -> Parsed<String> {
        self.work.consume(bytes.len())?;
        String::from_utf8(bytes).map_err(|_| {
            self.syntax(
                ErrorCode::UnexpectedCharacter,
                position,
                "Source token would produce an invalid UTF-8 annotation field",
            )
        })
    }
    fn string(&mut self, token: Token) -> Parsed<String> {
        let mut bytes = Vec::new();
        self.append_token(&mut bytes, token)?;
        self.finish(bytes, token.start)
    }
    fn append_token(&mut self, bytes: &mut Vec<u8>, token: Token) -> Parsed<()> {
        let text = self.text(token);
        self.work.append(bytes, text)?;
        Ok(())
    }
    fn lookahead(&self) -> Lexer {
        Lexer::new(self.current.map_or(self.lexer.position, |t| t.start))
    }
    fn next(&mut self, lexer: &mut Lexer) -> Parsed<Token> {
        Ok(lexer.next(self.input.as_bytes(), &mut self.work)?)
    }
    fn peek(&mut self, lexer: &mut Lexer) -> Parsed<Token> {
        Ok(lexer.peek(self.input.as_bytes(), &mut self.work)?)
    }
    fn increment_depth(depth: &mut usize) -> Result<()> {
        if *depth >= MAX_PROFORMA_PARSE_DEPTH {
            return Err(invalid("ProForma parse nesting limit exceeded"));
        }
        *depth += 1;
        Ok(())
    }
    fn skip_bracket(&mut self, lexer: &mut Lexer) -> Parsed<bool> {
        self.next(lexer)?;
        let mut depth = 1;
        while depth > 0 {
            let token = self.next(lexer)?;
            match token.kind {
                Kind::Lb => Self::increment_depth(&mut depth)?,
                Kind::Rb => depth -= 1,
                Kind::End => return Ok(false),
                _ => (),
            }
        }
        Ok(true)
    }
    fn integer(&mut self, token: Token, code: ErrorCode, message: &'static str) -> Parsed<i32> {
        self.work.consume(token.end - token.start)?;
        let bytes = self.text(token);
        let mut end = usize::from(bytes.first().is_some_and(|x| *x == b'+' || *x == b'-'));
        while bytes.get(end).is_some_and(u8::is_ascii_digit) {
            end += 1;
        }
        std::str::from_utf8(&bytes[..end])
            .ok()
            .and_then(|s| s.parse::<i32>().ok())
            .ok_or_else(|| self.syntax(code, token.start, message))
    }
    fn signed_integer(&mut self, message: &'static str) -> Parsed<i32> {
        let sign = if self.take(Kind::Plus)? {
            1
        } else if self.take(Kind::Minus)? {
            -1
        } else {
            1
        };
        let token = self.expect(Kind::Number, "charge value")?;
        let value = self.integer(token, ErrorCode::InvalidCharge, message)?;
        value
            .checked_mul(sign)
            .ok_or_else(|| self.syntax(ErrorCode::InvalidCharge, token.start, message))
    }
    fn float(&mut self, bytes: &[u8], position: usize, message: &'static str) -> Parsed<f64> {
        self.work.consume(bytes.len().saturating_mul(2))?;
        let value = std::str::from_utf8(bytes)
            .ok()
            .and_then(|s| s.parse::<f64>().ok());
        match value {
            Some(value)
                if value.is_finite()
                    && (value != 0.0 || !bytes.iter().any(|b| matches!(b, b'1'..=b'9'))) =>
            {
                Ok(value)
            }
            _ => Err(self.syntax(ErrorCode::InvalidMassValue, position, message)),
        }
    }
}

impl Parser<'_> {
    fn ion(&mut self) -> Parsed<PeptidoformIon> {
        let mut ion = PeptidoformIon::default();
        let chain = self.chain_with_charge(false)?;
        self.work.push(&mut ion.chains, chain)?;
        while !self.check(Kind::End)? {
            if self.take(Kind::Slash)? {
                if !self.take(Kind::Slash)? {
                    ion.charge = Some(self.charge()?);
                    break;
                }
                let chain = self.chain_with_charge(false)?;
                self.work.push(&mut ion.chains, chain)?;
            } else if self.take(Kind::Plus)? {
                ion.is_chimeric = true;
                let chain = self.chain_with_charge(true)?;
                self.work.push(&mut ion.chains, chain)?;
            } else {
                break;
            }
        }
        self.require_end("Unexpected characters after peptidoform ion")?;
        Ok(ion)
    }
    fn chain_with_charge(&mut self, chimeric: bool) -> Parsed<Peptidoform> {
        let mut chain = self.chain()?;
        if self.check(Kind::Slash)? {
            let mut look = self.lookahead();
            self.next(&mut look)?;
            let next = self.peek(&mut look)?;
            if next.kind == Kind::Slash {
                return Ok(chain);
            }
            if matches!(
                next.kind,
                Kind::Plus | Kind::Minus | Kind::Number | Kind::Lb
            ) {
                if next.kind == Kind::Lb {
                    self.skip_bracket(&mut look)?;
                } else {
                    if matches!(next.kind, Kind::Plus | Kind::Minus) {
                        self.next(&mut look)?;
                    }
                    if self.peek(&mut look)?.kind == Kind::Number {
                        self.next(&mut look)?;
                    }
                }
                let after = self.peek(&mut look)?.kind;
                if after == Kind::Plus || (after == Kind::End && chimeric) {
                    self.advance()?;
                    chain.charge = Some(self.charge()?);
                }
            }
        }
        Ok(chain)
    }
    fn chain(&mut self) -> Parsed<Peptidoform> {
        let mut chain = Peptidoform::default();
        let mut names = Vec::new();
        let name_position = self.current()?.start;
        while self.check(Kind::Lp)? {
            let mut look = self.lookahead();
            self.next(&mut look)?;
            if self.peek(&mut look)?.kind != Kind::Ra {
                break;
            }
            self.advance()?;
            while self.take(Kind::Ra)? {}
            let mut name = Vec::new();
            let mut depth = 0;
            while !self.check(Kind::End)? {
                let token = self.current()?;
                match token.kind {
                    Kind::Lp => Self::increment_depth(&mut depth)?,
                    Kind::Rp => {
                        if depth == 0 {
                            break;
                        }
                        depth -= 1;
                    }
                    _ => (),
                }
                self.append_token(&mut name, token)?;
                self.advance()?;
            }
            self.expect(Kind::Rp, "')' to close gene/protein prefix")?;
            if !names.is_empty() {
                self.work.append(&mut names, b" / ")?;
            }
            self.work.append(&mut names, &name)?;
        }
        if !names.is_empty() {
            chain.name = Some(self.finish(names, name_position)?);
        }
        while self.take(Kind::La)? {
            let entry = self.global_entry()?;
            self.work.push(&mut chain.global_mods, entry)?;
            while self.take(Kind::Comma)? {
                let entry = self.global_entry()?;
                self.work.push(&mut chain.global_mods, entry)?;
            }
            self.expect(Kind::Ra, "'>' to close global modifications")?;
        }
        chain.unlocalised_mods = self.unlocalised()?;
        while self.take(Kind::Lc)? {
            let modification = self.modification()?;
            self.work.consume(modification.alternatives.len())?;
            if modification
                .alternatives
                .iter()
                .any(|(_, label)| label.is_some())
            {
                return self.fail(
                    ErrorCode::UnexpectedCharacter,
                    "Labels are not allowed on labile modifications",
                );
            }
            self.expect(Kind::Rc, "'}'")?;
            self.work
                .push(&mut chain.labile_mods, LabileModification { modification })?;
        }
        if self.check(Kind::Lb)? && self.n_terminal_pattern()? {
            chain.n_term_mods = self.modifications()?;
            self.expect(Kind::Minus, "'-' after N-terminal modification")?;
        }
        chain.sequence = self.sequence()?;
        if chain.sequence.is_empty() {
            return self.fail(ErrorCode::EmptySequence, "Empty sequence");
        }
        if self.take(Kind::Minus)? {
            chain.c_term_mods = self.modifications()?;
            if chain.c_term_mods.is_empty() {
                return self.fail(
                    ErrorCode::UnexpectedCharacter,
                    "Expected modification after '-' for C-terminal modification",
                );
            }
        }
        Ok(chain)
    }
    fn n_terminal_pattern(&mut self) -> Parsed<bool> {
        let mut look = self.lookahead();
        if self.peek(&mut look)?.kind != Kind::Lb {
            return Ok(false);
        }
        while self.peek(&mut look)?.kind == Kind::Lb {
            if !self.skip_bracket(&mut look)? {
                return Ok(false);
            }
        }
        Ok(self.peek(&mut look)?.kind == Kind::Minus)
    }
    fn global_entry(&mut self) -> Parsed<GlobalModEntry> {
        if self.take(Kind::Lb)? {
            let modification = self.modification()?;
            self.expect(Kind::Rb, "']'")?;
            self.work.consume(modification.alternatives.len())?;
            if modification
                .alternatives
                .iter()
                .any(|(_, label)| label.is_some())
            {
                return self.fail(
                    ErrorCode::UnexpectedCharacter,
                    "Labels are not allowed on global modifications",
                );
            }
            self.expect(Kind::At, "'@' for global modification locations")?;
            let mut locations = Vec::new();
            if self.check(Kind::Id)? {
                let location = self.location()?;
                self.work.push(&mut locations, location)?;
            }
            while self.check(Kind::Comma)? {
                let mut look = self.lexer.clone();
                let after = self.peek(&mut look)?;
                if matches!(after.kind, Kind::Lb | Kind::Number) {
                    break;
                }
                self.advance()?;
                if self.check(Kind::Id)? {
                    let location = self.location()?;
                    self.work.push(&mut locations, location)?;
                }
            }
            Ok(GlobalModEntry::GlobalModification(GlobalModification {
                modification,
                locations,
            }))
        } else {
            let position = self.current()?.start;
            let mut isotope = Vec::new();
            if self.check(Kind::Number)? {
                let token = self.advance()?;
                self.append_token(&mut isotope, token)?;
            }
            if self.check(Kind::Id)? {
                let token = self.advance()?;
                self.append_token(&mut isotope, token)?;
            }
            if isotope.is_empty() {
                return self.fail(ErrorCode::InvalidFormula, "Expected isotope specification");
            }
            Ok(GlobalModEntry::IsotopeReplacement(IsotopeReplacement {
                isotope: self.finish(isotope, position)?,
            }))
        }
    }
    fn location(&mut self) -> Parsed<String> {
        let first = self.advance()?;
        let mut text = Vec::new();
        self.append_token(&mut text, first)?;
        if self.take(Kind::Minus)? {
            let token = self.expect(Kind::Id, "term identifier")?;
            self.work.append(&mut text, b"-")?;
            self.append_token(&mut text, token)?;
        }
        if self.take(Kind::Colon)? {
            let token = self.expect(Kind::Id, "location suffix")?;
            self.work.append(&mut text, b":")?;
            self.append_token(&mut text, token)?;
        }
        self.finish(text, first.start)
    }
    fn unlocalised(&mut self) -> Parsed<Vec<UnlocalisedMod>> {
        let mut output = Vec::new();
        while self.check(Kind::Lb)? {
            let mut look = self.lookahead();
            let mut question = false;
            loop {
                if self.peek(&mut look)?.kind != Kind::Lb {
                    break;
                }
                self.skip_bracket(&mut look)?;
                if self.peek(&mut look)?.kind == Kind::Caret {
                    self.next(&mut look)?;
                    if self.peek(&mut look)?.kind == Kind::Number {
                        self.next(&mut look)?;
                    }
                }
                if self.peek(&mut look)?.kind == Kind::Question {
                    question = true;
                    break;
                }
                if self.peek(&mut look)?.kind != Kind::Lb {
                    break;
                }
            }
            if !question {
                break;
            }
            while self.take(Kind::Lb)? {
                let mut entry = UnlocalisedMod::default();
                let modification = self.modification()?;
                self.work.push(&mut entry.modifications, modification)?;
                self.expect(Kind::Rb, "']'")?;
                if self.take(Kind::Caret)? {
                    let token = self.expect(Kind::Number, "occurrence count")?;
                    entry.occurrence = Some(self.integer(
                        token,
                        ErrorCode::InvalidMassValue,
                        "Invalid occurrence count",
                    )?);
                }
                self.work.push(&mut output, entry)?;
                if self.check(Kind::Question)? {
                    break;
                }
            }
            self.expect(Kind::Question, "'?' for unlocalised modification")?;
        }
        Ok(output)
    }
    fn sequence(&mut self) -> Parsed<Vec<SequenceSection>> {
        let mut sequence = Vec::new();
        loop {
            let token = self.current()?;
            if token.kind == Kind::Lp {
                let mut look = self.lookahead();
                self.next(&mut look)?;
                let ambiguous = self.peek(&mut look)?.kind == Kind::Question;
                self.advance()?;
                if ambiguous {
                    self.expect(Kind::Question, "'?' for ambiguous region")?;
                }
                let mut elements = Vec::new();
                while !self.check(Kind::Rp)? && !self.check(Kind::End)? {
                    let token = self.current()?;
                    if token.kind == Kind::Id {
                        self.advance()?;
                        for (i, byte) in self.text(token).iter().copied().enumerate() {
                            if !byte.is_ascii_alphabetic() {
                                return Err(self.syntax(
                                    ErrorCode::InvalidAminoAcid,
                                    if ambiguous {
                                        token.start
                                    } else {
                                        token.start + i
                                    },
                                    if ambiguous {
                                        "Invalid amino acid in ambiguous region"
                                    } else {
                                        "Invalid amino acid in range"
                                    },
                                ));
                            }
                            let modifications =
                                if i + 1 == token.end - token.start && self.check(Kind::Lb)? {
                                    self.modifications()?
                                } else {
                                    Vec::new()
                                };
                            self.work.push(
                                &mut elements,
                                SequenceElement {
                                    amino_acid: char::from(byte),
                                    modifications,
                                },
                            )?;
                        }
                    } else if !ambiguous && token.kind == Kind::Lb && !elements.is_empty() {
                        let mods = self.modifications()?;
                        for modification in mods {
                            self.work.push(
                                &mut elements.last_mut().unwrap().modifications,
                                modification,
                            )?;
                        }
                    } else {
                        break;
                    }
                }
                self.expect(Kind::Rp, "')'")?;
                let section = if ambiguous {
                    SequenceSection::AmbiguousRegion(AmbiguousRegion { elements })
                } else {
                    if elements.is_empty() {
                        return self
                            .fail(ErrorCode::UnexpectedCharacter, "Empty range is not allowed");
                    }
                    let modifications = if self.check(Kind::Lb)? {
                        self.modifications()?
                    } else {
                        Vec::new()
                    };
                    SequenceSection::ModifiedRange(ModifiedRange {
                        elements,
                        modifications,
                    })
                };
                self.work.push(&mut sequence, section)?;
            } else if token.kind == Kind::Id {
                self.advance()?;
                for (i, byte) in self.text(token).iter().copied().enumerate() {
                    if !byte.is_ascii_alphabetic() {
                        return Err(self.syntax(
                            ErrorCode::InvalidAminoAcid,
                            token.start,
                            "Invalid amino acid character",
                        ));
                    }
                    let modifications =
                        if i + 1 == token.end - token.start && self.check(Kind::Lb)? {
                            self.modifications()?
                        } else {
                            Vec::new()
                        };
                    self.work.push(
                        &mut sequence,
                        SequenceSection::Element(SequenceElement {
                            amino_acid: char::from(byte),
                            modifications,
                        }),
                    )?;
                }
            } else {
                break;
            }
        }
        Ok(sequence)
    }
    fn modifications(&mut self) -> Parsed<Vec<Modification>> {
        let mut output = Vec::new();
        while self.take(Kind::Lb)? {
            let modification = self.modification()?;
            self.work.push(&mut output, modification)?;
            self.expect(Kind::Rb, "']'")?;
        }
        Ok(output)
    }
    fn modification(&mut self) -> Parsed<Modification> {
        let mut modification = Modification::default();
        loop {
            let alternative = if self.check(Kind::Hash)? {
                (
                    ModificationTag::InfoTag(InfoTag::default()),
                    Some(self.label()?),
                )
            } else {
                let tag = self.tag()?;
                let label = if self.check(Kind::Hash)? {
                    Some(self.label()?)
                } else {
                    None
                };
                (tag, label)
            };
            self.work
                .push(&mut modification.alternatives, alternative)?;
            if !self.take(Kind::Pipe)? {
                break;
            }
        }
        Ok(modification)
    }
    fn tag(&mut self) -> Parsed<ModificationTag> {
        let token = self.current()?;
        if matches!(token.kind, Kind::Plus | Kind::Minus | Kind::Number) {
            return Ok(ModificationTag::MassDelta(self.mass()?));
        }
        if token.kind != Kind::Id {
            return self.fail(ErrorCode::UnexpectedCharacter, "Expected modification");
        }
        match self.text(token) {
            b"Formula" => {
                self.advance()?;
                self.expect(Kind::Colon, "':' after Formula")?;
                Ok(ModificationTag::FormulaTag(self.formula()?))
            }
            b"Glycan" => {
                self.advance()?;
                self.expect(Kind::Colon, "':' after Glycan")?;
                Ok(ModificationTag::GlycanComposition(self.glycan()?))
            }
            b"INFO" | b"info" | b"Info" => {
                self.advance()?;
                self.expect(Kind::Colon, "':' after INFO")?;
                Ok(ModificationTag::InfoTag(self.info()?))
            }
            b"Position" | b"position" | b"POSITION" => {
                self.advance()?;
                self.expect(Kind::Colon, "':' after Position")?;
                Ok(ModificationTag::PositionConstraint(self.position()?))
            }
            b"Cation" => {
                self.advance()?;
                self.expect(Kind::Colon, "':' after Cation")?;
                let position = self.current()?.start;
                let mut name = Vec::new();
                self.work.append(&mut name, b"Cation:")?;
                if self.check(Kind::Id)? {
                    let token = self.advance()?;
                    self.append_token(&mut name, token)?;
                }
                if self.take(Kind::Lb)? {
                    self.work.append(&mut name, b"[")?;
                    while !self.check(Kind::Rb)? && !self.check(Kind::End)? {
                        let token = self.advance()?;
                        self.append_token(&mut name, token)?;
                    }
                    self.expect(Kind::Rb, "']' for cation charge")?;
                    self.work.append(&mut name, b"]")?;
                }
                Ok(ModificationTag::NamedMod(NamedMod {
                    name: self.finish(name, position)?,
                    cv_hint: None,
                }))
            }
            b"Obs" => {
                self.advance()?;
                self.expect(Kind::Colon, "':' after Obs")?;
                let mut mass = self.mass()?;
                mass.source = MassDeltaSource::Obs;
                Ok(ModificationTag::MassDelta(mass))
            }
            b"UNIMOD" | b"MOD" | b"RESID" | b"XLMOD" | b"GNO" | b"unimod" | b"mod" | b"resid"
            | b"xlmod" | b"gno" | b"xlink" => Ok(ModificationTag::CvAccession(self.cv()?)),
            id => {
                let hint = match id {
                    b"U" => Some((CvDatabase::Unimod, MassDeltaSource::U)),
                    b"M" => Some((CvDatabase::Mod, MassDeltaSource::M)),
                    b"R" => Some((CvDatabase::Resid, MassDeltaSource::R)),
                    b"X" => Some((CvDatabase::Xlmod, MassDeltaSource::X)),
                    b"G" => Some((CvDatabase::Gno, MassDeltaSource::G)),
                    _ => None,
                };
                if let Some((database, source)) = hint {
                    let mut look = self.lexer.clone();
                    if self.peek(&mut look)?.kind == Kind::Colon {
                        self.advance()?;
                        self.advance()?;
                        if matches!(
                            self.current()?.kind,
                            Kind::Plus | Kind::Minus | Kind::Number
                        ) {
                            let mut mass = self.mass()?;
                            mass.source = source;
                            return Ok(ModificationTag::MassDelta(mass));
                        }
                        let mut named = self.named()?;
                        named.cv_hint = Some(database);
                        return Ok(ModificationTag::NamedMod(named));
                    }
                }
                Ok(ModificationTag::NamedMod(self.named()?))
            }
        }
    }
    fn named(&mut self) -> Parsed<NamedMod> {
        let position = self.current()?.start;
        let mut name = Vec::new();
        let (mut parens, mut brackets) = (0, 0);
        loop {
            let token = self.current()?;
            match token.kind {
                Kind::Id | Kind::Number => (),
                Kind::Minus => {
                    let mut look = self.lexer.clone();
                    let next = self.peek(&mut look)?;
                    if next.kind == Kind::Ra {
                        self.work.append(&mut name, b"->")?;
                        self.advance()?;
                        self.advance()?;
                        continue;
                    }
                    if !matches!(next.kind, Kind::Id | Kind::Number | Kind::Lp | Kind::Lb)
                        && parens == 0
                        && brackets == 0
                    {
                        break;
                    }
                }
                Kind::Lp => {
                    if name.is_empty() && parens == 0 && brackets == 0 {
                        break;
                    }
                    Self::increment_depth(&mut parens)?;
                }
                Kind::Rp => {
                    if parens == 0 {
                        break;
                    }
                    parens -= 1;
                }
                Kind::Lb => {
                    if name.is_empty() && parens == 0 && brackets == 0 {
                        break;
                    }
                    Self::increment_depth(&mut brackets)?;
                }
                Kind::Rb => {
                    if brackets == 0 {
                        break;
                    }
                    brackets -= 1;
                }
                _ => break,
            }
            self.append_token(&mut name, token)?;
            self.advance()?;
        }
        if name.is_empty() {
            return self.fail(ErrorCode::UnexpectedCharacter, "Expected modification name");
        }
        Ok(NamedMod {
            cv_hint: None,
            name: self.finish(name, position)?,
        })
    }
    fn cv(&mut self) -> Parsed<CvAccession> {
        let token = self.expect(Kind::Id, "CV prefix")?;
        let database = match self.text(token) {
            b"UNIMOD" | b"unimod" => CvDatabase::Unimod,
            b"MOD" | b"mod" => CvDatabase::Mod,
            b"RESID" | b"resid" => CvDatabase::Resid,
            b"XLMOD" | b"xlmod" | b"xlink" => CvDatabase::Xlmod,
            b"GNO" | b"gno" => CvDatabase::Gno,
            _ => {
                return Err(self.syntax(
                    ErrorCode::InvalidCvPrefix,
                    token.start,
                    "Invalid CV prefix",
                ));
            }
        };
        self.expect(Kind::Colon, "':' after CV prefix")?;
        let accession = if matches!(database, CvDatabase::Unimod | CvDatabase::Mod) {
            let token = self.expect(Kind::Number, "accession number")?;
            self.string(token)?
        } else {
            let position = self.current()?.start;
            let mut bytes = Vec::new();
            let mut depth = 0;
            loop {
                let token = self.current()?;
                match token.kind {
                    Kind::Id | Kind::Number => (),
                    Kind::Lb => Self::increment_depth(&mut depth)?,
                    Kind::Rb => {
                        if depth == 0 {
                            break;
                        }
                        depth -= 1;
                    }
                    _ => break,
                }
                self.append_token(&mut bytes, token)?;
                self.advance()?;
            }
            if bytes.is_empty() {
                return self.fail(ErrorCode::InvalidCvAccession, "Expected accession");
            }
            self.finish(bytes, position)?
        };
        Ok(CvAccession {
            database,
            accession,
        })
    }
    fn mass(&mut self) -> Parsed<MassDelta> {
        let position = self.current()?.start;
        let mut text = Vec::new();
        if self.take(Kind::Plus)? {
            self.work.append(&mut text, b"+")?;
        } else if self.take(Kind::Minus)? {
            self.work.append(&mut text, b"-")?;
        }
        if self.check(Kind::Number)? {
            let token = self.advance()?;
            self.append_token(&mut text, token)?;
        } else {
            return self.fail(ErrorCode::InvalidMassValue, "Expected mass value");
        }
        let error_position = self.current()?.start;
        let mass = self.float(&text, error_position, "Invalid mass value format")?;
        Ok(MassDelta {
            source: MassDeltaSource::None,
            mass,
            original_text: self.finish(text, position)?,
        })
    }
    fn formula(&mut self) -> Parsed<FormulaTag> {
        let position = self.current()?.start;
        let mut text = Vec::new();
        let mut brackets = 0;
        let mut charge = None;
        loop {
            let token = self.current()?;
            match token.kind {
                Kind::Id | Kind::Number | Kind::Lp | Kind::Rp => (),
                Kind::Lb => Self::increment_depth(&mut brackets)?,
                Kind::Rb => {
                    if brackets == 0 {
                        break;
                    }
                    brackets -= 1;
                }
                Kind::Minus => {
                    let mut look = self.lexer.clone();
                    if brackets == 0 && self.peek(&mut look)?.kind != Kind::Number {
                        break;
                    }
                }
                Kind::Colon => {
                    let mut look = self.lexer.clone();
                    let next = self.peek(&mut look)?;
                    if next.kind == Kind::Id && self.is_text(next, b"z") {
                        self.advance()?;
                        self.advance()?;
                        charge = Some(self.signed_integer("Invalid charge value")?);
                    }
                    break;
                }
                _ => break,
            }
            self.append_token(&mut text, token)?;
            self.advance()?;
        }
        if text.is_empty() {
            return self.fail(ErrorCode::InvalidFormula, "Empty formula");
        }
        Ok(FormulaTag {
            formula_string: self.finish(text, position)?,
            charge,
        })
    }
    fn glycan(&mut self) -> Parsed<GlycanComposition> {
        let mut glycan = GlycanComposition::default();
        while self.check(Kind::Id)? {
            let token = self.advance()?;
            let name = self.string(token)?;
            let count = if self.check(Kind::Number)? {
                let token = self.advance()?;
                self.integer(
                    token,
                    ErrorCode::InvalidMassValue,
                    "Invalid monosaccharide count",
                )?
            } else {
                1
            };
            self.work
                .push(&mut glycan.components, (GlycanComponent::Name(name), count))?;
        }
        if glycan.components.is_empty() {
            return self.fail(ErrorCode::UnknownMonosaccharide, "Empty glycan composition");
        }
        Ok(glycan)
    }
    fn info(&mut self) -> Parsed<InfoTag> {
        let position = self.current()?.start;
        let mut text = Vec::new();
        loop {
            let token = self.current()?;
            if matches!(
                token.kind,
                Kind::End | Kind::Rb | Kind::Pipe | Kind::Hash | Kind::Comma
            ) {
                break;
            }
            self.append_token(&mut text, token)?;
            self.advance()?;
        }
        Ok(InfoTag {
            text: self.finish(text, position)?,
        })
    }
    fn position(&mut self) -> Parsed<PositionConstraint> {
        let mut position = PositionConstraint::default();
        loop {
            let token = self.current()?;
            if matches!(
                token.kind,
                Kind::End | Kind::Rb | Kind::Pipe | Kind::Hash | Kind::Rc
            ) {
                break;
            }
            if token.kind == Kind::Id && matches!(self.text(token), b"N" | b"C") {
                let mut look = self.lookahead();
                self.next(&mut look)?;
                if self.peek(&mut look)?.kind == Kind::Minus {
                    self.next(&mut look)?;
                    let next = self.peek(&mut look)?;
                    if next.kind == Kind::Id && self.is_text(next, b"term") {
                        if self.is_text(token, b"N") {
                            position.n_term = true;
                        } else {
                            position.c_term = true;
                        }
                        self.advance()?;
                        self.advance()?;
                        self.advance()?;
                        self.take(Kind::Comma)?;
                        continue;
                    }
                }
            }
            if self.take(Kind::Comma)? {
                continue;
            }
            if token.kind == Kind::Id {
                for byte in self.text(token) {
                    if !byte.is_ascii_alphabetic() {
                        return Err(self.syntax(
                            ErrorCode::InvalidAminoAcid,
                            token.start,
                            "Invalid amino acid in position constraint",
                        ));
                    }
                    self.work.push(&mut position.residues, char::from(*byte))?;
                }
                self.advance()?;
            } else {
                return self.fail(
                    ErrorCode::UnexpectedCharacter,
                    "Expected amino acid residues or terminal position after Position:",
                );
            }
        }
        if position.residues.is_empty() && !position.n_term && !position.c_term {
            return self.fail(
                ErrorCode::UnexpectedCharacter,
                "Position constraint requires at least one residue or terminal position",
            );
        }
        Ok(position)
    }
    fn label(&mut self) -> Parsed<Label> {
        self.expect(Kind::Hash, "'#'")?;
        let position = self.current()?.start;
        let mut name = Vec::new();
        if self.check(Kind::Id)? {
            let token = self.advance()?;
            self.append_token(&mut name, token)?;
            if self.check(Kind::Number)? {
                let token = self.advance()?;
                self.append_token(&mut name, token)?;
            }
        } else if self.check(Kind::Number)? {
            let token = self.advance()?;
            self.append_token(&mut name, token)?;
        } else {
            return self.fail(ErrorCode::UnexpectedCharacter, "Expected label identifier");
        }
        let identifier = self.finish(name, position)?;
        let label_type = if identifier == "BRANCH" {
            LabelType::Branch
        } else if identifier.starts_with("XL") {
            LabelType::Crosslink
        } else {
            LabelType::Ambiguous
        };
        let score = if self.take(Kind::Lp)? {
            let token = self.expect(Kind::Number, "score value")?;
            let bytes = self.text(token);
            let value = self.float(bytes, token.start, "Invalid score value")?;
            self.expect(Kind::Rp, "')'")?;
            Some(value)
        } else {
            None
        };
        Ok(Label {
            label_type,
            identifier,
            score,
        })
    }
    fn charge(&mut self) -> Parsed<ChargeState> {
        if self.take(Kind::Lb)? {
            let mut adducts = Vec::new();
            let adduct = self.adduct()?;
            self.work.push(&mut adducts, adduct)?;
            while self.take(Kind::Comma)? {
                let adduct = self.adduct()?;
                self.work.push(&mut adducts, adduct)?;
            }
            self.expect(Kind::Rb, "']'")?;
            Ok(ChargeState::Adducts(adducts))
        } else {
            let sign = if self.take(Kind::Plus)? {
                1
            } else if self.take(Kind::Minus)? {
                -1
            } else {
                1
            };
            if !self.check(Kind::Number)? {
                return self.fail(ErrorCode::InvalidCharge, "Expected charge value");
            }
            let token = self.current()?;
            let value = self.integer(token, ErrorCode::InvalidCharge, "Invalid charge value")?;
            let charge = value.checked_mul(sign).ok_or_else(|| {
                self.syntax(
                    ErrorCode::InvalidCharge,
                    token.start,
                    "Invalid charge value",
                )
            })?;
            self.advance()?;
            Ok(ChargeState::Simple(charge))
        }
    }
    fn adduct(&mut self) -> Parsed<AdductIon> {
        let position = self.current()?.start;
        let mut text = Vec::new();
        while !self.check(Kind::End)? {
            if self.check(Kind::Colon)? {
                let mut look = self.lookahead();
                self.next(&mut look)?;
                let next = self.peek(&mut look)?;
                if next.kind == Kind::Id && self.is_text(next, b"z") {
                    break;
                }
            }
            let token = self.current()?;
            if matches!(token.kind, Kind::Rb | Kind::Comma) {
                return self.fail(
                    ErrorCode::UnexpectedCharacter,
                    "Unexpected token in adduct formula, expected ':z'",
                );
            }
            self.append_token(&mut text, token)?;
            self.advance()?;
        }
        if text.is_empty() {
            return self.fail(ErrorCode::UnexpectedCharacter, "Expected adduct formula");
        }
        self.expect(Kind::Colon, "':'")?;
        self.expect(Kind::Id, "'z'")?;
        let charge = self.signed_integer("Invalid adduct charge value")?;
        let occurrence = if self.take(Kind::Caret)? {
            let token = self.expect(Kind::Number, "occurrence count")?;
            Some(self.integer(
                token,
                ErrorCode::InvalidMassValue,
                "Invalid adduct occurrence count",
            )?)
        } else {
            None
        };
        Ok(AdductIon {
            formula: self.finish(text, position)?,
            charge,
            occurrence,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ion_chains_share_the_node_budget_until_final_publication() {
        let mut parser = Parser::new("A+B").unwrap();
        parser.work.nodes = 3;
        assert!(matches!(parser.ion(), Err(ParseFailure::Resource(_))));
        assert_eq!(parser.work.nodes, 0);
        let mut parser = Parser::new("A+B").unwrap();
        parser.work.nodes = 4;
        assert_eq!(parser.ion().unwrap().chains.len(), 2);
        assert_eq!(parser.work.nodes, 0);
    }

    #[test]
    fn lookahead_rescans_consume_the_same_remaining_work() {
        let mut parser = Parser::new("[a_long_annotation]-A").unwrap();
        assert!(parser.n_terminal_pattern().unwrap());
        let remaining = parser.work.remaining;
        assert!(parser.n_terminal_pattern().unwrap());
        assert!(parser.work.remaining < remaining);
        parser.work.remaining = 1;
        assert!(matches!(
            parser.n_terminal_pattern(),
            Err(ParseFailure::Resource(_))
        ));
    }

    #[test]
    fn output_growth_and_error_diagnostics_are_precharged() {
        let mut work = Work::new();
        work.bytes = 0;
        let mut text = Vec::new();
        assert!(work.append(&mut text, b"annotation").is_err());
        assert_eq!(text.capacity(), 0);
        let mut parser = Parser::new("[").unwrap();
        parser.work.bytes = 0;
        assert!(matches!(
            parser.syntax(ErrorCode::UnexpectedCharacter, 0, "Expected residue"),
            ParseFailure::Resource(_)
        ));
    }

    #[test]
    fn late_work_exhaustion_does_not_publish_a_partial_ion() {
        let text = "(>one)A[Oxidation]/2+(>two)B[+3.25]/3";
        let mut complete = Parser::new(text).unwrap();
        let before = complete.work.remaining;
        assert_eq!(complete.ion().unwrap().chains.len(), 2);
        let consumed = before - complete.work.remaining;
        let mut bounded = Parser::new(text).unwrap();
        bounded.work.remaining = consumed - 1;
        assert!(matches!(bounded.ion(), Err(ParseFailure::Resource(_))));
    }
}
