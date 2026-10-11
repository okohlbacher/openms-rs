// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The syntax check of `nlohmann::json` 3.11.2's `json::parse`, with its
//! exception texts, as `ParamJSONFile::load` meets them
//! (`ParamJSONFile.cpp:63` and `162-165` at core `bc9cc12`).
//!
//! `serde_json` decides the values; this module decides whether the document
//! is JSON at all, where it ends, and what the source's exception says when it
//! is not. It follows the vendored library at the same pin
//! (`src/openms/extern/nlohmann_json/include/nlohmann/detail/input/lexer.hpp`
//! and `parser.hpp`) token for token: the scanner's position bookkeeping, the
//! bytes it remembers for `last read: '…'`, and the parser's states. The
//! values themselves are not built. Its answer:
//!
//! * the byte offset at which the document ends: its end, or its first NUL
//!   byte after a complete value, which the source's scanner takes for the end
//!   of input (`lexer.hpp`, `case '\0'`), so bytes after it are never read;
//!   and the offsets of the number tokens that are exactly `-0`, which the
//!   source reads as the integer 0 (`strtoll`) where `serde_json` reads the
//!   float -0.0 ([`Checked`]);
//! * or the source exception's `what()`: a `parse_error.101` with its line,
//!   column and message, or the `out_of_range.406` a number beyond `double`
//!   raises when the parser reaches it (`parser.hpp`, `value_float`).
//!
//! The source's comment support (`ignore_comments`) is off in
//! `json::parse(std::ifstream)` and is not ported.

/// The scanner's token types (`lexer_base::token_type`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Token {
    Uninitialized,
    LiteralTrue,
    LiteralFalse,
    LiteralNull,
    ValueString,
    ValueUnsigned,
    ValueInteger,
    ValueFloat,
    BeginArray,
    BeginObject,
    EndArray,
    EndObject,
    NameSeparator,
    ValueSeparator,
    ParseError,
    EndOfInput,
    LiteralOrValue,
}

/// `lexer_base::token_type_name`.
fn token_name(token: Token) -> &'static str {
    match token {
        Token::Uninitialized => "<uninitialized>",
        Token::LiteralTrue => "true literal",
        Token::LiteralFalse => "false literal",
        Token::LiteralNull => "null literal",
        Token::ValueString => "string literal",
        Token::ValueUnsigned | Token::ValueInteger | Token::ValueFloat => "number literal",
        Token::BeginArray => "'['",
        Token::BeginObject => "'{'",
        Token::EndArray => "']'",
        Token::EndObject => "'}'",
        Token::NameSeparator => "':'",
        Token::ValueSeparator => "','",
        Token::ParseError => "<parse error>",
        Token::EndOfInput => "end of input",
        Token::LiteralOrValue => "'[', '{', or a literal",
    }
}

/// The scanner's messages for an unescaped control character in a string
/// (`lexer.hpp`, `scan_string`, cases `0x00` to `0x1F`).
const CONTROL_MESSAGES: [&str; 32] = [
    "invalid string: control character U+0000 (NUL) must be escaped to \\u0000",
    "invalid string: control character U+0001 (SOH) must be escaped to \\u0001",
    "invalid string: control character U+0002 (STX) must be escaped to \\u0002",
    "invalid string: control character U+0003 (ETX) must be escaped to \\u0003",
    "invalid string: control character U+0004 (EOT) must be escaped to \\u0004",
    "invalid string: control character U+0005 (ENQ) must be escaped to \\u0005",
    "invalid string: control character U+0006 (ACK) must be escaped to \\u0006",
    "invalid string: control character U+0007 (BEL) must be escaped to \\u0007",
    "invalid string: control character U+0008 (BS) must be escaped to \\u0008 or \\b",
    "invalid string: control character U+0009 (HT) must be escaped to \\u0009 or \\t",
    "invalid string: control character U+000A (LF) must be escaped to \\u000A or \\n",
    "invalid string: control character U+000B (VT) must be escaped to \\u000B",
    "invalid string: control character U+000C (FF) must be escaped to \\u000C or \\f",
    "invalid string: control character U+000D (CR) must be escaped to \\u000D or \\r",
    "invalid string: control character U+000E (SO) must be escaped to \\u000E",
    "invalid string: control character U+000F (SI) must be escaped to \\u000F",
    "invalid string: control character U+0010 (DLE) must be escaped to \\u0010",
    "invalid string: control character U+0011 (DC1) must be escaped to \\u0011",
    "invalid string: control character U+0012 (DC2) must be escaped to \\u0012",
    "invalid string: control character U+0013 (DC3) must be escaped to \\u0013",
    "invalid string: control character U+0014 (DC4) must be escaped to \\u0014",
    "invalid string: control character U+0015 (NAK) must be escaped to \\u0015",
    "invalid string: control character U+0016 (SYN) must be escaped to \\u0016",
    "invalid string: control character U+0017 (ETB) must be escaped to \\u0017",
    "invalid string: control character U+0018 (CAN) must be escaped to \\u0018",
    "invalid string: control character U+0019 (EM) must be escaped to \\u0019",
    "invalid string: control character U+001A (SUB) must be escaped to \\u001A",
    "invalid string: control character U+001B (ESC) must be escaped to \\u001B",
    "invalid string: control character U+001C (FS) must be escaped to \\u001C",
    "invalid string: control character U+001D (GS) must be escaped to \\u001D",
    "invalid string: control character U+001E (RS) must be escaped to \\u001E",
    "invalid string: control character U+001F (US) must be escaped to \\u001F",
];

const ILL_FORMED_UTF8: &str = "invalid string: ill-formed UTF-8 byte";
const HEX_DIGITS: &str = "invalid string: '\\u' must be followed by 4 hex digits";
const HIGH_SURROGATE: &str =
    "invalid string: surrogate U+D800..U+DBFF must be followed by U+DC00..U+DFFF";

/// The scanner (`detail::lexer`), without the values it would build.
struct Lexer<'a> {
    input: &'a [u8],
    /// The next byte `get` reads.
    next: usize,
    /// The last byte read, `None` at the end of the input (`eof()`).
    current: Option<u8>,
    next_unget: bool,
    chars_read_total: usize,
    chars_read_current_line: usize,
    lines_read: usize,
    /// The raw bytes read since the last token reset (`token_string`).
    token_string: Vec<u8>,
    /// A number's characters (`token_buffer`); strings are not kept.
    number: Vec<u8>,
    error_message: &'static str,
    /// A floating-point token's value (`value_float`).
    value_float: f64,
    /// The offsets of the number tokens `-0` scanned so far.
    negative_zeros: Vec<usize>,
}

impl<'a> Lexer<'a> {
    fn new(input: &'a [u8]) -> Self {
        Self {
            input,
            next: 0,
            current: None,
            next_unget: false,
            chars_read_total: 0,
            chars_read_current_line: 0,
            lines_read: 0,
            token_string: Vec::new(),
            number: Vec::new(),
            error_message: "",
            value_float: 0.0,
            negative_zeros: Vec::new(),
        }
    }

    /// `get`: the position counts even the end of input.
    fn get(&mut self) -> Option<u8> {
        self.chars_read_total += 1;
        self.chars_read_current_line += 1;
        if self.next_unget {
            self.next_unget = false;
        } else {
            self.current = self.input.get(self.next).copied();
            if self.current.is_some() {
                self.next += 1;
            }
        }
        if let Some(byte) = self.current {
            self.token_string.push(byte);
        }
        if self.current == Some(b'\n') {
            self.lines_read += 1;
            self.chars_read_current_line = 0;
        }
        self.current
    }

    /// `unget`: read the current byte again on the next `get`.
    fn unget(&mut self) {
        self.next_unget = true;
        self.chars_read_total = self.chars_read_total.saturating_sub(1);
        if self.chars_read_current_line == 0 {
            if self.lines_read > 0 {
                self.lines_read -= 1;
            }
        } else {
            self.chars_read_current_line -= 1;
        }
        if self.current.is_some() {
            self.token_string.pop();
        }
    }

    /// `reset`: a new token starts with the current byte.
    fn reset(&mut self) {
        self.number.clear();
        self.token_string.clear();
        if let Some(byte) = self.current {
            self.token_string.push(byte);
        }
    }

    fn error(&mut self, message: &'static str) -> Token {
        self.error_message = message;
        Token::ParseError
    }

    /// `get_token_string`: control characters as `<U+00XX>`, every other byte
    /// as it is.
    fn token_string(&self) -> Vec<u8> {
        let mut text = Vec::with_capacity(self.token_string.len());
        for &byte in &self.token_string {
            if byte <= 0x1F {
                text.extend_from_slice(format!("<U+{byte:04X}>").as_bytes());
            } else {
                text.push(byte);
            }
        }
        text
    }

    /// `position_string`: `line` is one-based, `column` the bytes read on it.
    fn position(&self) -> String {
        format!(
            " at line {}, column {}",
            self.lines_read + 1,
            self.chars_read_current_line
        )
    }

    fn skip_bom(&mut self) -> bool {
        if self.get() == Some(0xEF) {
            return self.get() == Some(0xBB) && self.get() == Some(0xBF);
        }
        self.unget();
        true
    }

    fn skip_whitespace(&mut self) {
        loop {
            self.get();
            if !matches!(self.current, Some(b' ' | b'\t' | b'\n' | b'\r')) {
                break;
            }
        }
    }

    fn scan(&mut self) -> Token {
        if self.chars_read_total == 0 && !self.skip_bom() {
            return self.error("invalid BOM; must be 0xEF 0xBB 0xBF if given");
        }
        self.skip_whitespace();
        match self.current {
            Some(b'[') => Token::BeginArray,
            Some(b']') => Token::EndArray,
            Some(b'{') => Token::BeginObject,
            Some(b'}') => Token::EndObject,
            Some(b':') => Token::NameSeparator,
            Some(b',') => Token::ValueSeparator,
            Some(b't') => self.scan_literal(b"true", Token::LiteralTrue),
            Some(b'f') => self.scan_literal(b"false", Token::LiteralFalse),
            Some(b'n') => self.scan_literal(b"null", Token::LiteralNull),
            Some(b'"') => self.scan_string(),
            Some(b'-' | b'0'..=b'9') => self.scan_number(),
            Some(0) | None => Token::EndOfInput,
            Some(_) => self.error("invalid literal"),
        }
    }

    fn scan_literal(&mut self, text: &[u8], token: Token) -> Token {
        for &expected in &text[1..] {
            if self.get() != Some(expected) {
                return self.error("invalid literal");
            }
        }
        token
    }

    /// `get_codepoint`: four hexadecimal digits, or `None`.
    fn codepoint(&mut self) -> Option<u32> {
        let mut codepoint = 0;
        for factor in [12, 8, 4, 0] {
            let digit = match self.get() {
                Some(byte @ b'0'..=b'9') => u32::from(byte - b'0'),
                Some(byte @ b'A'..=b'F') => u32::from(byte - b'A') + 10,
                Some(byte @ b'a'..=b'f') => u32::from(byte - b'a') + 10,
                _ => return None,
            };
            codepoint += digit << factor;
        }
        Some(codepoint)
    }

    /// `next_byte_in_range`: each following byte within its inclusive range.
    fn next_bytes_in(&mut self, ranges: &[(u8, u8)]) -> bool {
        for &(low, high) in ranges {
            match self.get() {
                Some(byte) if (low..=high).contains(&byte) => {}
                _ => {
                    self.error_message = ILL_FORMED_UTF8;
                    return false;
                }
            }
        }
        true
    }

    fn scan_string(&mut self) -> Token {
        self.reset();
        loop {
            let Some(byte) = self.get() else {
                return self.error("invalid string: missing closing quote");
            };
            let well_formed = match byte {
                b'"' => return Token::ValueString,
                b'\\' => match self.get() {
                    Some(b'"' | b'\\' | b'/' | b'b' | b'f' | b'n' | b'r' | b't') => true,
                    Some(b'u') => {
                        let Some(first) = self.codepoint() else {
                            return self.error(HEX_DIGITS);
                        };
                        if (0xD800..=0xDBFF).contains(&first) {
                            if self.get() == Some(b'\\') && self.get() == Some(b'u') {
                                let Some(second) = self.codepoint() else {
                                    return self.error(HEX_DIGITS);
                                };
                                if !(0xDC00..=0xDFFF).contains(&second) {
                                    return self.error(HIGH_SURROGATE);
                                }
                            } else {
                                return self.error(HIGH_SURROGATE);
                            }
                        } else if (0xDC00..=0xDFFF).contains(&first) {
                            return self.error(
                                "invalid string: surrogate U+DC00..U+DFFF must follow U+D800..U+DBFF",
                            );
                        }
                        true
                    }
                    _ => {
                        return self.error("invalid string: forbidden character after backslash");
                    }
                },
                0x00..=0x1F => return self.error(CONTROL_MESSAGES[usize::from(byte)]),
                0x20..=0x7F => true,
                0xC2..=0xDF => self.next_bytes_in(&[(0x80, 0xBF)]),
                0xE0 => self.next_bytes_in(&[(0xA0, 0xBF), (0x80, 0xBF)]),
                0xE1..=0xEC | 0xEE..=0xEF => self.next_bytes_in(&[(0x80, 0xBF), (0x80, 0xBF)]),
                0xED => self.next_bytes_in(&[(0x80, 0x9F), (0x80, 0xBF)]),
                0xF0 => self.next_bytes_in(&[(0x90, 0xBF), (0x80, 0xBF), (0x80, 0xBF)]),
                0xF1..=0xF3 => self.next_bytes_in(&[(0x80, 0xBF), (0x80, 0xBF), (0x80, 0xBF)]),
                0xF4 => self.next_bytes_in(&[(0x80, 0x8F), (0x80, 0xBF), (0x80, 0xBF)]),
                _ => return self.error(ILL_FORMED_UTF8),
            };
            if !well_formed {
                return Token::ParseError;
            }
        }
    }

    fn digits(&mut self) {
        while let Some(byte @ b'0'..=b'9') = self.get() {
            self.number.push(byte);
        }
    }

    /// `scan_number`: the state machine of RFC 8259's grammar, then the value
    /// of a floating-point token, which the parser checks.
    fn scan_number(&mut self) -> Token {
        // `scan` saw a byte, so it is the one before `next`.
        let start = self.next.saturating_sub(1);
        self.reset();
        let mut kind = Token::ValueUnsigned;
        let first = self.current.unwrap_or(b'0');
        self.number.push(first);
        if first == b'-' {
            kind = Token::ValueInteger;
            match self.get() {
                Some(byte @ b'0'..=b'9') => self.number.push(byte),
                _ => return self.error("invalid number; expected digit after '-'"),
            }
        }
        // `zero` reads no further digits; `any1` does.
        let leading_zero = self.number.last() == Some(&b'0');
        if leading_zero {
            self.get();
        } else {
            self.digits();
        }
        if self.current == Some(b'.') {
            self.number.push(b'.');
            kind = Token::ValueFloat;
            match self.get() {
                Some(byte @ b'0'..=b'9') => self.number.push(byte),
                _ => return self.error("invalid number; expected digit after '.'"),
            }
            self.digits();
        }
        if let Some(e @ (b'e' | b'E')) = self.current {
            self.number.push(e);
            kind = Token::ValueFloat;
            match self.get() {
                Some(sign @ (b'+' | b'-')) => {
                    self.number.push(sign);
                    match self.get() {
                        Some(byte @ b'0'..=b'9') => self.number.push(byte),
                        _ => {
                            return self
                                .error("invalid number; expected digit after exponent sign");
                        }
                    }
                }
                Some(byte @ b'0'..=b'9') => self.number.push(byte),
                _ => {
                    return self
                        .error("invalid number; expected '+', '-', or digit after exponent");
                }
            }
            self.digits();
        }
        self.unget();
        let text = String::from_utf8_lossy(&self.number);
        let integral = match kind {
            Token::ValueUnsigned => text.parse::<u64>().is_ok(),
            Token::ValueInteger => text.parse::<i64>().is_ok(),
            _ => false,
        };
        if integral {
            if self.number == b"-0" {
                self.negative_zeros.push(start);
            }
            return kind;
        }
        // `strtod` in the "C" locale: the token is always a complete decimal
        // number here, which Rust reads to the same correctly rounded value,
        // infinity on overflow included.
        self.value_float = text.parse::<f64>().unwrap_or(f64::INFINITY);
        Token::ValueFloat
    }
}

/// The parser (`detail::parser`) over [`Lexer`].
struct Parser<'a> {
    lexer: Lexer<'a>,
    last: Token,
}

impl Parser<'_> {
    fn next_token(&mut self) -> Token {
        self.last = self.lexer.scan();
        self.last
    }

    /// `exception_message`, framed as `parse_error::create(101, …)`.
    fn syntax_error(&self, expected: Token, context: &str) -> Vec<u8> {
        let mut what = format!(
            "[json.exception.parse_error.101] parse error{}: syntax error ",
            self.lexer.position()
        )
        .into_bytes();
        if !context.is_empty() {
            what.extend_from_slice(format!("while parsing {context} ").as_bytes());
        }
        what.extend_from_slice(b"- ");
        if self.last == Token::ParseError {
            what.extend_from_slice(self.lexer.error_message.as_bytes());
            what.extend_from_slice(b"; last read: '");
            what.extend_from_slice(&self.lexer.token_string());
            what.push(b'\'');
        } else {
            what.extend_from_slice(format!("unexpected {}", token_name(self.last)).as_bytes());
        }
        if expected != Token::Uninitialized {
            what.extend_from_slice(format!("; expected {}", token_name(expected)).as_bytes());
        }
        what
    }

    /// `sax_parse_internal`: one value, arrays and objects as a stack.
    fn value(&mut self) -> Result<(), Vec<u8>> {
        // true: an array, false: an object.
        let mut states: Vec<bool> = Vec::new();
        let mut skip_to_state_evaluation = false;
        loop {
            if !skip_to_state_evaluation {
                match self.last {
                    Token::BeginObject => {
                        if self.next_token() != Token::EndObject {
                            if self.last != Token::ValueString {
                                return Err(self.syntax_error(Token::ValueString, "object key"));
                            }
                            if self.next_token() != Token::NameSeparator {
                                return Err(
                                    self.syntax_error(Token::NameSeparator, "object separator")
                                );
                            }
                            states.push(false);
                            self.next_token();
                            continue;
                        }
                    }
                    Token::BeginArray => {
                        if self.next_token() != Token::EndArray {
                            states.push(true);
                            continue;
                        }
                    }
                    Token::ValueFloat => {
                        if !self.lexer.value_float.is_finite() {
                            let mut what =
                                b"[json.exception.out_of_range.406] number overflow parsing '"
                                    .to_vec();
                            what.extend_from_slice(&self.lexer.token_string());
                            what.push(b'\'');
                            return Err(what);
                        }
                    }
                    Token::LiteralFalse
                    | Token::LiteralNull
                    | Token::LiteralTrue
                    | Token::ValueInteger
                    | Token::ValueString
                    | Token::ValueUnsigned => {}
                    Token::ParseError => {
                        return Err(self.syntax_error(Token::Uninitialized, "value"));
                    }
                    Token::Uninitialized
                    | Token::EndArray
                    | Token::EndObject
                    | Token::NameSeparator
                    | Token::ValueSeparator
                    | Token::EndOfInput
                    | Token::LiteralOrValue => {
                        return Err(self.syntax_error(Token::LiteralOrValue, "value"));
                    }
                }
            } else {
                skip_to_state_evaluation = false;
            }
            let Some(&array) = states.last() else {
                return Ok(());
            };
            if array {
                if self.next_token() == Token::ValueSeparator {
                    self.next_token();
                    continue;
                }
                if self.last == Token::EndArray {
                    states.pop();
                    skip_to_state_evaluation = true;
                    continue;
                }
                return Err(self.syntax_error(Token::EndArray, "array"));
            }
            if self.next_token() == Token::ValueSeparator {
                if self.next_token() != Token::ValueString {
                    return Err(self.syntax_error(Token::ValueString, "object key"));
                }
                if self.next_token() != Token::NameSeparator {
                    return Err(self.syntax_error(Token::NameSeparator, "object separator"));
                }
                self.next_token();
                continue;
            }
            if self.last == Token::EndObject {
                states.pop();
                skip_to_state_evaluation = true;
                continue;
            }
            return Err(self.syntax_error(Token::EndObject, "object"));
        }
    }
}

/// What [`check`] learned about a document the source accepts.
#[derive(Debug)]
pub(super) struct Checked {
    /// The offset of the byte that ended the input: the length, or a NUL byte
    /// after the value.
    pub(super) end: usize,
    /// The offsets, all below `end`, of the number tokens that are exactly
    /// `-0`: the source's `strtoll` makes each the integer 0.
    pub(super) negative_zeros: Vec<usize>,
}

/// `json::parse(input)` in strict mode, as far as it succeeds or fails: what
/// [`Checked`] holds, or the `what()` of the exception the source throws. The
/// text is the source's bytes; a byte of `last read: '…'` may be any byte the
/// document holds.
pub(super) fn check(input: &[u8]) -> Result<Checked, Vec<u8>> {
    let mut parser = Parser {
        lexer: Lexer::new(input),
        last: Token::Uninitialized,
    };
    parser.next_token();
    parser.value()?;
    if parser.next_token() != Token::EndOfInput {
        return Err(parser.syntax_error(Token::EndOfInput, "value"));
    }
    let end = if parser.lexer.current.is_some() {
        parser.lexer.next - 1
    } else {
        input.len()
    };
    Ok(Checked {
        end,
        negative_zeros: parser.lexer.negative_zeros,
    })
}
