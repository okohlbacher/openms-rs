// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Naive source DateTime state, parsing and Gregorian arithmetic.
use crate::{Error, Result};
use chrono::{Datelike, Timelike};
use std::{
    fmt,
    hash::{Hash, Hasher},
    str::FromStr,
};

/// Maximum supplied date/time text, including bytes after a C-string NUL.
pub const MAX_DATETIME_INPUT_BYTES: usize = 1_048_576;
/// All source format strings, in source dispatch order.
pub const DATETIME_FORMATS: [&str; 7] = [
    "yyyy-MM-ddThh:mm:ss",
    "yyyy-MM-ddThh:mm:ss.zzz",
    "yyyy-MM-dd hh:mm:ss",
    "yyyy-MM-dd+hh:mm",
    "yyyy-MM-ddThh:mm:ssZ",
    "yyyy-MM-dd",
    "hh:mm:ss",
];

/// Fixed-size calendar fields with the source's independent validity flag.
/// A time-only value is valid despite a zero date. Arithmetic preserves validity.
/// Equality includes it; `source_less` intentionally does not.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DateTime {
    year: i32,
    month: i32,
    day: i32,
    hour: i32,
    minute: i32,
    second: i32,
    millisecond: i32,
    valid: bool,
}
impl DateTime {
    pub fn parse(value: &str) -> Result<Self> {
        let mut result = Self::default();
        result.set(value)?;
        Ok(result)
    }
    /// Automatic source parsing clears the value on syntax/calendar failure.
    /// Input-size rejection happens before clearing or scanning.
    pub fn set(&mut self, value: &str) -> Result<()> {
        bounded(value)?;
        self.clear();
        let mut d = Self::default();
        let mut parsed = false;
        if value.contains('.') && !value.contains('T') {
            let s = scan(value, "%d.%d.%d %d:%d:%d");
            d.assign(&s, &[2, 1, 0, 3, 4, 5]);
            parsed = s.count == 6;
        } else if value.contains('/') {
            let s = scan(value, "%d/%d/%d %d:%d:%d");
            d.assign(&s, &[1, 2, 0, 3, 4, 5]);
            parsed = s.count == 6;
        } else if value.contains('-') {
            if value.contains('T') {
                let selected = value.split_once('+').map_or(value, |(before, _)| before);
                let fraction = selected.contains('.');
                let s = scan(
                    selected,
                    if fraction {
                        "%d-%d-%dT%d:%d:%d.%d"
                    } else {
                        "%d-%d-%dT%d:%d:%d"
                    },
                );
                d.assign(&s, &[0, 1, 2, 3, 4, 5, 6]);
                parsed = s.count == if fraction { 7 } else { 6 };
                if parsed && fraction {
                    d.millisecond =
                        normalize_fraction(selected, d.millisecond).ok_or_else(parse_error)?;
                }
            } else if value.contains('Z') {
                let s = scan(value, "%d-%d-%dZ");
                d.assign(&s, &[0, 1, 2]);
                parsed = s.count == 3 && value.find('Z') == Some(value.len() - 1);
            } else if value.contains('+') {
                let s = scan(value, "%d-%d-%d+%d:%d");
                d.assign(&s, &[0, 1, 2, 3, 4]);
                parsed = s.count == 5;
            } else {
                let s = scan(value, "%d-%d-%d %d:%d:%d");
                d.assign(&s, &[0, 1, 2, 3, 4, 5]);
                parsed = s.count == 6;
            }
        }
        if !parsed {
            let s = scan(value, "%3s %3s %d %d:%d:%d %d");
            if s.count == 7 {
                if let Some(month) = month_abbrev(s.words[1]) {
                    d.year = s.ints[4];
                    d.month = month;
                    d.day = s.ints[0];
                    d.hour = s.ints[1];
                    d.minute = s.ints[2];
                    d.second = s.ints[3];
                    parsed = true;
                }
            } else {
                let s = scan(value, "%3s %3s %d %d");
                if s.count == 4 {
                    if let Some(month) = month_abbrev(s.words[1]) {
                        d.year = s.ints[1];
                        d.month = month;
                        d.day = s.ints[0];
                        parsed = true;
                    }
                }
            }
        }
        if !parsed || !valid_date(d.year, d.month, d.day) || !valid_time(d.hour, d.minute, d.second)
        {
            return Err(parse_error());
        }
        d.valid = true;
        *self = d;
        Ok(())
    }
    /// Source numeric setter, month/day/year followed by hour/minute/second.
    pub fn set_components(
        &mut self,
        month: u32,
        day: u32,
        year: u32,
        hour: u32,
        minute: u32,
        second: u32,
    ) -> Result<()> {
        let d = Self {
            year: signed(year)?,
            month: signed(month)?,
            day: signed(day)?,
            hour: signed(hour)?,
            minute: signed(minute)?,
            second: signed(second)?,
            millisecond: 0,
            valid: true,
        };
        if !valid_date(d.year, d.month, d.day) || !valid_time(d.hour, d.minute, d.second) {
            return Err(parse_error());
        }
        *self = d;
        Ok(())
    }
    pub fn set_date(&mut self, value: &str) -> Result<()> {
        bounded(value)?;
        let (s, order) = if value.contains('-') {
            (scan(value, "%d-%d-%d"), [0, 1, 2])
        } else if value.contains('.') {
            (scan(value, "%d.%d.%d"), [2, 1, 0])
        } else if value.contains('/') {
            (scan(value, "%d/%d/%d"), [1, 2, 0])
        } else {
            return Err(parse_error());
        };
        let mut d = *self;
        d.assign(&s, &order);
        if s.count != 3 || !valid_date(d.year, d.month, d.day) {
            return Err(parse_error());
        }
        d.valid = true;
        *self = d;
        Ok(())
    }
    pub fn set_date_components(&mut self, month: u32, day: u32, year: u32) -> Result<()> {
        let (m, d, y) = (signed(month)?, signed(day)?, signed(year)?);
        if !valid_date(y, m, d) {
            return Err(parse_error());
        }
        self.year = y;
        self.month = m;
        self.day = d;
        self.valid = true;
        Ok(())
    }
    pub fn set_time(&mut self, value: &str) -> Result<()> {
        bounded(value)?;
        let s = scan(value, "%d:%d:%d");
        if s.count != 3 || !valid_time(s.ints[0], s.ints[1], s.ints[2]) {
            return Err(parse_error());
        }
        self.hour = s.ints[0];
        self.minute = s.ints[1];
        self.second = s.ints[2];
        self.valid = true;
        Ok(())
    }
    pub fn set_time_components(&mut self, hour: u32, minute: u32, second: u32) -> Result<()> {
        let (h, m, s) = (signed(hour)?, signed(minute)?, signed(second)?);
        if !valid_time(h, m, s) {
            return Err(parse_error());
        }
        self.hour = h;
        self.minute = m;
        self.second = s;
        self.valid = true;
        Ok(())
    }
    /// Source-order signed components. C++ casts negative years to UInt here.
    pub const fn components(self) -> (i32, i32, i32, i32, i32, i32) {
        (
            self.month,
            self.day,
            self.year,
            self.hour,
            self.minute,
            self.second,
        )
    }
    pub const fn date_components(self) -> (i32, i32, i32) {
        (self.month, self.day, self.year)
    }
    pub const fn time_components(self) -> (i32, i32, i32) {
        (self.hour, self.minute, self.second)
    }
    pub const fn millisecond(self) -> i32 {
        self.millisecond
    }
    pub const fn is_valid(self) -> bool {
        self.valid
    }
    pub const fn is_null(self) -> bool {
        !self.valid
    }
    pub fn clear(&mut self) {
        *self = Self::default();
    }
    pub fn source_less(&self, other: &Self) -> bool {
        self.sort_key() < other.sort_key()
    }
    fn sort_key(self) -> [i32; 7] {
        [
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
            self.millisecond,
        ]
    }
    /// Formatted text is empty for any invalid value, even for an unknown format.
    pub fn format(self, format: &str) -> Result<String> {
        if !self.valid {
            return Ok(String::new());
        }
        let index = DATETIME_FORMATS
            .iter()
            .position(|f| *f == format)
            .ok_or_else(|| Error::InvalidValue("unknown DateTime format".into()))?;
        Ok(self.render(index))
    }
    pub fn iso_string(self) -> String {
        if self.valid {
            self.render(0)
        } else {
            String::new()
        }
    }
    pub fn get(self) -> String {
        if self.valid {
            self.render(2)
        } else {
            "0000-00-00 00:00:00".into()
        }
    }
    pub fn date_string(self) -> String {
        if self.valid {
            self.render(5)
        } else {
            "0000-00-00".into()
        }
    }
    pub fn time_string(self) -> String {
        if self.valid {
            self.render(6)
        } else {
            "00:00:00".into()
        }
    }
    fn render(self, index: usize) -> String {
        let Self {
            year: y,
            month: m,
            day: d,
            hour: h,
            minute: min,
            second: s,
            millisecond: ms,
            ..
        } = self;
        match index {
            0 => format!("{y:04}-{m:02}-{d:02}T{h:02}:{min:02}:{s:02}"),
            1 => format!("{y:04}-{m:02}-{d:02}T{h:02}:{min:02}:{s:02}.{ms:03}"),
            2 => format!("{y:04}-{m:02}-{d:02} {h:02}:{min:02}:{s:02}"),
            3 => format!("{y:04}-{m:02}-{d:02}+{h:02}:{min:02}"),
            4 => format!("{y:04}-{m:02}-{d:02}T{h:02}:{min:02}:{s:02}Z"),
            5 => format!("{y:04}-{m:02}-{d:02}"),
            _ => format!("{h:02}:{min:02}:{s:02}"),
        }
    }
    /// Unsupported formats or malformed/calendar-invalid input return an invalid
    /// value. Only the input bound returns an operation error.
    pub fn from_format(value: &str, format: &str) -> Result<Self> {
        bounded(value)?;
        let Some(index) = DATETIME_FORMATS.iter().position(|f| *f == format) else {
            return Ok(Self::default());
        };
        let (pattern, expected) = match index {
            0 => ("%d-%d-%dT%d:%d:%d", 6),
            1 => ("%d-%d-%dT%d:%d:%d.%d", 7),
            2 => ("%d-%d-%d %d:%d:%d", 6),
            3 => ("%d-%d-%d+%d:%d", 5),
            4 => ("%d-%d-%dT%d:%d:%dZ", 6),
            5 => ("%d-%d-%d", 3),
            _ => ("%d:%d:%d", 3),
        };
        let s = scan(value, pattern);
        let mut d = Self::default();
        if s.count != expected {
            return Ok(d);
        }
        if index == 6 {
            d.assign(&s, &[3, 4, 5]);
        } else {
            d.assign(&s, &[0, 1, 2, 3, 4, 5, 6]);
        }
        if index == 1 {
            let Some(ms) = normalize_fraction(value, d.millisecond) else {
                return Ok(Self::default());
            };
            d.millisecond = ms;
        }
        if (index != 6 && !valid_date(d.year, d.month, d.day))
            || !valid_time(d.hour, d.minute, d.second)
        {
            return Ok(Self::default());
        }
        d.valid = true;
        Ok(d)
    }
    /// Naive Gregorian UTC normalization; does not apply local timezone or DST.
    /// Year overflow rejects atomically. Milliseconds and validity are preserved.
    pub fn add_seconds(&mut self, seconds: i32) -> Result<&mut Self> {
        let months = i64::from(self.year) * 12 + i64::from(self.month) - 1;
        let year = months.div_euclid(12);
        let month = months.rem_euclid(12) + 1;
        let days = day_number(year, month, i64::from(self.day));
        let total = days * 86_400
            + i64::from(self.hour) * 3_600
            + i64::from(self.minute) * 60
            + i64::from(self.second)
            + i64::from(seconds);
        let (year, month, day) = civil_date(total.div_euclid(86_400));
        let year = i32::try_from(year)
            .map_err(|_| Error::InvalidValue("DateTime year exceeds i32".into()))?;
        let time = total.rem_euclid(86_400);
        self.year = year;
        self.month = month as i32;
        self.day = day as i32;
        self.hour = (time / 3600) as i32;
        self.minute = (time / 60 % 60) as i32;
        self.second = (time % 60) as i32;
        Ok(self)
    }
    pub fn now() -> Self {
        Self::clock_fields(chrono::Local::now())
    }
    pub fn now_utc() -> Self {
        Self::clock_fields(chrono::Utc::now())
    }
    fn clock_fields<T: chrono::TimeZone>(value: chrono::DateTime<T>) -> Self {
        Self {
            year: value.year(),
            month: value.month() as i32,
            day: value.day() as i32,
            hour: value.hour() as i32,
            minute: value.minute() as i32,
            second: value.second() as i32,
            millisecond: 0,
            valid: true,
        }
    }
    // Scan integer slots map onto year/month/day/hour/minute/second/millisecond.
    fn assign(&mut self, s: &Scanned<'_>, order: &[usize]) {
        for (slot, &value) in order.iter().zip(&s.ints) {
            match slot {
                0 => self.year = value,
                1 => self.month = value,
                2 => self.day = value,
                3 => self.hour = value,
                4 => self.minute = value,
                5 => self.second = value,
                _ => self.millisecond = value,
            }
        }
    }
}
impl FromStr for DateTime {
    type Err = Error;
    fn from_str(value: &str) -> Result<Self> {
        Self::parse(value)
    }
}
impl fmt::Display for DateTime {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.get())
    }
}
impl Hash for DateTime {
    fn hash<H: Hasher>(&self, state: &mut H) {
        if self.valid {
            self.render(1).hash(state)
        } else {
            "".hash(state)
        }
    }
}
fn bounded(value: &str) -> Result<()> {
    if value.len() > MAX_DATETIME_INPUT_BYTES {
        Err(Error::InvalidValue(
            "DateTime input exceeds byte limit".into(),
        ))
    } else {
        Ok(())
    }
}
fn parse_error() -> Error {
    Error::InvalidValue("invalid DateTime input or calendar fields".into())
}
fn signed(value: u32) -> Result<i32> {
    i32::try_from(value).map_err(|_| parse_error())
}
fn valid_time(h: i32, m: i32, s: i32) -> bool {
    (0..24).contains(&h) && (0..60).contains(&m) && (0..60).contains(&s)
}
fn valid_date(y: i32, m: i32, d: i32) -> bool {
    if y < 1 || !(1..=12).contains(&m) || d < 1 {
        return false;
    }
    let mut days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][(m - 1) as usize];
    if m == 2 && (y % 4 == 0 && y % 100 != 0 || y % 400 == 0) {
        days = 29;
    }
    d <= days
}
fn month_abbrev(value: &[u8]) -> Option<i32> {
    [
        b"Jan", b"Feb", b"Mar", b"Apr", b"May", b"Jun", b"Jul", b"Aug", b"Sep", b"Oct", b"Nov",
        b"Dec",
    ]
    .iter()
    .position(|m| m.as_slice() == value)
    .map(|i| i as i32 + 1)
}
fn normalize_fraction(value: &str, mut ms: i32) -> Option<i32> {
    if let Some(dot) = value.find('.') {
        let digits = value.as_bytes()[dot + 1..]
            .iter()
            .take_while(|c| c.is_ascii_digit())
            .count();
        for _ in digits..3 {
            ms = ms.checked_mul(10)?;
        }
        ms %= 1000;
    }
    Some(ms)
}
fn whitespace(c: u8) -> bool {
    c == b' ' || (b'\t'..=b'\r').contains(&c)
}
struct Scanned<'a> {
    ints: [i32; 7],
    words: [&'a [u8]; 2],
    count: usize,
}
// Only source-used %d and %3s conversions. Literal mismatch after the last
// assignment leaves its successful count, exactly like scanf. No heap buffers.
fn scan<'a>(value: &'a str, pattern: &str) -> Scanned<'a> {
    let bytes = value.as_bytes();
    let bytes = &bytes[..bytes.iter().position(|&c| c == 0).unwrap_or(bytes.len())];
    let p = pattern.as_bytes();
    let (mut i, mut j, mut ni, mut nw) = (0, 0, 0, 0);
    let mut out = Scanned {
        ints: [0; 7],
        words: [&[]; 2],
        count: 0,
    };
    while j < p.len() {
        if whitespace(p[j]) {
            while i < bytes.len() && whitespace(bytes[i]) {
                i += 1;
            }
            j += 1;
            continue;
        }
        if p[j] != b'%' {
            if bytes.get(i) != Some(&p[j]) {
                break;
            }
            i += 1;
            j += 1;
            continue;
        }
        j += 1;
        while i < bytes.len() && whitespace(bytes[i]) {
            i += 1;
        }
        if p[j] == b'3' {
            let start = i;
            while i < bytes.len() && i - start < 3 && !whitespace(bytes[i]) {
                i += 1;
            }
            if i == start {
                break;
            }
            out.words[nw] = &bytes[start..i];
            nw += 1;
            j += 2;
        } else {
            let negative = bytes.get(i) == Some(&b'-');
            if negative || bytes.get(i) == Some(&b'+') {
                i += 1;
            }
            let start = i;
            let mut n = 0i64;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                n = n * 10 + i64::from(bytes[i] - b'0');
                if n > i64::from(i32::MAX) + i64::from(negative) {
                    return out;
                }
                i += 1;
            }
            if i == start {
                break;
            }
            out.ints[ni] = if negative { -n as i32 } else { n as i32 };
            ni += 1;
            j += 1;
        }
        out.count += 1;
    }
    out
}
// March-based proleptic Gregorian eras; all intermediates fit i64 throughout
// the i32-year state domain, even after adding any signed i32 second interval.
fn day_number(mut y: i64, m: i64, d: i64) -> i64 {
    y -= i64::from(m <= 2);
    let era = y.div_euclid(400);
    let yo = y - era * 400;
    let mp = m + if m > 2 { -3 } else { 9 };
    era * 146097 + yo * 365 + yo / 4 - yo / 100 + (153 * mp + 2) / 5 + d - 1
}
fn civil_date(z: i64) -> (i64, i64, i64) {
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yo = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yo + era * 400;
    let doy = doe - (365 * yo + yo / 4 - yo / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = mp + if mp < 10 { 3 } else { -9 };
    (y + i64::from(m <= 2), m, d)
}
