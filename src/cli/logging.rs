// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The `-log` file and the debug levels: the native form of `TOPPBase`'s
//! `enableLogging_`, `writeLogInfo_`, `writeLogWarn_`, `writeLogError_` and
//! both `writeDebug_` overloads (`TOPPBase.cpp:2473-2538`), and of the
//! destructor that removes an empty log file (`TOPPBase.cpp:131-139`).
//!
//! As in the Release build, debug text reaches the log file only: the source
//! also streams it to `OPENMS_LOG_DEBUG`, which a Release build compiles out,
//! so `-debug` without `-log` prints nothing (oracle `debug1_run` …
//! `debug10_run` in `../oracle/toppbase-completion`). Info, warning and error
//! text goes to the console stream the caller passes *and* to the log file.
//!
//! A log line is `<YYYY-MM-DD hh:mm:ss> <ini location>: <text>`, in local
//! time; a parameter dump is framed by [`LOG_SEPARATOR`] lines. The file is
//! opened for appending on the first line written, and only once a
//! destination is known; a destination that cannot be opened writes nothing
//! and is not reported, as the source's unchecked `std::ofstream` does.

use crate::data_structures::DateTime;
use crate::param::Param;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::sync::Mutex;

/// The separator line around a debug parameter dump (source `LOG_SEPARATOR`).
pub const LOG_SEPARATOR: &str =
    " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ";

#[derive(Debug, Default)]
struct State {
    /// The source's `getIniLocation_()`.
    location: String,
    /// The source's `debug_level_`; -1 until the command line is parsed.
    debug_level: i64,
    /// The `log` parameter of the current parameter set, if any.
    destination: Option<String>,
    /// The open log file.
    file: Option<File>,
    /// The destination the file was opened for.
    opened: Option<String>,
}

/// The log file of one tool run (source members `log_` and `debug_level_`).
#[derive(Debug)]
pub struct ToolLog {
    state: Mutex<State>,
}

impl Default for ToolLog {
    fn default() -> Self {
        Self::new()
    }
}

impl ToolLog {
    /// A log with no destination and the source's initial debug level, -1.
    pub fn new() -> Self {
        Self {
            state: Mutex::new(State {
                debug_level: -1,
                ..State::default()
            }),
        }
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, State> {
        match self.state.lock() {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        }
    }

    /// Set the INI location that prefixes every line.
    pub fn set_location(&self, location: &str) {
        self.lock().location = location.to_owned();
    }

    /// Set the debug level, as the source assigns `debug_level_`.
    pub fn set_debug_level(&self, level: i64) {
        self.lock().debug_level = level;
    }

    /// The current debug level.
    pub fn debug_level(&self) -> i64 {
        self.lock().debug_level
    }

    /// Set the `log` value of the parameters in force: the source's
    /// `enableLogging_` reads `param_`, which is the command line at first and
    /// the resolved parameters after the INI merge. `None` or an empty value
    /// means no destination.
    pub fn set_destination(&self, destination: Option<&str>) {
        self.lock().destination = destination.filter(|d| !d.is_empty()).map(str::to_owned);
    }

    /// Source `enableLogging_`: open the destination once. Returns the
    /// `Writing to '<file>'` notice the source prints on standard output when
    /// the debug level is at least 1, for the caller to write there.
    fn enable(state: &mut State) -> Option<String> {
        if state.opened.is_some() {
            return None;
        }
        let destination = state.destination.clone()?;
        state.file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&destination)
            .ok();
        state.opened = Some(destination.clone());
        if state.debug_level >= 1 {
            let notice = format!("Writing to '{destination}'");
            let line = format!("{} {}: {notice}\n", DateTime::now().get(), state.location);
            if let Some(file) = state.file.as_mut() {
                let _ = file.write_all(line.as_bytes());
            }
            return Some(notice);
        }
        None
    }

    /// Append one line to the file, opening it first if needed.
    fn append(state: &mut State, text: &str) -> Option<String> {
        let notice = Self::enable(state);
        let line = format!("{} {}: {text}\n", DateTime::now().get(), state.location);
        if let Some(file) = state.file.as_mut() {
            let _ = file.write_all(line.as_bytes());
        }
        notice
    }

    /// The file part of `writeLogInfo_`, `writeLogWarn_` and `writeLogError_`:
    /// the caller writes `text` to its console stream. Returns the
    /// `Writing to` notice when this opened the file at debug level 1 or more.
    pub fn line(&self, text: &str) -> Option<String> {
        Self::append(&mut self.lock(), text)
    }

    /// Source `writeDebug_(text, min_level)`: a log line when the debug level
    /// is at least `min_level`. Returns the `Writing to` notice as
    /// [`line`](Self::line) does.
    pub fn debug(&self, text: &str, min_level: u32) -> Option<String> {
        let mut state = self.lock();
        if state.debug_level >= i64::from(min_level) {
            Self::append(&mut state, text)
        } else {
            None
        }
    }

    /// Source `writeDebug_(text, param, min_level)`: the text and the
    /// parameters between two separator lines. Returns the `Writing to` notice
    /// as [`line`](Self::line) does.
    pub fn debug_param(&self, text: &str, param: &Param, min_level: u32) -> Option<String> {
        let mut state = self.lock();
        if state.debug_level < i64::from(min_level) {
            return None;
        }
        let notice = Self::enable(&mut state);
        let block = format!(
            "{LOG_SEPARATOR}\n{} {} {text}\n{}{LOG_SEPARATOR}\n",
            DateTime::now().get(),
            state.location,
            param.to_text().unwrap_or_default()
        );
        if let Some(file) = state.file.as_mut() {
            let _ = file.write_all(block.as_bytes());
        }
        notice
    }

    /// The source's destructor: close the file, then remove the file the
    /// `log` parameter names when it is an empty regular file. As in the
    /// source, that holds whether or not this run wrote to it.
    pub fn finish(&self) {
        let mut state = self.lock();
        if let Some(mut file) = state.file.take() {
            let _ = file.flush();
        }
        if let Some(path) = state.destination.clone() {
            if std::fs::metadata(&path).is_ok_and(|m| m.is_file() && m.len() == 0) {
                let _ = std::fs::remove_file(&path);
            }
        }
    }
}
