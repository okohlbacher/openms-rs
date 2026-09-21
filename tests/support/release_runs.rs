// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The retained runs of `../oracle/topp-exception-exits`: how the C++ Release
//! build (core `bc9cc12`, cli `c19e494`, topp `174b576`, on ibminode06) ends a
//! tool run that returns an exit code and one that throws, and what reaches the
//! `-log` file (`tests/data/topp_exception_exits/`, provenance in
//! `tests/data/topp_exception_exits_provenance.json`).
//!
//! Each case directory holds the Release run's `argv.txt`, `env.txt`,
//! `exit_code.txt`, `stdout.txt`, `stderr.txt`, `tree.txt` and, where the run
//! left one, its `log.txt`. [`ReleaseRun::replay`] runs the same command line
//! through [`run_with_registry`] in a fresh directory, and
//! [`ReleaseRun::release`] gives the Release streams with the Release run's
//! paths mapped to the replay's. They are compared as recorded, except for
//! three things of the environment: the Linux `stty` probe line before a usage
//! text is dropped (a tool driven in process writes to explicit streams,
//! which are not probed), the figures of the closing `<tool> took …` line are
//! checked for shape rather than value (`support/took_line.rs`), and the log
//! file's time stamps are masked.

#![allow(dead_code)]

use openms::cli::{ExitCode, Tool, ToolHandler, ToolRegistrySources, run_with_registry};
use openms::system::file::TempDir;
use std::fs;
use std::path::{Path, PathBuf};

/// Where the Release run kept its case directories.
pub const ORACLE_RESULTS: &str = "/scratch/kohlbach/topp-exception-exits/run2/results";
/// Where the Release run read its inputs.
pub const ORACLE_INPUTS: &str = "/scratch/kohlbach/topp-exception-exits/inputs";
/// The Release installation prefix.
pub const RELEASE_PREFIX: &str =
    "/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576";
/// The terminal-width probe line the Linux Release build prints before usage
/// text when standard input is not a terminal.
pub const STTY_LINE: &str = "stty: 'standard input': Inappropriate ioctl for device\n";

/// The repository copy of each Release input, by the name the oracle gave it;
/// the oracle's `inputs.sha256` (retained beside the cases) records their
/// hashes, which equal these files'. An empty path is an input the test
/// derives itself.
const INPUTS: [(&str, &str); 12] = [
    (
        "FileInfo_1_input.dta",
        "file_info/inputs/FileInfo_1_input.dta",
    ),
    (
        "FileInfo_8_input.notype",
        "mzml_mobility/FileInfo_8_input.notype",
    ),
    ("x_1.txt", "fuzzy_string_comparator/fuzzydiff/x_1.txt"),
    ("x_1001.txt", "fuzzy_string_comparator/fuzzydiff/x_1001.txt"),
    (
        "PeakPickerHiRes_6_noforce.ini",
        "topp_peak_picker_hi_res/PeakPickerHiRes_6_noforce.ini",
    ),
    (
        "PeakPickerHiRes_6_input.mzML",
        "peak_picking/PeakPickerHiRes_6_input.mzML",
    ),
    (
        "PeakPickerHiRes_input.mzML",
        "peak_picking/PeakPickerHiRes_input.mzML",
    ),
    ("empty.mzML", "topp_peak_picker_hi_res/empty.mzML"),
    (
        "baseline_filter_tool_input.mzML",
        "baseline_filter_tool_input.mzML",
    ),
    (
        "SimpleSearchEngine_1.mzML",
        "topp_feature_finder_centroided/SimpleSearchEngine_1.mzML",
    ),
    (
        "FeatureFinderCentroided_1_input.mzML",
        "mzml_mobility/FeatureFinderCentroided_1_input.mzML",
    ),
    ("profile/FeatureFinderCentroided_1_input.mzML", ""),
];

fn repository(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn text(path: impl AsRef<Path>) -> String {
    path.as_ref().to_string_lossy().into_owned()
}

/// The repository file standing for the Release input `name`, or `None` for an
/// input the test derives itself.
pub fn input(name: &str) -> Option<PathBuf> {
    INPUTS
        .iter()
        .find(|(oracle, _)| *oracle == name)
        .and_then(|(_, relative)| {
            (!relative.is_empty()).then(|| repository(&format!("tests/data/{relative}")))
        })
}

/// One retained Release run.
pub struct ReleaseRun {
    pub name: &'static str,
    dir: PathBuf,
}

/// What a replay produced.
pub struct Replay {
    pub code: ExitCode,
    /// Standard output without the closing `<tool> took …` line.
    pub out: String,
    /// The closing line, when there was one.
    pub took: Option<String>,
    pub err: String,
    /// The `-log` file, time stamps masked; `None` when no file was left.
    pub log: Option<Vec<String>>,
    /// The case's working directory.
    pub cwd: PathBuf,
    _work: TempDir,
}

/// The Release side of a case, with its paths mapped to a replay's.
pub struct Released {
    pub exit: i32,
    /// Standard output without the closing `<tool> took …` line.
    pub out: String,
    /// The closing line, when there was one.
    pub took: Option<String>,
    /// Standard error without the `stty` probe line.
    pub err: String,
    /// The `-log` file, time stamps masked; `None` when the run left none.
    pub log: Option<Vec<String>>,
}

impl ReleaseRun {
    pub fn new(name: &'static str) -> Self {
        let dir = repository("tests/data/topp_exception_exits").join(name);
        assert!(dir.is_dir(), "no retained Release run {name}");
        Self { name, dir }
    }

    fn file(&self, name: &str) -> Option<String> {
        fs::read_to_string(self.dir.join(name)).ok()
    }

    /// The oracle input paths and what they are here; `derived` names the
    /// inputs the test derives itself, which take the place of the table's.
    fn input_map(&self, derived: &[(&str, &Path)]) -> Vec<(String, String)> {
        let mut map: Vec<(String, String)> = INPUTS
            .iter()
            .filter_map(|(name, _)| {
                let local = derived
                    .iter()
                    .find(|(oracle, _)| oracle == name)
                    .map(|(_, path)| path.to_path_buf())
                    .or_else(|| input(name))?;
                Some((format!("{ORACLE_INPUTS}/{name}"), text(local)))
            })
            .collect();
        // The longest name first, so `profile/X` is not taken for `X`.
        map.sort_by_key(|(oracle, _)| std::cmp::Reverse(oracle.len()));
        map
    }

    /// A Release file with the case's working directory, the inputs and the
    /// installation's `bin/` mapped, and the `stty` probe line dropped.
    fn mapped(&self, file: &str, cwd: &Path, derived: &[(&str, &Path)]) -> Option<String> {
        let mut content = self.file(file)?;
        content = content.replace(&format!("{ORACLE_RESULTS}/{}/cwd", self.name), &text(cwd));
        for (oracle, local) in self.input_map(derived) {
            content = content.replace(&oracle, &local);
        }
        content = content.replace(&format!("{RELEASE_PREFIX}/bin/"), "");
        Some(content.replace(STTY_LINE, ""))
    }

    /// The Release run of `tool`, mapped to `replay`'s directory.
    pub fn release(&self, tool: &str, replay: &Replay, derived: &[(&str, &Path)]) -> Released {
        let stdout = self.mapped("stdout.txt", &replay.cwd, derived).unwrap();
        let (out, took) = super::took_line::split_took_line(tool, &stdout);
        Released {
            exit: self.file("exit_code.txt").unwrap().trim().parse().unwrap(),
            out,
            took,
            err: self.mapped("stderr.txt", &replay.cwd, derived).unwrap(),
            log: self
                .mapped("log.txt", &replay.cwd, derived)
                .map(|content| masked_log(&content)),
        }
    }

    /// Run the Release command line through `T` in a fresh working directory
    /// that `setup` prepares as the oracle prepared the case's, with the
    /// registry of an uninstalled executable (the built-in manifest, no
    /// `OPENMS_TOOL_PREFIX_PATH`, no internal tools) and the case's own user
    /// directory, as the oracle ran each case with its own `HOME`.
    pub fn replay<T: Tool>(&self, derived: &[(&str, &Path)], setup: impl FnOnce(&Path)) -> Replay {
        let work = TempDir::new_in(std::env::temp_dir(), false).unwrap();
        let home = work.path().join("home");
        let cwd = work.path().join("cwd");
        fs::create_dir_all(&home).unwrap();
        fs::create_dir_all(&cwd).unwrap();
        setup(&cwd);
        let arguments: Vec<String> = self
            .mapped("argv.txt", &cwd, derived)
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect();
        assert_eq!(arguments[0], T::NAME, "{}", self.name);
        let mut sources = ToolRegistrySources::from_environment().unwrap();
        sources.prefixes.clear();
        sources.internal_tools_path = None;
        sources.ttd_internal_path = None;
        sources.file_context.user_override = Some(home);
        let registry = ToolHandler::new(sources);
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = run_with_registry::<T>(&arguments, &mut out, &mut err, &registry);
        let (out, took) =
            super::took_line::split_took_line(T::NAME, &String::from_utf8_lossy(&out));
        let log = fs::read_to_string(cwd.join("log.txt"))
            .ok()
            .map(|content| masked_log(&content));
        Replay {
            code,
            out,
            took,
            err: String::from_utf8_lossy(&err).into_owned(),
            log,
            cwd,
            _work: work,
        }
    }

    /// Replay the case and compare everything with the Release run: the exit
    /// code, standard output with the closing line present exactly where the
    /// Release build printed one, standard error, and the log file line for
    /// line (or its absence).
    pub fn assert_replayed<T: Tool>(&self, derived: &[(&str, &Path)], setup: impl FnOnce(&Path)) {
        let replay = self.replay::<T>(derived, setup);
        let release = self.release(T::NAME, &replay, derived);
        let name = self.name;
        assert_eq!(
            replay.code.as_i32(),
            release.exit,
            "{name}: exit\nstdout:\n{}\nstderr:\n{}",
            replay.out,
            replay.err
        );
        assert_eq!(replay.out, release.out, "{name}: stdout");
        assert_eq!(
            replay.took.is_some(),
            release.took.is_some(),
            "{name}: the closing line: port {:?}, Release {:?}",
            replay.took,
            release.took
        );
        assert_eq!(replay.err, release.err, "{name}: stderr");
        assert_eq!(replay.log, release.log, "{name}: log file");
    }
}

/// Compare a debug-level log file with the Release build's: the framework's
/// lines up to its read of `-no_progress` in order, then the rest as a
/// collection, because the tool body's `Value of … option` lines follow the
/// order in which the port reads its options and the input and output checks
/// run before the body (`docs/TOPP_CLI_SUPPORT.md`, *Log file and debug
/// levels*). `blocks` are runs of lines of the Release log, such as a
/// parameter dump, that must also appear in the port's log unbroken.
pub fn assert_debug_log(name: &str, actual: &[String], expected: &[String], blocks: &[&[String]]) {
    let framework = 1 + expected
        .iter()
        .position(|line| line.contains("Value of string option 'no_progress'"))
        .unwrap_or_else(|| panic!("{name}: no -no_progress line"));
    assert_eq!(
        actual[..framework],
        expected[..framework],
        "{name}: framework lines"
    );
    for block in blocks {
        assert!(!block.is_empty(), "{name}");
        assert!(
            expected.windows(block.len()).any(|w| w == *block),
            "{name}: the block is not the Release build's"
        );
        assert!(
            actual.windows(block.len()).any(|w| w == *block),
            "{name}: block missing or broken up:\n{}",
            block.concat()
        );
    }
    let (mut rest_actual, mut rest_expected) =
        (actual[framework..].to_vec(), expected[framework..].to_vec());
    rest_actual.sort();
    rest_expected.sort();
    assert_eq!(rest_actual, rest_expected, "{name}: the tool body's lines");
}

/// The lines of `log` from the one that ends with `header` through the
/// separator line that closes its parameter dump.
pub fn dump_block(log: &[String], header: &str) -> Vec<String> {
    let start = log
        .iter()
        .position(|line| line.trim_end().ends_with(header))
        .unwrap_or_else(|| panic!("no dump {header}"));
    let end = start
        + log[start..]
            .iter()
            .position(|line| line.starts_with(" - - - "))
            .unwrap_or_else(|| panic!("dump {header} not closed"));
    log[start - 1..=end].to_vec()
}

/// A log file with its `YYYY-MM-DD hh:mm:ss` time stamps masked, one entry
/// per line.
pub fn masked_log(content: &str) -> Vec<String> {
    content
        .split_inclusive('\n')
        .map(|line| {
            let bytes = line.as_bytes();
            let stamped = bytes.len() > 20
                && bytes[4] == b'-'
                && bytes[7] == b'-'
                && bytes[10] == b' '
                && bytes[13] == b':'
                && bytes[16] == b':'
                && bytes[19] == b' '
                && bytes[..4].iter().all(u8::is_ascii_digit);
            if stamped {
                format!("<time> {}", &line[20..])
            } else {
                line.to_owned()
            }
        })
        .collect()
}
