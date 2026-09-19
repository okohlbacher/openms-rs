# FileInfo tool support

Native coverage of the FileInfo TOPP tool, `OpenMS4-topp/src/FileInfo.cpp` at
topp `174b576e244e100f2345ca57a8e79aaa607156df`, on top of the TOPPBase
lifecycle (`docs/TOPP_CLI_SUPPORT.md`) and the FileInfo library
(`docs/FILE_INFO_SUPPORT.md`). Work package A5-FILEINFO-TOOL of the early TOPP
bundle, stage 1 of the FileInfo preview; `-i`, `-d` and `-c` are A6, which has
landed — see [FILE_INFO_CHECKS_SUPPORT](FILE_INFO_CHECKS_SUPPORT.md) — the
consensusXML, identification and FASTA branches A7, and `-v`, mzXML, mzData and
trafoXML A8.

| Artifact | Path |
| --- | --- |
| Tool | `src/cli/tools/file_info.rs` |
| Executable | `src/bin/FileInfo.rs` (`[[bin]]`, `required-features = ["mzml", "paramxml", "featurexml"]`) |
| Report branches | `src/format/file_info/` (A4), with `text_format.rs` (A2) |
| Tests | `tests/topp_file_info.rs` (33 cases, none ignored) |
| Fixtures | `tests/data/topp_file_info/` (`inputs/`, `expected/`), plus A4's `tests/data/file_info/` read in place |
| Manifest | `tests/data/topp_file_info_provenance.json` |
| Oracle drivers | `../oracle/topp-early-bundle/` (C1) and `../oracle/topp-file-info-tool/` (this package), outside this repository |

The tool is thin, as the source is: it registers the parameters, resolves the
input type, maps the flags onto `format::file_info::model::Options`, and routes
the library's two reports. No report branch lives here. It adds no module edge:
`cli -> format` and `cli -> system` already exist.

## Capability table

`-in` accepts, and `-in_type` names, the source's seventeen types. What each one
does here:

| Input | Flags | Result |
| --- | --- | --- |
| `dta`, `dta2d`, `mzML` | none, `-m`, `-p`, `-s` (any combination) | full report to `-out` or the output stream, TSV to `-out_tsv` |
| `featureXML` | none, `-m`, `-p`, `-s`, and `-d`/`-c`, which the source ignores there | full report and TSV |
| `dta`, `dta2d`, `mzML` | `-d`, `-c`, in any combination with the above | the detailed listing and the corrupt-data check in the report, A6 |
| any | `-v` | exit 11, `... schema and semantic validation (-v) is not ported` (A8) |
| `mzML` | `-i`, valid index | the index line, then the content; exit 0 (A6). Its counts are this port's decoder's — native difference 9 |
| `mzML` | `-i`, no index | the failure text and nothing after it; exit 6 (A6). Which files have no index is also this port's decoder's answer — native difference 9 |
| non-mzML | `-i` | exit 6 with the source's message and the usage text, before the library runs |
| `mzXML`, `mzData`, `mgf`, `sqMass`, `fid` | any | exit 11, `... peak-file branch for <type> input is not ported` |
| `consensusXML`, `idXML`, `mzid`, `pepXML`, `mzTab`, `trafoXML`, `fasta`, `pqp` | any | exit 11, `FileInfo <type> branch is not ported` |
| undetermined type | any | exit 10, `Error: Could not determine input file type!` |

Every refusal is explicit and writes no report: the C++ tool reports these
branches and exits 0 (decisions D3 and D5(c) of the early TOPP bundle map
`Error::Unsupported` to `INCOMPATIBLE_INPUT_DATA`).

## API mapping

Every member of the source `TOPPFileInfo` and the behaviour it carries.

| C++ | Rust |
| --- | --- |
| `class TOPPFileInfo : public TOPPBase` | `cli::tools::FileInfo`, a unit struct implementing `cli::Tool` |
| `TOPPFileInfo()` with name and description | `Tool::NAME`, `Tool::DESCRIPTION` |
| `registerOptionsAndFlags_()` | `Tool::register`, with the seventeen types in `INPUT_TYPES` (file-private) |
| `main_(int, const char**)` | `Tool::run_io`, which the framework calls; `Tool::run` forwards to it with the process streams |
| `outputTo_(ostream&, ostream&)` | the body of `run_io` after the two destinations are opened |
| `ofstream os` / `os_filt` | `-out` opened with `File::create`, or the `out` writer of `run_with` |
| `ofstream os_tsv` / `boost::iostreams::null_sink` | `-out_tsv` opened the same way, or nothing |
| `getGlobalLogInfo()` as the default report sink | the `out` writer, so a `run_with` caller captures the report |
| `FileTypes::nameToType(getStringOption_("in_type"))` | `FileType::from_name` |
| `FileHandler::getType(in)` | `detect_type`, `FileHandler::get_type` with a directory reported as unknown |
| `writeDebug_("Input file type: …", 2)` | not ported; `docs/TOPP_CLI_SUPPORT.md` records the debug helpers as unported |
| `writeLogError_` | `writeln!(err, …)` |
| `printUsage_()` | `cli::print_usage` with this tool's spec |
| `FileInfo::Options` assignment (eight members plus `log_type_`) | the `Options` literal in `run_io`, plus the native `source_dangling_references` |
| `FileInfo fi; fi.run(in, opt)` | `format::file_info::report::FileInfo::new().run(&input, &options)` |
| `os << FileInfo::toText(r)`, `os_tsv << FileInfo::toTSV(r)` | `write_all` of `to_text` and `to_tsv` |
| `if (opt.check_index && r.validation.index_checked && !r.validation.index_valid) return ILLEGAL_PARAMETERS;` | the same condition, reached since A6; TOPP_FileInfo_11's `WILL_FAIL` exit code is reproduced |
| `PARSE_ERROR`, `ILLEGAL_PARAMETERS`, `EXECUTION_OK` | `ExitCode::ParseError`, `IllegalParameters`, `ExecutionOk` |
| `throw Exception::FileNotWritable(...)` in `main_` | `open_output`'s `UNKNOWN_ERROR` with the source's `BaseException` wording; the source path is unreachable after `outputFileWritable_`, this one is reached for a directory (oracle `out_is_directory`) |
| `int main(int argc, const char** argv)` | `src/bin/FileInfo.rs`, three lines around `cli::run::<FileInfo>()` |

## Preserved source conventions

- **Registration** (`FileInfo.cpp:83-100`): `-in` required with the seventeen
  valid formats; `-in_type` a free string restricted to the same seventeen
  names; `-out` optional with format `txt`; `-out_tsv` optional, advanced, with
  format `tsv`; the flags `m`, `p`, `s`, `d`, `c`, `v`, `i` with their verbatim
  descriptions, in that order, before the common TOPP options. `--help`,
  `--helphelp` and `-write_ini` reproduce the C++ text and file
  (`usage_text_matches_the_cpp_text`, `write_ini_matches_the_cpp_file`).
- **Output order** (`main_`): `-out` and then `-out_tsv` are opened before
  anything is read, so an early failure leaves them empty and truncates an
  existing file (`c1_unknown_type_exits_10`,
  `an_existing_out_is_truncated_by_a_failing_run`, `c1_truncated_inputs_exit_3`).
- **Report routing**: without `-out` the report goes to the info log, here the
  `out` writer, and the TSV is discarded without `-out_tsv`
  (`without_out_the_report_goes_to_the_output_stream`, `out_tsv_without_out`).
- **Type resolution** (`outputTo_:105-118`): `-in_type` first, then
  `FileHandler::getType`, which reads the file name and then the content, so a
  featureXML map named `.tmp` runs the featureXML branch with or without
  `-in_type` (`a_featurexml_map_under_a_tmp_name`). An undetermined type is
  `PARSE_ERROR` (10) after the framework's `Warning: Could not determine format
  of input file '<path>'!`.
- **`-i` outside mzML** (`:119-124`): the error line, then the usage text on the
  error stream, then `ILLEGAL_PARAMETERS` (6), before the library runs.
- **The resolved type is forced** into the library options, so the loader
  accepts only that type.
- **`log_type`** is `CMD` unless `-no_progress`, as `TOPPBase::log_type_`; the
  native loaders report no progress, so it changes no output.
- **The timing line** (`TOPPBase::main`) never reaches `-out`: it is not ported
  at all, so the output stream carries the report alone.
- **Serial**: FileInfo has no parallel section in the source and none here.
  `-threads` is accepted and changes nothing (`threads_do_not_change_the_reports`).

## Native differences

1. **Unported branches and flags are refused** with
   `Error::Unsupported`, exit 11, instead of the source's report and exit 0; see
   the capability table. The message names the branch or the flag. Nothing is
   written to the reports, and `-out` stays empty
   (`unported_branches_are_refused_explicitly`).
2. **TOPP_FileInfo_9's registered input does not load.** Three mzML reader gaps
   outside this package (duplicate spectrum-level `userParam name`, a
   `dataProcessingRef` on primary arrays, a 64-bit float `charge array`;
   `docs/FILE_INFO_SUPPORT.md`, *Known reader gaps*) make the tool exit 11 where
   C++ exits 0. The registration is therefore run on the derived
   `FileInfo_9_strict_reader.mzML`, whose C++ report is identical apart from the
   file name, and the original input's refusal is asserted as a tripwire
   (`topp_file_info_9_registered_input_is_refused_by_the_mzml_reader`).
   FileInfo_12's input hits the same `charge array` gap. Since A6 its *index*
   is checked and parses with the counts the C++ reports, but the content after
   it still cannot be read, so TOPP_FileInfo_12's exit code stays out of reach
   and the loadable `-i` cases use the core `IndexedmzMLFile_1` fixture.
3. **A forced type the file contradicts.** The source's loaders detect the type
   themselves and throw when it is not the forced one: `ParseError`, exit 3, from
   `loadExperiment` (oracle `forced_dta_on_featurexml`,
   `forced_mzml_on_featurexml`) and `InvalidFileType`, exit 8, from `loadFeatures`
   (oracle `forced_featurexml_on_dta`). The library maps the first to
   `Error::InvalidValue`, exit 6, and the featureXML branch, which does not pass
   the forced type to the feature loader, reaches the loader's own refusal,
   exit 11. Both are recorded in `a_forced_type_the_file_contradicts_is_refused`;
   closing them belongs to the library (integrator request 2).
4. **Validation order.** The framework checks every option before the tool body,
   where the source reads `-in` inside `outputTo_`, after `-out` was opened. A
   run with a missing or empty `-in` therefore leaves no `-out` file here, while
   the C++ tool leaves an empty one (oracle `zero_byte_input`, C1
   `FileInfo_missing_in`); the exit codes agree (4 and 7).
5. **The dangling-reference warning.** The tool reads mzML with the
   source-compatible dangling-reference policy (decision D10), which writes one
   warning per distinct dangling ID to the crate's warning log stream, that is
   the process's standard error, not the tool's error stream. The C++ reader is
   silent. `c1_empty_mzml_with_a_dangling_reference` shows the report is the
   C++ one.
6. **Usage text.** As for every ported tool: no colours, no console-width
   shaping (the source's `stty` probe line is absent) and the version line names
   this crate's pinned core revision (`docs/TOPP_CLI_SUPPORT.md`).
7. **`-out` that cannot be opened** is `UNKNOWN_ERROR` (8) with the source's
   `BaseException` wording, matching the executed C++ for a directory; the
   source's own `FileNotWritable` throw is unreachable there because
   `outputFileWritable_` ran first.
8. **`writeDebug_`** of the detected type at debug level 2 is not ported.
9. **`-i` answers with this port's index decoder, which departs from the source
   in two measured places.** Whether a file has an index, and the spectrum and
   chromatogram counts printed on the line above the content, come from
   `src/format/indexed_mzml.rs`, not from a re-implementation of
   `IndexedMzMLDecoder`. Two boundaries of that decoder differ, and both
   departures are the owning package's own reviewed decisions, recorded in
   `docs/INDEXED_MZML_SUPPORT.md`:
   *a file shorter than 1023 bytes*, where the C++ seeks back past the start of
   the file and its regex searches an uninitialised heap buffer
   (`IndexedMzMLDecoder.cpp:165-168`, `:179-181`), so it reports no index and
   this tool's C++ counterpart exits `ILLEGAL_PARAMETERS`, while the port
   searches the file it has and finds the index that is there; and
   *an `<index>` element written without whitespace*, whose first `<offset>`
   the C++ DOM walk skips (`:280-282` and `:290-293`, upstream `CPP-337`), so
   every section is counted one short — a section holding a single offset is
   dropped entirely — while the port keeps it. Reproducing either would mean
   changing a decoder every index reader in the crate shares, and the second
   would break random access, because a dropped offset is a record that cannot
   be found. The departure therefore stands, and `-i`'s counts are not
   bit-identical to the source's on those two classes of file. No file an
   OpenMS writer produces is in either class: both indexed inputs this package
   runs on, `IndexedmzMLFile_1.mzML` and `FileInfo_12_input.mzML`, write each
   `<offset>` on its own indented line and are far longer than the search
   window, which is why the source's two boundaries stay invisible in practice.
   Measured against the Release build and pinned by the oracle cases
   `i_window_below`, `i_window_above` and `i_offsets_unspaced`; native
   differences 11 and 12 of `docs/FILE_INFO_CHECKS_SUPPORT.md` carry the full
   record, with `src/format/indexed_mzml.rs` named as the owner.

## Checked boundaries and evidence

Tier 1, retained upstream outputs (test-data `0cb15f2`,
`topp/CMakeLists.txt:881-904`): TOPP_FileInfo_1, _2, _3 and _9 run through the
tool with their registered arguments into a per-case temporary directory, and
the produced `-out` file compared with the retained expectation through C3's
`FuzzyStringComparator` with the pinned `FuzzyDiff.ini` (ratio 1.01, absdiff
0.01) and the registered whitelist `File name`, as the `_out1` comparisons do.
_9 runs on the derived input (native difference 2).

Tier 1, executed differential (product-SDK FileInfo, Debug, core `4fdec46`,
decision D7). The C1 oracle (`../oracle/topp-early-bundle`, manifest
`64e98543…`, run1) supplies the reports of FileInfo_1/2/3/9 with `-out_tsv`, the
empty featureXML and mzML, both FeatureFinderCentroided_1 outputs and the
`rep.dat` run, plus the exit codes and diagnostics of the unknown type, `-i` on
DTA, an invalid `-in_type`, wrong `-out`/`-out_tsv` extensions, an unwritable
`-out`, a missing input file, a missing `-in`, the two truncated inputs and the
bare invocation. This package's oracle (`../oracle/topp-file-info-tool`, both
runs reproduced) adds `--help`, `--helphelp`, `-write_ini`, an INI written by
C++ and fed back with the flags set, the two stdout routings, three forced-type
mismatches, a directory as `-out` and as `-in`, a `.tmp` featureXML forced and
detected, `-threads 4`, `-d -c` on featureXML, a zero-byte input and the
truncation of a pre-filled `-out`. Text and TSV are compared byte for byte with
only the `File name` text line and the `general: file name` TSV line normalised,
because they hold the path the run was given.

Tier 4: the refusal, message and empty `-out` of every unported branch and flag,
including all upstream registrations A7 and A8 will close; the executable's exit
status and its report on the process's standard output; the strict library
default behind the tool's source-compatible mzML reading.

The `File name` normalisation is the only tolerance used on a report; no
numeric tolerance is applied outside the registered FuzzyDiff comparisons.

## Deferrals

- `-i`, `-d` and `-c` are no longer deferred: A6 ported them and closed
  TOPP_FileInfo_11 and _19. TOPP_FileInfo_12 stays open on the `charge array`
  reader gap, not on the flag.
- consensusXML, idXML/mzIdentML and FASTA (A7): TOPP_FileInfo_7, _10, _13, _17,
  _18 and _20.
- `-v`, mzXML, mzData and trafoXML (A8): TOPP_FileInfo_4, _5, _6, _14, _15 and
  _16. TOPP_FileInfo_15's input is not retained here; it differs from _14's only
  inside the document and takes the same `-v` refusal.
- The three mzML reader gaps of native difference 2, and the forced-type
  refusal classes of native difference 3, both outside this package.
- `-log`, `-instance` and the timing line stay as `docs/TOPP_CLI_SUPPORT.md`
  records them for every tool.
