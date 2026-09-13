# Calling R scripts

`system::r_wrapper` covers `SYSTEM/RWrapper.h` and `SYSTEM/RWrapper.cpp` at SDK
`bc9cc12`. Three static members: find the interpreter, find the script, run it.

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class RWrapper` | the module itself | Three static methods become free functions. |
| `static std::string findScript(const std::string& script_file, bool verbose = true)` | `r_wrapper::find_script(context, script_file) -> Result<PathBuf>` | `Exception::FileNotFound` → `Error::Io` with `ErrorKind::NotFound`, this crate's counterpart. |
| `static bool findR(const std::string& executable = "Rscript", bool verbose = true)` | `r_wrapper::find_r(context, executable) -> Result<RCheck>` | The `bool` is `RCheck::found`. |
| `static bool runScript(const std::string& script_file, const std::vector<std::string>& cmd_args, const std::string& executable = "Rscript", bool find_R = false, bool verbose = true)` | `r_wrapper::run_script(context, script_file, cmd_args, executable, find_r_first) -> Result<ScriptOutcome>` | The `bool` is `ScriptOutcome::success`. |
| `executable = "Rscript"` (default argument, twice) | `r_wrapper::DEFAULT_EXECUTABLE` | Rust has no default arguments; callers pass the constant. |
| `bool verbose` (three times) | — | **Not ported**: it only chose whether to write to `OPENMS_LOG_INFO` / `OPENMS_LOG_ERROR`. The text is always returned in `RCheck::message` / `ScriptOutcome::message`, or in the error for `find_script`, and the caller decides where it goes. |

Native additions:

| Rust item | Why |
|---|---|
| `RCheck` with `found`, `output`, `message`, `executable` | The merged interpreter output the source assembles for its error message, and the path it resolved and threw away. |
| `ScriptOutcome` with `success`, `script`, `stdout`, `stderr`, `message` | The two drained buffers the source builds and only prints. `script` says which file was actually run. |
| `reassemble_lines` | The source's drain loop as a pure function; see *Preserved source conventions*. |
| `PROBE_ARGUMENTS`, `SCRIPT_ARGUMENTS`, `SCRIPTS_DIRECTORY`, `MAX_SCRIPT_ARGUMENTS` | The source's literals, named, plus the new ceiling. |
| `context: &FileContext` parameter | Resolution and shared-data lookup through the crate's caller-owned runtime locations. |

## Preserved source conventions

**The two argument vectors are exact**: `--vanilla -e sessionInfo()` for the
interpreter probe (`RWrapper.cpp:136`), and `--vanilla --quiet <script> <args…>`
for a script run (`RWrapper.cpp:64`–`70`). Both are `pub const`, so a reader can
see them without opening the implementation.

**`findScript` resolves the OpenMS data path first.** The source is
`File::getOpenMSDataPath()`, `ensureLastChar('/')`, then
`File::find(script_file, [path + "SCRIPTS"])` (`RWrapper.cpp:213`), all inside
one `try` whose `catch (...)` rethrows `FileNotFound`. That order is kept: an
existing absolute script still fails when the shared data cannot be located,
because the port asks for the data path before it asks for the file — which is
what the source does, and changing it would change which callers succeed.

**`runScript` swallows the not-found error.** `findScript` is the only entry
point that fails loudly; `runScript` catches it and returns `false`
(`RWrapper.cpp:52`–`57`). The class test pins exactly that contract, and
`run_script_degrades_to_false_without_raising` reproduces it.

**The check order is the header's**: optional `findR`, then `findScript`, then
the run. A bogus interpreter with `find_r_first` never reaches the script
lookup, so `ScriptOutcome::script` stays `None`.

**Both pipes are drained before the wait**, which is why the source starts a
`std::thread` for standard error: a child that fills a pipe buffer while the
parent waits would deadlock. This port's reader threads serve the same purpose.

**Output is reassembled the way the source drains it.** The source reads with
`std::getline` and appends `line + "\n"` (`RWrapper.cpp:90`, `:158`), so a final
line without a newline gains one and an empty stream stays empty.
`reassemble_lines` reproduces that, so `RCheck::output` and
`ScriptOutcome::stdout` / `stderr` hold the strings the C++ assembled.

**`findR` merges the two streams in the source's order**, standard output then
standard error, imitating Qt's `MergedChannels` as the comment says.

**The failure messages are transcribed**, including the layout of the probe
failure (`Output was:` / `------>` / `<------`) and of the script failure
(`--- ERROR MESSAGES ---`, `--- OTHER MESSAGES ---`, and the trailing
`Script failed. See above for an error description. ` with its trailing space).

**There is no time budget.** Unlike `JavaInfo` and `PythonInfo`, neither
`findR` nor `runScript` uses `Internal::waitForProcess`; a plotting script may
legitimately run for a long time. That is preserved, and the doc comments say
so, with `Invocation::timeout` available to a caller that wants one.

## Native differences

**The failure message names the interpreter that was actually requested.** The
source hardcodes the literal `Rscript` in `findR`'s error text — both in
`'Rscript' executable returned with error` and in `(command: 'Rscript …')` —
even when `executable` is something else. This interpolates `executable`, so
the diagnosis matches the command that ran. The sentence structure is otherwise
unchanged.

**Resolution goes through `FileContext`**, not the process `PATH` directly, so
both `find_r` and `run_script` can be exercised against a recorded stand-in
interpreter with no R installed. `run_script` resolves the interpreter itself
even when `find_r_first` is false, reaching the same
`Could not run '<exe>'. Is it installed and in PATH?` message the source
produces from a caught `process_error`.

**Arguments are refused above `MAX_SCRIPT_ARGUMENTS`** before the interpreter is
looked up, and `Invocation::preflight` additionally rejects an interior NUL. The
source passes the caller's vector through unbounded.

**Nothing is printed**; every diagnosis is returned.

**Threading.** `RWrapper.cpp` has no `#pragma omp`; its one `std::thread` exists
to avoid a pipe deadlock, not for parallelism, and the port's two reader threads
serve the same purpose.

## Checked boundaries and evidence

| Boundary | Behaviour |
|---|---|
| `find_r` with an interpreter that does not exist | `found == false`, `executable == None`, message names `FailedToStart`; never an error. |
| `find_r` with a present interpreter that exits non-zero | `found == false`, message quotes the merged output. |
| `find_script` with a missing script | `Error::Io` / `ErrorKind::NotFound`. |
| `find_script` with no locatable OpenMS shared data | the data-path error, as the source's `catch (...)` → `FileNotFound`. |
| `run_script` with `find_r_first` and a bogus interpreter | `success == false`, `script == None`, no error. |
| `run_script` with a missing script | `success == false`, message `Could not find R script '<x>'!`. |
| `run_script` with a script that exits non-zero | `success == false`, message holds both drained buffers. |
| `cmd_args.len() > MAX_SCRIPT_ARGUMENTS` | `Error::InvalidValue`; nothing is looked up or spawned. |
| An argument containing spaces | passed as one argument; asserted against the stand-in interpreter. |

Evidence is **tier 3** for the class test's assertions — no exception escapes
any of the three entry points, `findR("this_is_not_a_real_R_interpreter_xyz")`
is `false`, `findScript("definitely_nonexistent_script_qwerty.R")` raises
not-found, and both `runScript` degradation paths return `false` — and for the
transcribed argument vectors and messages. It is **tier 4** for the stand-in
interpreter runs, which derive the expected argument vector from the source
rather than from a running R, and for the ceilings. No C++ was executed; the
header has no retained output fixture.

Tests: `tests/system_process.rs` (`find_r_reports_a_bogus_interpreter_cleanly`,
`find_r_accepts_a_recorded_stand_in_interpreter`,
`find_script_reports_a_missing_script_as_not_found`,
`find_script_resolves_a_bundled_script_under_the_scripts_directory`,
`run_script_degrades_to_false_without_raising`,
`run_script_passes_the_sources_argument_vector_intact`,
`captured_output_is_reassembled_the_way_the_source_drains_it`,
`run_script_refuses_an_oversized_argument_list`) and
`src/system/r_wrapper.rs::tests`.
