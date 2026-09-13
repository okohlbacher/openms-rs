# Detecting Python

`system::python_info` covers `SYSTEM/PythonInfo.h` and `SYSTEM/PythonInfo.cpp`
at SDK `bc9cc12`. Three static members.

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class PythonInfo` | the module itself | Three static methods become free functions. |
| `static bool canRun(std::string& python_executable, std::string& error_msg)` | `python_info::can_run(context, python_executable) -> Result<PythonCheck>` | The `bool` is `PythonCheck::can_run`. |
| `std::string& python_executable` (in/out) | `PythonCheck::executable` | The source rewrites the caller's argument to the resolved absolute path; the port returns it instead of mutating an input. |
| `std::string& error_msg` (out) | `PythonCheck::message` | Carries the informational "resolved to" note even on success, as the source does. |
| `static bool isPackageInstalled(const std::string& python_executable, const std::string& package_name)` | `python_info::is_package_installed(context, python_executable, package_name) -> Result<bool>` | |
| `static std::string getVersion(const std::string& python_executable)` | `python_info::version(context, python_executable) -> Result<String>` | Empty string on every failure, as the source. |

Native additions:

| Rust item | Why |
|---|---|
| `PythonCheck` fields `exit_code`, `timed_out` | The source conflates "no Python" and "too slow" into one `false`. |
| `version_from_output(stdout, stderr)` | The banner rule as a pure function, so it can be tested on a recorded string rather than a live interpreter — which is the only way to test it without Python installed. |
| `validate_package_name` | See *Native differences*; this is the security hardening. |
| `VERSION_ARGUMENT`, `COMMAND_ARGUMENT`, `MAX_PACKAGE_NAME_BYTES` | The source's literals `"--version"` and `"-c"`, named, and the new ceiling. |
| `context: &FileContext` parameter | Resolution through the crate's caller-owned runtime locations; see *Native differences*. |

## Preserved source conventions

**`canRun` resolves first, then runs.** The source calls
`File::findExecutable(python_executable)` (`PythonInfo.cpp:50`), which rewrites
the name to an absolute path, and only then executes `--version`. This keeps
that order, so a name that resolves nowhere never reaches a spawn.

**The message text is transcribed**, with its two-space indentation:

- `  Python not found at '<exe>'!` / `  Make sure Python is installed and this location is correct.`
- then, only for a relative name, the `PATH` advice with the current `PATH`, and
  on macOS the application-bundle paragraph;
- `  Python found at '<exe>' but failed to run!` / `  Make sure you have the rights to execute this binary file.`
- `  Python was found at '<exe>' but the process timed out …`

**`error_msg` is not an error channel.** The source appends
`Python executable ('<given>') resolved to '<resolved>'` whenever resolution
changed the string, and assigns `error_msg` on the success path too. That is
reproduced: `PythonCheck::message` can be non-empty while `can_run` is `true`.

**The "failed to run" branch is reachable exactly where the class test reaches
it**: an existing file that is not executable resolves, then fails to spawn.

**The banner rule is exact.** `getVersion` reads `--version` line by line,
appending each line to one string **without a separator**, standard output
first and then standard error — some interpreters print the version on
stderr — and then trims `' '`, `'\t'`, `'\n'` and `'\r'` from both ends
(`PythonInfo.cpp:155`–`163`, `StringUtils::trim`). The missing separator means a
multi-line banner comes back with its lines run together; `version_from_output`
reproduces that and the test asserts it, because a caller matching the string
would otherwise see a different value than the C++ produced.

**The 30 second budget applies to all three probes**, and
`isPackageInstalled` pipes both of the child's streams to the null device.

## Native differences

**`python -c "import <name>"` no longer interpolates arbitrary text.** The
source builds `bp::args({"-c", "import " + package_name})`
(`PythonInfo.cpp:117`). The argument vector is safe from the *shell*, but the
second argument is Python **source code**, so a caller-supplied
`package_name` is executed: `isPackageInstalled(py, "os; import os; os.system('…')")`
runs that command with the user's privileges. This port validates the name
first — a dotted sequence of ASCII identifiers, at most
`MAX_PACKAGE_NAME_BYTES` — and a name that fails validation returns `false`
**without starting an interpreter**. That is the same answer the source
reaches, since a string that cannot be a module path cannot name an installed
package; it simply reaches it without executing anything. Callers that want the
distinction call `validate_package_name`, which returns
`Error::InvalidValue`. Non-ASCII identifiers, which Python 3 permits, are
refused: the probe answers a yes/no question about ordinary package names, and
the narrower rule is deliberate. Recorded in `OpenMS_CPP_ISSUES.md` as
*`PythonInfo::isPackageInstalled` executes its argument*.

**The resolved path is executed directly.** After `File::findExecutable` made
the name absolute, the source hands that absolute path back to
`bp::search_path`, which only ever walks `PATH` entries — so a perfectly good
absolute interpreter fails the probe whenever `PATH` is empty or unset, and the
class test cannot notice because its positive block is guarded by `canRun`
itself. This port runs the resolved path. See `OpenMS_CPP_ISSUES.md` entry
*empty `PATH` defeats an absolute executable*.

**Resolution goes through `FileContext`**, not the process `PATH` directly, for
the same reason as `system::java_info`: it is what makes the module testable
with no Python installed.

**`PATH` is read on every call, and an unset `PATH` is empty**, where the source
caches `getenv("PATH")` in a function-local `static` and constructs a
`std::string` from a null pointer if it is unset (`PythonInfo.cpp:59`).

**Nothing is printed.** The diagnosis is returned.

## Checked boundaries and evidence

| Boundary | Behaviour |
|---|---|
| Name that resolves nowhere | `can_run == false`, message contains `Python not found at`. |
| Existing but non-executable file | `can_run == false`, message contains `failed to run`. |
| Interpreter present | `can_run == true`, `executable` absolute and existing, message holds the "resolved to" note. |
| Probe exceeds 30 s | child killed, `timed_out == true`, `can_run == false`. |
| `is_package_installed` with a name that is not a module path | `Ok(false)`, no interpreter started. |
| `is_package_installed` with an absent interpreter | `Ok(false)`. |
| `version` with an absent, unstartable, timed-out or failing interpreter | `Ok(String::new())`. |
| Executable name with a NUL byte, or over the path ceiling | `Error::InvalidValue` before anything is searched or spawned. |

Evidence is **tier 3** for the class test's four literals — `canRun` false with
`"Python not found at"`, `canRun` false with `"failed to run"`,
`isPackageInstalled(py, "veryWeirdPackage___@@__@") == false` and
`isPackageInstalled(py, "math") == true` — and for the transcribed message and
banner rules. It is **tier 4** for the banner cases derived from
`std::getline` + `StringUtils::trim` semantics rather than transcribed
(`"Python 3.11.4\ntrailing\n"` → `"Python 3.11.4trailing"` is derived, not
observed), for the validation rule, and for the stand-in-interpreter tests. No
C++ was executed.

Tests: `tests/system_process.rs`
(`python_can_run_reproduces_the_two_message_literals`,
`python_can_run_resolves_a_bare_name_to_an_absolute_path`,
`python_package_detection_answers_the_class_tests_two_values`,
`python_package_detection_uses_the_stand_in_interpreter`,
`python_version_parsing_follows_the_sources_concatenate_then_trim`,
`python_version_is_read_from_the_stand_in_and_empty_when_absent`) and
`src/system/python_info.rs::tests`.
