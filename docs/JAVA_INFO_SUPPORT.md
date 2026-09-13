# Detecting Java

`system::java_info` covers `SYSTEM/JavaInfo.h` and `SYSTEM/JavaInfo.cpp` at SDK
`bc9cc12`. The header declares one static member.

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `class JavaInfo` | the module itself | A class with a single static method becomes a free function. |
| `static bool canRun(const std::string& java_executable, bool verbose_on_error = true)` | `java_info::can_run(context, java_executable) -> Result<JavaCheck>` | The `bool` is `JavaCheck::can_run`. |
| `bool verbose_on_error` | — | **Not ported**: it only chose whether to write the diagnosis to `OPENMS_LOG_ERROR`. There is no global log stream here, so the text is always returned in `JavaCheck::message` and the caller decides. Rust has no default arguments, so the source's `= true` has nothing to carry. |

Native additions:

| Rust item | Why |
|---|---|
| `JavaCheck` with `can_run`, `message`, `exit_code`, `timed_out`, `output` | Returns what the source logs and discards. `exit_code` is the number the source interpolates into its own message; `timed_out` separates "no Java" from "busy machine", which the source's single `false` conflates; `output` is the banner both pipes collected. |
| `VERSION_ARGUMENT` | The source's literal `"-version"`, named. |
| `context: &FileContext` parameter | See *Native differences*. |

## Preserved source conventions

**The command is `java -version` with both pipes attached**, under the 30 second
budget of `Internal::waitForProcess`, and success is exit code 0.

**The diagnosis text is transcribed**, including its structure of a
`Java-Check:` header line followed by two-space-indented lines, and including
all three branches:

- not found → `Java not found at '<exe>'!` / `Make sure Java is installed and this location is correct.`
- then, for a *relative* name, the `PATH` advice and the current `PATH`; for an
  *absolute* one, `You gave an absolute path to Java. Please check if it's correct.`
- timed out → `Java was found at '<exe>' but the process timed out …` with the
  `'force' flag` advice.
- non-zero exit → `Error executing '<exe>'!` / `Java returned a non-zero exit code (<n>).`

**The macOS-only paragraph about application bundles changing the system `PATH`
is kept**, under `cfg!(target_os = "macos")` where the source has `#ifdef __APPLE__`.

**Relative versus absolute is decided by `PathUtils::to_path(…).is_relative()`**,
the same predicate the source uses, through `system::path_utils::to_path`.

## Native differences

**Resolution goes through `FileContext`.** The source calls
`bp::search_path(java_executable)` directly, which reads the process `PATH` and
nothing else. `can_run` takes a `&FileContext` and uses
`FileContext::find_executable`, which is the crate's established caller-owned
runtime-locations pattern and is what makes this module testable with no JVM
installed: a test supplies a `search_path` pointing at a recorded stand-in.

Two behavioural consequences, both improvements, both documented at the item:

- an *absolute* Java is found even when `PATH` is empty or unset, where the
  source's `search_path` walks only `PATH` entries and therefore fails;
- resolution accepts a relative path against the current directory, matching
  the header's own "can be absolute, relative or just a filename".

**`PATH` is read on every call and an unset `PATH` is empty.** The source has
`static std::string path; if (path.empty()) path = getenv("PATH");`
(`JavaInfo.cpp:98`), which constructs a `std::string` from a null pointer when
`PATH` is unset — undefined behaviour — and caches the first value it ever saw.
Reported for the shared [C++ issue ledger](../OpenMS_CPP_ISSUES.md) as
*`getenv("PATH")` assigned into `std::string`*; the integrator owns that file
and assigns the `CPP-` number, so this package cites the defect by title rather
than by an identifier it cannot mint.

**The message is returned, not printed.** Nothing in this crate writes to a
global log stream.

**A timeout is distinguishable.** The source returns `false` for both "not
there" and "too slow"; `JavaCheck::timed_out` separates them while `can_run`
keeps the source's single answer.

## Checked boundaries and evidence

| Boundary | Behaviour |
|---|---|
| `can_run(ctx, "")` | `can_run == false`, message `Java not found at ''!` — the class test's only assertion. |
| Relative name not on the search path | `false`, message names the current `PATH`. |
| Absolute name that does not exist | `false`, message gives the absolute-path advice and does *not* name `PATH`. |
| Java present, exit code 0 | `true`, empty message, `output` holds the banner (Java writes it to standard error). |
| Java present, non-zero exit | `false`, message quotes the code. |
| Probe exceeds 30 s | child killed, `timed_out == true`, `can_run == false`. |
| Executable name with a NUL byte, or over the path ceiling | `Error::InvalidValue` before anything is searched or spawned. |

Evidence is **tier 3** for `canRun("") == false`, the one literal the class test
offers, and for the transcribed message text; **tier 4** for everything else,
including the stand-in-interpreter tests, which derive the expected answers from
the source's own control flow rather than from a running JVM. No C++ was
executed.

Tests: `tests/system_process.rs` (`java_cannot_run_from_an_empty_executable_name`,
`java_is_probed_through_a_recorded_stand_in_on_the_search_path`) and
`src/system/java_info.rs::tests`.
