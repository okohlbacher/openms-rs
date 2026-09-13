# Lexical path helpers

`system::path_utils` covers `SYSTEM/PathUtils.h` at SDK `bc9cc12`. The header is
54 lines, inline-only and has no `.cpp`, so the header *is* the implementation.
It exists for one reason: to give a caller a basename and a UTF-8 safe path
conversion without dragging in `File.h`.

## API mapping

Every public member of the header, with its Rust counterpart.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `PathUtils::basename(const std::string& file)` | [`system::path_utils::basename`](../src/system/path_utils.rs) | Delegates to `system::file::basename`, which the source itself duplicates from `File::basename`. Borrowed `&str` in, borrowed `&str` out — no allocation. |
| `to_path(const std::string& s)` | [`system::path_utils::to_path`](../src/system/path_utils.rs) | Returns `Result<PathBuf>`; the extra failure paths are the port's own bounds, not a source behaviour. |
| `MAX_PATH_BYTES` | `system::path_utils::MAX_PATH_BYTES` | Native only; aliases `system::file::MAX_PATH_BYTES` so the two agree. |

There is no unported member.

## Preserved source conventions

Both entry points stay **purely lexical**. Neither touches the filesystem,
normalises `.`/`..` segments, collapses repeated separators, or resolves
symbolic links, and neither cares whether the path exists.

`basename` reproduces all three outcomes of
`file.substr(file.find_last_of("\\/") + 1)`:

| Input | Result | Why |
|---|---|---|
| `"/data/run1.mzML"` | `"run1.mzML"` | everything after the last separator |
| `"run1.mzML"` | `"run1.mzML"` | no separator: the source reaches this through `npos + 1` wrapping to zero, a quirk its own comment calls out |
| `"/data/"`, `"/"`, `""` | `""` | a trailing separator leaves nothing after it |

Both slash styles are accepted on every platform, not only on Windows, and a
path that mixes them splits at the last one of either.

## Native differences

**`to_path` solves a problem Rust does not have.** The source constructs from
`std::u8string` because `std::filesystem::path(std::string)` interprets its
argument in the active Windows code page rather than as UTF-8, and it catches
`std::system_error` to fall back to that code page when the bytes are not valid
UTF-8 — a filename taken from `argv` under an ANSI locale. A Rust `&str` is
UTF-8 by construction and `PathBuf` stores it as UTF-8 on Unix and WTF-8 on
Windows, so the re-encoding is a no-op and **the fallback branch is
unreachable**: the byte sequence that triggers it cannot be held in a `&str`.
What survives is the invariant the source test actually asserts, that the
conversion round-trips byte for byte, and that is tested here.

**`basename` is not duplicated.** The source keeps two copies of the same
function so that `PathUtils.h` need not include `File.h`; Rust modules have no
such cost, so `path_utils::basename` calls `file::basename` and the two cannot
drift apart.

**`to_path` returns `Result`.** The source cannot fail outside its Windows
branch. This port refuses a path longer than `MAX_PATH_BYTES` (1 MiB) or
containing a NUL byte, which is the same policy `system::file` applies: a NUL
cannot survive any platform filesystem call, so one checked error at the
conversion replaces a failure at every later use.

## Checked boundaries and evidence

| Boundary | Value | Source behaviour |
|---|---|---|
| Path length | 1 MiB (`MAX_PATH_BYTES`) | unbounded |
| Interior NUL | refused | unchecked |

The class test has a **single section**, `to_path`, and all six of its
cross-platform round-trip literals — ASCII, `U+00E4`, the CJK string 日本語, a
name with spaces and parentheses, a 300-character name and the mixed
`测试 äö` — are transcribed into `tests/path_utils.rs`. Its three `_WIN32`-only
cases (`PathUtils_test.cpp:74-127`) are not ported: two exercise the ANSI
fallback for the lone byte `0xE4`, which no `&str` can hold — once on a bare
filename and once inside a UNC path — and the third is a plain ASCII UNC path
the port accepts with no special handling. `basename` has no section at all, so
its cases are derived from the source expression and are tier 4.

No C++ execution is claimed. Source hashes, line anchors and the class-test
review are in [the provenance record](../tests/data/path_utils_provenance.json).
