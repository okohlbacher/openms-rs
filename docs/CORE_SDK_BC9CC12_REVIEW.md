# Core SDK refresh to bc9cc12

The target advances from `82ce5b373c97f934ffd9b1ffd80215ca66473d0b` to
[`bc9cc12514c768385ce121d6ca4bb710fe1983c4`](https://github.com/okohlbacher/OpenMS4-core/tree/bc9cc12514c768385ce121d6ca4bb710fe1983c4),
fetched from the repository's default `codex/package-split` branch into a clean
detached worktree on 2026-09-11. The previous checkouts and all original fixture
pins are retained.

The [upstream comparison](https://github.com/okohlbacher/OpenMS4-core/compare/82ce5b373c97f934ffd9b1ffd80215ca66473d0b...bc9cc12514c768385ce121d6ca4bb710fe1983c4)
contains eight commits and eleven changed paths. Six commits are packaging and
CI only: a Homebrew tap and formula for the Core SDK, exporting Homebrew OpenMP
to consumers, testing the exported SDK target, and reporting the real formula
version. Two touch source.

## Runtime compatibility

Exactly **two** files inside the comparable scientific roots changed, both in
the Parquet reader:

- `src/openms/include/OpenMS/FORMAT/ParquetFile.h`
- `src/openms/source/FORMAT/ParquetFile.cpp`

`ParquetFile::readTable(const std::string&)` now delegates to the existing
`RandomAccessFile` overload and closes the file it owns before returning,
including on the error path, where it preserves the original exception. The
header documents the new ownership split: the filename overload closes its
input, the `RandomAccessFile` overload does not. The declared failure mode
widens from "if reading fails" to "if opening, reading or closing fails", and a
failed close raises `Exception::InvalidValue`.

The motivation is a Windows file-lifetime race: Arrow's asynchronous read tasks
can retain shared ownership of the input after `ReadTable` returns, so the file
could not immediately be replaced. The chunk-combining logic moved into the
`RandomAccessFile` overload unchanged in behaviour — it still avoids
`CombineChunks()` when every column already holds a single chunk.

**No ported surface is affected.** Every Parquet and Arrow header is `unmapped`
in the coverage ledger — 21 of them, including `ParquetFile.h` — and the crate
has no Arrow or Parquet code. `docs/REPOSITORY_ANALYSIS.md` keeps columnar
formats outside the portable Rust core pending a separate dependency decision.

## Test and reference impact

Two class tests changed, `ParquetFile_test.cpp` and `StopWatch_test.cpp`, the
latter making its timing assertions portable. Neither is a pinned reference
path. `StopWatch.h` and `StopWatch.cpp` — which *are* pinned, by
`tests/data/progress_logger_provenance.json` — are unchanged, so that manifest's
hashes still hold at the new target.

No carried-forward reference path changed bytes, so all 220 keep their recorded
hashes and none needs a new behaviour review. The class-test corpus stays at 703
files; its line count moves from 229,849 to 229,869 through the two edits above.

## Scope

Registered public headers remain **786**, and the whole registration union is
untouched: no `sources.cmake`, `CMakeLists.txt`, `includes.cmake` or
`OpenSwathAlgoFiles.cmake` changed, verified by
`tools/core_sdk_retarget.py`, which refuses to run when a registration input
moves and re-hashes all 125 registration-evidence files to prove it. The
scientific inventory stays at 1,578 files.

Source inspection and hash checks do not establish a native C++ build or
differential numerical parity. Rust build and test results are recorded
separately in [Validation](VALIDATION.md).
