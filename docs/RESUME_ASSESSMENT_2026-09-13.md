# Resume assessment and proposed next waves — 2026-09-13

This assessment checks local Git/worktrees, the current source ledger and code,
retained review evidence, and the already-running kim validation. It does not
start another conversion wave or claim to refresh GitHub. `origin/main` below
means the local remote-tracking ref. No implementation changes were made for
this assessment.

## Current checkpoint

- HEAD and local main: `45c424a97b504295645e33b77b3e7dcec74461e5`.
  Local `origin/main` points to the same commit.
- Active branch: `codex/sqlite-s1-resume`, with 28 modified tracked files and
  20 untracked files before adding this report. Those files are not committed
  or pushed by the resumed work.
- All 42 retained Claude worktrees are clean and their heads are ancestors of
  HEAD. There are no unmerged worktree results to recover.
- Core remains pinned to `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.
  CLI and TOPP have independent pins in `tests/data/topp_cli_provenance.json`:
  CLI `c19e49414bcd9ebdea42f89b3f74d2823205892c`, TOPP
  `174b576e244e100f2345ca57a8e79aaa607156df`.

Claude preserved the interrupted SQLite implementation in `068a1f5`, integrated
MATH/SYSTEM/COMPARISON waves A and B in `3b8ca94` and `4e8934b`, replaced the
hand-written FFT with rustfft's scalar planner in `d396f42`, and recorded the
crate-first policy in `45c424a`. Default Rayon parallelism is intentional.
These integrations should be kept; restarting or cherry-picking the old
worktrees would duplicate completed work.

The uncommitted resume adds SQLite limits/rollback and schema-guard corrections,
invalid Numpress-coordinate rejection, activation/intensity-unit metadata
transport, focused regressions and evidence reconciliation. New stable C++ issue
entries CPP-219–229 distinguish source review, native tests and unconfirmed
candidates. The main sqMass handler remains partial because its full metadata
snapshot inherits the native mzML subset.

## What the inventory actually says

| State | Committed baseline | Current working ledger |
|---|---:|---:|
| Complete | 53 | 53 |
| Native equivalent | 89 | 90 |
| Partial | 48 | 49 |
| Evidence requires review | 164 | 162 |
| Unmapped | 432 | 432 |
| Total registered core headers | 786 | 786 |

The working ledger therefore closes 143 headers, compared with 142 committed.
These are reviewed API states, not a scientific-parity completion percentage.
The working change closes the narrow SWATH handler and records the main SQLite
handler as partial. It is still subject to checkpoint integration review.

Five of 146 TOPP workflows have retained-output comparisons: DTAExtractor,
BaselineFilter, MapNormalizer, SpectraFilterWindowMower and MzMLSplitter.
A native TOPPBase-style framework exists, but OpenMS4-cli is not complete.
There is no comparable whole-package CLI completion ledger yet.

The current per-header ledger, rather than commit prose, is authoritative:
COMPARISON has 13/13 closed, MATH 13/17, SYSTEM 7/14. Some wave commit summaries
claim higher domain totals. Large remaining bodies include ANALYSIS (174
unmapped), FORMAT (104), FEATUREFINDER (43), ML (25), and QC (18).

## Validation and unresolved risks

Fresh local provenance, coverage, documentation-floor and diff-whitespace checks
pass. Public rustdoc coverage is 3,779/5,290 items (71.4%); a successful rustdoc
build does not mean all public APIs are documented or that their claims are true.

The existing remote run is under
`/ceph/ibmi/abi/oliver/openms-rs/results/sqlite-s1-resume-20260913-121422`.
It uses kim, node-local scratch, 32 compilation jobs and a separate target directory.
At inspection, these checks had completed successfully:

- Current Rust all features/all targets: 4,434 tests, no ignored tests.
- Current Rust no-default cargo test: 3,144 passing tests including doctests.
- All-feature doctests: 57 passing.
- Formatting, strict full Clippy and strict rustdoc.

The full Rust 1.85 all-target run subsequently failed (exit 101) in
`tests/system_process.rs`: 30 tests passed and four failed. Failures are
`constructed_callbacks_receive_the_two_streams_separately`,
`what_the_callbacks_collected_outlives_the_process_value`,
`a_budget_ends_the_call_when_a_descendant_still_holds_the_pipes`, and
`python_version_is_read_from_the_stand_in_and_empty_when_absent`.
Observed symptoms are failed process starts, missing captured output and an
empty interpreter version. Root cause is not established: do not call this a
compiler incompatibility or harmless flake yet. The suite's assertions do not
print the launch error in these cases. Retain the failing log (SHA-256
`883d9718fcbf52750463fd8202291a5ca475bf58d5a967462395453e352d60ea`),
capture the actual OS error, and compare isolated/parallel runs before changing
behavior. Full MSRV doctests passed; other matrix selections were still running.

**Update, same day (Claude, resuming this plan).** Root cause established, and it
is neither a compiler incompatibility nor Rust 1.85 specific. A throwaway copy on
kim with the `spawn()` error printed showed every failure to be `ETXTBSY` (errno
26, "Text file busy") on exec of a stand-in script written by the test moments
earlier: a child forked by a concurrently running test inherits that script's
write handle until its own exec. 18 of 25 parallel runs failed that way under
Rust 1.85 *and* 18 of 25 under 1.96; single-threaded runs failed 0 of 20. The
race is in the test harness, not in `ExternalProcess`. The tests in
`tests/system_process.rs` now run serially behind one lock; after the change 0 of
25 parallel runs failed under each toolchain. The earlier green 1.96 run had
passed by chance.
Focused leaf evidence records 31 handler tests and 10 activation tests passing
on Rust 1.85. This does not make the combined checkpoint green. The remote snapshot predates final
local documentation/provenance/ledger edits; compare runtime hashes and rerun
metadata gates on the final snapshot before publishing a checkpoint.

Most scientific groups still depend on source literals and native invariants,
not execution against a linked C++ SDK. Historical commit test counts and review
claims should remain historical unless matching retained logs/reports can be
located. Some available numerical review attempts ended at usage limits and do
not establish a completed independent review.

The requested S1 Fable review has not executed. Remote OAuth failed; local packet
export was rejected by automatic approval review pending explicit authorization
of the payload and destination. The renewed authorization request remains
unanswered. Codex review and passing tests must not be labeled Fable approval.

**Update, same day.** Claude ran the requested review after taking over the plan:
four Claude Fable 5.1 reviewers, one per area (SQLite handlers, mzML metadata
transport, the C++ issue log, ledger/documentation/CI), with every finding
re-checked by an independent agent told to refute it. Eight findings were
confirmed (seven distinct), none uncertain, one refuted. The two major ones were
a schema-invalid element order in metadata-only `<activation>` blocks, now fixed
with an XSD-validated regression, and the stale SWATH execution record, which the
final-snapshot validation run replaces. The record and every disposition are in
`tests/data/sqlite_s1_resume_validation/` (`fable-review.json`,
`review-dispositions.json`).

Concrete issues to address during stabilization:

- `docs/TOPP_CLI_SUPPORT.md` still describes a serial core and an unimplemented
  UpdateCheck, although the system implementation and default Rayon now exist.
  `src/cli.rs` stores `-threads` in context; it does not establish the promised
  execution policy. Logging and INI instance handling also remain incomplete.
- Network's exists-then-create sequence retains a concurrent overwrite race
  (CPP-227); some API/ledger text still promises no overwrite unconditionally.
- The Gumbel squared-negative-log-likelihood concern (CPP-228) needs an actual
  counterexample/oracle investigation. Preserve its candidate status; matching
  source behavior and proving statistical correctness are different checks.
- The C++ reference's documented TOPP launch failure is a historical observation.
  Recheck current binaries, package identity and run paths before treating it as
  a live blocker. Existing retained TOPP outputs remain usable independently.

## Revised execution plan

The user subsequently prioritized early FileInfo and FeatureFinderCentroided
builds. The [early TOPP build plan](EARLY_TOPP_BUILD_PLAN.md) replaces the previous
storage-first order. The status and evidence observations above are historical
observations from this assessment, not newly completed implementation.

1. Preserve the recovered checkpoint (the process-test failures are diagnosed and
   fixed; see the update above) alongside focused tool work; do not wait for full
   SQLite or SDK closure.
2. Build the reusable FileInfo reporter and CLI on existing mzML/featureXML
   readers; validate text/TSV, metadata, processing and statistics on source cases.
3. In parallel, port the picked-feature helper/trace-fitter/algorithm chain and
   connect FeatureFinderCentroided, starting with a tested non-FAIMS mzML preview
   and then closing seeds, fitting modes, FAIMS and output metadata contracts.
4. Add PeakPickerHiRes using the existing picker, then validate the complete
   inspection -> centroiding -> feature finding -> inspection workflow. Keep
   the five existing tools green; target an initial eight-executable bundle.
5. Resume remaining FileInfo formats, SQLite S2/S3, full CLI and other domain/tool
   waves after the early bundle, pulling forward only genuine target dependencies.

The linked plan records exact source paths, missing dependency groups, source
workflow cases, parallel ownership and acceptance criteria. No implementation,
commit or push is claimed by this planning revision.
