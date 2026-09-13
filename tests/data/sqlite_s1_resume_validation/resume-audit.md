# Resume reconciliation at 45c424a

Read-only audit on 2026-09-13. No repository edits, test execution, remote actions, Claude invocation, or source export. Checked Git history, retained worktree status, central/leaf docs, manifests, and bounded existing /tmp review-log names.

## Authoritative checkpoint

- Local main and origin/main were both `45c424a97b504295645e33b77b3e7dcec74461e5`; worktree initially clean.
- All **42** retained `.claude/worktrees` are clean (including untracked-file status) and their HEADs are ancestors of current main. Nothing discovered needs cherry-picking. Preserve directories until cleanup is separately wanted.
- `068a1f5` preserved interrupted S1 handler/SWATH code; merge `35f3fd2` already put it on main. Do not redo or overwrite these ports from the older c17e3a2 summary.
- `525e818` introduced default `parallel`/Rayon and the thread-count determinism contract.
- `3b8ca94` integrated MATH/SYSTEM/COMPARISON leaf wave A and audit fixes.
- `4e8934b` integrated dependent wave B and audit fixes (including eliminating duplicate comparison implementations, restoring source f32 arithmetic, and withdrawing unsupported SIMDe closure).
- `d396f42` replaced hand-written FFT with rustfft 6.4.1 scalar planner. `45c424a` records crate-first policy. Preserve both decisions.
- Core source pin remains `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

`docs/CORE_SDK_COMPLETION.md` and its JSON inputs give the current inventory: **786 registered headers; 53 complete, 89 native_equivalent, 48 partial, 164 evidence_requires_review, 432 unmapped**. Thus 142 are closed under that ledger; not a verified SDK completion percentage. Five of 146 TOPP workflows remain validated. Wave B commit records COMPARISON 13/13, MATH 15/17, SYSTEM 12/14, but per-header JSON is the authority.

## Validation: distinguish existing claims from retained execution

Git commit messages contain these historical Linux-node claims:

| Commit | all features | no default | other stated checks |
| --- | ---: | ---: | --- |
| 068a1f5 S1 preservation | 3,899 | 2,640 | clippy, MSRV1.85, rustdoc, fmt, six Python gates |
| 3b8ca94 wave A | 4,201 | 2,942 | clippy, MSRV1.85, rustdoc, fmt, seven Python gates |
| 4e8934b wave B | 4,417 | 3,101 | network feature3,158; quality/MSRV/seven Python gates |
| d396f42 rustfft | 4,418 | 3,102 | clippy, MSRV1.85, rustdoc, fmt, doc coverage, cycles, feature graph |

These are **historical author-reported outcomes**, not a newly executed audit. The commits do not add a frozen full-run summary with runtime hashes and retained remote log paths analogous to `tests/data/sqlite_connector_validation/summary.json`. I did not find newer retained full-run logs in the changed tracked paths or bounded relevant /tmp names. Parent should run its fresh remote baseline and retain commands, outcomes, snapshot hashes and logs instead of promoting commit prose into newly verified evidence.

Existing durable S0 and FORMAT validation stays valid for its own historical hashes. Do not rewrite it to describe later code.

S1 SWATH has stronger retained leaf evidence: `tests/data/mzml_sqlite_swath_provenance.json::native_validation` records 18 tests per current/MSRV sqlite-only selection on kim, strict clippy/rustdoc, source/test SHA256, and four local logs under `tests/data/sqlite_s1_validation/`. The exact-pinned adapted C++ probe reproduces CPP-196/197/202; it is not a full SDK or main-handler oracle.

Main-handler provenance lacks a final native_validation block. `docs/MZML_SQLITE_HANDLER_SUPPORT.md` explicitly remains pending integration, final hashes/review, and an activation-metadata fix. `handler_peer_review.json` retains two initial source-trace findings (expanded lossless size; sampled-noise snapshot self-read); it does not itself document dispositions. Parent/storage worker should connect fixed tests to each finding.

## Independent-review caveats

SWATH manifest `independent_review` says `claude-fable-5-1`, `blocked_not_executed`: remote OAuth expired; local escalation rejected pending concrete packet/destination approval. Support doc repeats it. This is not review approval; do not relabel it after native tests pass.

The Wave A commit says five numerical packages were reviewed by Codex Astra, no blockers. However the available `/tmp/cx-comparison-base.log`, `cx-fitters.log`, `cx-scorers-core.log`, `cx-scorers-advanced.log`, and `cx-fft-stats.log` all end in usage-limit errors, while corresponding `.md` files are empty. These attempts are not successful review reports. There might have been other execution not retained in the inspected scope; do not claim it either way without evidence. Claude audit-fix commits themselves are real preserved work and can be discussed as reviewed changes without asserting a missing external model approval.

No Claude invocation or source packet export was attempted during this audit.

## Concrete central-document updates

1. `docs/PORTING_WAVES.md`: add completed integration/preservation checkpoints for S1 and waves A/B/rustfft; make clear S1 code is merged but its contract/evidence closeout remains open. Keep OSW/OMS/S2 pending. Do not call all of S1 complete just because merged.
2. `README.md`: replace connector-only storage paragraph with optional sqlite connector + SWATH lookups and sqmass handler availability, linked to honest support limits. State SqMassFile/consumer/access/OSW/OMS still outstanding. Add one compact paragraph or table entries linking new numerical, comparison, process/network APIs; mention optional network and default parallel rather than listing dozens of headers.
3. `docs/PORTING_STATUS.md`: replace 'SQLite-family formats remain the next planned wave' with merged S1 implementation awaiting closeout and pending adapters; add wave A/B checkpoint links.
4. `docs/VALIDATION.md`: insert a **new** fresh remote resume checkpoint when parent baseline completes. Until then say later commit messages report counts but no new full-run artifact was found. Keep old S0 'next stages remain unported' statement explicitly historical, or replace it with a dated crosslink so readers do not infer current absence of S1.
5. `docs/SQLITE_STORAGE_PLAN.md`: its opening 'No Rust implementation...' describes plan creation, but should be labeled historical preparation and crosslinked to current S1 state.
6. `docs/FFT_SUPPORT.md`: implementation changed to rustfft, but API table still says recursion is flattened into equivalent stage loop; that describes removed hand implementation. Also scalar planner avoids CPU-dispatched kernels, but 'same bits on every machine' is a stronger claim than the retained same-node scalar-vs-SIMD measurement establishes. Record it as intended determinism contract and state actual tested platforms, unless cross-machine bitwise evidence is available. No numerical bug alleged here.
7. Shared C++ issue log ends at CPP-218 and did not change across Claude waves A/B. Several new manifests contain source defect lists, e.g. distribution fitters and network. Integrator should reconcile those candidates into stable IDs after source/evidence review rather than treating prose `cpp_issues` arrays as already centrally logged. Ordinary divergences must stay separate.

## Next-wave rationale

First finish **S1 contract/evidence closeout** against the newly merged numerical/system baseline. The activation metadata gap can silently lose source-supported RUN_EXTRA metadata, so it is a prerequisite for trustworthy sqMass storage. Retain source/test mappings, limits/rollback regressions and explicit review disposition before raising coverage.

Then follow already prepared `docs/SQLITE_CONSUMER_PLAN.md` for **S2**: SqMassFile, MSDataSqlConsumer, SpectrumAccessSqMass and the minimal OpenSwath ISpectrumAccess/BinaryDataArray/Spectrum/Chromatogram/metadata prerequisites. This exploits the existing handler and closes useful actual storage workflows. The plan accounts for six SqMassFile sections (87 in-section macros+8 helpers), seven access sections (39 macros), and no dedicated consumer test. Its CPP-204..217 findings are preparation, not runtime closure.

S2 needs view order/duplicate reconstruction over SQL IDs, explicit finish/error handling, settings preservation, run switching, bounded batches/atomic writes and Parquet declarations accounted for. OSW/OMS can proceed as independent mapped storage groups only with their real dependencies. FileHandler dispatch and CLI/TOPP porting should build on these concrete contracts, not assume that names or merged source establish feature completeness.
