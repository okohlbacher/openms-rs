# Porting instructions

Whenever repository analysis, porting or testing identifies an issue in the
original C++ SDK, record it in [OpenMS_CPP_ISSUES.md](OpenMS_CPP_ISSUES.md).
Include the source revision, affected files and functions, a concrete trigger,
explanation of the incorrect behavior, proposed C++ fix, evidence and the Rust
port's handling. Keep stable issue IDs and update existing entries instead of
duplicating them.

Distinguish an executed C++ reproduction from source review, independently
derived expectations and Rust-only tests. Retain durable evidence where available.
Do not claim that a proposed fix has been applied upstream. Unconfirmed candidates
must be labeled as such; ordinary API differences are not automatically defects.

When multiple agents work in parallel, send findings to the integrating agent,
which owns this shared log, to avoid conflicting edits and duplicate IDs.
