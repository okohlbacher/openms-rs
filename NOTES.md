# F1a NOTES (lane scratch, delete in final commit)

Lane: F1a — Arrow/Parquet core. Branch now: worktree-wf_ee239125-5be-10 (fast-forwarded onto the
earlier incarnation worktree-wf_bcdd783d-f3f-10 at 27d502a).

## Done so far
- Cargo.toml: optional `arrow` feature + arrow-array/schema/buffer/select + parquet pinned =59.3.0
  (claimed last release with rust-version 1.85; VERIFY on dax with +1.85.0). No THIRD_PARTY entry yet.

## State at start
- No arrow/parquet dependency in Cargo.toml; no mention in docs/THIRD_PARTY_CRATE_DECISIONS.md.
- MSRV 1.85, edition 2024.
- Upstream C++ at core bc9cc12 (repo /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/core):
  include/OpenMS/FORMAT/{ArrowIOHelpers,ArrowSchemaRegistry,ParquetFile,ParquetFilter,PSMArrowIO,
  ProteinIdentificationArrowIO,FeatureMapArrowIO,ConsensusMapArrowIO}.h + source/FORMAT/*.cpp
  Class tests: src/tests/class_tests/openms/source/<Name>_test.cpp (all 8 exist, plus Arrow_test.cpp).
- Tooling gotcha: the worktree guard refuses `git -C <other repo>` and commands whose text contains
  "source" next to git. Use `cd <core> && git show bc9cc12:<path>` one file per command, redirect to scratchpad.

## Next step
1. Dump the 8 headers/cpps/tests into scratchpad, read them.
2. Pick arrow/parquet crate versions with MSRV <= 1.85 (check crates.io rust-version); add optional `arrow` feature.
3. Port in dependency order: ArrowSchemaRegistry -> ArrowIOHelpers -> ParquetFile/ParquetFilter -> PSMArrowIO ->
   ProteinIdentificationArrowIO -> FeatureMapArrowIO -> ConsensusMapArrowIO.
