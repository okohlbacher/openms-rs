# Lane A4 notes (DEBT-05, EVID-03, EVID-04)

Branch worktree-wf_ee239125-5be-7 (fast-forwarded onto the earlier incarnation's
notes commit). Scratch scripts in the session scratchpad (parse_cpp.py,
class_*.tsv = the manual classification of all 354 entries, build_index.py,
splice.py, declare_retained.py = the SOURCE_PROVENANCE block for the lead).

## Done (commit 2c7278d-ish, "[DOC,TEST] Index every C++ issue entry ...")
- OpenMS_CPP_ISSUES.md: header corrected (001-060 at 82ce5b3, rest at bc9cc12
  + cli/topp/test-data; 38 name the 4fdec46 SDK); index of all 354 entries.
- check_source_citations.py: --verbose locations; UNCHECKABLE record enforced
  when every pin/retained file is reachable; tests added (69 pass).
- ../oracle/A4-retained-third-party/{eigen,boost,libstdcxx} fetched from
  ibminode06 Release install, sha256.txt, retain.sh, PROVENANCE.txt.
  With the lead's SOURCE_PROVENANCE block: 0 unreachable, 2 bare (Rust lines).

## Incarnation 3 findings (reproducibility)
- tools/measure_reproducibility.py (committed WIP 00f4474) runs in ~2 min:
  16860 digest occurrences / 11013 distinct; 75.6% / 70.1% third-party
  recomputable (in_repo match + upstream match + Cargo.lock).
- check_core_sdk.py sets: 4832 hash assertions (3056 distinct pairs): 758
  in-repo fixtures recomputed by CI; 310 ../oracle external shape-only;
  3757 C++ (current_sources 1813, graph 21, inventory 1578, carried 220,
  registration 125) shape-only in CI -- ALL match core bc9cc12 git objects
  (scratch A4/sdk_sets.py). Core is public: github.com/okohlbacher/OpenMS4-core.
- 428 in_repo 'mismatch' are mostly baseline_/previous_/source_/base_ keys
  (older versions of in-repo files): add a git-history lookup family.

## Incarnation 4 findings
- Merge prerequisite already an ancestor (1165047). Re-measured: same numbers.
- PUBLIC on GitHub (gh api): okohlbacher/openms-rs (this repo), OpenMS4-core
  @bc9cc12, OpenMS4-cli @c19e494, OpenMS4-topp @174b576, OpenMS4-test-data
  @0cb15f2 -- all four pins resolve. So the upstream family is third-party
  recomputable (not checked by CI; CI runs check_core_sdk.py without --source).
- No Rust code computes sha256 (no sha2 dep, no SHA IV constant): Rust digest
  literals are prose. CI recompute = Python gates (rust.yml:150-176,223) +
  cargo --locked crate checksums. CI gates read nothing outside the checkout,
  so running them in place in this worktree is CI-faithful.
- Plan: trace hashlib in each CI gate (scratch A4/trace_digests.py), intersect
  computed digests with recorded ones -> "recomputed by CI" figure.
- Plan: add git-history lookup for in_repo mismatches (older blob of same path).
- DONE 12e1f3b: measure tool has content index (tree/history/pinned) + --ci.
  Lab: 16860 occ / 10167 distinct; bound 82.7%/78.5%, content-only 4.0%,
  lab-only 8.5%, none 4.9%; CI recomputes 1902 occ (11.3%) / 1534 distinct.
- DONE 6a14273: citation checker reads core from packages/core git objects
  (third party had 158 false problems at f0990c0); dedupes identical files.
- THIRD PARTY on spock (/scratch/kohlbach/agents/A4-3p, scripts in scratch A4/
  third_party*.sh, runs run1/run2 at f0990c0, brun1/brun2 at 6a14273 bundle):
  public GitHub clones only; check_core_sdk.py --source PASSES (2107 files);
  measure tiers identical to lab; citations 4112 checked/47 skipped/17 unres
  (lab 4139/34/3; diff = retained libstdc++ under ../oracle). generate_enzymes
  needs rustfmt on PATH (passes with ~/.cargo/bin).
- WARNING: tools/check_schema_feature_graph.py runs cargo -- never run locally
  (I did once by accident, deleted the 694M target/).
- In progress: measure binding fixes (prefixed keys, log family, prose adopts
  manifest binding) -> rerun, then write docs/REPRODUCIBILITY.md.

## Next step
1. docs/REPRODUCIBILITY.md: measure sha256 assertions (all 64-hex in repo JSON/
   MD/TSV/RS), classify by what the digest names: in-repo file (recompute!),
   pinned upstream source (needs public OpenMS4-core at bc9cc12), ../oracle or
   /ceph or product-sdk (not reproducible). Read tools/check_core_sdk.py --source.
2. Lane report .planning/lane-reports/A4.md + ledger-fragments/wave4/A4.json.
3. Delete NOTES.md in the final commit.
