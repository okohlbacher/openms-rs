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

## Next step
1. docs/REPRODUCIBILITY.md: measure sha256 assertions (all 64-hex in repo JSON/
   MD/TSV/RS), classify by what the digest names: in-repo file (recompute!),
   pinned upstream source (needs public OpenMS4-core at bc9cc12), ../oracle or
   /ceph or product-sdk (not reproducible). Read tools/check_core_sdk.py --source.
2. Lane report .planning/lane-reports/A4.md + ledger-fragments/wave4/A4.json.
3. Delete NOTES.md in the final commit.
