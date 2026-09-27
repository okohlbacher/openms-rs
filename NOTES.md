# Lane A1 — EVID-01 notes (restart-safety file; deleted in the final commit)

## Task
Re-run the three mislabelled oracles on ibminode06 against the Linux Release
install `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`:
- `../oracle/topp-early-bundle` (run.sh -> cases/oracle.py; macOS product-sdk Debug core 4fdec46)
- `../oracle/topp-file-info-tool` (oracle.py)
- `../oracle/topp-peak-picker-tool` (oracle.py, derived/)
New oracle goes in `../oracle/A1-release-rerun/` (host+build label honest), run twice, DONE marker.
Then re-check every test that cites those oracles (grep the three names under tests/),
fix port (src/cli/tools/<tool>.rs) where Release differs; framework (src/cli.rs,
src/cli/context.rs) divergences -> report to lead (B8c owns).
Also fix tests/topp_peak_picker_hi_res.rs:~1179 which calls a macOS SDK oracle "Release".

## Known divergence
directory as `-in`: Linux Release exits 12, port exits 10.

## State
- Merged worktree-wf_ee239125-5be-1 (fast-forward to 9fcb962). No earlier incarnation.

## Next step
1. Read the three oracle.py drivers; port them into ../oracle/A1-release-rerun with a
   Linux runner (timeout instead of perl alarm; tool dir = release install bin).
2. rsync to ibminode06, run detached twice, pull results back.
3. Diff exit codes/stdout/stderr hashes against the old manifests; list divergences.
4. Grep tests for the three oracle names; re-check expectations; fix port.
