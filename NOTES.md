# INTEG-O notes (integrator for group O, branch worktree-wf_97065d91-48e-9)

Lane branches (merge order, all verified to exist):
OA 0a7-27 (wf_cd066dab), O1a 0a7-30, O1b 0a7-44, O1c 0a7-45, O2 0a7-47,
O3a wf_97065d91-48e-1, O3b wf_97065d91-48e-6 (also merges FF1 0a7-33),
O4a wf_cd066dab-0a7-49, O4b wf_a10fdd71-2d2-8, O5a wf_97065d91-48e-7,
O5b wf_cd066dab-0a7-46, O5c wf_97065d91-48e-5 (merges OA, O5b, F7b a10fdd71-2d2-5).

Known shared-file conflicts: src/analysis/openswath.rs (pub mod lines, keep both, cargo fmt),
src/format/mod.rs; O3a: register tests/data/openswath_peak_picking_provenance.json;
O3b: tests/data/mrm_feature_provenance.json under current_sdk_reference_manifests.
O3a has no lane report file: content is in the brief (4 of 5 headers complete,
MRMTransitionGroupPicker.h partial: fit_EMG pickApex, PickerPoolShared, lin_interpolate).

## Status
- [ ] 1 merge all lanes
- [ ] 2 cargo fmt + dedupe + sweep on dax (slot integ-O)
- [ ] 3 attack complete claims
- [ ] 4 attack evidence
- [ ] 5 INTEG-O.json + INTEG-O.md

## Next step
Merge lanes in order above; commit each merge; do NOT push, NOT rebase.
