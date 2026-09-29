# INTEG-B notes (integrator group B, wave 4)

Branch: worktree-wf_a10fdd71-2d2-11, base main f0990c06.

Lanes to merge (commits ahead of main):
- B0  worktree-wf_ee239125-5be-1   (112)
- B6  worktree-wf_cd066dab-0a7-1   (15)
- B7  worktree-wf_cd066dab-0a7-37  (124, contains B0 + -33)
- B8a worktree-wf_cd066dab-0a7-38  (127)
- B8b worktree-wf_cd066dab-0a7-39  (132)
- B8c worktree-wf_cd066dab-0a7-40  (124, contains B0 + -36)

Constraints: no local cargo (dax gate slot integ-B, sweep slot integ-B-sweep);
git commands must be plain single commands (worktree isolation guard).

## Done
- All six merged, no textual conflicts (a29ee5f9 B0, 95ddd81d B7, 4f2bc804 B8a,
  00c5683a B8b, 86c94e59 B8c, 9747db6a B6).
- Logs: scratchpad/integB/ (check1.log = cargo check --all-features --all-targets on dax).

## Complete claims to attack
B0: MSExperiment.h, LogStream.h. B6: NLargest, ThresholdMower, Normalizer, RankScaler,
SqrtScaler, WindowMower, SNEstimatorMedian, TheoreticalSpectrumGenerator, SpectrumAnnotator,
NucleicAcidSpectrumGenerator, DefaultParamHandler. B7: TOPPExternalToolBase, INIUpdater,
TOPPBase_defs. B8a: FeatureXMLFile. B8b: MathFunctions, PeptideIdentification, USI,
DecoyGenerator. B8c: Param.h.
Overlap: TOPPBase.h + ParamJSONFile.h in B0 and B8c fragments (B8c is later).

## Findings/fixes so far
- 36a9c871 collapse B8c list_item_text onto concept::numeric_formatting::ostream_g.
- 1e1124ac CLI to_int32/to_double -> string_utils (1e-400 now refused like Release).
- ca214457 ParamXML ZIP test + oracle ../oracle/integ-B-paramxml-zip (ParamXMLFile.h last gap closed -> complete).
- Complete claims checked OK: TOPPExternalToolBase, INIUpdater, TOPPBase_defs, FeatureXMLFile,
  MathFunctions, PeptideIdentification, USI, TSG (B6 chem trio), NLargest/ThresholdMower/Normalizer.
- Ledger: TOPPBase.h + ParamJSONFile.h take B8c's entry (superset of B0's).
- ZIP (B8a) also touches ledger scopes of ConsensusXMLFile, XMLFile, IdXMLFile, FileHandler,
  *Handler.h: amend text, status unchanged.
- Proposed CI: scratchpad/integB/jobs_patch.py, new_lines.tsv, run_lines.sh (slot integ-B).
- DefaultParamHandler unknown-key warning wording differs from Release (lead decision).

- 5c39075b InfoStream (B7) now wraps concept::log_stream::LogStream (dup cache removed).
- c-manifest: tests/data/topp_cli_base_classes_provenance.json (B7 had none).
- Ledger written: .planning/ledger-fragments/wave4/INTEG-B.json (56 entries) by
  scratchpad/integB/make_ledger.py; validate_ledger.py runs core_sdk_coverage.build (ok, complete 132->148).
- new_lines1 (merged tree): 13 proposed slice lines all PASS, counts in new_lines1.out.
- sweep1 slot integ-B-sweep (merged tree before fixes); sweep2 slot integ-B-sweep2 (final tree);
  new_lines2 at final tree (slot integ-B).

- e76f1ce1 zip_input multi-entry ParseError message = full what() (" in: " tail).
- Report draft: scratchpad/integB/INTEG-B.md (fill SWEEP_RESULT), copy to .planning/lane-reports/INTEG-B.md.
- sweep1 (merged, pre-fix): test-job lines 1-23 all pass (5901 all-features); stopped at 24 to free dax.
- sweep3 = final sweep, slot integ-B-sweep2, log dir scratchpad/integB/sweep3.

- sweep3 (slot integ-B-sweep2, tree = e76f1ce1): the ssh died with the previous incarnation
  after line 46. Lines 1-46 ALL PASS (remote logs /scratch/kohlbach/agents/integ-B-sweep2/.ci-line.N.log,
  summarised by scratchpad/integB/summ.sh). Lines 47-54 relaunched detached on dax via
  .rest.sh (DONE marker .rest.DONE, output .rest.out in the slot).

- sweep3 COMPLETE: 54/54 pass at e76f1ce1 (47-54 via .rest.sh: 47:11, 48:15, 49:12, 50:1 passed; fmt/fma-check/clippy/doc rc=0).
- Digestion attack (B8b finding): C++ tokenize_ yields a duplicate 0 for a leading zero-width cut
  (Asp-N/Lys-N on D/K-leading seq) -> ProteaseDigestion::digest (complete on main) emits an extra
  discarded empty product and shifted MC products; Rust dedups. Oracle launched:
  ../oracle/integ-B-digest-leading (run.sh pull). Plan: downgrade ProteaseDigestion.h in INTEG-B.json
  with named gap unless fixed; EnzymaticDigestion.h (partial) scope names countInternalCleavageSites.
- Still to member-check: DefaultParamHandler, WindowMower, RankScaler, SqrtScaler, SNEstimatorMedian,
  SpectrumAnnotator, NASG, DecoyGenerator (B8b, missing from report table), MSExperiment, LogStream.

- a54ea6a8 digestion leading-site oracle + tests/digestion_leading_site.rs (7 pass stable/1.85, clippy ok).
  Result: recorded native difference, ProteaseDigestion.h stays complete; INTEG-B.json amended by
  scratchpad/integB/ledger_digest.py (58 entries; do NOT rerun make_ledger.py, it would drop them).
- MSExperiment.h / LogStream.h name-level member check vs their docs: hold.

- facf18ef JSON -0 scanner collapse (param_json tests pass on dax, clippy ok).
- Report draft updated (items 8, 9; duplicates; table). sweep4 was stopped (pre-refactor tree).

## Next step
1. Final sweep sweep5 (slot integ-B-sweep2, log scratchpad/integB/sweep5) on tree facf18ef.
   If the ssh dies, remote logs stay in /scratch/kohlbach/agents/integ-B-sweep2/.ci-line.N.log
   (summarise with scratchpad/integB/summ.sh; resume with .rest.sh FROM TO).
2. Fill SWEEP_RESULT in scratchpad/integB/INTEG-B.md; cp to .planning/lane-reports/INTEG-B.md.
3. Delete NOTES.md in final commit.
