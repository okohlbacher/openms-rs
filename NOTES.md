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

## Next step
1. Poll sweep1 (slot integ-B-sweep) and new_lines1.out; fix failures.
2. Re-sweep after fixes; write INTEG-B.json + INTEG-B.md.
