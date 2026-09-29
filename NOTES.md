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

## Next step
1. Wait for check1; fix semantic merge breaks.
2. Start sweep integ-B-sweep; attack complete claims meanwhile.
