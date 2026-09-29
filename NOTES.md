# INTEG-A notes (integrator, group A)

Worktree branch: worktree-wf_a10fdd71-2d2-12, base f0990c0.

## Lane branches
- A1  worktree-wf_cd066dab-0a7-41 (2eb259fc)  EVID-01 oracle relabel
- A2a worktree-wf_ee239125-5be-3  (6dd92a78)  ledger only (net diff empty)
- A2b worktree-wf_ee239125-5be-4  (eade1684)  ledger only; report in its worktree
- A3a worktree-wf_cd066dab-0a7-2  (ca39a2c6)  DEBT-04 source order, code + tests
- A3b worktree-wf_ee239125-5be-6  (ea4b19c2)  docs/OPEN_ITEMS.md
- A4  worktree-wf_ee239125-5be-7  (60573e2e)  REPRODUCIBILITY, citation checker; report in its worktree
- A5a worktree-wf_cd066dab-0a7-3  (fad3f5fa)  chem/enzymes/unimod, feature unimod-xml
- A5b worktree-wf_cd066dab-0a7-4  (865a19f8)  IdentificationData etc.

## State
- ALL 8 lanes merged (HEAD 0da9ad40). A1 is built on B0 (worktree-wf_ee239125-5be-1 = 9fcb962a,
  Phase-B integration), so integ-A carries B0 too; integ-B (worktree-wf_a10fdd71-2d2-11) also merges B0.
- Conflicts: LOG_STREAM_SUPPORT.md (kept B0 text + A5a Colorizer pointer); 4 module roots keep-both.
- check_module_cycles: 4 new acyclic edges need --write (docs/module-cycles.json is not a lead file).
- A2a.md, A2b.md, A4.md reports are missing (worktree copies deleted); summaries are in the brief.
- Logs: scratchpad/logs/ (build1.log = test --no-run on dax slot integ-A)

## Findings
- build1 (test --no-run all-features) on dax: OK. sweep1 running (scratchpad/logs/sweep1).
- Member checker: scratchpad/members.py (headers exported to scratchpad/core-bc9cc12). STRICT=1 mode.
- Overlaps: ChromatogramTools (A2a complete vs A3a native_eq -> complete w/ A3a scope);
  StatisticFunctions (A3a complete vs A2b partial gap 2: Pearson/Matthews explicit f64::NAN) -> FIXING
  via x86_64 ops + collapse source_pearson (feature_finder_picked/scoring.rs) into pearson_correlation_coefficient;
  oracle ../oracle/integ-A-correlation (ibminode06) launched;
  MSSpectrum, IMDataConverter both partial -> merge scopes; SpectrumNativeIDParser A5b complete supersedes A3b.
- Ribonucleotide complete vs Element partial: same immutable-record design; Ribonucleotide has from_record (ok).

## Next step
1. Wait build1; fix compile errors.
2. Full sweep: OPENMS_GATE_HOST=dax OPENMS_GATE_JOBS=96 ~/.local/bin/openms-ci-sweep.sh integ-A <scratch>/sweep1
3. Attack complete claims (A5a, A5b, A2a, A2b, A3a StatisticFunctions).
4. Cross-group dup check vs integ-B (numeric_formatting, string_conversion, USI).
