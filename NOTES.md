# O5b lane notes (Library and inference)

Merged prerequisite worktree-wf_cd066dab-0a7-30 (O1a, fast-forward). No earlier incarnation notes.

## Scope (all missing in Rust; coverage entries at docs/core-sdk-coverage.json ~1816, 2002-2145)
C++ at core bc9cc12, src/openms/{include,source}/.../ANALYSIS/OPENSWATH/:
LevelContextInference, OpenSwathGeneInference, OpenSwathInferenceConfig, OpenSwathInferenceData,
OpenSwathLibraryIDNormalizer, OpenSwathPeptideInference, OpenSwathPeptidoformInference,
OpenSwathPercolatorScoring, OpenSwathProteinInference.
Class tests exist only for: LevelContextInference, OpenSwathLibraryIDNormalizer, OpenSwathPercolatorScoring.
Rust target dir: src/analysis/openswath/ (append pub mod in src/analysis/openswath.rs).

## Next step
1. Read the headers/cpp; order by dependency (Config/Data first).
