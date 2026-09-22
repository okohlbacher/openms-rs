# Lane NOTES (TOPPBase finish) -- delete before the final commit

Branch: worktree wf_673154fa-1e1-1 off main a660a0c.
Ledger: headers/packages/cli/include/OpenMS/APPLICATIONS/TOPPBase.h = partial (see docs/core-sdk-reviewed-apis.json scope).
Rust: src/cli.rs, src/cli/{context,defs,logging,processing,spec,usage,param_ctd,tool_handler,console,outcome}.rs, src/cli/tools/**
serde_json =1.0.150 float_roundtrip optional; features rna-json, proforma-json use it. paramxml = quick-xml + file-compression.

Gaps (ranked):
(a) JSON INI ParamJSONFile::load -> wire serde_json into paramxml feature; add THIRD_PARTY_CRATE_DECISIONS entry.
(b) InvalidValue exit 6->8; first per-tool ILLEGAL_PARAMETERS returns (MzMLSplitter, MapNormalizer, FileInfo, SpectraFilterWindowMower, lifecycle).
(c) FeatureFinderCentroided CITATIONS.
(d) tool bodies' error/warn lines -> cli::log_error/log_warning (file_info.rs:185, peak_picker_hi_res, feature_finder_centroided).
(e) UnregisteredParameter/WrongParameterType exit 12.
(f) Windows COLUMNS -- record only.

## Findings
(none yet)

## Next step
1. Find ParamJSONFile at core pin bc9cc12 and where TOPPBase calls it (cli pin c19e494); find Rust INI load path.
