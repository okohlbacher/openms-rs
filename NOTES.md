# B7 lane notes (CLI-01..05)

Merged prerequisite worktree-wf_ee239125-5be-1 (fast-forward to 9fcb962).

## Scope (cli pin c19e494, `git -C /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages/cli show c19e494:<path>`)
(The sandbox refuses `for` loops / `cd` around git: run one plain `git -C <abs> show` per call.)
- TOPPExternalToolBase.h/.cpp; test tests/source/TOPPExternalToolBase_test.cpp (ls -ld . / ls -0 / missing exe)
- MapAlignerBase.h/.cpp; tests MapAlignerBase_test.cpp, MapAlignerDefaults_test.cpp
- SearchEngineBase.h/.cpp; test SearchEngineBase_test.cpp (fixtures tests/data/Single_MS2_profileMode.mzML, degenerate_cases/{empty.fasta,empty_spec.mzML})
- INIUpdater.h/.cpp; test INIUpdater_test.cpp
- TOPPBase_defs.h: Citation already in src/cli/defs.rs; B0 closed the accessor-misuse arms (ToolContext::escaped_refusal, test accessor_misuse_ends_as_the_source_catch_blocks_end_it in tests/topp_cli_lifecycle.rs). Ledger entry in docs/core-sdk-reviewed-apis.json is stale (still "partial").

## Findings
- TOPP tools THROW the three defs exceptions directly: IDSplitter.cpp:98, IsobaricWorkflow.cpp:508, LuciphorAdapter.cpp:447 (RequiredParameterNotGiven), MapAlignerIdentification.cpp:282,578 (WrongParameterType), ProteinQuantifier.cpp:1140. -> add ToolError::{required_parameter_not_given, unregistered_parameter, wrong_parameter_type} with the catch texts (TOPPBase.cpp:466-493; RequiredParameterNotGiven quotes what() unless it starts with ').
- TOPPExternalToolBase: ep.run(..., verbose=true, READ_WRITE, env). callbacks append + OPENMS_LOG_INFO<<chunk at debug>=4; on failure with debug<4 writeLogError_("Standard output: "+out) and ("Standard error: "+err). Rust: src/system/external_process.rs ExternalProcess::with_callbacks/run(&Invocation, verbose)->Result<RunReport>.
- ToolSpec lacks registerTOPPSubsection_ (map assign, overwrites) -> add ToolSpec::register_topp_subsection (needed by MapAlignerBase).
- MapAlignerBase::getModelDefaults = TransformationModelDefaults::getDefaults, being ported by lane worktree-wf_ee239125-5be-30 as analysis::transformations::params::transformation_model_defaults. Do NOT duplicate: named gap.
- MapAlignmentTransformer: src/analysis/alignment_transformer.rs transform_experiment.
- SearchEngineBase deps: format::mzml::centroid_info(_with_options), system::file find_database (FileContext method), analysis::peptide_indexing::PeptideIndexing::run, FileHandler::get_type.
- INIUpdater: std::map<TDI,TDI> ordered/looked-up by operator< == sort_key string equality. list::create("") -> empty. ToolHandler::get_topp_tool_list(&self).
- Oracle: build like ../oracle/toppbase-completion/ctd_driver.sh (add -lOpenMS_CLI), run on ibminode06 via ssh -J kohlbach@sshgw.cs.uni-tuebingen.de kohlbach@ibminode06.

## Order
1. defs ToolError ctors  2. external_tool.rs  3. ini_updater.rs  4. search_engine.rs  5. map_aligner.rs  6. oracle B7-cli  7. docs/ledger/report

## Progress
- DONE (commit 70710a3, clippy all-features green on dax): src/cli/{external_tool,ini_updater,map_aligner,search_engine}.rs, defs.rs ToolError ctors, ToolSpec::register_topp_subsection.
- Gaps: SearchEngineBase::reindex_ (native PeptideIndexing has no Param surface/ExitCodes); MapAlignerBase::getModelDefaults (lane 5be-30 transformation_model_defaults); transformSpectraFiles_ progress logging.
- Oracle ../oracle/B7-cli (driver.cpp, cases.tsv, run.sh) shipped to ibminode06:/scratch/kohlbach/oracle/B7-cli, started detached; poll results/DONE:
  ssh -J kohlbach@sshgw.cs.uni-tuebingen.de kohlbach@ibminode06 ls /scratch/kohlbach/oracle/B7-cli/results
  then rsync results back to ../oracle/B7-cli/results and copy retained cases into tests/data/topp_cli_base_classes/.

## Next step
1. Fetch oracle results; write tests/topp_cli_base_classes.rs replaying cases (template tests/topp_b0_release.rs) + INIUpdater/class-test transcriptions.
2. Docs: TOPP_CLI_SUPPORT.md "Tool base classes" section; provenance json; doc coverage; ledger fragment; report.
