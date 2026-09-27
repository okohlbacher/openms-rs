# B7 lane notes (CLI-01..05)

Merged prerequisite worktree-wf_ee239125-5be-1 (fast-forward to 9fcb962).

## Scope (cli pin c19e494, `git -C ../OpenMS4-tests/packages/cli show c19e494:<path>`)
- include/OpenMS/APPLICATIONS/TOPPExternalToolBase.h + source/..cpp; test tests/source/TOPPExternalToolBase_test.cpp
- MapAlignerBase.h/.cpp; tests MapAlignerBase_test.cpp, MapAlignerDefaults_test.cpp
- SearchEngineBase.h/.cpp; test SearchEngineBase_test.cpp (fixtures tests/data/Single_MS2_profileMode.mzML, degenerate_cases/*)
- INIUpdater.h/.cpp; test INIUpdater_test.cpp (needs ToolHandler = src/cli/tool_handler.rs)
- TOPPBase_defs.h (Citation + 3 exceptions); ExitCode already in src/cli/parameter.rs
- Plan doc: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-R/.planning/CLI-PACKAGE.md sec 2,3

## Order
1. TOPPBase_defs (Citation, exceptions)  2. TOPPExternalToolBase  3. SearchEngineBase  4. MapAlignerBase  5. INIUpdater

## Next step
1. Read the C++ headers/tests above and src/cli/{defs,tool_handler,context,spec}.rs, src/system/external_process.rs.
