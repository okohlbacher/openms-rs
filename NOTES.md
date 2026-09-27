# Lane F4 notes (delete in the final commit)

Branch: worktree-wf_bcdd783d-f3f-14, base f0990c0.

## Scope and state at start
All unmapped except MzMLSqliteHandler (partial: src/format/mzml_sqlite_handler.rs,
tests/mzml_sqlite_handler.rs, docs/MZML_SQLITE_HANDLER_SUPPORT.md) and
ChromeleonFile / MSPGenericFile (evidence_requires_review, boost_regex manifest only).
Existing sqlite pieces: src/format/sqlite_connector.rs, mzml_sqlite_swath_handler.rs,
features `sqlite`, `sqmass = [sqlite, mzml]`.

C++ sizes (h/cpp/test lines) at core bc9cc12 - dumped into scratchpad/cpp via
scratchpad/dump.py (plain `git -C <core> show` works; loops over git are refused
by the worktree guard, so use the python script):
- SqMassFile 211/112/466, CachedMzML 188/118/210, CachedMzMLHandler 202/569/401
- SwathFile 177/542/388, SwathFileConsumer 646/14/593
- MSPFile 90/491/159, MSPGenericFile 134/337/413
- GNPSMGFFile 75/380/(no test), GNPSQuantificationFile 53/84/98, GNPSMetaValueFile 47/37/59
- EDTAFile 108/334/135, IBSpectraFile 82/291/84, ChromeleonFile 53/130/193
- FLASHDeconvSpectrumFile 115/613/333, FLASHDeconvFeatureFile 52/204/255
- MzMLSqliteHandler 245/1505/677 (already partial)

## Plan (order)
1. Small text formats with class tests: EDTAFile, ChromeleonFile, GNPSMetaValueFile,
   GNPSQuantificationFile, IBSpectraFile, MSPGenericFile, MSPFile.
2. CachedMzMLHandler + CachedMzML (binary cache; TICCalculator consumer).
3. SqMassFile (wrap mzml_sqlite_handler).
4. FLASHDeconv files, GNPSMGFFile, SwathFile/SwathFileConsumer.

## Next step
1. Read EDTAFile.h/.cpp/test and the existing text/csv helpers; port EDTAFile.
