# Lane A5a NOTES (restart-safety scratch; deleted in final commit)

Brief placeholder `__A5A__` was not expanded. Authoritative list: `.planning/a5-lists.json` key `A5a`
(main checkout, not worktree): 36 headers, all `partial` except UnimodXMLDataProvider (`evidence_requires_review`).
The brief's examples (MultipleTesting, SysInfo, PythonInfo, RWrapper, JavaInfo, GumbelDistributionFitter,
USI, MetaInfoRegistry) are all in the A5b list -> NOT this lane.

A5a list: XLMS/OPXLDataStructs, XLMS/OPXLHelper, APPLICATIONS/ConsoleUtils, CHEMISTRY/{AASequence,
DigestionEnzymeDB, DigestionEnzymeProtein, Element, ElementDB, EmpiricalFormula, EnzymaticDigestion,
ISOTOPEDISTRIBUTION/{CoarseIsotopePatternGenerator, IsoSpecWrapper, IsotopeDistribution}, RNaseDB, Residue,
ResidueDB, ResidueModification, Ribonucleotide, UnimodXMLDataProvider}, CONCEPT/{Colorizer, CommonEnums,
Exception, HashUtils, Macros, PrecisionWrapper, Types, UniqueIdIndexer, VersionInfo}, DATASTRUCTURES/{ConstRefVector,
ConvexHull2D, DBoundingBox, DataValue, FASTAContainer, IsotopeCluster, ListUtilsIO, MatchedIterator}

Ledger dump: scratchpad a5a_ledger.txt (regenerate from docs/core-sdk-coverage.json `headers` list, key `header`).

## Status per header (fill in as re-read)
(none yet)

## Next step
1. Re-read each scope vs C++ header at pin bc9cc12 (core). Small CONCEPT/DATASTRUCTURES first (cheap closes).
2. Record true remainder in ledger fragment .planning/ledger-fragments/wave4/A5a.json (main checkout path).
3. Close cheap gaps; commit per header group.
