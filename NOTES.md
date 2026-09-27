# Lane A5a NOTES (restart-safety scratch; deleted in final commit)

Incarnation 2 (brief fixed): the brief now names the 36 headers explicitly; they match `.planning/a5-lists.json` key `A5a`.
All `partial` except UnimodXMLDataProvider (`evidence_requires_review`). Prereq branch worktree-wf_bcdd783d-f3f-8 merged (fast-forward).
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
- WIP inherited (NOT yet wired into src/concept/mod.rs, not gated): src/concept/common_enums.rs, hash_utils.rs, types.rs.
  Next: review them against headers, wire `pub mod` lines, gate on dax.
- GumbelDistributionFitter / USI / MetaInfoRegistry / MultipleTesting etc. are A5b's list, not ours -> leave.

## Next step
1. Re-read each scope vs C++ header at pin bc9cc12 (core). Small CONCEPT/DATASTRUCTURES first (cheap closes).
2. Record true remainder in ledger fragment .planning/ledger-fragments/wave4/A5a.json (main checkout path).
3. Close cheap gaps; commit per header group.
