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
Done in code (gate g2 pending, then docs/ledger):
- CommonEnums, Types, HashUtils -> src/concept/{common_enums,types,hash_utils}.rs, tests/concept_basics.rs (all members)
- Macros -> src/concept/macros.rs (thread_critical, openms_precondition!/postcondition!), tests/concept_macros.rs
- NumericFormatting (support) -> src/concept/numeric_formatting.rs; text_format.rs now delegates (behaviour identical)
- PrecisionWrapper -> src/concept/precision_wrapper.rs; ListUtilsIO -> text_format.rs (VecLowPrecision, StringListPush,
  to_str_low_precision*), tests/list_utils_io.rs
- IsotopeCluster -> src/data_structures/isotope_cluster.rs (moved ChargedIndexSet etc. from ffp defs.rs, re-exported)
- DBoundingBox -> src/data_structures/dbounding_box.rs generic, tests/dbounding_box.rs
New acyclic module edges: data_structures->concept (report; lead may --write docs/module-cycles.json)
- GumbelDistributionFitter / USI / MetaInfoRegistry / MultipleTesting etc. are A5b's list, not ours -> leave.

## Next step
1. Poll logs/g2.log (scratch A5a/logs). Fix failures. Commit milestone.
2. Next headers: Ribonucleotide (<<, hash; needs EmpiricalFormula <<+hash, Element hash, IsotopeDistribution hash),
   ConvexHull2D (getMapPoints via DBoundingBox1, hash), MatchedIterator, UniqueIdIndexer, Colorizer, ConsoleUtils API,
   VersionInfo, DataValue gaps, then CHEMISTRY re-reads.
3. Support docs: docs/CONCEPT_BASICS_SUPPORT.md, DBOUNDING_BOX_SUPPORT.md, ISOTOPE_CLUSTER_SUPPORT.md (referenced in rustdoc!).
4. Ledger fragment + report at end.
