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

## Status per header
Code done + gated (tests pass on dax): CommonEnums, Types, HashUtils, Macros, NumericFormatting(support), PrecisionWrapper,
ListUtilsIO, IsotopeCluster, DBoundingBox (data_structures/dbounding_box.rs), ConvexHull2D (map_points, source_hash64),
MatchedIterator (data_structures/matched_iterator.rs; comparison::matched_alignment now uses it), UniqueIdIndexer
(concept/unique_id_indexer.rs), Colorizer (concept/colorizer.rs), ConsoleUtils (cli::ConsoleUtils; Windows width = gap),
VersionInfo (concept/version_info.rs; VersionDetails moved from system::update_check, alias kept; 6 build statics = gap).
Planned statuses: complete for all above except ConsoleUtils, VersionInfo (partial).
Clippy fix pending (type_complexity in macros.rs) -> logs/clippy2.log.
docs/doc-coverage.json is lead's: report "--write needed". module-cycles: new acyclic edges analysis->data_structures,
data_structures->concept, kernel->data_structures, system->concept (report; lead --write).

## Next step
1. Add ToStr for char/bool (ListUtilsIO completeness). Check clippy2.
2. Re-read remaining scopes vs C++/Rust: DataValue, ConstRefVector, FASTAContainer, Exception, all CHEMISTRY, XLMS -> record remainders.
3. Support docs (10) + provenance manifests (script) + ledger fragment + report.
