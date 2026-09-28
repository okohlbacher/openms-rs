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

## Incarnation 3 (restarted; same worktree)
Found uncommitted ToStr for char/bool + test (committed as WIP). clippy2.log (scratchpad A5a/logs) FAILED on tests:
list_utils_io.rs:81,92,138,144,161 (excessive_precision, inconsistent_digit_grouping, approx_constant),
concept_basics.rs:112-113 (assertions_on_constants), concept_macros.rs:76 (const_is_empty).
Scratch: $S=/private/tmp/claude-501/-Users-kohlbach-Claude-OpenMS-OpenMS4-Exploration-OpenMS4-R/c1467273-2c73-458b-a0aa-f084054684de/scratchpad/A5a (cpp/i/t symlinks to C++ pin, logs/).

## Incarnation 3 plan (cheap closes, in order; commit each)
- EmpiricalFormula: from_element, hydrogen, water, is_charged, has_element, to_map, elements() iter, lightest_isotope_mass,
  theoretical_isotopes_number, add_charge_adduct, to_source_string(toString), source_stream_text(operator<<: omits count<=1,
  charge -2 prints "--2" = C++ defect candidate), Ord via graph_cmp, derive Hash, source_hash64. Tests: EmpiricalFormula_test.cpp
  lines 115-122 (<), 203-207, 243-290 (addChargeAdduct), 375-399, 639-691.
- Ribonucleotide: Display ("Ribonucleotide '<code>' (<name>, <formula operator<<>)"), Hash, source_hash64 (Ribonucleotide.h:191-215).
- IsotopeDistribution: merge (cpp: sortByMass, trimLeft, trimRight, bins; IllegalArgument if more points), Ord (size then (mz,int)), Hash, Sorted enum.
- Element/ElementDB: by atomic number, has_element, name/symbol/number maps; Element Ord/Display/Hash. C++ ElementDB also holds
  isotope elements "(13)C" (name "(13)Carbon") + D/T symbol aliases (ElementDB.cpp:570-660) -> Rust gap unless added.
- EnzymaticDigestion filterByMissedCleavages + Specificity values; DigestionEnzymeDB hasEnzyme(ptr); Coarse generator setters.
- Then: support docs, provenance, ledger fragment (all 36), report.

## Next step
1. Add ToStr for char/bool (ListUtilsIO completeness). Check clippy2.
2. Re-read remaining scopes vs C++/Rust: DataValue, ConstRefVector, FASTAContainer, Exception, all CHEMISTRY, XLMS -> record remainders.
3. Support docs (10) + provenance manifests (script) + ledger fragment + report.

## Incarnation 4 (new worktree wf_cd066dab-0a7-3, merged worktree-wf_ee239125-5be-8 = incarnations 1-3 + lead's save 920f8a5)
Lead's save 920f8a5 holds uncommitted EmpiricalFormula work (+432 lines in src/chemistry/empirical_formula.rs) — ungated.
Old scratch dir is gone; new scratch: $S=<scratchpad>/A5a (recreate cpp symlinks).
Next step:
1. Review 920f8a5 EmpiricalFormula diff; gate check (lib + tests) on dax; fix clippy leftovers noted in incarnation 3.
2. Continue incarnation 3 plan from Ribonucleotide onward; then docs/ledger/report.
