# IDMapper, IonIdentityMolecularNetworking and SiriusFragmentAnnotation

Ports of three OpenMS 4 core headers (core `bc9cc12`):

| Header | Rust | Tests | Provenance |
|---|---|---|---|
| `ANALYSIS/ID/IDMapper.h` | `src/analysis/id/id_mapper.rs` | `tests/id_mapper.rs` | `tests/data/id_mapper_provenance.json` |
| `ANALYSIS/ID/IonIdentityMolecularNetworking.h` | `src/analysis/id/ion_identity_molecular_networking.rs` | `tests/ion_identity_molecular_networking.rs` | `tests/data/ion_identity_molecular_networking_provenance.json` |
| `FORMAT/DATAACCESS/SiriusFragmentAnnotation.h` | `src/format/sirius_fragment_annotation.rs` | `tests/sirius_fragment_annotation.rs` | `tests/data/sirius_fragment_annotation_provenance.json` |

## IDMapper

### API mapping

| Source member | Rust |
|---|---|
| `enum Measure {MEASURE_PPM, MEASURE_DA}` | `Measure::{Ppm, Da}` |
| `IDMapper()` / `IDMapper(const IDMapper&)` / `operator=` | `IDMapper::default()`, `Clone` |
| `DefaultParamHandler` (`getDefaults`, `getParameters`, `setParameters`) | `IDMapper::defaults`, `to_param`, `from_param`, `from_param_with_warnings` |
| parameter `rt_tolerance` / member `rt_tolerance_` | field `rt_tolerance` |
| parameter `mz_tolerance` / `mz_tolerance_` | field `mz_tolerance` |
| parameter `mz_measure` / `measure_` | field `measure` |
| parameter `mz_reference` (read from `param_` at each use) | field `mz_reference: MzReference` |
| parameter `ignore_charge` / `ignore_charge_` | field `ignore_charge` |
| `annotate(AnnotatedMSRun&, const PeptideIdentificationList&, const std::vector<ProteinIdentification>&, bool clear_ids, bool map_ms1)` | `annotate_spectra(&MSExperiment, &mut SpectraAnnotation, &[PeptideIdentification], &[ProteinIdentification], clear_ids, map_ms1)` |
| `annotate(AnnotatedMSRun&, const FeatureMap&, bool, bool)` | `annotate_spectra_from_features` |
| `annotate(FeatureMap&, ids, protein_ids, use_centroid_rt, use_centroid_mz, const PeakMap& spectra)` | `annotate_features` |
| `annotate(ConsensusMap&, ids, protein_ids, measure_from_subelements, annotate_ids_with_subelements, const PeakMap& spectra)` | `annotate_consensus` |
| `struct PeptideIdentificationListState` | `PeptideIdentificationListState` |
| static `mapPrecursorsToIdentifications(spectra, ids, mz_tol = 0.001, rt_tol = 0.001)` | `map_precursors_to_identifications` (explicit tolerances), `map_precursors_to_identifications_default` |
| protected `updateMembers_` | `from_param` (via `DefaultParamHandler::set_parameters_with`) |
| protected `getAbsoluteMZTolerance_` | `absolute_mz_tolerance` (public; the class test reaches it through a subclass) |
| protected `isMatch_` | `is_match` (public, same reason) |
| protected `checkHits_` | private `check_ids` |
| protected `getIDDetails_` | private `id_details` |
| protected `increaseBoundingBox_` | private `increase_box` |
| protected `checkMassType_` | `uses_average_mass` (public, warnings returned) |
| protected static `addIdentificationDataProcessing_` | `add_identification_data_processing` (public) |
| file-local `checkTMTType`, `NATIVE_ID_TYPE` | private `isobaric_type` |
| `OPENMS_LOG_INFO` / `OPENMS_LOG_WARN` statistics and warnings | returned `SpectraMappingStatistics` / `MapMappingStatistics` |
| — | `IDMapperLimits` (native resource ceilings), `SpectraAnnotation` (see below), `ID_MAPPER_NAME`, `UNKNOWN_SEARCH_RUN_IDENTIFIER` |

`METADATA/AnnotatedMSRun.h` is not ported (status `unmapped`, consumer `IDSplitter`). The
spectrum overloads therefore take the experiment read-only and write its identification half
into `SpectraAnnotation` (one peptide identification per spectrum plus the protein runs), which
an `AnnotatedMSRun` port can wrap. `SpectrumLookup` is used only for `findByNativeID`, which
is a native-ID map built in place (last spectrum wins, as `SpectrumLookup::addEntry_`).

### Preserved source conventions

- The spectrum RT sweep is the source's two-pointer walk over both RT-sorted lists (multimap
  order: by RT, ties in insertion order), with its strict comparisons: an identification at
  exactly a spectrum's RT is not mapped with a zero tolerance, while `is_match` accepts that
  distance. Hits are *merged* into the one identification per spectrum, which therefore has no
  identifier or score type of its own.
- Only identifications that are not `PeptideIdentification::empty()` in the source sense
  (identifier, hits, threshold, score type and higher-is-better all default) take part in
  spectrum mapping; feature and consensus mapping skip identifications without hits.
- Feature boxes, the one-second RT hash table, the per-mass-trace hull check, centroid
  fallback when any feature lacks a hull, the charge rules (`charges` contain the feature
  charge for precursor m/z; per-hit charge for peptide m/z; `|z|` equality for unidentified
  precursors) and the average-mass decision from `FeatureFinder` processing records are those
  of the source.
- Consensus mapping: `0` always matches as a charge for identifications; a subelement match
  stores `map_index` on the copy; unidentified precursors are pushed once per matching handle,
  and the one reused empty identification keeps the last `map_index` it was given.
- The identifier of empty precursor identifications is the map's first protein run, else its
  first unassigned identification, else `UNKNOWN_SEARCH_RUN_IDENTIFIER`, in which case a run
  with that identifier and the current date is appended.
- Isobaric (TMT/iTRAQ) consensus maps are mapped by spectrum reference exactly as the source
  does, including `lookForScanNrsAsIntegers` being decided from the *first* identification's
  reference.
- Statistics keep the source's counting quirks (for example a centroid-m/z precursor match on
  a feature map is counted as unassigned; "unmapped (empty)" also counts native-ID mappings).

### Native differences

- **Atomicity.** Every `annotate_*` plans its assignments first and commits only on success;
  the source mutates as it goes and leaves a half-annotated map on an exception.
- **Resource ceilings** (`IDMapperLimits`): RT bins, hash entries, comparisons and copied
  identifications are bounded and checked before the map is touched.
- **Non-finite input** is refused (`Error::InvalidValue`): identification RT/m/z, feature
  coordinates, and spectrum RT or precursor m/z of spectra with precursors. The source only
  checks `isnan` of identifications; a NaN or huge feature RT makes it index its hash table out
  of range.
- A feature that has convex hulls none of which contains a point is refused; the source casts
  `floor(DBL_MAX)` to an index.
- An isobaric consensus feature without handles is refused; the source dereferences the empty
  handle set. A missing column header gives an empty file name without inserting a header; the
  source's `std::map::operator[]` inserts a default one into the map.
- A `spectrum_reference` meta value that is not a string is `Error::InvalidValue`; the source's
  implicit `DataValue` conversion throws `ConversionError`.
- Run dates are text in this crate. `add_identification_data_processing` parses them; an
  unparsable or unset date gives no completion time (the source copies an unset `DateTime`).
- Log output is returned as statistics and warning strings.

### C++ defect candidates

1. `IDMapper.cpp:388-396`: `boost::regex scanregex{""}` followed by `if (lookForScanNrsAsIntegers
   && scanregex.empty())`. A successfully compiled empty pattern is not `empty()`, so
   `getRegExFromNativeID` is never consulted and `extractScanNumber(cf_scan_id, scanregex,
   false)` runs the empty pattern, which has no group 1: its tokens are empty strings, the
   conversion fails and `ParseError` is thrown whenever a scan ID has no exact match and the
   first identification's reference is an integer (Sage-style scan numbers). The port executes
   the same pattern through `concept::boost_regex` and returns `Error::Parse`. Confirmed by an
   executed C++ oracle (`../oracle/I4-idmapper-tmt`, Release build on ibminode06): `boost::regex{""}.empty()` is 0 and the Sage-style case throws `ParseError` "Could not extract scan number".
2. `IDMapper.cpp:351`: the integer check reads `ids[0].getSpectrumReference()` inside the loop
   over `ids[i]`; an unannotated first identification disables scan-number matching for all.
3. `IDMapper.cpp:963-969`: the centroid branch of the unidentified-precursor loop breaks before
   incrementing `matching_features`, so those matches are reported as unassigned.
4. `IDMapper.cpp:386`: `map.getColumnHeaders()[first_channel.getMapIndex()]` inserts a column
   header for an unknown map index; `*cf.getFeatures().begin()` is undefined for an empty
   consensus feature.

### Checked boundaries and evidence

- Tier 1: the consensusXML files the class test compares with (`IDMapper_3_out1`, `_3_out2`,
  `_5_out1`, `_6_out1..3`) and the retained outputs of `TOPP_IDMapper_1..6` (test-data
  `0cb15f2`, `topp/CMakeLists.txt:2963-2992`) are reproduced identification by identification
  (RT, m/z, hit sequences and charges, `map_index`, `spectrum_index`, `spectrum_reference`),
  with the TOPP tool's option defaults (`mz_reference=peptide`, `feature:use_centroid_mz=true`,
  topp `174b576` `src/IDMapper.cpp`) and, for consensus maps, its
  `sortPeptideIdentificationsByMapIndex`. Run identifiers and dates are not compared.
- Tier 1 (executed): the isobaric branch against `../oracle/I4-idmapper-tmt` (four cases).
- Tier 3: every in-memory literal of `IDMapper_test.cpp`.
- Tier 4: the disabled FeatureMap-to-spectra section's scenario, the zero-tolerance sweep,
  precursor partitioning, atomicity and limits.
- Reader defects met on the way (reported, not fixed here): the idXML/featureXML/consensusXML
  readers refuse `0000-00-00T00:00:00` (the source's own text for an unset date), the flanking
  residue `*`, and a repeated `SearchParameters` UserParam. `tests/id_mapper.rs::fixture_text`
  rewrites those in memory; none of them reaches the mapping.

## IonIdentityMolecularNetworking

| Source member | Rust |
|---|---|
| static `annotateConsensusMap(ConsensusMap&)` | `annotate_consensus_map(&mut ConsensusMap) -> Result<()>` |
| static `writeSupplementaryPairTable(const ConsensusMap&, const std::string& output_file)` | `write_supplementary_pair_table(&ConsensusMap, impl Write) -> Result<bool>`, `store_supplementary_pair_table(&ConsensusMap, path) -> Result<bool>` |
| file-local `VertexLabel`, `UndirectedOSMIdGraph` (Boost.Graph) | a vertex list and adjacency vectors; components numbered in first-vertex order as Boost's `connected_components` numbers them over `vecS` |
| — | `MAX_GRAPH_ELEMENTS` |

Preserved: row IDs on every feature, network numbers only on features with linked groups,
partner lists in ascending feature order, groups compared as text, linked groups removed; the
pair table's header, `MS1 annotation`, score `partners(a) + partners(b) - 2` in wrapping
unsigned arithmetic, `default` best ion, `dm/z=` via `StringUtils::toStr(double)`, reverse
edges dropped as written, comma replaced by `_` without quoting. Native: the map is unchanged
on error; a non-list linked-groups value or a non-numeric group (`std::stoull`) is
`Error::InvalidValue`; an empty map writes nothing (the source reads element 0 unchecked); a
partner list naming a row outside the map is an error (the source indexes out of range); the
result reports whether anything was written. Boost.Graph needs no crate here: the one
algorithm used is a breadth-first labelling.

## SiriusFragmentAnnotation

| Source member | Rust |
|---|---|
| nested class `SiriusTargetDecoySpectra` (`target`, `decoy`, constructors) | `SiriusTargetDecoySpectra { target, decoy }`, `new`, `Default` |
| static `extractAndResolveSiriusAnnotations(subdirs, score_threshold, use_exact_mass, decoy_generation)` | `extract_and_resolve_sirius_annotations` |
| static `extractAnnotationsFromSiriusFile(path, max_rank = 1, decoy = false, use_exact_mass = false)` | `extract_annotations_from_sirius_file` (all four arguments explicit) |
| static `extract_columnname_to_columnindex(const CsvFile&)` | `extract_columnname_to_columnindex(&CsvFile)` |
| protected `extractConcatNativeIDsFromSiriusMS_` | `extract_concat_native_ids` |
| protected `extractConcatMIDsFromSiriusMS_` | `extract_concat_m_ids` |
| protected `extractFeatureIDFromSiriusMS_` | `extract_feature_id` |
| protected `extractCompoundRankingAndFilename_` | `extract_compound_ranking_and_filename` |
| protected `extractCompoundRankingAndScore_` | `extract_compound_ranking_and_score` |
| — | `MAX_FILE_BYTES` |

Preserved: file layout and names, `formula_candidates.tsv` read with the crate's `CsvFile`
(tab, source comment handling), adduct whitespace removal, first row winning per rank, the
`_i` suffix for ranks above 1, the meta values and data-array layout (the ionization array stays
unnamed), `use_exact_mass` forced off for decoys, an empty extra-mass array for decoys, the
identifier-only decoy filler, resolution by native ID with ties going to the later directory,
output ordered by native ID. Native: a missing or empty `formula_candidates.tsv` yields an
empty result (the source dereferences the end of its empty map), a missing rank, a short or
empty fragment row and a non-numeric field are `Error::InvalidValue` (the source throws
`std::out_of_range`, reads out of bounds, or throws `ConversionError`), files are capped at
`MAX_FILE_BYTES`, the `decoy` meta value is the integer 0/1 the source's `bool` becomes, and
numbers are parsed by Rust's parser with `StringUtils`' whitespace, `+` and range rules rather
than by `std::from_chars`. Evidence: the class test's three sections (tier 3, its default
`TEST_REAL_SIMILAR` tolerance kept) plus exact values read off the fixture (tier 4).

## MetaboliteSpectralMatching

Rust: `src/analysis/id/metabolite_spectral_matching.rs`; tests:
`tests/metabolite_spectral_matching.rs`; provenance:
`tests/data/metabolite_spectral_matching_provenance.json`. Status **partial**: the
default `merge_spectra = true` needs `PROCESSING/SPECTRAMERGING/SpectraMerger.h`
(`mergeSpectraPrecursors`), which is not ported; `run` returns `Error::Unsupported` for it.

| Source member | Rust |
|---|---|
| class `SpectralMatch` (17 fields, getters/setters, copy/assignment) | `SpectralMatch` with public fields, `Default` (CCS = `ImTypes::DRIFTTIME_NOT_SET`), `Clone` |
| `MetaboliteSpectralMatching()` + parameters `prec_mass_error_value`, `frag_mass_error_value`, `mass_error_unit`, `report_mode`, `ionization_mode`, `merge_spectra`, `ccs_error_percent` | `MetaboliteSpectralMatching` fields, `defaults`, `to_param`, `from_param`, `from_param_with_warnings`; enums `MassErrorUnit`, `ReportMode`, `IonizationMode` |
| static `computeHyperScore(double, bool, const MSSpectrum&, const MSSpectrum&, double)` | `compute_hyper_score` |
| static `computeHyperScore(..., std::vector<PeptideHit::PeakAnnotation>&, double)` | `compute_hyper_score_with_annotations` |
| `run(PeakMap&, PeakMap&, MzTab&, std::string& out_spectra)` | `run(&mut MSExperiment, &mut MSExperiment, &mut MzTab, &str) -> Result<Vec<String>>` (warnings returned) |
| protected `updateMembers_` | `from_param` |
| protected static `computeHyperScore_` | private `hyper_score` |
| private `exportMzTab_` | private `export_mztab` |
| file-local `extractCCS`, `PrecursorMZLess`, `SpectralMatchScoreGreater` | private `extract_ccs`, stable sorts |
| base `ProgressLogger` | not ported (no progress output) |

Preserved: library sorted by precursor m/z in place; sliding `WindowMower` (20 Th, 5 peaks) on
the experiment in place, the thinned spectra stored as mzML when `out_spectra` is given; the
precursor window in the fragment unit; charge-sign skip; CCS from CCS drift time, 1/K0 via
`ImTypes::one_over_k0_to_ccs`, or a numeric `CCS` meta value, checked on the spectrum then the
first precursor; CCS filtering only when both are known; hyperscore `ln(dot) + ln(min(N,170)!)`
clamped at 0 with fewer than three matches scoring 0; RT floored to 0.1 s; identifier fallback
GNPS_Spectrum_ID → Massbank_Accession_ID → Metabolite_Name → native ID; top3/best; every mzTab
cell and optional column of `exportMzTab_` (`opt_ppm_error` and CCS columns through
`StringUtils::toStr(double)`, `opt_match_score` through a default `std::ostream`, i.e. `%g`
at precision 6). Native: `std::sort` ties (library entries with equal precursor m/z, equal scores)
keep input order instead of being unspecified; a library spectrum without precursor, unsorted
spectra, and short annotation arrays are errors (source UB / `std::out_of_range`); a scored
candidate lacking one of the string meta values `HMDB_ID`, `Sum_Formula`, `Inchi_String`,
`SMILES_String`, `Precursor_Ion` is `Error::InvalidValue` as the source's `ConversionError`;
`mztab_out` is untouched on error. `boost::math::factorial<double>` is reproduced by exact integer
arithmetic rounded to nearest (Boost's table holds the correctly rounded values); not executed
against Boost. Evidence: tier 3 for every class-test section (library, parameters and CCS
scenarios transcribed), tier 4 for the exported cells, hyperscore values and error cases
derived from the source formula.
