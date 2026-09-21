# mzIdentML support

Ported headers: `FORMAT/MzIdentMLFile.h` (107 lines) and
`FORMAT/HANDLERS/MzIdentMLHandler.h` (375 lines), with their implementations
`MzIdentMLFile.cpp` (115 lines) and `MzIdentMLHandler.cpp` (2302 lines).

Rust: [`src/format/mzidentml.rs`](../src/format/mzidentml.rs) — one module,
gated on the `idxml` feature, because it needs the identification stack
(`ProteinIdentification`, `PeptideIdentification`, `PeptideHit`,
`PeptideEvidence`, `SearchParameters`) and `quick-xml`, which is exactly what
`idxml` pulls in. No Cargo dependency or feature was added.

Tests: [`tests/mzidentml.rs`](../tests/mzidentml.rs).
Provenance and hashes: [`tests/data/mzidentml_provenance.json`](../tests/data/mzidentml_provenance.json).

`MzIdentMLFile::load` does **not** use the handler in the owned header: it
delegates to `Internal::MzIdentMLDOMHandler` (`MzIdentMLDOMHandler.h/.cpp`,
3258 lines), a separate, unowned header. That file is the read path this module
reproduces, so it was read in full and is hashed in the manifest with its role
recorded; its own ledger entry belongs to a later stage. The handler in the
owned header is the **write** path, plus a stream read path that is dead code
(see the API table).

The read path also finishes a cross-linking document with six functions from
`ANALYSIS/XLMS/OPXLHelper.h` (1474 lines), another unowned header. Those six
are reproduced as private helpers of this module — the read path is not
meaningful without them — and the file is hashed in the manifest with that role
recorded; the rest of `OPXLHelper.h` is a search algorithm and keeps its own
ledger entry. Section 5 has the whole cross-linking picture.

---

## 1. API mapping

### `FORMAT/MzIdentMLFile.h`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `MzIdentMLFile()` | [`SCHEMA_VERSION`], [`ReadOptions::default`], [`WriteOptions::default`] | The constructor's only effect is `XMLFile("/SCHEMAS/mzIdentML1.3.0.xsd", "1.3.0")`. There is no handle type to construct; the schema version is a constant and the ceilings are option structs. |
| `~MzIdentMLFile()` | not ported: `= default` in the source, nothing to release | |
| `void load(filename, poid, peid)` | [`load`], [`load_with_options`], [`load_with_registry`], [`load_into`] | Returns an owned [`MzIdentMLDocument`] instead of filling two out-parameters. `load_into` is the closest analogue and replaces the destination atomically; the source clears both containers first because its DOM handler only appends. |
| `void store(filename, poid, peid) const` | [`store`], [`store_with_options`], [`store_with_registry`] | The `.mzid` extension check is kept and happens before any output. |
| `bool isSemanticallyValid(filename, errors, warnings)` | **not ported**: needs `share/OpenMS/MAPPING/mzIdentML-mapping.xml`, which is not an embedded crate resource, and `FORMAT/VALIDATORS/MzIdentMLValidator.h`, a separate unowned header. | |
| `bool isValid(filename, os, used_version)` | `is_valid`, `is_valid_with_options`, with the optional `xml-schema` feature | [`detect_version`] picks the version and the matching bundled, unchanged `mzIdentML1.{0,1,2,3}.0.xsd` validates the file (1.0.0 with `FuGElightv1.0.0.xsd`, composed in memory). `used_version` is `report.schema.version()`. The retained TOPP_FileInfo_14/15 verdicts (1.1.0, valid and invalid at line 327) are ported; see [XML schema validation](XML_SCHEMA_SUPPORT.md). |
| `std::string detectVersion(filename) const` | [`detect_version`], [`detect_version_from_reader`] | Full port, including the 15-line header window, the `version="x.y.z"` preference over the `mzIdentML/x.y` namespace, and the fallback to the adapter default. |
| inherited `Internal::XMLFile::getVersion` / `isValid` | [`SCHEMA_VERSION`]; `isValid` as above | `MzIdentMLFile` overrides the inherited `isValid` with the version-detecting one above. |
| inherited `ProgressLogger` | nothing to port: the source only inherits it to hand a reference to the handler, which makes no progress call on load or store; the Release build makes none (`mzid_load`, `mzid_store` in `tests/progress_format_readers.rs`). A caller's logger type is `crate::concept::progress_logger::ProgressLogType` | |

### `FORMAT/HANDLERS/MzIdentMLHandler.h` — `Internal::IdentificationHit`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `IdentificationHit()` | [`IdentificationHit::new`] | Keeps the source defaults: `pass_threshold` true, every number 0. `Default` gives `pass_threshold` false, which is Rust's derive, so `new()` is the source-faithful constructor. |
| copy/move constructors, copy/move assignment | `Clone` | |
| `virtual ~IdentificationHit()` | not ported: `= default` | |
| `operator==`, `operator!=` | `PartialEq` | |
| `setId`, `getId` | field `id` | |
| `setCharge`, `getCharge` | field `charge` | |
| `setCalculatedMassToCharge`, `getCalculatedMassToCharge` | field `calculated_mass_to_charge` | |
| `setExperimentalMassToCharge`, `getExperimentalMassToCharge` | field `experimental_mass_to_charge` | |
| `setName`, `getName` | field `name` | |
| `setPassThreshold`, `getPassThreshold` | field `pass_threshold` | |
| `setRank`, `getRank` | field `rank` | 0-based, as the source stores `rank - 1`. |
| inherited `MetaInfoInterface` | field `metadata` ([`MetaInfo`]) | |

### `MzIdentMLHandler.h` — `Internal::SpectrumIdentification`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `SpectrumIdentification()` | `Default` | |
| copy/move constructors, copy/move assignment | `Clone` | |
| `virtual ~SpectrumIdentification()` | not ported: out-of-line but empty | |
| `operator==`, `operator!=` | `PartialEq` | |
| `setHits`, `getHits` | field `hits` | |
| `addHit` | [`SpectrumIdentification::add_hit`] | |
| protected `id_` | field `id` | The source has no accessor for it, which makes the member unreachable; exposing it costs nothing and the type is otherwise write-only. |
| inherited `MetaInfoInterface` | field `metadata` | |

### `MzIdentMLHandler.h` — `Internal::Identification`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `Identification()` | `Default` | |
| copy/move constructors, copy/move assignment | `Clone` | |
| `virtual ~Identification()` | not ported: out-of-line but empty | |
| `operator==`, `operator!=` | `PartialEq` | |
| `setCreationDate`, `getCreationDate` | field `creation_date` | `Option<String>` of the `xs:dateTime` text rather than a parsed `DateTime`, so an unparsable date is preserved instead of silently becoming the epoch. |
| `setSpectrumIdentifications`, `getSpectrumIdentifications` | field `spectrum_identifications` | |
| `addSpectrumIdentification` | [`Identification::add_spectrum_identification`] | |
| protected `id_` | field `id` | As above: no accessor in the source. |
| inherited `MetaInfoInterface` | field `metadata` | |

All three types are ported as value types and **nothing in this module consumes
them**, because the only code that fills them is the dead stream read path
below. They are part of the owned header's public API, so they are ported; they
are not part of the load/store contract.

### `MzIdentMLHandler.h` — `Internal::MzIdentMLHandler`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `MzIdentMLHandler(const vector<ProteinIdentification>&, const PeptideIdentificationList&, filename, version, logger)` (write) | [`write`], [`write_with_options`], [`write_with_registry`] | The constructor plus `writeTo` is one function here; there is no handler object to hold the inputs. |
| `MzIdentMLHandler(vector<ProteinIdentification>&, PeptideIdentificationList&, filename, version, logger)` (read) | **not ported** | This is the read-mode constructor of the dead stream path; the reader is [`read_with_registry`], which reproduces the DOM handler instead. |
| `~MzIdentMLHandler()` | not ported: `= default` | |
| `void onStartElement(qname, attributes)` | **not ported** | Vestigial. It parses `Peptide`, `Modification` and `SpectrumIdentificationItem` into `current_id_hit_`/`actual_peptide_`, and `onEndElement` moves the finished hit into `current_spectrum_id_`, which nothing ever stores; `id_`, `pro_id_` and `pep_id_` are never written. The header documents this class as appending identifications in read mode, which it does not do. Recorded in `OpenMS_CPP_ISSUES.md`. |
| `void onEndElement(qname)` | **not ported**, as above | |
| `void onCharacters(chars, length)` | **not ported**, as above | |
| `void writeTo(std::ostream&)` | [`write_with_registry`] | Full port of the element structure, with the id scheme, the number formatting, the C-terminal modification location and the `userParam` type attribute diverging (section 3). |

#### Protected and private members of `MzIdentMLHandler`

Not public API, but they define the write behaviour, so they are accounted for
here too.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `logger_` | nothing to port | The handler never calls it; see `ProgressLogger` above. |
| `cv_` | [`ControlledVocabulary::psi_ms`] | Loaded once from the embedded `psi-ms.obo` rather than `File::find("/CV/psi-ms.obo")`. |
| `unimod_` | not ported | The crate has no UniMod OBO; UniMod accessions and names come from the `ModificationsDB` record itself, which carries the same pair. |
| `tag_`, `open_tags_` | not ported | Parser state of the dead stream path. |
| `id_`, `pro_id_`, `pep_id_`, `cid_`, `cpro_id_`, `cpep_id_` | function parameters | |
| `current_spectrum_id_`, `current_id_hit_` | not ported | State of the dead stream path. |
| `handleCVParam_` | not ported | Only reachable from the dead stream path, where it validates a UNIMOD `Modification` cvParam and otherwise does nothing. |
| `handleUserParam_` | not ported | Declared in the header and **not defined anywhere** in the source. |
| `writeMetaInfos_` | `write_meta` | Writes the `type` attribute the schema defines instead of the source's `unitName` (section 3). |
| `getChildWithName_` | not ported | Declared in the header and not defined anywhere in the source (the DOM handler has its own copy). |
| `writeEnzyme_` | `write_enzyme` | Same CV fallback chain: the enzyme name, then `NoEnzyme` for "no cleavage", then `cleavage agent details`. |
| `writeModParam_` | `write_mod_params` | Emits protein-terminal specificity rules as well (section 3). |
| `writeFragmentAnnotations_` | `write_fragmentation`, `split_annotation` | The regex is replaced by an explicit parser with the same accepted shape; the `is_ppxl` flag is the `crosslinking` parameter and writes the same `cross-link_chain` / `cross-link_ioncategory` arrays. |
| `trimOpenMSfileURI` | `trim_file_uri` | |
| `writePeptideHit` | `write_item`, `write_modifications`, `write_score` | |
| `writeXLMSPeptideHit` | `plan_crosslinks`, `write_crosslink_peptide`, `write_crosslink_item`, `crosslinker_term`, `crosslink_skip` | Full port of the cross-linking output path (section 5). The heavy half of a labelled pair is built directly rather than by substituting strings in the light one, and the location of a terminal link follows the hit's specificity rather than whether `CrossLinksDB` happens to hold a terminal record of that mass. |
| `is_ppxl` (local of `writeTo`) | `Plan::crosslinking`, `is_crosslinking_run` | Same test: `is_cross_linking_experiment` or `SpectrumIdentificationProtocol == MS:1002494` on a run. |
| `initCvCaches_` | `WriteContext::details` | `getAllChildTerms("MS:1001143")`, descendants only. |
| `pep_sequences_`, `pp_identifier_2_sil_`, `sil_2_sdb_`, `sil_2_sdat_`, `ph_2_sdat_`, `sil_2_sip_`, `peptide_result_details_` | `Plan`, `WriteContext` | The same id and reference bookkeeping, with positional ids. |
| `actual_peptide_`, `current_mod_location_`, `actual_protein_` | not ported | State of the dead stream path. |
| private default constructor, copy constructor, copy assignment | not applicable | Deleted-by-privacy idiom; Rust functions need no such guard. |

### Read path reproduced from `MzIdentMLDOMHandler` (unowned header)

Listed because `MzIdentMLFile::load` is only meaningful through it. Every entry
is reproduced unless stated.

| C++ member | Rust counterpart |
|---|---|
| `readMzIdentMLFile` | [`read_with_registry`] |
| `writeMzIdentMLFile` | not ported (the source's `store` uses the stream handler; the DOM writer is commented out at `MzIdentMLFile.cpp:51`) |
| `parseParamGroup_`, `parseCvParam_`, `parseUserParam_` | `param_group`, `CvParam::from_node`, `user_value` |
| `parseAnalysisSoftwareList_` | `read_software` |
| `parseDBSequenceElements_`, `parsePeptideElements_`, `parsePeptideEvidenceElements_` | `read_sequence_collection`, `read_peptide` |
| `parseSpectrumIdentificationElements_` | `read_runs` |
| `parseSpectrumIdentificationProtocolElements_` | `read_protocols`, `modification_params` |
| `parseInputElements_` | `read_inputs` |
| `parseSpectrumIdentificationListElements_` | `read_lists`, `apply_result_params` |
| `parseSpectrumIdentificationItemElement_` | `read_item`, `select_score`, `merge_target_decoy`, `attach_evidences` |
| `parseSpectrumIdentificationItemSetXLMS` | `read_crosslink_result`, `read_crosslink_group`, `read_crosslink_fragmentation`, `crosslink_user_value` |
| XL half of `parsePeptideSiblings_` | `read_crosslink_modification`, `apply_crosslink_residue_modification` |
| `parseProteinDetectionListElements_`, `parseProteinAmbiguityGroupElement_`, `parseProteinDetectionHypothesisElement_` | `read_protein_detection` |
| `inferModificationLocation_` | `infer_modification_location` |
| `parsePeptideSiblings_` | `read_peptide`, `apply_modification` |
| `findSearchParameters_` | `additional_search_params` |
| `initScoreTermCaches_` | `ScoreTerms::new` |
| `toDoubleOrNaN_`, `toDoubleOrZero_` | `optional_value`, `score_value` |
| `buildCvList_` and the other `build*` writers | not ported (they belong to the commented-out DOM writer) |

### `ANALYSIS/XLMS/OPXLHelper.h` (unowned header), post-processing only

`readMzIdentMLFile` finishes a cross-linking document with six calls into this
header. The functions are reproduced as private helpers of this module, because
the read path is not meaningful without them; `OPXLHelper.h` as a whole - its
candidate enumeration, database digestion and scoring - stays unported and
keeps its own ledger entry.

| C++ member | Rust counterpart |
|---|---|
| `addProteinPositionMetaValues` | `add_protein_position_meta_values` |
| `addBetaAccessions` | `add_beta_accessions` |
| `addXLTargetDecoyMV` | `add_crosslink_target_decoy` |
| `removeBetaPeptideHits` | `merge_crosslink_hits` |
| `computeDeltaScores` | `compute_delta_scores` |
| `addPercolatorFeatureList` | `add_percolator_features`, `PERCOLATOR_FEATURES` |
| every other member (`enumerateCrossLinksAndMasses`, `digestDatabase`, `buildCandidates`, `buildPeptideIDs`, `combineTopRanksFromPairs`, …) | not ported: not on the read or write path |

---

## 2. Preserved source conventions

These are deliberate and observable, and several are surprising:

* **Elements are collected document-wide by tag name.** The source uses
  `DOMDocument::getElementsByTagName`, which ignores nesting, so an element
  under the wrong parent is still read. `gather` reproduces that, in document
  order.
* **`cvParam`s are keyed by accession in a sorted map, and that order decides
  the score.** The PSM score type is the first matching accession in
  *lexicographic* order, not document order: an MS-GF+ item that writes
  `MS:1002049`, `MS:1002050`, `MS:1002052` and `MS:1002053` yields
  `MS-GF:RawScore` because `MS:1002049` sorts first among the `MS:1001143`
  descendants. Pinned by `the_score_scan_follows_lexicographic_accession_order`.
* **Score selection order**: q-value child terms (skipping `MS:1002055`,
  without stopping the scan), then `MS:1001143` descendants (score type = the
  file's `name` attribute, orientation from the CV), then e-value child terms
  ("E-value"), then the bare `MS:1001143` (higher is better, assumed). An item
  with none of these yields **no hit**, which is why a
  `SpectrumIdentificationResult` can legitimately produce an empty
  `PeptideIdentification`.
* **`SearchParameters::db` is the `SearchDatabase` *location*,** not the
  `DatabaseName`. `whole.mzid` gives `MSDB`, `msgf_mini.mzid` `database.fasta`.
* **`missedCleavages` below zero becomes 1000** ("assuming unlimited").
* **Tolerances take the numerically greater of the plus/minus bounds,** and are
  ppm when the unit *name* is `parts per million`.
* **A threshold is the first `MS:1002482` descendant in accession order,** and
  `MS:1001494` ("no threshold") stops the scan without setting one.
* **`MinCharge`/`MaxCharge` win over an explicit `charges` list,** and are
  rendered as `min-max`.
* **mzIdentML rank 0 is treated as rank 1** (PMF data), then converted to
  OpenMS's 0-based rank.
* **Every `PeptideEvidence` of the referenced `Peptide` is attached to the
  PSM,** not just the ones its own `PeptideEvidenceRef` children name; the
  source calls those references redundant. A PSM therefore sees every protein
  the peptide maps to.
* **`isDecoy` merges into `target_decoy`** as target / decoy / target+decoy
  across a PSM's evidences.
* **`ProteinDetectionList` hypotheses are appended to the last run,** whichever
  run the list belongs to, and without checking for an existing accession. The
  element carries no reference that would identify the right run, so any other
  choice would be invented.
* **Protein hits are ordered by `(score, accession)`,** descending when higher
  scores are better: that is `ProteinHit::ScoreMore`, and it is why
  `sp|P0A9K9|SLYD_ECOLI` precedes `sp|P0A786|PYRB_ECOLI` in `msgf_mini` even
  though the latter is inserted first.
* **Modification-location inference** for a `Modification` with no `location`:
  N-terminus if the modification has any N-terminal variant, `length + 1` if it
  is exclusively C-terminal, nothing at all when `residues` names a concrete
  amino acid or the modification may also be internal.
* **`ModificationsDB::getModification`'s two-stage lookup** is reproduced in
  `resolve_modification`: with a residue and no requested specificity, an
  `ANYWHERE` lookup is tried first, "to avoid ambiguities (e.g.
  `Carbamidomethyl (N-term)`/`Carbamidomethyl (C)`)".
* **The writer's `SpectraData` format map** has four entries (mzML, mzXML,
  mzData, MGF) and defaults to mzML, including for its own `UNKNOWN`
  placeholder.
* **The writer's software CV fallback**: the OMSSA / Mascot / X\!Tandem /
  Sequest / MS-GF+ / Percolator / OpenPepXL aliases, then the engine name if
  PSI-MS has a term with it, then `analysis software`. A search engine PSI-MS
  does not know therefore comes back as `analysis software` after a round trip.
* **`passThreshold`** comes from the hit's `pass_threshold` meta value, else
  from comparing the score against the run's significance threshold when one
  was set, else true.
* **An empty `spectra_data`** is written as `location="UNKNOWN"`, and an
  OpenMS file URI keeps the source's `[`/`]` trimming and backslash
  normalisation.
* **A missing spectrum reference** falls back to `MZ:<mz>@RT:<rt>`.

### Cross-linking conventions

* **`MS:1002494` anywhere in any `AdditionalSearchParams` switches the whole
  document** to the cross-linking path, before anything else is read, and tags
  every run with the term - which is also what the writer tests.
* **Items of one result are grouped by the value of their `MS:1002511`
  cvParam,** and each group is one match. A result with no such cvParam is one
  group over its **first** item only ("fix for label-free mono-links").
* **The light half is the item with the smallest experimental m/z**; a match
  counts as labelled only when the m/z differ *and* the `spectrumID` names more
  than one spectrum, so a single spectrum whose items merely disagree on m/z is
  not mistaken for a label pair.
* **The alpha chain is the peptide carrying the donor modification**; a match
  with as many beta items as alpha ones is a `cross-link`, one whose alpha
  peptide carries the donor and the acceptor of the same link is a
  `loop-link`, anything else a `mono-link`.
* **A group with no donor at all falls back to the linear item path,** so a
  noncovalent association is read as ordinary candidates.
* **`rank` and `chargeState` come from the first item of the group** that
  declares a nonzero one; the score from the last `MS:1002681` or `MS:1003024`.
* **Fragment annotations come from the first item of the group** that has a
  `Fragmentation` block, and are rebuilt as
  `[<chain>|<category>$<series><index><loss>]`.
* **A cross-linker cvParam is not applied to the sequence.** `DSS` and the
  other XLMOD names are not in `ModificationsDB`, so the source's `has` guard
  skips them with a warning; the position and mass live on in the
  `xl_pos1`/`xl_mass` user parameters instead.
* **`userParam`s of an XL item are collected from its descendants,** which is
  how the `cross-link_chain` and `cross-link_ioncategory` arrays of its own
  `IonType`s become hit metadata.
* **One identification per spectrum reference** comes out of the merge step,
  with the beta chain folded into the alpha hit's `BetaPepEv:` values, and the
  `OpenPepXL:score` score type.
* **The first run gets the Percolator feature list** (`feature_extractor`,
  `extra_features`) whether or not it carries a cross-link.

---

## 3. Native differences

Each one is documented at the item in the module as well.

| Difference | Source | This port |
|---|---|---|
| **XML entity references and CDATA** | Xerces resolves both into the element's text | quick-xml reports each `&…;` as its own event and CDATA as another; both are resolved into the text, and an external entity - which would need the DTD this reader refuses anyway - is an [`Error::Unsupported`]. The sibling Mascot XML reader ignored those events and so deleted the character they stand for, turning a score of `1.5` into `15`. Pinned by `entity_references_and_cdata_are_character_data` and `an_external_entity_reference_is_refused`. |
| **Cross-link group index** | Groups by an index that counts *every* element child of the result, then looks the item up in the result's item list, so a result with a non-item element child before its items groups the wrong item | Groups by the item's own ordinal |
| **A result with no items** | Registers the group `(0 → item 0)` unconditionally and then dereferences a null `item(0)` | No identification, no crash |
| **A cross-link group with no experimental m/z** | `min_element` over a vector of NaN, then `light[0]` on an empty index vector | The group is skipped. Pinned by `a_crosslink_group_without_an_experimental_mz_is_skipped` |
| **The result's own parameters** | Applied to `pep_id_->back()`, so with several groups per result only the last one gets the retention time, and with no group at all they land on an unrelated earlier identification | Applied to every identification the result produced |
| **A loop-link's second position** | Set from the acceptor and then overwritten with the mono-link placeholder `"-"` two branches later, which makes the source's own loop-link branch dead and drops the second half of every loop-link on a store | Kept, so a loop-link round-trips. Pinned by `a_loop_link_keeps_both_positions_on_one_chain`. `OpenMS_CPP_ISSUES.md` |
| **`BetaPepEv:start` / `:end` (preserved convention)** | `StringUtils.h` provides the numeric `std::string += Int` overload, which appends decimal digits | Decimal numbers, matching the source; the former character-narrowing claim was incorrect |
| **Hits of a non-cross-linked identification** | `removeBetaPeptideHits` keeps only the first hit of every identification, so a spectrum read through the fallback path loses every candidate but the best | Only a folded beta chain is removed. Pinned by `noncovalent_association_resolves_every_sequence` |
| **The merged identification's run link** | Left empty, so every cross-linking PSM ends up unlinked from its protein run and the writer has to fall back to the first list | Kept |
| **A terminal cross-link's location** | Written at `location=0` (or the position plus two) only when `CrossLinksDB` also holds a terminal record of that mass, and otherwise emits an attribute list with no element name in front of it | The location follows the hit's terminal specificity, so `N_TERM` and `C_TERM` round-trip. Pinned by `terminal_crosslink_positions_round_trip` |
| **The beta chain's C-terminal location** | Written at the peptide length plus two, which its own reader reads as an internal position one past the end | The position plus two, i.e. the schema's C-terminus. `OpenMS_CPP_ISSUES.md` |
| **An ambiguous cross-linker mass** | Falls back to `getModification(one-letter code, full id, ANYWHERE)` - the name and residue arguments the wrong way round, which throws | The first record of that mass, as the intent reads |
| **`xl_pos1_protein` / `xl_pos2_protein`** | `evidence start + link position + 1`, where the source's own reader stored the file's 1-based `start`, so the protein coordinate is one too high | The same formula on the 0-based start this reader stores, i.e. the coordinate the source's comment describes |
| **XML ids** | `UniqueIdGenerator` values, so two stores of one document differ; elements are emitted from `std::set<std::string>`, i.e. sorted by their own XML text | Positional (`SIL_0`, `PEP_3`, `PEV_7`), emitted in a defined order. `store` is byte-reproducible, asserted in `store_and_reload`. The upstream FuzzyDiff whitelist exempts `id=` for exactly this reason. |
| **Run identifier** | A fresh `UniqueIdGenerator` value per run (with the source's own `TODO setIdentifier to xml id?`) | The `SpectrumIdentification` element's `id`, so a load is reproducible and peptide-to-run links survive. |
| **Wall-clock stamps** | An absent `activityDate` becomes `DateTime::now()` on read, an invalid run date becomes `DateTime::now()` on write, and the root `creationDate` is always `DateTime::now()` | An absent date stays `None` on read and the attribute is omitted on write; `WriteOptions::creation_date` sets one explicitly. The upstream `load` section asserts a nonzero date, which only holds because of the substitution; the Rust test asserts `None` and comments why. |
| **C-terminal modification location** | Written at `location = size`, which the source's own reader then applies to the last residue | Written at `size + 1`, the schema's C-terminus. Pinned by `a_c_terminal_modification_is_written_at_length_plus_one`. `OpenMS_CPP_ISSUES.md`. |
| **`PeptideEvidence` positions** | Read as the file's 1-based value into a member documented as 0-based, and written as position + 1, so each store/load cycle shifts them | Converted on both sides; `start="115"` reads as 114 and writes back as 115. `OpenMS_CPP_ISSUES.md`. |
| **`userParam` type** | The XSD type goes into `unitName`, which no reader reads | The schema's `type` attribute, so integers and doubles survive a round trip. `OpenMS_CPP_ISSUES.md`. |
| **`pass_threshold` storage** | `setMetaValue(String, bool)` promotes the bool to an integer, which `DataValue::toBool` then rejects | The crate's documented boolean strings `"true"`/`"false"`, so `MetaValue::to_bool` works. |
| **Number formatting** | `String(double)`, six significant digits | Shortest round-trip decimals, so an m/z or a score is not truncated by writing it. |
| **`nan` in required attributes** | Writes the literal `nan` for an absent m/z, which is not a valid `xsd:double` | Writes `NaN`, the schema's lexical form. An absent RT or m/z is `None`, never a NaN sentinel. |
| **Protein-terminal search modifications** | `writeModParam_` writes specificity rules for the peptide termini only, leaving a `@TODO: handle protein C-term/N-term`, so a protein-terminal modification cannot be resolved unambiguously on re-read | All four rules are written (`MS:1001189`, `MS:1001190`, `MS:1002057`, `MS:1002058`), so `Xlink:DTSSP[88] (Protein N-term)` round-trips. |
| **Ambiguous modification names** | `searchModifications` returns a set and every match is written, multiplying the element | Refused as ambiguous by the registry, so the output cannot silently gain modifications. |
| **`NumTolerableTermini` / `taxonomy`** | Consumed on read and dropped by the writer, losing the enzyme specificity | Written back as `userParam`s, so both survive a round trip. |
| **`AdditionalSearchParams` ordering** | Assigns the parsed block *over* the accumulated parameters, so a document that places `AdditionalSearchParams` after `ModificationParams`, `Enzymes` or the tolerances silently loses them | Seeded from `AdditionalSearchParams` first, then the other children, so the result does not depend on child order. The schema fixes that order, so no conformant file is affected. |
| **Namespaces** | `setDoNamespaces(false)`: tag names are matched as written, so a document that binds the mzIdentML namespace to a prefix is read as empty and then rejected for having no `SpectraData` | The default namespace and prefixed names both resolve; the root's namespace must be an mzIdentML one, and an element in a foreign namespace is refused. Pinned by `a_prefixed_namespace_is_resolved` and `a_foreign_root_namespace_is_refused`. |
| **Malformed XML** | A Xerces parse failure is logged and the partial document is used | Any XML error is an [`Error::Parse`]. |
| **Peptide sequence failures** | An unknown modification, a bad substitution location or an empty `PeptideSequence` ends in a logged warning and an *empty* `AASequence` in the peptide map, so the affected PSMs silently lose their sequence; one of those paths is an out-of-bounds write | Every one of them is an error. |
| **Protein hit order after load** | `ProteinIdentification::sort` with the `(score, accession)` key | Reproduced inside the reader, because the crate's `sort()` compares the score alone. Recorded as a gap against `src/identification.rs`. |
| **Evidence deduplication on write** | One `Peptide` and one evidence list per *sequence string*; a later PSM with the same sequence reuses the first PSM's evidence list and its own evidences are dropped | The `Peptide` element is still shared, but evidences are deduplicated by their full content, so no evidence is lost and the round trip is stable. |
| **A PSM whose evidence names an undeclared protein** | Logged and skipped, dropping the evidence | [`Error::InvalidValue`] before any output. |
| **A peptide identification not linked to a run** | Logged and dropped, losing the whole spectrum | [`Error::InvalidValue`] before any output. |
| **An empty `SpectrumIdentificationResult`** | Emitted, which the schema forbids | [`Error::MissingInformation`]. |
| **Inferred search modifications** | `ModificationDefinitionsSet::inferFromPeptides` when both lists are empty | Not inferred: writing search parameters the input never declared invents metadata, and the modifications are already in the `Peptide` elements. |

### Reference-resolution table

The behaviour the work package asks about, case by case.

| Case | Source | This port |
|---|---|---|
| duplicate element id (same kind) | `std::map::insert` keeps the **first** silently, so the second element's data is lost | [`Error::Parse`] naming the kind and id |
| dangling `peptide_ref` on a `SpectrumIdentificationItem` | `operator[]` default-constructs an empty `AASequence`: the PSM keeps its score and loses its sequence | [`Error::Parse`] |
| dangling `dBSequence_ref` on a `PeptideEvidence` | `operator[]` default-constructs: the PSM gets an empty protein accession and the run gains a `ProteinHit` with an empty accession | [`Error::Parse`] |
| `DBSequence` with an empty `accession` | Not registered at all, so every evidence pointing at it behaves as the case above | [`Error::Parse`] at the `DBSequence` |
| `SpectrumIdentificationList` no `SpectrumIdentification` references | `operator[]` yields index 0: its results are attributed to the **first** run | [`Error::Parse`] |
| two `SpectrumIdentification`s referencing one list | `insert` keeps the first, so the second run silently gets no spectra | [`Error::Parse`] |
| dangling `spectrumIdentificationProtocol_ref` | No protocol matches: the run keeps no search engine or parameters | **Tolerated**, same result — upstream `msgf_mini` does this |
| dangling `searchDatabase_ref` (`SpectrumIdentification` or `DBSequence`) | Empty database location and version; on a `DBSequence` the reference is never consulted | **Tolerated**, same result — upstream `msgf_mini` does this |
| dangling `spectraData_ref` | Empty `spectra_data` entry, later written as `location="UNKNOWN"` | **Tolerated**, same result |
| unknown `analysisSoftware_ref` | Empty search engine and version | **Tolerated**, same result |
| `AnalysisSoftware` without a name or version | Dropped with a logged error | **Tolerated**, dropped |
| `PeptideEvidenceRef` children | Ignored; evidences come from the `Peptide` id instead | Ignored, same |

---

## 4. Checked boundaries and evidence

### Resource ceilings

[`ReadOptions`] is checked before anything is allocated, and every stage draws
from one shared budget, so a refusal leaves the caller's document untouched:

| Field | Default | Bounds |
|---|---|---|
| `max_xml_bytes` | 64 MiB | Encoded input, BOM included; the reader takes `limit + 1` bytes and refuses on overflow |
| `max_elements` | 2,000,000 | Total elements, `cvParam`/`userParam` included |
| `max_depth` | 64 | Element nesting |
| `max_payload_bytes` | 256 MiB | Decoded tree plus the identifications built from it |
| `max_work` | 50,000,000 | Parser events, attribute normalisation and reference resolution |
| `max_list_items` | 1,000,000 | One library, evidence, hit or parameter list |

[`WriteOptions`] bounds the output (`max_output_bytes`, 64 MiB) and the element
count (`max_records`, 1,000,000). `max_records` is enforced twice: while the
plan is still growing - every planned element becomes at least one emitted one,
so the ceiling applies before the plan can outgrow it - and again on each
element as it is written. The whole document is built in memory and handed to
the writer only on success, so a refused write produces no bytes; `store`
publishes atomically through `path_io::write_plain`. [`MAX_ITEMS`] (1,000,000)
caps the identifications a single write accepts.

### No panics on file-derived data

* No indexing or slicing of file-derived text: `location`, `start`, `end` and
  every other numeric attribute goes through `checked_sub`/`try_from` and is
  compared against the sequence length before use. A single-character attribute
  is read with `chars()`, never `text[0]`, so `pre=""` is an error rather than
  a null character, and a multi-byte character is not split.
* Non-ASCII input is tested: `non_ascii_text_survives_the_round_trip` reads and
  writes an accession containing `日本語`.
* Every arithmetic step on a count is `checked_*` or `saturating_*`.
* The cross-linking path adds no index into file-derived data either: a group's
  item indices come from the item list's own range, a link position is compared
  against the chain's length before it becomes a residue index, and the five
  parallel `BetaPepEv:` lists are walked to the length of the shortest instead
  of by the first one's length, which is what the source indexes them by.
* `unsafe` is forbidden crate-wide; nothing here needs it.

### Test-to-section mapping

`MzIdentMLFile_test.cpp` has 15 `START_SECTION`s. All 15 are ported; none is
mapped and none is unaccounted.

| # | Upstream section | Assertion macros | Status | Rust test |
|---|---|---|---|---|
| 1 | `MzIdentMLFile()` | 1 | ported | `adapter_defaults_match_the_source_constructor` (`SCHEMA_VERSION == "1.3.0"`) |
| 2 | `~MzIdentMLFile()` | 0 | ported | same test; nothing to release |
| 3 | `void load(...)` | 41 | ported | `load_msgf_mini_matches_upstream_literals`, `load_replaces_rather_than_accumulates` (`peptide_ids[0]` score 195 with score type `MS-GF:RawScore`) |
| 4 | `[EXTRA]` modification without `location` (#5443) | 16 | ported | `missing_modification_location_is_inferred_or_skipped` (`PEPTIDEC(Carbamidomethyl)K`) |
| 5 | `void store(...)` | 69 | ported | `store_round_trip_preserves_whole_document`, `store_round_trip_preserves_modified_peptides`, `store_rejects_a_foreign_extension` (`variable_modifications.back() == "Acetyl (N-term)"`) |
| 6 | `[EXTRA] multiple runs` | 5 | ported | `three_runs_round_trip` (`precursor_tolerance == Ppm(20.0)` after the round trip) |
| 7 | `[EXTRA] thresholds` | 9 | ported | `thresholds_round_trip` (`significance_threshold == 0.5`, `pass_threshold == "false"` for both hits of spectrum 17) |
| 8 | `[EXTRA]` regression load of the example files | 0 | ported | `every_upstream_non_crosslinking_fixture_loads` plus the six cross-linking fixtures below |
| 9 | `[EXTRA] compability issues` | 0 (entirely commented out) | ported | `misplaced_elements_in_a_param_group_are_ignored`, `a_psm_without_a_recognised_score_yields_no_hit`, `a_psm_without_peptide_evidence_still_loads`, `an_identification_without_rt_keeps_no_coordinate`, `evidence_without_positions_keeps_them_unknown` — one test per condition its comments enumerate |
| 10 | `[EXTRA] XLMS data labeled cross-linker` | 40 | ported | `labelled_crosslinks_load_and_round_trip` (`xl_pos1` 3 / `xl_pos2` 4, `sequence_beta == "SAVIKTSTR"`, `spec_heavy_RT == 2125.5966796875`, score `-0.190406834856118`, annotation `[alpha|xi$b4]` at index 8, ten identifications, the 32-feature `extra_features` list) |
| 11 | `[EXTRA] XLMS data unlabeled cross-linker` | 40 | ported | `unlabelled_crosslinks_round_trip` (three identifications, `mono-link`/`cross-link`/`mono-link`, `KNVPIEFPVIDR` × `LGCKALHVLFER` at 0/3, `xl_mod == "DSS"`, five annotations with charges 1/1/1/1/2, `VEPSWLGPLFPDK(Xlink:DSS[156])TSNLR` at 12) |
| 12 | `[EXTRA]` mzIdentML 1.3 crosslinking scores and thresholds | 8 | ported | `crosslinking_v1_3_round_trips_at_schema_version_1_3_0` (every identification has a spectrum reference and a sequence; the stored file declares `version="1.3.0"` and loads again) |
| 13 | `[EXTRA]` mzIdentML 1.3 noncovalent association | 4 | ported | `noncovalent_association_resolves_every_sequence` (both candidates of the association, with their modifications resolved) |
| 14 | `[EXTRA]` mzIdentML 1.3 EDC crosslinking | 3 | ported | `edc_crosslinking_reads_linked_and_linear_peptides` (16 identifications from a file that also carries a CDATA `SiteRegexp`) |
| 15 | `[EXTRA]` mzIdentML 1.3 multiple spectra per identification | 3 | ported | `multiple_spectra_keep_distinct_references` |

Sections 3, 4, 5, 7, 10, 11 and 12 exceed five assertion macros; all seven are
ported rather than mapped, with their literals transcribed into the Rust tests.

Native tests beyond the upstream suite: the reference-resolution table above
(six refusal tests plus `dangling_metadata_references_stay_tolerated`), the
version-detection branches, the parser boundaries (non-ASCII, prefixed
namespace, foreign namespace, entity references, CDATA, an external entity, a
`DOCTYPE`, every ceiling, empty `PeptideSequence`, out-of-range substitution
location, negative rank), the score-order test, the C-terminal modification
location, the `ProteinDetectionList` rule, the `Fragmentation` round trip, the
three write refusals, the write ceilings including the planning-time record
ceiling, and four cross-linking cases no fixture reaches (a terminal link, a
loop-link, a group with no experimental m/z, a `userParam` typed by either
attribute). 55 tests in total.

### Evidence tier

Tier 3 (source review) for everything transcribed from the class test and the
nine upstream fixtures; tier 4 (independently derived) for the synthetic
documents, the ceilings and the error-variant choices. No C++ mzIdentML reader or writer was built or
executed and no differential output from those operations was retained, so this is **not** a tier 1
differential. Upgrading it needs an oracle driver under `../oracle/` that links
`libOpenMS` and prints the loaded identifications; the natural first comparison
is `msgf_mini`, whose values are all pinned above.

---

## 5. The cross-linking path

`MzIdentMLFile::load` takes a completely different path when any
`AdditionalSearchParams` declares `MS:1002494` ("crosslinking search"), and so
does this port.

### Reading

1. The `Peptide` elements are read with the XL branch of
   `parsePeptideSiblings_` (`read_crosslink_modification`): a `MS:1002509`
   ("crosslink donor") cvParam registers the peptide as an alpha chain under
   the link value the cvParam carries, together with the modification's
   `monoisotopicMassDelta` and the first UNIMOD or XLMOD name beside it; a
   `MS:1002510` ("crosslink acceptor") registers a beta chain. A peptide whose
   only cross-linking evidence is an `Xlink…`/`XLMOD:…` modification is
   registered as its own mono-link. Positions are `location - 1`, so an
   N-terminal link is `-1` and a C-terminal one the peptide length.
2. `parseSpectrumIdentificationItemSetXLMS` (`read_crosslink_result`,
   `read_crosslink_group`) groups the two to four
   `SpectrumIdentificationItem`s of one match by the value of `MS:1002511` and
   builds one identification per group: the light item's m/z and retention
   time, the alpha chain as the first `PeptideHit` and the beta chain as the
   second, the OpenPepXL scores and user parameters, the fragment annotations,
   and `xl_pos1`/`xl_pos2` with their terminal specificities.
3. The six `OPXLHelper` steps then run in the source's order
   (`finish_crosslinks`): protein-coordinate positions, the beta accessions,
   the per-chain target/decoy state, the merge that folds the beta chain into
   the alpha hit and collects the identifications of one spectrum, the delta
   scores, and the Percolator feature list on the first run.

The result is one identification per spectrum reference with one hit per
cross-link, which is what upstream sections 10 and 11 assert.

### Writing

A run that declares the term (or `is_cross_linking_experiment`) selects the
output path (`is_crosslinking_run`). Every chain of every hit becomes its own
`Peptide` element with the cross-linker `Modification` and the donor or
acceptor cvParam that pairs them — the source makes its peptide identity key
unique per hit by appending the link id, which leaves its own de-duplication
branch dead — the beta chain's `PeptideEvidence` elements are rebuilt from the
`BetaPepEv:` values the merge step left behind, all items of one spectrum share
a single `SpectrumIdentificationResult`, and a labelled match writes its heavy
half as a second item at the heavy m/z, the heavy retention time and the
calculated m/z shifted by `cross_link:mass_isoshift`.

The cross-linker accession comes from the bundled XLMOD vocabulary through
[`CrossLinksDB::global`] — the singleton the source uses too — searched by mass
difference and preferring the record whose full id names the match's own
`xl_mod`; when nothing matches, the source's `XLMOD:XXXXX` placeholder is
written with the name the hit carries.

### What is still not ported

`ANALYSIS/XLMS/OPXLHelper.h` as a whole: its candidate enumeration, database
digestion, spectrum matching and Percolator plumbing are a search algorithm,
not a file format, and they keep their own ledger entry. Only the six
post-processing functions `readMzIdentMLFile` calls are reproduced here, as
private helpers, because the read path is not meaningful without them. The
divergences that reproduction introduces are in section 3 — the loop-link
position the source overwrites and the first-hit-only merge — and each of them only keeps data the
source drops.

[`SCHEMA_VERSION`]: ../src/format/mzidentml.rs
[`ReadOptions`]: ../src/format/mzidentml.rs
[`ReadOptions::default`]: ../src/format/mzidentml.rs
[`WriteOptions`]: ../src/format/mzidentml.rs
[`WriteOptions::default`]: ../src/format/mzidentml.rs
[`MAX_ITEMS`]: ../src/format/mzidentml.rs
[`MzIdentMLDocument`]: ../src/format/mzidentml.rs
[`IdentificationHit::new`]: ../src/format/mzidentml.rs
[`SpectrumIdentification::add_hit`]: ../src/format/mzidentml.rs
[`Identification::add_spectrum_identification`]: ../src/format/mzidentml.rs
[`load`]: ../src/format/mzidentml.rs
[`load_with_options`]: ../src/format/mzidentml.rs
[`load_with_registry`]: ../src/format/mzidentml.rs
[`load_into`]: ../src/format/mzidentml.rs
[`store`]: ../src/format/mzidentml.rs
[`store_with_options`]: ../src/format/mzidentml.rs
[`store_with_registry`]: ../src/format/mzidentml.rs
[`read_with_registry`]: ../src/format/mzidentml.rs
[`write`]: ../src/format/mzidentml.rs
[`write_with_options`]: ../src/format/mzidentml.rs
[`write_with_registry`]: ../src/format/mzidentml.rs
[`detect_version`]: ../src/format/mzidentml.rs
[`detect_version_from_reader`]: ../src/format/mzidentml.rs
[`MetaInfo`]: ../src/metadata/value.rs
[`CrossLinksDB::global`]: ../src/chemistry/cross_links.rs
[`ControlledVocabulary::psi_ms`]: ../src/format/controlled_vocabulary.rs
[`Error::Parse`]: ../src/error.rs
[`Error::Unsupported`]: ../src/error.rs
[`Error::InvalidValue`]: ../src/error.rs
[`Error::MissingInformation`]: ../src/error.rs

### XML transport limits

The reader accepts UTF-8 bytes with an absent or UTF-8 encoding declaration;
other declared encodings return `Error::Unsupported`, including ASCII-only
documents declaring ISO-8859-1. Entity and structural checks in the reader are
not a general XML/XSD validation service; `is_valid` is.
