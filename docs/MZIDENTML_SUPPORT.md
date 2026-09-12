# mzIdentML support

Ported headers: `FORMAT/MzIdentMLFile.h` (107 lines) and
`FORMAT/HANDLERS/MzIdentMLHandler.h` (115 lines), with their implementations
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
| `bool isValid(filename, os, used_version)` | **partially ported**: the version-detection half is [`detect_version`]; XSD validation is not ported. | The crate's only schema validator is behind the optional `mzml-schema` feature (libxml) and is mzML-specific. |
| `std::string detectVersion(filename) const` | [`detect_version`], [`detect_version_from_reader`] | Full port, including the 15-line header window, the `version="x.y.z"` preference over the `mzIdentML/x.y` namespace, and the fallback to the adapter default. |
| inherited `Internal::XMLFile::getVersion` / `isValid` | not ported here | `XMLFile` is a separate header. The version this adapter declares is [`SCHEMA_VERSION`]. |
| inherited `ProgressLogger` | **not ported**: the crate has no progress-logger plumbing in the format adapters, and the source only inherits it to hand a reference to the handler, which never logs progress in the mzIdentML path. | |

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
| `logger_` | not ported | See `ProgressLogger` above. |
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
| `writeFragmentAnnotations_` | `write_fragmentation`, `split_annotation` | The regex is replaced by an explicit parser with the same accepted shape. |
| `trimOpenMSfileURI` | `trim_file_uri` | |
| `writePeptideHit` | `write_item`, `write_modifications`, `write_score` | |
| `writeXLMSPeptideHit` | **not ported** | See the cross-linking deferral, section 5. |
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
| `parseSpectrumIdentificationListElements_` | `read_lists` |
| `parseSpectrumIdentificationItemElement_` | `read_item`, `select_score`, `merge_target_decoy` |
| `parseProteinDetectionListElements_`, `parseProteinAmbiguityGroupElement_`, `parseProteinDetectionHypothesisElement_` | `read_protein_detection` |
| `inferModificationLocation_` | `infer_modification_location` |
| `parsePeptideSiblings_` | `read_peptide`, `apply_modification` |
| `findSearchParameters_` | `additional_search_params` |
| `initScoreTermCaches_` | `ScoreTerms::new` |
| `toDoubleOrNaN_`, `toDoubleOrZero_` | `optional_value`, `score_value` |
| `parseSpectrumIdentificationItemSetXLMS`, `buildCvList_` and the other `build*` writers | not ported |

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

---

## 3. Native differences

Each one is documented at the item in the module as well.

| Difference | Source | This port |
|---|---|---|
| **Cross-linking documents** | Read through a separate XL path plus six `OPXLHelper` post-processing steps | Refused with [`Error::Unsupported`] naming `MS:1002494`. Reading them as linear PSMs would split every cross-link into unrelated candidates. Section 5. |
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
count (`max_records`, 1,000,000). The whole document is built in memory and
handed to the writer only on success, so a refused write produces no bytes;
`store` publishes atomically through `path_io::write_plain`. [`MAX_ITEMS`]
(1,000,000) caps the identifications a single write accepts.

### No panics on file-derived data

* No indexing or slicing of file-derived text: `location`, `start`, `end` and
  every other numeric attribute goes through `checked_sub`/`try_from` and is
  compared against the sequence length before use. A single-character attribute
  is read with `chars()`, never `text[0]`, so `pre=""` is an error rather than
  a null character, and a multi-byte character is not split.
* Non-ASCII input is tested: `non_ascii_text_survives_the_round_trip` reads and
  writes an accession containing `日本語`.
* Every arithmetic step on a count is `checked_*` or `saturating_*`.
* `unsafe` is forbidden crate-wide; nothing here needs it.

### Test-to-section mapping

`MzIdentMLFile_test.cpp` has 15 `START_SECTION`s. Nine are ported, six are
unaccounted.

| # | Upstream section | Assertion macros | Status | Rust test |
|---|---|---|---|---|
| 1 | `MzIdentMLFile()` | 1 | ported | `adapter_defaults_match_the_source_constructor` (`SCHEMA_VERSION == "1.3.0"`) |
| 2 | `~MzIdentMLFile()` | 0 | ported | same test; nothing to release |
| 3 | `void load(...)` | 41 | ported | `load_msgf_mini_matches_upstream_literals`, `load_replaces_rather_than_accumulates` (`peptide_ids[0]` score 195 with score type `MS-GF:RawScore`) |
| 4 | `[EXTRA]` modification without `location` (#5443) | 16 | ported | `missing_modification_location_is_inferred_or_skipped` (`PEPTIDEC(Carbamidomethyl)K`) |
| 5 | `void store(...)` | 69 | ported | `store_round_trip_preserves_whole_document`, `store_round_trip_preserves_modified_peptides`, `store_rejects_a_foreign_extension` (`variable_modifications.back() == "Acetyl (N-term)"`) |
| 6 | `[EXTRA] multiple runs` | 5 | ported | `three_runs_round_trip` (`precursor_tolerance == Ppm(20.0)` after the round trip) |
| 7 | `[EXTRA] thresholds` | 9 | ported | `thresholds_round_trip` (`significance_threshold == 0.5`, `pass_threshold == "false"` for both hits of spectrum 17) |
| 8 | `[EXTRA]` regression load of the example files | 0 | ported | `every_upstream_non_crosslinking_fixture_loads` |
| 9 | `[EXTRA] compability issues` | 0 (entirely commented out) | ported | `misplaced_elements_in_a_param_group_are_ignored`, `a_psm_without_a_recognised_score_yields_no_hit`, `a_psm_without_peptide_evidence_still_loads`, `an_identification_without_rt_keeps_no_coordinate`, `evidence_without_positions_keeps_them_unknown` — one test per condition its comments enumerate |
| 10 | `[EXTRA] XLMS data labeled cross-linker` | 40 | **unaccounted** | needs `OPXLHelper`; `crosslinking_documents_are_refused_explicitly` pins the refusal instead |
| 11 | `[EXTRA] XLMS data unlabeled cross-linker` | 40 | **unaccounted** | as above |
| 12 | `[EXTRA]` mzIdentML 1.3 crosslinking scores and thresholds | 8 | **unaccounted** | as above; this section's fixture is the one retained for the refusal test |
| 13 | `[EXTRA]` mzIdentML 1.3 noncovalent association | 4 | **unaccounted** | the fixture declares `MS:1002494` |
| 14 | `[EXTRA]` mzIdentML 1.3 EDC crosslinking | 3 | **unaccounted** | the fixture declares `MS:1002494` |
| 15 | `[EXTRA]` mzIdentML 1.3 multiple spectra per identification | 3 | **unaccounted** | the fixture declares `MS:1002494` (three times) |

Sections 3, 4, 5, 7, 10, 11 and 12 exceed five assertion macros; the four of
them that are ported (3, 4, 5, 7) are ported rather than mapped, with their
literals transcribed into the Rust tests.

Native tests beyond the upstream suite: the reference-resolution table above
(six refusal tests plus `dangling_metadata_references_stay_tolerated`), the
version-detection branches, the parser boundaries (non-ASCII, prefixed
namespace, foreign namespace, every ceiling, empty `PeptideSequence`,
out-of-range substitution location, negative rank), the score-order test, the
C-terminal modification location, the `ProteinDetectionList` rule, the
`Fragmentation` round trip, the three write refusals and the write ceilings.
43 tests in total.

### Evidence tier

Tier 3 (source review) for everything transcribed from the class test and the
four upstream fixtures; tier 4 (independently derived) for the synthetic
documents, the ceilings and the error-variant choices. No C++ was built or
executed and no C++ output was retained, so this is **not** a tier 1
differential. Upgrading it needs an oracle driver under `../oracle/` that links
`libOpenMS` and prints the loaded identifications; the natural first comparison
is `msgf_mini`, whose values are all pinned above.

---

## 5. The cross-linking deferral

`MzIdentMLFile::load` takes a completely different path when any
`AdditionalSearchParams` declares `MS:1002494` ("crosslinking search"):

1. `parseSpectrumIdentificationItemSetXLMS` (643 lines) groups the two to four
   `SpectrumIdentificationItem`s of one cross-link spectrum match by the value
   of `MS:1002511`, and the XL half of `parsePeptideSiblings_` (about 200 more)
   reads donor/acceptor positions, the cross-link mass and the reagent name
   from `MS:1002509`/`MS:1002510`.
2. It then post-processes the result with six `OPXLHelper` functions —
   `addProteinPositionMetaValues`, `addBetaAccessions`, `addXLTargetDecoyMV`,
   `removeBetaPeptideHits`, `computeDeltaScores` and
   `addPercolatorFeatureList` — from `ANALYSIS/XLMS/OPXLHelper.h`, which is
   **not ported**.

`removeBetaPeptideHits` alone changes the hit counts that upstream sections 10
and 11 assert (`peptide_ids2[1].getHits().size() == 1`), so even a complete port
of the two handler branches could not reproduce those sections without
`OPXLHelper`. Reading such a document through the linear path would instead
produce one unrelated PSM per item, with the alpha and beta peptides as separate
candidates: silent corruption of the result rather than a missing feature.

So the reader refuses, with a message naming both the accession and what is
missing, and the six sections are reported as unaccounted. Re-scoping this is a
follow-up in two steps: port `ANALYSIS/XLMS/OPXLHelper.h` (and the
`CrossLinksDB`/`ProteinCrossLink` surface it needs, both of which already exist
in `src/chemistry/`), then the two XL branches and `writeXLMSPeptideHit`.
Nothing else in this module needs to change: the refusal is one call in
[`read_with_registry`], and the library, protocol and result readers are shared.

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
[`ControlledVocabulary::psi_ms`]: ../src/format/controlled_vocabulary.rs
[`Error::Parse`]: ../src/error.rs
[`Error::Unsupported`]: ../src/error.rs
[`Error::InvalidValue`]: ../src/error.rs
[`Error::MissingInformation`]: ../src/error.rs
