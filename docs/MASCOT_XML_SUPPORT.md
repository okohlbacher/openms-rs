# Mascot XML support

Port of `FORMAT/MascotXMLFile.h` / `.cpp` and the SAX handler it drives,
`FORMAT/HANDLERS/MascotXMLHandler.h` / `.cpp`, at source revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| | |
|---|---|
| Rust | `src/format/mascot_xml.rs` |
| Tests | `tests/mascot_xml.rs` |
| Provenance | `tests/data/mascot_xml_provenance.json` |
| Fixtures | `tests/data/MascotXMLFile_test_{1,2,3}.mascotXML`, `tests/data/MascotXMLFile_test_out_3.idXML` |
| Feature | `idxml` — see "Feature gating" |

Mascot XML is the result format of a Mascot database search: a `<header>`, the
search parameters, the protein/peptide `<hits>`, the `<unassigned>` peptides,
and optionally a `<queries>` section. Retention times are *not* part of the
format in general — they are recovered from each query's `<pep_scan_title>` or
`<StringTitle>` through a lookup the caller supplies, which is why the source
signature takes that helper separately from the file name.

## Feature gating

`mascot_xml` is behind the `idxml` feature. The module needs an XML parser, and
`quick-xml` is pulled in by `idxml`, `mzml`, `paramxml`, `featurexml`,
`consensusxml` and `cv-mapping`; adding a feature of its own is a `Cargo.toml`
change outside this package's scope. `idxml` is the right one to share: Mascot
results exist to be converted to idXML, and `IDFileConverter` is one of the two
direct TOPP consumers of this header.

## API mapping — `MascotXMLFile.h`

| C++ member | Rust counterpart |
|---|---|
| `class MascotXMLFile` | `MascotXmlFile`, a unit struct |
| base `Internal::XMLFile` | **not ported**: for Mascot XML the base's schema location and version are both empty, so it contributes no state. the general `XMLFile` base contract remains a separate unimplemented surface |
| `MascotXMLFile()` | `MascotXmlFile::new` (also `Default`) |
| `void load(const std::string& filename, ProteinIdentification&, PeptideIdentificationList&, const SpectrumMetaDataLookup&)` | `MascotXmlFile::load` and free `load`, returning `MascotXmlResult`; the three out-parameters become its three fields |
| `void load(const std::string& filename, ProteinIdentification&, PeptideIdentificationList&, std::map<std::string, std::vector<AASequence>>&, const SpectrumMetaDataLookup&)` | `MascotXmlFile::load_with_peptides` and free `load_with_peptides`. The `[in,out]` map is taken by shared reference: the source never writes to it |
| `static void initializeLookup(SpectrumMetaDataLookup&, const PeakMap&, const std::string& scan_regex = "")` | `MascotXmlFile::initialize_lookup`, returning the lookup plus the warnings the source logs. A non-empty `scan_regex` is `Error::Unsupported` — see "Native differences" |

Rust-only surface, all native: `read`, `read_with_peptides`,
`read_with_options`, `read_with_registry`, `ReadLimits`, `MascotXmlResult`,
`SpectrumTitleLookup`, `SpectrumMetaData`, `TitleReferenceFormat`, and the
constants `EVALUE_KEY`, `HOMOLOGY_THRESHOLD_KEY`, `IDENTITY_THRESHOLD_KEY`,
`SCORE_TYPE`.

## API mapping — `FORMAT/HANDLERS/MascotXMLHandler.h`

Not one of this package's owned headers, but `MascotXMLFile::load` is a thin
shell around it, so its behaviour is reproduced and its members are listed here.
This does not implement the general handler contract: its `XMLHandler` base, its
`fatalError`/`error`/`warning` reporting contract and its reuse by other
handlers are not ported.

| C++ member | Rust counterpart |
|---|---|
| `MascotXMLHandler(ProteinIdentification&, PeptideIdentificationList&, const std::string&, std::map<...>&, const SpectrumMetaDataLookup&)` | private `Handler::new` |
| `void onStartElement(const char16_t*, const XMLAttributes&)` | private `Handler::start_element` |
| `void onEndElement(const char16_t*)` | private `Handler::end_element` |
| `void onCharacters(const char16_t*, Size)` | the `Event::Text` and `Event::GeneralRef` arms of `Handler::run`, including the `tag_.empty()` guard. Predefined and numeric references are expanded; CDATA remains unsupported — see "Native differences" |
| `static std::vector<std::string> splitModificationBySpecifiedAA(const std::string&)` | private `Handler::split_modification`. Kept private because it needs the modification registry the handler holds; expose it if a caller ever needs it |
| private `protein_identification_`, `id_data_`, `actual_protein_hit_`, `actual_peptide_hit_`, `actual_peptide_evidence_`, `peptide_identification_index_`, `tag_`, `date_`, `date_time_string_`, `actual_query_`, `search_parameters_`, `identifier_`, `actual_title_`, `modified_peptides_`, `tags_open_`, `character_buffer_`, `major_version_`, `minor_version_`, `remove_fixed_mods_`, `lookup_`, `no_rt_error_` | the corresponding `Handler` fields. `actual_title_` is written but never read by the source and is not ported; `minor_version_` is parsed and never used, so the Rust reader only requires `majorVersion` |
| `XMLHandler::fatalError`/`error`/`warning` | `Error::Parse` for the fatal case; `MascotXmlResult::warnings` for the non-fatal ones |

## `SpectrumMetaDataLookup`, in subset

`METADATA/SpectrumLookup.h` and `METADATA/SpectrumMetaDataLookup.h` are separate
unported headers. `SpectrumTitleLookup` implements the part `MascotXMLFile`
reaches:

| C++ | Rust |
|---|---|
| `SpectrumLookup::empty()` | `SpectrumTitleLookup::is_empty` |
| `SpectrumMetaDataLookup::readSpectra(spectra, scan_regexp, get_precursor_rt)` | `SpectrumTitleLookup::read_spectra`, with the source's default `=(?<SCAN>\d+)$` scan expression and without the precursor-RT pass, which `MascotXMLFile` never requests |
| `SpectrumLookup::addReferenceFormat(regexp)` | `SpectrumTitleLookup::add_reference_format(TitleReferenceFormat)` |
| `SpectrumMetaDataLookup::getSpectrumMetaData(ref, meta, flags)` | `SpectrumTitleLookup::spectrum_meta_data(title, want_mz)`; the flag set is exactly the two the handler ever asks for |
| `SpectrumLookup::findByScanNumber` / `findByRT` / `findByNativeID` | private, reached through `spectrum_meta_data` |
| `SpectrumLookup::rt_tolerance` | `SpectrumTitleLookup::rt_tolerance`, default 0.01 s |
| `SpectrumMetaDataLookup::SpectrumMetaData` | `SpectrumMetaData`, with `Option<f64>` where the source uses a quiet NaN |
| everything else on both classes | **not ported** here |

A default-constructed lookup holds no spectra and no reference formats, so it
extracts nothing and leaves every retention time unset. That is exactly what the
upstream class test passes, and why its retained expected output carries no `RT`
attribute.

## Preserved source conventions

Each is covered by a named test in `tests/mascot_xml.rs`.

1. *Character accumulation.* `onCharacters` appends to one shared buffer and
   ignores text while no element is open, so text after a child's end tag is
   dropped and a parent's text is prepended to its first child's. The end
   handler trims the buffer and clears it. Reproduced event for event.
2. *`<NumQueries>` sizes the identification vector*, and `<peptide query="N">`
   indexes it at `N - 1`. Entries no peptide ever reaches stay empty and are
   dropped at the end, which is why this port records the count as an index
   space and materialises an entry only when a query references it; a repeated
   element with a *smaller* count truncates what was read, as `resize` does
   (`a_declared_query_count_is_not_materialised_up_front`).
3. *Score type and search engine* are the literal `Mascot`, stamped on the
   protein identification at each `</protein>` and on a peptide identification
   at its first inserted hit.
4. *The run identifier* is `Mascot_<date> <time>` from `<Date>`, split on `T`
   with any trailing `Z` removed, and is copied onto every peptide
   identification.
5. *`<prot_score>` is read with `toInt32`*, so a fractional protein score is a
   conversion error (`a_fractional_protein_score_is_a_conversion_error`).
6. *Significance thresholds.* `<pep_homol>` sets the threshold; `<pep_ident>`
   records both as hit meta values and replaces the threshold with the identity
   value when the homology value is larger or exactly zero — Matrix Science's
   rule that the homology threshold counts only when it exists and is smaller
   (`the_homology_threshold_only_replaces_a_larger_identity_threshold`).
7. *Fixed modifications are applied to `<pep_seq>`* only when no
   modified-peptide map was supplied. `Carboxymethyl (C)` modifies every C;
   `X (C-term)` and `X (Protein C-term)` the C terminus; `X (C-term R)` the C
   terminus only when the last residue is R; likewise for N. A name with fewer
   than two or more than three space-separated parts is a warning, not an error.
8. *`<pep_var_mod_pos>` is three `.`-separated fields*: the N-terminal slot, one
   digit per residue, the C-terminal slot. A non-`0` digit is the one-based
   index into `variable_modifications` — which is why that list must not be
   expanded before the document ends.
9. *Modification lists come from `<fixed_mods>`/`<variable_mods>` when they
   produced entries and from `<MODS>`/`<IT_MODS>` otherwise*
   (`mods_and_it_mods_are_read_only_without_the_dedicated_sections`). The guard
   is `search_parameters_.fixed_modifications.empty()` — list emptiness, not
   section presence — so an *empty* `<fixed_mods/>` section still lets a later
   `<MODS>` populate the fixed list. An earlier revision of this list said
   "when present"; that is only true for a section that contributed a name.
   Mascot XML 1.x has `<name>` only inside `<variable_mods>`; from 2.1 both
   sections have one, so `majorVersion` plus the enclosing element decides.
10. *Specificity groups are expanded* by `splitModificationBySpecifiedAA`:
    `Phospho (ST)` becomes `Phospho (S)` and `Phospho (T)`, each checked against
    the modification database; terminal specifications and anything that is not
    exactly `name (residues)` pass through. Fixed modifications are expanded as
    they are read; variable ones named in `<variable_mods>` only at
    `</mascot_search_results>`, to keep the index space `<pep_var_mod_pos>`
    refers to intact. The `<IT_MODS>` *fallback* is the exception: it expands
    immediately, as the source does, so "variable groups expand only at the
    document end" holds for the section and not for the fallback.
11. *A `<warning>` naming a modification that "can only be used as a variable
    modification" removes it from the fixed list*
    (`a_warning_element_removes_a_modification_from_the_fixed_list`).
12. *`</peptide>` de-duplicates by sequence.* A hit whose sequence is already
    stored for that query contributes only another protein accession as peptide
    evidence; `</u_peptide>` and `</q_peptide>` always insert.
13. *`<pep_scan_title>` resolution.* Retention time is always requested; the
    precursor m/z only when the identification has none yet. The first
    registered reference format that matches wins; a title matching none is not
    an error.
14. *`<StringTitle>` does two things*: replaces hit sequences from the
    modified-peptide map (every Mascot hit against every supplied sequence,
    quadratic in a always-small hit count, with a warning when the counts
    disagree), then, if the query still has no retention time, reads one from a
    `<m/z>_<RT>` title.
15. *`<RTINSECONDS>` inside `<queries>`* sets the query's retention time
    directly (`a_string_title_supplies_a_retention_time_when_the_query_has_none`).
16. *Post-processing in `load`.* Identifications with no hit at all are dropped
    silently; one with exactly one sequence-less hit is dropped and counted; the
    count of retention-time-less identifications is reported; a non-empty lookup
    that resolved nothing is an error; and a repeated first hit — equal score,
    sequence and charge — is collapsed
    (`identifications_without_a_sequence_are_dropped_with_a_count`,
    `a_repeated_first_hit_is_collapsed`).
17. *`initializeLookup` format set.* With raw data the scan-number and
    DTA-file-name formats are registered as well as the m/z-underscore-RT one;
    without it only the last, because it needs no spectrum to resolve
    (`initialize_lookup_reads_the_spectra_and_registers_the_default_formats`,
    `the_title_lookup_resolves_the_three_default_formats`).
18. *The scan-number matcher needs no backtracking.* The hand-coded equivalent
    of `[Ss]can( [Nn]umber)?s?[=:]? *(?<SCAN>\d+)` commits to the optional
    ` Number` and trailing `s` once it sees them, where a real engine would
    backtrack. That cannot change the outcome: if ` Number` or `s` is present
    and consuming it fails, *not* consuming it requires a digit at the `N` or
    the `s` itself, which is impossible. The groups are disjoint from `\d`, so
    the leftmost match is the same in both engines. The same disjointness holds
    for the DTA form. Source-reviewed against Boost's matcher; the Rust side of
    the table — `Scans followed later by 5`, `Scan Numbers 5`,
    `Scan Numbers =5`, `Scan Numberx5`, `Scan Number scan=7`, `scan=1 scan=2`,
    `xscan=7`, `scan=00012`, `scan=+5` and `scanx=3` — is asserted in
    `the_scan_number_matcher_needs_no_backtracking`. Note that *leftmost* wins
    here, unlike the MGF writer's accession table, where the token iterator
    takes the last match.

## Native differences

| Difference | Why |
|---|---|
| **A caller-supplied `scan_regex` is refused.** `initialize_lookup` accepts `None` (or an empty string) and registers the default formats; a non-empty expression returns `Error::Unsupported`. | The source compiles a Boost regular expression with named groups. This crate has no regular-expression engine and this package may not add a dependency, so the three default formats are hand-coded as `TitleReferenceFormat` variants. A caller needing another form registers a variant directly; extending the enum is the path to a fourth format. |
| **`<NumQueries>` is bounded *and* not materialised.** The source calls `resize` with the converted value, so 2,000,000,000 commits the memory, a negative value wraps to an enormous `size_t`, and a header repeating the element constructs and drops the vector again; `MascotXMLFile::load` then reserves the unfiltered count a second time for its filtered output. | `ReadLimits::max_queries` (5,000,000 by default) is checked before any allocation, a negative count is a parse error, the declared count is only an index space, and the filtered vector reserves the number of survivors. A five-million-query header with four repetitions went from 1.59 s to 0.02 s in release mode on the gate node. Recorded as `MXML-02` and `MXML-09` (`a_hostile_numqueries_is_bounded_rather_than_allocated`, `a_declared_query_count_is_not_materialised_up_front`). |
| **Index bounds are checked.** The source's guard is `peptide_identification_index_ > id_data_.size()`, so `query == size + 1` indexes one past the end and a document with no `<NumQueries>` indexes an empty vector. `<peptide query="0">` is *not* one of those cases: the member is a `UInt`, so `0 - 1` wraps to 4294967295 and the guard does catch it. The genuinely unchecked index is `id_data_[actual_query_ - 1]` in the `<queries>` branches, where `actual_query_` is also unsigned: `<query number="0">` followed by `<StringTitle>` or `<RTINSECONDS>` reads at index 4294967295. | All of them are refused, the `<peptide>` cases with the source's own "show_header=1" message. Recorded as `MXML-01` and `MXML-03` (`a_missing_numqueries_header_is_refused_rather_than_read_out_of_bounds`). |
| **A `<pep_var_mod_pos>` digit beyond the declared list, or a position beyond the sequence, is a parse error.** The source uses `vector::at` for the former, an uncaught `std::out_of_range`, and `AASequence::setModification` for the latter, which throws `IndexOverflow` when the index is not below the peptide length. | Both become `Error::Parse`, which is the same accept/reject boundary — only the exception type differs. An earlier revision of this row called the residue index unchecked; it is checked one level down (`AASequence.cpp`, `setModification`). No upstream fixture reaches either: all 645 aligned pairs in the two large fixtures are in range (`a_variable_modification_index_beyond_the_list_is_refused`). An empty N- or C-terminal slot (`.00.0` split into three fields with an empty first one) is skipped here, where the source reads `temp_string[0]`, gets the string's NUL terminator and fails its conversion. |
| **An empty `<pep_seq>` is handled.** The source's `(C-term X)` / `(N-term X)` branches dereference `end() - 1` and `begin()` without checking. | Recorded as `MXML-07`. |
| **Mascot's `-` flanking marker becomes a terminus marker.** The source stores the character verbatim, and `IdXMLFile` then writes `aa_before="-"`, which is neither of OpenMS's own `[`/`]` markers and which nothing downstream interprets — this crate's `FlankingResidue` cannot represent it at all. `-` before the peptide becomes `NTerminus`, after it `CTerminus`. | It is the information the character carries, and without the mapping `MascotXMLFile_test_2.mascotXML` (which uses `-` 
for both) could not be read. Recorded as `MXML-08`. |
| **The unresolved-retention-time warning is reported.** The source's guard is `if (!id_data_[i].getRT())`, which is false for the NaN it has just assigned, so an unresolved title reports nothing while a title legitimately encoding retention time 0 reports an error. | The port reports the unresolved case, which is the one a caller can act on. Recorded as `MXML-04`. |
| **Non-finite numbers are refused** wherever the source's `toDouble` would accept `inf`/`nan`. | Every consumer of an identification rejects them. An overflowing decimal literal is refused as a *conversion* error rather than becoming an infinity, which is what `std::from_chars` reports and what `toDouble` therefore throws; a second `+` (`++1`) is refused for the same reason. |
| **Predefined entities and numeric character references are expanded; CDATA and DTDs are refused.** Xerces expands references before the C++ handler sees character data. | quick-xml emits references separately, so the reader resolves them before appending to element text. Thus `<pep_score>1&#46;5</pep_score>` reads as 1.5, and `&amp;` in descriptions is supported. Undeclared references fail even in text the handler ignores, and references outside the root are malformed XML. CDATA is an explicit native limitation. Not a C++ defect (`entity_references_are_expanded_and_unresolvable_shapes_refused`, `references_outside_the_root_or_in_ignored_text_are_still_checked`). |
| **A document must be one complete element tree.** Xerces rejects a truncated document, content before the root and a second root; quick-xml's `check_end_names` only pairs the tags it sees. | Exactly one root element, closed, and no non-whitespace character data outside it, checked in the event loop. Without it a download truncated mid-export loaded as a valid partial result and never ran the root-close post-processing (`an_incomplete_or_multi_root_document_is_refused`). |
| **A reference format that matched but whose captured value does not convert is an error, not a miss.** `getSpectrumMetaData` returns after the first matching expression, and the `toInt32`/`toDouble` it calls inside throw, which the handler catches as a warning. | Falling through to the next format invents a successful association from a different part of the title: `500_12 scan=9223372036854775808` would resolve to RT 12 and m/z 500 after the scan-number format had already claimed it (`a_matched_format_whose_value_does_not_convert_does_not_fall_through`). |
| **`^` and `$` in the default formats are line anchors, and scan numbers are 32-bit.** Boost's perl syntax compiles them to `syntax_element_start_line`/`..._end_line` unless `no_mod_m` is set, and `initializeLookup` sets no flags; `SpectrumLookup` converts with `toInt32`. | A wrapped title matches on its later lines and a native ID with a trailing annotation line still yields its scan number; a digit run too long for `toInt32` is the source's `-1`, i.e. no scan-number entry and a warning (`the_title_anchors_are_line_anchors_and_scan_numbers_are_32_bit`). |
| **UTF-8 only.** The source hands the bytes to Xerces, which honours the declaration's encoding. | Every Mascot export in the upstream suite is ASCII; another encoding is an explicit error rather than mangled text. |
| **Namespace prefixes are stripped.** The source matches the Xerces qname, which for a default-namespace document — every Mascot export — is the local name; a prefixed document would match nothing at all. | Stripping is strictly more permissive. |
| **Resource ceilings.** `ReadLimits` bounds bytes, events, depth, queries, hits, per-element text and modification-list length. | The source has none (`resource_ceilings_bound_bytes_events_depth_and_text`). |
| **`ProteinIdentification::date_time` is a string.** The source holds a `DateTime`; this crate stores the serialized form and validates it with `CompletionTime`. | Matches the rest of this crate's identification records. The stored form is `DateTime::get()`'s `YYYY-MM-DD HH:MM:SS`, which differs from idXML's `T` separator. |
| **Warnings are returned, not logged.** | The crate has no global log stream. |
| **No OpenMP.** | The handler is serial upstream too. |

## Checked boundaries and evidence

**Evidence tier 1 (executed differential) for `MascotXMLFile_test_3`.**
`tests/data/MascotXMLFile_test_out_3.idXML` is output the pinned C++ produced:
the upstream section loads `MascotXMLFile_test_3.mascotXML` through
`MascotXMLFile::load`, writes it with `IdXMLFile::store` and fuzzy-compares the
two files with `FuzzyStringComparator::setAcceptableAbsolute(0.0001)`. That
retained file is committed here and compared record by record against this
port's own load of the same input:

- search parameters — database, version, taxonomy, charges, mass type, missed
  cleavages, both tolerances with their ppm flags, the one fixed and five
  variable modifications, and the enzyme (case-insensitively: the idXML writer
  lower-cases it);
- the run date, search engine and search-engine version;
- 7 protein hits, accession, score and empty sequence, in document order;
- 577 peptide identifications: score type, higher-score-better, significance
  threshold, precursor m/z (to 1e-4), and the absence of a retention time;
- 586 peptide hits: score, charge, the rendered modified sequence, every
  peptide evidence's accession and both flanking residues in order, and the
  `EValue` user parameter.

The retained file is read with a line-oriented attribute extractor in the test
rather than with `crate::format::idxml`, because that reader's fixed
50,000,000-unit work budget (`XmlLimits::default`) is exhausted after roughly
forty modified peptides: `AASequence::parse_with_budget` charges about 1.3
million units per annotated sequence, so 586 of them need about 790 million.
`IdXMLFile::store` writes exactly one element per line with plain
double-quoted attributes, so the extractor is exact for this fixture. Raising
that budget belongs to the package that owns `src/format/identification_xml.rs`;
it is listed under "Deferred".

**Evidence tier 3 (source review) for `MascotXMLFile_test_1` and `_2`.** Their
expectations are the class test's literals: for 1.0 the four fixed and three
variable modifications, three identifications at m/z 789.83 / 135.29 / 982.58,
two protein hits `AAN17824` (619) and `GN1736` (293), score type `Mascot`, date
`2006-03-09 11:31:52`, significance threshold 31.8621, hit scores 33.85 / 33.12
/ 43.9, the accession sets, and the three modified sequences; for 2.1 seven
missed cleavages, the modification lists, 1112 identifications with m/z
304.6967 / 314.1815 / 583.7948 at indices 0 / 1 / 1111, 66 protein hits,
`IPI00745872` and `IPI00908876` at score 122, date `2011-06-24 19:34:54`,
significance threshold 5, the five accessions of query 35, hit scores 5.34 /
14.83 / 17.5, and the sequences `VVFIK`, `LASYLDK`, `(Acetyl)AAFESDK`,
`(Acetyl)GALM(Oxidation)NEIQAAK` and `SHY(Phospho)GGSR`.

The class test writes the 1.0 sequences as PSI-MOD accessions —
`LHASGITVTEIPVTATN(MOD:00565)FK(MOD:00445)`,
`MRSLGYVAVISAVATDTDK(MOD:00445)` and `HSK(MOD:00445)LSAK(MOD:00445)`. This
crate's modification table carries UniMod accessions, so each is asserted
residue by residue (name and position) and then through the UniMod rendering
`LHASGITVTEIPVTATN(UniMod:7)FK(UniMod:52)` and so on. UniMod 7 is deamidation
and UniMod 52 guanidination, the same two modifications MOD:00565 and MOD:00445
name, so the two vocabularies agree on the result.

All 4 upstream `START_SECTION`s are ported — none is merely mapped.

**Independently derived (tier 4).** The resource ceilings, the bound checks, the
malformed and non-ASCII documents, the entity-reference expansion, CDATA and DTD refusal,
the one-complete-document rule, the lazily materialised query vector, the
`-` flanking mapping, the title-lookup format tests (including the line
anchors, the 32-bit scan width and the matched-but-unconvertible case), the
fractional-protein-score rejection and the threshold table are Rust-only checks
derived from reading the implementation.

**Second-model review.** An adversarial review by another model found that the
catch-all arm of the event loop *deleted* every XML entity reference from
element text, so `<pep_score>1&#46;5</pep_score>` was read as the score 15 —
silent corruption of a scientific value. Predefined and numeric references
are now expanded; undeclared references, CDATA and DTDs remain refused. The same review supplied the truncated-document,
NumQueries-amplification, matched-but-unconvertible, line-anchor and 32-bit
scan-width findings above, and corrected three claims in this document: the
`<MODS>`/`<IT_MODS>` guard tests list emptiness rather than section presence,
the `<IT_MODS>` fallback expands specificity groups immediately, and the
residue index of `<pep_var_mod_pos>` is checked one level down by
`AASequence::setModification` rather than being unchecked. Each behavioural
finding has a named regression test.

**Tolerances.** 1e-4 absolute wherever the upstream comparison uses
`setAcceptableAbsolute(0.0001)` or `TOLERANCE_ABSOLUTE(0.0001)`; 1e-9 absolute
for the scores and thresholds the class test compares with the tighter
`TOLERANCE_ABSOLUTE(0.00001)` default; exact for counts, charges, accessions,
flanking residues, sequences and every string.

## Deferred

- `src/format/identification_xml.rs`'s fixed 50,000,000-unit work budget cannot
  read an idXML with more than roughly forty modified peptides, so
  `idxml::load("MascotXMLFile_test_out_3.idXML")` fails with "peptide parent
  parsing resource limit exceeded". The per-sequence charge in
  `AASequence::parse_with_budget` scales with the whole modification database
  for every annotation. That file is outside this package's scope; until it is
  raised, the differential uses its own extractor.
- `METADATA/SpectrumLookup.h`, `METADATA/SpectrumMetaDataLookup.h` and
  `METADATA/SpectrumNativeIDParser.h` stay unported; only the subset above is
  reproduced, and a caller-supplied `scan_regex` cannot be honoured until a
  regular-expression facility exists in the crate.
- `FORMAT/XMLFile.h` and `FORMAT/HANDLERS/XMLHandler.h` stay unported; the
  handler's error-reporting contract is replaced by `Result` plus a warning
  list.
- No writer: Mascot XML is a search-engine output and the source has no `store`.
- `FORMAT/PepXMLFileMascot.h` is a separate header and is not touched here,
  although it produces the modified-peptide map `load_with_peptides` consumes.
- `ProgressLogger` is not involved in the source class either.
