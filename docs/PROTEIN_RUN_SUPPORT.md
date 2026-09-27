# Protein run metadata and export helpers

The native `ProteinIdentification` covers the source run record member by
member, including the members it inherits from `MetaInfoInterface`, against
Core SDK `bc9cc12`. It extends the existing [identification records](IDENTIFICATION_SUPPORT.md).
Expected values come from `ProteinIdentification_test.cpp` and, where the class
test is silent, from the Release oracle `protein-identification-closure`: a
driver built against the Linux x86_64 Release install on ibminode06, every
case run twice in a fresh process with identical records
([records](../tests/data/protein_identification_oracle.tsv),
[provenance](../tests/data/protein_identification_oracle_provenance.json)).
The header remains **partial**; the open gaps are ranked at the end.

## Equality

`SearchParameters::source_equal` and `ProteinIdentification::source_equal` are
the source `operator==`. The source settings comparison covers the fourteen
typed members and leaves out the settings' own metadata, so two runs that differ
only in `SE:` keys, per-engine settings or `InferenceEngine` are equal although
their `inference_engine()` differs (oracle case `sp_equality`). The run's own
metadata and both run-path lists do count, as they are run metadata in the
source. The derived `PartialEq` stays structural and compares everything; the
round-trip tests rely on that. The source `std::hash` specializations are not
ported: the oracle shows them inconsistent with the source equality (equal
settings, different hashes), and no source tool keys a container by them.

`ProteinGroup::source_less` is the source `operator<` ("best first": higher
probability, then fewer accessions, then the lexically smaller list), including
its NaN behaviour, which falls through to the accessions. `ProteinGroup::sort`
remains the validated sort. The three `get*DataArrayByName` pairs are
`*_data_array_by_name` and `*_data_array_by_name_mut`; a missing name is an
`InvalidValue` error where the source throws `ElementNotFound`.

`PeakMassType::NAMES`, `name` and `all_names` are `NamesOfPeakMassType` and
`getAllNamesOfPeakMassType`. `find_hit_mut` is the mutable `findHit`.

## Charge ranges

`SearchParameters::charge_range` is `getChargeRange` rule for rule and returns
`(min, max)`, checked against 71 strings: the whole string as one `toInt32`
integer; else a comma list of trimmed integers; else one or two `:` sides; else
a split at the first `-` (one minus) or the second (two or three); else `(0, 0)`.
A side drops every `-` and is negative if it had one, otherwise drops every `+`.
Reversed ranges stay reversed and odd strings keep the source value (`1-2-3` is
`(-12, 3)`); `2+ and 3+`, `abc` and an out-of-range number give `(0, 0)`.
Errors mirror the source exceptions: `InvalidValue` for `ConversionError`,
`MissingInformation` for three `:` sides.

Charges are no longer validated. The source never checks the string, and the
former native check rejected settings the source reads: test-data
`topp/SequestAdapter_4_output.idXML` stores `charges="+1+3"` (source `(0, 0)`),
which the native idXML reader refused. `mergeable` and `peptide_ids_mergeable`
validate nothing either, as in the source: an unreadable charge string is
compared as text (oracle case `mergeable_invalid_charges`), and a NaN tolerance
equals nothing. Both keep their `Result`, which is now always `Ok`.

## Coverage and modifications

`compute_coverage` and `compute_modifications` read and check only what the
computation uses: protein sequences and the evidence of the run's own
accessions. Scores, metadata and evidence for other accessions play no part
(oracle cases `coverage` and `validation_free_helpers`: NaN protein and peptide
scores, reversed evidence of a foreign accession). Unknown, reversed or
out-of-protein evidence of an own accession is a `MissingInformation` error, as
the source throws; an endpoint equal to the protein length passes the source
check and writes past its vector, and is an error here. Modifications keep the
native rules for unknown positions (the source would store a wrapped `Size`)
and for positions outside the protein. Both `ConsensusMap` overloads exist as
`compute_coverage_from_consensus_map` and
`compute_modifications_from_consensus_map`: every feature's identifications in
map order, then the unassigned ones on request.

## Primary MS run paths

The source keeps the mzML and raw-file run paths in the run metadata
(`spectra_data`, `spectra_data_raw`); the native record keeps them in
`primary_ms_run_paths` and `raw_ms_run_paths`. `primary_ms_run_path(raw)`,
`nr_primary_ms_run_paths(raw)`, `set_primary_ms_run_path(paths, raw)` and
`add_primary_ms_run_path(paths, raw)` follow the source lists exactly,
duplicates included (oracle case `run_paths`). The source warnings for empty
lists and non-mzML names are not logged. Two source states have no native form:
an absent key (which leaves the caller's output list untouched) and a
non-list `spectra_data` value (which throws `ConversionError`).

`set_primary_ms_run_path_from_experiment(paths, experiment)` is the
filesystem-aware `setPrimaryMSRunPath(const StringList&, MSExperiment&)`
(oracle case `run_paths_experiment`, 16 scenarios). With exactly one source file
path (`path + separator + name`, as `MSExperiment::getPrimaryMSRunPath` builds
it): an mzML name by the `getTypeByFileName` rule (any case, `.gz`/`.bz2`/`.zip`
peeled) that exists on disk becomes the only mzML path and nothing else changes;
a raw name becomes the raw path (existence unchecked) and `paths` the mzML
paths. Otherwise `paths` become the mzML paths. A directory counts as existing;
a `file://` path does not. The name rule is a restriction of
`format::file_types::type_by_file_name`, which this module cannot call because
`format` depends on it; the tests pin the two against each other.

## Engine, inference and export helpers

`inference_engine()` and `inference_engine_version()` return borrowed strings.
Explicit search-parameter metadata wins, including an explicitly empty string.
An existing value of another type errors, preserving the source's strict
DataValue conversion. Without explicit metadata, Fido, BayesianProteinInference,
Epifany and ProteinInference count as inference engines. Percolator counts only
when an indistinguishable group exists. Version falls back to the search-engine
version only when inference data exists.

`original_search_engine_name()` returns the stored engine unless its name contains
`Percolator` or `ConsensusID`. Otherwise it returns the first `SE:` key suffix
whose full key does not contain lowercase `percolator`, or `Unknown` if absent.
The source's "first" is registry order: the order in which the process first
registered each metadata name - for a run read from idXML its document order,
unless an earlier object registered a name first (oracle `engine_global`). The
native map is lexical, so with two or more `SE:` keys the engines can differ
(`SE:XTandem`, `SE:Comet`: source `XTandem`, native `Comet`). This changes tool
output: on two runs that differ only in key order, the Release PeptideIndexer
reads `SE:MSGFPlus` first as MS-GF+, switches to Trypsin/P and matches a peptide
after `RP` that it leaves unmatched with `SE:Comet` first. The native PeptideIndexer
repeats the lexical rule. The same order decides the order of metadata-derived
`search_engine_settings_as_pairs` pairs (oracle `settings_order`).

`search_engine_settings_as_pairs(engine)` emits the twelve standard fields in
source order when the requested engine is empty, or matches a non-Percolator,
non-ConsensusID run. Other requests select literal metadata prefixes and strip
the prefix plus one byte, normally a colon; a boundary inside a UTF-8 scalar
errors. `fill_indistinguishable_groups_with_singletons()` appends one group per
ungrouped accession in hit order. `copy_metadata_only(source)` copies the run
header, score configuration, search parameters, date, both run-path lists and
metadata, and keeps the target's hits and groups. Operations creating owned
output preflight a conservative one-million-descriptor and 64 MiB allowance.

## Inherited `MetaInfoInterface`

Both `ProteinIdentification` and `SearchParameters` inherit `MetaInfoInterface`;
both carry `metadata: MetaInfo`, whose name-keyed operations are the map's own
(see `METADATA/MetaInfoInterface.h`). Three differences are specific to this
record: `getKeys` order (above), the run paths, which live in dedicated lists
and not under `spectra_data`/`spectra_data_raw` in `metadata`, and the
`UInt`-index overloads and `metaRegistry()`, which are not ported.

## Remaining gaps, most consequential first

1. **Metadata key order.** `getOriginalSearchEngineName` and the metadata branch
   of `getSearchEngineSettingsAsPairs` follow `MetaInfo` registry order; the
   native map is lexical. Proven to change Release PeptideIndexer output; also
   reaches ConsensusID `-keep_old_scores` names and the mzTab engine list.
   Needs an order-preserving `MetaInfo` or registry (`METADATA/MetaInfo.h`).
2. **Specificities `no-cterm` and `no-nterm`.** `EnzymeTermSpecificity` has
   four of the six source values; the source exports `no-cterm`/`no-nterm`
   (oracle `names`). Adding them needs the idXML and mzIdentML writers and the
   PeptideIndexer mapping, which match the enum exhaustively.
3. **NaN singleton scores.** `fill_indistinguishable_groups_with_singletons`
   errors on a NaN score; the source appends a group with a NaN probability
   (oracle `singletons`).
4. **Hashes.** The three `std::hash` specializations are not ported.
5. **Typed members.** The date is its serialized string, compared as text
   rather than as a `DateTime`, and the digestion enzyme is its name and
   cleavage expression rather than the full `Protease` record.

[Tests](../tests/protein_identification_oracle.rs) compare every oracle record
listed above; [the earlier tests](../tests/protein_run.rs) cover the branch
tables, the mzTab crash regression, grouping, ownership and resource errors.
