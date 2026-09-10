# Modification registries and OBO providers

`ModificationsDB` stores immutable `Arc<ResidueModification>` records. The shared
global registry contains 3,127 specificity records, in provider order: the
unchanged 3,035-row UniMod/OpenMS TSV followed by 92 XLMOD mono-link records.
[CrossLinksDB](CROSS_LINKS_SUPPORT.md) separately loads 56 XLMOD cross-link
records. The byte-identical pinned ontology, CC-BY-3.0 attribution and hash are
in [resources/modifications](../resources/modifications/README.md).
Historical PSI-MOD is not bundled: its original redistribution grant has not
been confirmed. The reader supports caller-supplied PSI-MOD streams.

## Ownership and identity

Borrowed lookup methods (`find`, `get_modification`, `search_by_mass`,
`best_by_mass`) retain their existing signatures. `entries()` now exposes
`&[Arc<ResidueModification>]`; the corresponding `find_handles`,
`get_modification_handle`, `search_by_mass_handles` and `best_by_mass_handle`
clone shared handles. Sequences, definitions and generator outputs can retain
caller records after the database is dropped, without leaking allocations.
Caller registries may be extended atomically; the global registry is immutable.

`record_id()` returns `Option<u32>`. `unimod_accession()` is independently
optional; `obo_accession()` and `synonyms()` retain vocabulary identity.
`accession()` chooses the UniMod accession, then OBO accession, then an empty
string. Short IDs, full IDs, full names, available accessions and synonyms are
indexed. Only the UniMod prefix is case-insensitive. Ambiguous distinct full IDs
remain checked errors; duplicate full IDs resolve to the earliest stored record,
as in the previous native API. Registry enumeration and mass ties retain provider
and record order. Lookup matches are returned in registry record order.

Create a mutable `ModificationRecord` descriptor with `Default`, fill its public
identity/chemistry fields, and call `ResidueModification::from_record`. This
checks finite masses, uppercase residue origins and control-free names before
freezing the record. An empty `full_id` derives it from name and specificity.
Terminal wildcard `X` is represented as `None`, consistently with the existing
TSV API. `NeutralLoss::new` validates caller-supplied loss masses. Descriptors
need no invented UniMod ID and can carry `absolute_formula` independently of
`diff_formula`. `ModificationsDB::from_records` and `extend_records` accept owned
records; these append complete records rather than interpreting their accessions
as provider aliases.

## OBO loading and source conventions

`ModificationsDB::from_obo(reader, &OboReadOptions)` loads a buffered byte stream;
`extend_obo` returns counts for added records, added alias bindings and unresolved
alias specificity records. Parsing and append checks complete before publication.
A failed read, syntax error or resource limit leaves the previous registry and
all outstanding handles intact. Callers control provider order by appending
streams in order; the source's full provider sequence is UniMod, custom XML,
PSI-MOD, XLMOD. Loading PSI-MOD after the native global snapshot instead appends
its new standalone records after XLMOD.

The implementation follows pinned `OBODataProvider.cpp` and
`ModificationsDB::loadFromProviders_`:

- Terms are emitted by lexical OBO accession, then sorted unique residue origin,
  then N/C terminal expansions. The last term is flushed at EOF.
- PSI `id` remains the short ID; XLMOD uses the term name. Full IDs use the term
  name plus specificity. Synonyms remain separate aliases; PSI-MOD-label is
  retained as its synonym without changing the native short-ID meaning of
  `name()`.
- UniMod references in `def` bind the OBO accession to every already-loaded
  specificity of that UniMod record. The OBO alias's origin/terminus does not
  narrow this mapping, and its other names/synonyms are not installed. Missing
  targets are counted and omitted. Later providers do not retroactively resolve
  them. Standalone terms retain their own chemical fields.
- `DiffMono`, `DiffAvg`, `DiffFormula`, `Formula`, `MassMono`, `MassAvg`, `Origin`,
  `Source`, and `TermSpec` are read, along with XLMOD `monoisotopicMass`,
  `reactionSites` and `specificities`. Quoted `none` properties are skipped.
  Source classification preserves the source's unknown-to-empty fallback.
- Origins B/J/Z are excluded. Unrestricted X records are excluded; terminal X
  requires a nonzero declared monoisotopic delta unless it is a UniMod alias.
  Protein N/C origins expand to protein termini for mono-link loading and peptide
  termini for cross-link loading. Paired XLMOD specificities become their union;
  the registry does not infer which pairs are chemically permitted.
- Mono-link loading excludes exactly `reactionSites=2`; cross-link loading
  excludes exactly `reactionSites=1`. Missing or other counts do not imply
  exclusion. This preserves the actual source filter, including its limits.

Unknown nonchemical fields and ontology relations are ignored, as in the source.
Unknown properties still require the source's quoted-value syntax. This is a
modification provider, not a general ontology graph or lossless OBO writer.
Non-Term stanzas are isolated rather than leaking their fields into a previous
term. Protein N/C `TermSpec` spellings are accepted after whitespace removal,
correcting an unreachable source setter branch. Invalid UTF-8, control bytes,
nonfinite numbers, malformed formulas, missing published names and invalid
residue origins return checked errors. Quoting follows the source's literal
quote splitting, not a general escaped-string OBO grammar.

`Formula` is an absolute **free-residue** composition. The literal empty string
means absent, while nonempty zero/charged compositions remain explicitly present.
Residue chemistry only changes when a delta mass is nonzero or the delta formula
has atoms; absolute-only records with zero differences are source no-ops. With a
change, delta formula takes precedence, then absolute formula, then absolute or
delta mass arithmetic. Terminal chemistry uses declared delta fields. See
[sequence support](SEQUENCE_SUPPORT.md) for checked unavailable composition and
average-mass behavior and [modification support](MODIFICATION_SUPPORT.md) for
notation. Absolute formula is never silently treated as a delta.

## Bounds and validation

Defaults are 16 MiB input, 64 KiB per line (including newline), 100,000 terms,
200,000 records, 1,000,000 alias bindings and 128 MiB estimated registry payload.
Line/input bounds are checked before buffer growth. Expanded record counts and estimated payload
are checked before cloning each expanded record. The alias limit also bounds attempted UniMod-target
visits, including duplicates that would be discarded. Complete existing records,
index entries and new payload are checked before a registry snapshot is cloned.
`max_registry_bytes` conservatively charges strings, formula/map/loss storage and
indices; allocator capacity rounding, bookkeeping and temporary bounded parse
buffers are additional. Lookup uses standard ordered maps and sorted index
vectors. Records remain shared during transactional snapshots. Legacy `from_tsv`
keeps its existing in-memory format/parser; the new streaming limits apply to OBO
and to owned-record registry appends using default limits.

`tests/modification_registry.rs` covers validation, handle lifetimes, resource
limits, errors and atomic extension. Independent
`tests/registry_reference.rs` compares all 148 XLMOD specificity records, exact
mass bits, synonyms and ordering to pinned-source projections and checks
synthetic PSI alias/absolute-formula branches. These tests inspect source behavior
without executing or building C++.
