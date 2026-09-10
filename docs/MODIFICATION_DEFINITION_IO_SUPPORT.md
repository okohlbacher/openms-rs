# Portable modification definitions

The native `format::modification_definitions` module ports the public operations of `ModificationDefinitionIO` from OpenMS4-core revision `6bfc0e4711105f4eda2fea86812a83af7c7e791f`. It is available without XML features. The implementation uses owned registries and shared immutable records; it never mutates the global registry.

## Public operations

- `is_definition` selects named records with `ModificationProvenance::Defined`. `Defined`, `Cv`, and `MassOnly` correspond to source provenance values. Provenance is deliberately excluded from chemical equality and ordering. Bundled TSV and OBO providers produce `Cv`; caller-created `ModificationRecord` defaults to `Defined`.
- `collect` gathers definitions referenced by protein search-space names and peptide hits, grouped by run identifier. Bare search-space names collect every matching specificity. Unused runs do not create entries. `collect_iter` accepts a borrowed iterator so feature, subordinate-feature, consensus-feature, and unassigned identifications can share this operation without copying peptide records.
- `collect_feature_map` includes every nested subordinate and unassigned identification. `collect_consensus_map` includes every consensus feature and unassigned identification. Both include the runs' search-space names and have shared-budget variants, without requiring XML features. Traversal charges empty nodes and temporary references as well as scientific definitions.
- `encode` sorts and deduplicates serialized definitions. `attach` unions them with existing `modification_definitions` String metadata. Empty additions do not create a key; an existing empty key remains present. `encode_by_run` retains the source's last-value behavior for repeated run identifiers.
- `register_from` decodes a definition block into a caller-owned `ModificationsDB`. It returns the number of accepted records, including idempotent duplicates, and indexed diagnostics for malformed records. Chemical conflicts abort the operation atomically. `register_search_parameters` treats malformed embedded definitions as errors, appropriate for file interchange that must not silently lose chemistry.
- Shared-budget variants of collection, attachment, run encoding and search-parameter registration let XML adapters carry one cumulative work and payload allowance through the whole file.

`ResidueModification::{to_definition_string,from_definition_string,split_definition_records}` expose the source field codec. `ModificationsDB::{register_definition,has_defined_modification}` support direct typed registration. Direct registration retains complete record chemistry; a record may therefore be valid in a registry while being unrepresentable in the portable version-one format.

## Version-one record

The ten pipe-separated fields are version `1`, short name, full identifier, full name, origin, terminal specificity (`none` for Anywhere), delta formula, explicit delta monoisotopic mass, explicit delta average mass, and comma-separated neutral-loss formulas. Backslash escapes backslash, pipe and semicolon. Records are separated by unescaped semicolons. Field splitting removes an escape before any next character; a final lone backslash remains literal. Empty record-list entries are ignored.

Nine-field historical records are accepted without losses. An empty full identifier is derived from the name/site. Empty numeric fields mean zero. Explicit delta masses remain independent of the delta formula. Neutral-loss masses are reconstructed from their formulas. Finite doubles use shortest round-trip fixed/scientific text with the source exponent sign and minimum two digits; tie lengths prefer fixed notation. Tests assert scientific examples and exact decoded mass bits, without claiming differential C++ runtime verification.

The low-level encoder is explicitly a field projection, like C++. The interchange-level `encode` additionally decodes and compares the result before returning it. Absolute formulas/masses, vocabulary accessions, synonyms, hidden/classification metadata, independent neutral-loss mass overrides, and an atom-empty charged delta formula cannot be silently discarded. Such records return `Unsupported` when the projection changes their value. Named definitions with ordinary default fields, delta formulas or mass-only deltas are supported.

## Checked boundaries

Native records require finite masses, uppercase ASCII origins, and names without control characters. Unknown record versions, extra fields, malformed formulas and unsupported specificity names are errors. This is stricter than source decoding that ignores fields after the optional tenth field. Definition metadata must be a unitless String. Source string conversion of other metadata alternatives is not used.

A full-identifier conflict returns an error rather than the source warning followed by keeping the first registered chemistry. Equal records from different allocations are deduplicated by complete value, replacing pointer-dependent set identity. Diagnostic recovery is available for explicit `register_from` calls; identification-file readers use strict registration. Successful file input returns sequences owning their records after the temporary registry has been dropped.

Definition blocks are bounded to 16 MiB and 100,000 records. Collection and registration also limit work to 50 million units; default registry limits bound aggregate records, aliases and storage. XML variants charge the caller's remaining work and payload before copying definitions or registry indices. Each resolved protein accession and run identifier is charged before copying, so short repeated XML references cannot bypass the payload limit. Additional built-in per-operation caps remain active. Collection charges unmodified residues and absent-name lookups as well as selected modifications. These are conservative logical allocation allowances, not allocator RSS guarantees.

## Evidence

[modification_definition_io.rs](../tests/modification_definition_io.rs) covers the source's all-specificity collection, CV exclusion, deterministic union, duplicate registration count, named-versus-anonymous distinction, and literal `AEADNLDDK(TestIO:FromBlob)K` formula `C54H86N15O28P1` and mass `1423.5504442334`. The source uses a rounded mass assertion; Rust checks that literal within `1e-8` Da. Additional tests independently cover escapes, explicit signed-zero mass preservation, omitted fields, provenance equality, conflicts, budgets, and typed registration.

[idxml_definitions.rs](../tests/idxml_definitions.rs) checks self-contained definitions after a caller registry is dropped, search-only definitions, terminal and residue modifications, anonymous tags, exact second-write stability, unchanged global state, conflict rollback, UTF-16 and Latin-1 decoding, and strict attribute separators. Existing registry/chemistry/idXML tests remain active. Source hashes, locations and the distinction between literal and native cases are recorded in [modification_definition_io_provenance.json](../tests/data/modification_definition_io_provenance.json). No C++ code was built or executed.

The source ProForma Formula/INFO test is represented by the native anonymous mass-tag versus named-record distinction. This increment does not implement or claim ProForma parsing.
