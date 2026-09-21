# Native consensusXML interchange

`format::consensusxml` targets the Core SDK `6bfc0e4` handler and consensusXML
1.7. Its optional `consensusxml` feature is enabled by default. It supports owned
streams, transactional read-into, file loading, and atomic plain/gzip/bzip2 file
output. XML text uses UTF-8 output and the shared identification XML decoder.

The model preserves document identity, experiment type, typed map/column/feature
metadata, processing history, protein runs, peptide evidence, portable custom
modification definitions, assigned/unassigned peptide identifications, feature
centroids and handles. XML protein/run references are rebuilt explicitly. A
reserved native metadata entry retains existing run identifiers during native
round trips; source files without it receive new process-local run identifiers.
Path loading sets the loaded-file path and type, which are provenance fields
rather than serialized document content.

RT, m/z and intensity read filters use the source's half-open ranges. Filters do
not affect writing. The source applies centroid attributes only while reading a
handle element; an empty grouped-element list therefore retains a zero centroid.
Invalid nonnumeric unique-ID suffixes become zero, matching source behavior;
numeric overflow is rejected. The source's three supported experiment types are
retained, even though its XSD lists additional historical types.

## Protein-group quantities

The private quantity codec writes typed arrays under the source group-metadata
keys, including an accession-based ownership guard. Reading restores canonical
float, integer and string slots, then arbitrary named arrays in lexical order.
Legacy `abundances` arrays restore the two omitted zero-count arrays. Empty
arrays and all-zero `psm_count`/`distinct_peptides` arrays follow the source's
omission rules; serialized group-array layouts are normalized, not byte-identical
snapshots of arbitrary native vector layouts.

An ownership mismatch fails by default. Setting
`ReadOptions::discard_mismatched_quantities` opts into source-compatible discard;
`read_report` returns the associated diagnostic. Quantity metadata without an
ownership guard is discarded, as in the pinned handler. Regeneration removes
stale owned metadata, including old base entries, so filtered groups cannot be
resurrected. Duplicate array names across value types and the reserved guard name
are checked errors instead of silently overwriting quantities.

## Checked representation boundaries

Orphan peptide/run references, unknown protein accessions, duplicate column IDs,
and malformed input fail transactionally. Map/handle reference inconsistencies
are retained as in the source; `read_report` reports them and
`ConsensusMap::validate_consistency` checks them explicitly. Width fields,
protein modification positions, metadata units, software CV annotations, missing
required timestamps, and nonzero centroids without handles cannot be represented
and are rejected before writing. Scientific sequences and custom definitions use
the shared identification XML representation checks.

All public streams have byte, record, list, payload and work ceilings. Scientific
identification work shares the operation's allowance across runs and features.
Compressed input is decoded as a stream; output is published only after complete
serialization and compression. A caller-owned writer can still contain partial
bytes after an underlying I/O failure.

The adapter's syntax/representation checks are not a generic XSD validator.
Inherited `XMLFile::isValid` is `consensusxml::is_valid`, with the optional
`xml-schema` feature: real XSD validation against the bundled, unchanged
`ConsensusXML_1_7.xsd`, with the class test's verdicts for both fixtures, a
stored map and a stored map with protein-group quantities ported
([XML schema validation](XML_SCHEMA_SUPPORT.md)). The handler implementation
does not certify a TOPP tool or executed C++ differential parity.

## Progress

`ConsensusXMLFile` and its handler derive from `ProgressLogger`, and the file
hands the handler only its log type (`ConsensusXMLFile.cpp:76-78`, `:90-92`).
`consensusxml::load_with_progress` and `store_with_progress` therefore make the
handler's calls on a copy of the caller's logger, made as
`ProgressLogger::clone` makes one; a backend installed with `set_logger`
receives nothing, as a source file's does. Both sections are zero-width, so the
command backend prints one dot per call. Loading: `startProgress(0, 0, "loading
consensusXML file")`, then `setProgress(1)`, `setProgress(2)`, … for the root
and for every `map`, `consensusElement`, `IdentificationRun`, `ProteinHit`,
`PeptideHit` and `dataProcessing` element (`ConsensusXMLHandler.cpp:149`,
`:173`, `:254-256`, `:334`, `:424`, `:485`, `:582`), and `endProgress()` at
`</consensusXML>` (`:130-133`). Storing: `startProgress(0, 0, "storing
consensusXML file")`, `setProgress(1)` … `setProgress(5 + identification runs
+ column headers + consensus features)`, `endProgress()` (`:606-837`), with the
destination opened before the first call as `XMLFile::save_` opens it. These are
the Release build's calls, call for call (tier 1,
`tests/progress_format_readers.rs`). This reader parses the whole document
before converting any of it, so it makes the loading calls once the parse
succeeded: a document that is not well-formed makes none, where the Release
build has made the calls for the elements before the defect (the replay asserts
this against a truncated document), and one the conversion refuses has made
every set and no end. See `docs/PROGRESS_LOGGER_SUPPORT.md#format-readers`.

## Evidence

[The source/fixture manifest](../tests/data/consensusxml_provenance.json) pins the
handler, file wrapper, options, unique-ID rules, source class tests, two unchanged
source files and original schema. [The tests](../tests/consensusxml.rs) cover the
source's six-feature fixture, original range examples, typed metadata/ID/quantity
round trips, compressed transport, ownership and legacy reconstruction, malformed
references, bounds and atomic failures.

The source file loader also accepts ZIP archives. ZIP input is not yet supported
by the native shared path transport and remains an explicit completion gap.
