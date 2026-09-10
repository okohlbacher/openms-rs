# File types and shared experiment handling

`format::FileType` represents all 73 formats in Core SDK `6bfc0e4`.
`FileTypeList` preserves source annotation order, property intersections, duplicate
entries, dialog filter strings and exact filter lookup. Names are ASCII case
insensitive, with the source `pqt` alias. Directory types and source-file mzML CV
names retain source values. Enum values intentionally have no C++ numeric ABI.

The registry describes the **source SDK's** capabilities. For example, knowing
that `raw` denotes Thermo data does not provide a Rust RAW reader.
`FileHandler::can_read_experiment` and `can_write_experiment` report the actual
native adapters compiled into the crate.

## Filename behavior

The four lexical helpers do not access the filesystem. Both path separators,
special double extensions, unknown-extension acceptance and repeated compression
suffixes follow the pinned source. The source's case-sensitive `.pep.xml` alias
and unusual stripping results are retained: `sample.pep.xml` strips to
`sample.pep`, and `sample.pep.xml.gz` strips to `sample.pep.xml`.
Compression suffixes are peeled iteratively; an extension match at byte zero is
checked rather than underflowing. `consistent_output_type` returns `Unknown` for
conflicting recognized suffixes and explicit format choices.

## Experiment dispatch

`FileHandler::read_experiment` and `write_experiment` dispatch typed streams.
`load_experiment` uses the filename, then common content signatures when the
suffix is unknown, and enforces an optional allowed-format list. Supported input
formats are DTA, DTA2D, MS2, MGF and optional mzML. Output supports DTA, DTA2D,
MS2, MGF and optional mzML. MS2 writing is a native extension; the source registry
correctly retains the C++ adapter's read-only property. Adapter-specific metadata and scientific limits still
apply. DTA output requires exactly one spectrum and rejects experiment-level
metadata and chromatograms, avoiding the source dispatch's potential data loss.

`store_experiment` stages output in a newly created sibling file, flushes it,
then renames it into place. Validation or serialization failures preserve an
existing destination. I/O failures on a caller-owned writer can leave partial
stream output. Gzip/bzip2 transport uses the `file-compression` feature, enabled
by mzML, featureXML and consensusXML. Both decoders support concatenated members.
The bzip2 dependency uses its pure Rust backend. ZIP containers remain unsupported.

Content recognition reads at most 64 KiB and examines five lines, or up to 512
for IMS markers. It recognizes common XML and text signatures and replays all
preview bytes to the real parser. It is a heuristic, not a validator or a complete
port of every source content-detection rule. Directory probing, other formats,
identification/transition dispatch, source-file bookkeeping and
SHA1 helpers remain open. No adapter is substituted for an unsupported format.

## Feature and consensus dispatch

`read_feature_map`/`write_feature_map` and `read_consensus_map`/`write_consensus_map`
connect the optional featureXML and consensusXML stream adapters. Their
`load_*`/`store_*` variants share content recognition, allowed input types,
compressed transport and atomic publication. Loading records the source path
and file type on the map. Other feature/consensus formats remain explicit errors.
See [featureXML](FEATUREXML_SUPPORT.md) and [consensusXML](CONSENSUSXML_SUPPORT.md).

## Evidence

`tests/file_types.rs` covers the pinned registry count, the source's 50-readable
expectation, filter literals, directory types, filename vectors, compression and
alias edge cases. `tests/file_handler.rs` covers dispatch, unknown suffixes,
allowed types, feature gating, gzip, invalid output preserving existing files,
temporary cleanup and bounded sniffing. Source paths and hashes are recorded in
`tests/data/file_handling_provenance.json`; these are source-derived tests, not
executed differential C++ tests.
