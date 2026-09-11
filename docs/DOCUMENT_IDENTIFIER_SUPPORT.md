# Document identity and file provenance

`metadata::DocumentIdentifier` covers the source document value and operations at
SDK `82ce5b3`. The owned public fields are `identifier`, `loaded_file_path` and
`loaded_file_type`. Defaults are empty strings and `FileType::Unknown`. Equality
compares **only the identifier**, even when the paths or types differ; cloning,
assignment and `std::mem::swap` preserve all fields. Ordinary field mutation and
clone retain standard Rust costs. There are no accessor wrappers for plain fields.

`set_loaded_file_path` preserves absolute UTF-8 input exactly, including case,
dot segments, repeated separators and trailing separators. A leading slash is
absolute on every platform, matching the source Windows exception. Relative paths
are prefixed with the current directory and use generic slash separators on Windows;
empty input records the current directory. It does not check existence, resolve
symlinks, canonicalize paths or change the stored type. The native limit is 1 MiB
for input and resulting paths; NUL and non-UTF-8 derived paths are checked errors.
A failure preserves all fields.

`set_loaded_file_type` inspects the supplied file's contents independently of its
extension and the stored path. It reads at most 64 KiB of plain or decompressed
bytes, then uses the existing [bounded content heuristic](FILE_HANDLING_SUPPORT.md).
The type is published only after reading succeeds. Empty/unrecognized input becomes
Unknown; inaccessible input leaves the old type unchanged. gzip/bzip2 require
`file-compression`; ZIP returns Unsupported. Reaching the preview bound does not
promise validation of the unread tail or compressed trailer.

The shared detector now recognizes all five previously missing tabular signature
families: mzTab, msInspect TSV, specArray pepList, Kroenik and Percolator PSM output.
Marker order follows their source branches after MS2; existing higher-priority
scientific markers still win. This extends recognition, without adding parsers for
these formats. The native heuristic trims line indentation, so pepList and PSMS
recognize forms outside the exact source first-line predicates. Existing differences
also remain: a 64 KiB preview/512-line IMS bound, stricter PNG magic, finite numeric
DTA checks, UTF-8 lossy sniffing and some marker-order/whitespace differences. This
operation coverage is not an assertion of complete FileHandler parity.

The six direct tests cover source identity/empty-file literals, independent
file-insensitive equality and whole-record copying/swapping, exact absolute paths,
relative/empty paths, atomic bounds and I/O errors, all five source tabular markers
and their precedence, preview truncation and enabled/disabled compression. Source
hashes and the reference boundary are in
[the provenance record](../tests/data/document_identifier_provenance.json).
No C++ execution is claimed for this group.

This standalone value is a prerequisite for complete ExperimentalSettings.
MSExperiment and existing feature-map provenance migration remain separate work;
no duplicate document fields are added to those aggregates here.
