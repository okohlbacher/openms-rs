# Identification path loading and dispatch

The `idxml` feature exposes `format::idxml::{load, load_with_options,
load_with_registry, load_into, store, store_with_options, store_with_registry}`.
These compose the existing bounded idXML adapter with file transport. Returned
`IdXmlDocument` records include the source document identifier and unreferenced
search-parameter blocks; passing caller-owned modification registries keeps the
existing portable-chemistry behavior.

Input uses file magic, independently of suffix. Plain files always work;
gzip/bzip2 require `file-compression`. The idXML feature alone does not imply
compression. Decoded bytes are subject to the supplied `ReadOptions` limits.
`load_into` replaces its destination only after a successful complete load;
custom-options callers can assign the successful owned `load_with_registry`
result in the same way. ZIP remains unsupported.

Direct stores reject recognized non-idXML extensions, including compressed
variants; unknown extensions remain allowed. This source extension check runs
before creating the temporary file.

Source `IdXMLFile::store` opens an ordinary output stream rather than calling
`XMLFile::save_`. These native path writers therefore emit **plain XML regardless
of suffix**, including `.gz`, `.bz2` and `.zip`. They validate and write into a
sibling temporary file, flush/sync it, then rename it over the destination.
Serialization or publication failures preserve existing output and clean up the
temporary file. This is a native atomic-publication improvement. Caller-owned
`write` streams retain their existing possible partial-write behavior on I/O
failure. No compression API is implied by a plain idXML filename.

`FileHandler` additionally supplies compiled-adapter capability predicates and,
when idXML is enabled, identification read/write/load/store dispatch. Only idXML
is available; mzIdentML, OMS, other identification formats and graph persistence
remain explicit missing adapters. Stream dispatch uses an explicit `FileType`.
Path loading gives a recognized filename extension precedence, and otherwise
replays a bounded 64-KiB content preview to the parser. A nonempty input allowlist
must contain the recognized format.

Identification path output follows the source allowlist rule: a recognized
extension determines the type; an unknown extension is accepted only when a
single allowed type selects the adapter. A nonempty allowlist must contain the
selected type. Unlike the native experiment/map writer's optional requested-type
argument, this method accepts a source-style list of allowed output formats.
Unsupported or contradictory choices fail before touching output. FileHandler
returns the full native document, retaining its identifier rather than discarding
it through source vector-only dispatch.

These additions do not alter the idXML record parser/writer, add a C++ fallback,
implement source logger callbacks/schema-validation inheritance, or certify full
IdXMLFile/FileHandler parity. Existing record, metadata and sequence boundaries
are documented in [idXML support](IDXML_SUPPORT.md). Native rank persistence and
other existing interchange conventions are unchanged.

[Path tests](../tests/identification_paths.rs) cover original-source document
loading, registry/default options, magic/suffix disagreement, decoded size limits,
atomic destination preservation, source plain-output suffix rules, input/output
allowlists, content replay, unsupported formats, gzip/bzip2 corruption and the
compression-disabled build. Existing idXML and custom-definition suites validate
the reused serialization and chemistry paths. [Provenance](../tests/data/identification_paths_provenance.json)
records inspected sources at SDK `54a232f`; no C++ execution was used.
