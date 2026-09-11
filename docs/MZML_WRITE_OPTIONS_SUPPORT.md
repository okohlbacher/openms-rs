# mzML source writer options and indexed output

This group implements the source mzML writer's `PeakFileOptions` execution for the currently represented native experiment model. It adds `format::mzml::{write_with_peak_options, write_with_peak_options_and_limits, store_with_peak_options, store_with_peak_options_and_limits}`, borrowing the existing options and experiment. Stream/path forms return `PeakWriteReport { binary, indexed, xml_bytes }`; the byte length addresses uncompressed UTF-8 XML. No additional scientific option owner is introduced.

The source pin is [82ce5b373c97f934ffd9b1ffd80215ca66473d0b](https://github.com/okohlbacher/OpenMS4-core/tree/82ce5b373c97f934ffd9b1ffd80215ca66473d0b). [The manifest](../tests/data/mzml_writing_provenance.json) identifies the original handlers, options, class tests, schema, and independent references. This closes the bounded writer-option/index group, not all `MzMLFile` behavior or all source metadata domains. Alternate primary detector roles and independent noise arrays remain separate work; existing precise metadata/header loss guards continue to apply.

## Consumed options and precision

| Existing option | Source default | Behavior |
|---|---:|---|
| `write_index` | true | Indexed wrapper, one or two nonempty index sections, actual SHA-1 footer. |
| `mz_32_bit` | false | Requested ordinary spectrum m/z and chromatogram time precision. |
| `intensity_32_bit` | true | Requested ordinary spectrum/chromatogram intensity precision. |
| `zlib_compression` | false | Binary-array zlib, including Numpress-then-zlib. Independent of filename compression. |
| Mass/time Numpress | None | Coordinate codec; any non-None request also forces both primary ordinary fallback precisions to f64. |
| Intensity Numpress | None | Intensity codec; rejected/empty codec falls back to the source preparation precision. |
| Auxiliary-float Numpress | None | Aligned auxiliary float codec, with ordinary f32 fallback. |
| `force_tpp_compatibility` | false | Suppresses precursor isolation windows and emits even zero charge. Applies to both record kinds. |

For each primary dimension the source preparation rule is `prepared_f32 = requested_f32 && mass_time.compression == None`. Accepted Numpress always declares f64. Thus a mass/time codec request can force f64 intensity output even when intensity Numpress is disabled, and even when the coordinate codec fails its accuracy check. [MzMLHandler.cpp 5614–5658, 5708–5745](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L5614).

Stored coordinates are f64. Ordinary f32 output narrows once and checks finite representability before external output; signed zero survives. Stored intensities and auxiliary floats are f32, so f64 output promotes their existing value without inventing precision. Auxiliary arrays retain their canonical CV identity and type restrictions: a canonical f32-only role falls back to ordinary f32 when Numpress would require f64. Integer and ASCII-string arrays retain existing ordinary/zlib rules, including canonical charge-array i32 and other integer i64 encoding. All three Numpress configurations reuse the existing checked codec and fallback/report contract. See [Numpress support](MZML_NUMPRESS_SUPPORT.md).

The source writer does not apply read filters or sort input. These flags are deliberately ignored, even when their inactive values would be invalid for loading: metadata-only, fill-data, MS levels and all ranges, skip-chromatograms/XML-checks, sort flags, pool size, always-append-data and precursor selection mode. `force_mq_compatibility` and `write_supplemental_data` have no source mzML writer accesses. The potentially large MS-level vector is neither inspected nor cloned. Only fixed-size writing fields are copied.

Existing `write`/`write_with_options` and `store`/`store_with_options` remain ordinary unindexed output with their existing f64-coordinate/f32-intensity defaults. Existing `write_with_numpress` keeps its fixed legacy fallback precisions. Their ordinary bytes remain unchanged except the separately required selected-ion metadata correction below. New defaults use `PeakFileOptions`, hence indexed output. Explicit `write_index=false` permits an empty ordinary document; indexed output requires at least one spectrum or chromatogram. An empty record still receives an index entry.

## Precursor values and TPP

Every writer consumes a unit-free numeric `precursor.cv_terms.metadata["selected ion m/z"]` as the selected-ion CV and excludes it from activation user parameters, following [source 4587–4591 and 4736](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L4587). This is the existing owner used by isolation-target loading, not another metadata map. Consumed strings/lists/units/nonfinite/negative values fail before output.

When that discriminator is present, the effective precursor m/z is written as a separate isolation target even without offsets or `isolation_target_mz`. This preserves target250/selected999 and source target500/selected499.5 through default and isolation-mode reading. The native writer also preserves an explicitly represented zero target in this case. Existing typed `isolation_target_mz` takes precedence when present. See [isolation support](MZML_ISOLATION_SUPPORT.md) and its actual four-writer/two-reader-mode regression.

Explicit TPP mode intentionally removes the isolation window, retains selected m/z, and writes zero charge. Product isolation windows remain independent and are retained. The established native writer already emits one selected ion and its intensity even at zero; that represented shape remains. Activation, possible charge states, mobility, external IDs and supported scalar metadata retain their existing serialization. This is an explicit compatibility-loss option, not reversible isolation transport. [Source 4552–4594](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L4552).

## Indexed bytes, checksum, and paths

Each spectrum/chromatogram offset points exactly to the `<` of its emitted opening element, using UTF-8 bytes rather than characters. IDs are the same escaped IDs written on the records, including existing fallback IDs. Index sections preserve input order and contain every record once. `indexListOffset` points exactly to `<indexList`; the source helper captures before the preceding newline, which its decoder tolerates. This native layout follows the schema's element-offset description rather than claiming whitespace-byte equivalence.

The checksum is real incremental SHA-1 over bytes from byte zero through the final `>` of the opening `<fileChecksum>` tag. It excludes the digest and closing tags. The implementation hashes and counts only bytes actually accepted by the underlying writer, including partial writes. This corrects source placeholder checksum behavior **CPP-049**. Indexed empty output fails before requesting/opening a path or writing external bytes, correcting source invalid empty indexes **CPP-050**. See [the original C++ issue log](../OpenMS_CPP_ISSUES.md). Neither correction claims a security/authentication property.

Path storage uses the existing sibling temporary and atomic publication transport. `.gz` and `.bz2` select outer compression; array zlib remains independent. Indexes address the decompressed XML byte space. The existing plain-file `IndexedMzMLDecoder`/`has_index` remains a plain-byte seeking API and is not promoted to compressed random access. Unknown filename suffix behavior and ZIP rejection follow the existing path transport.

## Checked resource and error contract

`PeakWriteLimits` holds three independent resource domains: its configurable markup/index work and scratch allowance, the existing cumulative `NumpressCoderLimits` for preparation/attempts/fallback/retained binary text, and the existing fixed header allowance. Defaults are 512 MiB XML, 2 billion markup/index work units, and 256 MiB logical markup/index scratch; header and binary limits are not silently reset per record or array.

One header plan and one encoded binary payload are prepared. The writer then runs the same markup logic into a counting sink and finally into the external writer. Both markup traversals, record/metadata visits, escaped/scalar temporary text, offset storage, and checksum byte work are conservatively precharged before external output. Already prepared base64/header text is borrowed during both passes. No complete XML buffer, second parser, seek, or repeated binary encoding is used.

The measurement pass fixes all offsets and exact XML length using a 40-character checksum placeholder. Final emission observes only immutable borrowed data. Output limits use checked arithmetic; configured XML caps above `u64::MAX / 8` reject, keeping the hash's 64-bit bit-length and decoder's signed-offset domain representable. The work default fits a 32-bit `usize`; allocation/index arithmetic remains checked. Caller-controlled stream fragmentation is ordinary external I/O work, not an extra internally allocated document.

Validation, unsupported representation, conversion, configured-resource and binary preparation failures occur before external output; the input experiment/options remain unchanged. An external stream write/flush failure can leave its accepted prefix, as with existing stream APIs. Atomic path publication preserves the old destination on any failure. Ordinary allocator/runtime failures outside Rust's fallible reservation contract are not promoted to recoverable guarantees.

## Evidence and dependency

[The 17 new tests](../tests/mzml_write_options.rs), plus the updated isolation writer regression, cover the source precision matrix and f32 signed zero; source cross-option fallback; all three original Numpress base64 literals with and without zlib; canonical f32 fallback; both record kinds; TPP and selected-ion ownership; exact old ordinary output; ignored read options; complete mixed/escaped/multibyte index identities; XML and cumulative resource limits; finite narrowing; unchanged destinations; path compression; and partial-write I/O.

[The precision fixture generator](../tools/mzml_writing/precision_oracle.py) uses independent Python `struct.pack`/base64 expressions. It does not derive expected values from Rust. [The output checker](../tools/mzml_writing/check_output.py) independently uses Python `hashlib`, raw-byte locations and XML ID decoding, checking a complete ordered bijection between emitted records and index rows. Checksum prefix remainders 55, 56, 63 and 0 modulo 64 and successful writer chunk lengths 1, 7, 63 and 64 are exercised. Actual `xmllint` validates the unchanged pinned indexed XSD; schema success alone is not checksum/offset evidence. Source Numpress literals remain original C++ class-test constants; no C++ build is claimed for this group.

The optional exact `sha1 = 0.10.7` dependency uses `default-features=false` and `force-soft`, activated only with mzML. It is RustCrypto's portable Rust implementation, licensed MIT OR Apache-2.0; its original archives retain both notices. [Versioned API/source](https://docs.rs/sha1/0.10.7/sha1/) and [dependency license reference](MZML_WRITER_DEPENDENCY.md) identify it. No native crypto library, C/assembly build, or runtime network resource is added. Rust 1.98/all-features and Rust 1.85/mzML-only focused tests and strict lint are the validation targets; actual final commands/counts are recorded in the integration handoff. Cross-platform execution beyond the current host is left to CI, not inferred from these tests.
