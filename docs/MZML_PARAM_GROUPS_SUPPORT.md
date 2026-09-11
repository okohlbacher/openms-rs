# mzML referenceable parameter groups

The mzML reader supports `referenceableParamGroupList`, `referenceableParamGroup`, and `referenceableParamGroupRef` throughout the parameter contexts in the mzML 1.1 schema. Each definition stores its ordered CV and user parameters. At a reference, the reader applies those parameters through the same handlers used for inline parameters, including unit conversion, duplicate detection, precursor settings and binary-array encoding checks.

This covers run, spectrum, chromatogram, scan list, scan, scan window, selected ion, precursor/product isolation window, activation, binary array, and all schema header parameter contexts. Header references can precede the definitions because `fileDescription` comes before the group list in the schema; they are resolved before entering the run. Empty groups and repeated references to empty groups are supported. Nonempty repeated groups remain subject to ordinary duplicate-field rules. IDs support XML 1.0 Unicode NCNames and schema whitespace normalization. Namespace-prefixed elements and escaped attribute values work identically to inline input.

Group expansion does not change the existing [metadata preservation boundary](MZML_SUPPORT.md#metadata-preservation-boundary). For example, ordinary scan userParams and instrument inventories are still outside the native experiment representation. References in these contexts are resolved and charged against the limits even when their metadata is not retained. Group identity itself is not stored; the writer emits equivalent inline parameters for represented data.

## Validation and bounds

The reader rejects missing or duplicate group IDs, dangling references, misplaced or repeated group lists, count mismatches, invalid IDs, nested group references, non-parameter children and non-whitespace text in parameter elements. The schema allows CVs followed by userParams in a definition, but does not allow a definition to contain references. An empty group list is invalid; an empty group inside a nonempty list is valid. The supported scientific validators remain authoritative after expansion: an unsupported codec, invalid unit, duplicate MS level or binary-array userParam produces the same error as the equivalent inline input.

`ReadOptions` adds three limits:

- `max_param_groups`: 100,000 definitions, including empty or unused definitions.
- `max_total_params`: 10 million group definitions, parameter definitions, inline/expanded parameter applications, and reference occurrences combined. Even references to empty groups consume this budget.
- `max_param_bytes`: 512 MiB of conservative cumulative parameter storage and expansion accounting. It includes map/vector/string overhead and both parsed and retained text. Inline/definition attributes are charged before decoding/copying their keys and values; each reference application is charged before copying retained metadata. Definition storage and ignored-header expansion are included.

The byte accounting is deliberately conservative, not an exact allocator measurement. XML input and binary-array limits remain separate and continue to apply. These limits bound repeated-reference amplification independently of the input XML size. All read errors discard the local experiment; no partially parsed experiment is returned.

The reader remains a supported-subset parser, not an XSD or PSI semantic validator. In particular it does not validate every unrelated header inventory, arbitrary CV placement, all XML ordering constraints, or index checksums.

## Source evidence and independent tests

The source is OpenMS4-core `6bfc0e4711105f4eda2fea86812a83af7c7e791f`. `MzMLHandler.cpp` lines 1119–1127 expands each stored CV parameter in the reference's parent context; lines 1767–1775 store those CV definitions. The pinned handler does not store grouped userParams and uses an empty map entry for unresolved references. The Rust implementation deliberately follows the schema here: grouped userParams work through the inline handler, and unresolved references are errors. The schema also permits forward header references through its document-wide key/keyref constraints; the C++ streaming lookup does not resolve those later.

[The reference suite](../tests/mzml_param_groups.rs) factors an independently encoded existing fixture into groups, asserts literal decoded values, checks precursor settings and auxiliary float/integer/string arrays, covers every schema parameter context, and exercises malformed definitions, conflicts and resource amplification. The schema-valid [source projection](../tests/data/mzml_param_groups_source.mzML) retains the original header, two group definitions and first 15-peak spectrum from `MzMLFile_1.mzML`. The projection changes the declared encoding, retains only that spectrum with an updated list count, omits run userParams, and removes two unsupported binary-array userParams. The original base64 arrays decode independently to m/z 0 through 14 and intensity 15 through 1. This is a projection, not a claim that every field in the full upstream file is supported.

[Provenance](../tests/data/mzml_param_groups_provenance.json) records source/schema/fixture hashes, exact projection steps and source line ranges. The fixture was independently decoded with Python's standard library, and its structure was validated with `xmllint` against the pinned schema. No C++ code was compiled or executed. No new dependency is required. [Numpress](MZML_NUMPRESS_SUPPORT.md) is also supported through these groups. Additional XML encodings and unrepresented acquisition metadata remain open.
