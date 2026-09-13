# Native idXML interchange

`format::idxml` reads and writes identification results using the existing Rust identification types. It targets the idXML 1.5 representation, including portable modification definitions, in OpenMS4-core revision `6bfc0e4711105f4eda2fea86812a83af7c7e791f`. Original source fixtures remain pinned to revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. The optional `idxml` feature is enabled by default and shares the `quick-xml` dependency with mzML. It calls no C++ code.

```rust
use openms::format::idxml;
use std::fs::File;
use std::io::{BufReader, BufWriter};

# fn example() -> openms::Result<()> {
let document = idxml::read(BufReader::new(File::open("results.idXML")?))?;
println!("{} peptide identifications", document.peptide_identifications.len());
idxml::write(BufWriter::new(File::create("copy.idXML")?), &document)?;
# Ok(())
# }
```

`IdXmlDocument` contains `document_id`, `protein_identifications`, `peptide_identifications`, and `unreferenced_search_parameters`. Native run identifiers link the two identification collections. Unreferenced search parameter blocks remain available instead of being discarded. XML `SP_*` and `PH_*` IDs are document-local transport references; output generates fresh IDs.

## Supported information

| idXML information | Native representation |
|---|---|
| IdentificationRun date, search engine/version, search parameter reference | `ProteinIdentification` run fields and owned `SearchParameters` |
| Protein/peptide score type, score direction and significance threshold | Corresponding identification fields |
| Protein accession, sequence, score and coverage | `ProteinHit`; absent coverage and the source `-1` sentinel become `None` |
| Protein group and indistinguishable group encodings | `ProteinGroup` probability and accession lists |
| Peptide sequence, charge and score | `PeptideHit`; sequence parsing preserves B/Z/X, named modifications and numeric mass annotations; numeric terminal tags use leading placement or explicit terminal markers |
| Peptide RT, MZ and spectrum reference | Optional coordinates and `spectrum_reference` metadata |
| Protein references, flanking residues, start/end lists | Ordered `PeptideEvidence` records; zero-based inclusive positions, `-1` becomes `None` |
| Database/version, taxonomy, charges, mass type, enzyme, missed cleavages and tolerances | `SearchParameters`, including absolute/ppm tolerance units |
| Fixed and variable modification names | Ordered search parameter name lists |
| `EnzymeTermSpecificity` UserParam | Typed unknown/full/semi/none field |
| `spectra_data` and `spectra_data_raw` UserParams | Typed primary/raw MS run path lists |
| `fragment_annotation` UserParam | Typed peak annotations, including quoted separators and escaped quotes/backslashes |
| `_ar_<index>_*` UserParams | Typed analysis results, main score, score direction and named sub-scores |
| UserParam `string`, `int`, `float`, `stringList`, `intList`, `floatList` | `MetaValue` alternatives; integers use native i64, floats must be finite |

Metadata is supported on search parameters, protein identifications/hits, and peptide identifications/hits. The source bracketed list syntax is retained. Within string-list entries, the source `\|` escape represents a comma. Literal XML attribute whitespace follows XML 1.0 normalization; explicit character references such as `&#10;` retain the referenced value. Text, attribute names stored as values, and references are escaped when writing.

Evidence columns follow the source loader: the longest supplied list determines the evidence count, and shorter optional lists fill a prefix. Repeated protein accessions remain repeated evidence records. The writer retains the correspondence between each accession, flank and position. Empty accessions are allowed only after all nonempty accessions because XML IDREFS cannot encode an empty slot in the middle; an unrepresentable ordering fails before output. Unknown-only evidence records are preserved by writing their flanking markers.

Calendar dates use `YYYY-MM-DDTHH:MM:SS`, optionally followed by fractional seconds and `Z` or an explicit `±HH:MM` offset. Dates are retained without timezone conversion. Each run requires a date on writing. Empty runs and peptide identifications with no hits are retained.

## Differences chosen to preserve data

Reading an upstream file generates deterministic run identifiers from engine, date and run index, rather than using the C++ process-global unique-ID generator. These identifiers are unique within the document, not across independently loaded documents. Callers combining documents must assign distinct identifiers and update linked peptide IDs. Native identifiers round-trip through the reserved `openms-rust:run_identifier` UserParam on ProteinIdentification. Nonzero native hit ranks, which have no schema attribute, round-trip through `openms-rust:rank`. These are ordinary schema-valid UserParams; other implementations may retain them as generic metadata.

The writer retains hit order and empty peptide identifications. C++ sorts peptide hits and omits empty identifications. Output groups peptide identifications by protein-run order, as required by the XML hierarchy. Search parameter blocks are written separately for each run, avoiding C++ parameter deduplication that can collapse distinct metadata. Unused parameter blocks follow these run parameter blocks. Input shared-parameter references become equivalent owned values; original XML reference names and sharing are not retained.

Enzyme names retain their input spelling, including unknown names; the adapter does not silently substitute another enzyme. Converted special UserParams move to their typed fields rather than remaining duplicate generic metadata. Reserved metadata collisions fail explicitly. Peak annotation order and analysis score direction are preserved; the pinned C++ writer sorts peak annotations, and its loader omits the encoded analysis `higher_is_better` flag. Protein group references must resolve within the same run, and group/analysis indices must be contiguous from zero; analysis indices use canonical unsigned decimal notation. Unknown references and malformed lists are errors.

## Explicit limits

The adapter rejects constructs or native state it cannot represent faithfully:

- More than one ProteinIdentification block within one IdentificationRun. The native flat run model and the pinned loader do not preserve this schema-permitted arrangement correctly.
- Metadata attached to a FixedModification or VariableModification element; the native search model stores names only.
- Named definitions whose version-one projection cannot preserve complete chemistry or metadata, conflicting same-name definitions, or peptide modifications unsupported by the native sequence parser. No custom definition is replaced by a similarly named built-in record.
- Typed peptide modifications whose displayed sequence cannot be parsed back into the same attachment and chemistry. Some source modified-peptide generator paths put terminal records on residue slots or mismatched termini. Writing rejects those states before output instead of silently relocating or losing the annotation.
- ProteinHit modification positions, protein-group sample data arrays, and custom `digestion_regex` values. These lack an implemented idXML encoding.
- Metadata units and `MetaValueData::Empty`, which the schema cannot distinguish from an ordinary empty string without an extension. Empty strings and empty lists are supported separately.
- A literal `\|` inside a string-list entry, or a list containing exactly one empty string. The source encoding cannot distinguish these values from an escaped comma or empty list.
- Duplicate protein accessions within one run when writing, unmatched peptide run identifiers, missing evidence/group protein accessions, and evidence positions outside the schema's signed 32-bit range.
- Unknown elements/attributes, foreign element namespaces, invalid boolean/numeric/text values, malformed references, and unsupported reserved encodings. No unknown subtree is silently skipped.

The reader accepts UTF-8, ASCII, Latin-1 and UTF-16 XML 1.0 (with checked byte-order marks/declarations) and compatible idXML version declarations from 1.0 through 1.5; output is always 1.5. Historical files still need to use supported fields and sequence syntax. DTDs and external entities are not accepted. The shared XML tree supports CDATA and text geometry for map formats; idXML structural checks reject non-whitespace element text. Character references and the five predefined entities are resolved by quick-xml's `BytesRef::resolve_char_ref` and `escape::resolve_xml_entity`, so element text and attribute values follow one grammar: XML 1.0's `CharRef`, with a lowercase `x` for hexadecimal, no sign and no NUL. Until 2026-09-13 the shared reader resolved text references itself and read a signed number such as `&#+46;` or `&#x+2E;` as the character it names. Those are now `Error::Parse`. `&#X2E;`, `&#0;` and `&#x0;` were refused before and still are. Any other entity name is `Error::Unsupported`, because it would need the DTD this reader refuses. featureXML, consensusXML and transformation XML read through the same tree, so the change applies to them too. The legacy `userParam` spelling is accepted as `UserParam`. XML declarations enforce version/encoding/standalone field order, uniqueness, separators and allowed values. Processing instructions require a valid XML Name and cannot use the reserved case-insensitive `xml` target. Valid comments and processing instructions are ignored; stylesheet/schema URLs are never fetched. The parser enforces the represented structure and attributes, but it is not a general XSD validator. Run date validation covers positive four-digit Gregorian years, ordinary seconds `00..59`, and timezone offsets through 14 hours; leap seconds, year zero, and expanded/negative years are outside this subset.

## Bounds and failure behavior

`read_with_options` defaults to 64 MiB of XML, one million XML elements, and one million entries in any metadata, protein-group accession, fragment-annotation or evidence list. The shared tree permits up to 260 XML levels, while idXML structural validation still enforces its own schema layout. A shared 50-million work allowance and 256-MiB conservative payload allowance cover decoded input, tree storage, typed conversion, registry indices and sequence reconstruction. Resource limits can bind before the raw byte limit. The byte limit includes attributes, comments and processing instructions. Parsing creates a bounded XML tree and then owned native records; this is not a streaming API and uses multiple representations of some text during conversion. Returned records appear only after successful parsing and validation. Parse errors use line zero when no line number is available.

`write_with_options` defaults to 64 MiB of serialized output and one million output elements. It validates the document and stages its XML tree and complete byte output before touching the destination. A shared work/payload preflight also bounds native metadata, sequences, definitions and registry copying before serialization. These conservative logical allowances are not a bound on the Rust allocator's complete memory usage. Invalid or unsupported records, references, metadata, or output limits leave the supplied writer untouched. The final write and flush propagate I/O errors; an external writer can already contain partial output after such an I/O error. In particular, constructing a `File::create` writer truncates that file before this function receives it; callers needing filesystem replacement atomicity should stage a separate file themselves.

## Verification and provenance

`tests/idxml.rs` checks the exact pinned `IdXMLFile_whole.idXML` and `IdXMLFile_no_proteinhits.idXML` fixtures against source-test expectations, complete native round trips, ordered and repeated evidence, modified peptide sequences, typed metadata, annotations, independent run linking, calendar dates, XML normalization, character references in element text, malformed input, reserved-state errors, limits and writer error behavior. Fixture bytes and the schema are unchanged upstream copies. [idxml_provenance.json](../tests/data/idxml_provenance.json) records their source paths and SHA-256 hashes.

The [modification-generation workflow](../tests/modification_generation_workflow.rs)
also checks exact normal variant/search-definition round trips and verifies
zero writer calls when source terminal placements cannot be represented.

When `xmllint` is available, the Rust test independently validates generated XML against the exact copied `IdXML_1_5.xsd`. It passed on this development machine. The test explicitly reports a skip if the external validator is absent; ordinary parser/round-trip tests require only Rust. Source fixtures can also be checked without network access:

```sh
xmllint --nonet --noout --schema tests/data/IdXML_1_5.xsd tests/data/idxml_upstream_whole.idXML
cargo test --offline --test idxml
cargo test --offline --no-default-features --features idxml --test idxml
```

Reference: [IdXMLFile.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/FORMAT/IdXMLFile.cpp), [XMLHandler.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/FORMAT/HANDLERS/XMLHandler.cpp), [IdXMLFile_test.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/IdXMLFile_test.cpp), and [IdXML 1.5 schema](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/share/OpenMS/SCHEMAS/IdXML_1_5.xsd). Original source and fixture licensing is BSD-3-Clause; attribution is retained in the source headers and provenance manifest. No C++ build was performed.

## Caller-owned modification registries

`read_with_registry(reader, options, registry)` resolves sequence annotations
against an explicit `ModificationsDB`. Returned sequences own shared records and
survive dropping the registry. `write_with_registry(writer, document, options,
registry)` verifies exact typed reconstruction against that registry before any
writer call. Same-name records with different chemistry are errors. Writers collect named `Defined` modifications from search-space entries and all peptide attachments, then embed their portable records in each run's `modification_definitions` search metadata. Readers register these definitions in an owned copy of the supplied or global database before parsing peptide sequences. Portable custom records therefore round-trip through the default reader without an external registry. Input records and the immutable global database remain unchanged. Anonymous tags retain their normal bracket syntax.

The complete source codec, checked omissions and collection/registration APIs are described in [MODIFICATION_DEFINITION_IO_SUPPORT.md](MODIFICATION_DEFINITION_IO_SUPPORT.md). The shared private XML helpers also serve featureXML and consensusXML, while their run/group dialects remain separate. Early map-header modes stop consuming the byte stream at the requested root-child opening tag, including for UTF-16 input; ordinary idXML reads validate the entire input.

## File paths and dispatch

[Native path APIs](IDENTIFICATION_PATH_SUPPORT.md) compose this adapter with
bounded plain/gzip/bzip2 input and atomic plain output. FileHandler now dispatches
identification documents to idXML with source extension/allowlist rules.
