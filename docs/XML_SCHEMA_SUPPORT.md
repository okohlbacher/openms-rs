# XML schema validation

`format::xml_schema` ports `VALIDATORS/XMLValidator.h` and the `isValid` that
every `Internal::XMLFile` inherits (`XMLFile.h`), for every ported file format
whose source class registers a schema. It sits behind the optional, default-off
`xml-schema` feature and uses the same engine as
[mzML schema validation](MZML_SCHEMA_SUPPORT.md): the exact registry binding
**libxml 0.3.14** over the installed libxml2. `mzml-schema` now implies
`xml-schema`; the mzML path's behaviour, API and tests are unchanged, and its
engine, report types and preflight now live in `xml_schema.rs`.

The schemas are the ones the source ships in `share/OpenMS/SCHEMAS` at the SDK
pin `bc9cc12`, byte for byte, compiled into the crate from
`resources/schemas` ([notices and hashes](../resources/schemas/NOTICE.md)). No
schema is looked up at run time, no schema hint in a document is followed, and
nothing is fetched.

## API mapping

| Source | Rust | Notes |
| --- | --- | --- |
| `XMLValidator()` | none | Stateless: the free functions below take everything per call. |
| `XMLValidator::isValid(filename, schema, os)` | `xml_schema::validate_against(path, schema)`, `validate_against_with_options`, `validate_reader_against` | The caller's schema must be self-contained; see *Native differences*. |
| `XMLValidator::logError_` | `SchemaDiagnostic` | Level, message, filename, line, column, libxml2 domain and code, instead of one formatted line on `os`. |
| `valid_`, `filename_`, `os_` | no counterpart | The verdict is `SchemaValidationReport::is_valid`, the messages are its `diagnostics`. |
| `XMLValidatorErrorHandler_` | no counterpart | Its warning/error/fatal-means-invalid rule is `is_valid`'s. |
| `Internal::XMLFile::isValid(filename, os)` | `xml_schema::validate(kind, path)`, `validate_with_options`, `validate_reader` | `kind` is the schema the class registers. |
| `XMLFile(schema_location, version)` | `SchemaKind::location`, `SchemaKind::version` | The `schema_location_` and `schema_version_` each constructor passes. |
| `XMLFile::getVersion()` | `SchemaKind::version` | Each format module also keeps its own version constant. |
| `FeatureXMLFile::isValid` | `featurexml::is_valid` | `FeatureXML_1_9.xsd` |
| `ConsensusXMLFile::isValid` | `consensusxml::is_valid` | `ConsensusXML_1_7.xsd` |
| `IdXMLFile::isValid` | `idxml::is_valid` | `IdXML_1_5.xsd` |
| `ParamXMLFile::isValid` | `paramxml::is_valid` | `Param_1_8_0.xsd` |
| `TransformationXMLFile::isValid` | `transformation_xml::is_valid` | `TrafoXML_1_1.xsd` |
| `MzDataFile::isValid` | `MzDataFile::is_valid` | `mzData_1_05.xsd` |
| `MzXMLFile::isValid` | `MzXMLFile::is_valid` | `mzXML_idx_3.1.xsd` and the three schemas it includes |
| `MzIdentMLFile::isValid(filename, os, used_version)` | `mzidentml::is_valid`, `is_valid_with_options` | `detectVersion` picks `mzIdentML1.{0,1,2,3}.0.xsd`; `used_version` is `report.schema.version()`. |
| `PepXMLFile::isValid` | `pepxml::is_valid` | `pepXML_v114.xsd` |
| `MzMLFile::isValid` | `mzml::validate_schema` | Unchanged; [MZML_SCHEMA_SUPPORT.md](MZML_SCHEMA_SUPPORT.md). |
| `ImzMLFile::isValid` | `ImzMLFile::is_valid` | Unchanged; `mzml-schema`. |

The per-format entry points are compiled when both `xml-schema` and the
format's own feature are enabled. `xml_schema::validate` needs only
`xml-schema`, so any bundled schema can be used without its reader.

## Preserved source conventions

- **Verdict.** Any warning, error or fatal diagnostic makes a document invalid;
  only libxml2's informational level does not. This is the source's error
  handler, which reports all three (`VALIDATORS/XMLValidator.cpp:31-33`)
  through `logError_`, which sets `valid_ = false`
  (`VALIDATORS/XMLValidator.cpp:107-111`).
- **One registered schema.** `XMLFile::isValid` validates against the schema
  its class registered, whatever the document is (`XMLFile.cpp:393-401`). A
  consensusXML file checked with `SchemaKind::FeatureXML` is a report whose
  `is_valid` is false, not an error. Only mzML chooses between two schemas,
  as `MzMLFile::isValid` does.
- **mzIdentML version.** `mzidentml::is_valid` runs the ported
  `detect_version` first, as `MzIdentMLFile::isValid` does
  (`MzIdentMLFile.cpp:89-100`), and validates against the matching bundled
  schema. The TOPP FileInfo `-v` output reports the same version.
- **mzXML.** The indexed 3.1 schema is used for every mzXML document, indexed or
  not, as `MzXMLFile.cpp:18` registers it. An mzXML 2.1 document is invalid
  against it, as it is in the source.
- **Schema hints are not followed.** A document's `xsi:schemaLocation` or
  `xsi:noNamespaceSchemaLocation` plays no part: the schema is fixed before the
  document is read. The source likewise validates against the grammar it
  preloaded, and `share/OpenMS/SCHEMAS/README.md` notes that the hint URLs in
  OpenMS files were never fetched.

## Native differences

- **Engine.** libxml2's XSD 1.0 implementation, not Xerces-C. Messages, codes
  and columns are libxml2's and vary with its version; libxml2 often reports a
  line without a column. The two engines can disagree on corner cases of the
  standard; the indexed mzML identity-constraint defect CPP-054 is one such
  place, documented with the mzML schema.
- **Malformed XML is an error, not `false`.** A document that is not
  well-formed is `Error::Parse`, where the source reports Xerces' fatal error
  and returns `false`; `TransformationXMLFile_3.trafoXML` and
  `XMLValidator_syntax.xml` are both tested this way. Unsupported encodings, a
  DTD and exceeded limits are errors as well. A well-formed document that
  violates the schema is always a report.
- **A missing file** is `Error::Io` with `NotFound`, where the source throws
  `Exception::FileNotFound`.
- **No schema bound.** A default-constructed source `XMLFile` has no schema and
  `isValid` throws `Exception::NotImplemented` (`XMLFile_test.cpp:45-48`). The
  only kind without a bundled schema is `SchemaKind::External`, and
  `validate(SchemaKind::External, ..)` is `Error::InvalidValue`.
- **Compressed input** is recognised by its bytes and decompressed (gzip,
  bzip2), as every reader in the crate does. The source's `XMLValidator` hands
  the raw file to Xerces.
- **A caller's schema must be self-contained.** `validate_against` refuses a
  schema with `xs:include`, `xs:import`, `xs:redefine` or `xs:override`
  (`Error::Unsupported`), before libxml2 sees it, because resolving one would
  read other files or the network on the caller's behalf; the source leaves
  that to Xerces' default resolver. A DTD in the schema is refused as it is in
  a document, and a schema libxml2 cannot compile is `Error::InvalidValue`,
  where the source reports its errors and returns `false`.
- **Included schemas are composed in memory.** `mzXML_idx_3.1.xsd` and
  `mzIdentML1.0.0.xsd` reach four other bundled files through `xs:include`,
  and libxml2 resolves an include only through its process-wide loaders:
  against a memory buffer the relative location would name a file in the
  working directory. Rather than install a process-wide input callback, the
  engine composes the main schema and its includes into one document, as XSD
  1.0 section 4.2.1 defines an include: each included top-level component is
  copied with the namespace declarations of its own document, and a document
  without a default namespace gets `xmlns=""`, or the including target
  namespace for a chameleon include (`general_types_1.0.xsd`,
  `FuGElightv1.0.0.xsd`). Composition refuses a differing target namespace or
  differing `elementFormDefault`, `attributeFormDefault`, `blockDefault` or
  `finalDefault`, since those govern the local declarations it moves. The
  files in `resources/schemas` stay unchanged.
- **One validation at a time.** libxml2's schema contexts are not safe to use
  from several threads at once (the binding says so for libxml2 2.12 and
  later), so every engine call in the crate, mzML's included, holds one
  process-wide lock. The Rust preflight runs outside it; compilation,
  validation and the collection of diagnostics run inside.
- **Limits.** `SchemaValidationLimits` bound the Rust-side decoding, lexical and
  namespace preflight and the returned diagnostics, for the document and, with
  `validate_against`, the caller's schema on the same budget. They do not bound
  libxml2's own DOM, identity tables or run time. `mzidentml::is_valid` reads
  its version header within the same byte limit.

## Checked boundaries and evidence

Expected verdicts come from the C++ class tests, or from retained C++ output,
never from observing this port. `tests/data/xml_schema_provenance.json` lists
each one with its source line.

| Format | Schema (sha256) | Expectations from | Tests |
| --- | --- | --- | --- |
| XMLValidator | caller's `XMLValidator.xsd` | `XMLValidator_test.cpp:37-53`: valid, missing element, missing attribute, syntax, valid again, missing file | `tests/xml_schema.rs` |
| XMLFile | none | `XMLFile_test.cpp:45-48`: no schema bound | `tests/xml_schema.rs` |
| featureXML | `FeatureXML_1_9.xsd` (9c66e9ab…99a3) | `FeatureXMLFile_test.cpp:385-405`: both fixtures, a stored empty map, a stored loaded map | `tests/xml_schema.rs`, `tests/xml_schema_formats.rs` |
| consensusXML | `ConsensusXML_1_7.xsd` (37e05756…70fb) | `ConsensusXMLFile_test.cpp:274-289` and `331-422`: both fixtures, a stored loaded map, a stored map with protein-group quantities | both |
| idXML | `IdXML_1_5.xsd` (105c9f0c…6379) | `IdXMLFile_test.cpp:268-291`: a stored loaded document with three meta values. The stored empty document is not ported; see *Remaining gaps*. | `tests/xml_schema_formats.rs` |
| paramXML | `Param_1_8_0.xsd` (d672b6b9…941a) | `ParamXMLFile_test.cpp:63-243`: four stored parameter sets | `tests/xml_schema_formats.rs` |
| trafoXML | `TrafoXML_1_1.xsd` (07f16a75…6e8d) | `TransformationXMLFile_test.cpp:38-46`: fixtures 1, 2 and 4 valid, 3 not | both |
| mzData | `mzData_1_05.xsd` (c5c6ad63…4b7f) | `MzDataFile_test.cpp:830-847`: a stored empty and a stored loaded experiment | `tests/xml_schema_formats.rs` |
| mzXML | `mzXML_idx_3.1.xsd` (03ead7c7…a16b) with `mzXML_3.1_mod.xsd`, `separation_technique_1.0.xsd`, `general_types_1.0.xsd` | `MzXMLFile_test.cpp:586-600`: a stored loaded experiment | `tests/xml_schema_formats.rs`, `tests/xml_schema_includes.rs` |
| mzIdentML | `mzIdentML1.1.0.xsd` (8d12337d…f513); 1.0.0 with `FuGElightv1.0.0.xsd`, 1.2.0, 1.3.0 | TOPP_FileInfo_14 and TOPP_FileInfo_15 (test-data `topp/CMakeLists.txt:913-918`, `0cb15f2`): the retained C++ outputs say version 1.1.0, valid and invalid at line 327 | both |
| pepXML | `pepXML_v114.xsd` (64a81531…e2fb) | none: no source test asserts a pepXML verdict, so only the registration (`PepXMLFile.cpp:333`) is checked | `tests/xml_schema_formats.rs` |

The store-then-validate sections validate what the C++ writer stored. Here the
same content is stored by this port's writer and the same assertion made of
it: the class test's contract for a writer, not a comparison with C++ output.

Beyond the ported verdicts:

- `tests/xml_schema_includes.rs` checks the composition against libxml2's own
  include processing of the same files on disk: every verdict and every
  diagnostic, message and line, over eleven mzXML and ten mzIdentML
  documents, most of which violate the grammar in many places. Breaking the
  chameleon rewrite makes the mzXML schema fail to compile, which the suite
  catches.
- `tests/xml_schema_offline.rs` registers a catch-all libxml2 input observer
  before any validation, in its own process, and shows that no validation, the
  composed ones included, loads anything: schema hints pointed at `file:`,
  `http:` and bare include names, stylesheet instructions and XInclude produce
  no load, and a DTD fails before libxml2. A direct include from a memory
  schema, outside the public API, does reach the observer.
- Every bundled schema compiles without a warning, and a root no schema
  declares is a report for each of them.
- The unit tests in `xml_schema.rs` check that composition carries every
  component once and refuses what it would not reproduce, and that the
  verdict, limit and namespace rules the mzML path had still hold.

## Remaining gaps

- **The empty idXML store.** `IdXMLFile_test.cpp:274-277` stores a document
  with no protein run and asserts that it validates. The source writer emits a
  placeholder `SearchParameters` and `IdentificationRun`
  (`IdXMLFile.cpp:186-193`, `IdXMLFile.cpp:416-418`); this port's idXML writer
  refuses the document instead, so there is nothing to validate. That is an
  idXML writer gap, not a validation one.
- **Caller schemas that include or import others** are refused rather than
  resolved. Supporting them needs a resolver that reads only what the caller
  allows, which libxml2 before 2.14 offers only process-wide.
- **Schemas of formats not yet ported** are not bundled: `TraML1.0.0.xsd`,
  `protXML_v6.xsd`, `xQuest_1_0.xsd`, `ToolDescriptor_1_0.xsd`, and
  `Param_1_7_0.xsd`, which the CTD/CWL export uses. Each arrives with its
  format; none is re-scoped out. All of them compile under libxml2 2.9.13.
- **FileInfo `-v`** still refuses (`src/format/file_info/report.rs`,
  `check_flags_supported`). The validation it needs now exists; wiring it,
  with the semantic validation that follows for mzML and mzData, is separate
  work, and TOPP_FileInfo_14/15 would then become a byte-level differential.
- **Diagnostic wording** is libxml2's. Only verdicts and lines are compared
  with the source.

## Feature and CI

`xml-schema = ["dep:libxml", "dep:quick-xml", "file-compression"]`, off by
default; `tools/check_schema_feature_graph.py` still proves the default and
no-feature graphs select no libxml2. The suites need these lines:

```text
cargo test --locked --no-default-features --features xml-schema --test xml_schema --test xml_schema_offline --test xml_schema_includes
cargo test --locked --no-default-features --features "xml-schema featurexml consensusxml idxml paramxml mzml" --test xml_schema_formats
```

The all-features lines already run all four. Platform requirements are those
of [the mzML schema backend](MZML_SCHEMA_SUPPORT.md#optional-dependency-and-platforms).
