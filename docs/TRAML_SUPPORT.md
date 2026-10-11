# TraML support (`TraMLFile`, `TraMLHandler`)

Ports `FORMAT/TraMLFile.h` and `FORMAT/HANDLERS/TraMLHandler.h` (with their
`.cpp`) at OpenMS4-core `bc9cc12`. Rust: `src/format/traml.rs` (the file
adapter) and `src/format/traml_handler.rs` (the reader and writer), behind the
`traml` feature (`dep:quick-xml`, `file-compression`; in `default`).
`is_semantically_valid` additionally needs `semantic-validation`.

## API mapping

### `TraMLFile`

| C++ member | Rust | Notes |
|---|---|---|
| `TraMLFile()` / `~TraMLFile()` | free functions in `format::traml` | the class holds only the schema location and version |
| `load(filename, exp)` | `traml::load(path) -> Result<TargetedExperiment>` | refuses data the source drops; see below |
| (same, source behaviour) | `traml::load_source(path) -> Result<TraMLLoad>` | drops as the source does and returns its log messages |
| (same, explicit limits) | `traml::load_with_options(path, &ReadOptions)` / `read_with_options(reader, ..)` | |
| `store(filename, exp)` | `traml::store(path, &exp)` / `traml::write(writer, &exp)` | `.gz` / `.bz2` compressed as `XMLFile::save_` does; atomic publication |
| `isSemanticallyValid(filename, errors, warnings)` | `traml::is_semantically_valid(path) -> Result<ValidationReport>` | `report.errors` / `report.warnings` / `report.is_valid()`; mapping `resources/cv/TraML-mapping.xml` (unchanged copy of `share/OpenMS/MAPPING/TraML-mapping.xml`) |
| `XMLFile` schema `"/SCHEMAS/TraML1.0.0.xsd"`, version `"1.0.0"` | `traml::SCHEMA_VERSION` | XSD validation (`isValid`, inherited) is not wired: `format::xml_schema::SchemaKind` has no TraML entry yet (backlog: embed `TraML1.0.0.xsd`) |
| `ProgressLogger` base | not ported | the port reports no progress |

### `Internal::TraMLHandler`

| C++ member | Rust | Notes |
|---|---|---|
| write-only / read-only constructors | `write_document(out, &exp, &db)` / `read_document(input, &options, cv)` | the handler state is private to one call |
| `onStartElement`, `onEndElement`, `onCharacters` | private `Handler::start`, `end`, `characters` | same tag tables and routing |
| `writeTo(os)` | `write_document` | |
| `handleCVParam_` | private `Handler::handle_cv_param` (+ `check_cv_term`) | |
| `handleUserParam_` | private `Handler::handle_user_param` | |
| `writeUserParam_`, `writeUserParams_` | private `Writer::user_params` | |
| `writeCVParams_` (both), `writeCVList_` | private `Writer::cv_params` | |
| `writeTarget_`, `writeRetentionTime_`, `writeProduct_`, `writeConfiguration_` | private `Writer::target`, `retention_time`, `product`, `configuration` | |
| `getChildWithName_` | not ported | declared, never called by the handler |
| `cv_` (psi-ms.obo as `"PI"`) | `cv: Option<&ControlledVocabulary>`; `traml::read_with_options` passes `ControlledVocabulary::psi_ms()` | only `MS:` accessions are checked, as only psi-ms is loaded upstream |
| `XMLHandler::fromXSDString` | `traml_handler::from_xsd_string` | |
| `XMLHandler::writeXMLEscape` | `traml_handler::write_xml_escape` | |
| `XMLHandler::error` / `warning` (LOAD) | `TraMLLoad::warnings` or `Error::Unsupported` | see "Dropped data" |
| `ProductListType`, `ConfigurationListType` | `Vec<TraMLProduct>`, `Vec<Configuration>` | |

## Preserved source conventions

- Element routing by parent and grandparent tag, including the surprising
  ones: a `<Product>`'s `userParam`s land on the **transition**; a
  `<Precursor>`'s `userParam`s land on the owning transition or target, not on
  its precursor list; `IntermediateProduct` `userParam`s are "unhandled".
- TraML modification locations are 1-based on disk and 0-based in memory (`-1`
  N-terminus, sequence length C-terminus); a missing mass delta is `0`.
- A `cvParam` whose accession starts with `UniMod:` (case-insensitive) inside a
  `Modification` sets the UniMod id; on write the name comes from
  `ModificationsDB` for the residue and terminus, first match wins.
- Default transition flags (detecting, not identifying, quantifying) are not
  written; `true`/`1` are the true spellings on read.
- Retention-time unit: a later term without a unit does not overwrite an
  earlier unit; an iRT term without a unit reads as `RTUnit::Unknown`.
- Writer number formats: `os << double` at precision 15
  (`file_info::text_format::ostream_g(x, 15)`) everywhere, except the transition
  precursor m/z, which is `precisionWrapper` (`StringUtils::toStr(double)`,
  `param::value::format_float(x, true)`).
- The two default `cv` entries (MS, UO) are written when the experiment has
  none; the target's own CV terms and userParams are not written ("TODO"
  upstream).
- CV terms with an empty accession are read and written (the C++
  `OpenSwathDecoyGenerator` output carries 33 of them); they are stored through
  the crate-internal `CVTermList::add_source`, and `CVTermList::validate`
  still reports them.

## Native differences

- **Dropped data is refused by default.** Unknown elements, `RetentionTime` /
  `Configuration` / `Target` under an unexpected parent, a `cvParam` or
  `userParam` in an unhandled place, a `Precursor` under neither `Transition`
  nor `Target`, an obsolete CV term whose value fails its type check (the
  source returns early), and an interpretation ordinal/rank outside 0..=255
  (the source wraps into `unsigned char`) are `Error::Unsupported` unless
  `ReadOptions::source_compatible` is set; then the data is dropped (or
  wrapped) and the source message is returned as a warning.
- Compound `MS:1002476` (drift time) sets the **compound's** drift time. The
  source sets the pending *peptide's* drift time, which the next `<Peptide>`
  start discards (C++ defect candidate).
- The value-type checks of `handleCVParam_` sit inside a second
  `if (term.obsolete)` upstream and so only run for obsolete terms; this is
  reproduced, and listed as a C++ defect candidate.
- Product charge (`MS:1000041`) is converted through `toDouble` into an `int`
  upstream; the port truncates the same way and refuses a non-finite or
  out-of-range value.
- `<Sequence>` text is accumulated; Xerces delivers text in chunks and the
  source *assigns* each chunk, so a sequence split by an entity keeps only its
  last chunk upstream.
- Writer: strings the source writes unescaped (cvParam name/value/cvRef,
  peptide sequence, group label, formula, SMILES, UniMod name, protein
  sequence text) are escaped minimally (`&`, `<`, `"`), so every document the
  source writes well formed is written byte-identical and the others become
  well formed. A retention time without a value is refused (the source writes
  a stray `/>` or unit text).
- userParam order: the source follows `MetaInfoRegistry` order, this port the
  metadata key order (established crate divergence, see
  `docs/MZML_WRITER_CPP_PARITY.md`).
- Resource ceilings: `ReadOptions::{max_input_bytes, max_elements,
  max_depth}`; DTDs and external entities are refused.
- Numbers: `StringUtils::toDouble`/`toInt32` are approximated by trimmed Rust
  `parse` (accepts `inf`/`nan`, `+`); `cli::context::to_double` is the exact
  port but `format` may not depend on `cli`.

## Checked boundaries and evidence

`tests/traml.rs` (13 tests; `--features semantic-validation --test traml`
on dax: 13 passed):

- Tier 3: `TraMLFile_test.cpp` sections `store` (three round trips), `equal`,
  `assign`, `add`, `isSemanticallyValid` (empty document).
- Tier 1: every retained C++-written TraML output
  (`ConvertTSVToTraML_{,2,3,5,6}_output`, `TargetedFileConverter_{3,4,13}_output`,
  `FileMerger_9_output`, `OpenSwathDecoyGenerator_output`, test-data `0cb15f2`)
  is loaded and stored again and reproduces the retained file line for line,
  numbers compared as identical doubles (the upstream contract is FuzzyDiff,
  `topp/CMakeLists.txt:59`); only runs of consecutive `userParam` lines are
  sorted. `FileMerger_9` inputs merged plus the tool's `Software` record equal
  the retained output (`FileMerger.cpp:385-406` at topp `174b576`).
- Tier 4: unknown elements refused / reported, Product userParam routing and
  charge truncation, limits, DTD refusal, `from_xsd_string`, escaping, and an
  unset retention time refused with no file left behind.

## C++ defect candidates

1. `TraMLHandler.cpp:1357` sets `actual_peptide_.setDriftTime` for a compound.
2. `TraMLHandler.cpp:1112` and `:1125`: the value-type checks sit in a second
   `if (term.obsolete)`; the obsolete warning is logged twice.
3. `TraMLHandler.cpp:1502` and `:1517`: product charge through `toDouble`.
4. `writeRetentionTime_` (`TraMLHandler.cpp:880`) writes a malformed element
   for an unset value.
5. `writeCVList_` (`TraMLHandler.cpp:1723`) and several peptide/compound
   attributes are not escaped.
