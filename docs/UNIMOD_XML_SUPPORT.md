# Unimod XML: UnimodXMLFile and UnimodXMLDataProvider

`FORMAT/UnimodXMLFile.h`, `FORMAT/HANDLERS/UnimodXMLHandler.h` and
`CHEMISTRY/UnimodXMLDataProvider.h` at core `bc9cc12` are ported in
`src/format/unimod_xml.rs`, behind the non-default feature `unimod-xml`
(`quick-xml`, already a dependency of the XML formats). The source provider is
how `ModificationsDB` reads `unimod.xml` and `custom_mods.xml`; the Rust global
registry reads the compiled table `tools/generate_modifications.py` derives from
the same two files, and this reader serves caller-owned registries
(`ModificationsDB::extend_records`) and custom Unimod-style files.

## API mapping

| Source | Rust |
|---|---|
| `UnimodXMLFile()`, `~UnimodXMLFile()` | `UnimodXmlFile` (unit struct) |
| `UnimodXMLFile::load(filename, std::vector<ResidueModification*>&)` | `UnimodXmlFile::load(path) -> UnimodXmlLoad`, `load_from_reader` |
| `Internal::UnimodXMLHandler(mods, filename)`, `onStartElement`, `onEndElement`, `onCharacters` | the private `Handler` driven by quick-xml events |
| the handler's `warning(LOAD, ...)` for an unknown position | `UnimodXmlLoad::warnings` |
| `UnimodXMLDataProvider(const std::string& filename)` | `UnimodXmlDataProvider::new`, `filename()` |
| `UnimodXMLDataProvider::loadModifications()` (with `setFullId()`) | `UnimodXmlDataProvider::load_modifications`; `records_to_modifications` |
| `ResidueModification::SourceClassification`, `setSourceClassification(string)`, `getSourceClassification`, `getSourceClassificationName` | `chemistry::modification_classification::SourceClassification` (`from_text`, `name`, `value`), `ResidueModification::source_classification`, `source_classification_name` |

## Preserved source conventions

- Tags match by qualified name, `umod:X` or `X`.
- One record per `<umod:specificity>`; a multi-character site is origin `X`;
  `Anywhere`, `Any N-term`, `Any C-term`, `Protein N-term`, `Protein C-term`
  map to the term specificities and anything else warns and reads `Anywhere`.
- The delta's masses and element counts are the difference masses and
  formula; digits in an element symbol are its isotope (`13C` is `(13)C`);
  elements outside `<umod:delta>` and `<umod:NeutralLoss>` are skipped.
- A neutral loss without elements is dropped.
- Every site of a modification carries the classification of its **last**
  specificity (the handler's own TODO); the compiled table keeps each
  specificity's own. Text the source does not know (`Synth. pep. protect.
  gp.`) is kept by `classification()` and reads back empty through
  `source_classification_name()`, as `getSourceClassificationName` does.
- Absolute masses stay zero; provenance is `Cv`; the full identifier is derived
  as `setFullId()` derives it (`Oxidation (M)`, `Acetyl (Protein N-term)`).

## Native differences

- `hidden` is read from each specificity; the source has no such field.
- The path is opened as given (no `File::find` against the data path).
- A missing required attribute, an unparsable number, element or record id is
  an error; unbalanced or malformed XML is an error.
- Input is limited to `MAX_UNIMOD_XML_BYTES` (64 MiB) and
  `MAX_UNIMOD_RECORDS` site records.

## Checked boundaries and evidence

Tier 1: `tests/data/unimod_provider_release.tsv` (sha256
`e396582a3ddca5429846d4609aea9a502e48de57081dacba873919f2a6cdbdfd`) is
`loadModifications()` of the executed Release build (`../oracle/A5a-unimod`,
driver sha256 `1295df11d935d1388742df58d42c50b461c9d504c7091486eb12db4558ac7772`,
ibminode06, Linux x86_64 Release bc9cc12/c19e494/174b576, run twice,
identical) on the installed `unimod.xml` (2,791 records) and `custom_mods.xml`
(244), which are byte-identical to `resources/modifications/`.
`tests/unimod_xml.rs` compares every record field by field, checks agreement
with the compiled table (all fields but classification, which differs exactly
as described), and the handler rules on small documents.
