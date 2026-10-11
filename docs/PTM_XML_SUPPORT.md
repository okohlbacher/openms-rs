# PTMXML

`format::ptm_xml` (feature `kb-xml` or `idxml`) ports `FORMAT/PTMXMLFile.h` and
`FORMAT/HANDLERS/PTMXMLHandler.h` from OpenMS4-core
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Parsing runs on the shared bounded
SAX driver `src/format/sax_xml.rs` (quick-xml), which stands in for the
Xerces-based `Internal::XMLFile` / `XMLHandler` base classes and delivers
character data where the Release build's Xerces does. The Inspect and Sequest
adapters (`format::inspect_infile`, `format::sequest_infile`, `idxml`) read
their modification tables through it; lane F3b's separate
`format::ptm_xml_file` was folded into this module at wave-4 integration.

## API mapping

| Source member | Rust |
| --- | --- |
| `PTMXMLFile()` | none needed: free functions |
| `load(filename, ptm_informations)` | `load(path) -> PtmInformation`, `load_with_limits`, `read` / `read_with_limits` over any `BufRead` |
| `store(filename, ptm_informations) const` | `store(path, &PtmInformation)`, `write` over any `Write` |
| `std::map<std::string, std::pair<std::string, std::string>>` | `PtmInformation = BTreeMap<String, (String, String)>` (composition, possible amino acids) |
| `PTMXMLHandler(map, filename)`, `~PTMXMLHandler` | private `Handler`, one field per source member (`tag_`, `name_`, `composition_`, `open_tag_`) |
| `PTMXMLHandler::onStartElement` / `onEndElement` / `onCharacters` | `Handler` methods driven by `sax_xml` |
| `PTMXMLHandler::writeTo` | `write` |
| inherited `XMLFile::parse_` compression detection, `save_` suffix compression | gzip/bzip2 by leading bytes on read; `.gz` / `.bz2` suffix on `store` |
| — (native) | `ReadLimits` |

## Preserved source conventions

- The handler's state machine, including its quirks: text counts only
  directly inside `<name>`, `<composition>` and `<possible_amino_acids>`,
  nesting is not checked, `name_` and `composition_` carry over into later
  `<PTM>`s that lack them, an empty `<possible_amino_acids/>` stores nothing,
  a later entry replaces an earlier one of the same name, and every text is
  trimmed of space, tab, CR and LF (`StringUtils::trimmed`).
- Character data split by CDATA, a comment or a processing instruction
  arrives in pieces, as Xerces delivers it, and each piece overwrites the last.
- `writeTo` layout: no XML declaration, tab indentation, `\n` line ends,
  entries in map key order.

## Native differences

- `load` returns a new map, and only on success; the source clears the
  caller's map first, so a failed load leaves it empty.
- The writer escapes `& < > " '` (the source's `writeTo` does not, so a name
  with markup would not round-trip there) and refuses characters XML 1.0
  cannot carry.
- `ReadLimits` bound the input bytes, depth, events, text run, entry count and
  total result bytes (the carried-over composition could otherwise amplify a
  small document); exceeding one is `Error::InvalidValue`.
- A document type declaration or a non-UTF-8 encoding is
  `Error::Unsupported`; malformed XML is `Error::Parse` (the source's
  `ParseError`); a missing file is `Error::Io` (`FileNotFound`).

## Checked boundaries and evidence

`tests/ptm_xml.rs` (6 tests, `--features kb-xml`): the `load` and `store`
sections of `PTMXMLFile_test.cpp` over the byte-identical fixture
`tests/data/PTMXMLFile_PTMs.xml` (tier 3); multi-entry writer order from
`PTMXMLHandler.cpp:28-35`, escaping, compressed store/load and failure
atomicity (tier 4). Two tests compare against the executed Linux x86_64
Release build (`../oracle/integ-F-ptmxml`, ibminode06,
`openms4-release-bc9cc12-c19e494-174b576`, run twice; records in
`tests/data/ptm_xml_release/cpp_release_records.txt`, tier 1): the source's
`store` writes a name holding `<`, `>`, `&`, `"` or `'` unescaped, which its own
`load` then refuses, and reading the escaped name keeps only the last character
run, because Xerces splits the text at every reference and
`PTMXMLHandler::onCharacters` assigns (`PTMXMLHandler.cpp:56-66`); `A&amp;B`
reads as `&B`. This writer escapes (the documented native difference), so a
name holding one of those characters does not survive a round trip here either.
`tests/search_engine_ptm_infiles.rs` loads and stores the same fixture through
this module under `idxml` alone. Provenance: `tests/data/ptm_xml_provenance.json`.
