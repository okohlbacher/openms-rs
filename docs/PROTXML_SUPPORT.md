# ProtXMLFile

`FORMAT/ProtXMLFile.h` (core `bc9cc12`), ported by lane F3a in wave 4 as
`src/format/protxml.rs` (feature `idxml`) and salvaged from its worktree at
wave-4 integration, where the port was uncommitted; this document was written
at integration. It reads ProteinProphet protXML into a protein identification
(hits, protein groups, indistinguishable groups) and a peptide identification,
for `IDFileConverter`.

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `typedef ProteinIdentification::ProteinGroup ProteinGroup` | `crate::identification::ProteinGroup` | |
| `ProtXMLFile()` | `ProtXmlFile::new`, `Default` | stateless |
| `void load(filename, ProteinIdentification&, PeptideIdentification&)` | `ProtXmlFile::load`, `load_with_limits`, `read`, `read_with_limits` | returns a `ProtXmlDocument` (both identifications and the source's log); a caller-owned `UniqueIdGenerator` replaces the singleton the identifier is drawn from |
| `void store(filename, const ProteinIdentification&, const PeptideIdentification&, document_id = "")` | `ProtXmlFile::store` | always `Error::Unsupported`, the source's `NotImplemented` |
| protected handler members (`onStartElement`, `onEndElement`, `registerProtein_`, `matchModification_`, `resetMembers_`, ...) | the private handler | driven by the crate's Xerces-faithful SAX driver, `format::sax_xml`, with `SaxEncoding::Declared` |
| (none) | `SCORE_TYPE` (`ProteinProphet probability`), `MODIFICATION_TOLERANCE` (0.001 Da), `ReadLimits` | the source's inline constants |

## Preserved source conventions

- The identifier of both identifications is the next unique ID in decimal,
  one per `<program_details>`; the search engine, date and enzyme come from
  the header; `sample_enzyme` is looked up by name, lower-cased name or
  synonym, the keys the source's `DigestionEnzymeDB` registers (shared with the
  pepXML reader, `pepxml::find_enzyme`); an unknown enzyme is the source's
  `ElementNotFound`.
- A protein with probability 0 is filtered (`[EXTRA]` section of the class
  test); protein groups and indistinguishable groups are built as the source
  builds them; a degenerate peptide appears once per group.
- Attributes are read with the source's rules: `attributeAsDouble_` through
  `StringUtils::toDouble`, `attributeAsInt_` through Xerces'
  `XMLString::parseInt` (`sax_xml::xerces_parse_int`); a missing required
  attribute is the source's `fatalError` message.
- Modifications are matched by mass within 0.001 Da among the registry's
  modifications of the residue, with the source's warnings for a missing,
  unparsable or ambiguous modification, word for word.

## Native differences

- **PSI-MOD.** The Release registry also holds PSI-MOD, so its ambiguity
  warnings can list PSI-MOD names (for example the three
  `oxidation to L-methionine ... (M)` records beside `Oxidation (M)`); the crate
  does not bundle PSI-MOD (`docs/MODIFICATION_REGISTRY_SUPPORT.md`), so those
  names, and a warning left with a single candidate, do not appear here. The
  resolved sequences are the same.
- **Undefined behaviour is an error.** An `<indistinguishable_protein>`,
  `<peptide>` or `<mod_aminoacid_mass>` without an enclosing protein or
  peptide (an empty-vector or dangling-pointer access in the source), a charge
  outside 32 bits (the source's `(int)` cast wraps it) and a non-finite
  probability, coverage or score are refused.
- A load returns nothing on error and is bounded by `ReadLimits`; log lines
  past `max_log_lines` are counted, not kept.

## Checked boundaries and evidence

- **Tier 1, executed differential.** `../oracle/F3a-protxml`: the Release
  TOPP `IDFileConverter -test` of `openms4-release-bc9cc12-c19e494-174b576`
  converted four upstream and thirteen generated protXML inputs to idXML on
  ibminode06 (AMD EPYC 7763), twice, identical. `tests/protxml.rs` reads each
  Release idXML back and compares it field by field with this port's load, the
  tool's stderr line by line with this port's log (with the PSI-MOD
  difference above modelled by name), and the exit code and message of the
  six refused loads.
- **Tier 1, retained upstream output.** `TOPP_IDFileConverter_3` (test-data
  `0cb15f2`), whose retained idXML predates `ProteinHit`'s double score.
- **Tier 3.** Every `START_SECTION` of `ProtXMLFile_test.cpp`, with its
  literals.
- **Provenance.** `tests/data/protxml_provenance.json`.
