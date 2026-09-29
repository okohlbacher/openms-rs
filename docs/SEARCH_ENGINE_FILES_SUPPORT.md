# Search-engine input and output files

The adapters for the input (parameter) files and result files of the Sequest,
Inspect, PepNovo, OMSSA and X! Tandem search engines, and the three text feature
formats written by SpecArray, msInspect and Kroenik, port these OpenMS4-core
`bc9cc12` headers, each with its implementation:

| Source header | Rust module | Feature |
| --- | --- | --- |
| `FORMAT/SequestInfile.h` | `format::sequest_infile` | `idxml` |
| `FORMAT/SequestOutfile.h` | `format::sequest_outfile` | none |
| `FORMAT/InspectInfile.h` | `format::inspect_infile` | `idxml` |
| `FORMAT/InspectOutfile.h` | `format::inspect_outfile` | `mzml` |
| `FORMAT/PepNovoInfile.h` | `format::pepnovo_infile` | none |
| `FORMAT/PepNovoOutfile.h` | `format::pepnovo_outfile` | none |
| `FORMAT/OMSSAXMLFile.h` | `format::omssa_xml_file` | `idxml` |
| `FORMAT/OMSSACSVFile.h` | `format::omssa_csv_file` | none |
| `FORMAT/XTandemXMLFile.h` | `format::xtandem_xml_file` | `idxml` |
| `FORMAT/XTandemInfile.h` | `format::xtandem_infile` | none |
| `FORMAT/SpecArrayFile.h` | `format::spec_array_file` | none |
| `FORMAT/MsInspectFile.h` | `format::ms_inspect_file` | none |
| `FORMAT/KroenikFile.h` | `format::kroenik_file` | none |
| `FORMAT/PTMXMLFile.h`, `FORMAT/HANDLERS/PTMXMLHandler.h` | `format::ptm_xml_file` | `idxml` |

`PTMXMLFile` is not in the lane's list; it is ported because `InspectInfile`
and `SequestInfile` look modification names up in PTMXML tables. The `idxml`
modules need the crate's XML reader; `InspectOutfile` needs the mzXML and mzData
readers of `mzml`. No dependency was added.

Each module's rustdoc states the file grammar, every error and each divergence
from the source at the item it concerns. This page gives the mapping, the shared
machinery and the evidence.

## Shared machinery

- `data_structures::string_utils` ports `StringUtils::toDouble`, `toInt32`,
  `split` and `trim` with their acceptance rules. `cli::context` now re-exports
  its `to_double`/`to_int32` instead of keeping a private copy.
- `format::search_engine_text` holds what the text adapters share. It turns a
  sequence into the crate's strict `AASequence::parse` input the way the source's
  default permissive `AASequence::fromString` reads it (`*`, `#` and `+` become
  `X` outside a modification bracket, spaces are skipped). It also has a bounded
  whole-file reader that splits lines where `std::getline` does, the
  `std::string` search helpers with `npos` arithmetic, and `flanking_residue`,
  which the Sequest and Inspect result readers share.
- The X! Tandem and OMSSA readers are driven by the crate's one Xerces-faithful
  SAX driver, `format::sax_xml`, which every `XMLHandler` port shares; lane F3b's
  own driver was folded into it at wave-4 integration, and
  `format::search_engine_xml` keeps only the public `XmlReadLimits`. The driver
  delivers character data where the Release build's Xerces delivers it, split at
  every reference and ended by markup, a comment or a processing instruction,
  with white space inside the root element reported (the OMSSA handler resets
  state on the white space between elements; both handlers act on the first
  `characters` call of a text node, as the source does). It refuses a DTD. X!
  Tandem documents are read as Latin-1, as the source's
  `XMLFile::enforceEncoding_` does. The splitting is recorded from the Release
  build by `../oracle/integ-F-xerces-sax` (`tests/data/xerces_sax_release`).
- `inspect_infile::parse_modification_line` is the `handlePTMs` grammar that
  `InspectInfile` and `SequestInfile` share.
- `resources/search_engines/OMSSA_modification_mapping` is the source's
  `share/OpenMS/CHEMISTRY/OMSSA_modification_mapping`, embedded byte for byte.

## API mapping

C++ getters `getX()` are Rust `x()`, and setters keep the `set_` prefix. A C++
output parameter becomes a return value, or a result struct when there are
several. Copy construction, assignment and `operator==` become `Clone` and
`PartialEq`. The derived C++ classes whose parent only supplies a file name
(`Internal::XMLFile`) take a path, or a reader through `load_reader`.

| Source | Rust |
| --- | --- |
| `KroenikFile`, `MsInspectFile`, `SpecArrayFile::load(filename, FeatureMap&)` | `load(path) -> FeatureMap`, `load_reader` |
| their `store` templates (throw `NotImplemented`) | `store` returning `Error::Unsupported` |
| `OMSSACSVFile::load(filename, ProteinIdentification&, ids&)` | `load(path) -> Vec<PeptideIdentification>`; the source never touches the protein identification |
| `OMSSAXMLFile::load(filename, protein, ids, load_proteins, load_empty_hits)` | `load(path, OMSSALoadOptions) -> OMSSAXmlResult` |
| `OMSSAXMLFile::setModificationDefinitionsSet` | `set_modification_definitions` |
| private `readMappingFile_` | `new()` (bundled mapping), `with_mapping(text)`, `mapped_modifications` |
| `XTandemXMLFile::load(filename, protein, ids, mod_def_set)` | `load(path, &ModificationDefinitionsSet) -> XTandemXmlResult` |
| `XTandemInfile` setters and getters, `ErrorUnit`, `MassType` | same names, `ErrorUnit`, `MassType` enums |
| `XTandemInfile::write(filename, ignore_member_parameters, force_default_mods)` | `write(path, ..)`, `write_to(writer, ..)` |
| `SequestOutfile::load` | `load(path, &mut ProteinIdentification, p_value_threshold, pvalues, ignore_proteins_per_peptide) -> Option<PeptideIdentification>` |
| `getColumns`, `getSequences`, `getACAndACType`, `readOutHeader` | `get_columns`, `get_sequences`, `get_ac_and_ac_type`, `read_out_header -> SequestOutHeader` |
| `SequestInfile::store`, `getEnzymeInfoAsString`, `setEnzyme`, `handlePTMs` | `store`, `to_params_string`, `enzyme_info_as_string`, `set_enzyme`, `handle_ptms` |
| `SequestInfile` setters and getters (`getMaxAAPerModPerPeptide`, …) | same names (`max_aa_per_mod_per_peptide`, …) |
| protected `addEnzymeInfo`, `getEnzymeInfo_` | `add_enzyme_info`, `enzyme_info` |
| `PepNovoInfile::setModifications`, `getModifications`, `store` | `set_modifications`, `modifications`, `modification_definitions`, `lines`, `write`, `store` |
| `PepNovoOutfile::load(result_filename, ids, protein, score_threshold, index_to_precursor, fnames, modkey_to_mod_id)` | `load(path, score_threshold, &IndexPosMapping, &modkey_to_mod_id) -> PepNovoResult` |
| `PepNovoOutfile::getSearchEngineAndVersion` | `get_search_engine_and_version` |
| `InspectInfile::store`, `handlePTMs` and every accessor | `store`, `write`, `handle_ptms`, same-name accessors (`set_max_ptm_size` for `setMaxPTMsize`) |
| `InspectOutfile::load` | `load`, `load_at` (fixed run date) |
| `getWantedRecords`, `compressTrieDB`, `generateTrieDB` | `wanted_records`, `compress_trie_db`, `generate_trie_db` |
| `getACAndACType`, `getPrecursorRTandMZ`, `getLabels` | `ac_and_ac_type`, `precursor_rt_and_mz`, `labels -> DatabaseLabels` |
| `getSequences`, `getExperiment`, `getSearchEngineAndVersion`, `readOutHeader` | `sequences`, `experiment -> (MSExperiment, type name)`, `search_engine_and_version`, `read_out_header -> InspectOutHeader` |
| `PTMXMLFile::load`, `store` (and `PTMXMLHandler`) | `load`, `load_reader`, `store`, `write` |

Two declared C++ members have no Rust counterpart, because the pin gives them no
behaviour: `XTandemInfile::getNoiseSuppression`/`setNoiseSuppression` are
declared but never defined, and no member holds the value.

## Divergences, in brief

Each is documented in full at the Rust item.

- Untrusted input never panics. Every whole-file read, line, record, XML depth,
  attribute and text node is bounded, and the bound is checked before storage
  grows. Where the source has undefined behaviour, the port returns an error or
  a documented value instead. Examples: X! Tandem `aa` outside a group or before
  a domain, an OMSSA CSV first data line numbered -1, a PepNovo title past the
  line end, and a Sequest iterator incremented past the end.
- On error a loader returns nothing. The source has often already changed the
  caller's containers before it throws.
- `PTMXMLFile` escapes `&`, `<` and `>` when writing. The source writes the text
  raw, so it can produce a malformed document.
- Source behaviour that looks like a defect but is kept for parity, and named as
  a defect candidate in the rustdoc: `InspectOutfile::load` gives each
  identification the scan of the next spectrum. `XTandemXMLFile` keeps
  protein uids across loads in C++; each Rust load starts afresh, because a
  reader holds no state.

## Evidence

| Evidence | Tests |
| --- | --- |
| Upstream class tests (`src/tests/class_tests/openms/source/<Name>_test.cpp` at `bc9cc12`), fixtures copied unmodified | `tests/search_engine_feature_files.rs`, `search_engine_id_files.rs`, `search_engine_infiles.rs`, `search_engine_ptm_infiles.rs`, `search_engine_xml_files.rs`, `search_engine_inspect_outfile.rs` |
| Executed C++ Release build (ibminode06, `openms4-release-bc9cc12-c19e494-174b576`), two runs, byte-identical after masking the run date | `tests/search_engine_oracle.rs`, `search_engine_infiles.rs`, `search_engine_ptm_infiles.rs`, `search_engine_inspect_outfile.rs` against `tests/data/search_engine_oracle/` |

`tests/data/search_engine_files_provenance.json` records the sha256 of every
source file, fixture and oracle output. The oracle's driver is
`tests/data/search_engine_oracle/driver.cpp.txt`. Its host, compiler and hashes
records sit beside the outputs.

A class-test check that cannot run at the pin is replaced by an independently
derived or oracle expectation, and the test says so. Examples:
`InspectOutfile_test`'s `getSequences` section guards its checks with
`found_sequences.size() == 1` after asking for two records, and
`TEST_FILE_EQUAL` compares trimmed lines, so it does not see the trailing line
feed of the `.trie` fixtures.
