# UniProtKB XML

`format::uniprot_xml` (feature `kb-xml`) ports `FORMAT/UniProtXMLFile.h`
(with its records `UniProtEntry` and `UniProtFeature`) and
`FORMAT/HANDLERS/UniProtXMLHandler.h` from OpenMS4-core
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. It reads the fields the UniPEFF
TOPP tool needs to write a PEFF description line; parsing runs on the shared
bounded SAX driver `src/format/sax_xml.rs`.

## API mapping

| Source member | Rust |
| --- | --- |
| `UniProtFeature` (`type`, `description`, `has_position`, `has_range`, `position`, `begin`, `end`, `original`, `variation`) | `UniProtFeature` public fields; `type` is `feature_type`; coordinates `usize` |
| `UniProtEntry` (`dataset`, `accession`, `alt_accessions`, `name`, `full_name`, `primary_gene`, `ncbi_tax_id`, `tax_name`, `protein_existence`, `sequence_version`, `entry_version`, `sequence`, `features`) | `UniProtEntry` public fields |
| `UniProtXMLFile()` / `~UniProtXMLFile()` | none needed: free functions |
| `load(filename, entries)` | `load(path) -> Vec<UniProtEntry>`, `load_with_limits`, `read` / `read_with_limits` over any `BufRead` |
| `loadStreaming(filename, callback)` | `load_streaming(path, FnMut(UniProtEntry) -> Result<()>)`, `load_streaming_with_limits`, `read_streaming`, `read_streaming_with_limits` |
| `UniProtXMLHandler` (constructor, `onStartElement`, `onEndElement`, `onCharacters`, `CaptureTarget`, `resetEntry_`, `parsePosition_`, `stripWhitespace`) | private `Handler`, `Capture`, `position`, `strip_whitespace`, member for member |
| inherited gzip/bzip2 detection of `XMLFile::parse_` | by leading bytes on every read |
| — (native) | `ReadLimits` |

## Preserved source conventions

The handler's scoping rules, which follow the upstream C# UniPEFF tool: only
elements inside an `<entry>` count; the `<comment type="alternative
products">` subtree is skipped; the first `<accession>` is the accession and
later ones are alternatives; the entry mnemonic is the first `<name>` directly
under `<entry>` outside `<gene>`, `<organism>` and `<recommendedName>`; gene
and organism names are scoped so they never leak into the mnemonic; only the
canonical `<sequence>` (with a `length` attribute) is read, white space
removed; inside a `<feature>` only `<position>`, `<begin>`, `<end>`,
`<original>` and the first `<variation>` are read, and the feature is appended
whole at its end tag; text is kept only when no child element starts inside
it; a missing attribute reads as the empty string. `parsePosition_` is
`strtol` in base 10 over the whole attribute, with a negative, partial or
empty parse meaning unknown (0).

## Native differences

- `load` returns the entries, and none if the read fails; the source clears
  its output vector and streams into it, leaving a partial result on error.
- The streaming callback returns `Result`; an error stops the read and is
  returned unchanged. The source's `std::function<void(UniProtEntry&&)>`
  cannot stop it.
- A position beyond `int` is 0 (unknown); the source's `static_cast<int>`
  keeps the low 32 bits of `strtol`'s result, and turns an overflowing
  `LONG_MAX` into -1.
- `ReadLimits` bound input bytes (the default admits a full TrEMBL release
  when streamed), depth, events, one text run, features and accessions per
  entry, and the entries and text `read` collects; exceeding one is
  `Error::InvalidValue`. A DTD or a non-UTF-8 encoding is `Error::Unsupported`.

## Checked boundaries and evidence

`tests/uniprot_xml.rs` (4 tests, `--features kb-xml`): the `load` and
`loadStreaming` sections of `UniProtXMLFile_test.cpp` over the byte-identical
fixture `tests/data/UniProtXMLFile_test.xml` (tier 3); gzip input detected by
content, and failures and bounds (tier 4). Provenance:
`tests/data/uniprot_xml_provenance.json`.
