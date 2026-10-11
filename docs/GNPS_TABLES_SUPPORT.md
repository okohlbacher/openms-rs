# GNPS tables support

Ports of `FORMAT/GNPSQuantificationFile.h` and `FORMAT/GNPSMetaValueFile.h`
(core `bc9cc12`), the feature-quantification table (required) and the
meta-value seed table (optional) of GNPS Feature-Based Molecular Networking,
both produced by the `GNPSExport` TOPP tool. Rust:
`src/format/gnps_quantification.rs`, `src/format/gnps_meta_value.rs`. Tests:
`tests/gnps_tables.rs` (8; the TOPP differentials need `consensusxml`).
Provenance: `tests/data/gnps_tables_provenance.json`. The tool's MGF output
is `docs/GNPS_MGF_FILE_SUPPORT.md`.

## API mapping

| Source member | Rust |
|---|---|
| `GNPSQuantificationFile::store(const ConsensusMap&, const std::string& output_file)` (static) | `GNPSQuantificationFile::store(&ConsensusMap, path)`; `write(&ConsensusMap, writer)` |
| `GNPSMetaValueFile::store(const ConsensusMap&, const std::string& output_file)` (static) | `GNPSMetaValueFile::store(&ConsensusMap, path)`; `write(&ConsensusMap, writer)` |

## Preserved source conventions

- Quantification table: `#MAP` and `#CONSENSUS` header lines, one `MAP` row per
  column header, one `CONSENSUS` row per feature with RT, m/z, intensity,
  charge, width, quality; the IIMN columns (`row ID`, `best ion`, `partners`,
  `annotation network number`) exactly when the first feature carries
  `row ID`; per-map columns `rt_<i>`..`width_<i>` for `i` in `0..n` headers,
  whatever the header keys; when two handles share a map index the later one
  is written.
- Meta-value table: header `<TAB>filename<TAB>ATTRIBUTE_MAPID`, rows
  `<i><TAB><basename><TAB>MAP<i>` over the primary MS-run paths; basename after
  the last `/` or `\`.
- Numbers use `StringUtils::toStr`; meta values stream as a `DataValue`;
  fields go through the `SVOutStream` rules without quoting (a tab becomes
  `_`).

## Native differences

- A destination that cannot be opened is `Error::Io`; the source opens the
  stream unchecked and writes nothing.
- A line break inside a text field is `Error::InvalidValue`, where the
  source's `SVOutStream` throws `IllegalArgument`.

## Checked boundaries and evidence

- Tier 1: `TOPP_GNPSExport_1` (`-out_quantification`,
  `GNPSExport_1_out_quant.txt`) and `TOPP_GNPSExport_4`
  (`-out_meta_values`, `GNPSExport_4_out_meta_values.tsv`), test-data
  `topp/CMakeLists.txt`, reproduced from the same consensusXML inputs. The
  meta-value table matches byte for byte; the quantification table matches
  field by field, text exactly and numbers as identical `f64`: the retained
  file spells 52 values (44 lines) inexactly in the last digits (an older
  formatter), a count derived independently with exact 15-digit rounding.
- Tier 3: the single `START_SECTION` of each class test with its literals.
