# EDTAFile support

Port of `FORMAT/EDTAFile.h` / `EDTAFile.cpp` (core `bc9cc12`), Enhanced DTA
feature tables. Rust: `src/format/edta.rs`, with the source's
`StringUtils::toDouble`/`toInt32` in `src/format/string_utils.rs` (moved there
from `src/cli/context.rs` so the command line and the formats share one copy).
Tests: `tests/edta_file.rs` (16; the feature-map round trip needs
`featurexml`). Provenance: `tests/data/edta_file_provenance.json`.

## API mapping

| Source member | Rust |
|---|---|
| `EDTAFile()`, `virtual ~EDTAFile()` | `EDTAFile::new()` / `Default`, drop |
| `load(filename, ConsensusMap&)` | `load(path, &mut UniqueIdGenerator)`, `load_with_limits`, `read(reader, filename, ids)`, `read_with_limits` |
| `store(filename, const ConsensusMap&) const` | `store_consensus(path, &ConsensusMap)`, `write_consensus(writer, &ConsensusMap)` |
| `store(filename, const FeatureMap&) const` | `store_features(path, &FeatureMap)`, `write_features(writer, &FeatureMap)` |
| private `checkedToDouble_`, `checkedToInt_` | private `checked_double`, `checked_int` |
| anonymous `TYPE_UNDEFINED/OLD_NOCHARGE/OLD_CHARGE/CONSENSUS` | private enum |

The overloaded `store` becomes two names. The unique-ID singleton the source's
`setUniqueId()` draws from is a caller-owned `UniqueIdGenerator`.

## Preserved source conventions

- The separator is chosen from the first line (tab, else space, else comma,
  else space); the three variants (RT/m/z/intensity; plus charge and meta
  columns; consensus with `RT1` sub-feature quadruplets) are told apart as the
  source does; the header is optional only for three columns.
- `NA` or an absent column is -1; a sub-feature is kept only when RT and m/z
  differ from -1; a fourth cell of a three-column file is still converted.
- Unique IDs are drawn in the source's order (map, then each consensus
  feature and each quadruplet, skipped ones included).
- One `EDTA_Map <i>` column header per sub-feature column.
- The consensus writer writes `RT/m/z/intensity/charge` and the numbered
  quadruplets padded with `NA`; the feature writer the four columns only;
  numbers in `StringUtils::toStr` (intensities as `float`). The header calls
  the consensus overload "NOT IMPLEMENTED"; the implementation writes, and so
  does this port. Map indices, unique IDs, meta data, widths, qualities and
  identifications are not written, as in the source.
- Only the consensus writer checks the extension.

## Native differences

- An empty input is `Error::Parse`; the source dereferences the missing first
  line.
- NaN, infinite or `f32`-overflowing values are `Error::Parse`; the source
  stores them.
- Lines end in LF on every platform (the source's text stream writes CRLF on
  Windows).
- Resource ceilings: the `text::Limits` of the plain-text reader.

## Checked boundaries and evidence

- Tier 3: every `START_SECTION` of `EDTAFile_test.cpp` with its literals over
  the byte-identical upstream fixtures `EDTAFile_test_1..4.edta` and
  `EDTAFile_test_out_1.featureXML`.
- Tier 4: separator choice, mandatory header, `NA` handling, the converted
  fourth cell, meta-data column names, empty lines, non-finite refusal, both
  writer layouts, the extension check, unique-ID order.
