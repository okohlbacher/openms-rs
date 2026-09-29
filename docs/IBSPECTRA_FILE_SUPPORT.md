# IBSpectraFile support

Port of `FORMAT/IBSpectraFile.h` / `IBSpectraFile.cpp` (core `bc9cc12`), the
export of isobaric consensus maps into the IBSpectra table read by the R
package isobar. Rust: `src/format/ibspectra.rs`. Tests:
`tests/ibspectra_file.rs` (7; the real example needs `consensusxml`).
Provenance: `tests/data/ibspectra_file_provenance.json`.

## API mapping

| Source member | Rust |
|---|---|
| `IBSpectraFile()` | `IBSpectraFile::new()` / `Default` |
| `IBSpectraFile(const IBSpectraFile&)`, `operator=` | `Clone`, `Copy` (the class has no members) |
| `store(filename, const ConsensusMap&)` | `store(path, &ConsensusMap)`; `write(writer, &ConsensusMap)` |
| private `guessExperimentType_` | private method; returns the channel centres of iTRAQ 4-plex, TMT 6-plex or iTRAQ 8-plex |
| private `constructHeader_` | private method |
| private `getModifString_` | private method |
| file-local `IdCSV` | private `IdCsv` |

The `ANALYSIS/QUANTITATION` method classes are not ported; only their channel
centres are needed, transcribed from the `ItraqFourPlexQuantitationMethod`,
`TMTSixPlexQuantitationMethod` and `ItraqEightPlexQuantitationMethod`
constructors.

## Preserved source conventions

- Experiment type `labeled_MS2` or `itraq`, and 4/6/8 column headers, choose
  the method; the source's two `InvalidParameter` messages are kept.
- Unidentified features write one `UNIDENTIFIED_PROTEIN` /
  `UNIDENTIFIED_PEPTIDE` line; identified ones one line per accession of the
  first hit that the first protein identification contains (others skipped);
  non-unique peptides are written (`allow_non_unique = true`); intensity 0
  writes nothing.
- Modification string, theoretical mass (full monoisotopic + charge protons),
  channel matching by truncated integer m/z (last handle wins), the 1e-5
  intensity threshold and `StringUtils::toStr` numbers are the source's.
- The table is built completely before the file is written (the source's
  `TextFile::store`).

## Native differences

- An identified feature whose first peptide identification has no hit is
  `Error::InvalidValue`; the source reads `getHits()[0]` unchecked (undefined
  behaviour).
- A skipped accession is not reported on standard error.

## Checked boundaries and evidence

- Tier 1 retained output: the class test's real example compares with the
  retained `IBSpectraFile.ibspectra.csv` through `TEST_FILE_SIMILAR`, i.e. the
  fuzzy comparator with the class-test defaults (relative 1 + 1e-5, absolute
  1e-5, `ClassTest.cpp:35-38`).
- Tier 3: the remaining `START_SECTION`s of `IBSpectraFile_test.cpp`.
