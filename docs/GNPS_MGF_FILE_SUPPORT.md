# GNPSMGFFile support

Port of `FORMAT/GNPSMGFFile.h` / `GNPSMGFFile.cpp` (core `bc9cc12`), the
consensus-spectrum MGF writer for GNPS Feature-Based and Ion-Identity
Molecular Networking, used by the `GNPSExport` TOPP tool. Rust:
`src/format/gnps_mgf.rs` (features `consensusxml` and `mzml`). Tests:
`tests/gnps_mgf_file.rs` (9). Provenance:
`tests/data/gnps_mgf_provenance.json`. The tool's two other outputs are
`docs/GNPS_TABLES_SUPPORT.md`.

## API mapping

| Source member | Rust |
|---|---|
| `GNPSMGFFile()` | `GNPSMGFFile::new()` / `Default` |
| `DefaultParamHandler` defaults (`output_type`, `peptide_cutoff`, `ms2_bin_size`, `merged_spectra:cos_similarity`) | `defaults`, `to_param`, `from_param`, `from_param_with_warnings`; fields `output_type: OutputType`, `peptide_cutoff: i32`, `ms2_bin_size: f64`, `cos_similarity: f64` |
| `ProgressLogger` base | not ported; no progress is reported |
| `store(consensus_file_path, mzml_file_paths, out)` | `store(consensus_path, &[mzml paths], out)`; `write(&ConsensusMap, &[mzml paths], writer)` for an already loaded map |
| private `DEF_COSINE_SIMILARITY`, `DEF_MERGE_BIN_SIZE`, `DEF_PEPT_CUTOFF`, `DEF_MSMAP_CACHE` | public constants of the same names (`DEF_MERGE_BIN_SIZE` is the `float` 0.02 widened, as in the source) |
| file-local `binPeaks_`, `flattenAndBinSpectra_`, `writeMSMSBlockHeader_`, `writeMSMSBlock_`, `sortElementMapsByIntensity_`, `getElementPeptideIdentificationsByElementIntensity_` | private `flatten_and_bin`, the header/peak writes in `write`, `element_peptides` |

## Preserved source conventions

- Feature charge is the maximum of the feature's and its handles' charges;
  handles are ordered by intensity; each handle's map contributes the first
  peptide identification with matching integer `map_index` and a
  `spectrum_index`.
- `most_intense` keeps one annotation; `merged_spectra` keeps
  `peptide_cutoff` (-1 all) and adds every kept spectrum, the best one again
  included, whose binned contrast angle with the best spectrum reaches the
  threshold (bins of `ms2_bin_size` Da, spread 1, offset 0).
- A best spectrum without peaks skips the feature.
- Peaks are sorted by m/z and binned by the source's run rule (a bin closes
  when a peak lies more than the bin width above the bin's first m/z; closed
  bins with non-positive intensity sum are dropped, the last never is).
- Header numbers use `StringUtils::toStr`; peaks `std::fixed` with four
  decimals; blocks end `END IONS` and an empty line.

## Native differences

- A feature without a usable annotation (or `peptide_cutoff` 0) is
  `Error::MissingInformation`; the source reads `pepts[0]` of an empty vector.
- A `map_index` naming no input file, or a `spectrum_index` negative or past
  the file's spectra, is `Error::InvalidValue`; the source indexes unchecked.
- A non-integer `map_index`/`spectrum_index` is `Error::InvalidValue`, as the
  source's `DataValue` `(int)` conversion throws `ConversionError`
  (`DataValue.cpp:519-527`).
- A non-indexed mzML is an error from the on-disc reader; the source ignores
  `openFile`'s result.
- `std::sort` leaves equal intensities and equal m/z values in unspecified
  order; this port sorts stably (input order kept).
- `store` returns `Error::Io` when the output cannot be written; the source
  writes to an unchecked stream.

## Checked boundaries and evidence

- Tier 1: `TOPP_GNPSExport_1/2/3` (test-data `topp/CMakeLists.txt:989-1001`)
  reproduced with the parameters of `GNPSExport_1_mostint.ini`,
  `GNPSExport_2_merged.ini`, `GNPSExport_3_binsize.ini`, the same
  consensusXML and mzML inputs, against the retained C++ outputs. Every line
  is equal except 14 header values per file (8 `PEPMASS`, 6 `RTINSECONDS`)
  that the retained files spell inexactly in the last digit, e.g.
  `242.627651999999984` for the double `242.62765199999998344...`; each reads
  back as the identical `f64`, which is FuzzyDiff's criterion with no
  tolerance (upstream compares with ratio 1.01 / absdiff 0.01). The count 14
  was derived independently with Python's `decimal`. All peak lines match
  byte for byte, including both merge modes.
- Tier 4: defaults and restrictions of the constructor, parameter round trip,
  and the three refusals above.
- There is no upstream class test for this header.
