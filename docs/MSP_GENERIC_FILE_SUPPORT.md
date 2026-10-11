# MSPGenericFile support

Port of `FORMAT/MSPGenericFile.h` / `MSPGenericFile.cpp` (core `bc9cc12`),
metabolite spectral libraries in MSP text. Rust: `src/format/msp_generic.rs`
(the source's Boost expressions through the crate's `BoostRegex` facade;
`std::stod`/`std::stoul` and `StringUtils::toDouble` in
`src/format/string_utils.rs`). Tests: `tests/msp_generic_file.rs` (9).
Provenance: `tests/data/msp_generic_file_provenance.json`. Peptide MSP
libraries are `MSPFile`, `docs/MSP_FILE_SUPPORT.md`.

## API mapping

| Source member | Rust |
|---|---|
| `MSPGenericFile()` | `MSPGenericFile::new()` / `Default` |
| `MSPGenericFile(filename, MSExperiment& library)` | `MSPGenericFile::new()` then `load(filename)` |
| `~MSPGenericFile()` | drop |
| `getDefaultParameters(Param&)` | `MSPGenericFile::defaults()`; also `to_param`, `from_param`, `from_param_with_warnings` |
| private `updateMembers_` | `from_param` (reads `synonyms_separator`) |
| `load(filename, MSExperiment& library)` | `load(filename) -> MSExperiment`, `load_with_limits`, `read`, `read_with_limits` |
| `store(filename, const MSExperiment& library) const` | `store` (round-trip numbers), `store_with_options` (`StoreOptions::source()` for the source's `%g`), `write_with_options` |
| private `addSpectrumToLibrary` / `MSPGenericFile_friend::addSpectrumToLibrary` | `add_spectrum_to_library` (public, as the friend class exposes it to the class test) |
| private `loaded_spectra_names_`, `synonyms_`, `synonyms_separator_` | private fields; `synonyms_separator()` accessor |
| `is_valid` meta key | `IS_VALID` |

## Preserved source conventions

- Line classification in the source's order and with its expressions
  (peaks, synonyms, `Name`, `Num Peaks`, retention time, `PRECURSORMZ`,
  `CAS#/NIST#`, `INCHIKEY`, `SMILES`, `FORMULA`, `PRECURSORTYPE`, `CCS`, any
  other `key: value`), including which expressions use `no_mod_s` and so keep
  a trailing CR on CRLF files.
- `addSpectrumToLibrary`'s checks and their order (name, duplicate,
  `Num Peaks`), the in-place changes to the offered spectrum, and a missing
  retention time replaced by the library index.
- The writer layout (five `<mz>:<intensity> ` points per line, `Synon` lines,
  `CAS#/NIST#`).

## Native differences

- The source removes the last character after joining synonyms, wrong for a
  separator longer than one character; this port joins with the separator
  (identical for the default `|`).
- The source writes meta values in `MetaInfoRegistry` registration order
  (process-dependent); this port sorts by key.
- `store` defaults to numbers that read back exactly; the source's six-digit
  `%g` (which writes `1.5e+06` that its own reader cannot read back) is
  `StoreOptions::source()`.
- All spectra are checked for a name before the file is created; the source
  throws after partial output. A line break in a written field, and NaN or
  infinite numbers, are refused.
- `std::stod`'s escaping `std::invalid_argument`/`out_of_range` are
  `Error::Parse`; non-finite values and non-UTF-8 fields are refused.
- Resource ceilings: `ReadLimits` (line length, input size, 200 million
  peaks, 10 million spectra).

## Checked boundaries and evidence

- Tier 3: every `START_SECTION` of `MSPGenericFile_test.cpp` with its
  literals, including the friend-class section, over the byte-identical CRLF
  fixture `MSPGenericFile_input.msp`.
- Tier 4: line classification order and special fields, duplicates and
  missing retention times, both store number formats, synonym separators,
  refused line breaks. The read limits are not exercised by a test.
