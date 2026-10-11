# TransitionTSVFile support

Ports `ANALYSIS/OPENSWATH/TransitionTSVFile.h` (with its `.cpp`) at
OpenMS4-core `bc9cc12`. Rust: `src/analysis/openswath/transition_tsv_file.rs`
(no feature gate; its tests need `traml`).

## API mapping

| C++ member | Rust | Notes |
|---|---|---|
| `TransitionTSVFile()` / dtor | `TransitionTSVFile::new()` / `Default` | |
| `DefaultParamHandler` (`retentionTimeInterpretation`, `override_group_label_check`, `force_invalid_mods`) | `defaults()`, `to_param()`, `from_param()`, `from_param_with_warnings()`; fields `retention_time_interpretation`, `override_group_label_check`, `force_invalid_mods` | the `peak_picking.rs` pattern; values, descriptions, valid strings and the `advanced` tag as upstream |
| `updateMembers_` | `from_param` (private `from_complete_param`) | |
| `TSVTransition` (protected) | `TsvTransition` (public), `is_peptide()` | every field; `-1` sentinels kept |
| `header_names_` | `HEADER_NAMES` | |
| `strarray_` | not ported | declared, not defined or used |
| `getTSVHeader_` | `TransitionTSVFile::tsv_header(line)` | |
| `readUnstructuredTSVInput_` | `read_transitions(reader, FileType)` | returns the warnings it prints |
| `spectrastRTExtract`, `spectrastAnnotationExtract` | private `spectrast_rt_extract`, `spectrast_annotation_extract` | |
| `cleanupTransitions_` | private `cleanup_transition` | |
| `TSVToTargetedExperiment_` | `to_targeted_experiment(&mut [TsvTransition], &mut warnings)` | public, as `TransitionPQPFile` reuses it |
| `resolveMixedSequenceGroups_` | `resolve_mixed_sequence_groups` | |
| `createTransition_`, `createProtein_`, `createPeptide_`, `createCompound_`, `addModification_` | private `create_transition`, `create_protein`, `create_peptide`, `create_compound`, `modification` | |
| `interpretRetentionTime_` | `interpret_retention_time(rt) -> RetentionTime` | |
| `convertTransition_` | crate-internal `convert_transition` | reused by `TransitionPQPFile` |
| `writeTSVOutput_` | `write_targeted_experiment(writer, &exp)` | |
| `convertTargetedExperimentToTSV` | `convert_targeted_experiment_to_tsv(path, &exp)` | |
| `convertTSVToTargetedExperiment` (TraML structure) | `convert_tsv_to_targeted_experiment(path, FileType)`, `..._with_warnings`, `read_targeted_experiment(reader, FileType)` | |
| `validateTargetedExperiment` | `validate_targeted_experiment(&exp)` | |
| `convertTSVToTargetedExperiment` (`LightTargetedExperiment`), `streamTSVToLightTargetedExperiment_` | `convert_tsv_to_light_targeted_experiment(path, type)` / `read_light_targeted_experiment(reader, type)` (group O integration). The lines are parsed by `read_transitions` and the experiment is built with the streaming reader's rules: the retention time as read (`retentionTimeInterpretation` is not applied), modifications parsed from the sequence (a failure leaves them empty), a protein's UniProt accession at its own position, the label group check inline and in file order | the source streams; the port reads the transitions first (bounded by the reader's limits) |
| `convertLightTargetedExperimentToTSV` | `convert_light_targeted_experiment_to_tsv(path, exp)` / `write_light_targeted_experiment(writer, exp)` | the same 29 columns; `NA` for an empty group label, gene, sequence, full name, annotation and protein ids; collision energy -1 |
| `ProgressLogger` base | not ported | no progress reporting |

## Preserved source conventions

- Delimiter detection (tab, `;`, `,`; four fields), header and field cleanup
  (`"`, `'`, `,` removed), the synonym order of every column, `NA` handling,
  `1`/`0`/`TRUE`/`FALSE` flags, `TextFile::getLine` line endings, SpectraST MRM
  columns and annotation parsing, generated transition names (line number) and
  group ids (`<sequence>_<charge>`).
- A field count must equal the number of *distinct* header names.
- Mixed peptide label groups get the transition group id as label unless
  `override_group_label_check`.
- Modifications: N-terminal first, then C-terminal, then per residue; unknown
  UniMod record id is `-1`.
- Writer: numbers as `StringUtils::toStr(double)` (`500.0`,
  `628.434999999999945`), flags `0`/`1`, `NA` defaults as upstream, the
  `b-H20`/`y-H20` spelling (zero, not O).

## Native differences

- An unreadable file is `Error::Io`; the source reads it as empty and fails
  on the header.
- A missing `PrecursorMz` column is refused. The source reads it through
  `header_dict["PrecursorMz"]`, whose `operator[]` inserts index 0, so it
  parses the first column and then fails the next line's field count
  (`TransitionTSVFile.cpp:267`).
- Fragment number outside 0..=255 is refused; the source wraps it into the
  `unsigned char` ordinal.
- A mass-tag modification without a registry average mass
  (`PEPTIDEK[+8.0]`-style user modifications) is `Error::Unsupported`; the
  source creates a user-defined `ResidueModification`.
- `StringUtils::toDouble`/`toInt32` are approximated by trimmed Rust `parse`.
- Resource ceilings `ReadLimits::{max_bytes, max_transitions}`.
- Warnings are returned instead of logged.

## Checked boundaries and evidence

`tests/transition_tsv_file.rs` (10 tests; `--features semantic-validation
--test transition_tsv_file` on dax: 10 passed):

- Tier 3: `TransitionTSVFile_test.cpp` heavy sections
  (`convertTargetedExperimentToTSV`, `convertTSVToTargetedExperiment`) on
  `MRMAssay_detectingTransistionCompound_input.TraML`.
- Tier 1: the `TargetedFileConverter` heavy path (openswath `e29e7ee`,
  `TargetedFileConverter.cpp:229-276`) reproduced against retained C++ outputs
  at test-data `0cb15f2` (`topp/CMakeLists.txt:1736-1864`): ConvertTSVToTraML 1
  to 6 (2 with `minutes`, 3 with `seconds`, 4 the SpectraST MRM list),
  TargetedFileConverter 13 (TSV to TraML), 5 and 6 (TraML to TSV) and 8 (TSV to
  TraML to TSV with `override_group_label_check`, equal to its input). Lines are
  compared token by token with numbers as identical doubles (FuzzyDiff is the
  upstream contract).
- Tier 4: header detection, short lines, missing columns, label-group
  resolution both ways, invalid references refused.

## C++ defect candidates

1. `TransitionTSVFile.cpp:1145` stores `fragment_nr` into the
   `unsigned char` ordinal without a range check.
2. `TransitionTSVFile.cpp:195` never checks that the input opened.
3. `extractName<bool>` (`TransitionTSVFile.cpp:70-86`) treats an invalid flag
   as absent and silently tries the next synonym.
4. `TransitionTSVFile.cpp:1648`/`1651` write `b-H20`/`y-H20` (digit zero).

## Light path evidence (group O integration)

Tier 1: `../oracle/integ-O-tsvlight/driver.cpp` (ibminode06, Release install
`openms4-release-bc9cc12-c19e494-174b576`, g++ 14.4.0, run twice byte-identical) loaded ten cases with
`convertTSVToTargetedExperiment(Light)` (the TSV inputs `ConvertTSVToTraML_1/2/5/6`,
`TargetedFileConverter_8/12`, the SpectraST list `ConvertTSVToTraML_4_input.mrm`, and
`retentionTimeInterpretation=minutes`, `override_group_label_check=true` and `force_invalid_mods=true`
variants) and wrote each with `convertLightTargetedExperimentToTSV`. `light_path_matches_the_release_build`
replays the ten cases and compares the dump of the light experiment (every field, `%.17g`) and the TSV text.
Tier 3: the class test section "[EXTRA] Light path TSV roundtrip"
(`class_test_light_path_tsv_roundtrip`); the other two light sections of the class test go through
`TransitionPQPFile`'s light members.
