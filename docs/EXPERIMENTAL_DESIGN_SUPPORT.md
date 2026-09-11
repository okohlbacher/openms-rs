# Experimental design records, groupings and the TSV reader

This operation group ports `METADATA/ExperimentalDesign.h` and
`FORMAT/ExperimentalDesignFile.h` at SDK
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. Both public classes and every
public operation are implemented natively; no TOPP workflow that consumes a
design is certified port-ready by this group alone.

An experimental design maps quantitative values to the biological material they
were measured from. One row of the MS file section describes one channel of one
MS file; the sample section is a free-form table of factors, one named row per
sample. The complete model, glossary, worked examples and pitfalls are in the
source header and are not restated here.

## API

`metadata::{ExperimentalDesign, SampleSection, MSFileSectionEntry}` and
`format::experimental_design_file`.

| Source | Native |
| --- | --- |
| `ExperimentalDesign(fs, ss)` | `ExperimentalDesign::from_sections` — sorts, then validates |
| `get/setMSFileSection`, `get/setSampleSection` | `ms_file_section`, `set_ms_file_section`, `sample_section`, `set_sample_section` |
| `getFractionToMSFilesMapping` | `fraction_to_ms_files_mapping` |
| `getUniqueSampleRowToSampleMapping` | `unique_sample_row_to_sample_mapping` |
| `getSampleToPrefractionationMapping` | `sample_to_prefractionation_mapping` |
| `getConditionToSampleMapping` | `condition_to_sample_mapping` |
| `getSampleToConditionMapping` | `sample_to_condition_mapping` |
| `getConditionToPathLabelVector` | `condition_to_path_label_vector` |
| `getPathLabelTo{Sample,Fraction,FractionGroup,Prefractionation,Condition}Mapping` | `path_label_to_*_mapping(basename_only)` |
| `getNumberOf{Samples,Fractions,Labels,MSFiles,FractionGroups}` | `number_of_*` |
| `getSample`, `isFractionated`, `sameNrOfMSFilesPerFraction`, `filterByBasenames` | `sample`, `is_fractionated`, `same_nr_of_ms_files_per_fraction`, `filter_by_basenames` |
| `from{ConsensusMap,FeatureMap,Identifications}`, `annotateColumnHeaders` | `from_consensus_map`, `from_feature_map`, `from_identifications`, `annotate_column_headers` |
| `SampleSection::get{Samples,Factors,FactorValue,FactorColIdx,SampleName,SampleRow,ContentSize}`, `addSample`, `hasSample`, `hasFactor` | `samples`, `factors`, `factor_value`, `factor_value_by_row`, `factor_column_index`, `sample_name`, `sample_row`, `len`, `add_sample`, `has_sample`, `has_factor` |
| `ExperimentalDesignFile::load` | `experimental_design_file::{load, load_with_warnings, load_text}` |

`ColumnHeader::label_as_uint` ports `ConsensusMap::ColumnHeader::getLabelAsUInt`:
the annotated `channel_id` plus one, or 1 when absent. The `label` string is not
parsed, as in source.

## Preserved source conventions

Row sorting is `(fraction_group, fraction, label, sample, path)`. The validator
runs on construction and enforces unique `(fraction group, fraction, label)` and
`(path, label)`, fraction groups that are consecutive integers from 1, and one
sample per `(fraction group, label)` in a design with a single distinct label.
An empty file section is valid and unchecked.

Condition grouping ignores `Sample` and every column whose **name** contains
`replicate` or `Replicate`; `Donor` or `Rep` still split conditions. Condition
numbers are the lexicographic rank of the retained factor-value tuple, which is
built in the alphabetical order of the column names — not a reference level. The
weaker prefractionation rule ignores only `Sample` and keeps replicate columns.

A factor-less sample section — what `from_consensus_map`, `from_feature_map` and
`from_identifications` build — keeps the source's disagreement:
`sample_to_condition_mapping` and `sample_to_prefractionation_mapping` give every
sample its own group, while `condition_to_sample_mapping`,
`condition_to_path_label_vector` and `path_label_to_condition_mapping` collapse
all samples into one condition. Both `Sample*` mappings are keyed by sample
**name**, never by the stringified row index.

Both file layouts are auto-detected by the source detector, which scans every
line rather than only headers, does not trim cells and does not skip comment
lines. The one-table parser accepts unknown columns as sample metadata and
appends an absent `Label` or `Sample` column; without a `Sample` column the
fraction-group value becomes the sample name, and `Label > 1` is then rejected.
The two-table file section rejects unknown headers. Cells are trimmed, lines
starting with `#` are ignored, and a relative `Spectra_Filepath` is resolved
against the design file's directory, then the working directory, and is
otherwise kept as written.

`add_sample` keeps the source behavior of interning a repeated name to its first
row while still appending a content row, so `len` can exceed the number of
distinct names. `sample_name` reads the `Sample` column when the section has one
and otherwise the name store, so an inferred design answers too.

## Native differences

Every source throw becomes a typed `Error`; parse failures carry a one-based
line number. Four source behaviors are tightened rather than reproduced:

- A sample-section or data row shorter than its column map is rejected before
  access. The source reads past the end of the row vector ([CPP-059](../OpenMS_CPP_ISSUES.md)).
- A negative `Fraction_Group`, `Fraction` or `Label` is rejected. The source
  parses to `int32` and assigns to `unsigned`, wrapping `-1` to 4294967295
  ([CPP-060](../OpenMS_CPP_ISSUES.md)).
- A two-table file section naming a sample the sample section does not contain
  returns a typed error instead of a bare `std::out_of_range`.
- A nonnumeric or negative `fraction`/`fraction_group` column-header meta value
  is rejected instead of taking the source's unchecked `DataValue` cast
  ([CPP-058](../OpenMS_CPP_ISSUES.md)).

The source writes two advisory diagnostics to its global log stream: mismatched
factors for one sample in the one-table parser, and the removed-file count and
empty-design notice in `filterByBasenames`. This port has no global log stream,
so `load_with_warnings` returns the first as owned strings and
`filter_by_basenames` returns the removed count; an emptied design is visible
through `ms_file_section`. `load` discards the warnings.

`getUniqueSampleRowToSampleMapping` has a source debug assertion on a factor-less
section. This port takes the release behavior: every sample groups under the
empty tuple.

## Checked boundaries and evidence

`ExperimentalDesign::MAX_ROWS` (1,000,000) and `MAX_BYTES` (256 MiB, a
conservative logical estimate rather than measured allocation) are precharged
before any section is stored. Reading uses the existing `format::text` line
limits through `ReadOptions::limits`.

The nine source design fixtures are retained byte-for-byte under
`tests/data/experimental_design_*.tsv`; their hashes and those of the source
header, implementation and class test are in
[experimental_design_provenance.json](../tests/data/experimental_design_provenance.json).
Expected values are the literal `ExperimentalDesign_test.cpp` assertions plus
independently derived cases for the condition/prefractionation mappings, the
parser error paths and the resource bounds. No C++ executable was run for this
group. The consensus and feature fixtures of the source test are consensusXML
and featureXML documents; this port builds the equivalent column headers and
`spectra_data` metadata directly instead of reloading them.
