# Splitting identifications by file origin

`analysis::id_ripper` ports `IDRipper` from OpenMS4-core revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. `IDRipper::rip` takes borrowed protein/peptide identification slices and returns owned `RippingResult` partitions. It opens and writes no files and leaves both inputs unchanged.

## Origin conventions

Every peptide identification, including empty ones, must contain exactly one common origin annotation:

| Annotation | Partition identity and path |
| --- | --- |
| `file_origin` | A string path. Distinct values receive indices in first-appearance order across all peptides. |
| `map_index` | A nonnegative signed-32-bit index into the matching protein run's `primary_ms_run_paths` (`spectra_data` in C++). |
| `id_merge_index` | Same indexed-path rule, using the newer source annotation key. |

`detect_origin_annotation_format` exposes this check separately. Mixed conventions, multiple origin keys, missing annotations, the source `unknown` annotation and an empty peptide list cannot be autodetected and return an error. Indexed metadata follows the source text-to-integer conversion: typed integers and numeric strings are accepted, as are integral floats whose string form is an integer. Negative, fractional, overflowing and out-of-range values fail explicitly.

`split_identification_runs` preserves the zero-based input protein-run index in each `RipFileIdentifier`; otherwise it is `None` and origins combine runs. Files are returned sorted by run index and then origin index. Within a partition, protein runs follow their first peptide occurrence; peptide records retain input order. Protein hits are added once per accession, sorted within each newly encountered peptide's accession set.

`origin_fullname` retains the complete path. `output_basename` removes the last filename extension using the host's native path semantics (`sample.raw.mzML` becomes `sample.raw`). With `numeric_filenames=false`, different numeric output identities may not have the same basename. Setting it true bypasses this name collision check; callers can use the numeric identity fields to choose filenames. This library does not invent an output directory or write those files.

The selected origin metadata is removed from output peptide and protein-run copies. In indexed modes, each output run's primary MS paths shrink to its one origin path. `file_origin` mode preserves the original primary-path list. Raw MS paths remain unchanged in every mode, matching source metadata copying.

## Protein records and skipped identifications

Each output run receives only the proteins referenced by its retained peptide evidence. Scores, sequence, rank, modification positions, search parameters, dates, typed metadata and other represented fields retain their values. Coverage likewise retains its original scope; splitting does not recalculate per-file coverage or confidence. A protein accession is resolved **within its own identification run**. This deliberately fixes the source's global accession table, where the last run's protein could overwrite another run's sequence, scores or metadata.

Both protein-group vectors are omitted from output copies, following `ProteinIdentification::copyMetaDataOnly` used by IDRipper. `groups_not_copied` counts the group records omitted across those output copies. Input groups remain available unchanged; inference or group cleanup can be applied explicitly to the partitions if appropriate. Group probabilities are not silently recomputed from a subset of members.

The source skips empty peptide identifications and IDs whose hits have no protein accessions. The native result additionally retains these as `skipped_peptide_identifications`, with their original metadata. This makes every input peptide record accountable. Origin/run/path and basename checks still occur before that skip, matching source order.

## Checked source differences

- Run identifiers and accessions within one run must be unique. Identical accessions across different runs are supported with independent records.
- Every referenced protein must exist in the matching run. Partial or completely missing references return an error instead of producing dangling evidence or silently omitting the peptide.
- Without run splitting, the same numeric origin index cannot refer to different full paths in different runs, even when numeric filenames are enabled. The C++ key ignores the path and can combine unrelated files and relabel them with the first path. Choose run splitting to preserve these distinct origins.
- Empty origin paths and unusable nonnumeric basenames are errors.
- Failures return no partial partition result and never mutate input metadata. The C++ API removes origin annotations from caller-owned inputs while processing and may fail after partial mutation.

The implementation uses standard-library ordered maps/sets and owned output copies. Memory grows with the supplied records and number of resulting run copies; this is an in-memory operation rather than a streaming file splitter.

## Evidence

The pinned `IDRipper_test.cpp` contains constructor checks and an unimplemented ripping test, so it provides no partition goldens to claim as executed reference coverage. `tests/id_ripper.rs` instead checks independently hand-calculated partitions derived from the implementation: all three origin conventions, first-appearance and numeric ordering, basename extraction/collisions, indexed path reduction, hit union order, metadata/group behavior, multiple runs, retained skipped records, and checked reference/index errors. It explicitly tests the native corrections for run-specific protein values and differing paths sharing one numeric index.

Pinned references: `src/openms/source/ANALYSIS/ID/IDRipper.cpp`, its corresponding public header, `src/openms/source/METADATA/ProteinIdentification.cpp`, and `src/tests/class_tests/openms/source/IDRipper_test.cpp`. They are recorded in the [source inventory](source-inventory.json). No C++ reference build or execution occurred; derived code and tests retain BSD-3-Clause attribution. The file-writing TOPP tool is outside this core adapter.
