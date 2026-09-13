# MzMLSqliteHandler support

Status: the installed public handler API is implemented with the checked native
policies below. All 13 source class-test sections have an assertion-family mapping
and focused remote native tests. Integration-wide validation is tracked separately;
there is no executed full-handler C++ differential. The source is
`OpenMS4-core` revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.
The exact files, fixture hashes, source sections and literal expectations are in
[`mzml_sqlite_handler_provenance.json`](../tests/data/mzml_sqlite_handler_provenance.json).

## Public API inventory

`FORMAT/HANDLERS/MzMLSqliteHandler.h` is installed even though its class lives in
`OpenMS::Internal`. Its protected SQL helpers and stored counters are
implementation context, not additional public API targets. The public mapping
below is present in the native implementation. Metadata transport still has the
explicit native mzML subset boundary described below.

| Pinned C++ operation | Native operation | Contract |
| --- | --- | --- |
| Constructor `(filename, run_id)` | `new(path, run_id)` | Store path and writing run ID; source masks the sign bit. Reading uses the database run ID. |
| Implicit destructor | Rust ownership/drop | No persistent borrowed SQLite handle is exposed. |
| `setConfig` | `set_config` | Full metadata, lossy encoding, absolute m/z accuracy and batch size. |
| `setRunId` | `set_run_id` | Change subsequent writing run ID. |
| `getRunID` | `run_id` | Exactly one RUN row; both missing and multiple rows fail in the source. |
| `getNrSpectra`, `getNrChromatograms` | `nr_spectra`, `nr_chromatograms` | Count SQL records. |
| `readExperiment` | `read_experiment(meta_only)` | Full snapshot or SQL metadata projection, with optional primary data hydration. |
| `readSpectra`, `readChromatograms` | `read_spectra(ids, meta_only)`, `read_chromatograms(ids, meta_only)` | IDs are database record IDs, not positions in a returned vector. |
| `getSpectraIndicesbyRT` | `spectra_indices_by_rt(rt, delta, ids)` | Optional ID restriction; positive delta is an inclusive interval. Nonpositive source delta uses `RT >= target` and LIMIT 1 without ordering. |
| `createTables` | `create_tables` | Explicitly destructive database recreation. |
| `writeExperiment` | `write_experiment` | Run information, chromatograms and spectra. |
| `writeSpectra`, `writeChromatograms` | `write_spectra`, `write_chromatograms` | Append primary arrays and SQL metadata. |
| `writeRunLevelInformation` | `write_run_level_information` | RUN row and optional compressed mzML metadata snapshot. |

The current native defaults match lossy encoding, full metadata and m/z accuracy
0.0001, and initialize batch size to 500. The C++ constructor leaves batch size
uninitialized (CPP-190). Native selected
reads require unique nonnegative IDs and return ascending ID order. The default
`LossPolicy::Reject` checks whether the selected storage path can retain input
metadata; an explicit source-compatible loss policy permits documented losses.
These are native API policies, not source behavior inferred from tests.
Nonpositive finite accuracy, including the SqMassFile default −1 sentinel, uses
ordinary Numpress fixed-point estimation. Invalid configuration changes leave
the previous configuration intact.

## Native policy and resource boundaries

The public native additions are `HandlerConfig`, `HandlerLimits`, `config()`,
`set_loss_policy()` and the mutable `limits` field. SQLite errors retain their
cause through `Error::Io`; malformed values use checked errors, and unsupported
metadata uses `Error::Unsupported`. Source C++ exception class distinctions are
not preserved.

Each public read uses one SQLite read transaction. Missing files are opened
read-only and are not created. Record selection uses nonnegative 64-bit SQL IDs;
empty selected reads, duplicate requested IDs and missing IDs are errors. Positive
RT windows return IDs in ascending ID order. Nonpositive delta chooses the
earliest RT at or above the target, with ID as a tie-breaker, fixing the source's
unspecified LIMIT 1 ordering.

RT bounds (`rt - delta`, `rt + delta`, and `rt` itself for a nonpositive delta)
are bound as full binary `f64`. The source instead writes them into the SQL text
through its numeric formatter: 15 fractional digits with trailing zeros trimmed
for 0.01 <= |v| < 1e4, and the shortest round-trip scientific form otherwise. A
row whose retention time lies within that rounding of a bound can therefore be
selected differently, including in databases the source wrote. This is the same
deliberate precision correction the SWATH handler makes;
`small_rt_boundaries_use_full_binary_precision` pins it.

Each write operation commits atomically, and counters advance only after commit.
`create_tables` stages the seven-table schema and correct indexes in a sibling
file, then replaces the destination and resets counters. Existing SQLite sidecar
files prevent recreation. A newly constructed handler starts its append counters
at zero; reopening an already populated file does not resume appending.
Run-level information can precede record writes in a builder sequence. Such an
intermediate database can have an incomplete snapshot; a full read rejects
snapshot/SQL identity or count disagreement until the sequence is finished.
With Reject policy, each append checks its descriptive metadata against the
corresponding expected snapshot records, allowing rich metadata only when that
snapshot preserves it. A mismatch leaves both rows and writer IDs unchanged.
Individual low-level writes remain separate transactions. The source uses OpenMP
for spectrum/chromatogram array encoding; the native handler is serial. Batch
calls are recommended, as in the source header, with a normal chunk size of at
least 500 records; smaller positive chunks remain supported.

SQL values are bound, so record IDs and peptide strings containing punctuation
are literal data. Raw arrays use explicit little-endian f64 bytes; intensities
are checked when converted to native f32. Decoding requires one coordinate and
one intensity role of equal length, finite values, valid compression and correct
object ownership. Duplicate precursor/product rows and invalid activation values
are checked errors. Metadata-only reads avoid decoding primary arrays. Independent sampled-noise
grids carried in record metadata remain in a full snapshot; those grids are
separate from aligned auxiliary arrays.

Default limits allow one million combined records, ten million values per array,
twenty million combined values, 64 MiB per array BLOB or metadata snapshot,
512 MiB of logical allocation, a 2 GiB database and 500 million metered work
units. The allocation budget bounds what one operation materializes, not what
the file stores: stored payload is bounded by the database size limit, so a
database holding more than 512 MiB of array data still opens, counts and serves
metadata-only reads. These are configurable materialization and traversal ceilings, not a hard
process-memory or SQLite execution-time guarantee. Reads and writes check the
physical input file and logical page size; every write rechecks logical database
size and schema before committing, and recreation checks the staged database
before replacement. Fixed-schema scans reject views and virtual tables. A forged
ordinary table named `pragma_table_list` causes a checked query error instead of
bypassing these checks. A per-connection SQLite row-length limit supplements the
native byte checks; it applies to an entire row as well as individual fields.

Selections, native output records, metadata copies, array codec temporaries, and
snapshot record slots are charged before their allocations. Nested mzML parameter
and array allowances are reserved separately from the operation budget. The mzML
header registry retains its own fixed 256 MiB and 50 million work-unit limits;
these are separate from `HandlerLimits`, so the aggregate budget does not claim
to cover every nested parser allocation. SQLite's own allocations and query
planning are likewise outside that logical allocation counter.

## Source storage and metadata scope

Seven tables hold RUN, RUN_EXTRA, SPECTRUM, CHROMATOGRAM, PRECURSOR, PRODUCT and
DATA rows. The SQL metadata projection includes spectrum native ID, MS level,
retention time and polarity; precursor charge, target and offsets, drift time,
one activation method and energy, and optional peptide sequence; product target
and offsets; and chromatogram native ID plus precursor/product information.
Spectrum writing keeps only the first precursor and product and first activation
method. Chromatogram writing likewise stores one activation method. Spectrum
polarity is reduced to a positive/not-positive flag, so unknown polarity is
written as negative. Metadata numbers are rendered with stream precision 11.
The source schema has no auxiliary-array table.

When requested, RUN_EXTRA holds zlib-compressed mzML made from a metadata copy of
the experiment. The copy retains settings and per-record descriptive metadata,
but the source calls `clear(false)` on every spectrum and chromatogram before
serialization. That removes auxiliary arrays along with primary points. The
source header's promise of complete recovery therefore needs qualification
(CPP-201). The Rust full writer uses the documented [mzML metadata subset](MZML_SUPPORT.md).
It rejects peptide identifications and other metadata its serializer cannot
transport; `LossPolicy::Source` does not bypass those serializer checks.
Auxiliary arrays are rejected by default or deliberately dropped with Source
policy. SQL-only writes compare descriptive metadata with the SQL projection
and require Source policy when that projection would discard information.
A full-metadata flag alone is not evidence of full input preservation.

The C++ full reader uses the snapshot when configured and present, falling back
to SQL metadata if it is unavailable or the loaded experiment is empty. A
metadata-only read returns before primary data hydration. Selected reads use the
SQL projection rather than slicing the full metadata snapshot. Source SQL joins
can duplicate records for multiple precursor/product rows; missing
precursor/product rows can remove chromatograms from the projection. Native reads validate identity and cardinality rather than reproducing these
join-order and multiplicity failures. Full snapshot reads also check native IDs,
MS levels and retention times against SQL metadata.

## Compression and numeric comparisons

The source comments enumerate compression codes 0–7, but its actual decoder
accepts only 1, 5 and 6:

| Code | Implemented C++ decoding | Writer use |
| --- | --- | --- |
| 1 | zlib-compressed raw doubles | Lossless coordinates and intensities |
| 5 | zlib-compressed Numpress linear | Lossy m/z or chromatogram RT |
| 6 | zlib-compressed Numpress short logged float | Lossy intensities |

The native writer uses lossless code 1 for coordinate arrays shorter than three
values, even in lossy mode, to avoid an unusable linear fixed point. Nonempty
encoded arrays must have a finite positive fixed point. Linear encoding also
requires that its first two quantized coordinates fit unsigned 32-bit storage;
the raw source codec otherwise keeps only their low bits. Longer arrays with
unusable factors (for example three coordinates of 1e12) and negative initial
coordinates that would wrap fail atomically. Callers can explicitly select
lossless mode for those finite inputs. These are native corrections, not
additional source compression codes.

Other listed codes are rejected by the pinned decoder. Raw-double source I/O
uses host memory representation, so endian and byte-length behavior must be
explicit in the native implementation. Lossy m/z uses configured absolute
accuracy; chromatogram RT uses a fixed 0.05 seconds. The source disables Numpress
error checking during encoding. It does not safely validate paired array lengths
and distinct coordinate/intensity roles (CPP-191).

The paired retained fixtures contain two spectra (19,914 and 19,800 points) and
one chromatogram (48 points), with run ID 12345. Source full-data comparisons use
absolute OR relative acceptance, not both requirements simultaneously:

| Values | Absolute tolerance | Multiplicative relative tolerance |
| --- | --- | --- |
| Intensity | 0.0001 | 1.001 |
| m/z | 0.00001 | 1.000001 |
| Chromatogram RT | 0.05 seconds | 1.000001 |

Native reads carry the SQL run ID in `MSExperiment.sql_run_id`. The paired mzML
fixture contains a legacy `sqMassRunID` metadata key; the retained-file settings
comparison removes that key from the paired mzML expectation rather than
inventing a duplicate metadata owner on SQL reading.

Source tests also retain literal examples: spectrum 0 point 100 has m/z 204.817
and intensity 3857.86; chromatogram 0 point 20 has RT 0.200695 and intensity
147414.578125. Spectrum RTs are 0.2961 and 0.4738. These are source-test literals,
not output regenerated from C++ during this wave.

The source framework keeps numeric tolerances globally across sections. The
`cmpDataRT` calls leave absolute tolerance **0.05** and relative tolerance
**1.000001** active. Consequently, selected-spectrum RT assertions, subsequent
spectrum writer literals, and the first lossless chromatogram writer literals
inherit those values. The lossy chromatogram pass changes only the relative
value to **1.0002**; both its RT and intensity literals still allow absolute
error 0.05. The final reset occurs after these sections. The manifest retains
hashes of the test-framework definitions and this sequential tolerance ledger.
Stricter native assertions are additional checks, not the original thresholds.

## Class-test accounting and evidence

There are 13 START_SECTION blocks and 128 static assertion macros within them.
Eight additional assertion macros appear in three comparison helpers, for 136
static occurrences in the file. The helper checks execute in loops, and their
floating-point assertion macros run only after a similarity predicate fails;
these counts are not runtime assertion counts.

| Source section | In-section macros | Native test mapping |
| --- | ---: | --- |
| Constructor | 1 | `constructor_defaults_and_drop_do_not_create_a_file` |
| Destructor | 0 | `constructor_defaults_and_drop_do_not_create_a_file` |
| Run ID | 1 | `retained_run_id_and_counts` |
| Full and metadata-only experiment reading | 16 | `retained_experiment_metadata_and_all_numeric_values` |
| Spectrum count | 1 | `retained_run_id_and_counts` |
| Chromatogram count | 1 | `retained_run_id_and_counts` |
| Selected spectra | 21 | `selected_spectrum_class_test_literals` |
| Selected chromatograms | 15 | `selected_chromatogram_class_test_literals` |
| RT selection | 16 | `retention_time_class_test_literals_and_stable_nearest` |
| Experiment writing and recreation | 18 | `full_experiment_write_recreate_and_lossy_source_tolerances` |
| Loaded-path SQL injection regression | 2 | `run_level_loaded_file_path_is_bound_and_masked` |
| Spectrum append and recreation | 14 | `repeated_spectrum_appends_preserve_source_literals` |
| Chromatogram append, recreation and compression | 22 | `repeated_chromatogram_appends_lossless_and_lossy` |

`setConfig`, `setRunId`, `createTables` and run-level writing are additionally
covered by configuration atomicity, run-ID masking, destructive-recreation,
metadata-builder and rollback regressions. The source tolerance-state ledger is
applied separately to lossless and lossy chromatogram literals.

The mzML and sqMass files are exact retained upstream inputs. The existing sqMass
copy in `sqlite_s1_source_review` is referenced rather than duplicated. No handler
C++ class test or full SDK differential has been executed in this wave. Native
test names are mapped in the table. The separate exact raw Numpress probe ran 42
cases on kim and supports CPP-218's short-coordinate defect; it did not execute
MSNumpressCoder or this SQL handler. All probe sources, binaries and logs are
retained outside this repository and hashed in the manifest. Independent SQL
observations in the separate S1 review manifest are not an executed C++ handler
oracle.

## C++ findings

The shared [issue log](../OpenMS_CPP_ISSUES.md) owns stable entries. Relevant
handler findings are CPP-190 (uninitialized batching), CPP-191 (array pairing),
CPP-192 (hydration order), CPP-193 (unescaped record metadata), CPP-194 (partial
writes and counters), CPP-195 (wrong index tables), CPP-198 (recreation counters),
CPP-199 (negative activation value), CPP-201 (auxiliary-array loss), CPP-202
(reads create missing files), CPP-218 (unusable linear fixed points), and CPP-220
(negative initial coordinate quantization). Native handling and focused remote
regressions are recorded in the provenance manifest. The long-array extension to
CPP-218 and the negative-coordinate finding remain source-reviewed C++ cases;
the earlier 42-case raw-codec probe does not execute those triggers.
The loaded-path injection source regression tests an existing upstream fix and
must not be reported as an unpatched C++ path-injection defect.

## Integration review still in progress

The accompanying [precursor-activation update](MZML_PRECURSOR_ACTIVATION_SUPPORT.md)
handles source-supported collision-energy and supplemental-activation metadata
in RUN_EXTRA. That inherited mzML correction has its own source mapping and
focused tests. Combined feature-graph, full-suite and external-review status is
owned by the integrating wave record; passing these handler-focused checks is
not a claim of complete SDK or TOPP-tool validation.
