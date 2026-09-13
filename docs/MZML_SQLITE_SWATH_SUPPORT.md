# sqMass SWATH lookups

`format::mzml_sqlite_swath_handler` implements every public operation of the installed `FORMAT/HANDLERS/MzMLSqliteSwathHandler.h` at OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The module requires `sqlite`; it does not require mzML, Numpress, or the full sqMass reader/writer. Upstream places this helper in `Internal` and describes it as outside its stable public API. Its installed header remains part of the SDK port scope.

## API mapping

| Source member | Rust counterpart | Contract |
|---|---|---|
| constructor(filename) | `MzMLSqliteSwathHandler::new(path)` | Remembers the path without opening the file. |
| implicit destructor | Rust ownership | The handler owns a path; each query releases its connection and transaction before returning. |
| `readSwathWindows()` | `read_swath_windows()` | Distinct `(center, lower, upper)` tuples from precursor rows associated with MS2 spectrum IDs. |
| `readMS1Spectra()` | `read_ms1_spectra()` | All spectrum IDs whose MS level equals 1. |
| `readSpectraForWindow(map)` | `read_spectra_for_window(&window)` | IDs of precursor rows with isolation targets inside inclusive `center ± 0.01` m/z. |
| native extension | `with_limits(path, SwathReadLimits)` | Configure the maximum scanned rows; `new` uses 1,000,000. |

There are no other explicit public methods in the source header. Protected path/counter members are implementation details. The source counters are unused by these read operations.

`SwathWindow` exposes only the three fields actually populated by discovery: `center`, `lower`, and `upper`, all `f64` m/z values. `lower` and `upper` are absolute boundaries, calculated from the precursor offsets. This is not a full port of `OpenSwath::SwathMap`: the source leaves `ms1=false`, ion-mobility limits at `-1`, and its spectrum-access pointer empty. Those unpopulated fields and the spectrum-access association have no counterparts in this value. The original class test's `ms1 == false` assertion is accounted for by this documented output mapping; no Rust boolean assertion or full SwathMap coverage is claimed.

## Preserved behavior

Window discovery associates PRECURSOR.SPECTRUM_ID with SPECTRUM.ID at MS level 2. Repeated identical boundaries collapse; equal centers with differing boundaries remain distinct. Finite negative offsets retain the source arithmetic and can therefore produce inverted boundaries. Negative and positive zero compare equal for deduplication.

Center lookup uses only the supplied center. It ignores lower/upper even if they are nonfinite, does not filter MS level, does not require that a precursor ID appears in SPECTRUM, and retains repeated IDs when several matching precursor rows exist. A NULL isolation target does not match. These details follow the source PRECURSOR-only query. MS1 lookup excludes NULL MS levels.

Results are database record IDs, not ordinal spectrum positions. Rust retains nonnegative signed 64-bit SQLite IDs, including values above the source's 32-bit `int` range. All operations are serial, as are these upstream queries.

## Native corrections and boundaries

All accessors open the database explicitly read-only. A missing path fails without creating a file, correcting [CPP-202](../OpenMS_CPP_ISSUES.md#cpp-202--read-accessors-create-an-empty-database-when-the-input-path-is-missing). Each operation uses a read transaction so schema validation and all row scans share a snapshot. Statements and transactions close on success and error; no partial vector is returned.

Fallible row iteration reports SQLite errors even after a valid row. In the C++ loops, unchecked `sqlite3_step` can turn a late evaluation error into a successful prefix ([CPP-203](../OpenMS_CPP_ISSUES.md#cpp-203--swath-accessors-can-return-partial-success-after-sqlite-step-errors)). The native regression uses a generated MSLEVEL column with `abs(i64::MIN)` to trigger an actual later SQLite error. The C++ consequence is source-reviewed; that trigger has not been executed through C++ here.

Center lookup skips precursor rows whose spectrum ID is NULL, including chromatogram precursors. C++ incorrectly treats such a row as the end of the entire query ([CPP-196](../OpenMS_CPP_ISSUES.md#cpp-196--swath-selection-stops-at-a-matching-chromatogram-precursor-null-id)). Window discovery documentation describes distinct tuples rather than the source header's inaccurate distinct-center promise ([CPP-197](../OpenMS_CPP_ISSUES.md#cpp-197--swath-window-docs-promise-distinct-centers-but-query-deduplicates-full-bounds-tuples)).

Native center endpoints retain full binary f64 precision. The source serializes SQL endpoints through its numeric formatter (15 fractional decimal digits in the fixed-format range), so extremely small boundary differences can select a different row. A regression at center 0.12345678901234567 records this intentional precision difference; the adapted C++ probe does not validate original formatter rounding.

Source query order is unspecified. Rust sorts IDs ascending and windows lexicographically by center, lower, upper. Scans use the ordinary main SPECTRUM/PRECURSOR tables, then bounded Rust membership, deduplication, and sorting. Views and virtual tables are rejected before evaluation. The schema check invokes the table-valued pragma explicitly; a database table shadowing its name causes a SQL error rather than authorizing a forged view. Discovery requires both tables; MS1 lookup only SPECTRUM; center lookup only PRECURSOR. A schema without required columns fails. This adapter does not claim to validate the complete seven-table sqMass schema.

`SwathReadLimits::max_records` bounds the total source rows visited by one operation, including nonmatching and chromatogram rows. Discovery shares one budget across both tables. The other operations scan one table. A zero ceiling accepts only empty scanned tables. At most one additional row is stepped to detect a ceiling violation, before its fields are read or retained. Memory is O(limit), and native sorting is O(limit log limit); no SQL join or sort is performed before the ceiling check. SQLite page I/O, schema processing, and expressions in externally defined generated columns are not a wall-clock instruction budget.

All consumed SQL IDs and levels must have integer storage types, and consumed coordinates must be numeric; matching IDs must be nonnegative. Discovery rejects negative precursor IDs even when they would not join. Required discovery coordinates cannot be NULL. Center lookup skips NULL spectrum IDs before reading their other columns, and skips NULL targets. Nonfinite non-NULL targets on spectrum precursor rows are rejected even when they are outside the requested window; discovery checks coordinates of matching MS2 precursors. Supplied center and arithmetic-derived boundaries must be finite. SQLite can convert an inserted NaN to NULL, so such database input follows the stated NULL rules rather than recovering an erased NaN payload.

Underlying SQLite open/prepare/step/type failures map to `Error::Io`, retaining the rusqlite cause. Missing or nonordinary required tables, negative IDs, nonfinite values, and record-limit violations map to `Error::InvalidValue`. These error and validation policies intentionally replace C++'s unchecked step status, NULL sentinels, numeric coercions and narrowing.

## Evidence and tests

[The provenance manifest](../tests/data/mzml_sqlite_swath_provenance.json) pins the header, implementation, class test, header registration, output-value source, and original input hashes. All five source class-test sections are mapped: construction, destruction, discovery, MS1 lookup, and center lookup.

The C++ class test first loads `SwathFile.mzML` and writes a database through SqMassFile. The native query test instead inserts a reviewed projection of all 114 spectra's IDs, MS levels, and isolation center/offsets from that original mzML into narrow SQLite tables. The complete projection is retained in the manifest. This tests the SWATH query layer independently of the full sqMass writer; it is not an executed C++ output fixture or evidence of writer parity.

The original expectations are five windows, 19 MS1 spectra (IDs 0 through 108), and 19 spectra in each tested window (first IDs 1/109, second 2/110). The first window is 400–425 m/z centered at 412.5, the second has bounds 425–450, and the last 500–525. Tests compare these exactly: the source permits relative tolerance 1.0005, but these half-integer boundaries are exactly representable. A separate retained `SqliteMassFile_1.sqMass` input tests its two MS1 spectra and absence of SWATH spectra without regenerating the file.

Native tests cover NULL precursor rows, duplicate and orphan matches, distinct equal-center/different-width windows, inclusive tolerance endpoints, insertion-independent ordering, signed 64-bit IDs, invalid numeric types and values, zero/exact/exceeded shared budgets, schema/view rejection, missing paths, invalid database bytes, lock errors and connection cleanup. These are tier 4 invariants. Source literals and projected original input are tier 3 evidence. Executed C++ probes and remote validation are recorded separately when available; no unexecuted command is counted as passing evidence.

An adapted C++ probe executed on kim with g++ 13.3.0/C++20 and host SQLite 3.45.1 reproduces CPP-196/197/202. It compiles the exact pinned SWATH and connector translation units/headers, with isolated substitutes for OpenMS exception/types/string helpers and a pointer-only spectrum-access dependency. With a matching chromatogram precursor first it returns no IDs; in the middle it returns only the preceding spectrum; without it the control returns both spectra. Equal centers with differing offsets return both windows. A missing path remains absent after construction but a failed read creates a zero-byte file. Exact build/run logs, driver, adapters, binary, and input/output hashes are listed under `external_reference_artifacts` in the provenance manifest and remain outside this repository under `../oracle/sqlite-swath-s1-probe/`. This is tier 2 adapted-probe evidence, not a full SDK, full sqMass handler, or upstream class-test execution. Its string-formatting adapter also limits what can be inferred about arbitrary tolerance endpoint rounding.

Native validation on kim passed all 18 tests with current Rust and Rust 1.85.0 under the sqlite-only feature selection. The corresponding clippy gate and rustdoc build passed with warnings denied. The requested independent Fable 5.1 review has not executed: remote OAuth is expired and the local review packet awaits explicit approval after automatic approval review rejected the export. No Fable review conclusion is claimed.
