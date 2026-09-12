# SQLite connector

Source pin: `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.
This S0 increment covers the installed `FORMAT/SqliteConnector.h` interface:
three open modes, connection ownership, table and column existence, row counts,
SQL execution and byte-string binding. The implementation is
`src/format/sqlite_connector.rs`, behind the optional `sqlite` feature.

All public operations have native counterparts. The source, native implementation
and 19 native tests have been reviewed. Remote build and test results are recorded
separately in [`VALIDATION.md`](VALIDATION.md); a mapped operation is not by itself
evidence that a test was executed.

The pinned dependency is `rusqlite = 0.40.2` with bundled SQLite. This is a
native Rust wrapper around SQLite, with no OpenMS C++ library or SDK dependency.
Bundled SQLite is C and needs a C toolchain when this feature is enabled; the
feature is outside the crate defaults.

## Public API mapping

| Pinned C++ member | Native counterpart | Contract |
|---|---|---|
| `SqlOpenMode::READ_ONLY` | `SqlOpenMode::ReadOnly` | Existing database, no writes. |
| Conditional `SqlOpenMode::READONLY` alias | `SqlOpenMode::ReadOnly` | A source spelling alias, not another mode. |
| `SqlOpenMode::READWRITE` | `SqlOpenMode::ReadWrite` | Existing database, reads and writes. |
| `SqlOpenMode::READWRITE_OR_CREATE` | `SqlOpenMode::ReadWriteOrCreate` | Creates a missing database. The native mode default. |
| Deleted `SqliteConnector()` | No connector `Default` | Opening requires a path. |
| `SqliteConnector(filename, mode = READWRITE_OR_CREATE)` | `SqliteConnector::new(path)` and `open(path, mode)` | `new` selects the source default mode. |
| `~SqliteConnector()` | Owned connection and Rust drop | No copying of an owned native handle. |
| `tableExists(tablename)` | `table_exists(table)` | A fallible table-existence query. |
| `columnExists(tablename, colname)` | `column_exists(table, column)` | A fallible column-existence query. |
| `countTableRows(table_name)` | `count_table_rows(table)` | Returns a checked `usize` count; an unknown table fails. |
| `executeStatement(statement)` | `execute_statement(sql)` | Executes SQL, including a batch and caller-controlled transactions. |
| `executeBindStatement(prepare_statement, data)` | `execute_bind_statement(sql, data: &[&[u8]])` | Binds byte strings as BLOB values, using positions starting at one. |

The protected `openDatabase_` and `void* db_` are implementation details, not
public native extension points. `Internal::SqliteConnectorFriend` is a forward
declaration granting the private implementation access; no native handle is
exposed by this port.

## Pinned source behavior

The constructor calls `sqlite3_open_v2` with exactly the selected mode and
throws `SqlOperationFailed` on failure. `executeStatement` calls `sqlite3_exec`
with no row callback, so it accepts multiple statements and discards rows from
queries. It does not begin an implicit transaction: if a later statement fails,
earlier statements can already have changed the database. Transaction control
belongs to the SQL caller.

The binding operation prepares the first SQL statement and ignores its tail.
Each `std::string` is bound with `sqlite3_bind_blob` and its explicit length;
these are binary values, even when their bytes spell text. Unbound parameters
remain SQL `NULL`; excess supplied bindings cause an error. It steps once and
requires `SQLITE_DONE`, so it is intended for statements that do not return rows.
The source does not verify a read-back value in its binding test: that section
checks only that the row count becomes one.

`tableExists` asks `sqlite_master` for an exact table name and `type='table'`.
It does not search temporary tables or include views. `columnExists` iterates
`PRAGMA table_info` and compares column names with `strcmp`; no column matches
on a missing table. The public header says that the table must exist, but the
implementation returns false for that case.

The header says an unknown table makes `countTableRows` throw
`SqlOperationFailed`. Its actual prepare helper throws `IllegalArgument`, and
the source class test explicitly expects `IllegalArgument` for `UNKNOWN`.
The native error category is documented against both below rather than implying
that these source declarations agree.

## Native differences and operational boundary

The native connection owns its resources through `rusqlite`; it offers no
`Clone`, raw connection accessor or statement pointer. Table names are treated
as literal identifiers and quoted, and query values are bound. Thus an apostrophe
or SQL punctuation in a table name does not become executable query syntax.
Native query results check the SQLite operation outcome before inspecting rows.

The integrating S0 design preserves raw SQL execution and caller-controlled
transactions. It intentionally adds no wrapper work, memory or time budgets and
no SQL parser, implicit transaction, progress callback or query planner. A short
SQL string can perform substantial work; SQL byte length would not bound that
work. Higher-level format adapters must define their own input, schema and
transaction policies. SQLite's own engine constraints still apply.

The chosen native error mapping wraps SQLite operational failures in
`Error::Io`, preserving the underlying `rusqlite` error as an error source.
Explicit native validation failures use `Error::InvalidValue`. This does not
reproduce the distinction between C++ `IllegalArgument` and
`SqlOperationFailed`; callers use the crate error and its underlying cause.

The binding operation prepares and executes **only the first statement**.
Trailing SQL is ignored without preparing it, including malformed text or a
`PRAGMA` that could change connection state during preparation. This preserves
the pinned source behavior. Empty or comment-only binding input is a checked
`Error::InvalidValue`. A statement that returns a row is an SQLite operation
error, but it may already have effects; this wrapper does not roll them back.

Both SQL execution methods reject embedded NUL in the whole input before
executing any prefix. The source passes SQL through C strings and can execute
only the prefix. This check applies to SQL text, not BLOB data: embedded NUL,
non-UTF-8 bytes and zero-length BLOBs are preserved. Missing bindings remain
NULL; excess bindings and constraint failures are errors.

`execute_statement` consumes all result rows of every statement, so an error
while evaluating a later row is reported. It uses `rusqlite::Batch` and row
iteration rather than a convenience call that stops after the first row.
The busy timeout is explicitly zero, as with the source connection, so lock
conflicts return an error without an added retry delay. The constructor also
passes through SQLite special filenames such as `:memory:` and an empty
filename; the latter creates a private temporary database.

## Source test accounting and evidence

`SqliteConnector_test.cpp` is registered in
`src/tests/class_tests/openms/executables.cmake:290`. It deliberately includes
only the public header. All data is constructed through SQL in a temporary
database; there is no upstream input file or retained expected database to copy.

| Source section | Assertion macros | Literal or checked outcome | Native test |
|---|---:|---|---|
| Constructor | 2 | Non-null connector; `liveness_check` exists after `CREATE TABLE IF NOT EXISTS`. | `construction_creates_a_usable_database` |
| Destructor | 0 | Deletes the connector; contains no assertion. | `destruction_closes_the_connection` |
| `executeStatement` | 2 | Table `T` exists after insert `(1, 'a')`; `THIS IS NOT SQL` throws `IllegalArgument`. | `execute_statement_inserts_and_reports_malformed_sql` |
| `tableExists` | 2 | `T` true; `DOES_NOT_EXIST` false. | `table_exists_matches_exact_main_table_names` |
| `columnExists` | 3 | `ID` true; `NAME` true; `MISSING` false. | `column_exists_matches_exact_column_names` |
| `countTableRows` | 3 | Empty `T` has 0 rows; inserts `(1, 'a')`, `(2, 'b')` give 2; `UNKNOWN` throws `IllegalArgument`. | `count_table_rows_preserves_the_source_zero_and_two_counts` |
| `executeBindStatement` | 1 | Bind `bound_value` at `?1`; count becomes 1. | `execute_bind_statement_inserts_the_source_bound_value` |
| Open modes | 4 | Read-only `RO` exists and has 3 rows; missing file fails under `READ_ONLY` and `READWRITE`. | `open_modes_preserve_read_only_and_missing_file_behavior` |

There are **eight sections and 17 assertion macros**, not eight sections with
assertions. The constructor's pointer assertion maps to successful native
construction, and deletion maps to ownership/drop rather than pointer equality.
Neither substitution proves every operation or failure path.

Evidence is **tier 3, source review** for these transcribed outcomes. Native-only
boundary and resource-lifetime tests are independently derived, tier 4. No
full OpenMS C++ SDK was built or executed for S0. A separate adapted C++ probe
was executed and retained, as described below; it does not execute the upstream
class tests. Running SQLite through `rusqlite` alone is not an executed OpenMS
C++ differential.
Execution results belong in the integration validation record after the remote
gates finish. The manifest
[`sqlite_connector_provenance.json`](../tests/data/sqlite_connector_provenance.json)
pins the reviewed sources and section accounting.

## C++ findings

The integrating agent owns `OpenMS_CPP_ISSUES.md`. S0 source-reviewed candidates
are recorded as CPP-181 through CPP-187, plus the CPP-189 exception-documentation
mismatch. CPP-184, CPP-185 and CPP-189 have the adapted execution evidence below;
the other entries remain source-reviewed candidates.

| Issue | Source concern | S0 boundary |
|---|---|---|
| CPP-181 | Implicit copies duplicate the owned database pointer. | Owned private connection, no `Clone`. |
| CPP-182 | Failed-open constructor throws without closing a returned handle. | `rusqlite` owns opening and cleanup; no allocation-leak measurement claimed. |
| CPP-183 | Bind/step error branches bypass statement finalization. | RAII statement finalization; bind, constraint and returned-row errors have regression tests. |
| CPP-184 | Table-name helpers interpolate SQL syntax. | Bound table-name value and quoted identifiers; punctuation and expression cases tested. |
| CPP-185 | Query helpers ignore step errors and can leak on error. | Checked query/row outcomes; warmed-schema exclusive-lock regression covers all three helpers. |
| CPP-186 | Private string extraction truncates at embedded NUL. | Private extraction helpers are outside S0's public API. |
| CPP-187 | Private integer-to-string extraction uses a 32-bit read. | Private extraction helpers are outside S0's public API. |
| CPP-189 | Unknown-table exception documented as `SqlOperationFailed`, implemented and tested as `IllegalArgument`. | Native `Error::Io` retains the SQLite cause; both source contracts are documented. |

## Adapted C++ execution

A bounded probe compiled the exact pinned `SqliteConnector.cpp`, its exact
public/private headers and the pinned SQLite declaration header on `kim`, with
GCC 13.3 and the host's SQLite 3.45.1 shared library. The full SDK was not linked.
Three support headers were substituted: fixed-width type aliases/export macros,
the two exception classes, and numeric `StringUtils::toStr` used in diagnostics.
The exception adapters preserve the two distinct throw/catch paths and messages,
not the OpenMS exception ABI, inheritance hierarchy or file/line metadata.

These are **tier 2 adapted/isolated observations**. Sources, adapters, driver,
compiler/output logs, executable and the adapter manifest are hashed under
`external_reference_artifacts` in the provenance manifest and retained outside
this repository in `../oracle/sqlite-connector-s0-probe/`.

| Case | Observed pinned C++ behavior through adapters | Native regression |
|---|---|---|
| Actual table name containing an apostrophe | `tableExists` throws during SQL preparation. | `quoted_identifiers_are_literal_and_cannot_change_the_query` checks literal punctuation-bearing names. |
| Nonexistent name containing an `OR 1=1` expression | `tableExists` returns true. | The native expression-bearing name returns false. |
| Actual table name containing a space | `countTableRows` and `columnExists` throw during preparation. | Literal identifier quoting is checked. |
| Warm schema, second connection holds `BEGIN EXCLUSIVE` | Direct prepare returns `SQLITE_OK`; step returns `SQLITE_BUSY`. Both existence helpers nevertheless return false. | `locked_query_errors_are_not_reported_as_absence` requires errors from all three queries. |
| Row count under the same lock | Throws `SqlOperationFailed` and leaves one outstanding statement, counted through `sqlite3_next_stmt`. After rollback, a new count returns 3. | The native lock test checks error then recovery to count 3; it does not measure C++ or native heap allocations. |
| `countTableRows("UNKNOWN")` | Throws adapted `IllegalArgument`. | Native count test requires an `Error::Io` carrying `rusqlite::Error`. |

The eleven further native tests cover binary type/content and NULL defaults,
literal identifiers, batch failure effects and caller rollback/commit, errors
after a first result row, failed binding/constraint recovery, ignored binding
tails without preparing a `PRAGMA`, NUL SQL rejection, lock errors and SQLite's
temporary filename, externally created non-UTF-8 schema names, and bound
transaction/RETURNING behavior even when a downstream dependency enables
`rusqlite/extra_check`. The source destructor section is strengthened by dropping
an open transaction and checking that reopening sees zero uncommitted rows.
These native tests do not reproduce private `extractString` or
`extractValueIntStr` helpers.

## Scope after S0

`SqliteConnector_impl.h` is explicitly non-installed, as recorded in
`FORMAT/sources.cmake:166-175`. `SqliteHelper` statement functions,
`clearSignBit`, `SqlState`, extraction templates and throwing getters are
private implementation context. Their availability is not claimed through a
public Rust helper API. Later format modules can use safe SQLite primitives
internally as their actual requirements become known.

This connector alone does not implement `MzMLSqliteHandler`,
`MzMLSqliteSwathHandler`, `SqMassFile`, the SQLite consumer/access adapters,
`OSWFile`, `OMSFileStore`, `OMSFileLoad`, `OMSFile`, their schemas or
`FileHandler` dispatch. Those remain the subsequent SQLite stages in
[`PORTING_WAVES.md`](PORTING_WAVES.md).
