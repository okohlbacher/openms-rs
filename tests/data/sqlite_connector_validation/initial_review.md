**Verdict:** no blocking correctness bugs found in the connector or the 17 tests. The three open modes and five methods match the pinned source where intended and diverge only where documented. Findings below are ordered by severity; none require new abstractions or scope growth.

**1. `column_exists` decodes column names as `String` (minor fidelity gap).**
`row.get::<_, String>(1)` fails with a UTF-8 conversion error if a column name stored in the schema is not valid UTF-8. The source does a byte `strcmp` and would simply return `false`. Since `column_name: &str` can never equal non-UTF-8 bytes, the correct result is `false`, not `Err(Io)`. One-line fix: compare `row.get_ref(1)?.as_bytes()?` against `column_name.as_bytes()`. Also avoids an allocation per row. Reachable only via schemas created outside this API, so low priority.

**2. `raw_execute` behavior depends on rusqlite's `extra_check` feature (verify Cargo features).**
Without `extra_check`, `raw_execute` steps once and reports `ExecuteReturnedResults` on a row, which is exactly the source's `SQLITE_DONE` check and matches the doc sentence "a statement producing rows may already have effects before it is rejected". With `extra_check` enabled, rusqlite rejects before stepping any statement with `column_count() > 0` or where `sqlite3_stmt_readonly` is true, which includes `BEGIN`, `COMMIT`, and flag pragmas. That would break `execute_bind_statement("BEGIN", &[])` and falsify the doc. Confirm the feature is off; no code change needed if so.

**3. `destruction_closes_the_connection` cannot actually detect a leaked handle on Unix (test strength, queued test).**
The abandoned `BEGIN; INSERT` holds only a RESERVED lock. A leaked first connection would still let the reopened connection read `T` and count 0, and `remove_file` succeeds on open handles on macOS/Linux. Add one write on the reopened `ReadWrite` connection, for example `INSERT INTO T VALUES(2)`; it fails with `SQLITE_BUSY` if the first connection was not closed. The rollback-on-close assertion (count 0 after CREATE was autocommitted) is correct as written.

**4. `locked_query_errors_are_not_reported_as_absence` (queued test) relies on two SQLite properties; both hold with rusqlite's bundled SQLite.**
- `PRAGMA table_info` emits `sqlite3CodeVerifyNamedSchema`, so `sqlite3_step` opens a read transaction and hits `SQLITE_BUSY` under the locker's EXCLUSIVE lock even though the reader's schema is warm. If a future SQLite dropped that, the assertion would see `Ok(true)`. It has been present since 2011; expect pass.
- The file is in default DELETE journal mode. In WAL mode readers would not be blocked and all three `unwrap_err()` calls would fail. Nothing in the test changes journal mode, so this is a note, not a bug.
The warm-up is well placed: without it, prepare rather than step would fail, which still yields `Io` but for the wrong reason.

**5. Doc precision, `open`.** "SQLite special filenames are passed through" reads as including `file:` URIs. `SQLITE_OPEN_URI` is not set (same as the source), so URIs are treated as literal paths unless SQLite was built with `SQLITE_USE_URI`. Say ":memory: and the empty name".

**Verified as correct (details reviewers may want):**
- `Batch` + `raw_query` drain reproduces `sqlite3_exec` ordering: each statement is finalized before the next is prepared, so `CREATE TABLE T; INSERT INTO T` works and a failing later statement leaves prior effects. The post-first-row overflow test is a valid probe because `ORDER BY rowid` emits row 1 before evaluating row 2.
- `execute_bind_statement` prepares only from offset 0 and never calls `batch.next()` again, so the tail is never handed to `sqlite3_prepare`. The `foreign_keys=OFF` tail test is a real detector: flag pragmas mutate `db->flags` at compile time.
- Empty `&[]` binds via rusqlite's `zeroblob(0)` path, giving `typeof = 'blob'`, matching `sqlite3_bind_blob` with a non-null `c_str()` and length 0.
- Over-binding yields `SQLITE_RANGE` before any step, so no partial effects; unbound trailing parameters stay NULL.
- Busy timeout: rusqlite installs a 5 s handler on open; `busy_timeout(Duration::ZERO)` clears it, restoring the source's immediate `SQLITE_BUSY`.
- Leading `;` or comments before the first statement are skipped by SQLite's own parser in both implementations, so `Batch`'s null-statement skip introduces no divergence. Comment-only input maps the source's `SQLITE_MISUSE` throw to `InvalidValue`, which the doc states.
- `query_row` calls `check_no_tail`; the two literal queries carry no trailing `;`, so no spurious `MultipleStatement` error.

**Source quirks vs intended native differences (all documented in-code):**
- Quirk, not ported: `countTableRows` and `columnExists` interpolate the name raw, allowing `schema.table` or injection; `tableExists` breaks on `'`. Rust quotes or binds. Known OpenMS callers pass plain names.
- Quirk, not ported: `countTableRows` ignores the `sqlite3_step` return code; source header promises `SqlOperationFailed` for a missing table but `prepareStatement` throws `IllegalArgument`. Rust returns `Io`.
- Quirk, not ported: `executeBindStatement` leaks the statement on bind/step failure. Rust finalizes via `Drop`.
- Intended difference: embedded NUL is rejected whole rather than silently truncated at `c_str()`.
- Intended difference: `columnExists` is documented to require an existing table; both implementations actually return `false`.

Not reviewed: `docs/SQLITE_CONNECTOR_SUPPORT.md`, `src/format/mod.rs` gating, and `Cargo.toml` features were not supplied.
