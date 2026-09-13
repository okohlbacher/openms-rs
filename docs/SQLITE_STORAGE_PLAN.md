# S1 sqMass handler plan (read-only preparation)

Prepared 2026-09-13 against OpenMS4-core bc9cc12514c768385ce121d6ca4bb710fe1983c4, after Rust checkpoint a477773. This historical preparation did not implement Rust or execute C++. The handler and SWATH ports were subsequently integrated; see [current waves](PORTING_WAVES.md), [handler support](MZML_SQLITE_HANDLER_SUPPORT.md) and [SWATH support](MZML_SQLITE_SWATH_SUPPORT.md). Source paths below are relative to `.reference/openms4-core-bc9cc12`; source hashes and issue details are in [`the S1 source-review manifest`](../tests/data/sqlite_s1_review.json).

## Minimum dependency closure

Both handler headers are registered in `FORMAT/HANDLERS/sources.cmake` and installed/exported despite their `Internal` namespace and recommendation to prefer SqMassFile. That stability warning is not an SDK scope exclusion. All read/count/lookup paths should explicitly open read-only; the source default-creating connection mutates missing input paths before failing a table query.

Reuse `MSExperiment` (including `sql_run_id`), `MSSpectrum`, `MSChromatogram`, precursor/product/settings types, `format::numpress` raw codecs, `MSNumpressCoder` configuration/status checks, and the mzML reader/writer. No replacement kernel types or general SQL abstraction are needed. S0's `SqliteConnector` owns a private rusqlite connection. S1 needs one narrow crate-private accessor for typed rows, nullable values, mixed TEXT/INTEGER/REAL/BLOB parameters, reused statements, and a transaction. Do not broaden S0's public surface.

Full RUN_EXTRA support requires mzML and zlib as well as SQLite. Prefer a small `sqmass = ["sqlite", "mzml"]` feature (integration decision); alternatively gate the handler on `all(sqlite,mzml)`. Keep the SQLite connector independently usable. SWATH queries require only sqlite and a value record. There is no existing native `SwathMap` identified: introduce the minimal window value (`center`, `lower`, `upper`, source default fields if exposing a SwathMap mapping), without pulling in the OpenSwath spectrum-access interface just to transport three numbers. Inventory the upstream SwathMap header in provenance if any fields are mapped.

Use existing bounded compression internals through a small crate-private shared helper; avoid base64 encode/decode round trips and duplicating MSNumpressCoder. Preserve its `Encoded`/`EmptyInput`/`Disabled`/`Rejected` distinctions: a rejected encoding is an error, never an empty blob advertised as valid data. Apply checked arithmetic, aggregate decoded-byte/array/record budgets, and decompression ceilings before allocations. Reuse the current format error conventions.

## Public operation mapping

| C++ API | Native shape / behavior to retain |
|---|---|
| constructor(filename, run_id) | path-owning handler, explicit initialized defaults; source masks bit 63 of write run ID |
| setConfig | full metadata bool, lossy bool, linear absolute m/z accuracy default 0.0001 Th, batch limit default 500; validate finite positive accuracy and positive batch |
| setRunId | update masked write ID; reads derive actual RUN.ID |
| getRunID | exactly one RUN record required, checked integer conversion |
| getNrSpectra / getNrChromatograms | exact checked row counts |
| readExperiment(meta_only) | complete experiment or metadata only; RUN_EXTRA when enabled and present, otherwise table reconstruction; validate snapshot/table consistency |
| readSpectra / readChromatograms(indices,meta_only) | selected SQL record IDs, not ordinal positions; missing/negative/duplicate IDs rejected per source behavior; explicit deterministic result order |
| getSpectraIndicesbyRT | inclusive RT +/- delta when delta > 0; source delta <= 0 selects RT >= query, at most one; specify stable RT,ID ordering for nearest-at-or-after case |
| createTables | explicit destructive operation, correct seven-table schema/indexes, reset counters after success; fail before writes if existing file cannot be removed |
| writeExperiment | run + optional snapshot + chromatograms + spectra, one atomic public operation |
| writeSpectra / writeChromatograms | append with checked counters, prepared bound writes in bounded chunks; whole call atomic and counters advance only after commit |
| writeRunLevelInformation | parameterized RUN and optional RUN_EXTRA inserts atomically; docs explain this is low-level |
| readSwathWindows | distinct center/lower/upper tuples from MS2 precursors, source defaults for remaining fields |
| readMS1Spectra | MSLEVEL == 1 record IDs |
| readSpectraForWindow | center inclusive +/- 0.01 m/z; lower/upper ignored, exclude chromatogram NULL IDs; source does not restrict MS level |

Do not port protected template/query helpers as public APIs. Preserve stable source literal results while documenting deterministic ordering where C++ only offers an unspecified natural order. Do not guess that vector index equals a database ID. Keep append-to-existing-file unsupported unless deliberately implemented and tested; upstream explicitly documents that constructor counters start at zero.

## Exact storage schema and encoding

Seven tables, no schema version or foreign-key constraints in source:

| Table | Columns |
|---|---|
| DATA | SPECTRUM_ID INT, CHROMATOGRAM_ID INT, COMPRESSION INT, DATA_TYPE INT, DATA BLOB NOT NULL |
| SPECTRUM | ID INT PRIMARY KEY NOT NULL, RUN_ID INT, MSLEVEL INT, RETENTION_TIME REAL, SCAN_POLARITY INT, NATIVE_ID TEXT NOT NULL |
| RUN | ID INT PRIMARY KEY NOT NULL, FILENAME TEXT NOT NULL, NATIVE_ID TEXT NOT NULL |
| RUN_EXTRA | RUN_ID INT, DATA BLOB NOT NULL |
| CHROMATOGRAM | ID INT PRIMARY KEY NOT NULL, RUN_ID INT, NATIVE_ID TEXT NOT NULL |
| PRODUCT | SPECTRUM_ID INT, CHROMATOGRAM_ID INT, CHARGE INT, ISOLATION_TARGET REAL, ISOLATION_LOWER REAL, ISOLATION_UPPER REAL |
| PRECURSOR | SPECTRUM_ID INT, CHROMATOGRAM_ID INT, CHARGE INT, PEPTIDE_SEQUENCE TEXT, DRIFT_TIME REAL, ACTIVATION_METHOD INT, ACTIVATION_ENERGY REAL, ISOLATION_TARGET REAL, ISOLATION_LOWER REAL, ISOLATION_UPPER REAL |

DATA roles: 0=m/z, 1=intensity, 2=retention time. Exactly one coordinate and one intensity array per object with identical length, including zero. Do not trust row count >= 2 or row iteration order. Explicitly associate blobs by ID, validate roles and reject ambiguous ownership (both/neither spectrum/chromatogram), duplicate roles, or invalid references.

C++ decoder accepts only compression codes 1 (zlib doubles), 5 (zlib Numpress linear), 6 (zlib Numpress SLOF). Its comments enumerate 0..7 but that is not implemented support. Reject other codes at this checkpoint and document it. Lossless source stores host-native doubles, including intensities promoted to double. The retained fixture is little-endian; use explicit byte decoding and document the architecture assumption, never pointer reinterpretation. Lossy m/z uses estimated linear fixed point with configured absolute accuracy; RT independently uses 0.05 seconds; intensity uses SLOF. Source disables its Numpress error-tolerance verification in both lossy paths. Rust may fail explicit rejected encodings safely but must not widen tolerances to hide discrepancies.

RUN_EXTRA is zlib-compressed mzML, not base64; original retained file has indexedmzML and ISO-8859-1 declaration. Its source construction copies descriptive metadata and clears peaks *and auxiliary float/string/integer arrays*. Consequently full metadata does not guarantee full original experiment recovery. Minimal SQL also stores only the first spectrum precursor/product/activation, loses drift-time units and unknown-polarity distinction (writer maps unknown to negative), and selected reads never consult RUN_EXTRA. These boundaries require accurate docs and fixtures. Decide explicitly whether auxiliary-array input is rejected or accepted with a documented source-loss policy; do not claim lossless experiment round trips for it.

Some source numeric metadata is streamed with precision 11; native typed REAL parameters preserve f64 and should be documented as a correction rather than reproducing decimal rounding. Current source binds RUN filename/native ID as blobs via the BLOB-only helper, though schema columns say TEXT; retained fixture uses TEXT. Accept valid UTF-8 TEXT or BLOB for those historical RUN columns without accepting arbitrary invalid bytes.

Indexes: DATA by each object ID; SPECTRUM by RT, MS level, run; RUN_EXTRA by run; CHROMATOGRAM by run; PRODUCT and PRECURSOR by spectrum/chromatogram ID. Source's last four indexes accidentally target DATA; do not copy that mistake.

## Useful implementation sequence and parallel boundaries

1. Freeze S0 connector first. Root integrates private connection access and feature/module changes. One owner implements handler schema/config/typed query helpers and table-derived lossless reading; a second independently ports SWATH queries/value output and prepares fixtures/corruption tests; a third audits source/provenance/docs. Share only agreed schema/query helpers, not concurrent edits to central files.
2. Add bound transactional lossless writes and full metadata snapshot read/write; validate table/snapshot IDs and counts. Round-trip zero records, zero-length arrays, quoted IDs, and repeated create/write using the same instance. Avoid whole-input pre-encoded copies: use bounded chunking within a single transaction.
3. Add configured lossy writing with existing codecs, retained binary fixture reading, and full class-test numeric comparisons. Then independent focused Claude/Fable 5.1 review with bounded tool-free patch/context prompts and a written final review. Iterate actionable findings.
4. Run remote kim tests in unique slots using the inspected gate helper; no local full builds. S1 is complete only once every mapped public operation, malformed-input boundary, feature combination, docs/provenance checks, and review have passed. Commit/push is parent-owned.

Later consumer S2 `SqMassFile` and streaming consumer work must use actual database IDs (or an explicit contiguous-ID invariant). Its current C++ transform produces an extra empty selected-read batch for zero/exact-multiple-of-500 counts; a candidate entry is included because source review discovered it now.

## Retained fixtures and exact expected checks

| Source fixture | SHA-256 | Purpose |
|---|---|---|
| SqliteMassFile_1.sqMass (94,208 bytes) | de5bc7b1f695b3da53ade6e6113a58f6dc9e3609b1efc78ebffc0bb10a0da55f | Original schema/codecs, 2 spectra + 1 chromatogram, RUN ID 12345, full metadata |
| MzMLSqliteHandler_1.mzML (667,674 bytes) | c57ce081bc065a5e00193ee939ee1336f2b7dc5a2a1cbf29dee840432a4dc4a7 | Canonical decoded comparison, 19,914 and 19,800 spectrum peaks, 48 chromatogram points |
| SwathFile.mzML (292,637 bytes) | a92df088d6e4b3d6c750da6aae179baedf085a88bbe97a51020dc7159cd65675 | 19 cycles of MS1 plus 5 windows (114 spectra) |

Source class expectations: metadata-only has no peaks; full read restores 2 spectra/1 chromatogram. Spectrum 0 peak 100 is m/z 204.817, intensity 3857.86. Chromatogram point 20 is RT 0.200695, intensity 147414.578125. Copy exact class-test comparison policies: intensity absolute 1e-4 OR ratio <=1.001; m/z absolute 1e-5 OR ratio <=1.000001; RT absolute .05 OR ratio <=1.000001; dedicated lossy chromatogram comparison uses 1.0002 relative. Confirm these against executable oracle output before declaring tier 1 parity; never substitute a single blanket epsilon.

RT literals: (.4738,.1)->[1], (.296,.1)->[0], (.296,1.1)->[0,1]; restricting to {1}/{0} yields that record; (0,.1)->[]; (.3,-.1)->[1]; (0,-.1)->[0]. SWATH windows first center 412.5 with 400..425 bounds, second 425..450, last 500..525. MS1 IDs first0/last108 count19; first window IDs1/109 count19; second2/110 count19.

Regression matrix: all unsupported compression codes; bad/truncated zlib and Numpress; odd lossless byte count; duplicate/missing roles; unequal lengths in both row orders; empty arrays; corrupt/oversized blobs; arbitrary DATA insertion order; IDs with gaps and duplicate selections; RUN count 0/2; invalid metadata snapshot/count/identity; SQL NULLs/invalid enum values; NUL/apostrophe text and run-path injection fixture; write rollback after later statement failure and unchanged counters; recreate same object; invalid configuration; SWATH chromatogram precursor NULL rows before/between matching spectra; repeated center with different widths. Use SQL statement failure injection (e.g. trigger abort) to test transaction cleanup without filesystem races.

## Evidence and unresolved decisions

The [retained Python probe](../tests/data/sqlite_s1_source_review/probe.py) and [JSON output](../tests/data/sqlite_s1_source_review/result.json) performed read-only original-fixture inspection plus independent in-memory SQLite queries. They confirm the wrong-index fixture and demonstrate SQL row orders/NULL rows underlying the reviewed C++ control-flow defects. They do not execute C++ and do not upgrade evidence to tiers 1/2. Their hashes and the inspected unchanged sqMass fixture are retained in the review manifest. For executed C++ parity, create the oracle driver outside the repository under ../oracle, hash source/compiler/results and verify it links the requested pin; do not claim the existing prebuilt library represents this revision without verification.

Before declaring support, settle: SwathMap value mapping/module placement; feature name; native selected-read ordering; explicit policy for auxiliary arrays and additional precursor/product metadata; portable lossless-byte endianness scope; strictness for snapshots inconsistent with tables. Source unaligned reinterpret_cast, unchecked sqlite3_step paths and unordered RT LIMIT are additional review questions, not newly executed defects in this plan. Existing CPP181-187 cover helper-level SQLite issues; findings JSON deliberately avoids duplicating those.
