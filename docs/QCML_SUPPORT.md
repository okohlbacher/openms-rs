# qcML support

Port of `src/openms/include/OpenMS/FORMAT/QcMLFile.h` (206 lines) with
`src/openms/source/FORMAT/QcMLFile.cpp` (2137 lines), at core SDK revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

- Rust: `src/format/qcml.rs` (one module, gated on the `paramxml` feature)
- Tests: `tests/qcml.rs` (60 cases; one per upstream `START_SECTION`, plus
  native hardening)
- Fixtures: `tests/data/QcMLFile_reload_A.qcML`,
  `tests/data/QcMLFile_reload_B.qcML` (both unmodified upstream),
  `tests/data/QcMLFile_store_shape.qcML`,
  `tests/data/QcMLFile_unicode.qcML`, `tests/data/QcMLFile_latin1.qcML`,
  `tests/data/QcMLFile_stylesheet.qcML`,
  `tests/data/QcMLFile_report_sheet.xsl`
- Provenance: `tests/data/qcml_provenance.json`
- Status: **partial** — every member of the header is ported except
  `collectQCData`, which is a 926-line QC-metric computation over the kernel
  and identification types rather than a format concern. See **Deferred**.

## What qcML is

A qcML document reports quality control for a mass-spectrometry run. Its root
`<qcML>` holds `<runQuality>` elements, one per acquisition, and
`<setQuality>` elements, one per group of runs, each identified by an `ID`
attribute. Inside either:

- `<qualityParameter>` is one CV-identified scalar: a `name`, an `ID`, a
  `cvRef`/`accession` pair, an optional `value`, an optional unit and an
  optional `flag`. It is always an empty element.
- `<attachment>` carries bulk evidence for a parameter it names through
  `qualityParameterRef`, as either an opaque `<binary>` payload or an inline
  `<table>` with a space-delimited `<tableColumnTypes>` header and one
  space-delimited `<tableRowValues>` per row.

Two accessions are structural rather than informational. `MS:1000577` inside a
`<runQuality>` names the run; inside a `<setQuality>` it names one of the set's
members instead and is not kept as a parameter. `QC:0000058` inside a
`<setQuality>` names the set. A run or set with neither is named after its own
identifier.

The document closes with a fixed three-entry `<cvList>`, and — when the writer
finds `share/OpenMS/XSL/QcML_report_sheet.xsl` — an embedded XSL stylesheet plus
the processing instruction and DOCTYPE that make a browser render the report as
HTML.

## API mapping — `QcMLFile.h`

Every public member of the header, in declaration order. Protected members and
base classes follow.

### `class QualityParameter` → [`QualityParameter`]

| C++ member | Rust | Notes |
|---|---|---|
| `std::string name` | `QualityParameter::name` | written as `name`, required |
| `std::string id` | `QualityParameter::id` | written as `ID`, required |
| `std::string value` | `QualityParameter::value` | omitted when empty |
| `std::string cvRef` | `QualityParameter::cv_ref` | written as `cvRef`, required |
| `std::string cvAcc` | `QualityParameter::cv_acc` | written as `accession`, required |
| `std::string unitRef` | `QualityParameter::unit_ref` | see **unit spellings** |
| `std::string unitAcc` | `QualityParameter::unit_acc` | see **unit spellings** |
| `std::string flag` | `QualityParameter::flag` | documented upstream as "cv accession of the unit", a copy-paste slip; it is a boolean marker |
| `QualityParameter()` | `QualityParameter::default()` | all eight fields empty |
| `QualityParameter(const QualityParameter&)` | `Clone` | source is `= default` |
| `operator=` | `Clone` + assignment | |
| `operator==` | `QualityParameter::same_name` | the source compares `name` alone; `PartialEq` here is structural |
| `operator<` | `PartialOrd`/`Ord` | `name`-primary, then the remaining fields |
| `operator>` | `PartialOrd`/`Ord` | as above |
| `toXMLString(UInt)` | `QualityParameter::to_xml_string`, `::to_xml_string_with_options` | indentation is bounded by `MAX_INDENTATION` |

### `class Attachment` → [`Attachment`]

| C++ member | Rust | Notes |
|---|---|---|
| `std::string name` | `Attachment::name` | written as `name`, required |
| `std::string id` | `Attachment::id` | written as `ID`, required; documented upstream as "Name" |
| `std::string value` | `Attachment::value` | omitted when empty |
| `std::string cvRef` | `Attachment::cv_ref` | required |
| `std::string cvAcc` | `Attachment::cv_acc` | written as `accession`, required |
| `std::string unitRef` | `Attachment::unit_ref` | see **unit spellings** |
| `std::string unitAcc` | `Attachment::unit_acc` | see **unit spellings** |
| `std::string binary` | `Attachment::binary` | opaque; neither side decodes base64 |
| `std::string qualityRef` | `Attachment::quality_ref` | `qualityParameterRef`; optional here, required by the source's reader |
| `std::vector<std::string> colTypes` | `Attachment::col_types` | bounded by `MAX_COLUMNS` |
| `std::vector<std::vector<std::string>> tableRows` | `Attachment::table_rows` | bounded by `MAX_ROWS` and `MAX_TABLE_CELLS` |
| `Attachment()` | `Attachment::default()` | |
| `Attachment(const Attachment&)` | `Clone` | source is `= default` |
| `operator=` | `Clone` + assignment | |
| `operator==` | `Attachment::same_name` | as `QualityParameter` |
| `operator<` / `operator>` | `PartialOrd`/`Ord` | as `QualityParameter` |
| `toXMLString(UInt)` | `Attachment::to_xml_string`, `::to_xml_string_with_options` | |
| `toCSVString(const std::string&)` | `Attachment::to_csv_string` | rejects an empty separator |
| — | `Attachment::has_table` | native predicate for the source's `!colTypes.empty() && !tableRows.empty()` condition, which decides both serialisers' branch |

### `class QcMLFile` → [`QcMLFile`]

| C++ member | Rust | Notes |
|---|---|---|
| `QcMLFile()` | `QcMLFile::new`, `Default` | |
| `~QcMLFile()` | `Drop` | nothing to release |
| `map2csv(table, separator)` | `qcml::map_to_csv` | free function; the source member ignores `this` |
| `exportIDstats(filename)` | `QcMLFile::export_id_stats` | `Option<String>`; `None` for the source's `""` |
| `registerRun(id, name)` | `QcMLFile::register_run` | `Result`; rejects an empty id or name |
| `registerSet(id, name, names)` | `QcMLFile::register_set` | `Result`; `names` is a `&BTreeSet<String>` |
| `addRunQualityParameter(r, qp)` | `QcMLFile::add_run_quality_parameter` | `Result`; refuses an unregistered run instead of dropping |
| `addRunAttachment(r, at)` | `QcMLFile::add_run_attachment` | `Result`; creates the run bucket, as the source does |
| `addSetQualityParameter(r, qp)` | `QcMLFile::add_set_quality_parameter` | |
| `addSetAttachment(r, at)` | `QcMLFile::add_set_attachment` | |
| `removeAttachment(r, ids, at = "")` | `QcMLFile::remove_attachments_by_quality_ref` | `at` becomes `Option<&str>`; the default empty string is `None` |
| `removeAttachment(r, at)` | `QcMLFile::remove_attachments_by_accession` | distinct name for the overload |
| `removeAllAttachments(at)` | `QcMLFile::remove_all_attachments` | run-map-only reach preserved |
| `removeQualityParameter(r, ids)` | `QcMLFile::remove_quality_parameters` | `ids` is `&[String]`, not a mutable vector |
| `merge(addendum, setname = "")` | `QcMLFile::merge` | `setname` becomes `Option<&str>`; dedup policy is explicit |
| `collectSetParameter(setname, qp, ret)` | `QcMLFile::collect_set_parameter` | returns the values; `&self`, so no phantom entries |
| `exportAttachment(filename, qpname)` | `QcMLFile::export_attachment` | `Result<Option<String>>` |
| `exportQP(filename, qpname)` | `QcMLFile::export_quality_parameter` | `Option<&str>`; the source's `"N/A"` is `qcml::NOT_FOUND` |
| `exportQPs(filename, qpnames)` | `QcMLFile::export_quality_parameters` | byte-identical to the source, trailing comma included |
| `getRunIDs(ids)` | `QcMLFile::run_ids` | borrowed iterator replaces the out-parameter |
| `getRunNames(ids)` | `QcMLFile::run_names` | |
| `existsRun(filename, checkname = false)` | `QcMLFile::exists_run` / `::exists_run_or_name` | the boolean parameter becomes two names |
| `existsSet(filename, checkname = false)` | `QcMLFile::exists_set` / `::exists_set_or_name` | |
| `existsRunQualityParameter(filename, qpname, ids)` | `QcMLFile::exists_run_quality_parameter` | returns the ids; matches `cvAcc`, as the source does |
| `existsSetQualityParameter(filename, qpname, ids)` | `QcMLFile::exists_set_quality_parameter` | |
| `collectQCData(prot_ids, pep_ids, feature_map, consensus_map, inputfile_raw, remove_duplicate_features, exp)` | **not ported** | see **Deferred**; its `@param` constraints are recorded there |
| `store(filename)` | `QcMLFile::store`, `::store_with_options`, `::to_xml_string`, `::to_xml_string_with_options` | published atomically |
| `load(filename)` | `qcml::load`, `::load_with_limits`, `::read`, `::read_with_limits` | returns a new document rather than clearing in place |

### Protected members

| C++ member | Rust | Notes |
|---|---|---|
| `onStartElement(qname, attributes)` | not ported as an item | folded into the module-private reader; the source's SAX dispatch is not an API |
| `onEndElement(qname)` | not ported as an item | as above |
| `onCharacters(chars, length)` | not ported as an item | as above; the port assembles an element's complete text before splitting it |
| `runQualityQPs_` | `QcMLFile::run_quality_parameters` (read), `add_run_quality_parameter` (write) | a Rust type has no protected access, so the collections are reachable through borrowing accessors |
| `runQualityAts_` | `QcMLFile::run_attachments`, `add_run_attachment` | |
| `setQualityQPs_` | `QcMLFile::set_quality_parameters`, `add_set_quality_parameter` | |
| `setQualityAts_` | `QcMLFile::set_attachments`, `add_set_attachment` | |
| `setQualityQPs_members_` | `QcMLFile::set_members`, `register_set` | |
| `run_Name_ID_map_` | `QcMLFile::run_names`, `run_id_for_name` | |
| `set_Name_ID_map_` | `QcMLFile::set_names`, `set_id_for_name` | |
| `tag_`, `progress_`, `qp_`, `at_`, `row_`, `header_`, `name_`, `run_id_`, `names_`, `qps_`, `ats_` | not ported | transient parser state held on the object because the SAX callbacks cannot pass arguments. The port keeps it in a parser value that lives only for one read, which is what removes the state leaks recorded below |

### Base classes

| C++ base | Rust | Notes |
|---|---|---|
| `Internal::XMLHandler` | not ported | the port does not implement OpenMS's SAX handler interface |
| `Internal::XMLFile` | `qcml::VERSION` | constructed with an empty schema path and version `"0.7"`; neither side validates a schema |
| `ProgressLogger` | not ported | see **Deferred** |

### Native additions

`Limits`, `MergeOptions`, `WriteOptions`, `Stylesheet`, `QcMLFile::is_empty`,
`QcMLFile::set_ids`, `QcMLFile::set_names`, the `MAX_*` ceilings on all three
types, and the constants `qcml::VERSION`, `NOT_FOUND`, `RUN_NAME_ACCESSION`,
`SET_NAME_ACCESSION` and `SET_MEMBER_ACCESSION` — the last three name accessions
the source spells inline at five separate call sites.

## Unit spellings

The source's writer emits `unitRef=` and `unitAcc=` (`QcMLFile.cpp:88-95`,
`199-206`), while its reader parses `unitAccession` and `unitCvRef`
(`825-826`, `854-855`). A unit written by OpenMS therefore never survives being
read back by OpenMS.

This port's reader accepts **both** spellings, preferring the schema-shaped
`unitAccession`/`unitCvRef` when a document carries both. Its writer emits the
schema spellings by default, so a unit round-trips, and the source spellings
under `WriteOptions::source_unit_attributes`, so output can still match OpenMS
byte for byte.

## Preserved source conventions

These are reproduced deliberately, defects included, because they are what a
qcML consumer or the report stylesheet sees.

- **Element layout.** `<attachment ` is followed by ` name=`, so two spaces
  separate the tag from its first attribute. The `<table>` wrapper is written
  with no indentation and no newline, and `</table>` is followed directly by the
  indented `</attachment>`. Attribute order is `name`, `ID`, `cvRef`,
  `accession`, `value`, unit reference, unit accession, then `flag` or
  `qualityParameterRef`.
- **The fixed `cvList`.** Three hard-coded `<cv>` entries with fixed URIs and
  versions (PSI-MS 3.41.0, QC-CV 0.1.1, unit 1.0.0), independent of what the
  document actually references.
- **Lexical output order,** over the union of the parameter and attachment map
  keys. A run known only through `add_run_attachment` is written even though
  `exists_run` and `run_ids` do not see it, because both read the parameter map.
- **`exportQP`'s asymmetry.** A run parameter is matched on `accession`, a set
  parameter on `name` (`QcMLFile.cpp:638` and `659`). Undocumented upstream and
  almost certainly accidental, but it decides which lookups succeed.
- **`exportQPs`' trailing comma** and its `"N/A"` placeholder, so the field
  count stays stable.
- **`map2csv`'s misalignment.** Columns come from the first row's keys only, and
  a row missing one emits neither the cell nor its separator. `export_id_stats`
  builds exactly such a table, because the `QC:0000043`-`47` and
  `QC:0000053`-`57` CV names differ, so its `ms2` row is normally its label
  alone. Both are documented at the item with a `# Warning`.
- **`remove_all_attachments`' reach.** It iterates the run attachment map only,
  so a set is cleaned only when its identifier is also a run identifier with an
  attachment list, despite the header comment saying "all runs/sets".
- **`register_run`/`register_set` reset** the entry they name, so registering an
  existing identifier discards what it held.
- **A nameless entry is named after its own identifier**, and the name a
  registration supplies is stored in no element of its own.
- **`Attachment::to_csv_string`'s substitution**: the separator becomes `_`, or
  `$` when the separator is itself `_`, and each assembled line is trimmed of
  space, tab, carriage return and line feed.
- **Rows are written with exactly the cells they hold.** Neither side pads or
  checks a row against the column count.

## Native differences

Each is documented at the item that carries it. The rule throughout: where the
source silently discards information, the native default refuses and
`WriteOptions::source()`, `MergeOptions::source()` or
`WriteOptions::drop_unrepresentable` selects the source behaviour.

1. **Structural equality.** `operator==` on both records compares `name` alone.
   Rust's `==` is structural and `same_name` is the source predicate. `Ord` is
   `name`-primary with deterministic ties, which also makes `merge` deterministic
   where the source's unstable `std::sort` is not.
2. **`merge` deduplication.** The source's `std::unique` with a name-only
   comparator collapses distinct parameters that share a name.
   `MergeOptions::default()` collapses only exact duplicates;
   `MergeOptions::source()` reproduces the collapse.
3. **The `flag` value.** The source writes the literal `flag="true"` for any
   non-empty flag and reads the attribute back verbatim, so the value is lost.
   The default writer emits the value; `source_flag_literal` emits `"true"`.
4. **Unit spellings**, as above.
5. **`qualityParameterRef` is optional on read.** The source omits it when empty
   but reads it with `attributeAsString_`, which raises a fatal parse error on a
   missing attribute, so C++ cannot re-read its own attachment that hangs off the
   run. This port treats it as optional.
6. **Table cells.** The source substitutes `' '` with `'_'` in column types and
   *intends* to in row values, but concatenates the unsubstituted original
   (`QcMLFile.cpp:236-242`), so a cell containing a space silently becomes two
   cells on reload. The default writer refuses an empty cell or one containing
   XML whitespace, which a space-delimited table cannot represent;
   `source_table_text` reproduces the source, defect and all.
7. **Incomplete and over-complete attachments.** The source returns the empty
   string for an attachment with neither a binary payload nor a complete table,
   so `store` writes nothing for it, and writes only the binary when both are
   present, dropping the table. The default refuses both;
   `drop_unrepresentable` restores the source behaviour.
8. **Unpersistable names.** A run or set name that differs from its identifier
   survives a write only inside an `MS:1000577` or `QC:0000058` parameter. The
   default writer refuses to write a name that no parameter carries;
   `drop_unrepresentable` drops it, as the source does.
9. **Set membership recovery.** `store` documents each member as a `QC:0000005`
   parameter whose `ID` is the member's run identifier, but the source's reader
   only records members from `MS:1000577`, so C++ loses membership on every
   round trip. This port also reads `QC:0000005`, keeping the parameter as well.
   Its writer resolves a member through the name map and refuses one it cannot
   resolve; `drop_unrepresentable` reproduces the source's identifier-only
   lookup and its silent skip.
10. **XML escaping.** The source concatenates attribute values verbatim, so a
    value containing `&`, `<`, `>` or `"` produces a document no parser accepts.
    This port always escapes, and writes tab, line feed and carriage return as
    character references so attribute-value normalization cannot eat them. A
    character XML 1.0 cannot represent is `Error::InvalidValue`.
11. **Required attributes are checked on write.** An empty `name`, `ID`, `cvRef`
    or `accession` is `Error::MissingInformation`; the source writes `name=""`,
    which its own reader then accepts as an empty name.
12. **No phantom entries.** The source's removal methods and
    `collectSetParameter` reach their maps with `operator[]`, registering an
    empty run and set for an unknown key. The port touches only existing lists,
    and `collect_set_parameter` takes `&self`.
13. **Misplaced elements are refused.** The source decides a parameter's role by
    `parent_tag == "runQuality"` with an unconditional `else`, so a
    `qualityParameter` outside any entry is treated as a set parameter and
    attached to whichever entry closes next; a stray `tableRowValues` lands in
    the next attachment's table. The port makes all of these `Error::Parse`, and
    also rejects a nested or duplicated `runQuality`/`setQuality`.
14. **No parser state leaks.** The source never clears `names_`, so every set
    inherits the members of the sets before it, in the same file and across a
    second `load`; and it leaves `qp_` populated after a set-member parameter, so
    the next parameter's absent optional attributes keep the previous values. The
    port's parser state lives for one element or one entry.
15. **Complete character data.** The source uses the first non-empty character
    notification for a row and `StringUtils::split` clears its output, so a
    second chunk overwrites the first and a row containing an entity reference
    loses everything before it. The port assembles an element's whole text first.
    Repeated `<binary>` elements still concatenate, as the source does.
16. **Atomic writes and reads.** `store` serialises fully and publishes through
    the crate's atomic writer; the source streams into an `ofstream` and a
    failure part-way leaves a truncated report. `load` is a function returning a
    new document: the source's member clears its maps first, so a throwing parse
    leaves the object empty.
17. **`Option` instead of sentinels.** `export_attachment` returns `None` for
    "no match" and `Some("")` for "matched, but no table", which the source's
    single `""` cannot distinguish; `export_quality_parameter` returns `None`
    rather than `"N/A"`, and `export_quality_parameters` still emits the
    sentinel.
18. **No progress logging.** The source drives its `ProgressLogger` base during
    `load`.

## Checked boundaries and evidence

Resource ceilings, all refusals checked before the matching allocation. The
source has none: `parse_()` hands the whole file to the XML parser and the
handler appends to `std::vector` and `std::map` until the allocator fails.

Three of these once held only on the way out of the reader, which is what the
first release of this port claimed they did not do. `Attachment::MAX_COLUMNS`
was read off the column vector `split_cells` had already built;
`Attachment::MAX_TABLE_CELLS` and `Attachment::MAX_TEXT_BYTES` were charged
only when the finished attachment was committed into the document, by
`Attachment::preflight_table` and `check_text`, so a document could assemble an
oversized table or accumulate repeated `<binary>` elements in memory first and
be refused afterwards. All three are now decided on the captured text of the
element being read - already bounded by `Limits::max_text_bytes` - before any
cell `String`, row `Vec` or payload append exists, with the running cell count
carried on the parser and reset for each `<attachment>`. The commit-time and
write-time checks remain as the in-memory API's own guard.

| Ceiling | Default | Guards |
|---|---|---|
| `Limits::max_input_bytes` | 256 MiB | decoded document size |
| `Limits::max_depth` | 100 | nesting, which an embedded XSL stylesheet uses up |
| `Limits::max_elements` | 4,000,000 | markup nodes: elements, and the comments and processing instructions the port discards, each of which costs the parser the same walk as an element |
| `Limits::max_text_bytes` | 64 MiB | one `<binary>` or table element's character data, and one attribute value |
| `Limits::max_doctype_bytes` | 4 KiB | the internal DOCTYPE subset the source's own writer emits |
| `MAX_ATTRIBUTES` (module-private) | 64 | attributes on one element; the duplicate-name scan is quadratic in this count and no qcML element declares more than nine |
| `QcMLFile::MAX_ENTRIES` | 1,000,000 | runs, and separately sets |
| `QcMLFile::MAX_PARAMETERS_PER_ENTRY` | 1,000,000 | parameters per run or set |
| `QcMLFile::MAX_ATTACHMENTS_PER_ENTRY` | 1,000,000 | attachments per run or set |
| `QcMLFile::MAX_SET_MEMBERS` | 1,000,000 | member names per set |
| `QcMLFile::MAX_OUTPUT_BYTES` | 1 GiB | any serialisation or export, preflighted from field sizes |
| `Attachment::MAX_COLUMNS` | 100,000 | column types per table, counted in the captured text before the cells are built |
| `Attachment::MAX_ROWS` | 4,000,000 | rows per table |
| `Attachment::MAX_TABLE_CELLS` | 16,000,000 | cells per table, accumulated row by row as they are read |
| `Attachment::MAX_TEXT_BYTES` | 64 MiB | one cell, scalar field or binary payload, charged per `<binary>` element because repeated ones accumulate |
| `QualityParameter::MAX_TEXT_BYTES` | 16 MiB | one field |
| `QualityParameter::MAX_INDENTATION` | 64 | the source builds `std::string indent(level, '\t')` from an unchecked `UInt` |
| `Stylesheet::MAX_BYTES` | 16 MiB | an injected report stylesheet body |

Input policy, all independently derived:

- **DTDs.** A DOCTYPE is accepted only when it declares no `ENTITY` and names no
  external subset, which is exactly the shape `store` emits alongside a report
  stylesheet, so the port can read its own and OpenMS's stylesheet output
  without ever expanding an entity. Anything else is `Error::Unsupported`.
- **Encodings.** UTF-8 with or without a byte-order mark; UTF-16 of either
  endianness, checked against the declaration; US-ASCII; and ISO-8859-1, which
  the source's writer declares while emitting `std::string` bytes unchanged. A
  document declaring ISO-8859-1 whose bytes are valid UTF-8 is read as UTF-8,
  because that is what the source produces from UTF-8 strings; only bytes that
  are not valid UTF-8 are read as Latin-1. Any other declared encoding, and a
  declaration whose version is not 1.0, is `Error::Unsupported`.
- **No string is byte-sliced.** Splitting, trimming and prefixing go through
  `str` operations, and every index into decoded bytes goes through `get`. The
  non-ASCII fixture exists because an audit in this project found a reachable
  panic from byte-slicing a path with a non-ASCII component.

Parse cost is linear in the document, not quadratic. The reader reports a line
number with every diagnostic, and the first release of this port derived it by
counting the newlines of the whole byte prefix on **every** event, so a
document of many small nodes cost O(nodes x document bytes). Measured in
release on the gate node, a `<runQuality>` of self-closing
`<qualityParameter/>` elements took

| Elements | Document | Before | After |
|---|---|---|---|
| 10,000 | 777,868 B | 2,677.178 ms | 10.373 ms |
| 20,000 | 1,577,868 B | 10,781.590 ms | 18.468 ms |
| 40,000 | 3,177,868 B | 43,533.878 ms | 39.464 ms |

— twice the input for four times the time before, and for twice the time after.
`Parser::line_at` now carries the line number forward from the previous event's
byte position, so every byte of the document is examined once over the whole
parse. Comments drove exactly the same scan while being charged to no ceiling
at all, which is why they now count against `Limits::max_elements` alongside
processing instructions. `tests/qcml.rs` pins both: it times the same three
sizes, with and without a comment before every element, and refuses a
fourfold-larger document that costs more than eight times as much.

Evidence, in full in `tests/data/qcml_provenance.json`: **tier 3 source review**,
no tier 1 differential. No C++ was built or executed and no C++ output was
retained — `QcMLFile` has no TOPP test with a retained qcML output upstream, and
the class test's two fixtures are inputs, not captured outputs. All 43
`START_SECTION`s are ported, one Rust test each, named `section_*`:

- The **12 sections that assert at least one value** are transcribed literally
  and marked `upstream:` at the assertion. Exactly one of them, `load` with 11
  assertion macros, is above the five-macro threshold that forbids mapping, and
  it is ported in full. Among them: `existsRun("abc")` and
  `existsRun("somerun", true)` after `registerRun("abc","somerun")`;
  `existsSet("def")` and `existsSet("someset", true)` after
  `registerSet("def","someset",{"somerun1","somerun2"})`; `getRunIDs` yielding
  `["123","456"]`; `getRunNames` yielding `["testrun1","testrun2"]`; the
  `"somename"`/`"id"`/`"somevalue"` parameter surviving a copy and an
  assignment; name-only equality; `"somename" < "tomename"`;
  `"somename" > "romename"`; and all ten assertions of the `load` section on
  `QcMLFile_reload_A/B` — one run name `runAlpha`, then `runBeta` with A's names
  no longer resolving.
- The **29 `NOT_TESTABLE` sections**, plus the **2 destructor sections** whose
  bodies are only `delete ptr` and the comment "uh, twice?! No!", carry no
  upstream expected value at all. Each of those 31 is ported with an expectation
  derived from `QcMLFile.cpp` and marked `derived:` with the line it comes from. The byte layout asserted for
  `toXMLString` and `store` is transcribed from `QcMLFile.cpp:78-255` and
  `1972-2135`.
- `sections_mapped_with_evidence` is therefore **0**: no section is claimed
  against another test's coverage.

`tests/data/QcMLFile_store_shape.qcML` is transcribed writer output, which makes
it tier 3 and not tier 1; it exists because no upstream fixture contains an
attachment, so nothing upstream exercises tables, binary payloads, units or set
membership. The resource ceilings, the DOCTYPE and encoding policies, the
non-ASCII and Latin-1 fixtures, the malformed-input corpus and every
error-variant choice are tier 4, independently derived.

## Deferred

- **`collectQCData`** (`QcMLFile.cpp:1042-1967`, 926 lines) and its file-static
  helper `calculateSNmedian` are not ported, which is why the header's status is
  `partial`. Three reasons, and all three have to be dealt with together:
  1. It is a QC-metric computation, not a format concern: it reads a
     `FeatureMap`, a `ConsensusMap`, an `MSExperiment` with chromatograms and
     acquisition info, a `std::vector<ProteinIdentification>` with search
     parameters and a `PeptideIdentificationList`, and writes ~30 parameters and
     ~8 attachments into the document.
  2. It resolves every parameter's display name through
     `ControlledVocabulary::getTerm` over `CV/psi-ms.obo`, `CV/qc-cv.obo` and
     `CV/qc-cv-legacy.obo`. This crate bundles only `psi-ms.obo`, so the other
     two are new resources. Every lookup has a `catch (...)` fallback literal,
     so a port without them would diverge silently on every term that exists.
  3. It carries at least seven defects, several of them memory-unsafe, so a port
     has to decide case by case what to reproduce. All are recorded with line
     anchors in the provenance manifest: the integer-division slump percentages
     that are identically zero above 100 spectra; the duplicate-feature loop
     that indexes one past the end and then discards its result; the
     unchecked `getPrecursors().front()` and `exp.begin()`; m/z accumulated in
     `UInt`; the RIC-drop test that fires on nearly every spectrum; the
     missed-cleavage loop that is undefined for an empty sequence; and the
     consensus attachment that declares eight columns and writes seven values
     under a placeholder accession.

  Its header `@param` constraints are recorded here so they are not lost:
  `prot_ids` are the protein identifications from the ID file, `pep_ids` the
  peptide identifications, `feature_map` comes from a featureXML,
  `consensus_map` from a consensusXML, `inputfile_raw` is the mzML input file
  name (only its stem is used, as the run identifier *and* name),
  `remove_duplicate_features` removes duplicates in a set of merged features,
  and `exp` **requires `sortSpectra()` and `updateRanges()` to have been called
  first**.
- **The report stylesheet** is not bundled. `store` injects
  `share/OpenMS/XSL/QcML_report_sheet.xsl` (23,513 bytes) so a browser renders
  the qcML as HTML, and warns "No qcml stylesheet found, result will not be
  viewable in a browser!" when it cannot find it. This port reproduces the
  stylesheet-absent path by default and accepts the sheet through
  `WriteOptions::stylesheet` with `Stylesheet::from_file_text`, which applies the
  source's first-line strip, so a caller holding the file gets the same output.
- **`ProgressLogger`** is a base class the source drives during `load`. Wiring
  the crate's `src/concept/progress_logger.rs` into a parse would add a logging
  side effect to a library read, which no other format module here does.
- **The feature gate.** The module needs the crate's `quick-xml` dependency and
  is gated on `paramxml`, the smallest existing feature that provides it, because
  adding a feature is a `Cargo.toml` change outside a port package's scope. A
  dedicated `qcml` feature would be the honest gate.
- **Base64.** An attachment's `<binary>` is carried opaquely, as the source does.
  The crate's `base64` dependency is optional behind `numpress`, so a decode
  helper here would either widen the feature graph or be conditionally absent.
- **No OpenMP gap.** `QcMLFile.cpp` carries no `#pragma omp`; the source is
  serial and so is the port.

[`QualityParameter`]: ../src/format/qcml.rs
[`Attachment`]: ../src/format/qcml.rs
[`QcMLFile`]: ../src/format/qcml.rs
