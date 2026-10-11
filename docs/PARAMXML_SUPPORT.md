# OpenMS parameter XML (INI)

`format::paramxml` reads and writes the native `param::Param` hierarchy using the OpenMS parameter XML format through version 1.8.0. The `paramxml` feature uses the existing optional quick-xml dependency and is enabled by default. It can also be enabled alone.

`read`/`read_with_limits` and `load` return owned parameters. `read_into`/`load_into` use the source's accumulation semantics: absent keys survive, repeated keys update their value and tags, an empty incoming description preserves an existing nonempty description, and restrictions absent from the file remain in place. Unlike the source's incremental externally visible mutations, failure leaves the caller's entire tree unchanged. An internal builder shares a work budget across all records and validates the resulting tree once, avoiding a full tree copy per XML entry.

`write`/`write_with_limits` serialize and validate before writing to a stream, then flush and report errors. `store` validates before opening/truncating a file; `-` writes to stdout as in C++. I/O failure after writing begins can leave a partial external file.

## Source behavior covered

- Scalars: int, float/double, string, bool, input-file, output-file, output-prefix. Integer lists, float/double lists, string lists, and input/output-file lists include the empty-list case.
- Ordered entries and nested NODE sections, descriptions with the source `#br#` newline convention, XML escaping and tabs, and file/tag metadata.
- Required and advanced flags are independent and only the literal attribute value `true` adds the corresponding source tag. The reader does not inject schema defaults.
- Integer and floating restrictions use `min:max`, open endpoints, and the legacy `min-max` fallback. Malformed range *shape* is ignored as in the source; malformed numeric endpoints are errors. Numeric overflow and nonzero values underflowing to zero are checked errors, including the scalar cases where the C++ XML helper logs and substitutes zero. Default minima are `-i32::MAX` and `-f64::MAX`, not `i32::MIN` or negative infinity.
- String restrictions retain order and empty fields. File supported formats override restrictions. The native reader also retains supported formats for output-prefix values and equivalent legacy file restrictions, correcting cases the C++ reader loses on its own round trip.
- A bool input becomes a string constrained to `true,false`. Output uses `type="bool"` only for the source flag case: scalar value `false`, valid strings exactly `["true", "false"]`, and no overriding file type.
- NaN, infinities, signed zero and full finite f64 precision are preserved by the default writer. NaN is written as `NaN`. Unlike the source's limited decimal stream precision, the default writer's numbers use a round-trip representation. `WriteOptions::source()` writes the source's numbers instead, byte for byte (see *Writer options*).
- UTF-8, declared ISO-8859-1 and ASCII, and UTF-16 little/big endian input are supported. Output declares and emits UTF-8 by default; see *Writer options* for the source's ISO-8859-1 declaration. XML 1.0 literal line endings/attribute whitespace are normalized, while character references preserve tabs/newlines/carriage returns.

The source's nonsemantic schema attributes `short_description`, `position`, and scalar `default` are accepted and ignored, as in ParamXMLHandler. No remote schema, external entity or DTD is fetched. The inherited `XMLFile::isValid` is `paramxml::is_valid`, with the optional `xml-schema` feature: real XSD validation against the bundled, unchanged `Param_1_8_0.xsd`, with the four stored parameter sets of `ParamXMLFile_test.cpp:63-243` ported ([XML schema validation](XML_SCHEMA_SUPPORT.md)).

## Writer options

`write_with_options(output, param, limits, options)` and
`store_with_options(path, param, options)` take a `WriteOptions`, whose
`encoding` is an `OutputEncoding` and whose `numbers` is a `NumberText`.
`write`, `write_with_limits` and `store` use the default and are unchanged.

| `OutputEncoding` | Declaration | Bytes |
| --- | --- | --- |
| `Utf8` (default) | `encoding="UTF-8"` | every character as UTF-8 |
| `Latin1` (`WriteOptions::source()`) | `encoding="ISO-8859-1"`, as the source `ParamXMLFile::store` | a character up to U+00FF as its ISO-8859-1 byte, any other as a hexadecimal character reference (`&#x20AC;`) |

| `NumberText` | A floating-point item or restriction | A floating-point list item | An integer item |
| --- | --- | --- | --- |
| `RoundTrip` (default) | the shortest text that reads back to the same `f64` (`4`, `100000000000000000000`), `INF`, `-INF`, `NaN` | the same | refused beyond 32 bits, which the source's reader cannot read back |
| `Source` (`WriteOptions::source()`) | `ParamValue::toString()` / `StringUtils::toStr(double)` (`ParamXMLFile.cpp:110, 251, 256`): 15 fractional digits in `[1e-2, 1e4)` trimmed to one (`4.0`, `99.632000000000005`), the shortest round trip in scientific notation outside it (`1.0e20`, `1.5e-05`), `inf`, `-inf`, `NaN` | `os << value` at precision 15 (`ParamXMLFile.cpp:64, 328`), `%.15g`: `4`, `1e+20`, `0.1`, glibc's `nan`, `-nan`, `inf` | the 64-bit value, `std::to_string` (`ParamXMLFile.cpp:106`), as the source writes one a JSON `-ini` set beyond 32 bits although its reader then refuses it |

With `WriteOptions::source()` the source writer golden
`ParamXMLFile_test_writeXMLToStream.xml` is reproduced byte for byte, declaration
included, and so is every retained C++ `-write_ini` file of a ported tool
(`tests/topp_write_ini_text.rs`). The TOPP `-write_ini` option uses it
(`docs/TOPP_CLI_SUPPORT.md`). `Source` keeps 15 significant digits, as the
source does, so a value with more is not read back unchanged.

Native difference: the source declares ISO-8859-1 but copies the UTF-8 bytes of
its strings, so a non-ASCII character it writes reads back as several
different characters, in this reader as in the source's. `Latin1` keeps the
bytes consistent with the declaration, so the same tree reads back unchanged;
for ASCII content the bytes are the source's. Everything else (validation,
limits, escaping, atomic serialization before the destination is touched) is
shared with the default writer.

[Tests](../tests/paramxml.rs) (`source_declaration_writer_option`) compare the
golden byte for byte, check the declaration of both encodings, the ISO-8859-1
bytes and character references of values, descriptions and list items, and the
round trip through `read`, `store_with_options` and `load`.

## Checked native adaptations

Unknown types/elements/attributes, newer format versions, unsupported encodings, malformed XML and entity declarations return errors rather than silently dropping records. The source schema advertises additional types such as int-pair and input-prefix that its own ParamXMLHandler cannot load; those remain explicit unsupported inputs here.

The writer rejects values that the source representation would silently lose: Empty parameters, singleton empty string restrictions, comma-containing restrictions/tags, empty tags, literal `#br#` in descriptions, invalid XML characters, and, with the default `NumberText`, scalar integers outside the source reader's i32 range, nondefault restrictions on an inactive value type, and empty/colon-containing/duplicate sibling draft names. Native ParamValue itself still supports i64. Empty sections are retained by the native writer, including trees with no entries; this avoids source iterator-trace corner cases. These checks do not require the current value to satisfy its configured restriction: source parameter files can store defaults or editable values that have not yet passed parameter validation.

Limits default to 64 MiB raw/decoded/output XML, one million XML elements, one million items per list or comma-separated attribute, 128 XML nesting levels, and 65,536 bytes per fully qualified parameter path. The parameter layer's shared work/allocation ceilings also apply; these can reject a pathological tree before those XML ceilings. Attribute field slots and bytes are checked before splitting. Output escapes strings, tags and description newline markers incrementally within the byte budget. Reads and writes fail before publishing a partially parsed tree or validation-failing output.

## Evidence

[Tests](../tests/paramxml.rs) load four immutable upstream fixtures, compare the entire source writer golden with only its encoding declaration changed to UTF-8, exercise all source list/scalar alternatives, legacy tags, optional-required regressions, restrictions, special floats, Unicode/encodings, malformed XML, bounded failures and stream errors. A native workflow loads a threshold from INI, filters a DTA2D experiment, writes it, and reloads the expected retained peaks. It is a library integration test, not a full TOPP executable or C++ differential certification.

[Provenance](../tests/data/paramxml_provenance.json) pins the original headers, handlers, class tests, schema and fixtures to `6bfc0e4711105f4eda2fea86812a83af7c7e791f`. No generated Rust result is used as an expected scientific oracle.

## Remaining input transport

The native path loader recognizes plain, gzip and bzip2 input by content.
`load_with_limits` and `load_into_with_limits` bound decoded bytes; failed reads
preserve the existing parameter tree. A one-entry ZIP archive is read through
its file entry (`src/format/zip_ifstream.rs`), as `XMLFile::parse_` reads it
(`ParamXMLFile.cpp:357-361` calls `parse_`). `tests/paramxml_zip_input.rs` loads
stored, deflated and bzip2 entries, an entry below a directory entry, two entries
and an empty archive, and compares each with the Linux x86_64 Release build's
`ParamXMLFile::load` followed by `writeXMLToStream`, byte for byte
(`../oracle/integ-B-paramxml-zip`, ibminode06, two identical runs;
`tests/data/paramxml_zip_input_provenance.json`). Source
ParamXML output is plain, including its stdout option, regardless of suffix.
