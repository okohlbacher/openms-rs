# TextFile and CsvFile

`format::TextFile` and `format::CsvFile` implement the complete operational surface of the retained C++ helpers using owned UTF-8 strings and standard Rust readers, writers and filesystem paths. They do not require an optional feature or a new dependency. This is the source's literal CSV helper, not an RFC 4180 parser.

## TextFile

`new`/`default` creates an empty buffer. `from_reader` and `from_path` construct a buffer; `load_reader` and `load` replace it. These functions take `text::ReadOptions`, whose defaults match the source: no trimming, `first_n = -1`, retain empty lines, and no comment prefix. A positive `first_n` stops after that many retained lines without reading later lines. Zero and every negative value read until EOF; zero does **not** request an empty buffer.

Reading recognizes LF, CR and CRLF independently, including mixed endings and a final unterminated line. Terminators are removed. A terminal newline does not create an extra empty line. Optional trimming removes exactly ASCII space, tab, CR and LF. Empty-line skipping and comment-prefix filtering happen after trimming; skipped lines do not count toward `first_n`. An empty comment prefix disables comment filtering. UTF-8 BOMs and other characters are not automatically removed.

`get_line` exposes the single-line operation, returning a boolean for line/EOF and writing into a supplied String. `get_line_with_limits` applies limits to one call; the materializing reader instead shares counters across the whole requested input. At EOF the supplied String is cleared. Native errors preserve the previous supplied String rather than exposing the source's partially consumed line.

`add_line` is the native replacement for source `addLine` and `operator<<`. `lines`, `iter`, `iter_mut`, and borrowed `IntoIterator` provide immutable/mutable forward and reverse traversal. Mutable access can change line lengths; append and write operations revalidate the affected buffer's storage afterward. `clear`, `len`, and `is_empty` provide ordinary native buffer operations.

`write` and `store` preserve source suffix handling: a terminal CRLF is replaced by the platform newline, a terminal LF is retained through platform text translation, and entries without a terminal LF get one newline appended. A bare terminal CR is not stripped, and embedded line endings are not normalized individually. Windows text output translates each LF to CRLF as the source text-mode stream does; other targets emit LF. Buffer entries can contain embedded newlines, so storing and loading them may change the number of entries.

## CsvFile

CSV options select a single-byte separator (default comma), `item_enclosed` (default false), `first_n` (default -1), and text limits. `new`/`default` and `with_options` create an empty CSV buffer. File/reader construction and loading deliberately differ:

- `from_path` and `from_reader` mirror the source filename constructor: lines are untrimmed.
- `load` and `load_reader` mirror source `load`: lines are trimmed before testing comments.

Both always discard lines beginning with `#` at the time of that test. This prefix cannot be disabled. Thus the constructor retains `  #x`, while `load` discards it. Empty lines remain rows, and comments do not count toward `first_n`.

`row_count`, `clear`, `add_row`, `get_row`, `row`, `write`, and `store` cover the source row/storage operations. Clearing retains separator and enclosure settings. `get_row(index, output)` returns the source split-success boolean separately from checked errors; `row(index)` returns the same boolean with an owned vector:

- An empty row returns `false` and an empty vector.
- A nonempty row without the delimiter returns `false` and one unchanged field, including any enclosing characters.
- Otherwise every literal delimiter splits, including delimiters inside quotes, and empty fields are retained.
- If enclosure is enabled after a successful split, the first and last **byte** of each field are removed, without checking that they are quotes. A one-byte field becomes empty. An empty field causes a checked error.

`add_row` performs a literal join. Enclosure wraps each item in double quotes without escaping quotes, separators or newlines already present in the item. This can produce data that splits differently when read; the port retains that source behavior. An empty item list adds an empty row. `store` writes the raw buffer through TextFile, rather than parsing or reconstructing its fields.

## Checked native boundaries

Source strings can contain arbitrary bytes; these native APIs use UTF-8. Invalid UTF-8 input and byte splits/enclosure removal that break a Unicode character are errors. A non-ASCII separator remains selectable for source-compatible byte searches, but joining multiple fields with it is rejected because that single inserted byte cannot form valid UTF-8. ASCII separators, including NUL or newline, retain their literal behavior.

Out-of-range row indices, malformed empty enclosed fields, allocation failures, and I/O failures return `Result` errors. Loading replaces both contents and CSV options only after success. `get_row` likewise preserves its previous destination on error. The source can clear or partially update destinations before failing; that partial state is not reproduced.

`text::Limits` defaults to 256 MiB consumed input, 16 MiB per line/entry excluding input terminators, one million consumed or retained lines, 256 MiB conservative stored payload/slot accounting, and 512 MiB output. Limits can be lowered from those hard ceilings. Input limits include comments and empty lines even if skipped. The temporary line buffer is separately bounded by the line limit. Stored String capacity and vector slots are accounted conservatively; trimmed input is shrunk before retention. CSV row parsing/appending additionally caps fields at one million and checks field storage before allocation. Appending accounts for existing retained entries. Direct caller mutations through mutable text iterators are checked before subsequent append or output.

Output validates the complete buffer and output size before writing bytes; filesystem storage validates before creating/truncating the destination. Writers flush and report failure. An underlying I/O error can still leave a partial external file after writing begins. Distinct C++ exception classes and iostream flags map to Rust errors and the line/EOF boolean rather than an emulated stream state machine.

## Reference coverage

[Text tests](../tests/text_file.rs) cover all upstream text fixture positions and counts, source store literals, mixed line endings across one-byte reader buffers, comments/trimming/first_n interactions, early stopping, mutable iteration, resource limits and atomic failures. [CSV tests](../tests/csv_file.rs) cover all source `hello/world`, `the/dude`, `spectral/search`, `first/second/third`, and `4/5/6` rows, constructor/load differences, source non-RFC enclosure branches, Unicode boundaries, joins, field counts and storage failures.

The four source fixture files are copied byte-for-byte. [Provenance](../tests/data/text_csv_provenance.json) records the headers, implementations, class tests, StringUtils dependency and original fixture hashes at OpenMS4-core `6bfc0e4711105f4eda2fea86812a83af7c7e791f`. Source CSV test constructor/load exception sections enclosed in `#if 0` are not claimed as active upstream tests. Additional boundary expectations are independently derived from source branches; no Rust output was used to manufacture goldens and no C++ code was compiled or executed.
