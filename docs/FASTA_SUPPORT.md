# FASTA files and streaming

`format::fasta` implements the complete public FASTAFile operation set from
OpenMS4-core `54a232fe2cae9c590d5c997fa49d20e7769860fb`. Its header, implementation,
class test and original data file are byte-identical to the previous `6bfc0e4`
reference. Rust ownership replaces C++ stream references and move operations;
checked boundaries below remain deliberate native differences.

| Source operation | Native API |
| --- | --- |
| Entry construction/copy/equality | `FASTAEntry::new`, `Default`, `Clone`, `PartialEq`, public fields |
| `headerMatches`, `sequenceMatches` | `header_matches`, `sequence_matches` |
| `readStart`, `readNext` | `FASTAFile::read_start`, `read_next` |
| `position`, `setPosition`, `atEnd` | `position`, `set_position`, `at_end` |
| `readStartWithProgress`, `readNextWithProgress` | `read_start_with_progress`, `read_next_with_progress` |
| `writeStart`, `writeNext`, `writeEnd` | `write_start`, `write_next`, `write_end` |
| `load`, `store` | `load` returning owned entries, `load_into`, `store` |
| Inherited progress controls | Public `FASTAFile::progress: ProgressLogger` |

`FastaReader<R: BufRead>` provides the same parser over caller-owned streams:
`next_entry`, atomic `read_next`, iteration, `at_end`, `entries_read` and
`into_inner`. Positioning additionally requires `Seek`. `FastaWriter<W: Write>`
provides `write_entry` and consuming `finish`, which returns the writer after a
checked flush. Existing free `read`/`write` remain available, with explicit
`read_with_options`/`write_with_options` variants. All three owners accept
`FastaOptions`; file options are captured when each stream opens.

## Source parsing and writing

At initial open, the reader skips whitespace and leading `#` PEFF comment lines.
It retains the prior native extension accepting one UTF-8 BOM at stream start.
Headers begin with `>`; identifiers ignore leading spaces/tabs and CR, then stop
at a space, tab or LF. Descriptions remove tabs/CR but preserve ordinary spaces,
including additional spaces after the single identifier separator.

Sequences contain arbitrary bytes except space, tab, CR and LF, which are removed.
Consequently digits, asterisks, punctuation, modified sequence expressions,
semicolon lines, vertical tabs and non-ASCII characters are not interpreted as
an amino-acid alphabet. Each completed field must be valid UTF-8 for Rust strings.
The original `(ICPL:13C(6))` record therefore passes directly to `AASequence::parse`.

A new record is recognized only when the sequence loop consumes LF and the
immediately following byte is `>`. This preserves source edge cases: an indented
header belongs to the preceding sequence; `>a\n>b\nSEQ` is one record whose
sequence is `>bSEQ`; CR alone does not delimit headers. Semicolons are not comments.
Headers without a terminating LF and records without sequence data are errors.

The writer always emits `>identifier description\n`, including the separator
space when the description is empty. Sequence lines are **80 bytes**, with a final
LF after every nonempty chunk. Empty sequences write a header alone, although the
reader cannot read that as a complete standalone entry. Sequence whitespace is
written verbatim and removed if subsequently read. UTF-8 characters may straddle
an 80-byte line break; parsing rejoins the bytes before checking field encoding.

## Positioning and lifecycle

The reader leaves a recognized next header unconsumed. `position` returns its byte
position, or `None` for the source's `-1` EOF sentinel. Consuming a last record can
set EOF while successfully returning that record. `at_end` peeks and sets EOF;
on an empty stream it causes a later read to return none, whereas reading the
empty stream directly reports the source's first-record parse error. A leading
`#` comment ending at physical EOF also returns no records; the same comment
terminated by LF is followed by a first-record parse error.

Seeking clears EOF and native fused-error state, keeps cumulative read/record/work
counters, and does not repeat PEFF/BOM initialization. Seeking beyond the physical
stream length returns false without changing the logical position. A seek to the
exact end succeeds; the next direct read fails unless `at_end` first detects EOF.
Negative source stream positions are not representable by the `u64` API, and
actual seek failures return errors rather than source `setPosition`'s unchecked
success. Diagnostic line numbering restarts at one after a random seek.

A `FASTAFile` can read and write different files simultaneously. Aggregate `load`
and `store` use separate streams and leave those resident sessions intact.
Opening an already-open output is a checked error; `write_end` closes and flushes
it, and repeated `write_end` is harmless. Input/output operations before their
corresponding open return errors. Destruction makes a best-effort output flush;
call `write_end` or `finish` to observe I/O failures. The adapter adds no output or
flush replay after an observed streaming writer I/O error. An underlying writer
may have its own destructor behavior.

Path operations use **plain files**, as the source does; a `.fasta.gz` name does
not activate compression. Output permits FASTA and unknown extensions, rejecting
recognized other formats. Generic streams can be supplied by a caller's decoder,
which generally has no byte-seek operation. `store` preflights before truncating a
file, but transport failures can leave partial output; it is not an atomic rename.

Source progress labels, ranges and dispatch calls are retained: ordinary
`read_next` updates progress, its progress wrapper updates a successful read
again and ends on EOF, and final successful reads report `-1`. Inherited throttling
can suppress duplicate calls within one second. Aggregate load reports `0..1`;
store reports `0..record_count` and calls `next_progress` per entry. Repeated EOF
calls repeat `end_progress`, and failures can leave a started progress session
unfinished, as in the source; callers can explicitly manage recovery through
`progress`.

## Checked native boundaries

Limits are cumulative for each stream, including replay after seeks. Defaults:

| Limit | Default |
| --- | --- |
| Consumed input bytes | 512 MiB |
| Emitted output bytes | 512 MiB |
| Identifier + description + sequence bytes per entry | 16 MiB |
| Completed records | 1,000,000 |
| Work units | 4,000,000,000 |

Work charges input peeks, retained-byte append/UTF-8 scans, seeks, and conservative
output validation/emission visits. Entry limits are checked before growing field
vectors; no unbounded line is materialized. Read-ahead belongs to the supplied
buffered reader; file owners use the standard bounded `BufReader`. Aggregate
reads retain bounded complete records; streaming retains only the current entry.

Malformed or invalid UTF-8 reads fuse iteration until an explicit seek. The
caller-provided entry is unchanged on parser failure, and `load_into` commits only
a completely successful file. This improves the source's immediately cleared
aggregate destination. Writer preflight preserves the existing native rejection
of empty identifiers, whitespace/control characters in identifiers, `>` in an
identifier, and control characters in descriptions. The source's permissive
header injection is deliberately not reproduced. Sequence alphabet validation
was removed because it rejected legitimate source fixtures. Semantic validation
and limit failures occur before emitting the affected entry; aggregate writes
preflight every entry before emitting any bytes. Already emitted streaming
entries and external I/O effects cannot be rolled back.

## Reference evidence

[Tests](../tests/fasta.rs) check all **14 literal field assertions** across the
original **five-record fixture**, source PEFF/whitespace/asterisk cases, the source
seek-and-replay pattern and independently traced boundary cases. The second
record's description has no literal class-test assertion; round trips check its
preservation. Existing [format tests](../tests/formats.rs) retain BOM and malformed
header checks and now use the source's accepted sequence syntax.

[Provenance](../tests/data/fasta_provenance.json) records seven source hashes,
the unchanged-file comparison and both packaged fixture hashes. The field TSV
concatenates C++ literal strings only, preserving their original spelling and
starting line. No C++ was built or executed and no Rust-generated output was used
as a source golden. Progress callbacks, error/limit tests and platform-independent
seek failures are native or independently derived tests, not claimed upstream
runtime results.

Source: [FASTAFile.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/54a232fe2cae9c590d5c997fa49d20e7769860fb/src/openms/source/FORMAT/FASTAFile.cpp),
[FASTAFile.h](https://github.com/okohlbacher/OpenMS4-core/blob/54a232fe2cae9c590d5c997fa49d20e7769860fb/src/openms/include/OpenMS/FORMAT/FASTAFile.h),
[class tests](https://github.com/okohlbacher/OpenMS4-core/blob/54a232fe2cae9c590d5c997fa49d20e7769860fb/src/tests/class_tests/openms/source/FASTAFile_test.cpp).
