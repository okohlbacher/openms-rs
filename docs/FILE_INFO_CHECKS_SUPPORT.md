# FileInfo `-i`, `-d` and `-c` (A6-FILEINFO)

The three sections of `OpenMS::FileInfo::report_` that inspect a file rather
than summarise it: the indexed-mzML check, the detailed listing and the
corrupt-data check.

- Source: `FORMAT/FileInfo.h` and `FORMAT/FileInfo.cpp` at core `bc9cc12`,
  blocks `:827-846`, `:1779-1795`, `:1799-1848` and `:1851-1964`.
- Tool: `OpenMS4-topp/src/FileInfo.cpp` at topp `174b576`, `:118-123` and
  `:144-147`.
- Rust: `src/format/file_info/checks.rs`, wired into
  `src/format/file_info/report.rs` (the `-i` block and its early return) and
  `src/format/file_info/peaks.rs` (the two `-d` blocks and `-c`).
- Tests: `tests/file_info_checks.rs` (51) and the unit tests of `checks.rs`
  (11).
- Manifest: `tests/data/file_info_checks_provenance.json`.
- Oracle: `../oracle/a6-fileinfo`, 59 cases against the **Release** C++ build
  `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
  ibminode06, run twice and reproduced.

`FileInfo.h` stays `partial`: `-v` and the consensusXML, identification, FASTA,
mzTab, trafoXML and PQP branches belong to A7 and A8.

## API mapping

Every source member these three blocks touch, and where it went. `checks.rs`
exports no public item: the three writers are `pub(crate)`, as
`peaks::report` and `features::report` are, and the observable API stays
`FileInfo::run` with `Options`.

| Source | Rust | Notes |
|---|---|---|
| `Options::check_index` (`-i`) | `Options::check_index` | runs for every type, before the content |
| `Options::detailed` (`-d`) | `Options::detailed` | peak-file branch only |
| `Options::check_corrupt` (`-c`) | `Options::check_corrupt` | peak-file branch only |
| `report_` index block, `FORMAT/FileInfo.cpp:827-846` | `checks::write_index_check` | returns `false` for the source's early `return` |
| `report_` transition listing, `:1779-1795` | `checks::write_detailed_chromatograms` | |
| `report_` spectrum listing, `:1799-1848` | `checks::write_detailed_spectra` | |
| `report_` corrupt-data block, `:1851-1964` | `checks::write_corruption_check` | |
| `ValidationInfo::index_checked` | `ValidationInfo::index_checked` | set whenever `-i` ran |
| `ValidationInfo::index_valid` | `ValidationInfo::index_valid` | |
| `ValidationInfo::indexed_spectra` | `ValidationInfo::indexed_spectra` | `IndexedMzMLHandler::getNrSpectra()` |
| `ValidationInfo::indexed_chromatograms` | `ValidationInfo::indexed_chromatograms` | `getNrChromatograms()` |
| `CorruptionInfo::{performed,errors,warnings}` | `CorruptionInfo::{performed,errors,warnings}` | declared, never filled — see below |
| `DetailInfo::{performed,lines}` | `DetailInfo::{performed,lines}` | declared, never filled — see below |
| `Internal::IndexedMzMLHandler::openFile` + `parseFooter_` | `checks::parse_index` | reduced to what `parsing_success_` needs |
| `IndexedMzMLDecoder::findIndexListOffset` | `IndexedMzMLDecoder::find_index_list_offset` | already ported |
| `IndexedMzMLDecoder::parseOffsets` | `IndexedMzMLDecoder::parse_offsets` | already ported |
| `IndexedMzMLHandler::getParsingSuccess` | the `Option` `parse_index` returns | a handler that exists has parsed |
| `InstrumentSettings::NamesOfScanMode` | `ScanMode::name` | same 15 labels |
| `MSSpectrum::getDriftTimeUnitAsString` | `DriftTimeUnit::name` | `NamesOfDriftTimeUnit` |
| `Precursor::NamesOfActivationMethodShort` | `ActivationMethod::short_name` | |
| `Precursor::NamesOfActivationMethod` | `ActivationMethod::name` | |
| `MSExperiment::isSorted(false)` | `checks::nondescending` over the retention times | not the kernel predicate — see native difference 10 |
| `MSSpectrum::isSorted()` | `checks::nondescending` over the m/z values | not the kernel predicate — see native difference 10 |
| `ChromatogramSettings::getComment()` | not ported | no native field; always empty on this path — see below |
| `MSChromatogram::getName()` | `MSChromatogram::name` | |
| `MSChromatogram::getPrecursor().getMZ()` | `MSChromatogram::precursor.mz` | |
| `MSChromatogram::getProduct().getMZ()` | `MSChromatogram::product.mz` | |
| `TOPPFileInfo::outputTo_` `-i` guard, topp `:118-123` | `src/cli/tools/file_info.rs:184-190` | shipped by A5 |
| `TOPPFileInfo::outputTo_` exit, topp `:144-147` | `src/cli/tools/file_info.rs:219-222` | shipped by A5 |

## Preserved source conventions

- **Section order and every byte of it.** `-i` after the general header and
  before the content; the transition listing inside the chromatogram section,
  after the per-type counts; the spectrum listing after that section; `-c`
  last, before `-m`, `-p` and `-s`. The line texts, their leading and trailing
  spaces and their blank lines are the source's, including
  `" -- Detailed chromatogram listing -- "` with a space on each side and
  `"  activation methods: "` with a trailing one.
- **The early return of a failed `-i`.** `FORMAT/FileInfo.cpp:844` returns from
  `report_`, so the report ends with the failure text: no content, no `-m`,
  `-p` or `-s`, and not even the two trailing newlines of `:2443-2444`. Every
  other flag the caller set is skipped. The FileInfo tool turns that state
  into `ILLEGAL_PARAMETERS`, which is what upstream `TOPP_FileInfo_11` asserts
  with `WILL_FAIL 1`.
- **`-i` is not restricted to mzML in the library.** `report_` runs the check
  on whatever `-in` names; only the tool refuses a non-mzML input first. The
  port keeps that split, and refuses a branch it does not run only *after* the
  index check, so an unparsable index is reported in full even for a branch
  that is otherwise unsupported.
- **The m/z extent of the `-d` listing is the first and last stored peak**,
  `begin()` and `rbegin()`, not the smallest and largest. In a DTA or DTA2D
  scan, which no loader sorts, those differ.
- **An empty spectrum leaves the `m/z:` line unterminated.** The source writes
  the newline inside `if (!spectrum.empty())`, so the following
  `Precursors:  0` continues the same line. The port does not tidy it.
- **A blank line closes every precursor**, including the last.
- **`empty()` asks about the spectra alone**, so a file of chromatograms only
  reaches the transition listing and never the spectrum listing.
- **`-d` and `-c` are guarded inside the peak-file branch**, so a featureXML
  map ignores both and its report is identical to one written without them.
- **One name set for all three data-array kinds.** The source's single
  `std::map` means a float array and an integer array of the same name collide;
  the port uses one `BTreeSet` for the same reason.
- **A repeated value is reported once per repetition.** Both the retention-time
  and the m/z duplicate loops compare neighbours in the sorted vector, so a
  value stored three times gives two lines.
- **The `-c` header is written whatever the outcome**, so a clean file produces
  the header and nothing else.
- **Stream precision.** Nothing before these blocks changes it, so every
  `double` is written at the default 6 through
  `text_format::ostream_g`. Intensities are `float` promoted to `double`, as
  `operator<<` promotes them.

## Native differences

Each is documented at the item in `checks.rs` as well.

1. **`CorruptionInfo` and `DetailInfo` stay at their defaults.**
   `FileInfo.h:206-217` declares both with a `performed` flag and pre-rendered
   message lines, and `report_` never assigns to either at core `bc9cc12`: the
   text goes only into the stream. A caller using the library from Rust or
   pyOpenMS therefore cannot read the `-d` or `-c` findings as data, only as
   text. The port reproduces that rather than filling them, because filling
   them is a different API; the defect is recorded for `OpenMS_CPP_ISSUES.md`
   as `CPP-335`.
2. **An empty selected-reaction-monitoring chromatogram is refused.** The
   source reads `ms.front()` and `ms.back()` without checking, which is
   undefined on an empty container, so there is no behaviour to reproduce
   (wave 5's rule: refuse exactly where the source is out of bounds). The port
   returns `Error::InvalidValue` naming the chromatogram (upstream `CPP-336`).
   Unreachable through
   `FileInfo::run`: the mzML reader gives every chromatogram its points and
   `ChromatogramTools::convert_spectra_to_chromatograms` builds one point per
   source spectrum.
3. ~~**A NaN the source would sort is refused by `-c`.**~~ **CLOSED** under
   lead decision **D18**, which applies decision D16's machinery here. The
   entry is kept with its reasoning, because that is what makes the closure
   checkable.

   *What it said.* Such a value enters a `std::sort` whose comparator is then
   not a strict weak ordering, so the source is undefined; the port refused with
   `Error::InvalidValue` before the header was written.

   *Why that no longer follows.* Both calls are the unqualified
   `sort(v.begin(), v.end())` on a `std::vector<double>` — `ms1_rts` declared at
   `:1863` and sorted at `:1927`, `mzs` declared at `:1942` and sorted at
   `:1956`, with `using namespace std;` at `:47` and no comparator on either —
   which is exactly the call `crate::math::source_sort` reproduces comparison by
   comparison, tier 1 over 2,272 oracle inputs. D16 puts reproducing a
   permutation the standard leaves unspecified in scope, so `-c` now computes
   the answer the Release build computes instead of declining to.

   *That it is observable.* For MS1 retention times `{5.0, NaN, 5.0}` a
   three-element range is below `_S_threshold`, libstdc++ runs one
   `__insertion_sort` pass, every comparison involving the NaN is false and
   nothing moves — so the duplicate scan finds no equal neighbours and the
   Release build prints no duplicate line. A `f64::total_cmp` sort leaves
   `{5.0, 5.0, NaN}` and prints one. `a_nan_retention_time_keeps_the_release_builds_order`
   and `a_nan_mz_keeps_the_release_builds_order` pin both sorts.

   *What the port pays for it.* Nothing on ordinary data. The faithful sort
   allocates a permutation and an owned copy, and `-c` sorts every spectrum's
   m/z values, so `sort_as_the_source_does` takes it only where the permutation
   can be observed. Without a NaN, `operator<` on `f64` is a strict weak
   ordering, and two elements it calls equivalent are numerically equal and so
   bit-identical — with the single exception of `-0.0` against `+0.0`. Where
   every equivalence class is a set of bit-identical values the output
   *sequence* is a function of the multiset alone, so the permutation is
   unobservable and any correct sort writes the same bytes. The guard is
   therefore one O(n) pass for a NaN and for both zero spellings; everything
   else sorts in place with `f64::total_cmp` and allocates nothing.
   `both_paths_agree_wherever_the_fast_path_is_taken` runs both paths over the
   same inputs and asserts bit-identical output, rather than leaving the
   argument as an assertion.

   Infinities were never refused — `<`, `>` and `==` are defined on them.
   Unreachable through `FileInfo::run` either way: every loader on this path
   validates its coordinates, so the shapes above are built by hand in the unit
   tests.
4. **The chromatogram comment is always empty.**
   `ChromatogramSettings::getComment()` has no counterpart in
   `MSChromatogram`, and nothing in the pinned source calls `setComment` on a
   chromatogram — not `MzMLHandler`, not `ChromatogramTools`. Every
   chromatogram FileInfo can see therefore carries the empty default in C++ too,
   so the transition line ends in the same two spaces. `MzXMLHandler`,
   `MzDataHandler` and `XMassFile` set a *spectrum* comment, which this line
   does not print.
5. **No `std::cerr` diagnostics.** `IndexedMzMLDecoder` prints
   `findIndexListOffset Error: ...` and a dump of the searched bytes to
   `std::cerr` when it finds no footer (`IndexedMzMLDecoder.cpp:198-200`, the
   `else` of the match at `:183`). The port writes nothing there. No
   report line and no exit code depends on it; the oracle case
   `i_truncated_index` shows the two agree on both.
6. **Spectrum and precursor numbers do not wrap.** The source counts in `UInt`,
   which wraps above 4294967295; the port counts in `usize`. No loaded
   experiment reaches that count on any supported target.
7. **The footer is read before the first `-i` line is written.** The source
   writes `Checking mzML file for valid indices ... ` and then parses, letting
   an exception carry the line away. Nothing observes the difference: a run
   that cannot read the footer returns the error and no report.
8. **Without the `mzml` feature, `-i` is refused** with `Error::Unsupported`;
   there is no index decoder in that build.
9. **`Error::Io` and `Error::Parse` from the decoder propagate; every other
   decoding failure is the invalid-index outcome.** That is the source's own
   split, and the two errors carry its two exception classes.
   `Exception::FileNotFound`, `FileNotReadable` and `IOException` become
   `Error::Io` (`IndexedMzMLDecoder.cpp:79-87` in `parseOffsets`, `:152-160` in
   `findIndexListOffset`); the `Exception::ConversionError` of a footer offset
   that is not a non-negative 63-bit integer becomes `Error::Parse` (`:51-61`,
   the two throws of `stringToStreampos`) — an over-64-bit `<indexListOffset>`
   reaches it, and the tool then prints no report at all on either side.
   Everything else is the source's `return -1`: an offset outside the file
   (`:97-99`), a failed allocation (`:115-118`) and a malformed index (`:332`).
   The port's index byte and offset-count ceilings land in the
   failed-allocation bucket, which is where the source puts an index too large
   to hold.

10. **The two sortedness tests do not use the kernel predicates.**
    `MSExperiment::is_sorted` and `MSSpectrum::is_sorted` in this crate refuse a
    non-finite coordinate and therefore report a container holding an infinity
    as unsorted. The source compares neighbours with `>` alone, for which an
    infinity is ordinary, so `checks::nondescending` applies the source's
    comparison and the port writes no line the C++ would not. A NaN needs the
    same treatment for the same reason, and since item 3 closed it needs it in
    earnest: a NaN is never greater under `>`, so it never makes a container
    unsorted, while the kernel predicates would report one. Every other caller
    in the crate keeps the kernel predicates. Caught by the unit test
    `an_infinite_mz_is_checked_like_any_other`, which asserted only the
    duplicate line until it was strengthened to assert the absence of the
    unsorted one.

11. **A file shorter than 1023 bytes is read where the source reads nothing.**
    `findIndexListOffset` allocates `new char[buffersize + 1]`, seeks
    `-buffersize` from the end and reads (`IndexedMzMLDecoder.cpp:165-168`).
    On a shorter file the seek fails, the read writes nothing and the regex at
    `:179-181` searches a buffer nothing initialised, so what the source finds
    there is *indeterminate*: wave 5's rule does not ask for it to be
    reproduced. The port reads `min(length, 1023)` bytes and searches the whole
    file, which `find_index_list_offset_with_buffer_size`
    (`src/format/indexed_mzml.rs:74`) has always done. Measured: the Release
    build reports no index for the 967-byte `index_window_below.mzML` and exits
    `ILLEGAL_PARAMETERS`; the port finds its index and writes the report the
    C++ writes for `index_window_above.mzML`, the same file padded to 1217
    bytes, which both implementations agree on byte for byte. The two runs of
    the below-window case differ *only* in the bytes the source dumps to
    `std::cerr` (`:198-200`) — the uninitialised buffer itself — and the
    manifest records both raw hashes under `indeterminate_stderr`. **Owner:**
    `src/format/indexed_mzml.rs`, which every index reader in the crate shares,
    not this package. The departure is that package's own reviewed decision and
    is already recorded in `docs/INDEXED_MZML_SUPPORT.md` ("small files search
    their full contents"; "checked small-file reads … replace unchecked or
    platform-dependent source cases"). What is new here is that `-i` prints its
    consequence, and that the two oracle cases pin where it starts.

12. **An index whose first child is an `<offset>` keeps the offset the source
    drops.** `domParseIndexedEnd_` walks the children of each `<index>` as
    `iter = getFirstChild(); while (iter != lastChild) { iter = getNextSibling(); … }`
    (`:280-282` sets `iter`, and the walk itself is `:290-293`), advancing
    before it reads, so the first child is never looked at. A text node —
    any whitespace — immediately after the opening `<index …>` tag takes that
    place and the walk loses nothing, which is why every index an OpenMS
    writer produces parses. Only that one position matters: whitespace between
    the offsets or before `</index>` does not save the first offset, so the
    affected class is "first child is an `<offset>`", not "written without
    whitespace". Measured: the Release build counts one spectrum in the
    two-offset `index_offsets_unspaced.mzML`, the port counts two. Reproducing it would
    mean re-implementing the source's DOM walk in `parse_offsets`, which every
    other index reader in the crate calls — and it would break random access,
    because a dropped offset is a record that cannot be found. **Owner:**
    `src/format/indexed_mzml.rs`, not this package. That package already states
    the decision in `docs/INDEXED_MZML_SUPPORT.md` ("the native parser corrects
    the source DOM sibling loop that skips an offset when it is the first child
    without preceding whitespace"); the source defect is recorded for
    `OpenMS_CPP_ISSUES.md` as CPP-337.

## Checked boundaries and evidence

### Resource bounds

- `-c` reserves the MS1 retention-time buffer and, per spectrum, the m/z buffer
  with `try_reserve`, and refuses with `Error::InvalidValue` rather than
  aborting. Both are bounded by the experiment already in memory; the m/z
  buffer is reused across spectra.
- `-i` inherits `IndexReadLimits`: 1 MiB of footer, 16 MiB of index, 1,000,000
  offsets, 64 KiB of identifiers. Exceeding any of them makes the index invalid,
  as a failed allocation does in the source.
- `-d` allocates nothing beyond the report itself.

### Evidence

**Tier 1, executed differential.** The Release C++ FileInfo at all three pins
this port follows ran 59 cases on ibminode06, twice, reproduced
(`../oracle/a6-fileinfo/manifest.json`, `"reproduced": true`). 38 cases have
both reports compared byte for byte in `tests/file_info_checks.rs`, one
(`i_truncated_index`) its text, one (`i_window_below`) its report against the
padded file's because that is exactly the native difference, and the rest carry
exit codes, counts and refusals.

Three masks, all of them on stdout or stderr and never in `-out` or `-out_tsv`:
the `FileInfo took ...` footer; the `ProgressLogger`'s `-- done [took ...] --`
lines, which only the case that reproduces `TOPP_FileInfo_19` verbatim emits,
because that upstream test passes no `-no_progress`; and, for
`i_window_below` alone, the bytes `findIndexListOffset` dumps to `std::cerr`
after `Maybe this is not a indexedMzML.`, which for a file shorter than its
search window are the uninitialised buffer and differ from run to run. The
mask is listed per case in `indeterminate_stderr` with both runs' raw hashes;
every case whose dump is the file's real tail — `i_invalid_11`,
`i_truncated_index` — is left unmasked and reproduces exactly.

Covered by at least one compared case: a valid index with and without the other
flags; an absent index, a truncated index, an index on an empty mzML and on a
DTA; the transition listing on pure and mixed SRM input; the spectrum listing
on mzML, DTA, DTA2D, FAIMS and an experiment with no spectrum, with and without
an ion-mobility line, with and without precursors, and over an empty spectrum;
and every branch of `-c`:

| `-c` line | reached by |
|---|---|
| retention times not sorted | `c_scans`, `c_dta2d` |
| MS level 0 | `c_scans` |
| no peaks in spectrum | `c_scans`, `c_9_strict` |
| duplicate data-array name | `c_arrays`, `c_arrays_mixed` (C++ side), unit test (port side) |
| duplicate MS1 retention time | `c_scans`, `c_dta2d` |
| peak m/z not sorted | `c_unsorted_dta`, `c_unsorted_dta2d` |
| negative peak intensity | `c_peaks`, `c_unsorted_dta`, `c_unsorted_dta2d` |
| duplicate peak m/z | `c_peaks`, `c_unsorted_dta`, `c_unsorted_dta2d` |
| nothing (clean) | `c_clean`, `c_srm`, `c_faims`, `c_empty`, `c_dta` |

**Tier 1, retained upstream definition.** `TOPP_FileInfo_11`
(`topp/CMakeLists.txt:908-909` at test-data `0cb15f2`, `WILL_FAIL 1`) and
`TOPP_FileInfo_19` (`:928-930`) are reproduced, the first including its exit
code. `TOPP_FileInfo_12` (`:910-911`) and `TOPP_FileInfo_6` (`:896-898`) are
not, for the reader reasons below. The definitions are in the **test-data**
package, not in the topp pin: `174b576`'s own `CMakeLists.txt` is 63 lines and
names no FileInfo test.

**Tier 4.** The two undefined places above, the aggregates the source never
fills, the `-v` refusal that remains, and the infinity that is *not* refused.

### What is out of reach, and why

None of these is a gap in this package. Each is recorded twice with its owner:
in `known_gaps` of the reference manifest
`tests/data/file_info_checks_provenance.json`, and in `known_gaps` of the
oracle manifest `../oracle/a6-fileinfo/manifest.json`, where each entry also
carries its kind (`unreachable` or `divergence`) and the cases that show it.
Native differences 11 and 12 are in both lists too, because their owner is
`src/format/indexed_mzml.rs` rather than this package.

- **The `-c` duplicate-data-array-name line cannot be reached through the
  port's mzML reader.** `Record::check_array_kind` (`src/format/mzml.rs:1038`)
  refuses a repeated auxiliary array name with
  `Error::Parse("duplicate auxiliary array name")`, where the C++ reader loads
  the file and leaves the report to `-c` — which is what the flag exists for.
  The oracle's `c_arrays` and `c_arrays_mixed` record what the C++ writes, a
  test asserts the port's refusal, and a unit test of `checks.rs` covers the
  port's rendering of the line on an experiment built in memory. **Open for
  the lead:** a flag whose purpose is to report corrupt files cannot report
  this class of corruption while the reader refuses it.
- **The `-c` MS-level-0 line is reachable only for an optical scan mode.**
  `MSSpectrum::validate_scalars` (`src/kernel.rs:810-821`) refuses `ms_level`
  0 unless the scan mode is `ElectromagneticRadiation`, `Emission` or
  `Absorption`; the C++ loads an MS-level-0 *mass* spectrum.
  `corrupt_scans.mzML` gives its level-0 scan `MS:1000806 absorption spectrum`
  so both implementations load it and the line is exercised. Same shape of
  finding as the one above, at the kernel invariant.
- **`TOPP_FileInfo_12`'s own input cannot be loaded**: its `charge array` is
  stored as 64-bit float, which the strict mzML reader refuses
  (`canonical auxiliary array binary type`); A4 already records this. Its
  *index* is checked here and parses with the three spectra and no chromatogram
  the C++ reports; the loadable cases use the core `IndexedmzMLFile_1` fixture,
  whose `-i` run exits 0 in both implementations.
- **`TOPP_FileInfo_6`** runs `-d` on mzData, a peak format with no native
  reader on this path. The branch is refused as before this package; it belongs
  to A8.

### Fixtures

`corrupt_peaks.mzML`, `corrupt_scans.mzML`, `corrupt_arrays.mzML`,
`corrupt_arrays_mixed.mzML`, `clean_pair.mzML`, `unsorted_peaks.dta` and
`unsorted_peaks.dta2d` are generated by
`../oracle/a6-fileinfo/scripts/make_corrupt_fixtures.py`, whose rules the
manifest records. `index_window_below.mzML`, `index_window_above.mzML` and
`index_offsets_unspaced.mzML`, the three that pin native differences 11 and 12,
are generated by `scripts/make_index_window_fixtures.py`: the smallest mzML
payload both readers load, with an index written by hand so that the file's
size and the whitespace inside `<index>` are under the test's control. No value
in any of them comes from a Rust output: they are hand-built valid mzML, DTA
and DTA2D whose degeneracy is in the data, and every expected report is the
Release C++ output.
