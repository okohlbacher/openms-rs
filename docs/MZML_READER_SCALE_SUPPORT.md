# mzML reader: size-derived ceilings and timestamp sentinels

Two benchmark blockers are fixed here. Both were found by the OpenMS4 smoke
benchmark on `ibminode05`/`ibminode06` (results directory
`/ceph/ibmi/abi/oliver/bench/openms4/results/2026-09-14-smoke`), where **60 of
60 Rust tool executions on the real inputs failed before doing any work**:

1. `run/@startTimeStamp="-infinity"`, which ProteoWizard writes when the vendor
   file carries no acquisition date, made every tool exit 6 on the PXD001819
   `50amol_R1.mzML` input (`src/data_structures/datetime.rs:435`,
   `invalid DateTime input or calendar fields`).
2. The reader's fixed cumulative ceilings rejected every real input: 512 MiB of
   XML, 10 M peaks, 20 M array elements, 512 MiB of parameter storage.

## 1. Timestamp sentinels

### What the C++ does, executed

`../oracle/mzml-reader-scale/datetime_sentinel_cpp.sh` ran the C++ **Release**
build (core `bc9cc12`, prefix
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`) on a
three-spectrum slice of the benchmark input, once per timestamp spelling, with
`FileInfo`, `FileConverter` and `MzMLSplitter`. Results are in
`tests/data/mzml_reader_scale/cpp_timestamp_oracle.tsv`:

| `startTimeStamp` | `FileInfo` | `FileConverter` | stderr | written back as |
|---|---|---|---|---|
| `-infinity` | exit 0 | exit 0 | `Non-fatal error while loading '…': DateTime conversion error of "-infinity"` | no `startTimeStamp` |
| `infinity` | exit 0 | exit 0 | the same, with `"infinity"` | no `startTimeStamp` |
| `not-a-date-time` | exit 0 | exit 0 | the same, with `"not-a-date-time"` | no `startTimeStamp` |
| (empty) | exit 0 | exit 0 | silent | no `startTimeStamp` |
| `2014-01-01T00:00:00` | exit 0 | exit 0 | silent | `2014-01-01T00:00:00` |

A `processingMethod` completion time (`MS:1000747`) spelled `not-a-date-time`
was dropped without any Release-build output, and the written file carries only
the converter's own completion time.

The Release build's silence there is a build switch, not source behaviour.
`XMLHandler::warning` (`XMLHandler.cpp:88-107`) writes to `OPENMS_LOG_WARN`
only under `OPENMS_ASSERTIONS` and to `OPENMS_LOG_DEBUG` otherwise, and
`cvParamToValue` reaches that warning for an `xsd:dateTime` CV value
`DateTime::set` rejects. Executed on the **Debug** product SDK (the development
oracle, core `4fdec46`) with
`../oracle/mzml-reader-scale/debug_sdk_completion_time.sh`, log
`debug_sdk_completion_time.log`:

| completion time | `FileInfo` | `FileConverter` | stderr | written back as |
|---|---|---|---|---|
| `not-a-date-time` | exit 0 | exit 0 | `While loading '…': The CV term 'MS:1000747 - completion time' used in tag 'processingMethod' must be a valid date. The value is 'not-a-date-time'.` | only the converter's own `MS:1000747` |
| `-infinity` | exit 0 | exit 0 | the same, with `'-infinity'` | only the converter's own `MS:1000747` |
| `2014-01-01T00:00:00` | exit 0 | exit 0 | silent | kept, plus the converter's own |

This port has one warning level for both dropped values, so it reports the
completion time the way the Debug build does and the run timestamp the way both
builds do.

### Where the behaviour belongs

**Not in `DateTime`.** Source `DateTime::set` throws `Exception::ParseError`
for all three sentinels (`DateTime.cpp:332-335`), exactly as
`DateTime::parse` returns `Error::InvalidValue` here: the type has no notion of
a sentinel. The leniency lives at the mzML boundary, in two different places:

| Source | Rust |
|---|---|
| `XMLHandler::asDateTime_` (`XMLHandler.h:359-377`) trims, truncates to 19 characters, calls `DateTime::set`, catches `ParseError`, calls `error(LOAD, …)` and returns an unset `DateTime`. Used for `run/@startTimeStamp` (`MzMLHandler.cpp:1236`) | `mzml_header::read::run_timestamp`, under `ReadOptions::source_invalid_timestamps` |
| `XMLHandler::cvParamToValue` (`XMLHandler.cpp:232-243`) rejects an `xsd:dateTime` CV value that `DateTime::set` refuses, warns and returns `DataValue::EMPTY`, so `handleCVParam_` returns before `setCompletionTime` (`MzMLHandler.cpp:1539`, `:3246-3249`) | `mzml_header::read::processing_param`, under the same switch |
| `error(LOAD, …)` prints `Non-fatal error while loading …` on `OPENMS_LOG_ERROR` (`XMLHandler.cpp:71-87`) | one line per dropped value on the crate's warning log stream |

### The default is lenient, unlike the dangling-reference switch

The executed source is lenient here **unconditionally**: there is no strict mode
to opt into, and no tool passes one. `ReadOptions::source_invalid_timestamps`
therefore defaults to `true`, and `ReadOptions::source()` only adds
`source_dangling_references` on top.

That is a deliberate exception to decision D10 (*library defaults strict, tool
paths opt into source compatibility*), made because the strict default did not
cost a value, it cost the file: `FileHandler::load_experiment` passes
`ReadOptions::default()`, so with the switch off every TOPP tool exited 6 on the
PXD001819 input — the smoke benchmark's 60-of-60 failure. The two switches
differ in what the strict policy protects. A dangling `softwareRef` loses a
reference the file still contains, so refusing it protects data the reader could
otherwise silently drop, and only a tool that wants source behaviour turns it
off. An unparseable `startTimeStamp` carries no date at all: `-infinity` *is*
ProteoWizard's spelling of "this vendor file has no acquisition date", so
leniency loses nothing the file ever had. A caller who wants the strict policy
sets `source_invalid_timestamps: false`, which is what the regression test does.

### Native differences

- **The strict policy is available, and off.** The source has no strict mode;
  this port has one, unused by default.
- **The raw text is dropped, not retained.** The source keeps
  `mzml_start_time_stamp` when the attribute is longer than 19 characters, even
  when it did not parse, and its writer then omits it
  (`MzMLHandler.cpp:5209-5229`). This port's writer refuses that key without a
  valid date-time, so a rejected timestamp keeps nothing. The written run
  element is identical to the C++ one in both cases.
- **Milliseconds.** The source truncates to 19 characters before parsing; this
  port parses the whole text and keeps represented milliseconds (CPP-027,
  already documented in `MZML_HEADER_SUPPORT.md`). Only text the source would
  have parsed *differently* is affected, never a sentinel.
- **Warning.** Each dropped value writes one line
  (`Warning: mzML run startTimeStamp '-infinity' is not a date-time; …`). The
  source prints a non-fatal error for the run attribute in every build, and for
  the completion time a warning that only a Debug build shows (table above);
  this port has one level for both.
- **Empty text** is silently unset by default and an error under the strict
  policy; the source is silent in both cases.
- **Not covered.** Other `xsd:dateTime` CV terms that could appear in a header
  (`MS:1002435` data processing start time) still fail the read under both
  policies, as do non-date CV values of the wrong type, which the source also
  drops with a warning. That is a separate leniency, not part of this lane.

## 2. Size-derived ceilings

### The rule

Every cumulative allowance is now `floor + units * (consumed / per_bytes)`,
credited as the reader consumes XML bytes, and applies **together with** an
absolute ceiling:

```text
charge fails  <=>  charge > min(absolute ceiling, floor + rate * consumed)
```

`per_bytes` is 1 for the quantities that grow with every byte (peaks, decoded
bytes, elements, parameters, work) and larger for the counts that each need
their own start tag: one record per 512 consumed bytes, one binary array per
256, one parameter group per 4,096. `Allowance::after` recomputes the ceiling
from the total consumed rather than crediting each event, so a rate below one
unit per byte is exact however the input happens to be chunked.

`src/format/mzml_scaling.rs` holds `Allowance`, `InputScaling` and the `Ledger`
that reconciles the reader's plain `usize` counters after every XML event, so
the charge sites are unchanged. The absolute ceilings keep their `ReadOptions`
fields and now default to unbounded; a caller that sets one keeps it exactly.
`InputScaling::fixed()` restores the former fixed ceilings, which is what the
regression test uses to show that they reject a realistic document.

Quantities that need their own XML start tag — records, binary arrays,
parameter groups — are bounded by the input size too, but a record is the most
memory-amplifying thing a document can declare per XML byte (an empty
`<spectrum>` of about 240 bytes costs roughly a kilobyte of retained record),
so they get their own size-derived allowances rather than the input size alone.

### Where the defaults come from

Measured with an instrumented reader (no ceilings) on every staged benchmark
input on `ibminode06`, recorded in
`../oracle/mzml-reader-scale/probe_unlimited_ibminode06.log`
(the patch is `probe_instrumentation.diff`). The ratios are charge per consumed
XML byte, the maximum over all prefixes of the six mzML inputs (0.5 GB to
2.3 GB, uncompressed, zlib and mixed):

| allowance | floor (the former fixed ceiling) | per byte | largest measured | headroom |
|---|---|---|---|---|
| `peaks` | 10,000,000 | 8 | 0.087 | 92x |
| `array_bytes` | 512 MiB | 64 | 1.19 | 54x |
| `array_elements` | 20,000,000 | 16 | 0.17 | 94x |
| `records` | 1,000,000 | 1 per 512 bytes | 1 per 3,951 bytes | 7.7x |
| `arrays` | 1,000,000 | 1 per 256 bytes | 1 per 3,559 bytes | 13.9x |
| `param_groups` | 100,000 | 1 per 4,096 bytes | 1 per 1.2 GB | enormous |
| `params` | 10,000,000 | 4 | 0.007 | 570x |
| `param_bytes` | 512 MiB | 256 | 19.4 | 13x |
| `metadata_work` | 50,000,000 | 64 | 0.11 | 600x |
| `metadata_bytes` | 256 MiB | 64 | 0.11 | 600x |
| `selection_work` | 500,000,000 | 512 | not instrumented; ~180 units per selected peak, at most 0.37 peaks per byte | ~8x |
| `selection_bytes` | 256 MiB | 64 | 24 bytes per selected peak | ~7x |
| `count_work` | 50,000,000 | 16 | one unit per XML event | large |

The record, array and group densities, and the largest single array, were
counted directly from the staged inputs on `ibminode06`
(`../oracle/mzml-reader-scale/array_scan_ibminode06.log`), independently of the
instrumented reader:

| input | size | records | bytes/record | arrays | bytes/array | groups | largest array |
|---|---|---|---|---|---|---|---|
| `20120210_…_int20000_filtered.mzML` | 5,203,440 | 731 | 7,118 | 1,462 | 3,559 | 0 | 539 |
| `Zeitz_SIP_13-II_020_small_knubbel.mzML` | 7,622,518 | 1,929 | **3,951** | 1,538 | 4,956 | 0 | 289 |
| `UK222_picked.mzML` | 547,283,125 | 40,857 | 13,395 | 81,714 | 6,697 | 0 | 8,174 |
| `50amol_R1.mzML` | 1,197,928,823 | 43,746 | 27,383 | 87,492 | 13,691 | 1 | 43,745 |
| `20100219_SvNa_SA_Ecoli_preccorrected.mzML` | 1,506,847,945 | 34,895 | 43,182 | 69,790 | 21,591 | 0 | 34,894 |
| `UK222.mzML` | 2,317,975,830 | 40,857 | 56,733 | 81,714 | 28,366 | 1 | **53,824** |

The floors alone (1,000,000 records, 1,000,000 arrays) already hold every one of
these inputs 12 to 24 times over; the rates only matter above a million records,
which no real file approaches. A document of 2,000,000 empty records in 473 MiB
of XML (236 bytes each) is refused, where the unbounded ceiling accepted it at
4.4 times the input size in resident memory.

The margins also cover encodings the benchmark inputs do not use. The densest
realistic encoding of a peak is Numpress linear plus short-logged-float with
zlib, about 2.7 base64 characters per peak, which is 0.37 peaks and 6 decoded
bytes per input byte — still 20x inside `peaks` and 10x inside `array_bytes`.

Two ceilings stay absolute because they are per item, not cumulative:

- `max_xml_bytes`, now 1 TiB. It only stops an unbounded stream; a file ends on
  its own, and every other allowance scales with what was consumed.
- `max_array_bytes`, **unchanged at 64 MiB** per array: 8 million `f64` values.
  The largest single array in any benchmark input is the `UK222.mzML` TIC
  chromatogram, 53,824 elements and 431 KB decoded (table above), so the ceiling
  is 155 times the largest real array. It bounds the one transient value vector
  a decoded array allocates, and it is the only thing that stops a small
  document whose few arrays each inflate to hundreds of megabytes: a 203 KiB
  document with one spectrum, `defaultArrayLength="9999999"` and two
  zlib-compressed 64-bit arrays stays inside the cumulative peak, element and
  byte allowances a document of any size is granted. An eightfold increase to
  512 MiB, tried in the first round of this lane and reverted here, let exactly
  that document through at 1,547 times its own size in resident memory.

### Two allowances that were cumulative and should not have been

- **The Numpress coder's work and allocation allowance** was one budget for the
  whole document (500,000,000 work units, 512 MiB), and
  `numpress_coder::decode_text` charges 64 units per encoded byte. Any Numpress
  mzML beyond about 8 MB of encoded payload was therefore unreadable. The
  reader now builds the coder's `Work` per array, from the coder defaults plus a
  multiple of that array's own encoded text and declared value count. The
  cumulative element and byte allowances still bound the document. The private
  test `numpress_transport::tests::decoder_counters_span_arrays_and_fail_before_resetting`,
  which asserts that a caller-supplied `Work` spans arrays, is unchanged and
  still passes: only the reader's choice of lifetime changed.
- **The scientific selection budget** (`LoadOptions::max_selection_*`) is still
  cumulative, including the per-record scratch that
  `selection_budget_is_cumulative_across_records_and_independent_of_binary_caps`
  pins, but its ceilings are now size-derived, so the
  FeatureFinderCentroided load path reads a 2.3 GB run.

### What still refuses an adversarial document

`tests/mzml_reader_scale.rs` checks all of these with the defaults:

- a record declaring 1,000,000,000 points with a three-point payload
  (`peak count exceeds configured limit`, before any allocation);
- a zlib payload inflating 4,000-fold, declared as 50 million and as 9 million
  points (refused by the peak allowance and by the decoded-byte allowance);
- a 203 KiB document with one spectrum, `defaultArrayLength="9999999"` and two
  zlib-compressed 64-bit arrays, which no cumulative allowance stops: the
  per-array ceiling refuses it;
- a 2 MB document whose every spectrum expands a 20,000-parameter group,
  charging about 17,000 storage bytes per input byte against an allowance of
  256 (refused after about a megabyte, whatever the document's total size);
- a tightened `records`, `arrays` or `param_groups` allowance refusing a
  document the defaults read, at the count it tightens;
- 200 levels of nesting (the reader's fixed 128-level limit, unchanged).

The two documents too large for the test suite were run through the shipped
`MapNormalizer` on `ibminode06` instead, with its own default load options
(`../oracle/mzml-reader-scale/tool_probe_ibminode06.log` and
`tool_probe2_ibminode06.log`):

| document | size | result | peak RSS |
|---|---|---|---|
| one spectrum, `defaultArrayLength="9999999"`, two zlib 64-bit arrays | 219,069 B | refused, `binary array exceeds configured byte limit` | 4,096 KiB |
| 2,000,000 minimal empty records, 252 bytes each | 503,778,613 B | refused, `record count exceeds configured limit` | 1,723,392 KiB |
| 1,000,000 of the same records | 250,778,613 B | read (the tool then reports no intensities) | 884,736 KiB |

The zlib bomb is refused before any allocation. The record case shows what the
allowance buys and what it does not: the reader stops at the ceiling rather than
before it, so 2,000,000 records still cost 1.7 GB before the refusal, 3.4 times
the document. That constant is the record rate: a finer one would refuse real
files, and the floor deliberately admits a million records from a 240 MB
document. What changed is that the count is bounded at all — the first round of
this lane read all 2,000,000.

### Measured on the real inputs

`ibminode06`, release build, inputs staged on node-local `/scratch`, one
process per input (`../oracle/mzml-reader-scale/hpc_scale_06.sh` and
`hpc_scale_06_round2.sh`, logs `hpc_scale_ibminode06.log` for the first round
and `hpc_scale_ibminode06_round2.log` for the current defaults). "read" is the
`ReadOptions::source()` stream read; the process also counts the same file again
with `load_size_with_options`, which is the rest of the wall time. The node
carried a foreign load of 30 to 40 runnable processes throughout, so the times
are indicative and the peak RSS is not.

| input | size | spectra | peaks | read | process wall | peak RSS |
|---|---|---|---|---|---|---|
| sanity profile + sanity centroid | 5.2 MB + 7.6 MB | 730 + 1,929 | 252,018 | 0.3 s + 0.1 s | 0.8 s | 31,744 KiB |
| `UK222_picked.mzML` | 547 MB | 40,856 | 22,776,198 | 6.2 s | 12.4 s | 747,520 KiB |
| `50amol_R1.mzML` (zlib) | 1.20 GB | 43,745 | 88,434,492 | 15.3 s | 35.5 s | 1,946,624 KiB |
| `20100219_SvNa_SA_Ecoli_preccorrected.mzML` | 1.51 GB | 34,894 | 85,315,432 | 11.3 s | 25.7 s | 1,653,760 KiB |
| `UK222.mzML` | 2.32 GB | 40,856 | 197,765,338 | 17.8 s | 38.7 s | 3,491,848 KiB |
| `UK222.mzML`, FeatureFinderCentroided load options | 2.32 GB | 6,911 (MS1) | 96,433,834 | 16.7 s | 17.1 s | 1,972,736 KiB |

For comparison, the C++ Release `BaselineFilter` peaked at 3.2 GB on
`UK222.mzML` in the smoke benchmark, and `MapNormalizer` at 1.5 GB on
`50amol_R1.mzML`. Every spectrum count matches `inputs/MANIFEST.json`. Peak RSS
is within 1% of the first round's, which used a 512 MiB per-array ceiling and no
record or array allowance: no real input comes near either bound, which is the
point of the table above.

## Evidence

- **Tier 1, executed C++**: the timestamp table above, produced by
  `../oracle/mzml-reader-scale/datetime_sentinel_cpp.sh` on the Release build,
  and the completion-time table, produced by `debug_sdk_completion_time.sh` on
  the Debug product SDK (log `debug_sdk_completion_time.log`). Input, per-case
  inputs and written outputs are hashed in
  `tests/data/mzml_reader_scale_provenance.json`.
- **Tier 1, pinned source**: `XMLHandler::warning` writing to `OPENMS_LOG_WARN`
  only under `OPENMS_ASSERTIONS` (`XMLHandler.cpp:88-107`), which is why the
  Release build is silent for the dropped completion time and the Debug build is
  not.
- **Tier 4, measured**: the per-byte ratios, the direct record/array/group
  counts (`array_scan.sh`, log `array_scan_ibminode06.log`) and the HPC scale
  runs. Every staged benchmark mzML input is read by the `#[ignore]`d tests in
  `tests/mzml_reader_scale.rs`; spectrum counts are compared against
  `inputs/MANIFEST.json`, which was produced by a Python header scan, not by
  this port. Wall time and peak RSS per input are in
  `../oracle/mzml-reader-scale/hpc_scale_ibminode06_round2.log`, and the
  tool-level runs in `tool_probe_ibminode06.log` and
  `tool_probe2_ibminode06.log`.
- **Tier 4, Rust-only**: the synthetic-document tests, the adversarial
  documents, the explicit-ceiling precedence and the `Ledger` unit tests in
  `src/format/mzml_scaling.rs`.

`/usr/local/bin/cc` on `ibminode06` is an unrelated shell script that shadows
the real compiler, so `rustc` there links without producing a binary and every
build fails later with `build-script-build (never executed)`. The measurement
build passes `-C linker=/usr/bin/cc` around it
(`../oracle/mzml-reader-scale/build_06_round2.sh`); the gates run on `dax`,
which is unaffected.

## Tool wiring

`FileHandler::load_experiment` and `load_experiment_with_options` pass
`mzml::ReadOptions::default()`, which is why the timestamp default had to be the
source behaviour: both benchmark blockers are now fixed through that path, with
no change to `src/format/file_handler.rs`, which this lane does not own. The
executed proof is in `../oracle/mzml-reader-scale/tool_probe2_ibminode06.log`:
on the three-spectrum slice of the PXD001819 input whose `run/@startTimeStamp`
is `-infinity`, the shipped `MapNormalizer` and `DTAExtractor` now exit 0, print
`Warning: mzML run startTimeStamp '-infinity' is not a date-time; …`, and write
the run element back without `startTimeStamp`, exactly as the C++ oracle did.
The same binaries built from the previous round exited 6 with `invalid DateTime
input or calendar fields`.

On the **full** 1.2 GB input (`tool_probe_ibminode06.log`), `MapNormalizer` now
reads all 43,745 spectra in 18 s at 1,946,624 KiB and fails afterwards in the
**writer**, with `data array description resource limit exceeded`
(`mzml_header/write.rs`), which is a separate lane's blocker (the writer budget
that is cumulative where it should be per record). Both reader blockers are out
of the way of that path.

What still needs the integrator is the **other** switch:
`source_dangling_references` stays `false` by default (D10), so a tool that must
read a file with a dangling `softwareRef` still needs `ReadOptions::source()` on
the `FileHandler` load paths. That wiring belongs to the lane that owns
`file_handler.rs`.
