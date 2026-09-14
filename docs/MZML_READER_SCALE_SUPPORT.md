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
was dropped silently, and the written file carries only the converter's own
completion time.

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

`ReadOptions::source()` turns on this switch together with
`source_dangling_references`, following decision D10: **library defaults stay
strict, tool paths opt into source compatibility.** The default still refuses
an unparseable timestamp, because the value cannot be kept.

### Native differences

- **Strict by default.** The source has no strict mode.
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
  (`Warning: mzML run startTimeStamp '-infinity' is not a date-time; …`) where
  the source prints a non-fatal error for the run attribute and nothing at all
  for the completion time in a Release build.
- **Empty text** is silently unset under the option and an error under the
  default; the source is silent in both cases.
- **Not covered.** Other `xsd:dateTime` CV terms that could appear in a header
  (`MS:1002435` data processing start time) still fail the read under both
  policies, as do non-date CV values of the wrong type, which the source also
  drops with a warning. That is a separate leniency, not part of this lane.

## 2. Size-derived ceilings

### The rule

Every cumulative allowance is now `floor + per_byte * consumed`, credited as the
reader consumes XML bytes, and applies **together with** an absolute ceiling:

```text
charge fails  <=>  charge > min(absolute ceiling, floor + per_byte * consumed)
```

`src/format/mzml_scaling.rs` holds `Allowance`, `InputScaling` and the `Ledger`
that reconciles the reader's plain `usize` counters after every XML event, so
the charge sites are unchanged. The absolute ceilings keep their `ReadOptions`
fields and now default to unbounded; a caller that sets one keeps it exactly.
`InputScaling::fixed()` restores the former fixed ceilings, which is what the
regression test uses to show that they reject a realistic document.

Quantities that need their own XML start tag — records, binary arrays,
parameter groups — are already bounded by the input size (at least ~40 bytes
each), so they keep only the absolute ceiling.

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
| `params` | 10,000,000 | 4 | 0.007 | 570x |
| `param_bytes` | 512 MiB | 256 | 19.4 | 13x |
| `metadata_work` | 50,000,000 | 64 | 0.11 | 600x |
| `metadata_bytes` | 256 MiB | 64 | 0.11 | 600x |
| `selection_work` | 500,000,000 | 512 | not instrumented; ~180 units per selected peak, at most 0.37 peaks per byte | ~8x |
| `selection_bytes` | 256 MiB | 64 | 24 bytes per selected peak | ~7x |
| `count_work` | 50,000,000 | 16 | one unit per XML event | large |

The margins also cover encodings the benchmark inputs do not use. The densest
realistic encoding of a peak is Numpress linear plus short-logged-float with
zlib, about 2.7 base64 characters per peak, which is 0.37 peaks and 6 decoded
bytes per input byte — still 20x inside `peaks` and 10x inside `array_bytes`.

Two ceilings stay absolute because they are per item, not cumulative:

- `max_xml_bytes`, now 1 TiB. It only stops an unbounded stream; a file ends on
  its own, and every other allowance scales with what was consumed.
- `max_array_bytes`, now 512 MiB per array (eight times the former ceiling):
  67 million `f64` values, far beyond any instrument's single spectrum. It also
  bounds the one transient value vector a decoded array allocates, which is why
  it is not unbounded; the cumulative `array_bytes` allowance bounds each array
  as well.

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
- a 2 MB document whose every spectrum expands a 20,000-parameter group,
  charging about 17,000 storage bytes per input byte against an allowance of
  256 (refused after about a megabyte, whatever the document's total size);
- 200 levels of nesting (the reader's fixed 128-level limit, unchanged).

### Measured on the real inputs

`ibminode06`, release build, inputs staged on node-local `/scratch`, one
process per input (`../oracle/mzml-reader-scale/hpc_scale_06.sh`, log
`hpc_scale_ibminode06.log`). "read" is the `ReadOptions::source()` stream read;
the process also counts the same file again with `load_size_with_options`,
which is the rest of the wall time. The node carried a foreign load of about
0.2 per core throughout, so the times are indicative and the peak RSS is not.

| input | size | spectra | peaks | read | process wall | peak RSS |
|---|---|---|---|---|---|---|
| sanity profile + sanity centroid | 5.2 MB + 7.6 MB | 730 + 1,929 | 252,018 | 0.2 s + 0.1 s | 0.6 s | 31 MiB |
| `UK222_picked.mzML` | 547 MB | 40,856 | 22,776,198 | 6.8 s | 15.2 s | 731 MiB |
| `50amol_R1.mzML` (zlib) | 1.20 GB | 43,745 | 88,434,492 | 14.7 s | 26.0 s | 1.85 GiB |
| `20100219_SvNa_SA_Ecoli_preccorrected.mzML` | 1.51 GB | 34,894 | 85,315,432 | 9.4 s | 22.3 s | 1.58 GiB |
| `UK222.mzML` | 2.32 GB | 40,856 | 197,765,338 | 16.1 s | 42.0 s | 3.33 GiB |
| `UK222.mzML`, FeatureFinderCentroided load options | 2.32 GB | 6,911 (MS1) | 96,433,834 | 16.9 s | 17.1 s | 1.89 GiB |

For comparison, the C++ Release `BaselineFilter` peaked at 3.2 GB on
`UK222.mzML` in the smoke benchmark, and `MapNormalizer` at 1.5 GB on
`50amol_R1.mzML`. Every spectrum count matches `inputs/MANIFEST.json`. A first
run of the same tests, before the per-array ceiling was tightened from 2 GiB to
512 MiB, gave the same counts and RSS within 1% (`probe` and `hpc_scale` logs in
the oracle directory).

## Evidence

- **Tier 1, executed C++**: the timestamp table above, produced by
  `../oracle/mzml-reader-scale/datetime_sentinel_cpp.sh` on the Release build.
  Input, per-case inputs and written outputs are hashed in
  `tests/data/mzml_reader_scale_provenance.json`.
- **Tier 4, measured**: the per-byte ratios and the HPC scale runs. Every staged
  benchmark mzML input is read by the `#[ignore]`d tests in
  `tests/mzml_reader_scale.rs`; spectrum counts are compared against
  `inputs/MANIFEST.json`, which was produced by a Python header scan, not by
  this port. Wall time and peak RSS per input are in
  `../oracle/mzml-reader-scale/hpc_scale_ibminode06.log`.
- **Tier 4, Rust-only**: the synthetic-document tests, the adversarial
  documents, the explicit-ceiling precedence and the `Ledger` unit tests in
  `src/format/mzml_scaling.rs`.

## Tool wiring (integrator)

`FileHandler::load_experiment` and `load_experiment_with_options` still pass
`mzml::ReadOptions::default()`, so a TOPP tool run on
`50amol_R1.mzML` still exits 6. The one-line change the tools need is
`ReadOptions::source()` on the tool load paths, alongside the D10 wiring already
in progress for `source_dangling_references`.
