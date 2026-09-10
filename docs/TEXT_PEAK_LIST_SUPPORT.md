# MS2 and DTA2D peak lists

`format::ms2` and `format::dta2d` provide buffered `read`, `read_with_options`,
atomic `read_into`/`read_into_with_options`, `load`/`load_with_options`, and
`write`/`write_with_options`/`store`/`store_with_options`. The path wrappers delegate
to the same adapters. They do not add loaded-file path/type metadata: the native
experiment does not have OpenMS's separate DocumentIdentifier state.

The pinned [MS2File header](https://github.com/okohlbacher/OpenMS4-core/blob/6bfc0e4711105f4eda2fea86812a83af7c7e791f/src/openms/include/OpenMS/FORMAT/MS2File.h)
is an **input-only** adapter. Its Rust writer is an explicit native extension.
[DTA2DFile](https://github.com/okohlbacher/OpenMS4-core/blob/6bfc0e4711105f4eda2fea86812a83af7c7e791f/src/openms/include/OpenMS/FORMAT/DTA2DFile.h)
has source load, store and TIC-store operations, all represented here. No C++
execution, native dependency or external converter is involved.

## MS2 conventions

The source reader trims ASCII space/tab/CR/LF and collapses that whitespace
between fields. A line starting `S` must contain four fields; the fourth is the
precursor m/z. The scan-number fields are not interpreted. Each scan gets MS
level 2, one precursor with unknown charge/zero intensity, RT -1 and an
`index=N` native ID. Empty scans survive, and peak order is unchanged.

Lines starting `H`, `I`, `Z` or `D` are ignored regardless of their remaining
text, exactly as in the pinned adapter. In particular, **Z charge and I RTime
are not loaded**. Valid peak lines before the first S record are parsed and
counted against limits but discarded; a file with no S records has no spectra.
Other nonempty lines require exactly two finite numeric fields.

The native writer emits S records and peaks only. It requires MS2, unset RT and
exactly one precursor without charge, intensity or acquisition metadata. It
rejects populated names, spectrum type, annotations, identifications and
experiment metadata/chromatograms before producing output. Empty and ordinary
`index=digits` IDs are accepted as source-generated bookkeeping and regenerated;
other IDs are rejected. These restrictions make the stored scientific fields
readable by the source-compatible reader without silent loss.

## DTA2D conventions and filters

Data rows have exactly three fields: RT in seconds, m/z and intensity by default.
If a line contains a tab, the source uses tab delimiters; otherwise it uses single
spaces. Repeated delimiters create empty fields and are not silently collapsed.
Trailing whitespace is trimmed. Every `#` line is a header, not a free comment.
The first three header fields select any order of `SEC`/`MIN`, `MZ`, `INT`,
case-insensitively. Deprecated `RT`, `RETENTION_TIME`, `MASS-TO-CHARGE`, `IT` and
`INTENSITY` remain accepted. Additional header fields are ignored. A short or
invalid header returns a checked parse error.

Minutes multiply RT by 60. The source's minutes flag is sticky: after a MIN
header, a later SEC header does not reset it. A new spectrum begins only when
`abs(row_rt - spectrum_anchor_rt) > 0.0001`; the anchor is the first RT in the
group, not the preceding row. The initial anchor is the source's -1 sentinel.
Rows within its tolerance retain RT -1 and an empty native ID until a real RT
change. Later IDs count all raw groups, including those removed by filtering.
Spectra are MS1, unsorted input order is preserved, and empty filtered groups
are omitted.

`dta2d::ReadOptions` exposes `rt_range`, `mz_range`, `intensity_range` as optional
Rust ranges, plus `limits`. All three ranges include their start and exclude
their end, following OpenMS DRange. Intensity filtering uses the parsed f32
value widened to f64. Empty ranges select nothing; reversed/nonfinite endpoints
are errors. These are the only PeakFileOptions fields consumed by the C++
adapter. Generic MS-level, metadata-only, compression and precision flags have
no DTA2D behavior and are not presented as working filters.

Full storage emits `#SEC\tMZ\tINT` and round-trippable decimal rows. It requires
nonempty MS1 spectra without precursors and the same metadata restrictions as
MS2. Adjacent spectra within the grouping tolerance, or a first RT that would
be changed by the -1 sentinel, are rejected. Empty experiments are valid and
produce the header. Output uses Rust's shortest round-trip float spelling;
scientific values are preserved, rather than C++ formatter byte identity.

## TIC projection

`write_tic`, `write_tic_with_options`, `store_tic` and `store_tic_with_options`
explicitly project MS1 spectra into RT/zero-m/z/intensity rows. They use f32
intensity accumulation in peak order, include zero-intensity points for empty
MS1 spectra, and preserve scan order and duplicate RTs. They do not resample or
use stored chromatograms. Higher MS levels and unrelated metadata are
unconsumed, as in `calculateTIC()` with source defaults. Duplicate-RT rows can
be grouped if the resulting file is read as DTA2D spectra.

## Validation and limits

Finite signed coordinates/intensities are supported. Intensity tokens parse
directly to f32, avoiding an extra f64 rounding step. Nonfinite values, nonzero
values that underflow to zero, minute-conversion overflow and nonfinite TIC sums
are checked errors. Parsing returns owned data only on success; replacement
leaves the destination unchanged on parse, stream or limit errors.

Both adapters expose the same `Limits` type. Defaults and hard ceilings are
256 MiB input, 1 MiB per line, 20 million lines, 100,000 raw spectra/groups,
10 million raw peaks and 512 MiB output. Limits can be lowered. Ignored headers,
pre-scan peaks and filtered data still count. A bounded line reader uses at most
one extra byte to detect an exceeded input/line limit; vectors have checked
allocation. Peak/scan limits bound retained storage independently of input size.

Writers validate the entire represented payload and count formatted output
before writing any bytes or creating/truncating a path. All writers flush and
report flush errors. An actual destination I/O failure can still leave partial
output; this is not a transactional filesystem writer.

[Tests](../tests/text_peak_lists.rs) cover four exact source fixtures, source
filter results, scientific rounding, malformed inputs, atomic replacement,
limits, metadata-loss guards and flush failures. The
[provenance manifest](../tests/data/text_peak_lists_provenance.json) records
14 source hashes, four fixture hashes and derived versus literal expectations.
The first fixture TIC is exactly f32 `141649.671875`; the source class test's
short printed expectation is `141650`.
