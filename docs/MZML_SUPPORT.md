# mzML support in the Rust port

The `mzml` Cargo feature provides a **bounded mzML 1.1 reader/writer for common peak data**. It uses native Rust XML, base64, and zlib libraries; it does not call OpenMS C++. It is enabled by default and can be omitted with `--no-default-features`.

```rust
use std::fs::File;
use std::io::{BufReader, BufWriter};
use openms::format::mzml::{self, WriteOptions};

fn main() -> openms::Result<()> {
let experiment = mzml::read(BufReader::new(File::open("input.mzML")?))?;
mzml::write_with_options(
    BufWriter::new(File::create("output.mzML")?),
    &experiment,
    &WriteOptions { zlib_compression: true },
)?;
Ok(())
}
```

The API accepts `BufRead`/`Write`; file handling and outer gzip decoding are the caller's responsibility. `read` processes XML events, then returns an owned `MSExperiment` containing all records. It is not a lazy spectrum iterator or an on-disk random-access reader.

## Supported data

| Area | Reader | Writer |
| --- | --- | --- |
| XML | XML 1.0 encoded as UTF-8 or US-ASCII; mzML namespace, including namespace prefixes | UTF-8, plain mzML 1.1.0 |
| Indexed mzML | Sequentially reads the inner mzML; ignores index offsets/checksum | No index generation |
| Peak arrays | Little-endian IEEE f32 or f64 coordinate/intensity arrays, base64 with optional whitespace, uncompressed or zlib | f64 coordinates and f32 intensities; uncompressed by default, optional zlib |
| Named auxiliary arrays | `MS:1000786` names with f32/f64, signed i32/i64, or NUL-terminated ASCII strings; uncompressed or zlib | Native f32, signed integer annotations encoded as i64, and NUL-terminated ASCII strings; explicit `arrayLength` |
| Spectra | Native ID, MS level, profile/centroid flag, scan retention time, m/z and intensity | Same fields; unset retention time `-1` omitted |
| Chromatograms | Native ID, time and intensity arrays, one precursor | Same fields |
| Time units | Explicit seconds or minutes; minutes converted to seconds for scan time and chromatogram coordinates | Seconds |
| Precursors | One selected ion per precursor; selected m/z, charge/intensity, isolation target/offsets, activation methods/energy, possible charges, mobility and spectrum reference | Same supported fields; native fields outside this subset are rejected |
| User parameters | Direct `run`, `spectrum`, and `chromatogram` userParam names/values as strings | Stored maps written as string userParams |
| Container names | Reserved record userParam `openms-rust:name` | Same reserved userParam |

Intensity values are converted to the kernel's `f32` type. A f64 intensity outside the finite f32 range is an error; ordinary f64-to-f32 rounding is expected. All decoded coordinates and intensities must be finite. Peak order is preserved. Missing MS level uses the kernel default of 1; missing spectrum representation remains `Unknown`.

Auxiliary arrays preserve names and order within each native type on both spectra and chromatograms. Float values convert to finite f32; signed integer values must fit native i32, without rounding through floating point. A nonempty auxiliary array must have the record's peak count. Explicit `arrayLength="0"` preserves empty annotation placeholders independently of that count. Names must be nonempty and unique across all three auxiliary types.

The string encoding is `MS:1001479`: every ASCII element, including an empty string, ends in exactly one NUL. Empty elements remain aligned with peaks; unlike the source decoder's empty-token omission, they are never dropped. Embedded NUL and non-ASCII string values are rejected before writing, and missing final terminators, non-ASCII bytes or mismatched element counts are read errors. XML array names and record metadata still support UTF-8. This representation does not label UTF-8 binary string values as ASCII.

Empty records with `defaultArrayLength="0"` can have empty primary arrays; the reader also permits both primary arrays to be absent for an empty record. Nonempty records must have exactly one matching coordinate array and one intensity array. The writer generates `index=N` for an empty spectrum native ID and `chromatogram=N` for an empty chromatogram native ID; generated IDs become populated when read back. Nonempty spectrum IDs must satisfy the mzML `key=value` shape.

## Precursor acquisition fields

Precursors attach directly to `MSSpectrum` and `MSChromatogram`. Isolation target (`MS:1000827`) and lower/upper offsets (`MS:1000828`/`1000829`) use m/z units. A missing selected m/z falls back to the isolation target. When target and selected m/z agree, reading normalizes `isolation_target_mz` to `None`; the numeric target is retained implicitly, without preserving whether the source explicitly wrote it. A different target remains `Some(value)`. Writing rejects a negative effective target whenever an isolation window is emitted, including fallback from a signed selected m/z. Algorithms follow the source's selected-m/z convention when defining purity windows.

All 19 native activation methods are mapped to their PSI-MS accessions. Supplemental CID/HCD terms (`MS:1002679`/`1002678`) map to `Etcid`/`Ethcd`, losing the original alias spelling. Activation energy (`MS:1000509`) is represented in electronvolts. The distinct collision-energy term (`MS:1000045`) and supplemental collision energy are outside this subset and are ignored on read; they are not relabeled as activation energy. Possible charge states (`MS:1000633`) preserve order, duplicates and zero independently of the selected charge.

Selected-ion mobility supports drift time in milliseconds (`MS:1002476`, `UO:0000028`), inverse reduced mobility (`MS:1002815`, `MS:1002814`), FAIMS compensation voltage (`MS:1001581`, `UO:0000218`) and collision cross section (`MS:1002954`, `UO:0000324`). Finite signed values, including negative FAIMS voltages, are preserved. Specified units must match; absent units use the accession's stated native unit. No conversion is inferred.

`spectrumRef` is retained as `Precursor::spectrum_reference`. Reading permits unresolved references for filtered input; parent lookup can then use acquisition-order fallback. Writing requires each reference to identify a spectrum in the output document before any bytes are written. It does not enforce parent MS level or order as an XML constraint. Duplicate isolation/activation containers, duplicate supported scalar fields, conflicting mobility quantities and misplaced containers are rejected.

Drift-window offsets and arbitrary precursor CV terms/ordinary CV-list metadata have no supported writer representation and are rejected. A drift value requires a known unit, and a unit requires a value. Activation containers use an explicit placeholder only when neither a method nor an energy is available. These limitations apply equally to chromatogram precursors.

## Explicit limits and unsupported features

`read_with_options` accepts `ReadOptions`. Defaults are 512 MiB XML input, 64 MiB compressed or decoded binary bytes per array, 10 million total peaks, and 1 million total spectrum/chromatogram records. Additional cumulative limits are 512 MiB decoded array bytes, 20 million array elements and 1 million binary arrays. These include primary arrays; empty string elements count toward the element limit, while empty placeholders count toward the array limit. These bounds cover annotation storage independently of peak count. Parameter definitions and expansion have additional defaults of 100,000 groups, 10 million combined groups/parameters/reference occurrences, and 512 MiB cumulative conservative parameter bytes; these include unused definitions and ignored-header expansion.

XML nesting is limited to 128 levels. Base64 expansion has a corresponding bound. Zlib decoding uses bounded 8 KiB chunks with overrun detection and requires a complete stream with matching checksum and full input consumption; trailing or concatenated compressed data is rejected. Numeric bytes must match the declared count times precision exactly. Variable-size string data is checked against per-array and cumulative byte budgets, then against the exact declared element count. Adjust limits deliberately for larger experiments; the returned experiment still needs memory for all peaks, string objects and retained metadata.

The reader supports [referenceable parameter groups](MZML_PARAM_GROUPS_SUPPORT.md) through the same handlers as inline parameters, including grouped CV and user parameters. Missing references and malformed groups are errors. The reader rejects Numpress, integer/string primary peak arrays, unimplemented semantic auxiliary CV types, other binary CV terms or codecs, other time units, nonfinite values, unsupported XML encodings, DTDs, CDATA, and entity references in element text. Auxiliary units and binary-array userParams have no native storage and are rejected, including unit attributes attached to precision or compression terms. Standard XML escaping and numeric character references in attribute values are supported. Multiple selected ions within one precursor and duplicate userParam keys cannot be represented and are rejected. Conflicting or duplicate supported scientific CV values, duplicate records, missing arrays, incorrect array lengths, and relevant record/array/precursor/selected-ion/scan list counts are also errors.

Validation checks XML well-formedness and the scientific structures used by this subset. It does **not** run XSD or PSI controlled-vocabulary semantic validation while reading. Header inventories, arbitrary CV placement, all schema ordering constraints, and index/checksum validity are not verified. Parse errors use line zero when an exact XML line number is unavailable.

The writer rejects nonfinite data, duplicate IDs, invalid XML characters, reserved metadata keys, ambiguous auxiliary names, misaligned arrays and unencodable strings before producing output. It validates every record and array before the first write; an underlying stream error can still leave a partially written document. Schema-required software, instrument, processing, and activation containers use explicit placeholders for unavailable information. No instrument identity or activation method is invented.

## Metadata preservation boundary

This adapter is **not a lossless mzML archival converter**. Source files, instruments, scan windows, polarity, original data-processing history, arbitrary controlled-vocabulary lists, units on userParams, typed userParam semantics, product ions, acquisition CV fields outside the supported precursor subset, and nested userParams are not retained by this adapter. Multiple selected ions within one precursor remain unsupported. The only automatic name mapping is the reserved Rust name userParam; other title conventions are not inferred.

Stored metadata maps are string-valued. The `openms-rust:name` key is reserved and rejected in user maps; it is used only for the spectrum/chromatogram `name` fields. Unknown acquisition metadata outside the supported model is ignored during reading. Unsupported binary arrays are rejected, rather than discarded.

Use the documented peak/metadata subset for analysis and interchange. Preserve the original source file when its acquisition metadata is needed.

## Verification and provenance

The mzML tests cover independent mixed-precision input, both compression modes, spectra and chromatograms, minute-to-second conversion, precursor fields, UTF-8 and escaped attributes, empty records, malformed XML, inconsistent counts, truncated/corrupt/trailing zlib bytes, nonfinite/overflowing values, resource limits, and writer preflight errors. Independent review added regression tests for conflicting scientific CV fields, the reserved metadata key, and forbidden XML characters.

Independent auxiliary-array tests cover both float and signed integer widths, exact native integer limits, negative zero, empty strings/placeholders, ASCII controls, malformed encoding and cumulative resource boundaries. Writer review checks exact decoded bytes and validates populated and empty-array documents against the pinned XSD in both compression modes. The [chromatogram workflow](../tests/chromatogram_workflow.rs) connects peak picking, exact raw-sample integration and native mzML annotation interchange.

[Precursor workflow tests](../tests/precursor_workflow.rs) verify all 19 activation methods, four mobility quantities, isolation targets and widths, parent references, both compression modes, parser conflicts and output atomicity. Populated acquisition output is checked against the pinned XSD. The original `PrecursorPurity_input.mzML` is also bundled byte-for-byte: it is ASCII with a Latin-1 declaration, which the test changes to UTF-8 only in memory. All 3,872 original MS1 peaks are read, and seven scalar/map cases match independent source-derived numerical references. See [purity provenance](../tests/data/precursor_purity_provenance.json).

Two earlier small fixtures are copied byte-for-byte from the pinned OpenMS source: `MzMLFile_2_minimal.mzML` and `FIAMS_output/SerumTest_picked_10.mzML`. The indexed Serum fixture declares ISO-8859-1 but consists entirely of ASCII bytes. Its test explicitly changes only the declaration to UTF-8 in memory before parsing; the original fixture is unchanged. Its three decoded peaks are asserted against fixed values. The `mzml_independent.mzML` fixture was authored separately with Python's `struct`, `base64`, and `zlib` libraries, using explicit expected values and mixed precisions/codecs; it is not generated by the Rust writer.

The pinned `mzML_1_10.xsd` has identical schema content after newline normalization: the stored fixture uses CRLF and the immutable source uses LF. Both byte hashes are recorded separately. `xmllint --nonet` validated writer output for both a populated compressed experiment and an empty experiment on the development host. The schema test runs when `xmllint` is available and explicitly reports when it is unavailable. This provides independent **structural** validation, not complete PSI semantic or cross-application compatibility certification. No C++ build or C++ runtime comparison was performed.

Copied fixture/schema source paths and SHA-256 hashes are recorded in [mzml_provenance.json](../tests/data/mzml_provenance.json). Source revision: [OpenMS4-core `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`](https://github.com/okohlbacher/OpenMS4-core/tree/7c029e8cdba6abab503708ecdd56f6ab55e38ce4). The implementation was informed by the pinned [MzMLHandler](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp), [mzML schema](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/share/OpenMS/SCHEMAS/mzML_1_10.xsd), and local quick-xml 0.39.4 APIs.
