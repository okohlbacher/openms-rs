# MSNumpressCoder wrapper

`format::numpress_coder::MSNumpressCoder` implements the public wrapper operations
at OpenMS4-core revision `54a232fe2cae9c590d5c997fa49d20e7769860fb`, using the
[raw Numpress codecs](MSNUMPRESS_SUPPORT.md) and the existing
`peak_options::{NumpressConfig, NumpressCompression}` types. This increment does
not change mzML parsing/writing, `WriteOptions`, Numpress CV handling or the
other peak-file adapters.

The independent `numpress` Cargo feature enables the wrapper using the already
available optional `base64` and `flate2` dependencies. `mzml` enables `numpress`;
users can select `--no-default-features --features numpress` without XML. The raw
`format::numpress` numerical codecs remain unconditional. No new dependency or
lockfile change is required.

```rust
use openms::format::numpress_coder::{
    MSNumpressCoder, NumpressConfig, NumpressCompression,
};
let coder = MSNumpressCoder::default();
let config = NumpressConfig {
    compression: NumpressCompression::Linear,
    ..Default::default()
};
let result = coder.encode(&[100.0, 200.0, 300.00005], false, &config)?;
assert!(result.is_encoded());
let decoded = coder.decode(&result.output, false, &config)?;
assert_eq!(decoded.len(), 3);
# Ok::<(), openms::Error>(())
```

## API and outcomes

| Source operation | Native API |
| --- | --- |
| Constructor/destructor | `MSNumpressCoder::default()`, owned destruction |
| Config defaults/string mapping | Existing `NumpressConfig`, `set_compression`, `NumpressCompression::ALL`/`name` |
| `encodeNPRaw(double vector, result, config)` | `encode_raw`; `encode_raw_into` for source destination behavior |
| `encodeNP(double vector, result, zlib, config)` | `encode`; `encode_into` |
| `encodeNP(float vector, result, zlib, config)` | `encode_f32` |
| `decodeNPRaw(bytes, output, config)` | `decode_raw` |
| `decodeNP(base64, output, zlib, config)` | `decode`; atomic native `decode_into` |

All encoders return `Result<NumpressEncodeReport<T>>`, with bytes or text in
`output`, an explicit status, the actual fixed point when relevant, and a boolean
recording use of the maximal-factor fallback. Status is `Encoded`, `EmptyInput`,
`Disabled`, or `Rejected`. Rejection is either a raw codec/estimation/verification
message or a typed accuracy failure carrying the highest failing index and its
original/decoded values. Output is empty unless encoding succeeds. There is no
implicit stderr and **no fallback to ordinary uncompressed IEEE bytes**.

The source raw encoder returns without changing the caller's destination for
empty input, NONE mode or failure. `encode_raw_into` retains that behavior. The
source base64 wrapper clears its destination before encoding, so `encode_into`
clears on a completed skip or rejected encoding. Native resource/transport errors
leave destinations unchanged. All decode results are owned and `decode_into`
commits only success, deliberately replacing the source clear-before-decode
behavior with atomic error handling. `encode_f32` returns a report; it performs
exact f32-to-f64 promotion before entering the f64 encoder, including for NONE.
The source has no float-output decoder; none is invented.

## Factor selection and verification

Defaults come from the existing configuration: fixed point zero, tolerance
`0.0001`, mode NONE, estimation enabled and desired linear mass accuracy `-1`.
Names are exactly `none`, `linear`, `pic`, `slof`, case-sensitive. Safe remains a
raw-only codec because it is not a wrapper configuration mode.

Linear estimation uses the desired-accuracy helper only when
`linear_fp_mass_acc > 0`. A negative returned factor triggers the maximal-factor
helper, with the fallback recorded in the report. Otherwise linear estimation
uses the maximal factor directly. SLOF uses its own helper; PIC ignores all
factor/estimation settings. Explicit factors are used when estimation is disabled.
Unused scalar settings remain ignored, including nonfinite values on branches
where the source does not consume them.

A source corner is preserved: the desired-linear-accuracy helper returns zero
for fewer than three values. The wrapper does not treat zero as a fallback
signal. With verification enabled the resulting zero-factor stream is rejected;
with checking disabled it may be returned successfully, though decoding it
later fails the native finite-value checks. This is a reportable source behavior,
not an automatic switch to a different codec.

Verification runs only when `error_tolerance > 0`. Zero, negative values and NaN
disable it. Positive infinity enables the checks but makes finite-error bounds
unrestrictive. The comparisons scan from the final input back toward the first:

- PIC rejects nonfinite decoded values or absolute error **greater than or equal
  to one**. The magnitude of any positive configured tolerance is ignored.
- Linear/SLOF require finite original and decoded values. If original is zero,
  compare `abs(decoded)` to tolerance; if decoded is zero, compare
  `abs(original)`. Otherwise compare `abs(1 - original / decoded)` to tolerance.
  Rejection uses strict `>`; exact equality passes. This ratio is asymmetric.

The native raw codecs can reject nonfinite or invalid intermediate data before
a verification array exists. Such a failure is reported as a codec rejection
rather than an invented source failing index. As in the source catch-all encode
policy, errors returned from raw allocation/encoding/verification are explicit
rejections. Wrapper resource preflight failures and base64/zlib transport failures
are `Err`. Raw numerical checked boundaries remain those documented by the raw
codec module.

## Transport and checked corrections

Encoding runs Numpress, optionally zlib, then standard padded base64. It adds no
terminating NUL. Decoding reverses that order. Zlib uses the existing default
compression implementation. Compressed bytes may differ from source zlib while
the payload is equivalent; no exact compressed-writer-byte identity is claimed.

Source `Base64::decodeSingleString` ignores text shorter than four **bytes**,
even malformed text or a zlib request. The wrapper preserves this short-input
empty-result behavior. Longer input must have valid standard base64 length,
alphabet, padding and trailing pad bits. Whitespace and URL-safe variants are
rejected. This is an explicit checked correction to the source's permissive
SIMD byte transformation, not a claim of full `Base64` class parity.

Zlib decoding checks progress, CRC/structure and a complete single zlib stream.
Trailing bytes and concatenated members are rejected; the source iterative
inflater stops at the first stream end and can ignore trailing input. Gzip/raw
Deflate are not silently substituted for zlib. Transport decoding still happens
for NONE; only the subsequent raw decoder skips data. Raw empty input itself
returns empty without asking the raw codec to parse a header.

## Shared limits

`NumpressCoder { limits: NumpressCoderLimits { ... } }` supplies independent
resource settings. Defaults reuse the raw caps of 10,000,000 values, 128 MiB of
binary data and 500,000,000 work units, plus 192 MiB of base64 text and 512 MiB of
cumulative logical allocation accounting.

One shared counter covers promotion, estimation, optional fallback, encoding,
verification, error comparison, transport and output copies. Raw operations are
precharged at the same declared costs as the raw module; their local bounds do
not reset the wrapper counter. Verification reserves its known input count.
Standalone wrapper decoding conservatively charges the source-style decoded
upper bound, clamped to the value limit, before calling the raw decoder's exact
counting pass. This allocation accounting can reject a small decoded result
before its actual payload reaches the cap. Binary limits count every compressed
or decompressed buffer, not only base64 input size.

Base64 decoding preflights the exact padding-adjusted byte length. Zlib processes
bounded 16 KiB input/output chunks, charges every call before execution, and
checks each output extension. Output capacity grows geometrically; cumulative
new capacities and prior-buffer copy work are precharged before reservation.
A conservative 1 MiB logical allowance is charged before codec state creation.
This is logical payload/state accounting, not a process-RSS or allocator-overhead
guarantee. Decompression expansion and late errors expose no partial output.

## Evidence

[Tests](../tests/numpress_coder.rs) preserve all three source exact base64 strings
and the source 100-value length/accuracy checks. Raw source bytes are reused by
hash from the preceding raw-codec fixture. The
[transport fixture](../tests/data/numpress_coder_transport.tsv) was independently
compressed using Python zlib; it is a derived projection, not an upstream literal
or a native-writer-generated expectation. Its
[regeneration script](../tools/generate_numpress_coder_reference.py) uses only
Python's standard library.

The remaining tests cover fallback and zero-factor behavior, typed highest-index
rejections, exact tolerance boundaries, NONE/empty semantics, f32 promotion,
source short-base64 behavior, strict malformed input, multi-chunk zlib, CRC and
trailing data, expansion caps, cumulative raw/verification/transport allocations,
and atomic destinations. An independent read-only review checked source
arithmetic, resource precharges and transport boundaries. No C++ wrapper was
compiled or executed; prior raw C++ differential evidence applies only to the
separate raw-codec stage. [Provenance](../tests/data/numpress_coder_provenance.json)
records source, reused fixture and new projection/script hashes.
