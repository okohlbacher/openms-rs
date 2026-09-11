# Raw MSNumpress codecs

`format::numpress` implements the raw public numerical operation surface in
OpenMS4-core `FORMAT/MSNUMPRESS/MSNumpress.h` at revision
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. It has no additional dependency and is
available without optional XML features. This is separate from the
`MSNumpressCoder` wrapper, base64/zlib compression, mzML CV handling and file I/O.
Those transport operations are not added by this increment.

```rust
use openms::format::numpress::{
    optimal_linear_fixed_point, encode_linear, decode_linear,
};
let positions = [100.0, 200.0, 300.00005, 400.00010];
let fixed_point = optimal_linear_fixed_point(&positions)?;
let bytes = encode_linear(&positions, fixed_point)?;
let decoded = decode_linear(&bytes)?;
assert_eq!(decoded.len(), positions.len());
# Ok::<(), openms::Error>(())
```

## Public API

| Source operation | Rust operation |
| --- | --- |
| `optimalLinearFixedPoint` | `optimal_linear_fixed_point` |
| `optimalLinearFixedPointMass` | `optimal_linear_fixed_point_mass` |
| `optimalSlofFixedPoint` | `optimal_slof_fixed_point` |
| `encodeLinear`, `decodeLinear` | `encode_linear`, `decode_linear` |
| `encodePic`, `decodePic` | `encode_pic`, `decode_pic` |
| `encodeSlof`, `decodeSlof` | `encode_slof`, `decode_slof` |
| `encodeSafe`, `decodeSafe` | `encode_safe`, `decode_safe` |
| Pointer/count and vector overloads | Borrowed slices and owned result vectors |

Every operation also has a `_with_limits` variant taking a final
`&NumpressLimits` argument. All return `Result`. Results are fully computed before
being returned: callers can assign the result after `?` without exposing a
partially changed destination. No additional mutable-output wrapper is needed.

## Source arithmetic and bytes

Linear and SLOF headers contain the fixed-point IEEE binary64 bits in **big-endian**
order. Linear's first two quantized values are stored as little-endian 32-bit
words; SLOF values use little-endian 16-bit words. Safe stores each initial value
or double residual as big-endian binary64. No native-endian pointer reinterpretation
is used.

Linear quantization is the source `trunc(value * fixed_point + 0.5)` to signed
64-bit, including its asymmetric rounding for negative inputs. The first two
values serialize only their low 32 bits. Their decoder reads unsigned 32-bit
words. This defined source truncation is retained: a caller-selected factor can
lose high bits or turn a negative initial integer into a large positive decoded
value. The optimal helper is intended for ordinary positive m/z or RT inputs;
it is not a universal validation guarantee. Subsequent prediction uses
`current + (current - previous)`, with signed 32-bit residuals interpreted using
the SDK's two's-complement representation. Checked arithmetic preserves the
source operation order and rejects signed overflow.

The shared integer code removes leading zero or `f` hexadecimal digits and puts
their count in a leading nibble. Remaining digits are emitted least-significant
first. The first nibble occupies the high half of the next byte. A final unused
low nibble is zero padding. Decoding accepts complete nonminimal encodings, as
the source does. It rejects truncated headers/bodies and incomplete integer
codes. A zero low nibble is padding only at the final half-byte position; this
rule also permits a final zero integer encoded by nibble `8`.

PIC preserves the **compiled source guard**, which requires `value >= -0.5` and
`value + 0.5 <= INT_MAX`. Its header advertises a wider unsigned range than its
implementation accepts. Decoding supports every unsigned 32-bit pattern,
including values above that encoder limit. Negative inputs from `-0.5` through
zero round to zero.

SLOF uses literal `log(value + 1) * fixed_point`, rounds by adding `0.5`, and stores
a 16-bit unsigned integer. It does not substitute `ln_1p`. Its inverse is literal
`exp(word / fixed_point) - 1`, not `exp_m1`. The helper uses `floor(65534 / maximum)`
and initializes the maximum logged value to one. For finite input outside the
positive log domain, the source maximum comparison may ignore a nonfinite log;
the helper preserves that behavior, without promising the input is encodable.
Negative finite factors and small negative scaled values remain supported where
the source conversion truncates into the unsigned range. Negative-zero factors
also retain finite source outcomes: a positive decoded word divided by negative
zero yields negative infinity, whose exponential gives decoded `-1`.

Safe is the source codec name, **not a guarantee of bit-exact roundtrip for every
f64 sequence**. It stores the first two doubles and then double residuals from
linear prediction. Ordinary floating-point subtraction/addition can change the
last bits during reconstruction. The implementation keeps the literal source
operation order.

The linear accuracy helper returns zero for fewer than three inputs without
examining their values or the requested accuracy. Otherwise it returns
`0.5 / accuracy`, or the source sentinel `-1` when this exceeds the maximum
fixed point. It does not reject all finite negative accuracy values. Optimal
helpers return zero for empty data. An otherwise selected nonfinite result,
including the linear optimum for a one- or two-element all-zero array, is a
checked native error.

Empty linear/SLOF encoding is exactly the eight-byte factor header, and an
exactly eight-byte stream decodes to no values. Those unused factor bits are
retained even when nonfinite. Empty PIC/Safe encoding and decoding returns empty
vectors. Empty Safe decoding is a deliberate checked extension: the source
pointer decoder reads past an empty buffer.

## Checked bounds

Consumed numeric inputs, fixed points and returned/stored results must be finite,
except intermediate values that the documented source branches discard. Invalid
floating-to-integer conversions, i64 prediction overflow and out-of-range signed
residuals return errors. The C++ encoder's upper i64 test rounds `LLONG_MAX` to
`2^63` when comparing doubles; Rust explicitly rejects `2^63` before casting.
SLOF rejects odd payload lengths, which would read beyond the source buffer.
Safe rejects lengths not divisible by eight. None of the source pointer or
vector-overload out-of-bounds behavior is reproduced.

Default limits are 10,000,000 input/output values, 128 MiB of encoded bytes and
500,000,000 work units per operation. Limits are independent of mzML reader
limits and can be supplied explicitly. Encoder allocation is conservatively
preflighted at `8 + 5*n` for linear, `5*n` for PIC, `8 + 2*n` for SLOF and `8*n`
for Safe. A byte limit may therefore reject an encoding before its shorter final
compressed length is known. Work is precharged at 64 units per input value for
encoders and 64 per encoded input byte for decoders; helpers charge 16 per value.
Work bounds may bind before the value/byte limits.

Decoders first validate and count the **exact** output without allocating a
result vector, then reserve that count and materialize in a second pass. Both
passes are included in the shared work precharge. This avoids pessimistic
`2 * encoded_bytes` output caps rejecting a legitimate small decoded array.
Length arithmetic and reservations are checked. Input slices remain borrowed;
allocation overhead and process RSS are not promised by the logical bounds.

## Independent and executed evidence

[Tests](../tests/numpress.rs) combine three kinds of evidence:

- The wrapper class test contains three exact base64 strings for the four-value
  source vector. [The literal fixture](../tests/data/numpress_source_bytes.tsv)
  preserves the original text and source line, with raw bytes independently
  decoded using Python's standard library. The raw MSNumpress class test itself
  is only a constructor placeholder and supplies no scientific byte golden.
- Independent hexadecimal-digit/nibble oracles cover leading-zero/one patterns,
  full unsigned decoding, little-/big-endian fields and malformed streams.
  Source large-series length/accuracy assertions, selected numeric boundaries,
  Safe rounding and a long overflowing linear recurrence are checked separately.
- [295 differential cases](../tests/data/numpress_cpp_differential.tsv) were
  actually encoded using the **unmodified pinned C++ source**, linked only to a
  [temporary driver](../../oracle/probes/numpress_probe/probe.cpp). The driver also ran
  the decoders and all three helpers. Eight empty Safe cases skip its undefined
  empty decoder and use the native empty-result extension; all other 287 cases
  execute a C++ decoder. Rust checks every encoded byte and helper result;
  selected nonfinite helper results are expected native errors. Decoded values
  are compared exactly except SLOF, where host libm variation permits `2e-14`
  relative/absolute tolerance. These are generated differential fixtures, not
  upstream class-test literals.

The probe used Apple clang 21.0.0, C++17, `-O0`, on arm64 macOS. It required no
whole-project build, external library, source modification or dependency installation.
The original source's `LLONG_MAX`-to-double warning was retained. Exact command,
compiler/source/driver hashes, seed and fixture hashes are in
[provenance](../tests/data/numpress_provenance.json). The
[regeneration script](../tools/generate_numpress_reference.py) requires an
explicit source checkout and explicitly compiled probe; normal Rust tests run
only the packaged fixtures and never execute C++.

The original Johan Teleman/Lund University BSD-3-Clause notice is retained in the
Rust module. This covers raw numerical APIs, not the higher-level coder's error
fallback policy, estimated-factor defaults, base64/zlib, float-vector adaptation,
Numpress CV terms or mzML transport. Those remain a separate integration step.


The source [MSNumpressCoder wrapper](MSNUMPRESS_CODER_SUPPORT.md) is now available
under the independent `numpress` feature. It supplies base64/zlib transport and
source encoding-policy diagnostics; mzML transport integration remains separate.
