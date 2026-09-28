# SIMDe portability header

`SYSTEM/SIMDe.h` (core bc9cc12) declares no OpenMS API. Its whole content is:

- `#include <simde/x86/ssse3.h>`, which brings in the vendored SIMDe library
  (`src/openms/extern/simde`). SIMDe maps the SSE/SSSE3 intrinsics to native
  instructions where the target has them and to portable code elsewhere.
- Three free operators on `simde__m128i` — `|`, `|=` and `&` — defined only
  under `_MSC_VER`, because MSVC does not provide them for its vector type.
- An empty `namespace OpenMS {}`.

The header therefore has no member to port. What it provides is a build-time
facility for three implementation files, and each of those has a Rust
counterpart that does not need it:

| Source consumer at bc9cc12 | What the SIMD code does | Rust counterpart |
|---|---|---|
| `FORMAT/Base64.cpp` | 128-bit masks to split and join Base64 sextets when encoding and decoding binary arrays | the `base64` crate (`=0.22.1`, behind the `numpress` feature that `mzml` enables), used by the mzML, mzXML, mzData, imzML and numpress codecs; the decoded bytes are fixed by the format, not by the instruction set |
| `DATASTRUCTURES/StringUtils.cpp` (lines 450-498) | find the first or last byte that is or is not one of space, tab, newline, carriage return, 16 bytes at a time | byte-wise scans in the text parsers; the answer is the same position for any input |
| `FORMAT/HANDLERS/StringManager.cpp` (lines 147-207) | Xerces UTF-16 helpers: find the terminating zero, test whether eight UTF-16 units are ASCII, and pack them into bytes | none needed: there is no Xerces UTF-16 layer; `quick-xml` hands the parsers UTF-8 text, which is already the byte string these helpers build |

The build report's SIMD list — `BuildInfo::getActiveSIMDExtensions`, which in
the source lists the SIMDe probes that were defined — is ported as
[`system::build_info::active_simd_extensions`](../src/system/build_info.rs),
reading the Rust target features in the source's order and with its labels;
see [tests/build_info.rs](../tests/build_info.rs).

Explicit SIMD would be an optimisation and would need a determinism argument
(results must be bit-identical to the serial scalar path), so none is written.
If a profile ever shows the whitespace scans or the Base64 path to be hot, the
candidates are `std::arch` behind a runtime feature check (which needs
`unsafe`, forbidden in this crate) or a safe portable-SIMD crate such as
`wide`; either would get a `THIRD_PARTY_CRATE_DECISIONS.md` entry first.
