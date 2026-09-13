# SIMDe portability shim

`SYSTEM/SIMDe.h` at SDK `bc9cc12` (sha256 `f3e1180f0bfb78b34a2d88559a6b7f4d0ac1d7759405cac7ae60aba9679c7aab`,
pinned by [system_process_provenance.json](../tests/data/system_process_provenance.json))
is 43 lines, has no `.cpp`, no class test and **no Rust counterpart**. Nothing
in this port implements it or stands in for it. This document records that, and
why, so the header is accounted for rather than silently unmapped.

## Ledger status

**No review entry.** `docs/core-sdk-reviewed-apis.json` records a review only as
`complete`, `partial` or `native_equivalent`, and each of those has to cite the
Rust files and the tests that carry it (`tools/core_sdk_coverage.py` asserts
that every cited path exists and that `tests` is non-empty). This package
produced no Rust file and no test for this header, because there is nothing here
to implement, so it claims none of the three.

`docs/core-sdk-coverage.json` is generated. With no review entry, and with this
package's manifest citing the header, the regenerated row reads
`status: evidence_requires_review` with a null `review`, no
`candidate_rust_files`, and `tests/data/system_process_provenance.json` as its
only `reference_manifests` entry. That is the mechanical meaning of the bucket —
evidence exists, no accepted API review stands against it — and it is the state
this document describes. Read the bucket name literally and it promises more
than is coming: the review has been done, and its result is that no API review
is possible, because there is no API. This page is that result.

In particular the status is **not** `native_equivalent`. That status says
Rust's own types and standard library already provide the header's operations;
see [NATIVE_EQUIVALENTS.md](NATIVE_EQUIVALENTS.md). Nothing in Rust provides
`simde__m128i` operators or a single include point for a C SIMD level, because
Rust has neither problem. "The problem does not exist here" is not the same
claim as "the standard library already solves it", and only the first one is
true.

## What the header is

It is a C++ *include discipline*, not an API. Its own comment states the two
jobs:

1. include `<simde/x86/ssse3.h>` from exactly one place, so every translation
   unit sees the same SIMD level the build system enabled (`-mssse3` and
   friends, configured in `cmake/compiler_flags.cmake`);
2. define `operator|`, `operator|=` and `operator&` for `simde__m128i` under
   `_MSC_VER`, because GCC and clang provide them for their vector types and
   MSVC does not.

The header closes with an empty `namespace OpenMS {}`. It declares no type, no
function and no constant of its own, and its own comment forbids including it
from another header.

## Declarations, and what became of each

Every line of the header that declares or includes anything, and the
counterpart it has. There are five, and **the counterpart column is empty in
every row** — this is a not-ported inventory, not an API mapping, and it is
listed in full so that no declaration is left unaccounted for.

| C++ declaration | Rust counterpart | Why there is none |
|---|---|---|
| `#include <simde/x86/ssse3.h>` | none | A third-party C header. Nothing here includes it, directly or transitively. |
| `inline simde__m128i operator\|(const simde__m128i&, const simde__m128i&)` | none | MSVC-only, and it papers over a compiler gap Rust does not have. |
| `inline simde__m128i& operator\|=(simde__m128i&, const simde__m128i&)` | none | As above. |
| `inline simde__m128i operator&(const simde__m128i, const simde__m128i&)` | none | As above. |
| `namespace OpenMS {}` | none | Empty in the source; declares nothing. |

## Why there is no counterpart

Rust's portable SIMD story is different in kind. `std::arch` intrinsics are
`unsafe` and this crate forbids `unsafe` crate-wide; `std::simd` is nightly-only
and the crate's MSRV is 1.85. The port therefore contains no hand-written SIMD
at all, so there is nothing for a portability shim to make portable, and no
operator gap to fill — Rust has no per-compiler operator availability.

A future Rust SIMD path would not resurrect this header. It would be
`#[cfg(target_feature = …)]` on ordinary safe code, or a dependency, and either
way the decision belongs to whoever takes it, not to a shim ported ahead of a
need.

## Evidence

None is required and none is claimed: there is no behaviour to reproduce. What
this package did produce, and all that is cited here, is in
[system_process_provenance.json](../tests/data/system_process_provenance.json):
the header's pinned sha256, the source anchor recording what the header does and
that it declares no OpenMS member, and the `native_boundaries` note stating that
it has no counterpart. The header carries **0** class-test sections — no
`SIMDe_test.cpp` exists in the SDK — so recording it this way leaves no section
unmapped.

Related but **not** evidence for this decision:
`system::build_info::active_simd_extensions` reports which SIMD level is in
effect, deriving it from Rust target features where the source derives it from
the SIMDe arch macros. That function, its tests and
[BUILD_INFO_SUPPORT.md](BUILD_INFO_SUPPORT.md) belong to the `SYSTEM/BuildInfo.h`
port, which was delivered separately and is unchanged by this package; it
reports a SIMD level rather than providing SIMD portability, and it is named
here only so a reader looking for "where did the SIMD level go" finds it.
