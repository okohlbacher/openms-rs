# SIMDe portability shim

`SYSTEM/SIMDe.h` at SDK `bc9cc12` is 43 lines, has no `.cpp`, no class test and
no Rust counterpart. This document records why, so the header is accounted for
rather than silently unmapped.

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

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `#include <simde/x86/ssse3.h>` | — | Not ported: a third-party C header. |
| `inline simde__m128i operator\|(const simde__m128i&, const simde__m128i&)` | — | Not ported: MSVC-only papering over a compiler gap that does not exist in Rust. |
| `inline simde__m128i& operator\|=(simde__m128i&, const simde__m128i&)` | — | Not ported: as above. |
| `inline simde__m128i operator&(const simde__m128i, const simde__m128i&)` | — | Not ported: as above. |
| `namespace OpenMS {}` | — | Empty in the source. |

No public member of this header has a counterpart, and none is missing one.

## Why there is no counterpart

Rust's portable SIMD story is different in kind. `std::arch` intrinsics are
`unsafe` and this crate forbids `unsafe` crate-wide; `std::simd` is nightly-only
and the crate's MSRV is 1.85. The port therefore contains no hand-written SIMD
at all, so there is nothing for a portability shim to make portable and no
operator gap to fill — Rust has no per-compiler operator availability.

What the header exists to *enable* is visible in one place in this port:
`system::build_info::active_simd_extensions` reports the SIMD level in effect.
The source derives that from the SIMDe arch macros; the port derives it from
Rust target features, keeping the source's order and labels. See
`docs/BUILD_INFO_SUPPORT.md`.

A future Rust SIMD path would not resurrect this header. It would be
`#[cfg(target_feature = …)]` on ordinary safe code, or a dependency, and either
way the decision belongs to whoever takes it, not to a shim ported ahead of a
need.

## Evidence

None is required and none is claimed: there is no behaviour to reproduce. The
header carries **0** class-test sections — no `SIMDe_test.cpp` exists in the
SDK — so no section is left unmapped by recording it this way.

Ledger status: `native_equivalent`, on the grounds that the port achieves the
header's purpose — a consistent, build-system-backed SIMD level — without the
mechanism, because the mechanism addresses a C++ problem that Rust does not
have.
