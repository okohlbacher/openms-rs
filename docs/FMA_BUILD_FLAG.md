# The x86_64 FMA build flag, and the processor check that goes with it

`.cargo/config.toml` builds x86_64 with `-C target-feature=+fma`. This document
says why, what it costs, what the startup check does, and what was measured.

Short version for a user: **an x86_64 binary built from this checkout needs an
FMA3-capable processor — Intel Haswell (2013), AMD Piledriver (2012) or newer.**
On anything older it prints what it needs and exits 12 instead of crashing. To
build for an older processor:

```sh
RUSTFLAGS="-C target-feature=-fma" cargo build --release --locked
```

The results are the same either way; only the speed differs.

## 1. Why the flag

The port reproduces the reference build's glibc `exp`, `log` and `powf` bit for
bit (lead decisions D5 and D10 of wave 5), which means fusing multiply-adds
exactly where glibc does. In Rust that is `f64::mul_add`, which is a fused
multiply-add by contract and never a multiply followed by an add. On a baseline
x86-64 target the processor has no such instruction, so every `mul_add` becomes
an out-of-line call that dispatches again through a CPUID-set function pointer:
**41 such call sites in the hot math**, 9 in `exp`, 23 in `log` and 9 in `powf`.

[BENCHMARKS](BENCHMARKS.md) section 4 measured what that costs, on dax, run
`fma-ffc-dax-1` in
`/ceph/ibmi/abi/oliver/bench/openms4/results/2026-09-17-fma/`:

| pair | 1 thread | 32 threads |
|---|---:|---:|
| `rust-ffap` / `rust-main` — the cost of the bit-exact math | 1.214 | 0.917 |
| `rust-ffap-fma` / `rust-ffap` — what the flag removes | **0.715** | 0.922 |
| `rust-ffap-fma` / C++ Release | **1.055** | **0.894** |
| `rust-ffap` / C++ Release | 1.475 | 0.969 |

Output was `bitwise_equal` with and without the flag at both thread counts.
That run left the decision open; the user has since decided to enable it by
default, and this document records what that entails.

## 2. What the flag is scoped to

```toml
[target.'cfg(target_arch = "x86_64")']
rustflags = ["-C", "target-feature=+fma"]
```

Cargo reads `.cargo/config.toml` from the working directory of the `cargo`
invocation upwards, **not** from a package's own directory. So the flag applies
to builds run inside this checkout — the gates, CI, the benchmark harness — and
**not** to a project that depends on `openms` by path from somewhere else. Such
a consumer gets a baseline build, which is correct but slower, and its binaries
have no processor requirement and no check.

`aarch64` must not see this flag: fused multiply-add is in its baseline and the
target does not accept the feature name. The `cfg` scope is what keeps it away.

**Measured** (`cargo check --lib --no-default-features --target T -v`, counting
rustc invocations, macOS host, 2026-09-18):

| target | units compiled | units with `target-feature=+fma` |
|---|---:|---:|
| `x86_64-apple-darwin` | 24 | 24 |
| `x86_64-pc-windows-msvc` | 23 | 23 |
| `aarch64-apple-darwin` | 24 | **0** |

and on kim, a native `x86_64-unknown-linux-gnu` release build of all binaries,
**73 of 73 rustc invocations** carry the flag.

## 3. What the flag actually turns on

`fma` is not a feature a processor can have on its own: FMA3 instructions are
VEX-encoded, so rustc implies `avx`, `sse3`, `ssse3`, `sse4.1` and `sse4.2` from
it. Measured on kim with `cargo rustc --lib -- --print cfg`:

| build | `target_feature` |
|---|---|
| default | `avx fma fxsr sse sse2 sse3 sse4.1 sse4.2 ssse3` |
| `RUSTFLAGS="-C target-feature=-fma"` | `fxsr sse sse2` |

So the whole crate is compiled with VEX encoding and 256-bit auto-vectorisation,
not only the `mul_add` sites. In the release binaries on kim:

| binary | `vfmadd`/`vfmsub` | VEX (v-prefixed) mnemonics |
|---|---:|---:|
| `FeatureFinderCentroided`, default | 41 | 41,075 |
| `FeatureFinderCentroided`, opt-out | 2 | 442 |
| `FileInfo`, default | 0 | 23,738 |
| `FileInfo`, opt-out | 0 | 440 |

Two things to read out of that table. The 41 fused multiply-adds are exactly the
41 call sites the flag was for. And `FileInfo`, which never touches the ported
glibc math, still gains 23,738 VEX instructions — **the processor requirement is
not confined to the tools that do fused arithmetic.** The 2 `vfmadd` in the
opt-out binary are inside the standard library's own runtime-dispatched `fma`,
which is what a baseline build calls instead.

## 4. The opt-out, and why that one command is enough

`RUSTFLAGS` **replaces** the `rustflags` of `.cargo/config.toml` rather than
adding to them: Cargo takes its extra flags from the first of four sources that
is set, and the environment variable comes before the config file. `+fma` is the
only entry in that config, so one variable removes the whole requirement.

Verified on kim: `RUSTFLAGS="-C target-feature=-fma" cargo build --release
--locked --bins` succeeds, the resulting crate is compiled with the baseline
feature set (section 3), the binaries run, and their output is byte for byte the
output of the default build (section 7).

Anyone adding a second entry to that config must revisit the message in
`src/system/cpu_features.rs`, which quotes this command — and will be made to:
`tests/fma_build_flag.rs` reads `.cargo/config.toml` and fails if it holds any
table but the x86_64 one or any `rustflags` line but this one, because a second
entry would be silently dropped by the very command the crate prints.

## 5. The startup check

`src/system/cpu_features.rs` decides, and `src/cli.rs` acts:

* `build_requires_fma()` is `cfg!(all(any(target_arch = "x86", target_arch =
  "x86_64"), target_feature = "fma"))` — a compile-time fact about this crate;
* `cpu_provides_fma()` executes `cpuid` and reads the architectural FMA bit;
* `unsupported_cpu()` returns the message when the first is true and the second
  is false, and `None` otherwise;
* `cli::run` calls it as its **first statement** and, on `Some`, writes the
  message to standard error and returns `ExitCode::InternalError` (12). The
  table in [TOPP_CLI_SUPPORT](TOPP_CLI_SUPPORT.md) carries that row.

No `unsafe` anywhere in this crate; `#![forbid(unsafe_code)]` is unchanged.

The message names the requirement, the processor generations that satisfy it and
the rebuild command verbatim, and it is a constant, so printing it neither
allocates nor computes. It goes out through `Write::write_all` rather than the
formatting machinery, and `report_unsupported_cpu` is `#[cold]` and
`#[inline(never)]`, because this code runs on a processor that cannot execute
everything the build emitted and the less of the build it uses the better.

## 6. The standard library cannot answer this

The obvious implementation is `std::arch::is_x86_feature_detected!("fma")`. It
does not work, and the first version of this change shipped a guard that was not
in the binary at all.

`is_x86_feature_detected!` is documented to answer `true` **without asking the
processor** whenever the feature is already enabled at compile time. That is
exactly this case. Measured on kim (rustc 1.96.0), with a two-line probe whose
only body is the macro, compiled `-O` and marked `#[inline(never)]`:

| build | `cfg!(target_feature="fma")` | macro value | the probe function in the binary |
|---|---|---|---|
| baseline | `false` | `true` | present: loads the `std_detect` cache, calls the initialiser if empty |
| `-C target-feature=+fma` | `true` | `true` | **eliminated — folded to a constant** |

The same thing happened in the tool binaries. With the macro, `FileInfo::main`
disassembled to a bare `jmp` into the body of `run`, with no check in it:

```
<FileInfo::main>:
  jmp    <openms::cli::run_from_environment>
```

So the question goes to `raw-cpuid` (`=11.6.0`, MIT, x86 targets only, no
dependency of its own beyond `bitflags`, which the graph already had), which
always executes `cpuid`. Reading CPUID by hand means `core::arch::__cpuid`,
which is `unsafe`; the crate keeps that out of this one. The unit test
`the_detector_agrees_with_what_the_kernel_reports` compares its answer with the
`flags` line of `/proc/cpuinfo` on Linux, so the detector is checked against
something outside the build.

## 7. The flag changes no result

Measured on kim, release binaries built with and without the flag from the same
commit, on `tests/data/topp_feature_finder_centroided/FileConverter_31_output.mzML`:

* `FileInfo -in … -test`: both exit 0, and the reports are byte for byte equal;
* `FeatureFinderCentroided -in … -out … -test`: both exit 0, and the featureXML
  files are byte for byte equal (sha256 `cf3ac065a08814c0…`).

The whole test suite was run with the flag on; section 9 records the counts
against the recorded main run.

## 8. What the check does and does not guarantee

**It is exact for a processor with AVX but without FMA** — Sandy Bridge, Ivy
Bridge, Bulldozer. There, every instruction in the binary except the 41 fused
multiply-adds is legal, the check runs, the message is printed and the process
exits 12.

**It is best effort on a processor without AVX either** — Nehalem, Westmere,
Core 2 and older, all pre-2011. `fma` implies `avx`, so such a processor can
fault on a VEX instruction that is not a fused multiply-add, before the check is
reached. Two things were done about that, and neither is a guarantee:

1. the check is the first statement of `cli::run`, and everything after it is in
   `run_from_environment`, which is `#[inline(never)]`, so no instruction of the
   body can be hoisted above the check. This was not cosmetic: with the check
   one level lower, in `run_with`, `cli::run`'s own prologue emitted `vpxor`,
   `vmovdqu` and `vmovups %ymm0` while building the argument vector, all of them
   before the check;
2. the message path avoids the formatting machinery, as section 5 describes.

What remains on the path, measured in the release `FileInfo` on kim, is this
whole prologue — `cli::run` inlined into the tool's `main`, and not one VEX or
SSE instruction before the check:

```
<FileInfo::main>:
  push   %rax
  call   *<got>                      # unsupported_cpu
  test   %rax,%rax
  jne    <the message path>
  call   <openms::cli::run_from_environment>
  pop    %rcx
  ret
```

`main` itself, ahead of it, is six integer instructions and the call into
`std::rt::lang_start`, which is standard-library code compiled for the baseline.
`report_unsupported_cpu` is likewise integer-only: it locks standard error,
makes two `write_all` calls, drops the lock and returns 12. A different compiler
version may emit different code; this is a measurement, not a property of the
design.

**A build machine without FMA is a separate problem.** In a native build, Cargo
applies `target.<cfg>.rustflags` to host artifacts too: on kim all 13
build-script units were compiled with `+fma`. A build script that does
floating-point work would then need FMA on the *build* machine, before any of
this runs. The same `RUSTFLAGS` opt-out covers that case, because it replaces
the flags for every unit.

**C code is not affected.** `flate2`, `bzip2` and `sha1` are pure Rust in this
configuration, and `rusqlite`'s and `libxml`'s C is compiled by `cc` with its own
flags, which `RUSTFLAGS` does not reach.

**Unit tests are not covered.** The check is in the CLI entry point; a test
binary does not go through it, so on a processor without FMA `cargo test` would
still fault. That is acceptable — a machine that cannot run the tools cannot run
their tests either — but it is the reason CI's processor matters (section 10).

## 9. Evidence

Everything in this document was measured on the gate host kim (AMD EPYC 9654,
rustc 1.96.0) except the cross-target scoping in section 2, measured on the
macOS arm64 host, and the timings in section 1, which are from the recorded dax
run. The scripts and their output are under `/scratch/kohlbach/fma-default*` on
kim; the lane's scratchpad note `lane-notes/fma-default.md` names them.

**No emulator was available.** None of the five IBMI gate hosts — kim, dax,
spock, ibminode05, ibminode06 — has `qemu-user`, Intel SDE or an environment
module system that could provide one, so an FMA binary could not be run under a
`-cpu Nehalem` or `-cpu Westmere` guest. The end-to-end behaviour was instead
recorded on a **scratch copy of the checkout whose `cpu_provides_fma` was edited
to return `false`**, built and run on kim. Nothing in the repository was
weakened: the patch lives only in `/scratch/kohlbach/fma-default3/patched` and
its anchor is asserted, so it fails loudly rather than silently not applying.
That copy, with everything else identical:

| command | exit | stdout | stderr |
|---|---:|---|---|
| `FileInfo -in … -test` | **12** | 0 bytes | the message, verbatim |
| `FeatureFinderCentroided -in … -out … -test` | **12** | — | the message; **no output file was written** |
| `FeatureFinderCentroided --help` | **12** | — | the message |

`--help` refusing is deliberate, not an oversight: the check precedes every
source phase, so a processor that cannot run the tool is told that rather than
shown usage it cannot act on.

What that copy does *not* cover is the CPUID read itself, which is `raw-cpuid`'s
and is cross-checked against `/proc/cpuinfo` by
`the_detector_agrees_with_what_the_kernel_reports`. The decision between the two
answers is covered on any host by
`only_a_requiring_build_on_a_processor_without_fma_is_refused`, which injects
both.

Gates: see the lane record. The full `--all-features --all-targets` run with the
flag is compared against main `e1c3115`'s recorded 5,297 passed / 0 failed / 21
ignored.

## 10. Continuous integration

Every job now builds with the flag, so every x64 runner must have FMA.

* `macos-latest` and the other macOS runners are arm64 (GitHub documents
  `macos-latest` as "3 (M1), arm64"), so the flag is not applied there at all —
  the `cfg` scope, measured in section 2, is what makes that true.
* `ubuntu-latest` and `windows-latest` are x64. GitHub does not publish the
  processor model. Third-party measurement of the current fleet reports AMD EPYC
  7763 (Zen 3) and AMD EPYC 9V74 (Zen 4) for the x64 runners, both of which have
  FMA3; the ARM64 runners are Neoverse-N2. This is the best confirmation
  available without running a job, and it is a claim about a fleet that can
  change. If a runner ever lacks FMA, the tool binaries say so and exit 12 while
  the unit tests fault (section 8).

## 11. The benchmark harness

The harness builds with `cargo build --release --locked --offline`, inside the
checkout, so it now picks the flag up from `.cargo/config.toml` without being
told. Two consequences:

* a `rust-*` cell built from a commit at or after this change is an FMA build,
  and is **not** comparable with the `rust-main` and `rust-ffap` cells of
  BENCHMARKS section 4, which were baseline builds. The comparable cells there
  are `rust-ffap-fma` and `rust-main-fma`;
* a cell that deliberately wants a baseline build must now say so, with
  `RUSTFLAGS="-C target-feature=-fma"`, where before it got one by saying
  nothing.
