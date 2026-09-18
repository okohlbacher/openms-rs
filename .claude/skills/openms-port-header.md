---
name: openms-port-header
description: 'Port one OpenMS C++ header/implementation pair to Rust in the OpenMS4-R crate, or backport its Doxygen documentation onto Rust code that already exists. Use when porting any OpenMS class, algorithm or file format to Rust, when writing or reviewing docs for a ported module, or when asked to "port <Header>.h", "document <module>.rs", "work through the kernel", or continue the TOPP tool port. Covers the eight-artifact deliverable, the Doxygen-to-rustdoc mapping, evidence tiers and the differential test against retained C++ output.'
---

# Porting an OpenMS header to Rust

The crate is `OpenMS4-R`: a from-scratch Rust port of OpenMS 4, no C++ bindings,
no C++ build dependency. The goal is full replacement of `openms4-core`,
`openms4-cli` and `openms4-topp`.

Two jobs share this skill because they share the same source reading:

- **Port** a header/impl pair that has no Rust yet.
- **Backport documentation** onto Rust that already exists. 41% of public items
  were documented when this skill was written; the C++ documents nearly
  everything. Undocumented Rust is unfinished Rust.

Read `AGENTS.md` and `docs/DIFFERENTIAL_VALIDATION.md` in the repo first — they
are short and they own the evidence rules.

---

## Phase 1 — Read the source before writing anything

Find the pinned checkout from `openms::CORE_SDK_REVISION`; it is
`.reference/openms4-core-<short>/`. Never read the upstream working tree, which
moves.

Read, in this order:

1. The header. Every Doxygen block is a specification — `@note` and `@exception`
   especially, because they record behaviour you cannot infer from the code.
2. The `.cpp`. The header states intent; the implementation states what actually
   happens, including the bugs. **60 of the 794 headers have no `.cpp`** — they
   are header-only, usually templates. For those the header *is* the
   implementation; read every inline body, not just the declarations.
3. `src/tests/class_tests/openms/source/<Class>_test.cpp`. **This is the oracle.**
   Its literals are the expected values, and they are usually the only ones you
   will get without running C++.
4. `docs/core-sdk-coverage.json` for this header's status and its
   `direct_topp_consumers`, which tells you who needs it and how much surface is
   actually reachable.

While reading, write down every behaviour that is surprising. Those become
`@note`-equivalents in the Rust docs and, when they are defects, entries in
`OpenMS_CPP_ISSUES.md`.

---

## Phase 2 — Map the API

Rust idiom wins over transliteration; source *behaviour* wins over Rust
aesthetics. Where they conflict, keep the behaviour and document the difference.

| C++ | Rust |
|---|---|
| `getX()` / `setX(v)` | public field, or `x()` / `set_x(v)` when it validates |
| Method returning `bool` success | `Result<()>` |
| Method throwing `Exception::X` | `Result<T>` with the matching `Error` variant |
| Out-parameter `void f(T& out)` | return the value |
| `Size` / `UInt` index | `usize`, or `u32` when it is a stored record field |
| `DataValue` / `ParamValue` | `MetaValue` / `ParamValue` — never stringify |
| Overload set | distinct names — the real ones are `find_nearest`, `find_nearest_with_tolerance`, `find_nearest_in_window` |
| Global *mutable* singleton | caller-owned struct passed in |
| Immutable shared table | `OnceLock` global plus a caller-owned override; see `ModificationsDB::global` |
| `std::vector<T>&` accessor | `&[T]` / `&mut Vec<T>` |
| Iterator pair `begin()/end()` | `impl Iterator`, borrowed |
| Constructor / `operator=` / copy ctor | `new`, `From`, `Default`, `Clone` |
| `operator[]`, `operator==`, `operator<` | `Index`, `PartialEq`, `PartialOrd` |
| Static member function | associated function |
| `template<class T> class` | generic struct with trait bounds |
| Template specialisation | a separate `impl`, or a distinct type |
| `Real` / `DoubleReal` | `f32` for stored intensities, `f64` for coordinates |
| `std::map` / `std::set` | `BTreeMap` / `BTreeSet`, for deterministic order |
| `OPENMS_DLLAPI` | nothing; `pub` carries visibility |
| Conditional compilation | `#[cfg(feature = "...")]` |

Hard rules that are not negotiable:

- **No `unsafe`.** The crate forbids it.
- **No silent lossy behaviour.** If the source discards information, the native
  default refuses and an explicit option selects the source behaviour — see
  `dta::WriteOptions::source()` for the established pattern.
- **Bounded work.** Any operation whose cost scales with input takes a
  `MAX_ITEMS`/`MAX_BYTES`-style ceiling, checked in a `preflight` before
  anything is allocated or mutated, so a failure leaves the input unchanged.
  Copy `src/identification/run_mapping.rs`.
- **Atomicity.** An error must not leave a half-updated value. Build into a
  temporary, then commit.
- **Serial by default.** 36 source files carry `#pragma omp`. The port is
  deliberately serial — `docs/REPOSITORY_ANALYSIS.md` requires "serial reference
  behaviour first; optional Rust parallel iterators after deterministic tests".
  Do **not** introduce threads to match OpenMP. Record in the support doc that
  the source parallelises and the port does not, so the performance gap is
  stated rather than discovered.

**Splitting.** One C++ header does not have to become one Rust module. 194
headers are templated and some exceed a thousand lines. Split by type family or
by functional layer when a module would otherwise sprawl — `src/format/mzml*`
is ~20 files for one header — and say in the support doc which Rust files cover
the header, because the ledger records that mapping.

---

## Phase 3 — Backport the documentation

This is the half that gets skipped. Do not skip it. The target is *the same
level of detail and correctness as the C++*, not a one-line summary.

### Tag mapping

| Doxygen | rustdoc |
|---|---|
| `@brief` | first line of the doc comment, one sentence, no trailing period needed |
| body text | following paragraphs |
| `@param[in] x` | `# Arguments` when any parameter carries a constraint, unit or default; otherwise fold into prose. **Every `@param`'s semantic content must survive somewhere** — only the type restatement may be dropped. There are 5,789 of these; a "non-negative tolerance" or "in seconds" lost here is lost for good |
| `@param[out]` / `@param[in,out]` | describe in the return or the `&mut` receiver's prose |
| `@return` | prose in the first or last paragraph; `# Returns` only when non-obvious |
| `@throws` / `@throw` / `@exception` | `# Errors` section naming the condition and the `Error` variant. The singular `@throw` occurs 237 times — grep for both |
| `@note` | plain paragraph, or `# Notes`. **Never drop these** — they carry preconditions |
| `@warning` | `# Warning` section |
| `@pre` | state it in `# Errors` (checked) or as an explicit precondition sentence |
| `@see X` | intra-doc link ``[`X`]`` **only if `X` exists in the crate**; otherwise plain ``` `X` ``` and note the gap in the support doc. CI runs `cargo doc` with `-D warnings`, so an unresolved link fails the build, and many of the 179 `@see` targets are not ported yet |
| `@ref Foo "text"` | `[text](...)` to the support doc, or drop if it names a Doxygen group |
| `@ingroup` | drop; module structure carries it |
| `@name` group | `///` section headings or module grouping |
| `@c X` / `@p X` | `` `X` `` |
| `@b X` / `@em X` | `**X**` / `*X*` |
| `@code`…`@endcode` | **C++ example:** rewrite to the Rust API as a ```` ```rust ```` doctest and make it compile. **Any other language:** several blocks are Python for pyOpenMS (`MSSpectrum::rasterizeIMFrame`); render as ```` ```text ```` or drop. Never mark something `rust` that cannot compile — CI runs `cargo test --doc` |
| `@deprecated` | `#[deprecated(note = "...")]` when the Rust API keeps it, otherwise do not port the item and say so in the support doc |
| `@author` | drop; authorship lives in `AUTHORS` |
| `@f[`…`@f]` and inline `@f$`…`@f$` | keep the formula as text or a fenced block; do not invent MathJax |
| `@a` / `@li` / `@verbatim` / `@sa` | inline code / list item / fenced text / see-also. Any unlisted `@`-tag: carry its content, drop the markup |
| `@tparam` | generic parameter prose |
| `@todo` | drop from rustdoc; if still true, note it in the support doc |

### Rules

- **Translate, don't transcribe, but never silently drop.** Every `@note` and
  `@exception` must be *accounted for* in one of three ways: carried across,
  restated as the Rust equivalent, or explicitly neutralised — "the source warns
  that `Size` overflows here; `usize` arithmetic is checked, so this cannot
  occur". Silence is the one option that is not allowed, because a reader cannot
  tell it from an oversight.
- **Every divergence from source behaviour is documented at the item**, with one
  sentence saying what the source does and why this differs. These sentences are
  the port's real value.
- **Cite the source when the reason is non-obvious**: `` /// Source
  `MSSpectrum::findNearest` returns -1 when empty; this returns `None`. ``
- **Never claim tested behaviour that is not tested.**
- Module-level `//!` names the source header(s) and points at the
  `docs/*_SUPPORT.md`.
- Prefer a short sentence at the field over a long block at the struct.

### Documentation is not enforced by the compiler

`missing_docs` is not enabled — it cannot be, at 42% coverage. Run the ratchet
instead, which fails when any module loses documentation:

```bash
python3 tools/check_doc_coverage.py            # check
python3 tools/check_doc_coverage.py --report   # per-module, worst first
python3 tools/check_doc_coverage.py --write    # record the new floor after improving
```

**A module you touch leaves at 100%.** Pick your next target with `--report`.

### Worked example

This is the *target* state, not what the crate currently says. The existing Rust
one-liner is exactly the deficit this skill exists to close.

```cpp
/**
  @brief Binary search for the peak nearest to a specific m/z

  @param[in] mz The searched for mass-to-charge ratio searched
  @return Returns the index of the peak.

  @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.

  @exception Exception::Precondition is thrown if the spectrum is empty (not only in debug mode)
*/
Size findNearest(CoordinateType mz) const;
```

```rust
/// Binary search for the peak nearest to a specific m/z.
///
/// Midpoint ties choose the lower m/z; exact duplicates yield the first. That
/// tie rule is native detail, not a source statement, and is documented because
/// callers can observe it.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `mz` is not finite, and
/// [`Error::UnsortedData`] when the peaks are not sorted by m/z. The source
/// `@note` documents unsorted input as undefined and does not check it; this
/// port checks, because the cost is one pass and the alternative is a silently
/// wrong index.
///
/// An empty spectrum yields `Ok(None)` rather than the source's
/// `Exception::Precondition`, so absence is not an error path.
pub fn find_nearest(&self, mz: f64) -> Result<Option<usize>> {
```

The standard, item by item: the `@brief` becomes the summary line; the `@note`
is carried across *and* marked as now-checked; the `@exception` is mapped to its
Rust outcome; the extra finiteness check that the source lacks is stated; and
native-only detail is labelled as such so a reader never mistakes it for source
behaviour.

---

## Phase 4 — Evidence

Tiers, from `docs/DIFFERENTIAL_VALIDATION.md`:

1. **Executed differential** — an oracle driver linking the prebuilt C++, or a
   retained C++ output fixture. The only tier that catches a misread algorithm.
2. **Executed probe** — an extracted translation unit compiled and run.
3. **Source review** — transcribed class-test literals.
4. **Independently derived / Rust-only** — computed from the spec, or native
   invariants and resource bounds.

**Reach for tier 1 whenever a retained C++ output exists.** Concretely:

1. `grep -n "<Tool>" OpenMS4-tests/packages/test-data/topp/CMakeLists.txt` — each
   `add_test` names the exact invocation, and the following `${DIFF}` line names
   the produced file and the retained expected output.
2. Copy input and expected output into `tests/data/` under the group's prefix.
3. Reproduce the invocation through `run_with::<T>(...)` in an integration test,
   into a per-case temp directory. Cases that produce the same file name must
   not share one — they run in parallel.
4. Compare per the table in `docs/DIFFERENTIAL_VALIDATION.md` — read it, do not
   guess the tolerances.
5. Record input and expected-output hashes in the provenance manifest, with
   `evidence_tier` naming tier 1 and the upstream `CMakeLists.txt` line range.

`tests/topp_dta_extractor.rs` (byte-exact, text format) and
`tests/topp_mzml_splitter.rs` (canonical, mzML) are the two worked shapes.

**Never duplicate the tool body into its test.** A `[[bin]]` target cannot be
imported by an integration test, so the obvious move is to copy it — and the
copy drifts: the first two tools did exactly this and their registered parameter
descriptions had already diverged from the shipped binaries. Put the `Tool` impl
in a library module (`src/cli/tools/<tool>.rs`), keep `src/bin/<Tool>.rs` a
three-line `main`, and have the test import the library module. One definition,
tested and shipped.

Byte equality is **not** the default contract: the upstream suite compares with
FuzzyDiff. Use it only for plain-text formats where the source writes exact
bytes. For anything else compare decoded content.

**Do not try to run the prebuilt TOPP executables.** All 130 currently reject
their own tool manifest; `docs/DIFFERENTIAL_VALIDATION.md` records the state.
Tier 1 for a tool means the *retained* outputs above, not executing C++.

The upstream `${DIFF}` is FuzzyDiff with a whitelist —
`CMakeLists.txt:360` exempts `id=`, `href=`, `completion_time=` and `version=`.
For an XML-producing tool those four are exactly what you must not compare.

For a core algorithm with no retained output, tier 1 needs an oracle driver: a
small `main` in `../oracle/drivers/<group>.cpp` linking the prebuilt
`../core-build` library, reading a JSON case file and printing full-precision
results. `../core-build/lib/libOpenMS.dylib` links today and 710 class-test
binaries run.

**No C++ is committed to this repository.** Probe sources and oracle drivers
live in `../oracle/`; manifests record them under `external_reference_artifacts`
with their sha256, and `tools/check_core_sdk.py` enforces that shape.

---

## Phase 5 — The eight artifacts

A group is not done until all of these exist. This is the repo's actual
convention, measured across its commits:

1. Rust module(s) under `src/<domain>/`
2. `tests/<group>.rs` integration binary
3. `docs/<GROUP>_SUPPORT.md` — four required *contents*, not four literal
   headings: **API mapping** table / **preserved source conventions** / **native
   differences** / **checked boundaries and evidence**. The API mapping table lists *every* public source member and its
   Rust counterpart, or records it as not ported. A member missing from that
   table is the usual way an unported API goes unnoticed.
4. `tests/data/<group>_provenance.json`. Copy the shape of a neighbour rather
   than inventing keys — core groups use `source_revision`, `sources[{path,sha256}]`,
   `fixtures`, `native_implementation`, `native_tests`, `support_document` and
   `target_verification` (`tests/data/datetime_provenance.json`); TOPP groups use
   `package_revisions`, `method`, `evidence_tier` and `upstream_test_definition`
   (`tests/data/topp_cli_provenance.json`), because the cli/topp packages are
   versioned separately from the core SDK and must **not** be listed in
   `current_sdk_reference_manifests`
5. Ledger entry in `docs/core-sdk-reviewed-apis.json` with `status`, `rust`,
   `tests`, `documentation`, `scope`; then `python3 tools/core_sdk_coverage.py --write`
6. `docs/VALIDATION.md` checkpoint recording the executed checks
7. `README.md` row and `CHANGELOG.md` line
8. CI wiring, and it differs by kind of group:
   - **Library group:** add `--test <name>` to the `minimum-rust` job, on the
     line whose feature selection matches yours.
   - **TOPP tool:** add a `[[bin]]` with `required-features` to `Cargo.toml`,
     *and* append `--test topp_<name>` to the feature-sliced line in the **`test`**
     job — `rust.yml:25`, not `minimum-rust`. Both existing tool tests live
     there because they need `mzml paramxml`; wiring a tool into `minimum-rust`
     silently tests nothing.

Status in the ledger is earned, not asserted: `complete` needs every public API
reviewed *and* documented; `partial` means started; do not promote on hope. The
coverage script does not and cannot check this — it is the one honest judgement
the skill asks for, and `evidence_requires_review` is the correct default when
in doubt.

Every module you touch also leaves at **100% rustdoc coverage**, recorded by
`python3 tools/check_doc_coverage.py --write`. This is checkable, so it is
checked.

---

## Phase 6 — Verify

All of these, every time:

```bash
cargo fmt --all -- --check
cargo clippy --locked --all-features --all-targets -- -D warnings
cargo test --locked --all-features --all-targets
cargo test --locked --all-features --doc          # CI runs this; doctests must compile
cargo test --locked --no-default-features         # the portable core must still build
RUSTDOCFLAGS="-D warnings" cargo doc --locked --all-features --no-deps
for t in check_core_sdk test_core_sdk core_sdk_coverage test_core_sdk_coverage \
         check_schema_feature_graph check_doc_coverage; do
  python3 tools/$t.py && echo "$t OK" || echo "$t FAIL"
done
```

If the group touches chemistry, enzymes, modifications, RNA, monosaccharides or
any controlled vocabulary, it almost certainly regenerates an embedded table.
CI's quality job runs ten `--check` generators plus
`semantic_validator/projection.py` and `probes/ims_witness_source_oracle.py`;
mirror that job rather than guessing which apply:

```bash
grep -n 'python3 tools/' .github/workflows/rust.yml
```

Then run **your own test under the feature set CI will use for it**, and add
that line to the `minimum-rust` job. `.github/workflows/rust.yml` pins Rust
1.85 and runs ~20 distinct feature selections; a module that compiles only with
`--all-features` passes locally and fails there. `check_schema_feature_graph.py`
separately proves the default and no-default builds pull no C-dependent crate.

When a test fails, find the root cause before touching the expectation. If the
expectation really was wrong — because the port was *stricter than the source*
and the source's own fixtures violate it — change it and write the evidence into
the test as a comment. That has happened twice and both were real findings.

---

## Failure modes seen in practice

- **Documentation skipped** because the code compiled. The most common one.
- **`@note` dropped.** They encode preconditions the code does not check.
- **Stricter than source.** The port rejected an mzML whose list `count`
  disagreed with the child count; the upstream fixture does exactly that and C++
  loads it, so the port could not read its own reference data.
- **A lossy guard blocking parity.** The DTA writer refused to discard metadata
  that `DTAFile::store` silently drops. Resolution: guard stays the library
  default, explicit source option, tool opts in.
- **Byte comparison on a format that never matched byte-for-byte.** Check what
  the upstream test uses before asserting equality.
- **Ledger promoted without review.** `complete` means every public API was
  read; `evidence_requires_review` is the honest default.
- **Reading the upstream working tree** instead of the pinned checkout.
- **A `@code` block that is not C++.** Some are Python for pyOpenMS. Marking one
  ```` ```rust ```` breaks `cargo test --doc`.
- **Threading to match OpenMP.** The port is serial on purpose; state the gap
  instead.
- **Claiming tier 1 from a source-review fixture.** Tier 1 needs output the C++
  actually produced, hashed in the manifest.
- **A runtime CPU check the build flag folds away.**
  `is_x86_feature_detected!("fma" | "avx" | "sse3" | "ssse3" | "sse4.1" |
  "sse4.2")` is a compile-time `true` on x86_64 since `.cargo/config.toml` sets
  `-C target-feature=+fma`, and `-O` then deletes the guard, with no diagnostic
  and no clippy lint. Ask `system::cpu_features::cpu_provides_fma()` when the
  question is about the **processor**; `is_x86_feature_detected!` and
  `cfg!(target_feature = ...)` answer a question about the **build**.
