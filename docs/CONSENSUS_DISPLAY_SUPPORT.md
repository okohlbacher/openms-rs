# Consensus feature display support

Native coverage of the last two public members of `KERNEL/ConsensusFeature.h` at
Core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`: the stream operator and the
ratio description.

| Artifact | Path |
| --- | --- |
| Implementation | `src/kernel/consensus_display.rs` |
| Tests | `tests/consensus_display.rs` (9 tests) |
| Manifest | `tests/data/consensus_display_provenance.json` |

The feature-identification package recorded its residual verbatim as *"All
members covered except two: `Ratio::description_` has no Rust field and
`operator<<` has no Display"*. `Ratio::description` was added to
`src/kernel/features.rs` afterwards; this package supplies the operations over
it and the `Display` implementation, and `KERNEL/ConsensusFeature.h` has no
member left unported.

## API mapping

`docs/FEATURE_IDENTIFICATION_SUPPORT.md` holds the mapping for the other 30-odd
members of this header. The two rows below are the ones it records as not
ported, restated with what they now map to, together with every member of the
nested `Ratio` struct, which this package is the first to cover completely.

| C++ member | Rust |
| --- | --- |
| `operator<<(std::ostream&, const ConsensusFeature&)` | `impl Display for ConsensusFeature` in `src/kernel/consensus_display.rs` — **new here** |
| `struct ConsensusFeature::Ratio` | `kernel::features::Ratio` |
| `Ratio::Ratio()`, copy constructor, `operator=`, `virtual ~Ratio()` | derived `Default`, `Clone`; Rust needs no assignment operator and has no vtable to keep |
| `double Ratio::ratio_value_` | public field `Ratio::ratio_value`, checked by `Ratio::validate` |
| `std::string Ratio::denominator_ref_` | public field `Ratio::denominator_ref` |
| `std::string Ratio::numerator_ref_` | public field `Ratio::numerator_ref` |
| `std::vector<std::string> Ratio::description_` | public field `Ratio::description`, with `Ratio::description`, `Ratio::add_description`, `Ratio::set_description` and `Ratio::validate_description` as the checked path over it — **new here** |
| `// TODO ratio cv info` | not ported: a source `TODO`, not a member |

The header's `@TODO` on `Ratio` — *"members are public, names shouldn't end in
underscores"* — is satisfied by the port, whose fields carry the same names
without the trailing underscore.

## Preserved source conventions

* **The stream layout is reproduced exactly**, including the details that look
  like typos and are not: `"Intensity "` and `"Quality "` carry no colon while
  every other label does; `"Grouped features: "` and `"Meta information: "` end
  in a space before the newline, because the space is inside the source's string
  literal; and the closing banner line ends `"----------------- "` with a
  trailing space. `tests/consensus_display.rs` pins the whole block as one
  string literal, so none of this can be tidied away by accident.
* **Handle order is the source's `HandleSetType` order** — `FeatureHandle::IndexLess`,
  i.e. map index then unique ID — because the port keeps the handles as a sorted
  unique slice. Printing follows storage, as the source's `begin()`/`end()` loop
  does.
* **Ratios are not printed.** The source's operator ignores `ratios_`; so does
  this one.
* **A description is a list, not a set.** `add_description` appends a repeated
  line again, as `push_back` would, and accepts an empty line, because the
  source places no constraint on the strings.
* **Nothing validates a `Ratio` in the source**, and `Ratio::validate` in the
  port deliberately checks only the ratio value, so `ConsensusFeature::add_ratio`
  and `set_ratios` do not measure a description assigned straight to the public
  field. `validate_description` is the explicit check for untrusted input, and
  `tests/consensus_display.rs` asserts the split rather than leaving it implied.
* **No parallelism to reproduce:** `KERNEL/ConsensusFeature.cpp` contains no
  `#pragma omp`.

## Native differences

* **Number formatting.** The source wraps intensity, quality and every handle
  coordinate in `precisionWrapper`, which sets the stream precision to
  `writtenDigits<T>` — `std::numeric_limits<float>::digits10` = 6 for the `f32`
  intensity and quality, `std::numeric_limits<double>::digits10` = 15 for the
  `f64` coordinates — and the position goes through `DPosition`'s own operator
  (`DPosition.h:412-420`), which is `precisionWrapper` per coordinate,
  space-separated. This `Display` writes Rust's shortest round-trip form
  instead, as every other kernel `Display` in the port does: `1.5` prints as
  `1.5`, not `1.50000000000000`. The layout is identical; only the digits
  differ, and a reader comparing against C++ output must expect that.
* **Meta value order.** The source's `getKeys` yields names in
  `MetaInfoRegistry` index order — the order in which the running process first
  registered each name, which is not reproducible from a document. The port's
  `MetaInfo` is a `BTreeMap`, so names print sorted. The set of pairs is the
  same; only their order differs.
* **Two implementations of one layout, for now.** `ConsensusMap`'s `Display` in
  `src/kernel/map_operations.rs` writes the same per-feature block through a
  private `write_consensus_feature` helper, because that module could not call a
  `Display` that did not exist. This package may not edit `map_operations.rs`, so
  the duplicate stays; `display_agrees_with_the_consensus_map_block` asserts the
  two agree character for character, which makes any future drift a test
  failure rather than a silent inconsistency. Collapsing the helper into
  `{feature}` is left to whichever package next owns that file.
* **The description ceilings are native.** The source vector is unbounded.

## Checked boundaries and evidence

| Boundary | Value | Where |
| --- | --- | --- |
| `Ratio::MAX_DESCRIPTION_LINES` | 65 536 lines | `add_description`, `set_description`, `validate_description` |
| `Ratio::MAX_DESCRIPTION_BYTES` | 1 048 576 bytes (1 MiB) of UTF-8 | the same three |
| description byte-sum overflow | `checked_add` | `description_bytes` |

Both ceilings are per-ratio and are checked before anything is stored, so a
refused `add_description` or `set_description` leaves the description
byte-identical; `ratio_description_ceilings_are_checked_and_atomic` asserts that
for both, and asserts that a value exactly at each ceiling is accepted. The byte
preflight re-sums the stored lines on every append, because a `Ratio` keeps no
cached total; the sum is bounded by `MAX_DESCRIPTION_LINES` additions and
allocates nothing. `Display` allocates nothing beyond the formatter's own
buffer and cannot fail except as the sink fails.

**Evidence tiers.** `ConsensusFeature_test.cpp` has no `START_SECTION` for
either member, so the layout expectations are tier 3 (source review) read off
the implementation at `ConsensusFeature.cpp:391-418` rather than off a class-test
literal, and the `precisionWrapper` claims are read off
`CONCEPT/PrecisionWrapper.h`, `CONCEPT/Types.h:178-195` and
`DATASTRUCTURES/DPosition.h:412-420`. One expectation is tier 3 by
cross-reference instead: the printed block asserted in
`display_reproduces_the_source_stream_layout` is the same string the earlier
`cm_display` test in `tests/map_operations.rs` already asserted for
`ConsensusMap`, so two independently written tests now pin one layout. The
description ceilings, their atomicity and the `validate`/`validate_description`
split are tier 4 (Rust-only invariants). No C++ was built or executed, and no
retained C++ output exists for this header.
`tests/data/consensus_display_provenance.json` hashes the five source files and
pins the source anchors.

## Class-test section accounting

`ConsensusFeature_test.cpp` has 39 `START_SECTION`s. All 39 are already ported in
`tests/feature_identification.rs`, one Rust test each, as
`docs/FEATURE_IDENTIFICATION_SUPPORT.md` records in its own section table (92
sections across `BaseFeature_test.cpp`, `Feature_test.cpp` and
`ConsensusFeature_test.cpp`). **None of them exercises `operator<<`,
`Ratio::description_`, `addRatio`, `setRatios` or `getRatios`**: the two members
this package ports have no upstream section at all, which is why their evidence
is read off the implementation.

Self-audit (`ConsensusFeature.h`, this package's two members): 0 sections newly
ported because 0 exist, 0 mapped-with-evidence, 0 mapped-without-evidence, 0
unaccounted; 39 of 39 sections accounted for, all in the earlier package.

## Deferrals

* `src/kernel/map_operations.rs` keeps its private copy of the layout; see
  **Native differences**. Not a porting gap — the same output, twice.
* `tests/data/consensus_display_provenance.json` is not yet listed in
  `SOURCE_PROVENANCE.json`, so `tools/check_core_sdk.py` does not verify its
  hashes. Registering it is the integrator's step.

## Validation

```sh
cargo fmt --all -- --check
cargo clippy --locked --all-features --lib --test consensus_display -- -D warnings
cargo nextest run --locked --all-features --test consensus_display
cargo nextest run --locked --no-default-features --test consensus_display
cargo test --locked --all-features --doc
RUSTDOCFLAGS="-D warnings" cargo doc --locked --all-features --no-deps
cargo +1.85.0 check --locked --all-features --lib --test consensus_display
python3 tools/check_doc_coverage.py --report | grep consensus_display
```

CI: `tests/consensus_display.rs` needs no features and belongs on the
`--no-default-features` kernel line of the `minimum-rust` job
(`.github/workflows/rust.yml:79`), appended as `--test consensus_display`.
