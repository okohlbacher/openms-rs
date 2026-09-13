# Isotope source precision and bounding-box geometry support

Package B2-ISO-GEOM of the early TOPP bundle adds three source contracts that
`FeatureFinderAlgorithmPicked` needs and the native port did not provide:

1. an explicit **source-precision** mode in which coarse isotope patterns are
   computed with the C++ `Peak1D` binary32 arithmetic
   (`ProbabilityPrecision::SourceSingle`), leaving the `f64` default unchanged;
2. the **source `trimLeft`**, which keeps every peak when none reaches the
   cutoff (`IsotopeDistribution::trim_left_source`);
3. `BoundingBox2D::intersects`, `width` and `height` from `DBoundingBox.h`.

Source revision: OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.
Rust files: `src/chemistry/isotopes.rs` (both isotope headers) and
`src/kernel/geometry.rs` (`DBoundingBox.h` for two dimensions). The native
isotope conventions stay documented in [ISOTOPE_SUPPORT.md](ISOTOPE_SUPPORT.md);
this document covers the source-precision contract, the source `trimLeft` and
the bounding-box predicates. Provenance and hashes:
`tests/data/isotopes_source_precision_provenance.json`.

## Why the precision matters

`FeatureFinderAlgorithmPicked.cpp:364-425` precalculates one averagine pattern
per 100 Da mass window with `CoarseIsotopePatternGenerator(max_isotopes)`, trims
it at `intensity_percentage_optional`, classifies optional peaks against
`intensity_percentage` and scales by the maximum. All of that runs on `float`
intensities. The native `f64` pattern differs from the source by about 1e-7
relative, which is enough to flip a threshold decision or a seed score in the
last bit. The source-precision mode reproduces the executed C++ SDK bit for bit
(tier 1, below); B6 still has to close the pattern comparison for every
FeatureFinderCentroided_1 window against the C2 oracle.

## API mapping

### `CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h`

| Source member | Rust | Notes |
|---|---|---|
| `CoarseIsotopePatternGenerator(Size max_isotope = 0, bool round_masses = false)` | `CoarseIsotopePatternGenerator::new(Option<usize>, CoarseMassMode)`, `Default` | `0` is `None`; `Some(0)` and more than `MAX_ISOTOPE_PEAKS` are rejected. Precision starts at `Double`. |
| `~CoarseIsotopePatternGenerator()` | implicit drop | |
| `setMaxIsotope` | not ported | construct a new generator; the limit is fixed at construction |
| `setRoundMasses` | not ported | construct with the wanted `CoarseMassMode` |
| `getMaxIsotope` | `max_peaks` | |
| `getRoundMasses` | `mass_mode` | `false` = `Approximate`, `true` = `Nominal` |
| `setIsotopeOverride(const Element*, const IsotopeDistribution&)` | `set_isotope_override(&str, IsotopeDistribution)` | symbol key; validated (see native differences); narrowed to `f32` when read in source precision |
| `clearIsotopeOverrides` | `clear_isotope_overrides` | |
| `getIsotopeOverrides` | not ported | no accessor; keys are private atom values |
| `run(const EmpiricalFormula&)` | `run` | honours `ProbabilityPrecision` |
| `estimateFromPeptideWeight` | `estimate_from_peptide_weight` | |
| `estimateFromPeptideMonoWeight` | `estimate_from_peptide_mono_weight` | |
| `estimateFromPeptideWeightAndS` | `estimate_from_peptide_weight_and_sulfur` | |
| `approximateFromPeptideWeight` (static) | `approximate_from_peptide_weight` (associated) | `f64` recurrence; not affected by precision |
| `approximateIntensities` (static) | `approximate_intensities` (associated) | `f64`, as the source |
| `estimateFromRNAWeight` | `estimate_from_rna_weight` | |
| `estimateFromRNAMonoWeight` | `estimate_from_rna_mono_weight` | |
| `estimateFromDNAWeight` | `estimate_from_dna_weight` | |
| `estimateFromWeightAndComp(w, C, H, N, O, S, P)` | `estimate_from_weight_and_comp(w, AveragineComposition)` | |
| `estimateFromMonoWeightAndComp` | `estimate_from_mono_weight_and_comp` | |
| `estimateFromWeightAndCompAndS(w, S, C, H, N, O, P)` | not a generator method | formula via `AveragineComposition::estimate_average_mass_with_sulfur`, then `run` |
| `estimateForFragmentFromPeptideWeight` | `estimate_fragment_from_weights(.., PEPTIDE)` | |
| `estimateForFragmentFromPeptideWeightAndS` | not ported | fixed-sulfur fragment wrapper, deferred in ISOTOPE_SUPPORT.md |
| `estimateForFragmentFromRNAWeight` | `estimate_fragment_from_weights(.., RNA)` | |
| `estimateForFragmentFromDNAWeight` | `estimate_fragment_from_weights(.., DNA)` | |
| `estimateForFragmentFromWeightAndComp` | `estimate_fragment_from_weights(.., composition)` | |
| `calcFragmentIsotopeDist` | `calc_fragment_isotope_dist` | honours precision |
| `operator=` | `Clone` | |
| `convolve` (public) | `convolve` | honours precision |
| `convolvePow_` (protected) | `convolve_power` (public) | honours precision; see native differences |
| `convolveSquare_` (protected) | private (a pattern convolved with itself) | truncation equivalence below |
| `correctMass_` (protected) | private `correct_masses` | |
| `calcFragmentIsotopeDist_` (protected) | inside `calc_fragment_isotope_dist` | |
| `fillGaps_` (protected) | private `dense_from_distribution`, `single_from_distribution` | |
| `max_isotope_`, `round_masses_`, `isotope_overrides_` | private `max_peaks`, `mass_mode`, `overrides` | |
| (new) | `ProbabilityPrecision`, `with_precision`, `set_precision`, `precision` | explicit source binary32 option |

### `CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h`

| Source member | Rust | Notes |
|---|---|---|
| `MassAbundance` (`Peak1D`) | `IsotopePeak` | both fields `f64`; source precision stores binary32 values |
| `ContainerType`, iterator typedefs | `peaks()` slice | |
| `enum Sorted` | not ported | declared, unused by the implementation |
| `IsotopeDistribution()` | `Default` (`(0, 1)`), `empty()` | |
| copy and move constructors, destructor, `operator=` | `Clone`, moves | |
| `set(const ContainerType&)`, `set(ContainerType&&)` | `from_peaks` | validated |
| `getContainer` | `peaks` | |
| `getMax`, `getMin` | `max_mass`, `min_mass` | `None` when empty (source `0`) |
| `getMostAbundant` | `most_abundant` | `None` when empty (source `Peak1D(0, 1)`) |
| `size`, `clear`, `resize` | `len`, `clear`, `resize` | |
| `trimIntensities` | `trim_intensities` | |
| `sortByIntensity`, `sortByMass` | `sort_by_probability`, `sort_by_mass` | stable sorts |
| `renormalize` | `renormalize`, `renormalize_with(ProbabilityPrecision)` | |
| `merge` | not ported | deferred in ISOTOPE_SUPPORT.md |
| `trimRight` | `trim_right` | identical, including removal of every peak |
| `trimLeft` | `trim_left_source` | source; `trim_left` keeps the native all-below removal |
| `averageMass` | `average_mass` | |
| `operator==`, `operator!=` | `PartialEq` | |
| `operator<` | not ported | |
| `begin`, `end`, `rbegin`, `rend` | `peaks().iter()`, `.rev()` | |
| `insert(mass, intensity)` | `insert(IsotopePeak)` | |
| `operator[]` | `peaks()[i]` | read-only; mutation only through validated methods |
| `sort_`, `transform_` (protected) | not ported | internal helpers |
| `std::hash<IsotopeDistribution>` | not ported | `f64` fields |

### `DATASTRUCTURES/DBoundingBox.h` (with `width`/`height` of `DIntervalBase.h`)

| Source member | Rust | Notes |
|---|---|---|
| `DBoundingBox<D>`, `DIMENSION`, `Base`, `PositionType`, `CoordinateType` | `BoundingBox2D`, `Point2D`, `f64` | two dimensions only: 0 = RT, 1 = m/z; `D = 1` not ported |
| `DBoundingBox()` (empty sentinel) | `Option<BoundingBox2D>::None` | |
| copy constructor, `operator=(DBoundingBox)` | `Copy` | |
| `operator=(const Base&)`, `operator==(const Base&)` | not ported | the `DIntervalBase` port is the separate `data_structures::dinterval` type |
| destructor | implicit | |
| `DBoundingBox(minimum, maximum)` | `BoundingBox2D::new` | rejects inverted or non-finite corners; the source normalises |
| `enlarge(PositionType)`, `enlarge(x, y)` | `union` with a point box | `Option` for the empty start |
| `operator==(const DBoundingBox&)` | `PartialEq` | |
| `encloses(PositionType)`, `encloses(x, y)` | `encloses(Point2D) -> Result<bool>` | inclusive |
| `intersects` | `intersects` | new; inclusive |
| `isEmpty` | not ported | see native differences |
| `operator<<` | not ported | `Debug` |
| `std::hash<DBoundingBox<D>>` | not ported | |
| `DIntervalBase::width`, `height` | `width`, `height` | new |

## Preserved source conventions (`ProbabilityPrecision::SourceSingle`)

- **Narrowing points.** Element abundances and isotope-override weights become
  `f32` when read (`ElementDB.cpp:666-679` stores `double` literals into
  `Peak1D(double, float)`; the override path stores `float` on insertion).
  Rust `as f32` rounds to nearest, ties to even, as a C++ `float` assignment
  under the default rounding mode. Every element abundance on the traced path
  narrows to the executed SDK's bits (H, C, N, O, S, P, Br, Se).
- **Accumulation.** `convolve` and `convolveSquare_`
  (`CoarseIsotopePatternGenerator.cpp:347-355, 439-445`) loop `i` descending,
  then `j` descending, and evaluate `p + l[i] * r[j]` as a binary32 product and
  a binary32 sum. The port performs the same two operations in the same order.
  A contracted fused multiply-add would round once; the executed SDK build does
  not contract (measured, below).
- **Element order.** `run` iterates `EmpiricalFormula`'s
  `std::map<const Element*, SignedSize>`, so the convolution order is pointer
  order. The executed SDK iterates in `ElementDB` construction order (ascending
  atomic number: H, C, N, O, P, S, ..., Br). The port sorts formula atoms by
  atomic number, each natural element before its labelled isotopes in ascending
  mass number, and sums `getLightestIsotopeWeight` in that order too. Natural
  elements are measured; the labelled-isotope position follows
  `ElementDB::buildElement_` (the natural element is allocated before
  `storeIsotopes_`) and is source review only.
- **Exponentiation.** `convolvePow_` returns its input for exponent 1, starts
  from the input (odd exponent) or the identity, and convolves each successive
  square whose bit is set, lowest bit first. Squares the source computes after
  the highest set bit are unused and are not computed.
- **Truncation.** The source keeps squares at `max_isotope + 1` bins and never
  truncates gap-filled inputs; the port truncates every temporary to
  `max_isotope` bins. A retained bin `k < max_isotope` receives the products
  `(i, k - i)` for `i` from `k` down to 0 in both layouts, and the extra bin
  only feeds bins at or beyond `max_isotope`, so every retained value is
  identical. The executed cases with limits 3, 5, 7, 10, 11, 20 and unbounded
  confirm this bit for bit.
- **Charge.** The source convolves `H^charge` also for charge zero; that is a
  convolution with the identity, which copies every binary32 value, and is
  skipped. Positive charge adds hydrogen atoms after the formula elements, as in
  the source.
- **Renormalization.** `IsotopeDistribution.cpp:176-192`: the binary32 weights
  are summed from the end into a `double`, and each weight becomes the `float`
  narrowing of `weight / sum`. `renormalize_with(SourceSingle)` exposes the same
  arithmetic on any distribution.
- **Fragment weights.** `calcFragmentIsotopeDist_` accumulates complementary
  weights into a `float` in ascending precursor index and multiplies the sum by
  the `float` fragment weight. `estimate_fragment_from_weights` keeps the
  generator's overrides and precision in the inner solver and anchors masses at
  the source-order lightest-isotope mass.
- **`trimLeft`.** `IsotopeDistribution.cpp:210-220` erases before the first
  weight at or above the cutoff and erases nothing if none reaches it.
  `trim_left_source` does the same. `trimRight` already removes every peak when
  none reaches the cutoff, and `trim_right` matches it.
- **Bounding boxes.** `DBoundingBox::intersects` (`DBoundingBox.h:157-166`)
  returns false as soon as the other box starts above or ends below this one in
  a dimension, checking RT then m/z, so touching borders intersect. `width` and
  `height` are `max - min` of dimensions 0 and 1 (`DIntervalBase.h:319-329`).

## Native differences

- **Precision is opt-in.** `CoarseIsotopePatternGenerator::new` and `Default`
  use `ProbabilityPrecision::Double`, whose results are unchanged by this
  package. Narrowing discards precision, so it is never the default.
- **Stored type.** Results stay `IsotopePeak { f64, f64 }`. Under source
  precision every weight is an exact binary32 value, so a consumer can narrow
  back to `f32` without loss.
- **Errors instead of undefined values.** A weight beyond the binary32 range,
  a binary32 overflow during convolution, a zero or non-finite renormalization
  sum and a negative charge return `Error::InvalidValue`; the source produces
  infinities, NaN or `Exception::Precondition`. Failed renormalization leaves
  the distribution unchanged.
- **Override validation.** `set_isotope_override` requires declared isotope
  mass numbers, a strictly increasing nominal layout that includes the lightest
  declared isotope, and a positive finite total. The
  `FeatureFinderAlgorithmPicked.cpp:163-179` override inserts into a
  default-constructed distribution that already holds `(0, 1)`; the executed C++
  keeps that stray peak and produces different, longer patterns (windows
  0/1/5/10 have 30/110/436/811 bins instead of 6/26/148/247). The port rejects
  that input rather than reproduce it. Whether FeatureFinderCentroided should
  reproduce the defect or use the intended two-isotope override is a B6/B10
  decision; see the C++ issue candidate below.
- **`convolve_power` for exponent 1** returns the gap-filled nominal layout,
  as for every other exponent; the protected source helper returns its input
  unchanged. Probabilities are equal.
- **Allocation and work limits** (`MAX_ISOTOPE_PEAKS`,
  `MAX_CONVOLUTION_PRODUCTS`) apply in both precisions; source precision counts
  the same products as the native path.
- **`trim_left` keeps its native default** (removes every peak when all are
  below the cutoff); only `trim_left_source` reproduces the source. The source
  comparison promotes `float` to `double`; `trim_left_source` compares the stored
  `f64`, which is the same comparison for source-precision results.
- **Poisson approximations** are not affected by `ProbabilityPrecision`;
  `approximateFromPeptideWeight` keeps a `float` running product in the source.
- **Sorting** is stable where the source `std::sort` is unstable.
- **`BoundingBox2D`** has no empty sentinel and no `isEmpty`: emptiness is
  `Option`, and a zero-width or zero-height box is an ordinary degenerate box
  that `intersects` and `encloses` treat inclusively (the source `isEmpty`
  would call it empty, but `intersects` ignores `isEmpty`). `new` rejects
  inverted corners instead of normalising them. `width` and `height` overflow to
  positive infinity for extreme finite corners, as in the source.

## Checked boundaries and evidence

### Tier 1: executed C++ SDK probe

`../oracle/b2-iso-source-precision/probe.cpp` links the product SDK
`libOpenMS.dylib` (core 4fdec46, Debug, AppleClang 21, arm64; the four headers
used are hash-identical to the pin) and prints binary64 masses and binary32
`Peak1D` intensities as hexadecimal bit patterns. It ran twice in the fixed
oracle environment with byte-identical output. The 243-line output is the
repository fixture `tests/data/isotopes_source_precision/probe.tsv`, labelled
*oracle-generated (tier 1 executed differential)*. Hashes of the driver, run
script, cross-check, logs, binary and manifest are in the provenance file.

| Probe lines | Rust test | Contract |
|---|---|---|
| `element` (8) | `element_tables_narrow_to_the_source_peak1d_values` | table masses and `f32` abundances |
| `run` (10) | `run_matches_executed_cpp_bits` | glucose (0 and +2 charge), C6H14O6, C222N190O110, Br2, CBr2, C160, H2, rounded C100, unbounded C2 |
| `estimate` (10) | `averagine_estimates_match_executed_cpp_bits` | peptide/RNA/DNA at 100, 1000, 10000 Da (limit 3); peptide 1234.2 Da unbounded (317 bins) |
| `window`, `window_trim` (81 each) | `feature_finder_picked_windows_match_executed_cpp_bits_and_trims` | FeatureFinderAlgorithmPicked windows 0..80 (limit 20, width 100) and source `trimLeft`/`trimRight` at 0.1 % |
| `override`, `override_window` | `override_weights_narrow_like_source_insertion_and_match_cpp_windows` | 12C = 90 % override at limit 1020; the stray `(0, 1)` construction |
| `convolve` (2) | `raw_convolution_matches_executed_cpp_bits` | public convolution, identity case |
| `trim`, `renormalize` | `source_trim_left_keeps_every_peak_when_none_reaches_the_cutoff` | all-below cutoff; binary32 renormalization |
| `fragment` (3) | `fragment_weights_match_executed_cpp_bits` | peptide fragments, `calcFragmentIsotopeDist` on C1/C2 |
| `bbox` (21), `bbox_touch`, `bbox_extent` | `tests/geometry_bounding_box.rs` | intersects sequence, touching borders, width/height |

Every source-precision result is also checked to hold only exact binary32
values. Comparisons are bitwise. The probe ran on macOS arm64; the Rust tests
ran on Linux x86_64 (kim) and agree bit for bit, as expected for correctly
rounded IEEE-754 `+`, `*`, `/` and narrowing with no `libm` call on the path.

### Tier 3: upstream class tests (`TEST_REAL_SIMILAR` rule)

`tests/isotopes_source_precision.rs` implements ClassTest `isRealSimilar`
(`ClassTest.cpp:364-490`: absolute difference at most 1e-5 or ratio at most
1 + 1e-5) and asserts the class-test literals under source precision.

- `CoarseIsotopeDistribution_test.cpp`, 23 sections: the three constructors,
  setMaxIsotope (317 bins), convolve_, run (glucose, charged glucose, explicit
  adduct, negative charge), convolvePow_ (C222N190O110, Br2, CBr2),
  estimateFromWeightAndComp, estimateFromPeptideWeight (probabilities, masses,
  rounded masses), approximateFromPeptideWeight and approximateIntensities (KL
  below 0.05 against the source-precision truth), estimateFromPeptideWeightAndS,
  estimateFromRNAWeight, estimateFromDNAWeight, the three
  estimateForFragmentFrom*Weight sections and calcFragmentIsotopeDist. Not
  applicable: destructor, setRoundMasses (construction), the two NOT_TESTABLE
  getters. Not ported: estimateForFragmentFromPeptideWeightAndS.
- `IsotopeDistribution_test.cpp`, 25 sections: constructor, `==`, getMax,
  getMin, getMostAbundant, size, clear, trimRight, trimLeft (through
  `trim_left_source`), renormalize and `!=` are asserted; `operator<` is not
  ported; copy/assignment/set and the eight NOT_TESTABLE iterator sections need
  no source-precision variant (existing `tests/isotopes.rs` covers the native
  container).
- `DBoundingBox_test.cpp`, 18 sections: intersects (all 21 assertions), both
  encloses sections, both enlarge sections and equality (through `union` and
  `PartialEq`). Not ported: `D = 1` constructor, isEmpty, the `Base`
  assignment/equality, `operator<<` and hash. Constructors, destructors and
  copies need no test.
- `DIntervalBase_test.cpp`: the width (60) and height (75) sections.

The existing native ports in `tests/isotopes.rs` (both class tests) and every
other isotope or geometry consumer stay green: see the verification commands
below.

### Tier 4: independent derivation and native invariants

- `../oracle/b2-iso-source-precision/crosscheck.py` re-derives every `run`,
  `estimate` and `window` line with exact rational arithmetic rounded to
  binary32/binary64 in Python. 93 lines match only unfused accumulation, 8 match
  both, none match only a fused multiply-add and none match neither. With
  alphabetical element order, 78 of the 101 lines fail, so the element order is
  load-bearing.
- `intersects_equals_interval_overlap_on_an_exhaustive_grid` compares 10,000
  box pairs with the interval-overlap definition and checks symmetry.
- Binary32 overflow, weights beyond the binary32 range, zero sums, work and
  allocation limits, invalid cutoffs and the default precision are asserted in
  `source_precision_rejects_values_outside_binary32_and_keeps_limits` and
  `source_precision_is_explicit_and_the_default_is_unchanged`.

### Pending for B6

`FeatureFinderAlgorithmPicked` itself (window count from the MS1 maximum m/z,
optional-peak classification, maximum scaling, `trimmed_left`) is not part of
this package. B6 must show that `estimate_from_peptide_weight` under
`ProbabilityPrecision::SourceSingle`, followed by `trim_left_source` and
`trim_right`, equals the C2-ORACLE-PICKED `TheoreticalIsotopePattern` (item 5 of
`../oracle/featurefinder-picked`) bit for bit for every FeatureFinderCentroided_1
window (charge 2, width 100) and for the 12C = 90 % case. C2 had not landed when
this package was built; the windows 0..80 above cover FFC_1 masses up to
8050 Da, which includes every window whose maximum MS1 m/z times charge 2 stays
below that bound.

### Verification commands

Run from the worktree root on the kim build node through
`~/.local/bin/openms-kim-gate.sh` (Rust 1.96 and 1.85.0), plus local `cargo fmt`
and the repository checkers. The exact results are recorded in
`target_verification` of the provenance file.

```text
cargo fmt --all -- --check
cargo test --locked --all-features --lib --test isotopes_source_precision \
  --test geometry_bounding_box --test isotopes --test features --test fine_isotopes \
  --test fine_isotope_reference --test fine_isotope_stream --test fine_isotope_stream_review \
  --test theoretical --test rna_workflow --test adduct_workflow --test averagine_deisotoping \
  --test averagine_deisotoping_reference --test averagine_deisotoping_review --test deisotoping \
  --test feature_finding_metabo --test feature_identification --test feature_hypothesis \
  --test featurexml --test map_operations --test alignment_transformer --test ims_isotopes \
  --test ims_element --test ims_alphabet
cargo +1.85.0 test (same targets)
cargo test --locked --no-default-features --test isotopes_source_precision --test geometry_bounding_box
cargo clippy --locked --all-features --lib --tests -- -D warnings
RUSTDOCFLAGS='-D warnings' cargo doc --locked --all-features --no-deps
python3 tools/check_core_sdk.py
python3 tools/check_module_cycles.py
python3 tools/check_doc_coverage.py
```

## C++ issue candidates

- **FeatureFinderAlgorithmPicked abundance override keeps a stray `(0, 1)`
  peak.** `FeatureFinderAlgorithmPicked.cpp:163-179` builds the 12C and 14N
  overrides with `IsotopeDistribution isotopes; isotopes.insert(...)`, but the
  default constructor already holds `(0, 1)` (`IsotopeDistribution.cpp:33-36`).
  The override therefore places most weight 12 Da below carbon-12, and every
  non-default `isotopic_pattern:abundance_12C` or `abundance_14N` produces a
  wrong theoretical pattern (executed: 30 instead of 6 bins for the first
  window at 12C = 90 %). A `set()` of the two-isotope container, as before the
  override refactoring, would give the intended pattern.
- **Coarse pattern bits depend on `Element*` addresses.** `run` convolves in
  `std::map<const Element*, SignedSize>` order, which is heap-address order.
  It happens to follow `ElementDB` construction order in the executed build, but
  a different allocator could change the binary32 patterns (78 of 101 probe
  lines change under a different element order).

## Consumers

`FeatureFindingMetabo`, `Deisotoper`, theoretical spectra, fine isotopes, RNA
and adduct workflows construct generators through `new` and therefore keep the
`f64` default. `Feature::hull_bounding_box`, `ConvexHull2D::bounding_box` and the
featureXML reader and writer use `BoundingBox2D` without the new predicates.
