# Issues found in the original OpenMS C++ SDK

This is the ongoing upstream issue log for the Rust port. Findings below were
checked against **OpenMS4-core `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`**.
They are open in that source revision. **Fixes are proposed C++ changes, not
claims that the upstream repository has been patched.** Rust behavior is noted
separately. Findings first logged here: 2026-09-11.

“Executed” means a focused C++ reproduction ran. “Source-reviewed” means the
fault follows from the inspected source and stated input; it does not claim a
C++ runtime reproduction. A native regression alone is not C++ execution.
Intentional source behavior, ordinary API differences and unconfirmed suspicions
are not classified as confirmed defects here.

| ID | Issue | Evidence | Status |
| --- | --- | --- | --- |
| CPP-001 | DateTime silently replaces early dates after failed libc conversion | Executed, macOS | Open |
| CPP-002 | DateTime fractional normalization has signed integer overflow | Executed, UBSan | Open |
| CPP-003 | Stale system configuration discards repaired defaults | Source-reviewed | Open |
| CPP-004 | Isotope left-trimming retains every peak when all are below cutoff | Source-reviewed | Open |
| CPP-005 | Fragment isotope truncation writes beyond the result container | Source-reviewed | Open |
| CPP-006 | Non-preprocessed interpolation prepends unwanted zero coordinates | Source-reviewed | Open |
| CPP-007 | EMG evaluation produces NaN from finite tail inputs | Source-reviewed; independent numerical regression | Open |
| CPP-008 | ProForma modified ranges silently lose residue modifications | Source-reviewed | Open |
| CPP-009 | ProForma ambiguous mass validation ignores candidate modifications | Source-reviewed | Open |
| CPP-010 | IntegerMassDecomposer can hang on sorted positive weights | Source-reviewed; independent recurrence and arithmetic | Open |
| CPP-011 | Reusing MassTraceDetection retains stale array indices | Source-reviewed | Open |
| CPP-012 | Smoothed peak area mixes raw and smoothed intensities | Source-reviewed | Open |
| CPP-013 | FeatureFindingMetabo divides by zero total intensity | Source-reviewed | Open |
| CPP-014 | MassDecomposition addition can lower its maximum-count cache | Source-reviewed | Open |
| CPP-015 | Label-only cross-link endpoints make mass depend on chain order | Source-reviewed | Open |
| CPP-016 | mzML counting can decode binary arrays for spectra without RT | Source-reviewed | Open |
| CPP-017 | Skipping chromatograms suppresses unrelated mzML callbacks | Source-reviewed | Open |
| CPP-018 | Centroid inspection fails to restore reader options after an error | Source-reviewed | Open |
| CPP-019 | mzML writer overstates the number of processing records | Source-reviewed | Open |
| CPP-020 | ProForma try-mass validates a different resolved copy than it calculates | Source-reviewed; public AST trigger | Open |
| CPP-021 | CV XML output substitutes the wrong unit or dereferences an empty set | Source-reviewed | Open |
| CPP-022 | Legacy binary-type xrefs retain part of their prefix | Source-reviewed | Open |
| CPP-023 | Vocabulary stream output sends parent lines to global stdout | Source-reviewed | Open |
| CPP-024 | CV XML output leaves several attribute values unescaped | Source-reviewed | Open |
| CPP-025 | Processing-action fallback state leaks between mzML methods | Source-reviewed | Open |
| CPP-026 | mzML writer gives every processing step the same order value | Source-reviewed; schema contract | Open |
| CPP-027 | mzML writing discards processing completion seconds | Source-reviewed | Open |
| CPP-028 | Recognized software metadata reaches a missing mzML mapping path | Source-reviewed | Open |
| CPP-029 | Annotation-only brackets pass AASequence checks but fail conversion | Source-reviewed | Open |
| CPP-030 | Strict AASequence conversion silently drops terminal crosslinks | Source-reviewed; public AST trigger | Open |
| CPP-031 | Empty ambiguous regions shift AASequence attachment indices | Source-reviewed; public AST trigger | Open |
| CPP-032 | CV mapping namespace stripping rejects plain segments and loses attribute markers | Source-reviewed | Open |
| CPP-033 | Failed CV mapping loads contaminate later loads on the same reader | Source-reviewed | Open |
| CPP-034 | CV reference bulk assignment invalidates an aliased input iterator | Source-reviewed | Open |
| CPP-035 | Invalid CV mapping enum values silently change validation rules | Source-reviewed | Open |
| CPP-036 | XML compression sniffing reads uninitialized bytes from short files | Source-reviewed | Open |

## CPP-001 — DateTime ignores failed calendar conversion

**Affected file:** [`src/openms/source/DATASTRUCTURES/DateTime.cpp`, lines 48–83](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/DATASTRUCTURES/DateTime.cpp#L48), `addSecsToFields_`.

**Issue and reproduction:** `set("0001-01-01T00:00:00"); addSecs(0)` changes the
value to `1969-12-31 23:59:59` in the executed macOS probe. Adding zero should
preserve the valid year-one date. The helper accepts `timegm`'s failure result
and continues through `gmtime_r` without checking either conversion. A separate
libc diagnostic returned `time_t(-1)` with **errno still zero** for sampled input
years -1, 0, 1, 400, 1800 and 1899; 1900 and later sampled years converted.

**Evidence:** Unmodified source DateTime implementation/header compiled with
Apple clang 21, arm64-apple-darwin25.5.0, with UBSan. Of 301 ordinary probe rows,
35 expose this early-year arithmetic departure; the raw results are retained.
This is a host-dependent failure, not a claim that every libc rejects these dates.
[Durable source/probe evidence](tests/data/datetime_cpp_issues.json),
[raw source cases](tests/data/datetime_cpp_probe.tsv), and
[libc diagnostic](tests/data/datetime_libc_diagnostic.tsv) retain the results.

**Proposed fix:** Use checked, timezone-independent Gregorian arithmetic for
naive date addition. If retaining libc conversion, enforce and document its
supported range and check conversion results before publishing fields; checking
errno alone is insufficient on this host. Keep milliseconds and validity intact,
and leave the original value unchanged on failure. Add early-year, zero-date,
partial-date, leap-century and zero-increment regressions.

**Rust handling:** The native DateTime port uses deterministic checked Gregorian
normalization. All 35 differing rows are tested against independent calendar
month-stepping expectations; they are not counted as C++ matches.

## CPP-002 — Fractional seconds overflow signed int during normalization

**Affected file:** [`DateTime.cpp`, lines 230, 255 and 586](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/DATASTRUCTURES/DateTime.cpp#L230), the duplicated fractional parsing branches in `set` and `fromString`.

**Issue and reproduction:** Parsing `2000-01-01T00:00:00. 2147483647` assigns a
valid `int` through `%d`. The digit counter sees whitespace immediately after
the dot and counts zero digits; normalization then evaluates `2147483647 * 10`
in signed `int`, which is undefined behavior. The same input in the explicit
millisecond format reproduces the other branch; appending `+00:00` reproduces
the timezone-stripping branch.

**Evidence:** All three cases abort the unmodified-source UBSan probe with a
signed multiplication overflow diagnostic at the respective source lines.
[Exact inputs, UBSan stderr and source/compiler hashes](tests/data/datetime_cpp_issues.json)
are retained separately from the 301 ordinary probe rows.
[Probe instructions](tools/datetime_probe/README.md) describe the unmodified
source body and narrow dependency adapters.

**Proposed fix:** Share one checked normalization helper. Scale in `int64_t`
before applying the intended signed remainder, or reject an unsupported fraction
before multiplication. Also check conversion of out-of-range `%d` inputs.
For the whitespace case above, wide mathematical scaling followed by `% 1000`
gives zero; any rejection policy should be explicit and regression-tested.

**Rust handling:** Checked integer parsing/scaling rejects this input. Full `set`
retains source clear-first error behavior; partial setters remain atomic.

## CPP-003 — Configuration repair is computed and discarded

**Affected file:** [`src/openms/source/SYSTEM/FileConfig.cpp`, lines 128–160](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/SYSTEM/FileConfig.cpp#L128), `File::getSystemParameters`.

**Issue and reproduction:** Supply a readable `OpenMS.ini` containing only an
old `version` and a custom value. The function announces repair, builds defaults
in `p_new`, and merges existing settings into it, but returns `p`. Only the
version changes; missing defaults such as `threads` remain missing. The same
fault affects a missing version entry.

**Evidence:** The local source branch never assigns or returns `p_new` after
`p_new.update(p)`. The native compatibility regression
[`source_stale_config_returns_original_tree_with_only_version_updated`](tests/system_file.rs)
asserts the observed source-derived tree shape; no C++ runtime reproduction is
claimed for this finding.

**Proposed fix:** Return or move the repaired `p_new` into `p` after merging,
preserving the existing no-file-rewrite policy. Test both stale and missing
versions, retained custom values and presence of all five defaults.

**Rust handling:** Currently preserves this source behavior with warnings, as
documented in [system-file support](docs/SYSTEM_FILE_SUPPORT.md).

## CPP-004 — trimLeft does nothing when every peak is below cutoff

**Affected file:** [`src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.cpp`, lines 210–220](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.cpp#L210).

**Issue and reproduction:** For intensities `[0.1, 0.2]`, `trimLeft(0.5)` retains
the entire distribution. Erasure occurs only after finding an intensity at least
the cutoff, so the all-below case never erases anything. The correct retained
suffix is empty.

**Evidence:** Direct inspection of the complete loop. Native
[isotope regressions](tests/isotopes.rs) cover the all-below case and equality at
the cutoff; no C++ runtime reproduction is claimed.

**Proposed fix:** Find the first qualifying peak and erase the prefix up to that
iterator, allowing the iterator to equal `end()`. Add empty, all-below,
first-qualifying and exact-cutoff cases.

**Rust handling:** Clears the distribution when every intensity is below cutoff.

## CPP-005 — Truncated fragment distribution is indexed at full input length

**Affected file:** [`src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.cpp`, lines 450–520](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.cpp#L450), `calcFragmentIsotopeDist_`.

**Issue and reproduction:** Use a fragment distribution of three peaks, a
nonempty complementary distribution, and `max_isotope_ = 1`. The result is
resized to one peak, but the accumulation loop still visits all three input
peaks and accesses `result[1]` and `result[2]`. The final multiplication makes
the invalid access possible even when no precursor isotope matches.

**Evidence:** Result length is bounded by `r_max` at lines 461–472; the loop at
509–519 is instead bounded by `fragment_isotope_dist.size()`. This is a
source-reviewed out-of-bounds defect; a full C++ sanitizer run is not claimed.

**Proposed fix:** Bound every result-writing loop by `r_max` (or `result.size()`)
while retaining the existing complementary-index checks. Add sanitizer-backed
regressions where the configured length is shorter than the fragment input,
including empty precursor selections.

**Rust handling:** Bounds the full accumulation loop by the result length. See
[isotope support](docs/ISOTOPE_SUPPORT.md).

## CPP-006 — Interpolation overload resizes and then appends coordinates

**Affected file:** [`src/openms/source/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.cpp`, lines 209–234](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.cpp#L209), the `vector<pair<double,double>>` constructor with `preprocess = false`.

**Issue and reproduction:** With sorted anchors `(1,10), (2,20), (3,30)`, the
branch first resizes each coordinate vector to three zero entries, then appends
the three actual values. The resulting x vector is `[0,0,0,1,2,3]`, introducing
false anchors and duplicate coordinates into the interpolation backend.

**Evidence:** The same branch contains both `resize(data.size())` and
`push_back`. This source-reviewed defect concerns the alternate pair-vector
overload; it does not apply to the primary preprocessed annotated-point path.

**Proposed fix:** Use `reserve` followed by appending, or `resize` followed by
indexed assignment. Test both constructor overloads with preprocessing disabled
and assert identical anchor counts and interpolation results.

**Rust handling:** Uses the primary reviewed interpolation path and does not
reproduce this allocation defect. See [transformation support](docs/TRANSFORMATIONS_SUPPORT.md).

## CPP-007 — EMG tail expression overflows before its asymptotic branch

**Affected file:** [`src/openms/source/MATH/MISC/EmgGradientDescent.cpp`, lines 419–441](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/MATH/MISC/EmgGradientDescent.cpp#L419), `emg_point`; related derivative expressions in the same file also require review.

**Issue and reproduction:** For `x=-40`, `h=1`, `mu=0`, `sigma=1`, `tau=1`, the
positive-z branch is selected. It evaluates `exp(z*z)` with `z=41/sqrt(2)`, whose
exponent is 840.5 and overflows binary64. Multiplication by underflowed factors
then produces NaN although the finite EMG tail should underflow cleanly to zero.
The switch to the asymptotic branch at `z > 6.71e7` comes far too late to prevent
this intermediate overflow.

**Evidence:** Source expression inspection and the independent scalar/numerical
review in [EMG reference review](docs/EMG_REFERENCE_REVIEW.md), with the finite
input regression in [tests/emg.rs](tests/emg.rs). No executed C++ comparison is
claimed for this entry.

**Proposed fix:** Evaluate the positive-z expression using a scaled complementary
error function (`erfcx`) or a stable log-domain/asymptotic formulation, with
thresholds chosen for binary64 range. Apply the same stability review to loss
derivatives. Test ordinary tails, subnormal results and both branch transitions
against an independent high-precision reference.

**Rust handling:** Detects the nonfinite source-style result and returns an error;
it does not yet substitute a different scientific EMG formula.

## CPP-008 — ProForma modified ranges omit their residue annotations

**Affected file:** [`src/openms/source/CHEMISTRY/ProForma.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2113): resolution at 2113–2117, mass accumulation at 1958–1965, and mass issue checking at 2521–2532. Parser handling at 1145–1175 stores these inner annotations.

**Issue and reproduction:** `(M[UNIMOD:35]A)[+1]` retains the oxidation annotation
on the inner methionine in its parsed representation and text output. The three
scientific visitors process the range's own annotation but skip annotations on
`range.elements`. Oxidation is therefore neither resolved nor validated nor
included in the calculated mass, causing silent mass loss.

**Evidence:** Chemistry source review traced the same represented fields through
the parser, writer and all three omitted traversals. This is not yet an executed
C++ scientific reproduction.

**Proposed fix:** Visit each range element's modifications in resolution,
validation and accumulation before processing the range-level modifications.
Keep positions and cross-link counting consistent with ordinary elements. Add
modified-range examples with known, unknown and multiply annotated residues.

**Rust handling:** The native resolver and mass operations preserve this source
omission explicitly. [Mass support](docs/PROFORMA_MASS_SUPPORT.md) documents the
policy; [direct regressions](tests/proforma_mass.rs) retain the affected range.

## CPP-009 — Ambiguous mass checks ignore modifications on candidates

**Affected file:** [`ProForma.cpp`, lines 2506–2519](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2506), ambiguous-region mass validation; accumulation at 1949–1956.

**Issue and reproduction:** `(?I[UnknownMod999]L)` has equal unmodified candidate
masses, so mass validation reports no issue despite the unresolved annotation.
Accumulation uses only the first candidate and silently omits its unresolved
mass. `(?I[+10]L)` likewise passes the equal-base-mass check although its modified
candidate masses differ; the reported mass depends on which candidate is first.

**Evidence:** Source review shows that the validation loop compares residue
masses only and never calls the modification checker, while accumulation uses
only the first element. This is not yet an executed C++ scientific reproduction.

**Proposed fix:** Validate every candidate's modifications and compare complete
candidate masses. Report unresolved or mass-ambiguous candidates deterministically
instead of choosing a mass through candidate order. Add equal-base/different-mod,
unresolved-mod and truly equal modified-mass regression cases.

**Rust handling:** Resolution visits ambiguous-region annotations, matching the
source resolver. Native mass validation and accumulation preserve the source
omissions, with [direct regressions](tests/proforma_mass.rs) and an explicit
[compatibility policy](docs/PROFORMA_MASS_SUPPORT.md).

## CPP-010 — Integer mass decomposition can loop without progress

**Affected file:** [`src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IntegerMassDecomposer.h`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IntegerMassDecomposer.h#L355), witness creation at 355–367 and consumption at 420–430.

**Issue and reproduction:** Construct sorted positive weights `[10,16,25]` at
precision 1 and request the decomposition of 73. The table reports existence;
the exact composition is `[0,3,1]`, since `3*16+25=73`. However, the witness for
residue 3 is `(index 2, count 0)`. The single-result loop subtracts `0*25`, leaves
mass 73 and residue 3 unchanged, and repeats forever.

**Evidence:** An independently executed translation of the pinned table loops
produces final row `[0,41,32,73,64,25,16,57,48,89]` and the zero-count witness.
A separate exhaustive integer-composition calculation confirms the unique
solution. [Reproduction](tools/probes/ims_witness_source_oracle.py),
[results](tests/data/ims_witness_source_oracle.json) and
[source provenance](tests/data/ims_witness_source_oracle_provenance.json) are
retained. This is source-derived Python evidence, **not C++ execution**.

**Proposed fix:** Update the witness counter for each individually relaxed
residue; currently only the first counter advances before a loop can relax
several residues. Check against exhaustive small sorted alphabets. Also reject
zero-progress or invalid witnesses in `getDecomposition` to prevent a hang if
table invariants fail.

**Rust handling:** Checked witness traversal rejects nonprogressing state.
Earlier native regressions cover the equivalent unsorted 33-Da example; the
durable independent reproduction here also covers the sorted 73-Da case.

## CPP-011 — Mass trace detection reuses stale metadata-array state

**Affected files:** [`src/openms/source/FEATUREFINDER/MassTraceDetection.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FEATUREFINDER/MassTraceDetection.cpp#L80), discovery at 80–144, invocation at 482–485 and reads at 518–531; [`MassTraceDetection.h`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/FEATUREFINDER/MassTraceDetection.h#L181), member initialization at 181–184 and 204–206.

**Issue and reproduction:** Reuse one detector. First process spectra with an
aligned `FWHM_ppm` array (`Constants::UserParam::FWHM_MZ_ppm`), then process
nonempty MS1 spectra with detectable apices and no float arrays. Discovery sets
availability flags when an array is found, but never clears missing arrays on a
later run. Its consistency check accepts zero matching spectra. The later apex
read therefore indexes an absent array using the stale index, causing an
out-of-bounds access. Centroid and FWHM ion-mobility arrays have the same path.

**Evidence:** Direct review of initialization, discovery, validation and all
three consuming reads. No C++ sanitizer execution is claimed.

**Proposed fix:** Reset all three flags and indices before every discovery, or
discover into fresh local state and publish it after validation. Test repeated
runs with arrays, without arrays and with changed array order for each name.

**Rust handling:** Derives fresh per-call array state and tests sequential
inputs; checked failures preserve detector/output state. See
[mass trace detection support](docs/MASS_TRACE_DETECTION_SUPPORT.md).

## CPP-012 — Smoothed area accumulation uses raw peak intensities

**Affected files:** [`src/openms/source/KERNEL/MassTrace.cpp`, lines 61–78](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/KERNEL/MassTrace.cpp#L61), `computeSmoothedPeakArea`; public description in [`MassTrace.h`, lines 241–247](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/KERNEL/MassTrace.h#L241).

**Issue and reproduction:** RTs `[0,1]`, raw intensities `[10,20]` and smoothed
intensities `[1,2]` produce area 10.5 instead of the smoothed trapezoidal area
1.5. Only the first previous intensity comes from the smoothed vector. The
current value and subsequent previous values come from raw peaks, while interval
inclusion still tests the smoothed value. This contradicts the operation's
documented smoothed-area purpose even when all intensities are positive.

**Evidence:** Direct expression evaluation; no C++ runtime reproduction claimed.
The source's separate policy of including only positive-current intervals is
not the defect asserted here.

**Proposed fix:** Use smoothed intensities for both ends of each included
interval and document the nonpositive-value policy. Guard empty input before
indexing element zero. Add distinct raw/smoothed, empty and nonpositive tests.

**Rust handling:** Preserves the finite source expression for compatibility,
with checked empty input, and documents it in
[mass trace support](docs/MASS_TRACE_SUPPORT.md). A scientific correction must
be explicit and must not silently change existing reference results.

## CPP-013 — Feature finding does not guard zero normalization intensity

**Affected file:** [`src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp#L759), scores at 759–778 and total intensity at 911–915.

**Issue and reproduction:** Call the public feature finder directly with a
nonempty zero-intensity mass trace, raw intensity selected and isotope filtering
disabled. Summed intensity is zero, and the mandatory singleton hypothesis
receives `0/0` as its score. The API should reject an undefined normalization or
return a documented zero-signal result, rather than publishing NaN scores.
This finding concerns direct public input; it does not claim that the normal
upstream detection chain generates this trace.

**Evidence:** Source review of accumulation and unconditional division. No C++
execution is claimed.

**Proposed fix:** Validate the total before hypothesis creation and specify the
zero-signal policy. Add direct zero-intensity tests and define handling for
nonfinite totals and cancellation if signed quantified input is supported.

**Rust handling:** Rejects an invalid total before publishing output, consuming
IDs or changing caller-visible processing state. See
[feature finding support](docs/FEATURE_FINDING_METABO_SUPPORT.md).

## CPP-014 — Addition can reduce the cached maximum residue count

**Affected file:** [`src/openms/source/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.cpp`, lines 177–200](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.cpp#L177), `operator+`; compare `operator+=` at 70–94 and the maximum getter at 203–205.

**Issue and reproduction:** `MassDecomposition("A1") + MassDecomposition("B10 C5")`
stores A1/B10/C5 but reports maximum 5. The corresponding `+=` reports 10.
For each new key, `operator+` compares against the original left operand's
maximum, then overwrites the result's running maximum. A later count can
therefore lower it. `A10 + A100 B11` similarly reports 11 despite containing A110.

**Evidence:** Source control-flow review; native regressions pin the finite
behavior. This does not claim C++ execution.

**Proposed fix:** Compare against `d.number_of_max_aa_` at line 186. Assert that
both addition forms agree and the cache equals the maximum stored count for
new-key and existing-key combinations.

**Rust handling:** `checked_add` currently preserves this source behavior;
`checked_add_assign` retains the source's running-maximum update. See
[mass decomposition support](docs/MASS_DECOMPOSITION_SUPPORT.md).

## CPP-015 — Cross-link mass depends on endpoint traversal order

**Affected file:** [`src/openms/source/CHEMISTRY/ProForma.cpp`, lines 1927–1938](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L1927), `addModMass`; modification fallback at 1889–1909 and shared cross-link set at 2589–2598.

**Issue and reproduction:** Compare `K[#XL1]//K[+138.068#XL1]` with the chains
reversed. The first form reserves XL1 at its label-only, zero-mass endpoint,
then skips the chemistry-bearing endpoint as already counted. Reversal includes
the linker contribution. The same two linked chains therefore differ by
138.068 Da solely through their order.

**Evidence:** Source review of label insertion, mass lookup and chain traversal;
no C++ mass execution is claimed.

**Proposed fix:** Determine whether an endpoint supplies chemistry before
reserving its ID. Prefer collecting and validating one definition per link,
then summing once. Preserve explicit zero-mass chemistry and diagnose conflicting
definitions. Test both endpoint orders and annotation-only carriers.

**Rust handling:** Native mass operations preserve first-label reservation and
its endpoint-order dependence, with a [direct regression](tests/proforma_mass.rs)
and an explicit [compatibility policy](docs/PROFORMA_MASS_SUPPORT.md).

## CPP-016 — Count-only mzML loading can still decode peak arrays

**Affected file:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L1401), normal spectrum enqueue at 1401–1427, binary population at 203–217, error propagation at 244 and accepted-RT count handling at 2341–2351.

**Issue and reproduction:** Use `loadSize` with an active MS1 filter, a passing
MS1 spectrum with no scan-start-time term, `fill_data=true`, pool size 1 and
malformed binary data. Accepted RT handling ordinarily increments the count and
marks the spectrum skipped. Without RT, it does neither; the spectrum follows
the normal enqueue/populate path and can throw while decoding peak data despite
producing no count. The final pool flush also exists; this is not a lost-pool
claim.

**Evidence:** Source review tracing the count mode, missing-RT branch, buffer
threshold and binary exception path. No C++ runtime reproduction is claimed.

**Proposed fix:** Gate normal record enqueue and data population to `LD_ALLDATA`.
Define missing-RT counting separately from binary decoding. Test malformed
binary arrays in count mode, with and without RT and with a one-record pool.

**Rust handling:** The dedicated count reader never decodes binary peak arrays.
It retains the reviewed missing-RT count convention. See
[mzML count support](docs/MZML_COUNTS_SUPPORT.md).

## CPP-017 — The chromatogram skip option activates a global callback guard

**Affected file:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L149), `setOptions` at 149 and start-element guard at 870.

**Issue and reproduction:** Set `skip_chromatograms=true` before loading a file
containing experiment headers and spectra. `setOptions` immediately sets
`skip_chromatogram_`, and the generic start callback returns whenever that flag
is true. Header and spectrum callbacks are consequently suppressed even outside
a chromatogram. A chromatogram-only filter should leave those records readable.

**Evidence:** Source review of option assignment and the unconditional shared
guard. No C++ runtime reproduction is claimed.

**Proposed fix:** Keep the option separate from the active-record skip state;
activate the latter only upon entering a chromatogram and reset it at that
record's end. Test files with spectra alone and with both record kinds while
skipping chromatograms, including header retention.

**Rust handling:** Record-scoped handling skips chromatograms while preserving
spectrum processing/counting. See [mzML count support](docs/MZML_COUNTS_SUPPORT.md).

## CPP-018 — Centroid inspection leaves reader options changed after failure

**Affected file:** [`src/openms/source/FORMAT/MzMLFile.cpp`, lines 233–283](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/MzMLFile.cpp#L233), `getCentroidInfo`.

**Issue and reproduction:** Set `FillData(false)`, then call `getCentroidInfo`
with a missing or malformed file and catch the parse/file exception. The method
saves the old setting and enables data loading, but restores the old value only
after `transform` returns normally. An exception leaves `FillData(true)` on the
reused reader, changing subsequent loads unexpectedly.

**Evidence:** Direct source review of option mutation, throwing call and
normal-return restoration. No C++ runtime reproduction is claimed.

**Proposed fix:** Use an isolated reader/options copy for inspection or an RAII
scope guard that restores the original value on every exit. Test success, missing
input and malformed input with both initial settings.

**Rust handling:** This centroid-inspection convenience API remains outstanding.
Existing native load options are borrowed and do not require temporary mutation.

## CPP-019 — Declared processing count includes records that are not written

**Affected file:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, lines 5160–5196](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L5160), `dataProcessingList` count and auxiliary-array history emission.

**Issue and reproduction:** Write one spectrum with one float auxiliary array
and empty processing histories. The ordinary empty history yields one processing
record. Every float array is also counted unconditionally, giving `count="2"`,
but an array's separate record is emitted only if its history is nonempty.
Only `dp_sp_0` is written, so the declared count is wrong.

**Evidence:** Source review of differing count and emission predicates. This
does not claim C++ execution or that the XSD itself enforces count equality.

**Proposed fix:** Count exactly the nonempty array histories that will be emitted,
or build one emission list and use its length. Cover empty and nonempty histories,
multiple arrays and the empty-experiment fallback.

**Rust handling:** Full processing-history header transport is being planned.
The current fixed minimal processing list is separate; complete transport must
derive counts from the emitted records.

## CPP-020 — Formula interning can separate mass validation from calculation

**Affected file:** [`src/openms/source/CHEMISTRY/ProForma.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2642), try-mass copy/validation/calculation at 2642–2647, issue-check copy/resolution at 2491–2492 and formula identity creation at 1661–1673.

**Issue and reproduction:** Construct a public AST containing two ordinary M
elements. Give the first a `NamedMod` named `M[Formula:Cl101]` and the second a
`FormulaTag` with formula `Cl101`, no charge or labels. Begin with that exact
full ID absent from the modification database. This is a directly constructed
AST, **not a claim that text parsing produces that NamedMod**.

`tryGetMonoWeight` first resolves copy A. Its first annotation remains unknown;
the later formula creates the database entry `M[Formula:Cl101]`. Issue checking
then makes and resolves copy B, where the earlier name now resolves successfully,
so validation reports no issue. Calculation nevertheless uses copy A and silently
omits the earlier unresolved mass. It returns base MM plus one formula delta;
a later calculation with the now-populated database includes two deltas.

**Evidence:** Direct review of database insertion, both resolution passes and
the object passed to accumulation. No C++ runtime reproduction is claimed.

**Proposed fix:** Resolve to a stable interpretation, then validate and calculate
that same resolved object without a hidden second resolving copy. If forward
formula aliases are supported, collect definitions before resolving their uses
or repeat resolution explicitly before both operations. Test a fresh registry,
repeated calls and a prepopulated registry with identical ASTs.

**Rust handling:** Native mass operations retain the source's observable pass
order with one checked registry transaction and a [dedicated regression](tests/proforma_mass.rs).
The behavior is documented in [mass support](docs/PROFORMA_MASS_SUPPORT.md)
until a deliberate scientific correction is adopted.

## CPP-021 — CV XML formatting ignores the value's actual unit

**Affected file:** [`src/openms/source/FORMAT/ControlledVocabulary.cpp`, lines 120–135](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/ControlledVocabulary.cpp#L120), typed `CVTerm::toXMLString`.

**Issue and reproduction:** Give a term allowed units `UO:0000010` (second) and
`UO:0000031` (minute), then serialize value 2 with a minute unit. The writer
checks only whether the value has a unit, then substitutes the lexically first
allowed unit: seconds. It neither reads the actual unit nor converts the value,
so it changes the physical meaning. If the term has no allowed units but the
value has a unit, dereferencing `units.begin()` is undefined behavior.

**Evidence:** Direct source inspection of the unconditional first-unit selection
at line 129. No C++ execution is claimed.

**Proposed fix:** Serialize the actual DataValue unit identity. Validate it
against allowed units separately if desired; do not silently relabel it.
Handle an empty constraint set safely. Test multiple allowed units, explicit
units with no constraints and unitless values.

**Rust handling:** The vocabulary port uses the actual `MetaValue` unit
identity, including its name, with checked XML rendering. Allowed-unit metadata
remains a separate vocabulary constraint.

## CPP-022 — Legacy binary xref parsing removes the wrong prefix length

**Affected file:** [`src/openms/source/FORMAT/ControlledVocabulary.cpp`, lines 465–478](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/ControlledVocabulary.cpp#L465), binary data type xref parsing.

**Issue and reproduction:** Read a term with
`xref_analog: binary-data-type: MS:1000523`. The branch accepts both `xref` and
`xref_analog` spellings but always removes 22 bytes after whitespace removal,
the length of the ordinary prefix. It stores `a-type:MS:1000523` instead of
`MS:1000523`, corrupting the referenced accession.

**Evidence:** Source prefix checks and fixed substring offset; no C++ execution
is claimed.

**Proposed fix:** Remove the actual matched prefix including its separator.
Test ordinary and legacy spellings with whitespace, escaped separators and
quoted descriptions.

**Rust handling:** The vocabulary parser explicitly corrects this prefix
handling while retaining the ordinary source branch semantics.

## CPP-023 — Vocabulary printing splits output between two streams

**Affected file:** [`src/openms/source/FORMAT/ControlledVocabulary.cpp`, lines 616–628](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/ControlledVocabulary.cpp#L616), output-stream operator.

**Issue and reproduction:** Write a vocabulary with a parent relationship to an
`std::ostringstream`. Term, ID and name lines reach that stream, but parent
`is_a` lines go to global `std::cout`. The requested output is incomplete and
the operation unexpectedly writes to the process console.

**Evidence:** The inner parent loop names `cout` instead of `os`. No C++ runtime
reproduction is claimed.

**Proposed fix:** Send every line to `os`. Capture the requested stream and stdout
separately in a regression, verifying complete output and no console write.

**Rust handling:** All vocabulary diagnostic output uses the caller's
requested destination or an owned returned string.

## CPP-024 — CV parameter rendering does not escape every XML attribute

**Affected file:** [`src/openms/source/FORMAT/ControlledVocabulary.cpp`, lines 108–135](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/ControlledVocabulary.cpp#L108), both `CVTerm::toXMLString` overloads.

**Issue and reproduction:** Pass a CV reference containing an ampersand, such
as `MS&custom`, to either public formatter. The emitted `cvRef` contains the raw
ampersand and is not well-formed XML. Term IDs and unit identifiers have the
same unchecked concatenation path; names and values already use escaping.
These fields are caller-supplied strings, so safe built-in accessions alone do
not establish the public formatter's correctness.

**Evidence:** Direct review of each attribute concatenation and selective
escaping. No C++ runtime reproduction is claimed.

**Proposed fix:** Escape every attribute value through one shared helper and
reject characters XML 1.0 cannot represent. Add ampersand, quote, less-than,
Unicode and forbidden-control cases for each consumed attribute.

**Rust handling:** Vocabulary XML rendering validates and escapes all
consumed attributes before returning output. Unrelated opaque ontology text
is not rejected merely because it is not suitable for XML.

## CPP-025 — A later processing method can omit its required action term

**Affected files:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, lines 3848–3942](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L3848), `writeDataProcessing_`; required per-method transformation rule in [`share/OpenMS/MAPPING/ms-mapping.xml`, lines 99–101](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/share/OpenMS/MAPPING/ms-mapping.xml#L99).

**Issue and reproduction:** Write two processing methods: the first has peak
picking as an action, and the second has no action. The `written` flag is
initialized outside the method loop and remains true after the first method.
The second therefore omits the intended fallback `MS:1000543` action term.
Its output fails the source mapping's per-method MUST requirement for a data
transformation term. An empty-action method's output incorrectly depends on
the methods preceding it.

**Evidence:** Source review of the flag's lifetime, fallback branch and semantic
mapping scope. No C++ runtime or validator execution is claimed.

**Proposed fix:** Initialize `written=false` inside each method iteration.
Test empty-action methods before, between and after methods with actions,
checking each emitted method against the required semantic rule.

**Rust handling:** The complete header writer will track emitted actions per
method and derive each fallback independently.

## CPP-026 — Processing step order is always written as zero

**Affected files:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, line 3852](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L3852), processing method output; order contract in [`share/OpenMS/SCHEMAS/mzML_1_10.xsd`, lines 460–462](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/share/OpenMS/SCHEMAS/mzML_1_10.xsd#L460).

**Issue and reproduction:** Write two distinct consecutive processing methods.
The loop varies each software reference by index but emits `order="0"` for
every method. The schema describes this required attribute as the way to place
consecutive steps in their correct order, so the output loses distinct ordering
keys. A consumer relying on those keys cannot recover the intended sequence.

**Evidence:** Source and schema inspection. The XSD permits duplicate values,
and OpenMS's reader ignores the order attribute and retains document encounter
order at lines 1264–1276. This is therefore an interoperability defect, **not a
demonstrated OpenMS self-roundtrip loss or schema-validation failure**. No C++
execution is claimed.

**Proposed fix:** Write the actual step index as `order`. Separately review the
reader's handling of external documents with nontrivial ordering keys; define
duplicate-order handling explicitly. Test multi-step output and ordered external
input using the attribute's documented contract.

**Rust handling:** The complete header writer will emit distinct sequential
indices. Source encounter-order reading remains an explicit compatibility
decision until external ordering behavior is separately reviewed.

## CPP-027 — mzML writing discards processing completion seconds

**Affected files:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, line 3947](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L3947), completion-time output; reader at lines 3246–3249; full timestamp contract in [`src/openms/include/OpenMS/METADATA/DataProcessing.h`, lines 111–120](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/METADATA/DataProcessing.h#L111).

**Issue and reproduction:** A processing record completed at `2001-02-03
04:05:37` is written as `2001-02-03+04:05`. Reloading that value reconstructs
`04:05:00`, silently losing the seconds in the stored `DateTime`. No
minutes-only restriction is documented on the completion-time field.

**Evidence:** Source review of the explicit `yyyy-MM-dd+hh:mm` output format,
the completion-time reader and the DateTime parser. This is not an executed
C++ round trip. The existing source fixture uses minute precision and therefore
does not expose the loss.

**Proposed fix:** Emit an accepted timestamp format retaining seconds and any
stored fractional precision. Add a seconds-bearing round trip, with separate
fractional-second coverage where supported.

**Rust handling:** The complete header writer under development will retain
the represented completion timestamp's precision and test the seconds-bearing
case. The published minimal writer still rejects such unsupported headers.

## CPP-028 — Recognized software metadata can throw during mzML writing

**Affected files:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, line 3787](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L3787), software metadata output and validation at lines 3595 and 3668–3673; [`src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp`, line 538](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp#L538); mapping in [`share/OpenMS/MAPPING/ms-mapping.xml`, line 94](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/share/OpenMS/MAPPING/ms-mapping.xml#L94).

**Issue and reproduction:** Add an ordinary software metadata key recognized
by the loaded vocabulary, for example `completion time`, then write mzML.
The writer passes `/mzML/Software/cvParam/@accession` to validation. The mapping
only defines `/mzML/softwareList/software/cvParam/@accession`. Once the known
key reaches `locateTerm`, `rules_.at(path)` throws `std::out_of_range` instead
of emitting an allowed CV term or falling back to a user parameter. Unknown
arbitrary keys bypass this lookup and are not the trigger.

**Evidence:** Direct source review of the known-name branch, validation call,
map lookup and pinned mapping path. No C++ execution is claimed.

**Proposed fix:** Pass the exact software mapping path. Test both a recognized
software term and a recognized term disallowed at that path, verifying CV
promotion or user-parameter fallback respectively without an exception.

**Rust handling:** The complete header writer under development uses the
pinned software mapping path and checked lookup. Header regressions will
cover recognized names and path-specific fallback.

## CPP-029 — Annotation-only brackets pass conversion checks but fail conversion

**Affected file:** [`src/openms/source/CHEMISTRY/ProForma.cpp`, lines 2194–2200](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2194), issue collection; representability at 2268–2270 and strict attachment at 2307–2311.

**Issue and reproduction:** For `M[INFO:note]`, the issue collector deliberately
excludes annotation-only brackets from unresolved chemistry. It reports no
conversion issues and `isRepresentableAsAASequence` returns true. Default
strict conversion nevertheless throws an unresolved-modification error because
the attachment loop rejects every null handle, including this annotation.

**Evidence:** Direct source review of the chemistry predicate and the stricter
attachment branch. No C++ execution is claimed.

**Proposed fix:** Apply the same chemistry predicate during attachment. Define
annotation-loss policy consistently in both the diagnostic and conversion APIs;
if annotations are intentionally droppable, skip them in both. Test INFO-only,
position-only and empty ordinary brackets alongside actual unresolved chemistry.

**Rust handling:** The conversion group under development preserves this
observable source inconsistency and tests it explicitly. A correction must be
documented as a deliberate change to the conversion policy.

## CPP-030 — Strict AASequence conversion silently drops terminal crosslinks

**Affected file:** [`src/openms/source/CHEMISTRY/ProForma.cpp`, lines 2220–2260](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2220), terminal issue collection; terminal attachment at 2337–2346. Ordinary-residue crosslink detection is at 2206–2213.

**Issue and reproduction:** Construct a public peptidoform containing M and
an N-terminal label-only modification whose empty INFO tag carries a CROSSLINK
label. The terminal issue loops check chemistry and alternatives, but never
inspect labels. The issue list is empty, and strict conversion returns M while
discarding the crosslink. The C-terminal path has the same omission. No claim
is made that this directly constructed AST has been obtained by parsing text.

**Evidence:** Source comparison of ordinary versus terminal label handling and
the first-resolved-terminal attachment loops. No C++ execution is claimed.

**Proposed fix:** Include terminal crosslink labels in conversion diagnostics,
so `FAIL_ON_LOSS` rejects them consistently with ordinary-residue links. Test
both termini, label-only brackets and resolved chemistry bearing a link label.

**Rust handling:** The conversion group under development preserves the source
terminal omission with an explicit regression and compatibility note.

## CPP-031 — An empty ambiguous region shifts the conversion attachment index

**Affected files:** [`src/openms/source/CHEMISTRY/ProForma.cpp`, lines 2288–2293](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2288), residue emission; attachment cursor at 2330–2331; bounds check in [`src/openms/source/CHEMISTRY/AASequence.cpp`, lines 1442–1449](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/AASequence.cpp#L1442).

**Issue and reproduction:** Directly construct an empty ambiguous region
followed by a single M with a resolved oxidation modification, then use
`BEST_EFFORT`. The emitted sequence contains one residue because the empty
region contributes none. The attachment pass still increments its cursor for
that region, then attempts to modify index 1 of the one-residue sequence.
`AASequence::setModification` throws `IndexOverflow`. This is a checked C++
exception, not an out-of-bounds memory-access claim. Strict mode rejects the
ambiguous region earlier.

**Evidence:** Source review of the two unequal cursor rules and the actual
AASequence bounds check. No C++ execution or parser-produced-empty-region
claim is made.

**Proposed fix:** Advance the cursor only when an ambiguous region emitted a
residue, or reject empty regions explicitly before constructing the sequence.
Test an empty region before and between ordinary modified residues under both
permissive policies.

**Rust handling:** The conversion group under development preserves the source
index rule with checked errors and a direct public-AST regression.

## CPP-032 — CV mapping namespace stripping mishandles path segments

**Affected files:** [`src/openms/source/FORMAT/CVMappingFile.cpp`, lines 68–102](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/CVMappingFile.cpp#L68), namespace stripping; [`src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h`, lines 603–609](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h#L603), delimiter-free split behavior.

**Issue and reproduction:** Load a mapping with element path
`/mzML/run/@accession` and `strip_namespaces=true`. Splitting an unprefixed
segment at `:` returns one item. The loader expects an empty vector for that
case, treats the one-item result as invalid and raises a parse error. The
option therefore fails on ordinary paths and paths mixing prefixed and plain
segments. Separately, a path containing only prefixed segments such as
`/p:root/@p:accession` becomes `/root/accession`, losing the attribute marker.

**Evidence:** Direct source review of `split` and every reconstruction branch.
No C++ execution is claimed. The unchanged handling of `scopePath` is a
separate source convention, not part of this demonstrated finding.

**Proposed fix:** Accept a single unprefixed segment unchanged. When removing
a namespace prefix, retain any leading `@` that identifies an attribute.
Reject genuinely ambiguous multiple-colon segments explicitly. Test plain,
fully prefixed and mixed element/attribute paths.

**Rust handling:** The CV mapping loader under development applies these
checked corrections, while preserving the source's slash normalization and
separate `scopePath` behavior.

## CPP-033 — Failed CV mapping loads contaminate later loads

**Affected files:** [`src/openms/source/FORMAT/CVMappingFile.cpp`, lines 29–42](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/CVMappingFile.cpp#L29), load/publication/cleanup and callbacks at 52–65; cleanup dispatch in [`src/openms/source/FORMAT/XMLFile.cpp`, lines 41–55 and 116–120](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/XMLFile.cpp#L41); empty inherited [`XMLHandler::reset`, lines 37–39](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/XMLHandler.cpp#L37).

**Issue and reproduction:** Reuse one `CVMappingFile`. First load
`<CvMapping><CvReference cvName="old" cvIdentifier="OLD"/><CvMappingRule/></CvMapping>`.
The reference is accumulated before the rule's missing required `id` raises
an error. Catch it, then load `<CvMapping/>` into a fresh destination. The
second load publishes the stale OLD reference even though its input is empty.
Completed rules and partially accumulated rule terms can survive similarly.

**Evidence:** Source review confirms cleanup is after `parse_`, the inherited
RAII cleaner calls an empty `reset`, and this class supplies no override.
No C++ execution is claimed.

**Proposed fix:** Use local per-load state and publish only on success, or
reset all accumulators, including the current rule, at entry and on every exit.
Test failure after a reference, a completed rule and a partial term sequence,
followed by an empty and a valid load using the same reader.

**Rust handling:** The mapping loader under development uses local parsing
state and atomic destination publication; its reader can be reused after error.

## CPP-034 — CV reference bulk assignment can invalidate its input iterator

**Affected files:** [`src/openms/source/DATASTRUCTURES/CVMappings.cpp`, lines 68–75](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/DATASTRUCTURES/CVMappings.cpp#L68), `setCVReferences`; public const-reference getter at lines 77–80.

**Issue and reproduction:** On a nonempty object call
`m.setCVReferences(m.getCVReferences())`. The input aliases the destination
vector. The loop appends to that vector while traversing it, invalidating its
iterator when capacity grows; later iterator use has undefined behavior.
The public API does not forbid this ordinary getter-to-setter call.

The method also appends rather than replaces references for nonaliased input;
that observable behavior is recorded separately as a compatibility convention,
not the basis for this memory-safety finding.

**Evidence:** Direct source review of the const-reference argument, returned
member reference and vector mutation. No sanitizer/C++ execution is claimed.

**Proposed fix:** Decide replacement versus append semantics explicitly. Build
a separate draft from the supplied range before mutating either stored index
or vector, then publish it atomically. Test self-assignment and duplicate IDs.

**Rust handling:** The owned-input native bulk method cannot alias its internal
reference slice. It preserves source append semantics; a separately named
replacement operation supplies actual replacement.

## CPP-035 — Invalid CV mapping enums silently select different rules

**Affected file:** [`src/openms/source/FORMAT/CVMappingFile.cpp`, lines 103–157](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/CVMappingFile.cpp#L103), requirement-level and combination-logic attribute parsing.

**Issue and reproduction:** Supply otherwise complete rule attributes with
`requirementLevel="MAYY"` and `cvTermsCombinationLogic="ANDD"`. The loader
silently converts these misspellings into MUST and OR. Both unknown-value
branches contain only unimplemented exception comments, so malformed settings
change the resulting validation constraints without reporting the input error.

**Evidence:** Source review of default initialization, accepted literal branches
and empty error branches. No C++ parser/validator execution is claimed.

**Proposed fix:** Raise a parse error for unknown enum values, identifying the
attribute and supplied value. Cover every valid literal, empty values and typos.

**Rust handling:** The mapping loader under development preserves these source
fallbacks with explicit tests and documentation. Stricter rejection would be
a deliberate compatibility change rather than a hidden parser difference.

## CPP-036 — XML compression sniffing reads uninitialized short-file bytes

**Affected file:** [`src/openms/source/FORMAT/XMLFile.cpp`, lines 141–159](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/XMLFile.cpp#L141), the shared XML file loader's compression-prefix probe.

**Issue and reproduction:** Load an empty file or a one-byte non-NUL file.
The loader declares `char tmp_bz[3]` without initialization, requests two bytes
without checking the short read, and initializes only byte 2 before constructing
`std::string(tmp_bz)`. Empty input leaves both prefix bytes uninitialized;
one-byte input leaves byte 1 uninitialized. Constructing the prefix reads that
indeterminate data before XML parsing can report the incomplete document.

**Evidence:** Direct source review of buffer initialization, unchecked read and
C-string construction. No C++ execution or sanitizer result is claimed. This
finding does not additionally claim that reading `std::string[size()]` is out
of bounds; that sentinel access is permitted.

**Proposed fix:** Initialize the prefix buffer and inspect only the bytes
actually read. Require at least two bytes before comparing compression magic.
Test zero/one-byte inputs, ordinary short XML and complete compression headers.

**Rust handling:** Native file readers use initialized buffers and bounded
length-aware input. The new CV mapping transport reuses that path handling;
the C++ defect is in the shared upstream XMLFile implementation, outside its
five-header conversion group.

## Maintaining this log

Add an entry whenever porting or testing identifies a new original C++ defect.
Include the checked source revision, affected files/functions/lines, a concrete
trigger, observed versus expected behavior, evidence level, proposed C++ fix,
regression coverage and the Rust port's handling. Keep stable IDs and update
status only after verifying the upstream fix. Retain resolved entries with the
fixing commit. Link durable probe data once it lands; temporary investigation
paths are not sufficient long-term evidence. Do not promote a suspicion to an
executed or confirmed finding without verification.
