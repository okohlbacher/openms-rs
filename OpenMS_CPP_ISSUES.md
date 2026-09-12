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
| CPP-037 | ProForma modification combination loses charged-formula mass | Source-reviewed; public AST trigger | Open |
| CPP-038 | ProForma crosslink spectra count resolved linker mass twice | Source-reviewed; public AST trigger | Open |
| CPP-039 | SemanticValidator rejects allowed descendant units | Source-reviewed; custom vocabulary trigger | Open |
| CPP-040 | Failed semantic validation leaves stale XML paths and rule counts | Source-reviewed; repeated-call trigger | Open |
| CPP-041 | mzML header strings are inserted into attributes without escaping | Source-reviewed; ordinary text trigger | Open |
| CPP-042 | XLMS linear suffix losses divide the mass by charge twice | Source-reviewed; charge-two trigger | Open |
| CPP-043 | XLMS precursor isotope companions omit precursor charge division | Source-reviewed; charge-two trigger | Open |
| CPP-044 | Missing-path CV lookup changes behavior after validation | Source-reviewed; call-order trigger | Open |
| CPP-045 | PeptideEvidence rejects valid first-residue limits and accepts invalid ranges | Source-reviewed; zero-based range trigger | Open |
| CPP-046 | Semantic date validation applies date-time syntax to xsd:date | Source-reviewed; XSD lexical contract | Open |
| CPP-047 | ProForma XLMS link positions ignore preceding flattened ranges | Source-reviewed; explicit AST trigger | Open |
| CPP-048 | Zero centroid-inspection limit underflows and depends on spectrum type | Source-reviewed; unsigned counter trigger | Open |
| CPP-049 | Indexed mzML output always writes a placeholder file checksum | Source-reviewed; literal footer and schema contract | Open |
| CPP-050 | Empty indexed mzML declares zero indices but emits a dummy entry | Source-reviewed; empty-offset-list branch | Open |
| CPP-051 | mzML validation reuses parameter groups from earlier documents | Source-reviewed; successful repeated-call trigger | Open |
| CPP-052 | Four-line schema detection rejects indexed mzML with a longer XML prolog | Source-reviewed; independently executed XSD checks | Open |
| CPP-053 | Chromatogram primary-array selection can relabel pressure values as flow | Source-reviewed; competing primary-array trigger | Open |
| CPP-054 | Indexed mzML schema does not enforce offset ID references | Source-reviewed; independently executed XSD checks | Open |
| CPP-055 | SIMD Base64 decoding accepts characters outside the Base64 alphabet | Source-reviewed; malformed numeric payload trigger | Open |
| CPP-056 | mzML writing drops spectrum mobility when acquisition records are empty | Source-reviewed; synthetic scan branch | Open |
| CPP-057 | Spectrum settings equality omits both ion-mobility state fields | Source-reviewed; distinct mobility-state trigger | Open |
| CPP-058 | Floating DataValue casts read an inactive union member for strings and lists | Source-reviewed; public numeric-cast trigger | Open |

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

**Rust handling:** [Native centroid inspection](docs/MZML_CENTROID_SUPPORT.md)
borrows caller options and forces population only in a bounded local copy. Tests
cover success, missing/malformed inputs, decoding and classification failures;
the original options and their vector allocation remain unchanged.

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

**Rust handling:** The [native header writer](docs/MZML_HEADER_SUPPORT.md) derives declaration counts from its emitted processing registry, including auxiliary histories and the mandatory empty-history placeholder. Native round trips and independent XSD checks cover both writers.

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

**Rust handling:** The [native header writer](docs/MZML_HEADER_SUPPORT.md) tracks actions per method and emits an independently reversible fallback for each empty action set. Multi-method round trips cover the correction.

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

**Rust handling:** The [native header writer](docs/MZML_HEADER_SUPPORT.md) emits distinct sequential indices. Reading retains source XML encounter order explicitly; tests cover both conventions.

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

**Rust handling:** The [native header writer](docs/MZML_HEADER_SUPPORT.md) retains seconds and milliseconds through DataProcessing DateTime. Native regressions cover mzML (both writers), FeatureXML and ConsensusXML timestamps.

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

**Rust handling:** The [native header writer](docs/MZML_HEADER_SUPPORT.md) uses the pinned software mapping path and checked lookup. Regressions cover recognized names and path-specific user-parameter fallback.

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

**Rust handling:** The [native conversion group](docs/PROFORMA_CONVERSION_SUPPORT.md) preserves this
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

**Rust handling:** The [native conversion group](docs/PROFORMA_CONVERSION_SUPPORT.md) preserves the source
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

**Rust handling:** The [native conversion group](docs/PROFORMA_CONVERSION_SUPPORT.md) preserves the source
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

**Rust handling:** The [native CV mapping loader](docs/CV_MAPPING_SUPPORT.md) applies these
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

**Rust handling:** The [native mapping loader](docs/CV_MAPPING_SUPPORT.md) uses local parsing
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

**Rust handling:** The [native mapping loader](docs/CV_MAPPING_SUPPORT.md) preserves these source
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

## CPP-037 — ProForma combination discards a formula's charge contribution

**Affected files:** [`src/openms/source/CHEMISTRY/ProForma.cpp`, lines 1734–1738](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L1734), formula combination and re-interning; [`src/openms/source/CHEMISTRY/EmpiricalFormula.cpp`, lines 47–54 and 258–276](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/EmpiricalFormula.cpp#L47), charge-dependent mass and charge-free canonical text.

**Issue and reproduction:** Directly construct M with two pre-resolved
modification records and empty alternative lists. Give the records oxygen
formulas with charges +1 and 0 respectively, and set each declared delta mass
equal to its formula's `getMonoWeight()`. Their summed formula and declared
mass agree, so combination takes the formula branch. `toString()` drops the
charge and re-interns neutral `O2`. The combined modification consequently
loses one `PROTON_MASS_U` contribution despite passing the agreement check.

**Evidence:** Source review of charge-dependent mass, formula summation,
canonical serialization and neutral re-parsing. This is a public-AST trigger,
not a claimed parser-produced annotation or executed C++ conversion result.

**Proposed fix:** Check charge before choosing the formula-interning branch.
If the destination cannot represent it, preserve the summed declared mass in
an anonymous mass-only record or report the unsupported conversion explicitly.
Do not validate one charged formula then silently store a neutral one. Test
nonzero, cancelling and zero summed charges with exact mass accounting.

**Rust handling:** The [native conversion group](docs/PROFORMA_CONVERSION_SUPPORT.md) preserves this source
consequence explicitly and includes a charged-formula regression. A scientific
correction remains a separately documented compatibility decision.

## CPP-038 — ProForma crosslink spectra count resolved linker mass twice

**Affected files:** [`src/openms/source/CHEMISTRY/ProForma.cpp`, lines 2820–2833](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2820), crosslink extraction and sequence conversion; modification attachment at 2305–2316; [`src/openms/source/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.cpp`, lines 920–925 and 967–971](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.cpp#L920), combined precursor mass.

**Issue and reproduction:** Directly construct two nonempty chains with matching
ordinary-residue crosslink labels. On the alpha endpoint, supply a resolved
modification with positive delta mass D and a first mass-delta alternative D
bearing the label. Give the beta endpoint only the matching label. For example,
use alpha AM linked at M and beta MA linked at M. Request precursor peaks with
`generateSpectrum(ion, 1, 1, "M", false, true)`.

`findCrossLink` supplies D as the separate linker mass. `toAASequence(BEST_EFFORT)`
also attaches the resolved endpoint modification to alpha, so its mass already
includes D. The XLMS backend then adds alpha mass, beta mass and linker mass D.
The precursor is heavier by D than the two unmodified chains plus one linker;
at charge z its m/z excess is D/z. Other linked fragments using this same
precursor mass are affected. The source mass-calculation API counts the shared
linker label once and does not introduce this extra copy.

**Evidence:** Direct review of extraction, unfiltered resolved-handle attachment
and both backend precursor-mass expressions. This is a public-AST trigger with
source-derived arithmetic; no executed C++ spectrum comparison or parser-produced
linker record is claimed.

**Proposed fix:** Build the XLMS input sequences with the selected linker
modification removed, retaining unrelated residue and terminal modifications,
then pass its mass once through `cross_linker_mass`. Resolve and validate both
endpoints before extracting that mass. Test the same linker represented on one
or both endpoints and check precursor and linked-fragment mass conservation.

**Rust handling:** The [native ProForma spectrum wrapper](docs/PROFORMA_SPECTRA_SUPPORT.md)
preserves this finite source behavior. Independent tests distinguish one-endpoint
2D and two-endpoint 3D contributions from the chemically expected single linker.
The standalone XLMS backend adds its supplied linker only once. No native or
upstream scientific correction is claimed.

## CPP-039 — SemanticValidator compares descendant units with the measured term

**Affected file:** [`src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp`, lines 326–348](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp#L326), descendant-unit fallback in `handleTerm_`.

**Issue and reproduction:** Enable unit checking and define a measurement term
whose allowed-unit set contains U_PARENT. Supply known unit U_CHILD, a descendant
of U_PARENT, on that measurement term. Keep the measurement outside the unit
hierarchy. Exact membership fails as intended, but the fallback callback compares
each descendant with `parsed_term.accession` (the measurement) instead of
`parsed_term.unit_accession` (U_CHILD). The valid descendant unit is rejected as
not allowed. The inverse confusion can also accept the wrong unit if a malformed
vocabulary places the measured term under an allowed unit.

**Evidence:** Direct review of the membership check, callback capture/comparison
and error branch. No C++ validation execution is claimed. This concerns the
optional unit check, disabled by default, and is distinct from CV XML writer
unit selection in CPP-021.

**Proposed fix:** Compare the descendant accession with
`parsed_term.unit_accession`. Test an exact allowed unit, an allowed child and
an unrelated known unit using a small independent vocabulary.

**Rust handling:** The [native SemanticValidator](docs/SEMANTIC_VALIDATOR_SUPPORT.md) compares descendants with the supplied unit accession. Dedicated exact/descendant/unrelated-unit cases cover the correction. Source-name checking remains separate.

## CPP-040 — Failed semantic validation contaminates later document paths

**Affected files:** [`src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp`, lines 96–102 and 111–121](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp#L96), per-call initialization and start callbacks; end callbacks at 141–224; inherited cleanup in [`src/openms/source/FORMAT/HANDLERS/XMLHandler.cpp`, lines 37–39](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/XMLHandler.cpp#L37).

**Issue and reproduction:** Configure a MUST rule at `/r/cvParam/@accession`,
then reuse one validator. First validate
`<r><cvParam accession="MS:1"/></r>`. The missing required name throws after both
tag names have been pushed. Catch the exception and validate `<r/>`.
Initialization clears only errors and warnings, while the open-tag stack retains
`r/cvParam`. The second document is checked beneath that stale prefix, so the
required `/r/cvParam/@accession` rule is missed and validation can return true.
A fresh validator correctly rejects the same second document for its missing
required term. Fulfilled-rule counters can likewise survive failed parsing.

**Evidence:** Direct review of callback order, limited initialization, successful
end-only stack/counter cleanup, and the empty inherited reset. No executed C++
reproduction is claimed. CPP-033 documents the analogous failure in the separate
mapping-file reader; this entry concerns semantic validation results.

**Proposed fix:** Use local per-document stack, fulfilled counts and diagnostics,
or clear all of them on entry and every exceptional exit. Publish results only
for the current document. Test failed parsing before and after fulfilled terms,
followed by valid and semantically invalid documents on the same validator.

**Rust handling:** The [native SemanticValidator](docs/SEMANTIC_VALIDATOR_SUPPORT.md) keeps its XML stack, rule counters and report local to each operation. Tests reuse one validator after missing-attribute and malformed-tail failures, verifying no contamination or partial public report.

## CPP-041 — The mzML header writer does not escape several string attributes

**Affected file:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, line 3765](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L3765), software version; checksum output at 3797/3801 and fraction identifier at 5226.

**Issue and reproduction:** Set a software version to `v"&1` or an experiment's
fraction identifier to `fraction "A&B"`, then write mzML. These ordinary string
fields are inserted directly inside double-quoted XML attributes. Quotes close
the attribute early and bare ampersands start invalid entity references, so the
result cannot be parsed as XML. The same unescaped interpolation is used for
caller-supplied checksum strings, although arbitrary text there may separately
be an invalid checksum.

**Evidence:** Direct review of the insertion expressions and comparison with
adjacent source name/path strings that use `writeXMLAttribute_`. No C++ execution
is claimed. This is distinct from the ControlledVocabulary rendering defect in
CPP-024: these output sites bypass that renderer.

**Proposed fix:** Apply `writeXMLAttribute_` to every externally supplied string
attribute, including version, fraction identifier and checksum text. Test quote,
ampersand, less-than and Unicode values through XML writing and reading; enforce
any checksum content restrictions separately from XML escaping.

**Rust handling:** The [native header writer](docs/MZML_HEADER_SUPPORT.md) uses the checked XML attribute encoder. Version and fraction quote/ampersand values round-trip through both writers and XSD validation. Invalid digest strings retain the separate checksum content guard. No upstream patch is claimed.

## CPP-042 — XLMS linear suffix losses divide the mass by charge twice

**Affected file:** [`src/openms/source/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.cpp`, lines 302–311](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.cpp#L302), linear suffix generation; loss helper at 560–599.

**Issue and reproduction:** Generate linear y ions for peptide KS with link
position 0, maximum charge 2 and neutral losses enabled. The suffix S permits
water loss. The suffix loop first computes m/z as `pos = mono_weight / charge`,
then passes `pos` to a helper whose argument is a charged mass. That helper
subtracts the neutral loss and divides by charge again. For charged mass M,
loss L and charge z, the emitted loss m/z is `(M/z - L)/z` instead of `(M-L)/z`.
The error affects linear x/y/z loss peaks at charges above one; the intact peak
and the prefix path do not share this extra division.

**Evidence:** Direct review of the caller argument and both water/ammonia helper
branches, with independent algebra. No C++ spectrum execution is claimed.

**Proposed fix:** Pass `mono_weight` rather than `pos` to
`addLinearIonLosses_` in the suffix branch. Test charge-one and charge-two suffix
losses alongside prefixes, asserting the intact-to-loss spacing is L/z.

**Rust handling:** The [native XLMS generator](docs/THEORETICAL_XLMS_SUPPORT.md)
retains the finite source expression, with a charge-two suffix-loss regression
and a separately calculated correct chemical expectation. The support document
makes this compatibility behavior explicit. No upstream correction is claimed.

## CPP-043 — XLMS precursor isotope companions omit charge normalization

**Affected file:** [`src/openms/source/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.cpp`, lines 609–624](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.cpp#L609), intact precursor and isotope companion; corresponding water/ammonia paths at 639–654 and 668–683.

**Issue and reproduction:** Enable precursor peaks and isotopes with
`max_isotope >= 2`, then generate a precursor of charge 2 through either XLMS
crosslink overload. The monoisotopic peak uses `mono_pos / charge`, where
`mono_pos` is neutral mass plus the proton contribution. Its isotope companion
uses `mono_pos + C13C12_MASSDIFF_U / charge`. Thus the charged mass has not been
divided for the companion. It appears near twice the precursor m/z at charge 2
instead of being separated by the isotope spacing divided by charge. The same
error affects isotope companions of the precursor water/ammonia losses.

**Evidence:** Direct review of all three pairs of assignments. The variable's
charged-mass meaning follows its immediately preceding initialization. No C++
spectrum execution is claimed.

**Proposed fix:** Use `(mono_pos + C13C12_MASSDIFF_U) / charge` in all three
companion branches, or derive each from its already normalized monoisotopic m/z.
Test the companion spacing at charges 1, 2 and 3 for intact and both loss peaks.

**Rust handling:** The [native XLMS generator](docs/THEORETICAL_XLMS_SUPPORT.md)
retains these finite source values. Tests distinguish all three source
expressions from the independently calculated isotope spacing. No native or
upstream scientific correction is claimed.

## CPP-044 — Unmapped CV lookup depends on earlier validation calls

**Affected file:** [`src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp`, line 147](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp#L147), insertion of empty rule lists during closing callbacks; analogous insertion at 279; const lookup at 538.

**Issue and reproduction:** Construct a validator with no mapping for
`/r/cvParam/@accession`. Calling `locateTerm` for this path throws
`std::out_of_range` because it uses `rules_.at(path)`. Validate `<r/>`, then call
the same query with the same term. Closing r inserted an empty rule list through
`rules_[path]`, so the second query returns false. Neither the mapping nor the
vocabulary changed, but an unrelated validation call changed the query's missing-
path behavior. Validation also retains these unused path entries across files.

**Evidence:** Direct review of mutable callback lookups and the const public
predicate. No C++ execution is claimed. This is independent of failed-parser
state contamination in CPP-040 and of the incorrect software path in CPP-028.

**Proposed fix:** Avoid inserting rules while performing lookups, and define one
missing-path contract for `locateTerm`. Returning false is natural for its
allowed-term predicate; a documented checked exception can also be consistent.
Test the same unmapped query before and after validation of different documents.

**Rust handling:** The [native validator](docs/SEMANTIC_VALIDATOR_SUPPORT.md) builds a fresh bounded mapping index per call. An absent locate_term path consistently returns a checked error before or after unrelated validation. This deliberately removes incidental source cache state.

## CPP-045 — PeptideEvidence misclassifies valid and invalid position limits

**Affected file:** [`src/openms/source/METADATA/PeptideEvidence.cpp`, lines 79–84](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/METADATA/PeptideEvidence.cpp#L79), `hasValidLimits`; zero-based inclusive endpoint contract in the corresponding header.

**Issue and reproduction:** Evidence with start 0 and end 0 describes a legitimate
one-residue peptide at the protein's N terminus. `hasValidLimits()` nevertheless
returns false because it rejects any end equal to `N_TERMINAL_POSITION` (0).
Conversely, start 5/end 2 returns true, as do negative coordinates other than
the exact unknown sentinel -1. These do not describe a valid inclusive interval.

**Evidence:** Direct source review of the constants and three-condition predicate.
The existing native identification documentation and test already cover the
first-residue correction. No C++ execution is claimed.

**Proposed fix:** Accept exactly known nonnegative endpoints with `start <= end`.
Keep unknown positions distinct from an actual zero coordinate. Test 0..=0,
a normal interval, unknown endpoints, reversed limits and other negative values.

**Rust handling:** The existing `Option<usize>` model distinguishes missing from
zero and excludes negative coordinates. `has_valid_limits()` accepts ordered
known intervals, including 0..=0; `validate()` rejects reversed endpoints. These
corrections are retained and explicitly exercised by the identification tests.

## Maintaining this log

Add an entry whenever porting or testing identifies a new original C++ defect.
Include the checked source revision, affected files/functions/lines, a concrete
trigger, observed versus expected behavior, evidence level, proposed C++ fix,
regression coverage and the Rust port's handling. Keep stable IDs and update
status only after verifying the upstream fix. Retain resolved entries with the
fixing commit. Link durable probe data once it lands; temporary investigation
paths are not sufficient long-term evidence. Do not promote a suspicion to an
executed or confirmed finding without verification.

## CPP-046 — Semantic date validation applies date-time syntax to xsd:date

**Affected files:** [`src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp`, lines 502–516](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp#L502), XSD_DATE value branch; [`src/openms/source/DATASTRUCTURES/DateTime.cpp`, lines 188–335](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/DATASTRUCTURES/DateTime.cpp#L188), general date-time parser.

**Issue and reproduction:** Define an allowed CV term with `xref: value-type:xsd\:date`
and validate its value `2001-02-03`. The validator calls `DateTime::set`, whose
plain ISO branch requires six date/time fields. It rejects this valid date and
reports a wrong xsd:date value. Conversely, `2001-02-03T04:05:06` passes although
it is a date-time value. [The W3C xsd:date lexical contract](https://www.w3.org/TR/xmlschema11-2/#date)
consists of year, month and day, followed by an optional timezone, without a time.
The legacy parser's acceptance of a trailing `Z` does not repair the ordinary
no-timezone case.

**Evidence:** Pinned-source inspection of the validator and every DateTime parse
branch, compared with the primary XSD specification. Native SemanticValidator
regressions preserve both the date-only rejection and timestamp acceptance.
No complete C++ SemanticValidator execution is claimed.

**Proposed fix:** Validate this branch with a dedicated xsd:date lexical/calendar
check, including optional timezone bounds and the applicable year rules. Keep
xsd:dateTime separate if the vocabulary declares it. Cover date-only, zoned dates,
invalid calendar dates and timestamp rejection; changing the general DateTime
parser alone would still accept invalid date-time lexemes here.

**Rust handling:** The [native SemanticValidator](docs/SEMANTIC_VALIDATOR_SUPPORT.md) retains the source value conversion for compatibility, with explicit date-only rejection and timestamp-acceptance tests. Its value checks do not claim complete XSD lexical conformance. No upstream patch has been applied.

## CPP-047 — ProForma XLMS link positions ignore preceding flattened sections

**Affected file:** [`src/openms/source/CHEMISTRY/ProForma.cpp`, lines 2017–2041](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2017), `findCrossLink`; flattened sequence construction at 2280–2288; two-chain spectrum construction at 2820–2835.

**Issue and reproduction:** Construct an alpha Peptidoform whose first section
is an unmodified ModifiedRange containing A and G, followed by a SequenceElement
M with the first alternative carrying crosslink label XL1. Give beta an M with
the matching label followed by A, with `is_chimeric=false`. Label-only endpoints
and zero linker mass isolate this from CPP-038. The two-chain spectrum check
accepts the matching labels. BEST_EFFORT conversion emits alpha sequence AGM,
but `findCrossLink` increments its position only for top-level SequenceElement
sections and reports alpha position 0, rather than the M position 2. XLMS then
uses incorrect fragment boundaries and K-linked eligibility. A preceding
nonempty AmbiguousRegion has the analogous one-residue discrepancy.

**Evidence:** Direct source review of the position counter, two-chain issue
checks, BEST_EFFORT flattening and crosslink construction. This is an explicit
public AST trigger, not a claimed parsed-text or executed C++ spectrum result.

**Proposed fix:** Share flattened output-position accounting between conversion
and link extraction: advance by range length and by one for nonempty ambiguous
regions. Define link selection inside those sections consistently as well, or
reject unsupported section shapes before generation. Test a range and an
ambiguous region preceding each chain's endpoint, independently of linker mass.

**Rust handling:** The [native ProForma spectrum wrapper](docs/PROFORMA_SPECTRA_SUPPORT.md)
preserves this finite source behavior with explicit range/ambiguity regressions
and separate correct flattened-position expectations. Checked backend bounds
still apply. No upstream fix has been applied.

## CPP-048 — Zero centroid-inspection limit depends on the first spectrum type

**Affected files:** [`src/openms/source/FORMAT/MzMLFile.cpp`, lines 233–268](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/MzMLFile.cpp#L233), `getCentroidInfo`; public limit contract in `src/openms/include/OpenMS/FORMAT/MzMLFile.h`, lines 180–193; unsigned `Size` definition in `src/openms/include/OpenMS/CONCEPT/Types.h`, line 97.

**Issue and reproduction:** Call `getCentroidInfo(filename, 0)` on a valid file.
The remaining counter starts at zero. If the first accepted spectrum is
centroided or profile, its branch decrements the unsigned counter to `SIZE_MAX`,
and the zero-stop test fails. Parsing normally continues through the file,
counting its recognized spectra despite the requested zero limit. If the first
spectrum remains unknown, its branch does not decrement the counter; the same
stop test instead ends parsing after that one unknown spectrum. There is no
documented special zero mode or positive-limit precondition. This inconsistent
limit behavior is separate from the option-restoration defect in CPP-018.

**Evidence:** Direct source control-flow review and the exact `typedef size_t
Size` declaration establish defined unsigned wraparound. The public comment
says only the requested number of non-unknown spectra is inspected. This is
source evidence, not an executed C++ reproduction or inferred native test result.

**Proposed fix:** Handle zero before constructing the consumer or parsing input:
either return an empty result to implement a zero inspection count, or reject
zero with a documented positive-limit precondition. Guard decrement/termination
so it cannot wrap. Test zero with both known-first and unknown-first files,
then positive limits with unknown spectra interleaved between recognized ones.

**Rust handling:** [Native centroid inspection](docs/MZML_CENTROID_SUPPORT.md)
rejects a zero limit with `InvalidValue` before even requesting the path. Tests
cover this boundary and retain source count/stop order for positive quotas,
including Unknown interleaving and multiple MS levels. No upstream fix is claimed.

## CPP-049 — Indexed mzML output writes a constant placeholder checksum

**Affected files:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp`, lines 127–134](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp#L127), indexed footer output; checksum contract in `share/OpenMS/SCHEMAS/mzML_idx_1_10.xsd`, lines 1185–1189.

**Issue and reproduction:** Write an indexed mzML document through the source
writer. The footer unconditionally assigns `sha1_checksum = "0"` and emits
`<fileChecksum>0</fileChecksum>`. It never calculates the digest of the preceding
file bytes. Consequently every indexed output has the same placeholder,
irrespective of its spectrum data, metadata or byte offsets. The pinned format
schema documents a SHA-1 checksum covering the file start through the end of
the opening `fileChecksum` tag. Its `xs:string` type does not enforce this
semantic checksum requirement, so schema validation alone cannot detect the
placeholder.

**Evidence:** Direct review of the literal assignment, adjacent TODO specifying
the digest boundary, and the retained schema annotation. No C++ writer execution
or comparison against an executed digest is claimed.

**Proposed fix:** Update an incremental SHA-1 state with the exact serialized
bytes, including the complete opening checksum tag. Emit the resulting digest
without feeding the digest text into its own calculation. Test a minimal file
and varied spectrum/header data against an independent digest over the emitted
prefix; include byte-offset and line-ending cases.

**Rust handling:** The [native indexed writer](docs/MZML_WRITE_OPTIONS_SUPPORT.md)
calculates the real checksum from successfully written bytes. Independent Python
hashlib tests cover UTF-8, escaping, partial writes and SHA-1 block boundaries.
No upstream fix has been applied.

## CPP-050 — Empty indexed mzML declares zero indices but emits a dummy index

**Affected files:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp`, lines 80–123](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp#L80), `writeFooter_`; index-count and offset contracts in `share/OpenMS/SCHEMAS/mzML_idx_1_10.xsd`, lines 1119–1154.

**Issue and reproduction:** Request indexed output with no spectra or
chromatograms. Both offset vectors are empty, so `indexlists` is zero and the
writer emits `<indexList count="0">`. It then emits one `index` named `dummy`
containing an offset of -1 and an `idRef` of `dummy`. The declared count disagrees
with its actual child count, and the fabricated offset cannot address a real
indexed record. This is separate from the constant checksum in CPP-049.

**Evidence:** Direct review of the empty-vector count expression and the
`indexlists == 0` branch. The pinned schema describes `count` as the number of
indices and offsets as pointers to identified elements; its index and offset
elements each require at least one occurrence. The count relation is not an XSD
constraint, so no actual schema-validation failure is asserted here. No C++
writer execution is claimed.

**Proposed fix:** Define a consistent empty-output policy before writing. When
the selected indexed schema cannot represent zero indexed records, reject that
combination with a clear error or explicitly select ordinary mzML. Do not
silently emit a nonexistent index target. Derive every declared index count
from the indices actually written. Test empty, spectrum-only, chromatogram-only
and mixed experiments with independent count and byte-target checks.

**Rust handling:** The [native indexed writer](docs/MZML_WRITE_OPTIONS_SUPPORT.md)
rejects an empty experiment before requesting a path or writing external bytes.
Tests verify unchanged output; explicit ordinary mzML still accepts empty input.
No upstream fix has been applied.

## CPP-051 — mzML semantic validation reuses parameter groups from earlier documents

**Affected files:** [`src/openms/source/FORMAT/VALIDATORS/MzMLValidator.cpp`, lines 47–84](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/VALIDATORS/MzMLValidator.cpp#L47), group definition/reference callbacks; `src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp`, lines 86–108, inherited `validate` initialization.

**Issue and reproduction:** Reuse one `MzMLValidator` for two documents. The
first defines a referenceable parameter group `g` containing an allowed term
required by a mapping rule. The second references `g` without defining it.
The inherited `validate` clears diagnostics but does not clear `param_groups_`.
The reused validator therefore applies the earlier document's term, potentially
accepting a document that a fresh validator reports as missing a required term.
This contamination occurs after successful validation, independently of the
exception-state problem in CPP-040. Repeated definitions can also accumulate
terms from previous documents.

**Evidence:** Direct review of the persistent member, definition and reference
callbacks, and inherited initialization. There is no per-document group reset.
The two-document scenario is derived from that control flow; no executed C++
reproduction or full SDK build is claimed.

**Proposed fix:** Keep parameter groups, the current group ID and binary-array
state local to each validation operation, or clear them before parsing with
exception-safe cleanup. Preserve repeated-group behavior within a document.
Test reuse after successful and failed parses against a fresh validator,
including missing references and repeated IDs.

**Rust handling:** The [native mzML validator](docs/MZML_VALIDATOR_SUPPORT.md)
uses operation-local group, binary and rule state. Direct two-document tests
verify fresh behavior after successful and failed parses, while within-document
duplicate and forward-reference semantics remain source-compatible. No upstream
fix has been applied.

## CPP-052 — mzML schema selection depends on the first four physical lines

**Affected files:** [`src/openms/source/FORMAT/MzMLFile.cpp`, lines 58–80](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/MzMLFile.cpp#L58), `isValid`; `src/openms/source/FORMAT/TextFile.cpp`, lines 24–81, and `src/openms/include/OpenMS/FORMAT/TextFile.h`, line 54, four-line loading with empty/comment skipping disabled.

**Issue and reproduction:** Insert three XML comment lines after the XML
declaration of an indexed mzML document. Its `indexedmzML` root now starts on
line five. `isValid` concatenates only the first four trimmed lines and searches
for the literal `<indexedmzML`, so it selects the ordinary mzML schema. That
schema has no declaration for the indexed root. Schema choice therefore
depends on an otherwise legal XML prolog rather than the document element.

**Evidence:** [Recorded probe](docs/mzml-schema-selection-probe.json) retains
source hashes, the exact transformation and actual `xmllint` results. Both the
original `MzMLFile_4_indexed.mzML` and the transformed document pass the pinned
indexed XSD. The transformed document fails the ordinary XSD. The four-line
predicate is false by direct source review and independent byte inspection;
the C++ method was not executed. The inserted comments leave index offsets and
checksum untouched, so this probe asserts XSD validity only, not indexed-file
integrity. No full SDK build is claimed.

**Proposed fix:** Select the schema from the parsed document element's expanded
name, allowing legal prolog whitespace, comments and namespace prefixes. Keep
XML parsing bounded, then validate against the matching schema. Test a root
beyond line four and equivalent prefixed/default-namespace forms.

**Rust handling:** The optional [XSD validator](docs/MZML_SCHEMA_SUPPORT.md)
selects the schema from the checked document root and namespace. Direct tests
cover delayed and prefixed roots, UTF-16 input and malformed namespaces. No
upstream fix is claimed.

## CPP-053 — Chromatogram array promotion mismatches values and their physical role

**Affected files:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, lines 634–650](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L634), chromatogram population; `src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp`, lines 268–281, `computeDataProperties_`; `MzMLHandler.cpp`, lines 5614–5619 and 5693–5695, output role selection.

**Issue and reproduction:** Supply a chromatogram with equally sized time,
pressure and flow-rate arrays, in that order, with one pressure value of 200
and one flow value of 3. Population renames both physical arrays to `intensity
array`, while repeatedly setting the record's `mzml intensity array` metadata.
The final selector is `flow`, but `computeDataProperties_` returns the first
matching intensity array, so the resulting peak contains 200. Writing uses the
stored selector and emits that pressure value as flow, with flow units. A
canonical intensity array followed by a pressure array has the same mismatch.
The other promoted arrays are also excluded from supplemental-array creation
because they now carry the primary name.

**Evidence:** Direct review of the promotion loop, first-match lookup, peak
assignment, supplemental-array exclusion and writer's selector lookup. The
values above are an independent control-flow example, not an executed C++
reproduction. No full C++ SDK build or upstream fix is claimed.

**Proposed fix:** Select the primary array deterministically before assigning its
role. Set the selector from that selected array only, and preserve other arrays
under their actual names and physical units. Alternatively, reject ambiguous
primary candidates explicitly. Test both input orders and a canonical intensity
array combined with pressure, flow or detector-signal arrays, checking values,
roles and units after a write/read cycle.

**Rust handling:** The native reader rejects competing primary intensity
candidates before publishing a record. Writer preflight also rejects auxiliary
names which would become competing primaries on reload. Direct typed transport
tests cover pressure/flow ambiguity, canonical intensity competition and all
writer families; see [typed transport](docs/MZML_TYPED_TRANSPORT_SUPPORT.md).

## CPP-054 — Shipped indexed mzML schema leaves index references unchecked

**Affected files:** [`share/OpenMS/SCHEMAS/mzML_idx_1_10.xsd`, lines 1193–1200](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/share/OpenMS/SCHEMAS/mzML_idx_1_10.xsd#L1193), `KEY_ID_IDX` and `FKNID`; lines 1147–1155, `OffsetType`.

**Issue and reproduction:** Change the original indexed mzML fixture's offset
reference from `index=19` to the nonexistent `does_not_exist`. Both the record
key and offset-reference selectors start with `.//dx:indexedmzML/...` inside
the declaration for that root. They therefore look for another indexed root
below the current root and select no normal records or offsets. The reference
field also names `@id`, although offsets use `@idRef`. A dangling index
reference passes this shipped schema.

**Evidence:** The [recorded probe](docs/mzml-index-schema-probe.json) preserves
source hashes, the exact single replacement and independent `xmllint` results:
the original and changed documents both pass the unchanged pinned XSD. No
C++ method execution is claimed. The transformation does not repair byte offsets
or checksum; this probe demonstrates missing ID-reference validation only.

**Proposed fix:** Make the key selectors relative to the declared root:
`dx:mzML/dx:run/dx:spectrumList/dx:spectrum` and the equivalent chromatogram
path. Select `dx:indexList/dx:index/dx:offset` for the key reference and use
`@idRef` as its field. Test valid references, missing targets, duplicate record
IDs and both indexed record kinds. Byte offsets and checksum still require
separate integrity checks beyond XSD.

**Rust handling:** The optional [XSD operation](docs/MZML_SCHEMA_SUPPORT.md)
retains the exact original schema bytes and explicitly tests this limitation.
Separate genuine key/keyref failures verify that the engine enforces other
identity constraints. The indexed writer's independent tests check actual
ID/offset correspondence and checksum. The upstream schema has not been patched.

## CPP-055 — Base64 SIMD decoding does not validate its alphabet

**Affected files:** [`src/openms/source/FORMAT/Base64.cpp`, lines 83–110 and 166–208](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/Base64.cpp#L83), `registerDecoder_` and `stringSimdDecoder_`; `src/openms/include/OpenMS/FORMAT/Base64.h`, lines 310–350, numeric uncompressed decoding; `src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp`, lines 138–195, mzML numeric decoding.

**Issue and reproduction:** Pass 16 exclamation marks to the uncompressed
numeric Base64 decoder. The length check accepts this multiple of four. The
SIMD decoder then classifies each exclamation mark as a number using only an
upper-bound comparison against `'9' + 1`; it has no lower-bound or invalid-byte
check. It produces bytes that the caller interprets as numeric data instead of
reporting malformed Base64. mzML's optional whitespace removal cannot correct
these non-whitespace invalid characters.

**Evidence:** Direct review of the public numeric decode path, SIMD masks and
output conversion. No executed C++ reproduction or particular resulting numeric
value is asserted. The trigger is invalid Base64 independently of the selected
floating-point type and byte order.

**Proposed fix:** Validate alphabet membership, padding position and decoded
length before converting SIMD output into numeric values. Retain an explicit
whitespace-normalization policy, but reject invalid remaining characters. Test
punctuation, misplaced padding and malformed tails in ordinary and compressed
numeric decoding; no partially decoded array should escape on error.

**Rust handling:** The mzML transport uses checked Base64 decoding and rejects
invalid alphabet bytes. [The normalization option](docs/MZML_NORMALIZATION_SUPPORT.md)
executes `skip_xml_checks` without disabling decoding checks or XML legality.
Direct tests cover malformed alphabets, padding, XML, every supported codec and
both option settings. No upstream fix has been applied.

## CPP-056 — Synthetic mzML scans omit spectrum ion mobility

**Affected files:** [`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, lines 5393–5435](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L5393), `writeSpectrum_`; lines 5468–5500, the fallback for empty acquisition information.

**Issue and reproduction:** Construct a spectrum with drift time 1.5 and
`DriftTimeUnit::MILLISECOND`, no precursors and empty acquisition information.
The writer emits mobility only inside the first iteration of the acquisition
loop. With no acquisition records, it creates a synthetic scan containing RT,
zoom state and scan windows, but no mobility CV term. The otherwise represented
spectrum drift time and unit are absent from the file. The same omission affects
inverse reduced mobility, collision cross section and signed FAIMS voltage.

**Evidence:** Direct review of the complete acquisition loop, its first-scan
mobility switch and the synthetic-scan branch at the pinned revision. The input
above isolates spectrum mobility from selected-ion mobility. No C++ execution,
full SDK build or upstream fix is claimed.

**Proposed fix:** Share first-scan RT/mobility emission between real and synthetic
scans. Emit the same supported mobility term and unit in both branches, keeping
signed FAIMS values and the unset-value rules. Test all four units with zero,
one and multiple acquisition records, checking that exactly the first scan
carries the spectrum mobility and that a store/load cycle retains it.

**Rust handling:** Spectrum mobility transport is being implemented separately.
The native writer will use the same first-scan path for explicit and synthetic
acquisitions. No implemented native correction is claimed at this checkpoint.

## CPP-057 — Spectrum equality ignores ion-mobility format and peak type

**Affected files:** [`src/openms/source/METADATA/SpectrumSettings.cpp`, lines 23–43](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/METADATA/SpectrumSettings.cpp#L23), equality and inequality; `src/openms/include/OpenMS/METADATA/SpectrumSettings.h`, lines 170–171, mobility fields; `src/openms/source/KERNEL/MSSpectrum.cpp`, lines 514–529, spectrum equality.

**Issue and reproduction:** Create two otherwise identical `SpectrumSettings`
values and set their IM peak types to `IM_CENTROIDED` and `IM_PROFILE`.
`operator==` returns true because it compares neither `im_peak_type_` nor
`im_type_`. Changing only the IM format is also invisible to equality.
`MSSpectrum::operator==` delegates to this comparison and does not compare these
fields separately, so spectra with different IM representation state also
compare equal. The header declares ordinary equality; unlike the documented
name/range exclusions in spectrum equality, these omissions are not documented
as intentional.

**Evidence:** Direct review of the complete comparisons, member defaults and
the public setters/getters at `SpectrumSettings.cpp:80–97`. This is a source
control-flow finding, not an executed C++ reproduction.

**Proposed fix:** Compare both `im_type_` and `im_peak_type_` in
`SpectrumSettings::operator==`. Keep inequality as its negation. Test settings
and spectra which differ in each field independently, plus equal copies and
ordinary spectrum-type differences.

**Rust handling:** The existing native `SpectrumSettings` derives value equality
including both mobility fields. The spectrum mobility port will also include
its represented fields in record equality and add direct regressions. No
upstream patch or completed spectrum mobility implementation is claimed here.

## CPP-058 — Nonnumeric DataValue casts read an inactive union member

**Affected files:** [`src/openms/source/DATASTRUCTURES/DataValue.cpp`, lines 452–491](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/DATASTRUCTURES/DataValue.cpp#L452), conversion to long double, double and float; `src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h`, lines 159–184 and 409–417, conversion contract and tagged union; `src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, lines 1400–1403, elution-time fallback.

**Issue and reproduction:** Construct `DataValue(std::string("42"))` and request
a floating conversion, such as `static_cast<double>(value)`. Each floating
conversion rejects only Empty and handles Integer separately; every other
alternative reads `data_.dou_`. String and list alternatives store pointers in
other union members. Reading the inactive double member is undefined C++
behavior, violating the documented `ConversionError` contract for a wrong type.
This is not a numeric-string parser. No specific resulting number is guaranteed.

The same public cast is reached when a spectrum without scan RT contains a
String-valued `elution time (seconds)` user parameter: the mzML close handler
passes that DataValue directly to `setRT`.

**Evidence:** Direct review of all three conversions, the tagged-union members,
the public exception documentation and the mzML caller. No C++ execution,
sanitizer result or predicted pointer-dependent numeric output is claimed.

**Proposed fix:** Switch explicitly on the stored type: convert Integer, return
the floating member for Double, and throw `ConversionError` for every other
alternative. Apply the same rule to all three floating conversions. Test Empty,
String and all list alternatives, positive/negative numeric values and the mzML
fallback with a nonnumeric metadata value.

**Rust handling:** `MetaValue` accessors match the stored variant and cannot read
inactive storage. The typed mzML elution-time fallback accepts Integer/Float
variants and rejects String/list/Empty values before publishing a record. Direct
tests cover consumed nonnumeric values, explicit scan RT precedence and the
source ordering before primary-array metadata merges. No upstream fix is claimed.

## CPP-059 — Short experimental-design rows are read past the end of the row vector

**Affected files:** [`src/openms/source/FORMAT/ExperimentalDesignFile.cpp`, lines 226–237](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/ExperimentalDesignFile.cpp#L226), one-table content rows; the same file, lines 411–419, two-table sample rows; [`src/openms/source/METADATA/ExperimentalDesign.cpp`, lines 940–962](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/METADATA/ExperimentalDesign.cpp#L940), both `SampleSection::getFactorValue` overloads.

**Issue and reproduction:** Load a one-table design whose header declares five
columns and whose first data row has three cells, for example
`Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample` followed by
`1\t1\ta.mzML`. `parseOneTableFile_` reads
`cells[fs_column_header_to_index["Label"]]`, `["Fraction"]` and
`["Fraction_Group"]` with `std::vector::operator[]` before it reaches the
`parseErrorIf_(n_col != cells.size(), ...)` guard, so an out-of-range index is
dereferenced rather than reported. The check is present but runs three reads too
late.

The two-table sample section has no such check at all. `SAMPLE_HEADER` assigns
`n_col = sample_columnname_to_columnindex_.size()` and `SAMPLE_CONTENT` never
compares it, so a sample row shorter than its header reads
`cells[sample_columnname_to_columnindex_["Sample"]]` out of range and then stores
the short row in `content_`. Every later `getFactorValue` compounds it: both
overloads bound-check the row with `content_.at()` and then index the row with
`sample_row[col_index]` unchecked, so a stored short row is read out of range on
each access. All three sites are undefined behavior, not a diagnosable parse
error. The two-table sample section is the common shape, because the file section
is checked and the sample section is not.

**Evidence:** Direct review of both parser state machines, the placement of the
one guard relative to the reads it protects, the absent sample-row guard, and
both `getFactorValue` overloads. No C++ execution or sanitizer result is claimed.

**Proposed fix:** Move the `n_col != cells.size()` check in
`parseOneTableFile_` above the first cell read, add the same check to the
`SAMPLE_CONTENT` branch of `parseTwoTableFile_`, and replace
`sample_row[col_index]` with `sample_row.at(col_index)` in both `getFactorValue`
overloads so a section constructed directly through the public
`SampleSection(content, ...)` constructor cannot get past it either. Test a short
data row and a short sample row in both layouts, and a directly constructed
section whose content rows are shorter than its column map.

**Rust handling:** `format::experimental_design_file` checks cell counts before
any indexed access in both layouts, and `SampleSection::from_table` rejects a row
shorter than its column map when the section is built. `factor_value` and
`factor_value_by_row` return a typed error instead of indexing. Tests cover a
short one-table data row, a short two-table sample row and a directly built
short-row section. No upstream fix is claimed.

## CPP-060 — Negative design indices wrap to large unsigned values

**Affected files:** [`src/openms/source/FORMAT/ExperimentalDesignFile.cpp`, lines 226–256](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/ExperimentalDesignFile.cpp#L226), one-table `Label`/`Fraction`/`Fraction_Group`; the same file, lines 386–392, the two-table equivalents.

**Issue and reproduction:** Load a design containing `Fraction` `-1`. Both
parsers read the cell with `StringUtils::toInt32`, which accepts the negative
value, and assign it to the `unsigned` members of `MSFileSectionEntry`, so the
row is stored with fraction 4294967295. Nothing downstream rejects it:
`getNumberOfFractions` counts it as a distinct fraction, `isFractionated`
reports the design as fractionated, and `getPathLabelToFractionMapping` publishes
it. A negative `Label` is worse, because the one-table guard
`parseErrorIf_(!has_sample && (label > 1), ...)` tests the signed value: `-1` is
not greater than 1, so a design without a `Sample` column passes the multiplex
check and then stores label 4294967295. A negative `Fraction_Group` is the only
one caught, and only indirectly, by the `isValid_` consecutive-from-1 rule.

**Evidence:** Direct review of the conversion and assignment in both parsers, the
member types in `ExperimentalDesign.h` lines 468–482, and the guard that tests
the signed value before the assignment. No C++ execution is claimed.

**Proposed fix:** Reject a negative value in the parsers, before the assignment
and before the `label > 1` guard, with the existing `parseErrorIf_` diagnostic;
or read the three columns as unsigned. Test `-1` in each of the three columns in
both layouts, including a one-table design with no `Sample` column.

**Rust handling:** The parser converts each index through `i32` and then a
checked `u32` conversion, so a negative `Fraction`, `Fraction_Group` or `Label`
is a typed parse error naming the column and line. Tests cover a negative
fraction. No upstream fix is claimed.

## CPP-061 — SpectrumHelper::makePeakPositionUnique discards the whole spectrum record

**Affected files:** `src/openms/include/OpenMS/KERNEL/SpectrumHelper.h`, line 196.

**Issue and reproduction:** `makePeakPositionUnique` ends with `std::swap(p_new, p)` where `p_new` is default constructed and only ever received merged peaks. The result therefore loses retention time, MS level, name, native identifier, precursors, instrument settings and metadata. The `OPENMS_LOG_WARN` at line 150 announces only that the data arrays are dropped, so a caller is told about a smaller loss than actually occurs. Trigger: any non-empty spectrum with a retention time set.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Copy `p` into `p_new` (or call `copySpectrumMeta`) before filling it, then clear the data arrays. Test a spectrum carrying RT, MS level, name and metadata through the call.

**Rust handling:** `make_peak_position_unique` keeps every record field by default; `UniquePositionOptions::reset_metadata` selects the source behaviour explicitly.

## CPP-062 — SpectrumRangeManager::byMSLevel(0) can only throw

**Affected files:** `src/openms/include/OpenMS/KERNEL/SpectrumRangeManager.h`, lines 82-101 and 124-127.

**Issue and reproduction:** `extend` and `extendUnsafe` document `ms_level = 0` as addressing the global ranges, and `byMSLevel` declares `UInt ms_level = 0` as its default argument. The global ranges live in the base subobject and are never inserted into `ms_level_ranges_`, so calling `byMSLevel()` with its own default always throws. The class test codifies the throw.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Return the base subobject for level 0, or remove the misleading default argument. Test `byMSLevel()` with no argument on a populated manager.

**Rust handling:** `by_ms_level(0)` returns `None` and `global()` reads the base; both are stated in the API table.

## CPP-063 — MSExperiment::updateRanges registers no per-level entry for MS level 0

**Affected files:** `src/openms/source/KERNEL/MSExperiment.cpp`, lines 698-699, with `SpectrumRangeManager.h` line 84.

**Issue and reproduction:** A spectrum at MS level 0 extends the global ranges twice and never creates a per-level entry, so `SpectrumRangeManager::getMSLevels()` omits level 0 although `MSExperiment::getMSLevels()` reports it. The two level lists disagree for the same run.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Insert a per-level entry for level 0 as for any other level, or document the exclusion in both accessors. Test a run holding a level-0 spectrum.

**Rust handling:** Replicated, with the disagreement recorded in `docs/RANGES_SUPPORT.md`.

## CPP-064 — RangeBase accepts NaN and infinity and breaks its own emptiness invariant

**Affected files:** `src/openms/include/OpenMS/KERNEL/RangeManager.h`, lines 113-124, 157-171 and 235-251.

**Issue and reproduction:** `setMin`, `setMax` and `extend` perform no finiteness check. `setMin(NaN)` leaves a range for which the documented `isEmpty()` equivalence `min > max` is neither true nor false, and `extend(NaN)` is a silent no-op because every comparison against NaN is false. Later `contains`, `clampTo` and `pushInto` calls then behave arbitrarily.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Reject non-finite input, or document the range as undefined once a non-finite value is stored. Test `setMin`, `setMax` and `extend` with NaN and both infinities.

**Rust handling:** Non-finite input is rejected with `Error::InvalidRange` before any state changes.

## CPP-065 — RangeUtils energy and isolation predicates contradict their own notes for MS1 spectra

**Affected files:** `src/openms/include/OpenMS/KERNEL/RangeUtils.h`, notes at lines 517-518, 570 and 615 against the bodies at 542, 557, 595 and 640.

**Issue and reproduction:** The notes for `IsInCollisionEnergyRange`, `IsInIsolationWindowSizeRange` and `IsInIsolationWindow` state that MS1 spectra, and spectra with no collision energy, return true. The bodies return false in those cases regardless of the `reverse` flag. The code is what the `remove_if` callers need; the documentation describes the opposite filter.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Correct the three notes to state that the predicate returns false, and that `reverse` does not invert the early return.

**Rust handling:** The code behaviour is preserved; each predicate's rustdoc states the early return and the mismatch.

## CPP-066 — BinnedSpectrum default construction leaves a null bin matrix

**Affected files:** `src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h`, line 93; `src/openms/source/KERNEL/BinnedSpectrum.cpp`, line 191.

**Issue and reproduction:** The default constructor leaves `bins_` as a null pointer. `getBinIntensity()` dereferences it without a check, so a default-constructed object is unusable and crashes rather than reporting the error.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Either allocate an empty matrix in the default constructor or check the pointer in every accessor. Test `getBinIntensity` on a default-constructed object.

**Rust handling:** A bin-less object cannot be constructed; the default constructor is deliberately not ported.

## CPP-067 — BinnedSpectrum::getBinIntensity mutates the spectrum it reads

**Affected files:** `src/openms/source/KERNEL/BinnedSpectrum.cpp`, line 191.

**Issue and reproduction:** The method is non-const and uses Eigen's `coeffRef`, which inserts an explicit zero coefficient for every miss. Reading an absent bin therefore changes `nonZeros()` and can change the result of `operator==`, so two equal spectra stop comparing equal after one of them is read.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Use a read-only coefficient accessor and make the method const. Test equality before and after reading an absent bin.

**Rust handling:** Bin lookup is read-only and takes `&self`.

## CPP-068 — BinnedSpectrum::operator== ignores the bin offset

**Affected files:** `src/openms/source/KERNEL/BinnedSpectrum.cpp`, line 125.

**Issue and reproduction:** Equality compares the bin size, the ppm flag and the bin contents but not `offset_`. Two spectra binned on different grid origins compare equal while their bins mean different m/z intervals.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Include `offset_` in the comparison. Test two spectra that differ only in offset.

**Rust handling:** The derived `PartialEq` includes the offset; the difference is stated in `docs/COMPARISON_SUPPORT.md`.

## CPP-069 — BinnedSpectrum left-boundary guard relies on unsigned wraparound

**Affected files:** `src/openms/source/KERNEL/BinnedSpectrum.cpp`, line 88.

**Issue and reproduction:** The guard `static_cast<int>(idx - j - 1) >= 0` computes an unsigned difference that wraps, then truncates it to `int`. For `idx >= 2^31` the truncation changes sign and the guard admits an out-of-range index.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Compare the indices before subtracting, or use a signed type throughout.

**Rust handling:** Uses `saturating_sub`, so the boundary cannot wrap.

## CPP-070 — FeatureHandle::asMutable casts away constness of a possibly const object

**Affected files:** `src/openms/include/OpenMS/KERNEL/FeatureHandle.h`, line 148.

**Issue and reproduction:** `asMutable` applies `const_cast` to `*this` and the header carries a TODO acknowledging it. When the referent is genuinely const, writing through the result is undefined behaviour.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Provide a non-const overload instead of casting, so the compiler rejects mutation of a const handle.

**Rust handling:** Replaced by ordinary `&mut` access and `ConsensusFeature::set_handles`; no cast exists.

## CPP-071 — DRange default constructor contradicts its documentation

**Affected files:** `src/openms/include/OpenMS/KERNEL/../DATASTRUCTURES/DRange.h`, lines 69-77.

**Issue and reproduction:** The default constructor is documented as creating a range with all coordinates zero, but it calls the base constructor, which produces the empty sentinel with minimum `+DBL_MAX` and maximum `-DBL_MAX`. A caller trusting the comment gets the opposite of an all-zero range.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Correct the comment to describe the empty sentinel.

**Rust handling:** `Default` is the empty sentinel, as the body does, and the rustdoc says so.

## CPP-072 — DRange::united of two empty ranges returns the universal range

**Affected files:** `src/openms/include/OpenMS/DATASTRUCTURES/DRange.h`, lines 178-195.

**Issue and reproduction:** Uniting two empty ranges passes their inverted sentinel corners through `setMinMax`, which normalises by swapping them. The result is the universal range from `-DBL_MAX` to `+DBL_MAX`: the union of two empty sets becomes everything.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Return an empty range when both operands are empty.

**Rust handling:** Transcribed and tested as `united_is_the_bounding_range`, with the quirk documented.

## CPP-073 — DRange::extend comment contradicts the collapse it performs

**Affected files:** `src/openms/include/OpenMS/DATASTRUCTURES/DRange.h`, line 310.

**Issue and reproduction:** The `@param` text states that resulting invalid minima and maxima are not fixed automatically, while the body collapses an inverted dimension to its centre. The class test asserts the collapse, so the comment is stale rather than the code.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Update the comment to describe the collapse.

**Rust handling:** The collapse is ported and documented.

## CPP-074 — DIntervalBase default constructor documented as infinite corners

**Affected files:** `src/openms/include/OpenMS/DATASTRUCTURES/DIntervalBase.h`, line 51.

**Issue and reproduction:** The constructor is documented as placing the corners at infinity; it uses the finite `numeric_limits` extrema. Code testing for infinity never matches.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Correct the comment to name the finite extrema.

**Rust handling:** The exact finite sentinel is ported and named in the rustdoc.

## CPP-075 — Class tests with unreachable or duplicated assertions

**Affected files:** `src/tests/class_tests/openms/source/FeatureHandle_test.cpp`, `StandardTypes_test.cpp`, `BinnedSpectrum_test.cpp`.

**Issue and reproduction:** The `FeatureHandle` `IndexLess` section sets `lhs.setUniqueId` twice and never sets `rhs`, so the ordering it claims to test is not exercised. `StandardTypes_test.cpp` constructs `PeakSpectrum` and `PeakMap` twice and never `Chromatogram`. The `BinnedSpectrum` constructor section title omits the `bool unit_ppm` parameter.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Set both operands in the ordering test, cover `Chromatogram`, and correct the section title.

**Rust handling:** The Rust tests set both operands, cover all three aliases and name every parameter.

## CPP-076 — rasterizeIMFrame clears the caller's image before the check that can throw

**Affected files:** src/openms/source/KERNEL/MSSpectrum.cpp:889-907 (total_pixels, std::fill, empty early return, then the IM/peak size Precondition)

**Issue and reproduction:** The output buffer is zero-filled at line 892, before the ion-mobility array's length is compared against the peak count at line 900. A spectrum whose IM array length differs from its peak count therefore throws Exception::Precondition *after* the caller's pre-allocated image has been wiped. A caller that catches the exception and keeps rendering shows a blank frame rather than the previous one. Trigger: any non-empty spectrum whose IM float array has a different number of entries than the peak list — reachable after a partial mzML load or after peaks are appended without extending the array.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Move the `im_data.size() != this->size()` check above the `std::fill`, next to the other precondition checks, so every throw happens before the buffer is touched.

**Rust handling:** `MSSpectrum::rasterize_im_frame` validates the grid, the IM array's presence and its length, and every peak and mobility value, before allocating the output `Vec<f32>`. A rejected call allocates nothing and mutates nothing. Asserted by `rasterize_im_frame_rejects_invalid_grids_and_arrays` in tests/spectrum_mobility.rs.

## CPP-077 — rasterizeIMFrame multiplies the bin counts without an overflow check

**Affected files:** src/openms/source/KERNEL/MSSpectrum.cpp:889 (`const Size total_pixels = im_bins * mz_bins;`), :892 (std::fill), :946 (pixel_idx)

**Issue and reproduction:** `im_bins * mz_bins` is unchecked `size_t` arithmetic. With bin counts whose product overflows, `total_pixels` wraps to a small value, `std::fill` clears only that prefix, and the per-peak write at line 946 computes `mz_bin * im_bins + im_bin` from the *unwrapped* bin counts and writes far outside whatever the caller allocated — a heap buffer overflow. Trigger: `spec.rasterizeIMFrame(buf, 1ull << 33, 1ull << 33, …)` on a 64-bit build, or any pair whose product exceeds SIZE_MAX; the bin counts come straight from a viewer's zoom level or a Python caller.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Check the product before using it, e.g. `if (mz_bins > std::numeric_limits<Size>::max() / im_bins) throw Exception::InvalidValue(...)`, and document a maximum raster size.

**Rust handling:** `ImFrameRaster::pixels` uses `checked_mul` and returns `Error::InvalidValue` on overflow; `ImFrameRaster::validate` additionally caps the product at `MSSpectrum::MAX_RASTER_PIXELS` (16 777 216, i.e. 64 MiB of f32). Both are asserted in `rasterize_im_frame_rejects_invalid_grids_and_arrays`.

## CPP-078 — rasterizeIMFrame casts a non-finite coordinate to Int64

**Affected files:** src/openms/source/KERNEL/MSSpectrum.cpp:925-941 (range filter and `static_cast<Int64>`)

**Issue and reproduction:** The range filter at line 928 uses `mz < min_mz || mz > max_mz || im < min_im || im > max_im`. Every one of those comparisons is false for NaN, so a NaN m/z or a NaN ion-mobility value is *not* skipped and reaches `static_cast<Int64>((mz - min_mz) * mz_scale)` at line 935. Converting a NaN to an integer type is undefined behaviour in C++; in practice it produces an unspecified value, and the only guard afterwards is the upper clamp, so a negative result indexes before the buffer. Trigger: any spectrum carrying a NaN m/z or a NaN entry in its ion-mobility array, which an mzML with a corrupt binary array or a preceding arithmetic step can produce.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Reject or skip non-finite values explicitly — `if (!std::isfinite(mz) || !std::isfinite(im)) continue;` before the range filter — and clamp the computed bin from below as well as above.

**Rust handling:** `rasterize_im_frame` checks every peak m/z, peak intensity and mobility value with the crate's `finite` helper before allocating, returning `Error::InvalidValue`. Rust's float-to-integer `as` cast is saturating rather than undefined, so even without the check nothing could index out of bounds. Asserted by the `nan_mobility` and `nan_mz` cases in `rasterize_im_frame_rejects_invalid_grids_and_arrays`.

## CPP-079 — sortByPositionPresorted handles an incomplete chunk list differently in its two branches

**Affected files:** src/openms/source/KERNEL/MSSpectrum.cpp:397-441 (the no-data-array stable_sort at 405-408 versus the chunk-driven path at 409-440)

**Issue and reproduction:** When the spectrum has no float, string or integer data arrays the function ignores `chunks` entirely and stable-sorts the whole peak list. When it has any data array it sorts only inside the given chunks and merges only `[chunks.front().start, chunks.back().end)`, then applies the resulting permutation with selectUnchecked — so any peak outside that span keeps its position and stays unsorted. The same chunk list therefore produces a fully sorted spectrum or a partially sorted one depending on whether a data array happens to be attached. The early return at line 401 (`chunks.size() == 1 && chunks[0].is_sorted`) compounds it: a single sorted chunk covering only part of the spectrum leaves the rest unsorted in both branches. Trigger: call `sortByPositionPresorted` with a chunk list built before the last batch of peaks was appended, once with and once without a data array attached.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Validate the chunk list up front — first chunk starts at 0, each chunk starts where the previous ended, none ends before it starts, and the last ends at size() — and throw Exception::Precondition otherwise; then use one code path so the branch cannot change the result.

**Rust handling:** `MSSpectrum::sort_by_position_presorted` performs exactly that validation and returns `Error::InvalidValue` for a gap, an overlap, a reversed run or a list that stops short of the peak count. With the tiling guaranteed, the two source branches produce the same permutation, so the port keeps only the chunk-aware one. Asserted by `presorted_sort_matches_plain_sort_and_rejects_bad_chunks`.

## CPP-080 — sortByPositionPresorted trusts is_sorted and feeds std::inplace_merge an unsorted range

**Affected files:** src/openms/source/KERNEL/MSSpectrum.cpp:415-438 (per-chunk stable_sort for !is_sorted, then the recursive inplace_merge)

**Issue and reproduction:** Chunks whose `is_sorted` flag is true are never sorted and never checked. `std::inplace_merge` requires both input ranges to be sorted by the comparator; passing an unsorted range is a precondition violation whose result is unspecified, and the function returns a spectrum silently in the wrong order — with the data arrays permuted to match, so the corruption is invisible to `checkDataArraySizes_` and to any later consistency check. Trigger: any caller that mis-tracks its runs, e.g. `MSSpectrum::Chunks` used across an intervening sort or insertion; `isSorted()` afterwards is the only way to notice.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Verify the claim before merging — `OPENMS_PRECONDITION(std::is_sorted(...))` at minimum, or an unconditional `std::is_sorted` check per chunk, which costs one linear scan against the merge's own O(n log k).

**Rust handling:** `sort_by_position_presorted` verifies every run marked `is_sorted` and returns `Error::UnsortedData` when the claim is false, before anything is mutated. Asserted by `presorted_sort_matches_plain_sort_and_rejects_bad_chunks`.

## CPP-081 — Chunks::add can record a chunk whose start exceeds its end

**Affected files:** src/openms/include/OpenMS/KERNEL/MSSpectrum.h:80-83 (`Chunks::add`), used at src/openms/source/KERNEL/MSSpectrum.cpp:415-419 and 433

**Issue and reproduction:** `add` records `{previous_end, spec_.size(), is_sorted}` from a `const MSSpectrum&` captured at construction. Nothing requires the spectrum to have grown: if peaks were removed (pop_back, erase, clear, select with a subset) between two `add` calls, `spec_.size()` is below the previous chunk's end and the new chunk has `start > end`. Both `Size` members are unsigned, so the inversion is not detectable by sign. `sortByPositionPresorted` then forms `select_indices.begin() + chunk.start` and `select_indices.begin() + chunk.end` and passes that inverted pair to `std::stable_sort` and `std::inplace_merge`, which is undefined behaviour; `chunk.start` can also exceed the vector's size outright, so the iterator is past the end. Trigger: `Chunks c(spec); …; c.add(true); spec.pop_back(); c.add(true); spec.sortByPositionPresorted(c.getChunks());`.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Throw Exception::Precondition in `add` when `spec_.size() < chunks_.back().end`, and validate `start <= end <= size()` for every chunk at the top of `sortByPositionPresorted`.

**Rust handling:** `Chunks::add` takes the spectrum as an argument (Rust cannot hold a shared reference across the mutation the source performs) and returns `Error::InvalidValue` when the spectrum has shrunk below the previous run's end. `sort_by_position_presorted` independently rejects `end < start`. Asserted by `presorted_sort_matches_plain_sort_and_rejects_bad_chunks`.

## CPP-082 — An empty ion-mobility array makes isSortedByIM report true and sortByIonMobility a silent no-op

**Affected files:** src/openms/source/KERNEL/MSSpectrum.cpp:383-395 (sortByIonMobility), :506-512 (isSortedByIM), :20-37 (checkDataArraySizes_, which permits an empty array)

**Issue and reproduction:** `containsIMData()` inspects only array *names*, so a float array named e.g. 'mean ion mobility array' with zero entries marks the spectrum as an IM frame. `isSortedByIM` then runs `std::is_sorted` over that empty range, which is vacuously true, and reports the spectrum sorted by ion mobility although no peak carries a mobility value. `sortByIonMobility` takes the same short-circuit and returns without touching anything. Neither relates the array to the peak count — `checkDataArraySizes_` deliberately exempts empty arrays, and in `sortByIonMobility` it only runs inside `sort()`, which the short-circuit skips. Trigger: build a spectrum with peaks and add an empty, correctly named ion-mobility float array (what an mzML reader produces when the binary array is absent or fails to decode); `isSortedByIM()` returns true.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Require the ion-mobility array to hold one value per peak in both functions — reuse `checkDataArraySizes_`-style logic on the IM array specifically, and throw Exception::Precondition when the sizes disagree — rather than short-circuiting on an empty range.

**Rust handling:** `MSSpectrum::is_sorted_by_im` and `sort_by_ion_mobility` both go through a private `checked_im_values`, which returns `Error::InvalidValue` when the ion-mobility array's length differs from the peak count, and additionally reject non-finite values (`std::is_sorted` with `<` also reports a NaN-containing array as sorted). Asserted by `ion_mobility_sorting_rejects_unusable_arrays`.

## CPP-083 — mergePeaks leaves the inherited range cache too narrow, undocumented

**Affected files:** src/openms/source/KERNEL/MSChromatogram.cpp, lines 548-565; contrast src/openms/include/OpenMS/KERNEL/MSChromatogram.h, lines 439-441

**Issue and reproduction:** mergePeaks replaces the peak vector without calling updateRanges() and without clearing the inherited RangeManager, so a chromatogram whose ranges were current afterwards reports a range that no longer contains its own points. select() documents the opposite, benign direction ('selecting a subset can leave them too wide -- call updateRanges()'); mergePeaks documents nothing. Trigger: a.updateRanges(); a.mergePeaks(b); a.getMaxRT() still returns a's pre-merge maximum although a now holds b's later points, and getMaxIntensity() misses every summed point. Any consumer that trusts getRange() after a merge, for example a plotting axis or a range-based filter, silently drops data.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Call updateRanges() at the end of mergePeaks, or clearRanges() so the stale values cannot be read as current, and add an @note matching the one on select().

**Rust handling:** There is no cache: MSChromatogram::range_manager() recomputes from the current points on every call, so the state cannot exist. Recorded under 'Native differences' in docs/CHROMATOGRAM_MERGE_SUPPORT.md and asserted by tests/chromatogram_merge.rs::source_update_ranges.

## CPP-084 — mergePeaks leaves the destination's data arrays mis-sized, so a later sort or select throws

**Affected files:** src/openms/source/KERNEL/MSChromatogram.cpp, lines 548-553; header @note at src/openms/include/OpenMS/KERNEL/MSChromatogram.h, line 472; checkDataArraySizes_ at MSChromatogram.cpp, lines 355-373

**Issue and reproduction:** mergePeaks assigns a new peak vector and never touches float_data_arrays_, string_data_arrays_ or integer_data_arrays_, which keep their pre-merge length. The header @note says peak-level metadata 'is not guaranteed to be correct after merging', which understates the consequence: the chromatogram now fails its own checkDataArraySizes_, so the next sortByPosition(), sortByIntensity(), sort(lambda) or select() throws Exception::Precondition. Trigger: a chromatogram with one 2-entry integer data array for 2 peaks, merged with a 1-peak chromatogram at a distinct RT; a.size() becomes 3 while the array stays at 2, and a.sortByPosition() then throws.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Clear the three data arrays inside mergePeaks, as clear() already does for the same reason, and change the @note to say the arrays are dropped. Alternatively add a parameter selecting drop-or-keep, and state in the @note that keeping them leaves the chromatogram in a state its own sort and select reject.

**Rust handling:** MergedDataArrays makes the choice explicit: Reject (the default) refuses a merge when either side carries a non-empty array, Drop removes them so the result validates, and Source reproduces the untouched arrays. tests/chromatogram_merge.rs::merge_annotation_array_policies asserts that the Source variant then fails validate() and sort_by_position().

## CPP-085 — setSumSimilarUnion has external linkage at global namespace scope

**Affected files:** src/openms/source/KERNEL/MSChromatogram.cpp, lines 474-515

**Issue and reproduction:** The merge helper is defined at global scope (outside namespace OpenMS, after a file-level `using namespace OpenMS;`) with no `static` and no anonymous namespace, so it has external linkage under the plain name ::setSumSimilarUnion. It is used only by MSChromatogram::mergePeaks in the same translation unit. Any other translation unit in the program that defines a global setSumSimilarUnion with a different body but the same signature is an ODR violation the linker will not diagnose, and the symbol is needlessly exported from libOpenMS. Trigger: link a second object file defining `OpenMS::MSChromatogram::Iterator setSumSimilarUnion(...)` at global scope; which definition mergePeaks calls is unspecified.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Mark it `static`, or move it into an anonymous namespace inside namespace OpenMS, or make it a private static member of MSChromatogram.

**Rust handling:** The equivalent union is written inline in MSChromatogram::merge_peaks_with_options in src/kernel/chromatogram_merge.rs; Rust has no global namespace and no ODR hazard.

## CPP-086 — mergePeaks takes a non-const reference to a chromatogram it only reads

**Affected files:** src/openms/include/OpenMS/KERNEL/MSChromatogram.h, line 479; src/openms/source/KERNEL/MSChromatogram.cpp, lines 548-564

**Issue and reproduction:** The signature is `void mergePeaks(MSChromatogram& other, bool add_meta = false)` and the @param tag marks `other` as [in,out], but the body only calls other.begin(), other.end() and other.getMZ() and never writes to it. A caller holding a `const MSChromatogram&` — for example one iterating an MSExperiment by const reference — cannot call mergePeaks without copying the whole chromatogram. Trigger: `void f(const MSExperiment& e, MSChromatogram& a) { a.mergePeaks(e.getChromatograms()[0]); }` does not compile.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Change the parameter to `const MSChromatogram& other` and the @param tag to [in]. The body needs no change beyond using const iterators.

**Rust handling:** merge_peaks takes `&MSChromatogram`. A documented side effect is that the aliased call `a.merge_peaks(&a, ..)`, which the source permits, is a compile error here.

## CPP-087 — The chromatogram stream operator prints an always-empty settings block

**Affected files:** src/openms/source/METADATA/ChromatogramSettings.cpp, lines 149-154, reached from src/openms/source/KERNEL/MSChromatogram.cpp, line 24

**Issue and reproduction:** `std::ostream& operator<<(std::ostream& os, const ChromatogramSettings& /*spec*/)` takes its argument unnamed and writes only '-- CHROMATOGRAMSETTINGS BEGIN --' and '-- CHROMATOGRAMSETTINGS END --'. MSChromatogram's operator<< streams the settings between its own banners, so every chromatogram dump advertises a settings section that has never contained a single setting: not the native id, the chromatogram type, the precursor, the product or any meta value. Trigger: build a chromatogram with a product m/z and a native id, stream it, and observe the two delimiters with nothing between them. A developer reading the dump concludes the record carries no settings.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Print the members the way the class test and MSSpectrum's settings dump lead a reader to expect (native id, chromatogram type, precursor, product, meta values), or delete the two delimiter lines so the empty section is not advertised.

**Rust handling:** Display for MSChromatogram reproduces the layout exactly, delimiters included, because the class test matches on the surrounding text; the reason is documented at the impl and in docs/CHROMATOGRAM_MERGE_SUPPORT.md, and tests/chromatogram_merge.rs::source_stream_layout pins the full string.

## CPP-088 — updateRanges warns that ranges were already up to date on its very first call

**Affected files:** src/openms/source/KERNEL/MSChromatogram.cpp, lines 517-546 (the OPENMS_ASSERTIONS blocks at 519-524 and 533-545)

**Issue and reproduction:** The debug-build check reads the old extrema as 0 whenever a dimension isEmpty(), computes the new extrema the same way, and logs 'Update ranges was called but ranges were already up-to-date' when the four values match. On a chromatogram that has never had updateRanges() called, both the before and after values are 0 whenever the result is also empty or all-zero, so the warning fires on a first, entirely necessary call. Trigger, in a build with OPENMS_ASSERTIONS: `MSChromatogram c; c.updateRanges();` warns, and so does a chromatogram holding a single peak at RT 0 with intensity 0.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Compare emptiness as well as the numeric extrema — for example skip the warning when any dimension was empty beforehand — so the check cannot confuse 'no range yet' with 'range unchanged'.

**Rust handling:** Not applicable: range_manager() has no cache, so there is no redundant-refresh condition to warn about. The source warning is recorded as neutralised in docs/CHROMATOGRAM_MERGE_SUPPORT.md.

## CPP-089 — BaseFeature::sortPeptideIdentifications comparator is not a strict weak ordering

**Affected files:** src/openms/source/KERNEL/BaseFeature.cpp:125-143 (the lambda; the empty branch is 128-131)

**Issue and reproduction:** The comparator returns true whenever its left argument is empty, regardless of the right one. For two empty PeptideIdentifications both comp(a,b) and comp(b,a) are true, so asymmetry — and therefore the strict weak ordering std::sort requires — is violated, which is undefined behaviour; libstdc++'s unguarded insertion sort relies on the ordering being sane and can walk past the end of the range. Trigger: a BaseFeature whose peptides_ holds two or more identifications without hits (exactly what the class test's own `ids.resize(n)` produces, and what featureXML yields for a PeptideIdentification element with no hits) and then sortPeptideIdentifications().

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Make emptiness a proper key: `if (p1.empty()) return !p2.empty();` keeping the existing `if (p2.empty()) return false;`.

**Rust handling:** sort_peptide_identifications orders on an explicit total key (hits present before hits absent, then best score) with a stable slice::sort_by, so several empty identifications are fine and keep their relative order. Asserted in tests/feature_identification.rs::sort_peptide_identifications_is_checked_and_atomic.

## CPP-090 — The same comparator mutates its arguments, so hits are sorted only where the sort happens to compare

**Affected files:** src/openms/source/KERNEL/BaseFeature.cpp:126-127

**Issue and reproduction:** The lambda takes `PeptideIdentification&` (non-const) and calls `p1.sort(); p2.sort();` inside the comparison. std::sort's comparator must not modify the objects it compares. The observable consequence is that hits are sorted only for elements the sort actually visits: a feature with exactly one attached identification is never compared, so sortPeptideIdentifications() leaves its hits unsorted while claiming to have sorted, and with more elements which hit ends up first can depend on the introsort's internal comparison sequence.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Sort every element's hits in a separate pass first, then call std::sort with a `const PeptideIdentification&` comparator.

**Rust handling:** All hits are sorted in a dedicated pass before the identifications are ordered. Asserted in bf_sort_peptide_identifications (hits [0.9, 0.5] inside the moved identification) and stated at the item's rustdoc.

## CPP-091 — Mixed isHigherScoreBetter flags make the sort comparator asymmetric

**Affected files:** src/openms/source/KERNEL/BaseFeature.cpp:136-142

**Issue and reproduction:** The score direction is read from the left operand only. If identification A has isHigherScoreBetter() true and B false, comp(A,B) and comp(B,A) can both be true, again violating the strict weak ordering and giving undefined behaviour. The header comment says the identifications are assumed to share a score type, but nothing checks it, and nothing prevents a caller from attaching identifications from two search engines to one feature.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Read the flag once (from the first non-empty element) and use it for every comparison, or reject/assert mixed score types explicitly.

**Rust handling:** sort_peptide_identifications returns Error::InvalidValue when the non-empty identifications disagree on higher_score_better, leaving the feature untouched. Asserted in sort_peptide_identifications_is_checked_and_atomic.

## CPP-092 — updateIDReferences / updateAllIDReferences lose identification matches on a throwing translation

**Affected files:** src/openms/source/KERNEL/BaseFeature.cpp:247-259 and src/openms/source/KERNEL/Feature.cpp:209-216

**Issue and reproduction:** updateIDReferences swaps id_matches_ into a local set first, emptying the member, then translates and reinserts one reference at a time. RefTranslator::translate throws Exception::MissingInformation for an unmapped reference when allow_missing is false (IdentificationData.cpp:1392-1399). If it throws on the k-th reference the feature keeps only the k-1 already translated ones and has silently dropped the remainder; primary_id_ may already have been overwritten too. Feature::updateAllIDReferences compounds this by updating each feature as it descends, so a failure deep in the subordinate tree leaves the tree half translated. A caller that catches the exception and continues — e.g. a map merge — then holds a feature whose annotations are partly gone with no indication.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Translate into a local set and a local optional primary ID and only assign/swap after the loop completes; in Feature::updateAllIDReferences, translate every subordinate's references into locals (or verify them all) before assigning anything.

**Rust handling:** BaseFeature::update_id_references builds both the translated primary ID and the translated match set into temporaries and commits only after every translation succeeded. Feature::update_all_id_references adds a read-only verification pass over the whole subordinate tree before the apply pass. Asserted in update_id_references_translates_atomically and update_all_id_references_covers_subordinates_or_changes_nothing.

## CPP-093 — ConsensusFeature::Ratio default constructor leaves ratio_value_ indeterminate

**Affected files:** src/openms/include/OpenMS/KERNEL/ConsensusFeature.h:92-111

**Issue and reproduction:** Ratio() is user-provided with an empty body and no member has a default member initialiser, so `double ratio_value_` is default-initialised, i.e. holds an indeterminate value. `ConsensusFeature::Ratio r; use(r.ratio_value_);` — or `ratios_.resize(n)` followed by reading the new elements — is undefined behaviour and in practice reads whatever was in that memory. The struct also carries a virtual destructor with no virtual functions while being stored by value in std::vector<Ratio>, which costs a vtable pointer per ratio and invites slicing.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** `double ratio_value_ = 0.0;` and `Ratio() = default;`; drop `virtual` from the destructor, since nothing derives from Ratio.

**Rust handling:** Ratio derives Default, so ratio_value is a defined 0.0 and there is no vtable; add_ratio and set_ratios additionally reject a non-finite value. Asserted in ratios_are_checked_and_replaced_atomically.

## CPP-094 — ConsensusFeature::setRatios takes a non-const lvalue reference for a copy-in setter

**Affected files:** src/openms/include/OpenMS/KERNEL/ConsensusFeature.h:265 and src/openms/source/KERNEL/ConsensusFeature.cpp:321-324

**Issue and reproduction:** The body is a plain copy (`ratios_ = rs;`) but the parameter is `std::vector<Ratio>&`. A temporary or a const vector therefore does not bind: `cf.setRatios(buildRatios());` and `void f(const std::vector<Ratio>& r) { cf.setRatios(r); }` both fail to compile, forcing callers to materialise a non-const named vector and implying, wrongly, that the setter may modify the argument.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** `void setRatios(const std::vector<Ratio>& rs);`, optionally with a `std::vector<Ratio>&&` overload.

**Rust handling:** set_ratios takes Vec<Ratio> by value, so a temporary is the natural call; it validates every ratio before replacing the stored list.

## CPP-095 — getAnnotationState reports MULTIPLE_SAME when only one identification actually has hits

**Affected files:** src/openms/source/KERNEL/BaseFeature.cpp:154-177

**Issue and reproduction:** The `peptides_.size() == 1` shortcut requires non-empty hits, but the fallback loop skips identifications without hits and then maps `seqs.size() == 1` to FEATURE_ID_MULTIPLE_SAME. A feature carrying one annotated identification plus one empty one therefore reports 'multiple IDs (identical)' although exactly one ID exists, which is inconsistent with the shortcut immediately above it. Trigger: `ids.resize(2); ids[0].setHits(std::vector<PeptideHit>(1, hit));` then getAnnotationState(). Consumers that render the state string or branch on SINGLE see the wrong category.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Count the identifications that contributed a sequence and return FEATURE_ID_SINGLE when that count is one, independent of how many empty identifications are attached.

**Rust handling:** Preserved deliberately, because tools may depend on the existing state, and documented at BaseFeature::annotation_state and in docs/FEATURE_IDENTIFICATION_SUPPORT.md. Asserted in annotation_state_counts_only_identifications_that_have_hits.

## CPP-096 — Feature::getConvexHull() is a const method that lazily mutates the object and hands out a mutable reference

**Affected files:** src/openms/include/OpenMS/KERNEL/Feature.h:99 and 175-179; src/openms/source/KERNEL/Feature.cpp:93-137

**Issue and reproduction:** getConvexHull() is declared const but recomputes and writes the mutable members convex_hull_ and convex_hulls_modified_, and returns a non-const ConvexHull2D&. Two consequences: (1) two threads calling it concurrently on the same Feature race on those members with no synchronisation, and the codebase parallelises feature loops with OpenMP (36 source files carry #pragma omp), so a const-looking read is not safe to share; (2) any caller holding a `const Feature&` can mutate the cached hull — `f.getConvexHull().clear();` compiles — leaving the object in a state that disagrees with convex_hulls_ until something sets the flag again.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Return `const ConvexHull2D&`, and either compute the overall hull eagerly in setConvexHulls/getConvexHulls or guard the lazy computation (std::call_once or a mutex) if the laziness must stay.

**Rust handling:** Feature::convex_hull() computes and returns an owned ConvexHull2D with no cache and no mutable-through-shared path, so neither the race nor the const mutation is expressible; convex_hulls_modified_ and convex_hull_ are recorded as not ported in the API mapping table. Covered by f_get_convex_hull and f_encloses.

## CPP-097 — MRMFeature lookup by unknown key mutates the map and returns the first feature

**Affected files:** src/openms/source/KERNEL/MRMFeature.cpp lines 82-85 (getFeature) and 125-128 (getPrecursorFeature), non-const overloads

**Issue and reproduction:** Both non-const accessors read `features_.at(feature_map_[key])`. `std::map::operator[]` default-inserts `key -> 0` when the key is absent, so a read accessor silently mutates the map and then returns feature 0 — an unrelated feature — or throws std::out_of_range when the list happens to be empty. Trigger: `MRMFeature f; Feature a; f.addFeature(a, "chromatogram1"); f.getFeature("typo");` returns the feature stored under "chromatogram1", and a subsequent getFeatureIDs() now reports two ids, one of which was never added. The const overloads use .at() and throw, so the same call is an error or a wrong answer depending only on the constness of the receiver.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Look the key up with find() and throw Exception::ElementNotFound (or return a pointer/optional) in both overloads, matching the const behaviour. If a mutable reference really is needed for a missing key, require an explicit insert call.

**Rust handling:** `MRMFeature::feature` / `feature_mut` / `precursor_feature` / `precursor_feature_mut` return `Error::MissingInformation` for an unknown key and never mutate. Covered by tests/mrm.rs::unknown_keys_are_missing_information.

## CPP-098 — MRMFeature::addFeature strands a feature when a key repeats, while the sibling class throws

**Affected files:** src/openms/source/KERNEL/MRMFeature.cpp lines 70-80 (addFeature) and 105-115 (addPrecursorFeature)

**Issue and reproduction:** Both push onto the vector first and then assign `feature_map_[key] = Int(size) - 1`. A repeated key re-points the map at the new element and leaves the previously keyed feature in the list with no key referring to it: the data is still present, counted by getFeatures().size(), but unreachable by any lookup and invisible to getFeatureIDs(). Trigger: `f.addFeature(a, "x"); f.addFeature(b, "x");` gives getFeatures().size() == 2 and getFeatureIDs() one entry. The sibling class in the same data model, MRMTransitionGroup::addTransition (MRMTransitionGroup.h:133), emplaces first and throws Exception::InvalidValue on exactly this condition, so the two halves disagree on whether a duplicate key is an error.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Use map::emplace and throw Exception::InvalidValue as addTransition does, or replace the feature in place at the existing index instead of appending.

**Rust handling:** `add_feature` / `add_precursor_feature` return `Error::InvalidValue` on a repeated key and store nothing; `add_feature_with(..., DuplicateKeyPolicy::SourceOverwrite)` reproduces the source's append-and-orphan exactly. Covered by tests/mrm.rs::duplicate_feature_keys_reject_by_default.

## CPP-099 — MRMTransitionGroup::isInternallyConsistent cannot report an inconsistent group in a release build

**Affected files:** src/openms/include/OpenMS/KERNEL/MRMTransitionGroup.h lines 291-297, with src/openms/include/OpenMS/CONCEPT/Macros.h line 91

**Issue and reproduction:** The function's three checks — equal transition/chromatogram counts, equal map sizes, and isMappingConsistent_() — are all OPENMS_PRECONDITION, which expands to nothing unless OPENMS_ASSERTIONS is defined. In an ordinary release build the body reduces to `return true;`, so the function declared `bool isInternallyConsistent() const` can never return false, and isMappingConsistent_() is never called. In an assertions build a violation throws Exception::Precondition instead of returning false, so there is no build in which the boolean result is meaningful. Trigger: a group with one addTransition and no addChromatogram returns true from any release build.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Evaluate the three conditions and return their conjunction — the declared return type already promises that — and keep the OPENMS_PRECONDITIONs alongside for the assertions build if the fail-fast behaviour is wanted there.

**Rust handling:** `is_internally_consistent` evaluates the three conditions and returns the answer; one pass over the chromatogram key map. Covered by tests/mrm.rs::internal_consistency_detects_every_source_condition and the ported class-test section group_is_internally_consistent.

## CPP-100 — IDScoresAsMetaValue writes the transition_names meta value twice

**Affected files:** src/openms/source/KERNEL/MRMFeature.cpp lines 142 and 152

**Issue and reproduction:** The function makes 43 setMetaValue calls covering only 42 distinct keys: `setMetaValue(id + "transition_names", idscores.ind_transition_names)` appears at line 142 and again, identically, at line 152. The second write is dead work (the value is the same), but its position — immediately before the block of `ind_*` keys starting at line 153 — suggests a copy/paste where a different OpenSwath_Ind_Scores member was intended, so a score may be silently absent from the output. Trigger: call IDScoresAsMetaValue and count the resulting meta values: 42, not 43.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Delete line 152, or establish which member the second write was meant to publish and write that instead. A static list of the 42 key suffixes would make the block auditable.

**Rust handling:** `id_scores_as_meta_value` writes each of the 42 distinct keys once; the suffixes are exposed as `OpenSwathIndScores::KEY_SUFFIXES` and the count is asserted. Covered by tests/mrm.rs::id_scores_as_meta_value_writes_the_whole_block.

## CPP-101 — getLibraryIntensity clamps entries the caller already had in the output vector

**Affected files:** src/openms/include/OpenMS/KERNEL/MRMTransitionGroup.h lines 319-333

**Issue and reproduction:** The function appends the transitions' library intensities to the caller's vector and then runs its negative-to-zero clamp over `result.size()` — the whole vector, not just the appended range. Any negative value the caller had put in the vector before the call is silently overwritten with zero. Trigger: `std::vector<double> r; r.push_back(-5.0); group.getLibraryIntensity(r);` leaves r[0] == 0.0 even for an empty group. The class test always passes a fresh vector, so this never shows up upstream.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Record the starting size before appending and clamp from that index, or build into a local vector and insert it at the end.

**Rust handling:** `library_intensity` returns a fresh `Vec<f64>` and clamps only its own entries; the difference is documented at the item and in docs/MRM_SUPPORT.md. Covered by the ported class-test section group_get_library_intensity.

## CPP-102 — subset re-keys precursor chromatograms by nativeID and throws on a group the class's own test builds

**Affected files:** src/openms/include/OpenMS/KERNEL/MRMTransitionGroup.h lines 355-359, with src/tests/class_tests/openms/source/MRMTransitionGroup_test.cpp lines 200-213

**Issue and reproduction:** subset copies every precursor chromatogram with `addPrecursorChromatogram(pc, pc.getNativeID())`, discarding the key it was actually stored under. Two precursor chromatograms whose nativeIDs coincide — which includes the very common case of both being empty, because nothing requires a chromatogram to carry a nativeID — then collide in the new group and addPrecursorChromatogram throws Exception::InvalidValue. Trigger: the class test's own getPrecursorChromatogram section stores chrom1 (nativeID never set, so empty) under "dummy1" and again under "dummy2"; calling subset() on that group throws, although the group is perfectly valid and hasPrecursorChromatogram("dummy1"/"dummy2") both answer true.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Carry the precursor chromatograms under their original keys by iterating precursor_chromatogram_map_ rather than the vector, which also makes the subset a faithful copy instead of a re-keying.

**Rust handling:** `subset` reproduces the source re-keying exactly and returns `Error::InvalidValue` on the collision instead of throwing; the quirk is documented at the item, in the preserved-conventions list, and pinned as a source anchor. Covered by tests/mrm.rs::subset_rebuilds_features_and_rekeys_precursors, which asserts the precursor lands under its nativeID and not under its stored key.

## CPP-103 — subsetDependent indexes chromatogram_map_.at() without the guard subset uses

**Affected files:** src/openms/include/OpenMS/KERNEL/MRMTransitionGroup.h line 397, compared with lines 348-351 in subset

**Issue and reproduction:** subset guards the chromatogram transfer with `if (this->hasChromatogram(tr.getNativeID()))`; subsetDependent calls `chromatograms_[chromatogram_map_.at(tr_it->getNativeID())]` unconditionally. A selected transition with no chromatogram under its native ID therefore raises an uncaught std::out_of_range, which escapes as a bare STL exception rather than an OpenMS Exception, so the TOPP top-level handler reports it without file, line or function context. Trigger: a group with addTransition(t, "t1") where t.getNativeID() == "t1" and no addChromatogram, then subsetDependent({"t1"}). This is never exercised upstream (see the next issue).

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Apply the same hasChromatogram guard as subset, or throw Exception::ElementNotFound with the offending native ID so the failure is diagnosable.

**Rust handling:** `subset_dependent` returns `Error::MissingInformation` naming the key, and builds into a temporary so the receiver is untouched. Covered by tests/mrm.rs::subset_dependent_keeps_whole_features_and_drops_precursors.

## CPP-104 — The subsetDependent class-test section tests subset, leaving subsetDependent with no coverage

**Affected files:** src/tests/class_tests/openms/source/MRMTransitionGroup_test.cpp lines 328-352

**Issue and reproduction:** START_SECTION(MRMTransitionGroup subsetDependent(std::vector<std::string> tr_ids)) builds its fixture and then calls `mrmtrgroupsub = mrmtrgroup.subset(transition_ids);` — subset, not subsetDependent. The section duplicates the preceding subset section with a two-element id list and asserts nothing that subsetDependent does differently, so the two behaviours that distinguish it (copying each MRMFeature whole rather than rebuilding it, and dropping the precursor chromatograms) and its unguarded chromatogram_map_.at() are all untested. The upstream suite reports the section as passing.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Call subsetDependent in that section, and add a case with a selected transition that has no chromatogram plus a case asserting that a copied MRMFeature still carries the sub-features of transitions outside tr_ids.

**Rust handling:** The section is ported exactly as written (against subset, so its transcribed literals stay meaningful) as tests/mrm.rs::group_subset_dependent_section, and subset_dependent itself is covered separately by tests/mrm.rs::subset_dependent_keeps_whole_features_and_drops_precursors, which asserts the whole-feature copy, the dropped precursors and the missing-chromatogram error.

## CPP-105 — IndexedMzMLHandler copy constructor silently drops both native-id maps

**Affected files:** src/openms/source/FORMAT/HANDLERS/IndexedMzMLHandler.cpp:76-88 (member list at IndexedMzMLHandler.h:57-64); consumer at src/openms/include/OpenMS/KERNEL/OnDiscMSExperiment.h:97-103

**Issue and reproduction:** The copy constructor's initialiser list copies filename_, spectra_offsets_, chromatograms_offsets_, index_offset_, spectra_before_chroms_, a freshly reopened filestream_, parsing_success_ and skip_xml_checks_, but omits spectra_native_ids_ and chromatograms_native_ids_, which are therefore default-constructed empty. Every getMSSpectrumByNativeId / getMSChromatogramByNativeId call on a copied handler throws Exception::IllegalArgument for every identifier, while getNrSpectra and index-based access still work, so the failure looks like a missing spectrum rather than a broken copy. The class is meant to be copied: the comment at line 82 says reopening rather than copying the stream 'is critical for parallel access to the same file', OnDiscMSExperiment's copy constructor copies indexed_mzml_file_ by value, and OnDiscMSExperiment.h:65-67 recommends '#pragma omp parallel for firstprivate(ondisc_map)'. Trigger: copy a handler (or an OnDiscMSExperiment) opened on any indexed mzML and call getSpectrumByNativeId with an id the original resolves.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Add spectra_native_ids_(source.spectra_native_ids_) and chromatograms_native_ids_(source.chromatograms_native_ids_) to the initialiser list, in declaration order. A class test that copies the handler and repeats the by-native-id assertions would have caught it: the existing copy-constructor section only compares arrays fetched by index.

**Rust handling:** Not reproduced. There is no copy constructor; an independent reader is an independent IndexedMzMLHandler::open, and the ordered identifier maps are ordinary owned state built in open_with_limits. tests/indexed_mzml_handler.rs::a_second_handler_on_one_file_reads_the_same_data asserts that the second handler resolves both native identifiers, which is exactly the assertion the source would fail.

## CPP-106 — IndexedMzMLHandler::openFile accumulates index state instead of replacing it

**Affected files:** src/openms/source/FORMAT/HANDLERS/IndexedMzMLHandler.cpp:20-62 (parseFooter_) and :92-101 (openFile)

**Issue and reproduction:** parseFooter_ push_backs into spectra_offsets_ and chromatograms_offsets_ and emplaces into both native-id maps without clearing any of them, and openFile does not clear them either — it only closes and reopens the stream. Opening a second valid indexed mzML on the same object leaves getNrSpectra() reporting the sum of both files, with the first file's offsets now interpreted against the second file's stream, so getMSSpectrumById returns garbage or throws from the decoder for the stale entries. spectra_before_chroms_ is likewise recomputed from a mixed vector. The class test does call openFile repeatedly (IndexedMzMLFile_test.cpp:100-106), but its first call throws FileNotFound out of findIndexListOffset and its second returns early because findIndexListOffset yields -1, so nothing is appended and the defect stays hidden. Trigger: handler.openFile(a); handler.openFile(b); with two valid indexed mzML files.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Clear spectra_offsets_, chromatograms_offsets_, spectra_native_ids_ and chromatograms_native_ids_ at the top of parseFooter_ (or in openFile before calling it), and reset index_offset_ and spectra_before_chroms_ on the early-return path too.

**Rust handling:** Not reproduced. There is no reopen: IndexedMzMLHandler::open always constructs a fresh handler, so no container can carry state across files. tests/indexed_mzml_handler.rs::opening_replaces_rather_than_accumulates asserts that a second open of the same two-spectrum fixture still reports exactly 2 spectra and no third offset.

## CPP-107 — Record read length is unchecked in both directions and the read result is never inspected

**Affected files:** src/openms/source/FORMAT/HANDLERS/IndexedMzMLHandler.cpp:139-168 (getChromatogramById_helper_) and :199-228 (getSpectrumById_helper_)

**Issue and reproduction:** Both helpers compute std::streampos readl = endidx - startidx from two values that come straight out of the file's own footer index, then call new char[readl + std::streampos(1)]. Nothing checks that endidx >= startidx, that either offset lies inside the file, or that the difference is a plausible record size. A decreasing pair makes the length negative, so the new-expression throws std::bad_array_new_length — not an OpenMS exception, so it escapes every Exception:: handler a caller installed. An oversized pair is an unbounded allocation driven by untrusted input. Worse, filestream_.seekg/read are issued without checking gcount() or the stream state: when the range runs past the end of the file the read stops early, the tail of the buffer keeps the indeterminate values left by new char[], only buffer[readl] is set to '\0', and std::string text(buffer) then scans that indeterminate memory. Trigger: any indexed mzML whose <indexList> offsets are corrupt, truncated or hostile — precisely the untrusted input this class exists to read.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Validate the pair before allocating: require startidx <= endidx, clamp endidx to the file size obtained once at openFile, and reject a span above an explicit ceiling with Exception::ParseError. Replace new char[] / delete[] with std::vector<char> (value-initialised, exception-safe) and check filestream_.gcount() == readl after the read, throwing ParseError otherwise.

**Rust handling:** Reproduced as checked errors. record_range rejects a decreasing pair (Error::Parse 'index offsets do not increase across the record'), an end past the file length (Error::Parse 'record byte range extends past the end of the file') and a span above RecordReadLimits::max_record_bytes (Error::InvalidValue), all before any allocation; read_range then uses try_reserve_exact and verifies the byte count actually read. open_with_limits additionally rejects an indexListOffset or index entry beyond the file. Covered by decreasing_index_offsets_are_rejected_before_allocating, offsets_past_the_end_of_the_file_are_rejected and an_oversized_record_range_is_refused_before_reading, which also asserts the handler is unchanged after the refusal.

## CPP-108 — Chromatogram out-of-range message reports the spectrum count

**Affected files:** src/openms/source/FORMAT/HANDLERS/IndexedMzMLHandler.cpp:132-137 (and the same off-by-one wording at :192-197)

**Issue and reproduction:** getChromatogramById_helper_ checks chromToGet >= (int)getNrChromatograms() but its Exception::IllegalArgument message reads 'id needs to be smaller than the number of spectra' and interpolates getNrSpectra() as 'maximal allowed'. On the class test's own fixture (2 spectra, 1 chromatogram) asking for chromatogram 1 is rejected with 'maximal allowed is 2', which is both the wrong noun and a larger number than the check accepts — actively misleading while debugging an index. Separately, both helpers say 'maximal allowed is N' where the largest accepted index is N-1.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Use 'number of chromatograms' and getNrChromatograms() in the chromatogram helper, and report getNrX() - 1 (or reword to 'must be less than N') in both.

**Rust handling:** Not reproduced. out_of_range formats '<kind> index {index} is not below the indexed count {count}' from the same RecordKind and the same vector the bound was taken from, so the noun and the number cannot disagree, and 'not below' states the relation rather than an off-by-one maximum.

## CPP-109 — Record XML parse errors are discarded, and the handler always feeds the parser ill-formed input for the last record

**Affected files:** src/openms/source/FORMAT/HANDLERS/MzMLSpectrumDecoder.cpp:536-556 (domParseString_); ranges produced at src/openms/source/FORMAT/HANDLERS/IndexedMzMLHandler.cpp:142-155 and :202-215

**Issue and reproduction:** domParseString_ constructs a XercesDOMParser and calls parse() without ever calling setErrorHandler, so well-formedness errors go to the default handler and are silently dropped; only a null document element is caught. This is not a latent risk but a permanent condition: the handler's byte range for the last record of each kind deliberately runs to the start of the other list or to <indexList>, so the text handed to the parser ends with '</spectrum></spectrumList><chromatogramList ...>' or '</chromatogram></chromatogramList></run></mzML>'. Every such fetch parses input with content after the root element. The same silence means a genuinely truncated, mis-indexed or corrupt record yields an empty MSSpectrum with no peaks and no error, indistinguishable from a legitimately empty scan. Trigger: getMSSpectrumById for the last spectrum of any indexed mzML, and any file whose index offsets do not land exactly on a record.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Install an error handler that throws Exception::ParseError for errors and fatal errors, and trim the byte range in the two helpers at the record's own closing tag before handing it over, so the parser never sees trailing content.

**Rust handling:** Reproduced as a loud failure. record_xml requires the range to start with the expected element and trims it at the matching closing tag, then the assembled single-record document goes through src/format/mzml.rs, whose errors are returned rather than discarded; an unterminated record is Error::Parse 'record is not closed inside its byte range'. Covered by the_last_record_is_trimmed_at_its_own_closing_tag (asserting the trimmed spectrum contains no '</spectrumList>' and the trimmed chromatogram no '</run>') and an_offset_that_is_not_a_record_start_is_rejected.

## CPP-110 — getMSChromatogramById(int) recomputes ranges twice

**Affected files:** src/openms/source/FORMAT/HANDLERS/IndexedMzMLHandler.cpp:278-284, with the inner call at :297-302

**Issue and reproduction:** The value-returning overload calls getMSChromatogramById(id, c), which already ends with c.updateRanges() at line 301, and then calls c.updateRanges() again on the result. updateRanges is a full pass over every point of the chromatogram, so every by-index chromatogram fetch pays for it twice. The spectrum counterpart at :246-251 does not do this, so the two families are also inconsistent. Minor: wasted work only, no wrong answer.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Delete the second c.updateRanges() at line 282, matching getMSSpectrumById.

**Rust handling:** Does not arise. src/kernel.rs computes ranges on demand rather than caching them behind an updateRanges() call, so there is no range refresh to perform once, twice or not at all.

## CPP-111 — MapConversion::convert(PeakMap) sorts and indexes past the end of its vector

**Affected files:** src/openms/source/KERNEL/ConversionHelper.cpp:24-44, src/openms/source/KERNEL/MSExperiment.cpp:744, src/openms/include/OpenMS/KERNEL/MSExperiment.h:171

**Issue and reproduction:** n is clamped against MSExperiment::getSize(), which sums the peaks of every spectrum plus every chromatogram point, but the vector that is then partially sorted comes from get2DData, which skips every spectrum whose MS level is not 1. With any MS2 spectrum or chromatogram present the middle iterator of std::partial_sort is past the end (undefined behaviour) and the following loop reads tmp[element_index] past the end. The default argument Size(-1) reaches this path on every call that does not pass n explicitly, so any DDA run converted with the default triggers it.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Clamp against tmp.size() after get2DData has filled it, not against getSize().

**Rust handling:** MapConversion::peak_map_to_consensus clamps against the number of points actually collected; tests/conversion_helper.rs::peak_map_to_consensus_clamps_against_collected_points asserts 12 elements for a 12-MS1-peak run that also holds an MS2 spectrum and a chromatogram (getSize() would be 14).

## CPP-112 — FeatureMap::swap and ConsensusMap::swap do not swap the meta values

**Affected files:** src/openms/source/KERNEL/FeatureMap.cpp:323, src/openms/source/KERNEL/ConsensusMap.cpp:393

**Issue and reproduction:** Both swap the elements, ranges, DocumentIdentifier, UniqueIdInterface, the UniqueIdIndexer, every record vector and id_data_, but never the inherited MetaInfoInterface. After a swap each map carries the other map's data with its own meta values. clear(true) and operator== both treat the meta values as part of the map, so the omission is inconsistent with the rest of the class rather than deliberate.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Add MetaInfoInterface::swap(from) (or std::swap on the meta pointer) to both swap bodies.

**Rust handling:** FeatureMap::swap and ConsensusMap::swap reproduce the omission exactly and say so at the item, pointing callers at std::mem::swap for a complete exchange; tests/map_operations.rs::fm_swap and ::cm_swap assert that a meta value stays with its original map.

## CPP-113 — ConsensusMap::split indexes its result vector with the map index

**Affected files:** src/openms/source/KERNEL/ConsensusMap.cpp:702, 780, 801

**Issue and reproduction:** fmaps is sized by column_description_.size() but indexed by the feature handle's and the identification's map_index. The two agree only when the column headers happen to be keyed 0..n-1. A consensusXML whose headers are keyed sparsely (which nothing forbids, and which setPrimaryMSRunPath can itself create) makes fmaps[it->first] read and write past the end of the vector.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Build a map index to position mapping from column_description_'s keys, or check the index against fmaps.size() and throw.

**Rust handling:** ConsensusMap::split maps every index through check_column and returns Error::InvalidValue for one that does not name a column; tests/map_operations.rs::cm_split_error_paths covers it.

## CPP-114 — ConsensusMap::split dereferences an empty map in the isobaric branch

**Affected files:** src/openms/source/KERNEL/ConsensusMap.cpp:741

**Issue and reproduction:** (*new_feats.begin()).second is taken without checking that the consensus feature had any feature handle. A consensus feature with peptide identifications but no handles — which nothing in the class prevents — dereferences end().

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Guard on new_feats.empty() before the dereference, or fold it into the existing min_index check.

**Rust handling:** The min_index != Some(0) guard rejects an empty handle set before any identification is routed, and the isobaric target lookup is a checked Option.

## CPP-115 — ConsensusMap::appendRows pairs column headers by position, not by column index

**Affected files:** src/openms/source/KERNEL/ConsensusMap.cpp:85-92

**Issue and reproduction:** After merging rhs's headers, the loop advances an iterator over the merged map and a second over rhs's map in lockstep and sets getColumnHeaders()[it->first].size = it->second.size + it2->second.size. The two iterators are at the same ordinal position, not at the same column index, so a merge adds the size of an unrelated column and renames only the first min(|merged|, |rhs|) headers to mergedConsensusXMLFile. The remaining headers keep a filename that no longer describes their contents.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Iterate rhs's headers and accumulate into the merged entry with the same key.

**Rust handling:** ConsensusMap::append_rows reproduces the positional zip and documents it; tests/map_operations.rs::cm_append_rows asserts that header 0 is renamed while header 1 keeps "m2".

## CPP-116 — ConsensusMap::setPrimaryMSRunPath writes by position through a default-inserting map

**Affected files:** src/openms/source/KERNEL/ConsensusMap.cpp:518-536

**Issue and reproduction:** The count check fires only when column_description_ is non-empty, and the write loop uses column_description_[i] with i counting from zero. std::map::operator[] default-inserts, so a map whose headers are keyed otherwise (say {5}) passes the count check for one path and then gains a second, empty column at key 0 instead of being renamed. An empty map silently accepts any number of paths and invents that many columns.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Write into the existing keys in order, or reject a map whose keys are not 0..n-1.

**Rust handling:** ConsensusMap::set_primary_ms_run_path reproduces the behaviour into a temporary and documents it; tests/map_operations.rs::cm_primary_ms_run_path asserts the sparse-key case produces two columns.

## CPP-117 — MSExperiment::getPrimaryMSRunPath joins using the unstripped file:/// path

**Affected files:** src/openms/source/KERNEL/MSExperiment.cpp:922-925

**Issue and reproduction:** actual_path strips a leading file:/// only to decide whether the separator should be a backslash or a slash; the location that is pushed is then built from the original path, so a source file whose path is file:///C:/data yields file:///C:/data/run.mzML. FeatureMap::setPrimaryMSRunPath and ConsensusMap::setPrimaryMSRunPath pass that string straight to File::exists, which cannot find it, so the fallback list is silently used instead of the experiment's own run path.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Build ms_run_location from actual_path.

**Rust handling:** The private usable_experiment_run_path helper reproduces the spelling so the existence test sees what the source tests, and records the defect in its doc comment and in tests/data/map_operations_provenance.json.

## CPP-118 — UniqueIdIndexer::resolveUniqueIdConflicts can loop forever

**Affected files:** src/openms/include/OpenMS/CONCEPT/UniqueIdIndexer.h:138-163

**Issue and reproduction:** while (uniqueid_to_index_.contains(unique_id)) { setUniqueId(); ... } has no bound. A UniqueIdGenerator that has been seeded identically in two places, or one whose state has been reset, can keep producing an already-used value and hang the merge. The counter invalid_uids also does not count the IDs it assigns to previously unassigned elements, so the number it returns (and that FeatureMap::operator+= logs) understates what changed.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Bound the redraw loop and either count every assignment or rename the return value.

**Rust handling:** The private resolve_unique_id_conflicts bounds each element at MAX_UNIQUE_ID_REDRAWS (64) and returns Error::InvalidValue past it, while keeping the source's counting rule; documented in map_operations.rs.

## CPP-119 — FeatureMap::getPrimaryMSRunPath does not clear its output argument

**Affected files:** src/openms/source/KERNEL/FeatureMap.cpp:452-464

**Issue and reproduction:** toFill is only assigned when the spectra_data meta value exists; otherwise the caller's previous contents survive, and the UNKNOWN placeholder is appended only when the result is empty. A caller that reuses one StringList across maps therefore reads the previous map's paths back as if they belonged to the current one, with no diagnostic.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Clear toFill on entry.

**Rust handling:** FeatureMap::primary_ms_run_path returns a fresh Vec, so the question cannot arise; the source behaviour is stated in the item's rustdoc.

## CPP-120 — mzML list `count` attribute drives an unvalidated container reserve

**Affected files:** src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp:965, :978, :996, :1012, :1017; src/openms/include/OpenMS/FORMAT/HANDLERS/XMLHandler.h:393-398, :538-543; src/openms/source/FORMAT/HANDLERS/StringManager.cpp:109-112; src/openms/source/KERNEL/MSExperiment.cpp:520-523

**Issue and reproduction:** `attributeAsInt_` returns a signed `Int` straight from `xercesc::XMLString::parseInt` with no range or sign validation, and MzMLHandler passes that value directly into container reservations from an untrusted file: `bin_data_.reserve(attributeAsInt_(attributes, s_count))` (:1017), `exp_->reserveSpaceSpectra(scan_count_total_)` (:978) and `exp_->reserveSpaceChromatograms(chrom_count_total_)` (:1012), where `reserveSpaceSpectra(Size)` forwards to `spectra_.reserve(s)`. A `count="-1"` converts to `SIZE_MAX` on the unsigned parameter and raises `std::length_error` (not an OpenMS exception, so it escapes the format-error path callers expect); a plausible-looking `count="2000000000"` causes an immediate multi-gigabyte reservation before a single child element has been parsed. A one-line edit to a list attribute in an otherwise valid mzML is enough. Source-reviewed, not executed.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Clamp or validate the declared count before spending it: read it with `asUInt_`/a checked conversion, reject negatives, and cap the reserve at a configured ceiling (or at a bound derived from remaining input size) rather than at the attribute's face value — the attribute is only a capacity hint, so clamping costs nothing.

**Rust handling:** The Rust reader parses the same attribute as `usize` (so negatives and non-numerics are rejected outright) and deliberately spends it on no allocation at all: records are bounded by `ReadOptions::max_records` and arrays by `max_total_arrays` as each child actually opens, and the productList/scanWindowList/scanList counts are checked against the remaining parameter budget before any parameter is allocated. `tests/mzml_list_counts.rs::declared_counts_remain_resource_ceilings_before_allocation` pins that a 1e9 declared spectrumList or binaryDataArrayList count is accepted while allocating nothing.

## CPP-121 — rasterizeRTMZ's SUM aggregation is not reproducible across thread counts

**Affected files:** src/openms/source/KERNEL/MSExperiment.cpp:355-518

**Issue and reproduction:** rasterizeRTMZ chooses between a single-threaded branch (`if (num_threads <= 1)` at :355, writing straight into the output) and an OpenMP branch (:401-517) from a heuristic over omp_get_max_threads(), a 4 MB memory budget and sqrt(num_spectra) (:344-352). For RasterAggregation::MAX the two agree, because a maximum is associative and commutative. For SUM they do not: each thread accumulates its own f32 partial sum into a per-thread buffer (:460) and the merge adds those partials into the pixel (:499), so the f32 summation order — and the rounded pixel value — depends on how the spectra happened to be distributed. The same input therefore rasterizes to different SUM images on machines with different core counts, and differently again from the single-threaded branch. Nothing in the header's documentation warns of this, and the class test does not exercise the parallel branch's SUM output against the serial one.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Either document SUM as thread-count dependent, or make the merge deterministic — accumulate in f64 per thread and narrow once at the end, or partition spectra to threads by a fixed rule and merge in thread-id order with a fixed traversal. Accumulating in double costs one extra buffer's width but removes the nondeterminism entirely.

**Rust handling:** The port is serial and reproduces the single-threaded branch exactly, so its SUM output is deterministic. The difference is stated rather than hidden: src/kernel/experiment_mobility.rs:41-50 and docs/EXPERIMENT_MOBILITY_SUPPORT.md:310-322 say the serial image matches the source only for Max, and name the two lines that make Sum order-dependent.

## CPP-122 — rasterizeRTMZ's two negative-bin guards are dead code that mask an unchecked precondition

**Affected files:** src/openms/source/KERNEL/MSExperiment.cpp:362, :376, :388, :425-428, :453, :472

**Issue and reproduction:** rasterizeRTMZ skips a whole spectrum when rt_bin < 0 (:362 serial, :425-428 parallel) and one peak when mz_bin < 0 (the `mz_bin >= 0` tests at :376 and :388 serial, :453 and :472 parallel). Neither can fire on the input the function documents: the spectrum window comes from RTBegin(min_rt)/RTEnd(max_rt) and the peak window from MZBegin(min_mz)/MZEnd(max_mz), all lower_bound/upper_bound calls, and both scales are positive because min >= max is rejected at :288-295. They can fire only when the run is not sorted by RT or a spectrum is not sorted by m/z — exactly the case in which those binary searches return meaningless positions. So the guards do not protect against anything they can actually see, while the real precondition (the header's "@note should be sorted by RT and m/z") stays unchecked, and on unsorted input the function silently drops arbitrary spectra and peaks instead of reporting the violation.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Check the precondition instead of guarding its symptom: verify the RT order of the visited range and the m/z order of each contributing spectrum up front and throw Exception::Precondition (or at least emit a warning), then drop the two unreachable branches. If the guards are kept as belt-and-braces, a comment should say they are unreachable on sorted input.

**Rust handling:** rasterize_rt_mz calls rt_begin/mz_begin, which check_sorted and return Error::UnsortedData before any binning (src/kernel.rs:894-909 and :665-676), so the branch is unreachable in the port too — but for a stated reason rather than by accident, and unsorted input is reported rather than silently thinned. Recorded at src/kernel/experiment_mobility.rs:654-673 and docs/EXPERIMENT_MOBILITY_SUPPORT.md:323-338, with the note that Rust's saturating as-usize cast would clamp a negative bin to 0 rather than skip, were it ever reached.

## CPP-123 — None found — no new C++ source defect

**Affected files:** n/a

**Issue and reproduction:** This was a documentation- and test-accuracy task; no upstream C++ behaviour was re-examined, so nothing new is reportable for OpenMS_CPP_ISSUES.md. The four defects docs/ON_DISC_EXPERIMENT_SUPPORT.md already records for OnDiscMSExperiment.h/.cpp are unchanged and were not re-litigated.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** n/a

**Rust handling:** n/a

## CPP-124 — .ibd array reads are sized from the declared count with no comparison against the file's length

**Affected files:** src/openms/source/FORMAT/HANDLERS/ImzMLHandlerHelper.cpp:122-186 (readMzArray), :190-262 (readFloatVector_)

**Issue and reproduction:** Both readers validate only the element count, against MAX_IBD_ARRAY_ELEMENTS = 100,000,000, and then size the output from that count. Nothing compares offset + count * element_width against the actual length of the .ibd, which is the one fact that makes the request answerable. A malformed or hostile IMS:1000103 therefore commits the memory first and discovers the truncation only when fread comes up short: 800 MB for a float64 array at the ceiling, and 1.2 GB for a float32 one, because readMzArray stages a std::vector<float> alongside the std::vector<double> it widens into. Offsets and lengths in imzML come from one file and index into another, so this is reachable from an untrusted dataset without any corruption of the XML being detectable first.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** In validateCount_, or in a new preflight beside it, stat or fseek/ftell the .ibd once when it is opened, keep its length on the handler, and reject offset + count * width > length before resize/fread. Two comparisons and no extra I/O.

**Rust handling:** ImzMLBinaryIO::preflight rejects an unknown data type, a count above max_array_elements, a count above usize, a stored size above max_array_bytes, an offset + bytes that overflows u64, and any range that leaves the length measured when the .ibd was opened — all before allocating, and the allocation then uses try_reserve_exact. Covered by a_range_past_the_end_of_the_ibd_is_refused, an_element_count_above_the_ceiling_is_refused_before_allocating, an_offset_that_overflows_the_range_is_refused, a_byte_length_above_the_configured_limit_is_refused and a_truncated_ibd_fails_the_preflight.

## CPP-125 — ImzMLMeta::mz_data_type and int_data_type stay empty on an index-only load

**Affected files:** src/openms/source/FORMAT/HANDLERS/ImzMLHandler.cpp:248-252; src/openms/include/OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h:64-70; src/openms/source/FORMAT/ImzMLFile.cpp (loadSpectraIndex path)

**Issue and reproduction:** The two fields are assigned inside the `if (handler_.ibd_ && handler_.decode_ibd_ && ...)` decode branch of ImzMLInterceptConsumer::consumeSpectrum. ImzMLFile::loadSpectraIndex passes decode_ibd = false, so an index-only load returns an ImzMLMeta whose mz_data_type and int_data_type are empty strings even though every spectrum's XML declared MS:1000521. The header documents them as "first occurrence, dataset-level summary" with no caveat, and the value is available from the parse, not from the decode: ArrayMeta::dt is already set at endElement("binaryDataArray"). A caller that reads the index to decide how to allocate cannot learn the element width without decoding a spectrum it does not want.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Record the first non-UNKNOWN dt in ImzMLHandler::onEndElement("spectrum"), where cur_mz_meta_ and cur_int_meta_ are already snapshotted, instead of in the consumer's decode branch. The decode-branch assignment then becomes redundant and can go.

**Rust handling:** read_index sets ImzMLMeta::mz_data_type / int_data_type from the first spectrum in document order that declares a type, so a metadata-only parse is informative; documented at the fields and under native difference 4 of docs/IMZML_HANDLER_SUPPORT.md. Both fixtures report Float32 either way, so this does not change any upstream expectation.

## CPP-126 — An unnamed external auxiliary array vanishes from the index with no trace

**Affected files:** src/openms/source/FORMAT/HANDLERS/ImzMLHandler.cpp:264-269 (warn at decode) and :374-390 (omit from the index entry)

**Issue and reproduction:** An external binaryDataArray that carries neither MS:1000786 nor a child of MS:1000513 is warned about when a spectrum is decoded in memory, and is separately skipped when the ImzMLSpectrumIndex entry is built. The warning comes from spec_ims_, which OnDiscImzMLExperiment does not retain; the index it does retain simply has one fewer aux entry. So an on-disc consumer inspecting getIndex(i).aux cannot tell that an array was present in the file and discarded, and an index-only load emits no warning at all because consumeSpectrum's aux block is inside the decode guard. The file's array count and the index's disagree silently.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Keep a count of the dropped arrays on ImzMLSpectrumIndex, or log the skip from endElement("binaryDataArray") where the name is known, so the omission is visible to whoever holds the index.

**Rust handling:** ImzMLSpectrumIndex::unnamed_aux counts them, ImzMLSpectrumIndex::aux still holds only named arrays so aux.len() matches the source, and ImzMLHandler::spectrum reports one AuxSkipReason::Unnamed per dropped array in DecodedSpectrum::skipped_aux. Covered by skippable_auxiliary_arrays_are_reported_and_the_spectrum_survives.

## CPP-127 — verifyIbdUuid_ carries two contradictory doc comments, one describing behaviour the function does not have

**Affected files:** src/openms/source/FORMAT/ImzMLFile.cpp:77-87

**Issue and reproduction:** Two comment blocks sit immediately above the same function. The first says the .ibd "must begin with the 16-byte UUID declared in the .imzML XML" and that a mismatch means the files do not belong together, so it should "reject loudly instead of silently decoding garbage offsets". The second, added later, says the check is advisory and is "reported as a loud WARNING rather than a hard error so that legacy / non-conformant datasets still load". The code does the second. A reader who stops at the first comment will believe a mismatched pair is rejected and will not add the check their own caller needs.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Delete the superseded first block. The surviving comment already states the policy and the reason for it.

**Rust handling:** ImzMLHandler::uuid_status returns UuidStatus::{Match, Mismatch, NotDeclared, IbdTooShort} and open* rejects none of them, which keeps the source's actual tolerance while making the verdict available; the rustdoc states that the source warns and loads anyway. Covered by both_fixtures_ibd_headers_match_their_declared_uuid and a_mismatched_or_missing_uuid_is_reported_not_rejected.

## CPP-128 — The declared .ibd checksums are parsed and mirrored but never verified on any read path

**Affected files:** src/openms/source/FORMAT/ImzMLFile.cpp:148-160 (attachImzMLMeta_); src/openms/source/FORMAT/HANDLERS/ImzMLHandler.cpp:703-704

**Issue and reproduction:** IMS:1000091 (SHA-1) and IMS:1000090 (MD5) are parsed into ImzMLMeta and copied onto the loaded MSExperiment as imzml:ibd_sha1 / imzml:ibd_md5 MetaValues, and no read path recomputes either digest. The UUID header check covers only the first 16 bytes, so a .ibd that is truncated after the header, or whose array region was corrupted in transfer, decodes into plausible-looking numbers with no indication that the file no longer matches what the .imzML describes. imzML declares these checksums precisely because the two files travel separately.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Add an opt-in verification — a PeakFileOptions flag or an explicit ImzMLFile method — that streams the .ibd through the declared digest. It is one sequential pass and it should not be automatic, but it should be reachable.

**Rust handling:** ImzMLHandler::verify_ibd_sha1 re-hashes the .ibd in 64 KiB chunks and returns ChecksumStatus, bounded by max_checksum_bytes and never automatic. The upstream processed fixture's declared digest 7e8fdb93053915d3edb51b70aa0619ac209964df is reproduced in the_declared_ibd_sha1_is_reproducible_and_md5_is_only_parsed. MD5 is parsed only, because this crate has no MD5 implementation and adding a dependency is a Cargo.toml change outside this package's scope; ImzMLMeta::ibd_md5 preserves the declared string for a later stage.

## CPP-129 — OnDiscImzMLExperiment::open re-derives the .ibd path in a branch that cannot be taken

**Affected files:** src/openms/source/KERNEL/OnDiscImzMLExperiment.cpp:238-246; src/openms/source/FORMAT/ImzMLFile.cpp:568-570

**Issue and reproduction:** open() sets ibd_path_ from the override or from meta_.ibd_file_path, then re-implements ImzMLFile::inferIbdPath_ inline for the case where that string is empty. It is never empty: loadImpl_ resolves the path as override-or-inferIbdPath_ and assigns meta_.ibd_file_path unconditionally before parsing. The inline copy is unreachable duplicated logic; a future change to inferIbdPath_ would leave it silently behind.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Delete the `if (pimpl_->ibd_path_.empty())` block and keep `ibd_path_ = pimpl_->meta_.ibd_file_path`, which loadImpl_ already guarantees to be the file that was read.

**Rust handling:** One derivation only: `infer_ibd_path` in src/format/imzml_handler.rs, called by `OnDiscImzMLExperiment::open`; `open_with_ibd`/`open_with_limits` take the path explicitly. Documented on `open` and in docs/ON_DISC_IMZML_SUPPORT.md.

## CPP-130 — IonImage allocates width * height from unvalidated file-supplied image dimensions

**Affected files:** src/openms/source/IMAGING/IonImage.cpp:21-28; src/openms/include/OpenMS/IMAGING/IonImageExtraction.h:89; src/openms/source/FORMAT/ImzMLFile.cpp:380-385

**Issue and reproduction:** extractIonImage constructs IonImage(geom.getWidth(), geom.getHeight()), and resize() assigns width*height doubles plus a parallel std::vector<bool> with no ceiling. Both dimensions come from the imzML's IMS:1000042 / IMS:1000043 (raised to the largest observed coordinate), so a corrupt or hostile header sizes the allocation directly: a declared 100000 x 100000 grid asks for 80 GB of doubles before a single peak is read. The dataset need contain only one spectrum.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Bound the pixel product in MSImagingGeometry::setDimensions or IonImage::resize, e.g. against a MAX_IMAGE_PIXELS constant, and throw Exception::InvalidValue rather than attempting the allocation.

**Rust handling:** MAX_IMAGE_PIXELS = 16,777,216 (4096 x 4096) is checked in `ImagingGeometry::set_dimensions`, `add_pixel`, `IonImage::new`/`resize` and `ImagingRegion::from_mask`, and the allocation uses try_reserve_exact. A file declaring a larger grid fails `open()` and commits nothing; tested by `a_declared_grid_above_the_ceiling_is_refused` and `an_image_above_the_pixel_ceiling_is_refused`.

## CPP-131 — MSImagingRegion::fromMask computes the far corner in wrapping UInt arithmetic

**Affected files:** src/openms/source/IMAGING/MSImagingRegion.cpp:52-53

**Issue and reproduction:** max_x_ = origin_x + width - 1 and max_y_ = origin_y + height - 1 are evaluated in UInt. For an origin near the top of the range the sum wraps, producing a bounding box whose maximum is below its minimum — the exact state rectangle() rejects. contains() then answers false everywhere, intersects() reports no overlap with anything, and area() (bbox path) is not reached only because the shape is Mask; getBBoxWidth() underflows instead. The region is silently inert rather than rejected.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Check the sum for overflow before assigning, and throw Exception::InvalidValue as the sibling validations do.

**Rust handling:** `ImagingRegion::from_mask` computes the corner with `checked_add` and returns `Error::InvalidRange` on overflow; tested by the `u32::MAX` origin case in `a_mask_covers_only_its_set_bits`.

## CPP-132 — ImzMLWriter::store cannot write a metadata-only continuous dataset

**Affected files:** src/openms/source/FORMAT/HANDLERS/ImzMLWriter.cpp:485-491 (applyStoreOptions_ clear(false)), :692-754 (spectraShareMz_ / isContinuousMode_), :1410-1415 (store)

**Issue and reproduction:** PeakFileOptions::setMetadataOnly(true) makes applyStoreOptions_ clear every peak with spectrum.clear(false). isContinuousMode_ then calls spectraShareMz_ over the emptied spectra, which cannot find a non-empty reference and returns false, so the explicit "continuous" branch throws Exception::InvalidParameter. Since loading any continuous imzML sets imzml:imaging_mode = "continuous" on the experiment, a metadata-only store of a continuous dataset is impossible, while the same store of a processed dataset succeeds. Source review only; not reproduced against running C++.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Skip the shared-axis check when every spectrum is empty, or honour the declared mode for a metadata-only store: with no peaks the two layouts are indistinguishable on disk.

**Rust handling:** Reproduced deliberately rather than quietly downgrading the declared mode, and pinned by tests/imzml_writer.rs::metadata_only_over_a_declared_continuous_dataset_is_refused_as_in_the_source. Documented as a # Warning on apply_store_options and referenced from store_with_options; removing imzml:imaging_mode lets the same store succeed as processed.

## CPP-133 — ImzMLWriter::store leaks the ProgressLogger recursion depth on every throw

**Affected files:** src/openms/source/FORMAT/HANDLERS/ImzMLWriter.cpp:1437 (startProgress), :1519 (endProgress)

**Issue and reproduction:** logger.startProgress(0, work.size() + 2, ...) is called at line 1437 and logger.endProgress() only on the success path at line 1519. Every throw in between -- UnableToCreateFile on either path, ParseError from the array writers, the UUID header write or the fflush -- leaves ProgressLogger's recursion depth incremented, so later progress output in the same process is indented wrongly and, at the depth cap, suppressed. Source review only.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Guard the range with a scope object whose destructor calls endProgress, or wrap the body in try/catch and end the range before rethrowing.

**Rust handling:** store_with_options ends the progress range before propagating, so the nesting depth is balanced on every path; stated in the # Arguments prose for the logger parameter.

## CPP-134 — A failed ImzMLWriter::store leaves a truncated .ibd with no .imzML

**Affected files:** src/openms/source/FORMAT/HANDLERS/ImzMLWriter.cpp:1440-1511 (the fopen'd .ibd block), :1517 (the .imzML written last)

**Issue and reproduction:** The .ibd is created and streamed before the .imzML, and before several failures can still occur: a ParseError from writeFloat32Array's MAX_IBD_ARRAY_ELEMENTS check, a short fwrite, a failed fflush, or an unwritable .imzML path. UniqueFile_'s destructor only closes the handle, so a partially written .ibd outlives the failure and a later run can mistake it for the companion of an older .imzML. Source review only.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Write to a temporary path and rename on success, or unlink the .ibd when store exits by exception.

**Rust handling:** Every check that does not need the filesystem runs in a preflight before either file is created, so no rejected experiment leaves a file behind; tests/imzml_writer.rs asserts the absence of both files on every refusal path, including all eleven ceilings. A genuine I/O failure mid-write still leaves the partial files, as in the source, and that is stated in the # Errors section.

## CPP-135 — ImzMLFile_1_Example_Continuous.imzML is not schema-valid mzML 1.1.0, so ImzMLFile::isValid would reject its own reference fixture

**Affected files:** src/tests/class_tests/openms/data/ImzMLFile_1_Example_Continuous.imzML (sha256 a358a80751cfd014dc86d5a18ee04c21b98063d3523c0f03f29c375f49c6b0cc); consumed by src/tests/class_tests/openms/source/ImzMLFile_test.cpp and ImzMLFile_all_modes_test.cpp; validated by ImzMLFile::isValid -> Internal::XMLFile::isValid against resources/schemas/mzML_1_10.xsd

**Issue and reproduction:** Validating the fixture against the mzML 1.1.0 schema that ImzMLFile's own constructor registers produces eleven errors: ten cvParam elements omit the required cvRef attribute, and scanSettingsList (line 35) precedes softwareList (line 57) while the schema sequence is softwareList, scanSettingsList, instrumentConfigurationList, dataProcessingList. The same misordering makes a single-pass mzML reader see the dataProcessing softwareRef="sw1" as a forward reference. The C++ suite never notices because its only isValid section validates a file it has just written, and Xerces validation is disabled on the read path (fgSAX2CoreValidation false, ImzMLFile.cpp:576).

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Add cvRef="MS" / cvRef="IMS" to the ten cvParam elements that lack it, and move softwareList before scanSettingsList so the top-level element order matches mzML_1_10.xsd. Optionally add a class-test section that calls isValid on the unmodified fixture, which would have caught both.

**Rust handling:** Not a blocker: openms::format::imzml_file and imzml_handler scan only the IMS vocabulary and read the fixture correctly, and all ported sections pass against it unmodified. ImzMLFile::is_valid reports the eleven diagnostics faithfully, and the misordering is one of the three reasons format::mzml cannot currently read the fixture (recorded in docs/IMZML_FILE_SUPPORT.md under 'the mzML metadata gap'). The fixture was not modified.

## CPP-136 — ImzMLFile_2_Example_Processed.imzML declares mzML version="1.1" instead of 1.1.0 and is ISO-8859-1 with Latin-1 bytes

**Affected files:** src/tests/class_tests/openms/data/ImzMLFile_2_Example_Processed.imzML (sha256 066671590f01d4820b7ac5bd1f669d23dca04a91e40c08a35f7d0c6d6c729480), line 2 root attribute version="1.1", XML declaration encoding="ISO-8859-1"; consumed by ImzMLFile_test.cpp and ImzMLFile_all_modes_test.cpp

**Issue and reproduction:** The root mzML element declares version="1.1" where the schema and every other OpenMS fixture use the three-component 1.1.0, and the document is ISO-8859-1 carrying three Latin-1 bytes (0xFC for u-umlaut in an MS:1000590 contact affiliation, 0xDF twice for sharp-s in the contact address). A strict version check refuses the file; a UTF-8-only reader either refuses it or silently mis-decodes the three bytes into replacement characters, corrupting contact metadata. C++ tolerates both because MzMLHandler warns rather than failing on an unexpected version and Xerces transcodes ISO-8859-1, so the defect is invisible upstream.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Change version="1.1" to version="1.1.0" and re-save the file as UTF-8 with encoding="utf-8" in the declaration, preserving the three characters. Both are edits to the fixture, not to the library.

**Rust handling:** Not a blocker for imzML: openms::format::imzml_handler accepts the non-UTF-8 declared encoding for the ASCII-only IMS values it reads and rejects a non-UTF-8 byte only inside a value it actually consumes, so every ported section passes against the fixture unmodified. It does block two other paths that are outside this package: format::mzml rejects the file as 'only mzML 1.1 is supported', and mzml_schema refuses it as 'requires UTF-8, UTF-16 or ASCII-compatible bytes'. The fixture was not modified.

## CPP-137 — ImzMLFile's two buildImagingGeometry overloads, documented as the same source of truth, disagree on the pixel-size condition and on the skip-reason ordering

**Affected files:** src/openms/source/FORMAT/ImzMLFile.cpp:240-369 (the MSExperiment overload) and :371-476 (the index overload); specifically :277 vs :397 and :363 vs :472. Both are declared static public in src/openms/include/OpenMS/FORMAT/ImzMLFile.h

**Issue and reproduction:** The header calls the index overload 'the source-of-truth path' and says the MetaValue path exists for experiments already loaded, implying they agree. Two differences: (1) the MetaValue overload copies the pixel size when imzml:pixel_size_x and imzml:pixel_size_y merely exist (:363) where the index overload copies it only when both are strictly positive (:472), and MSImagingGeometry::setPixelSize validates nothing (MSImagingGeometry.cpp:28-33), so for a dataset declaring IMS:1000046/047 as 0 the two builders return geometries whose getPixelSizeX() differs (0 versus the default) and any physical-distance computation differs with it; (2) the MetaValue overload tests x<1||y<1 before reading imzml:z (:277) while the index overload tests z!=1 first (:397), so a spectrum at (0,0,2) is warned about as a non-conformant coordinate on one path and skipped silently on the other. The second is cosmetic; the first changes returned data.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Make :363 require pixel_size_x > 0 && pixel_size_y > 0, matching :472, and move the z!=1 test ahead of the coordinate test at :277 so both overloads report the same reason. Alternatively have the MetaValue overload build an ImzMLSpectrumIndex vector and delegate to the index overload, which removes the divergence by construction.

**Rust handling:** Both differences are reproduced deliberately, because they are observable, and both are documented at build_imaging_geometry_from_experiment and in docs/IMZML_FILE_SUPPORT.md. The ordering difference is asserted in the_meta_value_geometry_builder_reports_every_kind_of_unplaceable_spectrum, which checks that the index builder reports other_plane_count 1 and non_positive_count 0 for the same entry the MetaValue builder reports as non-positive. the_geometry_builders_agree_on_the_upstream_fixtures pins that the two agree on width, height and every pixel for both real fixtures. ImagingGeometry::set_pixel_size rejects only a non-finite size, as the source rejects nothing.

## CPP-138 — ClassTest::isRealSimilar reports any two infinities as similar, including +inf against -inf

**Affected files:** src/testframework/source/CONCEPT/ClassTest.cpp:364-489 (quotient at :438, sign test at :439, ratio test at :467)

**Issue and reproduction:** For two infinite arguments both `absdiff = number_1 - number_2` and `ratio = number_1 / number_2` evaluate to NaN. Every subsequent comparison against NaN is false, so `ratio < 0.`, `ratio < 1.` and `ratio > ratio_max_allowed` all fail to fire and control falls through to the final else, which sets fuzzy_message = "ratio of numbers is small" and returns true. TEST_REAL_SIMILAR(inf, -inf) therefore passes silently.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Reject non-finite arguments explicitly alongside the existing NaN guards at :374-385, e.g. `if (!std::isfinite(number_1) || !std::isfinite(number_2)) { fuzzy_message = "...not finite"; return false; }` — the existing NaN checks already establish the precedent and the message channel.

**Rust handling:** is_real_similar (tests/imzml_file.rs:225-259) reproduces the quirk deliberately, because its contract is to be the C++ oracle, and the_oracle_reproduces_the_upstream_infinity_quirk pins it so it cannot be mistaken for a port bug. close(), the assertion the 38 call sites use, refuses it: it requires both operands to be finite before asserting anything.

## CPP-139 — ClassTest::isRealSimilar is asymmetric — the quotient underflows to -0.0 and defeats its own opposite-sign branch

**Affected files:** src/testframework/source/CONCEPT/ClassTest.cpp:438-451 (ratio = number_1 / number_2; if (ratio < 0.))

**Issue and reproduction:** The opposite-sign test is `ratio < 0.`, read off the quotient rather than the operands. When the magnitudes are far enough apart the quotient underflows to negative zero, and `-0.0 < 0.` is false, so the branch is skipped. `-0.0 < 1.` is then true, so the reciprocal step sets ratio = 1. / -0.0 = -inf, and `-inf > ratio_max_allowed` is false — so the function returns true for two numbers of opposite sign that differ by 600 orders of magnitude. isRealSimilar(1e-300, -1e300) is true while isRealSimilar(-1e300, 1e-300) is false, and isRealSimilar(1.0, -inf) is true while isRealSimilar(-inf, 1.0) is false. The same-sign pair isRealSimilar(1e-300, 1e300) is correctly false, which shows this is a defect and not a deliberate tolerance. Confirmed by compiling and running the transcribed control flow, not by reading it.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Decide the sign case from the operands, not from the quotient: replace `if (ratio < 0.)` with a test such as `if (std::signbit(number_1) != std::signbit(number_2))`, placed before the quotient is taken. Both numbers are already known non-zero at that point, so signbit is unambiguous. That is what the branch's own message ("numbers have different signs") already claims it tests.

**Rust handling:** is_real_similar reproduces the asymmetry faithfully and the_oracle_reproduces_the_upstream_negative_zero_asymmetry pins both directions. close() closes it the way the fix above would: it computes opposite_signs from left.is_sign_negative() != right.is_sign_negative() and requires the absolute difference to be within 1e-5 whenever the signs differ. Both close() guards only ever reject pairs upstream accepted, so no ported assertion is weaker than its C++ original.

## CPP-140 — ImzMLHandler's non-external peak fill is dead code: MzMLHandler throws first on a mixed spectrum

**Affected files:** src/openms/source/FORMAT/HANDLERS/ImzMLHandler.cpp:198-234 and :608-610; src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp:407-413

**Issue and reproduction:** ImzMLInterceptConsumer::consumeSpectrum is written to handle a spectrum with exactly one external peak array: it enters its decode block on `(ims->mz_meta.is_ext || ims->int_meta.is_ext)` (:198), reads the external side from the .ibd, and fills the non-external side from the peaks MzMLHandler decoded — `mz_vec[i] = s[i].getMZ()` (:214-218), `int_vec[i] = s[i].getIntensity()` (:228-232) — so the mismatch throw at :234 should be reachable only for genuine corruption. It never runs. For such a spectrum the inline side's <binary> decodes to N values while the external side's placeholder <binary/> decodes to none, and MzMLHandler::populateSpectraWithData_ compares those two DECODED sizes at :410 and calls fatalError(LOAD, "The length of m/z and integer values of spectrum '...' differ ... Not reading spectrum!"), which throws Exception::ParseError. ImzMLHandler::onEndElement's two mitigations do not help: zeroing default_array_length_ (:608-610) only affects the defaultArrayLength comparison at :415-429, which populateSpectraWithData_ repairs rather than throws on, and zeroing bin_data_.back().size (:550-553) touches the declared length, not the decoded vector. So the fill at :214-232 is latent in the in-memory loader. The on-disc loader cannot reach it either, for an unrelated reason: ImzMLSpectrumIndex carries no is_ext field (ImzMLHandlerHelper.h:102-139), so Impl::readMz_/readInt_ gate on `mz_length`/`int_length != 0` (OnDiscImzMLExperiment.cpp:180-199) and an inline array, having no IMS:1000103, is silently empty — which then trips the on-disc mismatch throw at :75. Not exploitable and not reachable from a conformant imzML 1.1.0, which stores both peak arrays externally; the cost is dead code that reads as a supported path and a confusing error (a base-class message about a spectrum length, for a file whose real problem is a missing IMS:1000101).

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Decide which behaviour is intended and make the code say so. Either (a) keep the fill and let it work: have ImzMLHandler suppress the base class's decoded-size comparison for a mixed spectrum the way it already suppresses the defaultArrayLength warning — e.g. drop the external array's BinaryData entry from bin_data_ instead of only zeroing its `size`, so computeDataProperties_ does not find it and populateSpectraWithData_ returns early at :386-394 — and add an is_ext pair to ImzMLSpectrumIndex so loadSpectraIndex and OnDiscImzMLExperiment can tell an inline array from a zero-length one; or (b) delete the fill at :214-232 and reject a mixed spectrum explicitly in ImzMLHandler with a message that names IMS:1000101, so the diagnosis is not left to a base-class array-length error. A class test for a single-external-array spectrum would have caught this at either end; ImzMLFile_test.cpp has none, because no fixture is non-conformant.

**Rust handling:** The Rust port implements the interception layer, not MzMLHandler, so it follows option (a) as :214-232 describes it: src/format/imzml_handler.rs retains the inline base64 of a non-external peak array during the index scan (ImzMLSpectrumIndex::mz_inline / int_inline) and decodes it in spectrum() and in mz_array/intensity_array where the .ibd read would otherwise be, keeping the length check as the corruption guard. It also carries the is_ext pair the C++ index lacks (mz_external / int_external), which is what lets both Rust read paths agree. The consequence — Rust accepts a non-conformant file both C++ loaders reject, in a different class and with a different message — is stated in docs/IMZML_HANDLER_SUPPORT.md:340-355, in the DecodedSpectrum::inline_peaks and ImzMLHandler::spectrum rustdoc, and as a gap in tests/data/imzml_handler_provenance.json, with source anchors at ImzMLHandler.cpp:198/214/228/608, MzMLHandler.cpp:410 and ImzMLHandlerHelper.h:102.

## CPP-141 — DataValue::operator double() reads the union's double member for a non-numeric value

**Duplicate tracking:** This is the imzML call-site instance of [CPP-058](#cpp-058--nonnumeric-datavalue-casts-read-an-inactive-union-member). The stable ID and original evidence are retained; it is not an additional distinct DataValue defect.

**Affected files:** src/openms/source/DATASTRUCTURES/DataValue.cpp:466-478 (operator double()), :480-491 (operator float()), :452-464 (operator long double()); reached from src/openms/source/FORMAT/HANDLERS/ImzMLWriter.cpp:356-378

**Issue and reproduction:** operator double() throws only for EMPTY_VALUE and converts INT_VALUE; for STRING_VALUE and the three list types it falls through to `return data_.dou_;`. For a STRING_VALUE the live union member is a std::string*, so this reinterprets a pointer's bit pattern as a double — undefined behaviour, and in practice a garbage number rather than the Exception::ConversionError that operator int() and operator unsigned int() raise for the same input. extractMeta_ reaches it for all four of imzml:pixel_size_x, imzml:pixel_size_y, imzml:max_dim_x and imzml:max_dim_y, so `exp.setMetaValue("imzml:pixel_size_x", "wide")` followed by ImzMLFile::store writes an imzML carrying a nonsense pixel size and a nonsense IMS:1000044/45 extent computed from it.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Give operator double(), operator float() and operator long double() the type guard their integer siblings already have: throw Exception::ConversionError unless value_type_ is DOUBLE_VALUE or INT_VALUE. The INT_VALUE conversion and the EMPTY_VALUE throw stay as they are, so no caller that passes a number changes behaviour.

**Rust handling:** Not reproducible in the port — MetaValue is a tagged enum and MetaValue::as_f64 returns Err for every non-numeric variant, so there is no union to misread. dataset_meta therefore returns Error::InvalidValue for a non-numeric size, which is the one place it is deliberately stricter than extractMeta_. Pinned by dataset_meta_reads_every_key_and_refuses_the_wrong_type and the_numeric_keys_stay_strict_and_a_size_accepts_an_integer, the latter also asserting that an integer size is still accepted because the source accepts one. Recorded as cpp_issue_candidate 4 in tests/data/imzml_writer_provenance.json and under Source findings in docs/IMZML_WRITER_SUPPORT.md; unconfirmed against running C++.

## CPP-142 — ImzMLWriter::store treats one misaligned FloatDataArray as skippable or fatal depending on an unrelated PeakFileOptions flag

**Affected files:** src/openms/source/FORMAT/HANDLERS/ImzMLWriter.cpp:454, 468, 475, 522-531; src/openms/source/KERNEL/MSSpectrum.cpp:20-38, 40-66, 444-460; src/openms/include/OpenMS/KERNEL/MSSpectrum.h:364-376

**Issue and reproduction:** appendAndWriteFloatDataArrays_ skips a FloatDataArray whose size() differs from the spectrum's and warns, under an explicit "Per-peak contract" comment — the writer's stated policy. But applyStoreOptions_ runs first and reaches MSSpectrum::sortByPosition (:454 inside the peak-filter branch and :475 in the else branch) and MSSpectrum::select (:468), all of which run checkDataArraySizes_ and throw Exception::Precondition for exactly such an array. So the same experiment with the same array stores with a warning under default options over an already-sorted spectrum, and aborts the whole store as soon as getSortSpectraByMZ() has an unsorted spectrum to sort or an m/z or intensity range actually trims a peak — decided by a flag that has nothing to do with the array.

**Evidence:** Source review at the pinned revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. No C++ execution, sanitizer run or upstream fix is claimed.

**Proposed fix:** Apply the per-peak-contract check before applyStoreOptions_ and drop the offending arrays there, warning once per array as appendAndWriteFloatDataArrays_ does, so skip-and-warn is the only outcome and neither sortByPosition nor select ever sees a misaligned array. The kernel preconditions stay as they are — they are correct to refuse a reorder under a mis-sized annotation array.

**Rust handling:** Was reproduced in the port and is now fixed at this package's own call site, since src/kernel.rs is out of scope and is right to refuse: apply_store_options lifts every misaligned float, integer and string data array off the spectrum around the sort and the peak filters via detach_misaligned_arrays / restore_misaligned_arrays, restoring each at its original index with its original values, so the array and its original length still reach StoreReport::skipped_float_arrays. validate_data_arrays, select and sort_by_position are unchanged. Pinned by a_misaligned_array_is_skipped_whatever_the_peak_file_options_say (three option sets, one outcome) and apply_store_options_puts_a_misaligned_array_back_in_place. Recorded as cpp_issue_candidate 5; unconfirmed against running C++.


## FORMAT wave evidence scope

The following entries were consolidated against the pinned source and existing
IDs during integration. Source review does not imply executed C++ reproduction.
The two unconfirmed candidates are explicitly labeled. Some native APIs retain
source behavior for compatibility; a recorded proposed fix is not an upstream
change or a claim that the Rust behavior was corrected.

## CPP-143 — Negative CHEMMOD identifiers fail modification-cell parsing

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MzTabModification::fromCellString` in [src/openms/source/FORMAT/MzTab.cpp:130](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTab.cpp#L130).

**Trigger:** Parse 8-CHEMMOD:-18.010565.

**Issue:** Splitting on every hyphen yields three fields; the parser requires exactly two and throws.

**Proposed C++ fix:** Split position/identifier once with CHEMMOD signed mass retained. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Preserves source rejection of negative CHEMMOD strings. No native grammar correction is claimed. See [src/format/mztab.rs](src/format/mztab.rs), [tests/mztab.rs](tests/mztab.rs) and [MZTAB_SUPPORT.md](docs/MZTAB_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-144 — Quoted commas split a MzTab modification-list entry

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MzTabModificationList::fromCellString` in [src/openms/source/FORMAT/MzTab.cpp:236](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTab.cpp#L236).

**Trigger:** A modification parameter includes quoted name "blabla, [bla]".

**Issue:** Comma protection excludes in_quotes, so later comma split separates a quoted field.

**Proposed C++ fix:** Protect all commas within parameters/quotes and keep quote/bracket states distinct. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Preserves source quoted-comma splitting behavior, including the defective case. No native tokenizer correction is claimed. See [src/format/mztab.rs](src/format/mztab.rs), [tests/mztab.rs](tests/mztab.rs) and [MZTAB_SUPPORT.md](docs/MZTAB_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-145 — MzTab parameter rendering fails to quote a bare comma

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MzTabParameter::toCellString` in [src/openms/source/FORMAT/MzTabBase.cpp:327](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTabBase.cpp#L327).

**Trigger:** Set parameter name a,b and serialize/reparse.

**Issue:** Only comma-space triggers quoting, producing an extra cell component for bare comma.

**Proposed C++ fix:** Quote any comma and escape embedded quotes consistently. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Preserves source comma-space-only quoting rule. Bare comma rendering remains source-compatible and non-round-trippable. See [src/format/mztab.rs](src/format/mztab.rs), [tests/mztab.rs](tests/mztab.rs) and [MZTAB_SUPPORT.md](docs/MZTAB_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-146 — MzTab score-by-run header order differs from row order

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MzTabFile::generateMzTabProteinHeader_` in [src/openms/source/FORMAT/MzTabFile.cpp:2033](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTabFile.cpp#L2033); `MzTabFile::generateMzTabSectionRow_ (protein overload)` in [src/openms/source/FORMAT/MzTabFile.cpp:2114](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTabFile.cpp#L2114).

**Trigger:** Protein row has scores1/2 and runs1/2 with four distinct values.

**Issue:** Header iterates run then score; row iterates score then run, interchanging off-diagonal cells.

**Proposed C++ fix:** Share score/run ordering between header and row. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native SectionLayout aligns header and cells; checks source fixtures and synthetic column families. See [src/format/mztab_file.rs](src/format/mztab_file.rs), [tests/mztab_file.rs](tests/mztab_file.rs) and [MZTAB_FILE_SUPPORT.md](docs/MZTAB_FILE_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-147 — MzTab PSM optional columns are dropped on load

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MzTabFile::load` in [src/openms/source/FORMAT/MzTabFile.cpp:1278](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTabFile.cpp#L1278).

**Trigger:** PSH declares opt_global_score with a PSM value.

**Issue:** Header matches exactly opt_ rather than its prefix, so named optional columns never get registered.

**Proposed C++ fix:** Use hasPrefix(cells[i], "opt_"). No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native reader recognizes named optional columns. See [src/format/mztab_file.rs](src/format/mztab_file.rs), [tests/mztab_file.rs](tests/mztab_file.rs) and [MZTAB_FILE_SUPPORT.md](docs/MZTAB_FILE_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-148 — MzTab column-unit metadata parses its key as an index

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MzTabFile::load` in [src/openms/source/FORMAT/MzTabFile.cpp:683](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTabFile.cpp#L683).

**Trigger:** Load MTD colunit-protein with a valid value.

**Issue:** extractBracketIndex receives colunit without brackets; downstream indexed vector assignment has no resize.

**Proposed C++ fix:** Parse nonindexed column-unit value into vector append; validate malformed keys. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native column-unit grammar uses the expected vector form. See [src/format/mztab_file.rs](src/format/mztab_file.rs), [tests/mztab_file.rs](tests/mztab_file.rs) and [MZTAB_FILE_SUPPORT.md](docs/MZTAB_FILE_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-149 — MzTab-M assay custom metadata is written under ms_run

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MzTabMFile::generateMzTabMMetaDataSection_` in [src/openms/source/FORMAT/MzTabMFile.cpp:222](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTabMFile.cpp#L222).

**Trigger:** Assay1 has custom1 parameter while ms_run1 denotes a different entity.

**Issue:** Writer emits ms_run[1]-custom[1], attributing assay metadata to run.

**Proposed C++ fix:** Emit assay index prefix. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Correct by default; source_assay_custom_key opts into source spelling. See [src/format/mztab_m.rs](src/format/mztab_m.rs), [tests/mztab_m.rs](tests/mztab_m.rs) and [MZTAB_M_SUPPORT.md](docs/MZTAB_M_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-150 — MzTab-M column-unit families share an incorrect output key

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MzTabMFile::generateMzTabMMetaDataSection_` in [src/openms/source/FORMAT/MzTabMFile.cpp:347](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTabMFile.cpp#L347).

**Trigger:** Metadata has different feature and evidence column units.

**Issue:** All loops emit colunit_small_molecule, losing family distinction.

**Proposed C++ fix:** Use distinct feature/evidence keys. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Correct by default; source_colunit_keys explicit compatibility flag. See [src/format/mztab_m.rs](src/format/mztab_m.rs), [tests/mztab_m.rs](tests/mztab_m.rs) and [MZTAB_M_SUPPORT.md](docs/MZTAB_M_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-151 — MzTab-M exporter changes metadata keys before lookup

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MzTabM::getFeatureMapMetaValues_` in [src/openms/source/FORMAT/MzTabM.cpp:103](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MzTabM.cpp#L103).

**Trigger:** Feature metadata contains raw key with space.

**Issue:** Collector changes key to underscore spelling, so later original metadata lookup yields absence/null.

**Proposed C++ fix:** Retain raw lookup keys; normalize output column names only. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** substitute_keys_before_lookup default false, preserving raw key lookup. See [src/format/mztab_m.rs](src/format/mztab_m.rs), [tests/mztab_m.rs](tests/mztab_m.rs) and [MZTAB_M_SUPPORT.md](docs/MZTAB_M_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-152 — pepXML end_scan mismatch check reads start_scan twice

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `PepXMLFile::readRTMZCharge_` in [src/openms/source/FORMAT/PepXMLFile.cpp:926](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/PepXMLFile.cpp#L926).

**Trigger:** spectrum_query start_scan=1 end_scan=2.

**Issue:** Both locals read start_scan so merged-scan diagnostic never fires.

**Proposed C++ fix:** Read end_scan for endscan. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native parser compares actual attributes and reports mismatch. See [src/format/pepxml.rs](src/format/pepxml.rs), [tests/pepxml.rs](tests/pepxml.rs) and [PEPXML_SUPPORT.md](docs/PEPXML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-153 — pepXML drops a uniquely resolved undeclared modification

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `PepXMLFile::onStartElement` in [src/openms/source/FORMAT/PepXMLFile.cpp:1725](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/PepXMLFile.cpp#L1725).

**Trigger:** A mod_aminoacid_mass not in header resolves to exactly one modification in registry.

**Issue:** mods nonempty branch appends only if size>1; exactly one match silently omitted.

**Proposed C++ fix:** Always append first resolved match; warn only if multiple. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native reader installs uniquely resolved match. See [src/format/pepxml.rs](src/format/pepxml.rs), [tests/pepxml.rs](tests/pepxml.rs) and [PEPXML_SUPPORT.md](docs/PEPXML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-154 — pepXML fixed protein C-terminal modification misses terminal branch

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `PepXMLFile::onEndElement` in [src/openms/source/FORMAT/PepXMLFile.cpp:2123](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/PepXMLFile.cpp#L2123).

**Trigger:** Implicit fixed modification has PROTEIN_C_TERM specificity.

**Issue:** C-terminal conditional repeats PROTEIN_N_TERM and falls through to internal-residue loop.

**Proposed C++ fix:** Replace second repeated enum with PROTEIN_C_TERM. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native terminal logic accounts for protein C terminus; source-specific warnings documented. See [src/format/pepxml.rs](src/format/pepxml.rs), [tests/pepxml.rs](tests/pepxml.rs) and [PEPXML_SUPPORT.md](docs/PEPXML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-155 — qcML loses units across its own store/load

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `QcMLFile::QualityParameter::toXMLString` in [src/openms/source/FORMAT/QcMLFile.cpp:81](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/QcMLFile.cpp#L81); `QcMLFile::onStartElement` in [src/openms/source/FORMAT/QcMLFile.cpp:820](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/QcMLFile.cpp#L820).

**Trigger:** QualityParameter has nonempty unitRef and unitAcc.

**Issue:** Writer uses unitRef/unitAcc; reader recognizes unitCvRef/unitAccession.

**Proposed C++ fix:** Use canonical names and optional legacy aliases. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Reader accepts both legacy and schema unit spellings. Default writer uses schema spelling; WriteOptions::source() selects legacy output. See [src/format/qcml.rs](src/format/qcml.rs), [tests/qcml.rs](tests/qcml.rs) and [QCML_SUPPORT.md](docs/QCML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-156 — qcML table writer discards its normalized row copy

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `QcMLFile::Attachment::toXMLString` in [src/openms/source/FORMAT/QcMLFile.cpp:227](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/QcMLFile.cpp#L227).

**Trigger:** One table cell contains a space.

**Issue:** Writer substitutes spaces in copy_row, then concatenates original row, splitting cell on reload.

**Proposed C++ fix:** Serialize copy_row, or implement an escaping grammar. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Default writer refuses empty or XML-whitespace table cells. source_table_text reproduces unnormalized source rows; it does not normalize them into a corrected value. See [src/format/qcml.rs](src/format/qcml.rs), [tests/qcml.rs](tests/qcml.rs) and [QCML_SUPPORT.md](docs/QCML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-157 — qcML removeAllAttachments omits set-only entries

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `QcMLFile::removeAllAttachments` in [src/openms/source/FORMAT/QcMLFile.cpp:511](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/QcMLFile.cpp#L511).

**Trigger:** Set has matching attachment and no run attachment map entry of same ID.

**Issue:** Method iterates runQualityAts_ only despite all-runs/sets contract.

**Proposed C++ fix:** Visit union of run and set identifiers. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Preserves run-map-only iteration. Sets are reached only through a matching run attachment-map key. See [src/format/qcml.rs](src/format/qcml.rs), [tests/qcml.rs](tests/qcml.rs) and [QCML_SUPPORT.md](docs/QCML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-158 — qcML map2csv emits misaligned rows when a column is missing

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `QcMLFile::map2csv` in [src/openms/source/FORMAT/QcMLFile.cpp:684](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/QcMLFile.cpp#L684).

**Trigger:** First map row has columns A/B; next has only B.

**Issue:** Missing A emits no separator/cell, so B shifts into A column.

**Proposed C++ fix:** Emit empty missing cells and choose union of columns, or reject inconsistent rows. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Preserves first-row column selection and missing-cell misalignment. Resource limits do not correct this layout behavior. See [src/format/qcml.rs](src/format/qcml.rs), [tests/qcml.rs](tests/qcml.rs) and [QCML_SUPPORT.md](docs/QCML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-159 — qcML writer and reader disagree on set-member CV accession

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `QcMLFile::onEndElement` in [src/openms/source/FORMAT/QcMLFile.cpp:940](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/QcMLFile.cpp#L940); `QcMLFile::store` in [src/openms/source/FORMAT/QcMLFile.cpp:2072](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/QcMLFile.cpp#L2072).

**Trigger:** Register set membership to a run with its name parameter and store/reload.

**Issue:** Writer emits QC:0000005 membership record; reader treats MS:1000577 as set member.

**Proposed C++ fix:** Use one membership representation and resolve names/IDs consistently. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Reader also recognizes QC:0000005 membership while retaining the parameter. Writer checks representability; source policy can select original drops. See [src/format/qcml.rs](src/format/qcml.rs), [tests/qcml.rs](tests/qcml.rs) and [QCML_SUPPORT.md](docs/QCML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-160 — qcML TIC slump percentage truncates before multiplication

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `QcMLFile::collectQCData` in [src/openms/source/FORMAT/QcMLFile.cpp:1293](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/QcMLFile.cpp#L1293).

**Trigger:** Run has200 spectra,100 below threshold.

**Issue:** (100/200)*100 integer arithmetic yields0 instead of50.

**Proposed C++ fix:** Compute floating percentage 100.0*count/size with explicit empty guard; same RIC formula. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** collectQCData remains unported; avoid reproducing when QC-computation wave begins. See [src/format/qcml.rs](src/format/qcml.rs), [tests/qcml.rs](tests/qcml.rs) and [QCML_SUPPORT.md](docs/QCML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-161 — Percolator loader requires FileName through an unchecked map lookup

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `PercolatorInfile::load` in [src/openms/source/FORMAT/PercolatorInfile.cpp:239](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/PercolatorInfile.cpp#L239).

**Trigger:** Valid rectangular PIN omits optional FileName column.

**Issue:** Filename map stays empty; .at(UNKNOWN) throws std::out_of_range rather than declared parse error. Caller may catch; process abort is NOT inevitable.

**Proposed C++ fix:** Register default filename or reject missing requirement with structured ParseError. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native Error::MissingInformation names absent column. See [src/format/percolator_infile.rs](src/format/percolator_infile.rs), [tests/percolator_infile.rs](tests/percolator_infile.rs) and [PERCOLATOR_INFILE_SUPPORT.md](docs/PERCOLATOR_INFILE_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-162 — Percolator enzyme features use unmapped protein-terminal markers

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `PercolatorInfile::stampPinFeaturesOnHits (feature tracking overload)` in [src/openms/source/FORMAT/PercolatorInfile.cpp:520](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/PercolatorInfile.cpp#L520).

**Trigger:** PeptideEvidence has [ before or ] after at protein terminus.

**Issue:** enzN/enzC are computed before terminal marker becomes recognized dash.

**Proposed C++ fix:** Normalize markers before both enzyme features and peptide text. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Source feature behavior deliberately reproduced and documented; source correction not applied upstream. See [src/format/percolator_infile.rs](src/format/percolator_infile.rs), [tests/percolator_infile.rs](tests/percolator_infile.rs) and [PERCOLATOR_INFILE_SUPPORT.md](docs/PERCOLATOR_INFILE_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-163 — Percolator writer and loader disagree on trailing protein-list width

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `PercolatorInfile::load` in [src/openms/source/FORMAT/PercolatorInfile.cpp:234](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/PercolatorInfile.cpp#L234); `PercolatorInfile::stampPinFeaturesOnHits (feature tracking overload)` in [src/openms/source/FORMAT/PercolatorInfile.cpp:544](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/PercolatorInfile.cpp#L544).

**Trigger:** One PSM maps to two protein accessions.

**Issue:** Writer separates proteins with tab; loader rejects row wider than header.

**Proposed C++ fix:** Read final Proteins column as variable-width list, retaining semicolon support for Sage. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native writer preserves tabs and reader explicitly reports source incompatibility. See [src/format/percolator_infile.rs](src/format/percolator_infile.rs), [tests/percolator_infile.rs](tests/percolator_infile.rs) and [PERCOLATOR_INFILE_SUPPORT.md](docs/PERCOLATOR_INFILE_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-164 — Mascot query index guard accepts one-past-end

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::MascotXMLHandler::onStartElement` in [src/openms/source/FORMAT/HANDLERS/MascotXMLHandler.cpp:49](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MascotXMLHandler.cpp#L49).

**Trigger:** NumQueries=1 and peptide query=2, followed by a peptide field.

**Issue:** Zero-based index1 is not greater than size1 and later indexes outside vector.

**Proposed C++ fix:** Check index>=size and validate positive query before subtraction. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native lazy query indexing checks declared range before access. See [src/format/mascot_xml.rs](src/format/mascot_xml.rs), [tests/mascot_xml.rs](tests/mascot_xml.rs) and [MASCOT_XML_SUPPORT.md](docs/MASCOT_XML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-165 — Mascot RT failure test treats zero as failure and NaN as success

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::MascotXMLHandler::onEndElement` in [src/openms/source/FORMAT/HANDLERS/MascotXMLHandler.cpp:108](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MascotXMLHandler.cpp#L108).

**Trigger:** Title lookup leaves RT NaN, or legitimately returns RT0.

**Issue:** Boolean negation of double skips missing NaN and flags valid zero.

**Proposed C++ fix:** Use finite/NaN or explicit found flag. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native lookup reports absent information explicitly; zero RT retained. See [src/format/mascot_xml.rs](src/format/mascot_xml.rs), [tests/mascot_xml.rs](tests/mascot_xml.rs) and [MASCOT_XML_SUPPORT.md](docs/MASCOT_XML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-166 — Mascot MGF loader carries precursor and RT fields between blocks

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MascotGenericFile::load (template)` in [src/openms/include/OpenMS/FORMAT/MascotGenericFile.h:85](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/include/OpenMS/FORMAT/MascotGenericFile.h#L85); `MascotGenericFile::load (template)` in [src/openms/include/OpenMS/FORMAT/MascotGenericFile.h:140](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/include/OpenMS/FORMAT/MascotGenericFile.h#L140).

**Trigger:** First BEGIN IONS block has RTINSECONDS; second block has its own PEPMASS but omits optional RTINSECONDS.

**Issue:** The second spectrum retains the previous RT because only selected state is cleared. Missing PEPMASS similarly retains precursor state, but optional RT omission is sufficient for this trigger.

**Proposed C++ fix:** Reset per-query spectrum metadata or construct a fresh spectrum each block. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** The default CarryOver::Reset starts each MGF block with fresh metadata. CarryOver::Source explicitly opts into the source carry-over behavior; this differs from default general MGF input. See [src/format/mascot_generic.rs](src/format/mascot_generic.rs), [tests/mascot_generic.rs](tests/mascot_generic.rs) and [MASCOT_GENERIC_SUPPORT.md](docs/MASCOT_GENERIC_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-167 — mzIdentML writer places C-terminal modification at last-residue location

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::MzIdentMLHandler::writePeptideHit` in [src/openms/source/FORMAT/HANDLERS/MzIdentMLHandler.cpp:1363](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MzIdentMLHandler.cpp#L1363).

**Trigger:** A peptide of lengthN has a C-terminal modification.

**Issue:** Writer emits location N; schema convention is N+1, so it targets final residue.

**Proposed C++ fix:** Emit sequence.size()+1 and test terminal round trip. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native writer uses terminal location convention and tests modified peptides. See [src/format/mzidentml.rs](src/format/mzidentml.rs), [tests/mzidentml.rs](tests/mzidentml.rs) and [MZIDENTML_SUPPORT.md](docs/MZIDENTML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-168 — mzIdentML reader dereferences missing PeptideSequence child

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::MzIdentMLDOMHandler::parsePeptideSiblings_` in [src/openms/source/FORMAT/HANDLERS/MzIdentMLDOMHandler.cpp:2466](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MzIdentMLDOMHandler.cpp#L2466).

**Trigger:** A Peptide includes empty PeptideSequence element.

**Issue:** getFirstChild returns null and getNodeType dereferences it.

**Proposed C++ fix:** Handle empty text node explicitly, reject or represent empty sequence according to schema. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native parser returns checked outcome; no null dereference. See [src/format/mzidentml.rs](src/format/mzidentml.rs), [tests/mzidentml.rs](tests/mzidentml.rs) and [MZIDENTML_SUPPORT.md](docs/MZIDENTML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-169 — mzIdentML substitution position is used as unchecked string index

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::MzIdentMLDOMHandler::parsePeptideSiblings_` in [src/openms/source/FORMAT/HANDLERS/MzIdentMLDOMHandler.cpp:2491](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MzIdentMLDOMHandler.cpp#L2491).

**Trigger:** SubstitutionModification location0 or >peptide length with replacementResidue present.

**Issue:** Signed location-1 converts into string index and writes out of bounds.

**Proposed C++ fix:** Validate1<=location<=sequence length and nonempty residue attributes before indexing. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native parser checks positions and rejects invalid substitutions. See [src/format/mzidentml.rs](src/format/mzidentml.rs), [tests/mzidentml.rs](tests/mzidentml.rs) and [MZIDENTML_SUPPORT.md](docs/MZIDENTML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-170 — mzData checks missing and short arrays after unsafe indexing

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::MzDataHandler::fillData_` in [src/openms/source/FORMAT/HANDLERS/MzDataHandler.cpp:518](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MzDataHandler.cpp#L518).

**Trigger:** Spectrum has one array, or intensity/auxiliary array shorter than m/z.

**Issue:** precisions0/1 accessed before missing-array guard; length mismatch only logs, then peak loop indexes short data.

**Proposed C++ fix:** Require primary arrays before reading precision and validate all decoded lengths before peak construction. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native parser checks primary/aligned-array shapes with structured error. See [src/format/mzdata.rs](src/format/mzdata.rs), [tests/mzdata.rs](tests/mzdata.rs) and [MZDATA_SUPPORT.md](docs/MZDATA_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-171 — mzData writer emits scan modes its reader does not recognize

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::MzDataHandler::writeTo` in [src/openms/source/FORMAT/HANDLERS/MzDataHandler.cpp:919](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MzDataHandler.cpp#L919); `Internal::MzDataHandler::cvParam_` in [src/openms/source/FORMAT/HANDLERS/MzDataHandler.cpp:1088](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MzDataHandler.cpp#L1088).

**Trigger:** Store ABSORPTION/EMC/TDF mode then load.

**Issue:** Writer spellings absent from reader map; fallback loses mode and may modify previous spectrum for MSlevel>=2.

**Proposed C++ fix:** Use shared bidirectional scan-mode table and update current spec in fallback. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Default writer refuses scan modes that cannot round-trip through the source-compatible mapping. Parser updates the current spectrum, not the previous spectrum, for unknown MSn modes. It does not add recognition of these three spellings. See [src/format/mzdata.rs](src/format/mzdata.rs), [tests/mzdata.rs](tests/mzdata.rs) and [MZDATA_SUPPORT.md](docs/MZDATA_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-172 — Streaming mzML consumer references header entries declared only for first record

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MSDataWritingConsumer::consumeSpectrum` in [src/openms/source/FORMAT/DATAACCESS/MSDataWritingConsumer.cpp:62](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/DATAACCESS/MSDataWritingConsumer.cpp#L62); `Internal::MzMLHandler::writeSpectrum_` in [src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp:5250](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L5250).

**Trigger:** Second streamed spectrum has non-default source file or processing history differing from first.

**Issue:** Header derives from dummy one-record map; second sourceFileRef/processing fallback points to undeclared entry. This is independent of numeric-overload claims.

**Proposed C++ fix:** Predeclare record dependencies or refuse new dependencies after header publication. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native consumer rejects dependencies absent from first header rather than emitting dangling references. See [src/format/ms_data_writing_consumer.rs](src/format/ms_data_writing_consumer.rs), [tests/ms_data_writing_consumer.rs](tests/ms_data_writing_consumer.rs) and [MS_DATA_WRITING_CONSUMER_SUPPORT.md](docs/MS_DATA_WRITING_CONSUMER_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-173 — mzXML release decode reads beyond short peak payload

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::MzXMLHandler::doPopulateSpectraWithData_` in [src/openms/source/FORMAT/HANDLERS/MzXMLHandler.cpp:1162](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MzXMLHandler.cpp#L1162).

**Trigger:** peaksCount2 with decoded payload containing only one mz/intensity pair, assertions disabled.

**Issue:** Only assert validates length; loop trusts declared count and indexes past vector.

**Proposed C++ fix:** Replace assert with runtime payload-length validation before loop. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native runtime shape check rejects mismatch in all builds. See [src/format/mzxml.rs](src/format/mzxml.rs), [tests/mzxml.rs](tests/mzxml.rs) and [MZXML_SUPPORT.md](docs/MZXML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-174 — mzXML precursor value and window width depend on SAX chunking

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::MzXMLHandler::onCharacters` in [src/openms/source/FORMAT/HANDLERS/MzXMLHandler.cpp:569](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/MzXMLHandler.cpp#L569).

**Trigger:** Precursor numeric text is delivered in two character callbacks, e.g. an entity/comment boundary.

**Issue:** Each chunk sets m/z and halves previous width again; split comments retain only last chunk.

**Proposed C++ fix:** Accumulate per-element text and parse/apply once at closing tag. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native parser accumulates text and preserves full value/window. See [src/format/mzxml.rs](src/format/mzxml.rs), [tests/mzxml.rs](tests/mzxml.rs) and [MZXML_SUPPORT.md](docs/MZXML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-175 — SVOutStream probe stream remains poisoned after non-newline manipulator

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `SVOutStream::operator<<(std::ostream& (*)(std::ostream&))` in [src/openms/source/FORMAT/SVOutStream.cpp:101](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SVOutStream.cpp#L101).

**Trigger:** Write field, apply std::ends then std::endl then write next field.

**Issue:** Probe string contains NUL plus newline, never equals newline; bookkeeping emits unwanted separator at next line start.

**Proposed C++ fix:** Reset probe buffer/error state per manipulator or identify line-end operations explicitly. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native named end_line/newline avoid manipulator-detection state. See [src/format/sv_out_stream.rs](src/format/sv_out_stream.rs), [tests/sv_out_stream.rs](tests/sv_out_stream.rs) and [SV_OUT_STREAM_SUPPORT.md](docs/SV_OUT_STREAM_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-176 — TrafoXML omits unsupported parameter types after a non-fatal diagnostic

**Status and source:** source-reviewed behavior; unconfirmed defect / API error-policy difference. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `TransformationXMLFile::onStartElement` in [src/openms/source/FORMAT/TransformationXMLFile.cpp:145](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/TransformationXMLFile.cpp#L145); `Internal::XMLHandler::error` in [src/openms/source/FORMAT/HANDLERS/XMLHandler.cpp:71](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/HANDLERS/XMLHandler.cpp#L71).

**Trigger:** Param name slope type bogus value1.

**Issue:** Unsupported type logs a non-fatal error and omits the parameter. This can cause missing/default model parameters, but non-fatal continuation is explicit source behavior; whether it violates the public contract is not established.

**Proposed C++ fix:** Throw ParseError for unsupported types. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. candidate requiring public error-handling contract review.

**Rust handling:** Native Error::Unsupported preserves explicit failure. See [src/format/transformation_xml.rs](src/format/transformation_xml.rs), [tests/transformation_xml.rs](tests/transformation_xml.rs) and [TRANSFORMATION_XML_SUPPORT.md](docs/TRANSFORMATION_XML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-177 — Linear transformation accepts symmetric_regression but never uses it

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `TransformationModelLinear::TransformationModelLinear` in [src/openms/source/ANALYSIS/MAPMATCHING/TransformationModelLinear.cpp:18](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/ANALYSIS/MAPMATCHING/TransformationModelLinear.cpp#L18).

**Trigger:** Fit non-collinear data with symmetric_regression=true versus false.

**Issue:** Constructor stores symmetric_ but fitting path never consults it; the documented regression on y-x versus y+x is not selected.

**Proposed C++ fix:** Implement declared symmetric regression or reject/remove option and correct docs. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Native fitting behavior remains same; setting retained as documented source limitation. See [src/format/transformation_xml.rs](src/format/transformation_xml.rs), [tests/transformation_xml.rs](tests/transformation_xml.rs) and [TRANSFORMATION_XML_SUPPORT.md](docs/TRANSFORMATION_XML_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-178 — MSstats missing design pair silently uses sample0

**Status and source:** source-reviewed unchecked lookup; unconfirmed reachable defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MSstatsFile::storeLFQ` in [src/openms/source/FORMAT/MSstatsFile.cpp:432](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MSstatsFile.cpp#L432).

**Trigger:** A consensus channel label/file pair absent from the design but passing storeLFQ filename-subset and one-label guards; exact input not yet constructed.

**Issue:** If a missing pair reaches operator[], it inserts zero and may read sample row zero; upstream guards exist and a concrete reachable malformed pair remains to be demonstrated.

**Proposed C++ fix:** Check design/run membership and throw MissingInformation. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. candidate requiring a reachable file/label mismatch despite upfront filename and single-label checks.

**Rust handling:** Native lookup checks file/label and run mappings and returns MissingInformation; exact reachable C++ missing-pair case remains pending. See [src/format/msstats.rs](src/format/msstats.rs), [tests/msstats.rs](tests/msstats.rs) and [MSSTATS_SUPPORT.md](docs/MSSTATS_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-179 — MSstats unknown summarization method writes zero intensities

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MSstatsFile::constructFile_` in [src/openms/source/FORMAT/MSstatsFile.cpp:174](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MSstatsFile.cpp#L174).

**Trigger:** Library caller sets method bogus with nonzero intensities.

**Issue:** No branch assigns initial0, but CSV still emitted.

**Proposed C++ fix:** Validate method before processing. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Typed RetentionTimeSummarization rejects unknown names. See [src/format/msstats.rs](src/format/msstats.rs), [tests/msstats.rs](tests/msstats.rs) and [MSSTATS_SUPPORT.md](docs/MSSTATS_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## CPP-180 — MSstats aggregation collapses equal intensities at distinct times

**Status and source:** source-reviewed defect; no C++ execution of this defect. Revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `MSstatsFile::constructFile_` in [src/openms/source/FORMAT/MSstatsFile.cpp:140](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/MSstatsFile.cpp#L140).

**Trigger:** Distinct RTs have intensities10,10,20; sum requested.

**Issue:** std::set yields10,20, so sum30 instead of40 and mean15 instead of13.333.

**Proposed C++ fix:** Keep multiset/vector for sample intensities, separate RT duplicate policy. No upstream patch is claimed.

**Evidence:** Source review, with pinned source hashes and line references in [the FORMAT review manifest](tests/data/format_wave_cpp_review.json). This entry has no executed C++ reproduction or sanitizer evidence. confirmed by source review; not executed.

**Rust handling:** Preserves source intensity deduplication before summarizing, explicitly sorting and deduplicating f32 intensities. See [src/format/msstats.rs](src/format/msstats.rs), [tests/msstats.rs](tests/msstats.rs) and [MSSTATS_SUPPORT.md](docs/MSSTATS_SUPPORT.md). Rust regression execution is recorded separately in the wave validation record.

## FORMAT review corrections and unresolved observations

- **MzTabModification position-byte corruption / mzML numeric IDs / OPXL BetaPepEv numeric positions:** Global numeric overloads invalidate claim only char overload is viable. Isolated executed probe emits decimal. Do not log as C++ defect; fix docs and retain no full-SDK-execution claim.
- **MzTabSpectraRef::setSpecRefFile duplicates setSpecRef:** Naming ambiguity/alias alone is not a defect; no distinct promised behavior established.
- **Accepting too-new TrafoXML version:** Warning plus forward-compatibility policy is not demonstrably wrong without concrete misread future data.
- **SVOutStream absence of std::ostream inheritance in Rust:** Native API design difference, not source defect.
- **Percolator short Sage path subtraction:** Unsigned subtraction verified but std::string substr count can clamp oversized length; do not call it memory unsafety/automatic process abort. Wrong derived path/error deserves further contract review.
- **Unbounded allocation, non-atomic file I/O and permissive count handling:** Need concrete defect/contract assessment; do not log every native resource-policy addition as C++ bug.
- **MzTabM find_if end dereference:** Unchecked dereference visible but exporter constructs both collections internally; malformed externally supplied IDs may not reach this loop. Need reachable graph trigger before confirmed defect.
- **MzTabDouble value-only comparison:** State conflation is real but comparison contract may deliberately be value-based. Require public contract review before confirmed defect.
- **Additional FORMAT source anchors:** This bounded pass prioritized concrete incorrect output, missing bounds and self-roundtrip defects; remaining named doc findings are pending deeper verification, not certified false or complete.
- **All other defect anchors in support documents:** Pending deeper source and contract verification; surveyed does not mean independently verified. Do not promote them automatically to confirmed defects.

The numeric-overload correction has isolated executed evidence under
`../oracle/format-numeric-overload/`, with hashes recorded in the affected
format manifests. It substitutes small formatting adapters and is not a full
SDK build. No new defect ID is assigned to the disproved narrowing claim.

## Additional FORMAT observations awaiting source review

The [retained review queue](tests/data/format_cpp_review_queue.json) gives stable
`FORMAT-OBS-001` through `FORMAT-OBS-134` identifiers to the remaining leaf-document
reports, including the reported affected files, behavior and proposed changes.
These are **unreviewed observations, not 134 additional defects**: the queue
includes policy differences, possible duplicates and incomplete hypotheses.
Known overlaps point back to existing CPP entries. Resolve duplicates and check
source contracts before promoting a report; retain its observation ID when
recording the disposition. Source file hashes establish provenance only.

## SQLite S0 planning candidates

These candidates were identified while reading the next wave's source, before
implementation. They do not count as completed Rust functionality.

## CPP-181 — SqliteConnector permits copying an owned database handle

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Implicit copy constructor/assignment and ~SqliteConnector` in [src/openms/include/OpenMS/FORMAT/SqliteConnector.h:40](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/include/OpenMS/FORMAT/SqliteConnector.h#L40). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Copy an open connector and destroy both copies. The implicit copy duplicates db_ without ownership transfer; both destructors call sqlite3_close_v2 on the same handle. Assignment also loses the previous handle.

**Proposed C++ fix:** Delete copy operations and provide explicit move ownership if needed. No upstream change is claimed.

**Evidence:** Direct source review during S0 planning. No executed C++ reproduction, sanitizer result or Rust test is claimed; error-path and public-contract confirmation remains part of S0.

**Rust handling:** SQLite connector port has not started. Planned implementation uses rusqlite ownership, bound values, quoted identifiers and checked row/column APIs; no implemented correction or native test is claimed.

## CPP-182 — SqliteConnector does not close a handle when opening fails

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `SqliteConnector::openDatabase_` in [src/openms/source/FORMAT/SqliteConnector.cpp:32](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L32). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Construct a connector for a missing read-only database or an inaccessible path. sqlite3_open_v2 may return an allocated error handle. The constructor stores it then throws; the object destructor is not run, leaving that handle unclosed.

**Proposed C++ fix:** Close any returned handle before throwing, or use an owning local guard until successful construction. No upstream change is claimed.

**Evidence:** Direct source review during S0 planning. No executed C++ reproduction, sanitizer result or Rust test is claimed; error-path and public-contract confirmation remains part of S0.

**Rust handling:** SQLite connector port has not started. Planned implementation uses rusqlite ownership, bound values, quoted identifiers and checked row/column APIs; no implemented correction or native test is claimed.

## CPP-183 — SQLite bound statements leak on bind or step errors

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::SqliteHelper::executeBindStatement` in [src/openms/source/FORMAT/SqliteConnector.cpp:156](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L156). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Supply more blobs than placeholders, or violate a UNIQUE constraint during the step. Both error branches throw before sqlite3_finalize; the source TODO acknowledges the leak. The statement may also retain borrowed SQLITE_STATIC buffers.

**Proposed C++ fix:** Use a statement RAII guard on all exits; bind owned data or ensure it outlives the statement. No upstream change is claimed.

**Evidence:** Direct source review during S0 planning. No executed C++ reproduction, sanitizer result or Rust test is claimed; error-path and public-contract confirmation remains part of S0.

**Rust handling:** SQLite connector port has not started. Planned implementation uses rusqlite ownership, bound values, quoted identifiers and checked row/column APIs; no implemented correction or native test is claimed.

## CPP-184 — SQLite table-name helpers interpolate names as SQL syntax

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `SqliteConnector::tableExists, columnExists and countTableRows` in [src/openms/source/FORMAT/SqliteConnector.cpp:118](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L118). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Create a table whose name contains an apostrophe or SQL punctuation, then query it by its actual name. tableExists inserts the name into a quoted literal without escaping; countTableRows and columnExists insert an unquoted identifier. Valid names can fail or alter the query meaning.

**Proposed C++ fix:** Bind values such as sqlite_master.name and correctly quote SQL identifiers with doubled double quotes. No upstream change is claimed.

**Evidence:** Direct source review during S0 planning. No executed C++ reproduction, sanitizer result or Rust test is claimed; error-path and public-contract confirmation remains part of S0.

**Rust handling:** SQLite connector port has not started. Planned implementation uses rusqlite ownership, bound values, quoted identifiers and checked row/column APIs; no implemented correction or native test is claimed.

## CPP-185 — SQLite query helpers ignore step errors and can leak statements

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `SqliteConnector::countTableRows; Internal::SqliteHelper::tableExists and columnExists` in [src/openms/source/FORMAT/SqliteConnector.cpp:77](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L77). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Cause sqlite3_step to return SQLITE_BUSY or another error while querying an existing table. The helpers inspect column types without checking the step result. Existence queries can conflate an execution error with absence; countTableRows throws on NULL before finalizing. Reachable lock/error behavior still needs execution.

**Proposed C++ fix:** Check SQLITE_ROW/DONE explicitly, report other statuses, and finalize through RAII. No upstream change is claimed.

**Evidence:** Direct source review during S0 planning. No executed C++ reproduction, sanitizer result or Rust test is claimed; error-path and public-contract confirmation remains part of S0.

**Rust handling:** SQLite connector port has not started. Planned implementation uses rusqlite ownership, bound values, quoted identifiers and checked row/column APIs; no implemented correction or native test is claimed.

## CPP-186 — SQLite string extraction truncates embedded NUL bytes

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::SqliteHelper::extractValue<std::string> / extractString` in [src/openms/source/FORMAT/SqliteConnector.cpp:218](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L218). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Extract a SQLite TEXT value containing bytes a, NUL, b. The helper constructs std::string from a C string and loses bytes after the first NUL.

**Proposed C++ fix:** Use sqlite3_column_bytes with the returned pointer and retain an explicit byte length. No upstream change is claimed.

**Evidence:** Direct source review during S0 planning. No executed C++ reproduction, sanitizer result or Rust test is claimed; error-path and public-contract confirmation remains part of S0.

**Rust handling:** SQLite connector port has not started. Planned implementation uses rusqlite ownership, bound values, quoted identifiers and checked row/column APIs; no implemented correction or native test is claimed.

## CPP-187 — SQLite integer-to-string extraction narrows to 32 bits

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::SqliteHelper::extractValueIntStr` in [src/openms/source/FORMAT/SqliteConnector.cpp:261](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L261). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Read the SQLite INTEGER value 4294967297 through extractValueIntStr. The helper uses sqlite3_column_int despite SQLite INTEGER being signed 64-bit; conversion to text follows an already narrowed value.

**Proposed C++ fix:** Use sqlite3_column_int64 before decimal formatting. No upstream change is claimed.

**Evidence:** Direct source review during S0 planning. No executed C++ reproduction, sanitizer result or Rust test is claimed; error-path and public-contract confirmation remains part of S0.

**Rust handling:** SQLite connector port has not started. Planned implementation uses rusqlite ownership, bound values, quoted identifiers and checked row/column APIs; no implemented correction or native test is claimed.

## CPP-188 — MSDataWritingConsumer class test is disabled and calls a nonexistent constructor

**Status and source:** Source-reviewed test-coverage issue at `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; no compiler reproduction claimed.

**Affected files and functions:** `src/tests/class_tests/openms/source/MSDataWritingConsumer_test.cpp:25–29`, constructor section; `src/tests/class_tests/openms/executables.cmake:239`, test registration; `src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h:76`, constructor declaration. All three files are hashed in [the consumer manifest](tests/data/ms_data_writing_consumer_provenance.json).

**Trigger and issue:** Re-enable the commented-out class-test registration. The test instantiates `MSDataWritingConsumer()` while the header requires a filename; the stale test cannot exercise the current API. In its present disabled state it supplies no executed regression coverage. This is a test-maintenance defect, not a claim that the production consumer necessarily fails.

**Proposed C++ fix:** Update the fixture to a concrete consumer with a temporary output path, implement the empty test sections against the current contract, and re-enable registration. No upstream patch is claimed.

**Rust handling:** The native consumer suite is enabled and exercises settings, counts, streaming and output checks. Its expectations are source-derived and independent Rust checks; the disabled C++ test is not an executed oracle.
