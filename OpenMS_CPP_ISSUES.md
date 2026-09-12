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
