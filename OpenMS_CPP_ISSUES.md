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
| CPP-027 | mzML writing discards processing completion seconds | Source-reviewed; executed (CLI-2 oracle) | Open |
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
Executed since (early-TOPP wave 2, CLI-2 `f886d90`): the product SDK at core
`4fdec46` writes the `-test` completion time `1999-12-31 23:59:59` of
SpectraFilterWindowMower's processing record as `1999-12-31+23:59`
(`../oracle/topp-cli-lifecycle/cli2/manifest.json`, the retained
`SpectraFilterWindowMower_1_output.mzML`), and `MzMLHandler.cpp:3947` is
unchanged at `bc9cc12`. The CLI-2 processing-record comparisons truncate the
port's seconds to match.

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

## CPP-172 — mzML writers emit `dataProcessingRef` and `sourceFileRef` the header never declares

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`, on `ibminode06`.

**Status:** Executed. (Was: source-reviewed, no C++ execution, and scoped to the streaming consumer. Both the reproduction and the wider scope below are new.)

**Affected files/functions:** `src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp:5060-5069`, `MzMLHandler::writeHeader_`, which builds `dps` by **content**; `:5252-5255` and `:5258-5272`, `MzMLHandler::writeSpectrum_`, which searches the same `dps` by **pointer**; `src/openms/source/FORMAT/DATAACCESS/MSDataWritingConsumer.cpp:53-95`, `MSDataWritingConsumer::consumeSpectrum`, which additionally calls `writeHeader_` on a one-record dummy map so that its `dps_` never grows past one entry.

**Trigger:** Two independent ones.

1. **The ordinary whole-document writer,** `MzMLFile::store`, on any experiment whose records carry `DataProcessing` vectors that are equal in content but held in distinct objects. **Every `FileMerger` output is such a file** — `FileMerger` writes with `FileHandler().storeExperiment` (`FileMerger.cpp:538` at the TOPP pin `174b576`), and each merged part contributes its own history object.
2. **The streaming consumer,** `MSDataWritingConsumer`, on any mzML whose records do not all share the first record's `sourceFile` and `dataProcessing`. Reached from every TOPP tool that uses it — `PeakPickerHiRes`, `PeakPickerIM`, `NoiseFilterGaussian`, `NoiseFilterSGolay` and `FileConverter` in their low-memory modes, plus `CometAdapter`, `OpenSwathMzMLFileCacher`, `OpenSwathWorkflow`, `SageAdapter` and `TICCalculator`.

**Issue:** The two halves of the writer disagree about what makes two processing histories the same.

`writeHeader_` deduplicates them by content: `already_present = OpenMS::Helpers::cmpPtrContainer(exp[s].getDataProcessing(), dps[j])` (`:5060-5069`), and `cmpPtrContainer` reduces to `cmpPtrSafe`, whose own comment reads "We are not interested whether the pointers are equal but whether the contents are equal" and whose body is `*a == *b` (`Helpers.h:35-51`). Content-equal histories therefore collapse into one declared `dataProcessing` entry.

`writeSpectrum_` then searches that same `dps` by pointer: `spec.getDataProcessing() != dps[0]` (`:5258`) and `spec.getDataProcessing() == dps[i]` (`:5265`) over `std::vector<std::shared_ptr<const DataProcessing>>` (`SpectrumSettings.h:165`), which compares the `shared_ptr`s, not the objects. A record whose history is content-equal to a declared entry but a distinct object matches nothing, and the search falls through to `dp_ref_num = s` — the record's own position in the stream — so the attribute is written as `dataProcessingRef="dp_sp_<s>"` naming an element that does not exist.

`sourceFileRef` has a second, simpler form of the same fault: `:5252-5255` writes `sourceFileRef="sf_sp_<s>"` whenever the record has a non-default source file, unconditionally, without consulting the header at all, so it dangles for every record after the first even when the source file *is* the first record's. The same numbering applies to a reference a binary data array carries, written as `dp_sp_<s>_bi_<m>` (`:5567`, `:5597`, `:5806` for a spectrum's three array kinds, `:5965` and `:5997` for a chromatogram's).

In the streaming consumer both faults are amplified but not caused: only the first record reaches `writeHeader_`, so `dps_` holds exactly one entry and never grows, and the source's own `// TODO writeSpectrum assumes that dps_ has at least one value -> assert this here` (`MSDataWritingConsumer.cpp:93`) marks the gap.

The result is not valid mzML. Both attributes are `xs:IDREF` on `SpectrumType` (`share/OpenMS/SCHEMAS/mzML_1_10.xsd:851`, `:856`) against `xs:ID` on `DataProcessingType` and `SourceFileType`, and `dataProcessingRef` additionally carries `xs:keyref KEYREF_DPREF`, whose `refer` is `KEY_DP_ID`, the `id` of a `dataProcessingList/dataProcessing` (`:1064-1071`, `:983-990`).

A chromatogram is the silent form of the same defect: `writeChromatogram_` (`MzMLHandler.cpp:5879`) writes `id`, `index` and `defaultArrayLength` and no reference at all, so a chromatogram whose history is not the first record's is written under the list's `defaultDataProcessingRef` and its own history is lost without a diagnostic.

**Proposed C++ fix:** Make the two halves agree. In `writeSpectrum_`, search `dps` with `Helpers::cmpPtrContainer` — the same comparison `writeHeader_` used to build it — instead of `operator==`; that alone removes the whole-document trigger, because every history the experiment holds is then found. `sourceFileRef` should likewise be resolved against the header and not renumbered when the record's source file is one the header declares. For the streaming consumer a second remedy is still needed, because there the wanted entry genuinely is not in the published header: either publish the whole `dataProcessingList` up front from `setExperimentalSettings`, or omit the attribute and let the record inherit `defaultDataProcessingRef`, which loses information but keeps the file valid, or throw. Whichever is chosen, `MSDataWritingConsumer`'s class documentation should say so and the `// TODO` at `.cpp:93` should be resolved. Add a class-test case with two records carrying different histories (see CPP-188: that test is currently disabled).

**Evidence:** Executed on the Release build on `ibminode06`. Drivers and logs under `../oracle/p4-lowmemory` (`closediff1_06.sh`, `closediff2_06.sh`, `closediff3_06.sh`; `logs/closediff1_06.log`, `logs/closediff2_06.log`, `logs/closediff3_06.log`) and, for the whole-document trigger and the content/pointer split, `../oracle/integ-w7` (`refcheck_06.sh`, `dupdp_06.sh`; logs under the wave-7 integration log directory).

* **Whole-document writer.** The ~9.3 MB, 110-spectrum file the Release `FileMerger` builds from 22 copies of `PeakPickerHiRes_input.mzML` declares `<sourceFileList count="22">` and **exactly one** `<dataProcessing id="dp_sp_0">`, and carries **106** record `dataProcessingRef`s of which **105 dangle**, `dp_sp_5` through `dp_sp_109`; records 1 to 4 carry none, because their history is the first record's, and record 0 carries the one that resolves. The single declared entry against 22 merged parts is the content deduplication; the 105 references are the pointer comparison. A two-part merge is the minimal case: one declared `dp_sp_0`, six references, five dangling (`dp_sp_5` … `dp_sp_9`). No streaming consumer is involved in either.
* **Content against pointer, isolated.** On a five-record fixture whose `dp_sp_0` and `dp_sp_1` are byte-identical apart from their `id` (the committed `PeakPickerHiRes_refs_input.mzML` with `dp_sp_1`'s `softwareRef` repointed from `so_dp_1` to `so_dp_0`), `PeakPickerHiRes -processOption inmemory` declares one `dataProcessing`, `dp_sp_0`, and still writes `dataProcessingRef="dp_sp_1"` and `"dp_sp_2"` on records 1 and 2 — two dangling references in a 12 KB file. The control is the same tool and mode on the unmodified fixture, where the histories differ in content: it declares `dp_sp_0` and `dp_sp_1` and dangles nothing.
* **Streaming consumer.** On the same 110-record file, `PeakPickerHiRes -test -no_progress -processOption lowmemory` exits 0 and writes all 110 records; its output declares `sf_ru_0` and `dp_sp_0` and reproduces the same 106 references and the same 105 dangling identifiers. On the five-record `refs` fixture the whole numbering rule is visible: records 0 to 4 come out as `sourceFileRef="sf_sp_0"` … `"sf_sp_4"` against a header declaring one record source file, and `dataProcessingRef` appears as `dp_sp_0`, `dp_sp_1`, `dp_sp_2` and then not at all. Records 1 and 2 share one `dataProcessing` in the input and still get two different dangling identifiers, and record 4's source file *is* the first record's and it still gets the dangling `sf_sp_4`.
* **The source's own reader accepts these files:** `FileInfo` exits 0 and prints `Error: unregistered source file reference sf_sp_1.` for the spectrum references (`MzMLHandler.cpp:899-906`) and says nothing at all about the dangling `dataProcessingRef`, which CPP-262 already covers: `processing_[ref]` resolves with `std::map::operator[]`, so an unknown id silently yields an empty history.

**Rust handling:** `ReferencePolicy`, in `src/format/ms_data_writing_consumer.rs`. `Checked`, the default, refuses the record with `Error::Unsupported` rather than writing a reference that cannot resolve. `SourceDangling` reproduces the source in the source's own `sf_sp_<s>`, `dp_sp_<s>` and `dp_sp_<s>_bi_<m>` spelling — this writer's own identifiers are zero-padded positions in one namespace, so reusing them would make the reference *resolve*, to the wrong entry, instead of dangling — and `PeakPickerHiRes -processOption lowmemory` selects it, because refusing stops that mode after five records on any `FileMerger` output. On the 110-record file the port's low-memory output carries the same 105 dangling identifiers as the C++ one and the same decoded content. It is **not** an exact reproduction, and is not claimed as one: the port has no pointer identity in its model and decides "differs from the first record's" by the text the history renders, so on the `dp_sp_0`/`dp_sp_1` fixture above it writes no `dataProcessingRef` where the source dangles one. The two agree on every input whose textually distinct histories are also distinct objects. See [src/format/ms_data_writing_consumer.rs](src/format/ms_data_writing_consumer.rs), [tests/ms_data_writing_consumer.rs](tests/ms_data_writing_consumer.rs), [MS_DATA_WRITING_CONSUMER_SUPPORT.md](docs/MS_DATA_WRITING_CONSUMER_SUPPORT.md) and native difference 12 of [TOPP_PEAK_PICKER_HI_RES_SUPPORT.md](docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md).

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

**Evidence:** Direct pinned-source review. This issue has no executed C++ reproduction or sanitizer result; native S0 handling and regression checks are distinct evidence, described below and in the connector provenance.

**Rust handling:** SqliteConnector owns a rusqlite Connection and implements neither Clone nor Copy. Dropping an uncommitted connection rolls back through SQLite ownership; no duplicate raw handle API is exposed.

## CPP-182 — SqliteConnector does not close a handle when opening fails

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `SqliteConnector::openDatabase_` in [src/openms/source/FORMAT/SqliteConnector.cpp:32](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L32). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Construct a connector for a missing read-only database or an inaccessible path. sqlite3_open_v2 may return an allocated error handle. The constructor stores it then throws; the object destructor is not run, leaving that handle unclosed.

**Proposed C++ fix:** Close any returned handle before throwing, or use an owning local guard until successful construction. No upstream change is claimed.

**Evidence:** Direct pinned-source review. This issue has no executed C++ reproduction or sanitizer result; native S0 handling and regression checks are distinct evidence, described below and in the connector provenance.

**Rust handling:** Connection::open_with_flags owns failed-open cleanup. Missing read-only/read-write files are rejected and remain absent. Native tests check failures and later successful opening; they do not measure allocator leaks.

## CPP-183 — SQLite bound statements leak on bind or step errors

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::SqliteHelper::executeBindStatement` in [src/openms/source/FORMAT/SqliteConnector.cpp:156](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L156). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Supply more blobs than placeholders, or violate a UNIQUE constraint during the step. Both error branches throw before sqlite3_finalize; the source TODO acknowledges the leak. The statement may also retain borrowed SQLITE_STATIC buffers.

**Proposed C++ fix:** Use a statement RAII guard on all exits; bind owned data or ensure it outlives the statement. No upstream change is claimed.

**Evidence:** Direct pinned-source review. This issue has no executed C++ reproduction or sanitizer result; native S0 handling and regression checks are distinct evidence, described below and in the connector provenance.

**Rust handling:** Prepared statements are owned RAII values. Excess bindings, constraint errors and returned rows report errors and release the statement. Native regression tests verify that subsequent queries and bindings remain usable.

## CPP-184 — SQLite table-name helpers interpolate names as SQL syntax

**Status and source:** Reproduced in an adapted C++ probe on kim at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; not a full SDK build.

**Affected files and functions:** `SqliteConnector::tableExists, columnExists and countTableRows` in [src/openms/source/FORMAT/SqliteConnector.cpp:118](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L118). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Create a table whose name contains an apostrophe or SQL punctuation, then query it by its actual name. tableExists inserts the name into a quoted literal without escaping; countTableRows and columnExists insert an unquoted identifier. Valid names can fail or alter the query meaning.

**Proposed C++ fix:** Bind values such as sqlite_master.name and correctly quote SQL identifiers with doubled double quotes. No upstream change is claimed.

**Evidence:** Exact pinned connector implementation and declarations were compiled with small substitute StandardTypes, Exception and StringUtils support headers against host SQLite 3.45.1 (GCC 13.3, C++20). The ordinary table T exists; ABSENT does not. The existing table named odd'name causes IllegalArgument during preparation. The absent literal name "ABSENT' OR 1=1 --" incorrectly returns true. The existing table named "odd table" causes preparation errors in both row-count and column-existence queries. Probe sources, adapters, executable and logs are hashed in [the connector provenance](tests/data/sqlite_connector_provenance.json) and retained externally under `../oracle/sqlite-connector-s0-probe/`. This does not execute the original OpenMS exception ABI or full SDK dependency closure.

**Rust handling:** table_exists binds its name as a value; count_table_rows and column_exists quote the entire identifier and double embedded quotes. Native tests cover punctuation-bearing names, injection-shaped input and NUL identifier rejection.

## CPP-185 — SQLite query helpers ignore step errors and can leak statements

**Status and source:** Reproduced in an adapted C++ probe on kim at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; not a full SDK build.

**Affected files and functions:** `SqliteConnector::countTableRows; Internal::SqliteHelper::tableExists and columnExists` in [src/openms/source/FORMAT/SqliteConnector.cpp:77](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L77). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Cause sqlite3_step to return SQLITE_BUSY or another error while querying an existing table. The helpers inspect column types without checking the step result. Existence queries can conflate an execution error with absence; countTableRows throws on NULL before finalizing. The adapted probe reproduced this with an exclusive lock after a warm read.

**Proposed C++ fix:** Check SQLITE_ROW/DONE explicitly, report other statuses, and finalize through RAII. No upstream change is claimed.

**Evidence:** Exact pinned connector implementation and declarations were compiled with small substitute StandardTypes, Exception and StringUtils support headers against host SQLite 3.45.1 (GCC 13.3, C++20). After a warm read returns 3 rows, a second connection holds BEGIN EXCLUSIVE. A direct trace records successful preparation (SQLITE_OK, 0) followed by SQLITE_BUSY (5) at step. Both tableExists(T) and columnExists(T, ID) return false, with database error code 5. countTableRows(T) throws the adapted SqlOperationFailed and leaves one outstanding statement, counted using sqlite3_next_stmt. After the second connection rolls back, the count is again 3. Probe sources, adapters, executable and logs are hashed in [the connector provenance](tests/data/sqlite_connector_provenance.json) and retained externally under `../oracle/sqlite-connector-s0-probe/`. This does not execute the original OpenMS exception ABI or full SDK dependency closure.

**Rust handling:** Every query checks rusqlite preparation and row-step results; errors propagate as Error::Io rather than false. Owned statements clean up on all exits. Native locked-query tests cover all three operations and successful reuse after unlocking.

## CPP-186 — SQLite string extraction truncates embedded NUL bytes

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::SqliteHelper::extractValue<std::string> / extractString` in [src/openms/source/FORMAT/SqliteConnector.cpp:218](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L218). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Extract a SQLite TEXT value containing bytes a, NUL, b. The helper constructs std::string from a C string and loses bytes after the first NUL.

**Proposed C++ fix:** Use sqlite3_column_bytes with the returned pointer and retain an explicit byte length. No upstream change is claimed.

**Evidence:** Direct pinned-source review. This issue has no executed C++ reproduction or sanitizer result; native S0 handling and regression checks are distinct evidence, described below and in the connector provenance.

**Rust handling:** The S0 public connector is implemented, but the private extraction helper remains unported. Binary-binding tests preserve NUL bytes; they are not evidence for a text-extraction replacement.

## CPP-187 — SQLite integer-to-string extraction narrows to 32 bits

**Status and source:** Source-reviewed candidate; not reproduced in running C++, revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `Internal::SqliteHelper::extractValueIntStr` in [src/openms/source/FORMAT/SqliteConnector.cpp:261](https://github.com/okohlbacher/OpenMS4-core/blob/bc9cc12514c768385ce121d6ca4bb710fe1983c4/src/openms/source/FORMAT/SqliteConnector.cpp#L261). Public/private declarations are pinned in [the SQLite review manifest](tests/data/sqlite_connector_review.json).

**Trigger and issue:** Read the SQLite INTEGER value 4294967297 through extractValueIntStr. The helper uses sqlite3_column_int despite SQLite INTEGER being signed 64-bit; conversion to text follows an already narrowed value.

**Proposed C++ fix:** Use sqlite3_column_int64 before decimal formatting. No upstream change is claimed.

**Evidence:** Direct pinned-source review. This issue has no executed C++ reproduction or sanitizer result; native S0 handling and regression checks are distinct evidence, described below and in the connector provenance.

**Rust handling:** The S0 public connector is implemented, but the private integer-to-string extraction helper remains unported. Row counts use a checked i64-to-usize conversion; no extraction-helper correction is claimed.

## CPP-188 — MSDataWritingConsumer class test is disabled and calls a nonexistent constructor

**Status and source:** Source-reviewed test-coverage issue at `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; no compiler reproduction claimed.

**Affected files and functions:** `src/tests/class_tests/openms/source/MSDataWritingConsumer_test.cpp:25–29`, constructor section; `src/tests/class_tests/openms/executables.cmake:239`, test registration; `src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h:76`, constructor declaration. All three files are hashed in [the consumer manifest](tests/data/ms_data_writing_consumer_provenance.json).

**Trigger and issue:** Re-enable the commented-out class-test registration. The test instantiates `MSDataWritingConsumer()` while the header requires a filename; the stale test cannot exercise the current API. In its present disabled state it supplies no executed regression coverage. This is a test-maintenance defect, not a claim that the production consumer necessarily fails.

**Proposed C++ fix:** Update the fixture to a concrete consumer with a temporary output path, implement the empty test sections against the current contract, and re-enable registration. No upstream patch is claimed.

**Rust handling:** The native consumer suite is enabled and exercises settings, counts, streaming and output checks. Its expectations are source-derived and independent Rust checks; the disabled C++ test is not an executed oracle.

## CPP-189 — SqliteConnector row-count documentation names the wrong exception

**Status and source:** Reproduced in an adapted C++ probe on kim at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; not a full SDK build.

**Affected files and functions:** `SqliteConnector::countTableRows`, `src/openms/include/OpenMS/FORMAT/SqliteConnector.h:75`; `src/openms/source/FORMAT/SqliteConnector.cpp:82,145–153`; `src/tests/class_tests/openms/source/SqliteConnector_test.cpp:105–106`. Files are hashed in [the connector provenance](tests/data/sqlite_connector_provenance.json).

**Trigger and issue:** Call `countTableRows("UNKNOWN")` on a database without that table. The installed header promises `SqlOperationFailed`, but preparing the SELECT throws `IllegalArgument`; the class test explicitly expects that latter type. A caller handling only the documented exception can miss the actual failure.

**Proposed C++ fix:** Correct the header to document `IllegalArgument` for an unknown table. Alternatively, deliberately translate the exception and update the class test as an API change. No upstream fix is claimed.

**Evidence:** Exact pinned connector implementation and declarations were compiled with small substitute StandardTypes, Exception and StringUtils support headers against host SQLite 3.45.1 (GCC 13.3, C++20). countTableRows(UNKNOWN) throws the adapted IllegalArgument from the exact pinned prepareStatement path, agreeing with the source class test and contradicting the installed header. Substitute exception classes preserve distinct throw/catch paths; the original exception inheritance, ABI and metadata formatting are not exercised. Probe sources, adapters, executable and logs are hashed in [the connector provenance](tests/data/sqlite_connector_provenance.json) and retained externally under `../oracle/sqlite-connector-s0-probe/`. This does not execute the original OpenMS exception ABI or full SDK dependency closure.

**Rust handling:** The native rustdoc and support document state the source mismatch and specify Error::Io with an underlying SQLite cause for an absent table. The source class-test expectation is mapped to a native error regression.
## CPP-190 — Default handler writing reads an uninitialized SQL batch size

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:269` (`MzMLSqliteHandler::MzMLSqliteHandler`); `src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h:121` (`MzMLSqliteHandler::setConfig`); `src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h:238` (`MzMLSqliteHandler::sql_batch_size_`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:1261` (`MzMLSqliteHandler::writeSpectra`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:1472` (`MzMLSqliteHandler::writeChromatograms`).

**Trigger:** Construct a handler, createTables(), then write a nonempty spectra/chromatogram vector without first calling setConfig().

**Issue:** sql_batch_size_ is an int with no in-class initializer and is absent from the constructor initializer list. setConfig is its only assignment; writes compare sql_it to the indeterminate member. Higher-level SqMassFile calls setConfig, but the direct public API and source class-test writeExperiment path do not require it.

**Proposed C++ fix:** Initialize sql_batch_size_ to 500 in-class or in the constructor, keeping setConfig override; validate positive overrides. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. Initialize a typed default configuration explicitly and validate overrides.

## CPP-191 — Array hydration lacks pair length and role validation

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:196` (`populateContainer_sub_`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:255` (`populateContainer_sub_`).

**Trigger:** A spectrum DATA coordinate array has 2 doubles and its intensity array has 1 (or 3); separately, provide two coordinate-role rows and no intensity-role row.

**Issue:** Only the first nonempty array determines container length. Every subsequent array is copied using container length and an unchecked decoded-data iterator: shorter data is read beyond end, longer data is truncated. cont_data only counts rows and accepts >=2, so duplicate roles can satisfy completeness while a required role is absent.

**Proposed C++ fix:** Decode into separate role slots; reject duplicate or unexpected roles, require exactly coordinate+intensity, verify equal lengths (including empty arrays), then construct peaks with bounded indexing. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. Validate roles and equal lengths before constructing any returned spectrum/chromatogram.

## CPP-192 — Blob hydration assigns objects by SQL row encounter order instead of record identity

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:122` (`populateContainer_sub_`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:559` (`MzMLSqliteHandler::populateSpectraWithData_`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:511` (`MzMLSqliteHandler::populateChromatogramsWithData_`).

**Trigger:** Store metadata records ID0/native s0 and ID1/native s1, but insert DATA rows for ID1 before ID0; use a legal query plan that returns DATA insertion order.

**Issue:** The first encountered DATA ID is mapped to containers[0] regardless of that container’s actual SQL ID. Metadata and blob queries have no matching ORDER BY. Differing natural order throws a false native-ID mismatch; identical native IDs can instead put data on the wrong metadata record.

**Proposed C++ fix:** Carry the actual record-ID-to-result-position map out of metadata preparation and use it for every blob row; ORDER BY alone is insufficient if independent metadata snapshots are used. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. The retained Python SQLite observations support the SQL mechanism, not an execution of the C++ control flow. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. Hydrate by explicit record identity; reject missing/extra references and snapshot disagreement.

## CPP-193 — Spectrum/chromatogram IDs and peptide sequence remain unescaped SQL values

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:1167` (`MzMLSqliteHandler::writeSpectra`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:1194` (`MzMLSqliteHandler::writeSpectra`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:1405` (`MzMLSqliteHandler::writeChromatograms`).

**Trigger:** Write a native ID such as scan'1 or a precursor peptide_sequence containing an apostrophe; crafted text may add SQL syntax to the concatenated batch.

**Issue:** Values are surrounded by apostrophes but never escaped or parameterized, making valid text fail SQL parsing and allowing supplied content to alter the statement. The pinned run-path binding fix at lines 897-900 does not cover these fields. This is distinct from CPP184 identifier-helper interpolation.

**Proposed C++ fix:** Bind every text, numeric and binary field through prepared statements; never interpolate record content into SQL. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. Use typed rusqlite parameters throughout; test apostrophes and embedded NULs.

## CPP-194 — Metadata failures leave committed DATA rows and advanced writer counters

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:1261` (`MzMLSqliteHandler::writeSpectra`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:1472` (`MzMLSqliteHandler::writeChromatograms`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:897` (`MzMLSqliteHandler::writeRunLevelInformation`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:873` (`MzMLSqliteHandler::writeExperiment`).

**Trigger:** Cause a metadata insert to fail after blob insertion (e.g. native-ID apostrophe, duplicate ID from reopening an existing database, or SQLite trigger abort).

**Issue:** DATA executeBindStatement calls run before BEGIN TRANSACTION, so their successful inserts autocommit; counters advance before the metadata transaction succeeds. Later metadata failure rolls back at most that metadata transaction and leaves orphan DATA. RUN and RUN_EXTRA are likewise separate commits, and writeExperiment sequences multiple independently committed operations.

**Proposed C++ fix:** Use one transaction for each public write operation, including all DATA/metadata and run snapshot rows; publish counter increments only after commit. Share internal transaction-taking helpers for writeExperiment to avoid nested independent commits. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. Perform prepared bounded writes inside one transaction and update counters after commit; inject later-statement failures in tests.

## CPP-195 — Product and precursor indexes are accidentally built on DATA

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:1047` (`MzMLSqliteHandler::createIndices_`).

**Trigger:** Create a sqMass database and inspect indexes, or query joins against a large PRODUCT/PRECURSOR table.

**Issue:** Four named product/precursor indexes target DATA instead of PRODUCT/PRECURSOR. This adds redundant DATA insertion/index storage overhead while leaving the intended join columns unindexed. The retained original sqMass fixture has the same four wrong targets.

**Proposed C++ fix:** Create product_chr_idx/product_sp_idx on PRODUCT and precursor_chr_idx/precursor_sp_idx on PRECURSOR. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. The retained Python SQLite observations support the SQL mechanism, not an execution of the C++ control flow. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. Create correct table indexes and assert sqlite_master targets; do not require performance differences in correctness fixtures.

## CPP-196 — SWATH selection stops at a matching chromatogram precursor NULL ID

**Status and source:** Reproduced with an adapted C++ probe; source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteSwathHandler.cpp:92` (`MzMLSqliteSwathHandler::readSpectraForWindow`).

**Trigger:** PRECURSOR contains a chromatogram row with NULL SPECTRUM_ID and isolation target 412.5 before matching spectrum rows, then readSpectraForWindow(center=412.5).

**Issue:** The query selects all matching precursor rows, including chromatograms. The loop uses column-0 NULL as end-of-results and therefore stops on the first chromatogram row, returning no or only a prefix of matching spectra. Source writeExperiment writes chromatograms before spectra, so this row ordering can arise naturally.

**Proposed C++ fix:** Filter SPECTRUM_ID IS NOT NULL (prefer join to SPECTRUM for valid identities) and iterate with checked SQLITE_ROW/SQLITE_DONE status instead of a column-value sentinel. No upstream fix is claimed.

**Evidence:** Executed on kim with GCC 13.3/C++20 and host SQLite 3.45.1, using the exact pinned handler, connector and SwathMap sources plus declared type, string-formatting, exception and pointer-only adapters. With center 500, a matching chromatogram row first returned no IDs; placing it between two spectrum precursors returned only ID 1. Removing that row returned IDs 1 and 2. This is adapted execution, not a full SDK or original class-test run. Hashed sources, adapters, executable, database and logs are listed in [the SWATH provenance manifest](tests/data/mzml_sqlite_swath_provenance.json); `../oracle/sqlite-swath-s1-probe/result.log` has SHA-256 `792d2506c10b59e5900518c1afe191e37eed72402a715c7ce03838268f0a4bcb`. The earlier Python SQL observations remain independent SQL evidence.

**Rust handling:** S1 uses fallible typed row iteration and skips NULL spectrum IDs; the remote-tested native regression retains all matching spectrum rows around chromatogram rows.

## CPP-197 — SWATH window docs promise distinct centers but query deduplicates full bounds tuples

**Status and source:** Reproduced with an adapted C++ probe; source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteSwathHandler.h:57` (`MzMLSqliteSwathHandler::readSwathWindows documentation`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteSwathHandler.cpp:30` (`MzMLSqliteSwathHandler::readSwathWindows`).

**Trigger:** Two MS2 precursors have center 412.5 but lower/upper offsets 12.5 and 10.

**Issue:** SELECT DISTINCT(ISOLATION_TARGET), lower, upper applies DISTINCT to all three output columns, returning two windows with the same center. The new header explicitly promises one per distinct center. The tuple result is meaningful for variable widths, so this is a documentation/contract mismatch, not proof that tuple behavior itself should change.

**Proposed C++ fix:** Document distinct (center,lower,upper) tuples and add a differing-width example; if unique centers are intended instead, define how conflicting bounds are handled before changing SQL. No upstream fix is claimed.

**Evidence:** Executed on kim with GCC 13.3/C++20 and host SQLite 3.45.1, using the exact pinned handler, connector and SwathMap sources plus declared type, string-formatting, exception and pointer-only adapters. Two spectrum precursors at center 500 with offsets 10 and 20 returned two windows, 490–510 and 480–520, demonstrating full-tuple distinctness. This is adapted execution, not a full SDK or original class-test run. Hashed sources, adapters, executable, database and logs are listed in [the SWATH provenance manifest](tests/data/mzml_sqlite_swath_provenance.json); `../oracle/sqlite-swath-s1-probe/result.log` has SHA-256 `792d2506c10b59e5900518c1afe191e37eed72402a715c7ce03838268f0a4bcb`. The earlier Python SQL observations remain independent SQL evidence.

**Rust handling:** S1 preserves distinct tuples, including equal centers with differing bounds, and describes that behavior explicitly. Native regression tests compare both windows.

## CPP-198 — Recreating tables on a used handler retains counters from the deleted database

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:939` (`MzMLSqliteHandler::createTables`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:269` (`MzMLSqliteHandler::MzMLSqliteHandler`); `src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h:223` (`writer counter contract`).

**Trigger:** Construct/configure handler, createTables, write two spectra, createTables again on the same object, write a new spectrum, then select database ID0.

**Issue:** createTables deletes the file and recreates an empty database but never resets spec_id_/chrom_id_. The new spectrum is assigned ID2 while count is1; ID0 selection fails. The existing class test recreates through the same handler but checks counts rather than IDs. The counters are documented as global to a particular database file, and downstream SqMassFile transforms assume IDs begin at zero.

**Proposed C++ fix:** Reset both counters only after successful destructive schema recreation; retain existing counters if recreation fails. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. Reset counters on successful create_tables and regression-test IDs across repeated creation.

## CPP-199 — Metadata readers accept invalid negative activation enum values below -1

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:708` (`MzMLSqliteHandler::prepareChroms_`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:852` (`MzMLSqliteHandler::prepareSpectra_`).

**Trigger:** Read a PRECURSOR row with ACTIVATION_METHOD=-2 (or a more negative value) and valid surrounding metadata.

**Issue:** The reader excludes only -1 and checks value < SIZE_OF_ACTIVATIONMETHOD, allowing all smaller negative values through the cast and insertion into the activation-method set. This is outside the named enum domain. No downstream crash is claimed without execution.

**Proposed C++ fix:** Check 0 <= value && value < SIZE_OF_ACTIVATIONMETHOD; handle -1 as unset and reject other invalid values explicitly. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. Use checked enum conversion with a deliberate NULL/-1 unset policy.

## CPP-200 — SqMassFile transform issues an extra empty selected-read batch

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/SqMassFile.cpp:49` (`SqMassFile::transform`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:390` (`MzMLSqliteHandler::readSpectra`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:413` (`MzMLSqliteHandler::readChromatograms`).

**Trigger:** Transform a file containing zero spectra/chromatograms, or exactly 500,1000,... of either record type.

**Issue:** The <= count/batch_size loop runs a final iteration with start=end and an empty indices vector. The selected-read APIs require nonempty input. In a precondition-enabled build this fails immediately; for nonempty exact multiples with preconditions disabled, empty indices mean unrestricted metadata and the subsequent size comparison throws. The normal empty side of a spectra-only/chromatogram-only file also violates the stated precondition.

**Proposed C++ fix:** Loop while idx_start < count or iterate nonempty chunks of actual record IDs; cache counts and skip zero-record types. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. S1 remains unported; address in subsequent SqMassFile consumer wave and add 0/1/499/500/501 boundary tests.

## CPP-201 — Full metadata recovery promise omits auxiliary-array loss

**Status and source:** Unconfirmed documentation/contract candidate at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h:116` (`MzMLSqliteHandler::setConfig documentation`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:904` (`MzMLSqliteHandler::writeRunLevelInformation`); `src/openms/source/KERNEL/MSSpectrum.cpp:203` (`MSSpectrum::clear`); `src/openms/source/KERNEL/MSChromatogram.cpp:335` (`MSChromatogram::clear`).

**Trigger:** Store an experiment with a nonempty auxiliary float/string/integer data array using write_full_meta=true, then load it.

**Issue:** The RUN_EXTRA snapshot invokes clear(false), which clears all auxiliary arrays as well as peaks. DATA writing supports only primary coordinate/intensity arrays, so auxiliary values have no storage path. This contradicts the header parenthetical allowing complete recovery of the input file. Array loss is source-confirmed; whether the format intentionally excludes these arrays is a contract question, so classified as a documentation/contract candidate rather than an unqualified algorithm defect.

**Proposed C++ fix:** At minimum qualify the recovery promise and detect/document unsupported arrays. For complete recovery, provide an explicit format-compatible auxiliary-data representation or reject such input instead of silently accepting it. No upstream fix is claimed.

**Evidence:** Pinned source review. No executed C++ reproduction or sanitizer run is claimed. Source hashes and supporting artifacts are retained in [the S1 review manifest](tests/data/sqlite_s1_review.json).

**Rust handling:** Not implemented in S0. Choose and document an explicit auxiliary-array policy before advertising full experiment recovery. No handling implemented yet.

## CPP-202 — Read accessors create an empty database when the input path is missing

**Status and source:** Reproduced with an adapted C++ probe; source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteSwathHandler.cpp:23` (`MzMLSqliteSwathHandler::readSwathWindows`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteSwathHandler.cpp:58` (`MzMLSqliteSwathHandler::readMS1Spectra`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteSwathHandler.cpp:84` (`MzMLSqliteSwathHandler::readSpectraForWindow`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:280` (`MzMLSqliteHandler::readExperiment`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:390` (`MzMLSqliteHandler::readSpectra`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp:413` (`MzMLSqliteHandler::readChromatograms`); `src/openms/include/OpenMS/FORMAT/SqliteConnector.h:58` (`SqliteConnector default open mode`); `src/openms/source/FORMAT/SqliteConnector.cpp:44` (`SqliteConnector::openDatabase_`).

**Trigger:** Call a read accessor with a nonexistent filename whose parent directory is writable.

**Issue:** Read methods construct SqliteConnector(filename_) without an explicit mode, selecting READWRITE_OR_CREATE. SQLite creates an empty file, then the read query fails because the required tables do not exist. Thus a read error mutates the filesystem and is reported as a missing-table query failure instead of a missing-input open failure. The same call sites unnecessarily request write access to existing databases.

**Proposed C++ fix:** Pass SqlOpenMode::READ_ONLY from all read/count/lookup accessors. Retain creating modes only for explicit create/write operations; test that failed reads leave no file. No upstream fix is claimed.

**Evidence:** Executed on kim with GCC 13.3/C++20 and host SQLite 3.45.1, using the exact pinned handler, connector and SwathMap sources plus declared type, string-formatting, exception and pointer-only adapters. The missing path remained absent after construction. Calling `readMS1Spectra` raised adapted `IllegalArgument` and left a zero-byte file. Other listed read paths remain source-reviewed only. This is adapted execution, not a full SDK or original class-test run. Hashed sources, adapters, executable, database and logs are listed in [the SWATH provenance manifest](tests/data/mzml_sqlite_swath_provenance.json); `../oracle/sqlite-swath-s1-probe/result.log` has SHA-256 `792d2506c10b59e5900518c1afe191e37eed72402a715c7ce03838268f0a4bcb`. The earlier Python SQL observations remain independent SQL evidence.

**Rust handling:** S1 SWATH accessors open with explicit ReadOnly mode; native missing-path tests verify noncreation. The main storage handler is being validated separately and is not covered by the adapted SWATH probe.

## CPP-203 — SWATH accessors can return partial success after SQLite step errors

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; the SQLite error trigger is executed in a Rust regression, not a C++ reproduction.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLSqliteSwathHandler.cpp`: `readSwathWindows` (lines 39–53), `readMS1Spectra` (69–79), and `readSpectraForWindow` (98–108).

**Trigger:** Create an ordinary SPECTRUM table with ID rows 1 and -9223372036854775808, then add a virtual generated MSLEVEL column evaluating `abs(ID)`. The MS1 query can produce its first row before the next step reports SQLite integer overflow. Locks, corrupt pages or other execution failures provide additional possible triggers, but are not reproduced by this test.

**Issue:** Each accessor discards both initial and subsequent `sqlite3_step` return codes and treats a NULL column as end of results. Execution failure can consequently become a successful empty result or a successful prefix. Ignoring `sqlite3_finalize` also loses the saved statement error. This shares the unchecked-step mechanism of CPP-185 but concerns the independent public SWATH loops; CPP-196 separately covers a valid NULL row stopping the scan.

**Proposed C++ fix:** Accumulate rows only while step returns SQLITE_ROW, accept only SQLITE_DONE as successful completion, and throw for every other result. Use owned statement cleanup on all exits; never publish a partial result on query error. No upstream fix is claimed.

**Evidence:** Exact pinned source SHA-256 `8eb697bd44daee056ac3af662bab64b02f2ba82869a2f8f2bf188ae4dc277597`. Native test `sql_step_error_after_a_valid_row_is_not_mistaken_for_end_of_results` in [the SWATH tests](tests/mzml_sqlite_swath_handler.rs) executes a genuine SQLite step error after a valid row and verifies cleanup. Remote test evidence is recorded in [SWATH provenance](tests/data/mzml_sqlite_swath_provenance.json). The separate adapted C++ probe covers CPP-196/197/202 only; no executed C++ result is claimed for CPP-203.

**Rust handling:** Every fallible row advance propagates SQLite failure as `Error::Io`; the locally accumulated vector is discarded on failure and the transaction/connection close through ownership.

## CPP-204 — MSDataSqlConsumer owning raw pointer leaks during failed construction and can be shallow-copied

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h` (`MSDataSqlConsumer implicit copy constructor/assignment`, `handler_ ownership`); `src/openms/source/FORMAT/DATAACCESS/MSDataSqlConsumer.cpp:16` (`MSDataSqlConsumer::MSDataSqlConsumer`).

**Trigger:** Construct with a path whose createTables throws after handler allocation, or copy a successfully constructed consumer and destroy both copies.

**Issue:** Constructor allocates handler before throwing operations and no RAII owner releases it on failure. Public class does not disable copy; implicit copy duplicates owning pointer, leading to two destructors using/deleting it.

**Proposed C++ fix:** Use std::unique_ptr<MzMLSqliteHandler>, initialize safely, and delete copy operations or implement ownership-preserving copy semantics. No upstream fix is claimed.

**Evidence:** Header raw pointer; allocation line18 precedes reserve/create lines22–26; delete line42. Full header and source reviewed. No C++ allocation/leak probe. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; recommend owned Rust handler, no Clone for mutable consumer.


## CPP-205 — MSDataSqlConsumer destructor lets I/O exceptions escape noexcept destruction

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/DATAACCESS/MSDataSqlConsumer.cpp:29` (`MSDataSqlConsumer::~MSDataSqlConsumer`, `flush`, `writeRunLevelInformation`); `src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h` (`~MSDataSqlConsumer override`).

**Trigger:** Leave buffered records or unwritten RUN metadata, then encounter SQLite/disk failure at destruction.

**Issue:** Destructor calls throwing flush/write methods without catch. Its implicit noexcept override terminates the process on escape; explicit handler deletion is skipped.

**Proposed C++ fix:** Add an explicit checked finalize operation; keep destructor nonthrowing and make ownership RAII. Report persistent failures through finalize, not destruction. No upstream fix is claimed.

**Evidence:** Lines31 and38 call operations known to throw, inherited IMSDataConsumer destructor is nonthrowing. No executed termination reproduction. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; propose finish()->Result and nonthrowing Drop.


## CPP-206 — MSDataSqlConsumer changes buffered records to the next run ID

**Status and source:** Unconfirmed buffered-run ownership contract candidate at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/DATAACCESS/MSDataSqlConsumer.cpp:45` (`MSDataSqlConsumer::addRun`, `MSDataSqlConsumer::setRunId`, `MSDataSqlConsumer::flush`).

**Trigger:** Consume a spectrum under runA while below flush threshold, call setRunId(B) or addRun(...,B), then flush.

**Issue:** Buffers contain only records, not their original run IDs; flush writes using the handler current ID B. Records accepted before the run switch are reassigned to the later run. The header says the ID applies to subsequent writes, so deferred flush semantics may be intentional; review caller expectations before changing behavior.

**Proposed C++ fix:** Flush both buffers successfully before changing run ID, or retain per-record run ownership and flush grouped by original run. No upstream fix is claimed.

**Evidence:** Lines48/57 change handler ID with no flush; lines64/71 later call handler writers that bind its current run_id_. Header says setRunId applies to subsequent writes; buffered consumption ambiguity should be stated in fix rationale. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; require run-transition tests and deliberate per-run behavior.


## CPP-207 — SpectrumAccessSqMass unchecked view positions permit out-of-bounds access

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.cpp:27` (`SpectrumAccessSqMass(parent,indices)`, `getSpectrumById`, `getSpectrumMetaById`).

**Trigger:** Parent subset contains one entry; construct child with position -1, or call getSpectrumById/getSpectrumMetaById with -1 or >=subset length.

**Issue:** Nested constructor checks upper bound only, then indexes with negative int converted to size_t. Individual getters index the subset with no bounds check. Undefined behavior precedes checked SQL access.

**Proposed C++ fix:** Validate 0 <= id < visible count before every vector access; use checked at() or explicit validation for all constructors and accessors. No upstream fix is claimed.

**Evidence:** Lines43–45,75,106 show the unchecked accesses. Public header documents out-of-range parent selection throws. No UB execution attempted. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; use checked usize view positions/Result and validate SQL ID mapping.


## CPP-208 — SpectrumAccessSqMass bulk read does not preserve the configured view order or duplicates

**Status and source:** Source-reviewed defect; ordering additionally lacks sql guarantee at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.cpp:121` (`getAllSpectra`, `getSpectrumById`, `SpectrumAccessSqMass(handler,indices)`).

**Trigger:** Configure subset [1,0] or [1,1] and compare individual visible getters with getAllSpectra.

**Issue:** Individual access honors vector position, but bulk passes it to SQL IN without ORDER BY. SQL order is not guaranteed to match view order. Duplicate IDs collapse in SQL, then handler size mismatch throws despite individual accesses and reported subset count permitting duplicates.

**Proposed C++ fix:** Fetch unique underlying IDs once, map records by SQL ID, then reconstruct configured ordered view including repeated positions; alternatively reject duplicate views explicitly at construction if API changed. No upstream fix is claimed.

**Evidence:** Constructor stores input unchanged, getter uses sidx_[id], bulk forwards sidx_; handler uses WHERE ID IN(...) and compares unique returned count to original index count. No current C++ SQL-order execution. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; S1 sorted unique read API requires wrapper reordering and duplicate handling.


## CPP-209 — SpectrumAccessSqMass metadata index stays zero for every spectrum

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.cpp:114` (`getSpectrumMetaById`, `getAllSpectra`); `src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h:156` (`OSSpectrumMeta::index/default constructor`).

**Trigger:** Request metadata for visible spectrum1 or later, or bulk-read several spectra.

**Issue:** Source fills id,RT,ms_level but leaves index at default0, contradicting SpectrumMeta index documentation as zero-based consecutive spectrum-list index.

**Proposed C++ fix:** Assign the visible position to index in individual and bulk metadata results. No upstream fix is claimed.

**Evidence:** Assignments at114–117/153–156 omit index; constructor initializes index0. No executed reproduction. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; add exact metadata-index assertions for full and subset views.


## CPP-210 — SpectrumAccessSqMass class test repeats one out-of-range index instead of testing 50

**Status and source:** Source-reviewed test defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/tests/class_tests/openms/source/SpectrumAccessSqMass_test.cpp:97` (`parent-subset constructor START_SECTION`).

**Trigger:** The second out-of-range parent-subset test constructs indices2=[50].

**Issue:** TEST_EXCEPTION passes indices ([1]) again rather than indices2, so the alleged50 case is never tested.

**Proposed C++ fix:** Pass indices2 to the second constructor invocation; also add negative-index coverage. No upstream fix is claimed.

**Evidence:** Source lines97–103 show separate indices2 creation then old indices passed. No class-test execution. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; native independent tests should check both1 and50 and negative values.


## CPP-211 — SpectrumAccessSqMass installed example omits required handler run ID

**Status and source:** Source-reviewed documentation defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h` (`SpectrumAccessSqMass_example documentation`); `src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h` (`MzMLSqliteHandler(filename,run_id)`).

**Trigger:** Compile the installed example constructing MzMLSqliteHandler handler(file).

**Issue:** Pinned handler constructor requires two arguments and has no default run_id, so example cannot compile.

**Proposed C++ fix:** Pass handler(file,0) for reading. No upstream fix is claimed.

**Evidence:** Exact public header signature compared with example. No compiler run. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Native documentation should use its actual new(path,run_id) signature.


## CPP-212 — OpenSwath drift filtering dereferences missing or misaligned mobility arrays

**Status and source:** Source-reviewed defect at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h` (`ISpectrumAccess::filterByDrift`, `ISpectrumAccess::getSpectrumById(id,drift_start,drift_end)`); `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.cpp` (`SpectrumAccessSqMass::getSpectrumById`).

**Trigger:** Call base drift-filtered access on SqMass spectrum (which supplies only m/z and intensity), or filter a spectrum whose drift/intensity arrays are shorter than m/z.

**Issue:** All guards are commented out. filterByDrift dereferences null mobility pointer and advances drift/intensity iterators to m/z length without checking lengths. Every normal SqMass-returned spectrum lacks mobility arrays.

**Proposed C++ fix:** Validate all primary/mobility arrays and equal lengths; return a checked missing-mobility error and reject malformed ranges/arrays. No upstream fix is claimed.

**Evidence:** ISpectrumAccess header obtains nullable getDriftTimeArray then directly reads im_arr->data; source SqMass constructs only two primary arrays. No unsafe execution attempted. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Shared interface not implemented yet; missing mobility must be an explicit Result error, never a panic/null dereference.


## CPP-213 — SqMassFile store documentation hides unconditional replacement of an existing database

**Status and source:** Unconfirmed replacement-documentation candidate at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/include/OpenMS/FORMAT/SqMassFile.h` (`SqMassFile::store documentation`); `src/openms/source/FORMAT/SqMassFile.cpp:27` (`SqMassFile::store`).

**Trigger:** Store to an existing sqMass path expecting the documented create-if-necessary behavior.

**Issue:** Documentation says creating file/tables if necessary, but implementation always calls createTables which removes the existing database before write. Overwriting is common for store APIs; this wording alone does not prove append/preservation was promised. The candidate is the lack of explicit replacement/failure documentation, not overwrite semantics by themselves.

**Proposed C++ fix:** Explicitly document replacement/destructive behavior and failure consequences, or implement deliberate atomic replacement/append semantics. No upstream fix is claimed.

**Evidence:** store calls createTables unconditionally; handler createTables removes file. No current C++ write execution. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; stage whole output and atomically replace only after successful store if native chooses safer behavior.


## CPP-214 — MSDataSqlConsumer full metadata discards supplied experimental settings and addRun suppresses accumulated snapshot

**Status and source:** Source-reviewed loss; contract candidate at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/DATAACCESS/MSDataSqlConsumer.cpp:45` (`addRun`, `setExperimentalSettings`, `~MSDataSqlConsumer`, `consumeSpectrum`, `consumeChromatogram`).

**Trigger:** Use full_meta=true, send nondefault experimental settings and record-level metadata, call addRun, consume records, finish/destroy.

**Issue:** setExperimentalSettings is a no-op. addRun writes an empty snapshot and sets wrote_any_run, so destructor never writes accumulated peak_meta. Full read falls back from empty snapshot to SQL projection, losing metadata outside that projection.

**Proposed C++ fix:** Persist supplied settings and final per-run metadata snapshot; separate RUN existence from snapshot finalization, or document full_meta limitation for explicit/multiple-run workflow. No upstream fix is claimed.

**Evidence:** Line107 no-op; addRun empty meta lines49–52; destructor skip lines36–40; source handler fallback readExperiment empty snapshot. No execution; full-meta contract needs explicit policy review. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; preserve settings, distinguish add-run registration from final snapshot, and state multi-run compatibility limits.


## CPP-215 — MSDataSqlConsumer accepts negative signed buffer size before unsigned allocation

**Status and source:** Unconfirmed input-validation candidate at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/source/FORMAT/DATAACCESS/MSDataSqlConsumer.cpp:16` (`MSDataSqlConsumer constructor`); `src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h` (`constructor buffer_size argument`).

**Trigger:** Pass buffer_size=-1.

**Issue:** Signed argument becomes size_t maximum and reserve attempts an excessive allocation/throws rather than a clear invalid-size error; ctor raw-pointer leak separately confirmed. Zero immediately flushes rather than buffering.

**Proposed C++ fix:** Validate positive buffer size before allocation; decide and document whether0 is supported immediate-flush mode. No upstream fix is claimed.

**Evidence:** Source conversion/reserve reviewed; header has no documented negative semantics. No allocation probe. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; validate before any file replacement or allocation.


## CPP-216 — SpectrumAccessSqMass zero-width RT documentation differs from delegated first-at-or-after behavior

**Status and source:** Source-reviewed documentation mismatch candidate at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h` (`getSpectraByRT documentation`); `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.cpp:166` (`getSpectraByRT`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp` (`getSpectraIndicesbyRT`).

**Trigger:** Call getSpectraByRT(0.3,0) against source fixture with spectrumRTs0.2961 and0.4738.

**Issue:** Header describes exact interval[RT-delta,RT+delta], but delegated nonpositive branch returns firstRT>=target, potentially0.4738. This behavior is needed by base nearest-spectrum convenience.

**Proposed C++ fix:** Document special delta0 semantics and preserve them deliberately, or change base nearest-search contract together. No upstream fix is claimed.

**Evidence:** Source header and delegation read; source handler source test confirms nonpositive behavior for delta-.1 but no exactzero accessor test. No current accessor execution. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Not implemented yet; native RT boundary tests and explicit zero semantics required.


## CPP-217 — OpenSwath getMultipleSpectra count parameter can return more than requested

**Status and source:** Source-reviewed algorithm/documentation candidate at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Affected files and functions:** `src/openswathalgo/source/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.cpp` (`ISpectrumAccess::getMultipleSpectra both overloads`).

**Trigger:** Request n=0 or n=2 around an interior valid closest spectrum.

**Issue:** Always pushes closest once; then adds left+right for i<=n/2. n0 yields1; n2 canyield3 despite comment describing n as sequence length. Negative requests also return one.

**Proposed C++ fix:** Validate positive odd n and document odd-window semantics, or cap selection to requested count and define even tie policy. No upstream fix is claimed.

**Evidence:** Both full implementation overloads read; no class-test expectation checked yet and no execution. Exact source hashes and review scope are retained in [the S2 source-review manifest](tests/data/sqlite_s2_review.json). No executed C++ reproduction or S2 native test is claimed.

**Rust handling:** Shared interface not implemented yet; choose explicit bounded count contract after interface test audit.

## CPP-218 — Positive-accuracy Numpress encoding destroys one- and two-point coordinate arrays

**Status and source:** Reproduced through the exact standalone raw C++ codec at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The MSNumpressCoder/sqMass calling paths are source-reviewed, not executed as a full SDK.

**Affected files and functions:** `src/openms/source/FORMAT/MSNUMPRESS/MSNumpress.cpp`: `optimalLinearFixedPointMass` (237–245), `optimalLinearFixedPoint` (262 onward), `encodeLinear` (326 onward), and `decodeLinear` (443, 458); `src/openms/source/FORMAT/MSNumpressCoder.cpp` (`encodeNPRaw`, 131–151); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp` (`writeSpectra` and `writeChromatograms` Numpress configuration, 1078–1085 and 1319–1326).

**Trigger:** Encode a nonempty array containing one or two coordinates, such as [100], [100,101] or chromatogram RTs [1,2], using positive target accuracy 0.0001 or 0.05. The sqMass writer uses this configuration and disables codec error verification.

**Issue:** The accuracy estimator returns zero for fewer than three values, incorrectly commenting that the first two points are encoded as floats. They are actually integer-quantized by the factor, so zero destroys the input. Decoding divides zero by zero and yields NaN. MSNumpressCoder falls back only for a negative estimated factor, so zero is not corrected. Using the ordinary estimator alone is also insufficient: its factor is infinite for short all-zero arrays.

**Proposed C++ fix:** Select a positive finite factor for nonempty short arrays, respecting accuracy and integer range, or fall back to an explicitly supported lossless storage mode. Treat short all-zero input deliberately; reject nonfinite/zero factors where raw nonempty linear encoding or decoding cannot interpret them. No upstream fix is claimed.

**Evidence:** A standalone driver compiled the exact pinned raw translation unit/header on kim with GCC 13.3/C++20. All 42 cases were checked: 12 short positive-accuracy cases decoded to NaN, four short all-zero ordinary estimates were infinite and skipped before encoding, and 26 finite controls passed (including six explicit-factor-one zero controls). A separate sanitizer diagnostic stopped at the invalid floating-point-to-integer conversion for the infinite zero estimate; no resulting encoded/decoded values are used. Sources, driver, binary and logs are hashed in [handler provenance](tests/data/mzml_sqlite_handler_provenance.json) and retained under `../oracle/sqlite-short-numpress-s1-probe/`. This is tier 2 raw-codec execution; no full handler or upstream class-test run is claimed.

**Rust handling:** The sqMass handler writes one- and two-point coordinate arrays using source-supported lossless code 1, including all-zero inputs. Intensity SLOF behavior remains unchanged. The public raw codec's documented source semantics remain separate. Native short-array tests check stored compression tags and finite exact coordinate round trips.

**Additional source-reviewed trigger (resume):** Three coordinates `[1e12,1e12,1e12]` also yield a zero maximal linear factor after the positive-accuracy estimator falls back. This case was not in the executed 42-case C++ probe. The native handler rejects the invalid factor atomically and permits an explicit lossless retry; the new regression covers spectra and chromatograms. The executed probe count above remains unchanged.

## CPP-219 — Chromatogram precursor reload drops supplemental activation metadata

**Status and source:** Source-reviewed read/write mismatch at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; no executed C++ reproduction.

**Affected files and functions:** `src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, spectrum activation CV branch (1990–2009), chromatogram activation CV branch (2033–2158), shared `writePrecursor_` (4736 onward).

**Trigger:** Store a chromatogram precursor carrying supplemental beam-type collision-induced dissociation (MS:1002678), supplemental collision-induced dissociation (MS:1002679), or supplemental collision energy (MS:1002680), then reload it.

**Issue:** The shared writer emits these values as CV parameters. The spectrum reader handles them, but the chromatogram branch has no corresponding routes, losing supplemental metadata and its derived activation methods on reload.

**Proposed C++ fix:** Share the activation CV decoder between spectrum and chromatogram precursors, preserving supplemental values and derived methods consistently. Add chromatogram round-trip cases for all three accessions. No upstream fix is claimed.

**Evidence:** Exact source hashes, CV definitions and native test mapping are recorded in [activation provenance](tests/data/mzml_precursor_activation_provenance.json). This conclusion follows source review; native regressions cannot establish executed C++ behavior.

**Rust handling:** The common precursor decoder handles these accessions for both spectrum and chromatogram records. Native writing preserves the documented typed metadata subset, with explicit handling of supplemental-method combinations; see [support](docs/MZML_PRECURSOR_ACTIVATION_SUPPORT.md).

## CPP-220 — Negative initial linear-Numpress coordinates can wrap during sqMass writing

**Status and source:** Unconfirmed C++ defect candidate from source review and independent arithmetic at revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; no executed C++ reproduction for this trigger.

**Affected files and functions:** `src/openms/source/FORMAT/MSNUMPRESS/MSNumpress.cpp` (`encodeLinear`, first two values at 326–341, and `decodeLinear` at 443–458); `src/openms/source/FORMAT/MSNumpressCoder.cpp` (`encodeNPRaw`); `src/openms/source/FORMAT/HANDLERS/MzMLSqliteHandler.cpp` (`writeSpectra`, `writeChromatograms`).

**Trigger:** Write finite coordinate values `[-100,-99,-98]` with positive accuracy 0.0001 and lossy compression.

**Issue:** The estimated factor is positive, but the first two negative quantized integers are serialized using their low 32 bits and decoded as unsigned values. They cannot recover the original negative coordinates. The sqMass writer disables error verification and can therefore commit a shifted array. This candidate concerns missing storage-domain validation; the raw codec's supported input domain should also be clarified.

**Proposed C++ fix:** Require the first two truncated quantized values to fit `[0, UINT32_MAX]` before committing, or explicitly choose a supported lossless representation. No upstream fix is claimed.

**Evidence:** Source hashes and the separately executed native regression `invalid_lossy_coordinate_quantization_rolls_back_without_advancing_ids` are recorded in [handler provenance](tests/data/mzml_sqlite_handler_provenance.json). The earlier C++ short-array probe does not reproduce this negative-coordinate case.

**Rust handling:** The storage adapter checks the initial quantized values, returns an error without changing rows or counters, and allows explicit lossless retry. The public raw codec's documented source behavior is unchanged.

## CPP-221 — Gumbel result declares eval without defining it

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed missing public definition; C++ link failure inferred, not executed.

**Affected files/functions:** `src/openms/include/OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h:51`, `GumbelDistributionFitter::GumbelDistributionFitResult::eval(double) const`; `src/openms/source/MATH/STATISTICS/GumbelDistributionFitter.cpp`.

**Trigger:** A client constructs `GumbelDistributionFitter::GumbelDistributionFitResult(0.0, 1.0)` and calls `eval(0.0)`.

**Issue:** The public header declares the non-inline member but the implementation defines only `log_eval_no_normalize`, constructors/configuration and fitting. A source-tree search found no definition of this class's `eval`. A conventional linked client therefore references an undefined symbol; this is distinct from a member omitted by the Rust port.

**Proposed C++ fix:** Define and export the result member using the existing residual's Gumbel density `(z * exp(-z)) / b`, where `z = exp((a-x)/b)`, and add a public-client link/evaluation test. Ensure the symbol is visible in shared-library builds, since the nested result currently lacks its own export annotation. No fix has been applied upstream.

**Evidence:** Exact declaration and absence of a definition in the pinned source search; `GumbelDistributionFitter.cpp:56–66` supplies the already-used residual model. `tests/data/distribution_fitters_provenance.json` records the original finding. No executed C++ reproduction.

**Rust handling:** `src/math/fitters/gumbel.rs::GumbelDistributionFitResult::eval` exposes that residual model as a native extension with finite/positive-scale checks; existing unit tests compare its peak and log-density consistency. The module preamble says the declaration is “not ported,” while the method documentation accurately explains the extension. Rust's evaluation is present; `fit_weighted` on this least-squares class remains absent.

## CPP-222 — Gumbel least-squares fitter advertises an undefined weighted fit

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed missing public definition; C++ link failure inferred, not executed.

**Affected files/functions:** `src/openms/include/OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h:72–81`, `GumbelDistributionFitter::fitWeighted`; corresponding `.cpp`.

**Trigger:** Invoke `GumbelDistributionFitter().fitWeighted(x, w)` with ordinary equally sized nonempty samples and weights.

**Issue:** The public declaration and Doxygen promise a weighted histogram followed by a fit, but no method definition exists. The identically named method on **GumbelMaxLikelihoodFitter** is a different class and scientific operation. The source Gumbel class test calls that other fitter and therefore does not exercise this missing public symbol.

**Proposed C++ fix:** Implement the documented weighted-histogram fitting contract, including explicit binning and weight validation, and test it through this class; alternatively remove/deprecate the unsupported declaration through an explicit API decision. Do not silently substitute maximum likelihood for the documented histogram fit. No upstream fix claimed.

**Evidence:** Header/implementation/source symbol search and the source class test's `gmlf.fitWeighted` calls; `tests/data/distribution_fitters_provenance.json`. No C++ link probe executed.

**Rust handling:** `src/math/fitters/gumbel.rs` explicitly documents this missing source operation; the separate `gumbel_max_likelihood` module supplies the separately defined maximum-likelihood class. The weighted-histogram public contract remains unimplemented, not fulfilled by that other module.

## CPP-223 — Distribution-fitter documentation names nonexistent gnuplot members

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed documentation defect.

**Affected files/functions:** Class Doxygen in `src/openms/include/OpenMS/MATH/STATISTICS/GaussFitter.h:29–30`, `GumbelDistributionFitter.h:28–29`, and `GumbelMaxLikelihoodFitter.h:26–27`; claimed `getGnuplotFormula` member.

**Trigger:** Follow the class documentation and call `getGnuplotFormula()` after a fit.

**Issue:** None of those three classes declares the advertised member. Client compilation would fail at member lookup. The source may print a formula under verbose compilation flags; that is not a public getter. Other unrelated trace-fitters' real methods do not satisfy these class promises.

**Proposed C++ fix:** Remove the obsolete promise and document constructing a formula from returned parameters, or implement a deliberate public getter with tests. Review debug-only call sites in `ANALYSIS/ID/IDDecoyProbability.cpp`, which still mention the missing methods, when enabling that optional debug macro. No enabled-debug compile was run here.

**Evidence:** Complete public headers and implementation search; distribution-fitters provenance. No compilation attempted.

**Rust handling:** The native fitter APIs expose parameter/result values and do not invent the missing getter; `docs/DISTRIBUTION_FITTERS_SUPPORT.md` explicitly records the stale claim. No native gnuplot getter or executed C++ compatibility is claimed.

## CPP-224 — Gumbel maximum-likelihood fitting reads beyond a short weight vector

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed memory-safety defect; out-of-bounds access not executed.

**Affected files/functions:** `src/openms/source/MATH/STATISTICS/GumbelMaxLikelihoodFitter.cpp:47–66`, anonymous `GumbelDistributionFunctor::operator()`, and public `fitWeighted` at 74–83; corresponding public header.

**Trigger:** `fitWeighted({1.0, 2.0}, {1.0})` with the ordinary finite default initial parameters.

**Issue:** The objective iterates until `m_data.cend()` while dereferencing and incrementing `wit` without checking the weight length. The second sample dereferences the weight end iterator. Public `fitWeighted` validates neither length before passing the vectors to the optimizer. Longer weights are ignored rather than diagnosed; only shorter weights cause this specific read past end.

**Proposed C++ fix:** Require `x.size() == w.size()` before constructing the objective, return a documented fitting/argument error otherwise, and add short/long-weight regressions. No upstream change claimed.

**Evidence:** Direct loop and public entry-point source trace, recorded in distribution-fitters provenance. No sanitizer or C++ execution.

**Rust handling:** `src/math/fitters/gumbel_max_likelihood.rs::fit_weighted` rejects unequal lengths with `Error::InvalidValue` before iteration and leaves initial parameters unchanged. Existing `mismatched_weights_are_refused` unit coverage is present; not rerun here.

## CPP-225 — Memory-usage delta overwrites its minus sign

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed reporting defect.

**Affected file/function:** `src/openms/source/SYSTEM/SysInfo.cpp:324–334`, `SysInfo::MemUsage::diff_str_`, used by `delta`.

**Trigger:** A recorded working set falls from 4096 KiB to 1024 KiB before `delta("release")` formats the difference.

**Issue:** The function appends `"-"` when memory decreases, then overwrites the entire string with the absolute whole-MiB magnitude plus `" MB"`. The result is `"3 MB"` instead of `"-3 MB"`, so a decrease is indistinguishable from growth. The trigger is an independently derived example, not a sampled process reproduction.

**Proposed C++ fix:** Append the magnitude to the existing sign (`s += ...`) or construct sign and magnitude together. Add a deterministic unit test for decreasing sampled counters, avoiding dependence on allocator release behavior. No upstream fix claimed.

**Evidence:** Exact assignment following the sign branch; `tests/data/sys_info_provenance.json`. Source-only.

**Rust handling:** `src/system/sys_info.rs::difference_string` computes a signed i128 difference and preserves the minus sign. Existing unit coverage and `tests/sys_info.rs::a_negative_delta_keeps_its_sign` assert the `-3 MB` report with deterministic sample values. No fresh test execution here.

## CPP-226 — Update notification announces the local version instead of the offered update

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed user-visible reporting defect.

**Affected file/function:** `src/openms/source/SYSTEM/UpdateCheck.cpp:136–142`, `UpdateCheck::run` update-available notice.

**Trigger:** Running library/tool version 2.0.0, a successfully parsed server response 3.0.0, and a check that reaches the update-available branch.

**Issue:** The branch correctly determines that `server_version` is newer, but concatenates the local `version` argument into “Version ... is available”. The user is told that 2.0.0 is available rather than 3.0.0. The local tool version and running library version can differ; neither is necessarily the offered update.

**Proposed C++ fix:** Render the validated server version in the notice, retaining the tool name and URL. Add a fake-response test with different local/server values. No upstream fix claimed.

**Evidence:** Direct source branch/string construction; `tests/data/network_provenance.json`. Its phrase “reproduced verbatim” concerns native behavior, not an executed C++ oracle.

**Rust handling:** `src/system/update_check.rs::run` deliberately preserves the local-version notice, documents it, and returns the server response for callers. Existing `tests/update_check.rs` assertions pin local-version messages. This defect is **not corrected in the native notice**; neither a fix nor C++ execution is claimed.

## CPP-227 — Download filename selection does not prevent concurrent overwrite

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed race violating the documented no-overwrite promise; no concurrent C++ reproduction executed.

**Affected files/functions:** `src/openms/source/SYSTEM/Network.cpp:36–45`, `saveFileName_`, and 60–62, `Network::downloadFile`; `src/openms/include/OpenMS/SYSTEM/Network.h:36–39`.

**Trigger:** Two successful downloads choose the same absent basename. Both finish the `exists` check before either opens the destination; one opens/writes it, then the other opens the same path.

**Issue:** Existence checking and file creation are separate. `std::ofstream(filename, std::ios::binary)` opens for output and truncates an existing file. The second writer can destroy or replace the first result despite the unconditional public promise that existing files are never overwritten. Returning/logging the chosen path does not prevent or reliably detect this race.

**Proposed C++ fix:** Reserve each candidate atomically with exclusive creation, retry an already-existing candidate, and write through the reserved handle. Preserve the guarantee for dangling symlinks and concurrent creators too; do not add another pre-open existence check. Add a controlled concurrent-candidate test. No upstream change claimed.

**Evidence:** Header promise plus source check/open sequence; `tests/data/network_provenance.json`. Interleaving is source-derived, not executed.

**Rust handling:** `src/system/network.rs::download_file_with` also calls `save_file_name` then `fs::File::create`, retaining the race. Its prose currently repeats the no-overwrite promise, while the provenance correctly admits the race. On subsequent write failure, native cleanup can also remove that raced destination. This needs a native correction/documentation disposition; it is **not already fixed** by returning `PathBuf`.

## CPP-228 — Squaring Gumbel negative log likelihood can change the optimum

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** **Unconfirmed scientific-defect candidate.** Objective structure is source-reviewed; an actual failed fit versus a valid independent MLE remains to be demonstrated.

**Affected file/functions:** `src/openms/source/MATH/STATISTICS/GumbelMaxLikelihoodFitter.cpp:53–66`, objective residuals, and 74–83, Levenberg–Marquardt invocation.

**Candidate trigger:** A narrow-scale, finite, nondegenerate sample (for example two separated observations near zero with positive weights) whose weighted negative log likelihood can be negative around its optimum and zero on another parameter contour. Exact optimizer outcome depends on initial parameters and termination and has not been executed here.

**Concern:** The residual vector is `[NLL, 0]`; least squares minimizes `NLL²`. This has the same ordering as NLL while NLL stays nonnegative, but favors a zero contour over a negative minimum when one is reachable. Continuous probability densities can exceed 1, so negative log likelihood is not intrinsically nonnegative. The existing provenance's definite attraction statement overstates demonstrated behavior without an actual fitted counterexample.

**Proposed investigation/fix:** Compare a narrow-scale sample against an independently implemented scalar-NLL optimizer and check stationarity/likelihood. If confirmed, minimize NLL directly using a suitable solver with a positive-scale parameterization rather than squaring it. Keep any behavior change explicit because source-matching downstream results may move.

**Evidence:** Source residual construction and independent mathematical reasoning only. No C++ fit, numerical experiment or new native test executed in this reconciliation.

**Rust handling:** `src/math/fitters/gumbel_max_likelihood.rs` deliberately retains the `[NLL, 0]` residual with the native LM solver. Its API documents driving NLL toward zero. The suspected optimization issue remains shared, not fixed.

## CPP-229 — Download collision-suffix counter lacks an overflow guard

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** **Unconfirmed extreme-input robustness candidate.** Unbounded suffix probing and signed increment are source-reviewed; overflow reachability under a real filesystem has not been reproduced.

**Affected file/function:** `src/openms/source/SYSTEM/Network.cpp:41–45`, `saveFileName_`.

**Candidate trigger:** The basename and every successive suffix through `INT_MAX` already exist, or concurrent creators keep occupying candidates until the signed counter reaches its limit.

**Concern:** A signed `int` is incremented without checking overflow; if reached, increment beyond `INT_MAX` is undefined behavior. Even below that extreme, scan work scales with existing collisions without a caller-visible bound. A normal finite directory with a gap terminates: the previous wording that the loop simply “does not terminate” should not be copied as a general claim.

**Proposed C++ fix:** Add a documented finite collision limit with an I/O error, using checked counter arithmetic; combine with atomic exclusive creation from the separate race entry.

**Evidence:** Source loop only; no huge-directory experiment or C++ overflow reproduction.

**Rust handling:** `src/system/network.rs::save_file_name` checks suffixes 0 through `MAX_NAME_SUFFIX` (10,000) inclusively and then returns `Error::InvalidValue`. This bounds suffix probing; it does not fix the separate selection/open race.

## CPP-230 — FuzzyStringComparator reports the last above-one ratio as the maximum

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (`src/testframework`). Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`, whose comparator source is byte-identical.

**Status:** Executed. Report defect; the verdict is unaffected.

**Affected file/function:** `src/testframework/source/CONCEPT/FuzzyStringComparator.cpp:665–683`, `compareLines_`.

**Trigger:** Numbers that differ by ratios above one but within the relative tolerance, reported at verbose level 2 or higher.

**Issue:** The recorded lines of the maximum are replaced whenever `ratio > ratio_max_`, but `ratio_max_` itself is raised only on the failure branch, just before `reportFailure_`. On a passing comparison it stays at 1, so every above-one ratio replaces the recorded lines. A PASSED report therefore prints `relative_max: 1`, and "Maximum relative error was attained at these lines" names the last line with any ratio above one, not the line with the largest ratio.

**Proposed C++ fix:** Raise `ratio_max_` together with the recorded lines whenever `ratio > ratio_max_`, independently of the tolerance test.

**Evidence:** Oracle case `rep_success_last_ratio_line` in `../oracle/fuzzy-string-comparator/manifest.json`; `tests/data/fuzzy_string_comparator_provenance.json`.

**Rust handling:** `tests/support/fuzzy_string_comparator.rs` is test support for FuzzyDiff parity and reproduces the source report; every oracle case except the hexadecimal one of CPP-233 matches in verdict and log bytes.

## CPP-231 — A reused FuzzyStringComparator keeps its raised ratio maximum

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (`src/testframework`). Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/testframework/source/CONCEPT/FuzzyStringComparator.cpp:805`, `compareStreams`; the raise at 678 in `compareLines_`.

**Trigger:** One instance compares `1` against `1.5`, which fails, then `1` against `1.2` with a relative tolerance of 1.01 and an absolute tolerance of 0.

**Issue:** `compareStreams` resets only `is_status_success_`. The failed comparison leaves `ratio_max_` at 1.5, and `compareLines_` tests the tolerance only when a ratio exceeds `ratio_max_`, so the second comparison accepts 1.2 although it is outside the tolerance.

**Proposed C++ fix:** Reset `ratio_max_`, `absdiff_max_` and the recorded maximum lines at the start of every comparison.

**Evidence:** Oracle case `reuse_ratio_max_carries` in `../oracle/fuzzy-string-comparator/manifest.json`.

**Rust handling:** Reproduced: the test support carries the maxima over to the next comparison on the same instance, as documented in `docs/FUZZY_STRING_COMPARATOR_SUPPORT.md`.

## CPP-232 — A negative ratio that underflows to -0.0 passes the sign and ratio tests

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (`src/testframework`). Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/testframework/source/CONCEPT/FuzzyStringComparator.cpp:645–660`, `compareLines_`.

**Trigger:** `-1e-300` compared with `1e300`.

**Issue:** The quotient `element_1_.number / element_2_.number` underflows to `-0.0`, so `ratio < 0` is false and the different-signs failure is skipped. `ratio < 1` then takes the reciprocal, `-inf`, which never exceeds `ratio_max_`, and the pair is accepted.

**Proposed C++ fix:** Compare the operands' signs (for example with `std::signbit`) instead of the quotient's, and refuse a zero or infinite quotient of two non-zero numbers.

**Evidence:** Oracle case `num_ratio_underflow_negative_zero` in `../oracle/fuzzy-string-comparator/manifest.json`.

**Rust handling:** Reproduced in the test support; the oracle case matches.

## CPP-233 — The libc++ number fallback accepts hexadecimal floats, so verdicts depend on the standard library

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (`src/testframework`). Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85` (AppleClang, libc++).

**Status:** Executed on libc++; the libstdc++ side comes from a verifier probe with g++ 13.3, not retained.

**Affected file/function:** `src/testframework/source/CONCEPT/FuzzyStringComparator.cpp:160–211`, `fromCharsFloat`, compiled when `_LIBCPP_VERSION` is defined.

**Trigger:** `0x10` compared with `16`; also `0X1P-2` with `0.25`.

**Issue:** On libc++ builds number tokens are parsed with `strtod`, which accepts hexadecimal floats that `std::from_chars` in general format rejects, so `0x10` equals `16` on Apple builds only. The fallback also writes `strtod`'s result into the reset number of a letter element, which changes the failure report. Its comment assumes `std::from_chars` reports `result_out_of_range` only on overflow, but libstdc++'s `std::from_chars` also rejects `1e-400`, so underflow is platform-dependent too. Identical inputs therefore get different verdicts and reports by platform.

**Proposed C++ fix:** Make the fallback reject what `std::from_chars` rejects (hexadecimal prefixes, and `nan(...)` sequences it does not consume), decide underflow explicitly on both paths, and leave a letter element's number untouched.

**Evidence:** Oracle case `tok_hex_vs_decimal` in `../oracle/fuzzy-string-comparator/manifest.json`; the C3-FUZZY verifier's adversarial cases `adv_hex_capital_prefix_letters` and `adv_nan_dash_inside_plus`.

**Rust handling:** The test support follows the `std::from_chars` contract the source states: no hexadecimal floats, underflow accepted. `tok_hex_vs_decimal` is asserted as a known divergence from the macOS oracle.

## CPP-234 — A missing second input is reported as the first input file

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (`src/testframework`). Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed. Diagnostic defect.

**Affected file/function:** `src/testframework/source/CONCEPT/FuzzyStringComparator.cpp:891–898`, `openInputFileStream_` (message at 896).

**Trigger:** `compareFiles` with a readable first file and a missing second file.

**Issue:** `openInputFileStream_` serves both inputs but always logs "Error opening first input file '<name>'".

**Proposed C++ fix:** Pass the input's position or label into `openInputFileStream_`.

**Evidence:** Oracle case `file_missing_second` in `../oracle/fuzzy-string-comparator/manifest.json`.

**Rust handling:** Reproduced: both open failures report "Error opening first input file".

## CPP-235 — A NaN relative tolerance accepts every ratio

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (`src/testframework`). Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/testframework/source/CONCEPT/FuzzyStringComparator.cpp:297–305`, `setAcceptableRelative`.

**Trigger:** `setAcceptableRelative(NaN)`.

**Issue:** `ratio_max_allowed_ < 1.0` is false for NaN, so NaN is stored unchecked, and `ratio > ratio_max_allowed_` is then false for every ratio: no relative difference can fail.

**Proposed C++ fix:** Refuse a non-finite or non-positive tolerance in the setter.

**Evidence:** Oracle case `num_ratio_nan_setter` in `../oracle/fuzzy-string-comparator/manifest.json`.

**Rust handling:** Reproduced in the test support; the oracle case matches.

## CPP-236 — FuzzyDiff -sort defeats the same-file check

**Source revision:** topp `174b576e244e100f2345ca57a8e79aaa607156df` (`src/FuzzyDiff.cpp`) with core `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (`src/testframework`). Executed with the product-SDK FuzzyDiff.

**Status:** Executed.

**Affected file/function:** `src/FuzzyDiff.cpp:136–196`, the `-sort` branch of `main_`; `src/testframework/source/CONCEPT/FuzzyStringComparator.cpp:864–866`, `compareFiles`.

**Trigger:** `FuzzyDiff -sort -in1 a.tsv -in2 a.tsv`.

**Issue:** With `-sort`, FuzzyDiff writes both inputs to temporary files with unique names and compares those, so `compareFiles`'s "first and second input file have the same name. That's cheating!" check never fires, and a file compared with itself passes.

**Proposed C++ fix:** Check the original input names before sorting.

**Evidence:** Oracle case `fd_sort_same_file` in `../oracle/fuzzy-string-comparator/manifest.json`.

**Rust handling:** `fuzzy_diff` sorts in memory with the same verdict and exit code, so the self-comparison passes as in the source.

## CPP-237 — FeatureFinderCentroided's intensity filter keeps zero and subnormal intensities

**Source revision:** topp `174b576e244e100f2345ca57a8e79aaa607156df` (`src/FeatureFinderCentroided.cpp`) with core `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (`DPosition.h`). Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/FeatureFinderCentroided.cpp:182–184`, `main_`; `src/openms/include/OpenMS/DATASTRUCTURES/DPosition.h`, which has no `std::numeric_limits` specialisation (`minPositive` is at 326–329).

**Trigger:** Any input spectrum with zero or subnormal intensities.

**Issue:** The comment says "filter out zero (and negative) intensities", but the range starts at `std::numeric_limits<DPosition<1>>::min()`. Without a specialisation the primary template returns a value-initialised `DPosition`, which is 0, so zero and subnormal intensities pass and only negative ones are dropped.

**Proposed C++ fix:** Start the range at `RP_TYPE::minPositive()`.

**Evidence:** `../oracle/a3-format-io/manifest.json`: `intensity_bounds.mzML` keeps 0.0, 1e-310 and `DBL_MIN` and drops only -1.0; `tests/data/mzml_mobility_provenance.json`.

**Rust handling:** No FeatureFinderCentroided wrapper exists yet. `FileHandler::load_experiment_with_options` applies the range a caller passes; a source-faithful wrapper (packages C5 and B10) must pass `[0.0, f64::MAX)`. FeatureFinderCentroided_1 loads 112 spectra and 3084 peaks with either lower bound.

## CPP-238 — A FAIMS voltage of -1 V is indistinguishable from an unset drift time

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed for the reader and format query; the FAIMSHelper consequence is source-reviewed.

**Affected file/function:** `src/openms/include/OpenMS/IONMOBILITY/IMTypes.h:93`, `DRIFTTIME_NOT_SET = -1.0`, and its readers, including `src/openms/source/IONMOBILITY/FAIMSHelper.cpp:50`.

**Trigger:** A spectrum with `MS:1001581` (FAIMS compensation voltage) value `-1` volt.

**Issue:** The sentinel -1 is a valid compensation voltage. Such a spectrum loads with drift time -1 and unit FAIMS_CV, but `determineIMFormat` reports no ion mobility, and `FAIMSHelper::getCompensationVoltages` erases -1 V as a missing voltage and warns. The writer still writes the value, because it tests the unit before the value.

**Proposed C++ fix:** Represent "not set" apart from the value range, for example by the unit `NONE` alone or an optional drift time.

**Evidence:** `../oracle/a3-format-io/manifest.json`, `scan_mobility.mzML` spectrum 7.

**Rust handling:** The port keeps the source's -1 sentinel (`metadata::ImTypes::DRIFTTIME_NOT_SET`), so it shares the collision; `tests/mzml_mobility.rs` matches the oracle row.

## CPP-239 — The precursor drift-time writer has no FAIMS case

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; no C++ write was executed.

**Affected file/function:** `src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp:4607–4628`, the selected-ion writer; the reader maps the term correctly at 1864.

**Trigger:** A precursor whose drift-time unit is `FAIMS_COMPENSATION_VOLTAGE`.

**Issue:** The switch handles `MILLISECOND`, `VSSC` and `CCS`. A FAIMS voltage falls into `default`, which warns "Precursor drift time unit not set, assume milliseconds" and writes `MS:1002476` in milliseconds, so the voltage reads back as a millisecond drift time.

**Proposed C++ fix:** Add a FAIMS case that writes `MS:1001581` with `UO:0000218`, and keep the millisecond fallback for `NONE` only.

**Evidence:** Source review; `tests/data/mzml_mobility_provenance.json`.

**Rust handling:** The precursor writer takes its term from the shared mobility table in `src/format/mzml_precursor.rs`, which maps the FAIMS unit to `MS:1001581` with `UO:0000218`, so it does not fall back to milliseconds.

## CPP-240 — Upstream FAIMS fixtures spell the volt unit UO:000218

**Source revision:** test-data `0cb15f23fccc6ea196bfafcfbbf020958f36c3a3`.

**Status:** Source-reviewed fixture defect.

**Affected files:**
- `topp/FAIMS_CV-60C_V-45_Interleaved.mzML`, lines 324, 364, 404, 444, 484, 524, 564, 604, 644, 684, 724 and 764: all 12 `MS:1001581` cvParams in the file, which has 12 spectra.
- `topp/FAIMS_test_data.mzML`, lines 171 and 213: both `MS:1001581` cvParams in the file, which has 2 spectra.

**Issue:** Every FAIMS compensation voltage cvParam in both files spells the `unitAccession` `UO:000218`, one digit short; neither file contains `UO:0000218`. The volt term is `UO:0000218`, which the pinned writer emits. A reader that checks unit accessions rejects or ignores the unit.

**Proposed fix:** Correct the accession in both fixtures, or regenerate them with the current writer.

**Evidence:** In the test-data package, `git show 0cb15f2:topp/FAIMS_CV-60C_V-45_Interleaved.mzML | grep -n 'unitAccession="UO:000218"'` lists the 12 lines and `git show 0cb15f2:topp/FAIMS_test_data.mzML | grep -n 'unitAccession="UO:000218"'` the 2; `grep -c 'accession="MS:1001581"'` gives the same 12 and 2.

**Rust handling:** The mzML reader accepts `UO:000218` as volts for the FAIMS voltage only (native difference 3 in `docs/MZML_MOBILITY_SUPPORT.md`); the writer emits `UO:0000218`.

## CPP-241 — computeIntensityProfile dereferences begin() of an empty MassTraces

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed undefined behaviour; not executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.cpp:196–198`, `MassTraces::computeIntensityProfile`.

**Trigger:** An empty `MassTraces`.

**Issue:** `trace_it = this->begin()` is dereferenced (`trace_it->peaks`) and incremented without a check for an empty collection. The shipped callers pass non-empty traces.

**Proposed C++ fix:** Return early when `this->empty()`.

**Evidence:** Source review; `tests/data/feature_finder_picked_helper_structs_provenance.json`.

**Rust handling:** `MassTraces::intensity_profile` returns an empty profile.

## CPP-242 — computeIntensityProfile never terminates on a NaN retention time

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Executed (2026-09-16).

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.cpp:210–236`, `MassTraces::computeIntensityProfile`.

**Trigger:** A NaN retention time in a later trace compared against a profile entry, or a NaN entry copied from the first trace that a later trace reaches.

**Issue:** None of `>`, `<` and `==` holds for NaN, so no branch advances either iterator and the `while` loop never ends.

**Proposed C++ fix:** Treat an unordered comparison as an error, or advance the profile iterator.

**Evidence:** FeatureFinderCentroided_1 with the retention time of scan 50 set to NaN, through `FeatureFinderAlgorithmPicked::run` on the Release build (`nonfinite_stage` cases `rt_nan_mid`, `rt_nan_mid_bins3`, `rt_nan_mid_unsorted`), and with scan 0's retention time NaN in an unsorted input (`sort_mobility_stage` case `v3_rt_nan_first_unsorted`): the runs did not return within 30 s (killed, twice each) and a gdb stack sample shows `MassTraces::computeIntensityProfile` under `GaussTraceFitter::setInitialParameters_` in the seed loop. `extendMassTrace_` adds the NaN-RT peak because its NaN overall score is not below 0.01 (`.cpp:1535`). A NaN RT at the first or last scan of a sorted input never joins a trace (8 features); in an unsorted input the introsort moves it (`v3_rt_nan_second_unsorted` and `v3_rt_nan_last_unsorted`: 8 features). Round 4 adds: with `write_debug` (`debug:pseudo_rt_shift 500`), scan 50's or scan 20's retention time NaN (Gaussian and EGH): killed after 30 s, twice each, identical; the process had flushed 1,016,234 or 1,065,423 bytes of `debug/log.txt` and written the files of plots 0-2 or 0-9 (`../oracle/ffap-complete-fix4`, `fix4_vfi`). Also `tests/data/feature_finder_picked_helper_structs_provenance.json`.

**Rust handling:** Returns `Error::InvalidValue`; a NaN that is only copied or appended passes through as in the source. Refused at exactly that merge (`MassTraces::intensity_profile`); the endless runs are replayed as refusals; a debug run records a `NeverReturns` termination with the hanging seed's plot number after the seed's log lines, so a caller writes only the flushed log, which equals the executed file byte for byte, as do the seed map and the feature files (`a_seed_loop_that_never_returns_keeps_what_the_executed_process_had_written`).

## CPP-243 — updateBaseline leaves the baseline indeterminate when no trace holds a peak

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; not executed, because the result is an indeterminate value.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.cpp:135–158`, `MassTraces::updateBaseline`; the constructor at 82–85 initialises only `max_trace`.

**Trigger:** One or more traces, none of which holds a peak.

**Issue:** The early return covers only an empty collection. With traces but no peaks the loop never assigns `baseline`, which the constructor never initialised, so it keeps an indeterminate value.

**Proposed C++ fix:** Initialise `baseline` to 0 in the constructor, or set it to 0 when no peak is seen.

**Evidence:** Source review; `tests/data/feature_finder_picked_helper_structs_provenance.json`.

**Rust handling:** The baseline starts at 0.0, and `update_baseline` leaves it unchanged when no trace holds a peak.

## CPP-244 — getCompensationVoltages lets a NaN into std::set<double>

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/IONMOBILITY/FAIMSHelper.cpp:45` (insert) and 50 (erase), `getCompensationVoltages`.

**Trigger:** A FAIMS spectrum whose drift time is NaN, from any reader that accepts `NaN` as a cvParam value.

**Issue:** NaN breaks the strict weak ordering `std::set<double>` requires. With the NaN on the first FAIMS spectrum, every later voltage compares equivalent to it and is dropped; `erase(DRIFTTIME_NOT_SET)` then removes the NaN, so the set comes back empty and a spurious missing-voltage warning is logged. With the NaN on a later spectrum, only the NaN is dropped.

**Proposed C++ fix:** Skip or reject NaN voltages when collecting.

**Evidence:** Oracle cases `nan_first`, `nan_middle` and `nan_last` of `../oracle/pte-faims-helper`; `tests/data/faims_helper_provenance.json`.

**Rust handling:** `FaimsHelper::get_compensation_voltages` returns `Error::InvalidValue` naming the spectrum index.

## CPP-245 — filterPeptidesByFAIMSCV silently accepts parameters that can never match

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/IONMOBILITY/FAIMSHelper.cpp:58–82`, `filterPeptidesByFAIMSCV` (comparison at 71).

**Trigger:** A `cv_tolerance` of zero, a negative or NaN tolerance, or a NaN `target_cv` at any tolerance, an infinite one included.

**Issue:** The strict test `std::abs(pep_cv - target_cv) < cv_tolerance` can then never succeed, so only unannotated identifications are returned, with no diagnostic. A NaN target keeps only the unannotated identifications even at an infinite tolerance, because every comparison with NaN is false. An infinite target also keeps no annotated identification, not even one annotated with the same infinity, because `inf - inf` is NaN. That follows from IEEE arithmetic on a voltage `getCompensationVoltages` can itself return, so the port treats it as source behaviour rather than as part of this defect.

**Proposed C++ fix:** Validate `target_cv` and `cv_tolerance`.

**Evidence:** Oracle cases `tolerance_zero`, `tolerance_negative`, `tolerance_nan`, `target_nan` and `target_nan_tolerance_infinite`; for infinite targets `target_positive_infinity`, `target_negative_infinity`, their `_tolerance_infinite` variants and the two `faims_filter_infinite_annotation` cases. All are in `../oracle/pte-faims-helper` (fix-round manifest `c374c328`); `tests/data/faims_helper_provenance.json`.

**Rust handling:** `FaimsHelper::filter_peptides_by_faims_cv` returns `Error::InvalidValue` for a NaN target or a NaN, zero or negative tolerance; a positive infinite tolerance stays valid. Since `f3c29cb` an infinite target is filtered as the source filters it and returns the C++ identifiers. `MetaValue` cannot hold an infinite annotation, so the infinite-annotation cases run on their representable subset.

## CPP-246 — PeakTypeEstimator's comment states an 80% threshold for a 75% test

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed documentation defect; behaviour is unaffected.

**Affected file/function:** `src/openms/include/OpenMS/FORMAT/PeakTypeEstimator.h:147`, `estimateType`.

**Issue:** The line reads `if (evidence_ratio > 0.75) // 80% are profile`.

**Proposed C++ fix:** Correct the comment to 75%.

**Evidence:** Source review; `docs/PEAK_TYPE_ESTIMATOR_SUPPORT.md`.

**Rust handling:** The port uses 0.75, and its documentation states 0.75.

## CPP-247 — FeatureFinderAlgorithmPicked's abundance override keeps a stray (0, 1) peak

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:163–179`, `run`; `src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.cpp:33–36`, the default constructor.

**Trigger:** Any non-default `isotopic_pattern:abundance_12C` or `isotopic_pattern:abundance_14N`.

**Issue:** The override inserts the two isotopes into a default-constructed `IsotopeDistribution`, whose constructor already holds `(0, 1)`. The override becomes `{(0, 1), (12, a), (13, 1 - a)}`, which puts most of the weight 12 Da below carbon-12. Executed at 12C = 90% with `max_isotopes` 1020: windows 0, 1, 5 and 10 have 30, 110, 436 and 811 bins instead of the intended 6, 26, 148 and 247.

**Proposed C++ fix:** `set()` the two-isotope container, or `clear()` before inserting.

**Evidence:** `../oracle/b2-iso-source-precision/manifest.json` (`probe.tsv`); `tests/data/isotopes_source_precision_provenance.json`; the sizes are asserted in `tests/isotopes_source_precision.rs`.

**Rust handling:** `CoarseIsotopePatternGenerator::set_isotope_override` rejects that construction. FeatureFinderAlgorithmPicked's seed stage (B6-FFAP-SEEDS, `80bbdf1`) refuses a changed abundance with `Error::Unsupported` by default and builds the intended two-isotope distribution only under `AbundanceOverride::Intended`; C2 shows the source's effect on FeatureFinderCentroided_1 (27-bin windows and 0 seeds at 12C = 90%). Since wave 5 the port computes the intended two-isotope override by default (`AbundanceOverride::Intended`, lead decision of 2026-09-15), which FeatureFinderCentroided uses; `AbundanceOverride::Refuse` is the opt-out. The intended result is pinned against an adapted Release replay (`intended_abundance.cpp`: FFC_1 with 12C 90 finds 18 seeds, 1 candidate and 1 feature; 12C 99 25/8/8; 14N 95 13/2/2), where the executed C++ finds 0/0/0 at 12C 90.

## CPP-248 — CoarseIsotopePatternGenerator::run gives different bits in different runs of one binary

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed 200 times on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h:66` (`MapType_`) and 341 (`formula_`); `src/openms/source/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.cpp:114–120`, `run`; `src/openms/source/CHEMISTRY/EmpiricalFormula.cpp:57–67`, `getLightestIsotopeWeight`; element allocation at `src/openms/source/CHEMISTRY/ElementDB.cpp:580–588` and 634–664.

**Trigger:** Repeated executions with the same formula.

**Issue:** `formula_` is a `std::map<const Element*, SignedSize>`, ordered by heap address. `run` convolves and `getLightestIsotopeWeight` sums in that order, and binary32 accumulation depends on it. `CoarseIsotopePatternGenerator(0).run(EmpiricalFormula("C1H1N1O1S1P1"))` iterated `H C N O P S` in 198 of 200 runs, with bin-2 intensity `0x3d3540d4`, and `H N C O P S` in runs 40 and 147, printing `0x3d3540d5`. The same two runs changed `estimateFromPeptideWeight` for every FeatureFinderAlgorithmPicked window from 150 to 8050 Da. Br, Na, He, B and labelled isotopes also move between runs.

**Proposed C++ fix:** Order the map by atomic number and isotope mass number instead of by pointer.

**Evidence:** `../oracle/b2-iso-element-order/manifest.json` (`probe.cpp`, `run.sh`, `tally.py`, `results/runs.sha256`); `tests/data/isotopes_source_precision/element_order.tsv`.

**Rust handling:** The port iterates in ascending atomic number, with labelled isotopes after their natural element, so `ProbabilityPrecision::SourceSingle` reproduces the SDK runs that use that order, the majority for natural elements.

## CPP-249 — ElementDB builds iridium from rhenium's tables

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/CHEMISTRY/ElementDB.cpp:512`, the element table construction; `iridium_abundance` and `iridium_mass` at 510–511 are unused.

**Issue:** The call is `buildElement_("Iridium", "Ir", 77u, rhenium_abundance, rhenium_mass)`, so every iridium mass and isotope pattern uses rhenium's isotopes. `Os3Ir3` has a lightest-isotope weight of 1106.716344 Da (three osmium-184 plus three rhenium-185) in every one of 200 runs, where iridium gives 1124.739252 Da.

**Proposed C++ fix:** Pass `iridium_abundance` and `iridium_mass`.

**Evidence:** Case `pair_Os3Ir3` of `../oracle/b2-iso-element-order/manifest.json`; `tests/isotopes_source_precision.rs`.

**Rust handling:** The port uses the declared iridium table (`docs/CHEMISTRY_SUPPORT.md`); a test asserts that its `Os3Ir3` bits differ from every SDK run.

## CPP-250 — -instance cannot be used

**Source revision:** cli `c19e49414bcd9ebdea42f89b3f74d2823205892c` (`source/APPLICATIONS/TOPPBase.cpp`) with core `bc9cc12514c768385ce121d6ca4bb710fe1983c4` (`Param.cpp`). Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `source/APPLICATIONS/TOPPBase.cpp:166` (registration), 2104 (`getDefaultParameters_` exclusion) and 339 (the strict `Param::update`).

**Trigger:** Any tool run with `-instance <n>`.

**Issue:** `instance` is registered, but `getDefaultParameters_` excludes it while the command-line `Param` keeps it. `Param::update` with `fail_on_unknown_parameters` set then rejects every run that passes `-instance` ("Unknown (or deprecated) Parameter 'instance' given in outdated parameter file!", exit 6). INI instance sections other than 1 are unreachable, and TOPPBase_test's instance 5 and 6 `getStringOption_` checks pass only because they read `param_` after the failed update.

**Proposed C++ fix:** Remove `instance` from the command-line parameters before the update, as the lifecycle already does for `ini` (333).

**Evidence:** Oracle case `instance_on_command_line` in `../oracle/topp-cli-lifecycle/manifest.json`; `tests/data/topp_cli_lifecycle_provenance.json`.

**Rust handling:** The port reproduces exit 6. `docs/TOPP_CLI_SUPPORT.md` records why the `-instance 5` class-test case is not transcribed.

## CPP-251 — common: values override instance values for subsection parameters

**Source revision:** cli `c19e49414bcd9ebdea42f89b3f74d2823205892c` with core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `source/APPLICATIONS/TOPPBase.cpp:305` (`param_common_ = param_inifile_.copy("common:", true)`) and 320–330 (merge order); `src/openms/source/DATASTRUCTURES/Param.cpp:1262–1281`, leaf-name matching in `Param::update`.

**Trigger:** An INI with both `<Tool>:1:section:name` and `common:<Tool>:section:name`.

**Issue:** `param_common_` keeps the nested key `<Tool>:section:name`, and `finalParam.merge` adds it next to the instance value `section:name`. `Param::update` finds the nested key by its leaf name after it has applied the instance value, so the common value wins, the reverse of the intended precedence.

**Proposed C++ fix:** Take tool-specific common values only from `common:<Tool>:` (already copied separately) and strip that prefix before merging, so every value is matched by its full name.

**Evidence:** Oracle case `ini_instance_and_common` in `../oracle/topp-cli-lifecycle/manifest.json`: the output has `peakcount=1` although the instance section says 2.

**Rust handling:** The port reproduces the source precedence.

## CPP-252 — A common:<tool>: value for a top-level parameter is rejected

**Source revision:** cli `c19e49414bcd9ebdea42f89b3f74d2823205892c` with core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `source/APPLICATIONS/TOPPBase.cpp:305` and 339; `src/openms/source/DATASTRUCTURES/Param.cpp:1171–1183` (`findFirst`) and 1262–1281.

**Trigger:** An INI with `common:<Tool>:threads`, a top-level parameter.

**Issue:** The nested copy `<Tool>:threads` has no exact match, and `findFirst` matches only names that end in `:threads`, which excludes the root-level `threads`. The strict update therefore fails with exit 6 ("Unknown (or deprecated) Parameter 'SpectraFilterWindowMower:threads' …").

**Proposed C++ fix:** The same as CPP-251: strip the `<Tool>:` prefix from tool-specific common values before merging.

**Evidence:** Oracle case `ini_common_top_level` in `../oracle/topp-cli-lifecycle/manifest.json`.

**Rust handling:** The port reproduces exit 6.

## CPP-253 — StringUtils::number silently cuts its text at 63 bytes

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`, and by a probe that compiles a verbatim copy of the pinned function.

**Status:** Executed.

**Affected file/function:** `src/openms/source/DATASTRUCTURES/StringUtils.cpp:526–531`, `StringUtils::number(double, UInt)`.

**Trigger:** A value whose `%.*f` text is longer than 63 bytes, such as a corrupt or huge intensity or retention time of about 1e60 or more at the digit counts FileInfo uses.

**Issue:** `number` formats into `char buf[64]` with `std::snprintf(buf, sizeof(buf), "%.*f", static_cast<int>(n), d)` and returns the buffer without checking the length `snprintf` reports. The text is cut to 63 bytes with no diagnostic, so the magnitude is wrong and the decimals are gone: `number(1e100, 0)`, `number(1e100, 1)` and `number(1e100, 2)` all return the 63-digit integer `100000000000000001590289110975991804683608085639452813897813275`. FileInfo prints every retention-time, m/z, ion-mobility and intensity range through `number`, in the text report and in the `general:` TSV lines (`FileInfo.cpp:108–540`).

**Proposed C++ fix:** Size the buffer from `snprintf(nullptr, 0, ...)`, or format with `std::to_chars` into a buffer large enough for the value.

**Evidence:** The `number(x, 0..2)` columns of the D rows for bits `54b249ad2594c37d` and `d4b249ad2594c37d` in `../oracle/file-info-text-format/results/driver.tsv`; the same rows in `results/pin_probe.tsv`, from a probe that compiles a verbatim copy of `StringUtils.cpp:526–531`; `tests/data/file_info_text_format_provenance.json`.

**Rust handling:** `format::file_info::text_format::fixed` refuses a value whose text the source would cut; `fixed_truncated` reproduces the source bytes for the tool path.

## CPP-254 — StringUtils::number turns a digit count of 2^31 or more into a negative precision

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`, and by a probe that compiles a verbatim copy of the pinned function.

**Status:** Executed.

**Affected file/function:** `src/openms/source/DATASTRUCTURES/StringUtils.cpp:529`, the `static_cast<int>(n)` in `StringUtils::number(double, UInt)`.

**Trigger:** A digit count `n` of 2147483648 or more.

**Issue:** The `UInt` digit count becomes the `int` precision of `%.*f`. From 2^31 on it is negative, and a negative precision argument counts as omitted, so the value prints with six decimals instead of the requested count, without a diagnostic: `number(0.125, 2147483648)` returns `0.125000`.

**Proposed C++ fix:** Range-check `n` before the cast.

**Evidence:** The F rows with digit count 2147483648 in `../oracle/file-info-text-format/results/driver.tsv` (`0.125000`, `2.500000`, `-0.000000`, `0.333333`, `123.456000`) and in `results/pin_probe.tsv`; `tests/data/file_info_text_format_provenance.json`.

**Rust handling:** `fixed` refuses such a digit count; `fixed_truncated` reproduces the source text.

## CPP-255 — StringUtils::toStr documents significant digits but writes fraction digits

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The behaviour is executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Documentation defect. The behaviour is executed; whether the comments or the code state the intended precision is unconfirmed.

**Affected file/function:** `src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h:118–121`, the `toStr(float, bool)` and `toStr(double, bool)` comments; `src/openms/source/DATASTRUCTURES/StringUtils.cpp:313–316`, the precision-mapping comment, with the code at 374–377 and 384–387; `src/common/include/OpenMS/CONCEPT/Detail/NumericFormatting.h:42` and 82.

**Issue:** The header says `full_precision` selects "6-digit" (float) or "15-digit" (double) output, and the implementation comment says 6 and 15 "significant digits (general)". `appendToStr` passes `writtenDigits<float>()` or `writtenDigits<double>()` to `appendNumeric`, which for magnitudes in [1e-2, 1e4) calls `std::to_chars` with `std::chars_format::fixed` and that precision, a count of fraction digits. `toStr(1234.5678)` is `1234.567800000000034`, 19 significant digits including binary noise, and `toStr(9.995)` is `9.994999999999999`. At 1e4 and above, and below 1e-2, the same function writes the shortest round-trip scientific text, so the precision differs across magnitudes.

**Proposed C++ fix:** Correct the comments, or use general or shortest formatting in the fixed range; the latter changes output and needs a review of the reference files.

**Evidence:** The `toStr(x)` column of the D rows for bits `40934a456d5cfaad` and `4023fd70a3d70a3d` in `../oracle/file-info-text-format/results/driver.tsv` and `results/pin_probe.tsv`; source review of the cited lines.

**Rust handling:** `to_str` and `to_str_f32` reproduce the source behaviour, not the comments.

## CPP-256 — SignalToNoiseEstimatorMedian's AUTOMAXBYPERCENT mode reads and writes out of bounds

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The crash is executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85` (C1).

**Status:** Executed (crash) on the product SDK; the out-of-domain behaviour is executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576` (`../oracle/sne-completion/probes`: 7 of 9 inputs SIGSEGV/SIGABRT, 2 read neighbouring memory); the defined domain is derived in `docs/SIGNAL_TO_NOISE_SUPPORT.md`.

**Affected file/function:** `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h:191–233`, `computeSTN_`, the `AUTOMAXBYPERCENT` branch.

**Trigger:** `auto_mode = 1` with estimation running, for example PeakPickerHiRes with `signal_to_noise > 0`.

**Issue:** `std::max_element` is called with the comparator `a.getIntensity() > b.getIntensity()` (line 208), so it returns the minimum intensity, not the maximum. `bin_size = maxInt / 100` (211) is then 0 for a spectrum with a zero intensity, and `++histogram_auto[(int)((peak.getIntensity() - 1) / bin_size)]` (216) indexes the 100-bin vector with no bounds check: a quotient far above 99, a negative index for intensities below 1, or a division by zero converted to `int`. An empty container dereferences `end()` (209). Beyond `INT_MAX` points the `int` counters at `:216` and `:228` can overflow, and `(int)(p * n / 100)` at `:220` is undefined from `2^31` on; the Release build's 32-bit `cvttsd2si` returns `INT_MIN` there and skips the walk (measured at `n = 2^31` and `n = 3,000,000,001`, `../oracle/sne-fix`).

**Proposed C++ fix:** Use `std::max_element` with `<` (or `getIntensity()` less), return early for an empty container, guard `bin_size > 0`, and clamp the bin index to `[0, 99]`.

**Evidence:** C1 (`../oracle/topp-early-bundle`) records PeakPickerHiRes with `auto_mode 1` ending in SIGBUS 138 or SIGSEGV 139, varying between attempts. The P1 oracle case `extra_auto_mode_percentile_sn0` (`../oracle/peak-picker-hires`) shows the mode is accepted when no estimation runs; `tests/data/peak_picking_provenance.json`.

**Rust handling:** Computed exactly on the defined domain (no empty container, every quotient in `(-1, 100)`, `int` counters within range) in both profiles; the `:220` `INT_MIN` is reproduced; `Error::Unsupported` names `:209`, `:216`, `:228` or `:365` elsewhere (`docs/SIGNAL_TO_NOISE_SUPPORT.md`).

## CPP-257 — SignalToNoiseEstimatorMedian converts an unbounded bin quotient to int

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`: 8 `cpp257_*` cases (`../oracle/sne-completion`), including `PeakPickerHiRes::pick`; the 32-bit `cvttsd2si` at `libOpenMS.so` `0x186c5d0`/`0x186c64d` gives `INT_MIN`, which the clamp sends to bin 0.

**Affected file/function:** `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h:297` and `:308`, `computeSTN_`.

**Trigger:** A manual `max_intensity` small against the data, so that `intensity / bin_size` exceeds `INT_MAX`; for example `max_intensity 1` with intensities above 2^31; also `auto_max_stdev_factor 0` with `bin_count 1,000,000` and one intensity at `f32(1e12)` among 10,000, and negative intensities, which `PeakPickerHiRes` accepts.

**Issue:** `(int)(intensity / bin_size)` is converted before `std::min<int>` clamps it to the last bin. A `double` outside the `int` range makes the conversion undefined. arm64 saturates to `INT_MAX` (the last bin); x86-64 `cvttsd2si` gives `INT_MIN`, which `std::max(..., 0)` clamps to bin 0, so the same input can land in the first or the last histogram bin depending on the platform. The automatic modes reach it too: `auto_max_stdev_factor = 0` with a large `bin_count` (`cpp257_stdev_bins`), and negative intensities, which pull the range down against the largest intensity (`cpp257_neg_factor0`, `cpp257_neg_default`).

**Proposed C++ fix:** Clamp in `double` before converting: `std::min(intensity / bin_size, double(bin_count_minus_1))`, then cast.

**Evidence:** Source review of lines 258, 297 and 308; the P1 verifier's adversarial case against `../oracle/peak-picker-hires` (arm64). Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`: the 8 `cpp257_*` cases of `../oracle/sne-completion`, twice each, including `PeakPickerHiRes::pick` (`cpp257_pick_manual`).

**Rust handling:** `BinIndexConversion::X86_64Release` (`PickingCompatibility::source()`, the TOPP tools) reproduces bin 0; the native profile's `ClampBeforeTruncation` clamps first (the last bin, arm64's and the parameter description's answer).

## CPP-258 — PeakPickerHiRes's FWHM bisection can loop forever

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed. A hang cannot be recorded by the oracle; the port detects the fixed point.

**Affected file/function:** `src/openms/source/PROCESSING/CENTROIDING/PeakPickerHiRes.cpp:396–408` (left) and `:421–434` (right), `pick_`, with `report_FWHM` set.

**Trigger:** A picked peak whose spline maximum is not positive, or any search whose midpoint reaches a bracket end without meeting the tolerance.

**Issue:** The `do`/`while (fabs(int_mid - fwhm_int) > threshold)` loops have no step limit. With `max_peak_int <= 0`, `threshold = 0.01 * fwhm_int` is not positive, so the condition never becomes false. Otherwise, once `mz_left` (or `mz_right`) and `mz_center` are adjacent doubles, the midpoint equals one of them and the bracket stops moving; if the spline there is not within the tolerance of half height, the loop never ends.

**Proposed C++ fix:** Skip FWHM for a non-positive maximum and stop when the midpoint equals a bracket end or after a fixed number of halvings.

**Evidence:** Source review; `docs/PEAK_PICKING_SUPPORT.md`, native difference 4.

**Rust handling:** The port detects the fixed point exactly and returns `Error::InvalidValue`, with a 4,096-halving guard that no terminating search reaches.

## CPP-259 — PeakPickerHiRes weights ion mobility with samples its spline discards

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/PROCESSING/CENTROIDING/PeakPickerHiRes.cpp:250–260`, `:290–299`, `:333–342` and `:444–446`, `pick_`.

**Trigger:** A profile spectrum with an ion-mobility array and two samples at the same m/z inside one peak.

**Issue:** The support is a `std::map<double, double>`, so `peak_raw_data[pos] = intensity` overwrites an equal key, while `weighted_im += im * intensity` adds every sample. The weighted mean `weighted_im / total_intensity` then divides a sum over all samples by the intensity total of the samples the map kept, and the reported mobility is no longer a weighted mean of any subset.

**Proposed C++ fix:** Accumulate the mobility weight in the map entry (or rebuild it from the final map), or reject duplicate positions before picking.

**Evidence:** Oracle cases `source_duplicate_apex`, `source_duplicate_apex_flank`, `source_duplicate_apex_nocheck` and `source_duplicate_extension` on inputs with an `Ion Mobility` array (`../oracle/peak-picker-hires`; `tests/data/peak_picking/synthetic.tsv`).

**Rust handling:** The default refuses duplicate positions; `PickingCompatibility::source()` reproduces the source result bit for bit.

## CPP-260 — Orphaned PeakPickerHiRes class-test fixtures no longer match the code

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Test-data defect; executed.

**Affected file/function:** `src/tests/class_tests/openms/data/PeakPickerHiRes_orbitrap_sn0_out.mzML`, `PeakPickerHiRes_ftms_sn0_out.mzML`, `PeakPickerHiRes_orbitrap_ppmax.mzML`, `PeakPickerHiRes_ftms_ppmax.mzML`, `PeakPickerHiRes_orbitrap_sn4_out_ppmax.mzML` and `PeakPickerHiRes_ftms_sn4_out_ppmax.mzML`.

**Trigger:** Using these files as expected outputs.

**Issue:** No class test references the six files. The current picker does not reproduce the `sn0` outputs: at `signal_to_noise 0` it gives 82, 112 and 89 (orbitrap) and 314 and 319 (FTMS) centroids per spectrum against the stored 679, 860 and 640 and 9,359 and 9,384. `PeakPickerHiRes_orbitrap_ppmax.mzML` holds 9,778 points against a 1,210-point input, so it is not an output of that input.

**Proposed C++ fix:** Remove the files, or regenerate them and add the class-test sections that use them.

**Evidence:** P1's runs of the product SDK on the class-test inputs (`../oracle/peak-picker-hires`); a search of `src/tests` and `src/openms` at the pin finds no reference.

**Rust handling:** Not used; `docs/PEAK_PICKING_SUPPORT.md` records them as unused.

## CPP-261 — PeakPickerHiRes FTMS class-test files repeat a spectrum id

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Test-data defect; source-reviewed.

**Affected file/function:** `src/tests/class_tests/openms/data/PeakPickerHiRes_ftms.mzML`, `PeakPickerHiRes_ftms_sn1_out.mzML`, `PeakPickerHiRes_ftms_sn4_out.mzML`, `PeakPickerHiRes_ftms_sn0_out.mzML`, `PeakPickerHiRes_ftms_ppmax.mzML` and `PeakPickerHiRes_ftms_sn4_out_ppmax.mzML`.

**Trigger:** Loading the files with a reader that requires unique spectrum ids.

**Issue:** Each file has two `<spectrum>` elements with `id="spectrum=1"`. mzML requires a spectrum id to be unique within the run; `MzMLHandler` loads the files anyway, which hides the duplicate from the class test.

**Proposed C++ fix:** Renumber the second spectrum in the input and in the outputs.

**Evidence:** `grep` of the pinned files; `PeakPickerHiRes_test.cpp:281–350` loads them.

**Rust handling:** The native mzML reader refuses duplicate record ids, so P1 commits copies with only the second id renamed (`*.unique_ids.mzML`, recorded as adapted in `tests/data/peak_picking_provenance.json`).

## CPP-262 — Dangling mzML software and data-processing references are silently default-constructed

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp:920`, `:924`, `:948`, `:952`, `:1034`, `:1264` and `:1288`, `startElement`.

**Trigger:** An mzML `softwareRef`, `processingMethod/@softwareRef`, `dataProcessingRef` or `defaultDataProcessingRef` naming no definition, for example the upstream TOPP fixture `PeakPickerHiRes_5_input.mzML` (softwareRef `so_in_0` with no `softwareList`, `defaultDataProcessingRef` `dp_sp_0` with no `dataProcessingList`).

**Issue:** The references are resolved with `std::map::operator[]`, which inserts a default entry, so the instrument or method silently gets an empty `Software` and the record an empty processing history. No diagnostic is written, while an unregistered spectrum `sourceFileRef` is checked with `contains()` and warned about (`:899–906`). A software list placed after the instrument list also resolves to empty software.

**Proposed C++ fix:** Look the IDs up with `find`; warn or throw `ParseError` for an unknown ID, as for `sourceFileRef`.

**Evidence:** `../oracle/p2-mzml-leniency` (load, metadata-only load and transform on the upstream fixture and four synthetic cases; exit 0, stderr without a reference diagnostic); `tests/data/mzml_header_leniency_provenance.json`.

**Rust handling:** The reader rejects these references by default; `mzml::ReadOptions::source_dangling_references` reproduces the source result and warns once per distinct ID.

## CPP-263 — MzMLSplitter's parts carry no processing record

**Source revision:** topp `174b576e244e100f2345ca57a8e79aaa607156df` (`src/MzMLSplitter.cpp`) with cli `c19e49414bcd9ebdea42f89b3f74d2823205892c` (`TOPPBase.cpp`). Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/MzMLSplitter.cpp:146–167`, `main_`; `source/APPLICATIONS/TOPPBase.cpp:594–605`, `addDataProcessing_(PeakMap&, ...)`.

**Trigger:** Any MzMLSplitter run.

**Issue:** `addDataProcessing_(part, getProcessingInfo_(DataProcessing::FILTERING))` runs on a copy of the experiment whose spectra and chromatograms were moved out, before the part's own spectra and chromatograms are added. `addDataProcessing_` attaches the record to each spectrum and chromatogram present, so no output part records the filtering step.

**Proposed C++ fix:** Call `addDataProcessing_` after the spectra and chromatograms are added to the part.

**Evidence:** Oracle case `mzml_splitter_1` in `../oracle/topp-cli-lifecycle/cli2/manifest.json`; the retained `MzMLSplitter_output_part1/2.mzML` (test-data `0cb15f2`).

**Rust handling:** The port reproduces the call order, so its parts carry no record either (`tests/topp_mzml_splitter.rs`).

## CPP-264 — ParamXMLFile declares ISO-8859-1 but writes UTF-8 bytes

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed.

**Affected file/function:** `src/openms/source/FORMAT/ParamXMLFile.cpp:66`, `writeXMLToStream`.

**Trigger:** A parameter name, value, description, tag or restriction containing a non-ASCII character, written with `-write_ini` or `ParamXMLFile::store`.

**Issue:** The writer emits `<?xml version="1.0" encoding="ISO-8859-1"?>` and then copies the UTF-8 bytes of its `std::string` values, escaping only XML markup. A reader honouring the declaration decodes each multi-byte character as several Latin-1 characters, so the text does not read back as written, in Xerces as in any other parser.

**Proposed C++ fix:** Declare `UTF-8`, or transcode to ISO-8859-1 and write characters above U+00FF as character references.

**Evidence:** Source review of the declaration and the escaping path (`XMLHandler::writeXMLEscape`); the upstream writer golden `ParamXMLFile_test_writeXMLToStream.xml` is ASCII, so it does not expose the mismatch.

**Rust handling:** `paramxml::WriteOptions::source()` writes the ISO-8859-1 declaration with bytes consistent with it: characters up to U+00FF as their Latin-1 byte, others as character references (`docs/PARAMXML_SUPPORT.md`).

## CPP-265 — XML files are opened twice, so a FIFO given as input blocks

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FORMAT/XMLFile.cpp:141–166`, `parse_`.

**Trigger:** An XML input read from a FIFO fed by a single writer, for example `-ini` or `-write_ini -ini` with a named pipe.

**Issue:** `parse_` opens the file with `std::ifstream` to peek at two bytes for the compression check, closes it, and opens it again through `LocalFileInputSource`. The first open consumes the writer's data; the second open waits for a writer that never comes, and the tool blocks without a message.

**Proposed C++ fix:** Peek and parse through one stream, for example by wrapping the already opened stream in a Xerces `InputSource`.

**Evidence:** Oracle observation `write_ini_ini_fifo_single_writer` in `../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json`: killed by the 10-second alarm, exit 142, nothing written.

**Rust handling:** The port reads the file once and completes; `docs/TOPP_CLI_SUPPORT.md` records this as a deliberate difference.

## CPP-266 — An INI that exists but cannot be opened is reported as an internal error

**Source revision:** cli `c19e49414bcd9ebdea42f89b3f74d2823205892c` with core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FORMAT/TextFile.cpp:36–47`, `load`; cli `source/APPLICATIONS/TOPPBase.cpp:495–499`.

**Trigger:** `-ini` naming a Unix socket, or `/dev/tty` in a process without a controlling terminal, on a run or with `-write_ini`.

**Issue:** The file exists and `File::readable` holds, but the open fails, so `TextFile::load` throws `IOException`. TOPPBase has no handler for it and reaches the `BaseException` arm: "Error: Unexpected internal error (IO error for file ...)", exit 8 (`UNKNOWN_ERROR`), where an unreadable input is exit 2 and a corrupt one exit 3.

**Proposed C++ fix:** Map `IOException` on an input file to `INPUT_FILE_NOT_READABLE` with a message naming the open failure.

**Evidence:** Oracle cases `ini_socket`, `write_ini_ini_socket`, `ini_tty` and `write_ini_ini_tty` in `../oracle/topp-cli-lifecycle/ini_read_failures/manifest.json`.

**Rust handling:** The port reproduces exit 8 and the message for an open failure other than NotFound and PermissionDenied; a read failure after a successful open stays exit 3.

## CPP-267 — FileInfo::Result declares fields that run never fills

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed.

**Affected file/function:** `src/openms/include/OpenMS/FORMAT/FileInfo.h:194–198` (`ValidationInfo::schema_version`, `detail`) and `:220–236` (`Result::experiment_meta`, `statistics`, `corruption`, `detail`); `src/openms/source/FORMAT/FileInfo.cpp`, `run` and `report_`.

**Trigger:** Any library or pyOpenMS caller reading the structured result after `-m`, `-s`, `-c`, `-d` or `-v`.

**Issue:** `FileInfo.cpp` never assigns `experiment_meta`, `statistics` (the `NamedStats` blocks), `corruption`, `detail`, `ValidationInfo::schema_version` or `ValidationInfo::detail`. The metadata, statistics, corruption and detail output goes only into the text report, so structured callers see empty records that look like results.

**Proposed C++ fix:** Fill the fields where the text is produced, or remove them from the public result.

**Evidence:** Source review: no assignment to these members exists in `FileInfo.cpp` at the pin.

**Rust handling:** `model::FileInfoResult` keeps the fields and leaves them empty, as the source does (`docs/FILE_INFO_SUPPORT.md`).

## CPP-268 — FileInfo computes the FAIMS compensation voltages twice per peak file

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed.

**Affected file/function:** `src/openms/source/FORMAT/FileInfo.cpp:1673` and `:1742`, `report_`.

**Trigger:** Any peak file, and in particular one with a FAIMS spectrum whose voltage is missing.

**Issue:** `FAIMSHelper::getCompensationVoltages(exp)` runs once to fill `PeakInfo::faims_cvs` and again to print the `IM (FAIMS_CV)` line. Each call scans every spectrum and, for a missing voltage, logs "FAIMS compensation voltage is missing for at least one spectrum!", so the scan runs twice and the warning appears twice.

**Proposed C++ fix:** Reuse `pk.faims_cvs` for the text line.

**Evidence:** Source review of both call sites.

**Rust handling:** The port scans once and records the warning once in `FileInfoResult::warnings`.

## CPP-269 — FileInfo_test copies to a fixed temporary file name

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Test defect; source-reviewed.

**Affected file/function:** `src/tests/class_tests/openms/source/FileInfo_test.cpp:133–143`, the forced-type section.

**Trigger:** Two concurrent runs of the class test sharing a temporary directory.

**Issue:** The section copies its input to `File::getTempDirectory() + "/test_forced_type.tmp"` and removes it at the end, so concurrent runs overwrite and delete each other's file.

**Proposed C++ fix:** Use `File::TempDir` or a unique name.

**Evidence:** Source review.

**Rust handling:** `class_test_forced_type_selects_the_parse_branch` uses its own `TempDir`.

## CPP-270 — Retained FileInfo outputs keep a stale intensity padding

**Source revision:** test-data `0cb15f23fccc6ea196bfafcfbbf020958f36c3a3` with core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Test-data defect; source-reviewed.

**Affected file/function:** test-data `topp/FileInfo_3_output.txt:11` and `topp/FileInfo_7_output.txt:14`; the writer at `src/openms/source/FORMAT/FileInfo.cpp:140` and `:189`.

**Trigger:** Comparing FileInfo's `intensity:` range line with the retained outputs.

**Issue:** The retained files have six spaces after `intensity:`; the source at the pin writes one. TOPP_FileInfo_3 and _7 pass only because FuzzyDiff ignores whitespace differences, so the expectations no longer describe the output.

**Proposed C++ fix:** Regenerate the two retained outputs.

**Evidence:** The two retained lines against the two source lines.

**Rust handling:** The port writes one space and compares the retained FileInfo_3 output under FuzzyDiff, accepting both spellings.

## CPP-271 — mass_trace:min_spectra = 1 silently finds nothing

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:48–49` (default and minimum 1), `:1112` (`updateMembers_`) and `:340` (`run`).

**Trigger:** `mass_trace:min_spectra` set to 1, its declared minimum.

**Issue:** `min_spectra_ = floor(1 * 0.5) = 0`, so every trace score is `0 / 0 = NaN` and every peak a local maximum. No overall score reaches a threshold: the run reports 0 seeds and 0 features, exits 0 and warns about nothing.

**Proposed C++ fix:** Set the minimum to 2, or reject values that give `min_spectra_ == 0`.

**Evidence:** B6 driver case `ffc1_min_spectra_1` in `../oracle/b6-ffap-seeds/manifest.json`; `tests/data/feature_finder_picked_provenance.json`.

**Rust handling:** `Error::InvalidValue`; whether the port follows the source is open for the lead (`docs/FEATURE_FINDER_PICKED_SUPPORT.md`).

## CPP-272 — The overall seed score depends on the platform's powf

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85` (macOS arm64). Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576` as well.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:506`, `run`.

**Trigger:** Any run; the score feeds `seed:min_score`.

**Issue:** `std::pow(float, float)` calls the C library `powf`, which is not correctly rounded on every platform. On macOS 99 of 30,840 executed overall scores are one binary32 step from the correctly rounded value (12 of 3,084 on FeatureFinderCentroided_1), checked against 60-digit decimal arithmetic. A score next to `seed:min_score` can therefore select or drop a seed depending on the platform. The Linux x86_64 Release build binds `powf@GLIBC_2.27` of glibc 2.39 (Ubuntu `2.39-0ubuntu8.9`), whose ifunc selects `__powf_fma` on the AMD EPYC 7763 reference node; it is one binary32 step from the correctly rounded value on 8 of the 30,840 scores (2 of 3,084 on FeatureFinderCentroided_1), other scores than Apple's.

**Proposed C++ fix:** Evaluate the cube root in `double` and round once to `float`.

**Evidence:** C2 `ffap_stages` and `../oracle/b6-ffap-seeds` score arrays on both builds; `tests/data/feature_finder_picked/overall_rounding.tsv` lists the 8 Linux Release differences; `../oracle/ffap-complete-fix1` (`powf_probe` through `dlsym`: the resolved function, a 42x42 special grid, four generated sets of `2^26` pairs, every binary32 base with the exponent `1.0f/3.0f`; `logs/powf_fma_disasm.txt`).

**Rust handling:** `overall_score` computes the reference build's `powf`: `glibc_powf.rs` ports `__powf_fma` from Arm optimized-routines (MIT), with the FMA fusion and special cases read from the glibc 2.39 disassembly, and equals the executed library on every probed pair (all `2^32` bases with the exponent `1/3` included), so every overall score, the 8 misrounded ones included, is the Release build's on every host (lead decision D5 of wave 5).

## CPP-273 — write_debug reads an undeclared parameter and throws

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Executed (Linux x86_64 Release, `bc9cc12-c19e494-174b576`).

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:2137`, `writeFeatureDebugInfo_`; the declaration at `:124`.

**Trigger:** `debug` output enabled (the TOPP parameter `write_debug`), once a feature is written.

**Issue:** `writeFeatureDebugInfo_` reads `param_.getValue("debug:pseudo_rt_shift")`, but the declared parameter is `advanced:pseudo_rt_shift`. The lookup throws `ElementNotFound`, which escapes the OpenMP region that calls it.

**Proposed C++ fix:** Read `advanced:pseudo_rt_shift`.

**Evidence:** Source review of both lines. `FeatureFinderCentroided -algorithm:write_debug` on FeatureFinderCentroided_1: SIGABRT, shell status 134, in 5 of 5 single-thread runs and 3 of 3 four-thread runs (`../oracle/ffap-instr-completion` tool cases `a4`, `b1`, `c1`). An INI that supplies `algorithm:debug:pseudo_rt_shift` is dropped by the tool (verifier, SIGABRT).

**Rust handling:** Ported: `PseudoRtShiftKey::Source` stops with `Error::Unsupported` at that seed, and the tool exits 8 after writing what the executed run wrote. `PseudoRtShiftKey::Declared` reads `advanced:pseudo_rt_shift`. A numeric shift, infinite and NaN included, writes the executed files byte for byte (glibc's `-nan` in the `.plot` formulas).

## CPP-274 — A single retention time or m/z makes the intensity score undefined

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK (Debug) at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85` and on the C++ Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed (2026-09-15; source-reviewed before that).

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:244–245` (bin widths) and `:1837–1838`, `intensityScore_`.

**Trigger:** An MS1 input whose spectra share one retention time, or whose peaks share one m/z.

**Issue:** `intensity_rt_step_` or `intensity_mz_step_` is 0, so `(rt - rt_min) / intensity_rt_step_` is `0 / 0 = NaN`, and `(UInt) std::floor(NaN)` is undefined behaviour before `std::min` clamps it.

**Proposed C++ fix:** Use one bin per dimension when the range is empty, or reject such input with a message.

**Evidence:** Executed on the Release build `openms4-release-bc9cc12-c19e494-174b576` (port/ffap-semantics, `../oracle/ffap-sem-completion`): 26 configurations of the driver `degenerate_stage` (every RT equal, every m/z equal, a subnormal and an overflowing RT extent, with the FFC_1 INI, the defaults and `seed:min_score 0`; three repetitions at 1 and 4 threads, identical) and the probe `iscore_probe` (57 positions on four grids). `libOpenMS.so` compiles the conversion as `cvttsd2si` into a 64-bit register, low 32 bits, unsigned `cmovbe` cap, so NaN and infinite positions select half-bin 0; the distances are `0/0` or `inf/inf`, every intensity score is the default NaN `0xfff8000000000000`, no seed is found and the map is empty (exit 0). The same holds for an infinite retention time (`nonfinite_stage`); an infinite m/z makes the step infinite too but ends at step 3.1 (see the Size-conversion entry, CPP-314). `FileFilter_44_input.mzML` holds two MS1 spectra, not four, and is a short input (see CPP-312).

**Rust handling:** `DegenerateBinStep::Source` (the default and the tool's) reproduces the Release outcome bit for bit, NaN bits included, on Linux x86_64 and macOS arm64; `DegenerateBinStep::Refuse` returns `Error::InvalidValue` exactly when a step is zero or infinite and the seed loop visits a scan. The wave-5 port therefore no longer diverges here, and the `#[ignore]`d divergence test is gone. See also CPP-312, which is the same conversion reached from a short input.

## CPP-275 — charge_low above charge_high wraps the charge count

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Executed (Linux x86_64 Release, `bc9cc12-c19e494-174b576`).

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:197`, `run`.

**Trigger:** `isotopic_pattern:charge_low` more than one above `isotopic_pattern:charge_high`, or `charge_low` 1 with `charge_high` `INT_MAX`.

**Issue:** `UInt charge_count = charge_high - charge_low + 1` wraps for `charge_low > charge_high + 1`, and the `UInt` array count `3 + 2 * charge_count` (`.cpp:196-221`) wraps too, also for the positive count `2^31 - 1`. A count of `-1` or `2^31 - 1` wraps to 1 array and the loop writes arrays 1 and 2 past the end; `-2` and `-3` wrap to `2^32 - 1` and `2^32 - 3` arrays per spectrum, which stay in bounds if the allocation succeeds; `-4` and below wrap to `2^32 + 3 + 2 * count` (at least 9) and the pattern loop writes up to index `2^32 + 2 + count`, past the end.

**Proposed C++ fix:** Validate `charge_low <= charge_high` in `updateMembers_` and throw `InvalidParameter`.

**Evidence:** `../oracle/ffap-complete-fix3` `fix3_stage`, two runs each: `charge_low`/`charge_high` 4/2, `INT_MAX`/1, `INT_MAX`/498 and 1/`INT_MAX` die with SIGSEGV; 5/2, 6/2, 7/2 and 2/`INT_MAX` (no wrap, `2^32 - 1` arrays) throw `std::bad_alloc`; 3/2 returns an empty map. `../oracle/ffap-complete-fix6` `node/run_wrap.sh`, 20 counts twice each: which of the two happens is decided by the memory the process may have, not by the count. `sizeof(MSSpectrum::FloatDataArray)` is 88 and the array vector's `max_size()` is `(2^63 - 1) / 88`, so `resize` never throws `std::length_error` here; under a 16 GB address space every count up to 40,000,003 arrays (3.3 GiB) dies with SIGSEGV and every count from 80,000,003 (6.6 GiB) throws `std::bad_alloc`, while under a 500 GB one 100,000,003, 166,000,001, 200,000,003 and 1,000,000,003 arrays die instead and only `2^32 - 5` arrays still throw. `../oracle/ffap-complete-min6` `node/run_native_wrap.sh`, the round-6 minors, the same driver binary with no `ulimit -v` at all, twice each: 12,201,611, 12,201,613, 20,000,003, 40,000,003, 80,000,003, 100,000,003, 166,000,001, 200,000,003 and 1,000,000,003 arrays all die with SIGSEGV, so the capped `std::bad_alloc` rows are a property of the cap and not of the reference platform. `node/run_mem.sh` measures what the run holds when it writes past the end: the maximum resident set is 2.88 GiB at 12,201,611 arrays and 21.85 GiB at 100,000,003, about 232 bytes per array, because the pattern loop names and `assign`s every in-bounds array before the first out-of-bounds index; `2^32 - 5` arrays would need about 928 GiB on that model, 93% of the shared reference node's entire memory, and was not run uncapped.

**Rust handling:** `Settings::charge_count` refuses every wrapping count with `Error::InvalidValue`, whatever the `Limits` (the out-of-bounds cases as undefined, `-2` and `-3` before their allocation); counts up to `2^31 - 2` stay behind the native `Limits::max_charges`; `charge_low == charge_high + 1` gives zero charges, as in the source; the counts `-2` and `-3` stay refused unconditionally (lead decision D13 of wave 5); where the source writes out of bounds the run records a `DebugTermination` at `TerminationPoint::ScoreArrays`, as long as the one allocation between the wrap and that write - the first spectrum's `(3 + 2n) mod 2^32` arrays of 88 bytes each - is at most the bytes of 1,000,000,003 arrays (a crate constant, not a `Limits` field, so no caller can move it: lead decision D12) - the largest count measured to die on the reference node with the memory that node has, raised in the round-6 minors from the 1 GiB that round 6 read off runs under the harness's 16 GB address-space cap; above that line the executed outcome depends on the memory the process may have and nothing is recorded, as lead decision D6 has it for the isotope windows (executed: a reused object's `debug/log.txt` is left at the first run's flushed prefix after 4/2, `INT_MAX`/498 and 2141382846/1, and complete after 7/2's `std::bad_alloc`) (`charge_count_wraps_are_refused_whatever_the_limits`, `boundary_cases_match_the_linux_release_build`, `a_wrapped_score_array_count_records_its_termination_up_to_the_documented_ceiling`, `a_reused_instance_leaves_the_executed_log_at_every_later_termination`).

## CPP-276 — isotopeScore_ stops trying shorter isotope tails after a better fit

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; the constructed case runs in the port's test. Whether the narrowing is intended is unconfirmed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:1757–1788`, `isotopeScore_`.

**Trigger:** A pattern where a better fit with more trailing isotopes is found for an early `b`, and a later `b` fits best with fewer.

**Issue:** The inner loop starts at `e = best_end`, and `best_end` is raised inside the loop whenever a better fit is found. For every later `b` the combinations with fewer trailing isotopes are never tried, even when they fit better: in the constructed case the skipped candidate has correlation 0.850 and the returned one 0.746.

**Proposed C++ fix:** Start the inner loop at the `best_end` computed before the search, kept in a separate variable.

**Evidence:** Source review; `isotope_score_narrows_later_candidates_after_a_new_best_fit` in `tests/feature_finder_picked_seeds.rs`.

**Rust handling:** `isotope_score` reproduces the narrowing.

## CPP-277 — FeatureFinderAlgorithmPicked's empty-input message names ranges

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Message defect; source-reviewed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:1068–1071`, `run`.

**Trigger:** An input map with no peaks.

**Issue:** The check is `input_map.getSize() == 0`, which counts peaks, but the exception says "FeatureFinder needs updated ranges on input map. Aborting.", which describes neither the check nor the cause.

**Proposed C++ fix:** Report that the input holds no peaks.

**Evidence:** Source review.

**Rust handling:** The port keeps the message verbatim on the same check.

## CPP-278 — splitByFAIMSCV groups carry no ranges, so FeatureFinderCentroided fails on FAIMS input

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/IONMOBILITY/IMDataConverter.cpp:42–69`, `splitByFAIMSCV`; the consumer `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:242`.

**Trigger:** Any input with FAIMS compensation voltages, split for FeatureFinderCentroided.

**Issue:** Voltage groups are filled with `addSpectrum` and `updateRanges` is never called, so `spectrumRanges().byMSLevel(1)` throws `InvalidValue` "No ranges for this MS level" on every group, and FeatureFinderCentroided exits 8 on every FAIMS input.

**Proposed C++ fix:** Call `updateRanges()` on each group before returning.

**The whole C++ FAIMS fix for FeatureFinderCentroided** is this entry's
`updateRanges()` **plus** CPP-282 and CPP-283. With only this one applied the
tool runs to the end, and on any input where at least one cross-voltage merge
actually fires -- the ordinary multi-voltage case under the default
`-faims_merge_features true` -- it then writes an **empty** feature map,
because `mergeFAIMSFeatures` erases every feature that still carries unique
id 0 (CPP-282), and by `FeatureFinderCentroided.cpp:294-299` every feature of
a FAIMS run carries `FAIMS_CV`, so none is held back in the untouched
non-FAIMS group. It does **not** write an empty map when no merge fires: with
`-faims_merge_features false` the merge is never called
(`FeatureFinderCentroided.cpp:309`), and on single-voltage FAIMS input the
callback returns `false` for every pair (`FeatureOverlapFilter.cpp:445-448`),
so `removed_uids` stays empty and the `erase` at
`FeatureOverlapFilter.cpp:277-281` removes nothing. With CPP-282 also applied
the tool writes features, but a cluster of three or more voltages is split
into two features whose intensities double-count one member (CPP-283). The
Rust port applies all three:
`docs/TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md`, *The FAIMS closure*.

*Basis of the composite claim.* The three-fix chain is read from the pinned
source (`FeatureOverlapFilter.cpp:263-271` inserts into `removed_uids` only
when the callback returns `true`; `:277-281` erases exactly those ids;
`FeatureFinderCentroided.cpp:328-329` assigns the unique ids only *after* the
merge) and confirmed on the Rust port, whose source-faithful merge
(`FeatureOverlapFilter::merge_faims_features`) is call-for-call the source's.
Two executed tests of package B11 pin the two non-merging cases:
`a_single_faims_voltage_finds_the_features_of_the_plain_input` asserts the
console line `FAIMS feature merge: 8 -> 8 features (merged 0)` on
single-voltage input, and `two_faims_voltages_reproduce_the_release_build_group_by_group`
asserts that `-faims_merge_features false` prints no `FAIMS feature merge:`
line at all. **No patched C++ build was ever executed**: the pinned build
exits 8 before the merge on every FAIMS input (this entry), so the post-fix
behaviour of the C++ tool is read from its source, not observed -- in
particular nobody has seen a C++ run produce the empty map. The empty-map
outcome also rests on CPP-282's own premise, that the features still carry
unique id 0 when the merge runs.

**Evidence:** All 35 voltage groups in `../oracle/im-data-converter` throw at `byMSLevel(1)`, as do C2 `faims_facts` 1a and 1b (`../oracle/featurefinder-picked`); `tests/data/im_data_converter_provenance.json`.

**Rust handling:** Ranges are computed on demand, so a voltage group has its own ranges the moment it holds spectra; the defect is not reachable in the port and is not emulated. Re-executed for package B11 against the Release build with eight runs on six FAIMS inputs (`../oracle/b11-faims`), all rc 8 after the first `Processing FAIMS CV group:` line. `FeatureFinderCentroided` now processes FAIMS input end to end (`docs/TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md`, native difference 1).

## CPP-279 — splitByFAIMSCV destroys the chromatograms of FAIMS input

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/IONMOBILITY/IMDataConverter.cpp:84`, `splitByFAIMSCV`.

**Trigger:** FAIMS input that also holds chromatograms.

**Issue:** The groups receive spectra only, and `exp.clear(true)` then destroys the input's chromatograms, and any skipped spectra, without a message.

**Proposed C++ fix:** Move the chromatograms into a group (or return them) and report what was dropped.

**Evidence:** Oracle cases `faims_test_data` and `settings_and_chromatograms_faims` in `../oracle/im-data-converter`.

**Rust handling:** `FaimsSplit::dropped_chromatograms` and `skipped_spectra` return them. `FeatureFinderCentroided` loads MS level 1 only and uses no chromatogram, so it discards them as the source does; no output of that tool can show the difference.

## CPP-280 — splitByFAIMSCV groups NaN voltages unpredictably

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/IONMOBILITY/IMDataConverter.cpp:31` and `:42–64`, `splitByFAIMSCV`.

**Trigger:** A spectrum with unit FAIMS_CV and a NaN drift time.

**Issue:** The NaN reaches `std::set<double>` through `FAIMSHelper::getCompensationVoltages` (CPP-244) and the keys of `std::map<double, MSExperiment>`, whose ordering NaN breaks. With the NaN first, the whole input is returned unsplit under the NaN key with a spurious missing-voltage warning; with the NaN in the middle or last, `find` places the NaN spectrum in the unrelated -50 V group and the MS2 spectrum after it is skipped.

**Proposed C++ fix:** Reject or skip NaN voltages before inserting them, with a warning.

**Evidence:** Oracle cases `nan_first`, `nan_middle` and `nan_last` in `../oracle/im-data-converter`.

**Rust handling:** NaN voltages are refused with `Error::InvalidValue` and the input unchanged; the recorded C++ groups are asserted.

## CPP-281 — splitByFAIMSCV's information message says "Not" for "No"

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Message defect; source-reviewed.

**Affected file/function:** `src/openms/source/IONMOBILITY/IMDataConverter.cpp:35`.

**Trigger:** Input without FAIMS compensation voltages.

**Issue:** The message reads "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN."

**Proposed C++ fix:** "No FAIMS compensation voltages found in the data. ..."

**Evidence:** Source review; the executed oracle prints the same text.

**Rust handling:** `ImDataConverter::NO_COMPENSATION_VOLTAGES_INFO` keeps the text verbatim.

## CPP-282 — mergeFAIMSFeatures removes every FAIMS feature without a unique ID

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp:200–281`, `filter`, as called by `mergeFAIMSFeatures` (`:384–527`).

**Trigger:** FAIMS features whose unique IDs are 0, as FeatureFinderAlgorithmPicked returns them before a tool assigns IDs.

**Issue:** Removal is recorded in `removed_uids` by `getUniqueId()`. Once one feature with ID 0 is merged, every feature with ID 0 is skipped as a querier and erased, so the result holds no FAIMS feature at all.

**Proposed C++ fix:** Record removal by index or pointer instead of unique ID.

**Evidence:** Oracle case `c2_uid0_wipe` in `../oracle/feature-overlap-filter`, and C2 `faims_facts` (0 features); `tests/data/feature_overlap_filter_provenance.json`.

**Rust handling:** Reproduced exactly by `FeatureOverlapFilter::merge_faims_features` and by the executed case `c2_uid0_wipe`. Package B11 added the corrected route the tool takes: `FeatureFinderCentroided` draws a unique id for every feature before it merges -- the ids its own `applyMemberFunction(setUniqueId)` overwrites a few lines later -- so removal keys on real ids, and `FaimsMergeFidelity::Corrected` refuses a map whose FAIMS features share an id instead of erasing them.

## CPP-283 — FeatureOverlapFilter merges already removed features again

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp:205–271`, `filter`; the merge callback at `:430–499`.

**Trigger:** Three or more mutually overlapping features.

**Issue:** The inner loop does not skip candidates that are already in `removed_uids`, so a removed feature is passed to the callback again for a later survivor and its intensity is added twice. In `mergeFAIMSFeatures` a survivor absorbs at most one feature, because the callback removes its `FAIMS_CV` after the first merge and then refuses: 1000, 900 and 800 at three voltages become 1900 and 1700.

**Proposed C++ fix:** Skip candidates already marked removed, and let a survivor keep absorbing features after its first merge.

**Evidence:** Oracle case `c2_three_cvs` in `../oracle/feature-overlap-filter`.

**Rust handling:** Reproduced exactly by `FeatureOverlapFilter::merge_faims_features` and by the executed case `c2_three_cvs` (1900 and 1700). Package B11 added `FaimsMergeFidelity::Corrected` beside it, which skips a candidate already marked removed and tests a candidate's voltage against the voltages the survivor already stands for (`merged_centroid_IMs`) rather than against a `FAIMS_CV` the first merge removed, so 1000, 900 and 800 at three voltages become one feature of 2700. `FeatureFinderCentroided` runs the corrected merge; the source one stays available and stays tested.

## CPP-284 — FeatureOverlapFilter's float boxes miss pairs within tolerance

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp:144–191`, the `getBox` lambdas and the quadtree extent; `src/openms/extern/Quadtree/include/Box.h:60–64`, `intersects`.

**Trigger:** A zero tolerance, a tolerance below the `float` spacing of the coordinates, or coordinates near 2^24.

**Issue:** Candidates come from `Box<float>` intersection, which is strict, while the final test is an inclusive `double` distance. Two features at the same position with tolerance 0 have zero-width boxes that never intersect, and near 2^24 the `float` conversion collapses boxes, so pairs that satisfy the tolerance test never merge.

**Proposed C++ fix:** Build the boxes in `double` with a small outward margin, or use inclusive intersection for the candidate query.

**Evidence:** Oracle cases `bound_zero_tolerance`, `bound_tiny_tolerance` and `bound_float_collapse_rt` in `../oracle/feature-overlap-filter`.

**Rust handling:** Reproduced exactly.

## CPP-285 — Trace-level overlap compares only start times and drops some traces

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp:59–87`, `getFeatureBounds`.

**Trigger:** `filter` in trace mode.

**Issue:** The end of a mass trace is searched from the end of the hull outline, which returns to the first scan, so `rt_max` equals the start time and trace overlap compares start times only. A trace whose first scan's lower m/z is at or below 0 gets `rt_min` after `rt_max` and is silently skipped.

**Proposed C++ fix:** Take `rt_max` from the largest retention time with positive m/z on the outline.

**Evidence:** Oracle cases `edge_trace_bounds_collapse_trace` (two features whose traces overlap for 1.5 s are kept) and `edge_trace_inverted_bounds_skipped` with its control in `../oracle/feature-overlap-filter`.

**Rust handling:** Reproduced exactly.

## CPP-286 — FeatureOverlapFilter aborts a Debug build on valid feature maps

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK (Debug) at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed (Debug only).

**Affected file/function:** `src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp:139` and `:171–191`; `src/openms/source/KERNEL/FeatureMap.cpp:260`, `updateRanges`; `src/openms/extern/Quadtree/include/Quadtree.h:141`.

**Trigger:** A feature with a zero-width hull box outside the extent margin, or a feature without a hull in the hull modes.

**Issue:** `FeatureMap::updateRanges` skips empty or zero-width hull boxes (`DBoundingBox::isEmpty`), so the quadtree extent can exclude a feature's box, and `assert(box.contains(mGetBox(value)))` aborts. A hull-less feature converts the `±DBL_MAX` sentinel box to `float` and aborts the same way; a Release build silently ignores it.

**Proposed C++ fix:** Compute the extent from the same boxes `getBox` returns, and reject features without a hull in the hull modes.

**Evidence:** Oracle cases `edge_zero_extent_hull_outside_margin` and `edge_hull_less_convex_hull` (exit 134) in `../oracle/feature-overlap-filter`, compared with a Release replica of the pinned source.

**Rust handling:** Follows the Release replica; a hull-less feature in the hull modes is `Error::MissingInformation`.

## CPP-287 — FeatureOverlapFilter leaves the map changed when it throws

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the product SDK at core `4fdec46b205459b92e7d3b9e56df5d8e912d5c85`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp:141` and `:197` (`filter`); `:405–527` (`mergeFAIMSFeatures`).

**Trigger:** A subordinate without a convex hull in trace mode; a `FAIMS_CV` or merged list the callback cannot convert.

**Issue:** `filter` sorts the caller's map before `getFeatureBounds` throws `MissingInformation`, so the map comes back reordered. `mergeFAIMSFeatures` moves every feature into two temporary maps first; when its callback throws `ConversionError`, the caller's map is left holding moved-from features stripped of their metadata.

**Proposed C++ fix:** Validate before sorting and moving, or restore the map in a catch block.

**Evidence:** Oracle cases `edge_trace_missing_sub_hull` and `edge_faims_cv_empty` in `../oracle/feature-overlap-filter`.

**Rust handling:** Every error leaves the map exactly as it was (a journal of the overwritten values, or a copy where a later failure is possible).

## CPP-288 — FeatureOverlapFilter has undefined behaviour on reachable inputs

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; not executed.

**Affected file/function:** `src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp:40–43` (`getFeatureBounds`), `:115–120` (`tracesOverlap`), `:232–233` and `:442–443` (the `FAIMS_CV` conversions).

**Trigger:** Trace mode with a feature that has more subordinates than hulls, or a candidate without trace bounds; a string or list `FAIMS_CV`.

**Issue:** `feat.getConvexHulls()[i]` indexes past the hull vector when a feature has more subordinates than hulls, and `points.front()` reads an empty hull. `tracesOverlap` dereferences `find(...)->second` without checking for `end()`, which a feature whose traces were all skipped reaches. `(double)` on a string or list `DataValue` reads the wrong union member instead of throwing.

**Proposed C++ fix:** Check the hull count and emptiness, check `find` against `end()`, and convert `FAIMS_CV` with a type check.

**Evidence:** Source review.

**Rust handling:** Each case is refused with `Error::InvalidValue` or `Error::MissingInformation` before the map changes.

## CPP-289 — GaussTraceFitter's copy leaves region_rt_span_ uninitialised

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; not executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/GaussTraceFitter.cpp:25-46`, the copy constructor and `operator=`.

**Trigger:** Any copy of a `GaussTraceFitter` that has fitted, followed by `checkMaximalRTSpan`.

**Issue:** Both copy `height_`, `x0_` and `sigma_` and call `updateMembers_()`, but neither copies `region_rt_span_`, which `setInitialParameters_` sets (`.cpp:259-260`) and which `checkMaximalRTSpan` divides the fitted span against (`.cpp:98-101`). The copy therefore reads an indeterminate value, and `FeatureFinderAlgorithmPicked` treats a `true` answer as "Invalid fit: Fitted model is bigger than 'max_rt_span'".

**Proposed C++ fix:** Copy `region_rt_span_` in both members, or declare them `= default` now that every member is copyable.

**Evidence:** Source review of the pinned file. The port's own `Clone` copies the whole model, so no executed case distinguishes them.

**Rust handling:** Every field of the model is cloned; `region_rt_span` is part of the fitted state.

## CPP-290 — GaussTraceFitter's Jacobian carries 0.125 in the sigma column

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/GaussTraceFitter.cpp:199`, `GaussTraceFunctor::df`.

**Trigger:** Every Gaussian trace fit.

**Issue:** The analytic derivative of `theoretical_int * height * exp(-(rt-x0)^2 / (2 sigma^2))` with respect to `sigma` is `theoretical_int * height * e * (rt-x0)^2 / sigma^3`; the code writes `0.125 * trace.theoretical_int * height * e * pow2(rt - x0) * inv_sig3 * weight`. The extra factor of 1/8 does not stop the fit converging on the tested data, but the Levenberg-Marquardt trust region, and with it every evaluation-budget boundary, follows a wrong derivative.

**Proposed C++ fix:** Drop the `0.125`, and re-derive the published class-test parameters, which were produced with it.

**Evidence:** The executed Jacobian columns of `../oracle/gauss-trace-fitter` carry the factor; comparing the recorded column against the analytic derivative at the same vectors reproduces it exactly.

**Rust handling:** Reproduced: `gauss_trace_fitter.rs` transcribes the `0.125`, with a comment naming this issue, because removing it would change every fitted parameter.

## CPP-291 — An all-zero intensity profile "fits" with a NaN sigma

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/GaussTraceFitter.cpp:253-291` (`setInitialParameters_`) and `src/openms/source/FEATUREFINDER/TraceFitter.cpp:122-130` (`optimize_`).

**Trigger:** A mass trace whose intensities are all zero, or whose smoothed maximum equals the baseline.

**Issue:** `height_ = smoothed[max_index] - traces.baseline` is 0, so `alpha = (left_height + right_height) * 0.5 / height_` is `0/0 = NaN`. `if (alpha >= 1)` is false for NaN, so the guard the line above provides for the degenerate case is skipped and `sigma_ = delta_x * 0.5 / sqrt(-2.0 * log(alpha))` is NaN. Every residual is then NaN, Eigen stops with `CosinusTooSmall` (status 5), and `optimize_` accepts every status above `ImproperInputParameters`, so `fit` succeeds with a NaN model, NaN bounds, NaN FWHM and a NaN area.

**Proposed C++ fix:** Test `alpha` for NaN as well as for `>= 1`, or throw `UnableToFit` when `height_` is not positive.

**Evidence:** Executed case `start.all_zero` in `../oracle/gauss-trace-fitter` (driver linked against the product SDK, run twice with identical output).

**Rust handling:** Reproduced under the default compatibility, and available as a refusal: the port records the same NaN model so the executed comparison holds, and documents it as a native boundary.

## CPP-292 — The trace fitters read empty traces and index theoretical peaks unchecked

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; not executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/GaussTraceFitter.cpp:229-262`, `src/openms/source/FEATUREFINDER/EGHTraceFitter.cpp:333-368` (both `setInitialParameters_`), and `src/openms/source/FEATUREFINDER/TraceFitter.cpp` `computeTheoretical`.

**Trigger:** `MassTraces` whose intensity profile is empty (no trace has a peak), or a `k` beyond `trace.peaks.size()`.

**Issue:** With an empty profile `N` is 0, so `smoothed` is empty and `smoothed[max_index]` with `max_index == 0` reads out of bounds; `std::advance(it, max_index)` then walks an empty list and `total_intensities.rbegin()->first` dereferences its end. `computeTheoretical` indexes `trace.peaks[k]` with no bound check. All of it is undefined behaviour reachable from the public API.

**Proposed C++ fix:** Throw `UnableToFit` when the profile is empty, and bound-check `k`.

**Evidence:** Source review of the pinned files.

**Rust handling:** `Error::InvalidValue` (`UnableToFit-FinalSet`) before anything is read, in both fitters and in `compute_theoretical`.

## CPP-293 — ParamEntry::isValid narrows a 64-bit integer before its range check

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Executed (the acceptance); the restriction bypass is source-reviewed.

**Affected file/function:** `src/openms/source/DATASTRUCTURES/Param.cpp:120-127`, `ParamEntry::isValid`.

**Trigger:** An `INT_VALUE` parameter outside the `int` range, with or without `setMinInt`/`setMaxInt`. FeatureFinderAlgorithmPicked reaches it through every integer parameter it declares.

**Issue:** `int tmp = value;` narrows the 64-bit `ParamValue` before the `min_int`/`max_int` comparison, so the check runs on the truncated value: `2^32 + 5` passes `setMinInt(1)` as `5`, and the error message prints the truncated number too.

The same narrowing runs on the reading side: `updateMembers_` reads the members through `ParamValue::operator int` / `operator unsigned int`, which keep the low 32 bits (`ParamValue.cpp:445-461`). The consequences for a FeatureFinderAlgorithmPicked caller, all executed: `intensity:bins = 2^32 + 10` runs with 10 bins (10 features on FFC_1); `2^32` is refused as `'0'` and `2^31` as `'-2147483648'`; `min_spectra_ = (UInt) floor(value * 0.5)` keeps the low 32 bits of `cvttsd2si`, so `2^33 + 22` gives 11 and `2^62 + 30` gives 0. A negative value whose low 32 bits pass the check throws `ConversionError` from `operator unsigned int` either half way through `updateMembers_`, leaving the earlier members updated (`mass_trace:max_missing`, `intensity:bins`), or at the start of `run_` (`fit:max_iterations`).

**Proposed C++ fix:** Compare the 64-bit value, and reject anything outside the `int` range explicitly.

**Evidence:** The executed probe `../oracle/gauss-trace-fitter/param-range` shows `TraceFitter` accepting a `max_iteration` beyond `int`; that parameter carries no restriction, so it shows the acceptance and not the bypass. `../oracle/ffap-complete-fix2` `fix2_driver bigint`, 21 cases, two runs each, identical, executes the bypass and every consequence listed above on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Rust handling:** `Param` refuses a value that cannot be converted to `i32` in the restriction check, so `TraceFitterParams::to_param`/`from_param` round-trip only within `i32`; `tests/trace_fitter.rs::to_param_and_from_param_disagree_beyond_i32` pins the difference and must flip when this is fixed. FeatureFinderAlgorithmPicked is the exception: `algorithm::check_parameters` and `Settings` reproduce all of it with the source's texts, because that header is pinned tier 1 against the Release build. The crate's shared `Param` check keeps its native refusal for every other handler.

## CPP-294 — TraceFitter documents both RT-span checks backwards and calls maxfev an iteration count

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed (documentation).

**Affected file/function:** `src/openms/include/OpenMS/FEATUREFINDER/TraceFitter.h:170-191` (the two span checks), `:31-33` and `:253` (`max_iteration`, `max_iterations_`), and `src/openms/source/FEATUREFINDER/TraceFitter.cpp:129`.

**Trigger:** Reading the header to implement or call a `TraceFitter`.

**Issue:** The header says `checkMinimalRTSpan` returns true "when the model spans at least `min_rt_span` of the search area" and `checkMaximalRTSpan` true "when the model does not exceed `max_rt_span`". Both implementations return the opposite (`GaussTraceFitter.cpp:98-106`), and `FeatureFinderAlgorithmPicked.cpp:2055` and `:2081` treat a `true` answer as the failure. `max_iteration` is documented as "maximum number of LM iterations" but `optimize_` assigns it to `lmSolver.parameters.maxfev` (`TraceFitter.cpp:121`), Eigen's maximum number of *function evaluations*, which is a different and much smaller budget. Finally the `UnableToFit` text is "Could not fit the gaussian to the data" for every model, including EGH.

**Proposed C++ fix:** Reword both `@return` clauses to the implemented sense, rename the parameter or document it as `maxfev`, and take the model name from the subclass.

**Evidence:** Source review of the pinned header and implementations, against the two executed call sites in `FeatureFinderAlgorithmPicked.cpp`.

**Rust handling:** The port documents the implemented sense, names the field `max_iteration` with the `maxfev` meaning stated, and reproduces the "gaussian" wording verbatim in both fitters so the executed message comparison holds.

## CPP-295 — EGHTraceFitter::setInitialParameters_ has no alpha >= 1 guard

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/EGHTraceFitter.cpp:395-412`, `setInitialParameters_`.

**Trigger:** A flat or single-scan profile, or a baseline above the smoothed half-height edges, so that `alpha = (left_height + right_height) * 0.5 / height_` is at least 1.

**Issue:** `log_alpha = log(alpha)` is then zero or positive, `tau_ = -1 / log_alpha * (B - A)` is infinite or has the wrong sign, and `sigma_ = sqrt(-0.5 / log_alpha * B * A)` is the square root of a negative number, i.e. NaN. `optimize_` accepts `CosinusTooSmall` at the start point, so `fit` succeeds with a NaN sigma and a meaningless tau, and the bounds, FWHM and area computed from them are NaN. `GaussTraceFitter` guards exactly this case (`GaussTraceFitter.cpp:284-287`).

**Proposed C++ fix:** Add the same `alpha >= 1` guard, or throw `UnableToFit`.

**Evidence:** Executed C2 cases `flat3_egh` and `short3_egh`, and this package's `alpha_above_one_baseline`, in `../oracle/egh-trace-fitter` and `../oracle/featurefinder-picked`.

**Rust handling:** Reproduced, so the executed comparison holds; documented as a native boundary.

## CPP-296 — EGHTraceFitter::setInitialParameters_ reads an empty trace

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; not executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/EGHTraceFitter.cpp:333-368`.

**Trigger:** `MassTraces` whose computed intensity profile is empty.

**Issue:** `smoothed` is empty, so `smoothed[max_index]` reads out of bounds; `std::advance(it, max_index)` walks an empty list and `total_intensities.rbegin()->first` dereferences its end. Undefined behaviour, unlike the `N <= LEN + 1` short path `GaussTraceFitter` has.

**Proposed C++ fix:** Throw `UnableToFit` before smoothing when the profile is empty.

**Evidence:** Source review of the pinned file.

**Rust handling:** `Error::InvalidValue` (`UnableToFit-FinalSet`).

## CPP-297 — EGHTraceFunctor drops the baseline in one branch and takes |sigma| in one of two places

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed, with executed rows.

**Affected file/function:** `src/openms/source/FEATUREFINDER/EGHTraceFitter.cpp:62-71` (`operator()`), `:88` (`df`) and `:277-288` (`getGnuplotFormula`).

**Trigger:** Any fit whose trial parameters make `2 sigma^2 + tau (t - t_R)` non-positive, or whose sigma is negative.

**Issue:** Where the denominator is not positive the residual model is `fegh = 0.0`, while the positive branch is `baseline + theoretical_int * H * exp(...)` and `getGnuplotFormula` writes `baseline + (cond ? ... : 0)`. The plotted curve therefore falls to the baseline where the residual falls to zero. Separately, `df` uses `fabs(x_map(2))` for sigma while `operator()` uses the signed value, so for a negative sigma the sigma column of the Jacobian has the wrong sign and points away from the descent direction.

**Proposed C++ fix:** Use the baseline in the else branch, and treat sigma consistently in both members (or project sigma to its absolute value once, before either is called).

**Evidence:** Source review of the pinned file, with the executed functor rows of `../oracle/egh-trace-fitter` covering both branches.

**Rust handling:** Reproduced exactly, including the sign asymmetry, so the executed residual and Jacobian comparison holds to the last place.

## CPP-298 — extendMassTraces_ compares a pattern index with a trace index

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; the resulting output is executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:1460-1476`, `extendMassTraces_`.

**Trigger:** An isotope of the pattern whose extended trace has fewer than three peaks.

**Issue:** `p` indexes the *isotope pattern*, while `traces.max_trace` indexes the *traces collected so far* and is still 0 until the maximum trace has been pushed. So an invalid trace at `p == 0` takes neither branch and is appended anyway; any later invalid trace before the maximum breaks out of the loop; and the `traces.clear()` branch the comment describes ("Missing traces in the middle of a pattern are not acceptable") is unreachable, because `p < traces.max_trace` cannot hold while `max_trace` is 0.

**Proposed C++ fix:** Compare `p` against the pattern's own maximum position, or compare the trace index the loop is filling.

**Evidence:** Source review of the pinned file; the resulting traces are part of the executed feature output compared in `../oracle/b7-ffap-features` and against the Release build's `TOPP_FeatureFinderCentroided_1`.

**Rust handling:** Reproduced exactly, because it decides which traces a feature keeps.

## CPP-299 — The better-seed search re-reads a moving m/z

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; the resulting output is executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:1420-1445`, `extendMassTrace_`.

**Trigger:** A seed with a better nearby maximum in an adjacent scan.

**Issue:** `mz` is captured before the loop, but the nearest-peak lookup inside it is `map_[spectrum_index].findNearest(map_[starting_peak.spectrum][starting_peak.peak].getMZ())`, which re-reads the m/z of the *current* starting peak after `starting_peak` has moved. The search target therefore drifts with the accepted candidates while the acceptance window `std::fabs(mz - ...) >= pattern_tolerance_` stays anchored to the original m/z, so the two disagree about what is being searched for.

**Proposed C++ fix:** Search from the captured `mz`, or move the window with the seed; either is consistent, the mixture is not.

**Evidence:** Source review of the pinned file; the seeds this loop produces are part of the executed comparison.

**Rust handling:** Reproduced exactly.

## CPP-300 — The monoisotopic m/z correction uses the proton mass and a trace index

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; the reported value is executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:780-784`, `run_`, the `feature:reported_mz == "monoisotopic"` branch.

**Trigger:** Any feature reported with `feature:reported_mz` set to `monoisotopic`.

**Issue:** The correction is `(Constants::PROTON_MASS_U / c) * (traces.getTheoreticalmaxPosition() + trimmed_left)`. The spacing between isotope peaks is `C13C12_MASSDIFF_U` (about 1.00335 u), not `PROTON_MASS_U` (about 1.00728 u), and `getTheoreticalmaxPosition()` is a position in the pattern, not an isotope number. The reported monoisotopic m/z is biased by roughly 4 mDa per isotope step at charge 1.

**Proposed C++ fix:** Subtract `C13C12_MASSDIFF_U / charge` per isotope step, counted from the monoisotopic peak.

**Evidence:** Source review of the pinned file; the reported value is part of the executed feature comparison.

**Rust handling:** Reproduced exactly, with the constant named in the module documentation.

## CPP-301 — The final feature intensity picks the isotope window by m/z, not by mass

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; the reported value is executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:790` against `:356`.

**Trigger:** Any feature of charge 2 or more.

**Issue:** The intensity is `fitter->getArea() / getIsotopeDistribution_(f.getMZ()).max`, but the isotope windows were precalculated over *mass*: `max_mass = maxMZ * charge_high` at `:356`, and `getIsotopeDistribution_` bins its argument by `mass_window_width_`. For a doubly charged feature the window used is the one for half the feature's mass, whose maximum abundance differs, so the reported intensity is scaled by the wrong factor.

**Proposed C++ fix:** Pass `f.getMZ() * charge` (the mass) to `getIsotopeDistribution_`.

**Evidence:** Source review of the pinned file; the reported intensity is part of the executed feature comparison.

**Rust handling:** Reproduced exactly.

## CPP-302 — aborts_ is written from inside the OpenMP region without synchronisation

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:1129-1140` (`abort_`), called at `:627`, `:640` and `:725` inside the `#pragma omp parallel for` opened at `:595`.

**Trigger:** Two threads aborting a seed at the same time.

**Issue:** `abort_` does `aborts_[reason]++` on a `std::map<String, UInt>` and, with `debug_`, `abort_reasons_[seed] = reason` on a second map, both without a `critical` section, while the surrounding loop is parallel. Every other shared write in that loop is guarded (`:651`, `:716`, `:798`, `:812`). Concurrent insertion into a `std::map` is undefined behaviour, and the counts the tool prints come from it.

**Proposed C++ fix:** Put both writes in a `critical` section, or accumulate per thread and merge after the loop.

**Evidence:** Source review of the pinned file. The C2 driver records the library's abort map single-threaded only, for this reason. Executed: at `OMP_NUM_THREADS=4`, `debug/log.txt` differed in each of three runs (`../oracle/ffap-instr-completion` tool cases `c1`, `c2`). A four-thread run with no seed (`c3`) is defined and equals the single-thread run.

**Rust handling:** The port returns `RunOutput::aborts` from the serial merge of per-seed results, so the counts are deterministic at every thread count (asserted byte-identical at `-threads` 1/2/4/8/0); `abort_reasons_` and the debug log are collected in seed order as well. Giving the single-thread result for this race is the one documented exception to wave 5's undefined-behaviour rule (lead decision D11), because the determinism contract requires parallel output to equal serial output.

## CPP-303 — MorphologicalFilter leaves the last output sample unwritten for a one-sample element

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the C++ Release build `openms4-release-bc9cc12-c19e494-174b576` (gcc 14.4, `-O3 -DNDEBUG`).

**Status:** Executed.

**Affected file/function:** `src/openms/include/OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h:324-427` (`applyErosion_`) and the matching `applyDilation_`.

**Trigger:** `struc_size == 1` (one data point, or a Thomson width narrower than the spacing) on a signal of more than five samples, so the simple fallback at `:341` is not taken.

**Issue:** With `struc_size_half == 0` the lower-margin loop is empty, the middle loop writes output indices `0 .. size - 2`, and the higher-margin block's loops are empty as well, so `output[size - 1]` is never assigned. The caller sees whatever the output buffer held. Through `filterRange` that gives, for a one-sample element: `erosion`, `dilation`, `opening` and `closing` return the input with the last sample zeroed; `tophat` and `bothat` return zero everywhere with the last sample kept; `gradient` returns zero everywhere with the last sample the negated stale buffer value (see CPP-304). The identity of a one-sample erosion or dilation is the input itself, so every one of these is wrong in exactly one sample.

**Proposed C++ fix:** Write the last output sample in the higher-margin block, or route `struc_size == 1` to `applyErosionSimple_`/`applyDilationSimple_`, which are correct there.

**Evidence:** Executed differential against the Release build over the edge shapes and an exhaustive element-length sweep, and at tool level through `BaselineFilter`: `../oracle/baseline-filter-edges`, with the 1,698 committed expectation rows in `tests/data/baseline_filter_edges_expected.tsv` and the 528 sweep rows.

**Rust handling:** Reproduced exactly, including through the `BaselineFilter` tool; `-method erosion_simple` and `dilation_simple` select the simple variants, which keep the last sample, as the source does.

## CPP-304 — MorphologicalFilter::filterRange's static scratch buffer leaks between calls

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the C++ Release build `openms4-release-bc9cc12-c19e494-174b576` (gcc 14.4, `-O3 -DNDEBUG`).

**Status:** Executed.

**Affected file/function:** `src/openms/include/OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h:171` (the function-local `static std::vector` in `filterRange`), used by the `gradient`, `tophat`, `bothat`, `opening` and `closing` branches at `:208-225`.

**Trigger:** A `gradient` filter with a one-sample structuring element, applied to more than one spectrum in the same process.

**Issue:** The buffer is `static` "only to avoid reallocation", but `applyErosion_` does not write its last sample for a one-sample element (CPP-303), so `buffer[size - 1]` still holds the value a *previous* call left there. `output_begin[i] -= buffer[i]` then subtracts it, and the last sample of a gradient depends on which spectra were filtered before it, in which order, in the same process. Two runs that differ only in the order of unrelated spectra give different output.

**Proposed C++ fix:** Clear the buffer before use, size it per call, or fix CPP-303 so nothing stale is ever read.

**Evidence:** Executed: the `gradient` rows of `../oracle/baseline-filter-edges` with a one-sample element, and the history fixture `tests/data/baseline_filter_edges_history.mzML`, which reproduces the dependence on the preceding spectrum.

**Rust handling:** Reproduced deliberately, so the executed comparison holds; the behaviour and its history dependence are documented in `docs/MORPHOLOGICAL_FILTER_SUPPORT.md`.

## CPP-305 — The indexedmzML indexListOffset addresses the byte before <indexList

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the C++ Release build `openms4-release-bc9cc12-c19e494-174b576` (gcc 14.4, `-O3 -DNDEBUG`).

**Status:** Executed.

**Affected file/function:** `src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp:88-125`, `writeFooter_`.

**Trigger:** Any `indexedmzML` file the C++ writer produces (`PeakFileOptions::getWriteIndex()`, the default).

**Issue:** `Int64 indexlistoffset = os.tellp();` is taken *before* `os << "\n";`, so the recorded offset points at the newline that precedes `<indexList`, one byte early. The mzML indexedmzML schema defines the value as the offset of the `<indexList>` element. A reader that seeks to it and expects a `<` finds `\n`.

**Proposed C++ fix:** Take `tellp()` after the newline is written.

**Evidence:** Executed on the Release build: its `MapNormalizer` output declares `<indexListOffset>5150584</indexListOffset>` while `<indexList` starts at 5150585; verified with an independent offset checker on the produced files (`../oracle/mzml-writer-scale-parity`, `tools/mzml_writing/check_output.py`).

**Rust handling:** This port's offsets address the opening `<indexList` exactly, and the difference is recorded as a container fact in `docs/MZML_WRITER_CPP_PARITY.md`, next to the placeholder checksum of CPP-049, so nobody later "fixes" the Rust offsets to match.

## CPP-306 — Converting a 32-bit time array from minutes narrows the result back to float

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the C++ Release build `openms4-release-bc9cc12-c19e494-174b576` (gcc 14.4, `-O3 -DNDEBUG`).

**Status:** Executed.

**Affected file/function:** `src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp:217-222`, `decodeBase64Arrays`; the unit is set at `:365-373`, `handleBinaryDataArrayCVParam`.

**Trigger:** Reading any mzML binary data array that carries `MS:1000595` (time array) with `unitAccession="UO:0000031"` (minute) and 32-bit precision — what ProteoWizard writes for the TIC chromatogram of a Thermo run.

**Issue:** The minute-to-second conversion is applied in place. The 64-bit branch (`:210-216`) multiplies a `std::vector<double>` and keeps the `double`; the 32-bit branch is `for (auto& it : bindata.floats_32) { it = it * unit_multiplier; }`, where `it` binds to `float&`, so the `double` product is narrowed back to `float` on assignment and the converted seconds keep only 32-bit precision. A Numpress array is forced to `PRE_64` before either branch (`:185`), so it is unaffected. Reading the same physical times as a 64-bit minute array therefore gives the C++ itself a different result from reading them as a 32-bit minute array: on the 40,856-point TIC time array of `profile_hr_qe_silac_uk222/UK222.mzML`, 38,107 of 40,856 times differ between the two, by up to 2.44e-4 s (half an `f32` ULP at 4,400 s), and 8,806 of 40,855 point spacings move by more than 1e-3 relative. Downstream that is not cosmetic: `PeakPickerHiRes` finds its apex by bisection on a cubic spline through those points, so the picked chromatogram of that file moves by up to 3.19e-3 s in retention time and 1.75e-3 relative in intensity.

**Proposed C++ fix:** Promote the array to `floats_64` before applying a multiplier other than 1.0, as the Numpress path already does, or accumulate the product in `double` and store it in a 64-bit array.

**Evidence:** Executed on the Release build, both sides, by the fixing lane and independently by its verifier. The first value of that file's time array, raw `0x1.0081c4p+0` minutes, is stored as `0x1.e0f35p+5` = 60.118804931640625 s from the 32-bit minute array and as `0x1.e0f34f8p+5` = 60.118803977966309 s from a 64-bit array holding the same times; a 32-bit **seconds** array reproduces the narrowed value exactly, which pins the loss to `f32` precision and nothing else. Seven cases over 32/64-bit against minute/second plus the empty, single-point and unsorted shapes were read and picked by the Release build; the driver, the case generator and the raw results are in `../oracle/picked-chromatogram/`, and the projected values are `tests/data/peak_picking/chromatogram_time_oracle.tsv`.

**Rust handling:** Reproduced deliberately, because it decides ordinary output: `mzml::ReadOptions::source_time_array_precision`, which `ReadOptions::source()` sets and therefore every tool path that reproduces source loading, narrows the product exactly as the source does; the library default keeps the `f64`. The same TSV pins both modes — the `min32` rows the source mode, the `min64` rows the default. One deliberate divergence: a finite 32-bit time whose product with 60 is not finite is refused with the existing `nonfinite binary value` parse error instead of being stored as the source's infinity. See `docs/MZML_SUPPORT.md` and `docs/PEAK_PICKING_SUPPORT.md`.

## CPP-307 — The XML writers print doubles at 15 significant digits, so their own output does not round-trip

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the C++ Release build `openms4-release-bc9cc12-c19e494-174b576` (gcc 14.4, `-O3 -DNDEBUG`).

**Status:** Executed.

**Affected file/function:** `src/openms/source/FORMAT/XMLFile.cpp:362` and `:369`, `save_` (both the compressed and the uncompressed stream), and `src/openms/source/FORMAT/MzMLFile.cpp:170`, `storeBuffer`; the precision comes from `writtenDigits<double>()` = `std::numeric_limits<double>::digits10` = 15 (`src/openms/include/OpenMS/CONCEPT/Types.h:188-192`).

**Trigger:** Any `double` an XML writer prints through the stream's default formatting — for mzML, every `scan start time` and every other scalar written as text rather than as a binary array.

**Issue:** 15 significant digits is `digits10`, the number of decimal digits a `double` is guaranteed to *carry*, not the `max_digits10` = 17 needed to *recover* it. The comment beside the call says "set high precision for floating point numbers", but the value chosen loses information: a store followed by a load does not return the value that was stored. This is not a comparison tolerance — it is the writer discarding bits that the reader then cannot restore.

**Proposed C++ fix:** Use `std::numeric_limits<double>::max_digits10` (17), or write the shortest round-tripping representation (`std::to_chars` with no precision argument).

**Evidence:** Executed on the Release build's own `PeakPickerHiRes` output for `profile_hr_qe_silac_uk222/UK222.mzML`. All 40,856 scan start times were extracted from the input and reduced with `60.0 * StringUtils::toDouble(s)`, then rendered with `os.precision(writtenDigits(double()))`: the probe's text equals the text the tool actually wrote, 40,856 of 40,856, and the C++'s own written text fails to reparse to its own stored `double` in 10,671 of 40,856, worst case stored `0x1.00054ab606b7ap+12`, written `4096.33074`, reparsed `0x1.00054ab606b7bp+12`, a difference of 9.09e-13 s. Verified independently by the `fix/picked-chromatogram` verifier, who also confirmed the two *readers* agree bit for bit on all 40,856 (so the residual is writer-side only).

**Rust handling:** This port writes the shortest round-tripping text and round-trips its own output 40,856 of 40,856, so a decoded comparison against the C++ shows the 1-to-2-ULP spread above on retention times while the binary arrays are bit-identical. The difference is recorded as a native difference in `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md` and in `docs/BENCHMARKS.md` §3.6, so nobody later "fixes" the Rust text to match.

## CPP-308 — MapNormalizer divides every MS1 intensity by an unguarded maximum

**Source revision:** TOPP `174b576e244e100f2345ca57a8e79aaa607156df`, the TOPP tree of the Release build `openms4-release-bc9cc12-c19e494-174b576` that the behaviour below was executed on. The source text was read from the two local TOPP checkouts `d0234cc` and `6f8eb94`, which agree line for line at the lines cited.

**Status:** Executed.

**Affected file/function:** `src/MapNormalizer.cpp:93-103`, `TOPPMapNormalizer::main_`.

**Trigger:** Any input whose combined maximum intensity — `MSExperiment::getMaxIntensity()` after `updateRanges()`, which includes chromatograms — is zero, negative, or an empty range.

**Issue:** The tool computes `exp.updateRanges(); double max = exp.getMaxIntensity() / 100.0;` and then, for every MS1 peak, `pk.setIntensity(pk.getIntensity() / max);` with no check on `max`. Three degenerate paths follow, all executed on the Release build:

- **all-zero intensities:** `max` is 0, the division is `0.0 / 0.0`, and the tool exits 0 having written **NaN into every MS1 peak**. The MS2 spectra are untouched, so the file is half NaN and half data and nothing in the output says so.
- **all-negative intensities:** `max` is negative, so every MS1 intensity **changes sign** and is scaled by 100. `[-10, -200, -3000]` is written back as `[1000, 20000, 300000]`. Exit 0.
- **an empty combined range:** `getMaxIntensity()` reads an empty `RangeBase`. Here the C++ behaves correctly and this entry does **not** extend to it: it throws `InvalidRange` and writes nothing (`exit 8`, "Empty or uninitialized range object. Did you forget to call updateRanges()?").

Severity is low — these are degenerate inputs — but the first two are silent-wrong-answer paths, not crashes.

**Proposed C++ fix:** Refuse a non-positive or non-finite scale with the tool's own error, as the empty-range path already refuses. Normalising by a non-positive maximum has no defined meaning.

**Evidence:** Executed on the Release build on dax, all four degenerate inputs, with the exit status, whether an output file was written, the decoded intensity arrays of every output and the `FileInfo` combined ranges of every input recorded: `../oracle/map-normalizer-divergence/degenerate.sh` and `degenerate.log`, hashed in `tests/data/topp_map_normalizer_provenance.json` and registered in `SOURCE_PROVENANCE.json`. The inputs come from a deterministic generator (`make_degenerate.py`) with no C++ involved.

**Rust handling:** The port refuses a non-positive or underflowed scale with exit 6 and reproduces the C++'s refusal on the empty range (that refusal was itself **executed** against the C++ after an earlier pass reasoned, wrongly, that the source treated it as a no-op). One residual difference is deliberate and recorded rather than hidden: the port reaches exit 6 through `Error::InvalidRange`, where `TOPPBase` reaches 8, so a pipeline that branches on TOPP exit codes sees a different code for the same refusal. Not an issue in this entry, and not a defect of the C++: **the combined-maximum semantics itself is correct and intended**. `updateRanges()` including chromatograms is self-consistent, and `FileInfo` prints the combined, spectrum, per-MS-level and chromatogram ranges separately, so the source knows exactly what it is doing. Whether normalising MS1 peaks against a chromatogram point is scientifically right is a question for the tool's maintainer; the port must match it either way, and now does — all 87,492 intensity arrays of the 1.2 GB benchmark run agree bitwise with the C++ tool's.

## CPP-309 — DTAFile::load and DTAFile::store use different proton masses

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source-reviewed; arithmetic verified against the cited lines. Not separately executed as a C++ store-then-load round trip.

**Affected file/function:** `src/openms/include/OpenMS/FORMAT/DTAFile.h:120` (`load`) and `:204` (`store`).

**Trigger:** Any DTA store followed by a load, or the reverse, with a precursor charge greater than 1.

**Issue:** The two directions do not use the same constant. `load` converts the singly protonated mass to m/z with

```
precursor.setMZ((mh_mass - Constants::PROTON_MASS_U) / charge + Constants::PROTON_MASS_U);
```

using `Constants::PROTON_MASS_U` = 1.00727646677, while `store` converts m/z back with the literal `1.0`:

```
os << ((precursor.getMZ() - 1.0) * precursor.getCharge() + 1.0);
```

A store-then-load round trip therefore shifts the precursor by `(charge - 1) x 7.276` mDa: 7.3 mDa at charge 2, 14.6 mDa at charge 3, and so on. At charge 1 and charge 0 the two agree, which is why the asymmetry is easy to miss. This is inside the mass-accuracy window of any modern instrument and is not a rounding artefact.

**Proposed C++ fix:** Use `Constants::PROTON_MASS_U` in `store` as well. If the `1.0` is kept for backward compatibility with files other tools have written, say so at both call sites and document that the pair is not a round trip.

**Evidence:** The two expressions above, at the cited lines of the pinned header. The port already names the asymmetry (`MassConvention::LegacyOpenMS` reproduces `store`'s `1.0`; `MassConvention::ExactProton` is the exact inverse of `load`), and its DTAExtractor output is byte-equal to the Release tool's over 36,443 files, so the port executes the source's arithmetic — but a C++-only round trip at charge > 1 was not run as a separate reproduction, and this entry does not claim it was.

**Rust handling:** `src/format/dta.rs`. The default `MassConvention::ExactProton` is the exact inverse of the reader; `MassConvention::LegacyOpenMS` reproduces the source and is what the `DTAExtractor` tool passes, so the tool's output matches the C++'s. Documented in the module header rather than silently corrected.

## CPP-310 — DTAFile::store writes two different 15-digit rules on one line, both past their type's precision

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the C++ Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/include/OpenMS/FORMAT/DTAFile.h:184` (`os.precision(writtenDigits<double>(0.0))`) and `:215` (`os << it->getPosition() << " " << it->getIntensity() << "\n"`), by way of `DPosition::operator<<` (`src/openms/include/OpenMS/DATASTRUCTURES/DPosition.h:412-420`), `precisionWrapper` (`src/openms/include/OpenMS/CONCEPT/PrecisionWrapper.h:75-80`), `StringUtils::toStr` and `Internal::NumericFormatting::appendNumeric` (`src/common/include/OpenMS/CONCEPT/Detail/NumericFormatting.h`).

**Trigger:** Every peak line of every `.dta` file the C++ writes.

**Issue:** One line carries two numbers formatted by two different rules, and both exceed what their type can represent.

- The **m/z** is a `DPosition<1>`, whose `operator<<` goes through `precisionWrapper` to `StringUtils::toStr(double, true)` and so to `std::to_chars(..., std::chars_format::fixed, 15)`. That is 15 digits *after the decimal point*. For a typical fragment m/z of a few hundred that is 18 to 19 **significant** digits — two or three past the 17 a `double` can carry, so the last digits are an artefact of the binary representation: `350.133800546540385`.
- The **intensity** is a bare `float`, promoted to `double` and written through the stream's default float field at `os.precision(15)`, i.e. 15 **significant** digits. A `float` needs 9 to round-trip; the remaining six are noise: `583.898498535156`.

The inconsistency is accidental rather than chosen — it follows from `DPosition::operator<<` taking the `precisionWrapper` route while a bare `float` takes the stream's. The consequence is a file roughly twice the size any value in it justifies. Not a correctness defect: nothing is lost, and both forms read back correctly.

**Proposed C++ fix:** One rule for both, at each type's `max_digits10` (17 for `double`, 9 for `float`), or the shortest round-tripping form for both. Either choice shrinks the file substantially and neither loses a bit.

**Evidence:** Executed. The C++ Release `DTAExtractor` over the 1.2 GB Velos benchmark run writes 36,443 files totalling exactly **836,505,793 bytes**, and this port — after it was changed to reproduce both rules exactly — writes the same 36,443 files and the same 836,505,793 bytes, byte for byte, at 1 and at 32 threads (`docs/BENCHMARKS.md` §3.6). Before the change the port used Rust's shortest round-trip text for both numbers, about 8 significant digits for an `f32` intensity, and wrote 658,901,890 bytes — **21.2 % less** for the same values, which is the size of what the source's two rules add.

**Rust handling:** Reproduced exactly, in `src/format/dta.rs`, from the two formatters `src/format/file_info/text_format.rs` already ports (`to_str` for the fixed-15-fraction rule, `ostream_g` for the 15-significant rule). The module header states which number takes which rule and what each recovers on a round trip. The cost is recorded rather than hidden: matching the source's text made `DTAExtractor` 24.5 % slower on the benchmark run, and it is the one tool in wave 4 that got slower.

## CPP-311 — MzMLSplitter writes parts whose precursor spectrumRef does not resolve

**Source revision:** TOPP `174b576e244e100f2345ca57a8e79aaa607156df`, executed on the Release build `openms4-release-bc9cc12-c19e494-174b576`. Source text read from the local TOPP checkouts `d0234cc` and `6f8eb94`, which agree.

**Status:** Executed.

**Affected file/function:** `src/MzMLSplitter.cpp`, `TOPPMzMLSplitter::main_`; the option that would have prevented it is commented out at `:67-68` (`// @TODO: // registerFlag_("precursor", "Make sure precursor spectra end up in the same part as their fragment spectra")`).

**Trigger:** Splitting any file in which an MS2 spectrum's precursor names an MS1 spectrum that lands in a different part — i.e. essentially every real DDA run split into more than one part.

**Issue:** mzML 1.1's `xs:keyref KEYREF_PRECURSOR_SPECTRUMREF` (`share/OpenMS/SCHEMAS/mzML_idx_1_10.xsd`) requires `precursor/@spectrumRef` to resolve to a spectrum id **in the same document**, and the schema text says the attribute is "for precursor spectra that are local to this document". The splitter cuts the spectrum list at a byte or count boundary and copies each spectrum's precursor unchanged, so systematically produces parts that violate that keyref. The schema provides `sourceFileRef` + `externalSpectrumID` for exactly this case and the tool uses neither.

**Proposed C++ fix:** Rewrite a crossing reference as an external one (`sourceFileRef` + `externalSpectrumID`), or drop it, or finish the `@TODO` flag above — and, whichever is chosen, document the behaviour.

**Evidence:** Executed on the Release build over the 1.2 GB Velos benchmark run split into four parts: part 2 of 4 carries 4 references to `scan=10929`, which is in part 1. The part files are the run's own retained outputs.

**Rust handling:** The port's mzML writer used to refuse a precursor reference that does not name a spectrum in the file being written, which is why `MzMLSplitter` could not process real input at all in wave 3 ("precursor spectrum reference does not name an output spectrum" on every repetition). Matching the executed C++ was the right call for a port, so the writer now emits the dangling reference as the source does, under a tool-side option — the library default stays strict and `FileHandler` does not enable it (lead decision of 2026-09-15, decision D10). The rule and the schema file are described in `docs/MZML_SUPPORT.md`. With that, all four parts are bitwise equal to the C++ tool's on every decoded array.

## CPP-312 — FeatureFinderAlgorithmPicked divides a zero-width retention-time range by intensity:bins

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on both the Debug product SDK and the Release build.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:244-245` in `run_`, with the consequence at `:1837-1838` in `intensityScore_`.

**Trigger:** An input with more than `2 * min_spectra_` MS1 scans whose spectra share one retention time or one m/z; separately, any input with fewer than `2 * min_spectra_` scans.

**Issue:** The intensity-score grid is built with

```
intensity_rt_step_ = (…getMaxRT() - rt_start) / (double)intensity_bins_;
intensity_mz_step_ = (…getMaxMZ() - mz_start) / (double)intensity_bins_;
```

and neither numerator is checked for being zero. With `maxRT == minRT` the step is `0.0`, and `intensityScore_` then evaluates `std::floor((rt - rt_min) / intensity_rt_step_ * 2.0)` = `floor(0.0/0.0)` = `floor(NaN)`, whose conversion to `UInt` is undefined behaviour before `std::min` ever sees it. The **Release** build computes non-finite bin bounds and **silently returns an empty feature map with exit 0**, which a caller cannot distinguish from a run without features.

`FileFilter_44_input.mzML` is not an instance of the zero-width case: it holds two MS1 spectra, and with the default `mass_trace:min_spectra 10` (`min_spectra_` 5) the seed loop `min_spectra_ .. n - min(min_spectra_, n)` (`.cpp:297`, `:493`) is empty for every input of at most 10 scans, so its empty map is fixed by its length. Its Debug exit 8 is not the zero-width range either: steps 2 and 3.2 call `startProgress(min_spectra_, n - min(min_spectra_, n))` (`.cpp:297-298`, `:493-494`), e.g. `startProgress(5, 0)`, and `ProgressLogger::startProgress` checks `begin <= end` only with `OPENMS_ASSERTIONS` (`OPENMS_PRECONDITION` at `ProgressLogger.cpp:235`), so a Debug build throws `Exception::Precondition` and exits 8 on every input shorter than `2 * min_spectra_` scans, while a Release build stores the inverted range and returns the empty map; the two build types disagree on the same valid input, and the documented contract ("Sets the progress range from begin to end") does not say which is intended.

**Proposed C++ fix:** Refuse a zero-width range in either dimension with a stated error, or collapse to a single bin when the range is zero (`intensity_bins_ = 1`) and say so in the log. For the inverted progress range: skip or clamp the two progress sections when `begin > end`, or make `startProgress`'s check unconditional and document it. Debug and Release should agree.

**Evidence:** Tool cases `zero_rt`, `zero_mz` (1 and 4 threads, `seed:min_score 0`), `zero_mz_control_min_score_0`, `filefilter_44_force` and `fileconverter_31` executed on the Release `FeatureFinderCentroided` (three repetitions). The **Debug** exit 8 is `FFC_FileFilter_44_force` (`debug_only`, `tests/data/topp_feature_finder_centroided_provenance.json`) - the short-input case, not the zero-width one. The Release range behaviour is `../oracle/progress-logger-release-range/driver.cpp` (`tests/data/progress_logger_release_range.tsv`, 60 calls, 3 identical runs); the FFAP Release event sequence of the short case is `S 5 0` (`../oracle/ffap-instr-completion`, progress `short`). The Release run is the reference behaviour, as always in this log — never the Debug exit code.

**Rust handling:** The tool reproduces the Release outcome (exit 0, the source's lines, an empty map); `ProgressLogger` accepts inverted ranges as the Release build does and FeatureFinderAlgorithmPicked passes them unchanged. The wave-5 port therefore no longer diverges here: `a_short_input_never_reaches_the_seed_loop_as_in_the_cpp_release_build` (no longer `#[ignore]`d), `a_zero_width_retention_time_range_follows_the_cpp_release_build`, `a_zero_width_mz_range_follows_the_cpp_release_build` and `the_progress_event_sequence_matches_the_release_build` pin it. See also CPP-274, the same conversion reached from a degenerate range rather than from a short input.

## CPP-313 — FeatureXMLHandler caps its feature reservation at 1e5 on a premise current data exceeds

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Release build.

**Status:** Executed (the premise), source-reviewed (the effect).

**Affected file/function:** `src/openms/source/FORMAT/HANDLERS/FeatureXMLHandler.cpp:318`, `startElement`.

**Trigger:** Reading any featureXML that declares more than 100,000 features.

**Issue:**

```
map_->reserve(std::min(Size(1e5), count)); // reserve vector for faster push_back, but with upper boundary of 1e5 (as >1e5 is most likely an invalid feature count)
```

The comment's premise — that a declared count above 1e5 is "most likely an invalid feature count" — is wrong by an order of magnitude for current data. The benchmark's own `MassTraceExtractor` output declares **826,019** features and the executed Release build reads it correctly. The effect is only a lost reservation, so the vector grows repeatedly instead of once; it is a performance and comment defect, not a correctness one.

**Proposed C++ fix:** Reserve `count` outright, or raise the cap to something a current instrument run can actually exceed, and remove the claim in the comment. If a cap is wanted as a guard against a hostile declared count, say that instead — the current comment says the opposite.

**Evidence:** Executed: the C++ Release `FileInfo` at core `bc9cc12` reads the 826,019-feature map and the 2.06 GiB benchmark map without complaint, so the premise is false on real data. The cost of the missed reservation was not measured and is not claimed.

**Rust handling:** The port's own ceiling on this path was a different and worse defect — a fixed ~12.5 MB decode limit formed by three limits combined with `min()`, which refused the featureXML the port's own `FeatureFinderCentroided` writes — and it is fixed in this window: `src/format/featurexml_scaling.rs` derives the ceilings from the document's size, and features are streamed rather than retained. Both benchmark maps (59.6 MiB / 42,789 features and 2.06 GiB) now load and round-trip. The source has no ceilings at all, so this port is deliberately stricter and says where. See `docs/FEATUREXML_SCALE_SUPPORT.md`.

## CPP-314 — FeatureFinderAlgorithmPicked converts NaN, infinite and huge doubles to Size without a check

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:357` (the step-2.5 window count) and `:1200`, `getIsotopeDistribution_`.

**Trigger:** An infinite, NaN or very large m/z anywhere in the input, reached through `run`.

**Issue:** Both sites convert a `double` to `Size` with no range check: `Size num_isotopes = std::ceil(max_mass / mass_window_width_) + 1` and `Size index = (Size) std::floor(mass / mass_window_width_)`. The Release build emits `comisd` against `2^63`, then `cvttsd2si`, then `btc` (`libOpenMS.so` `0x18e46f4`, `0x18dddd0`): `+inf` and values of `2^64` and above give 0 windows, so an infinite or `1e300` m/z ends in "the value '12' was used but is not valid; IsotopeDistribution not precalculated. Maximum allowed index is 0"; a NaN m/z asks for window `2^63` ("9223372036854775808"); a count above `vector::max_size()` = 164,703,072,086,692,425 makes `resize` throw `std::length_error` ("vector::_M_default_append"; executed at `2^62 - 512`, `2^62`, 164,703,072,086,692,448 and `1.5*2^63`), and a count at or below it allocates 56 bytes per window (164,703,072,086,692,416 windows: `std::bad_alloc` on the reference node).

**Proposed C++ fix:** Reject non-finite m/z values in `run`, and check the count before converting.

**Evidence:** `nonfinite_stage` cases `mz_posinf_last`, `mz_posinf_spectrum0`, `mz_huge_1e300`, `mz_huge_2p63`, `mz_nan_*`; `sort_mobility_stage` cases `v2_mz_2p64`, `v2_mz_below_2p64`, `v3_mz_count_above_max_size`, `v3_mz_count_below_max_size` (`../oracle/ffap-sem-completion`, `../oracle/ffap-complete-fix2`).

**Rust handling:** `x86_64::truncate_to_u64` reproduces the conversion; above `max_size()` the port returns the `length_error` text (`IsotopeWindows::precalculate_onto`, `LENGTH_ERROR_WHAT`); at or below it `Limits::max_isotope_windows` refuses counts above its ceiling, since the source's outcome there depends on memory (lead decision D6 of wave 5). A debug run has opened `debug/log.txt` and created `debug/features` before either failure and keeps both (executed: `../oracle/ffap-complete-fix2` `lenerr_single` and `lenerr_reuse` at m/z `1e19` and `2e18`). Because `std::length_error` and `std::bad_alloc` are no OpenMS exceptions, `FeatureFinderCentroided` reports them from `TOPPBase`'s outer `std::exception` handler ("Unable to initialize or run FeatureFinderCentroided: vector::_M_default_append" / "std::bad_alloc", exit 12; executed `tool_1e19` and `tool_2e18`); the port's tool does the same for the `length_error` and exits 8 with its own message at its window ceiling.

## CPP-315 — FeatureFinderDefs is defined in two headers, so a translation unit including both does not compile

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Compiled with the Release install's GCC 14.

**Status:** Executed (compile).

**Affected file/function:** `src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h:24-55` and `src/openms/include/OpenMS/FEATUREFINDER/FeatureFinderDefs.h:19-50`.

**Trigger:** Any translation unit that includes both headers.

**Issue:** `FeatureFinderDefs` (with `NoSuccessor`, `IndexPair`, `IndexSet` and `ChargedIndexSet`) is defined in full in both headers, with no include guard shared between them, so including both is a redefinition and the translation unit does not compile.

**Proposed C++ fix:** Delete one definition and include `FeatureFinderDefs.h` from the algorithm header.

**Evidence:** `../oracle/ffap-sem-completion/drivers/defs_both_headers.cpp`: `-fsyntax-only` fails with the Release install's GCC 14.

**Rust handling:** One definition in `src/analysis/feature_finder_picked/defs.rs`, whose `ChargedIndexSet` compares its index sets only, as the source's inherited `std::set` operators do (executed: `../oracle/ffap-complete-fix1/drivers/defs_eq_probe.cpp`).

## CPP-316 — A reused FeatureFinderAlgorithmPicked instance carries four kinds of state into its next run

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed for (a), (c), (d) and (e); (b) is source-reviewed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:1138` (`abort_reasons_`), `:361` and `:378` (`isotope_distributions_`), `:231` (`log_`), `:137`, `:844` and `:1062` (the caller's map).

**Trigger:** Calling `run` twice on one `FeatureFinderAlgorithmPicked` object, or passing a non-empty map.

**Issue:** Four members survive a run and are never cleared.

(a) `abort_reasons_` is never cleared; its only writer is `:1138`. A reused object's abort map reads earlier runs' spectrum and peak indices from the **new** input. Executed: wrong peaks on an input of the same shape, and SIGSEGV on a shorter one (9 of 9), which leaves `debug/log.txt` at the first run's flushed prefix (1,114,578 of 1,118,230 bytes) and the second run's seed maps.

(b) `abort_reasons_` is a `std::map<Seed, std::string>` and `Seed::operator<` compares intensity only, so equally intense seeds collapse into one entry (source-reviewed).

(c) `isotope_distributions_` is never cleared. Step 2.5 resizes the kept windows and appends to them (`:361`, `:378`), so a reused object's second run finds different features on the same input (executed).

(d) `log_` is opened on every debug run and never closed (`:231`). A second debug run of the same object fails to open it and writes no log (executed). The first run's unflushed tail (up to 8,191 bytes) reaches the file only when the object is destroyed, so a process that dies in any later run, with or without `write_debug`, leaves the first run's flushed prefix (executed: 1,163,782 of 1,165,129 bytes after SIGSEGV in the seed loop or at the score arrays, 139,387 of 141,951 after SIGFPE in step 4; the complete log after a return).

(e) `run()` never clears a non-empty caller map (`:137`, `:844`, `:1062`). Old features are resolved against new ones, re-sorted, and re-annotated with the new input's `spectrum_index` and native id (executed).

**Proposed C++ fix:** Clear all four at the start of `run`, and document whether a non-empty output map is appended to or replaced. For (b), give `Seed` a total order.

**Evidence:** `../oracle/ffap-instr-completion` (reuse of one object over three runs, `stale scaled` and `stale oob`), `../oracle/ffap-complete-fix4` `fix4_reuse` (four reuse scenarios) and `../oracle/ffap-complete-fix5` reuse cases; every case run twice and identical but for the abort map's unique id.

**Rust handling:** (a) is refused at that read (`debug::abort_map`) with a `DebugTermination` (`TerminationPoint::AbortMap`, `TerminationKind::OutOfBounds`) whose `log_file_bytes` is that prefix. (c), (d) and (e) are reproduced: the instance keeps its isotope windows, its never-closed stream's counts (`FeatureFinderAlgorithmPicked::debug_log_file`, and every termination reports that length through `DebugTermination::log_file_bytes`) and the caller's features. (b) follows from the ported `Seed::is_less_intense_than`.

## CPP-317 — FeatureFinderAlgorithmPicked's step 4 computes a charge remainder without a zero check

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:936` and `:945`, the overlap resolution in `run`.

**Trigger:** A caller's feature with charge 0 (the featureXML default) in a mixed-charge overlap, or the pair `INT_MIN` and `-1`.

**Issue:** `f2.getCharge() % f1.getCharge()` and `f1.getCharge() % f2.getCharge()` are evaluated with no zero check and no `INT_MIN`/`-1` check. A charge-0 feature traps (executed: SIGFPE in every repetition), and so does `INT_MIN % -1` (SIGFPE, 2 of 2). A debug run leaves `debug/log.txt` at the flushed prefix (1,163,782 bytes, the pair's `Intersection` line still buffered), the seed map and the feature files, and no abort map or input.

**Proposed C++ fix:** Skip the pair when either charge is 0, and guard the `INT_MIN`/`-1` case.

**Evidence:** `../oracle/ffap-instr-completion` (runs into caller maps, 11 prefilled and 5 overlapping features) and `../oracle/ffap-complete-fix5`; the `INT_MIN % -1` case is the round-1 verifier's, SIGFPE 2 of 2.

**Rust handling:** Refused at that pair (`resolution::resolve_overlaps`) with a `DebugTermination` (`TerminationPoint::OverlapResolution`, `TerminationKind::ArithmeticTrap`) for runs with and without `write_debug`.

## CPP-318 — DefaultParamHandler::setParameters assigns param_ before checkDefaults throws

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/DATASTRUCTURES/DefaultParamHandler.cpp`, `setParameters`, reached from every `DefaultParamHandler` subclass.

**Trigger:** Any `setParameters` call whose value fails a restriction.

**Issue:** `param_` is assigned before `checkDefaults` throws, so after a refused call `getParameters()` returns the **refused** set while the typed members keep their old values. A caller that catches the exception and reads the parameters back sees a state the object never used.

**Proposed C++ fix:** Validate into a local `Param` and assign only after `checkDefaults` returns.

**Evidence:** Executed through FeatureFinderAlgorithmPicked (`../oracle/ffap-instr-completion`, the parameter surface including a refused set).

**Rust handling:** The port's `set_parameters` reports the error and leaves the typed settings as the source leaves them; the observable state after a refused call is pinned for FeatureFinderAlgorithmPicked.

## CPP-319 — ParamValue::operator double() returns the union's double member for a string or list value

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/DATASTRUCTURES/ParamValue.cpp:397-408`, `operator double()`; reached from `FeatureFinderAlgorithmPicked.cpp`, `writeFeatureDebugInfo_`.

**Trigger:** Reading a `STRING_VALUE` or any list-valued parameter as a `double`, for example `advanced:pseudo_rt_shift` given as a string.

**Issue:** `operator double()` throws only for `EMPTY_VALUE` and converts `INT_VALUE`; for `STRING_VALUE` and the three list types it falls through to `return data_.dou_;`. For a `STRING_VALUE` the live union member is a `std::string*`, so this reinterprets a heap pointer's bit pattern as a `double` (`movsd 0x8(%rdi)` in `libOpenMS.so`). `writeFeatureDebugInfo_` then writes address-dependent numbers near RT 0 into the debug `.dta` files, so one run's debug output cannot be compared with another's. This is the same defect CPP-058 and CPP-141 record for `DataValue`'s floating conversions, in the other value class; `ParamValue` is a separate file with its own copy of it, and this is the first entry with an executed, address-dependent consequence.

**Proposed C++ fix:** Give `operator double()`, `operator float()` and `operator long double()` the type guard their integer siblings already have.

**Evidence:** `../oracle/ffap-instr-completion` `ffap_shift_band_driver`: 3 processes wrote 3 different `0.dta` files at RT 0, `1e-295` and `1e-300`, and identical files at `5e-275` and `1e-289` (75 files each for a string shift, a list shift, and a string shift with a scan at RT `5e-275` and at `1e-289`).

**Rust handling:** The port's debug writer reproduces every shift the value of which is reproducible, and the address-dependent band is refused, not invented: `RejectedParameters::Shown` is the default and the heap-address bound is scoped to the reference platform (lead decision D8 of wave 5). The 304 non-finite shift files outside that band are byte for byte.

## CPP-320 — The debug comment at FeatureFinderAlgorithmPicked.cpp:1047 names the wrong erased array

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed (the effect), source-reviewed (the comment).

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:1047`, `run`.

**Trigger:** `write_debug` enabled.

**Issue:** The comment reads "store input map with calculated scores (without overall score)", but the loop erases float data array 2, which is `local_max` (`.cpp:202`), not an overall score. The written `debug/input.mzML` therefore keeps every overall-score array and loses `local_max`. A reader who trusts the comment reads the wrong array names. Documentation defect; no numeric consequence.

**Proposed C++ fix:** Erase the intended array, or fix the comment to say `local_max`.

**Evidence:** The executed `debug/input.mzML` array names (`../oracle/ffap-instr-completion`, tool cases `a1`-`a3`, byte-identical to the port's).

**Rust handling:** The port writes the same arrays as the executed build, so `debug/input.mzML` is byte for byte; the port's own comment names `local_max`.

## CPP-321 — FeatureFinderAlgorithmPicked's step 1 silently drops scans with a non-finite drift time

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:261` with `MSExperiment::areaBeginConst` (`MSExperiment.cpp:562-571`) and `AreaIterator::nextScan_` (`AreaIterator.h:277-298`).

**Trigger:** A scan whose drift time is NaN or `+-inf`.

**Issue:** `areaBeginConst` sets the mobility range to `RangeMobility{}.getNonEmptyRange()` = `[lowest, max]`, and `nextScan_` skips every scan whose drift time the range does not contain - which is every non-finite one. Such a scan's intensities are left out of **every** intensity quantile, while its peaks are still scored against those quantiles and can become seeds. Nothing is reported.

**Proposed C++ fix:** Iterate without a mobility filter in step 1, or refuse non-finite drift times in `run`.

**Evidence:** `../oracle/ffap-complete-fix1` (`nonfinite_stage_dt`, two runs each, identical): a NaN, `+inf` or `-inf` drift time on scan 50 removes its intensities from the FFC_1 quantiles (8 features); sixteen NaN drift times change them further; every drift time NaN leaves every cell empty and finds 10 features instead of 8; `f64::MAX`, `f64::MIN` and finite drift times keep the scan (`v2_dt_*`, `v3_dt_*`).

**Rust handling:** `IntensityThresholds::compute` applies the same filter (`RangeBase::contains` on the full range); the mzML reader refuses non-finite drift times, so only library callers reach it.

## CPP-322 — FeatureFinderAlgorithmPicked dereferences an empty best isotope pattern when feature:min_isotope_fit is 0

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed for the empty-pattern case; the `size_t(-1)` sub-case has never been observed and its executed outcome is unknown.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:614-625` and `extendMassTraces_` (`:1347-1349`).

**Trigger:** `feature:min_isotope_fit` 0 - inside its valid range `[0, 1]` - and a seed for which `findBestIsotopeFit_` finds no placement.

**Issue:** `findBestIsotopeFit_` returns 0 and leaves `best_pattern` empty; `0 < min_isotope_fit_` is false, so the seed is not aborted, and `extendMassTraces_` reads `pattern.spectrum[0]` and `map_[...][pattern.peak[0]]` of empty vectors: an out-of-bounds read. The same read with a non-empty pattern whose first matched isotope has no peak reads `map_[spectrum][size_t(-1)]`, heap metadata just before a spectrum's peaks; that sub-case was never observed (6,076 refusals of a 19,200-run grid were all empty patterns) and its executed outcome is unknown.

**Proposed C++ fix:** Abort the seed when the pattern is empty, or give the parameter a positive minimum.

**Evidence:** `../oracle/ffap-complete-fix3` (two runs each): the stage cases `g_avg_trace0`, `g_iso0_seed0`, `p_ipo_100_seed0_iso0` and `p_ipo_nan_seed0_iso0` die with SIGSEGV (a gdb backtrace in `../oracle/ffc-numerics-v2` shows the fault in `extendMassTraces_` under `run_`), `g_avg_trace0_iso_tiny` (bound `1e-300`) returns 14 features; the negative-intensity library runs `neg_oob1`, `neg_oob_seed035` and `neg_none_avg0` die with SIGSEGV with and without `write_debug` and leave a debug log truncated at the file buffer; `FeatureFinderCentroided` case `avg0` dies with SIGSEGV, having written its console lines only up to the FAIMS line.

**Rust handling:** Refused at that read (`extension::EMPTY_PATTERN_WHAT`), both sub-cases; a debug run keeps the seed's log lines, and every run (since round 5 also without `write_debug`) records an `OutOfBounds` termination (its SIGSEGV label is established for the empty pattern), and the tool writes only the flushed log and exits 8.

## CPP-323 — Heavy averagine windows are silently empty

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:365-424` with `CoarseIsotopePatternGenerator::estimateFromPeptideWeight` and `IsotopeDistribution::renormalize`/`trimLeft`/`trimRight`.

**Trigger:** An averagine window centre from about 273,770 Da on (executed between 273,769.5 Da, which keeps one bin, and 273,770.5 Da, the first empty window at width 1; at width 100 window 2738, centre 273,850 Da, is the first empty one and windows up to 2737 keep a single bin).

**Issue:** All 20 binary32 bins of the estimate underflow to zero. `renormalize` divides zero by the zero sum, so every weight is NaN; `trimLeft` erases nothing and `trimRight` everything, so the window is empty with maximum 0 and every peak of such a mass scores 0, without a message.

**Proposed C++ fix:** Compute the averagine in `double` or in log space, or report the mass limit.

**Evidence:** `../oracle/ffap-complete-fix3`, `u_*` cases (m/z 136,850.5 to `1e6` at charge 2, charge 1000 on the first 20 scans; two runs each): the runs return (8, 13 and 0 features), with m/z 136,850 (no such window) and 100,000 as controls; `../oracle/ffap-complete-fix4` `vw_*` cases (isotope windows printed at widths 1, 7.3, 100 and 200 across the boundary).

**Rust handling:** The crate-private `CoarseIsotopePatternGenerator::estimate_from_peptide_weight_source` reports the case, and step 2.5 empties the window as the source does; the public generator functions keep their error.

## CPP-324 — A NaN isotopic_pattern:intensity_percentage_optional is accepted and empties every isotope window

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/DATASTRUCTURES/Param.cpp`, `ParamEntry::isValid` (both range comparisons are false for NaN), and `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:372-374`.

**Trigger:** `isotopic_pattern:intensity_percentage_optional` set to NaN; the other `double` parameters accept NaN the same way.

**Issue:** `trimLeft(NaN)` erases nothing and `trimRight(NaN)` erases everything, so every isotope window is empty and the run finds no seed and no feature, silently.

**Proposed C++ fix:** Reject NaN in `isValid`.

**Evidence:** `../oracle/ffap-complete-fix3` `p_ipo_nan`, `p_ipo_nan_min0`, `p_ipo_nan_seeds`, `p_ipo_nan_seeds_min0`, `p_ipo_nan_bins3`, `p_ipo_nan_egh` (0 features each, two runs).

**Rust handling:** Step 2.5 applies the source's comparisons for a NaN cutoff (empty windows); the shared `IsotopeDistribution` trims keep refusing a NaN cutoff for their other callers.

## CPP-325 — An unsorted input with a mis-sized data array leaves the caller's map partly sorted

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:1081-1086`, `run`; `MSExperiment::sortSpectra` and `sortChromatograms` (`MSExperiment.cpp:791-822`); `MSSpectrum::sort` and `MSChromatogram::sort` (`checkDataArraySizes_`).

**Trigger:** An unsorted input in which some spectrum's non-empty data array differs in length from its peaks.

**Issue:** `run` takes the map by rvalue reference and sorts it in place. `sortSpectra` reorders all spectra and sorts the peaks of the earlier unsorted ones before `MSSpectrum::sort` throws `Exception::Precondition` for the first such spectrum **in retention-time order** (then chromatograms), so the caller's object is left half-sorted; and which spectrum is reported depends on `std::sort`'s order of equal retention times.

**Proposed C++ fix:** Check every data array before sorting anything.

**Evidence:** `../oracle/ffap-complete-fix3` `a_*` cases (two runs each): "FloatDataArray[0] size (25) does not match spectrum size (24)" for the spectrum that sorts first, not the first in input order; string and integer variants; an empty array 0 before a mis-sized array 1; "does not match chromatogram size (5)"; sorted spectra and exact arrays run through.

**Rust handling:** The same order and text (`Error::InvalidValue`); `run` consumes the experiment, so the partial sort is not observable.

## CPP-326 — FeatureFinderAlgorithmPicked's step-1 progress range wraps modulo 2^32

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:241`, `startProgress(0, intensity_bins_ * intensity_bins_, ...)`.

**Trigger:** `intensity:bins` at or above 65,536.

**Issue:** The product of two `UInt` values wraps modulo `2^32`, so 65,536 bins announce the range `[0, 0]` and 100,000 bins `[0, 1410065408]`, while `setProgress` then reports `Size` values up to `bins^2`. Cosmetic: the progress display is wrong, the computation is not.

**Proposed C++ fix:** Multiply in `Size`.

**Evidence:** `../oracle/ffap-complete-fix3` `fix3_driver progress`, nine `bins` values including `2^32 + 65,536` (narrowed to 65,536), two runs each.

**Rust handling:** The same wrapped start event (`seeds::step_one_progress`, crate-private). This is an in-bounds wrap, so lead decision D12 of wave 5 has the port reproduce it rather than refuse it.

## CPP-327 — Large retention times silently give features of infinite width and intensity

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:742` (`f.setWidth(fitter->getFWHM())`) and `:790` (`f.setIntensity(fitter->getArea() / ...max)`), which narrow `double` to `float`.

**Trigger:** Finite but large retention times. On FeatureFinderCentroided_1 with every retention time scaled, the first infinite width appears at a scale of `6e36` (1 of 9 features; 3 of 9 at `8e36`, 7 of 9 at `1e37`, all from `1.5e37`; none at `4e36`), and every intensity is infinite from `1e36` on (finite at `1e33`).

**Issue:** The fitted sigma passes about `1.44e38` and the `float` FWHM overflows. The run returns features with an infinite width, FWHM meta value and intensity, and `FeatureFinderCentroided` writes `inf` into the featureXML, without a message.

**Proposed C++ fix:** Check the fitted width and area against the `float` range, or store them as `double`.

**Evidence:** `../oracle/ffap-complete-fix4` (two runs each, identical): stage cases `vy_rt_1e33` (finite), `vy_rt_1e36` (finite widths, infinite intensities), `vy_rt_1e37` (7 of 9 widths infinite), `vy_rt_1e38` to `vy_rt_1e150` (all), Gaussian and EGH, `vx_rt_1e150` to `vx_rt_1e300`; `../oracle/ffap-complete-fix5` `run_onset.sh` (`nb_rt_2e36` to `nb_rt_5e37` and three jittered `1e37` inputs); the Release `FeatureFinderCentroided` on FFC_1 with every scan start time written as `ve36` or `ve39`: exit 0, featureXML with `<intensity>inf</intensity>` and FWHM `inf`.

**Rust handling:** FeatureFinderAlgorithmPicked stores the same values (the width field, the crate-private `MetaValue::source_float`), bit for bit; `BaseFeature::validate` and the featureXML writer refuse such features, so the port's `FeatureFinderCentroided` exits 3 without an output file (TOPP native difference 16, `infinite_feature_values_are_refused_by_the_featurexml_writer`). Whether the writer should follow the C++ build and write `inf` is a separate, recorded task.

## CPP-328 — A trace of zero intensities terminates FeatureFinderAlgorithmPicked

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `intensityScore_` (`src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmPicked.cpp:1913-1937`), `extendMassTrace_` (`:1535`), `MassTrace::getAvgMZ` and `getIsotopeDistribution_` at `:790`, inside step 3.3's OpenMP region.

**Trigger:** Peaks of intensity zero in intensity cells whose quantiles are zero (for example one m/z band of an input zeroed), `feature:reported_mz` `maximum` or `monoisotopic`, and thresholds that let such a feature through.

**Issue:** The zero peak's intensity score is `0 / 0 = NaN`, which passes the `< 0.01` test, so the zero peaks form a trace. Its average m/z is `0 / 0`, and `getIsotopeDistribution_(NaN)` converts the NaN to index `2^63` and throws `Exception::InvalidValue` ("the value '9223372036854775808' was used but is not valid; IsotopeDistribution not precalculated. Maximum allowed index is 15"), which leaves the OpenMP region, so `std::terminate` kills the process (SIGABRT, OpenMS's fatal-exception block on stdout). With `reported_mz average` the other traces keep the sum finite and the run returns.

**Proposed C++ fix:** Score a zero intensity as zero, skip traces without intensity, check the reported m/z, and catch exceptions inside the parallel region.

**Evidence:** `../oracle/ffap-complete-fix5` `run_band.sh`, 16 cases twice each, identical (with and without `write_debug`: 12 SIGABRT, 2 returns with `reported_mz average`, 2 SIGSEGV of the empty best pattern at `feature:min_isotope_fit` 0); the debug runs leave the flushed log prefix, the seed map and the feature files of every plot up to the terminating seed's, whose `.plot` prints the all-zero trace's m/z as `-nan`.

**Rust handling:** Refused at that seed with a `DebugTermination` (`TerminationKind::Exception`, the executed `what()` as its message) after the seed's log lines and feature files, all equal to the executed ones (`a_step_3_3_5_termination_keeps_what_the_executed_process_had_written`); `MassTrace::avg_mz` follows the Release build's SSE NaN rules, so the `-nan` is printed on every host.

## CPP-329 — SignalToNoiseEstimatorMedian leaves its percentage members uninitialised and operator= keeps stale ones

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Source-reviewed.

**Affected file/function:** `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h:434-436` and `:80-143`.

**Trigger:** Calling `getSparseWindowPercent()` or `getHistogramRightmostPercent()` before `init`, or assigning one estimator to another.

**Issue:** Neither constructor sets `sparse_window_percent_` or `histogram_oob_percent_`; `updateMembers_` does not set them and the copy constructor does not copy them, so both getters read an indeterminate `double` before `init`. `operator=` clears `stn_estimates_` but keeps the target's old percentages, so a reused estimator reports the previous input's diagnostics.

**Proposed C++ fix:** Initialise both to 0, and copy or reset them in the copy constructor and `operator=`.

**Evidence:** Source review of the pinned header. Precedent: CPP-289, the same class of defect in `GaussTraceFitter`.

**Rust handling:** The port's `NoiseEstimates` is stateless: both percentages are returned with the estimates of the run that computed them, so there is nothing to read early or to carry over.

## CPP-330 — SignalToNoiseEstimatorMedian's empty-median-bin fallback is dead code and its comment is wrong

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Source-reviewed, with a proof in `docs/SIGNAL_TO_NOISE_SUPPORT.md`.

**Affected file/function:** `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h:350-353`.

**Trigger:** Never; the branch is unreachable.

**Issue:** The `else` branch is commented "only possible if the rightmost bin was hit while empty (already flagged above)". The bin counts sum to `elements_in_window`, which is at least `(elements_in_window + 1) / 2`, so the median walk always stops in a non-empty bin and the branch cannot be reached. The comment describes a state the code cannot be in.

**Proposed C++ fix:** Delete the branch, or turn it into an assertion and fix the comment.

**Evidence:** Source review of the pinned header; the counting argument is written out in `docs/SIGNAL_TO_NOISE_SUPPORT.md`.

**Rust handling:** Not ported as a branch; the port's median walk carries the same invariant and the proof is recorded beside it.

## CPP-331 — SignalToNoiseEstimatorMedian divides zero by zero for an empty container

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h:373-374` in every mode, and `SignalToNoiseEstimator.h:127` in `AUTOMAXBYSTDEV`.

**Trigger:** `init` on an empty container.

**Issue:** `window_count` is 0, so `sparse_window_percent_ = sparse_window_percent_ * 100 / window_count` and the `histogram_oob_percent_` line are `0 / 0`; both percentages become NaN (`0xfff8000000000000`), as does `max_intensity_` in mode 0. No error is reported.

**Proposed C++ fix:** Return early for an empty container, or report 0 percent.

**Evidence:** Executed on the Release build: cases `empty_stdev`, `empty_manual`, `empty_stdev_chrom` and `progress_empty` (`../oracle/sne-completion`).

**Rust handling:** `NoiseCompatibility::nan_for_empty_input` reproduces the NaN bit pattern under `PickingCompatibility::source()`; the native profile reports the empty input instead.

## CPP-332 — estimateNoiseFromRandomScans ignores its ms_level filter and can read out of bounds

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.cpp:21-53`.

**Trigger:** Any call; the out-of-bounds read needs a percentile of exactly 100, an empty drawn scan, or a percentile outside `[0, 100]`.

**Issue:** Four defects in one function. `:44` reads `exp[scan]`, not `exp[spec_indices[scan]]`, so the drawn spectrum can have any MS level and can be empty, and `:50` then reads `tmp[0]` of an empty vector. `:42` scales by `(spec_indices.size() - 1)`, so the last candidate is never drawn. A percentile of exactly 100 reads `tmp[size()]` at `:50`. And `Size idx = tmp.size() * percentile / 100.0` at `:48` is an unchecked `double`-to-`unsigned long` conversion, after which `tmp.begin() + idx` and `tmp[idx]` both compute `_M_start + 4*idx mod 2^64` (`lea (%r12,%rcx,4),%r15`, scale 4, wraps) - undefined pointer arithmetic, in bounds or not depending on the value. The seed is `time(nullptr)`, so no run is reproducible.

**Proposed C++ fix:** Index through `spec_indices`, scale by `size()`, clamp the percentile and the index, and take the seed as a parameter.

**Evidence:** 47 in-domain cases executed on the Release build with the seed set through an interposed `time()` (`../oracle/sne-completion`); `../oracle/sne-followup` pins the in-bounds wrap in 10 cases (twice each, byte-identically); an empty drawn scan and an ordinary negative percentile read out of bounds (the `-50%` probe read mapped heap 8,192 bytes before the buffer, twice).

**Rust handling:** `estimate_noise_from_random_scans` / `RandomScanNoise` takes an explicit seed and reproduces libstdc++'s `minstd_rand0`, `generate_canonical<double,53>` and GCC 14.4.0's `nth_element` bit for bit, together with the `exp[scan]` indexing, the `(size()-1)` scaling and the `:48` conversion. The pointer wrap is reproduced where it stays in bounds (the element used is `e = idx mod 2^62`) and refused exactly when `e >= size` (`Error::Unsupported`), which is where the source reads out of bounds.

## CPP-333 — SignalToNoiseEstimatorMedian's int counters overflow beyond INT_MAX points

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The `:220` conversion is executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Source-reviewed, with the `:220` path executed.

**Affected file/function:** `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h:123` (the point count) and `SignalToNoiseEstimatorMedian.h:216` (the histogram counts), `:220` (`(int)` of `p*n/100`), `:228` (`elements_seen`), `:310`/`:324` (`elements_in_window`) and `:365` (`window_count`).

**Trigger:** A container with more than `INT_MAX` points, or a window that holds that many.

**Issue:** Every one of these counters is an `int` counting points. Past `2^31` they overflow, which is undefined for a signed type; `(int)(auto_max_percentile_ * c.size() / 100)` at `:220` is an undefined `double`-to-`int` conversion from `2^31` on, where the Release build's 32-bit `cvttsd2si` returns `INT_MIN` and the median walk is skipped entirely.

**Proposed C++ fix:** Use `size_t` or `SignedSize` for all of them.

**Evidence:** Source review of the pinned header, with `:220` executed at `n = 2^31` and `n = 3,000,000,001` (`../oracle/sne-fix`, five AUTOMAXBYPERCENT cases, each needing 64-90 GB).

**Rust handling:** `:220`'s `INT_MIN` is reproduced. The signed-overflow sites stay refused with `Error::Unsupported` naming the line, as a stated cost/benefit decision of the lead: each needs more than `2^31` points per spectrum and 64-90 GB per evidence run. `docs/SIGNAL_TO_NOISE_SUPPORT.md` records that, and the beyond-`INT_MAX` full-path evidence lives outside CI in `../oracle/sne-fix/harness`; unit tests pin the helpers at `n = 2^31`.

## CPP-334 — Param accepts NaN for a bounded double

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/DATASTRUCTURES/Param.cpp:146`, `ParamEntry::isValid`.

**Trigger:** Setting any `DOUBLE_VALUE` parameter that carries `setMinFloat`/`setMaxFloat` to NaN.

**Issue:** The check is `tmp < min_float || tmp > max_float`. Both comparisons are false for NaN, so NaN passes every bounded restriction. `SignalToNoiseEstimatorMedian` then runs with `win_len = NaN` or `auto_max_stdev_factor = NaN`, and FeatureFinderAlgorithmPicked with a NaN `isotopic_pattern:intensity_percentage_optional` (CPP-324); the restriction the parameter declares is simply not enforced.

**Proposed C++ fix:** Reject a non-finite value explicitly before the range comparison.

**Evidence:** Executed on the Release build: the `win_nan` case (`../oracle/sne-completion`), and the `p_ipo_nan*` cases of `../oracle/ffap-complete-fix3` for the FeatureFinderAlgorithmPicked path.

**Rust handling:** `NoiseCompatibility::source_value_domain` accepts NaN where the source does, under `PickingCompatibility::source()`; the native profile refuses it. FeatureFinderAlgorithmPicked applies the source's comparisons in its own module (see CPP-324), and the crate's shared `Param` check keeps its native refusal for every other handler.

## CPP-335 — FileInfo declares CorruptionInfo and DetailInfo and never fills them

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`.

**Status:** Executed.

**Affected file/function:** `src/openms/include/OpenMS/FORMAT/FileInfo.h:206-217` and `src/openms/source/FORMAT/FileInfo.cpp:1779-1964`, `FileInfo::report_`.

**Trigger:** `FileInfo::run` with `Options::detailed` or `Options::check_corrupt` set.

**Issue:** The header declares `CorruptionInfo {performed, errors, warnings}` and `DetailInfo {performed, lines}` as part of the structured `Result`, and documents the latter as "kept as pre-rendered lines". `report_` never assigns to either: every `-d` and `-c` message goes straight into the text stream. The class's stated purpose is that "the result is consumable directly from C++ and pyOpenMS without any stream" (`FileInfo.h:36-38`), and for these two flags it is not — a caller must re-parse the rendered text. Both flags are also the ones whose findings a caller would most want as data.

**Proposed C++ fix:** Push each message onto the matching vector as it is written, and set `performed` where the flag ran.

**Evidence:** Source review of the pinned header and implementation (a grep for `r.corruption` and `r.detail` in `FileInfo.cpp` returns nothing), plus the executed Release runs of `../oracle/a6-fileinfo`, whose `-d` and `-c` reports are entirely in the text stream.

**Rust handling:** Reproduced. `FileInfoResult::corruption` and `::detail` stay at their defaults after a run that requested both flags, and `tests/file_info_checks.rs` asserts that for every compared case. `docs/FILE_INFO_CHECKS_SUPPORT.md` records it as native difference 1.

## CPP-336 — FileInfo's `-d` listing reads front() and back() of an SRM chromatogram unchecked

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

**Status:** Source review (the input that triggers it cannot be produced by any pinned loader).

**Affected file/function:** `src/openms/source/FORMAT/FileInfo.cpp:1792`, `FileInfo::report_`.

**Trigger:** `-d` on a peak file holding a selected-reaction-monitoring chromatogram with no points.

**Issue:** `ms.front().getRT()` and `ms.back().getRT()` are called on every SRM chromatogram without checking `ms.empty()`. Both are undefined on an empty container.

**Proposed C++ fix:** Skip an empty chromatogram, or print a placeholder.

**Evidence:** Source review. No pinned loader produces an empty chromatogram — the mzML reader gives every chromatogram its points and `ChromatogramTools::convertSpectraToChromatograms` builds one point per source spectrum — so the line is currently unreachable rather than latent.

**Rust handling:** Refused with `Error::InvalidValue` naming the chromatogram, because there is no defined behaviour to reproduce (native difference 2 of `docs/FILE_INFO_CHECKS_SUPPORT.md`), with a unit test in `src/format/file_info/checks.rs`.

## CPP-337 — IndexedMzMLDecoder skips the first child of every `<index>` element

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576` (`FileInfo` sha256 `5d82c8a7248c9ccb37c6f4e64a0ed2202f375418befd095261590935de1172dc`), on `ibminode06`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FORMAT/HANDLERS/IndexedMzMLDecoder.cpp`, `IndexedMzMLDecoder::domParseIndexedEnd_` (`:206`) — `:280-282` declares `firstChild`, `lastChild` and `iter`, and `:290-293` is the walk itself.

**Trigger:** Reading the index of an `indexedmzML` whose `<index>` element has **no text node immediately after its opening tag** — that is, whose first child is an `<offset>` element. Whitespace elsewhere inside `<index>`, between the offsets or before `</index>`, does not avoid it.

**Issue:** The child walk is

    DOMNode* firstChild = currentNode->getFirstChild();   // :280
    DOMNode* lastChild  = currentNode->getLastChild();    // :281
    DOMNode* iter = firstChild;                           // :282
    while (iter != lastChild)                             // :290
    {
      iter = iter->getNextSibling();                      // :292
      DOMNode* currentONode = iter;                       // :293
      ...

It advances before it reads, so the first child is never looked at, and only that one position matters. Every index an OpenMS writer produces is indented, which puts a text node first and hides the defect; an index whose first child is an `<offset>` silently loses that offset, and an index whose section holds exactly one offset parses as EMPTY. Nothing reports an error: `parseOffsets` returns 0 (`:339`, against its `return -1` failure sites at `:97-99`, `:115-118` and `:332`) and the caller sees a short offset vector, so `FileInfo -i` prints a spectrum count that is too low and still exits 0.

**Proposed C++ fix:** Iterate the children with `for (DOMNode* iter = firstChild; iter != nullptr; iter = iter->getNextSibling())` and keep the existing `ELEMENT_NODE` test. (The loop's own NOTE explains why `DOMNodeList::item` is avoided; the sibling walk is right, only its first step is.)

**Evidence:** Executed, three ways, all on the Release build named above.

* **Pinned oracle case:** `../oracle/a6-fileinfo` case `i_offsets_unspaced` on `index_offsets_unspaced.mzML` (1247 bytes, sha256 `fd9dff92db516c03…`), two `<offset>` children and no whitespace inside `<index>`: the Release `FileInfo` prints "Found a valid indexed mzML XML File with 1 spectra and 0 chromatograms.", exit 0. Its control `index_window_above.mzML`, the same generator with a newline between the index elements, prints 1 for its one offset.
* **Probe,** `ibminode06` `/scratch/kohlbach/a6-fileinfo-close/probe`, each case run twice with identical output, inputs built by the pinned generator `../oracle/a6-fileinfo/scripts/make_index_window_fixtures.py` (sha256 `cf8d29a579391472…`) through its own `build(pad, entries, separator)`:

      build(250, 1, "\n")  1217 B  sha256 f194cd1a3381368797e78f6731b162d45b04492b94c101c4e7a2f5eca6a84a8a
                                -> "... with 1 spectra and 0 chromatograms."   exit 0
      build(250, 1, "")    1213 B  sha256 542ddaa763e20ce1acb17b38f973e83eb7c602a522b8019f1d1c587473382c4c
                                -> "... with 0 spectra and 0 chromatograms."   exit 0
      build(250, 3, "")    1281 B  sha256 fdcbc198100cd5d048be807670c238edf7596ab1fecbf3bb1e132657d1644ca0
                                -> "... with 2 spectra and 0 chromatograms."   exit 0

  The first is byte-identical to the pinned fixture `index_window_above.mzML` (sha256 `f194cd1a33813687…`), which is what makes the other two re-derivable: exactly one offset is lost per `<index>`, and a single-offset section disappears entirely.
* **Which whitespace matters,** measured by the closing reviewer on the same build: with a newline only *after* the opening `<index …>` tag the C++ counts 2 of 2 offsets; with a newline only *between* the two offsets, or only *before* `</index>`, it counts 1. The affected class is therefore "first child is an `<offset>`", which is what the mechanism above predicts.
* **Source:** the pinned lines above.

**Rust handling:** NOT reproduced, and deliberately so by the OWNING package, whose support document already states the decision: `docs/INDEXED_MZML_SUPPORT.md`, "the native parser corrects the source DOM sibling loop that skips an offset when it is the first child without preceding whitespace." Reproducing the skip would break random access, because a dropped offset is a record that cannot be found. A6 records only the consequence for `-i`: native difference 12 of `docs/FILE_INFO_CHECKS_SUPPORT.md`, native difference 9 of `docs/TOPP_FILE_INFO_SUPPORT.md` and `known_gaps` of both manifests, with `src/format/indexed_mzml.rs` named as owner; `tests/file_info_checks.rs` pins both sides.

Note under the same entry, with no number of its own because the source's behaviour there is indeterminate rather than wrong: `IndexedMzMLDecoder::findIndexListOffset` (`:165-168`) does `new char[buffersize+1]`, `f.seekg(-buffersize, f.end)`, `f.read(...)`. On a file shorter than `buffersize` (1023 by default, `IndexedMzMLDecoder.h:80`) the seek fails, the read writes nothing, and the regex at `:179-181` searches uninitialised memory; the `else` branch at `:198-200` then prints that memory to `std::cerr` after "Maybe this is not a indexedMzML.". Executed: two oracle runs of the same 967-byte file produced two different stderr dumps and the same report and exit code. A file that short cannot hold a real spectrum, so no writer-produced mzML reaches it, but `if (length < buffersize) buffersize = length;` before the seek would remove the read of uninitialised memory. The Rust side already does exactly that, and `docs/INDEXED_MZML_SUPPORT.md` records it ("small files search their full contents"), so this too is an owning-package decision A6 only inherits and prints.

## CPP-338 — FeatureFinderCentroided crashes intermittently at `-threads 32`

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Observed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`, binary sha256 `2781dd7cab48f482118375754038644251705c8bafab1958e310e4063dffb323`.

**Status:** Observed once, on the Release build. **Not diagnosed** — no faulting code identified, no reproduction under a debugger attempted.

**Affected file/function:** Unknown. The process died before writing any output; nothing in this evidence localises the fault.

**Trigger:** `FeatureFinderCentroided -threads 32` with `OMP_NUM_THREADS=32` on a 4,000-spectrum centroided Velos mzML. Intermittent: one failure in six executions of that cell in the wave-6 benchmark run (2026-09-18, `ibminode05`) and none in six in the wave-4 run (2026-09-16, same node, same binary, same input, same INI), so 1 in 12 observed at 32 threads and 0 in 12 at one thread.

**Issue:** The process terminated with `SIGSEGV` (signal 11, launcher exit 139) after 7.094 s wall and 12.818 s user, with 64 live threads and a peak RSS of 340,628 KiB, and wrote 0 output files. The four repetitions of the same cell that succeeded ran the same binary sha256, the same INI sha256 `2869134aeb3f98ed181fb4b8e089fd955a339d77958fb149f39a0f9d98a39b50`, the same input sha256 `6d0c151853d717aca0ce8ae625d5021aa10774bb9bda5921d21f164d74b04006` and the same environment, and took 25.2 s each. The node was quiet: foreign CPU 0.0092 per core, the pre-cell gate value 0.0098, `majflt` 0, and the repetition was **not** load-flagged. The last stdout line was `Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.`; the last stderr line was the known non-fatal `Non-fatal error while loading '…first4000.mzML': DateTime conversion error of "-infinity"`.

**Proposed C++ fix:** None proposed — the fault is not localised. The next step would be a repeat run of this cell under a debugger or with a core dump enabled, which needs a separate run.

**Evidence:** `/ceph/ibmi/abi/oliver/bench/openms4/results/2026-09-18-w6refresh/w6-ffc-subset/results.jsonl`, the record with `tool = FeatureFinderCentroided`, `impl = cpp-release`, `threads = 32`, `rep = 2`, `status = nonzero_exit`, `signal = 11`. It is the only record with a status other than `ok` or `case_load` in either wave, the pilot sub-runs of both waves included. The census behind that statement covers 334 timing executions in wave 6 and 228 in wave 4 counting the four measured sub-runs of each wave only, and 337 and 244 with the pilot sub-runs added. Reported in `docs/BENCHMARKS.md` §3.9.

**Rust handling:** Not applicable — this is a defect observation against the pinned C++ Release build, not a portability decision. For the record, neither Rust build (`+fma` default and the `-fma` opt-out) failed in any of its 24 executions of this cell, and none of the port's 186 measured repetitions in the run (96 `rust-release` + 90 `rust-nofma`) failed. Of the run's 282 measured repetitions, 281 succeeded; the one that did not is this C++ crash.

## CPP-339 — MzMLFile::transform's first pass leaves a ProgressLogger running, so every low-memory TOPP path fails on an mzML with both spectra and chromatograms

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`, on `ibminode06`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FORMAT/MzMLFile.cpp:178-231`, `MzMLFile::transform` and `transformFirstPass_`; `src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp:997-1006`, `MzMLHandler::startElement` for `chromatogramList` under `LD_RAWCOUNTS`.

**Trigger:** Any `MzMLFile::transform` call with progress logging on, over an mzML that holds both a `spectrumList` and a `chromatogramList`. Through the TOPP tools: `PeakPickerHiRes -processOption lowmemory`, `NoiseFilterGaussian -processOption lowmemory`, `NoiseFilterSGolay -processOption lowmemory`, `FileConverter -process_lowmemory` and `PeakPickerIM -processOption lowmemory`. `-no_progress` avoids it; `-test` does not.

**Issue:** `transform` parses the file twice through one `ProgressLogger` (`*this`). In the first pass the handler runs with `XMLHandler::LD_RAWCOUNTS`. At `<spectrumList>` it calls `logger_.startProgress(0, scan_count_total_, "loading spectra list")` (`MzMLHandler.cpp:966`) and at `<chromatogramList>` `logger_.startProgress(0, chrom_count_total_, "loading chromatogram list")` (`MzMLHandler.cpp:997`); immediately after the second one, now holding both counts, it throws `EndParsingSoftly` (`MzMLHandler.cpp:1001-1006`). The matching `logger_.endProgress()` at `</chromatogramList>` (`MzMLHandler.cpp:1493-1498`) and the outer `pg_outer.endProgress()` at `</mzML>` (`MzMLHandler.cpp:1524`) are therefore never reached, and the shared logger's `StopWatch` is still running when the second pass calls `startProgress` again. The `is_running_` guard in `StopWatch::start()` then throws `Exception::Precondition` (`src/openms/source/SYSTEM/StopWatch.cpp:41-43`) — an unconditional runtime check, not the assertions-only `OPENMS_PRECONDITION` macro, which `src/openms/include/OpenMS/CONCEPT/Macros.h:91` defines as nothing when `OPENMS_ASSERTIONS` is off; that is why a Release build fails too, and the failure below was observed on one. The tool exits 3 with

    Error: Unable to read file (- due to that error of type Precondition failed in: .../StopWatch.cpp@43-void OpenMS::StopWatch::start())

having written nothing. An input with only one record kind never takes the early throw, so its `startProgress`/`endProgress` pair is balanced — which is why the upstream registrations `TOPP_PeakPickerHiRes_3` and `_4` do not see it: their inputs hold spectra only and chromatograms only.

**Proposed C++ fix:** End the progress before throwing `EndParsingSoftly` in the two `LD_RAWCOUNTS` early-exit branches, or do not start a progress at all when `load_detail_ == LD_RAWCOUNTS`, since the first pass reports no useful progress anyway.

**Evidence:** `../oracle/p4-lowmemory/logs/probe_06.log`, section `(a)`. Minimal reproduction `../oracle/p4-lowmemory/outputs/probe_both_input.mzML`, 450,674 bytes, sha256 `23cb36d87ace821443539ec7a391711482ee6aad7100fb2c0243ccfcc6077f78`, produced on the node by that build's own `FileMerger` from the two upstream inputs `PeakPickerHiRes_input.mzML` (5 spectra) and `PeakPickerHiRes_2_input.mzML` (5 chromatograms). `PeakPickerHiRes -in both.mzML -out x.mzML -processOption lowmemory` exits 3; with `-no_progress` it exits 0 and writes 66,243 bytes; with `-test` it exits 3. Controls on the same build: the spectra-only and chromatogram-only inputs exit 0 in the same mode, and the in-memory mode on the two-kind input exits 0. The 2.3 GB benchmark input `UK222.mzML` (40,856 spectra, one chromatogram) fails the same way (`logs/rss_06.log`).

**Rust handling:** The port's `mzml::transform` carries no shared progress logger — the crate has no progress logging on this path at all — so `PeakPickerHiRes -processOption lowmemory` completes on every one of those inputs, including the 2.3 GB one. Recorded as native difference 4 in `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md`.

## CPP-340 — MSDataWritingConsumer never checks its output stream, so a low-memory TOPP run whose output cannot be created exits 0 having written nothing

**Source revision:** `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Executed on the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576`, on `ibminode06`.

**Status:** Executed.

**Affected file/function:** `src/openms/source/FORMAT/DATAACCESS/MSDataWritingConsumer.cpp:19-34`, the constructor, whose `ofs_.open(filename.c_str(), std::ios::out | std::ios::binary)` at `:33` is never checked; `doCleanup_` at `:151-173`, which closes the same stream at `:172` without checking it; `src/PeakPickerHiRes.cpp:170-186` at the TOPP pin `174b576`, `doLowMemAlgorithm`, which returns `EXECUTION_OK` unconditionally at `:185`.

**Trigger:** Any low-memory TOPP run whose `-out` cannot be opened for writing but passes the framework's pre-check. The reachable case is an `-out` that names an existing **directory**: `TOPPBase`'s writability pre-check accepts it, and `std::ofstream::open` then fails.

**Issue:** The constructor never tests `ofs_.is_open()` or the stream state, and nothing on the writing path does either — `consumeSpectrum` and `doCleanup_` stream into a failed `ofstream`, which silently discards everything. `doLowMemAlgorithm` returns `EXECUTION_OK`, so the tool exits 0 with an empty standard error having produced no output at all. The in-memory path of the same tool reports the same situation properly, because `MzMLFile::store` throws `UnableToCreateFile`. Silent data loss, and the two process options of one tool disagree about whether the run succeeded.

**Proposed C++ fix:** Throw `Exception::UnableToCreateFile` from the constructor when the stream is not open, as `MzMLFile::store` does, and check `ofs_` again in `doCleanup_` so that a write failure partway through is reported rather than swallowed.

**Evidence:** Executed on the Release build, `../oracle/p4-lowmemory/logs/closediff1_06.log` section F and `logs/closediff3_06.log` section D. With `-out` naming a pre-existing directory:

| | `-processOption lowmemory` | `-processOption inmemory` |
| --- | --- | --- |
| C++ Release | **exit 0**, empty stderr, nothing written | exit 5, `Error: Unable to write file (the file '…' could not be created. )` |

Two controls that do **not** reach the consumer, because the framework's pre-check catches them, exit 5 with `Cannot write output file given from parameter '-out'!` in both process options: an `-out` naming an existing read-only file, and an `-out` under a directory that does not exist.

**Rust handling:** `MSDataWritingConsumer::create` returns `Error::Io` when the file cannot be created, and `run_low_memory` propagates it, so the port answers `Error: Unexpected internal error (Is a directory (os error 21))` with `UNKNOWN_ERROR` in **both** process options. This is the one row of that tool's low-memory divergence table the port deliberately does not reproduce, because reproducing it means swallowing an I/O error on the one file the run exists to produce; recorded as native difference 13 of `docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md`, with the two controls that agree exactly, and pinned by `an_out_that_names_a_directory_is_reported_in_both_modes`.
