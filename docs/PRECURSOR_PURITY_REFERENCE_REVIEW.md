# Precursor purity reference review

[precursor_purity_reference.rs](../tests/precursor_purity_reference.rs) checks
scalar isolation purity, experiment parent selection, fuzzy scan estimates,
interpolation, and SPS fragment matching. The
[provenance manifest](../tests/data/precursor_purity_provenance.json) records ten
pinned source hashes and the independent extraction method. No C++ program was
built or executed. Tests use only the installed Rust package and bundled TSV
fixtures; they do not require the source checkout or temporary extraction files.

## Fixture and reference strength

The source `PrecursorPurity_input.mzML` contains seven spectra: two MS1 scans
with 2,003 and 1,869 peaks, separated by five MS2 scans. The source stores m/z
as uncompressed little-endian binary64 and intensity as uncompressed
little-endian binary32. Python's standard XML, base64, and binary unpacking
facilities decoded the arrays independently of the native mzML reader.

The [peak fixture](../tests/data/precursor_purity_peaks.tsv) retains 72 exact
observations: 37 from the earlier scan and 35 from the later scan. For each
source precursor in each MS1 scan, it includes the isolation window expanded
by 20 ppm and two immediately adjacent source peaks on either side. Duplicate
indices are merged. Authoritative hexadecimal IEEE bits accompany decimal
values and original peak indices. The
[spectrum fixture](../tests/data/precursor_purity_spectra.tsv) retains all seven
native IDs, scan levels, retention times in seconds, and five selected-ion
m/z, charge, intensity, and isolation offsets. All isolation targets equal the
selected-ion m/z, and the source has no precursor spectrum reference.

A separate Python calculation of the scalar source algorithm checked that
all seven scalar/map cases give identical results on the full original MS1
arrays and the compact extraction. The compact data is sufficient for these
cases, but is not a replacement for arbitrary whole-scan operations.

A [byte-identical original mzML](../tests/data/precursor_purity_input.mzML) is
also retained for the separate end-to-end reader-to-purity tests. Its full
115,055 bytes have the same SHA-256 as the pinned source fixture. This transport
check complements the independent TSV references; their expected observations
were not decoded by the Rust reader under test.

The [score fixture](../tests/data/precursor_purity_scores.tsv) keeps source test
literals separate from independently calculated values. Source totals are
occasionally truncated, and source ratios have five decimal places; comparisons
allow an absolute 0.00001 for these literals. Exact decoded intensity totals
are compared exactly, and independently calculated binary64 ratios use an
absolute tolerance of 1e-14. The latter values are derived references, not
captured C++ output. Original indices also identify every interfering peak.

The [SPS fixture](../tests/data/precursor_purity_sps.tsv) contains the eight
literal fragment-matching cases from the class test. Its ppm example uses
20 ppm, matching the actual source argument despite a comment saying 10 ppm.
A separate upstream self-consistency case uses the compact theoretical helper;
its expected 32 values follow from eight lengths, two b/y series, and two
charges. This test is explicitly self-consistency, while the literal cases
and the analytical binary32 boundary check provide separate numerical evidence.

Upstream interpolation tests contain fallback comparisons and a range check,
but no literal numeric purity table. Synthetic fuzzy and interpolation cases
in this suite derive expected values from the source arithmetic and branches;
they are not presented as upstream numerical goldens.

## Scientific conventions checked

- Scalar isolation and mass-match endpoints are inclusive. Tolerance is
  doubled, with ppm evaluated as `mz * tolerance * 2 * 1e-6`. Equal-distance
  nearest ties favor the lower m/z, and each matched observation is removed
  before the next isotope search.
- Scalar spacing is the source C13/C12 mass difference, 1.0033548378 Da,
  divided by the absolute charge; zero charge means one. The scalar routine
  can match other isotopes when the monoisotopic peak is absent, contrary to
  its header description.
- Scalar totals accumulate binary32 observations in binary64. In contrast,
  the fuzzy algorithm accumulates and divides in binary32 before widening
  its ratio. Two intensity-one neighbors around an intensity-2^24 precursor
  therefore leave the fuzzy total unchanged after sequential rounding.
- Fuzzy spacing uses the source neutron mass, 1.00866491566 Da. Exact strict
  window endpoints contribute half intensity; exact outer fuzzy endpoints
  are excluded from the denominator's neighbor scans. The initial globally
  nearest precursor is included without a mass tolerance.
- Fuzzy isotope lookup compares `lower_bound(expected)` with its physical
  successor. It never compares the predecessor. A close isotope immediately
  below the expected position can therefore be ignored; the reference test
  pins this defined source behavior rather than substituting true nearest
  matching. A physical candidate beyond the logical fuzzy range can contribute
  target intensity while being excluded from the denominator; the source can
  therefore also produce a single-scan purity above one. The implementation
  owner's separate test pins this case, and native code does not clip it.
- Fuzzy interpolation uses absolute RT differences and does not clamp its
  output. A finite extrapolation can yield a value below zero. Invalid next
  indices, a next scan of the wrong level, or a zero RT denominator fall back
  to the earlier scan estimate.
- The batch scalar map uses the preceding parent only; the header's statement
  about combining both neighboring MS1 scans does not describe the actual
  implementation. Parent lookup searches backward for a reference with a
  matching native ID and one lower MS level, then falls back to the nearest
  preceding scan with one lower level. A reference to a later scan does not
  select that scan.
- The source class test inserts zero-valued MS1/random keys through C++ map
  indexing. Those keys are not returned by the batch calculation. The native
  reference test expects exactly the five actual MS2 keys.
- SPS matches each precursor separately, including duplicate m/z values, and
  ignores precursor charge when selecting theoretical fragment charges. It
  narrows both tolerance bounds to binary32 before inclusive comparison with
  sorted binary32 b/y masses. Zero tolerance can consequently accept a
  distinct binary64 m/z that rounds to the same binary32 value.

## Checked boundaries and implementation review

The source fuzzy lookup can dereference beyond the physical spectrum, while
negative charge can reverse isotope stepping and prevent termination. Source
empty experiments and certain missing-precursor cases also index unchecked
containers. Native validation and bounded work must handle these branches
without reproducing undefined behavior. Those checked policies are distinct
from the finite scientific conventions above.

The scalar implementation was inspected for source isotope initialization,
operation order, inclusive searches, lower-tie selection, greedy removal, and
separate binary64 accumulation. Its work accounting charges copied isolation
observations, isotope iterations, binary searches, and physical removals. SPS
reuses the previously tested compact theoretical helper and preserves
source division-first ppm conversion and binary32 matching bounds.

The fuzzy implementation preserves the physical successor search when the
successor exists, even outside the logical binary-search interval. At physical
end, native code uses a remaining lower-bound candidate or reports no match;
it does not dereference missing observations. It rejects negative charge,
nonfinite arithmetic, nonpositive totals, negative fuzzy boundaries, and
non-advancing isotope steps. Window planning and both interpolation calculations
share bounded work accounting. Fuzzy output preserves the source's finite
unclamped extrapolation and zero-width early-stop behavior.

Batch suitability failures return checked errors instead of source warnings
and an empty map. With missing parents explicitly permitted, an unresolved MS2
parent retains a default zero score. Source spectrum references are searched
only among preceding scans at the required MS level. Review confirmed that the
native kernel lookup preserves the source's reference-first, fallback-second
ordering.

Review found that calling the general spectrum validator repeatedly could
traverse large auxiliary-array, identification, or CV collections without
charging that work to the purity budget. The integrated fix validates the
numerical observations and selected precursor fields used by purity directly,
with bounded loops, and does not recursively inspect unrelated annotations.
General container validation remains available separately. The owner's
regression covers 10,000 auxiliary placeholders reused by 10,000 MS2 scans,
including irrelevant malformed annotation content. Fuzzy precursor validation
uses the same constant-size numerical checks before planning windows.

All thirteen independent reference tests passed in the integrated current-compiler
run. The final validation report records the complete compiler and feature
matrix after integration.

An independent acquisition-integration review also found that a signed selected
m/z could become a negative fallback isolation target during mzML output. The
writer now checks the effective target before writing, with spectrum and
chromatogram regression cases. Direct chromatogram round trips exercise every
activation method and mobility quantity through the shared precursor adapter.
The original complete mzML fixture reproduces all seven scalar/map cases after
an explicit in-memory encoding declaration change; acquisition writer output
passes independent XSD validation. These checks complement the numerical
reference suite and do not imply complete PSI semantic validation.
