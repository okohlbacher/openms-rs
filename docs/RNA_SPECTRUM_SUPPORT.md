# RNA theoretical spectra

`NucleicAcidSpectrumGenerator` implements both public generation operations of
the pinned OpenMS class. It generates all nine RNA fragment series, ambiguous
a-B alternatives, signed charge states and optional precursor peaks. It uses
owned [RNA records and sequences](RNA_SUPPORT.md); no registry lookup or external
file is needed after an `NASequence` has been constructed.

```rust
use openms::chemistry::{NASequence, NucleicAcidSpectrumGenerator};
use std::collections::BTreeSet;

let sequence = NASequence::parse("[m1A]UC[C*]AC[A*]Gp")?;
let generator = NucleicAcidSpectrumGenerator {
    add_metainfo: true,
    add_first_prefix_ion: true,
    add_c_ions: true,
    add_y_ions: false,
    add_b_ions: false,
    ..Default::default()
};
let spectrum = generator.generate(&sequence, -1, -2)?;
let multiple = generator.generate_multiple(
    &sequence, &BTreeSet::from([-1, -3, -5]), -1,
)?;
```

The complete configuration contains the source's 13 boolean options and ten
f64 intensities. Only `add_b_ions` and `add_y_ions` default true; all intensities
default to 1.0. The a-B fields are `add_a_minus_b_ions` and
`a_minus_b_intensity`. Other names follow the source snake-case names.
Configurations are ordinary owned values supporting clone, debug and equality.
There are no hidden parameter caches.

## Public operations

- `generate(sequence, min_charge, max_charge)` creates one new spectrum.
- `append_to(spectrum, sequence, min_charge, max_charge)` appends and stably sorts
  all old and new peaks, with atomic updates to the affected arrays.
- `generate_multiple(sequence, &BTreeSet<i32>, base_charge)` returns a map of
  spectra with the source's cumulative charge-state behavior.
- `replace_multiple(output_map, sequence, charges, base_charge)` replaces the map
  atomically. An empty requested set clears it, without inspecting unused input
  settings or sequence chemistry.

New spectra retain source defaults: MS level one, unknown spectrum type and no
acquisition precursor records. An optional M peak is a spectral peak; it does
not populate `MSSpectrum::precursors`. Existing settings, identifiers, attached
identifications and metadata are retained by append. A workflow that needs MS2
or centroid metadata must set it explicitly.

## Mass conventions and source arithmetic

Fragment masses use the records' **independently declared monoisotopic masses**,
not their formula-derived masses. Terminal contributions are their declared mono
mass minus the empirical mass of natural H. Linkage and ion offsets are empirical
formula masses. These choices preserve the source's rounding and custom-record
behavior and differ deliberately from `NASequence::mono_mass`.

Let `m[i]` be declared nucleoside mass; F/T the five/three-prime contribution;
K the mass of H-1PO2; and S[i] the mass of SO-1 after a residue code ending `*`,
or zero otherwise. A final starred residue has no following linkage.

Prefix masses begin with `m[0] + F`, then accumulate, in order, the previous
prefix, next nucleoside, K and the preceding sulfur correction. Suffix masses
begin with `m[last] + T`, then accumulate previous suffix, next nucleoside, K
and the sulfur correction on that newly included residue. No full-length
fragment is emitted.

| Series | Offset and sulfur correction |
| --- | --- |
| a | Prefix minus H2O |
| b | Prefix |
| c | Prefix plus K; add SO-1 when its final residue is starred |
| d | Prefix plus HPO3; same final-residue correction |
| w | Suffix plus HPO3; add SO-1 on the boundary preceding the suffix |
| x | Suffix plus K; same boundary correction |
| y | Suffix |
| z | Suffix minus H2O |

`add_first_prefix_ion=false` omits a1/b1/c1/d1/a1-B. It does **not** omit
w1/x1/y1/z1, despite the source option description listing those ions.

For a-B, use the current residue's base-loss formula mass. At the first residue,
subtract H4O2 and ignore the five-prime end. For later residues, add the preceding
prefix plus H-5P, then the preceding residue's sulfur correction. A record whose
code ends `?` emits two peaks separated by CH2, each at half intensity, with
identical aN-B labels. Codes ending `?*` emit one peak, matching the record's
source `is_ambiguous` predicate. The registry's two alternative records are not
used for this calculation.

Intensity is narrowed to f32 only for an emitted peak. Ambiguous a-B intensity
is halved in f64 before narrowing; a large original value is permitted if its
emitted half is finite f32. Finite signed intensities are retained. Unused
intensity settings are not validated.

### Configuration-dependent precursor mass

The source reuses whichever fragment arrays were actually constructed:

| Available arrays | Precursor expression | Omitted sulfur correction |
| --- | --- | --- |
| Both | first prefix + last suffix + K | First linkage |
| Prefix only | last prefix + last nucleoside + K + T | Last linkage |
| Suffix only | last suffix + first nucleoside + K + F | First linkage |
| Neither | `NASequence::mono_mass(Full, 0)` | None |

The arithmetic order and missing-link corrections are retained. The last row
uses formula chemistry, whereas the other three use declared masses. Enabling
a fragment series can therefore change precursor mass. Prefix construction on
short sequences also depends on `add_first_prefix_ion`. This implementation does
not silently replace all paths with a formula-derived precursor.

## Charge and multiple-spectrum differences

An ordinary charged peak uses `abs(uncharged_mass / charge + PROTON_MASS_U)`.
It does not use natural-H/electron formula mass conversion.

Single-spectrum generation selects negative mode only when both endpoints are
negative, rejects mixed nonzero signs, and swaps endpoints by charge magnitude.
It emits magnitudes in the requested inclusive range subject to the strict
additional bound `magnitude < sequence.len()`. A final-only precursor is emitted
only at the originally requested maximum magnitude; if that magnitude is excluded
by the length bound, no final-only precursor is emitted at all.

Multiple-spectrum generation follows these separate source rules:

- There is no charge-magnitude-versus-length restriction.
- The smallest requested key selects positive or negative mode. The body follows
  this rule even for mixed-sign sets; it does not enforce the header's uniform-sign
  prose. Negative mode converts a positive base charge to its negative.
- Targets below the base magnitude are skipped. With metadata enabled they remain
  as empty spectra with named arrays; without metadata their keys are absent.
- Fragments and all-charge precursor peaks accumulate across successive targets.
  Final-only precursor peaks are excluded from the data inherited by the next
  target. Copies retain source insertion order until the final stable sort.
- A final-only M peak uses the **post-loop charge counter**: target +1 gets M at
  charge +2, target -3 gets M at charge -4. Charges stores that actual counter.
- The positive-mode final-only M expression has no absolute value. A custom
  negative precursor mass can produce a finite negative m/z here. The negative
  branch and ordinary copied peaks retain their absolute-value operation.

Zero charges are rejected when a peak would actually be divided by zero. Defined
no-row cases remain allowed: for example, the zero-charge copy can exclude the
only precursor, then a later nonzero copy can emit it. Empty sequence with a
processed final-only multi precursor returns an error rather than dereferencing
an absent source row. Charge absolute-value/counter overflow and nonfinite mass
arithmetic return checked errors. No unsigned wrapping or overflowing sign
multiplication is used.

## Annotation and append policy

With `add_metainfo`, generated arrays are `IonNames` and signed integer `Charges`.
Labels are a1, y2, a3-B or M without charge suffixes. Duplicate peaks and duplicate
labels are retained. All peak/array ordering is stable by m/z.

Append follows the source's **first-array** behavior rather than searching for
arrays by name. Existing first integer/string array names are preserved. The
following checked adaptations keep native arrays aligned:

- If additions need annotations for old unannotated peaks, fill missing names
  with empty strings and missing charges with zero.
- If no peaks are added, empty annotation placeholders stay empty. Old peaks and
  populated arrays are still sorted together. Metadata-enabled calls create the
  named placeholders when their array lists were absent.
- Other populated arrays have no defined values for new theoretical peaks and
  cause an error. Empty placeholders are retained. With metadata disabled,
  populated arrays likewise cause an error if new peaks would make them short.
- Malformed nonempty array lengths and nonfinite old peak coordinates/intensities
  are rejected before commit. Unused retention time, acquisition metadata and
  attached identification graphs are not validated or cloned.

The source can mutate before its final sort detects a short array, and its
already-sorted shortcut can leave malformed arrays unchecked. Native append
stages affected peak and array data and commits only after every check succeeds.
Existing array names and unrelated metadata retain their allocations.

## Bounds and validation evidence

Fixed public constants in `chemistry::nucleic_acid_spectrum_generator` bound an
operation to one million residues, 100,000 combined single-spectrum peaks or
100,000 peaks summed over the returned map, 4,096 requested keys, and 1,024
arrays per appended spectrum. Uncharged and cumulative working buffers also
have the 100,000-peak bound.

Every operation shares 50 million charged work units, 256 MiB of conservative
cumulative allocation accounting and 16 MiB of copied label bytes. These cover
mass construction, charge spans (including empty copies), copied output maps,
name vectors/text and stable-sort work/scratch. The formula-only precursor uses
the same remaining work and byte counters through NASequence's internal helper;
it does not start a fresh formula allowance. Output counts are preflighted across
the entire map, not reset for each target. Large empty charge spans are bounded
without iterating pointlessly through them.

The [independent RNA processing review](RNA_PROCESSING_REFERENCE_REVIEW.md)
records 132 source ion literals, of which 126 are actually compared upstream,
from the ordinary and sulfur-linked eight-residue examples. Native literal tests
use an explicit 0.001 Da tolerance, not an assumed ClassTest macro tolerance.
The upstream multiple-spectrum test constructs modified options but never applies
them: its real comparison uses default b/y ions without metadata. Independent
native tests additionally apply all nine series with annotations and verify
multi/single agreement where no precursor and length-limit difference applies.

[Direct tests](../tests/nucleic_acid_spectrum_generator.rs) cover charge/metadata
quirks, signed custom chemistry, half-intensity conversion, atomic failure and
aggregate limits. The private shared-budget regression proves that a formula
fallback cannot reset the generator allowance. No C++ reference was built or run.
Isotopic RNA envelopes, search scoring and newer RNA identification graphs are
separate work; they are not capabilities of this source class.
