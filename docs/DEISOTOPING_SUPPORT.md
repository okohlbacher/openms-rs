# Isotope cluster removal

`processing::deisotoping::Deisotoper` ports the simple
`Deisotoper::deisotopeAndSingleCharge` algorithm at OpenMS4-core revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. The separate
`AveragineDeisotoper` ports the source Poisson/KL-divergence method; see its
[model, preprocessing and selection conventions](AVERAGINE_DEISOTOPING_SUPPORT.md).
The details below describe the simple method.

```rust
use openms::processing::{SpectrumFilter, deisotoping::Deisotoper};
use openms::chemistry::C13C12_MASSDIFF_U;
use openms::{MSSpectrum, Peak1D};

let mut spectrum = MSSpectrum::from_peaks((0..3).map(|i|
    Peak1D::new(200.0 + i as f64*C13C12_MASSDIFF_U/2.0, 10.0-i as f32)
).collect());
Deisotoper { annotate_charge: true, ..Default::default() }
    .filter_spectrum(&mut spectrum)?;
assert_eq!(spectrum.len(), 1);
assert_eq!(spectrum.integer_data_arrays[0].data, [2]);
# Ok::<(), openms::Error>(())
```

## Algorithm and source conventions

The input must have sorted finite m/z values and nonnegative finite intensity.
For each unassigned seed in ascending m/z order, the algorithm tries charges
from maximum to minimum and takes the first qualifying ladder. Each expected
isotope is `seed_mz + isotope_number * C13C12_MASSDIFF_U / charge`; its nearest
observed peak must lie inside the inclusive tolerance window. Ppm windows use
the seed's m/z throughout the ladder, and nearest-distance ties prefer lower
m/z. Input order is checked once; each lookup then uses binary search.

Defaults are 10 ppm, charges 1–3, 3–10 isotope peaks, a decreasing-intensity
model starting at isotope index 2, retention of unassigned peaks, conversion to
single charge, and no added annotations, summed intensity or shared isotope extensions. The first heavy
isotope may exceed the monoisotopic intensity under that default. Missing peaks
or a failed intensity comparison stop a ladder; an already long enough prefix
can still qualify. These choices follow the source except for the native
disjoint-membership default; enable `allow_shared_isotopes` for source sharing.

One precursor with a known charge and positive neutral mass limits candidate
fragment neutral masses. Unknown charge zero, nonpositive neutral mass, no
precursor or multiple precursors disable that constraint. The comparison uses
atomic mass units and the proton constant, retaining the fixes covered by the
source's unknown-charge and unequal-charge regressions.

Single-charge conversion is `mz * charge - (charge - 1) * PROTON_MASS_U`.
Optional summed intensity uses the original observed cluster intensities.
Accepted monoisotopic peaks and any retained unassigned peaks are selected and
sorted again; every existing annotation array follows both permutations.

`deisotope` returns a `DeisotopingResult` containing the new spectrum and accepted
`IsotopeCluster` values. Cluster indices refer to the original input, with
clusters ordered before charge conversion. `filter_spectrum` replaces its
argument only after success. The inherited experiment filter also commits
atomically and applies to spectra only.

Optional integer arrays are `charge`, `iso_peak_count` and `feature_number`.
Unassigned peaks receive 0, 1 and -1 respectively. Isotope cluster IDs start at
zero. Enabling an annotation whose name already exists in any input array is
an error. A valid empty input remains unchanged without adding empty arrays.

## Deliberate corrections and limits

- Peaks cannot be reused within a ladder. The simple method defaults to disjoint
  clusters; `allow_shared_isotopes: true` additionally supports the source's
  cross-cluster reuse of heavy-isotope peaks. Assigned seeds remain excluded.
  With sharing enabled, a peak appears in each returned membership and contributes
  to each cluster's optional intensity sum. The separate Poisson/KL method
  defaults to source sharing, which is necessary to reproduce its real fixture.
  Neither method reuses its seed when a tolerance window spans the isotope spacing.
- Isotope count is assigned only for an accepted ladder. The source updates
  count during failed hypotheses, which can label an unassigned peak with a
  partial count greater than one.
- Feature numbers are copied into output before selection. The source swaps
  its feature vector into the annotation and then continues indexing the
  emptied vector; the Rust implementation keeps discovery state valid.
- The source header claims zero-intensity removal for both methods, but its
  simple implementation does not remove unassigned zeros. This port follows
  the implemented behavior. Use an explicit filter to remove them if required.
- Tolerance must be finite, nonnegative and at most 100 ppm or 0.1 Da. Charges
  are positive `u8` values, minimum isotope count is at least two and cannot
  exceed maximum count. Negative peak values, invalid arrays, annotation
  collisions, numeric overflow and excessive work are errors without mutation.
- The configurable default work limit is 10,000,000 attempted hypotheses and
  isotope lookups combined. Storage is linear in input size and returned cluster
  memberships; shared memberships are bounded by the work budget. No C++ state
  is retained.

This heuristic assumes a C13-spaced ladder. It does not establish elemental
composition or identify a peptide. The Poisson/KL method adds model scoring
and top-N preprocessing through its separate options. Negative-ion deisotoping
and specialized NuXL/oligonucleotide models remain outside these APIs.

`tests/deisotoping.rs` covers source precursor regressions, charge priority,
missing/decreasing ladders, exact tolerance endpoints, annotations, conversion,
intensity conservation, disjoint membership, malformed/resource failures and
empty/zero behavior. `tests/workflows.rs` independently generates a modified
peptide's coarse fragment envelopes, collapses them, and aligns them to the
monoisotopic theoretical spectrum while accounting for the source's neutral-H
versus proton mass convention. The C++ code was inspected, not executed.
