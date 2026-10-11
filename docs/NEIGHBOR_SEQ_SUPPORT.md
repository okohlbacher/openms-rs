# NeighborSeq support

Source: `src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h` and
`src/openms/source/ANALYSIS/ID/NeighborSeq.cpp` at core `bc9cc12`.
Rust: `src/analysis/id/neighbor_seq.rs`. Tests: `tests/neighbor_seq.rs`.
Consumer: the TOPP tool `DecoyDatabase` (subset-neighbor search).

## API mapping

| C++ | Rust | Notes |
|---|---|---|
| `NeighborSeq(std::vector<AASequence>&&)` | `NeighborSeq::new(Vec<AASequence>) -> Result<Self>` | by-value ownership, same contract; errors for a non-`X` peptide without a known mass |
| `MSSpectrum generateSpectrum(const AASequence&)` | `generate_spectrum(&self, &AASequence) -> Result<MSSpectrum>` | b/y, charge 1, b1 included |
| `static bool isNeighborSpectrum(...)` | `NeighborSeq::is_neighbor_spectrum(...) -> bool` | |
| `static int computeSharedIonCount(...)` | `NeighborSeq::compute_shared_ion_count(...) -> usize` | |
| `bool isNeighborPeptide(...)` | `is_neighbor_peptide(&mut self, ...) -> Result<bool>` | |
| `struct NeighborStats` (four counters, `total`, `percentOfTotal_`, `unfindable`, `noNB`, `oneNB`, `multiNB`) | `NeighborStats` (fields, `total`, `percent_of_total`, `unfindable`, `no_nb`, `one_nb`, `multi_nb`) | counters are `usize` |
| `NeighborStats getNeighborStats() const` | `neighbor_stats(&self)` | |
| protected `createMassLookup_()` | private, inside `new` | sorted `(mass, positions)` vector in place of `std::map<double, vector<int>>` |
| protected `findCandidatePositions_(...)` | private `find_candidate_positions` | index range in place of an iterator pair |
| private `digested_relevant_peptides_` | private field, read through `relevant_peptides()` | |
| private `mass_position_map_`, `spec_gen_`, `x_residue_`, `neighbor_stats_` | private fields (`x_residue_` becomes a residue-letter test) | `-1` becomes `None` |

## Preserved source conventions

- Neighbors are decided by the Dice fraction `2 * shared / (n1 + n2)` being
  **strictly** greater than the threshold; bins are `trunc(mz / bin)`, merged
  over the two sorted peak lists exactly as the source's set-intersection
  lambda does (a matched peak of the first spectrum is consumed).
- The precursor window is closed at both ends (`lower_bound`/`upper_bound`).
- The ppm window calls `ppm_to_mass(mono_weight, tolerance)` with the source's
  swapped arguments (declared `(ppm, mz_ref)`), which gives the same product
  up to rounding; kept for bit-identical windows.
- `isNeighborPeptide` keeps iterating after the first match so every matching
  relevant peptide's counter is incremented.
- The constructor always logs `Skipped k/n peptides with unknown('X') amino
  acids.` as a warning, also when `k` is zero.

## Native differences

- Errors instead of undefined behaviour: a negative or non-finite tolerance
  is `Error::InvalidValue` (the source only `assert`s in debug builds); a
  candidate or relevant peptide without a known mass (for example `B`/`Z`) is an
  error, where the source takes a mass from its residue table.
- Counters of `is_neighbor_peptide` are updated only after all comparisons of
  the call succeeded, so an error leaves them unchanged.
- The bin index is an `i64` with Rust's saturating float conversion; the
  source's `int` conversion is undefined outside its range and identical inside.
- `percent_of_total` computes in 128-bit integers, so `count * 100` cannot
  overflow as the source's `int` product can.
- Theoretical ion masses come from the native `TheoreticalSpectrumGenerator`,
  which sums fragment formula masses; they agree with the source's residue
  sums to far below any bin width used here.

## Checked boundaries and evidence

- Tier 3: every section of `NeighborSeq_test.cpp` (lines 26-245) is
  transcribed in `tests/neighbor_seq.rs`: spectrum m/z literals
  (`TEST_REAL_SIMILAR`, relative 1e-5), the neighbor-spectrum and shared-count
  tables at bin sizes 2, 1, 0.1 and 0.05, the `isNeighborPeptide` scenario with
  its stats before and after the repeated query, ownership after move, and the
  percentage formatters including the zero-total regression.
- Tier 4: `X` peptides are unfindable; ppm and Da windows just below and above
  a 0.036 Da difference; zero tolerance on an identical mass; refused invalid
  tolerances leave counters unchanged; empty spectra are never neighbors.
- No C++ oracle was executed for this header.
