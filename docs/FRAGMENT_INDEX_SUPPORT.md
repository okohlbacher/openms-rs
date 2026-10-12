# FragmentIndex support

Source: `src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h` and
`src/openms/source/ANALYSIS/ID/FragmentIndex.cpp` at core `bc9cc12`.
Rust: `src/analysis/id/fragment_index.rs` (index, build, query, SNES build
and realization) and `src/analysis/id/fragment_index_parameters.rs` (the
`DefaultParamHandler` half). Tests: `tests/fragment_index.rs`. No direct TOPP
consumer at the pin (ProSE and its TOPP tool use it).

Status: **partial**. The one unported member is the SNES spectrum query
(`querySpectrumSNES_`, source lines 1667-2276): `query_spectrum` returns
`Error::Unsupported` in SNES mode. Everything else below is ported.

## API mapping

| C++ | Rust | Notes |
|---|---|---|
| `struct Peptide` (ctor, `protein_idx`, `mod_bitmask_`, `sequence_`, `precursor_mz_`) | `Peptide` (`new`, `protein_idx`, `mod_bitmask`, `sequence`, `precursor_mz`) | |
| `struct SpectrumMatch` | `SpectrumMatch` | `Default` is all-zero, as the source's value-initialisation |
| `struct SpectrumMatchesTopN` (`hits_`, `operator+=`, `clear`) | `SpectrumMatchesTopN` (`hits`, `AddAssign<&Self>`, `clear`) | |
| `FragmentIndex()` + `setParameters` | `FragmentIndex::new(FragmentIndexParameters)`, `from_param(&Param)` | |
| `setParameters` on an existing object | `set_parameters(FragmentIndexParameters)` | keeps a built index, as the source |
| `getParameters` / defaults | `parameters()`, `FragmentIndexParameters::{defaults, to_param, from_param, from_param_with_warnings, validate}` | |
| `isBuild`, `getPeptides`, `getNumFragments` | `is_build`, `peptides`, `num_fragments` | |
| `addSpecialPeptide` (`DEBUG_FRAGMENT_INDEX` only) | not ported | debug-build-only API |
| `build(fasta_entries)` | `build(&mut self, &[FASTAEntry]) -> Result<()>` | |
| `clear` | `clear` | |
| `getPeptidesInMassWindow` | `get_peptides_in_mass_window` | |
| static `isOpenSearchMode` | `FragmentIndex::is_open_search_mode` | |
| `SNES_KIND_BIT_MASK`, `SNES_SLOT_MASK`, `enum SnesAnchor`, `isSingleCMother`, `isSingleNMother`, `isSnesMode` | same names in snake case | |
| `struct Hit`, `query(peak, range, charge)` | `Hit`, `query(&self, &Peak1D, (usize, usize), u16)` | `&self`: the source's non-const has no reason |
| `querySpectrum(spectrum, sms)` | `query_spectrum(&self, &MSSpectrum, &mut SpectrumMatchesTopN) -> Result<()>` | SNES mode: `Error::Unsupported` |
| `querySpectrum(spectrum, fasta, sms)` | `query_spectrum_with_fasta` | fasta only needed by the unported SNES query |
| `reconstructModifiedSequence` | `reconstruct_modified_sequence -> Result<AASequence>` | |
| `realizeSNESLength` (`-1` = none) | `realize_snes_length -> Result<Option<usize>>` | |
| `reconstructRealizedSubSequence` | `reconstruct_realized_sub_sequence -> Result<AASequence>` | |
| protected `Fragment`, `fi_fragments_`, `bucketsize_` | `Fragment`, `fragments()`, `bucket_size()` | read-only exposure for the class-test invariants |
| protected `updateMembers_` | `FragmentIndexParameters::from_param` + `FragmentIndex::new` | |
| protected `generatePeptides`, `generateSNESMothers_`, `generateFragmentsLightweight_`, `generateFragmentsForSeries_`, `initModificationTables_`, `initResidueMassTable_` | private | |
| protected `buildModSlots_`, `ModSlot`, `MAX_MOD_SLOTS` | `build_mod_slots -> Vec<ModSlot>`, `ModSlot`, `MAX_MOD_SLOTS` | public for inspection |
| protected `computeSnesSigmaDeltaSet_`, `snes_sigma_delta_set_*` | `compute_snes_sigma_delta_set`, `snes_sigma_delta_set(SnesAnchor)` | |
| private `queryPeaks`, `searchDifferentPrecursorRanges`, `computeMassWindow_`, `isOpenSearchMode_` | private `query_peaks`, `search_different_precursor_ranges`; public `compute_mass_window`, `FragmentIndexParameters::is_open_search_mode` | |
| private `trimHits` | `trim_hits` | public |
| private `querySpectrumSNES_` | **not ported** | the gap |
| protected member variables | fields of `FragmentIndexParameters` and private state | integer parameters (`fragment:min_mz`, `peptide:min_mass`, …) stay `i32`, converted to `f32` where the source converts |

## Preserved source conventions

- Every parameter, default, description (including the truncated
  `"Unit of fragment m"` and `"Add peaks of  x-ions"`), restriction and
  section of the constructor, including the declared-only `decoys`,
  `annotate:PSM`, `report:top_hits`, `peptide:motif` and `precursor:isotopes`.
- Validation of `updateMembers_`: finite tolerances, a positive window
  width, SNES needs both b and y ions; open-search mode is logged.
- Peptides: digestion with X/B/Z products skipped and counted; mass = water +
  proton + fixed terminal deltas + residue internal masses + fixed residue
  deltas; variable modifications as bitmasks over at most 32 slots (N-terminal,
  residues left to right, C-terminal), skipping subsets over the maximum or
  with two slots on one position; peptides outside `[min_mass, max_mass]` (in
  `f32`) dropped; sort by `(precursor_mz, protein_idx)`.
- SNES mothers per unambiguous span, Single-N and Single-C, `max_size` 0 =
  unlimited, only the lower mass bound, b-only / y-only fragments with the
  free terminus' fixed modification omitted.
- Fragments at charge 1 from cumulative `f64` sums narrowed to `f32`, ion
  offsets from `Residue::getInternalTo*Ion`, ions with index at most
  `min_ion_index` skipped, `[min_mz, max_mz]` inclusive; sorted by `(mz,
  peptide)`, buckets of 4096 recording their first m/z, each bucket sorted by
  peptide (stable, as the source's radix sort).
- Queries in the source's float arithmetic: peak mass `mz * z - (z - 1) *
  PROTON_MASS_U` computed in `f64` and narrowed; ppm tolerances as
  `Math::ppmToMass<float>`; bucket walk from one bucket before `lower_bound`;
  half-open peptide ranges; counts per (charge, isotope) block emitted in
  ascending peptide order at `max(min_matched_ions, 1)`; isotope iteration
  collapsed to 0 in open search; unknown precursor charge tries
  `min_charge..=max_charge`; `trimHits` ranking (matches, |isotope|, isotope,
  charge, peptide) and threshold.
- `reconstructModifiedSequence` applies a fixed terminal modification that has
  a residue origin to every peptide terminus, whatever the residue, as the
  source does.

## Native differences

- Errors instead of undefined behaviour: FASTA entries with characters other
  than `A`-`Z` are refused (the source indexes a 128-entry table with them);
  entries over 65535 residues are refused as in the source; out-of-range
  peptide coordinates in the reconstruction functions are errors; integer
  parameters outside their member type (negative counts, charges above 65535,
  isotope errors outside `i16`) are refused where the source wraps; a
  negative or non-finite fragment tolerance and an unknown enzyme are refused
  at construction.
- Peptides with equal `(precursor_mz, protein_idx)` keep generation order (a
  stable sort); the source's `std::sort` leaves it to the library.
- With 32 variable slots the source's `1u << 32` is undefined; the port
  enumerates all 2^32 subsets bounded by the other limits.
- `peptide:min_size` 0 is passed to the digestion as 1, which selects the same
  (never empty) products; the native digestion refuses 0.
- Match counting uses an ordered map in place of the source's thread-local
  window tables; results are identical and there is no cross-call state.
- `build` runs peptide generation and fragment generation on `threads`
  workers; the index is identical for every count
  (`build_is_identical_across_thread_counts`).

## Checked boundaries and evidence

- Tier 3, `FragmentIndex_test.cpp` transcribed in `tests/fragment_index.rs`
  (39 tests): `build()` expected peptide sets for four parameterisations and
  both ordering invariants (lines 178-262); enzyme specificity counts
  (270-375); the 65535 guard (379-392); `clear` (395-406); self-hits at charges
  1-4 with and without known charge (410-439); isotope errors -3..3 (443-499);
  tolerance jitter (503-569); fragment counts with and without
  `min_ion_index` (573-611); modification slots (615-695); cross-validation of
  six peptides against `ModifiedPeptideGenerator` (700-845); symmetric and
  asymmetric windows, open mode, isotope iteration (849-1041); validation
  (1045-1086); the half-open range contract and its regression guard
  (1089-1261); SNES enumeration, dispatch, realization in both directions
  with asymmetric tolerances, size, fixed and subset modifications, Σ sets,
  the kind-bit mask, the b/y requirement and X/B/Z spans (1267-1541,
  1851-2091, 2609-2744).
- Not transcribed, because they exercise the unported SNES query: lines
  1544-1692, 1779-1848, 2094-2606, 2747-2879. The non-SNES determinism section
  (1693-1776) tests thread-local buffer reuse, which has no Rust counterpart.
- Tier 4: SNES query reported unsupported, parameter round trip, Σ sets and
  slot order, thread-count invariance.
- No C++ oracle was executed for this header.
