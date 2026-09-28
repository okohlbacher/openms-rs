# CometModification support

Source header, core pin `bc9cc12`: `src/openms/include/OpenMS/ANALYSIS/ID/CometModification.h`
(`CometModification.cpp`). Rust: `src/analysis/id/comet_modification.rs`.
Tests: `tests/comet_modification.rs`. Provenance: `tests/data/comet_modification_provenance.json`.

The header has no direct TOPP consumer in the ledger; upstream it serves the
CometAdapter's `variable_modXX` generation.

## API mapping

| Source | Rust |
|---|---|
| `struct CometModification` | `pub struct CometModification` (`Clone`, `Debug`, `PartialEq`) |
| `static constexpr double MASS_TOLERANCE = 1e-6` | `CometModification::MASS_TOLERANCE` |
| `double mass{0.0}` | `mass: f64` |
| `std::string residues` | `residues: String` |
| `int binary_group{0}` | `binary_group: i32` |
| `int max_mods_per_peptide{5}` | `max_mods_per_peptide: i32` |
| `int term_distance{-1}` | `term_distance: i32` |
| `int nc_term{0}` | `nc_term: i32` |
| `bool required{false}` | `required: bool` |
| `CometModification()` | `Default` |
| `CometModification(const ResidueModification*, int, int)` | `CometModification::from_modification(&ResidueModification, i32, i32)` |
| `bool isMergeableWith(const CometModification&) const` | `is_mergeable_with(&self, &Self) -> bool` |
| `void merge(const CometModification&)` | `merge(&mut self, &Self)` |
| `std::string toCometString(Size) const` | `to_comet_string(&self, usize) -> String` |
| `static std::vector<CometModification> mergeModifications(const std::vector<CometModification>&)` | `merge_modifications(&[Self]) -> Vec<Self>` |

## Preserved source conventions

- The merge rules of the header brief, as the `.cpp` implements them, including
  the asymmetric residue order of a merge (`K` + `n` gives `Kn`).
- A terminal and an amino-acid entry merging resets `term_distance`/`nc_term`
  to `-1`/`0`; two terminal or two amino-acid entries keep them.
- `to_comet_string`: `variable_mod` + `setw(2)`/`setfill('0')` index, the mass
  in the `std::ostream` default (six significant digits, via
  `concept::fuzzy_string_comparator::format_g`), integers in decimal,
  `required` as `1`/`0`, and a fixed `0.0` neutral-loss column.
- `merge_modifications` is greedy first-fit in input order and never drops an
  entry.
- The mass test is the source's `|a - b| > MASS_TOLERANCE` refusal, NaN
  included.

## Native differences

- The source constructor takes a pointer and does not check it for null; the
  port takes a reference.
- This crate's `ModificationsDB` stores the origin of an any-residue terminal
  modification as `None` where the source stores `'X'`; `from_modification`
  treats both as `'X'`, so the residue string becomes `n`/`c` exactly as upstream.

## Source defect candidate

After a peptide N-term entry has merged with an amino-acid entry it has
`term_distance == -1`, so `isMergeableWith` no longer sees it as peptide N-term,
and a later peptide C-term entry of the same mass merges into it: `n`, `K`, `c`
of one mass collapse into `nKc`, one line that allows the modification at both
termini and on K. The header states that peptide N-term and C-term entries never
merge. The port keeps the source behaviour; the test
`a_merged_terminal_entry_accepts_the_opposite_terminal` pins it.

## Checked boundaries and evidence

Tier 3: every section and literal of `CometModification_test.cpp` at `bc9cc12`
(construction from `Oxidation (M)`, `Acetyl (N-term)`, `Amidated (C-term)`; the
twelve mergeability cases; the nine merge cases; the five string cases; the
five list-merge cases). Tier 4: default members, protein-terminal construction,
negative-mass and zero-mass formatting and three-digit indices, and the
merged-terminal behaviour above.

No OpenMP in the source.
