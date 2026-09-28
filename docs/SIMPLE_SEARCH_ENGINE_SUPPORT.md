# SimpleSearchEngineAlgorithm support

Source: `src/openms/include/OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h`
and `src/openms/source/ANALYSIS/ID/SimpleSearchEngineAlgorithm.cpp` at core
`bc9cc12`. Rust: `src/analysis/id/simple_search_engine.rs`. Tests:
`tests/simple_search_engine.rs`. Consumer: the TOPP tool `SimpleSearchEngine`
(the tool itself is not part of this port).

## API mapping

| C++ | Rust | Notes |
|---|---|---|
| `SimpleSearchEngineAlgorithm()` + `DefaultParamHandler` | `SimpleSearchEngineAlgorithm::default`, `defaults`, `to_param`, `from_param`, `from_param_with_warnings` | the peak-picker parameter pattern |
| `enum class ExitCodes` | `ExitCodes` | `UnexpectedResult`, `UnknownError`, `IllegalParameters` are never returned: failures are `Err` |
| `ExitCodes search(in_spectra, in_db, prot_ids, pep_ids) const` | `search(&self, &str, &str, &mut Vec<ProteinIdentification>, &mut Vec<PeptideIdentification>) -> Result<ExitCodes>` | same asymmetric output contract |
| (none) | `search_experiment` | native in-memory entry point used by `search` |
| protected `updateMembers_()` | private `from_complete_param` | |
| protected `struct AnnotatedHit_` (fields, `hasBetterScore`) | `AnnotatedHit` (fields, `has_better_score`, `better_first`) | `sequence` is an `Arc<str>`, not a view; the padding layout has no Rust meaning |
| protected static `preprocessSpectra_(PeakMap&, double, bool)` | `preprocess_spectra(&mut MSExperiment, f64, bool) -> Result<()>` | public |
| protected `postProcessHits_(exp, annotated_hits, protein_ids, peptide_ids, top_hits, fixed, variable, max_var, mods_fixed, mods_variable, missed_cleavages, prec_tol, frag_tol, prec_unit, frag_unit, min_z, max_z, enzyme, database_name) const` | `post_process_hits(&self, exp, annotated_hits, protein_ids, peptide_ids, database_name) -> Result<()>` | the member-copy arguments are read from `self` |
| protected members `precursor_mass_tolerance_` … `report_top_hits_` | public fields of the same names | `*_unit_` strings become `MassToleranceUnit`; `peptide_enzyme_specificity_` is `DigestionSpecificity` |
| file-local `mapPrecursorMassesToScans` | private `map_precursor_masses_to_scans` | sorted vector in place of `std::multimap` |
| `ProgressLogger` base (`startProgress`/`setProgress`/`endProgress` in `search`) | `search_with_progress`, `search_experiment_with_progress` taking a `ProgressLogger` | the four source sections; the plain entry points report nothing |

## Preserved source conventions

- Parameter names, defaults, descriptions (including the truncated
  `"Unit of fragment m"`), restrictions and section descriptions follow the
  constructor; duplicate modifications are made unique through a sorted set,
  with the source's warning, only when duplicates exist.
- Preprocessing: `ThresholdMower` default (0.05, not zero as the source
  comment says), `Normalizer` to one, sort by RT, per spectrum sort, deisotope
  and decharge (charges 1-3, 3-10 peaks, shared isotope extensions as the
  source), jumping 100 Th `WindowMower` keeping 20, `NLargest` 400, sort.
  Deisotoping is skipped with a warning beyond 100 ppm / 0.1 Da.
- Precursor masses `z * mz - z * PROTON_MASS_U - i * C13C12_MASSDIFF_U` per
  isotope offset `i`, only for spectra with one precursor, a charge in range
  and at least `peptide:min_size` peaks; windows `m ± m * tol * 1e-6` (ppm) or
  `m ± tol`, closed at both ends.
- Candidates skip `X`/`B`/`Z`, require a full motif match, and are scored once
  per distinct peptide with ln(hyperscore) against b/y ions at charge 1
  including b1; zero scores are dropped; per-spectrum lists are pruned to the
  top N whenever they reach `2 * top_hits`.
- Annotations use a `SpectrumAlignment` with its default 0.3 Th tolerance,
  fragment charges up to `min(z - 1, 2)`, and the precursor error correction
  `isotope_error * PROTON_MASS_U / z` (not the C13 spacing the mass used).
- Decoys: `DECOY_` prefix, plain reversal for `none` specificity, reversal
  within enzymatic peptides otherwise, then targets and decoys shuffled with a
  default-seeded (5489) Mersenne twister and Boost's `uniform_int` mapping,
  through the crate's `DecoyRandom` (made generic over the element type for
  this). The protein run identifier is `SSE_<DateTime::get()>`, stamped on
  every peptide identification, including pre-existing ones.
- Peptide indexing with `DECOY_` prefix, the search enzyme and specificity and
  `missing_decoy_action` silent; with decoys and `FDR:PSM > 0`: q-values,
  score filter, decoy and empty-identification removal, unreferenced proteins
  removed. Without decoys the source's warning is logged.

## Native differences

- Errors are `Err` and leave both outputs unchanged; the source maps
  peptide-indexer failures to exit codes after writing its outputs. An empty
  database still returns `ExitCodes::InputFileEmpty` after writing, as the
  source's `DATABASE_EMPTY` path does.
- Count parameters (`Size` members) refuse negative values; the source wraps
  them.
- The per-spectrum top-N selection is a stable sort. Hits that tie on every
  key of `hasBetterScore` (one peptide matching one spectrum under two
  isotope offsets, which needs a Dalton tolerance wider than the isotope
  spacing) keep the single-threaded source's insertion order; the source's
  `partial_sort` leaves them unspecified.
- Scoring and preprocessing run on `threads` workers with results identical
  for every count (`results_are_identical_across_thread_counts`); the source's
  OpenMP result order depends on its final re-sort by scan index.
- The protein run's IM unit is read from the first spectrum only when there
  is one; the source reads `exp[0]` unguarded.
- Ordinals of ion names beyond `i32` saturate where the source's `toInt32`
  throws; matched-ion counts beyond `u16` saturate where the source's cast
  wraps.
- `getAllSearchModifications` is reproduced by a private helper; the same
  helper exists privately in `src/format/mascot_generic.rs` and in
  `fragment_index_parameters.rs` and should be collapsed onto
  `ModificationsDB` at integration.

## Checked boundaries and evidence

- Tier 1, retained C++ outputs: `TOPP_SimpleSearchEngine_1` and `_2`
  (`topp/CMakeLists.txt:3773-3789` at test-data `0cb15f2`) run through
  `search` with the tool's `Search:` parameter subtree taken from the
  retained `.ini` files, and compared as decoded idXML with the upstream
  FuzzyDiff contract (ratio 1.01 or absolute 0.01, text exact, `date` and
  `db` whitelisted): search parameters, the protein run (hits, raw and primary
  run paths by file name), every peptide identification and hit (sequence,
  charge, score, evidences, every meta value with its type, every fragment
  annotation). The largest relative deviation measured on dax is 2.1e-10
  (`fragment_mz_error_median_ppm`, `SimpleSearchEngine_2`), 6.9e-11 for
  `_1`: the port reproduces the retained output to double rounding, far inside
  the contract.
- Progress (source review, `SimpleSearchEngineAlgorithm.cpp:653-655`,
  `686-712`, `718-864`, `870-891`): the recorded sections of a decoy search are
  exactly the source's four, with the per-protein values 1..22 of a
  single-threaded run, and the result equals the silent search's. The values
  are reported after the parallel digestion rather than interleaved with it.
- Tier 4: thread-count invariance, parameter round trip and the
  `updateMembers_` duplicate handling, refused restrictions and negative
  counts, the shortlex tie-break, the decoy/FDR path keeping only target
  evidence, and a missing database reported as an error.
- `SimpleSearchEngineAlgorithm_test.cpp` has no expected values (constructor
  and `NOT_TESTABLE`).
