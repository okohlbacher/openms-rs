# Assays and decoys (MRMIonSeries, MRMDecoy, MRMAssay)

Ports of the OpenMS core pin `bc9cc12` headers

| Header | Rust | Status |
|---|---|---|
| `ANALYSIS/OPENSWATH/MRMIonSeries.h` (+ `.cpp`) | `src/analysis/openswath/mrm_ion_series.rs` | complete |
| `ANALYSIS/OPENSWATH/MRMDecoy.h` (+ `.cpp`) | `src/analysis/openswath/mrm_decoy.rs` | partial: the `LightTargetedExperiment` members wait for `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h` |
| `ANALYSIS/OPENSWATH/MRMAssay.h` (+ `.cpp`) | `src/analysis/openswath/mrm_assay.rs` | partial: the `LightTargetedExperiment` members wait for `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h` |
| `ANALYSIS/OPENSWATH/TransitionListEvidenceFilter.h` (+ `.cpp`) | not ported | untouched: needs `OpenSwath::SwathMap`, `LightTargetedExperiment`, `ChromExtractParams` (`OpenSwathWorkflow.h`) and the OpenSWATH spectrum access layer |

Tests: `tests/mrm_ion_series.rs`, `tests/mrm_decoy.rs`, `tests/mrm_assay.rs`,
`tests/mrm_assay_oracle.rs`. Provenance: `tests/data/mrm_assay_provenance.json`.
Shared code touched: `src/chemistry/theoretical.rs` (`residue_loss_formulas`,
the one residue neutral-loss table, now used by both the theoretical spectrum
generator and the ion series) and `src/concept/math_functions/random_shuffler.rs` (the
draw-limited, generic `shuffle_limited` of `Math::RandomShuffler`, used over indices).

The TOPP consumers `OpenSwathAssayGenerator` and `OpenSwathDecoyGenerator` read and
write TraML, transition TSV and PQP (`TraMLFile`, `TransitionTSVFile`,
`TransitionPQPFile`), which are backlog items that depend on the types ported here;
no format is re-scoped out.

## API mapping

### MRMIonSeries

| Source | Rust |
|---|---|
| `MRMIonSeries()`, `~MRMIonSeries()` | `MRMIonSeries` (unit struct, `Default`); the members are associated functions |
| `typedef std::map<std::string, double> IonSeries` | `IonSeries = BTreeMap<String, f64>` (same key order) |
| `getIon(IonSeries&, ionid)` | `MRMIonSeries::get_ion(&IonSeries, &str) -> (String, f64)` |
| `annotateIon(ionseries, product_mz, mz_threshold)` | `MRMIonSeries::annotate_ion` |
| `annotateTransitionCV(tr, annotation)` | `MRMIonSeries::annotate_transition_cv` (native parse) and `annotate_transition_cv_with(.., &FormulaLossParsing)` |
| `annotateTransition(tr, peptide, prec_thr, prod_thr, reannotation, types, charges, specific, unspecific, round_decPow)` | `MRMIonSeries::annotate_transition(tr, peptide, prec_thr, prod_thr, reannotation, &IonSeriesOptions)` |
| `getIonSeries(sequence, precursor_charge, types, charges, specific, unspecific, round_decPow)` | `MRMIonSeries::get_ion_series(&AASequence, usize, &IonSeriesOptions) -> Result<IonSeries>` |
| private `annotationToCVTermList_`, `annotationToCV_` | private `annotation_to_interpretation`, `annotation_to_cv` |
| anonymous `extractChargeFromAnnotation_`, `extractOrdinalFromAnnotation_` | private `charge_of_annotation`, `ordinal_of_annotation` |
| the separate `fragment_types`, `fragment_charges`, loss switches and `round_decPow` arguments | `IonSeriesOptions { fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, round_dec_pow, unresolved_modifications, max_ion_series_entries, formula_loss }` |
| the `static const EmpiricalFormula nl_formula` of `annotationToCVTermList_` | `FormulaLossParsing::{Parse, SourceStatic}` |

### MRMDecoy

| Source | Rust |
|---|---|
| `MRMDecoy()` (`DefaultParamHandler("MRMDecoy")`) | `MRMDecoy::new` / `Default`; `MRM_DECOY_NAME` |
| parameters `non_shuffle_pattern`, `keepPeptideNTerm`, `keepPeptideCTerm`; `updateMembers_` | `defaults`, `to_param`, `from_param`, `from_param_with_warnings`; fields `keep_const_pattern`, `keep_n`, `keep_c` |
| `ProgressLogger` base | not ported (no progress output) |
| `IndexType` | `IndexType = Vec<usize>` |
| `ProteinVectorType`, `PeptideVectorType`, `TransitionVectorType`, `PeptideTransitionMapType` | not needed: `Vec<Protein>`, `Vec<Peptide>`, `Vec<ReactionMonitoringTransition>`, `BTreeMap<&str, Vec<&ReactionMonitoringTransition>>` internally |
| `generateDecoys(exp, dec, method, fraction, switchKR, tag, attempts, identity, prec_shift, prod_shift, prod_thr, types, charges, specific, unspecific, round_decPow)` | `MRMDecoy::generate_decoys(&self, &exp, &mut dec, &DecoyOptions)`; `method` is `DecoyMethod` (`from_source_str`) |
| `generateDecoysLight` | not ported (`LightTargetedExperiment`) |
| `switchKR(peptide)` | `MRMDecoy::switch_kr(&mut Peptide)` |
| `AASequenceIdentity(sequence, decoy)` | `MRMDecoy::aa_sequence_identity -> Result<f32>` |
| `shufflePeptide(peptide, identity, seed = -1, attempts = 100)` | `MRMDecoy::shuffle_peptide(&self, &Peptide, f64, u32, usize) -> Result<Peptide>`; the time seed is `DecoyOptions::shuffle_seed = None` |
| static `reversePeptide(peptide, keepN, keepC, pattern)` | `MRMDecoy::reverse_peptide -> Result<Peptide>` |
| static `findFixedResidues(sequence, keepN, keepC, pattern)` | `MRMDecoy::find_fixed_residues` |
| static `reversePeptideLight`, `shufflePeptideLight`, `switchKRLight` | not ported as such (`LightModification`); the heavy methods run the same string-level algorithms |
| protected `hasCNterminalMods_` | `MRMDecoy::has_cn_terminal_mods` (public) |
| protected static `hasCNterminalModsLight_`, `pseudoreversePeptideLight_` | not ported (`LightModification`) |
| protected `findFixedResidues_`, `findFixedAndTermResidues_` | `fixed_residues`, `fixed_and_terminal_residues` (public) |
| protected `pseudoreversePeptide_`, `reversePeptide_` | `pseudoreverse_peptide`, `reverse_peptide_full` (public) |
| protected `getModifiedPeptideSequence_` | `modified_peptide_sequence` (public) |
| members `keep_const_pattern_`, `keepN_`, `keepC_` | public fields `keep_const_pattern`, `keep_n`, `keep_c` |
| `Math::RandomShuffler` (default seed) in `generateDecoys` | `DecoyRandom::seeded(5489).shuffle_slice` |
| `boost::mt19937` + `boost::uniform_int<>` + `variate_generator` | private `BoostUniformInt` over `rand_mt::Mt` (crate-internal, shared with `MRMAssay`) |

### MRMAssay

| Source | Rust |
|---|---|
| `MRMAssay()`, `~MRMAssay()`, `ProgressLogger` base | `MRMAssay` (unit struct); progress not reported |
| `ProteinVectorType`, `PeptideVectorType`, `CompoundVectorType`, `TransitionVectorType`, `PeptideTransitionMapType`, `CompoundTransitionMapType` | not needed (plain `Vec`s and a private `group_by` into `BTreeMap`) |
| `ModifiedSequenceMap`, `SequenceMapT`, `IonMapT`, `PeptideMapT`, `TargetDecoyMapT` | private type aliases over `BTreeMap` (swath keys as `usize`, so -1 sorts last as the source's `size_t`) |
| `FragmentSeqMap`, `IonSeries` | `FragmentSeqMap = Vec<(f64, String)>`, `AssayIonSeries = Vec<(String, f64)>` |
| `reannotateTransitions(exp, prec_thr, prod_thr, types, charges, specific, unspecific, round_decPow)` | `MRMAssay::reannotate_transitions(&mut exp, prec_thr, prod_thr, &IonSeriesOptions) -> Result<()>` |
| `restrictTransitions(exp, lower, upper, swathes)` | `MRMAssay::restrict_transitions(&mut exp, lower, upper, &[(f64, f64)], &IonSeriesOptions) -> Result<()>` (the options carry only the peptide-conversion policy) |
| `detectingTransitions(exp, min, max)` | `MRMAssay::detecting_transitions(&mut exp, i32, i32)` |
| `uisTransitions(exp, types, charges, specific, unspecific, ms2, mz_threshold, swathes, round_decPow, max_alt, seed, disable_decoys)` | `MRMAssay::uis_transitions(&mut exp, &UisOptions) -> Result<()>` |
| `filterMinMaxTransitionsCompound(exp, min, max)` | `MRMAssay::filter_min_max_transitions_compound` |
| `filterUnreferencedDecoysCompound(exp)` | `MRMAssay::filter_unreferenced_decoys_compound -> Result<()>` |
| `reannotateTransitionsLight`, `restrictTransitionsLight`, `detectingTransitionsLight`, `uisTransitionsLight` | not ported (`LightTargetedExperiment`) |
| protected `getMatchingPeptidoforms_` | `MRMAssay::matching_peptidoforms` (public) |
| protected `getSwath_`, `isInSwath_` | `MRMAssay::swath_index`, `MRMAssay::is_in_swath` (public) |
| protected `getRandomSequence_(size, pseudoRNG)` | crate-internal `random_sequence`; public `random_sequence_with_seed(size, seed)` |
| protected `nchoosekcombinations_` | `MRMAssay::n_choose_k_combinations -> Result` |
| protected `addModificationsSequences_` | `MRMAssay::add_modifications_sequences -> Result` |
| protected `generateTheoreticalPeptidoforms_`, `generateTheoreticalPeptidoformsDecoy_` | `MRMAssay::theoretical_peptidoforms`, `theoretical_peptidoforms_decoy` |
| protected `generateTargetInSilicoMap_`, `generateDecoySequences_`, `generateDecoyInSilicoMap_`, `generateTargetAssays_`, `generateDecoyAssays_` | private functions of the same names in snake case |
| protected Light helpers (`generateTargetInSilicoMapLight_` and three more), `TargetDecoyMapLightT` | not ported (`LightTargetedExperiment`) |
| anonymous `sortIonMap`, `sortAndDeduplicatePeptideMap` | private `sort_ion_map`, `sort_and_dedup_peptide_map` |

## Preserved source conventions

- **Residue neutral losses** (`ResidueDB.cpp`): D, E, S, T lose `H2O1`; K, N, Q
  lose `H3N1`; R loses `H3N1`, `C1H2N2` and `C1H2N1O1`. A modified residue carries
  only its modification's losses (for example Phospho on S/T `H3O4P1`, Oxidation on
  M `C1H4O1S1`). The four unspecific losses are `H2O1`, `H3N1`, `C1H2N2`,
  `C1H2N1O1`; every other declared loss is specific.
- **Ion-series keys** are `<type><ordinal>^<charge>` and
  `<type><ordinal>-<formula>^<charge>`; m/z is rounded with `Math::roundDecimal`.
- **Tie rules** of `annotateIon`: closer by more than 1% of the threshold wins,
  otherwise lower charge, then lower ordinal; keys are visited in `std::map` order.
- **Random numbers**: `boost::mt19937(seed)` words (the reference MT19937 of
  `rand_mt`) mapped by Boost's bucket division. `boost::uniform_int<>` defaults to
  `[0, 9]`, so the source's `pseudoRNG() % 17` only ever picks one of the first ten
  residues `ANDCEQGHIL`, and `pseudoRNG() % len` one of the first ten positions;
  both reproduced (the class-test random sequence `CHLNHHQQNE` and every seeded
  shuffle literal match).
- **switchKR** of a non-K/R terminus: FNV-1a over the sequence bytes modulo 17
  (bytes above 127 sign-extended as on x86-64).
- **Map orders**: transitions grouped by peptide or compound reference in
  ascending byte order; UIS swath index -1 (`size_t` max) iterates last.
- **Name formats** of UIS transitions:
  `<index>_UIS_{<forms>}_<precursor>_<product>_<RT>_<annotation>` with
  `StringUtils::toStr(double)`; decoy ones `UISDECOY`, and the decoy index is not
  advanced for a skipped overlapping transition.
- **`a`/`b`/`c` ions from the suffix** in `annotateTransition`, the **added** mass of
  a non-named negative loss, and the default rounding `-4` of the reannotation and
  in silico ion series, all as the source (documented at the items).
- **`int` limits cast to `size_t`** in `detectingTransitions` and
  `filterMinMaxTransitionsCompound` (negative minimum keeps nothing, negative
  maximum is no limit).
- **`n choose k` with `k > n`** yields one combination of all positions (the
  source's bitmask truncation), so such a peptidoform carries fewer modifications.
- **Decoy method `shift`** leaves the sequence unchanged, so each decoy precursor
  equals its target and is excluded, as in the source.

## Native differences

- **Formula losses are parsed per annotation** (`FormulaLossParsing::Parse`). The
  source parses into a function-local `static`, so the first formula loss a process
  parses fixes the mass of every later one: the Release build reports
  `y7-H2O1^2` with the ammonia mass -17.0265490957 when `H3N1` came first.
  `FormulaLossParsing::SourceStatic` reproduces it with a caller-owned cache; the
  oracle test uses it.
- **Unresolvable modifications** are refused by default
  (`UnresolvedModificationPolicy::Reject`); `SourceSkip` reproduces the source.
- **Seeds are explicit**: the source seeds every `generateDecoys` shuffle and the
  UIS decoys from `time(nullptr)` when given -1; `None` does the same here, `Some`
  makes the result reproducible.
- **Unknown decoy method strings are errors** (`DecoyMethod::from_source_str`); the
  source silently leaves each sequence unchanged, so every decoy equals its
  target and is excluded.
- **Checked instead of undefined**: charge 0 (division), ordinals above 255
  (`unsigned char` truncation), unequal lengths in `AASequenceIdentity`, repeated
  pattern letters that run the fixed-position erase/insert out of range, positions
  past the C-terminus in `addModificationsSequences_`, and an id containing `decoy`
  without `_decoy` (`std::out_of_range`) are `Error::InvalidValue`.
- **Resource ceilings** (checked before work starts): ion-series entries
  (`DEFAULT_MAX_ION_SERIES_ENTRIES`, 2^22), shuffle work (`MAX_SHUFFLE_WORK`),
  decoy transitions (`MAX_DECOY_TRANSITIONS`), modification combinations
  (`MAX_PEPTIDOFORM_COMBINATIONS`, 2^20; the source enumerates all of them before
  comparing with `max_num_alternative_localizations`) and UIS transitions
  (`MAX_UIS_TRANSITIONS`).
- **Atomicity**: every mutating function builds into temporaries and commits at the
  end; on error its target is unchanged.
- The protected helpers the class tests exercise are public here.

## Checked boundaries and evidence

- **Tier 1** (`tests/mrm_assay_oracle.rs`, 8 tests): oracle
  `../oracle/O1c-assays` on ibminode06 against the Linux x86-64 **Release** install
  `openms4-release-bc9cc12-c19e494-174b576`, every case twice in fresh processes and
  the batch twice, byte-identical. Inputs are the class-test TraML files, dumped by
  the oracle (TraML is not ported) and rebuilt by the test. Covered:
  `reannotateTransitions` (both class-test sets), `restrictTransitions` (with and
  without swathes), `detectingTransitions`, `filterMinMaxTransitionsCompound`,
  `filterUnreferencedDecoysCompound`, `uisTransitions` (all five class-test sets),
  `generateDecoys` (pseudo-reverse, reverse, fraction 0.5, switchKR). All fields
  compare exactly except the unrounded decoy precursor m/z.
- **Declared tolerance (a finding)**: m/z within 4 ulp relative. Only the decoy
  precursor m/z needs it: `Residue::getMonoWeight(Internal)` is the full-residue
  formula mass (summed in `const Element*` map order) minus water, while the ported
  `AASequence` sums the internal composition in symbol order, so some residue masses
  (F, D, K among them) differ in the last bit; the `decoy_masses` oracle case pins
  it. This belongs to the `AASequence`/`Residue` port, not to this group.
- **Tier 3**: every literal of `MRMIonSeries_test.cpp`, `MRMDecoy_test.cpp` (all
  seeded shuffles, reversals, identity, switchKR, terminal checks, the heavy
  duplicate-detection regression) and `MRMAssay_test.cpp` (matching peptidoforms,
  swath lookup, random sequence, n-choose-k, addModificationsSequences, both
  peptidoform enumerations).
- **Tier 4**: native error paths, atomicity, parameters, ceilings.
- The source contains no `#pragma omp`; the port is serial.
