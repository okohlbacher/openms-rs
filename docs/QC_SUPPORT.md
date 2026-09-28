# QC metrics and FlagSet support

Source headers (core `bc9cc12`): `QC/QCBase.h`, `QC/SpectrumCount.h`,
`QC/TIC.h`, `QC/FeatureSummary.h`, `QC/MissedCleavages.h`,
`QC/IdentificationSummary.h` and `DATASTRUCTURES/FlagSet.h`.

Rust: `src/analysis/qc.rs` (module root) with `src/analysis/qc/{qc_base,
spectrum_count, tic, feature_summary, missed_cleavages,
identification_summary}.rs`, and `src/data_structures/flag_set.rs`.
Tests: `tests/qc.rs` (class tests), `tests/mzqc_file.rs` (end to end through
`MzQCFile`). Provenance: `tests/data/qc_provenance.json`.

These are the metrics `MzQCFile::store` needs (see
`docs/MZQC_FILE_SUPPORT.md`). The source keeps them in `OpenMS/QC`; they live
under `analysis` because their one format consumer sits in `format`, which
`analysis` already reaches and which already reaches `analysis`: a top-level
`qc` module would have added a module cycle (`tools/check_module_cycles.py`).

The other twelve headers of `OpenMS/QC` are not ported: `Contaminants`,
`FragmentMassError`, `FWHM`, `MQEvidenceExporter`, `MQExporterHelper`,
`MQMsmsExporter`, `Ms2IdentificationRate`, `Ms2SpectrumStats`, `MzCalibration`,
`PeptideMass`, `PSMExplainedIonCurrent`, `RTAlignment`. They serve the
QualityControl TOPP tool and each is a backlog item for it.

## API mapping

### `FlagSet<ENUM>` -> `FlagSet<E: FlagEnum>`

| Source | Rust |
|---|---|
| `FlagSet()` | `FlagSet::new`, `Default` |
| `explicit FlagSet(const ENUM&)` | `FlagSet::from_flag`, `From<E>` |
| copy constructor, `operator=(const FlagSet&)`, destructor | `Clone`/`Copy`, drop |
| `operator=(const ENUM&) = delete` | no counterpart: Rust never converts on assignment |
| `operator==` | `PartialEq`/`Eq` |
| `operator&`, `operator&=` (enum, set) | `BitAnd`, `BitAndAssign` for `E` and `FlagSet<E>` |
| `operator\|`, `operator\|=` (enum, set) | `BitOr`, `BitOrAssign` for both |
| `operator+`, `operator+=` (enum, set) | `Add`, `AddAssign` for both (same as `\|`) |
| `operator-`, `operator-=` (enum, set) | `Sub`, `SubAssign` for both |
| `isSuperSetOf(const FlagSet&)` | `is_super_set_of` |
| `isSuperSetOf(const ENUM&)` | `is_super_set_of_flag` |
| `empty` | `is_empty` |
| `value` | `value` |
| private `getPow_` | private `bit_of`, via `FlagEnum::bit` |

### `QCBase`

| Source | Rust |
|---|---|
| `enum class Requires` (`NOTHING` .. `ID`, `SIZE_OF_REQUIRES`) | `Requires` (bit positions 0..6), `Requires::SIZE_OF_REQUIRES`, `Requires::ALL` |
| `names_of_requires[]` | `qc_base::NAMES_OF_REQUIRES`, `Requires::name` |
| `enum class ToleranceUnit` (`AUTO`, `PPM`, `DA`, `SIZE_OF_TOLERANCEUNIT`) | `ToleranceUnit`, `ToleranceUnit::SIZE_OF_TOLERANCE_UNIT` |
| `names_of_toleranceUnit[]` | `qc_base::NAMES_OF_TOLERANCE_UNIT`, `ToleranceUnit::name` |
| `SpectraMap()`, `SpectraMap(const MSExperiment&)`, `~SpectraMap` | `SpectraMap::new`, `SpectraMap::from_experiment`, drop |
| `SpectraMap::calculateMap`, `at`, `clear`, `empty`, `size` | `calculate_map`, `at`, `clear`, `is_empty`, `len` |
| `using Status = FlagSet<Requires>` | `type Status = FlagSet<Requires>` |
| virtual `getName`, `requirements` | trait `QCBase::name`, `QCBase::requirements` |
| `isRunnable` | `QCBase::is_runnable`; its log lines `QCBase::runnable_warnings` |
| static `isLabeledExperiment` | `qc_base::is_labeled_experiment` (shares `kernel::map_operations::has_isobaric_analyzer`) |
| template static `hasPepID` | `qc_base::has_pep_id_feature_map`, `has_pep_id_consensus_map` |

### Metrics

| Source | Rust |
|---|---|
| `SpectrumCount()`, destructor, `getName`, `requirements` | `SpectrumCount::new`, drop, `QCBase` impl (`SpectrumCount::NAME`) |
| `SpectrumCount::compute` | `SpectrumCount::compute -> BTreeMap<usize, u32>` |
| `TIC()`, destructor, `getName`, `requirements` | `Tic::new`, drop, `QCBase` impl (`Tic::NAME`) |
| `TIC::Result` and its `operator==` | `TicResult` (`PartialEq`) |
| `TIC::compute(exp, bin_size = 0, ms_level = 1)` | `Tic::compute(exp, bin_size, ms_level) -> Result<TicResult>` |
| `TIC::addMetaDataMetricsToMzTab` | `Tic::add_meta_data_metrics_to_mz_tab` |
| `TIC::getResults` | not ported: declared in the header, defined nowhere (a call cannot link) |
| `FeatureSummary()`, destructor, `getName`, `requirements` | `FeatureSummary::new`, drop, `QCBase` impl |
| `FeatureSummary::Result` and its `operator==` | `FeatureSummaryResult` |
| `FeatureSummary::compute` | `FeatureSummary::compute -> Result<FeatureSummaryResult>` |
| `MissedCleavages()`, destructor, `getName`, `requirements` | `MissedCleavages::new`, drop, `QCBase` impl |
| `MissedCleavages::compute(FeatureMap&)` | `compute_feature_map(&mut FeatureMap)` |
| `MissedCleavages::compute(prot_ids, pep_ids)` | `compute(&[ProteinIdentification], &mut [PeptideIdentification])` |
| `MissedCleavages::getResults` | `results` |
| private `get_missed_cleavages_from_peptide_identification_` | private `count_peptide_identification` |
| `IdentificationSummary()`, destructor, `getName`, `requirements` | `IdentificationSummary::new`, drop, `QCBase` impl |
| `IdentificationSummary::UniqueID` | `UniqueId` (default count 0, threshold -1) |
| `IdentificationSummary::Result` and its `operator==` | `IdentificationSummaryResult` |
| `IdentificationSummary::compute` | `compute(&[ProteinIdentification], &mut [PeptideIdentification]) -> Result<_>` |

## Preserved source conventions

- `Requires::Nothing` is bit 0, a real bit: a metric requiring it would run
  only when the caller's status holds it too.
- `isRunnable`'s note is kept word for word, one line per missing input in
  bit order: `Note: Metric '<name>' cannot run because input data '<file>' is
  missing!`.
- `SpectraMap`: a repeated native ID maps to the last spectrum carrying it.
- `hasPepID` and `MissedCleavages::compute(FeatureMap&)` look at top-level
  features and unassigned identifications only, never at subordinates, as
  `MapUtilities::applyFunctionOnPeptideIDs` visits them (features first).
- `TIC::compute`: one point per spectrum in spectrum order at the spectrum's
  RT with its `f32`-accumulated intensity sum (`MSSpectrum::calculateTIC`);
  binning only for a positive bin size, through `LinearResamplerAlign`;
  intensities truncated to integers, RTs narrowed to `f32`; relative
  intensities `i / max * 100` in `f64` narrowed to `f32`, all 0 when the
  maximum is 0; a jump is `current > previous * 10`, a fall
  `current < previous / 10` with integer division.
- `addMetaDataMetricsToMzTab`: the vocabulary label is `total ion current`
  (the term's name, not a vocabulary label), the key is `custom.size()`, so
  sparse keys can be overwritten, and RTs are written with
  `StringUtils::toStr(float, false)`.
- `FeatureSummary`: the mean RT deviation is accumulated and divided in `f32`.
- `MissedCleavages`: fragments minus one of a digest with no missed cleavages;
  an empty sequence digests to nothing and so counts 4294967295 missed
  cleavages (`UInt32(size_t(0) - 1)`). The class documentation promises
  `FWHM` and `mass` annotations the code never writes; neither does the port.
- `IdentificationSummary`: a PSM is a peptide identification that is not
  `empty()` in the source's sense (identifier, hits, threshold 0, score type,
  higher-score-better), which is narrower than the crate's
  `PeptideIdentification::is_empty`; unique peptides are the unmodified
  sequences of first hits; the FDR thresholds are read from the first
  identification only, and only for the score type `FDR`; missed-cleavage
  sums wrap in `u32` as the source's `UInt`; a mean over nothing is NaN.

## Native differences

- `FlagEnum::bit` of 64 or more contributes no bit; the source's precondition
  is checked in debug builds only and undefined otherwise.
- `TIC::Result::area` is a `u64` summed exactly; the source's `UInt` wraps
  modulo 2^32 once the area exceeds 4,294,967,295, which real LC-MS runs
  reach. The jump test multiplies exactly; the source's `UInt` product wraps
  once a bin exceeds 429,496,729 and then miscounts. Both agree with the
  source whenever the source does not wrap.
- A TIC intensity beyond `u32::MAX` saturates (undefined in the source); NaN
  becomes 0.
- `TIC::compute` returns `Result` for the binning path (invalid spacing, a
  trace the resampler refuses); the source has no error path.
- `SpectrumCount` saturates its counts, `FeatureSummary` refuses more than
  `u32::MAX` features; the source's `UInt` would wrap.
- `FeatureSummary::compute` returns an error for a non-numeric `rt_deviation`
  where the source's `DataValue` conversion throws.
- `IdentificationSummary::compute` refuses empty `pep_ids`
  (`MissingInformation`), where the source reads `front()` of an empty list.
- Log output is returned, not written: `QCBase::runnable_warnings`,
  `MissedCleavages::log`.
- The missed-cleavage digest counts products (`ProteaseDigestion::
  peptide_count`) instead of materialising them.

## Checked boundaries and evidence

- Tier 3, `tests/qc.rs`: every section of `FlagSet_test.cpp`,
  `QCBase_test.cpp`, `SpectrumCount_test.cpp`, `TIC_test.cpp`,
  `FeatureSummary_test.cpp`, `MissedCleavages_test.cpp` and
  `IdentificationSummary_test.cpp` with their literals (the `FlagSet`
  precondition sections become the no-bit behaviour above).
- Tier 4, `tests/qc.rs`: the warning text of `is_runnable`, `hasPepID`,
  `addMetaDataMetricsToMzTab` (value text derived from
  `StringUtils::toStr(float, false)`), binning on an aligned grid, the
  non-numeric deviation, the source emptiness of a peptide identification.
- Tier 1, `tests/mzqc_file.rs`: `SpectrumCount`, `TIC` (all levels, MS1 with
  jump/fall, MS2), `FeatureSummary`, `MissedCleavages` and
  `IdentificationSummary` reproduce every metric value of the retained
  `TOPP_QCCalculator_2` output.
