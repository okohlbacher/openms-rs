# FeatureFindingMetabo

The native [implementation](../src/analysis/feature_finding_metabo.rs) covers the
complete class-specific `FeatureFindingMetabo` operation and all nineteen
configuration fields from reduced SDK
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. It consumes native `MassTrace` records
and produces `FeatureMap` plus optional grouped chromatograms. The preceding
`MassTraceDetection` and `ElutionPeakDetection` operations remain separate.
There is no general SVM training/loading facade; the two fixed metabolite
predictors are documented in [METABO_PREDICTOR_SUPPORT.md](METABO_PREDICTOR_SUPPORT.md).

## API and configuration

`FeatureFindingMetabo::new()` and `with_options(options)` return checked results.
`options()` borrows the requested configuration; `set_options` validates and
prepares its element alphabet before replacing any state. `limits` and `logger`
are public. The default logger has source CMD behavior and can be changed to
NONE or a caller backend through the existing `ProgressLogger` API.

`run(&mut traces, &mut rng)` returns `FeatureFindingMetaboOutput` containing the
feature map, chromatogram groups, and explicit configuration diagnostics.
`run_into` replaces both existing outputs and returns their previous ownership
plus diagnostics. Previous output records are not inspected or destroyed inside
this operation. The caller owns the seeded `UniqueIdGenerator`; the source global
singleton is not used.

| Option | Default | Behavior |
| --- | --- | --- |
| `local_rt_range` | 10 | Candidate centroid RT window, including when RT scoring is disabled |
| `local_im_range` | 0.02 | Candidate centroid IM window when any input trace has IM |
| `local_mz_range` | 6.5 | Candidate m/z window and isotope offset count |
| `charge_lower_bound` / `charge_upper_bound` | 1 / 3 | Inclusive charge traversal |
| `chrom_fwhm` | 5 | Stored source parameter, unused by this algorithm |
| `report_summed_intensities` | false | Sum member trace intensities instead of the first trace |
| `enable_rt_filtering` | true | RT overlap and cosine scoring |
| `isotope_rt_overlap_reference` | Longer | Longer-width or shorter-width/apex overlap rule |
| `min_isotope_rt_overlap` | 0.7 | Validated inclusive range 0–1 |
| `isotope_filtering_model` | Metabolites5 | Metabolites2, Metabolites5, Peptides, or None |
| `mz_scoring_13c` | false | C13 mean spacing instead of the Kenar mean |
| `use_smoothed_intensities` | true | Intensity selection for scoring/model inputs |
| `report_smoothed_intensities` | true | Reporting mode; disabled if smoothed scoring is disabled |
| `report_convex_hulls` | false | Member trace hulls in each feature |
| `report_chromatograms` | false | Raw chromatograms grouped by accepted feature |
| `remove_single_traces` | false | Discard charge-zero singleton hypotheses |
| `mz_scoring_by_elements` | false | Element isotope mass-defect windows; disables isotope filtering |
| `elements` | CHNOPS | Empirical-formula alphabet; counts do not multiply element entries |

Diagnostics replace the two configuration warnings. Element scoring changes the
stored *effective* model to None on successful runs, including empty runs;
`effective_isotope_filtering_model()` exposes that state. A subsequent
`set_options` resets it to the requested model. `contains_im_data()` updates from
each successful nonempty input; an empty call preserves its previous value, as
in the source.

## Source behavior retained

Input is sorted by cached centroid m/z. Hypotheses start from each trace, charge,
and accepted isotope prefix; they are not every possible combination. Candidate
m/z, RT and IM windows are inclusive, while Gaussian cutoffs and element-window
interiors use their source strict boundaries. RT scoring merges raw peaks by RT,
retaining only groups with exactly two entries. Duplicate entries from a single
trace can therefore form a pair; groups of three or more are omitted. Shorter
reference overlap also tests the longer trace's configured raw/smoothed apex.

The source computes RT, m/z and peptide intensity scores before combining them,
even if an earlier score is zero. Peptide scoring takes **raw** intensities for
already selected members and the configured intensity for the new candidate.
Mean scoring retains separate logarithm/exponential variance arithmetic and
separate square roots. Element offsets retain integer division of isotope gaps,
and include the lightest listed isotope independently of abundance.

Hypotheses are accepted by descending score. Exclusion uses trace **labels**, so
separate traces with the same label collide. Singletons carry charge zero and
legal-pattern metadata -1. Model inputs cap neutral mass at 1000 and fill absent
isotope ratios with zero; model label 2 is the accepted pattern. Optional
chromatograms use raw RT/intensity and first-member precursor m/z, have source
BasePeak type/native IDs, and remain in acceptance order. The feature map is
subsequently sorted by m/z, so its order can differ from chromatogram groups.
Features preserve labels, counts, apex, member intensity/RT/mz/IM vectors,
isotope distances, legality and FWHM metadata. Feature/map IDs consume source
MT19937-64 draws in acceptance order, then one map draw for nonempty input even
when no hypothesis is accepted. Empty input consumes no RNG draw.

## Checked native boundaries

All fallible processing is staged. A failure preserves input order, previous
outputs, RNG state, effective model and IM flag. On success, only then is input
reordered and state published. Source OpenMP/unstable-sort score ties have no
portable ordering guarantee; native equal scores preserve serial sorted-seed,
charge and prefix order, and equal m/z/RT values preserve encounter order.
Progress backend side effects cannot be rolled back; backend failures still
prevent scientific state publication.

Consumed nonfinite coordinates, intensities or arithmetic, zero total intensity,
nonfinite normalized peptide patterns, overflowing f32 outputs, and invalid
consumed source integer conversions return errors. Unused fields are not
speculatively validated. Negative finite RT/IM windows and inverted charge ranges
retain their defined source no-match branches. Negative m/z range with a positive
charge reaches a source out-of-range floating-to-unsigned offset conversion and
is rejected; a zero charge can keep its defined zero-offset behavior. The isotope
offset is bounded to the source signed32 loop domain. Empty/short trace ranges
and unavailable smoothed values fail only when the selected branch needs them.

Peptide averagine reuses the existing native empirical-formula and coarse-isotope
operations. Their probabilities are accumulated in f64; the source stores f32
isotope peaks between convolutions. This inherited precision difference is
explicit, and peptide scores are not claimed bitwise equal to C++. Default
metabolite modes do not use averagine. No new numerical backend is introduced.

## Resources and validation

Defaults cap one call at 1,000,000 traces, 10,000,000 aggregate peaks, 1,000,000
hypotheses/features, 50,000,000 work units and 256 MiB of cumulative allocation
allowance. These are configurable via `FeatureFindingMetaboLimits`. Traversals,
comparisons, geometric vector growth, isotope scoring, labels, output metadata,
hulls, chromatograms, temporary permutations and RNG copies share the same
ledger. The fixed predictors charge their full support-vector evaluation against
it and allocate no heap memory. The existing bounded MassTrace/FeatureHypothesis
queries are conservatively precharged before calls; ordinary Rust Clone/drop and
caller progress backends are outside the operation's resource guarantee.

The [tests](../tests/feature_finding_metabo.rs) execute the complete default chain:
**83** C13 features, **81** Kenar features and **80** element-window features.
All 81 committed featureXML records' scientific fields, typed metadata and hull
points are compared. The source writer compresses a hull copy before emitting
XML; the test applies the same compression. The historical comparator permits
absolute error <=1, or same-sign nonzero reciprocal-symmetric ratio <=1.001;
generated IDs and XML presentation are not numerical oracles. The source IM
projection also exercises the complete chain and retained IM metadata.

Independent small cases cover raw/smoothed/summed reporting, RT duplicates,
strict boundaries, source integer-gap behavior, both metabolite models and the
mass cap, charge/label/tie branches, IM-state reuse, output ownership, RNG and
progress failures, and aggregate limits. The predictor's six private tests use
488 executed LIBSVM rows and bit-check every packaged model constant. They are
counted once through this module, not again through its isolated test harness.

[Source/fixture provenance](../tests/data/feature_finding_metabo_provenance.json)
records exact SDK inputs and the independent fixture generator. Expected values
are projected from committed source data or derived from source expressions;
they are not recorded from this Rust implementation.

The independent peptide reference additionally covers four 1000 Da scores with
source f32 convolution rounding and all 24 active-element orders (the source
formula map is pointer ordered). Its maximum order-related score spread is
below 3.8e-8; the tests allow absolute score error 1e-7 for the inherited f64
implementation. A public two-trace case verifies the mixed raw/configured
peptide candidate path against those independent values.

Focused validation uses Rust 1.98/all-features and Rust 1.85/no-default,
`--locked --offline -j2`, including this module's public/private tests and
adjacent MassTrace, MassTraceDetection, ElutionPeakDetection and isotope tests.
Strict scoped library/test Clippy and formatting are checked on both versions.
To reproduce the independent fixtures, run the two generator scripts with the
pinned SDK directory; the first also takes an output directory. Both read source
inputs only. `python3 tools/generate_metabo_isotope_models.py --check` verifies the
immutable model tables against the packaged original resources.

The [file example](../examples/find_metabolite_features.rs) loads the original
`FeatureFindingMetabo_input1.mzML` unchanged, runs the three default operations
and writes 81 features. Its output passes the original FeatureXML 1.9 schema;
all scientific scalars and typed metadata match the committed expected file,
including 1,053 numerical comparisons under the source comparator. Hull export
is off by default in this example and covered separately above. The mzML reader
accepts this file's ASCII-only Latin-1 declaration without transcoding; full
header transport and complete TOPP command behavior remain outside this example.
