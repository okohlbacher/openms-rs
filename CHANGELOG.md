# Changelog

## 0.1.0 — Ongoing native Rust SDK port

- Add complete native peak-file option state, metadata/Product hash traits, and source-compatible MassDecomposition count records with checked arithmetic.

- Update the SDK target to `54a232f`, retaining historical fixture pins and explicit changed-source compatibility reviews.
- Complete native FASTA parsing/file/stream/seek/progress lifecycle, including source modified sequences and PEFF prologue handling.
- Add bounded experiment aggregation and XIC extraction with all four source reducers and product m/z metadata.
- Preserve chromatogram product isolation/scalar metadata through mzML; add bounded indexed-mzML footer/offset parsing and index detection.
- Validate with Rust 1.98 and 1.85 and correct the newer compiler test-formatting check.

- Add owned logging with source routing, prefixes, duplicate caching and notifications, and replaceable progress backends with process CPU timing.
- Complete public modification collection for feature/consensus maps and nested subordinates; charge even empty search-name lookups.
- Add gzip/bzip2 INI loading with content detection, decoded limits and atomic parameter updates, preserving source plain output.

- Add featureXML and consensusXML adapters with typed feature metadata, processing history, assigned/unassigned identifications, checked custom chemistry, source filters, and protein-group quantity ownership guards.
- Add reusable gzip/bzip2 transport using Rust backends, atomic file output and feature/consensus FileHandler dispatch.
- Add portable modification definition records and owned-registry registration shared across all three identification XML formats.
- Add native filesystem helpers, explicit runtime/data/configuration paths and owned temporary resources, with documented platform conventions.
- Native API migration: feature/map/column metadata now uses typed MetaInfo; map records gain processing history and loaded-file path/type.

- Add native parameter values, hierarchical parameters, defaults/restrictions/update/copy/merge, forward traces and both command-line parsers, with source quirks and transactional errors documented.
- Add default-parameter lifecycle handling, pure typed-state callbacks, source warning policies and atomic leaf-key metadata export.
- Add OpenMS parameter XML/INI read/write, legacy fixtures, UTF-8/Latin-1/UTF-16 input, full-precision special floats and checked metadata preservation.
- Add complete source TextFile/CsvFile helpers and ListUtils/StringListUtils operations, using native streams, slices and standard containers.
- Document four reviewed standard-container/alias equivalents separately from scientific class coverage.
- Add an exhaustive completion ledger for 786 registered public SDK headers and direct dependencies from 146 TOPP sources, with reviewed coverage distinguished from names and source-reference evidence.
- Port all 73 file-type descriptors and lexical filename helpers; add native experiment dispatch, bounded format sniffing, gzip transport and atomic path output.
- Add source MS2/DTA2D readers, DTA2D filters/storage/TIC, and a checked native MS2 writer with original source fixtures.
- Add mzML referenceable parameter groups, forward header references, shared inline validation and bounded definition/reference expansion.
- Add graph referential cleanup with all five source switches and conditional predicate filtering; surviving IDs remain stable and removed IDs cannot be reused.
- Add parent/match grouping, legacy sequence/evidence converters and charged peptide fragment mass queries, with focused source-derived tests.

- Update the target to reduced Core SDK 4.0.0 at `6bfc0e4`; record retained/removed source scope and 220 unchanged historical reference paths, expose the exact target in Rust, and verify SDK/RNA resource consistency in CI.
- Add the native identification sequence/provenance graph: typed graph-owned references, input/software/search/score history, parent/peptide/oligo registration, translated merge/copy, and inclusive sequence coverage with checked atomic updates.
- Add both RNase graph digestion operations, preserving custom RNA chemistry, parent flanks and processing history with cumulative parsing/digestion limits. Add independent source references, cross-parent budget/rollback regressions and a provenance-aware RNA example.
- Add typed peptide fragment formula dispatch plus identification observations, compounds, adducts, molecule references and observation matches. Registration, source merge semantics, score and annotation history, best-match queries, graph-owned translation and atomic resource checks are covered by focused source-derived tests.

- Add all fourteen RNA enzymes and native digestion, fixed/variable modified-RNA enumeration, and both annotated RNA spectrum-generation operations, with source-derived processing references and digestion-to-mzML workflows.
- Add RNA nucleoside records, the complete pinned 378-record registry, bounded TSV and optional MODOMICS JSON providers, and owned nucleic-acid sequences with source terminal/linkage/slicing and fragment mass conventions.
- Add independent RNA source formulas, masses, slices and complete registry projections, negative-ion isotope/mzML workflows, data provenance and a regeneration tool. MODOMICS redistribution terms remain separately tracked for release.
- Add the complete native Tagger API for measured spectra and m/z arrays: source mass/charge/tolerance traversal, fixed/variable modifications, I/L alternatives, deterministic two-stage registry lookup and bounded atomic append.
- Add all six original tag-count and 120 membership assertions, independent edge/mass-table references, digestion/indexing/mzML workflows and an extraction example.
- Add native AdductInfo parsing, electron-aware mass conversion and mono/average shifts, signed formula compatibility, source integer/whitespace grammar and independent ion-composition/mzML workflows.

- Add native DecoyGenerator with all source reversal/shuffle operations, seeded MT19937-64 and Boost interval mapping, preserved cache/cleavage quirks, transactional state and shared work/output limits.
- Add all thirteen literal decoy reference cases, independently reproducible integer RNG checks, FASTA/indexing/FDR/idXML workflows and a FASTA decoy example; preserve the helper's Boost license notices.

- Added all three SpectrumAnnotator operations and four IonNaming helpers: matched spectrum arrays, hit peak annotations, matching statistics and bounded charge display/parsing. Preserved ppm duplicate conventions, source no-op flags and f32 errors, with atomic updates and checked undefined statistics.
- Shared actual isotope/alignment work and precharged fragment/loss work across candidate peptides; coarse convolution also shares its allowance across theoretical envelopes. Anonymous mass-tag spellings now use shared immutable ownership during peptide slicing without changing their text or equality.
- Added original annotation fixtures, independent error/ratio references, modified-digestion and XML workflows, cumulative-work regressions and an annotation example.

- Added charge and isoelectric point with four pKa scales; seven hydrophobicity scales, GRAVY, moving profiles and moments; ten AAindex scales and source residue indicators; and gas-phase basicity. Preserved parent-residue and terminal-annotation conventions with checked inputs and work limits.
- Gas basicity preserves ordinary source arithmetic, uses a stable energy-domain retry for overflow and corrects empty-sequence high-temperature cancellation. Added complete table fixtures, independent numerical references, digestion/idXML workflows and a peptide-property example.

- Added an owning fallible fine-isotope iterator, absolute/relative threshold selection, raw/log probability access, checked peak conversion and custom binary64 isotope populations for both streaming and materialized patterns. Iterator formula input follows the raw wrapper's charge-ignoring convention; the high-level generator retains its natural-H convention.
- Shared the existing bounded isotope search with incremental iteration, preserving deferred expansion, source rounding and cumulative theoretical-spectrum work limits. Added an enrichment/streaming example.

- Added native fine isotope configurations with absolute/relative thresholds and total-probability coverage, preserving source abundance/output rounding, fixed labels and natural-H charge handling without a C++ dependency.
- Added fine theoretical fragment, neutral-loss and precursor envelopes with shared work limits and atomic append. Isotope loss intensities now retain the full floating-point product until final f32 storage, including the existing coarse path.
- Native API: `TheoreticalIsotopeModel::Fine { unexplained_probability }` selects fine envelopes; this enum now implements `PartialEq` without `Eq` because the probability is floating point. IsoSpec layer order remains outside the implemented surface.

- Added scalar precursor purity, SPS fragment matching, fuzzy scan purity, RT interpolation and experiment scoring with source numerical conventions and shared work limits.
- Connected precursor acquisition fields directly to spectrum/chromatogram records, including isolation offsets, activation, mobility and spectrum references. Parent lookup honors earlier referenced scans before acquisition-order fallback.
- Expanded the mzML subset to preserve supported precursor acquisition fields and reject unrepresentable native fields before writing. DTA/MGF exports reject richer precursor metadata rather than discarding it.
- Native API migration: `Precursor` is now `Clone` rather than `Copy`; `Precursor::new` is no longer const. Complete struct literals need `..Precursor::default()`. `PrecursorInfo` owns only `peak: Precursor`; field access delegates to that record through `Deref`/`DerefMut`, with explicit conversions in both directions.

- Added internal b/a fragments, abundant immonium peaks, activation-method presets and the compact mass-only ladder helper, preserving source charge, terminal, loss, ordering and floating-point conventions with bounded atomic output.

- Added shared ownership for modification records, caller-owned registry parsing/setters/generation, and exact custom-registry idXML interchange without leaked storage.
- Added bounded native OBO loading, PSI-MOD alias semantics, pinned XLMOD monolinks and a separate crosslink lookup database with preserved specificity and mass conventions.
- Added anonymous modification definitions and inference, retaining spelling, source compatibility behavior, and full-residue/terminal absolute-mass anchors.
- Added absolute-formula replacement for changed residues, including restoration of unknown B/Z/X chemistry; unchanged-residue records keep their original formula and mass.
- Preserved distinct same-name custom chemistry in protein modification observations, sequence duplicate filtering, peptide identity keys and conflict resolution.
- Native API migration: registry UniMod IDs are optional, registry entries and known sequence annotations use `Arc<ResidueModification>`, and `to_unimod_string()` returns `Result<String>` with non-UniMod mass fallback. `to_accession_string()` retains vocabulary accessions; peptide identity keys now own complete sequence values.

- Added fixed and variable peptide-modification generation with bounded combinatorial output, source terminal/ordering conventions, preserved existing annotations and atomic append.
- Added modification definition sets, compatibility checks, inference from peptide identifications and delta/absolute mass matching, including explicit stored absolute masses on owned modification records.
- Added a digest-variant enumeration example and a workflow verifying chemical formulas, independent fragment shifts, inferred search definitions and idXML round trips. The identification example now uses the fixed-modification generator.
- Preserved custom formula-free modification masses through peptide generation and monoisotopic fragments, with explicit unavailable composition and source absolute-mass precedence.
- idXML writing now rejects typed modification placements that cannot round-trip through its sequence syntax before touching the destination.

- Added native iterative peak picking with HiRes seeds, original-index regions, source refinement/suppression conventions and aligned integrated-intensity/width annotations.
- Added sliding and jumping window filters with source window/duplicate rules, stable ties and atomic spectrum/experiment updates.
- Added iterative mean noise estimation with three clipping passes, source histogram arithmetic, sparse-window diagnostics and checked historical percentile behavior.
- Added an end-to-end profile noise → iterative centroiding → window filtering example and exact mzML annotation round-trip tests. Independent reviews cover fourteen refinement cases, 480 window-selection cases and 192 mean-noise configurations.

- Added native EMG fitting with source training selection, analytic gradients, iRprop+ optimization, best-fit diagnostics and bounded truncated-side reconstruction.
- Added optional EMG preprocessing for peak integration, baseline estimation and shape metrics, with typed parameters and explicit finite bounds.
- Added a cropped-peak fitting example, independent area/centroid identities and a fitted-trace mzML workflow. The scientific core now uses the pure Rust `libm` dependency for the complementary error function.

- Added native legacy/corrected chromatogram peak picking with independent seed/boundary noise settings, source overlap handling, exact original indices and aligned peak annotations.
- Added sampled peak integration, all source integration/baseline choices and shape metrics for spectra and chromatograms, preserving source floating-point and Simpson conventions.
- Added a chromatogram integration example and a workflow covering exact boundaries, raw intensity sums, time-weighted areas and mzML annotation interchange.
- Added named mzML float/integer/ASCII string arrays, preserved empty strings and annotation placeholders, and bounded cumulative decoded array storage with validation before output.

- Added source Poisson/KL deisotoping with threshold/top-N preprocessing, longest-cluster selection, original-index membership, source shared-isotope behavior, optional disjoint clusters and checked aligned annotations.
- Corrected the threshold filter default to source 0.05; preserved direct Poisson recurrence rounding with a stable overflow fallback.
- Added explicit isotope sharing to the simple deisotoper and retained exact source precursor arithmetic in both methods.

- Expanded `AASequence` to unresolved B/Z/X and owned numeric mass tags with source precision-dependent registry lookup; formula and mass APIs now return `Result`.
- Preserved numeric annotations through digestion, indexing, idXML, modification mapping and conflict/filter operations; monoisotopic fragments support known mass without invented composition.

- Added peptide-to-protein indexing with ambiguity/mismatch rules, enzyme-aware evidence, I/L handling, decoy inference and atomic run reconstruction.
- Added basic protein score aggregation, representative counts, indistinguishable groups and greedy resolution for vector, single-run and consensus-map inputs.
- Added feature/spectrum identification conflict resolution, file-origin partitions and a FASTA-to-inference/FDR/idXML workflow test.
- Added score categories and atomic score switching, HyperScore/Morpheus fragment scoring, and native peptide/protein identification filters.
- Added legacy/Basic target-decoy FDR, peptide/protein q-values, picked proteins/groups, posterior estimates and ROC with pinned-source fixtures.
- Added optional bounded idXML 1.5 read/write with typed metadata, run/evidence references, independent schema validation and a complete identification-processing round trip.
- Added typed metadata, CV terms and validated acquisition/settings records.
- Added peptide/protein identification records, evidence coverage, observed protein modifications and spectrum/feature attachments; peak writers reject unsupported identification loss.
- Added retention-time regression, interpolation and robust LOWESS models and atomic experiment/feature/consensus transformations with source-derived reference tests.
- Expanded native digestion to all 33 pinned enzymes with full/semi/nonspecific enumeration and checked validity/count operations.
- Added a synthetic peptide-identification workflow example with checked fragment matches, evidence and protein coverage.
- Added native spectra, chromatograms, experiments and basic precursor metadata.
- Added feature/consensus containers, checked maps and IDs, scan-envelope geometry, containment and consensus/decharge summaries.
- Added formula/element chemistry, modified peptide masses, b/y fragments and tryptic digestion.
- Added an embedded UniMod/custom modification registry with preserved data/source licenses.
- Added coarse isotope patterns, convolution, enrichment, averagine and conditional fragment estimates.
- Added configurable theoretical peptide spectra with terminal ion series, neutral losses, precursor peaks, coarse isotope envelopes and aligned annotations.
- Added spectrum alignment, sparse binning, and common spectrum similarity scores.
- Added Gaussian/Savitzky–Golay smoothing and morphological baseline correction.
- Added HiRes spectrum/chromatogram centroiding, natural cubic splines, histogram-median noise estimation, FWHM and peak boundaries.
- Added simple C13-spacing deisotoping with optional charge conversion, intensity summation and annotations.
- Added normalization, filtering, scaling and linear resampling.
- Added DTA, FASTA and MGF interchange plus an optional documented mzML subset.
- Added source inventory, coverage notes, fixtures, examples, tests and CI configuration.
- Added a modified-peptide mass, isotope and theoretical-spectrum example.
- Recorded deliberate differences and remaining C++ library scope.

This version is an initial API and may change. It has not been published.
