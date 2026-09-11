# Porting status: ongoing native Rust port

Metabolite feature finding, complete experiment settings, DateTime and ProForma
annotation parsing, writing, resolution, [mass calculation](PROFORMA_MASS_SUPPORT.md),
[AASequence conversion](PROFORMA_CONVERSION_SUPPORT.md) and
[spectrum generation](PROFORMA_SPECTRA_SUPPORT.md) are native. The [completion ledger](CORE_SDK_COMPLETION.md) records reviewed
operation groups and remaining work. Full SDK/TOPP certification remains open.

The target is a feature-complete reduced Core SDK suitable for porting TOPP tools, with an idiomatic Rust API. Spectra, chemistry and common processing were the starting priorities. This document describes the implemented surface rather than claiming parity for every method of a similarly named C++ class. The [completion ledger](CORE_SDK_COMPLETION.md) tracks all registered public headers and direct TOPP dependencies.

The current target is SDK 4.0.0 at `82ce5b3`; the [SDK update](CORE_SDK_UPDATE.md) records the exact source inventory and extracted product backends excluded from this port’s remainder. Historical scientific fixtures retain their original pins.

[CV mapping records and XML loading](CV_MAPPING_SUPPORT.md) cover all five class-specific source APIs, with atomic loads and explicit compatibility corrections. [General semantic validation](SEMANTIC_VALIDATOR_SUPPORT.md) now supplies the complete class-specific mapping/term validator with ordered diagnostics and bounded, reusable operations. Format-specific derived validators and XSD validation remain separate.

## Mass traces and additional SDK values

The complete [MassTrace surface](MASS_TRACE_SUPPORT.md) is available, including
all centroid updates, cached widths and borders, quantification, smoothing and
hulls. [SDK constants](CONSTANTS_SUPPORT.md) expose every numeric entry and metadata
key with exact source values. The [monosaccharide database](MONOSACCHARIDE_SUPPORT.md)
provides the complete built-in source lookup surface. Mass-trace detection is
now available. ProForma sequence parsing and mass calculation are implemented;
[AASequence conversion](PROFORMA_CONVERSION_SUPPORT.md) and [ordinary/XLMS spectrum generation](PROFORMA_SPECTRA_SUPPORT.md) are implemented.

## Mobility containers, array descriptions and path operations

[MobilityPeak1D/Mobilogram](MOBILOGRAM_SUPPORT.md) now represent the source
value and scientific container operations, with current ranges and explicit
source partial-swap/equality behavior. Inherited generic range algebra and
mobility-bearing experiment operations remain. Generic arrays retain owned
metadata and shared processing handles; the existing constructor is unchanged,
and struct literals require defaults for the two new fields. Native processing
preserves retained descriptions. [XML guards](DATA_ARRAY_XML_SUPPORT.md) prevent
loss at transports that cannot encode them.

[IMSWeights](IMS_WEIGHTS_SUPPORT.md) implements the complete standalone weight
utility, including source GCD and floating-rounding behavior. The
[mzML filesystem API](MZML_PATH_SUPPORT.md) applies scientific load settings,
magic-based compression detection and atomic replacement/publication.
Standalone IMS isotope, alphabet and decomposition APIs are also implemented;
[source-supported mzML headers](MZML_HEADER_SUPPORT.md) are transported; consumers remain outstanding. Loaded-file
bookkeeping is retained in the experiment settings.

[ControlledVocabulary](CONTROLLED_VOCABULARY_SUPPORT.md) supplies complete
term definitions, OBO loading, graph queries, XML values and all five original
providers. General CV mapping and source-supported header transport are implemented;
[general semantic validation](SEMANTIC_VALIDATOR_SUPPORT.md) and [streaming consumers/transforms](MZML_CONSUMER_SUPPORT.md) are implemented. Centroid inspection and remaining mzML format-specific validation/writer options remain separate groups.

## Capability mapping

| C++ source area | Rust surface | Status and boundaries |
| --- | --- | --- |
| DATASTRUCTURES/Param, ParamValue | `param` | Full native values, ordered tree/traces, descriptions/tags/restrictions, defaults, update/copy/merge and command-line parsing; [details](PARAM_SUPPORT.md) |
| DATASTRUCTURES/DefaultParamHandler | `param::DefaultParamHandler` | Defaults, subsections, validation, typed-state update callbacks and leaf-key metadata export; [details](DEFAULT_PARAM_HANDLER_SUPPORT.md) |
| DATASTRUCTURES/ListUtils, StringListUtils | `data_structures` | Defined source conversions/formatting, membership, tolerance, prefix/suffix searches and ASCII case operations; [details](LIST_UTILS_SUPPORT.md) |
| FORMAT/TextFile, CsvFile, ParamXMLFile | `format::{TextFile,CsvFile,paramxml}` | Source literal text/CSV helpers and optional INI transport through 1.8.0; [text/CSV](TEXT_CSV_SUPPORT.md), [INI](PARAMXML_SUPPORT.md) |
| KERNEL/Peak1D, ChromatogramPeak | `Peak1D`, `ChromatogramPeak` | Coordinates, intensity, construction, equality |
| METADATA/Precursor | `Precursor` | Selected m/z, charge/intensity, activation, isolation, mobility, possible charges, CV terms and parent spectrum reference directly on kernel records; `PrecursorInfo` is a compatibility wrapper; [details](METADATA_SUPPORT.md) |
| KERNEL/MSSpectrum | `MSSpectrum` | Owned peaks and basic metadata; named float/integer/string arrays; validation, stable sorting, selection/retention, TIC, base peak, nearest and bounded searches, current ranges |
| KERNEL/MSChromatogram | `MSChromatogram` | RT peaks, precursor and product ions, basic metadata and aligned arrays; sorting, selection, ranges and nearest search |
| KERNEL/MSExperiment | `MSExperiment` | Spectrum/chromatogram storage, validation, RT sorting/search/filtering, MS-level queries, TIC chromatogram, ranges, referenced/acquisition-order parent lookup and complete window aggregation/XIC operations; [aggregation](EXPERIMENT_AGGREGATION_SUPPORT.md) |
| KERNEL/BaseFeature, Feature, FeatureHandle | `BaseFeature`, `Feature`, `FeatureHandle` | Measured values, metadata, attached peptide IDs, subordinate features, mass-trace hulls and owned handle identities |
| KERNEL/ConsensusFeature, FeatureMap, ConsensusMap | `ConsensusFeature`, `FeatureMap`, `ConsensusMap` | Mean/monoisotopic/decharge summaries, checked handle union and IDs, stable sorting, selection, ranges and column consistency; [details](FEATURE_SUPPORT.md) |
| DATASTRUCTURES/DataValue, METADATA/MetaInfo/CVTerm/settings | `metadata` | Typed values/units, CV terms, acquisition/precursor/instrument/processing/settings models and validated merging; [details](METADATA_SUPPORT.md) |
| METADATA/PeptideEvidence, PeptideHit/Identification, ProteinHit/Identification | `identification` | Owned records, stable score sorting, evidence/coverage, protein modification mapping and kernel attachments; [details](IDENTIFICATION_SUPPORT.md) |
| METADATA/ID/IdentificationData sequence and provenance layer | `identification::graph` | Stable owned IDs, input/software/search provenance, score histories, peptide/RNA parent matches, coverage, observations, compounds, adducts, molecule dispatch, observation matches, match groups, parent groups, bounded legacy sequence/evidence conversion, referential cleanup/filtering, atomic registration/merge/copy and RNase integration; graph persistence and full conversion remain; [details](IDENTIFICATION_GRAPH_SUPPORT.md) |
| ANALYSIS/ID/Scores, IDScoreSwitcherAlgorithm | `analysis::scores` | Six categories and 29 source names, numeric score switching/backups, per-record category discovery and atomic map adapters; [details](SCORING_SUPPORT.md) |
| ANALYSIS/ID/PrecursorPurity | `analysis::precursor_purity` | Scalar and batch scores, SPS matching, fuzzy scan purity and RT interpolation with bounded work and checked undefined source cases; [details](PRECURSOR_PURITY_SUPPORT.md) |
| ANALYSIS/ID/HyperScore, MorpheusScore | `analysis::psm_scoring` | Fragment matching, source score/error conventions, charge-aware and intensity-support variants; [details](SCORING_SUPPORT.md) |
| PROCESSING/ID/IDFilter | `analysis::id_filter` | Score/rank/sequence/charge/mass/evidence filters, duplicates and run-aware protein/evidence/group cleanup; [details](ID_FILTER_SUPPORT.md) |
| ANALYSIS/ID/FalseDiscoveryRate, IDScoreGetterSetter | `analysis::false_discovery_rate` | Legacy and Basic target/decoy curves, PSM/peptide/protein application, picked groups, posterior estimates and ROC; [details](FDR_SUPPORT.md) |
| ANALYSIS/ID/PeptideIndexing, selected AhoCorasickAmbiguous/DecoyHelper behavior | `analysis::peptide_indexing` | Bounded matching with protein ambiguities/mismatches, enzyme/I-L rules, decoy inference, evidence and run-wise protein reconstruction; optimized trie remains performance work; [details](PEPTIDE_INDEXING_SUPPORT.md) |
| ANALYSIS/ID/BasicProteinInferenceAlgorithm, selected IDBoostGraph operations | `analysis::protein_inference` | Representative score aggregation, counts, default grouping and greedy evidence resolution; vector/single-run/consensus conventions; [details](PROTEIN_INFERENCE_SUPPORT.md) |
| ANALYSIS/ID/IDConflictResolverAlgorithm | `analysis::id_conflict_resolver` | Best-score, matching-sequence and rank aggregation, between-feature intensity selection and per-spectrum deduplication; [details](ID_CONFLICT_SUPPORT.md) |
| ANALYSIS/ID/IDRipper | `analysis::id_ripper` | All three source origin annotations, checked in-memory file/run partitions and preserved skipped records; [details](ID_RIPPER_SUPPORT.md) |
| ANALYSIS/MAPMATCHING/MapAlignmentTransformer | `analysis::alignment_transformer` | Atomic experiment/feature/consensus/peptide RT updates, original-RT preservation, subordinate/hull/handle handling; [details](TRANSFORMATIONS_SUPPORT.md) |
| ANALYSIS/MAPMATCHING/TransformationModel family | `analysis::transformations` | Linear weighted-coordinate regression, linear/natural-cubic interpolation, robust LOWESS, inversion/deviations/windows; [details](TRANSFORMATIONS_SUPPORT.md) |
| DATASTRUCTURES/ConvexHull2D, selected 2D bounds | `ConvexHull2D`, `Point2D`, `BoundingBox2D` | Source scan-envelope containment, outline preservation, compression and rectangular bounds; not a general convex-polygon implementation |
| CHEMISTRY/ElementDB, EmpiricalFormula | `element_table`, `EmpiricalFormula` | Embedded 84-element data, labeled isotopes, mono/average masses, protonation charge, checked algebra; coarse and fine isotope generation are exposed separately |
| CHEMISTRY/AASequence, ResidueDB | `AASequence` | 20 standard residues plus U/O/J/B/Z/X, named/numeric residue and terminal annotations, checked formula/masses, subsequences and m/z; [details](SEQUENCE_SUPPORT.md) |
| CHEMISTRY/AASequence fragment convenience | `AASequence::fragment_ions` | b/y series only, charges 1..=maximum; modified residue/terminal masses included; no intensities, losses or isotope peaks in this convenience API |
| CHEMISTRY/TheoreticalSpectrumGenerator | `TheoreticalSpectrumGenerator` | a/b/c/x/y/z and z variants, internal fragments, immonium peaks, activation presets, compact mass helper, losses/precursors/coarse and fine envelopes and aligned names/charges; checked generation/append; [details](THEORETICAL_SPECTRA.md) |
| CHEMISTRY/TheoreticalSpectrumGeneratorXLMS | `TheoreticalSpectrumGeneratorXLMS`, `ProteinProteinCrossLink` | Complete class-specific linear/single/pair append operations, all source options and atomic aligned annotations; [source behavior and limits](THEORETICAL_XLMS_SUPPORT.md). ProForma wrappers are implemented separately; other OPXL records remain. |
| CHEMISTRY/ModificationsDB | `ModificationsDB`, `ResidueModification` | 3,035 UniMod/custom and 92 XLMOD monolink specificity records, shared owned handles, bounded OBO loading, aliases, caller records and mass lookup; [details](MODIFICATION_SUPPORT.md) |
| CHEMISTRY/CrossLinksDB | `CrossLinksDB` | Separate 56-record pinned XLMOD crosslink view, caller-owned bounded loading and shared registry searches; [details](CROSS_LINKS_SUPPORT.md) |
| CHEMISTRY/ModifiedPeptideGenerator | `ModifiedPeptideGenerator` | Fixed placement and bounded variable combinations, source reverse-site and terminal-path conventions, deterministic alternative order and atomic append; [details](MODIFIED_PEPTIDES_SUPPORT.md) |
| CHEMISTRY/ModificationDefinition, ModificationDefinitionsSet | `ModificationDefinition`, `ModificationDefinitionsSet` | Owned named/anonymous definitions, fixed/variable partitions, compatibility, delta/absolute mass matching and inference from peptide IDs; [details](MODIFICATION_DEFINITIONS_SUPPORT.md) |
| CHEMISTRY/IsoelectricPoint, HydrophobicityProfile, AAIndex, Residue hydrophobicity | `IsoelectricPoint`, `ProteomicsPkaScale`, `HydrophobicityProfile`, `HydrophobicityScale`, `AAIndex`, `AAIndexScale` | Four pKa scales, charge/pI, seven hydrophobicity scales, GRAVY/profiles/moments, ten indices and gas basicity with stable extreme-temperature evaluation; [details](PEPTIDE_PROPERTIES_SUPPORT.md) |
| CHEMISTRY/Ribonucleotide, RibonucleotideDB and data providers, NASequence | `Ribonucleotide`, `RibonucleotideDB`, `NASequence` | Full record identity and stored mass fields, all 378 pinned entries, bounded TSV/optional JSON providers, source parsing/slicing and all finite fragment formulas/masses; [details](RNA_SUPPORT.md) |
| CHEMISTRY/RNaseDB, DigestionEnzymeRNA, RNaseDigestion | `RNaseDB`, `DigestionEnzymeRNA`, `RNaseDigestion` | All fourteen enzymes, owned overrides, modified-code patterns, sequence and graph digestion, length/missed-cleavage limits, coordinates and end gains; arbitrary regex/XML remain; [details](RNASE_SUPPORT.md) |
| CHEMISTRY/ModifiedNASequenceGenerator | `ModifiedNASequenceGenerator` | Both operations, source maximum-one specificity behavior, exact subset/alternative order, bounded atomic output; [details](RNA_MODIFICATION_SUPPORT.md) |
| CHEMISTRY/NucleicAcidSpectrumGenerator | `NucleicAcidSpectrumGenerator` | Both operations, all nine series, source declared-mass and charge/precursor/sulfur conventions, annotations and bounded atomic updates; [details](RNA_SPECTRUM_SUPPORT.md) |
| CHEMISTRY/Tagger | `Tagger`, `TaggerOptions` | Complete source constructor, vector/spectrum extraction, append and maximum-charge setter, fixed/variable mass dictionary and I/L expansion; [details](TAGGER_SUPPORT.md) |
| CHEMISTRY/AdductInfo | `AdductInfo` | Full parser, component construction, electron-aware mass/m/z and mono/average shifts, signed formula compatibility and equality; [details](ADDUCT_SUPPORT.md) |
| CHEMISTRY/DecoyGenerator | `DecoyGenerator` | Full protein/peptide reversal and shuffling surface, seeded MT19937-64, source cache and positional enzyme rules, bounded overlapping unspecific products and atomic state changes; [details](DECOY_GENERATION_SUPPORT.md) |
| CHEMISTRY/SpectrumAnnotator, IonNaming | `SpectrumAnnotator`, `ion_naming` | All three annotation/statistics entrypoints and four name helpers; source duplicate/flag/statistical conventions with atomic updates, safe small-list quartiles and cumulative generation/alignment budgets; [details](SPECTRUM_ANNOTATION_SUPPORT.md) |
| CHEMISTRY/ISOTOPEDISTRIBUTION | `IsotopeDistribution`, `CoarseIsotopePatternGenerator`, `FineIsotopePatternGenerator`, `FineIsotopeIterator` | Coarse convolution/enrichment/averagine/conditional fragments plus bounded fine configurations, source thresholds and coverage, ordered raw streams and custom binary64 populations; [coarse](ISOTOPE_SUPPORT.md) and [fine](FINE_ISOTOPE_SUPPORT.md) |
| CHEMISTRY/ProteaseDB, ProteaseDigestion, EnzymaticDigestion | `ProteaseDB`, `Protease`, `ProteaseDigestion` | All 33 pinned enzymes, metadata and native rules; full/semi/nonspecific modified or unmodified products, ranges/counts, validity and missed-cleavage checks; [details](DIGESTION_SUPPORT.md) |
| PROCESSING/SCALING | `Normalizer`, `SqrtScaler`, `RankScaler` | Maximum/TIC normalization, negative-clamping square root, source rank behavior |
| PROCESSING/FILTERING | `ThresholdMower`, `NLargest`, `processing::window_mower::WindowMower` | Inclusive threshold, top N, sliding/jumping local selection and aligned arrays; [window filtering](WINDOW_MOWER_SUPPORT.md) |
| PROCESSING/SMOOTHING, BASELINE | Gaussian, Savitzky–Golay, morphology | Spectrum/chromatogram smoothing and ten morphology operations; [details](SIGNAL_PROCESSING.md) |
| PROCESSING/CENTROIDING/PeakPickerHiRes | `PeakPickerHiRes`, `CubicSpline2d`, `SignalToNoiseEstimatorMedian` | Spectrum/chromatogram centroiding, S/N, spacing/missing-flank rules, FWHM, boundaries and experiment selection; [details](PEAK_PICKING_SUPPORT.md) |
| PROCESSING/CENTROIDING/PeakPickerIterative | `processing::iterative` | HiRes seeds, iterative raw-sample refinement, original-seed suppression priority, integrated intensities, exact regions and rounded source arrays; [details](ITERATIVE_PICKING_SUPPORT.md) |
| PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative | `processing::mean_noise` | Three-pass histogram clipping, source fixed denominator, manual/global-standard-deviation/checked legacy-percentile ranges and sparse diagnostics; [details](MEAN_NOISE_SUPPORT.md) |
| ANALYSIS/OPENSWATH/PeakPickerChromatogram | `processing::chromatogram` | Legacy/corrected picking, smoothing, independent seed/boundary S/N, source overlap adjustment, raw sums and exact input regions; [details](CHROMATOGRAM_PICKING_SUPPORT.md) |
| MATH/MISC/EmgGradientDescent | `analysis::emg` | Source flank selection, literal analytic gradients, iRprop+ optimizer and truncated-side extrapolation; typed best-fit diagnostics, bounded work and numerical errors; [details](EMG_SUPPORT.md) |
| ANALYSIS/OPENSWATH/PeakIntegrator | `analysis::peak_integrator` | Intensity-sum, trapezoid and nonuniform Simpson integration, endpoint baselines and sampled shape metrics for spectra/chromatograms with optional EMG preprocessing; [details](PEAK_INTEGRATION_SUPPORT.md) |
| PROCESSING/DEISOTOPING/Deisotoper | `Deisotoper`, `AveragineDeisotoper` | Simple decreasing-intensity and Poisson/KL models, single-charge conversion, intensity sums and annotations; [simple](DEISOTOPING_SUPPORT.md) and [model](AVERAGINE_DEISOTOPING_SUPPORT.md) conventions |
| COMPARISON/SPECTRA, KERNEL/BinnedSpectrum | `comparison` | Absolute/ppm alignment and scoring, sparse bins, cosine/shared/agreeing, precursor, Zhang and Stein/Scott scores; [details](COMPARISON_SUPPORT.md) |
| PROCESSING/RESAMPLING | `LinearResamplerAlign`, `resample_to_grid` | Absolute spacing and supplied grids for intensity redistribution; no ppm spacing |
| FORMAT/DTAFile | `format::dta` | Buffered read/write with exact or legacy precursor-mass writer convention |
| FORMAT/FileTypes, FileNameUtils | `format::file_types` | Complete 73-format registry, properties, names, lexical suffix helpers and source dialog filters; source properties do not imply Rust reader support |
| FORMAT/FileHandler | `format::FileHandler` | Native experiment/feature/consensus dispatch, allowed types, bounded common content detection, gzip/bzip2 transport and atomic path output; other source dispatch families remain |
| FORMAT/FeatureXMLFile, ConsensusXMLFile | `format::featurexml`, `format::consensusxml` | Recursive feature geometry, source filters, typed metadata, processing history, portable identification chemistry and consensus group quantities; [featureXML](FEATUREXML_SUPPORT.md), [consensusXML](CONSENSUSXML_SUPPORT.md); ZIP input and inherited XSD validation remain |
| FORMAT/ModificationDefinitionIO | `format::modification_definitions` | Source escaped definition records, provenance-aware collection and owned local registration; [details](MODIFICATION_DEFINITION_IO_SUPPORT.md) |
| SYSTEM/File | `system::file` | Native filesystem, resource/configuration discovery and owned temporary resources, with explicit platform differences; [details](SYSTEM_FILE_SUPPORT.md) |
| FORMAT/MS2File, DTA2DFile | `format::ms2`, `format::dta2d` | Source text parsing, DTA2D ranges/storage/TIC and checked native MS2 writer; [details](TEXT_PEAK_LIST_SUPPORT.md) |
| FORMAT/OPTIONS/PeakFileOptions | `format::PeakFileOptions` | Complete source option state and Numpress configuration values; adapter execution remains separate; [details](PEAK_FILE_OPTIONS_SUPPORT.md) |
| CHEMISTRY/MASSDECOMPOSITION/MassDecomposition | `chemistry::MassDecomposition` | Complete count-container API and source cache quirks; [details](MASS_DECOMPOSITION_SUPPORT.md) |
| METADATA/Product hashing | `metadata::Product` and typed metadata | Equality-compatible hashes with signed-zero normalization; inherited source unit-state limits remain; [details](METADATA_HASH_SUPPORT.md) |
| FORMAT/FASTAFile | `format::fasta` | Complete file/stream reader/writer lifecycle, source lexical rules, bounded seek/progress and PEFF prologue skipping; [details](FASTA_SUPPORT.md) |
| FORMAT/MascotGenericFile | `format::mgf` | Buffered collection or streaming spectra; basic fields and unique extra key/value metadata; no Mascot submission |
| FORMAT/IdXMLFile | `format::idxml` | Optional bounded native identification read/write, typed UserParam metadata and run/evidence references; [detailed support](IDXML_SUPPORT.md) |
| FORMAT/HANDLERS/IndexedMzMLDecoder | `format::indexed_mzml` | Bounded footer discovery and both offset vectors, with source duplicate/order behavior and checked XML; [details](INDEXED_MZML_SUPPORT.md) |
| FORMAT/MzMLFile | `format::mzml` | Optional bounded subset; see [detailed support](MZML_SUPPORT.md) |

## Deliberate behavior differences

Native API:

- Rust owns values; exceptions become errors, output arguments become return values, parameters become typed configuration, and names use snake_case. There is no C++ ABI or source compatibility.
- Empty searches return `None`. Invalid indices, duplicate selection indices, inconsistent annotation arrays, nonfinite inputs and unsorted data for sorted-search APIs produce errors. Search validation is O(n), followed by binary lookup; repeated-query optimized sorted views are future work.
- Peak-position and intensity sorts preserve equal-value input order. No cached range invalidation is required because ranges are recomputed. Feature and consensus containers include legacy peptide/protein identification records. The separate `IdentificationData` sequence/provenance graph has stable owned references and checked atomic updates; observations, compounds, observation matches/groups, parent groups, referential cleanup and the sequence/evidence conversion bridge are implemented. Graph persistence and full legacy conversion remain unported. General ion-mobility operations and on-disk experiments remain unported; RT-binned TIC is available through the experiment summary API.
- Infallible low-level numerical summaries such as TIC assume valid peak data; call `validate` after direct public-field edits. TIC retains f32 behavior, including possible overflow for extreme inputs.
- Hulls preserve the source's scan-interval semantics, including its exact-interior-scan containment fallback. Outline-only hulls cannot answer containment queries. Checked feature-map operations reject ambiguous assigned IDs and subordinate trees deeper than 128; see [feature support](FEATURE_SUPPORT.md).

Processing:

- Chromatogram picking uses smoothed HiRes seeds and raw or smoothed boundary detection according to the selected method, then sums original intensities. Exact original indices supplement source f32 boundary annotations; seed and boundary noise options remain independent. Overlap adjustment follows the source and can leave overlapping regions. Crawdad is an unported external backend; see [picker support](CHROMATOGRAM_PICKING_SUPPORT.md).
- Peak integration selects existing samples inclusively without endpoint interpolation. Trapezoid retains f32 pair addition; Simpson retains source nonuniform spacing, even-count neighboring samples outside the bounds and its exact -1 subarea exclusion. Finite negative areas remain valid; undefined/nonfinite results are errors. Baseline and shape measurements follow sampled endpoint/threshold conventions; see [integration support](PEAK_INTEGRATION_SUPPORT.md).
- `ThresholdMower` defaults to the source 0.05 intensity cutoff. Use an explicit zero threshold to retain nonnegative zero-intensity peaks. The cutoff compares promoted f32 observations against f64.
- Maximum normalization of all-zero peaks is a no-op. A zero divisor with nonzero signed peaks is an error. C++ would divide by zero; the Rust port avoids NaN results and checks output intensity overflow.
- Maximum normalization still accepts all-negative input and uses its maximum, matching C++. The square-root scaler clamps negative intensities to zero without printing a global stderr warning.
- Rank scaling preserves the source's unusual rank offset, including the all-zero N+1 case. `NLargest` preserves the input order if at most N peaks exist; otherwise it returns descending intensity, with stable ties.
- Mutating filters validate or compute replacement values before committing changes. Experiment filters use temporary owned results so failures leave the experiment unchanged; this costs temporary memory.
- Absolute linear resampling redistributes intensities, not pointwise interpolation. Explicit-grid redistribution folds outside peaks into the nearest boundary; bounded `raster_align` instead excludes outside peaks, matching the corresponding C++ overloads.
- Resampling rejects nonempty auxiliary arrays because their semantics cannot be inferred. It preserves metadata and empty array placeholders. Invalid spacing/order and excessive grid sizes are errors. Reversed explicit bounds are errors rather than clearing the spectrum.
- HiRes picking defaults to disabled S/N estimation and preserves the source spline/extension rules. It requires distinct increasing coordinates and nonnegative intensities, copies acquisition metadata, reports omitted profile arrays, and supports selected intensity-weighted mobility arrays. FWHM and spline searches have finite iteration limits. The median noise estimator supports manual and mean-plus-standard-deviation histogram ranges; the source's unsafe percentile branch is not exposed.
- Simple deisotoping assumes centroid input and tries charge hypotheses from high to low. Defaults use 10 ppm, charges 1–3, and conversion to single charge. Simple clusters are disjoint by default, with explicit shared-isotope support; existing annotations follow retained/reordered peaks. The separate `AveragineDeisotoper` applies source threshold/top-N preprocessing, tests all charges and selects the longest acceptable Poisson/KL ladder, with highest charge breaking ties. It allows shared isotope extensions by default, reproducing the complete 104-peak source fixture, and offers disjoint selection explicitly. Both retain original input indices and aligned arrays; see [model support](AVERAGINE_DEISOTOPING_SUPPORT.md).

- EMG fitting preserves absolute-coordinate initialization, the literal source model and gradients, and best-iteration selection. Optional integration preprocessing extends the selected span. Parameters are typed records rather than an unaligned four-entry data array; invalid numerical expressions return errors. See [EMG support](EMG_SUPPORT.md).

Chemistry:

- Tagger preserves input index order, gap-relative ppm tolerance, strict source mass lookup, I/L branching and joint append sorting. Modification resolution uses two source lookup stages with deterministic provider-order ambiguity handling. Signed scalar placeholder masses remain private to its dictionary; [tag support](TAGGER_SUPPORT.md) documents all source edge cases and checked limits.
- Decoy generation retains source seed/cache history and positional anchoring, including N-terminal enzymes and overlapping unspecific products. Its private RNG retains Boost license terms. Checked modifications/empty-input policies and cumulative stateful limits are documented in [decoy support](DECOY_GENERATION_SUPPORT.md).

- The [chemistry support document](CHEMISTRY_SUPPORT.md) specifies charge/count grammar, exact isotope data, modifications and enzyme coverage.
- Peptide properties use parent residue codes. Charge/pI suppress annotated terminal groups; hydrophobicity and gas basicity ignore annotations. Gas basicity retains source ordinary arithmetic with documented overflow and empty-sequence numerical corrections; [property support](PEPTIDE_PROPERTIES_SUPPORT.md) states accepted alphabets, defaults and limits.
- Spectrum annotation preserves source array replacement, ppm last-match behavior and matched-only duplicates. Statistics retain their source flags and sort order; undefined enabled ratios return atomic errors. Generator and aligner counters are shared across candidates; coarse isotope convolution now also shares its limit across theoretical envelopes. Anonymous mass-tag text is immutable and shared across owned peptide slices. See [annotation support](SPECTRUM_ANNOTATION_SUPPORT.md).
- Negative-charge m/z uses absolute charge in the denominator. Formula charge follows OpenMS protonation semantics.
- Iridium uses the declared correct iridium table; the pinned C++ initialization accidentally passes rhenium's table.
- `AASequence` preserves unresolved B/Z/X and owns anonymous numeric tags. Formula, mono mass and average mass return `Result`; a known mass never implies a known formula. Source B/Z zero/negative placeholder masses and anonymous-tag formula omissions are replaced with explicit errors. Known formula and declared terminal-mass conventions remain tested; see [sequence support](SEQUENCE_SUPPORT.md).
- Digestion counts respect specificity and length filters. Unrestricted products use start-then-length order consistently for modified and unmodified sequences; the source APIs differ. Product validity supports explicit N-terminal M/MX and D|P allowances and retains full protein context for missed-cleavage counting; arbitrary regexes and runtime enzyme overrides remain unsupported. See [digestion support](DIGESTION_SUPPORT.md).
- The theoretical-spectrum generator has a larger scope than the b/y convenience method. Its defaults omit the first prefix ion; it can add source neutral losses, precursor peaks and coarse envelopes. Coarse envelopes retain the source's neutral-hydrogen-adduct mass convention. Unsupported terminal-loss/isotope combinations and chemically impossible loss formulas are handled explicitly; see [theoretical spectra](THEORETICAL_SPECTRA.md).

Identification analysis and XML:

- Peptide indexing implements the source ambiguity/mismatch and enzyme conventions with a bounded native matcher, complete evidence and explicit target/decoy policies. It reconstructs matched protein records as the source does. Both raw matching and `AASequence` peptide records support ambiguous letters; chemical calculations remain checked when composition or mass is unknown. The optimized parallel trie remains performance work; see [indexing support](PEPTIDE_INDEXING_SUPPORT.md).
- Basic protein inference preserves representative counts, the mean behind source `sum`, PEP conversion, default PSM-neighbor grouping and greedy resolution. Distinct source entry-point conventions are explicit. Finite native records reject source undefined/nonfinite retained results; see [inference support](PROTEIN_INFERENCE_SUPPORT.md).
- Conflict resolution moves rejected annotations to unassigned while preserving feature measurements. Origin splitting returns owned partitions and skipped records without mutating inputs, checks references within each run and rejects paths that would collide under source numeric keys; see [conflict support](ID_CONFLICT_SUPPORT.md) and [ripping support](ID_RIPPER_SUPPORT.md).

- Score switching preserves source score-name lookup and relative-value conventions, with atomic failures and conflict checks. HyperScore and Morpheus retain their overload-specific matching precision, charge rules and mean-error denominators; see [score support](SCORING_SUPPORT.md).
- Identification filters provide explicit keep/remove policies and run-aware evidence/protein/group cleanup. Filtering does not infer a protein model or rebuild all relationships automatically; see [filter support](ID_FILTER_SUPPORT.md).
- FDR distinguishes the actual legacy D/T and Basic pseudocount formulas. Source tie and lookup quirks are documented; undefined source probability-vector accesses are replaced with checked running means. These methods need appropriate target/decoy or posterior-probability input, not just arbitrary raw scores; see [FDR support](FDR_SUPPORT.md).
- idXML preserves the represented identification records and rejects unsupported state before writing. XML transport IDs are regenerated, while native run identity is retained via an explicit UserParam. Parsing and serialization have size limits; parsing is owned and bounded, not streaming. See [idXML support](IDXML_SUPPORT.md) for metadata and sequence limitations.

Text formats:

- mzML now preserves named auxiliary float/integer/ASCII string arrays, including empty strings and empty placeholders. Exact integer range, byte/element/array budgets and preflight validation are checked. Semantic auxiliary CV types and array metadata/units beyond the native representation remain unsupported; see [mzML support](MZML_SUPPORT.md).

- DTA's default writer is the exact inverse of the reader using the proton mass constant. `MassConvention::LegacyOpenMS` reproduces the C++ writer's 1.0-Da approximation.
- DTA, MGF and mzML writers reject attached peptide identifications before writing, because the supported peak-file representations cannot round-trip those records.
- DTA rejects multiple precursors. DTA and MGF reject nondefault precursor acquisition metadata before output. DTA/MS2 names are not encoded by the format.
- MGF retains unique extra fields as string metadata and applies global fields. Repeated fields (including repeated SEQ), multiple charge hypotheses, RT ranges, annotations after peak pairs, invalid MS levels and nonfinite values are rejected. Metadata keys are canonicalized to uppercase on read.
- MGF represents one precursor per spectrum, regenerates `index=N` native IDs on read, and does not store auxiliary arrays or full acquisition metadata. Chromatograms are rejected on write.
- FASTA preserves UTF-8 sequence bytes after removing space, tab, CR and LF, including source annotations, digits and semicolons; alphabet validation happens in chemistry. Source header spacing and PEFF prologue handling are retained. Empty/malformed entries are errors. Header descriptions are not interpreted as structured PEFF metadata; see [FASTA support](FASTA_SUPPORT.md).

Iterative processing preserves strict-next-sample seed association, binary32 centroid storage, original seed priority and the source’s asymmetric recenter search. Window selection retains duplicate-position/equality semantics and last-window quota rules. Mean-noise clipping retains the original window denominator through all three passes; its legacy percentile mode is explicitly checked and is not a corrected statistical percentile. Native resource limits, deterministic tie policies and atomic failures are described in the linked support documents.

Modified-peptide generation preserves the source’s distinct terminal handling in its one-modification and general paths, including native typed placements that ordinary sequence text cannot reconstruct. idXML rejects those unrepresentable placements before output. Normal generated residue variants retain chemistry, fragments and search definitions through identification interchange.

## Remaining scope

The full-library port remains in progress. The historical PSI-MOD dataset is not bundled in the default registry; the OBO reader supports caller-supplied PSI-MOD records and aliases. Other unimplemented areas include the remaining OpenMS metadata families; `IdentificationData` persistence and remaining legacy converter APIs; IsoSpec layered traversal and backend performance hints; RNA enzyme XML/arbitrary regex support; other centroiding families; calibration; feature finding/grouping and broader feature/consensus processing; landmark discovery and alignment/scoring beyond the documented models and spectrum comparisons; database search; probabilistic protein inference; the complete IDBoostGraph API; optimized peptide-indexing throughput and additional identification scorers/FDR overloads; broader quantification and OpenSWATH workflows; imaging; QC; ML/solvers; mzIdentML/pepXML and other formats; Arrow/Parquet; and vendor raw readers.

These are real unimplemented areas, not hidden C++ fallbacks. The [repository analysis](REPOSITORY_ANALYSIS.md) gives a dependency-ordered roadmap and acceptance criteria. Full compatibility requires independent C++ differential runs and representative workflow benchmarks in addition to source-derived Rust tests.

## Runtime reporting and map helpers

- `concept::log_stream` supplies owned log routing, prefixes, duplicate suppression, callbacks and thread-local streams; see [logging support](LOG_STREAM_SUPPORT.md) for platform and shutdown boundaries.
- `concept::progress_logger` supplies all core progress modes/backends and CPU/wall timing; see [progress support](PROGRESS_LOGGER_SUPPORT.md).
- Public modification-definition collection now accepts complete feature/consensus maps, including nested subordinates. INI path input now recognizes gzip/bzip2; ZIP remains open.

## Identification file paths

Native [idXML path loading/publication and identification dispatch](IDENTIFICATION_PATH_SUPPORT.md)
are available, including caller-owned chemistry, document IDs, source output
extension checks and plain bytes regardless of compression suffix.

## Experiment summaries

[Checked experiment summaries](EXPERIMENT_SUMMARY_SUPPORT.md) add source TIC
binning, chromatogram and combined ranges, aligned chromatogram sorting, total
peak counts, level/zero-intensity queries and spectrum-only array clearing.

## Mass decomposition

The [native mass decomposition algorithm](MASS_DECOMPOSITION_ALGORITHM_SUPPORT.md)
supports all five source settings, atomic configuration/append and bounded
source-ordered composition enumeration. [Public IMS integer/real solvers](IMS_DECOMPOSER_SUPPORT.md) now provide
existence, single/all compositions, counts and constrained real-mass queries.

## Scientific mzML loading

[Explicit scientific loading](MZML_LOAD_OPTIONS_SUPPORT.md) now consumes the
supported PeakFileOptions filters/sorting and preserves aligned arrays with
26 canonical binary roles. Source-supported headers and metadata-only reading are implemented; consumer paths remain open; [Numpress transport](MZML_NUMPRESS_SUPPORT.md) is implemented.


## IMS foundations, peak traversal and raw compression (2026-09-11)

Native [IMS isotope distributions and elements](IMS_ISOTOPE_SUPPORT.md) and
[alphabets with replaceable text parsers](IMS_ALPHABET_SUPPORT.md) now cover their
reviewed public operation sets. Per-call configuration replaces isotope globals;
container ownership, checked failures and portable parsing are documented.

[Area traversal](AREA_ITERATION_SUPPORT.md), [peak indices](PEAK_INDEX_SUPPORT.md)
and [filtered bulk peak exports](PEAK_DATA_SUPPORT.md) provide borrowed selection
and native numeric output. Source MS-level narrowing and f64/f32 RT grouping are
retained. [Unfiltered plain/rich 2D import/export](EXPERIMENT_2D_SUPPORT.md) is now available.
Scan-mobility overloads and Feature mass-trace specialization remain separate.

[Raw Numpress](MSNUMPRESS_SUPPORT.md) implements linear, PIC, SLOF and Safe codecs
and fixed-point helpers. Its 295 executed C++ reference cases are narrow codec
evidence, not a full SDK comparison. The [base64/zlib wrapper](MSNUMPRESS_CODER_SUPPORT.md) also provides source
estimation, verification and fallback diagnostics. [mzML Numpress transport](MZML_NUMPRESS_SUPPORT.md)
accepts all six source CV modes and provides preflighted writing with ordinary
fallback. Full metadata and PeakFileOptions writer execution remain open; normal
builds require no C++ compiler.

## Unique IDs and two-dimensional values

[Unique IDs and UUID generation](UNIQUE_ID_SUPPORT.md) use an owned source-compatible
MT19937-64 stream and a reusable native value trait. [Plain, mobility and rich 2D
peaks](PEAK2D_SUPPORT.md) preserve source numeric values and dimension labels; rich
values own metadata and implement the shared ID trait. Full inherited metadata
and generic range/container behavior remain separately tracked.

## Sequence coverage

[SequenceCoverage](SEQUENCE_COVERAGE_SUPPORT.md) computes the union of every exact
unmodified peptide occurrence, including overlaps and repeats, as a protein
coverage percentage. Modifications do not alter symbol matching; this utility
does not use identification scores, digestion rules or peptide evidence positions.

## Acquisition records and chromatogram conversion

[ChromatogramTools](CHROMATOGRAM_TOOLS_SUPPORT.md) implements both source conversion
operations with explicit grouping, first-record metadata and removal behavior.
Spectra/chromatograms now own acquisition records and shared processing handles.
[Processing propagation](PROCESSING_ACQUISITION_SUPPORT.md) retains these fields;
[mzML guards](MZML_ACQUISITION_GUARDS.md) reject their unrepresented transport.
Full settings bridges, comments and scan mobility remain separate work.
