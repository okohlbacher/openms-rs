# OpenMS for Rust

A native Rust port of selected [OpenMS4-core](https://github.com/okohlbacher/OpenMS4-core/tree/bc9cc12514c768385ce121d6ca4bb710fe1983c4) functionality: spectra, features, chemistry, identification records and analysis, retention-time transformations, common processing, and file interchange.

The target is a feature-complete native Core SDK that TOPP tools can be ported against. The [completion ledger](docs/CORE_SDK_COMPLETION.md) accounts for every registered public SDK header and maps direct dependencies from 146 local TOPP sources. It distinguishes reviewed APIs from partial and unverified coverage; full SDK and tool parity are still outstanding.

**This port is in progress and does not yet replace the full OpenMS library.** The current reduced SDK contains 807 physical include-directory headers and about 468,000 lines of first-party runtime code. The [SDK update](docs/CORE_SDK_UPDATE.md) records the current target, `bc9cc12`, and the product backends removed from its scope. This crate has its own Rust API, no C++ bindings, and no C++ build dependency. Native [FASTA lifecycle](docs/FASTA_SUPPORT.md), [aggregation/XICs](docs/EXPERIMENT_AGGREGATION_SUPPORT.md), [TIC and experiment summaries](docs/EXPERIMENT_SUMMARY_SUPPORT.md), [mzML Product transport](docs/MZML_PRODUCT_SUPPORT.md) and [indexed-mzML offsets](docs/INDEXED_MZML_SUPPORT.md) are now available. Native [peak-file option values](docs/PEAK_FILE_OPTIONS_SUPPORT.md), [metadata hashing](docs/METADATA_HASH_SUPPORT.md) and [mass-decomposition records](docs/MASS_DECOMPOSITION_SUPPORT.md) and the [native solver](docs/MASS_DECOMPOSITION_ALGORITHM_SUPPORT.md) are also available. The [repository analysis](docs/REPOSITORY_ANALYSIS.md) explains the source architecture, dependencies, and path toward broader coverage.

Native [IMS isotope/element operations](docs/IMS_ISOTOPE_SUPPORT.md) and [alphabets/parsers](docs/IMS_ALPHABET_SUPPORT.md), [area traversal](docs/AREA_ITERATION_SUPPORT.md), [peak indices](docs/PEAK_INDEX_SUPPORT.md), [filtered bulk peak export](docs/PEAK_DATA_SUPPORT.md), and [raw Numpress codecs](docs/MSNUMPRESS_SUPPORT.md) with the [base64/zlib wrapper](docs/MSNUMPRESS_CODER_SUPPORT.md) are available. Raw Numpress has 295 executed C++ reference cases; [mzML Numpress transport](docs/MZML_NUMPRESS_SUPPORT.md) now supports automatic reading and configured writing with ordinary fallback.

Native [unique IDs and UUIDs](docs/UNIQUE_ID_SUPPORT.md), [2D peak values](docs/PEAK2D_SUPPORT.md), [plain/rich experiment conversion](docs/EXPERIMENT_2D_SUPPORT.md) and [public IMS integer/real decomposers](docs/IMS_DECOMPOSER_SUPPORT.md) and [sequence coverage](docs/SEQUENCE_COVERAGE_SUPPORT.md) are available.

Native [spectrum–chromatogram conversion](docs/CHROMATOGRAM_TOOLS_SUPPORT.md) retains source-selected acquisition records. Spectra/chromatograms now own acquisition settings and shared processing handles; public struct literals need the new fields or `..Default::default()`. [Processing](docs/PROCESSING_ACQUISITION_SUPPORT.md) preserves them, and [mzML guards](docs/MZML_ACQUISITION_GUARDS.md) reject settings that cannot yet be serialized.

The [FORMAT wave](docs/FORMAT_WAVE_SUPPORT.md) adds mzTab/mzTab-M, mzXML, mzData, pepXML, mzIdentML, Mascot, qcML, Percolator, MSstats and transformation XML adapters. Its remaining exporter, validator and lookup gaps are recorded individually.

The optional `sqlite` feature supplies the [connector](docs/SQLITE_CONNECTOR_SUPPORT.md) and [SWATH lookups](docs/MZML_SQLITE_SWATH_SUPPORT.md). The `sqmass` feature adds the [SQLite spectrum/chromatogram handler](docs/MZML_SQLITE_HANDLER_SUPPORT.md), including bounded reads and transactional writes. SqMassFile, streaming/access adapters, OSW and OMS remain in the [next storage stages](docs/PORTING_WAVES.md).

Recent integrations add numerical fitters, scalar rustfft transforms, spectrum comparisons and system/process utilities. HTTP support uses optional `network`; the default `parallel` feature retains the documented serial/parallel bitwise contract. See the [porting status](docs/PORTING_STATUS.md) and per-header ledger for evidence and remaining API gaps.

## What works

[MassTrace](docs/MASS_TRACE_SUPPORT.md) provides owned trace peaks, cached centroids,
raw/smoothed quantification, FWHM and hull calculations.
[Mass-trace detection](docs/MASS_TRACE_DETECTION_SUPPORT.md) now supplies both source
run operations, mobility-aware trace growth and all eleven settings.
[Elution-peak detection](docs/ELUTION_PEAK_DETECTION_SUPPORT.md) adds trace splitting,
smoothing, extrema, width filtering and noise/SNR calculations. The complete
[constants and metadata-key collection](docs/CONSTANTS_SUPPORT.md) is also available.
[mzML record settings](docs/MZML_SETTINGS_SUPPORT.md) preserve spectrum scan modes,
polarity, zoom, scan windows and Product lists, plus chromatogram types.
[Ordered acquisition records](docs/MZML_ACQUISITION_SUPPORT.md) retain spectrum
scan identifiers, combination methods, scalar metadata and source-file references.
The [monosaccharide database](docs/MONOSACCHARIDE_SUPPORT.md) contains the complete
source collection of 24 records and 12 aliases.
[ProForma annotation data and text writers](docs/PROFORMA_SUPPORT.md) retain every
annotation variant, with 160 executed C++ writer comparisons.
[Text parsing and structured errors](docs/PROFORMA_PARSER_SUPPORT.md) cover both
source grammars, with 476 executed C++ parser/output/error comparisons.
[Protein-run helpers](docs/PROTEIN_RUN_SUPPORT.md) cover inference metadata, settings
export, singleton groups and metadata-only copying.
[Run mapping](docs/RUN_MAPPING_SUPPORT.md) resolves identification results to
ordered source-file paths, including merged runs and legacy path metadata.
[Feature hypotheses](docs/FEATURE_HYPOTHESIS_SUPPORT.md) supply borrowed isotope groups,
trace summaries, hulls and chromatogram exports. Optional [ProForma JSON](docs/PROFORMA_JSON_SUPPORT.md)
preserves the complete annotation schema.

[mzML count-only reading](docs/MZML_COUNTS_SUPPORT.md) supplies declared and filtered
record counts without decoding peak arrays.

[Contact and chromatography metadata](docs/EXPERIMENT_METADATA_SUPPORT.md) provide
complete ContactPerson, HPLC and Gradient records and operations.
[Document identity](docs/DOCUMENT_IDENTIFIER_SUPPORT.md) retains source equality,
path spelling and independent bounded content-type selection.

| Area | Implemented |
| --- | --- |
| Spectra and experiments | Peaks, chromatograms, precursors, aligned annotation arrays, sorting, selection, checked nearest/bound searches, ranges, base peaks, TIC and RT filtering |
| Features and geometry | Feature/consensus containers, checked maps and IDs, scan-envelope hulls, containment, consensus means and decharge summaries, [the identification surface of the feature containers](docs/FEATURE_IDENTIFICATION_SUPPORT.md): annotation state, checked peptide-identification sorting, primary IDs, observation-match sets and atomic reference translation |
| Chemistry | 84 element/isotope tables, formulas and masses, coarse isotope patterns/averagine, [fine isotope configurations, streaming and custom abundances](docs/FINE_ISOTOPE_SUPPORT.md), [named/numeric peptide annotations and unresolved residues](docs/SEQUENCE_SUPPORT.md), [33-enzyme full/semi/nonspecific digestion](docs/DIGESTION_SUPPORT.md), [fixed/variable modification generation](docs/MODIFIED_PEPTIDES_SUPPORT.md) and [named/anonymous definition sets and mass matching](docs/MODIFICATION_DEFINITIONS_SUPPORT.md), caller-owned registries and OBO/crosslink lookup, [charge/pI, hydrophobicity, amino-acid indices and gas basicity](docs/PEPTIDE_PROPERTIES_SUPPORT.md), [terminal/internal spectra, immonium ions, activation presets and compact mass ladders](docs/THEORETICAL_SPECTRA.md) with losses, precursors and annotations; [opt-in source-precision coarse isotope patterns, the source `trimLeft` and bounding-box predicates](docs/ISOTOPE_SOURCE_PRECISION_SUPPORT.md) for FeatureFinderAlgorithmPicked |
| Configuration and utilities | [Owned logging streams](docs/LOG_STREAM_SUPPORT.md), [progress reporting and timing](docs/PROGRESS_LOGGER_SUPPORT.md), [Typed parameter values](docs/PARAM_VALUE_SUPPORT.md), [hierarchical defaults, restrictions, updates and CLI parsing](docs/PARAM_SUPPORT.md), [default-parameter lifecycles](docs/DEFAULT_PARAM_HANDLER_SUPPORT.md), [OpenMS INI read/write](docs/PARAMXML_SUPPORT.md), [literal text/CSV helpers](docs/TEXT_CSV_SUPPORT.md) and [list utilities](docs/LIST_UTILS_SUPPORT.md); [filesystem, runtime paths and owned temporary resources](docs/SYSTEM_FILE_SUPPORT.md) |
| Processing | Normalization, threshold/top-N/window filtering, scaling, linear resampling, Gaussian/Savitzky–Golay smoothing, morphological baseline correction, median/iterative-mean noise estimates, HiRes and iterative centroiding, simple and Poisson/KL deisotoping with charge conversion; [PeakPickerHiRes with the source parameter contract](docs/PEAK_PICKING_SUPPORT.md) and opt-in source acceptance of degenerate input; [both signal-to-noise estimator headers](docs/SIGNAL_TO_NOISE_SUPPORT.md) are complete, with all three histogram ranges - the percentile range computed on the domain where the source is defined and refused outside it - the source build's binning conversions, the three warnings, optional progress reporting and `estimate_noise_from_random_scans` with an explicit seed and libstdc++'s `nth_element`, with the spectrum loop parallelised behind the `parallel` feature and byte-identical output at every worker count (1.88x faster than the C++ picker at 32 threads, [BENCHMARKS](docs/BENCHMARKS.md)); a reusable [cubic-spline fitter](docs/CUBIC_SPLINE2D_SUPPORT.md) that keeps its scratch across the millions of splines a centroiding run constructs, bit-identical by construction; [feature overlap filtering and FAIMS feature merging](docs/FEATURE_OVERLAP_FILTER_SUPPORT.md) in the source mode, with its quadtree |
| Chromatogram processing | Legacy/corrected peak picking, exact sample boundaries, raw intensity sums, time-weighted integration, baseline estimates, sampled shape metrics and EMG reconstruction, [millisecond-bucket chromatogram merging](docs/CHROMATOGRAM_MERGE_SUPPORT.md) |
| Spectrum comparisons | Absolute/ppm alignment, alignment scores, sparse bins, cosine/shared/agreeing scores, precursor, Zhang and Stein/Scott scores |
| RNA chemistry | [Nucleotide records, full pinned registry, TSV/JSON providers and nucleic-acid sequences](docs/RNA_SUPPORT.md), terminal/sulfur linkages, all source fragment formulas and masses, owned custom chemistry and checked slicing |
| Molecular adducts | [Neutral mass/m/z conversion, monomer/dimer notation, electron-aware shifts and formula compatibility](docs/ADDUCT_SUPPORT.md) |
| Decoy sequences | [Whole-protein and peptide reversal, seeded peptide shuffling and deterministic variants](docs/DECOY_GENERATION_SUPPORT.md), with source enzyme/cache conventions and checked work limits |
| Sequence tags | [Residue-mass tags from measured spectra](docs/TAGGER_SUPPORT.md), charge hypotheses, fixed/variable modifications, I/L alternatives and bounded atomic append |
| Spectrum annotation | [Fragment labels, hit peak annotations and match statistics](docs/SPECTRUM_ANNOTATION_SUPPORT.md), source ion-name parsing/display, shared calculation budgets and atomic updates |
| Precursor purity | [Scalar isolation scores, SPS matching, fuzzy scan purity and RT interpolation](docs/PRECURSOR_PURITY_SUPPORT.md), parent-scan lookup and checked experiment scoring |
| Experimental design | [Design records, fraction/label/sample mappings, condition and prefractionation grouping, basename filtering and consensus/feature/identification derivation](docs/EXPERIMENTAL_DESIGN_SUPPORT.md), with the tab-separated reader in both source table layouts |
| Metadata and identifications | Typed values and CV terms, acquisition settings, peptide/protein evidence and scores, protein coverage/modifications, attachments to spectra and feature maps |
| Identification graph | [Owned sequence and provenance records](docs/IDENTIFICATION_GRAPH_SUPPORT.md), stable typed references, score histories, parent/match groups, coverage, atomic registration/merge/copy, [referential cleanup and filtering](docs/IDENTIFICATION_CLEANUP_SUPPORT.md), legacy sequence/evidence conversion and RNase integration |
| Identification analysis | Score categories/switching, HyperScore/Morpheus fragment scores, peptide/protein filters, target/decoy FDR and q-values, picked proteins and probability estimates; peptide-to-protein indexing, basic protein inference/grouping, feature/spectrum conflicts and file-origin splitting |
| Retention-time models | Linear regression with coordinate weights, linear/natural-cubic interpolation, robust LOWESS, inverse refits, deviations, residual windows and atomic application to experiments/feature maps |
| File interchange | [Complete file-type registry and shared native experiment dispatch](docs/FILE_HANDLING_SUPPORT.md), streaming FASTA/MGF, DTA — whose writer reproduces the source's own two 15-digit numeric rules, so a 36,443-file extraction is byte-identical to the C++ tool's — and [MS2/DTA2D](docs/TEXT_PEAK_LIST_SUPPORT.md), optional bounded mzML with [scientific loading filters and canonical arrays](docs/MZML_LOAD_OPTIONS_SUPPORT.md), [reusable parameter groups](docs/MZML_PARAM_GROUPS_SUPPORT.md), precursor acquisition and annotation arrays, [featureXML](docs/FEATUREXML_SUPPORT.md) with [size-derived rather than fixed read ceilings](docs/FEATUREXML_SCALE_SUPPORT.md), so a 2 GB, 800,000-feature map loads, [consensusXML with protein quantities](docs/CONSENSUSXML_SUPPORT.md), and [idXML with native file paths and identification dispatch](docs/IDENTIFICATION_PATH_SUPPORT.md) with portable custom chemistry; gzip/bzip2 file transport, [indexed-mzML random access](docs/INDEXED_MZML_HANDLER_SUPPORT.md) that reads one spectrum or chromatogram at its index offset without loading the file; [FileHandler type detection by name and content, with option-taking experiment and feature loaders](docs/MZML_MOBILITY_SUPPORT.md), and the source [peak-type estimator](docs/PEAK_TYPE_ESTIMATOR_SUPPORT.md) on raw peaks; [source-compatible reading of dangling mzML software and data-processing references](docs/MZML_HEADER_SUPPORT.md); a [FileInfo library preview](docs/FILE_INFO_SUPPORT.md) for DTA, DTA2D, mzML and featureXML with `-m`, `-p` and `-s` in text and TSV |
| RNA processing | [Fourteen RNases, owned enzyme records and atomic graph registration](docs/RNASE_SUPPORT.md), [fixed/variable RNA modification generation](docs/RNA_MODIFICATION_SUPPORT.md), and [single/multiple annotated RNA spectra](docs/RNA_SPECTRUM_SUPPORT.md), with source mass/charge and terminal conventions |
| Ranges and predicates | [Range algebra and on-demand container ranges](docs/RANGES_SUPPORT.md) for spectra, chromatograms, mobilograms and the three experiment roles, with per-MS-level lookup and mobility; [spectrum and peak predicates](docs/RANGE_UTILS_SUPPORT.md) with the source reverse flag |
| Spectrum helpers | [Data-array lookup by name, intensity rebasing, position-unique merging and metadata copy](docs/SPECTRUM_HELPER_SUPPORT.md), generic over spectra and chromatograms |
| Geometry primitives | [Const-generic DPosition, DIntervalBase and DRange](docs/DPOSITION_SUPPORT.md) with the source empty sentinel, normalisation, intersection classification and extension rules |
| Ion mobility | [Drift time, ion-mobility arrays and frame rasterisation on spectra](docs/SPECTRUM_MOBILITY_SUPPORT.md): the `-1` drift-time sentinel, name-based IM-array detection over the nine pinned PSI-MS terms, mobility sorting, presorted chunked merging and bounded frame rasterisation; [FAIMS compensation voltages and FAIMS-CV identification filtering](docs/FAIMS_HELPER_SUPPORT.md); [spectrum- and scan-level mobility, the representation reset and unit-bearing ion-mobility arrays in mzML](docs/MZML_MOBILITY_SUPPORT.md); [splitting an experiment by FAIMS compensation voltage](docs/IM_DATA_CONVERTER_SUPPORT.md) |
| Mass-spectrometry imaging (imzML) | [The two-file imzML format](docs/IMZML_HANDLER_SUPPORT.md): the `.imzML` index and geometry over a companion `.ibd` with checked byte ranges, both continuous and processed storage, a [dataset writer](docs/IMZML_WRITER_SUPPORT.md), the [file adapter](docs/IMZML_FILE_SUPPORT.md) and an [on-disc pixel reader](docs/ON_DISC_IMZML_SUPPORT.md) with ion-image extraction |
| Targeted proteomics (SRM/MRM) | [Peak groups and transition groups](docs/MRM_SUPPORT.md): per-transition and precursor features, OpenSwath score sets, the three parallel transition/chromatogram maps with their consistency checks, both subsets and library-intensity/quality summaries |
| Feature and consensus maps | [Container operations of the two map types](docs/MAP_OPERATIONS_SUPPORT.md): annotation statistics, map arithmetic with caller-owned unique-ID resolution, row and column appends, splitting a consensus map back into feature maps, the primary MS run path on both containers and the two stream layouts; [conversions between peak, feature and consensus containers](docs/CONVERSION_HELPER_SUPPORT.md) with the source's asymmetric column-header sizes preserved |
| TOPP tools | [The TOPP command-line framework and the first executable tool](docs/TOPP_CLI_SUPPORT.md): registration, INI resolution, usage text, validation and the source exit codes; eight executable tools — `DTAExtractor`, `MzMLSplitter`, `MapNormalizer`, `SpectraFilterWindowMower`, `BaselineFilter`, and now [`PeakPickerHiRes`](docs/TOPP_PEAK_PICKER_HI_RES_SUPPORT.md) (in-memory centroiding, the `algorithm` subsection and `-write_ini`), [`FileInfo`](docs/TOPP_FILE_INFO_SUPPORT.md) (peak-file and featureXML reports with `-m`, `-p`, `-s`) and [`FeatureFinderCentroided`](docs/TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md) (features end to end, including FAIMS input: one run per compensation voltage, the `FAIMS_CV` annotation and the cross-voltage merge, with the two defects of the C++ merge corrected) — each reproducing its upstream tests against retained C++ output, with algorithm subsections. [`-threads`](docs/TOPP_THREADS_SUPPORT.md) reaches the five wave-2 tool bodies through a real worker pool, with a non-positive count meaning every processor as the executed C++ does, and no thread count changes an output byte; `PeakPickerHiRes` is the first tool whose pool does real work, and it is the sixth tool the shared `-threads` suite checks. Five of the eight tools build the pool and leave it idle, which every 32-thread figure in [BENCHMARKS](docs/BENCHMARKS.md) is labelled with. All eight tools agree with the C++ Release build on the data of every one, seven of them measured on full-size instrument data and `FeatureFinderCentroided` on a documented 4,000-spectrum subset because neither implementation finishes the full run: three of them are byte-equal and the five mzML writers are bitwise equal on every decoded array, differing only in identifier spelling, processing provenance and XML indentation the port does not write. Every tool follows TOPPBase's exit codes, including a bare invocation and unreadable, directory or `/dev/null` INI files, checked against the executed C++; [FuzzyStringComparator and FuzzyDiff](docs/FUZZY_STRING_COMPARATOR_SUPPORT.md) are ported as shared test support, with a decoded featureXML/mzML comparator; [the numeric text FileInfo prints](src/format/file_info/text_format.rs) (`StringUtils::number` and `toStr`, stream `%g` output and vector output) follows the executed C++, including `std::to_chars` ties to even, with its platform differences documented; `-write_ini` files, usage text and input and output format checks equal the executed C++ tools', and the four mzML tools attach their DataProcessing records |
| Feature finding (picked) | [FeatureFinderAlgorithmPicked helper structures](docs/FEATURE_FINDER_PICKED_HELPER_STRUCTS_SUPPORT.md): seeds, mass traces, intensity profiles and isotope patterns, the first step of the FeatureFinderCentroided chain; [the whole of FeatureFinderAlgorithmPicked](docs/FEATURE_FINDER_PICKED_SUPPORT.md), complete since wave 5: parameters, input checks, intensity, trace and isotope-pattern scores, pattern precalculation, seed selection, mass-trace extension, the isotope fit, the quality checks, cropping, abort bookkeeping and overlap resolution, so `run()` produces features end to end, plus the reusable instance (caller maps, accumulated aborts, the parameter surface, progress logging) and `write_debug` with `writeFeatureDebugInfo_` byte for byte. It is pinned tier 1 against the Linux x86_64 Release build: every source sort follows libstdc++'s introsort and 14.4.0 `stable_sort`, and the seed score, both trace fits and the isotope windows use that build's glibc 2.39 `powf`, `exp` and `log`, ported licence-clean from Arm optimized-routines, so the fits are bit for bit on every platform. It refuses only where the source has no reproducible answer - an out-of-bounds read or write, a loop that never ends, a process that terminates, output that depends on a heap address - and each such refusal records where and how the C++ process ends; the one accepted exception is the multi-thread abort race, which gives the single-thread result because parallel output must equal serial output. On a default x86-64 baseline build that fidelity costs 21 % at one thread, which `-C target-feature=+fma` removes with bitwise-identical output ([BENCHMARKS](docs/BENCHMARKS.md) §4); x86_64 builds set that flag by default since `port/fma-default` and require an FMA3-capable processor, refusing to start on any other with the rebuild command ([the FMA build flag](docs/FMA_BUILD_FLAG.md)); the [Gauss and EGH trace fitters](docs/TRACE_FITTER_SUPPORT.md) both models fit with, over a transcription of Eigen's Levenberg-Marquardt that matches the Linux x86_64 Release build bit for bit |
| Boost-compatible regular expressions | [One facade over `fancy-regex`](docs/BOOST_REGEX_SUPPORT.md) that gives Boost.Regex's answer or refuses at construction with `Unsupported`, never a different answer: 6.5 million compared cases over 134,075 patterns with 0 mismatches against a driver compiled from Boost 1.92, and 0 refusals in all sixteen OpenMS-derived pattern families |
| Examples | Read/filter/normalize/write spectra; stream FASTA and calculate tryptic peptide masses; inspect modified-peptide masses, isotopes and theoretical fragments; identify a synthetic modified peptide and compute protein coverage; pick and integrate chromatograms; reconstruct a cropped EMG peak; estimate profile noise, centroid and retain local peaks; enumerate modified digest products; stream natural and enriched isotope configurations; calculate peptide physicochemical properties; annotate measured fragments and inspect matching statistics; generate decoy FASTA; extract and match sequence tags; inspect charged RNA formulas and isotope patterns; digest RNA, enumerate variants and generate annotated fragments |

See the [coverage and differences](docs/PORTING_STATUS.md), [chemistry support](docs/CHEMISTRY_SUPPORT.md), and [mzML subset](docs/MZML_SUPPORT.md) before using this as a substitute for a C++ workflow. [Third-party crate decisions](docs/THIRD_PARTY_CRATE_DECISIONS.md) record, with measurements, which external C++ libraries are replaced by crates and which stay ported code. [Feature containers](docs/FEATURE_SUPPORT.md), [HiRes peak picking](docs/PEAK_PICKING_SUPPORT.md), [theoretical spectra](docs/THEORETICAL_SPECTRA.md), [simple deisotoping](docs/DEISOTOPING_SUPPORT.md), and [Poisson/KL deisotoping](docs/AVERAGINE_DEISOTOPING_SUPPORT.md) each have a defined supported subset.

[Typed metadata](docs/METADATA_SUPPORT.md), [identification records](docs/IDENTIFICATION_SUPPORT.md), and [retention-time models](docs/TRANSFORMATIONS_SUPPORT.md) document the data and alignment capabilities. [Score handling and fragment scores](docs/SCORING_SUPPORT.md), [identification filtering](docs/ID_FILTER_SUPPORT.md), [FDR calculations](docs/FDR_SUPPORT.md), and [idXML interchange](docs/IDXML_SUPPORT.md) describe supported identification workflows and source conventions. [Peptide indexing](docs/PEPTIDE_INDEXING_SUPPORT.md) links identified sequences to FASTA proteins; [basic protein inference](docs/PROTEIN_INFERENCE_SUPPORT.md) aggregates evidence and resolves groups. [Conflict resolution](docs/ID_CONFLICT_SUPPORT.md) and [file-origin splitting](docs/ID_RIPPER_SUPPORT.md) handle competing and merged identifications.

The [chromatogram picker](docs/CHROMATOGRAM_PICKING_SUPPORT.md) preserves source smoothing, seed and boundary conventions. The [peak integrator](docs/PEAK_INTEGRATION_SUPPORT.md) supports spectra and chromatograms, including the source's nonuniform Simpson averaging and sampled shape metrics. Optional [EMG fitting](docs/EMG_SUPPORT.md) reconstructs cropped peaks with the source iRprop+ optimizer and exposes typed fit diagnostics.

The [iterative picker](docs/ITERATIVE_PICKING_SUPPORT.md) refines HiRes seeds and reports exact input regions alongside the source’s rounded centroid and boundary arrays. [Window filtering](docs/WINDOW_MOWER_SUPPORT.md) supports sliding and jumping windows; the [iterative mean noise estimator](docs/MEAN_NOISE_SUPPORT.md) preserves the source’s three-pass clipping conventions.

IsoSpec layered traversal; arbitrary RNA enzyme regexes and XML import; other peak-picker families; feature finding/grouping; database search; probabilistic protein-inference engines; broader quantification and OpenSWATH workflows; vendor RAW formats; and Arrow/Parquet remain unimplemented. Identification-graph groups, referential cleanup and the bounded legacy sequence/evidence conversion bridge are implemented; graph persistence and the remaining converter APIs are outstanding.

The [mobility containers](docs/MOBILOGRAM_SUPPORT.md) provide checked mobilogram
search, sorting, selection and summaries. Generic `DataArray` values now retain
metadata and shared processing descriptions; existing struct literals need
`..Default::default()` or `DataArray::new`. Unsupported XML projections reject
these descriptions before output. [Integer mass weights](docs/IMS_WEIGHTS_SUPPORT.md)
and [mzML file operations](docs/MZML_PATH_SUPPORT.md) expose the next native SDK
utilities, including direct filtered loading and atomic compressed output.

## Use locally

Rust 1.85 or newer is required.

**On x86_64 this repository builds with `-C target-feature=+fma`**
(`.cargo/config.toml`), so a binary built here needs an FMA3-capable processor:
**Intel Haswell (2013), AMD Piledriver (2012) or newer**. A tool run on a
processor that has AVX but not FMA (Sandy Bridge, Ivy Bridge, Bulldozer) prints
what it needs and exits 12 rather than dying on an illegal instruction. On a
processor older than AVX entirely (before 2011) that refusal is best effort
rather than a guarantee: `+fma` implies AVX, so such a machine can fault on a
vector instruction that is not a fused multiply-add. The path to the check was
measured to be free of them for this compiler, which is a measurement and not a
property of the design -- see `docs/FMA_BUILD_FLAG.md`. To build for an older processor -- which changes no result and only
costs speed -- use exactly:

```sh
RUSTFLAGS="-C target-feature=-fma" cargo build --release --locked
```

`RUSTFLAGS` replaces the configured flags rather than adding to them, which is
why that one command is the whole opt-out. `aarch64` is unaffected, and so is a
project that depends on this crate by path from its own checkout, because Cargo
reads `.cargo/config.toml` from the working directory of the `cargo` invocation.
See [the FMA build flag](docs/FMA_BUILD_FLAG.md).

Add this local crate to a consuming project:

```toml
[dependencies]
openms = { path = "/absolute/path/to/OpenMS4-R" }
```

The `mzml`, `idxml`, `paramxml`, `featurexml`, `consensusxml`, `cv-mapping` and `rna-json` features are enabled by default. The XML features use Rust XML parsing. The independent `numpress` feature adds the base64/zlib wrapper and is enabled by mzML; raw Numpress remains available without features. The `file-compression` feature supplies gzip and bzip2 using Rust backends and is enabled by mzML, parameter XML, featureXML and consensusXML. `rna-json` adds serde_json for caller-supplied MODOMICS JSON. The embedded RNA registry and TSV reader remain available without that feature. Runtime reporting uses `chrono` for local timestamps and `cpu-time` for process CPU timing on Unix/Windows. The scientific core uses the small pure Rust `libm` library for EMG's complementary error function. Disable default features to use the scientific core, identification analysis and text readers without the XML, JSON and compression dependencies:

```toml
openms = { path = "/absolute/path/to/OpenMS4-R", default-features = false }
```

[ProForma AASequence conversion](docs/PROFORMA_CONVERSION_SUPPORT.md) supports both directions, all policies and diagnostics with an explicit registry. [CV mapping records and XML loading](docs/CV_MAPPING_SUPPORT.md) are available; `cv-mapping` is included by default and can be enabled independently. General CV/mapping validation uses the optional `semantic-validation` feature.

Enable ProForma JSON independently with `default-features = false, features = ["proforma-json"]`.

Enable featureXML or consensusXML independently with `default-features = false, features = ["featurexml"]` or `["consensusxml"]`. Enable only parameter XML with `default-features = false, features = ["paramxml"]`. Enable only idXML with `default-features = false, features = ["idxml"]` to use identification XML without the mzML codecs.

```rust
use openms::{MSSpectrum, Peak1D};
use openms::chemistry::AASequence;
use openms::processing::{Normalizer, SpectrumFilter};

fn main() -> openms::Result<()> {
    let mut spectrum = MSSpectrum::from_peaks(vec![
        Peak1D::new(100.0, 10.0),
        Peak1D::new(200.0, 40.0),
    ]);
    Normalizer::default().filter_spectrum(&mut spectrum)?;
    assert_eq!(spectrum.peaks[0].intensity, 0.25);
    assert_eq!(spectrum.find_nearest(199.0)?, Some(1));

    let peptide = AASequence::parse("DFPIANGER")?;
    println!("Neutral mass: {:.6} Da", peptide.mono_mass()?);
    println!("Doubly protonated m/z: {:.6}", peptide.mz(2)?);
    Ok(())
}
```

Run the bundled examples without supplying any input files:

```sh
cargo run --locked --example process_spectra
cargo run --locked --example digest_fasta
cargo run --locked --example peptide_analysis
cargo run --locked --example identify_peptides
cargo run --locked --example integrate_chromatogram
cargo run --locked --example fit_emg_peak
cargo run --locked --example process_profile
cargo run --locked --example enumerate_peptides
cargo run --locked --example stream_isotopes
cargo run --locked --example peptide_properties
cargo run --locked --example annotate_spectrum
cargo run --locked --example generate_decoys
cargo run --locked --example extract_tags
cargo run --locked --example analyze_rna
cargo run --locked --example process_rna
cargo run --locked --example identify_rna
```

The spectrum example uses the original OpenMS 121-peak DTA fixture, retains 14 peaks at threshold 10, and normalizes their TIC to one. The chromatogram example picks two synthetic peaks and prints raw intensity sums, time-weighted areas, baseline estimates and their signed differences. The EMG example fits a cropped synthetic trace and compares its observed and reconstructed areas with the known complete synthetic area. The profile example estimates noise for 401 samples, produces four iterative centroids and retains the two strongest peaks with their integration and width annotations. The enumeration example digests a small protein, applies fixed cysteine alkylation and lists peptide variants with up to two methionine oxidations and their masses. The isotope example prints a five-configuration natural glucose prefix and a threshold-selected enriched-carbon population. The property example reports charge at pH 7, isoelectric point, GRAVY and gas basicity at 500 K and 100 K. The annotation example matches eleven original source measurements and prints fragment labels and matching statistics. The decoy example writes two target proteins and two reproducible shuffled variants per target as FASTA. The RNA processing example digests AGUACG with RNase_T1, enumerates uridine modifications and prints annotated negative b/y fragments. The [RNA identification example](examples/identify_rna.rs) registers digest products with inclusive parent positions and processing history, then reports parent coverage. To use your own files:

```sh
cargo run --locked --example process_spectra -- input.mgf output.mgf
cargo run --locked --example digest_fasta -- proteins.fasta
cargo run --locked --example integrate_chromatogram -- chromatograms.mzML
```

The spectrum example also accepts DTA and mzML input. The output argument explicitly writes or replaces that MGF file. The FASTA example writes a TSV report to standard output; sequence offsets are zero-based and half-open. The peptide example prints modified-peptide masses, a coarse isotope envelope and annotated theoretical fragments. The identification example ranks candidates from a small embedded FASTA against a fixed synthetic spectrum, retains evidence and ion annotations, and computes coverage from the accepted peptide; its score is not a confidence/FDR estimate. The [identification pipeline test](tests/identification_pipeline.rs) separately demonstrates idXML → score switching/filtering → target/decoy q-values → reference cleanup → idXML with hand-computable confidence values. The [protein workflow test](tests/protein_workflow.rs) connects FASTA indexing, shared peptide evidence, basic inference, protein/peptide q-values, coverage, idXML round trips and file-origin partitions. Examples are demonstrations, not full TOPP command-line replacements.

## Validation

```sh
cargo test --locked --all-features --all-targets
cargo test --locked --all-features --doc
cargo test --locked --no-default-features
cargo test --locked --no-default-features --features idxml
cargo clippy --locked --all-features --all-targets -- -D warnings
cargo fmt --all -- --check
cargo doc --locked --all-features --no-deps
```

Append `--offline` when dependencies are already cached. The lockfile fixes the tested dependency versions. Tests use upstream assertions and fixtures, numerical invariants, round trips, malformed-input cases, and complete small workflows. The unmodified pinned raw Numpress C++ implementation was compiled for 295 differential cases; 287 also executed its decoder. The full C++ SDK was not built, and other source-derived tests are not live cross-language equivalence results. Cross-platform CI is provided; only the platforms actually listed in [validation notes](docs/VALIDATION.md) have been run locally.

## Design

Coordinates and retention times use `f64`; peak intensities use `f32` as in C++. Retention times are seconds. Algorithms return `Result` for invalid inputs and maintain valid annotation alignment. Centroiding reports profile arrays that lack an aggregation rule; resampling rejects populated auxiliary arrays. Public fields make construction straightforward; checked operations validate their preconditions, and ranges are computed on demand. Sorted searches currently spend O(n) validating order before the binary lookup.

The crate uses no unsafe Rust. Core chemistry is immutable and embedded; no OpenMS resource-directory setup is necessary. Unsupported scientific features are documented rather than represented by placeholder classes.

## Provenance and license

Current SDK target: `okohlbacher/OpenMS4-core` at `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Historical implementation and fixture provenance retains `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`, verified against GitHub on 2026-09-10. See the [SDK update](docs/CORE_SDK_UPDATE.md), [source provenance](SOURCE_PROVENANCE.json), the [historical source inventory](docs/source-inventory.json), and [fixture provenance](tests/data/README.md).

Implementation: BSD-3-Clause, with original OpenMS copyright and attribution in [LICENSE](LICENSE) and [AUTHORS](AUTHORS). Bundled UniMod-derived modification data: Design Science License, with complete source data and notices included. See [component licenses](LICENSES.md) and [data provenance](resources/modifications/README.md). The [RNA data notices](resources/rna/README.md) separately record the unresolved MODOMICS redistribution terms. This local package has not been published to crates.io.

The native metabolite feature-finding chain combines mass-trace detection,
elution peak detection and FeatureFindingMetabo. With centroided input:

```sh
cargo run --release --features mzml,featurexml --example find_metabolite_features -- input.mzML output.featureXML
```

See [feature-finding support](docs/FEATURE_FINDING_METABO_SUPPORT.md),
[fixed isotope predictors](docs/METABO_PREDICTOR_SUPPORT.md), and
[sample and instrument values](docs/EXPERIMENT_VALUES_SUPPORT.md).
This example exercises SDK operations; complete TOPP command behavior remains
subject to the [completion ledger](docs/CORE_SDK_COMPLETION.md).

Original C++ defects discovered during porting are tracked in
[OpenMS_CPP_ISSUES.md](OpenMS_CPP_ISSUES.md), with affected source files,
reproductions, evidence level and proposed upstream fixes.

[ProForma modification resolution](docs/PROFORMA_RESOLUTION_SUPPORT.md) now
fills shared chemistry handles using a caller-owned registry with atomic failure
handling. [Mass and m/z operations](docs/PROFORMA_MASS_SUPPORT.md) include
issue reports, availability checks and optional results. [Sequence conversion](docs/PROFORMA_CONVERSION_SUPPORT.md) is available;
[ProForma spectrum generation](docs/PROFORMA_SPECTRA_SUPPORT.md) now covers ordinary and crosslinked peptides.

[DateTime support](docs/DATETIME_SUPPORT.md) includes source-compatible parsing,
all seven formats, local/UTC clocks and checked Gregorian arithmetic.

[Experiment settings](docs/EXPERIMENTAL_SETTINGS_SUPPORT.md) now own sample,
instrument, chromatography, date, provenance and typed run metadata. Callers
using `MSExperiment::metadata` must migrate to `experiment.settings.metadata`.
Processing preserves the complete settings. [mzML header transport](docs/MZML_HEADER_SUPPORT.md)
now preserves source-supported contacts, instruments, software, source files and
processing histories. `DataProcessing.completion_time` now uses `DateTime` to
retain milliseconds. [PeptideEvidence](docs/PEPTIDE_EVIDENCE_SUPPORT.md) also
supports complete native value and hash-key operations.

[Controlled vocabularies](docs/CONTROLLED_VOCABULARY_SUPPORT.md) now provide
complete term records, cumulative OBO loading, hierarchy queries, typed XML
values and all five pinned source providers. All 9,254 terms and 16,852 name
aliases are checked against an independent source projection.

[General semantic CV validation](docs/SEMANTIC_VALIDATOR_SUPPORT.md) is available
with `features = ["semantic-validation"]`. It checks mapping rules, term names,
values and optional units, returning ordered errors and warnings. It preserves
documented source value conventions and does not perform full XSD validation.

[Crosslink spectrum generation](docs/THEORETICAL_XLMS_SUPPORT.md) now provides
all three source XLMS append operations and their full options. The owned
crosslink record preserves sequence identity. ProForma spectrum wrappers compose
these backends; other XLMS analysis classes remain separate work.

[Streaming mzML consumers](docs/MZML_CONSUMER_SUPPORT.md) support setup counts,
batched callbacks, early stopping and optional retention of modified records.
The consumer trait is available without XML features. Metadata-only record
loading now supports `fill_data=false` while retaining descriptor checks.

[Spectrum type queries and mzML centroid inspection](docs/MZML_CENTROID_SUPPORT.md)
use stored type, processing history and optional signal estimation. File inspection
counts by MS level with an explicit recognized-spectrum quota and preserved caller options.

[Isolation-target loading](docs/MZML_ISOLATION_SUPPORT.md) selects and filters by
the isolation window while retaining differing selected-ion metadata.
[Explicit mzML semantic validation](docs/MZML_VALIDATOR_SUPPORT.md) is available
with `mzml-validation`; it checks CV rules and values using the pinned vocabulary.

[Source mzML writer options](docs/MZML_WRITE_OPTIONS_SUPPORT.md) now support indexed
output, real SHA-1 checksums, primary precision, Numpress, zlib and TPP compatibility.
Stream writes preflight the complete output; path writes publish atomically.

[Binary whitespace normalization](docs/MZML_NORMALIZATION_SUPPORT.md) follows
`skip_xml_checks` while retaining checked XML, payload and resource validation.

[Explicit mzML XSD validation](docs/MZML_SCHEMA_SUPPORT.md) is available with
`features = ["mzml-schema"]`. This optional feature uses the original ordinary and
indexed schemas through libxml2 and requires its development library plus
libclang at build time. Default builds do not enable this dependency. Schema,
controlled-vocabulary, binary and index-integrity checks have separate contracts.

Spectrum and chromatogram metadata now uses typed values. See the
[API migration guide](docs/RECORD_METADATA_MIGRATION.md) and
[mzML transport details](docs/MZML_TYPED_TRANSPORT_SUPPORT.md) for optical spectra,
pressure/flow chromatograms, primary metadata and independent noise grids.
