# OpenMS for Rust

A native Rust port of selected [OpenMS4-core](https://github.com/okohlbacher/OpenMS4-core/tree/82ce5b373c97f934ffd9b1ffd80215ca66473d0b) functionality: spectra, features, chemistry, identification records and analysis, retention-time transformations, common processing, and basic file interchange.

The target is a feature-complete native Core SDK that TOPP tools can be ported against. The [completion ledger](docs/CORE_SDK_COMPLETION.md) accounts for every registered public SDK header and maps direct dependencies from 146 local TOPP sources. It distinguishes reviewed APIs from partial and unverified coverage; full SDK and tool parity are still outstanding.

**This port is in progress and does not yet replace the full OpenMS library.** The current reduced SDK contains 807 physical include-directory headers and about 468,000 lines of first-party runtime code. The [SDK update](docs/CORE_SDK_UPDATE.md) records the current target, `82ce5b3`, and the product backends removed from its scope. This crate has its own Rust API, no C++ bindings, and no C++ build dependency. Native [FASTA lifecycle](docs/FASTA_SUPPORT.md), [aggregation/XICs](docs/EXPERIMENT_AGGREGATION_SUPPORT.md), [TIC and experiment summaries](docs/EXPERIMENT_SUMMARY_SUPPORT.md), [mzML Product transport](docs/MZML_PRODUCT_SUPPORT.md) and [indexed-mzML offsets](docs/INDEXED_MZML_SUPPORT.md) are now available. Native [peak-file option values](docs/PEAK_FILE_OPTIONS_SUPPORT.md), [metadata hashing](docs/METADATA_HASH_SUPPORT.md) and [mass-decomposition records](docs/MASS_DECOMPOSITION_SUPPORT.md) and the [native solver](docs/MASS_DECOMPOSITION_ALGORITHM_SUPPORT.md) are also available. The [repository analysis](docs/REPOSITORY_ANALYSIS.md) explains the source architecture, dependencies, and path toward broader coverage.

Native [IMS isotope/element operations](docs/IMS_ISOTOPE_SUPPORT.md) and [alphabets/parsers](docs/IMS_ALPHABET_SUPPORT.md), [area traversal](docs/AREA_ITERATION_SUPPORT.md), [peak indices](docs/PEAK_INDEX_SUPPORT.md), [filtered bulk peak export](docs/PEAK_DATA_SUPPORT.md), and [raw Numpress codecs](docs/MSNUMPRESS_SUPPORT.md) with the [base64/zlib wrapper](docs/MSNUMPRESS_CODER_SUPPORT.md) are available. Raw Numpress has 295 executed C++ reference cases; [mzML Numpress transport](docs/MZML_NUMPRESS_SUPPORT.md) now supports automatic reading and configured writing with ordinary fallback.

Native [unique IDs and UUIDs](docs/UNIQUE_ID_SUPPORT.md), [2D peak values](docs/PEAK2D_SUPPORT.md), [plain/rich experiment conversion](docs/EXPERIMENT_2D_SUPPORT.md) and [public IMS integer/real decomposers](docs/IMS_DECOMPOSER_SUPPORT.md) and [sequence coverage](docs/SEQUENCE_COVERAGE_SUPPORT.md) are available.

Native [spectrum–chromatogram conversion](docs/CHROMATOGRAM_TOOLS_SUPPORT.md) retains source-selected acquisition records. Spectra/chromatograms now own acquisition settings and shared processing handles; public struct literals need the new fields or `..Default::default()`. [Processing](docs/PROCESSING_ACQUISITION_SUPPORT.md) preserves them, and [mzML guards](docs/MZML_ACQUISITION_GUARDS.md) reject settings that cannot yet be serialized.

## What works

| Area | Implemented |
| --- | --- |
| Spectra and experiments | Peaks, chromatograms, precursors, aligned annotation arrays, sorting, selection, checked nearest/bound searches, ranges, base peaks, TIC and RT filtering |
| Features and geometry | Feature/consensus containers, checked maps and IDs, scan-envelope hulls, containment, consensus means and decharge summaries |
| Chemistry | 84 element/isotope tables, formulas and masses, coarse isotope patterns/averagine, [fine isotope configurations, streaming and custom abundances](docs/FINE_ISOTOPE_SUPPORT.md), [named/numeric peptide annotations and unresolved residues](docs/SEQUENCE_SUPPORT.md), [33-enzyme full/semi/nonspecific digestion](docs/DIGESTION_SUPPORT.md), [fixed/variable modification generation](docs/MODIFIED_PEPTIDES_SUPPORT.md) and [named/anonymous definition sets and mass matching](docs/MODIFICATION_DEFINITIONS_SUPPORT.md), caller-owned registries and OBO/crosslink lookup, [charge/pI, hydrophobicity, amino-acid indices and gas basicity](docs/PEPTIDE_PROPERTIES_SUPPORT.md), [terminal/internal spectra, immonium ions, activation presets and compact mass ladders](docs/THEORETICAL_SPECTRA.md) with losses, precursors and annotations |
| Configuration and utilities | [Owned logging streams](docs/LOG_STREAM_SUPPORT.md), [progress reporting and timing](docs/PROGRESS_LOGGER_SUPPORT.md), [Typed parameter values](docs/PARAM_VALUE_SUPPORT.md), [hierarchical defaults, restrictions, updates and CLI parsing](docs/PARAM_SUPPORT.md), [default-parameter lifecycles](docs/DEFAULT_PARAM_HANDLER_SUPPORT.md), [OpenMS INI read/write](docs/PARAMXML_SUPPORT.md), [literal text/CSV helpers](docs/TEXT_CSV_SUPPORT.md) and [list utilities](docs/LIST_UTILS_SUPPORT.md); [filesystem, runtime paths and owned temporary resources](docs/SYSTEM_FILE_SUPPORT.md) |
| Processing | Normalization, threshold/top-N/window filtering, scaling, linear resampling, Gaussian/Savitzky–Golay smoothing, morphological baseline correction, median/iterative-mean noise estimates, HiRes and iterative centroiding, simple and Poisson/KL deisotoping with charge conversion |
| Chromatogram processing | Legacy/corrected peak picking, exact sample boundaries, raw intensity sums, time-weighted integration, baseline estimates, sampled shape metrics and EMG reconstruction |
| Spectrum comparisons | Absolute/ppm alignment, alignment scores, sparse bins, cosine/shared/agreeing scores, precursor, Zhang and Stein/Scott scores |
| RNA chemistry | [Nucleotide records, full pinned registry, TSV/JSON providers and nucleic-acid sequences](docs/RNA_SUPPORT.md), terminal/sulfur linkages, all source fragment formulas and masses, owned custom chemistry and checked slicing |
| Molecular adducts | [Neutral mass/m/z conversion, monomer/dimer notation, electron-aware shifts and formula compatibility](docs/ADDUCT_SUPPORT.md) |
| Decoy sequences | [Whole-protein and peptide reversal, seeded peptide shuffling and deterministic variants](docs/DECOY_GENERATION_SUPPORT.md), with source enzyme/cache conventions and checked work limits |
| Sequence tags | [Residue-mass tags from measured spectra](docs/TAGGER_SUPPORT.md), charge hypotheses, fixed/variable modifications, I/L alternatives and bounded atomic append |
| Spectrum annotation | [Fragment labels, hit peak annotations and match statistics](docs/SPECTRUM_ANNOTATION_SUPPORT.md), source ion-name parsing/display, shared calculation budgets and atomic updates |
| Precursor purity | [Scalar isolation scores, SPS matching, fuzzy scan purity and RT interpolation](docs/PRECURSOR_PURITY_SUPPORT.md), parent-scan lookup and checked experiment scoring |
| Metadata and identifications | Typed values and CV terms, acquisition settings, peptide/protein evidence and scores, protein coverage/modifications, attachments to spectra and feature maps |
| Identification graph | [Owned sequence and provenance records](docs/IDENTIFICATION_GRAPH_SUPPORT.md), stable typed references, score histories, parent/match groups, coverage, atomic registration/merge/copy, [referential cleanup and filtering](docs/IDENTIFICATION_CLEANUP_SUPPORT.md), legacy sequence/evidence conversion and RNase integration |
| Identification analysis | Score categories/switching, HyperScore/Morpheus fragment scores, peptide/protein filters, target/decoy FDR and q-values, picked proteins and probability estimates; peptide-to-protein indexing, basic protein inference/grouping, feature/spectrum conflicts and file-origin splitting |
| Retention-time models | Linear regression with coordinate weights, linear/natural-cubic interpolation, robust LOWESS, inverse refits, deviations, residual windows and atomic application to experiments/feature maps |
| File interchange | [Complete file-type registry and shared native experiment dispatch](docs/FILE_HANDLING_SUPPORT.md), streaming FASTA/MGF, DTA and [MS2/DTA2D](docs/TEXT_PEAK_LIST_SUPPORT.md), optional bounded mzML with [scientific loading filters and canonical arrays](docs/MZML_LOAD_OPTIONS_SUPPORT.md), [reusable parameter groups](docs/MZML_PARAM_GROUPS_SUPPORT.md), precursor acquisition and annotation arrays, [featureXML](docs/FEATUREXML_SUPPORT.md), [consensusXML with protein quantities](docs/CONSENSUSXML_SUPPORT.md), and [idXML with native file paths and identification dispatch](docs/IDENTIFICATION_PATH_SUPPORT.md) with portable custom chemistry; gzip/bzip2 file transport |
| RNA processing | [Fourteen RNases, owned enzyme records and atomic graph registration](docs/RNASE_SUPPORT.md), [fixed/variable RNA modification generation](docs/RNA_MODIFICATION_SUPPORT.md), and [single/multiple annotated RNA spectra](docs/RNA_SPECTRUM_SUPPORT.md), with source mass/charge and terminal conventions |
| Examples | Read/filter/normalize/write spectra; stream FASTA and calculate tryptic peptide masses; inspect modified-peptide masses, isotopes and theoretical fragments; identify a synthetic modified peptide and compute protein coverage; pick and integrate chromatograms; reconstruct a cropped EMG peak; estimate profile noise, centroid and retain local peaks; enumerate modified digest products; stream natural and enriched isotope configurations; calculate peptide physicochemical properties; annotate measured fragments and inspect matching statistics; generate decoy FASTA; extract and match sequence tags; inspect charged RNA formulas and isotope patterns; digest RNA, enumerate variants and generate annotated fragments |

See the [coverage and differences](docs/PORTING_STATUS.md), [chemistry support](docs/CHEMISTRY_SUPPORT.md), and [mzML subset](docs/MZML_SUPPORT.md) before using this as a substitute for a C++ workflow. [Feature containers](docs/FEATURE_SUPPORT.md), [HiRes peak picking](docs/PEAK_PICKING_SUPPORT.md), [theoretical spectra](docs/THEORETICAL_SPECTRA.md), [simple deisotoping](docs/DEISOTOPING_SUPPORT.md), and [Poisson/KL deisotoping](docs/AVERAGINE_DEISOTOPING_SUPPORT.md) each have a defined supported subset.

[Typed metadata](docs/METADATA_SUPPORT.md), [identification records](docs/IDENTIFICATION_SUPPORT.md), and [retention-time models](docs/TRANSFORMATIONS_SUPPORT.md) document the data and alignment capabilities. [Score handling and fragment scores](docs/SCORING_SUPPORT.md), [identification filtering](docs/ID_FILTER_SUPPORT.md), [FDR calculations](docs/FDR_SUPPORT.md), and [idXML interchange](docs/IDXML_SUPPORT.md) describe supported identification workflows and source conventions. [Peptide indexing](docs/PEPTIDE_INDEXING_SUPPORT.md) links identified sequences to FASTA proteins; [basic protein inference](docs/PROTEIN_INFERENCE_SUPPORT.md) aggregates evidence and resolves groups. [Conflict resolution](docs/ID_CONFLICT_SUPPORT.md) and [file-origin splitting](docs/ID_RIPPER_SUPPORT.md) handle competing and merged identifications.

The [chromatogram picker](docs/CHROMATOGRAM_PICKING_SUPPORT.md) preserves source smoothing, seed and boundary conventions. The [peak integrator](docs/PEAK_INTEGRATION_SUPPORT.md) supports spectra and chromatograms, including the source's nonuniform Simpson averaging and sampled shape metrics. Optional [EMG fitting](docs/EMG_SUPPORT.md) reconstructs cropped peaks with the source iRprop+ optimizer and exposes typed fit diagnostics.

The [iterative picker](docs/ITERATIVE_PICKING_SUPPORT.md) refines HiRes seeds and reports exact input regions alongside the source’s rounded centroid and boundary arrays. [Window filtering](docs/WINDOW_MOWER_SUPPORT.md) supports sliding and jumping windows; the [iterative mean noise estimator](docs/MEAN_NOISE_SUPPORT.md) preserves the source’s three-pass clipping conventions.

IsoSpec layered traversal; ProForma; arbitrary RNA enzyme regexes and XML import; other peak-picker families; feature finding/grouping; database search; probabilistic protein-inference engines; broader quantification and OpenSWATH workflows; vendor RAW formats; and Arrow/Parquet remain unimplemented. Identification-graph groups, referential cleanup and the bounded legacy sequence/evidence conversion bridge are implemented; graph persistence and the remaining converter APIs are outstanding.

The [mobility containers](docs/MOBILOGRAM_SUPPORT.md) provide checked mobilogram
search, sorting, selection and summaries. Generic `DataArray` values now retain
metadata and shared processing descriptions; existing struct literals need
`..Default::default()` or `DataArray::new`. Unsupported XML projections reject
these descriptions before output. [Integer mass weights](docs/IMS_WEIGHTS_SUPPORT.md)
and [mzML file operations](docs/MZML_PATH_SUPPORT.md) expose the next native SDK
utilities, including direct filtered loading and atomic compressed output.

## Use locally

Rust 1.85 or newer is required. Add this local crate to a consuming project:

```toml
[dependencies]
openms = { path = "/absolute/path/to/OpenMS4-R" }
```

The `mzml`, `idxml`, `paramxml`, `featurexml`, `consensusxml` and `rna-json` features are enabled by default. The XML features use Rust XML parsing. The independent `numpress` feature adds the base64/zlib wrapper and is enabled by mzML; raw Numpress remains available without features. The `file-compression` feature supplies gzip and bzip2 using Rust backends and is enabled by mzML, parameter XML, featureXML and consensusXML. `rna-json` adds serde_json for caller-supplied MODOMICS JSON. The embedded RNA registry and TSV reader remain available without that feature. Runtime reporting uses `chrono` for local timestamps and `cpu-time` for process CPU timing on Unix/Windows. The scientific core uses the small pure Rust `libm` library for EMG's complementary error function. Disable default features to use the scientific core, identification analysis and text readers without the XML, JSON and compression dependencies:

```toml
openms = { path = "/absolute/path/to/OpenMS4-R", default-features = false }
```

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

Current SDK target: `okohlbacher/OpenMS4-core` at `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. Historical implementation and fixture provenance retains `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`, verified against GitHub on 2026-09-10. See the [SDK update](docs/CORE_SDK_UPDATE.md), [source provenance](SOURCE_PROVENANCE.json), the [historical source inventory](docs/source-inventory.json), and [fixture provenance](tests/data/README.md).

Implementation: BSD-3-Clause, with original OpenMS copyright and attribution in [LICENSE](LICENSE) and [AUTHORS](AUTHORS). Bundled UniMod-derived modification data: Design Science License, with complete source data and notices included. See [component licenses](LICENSES.md) and [data provenance](resources/modifications/README.md). The [RNA data notices](resources/rna/README.md) separately record the unresolved MODOMICS redistribution terms. This local package has not been published to crates.io.
