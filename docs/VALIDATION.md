# Validation of the ongoing Rust port

## SDK update and extraction operations (2026-09-11)

The port now targets Core SDK 4.0.0 at `54a232f`. This increment completes the
reviewed native FASTA lifecycle and indexed mzML offset decoder, adds experiment
aggregation and XIC extraction with mzML Product interchange, and implements
peak display/hash traits. The original provenance is retained for carried-forward
fixtures; changed source files receive explicit review records.
[Recorded results](extraction-validation.json) contain the commands and log hashes.

| Check | Result |
| --- | --- |
| Rust 1.98, all features and all targets | 1,356 tests passed; examples compiled |
| Documentation examples | Three passed |
| Strict Clippy, all targets | Passed on Rust 1.98 with all features and Rust 1.85 without defaults |
| Rustdoc, optimized library, Rustfmt | Passed; documentation warnings denied |
| Indexed mzML alone, Rust 1.85 | Nine tests passed |
| Source/fixture audit | 1,834 distinct source/registration/reference files verified |
| Source-review and coverage regressions | Three source-review tests and two ledger tests passed |

All added components also passed their targeted current/minimum compiler tests
before integration. These overlapping runs are not added to the full-suite count.
Independent review closed XML lexical validation, duplicate-attribute work limits,
FASTA byte-wrapping and mzML Product metadata issues. Existing EMG regression tests
cover the private fitted-trace storage adjustment required by the larger Product
representation.

The preceding published commit `8fb47cb` passed Linux, macOS, Windows and Rust 1.85
CI, including 1,310 tests in the Linux all-target suite. Its quality job found a
new Rust 1.98 test-expression lint; the corrected expression and entire combined
crate now pass strict local Rust 1.98 checks.

The completion inventory still records 766 headers requiring implementation or
review. Twenty headers have a reviewed complete implementation or native equivalent;
this does not establish full dependency or TOPP workflow parity. No C++ executable
was built or run and no TOPP workflow is certified.

## Logging, progress and remaining map helpers (2026-09-10)

Owned logging, replaceable progress reporting, public feature/consensus
modification collection and compressed INI loading are integrated. Real local
timestamps and process CPU timing use pinned safe Rust adapters. The optimized
all-feature library builds. [Recorded results](runtime-foundations-validation.json)
include each command, outcome and log hash.

| Check | Result |
| --- | --- |
| Combined runtime, XML, definition, filesystem and unit suites, current Rust | 206 tests passed |
| Selected runtime, definition, filesystem and unit suites, Rust 1.85 without defaults | 128 tests passed |
| INI feature alone plus filesystem integration, Rust 1.85 | 32 tests passed |
| Final unit and CSV regression checks | 73 tests passed |
| Strict Clippy, all targets | Passed on current Rust with all features and Rust 1.85 without defaults |
| Doctests, Rustdoc, optimized library, Rustfmt | Passed; three doctests and documentation warnings denied |
| Source/fixture audit | 1,827 distinct current source/registration/reference files; four unchanged logging fixtures |
| Windows-style Git checkout filters | All four tested source fixtures retain exact bytes |
| Completion ledger and regressions | All 786 headers accounted for; two regression tests passed |

These selected runs overlap and supplement the complete 1,276-test suite in the
preceding increment. They are not a new complete all-target run. Reviews closed
an empty-name work-accounting gap and a logging flush-failure path that allowed
later output. Logging and progress also passed their isolated staged tests on
both compilers before integration.

CI for the preceding commit passed macOS and Rust 1.85. Windows exposed CSV
fixture checkout conversion, and the quality job exposed a Clippy test-expression
diagnostic. This increment protects all scientific fixture bytes and uses a
portable mutable byte array; local regression checks passed. Updated CI results
are tracked separately from these local results.

The completion gate still fails for 770 headers requiring implementation or
explicit review. Public modification collection and the native progress API now
have complete reviewed mappings. ZIP input, some platform logging/filesystem
behavior, full format/API coverage and executed C++/TOPP parity remain open.


## Map interchange and runtime resources (2026-09-10)

The target remains Core SDK 4.0.0 at `6bfc0e4`. Native featureXML and consensusXML,
portable modification definitions, typed feature metadata, filesystem helpers,
and shared gzip/bzip2 transport are integrated. The complete current-Rust suite
passes **1,276 tests**. [Machine-readable results](map-interchange-validation.json)
record the commands, outcomes and log hashes.

| Check | Result |
| --- | --- |
| Current Rust, all features and all targets | 1,276 tests passed; examples compiled |
| Rust 1.85, consensusXML alone plus unit/dispatch suites | 86 tests passed |
| Rust 1.85, featureXML alone plus unit/dispatch suites | 85 tests passed |
| Rust 1.85, idXML alone plus unit/custom-definition suites | 100 tests passed |
| Rust 1.85, no default features, selected unit/filesystem/definition suites | 95 tests passed |
| Strict Clippy, all targets | Passed on current Rust with all features and Rust 1.85 without default features |
| Final strict idXML-only Clippy, Rust 1.85 | Passed after narrowing a map-only helper's feature gate |
| Doctests, Rustdoc, optimized library, Rustfmt | Passed; three doctests, documentation warnings denied |
| Source/fixture audit | 1,825 distinct current source/registration/reference files; 11 batch fixtures verified |
| Coverage inventory and its regressions | All 786 headers accounted for; two regression tests passed |

Focused test counts overlap and must not be added to the full-suite count. Both
map dialects were tested against original fixtures; consensusXML output also
passed the original XSD through the locally available `xmllint`. Metadata,
protein references, portable chemistry, compressed paths and atomic failures
have independent regression coverage. Reviews caught and closed reference
expansion and validation-order resource-limit gaps.

The completion review now records ZIP input as missing from both map XML
loaders, and compressed input as missing from the INI loader. The INI entry was
therefore corrected from complete to partial. Map-based public modification
collection helpers, platform-specific filesystem behavior, inherited XML schema
validation, and broad SDK parity remain open. The completion gate still fails
for 772 headers requiring implementation or explicit review; it is not a count
of wholly absent Rust classes. No C++ executable was built or run, and no TOPP
workflow is yet certified as port-ready.


## Historical SDK foundations increment (2026-09-10)

The target remains Core SDK 4.0.0 at `6bfc0e4`. This increment adds identification
cleanup, reusable mzML parameter groups, MS2/DTA2D, file dispatch, text/list
utilities and the parameter/INI lifecycle used by TOPP tools. The final optimized
all-feature library builds successfully. Results are recorded in
[sdk-completion-validation.json](sdk-completion-validation.json).

| Check | Result |
| --- | --- |
| Full all-feature/all-target suite before the additive configuration layer | 1,128 tests passed |
| Final unit, parameter, INI, text and list suites, current Rust | 137 tests passed |
| Same selected suites, Rust 1.85 without default features | 122 tests passed |
| INI feature alone, Rust 1.85 | 14 tests passed |
| Final file dispatch and unit regression checks | 67 tests passed |
| Selected file, graph and unit suites, Rust 1.85 without default features | 99 tests passed |
| Strict Clippy, all targets | Passed with all features on current Rust and no default features on Rust 1.85 |
| Doctests, Rustdoc, optimized library, Rustfmt | Passed; three doctests, documentation warnings denied |
| Pinned source/fixture audit | 1,813 distinct source, registration and reference files verified |
| Completion inventory and its two regression tests | Passed; all 786 registered SDK headers accounted for |

These overlapping test counts are separate runs, not an additive total. The
full 1,128-test run preceded the additive parameter/text layer; focused tests
cover that layer and final whole-crate lint, documentation and release checks
cover the combined implementation. An initial inventory check found a stale
generated ledger during development; regeneration and the final check passed.

The completion gate deliberately still fails: 772 headers need implementation
or further review. Fourteen headers have a reviewed complete implementation or
native equivalent. An unmapped header can have Rust functionality that still
needs explicit review; this is not a count of wholly missing classes. The
[completion ledger](CORE_SDK_COMPLETION.md) is an inventory, not a completion
percentage. Source fixtures, schema validation and native library workflows do
not establish executed C++ differential parity or certify a complete TOPP tool.

## Historical reduced-SDK graph increment (2026-09-10)

The port now targets Core SDK 4.0.0 revision `6bfc0e4`. The optimized all-feature
library builds successfully. The complete all-feature/all-target suite passed
1,034 tests on Rust 1.96 and Rust 1.85 before this additive observation/match
layer; the no-default current run passed 905 tests before the layer was added.
The additive graph/formula checks now pass 54 private tests, 12 graph-operation
tests, 12 graph-record tests, 9 source-reference tests, 8 observation-match
tests, and 8 peptide-formula tests on the current compiler. The corresponding
focused suites pass on Rust 1.85 with no default features. Strict focused Clippy,
Rustfmt, source-target inventory checks and the 442-hash/438-link audit pass.

This increment adds observation and compound records, graph adducts, typed
molecules, observation matches, best-match queries, translated ownership and
typed peptide fragment formulas. At that stage, graph groups, cleanup/persistence
and the legacy converter were outside the port. Later increments add groups,
cleanup and the sequence/evidence conversion bridge; graph persistence and full
conversion remain outstanding.

The machine-readable results for this increment are in
[graph-validation-results.json](graph-validation-results.json).

Validated locally on **2026-09-10, macOS ARM64**. The complete all-feature suite
passes **993 tests on Rust 1.96.0: 947 integration tests across 110 suites and
46 unit tests**. The declared minimum Rust **1.85.0** passes all 129 selected
RNA and unit tests with all features. Both compilers pass 124 selected RNA and
unit tests without default features, plus three all-feature doctests each.
Full-library parity remains in progress.

The [machine-readable results](validation-results.json) retain the earlier
complete 938-test runs on both compilers and identify the current checks with
`rna_processing_` names. The new complete 993-test run covers every implemented
chemistry, processing, kernel and format suite on Rust 1.96. The minimum-compiler
and no-default-feature runs in this increment are focused checks, not complete
993-test runs.

## Historical RNA processing checks

| Check | Result |
| --- | --- |
| All features and all targets, Rust 1.96 | 993 tests passed; all fifteen examples compiled |
| All features, Rust 1.85, all RNA suites and unit tests | 129 tests passed |
| No default features, both compilers, all RNA suites and unit tests | 124 tests passed on each compiler |
| Doctests, all features | Three passed on each compiler |
| Clippy, all features/targets, Rust 1.96, warnings denied | Passed for the complete crate, tests and fifteen examples |
| Clippy, no default features/all targets, Rust 1.85, warnings denied | Passed |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| RNA processing example | Executed; RNase_T1 digestion, uridine variants and annotated negative b/y fragments |
| RNA registry regeneration | Original JSON/custom TSV reproduce all 378 embedded records and 375 distinct codes |
| RNA enzyme regeneration | Original XML reproduces all fourteen embedded enzymes and 84 field values |
| Package file list | 428 files; no package build or publication; RNA data release terms tracked separately |
| Current provenance and documentation | 421 source/fixture hashes and 413 local documentation links verified; RNA processing source lines, registry fields and ion bits checked alongside all retained scientific fixtures |

This increment adds fifty integration tests in five suites: twelve RNase,
fourteen modification-generation, twelve spectrum-generation, nine independent
source-reference and three workflow tests. Five new private tests verify the
shared formula work and allocation allowance, including the spectrum generator's
formula-only precursor path. The selected matrix includes all earlier RNA
record/provider, sequence, source-reference and workflow tests to check the
shared sequence helper refactor on both compilers.

Independent references cover fourteen enzyme records, 6,048 pattern/code cases,
twelve original digestion cases with 38 products, all four modification counts
(7, 6, 27 and 432), and all 126 fragment masses actually compared in the source.
The workflows connect modified cleavage, positional products and variable
modifications to annotated spectra and both mzML compression modes. Tests also
cover retained metadata, source charge/sulfur conventions, resource exhaustion
and atomic failures. No expected scientific value was changed to fit an
implementation result. See [RNase support](RNASE_SUPPORT.md),
[RNA modification generation](RNA_MODIFICATION_SUPPORT.md),
[RNA spectra](RNA_SPECTRUM_SUPPORT.md) and the
[independent review](RNA_PROCESSING_REFERENCE_REVIEW.md).

The checks below are retained historical evidence, with their original scope.

## Retained complete RNA foundation checks

| Check | Result |
| --- | --- |
| All features and all targets, Rust 1.96 and 1.85 | 938 tests passed on each compiler: 897 integration tests across 105 suites and 41 unit tests |
| No default features and all targets, Rust 1.96 | 863 tests passed |
| No default features, Rust 1.85, original RNA suites and unit tests | 70 tests passed |
| JSON-only feature, Rust 1.96, original RNA suites and unit tests | 73 tests passed |
| Doctests, all features | Three passed on each compiler |
| Clippy, all features/targets, warnings denied | Passed for the complete crate, tests and fourteen examples |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| RNA example | Executed; charged formula/mass/m/z, sulfur-aware suffix and coarse isotope probabilities |
| RNA registry regeneration | All 378 embedded records and 375 distinct codes reproduced |
| Package file list | 403 files at that stage; no package build or publication |
| Provenance and documentation at that stage | 359 source/fixture hashes and 382 local documentation links verified |

The foundation added thirteen record/provider tests, nine sequence tests,
eight independent source-reference tests and three workflows. These cover all
fifteen original formula assertions, eighteen mono/five average mass assertions,
ten positive slice examples and every registry entry. The JSON dependency
exposed two test-only empty-array type-inference ambiguities, resolved with
equivalent `is_empty()` assertions. A narrow precision-lint allowance preserves
the original carbon-13 mass literal in one workflow; the affected workflows
were rerun on both compilers. See [RNA support](RNA_SUPPORT.md) and the
[foundation reference review](RNA_REFERENCE_REVIEW.md).

## Earlier Tagger-specific checks

| Check | Result |
| --- | --- |
| All features, Rust 1.96 and 1.85 | 58 tests passed on each compiler: all 41 unit tests plus 17 Tagger integration tests |
| No default features, Rust 1.96 and 1.85 | 57 tests passed on each compiler; the mzML workflow excluded |
| Clippy, all features/targets, warnings denied | Passed for the final crate, tests and thirteen examples |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| Tag extraction example | Executed; five measured peaks produce the exact sorted tags EP, EPT, PE, PEP, PEPT and PT |
| Package file list | 380 files, including Tagger, its example, source fixtures and documentation; no package build or publication |
| Provenance and documentation at that stage | 293 source/fixture hashes and 353 local documentation links verified; six Tagger count rows and 120 membership rows verified at their original source lines alongside retained scientific fixtures |

The independent suite reproduces all six original Tagger counts and all 120
membership assertions using the source's theoretical-spectrum settings. The
source input sizes, 357 and 180 peaks, are checked first. Derived tests isolate
strict mass-window boundaries, ties, exact collisions, modified I/L behavior,
signed/unsorted coordinates and zero/inverted settings. Twelve private mass-table
tests cover the 19 base residues, two-stage resolution, deterministic provider
order, empty short IDs and source free-residue arithmetic.

Three private traversal/sort tests include a 9,999-residue path without recursion
and failures after partial traversal that preserve old string allocations. The
workflows connect measured tags to exact target/decoy substring matching,
modified digestion fragments and both mzML compression modes. See
[Tagger support](TAGGER_SUPPORT.md) and its [independent review](TAGGER_REFERENCE_REVIEW.md).

## Retained decoy and adduct checks before Tagger

| Check | Result |
| --- | --- |
| All features, Rust 1.96 and 1.85 | 58 tests passed on each compiler: all 26 unit tests plus 32 new integration tests |
| No default features, Rust 1.96 and 1.85 | 56 tests passed on each compiler; both XML workflows excluded |
| idXML without mzML, Rust 1.96 | 57 tests passed; decoy idXML workflow included |
| Clippy, all features/targets, warnings denied | Passed for the final crate, tests and twelve examples |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| Decoy FASTA example | Executed; two targets and two deterministic variants per target, six records total |
| Package file list | 369 files, including full Boost license notices and the reproducible RNG oracle; no package build or publication |
| Provenance and documentation at that stage | 277 source/fixture hashes and 336 local documentation links verified; new source/Boost hashes and RNG oracle checked alongside all retained scientific fixtures |

Decoy reference tests replay all thirteen unique original sequence cases,
including their shared random/cache history. Nine private RNG tests cover four
complete cycles for three seeds, a published standard value, rejection draws,
shuffle order and state handling. Three additional private tests cover identity
and transactional failures. The packaged integer oracle regenerates independent
checkpoints, checksums and permutations without running C++.

Adduct tests preserve the original scalar examples, then distinguish their
rounded tolerance from exact atomic/electron arithmetic. They check parser
edge cases, signed formulas, finite negative masses and native limits. The
workflows compare sodium-adduct isotope masses with complete ion composition,
preserve charge/metadata through mzML, and connect decoy FASTA to indexing,
synthetic FDR values and idXML. See [decoy support](DECOY_GENERATION_SUPPORT.md),
[adduct support](ADDUCT_SUPPORT.md) and their linked reference reviews.

## Retained full baseline before decoy and adduct additions

| Check | Result |
| --- | --- |
| All features and all targets, Rust 1.96 | 815 integration tests and 14 unit tests passed; eleven examples compiled |
| Documentation examples/API encapsulation | 3 doctests passed |
| No default features | 747 integration tests, 14 unit tests and 3 doctests passed; mzML/idXML excluded |
| All features and all targets, Rust 1.85 | Same 815 integration tests and 14 unit tests passed |
| Rust 1.85 documentation tests | 3 passed |
| idXML without mzML | 772 integration tests, 14 unit tests and 3 doctests passed; no base64/zlib feature dependency |
| Clippy, all features/targets, warnings denied | Passed |
| Rustfmt check | Passed |
| Rustdoc, warnings denied | Passed |
| Spectrum processing example | 121 input peaks → 14 retained peaks; normalized TIC 1.000000 |
| FASTA digestion example | Executed; positional peptide/mass/m/z TSV produced |
| Modified peptide analysis example | Executed; masses, five-bin isotope envelope and annotated b/y fragments produced; invalid input independently checked |
| Synthetic identification example | Executed; selected AC(Carbamidomethyl)DMK with seven fragment matches and 26.32% protein coverage; precursor/evidence/metadata checked |
| Chromatogram integration example | Executed; two synthetic peaks, original sample boundaries, raw intensity sums, time-weighted areas and baseline estimates produced |
| EMG fitting example | Executed; 36 cropped samples → 80 fitted samples, 700 iterations; sampled area 2,121.476414 → 2,501.460284 versus known complete synthetic area 2,506.628275 |
| Profile processing example | Executed; 401 profile samples → four iterative centroids → two retained peaks with aligned integration and width annotations |
| Modified peptide enumeration example | Executed; two digest products → four alkylated/oxidized variants plus one unchanged product; masses and doubly charged m/z values printed |
| Isotope streaming/enrichment example | Executed; first five natural glucose configurations and two enriched-carbon configurations above absolute threshold 0.1, with raw/log probabilities |
| Peptide-property example | Executed; modified and unmodified peptides produce finite charge, pI, GRAVY and gas basicity at both 500 K and 100 K |
| Spectrum annotation example | Executed; eleven source IFSQVGK measurements annotated with source b/y labels, exact widened intensity 12.100000262260437, longest y-series six and finite matching statistics |
| mzML input processing example | Executed using the independent mixed-precision fixture |
| Independent XML XSD checks | Populated/empty mzML, native idXML, auxiliary-array and populated-precursor mzML writer output validated with xmllint and their pinned schemas |
| Generated modification data | All 3,035 specificity rows regenerate exactly; original XML matches pinned archive byte for byte |
| Generated enzyme data | All 33 source records and 602,316 independent Python regex contexts regenerate exactly; bundled XML is byte-identical to the pinned source |
| Package file list | 351 files; spectrum annotation/ion naming, shared-work regressions, peptide-property tables/numerical fixtures, original data, schemas, isotope streaming/custom inputs and fine isotope/precursor purity integration, theoretical extensions, XLMOD/OBO and prior scientific fixtures, examples and licenses included; local reference snapshot and target directory excluded |
| Source/fixture provenance and documentation | 252 source/fixture hashes and 307 local documentation links verified; 41 spectrum-annotation rows and source/derived bits, 633 peptide-property rows and their binary64 values, isotope streaming/custom inputs, fine isotope/precursor purity and theoretical-spectrum sources, all 148 XLMOD mass literals/IEEE values, byte-identical ontology, prior scientific fixtures and libm archive/lock checksum checked; earlier chromatogram/schema checks retained |

A constructor-level source review during validation corrected `precursor_in_ms2`
to integer 0/1 and the spectrum ppm flag to the string projection `"0"`/`"1"`.
The 21 affected annotation/reference/workflow tests were rerun on both compilers
with all features; the optional-feature matrix, lint, documentation and examples
completed after the correction.

The unchanged modification/enzyme generators and data retain their successful
checks from earlier on the same date; they were not rerun for these API edits.

The package-list check does not publish or build a package. Cargo reported only
optional missing documentation/homepage/repository metadata for this local port.

## Test coverage

| Suite | Tests | What it establishes |
| --- | ---: | --- |
| RNase records, registry and digestion | 12 | Fourteen source enzymes, complete record identity, replacement/alias semantics, modified-code cleavage, end gains, product order/coordinates and atomic limits |
| RNA modification generation | 14 | Fixed/variable selection, terminal-site and maximum-one source quirks, alternative identity/order, exact preflight counts/bytes and atomic failures |
| RNA spectrum generation | 12 | Nine ion series, single/multiple charge branches, source mass/sulfur/precursor conventions, aligned annotations and shared resource limits |
| Independent RNA processing references | 9 | All fourteen enzymes, 6,048 predicates, 38 source digestion products, modification counts/identities and 126 compared source ion literals |
| RNA processing workflows | 3 | Modified cleavage, positional digests and variant mass shifts, annotated spectra and both mzML compression modes |
| RNA formula and spectrum budget internals | 5 | Shared work/allocation depletion, preserved mass bits and input state, empty-sequence behavior and cumulative formula-only precursors |
| RNA records and providers | 13 | Complete fields, identity, ordered duplicate registry, ambiguity resolution, JSON/TSV differences, diagnostics and checked limits |
| RNA sequences | 9 | Source parser/display, custom ownership, terminal and sulfur slicing, all fragment variants, natural-H/electron mass convention and atomic failures |
| Independent RNA references | 8 | Fifteen original formulas, 23 mono/average mass references, ten slices, 378-record registry projection and source branch/identity cases |
| RNA workflows | 3 | Independent charged formulas, custom carbon-13 registry lifetime, coarse/fine isotopes and both mzML compression modes |
| Sequence tags | 8 | Defaults/setter, exact append ownership/order, finite source option edges, ignored spectrum fields, registry lifetime and atomic resource failures |
| Independent Tagger references | 6 | Six exact source counts, 120 membership assertions, source peak counts/rounded traces, strict bounds, nearest ties, collisions, I/L, signed/unsorted and no-op behavior |
| Tagger mass-table internals | 12 | Nineteen free-residue formulas, fixed/variable order, exact keys, two-pass provider resolution, terminal/wildcard/anonymous records, formula/mass precedence and construction limits |
| Tagger traversal/sort internals | 3 | A 9,999-edge heap traversal, failure after emitted paths, string-byte sort work and allocation-preserving atomic errors, strict lookup/ties |
| Tagger workflows | 3 | Exact substring target/decoy matching, modified digest fragment gaps and both mzML compression modes |
| Decoy generation | 7 | All registered enzymes, reversal/shuffle source rules, composition where applicable, modified/empty inputs, cache/resource boundaries and atomic failures |
| Independent decoy references | 7 | All thirteen literal source cases, ordered RNG/cache history, reseeding, cross-context reuse, zero attempts, short products and unspecific order |
| Decoy RNG internals | 9 | Three seeds over four cycles, standard 10000th word, rejection mapping, shuffle order, clone/reseed and no-draw boundaries |
| Decoy identity/state internals | 3 | Forward/reverse maximum identity, strict ties and late cache/work/allocation rollback |
| Decoy workflows | 3 | FASTA, target/decoy evidence and synthetic FDR, composition through digestion and idXML transport |
| Molecular adducts | 8 | Complete parser/components/getters/equality, electron and n-mer conversions, mono/average shifts, compatibility and checked limits |
| Independent adduct references | 5 | Literal source cases, tighter atomic constants, isotope/charge/whitespace grammar, signed containment and source numerical boundaries |
| Adduct workflows | 2 | Complete sodium-ion composition, unchanged isotope probabilities, charge-dependent spacing and both mzML compression modes |
| Chemistry | 12 | Source-derived mass/formula and digestion values, residue/isotope tables, formula grammar/algebra limits, b/y mass conservation |
| Charge and isoelectric point | 10 | Four source pKa scales, terminal overrides/suppression, parent PTMs, U/O/ambiguous residues, saturated finite pH, exact endpoint/midpoint behavior and resource/convergence guards |
| pI internal work accounting | 1 | Endpoint and midpoint evaluations consume one precharged budget |
| Amino-acid indices and gas basicity | 7 | Source accessions/indicators, ordinary split arithmetic, independent low-temperature values, tied maxima, empty/high-temperature identity, formula-free annotations and input limits |
| Hydrophobicity profiles | 7 | Seven scales, literal GRAVY/window/moment values, sliding order, all annotation types, window/angle semantics and preallocated work limits |
| Independent property references | 7 | 503 source constants including 42 sentinels, isolated pKas, 30 source scalar assertions, 420 singleton/pair GB expressions, 100 Decimal-derived extreme-temperature cases and independent moments |
| Peptide-property workflows | 4 | Digestion retains only original terminal caps, mass-only chemistry needs no formula, source generator attachment distinctions and typed property metadata through idXML |
| Ion naming | 9 | Every source charge/ordinal case, first-line and caret priority, field/overflow fallbacks, about 18,000 round trips, Unicode text, output limits and allocation-free parsing |
| Spectrum annotation | 10 | Three source operations, array replacement, ppm last-match/duplicate branches, no-op flags, final sorting, safe small-list statistics, finite source special cases, scoped validation and atomic failures |
| Independent annotation references | 7 | Original measured bits and 17 literal assertions, independently derived fragment errors, exact binary top-N padding/sample variance/ratios, source label grammars and precursor absolute tolerance in ppm mode |
| Annotation workflows | 4 | Modified digestion against independent fragment formulas, distinct chemistry with identical sequence text, peak annotations/statistics through idXML and aligned arrays/acquisition through both mzML compression modes |
| Shared generation/alignment work | 3 | Actual alignment initialization/cells and shared residue/loss/fine/coarse allowances survive successive calls; standalone calls start fresh budgets |
| Annotation grammar internals | 1 | Complete source label regex semantics, including commas, losses, ordinal boundaries and terminal-series distinctions |
| Sequence chemistry and numeric tags | 14 | Source integer/decimal registry lookup, absolute/internal/H/OH masses, unresolved residue representation, formula/mass availability, owned tags and shared immutable spelling through independently owned slices, stable attachment, slicing, atomic setters and numeric parser limits |
| Sequence identification operations | 8 | B/Z/X non-mass filters, known-mass precursor filtering, exact modification IDs, owned protein observations, distinct resolver keys and atomic failures |
| Sequence chemistry workflow | 8 | Ambiguous digestion/indexing, numeric-tag idXML round trips, independent fragment/loss shifts, formula-dependent append rejection and extreme-mass cancellation regressions |
| Modifications | 11 | Full pinned registry load, specificity/name/mass lookup, neutral losses, isotope labels, terminal notation, atomic setters, modified digestion/fragments |
| Modified peptide generation | 12 | Fixed/variable source goldens, stable alternative order, terminal overwrite/duplicates, existing anonymous annotations, bounded combinations, atomic errors and custom formula-free mass/absolute/no-op rules |
| Independent modification generation | 9 | Seven fixed and nineteen variable source cases, full weighted-site order, exact typed terminal/fragment-mass distinctions, append limits, empty idXML rejection and independent custom mass branches |
| Modification definitions | 13 | Owned definition identity, fixed/variable set semantics, compatibility, count non-enforcement, absolute/delta matching, inference, conflicting full-ID chemistry and atomic limits |
| Independent definition reference | 5 | Seven literal compatibility cases, partition/merged precedence, negative-delta tolerance endpoints, stored/fallback absolute masses and named/anonymous all-hit inference |
| Modified peptide workflow | 5 | Digestion → variants → formulas and independent fragment shifts → inferred search definitions → idXML; mass-only fragments, isotope rejection and serialization preflight |
| Owned modification registry | 7 | Caller record validation, optional vocabulary IDs, shared handles, bounded atomic OBO appends, formula fields and complete Eq/Ord contracts including signed zero |
| Independent OBO registry reference | 10 | All 148 XLMOD records in exact source order with identities, specificities, synonyms and mass bits; PSI aliases and absent targets, stanza/EOF handling, literal empty versus zero formula and alias-work limits |
| Crosslink lookup | 6 | Separate 56-record registry, source DSS/BS3/EDC mass/site goldens, reactive-side union, terminal conversion, search eligibility and bounded caller-owned extensions |
| Anonymous modification definitions | 6 | Owned exact-spelling annotations, all-hit inference, source empty-short-ID compatibility, full-residue/H/OH mass anchors, unresolved deltas and atomic errors |
| Caller-owned chemistry workflow | 7 | Registry lifetime release and sharing, numeric/name resolution, absolute-formula X rescue, no-change/delta/terminal precedence, non-UniMod mass export, full-ID-only records and exact custom-registry idXML |
| Custom chemical identity | 6 | Same text with distinct formulas/vocabularies remains distinct in protein observations, owned peptide keys, sequence duplicate filtering, rank/spectrum conflicts; atomic late errors |
| Internal theoretical fragments | 7 | Source interval and ten-residue boundaries, ordered annotations, numeric residue and terminal chemistry, loss replacement, isotope independence, intensities, atomic append and custom declaration/storage limits |
| Compact mass-only spectra | 7 | Full-length six-series ladders, source terminal/float operation order, observed masses, independent suffix accumulation, sorted no-ops and combined/precision limits with late-error atomicity |
| Activation presets and immonium ions | 6 | Source mass goldens, all activation enum cases, inferred precursor metadata, unmodified residue eligibility and L-only branch, charge/intensity/isotope independence and aligned atomic append |
| Independent theoretical extension review | 9 | Analytical internal interval/mass/loss oracles, exact immonium constants, source rounded CID/partial ECD tables and independently rounded f32 compact-helper values |
| Scalar purity and SPS matching | 9 | Source isolation totals/residuals, inclusive doubled tolerances, charge/nearest conventions, f32 fragment windows, input immutability and shared work/unused-annotation limits |
| Fuzzy and interpolated purity | 11 | Neutron-spacing source successor lookup, f32 sums/division, strict and half-weight boundaries, ratios above one, zero-window/empty-parent behavior, RT extrapolation/fallbacks and nonfinite-field checks |
| Independent precursor purity review | 13 | Exact decoded source peaks and seven scalar/map goldens, SPS float32 bounds and counts, independent fuzzy arithmetic and parent reference/acquisition-order rules |
| Precursor metadata and interchange workflow | 10 | Shared native precursor ownership, referenced parents, scalar batches, original complete mzML fixture, all activation/mobility quantities, malformed XML/loss guards and populated acquisition XSD validation |
| Native fine isotope enumeration | 8 | Small exhaustive configurations, fixed labels/natural zero abundances, charge adducts, source rounding/coverage, deterministic thresholds and checked atom/output/frontier resource limits |
| Independent fine isotope references | 7 | All 44 source counts, 14 fructose mass/probability rows, 6 bromine configurations, subnormal/full-support tails and the 19,615-state insulin f32 coverage boundary |
| Independent fine theoretical spectra | 7 | Literal 10/5/50/12 source counts and exact mass/intensity tables, independent CHNO oracles, neutral H/charge division, terminal/loss rules and single/shared-budget atomic failures |
| Fine isotope workflow | 3 | Modified peptide formulas and charge adducts, aligned spectrum selection, unknown-composition errors and both compressed/uncompressed mzML round trips with precursor acquisition |
| Fine isotope internals | 2 | Deterministic equal-probability state identity and cumulative work allowance across calls |
| Fine isotope stream and custom inputs | 8 | Materializer equivalence, original-mode thresholds, charge-ignoring raw adapter, owned lifetimes, fused errors, checked conversion, zero-count validation and large-support prefixes |
| Independent stream references | 6 | Original 14-row/2548-state fructose and 10000-state insulin ordered cases, threshold counts, independently reconstructed natural tables, custom f64 precision and direct multinomial products |
| Fine isotope stream boundary review | 10 | Extreme finite weights, log underflow, late mass overflow, getter-consistent cutoff equality and adjacent values, cumulative custom dimensions, deferred expansion and 103041-state full streaming |
| Custom isotope spectrum workflow | 3 | Owned raw-coverage prefixes and continuation, distinct equal-mass configurations through materialization, aligned spectrum selection and both mzML compression modes |
| Digestion registry and specificity | 12 | All 33 enzyme predicates vs independent regex hashes, full/semi/nonspecific order and small-sequence oracle, validity/missed cleavages, modifications and resource limits |
| Typed metadata and acquisition | 17 | Scalar/list/unit invariants, explicit string bridges, CV merging, source settings defaults/unify, drift/isolation conventions, date/checksum validation and atomic errors |
| Identification records | 18 | Stable scores, typed evidence/metadata, charge ranges, grouping, interval coverage vs independent oracle, observed modifications, kernel attachments and file loss rejection |
| Score categories and switching | 9 | All 29 names/six categories, source lookup/backup conventions, restoration, heterogeneous records, reserved labels and atomic map/slice errors |
| Peptide-spectrum scoring | 8 | Source HyperScore/Morpheus goldens, overload-specific precision/charge/boundary rules, ion annotations/ordinals, error means and checked failures |
| Identification filters | 16 | Source cutoff/top-N/dense-rank/tie behavior, modifications, exact/sequence duplicates, signed charge/precursor errors, run-aware references/groups and map atomicity |
| FDR and q-values | 13 | OMSSA1534 and XTandem source thresholds, picked goldens, actual legacy/Basic formulas, score ties/directions, peptide representatives, posterior estimates/ROC, reserved labels and atomic limits |
| idXML | 16 | Exact pinned fixtures, native metadata/run/evidence round trips, independent XSD validation, malformed XML/reserved encodings, byte/element/list guards and preflight/flush errors |
| Identification processing workflow | 1 | Independent XML → score switching/top-hit selection → hand-computable target/decoy q-values → threshold/reference cleanup → coverage → exact native idXML round trip |
| Peptide-to-protein indexing | 17 | 121 matcher/17 decoy/15 enzyme oracle cases, ambiguity/mismatch/I-L rules, original-engine recovery, per-run evidence/protein reconstruction, metadata placement and preflight/atomic limits |
| Basic protein inference | 21 | Source merger scores/groups, best/product/mean, charge/modification representatives, greedy graph ties/negative scores, score restoration, vector/single/consensus conventions and atomic limits |
| Identification conflicts | 11 | Source rank aggregation and spectrum reports, best/matching selection, modified-sequence/charge keys, intensity winners, metadata/subordinate retention and atomic map errors |
| Identification origin partitions | 7 | Three origin formats, source ordering/path reduction/hit union, basename/key collisions, run-specific protein values, skipped records, group omission and checked errors |
| Combined protein workflow | 1 | FASTA → modified/shared peptide evidence → basic inference → distinct protein/PSM q-values → global coverage → native idXML round trip → file-origin partitions |
| Retention-time models | 12 | Source weighted/unweighted fits, interpolation/extrapolation goldens, LOWESS sine/cars/original tables, inverses, diagnostics/windows and resource/atomic errors |
| Container RT transformations | 7 | Source experiment/feature/consensus/peptide values, original metadata, nested hulls/subordinates/handles, optional spectrum IDs, ordering and late-error atomicity |
| Isotope distributions | 13 | Coarse glucose/heavy/bromine goldens, convolution, enrichment, peptide/RNA/DNA averagine, exact Poisson recurrence rounding and overflow fallback, conditional fragment probabilities, limits |
| Theoretical spectra | 12 | All ion-series mass tables, losses, precursors, isotope conventions, modified peptides, aligned atomic append and resource limits |
| Kernel | 19 | Source search/TIC values, ties/ranges, metadata alignment, mutation checks and malformed inputs |
| Feature/consensus geometry | 18 | Scan-envelope hulls, source containment quirks, consensus means/charges/decharge, map identities/ranges, subordinates and checked mutations |
| Processing | 14 | Original 121-peak fixture, normalization, source 0.05 threshold default and adjacent-f32 cutoff, top-N/rank behavior, atomic errors and resampling conservation/boundaries |
| Smoothing | 10 | Gaussian goldens/trapezoids/ppm behavior; Savitzky–Golay goldens, asymmetric edges and polynomial reproduction |
| Baseline morphology | 8 | All ten source operations, boundary windows, signed bottom-hat, width conversion and aligned atomic updates |
| Peak picking | 12 | Orbitrap/FTMS centroid goldens, full historical noise values, natural spline values/derivatives, FWHM, mobility, boundaries and experiment selection |
| Iterative peak picking | 8 | HiRes seed/refinement noise separation, source centroid and integration conventions, width/spacing rules, aligned annotations, selected experiments and atomic limits |
| Internal iterative conventions | 2 | Fourteen independent refinement configurations and deterministic ties; strict-next-sample association, asymmetric recenter search, original seed priority and rounded centroid storage |
| Independent iterative reference | 4 | Complete synthetic result, original Orbitrap/FTMS raw-region sums and centroid rounding, omitted annotations, experiment behavior and checked failures |
| Iterative profile workflow | 2 | Noise estimation → centroiding → local peak filtering, independent raw sums, preserved metadata and exact mzML width/intensity array round trip |
| Window filtering | 10 | Source 56/30 counts, exact triangle indices, 480 independent cases, sliding/jumping boundaries and order, duplicate membership, stable ties and atomic resource failures |
| Iterative mean noise | 7 | All 2,526 historical source outputs, hand-computed fixed-denominator clipping, strict window edge, histogram ceiling, signed/empty inputs, checked legacy percentile and exact work limits |
| Independent mean-noise reference | 3 | Fresh histogram rescans for 192 configurations, signed global statistics, effective bin ceiling and legacy-percentile f32 arithmetic |
| Chromatogram picking | 9 | Legacy/corrected methods, Gauss/SG smoothing, independent seed/boundary noise, raw sums, exact index regions, aligned arrays, overlap and checked failures |
| Internal chromatogram conventions | 2 | Source right-side closest-sample ties/end sentinel and sequential overlap midpoint assignments |
| Peak integration | 12 | Inclusive sampled sum/trapezoid/Simpson, all baseline choices, signed/empty/singleton inputs, shape metrics, duplicates, bounds and numeric/resource failures |
| Independent chromatogram processing | 12 | All 146 source trace samples and literal integration/background/shape/picker goldens, nonuniform/even Simpson behavior, f32 arithmetic, -1 sentinel exclusion and source noise/overlap conventions |
| Chromatogram workflow | 2 | Raw picking → exact full-precision boundaries → intensity sum/time-weighted area/baseline/shape → mzML peak-annotation round trip; f32 boundary rounding cannot change integration |
| EMG fitting | 8 | Source cutoff fit, three model branches, metadata and omitted arrays, zero bounds, best iterations, extrapolation, shared evaluation limits and numerical failures |
| Internal EMG conventions | 3 | All four analytic gradients against central differences, source training collection order and iRprop sign/zero/rollback rules |
| Independent EMG reference | 8 | Seven complete source fits and four full raw loss goldens, true relative parameter tolerances, ordered training losses, scalar tails/subnormals, branch boundaries and inclusive budgets |
| EMG integration | 6 | Both container types, all integration/baseline methods, unchanged supplied shape inputs, cropped and expanded spans, literal zero bounds and preserved inputs on errors |
| EMG workflow | 2 | Independent continuous area/centroid identities and cropped fit → mzML → sampled integration without unaligned parameter annotations |
| Deisotoping | 8 | Unknown/unequal precursor-charge source regressions, charge priority, ladders, intensities, annotations, disjoint membership and error atomicity |
| Poisson/KL deisotoping | 11 | Source defaults and threshold/top-N preprocessing, original indices and aligned arrays, both sharing policies on both algorithms, low-mass noise, bounded work and transactional errors |
| Complete deisotoping fixture | 3 | All 5,407 source input peaks retained; all 104 expected output peaks bit-exact, independent seed/charge mapping, shared membership, exact 103-peak disjoint subset and retained metadata/arrays |
| Independent deisotoping review | 8 | Adjacent-f32 KL thresholds for sizes 2–7, mixed-precision accumulation, longest/highest-charge selection, selected counts/sums, inclusive ppm endpoints, nearest ties and exact precursor arithmetic on both methods |
| Spectrum comparison | 19 | Source alignment/scoring/binning goldens; 1,000 compact-DP/reference-map comparisons; directed ppm and weighted rounding boundaries |
| Text formats | 12 | DTA/FASTA/MGF parsing, streaming/round trips, invalid input, writer preflight/flush errors and malformed-input corpus |
| mzML | 14 | Independent/upstream fixtures, codecs/precision/units, metadata subset, limits, malformed XML/arrays and independent XSD checks |
| Independent format/processing review | 5 | Regressions for nonfinite filtering, empty resampling, full-range charges and metadata key loss |
| Independent mzML review | 3 | Conflicting scientific CV fields, reserved names and forbidden XML characters |
| Independent mzML auxiliary reader | 10 | Both float/signed integer widths, exact integer range, empty ASCII elements/placeholders, malformed/unsupported metadata, variable string decompression and cumulative byte/element/array limits |
| Independent mzML auxiliary writer | 4 | Validation before first output, exact source binary encodings, zero-length distinctions, type/name order and metadata, all four compression/empty-record XSD cases |
| Synthetic identification workflow | 1 | Embedded FASTA → fixed modification → digestion → literal observed fragment matching → peptide evidence/annotations → protein coverage/modifications; precursor consistency and retained metadata |
| Workflows | 4 | FASTA → digestion/fragments → MGF; DTA → processing → MGF; modified peptide → isotope envelopes → deisotoping → alignment; profile → centroid → consensus |

Spectrum annotation and ion naming received separate source, numerical and
integration reviews. The four new integration suites contain 30 tests, with one
additional sequence-ownership regression and four private grammar/work tests.
Original source fixtures contain thirteen measured peaks and seventeen scalar
or label assertions. Eleven independently derived masses/errors distinguish
rounded source references, including its zero-MSE assertion, from actual
nonzero error statistics. Exact binary examples independently check padding,
sample deviation, quartile indices and current ratios. See the
[annotation support](SPECTRUM_ANNOTATION_SUPPORT.md),
[reference review](SPECTRUM_ANNOTATION_REFERENCE_REVIEW.md) and
[fixture provenance](../tests/data/spectrum_annotation_provenance.json).

Native annotations preserve source last-match array behavior, duplicate
matched-only records, parameter no-ops, sort order and raw absolute precursor
tolerance in ppm mode. Small-list quartiles and undefined enabled statistics
have documented safe policies. Workflows verify modified fragments from
independent formulas, complete custom-chemistry identity, and XML transport.

Shared theoretical and alignment counters now persist across candidate hits;
coarse convolution also retains one allowance across theoretical envelopes.
Focused exhaustion tests verify that calls cannot restart these counters.
Anonymous modification strings are shared immutably, with independent slices
and setters retaining their owned lifetime, value and ordering semantics.

Peptide physicochemical utilities received independent source-table and numerical
reviews. Five integration suites add 35 tests, and a private pI test verifies its
shared evaluation budget. All 200 public AAindex values, 182 hydrophobicity cells
and 59 pKa constants are checked. Private gas-basicity tables feed independent
source expressions for every canonical singleton and all 400 ordered pairs.
Thirty original scalar assertions retain their source spelling and tolerance
context. A separate 100-case grid uses 90-digit Decimal arithmetic, including
low temperatures, tied maxima, subnormal products and the largest finite inputs.
See the [property support](PEPTIDE_PROPERTIES_SUPPORT.md),
[reference review](PEPTIDE_PROPERTIES_REFERENCE_REVIEW.md) and
[fixture provenance](../tests/data/peptide_properties_provenance.json).

Ordinary gas basicity preserves source evaluation order. Stable overflow retries,
the underflowed-product limit and the empty-sequence identity are documented
native numerical differences. Terminal annotation tests preserve source pI
suppression without inventing PTM-specific pKas; hydrophobicity and gas basicity
continue to use parent residues. Digestion and identification interchange tests
exercise those rules together. The current complete Rust 1.96 run includes all
110 integration suites and 46 private unit tests. Three doctests pass on both
compilers, alongside the focused RNA processing checks detailed above.

Owning isotope streaming and custom binary64 populations received independent
source and numerical reviews. The source ordered fructose support (2,548 states),
fourteen literal configurations and insulin's 10,000-state prefix pass alongside
threshold and custom natural-table materialization checks. Custom weights retain
extra binary64 precision; hand-derived products preserve separate equal-mass
configurations. Raw formula charge is ignored while the high-level generator's
hydrogen-adduct rules and all prior fine-spectrum references remain unchanged.

Boundary tests distinguish finite logarithms from underflowed raw probabilities,
check late overflow and fused errors, and validate every custom cell even for
zero-count populations. A 10,000-category case succeeds when coverage stops at
the first configuration, then fails with a checked work error when its raw stream
requests more. Raising its threshold avoids expanding an ineligible suffix.
A separate stream drains 103,041 configurations although materialization rejects
that full support. Threshold tests retain a returned raw probability exactly and
reject the next larger binary64 value; materialized logarithmic cutoff behavior
is preserved. Distinct equal-mass peaks and their configuration indices survive
selection and both mzML compression modes.

The materializer boundary assertion checks the analytic logarithmic condition
rather than assuming a particular platform's `ln`/`exp` rounding direction. Both
compilers pass it alongside the raw-value equality assertions in the full matrix.

The existing search is shared without a new dependency. Iterator cutoffs use
returned raw values, with logarithmic fallbacks for relative underflow or an
unrepresentable mode. Source layers, performance hints and untrimmed layer-based
membership remain outside the current surface. Streaming retains a frontier and
visited set, so atom, work, state and memory limits bound its lifetime.

Native fine-isotope enumeration received independent source, mathematical and
theoretical-spectrum integration reviews. All 44 selected class-test count
assertions, the 14-row high-precision fructose table and six bromine configurations
are checked separately from older loose/header illustrations. The 19,615-state
insulin boundary verifies accumulation of stored f32 probabilities; full-support
tests retain configurations whose stored probabilities underflow to zero. A small
exhaustive combinatorial oracle checks masses and probabilities independently of
the native heap algorithm, including separate fixed labels and natural H adducts.

The theoretical integration reproduces the source's 10-, 5-, 50- and 12-peak
fine-spectrum cases and checks retained terminal chemistry, source loss formulas,
isotope annotation alignment and natural-H mass/charge conventions. Loss intensity
products now remain f64 until final storage, avoiding an intermediate f32 rounding
in both fine and coarse envelopes; existing coarse regressions still pass. A
105-prefix test fails specifically at the shared fine-work limit while each
single envelope is valid and total possible output is only 105 peaks. Both this
failure and late atom-limit errors leave appended observations unchanged. Modified
fine spectra preserve their annotations and precursor acquisition data through
both mzML compression modes.

This validates materialized native configurations with explicit deterministic
ties, inclusive threshold equality and checked resource limits. It does not
claim identical IsoSpec layering or untrimmed layer-dependent selections, or
an executed C++ differential comparison. No dependency was added.

Precursor purity received independent numerical, source and acquisition-format
reviews. Exact binary64 coordinates and binary32 intensities from the original
mzML support seven source scalar/map cases; eight SPS cases and independent fuzzy
arithmetic are also checked. Rounded upstream literals are kept separate from
independently calculated results. The full original fixture reproduces the same
scores through the native mzML reader after an explicit in-memory encoding
declaration change; the bundled source bytes remain unchanged. All 19 activation
methods, four mobility quantities, distinct isolation targets and spectrum
references are covered, including populated writer output checked by xmllint.

Review found full spectrum validation repeatedly traversed unused annotations
outside the purity work budget. Purity now validates only consumed numerical
fields; a 10,000-child/10,000-placeholder regression verifies this bound. Scoped
precursor validation remains explicit before fuzzy empty-parent shortcuts.
Malformed container nesting and loss of CV-list ordinary metadata on export are
also rejected. The writer rejects a negative effective isolation target before
output, including selected-m/z fallback; spectrum/chromatogram regression cases
prevent producing a document that the reader would reject. Native source-undefined iterator/arithmetic cases return checked
errors or documented end-candidate fallbacks; finite source overcount and RT
extrapolation are retained. No C++ binary was executed for these comparisons.

Theoretical-spectrum extensions received independent source and integration
reviews. Analytical references cover internal fragment masses, counts, source
start/end conventions and first-residue loss omission. Literal source immonium
masses and rounded activation tables are checked separately; the compact helper
has an independent bit-exact f32 analytical oracle, rather than captured C++
output. Coarse-envelope independence, custom mass tags, retained termini,
annotation alignment and error atomicity are covered across the focused suites.

Review found repeated custom loss declarations could consume unbounded work
before deduplication and distinct losses could fill templates before peak
preflight. Declaration visits now receive a shared preflight estimate, and
unique loss-template entries are charged before formula cloning. Regression
tests verify both limits for ordinary and internal generation while preserving
existing output. Fine isotope support was added subsequently and is covered above.

Owned modification records, OBO loading and crosslink lookup received independent
source and consumer reviews. All 148 projected XLMOD records match exact source
accession/site order, names, synonym sets and f64 mass bits. Synthetic PSI-style
cases verify all-specificity UniMod aliases, absent targets and absolute formula
precedence. The historical PSI-MOD snapshot remains unbundled; its source hash
is verified separately from packaged fixtures.

Review corrected source empty-sequence export, literal empty versus explicit
zero-formula behavior and repeated alias expansion work accounting. A further
consumer review found same-name custom chemistry being merged in protein
observations, sequence duplicate filtering and conflict keys. Complete chemical
value ordering now preserves those distinctions and agrees with equality,
including signed zero. Peptide identity keys retain owned sequence values.
Source-intentional textual keys in other identification algorithms remain.
Caller-owned registry tests verify record release after the last owner is dropped,
absolute-formula restoration of unknown residues, checked vocabulary/mass export,
and exact idXML chemistry validation before output.

Modified-peptide generation and definition sets received independent source and
integration reviews. Literal cases include 26, 71 and 199 variants, supplemented
by a complete small weighted-site oracle. Source maximum-one/general terminal
placement and formula-derived versus declared mass distinctions are preserved.
Review found and corrected inference merging of chemically different records
sharing a full ID; conflicting records now fail atomically in either order.
Custom formula-free records retain declared/absolute monoisotopic mass through
generation and fragment calculations, with explicit unavailable composition.
The source no-change rule and formula precedence are checked separately.
Normal generated variants and search definitions round-trip through idXML;
unrepresentable typed terminal states are rejected before any writer call.

Iterative peak picking received an independent source and numerical review.
The source class test has no numerical assertions; fourteen derived refinement
configurations explicitly cover source control flow and mixed precision. Two
real input profiles independently verify centroid rounding and raw-region sums;
they are not treated as iterative output goldens. Window selection matches
literal upstream counts and 480 independent small cases. Iterative mean noise
matches all 2,526 historical outputs at the source tolerance and fresh histogram
rescans over 192 configurations. The combined profile workflow verifies that
filtering preserves the picker’s aligned integration and width arrays through
mzML. The implemented source quirks and native checked-error policies are
recorded in the corresponding support documents.

EMG fitting received an independent numerical review. All twelve branch
expressions across four gradients match the source arithmetic, with separate
central-difference checks. Seven source fits retain their literal minute/second
coordinates, f32 container intensities and true relative parameter checks;
four source loss assertions distinguish full raw-f64 loss from training loss.
Derived training-order fixtures and erfc tail/subnormal tests independently
check operation ordering and special-function behavior. The strict numerical
error policy is documented, including failures after a previous finite best.
An additional mathematical workflow verifies total EMG area and centroid, then
checks cropped reconstruction and exact fitted-trace mzML interchange.

Independent reviews inspected the numerical modules and their source conventions.
Chromatogram integration and picking preserve literal source reference values,
mixed f32/f64 operation order, sampled boundaries and shape conventions. Reviews
caught an uncharged Gaussian coefficient-table allocation and verified the
corrected work limit. The combined workflow additionally preserves exact original
sample indices when f32 boundary metadata cannot represent their coordinates.

Auxiliary mzML arrays received independent reader and writer reviews. Tests
verify binary bytes, signed integer limits, ASCII terminators and empty elements,
validation before any write/flush, and cumulative decoded-storage bounds. A
review regression rejects units attached to encoding/compression terms, preventing
unrepresented metadata loss. The provenance audit distinguished the mzML schema
fixture's CRLF bytes from the immutable source's LF bytes and verified identical
normalized schema content; both hashes are now recorded.

Poisson/KL deisotoping received separate numerical and data-flow reviews. The
complete real-spectrum fixture revealed two source clusters sharing a heavy
isotope; the native default preserves both, while an explicit disjoint policy
omits only the independently identified overlapping seed. Source binary arrays,
decimal fixture values and all 104 seed/charge mappings were independently
verified. Review also corrected precursor operation order and low-mass noise
handling, and checked mixed f32/f64 KL acceptance boundaries.

Sequence review checks the distinction between known mass and known formula,
annotation ownership, position-preserving round trips, and source numeric lookup
precision. Regressions cover two large-mass cancellation errors: independently
accumulated suffix fragments preserve a small alanine fragment beside a large
tag, and equal/opposite terminal deltas preserve the peptide residue mass.

Protein inference, peptide indexing, conflict resolution, origin splitting and
the combined protein workflow received independent scientific reviews. The
matcher oracle includes 121 raw cases, 17 decoy-inference cases and 15 enzyme
boundaries. Review regressions address undefined-negative greedy scores,
original search-engine recovery, search-parameter metadata placement and
pre-clone sequence limits.

Identification scoring/filtering/FDR and idXML received separate peer reviews.
Regressions reject target/decoy label replacement, aliased analysis-result
indices, malformed declarations/processing instructions and encoded group lists
over their configured limit. Other reviews checked digestion context, typed
metadata, evidence coordinates, transformation inverses and atomic container
edits. The synthetic identification precursor was corrected against an
independent mass calculation.

All 84 elements and 283 isotope data points were compared against the source
declarations. All 2,341 inventory entries were verified against immutable Git
objects at the pinned revision; one CRLF checkout conversion is recorded
separately. Checks establish the documented surface; they do not count fully
ported C++ classes.

## Limits of this evidence

The C++ reference was inspected, **not built or executed**. These are
source-derived reference tests, independently reconstructed numerical
expectations, fixture/schema checks and Rust invariants. They are not runtime
differential tests against OpenMS or proof of full workflow equivalence.

Linux, macOS, Windows and Rust 1.85 CI have passed for the preceding published
commit, as recorded above; newer local changes have separate validation records.
Linux CI installs xmllint; on other hosts the schema tests explicitly report when
that optional executable is unavailable. It was available and executed locally.
Rust 1.85 compatibility was also exercised locally.

There are no runtime performance/memory benchmarks, external application
interoperability certification or full PSI controlled-vocabulary validation.
The capability documents describe current exclusions and deliberate source
corrections. No C++ source, contrib tree or vendored dependency was modified.
