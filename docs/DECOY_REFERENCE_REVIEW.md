# Independent decoy-generation references

The reference suite checks the complete native DecoyGenerator surface against the pinned OpenMS implementation and its literal class-test outputs. The [fixture manifest](../tests/data/decoy_provenance.json) records nine pinned source hashes, the original assertion locations, and four separately identified supplemental Boost 1.90 header hashes. No C++ executable or native production algorithm generated the expected strings.

## Literal evidence

[DecoyGenerator_test.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/DecoyGenerator_test.cpp) contains 15 relevant string assertions representing 13 unique operation/input/variant cases. All 13 are preserved in [decoy_source_goldens.tsv](../tests/data/decoy_source_goldens.tsv). Duplicate variant-zero assertions retain both line numbers.

The three `shuffle_peptides` cases must run on one object seeded with 4711, in source order. Their expected outputs are `DIESETEPTTP`, `ETTSRTPEPRIED`, and `SPERTETTPRIED`. The second call stores the final `IDE` product's decoy, reused by the third call. Resetting the object for each row would compare a different stateful computation.

The outer `shuffle` method creates a fresh generator seeded with `4711 + variant` for every sufficiently long digested product. Its seven unique fixture cases therefore use the existing receiver only to verify that receiver state has no effect. Its literal top-down sequence contains O and U; these decoy operations consume parent letters and do not require an empirical formula or mass.

## Independent deductions

The additional exact strings and state checks in [decoy_reference.rs](../tests/decoy_reference.rs) are deductions from source control flow. They are marked separately from literal upstream assertions:

- Zero attempts stores an unchanged peptide in the cache. Subsequent calls, including after reseeding, reuse that value without random draws. Conversely, a cached shuffled value is returned even when a later call requests zero attempts.
- The isolated `TESTRPEPTR` result `ETEPTSRRTP` follows from the first product of the literal outer Trypsin golden. Caching it with `no cleavage` and later applying Trypsin to `TESTRPEPTRIDE` reuses that final-product choice in a nonfinal context. With zero attempts for the remaining product, the result is `ETEPTSRRTPIDE`.
- A two-letter final product `AG` cannot improve its maximum forward/reverse identity. A direct 100-attempt call nevertheless consumes 100 draws. A zero-attempt call has identical output/cache but different RNG state; reseeding restores exact state equality. For `AKR` with Trypsin, the shuffled ranges have length one and consume no draws.
- Final products are fully reversed even if they end in a cleavage residue: `AK` becomes `KA`. The same nonfinal last-letter anchor is used with an N-terminal enzyme. Asp-N gives `AC | DPEF | DGK`, hence `AC | EPDF | KGD`.
- Source unspecific AASequence digestion enumerates `A,B,C,AB,BC,ABC` in length-then-start order. Reversal yields `ABCABBCCBA`; zero-attempt shuffling yields `ABCABBCABC`. Outer shuffling passes through the five short products and redigests the last `ABC`, yielding 17 letters with a fixed 14-letter prefix and a final permutation of `ABC`. The test deliberately does not assert an unobserved random permutation.
- Outer short products and complete variants can remain unchanged or duplicate one another. No deduplication or guaranteed target/decoy difference is imposed.

These checks follow [DecoyGenerator.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/DecoyGenerator.cpp) and the AASequence overload in [ProteaseDigestion.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/ProteaseDigestion.cpp). They do not duplicate the production RNG or shuffle implementation in the test suite.

## Source boundaries and native checks

Whole-protein reversal preserves an empty sequence. Outer shuffling returns one empty sequence per positive variant count, or no outputs for zero. The source peptide methods index an absent final product for empty input; native checked errors replace that undefined access. Native unsigned attempt/factor arguments represent the source's nonpositive loop behavior with zero.

All four transformations reject known, numeric residue, and numeric terminal modifications, including zero-factor/zero-attempt requests. The implementation's `!isModified()` precondition controls this behavior; header prose saying modifications are discarded does not override it. Source builds that disable preconditions need not exhibit the same rejection.

The cache key contains only parent peptide text. It intentionally omits enzyme, final/nonfinal context, attempts, and seed. Exclusive mutable access replaces the source's unsynchronized shared RNG behavior; concurrent source schedules and time-based default seeds have no reproducible golden.

The scientific audit also checks that resource limits cover overlapping unspecific products, nested outer digestion, all variant outputs, attempts, rejected random draws, cache comparison/copy costs, sequence allocation, and existing plus staged cache contents. Stateful output construction uses a copied RNG and new-entry delta, then commits only after successful output construction. Separate implementation tests cover late work/allocation failures and persistent cache limits; the independent suite checks unchanged state at public errors and cache hits.

## RNG evidence and limitations

The pinned [MathFunctions.h](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/MATH/MathFunctions.h) selects `boost::mt19937_64`, descending Fisher–Yates swaps, and Boost's inclusive integer distribution. Local Boost 1.90 headers supply supplemental implementation evidence, with separate hashes and Boost Software License 1.0 attribution. They are not part of the pinned OpenMS snapshot and do not establish which Boost release originally generated the class-test literals.

The packaged [integer reference generator](../tools/generate_decoy_rng_reference.py) independently reproduces the private RNG constants in [decoy_random.rs](../src/chemistry/decoy_random.rs). Run `python3 tools/generate_decoy_rng_reference.py` from the repository root to print its JSON; it needs only Python's standard library and does not call the native implementation. Its append-only integer recurrence uses prior sequence indices, whereas the native engine twists a cyclic array in place. The oracle omits Boost's normalization of the initial word's redundant low bits because those bits are replaced before use; separate native tests check normalization and the all-zero-state repair.

For each seed **0, 4711 and 18446744073709551615**, the oracle computes **1,248 consecutive words**, covering four complete 312-word twist cycles. Twenty zero-based checkpoints bracket midpoint and cycle boundaries; a 64-bit FNV-1a checksum over all words' little-endian bytes also checks the intervening output. Exact integer bucket division with rejection and descending swaps then derives three permutations of `ABCDEFGHIJKLMNOPQRST`: `ONLTGBCHMEJQFPRIKASD`, `IHEDJOFCGMSTAKRPNLBQ`, and `KEJPLCMBGQSRFDHOITNA`, respectively. Each uses 19 draws. These words, checksums and permutations are derived references, not additional literal OpenMS assertions. The manifest records the generator hash, checkpoint indices, checksum convention and per-seed summaries; the complete checkpoint values remain in the Rust tests and reproducible generator output without a duplicate fixture.

A separate primary-source literal anchors the raw engine: the C++ working draft specifies **9981545732273789042** as the 10,000th output of default-seeded `mt19937_64` (seed **5489**). Both the script and the native test check that exact integer. [C++ working draft, rand.predef paragraph 4](https://eel.is/c++draft/rand.predef#4).

The **nine private RNG tests** cover these multi-cycle references, the standard literal, clone/reseed state, Fisher–Yates order and draw counts, no-draw empty/singleton ranges, bucket boundaries including the full-width range, forced rejection and exhausted allowances, partial private progress before failure, and Boost normalization. Test-only inverse tempering supplies adversarial raw words without adding an injectable public RNG. Every accepted or rejected draw consumes the shared allowance before advancing the engine; the outer generator stages state and output to preserve public atomicity.

The derived helper retains `Copyright Jens Maurer 2000-2001` and `Copyright Steven Watanabe 2010, 2011`, with `BSD-3-Clause AND BSL-1.0` identification. The complete Boost Software License and component attribution belong to the package's [license records](../LICENSES.md); the independent Python reference script remains BSD-3-Clause.

Passing these exact fixtures establishes the stated source examples and independently tested state/ordering cases. It does not establish equivalence for every random seed, every source Boost release, platform clock, concurrent schedule, or unbounded input.

## Validation status

All seven independent references pass in the final Cargo runs on Rust 1.96 and 1.85 with all features and without default features, and in the current-compiler idXML-only run. The combined commands also pass all nine RNG tests, three private state tests, seven direct decoy tests and the enabled workflows. Strict all-target Clippy, rustdoc and formatting pass. The TSV, all source hashes and packaged integer-oracle output are independently verified. No C++ code was built or executed.
