# Peptide indexing

The native `analysis::peptide_indexing::PeptideIndexing` maps modified peptide hits to FASTA proteins, replaces evidence and target/decoy annotations, and reconstructs protein hits for each identification run. It uses existing native identification records and the enzyme registry. No C++ library or additional dependency is required.

Source reference: [PeptideIndexing.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/ID/PeptideIndexing.cpp), [AhoCorasickAmbiguous.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/ID/AhoCorasickAmbiguous.cpp), and [EnzymaticDigestion.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/EnzymaticDigestion.cpp), pinned at `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. Source hashes and independent oracle methodology are in [fixture provenance](../tests/data/peptide_indexing_matcher_reference_provenance.json).

## API and defaults

`run(database, protein_identifications, peptide_identifications)` accepts a borrowed FASTA slice and mutable identification slices, returning an `IndexingReport`. Successful runs replace records in place. Validation, matching, resource limits, and error policies complete before either slice is changed. Warnings are returned as strings in the report; the library does not print messages.

| Setting | Default and behavior |
| --- | --- |
| `enzyme`, `specificity` | `None`: resolve from each protein run's search parameters; unknown values fall back to Trypsin/full with warnings |
| `max_ambiguities` | 3; allowed range 0–10 |
| `max_mismatches` | 0; allowed range 0–10 |
| `il_equivalent` | false |
| `allow_nterm_protein_cleavage` | true; accepts removal of M or MX from a protein starting M |
| `decoy_rule` | `Auto`; explicit nonempty `Prefix(String)` and `Suffix(String)` are available |
| `unmatched_action` | `Error`; alternatives `Warn` retain unmatched hits, `Remove` remove hits and retain their identification records |
| `missing_decoy_action` | `Error`; alternatives `Warn` and `Silent` |
| `write_protein_sequence`, `write_protein_description` | false; true copies original FASTA fields to reconstructed protein hits |
| `keep_unreferenced_proteins` | false; true preserves unmatched old protein hits within each run, clearing their obsolete target/decoy labels |

`find_matches(peptide, protein)` exposes the raw ambiguity/mismatch matcher without enzyme filtering. It accepts ASCII residue letters, removes `*` from both inputs, and returns ascending zero-based starts in the transformed protein. It does not calculate chemistry. An empty peptide after removing stops is an error; a peptide longer than the protein has no matches.

`resolve_decoy_rule(database)` returns the effective affix, prefix/suffix position and whether automatic inference succeeded. `ResolvedDecoyRule::is_decoy` performs the corresponding case-sensitive classification.

`run_with_progress(database, proteins, peptides, &mut ProgressLogger)` is `run` reporting progress as the source's `ProgressLogger` base does; the caller passes the logger for the call and selects its type on it, and `run` itself reports nothing, as a source object of the default type `NONE`. See [Progress](#progress).

## Matching and evidence

Every accepted occurrence is retained, including overlaps and repeated occurrences within one protein. Evidence order follows FASTA order, then ascending position. Peptide modifications, scores, ranks, charge, hit order and unrelated metadata remain intact. Matching uses unmodified residue sequences; it does not require a protein modification to correspond to a peptide modification.

Protein-side ambiguous residues have the following expansions:

| Protein code | Compatible peptide letters |
| --- | --- |
| B | D, N |
| J | I, L |
| Z | E, Q |
| X | All 22 canonical residues, including O and U; excludes B, J, Z and X |

Literal case-insensitive equality costs neither ambiguity nor mismatch budget, including literal B/J/Z/X equality. Each unequal compatible protein expansion consumes ambiguity budget if available, then mismatch budget. Other unequal pairs consume mismatch budget. Equivalently, if A is the number of compatible unequal expansions and M the number of other unequal pairs, a match requires `M + max(0, A - max_ambiguities) <= max_mismatches`. Expansion is directional: a canonical protein residue does not expand an ambiguous peptide letter. Insertions and deletions are not supported.

When I/L equivalence is enabled, protein L and J become I; peptide L becomes I, while peptide J stays J, matching the source implementation. This option is rejected for Chymotrypsin, Chymotrypsin/P and TrypChymo. Both inputs are uppercased before these substitutions and enzyme validation.

Evidence starts and inclusive ends refer to the protein **after removing every `*`**. Flanks also come from the uppercase, optionally I/L-normalized protein. Protein boundaries use the native NTerminus/CTerminus markers; X flanks use Unknown. Enzyme validation can temporarily extend an M/MX-clipped match to position zero, while the evidence retains the actual matched start and flanks.

When `write_protein_sequence` is true, output protein sequences retain the **original FASTA text**, including stops and original I/L characters, as in C++. Consequently evidence offsets may not index the original stored sequence directly when it contains `*`, and normalized flanks may differ from its original spelling. Apply the documented transformations before comparing positions to stored FASTA text.

## Enzyme and run behavior

Registered enzymes support full, semi and nonspecific terminal checks. Missed-cleavage counts do not limit indexing matches, matching `FoundProteinFunctor`. X!Tandem permits additional D|P termini, even for the no-cleavage enzyme. MS-GF+/MSGFPLUS changes Trypsin to Trypsin/P. Relevant `SE:` search metadata is recognized, including original-engine recovery for Percolator and ConsensusID using the pinned first-key rule. Arbitrary digestion expressions are unsupported; an automatic run whose stored expression disagrees with the registered enzyme errors. An explicit enzyme selection overrides that expression.

Each hit receives `target_decoy` as target, decoy or target+decoy. Unmatched retained hits have this key removed. `protein_references` is unique for one matching FASTA accession, non-unique for multiple accessions, or unmatched for none; several occurrences within a single accession remain unique.

Matching protein hits are constructed anew in FASTA order with default score/rank/coverage/modifications and new target/decoy metadata, as in C++. Previous matched protein scores and unrelated hit metadata are not retained. Optional old unmatched hits retain their order before the newly matched hits. Run metadata, search parameters, score direction, protein groups and indistinguishable groups are preserved. Groups can consequently reference proteins removed by indexing; use the separate identification cleanup APIs when those groups should be pruned or recomputed.

All ten source `PeptideIndexer:*` settings are written to `run.search_parameters.metadata`, including effective enzyme/specificity and resolved decoy naming. They are not written to the run's top-level metadata. The report provides effective per-run settings, target/decoy/mixed/unmatched and uniqueness counts, output evidence/protein counts, work and warnings.

## Decoy inference

Automatic inference follows the pinned DecoyHelper implementation, including its unusual fixed regular expressions. Prefixes use its ordered vocabulary of decoy, dec, reverse, rev, reversed, `__id_decoy`, xxx, shuffled, shuffle, pseudo and random, followed by optional underscores. Alternation order means `reversed_A` is recognized with prefix `reverse`. Suffix expressions require a leading underscore and generally repeat the final letter rather than underscores: `_deco` and `_decoyyy` are recognized, while `_decoy_` is not. This behavior is preserved, not generalized to a different naming heuristic.

Inference counts names case-insensitively, requires the source's inclusive 40% database/80% observed-affix thresholds, rejects equal aggregate prefix/suffix counts, and checks eligible prefixes before suffixes. Returned spelling comes from the last observed occurrence; actual classification then remains case-sensitive. Failed inference returns the `DECOY_` prefix and `run` records a warning. Supplying an explicit affix avoids inference.

The missing-decoy policy concerns **matched peptide hits**, not merely the presence of decoy entries in FASTA. At least one decoy or target+decoy peptide hit is required by the default policy. Empty hit collections are exempt.

## Checked changes and limits

The following are deliberate differences from the pinned implementation:

- Inputs remain unchanged on all returned errors. C++ can return an error after partially replacing records.
- Run identifiers must be unique and every peptide identification must reference an existing run. FASTA accessions must be nonempty and unique. Unknown runs are not silently assigned to the first run.
- Automatic enzyme/specificity resolution and engine exceptions apply separately to each run. C++ assumes common settings and can propagate an engine exception from one run to all others.
- Old unreferenced protein hits are determined separately per run. C++ uses a global raw-match lookup, which can remove a protein from an unrelated run.
- All letters are uppercased before matching, I/L substitution and digestion. C++ case folding happens after its uppercase-only I/L substitutions, and enzyme checks may see lowercase letters.
- After stop removal, punctuation, whitespace, NUL, non-ASCII characters and modification notation in raw sequences are rejected. The C++ matcher can skip invalid protein bytes while computing incorrect offsets for matches spanning those bytes.
- Stale unmatched target/decoy labels are cleared. Kept orphan proteins use native Unknown through absence of the key, rather than the source's empty-string label.
- Empty peptide collections and collections containing only empty hit lists consistently perform protein-hit cleanup and return a successful report. An actual empty peptide hit is rejected.

The matcher is a direct, bounded scalar implementation, with cached results per distinct unmodified peptide and run. It scans protein-major, as the source's Aho-Corasick loop does: every distinct peptide is normalized once, in first-occurrence order, then every peptide is searched in one protein before the next protein, so each peptide's mappings still come out in FASTA order, then position order. It does not reproduce C++ trie construction, parallel execution, streaming FASTA cache or X-run performance shortcuts. Worst-case raw search work scales with distinct peptide count, database residues and peptide length; each raw candidate also invokes the shared enzyme validator. The optimized ambiguity trie remains porting work for proteome-scale throughput, not a claim of completed performance parity.

| Limit | Default |
| --- | --- |
| `max_records` | 1,000,000 each for FASTA entries, runs and peptide IDs; also 1,000,000 combined input peptide/protein hits |
| `max_residues` | 100,000,000 combined raw FASTA sequence bytes and all input peptide/protein-hit sequence bytes, checked before identification record clones |
| `max_matches` | 1,000,000 each for raw candidates per peptide/protein pair, cached accepted mappings, total output evidences and total output protein hits |
| `max_work` | 200,000,000 normalization bytes, residue comparisons, protein visits and charged enzyme-validation residues |

All configurable limits must be positive and use checked accumulation. Raw candidates are bounded before enzyme filtering, so a low match limit can error even when enzyme checks would reject those candidates. The shared digestion validator also enforces its one-million-residue protein length limit. Metadata and other existing record fields are validated and copied transactionally; the residue budget is not a total byte limit for arbitrary metadata. Limits may be raised explicitly for a known workload, but the scalar algorithm remains slower than the source trie.

Peptide records accept native `AASequence` with B/Z/X and named or numeric annotations; indexing uses its bare residues and does not require a formula or mass. Ambiguity expansion remains protein-side only, including when peptide records contain ambiguous letters. Literal stops are outside `AASequence`; raw `find_matches` additionally supports stop preprocessing. Invalid numerical/metadata fields in existing records are errors, even when indexing would replace some associated fields.

## Progress

The source's in-memory `run` makes two progress sections (`PeptideIndexing.cpp:371-373`, `:453-576`), and `run_with_progress` makes the same calls:

1. `Load first DB chunk` (`LOAD_PROGRESS_LABEL`) over `0` to `1`, with no set. The source loads its first FASTA chunk here; an in-memory database needs no loading, but the section is made, as the source makes it.
2. `Aho-Corasick` (`SCAN_PROGRESS_LABEL`) over `0` to the number of database entries, or to `i64::MAX` when there are exactly `SOURCE_PROTEIN_CACHE_SIZE` (400,000) of them, the source's test for a first chunk of unknown total (`:453`). The value is set after each protein is scanned, from `1` to `n` (`:490-496`). The section is skipped when there is no peptide hit to search for, as the source returns before it (`:381-394`, `:433-437`).

The source sets progress only from OpenMP thread 0; the values reproduced are those of a single-threaded run, where thread 0 scans every protein. `tests/progress_consumers.rs` replays the Release build's calls for a three-protein database, for an empty peptide list, for identifications without hits and for an empty database (tier 1, `../oracle/progress-consumers`), and checks the `i64::MAX` range at 400,000 entries. The indexing result of `run_with_progress` is the one `run` returns.

Differences:

- An empty database is refused with the rest of the up-front validation, before any progress; the source makes the first section and then returns `DATABASE_EMPTY` (`:375-379`).
- A failure inside a section still ends it (see `ProgressReporter::section`); the source's run cannot fail inside one.
- The source writes `Merge took ...` and a memory-usage line straight to `std::cout` after the scan (`:577-579`); those are not progress output and are not reproduced.

The protein-major scan changes one thing a caller can observe: when an input would be refused for more than one reason, which refusal is reported. Every peptide is now normalized before any protein is scanned, and the whole scan runs before any evidence is built, so an invalid peptide residue is reported ahead of a work or match limit that an earlier peptide's scan would have reached, and a scan limit ahead of an evidence limit. Inputs with a single cause, and every successful run, are unaffected: the same operations are charged to the same budgets, in total, and every mapping list is the same.

## Validation

`tests/peptide_indexing.rs` covers 121 independent matcher oracle cases, 17 decoy naming/threshold cases and 15 enzyme boundary cases. It also tests modified sequences, overlaps, evidence ordering, original FASTA preservation, per-run reconstruction, engine recovery, settings metadata placement, source protein-hit score resets, unmatched/missing-decoy policies, resource limits and transactional failures. Fixtures combine literal pinned upstream expectations with independently hand-traced source branches; no C++ build or execution is claimed.

The module's 17 tests pass with all features on the current Rust toolchain and Rust 1.85. Targeted Clippy passes with warnings denied. The combined protein workflow additionally exercises indexing with FASTA input, protein inference/FDR, idXML and identification splitting.
