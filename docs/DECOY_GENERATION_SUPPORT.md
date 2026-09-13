# Native decoy generation

`chemistry::DecoyGenerator` ports the complete scientific surface of the pinned
OpenMS4-core `DecoyGenerator`: protein reversal, peptide reversal, stateful
peptide shuffling, and deterministic multi-variant shuffling. It accepts an
`AASequence` and the existing `Protease` enum. `Protease::from_name` resolves
the pinned enzyme names and synonyms.

All transformations require unmodified input, matching the source precondition.
Named modifications, terminal annotations and numeric mass tags are rejected.
This deliberately follows the checked source behavior rather than its header's
contradictory note about discarding modifications. Parent residue letters,
including U/O/J/B/Z/X, can be rearranged without requiring a molecular formula.

## API and source conventions

| Operation | Behavior |
| --- | --- |
| `with_seed(seed)` | Reproducible stateful generator |
| `Default::default()` | Time-seeded convenience; use an explicit seed for reproducibility |
| `set_seed(seed)` | Reseeds the random engine; retains the peptide cache |
| `reverse_protein(protein)` | Reverses every parent residue; empty input gives an empty sequence |
| `reverse_peptides(protein, enzyme)` | Digests with full specificity and zero missed cleavages, reverses all but the last residue of each non-final product, and fully reverses the final product |
| `shuffle_peptides(protein, enzyme, max_attempts)` | Stateful, cached peptide shuffling using the same positional anchoring rule |
| `shuffle(protein, enzyme, decoy_factor)` | Produces whole-sequence variants using a fresh generator seeded with `4711 + variant` for each longer digest product; leaves the receiver unchanged |

Peptide reversal and stateful shuffling reject empty input because the source
indexes a missing final product. Multi-variant shuffling instead returns one
empty sequence per requested variant for empty input, and zero variants returns
an empty vector. Zero shuffle attempts returns and caches the original products.
Native attempt and variant counts are unsigned.

Anchoring is positional, including for enzymes that cleave at the N terminus.
The last product is fully reversed or shuffled even if the protein ends at a
cleavage residue. Multi-variant shuffling processes each outer product in a
fresh, isolated inner call. Its final cleavage residue can therefore move.
Products of at most two residues pass through unchanged only in that outer
operation. Direct stateful shuffling still attempts these short products.

Each attempt shuffles the previous candidate in place. The objective is the
larger of the fractions matching the target forward and backward. Only strict
improvement replaces the best candidate. Non-final products can stop after an
improvement reaches `1 / length + 1e-6`; the final product stops at exactly zero.
Rejected candidates still advance the random engine. Equal-quality candidates
retain the earlier best result.

## Persistent cache and reproducibility

The cache key is the raw peptide text. It deliberately omits enzyme, position,
attempt limit and seed, as in the source. A cache hit consumes no random draws;
repeated products within the same call consult newly staged entries too.
Reseeding does not discard cached choices. An earlier final-product result can
be reused in a non-final context, and a zero-attempt call can keep a product
unchanged on later calls with more attempts.

Stateful shuffling requires `&mut self`. Independent owned instances provide
reproducible execution without shared cache or random-state races. The outer
multi-variant operation ignores the receiver's seed and cache and does not
deduplicate identical variants or guarantee a difference from the target.

The private random helper draws MT19937-64 words from `rand_mt::Mt64` and
implements the Boost integer-bucket rejection mapping used by the source's
descending Fisher–Yates shuffle. Seeded source strings and independently derived
integer references test the mapping; using the same engine with modulo indexing
or another shuffle distribution would change the results. The inspected Boost 1.90 headers are supplemental
evidence, separate from the pinned OpenMS snapshot. Its time seeds and concurrent
source schedules are not reproducible reference cases.

## Unspecific digestion

The source AASequence digestion overload overrides zero missed cleavages for
`unspecific cleavage`, emitting all substrings in length-then-start order.
The decoy adapter preserves that ordering without changing native digestion's
existing enumeration. These products overlap, so concatenating them increases
sequence length and does not preserve the input composition.

For a protein of length `n`, the direct peptide operations concatenate
`n*(n+1)*(n+2)/6` residues. Multi-variant shuffling redigests each longer product
and can grow further. Checked output and work limits apply to these operations;
ordinary non-overlapping enzyme digestion retains composition.

## Checked native boundaries

| Limit | Per operation unless stated otherwise |
| --- | ---: |
| Input residues | 1,000,000 |
| Combined output residues | 1,000,000 |
| Digest products across outer and inner calls | 100,000 |
| Variants | 1,000 |
| Shared work, including actual random draws | 50,000,000 units |
| Cumulative allocation payload | 256 MiB |
| Persistent cached products | 100,000 entries |
| Persistent cache key/value and node allowance | 16 MiB |

Input validation precedes sequence copying. Output checks include per-residue
annotation storage, all requested variants and overlapping products. Shared
work counters cover scans, comparisons, shuffle attempts and rejected random
draws across a complete operation, including nested shuffles. Persistent cache
limits include existing entries plus staged key/value payloads.

Stateful operations clone only the small random engine and stage new cache
entries. Results are committed after every checked operation and sequence
construction succeeds. A returned error leaves the receiver's random state and
cache unchanged. Caller input sequences are immutable.

The [reference review](DECOY_REFERENCE_REVIEW.md) and
[fixture provenance](../tests/data/decoy_provenance.json) distinguish literal
source expectations from independently derived boundary cases. The
[workflow tests](../tests/decoy_workflow.rs) exercise FASTA output, peptide
indexing, target/decoy metadata, simple independently calculated FDR values and
idXML round-trips. The scores in that test are synthetic inputs, not search
results or a validation of biological confidence.

Run `cargo run --locked --offline --example generate_decoys` for a small target
and decoy FASTA example. No C++ runtime or additional Cargo dependency is needed.
