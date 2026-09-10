# RNA digestion and enzyme records

`RNaseDigestion` implements sequence digestion and identification-graph registration
for all fourteen pinned RNA enzymes, including modified nucleotide codes, missed
cleavages, length limits and terminal gains. `DigestionEnzymeRNARecord` contains the eight independent source fields;
`DigestionEnzymeRNA::from_record` checks and freezes them. `RNaseDB::global()`
contains the complete [source XML](../resources/enzymes/Enzymes_RNA.xml), projected
into Rust by [the generator](../tools/generate_rnases.py). No runtime data file,
XML parser or additional dependency is needed.

```rust
use openms::chemistry::{NASequence, RNaseDigestion};

fn main() -> openms::Result<()> {
    let mut digestion = RNaseDigestion::new("RNase_T1")?;
    digestion.missed_cleavages = 2;
    let input: NASequence = "pAUGUCGCAG".parse()?;
    let products = digestion.digest_with_positions(&input)?;
    assert_eq!(products.len(), 6);
    assert_eq!(products[0].sequence.to_string(), "pAUGp");
    assert_eq!((products[0].start, products[0].end), (0, 3));
    Ok(())
}
```

`digest` returns sequences; `digest_with_positions` additionally returns each
zero-based, half-open source interval and its missed-cleavage count.
`digest_into` replaces a supplied sequence vector only after every new product
succeeds. `new(name)` and `set_enzyme(name)` use the global registries;
`with_enzyme` and `set_enzyme_with_registry` accept an owned enzyme handle and a
caller-supplied nucleotide registry. Configuration captures gain records once,
so subsequent digestion survives registry drop. Configuration failures preserve
the previous enzyme, gains and options.

The native default is a configured RNase_T1. The C++ default constructor does
not initialize its RNA gain pointers, so source clients must call `setEnzyme`
before digestion. Native construction avoids that invalid state.

## Cleavage behavior

`cuts_after` patterns match the residue codes to the left of a boundary, in
left-to-right order. `cuts_before` patterns match right-hand residues, also in
left-to-right order. Commas separate patterns for successive residues. The XML
comments reverse these labels; this implementation follows executable C++.
Empty whole pattern lists impose no context requirement. Empty members within
a comma-separated list match any code but still occupy a context position.

Matching searches complete raw record codes. It does not substitute their
origins or parse their display text. All sixteen nonempty registered expressions
are compiled into small native predicates. Literal substring patterns and
comma compositions are also supported. Unsupported regex syntax returns
`Error::Unsupported`; an arbitrary Boost regex engine is not included. The
inherited main `regex` field is checked using the same subset at configuration,
because the C++ base setter compiles it even though RNA digestion does not use it.
All pinned RNA definitions leave that field empty.

Several source details affect modified RNA:

- RNase_T1 matches `G(?!m)`: it cuts after `m1G`, but not `Gm`.
- RNase_U2's `G|A(?!m)` still matches `Gm`; the lookahead applies only to A.
  Colicin_E5's `G|Q(?!m)` has the analogous G/Q distinction.
- RNase_MC1's `.*(?!m)$` always matches an empty suffix. It does not exclude
  a preceding methylated code.
- Cusativin's `^[^C]+$` requires the complete right-hand code to contain no C.
- MazF checks three separate right-hand codes with immediate `m6`/`m5`
  lookbehinds. Those lookbehinds never cross residue boundaries. There is no
  cut before the first residue, even when the input begins ACA.
- RNase_4 contains a literal closing-bracket alternative, whereas RNase_4p/4c
  contain `m1Y`. Ordinary display brackets are absent from raw record codes.

For proper enzymes, products are ordered by start position, then by increasing
missed-cleavage count. Increasing endpoints that exceed the maximum length are
skipped without visiting still-longer products; retained order and contents are
unchanged. Protein digestion's different enumeration order is not reused.
`no cleavage` returns at most the full input. `unspecific cleavage` enumerates
all permitted substrings by start then length, ignoring missed cleavages.
RNase_H matches every boundary but remains a proper enzyme, so missed cleavages
control which contiguous products it returns.

`min_length == 0` means one residue. `max_length == 0` means the input length;
larger maxima are clamped. Empty input and inverted or impossible length ranges
return no products. The latter avoids unsigned underflow in the source's
unspecific branch. Inherited protein-specificity modes are not consulted by the
C++ RNA algorithm and do not appear as native RNA settings.

## Terminal chemistry and identity

Digestion first slices the input, retaining the source sulfur-aware slicing
behavior. At each newly created internal end it then replaces the corresponding
terminal record with the enzyme gain, including a null gain. This can remove an
inherited five-prime sulfur group. A missing required slicing context still
fails before that replacement; no global context is substituted.

Original outside termini retain their records. Five-prime gain `p` resolves to
`5'-p`; other nonempty five-prime codes are looked up literally. Three-prime `p`
and `c` resolve to `3'-p` and `3'-c`. Other nonempty three-prime strings are wrapped
in brackets **before registry lookup**, preserving the source's unusual rule.
A caller's `cap` gain therefore requires a record whose actual code is `[cap]`.

The registry accepts owned records with `from_records`. Canonical names are
indexed both as written and in ASCII lowercase; synonyms are literal. Lookup
does not normalize arbitrary mixed-case input. Later same-name entries replace
older ones, and source alias/main-regex erasure rules are retained, including
colliding-index erasure. Remaining entries iterate in stable provider order,
replacing C++ pointer order. Returned handles survive registry destruction.
Native equality/ordering/hash includes RNA patterns and gains; the source's
inherited base comparison omits them.

## Identification-graph registration

`digest_identification_data(&mut graph)` digests all registered RNA parents in
accession order and skips parents of other molecule types. The registry-aware
`digest_identification_data_with_registry(&mut graph, &registry)` parses parent
text using caller-owned nucleotide records. Configured enzyme gains retain the
records captured when configuring the digestion object.

Products become owned `IdentifiedOligo` records in
[`identification::graph`](IDENTIFICATION_GRAPH_SUPPORT.md). Existing oligos are
retained, and equal complete chemical sequences merge their parent matches.
Distinct chemistry with equal display text remains distinct. Registration
attaches the graph's current processing step without reordering prior history.

Graph `ParentMatch` intervals are zero-based and **inclusive**, unlike the
half-open `DigestedOligo` intervals above. A product spanning `[start, end)` stores
`Some(start)` and `Some(end - 1)`. Neighbor markers are `[` and `]` at the original
parent ends; interior neighbors use the first byte of the raw nucleotide code.
For example, a code beginning `m` contributes `m`, not its base origin or full
bracketed display. A non-ASCII first byte returns a checked unsupported error
because it cannot form the source's one-byte neighbor as valid UTF-8.

The whole batch uses one bounded graph snapshot and commits only after every
parent, product and reference succeeds. A late malformed parent, unsupported
neighbor or resource failure preserves the original graph and its IDs. RNase
work, total input residues, generated products and logical bytes are shared
across all parents; parsed-parent payload and the parent-ID list also consume
the byte allowance. Graph registration independently shares one graph work and
allocation allowance across the same batch. It does not reset that allowance
for each product.

## Bounds, evidence and remaining work

Default digestion ceilings are one million input residues, 100,000 products,
50 million cumulative work units and 256 MiB conservative retained output.
The corresponding public `max_*` fields can lower these ceilings. Planning
checks product counts and plan bytes before copying sequences. Per-residue
prefix accounting charges each fragment for its own logical record payload,
with a conservative allowance for original ends, slicing context and enzyme
gains. Shared records may be counted repeatedly. This is not an allocator or
resident-memory measurement; input storage, prefix/cut/position plans and the
staged output coexist temporarily. There is no silent truncation on limits.

Record text fields allow 65,536 bytes each and 1,024 synonyms. Registry input is
bounded to 100,000 records, 128 MiB accounted records/indices and 50 million
work units, including conservative tree-comparison and index-update charges.
Cut and prefix scratch arrays each use at most about 8 MiB on 64-bit hosts;
these are separate from the retained-output ceiling. The sequence-only return
vector can overlap the position-bearing staging vector. A configured side allows at most 1,024 per-code patterns. Lookup
and matching retain finite input/work limits.

[Direct tests](../tests/rnase.rs) cover custom ownership, predicate precedence,
configuration rollback, half-open positions, sulfur/end replacement, empty and
inverted bounds, literal source index collisions and cumulative failures.
[Independent references](RNA_PROCESSING_REFERENCE_REVIEW.md) retain the original
twelve digestion cases and 38 output strings, all fourteen enzyme records, and
6,048 per-code pattern checks. No C++ library was built or executed.

[Graph references](../tests/identification_graph_reference.rs) additionally retain
the source three-product registration case and repeated matches at positions
zero and four, with independently derived inclusive ends and neighbor markers.
They also check rollback after a later malformed RNA parent.

An enzyme XML input provider and arbitrary Boost regex syntax remain outside
this implemented digestion surface. Graph persistence and full legacy conversion
remain unimplemented; observations, compounds, observation matches, parent/match
groups, referential cleanup and the sequence/evidence converter are covered;
see [graph support](IDENTIFICATION_GRAPH_SUPPORT.md).
