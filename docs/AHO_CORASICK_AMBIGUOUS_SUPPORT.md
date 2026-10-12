# AhoCorasickAmbiguous support

Source header, core pin `bc9cc12`: `src/openms/include/OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h`
(`AhoCorasickAmbiguous.cpp`). Rust: `src/analysis/id/aho_corasick_ambiguous.rs`.
Tests: `tests/aho_corasick_ambiguous.rs` (with the Release oracle fixtures `tests/data/aho_corasick_ambiguous_{cases,oracle}.tsv`). Provenance:
`tests/data/aho_corasick_ambiguous_provenance.json`.

The ledger's direct TOPP consumer is `Digestor`; upstream the trie is the matcher
inside `PeptideIndexing`. The existing `src/analysis/peptide_indexing.rs` does not
use this module: it has its own bounded scalar matcher with the same ambiguity
and mismatch semantics (see `docs/PEPTIDE_INDEXING_SUPPORT.md`). Switching it
onto `ACTrie` is left to a later change, because its tests pin the current
matcher's refusal order.

## API mapping

| Source | Rust |
|---|---|
| `constexpr char const AAtoChar[28]` | `AA_TO_CHAR: [u8; 28]` |
| `constexpr char const CharToAA[256]` | `CHAR_TO_AA: [u8; 256]` (built by a `const` block from `AA_TO_CHAR`; the test checks the layout) |
| `struct Hit` (`T = uint32_t`, default and 3-arg ctor, `operator==`) | `Hit { needle_index, needle_length, query_pos }` (`u32`), `Hit::new`, `Default`, `PartialEq` |
| `struct AA` | `AA` (`Copy`, `Ord`, `Hash`) |
| `AA()` | `AA::default()` = `AA::INVALID` |
| `explicit AA(char)` | `AA::from_byte(u8)` |
| `uint8_t operator()()` | `AA::code` |
| `operator==`, `operator<=` | derived `PartialEq`, `PartialOrd` (code order) |
| `isAmbiguous`, `isValid`, `isValidForPeptide` | `is_ambiguous`, `is_valid`, `is_valid_for_peptide` |
| `operator++()`, `operator++(int)` | `AA::succ -> Option<AA>` (`None` past `?`) |
| `operator-(AA)` | `impl Sub for AA` (wrapping) |
| `toChar` | `AA::to_char` |
| `static AA fromIndex(size_t)` | `AA::from_index -> Result<AA>` |
| `static size_t unambiguousAACount()` | `AA::unambiguous_aa_count` (22) |
| `class Index` (`T = uint32_t`), default invalid | `Index`, `Default` invalid (`u32::MAX`) |
| `Index(T)` | `Index::new`, `From<u32>` |
| `isInvalid`, `isValid` | `is_invalid`, `is_valid` |
| `T operator()()`, `T pos() const` | `Index::pos` |
| `T& pos()` | `Index::pos_mut` |
| `operator==`, `std::hash<Index>` | derived `PartialEq`, `Hash` |
| `struct Bitset { uint32_t bits; }` | `Bitset { bits: u32 }` |
| `set`, `clear`, `test`, `pop_count` | same names (`u32` positions) |
| `operator<<`, `operator<<=` | `impl Shl<u32>`, `impl ShlAssign<u32>` |
| `struct ACNode` (`ACNode()`, `ACNode(AA, uint8_t)`) | `ACNode` (`Default`, `ACNode::new`) |
| `ACNode::DepthHits` (1-bit `has_hit`, 7-bit `depth`) | `DepthHits { has_hit: bool, depth: u8 }` |
| `ACNode::ChildCountType` | `ChildCountType = u8` |
| `children_bitset`, `suffix`, `first_child`, `edge`, `nr_children`, `depth_and_hits` | same field names |
| `struct ACScout` (no default ctor; 5-arg ctor) | `ACScout::new` (no `Default`) |
| `ACScout::it_query` (`const char*`) | `it_query: usize` (offset into the state's query) |
| `tree_pos`, `max_aaa_leftover`, `max_mm_leftover`, `max_prefix_loss_leftover` | same |
| `size_t textPos(const ACTrieState&) const` | `ACScout::text_pos(&self)` (the offset is the position) |
| `AA nextValidAA()` | `ACScout::next_valid_aa(&mut self, query: &str)` (the query is passed; there is no pointer) |
| free `AA nextValidAA(const char*&)` | `next_valid_aa(&[u8], &mut usize) -> AA` |
| `struct ACTrieState` | `ACTrieState` (`Default`) |
| `setQuery`, `textPos`, `getQuery`, `nextValidAA` | `set_query`, `text_pos`, `get_query`, `next_valid_aa` |
| `const char* textPosIt() const` | `text_pos_it(&self) -> &[u8]` (the rest of the query) |
| `hits`, `scouts` (`std::queue`), `tree_pos` | `hits: Vec<Hit>`, `scouts: VecDeque<ACScout>`, `tree_pos` |
| `class ACTrie` | `ACTrie` |
| `ACTrie(uint32_t max_aaa = 0, uint32_t max_mm = 0)`, `~ACTrie()` | `ACTrie::new`, `Default` (0, 0), drop |
| `addNeedle`, `addNeedles`, `addNeedlesAndCompress` | `add_needle`, `add_needles`, `add_needles_and_compress` (`Result`) |
| `size`, `compressTrie`, `getNeedleCount` | `size`, `compress_trie`, `get_needle_count` |
| `setMaxAAACount`, `getMaxAAACount`, `setMaxMMCount`, `getMaxMMCount` | `set_max_aaa_count`, `get_max_aaa_count`, `set_max_mm_count`, `get_max_mm_count` |
| `bool nextHits(ACTrieState&) const` | `next_hits(&self, &mut ACTrieState) -> Result<bool>` |
| `void getAllHits(ACTrieState&) const` | `get_all_hits(&self, &mut ACTrieState) -> Result<()>` |
| private `nextHitsNoClear_`, `add_`, `addHits_`, `addHitsScout_`, `follow_`, `followScout_`, `stepPrimary_`, `stepScout_`, `createScouts_`, `createSubScouts_`, `createMMScouts_`, `createMMSubScouts_`, `findChildNaive_`, `findChildBFS_` | private methods of the same meaning (`add_hits` and `addHitsScout_` share `collect_hits`) |
| `.cpp` `_ambiguitiesOf`, `countSetBitsUpTo` | private `AA::ambiguities`, `count_set_bits_up_to` |

## Preserved source conventions

- The code order of `AAtoChar` (unambiguous `A`..`V` at 0..21, then `B J Z X`,
  `$`, `?`), case-insensitive letters, and `?` for every other byte including
  128..255.
- Ambiguity ranges: `B` = `D..N`, `J` = `I..L`, `Z` = `E..Q`, `X` = `A..V`,
  `$` = `A..X`; mismatches try every residue from `A` to `X` except the
  ambiguity range and the query residue itself, in code order.
- The breadth-first layout with children in edge order, the two-pass suffix
  links (first plain, then compressed where the suffix has no hit and no extra
  child), and the hit flag inherited along suffix links.
- The primary and scout walk exactly as the source: scouts are processed FIFO
  after each primary step, each runs until it loses its spawning residue or the
  query ends, and may queue sub-scouts; a scout only reports hits that include
  its spawning residue.
- `next_hits` returns `true` only when the *primary* reports; scout hits
  accumulate in `state.hits` until then, and those found after the primary's
  last hit are left in `state.hits` with a `false` return (a `while nextHits`
  loop drops them, as upstream's does; `get_all_hits` keeps them).
- Hit start = end of match minus needle length, so it is off by the number of
  skipped invalid characters (`*`) inside the match (the class test documents
  this).
- A NUL byte ends a query (the source walks a C string).
- `add_needles` keeps the needles before a refused one, so `get_needle_count`
  identifies it (the class test asserts this).

## Native differences

- **Refused instead of undefined** (`Error::InvalidValue`): an empty needle (the
  source marks the root as a hit and its hit collection then loops forever); a
  needle longer than `MAX_NEEDLE_LENGTH` = 127 (the seven-bit depth wraps); adding
  a needle after compression (the source appends to its breadth-first array with
  cleared construction lists); searching before compression (a debug-only
  assertion upstream); a state that was never given a query; a query or node
  count beyond `u32`.
- **Atomic needle insertion:** a refused needle leaves the trie unchanged. The
  source has already created the nodes of the needle's valid prefix when it
  throws, which `size()` shows, and those orphan nodes are past the end of its
  sparse needle table, which `compressTrie` then reads out of bounds.
- **A second `compress_trie`** is a no-op. Upstream it rebuilds from the
  discarded construction lists, keeping only a root that still has its old child
  bits, and the next search reads out of bounds.
- **Budgets above 255** saturate when copied into a scout; the source narrows
  `uint32_t` to `uint8_t`, so `max_aaa = 256` becomes 0 inside scouts.
- **`AA::succ`** returns `None` past `?` (the source asserts in debug builds and
  reads out of bounds otherwise); `Bitset` positions of 32 or more are ignored
  (undefined shifts upstream).
- `AA::from_index` returns an error for index 22 and above (a debug-only
  precondition upstream).
- `compress_trie` does not log its node and compression counts (the source logs
  them at info level).
- `ACScout` and `ACTrieState` hold offsets, not raw pointers into the query.

## Source defect candidates

1. An empty needle makes `addHits_` loop forever (the root is its own suffix and
   has `has_hit` set).
2. Needles of 128 residues or more wrap the seven-bit depth, so hits report a
   wrong length and start.
3. A needle refused for an invalid character leaves orphan nodes that are past
   the end of `vec_index2needles_`; `compressTrie` then indexes it out of range.
4. Calling `compressTrie` twice corrupts the trie.
5. Budgets above 255 wrap to small numbers inside scouts.

## Checked boundaries and evidence

Tier 3: every section and literal of `AhoCorasickAmbiguous_test.cpp` at
`bc9cc12`: the 22 `testCase` hit multisets over 13 tries (exact, ambiguous,
mismatch, prefix/suffix, `*` offsets, the scout first-residue rule), the needle
refusals and counts, `setQuery`, the `nextValidAA` position sequences, and the
`AA` table assertions. Tier 4 (native): table layout, neighbour ranges, the small
types, `get_all_hits` against the `next_hits` concatenation, and the refusals
above.

Tier 1 (executed differential): the class test compares sorted hit lists, so it
cannot see hit *order*. `release_oracle_hit_order` replays 322 cases (the class
test's 22 tries and queries, plus 300 seeded random tries with ambiguous needle
residues and queries with ambiguous residues, lowercase letters, `*` and empty
queries, at `max_aaa` 0-4 and `max_mm` 0-3) against the output of the Release
libOpenMS (`bc9cc12-c19e494-174b576`, Linux x86_64, run on ibminode06): the trie
size after compression, every `nextHits` result in library order, the hits left
in the state after the final `false`, and `getAllHits`. All 1186 records match
exactly. The driver ran twice per invocation and the oracle was invoked twice;
all four outputs are byte-identical. Driver, generator and raw output are under
`../oracle/I3b-aho/`, hashed in the provenance manifest.

No OpenMP in the source; the trie is immutable after compression and one state
per thread is the source's threading model, which `&self` search preserves.
