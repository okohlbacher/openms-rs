# Protein digestion support

`openms::chemistry::digestion` implements the 33 enzyme definitions in the pinned
[OpenMS4-core enzyme registry](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/share/OpenMS/CHEMISTRY/Enzymes.xml).
The original `Protease`, `ProteaseDigestion` and `DigestedPeptide` exports remain
available from `openms::chemistry`.

## API and coverage

| API | Behavior |
| --- | --- |
| `Protease` | Typed enzyme selection, name/synonym parsing, metadata, cleavage positions and internal-site counts |
| `ProteaseDB` | Immutable registry with the source's name index (name, lower-cased name, synonyms), expression index, name listings and `writeTSV`; see [ProteaseDB](#proteasedb) |
| `DigestionEnzymeProtein` | Name, expression, description, synonyms, terminal gains, PSI/X!Tandem IDs and optional Comet/MS-GF+/OMSSA IDs |
| `ProteaseDigestion::digest` | Modified `AASequence` peptides with original protein ranges and actual missed-cleavage counts |
| `digest_ranges` | Range products without constructing peptide strings; endpoints are always zero-based start and exclusive end |
| `digest_unmodified` | Borrowed slices into an uppercase unmodified protein string |
| `peptide_count`, `peptide_count_unmodified` | Counts matching the same specificity and length constraints as generation, without allocating peptide strings or range outputs |
| `count_internal_cleavage_sites`, `count_missed_cleavages` | All internal sites or sites strictly inside a specified protein range |
| `is_valid_product`, `is_valid_product_unmodified` | Terminal specificity and optional missed-cleavage checks with explicit `ProductValidation` allowances |
| `DigestionSpecificity`, `NAMES_OF_SPECIFICITY`, `by_name`, `value`, `from_value` | The source's `Specificity` enumeration with all six members (`SPEC_NONE` 0, `SPEC_SEMI` 1, `SPEC_FULL` 2, `SPEC_UNKNOWN` 3, `SPEC_NOCTERM` 8, `SPEC_NONTERM` 9), `NamesOfSpecificity` and `getSpecificityByName`; `digest` refuses the last three, validation treats them as full |
| `filter_by_missed_cleavages` | `filterByMissedCleavages`: the predicate receives the source's `Int` count |
| `NO_CLEAVAGE`, `UNSPECIFIC_CLEAVAGE` | `EnzymaticDigestion::NoCleavage` and `UnspecificCleavage` |
| `digest_ranges_source`, `digest_unmodified_source`, `peptide_count_source`, `count_internal_cleavage_sites_source`, `is_valid_product_source` | The source's own arithmetic, ported literally: `ProteaseDigestion::digest` (pair overload) with its `wrong_size`, `EnzymaticDigestion::digestUnmodified` (pair overload), `peptideCount`, `countInternalCleavageSites` and `isValidProduct_` with its `int` position and length and its three warnings, over `EnzymaticDigestion::tokenize_`'s positions; see [Source arithmetic](#source-arithmetic) |

Full specificity requires both peptide ends to be protein termini or enzyme
cleavage sites. Semi specificity requires at least one such end. `None`
enumerates every substring in the requested length interval. Selecting
`UnspecificCleavage` also enumerates every substring, regardless of specificity
and missed-cleavage settings. Products at different positions remain distinct
even if their sequences are identical.

The default remains fully specific trypsin, zero missed cleavages, minimum
length one and no maximum length. `Some(0)` is an invalid maximum; use `None` for
no limit. Added settings use defaults through Rust struct update syntax:

```rust
use openms::chemistry::{AASequence, DigestionSpecificity, Protease, ProteaseDigestion};

let protein: AASequence = ".(Acetyl)AC(Carbamidomethyl)KR.(Amidated)".parse()?;
let digestion = ProteaseDigestion {
    enzyme: Protease::LysC,
    specificity: DigestionSpecificity::Semi,
    missed_cleavages: 1,
    min_length: 2,
    max_length: Some(20),
    ..Default::default()
};
let products = digestion.digest(&protein)?;
# Ok::<(), openms::Error>(())
```

Cleavage depends on the unmodified residue letters, matching the source.
Modified residues remain modified in each product. Original terminal
modifications survive only when a product includes that original protein end.
Peptide formulas follow `AASequence::subsequence`; enzyme terminal-gain metadata
is not applied a second time. This also matches the source's peptide digestion.

## Registry details

The registry contains these exact names:

| Enzyme | Native variant |
| --- | --- |
| Trypsin | `Trypsin` |
| Arg-C | `ArgC` |
| Arg-C/P | `ArgCP` |
| Asp-N | `AspN` |
| Asp-N/B | `AspNB` |
| Asp-N_ambic | `AspNAmbic` |
| Chymotrypsin | `Chymotrypsin` |
| Chymotrypsin/P | `ChymotrypsinP` |
| CNBr | `CNBr` |
| Formic_acid | `FormicAcid` |
| Lys-C | `LysC` |
| Lys-N | `LysN` |
| Lys-C/P | `LysCP` |
| PepsinA | `PepsinA` |
| TrypChymo | `TrypChymo` |
| Trypsin/P | `TrypsinP` |
| V8-DE | `V8DE` |
| V8-E | `V8E` |
| leukocyte elastase | `LeukocyteElastase` |
| proline endopeptidase | `ProlineEndopeptidase` |
| glutamyl endopeptidase | `GlutamylEndopeptidase` |
| Alpha-lytic protease | `AlphaLyticProtease` |
| 2-iodobenzoate | `Iodobenzoate` |
| iodosobenzoate | `Iodosobenzoate` |
| staphylococcal protease/D | `StaphylococcalProteaseD` |
| proline-endopeptidase/HKR | `ProlineEndopeptidaseHKR` |
| Glu-C+P | `GluCPlusP` |
| PepsinA + P | `PepsinAPlusP` |
| cyanogen-bromide | `CyanogenBromide` |
| Clostripain/P | `ClostripainP` |
| elastase-trypsin-chymotrypsin | `ElastaseTrypsinChymotrypsin` |
| no cleavage | `NoCleavage` |
| unspecific cleavage | `UnspecificCleavage` |

Names are found as written and ASCII lower-cased (`trypsin`), synonyms as
declared; the lookup is otherwise case-sensitive. Examples include `Clostripain`
for Arg-C, `Glu-C` for glutamyl endopeptidase, `no_cut` for no cleavage and
`nonspecific` for unspecific cleavage. Identical cleavage expressions can belong
to multiple distinct entries, so `enzymes_by_regex` returns all matching names.
It accepts exact registered expressions only; it does not compile user regexes.
Unknown expressions return `Unsupported`.

Rules follow the actual source expressions, including explicit X/B/Z/J classes,
not simplified biological descriptions. For example, Formic_acid requires D/B/X
on both sides of the bond; proline endopeptidase requires an H/K/R/X followed by
P/X on the left and excludes P on the right. Iodosobenzoate cuts after W only,
while 2-iodobenzoate also includes X. Trypsin excludes all cleavage before P,
including WKP and MRP; this revision has no special exceptions for those motifs.

Unmodified string operations accept uppercase A–Z, including ambiguous codes.
`AASequence` also preserves B/Z/X and numeric annotations; its digestion products
do not require calculable chemistry. String operations reject modifications, spaces,
separators, lowercase letters and non-ASCII text. Sequence termini are always
cleavage boundaries; duplicate endpoint cuts and empty peptides are not emitted.
An empty input has the boundary list `[0]`, zero internal sites and zero products.

## ProteaseDB

`src/chemistry/digestion/protease_db.rs` ports `CHEMISTRY/ProteaseDB.h` and the
lookup half of its header-only base `CHEMISTRY/DigestionEnzymeDB.h` at core
`bc9cc12`. The source singleton loads the compiled-in
`BuiltInProteaseDataProvider`, then `share/OpenMS/CHEMISTRY/Enzymes.xml` when it
is found, whose definitions replace the built-in ones by name. At the pin both
carry the same 33 definitions in the same order, so the port replays both
registration passes (`DigestionEnzymeDB::addEnzyme_`, including its
remove-by-found-key rule) over the one compiled table to build the same name and
expression indices.

### API mapping

| Source member | Rust |
| --- | --- |
| `getInstance()` | `ProteaseDB::global()` |
| `ProteaseDB(std::vector<std::unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeProtein>>>)` | **not ported**: needs owned user-defined enzymes with an executable cleavage expression; the compiled `DigestionEnzymeProtein` has neither (see below) |
| `getEnzyme(name)` | `get_enzyme(name) -> Result<&'static DigestionEnzymeProtein>` |
| `getEnzymeByRegEx(regex)` | `get_enzyme_by_regex(regex)`; `enzymes_by_regex` (native) returns every enzyme with the expression |
| `getAllNames(names)` | `names()` |
| `hasEnzyme(const std::string&)` | `has_enzyme(name)` |
| `hasRegEx(regex)` | `has_regex(regex)` |
| `hasEnzyme(const DigestionEnzymeProtein*)` | `contains(&enzyme)`, by address as the source compares pointers |
| `beginEnzyme()` / `endEnzyme()` | `iter()`; `enzymes()` is the compiled table |
| `getAllXTandemNames` / `getAllCometNames` / `getAllOMSSANames` / `getAllMSGFNames` | `xtandem_names()` / `comet_names()` / `omssa_names()` / `msgf_names()` |
| `writeTSV(filename)` | `write_tsv(&mut impl Write)`, `write_tsv_file(path)` |
| `ConstEnzymeIterator`, `EnzymeIterator` typedefs | the `iter()` iterator |
| protected `addEnzyme_`, `loadFromProviders_`, the three index members | private `RegistryIndex` |

### Preserved source conventions

- Every enzyme is found under its name, `StringUtils::toLower(name)` (ASCII) and
  each synonym. When two enzymes claim one key the later registration wins.
- `getEnzymeByRegEx` returns the one enzyme the expression index keeps, the last
  registered: `(?<=[MX])` is `cyanogen-bromide`, not `CNBr`. An empty expression
  is never indexed.
- `writeTSV` writes the header `OpenMS_AllowedEnzymes` and one name per line.
- An enzyme without an engine identifier (`-1` in the source) is left out of that
  engine's listing; MS-GF+ identifier 0 counts.

### Native differences

- **Listing order.** The source keeps its enzymes in a
  `std::set<const DigestionEnzymeProtein*>`, so `getAllNames`, the engine
  listings, `beginEnzyme()` and `writeTSV` follow heap address order. On the
  Linux Release build that order starts `Arg-C/P, Asp-N, ..., Trypsin, Arg-C`:
  an allocator artefact that differs between platforms and allocators. The port
  lists in registration order (the `Enzymes.xml` order). The members are the
  same; the executed oracle checks them.
- `get_enzyme` and `get_enzyme_by_regex` return `Error::InvalidValue` where the
  source throws `ElementNotFound` and `IllegalArgument`.
- `write_tsv_file` reports a file that cannot be written; the source's unchecked
  `std::ofstream` writes nothing and reports nothing.
- `Protease::from_name` and `ProteaseDigestion::set_enzyme` resolve through the
  same index, so they accept a lower-cased name as source
  `ProteaseDigestion::setEnzyme` does.

### Checked boundaries and evidence

`tests/protease_db.rs` replays the executed C++ oracle `../oracle/B8b-dbs`
(Linux x86_64 Release install `bc9cc12-c19e494-174b576`, run twice on
ibminode06, byte-identical): all 33 definitions field by field, 125 name keys
(every name, lower- and upper-cased name, synonym, and unknown keys) through
`hasEnzyme` and `getEnzyme`, all 31 registered expressions plus two unknown ones
through `hasRegEx` and `getEnzymeByRegEx`, pointer membership, the five listings
(as sets, and in registration order natively) and the `writeTSV` file (the same
lines). `tests/data/protease_db_provenance.json` records the hashes. The
class-test literals of `ProteaseDB_test.cpp` are transcribed (tier 3).
`tests/digestion.rs` previously asserted that `trypsin` is unknown; the source
registers it and the oracle resolves it, so that assertion now states the
source behaviour.


## Ordering and compatibility choices

- Fully specific products follow increasing missed-cleavage count, then protein
  position, matching `ProteaseDigestion::digest`.
- Semi-specific products append after the fully specific products. They follow
  the source's alternating scan from both protein ends, extending to successive
  cleavage sites up to the missed-cleavage ceiling. Each range occurs once.
- Unrestricted products follow increasing start position, then length, matching
  `EnzymaticDigestion::digestUnmodified`. Both modified and unmodified native APIs
  use this order. The C++ `ProteaseDigestion` unspecific path instead orders by
  length; its `SPEC_NONE` path throws, despite the unmodified API supporting it.
- Counts use all generation settings. The C++ `peptideCount` method counts only
  the enzyme/missed-cleavage combinations, ignoring semi-specificity and length
  filtering. The native count therefore agrees with the actual native output.
- `ProductValidation` defaults to ignoring missed cleavages, as in the source.
  When enabled, the N-terminal allowance treats a start at residue 1 or 2 as
  protein start 0 if the protein begins with M; missed-cleavage counting includes
  that restored prefix. The D|P allowance accepts those bonds as terminal sites
  but does not add them to the enzyme's internal missed-cleavage count.
- Validation checks termini and optional missed cleavages, not the digestion
  length filter. As in the source, `None` specificity can explicitly check missed
  cleavages during validation even though unrestricted generation ignores the
  ceiling. `UnspecificCleavage` validation always ignores that ceiling.
- Missed-cleavage counts use the complete protein context. This preserves
  multi-residue recognition at range edges; it avoids the source's context loss
  when `SPEC_NONE` validity tokenizes an isolated substring.
- Empty, reversed or out-of-bounds validity ranges return false safely, including
  with the N-terminal allowance. Explicit range-counting operations reject such
  ranges. No pointer indexing precedes boundary validation.

### Leading cleavage site

Asp-N `(?=[DBX])`, Asp-N_ambic `(?=[DBEZX])` and Lys-N `(?=[KX])` are
lookahead-only (`share/OpenMS/CHEMISTRY/Enzymes.xml:43`, `:67`, `:135` at core
`bc9cc12`), so on a protein or peptide that starts with one of those residues
the expression matches with zero width at position 0. The source's
`EnzymaticDigestion::tokenize_` (`EnzymaticDigestion.cpp:92-116`) then yields
position 0 twice. On the Linux Release build
(`../oracle/integ-B-digest-leading`, ibminode06, two identical runs; retained in
`tests/data/digestion_leading_site_oracle.tsv`) this has these effects, and the
port, which keeps each boundary once, differs from the source exactly there:

- `countInternalCleavageSites` counts one more site than
  `count_internal_cleavage_sites` (`DFPIANGER` under Asp-N_ambic: 2 in the
  source, 1 here). IDFilter's `PeptideDigestionFilter` reproduces the source's
  count through a crate-private `Protease::source_token_cleavage_count`.
- Fully specific `digest` discards one empty product and, with missed cleavages
  allowed, repeats the first product as a missed-cleavage product (`DFPIANGER`,
  Asp-N, one missed cleavage: `0-9, 0-9` in the source). The port returns the
  source's list with the repeats removed, in the source's order.
- Semi-specific `digest` also loses products the source's own `isValidProduct`
  accepts (`DFPIANGER`, Asp-N, no missed cleavage: `1-9` ... `8-9` are missing).
  The port returns exactly the ranges `isValidProduct` accepts, as it does for
  every other input.
- `peptideCount` counts the discarded empty product and the repeats.

`isValidProduct` is not affected (it counts only sites strictly inside a
product), and the port matches all 17,376 of its oracle answers. Without a
leading site the port matches every oracle record (`countInternalCleavageSites`,
both digest paths in order, fully specific `peptideCount`).
`tests/digestion_leading_site.rs` asserts all of this.

### Source arithmetic

The native operations above keep each boundary once and count by the
generation settings. The source's `ProteaseDigestion::digest` and `peptideCount`
loop over `tokenize_`'s positions as Boost's token iterator yields them, with
`Size` arithmetic, so their results differ from the native ones in the ways
listed above and in these:

- `unspecific cleavage` has the missed-cleavage ceiling `numeric_limits<Size>::max()`
  (`ProteaseDigestion.cpp:106`), so `digest` emits every product of one length
  before the next (by length, then position); the native digest orders by start,
  then length, as `EnzymaticDigestion::digestUnmodified` does.
- `peptideCount` is `n (n + 1) / 2` for `unspecific cleavage`, and otherwise the
  token count plus the missed-cleavage levels, without the length filter and
  without semi-specific products.
- `digest` with specificity `none` throws `Exception::InvalidValue`
  (`ProteaseDigestion.cpp:93-98`); a semi-specific digest of the empty sequence
  throws it too (`EnzymaticDigestion.cpp:121-126`).
- `countInternalCleavageSites` of the empty sequence is `0 - 1` in `Size`
  arithmetic, every enzyme but `no cleavage` (whose expression `()` is not
  split and yields position 0).
- `isValidProduct_` (`EnzymaticDigestion.cpp:211-303`) applies the N-terminal
  methionine allowance before it looks at the length, so a position of 1 or 2
  with length 0 on a protein that starts with `M` is the fragment `0..position`;
  the native range API answers `false` for an empty range. It accepts a negative
  length, tokenises only the fragment's window (`tokenize_(sequence, pos - 2,
  end + 2)`, so a lookbehind cannot see the residue before the window), and
  writes three `OPENMS_LOG_WARN` lines (start beyond the sequence, end beyond
  it, empty fragment or sequence), which `is_valid_product_source` returns in
  `SourceValidity::warnings`. Where the source is undefined it is refused with
  an error: the random D|P allowance with a start or end before the sequence
  (a read of `sequence[-1]`) and, for specificity none without the
  missed-cleavage shortcut, a negative length whose window ends before it
  starts (the Release build crashes with SIGSEGV).
- The semi-specific scan (`semiSpecificDigestion_`) computes a length with
  `Size` arithmetic, so a candidate whose end precedes its start wraps and is
  discarded. With a leading zero-width site the forward cursor is parked on the
  duplicated position 0 for the whole scan and every start `1 .. n-1` is paired
  with it, which is how the source loses the products its own `isValidProduct`
  accepts.

`digest_ranges_source`, `digest_unmodified_source`, `peptide_count_source`,
`count_internal_cleavage_sites_source` and `is_valid_product_source` port exactly
these loops over the
crate-private `Protease::source_tokenize`, so a caller that must reproduce the
Release build's output, defects included, selects them by name.
`digest_unmodified_source` is `digestUnmodified`'s pair overload, which
differs from `digest` in three things it does before tokenising: it returns
nothing, and counts nothing, for a sequence shorter than `min_length`; it
enumerates every substring ordered by start, then length, for `unspecific
cleavage` and specificity `none`; and it digests a specificity that is none of
full, semi and none (`SPEC_UNKNOWN`, `SPEC_NOCTERM`, `SPEC_NONTERM`) as full
where `digest` refuses it. Evidence: tier 1,
`../oracle/integ-B-digest-source` on ibminode06 (Linux x86_64 Release install
`openms4-release-bc9cc12-c19e494-174b576`, conda g++ 14.4.0, `-O3 -DNDEBUG`,
every enzyme in its own process, two byte-identical runs), 14,528 records over
eight enzymes (Asp-N, Asp-N_ambic, Lys-N, Trypsin, Trypsin/P, Chymotrypsin,
`unspecific cleavage`, `no cleavage`), sixteen sequences including the empty
one, specificities full, semi and none, 0, 1, 2 and 5 missed cleavages and four
length filters, retained as `tests/data/digestion_source_oracle.tsv` and
asserted record by record by `tests/digestion_source_tokenization.rs`:
`digest` (6,144 records, products, order and discarded count, and the refusals),
`digestUnmodified`'s pair overload (6,144 records, as `(start, length)`),
`peptideCount` (1,536), `countInternalCleavageSites` (128) and the
`isValidProduct` grid (576 lines, `pos` and `length` from -1 to n + 1 for two
sequences, with and without the missed-cleavage check, the N-terminal allowance
and the D|P allowance, each cell the answer and the warning kind). The specificity
enumeration and `filterByMissedCleavages` are tier 3 (`EnzymaticDigestion.h:41-49`
and `EnzymaticDigestion.cpp:20, 205-209`; the `Size` to `Int` conversion of the
empty sequence's count is the implementation-defined modulo wrap of the Linux
x86_64 build and was not executed).

## Resource limits and remaining capabilities

Input sequences are limited to 1,000,000 residues. Defaults permit at most
100,000 accepted products and 10,000,000 sequence/candidate work units; both
limits are configurable through `max_products` and `max_work`. Length-filtered
candidate attempts still consume work. Modified peptide construction additionally
limits the total copied residues to 10,000,000, checked before any peptide output
is constructed. Range and borrowed-slice APIs avoid those sequence copies.
Counts obey the same work/product limits as range generation.

Custom runtime enzyme insertion, provider/XML overrides (the public
`ProteaseDB(providers)` constructor), arbitrary regular expressions, `no-cterm`/`no-nterm` specificity modes and the C++ discarded-by-length
counter are not exposed. Registry gains and search-engine IDs are read-only;
empty gain formulas and absent IDs remain empty/`None` exactly as declared.
The library performs no runtime XML or filesystem lookup.

## Provenance and verification

Source SHA-256 values at revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`:

| Source | SHA-256 |
| --- | --- |
| `share/OpenMS/CHEMISTRY/Enzymes.xml` | `2f160f3ec32db6257cb4eee43fd594b48cb36398a17b7af21bbba76bc8b16995` |
| `src/openms/source/CHEMISTRY/BuiltInProteaseDataProvider.cpp` | `97075b2157f09bb6cc7b86913483731a166a02cca2a794082ea9e970f85ffaca` |
| `src/openms/source/CHEMISTRY/ProteaseDigestion.cpp` | `2f9f20d2ac22d83cc1f91635b4355f0a6382f58f12bb1f6c0d3016cca1cc767a` |
| `src/openms/source/CHEMISTRY/EnzymaticDigestion.cpp` | `d9205b1fcb7bca227f3c941f76341f517809e14eebdfc82d82a7095a4121040b` |
| `src/openms/source/CHEMISTRY/DigestionEnzymeProtein.cpp` | `f7c838edd087b323d866dc64aae17b6d8c5cea39d0ef101a98e6319602ee6d9a` |

`tools/generate_enzymes.py` translates the pinned XML to the native immutable
table and compiled boundary predicates. It separately applies Python's regex
engine to every uppercase two-residue bond and every three-residue context.
The fixture stores a count and FNV-1a64 fingerprint of the 18,252 boolean results
for each enzyme. Tests compare all 602,316 native decisions against these compact
independent regex fingerprints. The maximum source lookbehind width is two
residues and lookahead width one, so these cover every internal-bond context.

Additional tests port full/semi/nonspecific products, counts and validity cases
from the source suite, compare generation with a separate exhaustive small-range
oracle, check metadata and aliases, preserve modifications, and exercise resource
and boundary errors. The original chemistry and modification tests still pass.

```text
python3 tools/generate_enzymes.py --check
cargo test --offline --test digestion --test chemistry --test modifications
cargo clippy --offline --test digestion -- -D warnings
```

Regeneration uses the unchanged bundled `resources/enzymes/Enzymes.xml`, Python's
standard library and `rustfmt`; no reference checkout or network is required.
Normal library use needs none of these runtime resources.
No C++ build or live C++/Rust differential run was performed.
