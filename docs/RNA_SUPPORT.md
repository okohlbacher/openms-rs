# RNA records and sequences

The native chemistry module implements `Ribonucleotide`, the complete pinned
`RibonucleotideDB`, its in-memory/TSV/MODOMICS JSON inputs, and `NASequence`.
Sequences own immutable nucleotide handles and retain custom chemistry after
the creating registry is dropped. This covers the source foundation APIs;
[RNase digestion](RNASE_SUPPORT.md), [RNA spectrum generation](RNA_SPECTRUM_SUPPORT.md)
and [modification enumeration](RNA_MODIFICATION_SUPPORT.md) build on this foundation.
The [identification sequence/provenance graph](IDENTIFICATION_GRAPH_SUPPORT.md)
stores owned RNA parents and oligos, processing history, parent positions and
coverage. RNase digestion can register its products atomically in that graph.
Observation match groups, parent groups and legacy conversion remain separate work; the graph's observations, compounds and observation matches are implemented.

```rust
use openms::chemistry::{NAFragmentType, NASequence};

let sequence = NASequence::parse("pA[C*]Gp")?;
let formula = sequence.formula(NAFragmentType::Full, -2)?;
let ion_mass = sequence.mono_mass(NAFragmentType::Full, -2)?;
let mz = ion_mass / 2.0;
assert_eq!(sequence.suffix(1)?.to_string(), "*Gp");
```

The [example](../examples/analyze_rna.rs) prints this sequence's charged formula,
mass/m/z, sulfur-aware suffix and coarse isotope probabilities. The
[workflows](../tests/rna_workflow.rs) check independent elemental sums, custom
carbon-13 labeling, negative ion coordinates and both mzML compression modes.

## Records and registries

Construct an editable `RibonucleotideRecord`, then freeze it using
`Ribonucleotide::from_record`. Its ten independent fields are name, code, new
code, HTML code, formula, origin, declared monoisotopic mass, declared average
mass, terminal specificity and base-loss formula. Changing an input formula
does not recalculate either declared mass. The source default record has code
and origin `.`, name `unknown ribonucleotide`, empty formula, zero masses,
Anywhere specificity and base-loss formula `C5H10O5`.

`is_modified` compares the one-character code with its origin; `is_ambiguous`
tests for a final `?`. The latter deliberately does not recognize `?*`, although
the JSON provider's alternative-code parser does. Finite signed stored masses
are accepted. Nonfinite masses and empty codes are rejected; source operations
on an empty code would otherwise dereference its final character.

Record equality, ordering and hashing use complete values, including both
formulas and independent masses. Positive and negative zero compare and hash
equally. A shared code is not sufficient to establish chemical identity.

`RibonucleotideDB::global()` retains all 378 pinned records, loaded in the
source order: 333 accepted JSON records followed by 45 OpenMS custom records.
There are 375 distinct codes. `entries()` preserves duplicates and order;
`get(code)` returns the last record with that code. `get_prefix(text)` returns
the longest complete code prefix, with UTF-8-safe boundaries. It is independent
of sequence parsing, which does not use longest-prefix matching.

Use `from_records` for ordinary owned records or `from_entries` for entries
with alternative codes. The optional `[String; 2]` in `RibonucleotideEntry`
corresponds to the source's two alternatives; the first code, if nonempty, activates
deferred resolution. All entries are indexed before alternatives are resolved.
They consequently point to the final code winners. A later successful ambiguity
overwrites an earlier mapping; a later unresolved one leaves the earlier valid
mapping and adds a diagnostic. `alternatives(code)` returns the two handles.
The registry is immutable after construction; cloning shares the frozen records.

The complete [data provenance and notices](../resources/rna/README.md) include
the original files, independently checked projection and regeneration command.
The default registry requires no runtime files or JSON parser. MODOMICS data
has separate, unresolved redistribution terms for release; the software's BSD
license is not applied to it. The current package is local and unpublished.

## Input providers

`ribonucleotide_db::read_tsv(&str)` and, with the `rna-json` feature,
`read_modomics_json(&str)` return a `RibonucleotideLoadReport`. Its entries can be
combined in order and passed to `from_entries`. Ordinary invalid rows are skipped
with indexed diagnostics, matching the source's logging-and-continuing behavior.
Invalid whole documents, headers and configured resource-limit violations are
fatal. No provider discovers files or downloads data implicitly.

The JSON reader iterates object keys lexicographically; arrays retain their
input order. Required fields, moiety shape and optional value types follow the
source. One one-character moiety sets the origin; four moieties select X and
the source terminal-code rules. It accepts the first two alternatives for codes
ending `?` or `?*`. Empty codes are omitted. Explicit base-loss formulas override
code heuristics.

The TSV reader accepts the original nine-column header and optional additional
columns, skips leading `#` comments, and replaces Unicode PRIME with ASCII
apostrophe throughout each data row. It retains source `QtRNA` code shortening,
`preQ0base` origin mapping, terminal-code priority and base-loss rules. An
ambiguity column uses the text before the first and after the last space.
Terminal new-code branches bypass base-loss and ambiguity heuristics.

The two providers intentionally differ in missing-mass handling:

| Input field | JSON | TSV |
| --- | --- | --- |
| Missing mono mass (JSON: absent/null; TSV: empty/`None`) | Formula mono mass | Zero |
| Explicit mono zero | Retain zero | Formula mono mass when formula is nonempty |
| Missing average mass (JSON: absent/null; TSV: empty/`None`) | Zero | Zero |
| Explicit average zero | Retain zero | Formula average mass when formula is nonempty |

The base-loss default is `C5H10O5`; deoxy codes use `C5H10O4`, and ribose methyl
codes use `C6H12O5`. The literal unusual `C10H19O21P` formula for Ar(p)/Gr(p)
is retained. JSON recognizes those suffixes; TSV requires exact codes.

## Sequence text, ownership and slicing

`parse` uses the embedded registry; `parse_with_registry` uses a supplied one.
`from_records` and `from_records_with_registry` retain supplied `Arc` handles.
The second form also captures that registry's optional `5'-p*` record for later
sulfur-aware slices. `residues`, checked `get_residue`/`set_residue`,
`set_sequence`, end getters/setters, `clear`, `len` and `is_empty` replace the
C++ pointer-container accessors.

The source syntax is preserved:

- A first `p` or `*` sets the five-prime phosphate or phosphorothioate. With
  input length greater than one, a final `p` or `c` sets the three-prime phosphate
  or cyclic phosphate. Interior ASCII spaces are skipped; other whitespace is
  not normalized and case is significant.
- Bare ASCII characters resolve individually and become sequence residues,
  regardless of their record's specificity. Non-ASCII custom codes need brackets.
- Brackets resolve exact complete codes. A bracketed terminal record sets that
  end wherever it appears; a later one overwrites it. Typed end setters likewise
  retain explicitly supplied placements without applying extra specificity rules.
- Display abbreviates exact terminal codes `5'-p`, `5'-p*`, `3'-p`, `3'-c`
  as `p`, `*`, `p`, `c`; other ends and multibyte residue codes use brackets.

Source display can lose unusual custom placement or same-code chemistry.
Even a residue-free three-prime phosphate or cyclic phosphate displays as `p`
or `c`, which does not reconstruct that same terminal-only sequence.
`checked_string_with_registry` verifies full reconstruction in the chosen
registry and returns an error on information loss. Ordinary `Display` still
produces the source spelling.

`prefix(length)` and `suffix(length)` reject lengths greater than or equal to
sequence length. Zero length is valid on a nonempty sequence and can retain
end modifications. `subsequence(start, optional_length)` rejects `start >= len`,
including zero-length requests there; excessive length is clamped. Prefixes
retain only the five-prime end, suffixes only the three-prime end, and a
subsequence retains original ends only when touching the corresponding boundary.

A slice immediately after a residue code ending `*` installs the captured
`5'-p*` record. This also applies to `suffix(0)` after a final starred residue.
Slicing never consults a hidden global registry. A missing required handle is a
checked error; `with_phosphorothioate_end` can supply or remove that context.
Complete sequence identity includes this context only when a starred residue
can make it relevant. Ordinary and empty sequences do not compare unequal
because of unused registry context. Pointer-dependent C++ ordering is replaced
with consistent value ordering.

## Formulas and mass conventions

Each RNA record describes a complete nucleoside. Between successive records,
add `H-1PO2`; after a code ending `*`, add `H-1POS` instead. A final starred
record alone adds no extra sulfur linkage. Let the resulting sum be B, and F/T
be the five/three-prime terminal formula minus H, or empty when absent.
For supported fragments, add `q` natural hydrogen atoms to B, then:

| Fragment | Additional formula |
| --- | --- |
| Full | F + T |
| a / b | F - H2O / F |
| c / d | F + H-1PO2 / F + HPO3; add SO-1 when the final residue is starred |
| w / x | T + HPO3 / T + H-1PO2; add SO-1 when F equals HPO2S |
| y / z | T / T - H2O |
| a-B | F - H4O2 - last residue formula + last base-loss formula |

All 20 finite `NAFragmentType` variants are represented. Ten legacy variants
(Internal, FivePrime, ThreePrime, Precursor, b/y minus water or ammonia,
NonIdentified and Unannotated) retain the source fallback: B without end groups
or charge hydrogens. They are not silently treated as Full.

`mono_mass(fragment, q)` and `average_mass(fragment, q)` evaluate that formula
then subtract `q * ELECTRON_MASS_U`. They return **ion mass, not m/z**: callers
divide by charge magnitude when appropriate. Natural hydrogen addition and
electron subtraction preserve source isotope and rounding conventions; using
`EmpiricalFormula::with_charge` would substitute proton mass semantics.

An empty sequence returns an empty formula before considering fragment, charge
or ends. Its charged mass is therefore `-q * ELECTRON_MASS_U`. Empty source
residue formulas, including N, remain empty; a nonempty chain of such records
can contain linkage chemistry alone. Signed formulas and finite negative masses
are retained, while atom-count overflow and nonfinite results are checked.

## Resource bounds and evidence

Records allow 4,096 code bytes and 65,536 bytes in each other text field.
Registry construction allows 100,000 entries, 128 MiB of conservative logical
storage and 50 million charged key/work operations. Prefix lookups also have a
work allowance.

Providers allow 16 MiB of input, 100,000 rows and 128 MiB of accounted output and
diagnostics. TSV rows are limited to 1 MiB. JSON is checked before deserialization
for depth 64, 250,000 structural tokens and 393,216 raw bytes per string.
Consumed formula fields and cumulative formula work are bounded separately.

Sequences allow one million residues, 16 MiB of input/display text and 256 MiB
of conservative logical payload, counting complete record payload per occurrence.
Formula operations share 50 million work units and 256 MiB of cumulative map
allocation accounting. Bounds are checked before relevant copying or formula
growth; sequence replacement failures preserve previous values.

The [independent review](RNA_REFERENCE_REVIEW.md) records original formula,
mono/average mass and slice assertions, every registry field and additional
boundary cases. [Validation](VALIDATION.md) states the executed feature/toolchain
checks. No C++ reference was built or executed, and these tests do not establish
full RNA search or downstream spectrum-generator parity.
