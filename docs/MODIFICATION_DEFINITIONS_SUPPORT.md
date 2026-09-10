# Modification definitions

`chemistry::{ModificationDefinition, ModificationDefinitionsSet}` ports the
search-definition models, compatibility predicate, mass matching, and inference
from `CHEMISTRY/ModificationDefinition.{h,cpp}` and
`CHEMISTRY/ModificationDefinitionsSet.{h,cpp}` at OpenMS4-core revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. The implementation and tests are
BSD-3-Clause; bundled chemical data retains its notices under
`resources/modifications`.

## Native interface

```rust
use openms::chemistry::{
    AASequence, ModificationDefinitionsSet, ModificationMatchOptions,
};

let definitions = ModificationDefinitionsSet::from_names(
    &["Carbamidomethyl (C)"], &["Oxidation (M)"],
)?;
assert!(definitions.is_compatible(&AASequence::parse("AC(Carbamidomethyl)M")?)?);
let matches = definitions.find_matches(15.994915, &ModificationMatchOptions::default())?;
assert_eq!(matches[0].definition.modification_name(), "Oxidation (M)");
# Ok::<(), openms::Error>(())
```

A definition owns a validated `ResidueModification` copy or an anonymous
`MassTag`. `new(name)` defaults to fixed and zero maximum occurrences;
`with_options` supplies both settings.
`from_modification` also accepts a record from a custom native registry, without
borrowing that registry's lifetime. `Default` is unset, fixed, and count zero.
Its name is empty and its chemistry getter returns an error. An unset definition
cannot enter a set. `set_modification` resolves before replacing; `fixed` and
`max_occurrences` are mutable public fields on standalone definitions.

`from_annotation` accepts either kind of `SequenceModification`, and
`set_annotation` replaces that owned chemistry while retaining the search
settings. `mass_tag()` exposes an anonymous annotation, including exact spelling,
origin, specificity, and independently known masses. The existing `modification()`
getter returns a named record, or `Unsupported` for a mass tag. Name-based
constructors continue to resolve the immutable global registry: retaining an
anonymous annotation never registers a name or changes subsequent parsing.

The set stores two ordered maps, one per fixed/variable partition. An insertion
keeps the first definition with a particular full ID in that partition. The
same ID may appear in both partitions; `len` counts both. `fixed_modifications`
and `variable_modifications` expose immutable ordered iterators, and their name
getters return ordered string sets. `modifications` merges in full-ID order,
preferring a fixed definition when both partitions contain the same ID.

`set_names` and its comma-separated counterpart replace both partitions
atomically. Comma-separated tokens are not trimmed, matching source `ListUtils`.
`set_modifications` accepts a slice corresponding to the source's single input
set: the first full ID wins globally, even across different fixed flags. In
contrast, `add_modification` deduplicates within its selected partition.

Definition equality compares chemical contents, fixed status, and occurrence
count. This replaces C++ pointer identity with owned value equality. Hashing is
consistent with equality, including the unset default. Definitions intentionally
have no `Ord`, since the source's full-ID-only ordering disagrees with equality.
Set equality includes both partitions and `max_modifications`, but excludes its
native resource budget. Mutating a returned reference cannot invalidate a map
key or move an entry between partitions.

## Compatibility and inference

`is_compatible` preserves the actual source predicate:

- Each fixed definition requires every residue with that origin letter to carry
  a residue modification with the same short modification name.
- Every present residue or terminal modification must have a full ID in either
  partition.
- Stored `max_occurrences` and `max_modifications` are **not enforced**.

The first rule ignores terminal specificity. A generic fixed N-terminal Acetyl
has wildcard origin `X` and does not require an ordinary peptide's terminal slot
to be acetylated; it does require any actual `X` residue to be modified. A fixed
`Gln->pyro-Glu (N-term Q)` checks every Q residue annotation, so it can reject a
peptide carrying that modification solely in its N-terminal slot. These are
preserved source semantics, not a stronger fixed-site validation policy.

Source anonymous records have an empty short ID. The fixed-residue check thus
accepts any anonymous annotation at that origin before the separate full-ID
membership check. For example, a set containing fixed `K[+12.345678901]` and
variable `K[+13.345678901]` accepts either listed tag at K. An unmodified K or an
unlisted tag still fails. The numeric spelling is not used as a short name.

`infer_from_peptides` reads **all hits** in every input identification, pooling
modification observations separately by residue letter, N terminus, and C
terminus. A site category with exactly one nonempty modification and no
unmodified observation produces a fixed definition. Otherwise its observed
modifications become variable. Residues absent from a peptide add no observation;
every hit, including an empty sequence, contributes both terminal observations.
Inference resets per-definition occurrence counts to zero and preserves the
set's `max_modifications`. It consults neither ranks, scores, nor protein evidence.
Consequently, an inferred fixed terminal definition can still encounter the
source compatibility inconsistency above.

Unresolved B/Z/X residues and anonymous mass annotations participate in inference
without requiring unknown chemistry. Anonymous results retain their owned tags,
including exact full IDs and decimal spelling. Equal numerical values written
differently remain distinct full IDs, as in the source's dynamic names. Numeric
syntax that already resolved to a known registry entry remains named.

Inference uses native chemical value identity, so separate allocations of the
same chemical record count as one observation, unlike the source's pointer
identity. If one full ID appears with conflicting chemical contents (possible
through a generator supplied with custom owned records), inference returns an
atomic error. It does not pick one mass based on hit or allocator order.

## Monoisotopic mass matching

`ModificationMatchOptions` defaults to both partitions, no residue or terminal
filter, delta-mass mode, and an inclusive 0.01 Da tolerance. Returned matches
are ordered by increasing absolute error. Equal errors retain fixed entries
before variable entries and full-ID order within each partition. An ID present
in both partitions may therefore produce two matches.

Origin filtering uses the **first character of the supplied residue string**,
as in the source; it does not first resolve a full residue name. Empty, `.` and
`X` selectors, and wildcard-origin definitions, bypass origin filtering.
Terminal filtering is exact when supplied.

Delta mode uses the record's declared monoisotopic mass difference. Absolute
mode compares a positive declared absolute mass directly. If that value is
nonpositive and a residue was supplied, it uses:

```text
declared difference + (unmodified full-residue formula mass − water formula mass)
```

The parentheses and full-formula calculation retain source arithmetic order.
Residue aliases and their case sensitivity follow the pinned `ResidueDB`.
Without a residue selector, the stored absolute value is compared literally;
the record's origin is not substituted as a missing residue.

Bundled UniMod/custom records, and records parsed by the existing TSV loader,
have stored absolute masses zero. `ResidueModification::with_absolute_masses`
adds finite declared mono/average masses to an owned copy without changing its
formula, deltas, global registry, or file schema. Nonpositive finite values keep
the source's fallback convention. This exercises the genuine positive stored
mass branch without requiring an OBO record. The shared registry also supports
bounded caller-owned OBO loading; definitions can own those records directly.

Anonymous residue tags have a different absolute-mass anchor from the numeric
value displayed inside a peptide. Source `createUnknownFromMassString` stores
the **full modified residue mass**. A signed tag uses its delta plus the full
unmodified residue formula mass. An absolute internal tag first subtracts
`full residue mass − water`, then adds the full residue mass back. These
operation orders are retained explicitly. Consequently, definition delta
matching for an absolute tag can differ in its final floating-point bits from
the tag's cached native delta, whose base is the internal residue formula.
Neither cached tag values nor the peptide's chemistry are changed by matching.

N-terminal tags use the H mass anchor and C-terminal tags use OH. Signed tags
add that anchor to obtain their stored absolute mass; absolute terminal tags
use the written value directly. These source values then follow the same
positive-value/fallback rule as named records.

Absolute B/Z/X tags retain a known internal mass and therefore have a full
modified residue mass equal to that mass plus water. They do not establish an
unmodified residue mass: delta matching returns `Unsupported` if such a record
passes the origin/terminal filters. Matching does not silently omit an eligible
record with unavailable chemistry or manufacture a delta from an empty formula.

Fallback on B/Z/X is `Unsupported`: the C++ empty-formula placeholders would
otherwise manufacture a negative internal residue mass. A positive stored
absolute mass needs no fallback and can therefore be compared with those
selectors. Unknown aliases (including `.`) error only when fallback is reached;
prior origin/terminal filters and positive stored masses can bypass lookup.
Non-ASCII residue selectors, nonfinite inputs, negative tolerances, disabled
fixed-and-variable selection, and arithmetic overflow are checked errors.

## Bounds and validation

`max_work` is a positive per-operation logical visit budget, default 1,000,000.
It is independent of the source's stored occurrence settings. Replacement
preflights input definition/name counts and identifier bytes before cloning or
database lookup; comma input also checks its bytes before splitting. Insertion
charges current set entries and the incoming full ID bytes. Compatibility
charges set entries, peptide sites and two termini, plus fixed-count times
peptide-length. Inference charges input identification count, hit count, and
each sequence length plus two termini and each present annotation's full-ID
bytes before allocating observation maps or cloning anonymous spellings.
Mass matching charges the complete stored entry count and selector bytes.

These are logical bounds, not exact CPU instruction counts. Ordered maps add
logarithmic factors; sorting at most `m` mass matches costs O(m log m). Result
record counts are bounded by stored definitions or observed sites; this is not
a byte/RSS allowance for custom chemical records. Standalone getters operate
on already-owned data. Immutable global registry initialization and
irrelevant caller-owned identification metadata are not charged. All fallible
set mutations construct a replacement before assignment.

`tests/modification_definitions.rs` covers source constructor/count/partition
expectations, all seven compatibility examples, pyro-Glu mass-match ordering,
the source inference example, owned stored masses and fallback arithmetic,
terminal quirks, anonymous annotations, equality/hash behavior, and atomic
resource/error paths. The independent reference suite and provenance are
maintained separately from these implementation tests. No C++ build or newly
executed C++ numerical result is claimed.

`tests/anonymous_modification_definitions.rs` adds source-derived residue and
terminal mass-anchor calculations, empty-short-ID compatibility, exact-spelling
inference, post-drop ownership, B/Z/X partial chemistry, and a checked inference
payload boundary. Existing custom-record collision tests use owned `Arc` handles
and require no leaked allocations.
