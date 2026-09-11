# ProForma / AASequence conversion

`Peptidoform` provides the complete pinned source AASequence conversion group.
The [ProForma header](PROFORMA_SUPPORT.md) now has native equivalents for every
public operation group, including separate [spectrum generation](PROFORMA_SPECTRA_SUPPORT.md).

| Operation | Native result |
|---|---|
| `aa_sequence_conversion_issues(&mut registry)` | `ConversionEvaluation<Vec<ConversionIssue>>` |
| `is_representable_as_aa_sequence(&mut registry)` | `ConversionEvaluation<bool>` |
| `to_aa_sequence(&mut registry)` | `ConversionEvaluation<AASequence>` using `FailOnLoss` |
| `to_aa_sequence_with_policy(policy, &mut registry)` | `ConversionEvaluation<AASequence>` |
| `Peptidoform::from_aa_sequence(&sequence)` | `Peptidoform` |

Every return is inside `Result`. Evaluation results contain `value` and owned
`warnings`. `ConversionWarning` distinguishes existing resolution warnings from
summed-formula versus declared-mass disagreements. Input ASTs and sequences stay
unchanged. Source has no ion overload for this group, so none is invented here.

```rust
use openms::chemistry::{ModificationsDB, proforma::{Peptidoform, WriteMode}};

let mut registry = ModificationsDB::global().clone();
let annotation = Peptidoform::parse("PEPM[Formula:O]TIDE")?;
let sequence = annotation.to_aa_sequence(&mut registry)?.value;
let restored = Peptidoform::from_aa_sequence(&sequence)?;
assert!(restored.to_text(WriteMode::Lossless)?.contains("Formula:"));
Ok::<(), Box<dyn std::error::Error>>(())
```

## Policies and diagnostics

The default `FailOnLoss` policy rejects reported conversion losses. Both
`DropUnlocalised` and `BestEffort` execute the same branches in the pinned C++
implementation; the native port preserves this behavior. They omit unlocalised,
labile and global modifications, keep the first ambiguous candidate, expand
range residue letters, and omit modifications inside ambiguous/range sections.
Global isotope replacements are ignored even by strict conversion and its issue
predicate. Neither policy synthesizes chemistry for unmatched mass deltas.

Source issue order is unlocalised, labile, first global-modification presence,
sequence sections, N-terminal brackets, then C-terminal brackets. A position
counts sequence sections rather than the flattened output residues. `None`
replaces source `SIZE_MAX` for unknown/terminal positions. Ordinary brackets
report multiple chemistry brackets, unresolved chemistry, genuine chemistry
alternatives and crosslink labels. INFO and position annotations are not
chemistry. Unknown residue letters are not checked by this predicate, although
AASequence construction can reject them. Source warning text and scientific
conditions are represented as typed owned data instead of global log output.

The [C++ issue ledger](../OpenMS_CPP_ISSUES.md) records preserved inconsistencies:

- **CPP-029:** `M[INFO:note]` has no conversion issue and passes the predicate,
  but strict attachment throws because its resolved handle is absent. A
  permissive conversion succeeds and drops the annotation.
- **CPP-030:** terminal crosslink labels are not checked by the predicate and
  disappear on conversion. A direct AST regression pins this behavior.
- **CPP-031:** an empty manually constructed ambiguous region emits no residue,
  but advances the attachment cursor. A subsequent modified M can therefore
  cause the source index error. Native index checks remain safe and occur at
  the corresponding attachment, before later sections perform more work.

These are source-reviewed reproductions, not executed C++ conversion tests.
Strict conversion also reports modified/ambiguous sections as unsupported rather
than attempting a standards-level interpretation of them.

## Modification selection and combination

The [explicit resolver](PROFORMA_RESOLUTION_SUPPORT.md) runs once on a bounded
copy. Lookup keeps its native deterministic first-provider policy in place of
allocation-dependent C++ pointer ordering. A single resolved residue handle is
attached directly, preserving owned record identity without a second name
lookup. Multiple handles are gathered in bracket order, including repeated
references to the same allocation. Empty alternatives retain a preexisting
handle; those handles can participate in combination without being counted as
chemistry by source diagnostics.

Combination sums declared diff mono masses in order. A sum with absolute value
at most `1e-6` leaves the residue unmodified. Nonempty diff formulas are also
summed. If every component supplied a formula, the sum is nonempty and its mass
agrees within `1e-3` Da, the canonical formula/site is interned. Otherwise an
anonymous mass-only record is interned using the summed declared mass. A
formula/mass disagreement warning is produced only when every component had a
formula. Absolute masses are populated from the source free-residue base.
Existing exact intern keys retain their previous record, as in the source.

Formula canonical text omits charge. Consequently, a charged summed formula can
pass the agreement check and then be reparsed and interned as neutral chemistry.
This scientific loss is tracked as **CPP-037**: two pre-resolved O formulas,
one carrying charge +1, combine into neutral O2 and lose one proton mass. Both
the output formula and the mass difference are explicitly tested; no C++
conversion execution is claimed. Source terminal attachment chooses
the first resolved bracket, including when a label-only bracket precedes it;
multiple terminal chemistry brackets are diagnosed rather than combined.

The returned object obeys existing [AASequence rules](SEQUENCE_SUPPORT.md):
unknown B/Z/X composition remains unknown unless its own supported absolute
chemistry supplies it, and cached construction requires nonnegative finite
calculated mass. Direct [ProForma mass](PROFORMA_MASS_SUPPORT.md) separately
allows finite signed results. These inherited native representation boundaries
are not changes to the underlying AASequence API.

## Reverse conversion and text precision

Reverse conversion visits residues then N/C termini and retains every known
modification's `Arc` in the resulting AST. Tag precedence follows source:

1. A UniMod accession produces a CV tag.
2. An anonymous record with a diff formula produces a Formula tag.
3. A named `Defined` record produces formula or mass chemistry followed by its
   name in INFO.
4. Other named records produce a NamedMod.
5. A remaining anonymous record produces a mass delta.

The source formula string omits charge. Fields not filled by the source retain
AST defaults. No registry is read or changed by reverse conversion.

Native anonymous `MassTag` values receive a detached owned anonymous record in
the returned AST. Their finite declared delta is retained. An absolute B/Z/X tag
has no delta in the existing native API; conversion locally reconstructs the
source zero-free-mass convention as **internal target mass + H2O**. This preserves
its explicit mass without changing the original tag's unknown delta/formula
accessors. Neutral H/OH terminal conventions remain unchanged.

Generated mass spelling is signed fixed decimal at round-trip precision;
both signed zeros render `+0`. Integral doubles require exact fixed zero
precision, because shortest general formatting can round their low integer
digits. For example, `100000000000000016384` must not become
`100000000000000020000`. Fractional values use shortest fixed decimal, including
subnormals. The existing canonical writer can independently round these tags;
lossless writing preserves their generated text.

## Transaction and limits

All forward operations use one private resolver session and shared budget.
A completed issue/predicate report publishes formula interning even when issues
exist. Successful conversion publishes new formula/combined records. Any `Err`
rolls back registry changes and publishes no sequence. This atomic behavior is
an explicit native correction to source singleton mutations surviving exceptions.

Per operation limits are **50 million work units**, **1 million consumed items**,
**256 MiB cumulative logical allocation allowance**, and **4 MiB per copied or
consumed text value**, including aggregated failure diagnostics. AST copies,
record traversal, formula sums, sparse map overhead, registry work, numeric
formatting, diagnostics, vector moves and destruction are precharged. Existing
AASequence rebuilding is conservatively charged against all potentially present
formula keys at every residue before its helper runs. A limit can therefore
reject work earlier than an implementation charging actual allocations; these
are logical allowances rather than physical allocator guarantees. Numeric or
count overflow is a checked error and is never swallowed as an unresolved name.

## Evidence

[Thirteen direct tests](../tests/proforma_conversion.rs) cover source examples,
all policies, diagnostics, source defects, identity, combinations, formula
round trips, anonymous and defined records, precision and transaction boundaries.
Four private tests force late resource failures after interning, cumulative
residue work, preallocation failure and fixed-decimal boundary behavior.
Adjacent mass/resolver/sequence/parser/writer/JSON tests remain separate evidence.

A [small reproducible probe](../tools/generate_proforma_conversion_mass_probe.py)
extracts the exact pinned `massDeltaText_` source helper and executes it with a
minimal exception stand-in. Its [1,078-row fixture](../tests/data/proforma_conversion_mass_text.tsv)
contains both signs of explicit boundaries and deterministic generated finite
binary64 values. Every row is checked through the complete native reverse API,
lossless writing and parsing. These are executed C++ **formatter-only** results;
no full SDK or C++ conversion algorithm was built or run. Source literals and
independent chemical tests remain labeled separately in the
[provenance manifest](../tests/data/proforma_conversion_provenance.json).
