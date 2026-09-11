# ProForma modification resolution

`Peptidoform::resolve_modifications(&mut self, &mut ModificationsDB)` implements
the complete source `ProForma::resolveModifications(Peptidoform&)` operation. It
returns `Result<Vec<ResolutionWarning>>`. The supplied registry owns any new
formula definitions; the peptidoform receives immutable `Arc` handles. Neither
the global registry nor a global logger is modified. There is no source ion
overload, and no additional ion-resolution API is introduced.

```rust
use openms::chemistry::{ModificationsDB, proforma::Peptidoform};

let mut registry = ModificationsDB::global().clone();
let mut peptide = Peptidoform::parse("EM[UNIMOD:35]K")?;
let warnings = peptide.resolve_modifications(&mut registry)?;
assert!(warnings.is_empty());
Ok::<(), Box<dyn std::error::Error>>(())
```

This operation fills handles; it does not rewrite annotation text, combine
brackets, convert to AASequence, calculate masses, validate cross-links or
generate spectra. [Mass/mz operations](PROFORMA_MASS_SUPPORT.md) use resolution
through their own atomic transaction, as does [AASequence conversion](PROFORMA_CONVERSION_SUPPORT.md).
[Spectrum generation](PROFORMA_SPECTRA_SUPPORT.md) also uses an explicit atomic transaction.

## Traversal and alternatives

The traversal is literal source order: sequence elements; ambiguous-region
elements; modified-range annotations; N-terminal annotations; C-terminal
annotations; unlocalised groups; labile annotations; then global modifications.
**Annotations on the elements inside a modified range are not visited**, matching
the source implementation. Global isotope replacements and stored locations,
occurrences, names, labels and charges are not consulted. An empty alternative
list retains its old resolved handle.

For each nonempty bracket, resolution first looks for an INFO annotation naming
a registered `Defined` record. It then attempts only the **first chemistry
alternative**. INFO and position constraints are annotations; every other tag,
including an unresolved glycan, counts as chemistry for this selection. Failure
of that first alternative does not try subsequent chemistry alternatives.
An annotation-only bracket resolves to no handle.

A resolving Defined INFO name wins over inline chemistry. Inline resolution
still runs and can create a formula record that remains in the registry even
when its handle is overridden. If both candidates resolve, equal nonempty diff
formulas mean agreement regardless of their stored masses. If either formula
is empty, absolute mass difference at most `1e-6` means agreement. Disagreement
returns `ResolutionWarning::DefinitionDisagreement` with both owned handles.
INFO text naming an ordinary vocabulary entry has no such overriding effect.

## Lookup behavior

All five CV prefixes are supported. An accession on a residue first searches
`Anywhere`, then falls back to unrestricted terminal specificity. A named tag
directly uses the supplied specificity; its CV hint is ignored. Lookup tests
the exact original alias before source UniMod-prefix normalization. The INFO
Defined gate uses the exact alias only. A narrow private index accessor retains
these rules without changing the ordinary native registry APIs.

The C++ name search selects by allocation-dependent pointer-set order. This
resolver deterministically selects the **first matching native provider entry**.
CV ambiguity returns an `AmbiguousAccession` warning while retaining that entry;
named lookup ignores the source ambiguity flag. It does not turn ambiguity into
the checked error used by `ModificationsDB::get_modification_handle`.
Missing alias keys return `ModificationNotFound` warnings. An existing alias
with incompatible residue/term simply stays unresolved. CV retry can produce
two missing-key warnings, as in source. Empty aliases omitted by the ordinary
native index are matched against source empty short/full-name or accession fields.

Source residue queries `X`, `.`, `?`, and no residue match non-X origins.
A named origin-X is a wildcard; an anonymous origin-X is an actual X except for
an unspecified query. Native residue fields are Unicode scalars rather than raw
C++ bytes: non-ASCII scalars are unknown to the formula-residue table, while
ordinary name matching retains scalar equality and wildcard behavior. Complete
UTF-8 modification names, INFO text and diagnostics are preserved.

Mass deltas search stored diff-mono masses with strict error `< 0.01` Da,
retaining the first provider tie. Zero-difference records are skipped unless
the query itself is zero. The mass source hint and preserved spelling do not
participate. Nonfinite consumed mass deltas are checked errors. Glycans, INFO
and position constraints do not independently produce a chemistry record.

## Formula definitions and transaction

A nonzero explicit formula charge, failed formula parsing, empty formula or
nonzero parsed charge remains unresolved. Resource failures are separate errors
and are never swallowed as parse failures. The existing native empirical-formula
grammar, checked `i32` counts, element data and mass accumulation remain in use;
this group does not claim bytewise parity with pointer-ordered C++ mass sums.

Neutral formulas are interned by site and canonical formula, with explicit
counts sorted by complete element-symbol spelling, including isotope labels.
Examples are `M[Formula:O1]`, `.n[Formula:C2H2O1]`, and
`.c[Formula:C2H2O1]`. Equivalent spellings reuse the same registry handle.
An already existing key wins without revalidating its chemistry or site.
This formula identity is distinct from a mass-bracket key.

New definitions have an empty short ID, `MassOnly` provenance, a signed
mass-bracket full name, and explicit diff formula/mono/average masses. Residue
absolute masses add the source free-residue masses; B/Z/X have zero base masses.
Unknown or lowercase formula sites remain unresolved. Terminal absolute mono
mass adds H for N-term and OH for C-term; terminal absolute average mass retains
the source zero default. Signed finite formulas and masses remain valid.

The implementation stages handle replacements without copying the whole AST.
It lazily clones the registry only for a new formula insertion, so later tags
see earlier additions. Publication of the staged registry and all handles occurs
only after every operation succeeds. A late invalid value or resource failure
leaves both original objects unchanged. Existing and newly returned `Arc`
records survive their registry through ordinary ownership.

## Limits and evidence

Each call shares **1,000,000 consumed items**, **50,000,000 work units**,
**256 MiB cumulative logical allocation allowance**, and **4 MiB per consumed
text value**. Lookup comparisons, collection walks, formatting, formula maps,
registry copies/extensions, diagnostics and replaced-record destruction are
precharged. Registry bounds also apply. Estimates are deliberately conservative:
repeated formula insertions charge repeated index copies and can exhaust work
before a batched implementation would. These are logical bounds, not exact
allocator bookkeeping or physical-memory ceilings. Ignored AST fields are
neither validated nor copied.

[Sixteen direct tests](../tests/proforma_resolution.rs) preserve source lookup,
formula and INFO examples and independently cover custom aliases, all CV kinds,
specificity/ambiguity, wildcard behavior, mass boundaries, all traversed groups,
formula identity, Unicode, diagnostics and rollback. Four private tests check
late work/byte/item exhaustion, invalid-formula versus resource failure, failing
lookup accounting and ignored state. Existing parser/writer/JSON/registry tests
remain separate regression evidence.

The [source manifest](../tests/data/proforma_resolution_provenance.json) records
the exact source hashes and relevant ranges. Primary implementations are
[`ProForma.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L1590)
and [`ModificationsDB.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ModificationsDB.cpp#L137).
The [source class tests](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/tests/class_tests/openms/source/ProFormaParser_test.cpp#L1415)
include resolution-only portions of conversion examples; this group does not
claim to execute the separate conversion APIs or an upstream C++ resolver.
Historical text/JSON/probe manifests are preserved.
