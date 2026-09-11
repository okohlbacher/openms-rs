# ProForma mass and m/z

Both `Peptidoform` and `PeptidoformIon` implement the complete source mass-operation
group: issue collection, availability checking, monoisotopic mass, m/z, and the
optional mass/mz variants. Every operation takes a caller-owned
`&mut ModificationsDB`. The AST is borrowed and remains unchanged.

| Native method | Returned value inside `Result` |
|---|---|
| `mass_calculation_issues` | `MassEvaluation<Vec<ConversionIssue>>` |
| `can_calculate_mass` | `MassEvaluation<bool>` |
| `mono_mass` | `MassEvaluation<f64>` |
| `mz` | `MassEvaluation<f64>` |
| `try_mono_mass` | `MassAttempt` |
| `try_mz` | `MassAttempt` |

`MassEvaluation<T>` contains `value` and owned resolution `warnings`.
`MassAttempt` contains `value: Option<f64>`, `issues`, and `warnings`, replacing
both source optional overloads and their mutable issue-output parameter.
An `Ok` attempt with no value means scientific unavailability. Invalid numerical
state or a resource limit returns `Err`. Chain `mz`/`try_mz` receive an explicit
`i32` charge before the registry argument; ion methods use the ion's charge.

```rust
use openms::chemistry::{ModificationsDB, proforma::PeptidoformIon};

let mut registry = ModificationsDB::global().clone();
let ion = PeptidoformIon::parse("PEPTIDE/2")?;
let result = ion.try_mz(&mut registry)?;
assert!(result.issues.is_empty());
let mz = result.value.expect("known peptide and nonzero charge");
assert!((mz - 400.6872).abs() < 0.001);
Ok::<(), Box<dyn std::error::Error>>(())
```

The [ProForma header](PROFORMA_SUPPORT.md) also supplies
[ordinary/XLMS spectrum generation](PROFORMA_SPECTRA_SUPPORT.md).
[AASequence conversion](PROFORMA_CONVERSION_SUPPORT.md) is available separately. These operations preserve the pinned SDK's behavior, including known
C++ defects; they are not a standards-corrected interpretation of all notation.

## Calculation and lookup

Resolution uses the [existing explicit registry rules](PROFORMA_RESOLUTION_SUPPORT.md),
including first native provider order in place of allocation-dependent C++ pointer
order. Stored resolved chemistry takes priority over tag mass, even with an empty
alternative list. Otherwise the first chemistry-carrying tag is selected; INFO
and position constraints do not count as chemistry. Mass deltas use their explicit
value only when resolution did not supply a record. A nearest registry record
within the source strict 0.01-Da window can therefore replace a rounded delta.

Unresolved Formula tags parse their formula text and use its mass, including
parsed charge. This fallback ignores the separate `FormulaTag.charge` field and
accepts an empty formula with mass zero, unlike the resolver's interning rules.
Invalid formulas are unavailable; a precharged resource failure remains an error.
Annotation-only brackets contribute zero. Unresolved names, accessions and glycans
have no mass. The native formula grammar, element data and deterministic atom
summation are inherited; no bitwise C++ pointer-order mass claim is made.

Each ordinary residue contributes its source free-residue mass minus H2O,
followed by its modifications. One water is added per chain, including an empty
chain. B/Z/X are source-known records with free mass zero, giving internal mass
−H2O. J/U/O are supported. Unknown, lowercase and non-ASCII residue scalars are
reported as unsupported. This local rule does not change AASequence's existing
unknown-composition policy.

After sequence sections and water, source order is N-terminal modifications,
C-terminal modifications, unlocalised modifications multiplied by their signed
occurrence (default one), labile modifications, and global modifications.
Global location matching counts only ordinary sequence elements, once each, and
only one-byte locations equal to that residue. Range/ambiguous elements and
terminal location names are not counted. Global isotope replacements are ignored.

m/z is `(mass + charge * PROTON_MASS_U) / abs(charge)`. Adducts supply only summed
`charge * occurrence`; their formula/name masses are ignored. Charge products and
sums must fit `i32`; `i32::MIN` uses a safe unsigned absolute value. Finite signed
mass/mz results remain valid. Consumed nonfinite inputs, arithmetic overflow or
nonfinite results return errors. Ignored label scores and other unused numeric
fields are not validated solely because a bounded copy contains them.

## Preserved C++ defects and diagnostic order

The [C++ issue ledger](../OpenMS_CPP_ISSUES.md) tracks these source defects; the
native initial port deliberately retains them:

- **CPP-008:** inner modifications on `ModifiedRange.elements` are omitted by
  resolution, validation and mass summation, although the parser/writer retain
  them. `(M[+16]A)[+1]` includes only the range's +1.
- **CPP-009:** ambiguity validation compares unmodified residue masses and ignores
  candidate modification availability. Calculation uses the first candidate,
  including its modifications. `(?I[+10]L)` and `(?LI[+10])` differ by 10 Da;
  an unresolved first-candidate modification can silently contribute zero.
- **CPP-015:** only the first alternative's crosslink label is inspected. The
  identifier is reserved before checking mass. A label-only first endpoint can
  suppress the later linker mass, so swapping crosslinked chains can change mass.
  The shared crosslink set applies to ordinary/first-ambiguous/range-wide/N/C
  annotations; unlocalised, labile and global modifications bypass it.

Issue order and source messages are retained, with zero-based positions and
`None` replacing the source unknown-position `SIZE_MAX`. Source message text that
embeds that sentinel retains native `usize::MAX` spelling. Ion issues receive
`Chain i: ` prefixes. The issue predicate does not reject chimeric ions: it can
return true even though scalar mass rejects the ion.

Strict scalar mass checks all chain issues before chimeric state. Optional ion
mass checks empty/chimeric state before resolving chains. An empty ion has mass
zero, even when marked chimeric. Missing/zero m/z charge is checked before any
mass work. These differences are intentional compatibility with the source
operation order, rather than one shared scientific-availability policy.

## Transaction, repeated passes and bounds

One private resolver session shares a lazy registry clone and one cumulative
budget across all chains and source passes. It publishes the registry only after
a public operation completes. Any `Err` leaves caller registry and AST unchanged.
An `Ok` issue/can/optional report is a completed operation: new formula definitions
are retained even when issues exist or its optional value is absent. This checked
rollback differs from source global mutation when an exception interrupts work.

The source's repeated resolution passes are preserved, including warning
repetitions. Strict mass validates a resolved copy, then calculates a freshly
resolved original copy. Chain optional mass resolves copy A, validates another
resolved copy B, then calculates A. Ion optional mass performs A/B validation
passes, then resolves fresh original copies for calculation. Formula interning
can make later passes resolve names that earlier passes could not, so collapsing
passes would change behavior. A direct AST regression names a not-yet-interned
`M[Formula:Cl101]` before a later formula annotation creates that exact key; the
optional method can validate B but omit the still-unresolved mass in A (**CPP-020**). This is
a source-reviewed defect, not a parser-produced spelling or an executed C++ probe.

Every call shares **1,000,000 consumed items**, **50,000,000 work units**,
**256 MiB cumulative logical allocation allowance**, and **4 MiB per copied or
consumed text value**. Fixed-depth AST vector/string payload is precharged before
copying, including ignored fields copied by source. Immutable chemistry handles
share ownership. Formula parsing, registry operations, issue/warning construction,
crosslink keys, numerical loops, vector moves and copied-value destruction use
the same ledger. Registry cloning occurs only when insertion is required.
Repeated source copies and conservative estimates can exhaust the budget earlier
than a different implementation; these are logical bounds, not physical allocator
or process-memory guarantees. Explicit early empty/charge/chimeric branches avoid
copying fields the source does not touch.

## Evidence

[Thirteen direct tests](../tests/proforma_mass.rs) cover all source mass-test
sections and independent structural/arithmetic cases. Source examples include
PEPTIDE, DFPIANGER, oxidation, phosphorylation, carbamidomethylation, terminal,
formula, unlocalised, labile, global, intra/interchain crosslinks, optional
results and charge errors. Rounded source literals are tested as rounded values.
Source comparisons to AASequence are retained with a small tolerance for stored
modification deltas rounded independently of formula-derived masses; the source
+138.068 try example compares against strict calculation, because nearest-mass
resolution selects stored DSS +138.06808.

Four private tests force late work/byte/item failure after interning, compare
completed-report versus scalar-error publication, check shared chain budgets,
and exercise early branches with oversized ignored AST state. Existing resolver
private/public and parser/writer/JSON/registry tests remain adjacent regression
evidence. The [manifest](../tests/data/proforma_mass_provenance.json) records exact
source hashes and class-test sections. No C++ mass execution or full SDK build is
claimed. Existing parser/writer/runtime probes retain their separate provenance.
