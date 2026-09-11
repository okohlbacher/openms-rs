# ProForma annotation data and text writing

`chemistry::proforma` provides the complete owned annotation data model and both
source text-serialization overloads from OpenMS4-core
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. Both source text grammars and
structured errors are also available; see [text parser support](PROFORMA_PARSER_SUPPORT.md).
All public operation groups of the pinned `ProForma` header now have native equivalents,
with the documented source behavior and checked resource/numerical boundaries.

The optional [JSON transport](PROFORMA_JSON_SUPPORT.md) provides both top-level
read/write operations with the source tagged schema. The
[modification resolver](PROFORMA_RESOLUTION_SUPPORT.md) fills chemistry handles
against a caller-owned registry. [Mass/mz operations](PROFORMA_MASS_SUPPORT.md)
include predicates, diagnostics and optional results. [AASequence conversion](PROFORMA_CONVERSION_SUPPORT.md)
provides both directions, all policies and diagnostics. [Spectrum generation](PROFORMA_SPECTRA_SUPPORT.md)
implements the six ordinary/crosslink operations using the existing native generators. Existing
[AASequence](SEQUENCE_SUPPORT.md), [modification records](MODIFICATION_SUPPORT.md),
and [MonosaccharideDB](MONOSACCHARIDE_SUPPORT.md) remain separately usable.

## Constructing and writing annotations

```rust
use openms::chemistry::proforma::{
    CvAccession, CvDatabase, Modification, ModificationTag, Peptidoform,
    SequenceElement, SequenceSection, WriteMode,
};

let peptide = Peptidoform {
    sequence: vec![SequenceSection::Element(SequenceElement {
        amino_acid: 'M',
        modifications: vec![Modification {
            alternatives: vec![(
                ModificationTag::CvAccession(CvAccession {
                    database: CvDatabase::Unimod,
                    accession: "35".into(),
                }),
                None,
            )],
            ..Default::default()
        }],
    })],
    ..Default::default()
};
assert_eq!(peptide.to_text(WriteMode::Lossless)?, "M[UNIMOD:35]");
Ok::<(), openms::Error>(())
```

`Peptidoform::to_text(mode)` and `PeptidoformIon::to_text(mode)` return
`Result<String>`. Both borrow their input and publish the complete string only
on success. `WriteMode::default()` is `Lossless`; the mode argument is explicit.

These methods serialize stored state, including unusual or grammatically invalid
state. They do not resolve accessions, check glycan names, parse formulas, escape
annotation text, infer label types, or validate linkage partners. Writing an AST
does not establish that its text will parse or round-trip through upstream.

## Complete data model

All source annotation variants and fields are retained in concrete public Rust
types, using owned `String`/`Vec`, signed `i32` counts and charges, and `Option` for
optional values. The source nested types are flattened inside this module.

| Source representation | Native representation |
| --- | --- |
| Conversion policy and issue vocabulary | `ConversionPolicy`, `ConversionIssueType`, `ConversionIssue`; `issue_type` replaces the field named `type` |
| Unknown issue position, source `SIZE_MAX` | `ConversionIssue.position: Option<usize>`; `None` is unknown |
| CV accession and optional named-modification hint | `CvDatabase`, `CvAccession`, `NamedMod` |
| Mass delta and source hint | `MassDelta`, `MassDeltaSource`, including exact `original_text` |
| Formula and glycan name/formula components | `FormulaTag`, `GlycanComposition`, `GlycanComponent::{Name, Formula}` |
| Free annotation and location constraints | `InfoTag`, `PositionConstraint`, all flags and residue order |
| Seven-way modification tag union | `ModificationTag`, with a variant for every source tag |
| Cross-link/branch/ambiguity labels | `LabelType`, `Label` with `label_type`, identifier and optional score |
| Alternative tags and resolved chemistry | `Modification.alternatives` and `resolved_mod: Option<Arc<ResidueModification>>` |
| Ordinary element, ambiguous region, modified range | `SequenceElement`, `AmbiguousRegion`, `ModifiedRange`, `SequenceSection` |
| Unlocalised, labile and position-specific global modifications | `UnlocalisedMod`, `LabileModification`, `GlobalModification` |
| Global isotope replacement | `IsotopeReplacement`, `GlobalModEntry` |
| Signed simple charge or adduct list | `ChargeState::{Simple, Adducts}`, `AdductIon` |
| Chain, chain group and linkage sites | `Peptidoform`, `PeptidoformIon`, `CrossLinkGroup` |

Ordinary `Clone` owns copies of annotation strings and vectors; immutable resolved
chemistry stays shared through `Arc`, so custom registry records outlive their
original registry. Writers neither traverse nor serialize resolved records.
`PartialEq` is a native convenience with ordinary IEEE float comparison and the
existing modification record's full-value equality; it does not imply pointer
identity. Stored defaults use explicit empty/zero values where supplied, replacing
uninitialized C++ aggregate scalar states.

## Source text behavior preserved

Both modes preserve vector order, repeated modifications, alternative order,
duplicate locations and signed/zero occurrence counts. The header describes
canonical sorting, but the implementation does not perform it.

* CV prefixes use `UNIMOD`, `MOD`, `RESID`, `XLMOD`, `GNO`; accession spelling is
  otherwise untouched. Named hints use `U`, `M`, `R`, `X`, `G`.
* Lossless masses copy a nonempty `original_text` exactly, even if it disagrees
  with the independent numeric mass. Otherwise masses use fixed four decimals
  and prepend `+` when the numeric value compares greater than or equal to zero.
  Source negative zero therefore becomes `+-0.0000` in this branch.
* Canonical label scores use fixed two decimals. Lossless label scores begin with
  classic-locale default stream formatting, six significant digits. Formatting
  flags persist: a previously formatted numeric mass changes subsequent lossless
  scores to fixed four decimals. A preserved mass spelling makes no such change.
  Every chain in an ion starts a fresh stream. The existing ParamValue stream
  formatter supplies the default floating representation.
* Label type is stored but ignored by the writer; identifier and score are used.
  Empty INFO plus a label writes `[#label]` inside square brackets, but a labile
  brace writes `{INFO:#label}`. Empty alternatives write `[]` or `{}`.
* Glycans emit each component followed by its count, including one, zero and
  negatives, without separators. Inline formula components retain their optional
  charge. Formula strings and glycan names receive no chemical interpretation.
* Chain names write `(>name)`. **The ion name is omitted.** A single-chain writer
  omits its charge; the ion writer emits chain charges only for a chimeric ion,
  joins chains with `+` or `//`, then emits an optional overall ion charge.
* An adduct list writes its entries followed by the absolute summed charge and
  sign, for example `/[Na:z+1^2,H:z-1]1+`. Empty lists write `/[]0+`. The source
  text parser does not consume that trailing aggregate-charge suffix; a general
  adduct text round trip is not promised.

Empty chains, empty ranges/ambiguities and otherwise unescaped source string
fields remain serializable. Unicode annotation text is preserved. Only individual
`char` fields (`amino_acid` and position-constraint residues) require ASCII in the
native writer; ASCII control characters remain literal source bytes.

## Checked limits and atomicity

Every call shares fixed bounds across all chains: **4 MiB output**, **1,000,000
consumed collection items**, and **50,000,000 work units**. Public constants name
these limits. Each collection is charged before iteration, each copied text is
charged by byte length, and output growth charges movement of existing bytes.
Geometric requested capacity is capped at the output limit. Each float reserves
2,048 work units before its bounded scalar formatting, covering the reused
formatter's temporary work as well as fixed-decimal conversion. These are explicit
logical limits, not a claim about allocator-internal bookkeeping or the size of
the caller's already-owned AST. Numeric-heavy output can reach the work limit
before the byte or item limit.

Consumed nonfinite masses/scores and non-ASCII character fields error. A mass
ignored in favor of exact lossless spelling remains unconsumed, as do ignored ion
names, chain charges and resolved records. Adduct multiplication and each running
sum use checked `i32` arithmetic; overflow errors replace source signed-overflow
undefined behavior. The absolute value of `i32::MIN` is printed through an
unsigned magnitude, a defined native extension over source `abs(INT_MIN)`.
All failures leave input ownership/state unchanged and return no partial string.

## Source compatibility boundaries

[Ordinary and crosslink spectrum generation](PROFORMA_SPECTRA_SUPPORT.md) completes
the pinned public operation groups. It composes the separate ordinary and XLMS
backends and retains their explicit numerical domains and known source defects.

Source itself has notable boundaries: parsing does not consult MonosaccharideDB,
glycan tags remain unresolved for source mass calculation, isotope-global entries
are ignored by that mass calculation, and adduct formula masses do not contribute
to source m/z. A future native standards extension must be distinguished from
those source behaviors. The full dependency analysis is retained in the parent
implementation handoff. Native header coverage does not imply complete ProForma
standard conformance or scientific correction of preserved source defects.

## Evidence

The [source header](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/CHEMISTRY/ProForma.h)
defines the model; [writer implementation lines 329–616](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L329)
defines the operation order. [Native tests](../tests/proforma.rs) directly construct
ASTs for the upstream writer assertions at source test lines 298, 548, 571, 635,
649, 1146–1167 and 1337–1407. Independently constructed examples cover every AST
variant, formatting state/rounding, all hints, charge extrema, owned custom
chemistry, Unicode and resource/error boundaries.

[The provenance manifest](../tests/data/proforma_provenance.json) records source
hashes for the original writer-only increment. Its historical grammar-fixture
status is preserved; the later [parser manifest](../tests/data/proforma_parser_provenance.json)
records the now-executed 176 positive and 22 negative grammar cases.

A separate [compiled source-writer probe](../tests/data/proforma_writer_probe_provenance.json)
extracts the exact public annotation declarations and complete private C++ writer
into a standalone translation unit. Only standard includes, DLL decoration,
namespace/class enclosure and an opaque unused modification-pointer declaration
are supplied. No scientific method is substituted and no C++ parser/backend is
linked. All 160 outputs match [native comparisons](../tests/proforma_writer_probe.rs):
twenty exact f64 bit patterns, both modes and four formatting/chain scenarios.
The cases include signed zeros, subnormals, extreme finite values and transitions
between default, fixed-four and fixed-two precision. This is an extraction probe,
not a full SDK build or complete parser-fixture parity.

The [text parser extraction probe](../tests/data/proforma_parser_probe_provenance.json)
adds 476 executed C++ outcome comparisons, including both text modes on accepted
inputs and code/position/message on rejected inputs. Its exception adapter and
scope are documented in [parser support](PROFORMA_PARSER_SUPPORT.md).
