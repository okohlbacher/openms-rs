# Peptide sequences and mass annotations

`chemistry::AASequence` stores uppercase peptide/protein residue sequences independently of whether their chemistry is fully known. It supports the canonical residues, U/O/J, ambiguous B/Z/X, named registry modifications and numeric monoisotopic annotations. The global modification registry remains immutable; resolved records use shared ownership and anonymous annotations are owned by the sequence and survive cloning, slicing, digestion and identification serialization without global registration.

The implementation is derived from [AASequence.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/AASequence.cpp), [ResidueModification.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/ResidueModification.cpp), and [AASequence_test.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/AASequence_test.cpp), pinned at `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. Integration test provenance records source hashes in [sequence workflow tests](../tests/sequence_workflow.rs). No C++ build or execution is required or claimed.

## Known chemistry is explicit

`formula() -> Result<EmpiricalFormula>`, `mono_mass() -> Result<f64>` and `average_mass() -> Result<f64>` expose separate questions. This is a native API change from earlier versions of this port, which returned these values directly.

| Sequence content | Formula | Monoisotopic mass | Average mass |
| --- | --- | --- | --- |
| Canonical/U/O/J residues and formula-bearing registry modifications | Known | Known | Known |
| Custom mass-bearing record without a delta or absolute formula | Unsupported | Known when base or explicit absolute residue mass is known | Unsupported |
| Changed residue with a supplied absolute formula | Known, including replacement of B/Z/X | Known | Known |
| Bare B, Z or X | Unsupported | Unsupported | Unsupported |
| Absolute B/Z/X numeric tag, e.g. `X[999]` | Unsupported | Known | Unsupported |
| Anonymous tag on a known residue or terminus | Unsupported | Known | Unsupported |
| Empty sequence | Empty formula | Zero | Zero |

`Error::Unsupported` means the requested chemistry is unavailable. `Error::InvalidValue` means the input or calculated chemistry is invalid. The parser does not invent average masses or formulas from a monoisotopic number. Bare B/Z/X never contribute fictitious zero or negative-water masses. J retains the existing isoleucine/leucine composition.

For example, `PEPTX[999]IDE` has the neutral mass of PEPTIDE plus 999 Da, while its formula and average mass remain unknown. `PEPBX[999]IDE` still has an unresolved B and therefore no complete monoisotopic mass. An absolute tag can supply an explicit internal mass for any of B/Z/X; a signed difference on any unresolved residue, including `X[+0]`, is rejected.

`subsequence`, `prefix` and `suffix` recalculate chemistry for their selected residues. Removing all ambiguous or anonymous components restores formula and average-mass access. Original terminal annotations survive only when their original terminus survives. Empty slices return the ordinary empty sequence.

## Fragment and charged formulas

`formula_for(PeptideFragmentType, charge: i32) -> Result<EmpiricalFormula>`
implements the full finite `ResidueType` surface of
[AASequence::getFormula](https://github.com/okohlbacher/OpenMS4-core/blob/6bfc0e4711105f4eda2fea86812a83af7c7e791f/src/openms/source/CHEMISTRY/AASequence.cpp#L383),
using the corrections in
[Residue.h](https://github.com/okohlbacher/OpenMS4-core/blob/6bfc0e4711105f4eda2fea86812a83af7c7e791f/src/openms/include/OpenMS/CHEMISTRY/Residue.h#L63).
The method treats **all supplied residues** as the requested fragment. Select
residues with `prefix`, `suffix` or `subsequence` first when a particular cut is
wanted. The existing `formula()` remains the cached full-formula convenience API.

| Fragment type | Terminal annotations retained | Correction to internal residue sum |
| --- | --- | --- |
| `Full` | Both | `H2O` |
| `Internal` | Neither | None |
| `NTerminal` | N | `H` |
| `CTerminal` | C | `OH` |
| `AIon`, `BIon`, `CIon` | N | `C-1O-1`, none, `NH3`, respectively |
| `XIon`, `YIon`, `ZIon` | C | `CO2`, `H2O`, `H-1N-1O`, respectively |
| `Zp1Ion`, `Zp2Ion`, `Precursor` | Neither | None: source fallback |
| `BIonMinusH2O`, `YIonMinusH2O`, `BIonMinusNH3`, `YIonMinusNH3` | Neither | None: source fallback |
| `NonIdentified`, `Unannotated` | Neither | None: source fallback |

The fallback variants deliberately return the internal formula. This source
method does not apply the radical or neutral-loss chemistry supported separately
by the theoretical spectrum generator.

Requested charge initializes the formula's charge metadata. It does not add
natural hydrogen atoms or divide the result by charge. Charges already carried
by retained modification formulas then add to that value. Empty sequences return
an empty, uncharged formula for every type and charge, including typed terminal
annotations retained by modification generation. This differs from RNA's
hydrogen-based charge convention; their fragment enums are separate types.

Retained anonymous/mass-only terminal annotations return `Error::Unsupported`;
discarded terminal annotations do not block formula access. Every included
residue still needs known composition. Custom absolute-formula replacement and
delta-formula precedence remain the same as for `formula()`.

```rust
use openms::chemistry::{AASequence, PeptideFragmentType};

let peptide = AASequence::parse("(Acetyl)PEPTIDE")?;
let b3 = peptide.prefix(3)?.formula_for(PeptideFragmentType::BIon, 2)?;
assert_eq!(b3.charge(), 2);
```

Formula construction adds N-terminal delta, C-terminal delta, internal residues,
then the type correction in source order. Atom and charge overflow return
`Error::InvalidValue`; signed finite atom counts in a fragment formula are
preserved as source algebra. No sequence or annotation is mutated on failure.
Each standalone call allows 50 million conservative map/traversal work units and
256 MiB of cumulative formula-node/scratch accounting, charged before allocation.
Graph callers can share these remaining allowances across formula operations.
The counters bound temporary maps as well as the returned formula; they are not
measurements of allocator RSS.

[Sequence formula tests](../tests/sequence_formula.rs) include the literal
`ACDEF` full/charged/b-ion formulas from `AASequence_test.cpp:341–346`, all 19
fragment kinds, terminal selection, charge extrema, slicing, unavailable
composition, custom absolute formulas, signed output and checked charge overflow.
Private tests exercise cumulative budgets and empty typed terminal states.

## Named and numeric syntax

Named annotations retain the existing registry syntax, including nested parentheses, UniMod accessions and explicit peptide/protein terminal names. For example, `(Acetyl)AC(Carbamidomethyl)M(Oxidation)K` remains supported. Existing named-annotation terminal fallback is retained.

Numeric brackets use these conventions:

| Input | Interpretation |
| --- | --- |
| `M[+15.99]` | A signed difference from M's internal mass; resolves to Oxidation at the written precision |
| `M[147.0354]` | An absolute internal residue mass; resolves to Oxidation at this precision |
| `M[147.035405]` | A more precise absolute mass; remains anonymous because the registry match lies outside its tolerance |
| `X[999.000]` | An explicit internal monoisotopic mass of 999 Da with original spelling retained |
| `[+42]PEPTIDE`, `[+42].PEPTIDE`, `n[+42]PEPTIDE`, `.[+42]PEPTIDE` | N-terminal mass difference |
| `PEPTIDEc[-1]`, `PEPTIDE.[-1]` | C-terminal mass difference |
| `n[43.0183900319]PEPTIDE` | Absolute N-terminal component mass; subtract the neutral H mass to obtain its delta |
| `PEPTIDEc[16.0187240319]` | Absolute C-terminal component mass; subtract the neutral OH mass to obtain its delta |

A leading plus **or minus** sign denotes a difference. An unsigned value denotes an absolute mass. Absolute residue annotations describe the residue's internal peptide mass, not its free amino-acid mass. Terminal absolute values describe only the terminal H or OH component. The neutral H and OH masses come from the same pinned element table as the formula model: approximately 1.0078250319 and 17.0027400319 Da.

A numeric annotation after a residue **always remains on that residue**. Only Anywhere registry modifications are considered there. Request a numeric terminal annotation before the sequence or with explicit n/c/dot notation. This deliberately differs from C++ retrying a boundary residue's numeric tag as a terminal modification. For example, native `Q[111]PEPTIDEK` preserves an anonymous Q internal mass; use `n[-17]QPEPTIDEK` for the source's approximate N-terminal pyro-Glu lookup. `C[143]PEPTIDEK` likewise remains a residue tag; `n[+40]CPEPTIDEK` requests the terminal Pyro-carbamidomethyl match.

The explicit attachment rule prevents chemical identity from changing when an internal anonymous tag becomes terminal through slicing or when a setter-created annotation is serialized and read again. C++ can rely on its mutable global registry to remember previously created anonymous annotations; the Rust registry stays immutable. No additional serialization syntax or hidden process history is introduced.

## Registry resolution and representation

The decimal spelling determines the registry lookup window, as in the pinned source:

- Without a decimal point, select the first source-ordered match within an inclusive 0.5 Da difference window.
- With a decimal point, select the closest match whose error is strictly less than `10^(-number of written decimal places)` Da. Ties keep source order. Written trailing zeros affect this lookup precision.
- Exclude zero-delta registry entries for a nonzero query.
- Use the actual residue and requested attachment specificity. Explicit numeric N/C lookups use peptide-terminal specificity; named annotations retain the existing protein-terminal fallback.
- If the base residue is B/Z/X, an unsigned tag stays an absolute owned annotation. No delta lookup is performed against C++'s placeholder mass or wildcard-X behavior.

Resolved annotations become `SequenceModification::Known(Arc<ResidueModification>)`.
`known()` borrows the record; cloning a peptide shares immutable chemistry without
copying its strings or leaking storage. `parse_with_registry(input, &registry)`
and all named/numeric `set_*_with_registry` variants use caller-supplied records
only during resolution. Dropping or extending that registry does not invalidate
an existing sequence or change its chemistry.

Display uses registry names, falling back to the full ID when a caller-supplied record has no short name. `to_accession_string()` uses UniMod or OBO accessions,
falling back to a custom name when no accession exists; reconstruction needs an
equivalent registry. `to_unimod_string()` now returns `Result<String>`: known
UniMod records retain their accessions, while other named records become absolute
internal-residue or H/OH-terminal mass brackets, following the source fallback.
No fictitious UniMod ID is generated. An unknown mass returns `Unsupported`; a negative terminal absolute component is rejected because a signed native bracket denotes a delta. An empty sequence exports an empty string, as in the source, including empty typed states carrying terminal records.
This export loses non-UniMod chemical identity and uses Rust's shortest
round-trippable decimal mass spelling. An unresolved numeric annotation becomes `SequenceModification::MassTag(MassTag)` with private, validated fields. Display and UniMod output preserve its exact original bracket content. Terminal Display uses `.[content]`, while the annotation's full identifier is `.n[content]` or `.c[content]`; residue full identifiers are `Q[content]`, `X[content]`, and so on.

`SequenceModification` exposes `known()`, `mass_tag()`, `name()`, `full_id()`, `record_id()`, `origin()` and `term_specificity()`. UniMod record IDs are optional; annotations without a UniMod entry report `None`. Anonymous names are the original numeric text, and their full IDs include attachment identity. Its `diff_formula()` and `diff_average_mass()` return Unsupported for anonymous annotations. `diff_formula()` also rejects a known custom record that carries a mass change without a delta formula (even when an absolute replacement formula is available); its declared mass differences remain accessible. `diff_mono_mass()` is available when a base mass is known; an absolute B/Z/X annotation has no known unmodified mass difference.

`MassTag` exposes its original `input()`, written `mass()`, `is_delta()`, optional `delta_mono_mass()` and optional `residue_mono_mass()`, `full_id()`, `origin()` and `term_specificity()`. The residue mass is the resulting internal residue mass and is absent for terminal tags. Accessors borrow annotations; consumers can clone them into owned records. No anonymous formula or mutable registry reference is exposed.

## Mass calculations and mutation

Formula-bearing residue modifications use their atom-formula deltas for monoisotopic mass. Known terminal modifications retain their declared source monoisotopic deltas. Average masses use the complete formula whenever it exists. The original declared terminal-mass correction remains in effect.

The [modified-peptide generator](MODIFIED_PEPTIDES_SUPPORT.md) can also apply
shared records from a caller-owned registry. When a record describes a nonzero mass or atom-formula difference, a delta
formula takes precedence. Otherwise an absolute formula replaces the free
residue formula, with water removed for the internal representation. This can
restore known composition and mass for B/Z/X. With neither formula, the record
keeps declared monoisotopic mass while composition and average mass remain
unavailable. For a residue, a nonzero declared absolute
mass is a free-residue mass; otherwise the declared delta is added to the
unmodified free residue, then water is subtracted. Formula-bearing records
retain formula precedence. A record with no formula or mass differences is
a no-op even when it stores absolute masses or formulas, matching the source’s residue-term
exception. Termini continue to use their declared deltas. Unlike the source’s
retained unmodified formula for mass-only records, the native formula API
reports missing composition explicitly.

`mz(charge)` requires a nonempty sequence and a known monoisotopic mass. The existing checked charge rules and proton mass apply. `fragment_ions(max_charge)` supports known monoisotopic numeric tags for b/y ions, propagating only the annotations retained by each fragment. Prefix and suffix masses accumulate independently to preserve a small complementary fragment beside a very large finite tag. Terminal differences combine before residue masses in full-mass calculation, preserving the peptide mass when large opposite terminal tags cancel.

Formula-dependent operations must obtain a known formula explicitly. The theoretical spectrum generator handles supported mass-only monoisotopic ions and documented neutral-loss cases, while rejecting formula-dependent isotope or unavailable loss calculations. See [theoretical spectra support](THEORETICAL_SPECTRA.md) for the precise options.

Existing `set_modification`, `set_n_terminal_modification` and `set_c_terminal_modification` operate on registry names; an empty name removes the annotation. New `set_mass_tag(index, text)`, `set_n_terminal_mass_tag(text)` and `set_c_terminal_mass_tag(text)` accept bracket **content** and replace the relevant annotation using the same numeric lookup. They do not cumulatively add to a previous annotation. All setters validate a staged copy before committing. Invalid syntax, unknown requested names, impossible known masses or out-of-range indices leave the sequence unchanged.

## Validation and remaining boundaries

Residues must be uppercase ASCII codes. Lowercase residues, literal stops, whitespace and inline protein modification syntax outside the supported grammar remain invalid. Numeric tag content has 1–1024 characters and accepts decimal digits, at most one decimal point and an optional leading sign. Scientific notation is deliberately unsupported: the source derives lookup tolerance by counting characters after a decimal point, so exponent notation would make that precision rule misleading. Empty numbers, NaN/infinity, overflow and nonzero decimal values that underflow to zero are rejected. Signed zero is retained as a delta and remains invalid on an unresolved base residue.

Duplicate annotations at the same attachment, unterminated brackets, excessive parenthesis nesting, residues after a C-terminal delimiter, and negative or nonfinite calculated known residue/peptide masses are rejected. Sequence construction rejects negative atom counts in fully known formulas; `formula_for` preserves signed fragment algebra as described above. Finite floating arithmetic still has ordinary machine precision; this is not arbitrary-precision mass arithmetic.

Numeric setters use replacement semantics. The source's cumulative `setModificationByDiffMonoMass` and arbitrary multi-modification merging are not provided by these setters. Anonymous formulas, average masses and undeclared isotope distributions remain unavailable unless a chemically defined record is supplied and the annotation explicitly replaced. The registry and known-modification coverage are described separately in [modification support](MODIFICATION_SUPPORT.md).

`tests/sequence_mass_tags.rs` contains 13 focused tests covering pinned integer/decimal modification lookups, historical source mass goldens with their appropriate numerical tolerance, fractional mass differences, H/OH terminal conventions, ambiguous-residue chemistry, slicing recovery, owned annotation identity, parser rejection and atomic setters. Independent `sequence_workflow` and `sequence_identification` tests exercise digestion, indexing, identification filtering/resolution, theoretical ions/losses, idXML round trips, and the large-mass cancellation regressions. Test provenance distinguishes source-derived expectations from deliberate native corrections.

`AASequence`, `SequenceModification` and `MassTag` have deterministic complete
value ordering consistent with equality. This prevents identification algorithms
from merging different custom chemistry that happens to share sequence text.
Ordering is native and independent of allocation addresses; it is not a promise
to reproduce C++ pointer order.

Anonymous mass tags retain their exact decimal spelling in immutable shared strings. Owned clones and fragment slices keep these strings alive independently; setters replace the affected tag without changing other sequences. Their text, mass values, equality and ordering are unchanged. Generation payload estimates continue to conservatively count the full logical ownership of these strings.
