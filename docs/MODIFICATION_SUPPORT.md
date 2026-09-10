# Peptide modification support

`AASequence` holds named residue and terminal modifications from an immutable
shared registry, alongside owned numeric mass annotations. The default `ModificationsDB`
contains 3,127 specificity records: 3,035 from pinned UniMod/OpenMS custom XML
plus 92 mono-link records from pinned XLMOD. Caller-owned registries can load
PSI-MOD or XLMOD OBO streams and validated custom records. Original XML, deterministic generation,
hashes and data licenses are in [resources/modifications](../resources/modifications/README.md).
The implementation has no C++ dependency or runtime resource lookup.

## Native API

```rust
use openms::chemistry::{AASequence, ModificationsDB, TermSpecificity};

let mut peptide = AASequence::parse("PEPC(Carbamidomethyl)M(Oxidation)K")?;
peptide.set_n_terminal_modification("Acetyl")?;
assert_eq!(peptide.as_str(), "PEPCMK");
assert!(peptide.is_modified());
let entry = ModificationsDB::global().get_modification(
    "UniMod:35", Some('M'), Some(TermSpecificity::Anywhere),
)?;
assert_eq!(entry.name(), "Oxidation");
# Ok::<(), openms::Error>(())
```

Names, full IDs, full names and case-insensitive `UniMod:` accession prefixes
resolve with explicit site/specificity filters. A query matching different full
IDs is an error; duplicate specificity rows with the same full ID use source
order. Peptide N/C termini are tried before protein N/C termini. Residue suffix
notation tries an internal modification first, then the applicable termini.
The single-residue sequence `A(Amidated)` correctly falls back to its C terminus.
This fallback applies to named annotations. Numeric brackets after a residue
always attach to that residue; a leading tag denotes the N terminus, and
`n`, `c`, or dot markers explicitly select a terminus. This keeps anonymous annotations stable after slicing
or setter operations. The C++ parser's numeric terminal retry is deliberately
not reproduced; see [sequence support](SEQUENCE_SUPPORT.md).

Examples accepted from the C++ tests include:

- `M(Oxidation)`, `C(UniMod:4)` and `K(Label:13C(6)15N(2))`;
- `(Acetyl)PEPTIDE`, `.(UniMod:1)PEPTIDE`, and `n(Acetyl)PEPTIDE`;
- `PEPTIDE(Amidated)`, `PEPTIDE.(UniMod:2)` and `PEPTIDEc(Amidated)`;
- `(UniMod:51)CPEPTIDE` for a protein-terminal specificity.

`as_str()` returns bare residues. `Display` includes modifications with explicit
terminal dots and uses full IDs for protein termini to preserve that specificity
when reparsed. Source generator paths can also create terminal records on
residue slots or on mismatched termini; these typed states remain available
natively, but ordinary sequence text cannot preserve every such placement.
`to_unimod_string() -> Result<String>` follows the OpenMS accession convention,
using absolute internal-residue or H/OH terminal masses when a record has no
UniMod ID; unavailable masses return errors. `to_accession_string()` preserves
OBO accessions when present. UniMod accessions alone cannot distinguish a protein specificity from a peptide
specificity with the same accession. Modification parentheses may nest to 128
levels; malformed and duplicate attachments are errors.

Setters replace a modification, or remove it when given an empty name. They
validate a temporary value before committing it. `subsequence`, `prefix`,
`suffix` and digestion retain residue modifications and retain a terminal
modification only if the corresponding original terminus remains. Digestion
uses the bare sequence cleavage rule, as in the C++ implementation.

## Mass and formula conventions

Sequence `formula()`, `mono_mass()` and `average_mass()` return `Result`.
The conventions below apply when the corresponding chemistry is known.
[Sequence support](SEQUENCE_SUPPORT.md) explains B/Z/X, numeric lookup precision,
absolute versus delta tags, and formula/average errors for mass-only annotations.
Sequence modification accessors borrow a `SequenceModification`; `known()`
returns a borrowed record owned through `Arc`, while `mass_tag()` exposes an
owned tag. Matching `SequenceModification::Known` exposes its shared handle;
registry handle methods clone it. Sequences and generated variants remain valid after a caller registry is dropped. `record_id()` returns
`Option<u32>` because OBO and custom records need not have UniMod IDs.

The neutral full formula contains all atom deltas, including isotope
substitutions. Residue modification masses are recalculated using the pinned
OpenMS element table, matching `Residue::setModification` when a delta formula
is present. Terminal mono masses retain the database's declared delta values,
matching `AASequence::getMonoWeight`. The average mass is evaluated from the
full formula. Consequently a terminally modified peptide's mono mass may differ
slightly from its formula's mono mass. Both are intentional source conventions.

For example, the pinned element table gives a `13C6/15N2` lysine label delta of
8.014200 Da; the XML declares 8.014199 Da. Tests distinguish these instead of
mixing tables or silently adjusting expected values. Coarse isotope patterns
operate on the formula and therefore use the formula mass anchor.

b/y fragments include only modifications on retained residues and the retained
terminus. Negative final atom counts, count overflow and nonfinite/negative
peptide mass are rejected. Zero and negative atom deltas remain valid inside
the modification database because they describe losses and substitutions.

Mass search uses declared monoisotopic deltas. `search_by_mass` includes exact
tolerance endpoints. Like the upstream best-match API, `best_by_mass` requires
strictly smaller error; tolerance zero returns no best hit. Ties use source
order. Neutral-loss formulas and declared masses are exposed per specificity.

The [modified peptide generator](MODIFIED_PEPTIDES_SUPPORT.md) applies fixed
modifications and enumerates bounded variable combinations. The
[definition-set API](MODIFICATION_DEFINITIONS_SUPPORT.md) manages fixed/variable
search definitions, compatibility, peptide-based inference and delta/absolute
mass matching. Owned modification records can additionally carry explicit
absolute masses; the embedded UniMod/custom records retain source defaults of
zero for these fields. OBO records can independently carry a declared absolute
free-residue formula. Formula precedence, alias resolution and bounded loading
are described in [registry support](MODIFICATION_REGISTRY_SUPPORT.md).
The generator supports mass-bearing custom records without atom formulas:
monoisotopic mass follows source declared/absolute precedence, while sequence
formula and average mass remain unavailable. Formula-dependent isotope
generation rejects those variants; ordinary monoisotopic fragments retain the
known mass shift. See [sequence support](SEQUENCE_SUPPORT.md).

## Remaining compatibility work

ProForma, Unimod XML alternative names/cross references, mutable global
databases and context validation against a parent protein remain unimplemented.
`from_tsv` retains its 11-field format. `AASequence::parse_with_registry` and
registry-aware setters use caller records without replacing the global database.
Historical PSI-MOD is not bundled because its original redistribution grant is
unconfirmed; callers can supply it through the bounded OBO loader.
Classification is retained per specificity, avoiding the source XML handler's
shared modification classification overwrite.

`tests/modifications.rs` contains source-derived formulas and named-notation
examples, data parser checks, declared/computed mass distinctions, neutral
losses, terminal round trips, atomic replacement, digestion preservation and
modified fragment mass checks. The C++ reference was inspected but not executed.
