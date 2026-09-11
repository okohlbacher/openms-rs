# PeptideEvidence value and key operations

The complete class-specific `METADATA/PeptideEvidence.h` surface at SDK
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b` is represented by
`identification::PeptideEvidence` and `FlankingResidue`. This closes the evidence
record, not the containing identification classes or all identification workflows.

| Source operation | Native operation |
| --- | --- |
| Default construction | `Default`, with empty accession, missing endpoints and unknown flanks |
| Full-field construction/accessors | Public fields; checked accession/range constructor `new` |
| Copy, move, assignment, destruction | `Clone`, ownership and ordinary Rust assignment/drop |
| Equality and inequality | All five fields through `PartialEq`/`Eq` |
| Less-than | Accession, start, end, before and after; missing positions first and source character ordering for valid flanks |
| `hasValidLimits` | `has_valid_limits`, with the documented inclusive-range correction |
| `std::hash` key use | Rust `Hash` for all fields, consistent with native equality |
| Unknown position -1 | `None` |
| N-terminal position 0 | `Some(0)` |
| Unknown/N-terminal/C-terminal flank X/[ /] | `Unknown`/`NTerminus`/`CTerminus`, round-tripped by `code`/`from_code` |

Coordinates are zero-based and inclusive. `positions()` returns a checked
inclusive range. The existing native predicate accepts a legitimate first-residue
mapping 0..=0 and rejects reversed intervals, correcting [CPP-045](../OpenMS_CPP_ISSUES.md).
Missing endpoints can be retained independently but do not form valid limits.
Negative coordinates cannot be represented; native `usize` coordinates can also
exceed the C++ signed Int range. Format-specific output constraints remain the
responsibility of each adapter.

Validated flanks use uppercase ASCII residues and the three source markers.
Directly constructed invalid public `Residue` variants remain distinguishable
under equality, hashing and total ordering; validation rejects them. Arbitrary
C++ char bytes are not part of the validated native model. Copy/move behavior
uses Rust ownership rather than observable C++ moved-from state.

Hashing supports native collections and considers the same stored fields as
equality. Its numeric bytes are not a C++ hash oracle or a stable interchange
identifier. Existing identification tests exercise inclusive coordinates, missing
values, ordering and downstream operations; hashed deduplication now retains
distinct flank records while collapsing repeated identical evidence.

[Provenance](../tests/data/peptide_evidence_provenance.json) pins the complete
header, implementation, source class test and hash helper. Most source class-test
sections are TODO placeholders; they are not claimed as executed or substantive
C++ reference cases. No C++ executable comparison was run for this value group.
