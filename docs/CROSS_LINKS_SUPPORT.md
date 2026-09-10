# Cross-link modification database

`chemistry::CrossLinksDB` composes the shared native modification registry with
the source's cross-link OBO selection rules. Source references are
`CHEMISTRY/CrossLinksDB.{h,cpp}`, `OBODataProvider.cpp`, and
`CrossLinksDB_test.cpp` at revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`.

The Rust implementation is BSD-3-Clause. The bundled XLMOD vocabulary is the
2016-07-13 release in the pinned source, created by Lutz Fischer and Gerhard
Mayer, and carries **Creative Commons Attribution 3.0 Unported** terms. Its
original attribution is retained in `resources/modifications/XLMOD.obo`; the
resource notices and independently checked provenance describe the data license
separately from the code license.

## Interface and ownership

```rust
use openms::chemistry::{CrossLinksDB, TermSpecificity};

let database = CrossLinksDB::global().database();
let linker = database.get_modification(
    "DSS", Some('K'), Some(TermSpecificity::Anywhere),
)?;
assert_eq!(linker.full_id(), "DSS (K)");
assert_eq!(linker.obo_accession(), Some("XLMOD:02001"));
assert_eq!(linker.diff_mono_mass(), 138.06807961);
# Ok::<(), openms::Error>(())
```

`global()` initializes the bundled vocabulary's 56 cross-link specificity
records once, without runtime filesystem or network access. It returns an
immutable shared instance. `database()` exposes
the existing lookup, mass-search, record iteration, and owned-handle APIs;
there is no second implementation of those operations.

`from_obo(reader, &OboReadOptions)` creates a caller-owned database using the same
bounded reader as `ModificationsDB`. It forces `cross_links_only=true` on a copy
of the options and preserves the caller's input, line, term, record, alias, and
registry payload limits. Options and input are not modified. `database_mut()`
permits checked extensions of an owned instance. As in the source's inherited
`addModification`, explicit custom additions need not be ontology cross-links;
this never mutates the global instance.

The registry's `Arc<ResidueModification>` handles remain valid after extension
or after a caller-owned database is dropped. Borrowed accessors remain useful
when a record is needed only while the database is borrowed. XLMOD accessions
are independent of optional UniMod record IDs; an ontology record does not gain
a fabricated UniMod ID.

`all_search_modifications()` returns full IDs with a nonempty OBO accession,
sorted ascending. The source stores PSI-MOD and XLMOD accessions in the same
field for this purpose. This is distinct from the general registry's
UniMod-based search listing. The source vector does not deduplicate names;
duplicate records, if explicitly added, retain duplicate output names. Custom
records without an OBO accession remain available through ordinary lookup but
are omitted from this particular list.

## Preserved source selection rules

- `reactionSites=1` is excluded. The implementation does not replace this rule
  with a stricter requirement that the field exist and equal two.
- The two specificity lists are unioned and deduplicated. Each ordinary origin
  creates a record; B/J/Z ambiguity origins are excluded. Generic X/Anywhere
  entries are omitted by the source's specificity filter.
- `Protein N-term` and `Protein C-term` become peptide N/C-terminal records in
  CrossLinksDB, with wildcard origin X represented natively as `None`. The
  general registry uses protein-terminal specificity for those same tokens.
- Zero-delta generic terminal entries fail the source's final specificity test.
- XLMOD names, ontology accessions, synonyms, source ordering and declared
  monoisotopic differences are retained. Negative cross-link mass differences
  are legitimate: EDC is `−18.01056027` Da.

The source does not preserve which reactive side a specificity came from.
For example, EDC produces eight records: D/E/K/S/T/Y and N/C termini. Those
records describe the union of possibilities; they do **not** assert that every
pair of listed sites can react together. Reaction-pair chemistry, spacer
geometry, cross-linked peptide graphs, and cross-link search/scoring are not
provided by this database class. Formulas absent from source records are not
inferred from their masses.

The wrapper shares the OBO reader's checked malformed-input and resource
policies. Input parsing completes before a caller-owned registry extension is
committed. See the shared registry documentation for provider alias resolution,
finite-mass validation, ordering, and bounded payload details.

## Verification

`tests/cross_links.rs` ports the upstream DSS/BS3/EDC expectations, including
eight EDC specificities and the two N-terminal matches near 138.068 Da. It also
checks ontology aliases, signed masses, the sorted search list, source filter
edges, preserved reader limits, caller-owned additions, and handle validity
after database drop. The independent registry reference suite compares the
complete pinned XLMOD projection and retains separate source/data provenance.
No C++ compilation or newly executed C++ result is claimed.
