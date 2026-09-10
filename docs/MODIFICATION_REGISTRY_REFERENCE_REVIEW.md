# Modification registry reference review

The reference tests cover the pinned XLMOD provider projection, independently
authored OBO examples, retained custom chemistry, and anonymous definition
inference. They run against packaged resources and do not require the original
C++ checkout. No C++ binary was built or executed.

## Evidence and redistribution

[registry_provenance.json](../tests/data/registry_provenance.json) records the
pinned revision, source and fixture SHA-256 hashes, and extraction method.
[registry_xlmod_reference.tsv](../tests/data/registry_xlmod_reference.tsv) contains
all 92 monolink and 56 cross-link specificity records projected from the bundled
2016-07-13 XLMOD release. Each row records the original accession, names, full ID,
origin, terminal specificity, decimal monoisotopic delta, its binary64 bits, and
synonyms. The projection is an independent translation of the source provider's
branches, rather than captured C++ output. Literal CrossLinksDB class-test
expectations additionally anchor DSS, BS3, EDC, and the 138.06807961 Da lookup.

XLMOD and the derived projection retain attribution to Lutz Fischer and Gerhard
Mayer and the pinned release's [CC BY 3.0 terms](https://creativecommons.org/licenses/by/3.0/).
The original notice remains in the bundled resource; transformation into rows is
identified explicitly. The older PSI-MOD snapshot is represented only by its
hash, version, and date in the manifest. It is not redistributed. The current
official PSI-MOD license does not by itself establish the historical terms of
that 2008 snapshot.

The three `registry_synthetic_*` fixtures are independently authored BSD-3-Clause
test data. They contain synthetic identifiers and descriptions, not copied
PSI-MOD terms. These fixtures exercise the PSI-style provider path without
depending on the unbundled ontology.

## Scientific behavior checked

[registry_reference.rs](../tests/registry_reference.rs) verifies complete
projection identity and order, with exact floating-point comparisons for every
pinned XLMOD mass. It also checks:

- Source ordering by ontology accession, then sorted distinct residue origins,
  followed by N-terminal and C-terminal expansions.
- B/J/Z origin filtering, wildcard-X filtering, duplicate-site collapse, and
  source union of both cross-linker reaction arms. The result describes allowed
  sites, not the pairing constraints of two chemical arms.
- Monolink exclusion of `reactionSites=2` and cross-link exclusion of
  `reactionSites=1`. Missing or other counts are not assumed to mean either one
  or two. CrossLinksDB forces its projection without modifying caller options.
- Protein terminal specificity for monolinks versus peptide terminal
  specificity for cross-links. Native terminal wildcard `None` corresponds to
  the source's origin `X`; the full identifiers are unchanged.
- A UniMod-linked PSI alias binding every already installed specificity of its
  UniMod target, regardless of the alias's own declared site. Only its accession
  becomes an alias; its display name and synonyms are not separately installed.
  Unresolved targets are counted and discarded, and provider order matters.
- Standalone accession/full-name/synonym lookup, last-record flushing without a
  final newline, and malformed suffixes leaving the prior registry unchanged.
- Bounded input, lines, terms, expanded records, aliases, and registry payload.
  Repeated aliases consume target-expansion work even when their retained
  results deduplicate; the test checks both the exact limit and its successor.

The synthetic chemistry examples test the source's free-residue conventions.
An absolute mass or formula with no mass/formula difference does not change the
residue. When a change is declared, a nonempty delta formula takes precedence;
otherwise an absolute formula replaces the free-residue formula, whose internal
form is obtained by subtracting water. Formula masses take precedence over
declared absolute masses. Terminal annotations use their declared differences,
not an absolute residue replacement. The custom registry can be dropped before
using a parsed peptide or generating variants from retained shared handles.
Literal empty absolute-formula text is absent. Nonempty `H0` or charge-only `+`
text remains an explicitly supplied formula, even though it contains no atoms;
it cannot silently become an ordinary mass-only change.

The earlier
[definition reference suite](../tests/modification_definitions_reference.rs)
now verifies successful inference of an owned anonymous `C[999]` annotation,
including its exact full identifier and mass-tag payload. The implementation
review also checked the source empty-short-ID compatibility rule and the
full-residue mass anchor used by anonymous absolute-mass matching. Anonymous
mass data does not establish an empirical formula or average mass.

## Explicit native boundaries

This is a modification-provider parser, not a general ontology graph parser.
Unknown nonchemical fields are ignored. Native unknown-stanza handling isolates
`[Typedef]` and other stanzas so that their fields cannot overwrite a preceding
modification; the source only recognizes `[Term]` as a boundary. Native parsing
also accepts the valid `Protein N-term`/`Protein C-term` property spellings,
which source whitespace removal otherwise makes unreachable. Invalid numeric
values, quoted properties, and terminal strings are checked rather than
producing incomplete scientific records.

The source's optional preferred PSI label is retained as a synonym. Native
`name()` continues to mean the short modification identifier; it does not add a
second preferred-display-label field. Ambiguous distinct specificity IDs require
caller disambiguation rather than choosing an unspecified source pointer order.

`to_unimod_string()` returns an empty string for an empty annotated sequence,
matching the source's early return. For a non-UniMod annotation it exports known
absolute mass brackets and intentionally loses vocabulary identity. Anonymous
tags retain their native original spelling. A negative absolute terminal
component is rejected even if the whole peptide mass is positive: emitting its
leading minus sign would instead mean a delta in the bracket grammar. The
accession representation and a caller-supplied registry provide a separate path
for retaining named custom chemistry.

## Review outcome

The ten registry reference tests and five updated definition reference tests
pass. Independent review identified the empty annotated export discrepancy and
quadratic alias-deduplication work. Both were corrected, with explicit tests for
empty export and bounded repeated-alias expansion. Validation results for the
final implementation are recorded
in the repository's current validation report; this review does not claim a
complete PSI-MOD/OBO graph implementation or a full OpenMS port.

## Identification identity review

A separate consumer review found that protein modification observations,
sequence-only duplicate filtering, peptide identity keys and conflict resolution
could merge different custom records sharing the same displayed identifiers.
Those paths now compare complete chemical values, with equality-consistent
ordering for named records, anonymous annotations and sequences. Signed-zero
masses remain equal. Existing position/full-ID and conflict tie ordering is
retained before the new complete-value tie-break. Shared record allocations take
a constant-time self-comparison path.

`tests/custom_modification_identity.rs` checks same-name records with different
formulas and vocabulary accessions, retained protein observations, rank/spectrum
conflicts, sequence duplicate policy, owned keys after source drop and atomic
late errors. String-based keys in protein inference, FDR, explicit sequence-list
filters and peptide indexing are retained because those source methods explicitly
request textual or unmodified sequences. The associated source hashes are in the
registry provenance manifest.
