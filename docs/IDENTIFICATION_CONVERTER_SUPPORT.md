# Identification sequence and evidence conversion

`identification::graph::IdentificationDataConverter` now provides the source
`importSequences` and `exportParentMatches` operations. It connects existing FASTA
records to graph parents and graph parent matches to legacy peptide evidence.
Both operations use native graph limits and commit atomically.

```rust
use openms::identification::graph::{
    IdentificationData, IdentificationDataConverter, MoleculeType,
};
use openms::format::fasta::FASTAEntry;

let mut graph = IdentificationData::new()?;
let entries = [FASTAEntry {
    identifier: "protein_1".into(),
    description: "example".into(),
    sequence: "PEPTIDE".into(),
}];
IdentificationDataConverter::import_sequences(
    &mut graph, &entries, MoleculeType::Protein, "DECOY_",
)?;
# Ok::<(), openms::Error>(())
```

`import_sequences` preserves identifier, description and sequence text exactly;
it does not parse peptide/RNA syntax or normalize case. The explicit Rust
arguments correspond to source defaults `Protein` and an empty decoy pattern.
A nonempty pattern matches any case-sensitive substring of an accession, rather
than only a prefix/suffix. All three molecule-type enum values are retained,
as in the source. Empty input is a no-op.

Parent registration determines repeated-accession behavior: existing nonempty
sequence/description must agree; missing fields can be filled; decoy flags combine
with logical OR; the first molecule type and coverage remain. Current processing
step is recorded through the normal registration path. A conflict or limit reached
late in an import rolls back the whole batch, preserving preexisting IDs.

`export_parent_matches(&graph, &matches, &mut hit)` appends to `hit.evidences`
and sorts the complete result by accession, inclusive start, inclusive end,
preceding character and following character. Duplicate evidences remain.
Unknown positions use `None` and sort before known positions, corresponding to
the source legacy `-1`. Empty flank strings become `X`; a nonempty flank uses
only its first byte. The remaining flank text and parent-match metadata are not
exported. Even an empty match map sorts preexisting evidences.

The exporter preserves the rest of the hit without cloning or validating unrelated
scores, annotations or sequence chemistry. It validates all supplied parent IDs,
including keys whose match set is empty. Invalid intervals, positions exceeding
legacy signed 32-bit range, or flank bytes that cannot be represented as an
uppercase amino acid or `X`/`[`/`]` marker fail before changing evidence. These
checks replace source narrowing or invalid legacy states; a modified RNA flank
such as `m6A` is therefore a checked error instead of an invented amino-acid code.

Each import shares the graph's complete transaction work/allocation limits.
Export checks combined existing/new evidence against `max_edges`, precharges
accession clones, vector storage and sorting work, and builds a replacement
before assignment. No new dependency is used.

[Focused tests](../tests/identification_converter.rs) reproduce the source's
five-parent FASTA assertion and independently check merge, decoy, append/sort,
reference and rollback rules. [Provenance](../tests/data/identification_converter_provenance.json)
records nine source hashes and the independent five-row fixture projection.
The original bridge test constructs `FASTAEntry` values directly from its TSV.
The [native FASTA reader](FASTA_SUPPORT.md) now also checks the unmodified source
file, including its PEFF prologue, annotated sequence and preserved description
spaces. Neither set of tests executes the C++ runtime.

Complete `importIDs`/`exportIDs` remains separate work, including run grouping,
score priority, metadata conventions and protein-group conversion. mzTab export
needs the mzTab data model, while feature/consensus converters need graph links
inside those map types. Graph persistence and referential cleanup are also not
implemented by this bridge. The existing [idXML adapter](IDXML_SUPPORT.md) operates
on legacy identification records.

Source: [IdentificationDataConverter.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/6bfc0e4711105f4eda2fea86812a83af7c7e791f/src/openms/source/METADATA/ID/IdentificationDataConverter.cpp#L764)
and [PeptideEvidence.cpp](https://github.com/okohlbacher/OpenMS4-core/blob/6bfc0e4711105f4eda2fea86812a83af7c7e791f/src/openms/source/METADATA/PeptideEvidence.cpp#L48).
