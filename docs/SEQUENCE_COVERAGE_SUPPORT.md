# Sequence coverage

`SequenceCoverage::get_coverage(&protein, &peptides) -> Result<f64>` implements
the complete public scientific operation of OpenMS Core SDK
`CHEMISTRY/SequenceCoverage.h` at revision
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. Inputs are native `AASequence` values;
the result is a percentage from 0 to 100, rather than a fraction.

```rust
use openms::chemistry::{AASequence, SequenceCoverage};

let protein = AASequence::parse("ACDEFGHIK")?;
let peptides = [AASequence::parse("ACD")?, AASequence::parse("FGH")?];
let percent = SequenceCoverage::get_coverage(&protein, &peptides)?;
assert_eq!(percent, 6.0 * 100.0 / 9.0);
# Ok::<(), openms::Error>(())
```

The example is explanatory Rust, not an additional crate doctest.

## Matching behavior

Every exact occurrence of every nonempty peptide contributes coverage. Matches
may overlap, including multiple occurrences of the same peptide: `AAA` covers
all five positions of `AAAAA`. The result counts the union of covered positions;
duplicated peptides, repeated matches and peptide order do not increase the
percentage beyond that union. Missing peptides and peptides longer than the
protein contribute nothing. Empty protein or an empty peptide list returns zero
immediately; individual empty peptides are skipped.

Matching uses the unmodified uppercase residue bytes from `AASequence::as_str()`.
Known modifications, anonymous mass tags and terminal annotations on either side
are ignored, just as the source calls `toUnmodifiedString()`. No formula or mass
availability is required. Symbols such as B, Z and X match themselves literally;
X is not a wildcard, and isoleucine and leucine are not treated as interchangeable.
Parsing and supported residue symbols follow the existing native
[AASequence contract](SEQUENCE_SUPPORT.md).

The final calculation retains source operation order:
`covered_count as f64 * 100.0 / protein_length as f64`. The native size limits
keep the integer count exactly representable. The utility neither digests a
protein nor consults peptide evidence coordinates, scores, registries or matches
in other proteins.

## Resource and failure behavior

The implementation borrows sequence strings and ignores their chemical metadata;
it does not clone or traverse annotation payloads. It uses one flag per protein
position and scans possible byte windows, avoiding an additional search framework.
The source's allocated unmodified strings and packed `vector<bool>` are replaced
by borrowed bytes and native boolean storage.

For a nonempty protein and peptide list, fixed checked limits are:

| Limit | Maximum |
|---|---:|
| Protein residues | 1,000,000 |
| Peptide entries | 1,000,000 |
| Sum of supplied peptide lengths, including duplicates and longer peptides | 1,000,000 |
| Cumulative work units | 50,000,000 |
| Requested coverage flags | One byte per protein residue, at most 1,000,000 bytes |

The empty-protein/list shortcuts happen before these limits and require no
allocation or traversal of unused inputs. For a substantive query, both peptide
passes, flag initialization and final coverage sum are charged before allocation.
Every eligible peptide additionally precharges
`(protein_length - peptide_length + 1) * (peptide_length + 1)` work: one visit
plus the full possible byte comparison at every position. Checked arithmetic
prevents overflow before allocating the flags. Successful matches then charge
their full span before marking it, even if earlier matches already covered it.

This comparison bound is **conservative**, not a count of actual byte comparisons.
A mismatching first byte may end equality testing immediately, while the entire
possible comparison length remains charged. Consequently a source query that
would run quickly can return a native resource error. For example, two copies of
a 5,000-residue `G` peptide against 10,000 `A` residues exceed the comparison
budget despite having no matches. The bound is explicit and reproducible; no
claim is made that the source itself imposes these limits. Requested flag storage
is a logical bound, not an exact allocator-overhead or process-memory ceiling.

Allocation uses a fallible reservation. Size, arithmetic, work or reservation
failure returns an error without changing any input. A late marking-budget failure
drops the local coverage flags and publishes no partial percentage. There are no
mutable options, global caches, output buffers or hidden modification lookups.

## Evidence

The [source-hash manifest](../tests/data/sequence_coverage_provenance.json) records
the pinned public header, implementation, complete class test and the consulted
AASequence unmodified-string implementation. All four source test cases are
retained: six of nine residues (approximately 66.6667%), an empty peptide list,
an empty protein, and overlapping `ABC`/`BCD` covering four of five positions.

Independent tests cover all binary A/G proteins through length seven with varied
peptides using a per-position interval oracle, overlapping repeat matches,
duplicates, annotation independence, literal unknown residues, no matches,
longer peptides, input ownership and conservative resource rejection. Private
tests check exact cumulative search/marking costs, a late failure after one
marked span, and the zero-work empty shortcuts. No C++ execution is used.
