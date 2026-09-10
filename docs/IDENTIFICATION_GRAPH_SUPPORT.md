# Identification sequence and provenance graph

`openms::identification::graph` adds the sequence/provenance layer of the current
Core SDK's `IdentificationData`. It coexists with the legacy peptide/protein
identification records. This is a native owned graph, with no C++ references,
database dependency or general-purpose graph framework.

```rust
use openms::chemistry::RNaseDigestion;
use openms::identification::graph::{IdentificationData, MoleculeType, ParentSequence};

let mut graph = IdentificationData::new()?;
let mut parent = ParentSequence::new("rna_1");
parent.molecule_type = MoleculeType::RNA;
parent.sequence = "pAUGUCGCAG".into();
let id = graph.register_parent_sequence(parent)?;
RNaseDigestion::default().digest_identification_data(&mut graph)?;
graph.calculate_coverages(true)?;
assert_eq!(graph.oligo_count(), 3);
assert_eq!(graph.parent(id)?.coverage, 1.0);
# Ok::<(), openms::Error>(())
```

Run the [RNA identification example](../examples/identify_rna.rs) with
`cargo run --locked --example identify_rna` to print oligos, inclusive parent
positions, neighbor markers, processing history counts and parent coverage.

## Records and references

The implemented records are `InputFile`, `ScoreType`, `ProcessingSoftware`,
`DBSearchParam`, `ProcessingStep`, `AppliedProcessingStep`,
`ScoredProcessingResult`, `ParentSequence`, `ParentMatch`, `IdentifiedPeptide`
and `IdentifiedOligo`. Metadata remains typed `MetaInfo`. Peptide and RNA nodes
own the existing `AASequence` and `NASequence` values and their shared immutable
chemical records.

Registered records have immutable lookup views. The eight distinct ID types
carry a private graph owner and stable slot; a reference from another graph or
an earlier `clear()` fails validation. A move retains IDs. `new()` and `clear()`
return `Result` so graph-owner counter exhaustion cannot wrap. `Default` is a
convenience that panics only on exhaustion of the u64 owner space.

`parents()`, `peptides()`, `oligos()` and the other collection iterators yield
`(id, &record)` in semantic key order. Corresponding singular lookups return
`Result<&record>`, and each family has a count method. Full custom chemical
values determine sequence identity, including relevant RNA sulfur context;
equal display strings are not a deduplication key.

Source pointer ordering is replaced by deterministic owner/slot ordering where
references are keys. RNA enzyme records use complete value identity; protein
enzyme handles refer to immutable uniquely named entries in the native fixed
registry. ParentMatch equality/ordering deliberately compares positions only;
the first neighbors and metadata for an equal-position entry survive insertion.

## Registration and successful merge rules

| Family | Registration key | Existing-key behavior |
|---|---|---|
| Input file | Nonempty name | Fill blank design ID, reject conflicting nonblank IDs, union primary files |
| Score type | CV accession and required name | Reject opposite orientation; retain original CV and metadata payload |
| Software | Name and version | Retain original software metadata and assigned-score priority list |
| Search parameters | Explicit settings, excluding metadata | Retain the first payload |
| Processing step | Timestamp, software, ordered input files, actions | Retain first payload; first optional search association wins |
| Parent sequence | Nonempty accession | Merge history/metadata, fill blank sequence/description or reject conflicts, OR decoy flag; retain original type and coverage |
| Peptide/oligo | Complete nonempty chemical sequence | Merge history/metadata and union parent matches |

All referenced scores, inputs, software, steps and parents must already belong
to the same graph. Parent references are checked even when their position sets
are empty; peptide parents must be proteins and oligo parents must be RNA.
Parentless identified sequences are permitted. Parent registration stores text
without parsing it. Finite coverage must lie in [0,1]. All finite signed search
tolerances and charges are retained as record values; no search is performed by
registering them.

ProcessingStep's optional validated timestamp represents unset source time as
`None`. The C++ automatic current-time convenience is not reproduced. ScoreType
permits a name-only CV term, as required by the source tests; it does not invent
an accession. DBSearchParam is distinct from the legacy SearchParameters record.

## Scores, history and updates

Applied processing steps are unique by optional step ID, including a score-only
`None` entry. Re-adding a step overwrites its matching scores without moving its
history position. Score lookup scans history in reverse insertion order.
Software priority comes before remaining scores in native ID order; duplicate
priority entries are retained when requesting all ordered scores. Empty priority
and history results use `Option`, rather than a source NaN/success-flag pair.

`set_current_processing_step` affects both new registration and re-registration
of parents, peptides and oligos. It adds an empty step when needed; it does not
reorder an earlier occurrence. Clearing the current step does not remove history.

Graph `add_parent_score`, `add_peptide_score` and `add_oligo_score` update the last
inserted step (or the score-only None step). They do not implicitly apply the
current step. The corresponding `set_*_meta_value` methods overwrite one typed
metadata key. Re-registration remains available for larger result updates.
All changes are staged and checked before replacing stored records.

## Graph merging and copying

`merge_from(&other)` returns a `ReferenceTranslator`. Its typed lookup methods
reject missing translations; they never leak a foreign ID. References are
translated in dependency order, including applied scores, parent positions and
step/search associations.

The destination keeps its current step and graph-level metadata. Its current
step can therefore be attached to merged/re-registered sequence results. Unlike
ordinary step registration, graph merge overwrites a step's search association
with the translated incoming link.

`try_clone_with_translation()` creates a new owner, copies graph metadata,
translates records without adding an extra active step, then restores the
translated source current step. No derived public Clone can accidentally copy
graph references into another owner. `clear()` removes records, metadata and
current step while retaining configured native limits.

Merges and the RNase graph operation use one bounded staged snapshot for the
whole call. Ordinary registration copies only an affected existing record when
it needs merging. A late conflict or resource error leaves every prior stored
record and its ID valid. This corrects source multi-index failure paths that can
partially mutate or erase a record. There is no unchecked-registration mode.

## Coverage and RNA digestion

Parent matches use inclusive positions stored as `Option<usize>`; `None`
represents the source unknown-position sentinel. Registration preserves unknown,
reversed and out-of-range intervals. Coverage ignores invalid positions and can
optionally require the identified sequence length to match the interval length.

`calculate_coverages(check_molecule_length)` parses referenced parent text and
unions valid intervals. Peptide and RNA parent lengths count residues, not text
bytes or bracket characters. `calculate_coverages_with_registries` supports
caller-owned peptide and RNA chemistry. Every parent's coverage is overwritten,
including zero for parents with no contributing intervals. No partial coverage
update survives a later parent parse failure.

The source has an observable empty-parent branch: when a molecule reaches an
empty-sequence parent, it breaks out of that molecule's parent loop. The native
graph retains this branch in deterministic parent-ID order. It does not silently
substitute a continue. Interval sorting/union replaces the source boolean array
without changing covered residue counts.

RNase graph digestion visits parents in accession order, skips non-RNA parents,
and appends/merges oligos with inclusive positions and source neighbor markers.
It retains the current processing step. See [RNA digestion](RNASE_SUPPORT.md)
for the registry-aware API, exact neighbor convention and shared batch bounds.

## Checked bounds

`GraphLimits` defaults and fixed maximum ceilings are 100,000 records summed
across families, one million graph edges, 256 MiB retained logical payload and
50 million work units per operation. Limits may be lowered, including to zero.
Each record also has a 64 MiB logical payload ceiling. Metadata keys/values,
CV terms, lists, history, matches and complete chemical comparison inputs are
measured before copying or comparing them. Shared chemistry may be counted many
times; these are conservative allowances, not resident-memory measurements.

Stable slot storage uses separately sorted slot indexes. Binary-search comparisons
are charged conservatively using complete payload sizes; sorted-index insertion
charges its linear movement. Cumulative allocation accounting shares the byte
ceiling across temporary maps, vector growth, snapshots and translated records.
Small operations on large records can consequently fail before the retained
graph reaches its nominal byte ceiling.

One private batch transaction shares its remaining work across every nested
registration. A new allowance is not granted for each RNA parent or product.
Coverage passes the same work and allocation counters into parent parsers.
Modified peptide parent parsing uses a conservative registry/annotation-dependent
preflight; plain parents avoid scanning the registry. RNA parsing charges residue
storage and payload incrementally before each reserve/copy, sharing the same
remaining counters; a zero byte allowance fails before sequence allocation.
No operation truncates data to fit a limit.

## Evidence and remaining graph layers

The source is `METADATA/ID` at Core revision
`6bfc0e4711105f4eda2fea86812a83af7c7e791f`; these files are byte-identical to the
historical `7c029e8` graph. [Direct graph tests](../tests/identification_graph.rs)
exercise registration, history, full chemical identity, translated copy/merge,
coverage and atomic failures. A private transaction regression distinguishes a
shared batch work allowance from independent single registrations. Independent
source references are in [identification_graph_reference.rs](../tests/identification_graph_reference.rs).
The [RNA workflow](../tests/rna_identification_workflow.rs) combines shared-parent
digestion, processing history, translated graph copy and generated spectra. It
also checks custom chemistry survival, raw-code neighbors and atomic product,
residue, parsing-byte and non-ASCII-neighbor failures.
No C++ build or execution is used for these references.

Observations, compounds, adduct registration, observation matches and groups,
parent groups, best-match queries, deletion/full cleanup, graph persistence and
IdentificationDataConverter remain unimplemented. No placeholder cleanup method
erases sequence results, and this sequence/provenance layer does not claim the
complete source IdentificationData API.
