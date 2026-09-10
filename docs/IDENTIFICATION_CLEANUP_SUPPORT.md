# Identification graph filtering and cleanup

`IdentificationData` implements the complete five-option `cleanup` operation and
both predicate-based removal helpers from the reduced Core SDK at `6bfc0e4`.
These operations remove records and repair references; they do not calculate new
identification scores or infer groups.

## API

- `cleanup(CleanupOptions) -> Result<CleanupReport>` performs the ordered cascade.
- `remove_observation_matches_if(|id, record| bool)` and
  `remove_parent_sequences_if(|id, record| bool)` remove selected records and run
  **default** cleanup only when the predicate selected at least one record.
- `try_remove_observation_matches_if` and `try_remove_parent_sequences_if` accept
  the same predicates returning `Result<bool>`. An error rolls back all graph
  changes. Predicate visits follow the graph's semantic iteration order; effects
  the caller performs outside the graph are not rolled back.

The option names and defaults match the source:

| Option | Default | Removal criterion |
| --- | --- | --- |
| `require_observation_match` | true | Molecules, observations and adducts lacking a surviving observation match |
| `require_identified_sequence` | true | Parents lacking a surviving identified peptide/oligo reference |
| `require_parent_match` | true | Peptides and oligos with no parent-reference entry |
| `require_parent_group` | false | Parents absent from every parent group |
| `require_match_group` | false | Observation matches absent from every match group |

An empty set of positions under a valid parent reference still counts as a parent
match. Compounds are not subject to the parent-match requirement. The options are
independent: all 32 combinations are tested.

## Source cascade and retained information

Cleanup runs once, in the source's order:

1. Optionally remove parents outside parent groups, then remove invalid parent
   links from all peptides and oligos.
2. Optionally remove parentless peptides and oligos. Always remove observation
   matches whose identified molecule was removed.
3. Optionally remove matches outside match groups.
4. Optionally remove molecules, observations and adducts unused by the surviving
   observation matches.
5. Optionally remove parents unused by the surviving sequences.
6. Repair both group kinds and remove empty groups. Parent grouping operations
   themselves remain, even when all their groups are removed.

Provenance records, input files, score definitions, graph metadata, processing
history and the current processing step remain. Cleanup does not recalculate
coverage or rescore groups. `CleanupReport` returns counts for removed records,
parent links and groups. Its `parent_group_scores_may_be_invalid` and
`match_group_scores_may_be_invalid` flags expose the source warnings when a
nonempty group's membership shrinks; existing scores are retained.

## Stable references and checked differences

Records are held in sparse standard-library maps, with a separate index for the
source's semantic key order. A successful deletion frees its record storage.
Surviving IDs retain their exact values. Removed IDs are rejected permanently,
even if an equal-key record is subsequently registered. Slot numbers advance
monotonically and are checked before allocation; there is no lifetime record
count ceiling or growing collection of deleted records. `clear` invalidates all
IDs by allocating a new graph identity. Copy/merge translates only live records.
Nested parent groups are addressed by their current position within a grouping
operation; those positional indices can change when preceding groups are removed.

Both source group-repair loops use Boost `multi_index::modify` and then dereference
the same iterator. When shrinking two keys makes them equal, modification can
erase the current element, making that dereference undefined. The native code
instead sorts the repaired keys stably, retaining the payload of the first group
in the **original semantic key order**, and removes later collisions. A removed
observation-match-group ID becomes invalid. No payload winner or numerical result
is claimed as source parity for that undefined case. Existing group payloads are
not merged during cleanup.

A complete cleanup/removal is staged in one bounded graph snapshot and committed
only after all checks pass. A predicate that selects nothing returns immediately
without a snapshot or cleanup. The graph's existing record, edge, byte and work
limits remain shared across the operation. Sparse table nodes, temporary reference
sets, group sorting, comparisons, payload destruction and snapshot copies are
charged conservatively. These limits describe retained logical payload and
cumulative operation allocation/work; they are not a process RSS guarantee.
Deletion releases retained record and edge allowance for later registrations.
No new dependency or source build is required.

## Evidence

[Focused tests](../tests/identification_cleanup.rs) preserve the literal source
cleanup counts (peptides 4→3→1; oligos 2→1→1), exercise all option combinations,
parent/match cascades, sparse-ID coverage/copy/merge, unchanged-history behavior,
key collisions, no-op filtering, repeated deletion/reinsertion and atomic error
paths. The other graph, group, match and converter suites protect existing APIs
through the sparse storage change.

[Provenance](../tests/data/identification_cleanup_provenance.json) records source
hashes and distinguishes source assertions, independently derived branch cases,
and native safety extensions. Persistence and the rest of the legacy converter
surface are separate remaining work; the cleanup surface itself has no stubbed
options.
