# Identification groups

`openms::identification::graph` stores parent groups and groups of related
observation matches. These are grouping results and scored provenance; the
records do not run a protein-inference algorithm or calculate group scores.

The implementation follows `ParentGroup.h`, `ObservationMatchGroup.h` and
`IdentificationData.cpp` at Core SDK revision
`6bfc0e4711105f4eda2fea86812a83af7c7e791f`.

## Parent groups

`ParentGroup` contains `parent_refs: BTreeSet<ParentId>` and
`scores: BTreeMap<ScoreTypeId, f64>`. It has no metadata or leader field, matching
the source. Empty groups and mixed parent molecule types are valid; all referenced
parents and score types must already belong to the graph. Scores must be finite,
but negative values are allowed.

`ParentGroupSet` contains a label, `Vec<ParentGroup>`, and `ScoredProcessingResult`.
`register_parent_group_set` always appends a new set, including repeated labels,
identical contents and empty sets. Its returned `ParentGroupSetId` is a native
convenience over the source's append-only vector.

The incoming vector is normalized into parent-reference-set order. Equal parent
sets retain the first complete score map, including an empty first map; scores
are neither merged nor aggregated. Stable sorting preserves that first-value
rule. All incoming groups are validated and measured before duplicates are
removed. Draft vectors can therefore represent invalid duplicates that C++'s
multi-index container would already have discarded; registration rejects them.

```rust
use openms::identification::graph::{IdentificationData, ParentGroup, ParentGroupSet, ParentSequence};
use std::collections::BTreeSet;

let mut graph = IdentificationData::new()?;
let first = graph.register_parent_sequence(ParentSequence::new("protein_a"))?;
let second = graph.register_parent_sequence(ParentSequence::new("protein_b"))?;
let mut groups = ParentGroupSet::new("inference results");
groups.groups.push(ParentGroup::new(BTreeSet::from([first, second])));
let id = graph.register_parent_group_set(groups)?;
assert_eq!(graph.parent_group(id, 0)?.parent_refs.len(), 2);
# Ok::<(), openms::Error>(())
```

`parent_group_sets()` iterates in registration order, rather than sorting labels.
`parent_group_set(id)` and `parent_group(id, index)` return checked immutable
views. `parent_group_set_count()` counts sets and `parent_group_count()` counts
nested groups. The current processing step is attached to each newly registered
set only when it is not already present; existing history is not reordered.

## Observation-match groups

`ObservationMatchGroup` contains `observation_match_refs: BTreeSet<ObservationMatchId>`
and a scored result. `register_observation_match_group` returns a stable
`ObservationMatchGroupId`. Equal member sets deduplicate, merging scores/history
and overwriting incoming metadata keys without changing membership. Registration
also applies the current processing step, preserving history order.

`observation_match_group(id)` returns a checked view; `observation_match_groups()`
iterates by member-reference set. Full native equality includes metadata, whereas
the source's custom group equality ignores metadata.

`match_group_all_same_molecule(id)` compares exact molecule variants and IDs.
Different adducts or observations can still refer to the same molecule.
`match_group_all_same_query(id)` compares exact observation IDs, rather than
coordinates or data-ID strings. Both return true for empty and singleton groups.
The corresponding record methods `all_same_molecule` and `all_same_query` accept
resolving closures; their source-compatible empty/singleton path does not call
the resolver.

## Translation, limits and failures

Copy/merge translates parent-group score IDs, parent IDs, both group families'
processing history, and observation-match member IDs. Parent group sets remain
append-only during merge. Equal translated observation-match group keys use the
ordinary scored merge rule. Copies restore the translated source current step;
merges retain and apply the destination current step.

The source's `IdentificationData::merge`, and therefore its copy constructor,
omits observation-match groups. The native implementation deliberately preserves
and translates them to avoid silent data loss. `ReferenceTranslator` provides
`parent_group_set` and `observation_match_group` methods. Foreign or stale IDs are
errors, and `clear()` invalidates both new ID families.

Existing [graph limits](IDENTIFICATION_GRAPH_SUPPORT.md#checked-bounds) cover all
group records. Nested parent groups count toward the total record ceiling,
including empty groups. Member references, scores and history count toward the
edge ceiling. Payload validation, stable sorting, duplicate comparisons,
translations and snapshots share bounded work/allocation counters. Registration
stages only the affected records; whole-graph merge/copy uses the existing bounded
transaction. A late invalid reference, nonfinite score, or resource limit leaves
the original graph unchanged. No groups are silently truncated to fit limits.

## Evidence and remaining scope

[Focused tests](../tests/identification_groups.rs) preserve the two literal source
registration cases and independently check duplicate labels, first score-map
retention, history, empty/singleton comparisons, mixed types, foreign references,
complete translated copies, native metadata equality and atomic resource failures.
[Provenance](../tests/data/identification_groups_provenance.json) records all six
source hashes and distinguishes source assertions from derived checks and native
corrections. No C++ executable was built or run.

Referential deletion/cleanup, grouping algorithms, graph persistence and the
remaining legacy conversion operations are separate porting work.
