# Identification graph records and score history

`identification::graph` represents provenance-aware peptide and RNA identifications
separately from the existing legacy `PeptideIdentification` and
`ProteinIdentification` APIs. Its record and score-history layer follows
OpenMS4-core [`6bfc0e4`](https://github.com/okohlbacher/OpenMS4-core/tree/6bfc0e4711105f4eda2fea86812a83af7c7e791f/src/openms/include/OpenMS/METADATA/ID).
These scientific source files are unchanged from the original `7c029e8` fixture
revision. This port uses owned native values and graph-specific IDs; it does not
store source iterators or raw pointers.

The implemented records are `InputFile`, `ScoreType`, `ProcessingSoftware`,
`DBSearchParam`, `ProcessingStep`, `AppliedProcessingStep`,
`ScoredProcessingResult`, `ParentSequence`, `ParentMatch`, `IdentifiedPeptide`
and `IdentifiedOligo`. `MoleculeType` includes Protein, Compound and RNA, and
`MassType` includes Monoisotopic and Average. The compound enum value is useful
for parent/search records; it does not imply that compound graph nodes have
already been implemented. Observations, identified compounds, match/group nodes
and the complete legacy converter remain separate graph work.

## Construction and ownership

Records have public editable fields. They are checked when passed to graph
registration or the checked record merge/update methods; immutable graph lookups
preserve registered keys and ownership. Empty draft names are representable,
but input-file names, score names and parent accessions must be nonempty when
registered. Identified peptide/oligo sequences must be nonempty. Parent sequence
text and descriptions may be empty; parent text is parsed only by operations
that consume it, such as RNA digestion or coverage calculation.

IDs are distinct types (`InputFileId`, `ScoreTypeId`, `ProcessingSoftwareId`,
`SearchParamId`, `ProcessingStepId`, `ParentId`, `PeptideId`, `OligoId`). Their
read-only `owner()` and `slot()` accessors expose opaque identity; external code
cannot construct or mutate them. The graph validates ownership and registration,
including references attached to empty parent-match sets. Record-level score
queries operate on supplied IDs and do not independently resolve their graph.

```rust
use openms::identification::graph::{
    IdentificationData, InputFile, ProcessingSoftware, ProcessingStep,
    ScoreType, ScoredProcessingResult,
};

let mut graph = IdentificationData::new()?;
let input = graph.register_input_file(InputFile::new("sample.mzML"))?;
let score = graph.register_score_type(ScoreType::new("search score", true))?;
let mut software = ProcessingSoftware::new("SearchTool", "1.0");
software.assigned_scores.push(score);
let software = graph.register_processing_software(software)?;
let mut step = ProcessingStep::new(software);
step.input_files.push(input);
step.date_time = Some("2026-09-10 12:00:00".parse()?);
let step = graph.register_processing_step(step, None)?;
let mut result = ScoredProcessingResult::default();
result.add_score(score, 42.0, Some(step))?;
assert_eq!(result.score(score), Some(42.0));
Ok::<(), openms::Error>(())
```

`ProcessingSoftware` reuses `metadata::Software`, including CV metadata, and adds
its ordered assigned-score list. `ProcessingStep` keeps ordered input-file IDs,
an action set, metadata, and an optional validated `CompletionTime`. `None`
explicitly represents an absent timestamp. Construction does not obtain the
machine's current time as the C++ convenience constructor does; callers provide
a timestamp when relevant.

`DBSearchParam` preserves the complete source search settings: molecule/mass
type, database/version/taxonomy, charge set, fixed/variable modification sets,
both mass tolerances and unit flags, enzyme, specificity, missed cleavages and
length limits. `None` specificity means source Unknown; `Some(DigestionSpecificity::None)` means nonspecific digestion. Tolerances must
be finite; signed finite values and otherwise unusual combinations are retained
because this record does not perform the search. Positive and negative floating
zero have the same registration key.

`GraphEnzyme` owns a pinned `DigestionEnzymeProtein` value or an
`Arc<DigestionEnzymeRNA>`. Protein names uniquely identify the immutable pinned
records; RNA keys include every record field, including synonyms and cleavage
and terminal-gain strings. Separately allocated equal RNA records have equal
keys; same-name records with different chemistry remain distinct. This replaces
allocation-dependent C++ enzyme-pointer identity with deterministic value identity.

## Registration keys and successful merge behavior

Registration keys are explicit, not derived from each record's full payload.
Most Rust records use full-value `PartialEq`, which is useful for detecting
changes; it must not be confused with the following registration keys.

| Record | Key | On a repeated key |
|---|---|---|
| InputFile | Name | Fill blank experimental-design ID, reject conflicting nonblank IDs, union primary files |
| ScoreType | CV accession and name | Require the same higher-better flag; retain the first other CV/metadata payload |
| ProcessingSoftware | Software name and version | Retain first software metadata and score priorities |
| DBSearchParam | Every search setting except metadata | Retain first metadata |
| ProcessingStep | Timestamp, software ID, ordered inputs, actions | Retain first metadata; search-link handling belongs to graph registration |
| ParentSequence | Accession | Merge history/metadata, fill blank sequence/description or reject conflicts, OR decoy flag, retain original coverage and molecule type |
| IdentifiedPeptide/Oligo | Complete native sequence value | Merge history/metadata and union parent-position matches |

Custom chemical annotations and RNA slicing context remain part of sequence
identity. Display strings are not keys. Merging a record leaves its own key
unchanged; the graph calls merges only after matching registration keys.

Merge and history-update errors are atomic in Rust: validation, work/payload
limits and conflicts are checked before committing a replacement. The C++
callback/container behavior can leave partial mutations after an exception;
that failure behavior is deliberately not reproduced.

## Score history

`ScoredProcessingResult::steps_and_scores` is an insertion-ordered vector unique
by optional processing-step ID. The `None` entry holds scores without provenance.
`add_processing_step`, `add_score` and `merge` update an existing step in place;
updating a step does not move it to the end. Incoming scores replace matching
values, and incoming metadata overwrites matching keys. Registration rejects
manually constructed vectors containing duplicate step IDs.

`score`, `score_at_step` and `score_and_step` return `Option` rather than the
source's NaN/success tuple. Cross-step lookup examines reverse insertion order,
not timestamps or the most recently updated entry. `steps_by_processing_step`
returns a borrowed view sorted by optional step ID without changing recency.
`number_of_scores` counts
score entries across all steps, including repeated types in different steps.
`most_recent_score` takes a resolver closure returning a borrowed software
priority list for a processing-step ID and selects the first score from the
latest step that has scores.

`AppliedProcessingStep::scores_in_order(priorities, primary_only)` keeps the
source's duplicate priority behavior: if a software priority list repeats an
assigned score, the complete result repeats it too. Unlisted scores follow in
deterministic typed-ID order. A score-only `None` step ignores supplied software
priorities. An empty score map returns an empty result. The fallback ID order
replaces source object-address order; no pointer-order compatibility is claimed.

A score CV may have a name without an accession (`ScoreType::new`). Its scoped
validator checks finite values/units and allows that source case without changing
`CVTerm::validate()` elsewhere. Nonempty graph CV accessions retain the existing
native whitespace restriction. Score values and floating metadata must be finite;
parent coverage must be finite and in `[0,1]`.

## Parent matches

`ParentMatch` uses optional **inclusive** start/end positions. `None` represents
the source unknown-position sentinel. Its equality and ordering compare only
positions; neighbors and metadata do not participate. Inserting an equal-position
match retains the first payload. Known positions sort before unknown ones.
`has_valid_positions(molecule_length, parent_length)` rejects missing, reversed
or overflowing intervals. A zero requested length skips that length check;
otherwise molecule length must match the inclusive span and end must lie inside
the parent. Registration can retain positions that are invalid for coverage,
matching source behavior. Unlike the C++ sentinel, an explicit
`Some(usize::MAX)` remains a supplied position: `[MAX,MAX]` has checked span one
when both length checks are disabled. `None` is the only unknown value. Graph
coverage supplies the actual parent length and rejects such out-of-range positions.

Neighbors are strings, with source constants `X` (unknown), `[` (left terminus)
and `]` (right terminus). `all_parents_are_decoys` requires at least one parent
and accepts a resolver closure; graph convenience methods validate and resolve
parent IDs. A parent entry with no positions still counts as a parent.

## Bounds and validation

One record is limited to 64 MiB of conservative logical payload, counting nested
metadata, CV values/units, score/history/edge storage and complete shared chemical
records. Traversal is charged before visiting arbitrary-size lists. Counting
chemical records per occurrence intentionally overestimates shared Arc memory,
but also bounds complete-value comparisons and same-ID custom-record keys.

Graph operations share the graph's configured work and allocation budgets;
standalone checked record operations use the default ceilings of 50 million work
units and 256 MiB of cumulative staging allocations. Payload bytes serve as a
conservative work estimate. Merge planning additionally charges map comparisons
and history scans before cloning, so an oversized or unusually expensive merge
can return a limit error before the nominal per-record byte ceiling. Returned
priority vectors and temporary history-index sets are bounded as well. These are
checked native resource limits, not OpenMS scientific constants.

The direct [record tests](../tests/identification_graph_records.rs) cover source
field/key behavior, CV name-only scores, priority duplicates, history recency,
atomic conflicts, signed zeros, custom chemistry, inclusive/unknown positions,
and work-limit errors. Graph registration, translation, coverage and RNA digestion
have separate integration/reference tests. Source contracts come from the
`METADATA/ID` headers and `IdentificationData_test.cpp`; no C++ build is used to
manufacture expected results.
