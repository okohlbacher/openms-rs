# Sequence tags from measured spectra

`chemistry::Tagger` ports the complete pinned Tagger surface: construction with
length, tolerance, charge and modification settings; extraction from m/z slices
or spectra; appending to an existing tag collection; and changing maximum charge.
It returns sorted, unique strings of parent amino-acid letters.

```rust
use openms::chemistry::{Tagger, TaggerOptions};
use openms::comparison::Tolerance;

let mut options = TaggerOptions::new(2, Tolerance::Absolute(0.02));
options.max_tag_length = 3;
let tagger = Tagger::new(options)?;
let tags = tagger.get_tags(&[150.0, 247.0527, 376.0953, 473.1480])?;
assert_eq!(tags, ["EP", "PE", "PEP"]);
```

The [executable example](../examples/extract_tags.rs) extracts six tags from a
small measured mass ladder. [Workflow tests](../tests/tagger_workflow.rs) use
tags for exact peptide/protein substring matching, compare a reversed decoy,
recover modified residues from digest fragments and preserve the tag set through
both mzML compression modes. Tagger itself does not score candidates or identify
peptides; database search remains outside this helper.

## Options and matching

`TaggerOptions::new(min_tag_length, tolerance)` defaults to maximum length 65535,
charges 1 through 1 and no fixed or variable modifications. Use
`Tolerance::Absolute(da)` or `Tolerance::Ppm(ppm)` explicitly. As in C++, finite
negative tolerance is normalized to its absolute value. Zero tolerance matches
no residues, even at exact mass equality.

For each starting peak and charge, Tagger explores later peak indices whose
mass difference times that charge matches a residue mass. Ppm tolerance is
computed as `(ppm / 1_000_000) * charged_gap`; it is relative to the gap rather
than either peak's m/z. The closest residue wins; equal distances retain the
lower mass. The tolerance is strict. A first candidate exactly at the lower
boundary rejects the query even when a later candidate lies inside, preserving
the source lookup behavior.

The base dictionary contains the 19 standard residues excluding I. An L match
branches into both L and I tags, including modified L. A separately installed I
mass emits only I. No marker distinguishes a modified residue in the returned
text. Exact mass collisions retain the most recently assigned residue letter.

Only m/z values are consumed. Spectrum intensities, precursors, acquisition
metadata and auxiliary arrays are ignored and remain untouched. Input index
order is preserved: unsorted finite coordinates are accepted, with the source's
early branch termination when a mass gap exceeds the largest allowed gap.
Signed and duplicate coordinates are accepted. Sorting an unsorted input can
therefore change the results; ordinary measured spectra should already be in
ascending m/z order.

Minimum and maximum lengths count residue edges, not peaks. Zero or inverted
length and charge bounds retain finite source behavior: zero minimum never
emits an empty tag, zero maximum prevents extension, and an inverted range
emits no new tags. Zero charge is accepted; it generally matches nothing but
can match with a sufficiently large absolute tolerance.

## Modification resolution

Use `fixed_mods` and `variable_mods` for registry names, full IDs or aliases.
`with_registry(options, &registry)` supports caller-owned records. Construction
stores owned scalar masses; registry lifetimes or subsequent changes cannot
alter a constructed Tagger.

The source resolves each request twice: first without residue/terminal filters,
then again using the selected short ID, selected origin and Anywhere
specificity. A terminal request can consequently succeed through an Anywhere
record with the same short ID. Its initial terminal chemistry is not applied.
Anonymous records with empty short IDs use a private matching path; the general
registry's public empty-name and ambiguity policies are unchanged.

C++ orders matching records by pointer value, and its fast second lookup returns
the last matching pointer despite a warning that says first. The native helper
uses the first matching record in stored provider order for both stages. This
is a deterministic replacement for allocation-dependent behavior; ambiguous
names are not guaranteed to choose the same record as a particular C++ process.

All fixed requests run in order, followed by variable requests. A fixed request
removes the first mass-sorted entry with its parent letter and assigns a mass
calculated from the fresh base residue. Repeated fixed requests do not accumulate
deltas. Variable requests add or overwrite mass keys without removing the base.
Fixed I does not remove unmodified L.

Masses preserve `Residue` operation order: calculate stored free-residue mass,
apply the modification's formula/absolute/delta precedence, then subtract water.
This can differ by a rounding bit from evaluating an internal formula directly.
The private scalar dictionary retains finite signed and zero masses. In
particular, source B/Z/X placeholders have free mass zero before water
subtraction. This does not supply masses for unresolved `AASequence` residues.

## Append behavior and checked limits

`get_tags` and `get_spectrum_tags` return a new collection. `append_tags` and
`append_spectrum_tags` sort and deduplicate the combined old and new strings;
old labels need not be valid amino-acid tags. The source early return when
minimum length exceeds the number of peaks leaves existing strings completely
untouched. Equal minimum length takes the normal sorting/deduplication path.
Every checked error also leaves existing output unchanged.

After the minimum-length early-return guard, each extraction permits at most
1,000,000 input peaks, 1,000,000 existing plus
emitted tags before deduplication, 64 MiB of existing plus emitted string data,
50,000,000 traversal/lookup/sort work units and 256 MiB of conservative staged
allocation accounting. Iterative traversal avoids recursive stack overflow.
Sorting charges actual string-byte comparisons; repeated paths consume budget
even when they later deduplicate. Existing strings are moved only after all
fallible work completes.

Inputs skipped by that early return are not validated, including unused
nonfinite coordinates or a peak count above the ordinary extraction limit.

Construction permits 10,000 total modification requests, 1 MiB of query bytes
including the second lookups, 200,000 registry records when modifications are
requested, and 50,000,000 conservative lookup/formula/table work units. Nonfinite
coordinates or arithmetic, unknown modification origins and missing resolution
return errors. Maximum charge `usize::MAX` is rejected because the source's
inclusive loop would overflow; large finite ranges remain subject to work limits.

The [independent review](TAGGER_REFERENCE_REVIEW.md) and
[source/fixture provenance](../tests/data/tagger_provenance.json) distinguish
literal C++ count/membership assertions from independently derived boundary
cases. The C++ reference was inspected without building or executing it.
