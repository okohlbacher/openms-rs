# Protein run metadata and export helpers

The native `ProteinIdentification` now includes the source inference-engine
queries/setters, original-engine recovery, run-level peptide-ID mergeability,
search-engine settings export, singleton-group completion and metadata-only
copying. These extend the existing [identification records](IDENTIFICATION_SUPPORT.md)
against Core SDK `82ce5b3`; the complete header remains under review.

`inference_engine()` and `inference_engine_version()` return borrowed strings.
Explicit search-parameter metadata wins, including an explicitly empty string.
An existing value of another type errors, preserving the source's strict
DataValue conversion. Without explicit metadata, Fido, BayesianProteinInference,
Epifany and ProteinInference count as inference engines. Percolator counts only
when an indistinguishable group exists. Version falls back to the search-engine
version only when inference data exists. An explicit version does not require
a valid engine value. Setters accept owned-convertible strings.

`original_search_engine_name()` returns the stored engine unless its name contains
`Percolator` or `ConsensusID`. Otherwise it returns the first `SE:` key suffix
whose full key does not contain lowercase `percolator`, or `Unknown` if absent.
Native metadata visits keys lexically; C++ visits numeric metadata-registry IDs.
This documented model difference can change the first winner in an ambiguous
multi-engine run and the order of metadata-derived export pairs. Names and case
tests otherwise follow the source exactly; metadata values are irrelevant here.

`peptide_ids_mergeable(other, experiment_type)` requires equal engine/version and
delegates to existing search-parameter mergeability. The `labeled_MS1` exception
for differing modification sets remains. Malformed native settings return an
error; source warning-stream output is represented by the boolean result.

`search_engine_settings_as_pairs(engine)` emits the twelve standard fields in
source order when the requested engine is empty, or matches a non-Percolator,
non-ConsensusID run. An empty request uses those fields even for merged runs.
Other requests select literal metadata prefixes and strip the prefix plus one
byte, normally a colon. As in source StringUtils::substr, a position past the
end is clamped; a boundary inside a UTF-8 scalar errors in the native API.
Numeric/list values use the already audited source formatting helpers, including
the upstream non-string mzTab export regression. Units are not appended.
Modification-list order and duplicates are retained; unrelated metadata is not
validated. The existing native four-value enzyme-specificity enum limits this
export; source no-N-terminal/no-C-terminal enum alternatives remain unrepresented.

`fill_indistinguishable_groups_with_singletons()` appends one group per previously
ungrouped accession, in hit order, and returns the number added. The first such
hit supplies its probability. Existing groups and quantity arrays remain intact.
Only newly consumed scores require finite values; a nonfinite score on an already
grouped hit is irrelevant. Duplicate and empty accessions follow ordinary exact
string identity. Staging preserves all group values when an error occurs.

`copy_metadata_only(source)` copies the run header, score configuration, search
parameters, date, both dedicated native run-path vectors and metadata. Source
paths reside inside its metadata, hence the two dedicated fields are copied too.
Target hits and both group vectors retain their original allocations and values.
Source protein results are not traversed, validated or cloned. Ordinary scalar
bits, including negative zero or unvalidated nonfinite state, are copied unchanged.

Operations creating owned output preflight a conservative one-million-descriptor
and 64 MiB logical payload allowance. Settings include cumulative formatting
scratch/output, metadata copies include sparse BTreeMap node capacity, and group
append includes final vector reallocation. Counts are conservative rather than a
promise to accept every object whose final payload is below 64 MiB. Existing
caller-owned group quantity arrays and source protein results are not copied.
Borrowed lookups and ordinary public field assignment retain normal Rust costs;
run mergeability retains the existing SearchParameters validation behavior.

[Tests](../tests/protein_run.rs) cover source branch tables, exact standard export
order, the original mzTab crash regression, all metadata value types, grouping,
ownership, UTF-8 boundaries and atomic resource errors. A separate private test
covers accounting overflow. [Provenance](../tests/data/protein_run_provenance.json)
pins the relevant source and helper files. No C++ execution is claimed for these
run helpers. Filesystem-aware primary-path selection and a full review of all
inherited ProteinIdentification/SearchParameters APIs remain open.
