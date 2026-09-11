# Identification-to-run mapping

`identification::IdentifierMSRunMapper` implements all public operations from the
Core SDK `82ce5b3` header. It maps protein-identification identifiers to ordered
primary MS run paths and complete path lists back to identifiers. Native borrowed
lookups, `Option`, `Result`, and ordinary cloning replace C++ output arguments,
exceptions and copy operations. `identifiers()` visits keys in lexical order.

`create()` replaces both maps. Repeated identifiers overwrite the forward entry;
distinct path lists still retain their reverse entries. Duplicate path lists
error even when their identifiers match, and empty lists are ordinary keys.
The source intentionally leaves its complete forward map and only the reverse
prefix preceding the first duplicate available after this error. The native
implementation preserves that state, including forward records after the
collision. This supports source-file resolution after a reverse-map ambiguity.
Resource errors instead occur before copying and retain the previous mapping.

`primary_ms_run_path()` selects the integer `id_merge_index` from a nonempty mapped
run, defaulting to zero when absent. A positive index outside the list falls back
to peptide `base_name`, or the empty string. Negative and noninteger consumed
indices error before fallback. Unknown and explicitly empty runs ignore the index
and use the legacy fallback directly. `base_name` uses the existing source-style
lenient scalar/list formatting; units are omitted. Paths and string fallback
values are borrowed, while numeric/list fallback text is allocated under the
existing formatting limits.

`validate_merge_index()` is a separate operation. Only runs containing two or more
files require an explicit, integer, nonnegative index inside the list. Unknown,
empty and single-file runs are exempt, including stale or wrongly typed index
metadata. Bounded native diagnostics include the PSM number and file count;
they do not reproduce the source exception text verbatim.

Paths retain exact spelling, order, duplicates and Unicode. The mapper performs
no filesystem access or normalization. It consumes the native dedicated
`primary_ms_run_paths` field; raw-path fields, protein results and unrelated
metadata are not traversed. This follows the existing native identification
model, whose path fields replace source metadata-backed storage.

Construction checks a cumulative one-million run/path descriptor limit and
64 MiB conservative allocation allowance before cloning any retained value.
Accounting covers both maps, sparse map nodes, path-vector slots and copied
identifier/path bytes, including forward overwrites and the reverse prefix.
These are conservative construction limits, not a bound on caller-owned input
or ordinary Rust clone/destruction costs. The type is dependency-free.

[Ten tests](../tests/run_mapping.rs) cover the source class-test lookup literals,
merged-file selection, complete forward publication on duplicate error, repeated
identifiers, empty paths, strict consumed metadata conversions, lenient fallback
types, path identity, ownership and atomic resource failures. Source and native
control flow were independently reviewed; no C++ execution is claimed.
[Provenance](../tests/data/run_mapping_provenance.json) pins the complete header,
implementation, class tests and conversion helpers. This reviewed header does not
establish full ProteinIdentification or TOPP workflow parity.
