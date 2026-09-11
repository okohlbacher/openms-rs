# ControlledVocabulary support

`format::controlled_vocabulary` implements the complete public operation group of OpenMS4-core `ControlledVocabulary.h` at `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. It requires no optional feature or additional dependency. `CVTermDefinition` is an ontology definition; the existing `metadata::CVTerm` remains a parameter instance.

The owned definition retains every source field: name, accession, description, parent/child IDs, obsolete flag, ordered synonyms, unparsed lines, all eleven XRefType alternatives, binary types and allowed units. Default/copy/assignment correspond to `Default`, `Clone` and Rust moves. Registry fields are private so borrowed lookup results cannot invalidate indexes.

The registry exposes exact accession/name lookup, the full ordered term map, parent/descendant queries, source-order early-stop callbacks, score direction, XML parameter rendering, FNV-1a hashing, diagnostic stream rendering, bounded OBO file/reader loading and the complete five-provider `psi_ms()` singleton. `first_child_with_name` is a direct convenience for the source mzML writer's descendant-name lookup. `psi_ms()` initializes once without installation paths or network access, and returns a recoverable error if its checked initialization fails.

## Loading and query behavior

The OBO reader implements every branch consumed by the source provider, preserving other in-term lines in `unparsed`. It is not an OWL reasoner or a complete OBO language interpreter: imports, Typedef semantics and unconsumed relationships do not acquire new behavior. Imports do not trigger network downloads.

Successful loads accumulate into the existing registry. Duplicate accessions replace definitions. Missing metadata headers retain previous values. Source child/name indexes are additive: old child links and aliases can survive a replacement, and a reused name resolves to the first source index entry. The optional description is tried only after an unsuccessful plain-name lookup, through literal `name + description` concatenation. Synonyms are retained but are not aliases.

Missing parent accessions create empty placeholder definitions. Their map key can differ from their empty `id` field. The provider preserves source ordered-map iteration, including whether a newly inserted placeholder is visited in the current pass. Such behavior is tested independently and is not silently replaced by a cleaned graph rebuild.

Parent/child sets iterate in lexical accession order. Descendant callbacks visit a shared DAG node once per source path, and stop immediately when the callback returns true. The root is excluded from ordinary descendant traversal; `add_all_child_terms` adds it explicitly. Acyclic self-is-child queries return false. An unknown queried parent may return false, while an unknown child errors, matching the source implementation rather than its broader exception comment.

Traversal is iterative. A traversed active-path cycle returns an error instead of recursive nontermination. Early success before reaching a later cycle still succeeds. Caller-set extensions and loads commit atomically; callback effects already performed before a later error cannot be rolled back.

Source trimming/removal recognizes only ASCII space, tab, CR and LF. Bracketed term stanza matching and warning-name comparisons are ASCII-insensitive; no Unicode case folding is added. Definitions/synonyms retain the source first-quote/second-quote slicing rather than applying a new OBO escape decoder. Known source oddities such as `X\!Tandem:hyperscore` therefore remain intact.

The extra BRENDA `DRV`/`part_of` parent interpretation runs only for the exact load name `brenda`. The source singleton calls this provider `BTO`, so those lines stay unparsed there. All eleven value types and their source aliases, xref/xref_analog forms, newer has_value_type relationships, list-type IDs, has_units and unknown-type diagnostics are covered. Warnings are returned as bounded `OboDiagnostic` records instead of being printed globally; message wording is concise native text, with source line and warning category retained.

## Explicit checked corrections

- Failed loads leave the old registry and header fields unchanged. Source can change them before an I/O/parse failure.
- Typed XML output uses the actual `MetaValue` unit identity. Source instead chooses the first allowed unit and can dereference an empty set. No implicit unit conversion or replacement is performed; the allowed-unit set remains vocabulary constraint metadata.
- All consumed XML attributes are escaped and XML 1.0 legality is checked. Accession/cvRef/unit fields are not exempt. Typed Empty omits `value`; a present empty string emits `value=""`, as the source overloads distinguish.
- The legacy `xref_analog:binary-data-type:` prefix uses its correct 29-byte length. Source incorrectly strips the ordinary 22-byte prefix length.
- Diagnostic term/id/name/parent lines all go to the selected writer. Source sends parent lines to process stdout. Diagnostic output is not a lossless OBO serializer.
- FNV-1a has an explicit 64-bit result. Source uses `size_t`, whose width depends on the target.

The source defects and fixes are recorded in the shared [C++ issue ledger](../OpenMS_CPP_ISSUES.md): CPP-021 (typed units), CPP-022 (xref_analog prefix), CPP-023 (selected output stream) and CPP-024 (XML attribute escaping). These entries distinguish direct source evidence from executed C++ results.

Typed scalar and list XML values use source full-precision formatting. The existing parameter float formatter is reused; no conversion through an integer-narrowing ParamValue copy occurs. Empty/list values are supported here even when a particular later mzML adapter has a narrower transport policy. Text rendering validates only the fields it consumes.

## Complete pinned providers

| Source load | Original file | Definitions loaded | Registry size after load | Name-index entries after load |
| --- | --- | ---: | ---: | ---: |
| MS | `psi-ms.obo` | 3,569 | 3,569 | 3,569 |
| PATO | `quality.obo` | 1,976 | 5,545 | 9,073 |
| UO | `unit.obo` | 285 | 5,780 | 10,384 |
| BTO | `brenda.obo` | 3,402 | 9,182 | 14,021 |
| GO | `goslim_goa.obo` | 72 | 9,254 | 16,852 |

These totals come from an independent Python transcription of source provider/index loops, and all native records and aliases are compared against it. The five files contain 9,304 raw Term stanzas; 50 UO accessions overwrite records already included by PSI-MS. The final registry has no unresolved-parent placeholders for these particular files. Final source getters are name `GO`, label `gene_ontology`, version `4.1.155`, and the retained PSI-MS URL.

Original resources total 2,340,432 bytes and remain byte-identical under `resources/cv`. Their Git attributes disable line-ending normalization. Tests check exact raw byte lengths and independently computed FNV values; the generator/provenance additionally verify SHA-256.

The BTO file is historical legacy text, not valid UTF-8: 335 bytes occur in invalid UTF-8 sequences. It is explicitly decoded as Windows-1252. General readers default to UTF-8 and expose `OboEncoding::Windows1252` as a caller-selected alternative; unsupported/undefined bytes are checked errors, never replacement characters. There is no encoding auto-detection.

Eight original NULs in the damaged `BTO:0002243` hypanthium description are preserved, along with its surrounding legacy text. No scientific definition is guessed. The raw file and explicit transcode policy preserve evidence; the description need not be representable as XML unless that field is actually consumed. Ontology data notices are separate from the code license: see [resource notices](../resources/cv/NOTICE.md).

## Limits and shared operations

Defaults are 16 MiB input, 1 MiB per physical line including delimiter, 100,000 registered terms, 1,000,000 parsed/derived entries per load, 1 billion conservative work units, 256 MiB cumulative logical allocation and 16 MiB rendered output. Limits are configurable through `VocabularyLimits`. Ordinary public value field access/Clone/equality/drop have ordinary Rust costs.

The work allowance includes conservative string-comparison bounds for all five cumulative loads and their repeated source index passes; it is not an estimate of CPU instructions or elapsed time. The checked full singleton uses 910,737,238 work units and 200,662,366 logical bytes on the validated 64-bit host, within the defaults. A raw whole-file slice and a small buffered file reader charge only bytes actually examined, preventing a quadratic whole-buffer overcharge for short lines.

All consuming operations precharge descriptors, strings, sparse map roots, geometric vector/string growth, parser warnings, registry drafts and graph frames. A map root reserves at least eleven full element slots even for one large definition. The default singleton uses one cumulative work/byte pair across all five providers. Shared crate-private lookup, clone, ancestry, descendant and first-name helpers accept a caller's remaining counters so XML consumers can retain one document-wide budget.

The public callback traversal charges its graph/ID operations; caller callback work is the caller's responsibility. The first-name helper charges name comparisons itself and avoids borrowing conflicts for consumers that need a metered source-order search. Bounded copies preflight the full owned payload before `Clone`.

## Validation and boundaries

The exact source class fixture is [controlled_vocabulary_source.obo](../tests/data/controlled_vocabulary_source.obo). The direct suite covers all source literals, all fields/type aliases, score direction and XML strings; independent tests cover cumulative aliases/children, placeholder order, source brenda gates, diamonds/cycles, deep traversal, encoding, typed units, XML legality, atomic failures and tiny buffered readers. Private tests cover the complete shared singleton budget, repeated-query exhaustion and sparse large map roots.

All 9,254 final definitions and 16,852 exact name-index entries are checked against [independent term projection](../tests/data/controlled_vocabulary_projection.tsv) and [alias projection](../tests/data/controlled_vocabulary_aliases.tsv). The generator reads only original source resources and does not invoke Rust to manufacture expected output. No executed C++ provider differential result is claimed.

This operation group supplies ontology values, the full source OBO projection and complete pinned providers. It does not itself implement CVMappingFile/CVMappings/SemanticValidator, arbitrary experimental-header transport, ontology edits, XML schema validation, or a general OBO roundtrip writer. Those capabilities remain separate consumer work.
