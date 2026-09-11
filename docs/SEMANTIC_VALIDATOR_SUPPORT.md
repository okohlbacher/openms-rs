# General semantic CV validation

The optional `semantic-validation` feature provides
`format::semantic_validator::{SemanticValidator, ParsedCVTerm, ValidationOptions,
ValidationLimits, ValidationReport}`. It depends on existing `cv-mapping`, with
no new library dependency. The feature is not included in default features.

This is the complete class-specific general `SemanticValidator.h` operation group
at SDK `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. It is a mapping/CV validator,
not an XSD validator, XPath evaluator, referenceable-parameter-group expander,
ontology downloader, or format-specific validator. It does not replace mzML's
separate pinned header predicate. Inherited general XMLFile APIs and C++ virtual
subclass extension points are outside this group.

## Public operations

`SemanticValidator::new(&mapping, &cv)` borrows immutable published CVMappings and
ControlledVocabulary instances. No complete rules or ontology records are cloned.
Its public `options` provides the source setter operations through direct fields:
`tag`, `accession_attribute`, `name_attribute`, `value_attribute`,
`unit_accession_attribute`, `unit_name_attribute`, `check_term_value_types` and
`check_units`. Defaults are `cvParam`, `accession`, `name`, `value`,
`unitAccession`, `unitName`, true and false, respectively. The separate `limits`
field controls native resource consumption. Options may be changed between calls.

`validate(path)` and `validate_reader(impl BufRead)` return an owned
ValidationReport. `report.errors` and `report.warnings` retain source order within
their respective vectors; `report.is_valid()` tests only that errors are empty.
Semantic failures return `Ok(report)`. I/O, malformed XML, missing required CV
parameter attributes, required graph lookup failures and resource exhaustion
return `Err`; a partial report is never published. Reuse after an error is safe.

`ParsedCVTerm` exposes all eight source fields, including the three presence
flags independently of their strings. Default, Clone and full equality are
ordinary owned Rust value operations. `locate_term(path, &parsed_term)` only
checks mapping selection against its accession; name, value, units and presence
flags do not participate. An absent path produces a stable checked error even
after validating a document containing that path. A present path whose rules
have no matching term returns false. An exact use_term match may succeed without
an ontology record; a requested descendant walk needs a valid ontology root.

## Source evaluation behavior

Rules are indexed by literal element path, retaining vector order. XML qualified
names are compared as written, without namespace-URI substitution. The source
empty ancestor path is `/`; a CV parameter at the document root is therefore
looked up at `//cvParam/@accession` with default names. Closing that root checks
`/cvParam/cvParam/@accession`, exactly as the source callback does.

Only matching parameter tags require accession and name attributes. Unknown
accessions emit one warning and skip all remaining semantic checks. Obsolete
known terms emit a warning and continue. Each rule stops at its first matching
mapped term: exact accession when use_term is true, then strict descendant
matching when allow_children is true. The provider's ordered child traversal is
used, preserving its documented stale-child behavior after additive OBO loads.
One incoming parameter can fulfill multiple rules. Counters share keys by path,
rule ID and mapped accession, preserving repeated rule-ID/accession interactions.

On element close, repeat violations are reported first in rule/term order,
followed by requirement violations in rule order. MUST AND requires all mapped
entries; MUST OR requires at least one; MUST XOR requires exactly one. MAY AND
allows none or all; MAY XOR allows at most one. The source's SHOULD AND/XOR
branches do not enforce either condition and are preserved. MAY/SHOULD OR always
passes. Empty rules, duplicate mapped accessions, and rules for elements that
never occur follow these same source branches; missing elements are not inferred.
The scope_path, use_term_name, term_name and cv_identifier_ref fields do not alter
this general class's evaluation.

The diagnostic sequence for a known term is obsolete warning, unit diagnostics,
location diagnostic, canonical-name error, then value error. Canonical names are
trimmed using only space/tab/CR/LF and compared case-sensitively. Optional units
are checked only when enabled. Exact permitted units or descendants are accepted;
missing/unknown/unrelated units error. Supplying a unit on a unitless term warns.
Unit-name presence remains distinct from accession presence; no unit-name match
is invented. The descendant comparison corrects CPP-039.

## Value conventions are source-compatible, not full XSD conformance

All eleven finite CV XRefType variants are handled. Absent or empty values error
for a required type, except that an explicitly present empty String is allowed.
Nonempty values on None are forbidden except for the source `PATO:` exception.
Integer and sign-restricted variants use source i32 conversion, including its
leading-plus/whitespace behavior and bounds. Decimal uses source f64 conversion,
including supported NaN/Inf and checked overflow/underflow. Boolean accepts
ASCII-trimmed, case-insensitive true/false plus 1/0. anyURI only requires a colon.

The source xsd:date branch calls DateTime::set. Thus ordinary date-only strings
such as `2001-02-03` are rejected while supported complete date-time strings are
accepted (CPP-046). This compatibility behavior is deliberate; no strict XSD
calendar parser was added. Existing checked DateTime and ListParse native
boundaries also apply. Disabling value checks skips all these value diagnostics.

## XML, transport and shared bounds

The existing CV mapping reader's bounded document decoder, XML lexical checks,
attribute handling and element stack now reside in one private `cv_xml` module.
A synchronous Start/End callback serves the mapping loader and this validator;
there is no owned XML tree or public parsing framework. Mapping record callbacks,
namespace option, limits, source behavior and atomic output publication remain
unchanged, and their existing fixture/projection tests remain required.

Input supports UTF-8 and UTF-16 LE/BE with checked BOM/declaration agreement.
US-ASCII and Latin-1-compatible declarations are accepted only for ASCII payload;
non-ASCII Latin-1 XML is explicitly unsupported. XML 1.0 character/name and
attribute-separator rules, duplicate attributes, entities, declarations, comments,
PI, CDATA, tag balance and single-root/text constraints are checked. Literal
attribute whitespace is normalized before reference expansion. DTD and external
entities are unsupported; schema hints cause no network access. Scientific text
is ignored only after lexical validation. This is the same bounded qname parser
policy as the CV mapping loader, without namespace-URI-based rule rewriting.

Paths reuse magic-based gzip/bzip2 transport, independent of filename suffix.
The complete decoded stream is consumed, including checksum/truncation failures.
ZIP is unsupported for this single-document operation. A stream reader receives
already decompressed XML. Zero/one-byte inputs are rejected safely.

Default caps are 16 MiB decoded input/normalized UTF-8, depth128, 1M XML elements,
1M encountered CV parameters, 100k mapping rules, 1M mapped terms, 100k total
diagnostics, 50M conservative work units and 128 MiB logical allocations/copies.
They are configurable resource limits, not wall-clock or physical-memory claims.

One operation counter pair covers mapping indexing, bounded input/decoding,
lexical scans, attributes, path construction, sparse map roots, ordered counters,
CV lookup/traversal comparisons, name/value validation and diagnostic formatting.
Graph calls receive the same counters rather than fresh provider budgets. Rule
and CV payload stays borrowed. Ignored XML and unknown terms consume applicable
input limits. Count/report growth and formatting are charged before allocation;
source u32 repetition overflow errors instead of wrapping. Ordinary Clone on
public value objects has ordinary Rust cost. CV construction is separately bounded
by the caller's vocabulary operations.

## Evidence and intentional source corrections

The unchanged source valid XML has 266 elements and 111 CV parameters, producing
zero errors/warnings. The unchanged corrupt XML has 271 elements and 119 CV
parameters, producing the class test's exact five errors and four warnings in
order. Tests use the original 738-term historical SemanticValidator CV plus the
published quality/unit/brenda/goslim_goa resources in source load order. The
modern PSI singleton is not substituted. Brenda's opaque historical bytes use
the existing Windows1252 decoding policy; all raw bytes stay unchanged.

Fifteen direct tests cover those literals and independent rule truth tables,
collisions, first-match traversal, stale children, custom names, root paths,
unit/value branches, error reuse, XML/transport and configured limits. Four private
regressions check diagnostic preflight, counter overflow, propagation of graph
resource errors, and cumulative path/count/report exhaustion after successful
terms. A separate Python projection verifies raw SHA-256 values, parses fixture
counts independently, and records all nine class-test diagnostic literals with
source line numbers. It never imports or executes Rust. No C++ executable oracle
is claimed. The source inventory includes 25 hashed source/dependency files.

The central [C++ issue ledger](../OpenMS_CPP_ISSUES.md) describes CPP-039 (unit
fallback target), CPP-040 (failed-parse state contamination), CPP-044 (state-dependent
absent-path locateTerm cache) and CPP-046 (xsd:date uses full DateTime parsing).
The first three are explicit native corrections; the date convention is retained.
The existing conservative XML, encoding, graph-cycle and numeric/resource
boundaries remain documented native adaptations.

Source implementation is BSD-3-Clause. Historical ontology notices and original
files are preserved; no new ontology license grant is inferred. Existing CV
resource licensing cautions still apply. Tests/data byte-preservation attributes
already cover these originals, so no new broad attribute rule was added.
