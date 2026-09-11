# CV mapping records and XML loading

`data_structures::{CVReference, CVMappingTerm, CVMappingRule, CVMappings,
RequirementLevel, CombinationsLogic}` provides the complete owned public value
operation group from the five pinned SDK headers. The optional `cv-mapping`
feature provides `format::cv_mapping::{CVMappingFile, ReadOptions}` and is included
in default features. The value records remain available without any features.

The source revision is `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`.
This is a mapping-file loader and record model. It does not perform semantic or
XSD validation, evaluate XPath, fetch vocabularies, or replace mzML's separately
reviewed header predicate. The source class has no mapping writer, and none is
claimed here. Its inherited general XMLFile operations are a separate API group.

## Public operations and ownership

The three leaf records expose every source field directly. `Default`, `Clone`,
assignment, field replacement and full equality replace C++ constructor,
copy/getter/setter and equality boilerplate. Rule terms and mapping rules are
ordered vectors; assigning a vector replaces it, and pushing appends one record.
RequirementLevel has Must=0, Should=1, May=2. CombinationsLogic has Or=0, And=1,
Xor=2. Rule defaults are Must/Or; all term booleans default false.

CVMappings keeps its ordered reference vector private so identifier queries
remain consistent. `cv_references()` borrows it, and `has_cv_reference()` tests
identifier presence. `add_cv_reference()` ignores an existing identifier and
returns false; the bool replaces process-global source stderr output. A unique
identifier is appended and returns true. Empty identifiers are ordinary keys.

Despite its name, source `setCVReferences` appends all input records, including
duplicates, and cannot clear the container with an empty argument. Native
`set_cv_references(Vec<CVReference>)` preserves this behavior. The separately
named `replace_cv_references` replaces both vector and index. Owned input avoids
source self-alias iterator invalidation: obtain a copy with `cv_references().to_vec()`
if the same records need appending. The source private map's last values are
completely determined by the reference vector, and only membership is publicly
queried. The native presence index therefore preserves all reachable source
public equality/query behavior without duplicate copies of the reference names.

Ordinary value mutation and Clone have normal Rust allocation costs. Resource
limits apply when an I/O operation consumes caller records or external bytes.

## Loading and exact source behavior

`CVMappingFile::default().load(path)` returns an owned CVMappings. `read(BufRead)`
accepts an already decoded stream. `load_into` and `read_into` atomically append
parsed references and replace rules in an existing destination, as a successful
source load does. Existing references absent from the new file remain; repeated
loads can append duplicates. On any I/O, XML, value or resource error, the
existing destination stays unchanged and subsequent calls start with clean
local parser state.

The parser follows source exact qnames: CvReference, CvMappingRule and CvTerm.
Other element names, container attributes and ordinary character content are
ignored. In particular `modelName`, `modelURI` and `modelVersion` are not part of
the five source records. XMLFile disables namespace processing, so a prefixed
mapping qname is also unrecognized; no namespace URI substitution is performed.
This is distinct from the path-string namespace option. Unrecognized content is
still checked for XML well-formedness and counted against input/resource limits.

Required attributes must be present but can be empty. Booleans require exact
lowercase `true` or `false`, without spaces or numeric aliases. A missing or empty
useTermName is false; a missing or empty isRepeatable is **true**, unlike a default
CVMappingTerm value. Unknown present requirementLevel values use Must and unknown
cvTermsCombinationLogic values use Or, preserving the source's explicit fallbacks.
Missing required values still error. No CV identifier existence check is added.

Within one document, the source has one current rule rather than a semantic
nesting stack. Rule starts replace its scalar fields without clearing terms;
CvTerm events append; rule ends publish and reset it. The native loader retains
this behavior even under unusual but well-formed wrapper/nested structures. It
does not pretend these accepted structures have been schema validated.

With strip_namespaces=false, paths remain verbatim. When true, empty slash
components are dropped, a leading slash is produced, and one namespace prefix is
removed from each cvElementPath segment. Unprefixed segments remain, and a
namespaced attribute keeps its `@` axis. Multiple colons in one segment error.
As in source, scopePath is unchanged. Empty paths stay empty. This corrects the
source option's ordinary-path rejection and attribute-axis loss.

## XML and file boundaries

Input supports UTF-8 with optional BOM, and UTF-16 LE/BE with BOM or the XML
signature. The declaration must agree with the bytes; only XML 1.0 is supported.
US-ASCII and ISO-8859-1-compatible declarations are accepted only for ASCII
payloads; non-ASCII Latin-1 bytes are an explicit unsupported encoding rather
than silently replaced characters. XML line endings and literal attribute
whitespace are normalized before numeric/predefined reference expansion, so a
character reference to LF remains LF. All consumed characters and references
must be legal XML 1.0. Duplicate attributes, adjacent unseparated attributes,
invalid names, malformed comments/PI/declarations, unmatched tags, multiple roots,
nonwhitespace outside the root and ordinary-text `]]>` are rejected.

DTD/external entities are explicitly unsupported; there are no network reads.
Schema hints are inert, as in this source loader. The referenced CvMapping XSD
is absent from the pinned SDK resources; no external schema is substituted.
Paths reuse shared magic-based gzip/bzip2 readers with existing optional
compression dependencies. ZIP is unsupported because this API expects one
scientific stream. Compression is independent of filename suffix, the entire
stream is read and corrupt/truncated compression is reported. Zero/one-byte
inputs are handled by length-aware byte slices and checked XML errors.

## Resource limits

ReadOptions keeps one operation's scientific option and explicit resource caps:
16 MiB raw decoded-stream bytes and normalized UTF-8 bytes; depth 128;
1,000,000 XML elements; 100,000 total output references; 100,000 parsed rule ends;
1,000,000 parsed terms; 50,000,000 conservative work units; 128 MiB logical copied
or allocated payload. These are configurable upper bounds, not physical-memory
or wall-clock guarantees.

The reader copies bounded input first and uses borrowed quick-xml events; it does
not build a full owned XML tree. One shared counter pair covers buffer growth,
character decoding, lexical scans, attribute copies, term/rule vectors, path
normalization, destination references and the identifier index. Replaced old
rules are not cloned. Duplicate/ignored records still consume applicable input
limits. The reference cap includes old destination references before any old
payload copy. At most one sentinel byte is inspected beyond the configured
input limit. Index preflight covers a sparse root with at least eleven slots.

## Independent evidence and source fixes

The complete class fixture has one reference, nine rules and 133 terms. The
original class-test first three rule term counts 14/32/46, first six detailed term
records, path/requirement/logic fields and all nine IDs are literal assertions.
Five additional unchanged inputs cover all shipped mapping resources and the
SemanticValidator fixture. An independent Python XML projection checks all 683
reference/rule/term records and every stored field in all six files. It verifies
all original SHA-256 values before projection and never imports or executes Rust
production. Projection rows are source-data expectations, distinguished from
literal assertions in the C++ class test; no C++ executable oracle is claimed.

Related source findings are documented in the central
[C++ issue ledger](../OpenMS_CPP_ISSUES.md): CPP-032 (namespace stripping), CPP-033
(reused failed-loader contamination), CPP-034 (bulk-set self-alias invalidation),
CPP-035 (invalid enum fallback), and CPP-036 (shared XMLFile short-read sniff).
Native namespace and parser-state corrections are deliberate; unknown enum
fallback and ordinary bulk append remain source-compatible. Dedicated regressions
exercise each relevant trigger and destination rollback. The small value/model
group introduces no ontology resource licensing changes; source/fixtures remain
under the pinned repository's BSD-3-Clause attribution, recorded with exact
source hashes in the provenance manifest.

[General SemanticValidator](SEMANTIC_VALIDATOR_SUPPORT.md) is now available through
the optional `semantic-validation` feature. The mapping loader and validator
share a private bounded XML event reader; all mapping-loader semantics and
fixtures remain unchanged. Full XSD and derived format validators remain separate.
