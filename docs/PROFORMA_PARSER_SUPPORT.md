# ProForma text parsing and structured errors

`chemistry::proforma` implements both complete text grammars and the structured
parse-error interface from OpenMS4-core `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`.
Parsing creates the [owned annotation AST](PROFORMA_SUPPORT.md). This adds text
parsing to the existing writers; the complete ProForma header still has separate
chemistry resolution, sequence conversion, mass and spectrum operations
that are not implemented here. [JSON transport](PROFORMA_JSON_SUPPORT.md) is
available separately with the `proforma-json` feature.

## Using the two grammars

```rust
use openms::chemistry::proforma::{Peptidoform, PeptidoformIon, WriteMode};

let peptide = Peptidoform::parse("M[UNIMOD:35]PEPTIDE")?;
assert_eq!(peptide.to_text(WriteMode::Lossless)?, "M[UNIMOD:35]PEPTIDE");
let ion = PeptidoformIon::parse("AC/2+DE/3")?;
assert!(ion.is_chimeric);
assert_eq!(ion.chains.len(), 2);
Ok::<(), Box<dyn std::error::Error>>(())
```

Both methods return `std::result::Result<Self, ParseFailure>` and publish only a
complete value. `Peptidoform::parse` accepts a single chain without charge or
chain separators; `PeptidoformIon::parse` additionally accepts `//`, `+` and charge
suffixes. It is not necessary to create a parser object or supply a registry.
The implementation uses a byte cursor and bounded lookahead rather than building
a second token collection.

The parser retains all source constructs: repeated chain-name prefixes, global
modifications and isotope replacements, unlocalised and labile modifications,
N/C-terminal modifications, individual residues, ambiguous regions and modified
ranges, alternative modification tags, labels and scores, formulas and glycans,
position constraints, CV accessions, named modifications and mass deltas, simple
charges and adduct lists. All source hints, optional fields and encounter order
are retained. Resolved modification handles stay `None` and cross-link groups
stay empty; neither is synthesized by the source parser.

Parsing performs syntactic work only. Unknown modification or glycan names,
accessions and unvalidated formula text remain available in the AST. It does not
invoke [modification registries](MODIFICATION_SUPPORT.md),
[MonosaccharideDB](MONOSACCHARIDE_SUPPORT.md), or the AASequence parser.

## Source grammar details

Whitespace is not skipped. ASCII letter runs form identifier tokens, including
lowercase letters and uncommon residue codes; sequence parsing accepts each ASCII
letter without a biochemical lookup. Other bytes become single identifier tokens.
Complete Unicode annotation names, names assembled from tokens, INFO text and
formula text are retained. Individual residue/position symbols use the source
ASCII grammar.

Numeric tokens accept decimal fractions with optional signs, including `.5`, but
not scientific notation. A dot joins a number only when a digit follows it.
Integer fields preserve the source conversion's integral prefix: a count of `2.5`
becomes `2`. Numeric CV accessions retain the original token text, including
leading zeros and fractional or signed spelling. Simple
and adduct charges preserve the tokenizer's separate-sign behavior, including
`--1` becoming `1`. Negative and zero counts are retained. Mass spellings remain
in `MassDelta.original_text` for lossless writing.

Repeated names join with `" / "`; `(>name)` and `(>>name)` both name the chain.
An empty name is omitted, and the ion-level name is never populated. Mixed `//`
and `+` separators are accepted; any `+` sets the single chimeric flag. A charge
before `+` belongs to the preceding chain. A final charge belongs to that chain
only when the chain followed `+`; otherwise it belongs to the ion. This contextual
source rule is retained instead of assigning every suffix to the last chain.

Some finite source behavior is permissive or asymmetric:

- An empty `(?)` ambiguous region counts as a sequence section; an empty `()`
  modified range is rejected. A completely absent section is an empty-sequence
  error. Parentheses inside ordinary named
  modifications are tracked without requiring a final zero nesting balance.
- Unlocalised bracket groups become separate entries. Global location lists may
  be empty or have a trailing comma. Labels are forbidden on global and labile
  modifications but are not required to have a matching cross-link partner.
- INFO stops at square-bracket, pipe, hash or comma delimiters; it does not treat
  a labile closing brace as an INFO terminator.
- The adduct writer emits a final aggregate-charge suffix that the source parser
  does not consume. Writing then parsing is not a general round-trip guarantee.

These choices follow the actual source parser, rather than adding stricter
ProForma-standard validation or inferring semantics from comments.

## Structured errors

`ParseFailure::Syntax(Box<ParseError>)` carries source-style diagnostics;
`ParseFailure::Resource(openms::Error)` carries checked input, work, nesting,
capacity or allocation failures. Resource errors are never disguised as an
additional source error-code variant.

`ErrorCode` exposes all 15 source variants and `as_str()`. `ParseError` provides
`new(code, position, input, message)`, `code()`, `position()`, `message()`,
`context_before()`, `context_after()`, `expected()`, `found()`, atomic
`set_expected_found(expected, found)`, and `formatted_message()`. Construction and
string-producing/mutating operations return the crate's checked `Result`.
`Display` and `std::error::Error` allow ordinary Rust error reporting.

Positions and context windows count **bytes**. The constructor clamps a supplied
position to the input length, keeps at most 20 bytes on each side and never copies
the whole input. Context getters return byte slices because a source window can
split a UTF-8 character. Human-readable formatting replaces only invalid display
fragments with U+FFFD and retains the source marker and ellipsis rules.

The source parser leaves expected/found strings empty even when its custom
message names an expected token. For example, `A[+1` yields
`UnexpectedCharacter` at position 4 with message `Expected ']'`, rather than the
header's illustrative `UnclosedBracket` diagnostic. Several declared source error
codes are never emitted by parsing; all remain constructible. The source C++
exception's file/line/function and global exception-handler effects are not a
Rust ABI interface.

## Checked boundaries and shared limits

The following public constants apply separately to each parse call:

| Constant | Limit |
| --- | ---: |
| `MAX_PROFORMA_PARSE_INPUT_BYTES` | 4 MiB |
| `MAX_PROFORMA_PARSE_NODES` | 1,000,000 appended AST collection items |
| `MAX_PROFORMA_PARSE_WORK` | 50,000,000 charged units |
| `MAX_PROFORMA_PARSE_BYTES` | 256 MiB cumulative requested buffer capacity |
| `MAX_PROFORMA_PARSE_DEPTH` | 256 bracket-lookahead or nested-name levels |
| `MAX_PROFORMA_DIAGNOSTIC_BYTES` | 64 KiB combined message/expected/found text |

One ledger covers all chains, copied token text, buffer growth, UTF-8 validation,
lookahead rescans and construction of syntax diagnostics. Limits are precharged;
checked capacity growth uses fallible standard allocations. Byte accounting is a
conservative sum of requested capacities, including old capacities when a vector
grows; it is not an exact physical-memory or allocator-bookkeeping ceiling.
There is no independent budget reset for another chain or lookahead. The node
limit counts collection entries, including alternatives and chain entries, rather
than claiming to count every scalar or aggregate field.

Integer conversions and separate-sign multiplication must fit `i32`. Numeric
mass and score values must be finite. Representable binary64 subnormals are kept,
while a nonzero decimal that rounds to zero is rejected; literal signed zero is
retained. This is a portable checked policy: C++ `std::stod` may report ERANGE even
for representable subnormals on some platforms. Decimal parsing is locale
independent, while source C character classification/conversion may follow its
process locale.

An individually published String must be valid UTF-8. For example, source parsing
of `A[Glycan:é]` would create two invalid one-byte glycan names; the native parser
returns a syntax error at that field instead. Intact Unicode assembled into any
other annotation field remains valid. This is a representation boundary, not a
blanket rejection of non-ASCII annotations. Other source overflow, undefined
conversion or excessive growth becomes a checked error without partial output.

## Evidence and remaining scope

The implementation follows the pinned
[source tokenizer and parser](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L40),
[structured errors](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/CHEMISTRY/ProForma.h#L556)
and [class tests](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/tests/class_tests/openms/source/ProFormaParser_test.cpp).
[Native parser tests](../tests/proforma_parser.rs) execute all 176 positive and 22
negative source fixture lines and all 62 component-grammar calls. Positive lines
use the source test's ion-dispatch heuristic. Negative lines use its single-chain
parser call; their success does not imply cross-link validation by the ion parser.
Source literal AST assertions and independent Unicode, numeric, diagnostic,
permissive-grammar and budget cases supplement these fixtures. Four private tests
exercise shared node/work/capacity limits and late error publication.

[The parser manifest](../tests/data/proforma_parser_provenance.json) records exact
source and copied-fixture hashes, case counts and native differences. These are
native source-derived tests; this parser group does not claim an executed C++
parser differential run or a complete SDK build. The earlier writer extraction
probe remains separate evidence in [the AST/writer support document](PROFORMA_SUPPORT.md).
[Modification resolution](PROFORMA_RESOLUTION_SUPPORT.md) is available as a
separate explicit operation. Conversion policies, mass/m/z and ordinary/XLMS
spectra remain unimplemented; no parser success claims those capabilities.

## Executed C++ text comparisons

A separate [extraction probe](../tests/data/proforma_parser_probe_provenance.json)
compiled the exact source AST, tokenizer, complete parser and writer, plus the
exact string-prefix helper. All 238 inputs are run through both grammars: the
198 unchanged source grammar lines plus 40 targeted branch inputs. The resulting
476 comparisons match native Rust: 382 accepted cases match both text writer
outputs, and 94 rejected cases match error code, byte position and original
message. [The regression](../tests/proforma_parser_probe.rs) reads only frozen
C++ outputs and does not require a C++ compiler at build or test time.

This is a source extraction, not a full SDK build. A small exception adapter
captures the parser's unmodified throw arguments; it does not run source exception
formatting or global-handler state. JSON, registry resolution and scientific
backends are absent. The native diagnostic tests above independently check
source formatting; the documented Unicode, numeric and resource boundaries
remain applicable. Probe generation and compiler identity are recorded with
all extracted block and fixture hashes. Expected data was not produced by Rust.
