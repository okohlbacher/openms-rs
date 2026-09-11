# Monosaccharide database

`chemistry::MonosaccharideDB` implements the complete public lookup surface of
Core SDK `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`'s MonosaccharideDB. It is a
small immutable built-in vocabulary for ProForma glycan symbols, containing all
**24 primary records and 12 synonyms** from the pinned source JSON. The source
header, implementation, class test and JSON are byte-identical to revision 54a.
See [source provenance](../tests/data/monosaccharide_provenance.json) and
[data notices](../resources/monosaccharides/README.md).

```rust
use openms::chemistry::MonosaccharideDB;

let db = MonosaccharideDB::global();
let sialic_acid = db.get_or_error("NeuAc").unwrap();
assert_eq!(sialic_acid.symbol, "Neu5Ac");
assert_eq!(sialic_acid.mass, 291.095416506);
assert_eq!(sialic_acid.formula, "C11H17N1O8");
```

The record's public owned fields are `symbol`, `name`, `mass`, `formula`, and
`synonyms`. Registry lookups borrow immutable records. Cloning a record gives the
caller independent editable strings and synonyms; it does not add aliases or
change the registry. The stored mass is a binary64 literal independent of the
formula text. No formula parsing, elemental recomputation, rounding to four
places, neutral-loss interpretation or glycan sequence parser is introduced.

## Public source mapping

| Source | Native |
| --- | --- |
| `getInstance()` | Thread-safe `MonosaccharideDB::global() -> &'static Self` |
| `hasSymbol(symbol)` | `has_symbol(&str) -> bool` |
| `getMonosaccharide(symbol)` | `get(&str) -> Option<&Monosaccharide>` |
| `getMonosaccharideOrThrow(symbol)` | `get_or_error(&str) -> Result<&Monosaccharide>` |
| `getAllSymbols()` | `all_symbols() -> Vec<&str>` |
| `getNumberOfMonosaccharides()` | `len()` and convenience `is_empty()` |

Symbols are exact and case-sensitive. The registry does not trim whitespace,
normalize Unicode or interpret a full descriptive name as a synonym unless that
name is explicitly listed. `NeuAc` and `Sialic Acid` resolve to the same borrowed
`Neu5Ac` record. `all_symbols` returns lexically sorted **primary symbols only**;
it borrows their strings instead of making source-style string copies. Unknown
queries return `None` or a native `InvalidValue` error, replacing null and
`ElementNotFound`. Error text deliberately does not copy a caller's possibly
large unknown symbol.

All lookup goes through the source synonym map, including primary names. Loading
follows lexical primary-key order, because the source JSON object uses an ordered
map. For each accepted primary, synonyms are inserted in their stored list order,
then that primary maps to itself. Repeated keys overwrite earlier alias targets.
Consequently, a later primary's synonym can shadow the spelling of an earlier
primary; lookup does not give primaries unconditional priority. This collision
case does not occur in the packaged data but is pinned by a private synthetic
builder test. Synonym lists themselves keep repetitions and input order.

## Built-in data and runtime boundary

The original [monosaccharides.json](../resources/monosaccharides/monosaccharides.json)
is packaged unchanged. A [standard-library Python generator](../tools/generate_monosaccharides.py)
produces the [embedded rows](../resources/monosaccharides/monosaccharides.rs),
storing exact binary64 bits and the original record strings. This allows identical
built-in values with all optional features disabled and requires no JSON parser
or new dependency at runtime.

The source JSON loader is private and discovers an installed file through
`File::find`. The native database uses this fixed compiled dataset: it does not
read an environment path, discover an installation, or allow a runtime file
replacement. There is no new public JSON provider or caller-registry constructor.
Updating the packaged vocabulary requires regenerating and rebuilding. This is
an explicit native data-discovery difference, not an assertion that source
installation overrides are supported.

Generation retains the source missing-mass skip and optional-field defaults:
missing names use the symbol, missing formulas are empty, and non-array synonyms
are ignored. Present text fields must have the expected type. The generator is
for the packaged scientific data, not a general replacement for every
nlohmann::json conversion corner; it rejects nonfinite masses and control
characters rather than emitting invalid or unbounded source code. Its input cap
is 1 MiB and 1,024 primary rows. The present data uses 24 rows and ordinary finite
positive masses. Runtime construction is bounded entirely by those fixed rows;
no caller query controls an allocation, traversal or copy of record payload.
`all_symbols` allocates only its 24 borrowed references. There is no configurable
execution-budget framework for this fixed vocabulary.

## Verification

The [native tests](../tests/monosaccharide_db.rs) cover every source class-test
accessor section, singleton identity, all 24 literal mass/formula pairs, exact
sorted symbols, all 12 aliases, owned clone independence, long/Unicode/whitespace
queries, and concurrent lookup. With `rna-json` enabled, an independent
serde_json decode of the untouched original compares every name, formula,
synonym and binary64 mass bit against the native registry. That decoder is used
only by the test; the module itself has no feature requirement.

A private synthetic test checks lexical load order, repeated synonyms,
primary/alias collision precedence, empty construction, and independent signed
stored mass/formula values. This is derived source-behavior evidence; it is not
an exposed custom registry interface. Existing historical fixture manifests are
unchanged. No C++ build or executable comparison was run.

Run `python3 tools/generate_monosaccharides.py --check` to verify deterministic
projection without a compiler or network. The integration handoff recommends
this check for the repository's generated-resource CI audit.
