# ProForma JSON transport

With the **`proforma-json`** feature, `chemistry::proforma` implements all four
public JSON operations from OpenMS4-core
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`: reading and writing `Peptidoform` and
`PeptidoformIon`. It reuses the existing pinned `serde_json` dependency without
coupling JSON availability to the unrelated `rna-json` feature. No default feature
is changed by this increment.

```toml
[dependencies]
openms = { path = "../openms-rs", features = ["proforma-json"] }
```

```rust
use openms::chemistry::proforma::Peptidoform;

let peptide = Peptidoform::parse("EM[UNIMOD:35]K")?;
let json = peptide.to_json()?;
let restored = Peptidoform::from_json(&json)?;
assert_eq!(restored, peptide);
Ok::<(), Box<dyn std::error::Error>>(())
```

Both types expose `to_json(&self) -> openms::Result<String>` and
`from_json(&str) -> openms::Result<Self>`. Calls return complete results or errors;
writing borrows the input and cannot partly modify it. JSON errors use the crate's
ordinary error type. Source JSON wrappers likewise use generic parse exceptions,
rather than the ProForma-specific [structured text errors](PROFORMA_PARSER_SUPPORT.md).
The exact C++ exception messages and nlohmann error numbers are not an API promise.

## Complete source schema

The implementation maps the [source schema helpers](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProFormaDataJson.h)
explicitly. It does not derive a Rust enum schema through serde. Objects use the
source field names and lexical key order; arrays retain encounter order.

| Structure | Source JSON form |
| --- | --- |
| Modification tag | `{"type":...,"value":...}` with `cv_accession`, `named_mod`, `mass_delta`, `formula`, `glycan`, `info`, or `position` |
| CV accession | `database`, `accession`; database is `UNIMOD`, `MOD`, `RESID`, `XLMOD` or `GNO` |
| Named modification | `name`, optional `cv_hint` |
| Mass delta | `source`, `mass`, `original_text`; source is `NONE`, `OBS`, `U`, `M`, `R`, `X` or `G` |
| Formula | `formula_string`, optional `charge` |
| Glycan | Array of `monosaccharide`/`count` objects; each monosaccharide is tagged `name` or `formula` |
| INFO / position | `text` / `residues`, `n_term`, `c_term` |
| Label | `type`, `identifier`, optional `score`; type is `CROSSLINK`, `BRANCH` or `AMBIGUOUS` |
| Modification | Array of alternative objects containing `tag` and optional `label` |
| Sequence section | Tagged `element`, `ambiguous_region` or `modified_range` |
| Sequence element | One-byte `amino_acid` string and `modifications` array |
| Ambiguous region / modified range | `elements` / `elements` and `modifications` |
| Global entry | Tagged `isotope_replacement` (`isotope`) or `global_modification` (`modification`, `locations`) |
| Unlocalised / labile modification | `modifications` plus optional `occurrence` / `modification` |
| Charge | Tagged `simple` integer or `adducts` array of `formula`, `charge`, optional `occurrence` |
| Chain | `global_mods`, `unlocalised_mods`, `labile_mods`, `n_term_mods`, `sequence`, `c_term_mods`, optional `name` and `charge` |
| Ion | `chains`, `is_chimeric`, optional `name` and `charge` |

JSON preserves chain charge and ion name, even though the source text writers omit
these in some contexts. It transports syntactically unusual stored AST values,
including empty sequences/regions/adduct lists and negative counts. Reading JSON
does not call the text parser, resolve formulas or glycans, check biochemical
residue codes, or validate cross-link partners.

`Modification.resolved_mod` is deliberately **not transported**. The writer does
not inspect, clone or serialize its external immutable chemistry. Every read
reconstructs it as `None`, matching the source pointer reset. Caller annotations
and label data survive; a resolved handle must be obtained separately when a
resolution API becomes available. `CrossLinkGroup` has private source ADL helpers
but is not a field reachable from either public JSON wrapper, so no extra native
standalone codec is invented for it.

## Input defaults and unusual source behavior

A chain may omit `global_mods`, which then becomes empty. Its other five vector
fields are required. An ion must have `chains`; omitted `is_chimeric` becomes
false. Omitted or explicit null optional names, charges, hints, scores and
occurrences become `None`. The writer omits these absent options. Position flags
may be omitted and default to false. Explicit null for either position flag or
`is_chimeric` is an error.

Ordinary vector fields require JSON arrays, even for empty values; explicit null
is not interchangeable with `[]`. Three source paths use generic JSON iteration:
**modification alternatives, glycan components and global modification locations**.
They accept arrays, object values in lexical key order, and null as an empty
collection. A primitive contributes one value, so a string used as `locations`
becomes a one-element location list. Primitive alternatives/components still fail
when the required item fields are absent. Writing these forms produces arrays.
Unknown fields are ignored, duplicate object keys retain the last value and all
variant/enum names remain case sensitive.

Source integer conversions accept JSON booleans as 1/0 and truncate floating
values toward zero. The native codec preserves this for charges, counts and
occurrences but requires the converted result to fit `i32`. Floating mass and
score fields reject booleans; boolean fields require actual booleans. Integer
conversion overflow becomes a checked error instead of implementation-dependent
narrowing or undefined float-to-integer conversion.

The source accepts a leading UTF-8 BOM. Bare JSON integer `-0` becomes integer
zero, while `-0.0` and `-0e0` retain negative floating zero. A bounded lexical
normalization preserves that distinction before serde decoding; it only matches
complete value tokens and never changes strings or exponent signs such as
`1e-0`. JSON exponent syntax and underflow to floating zero follow the JSON
library's numeric path, separately from the stricter decimal-only text grammar.

## Checked representation and resource limits

Ordinary String fields preserve valid Unicode, escaped control characters and
embedded NUL. `amino_acid` must contain exactly one ASCII byte, which may be
punctuation or NUL; the codec does not require an amino-acid letter. Position
residue characters use the same ASCII byte boundary. Non-ASCII source byte-vector
states are an explicit native representation restriction; annotation strings are
not restricted to ASCII. Invalid JSON Unicode escapes are rejected.

Stored floating values must be finite. Source nlohmann writing replaces NaN or
infinity with null, which cannot be read back as a required mass; native writing
returns an error instead. Finite binary64 values, representable subnormals and
floating signed zero are retained. Compact JSON numeric spellings may differ
between nlohmann and serde_json, particularly exponent formatting; semantic JSON
values and float round trips are the compatibility target, not universal
byte-identical floating text.

Public limits, shared across one complete read/write operation:

| Constant | Limit |
| --- | ---: |
| `MAX_PROFORMA_JSON_TEXT_BYTES` | 4 MiB input or output |
| `MAX_PROFORMA_JSON_ITEMS` | 1,000,000 charged structural/collection items |
| `MAX_PROFORMA_JSON_WORK` | 50,000,000 charged work units |
| `MAX_PROFORMA_JSON_BYTES` | 256 MiB cumulative logical capacity allowance |
| `MAX_PROFORMA_JSON_DEPTH` | 64 JSON container levels |

A bounded lexical preflight covers the entire input, including duplicates and
unknown fields, before serde allocates its document tree. It precharges
conservative sorted-map comparisons, structural storage and string scratch;
serde performs the actual JSON syntax/number/Unicode parsing. AST construction,
normalization, document building and serialization share the remaining ledger.
Each inserted writer field includes a full sparse BTree node allowance derived
from key/value size plus headroom, even for a one-field object. Output uses a
bounded sink, with fallible growth before writing each chunk.

Work estimates are deliberately conservative and may reject input before its
actual traversal would exhaust the allowance. Capacity accounting sums requested
storage, including growth and temporary representations; it is not an exact
physical-memory or allocator-overhead ceiling. Standard serde/map allocations
remain ordinary Rust allocations inside the preflight bounds. There are no
per-chain budget resets and no partly published AST or JSON text on failure.

## Evidence and remaining work

[Tests](../tests/proforma_json.rs) retain both source class-test JSON examples and
round-trip all 176 positive source grammar ASTs. Independent tests specify the
complete tagged schema, defaults, generic iteration, enum cases, coercions,
Unicode, duplicate keys, signed zero/exponents, resolved-handle omission and late
failures. Four private tests cover shared ledgers, sparse maps and depth/escape
preflight. The [manifest](../tests/data/proforma_json_provenance.json) distinguishes
source literals from derived assertions and records exact source/fixture hashes.
This group claims no executed C++ JSON differential run or full SDK build.

The [overall ProForma support](PROFORMA_SUPPORT.md) covers the pinned public operation groups.
[Chemistry resolution](PROFORMA_RESOLUTION_SUPPORT.md) is a separate explicit
operation; decoded handles initially remain unset. [Mass/mz operations](PROFORMA_MASS_SUPPORT.md)
and [AASequence conversion](PROFORMA_CONVERSION_SUPPORT.md) can resolve decoded
annotations using a supplied registry. [Ordinary/XLMS spectrum generation](PROFORMA_SPECTRA_SUPPORT.md)
is also available through a separate explicit operation.
