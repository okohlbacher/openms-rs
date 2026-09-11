# Native metadata and Product hashing

`Unit`, `MetaValueData`, `MetaValue`, `CVTerm`, `CVTermList` and `Product`
implement Rust `std::hash::Hash`. This closes the remaining declared hash
operation on the native Product record. The source reference is Core SDK
`54a232fe2cae9c590d5c997fa49d20e7769860fb`.

## Contract and fields

Hashing follows each existing native `PartialEq` implementation. No equality,
validation, unit, conversion, mutation or transport behavior changes in this
increment.

| Type | Fields included |
| --- | --- |
| `Unit` | Accession, name and vocabulary reference |
| `MetaValueData` | Alternative discriminator and complete typed content; list length, order and elements |
| `MetaValue` | Complete data and optional unit identity, including a unit attached to Empty |
| `CVTerm` | Accession, name, vocabulary reference and complete value/unit |
| `CVTermList` | All accession-map keys and ordered term vectors, plus ordinary typed metadata |
| `Product` | Target m/z, lower and upper isolation offsets, and the complete CVTermList |

Scalar doubles and every double-list element normalize -0 to +0 before their
bits are hashed. All other float bits are retained. Structs without direct
floats derive `Hash`, so every equality-significant field is included. The
standard library hashes strings, list boundaries, option variants and ordered
maps. BTreeMap insertion history does not affect the hash input, while duplicate
terms, their vector order and explicitly present empty accession buckets do.

Hashing does not call Display, round values, assign registry IDs, clone
metadata, allocate temporary collections or change input values. Work is linear
in the visited fields/elements and string bytes, with bounded call depth because
the model has no recursive value alternative. As an ordinary standard-library
trait it has no fallible work-budget argument; the caller owns input-size and
Hasher behavior. Any allocation or other work inside the caller's Hasher is its
own responsibility.

```rust
use openms::metadata::{MetaValue, Product};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

let mut product = Product::default();
product.mz = 500.25;
product.cv_terms.metadata.insert("label".into(), MetaValue::from("target"));
let mut state = DefaultHasher::new();
product.hash(&mut state);
let digest = state.finish();
```

## Native boundaries

The chosen Hasher determines the numeric digest. There is no promise of stable
values across releases, processes, architectures or languages. This does not
reproduce source FNV/golden-ratio mixing, DataValue decimal rounding/string
rendering, or numeric metadata-name registry IDs. No fingerprint should be used
as a collision-free or persistent chemical identity.

Native metadata equality already compares exact types, exact values and full
unit identity. The source scalar DataValue epsilon comparison and its rounded
hash are not introduced. In particular, nearby unequal native doubles are
hashed from their actual bits. Strings and typed scalar/list alternatives that
happen to have identical Display text remain different hash inputs.

No new `Eq` implementation is added. `Unit` retains its existing valid `Eq`.
`Product` has public floats and `MetaValueData` can contain unvalidated NaNs;
their equality remains partial. Hashing an invalid Product or raw data enum
does not invoke validation or turn NaN into a reflexive value. `MetaValue`
continues rejecting nonfinite floats at construction. Records without `Eq`
remain unsuitable as direct `HashMap`/`HashSet` keys: callers must explicitly
choose an appropriate validated or bitwise key model.

The existing native unit representation is unchanged: a CVTerm uses its
MetaValue's optional Unit, whereas source CVTerm and source DataValue have
independent unit fields. Empty-unit absence is `None`; checked Unit construction
requires a nonempty valid accession. Hashing preserves every state the native
model represents; it does not add the source's independent-unit states.

## Header coverage review

- `METADATA/Product.h`: all own public operations have native equivalents after
  this change: construction/copy, equality, scalar access/mutation and hashing.
  CV/ordinary metadata are preserved through the existing owned CVTermList.
  Native validation and hash-value differences remain as documented.
- `METADATA/CVTermList.h`: its class mutation/query/equality/hash operations are
  represented. The free source `hashCVTerms` / `hashCVTermList` helpers omit
  ordinary metadata; the native equivalent for that use is to hash `.terms()`
  directly. This is a native equality-compatible hash, not their source digest
  or exact field-omission fingerprint. Inherited numeric MetaInfo registry APIs
  remain a separate source API boundary.
- `METADATA/CVTerm.h`: own field/query/equality/hash operations are represented,
  but full source state equivalence is still limited by the deliberately unified
  unit representation described above. Do not label this an unrestricted full
  source-model port merely because its native value can now be hashed.
- `DATASTRUCTURES/DataValue.h`: remains partial. This hash addition does not
  supply the full ParamValue conversion, numeric-cast family, source list-size
  ordering, configurable source string/stream rendering, or independent numeric
  unit-code API. Existing documented exact/finite native metadata policies stay.
- `CONCEPT/HashUtils.h`, `MetaInfo.h` and `MetaInfoInterface.h`: this work does
  not implement or certify their full public headers, source hash primitives or
  process-wide numeric metadata-name registry.

## Verification

[`metadata_hash.rs`](../tests/metadata_hash.rs) uses the source Product hash-test
input shape and checks every equality-significant field through a recording
Hasher. It also covers all seven alternatives, scalar/list/nested signed zero,
unit presence and each unit field, list boundaries/order, duplicate/empty CV
buckets, ordinary metadata, insertion-order-independent map equality, lossy
Display collisions and NaN-capable public records. Tests do not pin a source or
Rust numeric digest. Existing [`metadata.rs`](../tests/metadata.rs) checks the
unchanged equality, validation and mutation behavior.

[`metadata_hash_provenance.json`](../tests/data/metadata_hash_provenance.json)
records all inspected source hashes. No C++ build or execution was used.
