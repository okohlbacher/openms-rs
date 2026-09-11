# Unique IDs and UUID generation

`concept::UniqueIdGenerator`, `UniqueId` and `HasUniqueId` provide native
counterparts of `UniqueIdGenerator.h` and `UniqueIdInterface.h` at SDK revision
`54a232fe2cae9c590d5c997fa49d20e7769860fb`.

```rust
use openms::concept::{HasUniqueId, UniqueId, UniqueIdGenerator};
let mut generator = UniqueIdGenerator::from_seed(546_666_321);
let mut value = UniqueId::default();
assert_eq!(value.ensure_unique_id(&mut generator), 1);
assert_eq!(value.unique_id(), 4_039_984_684_862_977_299);
assert_eq!(value.ensure_unique_id(&mut generator), 0);
value.set_unique_id_from_str("feature_00067890")?;
assert_eq!(value.unique_id(), 67890);
# Ok::<(), openms::Error>(())
```

## Generator mapping

| Source operation | Native operation |
| --- | --- |
| Once-initialized singleton | Caller-owned `UniqueIdGenerator::new` / `Default` |
| `getUniqueId` | `get_unique_id(&mut self)` |
| `getUUID` | `get_uuid(&mut self)` |
| `setSeed`, `getSeed` | `set_seed`, `seed`; deterministic `from_seed` constructor |
| Global synchronization | Caller-owned `Mutex<UniqueIdGenerator>` when shared draws are needed |

The default seed is system time in microseconds since the Unix epoch XOR the
process ID shifted left 32 bits. Explicit seeds reproduce the source MT19937-64
stream. One ID uses one raw engine word; no distribution mapping, collision
search or zero-redraw step is inserted. The source can theoretically emit zero,
which its ID interface considers invalid, despite its method descriptions saying
that generated IDs are valid. These random identifiers are not cryptographic
secrets and cannot guarantee uniqueness.

The existing private MT19937-64 implementation is reused, including its tested
Boost state normalization. Only crate-private visibility changes are required;
no recurrence, seed arithmetic or decoy shuffling operation is changed. Its
original Boost Software License notice remains in the shared implementation.
The source's six exact words for seed `546666321` independently check this reuse;
the engine's existing multi-cycle reference tests remain applicable.

A UUID consumes two consecutive words from the same owned generator. As in the
source `memcpy`, each word is laid out in the host's native byte order, then the
version-4 and variant-1 bits are applied. Output uses lowercase hexadecimal and
the standard 36-character hyphenated layout. Raw integer sequences are portable
across architectures; UUID strings for the same seed depend on byte order.
The native owned operation obtains both words without source-global interleaving;
a caller can lock one shared generator for the complete operation.

`Clone` deliberately copies native generator state and therefore repeats its
future sequence. The source singleton is not copied; this is an explicit native
ownership facility. Initialization and each draw use fixed-size state. UUID
formatting allocates a single 36-byte string. No dependency is added.

## Common value interface

`UniqueId(pub u64)` is the standalone value counterpart. It supports standard
copying, equality, ordering, hashing and zero-valued default construction.
`UniqueId::INVALID` is zero and `UniqueId::is_valid(value)` tests nonzero.
`HasUniqueId` supplies these operations to it, to a plain `u64`, and to other
records implementing two accessors over their existing ID field:

| Source operation | Native operation |
| --- | --- |
| `getUniqueId`, scalar `setUniqueId` | `unique_id`, `set_unique_id` |
| Valid/invalid count queries | `has_valid_unique_id`, `has_invalid_unique_id` return bool |
| `clearUniqueId` | `clear_unique_id` returns 1 for changed, 0 for already zero |
| `swap` | `swap_unique_id` exchanges ID fields only |
| Generating `setUniqueId()` | `assign_new_unique_id(&mut generator)` always returns 1 |
| `ensureUniqueId` | `ensure_unique_id(&mut generator)` draws only if currently zero |
| String `setUniqueId` | `set_unique_id_from_str` |
| Base-class mutable storage | `unique_id_mut` returns an exclusive `&mut u64` |

The parser takes the suffix after the last underscore, or the complete string
when no underscore exists. Only ASCII digits are accepted; leading zeros are
allowed. Empty or invalid suffixes clear the current ID. It uses defined wrapping
unsigned arithmetic just as C++, so `18446744073709551616` becomes zero and the
following integer becomes one. It does not trim spaces, accept signs or interpret
Unicode digits. A 1 MiB whole-input limit bounds the backward search and digit
scan; exceeding it returns an error before changing the ID. This resource limit
is the only parser behavior difference. Prefix contents otherwise have no meaning.

This maps the two classes' own operations. Other records implement the trait
individually; it does not establish complete inherited behavior for every feature,
map or peak class, nor does it emulate C++ virtual-method ABI or moved-from state.

## Validation

Six native tests include all source string-parser examples, ID mutation counts,
source RNG literals, cloned/reseeded multi-cycle state, independent little/big-endian
UUID byte references, the source 100,000-draw repeat/distinctness experiment,
1,000 UUID layouts, wrapping boundaries and four-thread caller-owned sharing
against the serial stream. Existing decoy and private engine tests run alongside.
No new C++ execution is claimed. Exact source paths and hashes are recorded in
[the manifest](../tests/data/unique_id_provenance.json).
