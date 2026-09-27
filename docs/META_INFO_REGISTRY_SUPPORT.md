# MetaInfoRegistry

`metadata::meta_info_registry` covers `METADATA/MetaInfoRegistry.h` and
`METADATA/MetaInfoRegistry.cpp` at SDK `bc9cc12`.

Tests: [`tests/meta_info_registry.rs`](../tests/meta_info_registry.rs), every
`START_SECTION` of `MetaInfoRegistry_test.cpp` (16, the two constructors and the
overload pairs gathered where they share a fixture), plus a unit test of index
exhaustion in the module.

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `MetaInfoRegistry()` | `MetaInfoRegistry::default()` | The thirteen reserved names, next index 1024. |
| copy constructor, `operator=`, destructor | `Clone`, `clone_from`, `Drop` | |
| `UInt registerName(name, description = "", unit = "")` | `register_name(&mut self, &str, &str, &str) -> Result<u32>` | See *Native differences* for the `Result`. |
| `void setDescription(UInt, ...)` | `set_description_at(&mut self, u32, &str) -> Result<()>` | |
| `void setDescription(const std::string&, ...)` | `set_description(&mut self, &str, &str) -> Result<()>` | |
| `void setUnit(UInt, ...)` | `set_unit_at(&mut self, u32, &str) -> Result<()>` | |
| `void setUnit(const std::string&, ...)` | `set_unit(&mut self, &str, &str) -> Result<()>` | |
| `UInt getIndex(const std::string&) const` | `index_of(&self, &str) -> Option<u32>` | `UInt(-1)` becomes `None`; the value is public as `UNREGISTERED_INDEX`. |
| `std::string getName(UInt) const` | `name_at(&self, u32) -> Result<&str>` | |
| `std::string getDescription(UInt) const` | `description_at(&self, u32) -> Result<&str>` | |
| `std::string getDescription(const std::string&) const` | `description(&self, &str) -> Result<&str>` | |
| `std::string getUnit(UInt) const` | `unit_at(&self, u32) -> Result<&str>` | |
| `std::string getUnit(const std::string&) const` | `unit(&self, &str) -> Result<&str>` | |
| private `next_index_` and the four maps | private fields | One map from name to index and one from index to (name, description, unit), which is the source's four maps with their shared key. `BTreeMap` rather than `unordered_map`; only lookups are exposed, so order is not observable. |

Native additions: `FIRST_REGISTERED_INDEX` (1024), `UNREGISTERED_INDEX`
(`u32::MAX`), `RESERVED_NAMES` (the constructor's thirteen entries), `len`,
`is_empty`, and `PartialEq`/`Eq`.

Every `Exception::InvalidValue` the header documents is `Error::InvalidValue`
with the source's message: "Unregistered index!" for an index, "Unregistered
name!" for the two setters by name, and "Unregistered Name!" — capital N, as
the source spells it — for the two getters by name.

## Preserved source conventions

- **Re-registering a name returns its first index and ignores the new
  description and unit.** Only the setters change them.
- **Indices count up from 1024 in registration order**; the reserved indices 1
  to 13 carry the source's descriptions verbatim and empty units.

## Native differences

- **The registry is caller-owned, not a process-wide singleton.** The source
  reaches one static registry through `MetaInfo::registry()` and serialises
  every member with `#pragma omp critical (MetaInfoRegistry)`. This crate has no
  mutable global state: lookups take `&self`, registration `&mut self`, and a
  registry shared between threads is put behind a `Mutex` (or `RwLock`) by its
  owner. The class test's `[EXTRA] multithreaded example` — one million
  concurrent registrations of the same new name, each answered with 1027 — is
  reproduced that way.
- **Index exhaustion is an error.** The source's `UInt` counter would wrap after
  `UINT_MAX`, handing out the sentinel, then 0 and the reserved indices again.
  `register_name` refuses with `Error::InvalidRange` once the next index would
  be `UNREGISTERED_INDEX`, leaving the registry unchanged.

## Source defect candidate

The throwing members raise `Exception::InvalidValue` from **inside** an
`omp critical` block (`setDescription`, `setUnit`, `getDescription(UInt)`,
`getUnit(UInt)`, `getName`). OpenMP requires a structured block to be left
normally, so in an OpenMP build an exception escaping the critical section is
not conforming, and it cannot be relied on to reach the caller's `catch`. This
was read from the source, not executed. The name-keyed getters avoid it
deliberately (their comment says the lookup "has to be outside the OpenMP
'critical' block"); the index-keyed ones do not. Not reproduced: the Rust
getters return the error.

## Evidence

Tier 3: every class-test literal is transcribed unchanged, including the
multithreaded section's `1027 * nr_iterations`. The error cases, the
ignored-arguments rule and the shared counter after a copy are derived from
`MetaInfoRegistry.cpp` and marked at their tests. No C++ was executed.
