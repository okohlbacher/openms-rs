# IsotopeCluster

`DATASTRUCTURES/IsotopeCluster.h` at core `bc9cc12` (header-only) is ported in
`src/data_structures/isotope_cluster.rs`. The source struct stores an isotopic
cluster: a group of peaks that likely come from one peptide with different
isotopic compositions, the scans they appear in, and a charge estimate when one
is available.

## API mapping

| Source | Rust |
|---|---|
| `struct IsotopeCluster` | `IsotopeCluster` |
| `IsotopeCluster()` | `IsotopeCluster::default()` |
| member `peaks` (`ChargedIndexSet`) | `IsotopeCluster::peaks` |
| member `scans` (`std::vector<Size>`) | `IsotopeCluster::scans` (`Vec<usize>`) |
| `typedef std::pair<Size, Size> IndexPair` | `IndexPair` (`(usize, usize)`) |
| `typedef std::set<IndexPair> IndexSet` | `IndexSet` (`BTreeSet<IndexPair>`) |
| `struct ChargedIndexSet : public IndexSet` | `ChargedIndexSet`: the base as the field `indices`, plus `Deref`/`DerefMut` |
| `ChargedIndexSet()` | `ChargedIndexSet::default()` (charge 0) |
| member `charge` (`Int`) | `ChargedIndexSet::charge` (`i32`) |

The three nested types are module items (Rust has no nested type
declarations); `analysis::feature_finder_picked::defs` re-exports them, since
`FeatureFinderDefs` names the same types.

## Preserved source conventions

- A default cluster is empty with charge 0, which the source documents as "no
  charge estimate".
- `IndexSet` orders by scan index, then peak index, as `std::set` of
  `std::pair` does.
- `ChargedIndexSet` compares and orders by its indices only; the charge takes
  no part, as the operators inherited from `std::set` do.

## Native differences

- `IsotopeCluster` has no equality operator, as the source has none; it
  derives `Clone` (the implicit copy constructor) and `Debug`.

## Checked boundaries and evidence

There is no `IsotopeCluster_test.cpp` at `bc9cc12`. `tests/isotope_cluster.rs`
checks the documented defaults (tier 3) and the ordering the typedefs name
(tier 4); the executed `FeatureFinderDefs` probe rows for `ChargedIndexSet`
are in `tests/feature_finder_picked_seeds.rs`. The type holds caller-built
data and performs no work of its own.
