# IMS public mass decomposers

This module implements the public scientific operations of OpenMS Core SDK
`MassDecomposer`, `IntegerMassDecomposer` and `RealMassDecomposer` at
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. It exposes arbitrary-weight decomposition
without requiring residue names or modification registries. The peptide-specific
`MassDecompositionAlgorithm` remains separate. These are bounded source-compatible
algorithms, including known finite source quirks; they are not an assurance that
an arbitrary unsorted alphabet yields every mathematical solution.

Use [IMSWeights](IMS_WEIGHTS_SUPPORT.md) when converting real masses to integers.
Its precision, rounding and optional GCD division remain explicit caller choices.
[IMSAlphabet](IMS_ALPHABET_SUPPORT.md) can supply named real masses, but is not a
required dependency of the decomposer implementation.

```rust
use openms::chemistry::{IMSIntegerMassDecomposer, IMSRealMassDecomposer, IMSWeights};

let integer = IMSIntegerMassDecomposer::from_integer_weights(&[4, 5, 6])?;
assert_eq!(integer.decompositions(10)?, vec![vec![0, 2, 0], vec![1, 0, 1]]);
assert_eq!(integer.decomposition(10)?, Some(vec![1, 0, 1]));

let weights = IMSWeights::from_masses(&[3.0, 5.0, 8.0], 1.0)?;
let real = IMSRealMassDecomposer::new(&weights)?;
assert_eq!(real.number_of_decompositions(8.0, 1.0)?, 2);
# Ok::<(), openms::Error>(())
```

This example is explanatory Rust, not an additional crate doctest.

## Public operation mapping

| Source operation | Native operation |
|---|---|
| Abstract `MassDecomposer<ValueType, DecompositionValueType>` | Object-safe `IMSMassDecomposer` trait |
| `IntegerMassDecomposer(const Weights&)` | `IMSIntegerMassDecomposer::new(&IMSWeights)` |
| Direct exact integer input | `from_integer_weights(&[u64])`, a native convenience |
| `exist(mass)` | `exists(u64) -> Result<bool>` |
| `getDecomposition(mass)` | `decomposition(u64) -> Result<Option<Vec<u32>>>` |
| `getAllDecompositions(mass)` | `decompositions(u64) -> Result<Vec<Vec<u32>>>` |
| Integer `getNumberOfDecompositions(mass)` | `number_of_decompositions(u64) -> Result<u32>` |
| `RealMassDecomposer(const Weights&)` | `IMSRealMassDecomposer::new(&IMSWeights)` |
| Real `getDecompositions(mass,error)` | `decompositions(f64,f64)` |
| Constrained overload | `decompositions_with_constraints(f64,f64,&BTreeMap<usize,(u32,u32)>)` |
| Real `getNumberOfDecompositions(mass,error)` | `number_of_decompositions(f64,f64) -> Result<u64>` |
| Copy/destruction | Owned `Clone` values and Rust destruction |

The concrete integer type also provides the four trait operations as inherent
methods. No C++ template instantiation or ABI is claimed: native masses use
portable `u64`, multiplicities and integer counts use `u32`, and real counts use
`u64`. This matches default C++ integer widths on LP64 Linux/macOS. C++
`unsigned long` is 32 bits on Windows; arbitrary source template widths are not
exposed. A native `usize` constraint index replaces C++ `unsigned int` indexing.

`None` denotes a missing single result. A zero-mass solution is present and has
one zero per alphabet entry. The immutable decomposer retains its own numeric
state; changing or dropping the original `IMSWeights` cannot change queries.

## Integer arithmetic, ordering and source quirks

Construction uses positive weights in caller order. It performs no sorting,
GCD reduction, normalization or duplicate elimination. Duplicate weights remain
distinct count positions. The first weight determines residue-table width.
Enumeration uses the source LCM residue-class traversal; the last coordinate
is not necessarily a simple ascending outer loop. For `[6,9,20]` at mass 60,
the order is `[10,0,0]`, `[7,2,0]`, `[4,4,0]`, `[1,6,0]`, `[0,0,3]`.

Existence consults the final residue table. The single-result method follows
its witness vector with the source strict/non-strict comparisons and counter
updates. It does not return the first enumerated solution: for `[4,5,6]` at
10 it returns `[1,0,1]`, while enumeration starts with `[0,2,0]`. The source
single-result loop also has an early break after adding a witness if that
witness exceeds the remaining mass. The native implementation preserves this
branch instead of substituting a different composition. A zero-progress
witness is a checked error in place of an unbounded loop.

A concrete nonprogressing source witness occurs for `[10,6,15]`. The final
residue row is `[0,21,12,33,24,15,6,27,18,39]`, but witness 3 is `(index 2,
count 0)`: the gcd-5 cache loop increments only its first residue's counter
before updating several residues. At mass 33, existence is true and enumeration
returns `[0,3,1]`, while source single-result traversal repeatedly subtracts
zero. Native single-result queries return an error for this case and masses
`33 + 10k`, preserving the other methods' results. The source rows and zero
witness have a dedicated regression; the exhaustive oracle is retained.

The source uses `first_weight * last_weight` as a finite infinity sentinel.
This can cause missing solutions and disagreement between methods. For
`[6,9,2]` at 13, enumeration returns `[0,1,2]` and count returns one, but
existence is false and the single-result method returns `None`. For
`[6,9,2,2]` at 13, independent exhaustive enumeration finds three solutions;
the source traversal returns only `[0,1,1,1]` and `[0,1,0,2]`, omitting
`[0,1,2,0]`. Tests retain both the mathematical oracle and the distinct source
result. Sorting an alphabet can change this behavior and the result order;
callers must also remap count indices if they choose to sort.

An empty or zero-weight alphabet is rejected before table construction.
A singleton has exact division-based existence, single, enumeration and count
operations. This is a native extension: source singleton enumeration works,
but source existence/single-result access a table that was never initialized.
Singletons need no sentinel or residue table and accept any positive `u64`.

Source unsigned arithmetic and narrowing can wrap. Native table sums, LCM
products (including multiplication before division), sentinel products,
witness counters and output multiplicities are checked. No silent wrapping,
partial-count truncation or arbitrary `2^53` restriction is introduced. The
real wrapper separately checks representable floating-point interval endpoints.

## Real intervals and constraints

For ordinary finite values the source calculation order is retained:

```text
start = ceil((1 + minimum_relative_rounding_error) * (mass - error) / precision)
end   = floor((1 + maximum_relative_rounding_error) * (mass + error) / precision)
search integer masses start <= m < end
```

Each candidate's real parent mass is accumulated in supplied alphabet order,
using the original masses and the emitted counts, then filtered by inclusive
absolute error. The exclusive integer upper bound is intentional. It can make
zero-tolerance searches empty even when a physical solution exists. It can also
exclude an exact target with a small positive tolerance: `[2.25,3.75]` at
precision 0.5 becomes `[5,8]`; target 2.25/error 0.1 has `start=end=5`.

The source real count method starts at one when `mass-error <= 0`. It is not
always equal to enumeration length, despite the source header's equivalence
claim. With masses `[2,3]`, precision 1, target 0.5 and error 0.5, enumeration
returns the all-zero vector while counting returns zero. At target zero/error
zero both return empty/zero. An endpoint such as `ceil(-1)` is rejected by
native enumeration rather than converted outside the unsigned range; the
source count overload can avoid that conversion by its start-at-one branch.

Finite signed target/error and finite signed nonzero precision are permitted
when their computed integer endpoints are valid. Negative error commonly
produces an empty interval or rejects every physical candidate; it is not
silently converted to its absolute value. Endpoints must be finite and in
`[0,2^64)` after the literal ceil/floor. Nonfinite stored/selected rounding
bounds and nonfinite parent masses are checked errors. A parent mass overflow
therefore errors instead of allowing an infinite candidate to reach the source
filter. The subtraction used for the final error comparison retains IEEE
behavior, so an infinite distance is rejected by a finite tolerance.

Constraints are inclusive and use original count indices. A lower bound larger
than its upper bound excludes all results. Out-of-range indices return a checked
error before traversal. The source header says these indices are skipped, but
the implementation indexes the vector unchecked; the native behavior is an
explicit safety boundary rather than a claim to reproduce that statement.

## Bounds and transactionality

Every construction/query has one cumulative work and logical allocation ledger:

| Bound | Limit |
|---|---:|
| Alphabet entries / maximum recursive depth | 128 |
| Residue-table cells | 4,000,000 |
| Work units per construction/query | 50,000,000 |
| Logical cumulative allocation bytes | 64 MiB |
| Retained decomposition rows | 100,000 |

The first dimension and full table product are checked before allocation.
Table initialization/copies, GCD loops, witness traversal, every recursive
state, count updates, real integer-range traversal, parent-mass calculation,
constraints and result copies consume shared work. Constraints rejecting every
candidate do not bypass that budget. Destination vectors and table/witness
scratch are charged before allocation; fallible reservations precede writes.
Byte accounting includes requested vector storage and cumulative replacement
capacity. It is a logical resource guard, not an exact operating-system heap or
allocator-overhead ceiling. Count/dimension caps separately bound retained size.

The source count implementations materialize temporary decomposition vectors,
despite the real header's opposite claim. Native counting visits the identical
candidate traversal without retaining vectors, and can exceed 100,000 results
while respecting work and checked integer counts. Enumeration uses one reusable
count vector instead of copying it at every recursion step. These allocation
adaptations do not change emitted order or arithmetic.

Queries return owned values only on complete success. Late output/work/allocation
failure drops staged results and cannot mutate the decomposer or publish a
partial result. Private tests force these failures after useful traversal and
verify that later queries still work. No append API, mutable solver configuration,
public generic solver engine or concurrent global cache is introduced.

## Evidence and integration boundary

The [provenance manifest](../tests/data/ims_decomposers_provenance.json) records
pinned headers, implementations and class tests. The upstream three class tests
exercise construction/destruction; their numerical operations are marked TODO
or not testable. The native suite reproduces the meaningful Natural19WithoutI
constructor setup at 0.01 precision after GCD division, then adds independent
Cartesian small-alphabet enumeration, dyadic real-mass filtering, duplicate and
non-coprime weights, witness/order distinctions, source sentinel omissions,
count-start differences, checked numeric boundaries and cumulative failure cases.
Derived source-loop examples are identified separately from upstream literals.
No C++ execution is used for this increment.

The private peptide solver remains unchanged. Its shared table/enumeration
conventions guided this implementation, but its residue alphabet, rounding
restrictions, output representation and lack of a witness vector differ.
A later narrow reuse can replace only that table/traversal after its existing
source tests and resource budgets are checked; this increment deliberately
avoids a cross-module refactor.
