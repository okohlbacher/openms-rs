# Mass decomposition solver

`chemistry::MassDecompositionAlgorithm` enumerates amino-acid count compositions
whose summed mass falls in a configured absolute tolerance. It implements the
public scientific operation of source `MassDecompositionAlgorithm` at SDK commit
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. It uses a private native integer solver;
there is no C++ library, FFI, external optimizer, or new dependency.

```rust
use openms::chemistry::{
    AASequence, MassDecompositionAlgorithm, MassDecompositionOptions,
    PeptideFragmentType,
};

let peptide = AASequence::parse("DFPIANGER")?;
let mass = peptide.mono_mass_for(PeptideFragmentType::Internal, 0)?;
let solver = MassDecompositionAlgorithm::from_options(MassDecompositionOptions {
    tolerance: 0.0001,
    ..Default::default()
})?;
let compositions = solver.decompositions(mass)?;
assert_eq!(compositions.len(), 842);
# Ok::<(), openms::Error>(())
```

These are compositions, not peptide sequences: residue position and fragment
chemistry are not reconstructed. The output value API is described in
[MASS_DECOMPOSITION_SUPPORT.md](MASS_DECOMPOSITION_SUPPORT.md).

## Configuration and ownership

`new()` creates the source defaults. `from_options(options)` uses the global
modification registry; `with_registry(options, &registry)` resolves a caller's
immutable registry. Both return checked `Result`s. The five source parameters
are represented directly by `MassDecompositionOptions`:

| Option | Default |
| --- | --- |
| `decomp_weights_precision` | `0.01` Da integer-grid scaling |
| `tolerance` | `0.3` Da absolute mass error |
| `fixed_modifications` | empty vector of names/full IDs/accessions |
| `variable_modifications` | empty vector of names/full IDs/accessions |
| `residue_set` | `DecompositionResidueSet::Natural19WithoutI` |

`options()`, `alphabet()` and `diagnostics()` return borrowed views. The alphabet
contains `(char, f64)` pairs in ascending symbol order. `set_options` and
`set_options_with_registry` build a complete replacement first; a failed update
preserves the old options, masses and solver. The resolved solver stores only
owned numeric masses, options and warnings, so the registry may subsequently be
dropped. Ignored-modification warnings are returned as diagnostics instead of
being printed to stderr.

The enum includes all eight source residue-set names and supports exact
case-sensitive parsing through `FromStr` and enumeration through `ALL`.
`Natural19WithoutI` and `Natural19WithoutL` have their indicated 19 residues;
`Natural20` has 20. Despite its name, source `Natural19J` also contains those same
20 I/L residues and no J. `AllNatural` includes O and U (22 residues).
`AmbiguousWithoutX` additionally includes B, J and Z (25); `Ambiguous` and `All`
include X as well (26).

Source B/Z/X have empty formulas and zero free-residue masses. Their initial
internal masses are consequently minus water. This port preserves the initial
values while resolving fixed modifications, then requires the completed
alphabet to have positive finite masses. Thus B/Z can be repaired by suitable
fixed records, while unrepaired ambiguous sets fail explicitly. X-origin
modifications are ignored by the source, so they cannot repair X. No unrecognized
set is silently substituted.

This is a typed Rust configuration/lifecycle API. It does not reproduce inherited
C++ `DefaultParamHandler` methods or ABI. Native `Param` and `DefaultParamHandler`
remain separate reusable facilities. The five solver settings are all supported.

## Source chemistry and ordering

Natural weights are evaluated as **full amino-acid formula mass minus water**,
retaining the source operation order rather than substituting an internal-formula
mass. Existing native elemental constants are used. Source formula accumulation
can depend on pointer order; exact binary identity with an executed C++ build is
not claimed. Both literal source solver counts are verified.

Named modifications resolve in first registry/provider order when ambiguous.
This deterministic native choice replaces the source's allocation-dependent
pointer selection and does not change the registry's ordinary checked lookup
policy. Fixed and variable definitions are separately keyed and sorted by full
ID, exactly as source `ModificationDefinition::operator<`; duplicate full IDs
keep the first selected record.

Fixed modifications operate in that order. A nonzero stored absolute mono mass
replaces the weight **directly**, without water subtraction. Otherwise a nonzero
declared delta is added to the current weight; repeated fixed deltas accumulate.
A missing origin key starts at zero, as in the source map. Terminal specificity,
formula and average-mass fields are not consulted by this algorithm. This differs
intentionally from applying a modification to an `AASequence`.

Variable definitions receive labels `a` through `z` in full-ID order. Each
consumes a label even if ignored for X origin or missing mass. A nonzero stored
absolute mono mass again becomes the direct weight; otherwise the declared delta
is added to the already fixed-modified origin weight. Source zero insertion of a
missing origin is preserved and rejected if the final alphabet remains invalid.
At most 26 distinct variable definitions are accepted, replacing source iterator
exhaustion with a checked error. Equal masses with different symbols remain
separate alternatives; there is no mass-based deduplication or I/L expansion.

## Numerical selection and enumeration

The solver computes `floor(mass / precision + 0.5)` for every weight. It divides
integer weights by their common GCD and multiplies the precision by that GCD;
weights are not recalculated after scaling. Relative minimum/maximum rounding
errors are computed from these final integers and original masses.

For query mass `m`, tolerance `e`, and effective precision `p`, the integer range
is exactly:

```text
start = ceil((1 + minimum_rounding_error) * (m - e) / p)
end   = floor((1 + maximum_rounding_error) * (m + e) / p)
start <= integer_mass < end
```

The upper integer endpoint is **exclusive**, as implemented by source
`RealMassDecomposer`, even though its surrounding description suggests an
inclusive range. Zero tolerance can therefore yield no result for an exact
integer-grid composition. An exact zero query with zero tolerance similarly
returns no compositions. Precision can affect such boundary cases; it is not
promised to change only cache usage.

The private extended residue table retains source GCD/cache-block construction,
its `first_weight * last_weight` infinity sentinel, and the recursive LCM/residue
traversal order. It is not reordered by mass. The unused single-witness bookkeeping
is omitted because it cannot affect this operation. The finite infinity sentinel
is a source limitation for some unusual alphabets: weights `[6, 9, 2, 2]` use
sentinel `12` and omit the valid mass-13 composition `[0, 1, 2, 0]`. This behavior is
preserved and regression-tested; arbitrary custom weights are not promised to
produce every mathematically possible composition. Candidate mass is accumulated
in symbol order as `sum(weight[i] * count[i])` in f64. The final real filter keeps
`abs(candidate_mass - m) <= e`. Filtering happens during generation to avoid a
large intermediate vector, preserving the source's retained order and duplicates.

`decompositions(mass)` returns owned results. `append_decompositions(&mut output,
mass)` preserves existing order and appends that same sequence. It stages new
results and commits only after all checks and final reservation succeed. A failed
append leaves existing values and their allocation unchanged.

## Checked limits

Each construction and each query has a fresh shared budget. The limits are
explicit errors, not truncation or partial results:

- 50 million charged name comparisons, formula/table operations, integer-range
  visits, recursive steps, candidate mass terms and result construction work.
- 64 MiB cumulative logical allocation accounting, including lookup
  temporaries, tables, result maps and vector storage. A query also limits the
  logical payload of **existing plus new** output to 64 MiB. Count-map estimates
  include 256 bytes per nonempty map plus 64 bytes per key and the value header;
  vector allowances include growth/staging. These conservative estimates exclude
  allocator bookkeeping, fragmentation and preexisting spare capacities, and are
  not a physical RSS guarantee.
- Four million u64 residue-table cells (32 MiB before ancillary storage), checked
  before allocation. The table width is the first scaled weight after GCD.
- 100,000 total existing plus new output values; 10,000 input modification names,
  one MiB of cumulative name bytes, and 200,000 registry records when looked up.
- Positive finite precision, nonnegative finite tolerance/query mass, positive
  finite final alphabet masses and nonzero rounded weights. Scaled bounds/weights
  must lie between zero and `2^53 - 1`, avoiding ambiguous float-to-integer casts.
- Checked u64 integer products/additions and source u32 intermediate counts.
  Fixed u64 arithmetic replaces the source platform-dependent unsigned-long width. Accepted
  output counts must also fit the source constructor's signed i32 token parser.
  Alphabet size bounds recursion to at most 52 public symbols (128 private cap).

A negative rounded search endpoint, nonfinite computation, too-fine precision,
zero integer weight or exhausted budget returns an error instead of a source
invalid cast, division by zero, overflowing arithmetic or unbounded work. Fixed
invalid intermediate masses may still be overwritten before final validation.

## Evidence and scope

[The provenance manifest](../tests/data/mass_decomposition_algorithm_provenance.json)
records exact source hashes and separates source literals from native analytical
checks. The class test's `DFPIANGER` counts are 842 at tolerance `0.0001` and 911
at `0.001`. The source `Weights` GCD example is also preserved. Direct tests cover
all residue sets, modification conventions, source boundary order, caller-owned
registry lifetime, atomic reconfiguration/append and explicit resource failures.
Private tests independently brute-force seven small integer alphabets for every
mass from 0 through 120, including equal weights and nontrivial GCDs.

The standalone public APIs of `IMSAlphabet`, `IMSElement`, `Weights`,
`IntegerMassDecomposer`, and `RealMassDecomposer` are not claimed as ported by this
module. In particular, the separate IMS constraint/count/single-witness overloads
are outside `MassDecompositionAlgorithm`'s public scientific operation. They are
now implemented in the separate [public IMS solver module](IMS_DECOMPOSER_SUPPORT.md).
No C++ execution or differential-runtime validation was performed for this
peptide solver.
