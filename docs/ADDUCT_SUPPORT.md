# Molecular adducts

`chemistry::AdductInfo` ports the pinned core's complete adduct helper: notation
parsing, checked component construction, mass/m/z conversion, monoisotopic or
average mass shift, formula compatibility, getters and record equality.

```rust
use openms::chemistry::{AdductInfo, EmpiricalFormula};

let glucose = EmpiricalFormula::parse("C6H12O6")?;
let adduct = AdductInfo::parse("2M+Na;1+")?;
let observed_mz = adduct.mz(glucose.mono_mass())?;
let monomer_mass = adduct.neutral_mass(observed_mz)?;
```

The atomic adduct formula is neutral. Its signed ion charge is stored separately
because ionization removes or adds electrons; `EmpiricalFormula` charge instead
represents protonation. The constructor rejects charged formulas, zero charge,
zero molecular multiplier and the unrepresentable magnitude of `i32::MIN` charge.

## Mass conventions

For monomer mass `M`, adduct atomic mass `A`, signed charge `q` and multiplier `n`:

```
observed_mz = ((M*n + A) - q*ELECTRON_MASS_U) / abs(q)
monomer_mass = ((observed_mz*abs(q) - A) + q*ELECTRON_MASS_U) / n
mass_shift = A - q*(PROTON_MASS_U + ELECTRON_MASS_U)
```

`mass_shift(false)` uses monoisotopic adduct mass; `mass_shift(true)` substitutes
average adduct mass while retaining the same proton/electron correction. The
source hydrogen atomic mass differs slightly from the separately tabulated
proton-plus-electron mass, so `M+H;1+` has a small nonzero mass shift. It is not
rounded to zero. Source operation order is retained.

Finite zero and negative scalar inputs are accepted, matching the source
implementation. Nonfinite inputs or intermediate overflow return errors. A
constructor name can be any text within the size bound; it need not be parsable
notation. All fields remain immutable after construction.

## Notation and identity

The text grammar uses one semicolon: `M+H;1+`, `M-H;1-`, `M+2K-H;1+`, or
`2M+CH3CN+Na;1+`. `FromStr` delegates to `parse`. Spaces, tabs, linefeeds and
carriage returns are removed; the resulting spelling becomes `name()`. Other
control or Unicode whitespace is not silently normalized. Bracketed forms such
as `[M+H]+` are outside this source grammar.

The final charge sign controls polarity even when the magnitude has a leading
sign. For example, `M;-2+` and `M;+-2+` both have charge +2. Molecular multipliers
and stoichiometric coefficients use the source's checked signed-i32 token range;
the typed constructor accepts the full nonzero u32 multiplier range.

Every plus/minus in the formula part separates operations. A numeric-only term
has an empty formula after its coefficient is removed: `M+H-1;1+` retains H.
Zero coefficients are valid, but the following formula must still parse.
Adjacent operators, leading/trailing operators, percent delimiters, unknown
elements and invalid integer tokens return errors.

Existing isotope notation and the D alias are accepted within fragments. Named
records compare normalized name, formula, charge and multiplier. Equivalent
chemistry with different retained spelling can therefore compare unequal.

## Candidate compatibility and limits

`is_compatible(candidate)` applies the source predicate
`candidate.contains(-adduct_formula)`. It ignores candidate charge and molecular
multiplier. A dimer losing two H still requires two H in the supplied candidate;
the method does not first double that candidate. Isotope labels are distinct
elements for this comparison. Signed candidate counts retain the native
formula's signed containment behavior.

Parsing and names are limited to 65,536 input bytes. Parsing accepts at most
4,096 additive terms and precharges at most 1,000,000 conservative text/formula
work units before repeated map operations. Existing native formulas use checked
i32 atom counts. Compatibility with an adduct loss of i32::MIN returns false,
because no native candidate can contain the required 2,147,483,648 atoms.

The [independent review](ADDUCT_REFERENCE_REVIEW.md) and
[provenance](../tests/data/adduct_provenance.json) retain literal source examples
and tighter checks from the original atomic constants. The
[workflow tests](../tests/adduct_workflow.rs) compare a doubly charged sodium
adduct with its complete ion composition and preserve masses, charge and adduct
metadata through mzML. No C++ runtime or new dependency is required.
