# IMS isotope distributions and elements

`IMSIsotopeDistribution` and `IMSElement` implement the public operation surfaces
of the corresponding OpenMS4-core IMS classes at
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. They use owned Rust values and explicit
operation settings. They are separate from the general
[isotope distribution API](ISOTOPE_SUPPORT.md) and fine isotope generators:
IMS bins store **mass defects**, with one abundance per nominal-mass offset.
This increment does not claim an IMSAlphabet or IMS molecule implementation.

```rust
use openms::chemistry::{
    IMSElement, IMSIsotopeDistribution, IMSIsotopeOptions, IMSIsotopePeak,
};

let hydrogen = IMSIsotopeDistribution::from_peaks(vec![
    IMSIsotopePeak { mass: 0.0078250319, abundance: 0.999885 },
    IMSIsotopePeak { mass: 0.01410178, abundance: 0.000115 },
], 1)?;
let settings = IMSIsotopeOptions { size: 3, abundances_sum_error: 0.0 };
let hydrogen_pair = hydrogen.pow(2, settings)?;
let element = IMSElement::from_distribution("H2", hydrogen_pair)?;
let isotope_zero_mass = element.mass(0)?;
let singly_ionized_mass = element.ion_mass(1)?;
# Ok::<(), openms::Error>(())
```

## Public operation mapping

| Source distribution surface | Native operation |
| --- | --- |
| `Peak`, peak/mass/abundance containers and iterator aliases | `IMSIsotopePeak`, `Vec`, immutable slices and their iterators |
| Empty/nominal constructor | `Default`, `new(nominal_mass: u32)` |
| Single-isotope constructor | `from_mass(f64)` |
| Peak-container constructor | `from_peaks(Vec<IMSIsotopePeak>, u32)` |
| Copy, assignment, equality/inequality, destruction | `Clone`, moves/assignment, `PartialEq`, normal Rust ownership |
| `SIZE`, `ABUNDANCES_SUM_ERROR` | `IMSIsotopeOptions { size, abundances_sum_error }`, passed by value |
| `size`, `empty` | `size(options)`, `is_empty()`; `stored_len()` exposes retained storage |
| Mass/abundance access | `mass(index)`, `abundance(index)`, `peaks()` |
| Nominal mass access/update | `nominal_mass()`, `set_nominal_mass(u32)` |
| `getAverageMass`, vector getters | `average_mass()`, `masses(options)`, `abundances(options)` |
| `normalize` | `normalize(options)` |
| Distribution `operator*=` | `convolve_assign(rhs, options)`, plus pure `convolve` |
| Observable source right-side padding | Explicit `convolve_assign_with_padded_rhs(&mut rhs, options)` |
| Integer `operator*=` | `pow_assign(power: u32, options)`, plus pure `pow` |
| Stream output | Checked `to_text(options)` |

| Source element surface | Native operation |
| --- | --- |
| Empty constructor | `Default` |
| Name plus nominal mass, single mass or distribution | `new(name, nominal)`, `from_mass(name, mass)`, `from_distribution(name, owned_distribution)` |
| Name and sequence access/update | `name`, `sequence`, `set_name`, `set_sequence` |
| Nominal, indexed, average, ion mass | `nominal_mass`, `mass(index)`, `average_mass`, `ion_mass(electrons: i32)` |
| Isotope distribution access/replacement | `isotope_distribution`, `set_isotope_distribution(owned_distribution)` |
| Electron constant | `IMSElement::ELECTRON_MASS_IN_U` |
| Copy, assignment, equality, destruction | Standard `Clone`, moves, `PartialEq`, owned destruction |
| Stream output | `to_text(options)` |

Source defaults for an indexed mass or ion calculation are supplied explicitly:
use index `0` or electron count `1`. Empty nominal construction may use nominal
`0`. The native element is a value, without the source virtual-destructor ABI or
C++ inheritance. All numerical public operations are available.

## Storage, truncation and arithmetic

Both source statics start at zero. `IMSIsotopeOptions::default()` therefore sets
`size = 0` and `abundances_sum_error = 0.0`. Settings belong to each operation;
there is no shared mutable configuration or hidden cross-thread state.

Construction retains every bin and does not normalize. `size(options)` returns
`min(stored_len, options.size)`, whereas `is_empty` checks actual storage. Thus a
single-isotope distribution can have accessible size zero and still be nonempty.
The raw indexed mass is evaluated in the source order:
`(stored_mass_defect + nominal_mass) + index`. Indexed access ignores truncation.
Average mass is the ordered sum of `mass(index) * abundance(index)` over **all**
stored bins, including a hidden tail; it does not divide by the abundance sum.
Equality likewise includes nominal mass and the complete stored tail.

`normalize` sums all stored abundances in input order. It scales only when
`sum > 0` and `abs(sum - 1) > abundances_sum_error`; equality with the error bound
does not scale. Finite negative abundances and negative error bounds remain
supported. Nonpositive sums remain unchanged. A running sum overflowing to
positive infinity gives scale zero and retains the source's finite signed-zero
results. A very small positive sum whose reciprocal makes a final abundance
nonfinite instead returns a checked error without modifying the distribution.

Convolution preserves the left-forward/right-backward pairing and floating-point
operation order of the source. Each destination abundance is a sum of pair
products. Its defect is the corresponding abundance-weighted defect sum divided
by that abundance, or exactly zero when abundance is zero. A nonfinite weighted
numerator discarded by this zero-bin branch does not cause an error. The output
has exactly `options.size` bins, adds the two nominal masses, and then applies
normalization to the resulting bins.

The source treats an empty right operand as a no-op and an empty left operand as
a complete copy of the right operand. These branches occur before configuration,
truncation, normalization and nominal addition. They retain long tails. Two
nonempty operands with size zero instead produce empty storage with the summed
nominal mass. The native methods preserve these distinctions, and validate only
settings actually consumed by the selected operation.

Source convolution physically pads both operands with zero bins, even though its
right parameter is declared const. The ordinary native operations preserve the
borrowed right operand and use virtual zero padding. Call the explicitly named
mutable-right operation when this representation side effect matters; it pads
the right side only after a successful nonempty fold, never truncates its tail,
and commits both outputs atomically. Self-folding is available through `pow` or
by borrowing the same value twice to pure `convolve`.

Powers follow the source binary, least-significant-bit-first multiplication
scheme and normalization after each fold. **Power zero is a source no-op**, as
is power one. It is not a mathematical unit distribution. Empty inputs and size
zero therefore have unusual but defined outcomes: an empty nominal-9 input at
power two becomes the default empty nominal-0 value; at power three it remains
empty nominal-9. A nonempty input at size zero similarly yields default empty at
power two and retains the original input at power three.

## Elements and text

A named constructor initially copies its name into the separate sequence field.
Changing either label later leaves the other untouched. Equality compares both
labels and the complete isotope distribution. Labels are arbitrary UTF-8,
including whitespace and NUL; no element-database lookup or chemical-name
validation is implied.

Ion mass subtracts the signed electron count times the literal IMS constant
`0.00054858` from isotope-zero mass. This intentionally preserves the older IMS
constant instead of substituting another native chemistry constant. It accepts
the full signed `i32` electron range and returns signed finite masses.

Text uses the existing ParamValue classic-locale stream formatter: six significant
digits, space-separated mass/abundance and one newline per accessible bin.
Element text includes the source `name:`, `sequence:`, and `isotope distribution:`
labels, followed by an extra newline. Formatting omits hidden bins. C++ stream
locale, caller-adjusted precision and stream flags are not emulated; `to_text` is
the deterministic default stream representation.

## Checked boundaries and limits

Stored mass defects and abundances must be finite; signed finite values are
accepted. Every returned numerical value and final stored fold result must also
be finite. An invalid index returns an error. Nominal mass is portable `u32`,
matching the source's unsigned-int nominal field on supported SDK platforms;
addition overflow returns an error. This is a deliberate change from the
source's **defined modulo-2^32 wrapping**, not a correction of undefined behavior.

Limits apply before expensive loops and allocation:

- At most 1,000,000 stored bins, and the same maximum requested size.
- At most 50,000,000 work units per operation, shared across all folds and copies
  of a power. Convolution precharges eight units per triangular bin pair, then
  output and normalization work.
- At most 64 MiB of cumulative **logical allocation accounting**, including
  intermediate vectors and copies. This bounds requested payload, not allocator
  bookkeeping or process RSS. Incoming owned vectors are not copied by constructors.
- At most 8 MiB of text, conservatively preflighted; scalar rendering additionally
  charges 1,024 work units and logical bytes per bin for both formatter calls.
  These shared limits can reject output before its final text alone reaches 8 MiB.
- Element name and sequence are each limited to 1 MiB. Their copies are bounded
  before allocation. Element rendering accounts both the intermediate isotope
  text and its copy into the final output.

Owned immutable storage means ordinary `Clone`, equality and destruction have
bounded input sizes; they retain ordinary Rust allocation behavior. Numerical
mutation, label setters and mutable-right convolution leave the original values
unchanged on checked errors. Allocation failure in explicit computational vector
and text reservations is returned as an error. No custom allocator or dependency
is introduced.

## Evidence and validation

The [provenance manifest](../tests/data/ims_isotopes_provenance.json) records six
pinned source hashes and separates literal source assertions from independent
arithmetic. `IMSElement_test.cpp` supplies the hydrogen mass defects and
abundances, oxygen mass `15.9994`, independent-label behavior, replacement and
ionic-mass assertions. `IMSIsotopeDistribution_test.cpp` contains construction
and destruction checks, but its numerical sections are TODO; it supplies no
numerical convolution goldens.

The native suites contain 13 distribution tests and five element tests. They
include all substantive element class-test assertions, independently derived
dyadic folds, nine small Cartesian species enumerations, binomial powers 2–12,
source empty/size-zero/power-zero cases, raw-tail behavior, padding, signed and
nonfinite boundaries, normalization overflow, nominal overflow, atomic errors,
and shared formatting/power resource limits. These tests verify native behavior
against source reading and independent arithmetic; no C++ build or executed C++
parity is claimed.
