# SDK constants and metadata keys

`openms::constants` (also `openms::concept::constants`) exposes all 38 numeric
entries and all 91 `user_param` strings from Core SDK
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. Numeric aliases and metadata key names
retain source spelling. The existing three chemistry mass names re-export the
same definitions, preserving their public paths and exact values.

The [implementation](../src/concept/constants.rs) retains the original source
literals and arithmetic expressions. These historical values are compatibility
data: for example, `M_PER_FOOT` remains the source's `3.048`, and Avogadro,
Boltzmann and other physical values are not replaced with newer SI definitions.
The isotope-spacing constant still differs slightly from rounded ElementDB
isotope masses. No existing scientific calculation changes as part of this port.

| Entries | Source units or meaning |
| --- | --- |
| `PI`, `E` | Mathematical constants |
| `ELEMENTARY_CHARGE`, `e0` | Coulombs |
| `ELECTRON_MASS`, `PROTON_MASS`, `NEUTRON_MASS` | Kilograms |
| Corresponding `_U` masses | Unified atomic mass units |
| `C13C12_MASSDIFF_U`, `ISOTOPE_MASSDIFF_55K_U` | Isotope spacing in unified atomic mass units |
| `AVOGADRO`, `NA`, `MOL` | Per mole |
| `BOLTZMANN`, `k` | Joules per kelvin |
| `PLANCK`, `h` | Joule seconds |
| `GAS_CONSTANT`, `R` | `NA * k` |
| `FARADAY`, `F` | `NA * e0` |
| `BOHR_RADIUS`, `a0` | Metres |
| `VACUUM_PERMITTIVITY` | Coulomb squared per joule metre |
| `VACUUM_PERMEABILITY` | `4 * PI * 1e-7` in joule second squared per coulomb squared metre |
| `SPEED_OF_LIGHT`, `c` | Metres per second |
| `GRAVITATIONAL_CONSTANT` | Newton metre squared per kilogram squared |
| `FINE_STRUCTURE_CONSTANT` | Dimensionless |
| `DEG_PER_RAD`, `RAD_PER_DEG` | Angular conversions |
| `MM_PER_INCH`, `M_PER_FOOT` | Source length-conversion literals |
| `JOULE_PER_CAL`, `CAL_PER_JOULE` | Energy conversions |

Source `EPSILON` is writable, so it maps to `epsilon()` and `set_epsilon(f64)`;
`DEFAULT_EPSILON` exposes its initial `1e-6`. Private atomic bit storage preserves
all values, including signed zero, infinities and NaN payloads, while preventing
data races on a mutable global. Each concurrent read sees a complete old or new
value. This is separate from `f64::EPSILON`, which means machine precision.
Algorithms with explicit tolerances retain those tolerances; this addition does
not turn them into consumers of shared epsilon. Any source algorithm that reads
that global must do so explicitly when ported.

All metadata strings borrow static storage and require no allocation or lookup
registry. They cover mobility/FWHM fields, peptide and fragment scores, target/
decoy annotations, cross-linking, SIRIUS, ion-identity networks and run metadata.
For example, `user_param::FWHM_MZ_ppm` is exactly `"FWHM_ppm"`, and
`user_param::SPECTRUM_REFERENCE` is exactly `"spectrum_reference"`.

[Two native tests](../tests/constants.rs) compare every numeric bit pattern and
string with [129 recorded values](../tests/data/constants_reference.tsv) emitted
by an actually compiled, unchanged source `Constants.h`. The
[C++ probe](../tests/data/constants_probe.cpp) uses only an empty `OpenMS/config.h`
include shim: this header contains no conditional scientific definitions and
uses no configuration macro. No source declarations or formulas are substituted.
This is a header-only reference execution, not a full C++ SDK build.
[Provenance](../tests/data/constants_provenance.json) records the source, probe,
fixture, compiler and command. Tests also cover the compatibility chemistry
exports and bit-preserving epsilon replacement/restoration.
