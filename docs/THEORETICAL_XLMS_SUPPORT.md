# XLMS theoretical spectra

This module implements the complete class-specific public operation group of
OpenMS4-core `TheoreticalSpectrumGeneratorXLMS` at
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`: all three generation overloads,
`LossIndex`, the 25 settings, construction and copying. It also implements the
complete `ProteinProteinCrossLink` record and its reaction type from
`OPXLDataStructs`. The rest of OPXLDataStructs and the ProForma spectrum wrappers
remain separate work. This is a native implementation; no C++ ABI or inherited
`DefaultParamHandler` facade is supplied.

Authoritative source: [generator header](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h),
[implementation](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.cpp),
[class test](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/tests/class_tests/openms/source/TheoreticalSpectrumGeneratorXLMS_test.cpp),
and [crosslink record](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h).

## API and source mapping

The public exports are `TheoreticalSpectrumGeneratorXLMS`, `XLMSOptions`,
`XLMSLimits`, `LossIndex`, `ProteinProteinCrossLink`, and
`ProteinProteinCrossLinkType` in `openms::chemistry`. Default/Clone replace the
source constructors/copy assignment. The options are read on every generation
call, replacing the source `updateMembers_` cache lifecycle.

| Source public operation | Native operation |
|---|---|
| `getLinearIonSpectrum` | `get_linear_ion_spectrum(&mut spectrum, &peptide, link_pos, frag_alpha, charge, link_pos_2)` |
| `getXLinkIonSpectrum` with precursor mass | `get_xlink_ion_spectrum(&mut spectrum, &peptide, link_pos, precursor_mass, frag_alpha, min_charge, max_charge, link_pos_2)` |
| `getXLinkIonSpectrum` with crosslink record | `get_crosslink_ion_spectrum(&mut spectrum, &crosslink, frag_alpha, min_charge, max_charge)` |
| `LossIndex` fields | `has_h2o_loss`, `has_nh3_loss` |
| `ProteinProteinCrossLink::getType` | `get_type()` |
| Record equality/hash | `Eq`/`Hash`, preserving sequence allocation identity |

All generation operations append to the supplied spectrum, sort the combined
peaks stably by m/z, and return `Result<()>`. Source default arguments are
explicit: linear charge is 1 and second link position is 0. A zero second
position means use the first position; the source does not enforce an ordered
loop-link pair. No MS level, spectrum type, precursor metadata, acquisition
metadata, identification data, or other unrelated fields are changed.

The crosslink record owns `Option<Arc<AASequence>>` for alpha/beta. Two equal
peptides in distinct allocations remain different record identities. Clone
shares these allocations and keeps them alive. Positions are `(isize, isize)`;
linker name, terminal specificities, and precursor correction are retained.
`new(mass)`, `set_cross_linker_mass(mass)`, and `cross_linker_mass()` enforce a
finite linker mass while permitting negative values. Signed zeros compare and
hash equally. Raw source NaN/Inf records are outside this native key domain;
C++ hash numbers, pointer addresses, and layout are not compatibility promises.
A nonempty beta classifies Cross; otherwise second position -1 classifies Mono;
all other states classify Loop, including the default null/null record. The
source enum sentinel is represented by the `COUNT` constant, not an invalid
native reaction type.

## Settings

| Setting | Default |
|---|---|
| add_isotopes | false |
| max_isotope | 2 |
| add_metainfo, add_charges | true |
| add_losses | false |
| add_precursor_peaks | true |
| add_abundant_immonium_ions | false |
| add_k_linked_ions, add_first_prefix_ion | true |
| add_a_ions, add_b_ions, add_y_ions | true |
| add_c_ions, add_x_ions, add_z_ions | false |
| a_intensity, b_intensity, c_intensity, x_intensity, y_intensity, z_intensity | 1.0 each |
| relative_loss_intensity | 0.1 |
| precursor_intensity, precursor_h2o_intensity, precursor_nh3_intensity | 1.0 each |

The source TODO flags `add_abundant_immonium_ions` and `add_first_prefix_ion`
remain inert. Isotopes add no companion when max_isotope<2 and exactly one
when max_isotope>=2, including values above 2. All ten intensity values are
validated as finite on each call, including inactive series. Negative finite
intensities remain supported; values that overflow the stored f32 are errors
when a peak is emitted.

## Scientific behavior and known source defects

Charges are visited first, then B,Y,A,X,C,Z series. Each ladder retains source
addition/subtraction order, terminal delta handling and ion corrections. The
private mass preparation uses free residue formula/mass minus water, matching
`Residue`; it does not substitute the ordinary TSG algorithm or AASequence's
internal-composition accumulation. Modification precedence is source-specific:
no-change gate, delta formula, absolute formula, declared absolute mono mass,
then delta mono mass. Native empirical formula accumulation retains the existing
deterministic element order, so exact C++ binary64 identity across pointer
orders is not claimed. B/Z/X have source zero free mass and therefore negative
internal mass; full/span queries reject bare X while direct residue ladder
access remains available. Known and numeric modifications use their retained
mass semantics.

Water losses occur for D/E/S/T and ammonia losses for K/N/Q/R, accumulated by
parent residue identity. The complete source loss alphabet is
`RHKDESTNQCUGPAVILMFYW`; modifications' own neutral-loss lists are ignored.
Linked losses also include the other peptide's flags. A linked ladder skips
an unavailable loss index while retaining its base/isotope output, including
the peptide-length endpoint. Empty or unsupported loss preprocessing returns a
checked error in place of source empty-vector access or map lookup failure.

Ordinary fragment peaks skip negative m/z but retain zero. Losses require
strictly positive neutral mass before division. K-linked peaks use intensity1,
subtract prefix-B then suffix-X mass, and annotate the opposite partner.
Precursor, water-loss and ammonia-loss peaks are always emitted together when
enabled, once at max_charge, independent of loss eligibility. Their negative
positions are retained. Isotope companions share their base annotation and
intensity; loss fragments have no isotope companions.

The backend deliberately preserves two finite source scientific defects. Tests
keep the literal source result separate from the independently derived chemical
expectation:

- **CPP-042:** linear suffix losses divide twice. For AS y1 at charge2,
  y1=53.528573578421; source y1-H2O=17.7590042573105, while the chemical value is
  44.523291046521. Source cpp311 passes divided m/z to helpers570/589.
- **CPP-043:** precursor companions add the isotope spacing to an undivided
  charged mass. For precursor1000 at charge2, the base is501.007276466771;
  source companion=1002.516230352442, while the chemical value is501.508953885671.
  Source cpp623/653/682 contains the three affected expressions.

These are source inspection plus independent scalar algebra, not executed C++
spectrum comparisons. The central [C++ issue log](../OpenMS_CPP_ISSUES.md)
records both defects. CPP-038 concerns the separate, unported ProForma wrapper;
this backend adds the supplied linker once and makes no compensation for that
wrapper's peptide modification behavior.

A missing alpha record returns immediately. An empty alpha suppresses both
alpha and beta ladders, while source K-linked and precursor branches remain
separate. Record terminal flags, name and precursor correction do not affect
this algorithm. Empty beta is omitted from precursor accumulation, preserving
`(alpha_full + linker_mass) + beta_full` order only when beta exists.

## Checked boundaries and atomic append

Finite signed charges are supported. Linear charge<=0 runs no charge iterations
but still performs applicable array creation/renaming/sorting. An inverted
linked interval may still emit precursors at max_charge. Zero charge fails only
when a chosen branch divides by it. Widened/checked charge-span arithmetic
avoids source loop overflow and 32-bit conversion wrap. Invalid consumed
positions, short C/X peptides and unrepresentable/nonfinite intermediate masses
return errors; safe endpoint and dormant-loop behavior remain available.
Missing-alpha no-op precedes configuration validation, as in source.

The first existing integer/string arrays are used and renamed `charge` and
`IonNames` when enabled. Their descriptions, metadata and shared processing
identities are moved intact. Missing/empty enabled arrays pad old peaks with
neutral annotations (0/empty string), a documented repair of invalid source
array lengths. Nonempty unrelated arrays cannot be extended with invented
values and cause an error; empty placeholders are retained. Disabled arrays
are still permuted if aligned and no new rows invalidate them.

Every fallible step precedes publication: mass preparation, generated rows,
final stable permutation and all array data. An error leaves the complete
spectrum unchanged. Existing acquisition/identification/metadata graphs are
preserved by ownership and are neither cloned nor destructed. Data and replaced
names are charged before copying/destruction.

Default limits are 4096 residues per peptide, 100000 combined peaks, 50 million
work units and 256 MiB cumulative logical allocation allowance. `XLMSLimits`
can adjust them. Both peptides, all charges, formatting, vectors/formulas and
sorting share one allowance. Charge loops are conservatively charged using
`charges * (residues + 1) * 16` before execution, even when negative peaks,
disabled series or empty-alpha ladder guards would avoid work. These estimates
are deterministic checked bounds, not measurements of physical heap usage or
RSS. Late count, work, allocation or finite-value errors never publish partial
results.

## Evidence and reproduction

[The direct generator tests](../tests/theoretical_xlms.rs) retain all three
upstream rounded mass arrays (52 values), six allowed-name sets (113 strings),
source settings/counts, charge distributions, small loop cases and isotope
saturation. These name sets are source membership expectations, not promises
that every allowed name is emitted. Independent tests cover terminal/numeric
and custom modifications, signed charges, the two logged defects, linked
endpoints, empty alpha, append ownership and atomic failures.
[Record tests](../tests/protein_cross_link.rs) cover every equality field,
reaction classification and owned allocation identity. Private tests verify
shared mass/charge budgets and late output rollback.

Run the dependency-free extractor against the exact source checkout:

```text
python3 tools/generate_xlms_reference.py /path/to/OpenMS4-core --check
```

It verifies the pinned class-test SHA before extracting unchanged decimals and
strings. It performs no C++ build/execution. The [provenance manifest](../tests/data/xlms_provenance.json)
records all 11 source hashes, fixture hashes and generator hash. The native
focused validation results are recorded in the integration handoff after both
compiler runs; source fixtures alone do not establish executed C++ parity.
