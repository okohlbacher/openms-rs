# Peptide physicochemical properties

The native chemistry API implements the scientific operations of the pinned
`IsoelectricPoint`, `HydrophobicityProfile` and `AAIndex` classes, including the
seven hydrophobicity scales exposed by `Residue`. All calculations use owned
`AASequence` inputs by reference and require no optional features or dependencies.
See the [source and independent numerical review](PEPTIDE_PROPERTIES_REFERENCE_REVIEW.md)
and [fixture manifest](../tests/data/peptide_properties_provenance.json).

## API

All types below are available under `openms::chemistry`.

| API | Behavior |
| --- | --- |
| `IsoelectricPoint::default()` | Lehninger pKa scale and tolerance `1e-4`; public `scale` and `tolerance` fields |
| `compute_charge(&sequence, ph)` | Fractional net charge at any finite pH |
| `compute_pi(&sequence)` | Bounded source bisection on pH 0–14 |
| `ProteomicsPkaScale` | `Lehninger`, `Emboss`, `Sillero`, `Bjellqvist` |
| `HydrophobicityScale::value(residue)` | Checked lookup for one uppercase canonical residue; `ALL` enumerates seven scales |
| `HydrophobicityProfile::compute_gravy(&sequence)` | Mean Kyte–Doolittle value |
| `compute_profile(&sequence, scale)` | One value per residue |
| `compute_windowed_profile(&sequence, window, scale)` | Moving arithmetic mean, with no padding |
| `compute_hydrophobic_moment(&sequence, window, angle_degrees)` | Normalized Eisenberg vector magnitude per window |
| `AAIndexScale::value(residue)` | Ten checked numeric scales, enumerated by `ALL`; `accession()` gives the source AAindex identifier |
| `AAIndex::aliphatic/acidic/basic/polar(residue)` | Literal source 0/1 indicators |
| `AAIndex::calculate_gb(&sequence, temperature_kelvin)` | Estimated gas-phase basicity with checked numerical evaluation |

Profile functions and gas basicity take explicit arguments. To reproduce source
defaults, use Kyte–Doolittle with window 7, moment window 11 with angle 100°, and
gas-basicity temperature 500 K. The ten AAindex accessions are KHAG800101,
VASM830103, NADH010106, NADH010107, WILM950102, ROBB760107, OOBM850104,
FAUJ880111, FINA770101 and ARGP820102. Hydrophobicity scales are Kyte–Doolittle,
Eisenberg, Hopp–Woods, Bull–Breese, Black–Mould, Guy and Eisenberg consensus.

Run the [peptide-property example](../examples/peptide_properties.rs) to print
charge, pI, GRAVY and gas basicity for three modified or unmodified peptides:

```sh
cargo run --locked --offline --example peptide_properties
cargo run --locked --offline --example peptide_properties -- 'PEPTIDER'
```

## Chemical conventions

These source models read parent residue codes. They do not require an empirical
formula or a known mass. Hydrophobicity and gas basicity ignore all annotations.
Charge and pI suppress an N- or C-terminal group whenever that separate terminus
has an annotation; sidechain modifications retain the parent residue's pKa.
A terminal-specific record stored on a residue does not suppress the separate
terminus. Thus these utilities do not predict PTM-specific property changes.

Charge and pI treat B/Z/X/J as neutral sidechains, accept U with acidic pKa 5.73,
and reject O. Bjellqvist retains its residue-specific terminal overrides. Every
numeric AAindex scale, hydrophobicity scale and gas-basicity calculation requires
the twenty canonical residues. Source `999` hydrophobicity sentinels become
errors. Indicator functions instead return zero for any unrecognized character;
their literal sets include F/G as aliphatic and W as basic. FAUJ880111 selects
H/K/R and is distinct from that basic indicator.

Charge and pI require a nonempty sequence. Bisection checks pH 0 before pH 14,
uses the tolerance both for endpoint charge and interval width, and selects the
endpoint with smaller absolute charge when no root is bracketed. Exact zero at
an interior midpoint advances the upper bound, preserving source behavior.
Tolerance must be finite and positive for pI; charge does not use it. Finite
extreme pH values give saturated finite charges.

Hydrophobicity also requires nonempty input. Windows must be positive and clamp
to sequence length, producing `n - min(window, n) + 1` entries. Moving means
subtract the outgoing value before adding the incoming value. Moment phase
starts at zero in every window. Finite negative or zero angles are valid; angle
conversion and phase products must remain finite. Eisenberg and Eisenberg
consensus are separate scales; moments always use the normalized former.

## Gas-basicity numerical behavior

Gas basicity preserves the source gas constant, backbone/sidechain pairing,
ordered accumulation and final division by `ln(2)`. The first sidechain is
omitted; later zero-energy sidechains still contribute a unit exponential.
Temperature must be finite and positive; the empty sequence is valid.

Ordinary finite evaluation follows source floating-point order. If it becomes
nonfinite, the implementation retries an algebraically equivalent shifted
energy sum, retaining finite low-temperature results and contributions from
tied maxima. If positive `R*T` rounds to zero, its maximum-energy limit is used.
Empty input uses the exact single-site identity at every temperature, correcting
the source's finite high-temperature cancellation. These are explicit numerical
differences; they do not change the underlying empirical model. Returned values
must be finite. Independent 90-digit references cover normal and extreme
temperatures; this is not a universal correctly-rounded guarantee.

## Resource limits and verification

Each operation accepts at most 1,000,000 residues and 50,000,000 work units.
Limits are checked before allocating outputs or traversing oversized input;
work accounting includes repeated passes. pI additionally limits bisection to
128 iterations, shares work across endpoint and midpoint evaluations, and
errors if floating-point interval narrowing stalls. Calculations return owned
results and leave inputs unchanged on success or error.

Focused tests cover literal source values, all scale entries and sentinels,
independently isolated pKas, numerical boundaries, rolling windows, moments and
resource limits. [Workflow tests](../tests/peptide_properties_workflow.rs) connect
digestion, terminal annotations, formula-free mass tags and typed identification
metadata through idXML. The C++ library was not built or executed; full core
parity remains in progress.
