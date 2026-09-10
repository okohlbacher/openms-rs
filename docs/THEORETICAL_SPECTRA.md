# Theoretical peptide spectra

`openms::chemistry::TheoreticalSpectrumGenerator` generates native `MSSpectrum`
values from `AASequence`. The implementation follows [OpenMS4-core revision
7c029e8cdba6abab503708ecdd56f6ab55e38ce4](https://github.com/okohlbacher/OpenMS4-core/tree/7c029e8cdba6abab503708ecdd56f6ab55e38ce4)
and uses the Rust residue, modification and isotope modules.

```rust
use openms::chemistry::{AASequence, TheoreticalSpectrumGenerator};

let peptide: AASequence = "AS(Phospho)K".parse()?;
let generator = TheoreticalSpectrumGenerator {
    add_metainfo: true,
    add_losses: true,
    ..Default::default()
};
let spectrum = generator.generate(&peptide, 1, 2, Some(3))?;
# Ok::<(), openms::Error>(())
```

## Supported behavior

| Setting or API | Behavior |
| --- | --- |
| `TheoreticalIonSeries` | a, b, c, x, y, z, z. (`ZPlusOne`) and z' (`ZPlusTwo`) terminal fragments |
| `generate` | Creates a centroid MS2 spectrum; inclusive positive fragment charge range; precursor metadata charge supplied explicitly or inferred as maximum fragment charge plus one |
| `add_first_prefix_ion` | Includes a1/b1/c1; otherwise prefix fragments start at ordinal 2. Suffix fragments always start at ordinal 1. Full-length fragments are excluded |
| `TheoreticalIonIntensities` | Configurable a/b/c/x/y/z intensities (z variants share z), plus independent intact/H2O/NH3 precursor intensities |
| `add_losses` | One peak per distinct residue loss formula present within a retained fragment; no combinations or repeated losses |
| `add_terminal_losses` | Adds water loss to prefix fragments and water/ammonia loss to suffix fragments; requires `add_losses` and no isotope envelopes |
| `TheoreticalIsotopeModel::Coarse` | Normalized low-mass isotope envelopes for intact fragments, neutral losses and precursors, with the configured maximum number of bins |
| `TheoreticalIsotopeModel::Fine` | Distinct exact-mass configurations covering `1 - unexplained_probability`, with source f32 probabilities and no normalization |
| `add_precursor_peaks` | Intact, water-loss and ammonia-loss peaks, independently of `add_losses`; maximum fragment charge only unless `add_all_precursor_charges` is set |
| `add_internal_fragments` | Both b/a internal fragments of length 2–10, independent of terminal series selection; source start/endpoint and loss rules described below |
| `add_abundant_immonium_ions` | Seven fixed singly charged, unit-intensity peaks for eligible unmodified residues |
| `generate_for_activation` | Source CID/HCID/HCD/ECD/ETD/ETciD/EThcD presets and charge selection |
| `append_mass_spectrum` | Source compact helper: sorted `f32` m/z values, charges from the supplied maximum down to one, first and full-length fragments included |
| `add_metainfo` | Aligned `IonNames` string and `Charges` integer arrays, one entry per peak |
| `sort_by_position` | Stable m/z ordering with aligned annotations; enabled by default |
| `append_to` | Validated atomic addition to an existing spectrum, preserving its other metadata and precursor entries |

Defaults select b/y ions with unit intensity, no first prefix ion, no losses,
no isotope envelopes, no precursor peaks and no annotation arrays. Precursor
metadata is still recorded for every nonempty peptide. An empty peptide yields a
default empty spectrum or leaves an appended spectrum unchanged. The pinned
source rejects c/x generation for a single-residue peptide; this behavior is
preserved. Other series yield no fragments for a monomer.

Unmodified D/E/S/T supply H2O loss; K/N/Q supply NH3 loss; R supplies NH3,
CH2N2 and CH2NO losses. A modified residue uses its modification's declared loss
formulas **in place of** the unmodified residue's losses, matching
`Residue::setModification`. Terminal modification loss declarations are not
added automatically. Phosphoserine loss and isotope-labeled lysine are tested.

## Mass and annotation conventions

Ion formulas apply these deltas to the full retained peptide, including H2O:

| Series | Formula delta |
| --- | --- |
| a | C-1 H-2 O-2 |
| b | H-2 O-1 |
| c | H1 N1 O-1 |
| x | C1 H-2 O1 |
| y | none |
| z | H-3 N-1 |
| z. | H-2 N-1 |
| z' | H-1 N-1 |

Prefix fragments retain the original N-terminal modification; suffix fragments
retain the original C-terminal modification. Both retain applicable residue
modifications. Monoisotopic m/z uses proton mass and preserves the independently
tabulated terminal modification mass deltas, as in OpenMS. Precursor loss peaks
use the loss-subtracted formula mass, which can differ slightly from subtracting
a loss from the intact precursor's declared mass.

Coarse and fine isotope generation follow the source's explicit **neutral hydrogen**
adduct formula: add `charge` H atoms, set formula charge to zero, calculate the
envelope, then divide mass by charge. Consequently the envelope retains the
hydrogen electron mass; it is about 0.0005486 m/z above the corresponding proton
calculation for an ordinary CHNOS peptide. Formula-based envelopes also exclude
the small difference between formula and declared terminal modification masses.
The coarse module uses the lightest supported isotope as its anchor and carbon-13
spacing for nominal bins. Fine mode enumerates actual isotope-count configurations
and computes their exact tabulated masses independently; it does not split coarse bins.
See [isotope support](ISOTOPE_SUPPORT.md) for probability and isotope-label details.

Each truncated coarse envelope is independently normalized. Its total intensity equals
the configured intact or loss intensity, subject to conversion to `f32`. Fine
envelopes remain unnormalized: source f32 abundances feed the configuration
probabilities, which are rounded to f32 before selection/output. Coverage sums
those widened f32 values in f64. The retained intensity therefore follows the
selected probability mass, not a forced sum of one. Both modes multiply the
configured intensity, relative loss factor when applicable, and isotope weight
in f64 before final f32 storage; there is no intermediate f32 loss-factor product.
Monoisotopic fragment names include charge markers, such as `b2++`. The source's
isotope intact names omit charge markers (`b2`); the separate charge array provides
the charge. Loss names include a canonical formula and charge markers, such as
`b2-H2O1++`. Precursor names include the adduct, such as `[M+2H-H2O]++`. Isotope
peaks repeat the same annotation without an isotope ordinal. With sorting disabled,
terminal series emit in b/y/a/c/x/z/z./z' order within each ascending charge,
followed by all internal groups, precursor groups and then immonium ions. Losses are deterministically ordered by canonical formula text.

## Validation and intentional differences

- Negative atom-count loss formulas are skipped in both monoisotopic and isotope
  modes. The source checks only its slower isotope path. This avoids generating
  chemically impossible a1(R)-CH2NO peaks in either mode.
- Configurations reject duplicate series, zero/reversed charge ranges,
  precursor charge below the maximum fragment charge, nonfinite/negative
  intensities, relative loss intensity above one, zero isotope bin count and
  nonfinite/out-of-range fine unexplained probabilities. Fine values must be in
  [0,1]; one produces empty envelopes while independent internal/immonium peaks
  remain enabled.
  `Some(0)` is invalid precursor charge; use `None` for automatic inference.
- Unsupported terminal-loss/isotope combinations return a clear error instead
  of silently ignoring the requested terminal losses. Empty formulas declared
  as modification neutral losses do not create duplicate intact peaks.
- Validation happens before mutation, including precursor charge validation.
  `append_to` works on a temporary result and replaces the target only on success.
- Existing annotation arrays are located by the exact names `IonNames` and
  `Charges`; duplicates are rejected. Existing unannotated peaks receive an
  empty name and charge zero when annotations are added. Either existing named
  array enables generation of both arrays. Empty unrelated arrays are preserved;
  populated unrelated arrays reject additions because their values for the new
  peaks are undefined. Existing array alignment and peak values must be valid.
- Per-peak intensities are `f32`; probability calculations are `f64`. Small
  rounding differences from the C++ float-based convolution are expected.

## Resource limits and deferred capabilities

Generation and combined append results are limited to 100,000 peaks
(`MAX_THEORETICAL_PEAKS`). Peptides may contain at most 4,096 residues
(`MAX_THEORETICAL_RESIDUES`), and `length² × enabled series` may not exceed
10,000,000 to bound repeated fragment construction. Peak-count preflight includes
potential terminal losses and the full requested isotope count, so it can
conservatively reject a request whose impossible losses would later be skipped.
Internal intact counts and up to seven immonium peaks are added after isotope
expansion. Internal losses observe the same peak limit as they are emitted;
only one start position (at most nine fragment lengths) is held at a time.
A separate 10,000,000 declaration-visit preflight includes repeated custom losses
before deduplication, with internal charge/type multiplicities. At most 100,000
distinct loss-template entries are retained cumulatively during a call, checked
before cloning their formulas. These conservative guards also count losses
whose formulas will later be rejected as impossible. Coarse envelopes also observe the isotope module's convolution limits. Fine
envelope sizes are checked as they are emitted; one shared 100,000,000-unit
fine work budget covers all fragments, losses and precursor charges in a
generation call. Per-envelope atom, frontier and storage caps apply as documented
in [fine isotope support](FINE_ISOTOPE_SUPPORT.md). Fragment charge
is a `u8` (1–255); precursor metadata charge is `u16` (1–65,535). Numerical
overflow and invalid intermediate formulas return errors.

Fine mode uses `TheoreticalIsotopeModel::Fine { unexplained_probability: 0.05 }`
to reproduce the source TSG default when its fine flag is selected. This is
95% requested coverage; the standalone fine generator defaults to 99%. IsoSpec
layer order is not reproduced; the separate native iterator exposes ordered
configurations and custom populations. The enum implements
`PartialEq`, without `Eq`, because this variant stores a floating probability.
Generation does not model fragmentation efficiencies or collision-energy effects;
configured intensities are deterministic weights.

## Internal fragments, immonium ions and presets

Internal groups emit b then a for each ascending charge, after all terminal
series. For each start, all intact lengths precede their neutral losses. The
source requires `start >= 1` and `start + 3 < peptide.len()`, so peptides of four
or fewer residues have none, and the last possible two-residue internal fragment
is omitted. Fragments exclude both peptide termini and stop at ten residues.
For `AGGGA`, the intact annotations are `GG`, `GGG`, `GG-CO`, `GGG-CO` at each
charge. Only loss annotations have charge markers. Residue modifications affect
mass and replace loss declarations; terminal modifications never enter these
fragments. Loss collection skips the first residue of each internal fragment,
matching the source loop. `add_terminal_losses` does not add internal losses.
Internal peaks remain monoisotopic even when terminal isotope envelopes are
selected. Formula-free numeric tags work when no terminal envelope or
formula-based precursor peaks are requested. As elsewhere in the native
generator, known impossible atom-count losses are skipped, and nonfinite or
negative final masses are errors; the source's fast internal path lacks these
checks.

Immonium peaks are emitted once per matching **unmodified** residue kind:

| Residue | m/z | Annotation |
| --- | ---: | --- |
| P | 70.0656 | `iP+` |
| C | 76.0221 | `iC+` |
| L | 86.09698 | `iL/I+` |
| H | 110.0718 | `iH+` |
| F | 120.0813 | `iF+` |
| Y | 136.0762 | `iY+` |
| W | 159.0922 | `iW+` |

They always have charge one and intensity one, independent of the requested
fragment charge range, series intensities and isotope settings. A modification
on one residue suppresses that slot's eligibility; a second unmodified residue
of the same kind restores it. Terminal modifications do not suppress eligibility.
Despite its `iL/I+` label, the source tests only L, so I alone emits no such peak.

`generate_for_activation(method, peptide, precursor_charge)` is an associated
function using fresh default settings with these series:

| Method | Series |
| --- | --- |
| CID | b/y |
| HCID, HCD | a/b/y |
| ECD, ETD | c/z./z' |
| ETciD, EThcD | all eight |

Other activation methods return `Unsupported`. Input charge zero follows the
source's default of two; inputs up to two generate charge-one fragments, and
larger inputs generate charges one and two. The source factory does not forward
the supplied charge to the spectrum generator: precursor metadata is inferred
as two or three from the fragment range. This convention is preserved even for
input charges one or five. Use `generate` when explicit precursor metadata is
required. The native unsigned parameter excludes negative charges; zero has no
logging side effect.

## Compact mass-only helper

`append_mass_spectrum(&mut values, &peptide, maximum_charge)` follows the source
`getPrefixAndSuffixIonsMZ` helper. It visits charges in descending order, appends
enabled ordinary a/b/c/x/y/z ladders including ordinal one **and the full length**,
then sorts the combined old and new `f32` values. Each ladder includes only its
corresponding terminal modification, including at full length. Radical z
variants and all other generator settings are ignored; repeated series act as
one enabled flag. A zero charge or empty peptide adds no peaks and still sorts
the existing values. This helper therefore has different endpoint and option
semantics from `generate`.

Residue masses are resolved once, and suffix ladders accumulate from the C end.
Protons, applicable terminal delta, ion conversion and residue masses are added
in the source's `f64` order before the final `f32` store. Known numeric masses
require no composition. Unresolved required masses, negative/nonfinite values
or `f32` overflow return an error before output changes. Limits cover 100,000
combined values, 4,096 residues while generating ions and 10,000,000 estimated
settings-scan, ladder and sort work units. Existing values are validated before
sorting even when no new ions are requested.

## Provenance and verification

Reviewed SHA-256 values at the pinned revision:

| Source | SHA-256 |
| --- | --- |
| `src/openms/source/CHEMISTRY/TheoreticalSpectrumGenerator.cpp` | `763d7e37041d3bed8f0b9fb122e43b5d6890fc699753e0f877ee3c50f2698325` |
| `src/openms/include/OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h` | `c0f8fc0213c7a9103f3e3ebeb305d498a79666afdb189b2b3b36c7d129cc81c9` |
| `src/tests/class_tests/openms/source/TheoreticalSpectrumGenerator_test.cpp` | `8fc12ab3cea2f804e1374fdc43368ad82588f465e01423753d03babb4c4a77e1` |

`tests/theoretical.rs` checks the original b/y and all-six-series mass tables,
loss names/counts, isotope masses/probabilities, precursor envelopes and
negative-loss regression counts. Independent tests cover z radical shifts,
modified-fragment and isotope mass anchors, modification loss replacement,
terminal losses, intensity scaling, source insertion order, aligned sorting,
empty/monomer inputs, invalid settings, limits and atomic append failures.

## Unresolved residues and numeric mass annotations

Sequences may preserve B/Z/X without defining a mass. Generation requires a
known monoisotopic mass for every residue and terminus, including the precursor
metadata, and errors for bare unresolved residues. An absolute numeric residue
tag such as `X[999]` supplies an internal mass; an unknown delta on an otherwise
known residue supplies its shift. These annotations support all monoisotopic ion
series without assigning an elemental formula. Named registry matches retain
their known formulas and isotope support.

Monoisotopic fragment losses follow the source's mass subtraction rules. An
anonymous modified residue supplies no neutral losses of its own and suppresses
its unmodified residue's losses; other retained residues still supply their
declared losses. If a fragment has a formula, impossible atom-count losses are
skipped as above. For a mass-only fragment, composition checks are unavailable;
calculated peak masses must still be finite and nonnegative. Terminal losses
remain an explicit option.

Coarse/fine isotope envelopes and the bundled `add_precursor_peaks` option require
the complete formula. The latter includes formula-based water/ammonia losses
even when fragment losses are disabled. Both return an error for mass-only
annotations, preserving an append target unchanged. This avoids the source's
silent omission of unknown tag masses in formula-based calculations. See
[sequence support](SEQUENCE_SUPPORT.md) and the independent cross-module checks
in [sequence_workflow.rs](../tests/sequence_workflow.rs).

```text
cargo test --offline --test theoretical
cargo clippy --offline --test theoretical -- -D warnings
```

The extensions have independent source-derived mass, ordering and boundary
checks in [theoretical extension review](THEORETICAL_EXTENSION_REFERENCE_REVIEW.md),
plus focused internal, preset and compact-helper tests.

No C++ build or live C++/Rust differential run was performed.

## Shared work in spectrum annotation

[Spectrum annotation](SPECTRUM_ANNOTATION_SUPPORT.md) invokes the same generator with shared residue, neutral-loss, fine-isotope and coarse-convolution counters across candidate hits. Coarse convolution now also shares its 50-million-product limit across all envelopes of an ordinary standalone theoretical spectrum; it is no longer restarted per envelope. Standalone calls otherwise start fresh counters, and all existing numerical/source conventions remain covered. Immutable anonymous mass-tag text is shared across owned fragment slices, avoiding repeated spelling copies.
