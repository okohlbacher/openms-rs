# Theoretical spectrum extension reference review

[theoretical_extension_reference.rs](../tests/theoretical_extension_reference.rs)
contains independent checks of internal fragments, immonium ions, activation
presets, and the compact floating-point mass helper. The
[provenance manifest](../tests/data/theoretical_extension_provenance.json)
records ten pinned source hashes, their links, the reference derivation, and
numerical tolerances. No C++ binary was built or executed. Tests require only
the Rust package and its bundled resources.

## Reference strength

The immonium values and CID/ECD values come directly from
`TheoreticalSpectrumGenerator_test.cpp`. The immonium comparison preserves the
seven exact decimal constants. CID's four-decimal values use an absolute
0.0001 Da tolerance. The ECD test checks the source's total count of seventeen
but only its eleven listed masses, using an absolute 0.000001 Da tolerance.

The source internal-fragment regression checks short-input termination and
increased output for `PEPTIDEK`; it does not supply a numerical mass table.
Our internal references are independently derived from the source loops and
the free-residue formulas and elemental masses in `ResidueDB.cpp`,
`ElementDB.cpp`, and `Constants.h`. They use an absolute 0.000000001 Da
tolerance for ordinary binary64 arithmetic differences. No call to the native
mass or spectrum generator supplies their expected values.

The compact-helper reference computes `AG` b1/b2/y1/y2 masses at charges one
and two, converts those independently calculated values to `f32`, and compares
their bits after sorting with existing output. This is an analytical reference,
not measured C++ output. The upstream `PrecursorPurity` test exercises the
helper through downstream self-consistency instead of a standalone literal
mass table.

## Source behaviors pinned by the independent suite

- Internal start positions obey `l >= 1 && l + 3 < n`, with retained lengths
  two through ten and neither peptide terminus. `AGGGA` therefore supplies only
  `GG` and `GGG`, each as b and a ions. The final possible two-residue interval
  is omitted. Inputs of length zero through four yield no internal ions.
- Fourteen alanines yield 124 internal peaks at one charge. Counts for retained
  lengths two through ten are 20, 20, 18, 16, 14, 12, 10, 8, and 6, including
  both series. The source's `PEPTIDEK` example yields 28 internal peaks.
- Internal labels retain raw sequence text, with `-CO` for a ions. Intact labels
  omit charge markers; loss labels contain them, while the aligned charge array
  records every charge. Repeated sequences are not merged.
- Each start emits intact lengths before loss lengths. Internal groups follow
  ordinary terminal-ion groups and precede precursor groups. Fixed immonium
  ions are last when sorting is disabled.
- Loss collection omits the first retained residue. `AKGGA` has no internal
  ammonia losses, while `AGKGA` supplies ammonia loss for both retained lengths.
  Intensity and neutral-loss mass subtraction are independently checked.
- Internal ions remain monoisotopic when the isotope option is selected, and
  do not acquire terminal-loss additions. Peptide terminal modifications do not
  change them. An absolute `X[999]` tag contributes its retained internal mass
  without establishing an empirical formula.
- The seven immonium constants remain charge one and unit intensity even when
  the requested fragment range begins at charge two. The separate preset tests
  also cover unmodified-residue eligibility, the source L-only I/L branch, and
  the activation factory's inferred precursor metadata charge.
- The compact helper includes first and full-length ordinary fragments,
  treats repeated series as boolean flags, and ignores radical z variants and
  the other generation settings. Empty input, charge zero, and no enabled
  ordinary series still sort valid existing values and do not require unknown
  residue chemistry.

## Implementation review

The internal implementation preserves proton-first accumulation, the a-ion CO
offset, the source loop bounds, and intact-before-loss ordering. The compact
helper preserves source accumulation through the final `f32` store: protons,
the retained terminal delta, ion conversion, and successive internal residue
masses. Suffix masses are accumulated from the C terminus directly. Activation
presets use the existing activation enum and default generator settings.

Native validation retains the documented restrictions on nonfinite values,
negative masses, impossible formula losses, and unsupported terminal-loss /
isotope combinations. These restrictions are separate from the source's
observable internal-fragment conventions. Materialized fine isotope patterns
were added subsequently; see [fine isotope support](FINE_ISOTOPE_SUPPORT.md).
Internal fragments retain their independent monoisotopic path. IsoSpec
layered traversal remains outside the implemented scope; native ordered streaming
and custom isotope populations are documented with fine support.

The review identified that repeated or distinct custom neutral-loss declarations
could be rescanned and cloned before the final peak-count check. The integrated
fix checks weighted declaration visits before scanning them and charges a
cumulative unique-loss-template budget before cloning formulas. Both ordinary
and internal fragments use these limits; repeated declarations cannot evade the
work bound merely by deduplicating to a small output. Independent inspection
confirmed the ordinary ordinal ranges and internal first-residue exclusion in
the accounting.

All nine reference tests pass on the current compiler with all features and
Rust 1.85 without default features. The final validation report records the
broader integrated checks. No additional scientific discrepancy was found in
the new internal, activation, immonium, or compact-helper implementations.
