# Chemistry support

The native `openms::chemistry` module implements molecular compositions and adduct mass conversion,
modified peptide chemistry, charge/pI and sequence properties, coarse and fine isotope patterns, configurable theoretical peptide spectra,
full/semi/nonspecific digestion with 33 enzymes, source-compatible decoy generation,
residue-mass tag extraction from spectra, RNA digestion/modification enumeration,
RNA spectrum generation, and RNA records/sequences with a full
pinned registry. Scientific calculations require no C++ library or runtime
resource files. Optional MODOMICS JSON import uses serde_json; embedded RNA
data and TSV import work without it. RNA digestion also registers owned oligos,
parent positions and processing history in the separate
[identification sequence/provenance graph](IDENTIFICATION_GRAPH_SUPPORT.md).

## Source and coverage

The historical chemistry reference is [OpenMS4-core at
7c029e8cdba6abab503708ecdd56f6ab55e38ce4](https://github.com/okohlbacher/OpenMS4-core/tree/7c029e8cdba6abab503708ecdd56f6ab55e38ce4).
The [current SDK update](CORE_SDK_UPDATE.md) verifies unchanged scientific sources
at `6bfc0e4`; historical fixture pins are retained.
The Rust implementation retains the upstream BSD-3-Clause attribution; the private decoy random helper also retains Boost Software License 1.0 terms. Embedded
UniMod data has its own Design Science License; see the [modification support](MODIFICATION_SUPPORT.md).

| Rust API | Implemented behavior | Reference |
| --- | --- | --- |
| `element_table`, `element`, `Element`, `Isotope` | All 84 declared natural-element tables, individual isotope masses and abundances, most-abundant-isotope and average masses | `src/openms/source/CHEMISTRY/ElementDB.cpp` |
| `EmpiricalFormula` | Parsing, isotope labels, signed atom counts, protonation charge, checked addition/subtraction/scaling, containment, mono/average mass and m/z | `src/openms/source/CHEMISTRY/EmpiricalFormula.cpp` |
| `AASequence` | Twenty canonical residues plus U/O/J and unresolved B/Z/X; named and numeric annotations, checked formula/masses/m/z and subsequences; [details](SEQUENCE_SUPPORT.md) | `src/openms/source/CHEMISTRY/ResidueDB.cpp`, `AASequence.cpp` |
| `IsoelectricPoint`, `HydrophobicityProfile`, `AAIndex` | Four pKa scales, seven hydrophobicity scales, profiles/moments, ten amino-acid indices and checked gas basicity; [details](PEPTIDE_PROPERTIES_SUPPORT.md) | `IsoelectricPoint.cpp`, `HydrophobicityProfile.cpp`, `Residue.cpp`, `AAIndex.h` |
| `fragment_ions` | b/y ions at every internal bond, charges 1 through the requested maximum, residue ordinals | `src/openms/include/OpenMS/CHEMISTRY/Residue.h` |
| `TheoreticalSpectrumGenerator` | a/b/c/x/y/z and z-variant spectra, intensities, neutral losses, precursors, coarse/fine isotope envelopes, aligned names/charges, atomic append; [details](THEORETICAL_SPECTRA.md) | `src/openms/source/CHEMISTRY/TheoreticalSpectrumGenerator.cpp` |
| `Ribonucleotide`, `RibonucleotideDB`, `NASequence` | Full nucleoside records and registry/providers; nucleic-acid parsing, owned identity, terminal/sulfur slicing, all finite source fragment formulas and ion masses; [details](RNA_SUPPORT.md) | `Ribonucleotide.cpp`, `RibonucleotideDB.cpp`, `NASequence.cpp`, TSV/JSON providers |
| `DigestionEnzymeRNA`, `RNaseDB`, `RNaseDigestion` | Fourteen pinned enzymes, caller-owned records/registry, modified-code matching, ordered length/missed-cleavage products and terminal gains, atomic graph registration with parent positions/history; [details](RNASE_SUPPORT.md) | `RNaseDigestion.cpp`, `DigestionEnzymeRNA.cpp`, `RNaseDB.cpp` |
| `ModifiedNASequenceGenerator` | Fixed replacement and bounded variable enumeration with complete record ownership, source fast-path/site order and atomic updates; [details](RNA_MODIFICATION_SUPPORT.md) | `ModifiedNASequenceGenerator.cpp` |
| `NucleicAcidSpectrumGenerator` | Both generation operations, nine ion series, signed charges, source precursor/sulfur branches and aligned annotations; [details](RNA_SPECTRUM_SUPPORT.md) | `NucleicAcidSpectrumGenerator.cpp` |
| `Tagger`, `TaggerOptions` | Complete measured m/z/spectrum tag extraction, length/charge/tolerance settings, fixed/variable masses, I/L branching and atomic append; [details](TAGGER_SUPPORT.md) | `Tagger.cpp`, `Residue.cpp` |
| `AdductInfo` | Molecular adduct parsing, charge/electron and n-mer conversions, mono/average mass shifts, compatibility and record identity; [details](ADDUCT_SUPPORT.md) | `AdductInfo.cpp` |
| `DecoyGenerator` | Protein and peptide reversal, stateful seeded shuffling and independent variants, source cache/digestion conventions and checked output/work; [details](DECOY_GENERATION_SUPPORT.md) | `DecoyGenerator.cpp`, `MathFunctions.h` |
| `SpectrumAnnotator`, `ion_naming` | Spectrum arrays, hit peak annotations, source match statistics and four charge/name helpers with shared work limits; [details](SPECTRUM_ANNOTATION_SUPPORT.md) | `SpectrumAnnotator.cpp`, `IonNaming.h` |
| `ModificationsDB` | Pinned UniMod/custom registry, specificity-aware name/accession and mass lookup, neutral losses | `ModificationsDB.cpp`, `UnimodXMLHandler.cpp` |
| `CoarseIsotopePatternGenerator`, `IsotopeDistribution` | Coarse distributions, convolution, enrichment overrides, averagine estimates and conditional fragment envelopes; [details](ISOTOPE_SUPPORT.md) | `CHEMISTRY/ISOTOPEDISTRIBUTION` |
| `ProteaseDB`, `Protease`, `ProteaseDigestion` | All 33 pinned enzymes and metadata; full/semi/nonspecific digestion, ranges, counts, validity and missed-cleavage checks; [details](DIGESTION_SUPPORT.md) | `ProteaseDigestion.cpp`, `EnzymaticDigestion.cpp`, `share/OpenMS/CHEMISTRY/Enzymes.xml` |

`PROTON_MASS_U`, `ELECTRON_MASS_U`, and `C13C12_MASSDIFF_U` use the upstream
`Constants.h` values. Formula masses use the declared ElementDB table, whose
rounded carbon isotope masses differ slightly from the dedicated spacing constant.

The reviewed local source files had these SHA-256 hashes:

| Source file | SHA-256 |
| --- | --- |
| `src/openms/source/CHEMISTRY/ElementDB.cpp` | `16bd77ddf39f9d374a06e70ac0e690358ea855184175a9a6ecedea64ae9b3e58` |
| `src/openms/source/CHEMISTRY/ResidueDB.cpp` | `e3fa85271263d2396f4ad647cc0fe634572e4a93437a0de2ef69ce8c8ef3675f` |
| `share/OpenMS/CHEMISTRY/Enzymes.xml` | `2f160f3ec32db6257cb4eee43fd594b48cb36398a17b7af21bbba76bc8b16995` |

## Native API decisions and differences

- Elements have private fields and no public constructor; use `element` or
  `element_table` to obtain validated immutable entries, and read their properties
  through `name()`, `symbol()`, `atomic_number()`, and `isotopes()`. This preserves
  the nonempty isotope-table invariant even when an element is copied.
- Formula charge means adding or removing protons, as in OpenMS. It does not
  represent electron-only ionization. `C6H12O6+2` is glucose plus two protons.
- Formula notation accepts `(13)C`, `D`/`T` aliases, negative counts, and signed
  terminal charges. `H-2` is a count of minus two hydrogen atoms; `H1-2` is one
  hydrogen atom with charge minus two. `with_charge` avoids this ambiguity.
  Display always includes atom counts for unambiguous round trips. Outer
  whitespace is trimmed; interior whitespace and parenthesized group multipliers
  are rejected. Bare trailing `-` is accepted as charge minus one.
- Arithmetic returns errors for `i32` count or charge overflow and never wraps.
  Zero-count entries disappear. Natural-element counts and explicit isotope
  counts are separate. Formula mass permits signed compositions for neutral
  loss arithmetic; m/z rejects a negative calculated mass.
- Negative-charge m/z divides by the absolute charge, returning physical positive
  m/z. C++ `AASequence::getMZ` divides by signed charge. Zero charge is an error.
  Empty sequences have mass zero and no fragment ions; their m/z is an error.
- Iridium uses the **declared iridium isotope table** (191 and 193). The pinned
  C++ source accidentally passes the rhenium tables to `buildElement_` for
  iridium (`ElementDB.cpp:512` at `bc9cc12`). This is an intentional correction,
  covered by a regression test. Iridium masses and isotope patterns therefore
  differ from the executed C++: every SDK run gives `Os3Ir3` a rhenium-based
  lightest-isotope weight of 1106.716344 Da, where iridium gives 1124.739252 Da
  (CPP-249).
- Peptides preserve uppercase B/Z/X without inventing a mass or composition.
  `formula`, `mono_mass` and `average_mass` return `Result`; unresolved chemistry
  returns `Unsupported`. Numeric annotations resolve to registry modifications
  when the source precision rule matches; otherwise an owned mass tag retains
  its mass with no invented atoms or average mass. J shares the composition of
  I/L. See [sequence support](SEQUENCE_SUPPORT.md) for grammar and available chemistry.
- `AASequence::fragment_ions` remains a lightweight b/y mass-list API, including
  every internal bond. `TheoreticalSpectrumGenerator` produces annotated
  `MSSpectrum` values with configurable series and intensities. Its source
  defaults use b/y ions and omit the first prefix ion; losses, envelopes,
  precursor peaks and annotations are opt-in. Configured intensities are
  deterministic weights, not a fragmentation-efficiency prediction.
- The generator carries retained residue and terminal modifications into each
  fragment. Modified residues use declared modification losses instead of the
  unmodified residue's losses. Impossible negative-count loss formulas are
  skipped, and unsupported terminal-loss/isotope combinations return errors.
  Coarse and fine fragment envelopes retain the source's neutral-hydrogen-adduct
  convention, including its electron-mass offset from proton-based m/z.
  See [theoretical spectra](THEORETICAL_SPECTRA.md) for formulas, charge limits,
  probability normalization, annotation conventions and resource bounds.
- Internal b/a fragments preserve source interval/loss rules; immonium peaks
  retain fixed source masses and unmodified-residue eligibility. Activation
  presets and the sorted `f32` mass-only helper have their own documented charge
  and endpoint conventions. [Fine isotope patterns](FINE_ISOTOPE_SUPPORT.md) retain distinct configurations, source f32 abundance/output rounding and natural-H charge handling; the owning raw iterator also supports custom binary64 populations and checked threshold selection.
- The immutable registry preserves all 33 pinned enzyme definitions, synonyms,
  terminal-gain metadata and search-engine IDs. Native predicates follow the exact
  source expressions. Trypsin excludes cleavage before P, including WKP/MRP.
  Unmodified string APIs accept uppercase A–Z for cleavage analysis; ambiguous
  B/Z/X can also be retained in `AASequence`; bare ambiguous residues have no
  calculable mass. Absolute residue tags can supply a monoisotopic mass.
- Digestion supports full, semi and unrestricted specificity while preserving
  retained residue/terminal modifications. Fully specific products use increasing
  missed-cleavage count then position; semi variants follow in source order;
  unrestricted products use start then length. Counts respect the same length and
  specificity settings as generation. Inclusive length bounds use `None` for an
  unbounded maximum. Range APIs avoid copying peptide sequences.
- Product validity offers explicit N-terminal M/MX and D|P allowances. It safely
  rejects invalid ranges and uses the full protein context for missed-cleavage
  counting. The [digestion support](DIGESTION_SUPPORT.md) documents unrestricted
  generation versus validity semantics, source differences and resource limits.

## Deliberately outside this increment

IsoSpec layered traversal and backend performance hints; RNA enzyme XML import; ProForma;
arbitrary enzyme regular expressions and protein-enzyme runtime registry overrides;
property prediction beyond the documented pI, hydrophobicity and AAIndex utilities;
and mutable element/residue databases are not implemented.
The separate identification graph covers observations, compounds, adducts,
observation matches, parent/match groups, referential cleanup and the bounded
legacy sequence/evidence conversion bridge. Graph persistence and full conversion
remain outstanding.

## Verification

`tests/chemistry.rs` includes golden cases from the upstream EmpiricalFormula,
AASequence, and ProteaseDigestion class tests. It also checks every ported element
and labeled isotope, signed-count and charge round trips, malformed input and
overflow, empty and single-residue sequences, fragment mass conservation,
proline boundaries, missed cleavages, and positional length filtering.
`tests/digestion.rs` checks all 33 compiled enzyme rules against 602,316
independent regex contexts, source digestion/validity examples, exhaustive small
range enumeration, modified products and resource boundaries.
`tests/modifications.rs`, `tests/isotopes.rs`, and `tests/theoretical.rs` extend
coverage to the embedded modification registry, isotope calculations, and
theoretical spectra. The latter includes source ion-mass tables, loss names and
counts, isotope probabilities and precursor envelopes, plus modified-fragment,
annotation-alignment, resource-limit and atomic-append checks. Each support
document records the relevant source provenance and remaining differences.

The historical `DFPIANGER` monoisotopic golden (1017.487958568 Da) differs from
its composition evaluated with the pinned element table (1017.4879641373 Da)
by 5.6 microdaltons. Tests retain the historical check with an appropriate
absolute tolerance and independently check the pinned formula mass tightly.

Run the chemistry tests without optional format dependencies:

```text
cargo test --offline --no-default-features --test chemistry --test modifications --test isotopes --test theoretical --test digestion
```

The chemistry implementation has been tested in Rust. No C++ build or live
C++/Rust differential execution was performed.
