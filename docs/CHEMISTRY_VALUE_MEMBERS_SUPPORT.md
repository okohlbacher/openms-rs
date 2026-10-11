# Chemistry value members: the remainder of eight CHEMISTRY headers

Lane A5a (wave 4) closed the members the ledger still listed as missing for
`CHEMISTRY/EmpiricalFormula.h`, `Element.h`, `ElementDB.h`,
`ISOTOPEDISTRIBUTION/IsotopeDistribution.h`,
`ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h`, `Ribonucleotide.h`,
`EnzymaticDigestion.h` and `DigestionEnzymeDB.h` at core `bc9cc12`. The earlier
parts of each header keep their own support documents
(`CHEMISTRY_SUPPORT.md`, `ISOTOPE_SOURCE_PRECISION_SUPPORT.md`,
`RNA_SUPPORT.md`, `DIGESTION_SUPPORT.md`); this document covers what was added.

| Header | Rust |
|---|---|
| `EmpiricalFormula.h` | `src/chemistry/empirical_formula.rs` (with the type in `src/chemistry/mod.rs`) |
| `Element.h`, `ElementDB.h` | `src/chemistry/element_db.rs` (with the table in `src/chemistry/mod.rs`, `elements.rs`) |
| `IsotopeDistribution.h` | `src/chemistry/isotopes/distribution_members.rs` |
| `CoarseIsotopePatternGenerator.h` | `src/chemistry/isotopes/generator_members.rs` |
| `Ribonucleotide.h` | `src/chemistry/ribonucleotide.rs` |
| `EnzymaticDigestion.h` | `src/chemistry/enzymatic_digestion.rs` |
| `DigestionEnzymeDB.h` | `src/chemistry/digestion.rs` (`ProteaseDB`), `src/chemistry/rnase.rs` (`RNaseDB`) |

## API mapping

### `EmpiricalFormula.h`

| Source | Rust |
|---|---|
| `EmpiricalFormula(SignedSize, const Element*, SignedSize charge = 0)` | `EmpiricalFormula::from_element` (an isotope element such as `(13)C` gives the labelled isotope) |
| `static hydrogen(int)`, `static water(int)` | `hydrogen`, `water` (checked: `Result`) |
| `getLightestIsotopeWeight` | `lightest_isotope_mass` |
| `calculateTheoreticalIsotopesNumber` | `theoretical_isotopes_number` |
| `addChargeAdduct(Int = 1, const std::string& = "H")` | `add_charge_adduct(count, adduct)` |
| `toString` | `to_source_string` (`Display` stays the round-trip form) |
| `toMap` | `to_map` |
| `isCharged`, `setCharge`, `hasElement`, `getNumberOf` | `is_charged`, `set_charge`, `has_element`, `number_of` |
| `begin`/`end` and the iterator typedefs | `entries()` yielding `FormulaEntry` |
| `operator<` | `Ord` |
| `friend operator<<` | `source_stream_text` |
| `std::hash<EmpiricalFormula>` | `Hash`; `source_hash64` is the source value |

### `Element.h`, `ElementDB.h`

| Source | Rust |
|---|---|
| `ElementDB::getInstance` | the free functions of `chemistry::element_db` |
| `getNames`, `getSymbols`, `getAtomicNumbers` | `names()`, `symbols()`, `atomic_numbers()` (`&'static HashMap`) |
| `getElement(const std::string&)`, `getElement(unsigned)` | `get_element`, `get_element_by_atomic_number` |
| `hasElement(const std::string&)`, `hasElement(unsigned)` | `has_element`, `has_atomic_number` |
| the isotope elements `(13)C`/`(13)Carbon` and the `D`/`T` symbols | built once by `element_db`; `Element::is_isotope` |
| `Element::getIsotopeDistribution` | `Element::isotope_distribution` |
| `Element::operator<` | `Ord` (atomic number, mono weight, symbol, name, average weight, isotopes) |
| `Element::operator<<` | `Element::source_stream_text` |
| `std::hash<Element>` | `Hash`; `Element::source_hash64` |
| `Element` constructors and the six setters | **not ported**: elements are database records (see below) |

### `IsotopeDistribution.h`, `CoarseIsotopePatternGenerator.h`

| Source | Rust |
|---|---|
| `enum Sorted` | `isotopes::Sorted` |
| `merge(double resolution, double min_prob)` | `IsotopeDistribution::merge` |
| `operator<` | `Ord` |
| `std::hash<IsotopeDistribution>` | `Hash`; `source_hash64` |
| mutable `operator[]`, mutable iterators, protected `transform_` | `set_peak`, `try_for_each_peak_mut` (validated) |
| `setMaxIsotope`, `setRoundMasses` | `set_max_peaks`, `set_mass_mode` |
| `getIsotopeOverrides` | `isotope_overrides()` as (symbol, distribution) pairs |
| `estimateForFragmentFromPeptideWeightAndS` | `estimate_fragment_from_peptide_weight_and_sulfur` |

### `Ribonucleotide.h`

`operator<<` is `source_stream_text` and `std::hash<Ribonucleotide>` is
`source_hash64`. `Display` and `Hash` already existed; `Display` writes the
formula's round-trip form.

### `EnzymaticDigestion.h`, `DigestionEnzymeDB.h`

| Source | Rust |
|---|---|
| `enum Specificity` (0, 1, 2, 3, 8, 9) and `SIZE_OF_SPECIFICITY` | `digestion::DigestionSpecificity` (also named `enzymatic_digestion::Specificity`), `SIZE_OF_SPECIFICITY` |
| `NamesOfSpecificity` | `NAMES_OF_SPECIFICITY` |
| `NoCleavage`, `UnspecificCleavage` | `NO_CLEAVAGE`, `UNSPECIFIC_CLEAVAGE` |
| `getSpecificityByName` | `DigestionSpecificity::by_name` |
| `getSpecificity`/`setSpecificity` with the source enum | the `specificity` field of `ProteaseDigestion` |
| `filterByMissedCleavages` | `ProteaseDigestion::filter_by_missed_cleavages` (the source's arithmetic, `count_internal_cleavage_sites_source`) |
| `DigestionEnzymeDB::hasEnzyme(const T*)` | `ProteaseDB::has_enzyme_record`, `RNaseDB::has_enzyme_record` |
| `DigestionEnzymeDB::getEnzymeByRegEx` | `ProteaseDB::enzyme_by_regex` (single, last registered); `enzymes_by_regex` returns all |
| the lowercase-name index of `addEnzyme_` | `ProteaseDB::get_enzyme`/`has_enzyme` (fixed; `RNaseDB` already had it) |

## Preserved source conventions

- An isotope element is its own entry everywhere: `(13)C` is counted, hashed
  and ordered apart from `C`, and `to_map`/`toString`/`operator<<` key by the
  source symbol text in byte order, so `(13)C` precedes `C` and `(2)H` precedes
  `C`.
- `operator<<` of a formula omits counts of one and every negative count, and
  writes a charge of `-2` as `--2` (`os << "-" << charge_`). `toString` writes
  every count, including `1` and negatives, and never the charge.
- `addChargeAdduct` ignores any charge written in the adduct and sets the
  charge to zero unconditionally.
- `merge` keeps the source's binning: `ceil((max - min) / resolution)` bins,
  index `round((m - min) / resolution)`, a peak landing on bin `n` dropped (so
  the heaviest peak is lost when the range is an exact multiple of the
  resolution), bins nothing landed in kept as `(0, 0)` when `min_prob` is zero,
  and a single remaining mass giving an empty distribution.
- `getEnzymeByRegEx` returns the enzyme registered last among those sharing an
  expression (`Clostripain/P`, `cyanogen-bromide`, `Glu-C+P`).
- Source hashes (`source_hash64`) use the `CONCEPT/HashUtils.h` recurrences
  over little-endian bytes and narrow `Peak1D` intensities to `float`, so they
  equal the SDK's `std::hash` values on x86_64.

## Native differences

- Elements are immutable database records: the `Element` constructors and
  setters are not ported, and no custom element can be created. This is the
  one gap left in `Element.h`.
- Iridium uses its declared table, not rhenium's (CPP-249); its isotope
  elements are `(191)Ir`/`(193)Ir` where the SDK has `(185)Ir`/`(187)Ir`, and
  its stream text and hash differ.
- `theoretical_isotopes_number` refuses a negative count on the binomial
  branch; the SDK converts it to unsigned and returns about `2.8e19` for
  `H-2O`. Boost's domain and overflow errors are errors.
- `water(0)` and `from_element(0, ..)` store no zero entries (the SDK stores
  `H0O0`, which is not empty); `water` errors where `2 * n` overflows.
- The multi-element `getLightestIsotopeWeight` sum follows this type's entry
  order; the SDK sums in element address order, so the last bit can differ.
- `IsotopeDistribution::merge` adds weights in `f64` (the SDK in `float`) and
  refuses an empty trimmed distribution, which the SDK reads past its end.
- A digestion holds `SPEC_UNKNOWN`, `SPEC_NOCTERM` and `SPEC_NONTERM` as the
  SDK does, and `ProteaseDigestion::digest` refuses them as the SDK does
  (`ProteaseDigestion.cpp:93-98`). The group-A `set_source_specificity` was
  removed at the wave-4 integration, when group B's `DigestionSpecificity`
  gained the three values.
- `estimate_fragment_from_peptide_weight_and_sulfur` refuses a fragment sulfur
  count above the precursor's, which the SDK subtracts as unsigned integers.

## Checked boundaries and evidence

Tier 1, executed Release oracles (ibminode06, Linux x86_64 Release install
`openms4-release-bc9cc12-c19e494-174b576`, `x86_64-conda-linux-gnu-g++`
14.4.0, libOpenMS core flags; each run twice, identical):

- `../oracle/A5a-chem` (`driver.cpp` sha256
  `939f86c72a993ab6cb4221a3a550d0ed94ceeae55fc7ef97ea5bd94b4ee50664`,
  `results/records.tsv` sha256
  `8de8db89310d5c776c0b908936d04e8394a338fb817f49f2e4bea50736f59263`):
  `tests/chemistry_members_oracle.rs` (element database sizes and lookups,
  element stream text, weights and hashes, formula text, hashes, isotope counts
  and lightest weights, ribonucleotide text and hashes, `merge`, distribution
  hashes and order, and the sulfur fragment estimate bit for bit under
  `ProbabilityPrecision::SourceSingle`), plus the lowercase-name rows of
  `tests/enzymatic_digestion_members.rs` and `tests/digestion.rs`.
- `../oracle/A5a-enzymes` (`driver.cpp` sha256
  `aabf6fcafac8babfba0c529589d3355183923a49c26df5d3a9899cd543237627`,
  `results/records.tsv` sha256
  `be7650be3d6cb9aa613250c2308ed01b4b41a2b954edba3d16c44c7a73659487`):
  `getEnzymeByRegEx` for all 30 protease expressions, the RNase name index.

Tier 2/3, class tests at `bc9cc12`: `EmpiricalFormula_test.cpp`,
`Element_test.cpp`, `ElementDB_test.cpp`, `IsotopeDistribution_test.cpp`,
`CoarseIsotopeDistribution_test.cpp` (lines 84-104, 355-397),
`Ribonucleotide_test.cpp`, `ProteaseDB_test.cpp`, in
`tests/empirical_formula_members.rs`, `tests/element_db.rs`,
`tests/isotope_distribution_members.rs`, `tests/coarse_generator_members.rs`,
`tests/ribonucleotide_source_forms.rs` and
`tests/enzymatic_digestion_members.rs`.

Bounded work: the element database is built once from the compiled table (a
few hundred entries); `merge` allocates at most as many bins as the trimmed
distribution has peaks; every other member is constant or linear in its
input, which is already bounded by its own type's limits.
