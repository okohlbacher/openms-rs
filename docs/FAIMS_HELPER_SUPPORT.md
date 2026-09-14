# FAIMS helper support

Native coverage of `IONMOBILITY/FAIMSHelper.h` and
`IONMOBILITY/FAIMSHelper.cpp` at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Artifact | Path |
| --- | --- |
| Implementation | `src/kernel/faims_helper.rs` |
| Tests | `tests/faims_helper.rs` |
| Fixtures | `tests/data/faims_helper/IM_FAIMS_test.mzML` (the pinned class-test file, byte-identical); `tests/data/faims_helper/oracle_cases.tsv` |
| Manifest | `tests/data/faims_helper_provenance.json` |
| Oracle driver | `../oracle/pte-faims-helper/` (outside this repository) |

The module lives in `kernel`, next to `spectrum_mobility` and
`experiment_mobility`, rather than in a new `ionmobility` top-level module. It
uses only module edges `kernel` already has (`concept`, `identification`,
`metadata`), so `tools/check_module_cycles.py` records no new edge.

Consumers in the pinned core: FileInfo calls `getCompensationVoltages`
unconditionally on every peak file, twice (`FileInfo.cpp:1673`, `:1742`);
`IMDataConverter::splitByFAIMSCV` groups spectra by it (`IMDataConverter.cpp:31`),
which FeatureFinderCentroided uses; `Biosaur2Algorithm.cpp:266` and
`SpectrumMetaDataLookup.cpp:240` call it too; and
`FeatureFinderIdentificationAlgorithm.cpp:561` calls `filterPeptidesByFAIMSCV`.

## API mapping

| C++ member | Rust |
| --- | --- |
| `class FAIMSHelper` | `kernel::faims_helper::FaimsHelper`, a zero-sized unit struct: the source class holds no state |
| `FAIMSHelper()` (implicit), `virtual ~FAIMSHelper()` | derived `Default`; no drop glue. The virtual destructor has no Rust counterpart: nothing derives from the class |
| `static std::set<double> getCompensationVoltages(const PeakMap& exp)` | `FaimsHelper::get_compensation_voltages(&MSExperiment) -> Result<CompensationVoltages>` |
| its `std::set<double>` return value | `CompensationVoltages::voltages: BTreeSet<CompensationVoltage>`; `CompensationVoltages::values()` iterates the voltages as `f64` in set order |
| its `OPENMS_LOG_WARN` record (`FAIMSHelper.cpp:52`) | `CompensationVoltages::warnings`, holding `FaimsHelper::MISSING_VOLTAGE_WARNING` verbatim |
| `static PeptideIdentificationList filterPeptidesByFAIMSCV(const PeptideIdentificationList& peptides, double target_cv, double cv_tolerance = 0.01)` | `FaimsHelper::filter_peptides_by_faims_cv(&[PeptideIdentification], f64, f64) -> Result<Vec<PeptideIdentification>>` |
| the default argument `cv_tolerance = 0.01` | `FaimsHelper::DEFAULT_CV_TOLERANCE`, passed explicitly |
| `PeptideIdentificationList` | a borrowed slice in, a `Vec` out |
| used: `Constants::UserParam::FAIMS_CV` | `concept::constants::user_param::FAIMS_CV` (ported earlier) |
| used: `IMTypes::DRIFTTIME_NOT_SET` | `metadata::ImTypes::DRIFTTIME_NOT_SET` (ported earlier) |
| used: `DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE` | `metadata::DriftTimeUnit::FaimsCompensationVoltage` (ported earlier) |
| native only | `CompensationVoltage` (`new`, `volts`, `TryFrom<f64>`, `From<CompensationVoltage> for f64`, total `Ord`, `Eq`, `Hash`); `FaimsHelper::MAX_SPECTRA`; `FaimsHelper::MAX_PEPTIDE_IDENTIFICATIONS` |

## The voltage key and its ordering

`std::set<double>` orders by `operator<` and treats two values as one element
when neither is less than the other. `CompensationVoltage` reproduces that for
every value it can hold:

- ascending numeric order, `-inf` first and `+inf` last, so FAIMS voltages come
  out `-65 < -55 < -45` (oracle `ascending_numeric_order`, `infinities`);
- `-0.0` and `+0.0` are one element, and the set keeps the sign of the zero it
  received first, as `std::set::insert` does (oracle `negative_zero_first`,
  `positive_zero_first`);
- no tolerance: adjacent doubles stay distinct (oracle
  `adjacent_doubles_stay_distinct`).

It is implemented as `f64::total_cmp` on the value with `-0.0` folded into
`+0.0`, with `Eq` and `Hash` defined on the same folded key. NaN is refused at
construction, so the order is a genuine total order.

The obvious alternative, `f64::to_bits`, is not this order. IEEE 754 bit
patterns place every negative value after every positive value and order
negative values by magnitude: `(-45.0f64).to_bits() < (-65.0f64).to_bits()`.
FAIMS data is mostly negative, so a bit-pattern key would reverse the list
FileInfo prints and the order in which `splitByFAIMSCV` creates its groups. The
test `voltages_order_numerically_and_not_by_bit_pattern` pins this.

## Preserved source conventions

- **Every spectrum is scanned.** A FAIMS spectrum anywhere in the experiment
  counts, not only the first (class-test section 4).
- **Only spectra whose unit is `FaimsCompensationVoltage` contribute**, whatever
  the drift time of the others, including the sentinel and NaN (oracle
  `sentinel_on_non_faims`, `unit_none_with_value`, `nan_on_non_faims`).
- **The sentinel `-1.0` is removed after collection**, by set equality, so
  exactly `-1.0` is removed, and a warning is produced exactly when something
  was removed (oracle `class_section_beyond_first_and_sentinel`,
  `only_sentinel`). The text is the source's.
- **One pass.** The source returns early for an experiment without spectra and
  pre-scans for any FAIMS spectrum before collecting. Collecting only from
  FAIMS spectra yields the empty set in both cases anyway; the results are
  identical.
- **Filter rule.** An identification is kept when its metadata has no
  `FAIMS_CV` key, or when `|cv - target_cv| < cv_tolerance`. The comparison is
  strict and has no epsilon (oracle `tolerance_half`: `-44.5` against `-45` at
  tolerance 0.5 is dropped). In binary arithmetic `-44.99` differs from `-45` by
  0.00999999999999801, so it is kept at the default tolerance (oracle
  `default_tolerance`). Only the exact key counts: a `FAIMS` key leaves an
  identification unannotated, and it is kept.
- **Infinite targets are filtered, not refused.** `FAIMSHelper.cpp` has no
  check on `target_cv`. Against `+inf` or `-inf`, `|cv - target_cv|` is
  infinite for every annotation and never less than the tolerance, an infinite
  tolerance included, so only the unannotated identifications are kept (oracle
  `target_positive_infinity`, `target_negative_infinity` and their
  `_tolerance_infinite` variants). The port returns the same identifiers. This
  keeps the two functions consistent: `get_compensation_voltages` returns
  infinite voltages (oracle `infinities`), and every voltage it returns is a
  valid target. The source keeps no annotated identification even when the
  annotation is the same infinity, because `inf - inf` is NaN (oracle
  `faims_filter_infinite_annotation`, where `+inf`, `-inf` and `-45`
  annotations are all dropped at an infinite tolerance). `MetaValue` refuses
  non-finite floats, so the port cannot build such an annotation. On the
  representable identifications it keeps what C++ keeps.
- **Integer annotations** convert to `f64` as `DataValue::operator double()`
  converts `INT_VALUE` (oracle `p4_int_-45`, kept).
- **Input order is preserved**, and kept identifications are copies; the input
  is not modified.

## Native differences

- **NaN voltages are refused.** A FAIMS spectrum with a NaN drift time makes
  `get_compensation_voltages` return `Error::InvalidValue` naming the spectrum
  index. The source inserts the NaN into its `std::set<double>`, which breaks
  the strict weak ordering the container requires. The executed oracle shows
  the result depends on insertion order. As the first element the NaN becomes
  the tree root, every later voltage compares equal to it and is dropped, and
  `erase(DRIFTTIME_NOT_SET)` then erases the NaN itself and logs the
  missing-voltage warning. The set comes back **empty with a spurious warning**
  (`nan_first`). Anywhere later, the NaN alone is dropped (`nan_middle`,
  `nan_last`). This is reported as a C++ issue candidate.
- **Warnings are returned, not logged.** This kernel module is not wired to
  `LogStream`; the caller decides where `CompensationVoltages::warnings` goes.
- **NaN filter parameters and non-positive tolerances are refused.** A NaN
  `target_cv`, or a `cv_tolerance` that is NaN, zero or negative, returns
  `Error::InvalidValue` before any identification is examined. The source
  accepts them and, because its strict comparison can never succeed, silently
  keeps only the unannotated identifications (oracle `tolerance_zero`,
  `tolerance_negative`, `tolerance_nan`, `target_nan`, and
  `target_nan_tolerance_infinite`, which shows that even an infinite tolerance
  cannot match a NaN target). A NaN target is refused for the reason
  `CompensationVoltage::new` refuses NaN, so no voltage
  `get_compensation_voltages` returns is refused as a target. An infinite
  target is not refused (see the preserved conventions above). A positive
  infinite tolerance stays valid and keeps every identification (oracle
  `tolerance_infinite`).
- **Non-numeric annotations are refused**, naming the identification index. An
  empty `DataValue` makes the source throw `Exception::ConversionError` (oracle
  `empty_datavalue_annotation`; C++ `MetaInfo` stores the empty value, and
  `metaValueExists` reports it). A string or list `DataValue` makes the source
  read an inactive union member (CPP-058, not executed).
- **NaN annotations cannot occur**: `MetaValue` refuses non-finite floats at
  construction, so every float annotation is finite.
- **Ceilings.** `MAX_SPECTRA` equals `ImTypes::MAX_SPECTRA` (100,000,000), and
  `MAX_PEPTIDE_IDENTIFICATIONS` is the same value. Both are checked before
  anything is collected or allocated. The filter's output allocation is
  fallible. The source has no ceiling.
- **`Result` return types.** The source functions cannot fail except through
  the conversion exception; the Rust functions return `Result` for the refusals
  above.

## Checked boundaries and evidence

| Evidence | Tier | What it covers |
| --- | --- | --- |
| `FAIMSHelper_test.cpp:32-104` literals | 3 (source review) | the constructor and destructor sections; 19 spectra and the voltages {-65, -55, -45} of `IM_FAIMS_test.mzML`; FAIMS detected beyond the first spectrum with the sentinel ignored; empty for non-FAIMS data |
| `oracle_cases.tsv`, `faims_file_*` | 1 (executed differential) | the C++ `MzMLFile` load of `IM_FAIMS_test.mzML`: native ID, MS level, drift time (hexadecimal) and unit of all 19 spectra, and the resulting voltage set in iteration order |
| `oracle_cases.tsv`, `faims_cvs` | 1 (executed differential) | 16 synthetic experiments: the two class-test experiments, empty, order, signed zeros, infinities, sentinel only, sentinel on another unit, unit `NONE`, duplicates, adjacent doubles, NaN first/middle/last and NaN on another unit. The Rust voltages match bit for bit, including the sign of zero, for the 13 cases without a FAIMS NaN; the three NaN cases are refused and their C++ results asserted |
| `oracle_cases.tsv`, `faims_filter*` | 1 (executed differential) | 15 `faims_filter` records: 13 parameter cases on nine identifications (double, integer, unannotated, other key, boundary values), including NaN and infinite targets at the default and at an infinite tolerance, one empty input, and an empty `DataValue` annotation. Kept lists match for the nine accepted cases, the four infinite-target cases among them. The five refused cases and the conversion exception are asserted against their recorded C++ outcomes. Two `faims_filter_infinite_annotation` records filter `+inf`, `-inf` and `-45` annotations by `±inf` at an infinite tolerance: C++ keeps only the unannotated identification, and the port keeps the same on the representable subset |
| oracle stderr (manifest) | 1 | the missing-voltage warning occurs three times: after the class-test sentinel case, after `only_sentinel`, and after `nan_first` (whose empty result shows `erase` removed the NaN root) |
| `tests/faims_helper.rs` native tests | 4 | ordering, equality and hashing of `CompensationVoltage`; NaN refusal with the spectrum index; parameter and annotation refusals; infinite targets accepted; strict boundary, order preservation and unchanged input; the ceilings |

The oracle is product-sdk (Debug, core `4fdec46`), accepted as a
development-time oracle. `git diff 4fdec46 bc9cc12` is empty for
`FAIMSHelper.{h,cpp}`, `FAIMSHelper_test.cpp`, `IM_FAIMS_test.mzML`,
`MzMLHandler.cpp`, `DataValue.cpp`, `MetaInfo.cpp`, `IMTypes.h` and
`Constants.h`, and the installed `FAIMSHelper.h` hashes equal to the pin. No
Debug-only precondition is on these paths.

## Class-test section accounting

| Section | Rust test | Status |
| --- | --- | --- |
| `FAIMSHelper()` | `constructor_and_destructor_sections` | passes |
| `~FAIMSHelper()` | `constructor_and_destructor_sections` | passes |
| `getCompensationVoltages(PeakMap& exp)` on `IM_FAIMS_test.mzML` | `get_compensation_voltages_section_on_the_values_the_cpp_reader_produced` | passes: the Rust reader loads the 19 spectra, each receives the drift time and unit the C++ reader produced, and the three literals and the oracle set hold |
| same, end to end through the Rust mzML reader | `get_compensation_voltages_section_through_the_mzml_reader` | **ignored until A3-FORMAT-IO**: the file's voltages are spectrum-level `MS:1001581` cvParams, which the Rust reader drops today. Measured failure: spectrum 0 keeps the sentinel `-1.0` (bits `0xBFF0000000000000`) where C++ reads `-55.0` |
| detects FAIMS beyond first spectrum and ignores sentinel | `get_compensation_voltages_detects_faims_beyond_the_first_spectrum_and_ignores_the_sentinel_section` | passes |
| returns empty for non-FAIMS | `get_compensation_voltages_returns_empty_for_non_faims_section` | passes |

Every section's expectations pass unchanged. The third section's reader half is
the only part that waits for another work package.

## Deferrals

- **End-to-end section 3** waits for A3-FORMAT-IO's spectrum- and scan-level
  drift time in the mzML reader. When it lands, remove the `#[ignore]` on
  `get_compensation_voltages_section_through_the_mzml_reader`; nothing else in
  this package changes.
- **Logging.** When kernel modules are wired to `LogStream`, the warning can be
  forwarded there; until then callers relay `CompensationVoltages::warnings`.

## Source defects observed

- **NaN compensation voltages break `std::set<double>`** (C++ issue candidate).
  `getCompensationVoltages` inserts drift times without a NaN check. With a NaN
  first, the result is empty and a spurious missing-voltage warning is logged;
  otherwise the NaN is silently dropped. The executed evidence is in
  `oracle_cases.tsv` (`nan_first`, `nan_middle`, `nan_last`). A NaN reaches this
  path from any reader that accepts `NaN` as a cvParam value.
- **Silent no-match parameters in `filterPeptidesByFAIMSCV`.** A tolerance of
  zero or less, or a NaN target or tolerance, keeps only unannotated
  identifications without any diagnostic, even at an infinite tolerance. This
  is a usability hazard rather than a crash; the port refuses these values. An
  infinite target also matches no annotation, but it is a voltage the source's
  own `getCompensationVoltages` can return, so the port filters it as the
  source does.
- **Non-numeric `FAIMS_CV` annotations** reach `DataValue::operator double()`,
  already logged as CPP-058.
