# Ion mobility types support

Native coverage of `IONMOBILITY/IMTypes.h` with its two implementation files,
and of the four ion-mobility members plus the three spectrum-type conversions of
`METADATA/SpectrumSettings.h`, at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Artifact | Path |
| --- | --- |
| Implementation | `src/metadata/im_types.rs`, `src/metadata/acquisition.rs` |
| Tests | `tests/im_types.rs` (24 tests) |
| Manifest | `tests/data/im_types_provenance.json` |

This work package closes the residual the spectrum-mobility package recorded
verbatim as *"Not covered: setIMFormat/getIMFormat/getIMPeakType/setIMPeakType,
which need MSSpectrum fields this package may not add"*. The four accessors are
now ported on `SpectrumSettings`, the class that declares them, together with
the whole of `IMTypes.h`, which had no ledger entry and no Rust module at all.
One part of that residual survives and is stated in **Deferrals**: the Rust
`MSSpectrum` flattens `SpectrumSettings` instead of inheriting it, and giving it
its own copy of the two fields means editing `src/kernel.rs`, which this work
package may not do.

## API mapping

### `IONMOBILITY/IMTypes.h`

| C++ member | Rust |
| --- | --- |
| `enum class DriftTimeUnit` | `metadata::DriftTimeUnit`, ported earlier in `src/metadata/acquisition.rs` |
| `DriftTimeUnit::NONE` / `MILLISECOND` / `VSSC` / `FAIMS_COMPENSATION_VOLTAGE` / `CCS` | `None` / `Millisecond` / `InverseReducedMobility` / `FaimsCompensationVoltage` / `CollisionCrossSection` |
| `DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT` | not ported: a count sentinel. `DriftTimeUnit::ALL.len()` is the count and `NAMES_OF_DRIFT_TIME_UNIT.len()` the array length |
| `extern const std::string NamesOfDriftTimeUnit[]` | `metadata::NAMES_OF_DRIFT_TIME_UNIT` |
| `DriftTimeUnit toDriftTimeUnit(const std::string&)` | `metadata::to_drift_time_unit` |
| `const std::string& driftTimeUnitToString(DriftTimeUnit)` | `metadata::drift_time_unit_to_string`; `DriftTimeUnit::name` and its `Display` are the same string |
| `enum class IMFormat` | `metadata::IonMobilityFormat`, ported earlier in `src/metadata/acquisition.rs` |
| `IMFormat::NONE` / `IM_PEAK` / `IM_SPECTRUM` / `UNKNOWN` | `None` / `PerPeak` / `PerSpectrum` / `Unknown` (the `Default`) |
| `IMFormat::SIZE_OF_IMFORMAT` | not ported: a count sentinel, as above |
| `extern const std::string NamesOfIMFormat[]` | `metadata::NAMES_OF_IM_FORMAT` |
| `IMFormat toIMFormat(const std::string&)` | `metadata::to_im_format` |
| `const std::string& imFormatToString(IMFormat)` | `metadata::im_format_to_string` |
| `enum class IMPeakType` | `metadata::IonMobilityPeakType`, ported earlier in `src/metadata/acquisition.rs` |
| `IMPeakType::IM_PROFILE` / `IM_CENTROIDED` / `UNKNOWN` | `Profile` / `Centroid` / `Unknown` (the `Default`) |
| `IMPeakType::SIZE_OF_IMPEAKTYPE` | not ported: a count sentinel, as above |
| `extern const std::string NamesOfIMPeakType[]` | `metadata::NAMES_OF_IM_PEAK_TYPE` |
| `IMPeakType toIMPeakType(const std::string&)` | `metadata::to_im_peak_type` |
| `const std::string& imPeakTypeToString(IMPeakType)` | `metadata::im_peak_type_to_string` |
| `class IMTypes` | `metadata::ImTypes`, a zero-sized unit struct: the source class holds no state and exists only to scope the statics |
| `IMTypes()`, `~IMTypes()` (implicit) | derived `Default`; no drop glue |
| `constexpr double IMTypes::DRIFTTIME_NOT_SET` | `ImTypes::DRIFTTIME_NOT_SET`. `kernel::spectrum_mobility` and `kernel::ranges` each hold a private copy of the same `-1.0`; `tests/im_types.rs` asserts the public constant agrees with `MSSpectrum::has_drift_time` |
| `constexpr double IMTypes::N2_BUFFER_GAS_MASS` | `ImTypes::N2_BUFFER_GAS_MASS` |
| `static IMFormat determineIMFormat(const MSExperiment&, int ms_level)` | `ImTypes::determine_im_format_for_ms_level` |
| `static IMFormat determineIMFormat(const MSSpectrum&)` | `ImTypes::determine_im_format` for the data half, and `ImTypes::determine_im_format_with_stored` for the complete function including the stored-format short circuit. Two functions because the Rust `MSSpectrum` carries no stored format; see **Native differences** |
| `static DIM_UNIT fromIMUnit(DriftTimeUnit)` | **not ported:** `DIM_UNIT` (`CONCEPT/CommonEnums.h`) has no Rust counterpart. `src/kernel/ranges.rs` documents that the port's dimension model is coarser than the source's `DIM_UNIT`, which separates the three ion-mobility units that `kernel::ranges::MSDim` merges into one mobility dimension. Nothing in the port can consume the return value, so the conversion would have no caller |
| `static double oneOverK0ToCCS(double, double, int, double buffer_gas_mass = N2_BUFFER_GAS_MASS)` | `ImTypes::one_over_k0_to_ccs` (N2 default) and `ImTypes::one_over_k0_to_ccs_with_buffer_gas` (explicit gas) |
| `static double ccsToOneOverK0(double, double, int, double buffer_gas_mass = N2_BUFFER_GAS_MASS)` | `ImTypes::ccs_to_one_over_k0` and `ImTypes::ccs_to_one_over_k0_with_buffer_gas` |
| anonymous-namespace `MASON_SCHAMP_CONSTANT` (`IMTypes.cpp:112`) | private `MASON_SCHAMP_CONSTANT` in `src/metadata/im_types.rs`, private as in the source |
| anonymous-namespace `reducedMass_` (`IMTypes.cpp:115`) | private `reduced_mass` |
| forward declarations `class MSExperiment`, `class MSSpectrum` | not applicable: Rust needs no forward declaration |

### `METADATA/SpectrumSettings.h`

The four ion-mobility members and the three spectrum-type conversions are this
work package's scope. Every other member is listed so that the header's public
surface is completely accounted for; those rows name where the member lives and
say plainly that this work package did not audit it. `docs/METADATA_SUPPORT.md`
is the METADATA package's own document for them.

| C++ member | Rust |
| --- | --- |
| `void setIMFormat(const IMFormat&)` | `SpectrumSettings::set_im_format` — **new here** |
| `IMFormat getIMFormat() const` | `SpectrumSettings::im_format` — **new here** |
| `void setIMPeakType(IMPeakType)` | `SpectrumSettings::set_im_peak_type` — **new here** |
| `IMPeakType getIMPeakType() const` | `SpectrumSettings::im_peak_type` — **new here** |
| `protected IMFormat im_type_ = IMFormat::UNKNOWN` | public field `SpectrumSettings::ion_mobility_format`, same default |
| `protected IMPeakType im_peak_type_ = IMPeakType::UNKNOWN` | public field `SpectrumSettings::ion_mobility_peak_type`, same default |
| `static const std::string NamesOfSpectrumType[]` | `SpectrumSettings::NAMES_OF_SPECTRUM_TYPE` — **new here** |
| `static StringList getAllNamesOfSpectrumType()` | `SpectrumSettings::all_names_of_spectrum_type` — **new here** |
| `static const std::string& spectrumTypeToString(SpectrumType)` | `SpectrumSettings::spectrum_type_to_string` — **new here** |
| `static SpectrumType toSpectrumType(const std::string&)` | `SpectrumSettings::to_spectrum_type` — **new here** |
| `enum class SpectrumType { UNKNOWN, CENTROID, PROFILE }` | `kernel::SpectrumType`, ported earlier; `SIZE_OF_SPECTRUMTYPE` is a count sentinel and is not ported |
| `SpectrumSettings()`, copy/move constructors, `~SpectrumSettings()`, copy/move assignment | derived `Default` and `Clone`; Rust moves need no declaration. Not audited here |
| `bool operator==` / `operator!=` | derived `PartialEq`. The source comparison omits `im_type_` and `im_peak_type_`; the port includes both, a divergence `docs/METADATA_SUPPORT.md` already records and `tests/metadata.rs` asserts. Not audited here |
| `void unify(const SpectrumSettings&)` | `SpectrumSettings::unify`. It leaves both ion-mobility fields local, as the source does — the source's `unify` never mentions them. Not otherwise audited here |
| `SpectrumType getType() const` / `void setType(SpectrumType)` | public field `SpectrumSettings::spectrum_type`. Not audited here |
| `getNativeID` / `setNativeID` | public field `native_id`. Not audited here |
| `getComment` / `setComment` | public field `comment`. Not audited here |
| `getInstrumentSettings` (const and mutable) / `setInstrumentSettings` | public field `instrument_settings`. Not audited here |
| `getAcquisitionInfo` (const and mutable) / `setAcquisitionInfo` | public field `acquisition_info`. Not audited here |
| `getSourceFile` (const and mutable) / `setSourceFile` | public field `source_file`. Not audited here |
| `getPrecursors` (const and mutable) / `setPrecursors` | public field `precursors`. Not audited here |
| `getProducts` (const and mutable) / `setProducts` | public field `products`. Not audited here |
| `setDataProcessing` / `getDataProcessing` (mutable) / `getDataProcessing` (const, constified copy) | public field `data_processing`, owned values rather than shared pointers. Not audited here |
| `protected type_`, `native_id_`, `comment_`, `instrument_settings_`, `source_file_`, `acquisition_info_`, `precursors_`, `products_`, `data_processing_` | the correspondingly named public fields. Not audited here |
| `operator<<(std::ostream&, const SpectrumSettings&)` | **not ported: no `Display`.** The source prints two banner lines and nothing else — `"-- SPECTRUMSETTINGS BEGIN --"` and `"-- SPECTRUMSETTINGS END --"` (`SpectrumSettings.cpp:185-190`), ignoring its argument entirely. Reported as a finding rather than reproduced |
| `std::hash<OpenMS::SpectrumSettings>` | **not ported.** `tests/metadata_hash.rs` covers `Product` and the CV term lists, not `SpectrumSettings`. Reported as a finding |

## Preserved source conventions

* **Name tables and their order.** The three arrays keep the source's order and
  spelling, `"<NONE>"` and `"1/K0"` included, and `tests/im_types.rs` asserts
  each array equals the enum's own `name()` sequence so the duplicated literals
  cannot drift from the single definition.
* **Exact, case-sensitive name matching.** The source conversions are
  `std::find` over the array, so `"centroid"` is not `"Centroid"` and
  `"IM_PEAK"` is not `"im_peak"`. The port matches byte for byte.
* **Determination precedence.** An ion-mobility float array wins over a
  per-spectrum drift time. A spectrum carrying both is `IM_PEAK` and is *valid*:
  the header claims `@throws Exception::InvalidValue` for that combination, but
  `IMTypesExperiment.cpp:66-72` only writes a debug record, and
  `IMTypes_test.cpp:183-186` asserts `IM_PEAK`. The implementation and the test
  are reproduced; the header comment is a source defect, recorded below.
* **The stored-format short circuit.** Any stored format other than `UNKNOWN` is
  returned without looking at the data — including a stored `NONE` on a spectrum
  that does carry ion mobility. `determine_im_format_with_stored` reproduces
  that exactly, and `tests/im_types.rs` pins all four cases.
* **`NONE` is erased before the per-level verdict.** A level holding both plain
  spectra and ion-mobility spectra takes the format of the ion-mobility ones; a
  level with none is `NONE`, and so is a level with no spectra at all.
* **The unreachable single-format guard.** The source's `occs.size() == 1`
  branch rejects anything but `IM_PEAK` and `IM_SPECTRUM` with *"subfunction
  returned invalid value(s)"*. That cannot happen — the per-spectrum function
  returns only `NONE`, `IM_PEAK` or `IM_SPECTRUM`, and `NONE` was just erased —
  and it is reproduced anyway, so a future change to the classification cannot
  silently return an unusable format.
* **The charge sign is ignored** in both cross-section conversions, and the
  rounded `28.0` is kept for N2 rather than the physical `28.006148`, because the
  calibration constant `1059.62245` was validated against that rounding.
* **Multiplication order.** `(C * |z|) / sqrt(mu) * (1/K0)` and
  `ccs * sqrt(mu) / (C * |z|)` are written exactly as the source evaluates them,
  so equal inputs give bit-identical results; the round trip through both is
  therefore accurate to one rounding rather than exact, as the source's own test
  acknowledges with `TEST_REAL_SIMILAR`.
* **No parallelism to reproduce.** Neither `IONMOBILITY/*` nor
  `METADATA/SpectrumSettings.cpp` contains a `#pragma omp`, so the port's serial
  policy costs nothing here.

## Native differences

* **No `SIZE_OF_*` sentinel, so four conversions cannot fail.**
  `driftTimeUnitToString`, `imFormatToString` and `imPeakTypeToString` throw
  `Exception::InvalidValue` only when handed the sentinel, and
  `spectrumTypeToString` likewise. A Rust enum cannot hold it, so those four are
  infallible and return `&'static str`. The `@throws` clauses are neutralised,
  not dropped.
* **`determineIMFormat(const MSSpectrum&)` is split in two.** The source reads
  the spectrum's own stored format first, through the `SpectrumSettings` base.
  The Rust `MSSpectrum` flattens those settings and has no stored format, so
  `determine_im_format(spectrum)` is the data half and
  `determine_im_format_with_stored(stored, spectrum)` is the whole function for
  a caller that holds the stored value on a `SpectrumSettings`. Passing
  `IonMobilityFormat::Unknown`, the source default, makes the two identical.
* **`determine_im_format_for_ms_level` classifies by the data half**, for the
  same reason. It agrees with the source for every spectrum whose stored format
  is the `UNKNOWN` default, which is every spectrum the port can build.
* **Non-finite input is rejected.** The source guards `value <= 0.0`, which a
  NaN passes — every comparison against NaN is false — so `oneOverK0ToCCS(NaN,
  …)` returns NaN. Both conversions require finite positive values here.
* **The buffer gas mass is checked.** The source never examines it, although a
  non-positive mass makes the reduced mass negative and the result NaN.
* **`i32::MIN` is a valid charge.** The source takes `std::abs` of an `int`,
  undefined for `INT_MIN`; the port uses `i32::unsigned_abs`.
* **`ms_level` is `i32`**, as in the source, and a negative level matches no
  spectrum — the source compares with `std::cmp_not_equal` against the
  spectrum's unsigned level, so the sign is preserved rather than wrapped.
* **The two log records are not emitted.** The source writes a debug record for
  a spectrum carrying both annotations and a warning for a drift time without a
  unit. No kernel or metadata module in this crate is wired to
  `concept::log_stream::LogStream`, so neither is written; both conditions stay
  directly observable at the call site (`spectrum.contains_im_data() &&
  spectrum.has_drift_time()`, and `spectrum.drift_time_unit ==
  DriftTimeUnit::None`), and the module documentation says so.
* **`all_names_of_spectrum_type` allocates only because the source's
  `StringList` does.** `NAMES_OF_SPECTRUM_TYPE` is the same list without the
  copy.
* **Equality includes the ion-mobility fields** where the source's does not.
  That divergence predates this work package; `docs/METADATA_SUPPORT.md` records
  it and `tests/metadata.rs` asserts it. It is restated on
  `SpectrumSettings::im_format` so that a reader of the accessor sees it.

## Checked boundaries and evidence

| Boundary | Value | Where |
| --- | --- | --- |
| `ImTypes::MAX_SPECTRA` | 100 000 000 spectra | `determine_im_format_for_ms_level` |

The ceiling matches `RangeManager::MAX_ITEMS` and
`MSSpectrum::MAX_MOBILITY_ITEMS`, so an experiment whose ranges can be computed
can also be asked for its ion-mobility format. It is checked before anything is
collected. Nothing in this module mutates its input: every function takes
`&self`-free arguments or a shared reference, so there is no partial-update path
to make atomic. The per-spectrum determination allocates nothing at all; the
per-level one allocates a `BTreeSet` that can hold at most four elements,
because that is the enum's cardinality.

The ceiling itself is not exercised by a test: allocating 100 000 001 spectra is
not a test. `tests/im_types.rs` asserts the constant's value and that the scan
works below it, which is the same choice `tests/spectrum_mobility.rs` makes for
`MSSpectrum::MAX_MOBILITY_ITEMS`.

**Evidence tiers.** Every expectation carrying a section number is tier 3
(source review): the literals of the 12 `START_SECTION`s of `IMTypes_test.cpp`,
the three ion-mobility sections of `MSSpectrum_test.cpp` and the
`getAllNamesOfSpectrumType` section of `SpectrumSettings_test.cpp`, transcribed
into `tests/im_types.rs`. The reserpine cross-section value `244.9402` and its
`TOLERANCE_ABSOLUTE(0.01)` come from `IMTypes_test.cpp:202-206`. The finiteness
and buffer-gas checks, the `i32::MIN` charge, the negative MS level, the name
table agreement and the ceiling are tier 4 (Rust-only invariants). No C++ was
built or executed, and no retained C++ output exists for these headers.
`tests/data/im_types_provenance.json` hashes the eleven source files and pins
the source anchors behind each documented quirk.

## Class-test section accounting

`IMTypes_test.cpp` has 12 `START_SECTION`s. All 12 are ported into
`tests/im_types.rs`, one Rust test each, none merely mapped.

| Source section | Rust test |
| --- | --- |
| `IMTypes()` | `im_types_default_constructor` |
| `~IMTypes()` | `im_types_destructor` |
| `DriftTimeUnit toDriftTimeUnit(const std::string&)` | `to_drift_time_unit_covers_every_name_and_rejects_others` |
| `const std::string& driftTimeUnitToString(const DriftTimeUnit)` | `drift_time_unit_to_string_covers_every_variant` |
| `IMFormat toIMFormat(const std::string&)` | `to_im_format_covers_every_name_and_rejects_others` |
| `const std::string& imFormatToString(const IMFormat)` | `im_format_to_string_covers_every_variant` |
| `IMPeakType string conversions` | `im_peak_type_string_conversions` |
| `static IMFormat determineIMFormat(const MSExperiment&, int)` | `determine_im_format_per_ms_level` |
| `static IMFormat determineIMFormat(const MSSpectrum&)` | `determine_im_format_per_spectrum` |
| `determineIMFormat returns IM_PEAK for centroided IM data` | `determine_im_format_for_centroided_im_data` |
| `static double oneOverK0ToCCS(double, double, int, double)` | `one_over_k0_to_ccs` |
| `static double ccsToOneOverK0(double, double, int, double)` | `ccs_to_one_over_k0_round_trips` |

`MSSpectrum_test.cpp` has 71 `START_SECTION`s. 68 are ported in
`tests/spectrum_mobility.rs` and the remaining three — the ones that package
could only map onto the bare enums — are ported here:

| Source section | Rust test | Substitution |
| --- | --- | --- |
| `void setIMFormat(IMFormat imf)` (`:1453`) | `set_im_format_round_trips` | asserted on `SpectrumSettings`, which declares the member, instead of on `MSSpectrum`, which inherits it in C++ |
| `IMPeakType getIMPeakType() const` (`:1463`) | `get_im_peak_type_defaults_to_unknown` | as above |
| `void setIMPeakType(IMPeakType)` (`:1470`) | `set_im_peak_type_round_trips` | as above |

`getIMFormat` has no section of its own upstream; `get_im_format_defaults_to_unknown_and_directs_to_the_determination`
covers it natively, including the header's note that an `UNKNOWN` value should be
resolved from the data.

`SpectrumSettings_test.cpp` has 33 `START_SECTION`s. One is ported here —
`static StringList getAllNamesOfSpectrumType()` as `all_names_of_spectrum_type`,
because this package added the members it exercises. None of its sections
touches the ion-mobility quartet. The other 32 exercise members the METADATA
package implemented; `docs/METADATA_SUPPORT.md` carries no section accounting,
so where each of those 32 is asserted is not established by any document. That
is a gap in that package's evidence and is reported as a finding, not silently
absorbed into this package's count.

Self-audit (`IMTypes.h` and the seven `SpectrumSettings.h` members): 16 sections
newly ported, 0 mapped-with-evidence, 0 mapped-without-evidence, 0 unaccounted.

### Fixture substitutions

`IMTypes_test.cpp` builds its `IMwithFDA` fixture by running
`IMDataConverter::reshapeIMFrameToSingle` over a one-spectrum experiment.
`IMDataConverter` is not ported, so `tests/im_types.rs` builds the reshaped
spectrum directly: one float data array named `"raw inverse reduced ion mobility
array"`, which is what `IMDataArrayUtils::setIMUnit` writes for the fixture's
`VSSC` unit (`IMDataArrayUtils.cpp:28`), and no per-spectrum drift time. Both
`determineIMFormat` overloads look only at `containsIMData()` and
`getDriftTime()`, so the substitution cannot change any asserted value.

## Deferrals

* **`MSSpectrum` still has no stored ion-mobility format or peak type.** In C++
  `MSSpectrum` inherits `SpectrumSettings` and so answers `getIMFormat()` and
  `getIMPeakType()` itself; twenty-odd call sites in `FORMAT/`,
  `PROCESSING/CENTROIDING/PeakPickerIM.cpp` and `FEATUREFINDER/` use them on a
  spectrum. The Rust `MSSpectrum` flattens the settings into named fields, and
  adding `ion_mobility_format` and `ion_mobility_peak_type` to it means editing
  the struct and its `Default` in `src/kernel.rs` (struct at `src/kernel.rs:255`,
  `Default` at `src/kernel.rs:283`), which this work package may not do. The
  whole surface those fields need is in place: the enums, the four accessors on
  the declaring class, and `ImTypes::determine_im_format_with_stored`, which
  takes the stored value as its first argument precisely so that the fields can
  be wired in without changing this module.
* **`IMTypes::fromIMUnit` is not ported**, because its return type `DIM_UNIT`
  is not ported; see the API mapping row.
* **`std::hash<SpectrumSettings>` and `operator<<(std::ostream&, const
  SpectrumSettings&)` are not ported.** Both are public members of
  `SpectrumSettings.h`, outside this package's residual, and neither has a Rust
  counterpart today.
* **`KERNEL/DimMapper.h` is deferred to the TOPPView package.** Basis: zero core
  includers; its six consumers are all `desktop/gui/source/VISUAL/*.cpp`
  (`PlotCanvas`, `LayerData*`). It is TOPPView code sitting in the `KERNEL`
  directory, so this is a package boundary, not a kernel gap. This is the
  KERNEL domain's one deliberate scope exclusion.
* `tests/data/im_types_provenance.json` is not yet listed in
  `SOURCE_PROVENANCE.json`, so `tools/check_core_sdk.py` does not verify its
  hashes. Registering it is the integrator's step.

## Source defects observed

* `IMTypes.h:100-107` documents `@throws Exception::InvalidValue if IM values are
  annotated as single drift time and float array` on
  `determineIMFormat(const MSSpectrum&)`. The implementation throws nothing for
  that input: `IMTypesExperiment.cpp:66-72` writes a debug record and returns
  `IM_PEAK`, and `IMTypes_test.cpp:183-186` asserts exactly that. The
  documentation is wrong, not the code.
* `IMTypesExperiment.cpp:38-44` guards a case that cannot occur, as described
  under **Preserved source conventions**. Harmless, but the error message
  *"subfunction returned invalid value(s)"* reports `occs.size()`, which is 1 on
  that path, so the message would be misleading if it ever fired.

## Validation

```sh
cargo fmt --all -- --check
cargo clippy --locked --all-features --lib --test im_types -- -D warnings
cargo nextest run --locked --all-features --test im_types
cargo nextest run --locked --no-default-features --test im_types
cargo test --locked --all-features --doc
RUSTDOCFLAGS="-D warnings" cargo doc --locked --all-features --no-deps
cargo +1.85.0 check --locked --all-features --lib --test im_types
python3 tools/check_doc_coverage.py --report | grep im_types
```

CI: `tests/im_types.rs` needs no features and belongs on the
`--no-default-features` kernel line of the `minimum-rust` job
(`.github/workflows/rust.yml:79`), appended as `--test im_types`.
