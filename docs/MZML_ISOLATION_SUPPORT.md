# Isolation-target precursor selection

The existing `LoadOptions.scientific.precursor_mz_selected_ion = false` now executes across scientific stream loading, path loading, retaining and nonretaining consumer loading. The existing count/setup reader already supports this source option; its scientific logic is unchanged. No new public type, dependency, parser, or writer option is introduced.

```rust
use openms::format::mzml::{self, LoadOptions, ReadOptions};
let mut options = LoadOptions::default();
options.scientific.precursor_mz_selected_ion = false;
let experiment = mzml::load_with_options("input.mzML", &options, &ReadOptions::default())?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Source events and representation

The implementation follows [MzMLHandler.cpp 1777–1810](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L1777) and [2157–2176](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L2157), with the original first-selected-ion rule for [user parameters](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L3534).

- An isolation target CV (`MS:1000827`) immediately overwrites `Precursor.mz`. For spectra it immediately tests the configured precursor range. Failure is sticky for the whole spectrum, including later targets or different precursors that pass. Ranges retain the source half-open interval `[min,max)` and existing unusual bound behavior.
- The first selected ion's `MS:1000744` never changes effective m/z or performs the isolation-mode range test. If its value differs from current `Precursor.mz`, it is stored in `precursor.cv_terms.metadata["selected ion m/z"]` as a finite, unit-free `MetaValueData::Float`. This uses existing typed precursor metadata, not the spectrum/chromatogram legacy String metadata maps.
- Equal values do nothing. In particular, an equal later CV does not erase an earlier differing stored value. Repeated differing selected-ion events replace that value. A target arriving after the first selected ion still overwrites effective m/z; no final-state recomputation replaces this encounter order.
- Without a target event, effective m/z remains zero and no target range test is synthesized. Missing selected-ion values create no extra metadata. The redundant `isolation_target_mz` field is cleared at precursor completion when equal to effective m/z, preserving the existing native representation convention.
- Later selected ions' scientific CV/userParam contents are ignored in isolation mode, as in source. Their markup, attribute syntax, list counts, known group references and cumulative parameter expansions remain checked. A later unconsumed `NaN` value need not be parsed; malformed XML or a missing referenced group still fails.
- Chromatogram precursors choose and retain values the same way, but precursor-range events do not exclude the entire chromatogram. Product isolation CVs do not become precursor filters.

Inline and referenced parameters follow the same path. The preexisting native selected-ion alias `MS:1000040` remains accepted; the source branch names only `MS:1000744`, so alias tests are native compatibility evidence.

The default `precursor_mz_selected_ion = true` behavior is unchanged, including existing native rejection of multiple selected ions and duplicate scientific quantities. This increment allows source later-ion ignoring and repeated target/selected-m/z CV events specifically in isolation mode. It does not claim full source reader or schema coverage. Other existing quantity, unit, formula, descriptor and model boundaries remain in force.

## Count and consumer behavior

The source count reader can count a spectrum at its scan RT event before a later precursor excludes it. Thus a target before RT gives count/filter parity, while a failing target after RT can leave `read_size_with_options` counting a record which full loading excludes. This source ordering is retained and directly tested; count parsing is not reimplemented for this option. [Count support](MZML_COUNTS_SUPPORT.md) describes early stops and missing-RT corrections.

Consumer setup still advertises raw expected counts. Its delivery pass uses the same effective-target reader and filtering, regardless of pool size, retaining mode, or `fill_data`. Input paths use the existing plain/gzip/bzip2 content-magic transport. Array population, consumer callback side effects and atomic retained destination publication follow [consumer support](MZML_CONSUMER_SUPPORT.md).

## Existing writer boundary

This group changes no writer. The existing writer serializes the new Float key as generic precursor activation `userParam` metadata; it does not yet choose that value for the selected-ion CV. It continues to emit effective `Precursor.mz` as selected ion. With no retained isolation target/offsets, it may emit no isolation window. A default-mode reread can retain the effective value and generic Float key, but this is not source selected-ion transport parity.

A concrete tested boundary is target250/selected999 without offsets: isolation loading returns effective250 plus Float999; current writing emits selected-ion250 and activation metadata999. Isolation-mode rereading then creates selected-ion metadata250 before encountering the duplicate activation key and returns a checked error. This limitation remains explicit until the separate source writer operation group consumes the key correctly. No round-trip guarantee is inferred from read support. [The independent fixture](../tests/data/mzml_isolation_target.mzML) supplies target500/selected499.5 for integration with that next group.

## Bounds, errors and evidence

Every effective spectrum target range event consumes the existing selection work allowance before execution. Every inline/reference application consumes the existing parameter count/byte ledger. A differing selected value additionally precharges a sparse typed metadata-map node and key before insertion. Parameter reuse cannot reset these counters; comparisons on later ignored scientific values are avoided. Existing XML, record, array, acquisition and consumer allowances remain separate and unchanged. Byte limits are logical conservative accounting, not physical RSS ceilings. Repeated values can consume storage allowance even when replacing an existing map entry.

Excluded records retain existing validation: a failing target does not turn later consumed invalid values or decoded malformed data into success. Errors publish neither a returned experiment nor a partially replaced/retained caller experiment; consumer callbacks already completed remain externally observable. New tests also close an existing native XML gap by rejecting literal `<` in raw attribute values before unescaping, including ignored/group content. Escaped `&lt;` remains valid. This is native malformed-input validation, not a discovered C++ defect.

[The 15 direct tests](../tests/mzml_isolation.rs) cover source literal500/499.5 values, absent/equal/repeated events, half-open endpoints, first/later ions, multiple precursors, referenced parameters, encounter order, chromatograms, default preservation, three path compression modes, both population modes, consumer pools, cumulative limits, atomic failure and the current writer boundary. The synthetic fixture uses two source class-test values from [MzMLFile_test.cpp 1420–1421](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/tests/class_tests/openms/source/MzMLFile_test.cpp#L1420); it is not an unchanged source XML fixture or executed C++ output. False-mode expectations are independently derived from the source event branches. [Provenance](../tests/data/mzml_isolation_provenance.json) preserves exact source hashes and fixture identity. No C++ execution is claimed.
