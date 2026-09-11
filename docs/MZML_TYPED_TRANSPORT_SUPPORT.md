# Typed records, primary roles and independent noise grids

This operation group ports the retained source MzMLHandler record metadata and additional primary-array routes at SDK `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. It builds on scientific loading, streaming consumers and the writer-options group. It does not claim the whole MzMLFile header complete or provide an XSD validator.

`MSSpectrum.metadata` and `MSChromatogram.metadata` now use the existing `MetaInfo`/`MetaValue` model. There is one metadata owner per record. Existing string values stay strings: `"42"` and `"[1, 2]"` are never implicitly parsed. All seven native value alternatives and units can be held in memory. See [the migration guide](RECORD_METADATA_MIGRATION.md).

## Source transport

Record user parameters use the existing XML scalar decoder. Integer, Float and String values, exact integer width, escaping and units are preserved; the private record-name marker requires a unit-free String. All ten spectrum CV-to-metadata routes and eight scalar scan routes are implemented. Source string-valued mass resolution, preset scan configuration and mass resolving power remain Strings. Elution time becomes a unit-free Float in seconds, multiplying minutes by 60. Without a scan RT event, source close-time fallback copies that record value to RT before primary metadata merges; it does not create an RT filtering event. Integer/Float fallback values are accepted; nonnumeric values are rejected rather than using the source inactive-union cast (CPP-058). An explicit scan RT leaves an unused fallback key untouched. Scan direction and scan law use source literal Strings. `spotID` becomes `maldi_spot_id`; export may use a typed user parameter instead of the original attribute. Scan-level CV metadata may likewise be exported as a record user parameter without duplicating acquisition metadata.

Primary-array metadata is moved into the same record map. For spectra the coordinate description merges before the intensity description, regardless of XML order. Chromatograms merge descriptions in XML order. Later values replace earlier values as in source MetaInfo addition. Primary processing history is rejected because the native primary values have no independent history owner. An unused wavelength array remains an ordinary aligned float array, including its description/history and its original encounter order.

Supported selectors are unit-free Strings:

| Record key | Value | Stored axis | Source CV and unit |
| --- | --- | --- | --- |
| `mzml coordinate array` | `wavelength` | spectrum coordinate | MS:1000617 / UO:0000018 |
| `mzml intensity array` | `absorption` | intensity | MS:1000515 / UO:0000269 |
| `mzml intensity array` | `pressure` | chromatogram intensity | MS:1000821 / UO:0000109 |
| `mzml intensity array` | `flow` | chromatogram intensity | MS:1000820 / UO:0000270 |
| `mzml intensity array` | `nonstandard` | chromatogram intensity | MS:1000786, value `detector signal` / UO:0000000 |

Without a selector, the existing m/z/time and intensity axes apply. Redundant or unknown selector strings are rejected on writing. m/z wins over wavelength regardless of encounter order. An empty unused wavelength is a valid empty auxiliary; a selected primary wavelength must have the declared peak count. MS level zero is valid only for the three typed non-mass modes ElectromagneticRadiation, Emission and Absorption. Unknown and mass-spectrometry modes still require a positive MS level.

Role values are not numerically rescaled. Absent or canonical additional-role units are accepted. Existing m/z/default-intensity unit acceptance remains unchanged by this group. A nondefault role's supplied unit is retained through the source `unit_accession` String key, except absorption, whose selector carries the role. Conflicting units are rejected rather than silently relabelled. Export supplies canonical role units, so reloading input which omitted units can add the canonical `unit_accession` key. Pressure/flow/nonstandard chromatogram type accessions MS:1003019, MS:1003020 and MS:1000626 use the existing `chromatogram type accession` String with the default typed Mass chromatogram kind; conflicting typed kinds are rejected.

## Independent noise

The spectrum keys `sampled noise m/z array`, `sampled noise intensity array`, and `sampled noise baseline array` each own a unit-free FloatList. Their lengths are independent of one another and of the peak count, including zero. Reading removes them from aligned auxiliary arrays before peak filtering/sorting. Their values are retained in f64 and remain unchanged by those operations. `fill_data=false` still validates descriptors but publishes no binary-derived primary metadata, arrays or noise lists.

All writers emit present noise keys, in the fixed order above, between the two primary arrays and the aligned auxiliaries. Each noise array has its own `arrayLength`, Float64 precision, and ordinary encoding; requested Numpress is suppressed for noise while ordinary zlib is retained. Noise does not consume an auxiliary description/header slot. Primary precision, indexing, TPP compatibility, checksums and selected-ion behavior remain those of [writer options](MZML_WRITE_OPTIONS_SUPPORT.md).

Noise descriptions or processing histories are rejected because the source FloatList representation cannot own them. Spectrum aligned arrays with reserved noise names are rejected before output, as are chromatogram aligned pressure/flow/detector arrays which would become competing primary intensities. The inverse domains remain aligned auxiliaries. A wavelength auxiliary cannot compete with a wavelength primary. Multiple promoted intensity arrays are rejected, avoiding source [CPP-053](../OpenMS_CPP_ISSUES.md) promotion/first-match mismatch.

## Checked boundaries and evidence

Generic metadata lists and Empty values are rejected on XML writing: source serializes them as strings and loses their type. The three noise FloatLists have the dedicated route above. Finite-value, XML, descriptor, duplicate-owner, array-length and existing codec checks remain enforced. Parsed metadata insertion/merge uses existing parameter/header ledgers; metadata copies and replacement destruction are precharged before mutation. Newly supported noise encoding is additionally precharged in the existing fixed metadata write allowance (50 million work units / 256 MiB conservative storage estimate). This is an accounting bound, not an allocator-level memory measurement. Reader binary and writer codec/selection limits remain independent and cumulative as before. Stream I/O failures may leave partial output; scientific/preflight failures precede external output.

Typed metadata is metered before full-record copies in the existing processing, theoretical append and EMG copy paths. Operations which move or ignore metadata do not gain extra copy/update accounting beyond their ordinary record validation. Annotation publishes Float tolerance and Integer ppm; alignment publishes Float `original_RT` and FloatList `original_rt`, preserving any preexisting value. MGF retains only unit-free String metadata; DTA now rejects any record metadata before writing. These are explicit transport boundaries.

The deterministic [source projection tool](../tools/generate_mzml_typed_reference.py) verifies 43 source hashes and projects 13 original source sections (1,092 lines) plus 26 complete CV stanzas. The [fixture](../tests/data/mzml_typed_transport_source.json) includes the original noise/PDA/pressure class-test literals and all metadata routing sections. [Direct transport tests](../tests/mzml_typed_transport.rs) combine those literals with independent ordering, typed-value, resource, filtering, callback and writer ownership checks. No C++ execution or byte-identical XML claim is made. [Record tests](../tests/record_metadata.rs) cover the model and flat formats.
