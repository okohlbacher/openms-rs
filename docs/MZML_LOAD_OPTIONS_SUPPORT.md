# mzML scientific loading options

`mzml::read_with_load_options(reader, &LoadOptions, &ReadOptions)` executes the represented `PeakFileOptions` loading choices. `LoadOptions.scientific` holds scientific choices; existing `ReadOptions` continues to bound XML, raw record/point counts, parameters and decoded arrays. Existing `read` and `read_with_options` retain their order-preserving behavior. This is a loading increment, not complete `MzMLFile` coverage.

## Executed choices

- Spectrum RT and explicit MS-level selection, including signed/duplicate MS-level option values and exact membership.
- Spectrum m/z and intensity selection, and chromatogram RT and intensity selection. Range endpoints use the source half-open predicate `!(value < min || value >= max)`. Empty/equal, inverted, infinite and NaN endpoint behavior follows that predicate; options are not normalized.
- Precursor selected-ion m/z selection using the actual source event rule. A selected value differing from the current precursor m/z triggers the range check; an equal isolation target/selected value, or selected zero with an initially zero target, bypasses it. Any failing precursor event excludes that spectrum. This applies to spectra carrying a precursor irrespective of MS level; chromatograms are unaffected.
- Source-default spectrum-m/z/chromatogram-RT sorting, independently disabled by their options. Equal coordinates retain input order, a deterministic native extension to source's unspecified ties.
- `skip_chromatograms` and native `LoadOptions.skip_spectra`.

Missing MS-level and scan-start CV parameters do not trigger their filters, matching the source event-driven checks. Explicit duplicate/conflicting CVs remain errors. Parameter-group references execute the same checks as inline CVs. Scan/chromatogram minute units are converted to seconds before range selection. Record order, identifiers, supported metadata, precursor/Product information and all retained annotation arrays remain intact.

Primary intensities and auxiliary floating arrays are selected while their decoded values are still `f64`; only retained values narrow to native `f32`. The aligned selection moves every nonempty float, integer and string annotation array with its points. Empty array placeholders retain their identity. Array names/string payloads are moved rather than cloned during permutation.

## Checked boundaries

Whole-record exclusion does not waive input validation. Skipped records are fully decoded, narrowed and validated, and still count against raw XML, record, point, binary and parameter limits. Thus a malformed/unsupported excluded record still errors, and a finite auxiliary or primary value exceeding native `f32` still errors when the entire record is skipped. This deliberately differs from source SAX early skipping. Conversely, a finite `f64` value on a discarded peak of a retained record need not fit `f32`. Nonfinite binary values, type violations and mismatched array lengths are always rejected before selection.

The function supports [metadata-only early stopping](MZML_HEADER_SUPPORT.md). It rejects `fill_data = false` for ordinary record loading, `skip_xml_checks = true`, and `precursor_mz_selected_ion = false` before reading input. Consumer/transform loading, isolation-target precursor output and the other outstanding `MzMLFile` operations remain separate work. The preexisting native selected-ion alias `MS:1000040` executes the same native filter as `MS:1000744`; source's selected-ion branch names the latter.

Write-only fields do not alter a read: `write_index` (including its true default), compression/precision/Numpress configurations, MQ/TPP compatibility and supplemental-data writing. `always_append_data` has no effect without a consumer. `max_data_pool_size` is a source batching hint; this native reader processes one record at a time and does not retain a batch. No unsupported codec is enabled by setting its write configuration.

## Canonical auxiliary identity

The reader/writer recognize all 26 non-primary binary-array descendants in the pinned PSI-MS vocabulary. [The extracted table](../tests/data/mzml_load_canonical_arrays.tsv) records their exact accessions, canonical names and declared binary-type masks. Primary m/z, time and intensity retain their existing dedicated roles; non-standard arrays use their explicit names.

Those 26 exact canonical names are reserved native array roles. Canonical input ignores a caller-supplied alternate `name` in favor of the vocabulary name, as source does; output writes the corresponding canonical accession. Non-standard input explicitly using a reserved canonical name errors because the name-only native array cannot retain that distinction. Case variants and aliases are ordinary distinct names, without inference. Duplicate resulting names are rejected. Canonical declared binary types are checked; terms with no declared type restriction remain permissive. Canonical charge arrays write signed 32-bit integers to preserve their declared type, while ordinary integer annotations retain the existing signed 64-bit writer representation. String/float/integer canonical type conflicts fail writer preflight before any bytes are emitted. [Header/reference transport](MZML_HEADER_SUPPORT.md) now retains auxiliary scalar metadata, represented units and processing references. Independent primary-array history and non-string or unit-bearing primary metadata remain checked errors.

This name-based policy preserves represented canonical identity without claiming a full controlled-vocabulary validator or a richer array metadata model.

## Resources and evidence

`max_selection_work` defaults to 500 million visits and `max_selection_bytes` to 256 MiB of cumulative index/permutation allocations. One counter covers every record, precursor CV application, MS-level membership scan, array check, filter, permutation and sort allowance; it is not reset between records or group references. Scratch capacity is charged before allocation and sorting. Raw decoding and output still obey the independent existing `ReadOptions` limits. The reader returns a new experiment only after full success.

[The 18 focused tests](../tests/mzml_load_options.rs) reproduce the source's four published restricted-load cases and compare all 65 primary source points against an independent Python base64/struct projection. They also cover double precision boundaries, all array kinds and 26 canonical roles, empty placeholders, minute scaling, CV presence/reference behavior, precursor bypass, unsupported option preflight, malformed excluded data and cumulative limits. Current and minimum-version checks also exercise the existing mzML, parameter-group and Product suites (49 tests total).

[Provenance](../tests/data/mzml_load_provenance.json) pins source code, source tests, original dataset, vocabulary and generated fixture hashes at `54a232fe2cae9c590d5c997fa49d20e7769860fb`. [The unmodified original file](../tests/data/mzml_load_source_original.mzML) is retained separately from [the scientific projection](../tests/data/mzml_load_source_projection.mzML). The projection strips unrelated unsupported metadata, binary-array processing references and user parameters, and corrects one declared auxiliary-array count to the actual four arrays; it is not represented as an unmodified source load golden. Published expected scalar/count literals remain distinct from decoded IEEE values and independent edge-case arithmetic. No C++ code was executed.

Independent peer review checked source filter events, permutation/resource accounting, deferred auxiliary narrowing and all canonical table entries. No unresolved review finding remains.
