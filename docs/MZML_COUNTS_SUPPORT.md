# mzML size counting

`format::mzml::{read_size, read_size_with_options, load_size, load_size_with_options}` implement the `MzMLFile::loadSize` operation against OpenMS4-core at 82ce5b373c97f934ffd9b1ffd80215ca66473d0b. They return `MzMLCounts { spectra, chromatograms }` without constructing peaks or an experiment. Stream APIs take `BufRead`; path APIs reuse content-magic plain/gzip/bzip2 loading. The `_with_options` forms take `&PeakFileOptions` and `&ReadOptions`. Existing full readers, scientific loaders, writers and their validation behavior are unchanged.

Default counting returns declared record-list counts. It does not verify the number of children, and it stops successfully at the second list's start once both counts are known. If only one list is present, it parses to EOF and reports zero for the absent list. Signed source 32-bit counts are retained internally and negative counts report zero; the exact `-1` sentinel prevents the usual early stop. Source `StringUtils` numeric syntax, including the unusual `+-1`, is reused. `metadata_only=true` stops before the first list's count is read and returns zero counts, after consuming its required `defaultDataProcessingRef`.

RT, MS-level or precursor-m/z filters select the source filtered-count algorithm. This is an event algorithm, so its result can differ from the length of a fully loaded, filtered experiment:

- A mismatching MS-level CV skips a spectrum. Missing MS-level metadata does not reject it.
- An accepted scan-start-time CV increments the spectrum count and skips subsequent XML events for that record. RT ranges are half-open; minutes are converted to seconds, while other unit labels follow the source seconds branch. Missing RT and elution-time-only metadata do not count.
- Precursor filtering only has an effect before the accepted RT. Only the first selected ion is used; an equal selected/isolation m/z bypasses the source selected-ion range check. `precursor_mz_selected_ion=false` instead checks the isolation target.
- One referenceable parameter-group expansion finishes all its CV terms even after the record becomes skipped. Multiple accepted RT terms in one group can therefore count a spectrum more than once. Separate inline RT events count at most once. This finite source quirk is preserved.
- Chromatograms are counted at their start; spectral filters do not reduce this count.

Peak m/z/intensity ranges do not select filtered counting or affect its result. Sorting, write settings, Numpress choices, `always_append_data`, and data-pool sizing do not affect these APIs. Neither `fill_data` nor `skip_xml_checks` enables binary decoding or disables XML validation here. Full-record metadata and header definitions unrelated to counting are not materialized or scientifically validated; no header roundtrip or consumer operation is claimed.

Two deliberate native corrections are included. A missing-RT spectrum never causes incidental binary decoding, whereas C++ may leave it in a normal data pool and decode it despite producing no count. `skip_chromatograms` skips chromatogram counting while preserving spectrum counting; the source initializes a shared skip flag that can suppress header and spectrum callbacks too. These corrections have source-derived tests, not C++ runtime differential claims.

The reader validates UTF-8 and XML 1.0 characters, matching markup, mzML element namespaces, duplicate attributes, consumed parameter group IDs/references, and required count-related fields. It does not apply a full mzML schema or validate skipped binary encodings, array counts or scientific metadata. Unknown/malformed consumed references and invalid/nonfinite consumed numbers are checked errors rather than source warnings or unsafe behavior. Unconsumed record fields remain unparsed, while malformed XML in those regions is still rejected. All declared and observed record counts are bounded, including ignored records; limits can therefore reject input the source would skip.

UTF-8, optional UTF-8 BOM, and US-ASCII are supported. ISO-8859-1 declarations are accepted only when the consumed bytes are ASCII, which permits the exact original class-test fixture. Non-ASCII ISO-8859-1 returns `Unsupported`; no Latin-1 transcoding is claimed. BOM handling is independent of input chunk size. DTDs and external entities are unsupported. The five predefined entities, valid character references and bounded CDATA are supported as text.

The implementation keeps quick-xml authoritative for markup, namespace scopes and tag matching. Between markup/reference events it discards ordinary text directly in chunks of at most 8192 bytes, validating UTF-8 and XML legality, including forbidden literal `]]>` across chunks. It stops before `<` or `&` so markup and references stay with quick-xml. This avoids buffering a large `<binary>` text event. The input wrapper tracks raw consumed XML bytes; direct reads intentionally do not use quick-xml's diagnostic offset. Errors use the existing non-positional mzML error convention.

Resource limits are cumulative:

| Bound | Counting use |
|---|---|
| `ReadOptions.max_xml_bytes` | All consumed decompressed XML bytes, including discarded text and BOM. |
| `max_records` | Declared list totals, observed physical records, and returned counts, including reference-induced multiple counts. |
| `max_param_groups` | Declared and actual stored parameter-group counts. |
| `max_total_params` | Every consumed start/empty-element descriptor plus each expanded parameter. |
| `max_param_bytes` | Temporary event/namespace allowances, attribute maps, stack/group/reference storage and repeated expansions. |
| `MAX_COUNT_EVENT_BYTES` |1 MiB maximum per buffered markup/comment/CDATA/reference event, further reduced when the remaining storage budget is smaller. |
| `MAX_COUNT_XML_DEPTH` |256 nested elements, including skipped payload. |
| `MAX_COUNT_WORK` |50 million semantic event, expanded-term, reference-ID and conservative MS-level membership visits. |

An event's conservative storage allowance is reserved before quick-xml reads or allocates it; unused allowance is returned, while actual costs accumulate. Array byte/element/peak limits and acquisition materialization mode are inapplicable because no arrays/acquisitions are produced. The fixed event buffer cap applies to comments and CDATA even when their content is irrelevant, while ordinary text can span the full permitted XML input.

Successful early return does not drain the XML or compressed tail. A trailing malformed body or gzip checksum can remain unexamined, unlike full-reader EOF validation. Underlying buffers/decompressors may prefetch; this is a logical parse/consume boundary, not a promise of zero physical reads beyond an XML token. ZIP files retain the shared path reader's explicit unsupported policy. Returning counts by value avoids partially published output on error; a borrowed input stream has naturally advanced.

`tests/mzml_counts.rs` covers the five literal source count pairs on the unchanged `MzMLFile_1.mzML`, independently derived count/parameter/precursor branches, checked native corrections, encoding and tiny-chunk parser state, large text under small storage limits, resource failures, path compression, and unconsumed I/O/checksum tails. `tests/data/mzml_counts_provenance.json` records source/fixture hashes and test provenance. No C++ execution was used to create expected values.

Consumer/transform methods, experiment-wide `ExperimentalSettings`, full metadata-only record loading, centroid-info consumers, and arbitrary header transport remain separate work. This operation does not advertise those capabilities.
