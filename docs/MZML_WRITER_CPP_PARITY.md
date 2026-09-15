# mzML writer: scale and C++ Release output parity

Work package `mzml-writer-scale-parity`. It removes the writer's scale blocker
and brings the written document close to the C++ Release writer's, measured
with the benchmark harness's decoded comparison (decision D6). It changes no
peak data: the decoded arrays of `MzMLSplitter`, `SpectraFilterWindowMower` and
`MapNormalizer` on the benchmark slice stay bitwise identical to the C++ output
(1200 of 1200 arrays each), before and after.

Sources are the pinned core `bc9cc12`
(`src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp`, `MzMLHandlerHelper.cpp`,
`include/OpenMS/FORMAT/OPTIONS/PeakFileOptions.h`). Executed evidence is the C++
Release build `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`
on the staged benchmark inputs, with the INIs of the 2026-09-14 smoke run.

## 1. Ceilings scale with the experiment

The writer preflights — the header plan, the settings validation and the
prepared peak-file writer's markup, index and binary budgets — each drew on one
fixed whole-document allowance. Realistic runs exhausted them: the smoke
benchmark measured 647 spectra passing and 648 failing on `UK222_picked`, and
on `integrate/wave2` `SpectraFilterWindowMower` already failed on the
600-spectrum slice with `controlled vocabulary resource limit exceeded`, while
the C++ writer stores the complete 40,856-spectrum runs.

Every such allowance is now multiplied by the experiment's *record shares*
(`mzml::writer_shares`: one for the experiment-level header plus one per
spectrum and chromatogram), so a ceiling is linear in what it has to cover and
cannot refuse a long run. Inside one document, amplification is still bounded,
and an explicit `PeakWriteLimits` value passed to
`write_with_peak_options_and_limits` or `store_with_peak_options_and_limits`
remains a whole-document ceiling that is never reset per record.
`PeakWriteLimits::for_experiment`, used by the entry points without explicit
limits, scales the same defaults and adds a fixed allowance per stored array
value (1 KiB of work and 256 bytes) plus per-array Numpress ceilings that grow
to the longest array. The source enforces no ceilings at all.

Measured on the staged inputs (release build, `dax`, shared node, one run;
`tests/mzml_writer_scale.rs`, the `#[ignore]`d HPC tests):

| Input | Spectra | Read | Store, indexed | Store, plain | Output |
|---|---:|---:|---:|---:|---:|
| `centroid_lcms_qe_silac_uk222_picked/UK222_picked.mzML`, 547 MB | 40,856 | 3.99 s | 3.65 s | 1.83 s | 533 MB |
| `profile_hr_qe_silac_uk222/UK222.mzML`, 2.32 GB | 40,856 | 11.59 s | 15.44 s | 4.76 s | 3.34 GB |

The output of the profile run is larger than its input because the input stores
32-bit m/z and both writers' defaults store 64-bit m/z (`mz_32_bit = false`).

## 2. Indexed mzML by default, in one pass

`mzml::write`, and therefore `FileHandler::store_experiment` and every TOPP
tool that stores mzML, now writes `indexedmzML` with record offsets and a
`fileChecksum`, as `MzMLFile::store` does with its default `PeakFileOptions`
(`write_index_ = true`, `PeakFileOptions.h:244`); the C++ Release tools write
`indexedmzML` on all benchmark inputs.

It is a single pass: each record's offset is taken as its `<spectrum` or
`<chromatogram` tag starts and the checksum is computed while the bytes go out,
so the markup is rendered once. The prepared two-pass writer
(`write_with_peak_options`) is unchanged and still produces the same bytes; the
plain layout that `MSDataWritingConsumer` splits per record stays available as
`write_with_options`. An experiment with neither spectra nor chromatograms has
nothing to index and stays plain, where the source writes an index with a
fabricated `-1` offset (CPP-050).

The checksum is a real SHA-1 over the bytes from the start of the document
through the opening `<fileChecksum>` tag, as the indexed schema specifies; the
source writes the constant `0` (CPP-049). It is what the extra store time above
costs: about 10.7 s of the 15.44 s on a 3.34 GB document, roughly 310 MB/s,
because `sha1` is pinned with `force-soft`. The C++ writer pays nothing for its
placeholder, so this is a real work difference in any timing comparison.

Checked on the C++ side: C++ `FileInfo` and C++ `MzMLSplitter` both read the
Rust `MapNormalizer` output of the benchmark slice (exit 0; 600 spectra,
186,536 peaks, split into 2 parts), and `tools/mzml_writing/check_output.py`
verifies independently, on the real tool outputs, that every offset addresses
its record's opening tag in order and that the digest is the SHA-1 of the
prefix.

## 3. No invented precursor intensity

The writer emitted `MS:1000042` for every precursor. The source writes it only
for a positive intensity (or when the meta value `peak intensity` exists), and
always with unit attributes, `MS:1000132` unless
`peak intensity unit accession` names another unit (`MzMLHandler.cpp:4596`).
On the benchmark slice the C++ output keeps the input's 149 terms while this
port wrote 198, adding `value="0"` to 49 precursors — the decoded comparison
reported exactly those 49 as `precursor term only in a`. A precursor without
the term reads back as `0.0`, so the default `+0.0` is now omitted; a negative
or `-0.0` intensity, and an explicitly united one, are still written, because
dropping them would change what reads back.

## 4. Metadata parity

Changed to follow the source writer:

| Item | Now | Source |
|---|---|---|
| Root element | `xsi:schemaLocation`, `accession` (even empty), `version`; no `id` | `MzMLHandler.cpp:4851` |
| `mzML@id` | written as the run `userParam name="mzml_id"` | `:1201`, `:5229` |
| `cvList` | the five pinned entries with the source's names, versions and URIs | `:4855` |
| `softwareList` | `so_in_0`, `so_configuration_<i>`, `so_default`, then the processing methods' software | `:5107` |
| `so_default` | the source's empty `Software()`: `version=""`, `MS:1000799 value=""` | `:5116`, `:3763` |
| `processingMethod@order` | always `0` | `:3852` |
| Run element | `id="ru_0"`, `defaultSourceFileRef` when the run has a source file | `:5211`, `:5219` |
| `dataProcessingRef` | omitted where the record's history is the list default | `:5257` |
| Valueless `cvParam` | no `value` attribute | `writeCV_`, `:3600` |
| Numbers in `cvParam`/`userParam` | `StringUtils::toStr(double)` text where it reads back unchanged, so an inherited `3.0`, `200.0` or `1.0e20` survives verbatim | `DataValue::toString`, `StringUtils.cpp:384` |
| Intensity array | carries `MS:1000131 number of detector counts` | `:5688` |

The empty-history placeholder that the schema forces is now its own software
entry, `so_default_empty_history`, so that `so_default` can be the source's
empty `Software()`; its name, version and marker userParam are unchanged, and
reading still normalises only that exact payload back to an empty history.

`float_text` keeps the source text only when it reads back as the same value.
The source's fixed branch writes 15 fraction digits, which loses precision
above 14 significant digits; those values keep Rust's shortest round-tripping
text rather than being truncated.

### Result of the decoded comparison

Rust versus C++ Release, same tool, INI and input
(`inputs/derived/sub_centroid_uk222_picked_first600.mzML`, and
`sub_profile_uk222_first600.mzML` for `BaselineFilter`), compared with the
benchmark harness's `equiv.py` (data and metadata verdicts, every difference
counted):

| Tool | Before | After |
|---|---|---|
| `MzMLSplitter` | data DIFFERENT (49 invented precursor terms), 2433 metadata differences in 25 categories | data equal within tolerance, 9 metadata differences in 7 categories |
| `SpectraFilterWindowMower` | did not run (`controlled vocabulary resource limit exceeded`) | data equal within tolerance, 13 metadata differences in 9 categories |
| `MapNormalizer` | data DIFFERENT (49), 2436 metadata differences in 25 categories | data equal within tolerance, 11 metadata differences in 7 categories |
| `BaselineFilter` | data DIFFERENT (49 terms plus 189 of 1200 arrays), 2433 metadata differences | data DIFFERENT (189 of 1200 arrays, a separate lane), 9 metadata differences in 7 categories |

Peak arrays stay bitwise identical for the three I/O tools throughout. The
remaining `equal within tolerance` is the precursor intensity, whose f32 value
the source prints with more digits (largest relative difference 5.2e-8).

### What still differs, and why

- **Generated header identifiers** (6 to 8 differences per file): this port
  writes `sf_<20 digits>`, `dp_<20 digits>` and `so_dp_<20 digits>_<20 digits>`
  where the source writes `sf_ru_0`, `dp_sp_0` and `so_dp_sp_0_pm_0`, so
  `sourceFile@id`, `dataProcessing@id`, `processingMethod@softwareRef`,
  `software@id`, `run@defaultSourceFileRef` and
  `spectrumList@defaultDataProcessingRef` differ in value while referring to
  the same elements. `MSDataWritingConsumer` splits the writer's records on
  the literal `dp_00000000000000000000`, so changing the scheme also changes
  `src/format/ms_data_writing_consumer.rs`, which this lane does not own; it is
  in the integrator requests.
- **`spectrum@dataProcessingRef` on the first spectrum** (1 difference): the
  source repeats the default reference on spectrum 0 and this port omits it
  everywhere it is the default, which also keeps the streaming consumer's
  per-record rendering byte-identical to the whole-document writer.
- **Software aliases**: for a tool name the pinned CV does not contain
  verbatim, the source retries `"<name> software"` and `"TOPP <name>"`, so
  `SpectraFilterWindowMower` becomes `MS:1002146 TOPP SpectraFilterWindowMower`
  while this port writes `MS:1000799` with the exact name. The alias does not
  read back as the name it was given; the exact-name transport is deliberate
  (see [header support](MZML_HEADER_SUPPORT.md)).
- **Container**: `indexListOffset` and `fileChecksum` differ by construction
  (real offsets and a real digest against the source's `0`). D6 reports these
  as container facts, not metadata differences.
- **Layout**, out of scope under D6: the XML declaration (`UTF-8` against the
  source's `ISO-8859-1` for ASCII-only bytes), indentation, line breaks and
  attribute order.
- **Element order inside a record**: this port writes the spectrum type after
  the MS level and sorts metadata by name, where the source follows its
  `MetaInfoRegistry` order. The decoded comparison keys parameters by
  accession or name within their record, so this is not a reported difference;
  restoring the source order would need an ordered `MetaInfo`.

## Evidence

- `tests/mzml_writer_scale.rs`: 50,000 realistic records through every default
  entry point (they failed before this package at about 650); the former fixed
  ceiling still refusing that experiment when passed explicitly; a 20,000-record
  file stored through `FileHandler` and read back equal; the indexed default
  verified by `tools/mzml_writing/check_output.py`; the header parity literals
  above; and two `#[ignore]`d HPC tests on the staged 547 MB and 2.32 GB inputs.
- The existing writer suites (`tests/mzml_write_options.rs`,
  `tests/mzml_header.rs`, `tests/mzml_settings.rs`, `tests/mzml_acquisition.rs`,
  `tests/mzml_precursor_activation.rs`, `tests/ms_data_writing_consumer.rs`)
  keep their assertions; the schema checks now validate indexed output against
  the pinned indexed XSD.
- The decoded comparison above was run on `dax` with the harness's own
  `equiv.py` and `tools.json` defaults; the driver is
  `../oracle/mzml-writer-scale-parity/`.
