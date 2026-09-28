# MzQCFile support

Source header (core `bc9cc12`): `FORMAT/MzQCFile.h`, implementation
`FORMAT/MzQCFile.cpp`; with it `FileHandler::computeFileHash` from
`FORMAT/FileHandler.h`.

Rust: `src/format/mzqc_file.rs` (feature `mzml`), the shared
`nlohmann::json` writer `src/format/nlohmann_json.rs`, and
`FileHandler::compute_file_hash` in `src/format/file_handler.rs`. The metrics
come from `crate::analysis::qc` (`docs/QC_SUPPORT.md`), the metric names from
PSI-MS and the embedded PSI QC vocabulary `resources/cv/qc-cv.obo`.
Tests: `tests/mzqc_file.rs`. Provenance: `tests/data/mzqc_file_provenance.json`.

The feature gate is `mzml` because it is the feature that brings the `sha1`
crate the input-file checksum needs, and an mzQC document always describes an
mzML input; no new feature or dependency was added.

## API mapping

| Source | Rust |
|---|---|
| `MzQCFile()` (defaulted) | `MzQcFile::new`, `Default` |
| `store(input_file, output_file, exp, contact_name, contact_address, description, label, feature_map, prot_ids, pep_ids)` | `MzQcFile::store(output_file, &MzQcInput { input_file, experiment, contact_name, contact_address, description, label, feature_map }, prot_ids, pep_ids) -> Result<MzQcReport>` |
| (the document text, built inside `store`) | `MzQcFile::text(input, prot_ids, pep_ids, creation_date) -> Result<(String, MzQcReport)>` |
| `FileHandler::computeFileHash` | `FileHandler::compute_file_hash` |
| `FileHandler::storeQC` | not ported: its qcML branch needs `QcMLFile::collectQCData`; its mzQC branch is exactly `MzQcFile::store` |

`MzQcInput` groups the seven non-identification parameters so the call has
four arguments; `prot_ids` is a shared slice (the source's non-const
reference is never written through), `pep_ids` stays mutable because the
missed-cleavage metric annotates the first hits, as in the source.

## Preserved source conventions

- Metric selection by `QCBase::Status`: raw-data metrics need a non-empty
  input file name, feature metrics a non-empty feature map, identification
  metrics non-empty protein *and* peptide identifications. The number of
  chromatograms, RT duration and m/z range are always written.
- Metric order, accessions and value shapes as the source writes them
  (table in `MzQcFile::text`); the RT mean shift is omitted when it is 0.
- A metric whose accession the vocabularies lack is skipped with the source's
  `<accession> not found in CV.` line (returned in `MzQcReport::log`).
- The two vocabularies are loaded as the source loads them:
  `psi-ms.obo` as `PSI-MS`, then `qc-cv.obo` as `QC`, into one vocabulary.
- Input-file record: absolute location (`File::absolutePath`, generic
  separators), base name, file format `MS:10000584` (sic: one digit too many
  for PSI-MS `MS:1000584`), completion time `getDate() + "T" + getTime()`
  (`0000-00-00T00:00:00` when unset), SHA-1, instrument name.
- Analysis software `MS:1009001` (no such PSI-MS term; the source's comment
  asks for one) named `QCCalculator`, `https://www.openms.de`.
- Controlled vocabularies: PSI QC 0.1.2 and PSI-MS 4.1.155, always.
- `computeFileHash`: an unreadable file hashes to the empty string, and a
  failing read ends the input.
- Text: `nlohmann::ordered_json::dump(2)` - insertion order, two spaces,
  shortest round-trip doubles with `.0` on integral values and two-digit
  exponents, `null` for non-finite values (a mean over nothing), no final
  newline; `qualityMetrics` is `null` if no metric was added. Plain bytes
  whatever the output extension, as `std::ofstream` writes them.

## Native differences

- Atomic output: the document is built before anything is written and
  published by rename, so a failing metric leaves an existing file unchanged;
  the source opens (truncates) the output first and leaves it empty.
- Ranges are computed on demand from the experiment; the source reads the
  cached ranges of `updateRanges()`, a documented precondition. Only the
  `sortSpectra()` precondition remains the caller's.
- RT duration and m/z range of an experiment without spectra (or peaks) are
  0 and `[4294967295, 0]`; the source converts infinite or `DBL_MAX` bounds
  to `UInt`, which is undefined.
- The TIC area is exact beyond 2^32 and the ten-fold jump test does not wrap
  (`docs/QC_SUPPORT.md`).
- The analysis-software version is the core SDK version this crate
  implements (`CORE_SDK_VERSION`, 4.0.0); the source writes its library
  version. The upstream test whitelists every `version` line.
- The creation date is the local time from the crate's `DateTime::now`;
  `MzQcFile::text` takes it as a parameter for reproducible text.
- `store` returns an error for an invalid spectrum (`MSExperiment::ranges`),
  for a non-numeric `rt_deviation` and for identifications without a usable
  enzyme; the source throws for the latter two and does not check the first.

## Checked boundaries and evidence

- Tier 1, `tests/mzqc_file.rs::qccalculator_2_reproduces_the_retained_mzqc`:
  TOPP_QCCalculator_2 (test-data 0cb15f2 `topp/CMakeLists.txt:5387-5389`)
  reproduced from its retained inputs through the library as
  topp 174b576 `src/QCCalculator.cpp` drives it, compared with
  `QCCalculator_2_output.mzQC` under `FuzzyDiff -ini FuzzyDiff.ini` (ratio
  1.01, absdiff 0.01) with the whitelist `creationDate`, `location`,
  `version`. Every line outside the whitelist is byte-identical except one
  number, see below.
- The experiment is read from `QCCalculator_input_strict_reader.mzML`, the
  retained mzML with the eight `dataProcessingRef` attributes on its m/z,
  intensity and time arrays deleted: the strict mzML reader refuses
  processing on a primary array (`docs/FILE_INFO_SUPPORT.md`), while the
  source's handler drops it. The input-file record is still written for the
  retained file, whose SHA-1 the output records.
- Oracle `oracle/F2-mzqc` (Linux x86_64 Release build
  `openms4-release-bc9cc12-c19e494-174b576`, ibminode06, gcc 14.4.0, run
  twice): see "Oracle" below.
- Tier 1, `text_is_the_release_builds_output_line_for_line`: the port's
  document equals the Release build's (`tests/data/mzqc_file/release/
  QCCalculator_2_release.mzQC`) line for line, the creation date included,
  except the input's absolute location, which names the directory each ran
  in.
- Tier 4: the whitelisted lines hold the given creation date, the absolute
  location, mzQC 1.0.0 and the SDK version; a document without inputs keeps
  the three unconditional metrics and logs the four notes in order; an
  unreadable file has an empty checksum; a failed store leaves the file.

## Oracle

`oracle/F2-mzqc` (beside the repository; `driver.cpp`, `remote_run.sh`,
`run.sh`, `results/hashes.txt`) builds QCCalculator's mzQC body against the
Linux x86_64 Release install `openms4-release-bc9cc12-c19e494-174b576` on
ibminode06 (gcc 14.4.0, `-O3 -DNDEBUG`), and runs it twice in fresh processes
on the retained inputs, once with the retained mzML and once with the
strict-reader derivative. Results (`results/differ.txt` empty):

- run 1 and run 2 agree apart from `creationDate`;
- the retained and the derived mzML give the same document apart from
  `creationDate`, so the derivative changes nothing the source reads;
- the RT mean shift of the pinned build is `0.916678249835968`
  (`0.91667824983596802` at 17 digits), the value the port writes. The
  retained upstream output, written by an older build (its software version
  line says 3.5.0), has `0.9166780114173889`; the relative difference,
  2.6e-7, is inside the upstream FuzzyDiff tolerance. No summation order of
  the six `float` deviations gives the retained value, so it predates the
  pinned `FeatureSummary` arithmetic or input.

The oracle also prints the six features in the order `sortByRT` leaves them,
with their deviations (`results/values.tsv`); summed in that order in `f32`
they give the pinned build's value, which is the port's.
