# Spectrum meta-data lookup

`metadata::spectrum_meta_data_lookup` covers `METADATA/SpectrumMetaDataLookup.h` and
`METADATA/SpectrumMetaDataLookup.cpp` at SDK `bc9cc12`. It builds on
[`SPECTRUM_LOOKUP_SUPPORT.md`](SPECTRUM_LOOKUP_SUPPORT.md) (`metadata::spectrum_lookup`).

Tests: [`tests/spectrum_meta_data_lookup.rs`](../tests/spectrum_meta_data_lookup.rs), every
`START_SECTION` of `SpectrumMetaDataLookup_test.cpp` with unchanged literals (the
`addMissingSpectrumReferences` section reads the same `MzMLFile_1.mzML` the class test does,
from `tests/data/mzml_validator/`), plus tests derived from `SpectrumMetaDataLookup.cpp` for what
the class test does not reach: `addMissingFAIMSToPeptideIDs` (no section upstream), the
precursor-RT and scan-number error lines, the index overflow, the uncaught `ElementNotFound` of
`addMissingIMToPeptideIDs`, `stop_on_error`, and reference-extraction details.

The two consumer-local subsets that predate this port (`mascot_xml::SpectrumMetaData`,
`pepxml::SpectrumMetaData`) are unchanged; they can move onto this type when their owners next
touch them.

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `typedef unsigned char MetaDataFlags`, `MDF_RT` … `MDF_ALL` | `MetaDataFlags = u8`, `MDF_*` constants | Same bit values. |
| `struct SpectrumMetaData` and its constructor | `SpectrumMetaData`, `Default` | Quiet NaN for the three times/m/z, charge and MS level 0, scan -1, empty native ID, as the source. |
| `SpectrumMetaDataLookup()`, destructor | `SpectrumMetaDataLookup::new()`, `Default`, `Drop` | |
| base class `SpectrumLookup` | `Deref<Target = SpectrumLookup>` | `find_by_*`, `is_empty`, `reference_formats`, `rt_tolerance` read through the base; `add_reference_format` and `set_rt_tolerance` mutate. The base's `read_spectra` is not reachable mutably (see below). |
| `readSpectra(spectra, scan_regexp = default_scan_regexp, get_precursor_rt = false)` | `read_spectra(&[MSSpectrum], &str, bool) -> Result<Vec<String>>` | Returns the lines the source logs. |
| `setSpectraDataRef(spectra_data)`, protected `spectra_data_ref` | `set_spectra_data_ref`, `spectra_data_ref()` | |
| `getSpectrumMetaData(Size, SpectrumMetaData&)` | `spectrum_meta_data(usize) -> Result<SpectrumMetaData>` | `IndexOverflow` becomes `Error::InvalidRange`. |
| static `getSpectrumMetaData(const MSSpectrum&, SpectrumMetaData&, const boost::regex& = boost::regex(), const std::map<Size,double>& = {})` | `spectrum_meta_data_of(&MSSpectrum, &mut SpectrumMetaData, Option<&BoostRegex>, &BTreeMap<usize, f64>) -> Result<Vec<String>>` | The empty default regex is `None`. Writes into the caller's record like the source. |
| `getSpectrumMetaData(const std::string&, SpectrumMetaData&, MetaDataFlags = MDF_ALL)` | `spectrum_meta_data_from_reference(&str, &mut SpectrumMetaData, MetaDataFlags) -> Result<()>` | Rust has no default argument; pass `MDF_ALL`. |
| static `addMissingRTsToPeptideIDs` | `add_missing_rts_to_peptide_ids(&mut [PeptideIdentification], &MSExperiment) -> Result<Annotation>` | `Annotation { success, log }` is the source's `bool` and its log lines. A missing RT is `None` or NaN. |
| static `addMissingIMToPeptideIDs` | `add_missing_im_to_peptide_ids` | |
| static `addMissingFAIMSToPeptideIDs` | `add_missing_faims_to_peptide_ids` | |
| static `addMissingSpectrumReferences(peptides, filename, stop_on_error = false, override_spectra_data = false, override_spectra_references = false, proteins = {})` | `add_missing_spectrum_references(&mut [PeptideIdentification], &str, bool, bool, bool, &mut [ProteinIdentification])` | Pass `&mut []` for the default proteins. See divergences. |
| private copy constructor and assignment | `Clone` | Harmless for a value that owns its data. |

## Divergences

- **Re-reading replaces.** The source's `readSpectra` clears neither `metadata_` nor the base
  class's retention-time, ID and scan maps, so after a second call index *i* still returns the
  first read's meta data and the first read's keys keep resolving to indexes of the new spectra.
  `read_spectra` replaces everything. A lookup read once is identical. Source defect candidate.
- **`addMissingSpectrumReferences` updates the proteins.** The source takes
  `std::vector<ProteinIdentification> proteins` by value and sets `spectra_data` on that copy, so
  `override_spectra_data` never reaches the caller. This port takes `&mut [ProteinIdentification]`
  and sets `primary_ms_run_paths` (the typed form of the `spectra_data` meta value, as the idXML
  reader stores it) to `["file://" + filename]`. Source defect candidate.
- **An empty scan expression.** The source's `setScanRegExp_` leaves the
  previous scan expression in place when given an empty one; with replacing re-reads, an empty
  expression simply disables scan numbers, as `SpectrumLookup` documents.
- **The stored ion-mobility format.** `addMissingIMToPeptideIDs` calls
  `IMTypes::determineIMFormat(spec)`, which returns a spectrum's stored `getIMFormat()` before
  looking at the data. This crate's `MSSpectrum` carries no stored format field, so the format is
  determined from the data (`ImTypes::determine_im_format`); for a spectrum whose stored format is
  the default `UNKNOWN`, which every spectrum built by a reader has, the result is the source's.
- **Numbers from references.** A charge is converted with `toDouble` and truncated to `Int`
  exactly as the source's implicit conversion, and refused outside `i32` where the source's
  conversion is undefined. A negative `LEVEL` is refused where the source wraps it to a huge
  `Size`.
- **Non-finite retention times** are refused while indexing, as `SpectrumLookup::read_spectra`
  does.
- **Log lines are returned**, not written to `OPENMS_LOG_*`; the wording is the source's.

## Evidence

Tier 3: `SpectrumMetaDataLookup_test.cpp` transcribed with unchanged literals; the other tests
are derived from `SpectrumMetaDataLookup.cpp` and cite its lines. No C++ was executed.
