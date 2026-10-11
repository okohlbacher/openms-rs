# FLASHDeconv output files support

Ports of `FORMAT/FLASHDeconvSpectrumFile.h` and
`FORMAT/FLASHDeconvFeatureFile.h` (core `bc9cc12`). Rust:
`src/format/flash_deconv_spectrum_file.rs`,
`src/format/flash_deconv_feature_file.rs`. Tests:
`tests/flash_deconv_files.rs` (5). Provenance:
`tests/data/flash_deconv_files_provenance.json`.

**Status: partial.** Only the header writers are ported. Every row writer
takes records of `ANALYSIS/TOPDOWN` (`DeconvolvedSpectrum.h`, `PeakGroup.h`,
`FLASHHelperClasses.h`: `MassFeature`, `PrecalculatedAveragine`), which have no
Rust counterpart; they are the named dependency of the remaining members.

## API mapping

| Source member | Rust |
|---|---|
| `FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(os, ms_level, detail, report_decoy)` | `FLASHDeconvSpectrumFile::write_deconvolved_masses_header(writer, ms_level, detail, report_decoy)` |
| `FLASHDeconvSpectrumFile::writeTopFDHeader(os, param)` | `FLASHDeconvSpectrumFile::write_top_fd_header(writer, &Param)` |
| `FLASHDeconvSpectrumFile::writeDeconvolvedMasses(dspec, os, file_name, avg, decoy_avg, tol, write_detail, record_decoy, noise_decoy_weight)` | not ported: needs `DeconvolvedSpectrum`, `PeakGroup`, `PrecalculatedAveragine` |
| `FLASHDeconvSpectrumFile::writeMzML(map, deconvolved_spectra, deconvolved_mzML_file, annotated_mzML_file, mzml_charge, tols)` | not ported: needs `DeconvolvedSpectrum` |
| `FLASHDeconvSpectrumFile::writeIsobaricQuantification(os, deconvolved_spectra)` | not ported: needs `DeconvolvedSpectrum` |
| `FLASHDeconvSpectrumFile::writeTopFD(dspec, os, filename, qval_threshold, min_ms_level, randomize_precursor_mass, randomize_fragment_mass)` | not ported: needs `DeconvolvedSpectrum`, `PeakGroup` |
| private `topFD_min_peak_count_` (3), `topFD_max_peak_count_` (500) | documented on the type; used only by `writeTopFD` |
| `FLASHDeconvFeatureFile::writeHeader(os, report_decoy = false)` | `FLASHDeconvFeatureFile::write_header(writer, report_decoy)` |
| `FLASHDeconvFeatureFile::writeTopFDFeatureHeader(os, ms_level)` | `FLASHDeconvFeatureFile::write_top_fd_feature_header(writer, ms_level)` |
| `FLASHDeconvFeatureFile::writeFeatures(mass_features, file_name, os, report_decoy)` | not ported: needs `FLASHHelperClasses::MassFeature` |
| `FLASHDeconvFeatureFile::writeTopFDFeatures(deconvolved_spectra, mass_features, scan_rt_map, file_name, os, ms_level)` | not ported: needs `DeconvolvedSpectrum`, `MassFeature` |

## Preserved source conventions

- Headers are spelled exactly as the source's stream statements, including
  its irregular trailing tabs: the MS1 summary header ends `Qscore2D\t\n`
  without decoys and `Qvalue\n` with them; every other variant keeps the tab
  before the newline.
- `ms_level == 1` selects the MS1 variant; every other level (0 included)
  the MS2 variant, as the source's `if/else`.
- `writeTopFDHeader` writes each entry's own name (the source's
  `ParamIterator` `name`, without section path) and the `ParamValue` stream
  text, in the tree's traversal order; `writeTopFDFeatureHeader` flushes.

## Native differences

- Write errors are returned as `Error::Io`; the source writes to unchecked
  streams.

## Checked boundaries and evidence

- Tier 3: the header `START_SECTION`s of `FLASHDeconvSpectrumFile_test.cpp`
  (`writeDeconvolvedMassesHeader`, `writeTopFDHeader`) and
  `FLASHDeconvFeatureFile_test.cpp` (`writeHeader`,
  `writeTopFDFeatureHeader`) with their literals; all eight spectrum-header
  variants and all feature headers compared in full with the text transcribed
  from `FLASHDeconvSpectrumFile.cpp:236-324` and
  `FLASHDeconvFeatureFile.cpp:18-44`.
- The remaining class-test sections (`writeIsobaricQuantification`,
  `writeMzML`, `writeTopFD`, `writeFeatures`, `writeTopFDFeatures`) are not
  reproduced; they need the unported records.
