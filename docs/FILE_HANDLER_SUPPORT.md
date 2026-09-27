# FileHandler: loadExperiment's arguments and the member ledger

`FORMAT/FileHandler.h` at core `bc9cc12` is the dispatch every TOPP tool loads
and stores through. This page records, member by member, what the port's
`format::file_handler::FileHandler` covers, with the evidence, and ranks what
remains. The type registry and the lexical helpers are in
[FILE_HANDLING_SUPPORT.md](FILE_HANDLING_SUPPORT.md); `getType` and the
option-taking loaders in [MZML_MOBILITY_SUPPORT.md](MZML_MOBILITY_SUPPORT.md).

## loadExperiment with all of its arguments

Source `FileHandler::loadExperiment(filename, exp, allowed_types, log,
rewrite_source_file, compute_hash)` (`FileHandler.h:164-179`,
`FileHandler.cpp:849-1052`) is `FileHandler::load_experiment_with_progress(path,
allowed, &ExperimentLoadOptions, &ProgressLogger)`. `ExperimentLoadOptions`
carries the handler's `PeakFileOptions` (`peaks`), the native mzML reader
limits (`mzml`, no source counterpart), `rewrite_source_file` and
`compute_hash`; its default is the source's default call.

**Log type.** The source hands each file object only `setLogType(log)`. The
port hands the reader a copy of the caller's logger made as
`ProgressLogger::clone` makes one: a fresh backend of the same type that shares
the clock, the nesting context and the GUI factory. A `ProgressLogger` rather
than a bare `ProgressLogType` is taken because the port's GUI factory lives on
the logger, not in a process global as the source's `make_gui_progress_logger`
does; `ProgressLogger::new()` with `set_log_type` is the plain log-type call.

| Detected type | Source reader and its calls | Port |
|---|---|---|
| mzML | `MzMLFile` (document, spectrum-list and chromatogram-list sections) | `mzml::load_with_progress` |
| mzXML | `MzXMLFile` | `mzxml::load_with_progress` |
| mzData | `MzDataFile` | `mzdata::load_with_progress` |
| DTA2D | `DTA2DFile` | `dta2d::read_with_progress` over the handler's (possibly decompressed) stream |
| DTA, MS2 | none: DTA gets no log type (`DTAFile().load`, `FileHandler.cpp:869-875`), and `MS2File::load` makes no call | none |
| MGF | `MascotGenericFile::load`: one section sized by the file | **none**: see divergence 1 |

**SRM conversion.** Every mzML load ends with
`ChromatogramTools().convertSpectraToChromatograms<PeakMap>(exp, true)`
(`FileHandler.cpp:910`, `KERNEL/ChromatogramTools.h:110-215`), and so do all
three of the port's experiment loaders for mzML (`load_experiment`,
`load_experiment_with_options`/`_with_read_options` and
`load_experiment_with_progress`): spectra of scan mode SRM are removed, and
those with one precursor and at least one peak become one SRM chromatogram per
precursor and product m/z, appended after the chromatograms the file carried.
Each SRM spectrum that cannot be converted writes the source's warning
(`ChromatogramTools.h:164-165`) to the warning log. An experiment without SRM
spectra is left untouched after one scan, exactly the source's result. The
conversion runs under `ChromatogramConversionLimits`: the defaults bound its
own output, and the record, work and byte ceilings are raised by the cost of
scanning and moving the records the reader already loaded, so a large file
with a few SRM spectra is not refused for its size.

**rewrite_source_file** (`FileHandler.cpp:1014-1051`). The first source file
is kept, or a default one; none warns `No source file annotated.`, more than
one warns `Expecting a single source file in mzML. Found N will take only first
one for rewriting.` Its name becomes the basename of the given path after
trailing `/` and `\` are stripped; its path `file://` plus the directory part
of `system::file::absolute_path` of that name, which is lexical as the source's
`std::filesystem::absolute` is (`./x.mzML` records `<cwd>/.`), with a third `/`
when the directory does not start with one; its file type the mzML term of the
detected type (`FileType::mzml_name`, `FileTypes.cpp:259-275`). Everything else
the record carried stays, the original checksum included when no hash is
computed. That record then replaces all source files.

**compute_hash** stores `FileHandler::compute_file_hash(path)` with type
SHA-1. `compute_file_hash` is `FileHandler::computeFileHash`
(`FileHandler.cpp:768-778`): the SHA-1 of the bytes as stored (a `.gz` file's
digest is that of its compressed bytes), 40 lowercase hex digits, the empty
string for a file that cannot be opened, and the digest of the bytes read so
far when a read fails, which makes a directory's digest that of no input
(`:182-185`). It uses the `sha1` crate already behind the `mzml` feature, so
it and `compute_hash` need that feature.

**mzXML and mzData** are now dispatched by every experiment loader, with the
peak options their readers consume (`FileHandler.cpp:886-902`); they were read
by the port but not reachable through `FileHandler`. They are not yet written
through it (`can_write_experiment` stays false for both).

### Divergences, each checked against the Release build

1. **MGF reports no progress.** The native MGF dispatch is the strict
   `format::mgf` adapter, not the `format::mascot_generic` port of
   `MascotGenericFile` that has the source's progress section. The loaded
   spectra counts and the rewritten source file agree with the Release build;
   the calls do not (Release: `S 0 131 "loading MGF"`, two sets, an end).
2. **A missing `.dta2d` fails before its section starts.** The source's
   `DTA2DFile::load` starts its section and then throws `FileNotFound`,
   leaving the nesting one level deeper; the port opens the file before it
   dispatches and makes no call.
3. **`compute_hash` without `rewrite_source_file` is refused**
   (`Error::InvalidValue`, before any I/O). The source guards it with an
   `OPENMS_PRECONDITION` (`FileHandler.cpp:853`), which only a debug build
   checks; the Release build ignores the flag, which its `mzml_hash_only` case
   shows.
4. A type outside the allowed list is `Error::InvalidValue` where the source
   throws `ParseError`; neither makes a progress call.

### Evidence

`tests/file_handler_load.rs`, 11 tests:

- the upstream class test `FileHandler_test.cpp:229-320` (mzML and DTA2D
  checksums, the DTA2D m/z range with and without the hash, the allowed-type
  refusal, the mzData and mzXML m/z-range loads, a non-ASCII path) and RFC
  3174's SHA-1 vectors;
- `tests/data/file_handler_load_release.tsv`, the C++ Release build at core
  `bc9cc12` on ibminode06 (`../oracle/file-handler-load`, provenance in
  `tests/data/file_handler_load_provenance.json`, two repetitions
  identical): 23 loads, each on a fresh `FileHandler` with log type GUI and a
  recording backend factory. For the 17 cases the port reproduces completely,
  every recorded progress call, every warning-level log line, every
  source-file field, the first twelve spectra and every chromatogram agree;
  the two relative-path cases agree against the working directory; the seven
  `computeFileHash` rows agree, the directory and the missing file included;
  the four divergences above are asserted as stated. A final test checks that
  every case in the fixture is exercised.

## Member ledger

| Source member | Rust | Status |
|---|---|---|
| `getType`, `getTypeByFileName`, `hasValidExtension`, `stripExtension`, `swapExtension`, `getConsistentOutputfileType`, `getTypeByContent` | `FileHandler::get_type`, `file_types::*`, `type_by_content` | ported ([MZML_MOBILITY_SUPPORT.md](MZML_MOBILITY_SUPPORT.md)) |
| `isSupported` | none; `can_read_experiment` and siblings report native availability | open |
| `getOptions`/`setOptions`, `getFeatOptions`/`setFeatOptions` | arguments (`ExperimentLoadOptions::peaks`, `load_feature_map_with_options`) | ported as arguments |
| `loadExperiment` | `load_experiment`, `load_experiment_with_options`, `_with_read_options`, `load_experiment_with_progress` | **partial**: DTA, DTA2D, mzML, mzXML, mzData, MGF, MS2 with every argument; SQMASS, XMASS, MSP, Bruker TDF, Thermo RAW open; imzML refused as in the source |
| `storeExperiment` | `store_experiment` | partial: DTA, DTA2D, MGF, mzML (MS2 a native extension); no log type; mzXML, mzData (with the source's chromatogram-to-spectra conversion), SQMASS, MSP open |
| `loadSpectrum`, `storeSpectrum` | none | open (DTA, XMASS) |
| `loadFeatures`, `storeFeatures` | `load_feature_map*`, `store_feature_map` | partial: featureXML only, no log type; TSV, PEPLIST, KROENIK, OMS, FEATUREPARQUET, EDTA open |
| `loadConsensusFeatures`, `storeConsensusFeatures` | `load_consensus_map`, `store_consensus_map` | partial: consensusXML only, no log type; EDTA, OMS, CONSENSUSPARQUET open |
| `loadIdentifications`, `storeIdentifications` | `load_identifications`, `store_identifications` | partial: idXML only, no log type; mzIdentML, OMS, xQuestXML, OMSSAXML, protXML, IDPARQUET open |
| `loadTransitions`, `storeTransitions` (TraML), `loadTransformations`, `storeTransformations` (trafoXML), `storeQC` (qcML, mzQC) | none | open |
| `computeFileHash` | `FileHandler::compute_file_hash` | **ported** |
| `applyPostLoadOptions_` (private) | none | only the TDF and RAW branches use it |

## Remaining gaps, ranked by the tools they block

1. **MGF through `mascot_generic`.** The source reads MGF with
   `MascotGenericFile`, whose progress section and carry-over semantics the
   port has (FileInfo already loads MGF that way); `FileHandler` still
   dispatches the stricter `format::mgf`. Switching closes MGF's log type; the
   Release rows for it are already in the fixture.
2. **The log type on the store and map entry points.** `storeExperiment`,
   `loadFeatures`/`storeFeatures`, `loadConsensusFeatures`/
   `storeConsensusFeatures` and the identification pair take a log type in the
   source; the port's readers and writers have `*_with_progress` for mzML,
   DTA2D, featureXML and consensusXML, but `FileHandler` does not wire them.
3. **mzXML and mzData output** through `store_experiment`, including the
   source's `convertChromatogramsToSpectra` before either write
   (`FileHandler.cpp:1184-1226`).
4. The identification, transition, transformation and QC formats the port has
   readers for (mzIdentML, TraML, trafoXML, qcML) but does not dispatch.
5. `loadSpectrum`/`storeSpectrum`, `isSupported`, and the formats with no
   native reader yet (SQMASS, XMASS, MSP, OMS, the parquet types, Bruker TDF,
   Thermo RAW), plus ZIP input.
