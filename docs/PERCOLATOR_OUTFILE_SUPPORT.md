# PercolatorOutfile support

`FORMAT/PercolatorOutfile.h` and `FORMAT/PercolatorOutfile.cpp` at `bc9cc12`,
ported to `src/format/percolator_outfile.rs`, together with the
regular-expression path of `METADATA/SpectrumLookup.h` and
`METADATA/SpectrumMetaDataLookup.h`, ported to
`src/metadata/spectrum_meta_data_lookup.rs` (group F's separate port of the regular-expression path was merged into it at the wave-4 integration). Tests are in
`tests/percolator_outfile.rs`; fixtures, hashes and source anchors are in
`tests/data/percolator_outfile_provenance.json`.

A Percolator PSM file (`.psms`) is one tab-separated header line,
`PSMId score q-value posterior_error_prob peptide proteinIds`, then one row per
PSM. The reader needs no Cargo feature; the test file is gated on `idxml`
because its differential reads the retained C++ outputs with the idXML reader.

## API mapping: `PercolatorOutfile.h`

Every public and private member of the header.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `enum class ScoreType { QVALUE, POSTERRPROB, SCORE, SIZE_OF_SCORETYPE }` | `ScoreType` | `SIZE_OF_SCORETYPE` is a count sentinel, not a value; `SCORE_TYPE_NAMES.len()` plays its role. |
| `static const std::string score_type_names[]` | `SCORE_TYPE_NAMES`, `ScoreType::name` | `["qvalue", "PEP", "score"]`, the `IDFileConverter -score_type` choices. |
| `static ScoreType getScoreType(std::string)` | `ScoreType::from_name`, `FromStr` | ASCII case folding as `StringUtils::toLower` in the C locale; `InvalidValue` for any other name, as the source. |
| `PercolatorOutfile()` | `PercolatorOutfile::new`, `Default` | Stateless. |
| `void load(filename, proteins, peptides, lookup, output_score = QVALUE)` | `PercolatorOutfile::load`, `PercolatorOutfile::read` | Returns a `PercolatorDocument` (the protein identification, the peptide identifications and the source's log lines) instead of filling out-parameters. `read` takes a stream. |
| `void getPeptideSequence_(std::string, AASequence&) const` (private) | `Rewrites::peptide` (private) | Flanking residues, `[unknown]`, `[UNIMOD:n]`, the missing `+` of a positive mass shift. |
| `void resolveMisassignedNTermMods_(std::string&) const` (private) | `Rewrites::resolve_misassigned_n_term_mods` (private) | The X! Tandem N-terminal shift rule, including the carbamidomethyl-on-C exception. |
| (none) | `HEADER`, `DEFAULT_REFERENCE_FORMATS`, `PercolatorDocument` | The header and the three reference formats the source spells inline in `load`. |

## Preserved source conventions

- **Default reference formats.** When the lookup has no reference format,
  `load` registers MS-GF+, Mascot and X! Tandem formats in that order, and the
  lookup keeps them afterwards, even when the load fails. The Mascot expression
  writes `[(scans:)(scan=)(spectrum=)]`, which Boost reads as a character class;
  it is kept verbatim because that is what the source matches.
- **An unresolvable reference is logged, not fatal**, and the hit keeps
  whatever the reference itself supplied before the look-up failed.
- **Scores.** All three scores are stored as `Percolator_score`,
  `Percolator_qvalue` and `Percolator_PEP`; `output_score` selects the hit's
  score and the identification's score type (`q-value`, `Posterior Error
  Probability`, both lower-is-better, or `Percolator_score`, higher-is-better).
- **Proteins.** One `ProteinHit` per distinct accession, in byte order (the
  source's `std::set<std::string>`); search engine `Percolator`, identifier
  `id`, the load time as the run date; fixed and variable modifications are
  inferred from the peptides with `ModificationDefinitionsSet::inferFromPeptides`.
- **Log text.** One error per unresolved reference, one warning per peptide
  whose `[unknown]` modifications were removed, the two "Created ..." counts,
  and one warning per kind of missing charge, retention time or m/z, returned in
  `PercolatorDocument::log` in the source's order.
- **`regex_replace` semantics.** The three rewrites replace every
  non-overlapping match left to right, as `boost::regex_replace` with default
  flags does, through `crate::concept::boost_regex`.

## Divergences (all refusals of undefined or unrepresentable input)

- A row with fewer than five fields is refused; the source reads `items[1]` to
  `items[4]` unchecked.
- A peptide shorter than two characters is refused; the source reads
  `peptide[len - 2]`.
- A non-finite score is refused: a hit's score and its meta values must be
  finite in this crate.
- An empty file is refused with the header message; the source's `getRow(0)`
  throws on it.
- A peptide whose first byte is not ASCII cannot take the misassigned-shift
  path and is refused with `Unsupported`.

## API mapping: `SpectrumLookup.h` and `SpectrumMetaDataLookup.h`

`SpectrumMetaDataLookup` folds the `SpectrumLookup` base into one type.
`mascot_xml::SpectrumTitleLookup` and `pepxml::SpectrumIndex` predate it and
still serve their own consumers.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `default_scan_regexp` | `DEFAULT_SCAN_REGEXP` | |
| `regexp_names_` (protected) | `REFERENCE_GROUP_NAMES` | |
| `reference_formats` | `reference_formats()`, `add_reference_format` | Read-only view; formats are added through the checked function. |
| `rt_tolerance` | `rt_tolerance` field | 0.01 s default. |
| `empty()` | `is_empty`, `len` | |
| `readSpectra(spectra, scan_regexp, get_precursor_rt)` | `read_spectra` | Replaces earlier contents; the source appends to its meta data and keeps stale index entries. Returns the log lines. |
| `findByRT`, `findByNativeID`, `findByIndex`, `findByScanNumber`, `findByReference` | `find_by_rt`, `find_by_native_id`, `find_by_index`, `find_by_scan_number`, `find_by_reference` | `ElementNotFound`/`IndexOverflow` become `InvalidRange`, `ParseError` becomes `InvalidValue`. |
| `addReferenceFormat` | `add_reference_format` | Nothing is registered on error. |
| `MetaDataFlags`, `MDF_*` | `MetaDataFlags` with `RT`, `PRECURSOR_RT`, `PRECURSOR_MZ`, `PRECURSOR_CHARGE`, `MS_LEVEL`, `SCAN_NUMBER`, `NATIVE_ID`, `ALL` | Bit values as the source. |
| `SpectrumMetaData` | `SpectrumMetaData` | Quiet NaN becomes `None`; the `-1` scan sentinel becomes `None`. |
| `getSpectrumMetaData(Size, SpectrumMetaData&)` | `meta_data` | |
| `getSpectrumMetaData(spectrum_ref, meta, flags)` | `meta_data_for_reference` | |

Not ported yet (the METADATA headers stay `partial`): the static
`getSpectrumMetaData(const MSSpectrum&, ...)` overload, `setSpectraDataRef`,
`addMissingRTsToPeptideIDs`, `addMissingIMToPeptideIDs`,
`addMissingFAIMSToPeptideIDs`, `addMissingSpectrumReferences`, the deprecated
static forwarders `extractScanNumber` (two overloads), `getRegExFromNativeID`
and `isNativeID`, and a transcription of `SpectrumLookup_test.cpp` and
`SpectrumMetaDataLookup_test.cpp`.

## Evidence

- **Tier 1, retained upstream output.** `TOPP_IDFileConverter_12`, `_13` and
  `_14` (`topp/CMakeLists.txt:1235-1245` at test-data `0cb15f2`) convert
  `IDFileConverter_12_input.psms` with `-score_type qvalue`, `PEP` and
  `score`. The retained outputs are read back with `crate::format::idxml` and
  compared record by record; the run date is excluded, as the upstream
  `-whitelist "IdentificationRun date"` excludes it.
- **Tier 3, class test.** The `getScoreType` and `load` sections of
  `PercolatorOutfile_test.cpp`; its other two sections only construct and
  destroy the reader.
- **Tier 3/4.** Reference formats, look-ups, peptide rewriting and refusals,
  derived from reading the source.
