# mzML record settings and spectrum products

Existing `mzml` readers and writers now transport spectrum scan mode, polarity,
zoom and ordered scan windows, ordered spectrum Product lists, and the nine
represented chromatogram types. This applies to ordinary, zlib and Numpress
arrays and to the existing path/dispatch entry points. Public options and binary
precision policies are unchanged.

## Scientific representation

All fourteen non-Unknown ScanMode values use the exact source CV mapping.
The fileContent header summarizes present spectrum modes in source order, with
the source Mass fallback for Unknown or no spectra. Chromatograms do not
contribute to this source header summary, so chromatogram-only files also use
the Mass fallback.
Positive/negative polarity is a spectrum CV. Zoom writes inside the scan; reading
also accepts its deprecated spectrum location. Repeated zoom terms are harmless.
Scan modes, polarity, chromatogram types and individual window bounds reject
repeated definitions, including through referenceable groups. The source usually
keeps the last definition; the native correction detects conflicting metadata.

Unknown scan mode emits no record scan-mode term, retaining the existing native
roundtrip. The source writer instead forces a mass-spectrum CV. A missing
chromatogram-type CV reads as Mass. The writer emits the appropriate type CV for
all nine represented types; source alias `MS:1001474` reads as SRM and writes the
canonical `MS:1001473`. Unknown chromatogram type is a checked output error,
because omission would reread as Mass. Source non-mass signal types requiring
unrepresented coordinate/intensity roles (`MS:1003019`, `MS:1003020`,
`MS:1000626`) are checked input errors. This is not general non-MS detector-data
support, and the existing positive MS-level rule remains.

Scan windows preserve order and finite signed binary64 endpoints, including
equal endpoints and signed zero. Missing endpoints use the source zero default;
the native model rejects inverted windows. Windows across multiple scans append
in encounter order. Existing multiple-RT rejection remains unchanged. A scan can
be written for windows/zoom without manufacturing an RT when native RT is -1.
The source normally emits its RT fallback even in that case.

The two bound CVs use m/z units by default. A nondefault unit is stored exactly in
`ScanWindow.metadata["unit_accession"]`, as in source. Both endpoints must agree
on the unit; mixed units are a checked correction. MS/UO accessions with seven
numeric digits are preserved without conversion or ontology lookup. The writer
supplies known names for m/z and the source-tested nanometer, otherwise omits the
optional unitName. XML still carries unitAccession and unitCvRef. A separately
supplied `unit_accession` userParam is rejected instead of overwriting the bound
identity. A native explicit-default `MS:1000040` metadata key is rejected because
the source reader removes that key; callers represent the default by absence.

Window metadata shares the Product scalar codec: strings, i64 integers, finite
f64 and optional MS/UO units are retained. Empty/list values and units on the
unit-accession control value are checked output errors. InstrumentSettings'
separate metadata map has no source writer location and remains rejected.
Unknown scan-window CVs retain the existing reader's ignored-metadata policy;
this increment interprets only the two represented bounds.

Products use the existing isolation-window implementation for target and
asymmetric offsets, scalar metadata and units. The spectrum owner now stores an
ordered vector instead of rejecting it. The singular chromatogram Product is
unchanged. Both owners use temporary parsing state, strict list/count/nesting
checks, and the same reference-expansion path. The writer explicitly emits
negative-zero offsets to preserve their sign; positive-zero offsets retain the
existing default omission. No numerical rounding is introduced.

## Boundaries and failure behavior

AcquisitionInfo, SourceFile and DataProcessing attachments still fail the O(1)
loss guard. Nondefault InstrumentSettings on **chromatograms** also fail: the
pinned source writer and standard chromatogram XML grammar have no scan-list
location. Consequently a spectrum-to-chromatogram conversion can produce native
settings that cannot be written losslessly by this bounded adapter. There is no
private metadata projection to hide that limitation. Chromatogram-to-spectrum
conversion can carry source SRM settings, precursor and products through mzML,
then convert back. The tests exercise this useful direction and the remaining
chromatogram rejection explicitly.

Array-description loss guards remain unchanged. All new unsupported-state guards
run before deep experiment validation. Existing path writes remain atomic on
serialization errors. Existing anonymous native-ID generation is unchanged;
exact object roundtrips require stable IDs supplied before writing.

Empty Product/window objects consume the shared parameter count and byte limits;
declared counts are checked before growth. Inline and referenced settings use
the same accounting and validation. Filtered-out records still receive the
reader's full-payload validation, including these settings. The writer charges
newly supported settings against a fixed cumulative 50-million-work / 256-MiB
logical-payload preflight before validation/rendering. These additional limits
are separate from binary and Numpress budgets. No new resource-option fields are
introduced, and they are not process-memory guarantees.

## Evidence

[Eight settings tests](../tests/mzml_settings.rs) cover all source mode/type
mappings, literal window/Product values and metadata, unit identity, legacy zoom,
reference expansion and count/byte limits, malformed structures, excluded-record
validation, signed zero, converter transport, and independent XSD validation of
both writer paths. The adjacent Product/guard/loading/Numpress/path suites remain
part of focused validation. Tests distinguish source literals from native
roundtrips and arithmetic examples; no expected numerical values are generated
by Rust production code and no C++ execution was performed for this increment.

[Provenance](../tests/data/mzml_settings_provenance.json) pins
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. The relevant handlers, metadata records,
source fixtures and schema are byte-identical to the previous 54a pin. The
[projection generator](../tools/generate_mzml_settings_reference.py) reuses the
already packaged original `MzMLFile_1.mzML` by hash. It preserves four spectra's
settings/decimal strings and parameter-group definitions, removes unrelated
payload and sets peak counts to zero. **The original Product list declares one
but contains two Products; the projection explicitly corrects only this count to
two**, matching the source class assertions. Original bytes remain unchanged.
The separate source PDA window 220–500 nanometer case is tested as a settings-only
case; unsupported PDA arrays and source MS level zero are not claimed.

See [Product transport](MZML_PRODUCT_SUPPORT.md),
[remaining acquisition guards](MZML_ACQUISITION_GUARDS.md) and
[general mzML support](MZML_SUPPORT.md) for adjoining interfaces and limitations.
