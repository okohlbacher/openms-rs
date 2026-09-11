# Remaining mzML acquisition guards

The [record settings extension](MZML_SETTINGS_SUPPORT.md) now transports spectrum
scan modes, polarity, zoom, windows and Product lists, plus represented
chromatogram types. The original blanket guard has been narrowed only for these
source-backed paths. [Acquisition transport](MZML_ACQUISITION_SUPPORT.md) further
retains informative spectrum scan lists and scalar acquisition metadata.

Before general validation or output, the shared ordinary/Numpress writer checks
unsupported state using scalar fields and container lengths:

- Nonempty **chromatogram** AcquisitionInfo entries, combination method or metadata.
- Any SourceFile strings, noncanonical file size, checksum type, CV terms or
  metadata. Only positive zero is the canonical default size; negative zero
  would otherwise silently lose its sign.
- Any DataProcessing handles, including a default pointed-to record.
- Spectrum InstrumentSettings metadata, and all nondefault chromatogram
  InstrumentSettings (the standard/source chromatogram grammar lacks that slot).
- Unknown chromatogram type, whose omission would reread as Mass.
- Auxiliary array description metadata or processing handles.

Empty metadata values count if the map has a key. An empty AcquisitionInfo vector
can still carry unsupported method/metadata on chromatograms. Checks do not allocate defaults,
clone handles, compare nested values or traverse shared processing payload.
Supported spectrum settings, acquisition records and Product lists then receive a bounded cumulative
preflight and ordinary value validation before XML begins. Unknown combination
methods, unrepresentable instrument-reference metadata and nonscalar acquisition
values fail during that bounded preflight.

Existing precursor isolation/activation/mobility, singular chromatogram Product,
record names/metadata, aligned arrays and binary precision remain on their
established paths. SourceFile and processing attachments remain at defaults in the supported reader
subset; spectrum AcquisitionInfo now follows its explicit Canonical/Source mode. Precursor-specific acquisition data
does not populate those separate attachments.

`mzml::write`, `write_with_options`, `write_with_numpress`, mzML path stores and
FileHandler dispatch share the guard. Existing plain/gzip/bzip2 destinations
remain unchanged on rejection. idXML, FeatureXML and ConsensusXML use separate
map/identification models and are unaffected by these raw-record fields.

[Seven guard tests](../tests/mzml_acquisition_guards.rs) cover the remaining
categories, early rejection before invalid nested/shared values, array
descriptions, negative-zero SourceFile size and destination preservation.
[Settings tests](../tests/mzml_settings.rs) cover the newly accepted categories.
The original [54a provenance](../tests/data/mzml_acquisition_guards_provenance.json)
is historical evidence for the blanket-guard increment; the
[82ce extension manifest](../tests/data/mzml_settings_provenance.json) records its
narrowing and unchanged source defaults. No C++ execution is claimed.
