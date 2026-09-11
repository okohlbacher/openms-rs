# Remaining mzML acquisition guards

The [record settings extension](MZML_SETTINGS_SUPPORT.md) now transports spectrum
scan modes, polarity, zoom, windows and Product lists, plus represented
chromatogram types. The original blanket guard has been narrowed only for these
source-backed paths. [Acquisition transport](MZML_ACQUISITION_SUPPORT.md) further
retains informative spectrum scan lists and scalar acquisition metadata.

Before general validation or output, the shared ordinary/Numpress writer checks
unsupported state using scalar fields and container lengths:

- Nonempty **chromatogram** AcquisitionInfo entries, combination method or metadata.
- Nondefault chromatogram SourceFile (mzML has no corresponding reference attribute).
- Source-file size and arbitrary CV payload unsupported by the source writer.
- Spectrum InstrumentSettings metadata, and all nondefault chromatogram
  InstrumentSettings (the standard/source chromatogram grammar lacks that slot).
- Unknown chromatogram type, whose omission would reread as Mass.
- Source-unsupported header values and nonscalar/unrepresentable metadata;
  represented source files, processing handles and auxiliary descriptions now
  use the [header/reference codec](MZML_HEADER_SUPPORT.md).

Empty metadata values count if the map has a key. An empty AcquisitionInfo vector
can still carry unsupported method/metadata on chromatograms. Cheap guards do not allocate defaults or traverse nested payload. The separate
header preflight validates and meters supported shared processing payload.
Supported spectrum settings, acquisition records and Product lists then receive a bounded cumulative
preflight and ordinary value validation before XML begins. Unknown combination
methods, unrepresentable instrument-reference metadata and nonscalar acquisition
values fail during that bounded preflight.

Existing precursor isolation/activation/mobility, singular chromatogram Product,
record names/metadata, aligned arrays and binary precision remain on their
established paths. SourceFile and processing attachments now resolve through the header registry;
spectrum AcquisitionInfo follows its explicit Canonical/Source mode. Precursor-specific acquisition data
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
