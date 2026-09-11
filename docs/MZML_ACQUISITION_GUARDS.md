# mzML guards for newly owned acquisition fields

The native spectrum/chromatogram records now carry instrument settings,
acquisition information, source-file information and processing handles.
Spectra also carry Product lists; chromatograms carry a typed chromatogram kind.
The current mzML adapter does not serialize these additional fields. It rejects
nondefault attached state before output instead of silently discarding it.

[The shared writer preflight](../src/format/mzml.rs) checks the kernel's O(1)
`has_acquisition_settings()` predicates before calling general experiment
validation. This applies to ordinary and Numpress output, including zlib arrays.
It also checks for DataArray description metadata/processing before general
validation can traverse those descriptions. Numpress retains its separate shared
binary work/allocation preflight; no codec or encoding behavior changes here.

## What is rejected

For both record types, any of the following is unrepresented:

- Instrument scan mode other than Unknown, zoom enabled, polarity other than
  Unknown, a nonempty scan-window vector, or instrument metadata.
- A nonempty acquisition vector, combination-method string, or acquisition
  metadata.
- Any source-file name/path/type/checksum/native-ID string, nonzero file size,
  non-Unknown checksum type, CV terms, or source-file metadata.
- A nonempty processing-handle vector, even if every pointed-to record is default.

A nonempty spectrum Product vector is rejected, including a vector containing
one default Product. Any chromatogram type other than Mass is rejected. Empty
metadata values still count when their map contains a key. A nonempty vector
of default acquisition/window objects is different from an empty vector.

These decisions cover all source default fields, not generic container `is_empty`
shortcuts. In particular, the source AcquisitionInfo vector may be empty while
its combination method or metadata remains nondefault. SourceFile includes its
CVTermList metadata as well as its CV terms. The native check inspects only
scalars and container lengths; it does not allocate default records, clone
processing handles, or compare arbitrary nested payloads. Source floating equality equates negative-zero source-file size with default
zero. The native guard deliberately rejects negative zero as well, because
dropping that attached field would erase its IEEE sign; only positive zero is
accepted as the canonical default size.

## Existing supported fields remain supported

Existing spectrum precursor records, precursor activation/isolation/mobility
fields, scalar spectrum/chromatogram metadata, names, spectrum representation,
RT, native IDs, aligned arrays, and the singular chromatogram Product remain on
their existing paths. The new guard deliberately distinguishes spectrum Product
lists from the singular chromatogram Product already implemented by mzML.

Readers leave the newly attached fields at their native/source defaults under
the existing supported reader subset. Existing precursor-specific acquisition
transport does not populate the new spectrum/chromatogram InstrumentSettings,
AcquisitionInfo, SourceFile or processing fields. This increment does not claim
new acquisition-metadata parsing or serialization.

All currently relevant XML paths are covered: `mzml::write`,
`write_with_options`, `write_with_numpress`, ordinary mzML path `store` operations
and `FileHandler` mzML output dispatch share the same rejection. Path stores
continue to preserve an existing destination on serialization failure, including
gzip and bzip2 destinations. File-backed output may create and remove its owned
temporary file; it does not replace the destination on rejection.

idXML, FeatureXML and ConsensusXML operate on separate identification/map models,
not `MSSpectrum` or `MSChromatogram`; their existing processing serialization is
unchanged. There is no additional current mzXML/mzData raw-record writer to patch.
Text formats are covered separately by the kernel/conversion migration.

## Source and validation

The source is OpenMS4-core `54a232fe2cae9c590d5c997fa49d20e7769860fb`.
[The manifest](../tests/data/mzml_acquisition_guards_provenance.json) records the
relevant headers/implementations and exact default/equality locations. The source
MzMLHandler has explicit acquisition, source-file, processing and spectrum
Product-list output; the current native adapter implements a smaller transport
surface, so checked rejection is an intentional limitation.

[Seven focused tests](../tests/mzml_acquisition_guards.rs) independently exercise
all nondefault field categories, all nine nondefault chromatogram types, default
versus nonempty-default containers, metadata-only SourceFile state, CV-only state,
the stricter negative-zero preservation check, invalid nested data rejected before deep validation,
shared-handle ownership, ordinary/Numpress/dispatch output preservation, and
plain/gzip/bzip2 destination preservation. Existing Numpress, mzML, loading,
parameter-group and Product tests are rerun. No numerical expectations were
created from production output and no C++ build/execution was performed.
