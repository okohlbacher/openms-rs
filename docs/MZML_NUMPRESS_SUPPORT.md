# mzML Numpress binary transport

This increment implements the linear, positive-integer (PIC), and short logged
float (SLOF) mzML binary transports against OpenMS4-core
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. The existing reader, filtered reader,
and their file-path callers accept these arrays automatically. A separate
writer operation provides the three source Numpress configurations without
changing `WriteOptions` or the behavior of ordinary writers.

## Public API

[The implementation](../src/format/mzml_numpress.rs) exports through `format::mzml`:

- `NumpressWriteOptions { binary: WriteOptions, mass_time: NumpressConfig,
  intensity: NumpressConfig, float_data_array: NumpressConfig,
  limits: NumpressCoderLimits }`.
- `write_with_numpress(writer, experiment, &options) -> Result<NumpressWriteReport>`.
- `NumpressWriteReport { encoded_arrays, fallback_arrays, ordinary_arrays }`.

The configurations default to disabled Numpress, as source `PeakFileOptions`
does. `binary.zlib_compression` applies to every array. Successful Numpress
arrays declare 64-bit floating data. Coordinates and chromatogram times use
`mass_time`, intensities use `intensity`, and named/canonical auxiliary floating
arrays use `float_data_array`. Integer and string arrays remain ordinary binary.

The writer first performs the existing whole-document metadata/value validation
and prepares every encoded array. Empty arrays, raw codec failures and failed
source error checks fall back to ordinary encoding. Original f64 coordinates and
f32 intensity/auxiliary values are preserved on fallback. Resource and
transport errors return before any XML is written; actual writer I/O failures
can leave partial output. Counts distinguish requested-but-fallback arrays from
those for which Numpress was disabled. They do not expose individual codec
rejection messages; use the standalone coder to inspect those in detail.

Default writer options produce the same ordinary bytes as `write_with_options`
for the tested fixtures, with and without zlib. This is not a promise of identical
compressed bytes across different zlib implementations or versions.

## Source behavior and checked boundaries

The reader accepts `MS:1002312`, `MS:1002313`, `MS:1002314` and the combined zlib
accessions `MS:1002746`, `MS:1002747`, `MS:1002748`. A standalone Numpress term and
one separate zlib term may appear in either order. Inline and referenced
parameters execute the same validation. Duplicate codec/compression declarations,
multiple modes, or `no compression` combined with Numpress are errors; source
CV handling instead overwrites fields in encounter order.

The source helper repairs missing Numpress precision to float64, promotes declared
float32 Numpress to float64, and repairs PIC integer declarations to float64.
These three repairs are retained. Other integer or string Numpress combinations
are rejected. Repairs apply only to Numpress. Canonical array identity/type checks
use the effective float64 type: for example, a charge array cannot silently
become a floating annotation. A canonical float32-only array (`mean charge array`)
is emitted as ordinary float32 when Numpress is requested, so its declared type
remains compatible; this is an explicit native correction to source emission.

Every decoded array must match its declared/default point count, including empty
auxiliary placeholders. The source Numpress branch omits the ordinary decoder's
length warning/repair. Invalid casts, nonfinite decoded values, malformed/truncated
nibbles, invalid base64 and truncated/corrupt/trailing zlib frames are checked
errors. Empty uncompressed Numpress text is accepted only when the count is zero;
empty compressed text is rejected. Nonempty text shorter than four characters is
rejected even though the standalone source-compatible wrapper ignores it.

Time unit scaling happens after decoding, before selection. Numpress values stay
f64 through the scientific loader's range checks, sorting and aligned selection;
only retained intensity and auxiliary values narrow to f32. Ordinary readers
narrow all represented values. Excluded whole records still undergo full native
decoding/validation, including f32 representability, as documented by the existing
loader. This does not claim source early-skip behavior.

The new entry point keeps the existing fixed ordinary precision choices; it does
not implement general `PeakFileOptions` writer migration. In particular, it does
not copy the source writer's cross-setting interaction in which enabling
mass/time Numpress may also select f64 ordinary intensity fallback. Independent
sampled noise grids, additional primary detector roles, indexed writing and other
unrepresented model fields remain separate work. [Header/array descriptions and
metadata-only headers](MZML_HEADER_SUPPORT.md), plus [consumer streaming and
disabled record population](MZML_CONSUMER_SUPPORT.md), are now implemented.
No new path writer overload is
introduced. The `numpress` feature is already enabled by `mzml`; no new dependency
or feature is required by this increment.

## Resource accounting

Readers preserve `ReadOptions` raw declared point, total array-element, per-array
byte and total decoded-byte caps. Numpress decoded storage is counted as eight
bytes per element, regardless of the precision label that needed repair. Compressed
and decompressed codec bytes are bounded by `max_array_bytes`. Additionally, one
shared coder session spans all Numpress arrays, with the standalone defaults of
500 million weighted work visits and 512 MiB cumulative logical allocations
(including temporary/retained decode buffers and a conservative zlib-state charge).
These additional session caps are currently fixed; they do not reset per record.
The ordinary reader's established resource contract is unchanged.

The new writer bounds counts and precharges binary validation work before the
general native validator visits any peak or auxiliary scalar. Array description
metadata is rejected through an O(1) check before deep validation.

The new writer uses its supplied coder limits cumulatively across all arrays,
including promotion, estimation, codec verification, fallback bytes, zlib,
base64 and retained encoded-array slots/text. Per-array element/binary/text limits
are checked too. These are binary transport bounds, not a new whole-document
metadata allocation budget. Existing unsupported metadata checks precede output.

## Evidence

[Tests](../tests/mzml_numpress.rs) cover all six accessions, separate zlib,
precision repairs, strict conflicts/counts, canonical identity, all three writer
configs, source fallback, referenced parameters, minutes, scientific filtering,
all aligned array kinds, shared budgets, pre-output failures and independent XSD
validation when `xmllint` is available. A private test checks the same decode
session exhausting work and allocation budgets across successive arrays.

[Provenance](../tests/data/mzml_numpress_provenance.json) hashes the source and
fixtures. The three exact class-test base64 strings are reused from the existing
[raw byte fixture](../tests/data/numpress_source_bytes.tsv); their zlib projections
are reused from the existing [wrapper fixture](../tests/data/numpress_coder_transport.tsv).

The [original upstream mzML](../tests/data/mzml_numpress_source_original.mzML) is
copied without modification. Its [projection](../tests/data/mzml_numpress_source_projection.mzML)
keeps all 18 chromatograms, 36 original binary array subtrees and their exact
encoded payloads. XML is reserialized with a UTF-8 declaration in place of the
original ISO-8859-1 declaration, with unrelated metadata/index sections omitted. The
[342-point table](../tests/data/mzml_numpress_source_values.tsv) is independently
derived by Python standard-library base64/zlib plus signed-nibble linear and SLOF
arithmetic. It is not a C++-executed golden or a Rust-generated expectation.
[The regeneration script](../tools/generate_mzml_numpress_reference.py) records the
transformation. No new C++ execution was performed for this transport increment;
the raw codec's prior differential evidence remains separately documented.
