# Native featureXML support

`format::featurexml` reads and writes featureXML 1.9 and the checked legacy forms
used by the pinned Core SDK. `FeatureFileOptions` is a passive owned options
record. Free functions provide `read`, `read_with_options`, `read_with_registry`,
`read_into`, `read_size`, `write`, `write_with_options`, `write_with_registry`,
`load`, `load_with_options`, `load_into`, `load_size`, `store`, and
`store_with_options`. The Cargo feature is `featurexml`.

## Fields and source conventions

The adapter transports document and unique IDs, typed map metadata, processing
history, identification runs and search parameters, protein hits, assigned and
unassigned peptide IDs, and recursive features. Each feature retains RT, m/z,
f32 intensity, overall and per-dimension qualities, signed charge, unique ID,
convex hulls, subordinates, and typed metadata. Hulls preserve the ordered XML
outline; both modern `pt` and legacy `hullpoint/hposition` syntax are accepted.
Writing compresses a copy of scan-backed hulls as the source does. Reading an
outline does not invent internal scan data.

Map `id` and obsolete `unique_id` are supported; the latter wins if both occur.
Nested feature IDs use the last underscore-delimited decimal component. Invalid
non-numeric IDs become unassigned zero; numeric overflow is rejected. Assigned
top-level IDs must be unique when writing. Old `description` and removed `model`
payloads are ignored, as in the source. The legacy `userParam` spelling is accepted.

`load` stores its path and `FileType::FeatureXml` on the returned map; stream
reading leaves loaded-file identity empty/unknown. Loaded-file identity is not
itself serialized. Path operations use the common plain/gzip/bzip2 transport;
compressed input is detected by content and output compression by suffix. Stores
validate the extension and replace the destination atomically.

Metadata retains scalar and list integer, floating, and string types. Width has
no XML element: `BaseFeature::set_width` stores matching numeric `FWHM` metadata.
Reading restores width from this metadata on top-level features only. Subordinate
FWHM metadata remains available, but its width field stays zero. Writing rejects
nonzero subordinate width and inconsistent top-level width/FWHM, preventing an
unreported change on the next read. A NaN width and a NaN `FWHM` agree, as
`setWidth` leaves them.

### Non-finite values

The source writes and reads them, so this dialect does too. `writeFeature_`
sends every scalar through `precisionWrapper` and `writeUserParam_` sends a
`DataValue` through the same conversion, and
`Internal::NumericFormatting::appendNumeric` answers before it formats
anything: `NaN` for a NaN of **either** sign, then `-inf` or `inf`
(`src/common/include/OpenMS/CONCEPT/Detail/NumericFormatting.h:29-35`; the
header records that `inf.0` was once produced and could not be read back).
`StringUtils::toDouble` reads all three: `tryParseNaN` takes any case of `nan`,
and `std::from_chars` takes any case of `inf`, `infinity` and `nan` with an
optional sign, after a leading `+` is skipped
(`StringUtils.cpp:35-58, 239-276`). The executed evidence is
`../oracle/featurexml-inf`, whose probes store every non-finite value in every
place featureXML has and read it back, and take twelve spellings through the
pinned reader.

So the writer writes `inf`, `-inf` and `NaN` for a non-finite position,
intensity, quality, overall quality, `float` meta value or `floatList` entry,
and the reader accepts `inf`, `+inf`, `-inf`, `infinity`, `Infinity`, `INF`,
`NaN`, `nan`, `NAN` and `-nan`. `FeatureFinderAlgorithmPicked` reaches infinite
widths, `FWHM` values and intensities on finite input, and the port's
`FeatureFinderCentroided` now writes what the Release build writes (this closed
TOPP native difference 16). `MetaValue::validate` and the public
`TryFrom<f64>` still refuse a non-finite value; the reader stores one through
the crate-private `MetaValue::source_float` and `source_float_list`, as the
algorithm does (lead decision D13).

Three refusals remain, all of them older than and separate from the spellings:

- **A non-finite hull point.** `ConvexHull2D::setHullPoints` validates nothing
  (`ConvexHull2D.cpp:119-123`), so the Release build writes and reads one; this
  port refuses it, because its hulls and the bounding boxes derived from them
  rest on finite coordinates — a scan ordering and `ConvexHull2D::bounding_box`
  are undefined for a NaN. No ported algorithm produces such a point: the hulls
  `FeatureFinderAlgorithmPicked` builds come from peak positions and stay finite
  at every retention-time scale measured, including the ones where its widths
  and intensities overflow. The refusal is a parse error, never a panic.
- **A finite value that `f32` cannot hold**, such as `1e50` in an `intensity`.
  The source narrows it to an infinity; this port refuses it.
- **A decimal literal that is merely out of range**, such as `1e999`. The
  pinned reader fails the load on one too — `attributeAsDouble_` lets the
  `ConversionError` out (`XMLHandler.h:401-406`), recorded as `LOAD_THREW` in
  `../oracle/featurexml-inf/results/spellings.tsv` — while an element's text
  would reach `asDouble_`, which logs a non-fatal line and keeps `0.0`. Only a
  spelled token is non-finite here.

A finite **negative** width also stays refused. `setWidth` stores any value
(`BaseFeature.cpp:87-94`), but no source path produces a negative FWHM, and the
mirror check would have to let it through silently.

Accepting such a document is what lets one reach the rest of the crate, and
the kernel's own finite invariants then apply. `FeatureMap::ranges` — the
port's `FeatureMap::updateRanges` — refuses a non-finite value, so the map
reads and its ranges are a checked error. The Release `FileInfo` on the
document above exits 0 and prints `retention time: -inf .. inf sec (inf min)`
and `Total ion current in features: nan`
(`../oracle/featurexml-inf/results/fileinfo_nonfinite.out`); this port's
`FileInfo` exits 6 with `Invalid parameter: invalid value: range value must be
finite`. That is a **checked refusal, never a panic**
(`a_nonfinite_map_reads_and_its_ranges_are_a_checked_error`); making the ranges
themselves non-finite is a kernel change with its own evidence, not a
featureXML one, and is open for the lead.

## Options

The defaults load convex hulls and subordinates, load all metadata and feature
payload, and apply no ranges. Setting `load_convex_hulls` or `load_subordinates`
false skips scientific interpretation of those subtrees. XML well-formedness
and document limits still apply. RT, m/z and intensity ranges include the lower
endpoint and exclude the upper endpoint. Empty/inverted ranges select nothing;
nonfinite bounds fail when used. Every feature is filtered independently, so a
rejected parent removes its entire subtree. These options never filter writing.

Metadata-only reading stops interpreting XML at the opening `featureList` and
returns the preceding map/identification data. `read_size` / `load_size` return
that element's declared count, without interpreting its feature payload. The
source's separate `FeatureFileOptions.size_only` flag is retained but does not
change normal reading. Metadata-only takes precedence over count-only reading,
returning count zero without interpreting the count attribute. The bounded prefix reader stops before
feature payload bytes are consumed or decoded; comments, quoted attributes and
UTF-16 input are handled before deciding that the root's featureList has begun.

## Identification transport and checked differences

The shared map/identification XML codecs handle protein references, peptide
positions and flanks, score metadata, native rank/run-identifier extensions, and
portable named modification definitions. Registries are owned locally; reading
a document does not mutate the global modification database. Custom sequence
chemistry must be representable by its portable definition or resolvable from
the explicitly supplied registry. Malformed definitions, chemical collisions,
and unsupported chemistry fail before caller state is changed.

The following native behavior is deliberate:

- Metadata following a subordinate peptide identification stays on the actual
  owning subordinate. The source resets its metadata pointer to the top-level
  feature and can misplace these values.
- Orphan peptide run identifiers and unknown protein accessions are errors. The
  source can warn and omit an orphan, or silently bind an accession to `PH_0`.
- Structured protein groups are rejected by this dialect. Source featureXML
  transports raw UserParams but does not encode its in-memory group structures;
  use consensusXML for group quantity arrays. Raw group-like UserParams are
  preserved without an invented group interpretation.
- Out-of-range f32 values, non-finite hull points, invalid dimensions,
  unsupported XML fields, malformed references, and unrepresentable native
  metadata are checked errors (see *Non-finite values* above for what is
  accepted instead). Required timestamps must be present rather than
  emitting an invalid schema timestamp.
- Metadata on a data-processing entry, an identification run, a protein hit or
  a peptide identification inside a featureXML still refuses a non-finite
  value: those are written through the shared map/identification codecs, which
  every dialect uses, and no source algorithm puts a non-finite value there.
  Only the map's own metadata and a feature's carry the source's spellings.

Modern IdentificationData graph attachments are not serialized by either source
map XML handler and are outside this transport. A complete public inherited
`XMLFile::isValid` schema-validation API is not claimed here. Schema checks of
written fixtures are separate validation evidence.

## Bounds and verification

`ReadOptions` contains `feature_options`, `Limits` and `InputScaling`;
`WriteOptions` contains `Limits` and `OutputScaling`. `Limits` holds the
absolute ceilings and the scalings hold their growth with the size of the work:
the effective ceiling is the smaller of the two, and every absolute ceiling but
`max_xml_bytes` (8 GiB) and `max_depth` (128 subordinate levels) defaults to
unbounded so the size-derived allowance decides. The reader's allowances are
`floor + rate * decoded bytes` with the former fixed ceilings as floors — one
million elements and list items, 50 million work units, 256 MiB of payload — so
a small document is bounded exactly as before while a real one is admitted. The
writer's grow with the counted size of the map instead: features, hull points,
identifications, identification hits and metadata entries. `Limits::former()`
with `InputScaling::fixed()` or `OutputScaling::fixed()` restores the former
fixed behaviour.

Work and payload are charged before geometry/metadata copies, registry copies
and identification conversion. Parsing, registry work, and conversion use shared
counters; limits apply to the whole operation, including all subordinates.
Native allocation and output checks precede materialization; an allocation the
host cannot satisfy is a checked error rather than an abort.

Features are converted as the parser finishes each one rather than after the
whole document is in memory, so the parse tree is the size of one feature. The
`featureMap` header is interpreted from the children that precede `featureList`,
which FeatureXML_1_9.xsd places last; a root child after `featureList` is
therefore refused (`featureMap content after featureList`) rather than read too
late to inform any feature.

`docs/FEATUREXML_SCALE_SUPPORT.md` records the ceiling that was there, the
measured charges the rates are derived from, and wall time and peak memory
against the C++ `FileInfo` on the 59.6 MiB and 2.06 GiB benchmark featureXML
files.

Reads produce an owned draft; `read_into` and `load_into` replace their target
only after success. Writers prepare and validate the complete output before
writing to a caller stream. A stream's own I/O failure can still leave bytes in
that external stream; atomic file replacement is provided by `store`.

[Focused tests](../tests/featurexml.rs) use three byte-identical upstream fixtures
and independent boundary cases. Two further fixtures are documents the Release
build itself wrote: `featurexml_nonfinite_release.featureXML` and
`featurexml_nonfinite_hull_release.featureXML`, which differ only in whether
their outline points are non-finite. [Provenance](../tests/data/featurexml_provenance.json)
pins source files, schemas and fixture SHA-256 hashes. The executed C++ `FileInfo` is used
as the scale reference in `docs/FEATUREXML_SCALE_SUPPORT.md`, which the
`#[ignore]`d `hpc_scale_benchmark_featurexml_files_load_and_round_trip` checks
against by reading the benchmark files by path, and the executed
`FeatureXMLFile` of `../oracle/featurexml-inf` is the oracle for the non-finite
spellings above. Everything else is read from the pinned source.

The source file loader also accepts ZIP archives. ZIP input is not yet supported
by the native shared path transport and remains an explicit completion gap.
