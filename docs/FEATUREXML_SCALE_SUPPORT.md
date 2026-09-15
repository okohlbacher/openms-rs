# featureXML reader and writer: size-derived ceilings

The OpenMS4 TOPP benchmark (results directory
`/ceph/ibmi/abi/oliver/bench/openms4/results/2026-09-15-w3threads`) ran
`FileInfo` on the 59.6 MiB benchmark featureXML at one and at 32 threads. All
twenty Rust executions failed, exit 3:

```
Unable to read file (parse error on line 0: identification XML byte limit exceeded)
```

The C++ `FileInfo` summarised the same file in about 0.8 s. The benchmark's
adversarial review bisected the port's real ceiling by truncating a
`FeatureFinderCentroided` output: it read 1,672,637 bytes and refused 3,196,930,
i.e. "the port's featureXML reader works to ~1.7 MB". The benchmark's second
featureXML input, 2.06 GiB with 826,019 features, was never attempted.

## The ceiling that was there, exactly

`Limits::default()` was 64 MiB of XML, one million elements and list items, 128
subordinate levels, 50,000,000 work units and 256 MiB of payload. Those three
cumulative ceilings were coupled inside the shared parser
(`identification_xml::parse_xml_with_budget`):

```rust
let decode_limit = max_bytes.min(*remaining_bytes / 8).min(*remaining_work / 4);
```

so the reader would not even *decode* more than

```
min(64 MiB, 256 MiB / 8, 50,000,000 / 4) = 12,500,000 bytes
```

and a document above that failed with `identification XML byte limit exceeded`.
Below it the same 50,000,000 work units were charged four per decoded byte
before the tree was built, so the last stretch under 12,500,000 was already
unreachable: a 12,500,000-byte document exhausted the work ceiling instead. The
review's 1.7 MB / 3.2 MB bracket is the payload ceiling biting first on
hull-dense content — the same coupling, a different term.

`tests/featurexml.rs::the_former_fixed_ceilings_still_refuse_exactly_what_they_refused`
pins both boundaries, `Limits::former()` and `InputScaling::fixed()` reproducing
the former ceilings exactly.

## What the source does

`FeatureXMLFile::load` (`src/openms/source/FORMAT/FeatureXMLFile.cpp:44-68`)
hands the file to `Internal::FeatureXMLHandler`, a Xerces SAX handler
(`FeatureXMLHandler.h:42-45`), through `XMLFile::parse_`
(`XMLFile.cpp:116`). There is no resource ceiling anywhere on that path: the
handler appends features until the document ends. Its only nod to size is a
reservation hint, not a refusal:

```cpp
map_->reserve(std::min(Size(1e5), count)); // reserve vector for faster push_back, but with
                                           // upper boundary of 1e5 (as >1e5 is most likely
                                           // an invalid feature count)
```
(`FeatureXMLHandler.cpp:318`) — and the executed C++ reads the 826,019-feature
map without complaint, so the comment is about the reservation, not about what
the reader accepts.

The port keeps bounded work as a principle rather than copying the absence of
ceilings, so the fix is to derive them from the input instead of fixing them.

## What changed

### 1. Size-derived allowances (`src/format/featurexml_scaling.rs`)

Modelled on `src/format/mzml_scaling.rs`, the established pattern. An
`Allowance` is `floor + units * (consumed / per_bytes)`, saturating.
`InputScaling` carries one per cumulative reader quantity, `OutputScaling` one
per writer quantity. `Allowance` is repeated rather than imported from
`mzml_scaling` because that module is compiled only with the `mzml` feature,
which `featurexml` does not imply.

The absolute ceilings in `Limits` stay, and the effective ceiling is the smaller
of the two, but their defaults become unbounded so the allowances decide:

| `Limits` field | former default | new default |
|---|---|---|
| `max_xml_bytes` | 64 MiB | 8 GiB |
| `max_records` | 1,000,000 | unbounded |
| `max_list_items` | 1,000,000 | unbounded |
| `max_depth` | 128 | 128 |
| `max_work` | 50,000,000 | unbounded |
| `max_payload_bytes` | 256 MiB | unbounded |

`max_xml_bytes` is the one quantity nothing can be derived from — on input it
bounds the decoded document, on output the rendered one — so it keeps a large
but finite default and stays the only absolute ceiling that does. `read_size`
and metadata-only reading stop at the opening `featureList`, so neither is
charged for the payload it never decodes.

The reader's allowance floors are the former fixed ceilings, so a small document
is bounded exactly as tightly as it was; the rates are:

| quantity | default | measured need |
|---|---|---|
| `records` (XML elements) | 1,000,000 + 1 per 8 B | densest input opens one per 42.5 B |
| `list_items` | 1,000,000 + 1 per 1,024 B | a list is a property of one value |
| `work` | 50,000,000 + 64/B | at most 8.08/B |
| `payload_bytes` | 256 MiB + 32/B | at most 16.9/B |

The `work` and `payload_bytes` charges are deliberately conservative multiples
of what is really spent and held — the parser charges eight payload units per
decoded byte for a buffer that is one byte per byte — so these are ceilings on a
proxy, not byte counts. Real peak memory is measured below. The payload rate is
deliberately the tightest of the four: because a streamed subtree's charge is
refunded, the payload counter measures what is *retained*, and 32 per byte sits
four times above the 8.2 a streamed featureXML charges and below the roughly 34
a whole retained tree costs.

Measured charges, `FileInfo` reading each benchmark input:

| input | decoded bytes | work charged | per byte | payload charged | per byte |
|---|---|---|---|---|---|
| sanity `LCMS-centroided.featureXML` | 193,180 | 1,560,064 | 8.08 | 3,259,742 | 16.87 |
| `UPS1_50amol_R1.featureXML` | 62,504,565 | 418,465,214 | 6.70 | 510,521,698 | 8.17 |
| `UPS1_500amol_R3.featureXML` | 2,207,668,458 | 15,724,807,079 | 7.12 | 17,831,447,048 | 8.08 |

The small file's higher ratios are its fixed header cost spread over few bytes,
which is what the floors exist for; the rates are set from the large ones, where
the ratio has settled.

The writer has no input to measure, so its allowances grow with the counted size
of the map instead: one unit per feature, per convex-hull point, per
identification, per identification hit and per metadata entry, summed over
subordinates to the same depth writing uses (`map_units`, `feature_units`). Hull
point counts come from `ConvexHull2D::point_count_bound`, so counting is linear
in the number of hulls rather than of points.

### 2. Features are streamed, not retained

Raising the ceilings alone was measured first, and it is not enough: with the
ceilings simply lifted, `FileInfo` read the 59.6 MiB map in 2.87 s at 980 MiB of
peak RSS and the 2.06 GiB map in 103 s at **35.9 GiB**, because the shared parser
builds a `Node` tree for the whole document and a `<pt x=… y=…/>` hull point
costs roughly 1.2 KB of tree for 46 bytes of text.

`identification_xml::Detach` lets a dialect name one root child whose children
are handed over as they close instead of being retained:
`featurexml::read_document` passes `container: "featureList", element: "feature"`
and converts each feature into the map as the parser finishes it. The tree in
memory is therefore one feature, not the whole list. The payload charged for a
detached subtree is refunded when the callback returns, keeping that ceiling a
statement about storage held; the work charged for it is not, because the effort
was really spent.

The header — everything the FeatureXML 1.9 sequence places before `featureList`:
`UserParam`, `dataProcessing`, `IdentificationRun`, `UnassignedPeptideIdentification` —
is converted on the first detached feature, from the root as parsed so far, or
after the parse when the list is empty. One consequence is checked and
deliberate: because a schema-valid `featureMap` ends with `featureList`, a root
child *after* it is now refused with
`featureMap content after featureList` rather than silently read too late to
inform any feature.

### 3. The decode no longer copies the document

`identification_xml::document` copied the byte buffer into a `String` and then
`str::replace`d line endings twice, so a document briefly existed three times
over. It now moves the buffer when the encoding allows (UTF-8, and ASCII-only
Latin-1 or US-ASCII, which is what every OpenMS writer produces), reserves the
exact worst case once for a real Latin-1 decode, and normalises line endings
only when a carriage return is present. `Read::read_to_end` is replaced by
`read_all_bounded`, which grows with `Vec::try_reserve`, so a document larger
than the host can hold is a checked error rather than an allocation abort — the
thing that makes an 8 GiB default ceiling safe to offer.

## Measured, against the C++ Release build

`FileInfo -in <file> -threads 1 -no_progress` on `dax` (AMD EPYC, 2.2 TiB RAM),
page cache warm, three repetitions, medians. C++ is the pinned Release build at
core `bc9cc12` (`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`);
Rust is this branch, `--release`.

| file | features | C++ wall | C++ peak RSS | Rust wall | Rust peak RSS | before this branch |
|---|---|---|---|---|---|---|
| `UPS1_50amol_R1.featureXML`, 59.6 MiB | 42,789 | 0.82 s | 55 MiB | 1.24 s | 150 MiB | refused, exit 3 |
| `UPS1_500amol_R3.featureXML`, 2.06 GiB | 826,019 | 21.5 s | 636 MiB | 44.2 s | 5.39 GiB | refused, exit 3 |

The Rust and C++ `FileInfo` reports are byte-identical on both files and on the
sanity featureXML, except for the C++ `FileInfo took …; Peak Memory Usage: …`
footer, which this port does not print. That includes every feature count, range
bound, total ion current, charge histogram and identification histogram.

The port is 1.5x (59.6 MiB) to 2.1x (2.06 GiB) the C++ wall time and 2.7x to
8.5x its peak memory. The remaining memory gap is the decoded document, which is
held as one `String` for the whole parse (2.06 GiB of the 5.39 GiB), plus the
per-feature `MetaInfo` maps of the retained `FeatureMap`. Removing the first
would mean decoding incrementally from the reader instead of up front; it is a
worthwhile follow-up and is not attempted here.

Two limits of this change are worth stating plainly. The writer still builds the
whole output tree before rendering, so writing a map of this size costs what
reading one used to; the HPC-scale test round-trips the 59.6 MiB map, not the
2.06 GiB one. And identification content — `IdentificationRun` and
`UnassignedPeptideIdentification` — is not streamed, because it is retained in
the `FeatureMap` anyway; a document whose bulk is identification rather than
feature data is therefore held as a tree, bounded linearly by the payload
ceiling rather than tightly.

## What is still refused

The ceilings are linear in the input, never sublinear, so nothing about this
change lets a document amplify a few bytes into unbounded work:

- a document larger than `max_xml_bytes` (8 GiB by default) is refused while
  decoding, with `identification XML byte limit exceeded`;
- a `featureMap` with no `featureList` is refused before any tree is built,
  because nothing in it could be streamed;
- an element density above one per 8 bytes, a work charge above 64 per byte or a
  payload charge above 32 per byte is refused. No streamed featureXML comes
  close; a flood of empty elements, and a document whose bulk is a single
  feature too large to stream, do not fit;
- subordinate nesting beyond `max_depth` (128) and XML nesting beyond
  `2 * max_depth + 16` are refused unchanged;
- every absolute ceiling a caller sets explicitly still wins over the
  size-derived one;
- an allocation the host cannot satisfy is a checked error, not an abort.

`tests/featurexml.rs` pins all of this:
`every_cumulative_ceiling_is_earned_by_the_documents_own_bytes`,
`an_absolute_ceiling_still_wins_over_the_size_derived_one`,
`streamed_features_are_not_retained_and_the_feature_list_must_come_last`,
`the_writer_ceiling_is_earned_by_the_map_it_is_given`, and the `#[ignore]`d
`hpc_scale_benchmark_featurexml_files_load_and_round_trip`, which reads both
benchmark files by path and round-trips the 59.6 MiB map through the writer.
