# Indexed-mzML offsets

With the `mzml` feature, `format::indexed_mzml::IndexedMzMLDecoder` implements
both public operations of the source `IndexedMzMLDecoder`: footer discovery and
index offset parsing. `format::mzml::has_index` exposes the source `MzMLFile`
convenience probe. These APIs operate on plain files and use no new dependency.

```rust,no_run
use openms::format::indexed_mzml::IndexedMzMLDecoder;

let decoder = IndexedMzMLDecoder::default();
if let Some(start) = decoder.find_index_list_offset("run.mzML")? {
    let offsets = decoder.parse_offsets("run.mzML", start)?;
    for (native_id, byte_offset) in offsets.spectra {
        println!("{native_id}: {byte_offset}");
    }
}
# Ok::<(), openms::Error>(())
```

`find_index_list_offset` searches the final 1023 bytes. The explicit
`find_index_list_offset_with_buffer_size` variant changes that window; small
files search their full contents and zero bytes produces `None`. The scanner
preserves the source's broad regular-expression behavior: namespace or other
prefixes ending in `indexListOffset` match, the first matching tag without
digits gives `None`, only consecutive ASCII digits are read, and an embedded
NUL terminates the probe. It is a linear byte scan without a regex dependency.

Finding a number does not prove valid XML, a valid index, a correct checksum,
or an in-file offset. `has_index` is deliberately the same discovery probe.
`parse_offsets` separately bounds the supplied index offset against file size,
reads the remaining suffix, wraps it in `indexedmzML`, and parses its index.
Returned spectrum and chromatogram vectors preserve order and duplicate IDs.
A repeated section replaces its earlier values, as in source; absent sections
are empty in the owned result. Errors return no partial vectors.

Offsets use nonnegative `u64` values up to `i64::MAX`, the portable signed
63-bit file-address range used by the source platforms. Invalid conversion and
I/O produce errors rather than diagnostic output plus a negative sentinel.
The parser accepts standard UTF-8 index structure, scalar numeric text,
CDATA, built-in/numeric XML references, comments and processing instructions.
DTD declarations, custom entities, misplaced index elements, malformed XML
names/attributes, duplicate attributes and trailing non-XML whitespace fail.
XML attribute CRLF normalization and escaped whitespace are preserved.

The native parser corrects the source DOM sibling loop that skips an offset
when it is the first child without preceding whitespace. Compact XML therefore
retains every offset, including singleton sections. Checked small-file reads,
atomic owned results, strict standard index element placement, UTF-8 decoding,
and bounded allocation replace unchecked or platform-dependent source cases.
These differences are part of the reviewed native API, not a claim of exact
behavior for malformed or undefined C++ inputs.

`IndexReadLimits` defaults to a 1 MiB footer search, 16 MiB index suffix,
1,000,000 parsed offsets, 65,536 bytes per ID, depth 64 and 64 attributes per
element. Replaced sections still consume the cumulative offset limit. Suffix
size is checked before allocation, failed allocations become errors, and
duplicate attribute detection uses a tree set rather than the XML reader's
quadratic default. Offset number text is limited to 64 bytes. These are
caller-configurable limits except for numeric text length; arbitrary settings
can increase memory/time consumption.

[Tests](../tests/indexed_mzml.rs) preserve the upstream indexed-mzML fixture
byte-for-byte and check its literal index position 667742, spectrum offsets
24146/345745 and chromatogram offset 665563. Independent cases cover scanner
semantics, compact XML, duplicate/replaced sections, Unicode/entity handling,
malformed structure and resource boundaries. The [source manifest](../tests/data/indexed_mzml_provenance.json)
pins the exact SDK revision and fixture. No C++ runtime or checksum comparison
has been executed. Indexed writing, indexed spectrum decoding and complete
file-backed experiments remain separate Core SDK work.
