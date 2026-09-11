# mzML Product transport

The `mzml` feature reads and writes `MSChromatogram::product` and the ordered `MSSpectrum::products` vector through the existing `mzml::read` and `mzml::write` APIs. This completes transport of the target m/z produced by native extracted-ion chromatograms. It also preserves asymmetric isolation offsets and representable Product metadata.

```rust
use openms::kernel::{MSChromatogram, MSExperiment};
use openms::metadata::Product;
use openms::format::mzml;

fn main() -> openms::Result<()> {
let chromatogram = MSChromatogram {
    native_id: "transition_1".into(),
    product: Product {
        mz: 250.125,
        isolation_window_lower_offset: 0.4,
        isolation_window_upper_offset: 0.75,
        ..Default::default()
    },
    ..Default::default()
};
let experiment = MSExperiment {
    chromatograms: vec![chromatogram],
    ..Default::default()
};
let mut bytes = Vec::new();
mzml::write(&mut bytes, &experiment)?;
Ok(())
}
```

The target uses `MS:1000827`; lower and upper offsets use `MS:1000828` and `MS:1000829`. Values are binary64. Explicit units must denote m/z (`MS:1000040`); missing units use this known quantity. A finite signed target is valid, as in the source Product model. Offsets must be finite and nonnegative. The source reader skips negative offsets; the native reader reports an error. No numerical rounding or interval symmetrization is introduced.

The writer always emits `<product><isolationWindow>`, including the target when it is zero. Positive offsets and negative-zero offsets are emitted; positive zero is represented by its default. No distinction is retained between missing and explicitly zero quantities. Product and precursor quantities remain independent. Peak arrays, array annotations, record names and record metadata retain their existing behavior under both uncompressed and zlib modes.

Product userParams are stored in `product.cv_terms.metadata`. Parsing follows `XMLHandler::fromXSDString`: double, float and decimal map to native f64; the source's six small-integer types use checked i32 conversion then widen to i64; long and arbitrary-integer types use checked i64. Their nominal signedness and narrower XSD ranges are not separately enforced by the source or this adapter. Unknown types, including boolean and absent type, remain strings. Floating metadata must be finite. The existing native numeric parser supplies checked decimal conversion; unsigned-long values above i64 cannot be represented.

Scalar metadata writes as `xsd:string`, `xsd:integer` or `xsd:double`. Empty strings are preserved. Optional MS/UO units retain accession, name and vocabulary reference; missing names remain empty, and a missing vocabulary reference is inferred from the accession prefix. Mismatched or unsupported vocabulary identities and unit attributes without an accession are errors. No ontology lookup, name rewriting or unit conversion is performed. XML escaping preserves Unicode and attribute tabs/newlines.

Both source and native Product models provide a generic CV list with ordinary metadata. The source mzML handler only transports the three isolation quantities and ordinary metadata. Arbitrary entries in that CV list, including present empty accession buckets, cannot be recovered through this source handler and are rejected before writing. Unknown Product isolation CVs are read errors. Native Empty and list metadata alternatives are also rejected before writing because the source serializes them as strings and loses their type. The writer preserves supported metadata as userParams rather than promoting names to vocabulary terms by ontology lookup.

Only one Product is represented per chromatogram. Duplicate containers, multiple isolation windows, repeated supported quantities, duplicate metadata names, and misplaced Product children are errors. Referenceable groups apply through exactly the same handlers, duplicate detection, parameter-count and byte accounting as inline values. Read errors return no partial experiment. Full writer preflight runs before the first output call; transport failures may still leave partial output. Spectrum Product vectors are supported by the [record settings extension](MZML_SETTINGS_SUPPORT.md), which reuses this same codec and preserves list order.

[Nine focused tests](../tests/mzml_product.rs) cover source literals 18.88/1/2 and `isolationwindow3`, ownership, both compression modes, typed metadata and units, exact i64 endpoints, grouped/direct conflicts, cumulative parameter bounds, malformed structure, source signed-target behavior, eleven late writer-rejection cases, and an actual two-region XIC extraction followed by mzML roundtrip. Writer output is independently checked against the bundled mzML XSD when `xmllint` is available. The source literals come from a spectrum Product in `MzMLFile_1.mzML`; the test explicitly routes the same values through the source's parallel chromatogram branch. It does not claim an upstream chromatogram fixture or C++ execution.

The original [Product provenance](../tests/data/mzml_product_provenance.json) pins its inspected source paths to `54a232fe2cae9c590d5c997fa49d20e7769860fb`. The [settings provenance](../tests/data/mzml_settings_provenance.json) records the later spectrum-list extension and signed-zero preservation at 82ce. The handler, Product definition, schema and data fixture are unchanged from the preceding SDK pin; the class test changed only an unrelated string loop variable. See [general mzML support](MZML_SUPPORT.md) and [parameter groups](MZML_PARAM_GROUPS_SUPPORT.md) for the remaining adapter boundaries and configurable limits.
