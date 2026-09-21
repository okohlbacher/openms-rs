# TransformationXMLFile support

`FORMAT/TransformationXMLFile.h` and `FORMAT/TransformationXMLFile.cpp` at
`bc9cc12`, ported to `src/format/transformation_xml.rs` and tested by
`tests/transformation_xml.rs`. Fixtures, hashes and source anchors are in
`tests/data/transformation_xml_provenance.json`.

A TrafoXML document holds exactly one `<Transformation>`: the fitted model's
name, the `<Param>` entries that model was fitted with, and the `<Pair>`
elements it was fitted from. The Rust module is gated on
`any(feature = "featurexml", feature = "consensusxml")`, which is what brings in
the shared bounded XML codec (`src/format/identification_xml.rs`) this module
reuses for parsing, escaping and rendering. No new Cargo feature was added.

## API mapping

Every public and protected member of the header, plus the two inherited members
the class test exercises.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `TransformationXMLFile()` | `VERSION`, `SCHEMA_LOCATION`, `ReadOptions::default()`, `WriteOptions::default()` | The constructor's whole effect is to pin the handler version to `1.1` and the schema to `/SCHEMAS/TrafoXML_1_1.xsd`. There is no handler object here; the two constants and the option defaults carry that state. |
| `void load(const std::string&, TransformationDescription&, bool fit_model = true)` | `load`, `load_with_options`, `read`, `read_with_options` | The out-parameter becomes the return value. `fit_model` is `ReadOptions::fit_model`. `load`/`read` take the default; the `_with_options` forms take the limits too. |
| `void store(const std::string&, const TransformationDescription&)` | `store`, `store_with_options`, `write`, `write_with_options` | `store` publishes the file only after the bytes are complete and flushed, through the crate's atomic-publication helper. |
| `void onStartElement(const char16_t*, const Internal::XMLAttributes&)` | not ported: a SAX callback of the source's own handler base. Its whole body is the element dispatch, which is reproduced inside `read_record_with_options`. |
| `Param params_` | `TransformationRecord::parameters` | Protected state in the source, public here: the untyped `(name, value)` map is exactly what the file contains, and is the only way to see the parameters of a model this port cannot fit. Typed as `BTreeMap<String, ParamValue>`, so document order is replaced by name order — the source's `Param` is also name-ordered. |
| `TransformationDescription::DataPoints data_` | `TransformationRecord::data` | `Vec<DataPoint>`, keeping the optional `note`. |
| `std::string model_type_` | `TransformationRecord::model_type` | |
| inherited `bool Internal::XMLFile::isValid(const std::string&, std::ostream&)` | `is_valid`, with the optional `xml-schema` feature: XSD validation against the bundled, unchanged `TrafoXML_1_1.xsd`; see [XML schema validation](XML_SCHEMA_SUPPORT.md). |
| inherited `const String& Internal::XMLFile::getVersion()` | `VERSION` | |
| (none) | `TransformationRecord::model_config` | The typed `ModelConfig` a record names. The source resolves the name inside `TransformationDescription::fitModel`; splitting it out is what lets a caller read a `b_spline` file's parameters without being able to fit it. |
| (none) | `TransformationRecord::description` | `fitModel` plus `setDataPoints`, as one call. |
| (none) | `TransformationRecord::from_description` | The inverse: what `store` serialises. |
| (none) | `read_record_with_options`, `load_record_with_options`, `write_record_with_options`, `store_record_with_options` | The record-level surface the source keeps private. |
| (none) | `ReadOptions`, `WriteOptions` | Explicit resource ceilings; the source has none. |

## Preserved source conventions

- **The `fit_model = false` path loses the model name.** `load` always calls
  `setDataPoints`, which resets the model type to `none` *even when it was
  `identity`*, and only then optionally calls `fitModel`. So a description read
  with `fit_model: false` reports `none`, exactly as the class test asserts.
  `read_record_with_options` is the way to see the name in that case.
- **`identity` is sticky.** `TransformationDescription::fit_model` returns early
  when the current model is `Identity`, so loading into a description that
  already holds one keeps it. That behaviour lives in
  `src/analysis/transformations.rs` and is unchanged here.
- **An empty model name is refused.** `store` throws
  `Exception::IllegalArgument` with "will not write a transformation with empty
  name"; `write_record_with_options` returns `Error::InvalidValue` carrying the
  same phrase. A record built from a description can never hit it, because
  `ModelConfig::None` is named `none`.
- **Empty parameters are skipped.** A `ParamValue::Empty` writes no `<Param>`.
- **The three written types.** `int`, `float` and `string`. The source funnels
  all three list types into `type="string"` holding the bracketed list text,
  which cannot be read back as a list; this port reproduces that on writing and
  therefore inherits the same one-way conversion.
- **`<Pairs>` only when there are pairs**, with `count` equal to their number.
- **`<Pairs count>` is advisory.** The source only passes it to
  `data_.reserve` and never compares it with the number of children, so a
  disagreement is accepted here too.
- **`note` is the only escaped attribute**, and is written only when non-empty.
- **Pre-3.0 weight tolerance.** An empty `x_weight`/`y_weight` string means
  unweighted, which is how TrafoXML files written before OpenMS 3.0 spell it.
- **Weight names.** `x`/`y` identity, `ln(x)`/`ln(y)`, `1/x`/`1/y`,
  `1/x2`/`1/y2`, mapping onto the crate's `WeightFunction`; anything else is
  refused, as `TransformationModel`'s constructor refuses it.
- **An unknown model name is an error**, as `fitModel` throws
  `Exception::IllegalArgument` for one.

## Native differences

Each of these is stated at the item in the rustdoc as well.

- **A document version above `1.1` is refused.** The source compares the
  version numerically and only warns, with the message "This might lead to
  undefined program behavior". A parser may not choose undefined behaviour, and
  `src/format/idxml.rs` already sets the precedent of a bounded accepted range,
  so `read_record_with_options` returns `Error::Unsupported`. Versions `1.0`
  and `1.1` are accepted; both upstream fixtures are `1.0`.
- **An unsupported `<Param type>` is an error.** The source reaches
  `XMLHandler::error`, which only writes to the error log and returns, so the
  parameter is silently dropped. This port returns `Error::Unsupported` naming
  the type.
- **An unknown element is an error.** The source reaches
  `XMLHandler::warning`, which logs at debug level in a release build and drops
  the element. This port returns `Error::Unsupported`, matching what
  `consensusxml`/`idXML` already do in this crate.
- **A non-finite coordinate or parameter is an error.** The source parses
  `from`, `to` and float parameters with `attributeAsDouble_` and accepts
  whatever `strtod` returns, including infinities. This port requires finite
  values, because the fitted models require them.
- **`isValid` needs the `xml-schema` feature.** It is XSD validation from the
  source's `Internal::XMLFile` base, ported as `is_valid` over libxml2. File 3
  is not even well-formed XML, carrying a second `</Transformation>` close tag,
  so it is `Error::Parse` where the source returns `false`. Without the
  feature, the structural reader reaches the same verdict on all four upstream
  fixtures.
- **A different parameter set is written for a data-fitted linear model.** In
  C++ the base `TransformationModel` constructor copies the caller's `Param`
  and the data-fitted branch merges the eight weight and datum defaults into
  it, so that path writes ten `<Param>` entries while the explicit-coefficient
  path writes two. This port writes `slope` and `intercept` always, and the
  weight parameters only when a weight is not the identity, so the file no
  longer depends on how the model was constructed. Both forms read back to the
  same model, and the class test's two-parameter assertion still holds.
- **Weighted models round-trip.** The source writes the weight parameters only
  when its `Param` happens to carry them, which after `invert()` it does and
  after a plain fit-from-coefficients it does not. Writing them whenever they
  are non-identity means a weighted model survives a store/load cycle instead
  of silently becoming unweighted.
- **Output bytes are not reproduced.** Attribute order within an element is the
  shared writer's (alphabetical) rather than the source's `type`, `name`,
  `value`; the source's two spaces after `<Param` are not emitted; indentation
  is two spaces per level rather than a tab. The upstream suite compares
  trafoXML with `FuzzyDiff`, not byte-for-byte.
- **`b_spline` is refused with `Error::Unsupported`.** The source implements
  it; `src/analysis/transformations.rs` does not. The file's parameters and
  pairs remain readable through `read_record_with_options`.
- **Interpolation types.** `polynomial` and `akima`, which the source's
  interpolated model accepts, return `Error::Unsupported`; `linear` and
  `cspline` map onto the crate's two variants.
- **No OpenMP gap.** Neither the header nor its implementation carries a
  `#pragma omp`.

## Checked boundaries and evidence

Ceilings, all in `ReadOptions`/`WriteOptions` and all checked before anything
is allocated or mutated:

| Ceiling | Default | What it bounds |
|---|---|---|
| `max_xml_bytes` | 64 MiB | undecoded input, or serialised output |
| `max_records` | 1,000,000 | XML elements, `<Pair>` included |
| `max_list_items` | 1,000,000 | entries in one decoded list |
| `max_payload_bytes` | 256 MiB | cumulative decoded-tree allocation |
| `max_work` | 50,000,000 | cumulative parser work units |
| `max_data_points` | 1,000,000 | `<Pair>` elements, checked against both the declared `count` and the actual child count before any pair is stored |

A zero ceiling is `Error::InvalidValue`. Reading builds the whole record and
only then constructs the description, so a failure leaves the caller's data
untouched; writing renders into a buffer and publishes the file only after the
bytes are complete.

No string is byte-sliced anywhere in the module: attribute values are compared
and parsed whole, and the only slicing is `str::parse`. The escaping round trip
is tested with `日本語 & <anchor>` and with `"quoted"` notes.

**Evidence: tier 3, source review.** The expectations in
`tests/transformation_xml.rs` are transcribed class-test literals — model names
`none`/`linear`/`interpolated`, parameter counts 0/2/2, slope π and intercept e
compared with `TEST_REAL_SIMILAR`'s relative tolerance (the fixture stores them
at six significant digits), the three pairs (1.2, 5.2) (2.2, 6.25) (3.2, 7.3),
`extrapolation_type` `two-point-linear` appearing in a reloaded interpolated
model, the `fit_model = false` path reporting `none`, the four fixtures' `isValid` verdicts (XSD
validation with `xml-schema`, structural outcomes without it), and `Exception::IllegalArgument` for the model name
`mumble_pfrwoarpfz`. Transcribed literals detect transcription drift but cannot
falsify a misread algorithm. No C++ was built or executed and no C++ output was
retained, so nothing here is tier 1 or 2. The escaping round trip, the
non-ASCII note, the coordinate-weight round trip, the version and
unknown-element refusals and the resource ceilings are independently derived
(tier 4), because no upstream fixture reaches them.

### Section accounting

All four `START_SECTION`s of `TransformationXMLFile_test.cpp` are ported. The
XSD `isValid` section is ported behind the optional `xml-schema` feature, in
`tests/xml_schema.rs` and `tests/xml_schema_formats.rs`; the structural
substitute below still runs without it.

| Section | Assertion macros | Rust test | One reproduced value |
|---|---|---|---|
| `TransformationXMLFile()` | 1 | `the_constructor_state_is_version_1_1_and_its_schema` | the handler version `1.1` |
| `[EXTRA] static bool isValid(const std::string&)` | 4 | `tests/xml_schema_formats.rs::transformation_xml_fixtures_validate_as_the_class_test_asserts` (`xml-schema`); structurally, `the_three_schema_valid_fixtures_read_and_the_invalid_one_does_not` | fixtures 1, 2 and 4 validate against `TrafoXML_1_1.xsd`; 3 is malformed XML |
| `void load(..., bool fit_model=true)` | 17 | `loading_fits_the_named_model_from_the_file` | `getModelType()` of file 4 is `interpolated`, and its second pair is (2.2, 6.25) |
| `void store(...)` | 18 (its `#if 0` b_spline block is excluded, as the compiler excludes it) | `storing_and_reloading_preserves_the_model_and_its_parameters` | a reloaded stored interpolated model has `params.size() == 2` with `extrapolation_type == "two-point-linear"` |

### C++ findings

Recorded for `OpenMS_CPP_ISSUES.md`; the integrator owns that file.

- **`symmetric_regression` is dead.**
  `TransformationModelLinear.cpp:34` reads it into `symmetric_` and nothing
  ever reads `symmetric_` again, while `TransformationModelLinear.h:24`
  documents that "Depending on parameter `symmetric_regression`, a normal
  regression (*y* on *x*) or ..." is performed. The documented alternative
  regression is not implemented. A TrafoXML carrying
  `symmetric_regression="true"` is accepted and ignored, by C++ and by this
  port alike.
- **Unconfirmed error-policy candidate (CPP-176): unsupported `<Param type>`.**
  `TransformationXMLFile.cpp:161` calls `XMLHandler::error`, which
  (`XMLHandler.cpp:71`) only logs. A file whose model parameters use any other
  type loads as a model with missing parameters, and `fitModel` then either
  falls back to defaults or fails for an unrelated reason. Non-fatal continuation
  is explicit source policy; a public-contract violation has not been established.
- **Forward-compatibility policy: a too-new document version is accepted.**
  `TransformationXMLFile.cpp:137` warns that this "might lead to undefined
  program behavior" and continues. This alone does not establish a source defect.
