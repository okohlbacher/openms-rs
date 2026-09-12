# FORMAT integration wave

This wave integrates the fifteen format modules left on Claude's reviewed branches,
then applies independent source, documentation and regression review. The source
pin is OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.
It adds library adapters; it does not add a TOPP executable or complete the SDK.

| Native module | Feature needed | Support and remaining boundary |
|---|---|---|
| `mztab` | none | [Cell and document records](MZTAB_SUPPORT.md); identification/consensus exporters and streams remain open |
| `mztab_file` | none | [Proteomics read/write](MZTAB_FILE_SUPPORT.md); exporter-dependent store overloads remain open |
| `mztab_m` | none | [Metabolomics records/export/write](MZTAB_M_SUPPORT.md); original OMS-input exporter comparison remains open |
| `ms_data_writing_consumer` | `mzml` | [Streaming consumer](MS_DATA_WRITING_CONSUMER_SUPPORT.md); indexed footer, semantic validation and progress remain open |
| `sv_out_stream` | none | [Separated-value output](SV_OUT_STREAM_SUPPORT.md) |
| `mzxml` | `mzml` | [mzXML](MZXML_SUPPORT.md); XSD validation, progress and FileHandler dispatch remain open |
| `mzdata` | `mzml` | [mzData](MZDATA_SUPPORT.md); XSD/semantic validation, progress and FileHandler dispatch remain open |
| `pepxml` | `idxml` | [pepXML](PEPXML_SUPPORT.md); caller-built spectrum lookup replaces automatic source-file loading |
| `qcml` | `paramxml` | [qcML records and files](QCML_SUPPORT.md); `collectQCData` metrics remain open |
| `msstats` | none; tests use `consensusxml` | [MSstats exports](MSSTATS_SUPPORT.md); converter executable is separate work |
| `percolator_infile` | none | [PIN load and feature preparation](PERCOLATOR_INFILE_SUPPORT.md) |
| `transformation_xml` | `featurexml` or `consensusxml` | [Transformation XML](TRANSFORMATION_XML_SUPPORT.md); XSD validation and resolving b-spline models remain open |
| `mascot_generic` | none | [Mascot generic files](MASCOT_GENERIC_SUPPORT.md); progress integration remains open |
| `mascot_xml` | `idxml` | [Mascot XML](MASCOT_XML_SUPPORT.md); arbitrary title-lookup regex remains open |
| `mzidentml` | `idxml` | [Linear and cross-link mzIdentML](MZIDENTML_SUPPORT.md); XSD/semantic validation and generic handler interface remain open |

The [public-header ledger](CORE_SDK_COMPLETION.md) records 21 header assessments
for this wave: five native equivalents and sixteen partial implementations.
The private `MzIdentMLDOMHandler` is mapped in its module's support document,
not counted as an additional public header. Structural XML checks are useful
regressions but do not implement the source XSD validators.

Each module has a `tests/data/*_provenance.json` manifest. Source literals,
independent Rust tests and comparisons with retained C++ output have different
strengths; see [the comparison policy](DIFFERENTIAL_VALIDATION.md). This
integration does not claim a newly executed full C++ SDK differential.
Uncopied large source fixtures are recorded separately from retained fixtures.
The shared [C++ issue log](../OpenMS_CPP_ISSUES.md) distinguishes source-reviewed
findings, unconfirmed candidates and a retained queue of unreviewed leaf reports.
MSstats uses a separately pinned tests package and its test-data gitlink.

Compilation and tests run on IBMI `kim` using node-local NVMe `/scratch`, with
separate checkout and build directories for concurrent workers. Final executed
checks and review evidence belong in [validation](VALIDATION.md).

The [next planned family](PORTING_WAVES.md) is SQLite-backed storage: connector, spectrum access,
sqMass, OSW and OMS. These remain core work. After that, continue the remaining
FORMAT and analysis dependencies from the ledger before extending CLI and TOPP
workflows; a format name alone is not evidence of complete method coverage.
