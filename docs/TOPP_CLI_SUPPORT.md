# The TOPP command-line framework

This group ports the `OpenMS4-cli` package — `TOPPBase` and the parameter
records around it — and the first TOPP tool built on it. It is the layer that
turns the library into executables, and it was previously outside the coverage
ledger entirely.

The framework lives in `src/cli.rs` and needs the `paramxml` feature, because
every TOPP tool supports `-ini` and `-write_ini`.

## API

| Source | Native |
| --- | --- |
| `TOPPBase` subclass with `registerOptionsAndFlags_` / `main_` | `trait Tool` with `register(&mut ToolSpec)` / `run(&ToolContext)` |
| `TOPPBase::ExitCodes` | `ExitCode`, same 15 variants and discriminants |
| `ParameterInformation` | `ParameterInformation`, same fields |
| `ParameterInformation::ParameterTypes` | `ParameterType`, same 16 variants |
| `register{String,Int,Double}Option_`, `registerFlag_` | `ToolSpec::register_*` |
| `registerInputFile_`, `registerOutputFile_`, `registerOutputPrefix_`, `registerOutputDir_` | `ToolSpec::register_*` |
| `register{String,Int,Double}List_`, `register{Input,Output}FileList_` | `ToolSpec::register_*_list` |
| `setValidStrings_`, `setValidFormats_`, `setMin/MaxInt_`, `setMin/MaxFloat_` | `ToolSpec::set_*` |
| `registerSubsection_`, `addText_`, `addEmptyLine_` | `ToolSpec::register_subsection`, `add_text`, `add_empty_line` |
| `get{String,Int,Double}Option_`, `get*List_`, `getFlag_`, `getParam_` | `ToolContext::{string,int,double,*_list,flag,param}` |
| `parseRange_` | `parse_range` |
| `main(argc, argv)` | `cli::run::<T>()`, or `run_with::<T>(args, out, err)` |
| `printUsage_` | internal, reached by `--help` / `--helphelp` |

Registration is a builder rather than protected methods on the tool, so the
registered set is inspectable without running the tool, and a tool cannot mutate
its own resolved parameters mid-run.

## Preserved source conventions

The eleven common parameters are registered in source order and with source
descriptions: `ini`, `log`, `instance`, `debug`, `threads`, `write_ini`,
`no_progress`, `force`, `test`, `-help`, `-helphelp`. As in source, the two help
flags are registered with a leading dash, so their command-line tokens carry two.

Resolution order is defaults, then the INI file section, then the command line.
`--help` and `--helphelp` short-circuit before validation, so usage prints
without a required parameter being supplied. `-write_ini` writes the resolved
tree and stops, and does not record the `write_ini` or `ini` request itself.
Advanced parameters are hidden from `--help` and listed by `--helphelp`, which
reports how many were hidden.

Validation covers required values, `setValidStrings_` sets, integer and float
ranges, registered file formats and input existence, and maps each failure to
the source exit code: `MISSING_PARAMETERS`, `ILLEGAL_PARAMETERS`,
`INPUT_FILE_NOT_FOUND`.

## Native differences

Errors are typed `Result` values mapped to exit codes at the boundary rather
than exceptions; `Error::Parse` becomes `PARSE_ERROR`, an I/O not-found becomes
`INPUT_FILE_NOT_FOUND`, a permission failure `CANNOT_WRITE_OUTPUT_FILE`, and
`Unsupported` `INCOMPATIBLE_INPUT_DATA`.

Not yet ported, and deliberately deferred: `ToolHandler` and the `.tools.tsv`
manifest discovery, `INIUpdater`, the CTD/CWL/JSON writers behind `-write_ctd`,
`-write_cwl` and `-write_json`, `SearchEngineBase`, `MapAlignerBase`,
`TOPPExternalToolBase`, and `UpdateCheck` — which should stay unimplemented. The
`log`, `instance` and `threads` values are accepted and exposed on the context
but do not yet redirect logging, select an INI instance section, or set a thread
count, because the port is serial.

## DTAExtractor and executed differential evidence

`src/bin/DTAExtractor.rs` is the first TOPP tool. Its three upstream tests are
reproduced in `tests/topp_dta_extractor.rs` against the retained C++ outputs:

| Upstream test | Arguments | Result |
| --- | --- | --- |
| `TOPP_DTAExtractor_1` | `-rt :61` | `DTAExtractor_RT60.0.dta` byte-identical |
| `TOPP_DTAExtractor_2` | `-level 1` | `DTAExtractor_RT60.0.dta` byte-identical |
| `TOPP_DTAExtractor_3` | `-level 2 -mz :1000` | `DTAExtractor_RT140.0_MZ5.0.dta` byte-identical |

These fixtures were produced by the C++ tool, so agreement is **tier 1 executed
differential evidence** under `docs/DIFFERENTIAL_VALIDATION.md` for the whole
chain: command line, parameter validation, mzML reading, the source number
formatter that names the output files, and DTA writing. This is the first
validated TOPP workflow in the port.

Two source behaviors had to be matched to get there, and both were real gaps:

1. **Header list counts are advisory on reading.** The port rejected an mzML
   whose `count` attribute disagreed with the number of children. The upstream
   fixture `DTAExtractor_1_input.mzML` declares `softwareList count="5"` with
   four entries and `dataProcessingList count="3"` with one, and C++ loads it —
   so the port could not read its own reference data. Reading now ignores the
   declared value; writing still emits the true count.
2. **`DTAFile::store` discards what DTA cannot represent.** The native writer
   refused to drop spectrum metadata, extra precursors, identifications or
   precursor acquisition metadata. Source writes the precursor mass and charge
   plus the peaks and silently ignores the rest, and uses the legacy proton mass
   (`(mz - 1.0) * charge + 1.0`), not the exact one. `dta::WriteOptions::source()`
   selects that behavior; the checked default is unchanged for library callers.

## MzMLSplitter

`src/bin/MzMLSplitter.rs` is the second tool, and the first whose output is
compared **canonically rather than byte for byte**. Both upstream
`TOPP_MzMLSplitter_*` invocations agree with the retained C++ parts on record
counts, native identifiers, MS levels, retention times and peak values. The two
mzML writers differ in serialisation detail, and the upstream test itself uses
FuzzyDiff rather than a byte comparison, so byte equality is not the contract.
Source conventions preserved: the part count derived from a file size in
KB/MB/GB base 1024, the remainder spread over the parts still to come, zero
padding to the width of the part count, and the refusal of `no_chrom` together
with `no_spec`.

## MapNormalizer and SpectraFilterWindowMower

`MapNormalizer` scales MS1 peak intensities to a percentage of the run maximum;
its upstream test reproduces the retained C++ output, and the most intense MS1
peak lands on 100. Higher MS levels are untouched and the source's commented-out
chromatogram branch is not ported.

`SpectraFilterWindowMower` is the first tool with an **algorithm subsection**,
the shape most remaining TOPP tools take. `Tool::subsection_defaults` ports
`getSubsectionDefaults_`: the algorithm's parameter tree is merged beneath the
tool's own defaults, so `-write_ini` emits it and an INI or command line can
override it. Its upstream test matches the retained output, and three further
tests cover the subsection itself — that the C++ `WindowMower` defaults
(`windowsize` 50, `peakcount` 2, `movetype` slide) reach the INI, that changing
`peakcount` changes the result, and that a value violating a registered
restriction is **ignored in favour of the default** rather than rejected, which
is what the source `Param::update` does.

One ordering detail the subsection support forced: a section description cannot
be set on a section that holds no entries, so `ToolSpec::to_param` no longer
writes subsection descriptions and the caller applies them after inserting the
algorithm defaults.

`BaselineFilter` removes the baseline by morphological filtering and reproduces
its retained output with zero difference across all 132 intensities. It takes
the filter's three parameters as ordinary options rather than a subsection,
exactly as the source does, and keeps the source's two refusals: a run holding
only chromatograms, and spectra that are not sorted by m/z.

The second DTA finding above is the first concrete instance of the port's "checked boundaries"
convention blocking C++ parity. The resolution pattern — keep the guard as the
library default, add an explicit source-behavior option, and have the tool opt
in — is the one to apply as further tools meet their own guards.
