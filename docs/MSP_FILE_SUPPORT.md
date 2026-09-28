# MSPFile support

Port of `FORMAT/MSPFile.h` / `MSPFile.cpp` (core `bc9cc12`), NIST peptide
spectral libraries in MSP text. Rust: `src/format/msp.rs`, with
`StringUtils::toFloat` added to `src/format/string_utils.rs`. Tests:
`tests/msp_file.rs` (12) and four unit tests in `msp.rs`. Provenance:
`tests/data/msp_file_provenance.json`. Metabolite MSP libraries are
`MSPGenericFile`, `docs/MSP_GENERIC_FILE_SUPPORT.md`.

## API mapping

| Source member | Rust |
|---|---|
| `MSPFile()` | `MSPFile::new()` (`Default` is the all-false value; `new` carries the source defaults) |
| `MSPFile(const MSPFile&)`, `operator=` | `Copy`/`Clone` |
| `~MSPFile()` | drop |
| `DefaultParamHandler` defaults (`parse_headers`, `parse_peakinfo`, `parse_firstpeakinfo_only`, `instrument`) | `MSPFile::defaults`, `to_param`, `from_param`, `from_param_with_warnings`; typed fields `parse_headers`, `parse_peakinfo`, `parse_firstpeakinfo_only: bool`, `instrument: InstrumentFilter` |
| `load(filename, PeptideIdentificationList& ids, PeakMap& exp)` | `load(filename) -> (Vec<PeptideIdentification>, MSExperiment)`, `load_with_limits`, `read(reader, limits)` |
| `load(filename, AnnotatedMSRun&)` | the same `load`; the pair `(identifications, experiment)` stands in for `AnnotatedMSRun`, which is not ported (`METADATA/AnnotatedMSRun.h`) |
| `store(filename, const AnnotatedMSRun&)` | `store(filename, &MSExperiment, &[PeptideIdentification])`, `write(experiment, identifications, writer)` |
| protected `parseHeader_` | private `parse_header` (applied with `parse_headers`) |
| `Constants::UserParam::IonNames` (read by `store`) | `ION_NAMES` |

## Preserved source conventions

- Lines end at LF only (`std::getline`); a CR stays in the line.
- `Name:` peptides lose every innermost `(...)` group (the source's
  ECMAScript `\((.*?)\)`, where `.` matches neither LF nor CR); the
  modifications come from `Mods=`. Both `Mods=` spellings are read
  (`n/pos,aa,name/...` and `n(pos,aa,name)(...)`), with the source's renaming
  of `Pyro-glu`, `CAM`, `AB_old_ICATd8`, `AB_old_ICATd0`, the N-/C-terminal
  attempt at positions 0 and `len-1` with residue fallback, and silent skipping
  of other positions (the source prints them to `cerr`).
- `Inst=` filtering drops the entry's identification and skips its peaks, but
  still advances the `index=<n>` counter; the pending spectrum's precursor and
  header meta values carry over to the next entry, as in the source.
- Peak lines are split on `\s{2,}|\t|\r` exactly as the source's
  `sregex_token_iterator`; m/z and intensity are converted with
  `StringUtils::toFloat` (one rounding to `f32`); annotations follow the
  source's unquote, `?`, first-space, comma, `/` and `^charge` rules.
- The line that ends a peak block is consumed unread.
- `MW:` sets `(MW + z * PROTON_MASS_U) / z` as the precursor m/z.
- `store` writes the source layout: only `M` with a difference formula of
  16 +- 0.01 Da is spelled `M(O)`; the N-terminal modification is position 0;
  the C-terminal modification is not written; `Inst=it`; intensities scaled to
  a maximum of 10000; `IonNames` entries or `"?"`; numbers in `std::ostream`'s
  six-digit `%g`.

## Native differences

- The source indexes without checks in several places; each is an
  `Error::Parse` here: a `Name:` line without a space, `/` or charge; `MW:`,
  `Comment:` or a peak block before the first `Name:`; a `Mods=` count above
  the listed entries or an entry with fewer than three fields; an annotation
  part without ion name (`"/x"`).
- An unquoted annotation field is `Error::Parse` (the source's
  `ConversionError` from `unquote`).
- `MW:` after a charge below 1 is `Error::InvalidValue`; the source divides by
  zero or by the unsigned wrap-around of a negative charge.
- `store` refuses experiment/identification slices of different lengths and an
  `IonNames` array shorter than its spectrum (source undefined behaviour).
- The source's terminal-modification fallback catches only
  `ElementNotFound`; the Rust sequence API reports every failure as an error,
  so any terminal failure falls back to the residue modification, whose own
  failure is returned.
- Resource ceilings: `ReadLimits` (line length, input size, 200 million
  peaks, 10 million entries), checked before each allocation.
- `ProgressLogger` is not a base class here; no progress is reported.

## Checked boundaries and evidence

- Tier 3: every `START_SECTION` of `MSPFile_test.cpp` (copy, assignment,
  load with and without the `qtof`/`it` filter, the store round trip) with its
  literals, over the byte-identical fixture `MSPFile_test.msp`.
- Tier 4: `f32` peak values and the `MW` precursor, `parse_headers`,
  `parse_firstpeakinfo_only = false` (the class test's commented-out
  expectation `y2-H2O`), `parse_peakinfo = false`, each refusal above, the
  consumed line after a peak block, both limits, and the writer layout of two
  entries (`M(O)`, `Mods=2/0,A,Acetyl/11,M,Oxidation`, scaling to 10000,
  `IonNames`).
- No C++ was executed for the writer's `MW:` digits; they follow
  `ostream_g(…, 6)`, the crate's shared `std::ostream` formatter.
