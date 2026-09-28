# B8a lane notes (BLOCK-06a: MzMLFile, FeatureXMLFile, StringUtils)

Branch: worktree-wf_ee239125-5be-34. Prereq worktree-wf_ee239125-5be-1 merged (fast-forward to 9fcb962).

## Gaps (from docs/core-sdk-coverage.json scopes)
- MzMLFile.h (partial): IMPeakType (MS:1003441) ion-mobility peak type; selected-ion drift time
  propagation to spectrum; remaining precursor activation metadata/unit routes; progress of
  loadSize, loadBuffer, storeBuffer and both transform overloads.
- FeatureXMLFile.h (partial): "source ZIP input" unsupported (check what C++ actually does).
- StringUtils.h (partial, only FileInfo scope in src/format/file_info/text_format.rs): not ported:
  toStr(float/double,false), appendToStrLowP, long double overloads, toStr(char), toStr(DataValue),
  toStr(ParamValue), numberLength, and the rest of the header (string utilities/conversions).

## Next step
1. Read skill ~/.claude/skills/openms-port-header/SKILL.md
2. Read C++ StringUtils.h at core bc9cc12 member by member; find existing Rust equivalents (grep).
3. MzMLFile.h member-by-member for the named gaps.
4. FeatureXMLFile ZIP input.

## Findings (session 1)
- C++ sources extracted to scratchpad/B8a/ (git -C <abs path to OpenMS4-tests/packages/core> show bc9cc12:...; only plain,
  literal commands pass the worktree guard — no loops/vars around git or the gate script).
- FeatureXML "ZIP input": XMLFile::parse_ (XMLFile.cpp:141-165) sends 'PK' files to CompressedInputSource ->
  ZipInputStream -> ZipIfstream (libzip): exactly one non-directory entry, else FileEmpty / ParseError
  ("ZIP archive contains N file entries; expected exactly 1: a b"). Rust refuses centrally in
  src/format/path_io.rs:36 -> fixing path_io::open closes it for featureXML AND mzML (MzMLFile::load uses parse_ too).
  Needs `zip` crate (2.4.2 MSRV 1.73; checking newer). Cargo.toml change + THIRD_PARTY_CRATE_DECISIONS entry.
- MzML IMPeakType: kernel MSSpectrum (src/kernel.rs:282) has NO im_peak_type field; SpectrumSettings (metadata/acquisition.rs)
  has one. Reader: MS:1003441 -> IM_CENTROIDED (MzMLHandler.cpp:1646); finish: containsIMData && UNKNOWN -> IM_PROFILE (:253);
  writer :5291. Adding a field to MSSpectrum touches ~376 literals (most use ..Default).
- StringUtils: many private duplicates (cli/context.rs:921/971 toInt32/toDouble, mascot_generic.rs, mascot_xml.rs,
  identification_xml.rs nan parsing). text_format.rs holds to_str/to_str_f32/fixed. Plan: new public
  src/data_structures/string_utils.rs porting the whole header; module cycle check needed (data_structures -> format?).

## Session 1, later
- ORACLES (outside repo): ../oracle/B8a-stringutils (DONE, 5068 records, ibminode06 Release, 2 runs identical)
  and ../oracle/B8a-zip-input (launched; `./run.sh pull`). Fixtures: tests/data/zip_input/*.zip (python zipfile).
- Release libzip links libbz2+libz+libcrypto: methods stored/deflate/bzip2 only (no lzma/zstd).
- Underflow-to-zero is ConversionError on Linux (libstdc++ fast_float); subnormals parse.
- Module graph: concept & param are sinks; data_structures -> param only; metadata -> data_structures.
  => numeric core goes to NEW src/concept/numeric_formatting.rs (NumericFormatting.h is CONCEPT/Detail),
     StringUtils to NEW src/data_structures/string_utils.rs, toStr(DataValue) as MetaValue::to_text in NEW
     src/metadata/value_text.rs; text_format.rs re-exports moved items (paths stay valid).
- zip crate: 7.2.0 is newest with MSRV<=1.85 (1.83); 8.x needs 1.88. Plan features deflate-flate2 + bzip2,
  inside file-compression feature. Owned entry reader: central directory via zip crate, then seek data_start,
  take(compressed_size), flate2 DeflateDecoder / bzip2 BzDecoder, CRC via flate2::Crc.
- Next: write numeric_formatting.rs, string_utils.rs, move text_format core; then path_io ZIP; then mzML gaps.

## Session 2 (branch worktree-wf_cd066dab-0a7-38; merged -34 = my session 1, and -1)
- Lead saved WIP 39ad0bc: numeric_formatting.rs, string_utils.rs (1910 lines), value_text.rs, zip_input.rs,
  path_io ZIP route, tests/string_utils.rs + oracle TSV, tests/zip_input.rs + oracle.tsv, Cargo zip dep.
  Compile state UNKNOWN at session start.
- Next step:
  1. Gate-check: cargo check --features <zip feature> and tests string_utils / zip_input on spock.
  2. Fix compile/test failures; commit StringUtils milestone; commit ZIP milestone.
  3. MzML IMPeakType + drift time + progress; 4. ledger + report.
