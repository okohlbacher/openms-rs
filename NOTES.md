# F7b lane notes (restart-safe)

Branch: worktree-wf_cd066dab-0a7-18 (incarnation 2 merged -13), base f0990c0. Slot F7b, gate host dax.

## Findings (took many calls)
- Headers/cpp read at core bc9cc12 (git -C ../OpenMS4-tests/packages/core show bc9cc12:<path>; the sandbox refuses loops/vars around git: run one literal git per call).
- Rust compressed input lives in src/format/path_io.rs `open()`: gzip (MultiGzDecoder), bzip2 (MultiBzDecoder), and ZIP "PK\x03\x04" is REFUSED (Error::Unsupported) = SHARED-02 gap. Coverage entries naming it: CVMappingFile, ConsensusXMLFile, FeatureXMLFile, FileHandler, ConsensusXMLHandler, FeatureXMLHandler, MzMLHandler, IdXMLFile, ParamXMLFile, XMLFile.
- C++ XMLFile::parse_ dispatches on 2 bytes: 'BZ', 1f 8b, 'PK' -> CompressedInputSource -> Bzip2/Gzip/ZipInputStream -> *Ifstream. FileHandler::getTypeByContent uses 'BZ', gzip, 'PK\3\4' and reads 8191 bytes from the *Ifstream.
- C++ semantics: GzipIfstream = zlib gzread (transparent for non-gzip input! concatenated members; trailing garbage ignored; closes when a request cannot be filled). Bzip2Ifstream = BZ2_bzRead SINGLE stream (closes when BZ_STREAM_END seen, even if returned == n). ZipIfstream = libzip, exactly one non-directory entry (names ending '/' skipped): 0 -> FileEmpty, >1 -> ParseError listing names, zip_open failure -> FileNotFound if missing else FileNotReadable; read(n==0) -> 0 no-op; read returning 0 closes.
- Class tests: GzipIfstream_test, Bzip2Ifstream_test, ZipIfstream_test (data in core src/tests/class_tests/openms/data: GzipIfStream_1.gz, GzipIfStream_1_corrupt.gz, Bzip2IfStream_1.bz2, Bzip2IfStream_1_corrupt.bz2, Deisotoper_input1.mzML(.zip), DirEntry.mzML.zip, Empty.zip, MultiEntry.mzML.zip). No tests for ZipArchiveFile, StringManager, XMLAttributes.
- ZipArchiveFile consumers: OpenSwath OSW parquet reader/writer, PercolatorScoring, TransitionParquetFile, FileHandler.cpp:968 (unzipDirectory).
- Crate: zip 6.0.0 (MSRV 1.83; required deps crc32fast, indexmap, memchr) cached locally; use default-features=false, features deflate-flate2 + bzip2. Plan: add to `file-compression` feature. Cargo.lock must be regenerated on dax (gate rsyncs up only; fetch lock back with ssh -J kohlbach@sshgw.cs.uni-tuebingen.de kohlbach@dax cat /scratch/kohlbach/agents/F7b/Cargo.lock).
- Error enum: InvalidValue, Parse{line,message}, Unsupported, Io, InvalidRange, MissingInformation. TempDir at src/system/file.rs:1440 (TempDir::new(keep), .path()).
- StringManager: pure UTF-16 helpers (convert both ways, strLength, parseInt (Xerces), isASCII, appendASCII, compress64_). XMLAttributes: view over Xerces Attributes: value(u16 qname)/value(narrow)/index (neg if absent)/valueByIndex/handle.

## Incarnation 3 findings
- Git: run git as a plain single command (no cd, no &&) or the sandbox refuses it.
- Gzip read error throws ParseError (header says ConversionError); after close() it calls gzerror(nullptr) -> NULL -> "unknown error (code: <gzerror_>)" and gzerror_ is UNINITIALISED in the filename ctor (defect candidate). Bzip2 BZ2_bzReadOpen failure -> ConversionError.
- ZipIfstream_test [EXTRA]: corrupt (text) zip -> FileNotReadable; FileHandler getTypeByFileName strips .zip (sample.mzML.zip -> MZML); end-to-end loadExperiment(Deisotoper_input1.mzML.zip) == plain (sizes, ms levels).
- ZipArchiveFile: zipDirectory (store-only, recursive, rel generic paths, removes existing), unzipDirectory (dir -> returned unchanged; !readable -> FileNotFound; TempDir/parquet_unpacked; abs path / escape -> InvalidValue), addOrReplaceFromFile (ZIP_CREATE; stored), listEntries (empty on failure), writeSidecarIndex (dir -> <path>.idx.json external; zip -> entry '<basename>.idx.json', JSON `{\n  "name": size,\n...\n}\n`, stored), extractEntryToTempFile (dir layout -> path, missing -> FileNotFound).

## Done (incarnation 3, committed)
- src/format/{gzip_ifstream,bzip2_ifstream,zip_ifstream,zip_archive_file,string_manager,xml_attributes}.rs; ZIP in path_io::open; zip =6.0.0 in file-compression; Cargo.lock fetched from dax.
- tests/compressed_streams.rs (20 pass on dax), tests/string_manager.rs, tests/mzml_counts.rs zip assertion now Io(InvalidData) (C++ FileNotReadable).
- StringManager class test is XMLHandler_test.cpp (compress64/isASCII/appendASCII/strLength) + MzMLFile_test.cpp UTF-8 section.
- Existing Release oracle ../oracle/B8a-zip-input (ibminode06): LZMA entry -> FileNotReadable (Release libzip has no LZMA), bad CRC -> ParseError, junk -> FileNotReadable, dir-only/empty -> FileEmpty.
- Own oracle ../oracle/F7b-streams launched on ibminode06 (run.sh push/pull).

## Incarnation 3, late state
- All tests green on dax; oracle ../oracle/F7b-streams pulled (results/records.tsv), differential tests tests/compressed_streams_oracle.rs + tests/zip_input_oracle.rs pass.
- Ledger fragment written (.planning/ledger-fragments/wave4/F7b.json, 16 entries); report draft in scratchpad F7b.md.
- CI sweep F7b-sweep running (rust.yml regenerated locally only for the sweep, then reverted).
- Remaining: sweep result -> write report, delete NOTES.md, final commit.

## Next step (older)
1. `bash ../oracle/F7b-streams/run.sh pull`; copy results/records.tsv into tests/data/compressed_streams/f7b_oracle_records.tsv (+ oracle data files) and write the differential test.
2. Provenance manifest tests/data/compressed_streams/provenance.json; THIRD_PARTY_CRATE_DECISIONS entry for zip.
3. Full CI sweep; MSRV +1.85.0 check; clippy; doc.
4. Report + ledger fragment; delete NOTES.md.
5. If time: IndentedStream public type (src/cli/usage.rs Indented).
