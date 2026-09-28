# F7b lane notes (restart-safe)

Branch: worktree-wf_cd066dab-0a7-13, base f0990c0. Slot F7b, gate host dax.

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

## Next step
1. Add zip dep; new module src/format/compressed_ifstream.rs (GzipIfstream, Bzip2Ifstream, ZipIfstream); wire ZIP into path_io::open.
2. Copy class-test data into tests/data (provenance manifest), port the 3 class tests.
3. ZipArchiveFile module. 4. StringManager + XMLAttributes. 5. IndentedStream public type. 6. Oracle on ibminode06 in ../oracle/F7b-streams/.
