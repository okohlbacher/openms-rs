# Compressed input streams, ZIP archives and the XML string helpers

Headers (core `bc9cc12`): `FORMAT/GzipIfstream.h`, `FORMAT/Bzip2Ifstream.h`,
`FORMAT/ZipIfstream.h`, `FORMAT/ZipArchiveFile.h`,
`FORMAT/HANDLERS/StringManager.h`, `FORMAT/HANDLERS/XMLAttributes.h`.

| Header | Rust | Feature |
|---|---|---|
| `GzipIfstream.h` | `src/format/gzip_ifstream.rs` | `file-compression` |
| `Bzip2Ifstream.h` | `src/format/bzip2_ifstream.rs` | `file-compression` |
| `ZipIfstream.h` | `src/format/zip_ifstream.rs`; ZIP detection in `src/format/path_io.rs` | `file-compression` |
| `ZipArchiveFile.h` | `src/format/zip_archive_file.rs` | `file-compression` |
| `HANDLERS/StringManager.h` | `src/format/string_manager.rs` | none |
| `HANDLERS/XMLAttributes.h` | `src/format/xml_attributes.rs` | none |

The `file-compression` feature is in every default build through `mzml`,
`featurexml`, `consensusxml`, `paramxml` and `cv-mapping`; it now also pulls
in `zip =6.0.0` (see `THIRD_PARTY_CRATE_DECISIONS.md`).

## ZIP input for the loaders (SHARED-02)

`path_io::open`, which every XML loader and `FileHandler` read through, used
to refuse a ZIP signature. It now recognises `PK\3\4` (a local file header)
and `PK\5\6` (an archive with no entries) and reads the archive's single
non-directory entry through `ZipIfstream`, as the source's
`XMLFile::parse_` -> `CompressedInputSource` -> `ZipInputStream` ->
`ZipIfstream` chain does. The source dispatches on the two bytes `PK`; four
bytes keep a text file that happens to start with "PK" readable. An archive
with zero or several file entries, an unreadable archive and an entry in an
unsupported method fail with the errors `ZipIfstream::open` documents.

`tests/zip_input_oracle.rs` replays lane B8a's Release oracle
(`../oracle/B8a-zip-input`, 52 records): `FeatureXMLFile::load`/`loadSize`
and `MzMLFile::load`/`loadSize` over stored, deflated, bzip2 and LZMA
archives, one with a directory entry, one with two entries, an empty and a
directory-only archive, one with a bad CRC and one with junk after the
signature. `tests/compressed_streams.rs` loads consensusXML, idXML, paramXML
and CV-mapping files from a ZIP archive and compares them with the plain
loads.

## Semantics kept, and how they were established

Evidence tiers: the upstream class tests (`GzipIfstream_test`,
`Bzip2Ifstream_test`, `ZipIfstream_test`, `XMLHandler_test` - which is the
StringManager test - and the UTF-8 section of `MzMLFile_test`) are ported in
`tests/compressed_streams.rs`, `tests/string_manager.rs` and the unit tests
of `string_manager.rs`. Two executed-oracle runs of the Linux x86_64 Release
build on ibminode06 (`../oracle/F7b-streams`, this lane; `../oracle/B8a-zip-input`,
lane B8a) are replayed in `tests/compressed_streams_oracle.rs`, which compares
byte counts, `isOpen`/`streamEnd`, the bytes and exception classes step by step.
The Release build links zlib 1.3.2, libbzip2 1.0.8 and libzip 1.11.2 (with
bzip2 and OpenSSL, without LZMA or Zstandard).

- **Gzip** (`gzread`): transparent for non-gzip input; concatenated members;
  trailing data after a member ignored, also when it starts with the gzip magic
  but has an unusable header; the stream closes when a request cannot be
  filled, so a request ending exactly at the end leaves it open; a failure is
  reported on the call after the bytes decoded before it.
- **Bzip2** (`BZ2_bzRead`): the reader loop is re-created over the `bzip2`
  crate's low-level decompressor, with libbzip2's 5000-byte input blocks, so a
  read closes as soon as the end-of-stream marker is decoded (also when it
  filled the request exactly); only the first stream is read; any failure
  discards the bytes of that call.
- **ZIP** (`zip_fread`): a read fills the request unless the entry ends; the
  stream closes only on a read that returns 0; a CRC-32 or size mismatch is
  reported after the bytes, on the next call.

## Divergences (each is asserted by a test)

| Case | Release build | Port | Why |
|---|---|---|---|
| Truncated gzip member | data so far, then a clean end | `Error::Parse`; `GzipIfstream::set_accept_truncated(true)` gives the source result | no silent data loss |
| gzip CRC-32 / ISIZE mismatch | silent when the member is decoded in one fetch and the read runs to the end; an error otherwise | always `Error::Parse`, after the bytes | the check is the only integrity guard |
| gzip header CRC mismatch; reserved header flags (first member) | empty stream; raw copy of the file | `Error::Parse` | both are corrupt input |
| corrupt deflate data | zlib decodes up to 1 MiB ahead and fails on the next call | fails where flate2 meets the damage | byte counts before the error differ |
| `ZipArchiveFile::addOrReplaceFromFile`, missing source | `InvalidValue` from `zip_close` (header says `FileNotFound`) | `Error::Io` `NotFound`, archive untouched | follows the header |
| `extractEntryToTempFile` on an unpacked layout with `../` | returns the escaping path | refused | path traversal |
| Escape check in `unzipDirectory` / `extractEntryToTempFile` | string-prefix comparison: `../<base>x/f` passes | component-wise | path traversal |
| JSON index names | written verbatim (`"` or `\` gives invalid JSON) | escaped | ordinary names give identical bytes |
| Archive entry order (`zipDirectory`, directory index) | file-system order | sorted | determinism; the oracle's order was the same |
| `StringManager::parseInt` beyond `i32` | `(int)` cast wraps (`4294967296` -> 0) | `Error::InvalidRange`; `parse_int_source` wraps | no silent wrap |
| `StringManager::convert` of an unpaired surrogate | pairs it with the next unit (`D800 0041` -> U+2441) | `Error::InvalidValue` | not valid UTF-16 |
| `appendASCII` with a non-ASCII unit | appends the low byte | `Error::InvalidValue`, string unchanged | precondition checked |
| Duplicate entry names in an archive | libzip counts both | the `zip` crate keeps the later one | crate model |

Not ported, by design: `shared_xerces_ptr`/`unique_xerces_ptr` (Rust
ownership frees the buffers), `CONST_XMLCH` (`str::encode_utf16`), and the
SIMD code paths of `StringManager` (scalar loops, identical results, no
`unsafe`). The source's `Exception::NotImplemented` for builds without libzip
has no counterpart, because the types exist exactly when the feature is on.

## C++ defect candidates

1. `GzipIfstream::read` calls `gzerror` after `close()` has set the handle to
   `nullptr`, and reports `gzerror_`, which the filename constructor never
   initialises: the Release build prints `unknown error (code: 1593851233)`,
   a different number per call site.
2. `GzipIfstream.h` documents `ConversionError`; the implementation throws
   `ParseError`.
3. The Release build's gzip reader accepts a member with a wrong CRC-32 or
   ISIZE when it is decoded in one fetch (oracle `g_bad_crc`, `g_bad_len`),
   and reads a truncated member as a shorter valid one.
4. `ZipIfstream.h` says a short read closes the stream; the implementation
   closes only when a read returns 0 (oracle `z_short`).
5. `ZipArchiveFile`'s escape check compares strings, so `../<base>x/...`
   escapes into a sibling directory; the unpacked-layout branch of
   `extractEntryToTempFile` has no check at all (oracle `extractDirEscape`).
6. `ZipArchiveFile::writeSidecarIndex` writes names unescaped into JSON.
7. `Xerces XMLString::parseInt` wraps `long` values into `int` silently.
8. `GzipIfstream`/`Bzip2Ifstream::read` cast `size_t n` to `unsigned int`;
   bzip2 then passes it as `int`, so requests of 2 GiB or more fail or wrap.
