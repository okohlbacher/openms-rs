# mzML filesystem operations

`format::mzml::{load, load_with_options, load_into, load_into_with_options,
store, store_with_options}` connect the existing mzML subset to the shared
filesystem transport. The target is Core SDK `54a232fe2cae9c590d5c997fa49d20e7769860fb`.
No new dependency is needed: the `mzml` feature already enables file compression.

`load` applies source-default `LoadOptions`, including stable m/z and RT sorting
with aligned annotations. `load_with_options(path, scientific, limits)` supplies
both scientific filters and the independent `ReadOptions` resource limits.
The existing `read` and `FileHandler::load_experiment` retain their established
stream-order behavior. Callers requiring source-default scientific loading can
use the new mzML path API explicitly. This does not add scientific options to
the other FileHandler formats.

Input compression is detected from bytes, even with an unknown or misleading
suffix. Plain, gzip and bzip2 input are supported; limits count decompressed XML
and decoded binary arrays. Invalid options fail before opening the input, and
corrupt compressed streams fail rather than returning a partial experiment.
ZIP input remains unsupported. `load_into` and `load_into_with_options` parse an
owned replacement before changing the destination, so every reported error leaves
the previous value intact. Source `MzMLFile::load` resets its destination first;
atomic failure is a deliberate native correction.

`store` and `store_with_options` publish through a sibling temporary file after
serialization, compression, flush and synchronization succeed. Existing files
survive any reported failure. Unknown and mismatched format suffixes are accepted,
as in direct source `MzMLFile::store` / `XMLFile::save_`; callers needing format
dispatch should use FileHandler. The shared transport recognizes `.gz` and `.bz2`
case-insensitively, whereas the source suffix checks are case-sensitive. ZIP
output is rejected explicitly; it is not written as plain XML under a ZIP name.
Binary-array zlib compression in `WriteOptions` is independent of outer filename
compression. The writer's existing representation checks remain in force.

This completes the file/owned-replacement entry points for the represented mzML
subset. It does not add source DocumentIdentifier loaded-path/file-type state,
full experimental metadata, XML schema/CV validation, `loadSize`, consumers,
transform passes, centroid inference, indexed output, or Numpress codecs.
`Read`/`Write` stream APIs already cover the native equivalent of buffer I/O.
See [scientific loading](MZML_LOAD_OPTIONS_SUPPORT.md) and
[mzML representation support](MZML_SUPPORT.md) for the remaining limits.

[Tests](../tests/mzml_paths.rs) cover source-derived projection filtering, sort
defaults and legacy stream order, aligned annotations, magic-based compression,
decoded limits, corrupt/missing files, independent outer/binary compression and
atomic replacement/output. [Provenance](../tests/data/mzml_paths_provenance.json)
records source hashes. No C++ runtime or full TOPP workflow has been executed.
