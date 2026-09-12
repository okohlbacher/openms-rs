# On-disc imzML experiment support

Source revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

| Item | Value |
|---|---|
| Ported header | `KERNEL/OnDiscImzMLExperiment.h` (266 lines) + `source/KERNEL/OnDiscImzMLExperiment.cpp` (362 lines) |
| Rust module | `src/kernel/on_disc_imzml_experiment.rs` |
| Test | `tests/on_disc_imzml_experiment.rs` (52 cases) plus 2 doctests |
| Provenance | `tests/data/on_disc_imzml_provenance.json` |
| Feature gate | `mzml` (existing; no Cargo change) |
| Builds on | `src/format/imzml_handler.rs` — see `docs/IMZML_HANDLER_SUPPORT.md` |
| Evidence | Tier 3 source review with one independent cross-check; tier 4 for the bounds and the native checks. Not a tier 1 differential. |
| Ledger status | `evidence_requires_review` |

## Scope

`OnDiscImzMLExperiment` is the kernel-level façade over the imzML reader that
`docs/IMZML_HANDLER_SUPPORT.md` describes. It owns four things the reader does
not: the dataset lifecycle (`open`/`close`/`isOpen`), the two-dimensional pixel
grid built eagerly at `open()`, coordinate-addressed spectrum access, and the
ion-image sweep. Everything below that — the IMS vocabulary, the `.imzML` index
scan, the `.ibd` range reads and the per-spectrum decode — is the reader's and
is called, not reimplemented.

### Three types that belong to another package

`getGeometry()` returns `MSImagingGeometry` and `extractIonImage()` returns
`IonImage`, both forward-declared in the header and both living in
`OpenMS/IMAGING/`, a directory no Rust package owns. The region-scoped
`extractIonImage` overload additionally needs `MSImagingRegion`, and the sweep
itself is the header-only template `Internal::extractIonImage` in
`IMAGING/IonImageExtraction.h`. The façade cannot be ported without them, so
this module reproduces them as `ImagingGeometry`, `ImagingRegion` and
`IonImage`.

That is **not** a claim on those headers. Their own class tests —
`MSImagingGeometry_test.cpp`, `MSImagingRegion_test.cpp`, `IonImage_test.cpp` —
were read, and the behaviours the façade depends on are checked here, but the
suites were not ported section by section and `IMAGING/MSImagingExperiment.h`,
the in-memory counterpart that shares the geometry, is entirely unported. A
later IMAGING work package owns those four headers and may replace these types;
the mapping tables below record exactly how much of each was reproduced, so it
knows what it is inheriting.

## API mapping — `KERNEL/OnDiscImzMLExperiment.h`

Every public member of the header, in declaration order.

| C++ member | Rust | Notes |
|---|---|---|
| `class OnDiscImzMLExperiment` | `OnDiscImzMLExperiment` | |
| `OnDiscImzMLExperiment()` | `OnDiscImzMLExperiment::new`, `Default` | The source allocates an empty `Impl`; the Rust default is closed, zero-length, with an empty grid and default metadata |
| `~OnDiscImzMLExperiment()` | `Drop` (derived) | The source's `Impl::~Impl` `fclose`s the `.ibd`; `File` closes itself, so there is no destructor to write. The source's out-of-line destructor exists only because `Impl` is incomplete in the header — a PIMPL artefact with no Rust counterpart |
| `OnDiscImzMLExperiment(const OnDiscImzMLExperiment&) = delete` | no `Clone` impl | **Not ported, deliberately.** The type owns an open file handle with one seek position; see *Move-only semantics* |
| `operator=(const OnDiscImzMLExperiment&) = delete` | no `Clone` impl | Same |
| `OnDiscImzMLExperiment(OnDiscImzMLExperiment&&)` | native move | **Not ported as an operation:** Rust moves by default. See *Move-only semantics* for the source's moved-from reset |
| `operator=(OnDiscImzMLExperiment&&)` | native move | Same. The source's `close()`-then-move-then-reset body has no counterpart |
| `void open(const std::string&, const std::string& ibd_path = "")` | `open`, `open_with_ibd`, `open_with_limits` | Rust has no default arguments, so the defaulted `ibd_path` becomes `open` (derive it) versus `open_with_ibd` (explicit). `open_with_limits` adds the reader's ceilings, which the source has no way to pass. Returns `Result<()>`; the source's `@throws FileNotFound` / `@throws ParseError` map to `Error::Io` / `Error::Parse` |
| `void close() noexcept` | `close` | Returns `()`; nothing here can fail either |
| `bool isOpen() const noexcept` | `is_open` | |
| `std::size_t getNrSpectra() const noexcept` | `len` | Plus `is_empty`, which the source expresses as `getNrSpectra() == 0` |
| `std::size_t size() const noexcept` | `len` | The header defines `size()` inline as a synonym for `getNrSpectra()`; one Rust method covers both |
| `const ImzMLSpectrumIndex& getIndex(std::size_t i) const` | `index` | `Result<&ImzMLSpectrumIndex>`; `Exception::IndexOverflow` becomes `Error::InvalidValue`. No `.ibd` read, as documented |
| `MSSpectrum getSpectrum(std::size_t i) const` | `spectrum` | `Result<MSSpectrum>`. Also `decoded_spectrum`, which returns the reader's `DecodedSpectrum` so the auxiliary arrays the source warns about and drops are visible rather than logged |
| `MSSpectrum operator[](std::size_t i) const` | **not ported as `Index`** | `std::ops::Index` must return a reference and cannot fail; `getSpectrum` returns a decoded value and can. `spectrum(i)` is the only form. The source's `operator[]` is a one-line forward to `getSpectrum`, so nothing is lost but the sugar |
| `MSSpectrum getSpectrumAtCoord(uint32_t x, uint32_t y, uint32_t z = 1) const` | `spectrum_at_coord`, `spectrum_at_pixel` | The defaulted `z` becomes the second name. Coordinates stay imzML-native 1-based. `Exception::ElementNotFound` becomes `Error::InvalidValue` |
| `IonImage extractIonImage(double mz, double tolerance_ppm) const` | `extract_ion_image` | |
| `IonImage extractIonImage(double mz, double tolerance_ppm, Size region_id) const` | `extract_ion_image_in_region` | Overload set becomes distinct names, per the crate's convention |
| `const ImzMLMeta& getImzMLMeta() const noexcept` | `meta` | Returns the reader's `ImzMLMeta`, already ported in full |
| `const MSImagingGeometry& getGeometry() const` | `geometry` | Returns `&ImagingGeometry`; see the support-type table |
| `MSImagingGeometry& getGeometry()` | `geometry_mut` | Rust cannot overload on receiver mutability, so the two overloads become two names |
| `uint32_t gridWidth() const noexcept` | `grid_width` | |
| `uint32_t gridHeight() const noexcept` | `grid_height` | |
| `struct Impl` / `std::unique_ptr<Impl> pimpl_` (private) | **not ported: not needed.** | The PIMPL exists so the header needs no expat/Xerces/`imzml::` declaration. Rust module privacy already hides implementation types, and nothing in this module's public signature names a parser type |
| `class MSImagingGeometry;` / `class IonImage;` (forward declarations) | `ImagingGeometry` / `IonImage` | Reproduced here; see below |
| `@see ImzMLFile, ImzMLMeta, ImzMLSpectrumIndex, OnDiscMSExperiment` | two of four resolve | `ImzMLMeta` and `ImzMLSpectrumIndex` are intra-doc links. `FORMAT/ImzMLFile.h` is unported, so the rustdoc names `ImzMLHandler` as the layer actually used; `KERNEL/OnDiscMSExperiment.h` is unported, so it names `format::indexed_mzml_handler` as the nearest analogue. Both gaps are stated at the item rather than left as dangling links, which `cargo doc -D warnings` would reject |

Members with no C++ counterpart, added because the source logs where this crate
returns:

| Rust | Why |
|---|---|
| `geometry_report() -> &GeometryReport` | The header documents the out-of-grid, `< 1` and duplicate-pixel cases as *warnings that still load*. `OPENMS_LOG_WARN` has no Rust equivalent here, so the verdict is retained |
| `uuid_status() -> Option<&UuidStatus>` | Same for the header's "UUID-header mismatch … reported as warnings". Mirrors `ImzMLHandler::uuid_status` |
| `imzml_path()`, `ibd_path()` | The source keeps `ibd_path_` privately and surfaces it only as `getImzMLMeta().ibd_file_path`; both are exposed because both appear in error messages |
| `is_empty()` | Clippy requires it next to `len()` |
| `decoded_spectrum()` | See `getSpectrum` above |
| `build_imaging_geometry()` (free function) | `ImzMLFile::buildImagingGeometry`, which `open()` calls; see below |

## API mapping — the support types

### `ImzMLFile::buildImagingGeometry(const std::vector<ImzMLSpectrumIndex>&, const ImzMLMeta&, MSImagingGeometry&)`

`FORMAT/ImzMLFile.h` is not this package's header, but `open()` calls this
static and the source's own comment calls it the source-of-truth builder shared
by the in-memory and on-disc paths. Ported in full as the free function
`build_imaging_geometry(index, meta) -> Result<(ImagingGeometry, GeometryReport)>`:
the out-parameter becomes a return value and the four `OPENMS_LOG_WARN`
branches become the report's counts and capped lists. The sibling overload
`buildImagingGeometry(const MSExperiment&, MSImagingGeometry&)`, which reads
`imzml:x` / `imzml:y` meta values instead of the index, is **not ported**: it
serves the in-memory `MSImagingExperiment` path, which is unported.

### `IMAGING/MSImagingGeometry.h` → `ImagingGeometry`

| C++ member | Rust | Notes |
|---|---|---|
| `struct Pixel { UInt x, y; Size spectrum_index; }` | `ImagingPixel` | Public fields, same names |
| `void setDimensions(UInt, UInt)` | `set_dimensions` | Now `Result<()>`: the product is checked against `MAX_IMAGE_PIXELS` |
| `UInt getWidth() const` | `width` | |
| `UInt getHeight() const` | `height` | |
| `void setPixelSize(double, double, const std::string& unit = "micrometer")` | `set_pixel_size` | Rust has no default arguments, so `unit` is explicit. Now `Result<()>`: non-finite extents are rejected |
| `double getPixelSizeX() const` | `pixel_size_x` | |
| `double getPixelSizeY() const` | `pixel_size_y` | |
| `const std::string& getPixelSizeUnit() const` | `pixel_size_unit` | `&str` |
| `void addPixel(UInt, UInt, Size)` | `add_pixel` | `Exception::InvalidValue` becomes `Error::InvalidValue`, plus the `MAX_IMAGE_PIXELS` count check |
| `bool hasPixel(UInt, UInt) const` | `has_pixel` | |
| `Size getSpectrumIndex(UInt, UInt) const` | `spectrum_index` | Returns `Option<usize>`: absence is the normal state of most of an image, not an error. The source throws `Exception::ElementNotFound` |
| `const std::vector<Pixel>& getPixels() const` | `pixels` | `&[ImagingPixel]`, insertion order preserved |
| `Size getNumberOfPixels() const` | `number_of_pixels` | Plus `is_empty`, native |
| `void clear()` | `clear` | Restores the documented `1.0`/`1.0`/`"micrometer"` default |
| `static constexpr Size NO_REGION` | `ImagingGeometry::NO_REGION` | Kept, because `add_region` refuses it; no accessor here ever returns it |
| `void addRegion(const MSImagingRegion&)` | `add_region` | Takes the region by value (moved in) rather than by const reference and copying. Plus the `MAX_REGIONS` count check |
| `void removeRegion(Size)` | `remove_region` | `ElementNotFound` becomes `Error::InvalidValue` |
| `void clearRegions()` | `clear_regions` | |
| `const std::vector<MSImagingRegion>& getRegions() const` | `regions` | |
| `const MSImagingRegion& getRegion(Size) const` | `region` | `Result<&ImagingRegion>`; kept as an error because `extract_ion_image_in_region` must tell an unknown region from an empty one. `has_region` answers the same question without an error |
| `Size getNumberOfRegions() const` | `number_of_regions` | |
| `Size regionOf(UInt, UInt) const` | `region_of` | `Option<usize>` instead of the `NO_REGION` sentinel. The source's "acquired pixel first" precondition is kept |
| `std::vector<Size> getRegionPixels(Size) const` | `region_pixels` | `Result<Vec<usize>>`; positions in `pixels()` |
| `std::vector<Size> getRegionSpectrumIndices(Size) const` | `region_spectrum_indices` | `Result<Vec<usize>>`; spectrum indices |
| `static UInt64 packKey_(UInt, UInt)` (private) | private `pack_key` | Layout `(y << 32) | x` kept; see *Native differences* for the container |

Nothing in the header is omitted.

### `IMAGING/MSImagingRegion.h` → `ImagingRegion`

| C++ member | Rust | Notes |
|---|---|---|
| `enum class Shape { Rectangle, Mask }` | `RegionShape` | Same variants |
| `static MSImagingRegion rectangle(Size, const std::string&, UInt, UInt, UInt, UInt)` | `ImagingRegion::rectangle` | `Result<Self>`; `Exception::InvalidValue` for an inverted box becomes `Error::InvalidRange`, which names the same condition more precisely |
| `static MSImagingRegion fromMask(Size, const std::string&, UInt, UInt, UInt, UInt, std::vector<bool>)` | `ImagingRegion::from_mask` | `Result<Self>`; adds the `MAX_IMAGE_PIXELS` cell check and the far-corner overflow check |
| `const std::string& getName() const` | `name` | |
| `bool contains(UInt, UInt) const` | `contains` | |
| `bool intersects(const MSImagingRegion&) const` | `intersects` | |
| `Shape getShape() const` | `shape` | |
| `Size getId() const` | `id` | |
| `UInt getMinX/getMinY/getMaxX/getMaxY() const` | `min_x`, `min_y`, `max_x`, `max_y` | |
| `UInt getBBoxWidth/getBBoxHeight() const` | `bbox_width`, `bbox_height` | |
| `const std::vector<bool>& getMask() const` | `mask` | `&[bool]`; empty for a rectangle, as documented |
| `Size area() const` | `area` | |

Nothing in the header is omitted. `std::vector<bool>` becomes `Vec<bool>`: the
C++ specialisation is bit-packed and Rust's is not, so a mask costs eight times
the memory. `MAX_IMAGE_PIXELS` caps that at 16 MB per mask, and nothing in this
port stores many.

### `IMAGING/IonImage.h` → `IonImage`

| C++ member | Rust | Notes |
|---|---|---|
| `IonImage() = default` | `Default` | Empty 0 x 0 image |
| `IonImage(UInt, UInt)` | `IonImage::new` | `Result<Self>`; the pixel count is checked |
| `void resize(UInt, UInt)` | `resize` | `Result<()>`, and atomic: a refused resize leaves the image unchanged |
| `UInt getWidth() const` | `width` | |
| `UInt getHeight() const` | `height` | |
| `bool hasPixel(UInt, UInt) const` | `has_pixel` | `false` out of bounds, as documented — the one accessor that does not report a range error |
| `double getIntensity(UInt, UInt) const` | `intensity` | `Result<f64>`; `Exception::IndexOverflow` becomes `Error::InvalidValue` |
| `void setIntensity(UInt, UInt, double)` | `set_intensity` | `Result<()>`; same mapping, plus a finiteness check |
| `void setMzRange(const RangeMZ&)` | `set_mz_range` | Takes `kernel::ranges::RangeBase`, the base type `RangeMZ` specialises |
| `const RangeMZ& getMzRange() const` | `mz_range` | Empty when never set, as documented |
| `const std::vector<double>& getData() const` | `data` | `&[f64]`, row-major |
| `const std::vector<bool>& getMask() const` | `mask` | `&[bool]`, indexed identically |
| `Size linearIndex_(UInt, UInt) const` (private) | private `linear_index` | |

Nothing in the header is omitted.

### `IMAGING/IonImageExtraction.h`

| C++ | Rust | Notes |
|---|---|---|
| `template <typename GetSpectrum> IonImage Internal::extractIonImage(const MSImagingGeometry&, double, double, const std::vector<Size>&, Size, GetSpectrum&&)` | private `OnDiscImzMLExperiment::extract` | The template exists to share one kernel between the in-memory and on-disc callers that differ only in how a spectrum is obtained. With the in-memory caller unported there is one caller, so the kernel is a private method rather than a generic free function. The `@note` about the lambda's explicit `-> const MSSpectrum&` return type is a C++ lifetime hazard with no Rust counterpart: `peaks_only` returns an owned value and the borrow checker rejects a dangling alternative. Reintroduce the generic form when the IMAGING package ports `MSImagingExperiment` |
| `const std::vector<Size>& pixel_indices` | `Option<&[usize]>` | `None` means the whole grid, so the source's `std::iota` allocation of one `Size` per pixel is not made |

## Preserved source conventions

1. **Both files, two-file discipline.** The `.ibd` path is derived exactly as
   `ImzMLFile::inferIbdPath_` derives it — a case-insensitive `.imzML` suffix
   replaced by `.ibd`, any other name gaining `.ibd` — through the reader's
   `infer_ibd_path`, and an explicit override is threaded through both the index
   load and the UUID check so neither targets a stale sibling. The upstream
   suite guards that case and so does `an_explicit_ibd_override_is_honoured`.
2. **Lazy peaks, eager grid.** `open()` parses the XML index and builds the
   grid; no peak array is touched until a decode asks. The grid is built at
   `open()` because it is an in-memory pass over an already-parsed index and
   because building it there surfaces a broken coordinate grid at `open()`
   rather than on the first pixel query. Both are the source's stated reasons.
3. **Non-fatal coordinate problems.** A pixel outside the declared
   `IMS:1000042` x `IMS:1000043` grid, a coordinate below 1, a duplicated pixel
   and a `z != 1` plane are all dropped from the grid and never abort the load;
   every such spectrum stays reachable by linear index. Of duplicates the first
   spectrum in document order keeps the pixel.
4. **The declared extent is a floor.** `gridWidth`/`gridHeight` report
   `ImzMLMeta::max_count_x`/`max_count_y`, which the reader raises to the
   largest observed coordinate; when the file declares neither, the grid derives
   `max_x + 1` / `max_y + 1` after the pixel loop. The two can therefore differ
   from `geometry().width()`/`height()`, and the port keeps that.
5. **Zero-based grid, one-based coordinates.** `getSpectrumAtCoord` takes
   imzML-native 1-based `(x, y, z)` and the grid stores 0-based `(x-1, y-1)`;
   only `z == 1` is addressable, because the grid is two-dimensional.
6. **Sorted peaks, aligned arrays.** `getSpectrum` attaches the auxiliary
   arrays first and sorts by m/z last, so the float data arrays stay aligned
   with peak order. `sort_by_position` reorders every parallel array, which is
   what makes the ordering safe.
7. **Inclusive extraction window.** `dm = mz * tolerance_ppm * 1e-6` and the sum
   runs from the first peak at or above `mz - dm` to the last at or below
   `mz + dm`, which is what `MZBegin`/`MZEnd` mean. A pixel with a spectrum but
   no peak in the window is valid with intensity 0; a pixel absent from the grid
   stays invalid. The image's m/z range records the window.
8. **Extraction ignores auxiliary arrays.** The sweep decodes only m/z and
   intensity, so a compressed auxiliary array fails `spectrum` and not
   `extract_ion_image`. The source splits `decodeSpectrum` from `decodePeaks`
   for exactly this and documents the asymmetry; the port keeps both halves.
9. **`close()` keeps the index.** It releases the `.ibd` and clears the grid,
   and `len`, `index` and `meta` keep answering afterwards.
10. **Compression is only `MS:1000576`.** Any other child of `MS:1000572` —
    zlib and every numpress variant — is refused without being enumerated. That
    check lives in the reader and is inherited unchanged.
11. **Serial.** Nothing in this header parallelises in C++ and nothing does here.
12. **Regions are a decoupled, disjoint overlay.** Membership is derived on
    demand from the footprint; no per-pixel region state is stored, and
    `add_region` refuses an overlap, so `region_of` has at most one answer.

## Native differences

Each is a deliberate divergence, documented at the item in rustdoc as well.

1. **`open()` is atomic.** The source replaces its `Impl` on the first line of
   `open()`, so a failed open discards whatever dataset was loaded. This builds
   into locals and commits only on success, per the crate's atomicity rule;
   `a_failed_open_leaves_the_previous_dataset_in_place` pins it.
2. **Warnings are returned, not logged.** `GeometryReport` carries the
   out-of-grid, `< 1`-coordinate, duplicate-pixel and other-plane counts with
   the first `MAX_LISTED_PROBLEMS` spectrum indices of each — the source's own
   cap of 20 — and `uuid_status()` carries the UUID verdict. Neither rejects
   anything, so the source's "still loads" policy is unchanged; a caller that
   wants strict conformance can now enforce it. This is the same choice
   `ImzMLHandler::uuid_status` made one layer down.
3. **Absence is not an error where it is normal.**
   `ImagingGeometry::spectrum_index` returns `Option`, and `region_of` returns
   `Option` rather than the `NO_REGION` sentinel. `region`, `region_pixels` and
   `region_spectrum_indices` keep the error, because an unknown region is a
   caller mistake rather than an empty answer.
4. **`operator[]` has no `Index` impl.** See the mapping table.
5. **Move-only semantics.** The source deletes copy construction and copy
   assignment and defines move construction and move assignment, each leaving
   the moved-from object holding a *fresh empty `Impl`* so that calling into it
   is defined rather than a null dereference. Rust needs none of that: a move is
   the default, a moved-from binding cannot be used at all, and there is nothing
   to reset. The type is deliberately **not `Clone`** — it owns an open `File`
   with one seek position, and two handles sharing it would interleave seeks.
   `meta()` and `geometry()` return clonable values for callers who need a copy.
6. **Every decode takes `&mut self`.** The source's methods are `const` and
   mutate the `FILE*` through the pointer, which is why the header carries a
   `@note` that the class is not thread-safe. Here that is a compile-time
   property.
7. **`close()` copies the index.** Keeping the source's post-close contract
   (`getNrSpectra` and `getIndex` keep answering) while actually releasing the
   file means moving the index out of the reader that owned it, so `close()` is
   O(spectra) where the source's is O(1).
8. **`IMPeakType` is not set.** `decodeSpectrum` marks a spectrum carrying
   per-peak ion mobility as `IMPeakType::IM_PROFILE`; `MSSpectrum` has no such
   field and this package may not add one. The same gap is recorded in
   `docs/SPECTRUM_MOBILITY_SUPPORT.md`. Observable consequence: a caller cannot
   distinguish profile from centroided mobility on an on-disc spectrum.
9. **Deterministic lookup containers.** `MSImagingGeometry` uses
   `std::unordered_map` for both its pixel key map and its region-id map;
   `BTreeMap` is used here, per the crate's determinism convention. The key
   layout `(y << 32) | x` is unchanged. Neither map is iterated in the source,
   so no observable order changes.
10. **Non-finite values are refused.** `set_pixel_size` and
    `IonImage::set_intensity` reject NaN and infinity; the source stores any
    `double`. A `.ibd` holding a non-finite intensity therefore produces an
    error rather than a pixel that compares unequal to itself.
11. **The far corner of a mask cannot wrap.** `MSImagingRegion::fromMask`
    computes `origin_x + width - 1` in wrapping `UInt` arithmetic, which for a
    large origin yields a box starting above where it ends. `from_mask` refuses
    it with `Error::InvalidRange`.
12. **The `iota` is not allocated.** The whole-image overload iterates the pixel
    range instead of materialising one `Size` per pixel, as the source does.
13. **One unreachable source error has no counterpart.**
    `decodePeaksInto_` raises `ParseError("failed to decode spectrum peaks from
    .ibd")` when a declared non-zero array decodes to nothing. The reader
    preflights every range against the `.ibd` length and reports a short read as
    an error, so a non-zero count cannot silently yield an empty array.
14. **`PeakFileOptions` and `ProgressLogger` are not threaded through.** The
    source's `open()` constructs an `ImzMLFile`, sets `ProgressLogger::NONE` and
    inherits that object's options. This façade reads no option and reports no
    progress. `SortSpectraByMZ` is the only option that would change an on-disc
    result, and this port always sorts, which is what the source does on this
    path.
15. **An inherited nuance of the `.ibd` derivation.** The reader's
    `infer_ibd_path` uses `Path::set_extension`, which for a file named exactly
    `.imzML` — a leading-dot name Rust treats as extensionless — yields
    `.imzML.ibd` where C++'s six-character strip yields `.ibd`. That function
    belongs to `src/format/imzml_handler.rs`; no real dataset is named that way.

## Checked boundaries and evidence

### Bounds

The pixel count is a product of two dimensions read from the file, and every
offset the sweep follows came out of the file and indexes into a second file.
Nothing is allocated before it is checked.

| Ceiling | Value | Guards | Source counterpart |
|---|---|---|---|
| `MAX_IMAGE_PIXELS` | 16,777,216 (4096 x 4096) | `ImagingGeometry::set_dimensions`, `add_pixel`, `IonImage::new`/`resize`, `ImagingRegion::from_mask` | None. `IonImage::resize` allocates `width * height` doubles plus a mask from the file's declared grid. The value matches `MSSpectrum::MAX_RASTER_PIXELS`, which bounds the same product-of-two-file-values shape, and sits far above the reader's own 5,000,000-spectrum ceiling |
| `MAX_REGIONS` | 4,096 | `ImagingGeometry::add_region` | None. Insertion tests the newcomer against every existing region |
| `MAX_LISTED_PROBLEMS` | 20 | `GeometryReport` list growth | `const Size max_listed = 20` in `buildImagingGeometry` |
| `ImzMLReadLimits` | reader's | every `.ibd` range, preflighted against the ceilings *and* the actual file length before any allocation | `MAX_IBD_ARRAY_ELEMENTS` only, with no file-length comparison |

Consequences, all tested:

- A file declaring a grid above the ceiling fails `open()` with
  `Error::InvalidValue` and commits nothing, where the source opens and then
  fails inside the first extraction. A legitimate dataset cannot reach this: it
  would need more pixels than the reader will index spectra.
- A hostile `IMS:1000102` offset (`u64::MAX - 8`) parses into the index and is
  refused at decode, before any allocation, on both the `spectrum` and the
  `extract_ion_image` path.
- A tightened `max_array_elements` reaches the `.ibd` through
  `open_with_limits` and refuses the decode while still allowing the index to
  be read.
- The region-intersection walk is bounded: two rectangles short-circuit, and
  any pair involving a mask walks at most the smaller bounding box, which for a
  mask *is* its mask and is therefore below `MAX_IMAGE_PIXELS`.
- Overflow is checked, not wrapped: the derived `max_x + 1`, a mask's
  `origin + extent - 1` and `u64 -> usize` grid offsets all use checked
  arithmetic.

### Error mapping

| Source exception | Rust |
|---|---|
| `Exception::FileNotFound` (`.ibd` not open) | `Error::Io` with `ErrorKind::NotFound`, naming the path |
| `Exception::FileNotFound` (file missing at `open()`) | `Error::Io`, from the reader |
| `Exception::ParseError` (malformed XML, bad IMS value, missing pixel coordinate) | `Error::Parse`, from the reader |
| `Exception::ParseError` (m/z / intensity length mismatch) | `Error::Parse`, message naming the `(x,y,z)` pixel as the source's does |
| `Exception::ParseError` (compressed external array) | `Error::Unsupported`, from the reader |
| `Exception::IndexOverflow` (spectrum index) | `Error::InvalidValue` |
| `Exception::IndexOverflow` (`IonImage` coordinate) | `Error::InvalidValue` |
| `Exception::ElementNotFound` (no spectrum at a coordinate) | `Error::InvalidValue`, matching `ImzMLHandler::spectrum_at_coord` |
| `Exception::ElementNotFound` (unknown region id) | `Error::InvalidValue` |
| `Exception::InvalidValue` (negative or non-finite `mz`/`tolerance_ppm`; pixel references a missing spectrum; duplicate or out-of-bounds pixel; duplicate, overlapping or `NO_REGION`-identified region) | `Error::InvalidValue` |
| `Exception::InvalidValue` (inverted region box) | `Error::InvalidRange` |

`ElementNotFound` maps to `InvalidValue` rather than `MissingInformation` for
consistency with the imzML reader, which made that choice first for the same
condition.

### Evidence

**Tier 3, source review.** `OnDiscImzMLExperiment` has no class test of its
own; it is constructed nineteen times in `ImzMLFile_test.cpp` and three times in
`ImzMLFile_all_modes_test.cpp`. The literals transcribed into
`tests/on_disc_imzml_experiment.rs`:

| Upstream section | Literal | Rust test |
|---|---|---|
| `mode 5: OnDiscImzMLExperiment random access` | `isOpen()`, 9 spectra, 3 x 3 grid, mode `continuous`, coordinate and index access agree | `the_continuous_fixture_opens_with_the_upstream_grid` |
| `const MSImagingGeometry& getGeometry() const` | `getWidth()==gridWidth()`, `getHeight()==gridHeight()`, `getNumberOfPixels()==getNrSpectra()`, `getSpectrumIndex(0,0)==0` | `the_grid_matches_the_dataset_and_maps_every_spectrum` |
| `mode 7: buildImagingGeometry` (the `MSExperiment` overload, not ported) | the grid maps 0-based `(2,2)` to the spectrum acquired at imzML `(3,3)`, which is index 8 | `the_grid_matches_the_dataset_and_maps_every_spectrum` |
| `processed imzML encoding` | mode `processed`, non-empty spectrum 0 | `the_processed_fixture_opens_and_decodes` |
| `[EXTRA] … honours an explicit .ibd path override` | a relocated `.imzML` opens only with the override | `an_explicit_ibd_override_is_honoured` |
| `OnDiscImzMLExperiment tolerates duplicate pixel coordinates by default` | `size()==2`, one pixel, `getSpectrumIndex(0,0)==0` | `a_duplicate_pixel_keeps_the_first_spectrum` |
| `void load sorts external peaks when getSortSpectraByMZ is true` | on-disc peaks sorted 121.0 / 131.0; `disc_img.getIntensity(0,0)==1310.0` | `the_extraction_is_faithful_to_the_upstream_sorted_peak_case` |
| `IonImage extractIonImage(double, double)` | `Exception::InvalidValue` for a negative `mz` and for a negative tolerance | `a_negative_or_non_finite_extraction_argument_is_refused` |
| `IonImage extractIonImage(double, double, Size)` | region image keeps the grid's dimensions; unknown region id throws | `a_region_extraction_covers_only_the_region` |
| `MSImagingGeometry_test.cpp` default/`addPixel`/`clear` sections | `1.0`/`1.0`/`"micrometer"`, insertion-order pixels, duplicate and bounds errors | `a_default_geometry_has_unit_micrometer_pixels`, `pixels_keep_insertion_order_and_reject_duplicates_and_strays`, `clear_restores_the_default_geometry` |
| `IonImage_test.cpp` | zeroed and masked-out fresh image, row-major indexing, out-of-bounds behaviour | `a_fresh_image_is_zeroed_and_fully_masked_out`, `writing_a_cell_marks_it_valid_and_is_row_major`, `image_coordinates_are_bounds_checked` |
| `MSImagingRegion_test.cpp` | inclusive boxes, mask validation, intersection | `a_rectangle_covers_its_inclusive_bounding_box`, `a_mask_covers_only_its_set_bits`, `intersection_is_symmetric_and_mask_aware` |

**One independent cross-check.** The C++ `extractIonImage` sections compare the
on-disc image against the in-memory `MSImagingExperiment` path pixel for pixel.
That path is unported, so `a_whole_image_extraction_sums_the_window_at_every_pixel`
and `a_region_extraction_covers_only_the_region` instead compare every pixel
against a sum computed in the test directly from `ImzMLHandler::spectrum` over
the same window — a second, independent summation over the same `.ibd` bytes,
which is the same shape of check the C++ section performs.

**Tier 4, independently derived.** `MAX_IMAGE_PIXELS`, `MAX_REGIONS`, the
hostile offset, the atomic-open behaviour, the non-finite rejections, the
`u32` far-corner overflow, the `close()`-retains-index contract and every
`Error`-variant choice. No upstream fixture reaches any of them.

**Not tier 1.** No C++ was built or executed and no C++ output was retained.

### Source defect recorded

`OnDiscImzMLExperiment::open` (`OnDiscImzMLExperiment.cpp:238-246`) re-derives
the `.ibd` path from `pimpl_->meta_.ibd_file_path` and falls back to its own
copy of `ImzMLFile::inferIbdPath_` when that string is empty. It never is:
`ImzMLFile::loadImpl_` (`ImzMLFile.cpp:568-570`) resolves the path as
override-or-`inferIbdPath_` and assigns `meta_.ibd_file_path` unconditionally
before parsing. The fallback is therefore unreachable duplicated logic, and a
future change to `inferIbdPath_` would silently leave this copy behind. Not a
behavioural bug today. The port has one derivation, in the reader.

## Section accounting

`OnDiscImzMLExperiment` is constructed in 22 distinct upstream sections: 19 in
`ImzMLFile_test.cpp` and 3 in `ImzMLFile_all_modes_test.cpp`. Every one is
accounted for below — 17 ported, 2 partial, 3 not ported. None is unaccounted.

| Upstream section | Status | Rust |
|---|---|---|
| `OnDiscImzMLExperiment random access` | ported | `the_continuous_fixture_opens_with_the_upstream_grid` |
| `const MSImagingGeometry& getGeometry() const` | ported | `the_grid_matches_the_dataset_and_maps_every_spectrum` |
| `IonImage extractIonImage(double mz, double tolerance_ppm) const` | ported | `a_whole_image_extraction_sums_the_window_at_every_pixel`, `a_negative_or_non_finite_extraction_argument_is_refused`. The section's in-memory comparison arm is replaced by an independent sum; see *Evidence* |
| `[EXTRA] OnDiscImzMLExperiment::open honours an explicit .ibd path override` | ported | `an_explicit_ibd_override_is_honoured` |
| `void store round-trip FloatDataArray ion mobility and non-standard` | **partial** | `an_auxiliary_array_reaches_the_decoded_spectrum` covers the on-disc half — two auxiliary arrays indexed, the ontology name, the `unit_accession`, `containsIMData()` — except the `getIMPeakType() == IM_PROFILE` assertion, which `MSSpectrum` cannot represent (native difference 8). The round trip itself needs `ImzMLWriter`, which is unported |
| `void load and OnDisc reject zlib-compressed external m/z and intensity` | ported | `a_compressed_peak_array_is_refused_on_both_paths` (`MS:1000574`) |
| `OnDiscImzMLExperiment rejects compressed zero-length aux arrays` | ported | `a_compressed_auxiliary_array_fails_a_decode_but_not_an_extraction` |
| `OnDiscImzMLExperiment tolerates duplicate pixel coordinates by default` | ported | `a_duplicate_pixel_keeps_the_first_spectrum` |
| `IonImage extractIonImage(double mz, double tolerance_ppm, Size region_id) const` | ported | `a_region_extraction_covers_only_the_region` |
| `void load sorts external peaks when getSortSpectraByMZ is true` | ported | `the_extraction_is_faithful_to_the_upstream_sorted_peak_case` covers the on-disc half, including the 1310.0 pixel |
| `void load and OnDisc reject numpress-compressed external m/z` | ported | `a_compressed_peak_array_is_refused_on_both_paths` (`MS:1002312`) |
| `void load skips one bad aux length and still returns all spectra` | ported | `an_auxiliary_array_of_the_wrong_length_is_skipped` |
| `void load and OnDisc drop zero-length aux without a ghost IM array` | ported | `a_zero_length_auxiliary_array_leaves_no_array` |
| `void store skips unnamed FloatDataArray` | **not ported: writer-driven.** | Its on-disc half asserts `getIndex(0).aux.size() == 0` for a file `ImzMLWriter` produced, and `ImzMLWriter` is unported. The reader's unnamed-array rule is `AuxSkipReason::Unnamed` and is tested at that level in `tests/imzml_handler.rs` |
| `void load drops phantom IntegerDataArray for integer-typed aux` | ported | `an_integer_typed_auxiliary_array_becomes_a_float_array` |
| `void load and OnDisc drop inline aux when peaks are external` | ported | `an_inline_auxiliary_array_is_not_decoded` |
| `void store skips FloatDataArrays named after the peak arrays` | **not ported: writer-driven.** | Its on-disc half asserts that the writer emitted no shadowing auxiliary array; nothing about it is a façade behaviour |
| `void store drops integer and string data arrays` | **not ported: writer-driven.** | Its on-disc half asserts empty integer and string arrays on a decoded spectrum, which `an_integer_typed_auxiliary_array_becomes_a_float_array` and `an_inline_auxiliary_array_is_not_decoded` already establish |
| `void load and OnDisc skip aux array without a supported binary data type` | ported | `an_auxiliary_array_without_a_supported_type_is_skipped` |
| `mode 5: OnDiscImzMLExperiment random access` | ported | `the_continuous_fixture_opens_with_the_upstream_grid` |
| `cross-mode consistency: full load vs on-disc vs RAM lookup` | **partial** | `every_pixel_of_the_upstream_grid_is_addressable` covers the on-disc arm and its agreement with the index. The full-load and in-RAM arms need `ImzMLFile::load` and `MSImagingExperiment`, both unported |
| `processed imzML encoding` | ported | `the_processed_fixture_opens_and_decodes` |

## Deferred

- **`IMAGING/MSImagingGeometry.h`, `IMAGING/MSImagingRegion.h`,
  `IMAGING/IonImage.h` and `IMAGING/IonImageExtraction.h`** stay unowned. The
  tables above record exactly which members were reproduced here and why; their
  class tests were read but not ported section by section, and the ledger status
  for them is untouched.
- **`IMAGING/MSImagingExperiment.h`**, the in-memory counterpart that shares the
  geometry, is unported, as is the `buildImagingGeometry(const MSExperiment&, …)`
  overload that serves it. Porting it should also restore the shared generic
  extraction kernel.
- **`MSSpectrum::im_peak_type`** — see native difference 8.
- **`PeakFileOptions` / `ProgressLogger`** — see native difference 14.
- **`docs/core-sdk-reviewed-apis.json` and `docs/core-sdk-coverage.json`** were
  not edited; both are outside this package's scope. `tools/core_sdk_coverage.py`
  already reports the coverage file stale from the imzML reader package's
  manifest alone, and its heuristic already maps
  `src/kernel/on_disc_imzml_experiment.rs` onto `KERNEL/OnDiscImzMLExperiment.h`
  and `IMAGING/IonImage.h`. The integrator adds the ledger entry and runs
  `--write`.
- **`SOURCE_PROVENANCE.json`'s `current_sdk_reference_manifests`** lists neither
  `tests/data/imzml_handler_provenance.json` nor
  `tests/data/on_disc_imzml_provenance.json`, so `check_core_sdk.py` does not yet
  hash-verify either. Adding both is the integrator's step.
- **`docs/doc-coverage.json`** was not rewritten.
  `src/kernel/on_disc_imzml_experiment.rs` measures 100.0% (88/88) and
  `src/kernel.rs` stays at 100.0% (70/70).
- **CI wiring**: append `--test on_disc_imzml_experiment` to a `--features mzml`
  line of the `minimum-rust` job in `.github/workflows/rust.yml`. The test passes
  under `--no-default-features --features mzml`, under `--all-features` and
  under `cargo +1.85.0`.
