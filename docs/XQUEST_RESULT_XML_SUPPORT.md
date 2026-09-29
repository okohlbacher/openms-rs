# XQuestResultXMLFile and XQuestResultXMLHandler

`FORMAT/XQuestResultXMLFile.h` and `FORMAT/HANDLERS/XQuestResultXMLHandler.h`
(core `bc9cc12`), ported by lane F3a in wave 4 as
`src/format/xquest_result_xml.rs` (feature `idxml`). xQuest result files
(`.xquest.xml`) hold cross-linking (XL-MS) identifications written by xQuest
and OpenPepXL. TOPP consumers: `IDFileConverter` (in and out), `XFDR` (in and
out), `OpenPepXL` and `OpenPepXLLF` (out).

Status: **partial**. Loading, storing, the getters and the handler's string
helpers are ported and differentially tested against the Release build; the
two `writeXQuestXMLSpec` overloads (the spec.xml side file) are not, because
their parameter types, `OPXLDataStructs::CrossLinkSpectrumMatch` and
`OPXLDataStructs::PreprocessedPairSpectra`, have no port yet (dependency:
`ANALYSIS/XLMS/OPXLDataStructs.h`).

## API mapping

### `XQuestResultXMLFile`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `XQuestResultXMLFile()` | `XQuestResultXmlFile::new`, `Default` | `n_hits_ = -1`; the score range is 0 (uninitialized in the source) |
| base `XMLFile("/SCHEMAS/xQuest_1_0.xsd", "1.0")` | `SCHEMA_LOCATION`, `SCHEMA_VERSION` | the schema is not validated against, as in the source |
| base `ProgressLogger` | not ported | the source reports no progress from these members |
| `void load(filename, pep_ids, prot_ids)` | `load`, `load_with_limits`, `read`, `read_with_limits` | returns `XQuestDocument` (`protein_identifications`, `peptide_identifications`, `log`); applies the five `OPXLHelper` steps (below) |
| `void store(filename, poid, peid) const` | `store` | extension check (`Exception::UnableToCreateFile` → `Error::InvalidValue`); renders fully before creating the file |
| `int getNumberOfHits() const` | `number_of_hits` | |
| `double getMinScore() const`, `double getMaxScore() const` | `min_score`, `max_score` | |
| `static void writeXQuestXMLSpec(...)` (labelled linkers) | **not ported** | needs `OPXLDataStructs::PreprocessedPairSpectra`, `CrossLinkSpectrumMatch` |
| `static void writeXQuestXMLSpec(...)` (label-free linkers) | **not ported** | needs `OPXLDataStructs::CrossLinkSpectrumMatch` |
| private `getxQuestBase64EncodedSpectrum_`, `wrap_` | **not ported** | only reached from `writeXQuestXMLSpec` |
| private `n_hits_`, `min_score_`, `max_score_` | private fields | |

### `Internal::XQuestResultXMLHandler`

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `static std::map<Size, std::string> enzymes` | `ENZYMES` | an unlisted number maps to the empty name, as `operator[]` does (the source also inserts it into the process-wide map) |
| `static std::map<std::string, UInt> months` | `MONTHS` | an unlisted abbreviation maps to 0 |
| reading constructor + `onStartElement` / `onEndElement` | `XQuestResultXmlHandler::parse` → `ParsedXQuest` | driven by `format::sax_xml` (`SaxEncoding::Declared`) |
| writing constructor + `writeTo(os)` | `XQuestResultXmlHandler::write_to(writer, proteins, peptides, date)` | `date` replaces `DateTime::now()` so output is reproducible; `store` passes the current local time |
| `getMinScore`, `getMaxScore`, `getNumberOfHits` | `ParsedXQuest::{min_score, max_score, n_hits}` | |
| `static StringList splitByNth(input, separator, n)` | `XQuestResultXmlHandler::split_by_nth` | returns `[String; 2]` |
| `static StringList splitByMiddle(input, separator)` | `XQuestResultXmlHandler::split_by_middle` | |
| private `extractDateTime_`, `addMetaValues_`, `getLinkPosition_`, `setPeptideEvidence_` | private `Reader::date_time`, the shared-meta merge, `link_position`, `Reader::evidences` | |
| private members (`decoy_string_`, `spectrum_index_light_`, …, `peptide_id_meta_values_`) | fields of the private `Reader` | `ms_run_path_` is collected but never used by the source, and not kept |
| (none) | `ReadLimits`, `XQuestDocument::suppressed_log_lines` | resource ceilings |

### `OPXLHelper` steps inside `load`

`ANALYSIS/XLMS/OPXLHelper.h` has no port. The five members `load` calls are
private functions of this module, statement for statement:
`addXLTargetDecoyMV`, `addBetaAccessions`, `combineTopRanksFromPairs`,
`removeBetaPeptideHits`, `computeDeltaScores`. The mzIdentML reader keeps
copies that diverge on purpose (`src/format/mzidentml.rs`); a later
`OPXLHelper` port should take these faithful ones.

## Preserved source conventions

- **OpenPepXL or xQuest.** `xquest_version` containing `OpenPepXL` selects the
  OpenPepXL flavour: `yyyy-MM-dd hh:mm:ss` dates, `enzyme_name`,
  `tolerancemeasure_ms1`, `scan_index_light`/`_heavy`, `charges` and a required
  (then ignored) `run_path`. Otherwise xQuest: `Fri Dec 18 12:28:42 2015`
  dates (an empty third field for a single-digit day), `enzyme_num`,
  `tolerancemeasure`, spectrum indices split out of the `spectrum` name with
  `splitByMiddle(spectrum, '_')`, `X` in `seq1`/`seq2` read as
  `M(Oxidation)`, and the searched charges collected from the results.
- **Enzymes** are looked up as `ProteaseDB::getEnzyme`: by name, lower-case
  name or synonym (`trypsin` → `Trypsin`); an unknown name is
  `ElementNotFound` (`x08`).
- **Variable modifications** are never found: the source reads the xQuest
  `<letter>,<mass>` form's mass through `double(DataValue(std::string))`,
  which returns the bits of the string's storage pointer, not the number, and
  only one-letter forms are looked up at all. The list is therefore empty
  whenever `variable_mod` is present (`x01`, `x03`). *C++ defect candidate.*
- **Mono-links** get the cross-linker's modification on the linked residue:
  of the registry's modifications within 0.01 Da of `xlinkermass` on that
  residue (`searchModificationsByDiffMonoMass`, registry order), the *last*
  whose id contains the cross-linker name; without one, the first, announced
  on standard output as `Set default mono-link: <id>` (`XQuestDocument::log`,
  `x02`).
- **Terminal link specificity** compares the link's residue letter with the
  cross-linker's residue lists, indexing the `seq1`/`seq2` *text* (modification
  annotations included), as the source does (`x09`).
- xQuest's U+2212 minus in `xlinkermw` and `xlinkermass` reads as `-` (`x10`).
- `decoy_prefix` is read as `istringstream >> bool` (`"0"` false, any other
  leading integer true, no integer false), absent means true; `decoy_string`
  replaces the default `decoy_` when present and not empty (`x15`).
- The alpha hit carries `xl_chain = MS:1002509`, the beta `MS:1002510`; both
  are removed by `removeBetaPeptideHits`, which folds the beta hit into
  `BetaPepEv:pre/post/start/end` (decimal positions, `-1` for unknown) and merges
  identifications by the byte order of their `spectrum_reference` text (`7,9`
  after `12`, `x01`).
- `combineTopRanksFromPairs` groups by the byte order of `spectrum_index`,
  sorts each group by descending alpha score as the Release build's
  `std::sort` over reverse iterators leaves it (`math::source_sort`; `x12` has
  21 tied identifications), drops consecutive equal identifications
  (`PeptideIdentification::operator==`, doubles within 1e-6: `x11`) and
  renumbers `xl_rank`.
- **Writing** follows `writeTo` byte for byte: doubles in the stream's
  15-significant-digit `%g` form or through `StringUtils::toStr` (15 fraction
  digits) exactly where the source uses each; list meta values with their
  brackets stripped (`monolinkmw`, `AArequired1`); attribute values unescaped;
  a new `<spectrum_search>` unless the hit's `spectrum_reference` is the
  *string* equal to the previous one (a loaded file stores integers, so every
  hit opens its own element, as in the recorded outputs); the monoisotopic
  weight as `AASequence::getMonoWeight` sums it (free-residue formula plus
  delta formula, minus water, per residue).

## Native differences

- **Meta value order.** The meta values written at the end of each
  `<search_hit>` come in byte order of their names; the source writes them in
  `MetaInfoRegistry` order, i.e. the order in which its process first
  registered each name, which depends on everything the process did before.
  The differential test compares that segment as a set.
- **Uninitialized source members.** `min_precursor_charge_` and
  `max_precursor_charge_` (xQuest flavour) start at 0 here; the source never
  initializes them, and the Release build read 0 for `x03` but 1 for `x04`, so
  the reported `precursor:min_charge` of an xQuest file is undefined there and
  0 here. Indices of a `<search_hit>` before any `<spectrum_search>` are 0.
- **Undefined indexing is refused.** Where the source indexes past a split
  (`rtsecscans`/`mzscans` without `:`, an empty `xlinkposition`, `charges`, a
  short date, a `spectrum` name without file part or indices), the port
  returns `Error::InvalidRange`.
- **Non-finite numbers.** The metadata model holds finite doubles only: a
  `score`, mass or tolerance of `nan`/`inf` is refused, and a delta score that
  is not finite (a zero score followed by a nonzero one) is left unset where
  the source stores the infinity.
- **Writer refusals.** No protein identification, a hit without a peptide
  evidence, or a string meta value where the source converts to a number
  (`double(DataValue)` reads the union's storage) are `Error::InvalidValue`
  rather than undefined behaviour.
- A `std::exception` for an unsupported link `type` becomes
  `Error::InvalidValue` carrying the source's log message (`x06`).

## Evidence

- Tier 1, executed differential: `../oracle/F3a-xquest` drove the Linux
  x86_64 Release libOpenMS (`openms4-release-bc9cc12-c19e494-174b576`) on
  ibminode06, twice, identically, over the class-test fixture, the upstream
  TOPP files `OpenPepXL_output.xquest.xml`, `XFDR_test_in2/3.xquest.xml`,
  `XFDR_test_out{1,2,3,4,5,7}.{xquest.xml,idXML}` and 16 generated cases,
  recording `load`, `store` and idXML → xquest.xml conversion, plus 14
  `splitByNth`/`splitByMiddle` calls. `tests/xquest_result_xml.rs` compares
  every field line by line. The six XFDR files over 700 kB are not retained in
  the repository; they passed the same tests in an ad-hoc run on dax.
- Tier 3: the literals of `XQuestResultXMLFile_test.cpp`
  (`class_test_store_and_reload`).
- Tier 4: resource ceilings and refusals.
