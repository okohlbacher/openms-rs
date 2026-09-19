# FileInfo: the consensusXML, identification and FASTA branches (A7)

Package **A7-FILEINFO** (early TOPP bundle, wave 8). It ports the three content
branches of `OpenMS::FileInfo` that A4, A5 and A6 left open, with their `-m`,
`-p` and `-s` arms:

| Branch | Source | Rust |
| --- | --- | --- |
| FASTA | `FORMAT/FileInfo.cpp:853-1076`, `:2001-2004`, `:2112-2114`, `:2377-2379` | `src/format/file_info/fasta.rs` |
| consensusXML | `:1146-1311`, `:1985-1989`, `:2101-2104`, `:2257-2372` | `src/format/file_info/consensus.rs` |
| idXML, mzIdentML | `:1312-1470`, `:1990-1996`, `:2105-2107`, `:2373-2376` | `src/format/file_info/identifications.rs` |

The TOPP tool (`OpenMS4-topp/src/FileInfo.cpp`, topp `174b576`) is a thin
wrapper; every line below comes from the library class at the core pin
`bc9cc12`. After this package `FileInfo.h` still has pepXML, mzTab, trafoXML and
PQP open (A8), together with `-v`.

Evidence and hashes: `tests/data/file_info_a7_provenance.json`.
Tests: `tests/file_info_a7.rs`, plus the `std::hash` unit tests in
`src/format/file_info/fasta.rs`. Two older files change with the scope:
`tests/file_info.rs` runs the consensusXML and FASTA class-test sections that
A4 had to leave as refusals, and `tests/topp_file_info.rs` reproduces
TOPP_FileInfo_7, _10, _13, _17, _18 and _20 through FuzzyDiff against the
retained upstream outputs instead of listing them as not ported. Each keeps a
tripwire that fails if a branch goes back to refusing.
Oracle: `../oracle/a7-fileinfo`, 60 cases against the Release C++ FileInfo of
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
ibminode06, run twice and reproduced. 40 of them have both reports compared
byte for byte in the test.

---

## 1. What each branch writes

### FASTA (`:853-1076`)

The entries load through `FASTAFile::load`. Before anything is counted, the
whole file is classified: `:900-912` walks every byte of every sequence and one
byte outside `"ACGTUNacgtunRYSWKMBDHVryswkmbdhv"` makes the file amino acid.
That single decision picks the labels (`nucleotide` / `amino acid`), the
alphabet the per-sequence ambiguity test uses, and which pair of totals ends the
block.

The report is, in order: the sequence count; the five-line length distribution,
but only when there is at least one sequence; the count of sequences carrying at
least one ambiguous residue; the duplicate-header and duplicate-sequence counts;
the residue total; one line per residue byte; and the two ambiguity totals.
Every percentage is `Math::percentOf(value, entries.size(), 2)`, which answers
`0.0` for a zero total rather than dividing.

Below three sequences the source does not ask for quartiles (`:994-998`): it
prints the minimum in place of the 25%ile and the maximum in place of the 75%ile,
and `:1051-1064` fills the structured statistics field by field with the same
substitution. Three or more go through `SummaryStatistics`.

The branch writes **nothing at all** to the TSV report, and the `-m`, `-p` and
`-s` arms for FASTA are empty (`:2001-2004`, `:2112-2114`, `:2377-2379`), so
`-m` and `-s` contribute only their titles and `-p` its title and the
no-information line.

### consensusXML (`:1146-1311`)

`ConsensusMap::updateRanges` first, then the size histogram in **descending**
size. Each row carries the number of consensus features of that size, the
sub-features they account for, and both again restricted to the consensus
features that carry at least one peptide identification. The size column is
right-aligned in `largest_size / 10 + 1` characters (`:1221`) — the source's own
expression, which is not a digit count: a largest size of 100 gives a width of
11.

A second block (`:1193-1209`) re-counts at the level of peptides: the same
sequence and charge seen in *n* maps counts once, and contributes every
sub-feature it was seen in. Then two total lines, the second of which pads the
second column with `field_width` spaces. The histogram, the peptide rows, the
totals and the range block are all skipped when the map holds no consensus
feature; two lines of their own replace them. The column headers and the
assigned and unassigned identification counts are written either way.

`-m` prints the document identifier with **no** TSV twin, unlike the featureXML
arm. `-s` writes eleven blocks and, again unlike featureXML, none of them to the
TSV.

### idXML and mzIdentML (`:1312-1470`)

Three TSV lines — database, version, taxonomy — come first, from the first
protein identification run. Then the search engines, deduplicated and ordered as
the source's `set<pair<string, string>>` orders them; the run, protein-hit and
non-redundant protein-hit counts; the matched-spectrum, peptide-sequence,
PSM-per-spectrum, peptide-hit, modified-top-hit and non-redundant peptide-hit
counts; and, when there are any, one line of modification counts that does not
end in a newline.

Two numbers are deliberately coarse and are reproduced as they are:

- `PSMs / spectrum` (`:1416`) is `Size / int`, so it truncates: a file with
  seven hits over three spectra prints `2`. The structured
  `IdentInfo::psms_per_spectrum` carries the real ratio, which is what the
  source's own `Result` records at `:1464`;
- the average peptide length is `Math::round` of the mean of the hit lengths,
  streamed at the report's precision, while the modified-top-hit percentage
  (`:1418`) goes the other way and is built as a `std::string`, so
  `StringUtils::appendToStr(double)` renders it — `80.0%`, not `80%`.

Modification counting uses **two different identities** (`:1353-1372`): a
terminal modification is counted under `ResidueModification::getId()`, a residue
modification under `getFullId()`. So an oxidised methionine is
`Oxidation (M)` while an N-terminal dimethylation is `Dimethyl`, with no origin
suffix. `getId()` is empty for a user-defined mass-only modification
(`ResidueModification.cpp:593`, `:631`, `:671` set the full identifier and not
the identifier), and the port maps that case to the empty key too.

---

## 2. Source defects this port reproduces

### 2.1 The consensusXML `-s` quality sample is twice as long as it should be

`:2263-2266` declares

```cpp
vector<double> qualities(size);   // size zero-initialised values
qualities.reserve(size);          // a no-op
```

and then appends to it in the loop, so the sample ends up `2 * size` long with
`size` leading zeros. `intensities` is declared empty and merely reserved, so it
has the length one would expect. `widths` has the same defect and is never
printed.

This is deterministic, in bounds, and explained by the executed instructions, so
decision D1 asks for it to be reproduced rather than refused — and the upstream
reference output says the same thing: `FileInfo_7_output.txt` records five
consensus features with `Intensities ... num. of values: 5` and
`Qualities ... num. of values: 10`. Asserted in
`consensus_upstream_7_with_all_flags`.

### 2.2 `-m`, `-p` and `-s` report an empty experiment for mzIdentML

None of the three sections has an `MZIDENTML` arm, so an mzIdentML input falls
into the trailing `else //peaks` arm of each (`:2005`, `:2115`, `:2384`) and is
reported off the `MSExperiment` the identification branch never loaded. The
result is the peak-file metadata layout with every field empty, a date of
`0000-00-00 00:00:00`, the peak-file data-processing arm finding an empty
experiment, and a peak-file `Intensities:` block over no values. Reproduced
verbatim; `identifications_mzidentml_and_its_peak_file_fall_through`.

The port calls `peaks::write_meta` with a default `MSExperiment` rather than
writing the constant text, so the two renderings cannot drift apart.

### 2.3 The FASTA duplicate buckets keep only the last index per hash

`:931` and `:949` **assign** each bucket a one-element vector,
`m_headers[id_hash] = { index };`, instead of appending to it. The bucket is
therefore `last index with this hash`, and a duplicate is reported when that one
index matches. Three identical entries still count two duplicates, because #1 is
compared with #0 and #2 with #1.

The consequence is that the result depends on `std::hash<std::string>`: two
different strings that collide hide a duplicate a collision-free hash would have
reported. Reproducing the reference build therefore needs *its* hash, so
`fasta.rs::string_hash` implements libstdc++'s `_Hash_bytes` — the 64-bit Murmur
variant with multiplier `0xc6a4a7935bd1e995`, seed `0xc70f6907` and a 47-bit
shift-mix. `oracle/a7-fileinfo/scripts/probe_std_hash.cpp` pins it on the
reference toolchain over the empty string, every `length % 8` tail case, bytes
above `0x7f` and an embedded NUL; 28 of its 60 values are asserted in a unit
test.

---

## 3. What this port refuses (decision D1)

D1 refuses exactly where the source's behaviour is an out-of-bounds access, a
data race, process termination or a loop that never ends. Three sites qualify,
all confirmed against the Release build.

### 3.1 A consensus sub-feature whose map index is outside the column headers

`:1176-1183` sizes one occurrence vector from `getColumnHeaders().size()` and
then indexes it with `FeatureHandle::getMapIndex()`, which is the file's `map=`
**id** and not a position. A file whose map ids are not exactly `0..n-1` is out
of bounds.

Measured on the Release build:

| input | map ids | headers | result |
| --- | --- | --- | --- |
| `a7_cons_no_headers.consensusXML` | 0 | 0 | **SIGSEGV** |
| `a7_cons_mapindex_high.consensusXML` | 0, 5 | 2 | exit 0, silent corruption |
| `a7_cons_ids_one_based.consensusXML` | 1, 2 | 2 | exit 0, wrong peptide row |
| `ConsensusID_3_input.consensusXML` (**upstream**) | 1, 2 | 2 | exit 0, wrong peptide row |

The last one matters most: it is a fixture of the upstream test suite, and the
Release build prints

```
  peptides (with different mod. and charge) observed in 1 maps: 2	 (features: 2 )
```

for peptides that are in **both** maps — the write at index 2 lands outside the
two-slot vector, so the count only ever sees slots 0 and 1. The port returns
`Error::InvalidValue` before any of the report is written. The same file without
an identification on the offending consensus feature never reaches the indexing
and is reported normally, exactly as the Release build reports it
(`c_mapindex_high_noid`).

### 3.2 An identification file with no protein run

`:1336-1341` reads `id_data.proteins[0]` unconditionally, while the structured
block at `:1451` guards the very same access with
`if (!id_data.proteins.empty())`. `a7_id_no_runs.idXML` segmentation-faults on
the Release build. (This crate's idXML reader refuses such a file first, for its
own reason, so the refusal is reached either way.)

### 3.3 A peptide identification with no hit

`:1354` reads `getHits()[0]` behind the guard `!id_data.peptides[i].empty()`,
but `PeptideIdentification::empty()` (`PeptideIdentification.cpp:210-217`) tests
for a *default-constructed object*, not for an empty hit list: an identifier, a
score type, a non-zero significance threshold or `higher_score_better == false`
each make it false on their own.

The guard is in fact **unreachable for any loaded file**: `IdXMLFile::load`
gives every `PeptideIdentification` the enclosing `IdentificationRun`'s
identifier, so `id_` is never empty. `a7_id_empty_hitlist_ok.idXML` clears the
score type and still segmentation-faults, which is what shows this. Any idXML or
mzIdentML holding a `<PeptideIdentification>` with no `<PeptideHit>` therefore
crashes the reference FileInfo.

---

## 4. Native differences

1. **A non-ASCII FASTA residue byte is refused.** The source counts and prints
   raw `char`s, so such a byte reaches the report as that single byte. This
   port's reports are Rust strings and the FASTA reader only yields such a byte
   as part of a multi-byte UTF-8 sequence, so no faithful rendering exists;
   `fasta.rs::reject_non_ascii` refuses before anything is written. No upstream
   fixture contains one.
2. **The residue table is keyed by `u8`, ordered as signed `char`.** The source
   iterates a `std::map<char, int>`, and `char` is signed on the reference
   build's x86_64 Linux target, so a byte at or above `0x80` sorts *before* `A`.
   `fasta.rs::signed_char_order` renders in that order; the structured
   `FastaInfo::residue_counts` is a `BTreeMap<u8, _>` and orders the other way.
   Difference 1 refuses every input that could tell the two apart.
3. **The duplicate warnings go to `FileInfoResult::warnings`.** The source
   writes them with `OPENMS_LOG_WARN`, not into either report; the library
   prints nothing and the caller decides. `fasta_duplicate_warnings_match_the_source_log`
   compares them with the tool's stderr, byte for byte and in order.
4. **A forced type that does not match the content is refused with this crate's
   message.** `FileHandler::loadIdentifications` refuses an idXML forced to
   mzIdentML with `type: idXML is not allowed for loading identifications` and
   the tool exits 3; the port returns
   `Error::InvalidValue("idXML is not an allowed input format")`, which is how
   `FileHandler` maps that refusal throughout this crate.

---

## 5. Known gap outside this package

The shared identification-XML reader cannot load an idXML with more than **14
modified peptide hits**, whatever the file size.
`AASequence::parse_with_budget` (`src/chemistry/sequence.rs:1147-1180`) charges
a per-modified-sequence preflight proportional to `ModificationsDB::len()`
against the single document-wide `XmlLimits::max_work = 50_000_000`
(`src/format/identification_xml.rs:411-425`). Measured: 14 modified hits load,
15 do not, while 800 *unmodified* hits are fine.

Two oracle cases are blocked by it — `FileFilter_25_input.idXML` (473 modified
hits) and `FalseDiscoveryRate_5_input.idXML` (75) — so they have no differential
here; the oracle records what the C++ prints for both. This is a limit of that
reader and of the sequence parsing budget, not of these branches, and it is
raised for the lead rather than worked around here.
