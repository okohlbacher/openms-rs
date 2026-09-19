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
Oracle: `../oracle/a7-fileinfo`, 72 cases against the Release C++ FileInfo of
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
ibminode06, run twice and reproduced. 50 of them have both reports compared
byte for byte in the test; two more are retained as evidence for section 5.2
rather than compared.

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
the Release build.

The port refuses the file, but **not from this branch**: the shared reader
refuses it first, and `identifications.rs`'s `data.proteins.first()` is defence
in depth that no accepted input reaches. Measured on `a7_id_no_runs.idXML`:
`Error::Parse`, `idXML needs at least one IdentificationRun`, and the tool exits
**3**, not the 6 the branch's own `Error::InvalidValue` would give. The other
reader answers the same shape: `mzidentml.rs:1596-1598` refuses a document with
no `SpectrumIdentification` element before it can produce a run-less
`MzIdentMLDocument`, and `identifications::report` has no entry that does not
load from a file. So the guard has no executed coverage on any input, and none
can be written while either reader keeps its own refusal; it stays because the
branch must not index an empty vector if a future reader ever hands it one.
`identifications_without_a_run_are_refused` pins the measured refusal, with its
message.

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
5. **Every NaN the FileInfo text layer prints is spelled `nan`, where glibc
   spells a NaN whose sign bit is set `-nan`.** This is a class of line, not one
   line: any of the seven value lines a `SummaryStatistics` block prints — the
   mean, the five order statistics and the variance — can carry one, and the
   consensusXML `-s` blocks are the first FileInfo path whose own arithmetic can
   produce one at all.

   *Scope of the claim, as measured.* It and native difference 6 are the only
   two classes of line on which the 53 compared oracle reports disagree, and
   this one is the only class that reaches more than one report. How many lines of it a report holds
   depends on the input: `c_zero_intensity_s` has one, `c_nan_one_s` nine and
   `c_nan_two_s` ten. The tests pin those counts and check each line against the
   class rather than against a single expected line.

   *Where it comes from.* `:2310` computes
   `it_ratio = element_intensity / (centroid_intensity > 0 ? centroid_intensity : 1)`
   and `:2312-2315` replaces every ratio below 1 by its reciprocal, so a
   sub-feature of intensity 0 under a centroid of positive intensity
   contributes `1 / 0 = +inf`, and one of intensity `-0.0` contributes
   `1 / -0.0 = -inf`. `Math::SummaryStatistics`
   (`StatisticFunctions.h:933-958`) then summarises `{1, +inf}`: the mean is
   `+inf`, and `Math::variance` (`:541-556`) adds `(1 - inf)^2 = +inf` to
   `(inf - inf)^2 = NaN`, so the `Relative intensity error` block reports a NaN
   variance. A `-inf` and a `+inf` in the same block make the `mean` and the
   `median` NaNs too, and `:2317` adding them makes a NaN that goes into the
   *sample* of the per-consensus-feature block rather than into a statistic
   summarised out of one — see section 5.2. Nothing is out of bounds, and both
   oracle runs agree, so D1 asks for all of it to be reproduced rather than
   refused. Oracle cases `c_zero_intensity{,_s,_all}` on
   `a7_cons_zero_intensity.consensusXML`, `c_nan_one{,_s,_all}` on
   `a7_cons_nan_one.consensusXML` and `c_nan_two{,_s}` on
   `a7_cons_nan_two.consensusXML`; all of them exit 0.

   *`variance` is not the only statistic these branches can make non-finite.*
   The frozen `c_zero_intensity_s` report itself carries `mean: inf` and
   `maximum: inf`, and the `c_nan_two_s` report carries a NaN on every value
   line of its block — the mean, all five order statistics and the variance. What is true of the oracle's cases is narrower and is what the tests
   assert: the `nan`/`-nan` class is the only way any of them disagrees.

   *The value is reproduced; only the text differs.* Measured on x86_64:
   `inf - inf` is `0xfff8000000000000`, SSE2's default NaN, whose sign bit is
   set — in the reference build (probe compiled on ibminode06, fed from `argv`
   so nothing is constant-folded, `printf` and `std::ostream` both `-nan`) and
   in this crate's `variance_with_mean` alike (the same probe in Rust on kim,
   optimised and unoptimised, `sign_negative=true`). The port computes the
   identical bits.

   *Why the text layer writes `nan` anyway.* `text_format`'s `nonfinite`
   ignores the sign of a NaN by design, and that design is not this package's:
   the sign of a *generated* NaN belongs to the hardware — AArch64's default
   NaN is the positive one, and this crate's `cross-platform` CI job runs the
   full suite on `macos-latest`, on tags and on `workflow_dispatch`
   (`.github/workflows/rust.yml:127`, `:134`) — and Apple libc writes `nan` for
   `0xfff8000000000000` regardless. `../oracle/file-info-text-format` measured
   exactly that bit pattern (`results/driver.tsv`, the `D fff8000000000000`
   row) and `tests/file_info_text_format.rs` asserts the `nan` it produced.
   Spelling the sign here would contradict that executed row and make every
   frozen expectation architecture-dependent. The A2 module note at
   `src/format/file_info/text_format.rs` states the rule.

   *Making the sign printable is its own wave, not this package's.* Spelling it
   honestly would first have to make the value host-independent: an
   x86_64-faithful `variance_with_mean` in `src/math/statistic_functions.rs`,
   which needs the x86_64 emulation promoted out of
   `analysis::feature_finder_picked::scoring` into shared math, plus A2's oracle
   row re-captured against the Linux Release build instead of the macOS SDK.
   That is a cross-cutting change to shared math which landed ports already
   consume, so it is carried forward for the lead rather than done here.

   *How it is pinned.* Every Release report is frozen whole. The bare ones match
   byte for byte through `check`; for the `-s` ones,
   `assert_report_but_the_nan_spelling` takes the number of lines of this class
   the report has and, for each differing line, asserts that the reference ends
   in `-nan`, that ours ends in `nan` and not `-nan`, and that putting the sign
   back reproduces the reference line character for character. A divergence
   outside the class, a change in how many lines carry a NaN, a `-nan` where the
   reference has a number, or a change of either spelling fails the test.
   `consensus_zero_intensity_sub_feature_makes_the_variance_a_nan` pins one such
   line and `consensus_nan_in_the_statistics_sample` pins nine, nine and ten.
6. **A sample holding both a negative and a positive zero is ordered by the
   IEEE-754 total order, where `std::sort` leaves it as it found it.** This is
   the same `std::sort` boundary as section 5.2 with **no NaN anywhere**, and
   unlike section 5.2 the port *accepts* the input: it exits 0, as the Release
   build does, and reproduces every line of the report except four.

   *Where it comes from.* `:2310-2311` pushes every intensity ratio into
   `it_delta_by_elems` **before** `:2312-2315` inverts the ones below 1, so a
   consensus feature with a sub-feature of intensity `-0.0` and one of `0.0`
   under a positive centroid contributes `-0.0` and `0.0` to the
   `Intensity ratios` sample. `Math::SummaryStatistics` hands that sample to
   `std::sort` (`StatisticFunctions.h:948`) with the default `operator<`, under
   which `-0.0 < 0.0` and `0.0 < -0.0` are **both false**. The two are therefore
   *equivalent*, the strict-weak-ordering precondition **holds** — nothing here
   is undefined — and every permutation is a conforming result. libstdc++
   leaves a range this size as it found it, `front()`, the quantiles and
   `back()` at `:952-956` are positional reads, and `ostream` writes `-0` for a
   negative zero, so the four lines are a property of the order the
   sub-features appear in the file. This port's `sort_ascending` uses
   `f64::total_cmp`, which orders `-0.0` before `0.0` deterministically, so it
   prints one of the two answers for both orders.

   *Measured, both orders.* `a7_cons_nan_one.consensusXML` and
   `a7_cons_zero_swapped.consensusXML` hold the **same consensus feature with
   its two sub-feature intensities exchanged**. On the Release build, each run
   twice and reproduced (oracle cases `c_nan_one_s` and `c_zero_swapped_s`,
   annotated `signed_zero_order` in the manifest):

   | line | `a7_cons_nan_one` | `a7_cons_zero_swapped` |
   | --- | --- | --- |
   | `minimum:` | `-0` | `0` |
   | `lower quartile:` | `-0` | `0` |
   | `upper quartile:` | `0` | `-0` |
   | `maximum:` | `0` | `-0` |

   The port prints the left-hand column for **both** files, so it agrees with
   the reference on `c_nan_one_s` and differs from it on `c_zero_swapped_s`
   on exactly those four lines. Each of the four lines the port prints is the
   Release build's own line for the other file, which the test asserts: nothing
   here is invented, the port simply always picks the permutation `total_cmp`
   puts first.

   *What is reproduced.* Everything else, including the `Ranges` line
   `intensity: -0.00 .. 100.00` against `intensity: 0.00 .. 100.00`, which the
   same equivalence produces through `std::min` in `updateRanges` and which the
   port **does** follow, because that path keeps the file's own order on both
   sides. The bare and `-out_tsv` reports of both files are reproduced byte for
   byte; `FileInfo.cpp:2257-2372` writes nothing to `os_tsv`, so only the text
   report carries the statistics at all.

   *Why it is a difference and not a refusal.* Refusing the shape would widen a
   refusal in `sort_ascending`, which every `SummaryStatistics` caller in the
   crate consumes, to an input the Release build handles in bounds, stably and
   with a well-defined precondition — a decision the lead has not taken (see
   section 5.2's open question, which this instance widens). Reproducing it
   means porting libstdc++'s `std::sort` permutation, which is the same
   shared-math wave `CPP-347` names; that wave has to cover **any** sample
   whose elements `std::sort` calls equivalent but `f64::total_cmp` orders, not
   only NaN-bearing ones.

   *How it is pinned.*
   `consensus_a_signed_zero_sample_is_ordered_by_the_total_order` asserts the
   two bare reports byte for byte, that the two frozen Release `-s` reports
   disagree on exactly those four lines plus the `Ranges` line, that the port's
   swapped report differs from its reference on exactly thirteen lines — nine
   of the class of native difference 5 and those four — and that each of the
   four equals the Release build's own line for the unswapped file. Before this
   test the divergence was real but unpinned: `docs/FILE_INFO_SUPPORT.md` item
   6 recorded that "a sample holding both zeros can print `-0` where the source
   prints `0`", and nothing in the suite would have noticed either a regression
   or a fix.

---

## 5. Known gaps outside this package

### 5.1 The identification-XML reader's modified-hit budget

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

### 5.2 A NaN next to a number in a `SummaryStatistics` sample

`:2310` divides and `:2312-2315` inverts, so a consensus feature with one
sub-feature of intensity `-0.0` and one of `0.0` under a positive centroid
contributes `(-inf) + (+inf)` at `:2317`, which `:2323` divides by `cm.size()`
and `:2327` pushes into `it_aad_by_cfs` — the *sample* of the `Average relative intensity error
within consensus features` block, not a statistic summarised out of one. The
sample then goes to `std::sort`.

Under `operator<` a NaN is incomparable with every value, itself included, so
`std::sort` is free to return **any** permutation of the elements it cannot
tell apart — and if the sample also holds two or more distinct numbers,
transitivity of incomparability fails (`1 ~ NaN` and `NaN ~ 3` while `1 < 3`)
and the call is undefined outright. Either way the question is the same: *is
the set of outputs the source may produce a singleton?*

Two shapes answer yes, and both are **reproduced**:

| sample | why the set of outputs has one member | oracle case |
| --- | --- | --- |
| one value | a one-element range has exactly one permutation | `c_nan_one_s` |
| every value a NaN | `std::sort` may permute freely, but every permutation prints the same eight lines | `c_nan_two_s` |

For the first the reference prints the NaN on the mean and on all five order
statistics, and `0` for the variance — the `n <= 1` substitution of
`StatisticFunctions.h:951`. For the second the variance is a NaN too, because
`n > 1` lets `Math::variance` run, so all seven value lines carry one and only
`num. of values` does not.
`SummaryStatistics::of_nan_sample` (`src/math/statistic_functions.rs`) computes
both without sorting, since there is nothing to order.

**A NaN next to a number is refused**, with
`Error::InvalidValue("statistics input must not contain NaN")`, where the
reference build exits 0. The reason is measured rather than assumed: in the
two-element samples measured here libstdc++ compares every pair involving the
NaN false and therefore moves nothing, so the `minimum`, quartile and `maximum`
lines the reference prints are positional reads of a range whose elements
`std::sort` was free to leave in any order. "Moves nothing" is a property of
this sample's **size and arrangement**, not of the NaN. `__introsort_loop` runs
only
above `_S_threshold`, which the headers this build was compiled with enumerate
as 16, so at 17 elements or more it is free to move the NaN, and for the
`{NaN, 2..n}` family measured here it does — a 20-element
sample `{NaN, 2..20}` prints `minimum: 2` and `median: -nan` — and even below
the threshold a block move can carry it, as `{3, NaN, 2}` sorting to
`{2, 3, NaN}` shows. Neither makes the order statistics any less a property of
the input order, or the refusal any less necessary; the full reading of the
headers is in
[the shared-math document](STATISTIC_FUNCTIONS_SUPPORT.md#nan-policy).
Oracle cases `c_nan_then_finite_s` and `c_finite_then_nan_s` hold **the same
two consensus features in opposite file order**, are each stable over three
runs, and disagree on exactly four lines:

```text
                    c_nan_then_finite_s      c_finite_then_nan_s
  minimum:          -nan                     2
  lower quartile:   -nan                     2
  upper quartile:   2                        -nan
  maximum:          2                        -nan
```

Reproducing those values means porting libstdc++'s `std::sort` permutation into
`sort_ascending`, which every `SummaryStatistics` caller in the crate consumes.
That is shared-math work of its own wave, so **this is a deferral, not a D1
refusal**, and it is raised for the lead. The manifest of
`../oracle/a7-fileinfo` records both reports under `unspecified_order` with the
reason, and `consensus_nan_in_the_statistics_sample` asserts the exact refusal,
that the run without `-s` still succeeds, and that the two retained reference
reports disagree on exactly those four lines.

**The same boundary without a NaN, which this port does not refuse.** `-0.0`
and `0.0` are also equivalent under `operator<`, so a sample holding both is
left in file order by `std::sort` for exactly the same reason — but there the
precondition *holds*, nothing is undefined, and this port accepts the input and
prints the order `f64::total_cmp` gives. That is native difference 6 of section
4, measured on `c_nan_one_s` against `c_zero_swapped_s` and pinned by
`consensus_a_signed_zero_sample_is_ordered_by_the_total_order`. It matters for
the lead's decision in two ways: the shared-math wave has to cover **any**
sample whose elements `std::sort` calls equivalent but `f64::total_cmp` orders,
not only NaN-bearing ones; and until this round it was the one instance of this
boundary that was *silent* — recorded in prose in `docs/FILE_INFO_SUPPORT.md`
item 6, but with no oracle case, no frozen report and no test, so neither a
regression nor a fix would have been noticed.
