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
Oracle: `../oracle/a7-fileinfo`, **75** cases against the Release C++ FileInfo of
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
ibminode06, run twice and reproduced. **55** of them have both their `-out` and
`-out_tsv` reports compared byte for byte in the test, and none carries an
exemption of any kind. The count was 46 before the NaN-spelling step of
2026-09-20 closed native difference 5; the nine reports that were compared line
for line under the `nan` / `-nan` exception are now compared like the rest, and
the two of them that section 5.2 once called "retained as evidence rather than
compared" are among them.

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
5. ~~**Every NaN the FileInfo text layer prints is spelled `nan`, where glibc
   spells a NaN whose sign bit is set `-nan`.**~~ **CLOSED** on 2026-09-20. It
   is kept here with its measurement, because the measurement is what makes the
   closure checkable. It was a class of line, not one line: any of the seven
   value lines a `SummaryStatistics` block prints — the mean, the five order
   statistics and the variance — can carry one, and the consensusXML `-s`
   blocks are the first FileInfo path whose own arithmetic can produce one at
   all.

   *Scope of the claim, as measured.* When native difference 6 closed, this was
   the **only** class of line on which any compared oracle report disagreed with
   the Release build. Nine reports carry a sign-bit NaN, and how many lines of it
   each holds depends on the input: `c_zero_intensity_s` and
   `c_zero_intensity_all` one each, `c_nan_then_finite_s` and
   `c_finite_then_nan_s` seven each, `c_nan_one_s`, `c_nan_one_all`,
   `c_zero_swapped_s` and `c_zero_swapped_all` nine each, and `c_nan_two_s` ten.
   (The earlier text of this paragraph said "nine reports" and then named eight,
   dropping `c_zero_intensity_all`; the provenance manifest had it right.) All
   nine are now compared byte for byte, and `assert_report_with_signed_nans`
   keeps the counts above as a tripwire: byte equality alone would still hold if
   both the port and a regenerated expectation stopped writing the sign, so the
   assertion additionally reads the counts out of the Release build's own frozen
   report, and asserts that neither report spells a NaN without its sign
   anywhere.

   *What closed it.* Two steps, in this order, because the spelling is only
   honest once the value is.

   - *The value.* `src/format/file_info/consensus.rs`'s `it_aad += it_ratio` —
     the `(-inf) + (+inf)` of section 5.2 — was plain Rust arithmetic, so the
     NaN it generated carried the host's sign bit, positive on arm64. It now
     goes through `crate::math::x86_64::add`, which answers an invalid
     operation with SSE2's default NaN `0xfff8000000000000`. A sweep of
     `src/format/file_info/` found no second instance: the consensusXML
     reader's `map.validate()` refuses a non-finite rt, m/z, intensity or
     width, so `rt_diff`, `mz_diff` and `it_ratio` are all operations on finite
     values, the `rt_aad`/`mz_aad` sums accumulate values already made
     non-negative and cannot cancel to a NaN, and the featureXML branch's
     `tic += intensity` runs after `RangeBase::extend_value` has already
     refused a non-finite intensity.
   - *The spelling.* `text_format::nonfinite` writes `-nan` for a NaN whose
     sign bit is set. All three of its call sites reach C `printf`;
     `StringUtils::toStr` does not come through it and still writes `NaN` for
     either sign, which is `NumericFormatting.h:29`.

   *What the spelling is measured on, and what is generalised.*
   `../oracle/a2-textfmt-linux` re-ran A2's `driver.cpp`, `cases.h` and
   `pin_probe.cpp` — byte-identical to the macOS oracle's — against the Linux
   x86_64 Release install on ibminode06. Of 1018 rows exactly one differs from
   the macOS capture, `fff8000000000000`, and it differs in the five `printf`
   columns and not in the `toStr` column. That corpus holds exactly three NaN
   bit patterns and exactly one sign-bit NaN row, so what is **measured** is
   `number(-NaN, n)` at `n` in `{0, 1, 2}` and `ostream(-NaN, p)` at `p` in
   `{6, 15}`; no sign-bit NaN `float` is pinned at all. Every other digit count,
   precision and the `float` overload are **generalised** from glibc writing the
   sign before `__printf_fp` dispatches on the class. A wider negative-NaN
   sweep, `../oracle/a2-textfmt-nan-sweep`, was captured by another lane while
   this one ran; it is named rather than cited, because it is not registered in
   this repository's manifests and is the lead's to fold in.

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

   *Why the text layer used to write `nan` anyway, and what changed.*
   `text_format`'s `nonfinite` ignored the sign of a NaN, on two grounds. The
   first was that the sign of a *generated* NaN belonged to the hardware —
   AArch64's default NaN is the positive one, and this crate's `cross-platform`
   CI job runs the full suite on `macos-latest`, on tags and on
   `workflow_dispatch` (`.github/workflows/rust.yml:127`, `:134`). The
   shared-math wave of 2026-09-19 and the value step above removed it: every
   NaN these branches can generate now carries the Release build's bits on any
   host. The second was that the pinned oracle row said `nan` — but it said so
   because it was captured with **Apple libc**, on the macOS product SDK, which
   is not the reference platform. `../oracle/a2-textfmt-linux` re-ran the same
   driver against the Linux Release install and row 92 reads
   `D fff8000000000000 -nan -nan -nan -nan -nan NaN`. That is the row
   `tests/file_info_text_format.rs` now embeds, and with both grounds gone
   `nonfinite` carries the sign.

   *How it is pinned.* Every Release report is frozen whole and every one of
   them — bare, `-s` and `-all` alike — is compared byte for byte through
   `check` or `assert_report_with_signed_nans`; nothing is exempted.
   `assert_report_with_signed_nans` is `assert_report` plus two tripwires that
   byte equality alone would not give: the number of lines of this class the
   *reference* report carries, and that neither report spells a NaN without its
   sign anywhere. So a regeneration that quietly dropped the sign on both sides,
   or a change that removed the NaN-bearing lines, still fails.
   `consensus_zero_intensity_sub_feature_makes_the_variance_a_nan` pins one such
   line in each of two reports, and `consensus_nan_in_the_statistics_sample`
   pins nine, nine, ten, seven and seven;
   `consensus_a_signed_zero_sample_keeps_the_release_builds_order` pins nine
   twice more.
6. ~~**A sample holding both a negative and a positive zero is ordered by the
   IEEE-754 total order, where `std::sort` leaves it as it found it.**~~
   **CLOSED** in the shared-math wave of 2026-09-19, under lead decision D16.
   It is kept here with its measurement, because the measurement is what makes
   the closure checkable.

   *Where it came from.* `:2310-2311` pushes every intensity ratio into
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
   sub-features appear in the file. This port's `sort_ascending` used
   `f64::total_cmp`, which orders `-0.0` before `0.0` regardless of input order,
   so it printed one of the two answers for both.

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

   *How it was closed.* `sort_ascending` is now
   `source_sort_by(&mut values, |a, b| a < b)` — `std::sort(begin, end)` itself,
   reproduced comparison by comparison by `crate::math::source_sort` from the
   GCC 14.4.0 libstdc++ headers the reference build was compiled with. A
   two-element range of equivalent values comes back in the order it went in, so
   the port prints the **left** column for `a7_cons_nan_one` and the **right**
   column for `a7_cons_zero_swapped`: both members of the measured pair, rather
   than one answer for both.

   *How the closure is pinned.*
   `consensus_a_signed_zero_sample_keeps_the_release_builds_order` (renamed from
   `consensus_a_signed_zero_sample_is_ordered_by_the_total_order`) asserts the
   two bare reports byte for byte, that the two frozen Release `-s` reports
   disagree on exactly those four lines plus the `Ranges` line, and then — the
   part that changed — that **both** the port's swapped report and its unswapped
   report match their own references through
   `assert_report_with_signed_nans`. Since native difference 5 closed that is
   byte-for-byte equality, with the nine sign-bit NaN lines of each report
   asserted as a tripwire rather than exempted. Zero order-statistic lines
   differ.
   `a_signed_zero_keeps_the_order_the_release_build_keeps` in
   `tests/statistic_functions.rs` pins the same behaviour at the
   `SummaryStatistics` level, in both input orders.

   *What is still reproduced the way it always was.* The `Ranges` line
   `intensity: -0.00 .. 100.00` against `intensity: 0.00 .. 100.00`, which the
   same equivalence produces through `std::min` in `updateRanges`; the bare and
   `-out_tsv` reports of both files, byte for byte.
   `FileInfo.cpp:2257-2372` writes nothing to `os_tsv`, so only the text report
   carries the statistics at all.

---

## 5. Known gaps outside this package, and one that closed

Section 5.1 is a gap. Section 5.2 is kept in place with its measurement because
the measurement is what makes its closure checkable, as native difference 6 is
kept in section 4.

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

### 5.2 A NaN next to a number in a `SummaryStatistics` sample — **CLOSED**

Closed in two steps: decision **D16** (shared-math wave, 2026-09-19) made the
shape *reproducible*, and the NaN-spelling step of 2026-09-20 made its text
byte-identical to the reference build's. The provenance manifest marks it
`RESOLVED`. Everything below is the record of how it was settled, not an open
item.

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

**That question is no longer the one that decides.** Lead decision **D16**
(shared-math wave, 2026-09-19) says that reproducing a permutation the standard
leaves *unspecified* is in scope, because the Release build runs one particular
algorithm deterministically and the port already reproduces it comparison by
comparison (`crate::math::source_sort`, tier 1 over 2,272 inputs). Since that
decision `sort_ascending` is `source_sort_by(&mut values, |a, b| a < b)` —
`std::sort(begin, end)` itself — and **every shape in this section is
reproduced**, not only the two whose output set is a singleton:

| sample | oracle case | before D16 | after |
| --- | --- | --- | --- |
| one value | `c_nan_one_s` | reproduced (one permutation exists) | reproduced |
| every value a NaN | `c_nan_two_s` | reproduced (all permutations print the same) | reproduced |
| a NaN then a number | `c_nan_then_finite_s` | **refused** | reproduced |
| a number then a NaN | `c_finite_then_nan_s` | **refused** | reproduced |

For the first the reference prints the NaN on the mean and on all five order
statistics, and `0` for the variance — the `n <= 1` substitution of
`StatisticFunctions.h:951`. For the second the variance is a NaN too, because
`n > 1` lets `Math::variance` run, so all seven value lines carry one and only
`num. of values` does not. `SummaryStatistics::of_nan_sample`, which used to
compute those two without sorting, is **gone**: they are ordinary sorts now.

**What made the last two hard, and what settles them.** The `minimum`, quartile
and `maximum` lines the reference prints are positional reads of a range whose
elements `std::sort` was free to leave in any order, and the two files prove it:
`c_nan_then_finite_s` and `c_finite_then_nan_s` hold **the same two consensus
features in opposite file order**, are each stable over three runs, and disagree
on exactly four lines:

```text
                    c_nan_then_finite_s      c_finite_then_nan_s
  minimum:          -nan                     2
  lower quartile:   -nan                     2
  upper quartile:   2                        -nan
  maximum:          2                        -nan
```

A genuinely sorted range cannot have that property, which is why `CPP-347`
stands as a C++ defect. But "which permutation" is a question with a measured
answer for this build: for a two-element range libstdc++ runs one
`__insertion_sort` pass, `2 < NaN` and `NaN < 2` are both false, so nothing
moves and the sample keeps its file order. The port reproduces that, and both
reports are now compared line for line by
`consensus_nan_in_the_statistics_sample` — seven differing lines each, all seven
of the `nan` / `-nan` class of native difference 5, nothing else.

"Moves nothing" is a property of the sample's **size and arrangement**, not of
the NaN, and the port reproduces that too. `__introsort_loop` runs only above
`_S_threshold`, which the headers this build was compiled with enumerate as 16,
so at 17 elements or more it is free to move the NaN and for the `{NaN, 2..n}`
family it does — `{NaN, 2..20}` prints `minimum: 2` and `median: -nan` — and
even below the threshold a block move can carry it, as `{3, NaN, 2}` sorting to
`{2, 3, NaN}` shows. `the_introsort_threshold_decides_where_a_nan_lands` in
`tests/statistic_functions.rs` pins the whole permutation for the 16-, 17- and
20-element cases; the full reading of the headers is in
[the shared-math document](STATISTIC_FUNCTIONS_SUPPORT.md#nan-policy).

The manifest of `../oracle/a7-fileinfo` records both reports under
`unspecified_order` with the reason; that annotation is now a statement about
the C++ standard's guarantee rather than about what the port does with them.

**The same boundary without a NaN is closed with it.** `-0.0` and `0.0` are also
equivalent under `operator<`, so a sample holding both is left in file order by
`std::sort` for exactly the same reason — and since D16 the port leaves it there
too. That was native difference 6 of section 4, measured on `c_nan_one_s`
against `c_zero_swapped_s`; it is now pinned as an *equality* by
`consensus_a_signed_zero_sample_keeps_the_release_builds_order`. Until the A7
round it was the one instance of this boundary that was *silent* — recorded in
prose in `docs/FILE_INFO_SUPPORT.md` item 6, but with no oracle case, no frozen
report and no test, so neither a regression nor a fix would have been noticed.

**What was left, and is now done.** At the end of the shared-math wave the five
frozen reports were still not byte-identical to the reference, for two reasons
outside that wave's scope. Both are closed:

- `src/format/file_info/consensus.rs` generated its NaN in plain Rust
  arithmetic — `it_aad += it_ratio` is the `(-inf) + (+inf)` of this very
  section — so that value was host-shaped even though every value
  `statistic_functions` produces was not. It goes through
  `crate::math::x86_64::add` now, and a sweep of the module found no second
  instance.
- A2's oracle row had to be re-captured against the Linux Release build rather
  than the macOS SDK. `../oracle/a2-textfmt-linux` did that, and
  `text_format`'s `nonfinite` now spells a sign-bit NaN `-nan` as glibc does.

All nine reports that carried native difference 5 are compared byte for byte,
with no exemption: the four `-s` reports of the shapes tabulated above
(`c_nan_one_s`, `c_nan_two_s`, `c_nan_then_finite_s`, `c_finite_then_nan_s`),
the signed-zero pair `c_zero_swapped_s` and `c_zero_intensity_s`, and the three
`-all` reports `c_nan_one_all`, `c_zero_swapped_all` and
`c_zero_intensity_all`. The `assert_report_but_the_nan_spelling` helper that
existed only for this difference is gone.

**The same closure for `FileInfo -c`, decision D18.** `-c` refused a NaN MS1
retention time or peak m/z on the same grounds this section once used: the
source's `std::sort` leaves the order undefined. `FileInfo.cpp:1927` and `:1956`
are the same unqualified `sort(v.begin(), v.end())` on a `std::vector<double>`
that D16 now reproduces — the vectors are declared at `:1863` and `:1942`,
neither call carries a comparator, and `:47` is `using namespace std;` — so the
refusal is closed with the same machinery. It is observable: for MS1 retention
times `{5.0, NaN, 5.0}` libstdc++ moves nothing and the Release build prints no
duplicate line, where a `f64::total_cmp` sort prints one. Item 3 of
*Native differences* in `docs/FILE_INFO_CHECKS_SUPPORT.md` carries the full
record.
