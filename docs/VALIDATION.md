# Validation of the ongoing Rust port

## Shared-math wave: the source's own arithmetic and its own `std::sort` (2026-09-19)

Decision **D16** was taken for this wave and is recorded in full, with its
provenance, in [the decision list below](#lead-decisions-d1-d17): **reproducing an
unspecified `std::sort` permutation is in scope**, because the port already does
it (`src/math/source_sort.rs`, tier 1 against two oracle drivers over 2,272
inputs), because doing so is D1-compliant by construction (it reproduces the
in-bounds measured behaviour and refuses exactly where the introsort reads
outside the vector), because refusing is worse (the signed-zero half is ordinary
finite data), and because the pin risk is managed (`source_sort` names the
sha256 of every libstdc++ header it reproduces). That decision closes the two
bullets wave 8 left open below.

The wave has three parts, all in the same two files plus the promotion.

**Part 1 — the promotion, behaviour-preserving.** `scoring::x86_64`,
`scoring::libstdcxx` and `analysis::feature_finder_picked::source_sort` are now
`src/math/x86_64.rs`, `src/math/libstdcxx.rs` and `src/math/source_sort.rs`.
They are not the picked feature finder's: they are the arithmetic and the
ordering the Release build's own results depend on, and shared math needs them.
Visibility is unchanged (the two emulation modules stay crate-private,
`source_sort` stays public), 28 path references across 14 Rust files follow, and
`RUSTDOCFLAGS=-D warnings cargo doc` is clean. That nothing moved but paths was
*checked*, not asserted: the extracted `x86_64` body is token-identical to the
block it came from, and `libstdcxx` differs only where rustfmt dropped a
trailing comma from a signature that fits on one line after the dedent.
`python3 tools/check_module_cycles.py` confirms the direction is cycle-safe —
`analysis → math` exists, `math → analysis` does not — at 64 cross-module edges
and 13 mutually-dependent pairs, unchanged. `glibc_libm.rs` was deliberately
left where it is and is the obvious candidate for the next promotion.

**Part 2 — the NaN bits.** IEEE 754 fixes every finite result of `+ - * / sqrt`
but not which NaN bit pattern an operation *generates* out of non-NaN operands.
The Release build's SSE2 answers `0xfff8000000000000`, which glibc spells
`-nan`; an arm64 host answers `0x7ff8000000000000`. Native difference 5 is that
spelling, and the *value* had to stop depending on the host before the spelling
could be argued about at all. Every function in
`src/math/statistic_functions.rs` whose own arithmetic can generate a NaN is now
built on `crate::math::x86_64`: `sum`, `mean`, `variance`, `variance_with_mean`,
`sd`, `sd_with_mean`, `covariance`, `mean_square_error`,
`root_mean_square_error`, `mean_absolute_deviation`, `mad`'s `fabs`, and the two
places an order statistic is *interpolated* rather than read —
`median_of_sorted`'s even-size average and `quantile`'s linear blend. The module
documentation states per function why it is or is not affected.

The invariant that makes this safe is that the helpers return the IEEE result
whenever it is not a NaN, so no finite value may move, and that is asserted
directly rather than hoped for: `the_x86_64_helpers_change_no_finite_result`
compares bit for bit against the plain-Rust arithmetic the module used before,
over a battery reaching subnormals, both zeros, `DBL_MAX` and ranges whose
squared deviations overflow to an infinity, across every function and every
equally long pair. **No frozen expectation in `tests/statistic_functions.rs`
moved.** `a_generated_nan_carries_the_release_builds_bits` then pins fifteen
bit patterns derived from the SSE2 rules of Intel SDM vol. 1 rather than from
this crate's output, including the *positive* NaN `andpd` leaves behind in
`mean_absolute_deviation` and the quieted payload a signalling NaN input keeps.

One correction to the brief this wave was given, recorded because it is the kind
of thing a reader would otherwise take on faith:
`variance_with_mean(&[+inf, -inf], 0.0)` does **not** produce a NaN. Both
squared deviations are `+inf` and they add, so the result is `+inf`; the test
pins that too. The shapes that do generate the indefinite NaN are the ones where
a *subtraction* cancels two infinities — `variance(&[1.0, +inf])`, which is the
`FileInfo` path itself, and `variance_with_mean(&[+inf, -inf], +inf)`.

**Part 3 — the permutation.** `sort_ascending` was
`values.sort_by(f64::total_cmp)`, and every entry point above it refused a NaN.
It is now `source_sort_by(&mut values, |a, b| a < b)`: `std::sort(begin, end)`
with the default `operator<` and no comparator, which is the call at
`MATH/StatisticFunctions.h:140` (`median`), `:244` (`quantile1st`), `:281`
(`quantile3rd`), `:189` (`MAD`, through its own `median`) and `:948`
(`SummaryStatistics`'s unqualified `sort(data.begin(), data.end())`), core
`bc9cc12`. `median`, `quantile1st`, `quantile3rd`, `mad` and
`SummaryStatistics::new` therefore **reproduce** a NaN-bearing sample instead of
refusing it, and the only refusal left in the sort path is D1's: an introsort
read outside the vector, which an asymmetric comparison cannot provoke and `<`
on `f64` keys — NaN keys included — is asymmetric.
`SummaryStatistics::of_nan_sample` and the two-shape special case it existed for
are gone; `new` reads its order statistics through the private `_of_sorted`
helpers, because `std::sort`'s own output is not ascending when a NaN is in it.

The evidence is the permutation itself, not the statistic.
`the_introsort_threshold_decides_where_a_nan_lands` pins the whole array for
`{NaN, 2..16}`, `{NaN, 2..17}` and `{NaN, 2..20}`, whose NaN lands at index 0, 8
and 10 — the three positions the "Where a NaN lands" section of
[STATISTIC_FUNCTIONS_SUPPORT](STATISTIC_FUNCTIONS_SUPPORT.md) measured against
the reference compiler, and the reason `{NaN, 2..20}` prints `minimum: 2` and
`median: -nan`. The port reproduces all three, and the 16/17 boundary is
`_S_threshold` (`bits/stl_algo.h:1806`, `:1880`, `:1899-1910`) doing exactly
what the doc says it does. `{3, NaN, 2}` sorting to `{2, 3, NaN}` — the block
move that carries a NaN without any comparison involving it being true — is
pinned too.

**The five frozen A7 oracle cases, which are the acceptance test.** All five are
now compared line for line against the retained Release reports, and every one
differs *only* in native difference 5's `nan` / `-nan` spelling:

| oracle case | before | after |
| --- | --- | --- |
| `c_nan_one_s` | compared, 9 NaN-spelling lines | unchanged, still 9 |
| `c_nan_two_s` | compared, 10 NaN-spelling lines | unchanged, still 10 |
| `c_nan_then_finite_s` | **refused**; retained as evidence, not compared | compared, 7 NaN-spelling lines, nothing else differs |
| `c_finite_then_nan_s` | **refused**; retained as evidence, not compared | compared, 7 NaN-spelling lines, nothing else differs |
| `c_zero_swapped_s` | compared, 9 NaN-spelling lines **plus 4 order-statistic lines** (native difference 6) | compared, 9 NaN-spelling lines, **0 order-statistic lines** |

`consensus_a_signed_zero_sample_is_ordered_by_the_total_order` asserted a
divergence and is renamed to
`consensus_a_signed_zero_sample_keeps_the_release_builds_order`, which asserts
equality with **both** members of the measured pair — the port reproduces the
swapped file's report and the unswapped file's report rather than collapsing
them onto one. `assert_report_but_the_nan_spelling` is a stronger comparison
than a line count: it fails if any line outside the NaN class differs at all, if
the number of NaN-spelled lines changes, if a `-nan` appears where the reference
has a number, or if the crate starts or stops writing the sign. **Native
difference 6 is closed.**

Byte-for-byte equality including the sign is not reachable from this wave and
was not attempted: `src/format/file_info/text_format.rs` was not touched, and
its `nonfinite` rule is part 3 of the promotion bullet, which still needs A2's
oracle row re-captured against the Linux Release build.

**One cost, measured rather than assumed — and larger than this wave reported.**
`sort_ascending` builds a permutation with an interpreted introsort and a
closure per comparison where it used to call `slice::sort_by`. This wave put
that at "2.3x at 1,000 values rising to 5.0x at 1,000,000" from a table with no
committed harness and no retained data, and read the blast radius as narrow
because `FileInfo`'s `summarize` works on "bounded samples". **Both were wrong.**
Re-measured with the committed harness of
[BENCHMARKS](BENCHMARKS.md) §8, the public entry point cost **19.1x** wall clock
and **2.9x** peak memory at ten million values; and
`src/format/file_info/peaks.rs:706` and `:714` hand `summarize` every MS1 peak
intensity in the file, bounded only by
`FileInfo::MAX_STATISTICS_VALUES = 1 << 27`, which is not a small sample but the
largest one the tool accepts. **Lead decision D17 closes it**, and the entry
point is now within 2 % of the library sort on any sample whose permutation
cannot be observed — which is all 809 of the repository's own corpus samples.

**What this wave did not close, and says so at the item.**

- `compute_rank` and `rank_correlation_coefficient` still refuse a NaN. Their
  `std::sort` (`:829-830`) is a lambda comparing `std::pair::second`, not the
  default `operator<`, and a NaN additionally defeats their *relative tie test*,
  whose two comparisons are both false against a NaN and which would therefore
  merge every block the NaN touches. That is a second, independent behaviour
  that no oracle row measures. For a NaN-free range the two sorts agree on the
  ranks anyway, because `operator<` and `total_cmp` differ only on `±0.0`, which
  the tie test makes one block either way.
- The public `_sorted` entry points (`median_sorted`, `quantile1st_sorted`,
  `quantile3rd_sorted`, `quantile`) still return `UnsortedData` for a NaN. They
  do not sort; they verify a *caller's* claim to have sorted, and there is no
  way to know which permutation a caller who asserts "already sorted" about a
  NaN-bearing range meant.
- `src/format/file_info/consensus.rs` still generates a NaN in plain Rust
  arithmetic: `it_aad += it_ratio` is `(-inf) + (+inf)` for the
  `a7_cons_nan_one` fixture. Every value `statistic_functions` produces is now
  host-independent; this one is not, and the spelling step has to take it with
  it. Part 2's brief scoped the rewrite to `statistic_functions.rs`, so this is
  reported rather than changed.
- `pearson_correlation_coefficient` and `matthews_correlation_coefficient` were
  left on plain arithmetic on purpose. Both substitute an explicit `f64::NAN`
  for a division the source actually performs — a divergence that predates this
  work and is documented at each item — so making only their *other* operations
  bit-faithful would leave that substituted NaN as the single host-shaped value
  in the result. The crate's bit-faithful Pearson is
  `analysis::feature_finder_picked::scoring::source_pearson`, which the picked
  feature finder measured against the Release build; that the two exist side by
  side is worth a reviewer's attention and converging them needs an oracle row.

**CPP-347** is rewritten as closed-by-reproduction: the port now reproduces the
defect rather than refusing it, and the underlying C++ defect —
`Math::SummaryStatistics` reading order statistics positionally out of a range
whose order the comparison did not determine — still stands and is still worth
raising with maintainers.

## Wave 8: the reader's round trip, the featureXML writer, the picker consumers, A7 and the citation checker (2026-09-19)

`integrate/wave8` merges `fix/reader-round-trip` (`42064c6`),
`fix/featurexml-nonfinite` (`511fafa`), `fix/picker-noise-consumers`
(`8d3e052`), `port/a7-fileinfo` (`7801a05`) and `tools/citation-checker`
(`f9c692f`) onto `main` `2ba9c1d`. All five merged without a conflict, and the
merged tree changes **255 files** against `main` (+29,435/−580), which is the
237 paths the five lanes touch plus 18 this pass adds: the 11 integrator-owned
records, and the seven files — one fixture and six frozen Release reports — the
A7 verifier's major finding required. See [FILE_INFO_A7_SUPPORT](FILE_INFO_A7_SUPPORT.md),
[MZML_HEADER_SUPPORT](MZML_HEADER_SUPPORT.md),
[MS_DATA_WRITING_CONSUMER_SUPPORT](MS_DATA_WRITING_CONSUMER_SUPPORT.md),
[FEATUREXML_SUPPORT](FEATUREXML_SUPPORT.md),
[CHROMATOGRAM_PICKING_SUPPORT](CHROMATOGRAM_PICKING_SUPPORT.md),
[ITERATIVE_PICKING_SUPPORT](ITERATIVE_PICKING_SUPPORT.md) and
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-8-status).

**Disjointness was checked, not assumed.** Three lanes touch `src/format/`, so
the name-only diffs of all five branches against `main` were listed and
compared pairwise: 237 paths, **zero** appearing in more than one lane.
`src/format/` splits cleanly — the reader lane owns `mzml.rs`,
`mzml_header/read.rs` and `ms_data_writing_consumer.rs`, the featureXML lane
owns `featurexml.rs` and `identification_xml.rs`, A7 owns `file_info.rs` and
`file_info/*`. Confirmed after the fact rather than left as an intention:
every one of the 237 appears in the merged diff, and the 18 paths that are not
a lane's are this pass's own.

### The five verdicts, and what this pass had to finish

| lane | rounds | closing verdict | carried into this pass |
|---|---|---|---|
| `fix/reader-round-trip` | 1 | **approve_with_notes**, 2 minors | nothing; both minors were applied on the branch |
| `fix/featurexml-nonfinite` | 2 (3 majors, 5 minors → 1 minor) | **approve_with_notes** | nothing |
| `fix/picker-noise-consumers` | 2 (3 majors, 4 minors → 2 minors) | **approve_with_notes** | nothing |
| `port/a7-fileinfo` | 3 | **changes_required**, 1 major + 5 minors | **all six applied here** |
| `tools/citation-checker` | 3 | **approve_with_notes**, 4 minors | notes only; recorded below |

A7 is the one lane that did not close clean. Its major was **not taken on the
verifier's word**: the fixture was generated from the lane's own oracle
generator, the oracle was re-run on `ibminode06` against the Release build at
the pins, and the port's side was measured from a build of the merged tree.

**A7's major — the same `std::sort` boundary, without a NaN, and silent.** The
A7 close round escalated to the lead a sample its own arithmetic can put a NaN
into, and told the lead the boundary was "pinned by
`consensus_nan_in_the_statistics_sample`, so either decision is a visible change
rather than a silent one". The verifier measured one fixture further out and
found a second instance that was *not* visible. `FileInfo.cpp:2310-2311` fills
the `Intensity ratios` sample **before** the inversion at `:2312-2315`, so one
consensus feature with sub-features of intensity `-0.0` and `0.0` puts both
zeros into that sample. `operator<` calls them equivalent — `-0.0 < 0.0` and
`0.0 < -0.0` are both false — so `std::sort`'s strict-weak-ordering precondition
**holds**, nothing is undefined, only the permutation is unspecified, and
libstdc++ leaves a range this size as it found it. `front()`, the quantiles and
`back()` at `:952-956` are positional reads, and `ostream` writes `-0` for a
negative zero. `sort_ascending` **then ordered** by `f64::total_cmp`, which puts
`-0.0` first, so the port printed one of the two answers for both file orders —
four lines the Release build does not print, on an input it accepts, with no
oracle case, no frozen report and no test. (Every sentence in this subsection is
wave 8's state of play; D16 closed it the next day, as the paragraph at the end
records.)

Measured here, both ways. `a7_cons_zero_swapped.consensusXML` is
`a7_cons_nan_one` with its two sub-feature intensities exchanged, generated by
`../oracle/a7-fileinfo/scripts/make_a7_fixtures.py`; the oracle gained
`c_zero_swapped`, `c_zero_swapped_s` and `c_zero_swapped_all` and was re-run on
`ibminode06` (75 cases, two runs, `reproduced: true`). On the Release build:

| line | `a7_cons_nan_one` | `a7_cons_zero_swapped` |
| --- | --- | --- |
| `minimum:` | `-0` | `0` |
| `lower quartile:` | `-0` | `0` |
| `upper quartile:` | `0` | `-0` |
| `maximum:` | `0` | `-0` |

The port **then printed** the left column for both. That the re-run is *additive* is
checked rather than asserted: **all 110 frozen expectations rebuild byte for
byte out of the new manifest**, of which 104 predate the re-run, and the
reference tool's sha256 `5d82c8a7…1172dc` is the same binary A6 recorded — the
same one A6 and the A7 lane ran, which is an independent check that all three
packages measured the same Release build.

Recorded as **native difference 6**, not as a refusal. Refusing would widen a
refusal in `sort_ascending` — which every `SummaryStatistics` caller in the
crate consumes — to an input the Release build handles in bounds, stably, and
with its own precondition satisfied, and that is the lead's decision rather than
an integrator's. `consensus_a_signed_zero_sample_is_ordered_by_the_total_order`
asserted the two bare reports byte for byte, the four lines on which the two
Release reports disagree, and that each of the port's four lines is the Release
build's own line for the unswapped file — so both spellings were measured and
none was derived from Rust output. The lead's open question about reproducing an
unspecified `std::sort` permutation then had **both** instances in front of it,
and `CPP-347` was rewritten around the general defect: order statistics read
positionally out of a range whose order the comparison did not determine. A fix
that only filters non-finite values leaves half of it in place.

**Closed the next day** by lead decision D16 and the shared-math wave: the four
lines are no longer a divergence, the test is renamed
`consensus_a_signed_zero_sample_keeps_the_release_builds_order` and asserts
equality with both retained Release reports, and native difference 6 is gone.
The wave section at the top of this document records the measurement.

A7's five minors were applied too. The one that matters beyond wording: the
claim "libstdc++ compares every pair involving a NaN false and therefore moves
nothing" was stated unconditionally in eight places and is false above
libstdc++'s insertion-sort threshold, where `__introsort_loop` does move the
NaN — a 20-element sample `{NaN, 2..20}` prints `minimum: 2` and `median: -nan`.
Scoped in all eight, including the oracle driver and its manifest, which was
re-emitted.

### Lead decisions D1-D17

The full text of D1-D13 is in
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-5-status), which is wave 5's
own list; every decision after it is recorded here in full. Two are new this
wave, and they collided: **both** `fix/reader-round-trip` and
`fix/picker-noise-consumers` numbered their decision D14. The reader lane's
number is cited in nine committed files — `MZML_HEADER_SUPPORT.md`, native
difference 12 of `TOPP_PEAK_PICKER_HI_RES_SUPPORT.md`,
`MS_DATA_WRITING_CONSUMER_SUPPORT.md`, the `ReadOptions` rustdoc in
`src/format/mzml.rs`, `tests/mzml_source_file_round_trip.rs`,
`tests/mzml_header_leniency.rs`, `tests/topp_peak_picker_hi_res.rs` and the two
manifests `tests/data/mzml_source_file_round_trip_provenance.json` and
`tests/data/mzml_header_leniency_provenance.json` — and the picker lane's only
in its handover text, so the reader lane keeps D14 and the
picker lane becomes **D15**. Nothing in the tree had to change.

- **D14** the mzML reader accepts a dangling `sourceFileRef` under
  `source_dangling_references`, as the Release build does, and keeps today's
  refusal as the strict default — because the port had otherwise introduced a
  file only the C++ reader would read.
- **D15** a picker's internal noise estimate reproduces the source
  **unconditionally** where the signal it reads is picker-generated or already
  validated (`PeakPickerChromatogram`, whose `snt_` reads the smoothed trace
  under `corrected`), and **follows the picker's own profile** where the
  estimator reads the caller's data (`PeakPickerIterative`, whose
  `snt.init(input)` reads the raw spectrum, as `PeakPickerHiRes` already did). A
  `PickingCompatibility` flag that the port cannot yet honour faithfully leaves
  its refusal in place in both profiles rather than returning different results
  under a flag that claims source behaviour.

**D16 and D17 postdate this wave** and are kept here rather than in a section of
their own, because this list is where the full text of every decision after
wave 5 lives. D16 belongs to the shared-math wave of the next day; D17 to the
repair round that followed it.

- **D16** (shared-math wave; taken by the lead on **2026-09-19**, in the session
  that briefed that wave). **Reproducing an unspecified `std::sort` permutation
  is in scope.** The question the `sort_ascending` bullet of
  [the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-5-status) and the "Still
  open" section of this document carried, and which the "Where a NaN lands"
  section of
  [STATISTIC_FUNCTIONS_SUPPORT](STATISTIC_FUNCTIONS_SUPPORT.md) and section 5.2
  of [FILE_INFO_A7_SUPPORT](FILE_INFO_A7_SUPPORT.md) left open — whether the
  port should reproduce a permutation the C++ standard leaves unspecified — is
  decided **yes**, on four grounds:
  1. **The port already does it.** `src/math/source_sort.rs` is a
     comparison-by-comparison, move-by-move port of the GCC 14.4.0 libstdc++
     introsort and `stable_sort`, validated tier 1 against two oracle drivers
     (`../oracle/ffap-instr-completion`, `../oracle/ffap-complete-fix1`) over
     2,272 inputs carrying ties, signed zeros, infinities and four NaN bit
     patterns. Using it in `sort_ascending` promotes executed evidence; it does
     not gamble on new behaviour.
  2. It is **D1-compliant by construction**: it reproduces the in-bounds,
     deterministic, measured behaviour and refuses exactly where the introsort's
     unbounded partition and final-insertion loops read outside the vector.
  3. **Refusing is worse.** The signed-zero half (native difference 6) is
     ordinary finite data the Release build summarises without complaint; a
     blanket refusal in `sort_ascending` would turn it away, and
     `sort_ascending` is what every `SummaryStatistics` caller in the crate
     consumes.
  4. The **pin risk is already managed**: `source_sort` names the sha256 of
     every libstdc++ header whose algorithm it reproduces, so a toolchain change
     is detectable rather than silent.

  Filed here rather than in the work packages, where the shared-math wave
  mistakenly appended it to wave 5's `D1-D13` list; that list is wave 5's own
  and now says so.

- **D17** (shared-math repair round, **2026-09-20**). **`sort_ascending` may
  take a proved-equivalent fast path.** D16's faithful sort measured **19.1x**
  wall clock and **2.9x** peak memory at n = 10,000,000 against the library sort
  it replaced ([BENCHMARKS](BENCHMARKS.md) §8), and
  `src/format/file_info/peaks.rs:706` and `:714` hand `summarize` **every MS1
  peak intensity in the file**, bounded only by
  `FileInfo::MAX_STATISTICS_VALUES = 1 << 27`. That is a real regression on
  `FileInfo -s` over a routine LC-MS run.

  The decision: `sort_ascending` sorts **in place with `f64::total_cmp` and no
  allocation at all** when the sample contains **no NaN** and **not both zero
  spellings**, and runs the libstdc++ permutation otherwise.

  This is not a fidelity compromise, and the reason is a theorem rather than a
  preference. Without a NaN, `operator<` on `f64` is a strict weak ordering
  whose incomparability relation is numeric equality; two numerically equal
  non-NaN doubles have the **same bits**, with `-0.0 == +0.0` the one exception
  in the format. So unless the sample holds both spellings of zero, every
  equivalence class is a set of bit-identical values, the sorted *sequence* is a
  function of the multiset alone, and `std::sort`'s choice among equivalents is
  unobservable — any correct sort writes what the Release build writes. The
  guard is therefore exactly one O(n) pass testing `is_nan()` and whether both a
  negative and a non-negative zero occur.

  The argument is not what the port rests on.
  `both_paths_agree_bit_for_bit_wherever_the_fast_one_is_taken` runs **both**
  paths over the same adversarial samples — both zeros, both infinities,
  subnormals, `DBL_MAX`, signalling and negative NaNs, heavy duplication,
  ascending, descending, organ-pipe and sawtooth shapes, and raw random bit
  patterns, at 21 lengths spanning libstdc++'s 16-element `_S_threshold` and its
  heapsort fallback — and compares the results bit for bit, NaN payloads
  included. Deleting either half of the guard makes it fail. Over every mzML
  fixture in `tests/data`, 809 of 809 statistics samples take the fast path
  (§8.4), so the cost D16 accepted was being paid on every sample and bought
  nothing on any of them.

**The lead's two decisions of this round**, both taken inside the A7 lane and
both recorded here because they shaped what landed:

1. **Native difference 5 takes route (b): document and pin the NaN spelling.**
   The FileInfo text layer spells every NaN `nan` where glibc spells a sign-bit
   NaN `-nan`. Spelling the sign honestly would first have to make the *value*
   host-independent, which needs the x86_64 emulation promoted out of
   `analysis::feature_finder_picked::scoring` into shared math, an
   x86_64-faithful `variance_with_mean` built on it, and A2's oracle row
   re-captured against the Linux Release build instead of the macOS SDK. That is
   a cross-cutting change to shared math which landed ports already consume, so
   it is carried forward as its own wave rather than done inside a package. The
   three parts and the owner are in the work packages.
2. **Under D1 the one-value NaN sample is reproduced, not refused.** The lane
   was told to accept the sample and print what the Release build prints. It did
   so for the two shapes whose set of possible outputs has exactly one member —
   a one-value sample and an all-NaN sample — and raised the rest rather than
   quietly extending the instruction. This pass agrees with that reading and
   says why in `docs/FILE_INFO_A7_SUPPORT.md` section 5.2: the instruction is
   well defined exactly where the permutation cannot be observed.

### Decision D14 in one paragraph, because it reverses an earlier one

Wave 7 taught `MSDataWritingConsumer` to reproduce the
`sourceFileRef="sf_sp_<s>"` the source's `writeSpectrum_` writes for a record
the streamed header cannot declare (`CPP-172`), while the reader kept refusing a
dangling `sourceFileRef`. On an input with per-record source files the port
therefore wrote a file neither it nor a strict reader would read back, while the
C++ reader read both its own output and the port's. That is a defect **this
project introduced**, not one it inherited. D14 reverses the
"`sourceFileRef` … stays strict" line drawn by work package P2-MZML-LENIENCY,
deliberately, and the record says so in `docs/MZML_HEADER_SUPPORT.md`, in
`tests/data/mzml_header_leniency_provenance.json` and in
`tests/mzml_header_leniency.rs`, where `sourceFileRef` moved out of the "stays
strict" list into an explicit strict-default / lenient-under-the-switch
assertion. What the Release build does was **measured with a reader probe**
rather than inferred from `FileInfo`'s exit code: a spectrum's dangling
reference leaves the record a default `SourceFile` and warns once per
occurrence; a chromatogram's default-constructs, inserts and sets an empty one
silently; a scan's and a precursor's set two present, empty metadata keys. The
round trip is executed both ways — twelve files, two implementations, three
inputs, both `-processOption` modes, every file read by both with exit 0.

### D1 applied both ways, this wave

**Refused under D1** (each measured on the Release build): the three
out-of-bounds `std::vector` accesses of the A7 branches — a consensus
sub-feature whose map index is at or beyond the column-header count
(`c_no_headers`, `c_mapindex_high`, `c_ids_one_based`, `c_cid3`; `CPP-341`), an
identification file with no protein run (`id_no_runs`; `CPP-342`) and a peptide
identification with no hit (`id_empty_hitlist`, `id_empty_hitlist_ok`;
`CPP-343`). Two of the four consensus cases and both identification cases end in
SIGSEGV; the other two exit 0 with a wrong number, which is why D1's
out-of-bounds clause rather than its measured-and-repeatable clause governs
them. Worth recording with `id_no_runs`: the refusal a user actually gets is the
shared idXML reader's ("idXML needs at least one IdentificationRun", tool exit
3), not the branch's own guard, which no input can reach.

**Reproduced under D1**, where the source is defined and in bounds: the
consensusXML `-s` quality sample that the source pre-sizes and then appends to,
so it is twice as long as it should be with a zero half, which upstream's own
`FileInfo_7_output.txt` records (`CPP-344`); the peak-file arm `-m`, `-p` and
`-s` give an mzIdentML input (`CPP-345`); the non-finite statistics of
`c_zero_intensity_s`; the two answerable NaN-sample shapes of `c_nan_one_s` and
`c_nan_two_s`; the featureXML writer's `inf`, `-inf` and `NaN`; the iterative
picker's division by a negative integrated intensity; and the chromatogram
picker's bin-index conversion for a histogram quotient outside `int` range.

**A NaN next to a number is refused, and that refusal is two different things
depending on the sample**, which is worth keeping apart because only one of them
needs a decision:

- *a D1 refusal*, for a sample holding a NaN and **two or more distinct
  numbers**. Transitivity of incomparability fails there — `1 ~ NaN` and
  `NaN ~ 3` while `1 < 3` — so `std::sort`'s strict-weak-ordering precondition
  is violated outright and the call is undefined. D1 refuses undefined
  behaviour, and no lead decision is needed.
- *a deferral that is explicitly NOT a D1 refusal*, for a NaN next to **at most
  one distinct number**. The precondition holds there, every element is
  equivalent, and only the permutation is unspecified; the Release build exits 0
  and prints a stable report, so D1 would have the port reproduce it. It did
  not, because the `minimum`, quartile and `maximum` lines are positional reads
  of a range `std::sort` was free to leave in any order: `c_nan_then_finite_s`
  and `c_finite_then_nan_s` hold the same two consensus features in opposite
  file order and disagree on exactly those four lines. Reproducing them means
  porting libstdc++'s `std::sort` permutation into `sort_ascending` in shared
  math, which every `SummaryStatistics` caller consumes. **The lead decided it:
  D16, in the shared-math wave of 2026-09-19 — reproducing the permutation is in
  scope, and both cases are now reproduced.** See `CPP-347`.

Both were refused by the same message, which is why the distinction is recorded
here rather than left to be read off the code. **Neither is refused any more**:
the first class was never a D1 refusal in the sense D1 means — `std::sort` runs
deterministically on such a range, it is only the *standard's guarantee* that is
void — and D16 says to reproduce the executed, measured behaviour and to refuse
only where libstdc++ reads outside the vector. The refusal that remains is that
one.

That deferral and **native difference 6** were the two halves of one question,
scoped differently on purpose. Both are samples whose elements `std::sort` calls
equivalent, in both the precondition holds, and in both the order statistics are
positional reads — but the NaN half was refused and the signed-zero half was
accepted with its divergence pinned. The difference was what refusing would
cost: refusing a NaN sample turns away an input no caller has a use for, while
refusing a signed-zero sample would turn away ordinary finite data the Release
build summarises without complaint. **The shared-math wave closed both**, and
`CPP-347` is written around what they share rather than around the NaN.

### The featureXML round trip, and one divergence in each direction

`NumericFormatting::appendNumeric` writes `NaN` for a NaN of either sign and
`inf`/`-inf` for an infinity (`CONCEPT/Detail/NumericFormatting.h:27-35`), and
`StringUtils::toDouble` reads all three back, so the featureXML dialect now does
the same for a feature's position, intensity, qualities, overall quality,
width/`FWHM` and every `float` and `floatList` meta value, through the
crate-private `MetaValue::source_float` (decision D13). The public
`TryFrom<f64>` and `MetaValue::validate` still refuse a non-finite value, and so
do the shared map and identification codecs, a hull point and a finite value
`f32` cannot hold. A literal `toDouble` cannot convert at all — `1e999`,
`banana`, `inf.0` — stays refused, which **agrees** with the source wherever the
literal is an attribute, whose `ConversionError` leaves the parse, and
**diverges** from it in an element's text, where `asDouble_` logs a non-fatal
line and keeps `0.0`, so the Release build loads a document this port refuses;
an underflowing literal diverges the other way. Both directions are recorded
with their executed evidence rather than one of them. Together with a failed
store taking the source's write-side arm (`CANNOT_WRITE_OUTPUT_FILE`,
`TOPPBase.cpp:430-435`, exit 5) instead of being announced as a read failure,
this closes **TOPP native difference 16**: `FeatureFinderCentroided` on a
retention-time-scaled input exits 0 and writes all 1,263 values the Release
build writes. `CPP-327`'s Rust-handling paragraph is corrected to match.

### The citation checker, and what a green run from it means

`tools/check_source_citations.py` resolves C++ source citations against the pins
this repository already declares and reads them back. At this integration head
it exits **0** over the whole tree: **3,344 citations resolved, 102 confirmed
against code quoted beside them, 192 ambiguous, 0 problems**, in about four
seconds. (The polish round below takes those to 3,346 / 104 / 196: two more
citations because a docstring there quotes its source verbatim, and four more
ambiguous because the guard learned to count two paths of one pin.) Getting
there took the six substitutions the lane asked for in `OpenMS_CPP_ISSUES.md`
— four in one `MSSpectrum.cpp` block whose lines had shifted by six, and two
single-line citations, one in `MzTabFile.cpp` and one in `SVOutStream.cpp`,
that had come to name blank lines — plus an editorial seventh, where a
`MzTabFile.cpp` citation named a real line but the wrong loop: the
`best_search_engine_score` loop rather than the score-outer row loop the entry
is about. Each new line number was read back at the pin, and the stale numbers
are deliberately not repeated here, because writing one beside its file name in
prose is exactly the defect this checker exists to catch.

**It earned its place on the way in.** Run against the merged tree it found a
citation defect A7 shipped and three review rounds missed:
`src/format/file_info/identifications.rs` cited `FileInfo.cpp:1354` for a
quotation that sits at `:1347`. Fixed by naming the three lines the sentence
means — `:1347` the guard, `:1352` the reference, `:1354` the read.

What a green run does **not** mean is stated by the tool itself and repeated
here, because a count that reads as more than it measures is the failure mode
this lane exists to prevent. On the tree as it now stands, only 104 of 3,346
citations are *confirmed*; the
rest are checked for existence, for a non-blank single line, and for annotation
blocks. Most citations in this repository paraphrase the source instead of
reproducing it, and a paraphrase cannot be read back. **The A6 defect that
motivated the lane is still not caught**: all eight places that write `:290-293`
for code at `:280-283` paraphrase, and the replay at the lane's final head
produces no finding. What the tool covers is the *class* — the transcribed block
in `CPP-337` that reproduces the same code is read line by line, and shifting it
fires. A further 196 citations name a file more than one pin carries with nothing
to tell them apart; one pin answered each, and that count is how often the pin
that answered may have been the wrong file of the right name. It is 196 rather
than the lane's 79 because this pass's own records cite those file names too,
and because the polish round counts a name one pin carries at two paths as
ambiguous as well. The lever that
would raise the confirmed fraction is a convention — quote the source verbatim
in the code span beside the citation — not a cleverer checker.

The lane's second task is separate and is fixed: `tools/check_core_sdk.py` gave
a different verdict per worktree, because `check_external_artifacts` asked
`ROOT.rglob` — the working directory — rather than the repository, so the
gitignored `.reference/` checkouts answered the "no C++ left behind" search.
Reproduced both ways in this pass: from the **main** worktree at `2ba9c1d` the
checker fails with `Exception.h still present in repo as [nine .reference/
paths]`, rc 1; with the merged tool and the same `ROOT` it passes,
`Verified 2092 distinct current source/registration/reference files`. The main
worktree was not modified — `git status --porcelain` was empty before and after.

### This pass's gates

All on `kim` through the gate script, slot `integ-w8`, one gate at a time on one
slot, detached and polled; logs in the session scratchpad under
`integ-w8-logs/final/`. Every figure below is summed from **all** the log's
`test result:` lines and cross-checked against an anchored recount of the
`... ok` and `... ignored` lines between them, so a dropped or interleaved line
cannot hide a binary. It earned its keep on the full run: the triple sum reads
5,477 and the anchored recount 5,475, because the shared ssh capture swallowed
two `... ok` lines, and taking the larger of each pair repairs exactly that.

| Gate | Result |
|---|---|
| `+1.85.0 check --locked --all-features --all-targets` | exit 0 (MSRV 1.85) |
| `clippy --locked --all-features --all-targets -- -D warnings` | exit 0 |
| `clippy --locked --no-default-features --all-targets -- -D warnings` | exit 0 |
| `doc --locked --all-features --no-deps`, `RUSTDOCFLAGS=-D warnings` | exit 0 |
| `test --locked --all-features --doc` | exit 0, 73 + 3 = **76 doctests**, unchanged from `main` |
| `test --locked --all-features --all-targets --no-fail-fast` | exit 0, **5,477 passed / 0 failed / 21 ignored** over 359 result lines |
| `test --locked --no-default-features --no-fail-fast` | exit 0, **3,660 passed / 0 failed / 3 ignored** over 335 result lines |
| `test --locked --no-default-features --features mzml --test mzml_source_file_round_trip` | exit 0, 5 passed |
| `test --locked --no-default-features --features consensusxml,idxml --test file_info_a7` | exit 0, 34 passed |
| `test --locked --no-default-features --test file_info_a7 --test picker_noise_consumers --test statistic_functions` | exit 0, 11 + 19 + 27 passed |

The doctest gate is run separately on purpose: `--all-targets` does not cover
doctests, so a battery without it passes vacuously on that slice.

**The full-suite difference from `main` is accounted for lane by lane**, not
merely noted. `main` is 5,401 passed / 0 failed / 21 ignored; this head is
5,477 / 0 / 21, a difference of +76 passed and **no change in
ignored**. The `.rs` diff against `main` adds exactly 76 `#[test]` functions
(5,430 → 5,506) and removes none, and they attribute cleanly: `fix/reader-round-trip` 7,
`fix/featurexml-nonfinite` 7, `fix/picker-noise-consumers` 19,
`port/a7-fileinfo` 42, this pass 1 (the signed-zero test), `tools/citation-checker`
0 Rust tests and 49 Python ones. Each lane's own reported full-suite figure
reconciles against the same base: 5,401 + 7 = 5,408 for the reader and the
featureXML lane, 5,401 + 19 = 5,420 for the picker lane, 5,401 + 42 = 5,443 for
A7 — the three numbers those lanes reported at their own heads — and the
citation lane reported 5,401 unchanged.

Two gates ran twice, and the reason is recorded rather than smoothed over. The
first battery's `doc` gate failed, rc 101, on a broken intra-doc link this pass
introduced: `[`sort_ascending`]` in the new "Signed zeros" section of
`src/math/statistic_functions.rs` names a private item, which rustdoc cannot
resolve and `-D warnings` therefore rejects. Fixed by writing the name in plain
backticks, and the whole battery re-run at the final head; no other gate was
touched by the change. `git diff --name-only a6c7f34..1e1833a` names three files
committed after the gate head and before the polish round: `docs/VALIDATION.md`,
`docs/EARLY_TOPP_WORK_PACKAGES.md` and `SOURCE_PROVENANCE.json`. Two are
Markdown; the third is a JSON record whose change is one prose sentence inside
`external_reference_note`, and `grep -rn SOURCE_PROVENANCE --include='*.rs'`
finds no Rust reader for it. So nothing the compiler reads moved after the gate
head — which is the claim that matters, and it is narrower than "Markdown
only".

Locally, from the integration worktree: `cargo fmt --all -- --check` exit 0;
`tools/check_doc_coverage.py` 4563/5918 = 77.1 %, recorded with `--write`;
`tools/check_module_cycles.py` 64 cross-module edges and 13 mutually-dependent
pairs, unchanged from `main`; `check_core_sdk.py`, `core_sdk_coverage.py`,
`test_core_sdk.py`, `test_core_sdk_coverage.py`,
`check_schema_feature_graph.py`, `test_source_citations.py` and all ten
`generate_*.py --check` / projection / probe checkers exit 0;
`check_source_citations.py` exits 0 with 3,344 citations resolved and 102
confirmed. Every changed JSON parses (14 files), and
`.github/workflows/rust.yml` parses as YAML. Determinism was re-checked here
rather than taken on report:
`parallel_determinism` 5 passed and `topp_threads` 7 passed, including
`sums_are_bit_identical_across_thread_counts` and
`map_collect_preserves_input_order_at_every_thread_count`.

No C++ was built or run in this pass except the A7 oracle re-run on
`ibminode06`. `ibminode05` was never contacted.

### The polish round on top of this head

Twelve minors were carried into a polish round after the integration head, and
each was checked against this head before anything was changed.

**Four were already applied** by the integrator while merging — the whole A7
group: the mechanism claim scoped in all eight places, integrator request 8's
refused class split in two, the oracle case note rewritten and its manifest
re-emitted, and the garbled test comment finished. They are named here so that
a re-reader does not go looking for them. The mechanism claim was nevertheless
sharpened again, because checking that the applied scoping was *true* turned up
a second condition it did not carry; that is the first bullet below. Three further findings turned out to
be about handover text rather than about the repository, and no committed file
carries them: the A7 handover's commit list, the citation lane's commit count,
and the ledger claim answered above.

**Eight needed work**, and two of those were not wording. **The chromatogram
`sourceFileRef` refusal is now executed** rather than only described, and
**the citation checker's ambiguity guard now counts per file rather than per
pin** and counts only where it counts `checked`, with `Pins.packaged` and both
directions of the subset relation pinned by tests that fail when the fix is
mutated away. The other six are the counts and wording below.

Three figures were re-measured rather than carried, and two of them moved
against what the round was told:

- The libstdc++ mechanism behind the NaN refusal was read out of the headers
  the Release build was compiled with, not inferred. `_S_threshold` is 16 and
  `__introsort_loop` runs only above it, so at sixteen elements or fewer
  `std::sort` is a single `__insertion_sort` pass — but a NaN can still be
  carried by a block move there, which the "property of the size" wording did
  not allow for. `{3, NaN, 2}` sorts to `{2, 3, NaN}` at three elements.
- The citation mutation replay was recomputed from scratch at both heads. The
  lane report's "caught 72, missed 2" does not reproduce; the four misses the
  verifier found do, at the lane head and here, and the count is a property of
  the shift chosen rather than a bound on the checker.
- The removed-assertion count was not six. It is **eight**; the two the audit's
  own grep missed are written as `outcome.assert_*` method calls.

**Gates at the polish head**, on `dax`, slot `w8-polish`, detached and polled,
through `/scratch/kohlbach/openms-rs-env.sh` for `pkg-config` and `libxml2`:

| Gate | Result |
|---|---|
| `fmt --all -- --check` | exit 0 |
| `clippy --locked --all-features --all-targets -- -D warnings` | exit 0 |
| `clippy --locked --no-default-features --all-targets -- -D warnings` | exit 0 |
| `+1.85.0 check --locked --all-features --all-targets` | exit 0 |
| `doc --locked --all-features --no-deps`, `RUSTDOCFLAGS=-D warnings` | exit 0 |
| `test --locked --all-features --doc` | exit 0, 73 + 3 = **76 doctests**, unchanged |
| `test --locked --all-features --all-targets --no-fail-fast` | exit 0, **5,478 passed / 0 failed / 21 ignored** over 359 result lines |
| `test --locked --no-default-features --no-fail-fast` | exit 0, **3,660 passed / 0 failed / 3 ignored**, unchanged |

The no-default suite does not move because the new test is behind the `mzml`
feature. 5,478 is 5,477 **+1**, and the one is
`a_chromatograms_source_file_ref_is_refused_under_both_policies`. The triple
sum and the anchored recount of `... ok` lines agree at 5,478 this time, with
0 `FAILED`, 0 `panicked`, 0 `failures:` and 0 lines starting `error`. Locally:
every `tools/*.py` checker exits 0, `core_sdk_coverage.py --write` leaves no
diff, `check_source_citations.py` reports 3,346 resolved / 104 confirmed / 196
ambiguous / 0 problems, and `test_source_citations.py` runs 52 tests. The round
adds one Rust test and eight Python ones, removes none, and removes no
assertion; `unsafe` and `#[ignore]` counts and `Cargo.toml`/`Cargo.lock` are
untouched.

The battery ran at `9cd037c` from a clean worktree. Everything committed after
it is Markdown — and this time that claim is the narrow one, checked the way
the finding above asks: `git diff --name-only 9cd037c..HEAD` returns
`docs/MZML_HEADER_SUPPORT.md` and this file, both Markdown, and nothing with an
extension the compiler reads.

The oracle directory under `../oracle/a7-fileinfo` was deliberately **not**
re-emitted, and the justification has to be narrower than it first read. The
note's first half is scoped to the two-element samples it describes and is
exact there. Its second half — "that 'moves nothing' is a property of the SIZE,
not of the NaN" — over-generalises in exactly the way the in-repo places did:
below the insertion-sort threshold a block move can carry a NaN too, as
`{3, NaN, 2}` sorting to `{2, 3, NaN}` shows, so the property is of the size
**and** the arrangement. The corrected statement lives in
[STATISTIC_FUNCTIONS_SUPPORT](STATISTIC_FUNCTIONS_SUPPORT.md) and
[FILE_INFO_A7_SUPPORT](FILE_INFO_A7_SUPPORT.md); the driver keeps the older
sentence because re-emitting it would move a sha256 already registered in
`SOURCE_PROVENANCE.json`, and the sentence changes no case, no fixture and no
recorded value. It is the next re-emission's to fix.

### Ignored tests

**28 `#[ignore]` attributes in the tree, all of them in `tests/` and none in
`src/` — byte-identical to `main`**, file for file and count for count, across
the same ten files (`featurexml` 1, `file_info` 6, `fuzzy_string_comparator` 1,
`gauss_trace_fitter` 1, `lm_budget_differential` 1,
`lm_eigen_path_differential` 1, `mzml_reader_scale` 8, `mzml_writer_scale` 3,
`topp_baseline_filter_edges` 3, `topp_threads` 3). No test was ignored, skipped,
deleted or weakened by this wave; no tolerance was widened; no expected value
was derived from Rust output. The `.rs` diff against `main` adds **0**
occurrences of `unsafe` and **0** new `#[ignore]`, and adds 76 `#[test]`
functions (5,430 → 5,506), which reconcile lane by lane: reader 7, featureXML 7,
picker 19, A7 42, this pass 1, the citation lane 0 Rust tests and 49 Python ones
(`tools/test_source_citations.py` 44, `tools/test_core_sdk.py` 3 → 8). The
polish round below adds one more of each kind, taking `#[test]` to 5,507 and
`test_source_citations.py` to 52. The gate
reports **21** ignored rather than 28 for the reason earlier waves recorded: the
platform-gated tests are not compiled on the Linux gate hosts.

Eight assertion sites were removed in this wave, and every one of them is a
specification change rather than a weakening. They group in three:

- **Four A7 scope tripwires**, which asserted a branch was *unported*: two in
  `tests/file_info.rs` (`class_test_run_consensusxml_is_pending` and
  `class_test_run_fasta_is_pending`, each asserting an
  `Error::Unsupported`) and two in `tests/topp_file_info.rs` (exit code 11 and
  an empty stdout). A7 ported those branches, so each is replaced by an
  assertion of what the branch now does, which is the stronger statement.
- **Three closing TOPP native difference 16**, all in
  `tests/topp_feature_finder_centroided.rs`: `assert_exit(InputFileCorrupt)`,
  `assert_err_contains("nonfinite feature value")` and the no-output-file
  check, replaced by `ExecutionOk`, an empty stderr and assertions on the
  values the Release build itself writes.
- **One lead-ordered**: `SummaryStatistics::new(&mut [NaN])` being refused,
  which the lead's instruction made false, replaced by two new refusal
  assertions for the mixed sample in both orders, by the untouched-input
  assertions around it, and by two whole new tests — `tests/statistic_functions.rs`
  goes from 85 assertion macros to 109.

The count is worth stating carefully, because the obvious way to take it
undercounts: `git diff main..HEAD -- '*.rs' | grep -cE '^-[[:space:]]*assert'`
returns five, missing the one written as an `unwrap_err()` argument and the two
written as `outcome.assert_*` method calls rather than as macros at the start
of a line.

Two hygiene facts, mechanically checked rather than asserted: `Cargo.toml` still
carries `unsafe_code = "forbid"` and `rust-version = "1.85"`, and no C++ entered
the repository — the wave's four oracle directories live under `../oracle/` and
are registered in `SOURCE_PROVENANCE.json` by sha256, all 29 artifact hashes
recomputed in this pass rather than copied from a report.

### Still open, and deliberately not closed here

- ~~**Two shared-math waves now point at the same file.**~~ **Closed** in the
  shared-math wave of 2026-09-19, run as one wave exactly as this pass
  recommended: the x86_64-faithful `variance_with_mean` and a
  libstdc++-faithful `sort_ascending` landed together in
  `src/math/statistic_functions.rs`, under lead decision D16. What A2's oracle
  row still owes — a re-capture against the Linux Release build rather than the
  macOS SDK — is the one part of the promotion bullet left, and
  `src/format/file_info/text_format.rs` was deliberately not touched until it
  arrives. See the shared-math wave section at the top of this document.
- **The whole-document mzML writer still deduplicates by content.**
  `mzml::write` declares one `dataProcessing` and writes no record reference
  where the C++ `MzMLFile::store` declares one and dangles two. Reproducing it
  would make every `mzml::write` caller emit references mzML forbids, by
  default, with nothing to opt into — a change to the crate's default output
  validity rather than a fidelity fix inside one policy, so it is the lead's
  call. Pinned by `the_whole_document_writer_still_deduplicates_by_content`.
- **A chromatogram's `sourceFileRef` is refused, not policed.** mzML 1.1 has no
  such attribute on `ChromatogramType`, so the reader refuses one with
  `Error::Unsupported` before any dangling-reference policy applies, while the
  C++ reads it and default-constructs an empty `SourceFile`. Nothing either
  writer emits produces one, so it is outside the round trip. The divergence
  itself is no longer only asserted in prose:
  `a_chromatograms_source_file_ref_is_refused_under_both_policies` executes it,
  with a reference the header declares and one it does not, under the strict
  default and under `source_dangling_references`.
- **The identification-XML reader's modified-hit budget.** An idXML with more
  than 14 modified peptide hits is refused by the shared reader, because
  `AASequence::parse_with_budget` charges a per-modified-sequence preflight
  against one document-wide 50-million work budget. Two A7 oracle cases
  (`FileFilter_25_input.idXML`, 473 modified hits, and
  `FalseDiscoveryRate_5_input.idXML`, 75) therefore have no differential. The
  oracle records what the C++ prints for both, so a later lane can close it
  without re-running the reference build.
- **The citation checker's item 7** — a bare range under a bare file name — was
  deliberately not attempted. It is the largest unchecked population left (486
  unresolved ranges at the lane's head, 50 here) and the largest false-positive
  surface in the tool; closing it needs its own measurement pass first.
- **`Pins.packaged()` has no test of its own.** The verifier showed that
  mutating it moves two citations on the real tree while the 44-test suite stays
  green. A note, not a blocker: the behaviour is measured, only the guard is
  missing.

### What this checkpoint does not claim

- **It does not claim the benchmark was re-measured.** Nothing ran on
  `ibminode05` in this wave or in this pass, and no cell of `BENCHMARKS.md`
  changed.
- **It does not claim three of the four oracles were re-executed here.** The
  reader, featureXML and picker oracles were run by their own lanes, twice each,
  and this pass registered their drivers by sha256 without re-running them. Only
  `../oracle/a7-fileinfo` was re-executed here, and only to add the three cases
  the A7 verifier's major finding asked for.
- **It does not claim the signed-zero divergence is surveyed.** It is measured
  on one purpose-built pair of fixtures that differ by one exchanged attribute.
  Whether real consensusXML files carry samples with both zeros is not measured,
  and the port's behaviour is deterministic either way.
- **It does not claim a green citation run means the citations are right.** It
  means 104 of 3,346 were read back against quoted code and none of the
  remaining 3,242 is impossible. The section above says what that leaves out.
  Nor does it mean that a wrong line number would always be caught: shifting
  each of the 96 unique confirmed citations by +40 in its own document catches
  92 of them, and the four it misses include two that a span the same unit
  cites around them answers for either way.
- **It does not claim `validated_topp_workflows` moved.** It is 8, unchanged.
  Three of the eight are the tools this wave changed, and all three were already
  validated; no lane added a tool with a `src/bin/*.rs` and a tier-1 manifest.
- **It does not claim the picker consumers moved the ledger.** The handover
  wrote that they "raise `evidence_requires_review` coverage". No count moves:
  `docs/core-sdk-coverage.json` reads complete 63 /
  evidence_requires_review 165 / native_equivalent 90 / partial 59 /
  unmapped 409 at `main` and the same five numbers at this head. What the lane
  did is add `tests/data/picker_consumers_provenance.json` as a second
  reference manifest to `PeakPickerChromatogram.h` and `PeakPickerIterative.h`,
  both of which were already at `evidence_requires_review` and stay there.

## Wave 7: the FileInfo checks, the low-memory picker and the wave-6 benchmark refresh (2026-09-19)

`integrate/wave7` merges `bench/wave6-refresh` (`3096196`), `port/a6-fileinfo`
(`0510382`) and `port/p4-lowmemory` (`4293aab`) onto `main` `36c26a0`. All three
merged without a conflict, and the merged tree changes 119 files against `main`
(+14,398/−663) before this pass's own records. The three lanes are disjoint by
construction: the benchmark lane owns one Markdown file, A6 owns
`src/format/file_info/*` and its tests and fixtures, P4 owns
`src/format/ms_data_writing_consumer.rs`, `src/cli/tools/peak_picker_hi_res.rs`
and the picker's fixtures. See
[FILE_INFO_CHECKS_SUPPORT](FILE_INFO_CHECKS_SUPPORT.md),
[TOPP_PEAK_PICKER_HI_RES_SUPPORT](TOPP_PEAK_PICKER_HI_RES_SUPPORT.md),
[MS_DATA_WRITING_CONSUMER_SUPPORT](MS_DATA_WRITING_CONSUMER_SUPPORT.md),
[BENCHMARKS](BENCHMARKS.md) §3 and
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-7-status).

### The three verdicts, and what this pass had to finish

| lane | round-2 verdict | closing verdict | carried into this pass |
|---|---|---|---|
| `bench/wave6-refresh` | changes_required (1 major, 5 minors) | **approve_with_notes**, 7 minors | 6 applied here, 1 was the integrator's own text |
| `port/a6-fileinfo` | approve_with_notes (3 minors) | **approve_with_notes**, 3 minors | 2 applied here, 1 declined with a reason |
| `port/p4-lowmemory` | changes_required (1 major, 2 minors) | **changes_required**, 2 majors + 5 minors | both majors and 4 minors applied here |

P4 is the one lane that did not close clean, so its two majors were applied in
the merged branch rather than carried. **Neither was taken on the reviewer's
word**: both were re-measured here against the C++ Release build at the pins on
`ibminode06`, and the measurements are registered as
`../oracle/integ-w7/refcheck_06.sh` and `dupdp_06.sh`.

**Major 1 — the dangling reference is not a streaming defect.** The reviewer
claimed the ordinary whole-document `MzMLFile::store` emits the same invalid
references, so P4's own 110-record `FileMerger` evidence file is already invalid
mzML before `PeakPickerHiRes` reads it. Measured: that file declares
`<sourceFileList count="22">` and **exactly one** `<dataProcessing id="dp_sp_0">`,
and carries **106** record `dataProcessingRef`s of which **105 dangle**,
`dp_sp_5` through `dp_sp_109`. A two-part merge is the minimal case — one
declared `dp_sp_0`, six references, five dangling. No streaming consumer is
involved in either. The root cause is confirmed at the pin: `writeHeader_`
deduplicates histories by content with `Helpers::cmpPtrContainer`
(`MzMLHandler.cpp:5060-5069`; `Helpers.h:35-51`) while `writeSpectrum_` compares
them by pointer over `std::vector<std::shared_ptr<const DataProcessing>>`
(`:5258`, `:5265`; `SpectrumSettings.h:165`). The single declared entry against
22 merged parts *is* the content deduplication; the 105 references *are* the
pointer comparison. CPP-172 is rewritten around that, with both triggers and a
fix that covers the non-streaming half.

**Major 2 — `SourceDangling` is not an exact reproduction, and three places said
it was.** The port decides "differs from the first record's" by the text the
history renders; the source decides it by pointer identity. Measured on the
committed five-record `refs` fixture with `dp_sp_1`'s `softwareRef` repointed
from `so_dp_1` to `so_dp_0`, so that `dp_sp_0` and `dp_sp_1` render identically
and only their `id` differs: the C++ low-memory output is **byte-identical** to
its output on the unmodified fixture (sha256 `7c75908440e1c154…`) and still
carries `dataProcessingRef="dp_sp_1"` and `"dp_sp_2"` against a header declaring
only `dp_sp_0`, while this port writes **no** `dataProcessingRef` at all. The
three claims are narrowed to "wherever the source's pointer comparison and
content equality agree", the divergence is recorded as such in native difference
12, and a new consumer test,
`a_history_equal_to_the_headers_by_content_is_not_renumbered`, fails loudly if
the decision ever becomes pointer-like. Reproducing the pointer rule means
carrying the input's own `dataProcessing` identifier through the reader into the
write decision, which is the mzML reader's model to change and is left open.

The same probe measured the figure P4's minor 3 disputed. **Both the lane and
its reviewer were wrong about it.** The lane wrote "three referenceable ids and
109 references"; the reviewer proposed "23 declared and 106 references". The
C++ low-memory *output* declares **two** referenceable ids, `sf_ru_0` and
`dp_sp_0`, and carries **106** references — the reviewer's 23 is the count for
the `FileMerger` *input*, not the tool's output. 106 is right and is also the
only self-consistent figure, since records 1 to 4 carry none and 105 dangle. The
document now says two and 106.

### What else was applied here rather than carried

- **A6, CPP-337's affected class.** The documents said an index loses its first
  offset when it is "written without whitespace". Only a text node *immediately
  after the opening `<index …>` tag* saves it — the walk sets `iter = firstChild`
  and advances before it reads, so only the identity of the first child matters,
  and whitespace between the offsets or before `</index>` does not help. That is
  the mechanism at the pin (`:280-282`, `:290-293`) and it is what the reviewer
  measured on the Release build. Corrected in five places.
- **The benchmark minors, each recomputed from the raw per-repetition records on
  dax rather than accepted.** User time is within 0.45 s in **fourteen** cells,
  not twelve (the six `-fma` tools at both thread counts, FileInfo counting
  twice). The widest 32-thread delta is **0.434 s**: the raw medians are
  26.339561 s and 25.905422 s, and the document's 0.435 was the difference of
  their three-decimal roundings. The worst peak RSS figures are **4,064 MiB** and
  **6,600 MiB**, not "4.06 GB" and "6.60 GB", which were MiB divided by 1000 and
  understated each by 4.9 %. The four last-digit RSS ratio shifts come from
  movements of at most **0.13 %**, which the sentence's own bound already said;
  "sub-tenth-of-a-percent" contradicted it. The `-fma` table's last column is a
  speed-up and is now headed `t1/t32`.
- **Caveat 12's dependency is now disclosed.** §3.1 and caveat 14 say the run's
  one load-flagged execution "enters no table", which is true; caveat 12's
  `-write_ini` figures are prose, and that execution is one of the five samples
  behind DTAExtractor's 3.133 ms median. Recomputed here: the five samples are
  3.133, 3.053, 2.216 (flagged), 3.933 and 3.955 ms, so the median over the four
  unflagged ones is 3.533 ms and the two ranges would read 3.5–4.2 ms and
  60.2–65.7 ms. The published figures keep the run's own filter, which is the
  basis every other median in that document uses, and caveat 12 now prints both.

**Declined, with the reason.** A6's reviewer asked for the CPP-337 correction to
be carried into `../oracle/a6-fileinfo/manifest.json` as well, which means
re-executing the oracle so the manifest is generated and not edited. The five
repository statements are corrected; the oracle manifest's one `known_gaps`
sentence still carries the narrower wording, and its sha256 is registered as it
stands. Re-running 59 Release cases to reword one sentence in a run artifact was
not worth the risk of moving 81 committed expectations in an integration pass.
It is listed as deferred.

### The wave-6 benchmark refresh

`docs/BENCHMARKS.md` gains §3, the wave-6 run of 2026-09-18 on ibminode05
(`2026-09-18-w6refresh`): eight TOPP tools at 1 and 32 threads, the port's
default x86_64 build — which now carries `-C target-feature=+fma` — against the
pinned C++ Release build `openms4-release-bc9cc12-c19e494-174b576`, plus a
`RUSTFLAGS="-C target-feature=-fma"` arm on eight of the nine cases. Output
equivalence under D6 is **unchanged from wave 4 on every tool and every count**,
at both thread counts: **three cases (two tools)** bitwise equal —
DTAExtractor on the Velos mzML and FileInfo on both its datasets — five equal
within tolerance with 0 arrays different, and FeatureFinderCentroided matching
all 4,076 features with 0 unmatched either way and metadata `equal` after the
algorithm was completed. All 52 repetition-determinism and all 26
thread-invariance checks are `bitwise_equal`, on all three implementations. All
sixteen `-fma`-against-`+fma` comparisons are `bitwise_equal` on data; metadata
is `equal` on ten of them and `not_applicable` on the **six comparisons whose
three cases** — DTAExtractor and the two FileInfo datasets — write DTA or plain
text and carry no metadata at all. Details and caveats in
[BENCHMARKS](BENCHMARKS.md) §3.

One tool moved outside the ~3 % cross-session drift band:
FeatureFinderCentroided, whose algorithm was completed between the waves and for
which the build flag is worth 38.8 % at one thread and 8.7 % at 32. The other
seven reproduce wave 4 inside the band. Peak RSS reproduces wave 4 on every tool
but FeatureFinderCentroided, whose Rust peak rose about 4 % — 310.3 → 322.8 MiB
at one thread and 306.4 → 317.8 at 32, against a C++ side that did not move, so
the RSS ratio went 0.86 → 0.89 and 0.81 → 0.84. That wave-4 comparison is the
major its round-2 review raised: the document had claimed the RSS was "unchanged
from wave 4 to the same three digits on every tool", which was false and which
printed no wave-4 number a reader could check. It now prints both waves side by
side. The lane ran nothing on ibminode05 in that round and re-timed nothing; the
correction came from the raw per-repetition records of both waves.

### The ledger

**A6 advances `FORMAT/FileInfo.h`, and it stays `partial`.** The flags now
covered are **`-i`, `-d` and `-c`**, which leaves exactly one flag refused,
`-v`; `partial` is now owed to A7 (consensusXML, identification, FASTA) and A8
(`-v`, mzXML, mzData, trafoXML) alone. `src/format/file_info/checks.rs` and
`tests/file_info_checks.rs` join its `rust` and `tests` lists, each in sorted
position — `checks.rs` sorts *before* `features.rs`, and
`tests/file_info_checks.rs` is *second* after `tests/file_info.rs` because `.`
precedes `_`; the lane's own placement hints had both the wrong way round.

**P4 advances two headers.** `FORMAT/DATAACCESS/MSDataWritingConsumer.h` records
the indexed footer (no longer an exception), `ReferencePolicy`, the `softwareList`
gap that the refs fixture exposed in `Checked`, and the content-versus-pointer
divergence above. `PROCESSING/CENTROIDING/PeakPickerHiRes.h` records
`-processOption lowmemory` as ported and drops it from the PARTIAL list, which
now holds the Mobilogram overloads, `pickExperiment` on `OnDiscMSExperiment` and
the `ProgressLogger` base.

**`validated_topp_workflows` does not move. It stays at 8.** That is the honest
answer, not a missing promotion. The counter is derived, not written: it counts
tools whose TOPP *package* provenance manifest declares tier 1 and names an
upstream test definition, and both `FileInfo` and `PeakPickerHiRes` were already
among the eight. This wave deepens those two workflows — three more FileInfo
flags, a second PeakPickerHiRes process option — without validating a ninth
tool. A6's new manifest, `tests/data/file_info_checks_provenance.json`, is a
**core SDK** reference manifest rather than a TOPP package one (it cites
`src/openms/` sources), so by construction it cannot move this counter either.
`core_sdk_coverage.py --write` confirms every count unchanged: 786 registered
public headers, 63 complete / 165 evidence_requires_review / 90
native_equivalent / 59 partial / 409 unmapped, 124 TOPP sources, 8 validated.

### C++ issues

Six new entries, `CPP-335` to `CPP-340`, and one rewrite. **The numbers are not
the ones the lanes proposed:** all three lanes independently claimed "CPP-335",
so they were assigned here after main's highest, `CPP-334`. A6's three keep
335–337, because `CPP-337` was the only number a lane had already written into
its own repository files; the benchmark crash is 338, and P4's two are 339 and
340. Every citation was re-read at the pins before the entries were written.

The dangling-reference finding is **not** a new number: `CPP-172` already covered
it as a source-review entry, so it is replaced in place, promoted to Executed,
and widened from "streaming consumer" to the `writeHeader_`/`writeSpectrum_`
asymmetry with its two triggers. All 340 entries are `##` headings; the lanes'
requests used `###`, which would have nested them one level too deep.

### Lead decisions of this wave, and where they landed

1. **Reproduce the source's dangling references rather than refuse them**
   (native difference 12). Discharged, and now qualified by the measured limit of
   that reproduction. The 110-record `FileMerger` file is measured end to end at
   the lane's final head and again here; it is pinned in the oracle rather than
   by a repository test, because the merged input is ~9.3 MB and not byte-stable
   (`FileMerger` runs without `-test` and stamps a time in). The committed
   five-record `refs` fixture pins all four cells of the rule, where the
   `FileMerger` file exercises one.
2. **`-i` answers with this port's index decoder** (native difference 9 of the
   FileInfo tool document). Discharged: reachable from both `-i` rows of the
   capability table and from the `ValidationInfo` row of the library document,
   with both boundaries, their pinned lines, their owner
   (`src/format/indexed_mzml.rs`) and why the departure stands.
3. **File CPP-337.** Discharged, as `CPP-337`, with three kinds of evidence.

### Gates at the integration head

Node discipline held: **nothing ran on ibminode05 at any point in this pass.**
The full suite ran on kim, every other gate on dax, and the two C++ probes on
ibminode06. Logs under the wave-7 integration log directory.

| gate | host | result |
|---|---|---|
| `fmt --all -- --check` | local | rc 0 |
| `clippy --locked --all-features --all-targets -- -D warnings` | dax | rc 0, zero warning lines |
| `+1.85.0 check --locked --all-features --all-targets` | dax | rc 0 |
| `doc --locked --all-features --no-deps`, `RUSTDOCFLAGS=-D warnings` | dax | rc 0, zero warnings |
| `test --locked --all-features --all-targets --no-fail-fast` | kim | **5,401 passed / 0 failed / 21 ignored** |
| `test --locked --no-default-features --all-targets --no-fail-fast` | dax | 3,565 passed / 0 failed / 3 ignored |
| `--no-default-features --features mzml,paramxml,featurexml --test topp_file_info --test topp_feature_finder_centroided --test file_info_checks` | dax | 51 + 45 + 34, rc 0 — CI job `test` |
| the same slice on **`+1.85.0`** with the five other targets of CI job `minimum-rust` | dax | 19 + 39 + 35 + 51 + 45 + 34, rc 0 |
| `--no-default-features --features mzml --test file_info_checks` | dax | 47 / 0 / 0 |
| `--no-default-features --features mzml,paramxml --test topp_peak_picker_hi_res --test ms_data_writing_consumer` | dax | 33 + 36, rc 0 |
| `--no-default-features --features mzml --test ms_data_writing_consumer --test mzxml --test mzdata` | dax | 36 + 64 + 61, rc 0 — CI job `minimum-rust` |
| `tools/check_core_sdk.py`, plain and `--source <absolute pin path>` | local | rc 0 both; 2,092 distinct source/registration/reference files verified at `bc9cc12` |
| `tools/check_doc_coverage.py --write` | local | 4,560/5,915 = 77.1 % |
| `tools/check_module_cycles.py` | local | 64 cross-module edges, 13 mutually-dependent pairs — unchanged |
| `tools/check_schema_feature_graph.py`, `tools/core_sdk_coverage.py`, `tools/test_core_sdk.py`, `tools/test_core_sdk_coverage.py` | local | rc 0 |
| YAML parse of `.github/workflows/rust.yml`; `json.load` of all six changed JSON files | local | all parse |

**The difference from `main` is fully accounted for.** `main` is
5,317 / 0 / 21; this head is 5,401 / 0 / 21, which is **+84**: A6 adds 63
(5,380 at `0510382`, the figure both the lane and its reviewer reported), P4
adds 20 (5,337 at `4293aab`, likewise), and this pass adds **one** —
`a_history_equal_to_the_headers_by_content_is_not_renumbered`, the test that
pins P4's second major. 63 + 20 + 1 = 84. The ignored count does not move.

The suite was summed by reconciliation rather than by trusting one number,
because the ssh capture on the gate path drops lines in two different ways: the
log is walked, each `test result:` line is compared with the `... ok` and
`... ignored` lines counted since the previous one, and the larger of each pair
is taken. The `--all-features` capture came through **clean** — 356 result lines
and **zero** disagreements — so its 5,401 is the same figure on either method.
The `--no-default-features` capture had exactly one damaged spot (a result line
reading 7 where 17 had been counted, an eaten `test result:` line whose window
then covered two binaries); the reconciliation repaired it, and no other line
disagreed. Neither log contains a `FAILED` marker, a `failures:` block or a
panic. The `--all-targets` form excludes doctests, which is why its
`--no-default-features` total is below the 3,626 A6 reported for the plain form.

### Ignored tests

**28 `#[ignore]` attributes in the tree, all of them in `tests/` and none in
`src/` — byte-identical to `main`**, file for file and count for count, across
the same ten files. No test was ignored, skipped, deleted or weakened by this
wave; no tolerance was widened; no expected value was derived from Rust output.
The `.rs` diff against `main` adds **0** occurrences of `unsafe` and **0** new
`#[ignore]`, and adds 85 `#[test]` functions (5,345 → 5,430). The gate reports
**21** ignored rather than 28 for the reasons earlier waves recorded: the
platform-gated tests are not compiled on the Linux gate hosts.

Two hygiene facts, mechanically checked rather than asserted: `Cargo.toml` still
carries `unsafe_code = "forbid"` and `rust-version = "1.85"`, and no C++ entered
the repository — the two integration probes are shell drivers under
`../oracle/integ-w7/`, registered in `SOURCE_PROVENANCE.json` by sha256, with no
copy inside the crate.

### Still open, and deliberately not closed here

- **The reader's other half.** This port's reader refuses an unregistered
  spectrum `sourceFileRef` under *either* dangling-reference policy
  (`src/format/mzml_header/read.rs:113-121`) where the source warns once and
  continues (`MzMLHandler.cpp:899-906`). Now that the writer reproduces the
  source's references, the port writes a low-memory output it will not read back
  on an input with per-record source files, and will not read the C++ output of
  the same run either, while the C++ reader reads both. The `FileMerger` case has
  only dangling `dataProcessingRef`s and round-trips on both sides. Changing it
  reverses an earlier lane's documented decision and belongs with the mzML
  reader, not this tool. **Closed in wave 8 as decision D14**
  (`fix/reader-round-trip`): the reader accepts what the Release build accepts
  under `ReadOptions::source_dangling_references`, and the refusal stays as the
  default strict profile.
- **Reproducing the pointer rule** in `SourceDangling`, above. **Closed in
  wave 8** for the streaming consumer, without the reader change this bullet
  anticipated: the identity was already in the model. What remains of it is in
  the whole-document writer and is listed in wave 8's own open items.
- **Two classes of corruption `-c` cannot report**, because this port refuses
  them before `-c` sees them: a repeated auxiliary array name
  (`src/format/mzml.rs:1038`) and an MS-level-0 mass spectrum
  (`src/kernel.rs:810-821`). The C++ loads both and lets `-c` do its job.
- **CPP-337's singleton consequence is a probe, not a pinned oracle case.**
  `build(250, 1, "")` makes the Release `FileInfo` print "0 spectra" and exit 0 —
  an entire index section vanishing with a success status. Pinning it costs a
  fixture, two expectations, a test and five recounts.
- **SpectraFilterWindowMower is unresolved at n = 3** in the benchmark, and has
  no `-fma` arm; the C++ FeatureFinderCentroided crash is one event, not a
  diagnosis. Both are labelled open in `BENCHMARKS.md` rather than resolved.

### What this checkpoint does not claim

- **It does not claim the benchmark was re-measured.** Nothing ran on
  ibminode05 in this wave's closing round or in this pass. No cell, table or
  measurement in `BENCHMARKS.md` changed; the raw per-repetition records were
  read from dax, and the figures this pass corrected are prose restatements
  recomputed from those same records.
- **It does not claim C++ evidence for the `dupdp` case beyond one fixture.**
  The content-versus-pointer split is measured on one purpose-built five-record
  file and on the 110-record `FileMerger` output. It is not a survey of how often
  real files carry content-equal, pointer-distinct histories.
- **It does not claim the oracle manifests were regenerated.** A6's oracle was
  re-executed by its own lane, twice, at its canonical path; this pass did not
  re-run it and did not edit it. The integration probes are new drivers with
  their own logs, not modifications of an existing oracle.
- **It does not claim anything about a non-FMA processor.** Every gate ran on
  kim, dax or ibminode06, all of which have FMA.

## Wave-6 FAIMS closure and the FMA build default (2026-09-18)

`integrate/wave6` merges `port/b11-faims` (`7921409`) and `port/fma-default`
(`86b0337`) onto `main` `e1c3115`. The two branches share **no file** — 13 files
against 10, `comm -12` of their name-only diffs empty, re-checked at those exact
heads — so both merges were conflict-free. That was verified rather than
assumed, and then confirmed after the fact: the merged tree changes exactly 23
files against `main`, which is 13 + 10. Every integrator-owned record was left
to this pass. See
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-6-status),
[TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT](TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md)
and [FMA_BUILD_FLAG](FMA_BUILD_FLAG.md).

Both lanes were adversarially verified and both came back
**approve_with_notes**: B11 in round 2 with four minors, `port/fma-default` with
three. Every minor was applied on its own branch before this pass, and each lane
reported one deviation from a verifier's *proposed fix*, neither of which
changed the finding:

- B11 declined the wording "a convenience wrapper of `mergeOverlappingFeatures`"
  for `mergeFAIMSFeatures`, because at the pin it does not delegate — it builds
  its own callback and calls `filter` itself (`FeatureOverlapFilter.cpp:384-527`,
  `:507-511`) — and took the verifier's alternative instead, quoting
  `mergeFAIMSFeatures`' own Doxygen block with per-line pins.
- `port/fma-default` could not route `atan_is_reference()` through
  `cpu_features::cpu_provides_fma()` as its verifier proposed, because
  `analysis -> system` closes a module cycle and
  `tools/check_module_cycles.py` — a CI step — refuses it. It reads the same
  architectural bit (leaf 1, `ECX` bit 12) through `raw_cpuid` directly, with the
  production copy named in the doc comment. The verifier's substance is met:
  `cpuid` is executed, so the answer is the processor's and not the build's.

### The ledger

**No header changes status this wave**, and that is the correct outcome rather
than a missing promotion. `PROCESSING/FEATURE/FeatureOverlapFilter.h` was
already `complete` and stays `complete` with its scope rewritten: the sentence
that deferred the FeatureFinderCentroided FAIMS closure to B11 is replaced by
what B11 actually added, one native mode (`FaimsMergeFidelity`) beside faithful
entry points that are unchanged. `IONMOBILITY/IMDataConverter.h` stays
`partial`, because B11 consumed `splitByFAIMSCV` and ported no further member;
its scope says so. Review state across the whole ledger is unchanged: complete
**63**, partial **59**, `native_equivalent` **90**,
`evidence_requires_review` 165, unmapped 409, 786 registered public headers.

**`validated_topp_workflows` stays 8, and that is a decision, not an oversight.**
The count is derived, not written: `tools/core_sdk_coverage.py`'s
`validated_workflows()` reads `SOURCE_PROVENANCE.json`'s
`topp_package_reference_manifests` and counts a tool whose manifest declares
tier 1 *and* names the upstream test definition it reproduces.
`FeatureFinderCentroided` was already in that set on `main` — checked, the eight
names are identical before and after — on the strength of its non-FAIMS
differential evidence, so the closure cannot raise the count. It must not lower
it either. The FAIMS evidence B11 adds is itself an executed differential
against retained C++ output: each compensation-voltage group written as its own
single-voltage mzML with the FAIMS cvParam removed, run through the C++ Release
build on `ibminode06`, and the port's features for that group compared with what
that run found. That strengthens the tier-1 claim. The **cross-voltage merge**
is the part with no C++ oracle, and the tool's tier-1 standing never rested on
it; the manifest's own `method` field says so in as many words.

### The FAIMS closure, and where its evidence stops

Decision **D5 is closed for `FeatureFinderCentroided` (2026-09-18)**. The tool
no longer refuses FAIMS input: it splits by compensation voltage, runs the
picked feature finder once per voltage on that voltage's seeds, annotates every
feature with its `FAIMS_CV` and merges across voltages under
`-faims_merge_features`. `IMDataConverter` stays `partial` for the members D5
left out.

Three C++ defects lie on that path, all already recorded: `CPP-278` (the split's
groups carry no ranges, so the C++ tool exits 8 on **every** FAIMS input),
`CPP-282` (the merge erases every feature whose unique id is still 0) and
`CPP-283` (a cluster of three or more voltages is split in two, double-counting
one member). The port answers all three. What this wave does **not** claim:

- **No patched C++ build was ever executed.** The composite statement appended
  to `CPP-278` — that fixing only the ranges makes the tool write an empty
  feature map — is read from the pinned source and confirmed on the port's
  call-for-call faithful merge. Nobody has seen a C++ run produce that empty
  map, and the outcome additionally rests on `CPP-282`'s own premise. The entry
  says this in its own "Basis of the composite claim" paragraph.
- **The statement is qualified, not general.** The empty map follows only on an
  input where at least one cross-voltage merge actually fires. With
  `-faims_merge_features false` the merge is never called, and on single-voltage
  FAIMS input the callback refuses every pair, so nothing is erased. Both
  non-merging cases are pinned by executed tests of this package.
- **The corrected merge has no C++ oracle and cannot get one from the pinned
  build.** It is pinned against the specification derived from the source's own
  parameter documentation and against hand-derived cases whose numbers are
  written out in `tests/feature_overlap_filter.rs`, each written next to the
  executed `c2_*` case that records what the source does instead.

One change outside the FAIMS path: the tool's `OPENMS_LOG_WARN` lines now go to
**stderr**, where the executed C++ writes them. No test asserted the old
destination.

### The processor guard: what it rests on, and what it does not cover

`.cargo/config.toml` sets `-C target-feature=+fma` for
`cfg(target_arch = "x86_64")`. This is a **breaking runtime change**: a binary
built from this checkout needs an FMA3-capable processor, Intel Haswell (2013)
or AMD Piledriver (2012) and newer.

What is measured, in [FMA_BUILD_FLAG](FMA_BUILD_FLAG.md) §7: the flag can change
only instruction encoding, element-wise vector width and `mul_add`'s
implementation. It cannot contract `a * b + c`, because Rust lowers that to
`fmul`/`fadd` with no fast-math flags, and it cannot auto-vectorise a
floating-point reduction, because that needs `reassoc`, which Rust never sets.
Both matter here and are not theoretical: `levenberg_marquardt.rs`'s
`lane_madd` reproduces Eigen's `pmadd` **unfused**, and `eigen_sum` hand-writes
Eigen's two-lane summation order. Five `rustc --emit=asm` probes back this, and
the probe source and both `rustc` commands are printed in the document, so the
check depends on nothing outside the file. Cross-check on the shipped binaries:
the crate has exactly 41 `f64::mul_add` call sites and the `+fma`
`FeatureFinderCentroided` on kim has exactly 41 `vfmadd`/`vfmsub`, so no fused
multiply-add in it came from anywhere but an explicit `mul_add`.

The guard's documented limits, each accepted by the lead rather than left
implicit:

- **Test binaries are not guarded.** They do not go through `cli::run`, so on a
  processor without FMA a test binary dies on `SIGILL` with no message, where a
  tool binary prints the requirement and exits 12. This is also the CI failure
  mode if a runner ever lacks FMA.
- **The flag stays scoped to x86_64.** A 32-bit x86 build does not get it, so
  the guard is inert there, which is consistent.
- **Pre-AVX processors get best effort**, measured rather than guaranteed. A
  separate baseline-built launcher is out of scope.
- **No runner's processor was measured at this integration.** Nothing in this
  pass executed a CPUID read, a `/proc/cpuinfo` read or a job on a GitHub
  runner. That the CI runners satisfy the requirement is an inference from the
  published images being far newer than 2013, not a measurement. The first CI
  run of this branch is the measurement.

### Standing hazard: `is_x86_feature_detected!` is a compile-time constant here

**From `port/fma-default` on, `is_x86_feature_detected!` answers a question
about the build, not about the processor, on x86_64.** `.cargo/config.toml`
builds x86_64 with `-C target-feature=+fma`, and the macro is documented to
answer `true` *without consulting the processor* for any feature the build
already enables: on `x86_64-unknown-linux-gnu` that is `fma`, `avx`, `sse3`,
`sse4.1`, `sse4.2` and `ssse3`, and on `x86_64-apple-darwin`, whose baseline
already has most of those, `fma`, `avx` and `sse4.2`. Under `-O` the guard
written around it is then deleted outright. **Nothing warns**: no compiler
diagnostic, no clippy lint, and in review the code reads exactly right while
being absent from the binary.

This is not hypothetical. It cost `port/fma-default` a shipped guard that was
not there — `FileInfo::main` disassembled to a bare `jmp` into the body of `run`
— and it had already, silently, turned the `atan_is_reference()` gate in
`src/analysis/feature_finder_picked/glibc_libm.rs` into a constant `true`, which
that lane fixes. The same short-circuit is built into the `cpufeatures` crate's
`new!` macro, which the lock already carries through `sha1`; it is unaffected
today only because `sha1` asks for `"sha"`, which `+fma` does not enable.

The rule, recorded here and as a bullet in
`.claude/skills/openms-port-header.md`: a question about the **processor** goes
to `system::cpu_features::cpu_provides_fma()`, or — from a module that may not
name `crate::system` without closing a module cycle — to `raw_cpuid` directly
for the same architectural bit. `is_x86_feature_detected!` and
`cfg!(target_feature = ...)` answer a question about the **build**, which on
x86_64 now has a known answer. Measurements, the affected sites and the
per-target feature table are in [FMA_BUILD_FLAG](FMA_BUILD_FLAG.md) section 6.

### Gates

All on **kim**, slot `integ-w6`. Every gate rc 0 on its first attempt —
**no gate exited 255, so none needed a rerun**. Driver and logs are in the
session scratchpad under `integ-w6-logs/`. Every figure below is summed from the
log's `test result:` lines, over **all** of them, so a dropped target cannot
hide.

The eight-gate battery ran at `288dfb6`, the last commit before this checkpoint
was written. Because the checkpoint itself changes two files, the six decisive
gates were then **re-run at `4c5c6ec`, the branch's final code-bearing head**,
rather than arguing that documentation cannot affect a build: MSRV check,
clippy, rustdoc, the full suite, the `--no-default-features` suite and the
doctests all exit 0 again with **every count identical** — 5317 / 0 / 21 over
355 lines, 3616 / 0 / 3 over 331 lines, 76 doctests. The only change after that
run is the paragraph you are reading.

| Gate | Result |
|---|---|
| `+1.85.0 check --locked --all-features --all-targets` | exit 0 (MSRV 1.85), non-vacuous: the log compiles `openms` itself |
| `clippy --locked --all-features --all-targets -- -D warnings` | exit 0, **0 warnings**, non-vacuous |
| `doc --locked --all-features --no-deps`, `RUSTDOCFLAGS=-D warnings` | exit 0 |
| `test --locked --all-features --all-targets` | **5317 passed, 0 failed, 21 ignored** over **355** result lines |
| `test --locked --no-default-features` | **3616 passed, 0 failed, 3 ignored** over **331** result lines |
| `test --locked --all-features --doc` | **76 passed, 0 failed** (73 + 3) |
| `+1.85.0` minimum-rust line :109 as changed (`--no-default-features`, + `fma_build_flag` + `build_info`) | 177 passed, 0 failed, 2 ignored over 10 result lines |
| `--all-features`, the two lanes' targets (`topp_feature_finder_centroided`, `feature_overlap_filter`, `fma_build_flag`, `build_info`, `topp_cli_lifecycle`) | 170 passed, 0 failed, 0 ignored: **45, 39, 3, 10, 73** |

Each log is internally consistent: for the full suite, 355 `Running` headers =
355 `running N tests` lines = 355 `test result:` lines, every status `ok`, and
the announced total 5338 equals 5317 + 0 + 21, so nothing was lost through the
pipe. The same three counts agree for every other gate.

An earlier single-gate probe, run before the battery and before any document was
written, had already proved the changed CI line at `+1.85.0`: rc 0.

The five per-target counts in the last row reproduce **exactly** what the B11
lane reported for them on its own branch (45, 39, 10, 73) with the FMA lane's
new target beside them, which is the check that the merge changed neither lane's
behaviour.

**Repair round (2026-09-18).** The audit of this checkpoint found that the
"does not claim" bullet above recorded an evidence gap that did not exist, and
that the gap had hidden wrong line citations. Correcting them changed
documentation, one JSON string and doc comments in
`src/cli/tools/feature_finder_centroided.rs` — nothing executable: with comment
lines stripped, that file is byte-identical to its parent at both repair
commits. The six decisive gates were re-run anyway, at `422233c`, the last
commit on this branch that changes anything the compiler reads; every commit
after it, including the one that adds this paragraph, is documentation only.
They ran with the worktree clean and equal to that commit for the whole run
(2112 tracked files, every blob hash recomputed, nothing untracked): MSRV check,
clippy, rustdoc, the doctests, the full suite and the `--no-default-features`
suite, **all exit 0 on the first attempt, none exited 255**, with every count
identical to the `4c5c6ec` run — **5317 / 0 / 21** over **355** result lines,
**3616 / 0 / 3** over **331**, and **76 doctests**. The full suite's target list
was compared run against run as well: 328 targets each, symmetric difference
empty, so no target silently dropped out.

One honest wrinkle in that log, because a future reader diffing it will hit it.
`x_test_all.log` has **four lines where the ssh transport interleaved a
`Running` header into the middle of another line**, and two of them swallow
numbers: at `:848` a `test result:` line lost its counts, and at `:7250` a test
line lost its trailing `ok`. A naive `awk` over the log therefore reports 5308
rather than 5317. Both reconcile exactly and independently: the 354 well-formed
result lines sum to 5308 passed + 21 ignored = 5329, and the announced total is
5338, leaving exactly the 9 tests of the one target whose result line was
mangled; and the unanchored `... ok` count is 5316, one short of 5317, which is
the line at `:7250`. `cargo test` also exits non-zero on any failure and this
gate exited 0, and the log contains no `FAILED` marker and no `failures:` block.
The `4c5c6ec` log has no interleaved line at all, which is why its `awk` total
was clean. **No test changed status; the difference is in the pipe, not the
suite.**

### The suite totals against `main`, target by target

`main`'s recorded figures were not taken on trust. `main` (`e1c3115`) was
checked out into its own worktree and run on the **same host, same toolchain,
same commands**, and it reproduces the wave-5 record exactly: **5297 passed, 0
failed, 21 ignored over 354 result lines** with all features, and **3599
passed, 0 failed, 3 ignored over 330 result lines** with none. The two runs were
then compared target by target, not just in total.

| | `main` `e1c3115` | `integrate/wave6` | delta |
|---|---|---|---|
| `--all-features --all-targets` | 5297 / 0 / 21 over 354 lines | **5317 / 0 / 21** over **355** lines | **+20 passed**, +1 line, ignored unchanged |
| `--no-default-features` | 3599 / 0 / 3 over 330 lines | **3616 / 0 / 3** over **331** lines | **+17 passed**, +1 line, ignored unchanged |

**Every one of those tests is accounted for, and no other target moved by a
single test.** With all features the +20 is exactly four targets:

| Target | `main` | wave 6 | delta | Lane |
|---|---:|---:|---:|---|
| `tests/topp_feature_finder_centroided.rs` | 39 | 45 | +6 | B11 |
| `tests/feature_overlap_filter.rs` | 33 | 39 | +6 | B11 |
| `src/lib.rs` unit tests | 358 | 363 | +5 | FMA (`cpu_features`) |
| `tests/fma_build_flag.rs` | — | 3 | +3 | FMA (new target) |

6 + 6 = **12** for `port/b11-faims` and 5 + 3 = **8** for `port/fma-default`,
which are exactly the deltas the two lanes reported on their own branches
(5309 and 5305 against the same 5297). The merge is additive to the test, which
is what two branches sharing no file must produce. The one extra result line is
`tests/fma_build_flag.rs`, the wave's only new test target.

Without default features the +17 is the same three feature-independent targets —
`src/lib.rs` +5, `feature_overlap_filter` +6, `fma_build_flag` +3 — plus **+3
doctests** (58 to 61), the `cpu_features` documentation examples;
`topp_feature_finder_centroided` is gated behind `mzml`, `paramxml` and
`featurexml` and contributes nothing there. The dedicated doctest gate shows the
same three: 73 + 3 = **76**, against `main`'s 73.

Locally on macOS arm64, all rc 0: `cargo fmt --all -- --check`; all seven
repository Python checkers — `check_core_sdk` (plain **and** with `--source`
against the pinned checkout, 2,092 distinct source/registration/reference files
verified at `bc9cc12`), `check_doc_coverage`, `check_module_cycles`,
`core_sdk_coverage`, `test_core_sdk`, `test_core_sdk_coverage` and
`check_schema_feature_graph`; and the ten generator `--check` scripts the
`quality` job runs. `.github/workflows/rust.yml` parses as YAML (6 jobs, 111
steps) and every changed JSON record round-trips through `json.load` — as does
every tracked `.json` in the repository, checked in passing.

Both ledger generators are idempotent: re-running
`core_sdk_coverage.py --write` and `check_doc_coverage.py --write` after the
edits leaves `docs/` clean. `docs/module-cycles.json` is **byte-unchanged**
against `main`, so neither lane added a cross-module edge.

The integrator pass touched **no** `.rs` file, no file under `tests/data/`, and
neither `Cargo.toml` nor `Cargo.lock`: its 15 files are the ledger, the
provenance, the C++ issue log, the crate and licence records, CI, the porting
skill and six documents. `Cargo.toml`'s `[lints.rust] unsafe_code = "forbid"` is
untouched, `src/` contains no `unsafe`, and **0** C++ files are tracked.

### CI audit

`.github/workflows/rust.yml` gains exactly two `--test` names, both on the same
`minimum-rust` line — the general `--no-default-features` slice that already
carried `feature_overlap_filter` — and one new step. `fma_build_flag` is this
wave's only new test target; before the change it ran at `1.85.0` only under
`--all-features --all-targets`, because `minimum-rust` has no bare
`cargo test --no-default-features` line the way the `test` job does.
`build_info` joins it for the same reason and because this wave changes what it
reports. The new `quality` step,
`RUSTFLAGS="-C target-feature=-fma" cargo check --locked --all-features
--all-targets`, keeps the opt-out this crate prints from silently rotting.

**No test binary in `tests/` is unrun.** There are **328** integration targets
on disk, one more than the 327 of wave 5, which is `tests/fma_build_flag.rs`.
Every `--test` name in the workflow resolves to one of them (**0 dangling**, 178
distinct names), and the `test --locked --all-features --all-targets` gate
below launched all 328. The other **150** are reached only by the
`--all-features --all-targets` steps, which is sufficient because all **80**
file-level `#![cfg(...)]` gates in `tests/` are cargo feature gates on features
declared in `Cargo.toml` — checked mechanically, **0** of them uses a
non-feature predicate — so `--all-features` satisfies every one.

### Ignored tests

**22 `#[ignore]` attributes in the tree, all of them in `tests/` and none in
`src/` — byte-identical to `main`**, file for file and count for count. No test
was ignored, skipped, deleted or weakened by this wave; no tolerance was
widened; no expected value was derived from Rust output. Across both branches
the `.rs` diff against `main` adds **0** occurrences of `unsafe` and **0** new
`#[ignore]`, and adds 20 `#[test]` functions. The gate below reports **21**
ignored rather than 22 for the reason wave 5 recorded: the macOS-arm64-gated
`macos_arm64_sdk_gap_report` in `tests/lm_eigen_path_differential.rs` is not
compiled on `kim`.

### What this checkpoint does not claim

- **It does not claim any C++ evidence for the cross-voltage merge.** See *The
  FAIMS closure, and where its evidence stops* above. No patched C++ build was
  executed, and the pinned one exits 8 before the merge on every FAIMS input.
- **It does not claim that the FMA flag was validated on a processor without
  FMA.** Every gate ran on `kim`, which has FMA, so the guard's refusal path was
  exercised by its unit tests and not by a real refusal on real hardware. That
  the gate counts are identical to each lane's pre-integration counts is the
  expected outcome of running on an FMA host, not evidence about a non-FMA one.
- **It does not claim anything about the CI runners' processors.** Nothing here
  executed a job on a GitHub runner.
- **It does not claim a performance result.** No benchmark was run at this
  integration. §4 of [BENCHMARKS](BENCHMARKS.md) is the wave-5 measurement on
  **dax**, and §4.5 now records the decision rather than new numbers.
- **It does not re-verify the pinned sources line by line — but every line
  reference this wave adds was checked against its pin.** The `CPP-278` core
  references were checked against the pinned checkout
  `.reference/openms4-core-bc9cc12`: `FeatureOverlapFilter.cpp:384-527` is
  `mergeFAIMSFeatures` and calls `filter` itself at `:507-511`, `:263-271`
  inserts into `removed_uids` only on a `true` callback, `:277-281` erases
  exactly those ids, `:445-448` refuses a same-voltage pair, and the function's
  own Doxygen block is `:156-180`, with the sentences quoted at `:157`,
  `:159-160` and `:175`. The `FeatureFinderCentroided.cpp` references were
  checked against the pinned `topp` revision `174b576`, which **is** locally
  reachable: `OpenMS4-tests/packages/topp` is the `OpenMS4-topp` repository, the
  pin is an ancestor of its `HEAD`, and `git show
  174b576:src/FeatureFinderCentroided.cpp` yields the pinned file (387 lines).
  Correct as written: `:258-281`, `:283-291`, `:294-299`, `:309`. **Five were
  wrong and are corrected here**: in
  `src/cli/tools/feature_finder_centroided.rs` the group log line is `:254` not
  `:253`, the combined-features line `:306` not `:307`, the merge log line
  `:313-314` not `:314-315`, and the merge's literal `5.0`/`0.05` arguments
  `:312` not `:313`; in `CPP-278` the unique ids are assigned at `:328-329`, not
  `:318-320`. Every correction is to a line number only — each claim's subject
  was right, including `CPP-278`'s, whose ids are indeed assigned after the
  merge call at `:312`. The wrong numbers match no local copy at any revision:
  the `topp-sdk-validation/source` copy is byte-identical to the pin, and the
  `topp` working tree at `HEAD`, which differs from the pin elsewhere, carries
  the same six FAIMS anchors at `254`, `306`, `312`, `313`, `318` and `328`.
  Checking those references opened the CLI pin `c19e494` as well, so every
  distinct `TOPPBase.cpp` citation in this repository was swept against it: 66
  when the sweep ran, and 67 in 100 occurrences across 26 files at this head,
  because the sentence below quotes a citation that no other file carried. One is wrong, and it is **not** this wave's:
  `TOPPBase.cpp:519-522`, which `main` already carries in three places, is blank
  lines and the head of `toolName_()` at the pin. The outer
  `catch (const std::exception&)` those three sentences mean — the one that
  reports `Unable to initialize or run ...` and returns `INTERNAL_ERROR` — is
  `:510-514`, and it is corrected here; the sibling citation `:505-508` for the
  `BaseException` handler was already right. No other `TOPPBase.cpp` citation
  mismatched.
  What is still unclaimed is the rest of those files: only the references this
  repository cites were read, and reading them is a source check, not an
  executed build.

## Wave-5 completion of FeatureFinderAlgorithmPicked and the noise estimators (2026-09-17)

`integrate/wave5` merges `port/ffap-complete` (`a11fc26`, itself the merge of
`port/ffap-instrumentation` `c79e66f`, `port/ffap-semantics` `511d29e` and
`port/progress-logger-release-range` `f89d5d4`, plus six combined fix rounds and a
minors pass) and `port/signal-to-noise` (`be70a98`) onto `main` `59e0e1c`. The two
branches share **no file**, so both merges were conflict-free and every
integrator-owned record was left to this pass. See
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-5-status) and
[BENCHMARKS](BENCHMARKS.md) §4.

This wave promotes three headers and closes the port's last documented divergences on
the feature-finder path:

| Header | Before | After |
|---|---|---|
| `FEATUREFINDER/FeatureFinderAlgorithmPicked.h` | `partial` | **`complete`** |
| `PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h` | `partial` | **`complete`** |
| `PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h` | `partial` | **`complete`** |
| `CONCEPT/ProgressLogger.h` | `native_equivalent` | `native_equivalent` (scope rewritten) |
| `FEATUREFINDER/EGHTraceFitter.h`, `GaussTraceFitter.h` | `complete` | `complete` (scope rewritten) |
| `CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h` | `partial` | `partial` (scope rewritten) |
| `FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h` | `partial` | `partial` (scope rewritten) |

Review state across the whole ledger: complete 60 -> **63**, partial 62 -> **59**,
`native_equivalent` 90 unchanged, 786 registered public headers unchanged.

**No status overstates its evidence.** Each promoted scope sentence states, in the
lane's own words, what is reproduced, what is refused and why, and the one accepted
exception. For `FeatureFinderAlgorithmPicked.h` the refusals are exactly lead decision
D1's classes — an out-of-bounds read or write, a loop that never ends, process
termination, or output that depends on a heap address — every one of them recording
where and how the executed process ends, and the single exception is the multi-thread
race on `aborts_`/`abort_reasons_`/`log_`, which gives the single-thread result under
**lead decision D11** because the determinism contract requires parallel output to equal
serial output.

### Lane evidence

- ffap-instrumentation (port/ffap-instrumentation, 8a53517/c79e66f): the FeatureFinderAlgorithmPicked instance and its side channels are ported and checked against the Linux x86_64 Release build. The oracle is ../oracle/ffap-instr-completion: drivers ffap_instr_driver (source 894fbfa9, binary d20242a3), ffap_progress_driver (bf54f670/a742aab8) and ffap_shift_band_driver (74489272/c8f85f74), plus tool_cases.py, run on ibminode06 with OMP_NUM_THREADS=1. Checked: reuse of one object over three runs, and runs into caller maps (11 prefilled and 5 overlapping features, NaN and infinite keys); the parameter surface, including a refused set; the full ProgressLogger event sequence of five cases under a counting time(); write_debug for tool cases a1, a2 and a3 (log.txt byte-identical, featureXML and mzML under D6), c3 (a1 at four threads, defined, equal to a1), and b1 (C++ SIGABRT; the port exits 8 with the same pre-termination files); 375 writeFeatureDebugInfo_ files byte-identical, plus 75 each for a string shift, a list shift, and a string shift with a scan at RT 5e-275 and at 1e-289 (3 processes, 3 addresses, identical). Refused where the source terminates or is not reproducible, and measured there: SIGABRT; SIGSEGV 7 of 7; SIGFPE 2 of 2; an address-dependent 0.dta in 3 of 3 processes at RT 0, 1e-295 and 1e-300; logs differing at 4 threads (c1, c2). Lane gates (dax, 8a53517): 5228 passed, 0 failed, 22 ignored over 353 test-result lines.
- port/ffap-semantics (FeatureFinderAlgorithmPicked semantics, evidence and records; fix round 1, 2026-09-16). The seed and feature stages were re-baselined on the Linux x86_64 Release build openms4-release-bc9cc12-c19e494-174b576 (Gaussian features bit for bit on glibc, EGH within 2.3038e-12; macOS arm64 recorded as a platform note). Degenerate intensity bins, short inputs and non-finite input follow that build where its outcome is measured and explained by the emitted instructions (cvttsd2si, cvttsd2si/btc, libstdc++'s probe sequence), pinned by 26 degenerate configurations, 228 probe positions and 189 non-finite cases; FeatureFinderDefs is ported against an executed probe; the intended abundance override is pinned against an adapted Release replay. Lane gates (dax, 884a2a5): 352 test binaries, 5,198 passed, 0 failed, 21 ignored.

  (F5 and F7 below supersede that lane's own list of remaining refusals.)

| Package | Merge | Evidence tier | Gates rerun by the approving verifier |
|---|---|---|---|
| `port/progress-logger-release-range` | f89d5d4 | Tier 1: the Debug-only `OPENMS_PRECONDITION(begin <= end)` (`ProgressLogger.cpp:235`) was a refusal in the port. It is removed from the wrapper and the command backend. An oracle driver linked to the Release install `openms4-release-bc9cc12-c19e494-174b576` ran 60 calls on ibminode06, 3 times (inverted, zero-width, nested and restarted ranges, end without start, NONE and default GUI), with every set/next forced past the same-second throttle and only the timing texts masked. The three runs are identical, and the port replays all 60 calls with equal outcome, depth and output bytes. The installed config.h leaves OPENMS_ASSERTIONS undefined, and libOpenMS.so has no copy of the range message | on kim: fmt, clippy `--all-targets`, rustdoc, `+1.85.0 check --all-targets`, 10 ProgressLogger-related targets 349/0/6, `--no-default-features --test progress_logger` 15/0, full `--all-features --all-targets` 5191 passed, 0 failed, 22 ignored |

The count in that row is **5191**, not the 5181 the progress lane first wrote; its
verifier recounted it.

- port/ffap-complete (merge 4709204, e562e68, 2415946; reconciliation 75cc49a; progress switch 04a14e3; doc link a74137b): the three branches merged with both lanes' tests kept (topp_feature_finder_centroided 36 passed, 0 ignored; feature_finder_picked 13; feature_finder_picked_seeds 32; feature_finder_picked_instrumentation 28; progress_logger 15); FeatureFinderAlgorithmPicked passes the inverted progress ranges of short inputs unchanged, and the_progress_event_sequence_matches_the_release_build compares every event exactly (S 5 0). Gates on kim (slot ffap-complete-gates, 04a14e3): 5,239 passed, 0 failed, 21 ignored over 353 test-result lines; rustdoc, clippy and +1.85.0 check again at a74137b.
- port/ffap-complete combined fix round 1 (60d050b, 7b40eed, ea6063d, e38534c, f6f9b58; lead decisions D1-D9 of wave 5): every std::sort FeatureFinderAlgorithmPicked reaches (spectra, chromatograms, step-1 cells, user seeds, seeds, feature map) follows libstdc++'s introsort and every std::stable_sort (the peaks of unsorted spectra and chromatograms) libstdc++ 14.4.0's, with the temporary-buffer halving, NaN and equal keys included; the overall seed score is the reference build's glibc 2.39 __powf_fma, ported from Arm optimized-routines (MIT) with the executed FMA fusion; step 1 skips scans with a non-finite drift time, as the area iterator does; step 2.5 returns std::length_error's text above vector::max_size(); ChargedIndexSet compares its index sets only; stream_number and the debug shift products print glibc's -nan with x86_64's NaN bits. Oracle ../oracle/ffap-complete-fix1 on ibminode06, every run twice and identical: sort_probe (2,272 inputs through libOpenMS.so's sortByPosition with and without a data array, sortSpectra, sortChromatograms and FeatureMap::sortByMZ, under a full, a partial and no temporary buffer), powf_probe (the dlsym-resolved powf, __powf_fma: a 42x42 grid, four sets of 2^26 pairs, every binary32 base with the exponent 1/3; the port is equal on all of them), nonfinite_stage_dt (52 stage cases: 43 returned, 8 threw, 1 never returned, all reproduced or refused at the merge; v2_rt_tie3_unsorted finds the executed 26 seeds), defs_eq_probe, and the non-finite shift runs (304 debug files byte for byte). The 9 NaN-sort cases of nonfinite_stage that the Release build returned are reproduced (170 of 170 returned runs). Gates on kim (snapshot ea6063d): fmt, clippy, eleven FFAP and progress targets, the library tests, the minimal-feature and no-default slices (3,560 passed, 0 failed, 3 ignored, counted from a piped log that probably lost a result line: the round-2 unpiped count is 3,568 with three more unit tests), +1.85.0 check, rustdoc, doctests; on dax the full --all-features --all-targets run, 5,245 passed, 0 failed, 21 ignored over 353 test-result lines; rustdoc again at f6f9b58.
- port/ffap-complete combined fix round 2 (checkpoint 23fdd7d, records 426a04c, gate record 4ce747d; the round-1 verifiers' findings): a debug run that fails in steps 1 to 2.5 keeps the opened debug/log.txt (its first line) and debug/features and leaves the stream open for the object's next run; FeatureFinderCentroided reports the step-2.5 std::length_error from TOPPBase's std::exception handler with exit 12; the parameter checks and typed members follow the source for 64-bit integer values (int narrowing with the source's InvalidParameter texts, operator unsigned int with its ConversionError half way through updateMembers_ or at the start of run_, the min_spectra cvttsd2si); getGnuplotFormula of both fitters computes its sums and products in the Release build's SSE operand order; a step-3.3.5 termination keeps the seed's log lines and feature files (source review); the introsort no longer panics on comparators that are not strict weak orderings, and its out-of-bounds guard is documented, with a proof, as unreachable for the algorithm's asymmetric comparisons; the instrumentation test compares fitted values bit for bit on Linux x86_64 (EGH within 2.4e-12), the TOPP debug input bit for bit, and the overall-score unit test against an executed row. Oracle ../oracle/ffap-complete-fix2 on ibminode06, every case twice and identical: fix2_driver bigint (21 cases), lenerr_single and lenerr_reuse (m/z 1e19: length_error; 2e18: bad_alloc; the reused object's files equal a fresh object's), formula (648 getGnuplotFormula texts), and the Release FeatureFinderCentroided on both huge inputs (exit 12). Gates (snapshot 426a04c): on kim fmt, clippy, twelve FFAP, progress and param targets (instrumentation 31, TOPP 37, feature_finder_picked 15, seeds 33), the library tests (14), the minimal-feature and no-default slices (whole no-default run 3,568 passed, 0 failed, 3 ignored over 329 result lines), +1.85.0 check, rustdoc, doctests; on dax the full --all-features --all-targets run, 5,252 passed, 0 failed, 21 ignored over 353 result lines (both totals counted from unpiped logs on the node).
- port/ffap-complete combined fix round 3 (checkpoints 6b9a5dc and 7f4b9c4, records f820e73 and e2a0661, clippy fix d562dff, gate record 7c1df72; the round-2 verifiers' findings; lead decisions D10-D12): both trace fitters call the reference build's glibc 2.39 exp and log, ported from Arm optimized-routines with the executed FMA fusion (glibc_libm), and the fit start values, EGH bounds and FWHM, the profile smoothing, cropping and quality arithmetic follow the Release build's SSE operand order, so every fit is bit for bit on Linux x86_64 and macOS arm64 (the former EGH bound 2.4e-12 and the macOS bounds 5.4e-13, 1.1e-3 and NONFINITE_FIT_GAP are gone; only the EGH area's atan on a host without glibc keeps a measured maximum, 0); step 2.5 empties windows whose binary32 bins all underflow and every window under a NaN intensity_percentage_optional, as the source's NaN weights and trimRight do; every wrapping UInt score-array count is refused whatever the Limits, and the step-1 progress range wraps as executed; the correlations divide by a zero denominator as Math::pearsonCorrelationCoefficient does; unsorted input with a mis-sized data array gives the source's Exception::Precondition text in introsort order; the empty best isotope pattern of extendMassTraces_ is refused as reachable, and a seed-loop refusal where the executed process dies or never returns records its termination with the seed's log lines. Oracle ../oracle/ffap-complete-fix3 on ibminode06, every case twice and identical: libm_probe (the host exp, log and atan: an 80-value grid and 13 sets of 2^26 inputs; exp and log equal to the port on every input on Linux and macOS, atan on Linux), gdb disassembly and table dumps of __ieee754_exp_fma and __ieee754_log_fma, fix3_stage (85 stage cases, boundary_stage.tsv.gz), vfi2_driver neg (three SIGSEGV runs with and without write_debug), fix3_driver progress (nine start events) and the Release FeatureFinderCentroided case avg0 (SIGSEGV). Gates (snapshot d562dff; e2a0661 changes only test prose, whose fmt and clippy were re-run): on kim fmt, clippy, fifteen FFAP, fitter, isotope, progress and param targets (feature_finder_picked 17, seeds 34, instrumentation 33, TOPP 38), the library tests (feature_finder_picked 18, isotopes 3), the minimal-feature and no-default slices (whole no-default run 3,572 passed, 0 failed, 3 ignored over 329 result lines), +1.85.0 check, rustdoc, doctests (67 + 3); on dax the full --all-features --all-targets run, 5,262 passed, 0 failed, 21 ignored over 353 result lines (both totals counted from unpiped logs on the node). One thread, FFC_1, best of 9: the algorithm takes 12.94 ms (symmetric) and 7.39 ms (asymmetric) on spock, against 9.61 and 6.84 ms before the round, because baseline x86_64 builds call an out-of-line fma for every fused multiply-add; macOS arm64 6.55 and 4.44 ms (6.38 and 4.27 before).
- port/ffap-complete combined fix round 4 (checkpoint 2c8e80e, records f8e86e4, gate record ddc35a7; the round-3 verifiers' findings): FeatureFinderAlgorithmPicked stores non-finite and negative FWHM, score and EGH values as setWidth and setMetaValue do (the crate-private MetaValue::source_float; the round-3 numerics verifier had shown an infinite float width reachable from large finite retention times, where the port refused the run; round 5 measured the onset on FFC_1 at a scale of 6e36); the host atan of the EGH area is used only on x86_64 Linux with glibc, and its departure elsewhere is recorded with measured rates; the reused object's isotope windows, the NaN-RT merge's NeverReturns debug output and 134 further stage cases are pinned against executed runs; the charge-wrap, empty-pattern and underflow records were corrected. Oracle ../oracle/ffap-complete-fix4 on ibminode06, every case twice and identical: fix4_stage (the round-3 numerics verifier's 124 cases, whose rows equal the verifier's capture, plus 10 width-boundary cases; extended_stage.tsv.gz), fix4_reuse (four reuse scenarios), fix4_vfi (the NaN-RT merge with and without write_debug, killed after 30 s) and the Release FeatureFinderCentroided on FFC_1 with its retention times scaled by 1e36 and 1e39 (exit 0, inf in the featureXML). Gates (snapshot f8e86e4): on kim fmt, clippy, fifteen FFAP, fitter, isotope, progress and param targets (324 passed, 0 failed, 3 ignored; feature_finder_picked 18, seeds 34, instrumentation 34, TOPP 39), the library tests (feature_finder_picked 18, isotopes 3, metadata 6), the feature slices, the whole no-default run (3,572 passed, 0 failed, 3 ignored over 329 result lines), +1.85.0 check, rustdoc, doctests (67 + 3); on dax the full --all-features --all-targets run, 5,265 passed, 0 failed, 21 ignored over 353 result lines (both totals counted from unpiped logs on the node). The extended replay runs its cases on up to eight worker threads; the feature_finder_picked target takes 187 s on kim and 129 s on macOS arm64 in the unoptimised test build (64 s on kim before the round).
- port/ffap-complete combined fix round 5 (checkpoints baa77ae and 76fd27c, records fd6f70a, gate record b5ec9e0; the round-4 verifiers' findings; lead decisions D13): the process-ending refusals outside the seed loop record their DebugTermination (the step-4 charge remainder of a caller's charge-0 feature, TerminationKind::ArithmeticTrap; a stale abort seed, OutOfBounds at TerminationPoint::AbortMap; a wrapped score-array count, OutOfBounds at ScoreArrays), for runs with and without write_debug (FeatureFinderAlgorithmPicked::termination); the instance keeps its never-closed log_ stream's counts (debug_log_file) and every termination states the length at which the executed process leaves debug/log.txt (DebugTermination::log_file_bytes), the flushed part of the debug run that opened the stream, in this run or an earlier one; the step-3.3.5 termination, found by a port-side search, is executed (a trace of zero intensities with reported_mz maximum or monoisotopic) and pinned with its debug side effects; MassTrace::avg_mz and MassTraces::intensity_profile follow the Release build's SSE NaN rules (the executed .plot prints -nan); the NaN-RT merge's seed maps, a reused instance's feature counts and the width-overflow onset (first at an RT scale of 6e36) are pinned; the stage tests compare bits under their bitwise tolerance, the sign of zero included; kernel::validate_given_finite_peaks is compiled with the mzML reader only and topp_threads' picking constants live in its Linux module, so every feature slice and macOS clippy --all-targets build without warnings. Oracle ../oracle/ffap-complete-fix5 on ibminode06, every case twice and identical but for the abort map's unique id: fix5_driver single, reuse and band, the unchanged ffap_instr_driver (stale oob and scaled) and fix4_stage (21 onset cases); 39 library cases (termination_digests.tsv.gz), 21 stage cases (width_onset_stage.tsv.gz). Gates (snapshot 76fd27c): on kim fmt, clippy, seventeen FFAP, fitter, mass-trace, isotope, progress and param targets (354 passed, 0 failed, 3 ignored; feature_finder_picked 19, seeds 35, instrumentation 37, TOPP 39), the library tests (feature_finder_picked 18, isotopes 3, metadata 6), the feature slices and their all-targets checks (no warning of this branch), the whole no-default run (3,572 passed, 0 failed, 3 ignored over 329 result lines, 0 warnings), +1.85.0 check, rustdoc, doctests (67 + 3); on dax the full --all-features --all-targets run, 5,270 passed, 0 failed, 21 ignored over 353 result lines (both totals counted from unpiped logs on the node); on macOS arm64 clippy --all-targets with all features and with none, and the changed targets.
- port/ffap-complete combined fix round 6 (checkpoint a969c27, gate record dc808a4; the round-5 verifiers' findings): where a wrapped score-array count ends the executed process is measured instead of assumed. sizeof(MSSpectrum::FloatDataArray) is 88 bytes and the array vector's max_size() is (2^63 - 1) / 88, so the only allocation between the wrap and the out-of-bounds write is one spectrum's arrays, and the same count decides differently with the memory the process may have: 100,000,003 arrays (8.2 GiB) throw std::bad_alloc under the 16 GB address space every oracle run of this branch uses and die with SIGSEGV under a 500 GB one. The port therefore records the ScoreArrays termination up to a documented 1 GiB of arrays - a crate constant no caller can move (lead decisions D6 and D12), where the earlier bound was 3 + 2 * Limits::max_charges and was falsified at 2005 arrays - and records nothing above it (the round-6 minors re-measured those counts with no address-space cap at all and moved the line; see the next paragraph, which supersedes this sentence); the debug side effects at that line are pinned (a reused object leaves the first run's flushed 1,163,782 bytes and its 79 files, a fresh debug run dies before debug/ exists). A debug run stopped by the port's own Limits::max_debug_bytes ceiling keeps its opened stream and debug/features, as the source leaves them, instead of losing its DebugOutput. MassTraces::update_baseline promotes its f32 with the emulated cvtss2sd, like MassTrace::avg_mz and MassTraces::intensity_profile, so a NaN baseline carries the executed bits into the stored score_fit and score_correlation and into the .plot formula; the remaining f32 promotions of the path use f64::from, measured equal on macOS arm64 for ten NaN patterns and compiled as cvtss2sd on x86_64 (a platform note). The EGH area's tolerance is recorded as equal to BITWISE, so it relaxes nothing. Oracle ../oracle/ffap-complete-fix6 on ibminode06, every case twice and identical: fix6_driver (= fix5_driver with the round-5 verifier's wrapLH variant and a sizes mode), node/run_wrap.sh with 20 wrapped counts from 1 to 2^32 - 5 arrays under a 16 GB and a 500 GB address space, node/run_cases6.sh with the three cases at 12,201,611 arrays; 20 rows (score_array_wraps.tsv) and 42 library cases (termination_digests.tsv.gz). Gates (snapshot a969c27; dc808a4 after it changes only documentation and the two FFAP manifests): on kim fmt, clippy (0 warnings), seventeen FFAP, fitter, mass-trace, isotope, progress and param targets (356 passed, 0 failed, 3 ignored; feature_finder_picked 19, seeds 35, instrumentation 39, helper_structs 26, TOPP 39), the library tests (feature_finder_picked 18, isotopes 3, metadata 6), the feature slices (158, 93, 41, 41) and their all-targets checks (no warning of this branch), the whole no-default run (3,572 passed, 0 failed, 3 ignored over 329 result lines, 0 warnings), +1.85.0 check, rustdoc, doctests (67 + 3); on dax the full --all-features --all-targets run, 5,272 passed, 0 failed, 21 ignored over 353 result lines (both totals counted from unpiped logs on the node); on macOS arm64 clippy --all-targets with all features and with none, rustdoc, the changed targets (feature_finder_picked 19 in 229 s, instrumentation 39, seeds 35, helper_structs 26, TOPP 39 and eight more) and cargo check --lib for four feature slices.
- port/ffap-complete, the round-6 minors (checkpoint 1a14793, gate record a11fc26; the round-6 verifiers' three minor findings): the score-array recording line is moved to what the reference platform measures, and the cvtss2sd promotion of MassTraces::update_baseline is pinned by executed bits instead of a bool. Round 6 had set SCORE_ARRAY_TERMINATION_CEILING_BYTES to 1 GiB on runs made under the 16 GB `ulimit -v` the oracle harness imposes; re-measured with round 6's own driver binary and no cap at all, every wrapped count from 12,201,611 to 1,000,000,003 arrays dies with SIGSEGV on the reference node, the five counts that threw std::bad_alloc under the cap included, so the line is now the bytes of 1,000,000,003 arrays and the conservative band round 6 documented is empty; the fixture carries the nine uncapped rows beside the capped ones (address_space_kib `none`) and the test asserts both halves - nothing recorded above the line, nothing measured to die left unrecorded - and fails at the old line. The same runs corrected the model: the pattern loop names and assigns every in-bounds array before the first out-of-bounds index, so the process holds about 232 bytes per array and not the 88 the array costs (maximum resident set 2.88 GiB at 12,201,611 arrays, 21.85 GiB at 100,000,003), which is why 2^32 - 5 arrays (about 928 GiB, 93% of the node) stays above the line and unrun. MassTraces::update_baseline's promotion is pinned against 90 executed baselines - 18 f32 patterns (quiet, signalling, negative and maximal-payload NaNs, both zeros, both infinities, both extremes, the two smallest subnormals, the smallest normal and an ordinary value) in five peak layouts - where the only NaN case in the repository had recorded is_nan(); the test states what it does not claim, namely that f64::from gives the same bits on the measured hosts, so it pins the value and not the instruction. Oracle ../oracle/ffap-complete-min6 on ibminode06, every case twice and identical: baseline_driver (node/run_baseline.sh, 90 rows), fix6_driver reused unchanged (node/run_native_wrap.sh, nine uncapped counts; node/run_mem.sh, two /usr/bin/time -v runs); 29 rows (score_array_wraps.tsv) and 90 rows (feature_finder_picked_helper_structs_update_baseline.tsv). Gates (snapshot 1a14793; a11fc26 after it changes only the three manifests' target_verification): on dax fmt, clippy (0 warnings), twelve FFAP, fitter, mass-trace and progress targets, the library tests, the feature slices, +1.85.0 check and rustdoc; on macOS arm64 fmt, clippy --all-targets with all features, the changed targets and the repository's Python checks (check_core_sdk, check_module_cycles, check_schema_feature_graph; doc coverage 4517/5872 = 76.9%, unchanged).

- **Signal-to-noise (port/signal-to-noise).** SignalToNoiseEstimatorMedian.h and SignalToNoiseEstimator.h are complete (docs/SIGNAL_TO_NOISE_SUPPORT.md).
  - The unmodified Linux x86-64 Release build openms4-release-bc9cc12-c19e494-174b576 ran on ibminode06, each case twice in a fresh process, byte-identically.
    - 143 cases: 89 estimator cases with every ratio, max_intensity_, both percentages, warnings and CMD progress output; 47 random-scan cases with the seed set through an interposed time(); 5 engine cases; 390 whole nth_element permutations; and PeakPickerHiRes::pick through libOpenMS's own instantiation. The fixtures are in tests/data/signal_to_noise/.
    - Five AUTOMAXBYPERCENT cases on 2^31 and 3,000,000,001 points (../oracle/sne-fix) all take the negative-range return, which pins the 32-bit cvttsd2si at :220. The port's full estimation matches them outside CI.
  - The estimator is a header template, so the driver instantiates it with libOpenMS's own compile flags. Its 88 floating-point instructions match libOpenMS's copies, identically for one instantiation and up to one memory displacement for the other. Both copies emit the same four 32-bit cvttsd2si.
  - The source profile reproduces two undefined-behaviour outcomes of that build, both measured conversions: CPP-257's bin 0 and :220's INT_MIN. It refuses exactly the out-of-bounds and signed-overflow sites.
  - The unchanged P1 driver re-run on Release printed records byte-identical to the arm64 Debug fixture.
  - The emulated :49-50 random-scan pointer wrap in 10 cases (../oracle/sne-followup, same driver
    binary 571d6f9b as sne-completion, each twice byte-identically): the drawn position wraps to an
    in-bounds element e = idx mod 2^62, which the port computes; the disassembly re-dump confirms the
    :48 conversion and the `lea (%r12,%rcx,4),%r15` pointer reused for the read.
  - Full suite: 5213 passed, 0 failed, 22 ignored.

### The six fix rounds and their verdicts

`port/ffap-complete` was verified after every round by **two independent lenses**, one
on numerics and one on instrumentation, each on its own detached checkout:

| Round | Head | Numerics lens | Instrumentation lens |
|---|---|---|---|
| 1 | `f6f9b58` | changes_required | changes_required |
| 2 | `4ce747d` | changes_required (averagine underflow, NaN `intensity_percentage_optional`, the EGH bound is a measured maximum) | approve_with_notes |
| 3 | `7c1df72` | approve_with_notes | changes_required (`DebugTermination` missing at two process-ending refusals outside the seed loop) |
| 4 | `ddc35a7` | approve_with_notes | changes_required |
| 5 | `b5ec9e0` | approve_with_notes | changes_required (the `score_arrays_overrun` bound was not supported by its evidence) |
| 6 | `dc808a4` | **approve_with_notes** | **approve_with_notes** |
| round-6 minors | `a11fc26` | (three minor findings applied; the ceiling re-measured uncapped) | |

`port/signal-to-noise`: round 1 changes_required (one major — the `> i32::MAX` refusal
was too broad), round 2 approve_with_notes, follow-up **approve**.
`port/progress-logger-release-range`: approve_with_notes, its one test-integrity minor
fixed by the lead in `f89d5d4`.

Every round-6 recommendation was to promote. This pass applied all three.

### Lead decisions D1-D17

The full text of D1-D13 is in
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-5-status); of D14-D17, in
[wave 8's own list](#lead-decisions-d1-d17) above. In one line each:

- **D1** reproduce a measured, repeatable, instruction-explained, **in-bounds** Release
  outcome; refuse out-of-bounds, races, termination and endless loops.
- **D2** libstdc++'s binary-search probes on NaN keys.
- **D3** every `std::sort` as introsort, every `std::stable_sort` as libstdc++ 14.4.0's.
- **D4** glibc's `-nan` in FFAP's and the trace fitters' formatters; the five others are
  reported, not changed.
- **D5** the reference build's glibc `powf`, ported licence-clean from Arm
  optimized-routines; reproduced exactly, so no fallback was needed.
- **D6** `std::length_error`'s text above `vector::max_size()`, native ceiling below.
- **D7** `ChargedIndexSet` equality by index sets.
- **D8** `RejectedParameters::Shown` default; the heap-address bound scoped to the
  reference platform.
- **D9** FAIMS out of scope (a FeatureFinderCentroided decision).
- **D10** `exp` and `log` ported with their FMA fusion; `atan` takes the fallback,
  because `__atan_fma` is LGPL-only with no licence-clean upstream.
- **D11** the abort race gives the single-thread result — **the one accepted exception
  to D1**.
- **D12** wraps that lead out of bounds are refused at the wrap, not behind a raisable
  ceiling; in-bounds wraps are reproduced. Round 6 made the recording line a crate
  constant; the minors moved it to the largest count measured to die **uncapped**.
- **D13** `MetaValue::source_float` accepted; the featureXML writer and CLI text split
  off — and since closed on `fix/featurexml-nonfinite`: the writer writes the source's
  `inf`/`-inf`/`NaN`, the reader takes them back, a failed store is the source's write
  failure with exit 5, and TOPP native difference 16 is closed; the longer test time
  accepted with no assertion dropped; charge counts `-2`/`-3` refused unconditionally;
  `atan` keeps the `libm` crate off the reference platform.
- **D14** the mzML reader accepts a dangling `sourceFileRef` under
  `source_dangling_references`, as the Release build does, and keeps today's refusal as
  the strict default (wave 8; the full paragraph is in that wave's section).
- **D15** a picker's internal noise estimate reproduces the source unconditionally where
  the signal it reads is picker-generated or already validated, and follows the picker's
  own profile where the estimator reads the caller's data (wave 8).
- **D16** reproducing an unspecified `std::sort` permutation is in scope, because the
  port already does it under tier-1 validation, because doing so is D1-compliant by
  construction, because refusing would turn away ordinary finite data, and because the
  libstdc++ header sha256s make a toolchain change detectable (shared-math wave,
  2026-09-19).
- **D17** `sort_ascending` may take a proved-equivalent fast path — an in-place
  `f64::total_cmp` sort that allocates nothing — where the sample holds no NaN and not
  both spellings of zero, because there every class `operator<` calls equivalent is a
  set of bit-identical values and the permutation is therefore unobservable (shared-math
  repair round, 2026-09-20).

### This pass's gates

All on `kim` through the gate script, slot `integ-w5`, detached, logs in the session
scratchpad under `integ-w5-logs/`. Every figure below is summed from the log's
`test result:` lines, over **all** of them, so a dropped target cannot hide.

| Gate | Result |
|---|---|
| `+1.85.0 check --locked --all-features --all-targets` | exit 0 (MSRV 1.85, no let-chains) |
| `clippy --locked --all-features --all-targets -- -D warnings` | exit 0, **0 warnings** |
| `doc --locked --all-features --no-deps`, `RUSTDOCFLAGS=-D warnings` | exit 0 |
| `test --locked --all-features --all-targets` | **5297 passed, 0 failed, 21 ignored** over 354 result lines, 327 integration targets + 27 unit-test binaries |
| `test --locked --no-default-features` | **3599 passed, 0 failed, 3 ignored** over 330 result lines |
| `+1.85.0` minimum-rust line :95 (`--no-default-features`, + `signal_to_noise`) | 80 passed, 0 failed, 0 ignored over 5 result lines |
| `+1.85.0` minimum-rust line :113 (`mzml paramxml`, + `signal_to_noise`) | 143 passed, 0 failed, 2 ignored over 6 result lines |
| `+1.85.0` minimum-rust line :114 (`mzml paramxml featurexml`, + `feature_finder_picked_instrumentation`) | 165 passed, 0 failed, 0 ignored over 5 result lines |

Locally on macOS arm64: `cargo fmt --all -- --check` exit 0, and all seven repository
Python checkers pass — `check_core_sdk` (plain **and** with
`--source .reference/openms4-core-bc9cc12`, 2,092 distinct source/registration/reference
files verified at `bc9cc12`), `check_doc_coverage`, `check_module_cycles`,
`core_sdk_coverage`, `test_core_sdk`, `test_core_sdk_coverage` and
`check_schema_feature_graph`. `.github/workflows/rust.yml` parses as YAML and every
changed JSON record round-trips through `json.load`.

No gate needed a rerun: none exited 255.

### CI audit

`.github/workflows/rust.yml` gains exactly three `--test` names, all in the
`minimum-rust` job, and each was confirmed on `+1.85.0` above: `signal_to_noise` on the
`--no-default-features` line (:95) and on the `"mzml paramxml"` line (:113), and
`feature_finder_picked_instrumentation` on the `"mzml paramxml featurexml"` line (:114).
The ProgressLogger lane needed none: `progress_logger` is already on the
`--no-default-features` line at :94.

**No test binary in `tests/` is unrun.** There are 327 integration targets on
disk. Every `--test` name in the workflow resolves to one of them (**0 dangling**), and
the `test --locked --all-features --all-targets` gate above launched
**327** of them, i.e. all of them. 176 of the 327 are
named on a `--test` line; the other 151 are reached only
by the three `--all-features --all-targets` steps, which is sufficient because all
80 file-level `#![cfg(...)]` gates in `tests/` are cargo feature gates on
features declared in `Cargo.toml`, so `--all-features` satisfies every one of them. The
`--test` lines exist to prove the feature slices, not to reach otherwise unreachable
targets.

### Ignored tests

22 `#[ignore]` attributes in the tree, all of them in `tests/` and none in `src/`, down
from 23 on `main`. The gate above reports **21** ignored rather than 22 because
`tests/lm_eigen_path_differential.rs`'s `macos_arm64_sdk_gap_report` is additionally
`#[cfg(all(target_os = "macos", target_arch = "aarch64"))]`, so it is not compiled on
`kim`; the two numbers agree once that is accounted for. **The one that went is the point of this wave**:
`tests/topp_feature_finder_centroided.rs`'s
`a_zero_width_retention_time_range_diverges_from_the_cpp_release_build`, whose reason
read "documented divergence: the port refuses a zero-width RT range that C++ Release
carries through to an empty feature map". The divergence is closed — the port now
follows the Release build — and three running tests replace it
(`a_zero_width_retention_time_range_follows_the_cpp_release_build`,
`a_zero_width_mz_range_follows_the_cpp_release_build`, and
`a_short_input_never_reaches_the_seed_loop_as_in_the_cpp_release_build`, which is no
longer ignored either). No test was ignored, skipped or weakened by this wave; **no new
`#[ignore]` was added anywhere on either branch**.

Every remaining one, with its reason and owner:

| Attribute | Reason on the attribute | Owner / why it is not a coverage gap |
|---|---|---|
| `tests/featurexml.rs:865` | HPC scale: reads the 59.6 MiB and 2.06 GiB benchmark featureXML files by path | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/file_info.rs:168` | mzML reader gaps outside A4: duplicate userParam 'name', processing on primary arrays, float charge array | mzML reader owner (outside package A4). A **real, named** coverage gap, each naming the exact construct and, where it applies, the C++ behaviour it does not yet match |
| `tests/file_info.rs:210` | mzML reader gap outside A4: the dangling spectrumList defaultDataProcessingRef 'dp_sp_0' is refused (unresolved dataProcessingRef); P2's source-compatible option (D10) | mzML reader owner (outside package A4). A **real, named** coverage gap, each naming the exact construct and, where it applies, the C++ behaviour it does not yet match |
| `tests/file_info.rs:325` | mzML reader gaps outside A4: duplicate userParam 'name', processing on primary arrays, float charge array | mzML reader owner (outside package A4). A **real, named** coverage gap, each naming the exact construct and, where it applies, the C++ behaviour it does not yet match |
| `tests/file_info.rs:351` | mzML reader gap outside A4: 'charge array' stored as 64-bit float is refused (canonical auxiliary array binary type); C++ converts it | mzML reader owner (outside package A4). A **real, named** coverage gap, each naming the exact construct and, where it applies, the C++ behaviour it does not yet match |
| `tests/file_info.rs:494` | mzML reader gap outside A4 (A3 request 5): C++ copies the selected-ion drift time 8.1 onto the MS2 spectrum (MzMLHandler.cpp:1871-1875), so its ion-mobility ranges end at 8.10; the Rust reader does not | mzML reader owner (outside package A4). A **real, named** coverage gap, each naming the exact construct and, where it applies, the C++ behaviour it does not yet match |
| `tests/fuzzy_string_comparator.rs:1197` | fills the 256 MiB log buffer | logging owner. A resource cost (256 MiB), not a behavioural gap |
| `tests/gauss_trace_fitter.rs:1593` | macOS-generated oracle against a solver that matches Linux x86_64 Release Eigen; prints a report, asserts no Rust value | the owning lane. Prints a report and asserts no ported value, so it cannot mask a regression |
| `tests/lm_budget_differential.rs:1257` | B3-LM gate report for the rejected levenberg-marquardt candidate; run with --ignored --nocapture | the owning lane. Prints a report and asserts no ported value, so it cannot mask a regression |
| `tests/lm_eigen_path_differential.rs:557` | measures the cost of matching Linux x86_64 Release on macOS arm64; asserts nothing | the owning lane. Prints a report and asserts no ported value, so it cannot mask a regression |
| `tests/mzml_reader_scale.rs:840` | reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/mzml_reader_scale.rs:846` | reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/mzml_reader_scale.rs:856` | reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/mzml_reader_scale.rs:867` | reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/mzml_reader_scale.rs:877` | reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/mzml_reader_scale.rs:906` | reads /ceph/ibmi/abi/oliver/bench/openms4/inputs on the IBMI nodes | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/mzml_writer_scale.rs:466` | HPC only: reads the 547 MB UK222_picked benchmark input from /ceph | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/mzml_writer_scale.rs:475` | HPC only: reads the 2.3 GB UK222 profile benchmark input from /ceph | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/topp_baseline_filter_edges.rs:293` | reads /ceph/ibmi/abi/oliver on the IBMI nodes | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/topp_baseline_filter_edges.rs:304` | reads /ceph/ibmi/abi/oliver on the IBMI nodes; multi-GB | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/topp_threads.rs:930` | HPC: reads /ceph/ibmi/abi/oliver/bench/openms4/inputs/derived; run on an IBMI node | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |
| `tests/topp_threads.rs:954` | HPC: reads multi-GB inputs under /ceph/ibmi/abi/oliver/bench/openms4/inputs | benchmark lane. Reads staged multi-gigabyte input on an IBMI node; unrunnable in CI by construction, and run by hand on the node |

**None of them hides a wave-5 behaviour.** 13 read staged multi-gigabyte inputs on
an IBMI node and cannot run in CI by construction; 5 are named mzML-reader gaps
outside package A4, each stating the exact construct it does not yet read and what the C++
does instead, so they are a visible backlog rather than a silent exemption; and the
remaining 3 print reports and assert no ported value; the fourth,
`verbose_3_log_bytes_are_bounded`, is skipped only for its 256 MiB memory cost and
does assert the ported log bound, as the table row above records. So none of them
could mask a regression.

Two of these files **were** touched by this wave — `tests/gauss_trace_fitter.rs` and
`tests/topp_threads.rs` — but neither diff adds, removes or edits an `#[ignore]`
attribute or the body of an ignored test: `gauss_trace_fitter.rs` gained trace-fitter
prose and `topp_threads.rs` moved its picking constants into its Linux module.

### What this checkpoint does not claim

- It does not claim a performance result on the wave-4 timing node. §4 of
  [BENCHMARKS](BENCHMARKS.md) ran on **dax**, not ibminode05, and its absolute numbers
  are comparable only with each other.
- It did not settle the `-C target-feature=+fma` question; nothing in the build
  configuration was changed by that wave. The user decided it on 2026-09-18 and
  `port/fma-default` carries it: x86_64 builds set the flag, output is unchanged,
  and a processor without FMA is refused with exit 12 rather than `SIGILL`
  ([FMA_BUILD_FLAG](FMA_BUILD_FLAG.md)).
- It does not claim the signal-to-noise signed-overflow sites are emulated. They stay
  refused as a stated cost/benefit decision; each needs more than `2^31` points and
  64-90 GB per evidence run.
- It does not claim `2^32 - 5` score arrays was executed uncapped: on the measured
  model that needs about 928 GiB, 93 % of the shared reference node, and it was
  deliberately not run.

## Wave-4 performance and correctness integration (2026-09-16)

`main` (`9a392fe`) carries, on top of wave 3's `fabd4b9`, 26 commits from eight
lanes: three performance lanes (`perf/peak-picker`, `perf/mzml-reader`,
`perf/spline-scratch`), one that removed dead work (`perf/validation`) and four
correctness fixes (`fix/map-normalizer`, `fix/tool-limits`,
`fix/featurexml-scale`, `fix/dta-precision`), plus the wave-3 leftovers
`fix/ffc-integration` and `fix/picked-chromatogram`. This shared-file pass on
`integrate/wave4-shared` records them in CI, the ledger, the provenance files,
the C++ issue log and the documentation, and replaces
[BENCHMARKS](BENCHMARKS.md) with the wave-4 run. See
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-4-status).

Tiers as before: tier 1 is an executed differential against a C++ oracle, tier 3
upstream class-test literals, tier 4 native derivation. The Release build
`openms4-release-bc9cc12-c19e494-174b576` (core `bc9cc12`, cli `c19e494`, topp
`174b576`) is the executed reference for everything on real data, and its
identity is in [BENCHMARKS](BENCHMARKS.md) §1. Every count below is from the
package's approving verifier, rerun on its own detached checkout through the
gate script on a Linux x86_64 IBMI node.

| Package | Merge | Evidence tier | Gates rerun by the approving verifier |
|---|---|---|---|
| `perf/peak-picker` (3 review rounds) | `2491afc` | Tier 4 with a bit-identity gate: the picking loop is parallel behind the `parallel` feature, and the output is **byte-identical at 1, 2, 4, 8 and 32 workers** and to the serial result. 1.66x at 32 workers on the lane's own 2.3 GB measurement, peak RSS −886 MB. The `-threads` policy reaches the picker through `ToolContext`, and `tests/topp_threads.rs` covers it as its sixth tool, on a profile input, holding it to **no** worker at `-threads 1` | on dax: `--all-features --all-targets`, `--no-default-features --features mzml,paramxml` picker targets, clippy `--all-targets`, rustdoc, fmt, `+1.85.0 check --all-targets` — all exit 0; three release builds of the binary; the pinned digest `bb13eecf…` reproduced |
| `perf/mzml-reader` | `821e783` | Tier 4, instruction-counted: buffer reuse across records, `decode_slice` in place of `decode_vec`'s zero fill, and a base64 accumulation that no longer pushes one `char` at a time, for **−47.2 %** of the instructions the load path executed. No decoded value changes; the pinned picker digest reproduces | callgrind on a fixed slice with a layout control; the full 2.3 GB input through three binaries; the lane's whole gate set on the branch tip |
| `perf/spline-scratch` | `e766311` | Tier 4 with a bit-identity argument: `CubicSpline2dFitter` replaces the eight heap vectors `with_max_points` allocated per spline. The recurrence is unchanged term by term and in its original evaluation order, so the coefficients are bit-identical **by construction**; `alloc::alloc` was 271,144,136 of the 590,918,793 instructions the construction cost. Measured wired into the picker: 33.19 s → 32.01 s at one thread, pick phase 17.24 s → 15.11 s, output sha256 unchanged | 11 interleaved A/B pairs on a quiet node; the branch's own gate set; the picker rebuilt with the fitter wired in and byte-identical at 1, 8 and 32 threads |
| `perf/validation` (3 rounds) | `4389b96` | Tier 4 with a mutation-checked premise: the reader's per-record `spectrum.validate()` is **removed as dead**, both the m/z and the intensity half, over all 197,765,338 points of the benchmark input, because the decoder already refuses a nonfinite value before a `Peak1D` or `ChromatogramPeak` is constructed. −44,121,898 instructions (15.0 Ir per peak); the branch as a whole −60,947,537 (−1.605 %) against its merge base, layout control 0 | on dax: fmt, `nextest --all-features --test kernel --test mzml --test peak_picking` 58/58 on 1.96 and 1.85, clippy `--all-targets`, rustdoc, `+1.85.0 check --all-targets` — all exit 0; 27 callgrind runs over 9 arms, reps agreeing to the instruction; the bit-identity gate reproduced five ways including from the merge base |
| `fix/map-normalizer` (2 rounds) | `c2ecade` | Tier 1 on 1.2 GB of real data: the tool normalised against the spectrum maximum where the source uses the combined maximum including chromatograms — 13.24x on 7,302 of 87,492 intensity arrays. All 87,492 now agree **bitwise** with the C++ Release tool's. The empty-range refusal is restored and was **executed** against the C++ rather than reasoned (the earlier reasoning was refuted by the execution) | on dax: the full 1.2 GB run against the C++ with one shared INI, 88,478,237 intensity and 88,434,492 m/z points, zero differences, port output identical at 1 and 32 threads; four degenerate inputs executed on the C++; the lane's gate set on both toolchains |
| `fix/tool-limits` | `5856c63` | Tier 1 on 1.2 GB of real data: `MzMLSplitter` and `SpectraFilterWindowMower` process full-size input for the first time. The writer refused precursor references that do not resolve inside a split part (`CPP-311`); the window mower applied its 1,000,000-point cap to the whole map instead of per spectrum. All four split parts and all 87,492 arrays bitwise equal to the C++'s | `real_data_differential` recorded in both tools' provenance manifests; the lane's gate set; `src/processing/window_mower.rs` documented at 100 % |
| `fix/featurexml-scale` | `23b6519` | Tier 1 with executed C++ at scale: the fixed ~12.5 MB decode ceiling (three limits combined with `min()`) is replaced by size-derived ceilings in the new `src/format/featurexml_scaling.rs`, and features are streamed. Both benchmark maps load and round-trip (59.6 MiB / 42,789 features at 150 MiB peak, 2.06 GiB at 5.39 GiB peak); the C++ Release `FileInfo` was executed on both and its report is byte-identical to the port's apart from the C++ timing footer | the lane's gate set on both toolchains; the shared `src/format/identification_xml.rs` changes proved additive by running idXML, consensusXML, mzIdentML and map_xml targets unchanged |
| `fix/dta-precision` | `db00dbf` | Tier 1, byte-equal on 836 MB: the writer reproduces the source's two 15-digit numeric rules, so `DTAExtractor` over the 1.2 GB Velos run writes 36,443 files and exactly 836,505,793 bytes, byte-equal to the C++ Release tool's. Before the fix it wrote 658,901,890 bytes (−21.2 %) | two full 836 MB output trees compared file by file on ibminode06; five interleaved repetitions per side; the lane's gate set on dax |
| `fix/ffc-integration`, `fix/picked-chromatogram` | carried from wave 3 | recorded in the wave-3 checkpoint below; the wave-4 benchmark confirms both at full size | — |

**Lead's results for the merged tree** (dax, detached): `test --locked
--all-features --all-targets` **5,189 passed, 0 failed, 22 ignored** on stable
and on `+1.85.0`; `+1.85.0 --no-default-features --all-targets` **3,500
passed, 3 ignored**; **70 doctests**; `clippy --all-targets -D warnings` and
rustdoc `-D warnings` clean; all six local Python checkers (`check_core_sdk`,
`check_doc_coverage`, `check_module_cycles`, `core_sdk_coverage`,
`test_core_sdk`, `test_core_sdk_coverage`) pass. **All four figures were
reproduced on kim by this pass** on its own tree — 5,189/0/22, 3,500/0/3 and
67 + 3 = 70 doctests, to the test — so the checkpoint quotes a measurement, not
a hand-over.

### Instrument-scale comparison

[BENCHMARKS](BENCHMARKS.md) is rewritten around the wave-4 run of 2026-09-16:
**all eight ported TOPP tools at 1 and 32 threads** — seven on full-size
instrument data, `FeatureFinderCentroided` on the documented 4,000-spectrum
subset, because neither implementation finishes the full run — on a quiet node, with 0 of 192 repetitions load-flagged. The wave-3
single-pair picker measurement it used to carry is superseded. In one line: at
one thread the port is faster on `SpectraFilterWindowMower` (0.73–0.77) and
`FileInfo`-on-mzML (0.803), level on `PeakPickerHiRes` (1.010, inside the ~3 %
band), and 1.10x–1.58x slower on the other six; at 32 threads it is 1.88x faster
on `PeakPickerHiRes` (0.531) and near-level on `FeatureFinderCentroided` (1.060)
while five tools ignore the flag; and it agrees with the C++ Release build on
the **data** of every tool — three byte-equal outputs and five mzML writers with
every decoded array bitwise identical.

Three things about that document are worth repeating here, because the
benchmark runner's own summary got them wrong and the adversarial review
corrected them:

- the 0.4–5.7 % output-byte deficit is **92–93 % XML indentation the port does
  not write** (the port's files contain zero tab characters), not metadata;
  `dataProcessingRef` appears **once per C++ file**, not once per spectrum;
- the C++ "32 threads" cells run 64 live threads, but the surplus is an **idle
  OpenBLAS pool sized by the harness's own `OMP_NUM_THREADS`** — with the
  variable unset the same binary peaks at 129 threads at `-threads 1` — so it is
  not a doubled compute budget and must not be framed as one;
- the SHA-1 `fileChecksum` cost shares are an **arithmetic projection** from a
  measured 793 MB/s throughput divided into the output size. No no-hash build
  was ever timed, so they are not a measured ablation.

The full-size `FeatureFinderCentroided` pilot is also stated asymmetrically on
purpose: the run's own log records the **C++** side killed at a 600 s cap and
holds no wall or exit record for the Rust side; the Rust result rests on the
reviewer's independent re-run at a 700 s cap, which was killed at the cap with
no output.

### The 22 ignored tests

All 22 on Linux with `--all-features` (23 `#[ignore]` attributes exist; the
macOS-only `macos_arm64_sdk_gap_report` is compiled out there), which matches
the lead's 22. **Exactly one is new in this window**, marked **NEW** below; the
other 21 are carried unchanged from wave 3 and their owners are unchanged.

| Test | Reason | Documented gap or external dependency | Owner |
|---|---|---|---|
| `featurexml.rs::hpc_scale_benchmark_featurexml_files_load_and_round_trip` **NEW** | reads the 59.6 MiB and 2.06 GiB benchmark featureXML files by absolute path | external dependency (`/ceph` staged inputs); the behaviour it exercises is documented in [FEATUREXML_SCALE_SUPPORT](FEATUREXML_SCALE_SUPPORT.md) | `fix/featurexml-scale` |
| `fuzzy_string_comparator.rs::verbose_3_log_bytes_are_bounded` | fills the 256 MiB log buffer | resource cost, not a gap; the bound it checks is asserted | C3-FUZZY |
| `lm_budget_differential.rs::levenberg_marquardt_crate_candidate_gate_report` | the gate report of the rejected crate candidate; asserts nothing | [THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md), Levenberg-Marquardt row | B3-LM |
| `gauss_trace_fitter.rs::solver_gap_probe_reports_the_known_gap` | a macOS-generated oracle against a solver that now matches Linux x86_64 Release Eigen; prints a report, asserts no Rust value | [TRACE_FITTER_SUPPORT](TRACE_FITTER_SUPPORT.md) "Known gap"; [DISTRIBUTION_FITTERS_SUPPORT](DISTRIBUTION_FITTERS_SUPPORT.md) §1 | B4-GAUSS / B3b |
| `lm_eigen_path_differential.rs::macos_arm64_sdk_gap_report` (macOS arm64 only, **not** in the 22) | measures the cost of matching Linux x86_64 Release on macOS arm64; asserts nothing | the same platform decision | B3b |
| `topp_feature_finder_centroided.rs::a_zero_width_retention_time_range_diverges_from_the_cpp_release_build` | documented divergence: the port refuses a zero-width RT range that the Release build carries through to an empty feature map | `CPP-274` and now `CPP-312`, which records the C++ side as executed; [TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT](TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md) native difference 2 | B10 |
| `file_info.rs::c1_file_info_9_mzml_mps`, `::a4_file_info_9_default_flags` | the strict mzML reader refuses `FileInfo_9_input.mzML` (a repeated spectrum userParam `name`, `dataProcessingRef` on the m/z and intensity arrays, a 64-bit float charge array) | [FILE_INFO_SUPPORT](FILE_INFO_SUPPORT.md) "Known reader gaps", decision D10 | mzML reader owner |
| `file_info.rs::a4_indexed_file_info_12_all_flags` | a 64-bit float `charge array` in `FileInfo_12_input.mzML` | the same | mzML reader owner |
| `file_info.rs::c1_empty_mzml_mps` | the dangling `defaultDataProcessingRef` of `empty.mzML` | D10; the option exists and A5's tool passes it, but `file_info::Options` still defaults strict | A6 |
| `file_info.rs::a4_mzml_file_1_all_flags` | the selected-ion drift time is not copied onto the MS2 spectrum | A3 request 5; the lead decided on 2026-09-15 to follow the executed source, lane still not opened | the lead |
| `mzml_reader_scale.rs` × 6 (`hpc_*`) | read `/ceph/ibmi/abi/oliver/bench/openms4/inputs` on the IBMI nodes | external dependency (multi-GB staged inputs) | `fix/mzml-reader-scale` |
| `mzml_writer_scale.rs::hpc_benchmark_centroid_uk222_picked_stores_through_the_tool_path`, `::hpc_benchmark_profile_uk222_stores_through_the_tool_path` | read the 547 MB and 2.3 GB benchmark inputs from `/ceph` | external dependency | `fix/mzml-writer-scale-parity` |
| `topp_baseline_filter_edges.rs::uk222_first600_matches_the_release_tool`, `::uk222_full_matches_the_release_tool` | read `/ceph/ibmi/abi/oliver` on the IBMI nodes; the second needs about 25 GB of memory | external dependency; the C++ side is a lane-private oracle directory | `fix/baseline-filter-last-point` |
| `topp_threads.rs::hpc_benchmark_slices_are_thread_invariant`, `::hpc_full_size_inputs_are_thread_invariant` (Linux + `parallel` only) | read the staged benchmark inputs under `/ceph` | external dependency | `fix/tool-threads` |

Every reason names a documented gap or an external dependency, and no ignore
hides an unexplained failure. Counted on the 22 Linux rows: **three reports**
meant to be read with `--ignored --nocapture` (one of them,
`verbose_3_log_bytes_are_bounded`, does assert its bound and is listed only
because it is ignored for its 256 MiB cost); **thirteen** need the IBMI `/ceph`
share (`mzml_reader_scale` × 6, `mzml_writer_scale` × 2,
`topp_baseline_filter_edges` × 2, `topp_threads` × 2 and the new `featurexml`
row); **five** name a reader gap (the `file_info` rows) and **one** a decided
divergence. Three plus thirteen plus five plus one is the 22 that
`--all-features --all-targets -- --ignored --list` prints. The wave-2 tripwire
`reader_gaps_behind_the_ignored_cases_are_still_present` still fails when a
reader gap closes.

The count moved 21 → 22 because of the one new row and nothing else: no ignore
was added to hide a failure, and none was removed. `fix/featurexml-scale` was
the only lane in this window that added a test needing an external input.

### CI

No lane in this window added a test **binary**, so no target was missing from
the feature-sliced lines. The audit that matters — *is any test binary in
`tests/` unrun?* — was redone from the tree rather than from the previous
answer: 325 test binaries exist, all 325 are built and run by `cargo test
--locked --all-features --all-targets`, which is the first step of both the
`test` and the `minimum-rust` job, and 174 of them are additionally named on a
reduced-feature line — 169 before this pass's four lines, 174 after. Nothing is unrun.

The audit did turn up three reduced-feature gaps that this window's work makes
worth closing, and four lines were changed or added:

| Line | Job | Why |
|---|---|---|
| `--no-default-features --features mzml … --test mzml …` (added `--test mzml`) | minimum-rust | `tests/mzml.rs` is gated `#![cfg(feature = "mzml")]` and was named on no reduced-feature line, so the rewritten reader and the removed validation loop were exercised only in the `--all-features` build |
| `--no-default-features --test kernel --test spline_math --test peak_picking --test window_mower` (new) | minimum-rust | the four library targets this window rewrote. They are not feature-gated, so the `test` job's bare `--no-default-features` line covers them on stable; the `minimum-rust` job has no such bare line, and covered them only through `--all-features` |
| `--no-default-features --features "mzml paramxml parallel" --test topp_threads --test topp_peak_picker_hi_res --test peak_picking_experiment` (new) | minimum-rust | the parallel picker is behind `#[cfg(feature = "parallel")]`, and `--all-features` was the **only** build that compiled it. The existing `mzml paramxml` line exercises the `not(parallel)` arm, so the two lines together now cover both |
| the same parallel line | test | the same slice on stable |

Each new or changed line was run on `+1.85.0` on kim; the results are in the
gate table below.

### The measurement hazard this window documented

One finding from `perf/validation` outlives its lane and applies to **every
future comparison of two commits of this port**, so it is recorded here and in
[EARLY_TOPP_WORK_PACKAGES](EARLY_TOPP_WORK_PACKAGES.md) rather than left in a
branch report.

`Record::finish` in `src/format/mzml.rs` has exactly **two code-generation
states** on the PeakPickerHiRes workload — 88,619,093 and 68,006,957
instructions of self cost, a quantum of **20,612,136** — and it flips between
them on source perturbations that have nothing to do with it. The two `main`
commits `e766311` and `8889ece` sit in different states, and the *only* source
file that differs between them is `src/cli/tools/map_normalizer.rs`, which the
picker never calls. `MSSpectrum::validate`'s self cost is identical in both, so
this is not an inlining transfer between those two symbols.

Three consequences, all of which cost this lane a measurement cycle:

1. **A per-function instruction diff is not attribution.** The same lane's
   `check_sorted` fusion showed −20.6 M in `Record::finish` and −5.3 M in
   `prepare_spectrum`, neither of which calls it. Only an ablation built on each
   arm attributes a saving to the pass it came from.
2. **Every arm must be rebuilt in the same batch.** Reusing a binary built in an
   earlier session moved a headline by 1.08 % — the reviewer rebuilt one arm
   from `git archive` and got a figure 261,553 instructions away from the lane's.
3. **A layout control is mandatory, and wall clock resolves nothing small
   here.** A control binary with two unrelated functions swapped in source
   order executes the *identical* instruction count, yet on this workload a
   plain rebuild of an identical tree is worth a couple of tenths of a second
   either way, and two `main` commits sit 0.100 s apart in the **opposite**
   direction to their instruction counts. Nothing under about 1 s is resolvable
   by wall clock on a loaded node, and about 0.2 s on a quiet one.

### Integration gates

All on `integrate/wave4-shared`, logs under the pass's scratch directory. Local
gates ran on the workstation (macOS arm64); every cargo gate ran on **kim**
through the gate script, on a detached node-local checkout of this branch.

| Gate | Result |
|---|---|
| `cargo fmt --all -- --check` (local) | exit 0 |
| every `tools/*.py` checker (local): `check_core_sdk`, `check_doc_coverage`, `check_module_cycles`, `check_schema_feature_graph`, `core_sdk_coverage`, `test_core_sdk`, `test_core_sdk_coverage`, and the ten `--check` generators (`generate_modifications`, `generate_enzymes`, `generate_ribonucleotides`, `generate_rnases`, `generate_monosaccharides`, `generate_controlled_vocabulary_reference`, `generate_cv_mapping_reference`, `generate_metabo_isotope_models`, `semantic_validator/projection.py`, `probes/ims_witness_source_oracle.py`) | **17 of 17 exit 0** |
| `check_core_sdk.py --source .reference/openms4-core-bc9cc12` (local) | exit 0; **2,092** distinct current source/registration/reference files verified at `bc9cc12` |
| YAML parse of `.github/workflows/rust.yml` (local) | parses; six jobs |
| `json.load` of every changed JSON (local) | 4 of 4 load: `SOURCE_PROVENANCE.json`, `docs/core-sdk-reviewed-apis.json`, `docs/core-sdk-coverage.json`, `docs/doc-coverage.json` |
| kim: `+1.85.0 check --locked --all-features --all-targets` | exit 0 |
| kim: `+1.85.0 test --locked --no-default-features --features mzml --test mzml …` (the changed mzML line, 12 targets) | exit 0: **324 passed, 0 failed, 0 ignored**; `tests/mzml.rs` runs for the first time on a reduced-feature line |
| kim: `+1.85.0 test --locked --no-default-features --test kernel --test spline_math --test peak_picking --test window_mower` (new) | exit 0: 24, 16, 15 and 11 = **66 passed, 0 failed** |
| kim: `+1.85.0 test --locked --no-default-features --features mzml,paramxml,parallel --test topp_threads --test topp_peak_picker_hi_res --test peak_picking_experiment` (new) | exit 0: 24, 20 and 11 = **55 passed, 0 failed, 2 ignored** (the two Linux-plus-`parallel` HPC thread tests, which this line is the first reduced-feature build to compile) |
| kim: `clippy --locked --all-features --all-targets -- -D warnings` | exit 0 |
| kim: `doc --locked --all-features --no-deps`, `RUSTDOCFLAGS=-D warnings` | exit 0 |
| kim: `test --locked --all-features --all-targets` | exit 0: **5,189 passed, 0 failed, 22 ignored** over 352 targets |
| kim: `+1.85.0 test --locked --no-default-features --all-targets` | exit 0: **3,500 passed, 0 failed, 3 ignored** |
| kim: `test --locked --all-features --doc` | exit 0: **67 + 3 = 70 doctests** |
| `check_doc_coverage.py --write` | floor **4,362/5,751 = 75.8 % to 4,400/5,770 = 76.3 %**; six modules move, three of them to 100 %: `src/format/featurexml.rs` 3/19 to 20/20, `src/format/featurexml_scaling.rs` new at 10/10, `src/format/identification_xml.rs` 0/2 to 2/2, plus `src/processing/peak_picking.rs` 23 to 27 items, `src/processing/spline/cubic.rs` 10 to 14 and `src/processing/window_mower.rs` 5/6 to 6/6 |
| `check_module_cycles.py`, then `--write` | 64 cross-module edges, 13 mutually-dependent pairs, exit 0 — and `--write` produces **no diff**. The `cli -> analysis` edge that three lanes asked the integrator to record was already recorded in the wave-3 pass, so there is **no new acyclic edge in this window** and none was invented to look like progress |

An independent adversarial audit of the pass reproduced every one of those
gates on its own detached worktree, recomputed all 125 provenance hashes (0
mismatched), confirmed both generators produce zero diff, ran the runnable
ignored tests against their recorded reasons and checked every `CPP-308`..`313`
citation in the pins. Verdict **approve with notes**: no blocker, no major, five
minors, all of them prose. All five are fixed in the follow-up commit — the C++
thread summary now carries the FileInfo/featureXML exception (33 peak threads,
not 64) and attributes the util ~5.3 figure to PeakPickerHiRes alone against a
per-cell range of 1.18 to 12.38; the work-package note no longer repeats the
picker lane's "does not use `in_thread_pool`" shorthand that this same pass
corrected elsewhere; the "all eight tools on full-size instrument data" claim is
qualified wherever it appeared (README, VALIDATION, PORTING_STATUS, CHANGELOG,
EARLY_TOPP_WORK_PACKAGES and the BENCHMARKS heading), because
`FeatureFinderCentroided` ran on the 4,000-spectrum subset; the reduced-feature
binary count is corrected from the pre-change 169 to 174; and the reader lane's
`base64-simd` request now has its row in the crate register rather than only a
work-package note.

One gate failed on its first attempt and the failure was in the invocation, not
the tree: the gate script interpolates its arguments unquoted, so
`--features "mzml paramxml parallel"` reached cargo as three words. Re-run as
`--features mzml,paramxml,parallel` it passes. The CI line keeps the quoted
form, which is correct inside a YAML `run:` step.

### Shared records changed by this pass

- **Ledger** (`docs/core-sdk-reviewed-apis.json`, then
  `core_sdk_coverage.py --write`): complete 60, evidence_requires_review 165,
  native_equivalent 90, **partial 61 to 62**, **unmapped 410 to 409**, over 786
  headers. The single status change is `FORMAT/DTAFile.h`, which enters at
  `partial` because `fix/dta-precision` made the port's writer the source's and
  produced byte-equal output over 836 MB, so there is now evidence to review; it
  is **not** `complete` (no claim about the filename overloads, the protected
  `default_ms_level_`, or whole-header completion). Eight further rows change
  scope without changing status, each for a stated reason:
  - `FORMAT/MzMLFile.h`: the reader rewrite (−47.2 % of the load path's
    instructions, no decoded value changed), the removal of the dead per-record
    validation with its written argument and mutation-checked test, the
    code-generation caveat, and the measured indentation fact behind the
    output-byte deficit.
  - `PROCESSING/CENTROIDING/PeakPickerHiRes.h`: the parallel entry points, the
    determinism evidence at five worker counts, and the wave-4 numbers. The
    wave-3 figures this row carried (38.7 s, 4,309 MiB) were a single
    unreplicated pair on a loaded build node and are now explicitly superseded.
    The `PARTIAL` clause gains the two scans that remain.
  - `MATH/MISC/CubicSpline2d.h`: `CubicSpline2dFitter`, marked as a **native
    addition**, not a port of a C++ member.
  - `KERNEL/MSSpectrum.h` and `KERNEL/MSChromatogram.h`: the wording is chosen
    so nobody reads the wave-4 change as a weakened `validate()`. No check was
    removed from `validate()`; the reader's dead **call** to it was.
  - `KERNEL/MSExperiment.h`: the spectra-only `ranges` against the source's
    combined `getMaxIntensity`, and the open `SummaryLimits::default().max_work`
    ceiling that refuses a file the port's own reader reads.
  - `FORMAT/FeatureXMLFile.h`: the size-derived ceilings, the executed C++
    `FileInfo` at both benchmark sizes, and the additive
    `identification_xml.rs` changes the other dialects share.
  - `APPLICATIONS/TOPPBase.h` (outside the registered union): the four tools
    that moved to executed full-size evidence, the picker as the sixth
    `-threads` tool, and four carried-open framework questions.

  **Validated TOPP workflows stay 8 of 124** — every ported tool already
  qualified in wave 3; what changed is the evidence behind four of them, not
  the count.
- **Provenance.** `SOURCE_PROVENANCE.json` registers **18 new oracle
  artifacts** across two directories: `map-normalizer-divergence` (ten, new in
  this window, four of them executed drivers and their logs — the evidence
  behind `CPP-308`) and `featurefinder-picked` (eight — C2, which the note in
  that key claimed was registered at W4.4 and was not). Every sha256 in the key,
  the **107 carried forward and the 18 added**, was recomputed from the file at
  this integration: **0 mismatched, 0 missing**. The other directories these
  lanes cite (`picked-chromatogram`, `mzml-reader-scale`,
  `topp-peak-picker-scale`, `baseline-filter-edges`, `lm-eigen-path`,
  `gauss-trace-fitter` including `solver-gap`, `b7-ffap-features`,
  `tool-threads`, `mzml-writer-scale-parity`, `topp-early-bundle`) were already
  registered and were re-verified unchanged rather than re-registered. No new
  in-repo manifest was produced: the six `tests/data/*_provenance.json` files
  this window changed were already listed. `tools/check_core_sdk.py` passes
  plain and with `--source .reference/openms4-core-bc9cc12`.
- **C++ issues.** `CPP-308` to `CPP-313`, each checked against the pinned source
  first (see below).
- **Crate register.** One new row, no new dependency. No lane in this window
  proposed taking a dependency, and two that could have are recorded as *not*
  taken: a SIMD base64 decoder, which the reader lane asked be written into
  [THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md) and now is, at
  *pending; measured, not adopted* (2.03x the current decoder on the hot path,
  worth about −0.9 s on a 2.3 GB input, against internal `unsafe` in the
  candidate and an `--offline` registry requirement); and any allocator change. The `mimalloc` measurement an
  earlier profiling lane reported (−1.18 s) is **withdrawn as not reproducible**
  — the library it `LD_PRELOAD`ed defines none of `malloc`/`free`/`calloc`/
  `realloc`/`posix_memalign`, so it initialised and intercepted nothing. Any
  future `mimalloc` decision needs a fresh measurement with an override build or
  a real `#[global_allocator]`.
- **Two stale integrator-owned documents, fixed here because the picker lane
  asked twice and neither pass acted.**
  [TOPP_THREADS_SUPPORT](TOPP_THREADS_SUPPORT.md) still said
  `outputs_are_byte_identical_across_thread_counts` covers "the five tools"
  (six since `3b943e4`) and stated the
  `every_executable_runs_its_body_on_the_requested_pool` contract without the
  picker's one exception; it now names the sixth tool, the `0`-workers-at-one
  exception and why the picker scopes its pool to the picking call.
  [TOPP_CLI_SUPPORT](TOPP_CLI_SUPPORT.md) still said the three wave-3a tools
  wire no pool, which has been false for `PeakPickerHiRes` since `3b943e4`.
  Both were verified against the code rather than against the lane report: the
  picker **does** call `ToolContext::in_thread_pool`, at
  `src/cli/tools/peak_picker_hi_res.rs:249`, around the picking call and not
  around the body, and skips it entirely at one worker — the lane's own
  shorthand ("does not use `in_thread_pool`") would have been wrong to copy.
- **Module graph and doc coverage.** Recorded in the gate table above: no new
  acyclic edge, and the doc-coverage floor ratchets from 75.8 % to 76.3 %.
- **Rewritten document.** [BENCHMARKS](BENCHMARKS.md).

### C++ issue candidates: what was logged and what was not

Logged `CPP-308` to `CPP-313`, each read in the pinned source before it was
written:

- `CPP-308` MapNormalizer's unguarded division (`src/MapNormalizer.cpp:93-103`
  of the pinned TOPP tree; the text was read in the two local TOPP checkouts
  `d0234cc` and `6f8eb94`, which agree line for line, and the three behaviours
  were **executed** on the pinned Release binary). Written from the lane's
  *amended* wording, not its first report: the entry deliberately does **not**
  extend to the empty-range case, where the C++ is correct and the port now
  matches it.
- `CPP-309` the `DTAFile` proton-mass asymmetry (`DTAFile.h:120` against
  `:204`). Marked **source-reviewed**, not executed: the arithmetic is verified
  against the cited lines, but a C++-only store-then-load round trip at charge
  > 1 was not run as a separate reproduction.
- `CPP-310` `DTAFile::store`'s two 15-digit rules. Executed, and the whole
  formatting chain was re-derived in the pinned source rather than taken from
  the report: `DTAFile.h:184` sets `os.precision(writtenDigits<double>())` = 15,
  `DPosition.h:412-420` routes the m/z through `precisionWrapper` to
  `NumericFormatting::appendNumeric(..., chars_format::fixed, 15)` — 15 digits
  *after the point* — while the `float` intensity takes the stream's default
  field at 15 *significant* digits.
- `CPP-311` MzMLSplitter's unresolvable `precursor/@spectrumRef`. Executed; the
  disabled flag that would have prevented it is at `MzMLSplitter.cpp:67-68`.
- `CPP-312` FeatureFinderAlgorithmPicked dividing a zero-width RT range
  (`FeatureFinderAlgorithmPicked.cpp:244-245`, consequence at `:1837-1838`).
  Executed on both builds; the Release behaviour — a silent empty feature map
  with exit 0 — is the reference, never the Debug precondition.
- `CPP-313` `FeatureXMLHandler.cpp:318`'s `min(Size(1e5), count)` reservation
  and the comment whose premise current data exceeds by an order of magnitude.

**Not logged, with the reason:**

- *The indexed-mzML `fileChecksum` placeholder.* Three separate lanes proposed
  it. It is already `CPP-049`, confirmed again here in
  `MzMLHandlerHelper.cpp:127-132`; `CPP-305` covers the `indexListOffset`
  beside it. No new entry.
- *The 32-bit time-array narrowing and the 15-significant-digit XML writer.*
  Already `CPP-306` and `CPP-307`, logged in the wave-3 pass.
- *The stale `TOPP_DTAExtractor_{1,2,3}_output.dta` reference files*, which the
  DTA lane reported as saying `120 100` where the pinned tool writes
  `120.0 100`. **Could not be checked against the pinned source**: the upstream
  TOPP test-data tree is not present in any local checkout, so only the tool's
  own output was available (which does write `120.0 100`). The candidate stays
  in the lane report until the test data can be read.
- *MapNormalizer's combined-maximum semantics.* Explicitly not a defect. The
  lane says so and this pass agrees: `updateRanges()` including chromatograms is
  deliberate and self-consistent, and `FileInfo` prints the combined, spectrum,
  per-level and chromatogram ranges separately. Whether normalising MS1 peaks
  against a chromatogram point is scientifically right is the tool maintainer's
  question; the port must match it either way.
- *The `startTimeStamp="-infinity"` warnings from both implementations.* An
  artefact of the benchmark input file, not of either program.
- *`src/cli.rs::run_failure` mapping `InvalidValue`/`InvalidRange` to exit 6
  where `TOPPBase` reaches 8.* A **port-side** decision, recorded at that
  function and mentioned in `CPP-308`'s Rust handling; not an upstream defect
  and not logged as one.

## Early TOPP bundle wave 3 integration (2026-09-15)

`integrate/wave2` (`864b295`) now carries, on top of the wave-2 tip `1c14d60`,
the two trace fitters, the wave-3a scaffold and its three tools, the
FeatureFinderAlgorithmPicked feature stage, the Boost.Regex facade, seven fix
lanes and the Levenberg-Marquardt rewrite. This shared-file pass on
`integrate/wave3-shared` records them in CI, the ledger, the provenance files,
the C++ issue log, the crate register and the documentation, and adds
[BENCHMARKS](BENCHMARKS.md). See
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-3-status).

The seventh fix lane, `fix/picked-chromatogram`, merged as `b1700de` (docs
`864b295`) while this pass was under audit, and this pass was rebased onto it.
Two consequences are recorded rather than smoothed over: the tree-wide gate
numbers below are the lead's at `4091665`, which is one merge behind that tip,
and the instrument-scale picked TIC chromatogram is now bit-identical, so every
statement about it here, in [BENCHMARKS](BENCHMARKS.md), in the ledger and in
[PORTING_STATUS](PORTING_STATUS.md) describes a closed finding. This pass's own
kim gates below were rerun on the rebased tree.

Tiers as in waves 1 and 2: tier 1 is an executed differential against a C++
oracle, tier 3 upstream class-test literals, tier 4 native derivation. Wave 3
adds a second C++ reference beside the Debug product SDK: the **Release build**
`openms4-release-bc9cc12-c19e494-174b576` (core `bc9cc12`, cli `c19e494`, topp
`174b576`, gcc 14.4, `-O3 -DNDEBUG`, no `-march`), whose identity and smoke
check are in [BENCHMARKS](BENCHMARKS.md) §1. Every count below is from the
package's approving verifier, rerun on its own detached checkout through the
gate script (Linux x86_64: spock, kim or dax); "1.96" is current stable and
"1.85" the minimum Rust.

| Package | Branch (merge) | Evidence tier | Gates rerun by the approving verifier |
|---|---|---|---|
| B4-GAUSS (fix round 3) | `2fc61e5` | Tier 3: 16 TraceFitter_test and 17 GaussTraceFitter_test sections; tier 1 against `../oracle/gauss-trace-fitter` and C2: start values, residuals, Jacobians and queries bit-identical, fits and the 1..500 budget sweep within 1e-9, statuses/`nfev`/`njev` and the budget boundaries exact; the 79-input solver-gap probe recorded, not asserted | `trace_fitter` 22 and `gauss_trace_fitter` 34 (including the `#[ignore]`d report) on 1.96 and 1.85, with and without default features; `--lib feature_finder_picked` 3; doctest 1; clippy `--all-targets` and rustdoc exit 0; the solver-gap oracle re-executed in a copy, all four result files byte-identical and every manifest hash verified; 43 in-repo and external sha256 recomputed, 0 mismatches |
| B5-EGH | `7d0c975` | Tier 3: the EGHTraceFitter_test sections; tier 1 against `../oracle/egh-trace-fitter` and C2: 80/80 start points, 515/515 parameter vectors and 11,023/11,023 functor values bit-identical when the oracle platform's own `exp`/`log` are used, fits and derived quantities within 1e-9 (measured maximum 1.98e-11) | on kim at `7d0c975`: `egh_trace_fitter` 32, `gauss_trace_fitter` 33 + 1 ignored, `trace_fitter` 22 on 1.96 and 1.85; `doc --all-features` exit 0 on both toolchains (the gate that failed 101 at `cf55e91`), with a negative control that re-injects one explicit link and fails 101; clippy exit 0; a 419-case before/after probe (44,246 lines) bit-identical except two error strings |
| B7-FFAP-FEATURES | `ba913da` | Tier 1: the feature stage replayed against `../oracle/b7-ffap-features` and, through the wrapper, against the Release build — the `FeatureFinderCentroided_1` family gives the Release build's own counts: 8 features, 30 hulls and 120 hull points in the default run, 24 seeds in the `-seeds` run and 1,054 hull points under `-debug 5`, with bit-identical m/z and hull points and rt/`score_fit`/`score_correlation` within 1e-9 (5.5e-13, 2.2e-10, 7.7e-12) | `test --locked --all-features --all-targets` 341 binaries, 4,960 passed, 0 failed; 14 feature-finding targets on 1.85, 271 passed; `--no-default-features --features mzml,paramxml,featurexml` 11 + 29; `build --no-default-features` ok; clippy and rustdoc exit 0; the verifier's own determinism harness fingerprinted every feature field, hull point, log line and abort entry at `Threads::serial()`, 2, 3 and 8 — byte-identical |
| wave-3a scaffold | `8f0bb3e` | Integrator-owned: the three tool registrations, their `[[bin]]` entries and `FileHandler::load_experiment_with_read_options` | covered by the wave-3a package gates below |
| P3-PICKER-TOOL | final commit (merge `4c2806d`) | Tier 1: the four registered workflows `TOPP_PeakPickerHiRes_1/_2/_5/_6` against the retained outputs (decoded, D6), the six parameter-failure registrations and `TOPPWRITEINI_OVERWRITE`, and the C1 oracle regressions plus `../oracle/topp-peak-picker-tool` | `--no-default-features --features mzml,paramxml --test topp_peak_picker_hi_res --lib` 17 + 243 on 1.96 and 1.85; the seven-target TOPP line 109 passed; `--all-features --all-targets` 338 ok lines, 0 failures; `+1.85.0 check --all-targets`, clippy, rustdoc exit 0; three release builds of the binary; 17 of the verifier's own product-SDK C++ reference runs |
| P4-PICKER-LOWMEM | `4293aab` (integrated at `integrate/wave7`) | Tier 1: the two registered low-memory workflows `TOPP_PeakPickerHiRes_3` and `_4` against the retained outputs (decoded, D6), and `../oracle/p4-lowmemory` — the C++ Release build at the pins in both modes over the upstream fixtures, over failing inputs of six kinds at two positions of a 110-record file, over a 110-record `FileMerger` output and a five-record reference fixture, and over the 2.3 GB `UK222.mzML`; that differential reproduced both retained files byte for byte, confirmed every divergence, established how each implementation ends a failing run and what it leaves on disc, and pinned the source's dangling-reference numbering record by record. Extended at integration: `../oracle/integ-w7` re-measured the same rule on the whole-document writer and isolated the content-versus-pointer split (CPP-172) | the lane's battery on kim, one gate at a time: `fmt --check`; clippy `--locked --all-features --all-targets -D warnings`; `+1.85.0 check --locked --all-features --all-targets`; `doc` with `RUSTDOCFLAGS=-D warnings`; the reduced slice `--no-default-features --features mzml,paramxml`; and `test --locked --all-features --all-targets` reconciled to 5,337 / 0 / 21. Re-verified by the reviewer on spock at the same head |
| A6-FILEINFO | `0510382` (integrated at `integrate/wave7`) | Tier 1: 59 cases of `../oracle/a6-fileinfo` against the Linux x86_64 Release build `openms4-release-bc9cc12-c19e494-174b576` on ibminode06, executed twice and reproduced; 38 of them compare `-out` and `-out_tsv` byte for byte and one its text, with only the lines that embed the input path normalised and three masks (the `FileInfo took` footer, the ProgressLogger `-- done [took …] --` lines, and the `std::cerr` byte dump of the one below-window case). `TOPP_FileInfo_11` (WILL_FAIL) and `_19` reproduced; `_12`'s exit code not, for a reader reason recorded in the checkpoint. Tier 4: an empty SRM chromatogram and a NaN entering either `std::sort` are refused | all eight green on dax at `0510382`, first attempt each: `fmt`; clippy `-D warnings`; `test --locked --all-features --all-targets` 5,380 / 0 / 21; the `mzml,paramxml,featurexml` slice 143 / 0 / 5; `--features mzml --test file_info_checks` 47 / 0 / 0; `--no-default-features` 3,626 / 0 / 3; `+1.85.0 check --locked --all-features --all-targets`; `doc` with `RUSTDOCFLAGS=-D warnings`. Re-run independently by the reviewer on spock with every figure reproduced |
| A5-FILEINFO-TOOL | `a4eb586` | Tier 1: `TOPP_FileInfo_1`, `_2`, `_3` and `_9` through FuzzyDiff against the retained outputs, and 17 executed product-SDK cases in `../oracle/topp-file-info-tool`, run twice | `--features mzml,paramxml,featurexml --test topp_file_info --test file_info` 33 and 58 (5 ignored) on 1.96 and 1.85; `--all-features` four targets 33/59/73/8; `--no-default-features --test file_info` 22; `+1.85.0 check --all-targets`, clippy, rustdoc exit 0; 66 + 2 doctests; a release build of the binary; the verifier reran all 40 C1 FileInfo cases, all 17 package cases and 47 further argv pairs against the C++ |
| C5-FFC-WRAPPER | `ba67aa9` | Tier 1: 25 cases against 29 executed C++ runs (C1 plus `../oracle/ffc-wrapper-c5`, each run twice and reproduced) | `--features mzml,paramxml,featurexml --test topp_feature_finder_centroided` 25 passed; the same on 1.85 and under `--all-features`; the whole suite under `--all-features`; clippy, rustdoc and `+1.85.0 check --all-targets` exit 0; a release build of all bins; the oracle re-executed from a copy, all 7 cases reproducing their exit codes and hashes |
| fix/ffc-integration | `4d53a7e` | Tier 1: the wrapper re-derived against the Release build now that B7 ports the algorithm — `FFC_1` decoded within 1e-9 relative on the fitted fields, byte-identical at `-threads` 1/2/4/8/0, and the two Debug-only cases re-expected from the Release build's exit 0 | whole suite 5,124 passed / 0 failed / 21 ignored on 1.96 **and** 1.85; `+1.85.0 --no-default-features --all-targets` 3,474 passed; 69 doctests; clippy, rustdoc, fmt exit 0; macOS arm64 26 passed / 1 ignored; six C++ Debug cases and seven Release cases rerun by the verifier; a probe merge onto the then-tip `4c2806d` green |
| crate/regex-facade | `2a29bd0` | Tier 1: a driver compiled against Boost.Regex 1.92 produces every compared answer. Full corpus 134,075 patterns / 6,531,674 compared cases / 0 answer mismatches / 0 compile or mark-count differences / 0 panics, 39,452 refusals in 44 categories; committed fixture 2,416 / 162,714 / 0, 512 refusals. `refusals.py`, an independent model built from the Boost headers, reproduces the refusal set exactly | `--test boost_regex` 27 and `--lib concept::boost_regex` 19 on 1.96 and 1.85, and on 1.85 without default features; clippy and rustdoc exit 0; the verifier added 41 family runs of its own (1,437,572 patterns, 160,773,228 cases) and `fuzz.py` at a fresh seed (80,000 expressions, 3,659,817 cases), 0 mismatches everywhere; the full corpus also clean in an overflow-checked debug-assertion build and with `MAX_AUTOMATON_ATOMS = 0`; the counter-check that the grids are not vacuous reproduces the document's 55,968 answer differences at the previous commit |
| bundle/B3b-LM-FIDELITY | `dc56a9f` | Tier 1: `../oracle/lm-eigen-path` traces every intermediate of Eigen 5.0.1's `minimizeOneStep`, `lmpar2`, `qrsolv` and `ColPivHouseholderQR` around the library's own trace functors over 141 fits. After the rewrite all 141 agree with the Linux x86_64 Release build bit for bit in every evaluation argument, the final parameters, the status, `nfev` and `njev`, and the port's fits equal that build's 141 of 141 | on spock (a third host): `--all-features --all-targets --no-fail-fast` 5,127 passed / 20 ignored / 6 failed, identical on 1.85 — the six are `tests/topp_feature_finder_centroided.rs` and were proven pre-existing on a pristine `aa0816e` worktree with byte-identical panic bodies; clippy and rustdoc exit 0; the six fitter/solver targets green without default features; native macOS arm64 runs of `lm_eigen_path_differential`; all 11 oracle artifacts re-hashed |
| fix/mzml-reader-scale | `f03ef85` | Tier 1: the C++ Release FileInfo, FileConverter and MzMLSplitter over five `startTimeStamp` sentinels and a control, plus the Debug SDK for the completion-time case; tier 4: the per-byte allowance measurement over all six benchmark inputs, with adversarial documents still refused | on kim: `test --locked --all-features` 320 ok lines, 4,940 passed on 1.96 and 1.85; fmt, clippy, rustdoc exit 0; a trial merge with the then-tip `a185bbb` conflict-free and green (326 ok lines, 5,032 passed); the six `#[ignore]`d HPC tests run one process each on ibminode06; an independent array scan of all four large inputs |
| fix/mzml-writer-scale-parity | `522ce8f` | Tier 1: the C++ Release tools' own output compared field by field after decoding (D6) on 600- and 5,000-spectrum real slices; the index and digest verified independently by `tools/mzml_writing/check_output.py` | on spock: fmt, clippy `--all-targets`, `test --locked --no-fail-fast`, the same with `--all-features`, `+1.85.0 test`, rustdoc — all exit 0; real-input comparison on ibminode06 at two slice sizes plus BaselineFilter; the full-size `#[ignore]`d writer tests; four adversarial documents refused with the destination untouched |
| fix/picker-scale | `6b55773` | Tier 1: the shipped tool end to end on the 2.3 GB benchmark input against the C++ Release tool, 22,776,198 centroids bit-identical in m/z and intensity | on spock: `+1.85.0 --no-default-features --all-targets` 343 targets / 3,476 tests (the blocker's proof: the example builds under no default features); 17 targets / 496 tests under `--all-features`; clippy and rustdoc exit 0; `--all-features --all-targets --no-fail-fast` 5,123 passed / 6 failed, the six proven pre-existing on a pristine `aa0816e` worktree; the memory cost of the replaced test measured at 118 MiB against 835 MiB |
| fix/tool-threads | `a442694` | Tier 1: `../oracle/tool-threads` — 14 executed C++ cases and 60 Rust-versus-C++ pairs, with `/proc/<pid>/task` sampled to count the threads each implementation starts; the executed C++ shows `-threads -1` and `-7` start `omp_get_num_procs()` threads, which `Threads::from_cli` had mapped to one | on spock: eight TOPP and parallel targets under `--all-features` on 1.96 and 1.85; `--no-default-features --features mzml,paramxml`; clippy, rustdoc, fmt and the doctests exit 0; the verifier's own thread probe over 49 further cases on ibminode06 |
| fix/baseline-filter-last-point | `4a597cf` | Tier 1: the C++ Release BaselineFilter over the edge shapes, an exhaustive element-length sweep, the tool and the full UK222 run (`../oracle/baseline-filter-edges`) | on kim: fmt, clippy `--all-targets`, `test --locked --all-features --all-targets`, twelve targets on 1.85, `--no-default-features`, rustdoc — all exit 0; the verifier built its own independent sweep oracle (37,730 forked cases, 1,199,520 samples) and its own mzML decoder, reproducing 119 of 119 committed tool rows; a negative control that reverts the four method arms fails |
| fix/picked-chromatogram | `e269586` (`b1700de`) | Tier 1 in both directions: `../oracle/picked-chromatogram` — the C++ Release build reads and picks seven cases built from the benchmark file's own TIC arrays (32/64-bit against minute/second, plus the empty, single-point and unsorted shapes) and every value is compared as an IEEE-754 bit pattern; the `min32` rows pin `ReadOptions::source_time_array_precision` and the `min64` rows pin the library default | on dax/spock: fmt, clippy `--all-targets`, rustdoc, `test --locked --all-features --all-targets --no-fail-fast` (5,141 passed; the only failures were the six `topp_feature_finder_centroided` tests that `fix/ffc-integration` then closed, proven pre-existing on the merge base), 21 mzML targets and 5 picking targets without default features on 1.96 and 1.85, 3 doctests; the verifier rebuilt the oracle from its own source and its own probe and got the committed TSV byte for byte (209 rows), re-ran the full 2.3 GB input through three binaries (C++, Rust before, Rust after) with its own comparators, and proved the tests non-vacuous by forcing the switch off (5 of the new tests fail) |

Lead's results for the merged tree at `4091665` (dax, detached): `test --locked
--all-features --all-targets --no-fail-fast` **5,141 passed, 0 failed, 21
ignored** on stable and on `+1.85.0`; `+1.85.0 --no-default-features
--all-targets` 3,486 passed, 3 ignored; 69 doctests; `clippy --all-targets -D
warnings` and rustdoc `-D warnings` clean.

### Instrument-scale comparison

The first Rust-against-C++ measurement on real data, and the reason
[BENCHMARKS](BENCHMARKS.md) exists. `PeakPickerHiRes` on
`profile_hr_qe_silac_uk222/UK222.mzML` (2,317,975,830 bytes, 40,856 spectra,
197,765,338 raw points), both tools driven by the INI the C++ tool writes with
`-write_ini`, `threads 1`:

| | C++ Release | Rust |
|---|---|---|
| wall | 26.80 s and 28.85 s (two runs) | 38.69 s |
| peak RSS | 3,977,240 KiB = 3,884 MiB | 4,412,352 KiB = 4,309 MiB |
| centroids | 22,776,198 | 22,776,198, **every m/z and every intensity bit-identical** |
| picked TIC chromatogram | 8,174 points | 8,174 points, **every retention time and every intensity equal** |

The chromatogram row is the one this pass had to correct. The run above showed
it differing (8,173 of 8,174 retention times by up to 3.19e-3 s, and 7,891 of
8,174 intensities by more than 1e-6 relative, worst 1.75e-3);
`fix/picked-chromatogram` then root-caused that in the **mzML reader**, not in
`pick_chromatogram` — the source narrows a unit-converted 32-bit time array
back to `f32` (`CPP-306`), which this input's TIC time array hits and the
picker's spline apex amplifies — and `ReadOptions::source_time_array_precision`
makes the two agree bit for bit over the whole file. Its verifier reproduced
both states from its own comparators on its own runs, before and after.

One difference remains, recorded rather than explained away: spectrum retention
times differ in 10,671 of 40,856 records by at most 9.09e-13 s, 1 to 2 ULP of
`f64`. That one is writer-side and the port's text is the correct side — the
C++'s own output fails to reparse to its own stored `double` in exactly those
10,671 records, because its writers print 15 significant digits (`CPP-307`);
the two readers agree bit for bit on all 40,856. The run was made by
`../oracle/topp-peak-picker-scale/bench_06.sh` on ibminode06 under foreign
load, once per implementation, so the wall times are indicative and carry no
median, IQR or confidence interval; the peak-RSS figures come from
`/usr/bin/time -v` and are not affected. BENCHMARKS §4 states the rest
of the caveats, including the harness's own peak-RSS measurement error, the
start-up floor, and that the Rust tools are serial where C++ uses OpenMP.

### The 21 ignored tests

All 21 on Linux with `--all-features` (22 `#[ignore]` attributes exist; the
macOS-only `macos_arm64_sdk_gap_report` is compiled out there). Three survive
`--no-default-features`, which matches the lead's 3.

| Test | Reason | Documented gap or external dependency | Owner |
|---|---|---|---|
| `fuzzy_string_comparator.rs::verbose_3_log_bytes_are_bounded` | fills the 256 MiB log buffer | resource cost, not a gap; the bound it checks is asserted | C3-FUZZY |
| `lm_budget_differential.rs::levenberg_marquardt_crate_candidate_gate_report` | the gate report of the rejected crate candidate; asserts nothing | [THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md), Levenberg-Marquardt row | B3-LM |
| `gauss_trace_fitter.rs::solver_gap_probe_reports_the_known_gap` | a macOS-generated oracle against a solver that now matches Linux x86_64 Release Eigen; prints a report, asserts no Rust value | [TRACE_FITTER_SUPPORT](TRACE_FITTER_SUPPORT.md) "Known gap"; [DISTRIBUTION_FITTERS_SUPPORT](DISTRIBUTION_FITTERS_SUPPORT.md) §1, the user's platform decision | B4-GAUSS / B3b |
| `lm_eigen_path_differential.rs::macos_arm64_sdk_gap_report` (macOS arm64 only) | measures the cost of matching Linux x86_64 Release on macOS arm64; asserts nothing | the same platform decision | B3b |
| `topp_feature_finder_centroided.rs::a_zero_width_retention_time_range_diverges_from_the_cpp_release_build` | documented divergence: the port refuses a zero-width RT range that the Release build carries through to an empty feature map | `CPP-274`, [TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT](TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md) native difference 2 | B10 |
| `file_info.rs::c1_file_info_9_mzml_mps`, `::a4_file_info_9_default_flags` | the strict mzML reader refuses `FileInfo_9_input.mzML` (a repeated spectrum userParam `name`, `dataProcessingRef` on the m/z and intensity arrays, a 64-bit float charge array) | [FILE_INFO_SUPPORT](FILE_INFO_SUPPORT.md) "Known reader gaps", decision D10 | mzML reader owner |
| `file_info.rs::a4_indexed_file_info_12_all_flags` | a 64-bit float `charge array` in `FileInfo_12_input.mzML` | the same | mzML reader owner |
| `file_info.rs::c1_empty_mzml_mps` | the dangling `defaultDataProcessingRef` of `empty.mzML` | D10; the option now exists and A5's tool passes it, but `file_info::Options` still defaults strict | A6 |
| `file_info.rs::a4_mzml_file_1_all_flags` | the selected-ion drift time is not copied onto the MS2 spectrum | A3 request 5; the lead decided on 2026-09-15 to follow the executed source, lane not yet opened | the lead |
| `mzml_reader_scale.rs` × 6 (`hpc_*`) | read `/ceph/ibmi/abi/oliver/bench/openms4/inputs` on the IBMI nodes | external dependency (multi-GB staged inputs) | fix/mzml-reader-scale |
| `mzml_writer_scale.rs::hpc_benchmark_centroid_uk222_picked_stores_through_the_tool_path`, `::hpc_benchmark_profile_uk222_stores_through_the_tool_path` | read the 547 MB and 2.3 GB benchmark inputs from `/ceph` | external dependency | fix/mzml-writer-scale-parity |
| `topp_baseline_filter_edges.rs::uk222_first600_matches_the_release_tool`, `::uk222_full_matches_the_release_tool` | read `/ceph/ibmi/abi/oliver` on the IBMI nodes; the second needs about 25 GB of memory | external dependency; the C++ side is a lane-private oracle directory (see the carried-forward note) | fix/baseline-filter-last-point |
| `topp_threads.rs::hpc_benchmark_slices_are_thread_invariant`, `::hpc_full_size_inputs_are_thread_invariant` (Linux + `parallel` only) | read the staged benchmark inputs under `/ceph` | external dependency | fix/tool-threads |

Every reason names a documented gap or an external dependency, and no ignore
hides an unexplained failure. The buckets, counted on the 21 Linux rows: **four
reports** meant to be read with `--ignored --nocapture` (three on Linux — the
fourth, `macos_arm64_sdk_gap_report`, is compiled out there — and one of the
four, `verbose_3_log_bytes_are_bounded`, does assert its bound and is listed
here only because it is ignored for its 256 MiB cost); **twelve** need the IBMI
`/ceph` share (`mzml_reader_scale` × 6, `mzml_writer_scale` × 2,
`topp_baseline_filter_edges` × 2, `topp_threads` × 2); **five** name a reader
gap (the `file_info` rows) and **one** a decided divergence (the
FeatureFinderCentroided zero-width RT range), each with a live tripwire or an
owner. Three plus twelve plus five plus one is the 21 that
`--all-features --all-targets -- --ignored --list` prints. The wave-2 tripwire
`reader_gaps_behind_the_ignored_cases_are_still_present` still fails when a
reader gap closes.

### CI

The important finding first, because it contradicts the brief this pass
started from. **`.github/workflows/rust.yml` did run the three new tool test
binaries.** The first step of both the `test` and the `minimum-rust` job is
`cargo test --locked --all-features --all-targets`, and none of
`topp_peak_picker_hi_res`, `topp_file_info` or `topp_feature_finder_centroided`
is gated out of it (their `#![cfg(...)]` headers ask for `mzml`, `paramxml` and
`featurexml`, all of which `--all-features` enables). That is exactly the
command B3b's verifier ran when it found the six FeatureFinderCentroided
failures on the merged tree, so CI would have gone red on them. One lane report
states the opposite in as many words ("CI does not run that binary, so CI will
not catch it") and is wrong; the `fix/baseline-filter-last-point` verifier had
already caught the same mistake in its own lane and made the fixer restate the
request. The real gap was narrower: fourteen targets had **no feature-sliced
line**, so nothing proved they build and pass outside `--all-features`.

Added in this pass, in the job's quoted space-separated feature style:

- `test` job, the `"mzml paramxml"` line: `topp_baseline_filter_edges`,
  `topp_peak_picker_hi_res`, `topp_threads`; and a new
  `"mzml paramxml featurexml"` line for `topp_file_info` and
  `topp_feature_finder_centroided`.
- `minimum-rust`: `lm_eigen_path_differential` on the existing
  no-default-features line, a new no-default-features line for `boost_regex`,
  `trace_fitter`, `gauss_trace_fitter`, `egh_trace_fitter` and
  `baseline_filter_edges`; `mzml_reader_scale` and `mzml_writer_scale` on the
  `--features mzml` line; `topp_peak_picker_hi_res`, `topp_threads` and
  `topp_baseline_filter_edges` on the `"mzml paramxml"` line; and
  `feature_finder_picked`, `topp_file_info` and `topp_feature_finder_centroided`
  on the `"mzml paramxml featurexml"` line.

Whole-job audit, mechanical over all 325 integration targets in `tests/`:
**no test binary is unrun by CI.** 169 now have an explicit `--test` line;
the other 156 run through a line with no target list. Nine of those 156 are
feature-gated and reach only the two `--all-features --all-targets` steps —
`mzml`, `mzml_auxiliary_review`, `mzml_auxiliary_writer_review`,
`mzml_param_groups`, `mzml_review`, `on_disc_experiment` (feature `mzml`) and
`network`, `network_get_request`, `update_check` (feature `network`, off by
default and not in any narrower line). No line was added for them: the mzML six
are already covered at the feature boundary by the 18-target `--features mzml`
line, and the network three would make CI depend on the internet. Recorded
here so the choice is visible rather than accidental.

### Shared records changed by this pass

- **Ledger** (`docs/core-sdk-reviewed-apis.json`, then
  `core_sdk_coverage.py --write`), from complete 56,
  evidence_requires_review 156, native_equivalent 90, partial 61,
  unmapped 423 to **60, 165, 90, 61 and 410** over 786 headers:
  - `TraceFitter.h`, `GaussTraceFitter.h` and `EGHTraceFitter.h` get real
    review entries at `complete` (every public and protected member mapped,
    tier 3 class-test sections, tier 1 against the C2 and package oracles).
    They had been pushed to `evidence_requires_review` in the last pass only by
    a manifest citation; that is now replaced by evidence. The lead's wave-2
    note to return them to `unmapped` is superseded by B4 and B5 merging.
  - `MorphologicalFilter.h` unmapped to `complete` (fix/baseline-filter-last-point).
  - `FeatureFinderAlgorithmPicked.h` stays `partial`, with the scope rewritten:
    B6 plus B7 port the whole algorithm and `run()` produces features end to
    end. It is **not** promoted to `complete`, for exactly the reason the lead's
    decision 5 keeps `SignalToNoiseEstimatorMedian.h` partial — `write_debug` is
    refused, so `writeFeatureDebugInfo_` and `abort_reasons_` are deliberately
    not ported. Promoting it is the same question, and is the lead's.
    **Superseded by the wave-5 checkpoint above (2026-09-17): `write_debug` and
    `writeFeatureDebugInfo_` are ported byte for byte, `abort_reasons_` is
    reproduced, and both headers are now `complete`.**
  - `PeakPickerHiRes.h` and `FileInfo.h` stay `partial` (P4, and A6/A7/A8) with
    their tools recorded; `MzMLFile.h` gains the reader-scale and writer-parity
    scope; the four `MATH/STATISTICS` fitter rows gain B3b's result and
    `tests/lm_eigen_path_differential.rs`.
  - After the rebase onto `864b295`, two of those scopes changed again without
    a status change: `PeakPickerHiRes.h` loses the picked-chromatogram
    divergence from its `PARTIAL` clause (it was never a picker defect) and
    states the instrument-scale memory in MiB converted from KiB by 1024, which
    the first draft of this pass had divided by 1000; `MzMLFile.h` gains
    `ReadOptions::source_time_array_precision` as the third
    source-compatibility switch, with `CPP-306` and the seven-case oracle.
  - **Validated TOPP workflows 5 → 8**: `PeakPickerHiRes`, `FileInfo` and
    `FeatureFinderCentroided` join, because their tier-1 manifests are now
    registered in `topp_package_reference_manifests`.
- **Provenance.** `SOURCE_PROVENANCE.json` registers **64 new oracle
  artifacts** across 16 oracle directories (`b7-ffap-features`, `gauss-trace-fitter`
  with its `solver-gap`, `exp-simulation` and `param-range` sub-runs,
  `egh-trace-fitter`, `lm-eigen-path`, `boost-regex`, `mzml-reader-scale`,
  `mzml-writer-scale-parity`, `baseline-filter-edges`, `topp-peak-picker-scale`,
  `picked-chromatogram`, `topp-peak-picker-tool`, `topp-file-info-tool`,
  `ffc-wrapper-c5`, `tool-threads`, `topp-early-bundle` and `release-build`;
  `picked-chromatogram`'s four came with the rebase onto `864b295`), each sha256
  recomputed from the file and equal to the value its own lane's manifest
  records. Four manifests join `current_sdk_reference_manifests`
  (`boost_regex`, `gauss_trace_fitter`, `egh_trace_fitter` and the wave-1
  omission `fuzzy_string_comparator`) and three tool manifests join
  `topp_package_reference_manifests`. `tools/check_core_sdk.py` passes plain
  and with `--source .reference/openms4-core-bc9cc12` (2,092 files verified,
  1,720 added source references).
- **C++ issues.** `CPP-289` to `CPP-307`, every one checked against the pinned
  source first (see below). `CPP-274` gains executed evidence at algorithm
  level and `CPP-265` loses a wrong citation.
- **Crate register.** The Boost.Regex row is closed (done, `2a29bd0`); the
  Levenberg-Marquardt row records B3b's result, and its two pins moved to
  `[dev-dependencies]` per the lead's decision.
- **Module graph.** One new acyclic edge, `cli -> analysis` (64 edges, 13
  mutually-dependent pairs, unchanged). **Doc coverage** floor 4,362/5,751 =
  75.8 %.
- **New document.** [BENCHMARKS](BENCHMARKS.md).

### C++ issue candidates: what was logged and what was not

Logged `CPP-289` to `CPP-307`, each verified against
`.reference/openms4-core-bc9cc12` before being written: the four B4 candidates
plus its `Param` finding (289-293), the shared `TraceFitter` documentation
defect (294), the four B5 candidates (294-297, the wording one shared), the
five B7 feature-stage candidates (298-302), the two morphology candidates
(303-304) and the writer's `indexListOffset` (305). `CPP-306` and `CPP-307` are
the two candidates `fix/picked-chromatogram` handed over, added after the
rebase onto its merge: the reader's `float&` narrowing of a converted 32-bit
time array (`MzMLHandlerHelper.cpp:217-222`, against the `double` of the
64-bit branch at `:210-216` and the forced `PRE_64` of the Numpress path at
`:185`), and the XML writers' `writtenDigits<double>()` = 15, which is
`digits10` and not `max_digits10`, so the C++'s own text fails to reparse to
its own stored `double` in 10,671 of 40,856 retention times of the benchmark
output.

Omitted, with the reason:

- **The zero-width retention-time range** the FeatureFinderCentroided lane
  raised is `CPP-274`, logged in wave 2. It is the same defect one call frame
  earlier, so it gained the lane's executed evidence — Debug dies in an
  `OPENMS_PRECONDITION` inside `ProgressLogger::init`, Release computes
  non-finite bin bounds and writes an empty feature map — instead of a new
  identifier.
- **Everything the Boost.Regex lane found** (`create_startmap`'s stale case
  flag, the `\<`/`\>` map sharing, the backstep blow-up, the fancy-regex
  optimizer rewrites and the regex-syntax prefix factoring). These are
  Boost.Regex and crate defects, not OpenMS ones, and no pinned OpenMS
  expression has any of the shapes. They are recorded in
  [BOOST_REGEX_SUPPORT](BOOST_REGEX_SUPPORT.md) and
  [THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md).
- **The C++ `fileChecksum` placeholder** is `CPP-049`, already logged; only the
  offset defect beside it is new.
- **The `-threads` semantics** (`-threads -1` starting `omp_get_num_procs()`
  threads) is the source behaving as documented in `TOPPBase.cpp:92-95`; the
  port's own `Threads::from_cli` was wrong, which is a port fix, not a C++
  issue.

### Integration gates

Run on kim through `~/.local/bin/openms-kim-gate.sh` (slot `integ-w3-shared`),
one at a time, on the Rust tree of this pass; logs under the session
scratchpad. Local checks ran on macOS arm64.

| Gate | Result |
|---|---|
| `+1.85.0 check --locked --all-features --all-targets` | exit 0 |
| `+1.85.0 test --locked --no-default-features --test boost_regex --test trace_fitter --test gauss_trace_fitter --test egh_trace_fitter --test baseline_filter_edges` (the new line) | exit 0: 27, 22, 33 + 1 ignored, 32 and 5 passed |
| `+1.85.0 test --locked --no-default-features` with the eight targets of the amended helper line | exit 0: 26, 16, 6, 30 + 1 ignored, 35, 7 + 1 ignored, 4 and 33 passed (`lm_eigen_path_differential` is the addition) |
| `+1.85.0 test --locked --no-default-features --features mzml` with the 18 targets of the amended mzML line | exit 0: 240 passed, 0 failed, 8 ignored (`mzml_reader_scale` and `mzml_writer_scale` are the additions, 6 + 2 of the ignores) |
| `+1.85.0 test --locked --no-default-features --features mzml,paramxml --test topp_cli_lifecycle --test peak_picking_experiment --test topp_peak_picker_hi_res --test topp_threads --test topp_baseline_filter_edges` | exit 0: 73, 19, 18, 7 and 4 + 2 ignored passed |
| `+1.85.0 test --locked --no-default-features --features mzml,paramxml,featurexml --test feature_finder_picked_seeds --test feature_finder_picked --test topp_file_info --test topp_feature_finder_centroided` | exit 0: 29, 11, 33 and 26 + 1 ignored passed |
| `clippy --locked --all-features --all-targets -- -D warnings` | exit 0, no warning or error line |
| `RUSTDOCFLAGS='-D warnings' doc --locked --all-features --no-deps` | exit 0, no warning or error line |
| `test --locked --all-features --all-targets -- --ignored --list` | exit 0: exactly 21 names, one per row of the table above |

Every row was rerun on the rebased tree (this pass's commit on `864b295`), so
the counts include `fix/picked-chromatogram`'s five new
`peak_picking_experiment` tests — that target is 19 here where the first run of
this pass recorded 14, and it is the only count the rebase moved. The 21
ignored names are unchanged by it.

The two changed `test`-job lines are the same commands with the stable
toolchain and the same target lists, so they are covered by the two
`mzml,paramxml` and `mzml,paramxml,featurexml` rows above.

Local: `cargo fmt --all -- --check`; the `--check` generators — the nine that
need no argument, including the two that run a probe
(`generate_metabo_isotope_models.py`, `tools/probes/ims_witness_source_oracle.py`),
and four more given the pinned source root
(`generate_mzml_typed_reference.py`, `generate_mzml_validator_reference.py`,
`generate_proforma_spectra_reference.py`, `generate_xlms_reference.py`), all
thirteen exit 0; the two that cannot run here are
`generate_metabo_predictor_reference.py`, which takes a prebuilt probe binary,
and `generate_proforma_conversion_mass_probe.py`, which compiles a C++ probe
and so hits the same Xcode licence gate as the feature-graph checker below;
`check_doc_coverage.py`,
`check_module_cycles.py`, `check_core_sdk.py` (also with
`--source .reference/openms4-core-bc9cc12`), `test_core_sdk.py`,
`core_sdk_coverage.py` and `test_core_sdk_coverage.py` with no argument, which
is their check mode; a YAML parse of `rust.yml`; and `json.load` of every
changed JSON file.

`check_schema_feature_graph.py` **was not run whole in this pass**, and this is
the one checker of the brief's "every `tools/*.py` checker" that is not covered
above. Run locally it exits 1 on its first feature selection, in `cargo check`,
with `linking with cc failed: exit status: 69 … You have not agreed to the
Xcode license agreements`; `~/.local/bin/openms-kim-gate.sh` executes only
`cargo $*` on the remote host, so it cannot drive a Python checker there
either. Its other half was run: the `cargo tree` assertions pass on all four
feature selections (default, `--no-default-features`, `+sqlite`, `+sqmass`) —
no `libxml`/`bindgen`/`clang-sys` in any of them, no `rusqlite`/`libsqlite3-sys`
outside the SQLite ones and both present inside them — and the same run
confirms that the `[dev-dependencies]` move keeps `levenberg-marquardt` and
`nalgebra` out of the normal build graph in all four. CI covers the whole
script: the `portable-feature-graph` job runs it on ubuntu on every push.
Giving the gate script a way to run a repository checker on the gate host would
close the hole; **the lead**.

## Early TOPP bundle wave 2 integration (2026-09-14)

`integrate/wave2` (`1c14d60`) merges eight verifier-approved package branches
with `--no-ff` onto the wave-2 scaffold `42c21e7`, which sits on `1d80eed`. The
shared-file pass on `integrate/wave2-shared` then records them in the module
graph, the ledger, the provenance files, the minimum-Rust CI job, the crate
register, the licences and the C++ issue log. Not included: B4-GAUSS (fix round
3), B5-EGH (waits for B4), lane B3b (the transcription's own first-step
divergence from Eigen, under investigation), the wave-3a packages P3, A5 and C5,
and the Boost.Regex facade. See
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-2-status).

Tiers as in wave 1. The counts below are from each package's approving verifier,
rerun on its own detached checkout through the gate script (spock for first
reviews, kim for the re-reviews of A4's and B9's fix rounds; Linux x86_64);
"1.96" is current stable and "1.85" the minimum Rust.

| Package | Branch (merge) | Evidence tier | Gates rerun by the approving verifier |
|---|---|---|---|
| B3-LM (partial: D2 fallback) | `f604b4a` (`14284fb`) | Tier 1: the kept transcription against C2 at 29,004 trace-fit budgets (58 fits at max_fev 1..500, 4 degenerate fits at 500), status, nfev and njev exact, x within 1e-9 relative with a 1e-12 floor (largest 6.36e-10); tier 4: Eigen's counting rule on the 16 distribution-fitter fits; the `levenberg-marquardt` candidate is a measurement only (ignored test) and failed the gate | `--no-default-features` and `--all-features` `lm_budget_differential`, `math_distribution_fitters`, `posterior_error_probability` 7 (1 ignored), 30 and 9 on 1.96; `--lib math::` 31; 1.85 `--no-default-features --lib` 225 plus the three targets; clippy `--all-targets` and rustdoc exit 0; the ignored gate report in debug, release and `minpack-compat` reproduced 29,003 of 29,004 and 3,811 failures; the fixture regenerated byte for byte; thread-vs-serial 0 mismatches |
| P1-PICKER-LIB | `1566a5f` (`6386e10`) | Tier 1: 93 product-SDK cases, 8,795 centroids, noise ratios and percentages bit for bit in the native default and `PickingCompatibility::source()`; defaults equal the SDK's stored `getDefaults`; retained class-test outputs; tier 4 refusals and limits | `--all-features --lib` plus 20 targets on 1.96 (lib 282, `peak_picking` 15, `peak_picking_experiment` 10), the same without three targets on 1.85; 1.85 `--no-default-features --features mzml,paramxml` 15 and 10; clippy and rustdoc exit 0; the oracle driver re-executed byte-identically; 26 of 28 new adversarial C++ cases bit-identical (the 2 others: a zero intensity total, open note); three mutation checks |
| P2-MZML-LENIENCY | `15a5f09` (`fabaea4`) | Tier 1: `../oracle/p2-mzml-leniency`, 68 lines on load, metadata-only and transform; tier 4: strict messages, scope, warnings, the padded-placeholder panic fix | `--all-features --lib` plus 47 binaries 1,141 passed on 1.96 and 1.85; 1.85 `--no-default-features --features mzml --lib` plus 6 targets 327; clippy and rustdoc exit 0; the oracle re-executed byte-identically and on 5 adversarial cases (3 match, 2 give the documented-gap findings) |
| CLI-2 | `f886d90` (`b4f4456`) | Tier 1: `../oracle/topp-cli-lifecycle` (34 CLI-2 cases run twice, 15 INI read-failure cases), all five `-write_ini` files, 10 help texts and the retained processing records; tier 3: the `-write_ini` part of `TOPPBase_test`; tier 4 run-phase mappings | `--no-default-features --features mzml,paramxml` 7 targets 108 passed (lifecycle 73) on 1.96 and 1.85, socket and terminal cases executed; 1.85 `check --all-features --all-targets`, clippy and rustdoc exit 0; 13 product-SDK probes |
| A4-FILEINFO-CORE (fix round) | `bb28715` (`93bd214`) | Tier 1: 34 text and 34 TSV reports byte for byte against the product-SDK FileInfo (C1 and `../oracle/file-info-core`), TOPP_FileInfo_1/2/3/9 through FuzzyDiff; tier 3: 7 of 9 FileInfo_test sections; tier 4 refusals | on kim: `--all-features --test file_info` 59 passed, 5 ignored on 1.96 and 1.85; `--no-default-features` 22 on both; `mzml,featurexml` 58 (5 ignored) on both; `mzml` 42 (5 ignored); `featurexml` 38; doctests 8; clippy and rustdoc exit 0; `-- --ignored` fails all 5 exactly as their reasons say; 12 new SRM reports against the SDK byte-identical |
| B6-FFAP-SEEDS | `80bbdf1` (`a532c31`) | Tier 1: settings, 21 quantiles, isotope windows, every per-peak score and the seed lists in 7 configurations (C2 and `../oracle/b6-ffap-seeds`) bit for bit, except 99 of 30,840 overall scores that are the correctly rounded value one binary32 step from the Apple `powf` oracle; the seed lists come from replicas of the source selection (adapted) | `feature_finder_picked_seeds`, `feature_finder_picked_helper_structs`, `isotopes_source_precision`, `geometry_bounding_box` 29, 26, 16 and 6 on 1.96 and 1.85; 1.85 `--no-default-features --features mzml,paramxml,featurexml` 29; clippy and rustdoc exit 0; 29 on local macOS arm64; a reviewer C++ probe with 61,127 per-peak scores identical in six configurations apart from 60 correctly rounded overall scores |
| B8-IMSPLIT | `c527c5f` (`09ce3d4`) | Tier 3: the constructor, destructor and three splitByFAIMSCV sections; tier 1: `../oracle/im-data-converter` (26 cases, 369 input records, 41 groups, log lines through a LogStream emulation) and the C2 split records; tier 4 | `--test im_data_converter --test faims_helper` 19 and 15 with all features and with `mzml` on 1.96 and 1.85; `--no-default-features` 16 on both; doctest 1; 1.85 `check --all-features --all-targets`, clippy and rustdoc exit 0; a randomized differential of 20,000 experiments against a transliteration of the source |
| B9-OVERLAP (fix round) | `f1c7785` (`1c14d60`) | Tier 3: 14 of 14 FeatureOverlapFilter_test sections; tier 1: 78 product-SDK cases and 16,917 callbacks replayed bit for bit, the pinned quadtree header executed; tier 2: a Release replica for the 4 Debug-only aborts; tier 4 atomicity | on kim: `--lib --test feature_overlap_filter` lib 285 and 33 with all features, 228 and 33 without default features, 251 and 33 with defaults, on 1.96, and 285/33 and 228/33 on 1.85; 1.85 `check`, clippy and rustdoc exit 0; 400 random maps each for the journal and snapshot paths equal; peak RSS 27.6 MB at 64,000 co-located features (2.35 GB at 16,000 before the fix); four mutations caught |

Lead's kim results for `integrate/wave2` at `1c14d60`: `+1.85.0 check --locked
--all-features --all-targets` exit 0; `test --locked --all-features
--all-targets --no-fail-fast` on dax (1.96): 334 result lines, 4,858 passed,
0 failed, 7 ignored. At `1d80eed` there was 1 ignored test (the 256 MiB
log-buffer case in `fuzzy_string_comparator`); the six new ignores are listed
below.

### New ignored tests

| Test | Reason | Owner |
|---|---|---|
| `tests/lm_budget_differential.rs` `levenberg_marquardt_crate_candidate_gate_report` | the gate report of the rejected crate candidate; it asserts nothing and runs with `--ignored --nocapture` | B3-LM, lane B3b |
| `tests/file_info.rs` `c1_file_info_9_mzml_mps`, `a4_file_info_9_default_flags` | the strict mzML reader refuses `FileInfo_9_input.mzML`: a repeated spectrum userParam `name`, `dataProcessingRef` on the m/z and intensity arrays, a 64-bit float charge array (C++ loads it) | mzML reader owner (held by no wave-2 package), decision D10 |
| `tests/file_info.rs` `c1_empty_mzml_mps` | the dangling `defaultDataProcessingRef` of `empty.mzML`; P2's `source_dangling_references` now exists, but FileHandler and `file_info::Options` do not pass it | integrator (FileHandler call site) and A5 (the Options field) |
| `tests/file_info.rs` `a4_indexed_file_info_12_all_flags` | a 64-bit float `charge array` in `FileInfo_12_input.mzML` | mzML reader owner, D10 |
| `tests/file_info.rs` `a4_mzml_file_1_all_flags` | the selected-ion drift time is not copied onto the MS2 spectrum (A3 request 5) | the lead |

Each reason names a documented gap: `docs/FILE_INFO_SUPPORT.md` "Known reader
gaps" and `docs/DISTRIBUTION_FITTERS_SUPPORT.md` section 8. The active tripwire
`reader_gaps_behind_the_ignored_cases_are_still_present` fails once a reader
gap closes, and A4's re-verifier ran the five with `--ignored`: each fails for
its stated reason. No ignore hides an unexplained failure.

### Shared records changed by this pass

- **Module graph.** Five acyclic edges recorded (63 edges, 13 mutual pairs):
  `analysis -> math` and `analysis -> param` (B6), `format -> math` (A4),
  `processing -> param` (P1), `processing -> concept` (B9). B3, P2, CLI-2 and
  B8 add none.
- **Ledger**, from complete 55, evidence_requires_review 158,
  native_equivalent 90, partial 55, unmapped 428 to 56, 156, 90, 61 and 423
  over 786 headers:
  - `FeatureOverlapFilter.h` unmapped to complete (B9: every public member,
    14/14 sections, 78 executed cases); `extern/Quadtree` recorded complete
    outside the registered union. FeatureFinderCentroided, FeatureFinderMetabo
    and FeatureFinderMultiplex lose it from their open SDK headers.
  - `FeatureFinderAlgorithmPicked.h` unmapped to partial (B6, front half only)
    and `FileInfo.h` unmapped to partial (A4 ports its model, `run`, `toText`,
    `toTSV` and the peak and featureXML branches).
  - `PeakPickerHiRes.h`, `SignalToNoiseEstimatorMedian.h`,
    `SignalToNoiseEstimator.h` (P1) and `IMDataConverter.h` (B8)
    evidence_requires_review to partial.
  - `GaussTraceFitter.h` and `EGHTraceFitter.h` unmapped to
    evidence_requires_review, only because B3's budget manifest cites their
    functors; nothing of either header is ported, and B4 and B5 own them.
  - `MzMLFile.h` (P2's option), `ParamXMLFile.h` (CLI-2's writer option) and the
    TOPPBase and ParameterInformation entries outside the union (CLI-2) are
    amended without a status change; `TraceFitter.h`, `FileHandler.h`,
    `MzMLHandler.h`, `SplineBisection.h` and `Constants.h` gain reference
    manifests; `ProcessingStep.h` gains `file_info/model.rs` as a name-collision
    candidate. Validated TOPP workflows stay at 5 of 124.
- **Provenance.** `SOURCE_PROVENANCE.json` lists the P1, A4, B8 and B9 group
  manifests under `current_sdk_reference_manifests`, registers 29 oracle
  artifacts of the eight packages and updates two CLI-1 rows for CLI-2, each
  sha256 recomputed and equal to its group manifest; B3's, P2's and B6's manifests are
  indexed through `origin_key` only. `distribution_fitters_provenance.json`
  gains the budget differential; `topp_cli_provenance.json` gains CLI-2's
  tools, tests and test definitions. `tools/check_core_sdk.py` passes plain and
  with `--source .reference/openms4-core-bc9cc12` (2,083 files verified).
- **CI.** The minimum-Rust job runs `mzml_header_leniency`,
  `lm_budget_differential`, `feature_overlap_filter`, `im_data_converter`,
  `file_info`, `peak_picking_experiment` and a new
  `--features "mzml paramxml featurexml"` line for
  `feature_finder_picked_seeds`.
- **Crate register.** The Levenberg-Marquardt row is reopened (measured, not
  adopted); the quadtree and the overall-score power are recorded as kept code.
- **C++ issues.** CPP-256 to CPP-288; CPP-027 gains executed evidence and
  CPP-247 B6's handling.
- **Licences.** The MIT notice of `extern/Quadtree` for `quadtree.rs`.

### Integration gates

Run on kim through `~/.local/bin/openms-kim-gate.sh` (slot `integ-w2-shared`),
one at a time, on the Rust tree of this pass (only documentation and JSON
changed after them); none returned 255:

| Gate | Result |
|---|---|
| `+1.85.0 check --locked --all-features --all-targets` | exit 0 |
| `+1.85.0 test --locked --no-default-features --features mzml` with the 16 targets of the changed mzML line | 224 passed, 0 failed (`mzml_header_leniency` 7) |
| `+1.85.0 test --locked --no-default-features --test feature_finder_picked_helper_structs --test isotopes_source_precision --test geometry_bounding_box --test fuzzy_string_comparator --test file_info_text_format --test lm_budget_differential --test feature_overlap_filter` | 26, 16, 6, 30 (1 ignored), 35, 7 (1 ignored) and 33 passed |
| `+1.85.0 test --locked --no-default-features --features mzml --test peak_type_estimator --test faims_helper --test mzml_mobility --test im_data_converter` | 11, 15, 12 and 19 passed |
| `+1.85.0 test --locked --no-default-features --features mzml,featurexml --test file_handler_type_detection --test file_info` | 8 passed; 58 passed, 5 ignored |
| `+1.85.0 test --locked --no-default-features --features mzml,paramxml --test topp_cli_lifecycle --test peak_picking_experiment` | 73 and 10 passed |
| `+1.85.0 test --locked --no-default-features --features mzml,paramxml,featurexml --test feature_finder_picked_seeds` | 29 passed |
| `clippy --locked --all-features --all-targets -- -D warnings` | exit 0, no warnings |
| `RUSTDOCFLAGS='-D warnings' doc --locked --all-features --no-deps` | exit 0, no warnings |

Local: `cargo fmt --all -- --check`, every `tools/*.py` checker in check mode,
a YAML parse of `rust.yml` and `json.load` of every changed JSON file pass.
Correction added by the wave-3 pass: "every checker" there did not include
the two generators that need a probe, `generate_metabo_isotope_models.py` and
`tools/probes/ims_witness_source_oracle.py`, which were not run in the wave-2
pass. Both were run at wave 3 and both pass (999 support vectors checked, and
the two source-recurrence rows with their independent exact compositions).

## Early TOPP bundle wave 1 and crate wave 1 integration (2026-09-14)

`integrate/wave1` merges nine verifier-approved branches onto main `3e171b2`.
The shared-file pass on `integrate/wave1-shared` then records them in the
ledger, provenance, CI workflow, licences and C++ issue log. `crate/digamma`
(`d94274c`) was dropped, and the unused `special` crate went with it; see
[THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md). The four
fix-round merges that followed, A2's first merge among them, are recorded below.
Not included: `crate/regex-facade`, still in review.

Tier 1 is an executed differential against a C++ oracle: the product SDK
(Debug, core `4fdec46`, accepted under decision D7) or, for the crate lanes, the
Boost or libc++ code the crate replaces. Tier 3 is upstream class-test literals;
tier 4 is native derivation. Every count below is from the package's approving
verifier, rerun on its own detached checkout through the kim gate (Linux
x86_64); "1.96" is current stable and "1.85" the minimum Rust.

| Package | Branch (merge) | Evidence tier | Gates rerun by the approving verifier |
|---|---|---|---|
| crate/quantile | `29bfbc4` (`823c8e9`) | Tier 1: the 40 reference quantiles against Boost 1.92 and a 100-digit reference; 97,489 grid points against Boost's double quantile; lfdr fixtures unchanged | `--all-features --lib --test multiple_testing --test comparison_scorers` 277/48/19 passed on 1.96 and 1.85; default features 245/48/19; 1.85 `--no-default-features --lib` 223; clippy `--all-targets -D warnings` and rustdoc `-D warnings` exit 0 |
| crate/mt64 | `1d0e9c8` (`f913cf5`) | Tier 1: Boost.Random `mt19937_64` and `uniform_int` plus libc++ `std::mt19937_64`, 19,649 output lines byte-identical | nine targets 159 passed on 1.96 and 1.85; `--lib` 274; 1.85 `--no-default-features` targets 83; clippy `--all-targets` and rustdoc exit 0 |
| crate/existing-crates | `be2bb15` (`db19767`) | Tier 4: integer and byte code only; 14,141,319 reference names, 655,360,000 calendar combinations and 9M URLs compared with the replaced code | nextest `--all-features --lib --tests` 4443/4443 on 1.96 and 1.85; clippy `--lib --tests` exit 0 (also no-default and `network`); rustdoc exit 0; doctests 55+2 on 1.96 and 57 on 1.85; 1.85 slices no-default, `idxml`, `featurexml,consensusxml` and `network` passed |
| C3-FUZZY | `30d82db` (`72650a5`) | Tier 1: 137 comparator cases against product-SDK `libOpenMSTestFramework.a` (one documented hexadecimal divergence) and FuzzyDiff runs; class-test literals tier 3 | `--test fuzzy_string_comparator` 33 passed, 1 ignored on 1.96 and 1.85; `--no-default-features` 30 passed, 1 ignored on both; `-- --ignored` 1 passed; clippy `--lib --tests` and `--all-targets` exit 0; fmt and rustdoc exit 0 |
| A3-FORMAT-IO | `ff44438` (`c8b0141`) | Tier 1: `../oracle/a3-format-io`, re-executed byte for byte; class-test literals and DTA2D/DTA/MGF routing tier 3; round trips tier 4 | full `--all-features --no-fail-fast` 4514 passed on 1.96; on 1.85 4513 passed and 1 failed (`ms_data_writing_consumer`: concurrent runs raced on a fixed temporary path; reproduced and fixed, see below; 28/28 in three reruns); `--no-default-features --features mzml` 19 on both; `mzml,featurexml` 8 on both; clippy `--all-targets` and rustdoc clean |
| B1-HELPERS | `1bd8686` (`09575c3`) | Tier 1 for the 13 computational members (164 oracle rows, bitwise); tier 4 for the accessors | `--test feature_finder_picked_helper_structs` 26 passed on 1.96 (default, no-default and all features) and 1.85 (default and no-default); clippy `--all-targets`, rustdoc and fmt exit 0 |
| A1-PTE-FAIMS | `cd7b83e` (`3d88650`) | Tier 1: `../oracle/pte-faims-helper` (4024 estimator, 1503 voltage and 1500 filter rows); class-test literals tier 3 | `--all-features --test peak_type_estimator --test faims_helper` 10 and 13 passed (1 ignored) on 1.96 and 1.85; `--no-default-features` 8 and 12 on both; clippy `--all-targets` and rustdoc exit 0; doctests 2 on both; the CI command `--all-features --all-targets` 4460 passed, 1 ignored |
| B2-ISO-GEOM | `799f610` (`3982561`) | Tier 1: `../oracle/b2-iso-source-precision` and 200 repeated runs in `../oracle/b2-iso-element-order` (the verifier added 400 more); class-test literals tier 3; work counts tier 4 | `--all-features` lib and 24 targets 645 passed on 1.96 and 1.85; `--no-default-features --test isotopes_source_precision --test geometry_bounding_box` 16 and 6 on both; clippy `--all-targets` and rustdoc exit 0; full 1.96 suite 4517 passed |
| CLI-1 | `5ae7826` (`f6bdd99`) | Tier 1 exit codes and diagnostics: `../oracle/topp-cli-lifecycle` (38 cases plus 6 INI read failures); TOPPBase_test literals tier 3; run-phase mappings tier 4 | six TOPP targets `--all-features --no-fail-fast` passed on 1.96 and 1.85 (lifecycle 62, 0 ignored); `-- --ignored` ran 0 tests on both; `--no-default-features --features mzml,paramxml` 62, 5 and 4 on both; clippy `--all-targets` and rustdoc exit 0 |

Shared records changed by this pass: the ledger (PeakTypeEstimator.h complete;
FAIMSHelper.h, FeatureFinderAlgorithmPickedHelperStructs.h, DBoundingBox.h,
CoarseIsotopePatternGenerator.h and IsotopeDistribution.h partial; TOPPBase.h
and ParameterInformation.h recorded outside the registered union), the module
ratchet (concept to chemistry removed, cli to concept and cli to metadata
added), `SOURCE_PROVENANCE.json`, the minimum-Rust CI job, `LICENSES.md` and
CPP-230 to CPP-252. The open follow-ups are listed in
[the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-1-status).

### Fix-round merges

Four `--no-ff` merges of approved fix rounds follow the audit repair `b0ac76a`
on `integrate/wave1`. The counts are from each approving verifier's rerun on its
own detached checkout through the kim gate.

| Package | Fix commit (merge) | What changed | Gates rerun by the approving verifier |
|---|---|---|---|
| B2-ISO-GEOM | `42f142f` (`05ca267`) | Documentation and provenance strings only. `ProbabilityPrecision::SourceSingle` reproduces the SDK bit for bit only in runs that iterate elements in the port's order, the majority order for natural elements, and never for a formula containing iridium, which `ElementDB.cpp:512` builds from rhenium's tables (CPP-249) | rustdoc `-D warnings` exit 0; `--all-features --test isotopes_source_precision` 16 passed on 1.96; fmt and `check_core_sdk` exit 0. Executed differential: the SDK's 84 element tables equal the port's except Ir, and 83 of 3,486 element pairs differ, every one containing Ir. 1.85 not rerun: comments and JSON strings only |
| A1-PTE-FAIMS | `f3c29cb` (`cab8976`) | `filter_peptides_by_faims_cv` refuses only a NaN target; infinite targets are filtered as C++ filters them, and the oracle gained 7 records (137 in all). The `estimate_type_with_limits` parity statement is limited to spectra without a stored type or data-processing records | `--all-features --test peak_type_estimator --test faims_helper` 11 and 14 passed (1 ignored) on 1.96 and 1.85; `+1.85.0 --no-default-features` 9 and 13; clippy `--all-targets` and rustdoc exit 0; the oracle rebuilt from a copy reran with stdout byte-identical and stderr identical after path normalisation; a C++ boundary probe of 60 target and tolerance cases matched the port on all 50 non-NaN cases |
| CLI-1 | `86d2733` (`9612872`) | A `-ini` that is neither a regular file nor a directory skips the readability precheck: `/dev/null` exits 3 and a mode-000 FIFO exits 2, before a run and with `-write_ini` (`ini_read_failures`, 10 oracle cases) | six TOPP targets `--all-features --no-fail-fast` on 1.96 and 1.85: lifecycle 65, BaselineFilter 2, DTAExtractor 5, MapNormalizer 2, MzMLSplitter 4, SpectraFilterWindowMower 4, 0 ignored; `+1.85.0 --no-default-features --features mzml,paramxml` lifecycle 65; clippy and rustdoc exit 0; the oracle driver rerun reproduced all 10 exit codes; three mutation runs |
| A2-TEXTFMT | `6cec951` (`7f0a8e0`) | The package's first merge, after two fix rounds: `StringUtils::number`, full-precision `toStr` with `std::to_chars` ties to even, integer and string `toStr`, vector `operator<<` and stream `%g` text for FileInfo, with the platform differences documented (Apple libc `%g` ties, NaN sign, the INT_MAX-byte `%f` band) | `--test file_info_text_format` 35 passed with default, all and no default features on 1.96 and 1.85; doctests `format::file_info::text_format` 6 on both; clippy `--all-targets` and rustdoc exit 0; an independent libOpenMS driver over 3,224,411 doubles with 0 mismatches; all 2^32 floats hashed identically on macOS and Linux; a 43,220-row `%g` probe on Apple libc and glibc |

The final shared-file pass records these merges. `FAIMSHelper.h` is `complete`;
`StringUtils.h` and `ListUtilsIO.h` are `partial` with FileInfo scope, and the
`Types.h` review maps `writtenDigits<float>` and `writtenDigits<double>`.
`FileInfo.h` is `unmapped` again, because A1's and A2's citations of
`FileInfo.cpp` are now `context_sources`. `SOURCE_PROVENANCE.json` registers the
A1, CLI-1 and A2 oracle artifacts, the CLI-1 group manifest moved to
`tests/data/`, the minimum-Rust job runs `--test file_info_text_format`,
CPP-253 to CPP-255 are new and CPP-245 is extended, and
[THIRD_PARTY_CRATE_DECISIONS](THIRD_PARTY_CRATE_DECISIONS.md) records a pending
candidate for a readability query that never opens the file.

The ledger delta of this pass, from the committed `docs/core-sdk-coverage.json`
(complete 54, evidence_requires_review 160, native_equivalent 90, partial 54,
unmapped 428) to 55, 157, 90, 55 and 429 over 786 headers:
- A2's merged manifest adds `tests/data/file_info_text_format_provenance.json`
  to the reference manifests of `Types.h`, `StringUtils.h` and `ListUtilsIO.h`.
- Respelling the `FileInfo.cpp` citations moves `FileInfo.h` from
  evidence_requires_review to unmapped, with no reference manifest left.
- The A2 reviews move `StringUtils.h` and `ListUtilsIO.h` from
  evidence_requires_review to partial and amend the `Types.h` review (rust,
  tests, scope), whose status stays partial.
- Moving the CLI-1 group manifest changes no row: it cites no `src/openms/`
  path.
- Promoting `FAIMSHelper.h` from partial to complete also removes it from the
  open SDK headers of FeatureFinderMetabo, FeatureFinderMultiplex and
  IonMobilityBinning. Validated TOPP workflows stay at 5 of 124.

### FAIMSHelper_test section 3 through the Rust reader

`get_compensation_voltages_section_through_the_mzml_reader` in
`tests/faims_helper.rs` waited for A3's spectrum-level `MS:1001581` read. With
A3 merged its `#[ignore]` is removed and no assertion changed. Through the kim
gate (slot `integ-final2`), `test --locked --no-default-features --features mzml
--test faims_helper`, the minimum-Rust CI line, and `--all-features` each pass
15 tests with 0 ignored, on 1.96 and on 1.85.0, and the test runs by name in all
four. Every `FAIMSHelper_test.cpp` section now passes, and the infinite and NaN
target contract of `f3c29cb` rests on the executed oracle, so `FAIMSHelper.h`
is `complete`.

### Unique temporary directories in tests and doctests

The audit repair fixed the one fixed-name path that raced, in
`tests/ms_data_writing_consumer.rs`, and listed eight more:
`tests/mascot_generic.rs:293` and `:1537`, `tests/sv_out_stream.rs:434`,
`tests/mztab_m.rs:958`, `:2425` and `:2531`, and the doctests at
`src/format/imzml_file.rs:580` and `src/format/imzml_writer.rs:890`. Each now
creates its own `TempDir::new_in(std::env::temp_dir(), false)`, removed on drop,
and no assertion changed. A re-grep of `src`, `tests`, `examples` and `benches`
finds no other written path under a fixed name. Kim gate results (slot
`integ-final2`, cargo 1.96.0 and 1.85.0):

| Gate | 1.96 | 1.85.0 |
|---|---|---|
| `test --locked --no-default-features --test mascot_generic --test sv_out_stream --test mztab_m` (the minimum-Rust CI line's feature set) | 46, 28 and 48 passed, 0 ignored | 46, 28 and 48 passed, 0 ignored |
| the same with `--all-features` | 46, 28 and 48 passed | 46, 28 and 48 passed |
| `test --locked --all-features --doc format::imzml_` | 3 passed, among them `imzml_file::ImzMLFile` and `imzml_writer::store` | 3 passed, the same doctests |

The other gates of this pass, on the same working tree: `+1.85.0 check
--locked --all-features --all-targets` exit 0; the changed minimum-Rust line
`+1.85.0 test --locked --no-default-features --test
feature_finder_picked_helper_structs --test isotopes_source_precision --test
geometry_bounding_box --test fuzzy_string_comparator --test
file_info_text_format` 26, 16, 6, 30 (1 ignored, the 256 MiB log-buffer case)
and 35 passed; `clippy --locked --all-features --all-targets -- -D warnings`
and `RUSTDOCFLAGS='-D warnings' doc --locked --all-features --no-deps` exit 0
with no warnings.

### The ms_data_writing_consumer failure under 1.85

The A3 verifier's one failure was `the_filename_constructor_creates_and_truncates`,
which panicked at `tests/ms_data_writing_consumer.rs:377:36` with `Io(NotFound)`
from `mzml::load`, right after a successful `finish()`. The cause is a race in
the test, not the toolchain and not the writer:
- The test wrote to the fixed directory
  `std::env::temp_dir()/openms_ms_data_writing_consumer_file` and removed it at
  the end.
- `TMPDIR` is unset on kim, so every agent slot shares `/tmp`.
- A concurrent run of the same binary could remove `out.mzML` between
  `finish()` and `load()`. `MSDataWritingConsumer::create` opens the path with
  `File::create`, so writing through the open handle still succeeds and only
  `load` fails.

Reproduced on kim, 50 runs or rounds per mode. Each campaign used its own
`TMPDIR` under `/tmp`, so no other agent's run could disturb it or be disturbed.

| Binary | Sequential, default threads | Sequential, `--test-threads=1` | 8 concurrent processes, whole binary | 8 concurrent processes, this test only |
|---|---|---|---|---|
| 1.85.0, before | 50/50 passed | 50/50 passed | 311 of 400 failed | 298 of 400 failed |
| 1.96.0, before | 50/50 passed | 50/50 passed | 307 of 400 failed | 297 of 400 failed |
| 1.85.0, after | 50/50 passed | 50/50 passed | 400/400 passed | 400/400 passed |
| 1.96.0, after | 50/50 passed | 50/50 passed | 400/400 passed | 400/400 passed |

All 1,213 failures before the fix carry the verifier's signature (`377:36`,
`NotFound`). The fix gives the test its own directory through the crate's port
of `File::TempDir`, `TempDir::new_in(std::env::temp_dir(), false)`, as
`tests/indexed_mzml.rs` already does. `+1.85.0 test --locked --all-features
--test ms_data_writing_consumer` passes 28/28.

### Integration gates

Run on kim (x86_64 Linux) at `ade680f`, the integrated tree without the
Boost.Regex facade, sequentially in one slot through
`~/.local/bin/openms-kim-gate.sh`:

| Gate | Result |
|---|---|
| `test --locked --all-features --all-targets --no-fail-fast` (1.96.0) | 327 binaries, 4675 passed, 0 failed, 1 ignored (the 256 MiB log-buffer case) |
| `+1.85.0 test --locked --all-features --all-targets --no-fail-fast` | 327 binaries, 4675 passed, 0 failed, 1 ignored |
| `+1.85.0 test --locked --no-default-features --all-targets --no-fail-fast` | 322 binaries, 3249 passed, 0 failed, 1 ignored |
| `test --locked --all-features --doc` | 65 passed, 0 failed |
| `clippy --locked --all-features --all-targets -- -D warnings` | exit 0, no warnings |
| `RUSTDOCFLAGS='-D warnings' doc --locked --all-features --no-deps` | exit 0, no warnings |

Local: `cargo fmt --all -- --check` and every `tools/*.py` checker in check
mode exit 0. The counts are summed over every `N passed` in the logs; cargo's
interleaved output hides one `test result:` line from a line-anchored count on
the stable log, which is why a naive count reads 4666 there. The per-test name
lists of the two toolchains agree.

## SQLite S1 resume checkpoint (2026-09-13, final snapshot)

GPT's uncommitted S1 resume, with the fixes from its Fable review and the
process-test harness fix, was validated on IBMI `kim` using node-local
`/scratch` sources, 32 compile jobs, rustc 1.96.0 and minimum
rustc 1.85.0. The run validated a snapshot hashed as
`7811316fb1153a3c3ebccd13f4bbd2493632495c2089b5af53a8fbc8448db744` (every file except `target/`); only this record and
the provenance manifests were edited afterwards.

| Check | Executed result |
|---|---|
| All features/all targets, current Rust | 4,437 passed |
| All features/all targets, Rust 1.85 | 4,437 passed |
| No default features, current Rust | 3,144 passed |
| All-feature doctests, each compiler | 57 and 57 passed |
| SQLite-only, both compilers, with and without `rusqlite/extra_check` | 20 SWATH tests in each of four selections |
| sqMass, both compilers, with and without `rusqlite/extra_check` | 33 handler and 11 activation tests in each of four selections |
| `tests/system_process.rs`, ten consecutive parallel runs per compiler | 10 and 10 of 10 runs passed |
| Formatting, full strict Clippy, rustdoc, focused SQLite/sqMass Clippy and rustdoc | passed |
| Provenance, coverage, documentation, module-cycle and feature-boundary gates | passed |

All 29 steps exited 0:

| Step | Exit | Seconds |
|---|---:|---:|
| `fmt` | 0 | 2.56 |
| `clippy` | 0 | 29.96 |
| `cargo_all_targets` | 0 | 109.25 |
| `portable` | 0 | 59.09 |
| `doctests` | 0 | 3.38 |
| `rustdoc` | 0 | 14.5 |
| `msrv_all_targets` | 0 | 130.25 |
| `msrv_doctests` | 0 | 8.04 |
| `sqlite_only` | 0 | 17.8 |
| `sqlite_only_extra_check` | 0 | 18.04 |
| `sqmass_only` | 0 | 21.95 |
| `sqmass_only_extra_check` | 0 | 21.88 |
| `msrv_sqlite_only` | 0 | 19.23 |
| `msrv_sqlite_only_extra_check` | 0 | 19.73 |
| `msrv_sqmass_only` | 0 | 23.96 |
| `msrv_sqmass_only_extra_check` | 0 | 23.88 |
| `swath_clippy` | 0 | 13.27 |
| `swath_rustdoc` | 0 | 3.48 |
| `sqmass_clippy` | 0 | 16.91 |
| `mzml_only_activation` | 0 | 20.96 |
| `msrv_system_process_x10` | 0 | 16.04 |
| `system_process_x10` | 0 | 15.99 |
| `check_core_sdk` | 0 | 0.11 |
| `test_core_sdk` | 0 | 0.05 |
| `core_sdk_coverage` | 0 | 0.2 |
| `test_core_sdk_coverage` | 0 | 0.22 |
| `check_doc_coverage` | 0 | 0.09 |
| `check_module_cycles` | 0 | 0.16 |
| `check_schema_feature_graph` | 0 | 47.71 |

The focused logs and `results.json` are retained in
[`tests/data/sqlite_s1_resume_validation/final`](../tests/data/sqlite_s1_resume_validation/final), each checked against the sha256 the run recorded;
the whole-tree logs remain at `/ceph/ibmi/abi/oliver/openms-rs/results/s1-resume-final-20260913-165118` and are identified by the sha256 values in
`results.json`. The review record and dispositions are in
[`tests/data/sqlite_s1_resume_validation`](../tests/data/sqlite_s1_resume_validation).
No macOS or Windows execution and no C++ differential is claimed.


## SQLite S0 checkpoint (2026-09-13)

The [public SQLite connector](SQLITE_CONNECTOR_SUPPORT.md) was validated on
IBMI `kim`, using node-local `/scratch` sources and build targets, 32 compile
jobs, Rust 1.96.0 and minimum Rust 1.85.0. SQLite remains optional.

| Check | Final executed result |
|---|---|
| All features, current Rust (`nextest` and all-target `cargo test`) | 3,854 passed |
| All features/all targets, Rust 1.85 | 3,854 passed |
| No default features, current Rust | 2,640 passed |
| SQLite-only, both compilers, with and without `rusqlite/extra_check` | 19 passed in each of four selections |
| Doctests, each compiler | 26 passed |
| Formatting, full strict Clippy and rustdoc | passed |
| SQLite-only strict Clippy with `extra_check` | passed in the focused worker run |
| Provenance, coverage and feature-boundary gates | passed |

All eight upstream connector class-test sections and 17 assertion macros are
mapped to named native tests. The 19 native tests also cover binary bindings,
literal identifiers, external non-UTF-8 schema names, lock errors and recovery,
transaction cleanup, late row errors, ignored binding tails and dependency
feature unification. The feature-boundary check on Rust 1.85 verifies that
default/no-default builds exclude SQLite and libxml, and SQLite-only excludes
libxml while selecting the bundled SQLite dependency.

The initial full run passed before review refinements. Claude Fable 5.1 then
reviewed the connector and re-reviewed the changes. Byte comparison removed
an unnecessary UTF-8 failure, a stronger destructor test proved release of the
write lock, and bound execution was made independent of the dependency's
optional precheck. An overstatement in the first review was corrected against
the pinned dependency source. Both reports and their dispositions are
[retained](../tests/data/sqlite_connector_validation/reviews.json).

A separate adapted C++ probe reproduced CPP-184, CPP-185 and the CPP-189
exception-documentation mismatch using the exact pinned connector implementation
and declarations, substitute support headers and host SQLite 3.45.1. Its
sources, binary and logs are retained outside this repository and hashed in
[the provenance manifest](../tests/data/sqlite_connector_provenance.json).
This is isolated adapted execution, not a full SDK build, original exception ABI
test or execution of the upstream class-test binary.

[The validation record](../tests/data/sqlite_connector_validation/summary.json)
contains final runtime hashes, commands, outcomes and retained remote log paths.
The final full run has 19 successful checks; the initial run has 17. Strong local
verification checked 2,007 distinct pinned source/registration/reference files.
No new macOS or Windows execution is claimed. The next storage stages remain
unported; their source findings and fixture inspection are recorded separately
in the [S1 plan](SQLITE_STORAGE_PLAN.md), without claiming native fixes.


## FORMAT integration checkpoint (2026-09-12)

The fifteen-module [FORMAT wave](FORMAT_WAVE_SUPPORT.md) was built and tested
on IBMI `kim`, with sources and target directories on node-local NVMe `/scratch`.
Rust 1.96.0 and minimum Rust 1.85.0 were used with 32 compile jobs. Independent
workers used separate scratch directories; final checks used frozen snapshots.

| Check | Executed result |
|---|---|
| All features, current Rust (`nextest` and `cargo test --all-targets`) | 3,835 passed |
| No default features, current Rust | 2,640 passed |
| All features/all targets, Rust 1.85 | 3,835 passed |
| Doctests, each compiler | 26 passed |
| Six minimum-feature FORMAT selections on Rust 1.85 | all passed |
| Full strict Clippy, rustdoc and formatting | passed after fixes |
| Six provenance/coverage/feature-graph checks | passed |
| Ten scientific-data/reference regeneration checks | passed |

The first combined run found a module-order formatting change, three test-style
Clippy diagnostics and a DTA documentation link that needed qualification.
Those corrections were followed by green checks, another full current-Rust test
run, the affected minimum-feature tests and an MSRV all-target build check.
Runtime source and Rust-test hashes match the final tested snapshot; later edits
record evidence and documentation. Strong local source verification checked
2,004 distinct pinned source/registration/reference files.

Five completed Claude Fable 5.1 review reports and their dispositions are
[retained with the validation artifacts](../tests/data/format_wave_validation/reviews.json).
Review found a Percolator line-limit regression and helped check order-sensitive
pepXML chemistry. Unsupported review claims were corrected against source and
regressions; a model's approval is not an SDK-wide correctness proof.

[The machine-readable record](../tests/data/format_wave_validation/summary.json)
contains commands, outcomes, snapshot hashes and remote log locations. Full logs
remain under `/ceph/ibmi/abi/oliver/openms-rs/results/format-final-20260912-214524`
and `format-final-20260912-215205`. Initial failed quality checks are retained,
not erased. This checkpoint adds no executed full-SDK C++ differential and makes
no new macOS/Windows validation claim. Public API gaps remain in the coverage
ledger and per-format support documents.

## Historical throughput measurement before FORMAT (2026-09-12)

The full sweep on the remote host was profiled after wave 1 rather than tuned by
assumption. Splitting `cargo test --all-features --all-targets` showed a rebuild
after touching `src/lib.rs` costs 10 s while running the already-built tests
costs 47 s: **82% of the wall clock was test execution**, not compiling.

`cargo test` runs each test binary in turn. With 225 integration binaries the
per-binary serialisation dominates, and no amount of build parallelism touches
it. `cargo-nextest` runs every test from every binary in one work-stealing pool:

| Runner | kim (384 c) | Mac (16 c) |
| --- | --- | --- |
| `cargo test --all-features --all-targets` | 48 s | 56 s |
| `cargo nextest run --all-features` | **9 s** | **18 s** |
| `cargo test --no-default-features` | 38 s | — |
| `cargo nextest run --no-default-features` | **8 s** | — |

At that earlier, smaller snapshot, the complete sweep — build, both test
selections, doctests, clippy, the MSRV 1.85 gate, rustdoc, fmt and the six
Python gates — took **42 s** on kim. This is not a timing claim for the later
FORMAT integration or for a cold build.

Two things did not help and are recorded so they are not retried:

- **More build jobs.** The default thread count already reaches 7.7 s; forcing
  384 gives 7.7 s. Compilation is 10 s incrementally, so raising
  `CARGO_BUILD_JOBS` past the current 96 changes nothing measurable.
- **A RAM-backed `target/`.** `/dev/shm` offers 1.2 TB, but a cold build there
  took 29 s against 22 s on node-local NVMe. The node's 2.2 TB of RAM already
  page-caches `/scratch`, so tmpfs only adds a copy.

This corrects an earlier judgement in this project's own plan, which stated that
nextest should not be installed because "full-suite time is compile/link of 225
test binaries, not test execution". That was asserted without measurement and is
wrong by a factor of four on both hosts.

Doctests are not a nextest feature and keep `cargo test --locked --all-features
--doc`. CI continues to use `cargo test`, so the runner change affects local and
remote development loops only; every test still runs in both. nextest runs each
test in its own process, which is stricter than `cargo test`'s shared-process
threads, and all 2,506 tests pass under it.

## imzML family: the five-header format, and what four audits found (2026-09-12)

The user withdrew an earlier deferral of `KERNEL/OnDiscImzMLExperiment.h` and then
set a standing rule that **all file formats stay in core**. The imzML family was
ported as a staged DAG on that basis: the reader core, then the writer and the
kernel-level facade in parallel, then `ImzMLFile` carrying the family's whole
class-test suite. It is not five headers but nine — the facade pulls in
`MSImagingGeometry`, `MSImagingRegion` and `IonImage`, and a tenth,
`IonImageExtraction.h`, is outside the registered union and recorded separately.

| Package | Headers | Class-test sections | Rust |
| --- | --- | --- | --- |
| Reader core | ImzMLHandlerHelper, ImzMLHandler | own tests; family suite is stage 3's | `src/format/imzml_handler.rs` |
| Writer | ImzMLWriter | round-trip against the reader | `src/format/imzml_writer.rs` |
| On-disc facade | OnDiscImzMLExperiment + 3 imaging headers | 17 ported, 22 mapped | `src/kernel/on_disc_imzml_experiment.rs` |
| File adapter | ImzMLFile | **50 ported, 0 mapped, 0 unaccounted** | `src/format/imzml_file.rs` |

All four audits returned zero blockers and `bounds_enforced: true`, which was the
property that mattered: an imzML dataset is two files, and every array offset and
length comes out of the XML and indexes into the companion `.ibd`.

**One finding should have been a blocker and two auditors found it independently.**
`infer_ibd_path` tested its suffix by byte-slicing `text[text.len() - 6..]`, which
panics whenever the sixth-from-last byte is a UTF-8 continuation byte. Every public
entry point in the family calls it first — `load`, `load_experiment`,
`load_into_consumer`, `load_spectra_index`, `store` and the facade — so
`dir/日本語.txt` aborted the process before a single validation ran. Reproduced,
then fixed with a char-boundary-safe `str::get`. The same rewrite closed a second
divergence: the suffix is now replaced by truncation, matching the source's
`p.substr(0, p.size() - 6) + ".ibd"`, where `PathBuf::set_extension` had treated a
name that is entirely `.imzML` as an extensionless hidden file and appended.

**Eight further majors were fixed rather than carried.** The reader returned `Err`
for a spectrum with exactly one external peak array, where ImzMLHandler.cpp:198-234
fills the non-external side from the inline peaks and succeeds; `extract_ion_image`
read arrays that `spectrum()` refused for the same pixel, because one path tested
the externality flags and the other did not; and the writer's float cvParam text
used Rust's `Display` where the source uses `std::to_chars` with `chars_format::fixed`
precision 15 inside [1e-2, 1e4) and shortest-round-trip scientific outside it
(`NumericFormatting.h:26-135`), its six vocabulary meta keys errored where
`DataValue::toString()` is deliberately lenient, and a misaligned auxiliary array was
skipped under default options but fatal under any sort or trimming filter.

**The test oracle was weaker than it claimed, and the re-run cleared the port.**
The ported `ClassTest::isRealSimilar` omitted the opposite-sign branch of
ClassTest.cpp:439-451, so `close(-1.0, 1.0)` returned true where C++ returns false,
and it backed 38 assertions. Ported faithfully, **all 38 still pass** — the weak
oracle was not masking a defect. The same package checked the three sibling suites
that carry a similar helper (`binned_spectrum`, `feature_handle`, `rich_peak2d`) and
established that their `|1.0 - ratio| <= 1e-5` form rejects negative ratios, so they
never had the defect; that was independently re-derived here before accepting it.

**One bound was closed by the integrator because no package owned the file.**
`ImzMLFile::preflight` sums `mz_length` from the index, but an inline array carries
its length in the XML, so inline peaks reached the caller uncounted once the reader
began decoding them — bounded only by the 512 MiB XML ceiling, roughly 96 million
`f32` peaks. A `PeakBudget` now charges every decoded spectrum against
`max_loaded_peaks` in both load loops, failing as soon as the ceiling is crossed.

Nineteen further C++ defects are recorded as CPP-124 to CPP-142.

Gates on the Linux node at the integration commit: nextest all-features **3,193
passed**, no-default-features 2,376 passed, doctests 16 passed, clippy `-D warnings`
clean, `cargo +1.85.0 check --all-features --all-targets` clean, rustdoc
`-D warnings` clean, `cargo fmt --check` clean, all six Python gates green.

A correction to the preceding checkpoint: the sweep it cites for kernel wave 3 was
killed by an ssh disconnect after the all-features run, so its no-default-features,
clippy, MSRV, rustdoc and Python gate results were never obtained. Those gates are
verified here, on a tree that contains that work.

## Kernel wave 3 completion and its audit fixes (2026-09-12)

The three packages a session limit had killed were re-run from the integrated
base and merged: MSExperiment/AreaIterator residuals, the MSSpectrum ion-mobility
quartet with ConsensusFeature's last two members, and the OnDiscMSExperiment
facade. Kernel headers closed or native-equivalent rise from 27 of 34 to 30 of 34,
and `IMTypes.h` closed with a first entry for `SpectrumSettings.h` — the ion-mobility
quartet turned out to be declared in METADATA, not KERNEL.

Three audits returned, none with a blocker, and one reported
`section_audit_honest: false`. Every finding was fixed in a follow-up wave rather
than merged as-is. Three are worth recording because each was a claim the code did
not support.

**The test-mapping rule was circumvented.** The MSExperiment package reported "63
sections mapped with a cited asserted value". Its support doc accounted for all 54
mapped `MSExperiment_test.cpp` sections with a bare list of nine test *files* — no
function, no value. The auditor counted assertion macros per section and found 22
above the five-macro threshold that mandates porting, one of them with 64 macros;
sections 5 and 6 (copy and move assignment, asserting `getMinMZ 5.0`, `getMaxMZ
10.0` and a moved-from size of 0) had no Rust evidence anywhere in the repo. All 22
are now ported. The honest accounting is 56 ported, 21 mapped with a named function
and value, and **4 unaccounted** — `set2DData<add_mass_traces=true>`, both
`getFirstProductSpectrum` overloads and `operator<<`, whose members are unported.
`MSExperiment.h` stays `partial` for exactly those four.

**A rustdoc claim misdescribed the C++ it cited.** `area_iteration.rs` said a
reversed low/high ion-mobility pair "silently selects nothing" upstream.
`AreaIterator.h:277` builds `RangeMobility{low_im_, high_im_}`, and
`RangeBase(min,max)` (RangeManager.h:48-52) *throws* `InvalidRange` when `min > max`.
The port's `Err` agreed with the source by accident, not by the stated reasoning.
Several line anchors had drifted and an OpenMP note claimed the serial rasterizer
"computes the same image" — true only for `Max`, since the parallel branch merges
per-thread f32 buffers and so differs for `Sum`.

**An audit found a defect in this port's own mzML reader.** A fixture substitution
in the OnDisc package was documented as forced by an unavailable upstream file. The
file is committed in-tree, byte-identical to the pinned copy; what blocked it was
this reader rejecting `binaryDataArrayList count="2"` with four children. Checking
every list handler in `MzMLHandler.cpp` showed upstream **never** compares a declared
count against the actual number of children: `binaryDataArrayList` feeds only
`bin_data_.reserve` (:1015), `selectedIonList`'s count only warns when above one
(:1371), and `precursorList`, `productList`, `scanWindowList` and
`referenceableParamGroupList` have no list handler at all. This port hard-errored at
four sites. All are advisory on reading now, in both the reader and the `loadSize`
counting path, while every declared count remains a **resource ceiling** enforced
before allocation and writing still emits the true count. The earlier header-list fix
(`src/format/mzml_header.rs:107-116`) had addressed only one instance of this defect
class; this is the general case, and it had blocked real data twice.

The last item cost a deliberate reversal. The count fix initially kept
`referenceableParamGroupList` strict so `read` and `read_size` would agree, and pinned
that with a test. Relaxing only the reader would have left `read` accepting a document
`read_size` rejects, so both paths were relaxed together and the test rewritten to
assert that the two readers still agree — the property the strict check had existed to
protect.

Gates at the integration commit, on the 384-core node: build 21 s, nextest
all-features 2,917 passed, no-default-features 2,362 passed, doctests green, clippy
`-D warnings` clean, `cargo +1.85.0 check --all-features --all-targets` clean, rustdoc
`-D warnings` clean, `cargo fmt --check` clean, all six Python gates green. Four further
C++ defects recorded as CPP-120 to CPP-123.

## Kernel wave 3, partial: the map containers (2026-09-12)

Wave 3 launched four work packages. **Three were killed mid-run by an account
session limit** (WP7 MSExperiment/AreaIterator residuals, WP11 residual closure,
WP12b the OnDisc facade); they produced no commits and are queued for a clean
re-run. WP10 completed and is integrated here. A concurrently launched imzML wave
died the same way before its first stage committed.

| Work package | Headers | Class-test sections | Rust |
| --- | --- | --- | --- |
| Map containers | FeatureMap, ConsensusMap, ConversionHelper | 74 ported | `src/kernel/map_operations.rs`, `src/kernel/conversion_helper.rs` |

FeatureMap.h and ConsensusMap.h had never been reviewed member by member — both
sat at `evidence_requires_review` with no entry in the reviewed-API ledger at all.
All 74 upstream sections (32 + 39 + 3) are ported with transcribed literals, tier 3
evidence under [the differential validation policy](DIFFERENTIAL_VALIDATION.md);
none is merely mapped. Nine further C++ defects are recorded as CPP-111 to CPP-119.

**The audit found three documentation defects, all fixed before merge.** Each was
a claimed equivalence the code does not have: `isMapConsistent` was documented as
a plain mapping to `validate_consistency()`, which is in fact strictly stricter —
the source checks only duplicate column descriptions and unregistered handle map
indices, while the Rust additionally rejects a bad `experiment_type`, non-finite
coordinates, duplicate unique IDs and invalid attached records, so a map the
source calls consistent can return `Err`. Both `updateRanges()` rows omitted that
`ranges()` is fallible where the source cannot fail, because it opens with
`self.validate()?`. And the five `FeatureMap` sorts are stable here but use
`std::sort` upstream, which is not — a divergence the document was meticulous
about elsewhere and silent about here. The range *content* was verified faithful
in both maps.

**One finding was a defect in the integrator's own instructions.** The package
added four lines to `src/kernel.rs` where the rule allowed two. It registered two
modules, and the rule had assumed one module per package; the wording is corrected
for the re-run rather than charged against the package.

Gates, independently re-run by the auditor at the package commit: MSRV
`cargo +1.85.0 check --all-features --all-targets` clean, nextest all-features
2,825 passed, no-default-features 2,309 passed, clippy `-D warnings` clean,
rustdoc `-D warnings` clean, and 100% rustdoc on both new modules.

## Kernel wave 2: ion mobility, chromatogram merging, feature identification, MRM and indexed mzML (2026-09-12)

Five work packages ported in parallel git worktrees, each followed by an
**independent adversarial audit in its own worktree** that re-ran every gate
itself rather than trusting the port's report. Kernel headers closed or
native-equivalent rise from 18 of 34 to 24 of 34.

| Work package | Headers | Class-test sections | Rust |
| --- | --- | --- | --- |
| Ion mobility | MSSpectrum (mobility surface) | 68 ported, 3 mapped | `src/kernel/spectrum_mobility.rs` |
| Chromatogram merging | MSChromatogram, Mobilogram | 91 ported | `src/kernel/chromatogram_merge.rs` |
| Feature identification | BaseFeature, Feature, ConsensusFeature | 92 ported | `src/kernel/feature_identification.rs` |
| MRM | MRMFeature, MRMTransitionGroup | 44 ported | `src/kernel/mrm.rs` |
| Indexed mzML | FORMAT/HANDLERS/IndexedMzMLHandler | 16 ported | `src/format/indexed_mzml_handler.rs` |

311 upstream class-test sections were ported and 3 mapped with cited evidence,
0 unaccounted. All transcribed literals, which is tier 3 evidence under
[the differential validation policy](DIFFERENTIAL_VALIDATION.md). No C++ was
executed. Thirty-five further C++ defects were recorded as CPP-076 to CPP-110,
including an `MRMFeature` lookup by unknown key that default-inserts into the map
and returns the first feature, a `BaseFeature::sortPeptideIdentifications`
comparator that is not a strict weak ordering and mutates its own arguments, and
an `IndexedMzMLHandler::openFile` that accumulates index state instead of
replacing it.

**No audit returned a blocker.** Every auditor independently re-ran
`cargo +1.85.0 check --locked --all-features --all-targets`, the full nextest
sweep, clippy, rustdoc and the doc-coverage report, enumerated each owned
header's public members by hand against the support doc's API table, and counted
`START_SECTION`s itself. Four verdicts were `accept_with_fixes` and one `accept`.

**Both major findings were claim defects, not code defects, and both were
fixed before merge.** The ion-mobility package proposed a ledger entry for
`MSSpectrum.h` that would have *replaced* the accumulated record of three earlier
waves — `documentation` and `scope` pointed only at the new support document.
Applied verbatim it would have erased the InstrumentSettings, AcquisitionInfo,
`getType(bool)` precedence and record-metadata review history. The integrator
merged instead: `rust` and `tests` unioned, the prior scope preserved and the new
wave's sentence appended naming its document. The indexed-mzML package documented
`setSkipXMLChecks` as "not ported" in three table rows and its provenance
manifest; the auditor traced `options_mut().skip_xml_checks` through
`mzml::read_with_load_options` to the Base64 whitespace strip at
`src/format/mzml.rs:2386`, which is the entire effect of the flag in the source
as well. All four claims were corrected to record the member as ported.

Three minor findings against the ion-mobility package are recorded rather than
fixed: an unreachable `DriftTimeUnit::None` branch documented as reachable, a
non-finite rejection whose stated rationale covers NaN but not infinities, and
the bottom-up chunk merge whose stability argument is untested because both
presorted fixtures use all-distinct m/z. The last is a real coverage gap — no
case has equal m/z spanning a chunk boundary.

Full sweep on the 384-core node at the merge commit: build 32 s, nextest
all-features 2,763 passed, doctests 8 passed, nextest no-default-features 2,247
passed, clippy `-D warnings` clean, `cargo +1.85.0 check --all-features
--all-targets` clean, rustdoc `-D warnings` clean, `cargo fmt --check` clean, and
all six Python gates green after the integrator regenerated the three generated
files.

## Kernel wave 1: ranges, predicates, helpers, gap-0 review and geometry (2026-09-12)

Five work packages ported in parallel git worktrees, then rebased onto wave 0,
merged and verified together. Kernel headers closed or native-equivalent rise
from 8 of 34 to 18 of 34.

| Work package | Headers | Class-test sections | Rust |
| --- | --- | --- | --- |
| Ranges | RangeManager, SpectrumRangeManager, ChromatogramRangeManager | 48 ported | `src/kernel/ranges.rs` |
| Predicates | RangeUtils | 33 ported | `src/kernel/range_utils.rs` |
| Helpers | SpectrumHelper | 9 ported | `src/kernel/spectrum_helper.rs` |
| Gap-0 review | DPeak, StandardTypes, RichPeak2D, FeatureHandle, BinnedSpectrum | 51 ported | `src/kernel/gap_closures.rs` and existing modules |
| Geometry | DPosition, DIntervalBase, DRange | 104 ported | `src/data_structures/{dposition,dinterval,drange}.rs` |

245 upstream class-test sections were ported with transcribed literals, which is
tier 3 evidence under [the differential validation policy](DIFFERENTIAL_VALIDATION.md).
No C++ was executed. Fifteen further C++ defects were recorded as CPP-061 to
CPP-075, including a `makePeakPositionUnique` swap that discards the whole
spectrum record while warning only about data arrays, and a `DRange::united` of
two empty ranges that returns the universal range.

**Ranges are computed on demand, deliberately.** The source caches ranges in a
mutable member refreshed by `updateRanges()`. Peak vectors are public here, so a
cache cannot be invalidated soundly; the algebra is ported as pure values and the
inherited container surface becomes `range_manager()` accessors. All 21 TOPP
`updateRanges` call sites were audited before choosing this: every one is a plain
update-then-read, and the mutable `getRange()` has no mutating caller in core, so
no ported tool's numbers change. The combined experiment role now folds in
chromatogram retention time, intensity and **product m/z**, which the previous
`MSExperiment::ranges()` omitted entirely; a ported `FileInfo` would have printed
a narrower m/z range.

**Two defects came from stale worktree bases.** Three of the five agents branched
from `a463e3e`, 26 commits behind, so their own green gates were green against a
tree without the typed record metadata, the drift-time fields or the 1.85 fixes.
Rebasing exposed both: two SpectrumHelper test assertions used the pre-migration
metadata API, and `copySpectrumMeta` documented a drift-time deferral that wave 0
had already made obsolete. The implementation was correct by construction; the
rustdoc, support document and API table were not, and nothing tested it. Each
branch was rebased onto wave 0 and re-verified before merging.

**The no-default-features CI line was broken by wave 0 and is now fixed.** Gating
`pub mod cli` on `paramxml` left the five `tests/topp_*.rs` files using
`openms::cli` without a gate of their own, so `cargo test --locked
--no-default-features` (rust.yml line 22) failed to compile. Verified after wave 0
were `--all-features` and the Python gates, not that line. Each file now carries
`#![cfg(all(feature = "mzml", feature = "paramxml"))]`.

Recorded checks on the integrated tree (16-core Apple Silicon, cargo 1.96):

| Check | Result |
| --- | --- |
| `cargo build --locked --all-features --all-targets` | clean |
| `cargo test --locked --no-default-features` | clean, after the gate fix |
| `cargo fmt --all -- --check` | clean |
| `python3 tools/check_core_sdk.py` | 1,060 added source references agree |
| `python3 tools/core_sdk_coverage.py --write` | complete 16, native-equivalent 70, partial 20, evidence-requires-review 159, unmapped 521 |
| `python3 tools/check_doc_coverage.py --write` | 1,764 of 3,401 public items = 51.9%, all new modules at 100% |

**A remote build host caught a defect the local integration missed.** The full
sweep was also run on an IBMI HPC node (`kim`, 384 cores). It reported two
unresolved intra-doc links in `spectrum_helper` that the local run had not been
repeated after merging. The cause was the wave-1 merge resolution itself: adding
an outer `///` doc comment on `pub mod spectrum_helper` in `src/kernel.rs` makes
rustdoc resolve that module's inner `//!` links in the parent `kernel` module,
so `[`PeakContainer`]` reported "no item named `PeakContainer` in module
`kernel`" although the trait exists in the module. One work package had warned
of this mechanism in its own notes. Both links now use full crate paths, and the
module records why.

### Remote build host

`kim` was provisioned as a second verification host: both toolchains in shared
CephFS home (install once, visible on every node), source tree and `target/` on
node-local NVMe `/scratch`, and libxml2 2.15.4 with pkg-config in a shared
micromamba environment under `/ceph/ibmi/abi/oliver/envs/rustbuild`. Three
environment facts had to be discovered and are recorded in
`/scratch/kohlbach/openms-rs-env.sh`: `pkg-config` must be on `PATH` for the
`libxml` build script, the node ships `libclang.so.1` without its resource
headers so bindgen needs GCC 13's include directory, and the default
`ulimit -n` of 1024 starves parallel `rustc`.

| Check | kim (384 c, 96 jobs) | Mac (16 c) |
| --- | --- | --- |
| `build --all-features --all-targets` | 22 s | 16 s |
| `test --all-features --all-targets` | **48 s**, 2,506 passed | **176 s**, 2,235 passed |
| `test --no-default-features` | 31 s, 2,028 passed | — |
| `clippy --all-features --all-targets` | 20 s | — |
| `cargo +1.85.0 check --all-features --all-targets` | 16 s | clean |
| six Python gates | all pass | all pass |

The test suite runs 3.7 times faster; the cold build does not, because the
toolchain is read over CephFS. `target/` reaches 23 GB there against 4 GB
locally, which node-local `/scratch` absorbs.


## Kernel wave 0: MSRV, build baseline and scaffold (2026-09-12)

Preparation for the parallel kernel port. Three findings are recorded because
each corrects a claim made earlier in this project.

**The minimum Rust version was broken.** `cargo +1.85.0 check --locked
--all-features --all-targets` failed on release `5688775` with `E0658` at five
let-chain sites: `src/cli.rs:128,345`, `src/metadata/experimental_design.rs:163,721`
and `src/format/experimental_design_file.rs:178`. Let-chains stabilised in Rust
1.88; edition 2024 accepts the syntax, so rustc 1.96 never reported it and the
crate's own `rust-version = "1.85"` was not enforced locally. CI `minimum-rust`
was red. All five are rewritten as nested `if`s; the 35 tests covering those
sites pass unchanged. An adversarial review found three of the five; the 1.85
compiler found the other two, and is now the gate — not a grep.

**Build baseline, the first recorded.** Before any change `target/` held 27 GB
(53.7 GiB of files) for 225 test binaries with full DWARF. With
`[profile.dev] debug = "line-tables-only"` (test semantics unchanged):

| Check | Result |
| --- | --- |
| `cargo test --locked --all-features --all-targets --no-run --timings` | 16 s, clean tree |
| `cargo test --locked --all-features --all-targets` | 2,231 passed, 176 s |
| `cargo test --locked --all-features --doc` | 4 passed (excluded by `--all-targets`; 2,235 total) |
| `target/` after the full build | 4.0 GB |
| `cargo +1.85.0 check --locked --all-features --all-targets` | clean |
| `cargo +1.85.0 check --locked --no-default-features` | clean |
| `cargo fmt --all -- --check` | clean |
| `check_core_sdk`, `core_sdk_coverage`, `test_core_sdk_coverage`, `check_schema_feature_graph`, `check_doc_coverage` | pass |

Machine: 16-core Apple Silicon, cargo 1.96. `cargo-timing.html` is kept outside
the repository.

**The struct-literal risk was phantom.** The plan feared ~400 struct literals
would break when fields were added. A regex over `MSSpectrum {` counted return
types, `impl` blocks and closure bodies; the compiler reports **zero** `E0063`
missing-field errors after adding `MSSpectrum::{drift_time, drift_time_unit}`,
`MSExperiment::sql_run_id`, `BaseFeature::{primary_id, id_matches}` and
`ConsensusFeature::ratios`. Every real literal already used
`..Default::default()`, the crate's existing convention. An automated fixer
built on the same regex was reverted in full rather than patched.

Scaffold additions: the fields above, `Ratio` (ports `ConsensusFeature::Ratio`:
`ratio_value`, `denominator_ref`, `numerator_ref`), and `Error::{InvalidRange,
MissingInformation}` mapped to `ILLEGAL_PARAMETERS` and `MISSING_PARAMETERS`.
`src/kernel.rs`, `src/kernel/features.rs` and `src/error.rs` remain at 100%
rustdoc coverage. Identification data is attached to maps **by reference**, a
deliberate divergence from `FeatureMap.h:294`, because the graph is not `Clone`
and embedding it would strip `Clone`/`PartialEq` from both map types.

## Indexed mzML writing and binary normalization (2026-09-11)

[Recorded checks](mzml-output-validation.json) cover complete represented
[source writer options](MZML_WRITE_OPTIONS_SUPPORT.md) and the
[binary whitespace option](MZML_NORMALIZATION_SUPPORT.md). The group adds 24 tests;
MzMLFile remains partial while typed/noise/detector transport and XSD validation
are implemented separately. Spectrum mobility/IMPeakType and additional precursor
activation metadata/unit routes also remain.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,139 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/minimal mzML-validation selection | 373 passed |
| Strict Clippy on both compilers; release library, Rustdoc, Rustfmt | Passed |
| Source audit, extraction integrity and two independent fixture regenerations | Passed |
| Completion checks and optional checksum dependency boundary | Passed |

All 15 integrated checks passed on their first run. Selected counts overlap.
Twenty-six final extraction files remain exact; Cargo.toml retains only the
previous semantic-feature registration in addition to the frozen writer change.
The final reader's normalization patch starts from the exact reviewed writer
and earlier centroid/semantic registrations. All 28 source hashes and nine
fixture/generator/reused-resource hashes were independently verified.

The writer reuses one prepared header/binary payload across two precharged markup
passes. Independent Python checks verify the actual SHA-1 prefix and every
record/index ID and byte offset, including partial writes and UTF-8. The unchanged
indexed XSD passes actual validation; this is separate from index integrity.
Original Numpress literals and independent precision/whitespace input fixtures
are retained. CPP-049/050 are corrected, and both normalization settings retain
checked malformed-Base64 rejection for CPP-055. No C++ method execution or full
C++ SDK build is claimed.

The [published ProForma checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34588912206)
passed CI. The ledger remains at 72 complete/native-equivalent headers and 714 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,945 distinct current files, 220 historical references, 21 graph references
and 952 added references. The source issue log contains 56 entries, including the newly source-reviewed
CPP-056 synthetic-scan mobility omission. Full SDK
completion remains outstanding.

## mzML inspection, isolation loading and semantic validation (2026-09-11)

[Recorded checks](mzml-operations-validation.json) cover
[spectrum type/centroid inspection](MZML_CENTROID_SUPPORT.md),
[isolation-target loading](MZML_ISOLATION_SUPPORT.md) and the complete native
[MzMLValidator specialization](MZML_VALIDATOR_SUPPORT.md). This checkpoint adds
47 tests and closes one public header; MzMLFile itself remains partial.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,115 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/minimal mzML-validation selection | 349 passed |
| Spectrum type queries without format features, Rust 1.85 | Seven passed |
| Strict Clippy on both compilers; release library, Rustdoc, Rustfmt | Passed |
| Source audit, source projection, extraction integrity and completion checks | Passed |

All 15 integrated checks passed on their first run. Selected test counts overlap.
Of 34 integrated extraction files, 33 remain exact; the shared reader merges
only the frozen isolation implementation and centroid/semantic registrations.
All 47 source hashes and 16 fixture/resource/tool/projection hashes were checked,
including the unchanged DTA files, original validator inputs, mapping and reused
ontologies. Native tests verify cumulative limits, source event order, caller
option preservation and fresh validator state across documents. No C++ method
execution or full SDK build is claimed.

The [published XLMS checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34586771858)
passed CI. The issue log now contains 55 entries. New findings concern mixed
chromatogram primary roles, ineffective indexed-schema ID references and unchecked
source Base64 alphabet bytes. The [index-schema probe](mzml-index-schema-probe.json)
records actual acceptance of a dangling reference by the unchanged XSD; it does
not claim byte-offset or checksum validity.

The ledger records 72 complete/native-equivalent headers and 714 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,942 distinct current files, 220 historical references, 21 graph references
and 924 added references. Indexed writing, complete writer options, typed/noise
transport, runtime XSD validation, whitespace-normalization options and broader
SDK work remain outstanding.

## ProForma spectrum generation (2026-09-11)

[Recorded checks](proforma-spectra-validation.json) cover all six
[ProForma spectrum operations](PROFORMA_SPECTRA_SUPPORT.md), completing the pinned
ProForma public header's native operation groups. Source compatibility and
checked boundaries remain explicit; this is not full notation-standard certification.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,068 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/no-default library and adjacent chemistry suites | 274 passed |
| Strict Clippy on both compilers; release library, Rustdoc, Rustfmt | Passed |
| Source audit, ProForma fixture regeneration and completion checks | Passed |

All 12 integrated checks passed on their first run. Selected counts overlap;
this group adds 21 tests. Twelve frozen extraction files remain exact; the
thirteenth has only a root module-documentation correction. All 30 source hashes
and both fixture/generator hashes were independently verified. The fixture
retains ten full source sections and 18 literal assertions, including six
generation calls; it does not invent numerical spectrum oracles. Separate tests
exercise real-backend composition, repeated resolution/warning order, shared
resource limits and atomic registry publication. CPP-038/047 finite source
behavior is retained beside independent chemical/position expectations.
No C++ spectrum execution or full SDK build is claimed.

The [published general-validator checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34585986090)
passed all CI jobs. The source issue log contains 52 entries, including a
cross-document mzML parameter-group state leak and four-line XSD schema-selection
defect. The latter has independent pinned-schema validation evidence in
[mzml-schema-selection-probe.json](mzml-schema-selection-probe.json). The ledger now
records 71 complete/native-equivalent headers and 715 requiring implementation
or review, with zero certified TOPP workflows. Source verification covers 1,933
distinct current files, 220 historical references, 21 graph references and 877
added references. Remaining SDK work is tracked in the completion ledger.

## Streaming mzML consumers (2026-09-11)

[Recorded checks](consumer-validation.json) cover the complete unconditional
MSDataConsumer interface, both source [transform operation groups](MZML_CONSUMER_SUPPORT.md)
and disabled scientific data population with `fill_data=false`.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,047 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/mzML-only library and adjacent selection | 236 passed |
| Rust 1.85/unconditional interface without format features | One passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, consumer fixture regeneration and completion checks | Passed |

All 14 integrated checks passed on their first run. Selected counts overlap;
the full suite adds 24 tests. All 15 frozen extraction files remain exact.
Independent regeneration reproduces both fixture files after five documented
source child-count repairs, retaining original scientific data. Both transform
modes reproduce the source test's four spectra, 40 peaks and TIC 350. Source
control-flow review and independent tests cover setup order, separate pools,
mutation/stop/error boundaries, descriptor validation and atomic destination
publication. No C++ consumer execution or full SDK build is claimed.

The [published header checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34584903921)
passed cross-platform, minimum-Rust and quality CI. The C++ issue log contains 50
entries, including the indexed writer's placeholder checksum and empty dummy index. The ledger records
70 complete/native-equivalent headers, 716 requiring implementation or review,
and zero certified TOPP workflows. Source verification covers 1,933 distinct
current files, 220 historical references, 21 graph references and 847 added
references. Remaining mzML source options/validation, centroid inspection,
indexed output, ProForma wrappers and broader SDK work are tracked separately.

## Crosslink spectrum generation (2026-09-11)

[Recorded checks](xlms-validation.json) cover the complete class-specific
[XLMS generator](THEORETICAL_XLMS_SUPPORT.md), all 25 source options, and the
ProteinProteinCrossLink record/reaction enum. Other OPXL records remain separate.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 2,023 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/no-default library and adjacent selection | 176 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, XLMS fixture extraction and completion checks | Passed |

All 12 integrated checks passed on their first run. Selected counts overlap;
the full suite adds 25 tests. All ten frozen extraction files remain exact.
The independent extractor reproduces 52 literal source masses and 113 allowed
annotation strings. Native tests additionally cover source branches, cumulative
limits, sequence identity and atomic aligned appending. Review corrected linked
loss endpoint guards, empty-alpha pair behavior and portable charge-span handling
before integration. No C++ spectrum execution or full SDK build is claimed.

The generator explicitly retains the finite source suffix-loss and precursor
isotope defects (CPP-042/043), with separately calculated chemical expectations.
The C++ issue log contains 48 entries. The ledger records 69 complete or
native-equivalent headers, 717 requiring implementation or review and zero
certified TOPP workflows. Source verification covers 1,932 distinct current
files, 220 historical references, 21 graph references and 835 added references.
ProForma wrappers, mzML consumers and broader XLMS analysis remain ongoing work.

## General semantic validation (2026-09-11)

[Recorded checks](semantic-validation.json) cover the complete class-specific
[SemanticValidator group](SEMANTIC_VALIDATOR_SUPPORT.md) and the shared XML reader
used by CV mappings. Enable the optional `semantic-validation` feature.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,998 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/semantic feature and adjacent selection | 163 passed |
| Rust 1.85/CV mapping feature and library selection | 127 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, fixture projections and completion checks | Passed |

All 16 integrated checks passed on their first run. Selected counts overlap;
the full suite adds 19 tests. All 14 reviewed extraction files remain
byte-identical after integration. The historical valid fixture produces no
diagnostics; the corrupt fixture reproduces all five error and four warning
messages in source order using its original 738-term vocabulary. Four raw
fixtures retain their source bytes. No C++ execution or full SDK build is claimed.

The native validator corrects descendant-unit lookup, failed-parse state leakage
and history-dependent missing-path lookup (CPP-039/040/044). It retains the source
date conversion behavior (CPP-046); complete XSD conformance and derived format
validators are separate work. The C++ issue log contains 47 entries, including
the newly recorded ProForma flattened-range position defect.

The ledger records 68 complete/native-equivalent headers and 718 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,930 distinct current files, 220 historical references, 21 graph
references and 824 added references. Streaming consumers and the separately
reviewed XLMS/ProForma spectrum work are outside this checkpoint.

## mzML headers and peptide-evidence keys (2026-09-11)

[Recorded checks](header-evidence-validation.json) cover the complete
[source-supported header/reference group](MZML_HEADER_SUPPORT.md), metadata-only
reading and [PeptideEvidence value/key operations](PEPTIDE_EVIDENCE_SUPPORT.md).
DataProcessing timestamps now use DateTime across mzML, FeatureXML and ConsensusXML.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,979 tests passed |
| Doctests, Rust 1.98 | Four passed |
| Rust 1.85/header and adjacent selection | 321 passed |
| Rust 1.85/native selection | 171 passed |
| Final Rust 1.85/path and precursor workflows | 14 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, reference regeneration and completion checks | Passed |

Selected totals overlap. The full suite adds 26 tests. Two old rejection tests
needed capability updates: metadata-only paths and scalar precursor metadata
now succeed. Their replacement assertions verify the new results and retain
unsupported-value/error checks. An initial positive path test used a historical
projection lacking a required processing reference; the final test uses a complete
writer-produced document. All failures and successful reruns are retained in the
record, without claiming first-run success.

All 35 reviewed header files remain byte-identical after integration. Six
fixtures and the formatted mapping table regenerate from the pinned C++ source;
116 reader/writer instrument pairs agree. Both rich-header writer outputs pass
actual independent XSD validation. No new C++ execution is claimed.

The [published parent checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34582838457)
passed all five cross-platform/minimum/quality jobs. The C++ issue log now has
46 entries. The completion ledger records 67 complete/native-equivalent headers,
719 requiring implementation or review, and zero certified TOPP workflows.
Source verification covers 1,926 distinct current files, 220 historical references,
21 graph references and 799 added references. Consumers, ProForma spectra and
separately staged semantic validation remain outside this checkpoint.

## ProForma conversion and CV mappings (2026-09-11)

[Recorded checks](conversion-mapping-validation.json) cover complete ProForma
AASequence conversion in both directions and the five-class CV mapping group,
adding 30 tests against SDK `82ce5b3`.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,953 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Rust 1.85/CV mapping and library selection | 126 passed |
| Rust 1.85/native and adjacent selection | 260 passed |
| Rust 1.85/JSON and mass selection | 142 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, independent projections and completion ledger | Passed |

The [previous published checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34580543352)
passed all Linux, macOS, Windows, minimum-Rust and quality jobs. This is distinct
from CI for the present increment.

All 21 integrated checks passed on their first run. Selected totals overlap.
Every one of 1,078 executed C++ formatter cases is exercised through native
reverse conversion, writing and parsing. Only the exact source formatter helper
was compiled; this is not an executed C++ conversion or spectrum oracle. Mapping
tests compare all 683 full records projected independently from six unchanged
source XML files. Complete operation and evidence review is recorded in the
[conversion](PROFORMA_CONVERSION_SUPPORT.md) and [mapping](CV_MAPPING_SUPPORT.md)
support documents.

The [C++ issue log](../OpenMS_CPP_ISSUES.md) contains 38 entries. New source-reviewed
findings document charge loss during formula combination and duplicate linker
mass in the still-unported ProForma XLMS spectrum wrapper. Proposed C++ fixes,
executed evidence and native compatibility behavior remain separate.

The ledger records 66 complete or native-equivalent headers and 720 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,925 distinct current files, 220 historical references, 21 graph references
and 776 added references. Full mzML headers/consumers, ProForma spectra and general
semantic validation remain ongoing groups; staged work is not certified here.

## ProForma mass and controlled vocabularies (2026-09-11)

[Recorded checks](vocabulary-mass-validation.json) cover the complete ProForma
mass operation group and ControlledVocabulary at SDK `82ce5b3`, adding 37 tests.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,923 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Rust 1.85/native unit and adjacent selection | 240 passed |
| Rust 1.85/ProForma JSON and mass selection | 138 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, independent projections and completion ledger | Passed |

All 17 integrated checks passed. Selected totals overlap. Root review covered
the complete new operation groups, including shared registry transactions,
source mass pass order, iterative graph behavior and cumulative allocation
limits. The vocabulary tests compare every field of all 9,254 final terms and
all 16,852 name aliases with an independent source-loop projection. Five raw
providers retain their original bytes and separate data-license notices.
No executed C++ mass or vocabulary differential result is claimed.

The [previous checkpoint's CI](https://github.com/okohlbacher/openms-rs/actions/runs/34578513534)
passed every Linux, macOS, Windows, minimum-Rust and quality job, confirming the
Windows resource checkout correction. That result is distinct from CI for
this later change.

The [C++ issue log](../OpenMS_CPP_ISSUES.md) now contains 36 stable entries,
including ten newly source-reviewed defects found during subsequent header,
conversion and mapping work. Executed evidence, source deductions, proposed
upstream fixes and native handling remain distinguished.

The completion ledger records 61 complete or native-equivalent headers and
725 still requiring implementation or review, with zero certified TOPP workflows.
Source verification covers 1,915 distinct current files, 220 historical references,
21 graph references and 734 added references. Full mzML headers/consumers,
ProForma sequence conversion/spectra and CV mapping/semantic validation remain
separate ongoing groups.

## Experiment settings, DateTime and ProForma resolution (2026-09-11)

[Recorded checks](settings-resolution-validation.json) cover the integrated
settings ownership migration, DateTime operations and ProForma resolver at
SDK `82ce5b3`, with 42 additional tests.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,886 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Rust 1.85/native unit and adjacent selection | 400 passed |
| Rust 1.85/JSON and RNA selection | 155 passed |
| Rust 1.85/mzML and settings selection | 103 passed |
| Final private-test repeats, current/minimum | 110/102 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed |
| Source audit, independent projections and completion ledger | Passed |

Selected totals overlap the full suite. The record retains an initial minimum
lint failure and its successful retry after correcting a test initializer.
The original-file metabolite example still produces 81 features; its output
passes the original schema and 1,458 scalar/typed-metadata comparisons with the
source expected file. Default hull omission remains covered by direct tests.

DateTime retains 301 executed unmodified-source probe rows: 266 native matches
and 35 deliberately corrected early-year calendar results. Independent Python
month stepping regenerates all 35 corrections exactly. Three separate UBSan
executions expose signed fractional-second overflow. The
[C++ issue log](../OpenMS_CPP_ISSUES.md) records these and source-reviewed defects,
with an independent sorted-weight reproduction of the IMS nonprogressing witness.
Proposed upstream fixes are distinct from the native compatibility policy.

The [previous checkpoint's CI](https://github.com/okohlbacher/openms-rs/actions/runs/34576479236)
failed on Windows because Git converted the original SVM resources to CRLF.
Explicit byte-preservation rules fix that checkout issue; all four resources
now remain identical under Windows-style conversion. Cross-platform CI for
this checkpoint is a separate check, not inferred from local success.

The ledger records 60 complete or native-equivalent headers and 726 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,910 distinct current files, 220 historical references, 21 graph references
and 702 added references. Full header transport, streaming consumers and the
remaining ProForma scientific operations are ongoing work.

## Metabolite feature finding and experiment values (2026-09-11)

[Recorded checks](metabo-values-validation.json) cover the integrated feature
finder, DocumentIdentifier, five sample/instrument value types and the mzML
ASCII-subset Latin-1 fix against SDK `82ce5b3`.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,844 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Native unit/adjacent tests, Rust 1.85/no defaults | 306 passed |
| Rust 1.85/JSON and RNA selection | 134 passed |
| Rust 1.85/mzML count and adjacent selection | 70 passed |
| Strict Clippy, release library, Rustdoc and Rustfmt | Passed; warnings denied |
| Source audit, model projections, ledger and regression checks | Passed |

The batch adds 41 tests; selected totals overlap the full suite. The feature
finder reproduces the source 83/81/80 counts and all 81 expected scientific
records, including typed metadata and compressed hulls. Both fixed classifiers
match 488 separately executed LIBSVM reference cases; all 7,755 model constants
are checked by bits. Four independent peptide scores cover source f32 rounding
and the declared native f64 precision boundary.

The file example reads the unchanged original mzML and writes 81 features.
Its featureXML passes the original schema; all scalar and typed metadata fields
match the source expected output, including 1,053 numerical comparisons. This
workflow exposed an ASCII-only Latin-1 declaration compatibility gap, now fixed
with explicit non-ASCII rejection. It does not certify complete TOPP behavior.

The ledger now has 58 complete or native-equivalent headers and 728 requiring
implementation or review, with zero certified TOPP workflows. Source verification
covers 1,908 distinct current files, 220 historical references, 21 graph references
and 660 added references. The [previous checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34574946683)
passed every GitHub job. Full experiment/header transport, generic SVM support
and ProForma scientific backends remain separate work.


## Experiment metadata, feature hypotheses, ProForma JSON and mzML counts (2026-09-11)

All four groups are integrated against SDK `82ce5b3`. [Recorded checks](annotation-counts-validation.json)
retain exact commands, outcomes, log hashes, source scope and review fixes.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,803 tests passed |
| Doctests, both compiler versions | Four passed per compiler |
| Native unit/adjacent tests, Rust 1.85/no defaults | 217 passed |
| Rust 1.85/JSON and RNA selection | 122 passed |
| Rust 1.85/mzML count and adjacent selection | 55 passed |
| Strict Clippy | All targets on current/all-features and minimum/no-defaults; minimum JSON/mzML targets also passed |
| Release library, Rustdoc and Rustfmt | Passed; documentation warnings denied |
| Source audit, projection, ledger and regression checks | Passed |

The batch adds 57 tests and one lifetime doctest. Selected totals overlap the full
suite. Source metadata quirks, ordered borrowing, all reachable JSON schema branches,
and the five source count pairs are covered, with independent resource and boundary
regressions. Reviews corrected retained allocation accounting, JSON exponent handling
and sparse map allowances, as well as XML declaration/attribute validation and tiny
input-chunk handling. The original scientific fixtures remain unchanged.

Source verification covers 1,894 distinct current files, 220 historical references,
21 graph references and 601 added references. The ledger records 51 complete or
native-equivalent headers, with 735 requiring implementation or review and zero
certified TOPP workflows. FeatureFindingMetabo orchestration, ProForma scientific
backends and full experiment/header/consumer support remain separate work. No new
C++ execution is claimed for this batch. The [elution and mapping checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34572998432)
passed all GitHub CI jobs.


## ProForma text parsing and structured errors (2026-09-11)

The complete source single-chain and ion grammars, all error codes and diagnostic
operations are integrated against SDK 82ce5b3.
[Recorded checks](proforma-parser-validation.json) retain commands, outcomes,
log hashes, staged tests and the C++ extraction boundary.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,746 tests passed |
| Doctests | Three passed |
| Native unit and adjacent tests, Rust 1.85/no defaults | 195 passed |
| Strict Clippy, all targets | Both compiler/feature configurations passed |
| Release library, Rustdoc and Rustfmt | Passed; documentation warnings denied |
| Source audit, projection, ledger and regression checks | Passed |

The 16 new tests include all 198 upstream grammar cases, 62 source component calls,
structured diagnostics and cumulative resource failures. A separately compiled
[exact-source probe](../tests/data/proforma_parser_probe_provenance.json) runs 238
inputs through both C++ grammars. All 476 comparisons match: 382 accepted cases
match both text modes, and 94 errors match code, byte position and original message.
The probe includes unchanged AST/tokenizer/parser/writer and prefix-helper blocks;
a capture-only exception adapter replaces SDK exception infrastructure. It does
not execute source exception formatting or a full C++ SDK. Frozen outputs add
no C++ dependency to native builds or CI.

Source verification covers 1,886 distinct current files, 220 historical references,
21 graph references and 556 added references. ProForma remains a partial SDK header:
JSON, resolution/conversion, mass/mz and spectrum methods are outstanding in this
checkpoint. The ledger remains at 48 complete/native-equivalent headers and 738
requiring implementation or review, with zero certified TOPP workflows.
The [earlier detection/acquisition checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34571964157)
passed every CI job, including Linux, macOS, Windows and minimum Rust.


## Elution-peak detection and identification run mapping (2026-09-11)

Both operation groups are integrated against SDK 82ce5b3.
[Recorded checks](elution-mapping-validation.json) include commands, outcomes,
log hashes, review scope and the resolved integration-formatting mismatch.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,730 tests passed |
| Doctests | Three passed |
| Native unit and adjacent tests, Rust 1.85/no defaults | 179 passed |
| Strict Clippy, all targets | Both compiler/feature configurations passed |
| Release library, Rustdoc and Rustfmt | Passed; documentation warnings denied |
| Source audit, projection, ledger and regression checks | Passed |

The 24 new tests cover all six elution options, source trace splitting and
smoothing, extrema and width/noise calculations, as well as complete run/path
mapping and merged-file selection. Elution also has 49 distinct passing staged
checks per compiler, overlapping the full suite. Source numerical fixtures and
independent small-case oracles remain unchanged. Both scientific results and
input trace updates roll back on operation failure; external progress output
has a separate documented scope.

Run mapping preserves the source's deliberate duplicate-error state: complete
forward mappings remain available while reverse mappings stop before the first
collision. Resource failures preserve the prior mapping. Source and native
implementations were independently reviewed. No C++ execution is claimed for
these two groups.

Source verification covers 1,886 distinct current files, 220 historical references,
21 graph references and 548 added references. The ledger now records 48 complete
or native-equivalent headers and 738 requiring implementation or review, with
zero certified TOPP workflows.


## Detection, ProForma writing, protein runs and mzML acquisitions (2026-09-11)

All four additions are integrated against SDK 82ce5b3. [Recorded checks](detection-acquisition-validation.json)
include exact commands, outcomes, log hashes, staged validation and review fixes.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,706 tests passed |
| Doctests | Three passed |
| Native unit and adjacent tests, Rust 1.85/no defaults | 155 passed |
| Staged mzML acquisition and adjacent tests | 101 current / 100 minimum passed |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Staged mzML-only strict Clippy, Rust 1.85 | Passed |
| Release library, Rustdoc, Rustfmt | Passed; documentation warnings denied |
| Source audit, projection, ledger and regression checks | Passed |

This group adds 48 tests; selected and staged totals overlap the full suite.
The [ProForma extraction probe](../tests/data/proforma_writer_probe_provenance.json)
actually compiled the exact C++ annotation declarations and complete writer,
with no scientific substitutions or backend linkage. All 160 output cases match
native output, covering both modes, twenty precise float bit patterns and four
formatting/chain scenarios. This is not a full C++ SDK build. The other three
increments use source fixtures, branch analysis and independent native oracles.

Mass-trace detection retains source growth/termination and metadata rules, with
explicit atomic failure and reusable-state corrections. Protein-run helpers
preserve target result ownership; native lexical metadata ordering is distinguished
from C++ registry ordering. The mzML tests include real independent XSD checks
for acquisition metadata combined with zoom and scan windows in both writers.
Zoom CV placement repairs a source ordering bug; all supported parameters precede
scan windows. Historical fixture bytes and pins are retained.

Source verification covers 1,883 distinct current files, 220 historical references,
21 graph references and 533 added references. The ledger records 46 complete or
native-equivalent headers and 740 requiring implementation or review; no TOPP
workflow is certified. The [preceding published checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34570080887)
passed all CI jobs, including Linux, macOS, Windows and minimum Rust.


## Mass traces, constants, monosaccharides and mzML settings (2026-09-11)

The four additions are integrated against SDK 82ce5b3.
[Recorded checks](trace-settings-validation.json) include commands, outcomes,
log hashes, staged scope and the executed C++ constants probe.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,658 tests passed |
| Doctests | Three passed |
| Native unit and adjacent tests, Rust 1.85/no defaults | 153 passed |
| Staged mzML settings and adjacent tests | 91 passed on Rust 1.98; 90 on Rust 1.85/mzML only |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Staged mzML-only strict Clippy, Rust 1.85 | Passed |
| Release library, Rustdoc, Rustfmt | Passed; documentation warnings denied |
| Source audit, projection and ledger/regression checks | Passed |

This group adds 29 tests. Selected and staged totals overlap the full suite.
All 38 numeric declarations and 91 metadata strings match an actually compiled,
unchanged C++ Constants header. Its unused configuration include needs only an
empty shim, with no scientific declarations replaced. This is not a full SDK
build. The other three ports do not claim C++ execution.

MassTrace tests preserve source numerical/cache quirks and use independently
computed expected values where broad upstream comparison constants are stale.
Monosaccharides retain every source field, literal mass and synonym precedence.
The mzML projection preserves scientific literals and explicitly repairs its
original Product count mismatch; ordinary and Numpress writer outputs pass
independent schema tests. Acquisition fields outside this expanded representation
remain guarded before output.

The [preceding acquisition checkpoint](https://github.com/okohlbacher/openms-rs/actions/runs/34568959479)
passed every CI job. The ledger now records 45 complete/native-equivalent headers,
741 requiring implementation or review, and no certified TOPP workflow.


## Acquisition, chromatogram conversion and SDK refresh (2026-09-11)

The acquisition fields, [ChromatogramTools](CHROMATOGRAM_TOOLS_SUPPORT.md),
[processing propagation](PROCESSING_ACQUISITION_SUPPORT.md) and
[mzML write guards](MZML_ACQUISITION_GUARDS.md) are integrated against
[SDK 82ce5b3](CORE_SDK_82CE5B3_REVIEW.md). [Recorded checks](acquisition-validation.json)
include exact commands, outcomes, log hashes and separate staged scope.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,629 tests passed |
| Doctests | Three passed |
| Native acquisition and adjacent tests, Rust 1.85/no defaults | 167 passed |
| Staged mzML guards and adjacent tests | 79 passed on Rust 1.98; 78 on Rust 1.85/mzML only |
| Final signed-zero guard regression | Seven passed on each Rust version |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Staged mzML-only strict Clippy, Rust 1.85 | Passed |
| Release library, Rustdoc, Rustfmt | Passed; documentation warnings denied |
| Source audit and ledger/regression checks | Passed; 1,872 distinct current files verified |

This group adds 27 tests. The full count also includes the preceding seven-test
SequenceCoverage addition; selected totals overlap. Review covered exact grouping
and encounter order, atomic conversion failures, preservation of owned metadata
and shared processing identity, cumulative nested copy budgets, and early XML
rejection for fields that cannot yet be represented. Existing numerical
processing operations are unchanged.

Both preceding published checkpoints passed every CI job:
[5a90169](https://github.com/okohlbacher/openms-rs/actions/runs/34567597493) and
[4d2b648](https://github.com/okohlbacher/openms-rs/actions/runs/34567824647).
These are separate from the local checks for this new group. No new C++ execution
or full SDK build is claimed. The ledger records 42 complete/native-equivalent
headers, 744 requiring implementation or review, and no certified TOPP workflow.


## Standalone sequence coverage (2026-09-11)

[SequenceCoverage](SEQUENCE_COVERAGE_SUPPORT.md) is integrated. Seven new tests
cover the complete source operation, independent positional enumeration and
bounded failures. [Recorded checks](sequence-coverage-validation.json) include
82 unit/selected tests on Rust 1.98/all features and 78 on Rust 1.85/no defaults,
strict scoped Clippy on both, and the release library build. Formatting, source
and ledger checks also pass; 1,860 distinct current source files were verified.

This standalone addition leaves prior operations unchanged. Its focused totals
overlap existing tests and supplement the preceding full 1,595-test/three-doctest
checkpoint. No new C++ execution or full combined suite is claimed here. The
ledger now has 41 complete/native-equivalent headers; 745 require implementation
or review, and no TOPP workflow is certified.


## Unique IDs, 2D conversion, IMS solvers and mzML Numpress (2026-09-11)

The four operation groups are integrated. [Recorded checks](transport-values-validation.json)
include exact commands, outcomes, log hashes and separate staged validation scope.

| Check | Result |
| --- | --- |
| Full suite, Rust 1.98/all features/all targets | 1,595 tests passed |
| Doctests | Three passed |
| Native/decoy tests, Rust 1.85/no defaults | 108 passed |
| mzML Numpress and adjacent staged tests, both Rust versions | 69 passed per configuration |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Scoped mzML-only strict Clippy, Rust 1.85 | Passed |
| Release library, Rustdoc, Rustfmt | Passed; documentation warnings denied |
| Source audit and ledger/regression checks | Passed; 1,859 distinct current files verified |

This group adds 48 tests. Selected totals overlap the full suite. Independent
schema tests validate both Numpress writer variants, and 36 unchanged upstream
binary payloads decode to 342 independently checked points. No new C++ execution
or full SDK build is claimed; the earlier raw-codec probes remain separate evidence.

Reviews covered ID word consumption and native-endian UUID layout, 2D grouping
and direct integer metadata conversion, source decomposition table/order quirks,
and mzML codec/type/fallback behavior. A source zero-witness loop now returns an
error; binary validation work is bounded before mzML scalar traversal. Ordinary
fallback retains its existing numeric precision and optional zlib setting.

The ledger records 40 complete/native-equivalent headers and 746 requiring
implementation or review. Full SDK parity and TOPP readiness remain open;
no TOPP workflow is certified.


## IMS foundations, peak traversal and Numpress (2026-09-11)

The isotope/element/alphabet, area/peak-export/index, raw Numpress and configurable
base64/zlib wrapper APIs are integrated. [Recorded checks](foundation-validation.json)
include commands, outcomes, log hashes and the scope of C++ reference execution.

| Check | Result |
| --- | --- |
| New/adjacent native suites, Rust 1.98/all features | 167 tests passed |
| Corresponding native suites, Rust 1.85/no defaults | 164 tests passed |
| Wrapper/raw/options/mzML, Rust 1.98/all features | 121 tests passed |
| Wrapper/raw/options, Rust 1.85/Numpress only | 104 tests passed |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Numpress-only strict Clippy, Rust 1.85 | Passed |
| Doctests, Rustdoc, release library, Rustfmt | Passed; three doctests and documentation warnings denied |
| Source audit and review/ledger regressions | Passed; 1,850 distinct current files verified |

This additive group introduces 81 integration tests. Selected totals overlap;
they supplement the preceding full 1,466-test suite rather than establish an
additive total. The [preceding commit](https://github.com/okohlbacher/openms-rs/actions/runs/34565942433)
passed every Linux, macOS, Windows, minimum-Rust and quality CI job. A missing
cached dependency interrupted the first documentation attempt; an isolated build
completed it and all remaining checks without production changes.

Raw Numpress tests compare 295 cases against an actually compiled, unmodified
pinned C++ implementation. Of those, 287 execute its decoder; eight empty Safe
cases avoid undefined source decoding. Encoded bytes and fixed-point helpers are
exact; SLOF decoding permits a documented host-math tolerance. Wrapper transport
fixtures are independently derived Python projections, not C++ wrapper runs.
The full C++ SDK has not been built or differentially validated.

Source reviews covered isotope convolution order, portable parser boundaries,
borrowed mutable traversal, source RT grouping, and compression rejection/fallback
semantics. The ledger records 33 complete/native-equivalent headers, with 753
requiring implementation or review. mzML Numpress wiring and full TOPP readiness
remain open; no TOPP workflow is certified.

## Mobility, array descriptions, weights and mzML paths (2026-09-11)

Mobilogram operations, IMSWeights, generic array descriptions and mzML file APIs
are integrated. [Recorded checks](mobility-validation.json) include commands,
results and log hashes.

| Check | Result |
| --- | --- |
| Rust 1.98, all features and all targets | 1,466 tests passed; examples compiled |
| Rust 1.85, selected native and adjacent processing suites | 137 tests passed |
| Rust 1.85, selected XML and path suites | 72 tests passed |
| Strict Clippy, all targets | Passed on Rust 1.98/all features and Rust 1.85/no defaults |
| Selected XML Clippy, Rust 1.85 | Passed |
| Documentation examples, Rustdoc, release library, Rustfmt | Passed; three doctests and documentation warnings denied |
| Source audit and review/ledger regressions | Passed; 1,842 distinct current files verified |

Focused counts overlap the full suite. Independent source review covered mobility
sorting/search and annotation alignment, weight quantization/GCD quirks, description
preservation and XML rejection before publication, and compressed path I/O. Array
struct literals now need the description fields or `..Default::default()`; existing
constructors retain their use. Full array-description XML transport remains open.

The preceding [scientific-operations commit](https://github.com/okohlbacher/openms-rs/actions/runs/34564499936)
passed all CI jobs on Linux, macOS, Windows, minimum Rust and quality checks.
The ledger records 25 complete/native-equivalent headers and 761 still requiring
implementation or review. No TOPP workflow is certified.

## Scientific loading, decomposition and experiment operations (2026-09-11)

The native mass-decomposition solver, mzML scientific filtering and canonical
array types, experiment summaries, and idXML filesystem/dispatch operations are
integrated. [Recorded checks](scientific-operations-validation.json) include the
commands, outcomes and log hashes.

| Check | Result |
| --- | --- |
| Rust 1.98, all features and all targets | 1,432 tests passed; examples compiled |
| Rust 1.85, mzML-only plus solver/summary/value/unit suites | 147 tests passed |
| Rust 1.85, idXML-only paths/definitions plus summaries | 39 tests passed |
| Strict Clippy, all targets | Passed on Rust 1.98 with all features and Rust 1.85 without defaults |
| Focused mzML/solver/summary Clippy, Rust 1.85 | Passed |
| Documentation examples, Rustdoc, optimized library, Rustfmt | Passed; three doctests and documentation warnings denied |
| Source/fixture audit | 1,840 distinct current source/registration/reference files verified |
| Source-review and completion-ledger regressions | Three source-review tests and two ledger tests passed |

Selected test counts overlap; they are not added to the full-suite count.
Independent source review covered decomposition order and finite residue-table
limitations, mzML raw-precision filtering and aligned annotations, summary
floating-point accumulation, and idXML plain-output/extension behavior. Native
regressions cover bounded work, corrupt input and atomic failure. A preexisting
Rust 1.85 test-expression lint was corrected without changing fixture bytes.

The preceding published [extraction increment](https://github.com/okohlbacher/openms-rs/actions/runs/34561162608)
and [value-API increment](https://github.com/okohlbacher/openms-rs/actions/runs/34561620258)
passed every CI job: Linux, macOS, Windows, minimum Rust and quality checks.

The completion ledger records 23 complete/native-equivalent headers and 763
requiring implementation or review. The solver's scientific API is represented;
standalone IMS utilities are separately tracked. Full mzML metadata, codecs,
consumer/transform behavior and many other SDK APIs remain. No C++ runtime
comparison or certified TOPP workflow is claimed.


## Peak-file, metadata and composition values (2026-09-11)

PeakFileOptions, equality-compatible metadata/Product hashing, and the complete
MassDecomposition count-container API are integrated. [Recorded checks](value-apis-validation.json)
include commands, outcomes and log hashes.

| Check | Result |
| --- | --- |
| Combined new/adjacent suites, Rust 1.98, all features | 58 tests passed |
| Same applicable suites, Rust 1.85 without default features | 49 tests passed |
| Strict Clippy, all targets | Passed on Rust 1.98 with all features and Rust 1.85 without defaults |
| Doctests, optimized library, Rustfmt | Passed; three doctests |
| Source audit and completion ledger | Passed |

These selected runs supplement the preceding full 1,356-test suite; their counts
overlap and are not an additive full-suite total. Independent reviews covered
source option activation/defaults and the distinct `+` versus `+=` cached-maximum
semantics in MassDecomposition. Native hash tests record every equality-significant
field and signed-zero normalization without assuming cross-language digest values.

At this checkpoint PeakFileOptions and MassDecomposition were value APIs. The
following scientific-operations increment adds loading execution and the solver. Product's own hash
operation is now present, while inherited CVTerm/DataValue independent-unit states
and numeric registry semantics remain under review. The inventory records 22
reviewed complete/native-equivalent headers and 764 still requiring work or review.
No C++ differential execution or certified TOPP workflow is claimed.

## SDK update and extraction operations (2026-09-11)

The port now targets Core SDK 4.0.0 at `54a232f`. This increment completes the
reviewed native FASTA lifecycle and indexed mzML offset decoder, adds experiment
aggregation and XIC extraction with mzML Product interchange, and implements
peak display/hash traits. The original provenance is retained for carried-forward
fixtures; changed source files receive explicit review records.
[Recorded results](extraction-validation.json) contain the commands and log hashes.

| Check | Result |
| --- | --- |
| Rust 1.98, all features and all targets | 1,356 tests passed; examples compiled |
| Documentation examples | Three passed |
| Strict Clippy, all targets | Passed on Rust 1.98 with all features and Rust 1.85 without defaults |
| Rustdoc, optimized library, Rustfmt | Passed; documentation warnings denied |
| Indexed mzML alone, Rust 1.85 | Nine tests passed |
| Source/fixture audit | 1,834 distinct source/registration/reference files verified |
| Source-review and coverage regressions | Three source-review tests and two ledger tests passed |

All added components also passed their targeted current/minimum compiler tests
before integration. These overlapping runs are not added to the full-suite count.
Independent review closed XML lexical validation, duplicate-attribute work limits,
FASTA byte-wrapping and mzML Product metadata issues. Existing EMG regression tests
cover the private fitted-trace storage adjustment required by the larger Product
representation.

The preceding published commit `8fb47cb` passed Linux, macOS, Windows and Rust 1.85
CI, including 1,310 tests in the Linux all-target suite. Its quality job found a
new Rust 1.98 test-expression lint; the corrected expression and entire combined
crate now pass strict local Rust 1.98 checks.

The completion inventory still records 766 headers requiring implementation or
review. Twenty headers have a reviewed complete implementation or native equivalent;
this does not establish full dependency or TOPP workflow parity. No C++ executable
was built or run and no TOPP workflow is certified.

## Logging, progress and remaining map helpers (2026-09-10)

Owned logging, replaceable progress reporting, public feature/consensus
modification collection and compressed INI loading are integrated. Real local
timestamps and process CPU timing use pinned safe Rust adapters. The optimized
all-feature library builds. [Recorded results](runtime-foundations-validation.json)
include each command, outcome and log hash.

| Check | Result |
| --- | --- |
| Combined runtime, XML, definition, filesystem and unit suites, current Rust | 206 tests passed |
| Selected runtime, definition, filesystem and unit suites, Rust 1.85 without defaults | 128 tests passed |
| INI feature alone plus filesystem integration, Rust 1.85 | 32 tests passed |
| Final unit and CSV regression checks | 73 tests passed |
| Strict Clippy, all targets | Passed on current Rust with all features and Rust 1.85 without defaults |
| Doctests, Rustdoc, optimized library, Rustfmt | Passed; three doctests and documentation warnings denied |
| Source/fixture audit | 1,827 distinct current source/registration/reference files; four unchanged logging fixtures |
| Windows-style Git checkout filters | All four tested source fixtures retain exact bytes |
| Completion ledger and regressions | All 786 headers accounted for; two regression tests passed |

These selected runs overlap and supplement the complete 1,276-test suite in the
preceding increment. They are not a new complete all-target run. Reviews closed
an empty-name work-accounting gap and a logging flush-failure path that allowed
later output. Logging and progress also passed their isolated staged tests on
both compilers before integration.

CI for the preceding commit passed macOS and Rust 1.85. Windows exposed CSV
fixture checkout conversion, and the quality job exposed a Clippy test-expression
diagnostic. This increment protects all scientific fixture bytes and uses a
portable mutable byte array; local regression checks passed. Updated CI results
are tracked separately from these local results.

The completion gate still fails for 770 headers requiring implementation or
explicit review. Public modification collection and the native progress API now
have complete reviewed mappings. ZIP input, some platform logging/filesystem
behavior, full format/API coverage and executed C++/TOPP parity remain open.


## Map interchange and runtime resources (2026-09-10)

The target remains Core SDK 4.0.0 at `6bfc0e4`. Native featureXML and consensusXML,
portable modification definitions, typed feature metadata, filesystem helpers,
and shared gzip/bzip2 transport are integrated. The complete current-Rust suite
passes **1,276 tests**. [Machine-readable results](map-interchange-validation.json)
record the commands, outcomes and log hashes.

| Check | Result |
| --- | --- |
| Current Rust, all features and all targets | 1,276 tests passed; examples compiled |
| Rust 1.85, consensusXML alone plus unit/dispatch suites | 86 tests passed |
| Rust 1.85, featureXML alone plus unit/dispatch suites | 85 tests passed |
| Rust 1.85, idXML alone plus unit/custom-definition suites | 100 tests passed |
| Rust 1.85, no default features, selected unit/filesystem/definition suites | 95 tests passed |
| Strict Clippy, all targets | Passed on current Rust with all features and Rust 1.85 without default features |
| Final strict idXML-only Clippy, Rust 1.85 | Passed after narrowing a map-only helper's feature gate |
| Doctests, Rustdoc, optimized library, Rustfmt | Passed; three doctests, documentation warnings denied |
| Source/fixture audit | 1,825 distinct current source/registration/reference files; 11 batch fixtures verified |
| Coverage inventory and its regressions | All 786 headers accounted for; two regression tests passed |

Focused test counts overlap and must not be added to the full-suite count. Both
map dialects were tested against original fixtures; consensusXML output also
passed the original XSD through the locally available `xmllint`. Metadata,
protein references, portable chemistry, compressed paths and atomic failures
have independent regression coverage. Reviews caught and closed reference
expansion and validation-order resource-limit gaps.

The completion review now records ZIP input as missing from both map XML
loaders, and compressed input as missing from the INI loader. The INI entry was
therefore corrected from complete to partial. Map-based public modification
collection helpers, platform-specific filesystem behavior, inherited XML schema
validation, and broad SDK parity remain open. The completion gate still fails
for 772 headers requiring implementation or explicit review; it is not a count
of wholly absent Rust classes. No C++ executable was built or run, and no TOPP
workflow is yet certified as port-ready.


## Historical SDK foundations increment (2026-09-10)

The target remains Core SDK 4.0.0 at `6bfc0e4`. This increment adds identification
cleanup, reusable mzML parameter groups, MS2/DTA2D, file dispatch, text/list
utilities and the parameter/INI lifecycle used by TOPP tools. The final optimized
all-feature library builds successfully. Results are recorded in
[sdk-completion-validation.json](sdk-completion-validation.json).

| Check | Result |
| --- | --- |
| Full all-feature/all-target suite before the additive configuration layer | 1,128 tests passed |
| Final unit, parameter, INI, text and list suites, current Rust | 137 tests passed |
| Same selected suites, Rust 1.85 without default features | 122 tests passed |
| INI feature alone, Rust 1.85 | 14 tests passed |
| Final file dispatch and unit regression checks | 67 tests passed |
| Selected file, graph and unit suites, Rust 1.85 without default features | 99 tests passed |
| Strict Clippy, all targets | Passed with all features on current Rust and no default features on Rust 1.85 |
| Doctests, Rustdoc, optimized library, Rustfmt | Passed; three doctests, documentation warnings denied |
| Pinned source/fixture audit | 1,813 distinct source, registration and reference files verified |
| Completion inventory and its two regression tests | Passed; all 786 registered SDK headers accounted for |

These overlapping test counts are separate runs, not an additive total. The
full 1,128-test run preceded the additive parameter/text layer; focused tests
cover that layer and final whole-crate lint, documentation and release checks
cover the combined implementation. An initial inventory check found a stale
generated ledger during development; regeneration and the final check passed.

The completion gate deliberately still fails: 772 headers need implementation
or further review. Fourteen headers have a reviewed complete implementation or
native equivalent. An unmapped header can have Rust functionality that still
needs explicit review; this is not a count of wholly missing classes. The
[completion ledger](CORE_SDK_COMPLETION.md) is an inventory, not a completion
percentage. Source fixtures, schema validation and native library workflows do
not establish executed C++ differential parity or certify a complete TOPP tool.

## Historical reduced-SDK graph increment (2026-09-10)

The port now targets Core SDK 4.0.0 revision `6bfc0e4`. The optimized all-feature
library builds successfully. The complete all-feature/all-target suite passed
1,034 tests on Rust 1.96 and Rust 1.85 before this additive observation/match
layer; the no-default current run passed 905 tests before the layer was added.
The additive graph/formula checks now pass 54 private tests, 12 graph-operation
tests, 12 graph-record tests, 9 source-reference tests, 8 observation-match
tests, and 8 peptide-formula tests on the current compiler. The corresponding
focused suites pass on Rust 1.85 with no default features. Strict focused Clippy,
Rustfmt, source-target inventory checks and the 442-hash/438-link audit pass.

This increment adds observation and compound records, graph adducts, typed
molecules, observation matches, best-match queries, translated ownership and
typed peptide fragment formulas. At that stage, graph groups, cleanup/persistence
and the legacy converter were outside the port. Later increments add groups,
cleanup and the sequence/evidence conversion bridge; graph persistence and full
conversion remain outstanding.

The machine-readable results for this increment are in
[graph-validation-results.json](graph-validation-results.json).

Validated locally on **2026-09-10, macOS ARM64**. The complete all-feature suite
passes **993 tests on Rust 1.96.0: 947 integration tests across 110 suites and
46 unit tests**. The declared minimum Rust **1.85.0** passes all 129 selected
RNA and unit tests with all features. Both compilers pass 124 selected RNA and
unit tests without default features, plus three all-feature doctests each.
Full-library parity remains in progress.

The [machine-readable results](validation-results.json) retain the earlier
complete 938-test runs on both compilers and identify the current checks with
`rna_processing_` names. The new complete 993-test run covers every implemented
chemistry, processing, kernel and format suite on Rust 1.96. The minimum-compiler
and no-default-feature runs in this increment are focused checks, not complete
993-test runs.

## Historical RNA processing checks

| Check | Result |
| --- | --- |
| All features and all targets, Rust 1.96 | 993 tests passed; all fifteen examples compiled |
| All features, Rust 1.85, all RNA suites and unit tests | 129 tests passed |
| No default features, both compilers, all RNA suites and unit tests | 124 tests passed on each compiler |
| Doctests, all features | Three passed on each compiler |
| Clippy, all features/targets, Rust 1.96, warnings denied | Passed for the complete crate, tests and fifteen examples |
| Clippy, no default features/all targets, Rust 1.85, warnings denied | Passed |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| RNA processing example | Executed; RNase_T1 digestion, uridine variants and annotated negative b/y fragments |
| RNA registry regeneration | Original JSON/custom TSV reproduce all 378 embedded records and 375 distinct codes |
| RNA enzyme regeneration | Original XML reproduces all fourteen embedded enzymes and 84 field values |
| Package file list | 428 files; no package build or publication; RNA data release terms tracked separately |
| Current provenance and documentation | 421 source/fixture hashes and 413 local documentation links verified; RNA processing source lines, registry fields and ion bits checked alongside all retained scientific fixtures |

This increment adds fifty integration tests in five suites: twelve RNase,
fourteen modification-generation, twelve spectrum-generation, nine independent
source-reference and three workflow tests. Five new private tests verify the
shared formula work and allocation allowance, including the spectrum generator's
formula-only precursor path. The selected matrix includes all earlier RNA
record/provider, sequence, source-reference and workflow tests to check the
shared sequence helper refactor on both compilers.

Independent references cover fourteen enzyme records, 6,048 pattern/code cases,
twelve original digestion cases with 38 products, all four modification counts
(7, 6, 27 and 432), and all 126 fragment masses actually compared in the source.
The workflows connect modified cleavage, positional products and variable
modifications to annotated spectra and both mzML compression modes. Tests also
cover retained metadata, source charge/sulfur conventions, resource exhaustion
and atomic failures. No expected scientific value was changed to fit an
implementation result. See [RNase support](RNASE_SUPPORT.md),
[RNA modification generation](RNA_MODIFICATION_SUPPORT.md),
[RNA spectra](RNA_SPECTRUM_SUPPORT.md) and the
[independent review](RNA_PROCESSING_REFERENCE_REVIEW.md).

The checks below are retained historical evidence, with their original scope.

## Retained complete RNA foundation checks

| Check | Result |
| --- | --- |
| All features and all targets, Rust 1.96 and 1.85 | 938 tests passed on each compiler: 897 integration tests across 105 suites and 41 unit tests |
| No default features and all targets, Rust 1.96 | 863 tests passed |
| No default features, Rust 1.85, original RNA suites and unit tests | 70 tests passed |
| JSON-only feature, Rust 1.96, original RNA suites and unit tests | 73 tests passed |
| Doctests, all features | Three passed on each compiler |
| Clippy, all features/targets, warnings denied | Passed for the complete crate, tests and fourteen examples |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| RNA example | Executed; charged formula/mass/m/z, sulfur-aware suffix and coarse isotope probabilities |
| RNA registry regeneration | All 378 embedded records and 375 distinct codes reproduced |
| Package file list | 403 files at that stage; no package build or publication |
| Provenance and documentation at that stage | 359 source/fixture hashes and 382 local documentation links verified |

The foundation added thirteen record/provider tests, nine sequence tests,
eight independent source-reference tests and three workflows. These cover all
fifteen original formula assertions, eighteen mono/five average mass assertions,
ten positive slice examples and every registry entry. The JSON dependency
exposed two test-only empty-array type-inference ambiguities, resolved with
equivalent `is_empty()` assertions. A narrow precision-lint allowance preserves
the original carbon-13 mass literal in one workflow; the affected workflows
were rerun on both compilers. See [RNA support](RNA_SUPPORT.md) and the
[foundation reference review](RNA_REFERENCE_REVIEW.md).

## Earlier Tagger-specific checks

| Check | Result |
| --- | --- |
| All features, Rust 1.96 and 1.85 | 58 tests passed on each compiler: all 41 unit tests plus 17 Tagger integration tests |
| No default features, Rust 1.96 and 1.85 | 57 tests passed on each compiler; the mzML workflow excluded |
| Clippy, all features/targets, warnings denied | Passed for the final crate, tests and thirteen examples |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| Tag extraction example | Executed; five measured peaks produce the exact sorted tags EP, EPT, PE, PEP, PEPT and PT |
| Package file list | 380 files, including Tagger, its example, source fixtures and documentation; no package build or publication |
| Provenance and documentation at that stage | 293 source/fixture hashes and 353 local documentation links verified; six Tagger count rows and 120 membership rows verified at their original source lines alongside retained scientific fixtures |

The independent suite reproduces all six original Tagger counts and all 120
membership assertions using the source's theoretical-spectrum settings. The
source input sizes, 357 and 180 peaks, are checked first. Derived tests isolate
strict mass-window boundaries, ties, exact collisions, modified I/L behavior,
signed/unsorted coordinates and zero/inverted settings. Twelve private mass-table
tests cover the 19 base residues, two-stage resolution, deterministic provider
order, empty short IDs and source free-residue arithmetic.

Three private traversal/sort tests include a 9,999-residue path without recursion
and failures after partial traversal that preserve old string allocations. The
workflows connect measured tags to exact target/decoy substring matching,
modified digestion fragments and both mzML compression modes. See
[Tagger support](TAGGER_SUPPORT.md) and its [independent review](TAGGER_REFERENCE_REVIEW.md).

## Retained decoy and adduct checks before Tagger

| Check | Result |
| --- | --- |
| All features, Rust 1.96 and 1.85 | 58 tests passed on each compiler: all 26 unit tests plus 32 new integration tests |
| No default features, Rust 1.96 and 1.85 | 56 tests passed on each compiler; both XML workflows excluded |
| idXML without mzML, Rust 1.96 | 57 tests passed; decoy idXML workflow included |
| Clippy, all features/targets, warnings denied | Passed for the final crate, tests and twelve examples |
| Rustdoc and Rustfmt | Both passed; documentation warnings denied |
| Decoy FASTA example | Executed; two targets and two deterministic variants per target, six records total |
| Package file list | 369 files, including full Boost license notices and the reproducible RNG oracle; no package build or publication |
| Provenance and documentation at that stage | 277 source/fixture hashes and 336 local documentation links verified; new source/Boost hashes and RNG oracle checked alongside all retained scientific fixtures |

Decoy reference tests replay all thirteen unique original sequence cases,
including their shared random/cache history. Nine private RNG tests cover four
complete cycles for three seeds, a published standard value, rejection draws,
shuffle order and state handling. Three additional private tests cover identity
and transactional failures. The packaged integer oracle regenerates independent
checkpoints, checksums and permutations without running C++.

Adduct tests preserve the original scalar examples, then distinguish their
rounded tolerance from exact atomic/electron arithmetic. They check parser
edge cases, signed formulas, finite negative masses and native limits. The
workflows compare sodium-adduct isotope masses with complete ion composition,
preserve charge/metadata through mzML, and connect decoy FASTA to indexing,
synthetic FDR values and idXML. See [decoy support](DECOY_GENERATION_SUPPORT.md),
[adduct support](ADDUCT_SUPPORT.md) and their linked reference reviews.

## Retained full baseline before decoy and adduct additions

| Check | Result |
| --- | --- |
| All features and all targets, Rust 1.96 | 815 integration tests and 14 unit tests passed; eleven examples compiled |
| Documentation examples/API encapsulation | 3 doctests passed |
| No default features | 747 integration tests, 14 unit tests and 3 doctests passed; mzML/idXML excluded |
| All features and all targets, Rust 1.85 | Same 815 integration tests and 14 unit tests passed |
| Rust 1.85 documentation tests | 3 passed |
| idXML without mzML | 772 integration tests, 14 unit tests and 3 doctests passed; no base64/zlib feature dependency |
| Clippy, all features/targets, warnings denied | Passed |
| Rustfmt check | Passed |
| Rustdoc, warnings denied | Passed |
| Spectrum processing example | 121 input peaks → 14 retained peaks; normalized TIC 1.000000 |
| FASTA digestion example | Executed; positional peptide/mass/m/z TSV produced |
| Modified peptide analysis example | Executed; masses, five-bin isotope envelope and annotated b/y fragments produced; invalid input independently checked |
| Synthetic identification example | Executed; selected AC(Carbamidomethyl)DMK with seven fragment matches and 26.32% protein coverage; precursor/evidence/metadata checked |
| Chromatogram integration example | Executed; two synthetic peaks, original sample boundaries, raw intensity sums, time-weighted areas and baseline estimates produced |
| EMG fitting example | Executed; 36 cropped samples → 80 fitted samples, 700 iterations; sampled area 2,121.476414 → 2,501.460284 versus known complete synthetic area 2,506.628275 |
| Profile processing example | Executed; 401 profile samples → four iterative centroids → two retained peaks with aligned integration and width annotations |
| Modified peptide enumeration example | Executed; two digest products → four alkylated/oxidized variants plus one unchanged product; masses and doubly charged m/z values printed |
| Isotope streaming/enrichment example | Executed; first five natural glucose configurations and two enriched-carbon configurations above absolute threshold 0.1, with raw/log probabilities |
| Peptide-property example | Executed; modified and unmodified peptides produce finite charge, pI, GRAVY and gas basicity at both 500 K and 100 K |
| Spectrum annotation example | Executed; eleven source IFSQVGK measurements annotated with source b/y labels, exact widened intensity 12.100000262260437, longest y-series six and finite matching statistics |
| mzML input processing example | Executed using the independent mixed-precision fixture |
| Independent XML XSD checks | Populated/empty mzML, native idXML, auxiliary-array and populated-precursor mzML writer output validated with xmllint and their pinned schemas |
| Generated modification data | All 3,035 specificity rows regenerate exactly; original XML matches pinned archive byte for byte |
| Generated enzyme data | All 33 source records and 602,316 independent Python regex contexts regenerate exactly; bundled XML is byte-identical to the pinned source |
| Package file list | 351 files; spectrum annotation/ion naming, shared-work regressions, peptide-property tables/numerical fixtures, original data, schemas, isotope streaming/custom inputs and fine isotope/precursor purity integration, theoretical extensions, XLMOD/OBO and prior scientific fixtures, examples and licenses included; local reference snapshot and target directory excluded |
| Source/fixture provenance and documentation | 252 source/fixture hashes and 307 local documentation links verified; 41 spectrum-annotation rows and source/derived bits, 633 peptide-property rows and their binary64 values, isotope streaming/custom inputs, fine isotope/precursor purity and theoretical-spectrum sources, all 148 XLMOD mass literals/IEEE values, byte-identical ontology, prior scientific fixtures and libm archive/lock checksum checked; earlier chromatogram/schema checks retained |

A constructor-level source review during validation corrected `precursor_in_ms2`
to integer 0/1 and the spectrum ppm flag to the string projection `"0"`/`"1"`.
The 21 affected annotation/reference/workflow tests were rerun on both compilers
with all features; the optional-feature matrix, lint, documentation and examples
completed after the correction.

The unchanged modification/enzyme generators and data retain their successful
checks from earlier on the same date; they were not rerun for these API edits.

The package-list check does not publish or build a package. Cargo reported only
optional missing documentation/homepage/repository metadata for this local port.

## Test coverage

| Suite | Tests | What it establishes |
| --- | ---: | --- |
| RNase records, registry and digestion | 12 | Fourteen source enzymes, complete record identity, replacement/alias semantics, modified-code cleavage, end gains, product order/coordinates and atomic limits |
| RNA modification generation | 14 | Fixed/variable selection, terminal-site and maximum-one source quirks, alternative identity/order, exact preflight counts/bytes and atomic failures |
| RNA spectrum generation | 12 | Nine ion series, single/multiple charge branches, source mass/sulfur/precursor conventions, aligned annotations and shared resource limits |
| Independent RNA processing references | 9 | All fourteen enzymes, 6,048 predicates, 38 source digestion products, modification counts/identities and 126 compared source ion literals |
| RNA processing workflows | 3 | Modified cleavage, positional digests and variant mass shifts, annotated spectra and both mzML compression modes |
| RNA formula and spectrum budget internals | 5 | Shared work/allocation depletion, preserved mass bits and input state, empty-sequence behavior and cumulative formula-only precursors |
| RNA records and providers | 13 | Complete fields, identity, ordered duplicate registry, ambiguity resolution, JSON/TSV differences, diagnostics and checked limits |
| RNA sequences | 9 | Source parser/display, custom ownership, terminal and sulfur slicing, all fragment variants, natural-H/electron mass convention and atomic failures |
| Independent RNA references | 8 | Fifteen original formulas, 23 mono/average mass references, ten slices, 378-record registry projection and source branch/identity cases |
| RNA workflows | 3 | Independent charged formulas, custom carbon-13 registry lifetime, coarse/fine isotopes and both mzML compression modes |
| Sequence tags | 8 | Defaults/setter, exact append ownership/order, finite source option edges, ignored spectrum fields, registry lifetime and atomic resource failures |
| Independent Tagger references | 6 | Six exact source counts, 120 membership assertions, source peak counts/rounded traces, strict bounds, nearest ties, collisions, I/L, signed/unsorted and no-op behavior |
| Tagger mass-table internals | 12 | Nineteen free-residue formulas, fixed/variable order, exact keys, two-pass provider resolution, terminal/wildcard/anonymous records, formula/mass precedence and construction limits |
| Tagger traversal/sort internals | 3 | A 9,999-edge heap traversal, failure after emitted paths, string-byte sort work and allocation-preserving atomic errors, strict lookup/ties |
| Tagger workflows | 3 | Exact substring target/decoy matching, modified digest fragment gaps and both mzML compression modes |
| Decoy generation | 7 | All registered enzymes, reversal/shuffle source rules, composition where applicable, modified/empty inputs, cache/resource boundaries and atomic failures |
| Independent decoy references | 7 | All thirteen literal source cases, ordered RNG/cache history, reseeding, cross-context reuse, zero attempts, short products and unspecific order |
| Decoy RNG internals | 9 | Three seeds over four cycles, standard 10000th word, rejection mapping, shuffle order, clone/reseed and no-draw boundaries |
| Decoy identity/state internals | 3 | Forward/reverse maximum identity, strict ties and late cache/work/allocation rollback |
| Decoy workflows | 3 | FASTA, target/decoy evidence and synthetic FDR, composition through digestion and idXML transport |
| Molecular adducts | 8 | Complete parser/components/getters/equality, electron and n-mer conversions, mono/average shifts, compatibility and checked limits |
| Independent adduct references | 5 | Literal source cases, tighter atomic constants, isotope/charge/whitespace grammar, signed containment and source numerical boundaries |
| Adduct workflows | 2 | Complete sodium-ion composition, unchanged isotope probabilities, charge-dependent spacing and both mzML compression modes |
| Chemistry | 12 | Source-derived mass/formula and digestion values, residue/isotope tables, formula grammar/algebra limits, b/y mass conservation |
| Charge and isoelectric point | 10 | Four source pKa scales, terminal overrides/suppression, parent PTMs, U/O/ambiguous residues, saturated finite pH, exact endpoint/midpoint behavior and resource/convergence guards |
| pI internal work accounting | 1 | Endpoint and midpoint evaluations consume one precharged budget |
| Amino-acid indices and gas basicity | 7 | Source accessions/indicators, ordinary split arithmetic, independent low-temperature values, tied maxima, empty/high-temperature identity, formula-free annotations and input limits |
| Hydrophobicity profiles | 7 | Seven scales, literal GRAVY/window/moment values, sliding order, all annotation types, window/angle semantics and preallocated work limits |
| Independent property references | 7 | 503 source constants including 42 sentinels, isolated pKas, 30 source scalar assertions, 420 singleton/pair GB expressions, 100 Decimal-derived extreme-temperature cases and independent moments |
| Peptide-property workflows | 4 | Digestion retains only original terminal caps, mass-only chemistry needs no formula, source generator attachment distinctions and typed property metadata through idXML |
| Ion naming | 9 | Every source charge/ordinal case, first-line and caret priority, field/overflow fallbacks, about 18,000 round trips, Unicode text, output limits and allocation-free parsing |
| Spectrum annotation | 10 | Three source operations, array replacement, ppm last-match/duplicate branches, no-op flags, final sorting, safe small-list statistics, finite source special cases, scoped validation and atomic failures |
| Independent annotation references | 7 | Original measured bits and 17 literal assertions, independently derived fragment errors, exact binary top-N padding/sample variance/ratios, source label grammars and precursor absolute tolerance in ppm mode |
| Annotation workflows | 4 | Modified digestion against independent fragment formulas, distinct chemistry with identical sequence text, peak annotations/statistics through idXML and aligned arrays/acquisition through both mzML compression modes |
| Shared generation/alignment work | 3 | Actual alignment initialization/cells and shared residue/loss/fine/coarse allowances survive successive calls; standalone calls start fresh budgets |
| Annotation grammar internals | 1 | Complete source label regex semantics, including commas, losses, ordinal boundaries and terminal-series distinctions |
| Sequence chemistry and numeric tags | 14 | Source integer/decimal registry lookup, absolute/internal/H/OH masses, unresolved residue representation, formula/mass availability, owned tags and shared immutable spelling through independently owned slices, stable attachment, slicing, atomic setters and numeric parser limits |
| Sequence identification operations | 8 | B/Z/X non-mass filters, known-mass precursor filtering, exact modification IDs, owned protein observations, distinct resolver keys and atomic failures |
| Sequence chemistry workflow | 8 | Ambiguous digestion/indexing, numeric-tag idXML round trips, independent fragment/loss shifts, formula-dependent append rejection and extreme-mass cancellation regressions |
| Modifications | 11 | Full pinned registry load, specificity/name/mass lookup, neutral losses, isotope labels, terminal notation, atomic setters, modified digestion/fragments |
| Modified peptide generation | 12 | Fixed/variable source goldens, stable alternative order, terminal overwrite/duplicates, existing anonymous annotations, bounded combinations, atomic errors and custom formula-free mass/absolute/no-op rules |
| Independent modification generation | 9 | Seven fixed and nineteen variable source cases, full weighted-site order, exact typed terminal/fragment-mass distinctions, append limits, empty idXML rejection and independent custom mass branches |
| Modification definitions | 13 | Owned definition identity, fixed/variable set semantics, compatibility, count non-enforcement, absolute/delta matching, inference, conflicting full-ID chemistry and atomic limits |
| Independent definition reference | 5 | Seven literal compatibility cases, partition/merged precedence, negative-delta tolerance endpoints, stored/fallback absolute masses and named/anonymous all-hit inference |
| Modified peptide workflow | 5 | Digestion → variants → formulas and independent fragment shifts → inferred search definitions → idXML; mass-only fragments, isotope rejection and serialization preflight |
| Owned modification registry | 7 | Caller record validation, optional vocabulary IDs, shared handles, bounded atomic OBO appends, formula fields and complete Eq/Ord contracts including signed zero |
| Independent OBO registry reference | 10 | All 148 XLMOD records in exact source order with identities, specificities, synonyms and mass bits; PSI aliases and absent targets, stanza/EOF handling, literal empty versus zero formula and alias-work limits |
| Crosslink lookup | 6 | Separate 56-record registry, source DSS/BS3/EDC mass/site goldens, reactive-side union, terminal conversion, search eligibility and bounded caller-owned extensions |
| Anonymous modification definitions | 6 | Owned exact-spelling annotations, all-hit inference, source empty-short-ID compatibility, full-residue/H/OH mass anchors, unresolved deltas and atomic errors |
| Caller-owned chemistry workflow | 7 | Registry lifetime release and sharing, numeric/name resolution, absolute-formula X rescue, no-change/delta/terminal precedence, non-UniMod mass export, full-ID-only records and exact custom-registry idXML |
| Custom chemical identity | 6 | Same text with distinct formulas/vocabularies remains distinct in protein observations, owned peptide keys, sequence duplicate filtering, rank/spectrum conflicts; atomic late errors |
| Internal theoretical fragments | 7 | Source interval and ten-residue boundaries, ordered annotations, numeric residue and terminal chemistry, loss replacement, isotope independence, intensities, atomic append and custom declaration/storage limits |
| Compact mass-only spectra | 7 | Full-length six-series ladders, source terminal/float operation order, observed masses, independent suffix accumulation, sorted no-ops and combined/precision limits with late-error atomicity |
| Activation presets and immonium ions | 6 | Source mass goldens, all activation enum cases, inferred precursor metadata, unmodified residue eligibility and L-only branch, charge/intensity/isotope independence and aligned atomic append |
| Independent theoretical extension review | 9 | Analytical internal interval/mass/loss oracles, exact immonium constants, source rounded CID/partial ECD tables and independently rounded f32 compact-helper values |
| Scalar purity and SPS matching | 9 | Source isolation totals/residuals, inclusive doubled tolerances, charge/nearest conventions, f32 fragment windows, input immutability and shared work/unused-annotation limits |
| Fuzzy and interpolated purity | 11 | Neutron-spacing source successor lookup, f32 sums/division, strict and half-weight boundaries, ratios above one, zero-window/empty-parent behavior, RT extrapolation/fallbacks and nonfinite-field checks |
| Independent precursor purity review | 13 | Exact decoded source peaks and seven scalar/map goldens, SPS float32 bounds and counts, independent fuzzy arithmetic and parent reference/acquisition-order rules |
| Precursor metadata and interchange workflow | 10 | Shared native precursor ownership, referenced parents, scalar batches, original complete mzML fixture, all activation/mobility quantities, malformed XML/loss guards and populated acquisition XSD validation |
| Native fine isotope enumeration | 8 | Small exhaustive configurations, fixed labels/natural zero abundances, charge adducts, source rounding/coverage, deterministic thresholds and checked atom/output/frontier resource limits |
| Independent fine isotope references | 7 | All 44 source counts, 14 fructose mass/probability rows, 6 bromine configurations, subnormal/full-support tails and the 19,615-state insulin f32 coverage boundary |
| Independent fine theoretical spectra | 7 | Literal 10/5/50/12 source counts and exact mass/intensity tables, independent CHNO oracles, neutral H/charge division, terminal/loss rules and single/shared-budget atomic failures |
| Fine isotope workflow | 3 | Modified peptide formulas and charge adducts, aligned spectrum selection, unknown-composition errors and both compressed/uncompressed mzML round trips with precursor acquisition |
| Fine isotope internals | 2 | Deterministic equal-probability state identity and cumulative work allowance across calls |
| Fine isotope stream and custom inputs | 8 | Materializer equivalence, original-mode thresholds, charge-ignoring raw adapter, owned lifetimes, fused errors, checked conversion, zero-count validation and large-support prefixes |
| Independent stream references | 6 | Original 14-row/2548-state fructose and 10000-state insulin ordered cases, threshold counts, independently reconstructed natural tables, custom f64 precision and direct multinomial products |
| Fine isotope stream boundary review | 10 | Extreme finite weights, log underflow, late mass overflow, getter-consistent cutoff equality and adjacent values, cumulative custom dimensions, deferred expansion and 103041-state full streaming |
| Custom isotope spectrum workflow | 3 | Owned raw-coverage prefixes and continuation, distinct equal-mass configurations through materialization, aligned spectrum selection and both mzML compression modes |
| Digestion registry and specificity | 12 | All 33 enzyme predicates vs independent regex hashes, full/semi/nonspecific order and small-sequence oracle, validity/missed cleavages, modifications and resource limits |
| Typed metadata and acquisition | 17 | Scalar/list/unit invariants, explicit string bridges, CV merging, source settings defaults/unify, drift/isolation conventions, date/checksum validation and atomic errors |
| Identification records | 18 | Stable scores, typed evidence/metadata, charge ranges, grouping, interval coverage vs independent oracle, observed modifications, kernel attachments and file loss rejection |
| Score categories and switching | 9 | All 29 names/six categories, source lookup/backup conventions, restoration, heterogeneous records, reserved labels and atomic map/slice errors |
| Peptide-spectrum scoring | 8 | Source HyperScore/Morpheus goldens, overload-specific precision/charge/boundary rules, ion annotations/ordinals, error means and checked failures |
| Identification filters | 16 | Source cutoff/top-N/dense-rank/tie behavior, modifications, exact/sequence duplicates, signed charge/precursor errors, run-aware references/groups and map atomicity |
| FDR and q-values | 13 | OMSSA1534 and XTandem source thresholds, picked goldens, actual legacy/Basic formulas, score ties/directions, peptide representatives, posterior estimates/ROC, reserved labels and atomic limits |
| idXML | 16 | Exact pinned fixtures, native metadata/run/evidence round trips, independent XSD validation, malformed XML/reserved encodings, byte/element/list guards and preflight/flush errors |
| Identification processing workflow | 1 | Independent XML → score switching/top-hit selection → hand-computable target/decoy q-values → threshold/reference cleanup → coverage → exact native idXML round trip |
| Peptide-to-protein indexing | 17 | 121 matcher/17 decoy/15 enzyme oracle cases, ambiguity/mismatch/I-L rules, original-engine recovery, per-run evidence/protein reconstruction, metadata placement and preflight/atomic limits |
| Basic protein inference | 21 | Source merger scores/groups, best/product/mean, charge/modification representatives, greedy graph ties/negative scores, score restoration, vector/single/consensus conventions and atomic limits |
| Identification conflicts | 11 | Source rank aggregation and spectrum reports, best/matching selection, modified-sequence/charge keys, intensity winners, metadata/subordinate retention and atomic map errors |
| Identification origin partitions | 7 | Three origin formats, source ordering/path reduction/hit union, basename/key collisions, run-specific protein values, skipped records, group omission and checked errors |
| Combined protein workflow | 1 | FASTA → modified/shared peptide evidence → basic inference → distinct protein/PSM q-values → global coverage → native idXML round trip → file-origin partitions |
| Retention-time models | 12 | Source weighted/unweighted fits, interpolation/extrapolation goldens, LOWESS sine/cars/original tables, inverses, diagnostics/windows and resource/atomic errors |
| Container RT transformations | 7 | Source experiment/feature/consensus/peptide values, original metadata, nested hulls/subordinates/handles, optional spectrum IDs, ordering and late-error atomicity |
| Isotope distributions | 13 | Coarse glucose/heavy/bromine goldens, convolution, enrichment, peptide/RNA/DNA averagine, exact Poisson recurrence rounding and overflow fallback, conditional fragment probabilities, limits |
| Theoretical spectra | 12 | All ion-series mass tables, losses, precursors, isotope conventions, modified peptides, aligned atomic append and resource limits |
| Kernel | 19 | Source search/TIC values, ties/ranges, metadata alignment, mutation checks and malformed inputs |
| Feature/consensus geometry | 18 | Scan-envelope hulls, source containment quirks, consensus means/charges/decharge, map identities/ranges, subordinates and checked mutations |
| Processing | 14 | Original 121-peak fixture, normalization, source 0.05 threshold default and adjacent-f32 cutoff, top-N/rank behavior, atomic errors and resampling conservation/boundaries |
| Smoothing | 10 | Gaussian goldens/trapezoids/ppm behavior; Savitzky–Golay goldens, asymmetric edges and polynomial reproduction |
| Baseline morphology | 8 | All ten source operations, boundary windows, signed bottom-hat, width conversion and aligned atomic updates |
| Peak picking | 12 | Orbitrap/FTMS centroid goldens, full historical noise values, natural spline values/derivatives, FWHM, mobility, boundaries and experiment selection |
| Iterative peak picking | 8 | HiRes seed/refinement noise separation, source centroid and integration conventions, width/spacing rules, aligned annotations, selected experiments and atomic limits |
| Internal iterative conventions | 2 | Fourteen independent refinement configurations and deterministic ties; strict-next-sample association, asymmetric recenter search, original seed priority and rounded centroid storage |
| Independent iterative reference | 4 | Complete synthetic result, original Orbitrap/FTMS raw-region sums and centroid rounding, omitted annotations, experiment behavior and checked failures |
| Iterative profile workflow | 2 | Noise estimation → centroiding → local peak filtering, independent raw sums, preserved metadata and exact mzML width/intensity array round trip |
| Window filtering | 10 | Source 56/30 counts, exact triangle indices, 480 independent cases, sliding/jumping boundaries and order, duplicate membership, stable ties and atomic resource failures |
| Iterative mean noise | 7 | All 2,526 historical source outputs, hand-computed fixed-denominator clipping, strict window edge, histogram ceiling, signed/empty inputs, checked legacy percentile and exact work limits |
| Independent mean-noise reference | 3 | Fresh histogram rescans for 192 configurations, signed global statistics, effective bin ceiling and legacy-percentile f32 arithmetic |
| Chromatogram picking | 9 | Legacy/corrected methods, Gauss/SG smoothing, independent seed/boundary noise, raw sums, exact index regions, aligned arrays, overlap and checked failures |
| Internal chromatogram conventions | 2 | Source right-side closest-sample ties/end sentinel and sequential overlap midpoint assignments |
| Peak integration | 12 | Inclusive sampled sum/trapezoid/Simpson, all baseline choices, signed/empty/singleton inputs, shape metrics, duplicates, bounds and numeric/resource failures |
| Independent chromatogram processing | 12 | All 146 source trace samples and literal integration/background/shape/picker goldens, nonuniform/even Simpson behavior, f32 arithmetic, -1 sentinel exclusion and source noise/overlap conventions |
| Chromatogram workflow | 2 | Raw picking → exact full-precision boundaries → intensity sum/time-weighted area/baseline/shape → mzML peak-annotation round trip; f32 boundary rounding cannot change integration |
| EMG fitting | 8 | Source cutoff fit, three model branches, metadata and omitted arrays, zero bounds, best iterations, extrapolation, shared evaluation limits and numerical failures |
| Internal EMG conventions | 3 | All four analytic gradients against central differences, source training collection order and iRprop sign/zero/rollback rules |
| Independent EMG reference | 8 | Seven complete source fits and four full raw loss goldens, true relative parameter tolerances, ordered training losses, scalar tails/subnormals, branch boundaries and inclusive budgets |
| EMG integration | 6 | Both container types, all integration/baseline methods, unchanged supplied shape inputs, cropped and expanded spans, literal zero bounds and preserved inputs on errors |
| EMG workflow | 2 | Independent continuous area/centroid identities and cropped fit → mzML → sampled integration without unaligned parameter annotations |
| Deisotoping | 8 | Unknown/unequal precursor-charge source regressions, charge priority, ladders, intensities, annotations, disjoint membership and error atomicity |
| Poisson/KL deisotoping | 11 | Source defaults and threshold/top-N preprocessing, original indices and aligned arrays, both sharing policies on both algorithms, low-mass noise, bounded work and transactional errors |
| Complete deisotoping fixture | 3 | All 5,407 source input peaks retained; all 104 expected output peaks bit-exact, independent seed/charge mapping, shared membership, exact 103-peak disjoint subset and retained metadata/arrays |
| Independent deisotoping review | 8 | Adjacent-f32 KL thresholds for sizes 2–7, mixed-precision accumulation, longest/highest-charge selection, selected counts/sums, inclusive ppm endpoints, nearest ties and exact precursor arithmetic on both methods |
| Spectrum comparison | 19 | Source alignment/scoring/binning goldens; 1,000 compact-DP/reference-map comparisons; directed ppm and weighted rounding boundaries |
| Text formats | 12 | DTA/FASTA/MGF parsing, streaming/round trips, invalid input, writer preflight/flush errors and malformed-input corpus |
| mzML | 14 | Independent/upstream fixtures, codecs/precision/units, metadata subset, limits, malformed XML/arrays and independent XSD checks |
| Independent format/processing review | 5 | Regressions for nonfinite filtering, empty resampling, full-range charges and metadata key loss |
| Independent mzML review | 3 | Conflicting scientific CV fields, reserved names and forbidden XML characters |
| Independent mzML auxiliary reader | 10 | Both float/signed integer widths, exact integer range, empty ASCII elements/placeholders, malformed/unsupported metadata, variable string decompression and cumulative byte/element/array limits |
| Independent mzML auxiliary writer | 4 | Validation before first output, exact source binary encodings, zero-length distinctions, type/name order and metadata, all four compression/empty-record XSD cases |
| Synthetic identification workflow | 1 | Embedded FASTA → fixed modification → digestion → literal observed fragment matching → peptide evidence/annotations → protein coverage/modifications; precursor consistency and retained metadata |
| Workflows | 4 | FASTA → digestion/fragments → MGF; DTA → processing → MGF; modified peptide → isotope envelopes → deisotoping → alignment; profile → centroid → consensus |

Spectrum annotation and ion naming received separate source, numerical and
integration reviews. The four new integration suites contain 30 tests, with one
additional sequence-ownership regression and four private grammar/work tests.
Original source fixtures contain thirteen measured peaks and seventeen scalar
or label assertions. Eleven independently derived masses/errors distinguish
rounded source references, including its zero-MSE assertion, from actual
nonzero error statistics. Exact binary examples independently check padding,
sample deviation, quartile indices and current ratios. See the
[annotation support](SPECTRUM_ANNOTATION_SUPPORT.md),
[reference review](SPECTRUM_ANNOTATION_REFERENCE_REVIEW.md) and
[fixture provenance](../tests/data/spectrum_annotation_provenance.json).

Native annotations preserve source last-match array behavior, duplicate
matched-only records, parameter no-ops, sort order and raw absolute precursor
tolerance in ppm mode. Small-list quartiles and undefined enabled statistics
have documented safe policies. Workflows verify modified fragments from
independent formulas, complete custom-chemistry identity, and XML transport.

Shared theoretical and alignment counters now persist across candidate hits;
coarse convolution also retains one allowance across theoretical envelopes.
Focused exhaustion tests verify that calls cannot restart these counters.
Anonymous modification strings are shared immutably, with independent slices
and setters retaining their owned lifetime, value and ordering semantics.

Peptide physicochemical utilities received independent source-table and numerical
reviews. Five integration suites add 35 tests, and a private pI test verifies its
shared evaluation budget. All 200 public AAindex values, 182 hydrophobicity cells
and 59 pKa constants are checked. Private gas-basicity tables feed independent
source expressions for every canonical singleton and all 400 ordered pairs.
Thirty original scalar assertions retain their source spelling and tolerance
context. A separate 100-case grid uses 90-digit Decimal arithmetic, including
low temperatures, tied maxima, subnormal products and the largest finite inputs.
See the [property support](PEPTIDE_PROPERTIES_SUPPORT.md),
[reference review](PEPTIDE_PROPERTIES_REFERENCE_REVIEW.md) and
[fixture provenance](../tests/data/peptide_properties_provenance.json).

Ordinary gas basicity preserves source evaluation order. Stable overflow retries,
the underflowed-product limit and the empty-sequence identity are documented
native numerical differences. Terminal annotation tests preserve source pI
suppression without inventing PTM-specific pKas; hydrophobicity and gas basicity
continue to use parent residues. Digestion and identification interchange tests
exercise those rules together. The current complete Rust 1.96 run includes all
110 integration suites and 46 private unit tests. Three doctests pass on both
compilers, alongside the focused RNA processing checks detailed above.

Owning isotope streaming and custom binary64 populations received independent
source and numerical reviews. The source ordered fructose support (2,548 states),
fourteen literal configurations and insulin's 10,000-state prefix pass alongside
threshold and custom natural-table materialization checks. Custom weights retain
extra binary64 precision; hand-derived products preserve separate equal-mass
configurations. Raw formula charge is ignored while the high-level generator's
hydrogen-adduct rules and all prior fine-spectrum references remain unchanged.

Boundary tests distinguish finite logarithms from underflowed raw probabilities,
check late overflow and fused errors, and validate every custom cell even for
zero-count populations. A 10,000-category case succeeds when coverage stops at
the first configuration, then fails with a checked work error when its raw stream
requests more. Raising its threshold avoids expanding an ineligible suffix.
A separate stream drains 103,041 configurations although materialization rejects
that full support. Threshold tests retain a returned raw probability exactly and
reject the next larger binary64 value; materialized logarithmic cutoff behavior
is preserved. Distinct equal-mass peaks and their configuration indices survive
selection and both mzML compression modes.

The materializer boundary assertion checks the analytic logarithmic condition
rather than assuming a particular platform's `ln`/`exp` rounding direction. Both
compilers pass it alongside the raw-value equality assertions in the full matrix.

The existing search is shared without a new dependency. Iterator cutoffs use
returned raw values, with logarithmic fallbacks for relative underflow or an
unrepresentable mode. Source layers, performance hints and untrimmed layer-based
membership remain outside the current surface. Streaming retains a frontier and
visited set, so atom, work, state and memory limits bound its lifetime.

Native fine-isotope enumeration received independent source, mathematical and
theoretical-spectrum integration reviews. All 44 selected class-test count
assertions, the 14-row high-precision fructose table and six bromine configurations
are checked separately from older loose/header illustrations. The 19,615-state
insulin boundary verifies accumulation of stored f32 probabilities; full-support
tests retain configurations whose stored probabilities underflow to zero. A small
exhaustive combinatorial oracle checks masses and probabilities independently of
the native heap algorithm, including separate fixed labels and natural H adducts.

The theoretical integration reproduces the source's 10-, 5-, 50- and 12-peak
fine-spectrum cases and checks retained terminal chemistry, source loss formulas,
isotope annotation alignment and natural-H mass/charge conventions. Loss intensity
products now remain f64 until final storage, avoiding an intermediate f32 rounding
in both fine and coarse envelopes; existing coarse regressions still pass. A
105-prefix test fails specifically at the shared fine-work limit while each
single envelope is valid and total possible output is only 105 peaks. Both this
failure and late atom-limit errors leave appended observations unchanged. Modified
fine spectra preserve their annotations and precursor acquisition data through
both mzML compression modes.

This validates materialized native configurations with explicit deterministic
ties, inclusive threshold equality and checked resource limits. It does not
claim identical IsoSpec layering or untrimmed layer-dependent selections, or
an executed C++ differential comparison. No dependency was added.

Precursor purity received independent numerical, source and acquisition-format
reviews. Exact binary64 coordinates and binary32 intensities from the original
mzML support seven source scalar/map cases; eight SPS cases and independent fuzzy
arithmetic are also checked. Rounded upstream literals are kept separate from
independently calculated results. The full original fixture reproduces the same
scores through the native mzML reader after an explicit in-memory encoding
declaration change; the bundled source bytes remain unchanged. All 19 activation
methods, four mobility quantities, distinct isolation targets and spectrum
references are covered, including populated writer output checked by xmllint.

Review found full spectrum validation repeatedly traversed unused annotations
outside the purity work budget. Purity now validates only consumed numerical
fields; a 10,000-child/10,000-placeholder regression verifies this bound. Scoped
precursor validation remains explicit before fuzzy empty-parent shortcuts.
Malformed container nesting and loss of CV-list ordinary metadata on export are
also rejected. The writer rejects a negative effective isolation target before
output, including selected-m/z fallback; spectrum/chromatogram regression cases
prevent producing a document that the reader would reject. Native source-undefined iterator/arithmetic cases return checked
errors or documented end-candidate fallbacks; finite source overcount and RT
extrapolation are retained. No C++ binary was executed for these comparisons.

Theoretical-spectrum extensions received independent source and integration
reviews. Analytical references cover internal fragment masses, counts, source
start/end conventions and first-residue loss omission. Literal source immonium
masses and rounded activation tables are checked separately; the compact helper
has an independent bit-exact f32 analytical oracle, rather than captured C++
output. Coarse-envelope independence, custom mass tags, retained termini,
annotation alignment and error atomicity are covered across the focused suites.

Review found repeated custom loss declarations could consume unbounded work
before deduplication and distinct losses could fill templates before peak
preflight. Declaration visits now receive a shared preflight estimate, and
unique loss-template entries are charged before formula cloning. Regression
tests verify both limits for ordinary and internal generation while preserving
existing output. Fine isotope support was added subsequently and is covered above.

Owned modification records, OBO loading and crosslink lookup received independent
source and consumer reviews. All 148 projected XLMOD records match exact source
accession/site order, names, synonym sets and f64 mass bits. Synthetic PSI-style
cases verify all-specificity UniMod aliases, absent targets and absolute formula
precedence. The historical PSI-MOD snapshot remains unbundled; its source hash
is verified separately from packaged fixtures.

Review corrected source empty-sequence export, literal empty versus explicit
zero-formula behavior and repeated alias expansion work accounting. A further
consumer review found same-name custom chemistry being merged in protein
observations, sequence duplicate filtering and conflict keys. Complete chemical
value ordering now preserves those distinctions and agrees with equality,
including signed zero. Peptide identity keys retain owned sequence values.
Source-intentional textual keys in other identification algorithms remain.
Caller-owned registry tests verify record release after the last owner is dropped,
absolute-formula restoration of unknown residues, checked vocabulary/mass export,
and exact idXML chemistry validation before output.

Modified-peptide generation and definition sets received independent source and
integration reviews. Literal cases include 26, 71 and 199 variants, supplemented
by a complete small weighted-site oracle. Source maximum-one/general terminal
placement and formula-derived versus declared mass distinctions are preserved.
Review found and corrected inference merging of chemically different records
sharing a full ID; conflicting records now fail atomically in either order.
Custom formula-free records retain declared/absolute monoisotopic mass through
generation and fragment calculations, with explicit unavailable composition.
The source no-change rule and formula precedence are checked separately.
Normal generated variants and search definitions round-trip through idXML;
unrepresentable typed terminal states are rejected before any writer call.

Iterative peak picking received an independent source and numerical review.
The source class test has no numerical assertions; fourteen derived refinement
configurations explicitly cover source control flow and mixed precision. Two
real input profiles independently verify centroid rounding and raw-region sums;
they are not treated as iterative output goldens. Window selection matches
literal upstream counts and 480 independent small cases. Iterative mean noise
matches all 2,526 historical outputs at the source tolerance and fresh histogram
rescans over 192 configurations. The combined profile workflow verifies that
filtering preserves the picker’s aligned integration and width arrays through
mzML. The implemented source quirks and native checked-error policies are
recorded in the corresponding support documents.

EMG fitting received an independent numerical review. All twelve branch
expressions across four gradients match the source arithmetic, with separate
central-difference checks. Seven source fits retain their literal minute/second
coordinates, f32 container intensities and true relative parameter checks;
four source loss assertions distinguish full raw-f64 loss from training loss.
Derived training-order fixtures and erfc tail/subnormal tests independently
check operation ordering and special-function behavior. The strict numerical
error policy is documented, including failures after a previous finite best.
An additional mathematical workflow verifies total EMG area and centroid, then
checks cropped reconstruction and exact fitted-trace mzML interchange.

Independent reviews inspected the numerical modules and their source conventions.
Chromatogram integration and picking preserve literal source reference values,
mixed f32/f64 operation order, sampled boundaries and shape conventions. Reviews
caught an uncharged Gaussian coefficient-table allocation and verified the
corrected work limit. The combined workflow additionally preserves exact original
sample indices when f32 boundary metadata cannot represent their coordinates.

Auxiliary mzML arrays received independent reader and writer reviews. Tests
verify binary bytes, signed integer limits, ASCII terminators and empty elements,
validation before any write/flush, and cumulative decoded-storage bounds. A
review regression rejects units attached to encoding/compression terms, preventing
unrepresented metadata loss. The provenance audit distinguished the mzML schema
fixture's CRLF bytes from the immutable source's LF bytes and verified identical
normalized schema content; both hashes are now recorded.

Poisson/KL deisotoping received separate numerical and data-flow reviews. The
complete real-spectrum fixture revealed two source clusters sharing a heavy
isotope; the native default preserves both, while an explicit disjoint policy
omits only the independently identified overlapping seed. Source binary arrays,
decimal fixture values and all 104 seed/charge mappings were independently
verified. Review also corrected precursor operation order and low-mass noise
handling, and checked mixed f32/f64 KL acceptance boundaries.

Sequence review checks the distinction between known mass and known formula,
annotation ownership, position-preserving round trips, and source numeric lookup
precision. Regressions cover two large-mass cancellation errors: independently
accumulated suffix fragments preserve a small alanine fragment beside a large
tag, and equal/opposite terminal deltas preserve the peptide residue mass.

Protein inference, peptide indexing, conflict resolution, origin splitting and
the combined protein workflow received independent scientific reviews. The
matcher oracle includes 121 raw cases, 17 decoy-inference cases and 15 enzyme
boundaries. Review regressions address undefined-negative greedy scores,
original search-engine recovery, search-parameter metadata placement and
pre-clone sequence limits.

Identification scoring/filtering/FDR and idXML received separate peer reviews.
Regressions reject target/decoy label replacement, aliased analysis-result
indices, malformed declarations/processing instructions and encoded group lists
over their configured limit. Other reviews checked digestion context, typed
metadata, evidence coordinates, transformation inverses and atomic container
edits. The synthetic identification precursor was corrected against an
independent mass calculation.

All 84 elements and 283 isotope data points were compared against the source
declarations. All 2,341 inventory entries were verified against immutable Git
objects at the pinned revision; one CRLF checkout conversion is recorded
separately. Checks establish the documented surface; they do not count fully
ported C++ classes.

## Limits of this evidence

The C++ reference was inspected, **not built or executed**. These are
source-derived reference tests, independently reconstructed numerical
expectations, fixture/schema checks and Rust invariants. They are not runtime
differential tests against OpenMS or proof of full workflow equivalence.

Linux, macOS, Windows and Rust 1.85 CI have passed for the preceding published
commit, as recorded above; newer local changes have separate validation records.
Linux CI installs xmllint; on other hosts the schema tests explicitly report when
that optional executable is unavailable. It was available and executed locally.
Rust 1.85 compatibility was also exercised locally.

There are no runtime performance/memory benchmarks, external application
interoperability certification or full PSI controlled-vocabulary validation.
The capability documents describe current exclusions and deliberate source
corrections. No C++ source, contrib tree or vendored dependency was modified.
