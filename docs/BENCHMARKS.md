# Benchmarks: this port against the C++ OpenMS4 Release build

What is compared, how, and what the numbers do and do not say. The current
figures are the **wave-4 run of 2026-09-16**: all eight ported TOPP tools, on
full-size instrument data, at 1 and 32 threads, five (or three, or ten)
repetitions per cell on a quiet node, with the equivalence of every output
judged separately for data and for metadata. They replace the wave-3 results
entirely. The port is faster on two tools and level on a third at one thread,
and 1.9x faster on the picker at 32 threads; it is slower on the rest and uses
more memory on most. Section 6 lists what the run does **not** establish, and
section 5 the caveats that qualify every figure in it.
Section 4 is a separate, narrower wave-5 run that answers one build-flag
question and is **not** comparable with the wave-4 tables in absolute terms.

Everything below is against **one** C++ reference: the optimised build at
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`. The
product SDK used as the correctness oracle elsewhere in this repository is a
*Debug* build of a different core revision and is never used for timing.

## 1. The C++ reference build

| | |
|---|---|
| Prefix | `/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` |
| Pins | core `bc9cc12`, cli `c19e494`, topp `174b576` (all three tree hashes match the local checkouts) |
| Compiler | conda-forge gcc 14.4.0, default target `-march=x86-64 -mtune=generic` |
| Flags | `CMAKE_BUILD_TYPE=Release`, `CMAKE_CXX_FLAGS` empty, so `-O3 -DNDEBUG` on 795 of 795 core, 7 of 7 cli and 122 of 122 topp translation units; `-mssse3` and `-ffp-contract=off` on 743 core units; `-fopenmp` throughout. No `-g`, `-O0`, `-march`, `_GLIBCXX_ASSERTIONS` or sanitizer. `config.h` has `#if (0) #define OPENMS_ASSERTIONS`, so `OPENMS_PRECONDITION` and `OPENMS_POSTCONDITION` expand to nothing |
| Manifest | `BUILD_MANIFEST.json` in the prefix, sha256 `2efe7501c5c7e3f6dd3feef78628ae6f1158e762b7fc04154de21cfd095b293a`; the same file is registered in `SOURCE_PROVENANCE.json` as `../oracle/release-build/manifest.json` |
| Libraries | binaries carry `DT_RPATH` (not RUNPATH), so staging must merge `deps/lib` into the staged `lib/`; `ldd` of the staged copy then resolves 0 libraries under `/ceph` |

It was checked against the recorded Debug oracle before being used: all 186
C1-ORACLE-TOPP cases were rerun, 184 matched, and the 2 that did not are
exactly the two the manifest flags `debug_only` (they exit 8 on a
precondition in Debug and 0 with an empty featureXML in Release — the Release
behaviour is the reference, **never** the Debug exit code). Of 143 compared
outputs, 27 are byte-identical to Debug, 34 identical after path substitution,
82 FuzzyDiff-equal and 0 different. The 82 FeatureFinderCentroided featureXML
files differ from Debug only in last floating-point digits (largest relative
difference 9.06e-11), so **bit-exact featureXML expectations taken from the
macOS Debug oracle do not hold against this build**; use FuzzyDiff or
regenerate.

## 2. The harness

`/ceph/ibmi/abi/oliver/bench/openms4/harness` — bash and the Python standard
library only, documented in its own `README.md`, every file hashed in
`run.json`.

- `bench.py` runs cases (tool × dataset × implementation × thread count) with
  one warm-up and R repetitions, default 5. Inputs and the C++ prefix are
  staged to node-local `/scratch`; outputs go there too; `TMPDIR` and
  `OPENMS_TMPDIR` are on `/dev/shm`. Both implementations get `-threads N`
  **and** `OMP_NUM_THREADS = RAYON_NUM_THREADS = N`.
- Measurement is `os.wait4` **through a small launcher**, so peak RSS is the
  tool's own: `/bin/true` through the same launcher records 1,024 KiB, and the
  naive `Popen`+`wait4` path records the harness's own high-water mark instead.
  Recorded per repetition: wall, user, sys, peak RSS, page faults, I/O, exit
  code or signal, and the output sha256. Every run records a `/bin/true` floor
  before and after the cases.
- **Load gate.** Before each case, a 2 s instantaneous sample of CPU-busy time
  and run queue per core; it waits above 0.25, refuses above 1.0 and flags in
  between. Foreign CPU per core is measured during every repetition; above 1.0
  the repetition is invalid and re-run. Load averages before and after are
  recorded as well — load average alone is a poor gate (under a real 0.14/core
  foreign load, `loadavg1` showed 0.067/core).
- `summarize.py` reports the **median** and IQR of wall, user, sys and peak
  RSS, and the rust/cpp ratio of medians with a 95 % percentile bootstrap CI
  (10,000 resamples, n ≥ 3). No means.
- `equiv.py` walks a ladder: byte-equal, then the C++ `FuzzyDiff` (ratio
  1.000001, absdiff 1e-9, C1-style whitelist), then a decoded comparison that
  base64/zlib-decodes the binary arrays, and otherwise DIFFERENT with the path
  of the first difference. Non-verdicts (`not_decidable`, `missing_output`,
  `no_counterpart`, `error`) are always recorded, never dropped.
- `node.json` records `lscpu`, `lscpu -e`, `numactl`, governor counts, boost,
  `amd_pstate`, kernel and THP per run.

**Inputs** are staged under `/ceph/ibmi/abi/oliver/bench/openms4/inputs`, 7.4 GB,
with `inputs/MANIFEST.json` recording source path, size, sha256, spectrum or
feature counts, per-MS-level centroid/profile counts from a full streaming
pass, instrument and compression. Nothing was downloaded; everything came from
ABI storage. The classes are: profile high-resolution (Q Exactive
`UK222.mzML`, 2.21 GB; LTQ Orbitrap XL E. coli, 1.44 GB), centroided LC-MS
(PXD001819 Orbitrap Velos `50amol_R1.mzML`, 1.14 GB; `UK222_picked.mzML`,
522 MB), featureXML (59.6 MB and 2.1 GB) and three small sanity files.

**Parameters.** Both implementations read the *same* INI, produced by the C++
tool's own `-write_ini` and hashed in the run record. That is the only way the
two are known to be doing the same work; where a plan runs a tool without
`-ini`, each side uses its own defaults and the comparison is not sound.

**Timing node.** `ibminode05`: AMD EPYC 7763, 1 socket, 64 cores × 2 SMT = 128
logical CPUs, 1 NUMA node, 1019 GB, governor `schedutil` with boost on, kernel
6.8, Ubuntu 24.04, `/scratch` on NVMe. One socket, so NUMA binding is
unnecessary. The node is **not** exclusive and its CPU frequency is **not**
pinned; what that did and did not cost the wave-4 run is in caveats 3 and 9.

## 3. The wave-4 run: eight tools, 1 and 32 threads

Run directory `/ceph/ibmi/abi/oliver/bench/openms4/results/2026-09-16-w4threads`,
finalised 2026-09-16T04:00:04+00:00. It replaces the wave-3 results, which are
superseded throughout and were, in two places, measuring something other than
what they claimed: DTAExtractor's apparent win compared two outputs of different
sizes, and MapNormalizer's ratios compared two programs computing different
answers. §3.7 is the tool-by-tool comparison.

All eight ported TOPP tools are measured on **full-size instrument data**. The
three tools whose Rust build refused real input in wave 3 — MzMLSplitter,
SpectraFilterWindowMower and FileInfo on featureXML — now accept it, so none of
wave 3's reduced substitute inputs is used. FeatureFinderCentroided is the one
exception and runs on the same documented 4,000-spectrum subset as wave 3; §5.6
says why, and what that costs the conclusion.

### 3.1 What ran

| item | value |
|---|---|
| Rust | `openms-rs` `9a392fee4aefd4c2a1d6baa022a6b81f8a0493e5`, 26 commits after wave 3's `fabd4b9`. `cargo build --release --locked --offline --bins`, default features (which include `parallel`), rustc 1.96.0, no `RUSTFLAGS`, x86-64 baseline, no `[profile.release]` override. Source tarball sha256 `4142c674969ca4a9…`; build manifest at `/scratch/kohlbach/bench/openms4/rust/9a392fee4aef/BUILD_MANIFEST.json` |
| C++ | the §1 Release build, byte-identical to the binary wave 3 used (sha256 verified per tool across the two runs), staged node-local, `ldd_probe.resolved_under_ceph` = 0 |
| node | ibminode05 only. Timing only; nothing was built there for the C++ side |
| parameters | one INI per tool, from the C++ `-write_ini`, one sha256 per tool across implementations and repetitions. Three of four INIs the reviewer diffed are byte-identical to a fresh `-write_ini`; FileInfo carries three symmetric documented edits (`m`/`p`/`s` = true) |
| command | `-ini <shared> -in … -out … -threads N -no_progress`, with `RAYON_NUM_THREADS = OMP_NUM_THREADS = N` |
| order | implementations interleaved within each round in a seeded random order (seed 20260916) |
| repetitions | 1 discarded warm-up + 5 measured rounds; cells whose warm-up is under 10 s are raised to 10; SpectraFilterWindowMower uses 3, because one execution takes 8.6 to 11.8 minutes |
| load | gate value 0.0163–0.0215 per core over the 36 per-cell samples; foreign CPU during the measured repetitions 0.0092–0.0184 per core, median 0.0168. **0 of 192 repetitions load-flagged, 0 refused, 0 retried, 0 failed**; `majflt` 0 in all 192 |
| peak-RSS floor | `/bin/true` through the launcher recorded 1,024 KiB in all four sub-runs, before and after the cases. The naive `Popen`+`wait4` method recorded 24,576–43,104 KiB in the same runs, i.e. the harness's own high-water mark — so the caveat that wave 1 raised is fixed and verified, and the numbers below are the tools' own |
| self-test | 28 harness tests, 0 failed |

The run and its arithmetic were **independently re-measured**: every median,
IQR, ratio, bootstrap CI, start-up-corrected value, peak RSS, thread count and
output-byte figure recomputes exactly from the 192 raw repetitions, and an
independent re-run of eight cells with a separate runner and a separate `/proc`
sampler landed within 1.80 % on the ratio (worst PeakPickerHiRes at 32 threads,
0.521 against 0.531). **No number in either table needed correcting.** Three of
the runner's *explanations* did, and those corrections are applied below.

### 3.2 One thread

Wall is the median of the measured repetitions with the interquartile range.
Ratio is rust/cpp (**< 1 = Rust faster**) with a 10,000-resample percentile
bootstrap 95 % CI of the ratio of medians. The start-up-corrected column
subtracts each side's own 10-record-slice median from both sides. Peak RSS is
GNU `time %M` of the tool itself, median over repetitions.

| tool | dataset | n (r,c) | rust wall med s [IQR] | cpp wall med s [IQR] | ratio rust/cpp | 95 % CI | start-up corr. | rust RSS MiB | cpp RSS MiB | RSS ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| DTAExtractor | Velos centroid 1.2 GB | 5,5 | 28.913 [28.896–29.003] | 26.318 [26.314–26.357] | 1.099 | [1.095, 1.102] | 1.100 | 1691.0 | 1516.0 | 1.12 |
| MzMLSplitter | Velos centroid 1.2 GB | 5,5 | 18.506 [18.472–18.520] | 14.522 [14.514–14.526] | 1.274 | [1.271, 1.280] | 1.278 | 4063.7 | 1524.0 | 2.67 |
| BaselineFilter | QE profile 2.3 GB | 5,5 | 24.964 [24.911–24.987] | 15.809 [15.762–15.833] | 1.579 | [1.557, 1.591] | 1.584 | 6600.0 | 3199.0 | 2.06 |
| MapNormalizer | Velos centroid 1.2 GB | 5,5 | 16.671 [16.655–16.716] | 14.749 [14.728–14.755] | 1.130 | [1.128, 1.136] | 1.133 | 1778.8 | 1518.0 | 1.17 |
| SpectraFilterWindowMower | Velos centroid 1.2 GB | 3,3 | 515.513 [515.469–527.383] | 705.732 [705.391–705.791] | **0.73–0.77** | [0.730, 0.765] | 0.73–0.77 | 3346.8 | 1516.0 | 2.21 |
| PeakPickerHiRes | QE profile 2.3 GB | 5,5 | 25.367 [25.326–25.420] | 25.114 [25.096–25.146] | 1.010 | [1.005, 1.013] | 1.011 | 3371.0 | 3882.8 | 0.87 |
| FileInfo | Velos centroid 1.2 GB | 5,5 | 12.807 [12.752–12.815] | 15.951 [15.944–15.967] | **0.803** | [0.797, 0.806] | 0.803 | 2427.9 | 2009.0 | 1.21 |
| FileInfo | featureXML 60 MB | 10,10 | 1.493 [1.486–1.497] | 1.005 [0.999–1.008] | 1.486 | [1.475, 1.498] | 1.569 | 151.0 | 79.0 | 1.91 |
| FeatureFinderCentroided | Velos 4,000-spectrum subset | 5,5 | 112.380 [112.264–112.389] | 90.086 [90.003–90.104] | 1.247 | [1.244, 1.250] | 1.249 | 310.3 | 360.1 | 0.86 |

The SpectraFilterWindowMower ratio is printed as a range on purpose: n = 3 and
one repetition is 4.6 % high (§5.5).

### 3.3 Thirty-two threads

| tool | dataset | n (r,c) | rust wall med s [IQR] | cpp wall med s [IQR] | ratio rust/cpp | 95 % CI | start-up corr. | rust RSS MiB | cpp RSS MiB | RSS ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| DTAExtractor | Velos centroid 1.2 GB | 5,5 | 28.992 [28.907–29.007] | 21.752 [21.725–21.800] | 1.333 | [1.309, 1.349] | 1.336 | 1691.0 | 1525.0 | 1.11 |
| MzMLSplitter | Velos centroid 1.2 GB | 5,5 | 18.562 [18.462–18.584] | 9.818 [9.817–9.820] | 1.891 | [1.876, 1.905] | 1.906 | 4063.8 | 1530.8 | 2.65 |
| BaselineFilter | QE profile 2.3 GB | 5,5 | 25.240 [25.175–25.629] | 14.406 [14.342–14.474] | 1.752 | [1.731, 1.793] | 1.762 | 6600.0 | 3198.7 | 2.06 |
| MapNormalizer | Velos centroid 1.2 GB | 5,5 | 16.682 [16.669–16.800] | 9.983 [9.953–9.989] | 1.671 | [1.653, 1.688] | 1.684 | 1778.9 | 1527.4 | 1.16 |
| SpectraFilterWindowMower | Velos centroid 1.2 GB | 3,3 | 516.581 [516.378–540.840] | 640.095 [639.899–640.981] | **0.80–0.88** | [0.804, 0.883] | 0.80–0.88 | 3347.0 | 1530.8 | 2.19 |
| PeakPickerHiRes | QE profile 2.3 GB | 5,5 | 12.592 [12.540–12.594] | 23.691 [23.647–23.734] | **0.531** | [0.527, 0.535] | 0.530 | 3367.0 | 3899.6 | 0.86 |
| FileInfo | Velos centroid 1.2 GB | 5,5 | 12.850 [12.807–12.877] | 11.351 [11.335–11.356] | 1.132 | [1.125, 1.137] | 1.137 | 2427.9 | 2016.7 | 1.20 |
| FileInfo | featureXML 60 MB | 10,10 | 1.492 [1.486–1.495] | 1.046 [1.043–1.056] | 1.427 | [1.409, 1.433] | 1.520 | 151.0 | 84.0 | 1.80 |
| FeatureFinderCentroided | Velos 4,000-spectrum subset | 5,5 | 26.772 [26.762–26.787] | 25.249 [25.240–25.272] | 1.060 | [1.057, 1.062] | 1.065 | 306.4 | 378.9 | 0.81 |

### 3.4 Where the port is faster, where it is slower

**Faster.**

- `SpectraFilterWindowMower`, both thread counts (roughly 0.73–0.77 at 1 thread
  and 0.80–0.88 at 32), on the full 1.2 GB Velos run. This is a real CPU win —
  512 s against 704 s of *user* time for bitwise-identical peak data — and it is
  the first time this comparison exists at all, because wave 3's Rust build
  refused the input.
- `FileInfo` on the 1.2 GB mzML at 1 thread, 0.803. The cleanest result in the
  run: same input, byte-identical 2,831-byte report, nothing written,
  reproduced independently at exactly 0.803, start-up 0.4–0.6 % of either side.
- `PeakPickerHiRes` at 32 threads, 0.531. The port's picking loop runs on 32
  workers and turns 25.367 s into 12.592 s (2.01×) while the C++ picker gains
  1.06×. A genuine parallelism result, not a start-up or work-volume artefact.

**Level.** `PeakPickerHiRes` at 1 thread, 1.010, CI [1.005, 1.013] — inside the
~3 % cross-session band of §5.4, so read it as *no difference measured*, not as
1 % slower.

**Slower.** Everything else, at both thread counts: `DTAExtractor` 1.099 /
1.333, `MapNormalizer` 1.130 / 1.671, `MzMLSplitter` 1.274 / 1.891,
`BaselineFilter` 1.579 / 1.752, `FileInfo` on the 60 MB featureXML 1.486 /
1.427, `FeatureFinderCentroided` 1.247 / 1.060, and `FileInfo` on the mzML at 32
threads 1.132 (the 1-thread win is gone).

**Memory.** The port uses more resident memory on six of the nine cases — worst
MzMLSplitter at 4.06 GB against 1.52 GB (2.67×), then BaselineFilter 6.60 GB
against 3.20 GB (2.06×) and SpectraFilterWindowMower 3.35 GB against 1.52 GB
(2.21×) — and less on PeakPickerHiRes (0.87×/0.86×) and FeatureFinderCentroided
(0.86×/0.81×). Peak RSS is essentially thread-count-independent on both sides.
There is no memory parity and none is claimed.

### 3.5 Thread behaviour, sampled rather than assumed

`peak thr` is the maximum thread count in `/proc/<tool>/task`, sampled every
5 ms. `CPU util` is (user + sys)/wall, i.e. the average number of CPUs the
process actually kept busy. A tool that asked for 32 threads and shows util
~1.0 did **not** use them, however many threads exist. The `busy thr` column
wave 3 printed is deliberately absent: its 100 ms per-thread sampling period is
longer than several of these cells, and the wave-3 review showed it reporting
1 busy thread beside a CPU utilisation of 5.86.

| tool | dataset | impl | t1 peak thr | t1 util | t32 peak thr | t32 util | t32 user s | t32/t1 wall | verdict at 32 |
|---|---|---|---|---|---|---|---|---|---|
| DTAExtractor | Velos 1.2 GB | rust | 2 | 1.00 | 33 | 1.00 | 25.9 | 1.00× | **stayed serial** |
| DTAExtractor | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | 6.12 | 130.3 | 1.21× | parallel |
| MzMLSplitter | Velos 1.2 GB | rust | 2 | 0.99 | 33 | 0.99 | 13.8 | 1.00× | **stayed serial** |
| MzMLSplitter | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | 12.38 | 118.7 | 1.48× | parallel |
| BaselineFilter | QE 2.3 GB | rust | 2 | 1.00 | 33 | 1.00 | 17.7 | 0.99× | **stayed serial** |
| BaselineFilter | QE 2.3 GB | cpp | 2 | 1.00 | 64 | 9.03 | 122.5 | 1.10× | parallel |
| MapNormalizer | Velos 1.2 GB | rust | 2 | 1.00 | 33 | 1.00 | 13.6 | 1.00× | **stayed serial** |
| MapNormalizer | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | 12.22 | 119.1 | 1.48× | parallel |
| SpectraFilterWindowMower | Velos 1.2 GB | rust | 2 | 1.00 | 33 | 1.00 | 513.6 | 1.00× | **stayed serial** |
| SpectraFilterWindowMower | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | 1.18 | 750.0 | 1.10× | parallel |
| PeakPickerHiRes | QE 2.3 GB | rust | 1 | 1.00 | 33 | 2.26 | 24.7 | 2.01× | parallel |
| PeakPickerHiRes | QE 2.3 GB | cpp | 2 | 1.00 | 64 | 5.82 | 132.9 | 1.06× | parallel |
| FileInfo | Velos 1.2 GB | rust | 1 | 1.00 | 1 | 1.00 | 11.1 | 1.00× | **no pool built** |
| FileInfo | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | 10.85 | 121.4 | 1.41× | parallel |
| FileInfo | featureXML 60 MB | rust | 1 | 1.00 | 1 | 1.00 | 1.4 | 1.00× | **no pool built** |
| FileInfo | featureXML 60 MB | cpp | 2 | 1.00 | 33 | 4.25 | 4.4 | 0.96× | parallel |
| FeatureFinderCentroided | Velos 4,000 subset | rust | 1 | 1.00 | 33 | 4.72 | 125.8 | 4.20× | parallel |
| FeatureFinderCentroided | Velos 4,000 subset | cpp | 2 | 1.00 | 64 | 4.22 | 106.2 | 3.57× | parallel |

**The Rust side uses its pool in two tools of eight.** PeakPickerHiRes builds a
33-thread pool (32 workers plus main) and reaches util 2.26 with a 2.01× wall
speed-up — the picking loop parallelises well, but the serial read and write
around it bound the total, and the pool is alive in only 282 of 2,632 samples.
FeatureFinderCentroided reaches util 4.72 and 4.20×. Five tools (BaselineFilter,
DTAExtractor, MapNormalizer, MzMLSplitter, SpectraFilterWindowMower) build the
33-thread pool from the flag and leave it idle: util 1.00, user time unchanged
from 1 thread. FileInfo builds no pool at all. At `-threads 1` the Rust tools
peak at 1 or 2 threads, so the pool is only constructed when more than one
thread is asked for. **This is the largest structural asymmetry in the run**:
on those five tools every 32-thread ratio measures a serial program against a
parallel one.

**The C++ side's extra threads are a library pool, not extra compute.** Every
C++ cell at `-threads 32` peaks at 64 threads except FileInfo on the featureXML,
which peaks at 33, and every cell at `-threads 1` peaks at 2 — but that surplus
is not an OpenMP team doing tool work, and it is wrong to read the 32-thread
column as "C++ was given twice the budget":

- at `-threads 1` the second thread is `jemalloc_bg_thd`, which accumulates no
  CPU at all — utilisation is exactly 1.00;
- on PeakPickerHiRes at `-threads 32` there are 63 tool-named threads but the
  whole process accumulates 136.0 s of thread CPU over 25.5 s of wall
  (util ~5.3), i.e. about five cores of work, not 64. That figure is one tool's:
  the per-cell utilisations in the table above range from 1.18
  (SpectraFilterWindowMower) to 12.38 (MzMLSplitter), so no single number
  describes the C++ column — what holds for every cell is only that the peak
  thread count is not the compute budget;
- the surplus tracks `OMP_NUM_THREADS`, which the harness sets for the C++ side,
  not the tool: with the variable unset the same binary peaks at **129** threads
  at `-threads 1` and **160** at `-threads 32`, with wall unchanged. Those are
  the OpenBLAS server pool sized at library load — the same mechanism this port
  documents at `src/cli/context.rs:164-168`.

So the tool bodies get comparable thread budgets. What differs is that the port
uses its pool in two tools of eight.

**High C++ CPU at 32 threads is largely libgomp spin-wait, not work.** The
clearest case is SpectraFilterWindowMower: 64 threads and 750 s of user CPU buy
66 s of wall (705.7 → 640.1 s). BaselineFilter spends 122.5 s of user CPU for a
1.10× gain and PeakPickerHiRes 132.9 s for 1.06×. Read `user s` next to wall
time; wall alone flatters the C++ 32-thread column.

### 3.6 Output equivalence, data and metadata judged separately

Neither comparison stops at the first difference; differences are counted per
category. Tolerances are per tool (`tools.json`): mzML default 0.001 ppm m/z and
1e-6 relative intensity; PeakPickerHiRes 0.01 ppm / 1e-5; FeatureFinderCentroided
0.01 ppm / 0.01 s / 1e-4 relative with a 10 ppm, 5 s match window.

| tool | dataset | data | metadata | what differs |
|---|---|---|---|---|
| DTAExtractor | Velos 1.2 GB | bitwise_equal | n/a | 36,443 files, 36,443 byte-equal, 836,505,793 B on both sides |
| MzMLSplitter | Velos 1.2 GB | equal_within_tolerance | **DIFFERENT** | 87,492 arrays over four parts, all bitwise; 18 metadata differences in 14 categories on part 1, 13 in 9 on the others |
| BaselineFilter | QE 2.3 GB | equal_within_tolerance | **DIFFERENT** | 81,713 of 81,714 arrays bitwise, 1 within tolerance (one `f32` ULP, ≤ 5.96e-8); 14 differences in 12 categories |
| MapNormalizer | Velos 1.2 GB | equal_within_tolerance | **DIFFERENT** | 87,492 of 87,492 arrays bitwise; 20 differences in 14 categories |
| SpectraFilterWindowMower | Velos 1.2 GB | equal_within_tolerance | **DIFFERENT** | 87,492 of 87,492 arrays bitwise over 43,745 spectrum pairs, 0 peak-count mismatches; 22 differences in 16 categories |
| PeakPickerHiRes | QE 2.3 GB | equal_within_tolerance | **DIFFERENT** | 81,714 of 81,714 arrays bitwise; 16 differences in 14 categories |
| FileInfo | Velos 1.2 GB | bitwise_equal | n/a | 1 file, byte-equal (2,831 B) |
| FileInfo | featureXML 60 MB | bitwise_equal | n/a | 1 file, byte-equal (1,949 B) |
| FeatureFinderCentroided | Velos 4,000 subset | equal_within_tolerance | **equal** | 4,076 of 4,076 features matched, 0 unmatched either way, 0 charge disagreements, 0.0 ppm m/z, 0.0 s RT, largest relative intensity difference 9.477e-8 |

The verdicts are identical at 1 and at 32 threads, and were reproduced
independently with a separate decoder, one verdict per tool.

**Every tool agrees on the data**, and both of wave 3's data findings are
closed: MapNormalizer's 7,302 intensity arrays differing by ~13× are gone (all
87,492 now bitwise), and DTAExtractor's 21.2 % smaller output at reduced
precision is gone (both sides write exactly 836,505,793 bytes across 36,443
byte-equal files).

**Metadata still differs on the five mzML writers**, in a narrow class:
identifier spelling and processing provenance, never numbers. The recurring
items are the `dataProcessing` / `software` / `sourceFile` id strings (Rust
`dp_00000000000000000000`, `so_dp_…`, `sf_…` against C++ `dp_sp_0`,
`so_dp_sp_0_pm_0`, `sf_ru_0`), a `spectrum@dataProcessingRef` attribute only
C++ writes, a chromatogram time array whose `unitCvRef` is `UO` in Rust and `MS`
in C++, a `chromatogram/precursor/activation` subtree present only in C++, a
`cvParam MS:1000543` against a Rust-only `userParam openms-rust:empty-processing-actions`,
and a software term that differs (`MS:1000799` in Rust against `MS:1002135` /
`MS:1002146` in C++).

**Container parity.** Both sides write `indexedmzML` with an index. C++ writes
`<fileChecksum>0</fileChecksum>` (`CPP-049`) and an `indexListOffset` one byte
early (`CPP-305`); this port computes a real SHA-1 and addresses the opening
`<indexList` exactly. §5.1 quantifies what that costs, and labels the figure.

**Determinism and thread-invariance.** All 36 repetition checks and all 18
thread-invariance checks are `bitwise_equal` at the **data** level, on both
implementations: decoded arrays, features and `.dta` peak lines are identical
across every repetition and between 1 and 32 threads. That is *not* the same as
byte-identical files, and wave 3's report overstated it. Counting distinct
output sha256 per cell: DTAExtractor, MzMLSplitter and FileInfo produce one
distinct sha256 across all repetitions on both implementations — their files
really are byte-identical. BaselineFilter, MapNormalizer, PeakPickerHiRes,
FeatureFinderCentroided and SpectraFilterWindowMower produce 3 to 5 distinct
sha256 out of 3 to 5 repetitions, on **both** sides, because they stamp a
processing completion time into the output; on the Rust side the SHA-1
`fileChecksum` then follows from it. The correct statement is: the numeric
content is deterministic and thread-invariant everywhere; the files carry a
timestamp in five of the eight tools. The port's own bit-identity gate uses
`-test`, which suppresses the timestamp, and reproduces exactly.

### 3.7 Against wave 3 (2026-09-15, Rust main `fabd4b9`)

The control is strong: the C++ binaries are byte-identical across the two waves,
no Rust binary was reused, and all 18 C++ cells reproduce between waves within
± 0.9 % under essentially identical foreign load. Against that fixed control the
Rust side moved −19 % to −24 % at 1 thread on four tools, −62 % on
PeakPickerHiRes at 32 threads, and +24.5 % on DTAExtractor.

| tool | w3 t1 | w4 t1 | w3 t32 | w4 t32 | rust wall t1 | rust RSS t1 | what changed |
|---|---|---|---|---|---|---|---|
| DTAExtractor | 0.880 | **1.099** | 1.076 | 1.333 | 23.220 → 28.913 s | 1900 → 1691 MiB | wave 3's win was an artefact: the port wrote 658.9 MB against C++'s 836.5 MB at ~8 significant intensity digits instead of the source's two 15-digit rules. It now writes the source's formats, both sides emit 836,505,793 bytes and the output is byte-equal — and the formatting costs 24.5 % more wall. The one tool that got slower, and it got correct |
| MzMLSplitter | FAILED | **1.274** | FAILED | 1.891 | — | — | wave 3 refused every repetition ("precursor spectrum reference does not name an output spectrum"). Now runs the full 1.2 GB Velos file; all four parts bitwise equal |
| BaselineFilter | 2.053 | **1.579** | 2.270 | 1.752 | 32.587 → 24.964 s | 6675 → 6600 MiB | the rewritten mzML reader; C++ unchanged (15.87 → 15.81 s) |
| MapNormalizer | 1.400 | **1.130** | 2.055 | 1.671 | 20.636 → 16.671 s | 1987 → 1779 MiB | faster reader, and the normalisation defect is fixed. Wave 3's ratio compared two programs computing different answers; this one does not |
| SpectraFilterWindowMower | FAILED | **0.73–0.77** | FAILED | 0.80–0.88 | — | — | the 1,000,000-point cap was applied to the summed map (88.4 M peaks, 88× over) so the tool refused every real input; it is now per spectrum. First real-size comparison |
| PeakPickerHiRes | 1.320 | **1.010** | 1.397 | **0.531** | 33.190 → 25.367 s | 4310 → 3371 MiB | the largest single move. At 1 thread the reader rewrite and the spline scratch buffers close the gap entirely; at 32 the new parallel spectrum loop gives 2.01× against the C++ picker's 1.06×. Peak RSS drops 939 MiB, consistent with the ~886 MB claimed for that change |
| FileInfo (1.2 GB mzML) | 1.037 | **0.803** | 1.445 | 1.132 | 16.551 → 12.807 s | 2639 → 2428 MiB | reader rewrite; output byte-equal as before |
| FileInfo (60 MB featureXML) | FAILED | **1.486** | FAILED | 1.427 | — | — | wave 3: "identification XML byte limit exceeded" — the port could not read the featureXML its own FeatureFinderCentroided wrote. It now reads the 60 MB, 42,789-feature file and its report is byte-equal |
| FeatureFinderCentroided | 1.257 | 1.247 | 1.066 | 1.060 | 112.955 → 112.380 s | 313 → 310 MiB | unchanged, as expected — no commit in this window touches it. The wave-3 review found this cell bimodal (2 of 5 repetitions at ~121 s); **it did not recur**, all five repetitions fall in 112.26–112.39 s |

Three of the eight tools could not be compared at real size in wave 3 and now
can. Four 1-thread ratios improved by 19–24 %. Both of wave 3's data-correctness
findings are closed. One tool got slower because it stopped writing less than it
should.

## 4. The wave-5 run: the FMA build-flag question (dax, 2026-09-17)

Run directory
`/ceph/ibmi/abi/oliver/bench/openms4/results/2026-09-17-fma/fma-ffc-dax-1`
(`summary.md`), plan `fma-ffc`
(`/ceph/ibmi/abi/oliver/bench/openms4/fma-2026-09-17/config/plan_fma_ffc.json`,
sha256 `46be874d7147ed45`), started 2026-09-17T05:34:38+00:00. Harness git
`962a87e3de8b`. The staged notes are in the session scratchpad as
`fma-benchmark.md`; that file was written *before* the run and its section 4
still says "not measured" — the numbers below supersede it.

**This run is not comparable with §3.** It ran on **dax** (AMD EPYC 9654, 2
sockets, 96 cores each, 384 logical CPUs), not on **ibminode05** (EPYC 7763),
which is the wave-4 timing node and the only node §3's absolute figures may be
read against. Nothing here may be put in the same table as a wave-4 second, and
no wave-4 ratio may be updated from it. The five cells below are comparable
**with each other**, because they are interleaved repetitions of one plan on one
node in one session.

It is also a single tool on a single dataset: `FeatureFinderCentroided` on
`sub_centroid_velos_50amol_r1_first4000`, the same documented 4,000-spectrum
subset §3 uses, with the same shared INI (sha256 `2869134aeb3f98ed`, the C++
`-write_ini` defaults, no edits). 1 warm-up and **5 interleaved measured
rounds**, seed 20260917, at 1 and 32 threads. **0 of the repetitions were
load-flagged** (gate: flag above 0.25 per core, 0.05 at threads=1).

### 4.1 What ran

| cell | build |
|---|---|
| `cpp-release` | the §1 C++ Release build, core `bc9cc12` / cli `c19e494` / topp `174b576`, gcc 14.4 `-O3`, no `-march`, staged node-local |
| `rust-main` | `openms-rs` main `59e0e1c`, `cargo build --release --locked --offline --bins`, default features, **no `RUSTFLAGS`** (x86-64 baseline) |
| `rust-ffap` | `port/ffap-complete` `ddc35a7` (fix round 4), same build, no `RUSTFLAGS` |
| `rust-ffap-fma` | `port/ffap-complete` `ddc35a7`, same build with `RUSTFLAGS='-C target-feature=+fma'` |
| `rust-main-fma` | main `59e0e1c` with `+fma` — the **control**, which isolates the code-generation effect of the flag from the effect of inlining the port's `mul_add` sites |

The branch measured is `ddc35a7`, not the round-6 head the rest of this wave
records. Fix rounds 5 and 6 add a few branches per peak (the SSE NaN rules of
`MassTrace::avg_mz`, the intensity profile and `MassTraces::update_baseline`);
they are not expected to move these figures and were not measured.

### 4.2 Medians

Wall-clock median of 5 measured rounds, in seconds.

| cell | 1 thread | 32 threads | peak RSS 1 thread, MiB |
|---|---:|---:|---:|
| `cpp-release` | 82.069 | 22.641 | 360.0 |
| `rust-main` | 99.702 | 23.939 | 305.0 |
| `rust-ffap` | 121.051 | 21.946 | 319.0 |
| `rust-ffap-fma` | **86.581** | **20.233** | 321.0 |
| `rust-main-fma` | 97.373 | 23.579 | 306.0 |

### 4.3 Ratios

| pair | 1 thread | 95% CI | 32 threads | 95% CI |
|---|---:|---|---:|---|
| `rust-main` / C++ | 1.215 | [1.213, 1.218] | 1.057 | [1.054, 1.062] |
| `rust-ffap` / C++ | 1.475 | [1.468, 1.521] | 0.969 | [0.956, 0.974] |
| `rust-ffap-fma` / C++ | **1.055** | [1.053, 1.057] | **0.894** | [0.889, 0.907] |
| `rust-main-fma` / C++ | 1.186 | [1.185, 1.190] | 1.041 | [1.038, 1.044] |
| `rust-ffap` / `rust-main` | **1.214** | [1.208, 1.252] | 0.917 | [0.905, 0.920] |
| `rust-ffap-fma` / `rust-ffap` | **0.715** | [0.693, 0.719] | 0.922 | [0.917, 0.945] |
| `rust-main-fma` / `rust-main` | **0.977** | [0.974, 0.980] | 0.985 | [0.982, 0.987] |

Read plainly:

- **The completed FeatureFinderAlgorithmPicked costs 21 % at one thread on a
  default build.** `rust-ffap` / `rust-main` is 1.214. The cost is the port's
  bit-exact numerics: on a baseline x86-64 build every `f64::mul_add` of the
  ported glibc `powf`, `exp` and `log` becomes two indirect calls into
  `compiler_builtins`' `fma` stub (41 call sites in the tool binary), instead of
  one instruction.
- **`-C target-feature=+fma` removes that cost and more.** `rust-ffap-fma` /
  `rust-ffap` is 0.715 at one thread; against C++ the completed port goes from
  1.475 to **1.055** at one thread and from 0.969 to **0.894** at 32.
- **The control gains only 2.3 %.** `rust-main-fma` / `rust-main` is 0.977, so
  almost all of the 28.5 % is the `mul_add` inlining, not the flag's other
  effects.
- **Output is unchanged by the flag.** `rust-ffap-fma` against `rust-ffap` is
  `bitwise_equal` (data) and `equal` (metadata) at **both** thread counts, 4,076
  features matched, 0 unmatched, 0 charge disagreements. Every Rust cell is
  `equal_within_tolerance` against C++ on the same counts, and every cell's five
  repetitions are `bitwise_equal` to each other.

### 4.4 What the flag actually turns on

`-C target-feature=+fma` is **not** limited to the `mul_add` sites. rustc
implies `avx`, `sse3`, `ssse3`, `sse4.1` and `sse4.2` from `fma`
(`rustc --print cfg -C target-feature=+fma`), so the whole crate is compiled
with VEX encoding and 256-bit auto-vectorisation. In the measured binaries the
tool goes from 469 to 44,645 VEX instructions on the branch, and from 467 to
40,259 on main, which has no `mul_add` at all.

Two consequences, both for the decision and neither of them a measurement:

- **It raises the minimum CPU.** The binaries then require an FMA-capable
  processor — Intel Haswell or AMD Piledriver and later — instead of baseline
  x86-64.
- **A narrower alternative exists and was not measured.**
  `#[target_feature(enable = "fma")]` on the three ported replica functions with
  runtime dispatch would confine the change to them.

### 4.5 The decision is open

**Nothing in the build configuration was changed by this wave.** No
`RUSTFLAGS`, no `[profile.release]` override and no `.cargo/config.toml` entry
was added, on any branch or in CI; the `+fma` cells above were built out of tree
by the harness. Whether release builds should adopt `-C target-feature=+fma`,
use narrower per-function runtime dispatch, or accept the 21 % at one thread is
**open with the user**, and is recorded as open in `docs/VALIDATION.md` and
`docs/EARLY_TOPP_WORK_PACKAGES.md`.

### 4.6 Which caveats of §5 apply

§5's caveats on cross-session drift (treat differences below about 3 % as not
established without a second session), on the unfixed CPU governor, on the
differing build flags of the two sides, and on the 4,000-spectrum subset all
apply here unchanged. The three ratios called out above — 1.214, 0.715 and
0.977 — are 21 %, 28.5 % and 2.3 %; the first two are far outside that band,
and **the 2.3 % control figure is inside it**, so read that one as "small, and
not separated from noise by this session alone". §5's node-specific figures and
§6's list are written for the wave-4 run on ibminode05 and are not restated
here.

## 5. Caveats that apply to every number above

1. **The port does work the C++ side skips on every indexed mzML it writes, and
   the size of that work is a projection, not a measurement.**
   `src/format/mzml_write_options.rs:410`/`:496` computes a real SHA-1 over the
   output and emits it as `fileChecksum`; the C++ emits the literal `0`
   (`CPP-049`). The `sha1` crate is pinned with `features = ["force-soft"]`, so
   no SHA-NI. Throughput of that exact crate and feature was **measured** on
   this node at 793 MB/s (three trials over 2 GiB: 790.9, 793.0, 794.0). The
   per-tool shares below are that throughput **divided into the output size**:
   an arithmetic projection that assumes the cost is purely additive. *No
   no-hash build was ever timed*, so these are estimates, not a measured
   ablation, and they must not be read as "the port would otherwise be level".

   | tool | output bytes | projected SHA-1 s | 1-thread gap s |
   |---|---|---|---|
   | BaselineFilter | 3.336 GB | ~4.21 | 9.16 |
   | MzMLSplitter | 1.605 GB | ~2.02 | 3.98 |
   | MapNormalizer | 1.605 GB | ~2.02 | 1.92 |
   | PeakPickerHiRes | 0.536 GB | ~0.68 | 0.25 |
   | SpectraFilterWindowMower | 0.253 GB | ~0.32 | −190 (port ahead) |

   C++'s zero checksum is arguably the non-conforming side, so this is a
   disclosure, not a correction — but four of the five mzML ratios include work
   only one implementation performs.

2. **The two sides do not write the same number of bytes, and the reason is
   indentation, not metadata.** Rust/C++ output bytes: BaselineFilter 0.9958,
   MapNormalizer 0.9906, MzMLSplitter 0.9906, PeakPickerHiRes 0.9747,
   SpectraFilterWindowMower 0.9434; FeatureFinderCentroided 1.0072 the other
   way; DTAExtractor and both FileInfo cases exactly 1.0000. A byte census of
   the run's own retained outputs attributes **92–93 % of the deficit to XML
   indentation the port does not write**: PeakPickerHiRes deficit 13,912,722 B
   against a C++ tab surplus of 12,885,704 B (92.6 %); MapNormalizer 15,160,297
   / 13,990,168 (92.3 %); SpectraFilterWindowMower 15,161,155 / 13,990,180
   (92.3 %); BaselineFilter 13,912,527 / 12,885,644 (92.6 %). The Rust files
   contain **zero** tab characters; C++ also writes ~258,358 more newlines,
   about 6.3 extra lines per spectrum. The metadata differences of §3.6 account
   for tens of bytes per file, not megabytes: `dataProcessingRef=` occurs
   **exactly once per C++ file**, not once per spectrum — about 30 bytes.

   The consequence is an output-fidelity fact worth stating on its own: *every
   line of the port's mzML differs textually from the source's*, so the two
   files are not diffable even where all 81,714 or 87,492 arrays are bitwise
   identical.

3. **CPU frequency is not pinned, but both sides ran at full boost.** Governor
   `schedutil` on all 128 CPUs, boost on. The busiest core's `scaling_cur_freq`
   snapshot after each execution has a median of 3,525 MHz for Rust and
   3,523–3,525 MHz for C++, at both thread counts, against a 3,530 MHz maximum.
   The node-wide mean over 128 mostly-idle CPUs is 1,587–1,864 MHz and is *not*
   the clock the tool ran at — wave 3's caveat quoted that mean and called it
   the largest uncontrolled factor, which its review corrected. Caveat on the
   caveat: the harness records the maximum over all 128 CPUs, an upper bound on
   the core actually running the tool, so this does not by itself prove the C++
   32-thread cells were never clocked down.

4. **One session only. Treat ratio differences under about 3 % as not
   established.** Two independent sessions on this node have previously moved
   single-case medians by up to 2.1 % while within-session IQRs were 0.1–0.5 %.
   The bootstrap CIs resample one session's repetitions and describe
   within-session noise only. Exactly one row falls inside that band:
   **PeakPickerHiRes at 1 thread, 1.010** — "no difference measured".

5. **SpectraFilterWindowMower has n = 3 and an outlier at both thread counts.**
   Raw Rust walls: 1 thread [515.424, 515.513, 539.252] (one +4.6 %), 32 threads
   [516.175, 516.581, 565.100] (one +9.4 %). The C++ triples are tight
   (705.9 / 705.1 / 705.7 at 1 thread). The median lands on the low mode and the
   bootstrap CI is one-sided around it. Quote this tool as **0.73–0.77 at 1
   thread and 0.80–0.88 at 32**, never to three digits — which is why the tables
   above print ranges.

6. **FeatureFinderCentroided is measured on a 4,000-spectrum subset, not the
   full run** (`inputs/derived/sub_centroid_velos_50amol_r1_first4000.mzML`,
   91 MB, manifest `MANIFEST_ffc_subset_4000.json`), the same subset wave 3
   used. Neither implementation finishes the full 43,745-spectrum run. The
   evidence is asymmetric and is stated as such: the run's own pilot log records
   the **C++** side killed at a 600 s cap, and holds *no* wall or exit record for
   the Rust side (47 s elapsed, no `RUST WALL` line), so that file does not
   establish the Rust result. An **independent re-run of the Rust side at a
   700 s cap** was killed at the cap with no output. The conclusion holds; the
   pilot log's `RSS 1024` and `EXIT 0` must not be quoted, because the launcher
   there measured `timeout`, not the tool. A single full-size cell would cost
   over an hour and the four-round matrix over eight. Scaling behaviour on the
   full run may differ from the subset's.

7. **Start-up is under 1 % of every full-size case, and the correction is
   smaller than the `-write_ini` figure.** The C++ binaries pay 63–70 ms on
   `-write_ini` alone against the port's 3.0–4.1 ms, a difference of 62–67 ms.
   But the start-up-corrected column subtracts the **10-record slice**, whose
   C++ surplus is 8.4 ms (BaselineFilter), 12.9 (DTAExtractor), 13.0
   (MzMLSplitter), 18.5 (MapNormalizer), 28.5 (FileInfo/mzML), 46.3
   (FileInfo/featureXML) and 57 ms (FeatureFinderCentroided). Either way the
   correction moves only the ~1 s FileInfo/featureXML row (1.486 → 1.569 and
   1.427 → 1.520) — *against* the port, since the port's start-up is the cheaper
   one. Nothing favourable to the port is a start-up artefact.

8. **This is a warm-cache CPU benchmark.** `majflt` is 0 in every one of the 192
   measured repetitions: inputs were staged node-local and resident. Fair to
   both sides, excludes I/O, not a cold-start measurement.

9. **The node is not reserved.** Another user's two single-threaded Python jobs
   ran throughout (about 2 cores of 128). 0 of 192 repetitions were
   load-flagged. Exclusive use would need an admin or a Slurm reservation.

10. **ibminode06 was not used.** It answers SSH from ibminode05 but rejects the
    key, and it has no `ssh-config` entry on the workstation, so its load could
    not be read; nothing in this run touched it. (The C++ Release prefix was
    originally *built* there in an earlier session; that is recorded in its build
    manifest and is not part of this run.)

11. **Build configurations differ and cannot be equalised.** C++ is gcc 14.4
    `-O3 -DNDEBUG -mssse3 -ffp-contract=off`, shared libraries, no LTO; Rust is
    opt-level 3, thin-local LTO, codegen-units 16, statically linked, no
    `[profile.release]` override, so debug-assertions and overflow-checks are
    off. Neither uses native CPU flags — both are x86-64 baseline, the one thing
    that *could* be equalised. This measures two shipped configurations, not two
    compilers. The C++ side is a genuine Release build: `libOpenMS.so` and the
    tools carry zero `__assert_fail` references and zero "Assertion" strings.

12. **n < 10 on every large case** (n = 5, and n = 3 for
    SpectraFilterWindowMower). None of the large-case ratios is a precision
    measurement.

13. **Instruction counts are exact; wall clock on this workload is not.** The
    validation lane established, with a layout control, that one function
    (`Record::finish`) has two code-generation states 20,612,136 instructions
    apart on this workload and flips between them on source edits to files it
    never calls — the only source file between two `main` commits that flipped
    it is `src/cli/tools/map_normalizer.rs`, which is not on the
    PeakPickerHiRes path. Separately, a rebuild of an identical tree is worth a
    couple of tenths of a second either way. Any future comparison of two
    commits of this port must carry a layout control and must not attribute a
    sub-1 % wall difference to a source change. See
    [EARLY_TOPP_WORK_PACKAGES](EARLY_TOPP_WORK_PACKAGES.md).

## 6. What this run does not establish

- **Anything about FeatureFinderCentroided at production scale.** The
  comparison is a 4,000-spectrum subset; on the full run neither side finishes
  (§5.6).
- **Byte-level output compatibility with the source for any mzML-writing
  tool.** The port writes no indentation, so every line differs textually even
  where every array is bitwise identical (§5.2).
- **Memory parity.** The port uses 1.1×–2.7× the peak RSS on six of nine cases
  (§3.4).
- **Anything at better than about 3 % across sessions** (§5.4), and nothing at
  all about cold-start or I/O-bound behaviour (§5.8).
- **A decomposition of the mzML write path.** The SHA-1 share is projected from
  a throughput measurement, not ablated (§5.1).
- **Thread scaling beyond 1 and 32**, and nothing about the five Rust tools that
  ignore `-threads` beyond the fact that they ignore it (§3.5).
- **Any claim that the C++ 32-thread column had a doubled compute budget.** It
  did not; the surplus threads are an idle library pool the harness itself sizes
  (§3.5).

## 7. Reproducing

```
# C++ reference: already built and staged; verify against its manifest
sha256sum /ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576/BUILD_MANIFEST.json

# the wave-4 matrix
H=/ceph/ibmi/abi/oliver/bench/openms4/harness
python3 $H/bench.py run --plan $H/plans/w4threads.json --run-id <date>-w4threads
python3 $H/bench.py summarize --run-id <date>-w4threads

# the port's own bit-identity gate for PeakPickerHiRes (needs -test; without it
# the output carries a wall-clock completion time and a SHA-1 over it)
PeakPickerHiRes -in UK222.mzML -out o.mzML -ini pph.ini -threads 1 -test
# 535,613,726 B, sha256 bb13eecfe092a272b08ddc71feec3780c7bc876e8847e8a45b173cda9d2dad52
```

Results live under `/ceph/ibmi/abi/oliver/bench/openms4/results/`, indexed by
`results/RUNS.tsv`; the wave-4 run is `2026-09-16-w4threads` and its
`summary.md` carries the full per-case detail. The oracle artifacts cited here
are registered with their sha256 in `SOURCE_PROVENANCE.json`.
