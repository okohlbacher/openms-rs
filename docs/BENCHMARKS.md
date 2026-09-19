# Benchmarks: this port against the C++ OpenMS4 Release build

What is compared, how, and what the numbers do and do not say. The current
figures are the **wave-6 run of 2026-09-18**: all eight ported TOPP tools, on
full-size instrument data, at 1 and 32 threads, five (or three, or ten)
repetitions per cell on a quiet reserved node, with the equivalence of every
output judged separately for data and for metadata. They replace the wave-4
results, whose figures are kept beside them as the comparison column.

Two things changed under the tables since wave 4, and §3 is built to keep them
apart. FeatureFinderCentroided is a **different program** — the
FeatureFinderAlgorithmPicked port was completed in the 122 commits between the
runs. And every x86_64 binary is now built with **`-C target-feature=+fma`**,
which `.cargo/config.toml` sets for the whole checkout, so the default build
changed underneath all eight tools and not only the one the flag was introduced
for. A third implementation, the same commit built with the documented opt-out,
is carried through eight of the nine cases so that the flag's effect is measured
within this session rather than inferred across two (§3.7).

Section 6 lists what the run does **not** establish, and section 5 the caveats
that qualify every figure in it. Section 4 is the earlier, narrower wave-5 run
that framed the build-flag question on a different node; it is **not** comparable
with the wave-6 tables in absolute terms and is kept for how the decision was
reached.

Everything below is against **one** C++ reference: the optimised build at
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576`, whose
per-tool binaries are byte-identical to the ones wave 4 measured. The product
SDK used as the correctness oracle elsewhere in this repository is a *Debug*
build of a different core revision and is never used for timing.

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
unnecessary. Its CPU frequency is **not** pinned; the node was reserved for
the benchmark lane for the whole of the wave-6 run, and what that did and did
not buy is in caveats 9 and 14.

This run did not modify the harness. Its plans and `impls.json` live outside it,
under `/ceph/ibmi/abi/oliver/bench/openms4/w6-2026-09-18/config/`, and are passed
with `--plan` and `--impls`; `run.json` records the harness's git head, its
working-tree status and the sha256 of every harness file, and all three are
byte-identical to wave 4's.


## 3. The wave-6 run: eight tools, 1 and 32 threads, on the FMA default build

Run directory `/ceph/ibmi/abi/oliver/bench/openms4/results/2026-09-18-w6refresh`.
It replaces the wave-4 tables, which are kept beside the new ones in §3.2 and
§3.3 as the comparison column, and discussed tool by tool in §3.8.

Two things changed since wave 4, and the run is built to tell them apart where
it can. The **FeatureFinderAlgorithmPicked port was completed**, so
FeatureFinderCentroided is a different program from the one wave 4 measured. And
**every x86_64 binary is now built with `-C target-feature=+fma`**, because
`.cargo/config.toml` sets it and the harness's build script reaches it by
`unset`ting `RUSTFLAGS` and building from the checkout. Wave 4's binaries were
baseline x86-64. So the default build changed underneath every one of the eight
tools, not only the one the flag was introduced for.

To keep that from becoming an unattributable difference, the run carries a
third implementation on eight of the nine cases: `rust-nofma`, the same commit
built with the documented opt-out `RUSTFLAGS="-C target-feature=-fma"`. It is
interleaved with the other two in the same rounds on the same node, so
`rust-nofma / rust-release` is a within-session measurement of the flag and is
not exposed to the ~3 % cross-session drift band. §3.7 reports it.

### 3.1 What ran

| item | value |
|---|---|
| Rust, `rust-release` | `openms-rs` `36c26a0a301d8280faf5403c5e824f952d4a5138` (main). `cargo build --release --locked --offline --bins`, default features, **no `RUSTFLAGS`** — so `.cargo/config.toml` applies `-C target-feature=+fma` on x86_64. rustc 1.96.0, built on ibminode05 in 52 s. Source tarball sha256 `8fd36a41950798fa…`, build manifest sha256 `38d6be60012b2fe3` |
| Rust, `rust-nofma` | the same tarball, the same command and the same node, with `RUSTFLAGS="-C target-feature=-fma"`, which **replaces** the config flags. Build manifest sha256 `bfa49688cc6bbb02`. Present on every case except SpectraFilterWindowMower (§5, caveat C-B) |
| C++ | the §1 Release build. Its per-tool binary sha256 are **byte-identical to wave 4's**, staged node-local, `ldd_probe.resolved_under_ceph` = 0 |
| node | ibminode05 only, reserved for this lane for the whole run |
| harness | `/ceph/…/harness`, **unmodified**: git HEAD `962a87e3de8b` and `files_sha256` are byte-identical to wave 4's `run.json`, and the working-tree status string matches character for character. This run's plans and `impls.json` live outside the harness, under `/ceph/…/w6-2026-09-18/config/`, and are passed with `--plan` and `--impls` |
| inputs | the **same files wave 4 used**, verified by `dataset_sha256` per dataset, staged node-local to `/scratch/kohlbach/bench-w6/run/inputs` |
| parameters | one INI per tool, written by the C++ `-write_ini` and shared by all three implementations. `used_sha256` equals `cpp_write_ini_sha256` unedited on five of the eight tools (BaselineFilter, FeatureFinderCentroided, MapNormalizer, PeakPickerHiRes, SpectraFilterWindowMower). Three carry the harness's documented edits from `tools.json`, each applied once to the shared file and therefore symmetric across implementations: FileInfo `m`/`p`/`s` = true, MzMLSplitter `parts` = 4, DTAExtractor `level` = 2 |
| command | `-ini <shared> -in … -out … -threads N -no_progress`, with `RAYON_NUM_THREADS = OMP_NUM_THREADS = N` |
| order | implementations interleaved within each round in a seeded random order (seed 20260918) |
| repetitions | 1 discarded warm-up + 5 measured rounds; cells under 10 s raised to 10; SpectraFilterWindowMower 3 |
| temp | `TMPDIR` and `OPENMS_TMPDIR` under `/dev/shm`, per sub-run |
| peak-RSS floor | `/bin/true` through the launcher recorded 1,024 KiB before and after the cases (limit 4,096); the naive `Popen`+`wait4` path recorded the harness's own high-water mark instead |
| size | four sub-runs (`w6-profile`, `w6-centroid`, `w6-ffc-subset`, `w6-sfwm`), **52 cells, 282 measured repetitions of which 281 succeeded** (the one that did not is §3.9), 333 successful timing executions counting the 52 warm-ups, and 783 counting the `-write_ini` and start-up-slice baselines too |
| load | pre-cell gate value 0.0016–0.0125 per core against a flag limit of 0.25 (0.05 at threads = 1); foreign CPU during the measured repetitions 0.0012–0.0199 per core, median 0.0091 — both ranges and the median are the same whether the failed repetition is counted or not. **0 of the 282 measured repetitions load-flagged, 0 refused, 0 retried**; `majflt` 0 in every one. Exactly one execution of the whole run carried a flag, and it enters no table here: a 2.2 ms `-write_ini` start-up baseline (DTAExtractor, `rust-release`, threads = 1, repetition 3) saw 0.0565 per core against the single-thread limit of 0.05 |
| failures | **one**, on the C++ side: `FeatureFinderCentroided` at 32 threads died of `SIGSEGV` in one repetition (§3.9). Every Rust execution of both builds completed |
| wall | 2026-09-18T19:40 to 2026-09-19T01:05 local, sequential, nothing else of this project's on the node |

### 3.2 One thread

Wall is the median of the measured repetitions with the interquartile range.
Ratio is rust/cpp (**< 1 = Rust faster**) with a 10,000-resample percentile
bootstrap 95 % CI of the ratio of medians. The wave-4 column is that run's
ratio **recomputed from its own raw repetitions by the method used here** —
the same four measured sub-runs, the same filter (`phase == "timing"`, no
warm-up, `status == "ok"`) and the same median, with the pilot sub-runs of both
runs excluded on both sides — so the two columns really do come from one
computation rather than a transcription. Recomputed that way it reproduces the
published wave-4 table to three digits on **all eighteen cells**, with no
exception. "moved" is the change in the ratio, taken from the full-precision
medians rather than from the three-digit columns, so recomputing it by hand
from the printed ratios can differ in the last digit; the drift band is about
3 % (§5.2), so a figure smaller than that is **not** a result.

| tool | dataset | n (r,c) | rust wall med s [IQR] | cpp wall med s [IQR] | ratio rust/cpp | 95 % CI | wave-4 ratio | moved | rust RSS MiB | cpp RSS MiB | RSS ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DTAExtractor | Velos centroid 1.2 GB | 5,5 | 29.122 [29.119–29.129] | 26.345 [26.344–26.397] | 1.105 | [1.101, 1.113] | 1.099 | +0.6 % | 1691.0 | 1517.0 | 1.11 |
| MzMLSplitter | Velos centroid 1.2 GB | 5,5 | 18.741 [18.713–18.766] | 14.483 [14.448–14.485] | 1.294 | [1.292, 1.301] | 1.274 | +1.5 % | 4064.0 | 1524.0 | 2.67 |
| BaselineFilter | QE profile 2.3 GB | 5,5 | 25.187 [25.152–25.211] | 15.961 [15.857–15.974] | 1.578 | [1.565, 1.599] | 1.579 | −0.1 % | 6600.0 | 3200.9 | 2.06 |
| MapNormalizer | Velos centroid 1.2 GB | 5,5 | 16.752 [16.737–16.777] | 14.670 [14.659–14.690] | 1.142 | [1.138, 1.145] | 1.130 | +1.0 % | 1778.8 | 1518.0 | 1.17 |
| SpectraFilterWindowMower | Velos centroid 1.2 GB | 3,3 | 564.434 [547.575–575.153] | 699.092 [698.780–699.305] | **0.76–0.84** | [0.759, 0.839] | 0.730 | see §3.8 | 3345.8 | 1516.0 | 2.21 |
| PeakPickerHiRes | QE profile 2.3 GB | 5,5 | 25.407 [25.339–25.412] | 25.275 [25.218–25.319] | 1.005 | [1.000, 1.011] | 1.010 | −0.5 % | 3371.0 | 3880.5 | 0.87 |
| FileInfo | Velos centroid 1.2 GB | 5,5 | 12.716 [12.677–12.717] | 15.976 [15.970–15.982] | **0.796** | [0.793, 0.800] | 0.803 | −0.9 % | 2427.9 | 2008.0 | 1.21 |
| FileInfo | featureXML 60 MB | 10,10 | 1.508 [1.503–1.520] | 1.006 [1.003–1.007] | 1.499 | [1.492, 1.515] | 1.486 | +0.9 % | 151.0 | 79.0 | 1.91 |
| FeatureFinderCentroided | Velos 4,000-spectrum subset | 5,5 | 101.063 [100.948–101.223] | 89.974 [89.973–90.023] | **1.123** | [1.118, 1.126] | 1.247 | **−10.0 %** | 322.8 | 361.7 | 0.89 |

### 3.3 Thirty-two threads

| tool | dataset | n (r,c) | rust wall med s [IQR] | cpp wall med s [IQR] | ratio rust/cpp | 95 % CI | wave-4 ratio | moved | rust RSS MiB | cpp RSS MiB | RSS ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DTAExtractor | Velos centroid 1.2 GB | 5,5 | 29.303 [29.292–29.449] | 21.705 [21.633–21.725] | 1.350 | [1.342, 1.427] | 1.333 | +1.3 % | 1691.0 | 1523.7 | 1.11 |
| MzMLSplitter | Velos centroid 1.2 GB | 5,5 | 18.848 [18.834–18.848] | 9.729 [9.723–9.761] | 1.937 | [1.925, 1.941] | 1.891 | +2.5 % | 4063.8 | 1531.0 | 2.65 |
| BaselineFilter | QE profile 2.3 GB | 5,5 | 25.212 [25.073–25.247] | 14.406 [14.387–14.421] | 1.750 | [1.711, 1.762] | 1.752 | −0.1 % | 6600.0 | 3198.7 | 2.06 |
| MapNormalizer | Velos centroid 1.2 GB | 5,5 | 16.618 [16.564–16.733] | 9.929 [9.916–9.961] | 1.674 | [1.663, 1.699] | 1.671 | +0.2 % | 1778.9 | 1526.6 | 1.17 |
| SpectraFilterWindowMower | Velos centroid 1.2 GB | 3,3 | 537.475 [532.947–556.211] | 641.690 [641.058–641.737] | **0.82–0.90** | [0.823, 0.898] | 0.807 | see §3.8 | 3345.0 | 1532.7 | 2.18 |
| PeakPickerHiRes | QE profile 2.3 GB | 5,5 | 12.734 [12.703–12.764] | 23.810 [23.790–23.871] | **0.535** | [0.532, 0.537] | 0.531 | +0.6 % | 3370.0 | 3895.0 | 0.87 |
| FileInfo | Velos centroid 1.2 GB | 5,5 | 12.701 [12.685–12.797] | 11.336 [11.314–11.344] | 1.120 | [1.118, 1.133] | 1.132 | −1.0 % | 2427.9 | 2017.0 | 1.20 |
| FileInfo | featureXML 60 MB | 10,10 | 1.507 [1.497–1.516] | 1.054 [1.044–1.057] | 1.430 | [1.417, 1.447] | 1.427 | +0.2 % | 151.0 | 84.0 | 1.80 |
| FeatureFinderCentroided | Velos 4,000-spectrum subset | 5,**4** | 23.276 [23.247–23.474] | 25.216 [25.214–25.230] | **0.923** | [0.921, 0.933] | 1.060 | **−12.9 %** | 317.8 | 378.8 | 0.84 |

**The C++ FeatureFinderCentroided cell has n = 4, not 5.** One of its five
measured repetitions died of `SIGSEGV`; §3.9 is that event on its own, because
it is the only execution in either wave that failed.

### 3.4 Where the port is faster, where it is slower

**One tool moved, and it is the one whose algorithm was completed.**
FeatureFinderCentroided goes from 1.247 to **1.123** at one thread and from
1.060 to **0.923** at 32 — the first 32-thread case other than the picker and
the window mower where the port is outright faster than the C++ Release build.
§3.7 shows that most of that is the build flag, not the algorithm: without it
the same commit is 1.560 and 1.003.

**Six of the remaining seven are where wave 4 left them.** Their ratios moved by
between −0.9 % and +1.5 % at one thread and −1.0 % and +2.5 % at 32 — all inside
the ~3 % cross-session band, so **nothing else moved** is the honest reading, not
"MzMLSplitter got 2.5 % worse". The largest single figure is MzMLSplitter at 32
threads, +2.5 %, still inside the band; at one thread the widest is the same tool
at +1.5 %.

**The seventh, SpectraFilterWindowMower, is not resolvable.** Its medians moved
from 0.730 to 0.807 at one thread and 0.807 to 0.838 at 32, which looks like a
move outside the band — but at n = 3 with one repetition 8–10 % above the other
two, the bootstrap intervals are [0.759, 0.839] against wave 4's [0.730, 0.765]
and [0.823, 0.898] against [0.804, 0.883]. They overlap at 32 threads and touch
at one. The tables therefore print ranges for this tool, as wave 4's did, and
§5.10 gives the raw triples. What both waves do agree on is the direction and
its size: every Rust repetition finishes well ahead of every C++ one.

**Faster than C++.** FileInfo on the 1.2 GB mzML at one thread (0.796),
PeakPickerHiRes at 32 threads (0.535), FeatureFinderCentroided at 32 threads
(0.923), and SpectraFilterWindowMower at both (§3.2/§3.3).

**Level.** PeakPickerHiRes at one thread, 1.005, CI [1.000, 1.011] — inside the
band, so *no difference measured*, exactly as in wave 4.

**Slower.** DTAExtractor 1.105 / 1.350, MapNormalizer 1.142 / 1.674,
MzMLSplitter 1.294 / 1.937, BaselineFilter 1.578 / 1.750, FileInfo on the 60 MB
featureXML 1.499 / 1.430, FeatureFinderCentroided at one thread 1.123, and
FileInfo on the mzML at 32 threads 1.120.

**Memory.** The port uses more resident memory on **seven of the nine cases**
— six of the eight tools, FileInfo being measured on two datasets — worst
MzMLSplitter 4.06 GB against 1.52 GB (2.67×), BaselineFilter 6.60 GB against
3.20 GB (2.06×); the two below parity are PeakPickerHiRes (0.87×) and
FeatureFinderCentroided (0.89× at one thread, 0.84× at 32). There is no memory
parity and none is claimed.

**Memory against wave 4: one tool moved, the same one.** The peak RSS of both
sides reproduces wave 4 on every tool **except FeatureFinderCentroided**, whose
algorithm was completed in the window (§5.4). Its Rust peak rose from 310.3 to
322.8 MiB at one thread (+4.0 %) and from 306.4 to 317.8 MiB at 32 (+3.7 %)
against a C++ side that did not move (360.1 → 361.7 and 378.9 → 378.8 MiB), so
its RSS ratio went from 0.86 to **0.89** at one thread and from 0.81 to **0.84**
at 32. On every other tool the wave-6 Rust median is within 0.09 % of wave 4's
and the C++ median within 0.13 %; twelve of those sixteen printed RSS ratios are
identical across the waves and the other four differ by one unit in the last
digit, from that sub-tenth-of-a-percent movement alone — DTAExtractor
1.12 → 1.11 at one thread, and at 32 threads MapNormalizer 1.16 → 1.17,
SpectraFilterWindowMower
2.19 → 2.18 and PeakPickerHiRes 0.86 → 0.87. The wave-4 figures quoted
here are
medians of the same four measured sub-runs under the same filter as §3.2's
wave-4 ratio column, so they are comparable cell for cell.

### 3.5 Thread behaviour, sampled rather than assumed

`peak thr` is the maximum thread count in `/proc/<tool>/task`, sampled every
5 ms. `CPU util` is (user + sys)/wall, the average number of CPUs the process
actually kept busy. A tool that asked for 32 threads and shows util ~1.0 did
**not** use them, however many threads exist.

| tool | dataset | impl | t1 peak thr | t1 util | t32 peak thr | t32 util | t32 user s | t32/t1 wall |
|---|---|---|---|---|---|---|---|---|
| DTAExtractor | Velos 1.2 GB | rust `+fma` | 2 | 1.00 | 33 | 1.00 | 26.3 | 0.99× |
| DTAExtractor | Velos 1.2 GB | rust `-fma` | 2 | 1.00 | 33 | 1.00 | 25.9 | 0.99× |
| DTAExtractor | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | 6.16 | 130.4 | 1.21× |
| MzMLSplitter | Velos 1.2 GB | rust `+fma` | 2 | 0.99 | 33 | 0.99 | 14.1 | 0.99× |
| MzMLSplitter | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | 12.47 | 118.7 | 1.49× |
| BaselineFilter | QE 2.3 GB | rust `+fma` | 2 | 1.00 | 33 | 1.00 | 17.6 | 1.00× |
| BaselineFilter | QE 2.3 GB | cpp | 2 | 1.00 | 64 | 9.03 | 123.3 | 1.11× |
| MapNormalizer | Velos 1.2 GB | rust `+fma` | 2 | 1.00 | 33 | 1.00 | 13.6 | 1.01× |
| MapNormalizer | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | 12.29 | 119.1 | 1.48× |
| SpectraFilterWindowMower | Velos 1.2 GB | rust `+fma` | 2 | 1.00 | 33 | 1.00 | 534.4 | 1.05× |
| SpectraFilterWindowMower | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | **1.18** | 752.0 | 1.09× |
| PeakPickerHiRes | QE 2.3 GB | rust `+fma` | 1 | 1.00 | 33 | 2.23 | 24.7 | 2.00× |
| PeakPickerHiRes | QE 2.3 GB | cpp | 2 | 1.00 | 64 | 5.85 | 133.9 | 1.06× |
| FileInfo | Velos 1.2 GB | rust `+fma` | 1 | 1.00 | 1 | 1.00 | 11.1 | 1.00× |
| FileInfo | Velos 1.2 GB | cpp | 2 | 1.00 | 64 | 10.88 | 121.5 | 1.41× |
| FileInfo | featureXML 60 MB | rust `+fma` | 1 | 1.00 | 1 | 1.00 | 1.4 | 1.00× |
| FileInfo | featureXML 60 MB | cpp | 2 | 1.00 | 33 | 4.23 | 4.4 | 0.95× |
| FeatureFinderCentroided | Velos 4,000 subset | rust `+fma` | 1 | 1.00 | 33 | **4.82** | 111.9 | **4.34×** |
| FeatureFinderCentroided | Velos 4,000 subset | rust `-fma` | 1 | 1.00 | 33 | 6.03 | 152.2 | 5.55× |
| FeatureFinderCentroided | Velos 4,000 subset | cpp | 2 | 1.00 | 64 | 4.22 | 106.2 | 3.57× |

Of the six tools measured without fused multiply-adds, **DTAExtractor's `-fma`
row is kept as the exemplar and the other five are left out**, because on all
six the `-fma` build is indistinguishable from the `+fma` one and the rows would
only repeat each other: peak threads identical in every cell, utilisation
identical to two decimals in every cell but PeakPickerHiRes at 32 threads (2.23
against 2.22), and user time within 0.45 s in all twelve cells — widest
DTAExtractor, 0.440 s at one thread and 0.435 s at 32, the latter being the
difference between the 26.3 s and 25.9 s its two rows print. That is itself the
point.
SpectraFilterWindowMower has no `-fma` row because it has no `-fma` arm
at all (§5.3). On FeatureFinderCentroided the two builds are **not**
indistinguishable, which is why both of its rows are shown: the `-fma`
build reaches util 6.03 and scales 5.55× where the `+fma` build reaches
4.82 and 4.34×. **That is
not better parallelism, it is more work to spread** — 152.2 s of thread CPU
against 111.9 s for the same 4,076 features. A speed-up ratio flatters the
slower build.

Otherwise the wave-4 picture is unchanged and reproduces closely: the port uses
its pool in two tools of eight (PeakPickerHiRes 2.00× and
FeatureFinderCentroided 4.34×), five build the 33-thread pool and leave it idle
(util 1.00, user time unchanged from one thread), and FileInfo builds no pool at
all. The C++ side peaks at 64 threads in every 32-thread cell but FileInfo on
the featureXML (33), and at 2 at `-threads 1`; as wave 4 established, that
surplus is a library pool the harness itself sizes through `OMP_NUM_THREADS`
and not a doubled compute budget — the per-cell utilisations run from **1.18 to
12.47**, not 64. The floor of that range is the clearest case: C++
SpectraFilterWindowMower holds 64 threads and spends 752.0 s of user CPU to turn
699.1 s of wall into 641.7 s, a utilisation of 1.18 for a 1.09× gain.

### 3.6 Output equivalence, data and metadata judged separately

Compared under **lead decision D6** — XML outputs are compared by decoded
content (FuzzyDiff's per-number rule, id exclusion, exact structure), never by
line layout — declared in
[EARLY_TOPP_WORK_PACKAGES](EARLY_TOPP_WORK_PACKAGES.md) and implemented by the
comparison policy of [DIFFERENTIAL_VALIDATION](DIFFERENTIAL_VALIDATION.md).
Neither comparison stops at the first difference. Tolerances are per tool
(`tools.json`): mzML default 0.001 ppm m/z and 1e-6 relative intensity;
PeakPickerHiRes 0.01 ppm / 1e-5; FeatureFinderCentroided 0.01 ppm / 0.01 s /
1e-4 relative with a 10 ppm, 5 s match window.

| tool | dataset | data | metadata | what the comparison counted |
|---|---|---|---|---|
| DTAExtractor | Velos 1.2 GB | bitwise_equal | n/a | 36,443 files, **36,443 byte-equal** |
| MzMLSplitter | Velos 1.2 GB | equal_within_tolerance | **DIFFERENT** | 87,492 arrays over four parts, **all bitwise**, 0 different; 57 metadata differences in 41 categories summed over the four parts |
| BaselineFilter | QE 2.3 GB | equal_within_tolerance | **DIFFERENT** | 81,713 of 81,714 arrays bitwise, 1 within tolerance, **0 different**; 14 differences in 12 categories |
| MapNormalizer | Velos 1.2 GB | equal_within_tolerance | **DIFFERENT** | 87,492 of 87,492 arrays bitwise, **0 different**; 20 differences in 14 categories |
| SpectraFilterWindowMower | Velos 1.2 GB | equal_within_tolerance | **DIFFERENT** | 87,492 of 87,492 arrays bitwise, **0 different**; 22 differences in 16 categories |
| PeakPickerHiRes | QE 2.3 GB | equal_within_tolerance | **DIFFERENT** | 81,714 of 81,714 arrays bitwise, **0 different**; 16 differences in 14 categories |
| FileInfo | Velos 1.2 GB | bitwise_equal | n/a | 1 file, byte-equal |
| FileInfo | featureXML 60 MB | bitwise_equal | n/a | 1 file, byte-equal |
| FeatureFinderCentroided | Velos 4,000 subset | equal_within_tolerance | **equal** | 4,076 of 4,076 features matched, **0 unmatched either way**, 0 charge disagreements, m/z 0.0 ppm, RT 0.0 s, largest relative intensity difference 9.477e-8, largest absolute quality difference 5.0e-7, same order |

The verdicts are **identical at 1 and at 32 threads**, and every count above is
identical to wave 4's. No tool's agreement with the C++ build changed, in either
direction, including the one tool whose algorithm was rewritten: the completed
FeatureFinderAlgorithmPicked still matches all 4,076 features with zero m/z and
zero RT difference, and its metadata is still `equal`.

**Metadata still differs on the five mzML writers**, in the same narrow class
wave 4 documented: identifier spelling and processing provenance, never
numbers — the `dataProcessing` / `software` / `sourceFile` id strings, a
`spectrum@dataProcessingRef` attribute only C++ writes, a chromatogram time
array whose `unitCvRef` is `UO` in Rust and `MS` in C++, a
`chromatogram/precursor/activation` subtree only C++ writes, a `cvParam
MS:1000543` against a Rust-only `userParam
openms-rust:empty-processing-actions`, and a differing software term.

**Container parity** is unchanged: C++ writes `<fileChecksum>0</fileChecksum>`
(`CPP-049`) and an `indexListOffset` one byte early (`CPP-305`); this port
computes a real SHA-1 and addresses the opening `<indexList` exactly. §5.1
quantifies what that costs and labels the figure a projection.

**Determinism and thread-invariance.** All **52** repetition checks (18
`rust-release`, 18 `cpp-release`, 16 `rust-nofma`) and all **26**
thread-invariance checks (9, 9 and 8) are `bitwise_equal` at the **data** level,
on all three implementations. Nine cases, not eight, because FileInfo is
measured on two datasets; `rust-nofma` has one case fewer because
SpectraFilterWindowMower has no `-fma` arm (§5.3). As in
wave 4 that is not the same as byte-identical files: counting distinct output
sha256 per cell, DTAExtractor, MzMLSplitter and FileInfo produce one distinct
sha256 across all repetitions, while the five tools that stamp a processing
completion time into the output produce as many distinct sha256 as they have
repetitions — on the C++ side too, and on **two** of its cells only four of
five, because two repetitions happened to land in the same second:
MapNormalizer at one thread and PeakPickerHiRes at 32. The C++ 32-thread
FeatureFinderCentroided cell also carries four distinct sha256, but out of the
four repetitions that wrote an output rather than five, for the reason §3.9
gives; that is as many distinct files as it has repetitions, not a collision.
The numeric content is deterministic and thread-invariant everywhere; five of
the eight tools carry a timestamp.

### 3.7 What `-C target-feature=+fma` is worth, measured

#### The binaries, before any timing

Both builds come from the same `git archive` of `36c26a0`, the same command and
the same toolchain on the same node; `rust-nofma` only exports
`RUSTFLAGS="-C target-feature=-fma"`, which **replaces** the
`.cargo/config.toml` flags rather than adding to them.

| tool | `vfma*` `+fma` | `vfma*` `-fma` | `%ymm` `+fma` | `%ymm` `-fma` |
|---|---:|---:|---:|---:|
| BaselineFilter | 0 | 0 | 10,209 | 162 |
| DTAExtractor | 0 | 0 | 8,663 | 162 |
| **FeatureFinderCentroided** | **41** | 2 | 16,147 | 162 |
| FileInfo | 0 | 0 | 10,463 | 162 |
| MapNormalizer | 0 | 0 | 10,095 | 162 |
| MzMLSplitter | 0 | 0 | 10,416 | 162 |
| PeakPickerHiRes | 0 | 0 | 10,766 | 162 |
| SpectraFilterWindowMower | 0 | 0 | 10,332 | 162 |

**Only FeatureFinderCentroided contains a fused multiply-add**, and exactly
**41** of them — the same 41 call sites `docs/FMA_BUILD_FLAG.md` counts in the
ported glibc `exp`, `log` and `powf`. The other seven have **zero** in *both*
builds, so for them the flag cannot remove an out-of-line `mul_add` call,
because there was none. (The 2 left in the `-fma` FeatureFinderCentroided are
not on the ported math path.)

**That does not make the flag a no-op for them.** rustc implies `avx`, `sse3`,
`ssse3`, `sse4.1` and `sse4.2` from `fma`, so the whole crate is recompiled with
VEX encoding and 256-bit auto-vectorisation: every tool goes from a flat 162
`%ymm` operands to between 8,663 and 16,147. Whatever the flag does to the seven
non-FFC tools it does through auto-vectorisation. Reading "0 `vfma*`, therefore
no effect" off that table would be an inference, not a measurement — which is
why the run measures it.

#### The measurement

`rust-nofma / rust-release` at the same thread count, interleaved in the same
rounds on the same node. **Above 1 means the flag is faster.** Both builds are
the same commit, so this ratio is not exposed to cross-session drift.

| tool | dataset | threads | `-fma` med s | `+fma` med s | `-fma` / `+fma` | 95 % CI | verdict |
|---|---|---|---|---|---|---|---|
| DTAExtractor | Velos 1.2 GB | 1 | 28.681 | 29.122 | 0.985 | [0.979, 0.993] | inside the band |
| DTAExtractor | Velos 1.2 GB | 32 | 28.901 | 29.303 | 0.986 | [0.933, 0.999] | inside the band |
| MzMLSplitter | Velos 1.2 GB | 1 | 18.599 | 18.741 | 0.992 | [0.987, 0.997] | inside the band |
| MzMLSplitter | Velos 1.2 GB | 32 | 18.809 | 18.848 | 0.998 | [0.989, 1.003] | inside the band |
| BaselineFilter | QE 2.3 GB | 1 | 25.563 | 25.187 | 1.015 | [0.987, 1.023] | inside the band |
| BaselineFilter | QE 2.3 GB | 32 | 25.331 | 25.212 | 1.005 | [1.001, 1.027] | inside the band |
| MapNormalizer | Velos 1.2 GB | 1 | 16.758 | 16.752 | 1.000 | [0.993, 1.003] | inside the band |
| MapNormalizer | Velos 1.2 GB | 32 | 16.743 | 16.618 | 1.008 | [0.990, 1.011] | inside the band |
| PeakPickerHiRes | QE 2.3 GB | 1 | 25.204 | 25.407 | 0.992 | [0.990, 0.998] | inside the band |
| PeakPickerHiRes | QE 2.3 GB | 32 | 12.697 | 12.734 | 0.997 | [0.994, 1.007] | inside the band |
| FileInfo | Velos 1.2 GB | 1 | 12.804 | 12.716 | 1.007 | [1.001, 1.011] | inside the band |
| FileInfo | Velos 1.2 GB | 32 | 12.833 | 12.701 | 1.010 | [0.997, 1.012] | inside the band |
| FileInfo | featureXML 60 MB | 1 | 1.497 | 1.508 | 0.993 | [0.981, 0.997] | inside the band |
| FileInfo | featureXML 60 MB | 32 | 1.496 | 1.507 | 0.993 | [0.984, 1.001] | inside the band |
| **FeatureFinderCentroided** | Velos 4,000 subset | 1 | 140.325 | 101.063 | **1.388** | [1.355, 1.409] | **established** |
| **FeatureFinderCentroided** | Velos 4,000 subset | 32 | 25.303 | 23.276 | **1.087** | [1.075, 1.091] | **established** |
| SpectraFilterWindowMower | Velos 1.2 GB | 1, 32 | — | — | **not measured** | — | §5, C-B |

**Where it is measurable it is large.** On FeatureFinderCentroided the flag is
worth **38.8 %** at one thread and **8.7 %** at 32. Against the C++ build that
is the difference between 1.560 and **1.123** at one thread, and between 1.003
and **0.923** at 32: the flag is what puts this tool ahead of the C++ Release
build at 32 threads, and it removes four fifths of its one-thread deficit.

**Where it is not measurable, say so plainly.** On the other seven cases all
fourteen ratios lie between 0.985 and 1.015 — inside the ~3 % band on every one
— and the sign is **not even consistent**: the flag is nominally faster on
BaselineFilter, MapNormalizer and FileInfo-on-mzML and nominally slower on
DTAExtractor, MzMLSplitter, PeakPickerHiRes and FileInfo-on-featureXML. The
correct statement is **no effect was measured on any tool without a fused
multiply-add**, not "the flag costs DTAExtractor 1.5 %". The AVX vectorisation
the flag also switches on is real in the binaries and did not show up in the
wall clock of these seven workloads, which are dominated by XML parsing,
base64 and I/O rather than by float arithmetic.

**The flag does not change any output.** All sixteen `-fma` against `+fma`
comparisons are `bitwise_equal` on **data**, at both thread counts, on every
case measured. On **metadata** ten of them are `equal` and the other six are
`not_applicable`: DTAExtractor and the two FileInfo cases write DTA or plain
text, which carries no metadata to compare — the same three cases §3.6 prints
as `n/a` against the C++ build. On FeatureFinderCentroided all four maximum
differences — m/z, RT, intensity and quality — are **exactly 0.0**, over all
4,076 features;
the two files differ by 8 bytes of timestamp. This reproduces the wave-5 finding
on a second node with a completed algorithm.

### 3.8 Against wave 4, tool by tool

The control is strong: the C++ binaries are byte-identical across the two waves
(`binary_sha256` per tool), the inputs are byte-identical (`dataset_sha256` per
dataset), the harness is byte-identical (`files_sha256` and even the
working-tree status string), and the node is the same. Against that fixed
control **the C++ medians reproduce between the two waves within ±0.96 % on all
eighteen cells** (widest: BaselineFilter at one thread, +0.96 %, and
SpectraFilterWindowMower at one thread, −0.94 %),
which is the run's own estimate of what a session boundary is worth and the
reason the drift band of §5.2 is set where it is.

| tool | w4 t1 | w6 t1 | w4 t32 | w6 t32 | what moved, and why |
|---|---|---|---|---|---|
| DTAExtractor | 1.099 | 1.105 | 1.333 | 1.350 | nothing established. `src/cli/tools/dta_extractor.rs` is untouched in the window, but 32 source files changed in it and the shared mzML reader, the progress logger, `cli.rs`, `kernel.rs` and `metadata/value.rs` are all on this tool's path, so this row is *not* an unchanged-code control; the flag gives it no `vfma*` |
| MzMLSplitter | 1.274 | 1.294 | 1.891 | 1.937 | nothing established. +2.5 % at 32 threads is the run's largest move and is still inside the band |
| BaselineFilter | 1.579 | 1.578 | 1.752 | 1.750 | nothing, to three digits, at both thread counts |
| MapNormalizer | 1.130 | 1.142 | 1.671 | 1.674 | nothing established |
| SpectraFilterWindowMower | 0.730 | 0.76–0.84 | 0.807 | 0.82–0.90 | **not resolvable at n = 3.** The port is substantially faster in both waves — every Rust repetition beats every C++ repetition — but the wave-6 Rust repetitions are more scattered than wave 4's, and the bootstrap intervals overlap at 32 threads and touch at one. This tool has no `-fma` arm (§5.3), so nothing here separates the build flag from the scatter |
| PeakPickerHiRes | 1.010 | 1.005 | 0.531 | 0.535 | nothing established, at either thread count, despite `peak_picking/noise.rs` and the new `noise_estimation.rs` changing in the window. Still level at one thread and still 1.87× faster at 32 |
| FileInfo (1.2 GB mzML) | 0.803 | 0.796 | 1.132 | 1.120 | nothing established; the one-thread win is intact |
| FileInfo (60 MB featureXML) | 1.486 | 1.499 | 1.427 | 1.430 | nothing established |
| **FeatureFinderCentroided** | 1.247 | **1.123** | 1.060 | **0.923** | the only tool that moved. Two causes, and they are separable: the FeatureFinderAlgorithmPicked port was completed (`instance.rs`, `source_sort.rs`, `scoring.rs`, `seeds.rs`, and the bit-exact `glibc_libm.rs` / `glibc_powf.rs`), which makes the tool *slower* — the `-fma` build of this commit takes 140.3 s where wave 4's incomplete algorithm took 112.4 s — and the build flag, which more than pays that back (§3.7) |

So the summary of the movement is short: **one of eight tools moved, and it is
the one that was rewritten and the one the flag was introduced for.** Six of the
other seven reproduce wave 4 inside the drift band, which is also the strongest
statement this run makes about its own repeatability; the seventh,
SpectraFilterWindowMower, is too noisy at n = 3 to say either way, and this run
did not buy it a `-fma` arm that might have narrowed the question.

### 3.9 The one execution that failed

Across wave 4's **228** and wave 6's **334** timing executions — the four
measured sub-runs of each run, with the `case_load` bookkeeping rows (36 and 52)
excluded and the warm-ups included; 244 and 337 if the pilot sub-runs of both
runs are counted as well — exactly **one** execution failed, and it is on the
C++ side:

| field | value |
|---|---|
| tool | `FeatureFinderCentroided`, C++ Release (`binary_sha256` `2781dd7cab48f482…`) |
| case | Velos 4,000-spectrum subset, `-threads 32`, `OMP_NUM_THREADS=32`, measured repetition 2 of 5 |
| outcome | `SIGSEGV` — signal 11, launcher exit 139, `Command terminated by signal 11` |
| when | 7.094 s wall, 12.818 s user, 64 live threads, peak RSS 340,628 KiB |
| output | none — 0 files written |
| load | foreign CPU 0.0092 per core, **not** load-flagged; the pre-cell gate value was 0.0098 per core, the same value every repetition of that cell ran under, including the four that succeeded |
| parameters | the same INI sha256 `2869134aeb3f98ed…` as the four repetitions that succeeded |
| last output | stdout `Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.`; stderr the known non-fatal `DateTime conversion error of "-infinity"` |

It is **intermittent and thread-count-dependent as far as this evidence goes**:
one crash in six executions at 32 threads in this run, none in six at 32 threads
in wave 4, and none in any one-thread execution of either wave — 1 in 12
observed at 32 threads, 0 in 12 at 1 thread. Neither Rust build crashed in any
of its 24 executions of this cell.

Two consequences, and no more than two. The C++ 32-thread median in §3.3 rests
on **n = 4**, which is why that cell is marked. And this is a defect
observation against the pinned C++ Release build that belongs in
`OpenMS_CPP_ISSUES.md`; the benchmark lane does not own that file, so the entry
is proposed to the integrator rather than written here. **One event is not a
diagnosis**: nothing here identifies the faulting code, and no attempt was made
to reproduce it under a debugger, which would need a separate run.

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

### 4.5 The decision: on by default since `port/fma-default`

**The user decided on 2026-09-18 to enable the flag by default.**
`.cargo/config.toml` now sets `-C target-feature=+fma` for
`cfg(target_arch = "x86_64")`, with a startup check that refuses a processor
without FMA instead of letting it die on `SIGILL`. `docs/FMA_BUILD_FLAG.md`
records the scope, the opt-out and everything that was measured. Until that
change nothing in the build configuration was set on any branch or in CI, and
the `+fma` cells above were built out of tree by the harness.

**For the harness this changes what "no `RUSTFLAGS`" means.**
`build/build_rust_orig.sh` extracts a tarball of the repository, `unset`s
`RUSTFLAGS` and runs `cargo build --release --locked --offline --bins` from the
extracted checkout, so from that commit on it picks the flag up from
`.cargo/config.toml` without being told. That the archive carries the file is
checked here: `git archive HEAD | tar -t` lists `.cargo/config.toml`, and
`.gitattributes` marks nothing under `.cargo/` `export-ignore`. The script
itself lives on the cluster and was not read at this integration; the FMA lane
reports that it warns about the `--offline` consequence below. Three
consequences. **`--offline` will fail until the shared `~/.cargo` cache on the
benchmark nodes holds `raw-cpuid 11.6.0`**, the one crate this change adds. A `rust-*` cell built at or after `port/fma-default` is an FMA build and is **not** comparable
with the `rust-main` and `rust-ffap` cells of this section, which were baseline
builds; the comparable cells are `rust-ffap-fma` and `rust-main-fma`. And a cell
that wants a baseline build must now say `RUSTFLAGS="-C target-feature=-fma"`,
where before it got one by saying nothing.

**Since 2026-09-18 the flag is not a proposal but the build.** Every x86_64
binary produced from this checkout carries it, including the eight measured in
§3, and `src/system/cpu_features.rs` refuses a processor that cannot run them.
What that costs and buys is no longer a single tool's figure on another node:
§3.7 measures the flag against its own opt-out on eight of the nine cases of the
wave-4 matrix, on the timing node, within one session. Read §3.7 for the
current numbers and this section for how the question was framed when it was
still open — the two are on different nodes and must not be put in one table.

### 4.6 Which caveats of §5 apply

§5's caveats on cross-session drift (treat differences below about 3 % as not
established without a second session), on the unfixed CPU governor, on the
differing build flags of the two sides, and on the 4,000-spectrum subset all
apply here unchanged. The three ratios called out above — 1.214, 0.715 and
0.977 — are 21 %, 28.5 % and 2.3 %; the first two are far outside that band,
and **the 2.3 % control figure is inside it**, so read that one as "small, and
not separated from noise by this session alone". §5's node-specific figures and
§6's list are written for the **wave-6** run on ibminode05 and are not restated
here; this section's own figures are dax's and stand on their own.


## 5. Caveats that apply to every number above

1. **The wave-4 column and the wave-6 column do not differ by one thing.**
   Between `9a392fe` and `36c26a0` there are 122 commits, and the default
   x86_64 build acquired `-C target-feature=+fma`. A wave-4-to-wave-6
   difference therefore mixes three causes: the source changes, the flag, and
   cross-session drift. Only one is separated here, and only within this
   session: `rust-nofma` is the same commit, the same command and the same node
   as `rust-release`, so `rust-nofma / rust-release` isolates the flag exactly
   (§3.7). Nothing isolates the source changes from drift. Where §3.8 says a
   tool "did not move", it means the two ratios agree inside the drift band —
   not that the code is unchanged.

2. **One session only. Treat ratio differences under about 3 % as not
   established.** Two independent sessions on this node have previously moved
   single-case medians by up to 2.1 % while within-session IQRs were 0.1–0.5 %.
   The bootstrap CIs resample one session's repetitions and describe
   within-session noise only. In §3.2 and §3.3 **every wave-4-to-wave-6 move
   that those tables print, except FeatureFinderCentroided's, is inside this
   band**, and so are all fourteen non-FFC entries of the §3.7 flag table. None
   of them is a result. The figures that *are* outside the band:
   FeatureFinderCentroided's −10.0 % and −12.9 % against wave 4, and its flag
   ratios 1.388 and 1.087. One move is outside the band but is not printed as a
   move: SpectraFilterWindowMower's point ratio goes from 0.730 to 0.807 at one
   thread (+10.5 %) and from 0.807 to 0.838 at 32 (+3.8 %), and at n = 3 its
   bootstrap intervals overlap wave 4's at 32 threads and touch at one, so §3.4
   and caveat 10 below treat it as unresolved and the tables print a range and
   "see §3.8" instead of a figure. Unresolved is not the same as unmoved.

3. **The flag's effect on SpectraFilterWindowMower was not measured.** That plan
   carries no `-fma` arm: a third implementation would have added about 69
   minutes of wall for one more cell, and the cell already carries the n = 3
   outlier caveat below, which would dominate any flag effect. §3.7's static
   table shows the tool has zero `vfma*` instructions and the same AVX change as
   the other six, so it is *expected* to behave like them — expectation, not
   measurement.

4. **FeatureFinderCentroided is not the same program wave 4 measured.** The
   FeatureFinderAlgorithmPicked port was completed in this window. Its wave-4
   ratio measured an earlier, incomplete algorithm. The two numbers share a row
   because they are the same tool on the same input, not because they are two
   measurements of one program.

5. **The processor-refusal path was not executed.** Every binary here ran on an
   EPYC 7763, which has FMA3, so the startup check took its passing branch in
   every execution and its refusing branch in none. The refusal is covered by
   unit tests over the decision function `message_for(build_requires_fma,
   cpu_provides_fma)` — including the refused combination and the exact rebuild
   command in the message — but no pre-2013 processor ran a binary from this
   build, here or anywhere in this repository's records.

6. **One C++ execution crashed, and one event is not a diagnosis.** §3.9 has it
   in full. It costs the C++ 32-thread FeatureFinderCentroided cell one
   repetition (n = 4). Nothing here identifies the faulting code; it was not
   reproduced under a debugger.

7. **The port does work the C++ side skips on every indexed mzML it writes, and
   the size of that work is a projection, not a measurement.**
   `src/format/mzml_write_options.rs` computes a real SHA-1 over the output and
   emits it as `fileChecksum`; the C++ emits the literal `0` (`CPP-049`). The
   `sha1` crate is pinned with `features = ["force-soft"]`, so no SHA-NI.
   Throughput of that exact crate and feature was measured on this node in
   wave 4 at 793 MB/s and was **not re-measured here**; the per-tool shares are
   that throughput divided into the output size, an arithmetic projection that
   assumes the cost is purely additive. *No no-hash build has ever been timed*,
   so these are estimates, not an ablation, and they must not be read as "the
   port would otherwise be level".

   | tool | output bytes | projected SHA-1 s | 1-thread gap s |
   |---|---|---|---|
   | BaselineFilter | 3.336 GB | ~4.21 | 9.23 |
   | MzMLSplitter | 1.605 GB | ~2.02 | 4.26 |
   | MapNormalizer | 1.605 GB | ~2.02 | 2.08 |
   | PeakPickerHiRes | 0.536 GB | ~0.68 | 0.13 |

   C++'s zero checksum is arguably the non-conforming side, so this is a
   disclosure, not a correction — but four of the five mzML ratios include work
   only one implementation performs.

8. **The two sides do not write the same number of bytes, and the reason is
   indentation, not metadata.** Rust/C++ output bytes, unchanged from wave 4:
   SpectraFilterWindowMower 0.9434 — the largest deficit of the set —
   PeakPickerHiRes 0.9747, MapNormalizer 0.9906, MzMLSplitter 0.9906,
   BaselineFilter 0.9958; FeatureFinderCentroided 1.0072 the other way;
   DTAExtractor and both FileInfo cases exactly 1.0000. Wave 4's byte census
   attributed 92–93 % of the deficit to XML indentation the port does not write:
   the Rust files contain **zero** tab characters, and C++ writes about 6.3 extra
   lines per spectrum. The metadata differences of §3.6 account for tens of bytes
   per file, not megabytes. The consequence is an output-fidelity fact worth
   stating on its own: *every line of the port's mzML differs textually from the
   source's*, so the two files are not diffable even where all 81,714 or 87,492
   arrays are bitwise identical.

9. **CPU frequency is not pinned.** Governor `schedutil` on all 128 CPUs, boost
   on. The harness records the maximum `scaling_cur_freq` over all 128 CPUs
   after each execution, which is an upper bound on the core that actually ran
   the tool, so this does not by itself prove no cell was clocked down.

10. **SpectraFilterWindowMower has n = 3.** Raw Rust walls: one thread [530.717, 564.434, 585.872], 32 threads [528.420, 537.475, 574.946] — one repetition 10.4 % and 8.8 % above the lowest. The C++ triples are tight (698.469 / 699.092 / 699.517 and 640.425 / 641.690 / 641.783), and foreign load during these cells was *lower* than in wave 4 (0.0012–0.0091 per core against a flat 0.0167), so the scatter is the port's and not the node's. Wave 4 carried the same caveat with a tighter Rust distribution ([515.424, 515.513, 539.252] at one thread). Quote this tool as **0.76–0.84 at one thread and 0.82–0.90 at 32**, never to three digits — which is why the tables print ranges rather than the bare medians 0.807 and 0.838.

11. **FeatureFinderCentroided is measured on a 4,000-spectrum subset, not the
    full run** (`inputs/derived/sub_centroid_velos_50amol_r1_first4000.mzML`,
    91 MB, manifest `MANIFEST_ffc_subset_4000.json`), the same subset waves 3
    and 4 used. Wave 4 established that neither implementation finishes the full
    43,745-spectrum run under a 600–700 s cap; that was **not retested here**,
    and the completed algorithm is slower at one thread than the one wave 4
    tested, so the full run is if anything further out of reach. A single
    full-size cell would cost over an hour. Scaling behaviour on the full run
    may differ from the subset's.

12. **Start-up is under 1 % of every full-size case.** The C++ binaries pay
    63.7–69.5 ms on `-write_ini` alone against the port's 3.1–4.2 ms, a
    difference of 60.6–65.7 ms per tool, which is under 0.7 % of every full-size
    case here — widest MapNormalizer at 32 threads, 0.65 % of its 9.93 s C++
    cell — and near-identical between the two Rust builds (within 0.45 ms on
    every tool; widest MzMLSplitter, 0.443 ms). It matters only on the ~1 s
    FileInfo/featureXML row, where it is 4.2 % of the Rust cell and 6.3 % of the
    C++ one, and there it works *against* the port, whose start-up is the
    cheaper one.

13. **This is a warm-cache CPU benchmark.** `majflt` is 0 in every measured
    repetition: inputs were staged node-local and resident. Fair to both sides,
    excludes I/O, not a cold-start measurement.

14. **The node was reserved for this lane but is not exclusive by
    construction.** ibminode05 was held for the benchmark for the whole run and
    no other lane ran on it. Foreign CPU during the 282 measured repetitions was
    0.0012–0.0199 per core (median 0.0091) and the pre-cell gate value
    0.0016–0.0125 per core, against a flag limit of 0.25 (0.05 at threads = 1) —
    the same figures §3.1 prints, and unchanged to the printed digit whether the
    one failed repetition is included or not. **0 of those 282 repetitions were
    load-flagged, 0 refused and 0 retried.** One execution outside them was:
    a 2.2 ms `-write_ini` start-up baseline at 0.0565 per core (§3.1). It
    feeds no table.
    Exclusive use would still need an admin or a Slurm reservation.

15. **Build configurations differ and can no longer be equalised at all.** C++
    is gcc 14.4 `-O3 -DNDEBUG -mssse3 -ffp-contract=off`, shared libraries, no
    LTO; Rust is opt-level 3, thin-local LTO, codegen-units 16, statically
    linked. Wave 4 could say "neither side uses native CPU flags — both are
    x86-64 baseline, the one thing that could be made equal". **That is no longer
    true**: the port's default build now targets FMA3/AVX and the C++ build does
    not. The `rust-nofma` arm is the closest this run comes to the old
    like-for-like comparison, and §3.7 reports it. This measures two shipped
    configurations, not two compilers.

16. **n < 10 on every large case** (n = 5, n = 4 on one cell after the crash, and
    n = 3 for SpectraFilterWindowMower). None of the large-case ratios is a
    precision measurement.

17. **Instruction counts are exact; wall clock on this workload is not.** The
    validation lane established, with a layout control, that one function
    (`Record::finish`) has two code-generation states 20,612,136 instructions
    apart on this workload and flips between them on source edits to files it
    never calls. Any future comparison of two commits of this port must carry a
    layout control and must not attribute a sub-1 % wall difference to a source
    change.

## 6. What this run does not establish

- **That the flag is what moved FeatureFinderCentroided between wave 4 and
  wave 6.** Three things changed at once — 122 commits, the build flag, and the
  session. §3.7 separates the flag *within this session*, and it is large enough
  to account for most of the move; it does not follow that the algorithm change
  contributed nothing, and the two are not separated against wave 4 (§5.1).
- **The flag's effect on SpectraFilterWindowMower.** Not measured (§5.3).
- **That the flag is free on the seven other cases measured.** "No effect measured at
  n = 5 on this workload" is not "no effect". These are XML- and I/O-bound
  workloads; a float-heavy one might differ (§5.2).
- **Anything about a processor without FMA3.** Every binary here ran on a CPU
  that has it; the refusal path is covered by unit tests over its decision
  function, not by execution (§5.5).
- **Anything about aarch64 or any non-x86_64 target.** The flag is scoped to
  `cfg(target_arch = "x86_64")` and those targets never see it; none was built
  or run here.
- **That `+fma` is the best available choice.** The narrower alternative —
  `#[target_feature(enable = "fma")]` on the three replica functions with
  runtime dispatch — is still unmeasured, as §4.4 said when the decision was
  open. This run measures the flag that was adopted, against its own opt-out,
  and nothing else.
- **Why the C++ FeatureFinderCentroided crashed**, or how often it does. One
  event, no faulting code identified, no reproduction attempt (§3.9, §5.6).
- **Anything about FeatureFinderCentroided at production scale.** The comparison
  is a 4,000-spectrum subset (§5.11).
- **Byte-level output compatibility with the source for any mzML-writing tool.**
  The port writes no indentation, so every line differs textually even where
  every array is bitwise identical (§5.8).
- **Memory parity.** The port uses 1.1×–2.7× the peak RSS on seven of the nine
  cases — six of the eight tools — and less on the other two (§3.4).
- **Anything at better than about 3 % across sessions** (§5.2), and nothing at
  all about cold-start or I/O-bound behaviour (§5.13).
- **A decomposition of the mzML write path.** The SHA-1 share is projected from
  a throughput measurement taken in wave 4, not ablated and not re-measured
  (§5.7).
- **Thread scaling beyond 1 and 32**, and nothing about the five Rust tools that
  ignore `-threads` beyond the fact that they ignore it (§3.5).
- **Any claim that the C++ 32-thread column had a doubled compute budget.** It
  did not; the surplus threads are a library pool the harness itself sizes
  (§3.5).

## 7. Reproducing

```sh
# C++ reference: already built and staged; verify against its manifest
sha256sum /ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576/BUILD_MANIFEST.json

# the two Rust builds of this run, from one clean export of 36c26a0.
# The build scripts are archived beside the run config, as wave 5's are (§4):
# /ceph/ibmi/abi/oliver/bench/openms4/w6-2026-09-18/build/ holds build_rust_orig.sh,
# build_rust_flags.sh, build_both.sh, flags.diff (the one added line:
# `if [ -n "${BUILD_RUSTFLAGS:-}" ]; then export RUSTFLAGS="$BUILD_RUSTFLAGS"; fi`),
# the two build logs and SHA256SUMS. /scratch is node-local and is not preserved.
BLD=/ceph/ibmi/abi/oliver/bench/openms4/w6-2026-09-18/build
B=/scratch/$USER/bench-w6
$BLD/build_rust_orig.sh  openms-rs-36c26a0a…tar.gz $B/default   # no RUSTFLAGS -> +fma from .cargo/config.toml
BUILD_RUSTFLAGS="-C target-feature=-fma" \
  $BLD/build_rust_flags.sh openms-rs-36c26a0a…tar.gz $B/nofma   # the documented opt-out

# the wave-6 matrix: the harness is unmodified, the config is this run's own
H=/ceph/ibmi/abi/oliver/bench/openms4/harness
C=/ceph/ibmi/abi/oliver/bench/openms4/w6-2026-09-18/config
R=/ceph/ibmi/abi/oliver/bench/openms4/results/2026-09-18-w6refresh
python3 $H/bench.py run --plan $C/w6_profile.json --impls $C/impls.json --results $R --run-id w6-profile
# …likewise w6_centroid, w6_ffc (with --manifest …/MANIFEST_ffc_subset_4000.json), w6_sfwm
python3 $C/consolidate_w6.py && python3 $C/tables_w6.py

# which build a binary is: 41 vfma* instructions is the +fma FeatureFinderCentroided,
# 2 is the opt-out; any tool with thousands of %ymm operands is a +fma build
objdump -d $B/default/36c26a0a301d/bin/FeatureFinderCentroided | grep -cE 'vfm(add|sub|nmadd|nmsub)'
```

Results live under `/ceph/ibmi/abi/oliver/bench/openms4/results/`, indexed by
`results/RUNS.tsv`. The wave-6 run is `2026-09-18-w6refresh`; each sub-run's
`summary.md` is the harness's own output and carries the full per-case detail,
and `w6_consolidated.json` carries the merged statistics the tables above print.
The oracle artifacts cited here are registered with their sha256 in
`SOURCE_PROVENANCE.json`.
