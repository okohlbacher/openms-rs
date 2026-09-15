# Benchmarks: this port against the C++ OpenMS4 Release build

What is compared, how, and what the first numbers do and do not say. Nothing
here is a claim about the port being faster or slower in general: one tool has
been measured on one instrument-scale input, once per implementation, on a
shared node. The method section exists so the next run is comparable; the
results section is deliberately small.

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

- `bench.py` runs cases (tool × dataset × implementation × threads 1/16/full)
  with one warm-up and R repetitions, default 5. Inputs and the C++ prefix are
  staged to node-local `/scratch`; outputs go there too; `TMPDIR` and
  `OPENMS_TMPDIR` are on `/dev/shm`. Both implementations get `-threads N`
  **and** `OMP_NUM_THREADS = RAYON_NUM_THREADS = N`.
- Measurement is `os.wait4`: wall, user, sys, peak RSS, page faults, I/O, exit
  code or signal, and the output sha256.
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
unnecessary. The node is **not** exclusive: another user's job and an idle
GitHub Actions runner were present during the 2026-09-14 smoke timing. CPU
frequency is not pinned.

## 3. What has been measured so far

### 3.1 PeakPickerHiRes on the 2.3 GB Q Exactive profile run

The only instrument-scale Rust-vs-C++ comparison that exists. Oracle directory
`../oracle/topp-peak-picker-scale` (manifest sha256
`d92c8993b35604fa15d06454c60fff46aefa6186dc2f549c259eca179ebd20e7`), driver
`bench_06.sh`.

Input `profile_hr_qe_silac_uk222/UK222.mzML`, sha256
`bd6f6e19…dc167`, 2,317,975,830 bytes, 40,856 spectra (6,911 MS1 profile,
33,945 MS2 centroid), 1 chromatogram, 197,765,338 raw points. INI from the C++
`-write_ini`, sha256 `dc0f7b60…83ff`, defaults: `threads 1`, `processOption
inmemory`, `algorithm:ms_levels` empty (automatic), `signal_to_noise 0.0`,
`report_FWHM false`. Rust binary built from `bundle/P3-PICKER-TOOL` merged with
`fix/picker-scale`, `--release --locked --no-default-features --features
mzml,paramxml`.

| | C++ Release | Rust |
|---|---|---|
| exit | 0 | 0 |
| wall | 26.80 s and 28.85 s in the two runs | 38.69 s |
| peak RSS | 3,977,240 and 3,977,188 KiB = **3,884 MiB** | 4,412,352 KiB = **4,309 MiB** |
| output | 549,528,678 B | 535,615,957 B |
| stdout | `MS-level 1: 6911 / 6911`, `MS-level 2: 0 / 33945` | identical (the C++ adds progress logging and a timing line this port does not write) |

So the port is about **1.4× slower** and uses about **11 % more** resident
memory on this input, at one thread.

**Equivalence, data and metadata judged separately** (the decoded D6
comparison, `decoded_compare_probe.rs`, reproduced independently by the
verifier):

- 40,856 spectra on both sides; native ids, MS levels and peak counts equal.
- **22,776,198 centroids, every m/z and every intensity bit-identical.** This
  is the result that matters, and it is exact, not within a tolerance.
- Spectrum retention times differ in 10,671 of 40,856 records, by at most
  9.09e-13 s — 1 to 2 ULP of `f64` on retention times of 60 to 4,400 s, after a
  round trip through each implementation's own writer. This one is writer-side
  and the port's text is the correct side: the two readers were shown to agree
  bit for bit on all 40,856, while the C++'s own output fails to reparse to its
  own stored `double` in exactly those 10,671 records, because its XML writers
  print `writtenDigits<double>()` = 15 significant digits, which is `digits10`
  and not the `max_digits10` = 17 a round trip needs (`CPP-307`).
- **The picked TIC chromatogram is bit-identical.** 8,174 points on both
  sides, every retention time and every intensity equal. This run first showed
  it differing — 8,173 of 8,174 retention times by at most 3.19e-3 s, and 7,891
  of 8,174 intensities by more than 1e-6 relative, 71 by more than 1e-3, worst
  1.75e-3 (252,284,752 against 252,727,072 at 4,393.5 s) — while every spectrum
  centroid agreed, and the retained workflow-2 fixture (five short
  chromatograms) was bit-exact, so only real data showed it. Lane
  `fix/picked-chromatogram` (merged as `b1700de`, after the run above)
  root-caused it in the **mzML reader**, not in the picker:
  `MzMLHandlerHelper::decodeBase64Arrays` applies the minute multiplier of an
  `MS:1000595` time array in place, and for a 32-bit array the element is a
  `float&`, so the `double` product is narrowed back to `float`
  (`MzMLHandlerHelper.cpp:217-222`; the 64-bit and Numpress paths are not —
  `CPP-306`). This input's TIC time array is exactly that case, 32-bit in
  minutes, and the picker's spline apex amplifies the at-most-2.44e-4 s
  per-point difference into the numbers above.
  `ReadOptions::source_time_array_precision`, which `ReadOptions::source()`
  sets and every tool load path therefore uses, reproduces the narrowing; the
  library default keeps the full `f64` product. With it the two picked
  chromatograms agree bit for bit over the whole file. The wall and peak-RSS
  numbers in the table above predate that change and were not re-measured; it
  touches one multiplication per time element of one 40,856-point array, so no
  re-measurement is claimed either way.
- The C++ `FuzzyDiff` never reaches the data on this pair: it fails at line 1,
  column 31 on the XML declaration encoding (`ISO-8859-1` against `UTF-8`),
  the documented container difference of decision D6.

**What this measurement is not.** It was run by `bench_06.sh`, not by
`bench.py`, and therefore:

- on **ibminode06**, the build node, not the idle timing node ibminode05, with
  foreign load of 24 to 35 runnable tasks throughout. The manifest says so
  itself: wall times are indicative; peak RSS is not affected by foreign load;
- **once per implementation**, so there is no median, no IQR and no confidence
  interval. The two C++ runs, 2 s apart in wall time, are the only spread
  information available;
- with `OMP_NUM_THREADS` unset. `-threads 1` reaches the OpenMP team through
  `TOPPBase`, but libgomp still spins up a thread per processor at start-up,
  which cost about 6 % of wall on a 1.8 s job in the review's own measurement.
  The `bench.py` runs set `OMP_NUM_THREADS = N`, so C++ numbers from the two
  paths are not directly comparable;
- output sizes differ by 13.9 MB because the two writers make different
  container choices, not because either dropped data (see the work-parity note
  below).

### 3.2 The 2026-09-14 smoke run: five tools, no valid Rust numbers

Recorded for completeness because it is the only run through the full harness.
On the full-size inputs, **60 of 60 Rust executions failed before doing any
work**: four tools exited 6 on `startTimeStamp="-infinity"` and BaselineFilter
exited 3 on the reader's peak-count ceiling. The C++ column of that run is
valid (DTAExtractor 57.5 s, MzMLSplitter 14.6 s, SpectraFilterWindowMower
274.2 s, MapNormalizer 15.1 s at `-threads 1`; BaselineFilter 15.9 s on
UK222), and it is the reason the reader, writer, thread and morphology lanes
of wave 3 exist. Both blockers are now fixed (`fix/mzml-reader-scale`,
`fix/mzml-writer-scale-parity`), and the run has **not** been repeated.

The sub-limit supplement of that day — ratios on 0.16 s to 1.9 s runs over
first-N-spectra slices — should not be quoted. Its review showed the headline
("Rust faster on 4 of 5 tools") was a start-up effect: with the fixed cost
subtracted (C++ about 65 ms before any mzML work, Rust about 3 ms), it becomes
2 faster, 1 equal, 2 slower.

## 4. Caveats that apply to every number here

1. **The harness peak-RSS figure was wrong for small processes, and is
   partly fixed.** `bench.py` starts tools with `subprocess.Popen`, which uses
   `vfork`; on `exec`, Linux copies the parent's peak RSS into the child's
   `ru_maxrss`. The recorded value was therefore `max(harness Python peak, tool
   peak)`: `/bin/true` reported 25,600 KiB against 1,024 KiB from
   `/usr/bin/time`, and after the harness process had peaked at 400 MiB,
   `/bin/true`, a 30 MiB Rust tool and a 72 MiB C++ tool all reported 430,340
   KiB. The smoke run's "Rust RSS 48,172 KiB" in all 48 failing runs was this
   floor; the real value was 30,720 KiB. C++ values at 1.5–3.2 GB are
   unaffected, and the §3.1 numbers come from `/usr/bin/time -v` directly, not
   from the harness. Any future harness run must measure RSS through a small
   intermediary (`/usr/bin/time -f %M`, a cgroup v2 `memory.peak`, or
   `/proc/<pid>/status` `VmHWM` before reaping) and record `/bin/true` through
   the launcher as a floor check.
2. **Start-up dominates short runs.** Fixed cost per process: C++ about 120 ms
   on a 10-spectrum input and 65 ms for `-write_ini` alone, Rust about 95 ms
   and 3 ms. Never headline a sub-second run, and report a start-up baseline
   next to every ratio.
3. **Rust tools are serial where C++ is OpenMP-parallel.** No ported tool
   passes a thread policy to its algorithm yet: `pick_experiment` takes none,
   and `FeatureFinderCentroided`'s ported stage is serial in the source as
   well. `-threads` is wired through `ToolContext::in_thread_pool` and is
   asserted to change no output byte, but it changes no work either. So the
   only like-for-like comparison today is `-threads 1`; any 16-thread or
   full-machine row must be labelled "Rust serial against C++ parallel". C++
   short jobs at 16 threads are in fact *slower* than at 1 thread (32 threads,
   CPU utilisation 4.6–12).
4. **Equivalence must be judged on data and metadata separately.** A verdict
   that stops at the first difference says nothing: on three of the smoke
   tools the harness reported DIFFERENT while every decoded peak array was
   bitwise identical, and the real difference was the mzML root element. The
   §3.1 comparison is split this way on purpose. Whole-line FuzzyDiff
   whitelists are worse than useless for this: a `startTimeStamp` whitelist
   drops the entire `<run>` line, and a featureXML `id=` whitelist drops every
   line containing `id=`. FuzzyDiff's 1e-6 relative default is 1 ppm on m/z,
   far too loose to call centroids equal — which is why §3.1 asserts bit
   identity instead.
5. **Work parity.** The two writers do not do the same work. This port writes
   `indexedmzML` with real record offsets and a real SHA-1 `fileChecksum`; the
   C++ writer writes the constant `0` (CPP-049) and an `indexListOffset` one
   byte early (CPP-305). The digest costs about 10.7 s on a 3.34 GB document at
   roughly 310 MB/s, because `sha1` is pinned with `force-soft`. Any writer
   timing must either enforce output parity or state this difference; §3.1 is
   a picker measurement whose output is written by both sides, so the digest is
   inside the Rust column.
6. **Build flags are not equalised beyond "no native CPU flags".** C++ is gcc
   `-O3 -mssse3 -ffp-contract=off`, shared libraries with PLT calls, no LTO,
   conda-forge zlib. Rust is the default release profile: opt-level 3,
   codegen-units 16, thin-local LTO, x86-64 baseline, `miniz_oxide`, a static
   binary. Neither uses `-march=native`.
7. **Uncertainty is understated by the within-session CI.** The bootstrap CI
   reflects noise inside one session only. An independent interleaved rerun
   moved medians by 0.3 % to 2.1 %, outside the reported IQRs of 0.1 % to
   0.5 %, while ratios reproduced to about 1 %. Treat effects below about 3 %
   as not measured.
8. **The node is shared and its frequency is not fixed.** Ask for exclusive
   use of ibminode05, or a Slurm reservation, before quoting anything at
   percent precision.

## 5. Not measured yet

- **FeatureFinderCentroided at scale.** The wrapper and the algorithm are
  ported and equivalent on `TOPP_FeatureFinderCentroided_1` (8 features,
  bit-identical m/z and hull points, fitted fields within 1e-9 relative
  against the Release build), but no instrument-scale run exists. A 1 GB
  centroided input may take tens of minutes per C++ execution, so size it with
  `--reps 1` first.
- **FileInfo at scale.** Only C1-sized inputs have been run.
- Every other tool: the five smoke tools need the full-size run repeated now
  that the reader and writer blockers are fixed, and
  `SpectraFilterWindowMower` additionally needs its own input-point ceiling
  lifted (`src/processing/window_mower.rs`), which stops it at 5,000 real
  spectra.
- Thread scaling of anything, for the reason in caveat 3.

## 6. Reproducing

```
# C++ reference: already built and staged; verify against its manifest
sha256sum /ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576/BUILD_MANIFEST.json

# the instrument-scale picker comparison of section 3.1
bash ../oracle/topp-peak-picker-scale/sync_build_06b.sh   # build the Rust tool on the node
bash ../oracle/topp-peak-picker-scale/bench_06.sh run2    # both tools under /usr/bin/time -v

# a harness run (medians, IQR, load gate, equivalence ladder)
H=/ceph/ibmi/abi/oliver/bench/openms4/harness
python3 $H/bench.py run --plan $H/plans/smoke.json --run-id <date>-smoke
python3 $H/bench.py summarize --run-id <date>-smoke
```

Results live under `/ceph/ibmi/abi/oliver/bench/openms4/results/`, indexed by
`results/RUNS.tsv`. The oracle artifacts of section 3.1 are registered with
their sha256 in `SOURCE_PROVENANCE.json`.
