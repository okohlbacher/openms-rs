# TOPP `-threads` support

Lane `tool-threads`. Rust: `src/cli/context.rs` (`ToolContext::thread_policy`,
`ToolContext::in_thread_pool`, `ToolContext::THREAD_CEILING`,
`ToolContext::WORKER_STACK_BYTES`) and the `Tool::run` of the five wave-2 tools
under `src/cli/tools/`, plus `PeakPickerHiRes`, which takes the policy through
`ToolContext::thread_policy` and scopes its own pool to the picking call rather
than wrapping its body (see "The picker is the sixth tool" below). Test:
`tests/topp_threads.rs`. Source: OpenMS4-cli
c19e494 `source/APPLICATIONS/TOPPBase.cpp` (sha256
`326b96f85b4041febec49e18252d7be2ad0a5450abe5249adecc28ca315fdc05`) and
`include/OpenMS/APPLICATIONS/TOPPBase.h`.

## API mapping

| Source | Rust |
|---|---|
| `registerIntOption_("threads", "<n>", 1, ...)` (`TOPPBase.cpp:168`) | `register_common` in `src/cli.rs` (unchanged) |
| `getParamAsInt_("threads", 1)` (`TOPPBase.cpp:408`) | `ToolContext::threads` |
| `static void setMaxNumberOfThreads(int num_threads)` (`TOPPBase.h:119`, `TOPPBase.cpp:84-98`) | `ToolContext::thread_policy`, the worker count |
| `omp_set_num_threads(num_threads)` before `main_` (`TOPPBase.cpp:408-415`) | `ToolContext::in_thread_pool`, called by every tool's `Tool::run` |
| `omp_get_num_procs()` | `Threads::all()`, that is `std::thread::available_parallelism` |

### How a tool uses it

```rust,ignore
impl Tool for MyTool {
    fn run(ctx: &ToolContext) -> Result<ExitCode> {
        ctx.in_thread_pool(|| Self::run_in_pool(ctx))?
    }
}
impl MyTool {
    fn run_in_pool(ctx: &ToolContext) -> Result<ExitCode> {
        // the body; parallel algorithms take ctx.thread_policy()
    }
}
```

A tool that overrides `Tool::run_io` computes inside the closure and writes to
`out` and `err` after `in_thread_pool` returns, because those streams are not
`Send`.

## Preserved source conventions

- A positive count `n` gives `n` workers, from `-threads` or from the INI
  `threads` value; the command line wins over the INI file, as for every
  parameter.
- Zero **and every negative count** give every available processor
  (`if (num_threads <= 0) num_threads = omp_get_num_procs()`). This replaces
  the earlier Rust behaviour, where `thread_policy` mapped a negative count to
  one worker through `Threads::from_cli`.
- `OMP_NUM_THREADS` does not size the tool body's workers. The source's
  `omp_set_num_threads` overrides it.
- The registered default stays 1, so a run without `-threads` is serial.
- The source `@note` that the setting only works when OpenMS is compiled with
  OpenMP carries over: without the `parallel` feature the body runs on the
  calling thread and every computation is serial.

## Native differences

- **Scoped pool instead of a process-wide limit.** The body runs on a rayon
  pool of exactly `thread_policy().get()` workers, built for the call and
  ended when it returns. Rayon inside the body reports that count, and rayon
  parallel iterators use this pool rather than rayon's global one. A pool is
  built for one worker too, so `-threads 1` bounds rayon to one worker.
- **Worker names.** Workers are named `openms-<index>`, visible as
  `/proc/<pid>/task/*/comm`.
- **Stack size.** Workers get 8 MiB (`WORKER_STACK_BYTES`), the main-thread
  stack on Linux and macOS, because the body moves off the main thread.
- **Ceiling.** A count above both `THREAD_CEILING` (1024) and the available
  processors is clamped to the larger of the two. The source has no bound;
  libgomp tries to create every requested thread and aborts when it cannot.
  Clamping never changes a result, by the determinism contract.
- **Thread-start failure.** When the operating system refuses the worker
  threads, the run ends with `Error::Io`, exit code `UNKNOWN_ERROR`, before any
  of the body ran.
- **cgroup quotas.** `available_parallelism` honours a Linux cgroup CPU quota;
  `omp_get_num_procs` counts only the affinity mask. Under a quota the port can
  start fewer workers for `-threads 0`.
- **`RAYON_NUM_THREADS`** is ignored, because the pool size is always explicit.
- **Start-up threads.** Without `OMP_NUM_THREADS`, the C++ tools start one
  extra thread per processor at library load, even for `-write_ini`. gdb traced
  them to OpenBLAS `blas_thread_init` in the conda-forge `liblapack.so.3` that
  the Release build links, not to OpenMS. They spin for about 0.1 s of CPU each
  and do no tool work. The port links no BLAS and starts no such threads.

## Determinism

The determinism contract of `src/concept/parallel.rs` applies to the tool
body: output bytes must not depend on the worker count. Under `-test` the
outputs are byte-identical at every count. Outside `-test`, the processing
record lists every parameter, so the `threads` value appears in the output as
given, next to the completion time. The C++ outputs differ across thread counts
in exactly those two values (benchmark smoke run 2026-09-14). A harness that
checks thread invariance byte for byte must therefore run with `-test` or
ignore those two values, for both implementations.

## Checked boundaries and evidence

Executed C++ (tier 1 probe): the Release build
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` on
ibminode06, with `/proc/<pid>/task` sampled every millisecond. The driver,
cases and results live outside the repository in
`../oracle/tool-threads/`:

| File | sha256 |
|---|---|
| `threads_probe.py` | `48dfaaef32adabdc42f157897fdfc13b0b01a6893d62b616012464ea77c4fde8` |
| `make_cpp_cases.py` | `55285cd7443389f6227774530700cae26d6c314f73807d38fffe1653243f54d1` |
| `cpp_cases.json` | `3e653c95ac2b6623555bad54ce4680686a76a4fea239d1d660b79649859ce770` |
| `cpp_results_ibminode06.jsonl` | `46f7d3ae8dd3b555bb0a016bc0925d33f7dd3aa0aa8d67522b6af56bbdcb2698` |
| `make_pair_cases.py` | `f3d620de3c8b19b11894e48d49341ff4fa417b51857726a165be3668517db5ab` |
| `pair_cases.json` | `4fe138416b14873ff3340e37dccc5c1a3c587d28c223e9cb91eb57c38c279a46` |
| `pair_results_dax.jsonl` | `5edd1b4131e8a7b9ac44965c816851810fd12f184a3048ee647b36aa261d6b48` |

MapNormalizer on `inputs/derived/sub_centroid_uk222_picked_first5000.mzML`,
ibminode06 (128 logical processors):

| Case | Tasks | Busy tasks | Reading |
|---|---|---|---|
| `-threads 1`, `OMP_NUM_THREADS=1` | 2 | 1 | one idle helper thread |
| `-threads 1`, no `OMP_NUM_THREADS` | 129 | 128 | OpenBLAS start-up threads at about 0.1 s CPU each |
| `-threads 4`, `OMP_NUM_THREADS=1` | 5 | 4 | `-threads` overrides the variable |
| `-threads 4`, `OMP_NUM_THREADS=16` | 20 | 19 | 15 OpenBLAS start-up threads plus a 4-thread body team |
| `-threads 0`, `OMP_NUM_THREADS=1` | 129 | 128 | `omp_get_num_procs()` |
| `-threads -1` and `-threads -7`, `OMP_NUM_THREADS=1` | 129 | 128 | negative means all processors |
| `-threads 0`, `taskset -c 0-3` | 5 | 4 | affinity mask honoured |
| `-write_ini`, no `OMP_NUM_THREADS` | 128 | 96 | start-up threads without any tool work |

Native tests (`tests/topp_threads.rs`):

- `policy_follows_toppbase_semantics`, `policy_clamps_only_absurd_counts`,
  `ini_threads_value_sets_the_policy`: the worker count for default, positive,
  zero, negative, INI and absurd values.
- `body_runs_inside_a_pool_of_the_policy_size` (feature `parallel`):
  `rayon::current_num_threads()` inside the body equals the count, and the body
  runs on an `openms-<index>` worker.
- `serial_build_runs_the_body_on_the_calling_thread` (no `parallel`).
- `body_errors_pass_through_the_pool`: exit codes and diagnostics of errors
  raised inside the body are unchanged.
- `outputs_are_byte_identical_across_thread_counts`: the **six** tools of
  `TOOLS` under `-test` at 1, 2, 8 and all workers. `PeakPickerHiRes` joined in
  `3b943e4`, and it is the first tool for which this assertion is about a
  parallel body rather than a serial one.
- `outside_test_mode_only_the_provenance_records_the_thread_count`.
- Linux, the real executables sampled through `/proc`:
  `every_executable_runs_its_body_on_the_requested_pool` (exactly `n`
  `openms-*` workers and `n + 1` tasks for `-threads 1`, `2` and `4`, with
  byte-identical outputs) — with **one stated exception**: `PeakPickerHiRes`
  starts `0` workers at `-threads 1`, by design, because its pool exists only
  around the picking call and is not built for a single worker. That is why
  `expected_workers` special-cases it and why the picker is sampled on a
  profile input (`picking_input`), so the sampled region is one the tool
  actually parallelises. Also
  `zero_and_negative_counts_start_every_available_processor`,
  `thread_environment_variables_do_not_size_the_pool`,
  `ini_threads_value_sizes_the_pool`.
- `#[ignore]`d, IBMI nodes only (`cargo test --release --test topp_threads --
  --ignored`): `hpc_benchmark_slices_are_thread_invariant` (the benchmark's
  600- and 5000-spectrum UK222 slices, all six tools at 1, 2 and 16 workers)
  and `hpc_full_size_inputs_are_thread_invariant` (PXD001819 50amol_R1, 1.2 GB,
  and UK222, 2.3 GB). Both check the pool size, and that exit status,
  diagnostics and outputs agree across counts. On the full-size inputs the
  tools still exit early (6 and 3) on the reader blockers other lanes own.

Rust release binaries against C++ on dax (384 logical processors, load about
0.15 per processor), median of 3 runs:

| Tool, input | `-threads` | Rust tasks | Rust wall (s) | C++ tasks | C++ wall (s) |
|---|---|---|---|---|---|
| DTAExtractor, first 5000 | 1 | 2 | 1.168 | 2 | 1.567 |
| DTAExtractor, first 5000 | 16 | 17 | 1.177 | 32 | 1.565 |
| DTAExtractor, first 5000 | 0 | 385 | 1.216 | 512 | 2.376 |
| MapNormalizer, first 600 | 1 | 2 | 0.164 | 2 | 0.169 |
| MapNormalizer, first 600 | 16 | 17 | 0.167 | 32 | 0.176 |
| MapNormalizer, first 600 | 0 | 385 | 0.214 | 512 | 0.404 |

The `-threads 1` and `16` rows set `OMP_NUM_THREADS` to the same value; the
`0` rows leave it unset. The pool costs about 3 ms at 16 workers and about
50 ms at 384. Those five tools are serial, so the extra workers stay idle — and
they still are: the wave-4 benchmark sampled all eight executables on full-size
data and found `BaselineFilter`, `DTAExtractor`, `MapNormalizer`, `MzMLSplitter`
and `SpectraFilterWindowMower` at CPU utilisation 1.00 with a 33-thread pool
alive at `-threads 32`, and `FileInfo` building no pool at all.

## The picker is the sixth tool

`PeakPickerHiRes` (`3b943e4`, wave 4) is the first tool whose pool does work.
It uses the same `ToolContext::in_thread_pool`, but **around the picking call
rather than around the whole body** (`src/cli/tools/peak_picker_hi_res.rs:249`),
and it skips the pool entirely when the policy is one worker (`:246-248`). Two
measured reasons:

- wrapping the body runs the whole tool, including the gigabyte-scale read and
  write, on a pool worker, which costs **0.68 s** on a 2.3 GB run through
  glibc's per-thread arenas;
- the picker has exactly one parallel region, and a grep of `src/format`,
  `src/kernel`, `src/metadata` and the other CLI tools confirms none of them
  contains `rayon`, so there is no stray global-pool `par_iter` for the body
  wrapper to bound.

A third reason is structural rather than measured: this tool overrides
`run_io`, which is what `run_with` calls, so a wrapper placed in `Tool::run`
would never execute for it. Inside the picking call the parallel work is
distributed by `BatchWorkers` in `src/processing/peak_picking.rs`, which exists
because `concept::parallel::map_collect` builds a `ThreadPoolBuilder`
unconditionally with no ambient-pool check.

The trade-off is a framework question, not a picker one, and it is carried
forward in [the work packages](EARLY_TOPP_WORK_PACKAGES.md#wave-4-status):
either `in_thread_pool` grows a documented one-worker fast path, or its doc
records that a tool scoping its pool to its parallel region may skip it.

The determinism contract is unchanged and is met: the picker's output is
byte-identical at 1, 2, 4, 8 and 32 workers and equal to the serial result. On
the 2.3 GB benchmark input it turns 25.367 s into 12.592 s (2.01x) where the
C++ picker gains 1.06x; see [BENCHMARKS](BENCHMARKS.md) §3.5.
