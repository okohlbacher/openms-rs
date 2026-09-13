# Wall-clock and process CPU timing

`system::stop_watch` covers `SYSTEM/StopWatch.h` and `SYSTEM/StopWatch.cpp` at
SDK `bc9cc12`. The class measures wall, user and kernel time of the current
process with start / stop / resume semantics, and formats a duration with only
the units it needs. Its five direct TOPP consumers are `Decharger`, `Epifany`,
`FeatureFinderLFQ`, `MetaboliteAdductDecharger` and `ProteinInference`.

## API mapping

Every public member of the header, plus the private `TimeDiff_` surface the port
exposes, with its Rust counterpart.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `StopWatch()` | `StopWatch::new`, `StopWatch::default` | Stopped, everything zero. |
| `StopWatch(const StopWatch&)`, `operator=` | `Clone`, `Copy` | Whole-value copy; the type owns no resources. |
| `~StopWatch()` | — | Implicit and empty on both sides; nothing to release. |
| `void start()` | `StopWatch::start() -> Result<()>` | `Exception::Precondition` → `Error::InvalidValue`. Clears first, so it cannot resume. |
| `void stop()` | `StopWatch::stop() -> Result<()>` | `Exception::Precondition` → `Error::InvalidValue`. |
| `void resume()` | `StopWatch::resume() -> Result<()>` | `Exception::Precondition` → `Error::InvalidValue`. |
| `void reset()` | `StopWatch::reset()` | Infallible here; see *Native differences*. |
| `void clear()` | `StopWatch::clear()` | |
| `double getClockTime() const` | `StopWatch::clock_time() -> f64` | Always available. |
| `double getUserTime() const` | `StopWatch::user_time() -> Option<f64>` | `None` where the platform reports no user/kernel split. |
| `double getSystemTime() const` | `StopWatch::system_time() -> Option<f64>` | `None` where the platform reports no user/kernel split. |
| `double getCPUTime() const` | `StopWatch::cpu_time() -> Option<f64>` | `user + system` where the split exists, else the platform aggregate. |
| `bool isRunning() const` | `StopWatch::is_running() -> bool` | |
| `bool operator==(const StopWatch&) const` | `PartialEq` (derived) | Compares accumulated interval, last start reading and running flag, as the source does. |
| `bool operator!=(const StopWatch&) const` | `PartialEq` (derived) | |
| `bool operator<(const StopWatch&) const` | `StopWatch::cpu_time_cmp -> Option<Ordering>` | `Ordering::Less`. Not `PartialOrd`; see *Native differences*. |
| `bool operator<=(const StopWatch&) const` | `cpu_time_cmp != Some(Greater)` | The source defines it as `!(other < self)`. |
| `bool operator>=(const StopWatch&) const` | `cpu_time_cmp != Some(Less)` | The source defines it as `!(self < other)`. |
| `bool operator>(const StopWatch&) const` | `cpu_time_cmp == Some(Greater)` | The source defines it as `other < self`. |
| `std::string toString() const` | `StopWatch::summary() -> String` | Same four components, same order, same separators. |
| `static std::string toString(double)` | `StopWatch::format_seconds(f64) -> Result<String>` | Errors only where the source's cast is undefined. |
| `TimeType` (`clock_t` / `UInt64`) | — | Not ported: the port stores microseconds, not platform ticks. |
| `SecondsTo100Nano_` | — | Not ported: Windows-only tick scale, unreachable here. |
| `cpu_speed_` (`sysconf(_SC_CLK_TCK)`) | `MICROS_PER_TICK` (private) | Fixed at `USER_HZ` = 100 for the `/proc` interface glibc also answers with. |
| `TimeDiff_` | `TimeSample` | Public here, so the arithmetic is testable in isolation. |
| `TimeDiff_::user_ticks`, `kernel_ticks` | `CpuSample::Split { user_micros, kernel_micros }` | Microseconds, not ticks. |
| `TimeDiff_::start_time`, `start_time_usec` | `TimeSample::wall_seconds`, `wall_micros` | Same seconds-plus-microseconds split. |
| `TimeDiff_::userTime()` | `TimeSample::user_time() -> Option<f64>` | |
| `TimeDiff_::kernelTime()` | `TimeSample::kernel_time() -> Option<f64>` | |
| `TimeDiff_::getCPUTime()` | `TimeSample::cpu_time() -> Option<f64>` | |
| `TimeDiff_::clockTime()` | `TimeSample::clock_time() -> f64` | |
| `TimeDiff_::operator-` | `TimeSample::difference` | |
| `TimeDiff_::operator+=` | `TimeSample::sum` | |
| `TimeDiff_::operator==` | `PartialEq` (derived) | |
| `TimeDiff_::ticksToSeconds_` | folded into the sampler | Ticks become microseconds at sample time. |
| `StopWatch::snapShot_()` | `TimeSample::now()` | Public here; it is the only injection point a test needs. |
| `accumulated_times_`, `last_start_`, `is_running_` | private fields of `StopWatch` | Same three, same meaning. |

Native additions: `MAX_FORMATTABLE_SECONDS`, `CpuSample::{user_micros,
kernel_micros, total_micros}`.

## Preserved source conventions

**The state machine is the source's.** `start` on a running watch, `stop` on a
stopped one and `resume` on a running one are all refused. `start` *clears
first*, so `start`, `stop`, `start` cannot resume a measurement — the header says
so explicitly and the port does the same. `reset` zeroes the accumulated total
and keeps a running watch running; `clear` zeroes it and stops.

**Reading order is the source's.** `clock_time` on a running watch converts the
accumulated total and the current interval to `f64` *separately* and then adds
them, rather than normalising first; floating-point addition is not associative,
so the order matters and it is kept. `cpu_time` is literally
`user_time() + system_time()`, each of which takes **its own clock sample**,
exactly as `getCPUTime()` calling `getUserTime()` then `getSystemTime()` does.

**The carry quirk is preserved.** `TimeDiff_::operator+=` carries a whole second
only while the microsecond remainder is *strictly greater* than one million, so
a remainder of exactly `1 000 000` is left unnormalised. It cannot change any
reported duration — a million microseconds contribute the same second either way
— but it *is* visible through equality, so a port that silently normalised would
make two watches compare equal that the source considers different. The borrow
rule in `operator-` (`while (usec < 0)`) is preserved in the same way. Both loops
are expressed as one integer division instead of one iteration per second, which
is the same result in constant rather than unbounded time.

**Equality is state equality.** `operator==` compares the accumulated interval,
the last start reading and the running flag — not elapsed times. Two watches
started at different moments are unequal even if they have accumulated nothing.

**`toString(double)` is transcribed.** Days are unpadded, hours / minutes /
seconds are zero-padded to two digits, and the fallback branch prints the
*unrounded* argument with `snprintf("%.2f")` — which is why `59.999` renders as
`60.00 s` while `60.0` renders as `01:00 m`. A negative duration cannot make any
unit branch positive, so it always takes the seconds form.

Rust's `{:.2}` replaces that `snprintf`. Both round the exact binary value to
nearest, ties to even, so they agree on the values where a naive
round-half-away implementation would not: `0.125` → `0.12`, `0.025` → `0.03`,
`2.675` → `2.67`, `1.005` → `1.00`, and `-0.004` → `-0.00` with the sign kept.
That was checked by running the two formatters side by side rather than
reasoned about.

## Native differences

**Wall time is monotonic.** The source samples `gettimeofday`, the realtime
clock: an NTP correction or an administrator changing the date is measured as
elapsed wall time, and can make an interval negative. This port samples
`std::time::Instant`, so an interval is never negative and never absorbs a clock
adjustment. Absolute readings are expressed relative to a process-wide monotonic
origin rather than the Unix epoch; only differences are reported, so the origin
is not observable.

**CPU time is an `Option`, and the split is Linux-only.**

| Platform | `clock_time` | `user_time` / `system_time` | `cpu_time` |
|---|---|---|---|
| Linux | real | real, `/proc/self/stat` `utime`/`stime` | `user + system` |
| macOS, other Unix, Windows | real | `None` | aggregate process CPU clock |
| anything else | real | `None` | `None` |

The source reaches the per-process user/kernel split through the POSIX `times`
system call, which works on macOS too. This port cannot call it without adding a
libc dependency, and dependencies are the integrator's decision, so on Linux it
reads `/proc/self/stat` and elsewhere it reports the split as **unavailable**
rather than as zero. `/proc/self` names the thread group, so `utime` and `stime`
already aggregate over every thread of the process — which is the contract the
header states. Where only an aggregate exists, `cpu_time` still answers, from
the `cpu-time` crate that is already a dependency of this crate.

**Ticks become microseconds at sample time.** The source accumulates integer
ticks and divides once by `sysconf(_SC_CLK_TCK)`. This port multiplies each tick
by 10 000 µs — exact integer arithmetic — accumulates microseconds and divides
once by 1e6. Both compute the correctly rounded `f64` nearest to the same
rational `n/100`, so the two agree **bit for bit**, not merely within a
tolerance. `USER_HZ` is fixed at 100 for the `/proc` interface regardless of the
kernel's timer frequency, and glibc answers `_SC_CLK_TCK` with 100 on Linux, so
the divisor is the same number on both sides.

**`PartialOrd` is deliberately absent.** The header documents `<`, `<=`, `>` and
`>=` as comparing "clock, user and system time"; the implementation compares
`getCPUTime()` alone. The implementation is what the port reproduces — but an
ordering that ignores wall time while equality does not would be inconsistent
with `PartialEq`, so it is exposed as the explicit `cpu_time_cmp` instead of a
trait. The mismatch between the header's wording and its own code is recorded
here rather than in `OpenMS_CPP_ISSUES.md`, because it is a documentation defect
with no behavioural consequence.

**`reset` cannot fail.** The source routes it through `start()`, which throws
when running; the preceding `clear()` is the only reason that never happens.
The port restarts directly, so there is no failure path to expose.

**`format_seconds` range-checks.** `(TimeType)time_in_seconds` is undefined
behaviour in C++ for NaN, infinity and any value outside the integer type, so
the port returns `Error::InvalidValue` for those. The day count keeps `i64`
width where the source narrows it to `int`. A component the platform cannot
report renders as `n/a` in `summary()`, a form the source cannot produce.

**Arithmetic saturates.** All `TimeSample` arithmetic saturates at the `i64`
bounds instead of wrapping. Reaching them needs about 292 billion years of
accumulated time, so it cannot occur; the port simply does not depend on that to
stay panic-free.

**Serial.** `StopWatch.cpp` carries no `#pragma omp` and neither does this
module, so there is no parallelism gap to record beyond the crate-wide one.

## Checked boundaries and evidence

| Boundary | Value | Source behaviour |
|---|---|---|
| `format_seconds` domain | finite, magnitude ≤ 9e18 s | unchecked cast (UB) |
| `/proc/self/stat` read | ≤ 64 KiB | n/a (uses `times`) |
| Sample arithmetic | saturating `i64` | wrapping / signed overflow |

All 21 class-test sections are mapped in `tests/stop_watch.rs`. The seven
`toString` literals are transcribed verbatim (tier 3) and re-derived from the
source expression at every unit boundary (tier 4). Fourteen sections are
`NOT_TESTABLE` or deferred upstream; each is still covered by the behaviour the
source test says is "tested below". The two ordering sections the source calls
untestable — because it does not control host CPU scheduling — are covered by
the one property that does not depend on it: the relation reads CPU time alone,
so two cleared watches compare `Equal` however they got there.

**The CPU clock is pinned by test, not taken on trust.** The source's `wait()`
is a **busy loop** (`StopWatch_test.cpp:21-29`), not a sleep, and its `bool
stop()` section then bounds four scale-sensitive quantities against the length
of that wait with `t_wait = 0.2 s`: `getCPUTime() > t_wait/2` (line 144),
`getUserTime() > t_wait/2` (line 152), `getUserTime() < t_wait*2` (line 154) and
`getSystemTime() < t_wait*2` (line 156). All four are transcribed with the
source's own half-to-double factors — deliberately generous, so a loaded CI host
does not turn them into flakes — and they are the only coverage that the tick
scale of `MICROS_PER_TICK` is right: a reading divided by the wrong `USER_HZ`
misses them by that whole factor. They hold only because the port's `wait`
helper busy-loops as the source's does; a sleeping thread accrues no CPU time,
so substituting a sleep would make every one of them vacuous. Because CPU time
is a property of the whole process and `cargo test` runs a binary's tests as
threads of one process, every test in the file that burns CPU on purpose takes a
process-wide lock, exactly as `tests/sys_info.rs` does for its working-set
readings.

The remaining timing assertions of that section are transcribed too:
`getClockTime()` above `0.1 s` and above `t_wait*0.95` and below `t_wait*3`, the
component-wise `<=` against the watch that never stopped (lines 164-165), and
`s_resume.getCPUTime() > (t_wait + t_wait_more)/2` across a stop/resume boundary
(line 175). One is not: `getClockTime() < 0.3` (line 126). Line 146 bounds the
same value at `t_wait*3 = 0.6` with the comment "be a bit more loose if e.g. a VM
is busy", so upstream loosened the wall-clock ceiling and left the tighter one
standing, which makes the loosening ineffective. The effective bound is the
looser of the two, and that is the one asserted here; the inconsistency is
recorded as a C++ defect candidate in the provenance record.

The frozen-while-stopped invariant is exact rather than statistical: a stopped
watch computes from stored integers, so repeated reads are bit-identical.

No C++ execution is claimed. Source hashes, line anchors and the class-test
review are in [the provenance record](../tests/data/stop_watch_provenance.json).
