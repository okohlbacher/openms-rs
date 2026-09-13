# Process and system memory reporting

`system::sys_info` covers `SYSTEM/SysInfo.h` and `SYSTEM/SysInfo.cpp` at SDK
`bc9cc12`. Its direct TOPP consumers are `FileInfo` and `TICCalculator`.

This header is the most platform-specific in the domain: the source implements
every reading three times, once each for Windows, macOS and Linux, and its own
documentation warns that outside Windows the numbers "might be very unreliable,
depending on operating system and kernel version". **This crate is Linux-first**,
and the port says so in its types rather than in a comment.

## API mapping

Every public member of the header, with its Rust counterpart.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `std::string bytesToHumanReadable(UInt64 bytes)` | `sys_info::bytes_to_human_readable(u64) -> String` | Portable; exact on every platform. |
| `static bool SysInfo::getProcessMemoryConsumption(size_t&)` | `sys_info::process_memory_consumption() -> Option<u64>` | KiB. `None` replaces the `false` + zeroed out-parameter. |
| `static bool SysInfo::getProcessPeakMemoryConsumption(size_t&)` | `sys_info::process_peak_memory_consumption() -> Option<u64>` | KiB. |
| `static bool SysInfo::getFreeSystemMemory(size_t&)` | `sys_info::free_system_memory() -> Option<u64>` | KiB. |
| `static Int64 SysInfo::getProcessId()` | `sys_info::process_id() -> u32` | `std::process::id`; portable. |
| `SysInfo::MemUsage` | `sys_info::MemUsage` | |
| `MemUsage::MemUsage()` | `MemUsage::new`, `MemUsage::default` | Both record the first time point, as the source's constructor does. |
| `MemUsage::mem_before` | `MemUsage::mem_before: Option<u64>` | |
| `MemUsage::mem_before_peak` | `MemUsage::mem_before_peak: Option<u64>` | |
| `MemUsage::mem_after` | `MemUsage::mem_after: Option<u64>` | |
| `MemUsage::mem_after_peak` | `MemUsage::mem_after_peak: Option<u64>` | |
| `void MemUsage::reset()` | `MemUsage::reset(&mut self)` | |
| `void MemUsage::before()` | `MemUsage::before(&mut self)` | |
| `void MemUsage::after()` | `MemUsage::after(&mut self)` | |
| `std::string MemUsage::delta(const std::string& event = "delta")` | `MemUsage::delta(&mut self, event: &str) -> Result<String>` | No default argument in Rust: pass `"delta"` for the source default. |
| `std::string MemUsage::usage()` | `MemUsage::usage(&mut self) -> String` | |
| `std::string MemUsage::diff_str_(size_t, size_t)` (private) | `difference_string` (private) | Emits the sign the source drops; see *Native differences*. |
| `read_off_memory_status_linux`, `read_available_memory_linux`, `statm_t` (file-local) | `status_field`, `meminfo_field`, `kibibytes` (private) | Different `/proc` files, same counters. |

Native additions: `MAX_EVENT_LABEL_BYTES`.

There is no unported public member.

## Preserved source conventions

**Units.** Every memory figure is in kibibytes, which the source calls "KB", and
every report renders a difference in mebibytes by integer division of the
kibibyte value by 1024 — truncating towards zero, so a 1535 KiB change prints as
`1 MB`.

**`bytesToHumanReadable` is reproduced exactly**, including three things that
look like bugs and are not worth changing:

- the unit is never pluralised, so `1 byte`, `2 byte`, `1000 byte`;
- the value is printed with `std::setprecision(4)` in the *general* format, so
  four significant digits with trailing zeros dropped — `2 KiB`, `1.5 KiB`,
  `45.34 MiB`, `1023 byte`. The divide-until-below-1024 loop stops at the first
  quotient, so the printed value is either exactly `0` or in `[1, 1024)`; over
  that domain the general format never selects scientific notation and the
  significant-digit count is fixed by the width of the integer part, which is
  how the port computes it. That equivalence was checked independently against
  the C `%.4g` conversion over 500 027 values covering the whole reachable
  domain, with no disagreement;
- from 2^60 bytes upwards the six-entry unit table is exhausted and the source
  returns its literal apology, `Congrats. That's a lot of bytes: <count>`, which
  is reachable for a `UInt64` and is therefore reproduced verbatim.

**`MemUsage` report assembly.** `delta` and `usage` record the second time point
themselves if it is missing. The peak component is appended only when the peak
reading is greater than zero — the source's own test for whether the platform
supports it — and the wording (`(working set delta)`, `(peak working set delta)`,
`(working set)`, `(peak working set)`) and the `, ` separator are the source's.

## Native differences

**Every reading is an `Option`, and only Linux answers.**

| Reading | Linux | macOS | Windows | other |
|---|---|---|---|---|
| working set | `VmRSS` from `/proc/self/status` | `None` | `None` | `None` |
| peak working set | `VmHWM` from `/proc/self/status` | `None` | `None` | `None` |
| free physical memory | `MemAvailable`, else `MemFree`, from `/proc/meminfo` | `None` | `None` | `None` |
| process id | real | real | real | real |

The source's macOS branch uses Mach `task_info` and `host_statistics64` and its
Windows branch uses `GetProcessMemoryInfo` and `GlobalMemoryStatusEx`; neither is
reachable from this crate without a platform binding it does not take. Returning
`None` is deliberate: the source's contract is "`false`, and the out-parameter is
set to 0", and a caller that ignores the `bool` sees a plausible zero. `None`
cannot be mistaken for a measurement.

**Two of the three Linux readings come from a different `/proc` file than the
source's, on purpose; the third comes from the same one.** The source computes
the working set as `statm.resident * sysconf(_SC_PAGESIZE) / 1024` and the peak
as `getrusage(...).ru_maxrss`. The kernel publishes both of those counters in
`/proc/self/status` as `VmRSS` and `VmHWM`, already in kibibytes — the source's
own quoted `proc(5)` excerpt says `resident` is "the same as VmRSS in
/proc/[pid]/status". Reading them there gives the identical number without a
`sysconf` call, which is what makes the module libc-free.

Free memory is the exception. The source's **primary** path already reads
`MemAvailable` from `/proc/meminfo` and returns on the first match
(`SysInfo.cpp:107-133`), which is exactly what this port does — same file, same
key, no divergence. Only the fallback differs: on a kernel too old to publish
`MemAvailable` the source falls back to
`sysconf(_SC_AVPHYS_PAGES) * sysconf(_SC_PAGESIZE) / 1024` and this port falls
back to `MemFree` in the same file. Both answer the narrower question —
physically free pages, without the reclaimable page cache `MemAvailable` counts
— but the port does not claim the two agree to the byte, because
`_SC_AVPHYS_PAGES` is answered by whichever libc the source was linked against.

**A negative delta keeps its sign.** `diff_str_` appends `"-"` to a local string
and then **assigns over that string** with the magnitude instead of appending to
it, so the sign is discarded and every source report reads as an increase. This
port emits the sign. The magnitude is unchanged: the source takes the absolute
value before dividing, which for truncation towards zero is the same number
either way. Recorded as a C++ defect candidate in
[the provenance record](../tests/data/sys_info_provenance.json); it belongs in
`OpenMS_CPP_ISSUES.md`, which the integrating agent owns.

**`Option` replaces the `0` sentinel.** The source treats `mem_after == 0` as
"not recorded yet", so a platform that reports a genuine zero is resampled on
every print, and "unsupported" and "zero" are the same value. The port keeps
`None` for "not recorded" and for "the platform does not report this", and a
reading that is missing renders as `n/a`, never as a number.

**`delta` is bounded and fallible.** The source takes any `event` string. This
port refuses one above `MAX_EVENT_LABEL_BYTES` (1 MiB) before building the
report, so a caller-controlled string cannot drive an unbounded allocation.
There is no default-argument syntax in Rust, so the source's `event = "delta"`
default becomes an explicit argument.

**`process_id` returns `u32`.** The source returns `Int64` only because `getpid`
and `_getpid` have different signed types; no platform identifier needs the sign
or the width.

**`/proc` parsing is strict.** Reads are capped at 256 KiB, and a matched line
must carry the `kB` unit — which the source does not check. Every key this
module asks for publishes kibibytes, so a line without the suffix is not the
record that was wanted.

**Serial.** `SysInfo.cpp` carries no `#pragma omp`; there is no parallel
behaviour to reproduce.

## Checked boundaries and evidence

| Boundary | Value | Source behaviour |
|---|---|---|
| `/proc` read | ≤ 256 KiB | unbounded `fopen`/`getline` |
| `MemUsage::delta` event label | ≤ 1 MiB | unbounded |
| `/proc` line unit | `kB` required | unchecked |

Both class-test sections are mapped in `tests/sys_info.rs`. All six
`bytesToHumanReadable` literals are transcribed (tier 3) and the unit boundaries
around them — `0`, `1`, `1000`, `1023`, `1024`, `1536`, `45.34 MiB`, `2^60`,
`u64::MAX` — are derived from the source's format (tier 4). The memory section
is reproduced without its 20 MB mzML fixture: a 64 MiB allocation is made
resident page by page and the same "grew by more than 10 000 KB" assertion is
made, which keeps a SYSTEM test free of a format dependency and gives more than
six times the margin the source asks for — 64 MiB is 65 536 KiB against a
10 000 KB threshold. On a platform where the readings are unavailable, the same
section asserts that they are all `None`.

No C++ execution is claimed. Source hashes, line anchors and the class-test
review are in [the provenance record](../tests/data/sys_info_provenance.json).
