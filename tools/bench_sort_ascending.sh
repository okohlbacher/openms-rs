#!/usr/bin/env bash
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
#
# The measurement behind section 8 of docs/BENCHMARKS.md: what
# math::statistic_functions::sort_ascending costs, per implementation and per
# sample size, in wall clock and in peak resident memory.
#
# Three implementations, on the same bytes:
#
#   library      slice::sort_unstable_by(f64::total_cmp) -- what the fast path
#                of lead decision D17 runs; in place, allocates nothing.
#   faithful     source_sort_by(|a, b| a < b) -- the libstdc++ introsort
#                reproduced comparison by comparison, which is what
#                sort_ascending runs when the permutation is observable.
#   entry_point  median(&mut sample) -- the public entry point, which picks the
#                path. On this sample (no NaN, one spelling of zero) it picks
#                the fast one, so this row is the regression the decision was
#                taken about.
#
# Each cell is one fresh process running one implementation at one size, so the
# peak RSS /usr/bin/time reports is that implementation's own and not the
# harness's. The test binary is built once, in release.
#
# --corpus instead reports how often the fast path is actually taken, over
# every mzML fixture in tests/data: one sample per file for the MS1 intensities
# of `FileInfo -s` and the MS1 retention times of `-c`, and one per spectrum for
# the peak m/z of `-c`.
#
# Usage:
#   tools/bench_sort_ascending.sh [--reps N] [--sizes "1000 10000 ..."] [--out FILE]
#   tools/bench_sort_ascending.sh --corpus
#
# Output is a markdown table on stdout and, with --out, the raw per-cell TSV.

set -euo pipefail

REPS=5
SIZES="1000 10000 100000 1000000 10000000"
RAW=""
CORPUS=0

while [ $# -gt 0 ]; do
  case "$1" in
    --reps)  REPS="$2";  shift 2 ;;
    --sizes) SIZES="$2"; shift 2 ;;
    --out)   RAW="$2";   shift 2 ;;
    --corpus) CORPUS=1;  shift 1 ;;
    -h|--help) sed -n '6,38p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

cd "$(dirname "$0")/.."

echo "# building the test binary in release" >&2
cargo test --release --test statistic_functions --no-run >/dev/null 2>&1 || {
  cargo test --release --test statistic_functions --no-run
  exit 1
}

# The most recently built statistic_functions test binary.
BIN=$(cargo test --release --test statistic_functions --no-run --message-format=json 2>/dev/null \
      | python3 -c '
import json, sys
path = None
for line in sys.stdin:
    line = line.strip()
    if not line.startswith("{"):
        continue
    record = json.loads(line)
    if record.get("reason") == "compiler-artifact" and record.get("executable"):
        target = record.get("target", {})
        if target.get("name") == "statistic_functions":
            path = record["executable"]
print(path or "")')

if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
  echo "could not locate the statistic_functions test binary" >&2
  exit 1
fi
echo "# binary: $BIN" >&2

if [ "$CORPUS" = "1" ]; then
  OPENMS_BENCH_IMPL=corpus "$BIN" how_often_the_fast_path_is_taken_over_the_mzml_corpus \
    --exact --nocapture | grep '^CORPUS'
  exit 0
fi

# GNU time reports peak RSS in KiB after "Maximum resident set size"; BSD/macOS
# /usr/bin/time -l reports it in bytes on a "maximum resident set size" line.
case "$(uname -s)" in
  Darwin) TIME_ARGS="-l" ;;
  *)      TIME_ARGS="-v" ;;
esac

TSV=$(mktemp)
trap 'rm -f "$TSV"' EXIT

for n in $SIZES; do
  for impl in library faithful entry_point; do
    err=$(mktemp)
    out=$(OPENMS_BENCH_N="$n" OPENMS_BENCH_IMPL="$impl" OPENMS_BENCH_REPS="$REPS" \
          /usr/bin/time $TIME_ARGS \
          "$BIN" sort_ascending_benchmark --exact --nocapture 2>"$err" || true)
    ns=$(printf '%s\n' "$out" | awk -F'\t' '$1=="BENCH" {print $5}')
    # Peak RSS, normalised to KiB.
    kib=$(awk '
      tolower($0) ~ /maximum resident set size/ {
        for (i = 1; i <= NF; i++) if ($i ~ /^[0-9]+$/) { value = $i }
      }
      END {
        if (value == "") { print "" }
        else if (unit == "bytes") { printf "%d\n", value / 1024 }
        else { print value }
      }' unit="$([ "$TIME_ARGS" = "-l" ] && echo bytes || echo kib)" "$err")
    rm -f "$err"
    if [ -z "$ns" ]; then
      echo "# cell n=$n impl=$impl produced no timing" >&2
      continue
    fi
    printf '%s\t%s\t%s\t%s\n' "$n" "$impl" "$ns" "${kib:-NA}" >>"$TSV"
    echo "# n=$n $impl ${ns}ns ${kib:-NA}KiB" >&2
  done
done

if [ -n "$RAW" ]; then cp "$TSV" "$RAW"; fi

echo
echo "| n | library (fast path) | faithful (libstdc++) | entry point | faithful / library | peak RSS library | peak RSS faithful |"
echo "| --- | --- | --- | --- | --- | --- | --- |"
python3 - "$TSV" <<'PY'
import sys

rows = {}
for line in open(sys.argv[1]):
    n, impl, ns, kib = line.rstrip("\n").split("\t")
    rows.setdefault(int(n), {})[impl] = (int(ns), kib)


def human(ns):
    if ns < 1_000:
        return f"{ns} ns"
    if ns < 1_000_000:
        return f"{ns / 1_000:.0f} us"
    if ns < 1_000_000_000:
        return f"{ns / 1_000_000:.1f} ms"
    return f"{ns / 1_000_000_000:.2f} s"


def mib(kib):
    return "NA" if kib == "NA" else f"{int(kib) / 1024:.0f} MiB"


for n in sorted(rows):
    cell = rows[n]
    lib = cell.get("library")
    fai = cell.get("faithful")
    ent = cell.get("entry_point")
    ratio = f"{fai[0] / lib[0]:.1f}x" if lib and fai and lib[0] else "NA"
    print(
        f"| {n:,} | {human(lib[0]) if lib else 'NA'} "
        f"| {human(fai[0]) if fai else 'NA'} "
        f"| {human(ent[0]) if ent else 'NA'} "
        f"| {ratio} "
        f"| {mib(lib[1]) if lib else 'NA'} "
        f"| {mib(fai[1]) if fai else 'NA'} |"
    )
PY
