# What a third party can recompute

The port's evidence rests on sha256 digests. Manifests pin the C++ source a
port was read from, the fixtures a test compares against and the oracle output
they came from, and prose repeats those digests. A digest is evidence only to a
reader who can recompute it. This page says how many of the recorded digests a
third party who clones this repository can recompute, how many CI recomputes,
and what it would take to raise either number.

Measured on 2026-09-27 at lane A4's commit `b874a5f` (off `main` `f0990c0`) with
`tools/measure_reproducibility.py`, on two machines:

- **the lab**: macOS arm64, this checkout, with `../oracle/` and `.reference/`
  beside it;
- **a third party**: spock (Linux x86_64, Python 3.12.3, git 2.43.0), with
  nothing but `git clone` of five public GitHub repositories (below) and this
  branch. No `../oracle/`, no `.reference/`, no `/ceph`.

The two machines agree on every count below, and two runs on the third-party
machine wrote byte-identical results. The one difference is expected: on the lab
machine all 1,430 digests that name an oracle artifact recompute against
`../oracle/`, which the third party does not have.

## The numbers

The tracked files record **16,860 digest occurrences**, **10,074 distinct**
(artifact, digest) pairs; 16,689 of the occurrences sit in JSON manifests, the
rest in prose, TSV sums files and `Cargo.lock`. An occurrence is one digest
written in one place.

| Tier | Occurrences | Share | Distinct | Share |
| --- | ---: | ---: | ---: | ---: |
| **Bound**: the text names bytes a third party can get, and they hash to it | 14,063 | 83.4% | 7,960 | 79.0% |
| **Content only**: the text does not bind it (or binds it to the wrong file), but a public file hashes to it | 544 | 3.2% | 503 | 5.0% |
| **Lab only**: an `../oracle/` artifact or a build, no public bytes | 1,445 | 8.6% | 814 | 8.1% |
| **None**: nothing anyone can recompute | 808 | 4.8% | 797 | 7.9% |
| **Recomputed by CI** on every push | 1,388 | 8.2% | 1,006 | 10.0% |
| CI, if it also ran `check_core_sdk.py --source` | 8,635 | 51.2% | 3,267 | 32.4% |

So **a third party can recompute 86.6% of the recorded digests** (83.4% bound
plus 3.2% by content), **CI recomputes 8.2%**, and one added CI step would take
that to 51.2%.

What the bound tier is made of:

| Family | Occurrences | Distinct | Recomputed from |
| --- | ---: | ---: | --- |
| C++ package file, digest matches | 8,066 | 2,867 | the public package repository, at the revision the manifest names or the declared pin |
| Tracked file, digest matches now | 5,509 | 4,735 | the clone |
| Tracked file, digest matches an earlier version of it | 345 | 215 | the clone's history (`git log -- <path>`) |
| `Cargo.lock` checksum | 143 | 143 | crates.io, by cargo on every fetch |

The empty-input digest deserves a word: 514 occurrences are the sha256 of zero
bytes. 472 of them name a tracked empty file (an expected-empty stderr, say) and
are counted as bound; the measurement never credits it by content or to CI,
because every empty file and every hash of nothing produces it.

## What CI recomputes, and what it only checks for shape

**No Rust test recomputes a digest.** The crate has no sha256 implementation
and no dependency that provides one; every digest literal in Rust source is
documentation.

**Cargo** verifies each crate's `Cargo.lock` checksum when it downloads it
(`--locked` everywhere in `.github/workflows/rust.yml`): 143 occurrences.

**The Python gates** are the rest. Of the 23 Python steps in `rust.yml`, 11 can
compute a digest at all (the gate or a `tools/` module it imports uses
`hashlib`); run in place with every sha256 they compute logged, they recompute
1,245 recorded occurrences. Most of that is `tools/check_core_sdk.py`, and the
way it is run matters:

`rust.yml` runs `python3 tools/check_core_sdk.py` with no `--source`. Over the
SDK record (`docs/core-sdk-update.json`, `SOURCE_PROVENANCE.json` and the 169
current-SDK reference manifests) that gate holds 4,832 digest assertions, 3,056
distinct (path, digest) pairs:

| Assertions | Count | Without `--source` (CI) | With `--source` |
| --- | ---: | --- | --- |
| In-repo fixtures | 758 | recomputed | recomputed |
| C++ scientific-file inventory | 1,578 | shape only | recomputed |
| Current-source references in manifests | 1,813 | shape only | recomputed |
| Identification-graph sources | 21 | shape only | recomputed |
| Carried-forward references | 220 | shape only | recomputed |
| Registration evidence | 125 | shape only | recomputed |
| `../oracle/` artifacts | 310 | shape, and no copy in the repo | same |
| Original bytes of changed carried references | 7 | shape only | shape only |

Every one of the 3,757 C++ assertions matches core `bc9cc12` byte for byte. On
spock, against a fresh `git clone` of the public core repository checked out at
the pin, `check_core_sdk.py --source` reports *Verified 2107 distinct current
source/registration/reference files* and passes.

## A correction to the earlier count

An earlier review of this codebase (in the project's planning notes, which are
not part of the repository) counted 3,676 distinct assertions and found 833 of
them, 22.7%, checkable by someone who clones the repository, filing the 2,300
C++ source assertions under "needs a 2.6 GB pinned checkout" of the gitignored
`.reference/`. That last step is wrong. `.reference/` is only the lab's
unpacked copy of repositories that are public on GitHub at every revision this
repository names:

- `github.com/okohlbacher/OpenMS4-core` at the pin `bc9cc12` and at the older
  revisions the issue log and the manifests name (`82ce5b3`, `54a232f`,
  `4fdec46`, `7c029e8`, `1b38a99`, `56b181f`);
- `OpenMS4-cli` at `c19e494`, `OpenMS4-topp` at `174b576`, `OpenMS4-test-data`
  at `0cb15f2` (and at `4a2d15e`, which `tests/data/msstats_provenance.json`
  pins);
- this repository, `github.com/okohlbacher/openms-rs`.

What *is* true is that CI does not check them: without `--source` the C++
digests are held to their shape. Like for like, over the 4,832 assertions of
the SDK record, a third party can recompute 4,515 (93.4%: the fixtures and every
C++ assertion) and CI recomputes 758 (15.7%). Over every digest in every tracked
file, the scope of the tables above, the shares are 86.6% and 8.2%.

## How to recompute it yourself

```sh
git clone https://github.com/okohlbacher/openms-rs.git
mkdir -p OpenMS4-tests/packages && cd OpenMS4-tests/packages
git clone https://github.com/okohlbacher/OpenMS4-core.git core
git clone https://github.com/okohlbacher/OpenMS4-cli.git cli
git clone https://github.com/okohlbacher/OpenMS4-topp.git topp
git clone https://github.com/okohlbacher/OpenMS4-test-data.git test-data
git -C core checkout --detach bc9cc12514c768385ce121d6ca4bb710fe1983c4
cd ../../openms-rs
python3 tools/check_core_sdk.py --source ../OpenMS4-tests/packages/core
python3 tools/check_source_citations.py --report
python3 tools/measure_reproducibility.py --ci --source ../OpenMS4-tests/packages/core
```

The package repositories go beside the clone, at
`../OpenMS4-tests/packages/<name>`; the tools read the other pins out of their
git objects, never out of a working tree, so only `core` has to be checked out
for `--source`, which requires a clean checkout at exactly the pin. The
generator gates call `rustfmt`, so the Rust toolchain has to be on `PATH`
(without it `generate_enzymes.py --check` fails on a missing `rustfmt`, not on
a digest). Nothing here needs `../oracle/`, `.reference/` or a C++ build.
`measure_reproducibility.py --ci` runs only the gates that can hash; it never
runs `tools/check_schema_feature_graph.py`, which builds the crate with cargo.

**Citations.** The same clones let a third party run the source-citation
checker, which reads core revisions out of the core package's git objects when
no `.reference/` checkout holds them. On spock it checks 4,112 C++ line
citations and passes; the lab checks 4,139. The difference is exactly the
third-party headers the lab keeps under `../oracle/` (Eigen 5.0.1, Boost 1.89.0
and GCC 14.4 libstdc++, admitted by digest): 27 citations the lab reads
against them are skipped or left unresolved on a machine that does not have
them (47 skipped and 17 unresolved there, against 34 and 3 here).

## What nobody outside the lab can recompute

**Lab only, 1,445 occurrences (8.6%).** 1,365 name an artifact under
`../oracle/` whose bytes are public nowhere: 751 files in 76 oracle
directories, 134 MB in all, the largest a 48 MB Boost.Regex corpus.
`../oracle/` sits outside the repository by policy (`check_core_sdk.py` refuses
a copy inside it). All of them recompute on the lab machine. (65 more oracle
digests equal a tracked, historical or pinned file and count as content only.)
The other 80 name a build: a driver binary, `libOpenMS`, an install under
`/ceph`, the macOS product SDK. Only an identical build reproduces those.

**None, 808 occurrences (4.8%).**

- **Console logs never kept, 552.** Validation records store `log_sha256` (and
  `stderr_sha256`) for each gate or probe they ran, but not the log (566 in all;
  14 of them happen to equal a tracked file).
- **In-repo digests of states no commit holds, 12.** They name a tracked file,
  but match neither the file now, nor any version of it in the history of
  `HEAD`, nor any public file: the validation ran on an uncommitted working tree.
  Five are in `tests/data/mzml_schema_provenance.json` (`native_baseline`
  names head `1c13a69`, where `src/format/mzml.rs` hashes differently), two
  each in `tests/data/mzml_precursor_activation_validation/summary.json` and
  `tests/data/mzml_sqlite_handler_provenance.json`, one each in
  `docs/proforma-spectra-validation.json`,
  `tests/data/progress_logger_provenance.json` and
  `tests/data/sdk_numeric_formatting_provenance.json`.
- **Unbound, 239.** A digest whose text names no artifact that could be found:
  a `snapshot_sha256`, `transformed_sha256` or `driver_sha256` with no path, a
  table cell, an abbreviated header path that no package holds. They may well
  be right; nobody can recompute them from what is written.
- **Upstream, not found, 5.** Header paths under `include/` that name an
  installed SDK rather than a source tree (`tests/data/mzml_sqlite_swath_provenance.json`).

**Content only, 544 occurrences (3.2%).** The text does not bind the digest to
the file it is of - a prefixed key with no path, prose naming a neighbour - but
a tracked file, a version in the history, or a pinned package file hashes to
it. A reader can recompute these by searching; the text does not tell them
where.

## What it would take to raise it

In order of gain for the effort:

1. **Run `check_core_sdk.py --source` in CI** against a clone of the public core
   at the pin: CI's share goes from 8.2% to 51.2% of the occurrences (10.0% to
   32.4% distinct). It costs one step and a 500 MB checkout; a shallow fetch of
   the one commit is enough, since `--source` needs only a clean working tree at
   the pin. (`rust.yml` is the lead's file.)
2. **Run the citation checker in CI** the same way, with the four packages
   cloned beside the checkout. CI runs only `test_source_citations.py` today;
   the checker itself needs the pins.
3. **Publish the oracle artifacts** the manifests name - 751 files, 134 MB - as
   an archive or a separate public repository with its own digests, outside this
   repository (C++ drivers are among them, and C++ never enters this one).
   That moves 1,365 occurrences, 8.1 points, from lab only to bound. The 80
   build digests stay: a build is not an artifact anyone can re-derive byte for
   byte.
4. **Keep the logs, or stop recording their digests**: 566 occurrences, 3.4
   points, are digests of text that was thrown away.
5. **Re-record the 12 uncommitted-state digests** at a commit, or mark them as
   working-tree snapshots.
6. **Bind every digest to a path**: the 544 content-only and 239 unbound
   occurrences. A manifest key that names what it hashes (`source_path` beside
   `source_sha256`) is enough.

With 1 to 6 done, what remains is the builds: about 0.5% of the occurrences.

## Measuring it again

```sh
python3 tools/measure_reproducibility.py                  # the tiers and families
python3 tools/measure_reproducibility.py --ci --json      # with what CI recomputes, as JSON
python3 tools/measure_reproducibility.py --list none      # every digest nobody can recompute
python3 tools/test_measure_reproducibility.py             # the binding rules it relies on
```

It measures and never fails a build. It takes about two minutes, most of it
hashing the pinned packages (1 GB) and every blob in the history of `HEAD` to
answer the content lookups.
