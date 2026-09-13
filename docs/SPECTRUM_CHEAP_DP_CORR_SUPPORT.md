# SpectrumCheapDPCorr support

Port of `src/openms/include/OpenMS/COMPARISON/SpectrumCheapDPCorr.h` and
`src/openms/source/COMPARISON/SpectrumCheapDPCorr.cpp` at Core SDK revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `02ed71fa5cc9080fb6e2651cc1b097022f6047211d1df10165da879f64d06607`,
`.cpp` sha256 `5bad4274f5b7bc69cbc9e8696535f4e718529d987367b9695760a839ea2bff4f`).

Rust: [`comparison::SpectrumCheapDPCorr`](../src/comparison.rs).
Tests: [`tests/comparison_scorers.rs`](../tests/comparison_scorers.rs).
Provenance: [`tests/data/comparison_scorers_provenance.json`](../tests/data/comparison_scorers_provenance.json).

## API mapping

Every public and private member of the header appears here.

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class SpectrumCheapDPCorr : public PeakSpectrumCompareFunctor` | `pub struct SpectrumCheapDPCorr` implementing [`PeakSpectrumCompareFunctor`](PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md) | composition of a `DefaultParamHandler` instead of inheritance |
| `SpectrumCheapDPCorr()` | `SpectrumCheapDPCorr::new() -> Result<Self>` | reproduces `setName("SpectrumCheapDPCorr")`, the three `defaults_.setValue` calls, `factor_ = 0.5` and `defaultsToParam_()` (`SpectrumCheapDPCorr.cpp:29-34`), in that order. No `Default` impl: construction is fallible, as the wave-A functors are |
| `SpectrumCheapDPCorr(const SpectrumCheapDPCorr& source)` | `Clone` | the C++ copy constructor copies the base, `lastconsensus_` and `factor_`, and **not** `peak_map_`, so a copy's `getPeakMap()` is empty while its `lastconsensus()` is not. `Clone` copies every field of the Rust struct, including the peak map, which is a strictly better-defined copy of state the source leaves inconsistent |
| `~SpectrumCheapDPCorr() override` | drop glue | C++ destructor is `= default` |
| `SpectrumCheapDPCorr& operator=(const SpectrumCheapDPCorr& source)` | `Clone::clone_from`, or assignment of a clone | same field set as the copy constructor, plus the self-assignment guard |
| `double operator()(const PeakSpectrum& a, const PeakSpectrum& b) const override` | `PeakSpectrumCompareFunctor::score` **and** `compare(&mut self, ..)` | the source's `const` method writes three `mutable` members; see "The `mutable` split" below |
| `double operator()(const PeakSpectrum& a) const override` | `PeakSpectrumCompareFunctor::self_score`, the trait default | the C++ override is `return operator()(spec, spec)` (`SpectrumCheapDPCorr.cpp:69-72`) |
| `const PeakSpectrum& lastconsensus() const` | `last_consensus(&self) -> &MSSpectrum` | recorded by `compare`, not by `score` |
| `std::map<UInt, UInt> getPeakMap() const` | `peak_map(&self) -> &BTreeMap<usize, usize>` | returns a borrow rather than a copy of the map, and `usize` rather than `UInt`; ascending key order is preserved |
| `void setFactor(double f)` | `set_factor(&mut self, f64) -> Result<()>` | `Exception::OutOfRange` becomes `Error::InvalidRange`; both bounds stay exclusive |
| `private: double dynprog_(const PeakSpectrum&, const PeakSpectrum&, int, int, int, int) const` | private `CheapDpRun::dynamic_program` | takes slices and a cell budget; emits its consensus peaks and peak-map entries into the run instead of into `mutable` members |
| `private: double comparepeaks_(double, double, double, double) const` | private `CheapDpRun::compare_peaks` | returns `Result<f64>`; an `int_cnt` outside `0..=3` is an error rather than `-1` |
| `private: static const std::string info_` | not ported | declared and never defined anywhere in the SDK; it has no observable effect |
| `private: mutable PeakSpectrum lastconsensus_` | field behind `last_consensus()` | see the split below |
| `private: bool keeppeaks_` | not ported as a member | **uninitialised upstream**; see "The shadowed flag" |
| `private: mutable double factor_` | field behind `factor()` / `set_factor()` | see the split below |
| `private: mutable std::map<UInt, UInt> peak_map_` | field behind `peak_map()` | see the split below |

`@htmlinclude OpenMS_SpectrumCheapDPCorr.parameters` documents three entries,
registered with the source's own defaults and description strings:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `variation` | `0.001` (float) | half-window as a fraction of the mean m/z of the pair, and the standard deviation of the Gaussian match term |
| `int_cnt` | `0` (integer) | `0` product, `1` sqrt of product, `2` sum, `3` agreeing intensity |
| `keeppeaks` | `0` (integer) | keep peaks without an alignment partner in the consensus |

## Preserved source conventions

- **The scan.** `operator()` walks both peak lists at once. A pair further apart
  than `variation = (mz_x + mz_y) / 2 * var` cannot be aligned, so the lower peak
  is consumed alone; otherwise the maximal run of mutually pairable peaks is
  measured and, when it is longer than one peak on both sides, handed to
  `dynprog_`. The run loop's four comparisons, its increment order and its two
  early breaks are transcribed one to one; only `!(a < b)` is written as
  `a >= b`, which is the same predicate for the finite coordinates this port
  guarantees.
- **The recurrence.** `dynprog_` is a global alignment with zero gap cost and a
  match term of `comparepeaks_`. The source's odd two-step comparison -
  `max(left, diagonal) > above`, then a second strict `diagonal > left` - fixes
  the tie order as *above*, then *left*, then *diagonal*, and that order is
  reproduced because it selects which consensus peaks and which peak-map entries
  a traceback emits.
- **Boost's normal density.** `comparepeaks_` builds
  `boost::math::normal_distribution<double>(0., variation)` and evaluates its pdf
  at `posa - posb` (`SpectrumCheapDPCorr.cpp:320-327`). The port transcribes
  Boost's pdf statement by statement - `exponent = x - mean`,
  `exponent *= -exponent`, `exponent /= 2 * sd * sd`, `result = exp(exponent)`,
  `result /= sd * root_two_pi` - and takes `root_two_pi` from Boost's own decimal
  literal rather than computing `(2 pi).sqrt()`, because the last bit of that
  constant enters every matched pair.
- **The scan's addition order.** `score` accumulates `dynprog_` results and
  `comparepeaks_` results in encounter order, which is the order reproduced here.
- **The consensus precursor.** One precursor, its m/z the mean of the two inputs'
  first precursors and its charge the **first** spectrum's
  (`SpectrumCheapDPCorr.cpp:94-96`).
- **The two consensus weightings disagree with each other, and both are kept.**
  The one-to-one branch weights the first spectrum by `1 - factor_`
  (`SpectrumCheapDPCorr.cpp:177-178`); the `dynprog_` traceback weights the
  *second* one by `1 - factor_` (`SpectrumCheapDPCorr.cpp:256-257`). Nothing in
  the source reconciles them.
- **`factor_` does not reach the score.** It appears only in consensus peaks, so
  the number `compare` returns is independent of it - asserted in the tests.

## Native differences

### The `mutable` split

`operator()` is `const` and writes `lastconsensus_`, `peak_map_` and `factor_`
through `mutable`. `PeakSpectrumCompareFunctor::score` takes `&self` and is
genuinely read-only here, so the two uses are separated:

* `score(&self, a, b)` returns the number and discards the rest.
* `compare(&mut self, x, y)` returns the same number, records the consensus and
  the peak map, and resets the factor to `0.5` as the source's last statement
  does.

Both call one private routine, so they cannot disagree; a test asserts their
equality on the fixture pair. A failure inside `compare` leaves the previously
recorded consensus, peak map and factor untouched, because the new values are
built in a temporary and committed only on success.

### The shadowed flag

`operator()` declares `bool keeppeaks_ = (int)param_.getValue("keeppeaks")`
(`SpectrumCheapDPCorr.cpp:82`) - a **local** that shadows the member of the same
name. The member is never initialised by either constructor, never assigned
anywhere, and is read by `dynprog_` at lines 275 and 286 to decide whether an
unpaired peak enters the consensus. Reading it is undefined behaviour, and the
consensus a `dynprog_` block contributes is therefore not defined by the source.
This port reads the registered `keeppeaks` parameter in both places, which is the
only defined reading and is plainly what was meant. The class test's
`lastconsensus().size() == 121` does not distinguish the two, because a
self-alignment leaves no unpaired peak.

### Other differences

- `int_cnt` outside `0..=3` is `Err(Error::InvalidValue)`. The source returns
  `-1` behind a `// TODO exception`, which a caller summing scores cannot tell
  from a real contribution.
- A `variation` of zero, or a pair whose mean m/z is zero, makes the Gaussian
  scale zero. Boost raises a domain error there under its default policy; this is
  `Err(Error::InvalidValue)`.
- Unsorted peaks are `Err(Error::UnsortedData)`. The scan assumes ascending m/z
  and mis-aligns silently otherwise; the source checks nothing.
- Negative intensities are refused, because `int_cnt == 1` takes the square root
  of the intensity product.
- `MAX_DP_CORR_CELLS = 1_000_000` bounds the dynamic-programming cells of one
  comparison, summed over every block. The source has no ceiling: a `variation`
  near its own documented maximum of `1` makes a single run span both spectra,
  and `dynprog_` then allocates `(n + 1) * (m + 1)` doubles *and* as many ints.
- A consensus intensity that does not fit `f32` is an error; the source's
  narrowing assignment is undefined behaviour there.
- The unreachable `else` branches of both peak-map writes are reproduced with
  their source semantics and marked unreachable in the code. The scan's version
  compares the two *indices* rather than the stored value, unlike the otherwise
  identical `dynprog_` code; neither can run, because the map is cleared per call
  and each key is written at most once.
- The traceback refuses an unset direction cell instead of looping forever. No
  input can reach it - every cell with `i >= 1` and `j >= 1` is written - but the
  source's `for(;;)` has no other exit.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| unsorted peaks | `Err(Error::UnsortedData)`; source mis-aligns silently |
| negative intensity | `Err(Error::InvalidValue)` |
| two empty spectra | `Ok(0.0)`; the scan never runs, as upstream |
| a pair at m/z zero | `Err(Error::InvalidValue)`; Boost's `check_scale` domain error |
| `variation <= 0` | `Err(Error::InvalidValue)` |
| `int_cnt` outside `0..=3` | `Err(Error::InvalidValue)`; source returns `-1` |
| `keeppeaks` inside `dynprog_` | the registered parameter; source reads an uninitialised member |
| dynamic-programming cells above `MAX_DP_CORR_CELLS` | `Err(Error::InvalidValue)` before either buffer is allocated; recorded state untouched |
| `set_factor` outside `(0, 1)` | `Err(Error::InvalidRange)`, both bounds exclusive, as upstream |
| consensus intensity beyond `f32` | `Err(Error::InvalidValue)` |

OpenMP: no `#pragma omp` in this header or its `.cpp`; the port is serial, as the
source is. The recurrence is sequential in both directions, so there is nothing
here a deterministic parallel reduction could take.

### Class-test sections

`src/tests/class_tests/openms/source/SpectrumCheapDPCorr_test.cpp`
(sha256 `a45c2b84a8437604f94e1b690f1d33ebb9718b021f501290b47601811a83dde0`),
**nine sections, all ported**. The fixtures are `Transformers_tests.dta` (121
peaks) and `Transformers_tests_2.dta` (93 peaks), retained byte-identical.

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `SpectrumCheapDPCorr()` | `cheap_dp_corr_construction_copy_and_assignment` | `name() == "SpectrumCheapDPCorr"`, the three defaults with their source values and types, `factor() == 0.5` | 4 |
| `~SpectrumCheapDPCorr()` | same | construction and drop | 4 |
| `SpectrumCheapDPCorr(const SpectrumCheapDPCorr&)` | same | equal parameters and name, the upstream pair | 3 |
| `SpectrumCheapDPCorr& operator=(const SpectrumCheapDPCorr&)` | same | equal parameters and name | 3 |
| `double operator()(const PeakSpectrum&, const PeakSpectrum&) const` | `cheap_dp_corr_reproduces_the_upstream_golden_scores` | `10145.4` and `12295.5` at the upstream `TOLERANCE_ABSOLUTE(0.1)`, and `10145.449278148695` / `12295.522100159595` to 1e-6; repeated on a second, fresh object as upstream does | 3 |
| `const PeakSpectrum& lastconsensus() const` | same | `121`, equal to the input peak count, plus the consensus precursor's m/z and charge | 3, plus 4 |
| `(Map<UInt, UInt> getPeakMap() const)` | same | `121` entries, every one of them `i -> i` | 3, plus 4 |
| `double operator()(const PeakSpectrum&) const` | same | `self_score(a) == 12295.5`, and `score` equal to `compare` exactly | 3 |
| `void setFactor(double f)` | `cheap_dp_corr_factor_is_range_checked_and_reset_after_a_comparison` | `0.3` accepted, `1.1`, `1.0` and `0.0` refused; the reset to `0.5` after a comparison; the factor's effect on the consensus and its absence from the score | 3, plus 4 |

The transcribed literals `10145.4`, `12295.5` and `121` are tier 3. All four were
reproduced to full precision, before any Rust existed, by an independent Python
model of the scan, of `dynprog_` and of Boost's pdf - `10145.449278148695` and
`12295.522100159595` - which is what makes the tight tolerances above
defensible. `121` is additionally derived: the fixture has 121 peaks, a
self-alignment pairs all of them, and with `keeppeaks` cleared only paired peaks
enter the consensus, so both the consensus length and the peak-map size must
equal the input length. The Rust test asserts that identity, not just the number.

Three further tests cover behaviour the class test never exercises:
`cheap_dp_corr_intensity_terms_and_kept_peaks` pins all four `int_cnt` branches
against closed forms and the `keeppeaks` consensus,
`cheap_dp_corr_refuses_undefined_input` pins the guards, and
`cheap_dp_corr_bounds_its_dynamic_programming_block` pins the cell ceiling and
the untouched recorded state after a refusal.
