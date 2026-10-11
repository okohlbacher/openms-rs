# The DefaultParamHandling trait

`openms::param::DefaultParamHandling` is the typed side of the source's `DefaultParamHandler` base class. [DEFAULT_PARAM_HANDLER_SUPPORT.md](DEFAULT_PARAM_HANDLER_SUPPORT.md) describes the handler object itself (current and default trees, subsections, the two check switches); this trait is what a class that *derives* from it becomes. An algorithm without it has no `Param` surface and cannot be driven by a TOPP tool's `-ini`.

## The contract

| Source | Rust |
|---|---|
| `DefaultParamHandler(name)`, `getName()` | `const NAME` |
| `subsections_` | `const SUBSECTIONS` (default none) |
| `check_defaults_`, `warn_empty_defaults_` | `const CHECK_DEFAULTS`, `const WARN_EMPTY_DEFAULTS` (default `true`) |
| the constructor's `defaults_.setValue(...)`, `getDefaults()` | `fn defaults() -> Result<Param>` |
| `param_`, `getParameters()` | `fn to_param(&self) -> Result<Param>` |
| `updateMembers_()` | `fn update_members(&Param) -> Result<Self>` |
| `setParameters(p)` on a new object | provided `from_param`, `from_param_with_warnings`, `from_param_with_handler` |
| `setParameters(p)` on an existing object | provided `set_parameters(&mut self, &Param)` |
| `defaultsToParam_()` at the end of the constructor | `Default`, whose `to_param()` equals `defaults()` |

An implementation writes four items - `NAME`, `defaults`, `to_param`, `update_members` - and, when the type carries members that are not parameters (resource ceilings, compatibility profiles), `keep_native_members`, so that `set_parameters` keeps them. `update_members` receives a tree that already holds every default and has passed the type and restriction check, so it reads each key without a fallback, exactly as the source's `updateMembers_` does.

The shape is the one `PeakPickerHiRes` established (`src/processing/peak_picking.rs`): an inherent `defaults()`/`to_param()`/`from_param()`/`from_param_with_warnings()` over a private `from_complete_param`. A type written that way adopts the trait with four delegating lines. The rustdoc of the trait carries a compiled example.

Helpers for the source's spellings live beside it in `openms::param::handler`: `set_flag` (a `"true"`/`"false"` string with those two valid strings, the source's boolean), `flag_parameter`, `unsigned_parameter` (the source's `(UInt)param_.getValue(key)`, which throws `ConversionError` for a negative value) and `int_parameter` (a count written back as the source's 32-bit `int`). `Param::replace_value` changes a value while keeping its description and restrictions.

## Types that implement it

| Type | Source class | Parameters | Notes |
|---|---|---|---|
| `processing::NLargest` | `NLargest` | `n` | `peakcount_` is the field `n`; a negative `n` is the source's `ConversionError`. |
| `processing::ThresholdMower` | `ThresholdMower` | `threshold` | The source re-reads `param_` on every call; the field is that value. |
| `processing::window_mower::WindowMower` | `WindowMower` | `windowsize`, `peakcount`, `movetype` | Keeps its native ceilings through `set_parameters`. |
| `processing::Normalizer` | `Normalizer` | `method` | The enum cannot hold a method outside the valid strings. |
| `processing::RankScaler`, `processing::SqrtScaler` | `RankScaler`, `SqrtScaler` | none | Empty defaults: every `setParameters` warns "No default parameters for DefaultParameterHandler '...' specified!", as the Release build does. |
| `processing::smoothing::GaussFilter` | `GaussFilter` | `gaussian_width`, `ppm_tolerance`, `use_ppm_tolerance`, `write_log_messages` | The width the active mode does not use is kept for `to_param`; `write_log_messages` is stored, the warning itself is not emitted. |
| `processing::smoothing::SavitzkyGolayFilter` | `SavitzkyGolayFilter` | `frame_length`, `polynomial_order` | An even frame is made odd as in the source, which keeps the given length in `param_`; `updateMembers_`'s order check is an error. |
| `processing::LinearResamplerAlign` | `LinearResamplerAlign` | `spacing`, `ppm` | `ppm` = `true` is `Error::Unsupported`: the ppm branch is not ported. |
| `processing::iterative::PeakPickerIterative` | `PeakPickerIterative` | nine, four with the source's trailing underscore | `sn_win_len_`/`sn_bin_count_` configure the typed noise estimator. |
| `processing::peak_picking::PeakPickerHiRes` | `PeakPickerHiRes` | seventeen: eight, and the median estimator's nine under the `SignalToNoise:` node | The reference shape; its inherent methods delegate to the trait. |
| `processing::peak_picking::SignalToNoiseEstimatorMedian` | `SignalToNoiseEstimatorMedian` | nine | |
| `processing::mean_noise::SignalToNoiseEstimatorMeanIterative` | `SignalToNoiseEstimatorMeanIterative` | nine | The three range parameters are kept whatever `auto_mode` selects. |
| `chemistry::TheoreticalSpectrumGenerator` | `TheoreticalSpectrumGenerator` | thirty | Intensities are `f64`, as the source's `double` members are; see below. |
| `chemistry::SpectrumAnnotator` | `SpectrumAnnotator` | eight | A negative `topNmatch_fragmenterrors` is the source's `ConversionError`. |
| `chemistry::NucleicAcidSpectrumGenerator` | `NucleicAcidSpectrumGenerator` | twenty-three | No restrictions on any intensity, as in the source. |

`TheoreticalSpectrumGenerator` held its intensities as `f32` before it gained the contract. The source's are `double` (`TheoreticalSpectrumGenerator.h:141-153`), multiplied in double precision and narrowed once into `Peak1D`. With `b_intensity` 0.7 and `relative_loss_intensity` 0.9, which `f32` cannot represent, the `f32` fields put intensities one ulp away from the Release build's; they are now `f64` and match it bit for bit.

## Native differences

- A refused tree changes nothing. The source assigns `param_` before `checkDefaults` can throw, so its `getParameters()` then shows the refused values while the members keep the old ones.
- Warnings are returned, not logged. The unknown-key warning reads `<Name>: unknown parameter '<key>'` where the source logs `Warning: <Name> received the unknown parameter '<key>'!`.
- An `updateMembers_` that throws leaves the source object half updated; `update_members` returns the error and nothing changes.
- `to_param()` rebuilds the tree from the typed members, so it carries the defaults' descriptions and restrictions. The source's `getParameters()` after `setParameters(p)` keeps each entry of `p` as given and any unknown key; `from_param_with_handler` returns a handler holding exactly that tree.
- Values a typed member cannot hold are errors where the source would wrap or truncate: a count above `i32::MAX` in `to_param`, a negative `max_isotope` for the coarse isotope model.

## Evidence

Tier 1, executed on ibminode06 against the Release install `openms4-release-bc9cc12-c19e494-174b576`, every case twice in fresh processes with identical output:

- `../oracle/B6-param-defaults` writes `getDefaults()` of each of the sixteen classes through `ParamXMLFile`, and `getParameters()` after `setParameters` with an unknown key, with an empty tree and with the values of six `setParameters` probes (five refused) and the `NLargest(7)` constructor. `tests/default_param_handling.rs` compares every name, value, type, description, restriction and the order, and the warnings and exception texts.
- `../oracle/B6-tsg-intensity` prints every peak of three `TheoreticalSpectrumGenerator` configurations set through `setParameters` (184 peaks), plus three `setParameters` probes of the chemistry contracts. `tests/theoretical_parameter_intensities.rs` requires bit-identical intensities, names and charges.

`tests/data/param_defaults_provenance.json` records the source, fixture, oracle and binary hashes. `SignalToNoiseEstimatorMeanIterative<MSChromatogram>` does not compile at the pin (the header calls `getMZ()` on `ChromatogramPeak`), so only the `MSSpectrum` instantiation was dumped.

## Open

- The m/z values of the `TheoreticalSpectrumGenerator` cases differ from the Release build's by at most two ulp: the port sums the masses in its own order. The intensities, names and charges are identical.
- `GaussFilter`'s `write_log_messages` warning is not emitted.
- Types outside this list that were written in the `PeakPickerHiRes` shape adopt the trait by adding the four delegating lines.
