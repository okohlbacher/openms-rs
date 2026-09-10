# Spectrum annotation and ion names

The native chemistry module implements all three scientific entrypoints of the
pinned `SpectrumAnnotator` and all four `IonNaming` helpers. It uses the existing
native theoretical spectrum generator, spectrum alignment and identification
records. No C++ library, regular-expression package or additional dependency is
required. See the [source reference review](SPECTRUM_ANNOTATION_REFERENCE_REVIEW.md)
and [fixture provenance](../tests/data/spectrum_annotation_provenance.json).

## API

Import `openms::chemistry::{SpectrumAnnotator, TheoreticalSpectrumGenerator}` and
`openms::comparison::SpectrumAlignment`. Each operation returns `Result<()>` and
commits its affected fields only after all calculations succeed.

| Method on `SpectrumAnnotator` | Result |
| --- | --- |
| `annotate_matches(&mut spectrum, &hit, &generator, &alignment)` | Sort measured peaks by m/z; replace all auxiliary arrays with ion labels, charges and match errors |
| `add_ion_match_statistics(&mut identification, &mut spectrum, &generator, &alignment)` | Calculate enabled statistics for each hit; retain the final hit's annotations on the measured spectrum |
| `add_peak_annotations(&mut hit, &spectrum, &generator, &alignment, include_unmatched_peaks)` | Replace the hit's peak annotations; leave the measured spectrum unchanged |

Set `generator.add_metainfo = true` for spectrum annotations and statistics.
The source's generator defaults disable metadata. If a match needs missing
arrays, `annotate_matches` returns an error. With no matches, it writes empty
labels and zero charges/errors. Peak-annotation output separately follows the
source named-array lookup and permits missing labels/charges as empty/zero.

Run the [example](../examples/annotate_spectrum.rs):

```sh
cargo run --locked --offline --example annotate_spectrum
```

It annotates the source's eleven IFSQVGK measurements, prints complete ion names
and charges, and reports match counts, intensity, longest series and mass errors.
The measured input is kept separate from generated theoretical peaks. This is
an annotation demonstration, not a search engine or a confidence estimate.

## Matching and data preservation

Generation uses fragment charges 1 through `min(hit.charge, 2)`; hit charge must
be positive. Precursor charge is left to the generator's ordinary inference.
All configured ion series, losses, isotopes, internal/immonium fragments and
owned modification chemistry remain available subject to the generator's
[documented constraints](THEORETICAL_SPECTRA.md). Generated peaks must actually
be sorted for alignment; disabling generator sorting is not itself an error
when the output happens to remain sorted.

Spectrum annotation replaces every float, integer and string array, with:

| Name | Stored values |
| --- | --- |
| `IonNames` | String label per measured peak; empty for unmatched peaks |
| `Charges` | Integer fragment charge; zero for unmatched peaks |
| `IonMatchError` | Absolute m/z error narrowed to f32; zero for unmatched peaks |

`IonNames` is the implemented source spelling, despite singular `IonName` in its
header prose. Errors are absolute Daltons even for ppm alignment. The spectrum
also receives string metadata `fragment_mass_tolerance` and
`fragment_mass_tolerance_ppm` (`"0"`/`"1"`, the source integer projected into
the native spectrum's string metadata). Its other metadata, acquisition
fields and attached identifications are retained.

Absolute alignment is one-to-one. Ppm alignment may assign several theoretical
peaks to one measured peak; the last alignment wins in spectrum arrays and in
the include-unmatched peak-annotation branch. Matched-only peak annotations
instead preserve one record per alignment, including repeated measured peaks.
Every annotation uses measured m/z and measured f32 intensity widened to f64.
This follows [native alignment conventions](COMPARISON_SUPPORT.md).

Statistics first sort by increasing intensity. If precursor statistics are
enabled and at least one precursor exists, they sort by m/z again. The final
hit determines the resulting peak order and arrays. Native sorts preserve ties,
whereas the C++ source does not specify equal-intensity order. Empty measured
spectra or identifications without hits return immediately without changes.

Validation covers the data actually used. Unrelated metadata, existing hit
annotations and evidence graphs are not deep-cloned or recursively validated;
replaced spectrum arrays need not be copied. Any later error leaves both mutable
inputs unchanged. File writers still apply their own full validation.

## Statistics options and conventions

All boolean options default to true; `top_n_fragment_errors` defaults to seven.

| Option | Metadata or behavior |
| --- | --- |
| `basic_statistics` | `matched_ions`, `matched_intensity`, `matched_ion_number`, `peak_number`, `sum_intensity` |
| `list_of_ions_matched` | Source option retained as a no-op; `basic_statistics` controls the ion list |
| `max_series` | `max_series_type`, `max_series_size` for consecutive a/b/c/x/y/z ordinals |
| `sn_statistics` | `sn_by_matched_intensity`, `sn_by_median_intensity` |
| `precursor_statistics` | Integer `precursor_in_ms2` (0/1); uses the raw alignment tolerance as absolute Daltons even in ppm mode |
| `top_n_fragment_errors` | Nonzero enables all five error statistics; zero disables that block |
| `fragment_error_statistics` | Source option retained as a no-op; the top-N value controls error statistics |
| `terminal_series_match_ratio` | `NTermIonCurrentRatio`, `CTermIonCurrentRatio` |

Identification metadata receives numeric `fragment_match_tolerance`, whose key
differs from the spectrum key. Existing unrelated hit/identification metadata
and hit peak annotations survive statistics updates.

Top-N errors are selected in descending measured-intensity order, then padded
with zeros to exactly N when fewer matches exist. Their mean, mean square and
sample standard deviation use that padded length. Median uses the whole error
list; IQR uses sorted indices `n/4` and `n/4 + n/2`. The native implementation
uses those intended indices safely for n=1–3, where the source's chained
selection subranges are invalid. Nonempty N=1 is an error because sample
deviation is undefined; an empty error list retains the source's five zeros.

Terminal-current and series-length classification follows the source's complete
label grammars, including loss-label distinctions. It does not use the broader
ion-name parser below. Longest-series ties retain lexical a/b/c/x/y/z order.
Source intensity/error accumulation occurs before rejecting an out-of-range
series ordinal from the ion list; list counts can therefore differ from other
matched contributions.

The two S/N statistics retain source f32/f64 arithmetic and special zero results
for all-matched or missing median groups. Enabled final ratios/statistics must
be finite. Undefined cases, such as zero matched intensity with terminal ratios
enabled, return an error; callers can disable the corresponding statistic.
No probability or identification confidence is inferred from these summaries.

## IonNaming

The functions live in `openms::chemistry::ion_naming`:

| Function | Behavior |
| --- | --- |
| `charge_suffix(i32) -> String` | Zero gives an empty suffix; magnitudes 1–8 give repeated signs; larger magnitudes give a sign and number |
| `charge_from_name(&str) -> i32` | Read source caret, trailing-sign or sign-number notation; unknown/malformed values give zero |
| `with_charge(&str, i32) -> Result<String>` | Keep an existing nonzero named charge; otherwise insert the supplied suffix at the end of the first line |
| `ordinal_from_name(&str) -> u32` | Read up to nine digits immediately after an ASCII letter; otherwise zero |

Charge parsing inspects only the first CR/LF-delimited line. The last valid
caret token has priority over suffix parsing; mass-delta or confidence fields
must not become numeric charges. `y3-H2O1+` has ordinal 3, and crosslink text
such as `[alpha|ci$y3]++` has no ordinal in the recognized position. UTF-8 free
text after the first line is preserved. These helpers implement the source
conventions and do not claim complete mzPAF parsing.

Parsers allocate no strings and use checked integer conversion, including i32
extrema. `with_charge` limits copied output to 1 MiB, including unchanged returns
and free-text lines. A suffix alone occupies at most eleven bytes. Immutable
anonymous modification spellings are also shared between owned peptide slices,
so fragment generation does not repeatedly copy those strings.

## Resource accounting and evidence

One annotation operation permits at most one million measured peaks, 1,000 hits,
one million top-N entries and 100,000 inspected precursors. Annotation stages
share 50 million work units, 64 MiB of copied labels and a conservative 128 MiB
cumulative allocation-payload allowance. Generated spectra and alignment retain
their own storage limits.

Every candidate shares 50 million theoretical residue-work units, ten million
neutral-loss declaration visits, 100 million fine-isotope work units, 50 million
coarse convolution products and 50 million alignment initialization/cells.
Precharged counters or actual inner-loop counters are consumed before the work.
Coarse convolution now also shares its allowance across all envelopes of a
standalone theoretical-spectrum call, closing the earlier per-envelope reset.
Ordinary standalone alignment retains its configured `max_cells`; standalone
generation and isotope operations start fresh allowances.

Focused tests cover source labels/statistics, ppm duplicate semantics, source
flag behavior, checked undefined statistics, sorting, limits and atomic failures.
Independent references retain original measurement bits and distinguish rounded
source assertions from directly derived nonzero errors. [Workflow tests](../tests/spectrum_annotation_workflow.rs)
connect modified digestion, distinct chemistry with equal sequence text, typed
identification metadata, ion-name display and both XML formats. The C++ library
was inspected but not built or executed. The full core port remains in progress.
