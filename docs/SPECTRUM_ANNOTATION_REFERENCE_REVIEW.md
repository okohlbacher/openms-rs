# Spectrum annotation reference review

The [independent reference suite](../tests/spectrum_annotation_reference.rs) covers SpectrumAnnotator and selected IonNaming behavior from OpenMS4-core revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. The [manifest](../tests/data/spectrum_annotation_provenance.json) records eighteen source hashes and the extraction method. No C++ code was built or executed. Tests use packaged literals and direct arithmetic, with no runtime dependency on the source checkout, Python or temporary preparation files.

## Source and derived evidence

| Fixture | Rows | Meaning |
| --- | ---: | --- |
| [Source peaks](../tests/data/spectrum_annotation_source_peaks.tsv) | 13 | The source's eleven IFSQVGK peaks in original input order, plus its two unmatched peaks. Original decimal m/z and intensity tokens, binary64 m/z bits, binary32 intensity bits and expected labels are retained. |
| [Source statistics](../tests/data/spectrum_annotation_source_statistics.tsv) | 17 | Literal numeric, boolean and string assertions, with source assertion kind and line number. |
| [Derived masses/errors](../tests/data/spectrum_annotation_derived_masses.tsv) | 11 | Independent calculations from source element masses, full residue formulas, proton mass, water conversion and prefix/suffix summation. These are not captured C++ output. |

The original experiment uses peptide IFSQVGK, charge two, b/y ions with metadata, and absolute tolerance 0.1 Da. Its eleven measured intensities are `1.1f`, so the exact widened total is `11 * f64(f32(1.1))`, not decimal 12.1. The two unmatched measurements have intensity `0.5f` and m/z 100/1000. Matching tests retain the expected sorted labels `y1+,y2+,b2+,y3+,b3+,y4+,b4+,y5+,b5+,b6+,y6+`.

The source metadata equality assertion at 12.1 is not treated as bitwise numeric evidence. Its `topN_MSEfragmenterror` assertion at zero is also weak: the measured peaks differ from predicted masses, and their actual mean squared error is positive. The suite preserves those literals with explicit tolerances, then separately checks f32 intensity widening and positive error statistics from independent mass calculations. Source formula iteration uses element pointers and can change the last binary64 rounding bits; independently predicted error statistics therefore use absolute tolerance `1e-9`, not a claim that Python reproduces C++ mass bits exactly. Missing ClassTest default relative/absolute behavior is not inferred.

A separate four-peak construction has exact binary errors `[0.125, 0.25, 0.5, 1]` and intensities `[1, 2, 4, 8]`. It gives independently calculable statistics:

- Top seven errors are `[1, 0.5, 0.25, 0.125, 0, 0, 0]`: mean `1.875/7`, MSE `1.328125/7`, and sample variance `(1.328125 - 1.875²/7)/6`.
- Median error is 0.375; the source quartile indices give IQR 0.75.
- All peaks matched gives match-based S/N zero. The median-intensity groups give S/N four.
- Prefix/suffix intensity ratios are `4/15` and `11/15`.

These exact binary inputs check padding, sample rather than population deviation, intensity order and group assignment without copying production results into an oracle.

## Observable annotation behavior

`annotate_matches` sorts measured peaks by m/z, replaces all three kinds of auxiliary arrays with `IonNames`, `Charges` and `IonMatchError`, and preserves unrelated metadata. `IonNames` is the actual source name despite the header's singular spelling. Errors are absolute Daltons narrowed to f32, even for ppm matching. The measured f32 intensity is widened for hit peak annotations; theoretical intensity is not substituted.

A large-ppm, one-measured-peak example independently produces two charge hypotheses. The last theoretical match wins in per-measured-peak arrays and in the include-unmatched annotation branch. Matched-only annotations instead retain both matches, including repeated measured m/z and intensity. This is intentional source behavior, not nearest-error selection or annotation concatenation.

Statistics first sort by intensity and leave the final hit's arrays on the measured spectrum. Equal-intensity source ordering is not specified by `std::sort`; the native implementation uses stable deterministic ties. If enabled precursor statistics has a precursor to inspect, it sorts by m/z again. The source passes the raw tolerance number to an absolute-Dalton precursor lookup even when alignment uses ppm; the independent test preserves that distinction. Spectrum metadata uses `fragment_mass_tolerance` and its ppm flag, while identification metadata uses `fragment_match_tolerance`.

Both `precursor_in_ms2` and `fragment_mass_tolerance_ppm` pass a C++ boolean to `setMetaValue(const DataValue&)`. There is no `DataValue(bool)` constructor, so integer promotion stores 0 or 1. Independent assertions require that exact integer type on peptide hits; the native spectrum's string metadata projects it as `"0"` or `"1"`. The fixture's `bool` kind preserves the source assertion literal, not a boolean or string storage type.

The stored `list_of_ions_matched` option does not control the source matched-ion list; `basic_statistics` does. The stored `fragment_error_statistics` option does not control the error block; nonzero top-N does. The binary-error test sets both unused options false and still verifies their source results.

## Parsing and statistical boundaries

SpectrumAnnotator uses complete-name regex matches, not IonNaming's broader charge/ordinal parser. Its N-terminal and C-terminal character classes contain literal commas as well as a/b/c and x/y/z. Only digits followed by zero or more plus signs qualify for the terminal-current ratios. Losses, caret notation and numeric charge suffixes do not qualify. The series-position regex additionally accepts a sign/comma run, a word-character suffix, and final plus signs; neutral-loss labels can therefore extend a series. The series map still requires a/b/c/x/y/z and resolves longest-run ties in lexical type order. Raw internal-fragment sequence labels and immonium labels are not conventional numbered terminal-series labels.

The source accumulates matched intensity/error before testing an out-of-range series ordinal, but excludes that label from the matched-ion list on its later `continue`. This can make list counts differ from accumulated matched contributions. The native grammar and small synthetic private-helper checks must preserve this ordering where the public generator can produce such a label; no test-only public injection API is added.

For matched counts one, two and three, the source's chained `nth_element` subranges are invalid. Native sorting of the intended order-statistic indices is explicitly a safe extension. The independent suite checks these intended values and does not call them upstream numerical goldens. Nonempty top-N one requests an undefined sample deviation; native code returns an atomic error. Enabled ratios that become nonfinite likewise require checked errors rather than invented finite scores. An error must leave both the identification and spectrum unchanged.

IonNaming is checked separately against source literals and conservative branch examples: last-caret priority, a valid zero-caret short circuit, malformed-caret fallback, same-sign suffixes even after slash/star fields, suppression of sign-plus-number fallback after those fields, CR/LF first-line boundaries, i32 extrema, ordinal length limits and loss digits that must not become ordinal digits. UTF-8 free text is preserved by insertion at the ASCII line boundary. These helpers cover the source naming conventions; they are not a complete mzPAF parser. ASCII character recognition is explicit rather than dependent on the process locale.

## Independent production review

The review reads scientific operations and scoped resource accounting without modifying production. Initial inspection identified two bookkeeping gaps while the implementation was still being written: allocation of the include-unmatched last-match map and sort scratch accounting by the wrong element type. Both were corrected by the implementation owner and verified on final review. The completed review also checks cumulative annotation allocations, label copying, scoped input handling, sorting and atomic commit boundaries. Shared generator/alignment hooks are invoked with the same operation counters across hits; their internal accounting has a separate integration review. A subsequent constructor-level review found that two boolean arguments had incorrectly been represented as strings; both were corrected to source integer semantics, and the independent tests now reject the previous representation. No remaining actionable discrepancy was found in this scope.

All seven reference tests and strict Clippy checks pass on the current compiler with all features and on Rust 1.85 with default features disabled. The fixture audit verifies all eighteen source hashes, three fixture hashes, forty-one rows, original measured input order and sorted labels, source assertion locations and recorded numeric bits.

This evidence establishes behavior for the stated source utilities and native checked boundaries. It does not validate a search score, identification confidence, or biological assignment, and does not claim that the full OpenMS core port is complete.
