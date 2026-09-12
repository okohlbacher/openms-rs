**Review scope:** the supplied diff only, no tools. No blocking defects in the guards themselves. Four findings, ordered by risk.

**1. mzData: trailing self-closing element is not covered by the multi-root guard (test gap, possible defect).**
The new check lives only in the `Event::Start` arm. quick-xml emits `Event::Empty` for `<x/>` unless the reader has `expand_empty_elements(true)`. The new test uses a self-closing root twice and passes, so either expansion is on or the Empty arm has its own guard. The diff shows neither, so the test does not justify the assertion for the general case.
- Trigger: `<mzData>…</mzData><x/>` or `<mzData/><x/>`.
- Correction: apply the same `seen_root && self.stack.is_empty()` check in the Empty arm (or confirm expansion is enabled in `Parser::new`) and add `format!("{root}<x/>")` to the loop in `one_xml_document_excludes_trailing_roots_and_outside_character_data`.

**2. Percolator: `max_lines = max_rows + 1` regresses valid inputs, most clearly for Sage siblings.**
`tab_options` is shared by the `.pin`, the Sage `.tsv`, and `matched_fragments.sage.tsv`. The fragment file holds one row per matched ion, routinely 10 to 50 rows per PSM. A caller that sizes `max_rows` to its PSM count now has the fragment file refused where the old independent CSV default accepted it. The same clamp also bites the `.pin` itself when the physical line count exceeds the row count: a trailing newline that the CSV layer counts as a line, or any comment line, refuses a file with exactly `max_rows` data rows.
- Trigger: `load(path, ReadOptions { sage_annotation: true, max_rows: <psm count>, .. })`; or a `.pin` at the row ceiling with one comment line.
- Correction: drop `max_lines` from the clamp. The existing `row > options.max_rows` checks already bound rows, and `max_input_bytes` now bounds staging, which is the stated goal. If a line cap is wanted for the `.pin`, pass a separate `csv::Limits` for the fragment file and use `max_rows + 2` or a comment-aware count for the others.

**3. Percolator: error type shifts at the ceiling.**
With the clamp, a `.pin` holding `max_rows + 1` data rows fails inside `CsvFile` with the CSV limit error before `read_csv` reaches its "row limit exceeded" `InvalidValue`. Callers that match on the variant or message change behaviour. Tests pass, so no current test pins this, which means the boundary is untested either way.
- Trigger: header plus `max_rows + 1` rows, `max_rows` small.
- Correction: resolved by finding 2. Otherwise add one test asserting the variant at `max_rows` and `max_rows + 1`.

**4. mzIdentML: the reversed `BetaPepEv` doc claim is not verifiable from the diff and is worded incorrectly.**
The new comment says a `std::string += Int` overload in `StringUtils.h` calls `std::to_string`. No overload can be added to `std::string`; only `OpenMS::String` carries `operator+=(Int)`. Whether the source appends digits or a code point therefore depends on the declared type of the accumulator in `OPXLHelper`/`MzIdentMLHandler`. If it is `std::string`, the deleted paragraph was correct and this revision introduces a wrong statement. Doc-only, no behavioural effect.
- Correction: verify the variable type at the call site and cite `String.h` rather than `StringUtils.h`; keep the divergence bullet if the type is `std::string`.

**Minor, not blocking.**
- Encoding checks in Mascot and mzIdentML accept any case of `UTF-8` and refuse the `utf8` alias. Refusing as `Unsupported` is a safe conservative choice.
- The Mascot encoding test lacks the `contains("ISO-8859-1")` guard the mzIdentML test has. If `document()` ever changes its declaration casing, the test still fails loudly rather than passing vacuously, so no change needed.
- mzData outside-root text trims only ASCII whitespace. If the reader does not strip a UTF-8 BOM before the first event, a BOM-prefixed file is now refused. quick-xml strips it in current versions, so this is a version pin concern only.

**Test justification summary.**
- Mascot: the three reference cases and the ISO-8859-1 case exercise exactly the new branches. Justified.
- mzIdentML: the encoding case is justified.
- mzData: the listed cases are justified; the Empty-element case is missing (finding 1).
- Percolator: the staging test proves the byte ceiling and the zero-option rejection. Nothing exercises the new line clamp or the Sage sibling path under it (findings 2 and 3).
