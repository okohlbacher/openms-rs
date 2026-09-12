**Revised dispositions**

1. **Line clamp — withdrawn (resolved).** The clamp is gone. `max_bytes` now bounds input, storage, and line size independently via `csv::Limits`, and the post-CSV `max_rows` checks remain. No further concern.

2. **max_rows vs. comment lines — resolved.** The regression accepting one data row with leading and trailing comments at `max_rows=1` covers the case I raised. Accept.

3. **Sage annotation rows limited by max_rows — retracted as a new finding.** The per-row check in the annotation loop predates this patch. Fragment files legitimately have many more rows than PSM files, so a shared `max_rows` can reject valid pairs, but that is a preexisting behavior and out of scope here. Non-blocking follow-up at most.

4. **mzData second root — resolved.** With `expand_empty_elements=true`, a trailing `<x/>` after the root becomes an `Event::Start`, which hits the "multiple mzData XML roots" check because the stack is empty and `seen_root` is set. The added test exercises exactly that path. Accept.

5. **Global `operator+=(std::string&, int)` — retracted.** My claim that such an overload cannot exist was wrong for the pinned source. The header block at revision bc9cc12 defines global nonmember `operator+=` overloads for `std::string` against every integral and floating type, plus `operator+` built on them. Two evidence tiers, kept distinct:
   - **Source review:** `StringUtils.cpp` implements `appendToStr` with `std::to_chars` into a stack buffer and appends the decimal characters. That is what rules out control-byte output.
   - **Isolated probe (tier 2):** confirms unqualified lookup from `OpenMS` member and `OpenMS::Internal` free-function contexts resolves to these globals and produces decimal text. It substitutes adapter bodies and has no SDK linkage, so it validates overload resolution only, not the production bodies. The provenance JSON already labels it correctly.

6. **Comment wording — accepted.** "appendToStr uses to_chars" matches the source.

**Remaining blocking issues in the revised patch:** none.

**Non-blocking notes**
- Keep the tier-2 label on the probe wherever it is cited so nobody later reads it as an executed SDK differential.
- `max_line_bytes = min(cap, max_bytes)` allows a single line to consume the whole byte budget. Fine, since input bytes are bounded by the same value, but worth one comment if anyone later raises `max_bytes` independently.
- The shared `max_rows` for PSM and Sage annotation files is the one preexisting sharp edge worth a follow-up ticket, not a change in this patch.
