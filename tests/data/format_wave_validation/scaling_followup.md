**Verdict:** the bound is sound for what it models, but two fallback paths still install something the bound never described, and one exclusion condition is too narrow. Nothing in the qcML change loses validation that the parser did not already own, with one thing to verify.

**Peptide fast path**

- **Terminal guard does not mirror the terminal setter.** The guard looks the record up with the record's own `term_specificity()` and the terminal residue. The setter presumably asks for an N-term or C-term specificity. A record resolved as `Anywhere` from a mass-only pepXML entry passes the guard, so its delta is subtracted, but the setter rejects it by name and installs the mass-only tag instead. That tag is resolved at install time against a formula that already contains every residue addition, not the ordered prefix. If it resolves to a different isobaric record with a larger negative part, the fast path accepts where the ordered path would have rejected or fallen back further. This is the original grouping bug surviving through the fallback hole. Fix: treat any terminal entry whose name lookup with the setter's specificity fails as ineligible, or replicate the setter's exact lookup here.

- **Residue name lookup may be stricter than the ordered setter.** The fast path resolves with `Some(TermSpecificity::Anywhere)`. If `set_modification_with_registry` accepts terminal-specific records on the first or last residue, the same plan yields a name in the ordered path but a mass tag plus a spurious diagnostic in the fast path. That is a visible output difference, not just a performance one. Verify which specificity the setter uses and pass the same value.

- **Absolute-formula exclusion is too narrow.** A record with a non-empty diff formula and an absolute replacement formula is not excluded. If chemistry applies the absolute formula when present, the effective delta is absolute minus residue formula, which can be more negative than the diff's negative part, and the bound under-subtracts. Exclude on `absolute_formula().is_some()` alone unless the chemistry provably ignores it when a diff exists.

- **Per-residue mass floor is not part of the bound.** The bound proves atom counts only. The ordered path also rejects a residue whose mass goes negative. This is safe only if bulk installation and the respelled parse run the same per-residue check and return an error, which the comment asserts but the excerpt cannot show. Worth one test: a record with a diff mass larger than glycine, applied to glycine, must reach the ordered path in both flows.

- **Mass-only fallback for a Known residue record can change resolution.** The fallback spells `diff_mono_mass()` as text and the parse resolves it by mass at that residue. If it resolves to a different record with a different formula, the bound loop over `result` correctly picks up the parsed record, so this case is covered. No finding, noted as verified.

- **Terminal mass tags disable the fast path entirely.** Any terminal `ResolvedModification::Mass`, known or unknown, returns `Ok(None)`. Correct, but every peptide with a terminal mass annotation pays the ordered cost. This is a documented limit, not a bug.

- **Duplicate slot entries are silently merged.** Two plan entries for the same residue overwrite one tag and increment `tag_count` twice, and a handle plus a tag on the same position would apply both. The planner promises first-per-slot, so this is a latent invariant rather than a live bug. A debug assertion would cost nothing.

- **Bracket-safe strings that do not parse.** Text like `+` or `1.2.3` passes `bracket_safe`, fails the parse, and falls to the ordered path, which rejects identically. No acceptance difference, only one wasted parse.

**qcML commit**

- **Nothing validated by the old single-child setters is skipped** as far as the excerpt shows, provided `register_run` and `register_set` still length-check the ID and name. Those are the only fields not passed through `check_text` here. Confirm they do, or add the check before registration.

- **Set membership is not validated.** Members are passed to `register_set` unchecked against known run IDs. If the old setter path never checked either, this matches the source. If it did, that check is now gone. Verify against the previous implementation.

- **Attachment references are not cross-checked** against quality parameter IDs. The C++ source does not check them either, so this is parity, not a regression.

- **Error mapping drops line numbers** for any error kind other than `InvalidValue` and `MissingInformation`. Cosmetic.

**Remaining limits**

The bound cannot prove anything about the ordered path's per-residue mass check or about records whose installed chemistry differs from the record subtracted. Both terminal findings above are instances of the second. Closing them requires the guard to call the same lookup the setter calls, then subtract the record that lookup returns rather than the planned one. Until then the fast path is correct only when every terminal name resolves identically under both lookups, which the code assumes but does not enforce.
