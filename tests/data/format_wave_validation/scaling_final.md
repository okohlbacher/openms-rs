**Verdict: no remaining concrete bug in the batching path. Earlier conditional findings are discharged.**

Discharged findings, with the excerpt that closes each:

- **Terminal lookup mismatch.** The guard now runs the same two-step specificity search as `resolve_terminal`, on the same first/last residue, and bounds the returned handle rather than the plan's record. Any `find` hit that fails to yield a handle falls back to the sequential path, where `set_terminal` would take the mass-tag fallback anyway. Consistent.
- **Anywhere leaking into a terminal slot.** Planner restricts terminal slots to N/C-term specificities, and `set_modification_with_registry` explicitly requests `Anywhere` for residues. The bulk handle lookup requests `Anywhere` too, so bulk and sequential resolve identically.
- **Absolute-formula records.** `bound_atom_contributions` refuses to batch when the diff formula is empty but an absolute formula exists. Both `internal_formula` and `internal_mass` prefer the diff formula when nonempty, so the bound models the same chemistry the rebuild uses.
- **Regrouping hiding a transient deficit.** `lower` equals the original formula plus every negative part. Any subset of deltas applied in any order is elementwise at least that, so a nonnegative `lower` proves every sequential intermediate nonnegative. The late-donor case yields a negative `lower`, forces the sequential path, and rejects there. The early-donor case is handled the same way, since terminal positives are never credited to `lower`. Matches the two tests you describe.
- **Mass tags disabling the composition check.** Once any residue mass tag is set, `rebuild_chemistry` yields no formula and the negative-atom check is off for the rest of the sequence, in both paths. The bulk path only ever falls back to sequential in that situation, never the reverse, so behavior is unchanged and the cost is only a perf fallback.
- **Terminal delta ordering.** Terminals are applied last but in plan order, and `rebuild_chemistry` sums both terminal deltas before residues. Mass overflow checks are therefore hit in the same order as the sequential path.
- **Overflow in `upper`.** Now rejects with `Ok(false)`, so an oversized positive sum can no longer batch.

One remaining observation, scoped rather than conditional. When `tag_count` is nonzero, `annotated` is rebuilt from bare residue characters, so any modifications already present on `sequence` are dropped, while `lower` is seeded from `sequence.formula()` which includes them. When `tag_count` is zero, `result` is a clone and keeps them. If the pepXML parser always hands in a plain peptide, this is unreachable and not a bug. If a caller can pass a pre-annotated sequence, the two branches diverge and the bound is seeded with atoms the result no longer has. A one-line `debug_assert!` that `sequence` carries no modifications would pin the assumption.

Minor, non-blocking: `length` counts residues as `size_of::<char>()` but the string is UTF-8, so the meter slightly over-charges for ASCII peptides. Harmless.

qcML: as you state, `register_set` never validated member-to-run membership before or after, so no regression there.
