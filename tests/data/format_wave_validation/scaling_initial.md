Write is disabled, so the findings go directly in this message. Review complete.

## Findings

**1. The batching gate excludes far more than atom-removing dependent records** (`src/format/pepxml.rs:2580-2589`)

- The negative-atom test rejects ordinary positive-delta modifications. Deamidated is `H(-1) N(-1) O`, and every SILAC label is `C(-6) 13C(6)` plus nitrogen variants (unimod.obo lines 324, 4256, 5008, 5144). Concrete trigger: a search with fixed Label:13C(6) on K sends every K-containing hit through the per-setter path, so cost becomes length times sites again. Real peptides stay under budget, but the quadratic is back for adversarial inputs spelled with those names.
- The `absolute_formula().is_some()` clause is order-independent and should go. Bundled TSV records always carry `None` (`modifications.rs:514`). Only OBO `Formula:` records get `Some` (`obo.rs:285`), and the rebuild consults it only when the diff formula is empty, where it replaces the base composition without touching other annotations (`sequence.rs:959-968, 977-997`). A caller loading PSI-MOD into a custom registry would lose batching for every record.
- The actual hazard is narrow. It exists only when residues are installed before terminals, the final formula is `Some`, and the residue-only intermediate is negative. Any mass tag makes the formula `None` (`sequence.rs:958, 1037`) and the negativity check at line 1061 never runs.

Simplest correct alternative: put the terminals into the same `with_resolved_modifications` call, which already accepts N- and C-terminal handles (`sequence.rs:678-683`). Resolve Known terminals with the same two-specificity lookup `resolve_terminal` performs (`sequence.rs:1490-1498`, both `find` and `get_modification_handle` are public, and the same loop already exists at `pepxml.rs:1805-1824`). Respell Mass terminals as a `.[text]` prefix or suffix, which the parser resolves through the identical `resolve_mass` call `set_terminal` makes (`sequence.rs:1309-1318` versus `640-649`). That yields one rebuild on the final state. If it errors, return `Ok(None)` and let the ordered path reproduce its per-setter fallback. The gate is then unnecessary, and the result matches the source, which installs pointers and never fails (`PepXMLFile.cpp:2081-2121`).

**2. An existing divergence the gate does not cover** (`pepxml.rs:2412-2423` with `2708-2733`)

A fixed header N-term Dimethyl with `variable="N"` plus three explicit Unknown:177 (D) annotations plans residues before the terminal. The third `set_residue` rebuild goes to H=-4, the named setter errors, and the residue becomes a mass tag with a diagnostic. The source produces three named residues. The fix in finding 1 covers this. Worth adding as a test beside `tests/pepxml.rs:1410`.

**3. Batched path errors hard where the ordered path falls back** (`pepxml.rs:2686` versus `2720-2727`)

`with_resolved_modifications(...)?` fails the whole read. The gate blocks negative-count cases, but `checked_add` overflow (`chemistry/mod.rs:353`) remains reachable with a custom registry. The fallback-to-ordered step in finding 1 closes it.

**4. Accounting: no bypass found**

- The gate pre-scan charge (line 2582) is sunk on bail, and `charge_rebuild` at line 2668 is charged even when the sequence is only cloned. Both over-charge, harmless.
- The fixed-residue scan charge (2443-2448) equals the real `contains` cost per declaration. `lookup_from_header` (2788-2790) charges the origin length only on a mass match, which is when it is scanned. `start_modification_info` charges one per candidate (1774).
- `tighten` runs after decode and only lowers (2913).
- Pre-existing and not introduced here: every registry mass search walks all entries (`modifications.rs:860-869`) but is charged 64 at line 1972, and again inside the respelled parse via `lookup_mass` (`sequence.rs:1416, 1419`) charged only `text.len()` at line 2661. Constant factor, bounded by registry size.

**5. Plan order versus source: no regression.** Explicit modifications keep encounter order (`PepXMLFile.cpp:2081`). The `taken` flags mirror `isModified()` and the terminal checks. The C-term copy-paste divergence (2425-2427 versus cpp line 2135) is deliberate and documented.

**6. Respell condition** (`pepxml.rs:2485`). Attributes with `.`, `n`, or `c` terminal tokens differ from `as_str()` and silently take the ordered path. Cost only, no correctness issue.

**7. qcML commit is equivalent** (`src/format/qcml.rs:2696-2745`). `check_parameter` and `check_attachment` (2045-2058) are exactly what `push_parameter` and `push_attachment` run (899, 919). Per-entry ceilings are enforced at 2632, 2660, and 2676 before each push. ID text and entry ceilings live in `register` (808-822). The source also registers first and then adds (`QcMLFile.cpp:956-995`). One nit: `register_set` clones the owned member set (2738 into 802); passing by value avoids the copy.

**8. Tests**

- `tests/pepxml.rs:1269` asserts wall-clock time under 60 seconds. The budget assertion is the real pin. Expect CI flakes under load if it stays.
- `tests/pepxml.rs:1410` pins only the terminal-first case. Add the reverse case from finding 2.

I could not create the plan file because the Write tool is disabled in this session.
