# Independent Tagger reference review

The [reference suite](../tests/tagger_reference.rs) preserves all six exact tag counts and all 120 tag-presence/absence assertions in the first section of the pinned [Tagger class test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/Tagger_test.cpp). The [manifest](../tests/data/tagger_provenance.json) records 14 source hashes, fixture hashes, source line numbers, literal trace values and separately identified derived checks. No C++ build or execution generated this evidence.

## Source integration goldens

The source theoretical spectrum enables a/b/y ions, first prefix ions, neutral losses and precursor peaks, without isotope expansion or metadata. The unmodified peptide `PEPTIDETESTTHISTAGGER` at fragment charges 1 and 2 produces 357 peaks. The D-oxidized peptide at fragment charge 2 produces 180. These peak counts are checked before Tagger assertions, so a dependency mismatch is visible.

| Spectrum | Tagger charge range | Modification table | Exact tag count |
| --- | --- | --- | ---: |
| Unmodified | 1 | Ordinary residues | 890 |
| Unmodified | 2 | Ordinary residues | 1006 |
| Unmodified | 1–2 | Ordinary residues | 1094 |
| D-oxidized | 1–2 | Ordinary residues | 545 |
| D-oxidized | 2 | Fixed Oxidation (D) | 667 |
| D-oxidized | 2 | Variable Oxidation (D) | 739 |

All use minimum/maximum tag lengths 2/5 and 10 ppm. [Count rows](../tests/data/tagger_source_counts.tsv) and [membership rows](../tests/data/tagger_source_membership.tsv) retain their original source locations. The presence/absence cases include broken modified-residue paths, restored modified paths, first-prefix omissions and charge-dependent false-positive suffix tags. The vector and spectrum overloads must agree exactly, and returned tags must be strictly lexicographically ordered.

These are native TSG-to-Tagger tests against published source assertions. The native-generated spectra are not mislabeled as saved C++ output, and matching the counts alone is not claimed to establish every possible Tagger path.

## Independent small and boundary cases

The source's small P/E traces retain their original f64 addition expressions and rounded anchors, including 0.015 and 0.025 Da offsets. At 0.02 Da, PE/EP/PEP are present in the expected traces and the outside-tolerance PE path is absent. A common offset cancels in later differences; it is not added to every edge.

Custom records independently construct exact internal masses 32, 34, 36 and 64 from an absolute **free-residue** mass plus a nonzero difference marker. The source free-mass-minus-water operation then yields those keys. They isolate behavior that broad proteomic tolerances cannot distinguish:

- An exact excluded lower-bound candidate causes the whole query to fail, even when a later exact match exists. This retains the source's immediate rejection after `lower_bound`.
- Equal nearest errors retain the lower-mass candidate. A successor examined after iterator increment cannot replace a closer valid entry when it lies outside the window.
- Zero tolerance produces no matches, including at an exact residue mass. Negative tolerance is normalized with absolute value.
- Exact mass collisions retain the last inserted origin according to the caller's modification list order.
- A selected modified L still produces both L and I branches. Two such gaps yield II, IL, LI and LL. A selected I produces only I.

The source only returns the nearest single table entry. These tests do not substitute an all-residues-within-tolerance matcher.

## Finite source behavior and deliberate native boundaries

Input m/z order is preserved, including unsorted and signed values. The independent unsorted trace places an overlarge early gap before a valid P edge and verifies that source early pruning suppresses the path. A separate zero-charge case with broad absolute tolerance verifies selection of G from zero neutral gaps. No intensity or metadata filtering is inferred from the fact that one overload accepts MSSpectrum.

The append tests distinguish `min_length > peak_count`, which leaves existing output untouched, from `min_length == peak_count`, zero maximum length and inverted charge ranges, which still sort/deduplicate existing strings. A later maximum-charge setter changes the search range without reconstructing the mass table.

ResidueDB constructs B/Z/X from empty free-residue formulas, giving zero free mass. Tagger's private residue table therefore has mass minus water for an unchanged modification at those origins. Independent descending-position tests retain this source geometry. It does not change general AASequence's explicit unavailable-mass policy or invent a composition for ambiguous peptide residues.

The source's two modification resolutions are preserved: unrestricted initial name lookup, then short name at the original origin and ANYWHERE. For ambiguous results, native provider order is the explicit deterministic replacement for allocation-dependent source pointer order. Public ModificationsDB ambiguity handling remains unchanged. A private empty-short-ID lookup scans actual anonymous records in provider order, without exposing empty aliases through the public registry. Unknown origin codes error instead of dereferencing a null source residue pointer.

Nonfinite arithmetic, allocation/work exhaustion and unsafe inclusive maximum-charge overflow are checked boundaries. An explicit DFS stack replaces unbounded recursion. Source output ordering, finite unusual inputs and the private signed-mass table remain expressible within those limits.

## Misleading source comments excluded from the contract

The implementation uses `(tolerance / 1e6) * neutral_gap_mass` for ppm and a literal dalton value for absolute mode. Some newer class-test prose describes absolute arguments of 10 or 100 as ppm scaled by fragment m/z. Those explanations are incorrect and are not carried into the native API. A synthetic trace also labels 103.0094 as I, though it is approximately cysteine; its literal numbers remain source data, not corrected residue annotations.

Free-residue formula mass followed by water subtraction is retained. Direct internal-formula mass or separate addition of a modification delta can differ in floating order. C++ formula maps use pointer keys, so their atom summation order is not a portable bitwise promise; the native helper retains deterministic formula arithmetic and the source's subsequent subtraction.

## Validation and production review

All six independent tests pass in the final unified commands on Rust 1.96 and 1.85, both with all features and without default features. Those commands also pass all 41 private unit tests, eight direct Tagger tests and the enabled Tagger workflows. Strict all-target Clippy, rustdoc and formatting pass. Both source input sizes, all six tag counts and all 120 membership assertions matched unchanged. Strict focused Clippy checks also pass for both compiler/feature combinations. All 14 pinned source hashes, both fixture hashes, source line references and local document links were verified.

Read-only review covered exact matcher bounds, explicit DFS/backtracking and I/L replay, cumulative path/text/work limits, source-ordered modification replacement, and atomic sorting/commit of existing and new output. No remaining scientific or atomicity defect was identified. The constructor checks modification-name allocation bounds and resolves owned numeric values without retaining the caller registry. Source count and membership goldens were not weakened to accommodate a dependency difference.
