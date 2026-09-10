# Peptide-property reference review

This review covers the native AAIndex, gas-phase basicity, hydrophobicity and isoelectric-point calculations against OpenMS4-core revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. Expected data come from pinned source declarations and independent arithmetic. No C++ library was built or executed, and no native production output generated a reference value.

The [reference suite](../tests/peptide_properties_reference.rs) reads five small TSVs. It has no dependency on the source checkout, temporary preparation files, Python, a network service, or an additional Rust dependency. The [manifest](../tests/data/peptide_properties_provenance.json) records sixteen source hashes, original source links, extraction details and fixture hashes.

## Source constants and assertions

| Fixture | Rows | Evidence and checks |
| --- | ---: | --- |
| [AAIndex and GB](../tests/data/peptide_properties_aaindex_gb.tsv) | 262 | All 200 public scale values and 62 protected GB lookup entries; original decimal spellings, binary64 bits and source lines. |
| [Hydrophobicity](../tests/data/peptide_properties_hydrophobicity.tsv) | 182 | Seven 26-letter source rows: 140 usable values and 42 `999` sentinels. Sentinels require errors, not numeric scores. |
| [pKa values](../tests/data/peptide_properties_pka.tsv) | 59 | Four nine-value base tables, twenty Bjellqvist N-terminal overrides, two C-terminal overrides and selenocysteine pKa 5.73. |
| [Source scalar assertions](../tests/data/peptide_properties_source_scalar_literals.tsv) | 30 | Three GB, fourteen pI/charge and thirteen hydrophobicity/profile/window/moment assertions. |
| [Derived GB arithmetic](../tests/data/peptide_properties_gb_decimal.tsv) | 100 | Separately derived, 90-digit Decimal calculations over ten sequences and ten temperatures. These are not additional upstream goldens. |

Extraction preserves original decimal tokens rather than replacing them with formatted native values. The hexadecimal column is the IEEE binary64 representation obtained by parsing that token. Source comments are removed only for token matching; line positions remain intact. Independent preparation checks compared all 181 literal AAIndex lookup assertions, seven alanine hydrophobicity assertions and the class test's separate twenty-value Eisenberg table with the extracted declarations. Regeneration was byte-identical, and every source/fixture hash and decimal/bit pairing was verified.

All 200 public AAIndex values and all 182 hydrophobicity cells are directly tested. The 62 protected GB entries remain private in the native API. Their values feed an independent source-expression evaluator used for every canonical singleton and all 400 ordered canonical pairs at 100,000 K. At that temperature weak terms remain numerically visible beside the strongest site; this is more discriminating than only testing arginine-rich peptides at 500 K. Additional ordinary-temperature cases distinguish first-sidechain omission and tied strong sites. Read-only production review also checks the source table column mapping.

Every pKa constant is checked through the public charge calculation at its pKa and at pKa ±0.75. Terminal groups are isolated by annotating the opposite terminus; `X` supplies a neutral sidechain when a base terminal value is needed without a Bjellqvist override. Sidechains are isolated by annotating both termini. Override checks include the separately calculated sidechain contribution. The independent formula uses acid/base concentration fractions, rather than copying the implementation's reciprocal-of-one-plus-power expression. Selenocysteine is checked in all four scales. The absolute comparison tolerance is `2e-14`.

The scalar fixture records explicit source `TOLERANCE_ABSOLUTE` context. The pinned first-party snapshot does not contain the ClassTest implementation needed to establish default relative tolerances or their combination with absolute tolerance. The Rust reference suite therefore states its own comparisons: GB absolute `0.01`; nonzero pI absolute `5e-4`, accounting for rounded references and the default bisection interval; zero-charge residual `1e-3`; hydrophobicity scalars `1e-8`. These are not claimed to reproduce unavailable test-macro defaults. The five EMBOSS literals are attributed by the source test to an independent calculation; this suite preserves those literals and attribution without claiming to have rerun EMBOSS.

## Numerical and scientific distinctions

AAIndex predicates preserve unusual source memberships: aliphatic includes F/G, basic includes W, and the positive-charge index FAUJ880111 instead selects H/K/R. Unknown/lowercase characters yield zero for predicates but errors for numeric scales. All seven hydrophobicity scales require canonical residues; pI accepts B/Z/X/J as neutral sidechains, handles U separately and rejects O. Properties use parent residue codes and do not need a sequence formula or mass.

The source GB expression retains its old gas constant, `R = (6.0221367e23 * 1.380657e-23) / 1000`, backbone/sidechain split pairing, sequential accumulation and final division by `ln(2)`. The first residue's sidechain is omitted; a later sidechain value of zero still contributes `exp(0)`. These choices are observable and are not replaced with a different chemical model. Ordinary finite evaluation is compared with the independently read source expression within four ULPs to allow platform math rounding.

The GB Decimal grid uses exact binary64 source energies and exact binary64 R/T inputs, high-precision multiplication, exponential/logarithm evaluation and `ln(2)` at ninety decimal digits. Shifted exponents at or below −10000 are omitted, with a total bound below the number of sites times `exp(-10000)`. The result is then rounded to binary64. The grid covers empty input, A/K, AK/KA, three source peptides, ARR and all twenty canonical residues, from the smallest positive subnormal temperature to the largest finite temperature. Tests allow four ULPs from that independent reference. This is numerical evidence over the stated grid, not a universal correctly-rounded guarantee.

Three native GB numerical policies deliberately extend or correct source evaluation:

- A nonfinite ordinary evaluation retries an algebraically equivalent energy-domain sum. This retains finite low-temperature results without exponential overflow; tied maxima contribute their entropy term.
- If positive `R*T` rounds to zero, the maximum-energy limit is returned. Under the input bound, the omitted entropy contribution is far below a result ULP.
- Empty input uses its exact single-site identity at every finite positive temperature. Source direct arithmetic otherwise loses energy at enormous temperatures even when its result remains finite.

These are documented numerical differences, not assertions of bitwise parity. The source's tests at 100 K assert only inequality from default-temperature literals and do not establish finite output values.

Hydrophobic moments are independently checked using the normalized Eisenberg fixture, including 100°, 160°, zero and negative angles. Window phase restarts at zero; multiplication precedes division in degree conversion; square/add/square-root behavior is retained. The seven normalized/unnormalized scales remain distinct. Annotation checks use mass-only sequences with unavailable composition to prove that property routines do not accidentally require empirical formulas. Hydrophobicity and GB ignore all annotations; pI suppresses any separately attached terminal annotation while leaving parent sidechain pKa unchanged. This does not model PTM-specific charge or hydrophobicity changes.

## Review scope and limitations

Read-only production review checks table mapping, source split/rolling/trigonometric addition order, terminal and sidechain selection, bisection endpoint priority, exact-zero midpoint handling, and bounded work before allocation or repeated traversal. Numeric temperature/pH/angle checks and errors for exhausted pI work or stalled bisection are native checked boundaries. The underlying scales and formulas are the pinned source's empirical models; matching these references is not independent validation of their biological applicability.

This increment covers these property utilities. It does not claim a complete Rust replacement for the full OpenMS core, new property scales, pKa predictions from modification chemistry, or universal numerical parity across C++ math libraries.

## Focused validation

All seven reference tests pass on the current compiler with all features and on Rust 1.85 with default features disabled. Strict Clippy passes for the reference target in both configurations. The fixture audit verifies sixteen source hashes, five fixture hashes, all 633 rows and their decimal/bit pairs, unchanged copies of the four source-extracted TSVs, exact native/source correspondence for all 62 private GB values, and the document's relative links. No source-derived expected value was changed to match production. The read-only implementation review found no actionable source-order, numerical-policy or resource-limit discrepancy in this scope.
