# ProForma spectrum generation

This implements all six spectrum operations in the [pinned ProForma header](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/include/OpenMS/CHEMISTRY/ProForma.h#L1047). Both `Peptidoform` and `PeptidoformIon` expose `can_generate_spectrum`, `spectrum_generation_issues`, and `generate_spectrum`. Each takes a caller-owned mutable `ModificationsDB` and returns a `ConversionEvaluation<T>` containing the result and ordered, owned warnings. This group adds no parser, resolver, conversion, mass, or generic spectrum-generator substitute: it composes the existing implementations.

`SpectrumGenerationOptions` replaces the five source configuration arguments with `min_charge = 1`, `max_charge = 1`, `ion_types = "by"`, `add_losses = false`, and `add_metainfo = true`. The ordinary and linked generator settings are intentionally distinct. The source implementation and helper branches are at [ProForma.cpp 2016–2091](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2016) and [2745–2855](https://github.com/okohlbacher/OpenMS4-core/blob/82ce5b373c97f934ffd9b1ffd80215ca66473d0b/src/openms/source/CHEMISTRY/ProForma.cpp#L2745).

```rust
use openms::chemistry::{ModificationsDB, proforma::{Peptidoform, SpectrumGenerationOptions}};
let peptide = Peptidoform::parse("PEM[UNIMOD:35]PTIDE")?;
let mut registry = ModificationsDB::global().clone();
let report = peptide.generate_spectrum(&SpectrumGenerationOptions::default(), &mut registry)?;
assert!(!report.value.is_empty());
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Advisory checks and traversal

Chain advisory checks call the existing AASequence representability collector, which resolves a local sequence copy. If the first pass reports problems, source spectrum collection calls the conversion issue collector again. Newly interned formula definitions are visible to that second pass, so the returned issues can differ. A single-chain ion delegates to the chain collector; generation then delegates again to the chain generator, preserving the extra pass. Warning duplication and order follow these passes.

Ion checks return the source issue text and position zero, in this order: no chains, chimeric input, single-chain delegation, unsupported chain count, missing crosslink label, differing labels. Exactly two linked chains receive only the shallow label checks at this stage. These predicates do not promise successful strict conversion, valid generator charge settings, or available residue chemistry. For example, INFO-only annotation can pass the predicate and fail strict attachment (inherited [CPP-029](../OpenMS_CPP_ISSUES.md)); bare X can pass conversion checks and fail ordinary mass generation.

The first crosslink search visits ordinary `SequenceElement` entries only. Within each residue it checks the first alternative of each modification and stops at the first crosslink label. Later alternatives, later labels, terminal labels, modified ranges and ambiguous subsequences do not participate. Raw `MassDelta.mass` takes precedence over an existing resolved handle's difference mass; otherwise an existing handle supplies the mass, or it is zero. Searches happen on the original input before the two local conversions. Label scores are ignored. Matching uses label identifiers. Chain/ion names, charges and adducts are not consumed by this source operation and are not copied or charged by its private sequence-only session.

## Scientific output and preserved source defects

Single chains use `FailOnLoss` conversion and the existing [ordinary generator](THEORETICAL_SPECTRA.md). Literal membership of `a,b,c,x,y,z,M,I` selects the six fragment ladders, precursor peaks and immonium ions. Case matters; unknown characters and repeated flags have no further effect. Source defaults remain: no first prefix ion, MS level 2, centroid type, automatically inferred precursor charge, and ordinary annotation names `IonNames` / `Charges`. Losses and annotation flags are forwarded. Empty direct chain ASTs return the backend's empty result; parsing an empty string remains a separate syntax error.

Two-chain generation converts alpha then beta using `BestEffort`, constructs an owned `ProteinProteinCrossLink`, calls the [XLMS backend](THEORETICAL_XLMS_SUPPORT.md) for alpha then beta, and performs a final stable position sort. It sets only the six ladder flags, precursor flag, losses and metainfo. `I` is ignored here; default K-linked ions and charge annotations remain enabled even when metainfo is disabled. Isotopes remain disabled. Duplicate precursor peaks from the two calls are retained. It does not substitute ordinary generator metadata for XLMS metadata.

Two verified source defects are preserved with separate chemical expectations in the tests:

- **CPP-038:** converted endpoints retain attached linker chemistry, and a separate linker mass is added again. One attached endpoint with difference mass D therefore contributes 2D in the linked precursor calculation; two attached endpoints can contribute 3D. At charge two, a one-endpoint D=100 example differs by 50 m/z from the one-linker chemical expectation. This wrapper does not compensate for the defect.
- **CPP-047:** the crosslink position counter advances only over ordinary sequence elements, although conversion flattens ranges and ambiguous sequences. A range AG before a linked M yields source index zero, rather than flattened index two. The native wrapper uses that finite source index and retains checked backend rejection if an index is out of range.

The separate linker mass uses alpha only when `alpha_mass > 0.001`; otherwise it uses beta. This is a strict signed threshold, not an absolute-value test. The pre-resolution timing can make a pre-resolved NamedMod produce a different separate linker mass from the same unresolved original AST. These behaviors are explicit source compatibility choices, not mathematical or chemical corrections. See [the C++ issue log](../OpenMS_CPP_ISSUES.md) and [conversion support](PROFORMA_CONVERSION_SUPPORT.md) for inherited representation, ambiguity and annotation boundaries.

## Transactions and limits

The input AST is immutable. One resolver session, warning stream, item allowance, work allowance and logical byte ledger cover all issue passes, both conversions, both generator calls and final publication. Formula records are staged only when insertion is needed; existing records remain shared through `Arc`. Completed advisory results can publish staged definitions even if issues are nonempty. Successful spectrum calls publish the registry and return the full spectrum together. Any checked conversion, numerical, resource or backend error leaves the caller's registry unchanged; no partial spectrum is returned.

Shared defaults are the existing ProForma limits: 50 million work units, 256 MiB of cumulative logical allocation estimates, one million traversed items and 4 MiB per consumed text field. AST traversal follows its fixed-depth container structure. These are conservative accounting limits, not physical RSS ceilings. The ordinary adapter precharges all wrapper-reachable chemistry, duplicate loss declarations and comparisons, slice/formula reconstruction, output and sorting before calling the existing generator. It deliberately uses full-peptide bounds even for shorter fragments and possible rows even when later filters remove them. This can reject a large input earlier than actual-work accounting would. Each ordinary generator's existing category limits still apply, including 4096 residues and 100,000 output peaks. XLMS preserves its 4096-residue, 100,000-peak, 50-million-work and 64-MiB logical per-call limits while also consuming the outer remaining counters across both sides.

The existing ordinary API represents fragment charges as `u8`; this wrapper accepts 1..=255 there and checks conversion from its signed options. XLMS retains finite source behavior for negative charges, inverted ranges and signed output; a consumed division by zero or nonfinite result is rejected. Out-of-range positions, count overflow, unrepresentable chemistry and allocation failure are checked errors. Unicode is preserved in consumed text; sequence symbols retain the existing conversion/AASequence boundaries. Standalone ordinary and XLMS public behavior and numerical code remain unchanged by the private shared-budget adapters.

## Evidence

[Direct tests](../tests/proforma_spectra.rs) reproduce the source spectrum class-test assertions and add independent branch/transaction cases. They cover all flags, ordinary/backend output equivalence, exact issue text/order, ignored context, pass-sensitive interning, alpha/beta warning ordering, resolved-handle timing, threshold endpoints, source defects 038/047, signed XLMS charge cases and registry rollback. Backend comparisons test composition/forwarding and are not independent C++ spectrum oracles. The two CPP-038 numerical checks instead use explicit free-residue formulas, water and proton constants and distinguish source output from the one-linker chemical expectation.

[The source projection](../tests/data/proforma_spectra_source.json) preserves ten complete sections and 18 assertions from `ProFormaParser_test.cpp:2998–3135`, including six generation calls, 15 parse inputs and 11 explicit resolution calls. Its six spectrum assertions retain the source's nonempty, minimum-count or comparative strength. It contains no generated peak masses or stronger invented count oracle. Reproduce it with:

```text
python3 tools/generate_proforma_spectra_reference.py /path/to/exact/SDK --check
```

[Provenance](../tests/data/proforma_spectra_provenance.json) lists the exact source hashes, extraction script and fixture. The source extraction reads text only; no C++ compiler or full SDK execution is claimed for this group. Existing ordinary and XLMS numerical fixtures retain their original separate provenance. Focused validation and exact extraction details are recorded in the handoff accompanying this group.
