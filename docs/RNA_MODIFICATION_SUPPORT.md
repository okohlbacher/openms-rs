# RNA modification generation

`chemistry::ModifiedNASequenceGenerator` implements both public operations of the pinned C++ `ModifiedNASequenceGenerator`: fixed replacement and combinatorial variable placement. It uses the native owned RNA records and sequences described in [RNA_SUPPORT.md](RNA_SUPPORT.md). No additional feature or dependency is required.

```rust
use openms::chemistry::{ModifiedNASequenceGenerator, NASequence, RibonucleotideDB};

fn main() -> openms::Result<()> {
    let database = RibonucleotideDB::global();
    let sequence: NASequence = "AUAUAUA".parse()?;
    let modifications = vec![database.get("m3U")?, database.get("s4U")?];
    let generator = ModifiedNASequenceGenerator::default();
    let variants = generator.variable_modifications(&modifications, &sequence, 3, true)?;
    assert_eq!(variants.len(), 27);
    Ok(())
}
```

`apply_fixed_modifications(&modifications, &mut sequence)` replaces one sequence atomically. `variable_modifications(&modifications, &sequence, maximum, keep_original)` returns a new vector. `apply_variable_modifications(&modifications, &sequence, maximum, &mut output, keep_original)` appends atomically, retaining existing entries first. All three return `Result`. Modification inputs are slices of `Arc<Ribonucleotide>`, so caller-created records work without registering names or keeping a registry alive.

## Compatibility and source behavior

The generator installs complete records; it does not add mass differences or combine modification formulas. It accepts the finite signed stored fields, empty formulas and explicit placement states already accepted by the RNA record/sequence APIs. Sequence formulas and masses remain separate operations. Ambiguous record codes are not expanded through the registry's alternatives map.

An original residue is protected when `is_modified()` is true: its code has byte length other than one, or its one-byte code differs from its own origin. Otherwise a candidate matches when its origin equals that one-byte code. `X` and `.` are literal origins, not wildcards.

Fixed modification placement has two effects:

- An empty five-prime or three-prime slot receives the **first** corresponding terminal candidate. The candidate's origin is ignored, including for empty sequences. An occupied end remains unchanged.
- An eligible ordinary residue receives the **last** matching Anywhere candidate. Compatibility uses the original residue throughout; an earlier replacement does not prevent a later matching candidate from winning.

Variable generation counts new placements only. Already-modified residues and occupied terminal slots do not reduce `maximum`. With no candidates, maximum zero, or no compatible sites, the only possible addition is the original sequence when requested.

The source's maximum-one shortcut differs from its general path. With `maximum == 1`, it visits ordinary residues right to left, checks origin, and **does not check the candidate's terminal specificity**. Thus a terminal-specific record can replace an ordinary residue. This shortcut never places a record in an end slot; an already occupied end does not prevent the unusual residue replacement. Empty input yields only the optional original.

With `maximum >= 2`, ordinary sites accept only Anywhere records. Each empty terminus is a separate site accepting records of its own specificity, independent of origin. Empty input can consequently yield terminal-only variants. No duplicate terminal compatibility entry is added by a second residue pass, unlike certain peptide-generator cases.

## Order, identity and serialization

C++ inputs are pointer-ordered sets. Native inputs retain **caller slice order**, removing only repeated Arc allocation identities. Clones of the same handle count once; separately allocated equal-value records remain separate alternatives. Records with identical codes but different chemistry also remain separate. Neither candidates nor output are deduplicated by code, formula, displayed text or value. Passing a canonical base record as its own modification can produce repeated unchanged values.

The optional original begins the appended block. Variants then appear by increasing placement count. Ascending site order is three-prime, five-prime, then ordinary residues from left to right; subsets are visited in reverse lexicographic order. For four ascending sites a,b,c,d, two-site order is **cd, bd, bc, ad, ac, ab**. Within each subset the first selected site's alternatives form the outer loop, and the last site's alternatives change fastest. Single-placement order is rightmost ordinary residue first, followed by five-prime and three-prime sites when the general path includes them.

Generated sequences retain the input's owned sulfur-slicing context. Introducing a `*` record into a sequence with no captured `5'-p*` context can cause a later sulfur-boundary slice to fail; generation does not silently substitute a global record.

Use `NASequence::checked_string_with_registry` before exchanging text. A terminal-specific record installed in a residue usually reparses into an end slot, and a custom same-code record may resolve to different chemistry. Those states are usable natively but cannot always round-trip through display text. Terminal-only empty sequences have the source's empty-formula behavior. Generated record values and lifetimes remain independent of display.

## Bounds and errors

Default configurable bounds are:

| Field | Default | Scope |
|---|---:|---|
| `max_residues` | 10,000 | Input sequence residues |
| `max_sites` | 1,024 | Supplied candidates before identity deduplication, and total retained compatible site/record entries |
| `max_work` | 10,000,000 | Cumulative scan, comparison, planning, validation and copy allowance |
| `max_outputs` | 100,000 | Existing plus newly emitted entries; fixed replacement counts as one output |
| `max_output_bytes` | 256,000,000 | Conservative logical payload of existing and newly emitted sequences |

Residue, site and work limits must be positive. Zero output limits are valid, including a successful append that has neither existing nor new entries. Existing output count and payload are checked even on no-op requests. Existing entries need not obey the new input's residue limit, but their traversal and payload consume the call's budgets.

Pointer deduplication uses a simple bounded quadratic scan; its complete comparison allowance is charged first. Compatibility storage is bounded, and a weighted-subset dynamic program checks counts and byte estimates before any combinatorial sequence cloning. Checked arithmetic rejects overflow. Iterative subset and alternative counters avoid recursive stack growth or a materialized powerset.

Payload includes sequence/container handles, logically referenced record strings and formula nodes, both ends, and the private sulfur context. New variants conservatively retain the input's payload estimate and add installed records without subtracting replaced ones. Shared records may therefore be counted many times. The estimate is intentionally not an allocator or resident-memory measurement. Planning allocations also have a separate cumulative 64 MiB ceiling; staged output and final vector growth can temporarily overlap, while remaining bounded by the checked output quantities. Allocation failures are returned where supported by the allocator.

Per-variant copy/validation work is precharged before staging. Existing sequence limits are still enforced on every final variant. A later failure—including a final record replacement that exceeds the sequence's text limit—leaves the input sequence or output entries unchanged. Caller-configured limits can reject otherwise finite C++ workloads; they do not silently truncate the generated set.

## Source and validation

The authoritative implementation is [ModifiedNASequenceGenerator.cpp at 7c029e8](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/ModifiedNASequenceGenerator.cpp); its maximum-one specificity behavior overrides the broader prose in the header. SHA-256: `a838584c26985f75c50b1b1922d24cc8ade0379cde17f4473d44954acc10ed7e`.

[The pinned class test](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/ModifiedNASequenceGenerator_test.cpp) supplies fixed `AUAUAUA → A[s4U]A[s4U]A[s4U]A` and variable counts 7, 6, 27 and 432, plus six literal single-placement strings. Test SHA-256: `fc9fd9eaf392d09f721c0f65ad0633bcb3a80ec25dcf2d3157738ec8b80b4ef6`.

[Direct native tests](../tests/modified_na_sequence_generator.rs) preserve those assertions and add independent subset order, typed terminal placement, caller ordering, Arc identity, signed records, ambiguity/no-op multiplicity, sulfur context, overflow, cumulative limits and late-failure atomicity checks. No C++ build or execution is used for these references.
