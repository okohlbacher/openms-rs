# Modification generation reference review

This independent review covers `ModifiedPeptideGenerator`,
`ModificationDefinition`, and `ModificationDefinitionsSet` at OpenMS4-core
revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. It checks actual source
branches and literal upstream assertions; no C++ code was built or executed.
It does not claim complete OpenMS chemistry or OBO registry coverage.

The [provenance manifest](../tests/data/modification_generation_provenance.json)
records SHA-256 hashes and pinned links for 18 source files and four fixtures.
The primary sources are
[`ModifiedPeptideGenerator.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/ModifiedPeptideGenerator.cpp),
[`ModificationDefinitionsSet.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/ModificationDefinitionsSet.cpp),
their headers and class tests, and the resolved-annotation paths in
`AASequence.cpp`, `Residue.cpp`, and `ResidueDB.cpp`.

## Literal and derived evidence

The fixtures reconstruct seven literal fixed-placement outputs, nineteen
variable-generation configurations, and seven compatibility assertions. The
variable cases include empty inputs, maximum zero/one/general paths, existing
annotations, rightmost-methionine-first order, two modification alternatives
per cysteine, and the source's 26-, 71-, four-, and 199-result cases.
[`modification_generation_expected.tsv`](../tests/data/modification_generation_expected.tsv)
holds literal sequence expectations where the source gives them. Cases with
source pointer-dependent order compare sorted identities or counts, as the
source test does; the fixture marks each comparison mode explicitly.

Additional tests are independently derived from the source statements, rather
than presented as measured C++ outputs. An eight-result `CC` example checks
every variant identity and the complete weighted-site order. Two alternatives
per site produce two right-site variants, two left-site variants and four
two-site variants. Within the last group the right site's alternative varies
fastest. This detects duplicate, missing, or misplaced variants that a matching
count alone would miss.

The independent suites are
[`modified_peptides_reference.rs`](../tests/modified_peptides_reference.rs) and
[`modification_definitions_reference.rs`](../tests/modification_definitions_reference.rs).
Implementation-owner tests and the separate
[workflow suite](../tests/modification_generation_workflow.rs) complement these
references; they are not used to manufacture expected fixture outputs.

## Terminal and chemical identity

Several unusual outcomes are finite, defined source behavior and are retained:

- Fixed placement has a terminal prepass that ignores origin, followed by a
  residue pass that may overwrite an existing terminal annotation. Already
  modified residue slots remain protected.
- Maximum-one variable generation visits residues from right to left and puts
  a matching residue-specific terminal record on the **residue slot**. Generic
  terminal records do not match ordinary residue letters in this path.
- General generation inserts vacant terminal alternatives without checking
  origin, then inserts matching boundary alternatives again. A matching
  residue-specific terminal record can therefore produce two identical
  variants. The second pass can also replace an occupied terminal slot.
- Protein N/C-terminal records are ignored by placement. They are not silently
  treated as peptide-terminal records.
- An empty sequence can receive a terminal record in the fixed/general paths.
  The source bypasses formula and mass accumulation for empty sequences; its
  empty formula and zero mass are retained alongside the typed annotation.

For `QAA` and `Gln->pyro-Glu (N-term Q)`, maximum one yields one residue-annotated
variant; a maximum greater than one yields two identical N-terminal variants.
The empirical formulas agree, but their monoisotopic masses need not be
identical. `Residue::setModification` calls `setFormula`, recomputing the
residue's mass from composition, while a true terminal annotation contributes
its declared rounded mass difference. Here the two deltas are approximately
−17.0265490957 and −17.026549 Da. Independent b/y-ion checks retain that
distinction: b ions receive the appropriate delta per charge, while y ions
excluding the Q remain unchanged.

Some such typed states cannot be reconstructed by the ordinary `AASequence`
parser. The generator preserves their positions; it does not weaken the
parser's positional checks or invent a serialization syntax. The idXML writer
preflights exact typed identity after parsing the displayed sequence and rejects
unrepresentable states before output. The independent empty-annotated-sequence
case complements the workflow's terminal-on-residue and foreign-origin terminal
cases. Ordinary generated variants retain lossless idXML interchange.

The native implementation replaces unspecified pointer/unordered alternative
ordering with full-ID order. This also makes competing fixed assignments
deterministic. Source undefined behavior, such as the null cached residue for a
generic terminal record applied to bare `X` in the maximum-one path, is a
checked error. Chemistry that produces negative atom counts or invalid masses
also fails atomically instead of bypassing existing chemical invariants.

## Definitions, inference and absolute masses

The source compatibility predicate does not enforce stored occurrence or total
modification maxima. It requires every residue matching a fixed definition's
origin to have the same short modification name, irrespective of terminal
specificity, and requires each present annotation's full ID to be allowed.
Consequently, generic fixed N-terminal Acetyl accepts ordinary unmodified `AAA`,
whereas an inferred fixed N-terminal pyro-Q can reject the same peptide carrying
that modification solely in its terminal slot. Both cases are explicitly tested.

Inference reads every hit, pools residue letters and terminal slots separately,
and distinguishes an exclusively modified site category from a category that
also has an unmodified observation. Replacement is atomic. Anonymous mass tags
are retained as owned annotations in inferred definitions; they do not become
named immutable registry records. Their named-record accessor returns
`Unsupported`, while their mass-tag accessor retains the exact spelling.
Conflicting chemical records sharing a full ID are rejected rather than letting
hit order select the inferred chemistry.

Mass matching preserves origin filtering before residue-alias resolution,
inclusive absolute tolerance, exact terminal filters, and fixed-before-variable
order for equal errors. Tests cover negative deltas, the adjacent floating-point
tolerance below an exact inclusion boundary, duplicate IDs in both partitions,
and source inference/compatibility inconsistencies.

Absolute mode first uses a positive stored monoisotopic mass. Otherwise, a
nonempty residue selector invokes the source fallback: declared difference plus
the mass of the **full residue formula minus water**. This is not replaced with
an independently summed peptide mass. Without a residue selector, even a zero
or negative stored value is compared literally. Native B/Z/X fallback rejects
unknown chemistry instead of manufacturing a negative mass from source empty
formula placeholders; a positive stored mass bypasses that fallback.

The source provider distinction matters:
[`OBODataProvider.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/CHEMISTRY/OBODataProvider.cpp)
reads `MassMono` and `MassAvg` into absolute fields. UniMod/custom XML passes
through `UnimodXMLDataProvider` and `UnimodXMLHandler`, which populate difference
masses. `ResidueModification` initializes absolute fields to zero; its record
loader does not change that convention. The native checked owned-record builder
exercises positive stored masses, including lookup bypass for B/Z/X, without
changing the embedded TSV schema. The subsequent OBO and XLMOD provider work is
covered by the [registry reference review](MODIFICATION_REGISTRY_REFERENCE_REVIEW.md).

## Formula-free custom records

The public generator also accepts shared owned records from a custom native
registry. A separate source-derived regression exercises the actual
`Residue::setModification` branch rules on eight small custom records:

- With an empty delta formula and a nonzero mono or average difference, a
  nonzero absolute mono value is a **free-residue** mass. The internal mass is
  that value minus water, independent of whether the original residue's mass
  is known. The test includes an absolute override on `X`.
- Without an absolute override, the operation is full-residue formula mass plus
  declared mono difference, then water subtraction. A binary64 assertion pins
  that order independently of a separately summed peptide mass.
- Absolute-only records with zero differences and no formula leave residue
  chemistry unchanged, despite retaining their typed annotation.
- A nonempty delta formula takes precedence over declared absolute and delta
  masses. True terminal slots continue to use the declared mono difference.
- An average-only difference still indicates missing composition, even though
  the calculated mono mass remains unchanged. Invalid negative internal mass
  returns an error before fixed placement changes the input peptide.

For a mass-bearing change without a formula, the native sequence exposes its
known mono mass but reports formula and average mass as unavailable. This is an
explicit correction to the source's retention of the old formula after a
formula-free mass change; that retained composition would misrepresent the new
residue. Fragment checks preserve the free-to-internal water subtraction and
leave complementary ions outside the modified residue unchanged. The separate
workflow checks theoretical monoisotopic ions, checked isotope rejection and
custom-record interchange boundaries. No placeholder composition is invented.

## Resource and mutation review

The generator counts weighted subsets and alternative multiplicities before
cloning combinatorial peptide outputs. Existing append contents count toward
output count and payload limits, even when the request generates no new variants.
The independent tests verify unchanged append vectors on both count and byte
failures. Production review also checked overflow arithmetic, site/work charging,
bounded subset planning, and final staged commit after chemistry validation.
The payload allowance is a documented conservative logical bound, not a process
memory measurement.

Definitions use checked work budgets and preflight replacement input before
cloning. Count settings remain source metadata rather than being repurposed as
resource limits. The public views cannot mutate map identity. Details and native
boundaries are documented in
[modified-peptide support](MODIFIED_PEPTIDES_SUPPORT.md) and
[definition support](MODIFICATION_DEFINITIONS_SUPPORT.md).

Focused validation runs both reference suites with all features on the current
toolchain and with/without optional formats on Rust 1.85. The all-feature suites
contain fourteen tests; thirteen run without the idXML feature. Strict Clippy,
formatting, fixture/source hashes and local documentation links are checked
separately. No new dependency, vendored source edit, or C++ execution is involved.
