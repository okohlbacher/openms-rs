# Modified peptide generation

`chemistry::ModifiedPeptideGenerator` applies fixed modifications and enumerates
variable variants using shared immutable modification records. Resolve a name
list once with `get_modifications`, or `get_modifications_with_registry` for a
caller-owned database, then reuse its `Arc<ResidueModification>` handles across
peptides and threads. Neither the registry nor the handle list must outlive the
generated peptides; the last owner releases the records. Prefer full specificity IDs such as `Oxidation (M)`;
ambiguous or unknown names return errors. Repeated full IDs are collapsed, and
records are ordered by full ID. The source uses pointer/unordered-map order,
which is unspecified across executions; the native order is deterministic.

`apply_fixed_modifications` changes one peptide only after successful placement
and chemistry checks. `variable_modifications` returns a new vector;
`apply_variable_modifications` appends to an existing vector atomically. Both
variable methods take a maximum number of newly placed modifications and a
`keep_unmodified` flag. The flag preserves the supplied peptide, including any
fixed, named or anonymous annotations it already carries. Existing modifications
do not consume the new-placement maximum.

## Fixed placement

The source has two passes, both retained here:

1. Each empty peptide terminal slot receives the first corresponding N-terminal
   or C-terminal record, **without checking its residue origin**.
2. Residues are visited from left to right. A residue already carrying any
   annotation, including an anonymous mass tag, is skipped. Matching Anywhere
   records replace that residue's modification. Matching terminal records on
   the first/last residue replace the corresponding terminal slot.

The second pass can overwrite an existing terminal annotation, including a mass
tag. It can also replace a terminal record installed in the first pass. Competing
Anywhere assignments on an initially unmodified residue are all visited; the
last record in native full-ID order wins. Competing first-pass terminal records
use the first record in that order, before any second-pass overwrite. Protein
N-terminal and protein C-terminal records are resolved but ignored by placement,
matching the pinned generator rather than the broader sequence parser.

## Variable placement and order

No supplied records or a maximum of zero produces only the optional unchanged
peptide. Otherwise the source uses distinct maximum-one and general algorithms.

For a maximum of **one**, unmodified residues are visited from right to left.
Origins must match exactly. Anywhere records and boundary-compatible
residue-specific terminal records are placed **on the residue slot**. Generic
terminal records have source origin `X`; they produce no variant on ordinary
residues in this path. A bare `X` matching such a generic terminal record would
dereference a null cached residue in the C++ implementation; the native port
returns a checked error instead.

For a maximum **greater than one**, empty terminal slots first receive all
corresponding peptide-terminal alternatives, regardless of residue origin. A
second pass over unmodified residues adds matching Anywhere alternatives and
boundary-compatible terminal alternatives. This second pass does not test
whether the terminal slot was already occupied. Thus it can introduce a terminal
overwrite, and it can add the **same terminal alternative twice** when both
passes match. Those duplicate output variants are retained; there is no final
deduplication. At most one alternative is selected for each actual site.

Sites are processed in reverse residue-index order, followed by N terminus and
C terminus. Subset groups follow the source's append order: for sites A then B,
the groups are unchanged, A, B, A+B. Within a group, alternatives for earlier
sites vary fastest. The unchanged peptide, when requested, is first. Alternative
ordering within each pass uses native full-ID order. This preserves the source
rightmost-methionine-first examples and the 26-, 71- and four-variant test cases.

## Typed states, chemistry and interchange

The unusual source placements above are represented directly as typed
`AASequence` annotations using an internal bulk application helper. Ordinary
public setters and parsing retain their existing positional validation. Inspect
`residue_modification`, `n_terminal_modification`, and `c_terminal_modification`
to distinguish a terminal-specific record placed on a residue from a record in
a terminal slot.

Some source-generated typed states cannot round-trip through ordinary sequence
text. For example, the maximum-one path can put a pyroglutamate terminal record
on a Q residue; ordinary parsing may resolve its displayed name into a terminal
slot. Similarly, the general path may place a residue-specific terminal record
at a terminus with a different residue. Display is descriptive for these states,
not a promise of lossless interchange. The idXML writer checks that displayed
sequence text reconstructs the exact typed sequence before writing, and rejects
such unrepresentable states before output. Normal generated variants round-trip.

Known formulas and masses are recomputed once per final variant. Negative atom
counts, negative calculated known residue/peptide masses, and numerical overflow
remain errors; a source pointer assignment does not bypass native chemistry
checks. Residue placements use the modification's formula-derived mass; terminal
slots use its declared mass difference. Thus the maximum-one and general
pyroglutamate variants can have identical formulas but masses differing by about
`9.57e-8` Da, matching the source's two arithmetic paths. Unresolved B/Z/X
chemistry and existing anonymous mass tags remain
representable without inventing formulas or average masses. Existing occupied
residue slots are protected in both paths.

Custom immutable registry records can also describe a mass change without a
delta formula. The source considers a residue changed when either declared mass
difference is nonzero or a delta formula has atoms. With atoms, formula-derived
mass retains precedence over declared absolute values. Without delta-formula atoms, an absolute formula takes precedence and replaces
the free-residue composition before removing water; it can resolve X/B/Z
chemistry. With neither formula, a nonzero declared absolute monoisotopic mass is interpreted as a **free-residue mass** and
water is subtracted to obtain the internal mass; otherwise the calculation is
`(unmodified full-residue formula mass + declared mono difference) - water`, in
that order. An absolute override can supply mass for an otherwise unknown X/B/Z
residue. A delta alone cannot resolve an unknown base mass.

Mass-bearing changes with neither formula make the sequence's composition and
average mass unavailable, and `SequenceModification::diff_formula` reports
`Unsupported`. Monoisotopic peptide/fragment processing remains available where
all residue masses are known. Terminal slots continue using the declared mono
difference, independently of any absolute value, and likewise lose known
composition when no formula describes their mass change. A record with zero
differences and no formula atoms is a source residue no-op even if it stores
absolute masses or an atom-empty formula charge; those absolute values, absolute formulas and
charge do not change the residue. Terminal formula handling retains the separate
source convention. Negative calculated masses and overflow still fail atomically.

An empty sequence has no residue sites. Its maximum-one path generates no new
variant. The source general and fixed paths can nevertheless attach terminal
records to an empty sequence. The port retains those typed records and the
source's empty formula and zero mass. This empty-state exception is also subject
to the interchange check above.

## Bounds and allocation

Defaults are 10,000 input residues, 1,024 compatible placement entries,
10,000,000 work units, 100,000 output peptides and 256,000,000 accounted output
bytes. Point/site/work limits must be positive; zero output count or byte limits
are useful for checking that no output is produced. `max_sites` limits both the
supplied record count and the number of compatible placement entries, including
duplicated terminal entries and successive fixed assignments. The input residue
limit applies to the peptide being processed. Existing appended records are
counted individually for work and output bytes.

A checked dynamic program counts weighted subsets before any combinatorial
sequence cloning. The complete output count and conservative payload allowance
are checked, including existing append contents even when no new variants are
requested. It then stores subset index lists and emits each final sequence with
one clone/rebuild, avoiding the source's repeated intermediate sequence copies.
There is no recursive traversal or unchecked binomial coefficient conversion.
Size arithmetic and work counters use checked integer operations.

The byte allowance includes the inline `AASequence` value, residue-string bytes,
allocated annotation slots, and **both owned strings** of every anonymous mass
tag. Formula storage is conservatively allowed for six standard-residue element
types plus every existing known modification's delta and absolute formula entries; an entry costs
`size_of((Atom, i32)) + 3 * size_of(usize)`. Each generated sequence additionally
allows the largest candidate formula-entry cost times the maximum selected
sites. Fixed placement adds all final assigned records' entry allowances.
Duplicate formula elements and overwritten anonymous strings are intentionally
overcounted. Known chemistry uses shared `Arc` handles; record strings/formulas
are not cloned per output and their shared allocations are not charged as
per-peptide payload. This is an explicit
logical payload convention, not an allocator-overhead or process-RSS bound; it
may reject output that a tighter accounting would accept.

Work units cover input/existing-output residue scans, a logarithmic record-sort
allowance, record/site compatibility checks, dynamic-program cells, formula-cost
scans, subset-list visits and copies, and a residue-plus-placement rebuild
allowance for all outputs before they are allocated. Variable generation remains
combinatorial, with work and output caps bounding enumeration. Staging temporarily
holds generated peptides alongside the existing output, while index groups use
storage proportional to their total selected-site entries. Any error leaves
the original peptide or append vector unchanged. The initial name-resolution
helper is separate from these generation limits; it does not clone peptides.

## Sources and validation

The implementation follows OpenMS4-core
[`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`](https://github.com/okohlbacher/OpenMS4-core/tree/7c029e8cdba6abab503708ecdd56f6ab55e38ce4),
`CHEMISTRY/ModifiedPeptideGenerator.cpp/.h`, its complete class test, and the
resolved-pointer, formula and mass paths of `AASequence.cpp`. Independent source
fixtures and hashes are recorded under `tests/data/modification_generation_*`.

Focused tests cover literal fixed examples, exact reverse-site order, 26 and 71
combinations, the four terminal/oxidation variants, competing assignments,
terminal multiplicity/overwrite, ignored protein termini, empty states,
unresolved and anonymous chemistry, atomic errors, checked overflow, and
existing-output limits. An independent suite checks the extracted source
fixtures and typed terminal cases. No C++ code was built or executed and no new
dependency was added.
