# PEFF (PSI Extended FASTA Format)

`format::peff` ports `FORMAT/PEFFFile.h` / `PEFFFile.cpp` from OpenMS4-core
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`: the eight record types, the reader
and writer, `isPEFFFile`, `toProForma`, and the `PEFFEntry` sequence members
that materialise, digest and enumerate proteoforms. It is always compiled (no
feature): the reader and writer are plain text, and the sequence members use
only `chemistry`.

| Rust file | Covers |
| --- | --- |
| `src/format/peff.rs` | records, `PEFFFile`, `read`/`write` over streams |
| `src/format/peff/read.rs` | `readStart` header parse, `readNext`/`readEntry_`, `parseAnnotations_`, tuple parsers, `isPEFFFile` |
| `src/format/peff/write.rs` | `formatHeader_`, database blocks, `formatEntry_`, `toProForma` |
| `src/format/peff/text.rs` | PEFF escapes, annotation identifiers, `CustomKeyDef` field-type validators |
| `src/format/peff/sequences.rs` | `getSequence` … `generatePeptides`, `enumeratePEFFModifications_` |

## API mapping

| Source member | Rust |
| --- | --- |
| `PEFFModification` (fields, both constructors, `operator==`, `Type`) | `PEFFModification` public fields, `new`, `with_details`, `PartialEq`; `PEFFModificationType::{PsiMod, Unimod, Generic}`, `of_accession`; field `type` is `modification_type` |
| `PEFFVariantSimple` | `PEFFVariantSimple`, `new`, `Default` (`variant_aa` `'\0'`) |
| `PEFFVariantComplex` | `PEFFVariantComplex`, `new` |
| `PEFFProcessedRegion` | `PEFFProcessedRegion`, `new` |
| `PEFFDisulfideBond` | `PEFFDisulfideBond`, `new` |
| `numeric_limits<UInt>::max()` sentinel | `ANNOTATION_ID_NOT_SET` |
| `PEFFEntry` fields, copy/move, `operator==` | `PEFFEntry` public fields, `Clone`, `PartialEq`; `std::map` → `BTreeMap` |
| `PEFFEntry::toFASTAEntry` / `fromFASTAEntry` | `to_fasta_entry` / `from_fasta_entry` |
| `PEFFEntry::getSequence` | `aa_sequence` |
| `PEFFEntry::getModifiedSequence` | `modified_sequence`, `modified_sequence_with_registry` (also returns the warnings) |
| `PEFFEntry::getVariantSequences(desc, seqs, include_complex)` | `variant_sequences(include_complex) -> PEFFPeptides` |
| `PEFFEntry::getProcessedSequence(region = "PEFF:0001021")` | `processed_sequence(region)`; default as `SIGNAL_PEPTIDE` |
| `PEFFEntry::digestWithVariants(digestor, desc, seqs, min, max, ref, var, mod)` | `digest_with_variants(digestor, &VariantDigestOptions)`, `digest_with_variants_with_registry` |
| `PEFFEntry::generatePeptides(digestor, desc, seqs, fixed, variable, max_var, min, max, ref, var, mod, cap)` | `generate_peptides(digestor, &PeptideGenerationOptions)`, `generate_peptides_with_registry` |
| private static `enumeratePEFFModifications_` | public associated `PEFFEntry::enumerate_peff_modifications` (returns the warning too) |
| `MAX_COMBINATORIAL_ELEMENTS`, `MAX_MODS` (function-local, both 20) | `MAX_COMBINATORIAL_ELEMENTS` |
| `PEFFCustomKeyDef` | `PEFFCustomKeyDef` |
| `PEFFDatabaseMetadata`, `SequenceType::{AA, NA}` | `PEFFDatabaseMetadata` (`Default` version `1.0`), `PEFFSequenceType::{Aa, Na}` |
| `PEFFFile()` / `~PEFFFile()` | `PEFFFile::new` / `Drop` (closes the output as `writeEnd`) |
| inherited `ProgressLogger` | public field `PEFFFile::progress` |
| `load(filename, entries, headers)` | `load(path) -> PEFFDocument` (entries, headers, warnings) |
| `store(filename, entries, header)` | `store(path, entries, header) -> warnings` |
| `store(filename, entries, headers)` | `store_multi` |
| `readStart`, `readNext`, `getHeaders`, `atEnd` | `read_start`, `read_next`, `headers`, `at_end` |
| `writeStart(filename, header)` / `(filename, headers)` | `write_start` / `write_start_multi` |
| `writeNext`, `writeEnd` | `write_next`, `write_end` |
| static `isPEFFFile` | `PEFFFile::is_peff_file` |
| static `toProForma` | `PEFFFile::to_proforma` |
| private `parseHeaderLine_`, `parseAnnotations_`, `parse*_`, `parseParenList_`, `formatHeader_`, `formatEntry_`, `readEntry_` | private functions in `read.rs` / `write.rs`, each citing its source |
| — (native) | `PEFFReader<R: BufRead>`, `read`, `read_with_limits`, `write`, `ReadLimits`, `PEFFFile::warnings`, `PEFFPeptides` |

Every public member of the header is mapped; nothing is left unported.

## Preserved source conventions

- Reading follows `readStart` / `readEntry_` byte for byte: `#` header lines,
  the file description block carried into every database block, the
  required-key warnings (`PEFFFile.cpp:1493-1527`), `>` entries whose sequence
  ends at a line end directly followed by `>`, carriage returns dropped, and
  the two `ParseError` messages of `readNext`.
- `parseAnnotations_`: the `\Key=value` separator is an unescaped backslash
  followed by a letter, digit or underscore; `\PName` in list form keeps a
  balanced `(…)` whole; `\ModRes*` comma positions expand into one
  modification each; the PEFF escapes `\\`, `\|`, `\(`, `\)` are reversed once,
  at the field that owns the text, so a tuple's escaped `|` never splits it.
- `CustomKeyDef` validation (field types `xsd:decimal`, `xsd:boolean`,
  `xsd:date`, `xsd:time`, …) produces the source's warnings, never errors.
- Writing emits the source's header order, 80-character sequence lines,
  `\ModResPsi` / `\ModResUnimod` / `\ModRes` by modification type, and the
  required-key warnings on write.
- `toProForma`: localised modifications as `[accession]` after the residue,
  unlocalised ones (position 0) as a `<?[…]>` prefix, `[name]` without an
  accession; disulfide bonds are not rendered.
- `getModifiedSequence` skips unknown and out-of-range positions silently and
  resolves each other modification by accession, then name, for its residue
  with anywhere specificity; an unresolved one is the warning
  `Could not resolve modification: <accession> (<name>)`.
- `getProcessedSequence`: for a signal peptide (accession `PEFF:0001021` *or*
  a region named `signal peptide`) the mature protein after the region, else
  the region itself; modifications are not applied; an empty sequence when no
  region yields one.
- `digestWithVariants` / `generatePeptides`: the digest uses the digestor's
  enzyme, specificity and missed cleavages with the caller's length window
  (`max_length` 0 = none), in the source's product order; stop-codon variants
  are never used; combinations with two variants on one residue, or with a
  modification that does not resolve on its (possibly substituted) residue,
  are dropped; more than 20 variants plus modifications in one peptide are cut,
  variants first, with the source's warning; descriptions `E4L(tag)`,
  `UNIMOD:35@4` joined by `;`, and `+varmod` for peptides a variable
  modification changed. `max_peff_mods_per_proteoform` k > 0 keeps the
  subsets of size at most k.

## Native differences

- **Out-parameters become return values.** `load` returns a `PEFFDocument`;
  the sequence members return `PEFFPeptides` (index-aligned `descriptions`
  and `sequences`, plus `warnings`). The source's many positional
  arguments of `digestWithVariants` / `generatePeptides` are the fields of
  `VariantDigestOptions` / `PeptideGenerationOptions`, whose `Default`s are
  the source's default arguments.
- **Warnings are returned, not logged.** Every `OPENMS_LOG_WARN` line of the
  reader, writer and sequence members is returned verbatim, in order.
- **Aligned outputs.** `getVariantSequences` and the reference branch of
  `digestWithVariants` push the description *before* `AASequence::fromString`;
  when that throws (for example a lowercase variant residue) the source's two
  vectors fall out of step. Here a sequence that does not parse drops its
  description too. See the defect candidates in the lane report.
- **Registries.** The source resolves against the process-wide
  `ModificationsDB`, which loads Unimod *and* PSI-MOD. `ModificationsDB::global()`
  holds the pinned Unimod only, so a `MOD:` accession resolves only through a
  caller registry extended with PSI-MOD (`ModificationsDB::extend_obo`, which
  registers PSI-MOD's Unimod aliases as the source does); every sequence
  member has a `_with_registry` form for this. With several matches the source
  warns and takes the first; this takes the first silently.
- **Unspecific cleavage order.** `ProteaseDigestion::digest_ranges` enumerates
  unspecific products by start; the source's `ProteaseDigestion::digest`
  (tokenize at every residue, then by missed-cleavage level) yields them by
  length, then start, so the ranges are re-sorted to that order.
  `DigestionSpecificity::None` is refused as the source's `InvalidValue`.
- **Output ceiling.** `max_peptides` (default 2^24) bounds the generated
  peptides; exceeding it is `Error::InvalidValue`, never a truncation. The
  source has no bound beyond its 2^20 combinations per peptide.
- **Resource bounds on reading** (`ReadLimits`): line, sequence, header-line,
  items-per-entry, entry and warning counts. Warnings past `max_warnings` are
  counted in `PEFFReader::suppressed_warnings`.
- **UTF-8.** Fields must be UTF-8 (`Error::Parse` otherwise); the source keeps
  bytes. A variant substitution that would split a multi-byte character is
  skipped (such text cannot parse as a sequence in either implementation).

## Checked boundaries and evidence

**Tier 1 (executed differential).** `tests/data/peff_sequences_oracle.tsv` is
the retained stdout of `../oracle/F6-peff-sequences/driver.cpp`, built against
the Linux x86_64 OpenMS4 **Release** install `bc9cc12-c19e494-174b576`
(conda-forge g++ 14.4, `-O2`) and run twice on **ibminode06**, byte-identical.
The `oracle_*` tests of `tests/peff_sequences.rs` compare every (description,
peptide) pair, in order, of 23 cases: `getVariantSequences` (including the
desynchronised vectors, which the port deliberately does not reproduce),
`digestWithVariants` under trypsin, chymotrypsin, two missed cleavages, semi
specificity, unspecific cleavage and no cleavage, with variants, conflicting
variants, modifications on substituted residues and name-only modifications,
and `generatePeptides` with fixed and variable modifications, missed
cleavages, no reference and the proteoform cap; plus `getModifiedSequence`
(PSI-MOD alias, and the two warnings of the retained stderr) and
`getProcessedSequence` (the `signal peptide` name rule, a region ending at the
sequence end, an absent region).

**Tier 3 (source review).** `tests/peff.rs` (32 tests) and the first 18 tests of
`tests/peff_sequences.rs` transcribe the sections of `PEFFFile_test.cpp` at
`bc9cc12`, over its eleven byte-identical `PEFFFile_*.peff` fixtures; the
proteoform tables (`REF_TABLE`, `K16R_TABLE`, `S14A_TABLE`, `CV_TABLE`, 52 rows
with an exact ProForma string and an exact `AASequence::toString` each) are
transcribed row for row, and the proteoform-cap counts 10 / 46 / 130 / 512 are
the section's binomial sums. The PSI-MOD assertions use a registry extended
with one independently authored OBO term (`MOD:00046`, origin S, linked to
Unimod 21).

**Tier 4.** The `ReadLimits` refusals, the header warnings of the minimal
fixture (derived from `PEFFFile.cpp:1493-1527`) and the unresolved-modification
warning without PSI-MOD (derived from `PEFFFile.cpp:586`).

The reader and writer have no executed oracle; their evidence is tier 3. See
`tests/data/peff_provenance.json`.
