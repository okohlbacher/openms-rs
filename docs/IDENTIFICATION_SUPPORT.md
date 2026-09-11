# Peptide and protein identification records

`openms::identification` ports owned identification records from OpenMS4-core revision `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. It connects modified peptide sequences, scores, protein evidence and typed metadata to spectra and feature maps. These records support storing and examining search results. Separate [score handling](SCORING_SUPPORT.md), [filtering](ID_FILTER_SUPPORT.md), and [FDR](FDR_SUPPORT.md) modules operate on them. [Peptide indexing](PEPTIDE_INDEXING_SUPPORT.md), [basic protein inference](PROTEIN_INFERENCE_SUPPORT.md), [conflict resolution](ID_CONFLICT_SUPPORT.md), and [file-origin splitting](ID_RIPPER_SUPPORT.md) provide further native operations. A database search engine, probabilistic inference engines and the separate `IdentificationData` reference graph remain outside this implemented surface.

## Native API

| Source model | Rust coverage |
| --- | --- |
| `PeptideEvidence` | Protein accession, optional zero-based inclusive endpoints, typed flanking residue/terminal markers, validation, ordering and native hashing; [complete evidence API](PEPTIDE_EVIDENCE_SUPPORT.md) |
| `PeptideHit` | Modified `AASequence`, score, rank, charge, evidence, peak annotations, analysis results, typed metadata, sequence/charge identity and target/decoy category |
| `PeptideIdentification` | Run identifier, candidate hits, score direction/type/threshold, optional RT/m/z, stable sorting, best-hit and accession lookup, spectrum reference and experiment label |
| `ProteinHit` | Score/rank, accession, raw protein sequence, optional coverage percentage, observed modifications, typed metadata, description and target/decoy category |
| `ProteinGroup` | Probability, ordered accessions and float/integer/string quantitative arrays; source comparison order |
| `ProteinIdentification::SearchParameters` | Database/version/taxonomy, charge declaration, mass type, fixed/variable modifications, missed cleavages, absolute/ppm tolerances, enzyme name/regex/specificity, typed metadata and mergeability |
| `ProteinIdentification` | Run/search-engine information, score settings, hits/groups, run paths, optional serialized date, coverage and observed-modification calculation |

Public fields use native names and can be edited directly. Call `validate` after such edits. Checked numerical and mutating operations validate their inputs and compute results before committing, so an error leaves prior records intact. Finite signed scores and coordinates are allowed. Protein sequences retain ASCII residue letters including ambiguous codes, `*`, `-` and `.`; they do not require an empirical formula. Peptide sequences preserve B/Z/X and named or numeric annotations without requiring calculable chemistry. The constructor trims outer accession/sequence whitespace. Raw date strings are retained without parsing or timezone inference.

The new records use [typed metadata](METADATA_SUPPORT.md). Existing measured kernel fields still use string metadata, with explicit conversion helpers. Rank and significance threshold are dedicated fields, not aliases that synchronize with entries named `rank` or `significance_threshold`. Run paths are likewise owned fields without hidden metadata aliases.

## Evidence and coverage

`PeptideEvidence` endpoints are **inclusive**: `0..=2` covers three residues. Missing endpoints are `None`. Flanking markers map `X` to unknown, `[` to the protein N terminus and `]` to the protein C terminus; explicit uppercase residue characters cover other flanks. Evidence ordering follows accession, start, end, before and after; missing positions sort first, and flanks use source character order. The native validity helper accepts the legitimate one-residue mapping `0..=0`, which the source helper excludes.

`compute_coverage` uses all evidence supplied by all peptide candidates, regardless of run identifier, matching the pinned implementation. It groups by protein accession and unions inclusive intervals, so overlap, duplicate evidence and adjacent regions are counted once. Every protein hit requires a nonempty sequence, even when it has no matching evidence; such hits receive zero coverage. Evidence for other protein accessions is ignored by the coverage calculation. It does not check sequence identity or require evidence length to equal peptide length, matching the source's interval-based calculation.

The implementation stores intervals rather than a protein-length bit array. Memory follows evidence count. Unknown or reversed intervals and an endpoint outside the sequence are errors. In particular, `end == protein_length` is rejected: the C++ bound check permits it before filling through `end + 1`, an out-of-bounds write. No partial percentages are assigned if a later protein fails validation. Coverage is `100 * union_length / protein_length`; `None` represents unknown coverage.

`compute_modifications` maps a modified peptide's residue at index `i` to `evidence.start + i`, its N-terminal modification to `start`, and its C-terminal modification to `end`. Duplicate `(position, complete modification value)` observations collapse. Distinct chemical records at the same position survive even when their displayed names and full IDs coincide. Skip sets accept modification names or full IDs. `ProteinModification` owns a cloned `SequenceModification`, so anonymous mass tags survive after the peptide records are dropped; known annotations retain shared owned `Arc` records. Known positions are required; arithmetic overflow, modified residues beyond the evidence endpoint and mapped modifications beyond a known protein sequence are errors. Proteins with no collected modifications retain their existing entries, including when every observation is skipped, matching the source. This method does not infer missing evidence or verify peptide/protein sequence identity.

## Scores, grouping and search settings

Peptide and protein sorting is stable and follows `higher_score_better`. Equal scores retain candidate order; sorting does not assign ranks. `best_hit` is a native convenience that leaves stored order intact and returns the first best-scoring candidate. Sequence/charge identity includes all modifications and excludes score and evidence. Peak annotations sort by m/z, charge, annotation text and intensity. Protein groups sort by descending probability, then ascending accession count, then the lexical accession list. Quantitative arrays represent samples, so their lengths need not equal the accession count; group equality includes the arrays.

Peptide target/decoy metadata accepts `target`, `decoy` and `target+decoy`; only `decoy` makes `is_decoy` true. Proteins allow the first two categories. Unknown status removes the metadata key. Getters accept case variations and reject invalid values. Spectrum-reference and experiment-label getters use lenient typed-value stringification as in source helpers. Setting an empty experiment label explicitly clears it; the source instead ignores the empty assignment. Native `is_empty` tests the candidate collection, rather than the C++ method's combined checks of identifier and score settings.

Charge declarations accept one signed integer, comma-separated signed charges, or two endpoints separated by `:` or `-`. A sign may precede or follow a charge. Lists return their minimum and maximum; empty declarations return `None`. Reversed ranges, malformed input, ambiguous delimiters and int32 overflow are errors. Pinned class-test cases `1,2,3`, `+2-+5`, `-1,-2,-3` and `2` are covered, as are signed int32 boundaries.

Search-parameter mergeability compares database basenames across either path separator, database version, taxonomy, exact charge declaration, tolerances including units, and enzyme name/regex/specificity. Fixed and variable modification lists compare as sets. `labeled_MS1` permits different modification sets. **Mass type and missed-cleavage count are ignored by this source operation**, and the native method preserves that behavior. Use full structural equality if those settings must also match. No filesystem paths are opened during these comparisons.

## Kernel attachment and file boundaries

`MSSpectrum` and `BaseFeature` own `peptide_identifications`. Both `FeatureMap` and `ConsensusMap` own `protein_identifications` and `unassigned_peptide_identifications`. Validation traverses these records, and ordinary sorting, filtering, selection and cloning preserve their attachment. Clearing a container with `clear(false)` retains its metadata-level records; `clear(true)` resets them. Clearing a feature map removes its features, including their attached IDs, while retaining map-level runs and unassigned IDs when requested.

No run-reference graph is inferred or repaired: an empty or nonmatching run identifier may be retained. Referenced accessions are ordinary strings. This preserves the ability to assemble partial identification results without pretending to implement `IdentificationData` or `IDMapper`.

DTA, MGF and the current mzML writer reject spectra with attached peptide identifications before emitting any bytes. The supported peak-file representations cannot round-trip those structured records. Existing acquisition and annotation limitations in each format remain documented separately. The optional [idXML adapter](IDXML_SUPPORT.md) reads/writes the supported identification model and explicitly rejects fields that cannot be retained. JSON, mzIdentML/pepXML, sequence database search, probabilistic protein inference, USI/UID construction and reference-graph merging remain unimplemented.

## Verification

`tests/identification.rs` exercises source defaults, stable score directions/ties, marker ordering, modified-sequence identity, accession references, typed metadata, target/decoy categories, annotation/group order, source charge-range cases and mergeability. Coverage is checked against an independent bit-array oracle over small overlapping intervals. Additional tests cover terminal/residue modification placement, duplicate and skipped observations, atomic malformed-evidence errors, nested kernel validation, attachment preservation and writer preflight rejection.

The [identification example](../examples/identify_peptides.rs) connects FASTA input, fixed cysteine modification, digestion, fragment matching, owned evidence and protein coverage. `tests/identification_workflow.rs` verifies the selected peptide, annotations, coordinates, mapped modification and metadata against independent expectations. `tests/identification_pipeline.rs` checks an independent idXML input through score switching, top-hit filtering, hand-computable q-values, protein/evidence/group cleanup, coverage and native idXML round trip. `tests/protein_workflow.rs` additionally verifies FASTA indexing, shared evidence, basic protein inference, separate protein/PSM q-values, complete coverage and file-origin partitioning.

The implementation was independently reviewed against `METADATA/{PeptideEvidence,PeptideHit,PeptideIdentification,ProteinHit,ProteinIdentification}.cpp` and corresponding headers/class tests in the immutable source snapshot. Their hashes are recorded in [source-inventory.json](source-inventory.json). The source C++ library was not built or run; these checks are source-derived expectations and native invariants, not runtime differential certification.

Custom modification identity is chemical value identity, not display text.
`PeptideHit::identity_key()` returns an owned `(AASequence, charge)` pair and
agrees with `same_sequence_and_charge`. Named records, anonymous annotations,
and sequences have deterministic equality-consistent ordering, including
signed-zero mass values. Protein observations retain their original position/
full-ID ordering, with complete chemical values breaking identifier ties.
