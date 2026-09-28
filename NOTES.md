# F6 lane notes (Knowledge-base XML) — restart state

Scratch (lane-private): $SCRATCH/F6/cpp  (git archive of core@bc9cc12 lane files; the
scratchpad root is shared with other lanes, so only use the F6/ subdir).

## What is missing (all unmapped at f0990c0 except where noted)
- FORMAT/PEFFFile.h (637 h / 3464 cpp) — test PEFFFile_test.cpp (2475 lines), 11 .peff fixtures
  in src/tests/class_tests/openms/data/PEFFFile_*.peff. TOPP consumer UniPEFF.
- FORMAT/UniProtXMLFile.h (117/29) + HANDLERS/UniProtXMLHandler.h (110/396) — test
  UniProtXMLFile_test.cpp, fixture UniProtXMLFile_test.xml. TOPP consumer UniPEFF.
- FORMAT/UnimodXMLFile.h (57/34) + HANDLERS/UnimodXMLHandler.h (89/234) — test
  UnimodXMLFile_test.cpp (uses share/OpenMS/CHEMISTRY/unimod.xml).
  Status evidence_requires_review; tests/data/modification_*_provenance.json reference it.
- CHEMISTRY/UnimodXMLDataProvider.h (36/42) — related gap; unimod.xml flattened offline today
  (tools/generate_modifications.py).
- FORMAT/PTMXMLFile.h (55/33) + HANDLERS/PTMXMLHandler.h (55/70) — test PTMXMLFile_test.cpp,
  fixture PTMs.xml.

## Order
1. PTMXMLFile + handler (small, self-contained) 2. UnimodXMLFile + handler + DataProvider
3. UniProtXMLFile + handler 4. PEFFFile (largest; UniPEFF consumer).

## Done so far (kb-xml slice green on dax: ptm_xml 6, unimod_xml 7, uniprot_xml 4, lib 12)
- Cargo feature `kb-xml` = quick-xml + file-compression (Cargo.toml, report it).
- src/format/sax_xml.rs: shared bounded SAX driver (qnames, attr normalisation, coalesced text).
- src/format/ptm_xml.rs + tests/ptm_xml.rs (fixture tests/data/PTMXMLFile_PTMs.xml).
- src/format/unimod_xml.rs (+ UnimodXMLDataProvider) + tests/unimod_xml.rs
  (Retain default = matches ModificationsDB::global(); ReadOptions::source() = C++ classification rule).
- src/format/uniprot_xml.rs + tests/uniprot_xml.rs (fixture tests/data/UniProtXMLFile_test.xml).
- PEFF fixtures copied to tests/data/PEFFFile_*.peff (11 files, unchanged).
- Gate: logs in $SCRATCH/F6/logs/*.log; slot F6 on dax.

## Incarnation 3 (branch worktree-wf_cd066dab-0a7-11)
Merged worktree-wf_ee239125-5be-16 (incarnation 2). It holds an UNGATED PEFF port saved by the lead:
src/format/peff.rs + peff/{read,text,write}.rs + tests/peff.rs (~4.1k lines). First: gate it on dax.

## PEFF gap found (incarnation 3)
PEFF I/O is ported; PEFFEntry sequence members are NOT: getSequence, getModifiedSequence,
getVariantSequences, getProcessedSequence, digestWithVariants, generatePeptides,
enumeratePEFFModifications_ (C++ PEFFFile.cpp 528-1300). Rust deps exist:
chemistry::sequence::AASequence (parse, set_modification), ModificationsDB::global().get_modification,
chemistry::digestion::ProteaseDigestion (digest_ranges; min/max in struct; C++ passes them),
chemistry::modified_peptides::ModifiedPeptideGenerator. C++ tests: PEFFFile_test.cpp 365-431,
882-1335, 1612-2350. Plan: src/format/peff/sequences.rs + tests/peff_sequences.rs.
C++ defect candidate: getVariantSequences / digestWithVariants push the description before
AASequence::fromString, so a throw desynchronises descriptions and sequences.

## Done (incarnation 3)
- PEFF compile fixes; PEFFEntry sequence members in src/format/peff/sequences.rs.
- tests/peff_sequences.rs: 18 class-test ports + 8 oracle tests (26 green on dax).
- Oracle ../oracle/F6-peff-sequences (ibminode06 Release, 2 identical runs) -> tests/data/peff_sequences_oracle.tsv.
- docs/{PEFF,PTM_XML,UNIMOD_XML,UNIPROT_XML}_SUPPORT.md, tests/data/*_provenance.json.
- Ledger fragment written (.planning/ledger-fragments/wave4/F6.json).
- MSRV 1.85 kb-xml slice: peff 32, peff_sequences 18 (pre-oracle), ptm 6, unimod 7, uniprot 4.

## Next step
1. Clippy all-features (running), CI sweep F6-sweep (running) -> fix anything.
2. Measure: no-default-features peff+peff_sequences counts; kb-xml lib unit tests count.
3. Write report .planning/lane-reports/F6.md; delete NOTES.md in final commit.
