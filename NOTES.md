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

## Done so far (uncompiled until the dax gate reports)
- Cargo feature `kb-xml` = quick-xml + file-compression (Cargo.toml, report it).
- src/format/sax_xml.rs: shared bounded SAX driver (qnames, attr normalisation, coalesced text).
- src/format/ptm_xml.rs + tests/ptm_xml.rs (fixture tests/data/PTMXMLFile_PTMs.xml).
- src/format/unimod_xml.rs (+ UnimodXMLDataProvider) + tests/unimod_xml.rs
  (Retain default = matches ModificationsDB::global(); ReadOptions::source() = C++ classification rule).
- src/format/uniprot_xml.rs + tests/uniprot_xml.rs (fixture tests/data/UniProtXMLFile_test.xml).
- PEFF fixtures copied to tests/data/PEFFFile_*.peff (11 files, unchanged).
- Gate: logs in $SCRATCH/F6/logs/*.log; slot F6 on dax.

## Next step
1. Get kb-xml slice green on dax (ptm_xml, unimod_xml, uniprot_xml, lib sax_xml).
2. Docs: docs/{PTM_XML,UNIMOD_XML,UNIPROT_XML}_SUPPORT.md + provenance JSONs; commit.
3. PEFF: read PEFFFile.cpp 1302-3464 (I/O) then 504-1300 (sequence methods); port to src/format/peff.rs.
4. C++ oracle on ibminode06 (../oracle/F6-kbxml/): dump unimod/uniprot/ptm/peff results.
