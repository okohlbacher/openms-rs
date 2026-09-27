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

## Next step
1. Read existing Rust XML patterns (src/format/transformation_xml.rs, paramxml.rs) + skill.
2. Port PTMXMLFile; commit.
