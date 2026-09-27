# F3a lane notes (Identification XML) — restart-safe scratch; delete in final commit

Branch: worktree-wf_bcdd783d-f3f-12 (off main f0990c0). Never push/rebase.

## Current coverage status (docs/core-sdk-coverage.json)
- PepXMLFile partial (src/format/pepxml.rs), PepXMLFileMascot UNMAPPED, ProtXMLFile UNMAPPED,
  IdXMLFile partial (identification_xml.rs, idxml.rs), MzIdentMLFile partial (mzidentml.rs),
  MascotXMLFile partial (mascot_xml.rs), XQuestResultXMLFile UNMAPPED, PercolatorOutfile
  evidence_requires_review (no rust), MzIdentMLHandler partial, MascotXMLHandler partial,
  XQuestResultXMLHandler UNMAPPED.

## Known bug to fix first (priority 1)
- AASequence::parse_with_budget (src/chemistry/sequence.rs:1147) charges per modified
  sequence a preflight of annotations * db.len()*(max_id+16)*16 + ... against ONE
  document-wide budget (identification_xml.rs:1276, 1599) -> >14 modified hits fail.
  Fix: amortise the registry-dependent cost (db scan once per document / cache
  max_id,max_atoms,keys) instead of per sequence. Two A7 oracle cases lack differential.

## Plan / order
1. Budget fix (unblocks idXML differentials, MascotXML fixture extractor workaround).
2. PercolatorOutfile (small, IDFileConverter consumer).
3. ProtXMLFile (IDFileConverter consumer).
4. XQuestResultXMLFile + handler (IDFileConverter, OpenPepXL, XFDR).
5. PepXMLFileMascot.
6. Gaps in partials if time.

## Next step
1. Read identification_xml.rs budget use sites and the A7 oracle note; fix budget.
