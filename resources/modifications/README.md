# OpenMS Rust Modification Table

This directory contains an explicitly transformed dataset named **OpenMS Rust
Modification Table**, generated on **2026-09-10**. It is derived from the UniMod
and custom modification XML distributed with OpenMS4-core revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`.

- `unimod.xml` is the complete unchanged source data from
  `share/OpenMS/CHEMISTRY/unimod.xml`: Copyright (C) 2002-2006 Unimod. This
  information may be copied, distributed and/or modified under the accompanying
  **Design Science License**, and comes without any warranty.
- `custom_mods.xml` is the complete unchanged OpenMS custom table from the same
  source directory, covered by the OpenMS BSD-3-Clause notice in the root LICENSE.
- `openms-rust-modifications.tsv` contains 3,035 specificity records flattened
  from 1,543 UniMod entries and 14 custom entries. The derivative UniMod data
  retains the Design Science License. The implementation code is BSD-3-Clause.
- `DESIGN_SCIENCE_LICENSE.txt` is the complete license text, copied from
  [GNU's license archive](https://www.gnu.org/licenses/dsl.html).
- `provenance.json` records the source and generated-data hashes.

The transformation flattens each modification into one row per specificity,
canonicalizes declared delta and neutral-loss elemental compositions, retains
declared mono/average mass differences, and drops XML-only descriptions, cross
references, ignored elements and brick definitions. The original XML is included
so the preferred source data accompanies the derivative. No original source
data is modified. Regenerate deterministically with:

```sh
python3 tools/generate_modifications.py
```

TSV fields are: numeric record ID, short name, full name, source site, normalized
terminal specificity, monoisotopic delta, average delta, empirical formula,
hidden flag, classification, and semicolon-separated neutral losses. Each loss
is `formula@mono_mass@average_mass`. Leading `#` lines are comments.

The runtime maps a terminal site with `anywhere` position to the wildcard `X`,
matching the pinned OpenMS XML reader; otherwise a terminal site accepts any
residue at that specified terminus. Names/accessions can therefore be ambiguous
without a residue and terminal specificity. These tables are pinned historical
data, not a claim to contain the latest UniMod release.

The registry also loads the unchanged `XLMOD.obo` ontology dated 2016-07-13,
with 92 mono-link specificity records appended after the 3,035 table records.
CrossLinksDB separately loads its 56 cross-link specificity records. XLMOD is
licensed **CC-BY-3.0**; see `XLMOD_LICENSE.md` and `xlmod_provenance.json` for
attribution, source, and exact hashes. This additional provider does not change
the 11-field TSV schema or its deterministic generator. Historical PSI-MOD is
not bundled; callers can supply an OBO stream through the bounded native reader.
