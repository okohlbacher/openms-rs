# Pinned enzyme source data

`Enzymes.xml` is an unchanged copy of OpenMS4-core
`share/OpenMS/CHEMISTRY/Enzymes.xml` at revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`.

Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin.
BSD-3-Clause; see the repository [LICENSE](../../LICENSE).

SHA-256: `2f160f3ec32db6257cb4eee43fd594b48cb36398a17b7af21bbba76bc8b16995`.

`python3 tools/generate_enzymes.py --check` verifies the generated native table and
independent regex test hashes from this source. Normal library use needs no XML
parser or runtime resource files. See [digestion support](../../docs/DIGESTION_SUPPORT.md).

## RNA enzymes

`Enzymes_RNA.xml` is the unchanged pinned `share/OpenMS/CHEMISTRY/Enzymes_RNA.xml`,
SHA-256 `4c90a3630a1553ef7c75bae00f2731e2f5bb690f3ee9434d54cc3315a95fe899`.
`python3 tools/generate_rnases.py --check` verifies all fourteen generated records
in `rna_enzymes.rs`. This is OpenMS enzyme data under the same software notice
above. [RNA digestion support](../../docs/RNASE_SUPPORT.md) describes native
per-code matching, custom records and the remaining XML-input/regex scope.
