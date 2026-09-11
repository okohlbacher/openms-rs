# Pinned ontology resource notices

These five files are unmodified copies from OpenMS4-core revision `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`, under `share/OpenMS/CV`. They are scientific data with their own provenance, separate from the Rust/OpenMS BSD-3-Clause implementation license. Original headers, creator lists, dates, citations and URLs are retained verbatim in the files.

| File | Included resource date/version | Byte handling |
| --- | --- | --- |
| psi-ms.obo | 4.1.155; 2024-06-03 | Original UTF-8 bytes |
| quality.obo | 2008-07-30 | Original UTF-8 bytes |
| unit.obo | 2011-10-12 | Original UTF-8 bytes |
| brenda.obo | 2008-05-01; revision 1.30 | Original legacy bytes; provider explicitly transcodes Windows-1252 |
| goslim_goa.obo | 2004-09-27 | Original ASCII bytes |

The provider does not update terms from the network. Read this table as the included historical resource version, not a claim that these are the latest ontology releases. Projection fixtures identify their derivation and the decoding/prefix corrections separately.

## License and attribution evidence

- **PSI-MS:** the included file explicitly provides the ontology under [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/). Its header lists the HUPO-PSI contributors and links [the ontology](http://purl.obolibrary.org/obo/ms/psi-ms.obo). Preserve those original attributions and this version when using the data.
- **BTO:** the [BRENDA Tissue Ontology project](https://github.com/BRENDA-Enzymes/BTO) states CC BY 4.0, with a [dedicated license text](https://github.com/BRENDA-Enzymes/BTO/blob/master/LICENSE.MD). Its requested citation is Gremse et al., “The BRENDA Tissue Ontology (BTO): the first all-integrating ontology of all organisms for enzyme sources,” Nucleic Acids Research39:D507–D513 (2011). The included 2008 file itself does not have a license statement; this project-policy evidence is recorded separately from its historical date. Do not conflate BTO terms with a license for the entire BRENDA enzyme database.
- **UO:** the [current owner repository](https://github.com/bio-ontology-research-group/unit-ontology) states CC BY 4.0. The [OBO Foundry catalog](https://obofoundry.org/ontology/uo.html) still lists CC BY 3.0. The included 2011 file contains neither statement. Both sources are retained as provenance evidence; this notice does not silently infer the exact historical license version from a current catalog.
- **PATO:** the [OBO catalog](https://obofoundry.org/ontology/pato.html) lists CC BY 3.0. The [owner repository license](https://github.com/pato-ontology/pato/blob/master/LICENSE) is a BSD-3-Clause software-style notice dated 2014. The included 2008 quality ontology file lacks an embedded license. These are distinct current project/data-license records, not proof that the historical file was retroactively relicensed as software.
- **GO:** the [Gene Ontology Consortium citation/license policy](https://geneontology.org/docs/go-citation-policy/) provides GO data/products under CC BY 4.0 and asks that the exact release be identified. This package contains the 2004 GOA-slim file, including its original citation, not current full GO. Retain the Consortium attribution and this historical file date.

Owner/catalog sources were checked 2026-09-11. Historical snapshot licensing details for PATO/UO remain explicit provenance limitations. These notices do not replace the named licenses, assert ownership of third-party definitions, or label all ontology data BSD-3-Clause. Keep original creator/citation headers and applicable license links with redistributed resources.

## Existing source damage

The original brenda.obo contains 335 bytes invalid under UTF-8 and eight NUL bytes in the `BTO:0002243` hypanthium definition. The raw bytes remain unchanged; the provider's explicit Windows-1252 decoding retains the damaged description as opaque data. No inferred scientific text replaces it. XML consumers validate only fields actually rendered. Exact SHA-256 hashes, offsets and source paths are in [the provenance manifest](../../tests/data/controlled_vocabulary_provenance.json).
