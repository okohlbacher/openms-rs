# Original schema resources

Every file here is copied byte for byte from the OpenMS core SDK's
`share/OpenMS/SCHEMAS`. The two mzML schemas come from OpenMS4-core revision
`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`; the rest from the current SDK pin,
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`, where the two mzML files are
unchanged.

| Local file | Original source | SHA-256 | Registered by |
| --- | --- | --- | --- |
| mzML_1_10.xsd | share/OpenMS/SCHEMAS/mzML_1_10.xsd | 8be24aab80f1f43a610d84745534e25abcd8734738c2699b0f67b4c242c154db | `MzMLFile`, `ImzMLFile` |
| mzML_idx_1_10.xsd | share/OpenMS/SCHEMAS/mzML_idx_1_10.xsd | b99b668f2d9faabbf62ab8cde1f778d126bf0213d083472590670ee64318f5b0 | `MzMLFile` (indexed) |
| FeatureXML_1_9.xsd | share/OpenMS/SCHEMAS/FeatureXML_1_9.xsd | 9c66e9abeead72ad60296ff442a8029fd04ec76c3ab80c2c0ff46063743b99a3 | `FeatureXMLFile` |
| ConsensusXML_1_7.xsd | share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd | 37e0575633930869c39aa67aad2d8bef28d495df105c4fb474717e2597b370fb | `ConsensusXMLFile` |
| IdXML_1_5.xsd | share/OpenMS/SCHEMAS/IdXML_1_5.xsd | 105c9f0cbc4ec275df8c1689fe31e409557a6535824e13d96ecef138a6986379 | `IdXMLFile` |
| Param_1_8_0.xsd | share/OpenMS/SCHEMAS/Param_1_8_0.xsd | d672b6b90b454ce9fad11524831e9d26ab52daba3f5b4d98beaed441aec2941a | `ParamXMLFile` |
| TrafoXML_1_1.xsd | share/OpenMS/SCHEMAS/TrafoXML_1_1.xsd | 07f16a75b039ea91a0c839c5ed15a3f09e4935ba4598b81b9651e5ff76138e6d | `TransformationXMLFile` |
| mzData_1_05.xsd | share/OpenMS/SCHEMAS/mzData_1_05.xsd | c5c6ad6343d8f87cfd6a81017b30ca30defbed3f6b9a2e14eadc3ab15c084b7f | `MzDataFile` |
| mzXML_idx_3.1.xsd | share/OpenMS/SCHEMAS/mzXML_idx_3.1.xsd | 03ead7c7798cb9f54ca27b8c056ee328093ab5ef94f0a03685cdf559820aa16b | `MzXMLFile` |
| mzXML_3.1_mod.xsd | share/OpenMS/SCHEMAS/mzXML_3.1_mod.xsd | f993a2654a0820a05b99dce5f0f08ca52d11740c89a10ab1a7efb271b69b9376 | included by `mzXML_idx_3.1.xsd` |
| separation_technique_1.0.xsd | share/OpenMS/SCHEMAS/separation_technique_1.0.xsd | d061c1dc60894c268ea9634ddcb671612a93c5f33b18ed036754cb928ec74132 | included by `mzXML_3.1_mod.xsd` |
| general_types_1.0.xsd | share/OpenMS/SCHEMAS/general_types_1.0.xsd | c937229c77b517d5334b30cd4007a74600903e977238fec8bd480b7c2fabd8de | included by `mzXML_3.1_mod.xsd` |
| mzIdentML1.0.0.xsd | share/OpenMS/SCHEMAS/mzIdentML1.0.0.xsd | eaf3cfd583834bb6940725bfeb63ffc8be779bb790c5f5f5bfdc68994494149e | `MzIdentMLFile` (detected 1.0.0) |
| FuGElightv1.0.0.xsd | share/OpenMS/SCHEMAS/FuGElightv1.0.0.xsd | cba9529f8b9c76446114e28f30a04bef5f92e49c2b8fd1ea2a71a8da2f4f9d35 | included by `mzIdentML1.0.0.xsd` |
| mzIdentML1.1.0.xsd | share/OpenMS/SCHEMAS/mzIdentML1.1.0.xsd | 8d12337d8d5abd50a30a68ee540a8180b62a177ee0d8d99d2d0ba30932f6f513 | `MzIdentMLFile` (detected 1.1.0) |
| mzIdentML1.2.0.xsd | share/OpenMS/SCHEMAS/mzIdentML1.2.0.xsd | beca6afe670394bc5edb810632712ea4fd434d08f62491dacc7da83c32ff7f7f | `MzIdentMLFile` (detected 1.2.0) |
| mzIdentML1.3.0.xsd | share/OpenMS/SCHEMAS/mzIdentML1.3.0.xsd | abac61f57e5dcd2eed76a547a3bbd71a47e241d378ac6963fee1f0c0390a6204 | `MzIdentMLFile` (default) |
| pepXML_v114.xsd | share/OpenMS/SCHEMAS/pepXML_v114.xsd | 64a81531831ae268c1fc421fd788d8de2c0ac6cf300362723b01239e3579e2fb | `PepXMLFile` |

Original PSI/OpenMS creator comments, declarations, byte-order marks and
encoding bytes are retained. `mzXML_3.1_mod.xsd` is the source's own locally
modified mzXML 3.1 schema, as its opening comment explains; it is carried as the
source ships it, not as the upstream sashimi file.

The OpenMS-authored schemas contain no separate embedded license grant; the
source repository's BSD-3-Clause LICENSE.md and its hash are retained as
provenance, without asserting a new license for third-party standards content.
The four mzIdentML schemas and `FuGElightv1.0.0.xsd` state that they are
distributed under the Creative Commons Attribution 2.0 license
(<http://creativecommons.org/licenses/by/2.0/>); that statement and their PSI
attribution are retained unchanged in each file. See the
[mzML schema evidence](../../tests/data/mzml_schema_provenance.json), the
[other schemas' evidence](../../tests/data/xml_schema_provenance.json) and the
[project notices](../../LICENSES.md).

The indexed mzML schema's original ID-reference limitation is retained and
documented as [CPP-054](../../OpenMS_CPP_ISSUES.md). No schema has been repaired
or reformatted. The two that `xs:include` others are composed in memory at
validation time (`src/format/xml_schema.rs`, `compose`); the files here stay
unchanged. XSD validity remains separate from index offsets and checksums.
