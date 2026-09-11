# Core SDK completion ledger

Target: `54a232fe2cae9c590d5c997fa49d20e7769860fb`. This ledger covers all **786 registered public headers** and direct includes from **146 TOPP source files**.

This is a work inventory, not a completion percentage. Matching declarations and source references remain unverified until each API and its behavior are reviewed. No TOPP workflow is yet certified as port-ready. Physical unregistered headers and product backends are tracked separately by the SDK source inventory.

| Review state | Headers |
| --- | ---: |
| complete | 10 |
| evidence_requires_review | 163 |
| native_equivalent | 10 |
| partial | 11 |
| unmapped | 592 |

## Highest fan-out open SDK dependencies

These counts show direct consumers; they do not establish full dependency closure.

| Header | Direct TOPP consumers | State |
| --- | ---: | --- |
| `OpenMS/FORMAT/FileHandler.h` | 118 | partial |
| `OpenMS/CONCEPT/LogStream.h` | 76 | partial |
| `OpenMS/KERNEL/MSExperiment.h` | 70 | partial |
| `OpenMS/METADATA/ProteinIdentification.h` | 56 | evidence_requires_review |
| `OpenMS/SYSTEM/File.h` | 53 | partial |
| `OpenMS/KERNEL/ConsensusMap.h` | 40 | evidence_requires_review |
| `OpenMS/KERNEL/FeatureMap.h` | 22 | evidence_requires_review |
| `OpenMS/FORMAT/MzMLFile.h` | 19 | partial |
| `OpenMS/CONCEPT/Constants.h` | 16 | evidence_requires_review |
| `OpenMS/PROCESSING/ID/IDFilter.h` | 14 | evidence_requires_review |
| `OpenMS/CHEMISTRY/ProteaseDB.h` | 13 | evidence_requires_review |
| `OpenMS/METADATA/PeptideIdentification.h` | 12 | evidence_requires_review |
| `OpenMS/CHEMISTRY/ModificationsDB.h` | 11 | evidence_requires_review |
| `OpenMS/FORMAT/FeatureXMLFile.h` | 11 | partial |
| `OpenMS/FORMAT/MzTabFile.h` | 11 | unmapped |
| `OpenMS/MATH/MathFunctions.h` | 11 | evidence_requires_review |
| `OpenMS/CONCEPT/Exception.h` | 10 | unmapped |
| `OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h` | 10 | unmapped |
| `OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h` | 9 | evidence_requires_review |
| `OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h` | 9 | unmapped |
| `OpenMS/FORMAT/ConsensusXMLFile.h` | 9 | partial |
| `OpenMS/FORMAT/IdXMLFile.h` | 9 | evidence_requires_review |
| `OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h` | 8 | unmapped |
| `OpenMS/FORMAT/ExperimentalDesignFile.h` | 8 | unmapped |
| `OpenMS/MATH/StatisticFunctions.h` | 8 | evidence_requires_review |
| `OpenMS/ANALYSIS/ID/PeptideIndexing.h` | 7 | evidence_requires_review |
| `OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h` | 7 | unmapped |
| `OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h` | 7 | unmapped |
| `OpenMS/FORMAT/MzTab.h` | 7 | unmapped |
| `OpenMS/FORMAT/PepXMLFile.h` | 7 | unmapped |
| `OpenMS/FORMAT/SVOutStream.h` | 7 | unmapped |
| `OpenMS/IONMOBILITY/IMTypes.h` | 7 | unmapped |
| `OpenMS/METADATA/ExperimentalDesign.h` | 7 | unmapped |
| `OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h` | 6 | unmapped |
| `OpenMS/CHEMISTRY/ProteaseDigestion.h` | 6 | evidence_requires_review |

## Completion requirements

1. Review every public API against the exact source revision and configuration; record a tested native implementation or an explicit standard-library equivalent.
2. Finish partial containers, formats, numerical routines and domain algorithms; preserve required metadata and error behavior.
3. Check transitive tool dependencies, resources and platform backends, then run ported TOPP workflows against C++ results.
4. Validate supported Rust versions, feature combinations, release builds and performance on realistic inputs.

The complete per-header and per-tool records are in [core-sdk-coverage.json](core-sdk-coverage.json). Explicit reviews are maintained in [core-sdk-reviewed-apis.json](core-sdk-reviewed-apis.json); TOPP source hashes and includes are in [topp-source-inventory.json](topp-source-inventory.json).

Regenerate with `python3 tools/core_sdk_coverage.py --write`. CI checks that the ledger agrees with the source inventory, current Rust declarations and evidence. To replace the local TOPP snapshot, pass `--topp-source /path/to/OpenMS/src/topp --write`.
