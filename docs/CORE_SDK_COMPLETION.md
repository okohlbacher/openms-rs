# Core SDK completion ledger

Target: `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. This ledger covers all **786 registered public headers** and direct includes from **124 TOPP source files**.

This is a work inventory, not a completion percentage. Matching declarations and source references remain unverified until each API and its behavior are reviewed. A TOPP workflow counts as validated only when an executed differential comparison against retained C++ output is recorded in its provenance manifest. Physical unregistered headers and product backends are tracked separately by the SDK source inventory.

Validated TOPP workflows: **8** of 124, each reproducing its upstream test against retained C++ output.

| Review state | Headers |
| --- | ---: |
| complete | 63 |
| evidence_requires_review | 165 |
| native_equivalent | 90 |
| partial | 59 |
| unmapped | 409 |

## Highest fan-out open SDK dependencies

These counts show direct consumers; they do not establish full dependency closure.

| Header | Direct TOPP consumers | State |
| --- | ---: | --- |
| `OpenMS/FORMAT/FileHandler.h` | 100 | partial |
| `OpenMS/CONCEPT/LogStream.h` | 69 | partial |
| `OpenMS/KERNEL/MSExperiment.h` | 55 | partial |
| `OpenMS/METADATA/ProteinIdentification.h` | 47 | partial |
| `OpenMS/DATASTRUCTURES/StringUtils.h` | 12 | partial |
| `OpenMS/FORMAT/MzMLFile.h` | 12 | partial |
| `OpenMS/CHEMISTRY/ProteaseDB.h` | 11 | evidence_requires_review |
| `OpenMS/METADATA/PeptideIdentification.h` | 11 | evidence_requires_review |
| `OpenMS/PROCESSING/ID/IDFilter.h` | 9 | evidence_requires_review |
| `OpenMS/CHEMISTRY/ModificationsDB.h` | 8 | evidence_requires_review |
| `OpenMS/FORMAT/FeatureXMLFile.h` | 8 | partial |
| `OpenMS/MATH/MathFunctions.h` | 8 | partial |
| `OpenMS/FORMAT/ConsensusXMLFile.h` | 7 | partial |
| `OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h` | 7 | partial |
| `OpenMS/FORMAT/MzTabFile.h` | 7 | partial |
| `OpenMS/FORMAT/IdXMLFile.h` | 6 | partial |
| `OpenMS/FORMAT/MzTab.h` | 6 | partial |
| `OpenMS/FORMAT/QcMLFile.h` | 6 | partial |
| `OpenMS/CHEMISTRY/ProteaseDigestion.h` | 5 | evidence_requires_review |
| `OpenMS/CONCEPT/VersionInfo.h` | 5 | evidence_requires_review |
| `OpenMS/FORMAT/OMSFile.h` | 5 | unmapped |
| `OpenMS/FORMAT/PepXMLFile.h` | 5 | partial |
| `OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h` | 5 | partial |
| `OpenMS/SYSTEM/StopWatch.h` | 5 | partial |
| `OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h` | 4 | evidence_requires_review |
| `OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h` | 4 | unmapped |
| `OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h` | 4 | unmapped |
| `OpenMS/ANALYSIS/ID/PeptideIndexing.h` | 4 | evidence_requires_review |
| `OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h` | 4 | unmapped |
| `OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h` | 4 | evidence_requires_review |
| `OpenMS/CONCEPT/Exception.h` | 4 | evidence_requires_review |
| `OpenMS/IONMOBILITY/IMDataConverter.h` | 4 | partial |
| `OpenMS/SYSTEM/JavaInfo.h` | 4 | partial |
| `OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h` | 3 | evidence_requires_review |
| `OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h` | 3 | unmapped |

## Completion requirements

1. Review every public API against the exact source revision and configuration; record a tested native implementation or an explicit standard-library equivalent.
2. Finish partial containers, formats, numerical routines and domain algorithms; preserve required metadata and error behavior.
3. Check transitive tool dependencies, resources and platform backends, then run ported TOPP workflows against C++ results.
4. Validate supported Rust versions, feature combinations, release builds and performance on realistic inputs.

The complete per-header and per-tool records are in [core-sdk-coverage.json](core-sdk-coverage.json). Explicit reviews are maintained in [core-sdk-reviewed-apis.json](core-sdk-reviewed-apis.json); TOPP source hashes and includes are in [topp-source-inventory.json](topp-source-inventory.json).

Regenerate with `python3 tools/core_sdk_coverage.py --write`. CI checks that the ledger agrees with the source inventory, current Rust declarations and evidence. To replace the local TOPP snapshot, pass `--topp-source /path/to/OpenMS/src/topp --write`.
