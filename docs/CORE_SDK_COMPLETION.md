# SDK completion ledger: core and cli packages

Targets: core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`, cli `c19e49414bcd9ebdea42f89b3f74d2823205892c`. This ledger covers all **794 registered public headers** of those packages (786 core and 8 cli) and direct includes from **124 TOPP source files**.

This is a work inventory, not a completion percentage. Matching declarations and source references remain unverified until each API and its behavior are reviewed. A TOPP workflow counts as validated only when an executed differential comparison against retained C++ output is recorded in its provenance manifest. Physical unregistered headers and product backends are tracked separately by the SDK source inventory.

Validated TOPP workflows: **9** of 124, each reproducing its upstream test against retained C++ output.

| Review state | core | cli | Headers |
| --- | ---: | ---: | ---: |
| complete | 125 | 0 | 125 |
| evidence_requires_review | 23 | 0 | 23 |
| native_equivalent | 91 | 0 | 91 |
| partial | 138 | 2 | 140 |
| unmapped | 409 | 6 | 415 |

## Highest fan-out open SDK dependencies

These counts show direct consumers; they do not establish full dependency closure.

| Header | Package | Direct TOPP consumers | State |
| --- | --- | ---: | --- |
| `OpenMS/APPLICATIONS/TOPPBase.h` | cli | 108 | partial |
| `OpenMS/FORMAT/FileHandler.h` | core | 100 | partial |
| `OpenMS/CONCEPT/LogStream.h` | core | 69 | partial |
| `OpenMS/KERNEL/MSExperiment.h` | core | 55 | partial |
| `OpenMS/METADATA/ProteinIdentification.h` | core | 47 | partial |
| `OpenMS/DATASTRUCTURES/StringUtils.h` | core | 12 | partial |
| `OpenMS/FORMAT/MzMLFile.h` | core | 12 | partial |
| `OpenMS/CHEMISTRY/ProteaseDB.h` | core | 11 | partial |
| `OpenMS/METADATA/PeptideIdentification.h` | core | 11 | partial |
| `OpenMS/PROCESSING/ID/IDFilter.h` | core | 9 | partial |
| `OpenMS/CHEMISTRY/ModificationsDB.h` | core | 8 | partial |
| `OpenMS/FORMAT/FeatureXMLFile.h` | core | 8 | partial |
| `OpenMS/MATH/MathFunctions.h` | core | 8 | partial |
| `OpenMS/FORMAT/ConsensusXMLFile.h` | core | 7 | partial |
| `OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h` | core | 7 | partial |
| `OpenMS/FORMAT/MzTabFile.h` | core | 7 | partial |
| `OpenMS/APPLICATIONS/TOPPExternalToolBase.h` | cli | 6 | unmapped |
| `OpenMS/FORMAT/IdXMLFile.h` | core | 6 | partial |
| `OpenMS/FORMAT/MzTab.h` | core | 6 | partial |
| `OpenMS/FORMAT/QcMLFile.h` | core | 6 | partial |
| `OpenMS/CONCEPT/VersionInfo.h` | core | 5 | partial |
| `OpenMS/FORMAT/OMSFile.h` | core | 5 | unmapped |
| `OpenMS/FORMAT/PepXMLFile.h` | core | 5 | partial |
| `OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h` | core | 5 | partial |
| `OpenMS/SYSTEM/StopWatch.h` | core | 5 | partial |
| `OpenMS/APPLICATIONS/MapAlignerBase.h` | cli | 4 | unmapped |
| `OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h` | core | 4 | partial |
| `OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h` | core | 4 | unmapped |
| `OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h` | core | 4 | unmapped |
| `OpenMS/ANALYSIS/ID/PeptideIndexing.h` | core | 4 | partial |
| `OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h` | core | 4 | unmapped |
| `OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h` | core | 4 | partial |
| `OpenMS/CONCEPT/Exception.h` | core | 4 | partial |
| `OpenMS/IONMOBILITY/IMDataConverter.h` | core | 4 | partial |
| `OpenMS/SYSTEM/JavaInfo.h` | core | 4 | partial |

## Completion requirements

1. Review every public API against the exact source revision and configuration; record a tested native implementation or an explicit standard-library equivalent.
2. Finish partial containers, formats, numerical routines and domain algorithms; preserve required metadata and error behavior.
3. Check transitive tool dependencies, resources and platform backends, then run ported TOPP workflows against C++ results.
4. Validate supported Rust versions, feature combinations, release builds and performance on realistic inputs.

The complete per-header and per-tool records are in [core-sdk-coverage.json](core-sdk-coverage.json). Explicit reviews are maintained in [core-sdk-reviewed-apis.json](core-sdk-reviewed-apis.json); the package inventories are [core-sdk-update.json](core-sdk-update.json) and [cli-sdk-update.json](cli-sdk-update.json); TOPP source hashes and includes are in [topp-source-inventory.json](topp-source-inventory.json).

Regenerate with `python3 tools/core_sdk_coverage.py --write`. CI checks that the ledger agrees with the package inventories, current Rust declarations and evidence. To replace the local TOPP snapshot, pass `--topp-source /path/to/OpenMS/src/topp --write`. To check the cli inventory against its pin, pass `--cli-source /path/to/OpenMS4-cli`, and add `--write` to regenerate it.
