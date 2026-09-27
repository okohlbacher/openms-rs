# F7a lane notes (restart-safety scratch; delete in final commit)

Lane: XML handlers + MS data consumers. Worktree branch: worktree-wf_ee239125-5be-17, base f0990c0.

## Coverage status at base (docs/core-sdk-coverage.json)
- DATAACCESS: MSDataAggregatingConsumer, MSDataCachedConsumer (TOPP FileConverter), MSDataChainingConsumer,
  MSDataSqlConsumer, MSDataStoringConsumer, NoopMSDataConsumer = unmapped;
  MSDataTransformingConsumer = evidence_requires_review (NovorAdapter); MSDataWritingConsumer = partial
  (src/format/ms_data_writing_consumer.rs; gaps: per-record semantic validation + progress).
- HANDLERS: ConsensusXMLHandler, FeatureXMLHandler, XMLHandler, MzML*, MzData, MzXML = partial (see review scopes).
- FORMAT: XMLFile partial, IndexedMzMLFileLoader partial (missing type + options_ member + OnDiscPeakMap typedef).

## Plan (dependency order)
1. Consumer trait family: find existing IMSDataConsumer equivalent (grep ms_data_writing_consumer.rs).
2. Port Noop, Storing, Chaining, Transforming, Aggregating, Cached, Sql consumers (small headers) -> src/format/dataaccess/ or alongside.
3. IndexedMzMLFileLoader type (join load_experiment + writing consumer).
4. Re-read handler partial scopes; close small gaps.

## Next step
1. Read src/format/ms_data_writing_consumer.rs for the consumer trait; read C++ headers of the 7 consumers at bc9cc12.
