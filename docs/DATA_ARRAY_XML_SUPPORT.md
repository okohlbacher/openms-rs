# Data-array description transport and guards

Native `DataArray<T>` retains `MetaInfoDescription` metadata and shared
`DataProcessing` descriptions separately from its name and values.
[mzML header/reference transport](MZML_HEADER_SUPPORT.md) now preserves auxiliary
scalar metadata, represented units and processing references for float, integer
and string arrays. Empty arrays and default processing objects remain meaningful
and round-trip through explicit descriptions. Scientific filtering and sorting
keep descriptions attached to the selected arrays.

Primary arrays have no independent native description owner. String metadata
without units merges into the record's string metadata in source order;
non-string/unit-bearing primary metadata and independent processing histories
remain checked errors. Nonscalar or unrepresentable auxiliary metadata also
fails before output.

Identification and map XML retain separate limits. Both protein-group families
in the shared identification XML preflight reject description metadata or
processing vectors, including empty/zero arrays and default pointed-to records.
This prevents quantity projection from silently dropping descriptions. The
check uses container lengths after shared array-slot precharge; no description
payload is cloned to decide rejection. FeatureXML's earlier structured-group
guard and idXML's quantity-array limit remain in force.

[Transport regressions](../tests/data_array_xml.rs) cover mzML round trips for all
three auxiliary types, and unchanged-output failures in consensusXML, idXML and
FeatureXML. [Header regressions](../tests/mzml_header.rs) additionally cover
reference identity, scalar units, primary metadata order and both mzML writers.
Source `MetaInfoDescription` and XML quantity conventions remain pinned in the
broader Mobilogram/DataArray and XML provenance manifests.
