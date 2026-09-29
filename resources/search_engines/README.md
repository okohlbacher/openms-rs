# Search-engine resources

`OMSSA_modification_mapping` is copied unmodified from
`share/OpenMS/CHEMISTRY/OMSSA_modification_mapping` of OpenMS4-core at
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` (BSD-3-Clause), sha256
`e7b8d09146a656f944a8f119ab147c5b0a2bef8968d8ec00e7f261f0d51cc5fb`, 4495 bytes.

It maps OMSSA's internal modification numbers to UniMod names, one
`number,description,name[,name...]` line each; `#` starts a comment. The source
`OMSSAXMLFile` reads it from the share directory on construction
(`File::find("CHEMISTRY/OMSSA_modification_mapping")`); the Rust
`format::omssa_xml_file` embeds it with `include_str!` and also accepts a
caller's own mapping text.
