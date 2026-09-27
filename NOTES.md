# B8c lane notes (BLOCK-05: TOPPBase remainder)

Merged prerequisite worktree-wf_ee239125-5be-1 (fast-forward to 9fcb962).
No earlier incarnation notes.

## Scope (from brief, ranked)
1. -write_ini numeric formatting: Release writes `4.0`, `1.0e20`, `0.0:`; Rust writes `4`,
   `100000000000000000000`, `0:`. Integer beyond 32 bits: C++ Param truncates, Rust refuses.
   Oracle: ../oracle/json-ini/run1/results/{w_pphr_int_huge,r_pphr_huge}. Files:
   src/format/paramxml.rs, src/param/** (integer width ONLY; B6 owns param otherwise).
2. JSON INI syntax-error text: nlohmann wording vs serde_json -> src/format/param_json.rs
3. parse error prints `parse error on line N: <msg>`; source prints `<msg> in: <file>`
4. Error::MissingInformation exits 7; source 8
5. FORMAT/ParamCWLFile.h
Windows console width: recorded, not chased.

## Next step
1. Read oracle cases w_pphr_int_huge / r_pphr_huge and the C++ writer (ParamXMLFile.cpp, Param value toString).
