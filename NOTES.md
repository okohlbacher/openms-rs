# F7b lane notes (restart-safe)

Branch: worktree-wf_cd066dab-0a7-13, base f0990c0.

## Headers and current coverage status
- FORMAT/GzipIfstream.h, Bzip2Ifstream.h, ZipIfstream.h, ZipArchiveFile.h: unmapped (no Rust).
- FORMAT/HANDLERS/StringManager.h, XMLAttributes.h: unmapped.
- FORMAT/IndentedStream.h: partial (src/cli/usage.rs `Indented`, crate-internal; gap = no public type, OpenMSInfo unported).
- FORMAT/Base64.h: partial (gap = private SIMD encoder/decoder pair; base64-simd measured, not adopted).
- SHARED-02 = ZIP-input gap named in five ledger entries (grep .planning for SHARED-02).

## Next step
1. Read C++ headers/cpp/class tests at core bc9cc12 for the eight headers.
2. grep SHARED-02 in .planning / docs to see what ZIP input must deliver.
3. Check Cargo.toml for flate2 / bzip2 / zip crates already present.
