# Indexed mzML checksum dependency

`sha1 = "=0.10.7"`, `default-features = false`, `features = ["force-soft"]`, optional and activated only by `mzml`, supplies the private incremental mzML checksum. Its API is safe; its dependency internals do not change this crate's `unsafe_code = "forbid"` policy. SHA-1 is used for the format's prescribed checksum, not authentication.

The source package declares `MIT OR Apache-2.0`, with unchanged `LICENSE-MIT` and `LICENSE-APACHE` in its published archive. See the [versioned source package](https://docs.rs/crate/sha1/0.10.7/source/), [MIT notice](https://docs.rs/crate/sha1/0.10.7/source/LICENSE-MIT), and [Apache notice](https://docs.rs/crate/sha1/0.10.7/source/LICENSE-APACHE). Binary redistributors retain the chosen dependency notices alongside the existing dependency notices. No vendored source is modified.

The isolated lockfile adds `sha1 0.10.7`, `digest 0.10.7`, `block-buffer 0.10.4`, `crypto-common 0.1.7`, `generic-array 0.14.7`, `typenum 1.20.1`, `version_check 0.9.5`, and `cpufeatures 0.2.17`; already present `cfg-if`/`libc` versions are reused. `force-soft` selects the portable software compressor. The package's README states MSRV 1.41; the integrated Rust 1.85 build is separately required and recorded. These are dependency references, not a relicensing of dependency code under OpenMS's BSD license.
