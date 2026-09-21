# mzML schema validation

The optional, default-off `mzml-schema` feature implements the class-specific `MzMLFile::isValid` operation with a real XSD engine. It uses exact registry **libxml 0.3.14**, backed by installed libxml2, and the unchanged ordinary/indexed mzML 1.1 schemas from SDK `82ce5b373c97f934ffd9b1ffd80215ca66473d0b`. It does not use a subprocess, an approximate schema checker, or application-provided schema URLs.

## API and source policy

Available through `format::mzml` and `format::mzml_schema`:

| Native operation | Behavior |
| --- | --- |
| `validate_schema(path)` | Default limits; shared gzip/bzip2 magic detection, independent of suffix. |
| `validate_schema_with_options(path, &options)` | Explicit native input/preflight/report limits. |
| `validate_schema_reader(reader, &options)` | Validate a caller-owned `BufRead`; the caller supplies uncompressed XML. |

`SchemaValidationReport` owns the selected `SchemaKind`, structured diagnostic severity, message, optional filename/line/column, and engine domain/code. `is_valid()` requires successful validation with no Warning, Error or Fatal diagnostics. Information (`None` engine severity) does not invalidate. A schema violation returns `Ok(report)` with false validity. Malformed or unsupported XML, I/O, resource limits, setup and recoverable engine failures return the existing typed `Error`. No caller output is partially replaced. Contexts and temporary documents are local to each call.

This preserves the warning policy of `VALIDATORS/XMLValidator.cpp:31–33,107–111`, which marks warning, error and fatal callbacks invalid. Diagnostic wording and codes come from the installed libxml2, not Xerces, and are not promised stable across C-library versions. The safe binding does not expose owned structured document-parser warnings. Schema compiler/validator warnings are drained even on successful calls, but incidental document-parser warnings cannot be reported through this API. The native lexical and namespace checks reject the tested malformed inputs before C; this is not a promise of identical arbitrary Xerces/libxml parser diagnostics.

The source test at `MzMLFile_test.cpp:1181–1201` asserts validity of a C++-written empty document, a C++-written loaded document, and the unchanged indexed fixture. Tests retain the original fixtures and distinguish their literal validity checks from separately generated native-writer documents. No C++ executable was run for this group.

## XML selection and offline behavior

Root selection uses the normalized expanded name `{http://psi.hupo.org/ms/mzml}mzML` or `indexedmzML`. Prefix spelling, BOMs and prolog line count do not change the selected schema, correcting [CPP-052](../OpenMS_CPP_ISSUES.md). The shared bounded XML decoder/scanner performs one lexical pass; a private callback scope stack checks namespace bindings, each QName component, undeclared prefixes and duplicate expanded attribute names. Unprefixed attributes have no default namespace. Empty-element and ordinary closing callbacks restore shadowed bindings.

This check is necessary because actual libxml2 2.9.13 `xmlReadMemory` returns a document after some namespace errors even with recovery disabled. Before the guard, an unused illegal `xmlns:xml` redefinition in the original fixture produced `engine_valid=true` with no schema diagnostics. Forbidden ASCII namespace-URI characters and malformed percent escapes are also rejected on normalized declarations, including unused ones. The before-fix observations are retained under `docs/mzml_schema/`; they are binding/parser observations, not additional OpenMS C++ defect claims. The lexical URI check is not a full URI grammar or resolver. Relative namespace references remain accepted (deprecated, not forbidden); associated parser warnings are unavailable. Unicode namespace names are retained. Namespace comparisons do not decode percent escapes or resolve relative paths. See [W3C namespace constraints](https://www.w3.org/TR/REC-xml-names/#ns-decl) and [URI syntax](https://www.rfc-editor.org/rfc/rfc3986).

Input supports the shared XML 1.0 decoder: UTF-8, BOM/declaration-detected UTF-16LE/BE, and ASCII-only payload with compatible Latin-1 declarations. Non-ASCII Latin-1 is explicitly unsupported. The already-validated declaration's encoding value is replaced with UTF-8 after decoding, preserving standalone, quote style and line layout (shared input line endings are normalized). No binding `encoding=Some` override is used. The raw schemas retain their own original encodings: ordinary schema comments contain Windows-1252 bytes; indexed schema is UTF-8.

Every DOCTYPE is rejected before C, including internal/general/parameter entities. Only the two static byte-pinned schemas are compiled; they contain no external include/import/redefine/DTD/entity dependencies. Parser settings explicitly disable recovery, external DTD loading and huge mode, and enable `NONET`. `no_def_dtd=false` is intentional: the binding's misleading option maps to flag4, which enables DTD loading for XML when true. No XInclude, stylesheet, arbitrary grammar or catalog processing is requested. No process-global resolver, error hook or cleanup routine is installed.

A separate integration-test process records all actual libxml input-loader attempts. File/HTTP/custom document schema hints, stylesheet instructions and XInclude produce zero loads; DTDs fail before C. A direct-binding arbitrary import used only as a positive control reaches the observer and receives empty bytes. The test also runs with a custom catalog environment. The public API never exposes that arbitrary import operation.

## What validity establishes

The engine checks real XSD sequence/content models, required attributes, numeric facets, identity keys and key references. Both schemas reject duplicate spectrum keys and dangling software references in the tests.

It does **not** establish CV semantics, binary decoding, declared-count equality, index offsets, checksum correctness, or all index references. The shipped indexed schema has the verified selector/field defect [CPP-054](../OpenMS_CPP_ISSUES.md). Exactly one original `idRef="index=19"` changed to `idRef="does_not_exist"` still validates with zero diagnostics. Raw schemas remain unchanged. The provenance records original/transformed hashes and this exact recipe; independent genuine key/keyref failures prevent confusing that schema defect with an incomplete engine. Use separate semantic/index/binary operations for their own contracts.

## Resource and failure boundary

Defaults cap encoded input and decoded/canonical UTF-8 at 16 MiB, depth at128, elements at1,000,000, shared native work at50,000,000 and logical native allocation at128 MiB. One meter spans reading/decoding, lexical events, attribute copies, namespace scopes/lookups/comparisons, encoding replacement and report publication. Input lengths are also checked against the binding's i32 range before C.

**These are not hard C-engine memory or time limits.** The DOM, identity tables, schema compilation, internal temporary data and validation runtime are opaque. The safe binding accumulates its own unbounded diagnostic vector and strings before returning. The default10,000 diagnostics/1 MiB message+filename limit is checked afterward, before publishing the report; native allocation accounting cannot retroactively bound that callback vector. No timeout, allocator hook or isolation process is provided. Allocator aborts and C faults are not recoverable `Error`s. Ordinary binding panics are caught; the normal process panic hook may still print. Native runtime shared-library loading failures can occur before Rust code starts.

## Optional dependency and platforms

The exact registry archive is SHA-256 `d3e969617eea4e728856e3629229ed5ee92aa0b0a8cc2b7b42171b8a9ba4916e`, upstream commit `95f02f36449d32d4b578cdf769370fe580098346`. Version0.3.14 uses edition2024 and is actually compiled/tested on Rust1.85; newer0.3.21 explicitly requires1.88. No third-party source is patched. The prior isolated spike verified archive/source bytes and safe parser/schema/error/initialization implementations; its conclusions are supplemented by the namespace regressions here.

Default and no-feature builds do not select libxml, bindgen, clang-sys or vcpkg. `tools/check_schema_feature_graph.py` checks both actual normal/build graphs and builds them with intentionally nonexistent native-library overrides. Cargo lockfile metadata resolution may still include optional packages. The OpenMS-owned build script uses std only and performs no probing; only enabled Windows/MSVC builds emit `bcrypt` and `ws2_32`, matching the upstream fix made after0.3.14.

Enabled builds require:

| Platform | Build requirements and linkage |
| --- | --- |
| Linux | `libxml2-dev`, `libclang-dev`/clang, pkg-config; `libxml2-utils` supports existing independent test tooling. Dynamically linked applications need the corresponding libxml2 runtime. |
| macOS | libxml2 development files, LLVM/libclang and pkgconf. CI uses Homebrew `libxml2 llvm pkgconf` with explicit pkg-config/libclang paths. |
| Windows MSVC | LLVM/libclang and vcpkg `libxml2:x64-windows-static-md`; explicit `VCPKG_ROOT`/`VCPKGRS_TRIPLET`. The `-md` triplet matches Rust's default dynamic CRT while linking the native package statically. The two SDK link libraries above require no dependency patch. MinGW is not the configured target. |

The workflow installs these packages before existing all-feature jobs, adds minimal schema tests, and adds Rust1.85 schema coverage on macOS/MSVC alongside the existing Linux minimum job. The portable graph job requires no native package setup. Local execution is macOS arm64, Rust1.85/1.98, libxml2 headers/runtime2.9.13; Linux/MSVC execution is **prepared CI, not locally observed evidence**. Native dependency versions/security maintenance remain the distributor's responsibility.

## Attribution and evidence

OpenMS code is BSD-3-Clause. The Rust libxml binding is MIT; its unchanged notice is [retained here](mzml_schema/libxml-MIT.txt). System libxml2 is a separate MIT-licensed C library whose installed distribution supplies its notices and transitive native dependencies. The exact PSI/OpenMS schemas retain their original creator comments and bytes; neither schema contains a separate embedded license grant. Preserve existing OpenMS/PSI attribution and the repository resource policy. `.gitattributes` disables text conversion for these raw schema files.

`tests/data/mzml_schema_provenance.json` records source paths/hashes, literal source assertions, original resources and projection recipes. Focused public tests cover originals, schema rules, CPP-052/054, namespace/encoding errors, compressed paths, contexts, limits and offline loading. Private tests cover information-versus-warning policy, declaration preservation and cumulative result accounting. Validation logs retain the initial namespace failures and subsequent fixes; only final successful checks count as completed validation.
