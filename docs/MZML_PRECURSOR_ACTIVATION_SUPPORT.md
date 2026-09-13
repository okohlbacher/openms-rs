# mzML precursor activation and intensity units

This closes the narrow activation metadata and selected-ion intensity-unit
transport routes in the native mzML reader/writer. It does **not** close the
whole `MzMLFile` or `MzMLHandler` API. Source revision:
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

Implementation: `src/format/mzml_precursor.rs`, the parameter dispatch in
`src/format/mzml_record.rs`, and the precursor reader/writer in
`src/format/mzml.rs`. Public entry points remain `mzml::read`,
`read_with_options`, `write`, and their existing load/streaming/encoding variants.
There is no new public API or dependency.

## Source routes and native representation

| Source route | Native representation |
|---|---|
| selectedIon `MS:1000042`, peak intensity | `Precursor.intensity` (`f32`) |
| selectedIon intensity `unitAccession` | String metadata `peak intensity unit accession`; missing/default `MS:1000132` is implicit |
| activation `MS:1000245`, charge stripping | String metadata `charge stripping = "true"` |
| activation `MS:1000045`, collision energy | Floating metadata `collision energy`, with its unit |
| activation `MS:1000412`, buffer gas | String metadata `buffer gas`, with any explicit unit |
| activation `MS:1000419`, collision gas | String metadata `collision gas`, with any explicit unit |
| activation `MS:1000138`, normalized collision energy | Floating metadata under the source's historical key `percent collision energy` |
| activation `MS:1000869`, collision gas pressure | Floating metadata `collision gas pressure`, with its unit |
| activation `MS:1002679`, supplemental collision-induced dissociation | String metadata of that name, also derives `ActivationMethod::Etcid` |
| activation `MS:1002678`, supplemental beam-type collision-induced dissociation | String metadata of that name, also derives `ActivationMethod::Ethcd` |
| activation `MS:1002680`, supplemental collision energy | Floating metadata `supplemental collision energy`, with its unit |
| activation `MS:1000509`, activation energy | Existing separate `Precursor.activation_energy`, in electronvolts |
| Other existing activation-method CVs | Existing `Precursor.activation_methods` set; unchanged |
| Scalar precursor userParam | Existing typed `cv_terms.metadata`; unchanged except conflicting intensity-unit sources now reject |

Collision energy metadata is distinct from the typed activation-energy field.
A unit retains its accession, CV reference and supplied name. Intensity-unit
metadata retains the accession only, matching the source; output resolves its
canonical name from the pinned controlled vocabulary. Both MS and UO accessions
are accepted if present in that vocabulary. This follows the source lookup;
it does not validate that a known term is physically appropriate for intensity.

Known activation metadata is emitted as a CV parameter when its native type
matches the source reader: floats for numeric terms and strings for text terms.
`charge stripping` is promoted only for the exact unitless string `"true"`.
Other scalar types/values remain userParams, preserving their identity and type.
Both spectrum and chromatogram precursors use the same native transport.

## Preserved conventions and deliberate differences

The pinned `MzMLFile_test.cpp:1420–1434,1484–1485` uses intensity 30 in
`MS:1000131` units, collision and supplemental collision energies 25 in
`UO:0000266`, and an empty supplemental beam-type dissociation term. The native
reader derives EThcD alongside ETD as the source spectrum reader does. An
explicit supplemental term plus the corresponding combined method emits only
the supplemental CV, avoiding a redundant combined-method CV.

A caller can construct supplemental metadata without a combined method. The
native writer keeps that value as a userParam to avoid adding an activation
method when the result is read back. The source writer promotes the metadata,
and its spectrum reader consequently adds a method. To request the CV form in
Rust, include the corresponding combined method in `activation_methods`.

The pinned source chromatogram activation branch omits the three supplemental
routes, although its shared precursor writer can emit them. The native reader
preserves them for chromatograms too. This is based on source inspection; it is
not an executed C++ finding; tracked as CPP-219 in `OpenMS_CPP_ISSUES.md`.

Missing and explicit default intensity units normalize to the same native
representation, matching the source reader. Explicit native metadata naming
`MS:1000132` is rejected on write because that entry would disappear on reload;
remove the entry to use the default. Zero intensity still writes a selected-ion
intensity CV, retaining a nondefault unit even at zero. The C++ writer normally
writes intensity only if positive or its legacy metadata key is present.

Unknown activation CV parameters retain the surrounding reader's established
policy; this change does not introduce generic CV preservation. Invalid or
nonfinite recognized numeric values return an error rather than the source's
warning-and-skip conversion behavior. Charge stripping follows the source's
boolean-presence convention, but rejects units that would otherwise be dropped.

## Checked boundaries and evidence

The existing total XML byte, parameter count and aggregate parameter-byte
limits apply before retention. Each newly retained metadata entry additionally
charges the shared conservative allocation budget. Referenceable parameter
groups use the same dispatch and budgets. No independent unbounded collection
or duplicate codec was introduced. Existing whole-experiment writer preflight
validates metadata before writing bytes; I/O failure after output begins retains
the public writer's existing behavior.

Duplicate recognized activation metadata and duplicate intensity values fail.
A CV unit and `peak intensity unit accession` userParam cannot both supply the
unit, in either input order, even if identical or the CV uses the implicit
default. This refuses source overwrite semantics explicitly. Unit CV-prefix
conflicts, attributes without an accession, unknown intensity accessions,
non-string native intensity-unit metadata, and units attached to the unit
accession metadata itself also fail.

`tests/mzml_precursor_activation.rs` distinguishes source-reviewed expectations
from native roundtrip/guard tests. Metadata maps, scalar types, unit identities,
activation sets and literal values compare exactly; there is no relaxed numeric
tolerance or byte-for-byte mzML claim. Provenance is in
`tests/data/mzml_precursor_activation_provenance.json`. No C++ executable was run
for this closure, and no Claude review is claimed.

## Native use

```rust
use openms::metadata::{ActivationMethod, MetaValue, Unit};
use openms::Precursor;

let mut precursor = Precursor { intensity: 30.0, ..Default::default() };
precursor.cv_terms.metadata.insert(
    "peak intensity unit accession".into(), "MS:1000131".into());
precursor.cv_terms.metadata.insert("collision energy".into(),
    MetaValue::try_from(25.0)?.with_unit(
        Unit::new("UO:0000266", "electronvolt", "UO")?)?);
precursor.activation_methods.insert(ActivationMethod::Ethcd);
precursor.cv_terms.metadata.insert(
    "supplemental beam-type collision-induced dissociation".into(), "".into());
# Ok::<(), openms::Error>(())
```

Attach the precursor to a spectrum or chromatogram and use `mzml::write`; read
it back through `mzml::read`. The corresponding integration tests exercise both
record types and compare the complete precursor.

Focused remote verification on kim: 35 tests across activation, acquisition and
isolation; 10 activation tests on Rust 1.85; 51 mzML-enabled doctests; strict
rustdoc; and strict clippy for the library and these three test targets passed.
The broader feature-sliced clippy found existing dead-code helpers in
`tests/data_array_xml.rs`; the full-snapshot format check saw concurrent SQLite
edits. Final whole-tree results belong to the integration checkpoint. Retained
logs and exact tested runtime hashes are in
`tests/data/mzml_precursor_activation_validation/summary.json`.
