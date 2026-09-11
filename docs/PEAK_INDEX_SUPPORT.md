# Peak and feature indices

`kernel::PeakIndex` represents the full public `PeakIndex` value/access surface
at Core SDK `54a232fe2cae9c590d5c997fa49d20e7769860fb`. It has public `peak` and
`spectrum` indices, ordinary copy/equality/hash traits, and the source invalid
sentinel `usize::MAX` in each default component. It owns no data or allocation.

`new(spectrum, peak)` follows the source argument order. `for_feature(peak)`
leaves the spectrum unset. `is_valid` checks only that the peak is not the
sentinel; it does not verify bounds and an unset spectrum does not make a
feature index invalid. `clear` resets both components.

`get_feature(&[T])` returns a borrowed feature-like value using the peak index,
ignoring spectrum. Pass `&map.features` for a native feature or consensus map.
`get_spectrum(&MSExperiment)` checks spectrum independently of peak validity.
`get_peak(&MSExperiment)` checks spectrum first, then the selected peak. These
accessors always return checked errors for invalid or out-of-bounds components;
the C++ preconditions are only enforced in debug builds. No pointer ownership,
copying of scientific values, or lifetime extension is introduced.

Indices can be constructed from an area iterator row's spectrum/peak indices.
An exhausted iterator has no row; `PeakIndex::default()` is the source-compatible
invalid value where an explicit sentinel is needed. No implicit remapping occurs
when a caller reorders a container.

[Tests](../tests/peak_index.rs) cover the source constructor and feature/peak-map
literals, borrowed reference identity, unset components, empty spectra, both
dimension boundaries and all 36 combinations of a small independent slice-access
oracle. [Source hashes](../tests/data/peak_index_provenance.json) identify the
header, empty implementation unit and class test. No C++ execution is claimed.
