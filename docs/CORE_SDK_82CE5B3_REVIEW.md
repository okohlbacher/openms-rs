# Core SDK refresh to 82ce5b3

The target advances from `54a232fe2cae9c590d5c997fa49d20e7769860fb` to
[`82ce5b373c97f934ffd9b1ffd80215ca66473d0b`](https://github.com/okohlbacher/OpenMS4-core/tree/82ce5b373c97f934ffd9b1ffd80215ca66473d0b),
verified against upstream HEAD and fetched into a separate clean checkout on
2026-09-11. The previous checkout and original fixture pins are retained.

The [upstream comparison](https://github.com/okohlbacher/OpenMS4-core/compare/54a232fe2cae9c590d5c997fa49d20e7769860fb...82ce5b373c97f934ffd9b1ffd80215ca66473d0b)
contains three commits and thirteen changed paths:

- `74526a8` shares numeric formatting between Core and its C++ test framework.
- `c6adde1` fixes Windows process/executable-path probes and their tests.
- `82ce5b3` strengthens Windows SDK acceptance diagnostics and reproduction.

## Runtime compatibility

The only changed file within the comparable scientific roots is
`src/openms/source/DATASTRUCTURES/StringUtils.cpp`. Its private `appendNumeric`
function moves to `src/common/include/OpenMS/CONCEPT/Detail/NumericFormatting.h`;
six callers now qualify the shared namespace. The complete function body,
including comments, is byte-identical after stripping indentation on each line.
Its normalized SHA-256 is
`1228c1db508dc1c02942fad91bb8c0e7c5cb8b8abf8064fc15cdbc60da5773d5`.

Float/double thresholds, precision, special values, trimming and exponent
formatting do not change. Native parameter/list formatting needs no numerical
change. The existing extended-long-double boundary remains explicit. The test
framework now uses this same formatter; no scientific class-test literal or
dataset changes. Core and test-framework build files add the private include
path. Remaining changes affect CI, runtime probes and installed acceptance tests,
including checking Arrow results rather than silently discarding them.

All source headers/implementations used by the next MassTrace, MonosaccharideDB
and mzML settings work remain unchanged. This update does not add or remove a
registered public scientific API or require deleting native implementations.

## Inventory and evidence

The comparable inventory still contains 1,578 scientific files, including
807 physical include-directory headers and 786 registered public headers.
Its physical line count is 468,228. The new 141-line shared private formatter is
outside those inventory roots and is explicitly tracked in
[the refresh manifest](../tests/data/sdk_numeric_formatting_provenance.json);
it is not counted as a new public SDK header. This avoids conflating a source
relocation with API removal or a scientific rewrite.

The manifest hashes all thirteen changed paths. The generated inventory refreshes
current file/registration hashes and preserves the original historical reference
hashes alongside explicit reviews. Per-fixture `target_verification` records the
new target without repinning the source from which expected values were obtained.
The prior [54a232f review](CORE_SDK_54A232F_REVIEW.md) remains applicable to earlier
changes. No new C++ execution or full SDK build is claimed by this refresh.
