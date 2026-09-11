# MassDecomposition values

`chemistry::MassDecomposition` ports the public count-container API from `CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h`. It represents symbol frequencies; it does not find compositions with a requested mass. The separate [native decomposition algorithm](MASS_DECOMPOSITION_ALGORITHM_SUPPORT.md) performs mass searches; it remains a distinct class.

```rust
use openms::chemistry::MassDecomposition;

fn main() -> openms::Result<()> {
    let mut composition = MassDecomposition::parse("C3 M4")?;
    composition.checked_add_assign(&MassDecomposition::parse("S2")?)?;
    assert_eq!(composition.to_text()?, "C3 M4 S2");
    assert_eq!(composition.to_expanded_string()?, "CCCMMMMSS");
    assert!(composition.contains_tag("SMC")?);
    assert!(composition.compatible(&MassDecomposition::parse("C2 S1")?));
    Ok(())
}
```

| Source API | Native API |
|---|---|
| Default constructor, copy, assignment | `new`, `Default`, `Clone`, ordinary owned assignment |
| String constructor | `parse`, `FromStr` |
| `operator+` | `checked_add` |
| `operator+=` | `checked_add_assign` |
| `toString` | `to_text` |
| `toExpandedString` | `to_expanded_string` |
| `getNumberOfMaxAA` | `number_of_max_aa` |
| `operator<` | `source_cmp`, `source_less` |
| Equality against a string | `equals_text` |
| `containsTag`, `compatible` | `contains_tag`, `compatible` |

## Source parsing and value identity

Tokens are separated by the literal ASCII space, without quote protection or general whitespace splitting. A token consists of one ASCII byte symbol followed by a signed-i32 decimal spelling. Count parsing reuses the source-derived native `ListParse` conversion, including trimming space/tab/CR/LF around the count and removing one leading `+`. Counts must be nonnegative. The empty string is an empty composition. Empty tokens, missing counts and numeric overflow are errors. Thus leading, trailing or repeated spaces ordinarily fail.

Everything from the first `(` onward is ignored. Only when this suffix is present does the source trim the remaining complete input before splitting. For example, `" C3 M4 (notes)"` is accepted while `" C3 M4"` fails. Ignored text still counts against the input byte limit. No balanced-parenthesis parser or amino-acid alphabet validation is added. Lowercase letters, digits, punctuation and ASCII control-byte symbols retain source behavior when they are syntactically representable. Non-ASCII symbols are rejected instead of introducing platform-dependent signed-char ordering or an arbitrary byte-string output API.

A duplicate symbol replaces its map count but does not reduce the separately cached maximum. Therefore `A9 A1` formats as `A1`, expands as `A`, and reports maximum 9. `equals_text("A9 A1")` is true, while `equals_text("A1")` is false. Rust `Eq` and `Hash` include both the count map and the maximum. Source comparison ignores the maximum and compares sorted symbol/count pairs lexicographically; this is exposed explicitly through `source_cmp` rather than an inconsistent Rust `Ord` implementation. Counts compare numerically, so `A2` precedes `A10`.

Both addition forms start with the left operand's cache and ignore the right operand's cache. `checked_add_assign` raises the running maximum from resulting counts. Source `operator+` has a further finite quirk: for each newly introduced key, it compares the count against the **original left** maximum, then assigns that count to the result cache. A later new key can lower an earlier maximum. Thus `A1 + B10 C5` reports maximum 5, while applying `+=` reports 10. An earlier existing-key sum can be lowered too: `A10 + A100 B11` reports 11 despite containing A110. `checked_add` preserves this source behavior. Addition order and choice of operator can therefore affect cached identity even when result maps agree. Copying a value preserves its entire state.

Zero counts retain their map keys. `Z0` expands to an empty string but formats as `Z0`; an empty composition is not compatible with `Z0`, because source compatibility requires the key to exist. An empty tag is always contained. Tags are multisets: their order is irrelevant, and spaces or modification spellings are not special syntax. `compatible` checks whether every key/count of its argument is present with sufficient count. Its false result does not print the source's unsolicited stderr diagnostic.

Both formatters use ascending ASCII byte order. Compact formatting includes zero counts and applies the source's final whitespace trim. An unusual leading tab/newline symbol can therefore disappear from compact text while remaining in expanded output. Compact text is not guaranteed to reconstruct duplicate-history or operator-specific cached maxima, or counts increased beyond the constructor's i32 range. These are source conventions, not a canonical serialization promise.

## Checked native boundaries

The input and tag cap is one MiB per operation. Expanded output is capped at 16 MiB. Formatting sums lengths with checked arithmetic before allocation or repetition. At most 128 ASCII map keys exist, so compact output, comparison, copying, compatibility and addition have small fixed structural bounds. `parse` and `contains_tag` use a 50-million-unit work budget, charging input traversal and a conservative full-map comparison allowance before each lookup/insertion. The tag budget can reject a large multiset test even below the byte limit.

Stored counts use `usize`. The source converts negative i32 counts to unsigned `Size`; this native API rejects them. Source unsigned addition wraps on overflow; native addition reports an error. `checked_add_assign` verifies every sum before changing any map entry or the cached maximum, including when a later key fails. No counts are silently saturated. Parsing and value-producing operations return owned results, so a failure cannot partially replace a caller's value.

## Verification

[Ten focused tests](../tests/mass_decomposition.rs) cover all distinct literal assertions in the upstream class test, plus duplicate history, both addition directions and operator-specific maximum behavior, map ordering, zero-key behavior, arbitrary ASCII symbols, suffix/whitespace grammar, i32 endpoints, checked usize overflow, atomic later-key failure, source formatting limitations, and input/work/output bounds. Expected literals are transcribed from the source; boundary expectations are independently derived from its operations and the stated native corrections.

[Provenance](../tests/data/mass_decomposition_provenance.json) records the exact header, implementation, class test and StringUtils dependencies at SDK revision `54a232fe2cae9c590d5c997fa49d20e7769860fb`. No C++ build or execution was used. The module has no new runtime dependencies and uses an owned standard-library `BTreeMap`.
