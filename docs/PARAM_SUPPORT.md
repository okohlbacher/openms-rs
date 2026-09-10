# Parameter trees

The native `param` module implements the retained Core SDK `Param`, `ParamEntry`, `ParamNode`, and forward parameter traversal. `ParamValue` is documented separately in [PARAM_VALUE_SUPPORT.md](PARAM_VALUE_SUPPORT.md). The implementation uses owned Rust strings, vectors, and sets; it does not depend on the C++ library.

## Public API mapping

| Source surface | Native API |
| --- | --- |
| Entry construction, public data, `isValid` | `ParamEntry::new`, public fields, `is_valid` / `validation_error` |
| Node construction, local/recursive lookup, `suffix`, `size`, both insertions | `ParamNode::new`, `find_entry`, `find_node`, `find_parent_of`, `find_entry_recursive`, `suffix`, `size`, `insert_entry`, `insert_node` |
| Tree construction/copy/move/assignment | `Param::new`, `Default`, `from_root`, `Clone`, `checked_clone`, Rust moves/assignment |
| `getEntry`, value/type/description getters, `exists`, `hasSection` | `entry`, `value`, `value_type`, `description`, `exists`, `has_section` |
| `setValue`, all tag and restriction methods | `set_value`, `tags`, `has_tag`, `add_tag`, `add_tags`, `clear_tags`, `valid_strings`, `set_valid_strings`, `set_min_int`, `set_max_int`, `set_min_float`, `set_max_float` |
| Section add, descriptions, tree size/empty/clear | `add_section`, `section_description`, `set_section_description`, `size`, `is_empty`, `clear` |
| Whole-tree insertion/removal/copy | `insert`, `insert_entry`, `remove`, `remove_all`, `copy`, `copy_subset`, `copy_subset_with_messages` |
| Defaults, checks, all update overloads, merge | `set_defaults`, `check_defaults`, `update`, `update_with_options`, `merge` |
| Both command-line parsers | `parse_command_line`, `parse_command_line_mapped` and `CommandLineOptions` |
| Iterator begin/end, copying/increment, dereference, name/trace | `iter`, `ParamIterator::new` / `Default`, `Clone`, standard fused `Iterator`, `ParamItem::{entry,key,trace}`, `end_trace` |
| `findFirst` / `findNext` | `find_first`, `find_next`; the latter accepts a borrowed prior entry |
| Stream output | `to_text`, `Display`, default classic-locale source value formatting |
| Source equality predicates | Explicit checked `source_equal`; ordinary Rust `PartialEq` compares complete values and order |

Entry fields retain all restrictions independently of the current value type. Default bounds are `-f64::MAX` / `f64::MAX` and `-i32::MAX` / `i32::MAX`; an unchanged bound is an inactive sentinel. NaN and infinities remain representable, and restriction comparisons follow the source. Scalar integer restriction checks narrow to `i32` with an error on overflow. Strings and string lists are type-strict. File tags bypass string restrictions; `output prefix` bypasses scalar strings only.

Entries precede child sections during iteration, preserving insertion order within each vector. Transitions describe opened and closed sections, suppress empty section pairs, and retain the final close transitions through `end_trace`. Each yielded row is valid independently; exhaustion yields `None` repeatedly. A foreign entry passed to `find_next` is rejected. Root-level entries are excluded by the historical `:leaf` suffix search.

## Preserved source behavior

- Inserting an existing entry replaces its value and tags, retains its restrictions, and replaces its description only when the old description is empty or the new one is nonempty. New entries retain the supplied restrictions.
- Intermediate sections may have the same name as a leaf. Only a collision at the final insertion target is rejected. Empty path components and colon-containing draft names are retained. Setter tags may contain commas, while `add_tag`, `add_tags`, and `set_valid_strings` reject them.
- `has_section` uses the source prefix-parent predicate, so a matching leaf or partial section name can return true. It strips one trailing colon. Section-description lookup does not strip that colon. An empty section query, undefined in the source, is a checked error.
- `copy` without a trailing colon selects local names by prefix. Removing the prefix can produce empty names. Copying an empty section with a trailing colon returns an empty tree because of the source parent lookup. `copy_subset` selects root entries by name and whole root sections; absent names can be returned as diagnostics.
- Deletion prunes only affected empty ancestors. Prefix deletion uses a linear stable `Vec::retain` implementation with the same successful result as the source's repeated vector erasure.
- Defaults preserve existing values and copy only restrictions relevant to each new value. Source prefix/section-description lookup expressions are retained, including the distinction between normalized and original prefixes.
- `check_defaults` preserves the source's two different prefixed lookup expressions. Unknown entries produce returned warnings; incompatible types or violated applicable restrictions produce errors.
- Updates preserve `:version` values and `:type` values having at least two colons. A missing full key can map to exactly one matching nested leaf. Equal values are not revalidated. `fail_on_unknown_parameters` takes precedence over adding unknown entries. A report with `success == false` retains successful changes, as in the source; resource and conversion errors roll back the entire operation.
- Merge adds missing values and retains existing entry payloads. Its original section-description lookup can replace an existing description when the duplicated prefix lookup is empty.
- Command-line input includes the executable at index zero, which is skipped. An option begins with `-` and has a nondigit second byte; negative numbers such as `-1.0` are values, while `-.3` is an option. Mapped parsing prioritizes multiple-argument, then flag, then one-argument registrations. Flags store the string `true`; unknown options and plain arguments accumulate separately. Repeated recognized options replace their previous values.

Source equality ignores descriptions and restrictions, and node equality ignores order. Its membership-based comparison can even be asymmetric for duplicate draft entries. That behavior is available only through `source_equal`; native `PartialEq` deliberately compares full values and preserves ordering. Floating-point NaN keeps its normal unequal-to-self semantics, so these types do not implement `Eq`.

## Checked boundaries

Checked tree operations cap depth at 128, live entries and nodes at 100,000 each, logical payload/allocation at 64 MiB, and cumulative work at 50 million units. Work includes text scans, lookups, comparison input, copied payloads, iteration paths, formatting, and vector moves. Nested defaults, updates, merge, and command-line parsing share one operation budget. A failed operation leaves the destination unchanged. Error and warning text conveys the source condition rather than reproducing the C++ logging stream byte for byte.

Ordinary atomic mutations use one bounded tree snapshot. Repeating them on a large tree can be quadratic. The XML adapter instead uses a crate-private owned builder with one shared budget and final validation; a failed builder cannot publish partial state. The insertion-order vectors still use linear sibling lookup, so the work budget can bind before the count limits for broad trees. No numeric restriction, logging flag, or unused field is silently normalized.

Public `ParamNode` is an owned draft. Checked construction validates depth before recursive processing, and its destructor drains children iteratively even for a rejected deep draft. Ordinary Rust `Clone`, direct field edits, and comparisons on caller-created drafts follow Rust collection conventions; use checked tree operations for resource-limited input. Rust strings represent UTF-8 text, including embedded NULs, rather than arbitrary invalid-UTF-8 C++ byte strings.

## Evidence

[tests/param.rs](../tests/param.rs) contains source-derived restrictions, insertion, query, iterator trace, copy/remove, defaults, update, and both command-line parser cases, plus source quirks and transactional boundary cases. Private module tests pin empty-name work charging, allocation failure, and shared builder failure. [param_provenance.json](../tests/data/param_provenance.json) pins the authoritative source files and test source. XML behavior and its original source fixture are covered separately by [PARAMXML_SUPPORT.md](PARAMXML_SUPPORT.md).
