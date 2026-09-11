# Contact and chromatography metadata

Native `metadata::ContactPerson`, `HPLC` and `Gradient` cover their complete
class-specific public operations at Core SDK `82ce5b3`. These are the first
additional value records needed by `ExperimentalSettings`; this group does not
yet add that aggregate to MSExperiment or serialize a complete mzML header.

ContactPerson owns all seven source strings and typed metadata. Public fields,
Clone and equality replace ordinary source accessors and copy operations.
`set_name()` retains the source's literal split rules: commas take precedence,
the first two comma-separated fields become last and first names after trimming
space/tab/newline/carriage-return, and further fields are ignored. Without a
comma, only the first two literal space-separated fields are used, including
empty fields and without trimming. Without either delimiter, only last_name
changes; first_name is retained even when the new name is empty. Unicode strings
remain intact. The operation limits input to 4 MiB and stages copied names before
publication. Direct field ownership and metadata follow ordinary Rust costs.

HPLC owns instrument, column and comment strings, temperature in degrees Celsius,
pressure in bar, flux in microliters per second, and an owned Gradient. Temperature
defaults to 21; pressure and flux default to zero. Native u32 fields preserve the
public source unsigned domain without its signed private-storage conversion.
All fields participate in equality and cloning; changing them does not impose
additional physical validation absent from the source record.

Gradient keeps ordered unique eluent names, strictly increasing signed i32
timepoints, and row-major u32 percentages. Names are exact and case-sensitive;
empty names are legal. Adding an eluent appends a zero-filled row; adding a time
appends a zero to each active row. Existing timepoints cannot be repeated, and
negative times are accepted when strictly increasing. Checked percentage queries
require known names/times; updates additionally require values no greater than 100.
`is_valid()` checks that every active timepoint totals 100. No active timepoints
is valid, while active timepoints with no eluents is invalid.

The source's clear-axis behavior is deliberately preserved: `clear_eluents()`
removes only names and `clear_timepoints()` removes only times. Percentage rows
and values remain visible, and subsequent additions can address these stale
values under new names/times. `clear_percentages()` is the explicit operation
that rebuilds the table to the current axes with zeros. If stale storage would
cause an out-of-bounds source access, native queries and updates return an error.
The percentage getter exposes an immutable view, so callers cannot create
arbitrary inconsistent storage directly.

Gradient allocation operations preflight cumulative axes, row descriptors and
stored cells against one million items and a 64 MiB conservative vector/string capacity
allowance, with fourfold growth headroom. Empty-axis retained capacities, stale cells and the
old table during reset count too; a small final table can therefore still exceed
the temporary allowance. Allocation occurs before publishing changed values.
Failed reserves can change private vector capacity, while all observable values
remain unchanged. Append uses geometric growth. Source undefined integer/index
behavior becomes checked errors; normal Clone, direct HPLC/contact assignment
and destruction retain ordinary Rust costs. No exact allocator-memory bound is
claimed.

[Tests](../tests/experiment_metadata.rs) retain source contact/table/default
literals, cover stale storage and repair, signed and unsigned boundaries,
resource failures and all 10,201 two-eluent percentage pairs. Source headers and
implementations were independently reviewed; no C++ execution is claimed.
[Provenance](../tests/data/experiment_metadata_provenance.json) pins the records,
class tests and string helper. Other experiment-header records, the complete
ExperimentalSettings aggregate and its XML/consumer transport remain open.
