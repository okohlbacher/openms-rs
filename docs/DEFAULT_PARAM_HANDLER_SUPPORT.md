# Default parameter handlers

`openms::param::DefaultParamHandler` provides the retained Core SDK's parameter/default trees, diagnostic name, delegated subsections, validation switches, and member-update lifecycle using the native [Param](PARAM_SUPPORT.md) API. It owns its configuration and borrows no external tree. No inheritance, C++ ABI, or stored callback is required.

## Construction and updates

`new(name)` creates empty current/default trees and subsections, with `check_defaults` and `warn_empty_defaults` enabled. Read accessors expose these six fields. `set_name`, `set_defaults`, and `set_subsections` configure them; the two Boolean switches have direct setters. Installing defaults does not change current parameters. Complete initialization with `defaults_to_parameters` or its callback variant.

`set_parameters(&Param)` starts from the supplied tree and fills missing default entries. It does not merge with the previous current tree: an omitted option returns to its default, and an omitted nondefault option disappears. Default checking retains `Param`'s source rules: unknown entries yield warnings; mismatched types and violated applicable restrictions yield errors. Warnings are returned as `Vec<String>` rather than sent to global logging streams. When checking is disabled, completion of missing defaults still occurs. The empty-default warning requires both switches to be enabled.

Each registered subsection is excluded from the validation copy, while its values remain in the committed parameter tree and the callback sees them. Values are retained verbatim, including duplicate/empty strings and trailing colons. As in the source, exclusion appends a colon; normally supply `"algorithm"`, not `"algorithm:"`. Entries such as `algorithm_extra:option` are not excluded by `algorithm`.

`defaults_to_parameters` fills missing current values and preserves existing values. It does not check restrictions, including restrictions on the default values themselves. It reports only the **first** missing default description in source iteration order; the source loop stops after that entry despite the plural warning text. These behaviors follow the implementation rather than a header comment implying that initialization simply copies all defaults.

## Keeping typed settings synchronized

`set_parameters_with` and `defaults_to_parameters_with` take a callback `FnOnce(&Param) -> Result<T>` and return `(T, Vec<String>)`. The callback constructs new owned typed settings from the complete staged tree. Validation and the callback must succeed before the handler commits. Install the returned settings after success:

```rust
use openms::param::{DefaultParamHandler, Param, ParamValue};

struct Settings { count: i32 }

fn settings(parameters: &Param) -> openms::Result<Settings> {
    Ok(Settings { count: parameters.value("count")?.to_i32()? })
}

fn main() -> openms::Result<()> {
    let mut defaults = Param::new();
    defaults.set_value("count", ParamValue::from(2), "Number of items", &[])?;
    defaults.set_min_int("count", 1)?;
    let mut handler = DefaultParamHandler::new("Example")?;
    handler.set_defaults(defaults)?;
    let (mut typed, _) = handler.defaults_to_parameters_with(settings)?;
    assert_eq!(typed.count, 2);

    let mut supplied = Param::new();
    supplied.set_value("count", ParamValue::from(3), "", &[])?;
    let (replacement, warnings) = handler.set_parameters_with(&supplied, settings)?;
    typed = replacement;
    assert_eq!(typed.count, 3);
    assert!(warnings.is_empty());
    Ok(())
}
```

The callback must not mutate existing external state: arbitrary callback side effects cannot be rolled back. Its own application work is outside the handler's resource budget. The noncallback methods use a no-op callback. A failed validation never invokes the callback, and a returned callback error leaves the handler unchanged. This deliberately corrects source `setParameters`, which assigns its current tree before checking restrictions and can leave member variables out of sync after an exception. Tests demonstrate both failure paths with an algorithm owning a handler and typed settings.

`checked_clone` copies the owned configuration with resource checks. Normal Rust `Clone` and moves replace the source copy/assignment operations; they do not invoke a derived callback. `source_equal` compares all handler policy fields and uses the historical name/value-only `Param` equality. Native `PartialEq` instead compares the complete stored trees, including their descriptions and restrictions.

## Metadata conversion and limits

`write_parameters_to_meta_values(parameters, metadata, prefix)` converts all seven [ParamValue](PARAM_VALUE_SUPPORT.md) alternatives to native metadata, widening integer lists from i32 to i64. A nonempty prefix receives a final colon if needed. The source uses each entry's **leaf name**, not its full parameter path: `a:option` and `b:option` both become `option`, or `prefix:option`. The last visited value wins, following source entry/section iteration order. Existing unrelated metadata, including its units, remains intact. No units are invented for replacement values.

Native [MetaValue](METADATA_SUPPORT.md) excludes nonfinite floating-point values. Exporting NaN or infinity therefore returns an error even though parameter storage itself permits them. Metadata conversion, all updates, and checked copies complete fallible work before committing; errors preserve the original destination.

The destination is the `ParameterMetaSink` trait rather than the `MetaInfo` type, and `MetaInfo` implements it, so that `param` names nothing in `metadata`. `param` is the lower module and `metadata` is the one that names it; a `param` that named `MetaInfo` would close a module cycle between the two, and before `metadata -> chemistry` was cut it closed the longer `param -> metadata -> chemistry -> param`. The conversion and the walk of the existing destination are therefore written in `metadata`, charged against the handler's own work and allocation limits through `ParamBudget`. Calls are unchanged — the destination type is inferred — and the behaviour above is the sink's contract, not an implementation detail of one destination.

Operations share `Param`'s 50-million-unit work and 64-MiB logical allocation ceilings across tree scans, default completion, validation, copies, warnings, subsection exclusion and metadata preparation. Subsection configuration additionally permits at most 100,000 strings. Parameter-tree depth/node and value-list limits also apply. Existing metadata contributes to the export payload budget; key comparison work includes long common prefixes. Logical accounting is conservative and includes transient copies, so these bounds do not promise that every object individually below 64 MiB can be processed. Borrowed accessors, Boolean setters, and ordinary Rust `Clone` retain their normal behavior.

## Evidence

[default_param_handler.rs](../tests/default_param_handler.rs) contains nine focused tests covering the source class-test literals, default and subsection behavior, warning switches, explicit source equality, typed callbacks, metadata collisions and all value alternatives, error atomicity, and resource rejection. [default_param_handler_provenance.json](../tests/data/default_param_handler_provenance.json) records six exact SDK source hashes and the source assertion/behavior locations. Additional native tests establish checked behavior without executing or building C++.
