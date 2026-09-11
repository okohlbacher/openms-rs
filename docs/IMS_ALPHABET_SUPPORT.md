# IMS alphabets and text parsers

`IMSAlphabet`, `IMSAlphabetParser` and `IMSAlphabetTextParser` provide native
counterparts for the three OpenMS4-core headers at
`54a232fe2cae9c590d5c997fa49d20e7769860fb`. The alphabet retains complete owned
[IMS elements and isotope distributions](IMS_ISOTOPE_SUPPORT.md), including
independent sequence labels and hidden isotope bins. No new dependency is used.
This covers alphabet loading, container operations and verbose output; it does
not expose another mass-decomposition engine.

```rust
use openms::chemistry::{IMSAlphabet, IMSIsotopeOptions};
use std::io::Cursor;

let mut alphabet = IMSAlphabet::read(Cursor::new(
    "# name and monoisotopic mass\nH 1.0078250319\nO 15.9949146196\n",
))?;
alphabet.push_mass("C", 12.0)?;
alphabet.sort_by_mass()?;
let masses = alphabet.masses(0)?;
let hydrogen = alphabet.mass_by_name("H")?;
let description = alphabet.to_text(IMSIsotopeOptions {
    size: 2,
    abundances_sum_error: 0.0,
})?;
# Ok::<(), openms::Error>(())
```

## Complete native operation mapping

| Source alphabet operation | Rust counterpart |
| --- | --- |
| Empty/copy/container constructors | `new`, `Default`, `Clone`, `from_elements(Vec<IMSElement>)` |
| Container, iterator and scalar aliases | `Vec`, immutable `elements()` slice, its iterators, `usize`, `f64`, `str` |
| `size`, native empty query | `len`, `is_empty` |
| Indexed/name `getElement` | `element(index)`, `get(name)` |
| `getName` | `name(index)` |
| Indexed/name `getMass` | `mass(index)`, `mass_by_name(name)` |
| `getMasses(isotope_index=0)` | `masses(isotope_index)`; pass `0` explicitly |
| `getAverageMasses`, `hasName` | `average_masses`, checked `has_name` |
| `setElement(name,mass,forced=false)` | `set_element(name,mass,forced)`; pass `false` explicitly |
| Both `push_back` overloads | `push(owned_element)`, `push_mass(name,mass)` |
| `erase`, `clear` | Same names |
| `sortByNames`, `sortByValues` | `sort_by_names`, `sort_by_mass` |
| Default file `load` | `load(path)`; `read(BufRead)` also constructs an owned result |
| File `load` with replaceable parser | `load_with_parser(path, &mut dyn IMSAlphabetParser)` |
| Stream output | `to_text(options)`, `write(writer, options)` |
| Assignment, destruction | Normal Rust ownership |

`IMSAlphabetParser` is a small trait with `parse(&mut dyn BufRead)`, `elements()`
and `elements_mut()`, plus a default plain-file `load(&Path)`. Its concrete
`IMSAlphabetTextParser` implementation also provides `new`, `Default`, direct
`parse(impl BufRead)` and `load(impl AsRef<Path>)` methods. Parsed output is a
`BTreeMap<String, f64>`.

This trait retains the source's replaceable-format capability without C++
virtual-method ABI or template instantiations for arbitrary containers, numeric
types or stream subclasses. A custom Rust parser can maintain additional state
and expose a name-to-`f64` map at the alphabet boundary. Native `BufRead`, owned
maps and their normal mutation operations replace source container/stream aliases.

The source has **no flat alphabet store method**. Its output operator prints
verbose element descriptions, not the two-column input format. Accordingly,
`write` is verbose output and is not advertised as a parse round-trip codec.

## Container semantics

Construction and append preserve input order and allow duplicate names. Name
lookup, replacement and erase operate on the **first** matching element.
Replacement constructs a fresh single-mass element: it resets the old sequence
label to the requested name, sets nominal mass to zero and uses one unit-abundance
bin. Later duplicate entries remain intact. If a name is absent and `forced` is
false, the operation is a no-op; an unused nonfinite mass is not validated in
that branch.

Indexed masses use the complete element's raw isotope access, independent of the
isotope display-size setting. Average masses use all stored bins and preserve
source accumulation order. Missing names, out-of-range indices and unavailable
isotope bins return checked errors. Nominal-only elements with empty isotope
storage can be retained in the alphabet; their average mass is zero. Mass sorting
zero or one element is a source no-op and does not ask for a missing mass.

Name sorting uses lexical UTF-8 byte order. Mass sorting uses finite isotope-zero
masses. Equal names and equal masses retain prior order; signed zeros compare
equal. This deterministic stable order replaces the source `std::sort` tie order,
which is not specified. Sorting stages only indices and scalar mass keys, then
moves complete elements in place after every fallible operation has succeeded.
Sequence labels, isotope tails and ownership are preserved.

Loading builds single-mass elements from the parser's lexical map order and then
sorts by mass. Thus equal-mass names loaded through a map retain lexical order.
An ordinary container constructed or appended directly retains duplicates; the
text parser's map does not.

## Built-in text format

The reader consumes LF-delimited lines, retaining a CR before LF for stream
whitespace processing. Its comment precheck matches the source exactly: skip
only leading **space and tab**, then ignore an empty line or a line beginning
with `#`. Name/mass token extraction uses the six classic ASCII whitespace
characters: space, tab, LF, CR, vertical tab and form feed. Rust's narrower
`is_ascii_whitespace` predicate is deliberately not used.

Each remaining line supplies an unquoted name followed by a decimal/scientific
mass prefix. Optional signs, leading/trailing decimal points with at least one
digit, and complete `e`/`E` exponents are accepted. Trailing columns, comments or
other text are ignored after the numeric prefix. The first duplicate name wins;
parsing a new input replaces the previous map. Finite signed numbers are allowed,
and decimal underflow may produce signed zero. Alphabet element mass access
subsequently retains its source addition of nominal zero, which can change a
negative zero's sign.

Malformed rows, missing numbers, incomplete exponents, overflow, nonfinite values
and invalid UTF-8 produce errors with one-based line numbers. The previous built-in
parser map and an existing alphabet remain unchanged. These are intentional
corrections to source parsing that clears first, ignores stream fail state and
may use an uninitialized, stale or implementation-dependent mass after failure.
A malformed duplicate row still errors, even though an already-parsed value
would win the map insertion.

The source delegates numeric extraction to the active C++ stream implementation
and locale. Hexadecimal floating syntax, locale decimal separators/grouping and
precise `num_get` consumption of adjacent letters vary. The built-in native parser
uses a portable decimal/scientific policy instead of claiming byte-identical
behavior across C++ libraries. Hexadecimal `0x`/`0X` prefixes are explicitly
rejected and can be implemented by a custom parser. Under the native prefix
policy `1.25suffix` gives `1.25`; some C++ libraries consume additional hexadecimal
letters and then signal failure for similar inputs. A CR-only or vertical-tab-only
line that passes the narrower comment precheck but has no name/mass is a checked
malformed line, replacing the source's failed extraction behavior.

## Custom parsers, atomicity and resource limits

The default parser and default alphabet load share a budget across line reads,
map construction, owned element construction, sorting and replacement. A custom
parser owns its computation and side effects, including what it does with its
input stream. The alphabet revalidates the entire returned map before committing:
count, name lengths, logical payload and every mass must be valid. Custom parser
state is not rolled back when later alphabet validation fails. `elements_mut`
retains ordinary Rust mutation semantics; returning a large or invalid map does
not bypass alphabet validation. A built-in parser also meters its previous map
before replacing it; direct edits can be cleared through the mutable map if they
exceed these limits.

The native limits are:

- 100,000 alphabet elements or parser map entries.
- 1 MiB per name or sequence (from `IMSElement`), and per input line including LF.
- 64 MiB input and 64 MiB cumulative logical allocation accounting per default
  operation. Retained alphabet payload is also limited to 64 MiB, including full
  isotope vectors, labels and element slots.
- 50,000,000 work units, including input bytes, map-key comparisons, name sorting,
  copies and replacement destruction. Prefix comparisons are charged by compared
  byte-length bounds. Fallible stable sorting measures actual comparison calls.
- 8 MiB verbose output, conservatively checked before formatting. Shared output
  work/allocation includes all nested element/isotope scalar formatting, temporary
  strings and the final combined string, so a resource limit may be reached before
  the destination alone reaches 8 MiB.

Map insertion uses a conservative node allowance and comparison bound rather
than claiming a precise allocator-specific memory ceiling. Incoming owned values
are moved; their preexisting allocation is not attributed to the current call.
`Clone`, equality, direct parser-map mutation and `clear` use normal Rust ownership
and destruction. Computational and I/O entry points return checked errors before
committing modified alphabet state. `write` fully validates/formats before its
first byte; an underlying writer failure can still leave a partial external
stream, as with ordinary Rust `Write`.

## Evidence

The [manifest](../tests/data/ims_alphabet_provenance.json) records the three source
headers, two implementations and three class tests separately from native
policies. Literal source tests supply hydrogen/oxygen/nitrogen masses `1/16/14`,
carbon `12`, peptide letters `A/R` with masses `71.03711/156.10111`, named/indexed
access, replacement/erase/sort cases and the custom-parser extension example.

Twelve native tests cover those assertions plus duplicate handling, sequence
replacement, deterministic ties, multiple isotopes, decimal prefixes and
underflow, vertical-tab extraction, malformed/late-I/O rollback, custom mutated
maps, exact verbose formatting and cumulative comparison/storage limits. No C++
execution is claimed. The previously frozen isotope/element implementation is a
stage dependency and is not part of this increment's extraction list.
