# PepXMLFileMascot

`FORMAT/PepXMLFileMascot.h` (core `bc9cc12`), ported by lane F3a in wave 4 as
`src/format/pepxml_mascot.rs` (feature `idxml`). It reads the pepXML Mascot
exports into peptide sequences with their modifications, keyed by spectrum
title. No TOPP tool calls it (`docs/core-sdk-coverage.json` lists no direct
consumer); full pepXML identifications are `PepXMLFile`
(`src/format/pepxml.rs`).

## API mapping

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `PepXMLFileMascot()` | `PepXmlFileMascot::new`, `Default` | stateless |
| `void load(filename, std::map<std::string, std::vector<AASequence>>& peptides)` | `PepXmlFileMascot::load`, `load_with_limits`, `read`, `read_with_limits` | returns `MascotPepXmlDocument`: `peptides` (a `BTreeMap`, the `std::map`'s byte order) and the handler's log |
| `@exception Exception::FileNotFound` | `Error::Io` | from opening the file |
| `@exception Exception::ParseError` | `Error::Parse` | malformed XML, a missing required attribute, an unparsable `peptide` |
| protected `onStartElement`, `onEndElement` | private `Handler::start`, `Handler::end` | driven by `format::sax_xml` with `SaxEncoding::Declared` (Xerces honours the declared encoding; the source enforces none) |
| protected `matchModification_(mass, description)` | private `Handler::match_modification` | first variable definition with `|mass difference| < 0.001`; `MODIFICATION_TOLERANCE` |
| protected members `actual_title_`, `actual_sequence_`, `actual_modifications_`, `peptides_`, `actual_aa_sequences_`, `fixed_modifications_`, `variable_modifications_` | fields of the private `Handler` | per read, never shared between reads |
| (none) | `ReadLimits`, `MascotPepXmlDocument::log`, `suppressed_log_lines` | resource ceilings; the source's `OPENMS_LOG_ERROR` output as data |

Every public and protected member of the header is accounted for above; the
class has no other API.

## Preserved source conventions

- Element names are matched exactly as written (no namespace resolution), as
  Xerces delivers them with namespace processing off.
- `<aminoacid_modification>` and `<terminal_modification>` with
  `variable="Y"` (exactly; `y` is fixed) are variable definitions (description
  and mass), any other value fixed (description only). Both read `variable`
  first; for a variable one the Release build (GCC, arguments right to left)
  reads `mass` before `description`, so with both missing the error names
  `mass` (oracle case `m14`).
- A `<mod_aminoacid_mass>` records its matched description (or an empty one)
  and its `position`, parsed as Xerces' `parseInt` and converted to unsigned:
  `-1` becomes 4294967295, and position 0 wraps to index 4294967295, both the
  source's `IndexOverflow` (`m08`, `m10`). Modifications collected outside a
  `<search_hit>` apply to the next hit (`m07`).
- At the end of a `<search_hit>` the `peptide` is parsed with
  `AASequence::fromString`'s permissive default, then each recorded
  modification applied in order: `C-term` or `N-term` in the description sets
  the whole description as the terminal modification; a description of exactly
  two words (`StringUtils::split` on one space) also sets its first word on the
  residue at the position (`m02`: `Acetyl (N-term)` gives
  `.(Acetyl)K(Acetyl)PEPTIDE`); any other description, including the empty one
  of an unmatched mass, is the non-fatal error
  `Cannot parse modification '<description>@<position>'` and the hit is kept
  (`m06`, `m16`, `m26`). Then every fixed modification of two words: `(C-term)`
  and `(N-term)` set the terminus, any other origin (parentheses removed) sets
  the name on every residue whose letter occurs in it, overwriting a variable
  modification there; other shapes log `Cannot parse fixed modification '...'`
  (`m04`, `m05`, `m22`).
- A registry miss (`Foo`, or `Amidated` on a residue) is the source's
  `InvalidValue` exception and ends the load (`m03`, `m11`, `m23`).
- At the end of a `<spectrum_query>` its hits are inserted under its title; a
  title seen before keeps its first entry (`std::map::insert`), and hits that
  occurred outside any query join the next one (`m07`).
- Log lines are `Non-fatal error while loading '<file>': <message>`, with the
  file name as passed (empty for `read`).

## Native differences

- **No state survives a read.** The source resets its members only after a
  successful parse, so a load after one that threw still matches masses
  against the failed file's variable modifications: the oracle's `m28_reuse`
  loads a throwing file defining `Oxidation (M)` and then a file defining only
  `Dioxidation (M)`, and gets `M(Oxidation)PEPTIDE`; a fresh object gets
  `MPEPTIDE` and the log line (`m28_fresh`). Every read here is a fresh object.
  C++ defect candidate (unnumbered): stale handler state after an exception.
- **A failed read returns nothing.** The source's out-parameter keeps the
  queries completed before the exception; here an error carries no partial
  map. (Every recorded failing case failed in its first query.)
- **Empty peptide with a terminal modification** (`m19`): the source stores a
  sequence that is only a terminal modification (`.(Acetyl)`); the crate's
  `AASequence` has no terminus without residues, so this is
  `Error::InvalidValue`.
- **Permissive parsing is emulated here.** The crate's `AASequence::parse` has
  no permissive mode, so outside round and square brackets `*`, `#` and `+`
  become `X` and spaces are dropped before parsing, as the source's
  `parseString_` does (`m20`). This belongs in `AASequence` itself; it is
  local until that port grows the option.
- **Printing a protein-terminal modification** is an `AASequence`
  difference, not this reader's: the source's `toString` prints every
  modification by its id (`ResidueModification.cpp:959-972`), so `m26`'s
  sequence is `.(Acetyl)SPEPTIDE`, while the crate's `Display` prints a
  protein-terminal one by its full id, `.(Acetyl (Protein N-term))SPEPTIDE`.
  The resolved modification is the same record. The test compares every
  sequence through the source's rule and checks `Display` against it for every
  sequence without a protein-terminal modification.
- **Error kinds.** `Exception::ConversionError` from `StringUtils::toDouble` is
  `Error::Parse` (it ends the parse either way), and so is a `peptide` the
  crate's `AASequence::parse` refuses (the source's `fromString`
  `ParseError`, `m21`); `IndexOverflow` is
  `Error::InvalidRange`; the registry's `InvalidValue` stays `InvalidValue`.
  A position outside the 32-bit range, which the source's `(int)` cast would
  wrap, is refused.
- **Resource ceilings** (`ReadLimits`), each checked before its work:
  document bytes (reported by the encoding-aware decoder as `Error::Parse`),
  XML events, depth, modification definitions (every observed mass is compared
  with every variable definition), modifications per hit, residues per peptide
  (each residue setting rebuilds the sequence), sequences held, and retained
  log lines (later ones are counted).

## Checked boundaries and evidence

- **Tier 1, executed differential.** `../oracle/F3a-pepxml-mascot`: a driver
  linking the Release `libOpenMS` of
  `openms4-release-bc9cc12-c19e494-174b576`, built and run on **ibminode06**
  (Linux x86_64, conda g++ 14.4.0, libOpenMS's own core flags), loaded the
  upstream fixture and 27 generated cases, twice, byte-identical. Its stdout,
  stderr and exit codes are `tests/data/pepxml_mascot/expected/`, compared by
  `release_differential`, `both_missing_names_mass_first` and
  `a_read_after_a_failed_read_starts_empty` in `tests/pepxml_mascot.rs`: every
  title, `toString` and `toUnmodifiedString` of every sequence, every
  `Non-fatal error` line, and every exception's kind.
- **Tier 3.** The class test's only functional section,
  `load(filename, peptides)`, literal by literal (`class_test_load`); its
  constructor and destructor sections only construct.
- **Tier 4.** Every ceiling refuses a value below the fixture's need, and
  the four item ceilings accept it exactly
  (`ceilings_refuse_instead_of_truncating`); log suppression
  (`log_lines_past_the_ceiling_are_counted`).
- Provenance: `tests/data/pepxml_mascot_provenance.json`.
