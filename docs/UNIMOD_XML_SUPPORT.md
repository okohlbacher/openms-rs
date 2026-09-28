# Unimod XML

`format::unimod_xml` (feature `kb-xml`) ports `FORMAT/UnimodXMLFile.h`,
`FORMAT/HANDLERS/UnimodXMLHandler.h` and `CHEMISTRY/UnimodXMLDataProvider.h`
from OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`: a runtime reader
for unimod.org modification files. The default build does not need it —
`ModificationsDB::global()` is embedded from the table
`tools/generate_modifications.py` flattens from the pinned
`resources/modifications/unimod.xml` offline — so this closes the gap for a
newer `unimod.xml`, OpenMS's `custom_mods.xml`, or a caller's own file:
`ModificationsDB::from_records(unimod_xml::load(path)?)`.

## API mapping

| Source member | Rust |
| --- | --- |
| `UnimodXMLFile()`, `~UnimodXMLFile()`, private copy/assign | none needed: free functions |
| `UnimodXMLFile::load(filename, std::vector<ResidueModification*>&)` | `load(path) -> Vec<ResidueModification>`, `load_with_options`, `read` / `read_with_options` over any `BufRead` |
| `UnimodXMLHandler(mods, filename)`, `~UnimodXMLHandler` | private `Handler`, one field per source member that affects the result |
| `onStartElement` / `onEndElement` / `onCharacters` | `Handler` methods driven by `sax_xml` |
| `UnimodXMLDataProvider(filename)` | `UnimodXMLDataProvider::new`, `with_options`, `filename`, `options` |
| `UnimodXMLDataProvider::loadModifications` | `load_modifications -> Vec<ResidueModification>` (owned values for `unique_ptr`) |
| per-specificity `classification` / `hidden` handling | `SpecificityAttributes::{Retain, Source}`; `ReadOptions::source()` reproduces the source |
| — (native) | `ReadLimits`, `ReadOptions` |

## Preserved source conventions

- Tags match with and without the `umod:` prefix, and only those two.
- Only `<umod:element>` children of `<umod:delta>` and `<umod:NeutralLoss>`
  build formulas; `<Ignore>`, `<aa>` and `<brick>` do not (OpenMS issue
  10030), so the difference formula weighs the recorded delta mass.
- A neutral loss without element children is dropped; its masses come from its
  own attributes.
- A site longer than one character (`N-term`, `C-term`) is the wildcard `X`;
  an unknown `position` is `Anywhere` (the source warns).
- Element symbols with leading digits become isotopes: `13C` × 6 is `(13)C6`.
- One record per `<umod:specificity>`, in document order, with provenance
  `Cv`, as `setProvenance(CV)`.
- `ReadOptions::source()`: `hidden` false and every specificity carrying the
  classification of the modification's *last* specificity, canonicalised as
  `setSourceClassification` + `getSourceClassificationName`
  (`UnimodXMLHandler.cpp:60-63`, a source `TODO`).

## Native differences

- The default `SpecificityAttributes::Retain` keeps each specificity's own
  classification and hidden flag — what the embedded registry records — so
  reading the pinned file reproduces `ModificationsDB::global()` record for
  record. `Source` is the opt-in for source-identical records.
- Records carry their full identifier (`Oxidation (M)`) from construction; the
  source leaves it to the data provider's `setFullId`. An empty `title` is
  `Error::MissingInformation`, the source's `setFullId` exception.
- `load` reads the path it is given; the source first resolves it with
  `File::find` (use `system::file::FileContext::find` for shared-data names).
- Masses must be finite. Structures the source mishandles through a dangling
  pointer (a specificity outside a modification, nested modifications) are
  refused rather than read.
- `ReadLimits` bound input bytes, depth, events, records, specificities per
  modification and losses per specificity; exceeding one is
  `Error::InvalidValue`. A DTD or non-UTF-8 encoding is `Error::Unsupported`.
- Results are returned, and only on success; the source appends raw pointers
  to the caller's vector as it goes.

## Checked boundaries and evidence

`tests/unimod_xml.rs` (7 tests, `--features kb-xml`) over
`resources/modifications/unimod.xml`, byte-identical to the pinned
`share/OpenMS/CHEMISTRY/unimod.xml` (sha256 `22231d32…91091`): the class-test
sections `load` and the two `[EXTRA]` formula/mass agreement invariants
(tier 3); the runtime reader reproducing the embedded registry, which an
independent Python path derives from the same file (tier 4 cross-check); the
source classification rule; the data provider; and limit refusals.
Provenance: `tests/data/unimod_xml_provenance.json`.
