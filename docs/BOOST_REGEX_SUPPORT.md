# Boost.Regex facade over fancy-regex

[`src/concept/boost_regex.rs`](../src/concept/boost_regex.rs) runs regular
expressions with the semantics of Boost.Regex, which is where OpenMS takes them
from. The engine is the `fancy-regex` crate (`=0.19.2`, features `std`,
`unicode`, `perf`), chosen in
[`THIRD_PARTY_CRATE_DECISIONS.md`](THIRD_PARTY_CRATE_DECISIONS.md): the enzyme and
RNase expressions need lookaround and the WIFF native-ID expression needs a
repeated group name, neither of which the `regex` crate has. Expressions without
such features still run on `fancy-regex`'s `regex-automata` delegate.

Tests: [`tests/boost_regex.rs`](../tests/boost_regex.rs). Manifest:
[`tests/data/boost_regex_provenance.json`](../tests/data/boost_regex_provenance.json).

This module owns no OpenMS header. It is the shared engine for the consumers that
follow in the next wave, each of which currently carries a hand-written scanner
or refuses arbitrary expressions with `Unsupported`: Mascot XML and generic
titles, pepXML, Percolator input, mzIdentML, mzTab, indexed mzML, digestion,
RNase digestion and peptide indexing. No consumer is switched in this change, so
no existing output changes.

## API mapping

| Boost.Regex | Rust |
| --- | --- |
| `boost::regex re(pattern)` | `BoostRegex::new(pattern)` |
| `boost::regex re(pattern, boost::regex::icase)` / `no_mod_s` | `BoostRegex::with_options(pattern, RegexOptions { icase, no_mod_s, .. })` |
| `regex_error` thrown by the constructor | `Error::InvalidValue` (Boost refuses the pattern too) or `Error::Unsupported` (valid Boost syntax the facade does not translate); see below |
| `basic_regex::str()` | `BoostRegex::as_str` |
| `basic_regex::mark_count()` | `BoostRegex::mark_count` |
| named sub-expressions of the pattern | `BoostRegex::capture_names` |
| `regex_search(s, m, re)` | `BoostRegex::search` -> `Option<Captures>` |
| `regex_search(s, re)` | `BoostRegex::is_search_match` |
| `regex_match(s, m, re)` | `BoostRegex::full_match` -> `Option<Captures>` |
| `regex_match(s, re)` | `BoostRegex::is_full_match` |
| `smatch::size()` | `Captures::len` |
| `smatch[0]` | `Captures::range` |
| `smatch[i]`, `sub_match::matched`, `first`/`second` | `Captures::get(i)`, `None` when the group did not take part |
| `smatch["NAME"]` (`named_subexpression`) | `Captures::name("NAME")` |
| `sregex_token_iterator(begin, end, re, subs)` | `BoostRegex::tokens(haystack, subs)` -> `Tokens`, an iterator of `Result<SubMatch>` |
| `ssub_match` (`first`, `second`, `matched`, `str()`) | `SubMatch { range, matched }`, `SubMatch::as_bytes` |
| `regex_error` (`error_complexity`, `error_stack`) thrown while matching | `Error::InvalidValue` from the backtrack limit (`RegexOptions::backtrack_limit`, default `DEFAULT_BACKTRACK_LIMIT`) |
| not in Boost | `MAX_PATTERN_BYTES`, `MAX_GROUP_DEPTH`, `MAX_REPEAT` |

Not provided, because no port target needs them yet: `regex_replace` (used by
`PercolatorOutfile`), `regex_iterator`, match flags (`match_not_bol`,
`match_prev_avail`, `match_partial`, `format_*`), searches over a sub-range with a
separate base iterator (a consumer slices the haystack instead, which is exact for
every pinned expression searched that way), the `no_mod_m`, `mod_x`, `nosubs` and
`collate` construction flags, the basic/extended/awk/grep/egrep syntaxes, `wregex`
and ICU `u32regex`. Offsets are byte positions; the haystack is `&[u8]`.

## Every Boost.Regex expression in the pinned sources

Every non-test file of core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4` that
uses Boost.Regex, and where its expressions are in the corpus:

| Source | Expressions | Call | Corpus family |
| --- | --- | --- | --- |
| `EnzymaticDigestion.cpp`, `ProteaseDigestion.cpp` | every `RegEx` of `Enzymes.xml` (29 distinct and `()`) | `sregex_token_iterator(-1)` | `ENZ` |
| `RNaseDigestion.cpp` | every comma-separated `CutsAfter`/`CutsBefore` member of `Enzymes_RNA.xml` (16 and the empty expression) | `regex_search` per nucleotide code | `RNA` |
| `SpectrumLookup.cpp`, `SpectrumMetaDataLookup.cpp` | `=(?<SCAN>\d+)$`; caller formats | `regex_search`, token iterator group 1 | `LOOKUP`, `CLASS` |
| `SpectrumNativeIDParser.cpp` (also reached from `IDMapper.cpp`, `Percolator.cpp`, `PercolatorInfile.cpp`, `USI.cpp`, `ArrowIOHelpers.cpp`) | the seven `getRegExFromNativeID` expressions and the WIFF `cycle=(?<GROUP>\d+)\s+experiment=(?<GROUP>\d+)` | token iterator `{1}` and `{1,2}` | `NATIVE` |
| `MascotXMLFile.cpp` | the three reference formats | `regex_search` with named groups | `TITLE` |
| `MzTabFile.cpp` | `^.*?\[(\d+)\].*$`, `^.*?\[\d+\].*?\[(\d+)\].*$` | token iterator group 1 | `MZTAB` |
| `IndexedMzMLDecoder.cpp` | `<[^>/]*indexListOffset\s*>\s*(\d*)` | `regex_search` | `IDX` |
| `MzIdentMLHandler.cpp` | `frag_regex_tweak` | `regex_match` | `FRAG` |
| `PepXMLFile.cpp` | `(.*?)([A-Z]+)(.*?)` on the first `)`-piece of each enzyme expression | `regex_match` | `PEPXML` |
| `FalseDiscoveryRate.cpp`, `FASTAContainer.h` | decoy prefix and suffix alternations | `regex_search` | `DECOY` |
| `SpectrumAnnotator.cpp` | the four ion-name grammars | `regex_match` | `ANNOT` |
| `MSPGenericFile.cpp` | fourteen line expressions, three with `no_mod_s` and one with `no_mod_s \| icase` | `regex_search`, repeated from the end of the previous match | `MSP` |
| `ChromeleonFile.cpp` | twelve header expressions with `no_mod_s` | `regex_match` | `CHROM` |
| `PercolatorOutfile.cpp` | three modification expressions | `regex_search` (its `regex_replace` is not provided) | `PERCOUT` |
| `MRMFeaturePickerFile.cpp`, `MRMFeatureQCFile.cpp`, `AbsoluteQuantitationMethodFile.cpp`, `XTandemInfile.cpp` | header and N-term expressions | `regex_search` | `MRM` |
| `IDMapper.cpp` | the empty expression | `regex_search` | `RNA` (the empty member) |
| `SimpleSearchEngineAlgorithm.cpp`, `ConsensusMapNormalizerAlgorithmMedian.cpp`, `MascotXMLFile::initializeLookup`, `SpectrumLookup::addReferenceFormat` | caller-supplied expressions | `regex_match`, `regex_search` | `SYN` and `FUZZ` probe the syntax |
| class tests `SpectrumLookup_test`, `SpectrumNativeIDParser_test`, `SpectrumMetaDataLookup_test`, `LogStream_test`, `LogConfigHandler_test` | their expressions | `regex_search`, `regex_match`, token iterator | `CLASS` |

`thirdparty/percolator/ScoreHolder.cpp`, `IDFilter.cpp` and `MSPFile.cpp` use
`std::regex` (ECMAScript), whose anchors and classes differ; they are not Boost
sites and this facade is not their replacement.

## Preserved Boost conventions

Each of these was measured against Boost 1.92 and is part of the corpus.

- **Bytes and ASCII classes.** Matching is by byte. `\d`, `\w`, `\s`, the POSIX
  classes and case folding follow Boost's `char` traits in the C locale: `\s` and
  `[:space:]` include the vertical tab, `[:blank:]` is space, tab and vertical tab,
  `\h` is space and tab, `\v` outside a bracket expression is `\n`, `\v`, `\f`,
  `\r`, and inside one `\v` is the vertical tab and `\b` is a backspace. The
  Arabic-Indic digit U+0663 is not `\d`. Every class is emitted as explicit byte
  ranges, so no engine class table is consulted.
- **The dot.** `.` matches every byte, line separators included, unless `(?-s)`
  or `no_mod_s` is in effect; then it excludes `\n`, `\r` and `\f` (Boost's
  `is_separator<char>`), not only `\n`.
- **Line anchors.** `^` matches at the start, after `\n` or `\f`, and after `\r`
  unless `\n` follows; `$` matches at the end, before `\r` or `\f`, and before
  `\n` unless `\r` precedes it. Neither matches inside a CRLF pair, so `.$` on
  `"\r\n"` matches the `\n`. `(?-m)` turns them into buffer anchors. `` \` `` and
  `\A` are the buffer start, `\'` and `\z` the buffer end.
- **Modifiers.** `(?imsx-imsx)` lasts to the end of the enclosing group and
  across `|`; `(?i:...)` is scoped. `(?)` is an error, `(?-)` is accepted.
- **Groups and names.** `(?<name>...)` and `(?'name'...)` capture; a name may be
  repeated. `smatch["name"]` is the first group of that name that took part; when
  none did, Boost returns its null sub-match at the end of the match. An unmatched
  group is positioned at the end of the searched range. `(?P...)` is an error.
- **Alternation and full matches.** Leftmost-first priority; `regex_match`
  backtracks into alternatives until a match spans the whole haystack.
- **The token iterator.** `-1` is the text before the match (matched only when
  non-empty), `-2` the text after it, an index past the groups Boost's null
  sub-match; after the last match the remainder follows when the first index is
  `-1`. After an empty match the next search carries `match_not_initial_null`:
  a match ending where the previous one did is refused but the search
  backtracks, so `(?=K)|K` on `AK` yields the empty match at 1 and then `K`.
  The facade reproduces this exactly with an anchored search that refuses empty
  matches, then an ordinary search one byte on.
- **Lookbehind width.** A lookbehind must have a single width under Boost's
  `calculate_backstep`: alternatives of equal width, `{n}` applied directly to a
  one-byte atom, no repeated group (`(?:a){2}` is refused), no backreference;
  nested lookarounds count zero. Variable-length lookbehind such as `(?<=K|RR)`,
  `(?<=[KR]+)` or `(?<=KR?)` is an error, as in Boost, although `fancy-regex`
  would accept some of it.
- **Quantifiers.** Nothing to repeat, a repeat of `^ $ \b \B \< \> \A \z` or of
  a repeat is an error; a `{` that does not form a bound is a literal (`a{,3}`,
  `a{x}`, `x{2}{`), spaces inside a bound are allowed (`a{ 2}`), `a{3,2}` is an
  error. A comment `(?#...)` is transparent to a following quantifier and may run
  to the end of the pattern. A repeated lookaround is tried at most once.
- **Escapes.** `\xH`, `\xHH` and `\x{H...}` up to `0x7F`; `\a \e \f \n \r \t`;
  punctuation escapes are literals. A backreference `\1`-`\9` must name an
  existing group, which may come later in the pattern. `(?=)` and `(?>)` are
  errors; `(?!)` and an empty lookbehind are accepted.
- **Bracket expressions.** `]` first is a literal, `[&&]` and `[a--b]` are
  literals, `-` is a literal first, last or as a range endpoint (`[!--]`), `[:^name:]`
  negates, an unknown POSIX name and a reversed range are errors.

## Native differences

### Byte matching through transcoding

`fancy-regex` 0.19.2 configures `regex-automata` with its default `utf8_empty`
even in `BytesMode::Ascii`. Measured consequences with a bare engine: Boost's
empty matches between the bytes of a multi-byte character are skipped (`""` and
`x*` token lists over `"é"` lose positions), and `(?s:.)*` on the one-byte
haystacks `0x85` and `0xA0` panics inside `regex-automata` ("reverse search must
match if forward search does"). The facade therefore runs the engine in Unicode
mode over a transcoded haystack: each byte above `0x7F` becomes its own
private-use code point `U+E000 + byte`, and every offset is mapped back. Patterns
are ASCII-only, so every accepted construct sees such a byte exactly as Boost
does (one character, no class, no case, non-word, distinct for backreferences).
A pure-ASCII haystack is used without copying; any other costs one pass and a
position table per call (per `tokens` iterator).

### Refused constructs

The facade refuses with `Error::Unsupported`, rather than approximate, every
construct below. Boost accepts them.

| Construct | Why |
| --- | --- |
| possessive quantifiers `*+ ++ ?+ {n}+` | not translated |
| `\Q...\E`, `\K`, `\G` | not translated |
| recursion `(?R)`, `(?1)`, `(?+1)`, `(?-1)`, `(?&name)`; conditionals `(?(...)...)`; branch reset `(?\|...)`; verbs `(*...)` | not translated |
| `\Z` | Boost's start map for `\Z` leaves out `\f`, so a leading `\Z` is never tried at a form feed (`\Z` on `"\f"` matches at 1, not 0); reproducing that bug is not worth it |
| octal `\0...`, `\c` control escapes, `\g` and `\k` backreferences, multi-digit backreferences `\10`, `\p`/`\P` properties, `\N`, `\R`, `\X`, `\C` | not translated |
| escape letters Boost reads as a literal or a locale class (`\y`, `\j`, `\l`, `\u`, `\L`, `\U`, `\E` ...), the same inside a bracket expression (`[\y]`, `[\0]`, `[a-\d]`), and `[\V]` | not translated |
| collating elements `[[.a.]]` and equivalence classes `[[=a=]]` | locale-dependent |
| `(?x)` | not translated |
| group names outside `[A-Za-z0-9_]` (Boost accepts `(?<n-x>...)` and `(?<>...)`) | not translated |
| non-ASCII pattern bytes and `\x` escapes above `0x7F` | a pattern byte above `0x7F` would have to split a UTF-8 character; on platforms where `char` is unsigned Boost also accepts `\x{80}`-`\x{FF}`, on signed-`char` platforms it rejects them |
| a capturing group inside a lookaround or an atomic group | Boost runs these as independent sub-matches and does not restore captures made inside them when the surrounding match backtracks past them: `(?!(a))b` or a failed branch after `(?>...()...)` leaves the group set. `fancy-regex` restores them. |
| a repeat with more than one repetition of a group that can match empty and captures | Boost records a final empty iteration: `(a*)+` on `"a"` gives group 1 = `1,1`, `fancy-regex` gives `0,1` |
| a repeat of a group that only asserts positions, such as `(?:$)+` | Boost never matches `(?:$)+` at all |
| a repeat of a modifier group that switches case sensitivity, such as `(?i)+a` | Boost undoes the switch when the empty repetition is abandoned, so `a` stays case-sensitive |
| nesting deeper than `MAX_GROUP_DEPTH` (48), repeat bounds above `MAX_REPEAT`, patterns longer than `MAX_PATTERN_BYTES` (64 KiB, `InvalidValue`) | work bounds |

None of these constructs occurs in a pinned OpenMS expression: every family
derived from the OpenMS sources compiles in full. How often each was hit in the
corpus is under Evidence.

### Other differences

- **Work bound.** A search may take `RegexOptions::backtrack_limit` backtracking
  steps (default `1_000_000`) before it returns `Error::InvalidValue`; the token
  iterator ends after such an error. Only expressions on the backtracking engine
  (lookaround, backreferences, atomic groups, line anchors, word boundaries)
  consume it. Boost's `error_complexity` limit is computed differently, so the
  inputs that exceed the two limits differ; below both, results are identical.
- **Error kinds and messages.** Boost's `regex_error` codes are not reproduced;
  a pattern Boost refuses is `InvalidValue` with the facade's own reason.
- **`tokens` with no sub-expression index** is `InvalidValue`; Boost indexes past
  the end of its vector.
- **Group-count safety net.** Construction fails with `Unsupported` if the
  engine program's group count ever differs from the pattern's. After the
  zero-repeat translation (`(?:X|(?!)){0}`, which keeps the groups of `X{0}`),
  no corpus pattern hits it.
- **Name hashing.** Boost looks names up by a hash and would conflate two names
  with the same hash; the facade compares names.
- **Speed.** Lookbehind-first enzyme scans run on the backtracking engine and were
  measured 4.8-5.7x slower than Boost in the survey; the generated protease table
  stays the fast path for the pinned enzymes.

## Checked boundaries and evidence

**Tier 1, executed differential.** The oracle lives outside the repository in
`../oracle/boost-regex/`: `driver.cpp` (compiled with Apple clang 21.0.0,
`-std=c++17 -O2`, against the header-only Boost.Regex 1.92 in
`/opt/homebrew/include`), `gen.py` (deterministic, seed 20260913, reads the pinned
checkout) and `build.sh`. The manifest records their sha256, the Boost headers'
sha256 and the full corpus's sha256; rerunning `gen.py` and the driver reproduced
the fixture hashes byte for byte.

The driver uses OpenMS's call shapes: construction with `perl` plus `icase` and
`no_mod_s`, `regex_search` and `regex_match` into an `smatch`, named lookups, and
`sregex_token_iterator` over index lists (`-1`, `0`, `1`, `1,2`, `-2`, `5`, `-3`,
`0,-2,-1`, `1,1,-1`). For every pattern the test requires the same compile outcome
and mark count, and for every case byte-identical output: match flag, every group
span (with Boost's position for unmatched groups), every named lookup (with the
null sub-match's position) and every token.

| | Patterns | Boost cases | Compared by the facade | Mismatches |
| --- | ---: | ---: | ---: | ---: |
| Full corpus (`../oracle/boost-regex/corpus_full.txt`) | 1,817 | 1,965,843 | 1,423,654 | 0 |
| Committed fixture (`tests/data/boost_regex_corpus.txt`) | 1,817 | 194,345 | 146,437 | 0 |

The fixture keeps every pattern, every fixed input (edge, structured and
source-derived) and every run except extra operations marked full-only, and cuts
seeded random inputs to a prefix per input set. It is 964,427 (236,307 corpus and 728,120 Boost output) bytes.
Compared cases exclude the patterns Boost refuses (which the facade must refuse
too) and those the facade refuses as listed above.

Per family in the committed fixture (the families other than `SYN` and `FUZZ` come from OpenMS sources):

| Family | Patterns | Refused by the facade | Refused by both | Compared cases |
| --- | ---: | ---: | ---: | ---: |
| `ANNOT` | 4 | 0 | 0 | 1,188 |
| `CHROM` | 12 | 0 | 0 | 1,344 |
| `CLASS` | 24 | 0 | 0 | 2,768 |
| `DECOY` | 2 | 0 | 0 | 516 |
| `ENZ` | 30 | 0 | 0 | 17,233 |
| `FRAG` | 1 | 0 | 0 | 338 |
| `FUZZ` | 1,111 | 141 | 422 | 8,768 |
| `IDX` | 1 | 0 | 0 | 231 |
| `LOOKUP` | 1 | 0 | 0 | 628 |
| `MRM` | 5 | 0 | 0 | 615 |
| `MSP` | 14 | 0 | 0 | 2,310 |
| `MZTAB` | 2 | 0 | 0 | 514 |
| `NATIVE` | 8 | 0 | 0 | 2,826 |
| `PEPXML` | 1 | 0 | 0 | 171 |
| `PERCOUT` | 3 | 0 | 0 | 354 |
| `RNA` | 17 | 0 | 0 | 6,407 |
| `SYN` | 578 | 84 | 79 | 99,308 |
| `TITLE` | 3 | 0 | 0 | 918 |

Adversarial inputs in every family: the empty string, `\n`, `\r`, `\r\n`, `\f`,
`\v`, tab, space, NUL, `\n\r`, `\r\r\n`, `\f\n`, `\r\f`, `é`, lone `0xFF`, `0xC3`,
`0x85` and `0xA0`, U+0663, NBSP, NEL and U+2028, plus structured inputs with
separators embedded (`scan=12\r\n`, `foo\f500_12`, `[1]\r\n[2]`).

Refusals over the whole corpus, by category:

| Construct | Patterns |
| --- | ---: |
| a capturing group inside a lookaround or atomic group | 111 |
| a repeat of a group that can match the empty string and captures | 34 |
| a repeat of a group that only asserts a position | 14 |
| `\Z` | 7 |
| an escape letter without a translated meaning | 7 |
| `\Q...\E` quoting | 4 |
| a named or relative backreference | 4 |
| a named-character, line-ending or grapheme escape | 4 |
| a non-ASCII byte | 4 |
| a possessive quantifier | 4 |
| a character property escape | 3 |
| a control-character escape | 3 |
| a hexadecimal escape above `0x7F` | 3 |
| a repeat of a modifier group that switches case sensitivity | 3 |
| an escape inside a character class without a translated meaning | 3 |
| `\K` | 2 |
| a backreference with more than one digit | 2 |
| a backtracking control verb | 2 |
| a conditional expression | 2 |
| an octal escape | 2 |
| `\G` | 1 |
| a branch-reset group | 1 |
| a collating element | 1 |
| an equivalence class | 1 |
| a group name outside `[A-Za-z0-9_]` | 1 |
| a recursive sub-expression | 1 |
| the `x` (extended) modifier | 1 |

84 of the 225 are syntax probes, listed one by one in `EXPECTED_UNSUPPORTED`;
141 are grammar-generated, 111 of them for a capture inside a lookaround or atomic
group, which the generator produces often. The test asserts both the list and
these counts.

**Tier 3.** `class_test_expressions` transcribes the regular-expression half of
`SpectrumNativeIDParser_test`, `SpectrumLookup_test` and
`SpectrumMetaDataLookup_test` (`spectrum=42` -> `42`, `1 2 3 42` -> last token
`42`, `rt=5.0,mz=1000.0` -> RT `5.0`, MZ `1000.0`).

**Tier 4.** Limits (`MAX_PATTERN_BYTES`, `MAX_GROUP_DEPTH`, zero backtrack limit,
empty index list), the backtracking limit as an error that also ends a token
iterator, token indices `-2` and past the groups, CRLF and form-feed anchors,
3,000 random patterns built from metacharacters, escapes and a non-ASCII
character that must compile or fail without panicking and match within bounds,
and `BoostRegex: Send + Sync`.

**Determinism.** Results are integer offsets from exact matching, with no floating
point. `memchr` and `aho-corasick` choose SIMD prefilters at run time, which
changes speed and not results; the survey measured identical `fancy-regex` output
on aarch64 macOS and x86_64 Linux. This change's tests ran on x86_64 Linux under
the current toolchain and Rust 1.85.0. No output can differ between machines.

**Crate issues found.** `fancy-regex` 0.19.2 leaves `regex-automata`'s
`utf8_empty` at its default in bytes mode, with the skipped empty matches and the
panic described above; and a repeated nullable capturing group reports the last
non-empty iteration where Perl and Boost report the final empty one. Both are
candidates for upstream reports.
