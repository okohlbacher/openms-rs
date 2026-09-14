# Boost.Regex facade over fancy-regex

[`src/concept/boost_regex.rs`](../src/concept/boost_regex.rs) runs regular
expressions with the semantics of Boost.Regex, which is where OpenMS takes them
from. The engine is the `fancy-regex` crate (`=0.19.2`, features `std`,
`unicode`, `perf`), chosen in
[`THIRD_PARTY_CRATE_DECISIONS.md`](THIRD_PARTY_CRATE_DECISIONS.md): the enzyme and
RNase expressions need lookaround and the WIFF native-ID expression needs a
repeated group name, neither of which the `regex` crate has. Expressions without
such features run entirely on `fancy-regex`'s `regex-automata` delegate.

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
| `regex_error` (`error_complexity`, `error_stack`) thrown while matching | `Error::InvalidValue` from the backtracking budget (`RegexOptions::backtrack_limit`, default `DEFAULT_BACKTRACK_LIMIT`, scaled with the haystack) or from the engine's fixed stack; see [Work bounds](#work-bounds) |
| not in Boost | `MAX_PATTERN_BYTES`, `MAX_GROUP_DEPTH`, `MAX_REPEAT`, `MAX_LOOKBEHIND_WIDTH`, `MAX_TRANSLATED_BYTES`, `BACKTRACK_LIMIT_BYTES`, `MAX_BACKTRACK_SCALE` |

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
| `SimpleSearchEngineAlgorithm.cpp`, `ConsensusMapNormalizerAlgorithmMedian.cpp`, `MascotXMLFile::initializeLookup`, `SpectrumLookup::addReferenceFormat` | caller-supplied expressions | `regex_match`, `regex_search` | `SYN`, `ADV` and `FUZZ` probe the syntax |
| class tests `SpectrumLookup_test`, `SpectrumNativeIDParser_test`, `SpectrumMetaDataLookup_test`, `LogStream_test`, `LogConfigHandler_test` | their expressions | `regex_search`, `regex_match`, token iterator | `CLASS` |

`thirdparty/percolator/ScoreHolder.cpp`, `IDFilter.cpp` and `MSPFile.cpp` use
`std::regex` (ECMAScript), whose anchors and classes differ; they are not Boost
sites and this facade is not their replacement.

## Preserved Boost conventions

Each of these is exercised by the corpus and agrees with Boost 1.92 on every
compared case.

- **Bytes and ASCII classes.** Matching is by byte. `\d`, `\w`, `\s`, `\h`,
  `\v` and the POSIX classes follow Boost's `char` traits in the C locale: `\s`
  and `[:space:]` include the vertical tab, `[:blank:]` is space, tab and vertical
  tab, `\h` is space and tab, `\v` outside a bracket expression is `\n`, `\v`,
  `\f`, `\r`, and inside one `\v` is the vertical tab and `\b` is a backspace. The
  Arabic-Indic digit U+0663 is not `\d`.
- **Bracket expressions.** The translator builds every set the way Boost's
  `basic_regex_creator::append_set` fills its byte map and emits the resulting
  bytes, so the engine's class tables and case folding are not consulted for a
  set. Under `icase` Boost lower-cases both range endpoints before ordering them
  and each input byte before the lookup, so `[a-Z]` is valid, `[Z-a]` and `[A-_]`
  are errors, and `[@-Z]` matches `[`, `\` and `]`. All negated classes of one
  expression are complemented once, together, so `[\S\D]` matches only bytes that
  are neither space nor digit. `]` first is a literal, `[&&]` and `[a--b]` are
  literals, `-` is a literal first, last or as a range endpoint (`[!--]`),
  `[:^name:]` negates, and an unknown POSIX name or a reversed range is an error.
- **Case-insensitive literals, backreferences and word boundaries.** These use
  the engine's Unicode tables, which agree with Boost's C-locale tables on every
  character a haystack can contain here: ASCII and the transcoded bytes described
  below.
- **The dot.** `.` matches every byte, line separators included, unless `(?-s)`
  or `no_mod_s` is in effect; then it excludes `\n`, `\r` and `\f` (Boost's
  `is_separator<char>`), not only `\n`.
- **Line anchors.** `^` matches at the start, after `\n` or `\f`, and after `\r`
  unless `\n` follows; `$` matches at the end, before `\r` or `\f`, and before
  `\n` unless `\r` precedes it. Neither matches inside a CRLF pair, so `.$` on
  `"\r\n"` matches the `\n`. `(?-m)` turns them into buffer anchors. `` \` `` and
  `\A` are the buffer start, `\'` and `\z` the buffer end.
- **Modifiers.** `(?ims-ims)` lasts to the end of the enclosing group and across
  `|`; `(?i:...)` is scoped. `(?)` is an error, `(?-)` is accepted. `x` is refused
  (below).
- **Groups and names.** `(?<name>...)` and `(?'name'...)` capture; a name may be
  repeated, and may be empty (`(?<>a)`, `(?''a)`), as in Boost. `smatch["name"]`
  is the first group of that name that took part; when none did, Boost returns its
  null sub-match at the end of the match. An unmatched group is positioned at the
  end of the searched range. `(?P...)` is an error.
- **Comments.** `(?#...)` ends at the first `)` or the pattern end, appends
  nothing, and is transparent to a following quantifier. Its content is skipped,
  non-ASCII bytes included (`(?#é)a` is accepted).
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
  a repeat is an error. A bound is read as Boost's `cpp_regex_traits::toi` reads
  it (`std::istream >> intmax_t`): an optional sign, decimal digits, leading zeros
  allowed (`a{0000000001}` is `a{1}`), spaces around the numbers (`a{ 2}`); a
  negative maximum is unbounded (`a{2,-1}` is `a{2,}`); a `{` that does not form a
  bound is a literal (`a{,3}`, `a{x}`, `x{2}{`, a value outside `intmax_t`);
  `a{3,2}` is an error. A repeated lookaround is tried at most once.
- **Escapes.** `\xH`, `\xHH` and `\x{H...}` up to `0x7F`; `\a \e \f \n \r \t`;
  punctuation escapes are literals. A backreference `\1`-`\9` must name an
  existing group, which may come later in the pattern. `(?=)` and `(?>)` are
  errors; `(?!)` and an empty lookbehind are accepted.

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
are ASCII-only outside comments, so every accepted construct sees such a byte
exactly as Boost does (one character, no class, no case, non-word, distinct for
backreferences). A pure-ASCII haystack is used without copying. Any other costs
one pass per call (per `tokens` iterator) and memory of about eleven bytes per
haystack byte, three for the transcoded character and eight for the offset table:
one search over 16 MB of `0xFF` peaks at 195 MB resident, against 19 MB for 16 MB
of ASCII (measured, release build, aarch64 macOS).

### Refused constructs

The facade refuses with `Error::Unsupported`, rather than approximate, every
construct below. Boost accepts them. Refusals that exist to keep work or memory
bounded are explained under [Work bounds](#work-bounds) and
[Construction cost](#construction-cost).

| Construct | Why |
| --- | --- |
| possessive quantifiers `*+ ++ ?+ {n}+` | not translated |
| `\Q...\E`, `\K`, `\G` | not translated |
| recursion `(?R)`, `(?1)`, `(?+1)`, `(?-1)`, `(?&name)`; conditionals `(?(...)...)`; branch reset `(?\|...)`; verbs `(*...)` | not translated |
| `\Z` | Boost's start map for `\Z` leaves out `\f`, so a leading `\Z` is never tried at a form feed (`\Z` on `"\f"` matches at 1, not 0); reproducing that bug is not worth it |
| octal `\0...`, `\c` control escapes, `\g` and `\k` backreferences, multi-digit backreferences `\10`, `\p`/`\P` properties, `\N`, `\R`, `\X`, `\C` | not translated |
| escape letters Boost reads as a literal or a locale class (`\y`, `\j`, `\l`, `\u`, `\L`, `\U`, `\E` ...), the same inside a bracket expression (`[\y]`, `[\0]`, `[a-\d]`), and `[\V]` | not translated |
| collating elements `[[.a.]]` and equivalence classes `[[=a=]]` | locale-dependent |
| the `x` modifier (`(?x)`) | not translated |
| group names with a character outside `[A-Za-z0-9_]` (Boost accepts `(?<n-x>...)`) | not translated |
| non-ASCII pattern bytes outside comments, and `\x` escapes above `0x7F` | a pattern byte above `0x7F` would have to split a UTF-8 character; on platforms where `char` is unsigned Boost also accepts `\x{80}`-`\x{FF}`, on signed-`char` platforms it rejects them |
| a capturing group inside a lookaround or an atomic group | Boost runs these as independent sub-matches and does not restore captures made inside them when the surrounding match backtracks past them: `(?!(a))b` or a failed branch after `(?>...()...)` leaves the group set. `fancy-regex` restores them. |
| a repeat allowing more than one repetition of a group that can match empty and captures | Boost records a final empty iteration: `(a*)+` on `"a"` gives group 1 = `1,1`, `fancy-regex` gives `0,1` |
| a counted repeat (`{n}`, `{n,}`, `{n,m}` with `n` or `m` above 1) of a group or backreference that can match the empty string | Boost ends a repeat as soon as an iteration matches the empty string, even below the minimum (`repeater_count::check_null_repeat`); the engine iterates to the bound. `(?:\b\|a){2}b` matches `ab` in the engine and not in Boost, and `(?:a{0}\b){999999999}` ran for 16.9 s on one byte |
| a repeat of a group that only asserts positions, such as `(?:$)+` | Boost never matches `(?:$)+` at all |
| a repeat of a modifier group that switches case sensitivity, such as `(?i)+a` | Boost undoes the switch when the empty repetition is abandoned, so `a` stays case-sensitive |
| a repeat allowing more than one repetition inside an atomic group or a negative lookaround, such as `(?>a+)b`, `(?!.*x)` or `(?<!(?=.*b)a)` | work bound: the engine discards the backtracking branches pushed there without counting them |
| a lookbehind wider than `MAX_LOOKBEHIND_WIDTH` (255) bytes | work bound |
| repeat bounds whose value exceeds `MAX_REPEAT` (999,999,999), and a repeat whose shortest match, multiplied out over nested repeats, exceeds it | work bound; the product also overflows `fancy-regex`'s analyzer, which panicked in builds with overflow checks on `(?:(?:a{999999999}){999999999}){999999999}` |
| nesting deeper than `MAX_GROUP_DEPTH` (48) and a translation longer than `MAX_TRANSLATED_BYTES` (512 KiB) | work and memory bounds; a pattern longer than `MAX_PATTERN_BYTES` (64 KiB) is `InvalidValue` |

None of these constructs occurs in a pinned OpenMS expression: every family
derived from the OpenMS sources compiles in full. How often each was hit in the
corpus is under Evidence. A translation the engine itself refuses (for example an
automaton over its size limit, such as `(?:(?:a{999}){999}){999}`) is
`Unsupported` too.

### Work bounds

`fancy-regex` fails a search that takes more backtracking steps than its budget,
or that needs more than the 1,000,000 branches of its fixed backtracking stack
(`MAX_STACK`, `RuntimeError::StackOverflow`). The facade reports both as
`Error::InvalidValue` and names the limit in the message ("exceeded the
backtracking limit of N steps", "exceeded the engine's backtracking stack of
1000000 entries"); a token iterator ends after such an error. An error is never a
different answer.

**The budget is per search and shared.** One search spends a single budget over
every start position it tries, so a scan that costs a few steps per position runs
out on a long haystack. The facade therefore scales the budget with the haystack:
a search over `n` bytes may take `backtrack_limit × s` steps, where `s` is the
smallest of 1, 4, 16 and 64 (`MAX_BACKTRACK_SCALE`) with
`n ≤ s × BACKTRACK_LIMIT_BYTES` (64 KiB), and 64 beyond. `fancy-regex` fixes a
budget when it compiles a program, so the programs for a larger budget are
compiled the first time a haystack needs them and kept. Each search of a token
iterator has its own budget, scaled by the whole haystack. Only expressions that
need the backtracking engine (lookaround, backreferences, atomic groups, line
anchors, word boundaries) spend budget; the others run entirely on an automaton,
in time proportional to the haystack length times the pattern length.

**The engine does not count the work between backtracking steps.** Its budget
counts branches popped from the stack, and three kinds of work pushed no branch or
discarded them. Repeated from every start position, each was quadratic in the
haystack without spending the budget:

- an automaton the engine runs anchored at a position (the trailing easy part of
  an expression, the body of a lookahead) can scan to the end of the haystack;
- a counted repeat pushes no branch for the iterations below its minimum;
- the branches pushed inside an atomic group or a negative lookaround are
  discarded, uncounted, once the body matches.

The facade spells every expression that needs backtracking so that the budget
pays for this work: the expression and every positive lookahead in it end in an
always-true `(?=)`, which makes the engine run the sub-expressions before it on
its backtracking machine instead of an automaton; every iteration of a repeat with
a minimum above one passes an empty alternative beside a never-matching one
(`(?:|(?!))`), which pushes a branch; and a repeat inside an atomic group or a
negative lookaround is refused. An expression that needs no backtracking keeps its
plain spelling and runs entirely on an automaton. With that, the engine does work
proportional to the translated pattern and its lookbehind widths between two
backtracking steps, and a backreference compares at most the length of its group,
as Boost's does, so every search's work is bounded by its budget times that.

Measured in a release build on aarch64 macOS, the facade before this spelling,
after it, and Boost 1.92 on the same haystack (`error_complexity` is Boost's
`regex_error`):

| Expression | Haystack | Before | After | Boost |
| --- | --- | --- | --- | --- |
| `(?:a{0}\b){999999999}` (review) | `a` | 16.9 s, match | refused at construction | 0.007 ms, match 0..0 |
| `(?:(?=A)A){999999}` | 64 KiB of `A` | 26.5 s, no match | 48 ms, budget | 186 ms, `error_complexity` |
| `(?:(?=A)A){999999}` | 1 MB of `A` | not measured | 39 ms, match 0..999999 | 16 ms, match 0..999999 |
| `(A)\1{999999}` | 64 KiB of `A` | 12.6 s, no match | 35 ms, budget | 213 ms, `error_complexity` |
| `(?:(?=A)A){1000}B`, token iterator | 1 MB of `A` | 12.6 s, 1 token | 781 ms, budget | 534 ms, `error_complexity` |
| `(?=A).*B` | 64 KiB of `A` | 543 ms, no match | 16 ms, budget | 0.03 ms, no match |
| `(?=A).*B` | 1 MB of `A` | over 14 s | 13 ms, stack | 0.28 ms, no match |
| `^(.+): (.+)` | `a\n` × 32,768 | 302 ms, no match | 20 ms, budget | 0.02 ms, no match |
| `^(.+): (.+)` | `a\n` × 500,000 | over 14 s | 10 ms, stack | 0.29 ms, no match |
| `(?=.*B)A` | 64 KiB of `A` | not measured | 14 ms, budget | 28 ms, `error_complexity` |
| `(?>A*$)B` | 100 KB of `A` | 29.9 s, no match | refused | 0.06 ms, no match |
| `(?!A*$)x` | 100 KB of `A` | 35.4 s, no match | refused | 2,580 ms, no match |
| `(?<!(?=.*b)\w)[ab](?<=a)b`, token iterator (found by fuzzing) | `ab` × 32,768 | 12.0 s | refused | not measured |
| `(?<=a{50000})b` | 64 KiB of `a` | 7.3 s, no match | refused | 399 ms, no match |
| `(?<=a{255})b` | 64 KiB / 1 MB of `a` | not measured | 27 ms / 424 ms, budget | 10 ms / 151 ms, no match |

The same limits on expressions shaped like OpenMS's:

| Expression | Haystack | Facade | Boost |
| --- | --- | --- | --- |
| `^x` | `a\n` × 250,000 (500 KB) | 29 ms, no match (the unscaled budget ends near 333 KB) | not measured |
| `^x` | `a\n` × 8,000,000 (16 MB) | 1,061 ms, no match | 10 ms, no match |
| `=(?<SCAN>\d+)$` | `=`, 100,000 digits, `x` | 10 ms, no match | 0.1 ms, no match |
| `=(?<SCAN>\d+)$` | `=`, 1,000,000 digits, `x` | 13 ms, stack | 1 ms, no match |
| `\bx` | 16 MB of `A` | 329 ms, no match | 4 ms, no match |
| `(?<=[KRX])(?!P)`, token iterator | 16 MB of `A` | 535 ms, 1 token | 171 ms, 1 token |
| `(?:\bA\|\BA){999999}` | 1 MB of `A` | 41 ms, stack | 14 ms, match 0..999999 |

A greedy repeat followed by something that needs backtracking pushes one branch per
byte it takes, so it fails with the stack error beyond about a million bytes where
Boost's growable stack answers, and a search that is quadratic in the haystack
stops with the budget where Boost's own optimizations sometimes answer at once
(`(?=A).*B`, `^(.+): (.+)`). Realistic data is not affected: a 10 MB random
protein gives 950,769 trypsin tokens in both engines.

The work-bound fuzzer (scratch `workfuzz`, 30,000 random expressions built from
lookarounds, atomic groups, anchors, word boundaries, backreferences, repeats with
bounds up to 999,999, wide lookbehinds and easy tails) searched every expression
that compiled with `regex_search`, `regex_match` and the first 64 tokens over
64 KiB of `a`, of `ab` and of `a\n`, and over 1 MB of `A` followed by `B`. Before
the negative-lookbehind refusal, 2,187 expressions compiled and 26,244 operations
ran in 265 s; the slowest took 12.0 s, all three of that kind a lookahead scanning
inside a negative lookbehind, and the next 1.3 s. After it, with a fresh seed and
40,000 random expressions, 2,784 compiled and 33,408 operations ran in 191 s: the
slowest took 2.2 s (the first 64 tokens of `\w*?^|` over 1 MB, answered), the next
1.2 s, and every other operation less than 0.7 s.

### Construction cost

The translation is longer than the pattern: a `^` becomes 40 bytes, a `$` 39, a
`.` 6, a literal 4 and a case-insensitive letter 9; the counted spelling adds 13
bytes per counted repeat, 4 per positive lookahead and 7 in all. The engine
compiles three programs per budget (search, full match, and the non-empty retry of
the token iterator) and holds an automaton per lookaround, so memory grows with
the translation. `MAX_TRANSLATED_BYTES` (512 KiB) caps it; peak resident set size
of constructing one expression near the cap, release build:

| Pattern | Pattern bytes | Construction | Peak RSS |
| --- | ---: | ---: | ---: |
| `^` × 12,190 | 12,190 | 58 ms | 66 MB |
| `^a\|` × 11,000 then `b` | 33,001 | 72 ms | 86 MB |
| `a` × 58,000, `icase` | 58,000 | 118 ms | 52 MB |
| `a` × 65,536 | 65,536 | 27 ms | 37 MB |
| `(?<=a)` × 10,922 | 65,532 | 7 ms | 19 MB |
| `^` × 12,190, then one search over 5 MB (scaled programs compiled) | 12,190 | 58 ms | 102 MB |

`a` × 65,536 with `icase` translates to more than 512 KiB and is refused before
the engine sees it. An expression that needs no backtracking is searched by an
automaton whose cost grows with the pattern length too: the 64 KiB literal above
took 10 s to find its match at the start of 5 MB of `a`.

### Other differences

- **Error kinds and messages.** Boost's `regex_error` codes are not reproduced;
  a pattern Boost refuses is `InvalidValue` with the facade's own reason.
- **`tokens` with no sub-expression index** is `InvalidValue`; Boost indexes past
  the end of its vector.
- **Group-count safety net.** Construction fails with `Unsupported` if any
  engine program (search, full match, the non-empty retry, for every budget
  compiled) numbers its groups differently from the pattern. After the
  zero-repeat translation (`(?:X|(?!)){0}`, which keeps the groups of `X{0}`), no
  corpus pattern hits it.
- **Name hashing.** Boost looks names up by a hash and would conflate two names
  with the same hash; the facade compares names.
- **Speed.** Expressions on the backtracking engine are slower than Boost (above).
  The counted spelling costs little on OpenMS's expressions; release build, same
  haystacks, facade at commit 854d996 and now: `(?<=[KR])(?!P)` tokens over a
  10 MB protein 255 and 258 ms, `(?<=[KRX])` 251 and 255 ms, `(?=[DBX])` 55 and
  211 ms (its lookahead body now runs on the backtracking engine),
  `=(?<SCAN>\d+)$` 2.4 and 2.5 ms, `^(?:Name|NAME): (.+)` with `no_mod_s` 2.6 and
  2.4 ms. The generated protease table stays the fast path for the pinned enzymes.

## Checked boundaries and evidence

**Tier 1, executed differential.** The oracle lives outside the repository in
`../oracle/boost-regex/`: `driver.cpp` (compiled with Apple clang 21.0.0,
`-std=c++17 -O2`, against the header-only Boost.Regex 1.92 in
`/opt/homebrew/include`), `gen.py` (deterministic, seed 20260913, reads the pinned
checkout), `build.sh <checkout or worktree>` (builds the driver, regenerates both
corpora, runs Boost and copies the fixture pair into `tests/data`) and
`refusals.py`. The manifest records their sha256, the Boost headers' sha256 and the
full corpus's sha256. Rerunning `gen.py` in a scratch directory, with a driver
built there from `driver.cpp`, reproduced both corpora and both Boost outputs byte
for byte.

The driver uses OpenMS's call shapes: construction with `perl` plus `icase` and
`no_mod_s`, `regex_search` and `regex_match` into an `smatch`, named lookups, and
`sregex_token_iterator` over index lists (`-1`, `0`, `1`, `1,2`, `-2`, `5`, `-3`,
`0,-2,-1`, `1,1,-1`). For every pattern the test requires the same compile outcome
and mark count, and for every case byte-identical output: match flag, every group
span (with Boost's position for unmatched groups), every named lookup (with the
null sub-match's position) and every token.

Which patterns the facade refuses, by category, and therefore how many cases it
compares, is derived on the oracle side. `refusals.py` re-implements the refusal
rules of this document in Python, independently of the Rust translator, applies
them to every corpus pattern and combines them with Boost's compile outcomes; the
test's `EXPECTED_UNSUPPORTED`, `EXPECTED_REFUSALS` and `EXPECTED_CASES` are its
output. On both corpora it agrees with the facade on every refusal and category,
and it finds no pattern Boost rejects that the rules would accept.

| | Patterns | Corpus cases | Boost answers | Compared by the facade | Mismatches |
| --- | ---: | ---: | ---: | ---: | ---: |
| Full corpus (`../oracle/boost-regex/corpus_full.txt`) | 2,068 | 1,983,043 | 1,639,737 | 1,402,806 | 0 |
| Committed fixture (`tests/data/boost_regex_corpus.txt`) | 2,068 | 197,058 | 172,185 | 147,209 | 0 |

Boost answers every case of a pattern it compiles; compared cases exclude the
patterns the facade refuses. The fixture keeps every pattern, every fixed input
(edge, structured and source-derived) and every run except extra operations marked
full-only, and cuts seeded random inputs to a prefix per input set. It is 990,331
bytes (253,972 corpus and 736,359 Boost output).

Per family in the committed fixture (the families other than `SYN`, `ADV` and
`FUZZ` come from OpenMS sources; `ADV` holds the probes of the independent review
and of this round's work bounds):

| Family | Patterns | Refused by the facade | Refused by both | Compared cases |
| --- | ---: | ---: | ---: | ---: |
| `ADV` | 251 | 61 | 9 | 6,767 |
| `ANNOT` | 4 | 0 | 0 | 1,188 |
| `CHROM` | 12 | 0 | 0 | 1,344 |
| `CLASS` | 24 | 0 | 0 | 2,768 |
| `DECOY` | 2 | 0 | 0 | 516 |
| `ENZ` | 30 | 0 | 0 | 17,269 |
| `FRAG` | 1 | 0 | 0 | 338 |
| `FUZZ` | 1,111 | 201 | 422 | 7,808 |
| `IDX` | 1 | 0 | 0 | 231 |
| `LOOKUP` | 1 | 0 | 0 | 628 |
| `MRM` | 5 | 0 | 0 | 615 |
| `MSP` | 14 | 0 | 0 | 2,310 |
| `MZTAB` | 2 | 0 | 0 | 514 |
| `NATIVE` | 8 | 0 | 0 | 2,826 |
| `PEPXML` | 1 | 0 | 0 | 171 |
| `PERCOUT` | 3 | 0 | 0 | 354 |
| `RNA` | 17 | 0 | 0 | 6,395 |
| `SYN` | 578 | 85 | 79 | 94,249 |
| `TITLE` | 3 | 0 | 0 | 918 |

The `ADV` family covers: counted repeats of nullable groups and backreferences,
with and without their nested forms, and the unbounded forms the engine agrees on;
repeats whose shortest match exceeds `MAX_REPEAT`; case-insensitive ranges
spanning non-letters, with and without `icase` and inline `(?i)`; bracket
expressions with several negated classes; repeat bounds as `toi` reads them (signs,
leading zeros, overflow); empty group names and non-ASCII comments; sets that match
nothing or every byte; repeats inside atomic groups and negative lookarounds and
their allowed neighbours; lookbehind widths at 255 and 256; and the counted
spelling of expressions that need backtracking.

Adversarial inputs in every family: the empty string, `\n`, `\r`, `\r\n`, `\f`,
`\v`, tab, space, NUL, `\n\r`, `\r\r\n`, `\f\n`, `\r\f`, `é`, lone `0xFF`, `0xC3`,
`0x85` and `0xA0`, U+0663, NBSP, NEL and U+2028, plus structured inputs with
separators embedded (`scan=12\r\n`, `foo\f500_12`, `[1]\r\n[2]`).

Refusals over the whole corpus, by category (`refusals.py`, equal to the facade's):

| Construct | Patterns |
| --- | ---: |
| a repeat inside an atomic group or a negative lookaround | 95 |
| a capturing group inside a lookaround or atomic group | 85 |
| a counted repeat of a group that can match the empty string | 43 |
| a repeat of a group that can match the empty string and captures | 29 |
| a repeat of a group that only asserts a position | 11 |
| a counted repeat of a backreference that can match the empty string | 10 |
| `\Z` | 7 |
| an escape letter without a translated meaning | 7 |
| a repeat whose shortest match is longer than `MAX_REPEAT` bytes | 5 |
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
| a lookbehind wider than `MAX_LOOKBEHIND_WIDTH` bytes | 2 |
| an octal escape | 2 |
| `\G` | 1 |
| a branch-reset group | 1 |
| a collating element | 1 |
| an equivalence class | 1 |
| a group name outside `[A-Za-z0-9_]` | 1 |
| a recursive sub-expression | 1 |
| a repeat bound above `MAX_REPEAT` | 1 |
| the `x` (extended) modifier | 1 |

146 of the 347 are probes outside `FUZZ`, listed one by one in
`EXPECTED_UNSUPPORTED`; 201 are grammar-generated. The test asserts the list, these
counts and the number of compared cases.

**Tier 1, transcribed.** `case_insensitive_ranges_follow_boost` and
`negated_classes_are_complemented_together` assert Boost's compile outcomes and
answers for the review's probes, transcribed from the oracle output: `[@-Z]` under
`icase` matches `[`, `\`, `]`, `@`, `a`, `_` and `` ` ``; `[A-z]` under `icase`
does not match `[`, `_` or `` ` ``; `[a-Z]` and `[k-K]` compile under `icase`;
`[Z-a]`, `[A-_]` and `[\x41-\x5b]` under `icase` and `[a-Z]` without it are
errors; `[\S\D]` does not match `1`, space or tab; `[[:^lower:][:^upper:]]` does
not match `a` or `A`.

**Tier 3.** `class_test_expressions` transcribes the regular-expression half of
`SpectrumNativeIDParser_test`, `SpectrumLookup_test` and
`SpectrumMetaDataLookup_test` (`spectrum=42` -> `42`, `1 2 3 42` -> last token
`42`, `rt=5.0,mz=1000.0` -> RT `5.0`, MZ `1000.0`).

**Tier 4.**

- `limits_are_errors`: `MAX_PATTERN_BYTES`, `MAX_GROUP_DEPTH`, a zero budget, an
  empty index list, token indices `-2` and past the groups, a lookbehind of exactly
  `MAX_LOOKBEHIND_WIDTH` bytes (matches, as Boost's `(?<=A{255})K` does at 255..256)
  and one byte wider (refused), and the `MAX_TRANSLATED_BYTES` boundary for
  case-insensitive letters.
- `nested_repeat_bounds_do_not_overflow_the_engine`: the 42-byte
  `(?:(?:a{999999999}){999999999}){999999999}` (the review's "43-byte" pattern)
  and its relatives are refused in
  the test profile, which has overflow checks on and in which they panicked inside
  the engine.
- `nullable_counted_repeats_are_refused_within_a_time_bound`: the review's blocker
  patterns and their nested forms are refused within a watchdog, and the unbounded
  forms answer or report a limit over 16 KiB within one.
- `uncounted_engine_work_is_bounded`: each family above stops with the budget on
  16 KiB, or is refused, within a watchdog.
- `backtracking_budget_scales_with_the_haystack`: the budget named in the error
  at each tier boundary (64 KiB, 64 KiB + 1, 256 KiB, 256 KiB + 1, 1 MiB,
  1 MiB + 1), a 500 KB line scan that needs the scaled budget, and the stack error
  beyond a million digits.
- CRLF and form-feed anchors, the token iterator after empty matches, 3,000 random
  patterns built from metacharacters, escapes and a non-ASCII character that must
  compile or fail without panicking, and `BoostRegex: Send + Sync`.

Beyond the tests, a debug-profile probe (overflow checks on) constructed 200,018
random patterns with large nested bounds, lookarounds and backreferences, 8,336 of
which compiled, without a panic.

**Determinism.** Results are integer offsets from exact matching, with no floating
point. The budgets are fixed per haystack length and the engine's limits are
constants, so which searches fail does not depend on the machine. `memchr` and
`aho-corasick` choose SIMD prefilters at run time, which changes speed and not
results; the survey measured identical `fancy-regex` output on aarch64 macOS and
x86_64 Linux, and this change's tests ran on x86_64 Linux under the current
toolchain and Rust 1.85.0.

**Crate issues found.** In `fancy-regex` 0.19.2:

- `regex-automata`'s `utf8_empty` stays at its default in bytes mode, with the
  skipped empty matches and the panic described above;
- a repeated nullable capturing group reports the last non-empty iteration where
  Perl and Boost report the final empty one;
- the analyzer multiplies a repeat's shortest match by its minimum without an
  overflow check, which panics in builds with overflow checks;
- the backtrack limit counts popped branches only, so a counted repeat below its
  minimum, an automaton delegate scanning ahead, and the branches an atomic group
  or negative lookaround discards all escape it, and a finite repeat of a nullable
  group runs all its empty iterations.

All are candidates for upstream reports.
