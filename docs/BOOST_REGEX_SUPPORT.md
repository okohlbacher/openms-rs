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
no existing output changes. A consumer that accepts caller-supplied expressions
should expect the limits under [Work bounds](#work-bounds) to stop some searches
Boost answers, and should report the error rather than treat it as no match: in
particular a positive lookahead used as a filter, such as `(?=[A-Z]*K)[A-Z]+R`,
can exhaust the budget on a few kilobytes, because the engine backtracks into
lookahead bodies.

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
| not in Boost | `MAX_PATTERN_BYTES`, `MAX_GROUP_DEPTH`, `MAX_REPEAT`, `MAX_AUTOMATON_ATOMS`, `MAX_LOOKBEHIND_WIDTH`, `MAX_TRANSLATED_BYTES`, `BACKTRACK_LIMIT_BYTES`, `MAX_BACKTRACK_SCALE` |

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
| a lazy repeat with a finite maximum at least two above its minimum (`{n,m}?` with `m >= n + 2`) of a one-byte atom (a literal, `.`, a class or bracket expression) that starts the expression: nothing precedes it but group starts and ends, flag groups that keep case sensitivity, `^ $ \b \B \< \> \A \z`, comments and whole lookarounds, the expression has no backreference, no `\|` outside a lookaround precedes it or separates the alternatives of a group around it (or of the expression), and no quantifier applies to a group around it | a Boost bug, refused rather than reproduced like `\Z`. `basic_regex_creator::probe_leading_repeat` marks such a repeat as leading. When the lazy repeat extends below its maximum, Boost records the position it got to (`unwind_char_repeat`, `unwind_short_set_repeat`, `unwind_fast_dot_repeat`, `unwind_slow_dot_repeat`), and a failed start position resumes the search behind it (`match_prefix`, `find_restart_*`), skipping start positions that match: `a{1,3}?\b` on `aaaa` matches `3..4`, and `.{3,5}?\b` on `ababab` does not match, where every start position gives `1..4` and `1..6`. A flag group that changes case sensitivity adds a `toggle_case` state and an alternation or a quantifier adds a state in front, which ends Boost's walk, so `(?i)a{1,3}?\b` and `x\|a{1,3}?\b` are translated |
| octal `\0...`, `\c` control escapes, `\g` and `\k` backreferences, multi-digit backreferences `\10`, `\p`/`\P` properties, `\N`, `\R`, `\X`, `\C` | not translated |
| escape letters Boost reads as a literal or a locale class (`\y`, `\j`, `\l`, `\u`, `\L`, `\U`, `\E` ...), the same inside a bracket expression (`[\y]`, `[\0]`, `[a-\d]`), and `[\V]` | not translated |
| collating elements `[[.a.]]` and equivalence classes `[[=a=]]` | locale-dependent |
| the `x` modifier (`(?x)`) | not translated |
| group names with a character outside `[A-Za-z0-9_]` (Boost accepts `(?<n-x>...)`) | not translated |
| non-ASCII pattern bytes outside comments, and `\x` escapes above `0x7F` | a pattern byte above `0x7F` would have to split a UTF-8 character; on platforms where `char` is unsigned Boost also accepts `\x{80}`-`\x{FF}`, on signed-`char` platforms it rejects them |
| a capturing group inside a lookaround or an atomic group | Boost runs these as independent sub-matches and does not restore captures made inside them when the surrounding match backtracks past them: `(?!(a))b` or a failed branch after `(?>...()...)` leaves the group set. `fancy-regex` restores them. |
| a repeat allowing more than one repetition of a group that can match empty and captures | Boost records a final empty iteration: `(a*)+` on `"a"` gives group 1 = `1,1`, `fancy-regex` gives `0,1` |
| a repeat allowing more than one iteration (`*`, `+`, `{n,}`, and `{n}` or `{n,m}` above 1, greedy or lazy) of a group that can match the empty string, and a bounded one (`{n}`, `{n,m}` above 1) of a backreference that can | Boost ends a repeat as soon as an iteration matches the empty string, even below the minimum: it accepts the iteration and takes the exit (`repeater_count::check_null_repeat`, `perl_matcher::match_rep`). The engine differs either way. A bounded repeat iterates to its bound: `(?:\b\|a){2}b` matches `ab` in the engine and not in Boost, and `(?:a{0}\b){999999999}` ran for 16.9 s on one byte. An unbounded one (`RepeatEpsilon`) fails the empty iteration and backtracks into the group's other alternatives: `(?:b?\|a)*` on `ba` matches `0..2` in the engine and `0..1` in Boost, `(?:(?:a\|b?)*?)+b` on `abab` `0..4` and `0..2`. A backreference matches the same text in every iteration, so its unbounded repeats agree and are translated |
| a repeat of a group that only asserts positions, such as `(?:$)+` | Boost never matches `(?:$)+` at all |
| a repeat of a modifier group that switches case sensitivity, such as `(?i)+a` | Boost undoes the switch when the empty repetition is abandoned, so `a` stays case-sensitive |
| a repeat allowing more than one repetition inside an atomic group or a negative lookaround, such as `(?>a+)b`, `(?!.*x)` or `(?<!(?=.*b)a)` | work bound: the engine discards the backtracking branches pushed there without counting them |
| a lookbehind wider than `MAX_LOOKBEHIND_WIDTH` (255) bytes | work bound |
| repeat bounds whose value exceeds `MAX_REPEAT` (999,999,999), and a repeat whose shortest match, multiplied out over nested repeats, exceeds it | work bound; the product also overflows `fancy-regex`'s analyzer, which panicked in builds with overflow checks on `(?:(?:a{999999999}){999999999}){999999999}` |
| nesting deeper than `MAX_GROUP_DEPTH` (48) and a translation longer than `MAX_TRANSLATED_BYTES` (512 KiB) | work and memory bounds; a pattern longer than `MAX_PATTERN_BYTES` (64 KiB) is `InvalidValue` |

None of these constructs occurs in a pinned OpenMS expression: every family
derived from the OpenMS sources compiles in full. How often each was hit in the
corpus is under Evidence. A translation the engine itself refuses would be
`Unsupported` too; no corpus pattern is. (`(?:(?:a{999}){999}){999}`, whose
automaton the engine refused before this round, is now searched on the
backtracking machine, see [Work bounds](#work-bounds).)

### Work bounds

A search runs either on an automaton or on `fancy-regex`'s backtracking machine.
Each has its own bound, and the facade chooses per search.

**Automaton searches.** An expression that needs no backtracking (no lookaround,
backreference, atomic group, line anchor or word boundary) can be handed whole to
`fancy-regex`'s `regex-automata` delegate. That delegate spends no budget and never
fails, and its slowest mode, the PikeVM it falls back to when its lazy DFA runs
out of cache, does work proportional to the haystack length times the size of the
automaton, not of the pattern. Repeat bounds multiply that size: at commit
c93fefd, `\w{0,9999}b` (11 bytes) took 8.0 s over 100 KB, and a 64 KiB literal
took 276 s to fail over 1 MB of `a` (Boost: 33.5 s). The facade therefore counts
the *atoms* of every expression: its one-byte atoms, assertions and
backreferences with every repeat written out, a repeat multiplying the atoms of
what it repeats by its maximum, or by its minimum (at least 1) when unbounded.
`\w{0,9999}b` has 10,000, the largest expression of the pinned sources that needs
no backtracking (`TransitionGroupPicker:PeakPickerChromatogram:(.+)`) 46. A search
runs on the automaton only when the expression has at most
`MAX_AUTOMATON_ATOMS` (4,096) atoms divided by the budget's scale factor below:
4,096 up to 64 KiB, 1,024 up to 256 KiB, 256 up to 1 MiB and 64 beyond. Every
other search runs on the backtracking machine with the counted spelling, where the
budget bounds it. So an automaton search does at most about 2.7 × 10⁸ atom steps
up to 4 MiB, and 64 per byte beyond. Measured with expressions that defeat the
lazy DFA, at each tier's limit, over random `a`/`b` haystacks:

| Expression | Atoms | Haystack | Search |
| --- | ---: | --- | ---: |
| `a.{0,2046}a.{0,2046}c` | 4,095 | 64 KiB | 1.83 s, no match |
| `(?:.{0,10}a){0,372}c` | 4,093 | 64 KiB | 1.64 s, no match |
| `a.{0,510}a.{0,510}c` | 1,023 | 256 KiB | 1.96 s, no match |
| `(?:.{0,10}a){0,93}c` | 1,024 | 256 KiB | 1.80 s, no match |
| `a.{0,126}a.{0,126}c` | 255 | 1 MiB | 1.40 s, no match |
| `(?:.{0,10}a){0,23}c` | 254 | 1 MiB | 1.44 s, no match |
| `a.{0,30}a.{0,30}c` | 63 | 16 MiB | 5.72 s, no match |
| `(?:.{0,6}a){0,9}c` | 64 | 16 MiB | 6.63 s, no match |
| `a.{0,2047}a.{0,2047}c` (one atom more) | 4,097 | 64 KiB | 10 ms, budget |
| `a.{0,511}a.{0,511}c` (one atom more) | 1,025 | 256 KiB | 40 ms, budget |

An expression the lazy DFA handles is much faster at the same size (`\w{0,4094}b`
over the same 64 KiB: 0.9 ms, match). The review's slow cases, at c93fefd, now, and in Boost 1.92
(`error_complexity` is Boost's `regex_error`):

| Expression | Haystack | c93fefd | Now | Boost |
| --- | --- | --- | --- | --- |
| `(?:.{0,999})*?b`, search | 1 MB of `a`, then `b` | 11.0 s, match 0..1000001 | refused at construction | 0.57 ms, match 0..1000001 |
| `(?:.{0,999})*?b`, full match | 1 MB of `a`, then `b` | 10.7 s, match | refused at construction | 0.05 ms, match |
| `(?:.{0,9999})*?b` | 100 KB of `a`, then `b` | 11.3 s, match | refused at construction | 0.02 ms, match |
| `\w{0,9999}b` | 100 KB of `a`, then `c` | 8.0 s, no match | 56 ms, budget | 80 ms, `error_complexity` |
| `(?:[^b]{999}.{0,3})*?a{1,999}`, full match | 1 MB of `a` | 12.2 s, match 0..1000000 | 33 ms, stack | 4.6 ms, match 0..1000000 |
| `a` × 65,535 then `b` | 1 MB of `a` | 276 s, no match | 1.19 s, no match | 33.5 s, no match |

The first three are unbounded repeats of a group that can match the empty string,
refused since this round for their answers (see [Refused constructs](#refused-constructs)).

**Backtracking searches.** `fancy-regex` fails a search that takes more
backtracking steps than its budget, or that needs more than the 1,000,000 branches
of its fixed backtracking stack (`MAX_STACK`, `RuntimeError::StackOverflow`). The
facade reports both as `Error::InvalidValue` and names the limit in the message
("exceeded the backtracking limit of N steps", "exceeded the engine's backtracking
stack of 1000000 entries"); a token iterator ends after such an error. An error is
never a different answer.

**The budget is per search and shared across start positions.** The engine runs
all start positions of one search in one pass and spends a single budget over
them, so a scan that costs a few steps per position runs out on a long haystack.
The facade therefore scales the budget with the haystack: a search over `n` bytes
may take `backtrack_limit × s` steps, where `s` is the smallest of 1, 4, 16 and 64
(`MAX_BACKTRACK_SCALE`) with `n ≤ s × BACKTRACK_LIMIT_BYTES` (64 KiB), and 64
beyond. `fancy-regex` fixes a budget when it compiles a program, so the programs
for a larger budget, and for the automaton limit that goes with it, are compiled
the first time a haystack needs them and kept. Each search of a token iterator
has its own budget, scaled by the whole haystack.

**The stack is fixed.** Every greedy iteration followed by something that needs
backtracking pushes one branch, and the stack holds 1,000,000, so such a run fails
beyond about a million bytes where Boost's growable stack answers. Measured
ceilings, the longest haystack that still answers:

| Expression | Haystack | Answers up to | Then |
| --- | --- | ---: | --- |
| `=(?<SCAN>\d+)$` | `=`, digits, `x` | 999,999 digits | stack |
| `^(.+): (.+)` | `a` × n (one line) | 999,998 bytes | stack |
| `(\d+)$`, full match | digits | 1,000,000 digits | stack |
| `(?=A).*B` | `A` × n | 1,411 bytes | budget (quadratic) |
| `^(.+): (.+)` | `a\n` × n | 1,992 bytes | budget (quadratic) |
| `[A-Z]+(?=R)` | `A` × n | 1,412 bytes | budget (quadratic) |

**The engine does not count the work between backtracking steps.** Its budget
counts branches popped from the stack, and three kinds of work pushed no branch or
discarded them. Repeated from every start position, each was quadratic in the
haystack without spending the budget:

- an automaton the engine runs anchored at a position (the trailing easy part of
  an expression, the body of a lookahead) can scan to the end of the haystack;
- a counted repeat pushes no branch for the iterations below its minimum;
- the branches pushed inside an atomic group or a negative lookaround are
  discarded, uncounted, once the body matches.

The facade spells every expression it searches on the backtracking machine so that
the budget pays for this work (the *counted spelling*), and refuses a repeat inside
an atomic group or a negative lookaround:

- `fancy-regex`'s `compile_concat` hands the longest trailing run of
  sub-expressions that need no backtracking to an anchored automaton. Every
  positive lookahead body of the counted spelling ends in an always-true `(?=)`,
  so that run is empty and the body runs on the backtracking machine, where every
  loop iteration pushes a branch.
- The whole expression needs two: `(?:X)(?=)(?=)`. Unless empty matches are
  refused, `fancy-regex` first rewrites an expression ending in a positive
  lookahead, `X(?=Y)`, into an explicit group 0 around `X` followed by `Y`
  (`optimize_trailing_lookahead`). With one `(?=)` that leaves `Y` empty, and an
  `X` that needs no backtracking is then handed to an automaton whole; an `X` that
  does still ran on the machine only because the empty `Y` made group 0 a middle
  child of the concatenation, which `compile_concat` visits as hard. With two, the
  rewrite takes the second, and group 0 keeps the first and needs backtracking.
- Every iteration of a repeat with a minimum above one passes an empty alternative
  beside a never-matching one (`(?:|(?!))`), which pushes a branch.

With that, the engine does work proportional to the translated pattern and its
lookbehind widths between two backtracking steps, and a backreference compares at
most the length of its group, as Boost's does, so every search's work is bounded
by its budget times that. The translated pattern can be long: `a` × 65,535 then
`b` spends one step per start position but compares up to 64 KiB at each (1.19 s
over 1 MB, above).

**Lookahead bodies are backtracked into.** The engine compiles a positive
lookahead as a saved position around its body, without discarding the branches the
body pushed. When what follows the lookahead fails, the engine backtracks into the
body and tries its other ways of matching, which cannot change the outcome (the
facade refuses captures inside lookarounds) but costs steps. Boost runs a lookahead
as an independent sub-match and never does this. Before this spelling, a
lookahead body that needed no backtracking ran on an anchored automaton, uncounted
and atomic; now it spends the budget, so a lookahead with an ambiguous repeat
(`(?:a|aa)*`, `a*a*a*`, `.*`, `[A-Z]*`) followed by a continuation that fails can
exhaust the budget on a short haystack. These are errors, not different answers.
Measured, with the longest haystack that still answers:

| Expression | Haystack | 854d996 | Now (as c93fefd) | Boost | Now answers up to |
| --- | --- | --- | --- | --- | ---: |
| `(?=(?:a\|aa)*)b` | `a` × 60, `c` | 0.02 ms, no match | 14 ms, budget | 0.06 ms, no match | 23 `a` |
| `(?=a*a*a*)b` | `a` × 300, `c` | 0.09 ms, no match | 16 ms, budget | 0.05 ms, no match | 67 `a` |
| `(?=.*a).*b` | `a` × 4,000 | 9.8 ms, no match | 8.5 ms, budget | 0.02 ms, no match | 142 bytes |
| `(?=[A-Z]*K)[A-Z]+R` | `A` × 12,000 | 19 ms, no match | 13 ms, budget | 57 ms, no match | 1,411 bytes |

A caller-supplied filter of this shape should expect the error on inputs of a few
dozen bytes to a few kilobytes.

Measured in a release build on aarch64 macOS, the facade before the counted
spelling (854d996, as recorded in the first round), with it (measured this round),
and Boost 1.92 on the same haystack:

| Expression | Haystack | 854d996 | Now | Boost |
| --- | --- | --- | --- | --- |
| `(?:a{0}\b){999999999}` (first review) | `a` | 16.9 s, match | refused at construction | 0.007 ms, match 0..0 |
| `(?:(?=A)A){999999}` | 64 KiB of `A` | 26.5 s, no match | 30 ms, budget | 186 ms, `error_complexity` |
| `(?:(?=A)A){999999}` | 1 MB of `A` | not measured | 24 ms, match 0..999999 | 16 ms, match 0..999999 |
| `(A)\1{999999}` | 64 KiB of `A` | 12.6 s, no match | 17 ms, budget | 213 ms, `error_complexity` |
| `(?:(?=A)A){1000}B`, token iterator | 1 MB of `A` | 12.6 s, 1 token | 438 ms, budget | 534 ms, `error_complexity` |
| `(?=A).*B` | 64 KiB of `A` | 543 ms, no match | 8.8 ms, budget | 0.03 ms, no match |
| `(?=A).*B` | 1 MB of `A` | over 14 s | 4.9 ms, stack | 0.28 ms, no match |
| `^(.+): (.+)` | `a\n` × 32,768 | 302 ms, no match | 11 ms, budget | 0.02 ms, no match |
| `^(.+): (.+)` | `a\n` × 500,000 | over 14 s | 4.1 ms, stack | 0.29 ms, no match |
| `(?=.*B)A` | 64 KiB of `A` | not measured | 9.2 ms, budget | 28 ms, `error_complexity` |
| `(?>A*$)B` | 100 KB of `A` | 29.9 s, no match | refused | 0.06 ms, no match |
| `(?!A*$)x` | 100 KB of `A` | 35.4 s, no match | refused | 2,580 ms, no match |
| `(?<!(?=.*b)\w)[ab](?<=a)b`, token iterator (found by fuzzing) | `ab` × 32,768 | 12.0 s | refused | not measured |
| `(?<=a{50000})b` | 64 KiB of `a` | 7.3 s, no match | refused | 399 ms, no match |
| `(?<=a{255})b` | 64 KiB / 1 MB of `a` | not measured | 17 ms / 257 ms, budget | 10 ms / 151 ms, no match |

The same limits on expressions shaped like OpenMS's (facade now):

| Expression | Haystack | Facade | Boost |
| --- | --- | --- | --- |
| `^x` | `a\n` × 250,000 (500 KB) | 18 ms, no match (the unscaled budget ends near 333 KB) | not measured |
| `^x` | `a\n` × 8,000,000 (16 MB) | 578 ms, no match | 10 ms, no match |
| `=(?<SCAN>\d+)$` | `=`, 1,000,000 digits, `x` | stack (above) | 1 ms, no match |
| `\bx` | 16 MB of `A` | 191 ms, no match | 4 ms, no match |
| `(?<=[KRX])(?!P)`, token iterator | 16 MB of `A` | 281 ms, 1 token | 171 ms, 1 token |
| `(?:\bA\|\BA){999999}` | 1 MB of `A` | 23 ms, stack | 14 ms, match 0..999999 |

A search that is quadratic in the haystack stops with the budget where Boost's own
optimizations sometimes answer at once (`(?=A).*B`, `^(.+): (.+)`). Realistic data
is not affected: a 10 MB random protein gives 950,769 trypsin tokens in both
engines.

The work-bound fuzzer (scratch `workfuzz`) builds random expressions from
lookarounds, atomic groups, anchors, word boundaries, backreferences, repeats with
bounds up to 999,999, wide lookbehinds and easy tails, and searches every
expression that compiles with `regex_search`, `regex_match` and the first 64 tokens
over 64 KiB of `a`, of `ab` and of `a\n`, and over 1 MB of `A` followed by `B`.
In the first round, before the negative-lookbehind refusal, the slowest operation
took 12.0 s (a lookahead scanning inside a negative lookbehind); after it, the
slowest of 33,408 took 2.2 s. The second review's fuzzer found
`(?:[^b]{999}.{0,3})*?a{1,999}` (above) as its slowest. This round the fuzzer also
drew `{1,999}`, `{0,999}`, `{0,9999}`, `\w{0,2000}`, `.{0,3}`, `[^b]`, `(?i:A)`,
`(?:a|aa)` and lazy bounded repeats: with 40,000 random expressions, 3,045
compiled and 36,540 operations ran in 213 s; the slowest took 2.24 s
(`.[^b]{3,7}?a[ab]*?(?:a|aa)\w{0,2000}\w{0,2000}`, a full match over 64 KiB of `a`,
answered on the automaton at 4,013 atoms), the next four 2.1 s (the same
`\w{0,2000}\w{0,2000}` shape), and every other operation 1.5 s or less.

### Construction cost

The translation is longer than the pattern: a `^` becomes 40 bytes, a `$` 39, a
`.` 6, a literal 4 and a case-insensitive letter 9; the counted spelling adds 13
bytes per counted repeat, 4 per positive lookahead and 12 in all. The engine
compiles three programs per budget (search, full match, and the non-empty retry of
the token iterator) and holds an automaton per lookaround and per run of
sub-expressions it delegates, so memory grows with the translation; an expression
searched on the automaton also builds that automaton, whose size grows with its
atoms. `MAX_TRANSLATED_BYTES` (512 KiB) caps the translation, and
`MAX_AUTOMATON_ATOMS` the automaton: an expression with more atoms is compiled only
for the backtracking machine, whose program grows with the translation, not with
the repeat bounds. Peak resident set size of constructing one expression near the
caps, release build:

| Pattern | Pattern bytes | Construction | Peak RSS |
| --- | ---: | ---: | ---: |
| `^` × 12,190 | 12,190 | 58 ms | 67 MB |
| `^a\|` × 11,000 then `b` | 33,001 | 69 ms | 86 MB |
| `a` × 58,252, `icase` (the longest that fits) | 58,252 | 42 ms | 25 MB |
| `a` × 65,536 | 65,536 | 17 ms | 18 MB |
| `(?<=a)` × 10,922 | 65,532 | 7 ms | 18 MB |
| `^` × 12,190, then one search over 5 MB (scaled programs compiled) | 12,190 | 56 ms | 102 MB |
| `a` × 65,536, then one search over 5 MB of `a` | 65,536 | 14 ms | 24 MB |

`a` × 58,253 with `icase` translates to more than 512 KiB and is refused before
the engine sees it. At c93fefd the 65,536-atom literal was compiled as an automaton
too (37 MB) and took 10 s to find its match at the start of 5 MB of `a`; it now
runs on the backtracking machine and takes 12 ms.

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
  haystacks, facade at commit 854d996, at c93fefd and now: `(?<=[KR])(?!P)` tokens
  over a 10 MB random protein 255, 258 and 263 ms, `(?<=[KRX])` 251, 255 and
  263 ms, `(?=[DBX])` 55, 211 and 213 ms (its lookahead body runs on the
  backtracking engine since c93fefd); at 854d996 and c93fefd, `=(?<SCAN>\d+)$`
  2.4 and 2.5 ms and `^(?:Name|NAME): (.+)` with `no_mod_s` 2.6 and 2.4 ms (not
  re-measured; this round changes their spelling only by a second trailing `(?=)`).
  The generated protease table stays the fast path for the pinned enzymes.

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
| Full corpus (`../oracle/boost-regex/corpus_full.txt`) | 39,690 | 3,742,129 | 3,398,823 | 2,185,378 | 0 |
| Committed fixture (`tests/data/boost_regex_corpus.txt`) | 2,142 | 199,344 | 174,471 | 147,830 | 0 |

Boost answers every case of a pattern it compiles; compared cases exclude the
patterns the facade refuses. The fixture keeps every pattern except the two
systematic probe grids `UNB` and `LEAD` of the full corpus, every fixed input
(edge, structured and source-derived) and every run except extra operations marked
full-only, and cuts seeded random inputs to a prefix per input set. It is 1,007,508
bytes (260,222 corpus and 747,286 Boost output). A scratch build with
`MAX_AUTOMATON_ATOMS` set to 0, which searches every pattern with the counted
spelling on the backtracking machine, also compares all 2,185,378 cases of the
full corpus with 0 mismatches.

Per family in the committed fixture (the families other than `SYN`, `ADV` and
`FUZZ` come from OpenMS sources; `ADV` holds the probes of both independent
reviews and of the work bounds):

| Family | Patterns | Refused by the facade | Refused by both | Compared cases |
| --- | ---: | ---: | ---: | ---: |
| `ADV` | 325 | 102 | 9 | 7,837 |
| `ANNOT` | 4 | 0 | 0 | 1,188 |
| `CHROM` | 12 | 0 | 0 | 1,344 |
| `CLASS` | 24 | 0 | 0 | 2,768 |
| `DECOY` | 2 | 0 | 0 | 516 |
| `ENZ` | 30 | 0 | 0 | 17,269 |
| `FRAG` | 1 | 0 | 0 | 338 |
| `FUZZ` | 1,111 | 215 | 422 | 7,584 |
| `IDX` | 1 | 0 | 0 | 231 |
| `LOOKUP` | 1 | 0 | 0 | 628 |
| `MRM` | 5 | 0 | 0 | 615 |
| `MSP` | 14 | 0 | 0 | 2,310 |
| `MZTAB` | 2 | 0 | 0 | 514 |
| `NATIVE` | 8 | 0 | 0 | 2,826 |
| `PEPXML` | 1 | 0 | 0 | 171 |
| `PERCOUT` | 3 | 0 | 0 | 354 |
| `RNA` | 17 | 0 | 0 | 6,395 |
| `SYN` | 578 | 86 | 79 | 94,024 |
| `TITLE` | 3 | 0 | 0 | 918 |

The `ADV` family covers: bounded repeats of groups and backreferences that can
match the empty string, with and without their nested forms; unbounded repeats of
such groups (refused: the second review's examples, where the engine answers
differently) and of such backreferences (compared); lazy leading repeats that
Boost restarts from (refused) and their neighbours that Boost does not mark
leading (compared); expressions that need no backtracking but exceed
`MAX_AUTOMATON_ATOMS`, and lookaheads whose bodies the engine backtracks into
(compared on short inputs); repeats whose shortest match exceeds `MAX_REPEAT`; case-insensitive ranges
spanning non-letters, with and without `icase` and inline `(?i)`; bracket
expressions with several negated classes; repeat bounds as `toi` reads them (signs,
leading zeros, overflow); empty group names and non-ASCII comments; sets that match
nothing or every byte; repeats inside atomic groups and negative lookarounds and
their allowed neighbours; lookbehind widths at 255 and 256; and the counted
spelling of expressions that need backtracking.

The full corpus adds two grids. `UNB` (19,183 patterns, the second review's probe)
repeats every body of two nullable or non-nullable parts, concatenated or
alternated, with `*`, `+`, `*?`, `+?`, `{1,}` and `{0,}`, before four suffixes,
over 19 inputs with search, full match and tokens: 15,346 refused, 218,709 cases
compared. `LEAD` (18,365 patterns) puts eight atoms under ten greedy, lazy, bounded
and unbounded quantifiers before ten followers, inside twenty prefixes and wrappers
(anchors, lookarounds, alternations before and after, capturing and non-capturing
groups, case-changing and case-keeping flag groups, a quantified lookahead, a
trailing backreference), with and without `icase`: 2,506 refused, 570,924 cases
compared. Over the full corpus the facade refuses 18,255 patterns: 14,843 for a
repeat of a group that can match the empty string and 2,533 for a leading lazy
repeat among them.

Adversarial inputs in every family: the empty string, `\n`, `\r`, `\r\n`, `\f`,
`\v`, tab, space, NUL, `\n\r`, `\r\r\n`, `\f\n`, `\r\f`, `é`, lone `0xFF`, `0xC3`,
`0x85` and `0xA0`, U+0663, NBSP, NEL and U+2028, plus structured inputs with
separators embedded (`scan=12\r\n`, `foo\f500_12`, `[1]\r\n[2]`).

Refusals over the whole corpus, by category (`refusals.py`, equal to the facade's):

| Construct | Patterns |
| --- | ---: |
| a repeat inside an atomic group or a negative lookaround | 95 |
| a capturing group inside a lookaround or atomic group | 85 |
| a repeat of more than one iteration of a group that can match the empty string | 73 |
| a repeat of a group that can match the empty string and captures | 28 |
| a lazy repeat with a finite maximum of a one-byte atom that starts the expression | 27 |
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

188 of the 403 are probes outside `FUZZ`, listed one by one in
`EXPECTED_UNSUPPORTED`; 215 are grammar-generated. The test asserts the list, these
counts and the number of compared cases.

**Tier 1, transcribed.** `case_insensitive_ranges_follow_boost` and
`negated_classes_are_complemented_together` assert Boost's compile outcomes and
answers for the first review's probes, transcribed from the oracle output: `[@-Z]`
under `icase` matches `[`, `\`, `]`, `@`, `a`, `_` and `` ` ``; `[A-z]` under
`icase` does not match `[`, `_` or `` ` ``; `[a-Z]` and `[k-K]` compile under
`icase`; `[Z-a]`, `[A-_]` and `[\x41-\x5b]` under `icase` and `[a-Z]` without it
are errors; `[\S\D]` does not match `1`, space or tab; `[[:^lower:][:^upper:]]`
does not match `a` or `A`. `nullable_group_repeats_are_refused_within_a_time_bound`
and `leading_lazy_repeats_are_refused` assert the second review's probes: each of
its examples is refused with its category, and the neighbours the rules translate
give Boost's answers (`(?:b?|a)?` on `ba` matches `0..1`, `(b?)\1*a` on `bab`
`0..2` with group 1 at `0..1`; `x|a{1,3}?\b`, `(?i)a{1,3}?\b` and
`(?=x)?a{1,3}?\b` on `aaaa` match `1..4`, `a{1,2}?\b` `2..4`, `a{1,}?\b` and
`(?:a{1,3}?)+\b` `0..4`, and `a{1,3}?\b(a)\1` does not match).

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
- `nullable_group_repeats_are_refused_within_a_time_bound`: the first review's
  blocker patterns and their nested forms, and the second review's unbounded
  forms, are refused within a watchdog, and the translated neighbours answer or
  report a limit over 16 KiB within one.
- `leading_lazy_repeats_are_refused`: the second review's leading lazy repeats
  are refused within a watchdog (see tier 1 above for the answers).
- `uncounted_engine_work_is_bounded`: each family above stops with the budget on
  16 KiB, or is refused, within a watchdog.
- `large_automata_spend_the_budget`: `\w{0,9999}b` over 100 KB stops with the
  scaled budget within a watchdog, and, with a budget of one step, `\w{0,4095}b`
  (exactly `MAX_AUTOMATON_ATOMS` atoms) answers over 64 KiB on the automaton and
  stops over 64 KiB + 1, and `\w{0,4096}b` stops on 64 bytes.
- `forced_backtracking_spends_the_budget`: `(?=A).*B`, `.{0,5000}B` and
  `(?:A|B){0,4096}C` over 64 KiB of `A` stop with the budget, which they would not
  if a crate update handed their trailing run or the whole expression to an
  automaton; the module test `bounded_spelling_runs_on_the_backtracking_machine`
  checks the same on the engine directly, including that one trailing `(?=)` is
  not enough.
- `lookahead_bodies_backtrack_within_the_budget`: the four lookahead shapes above
  answer as Boost or stop with the budget, within a watchdog.
- `backtracking_budget_scales_with_the_haystack`: the budget named in the error
  at each tier boundary (64 KiB, 64 KiB + 1, 256 KiB, 256 KiB + 1, 1 MiB,
  1 MiB + 1), a 500 KB line scan that needs the scaled budget, and the stack error
  beyond a million digits.
- CRLF and form-feed anchors, the token iterator after empty matches, 3,000 random
  patterns built from metacharacters, escapes and a non-ASCII character that must
  compile or fail without panicking, and `BoostRegex: Send + Sync`.

Beyond the tests, a probe built with overflow checks and debug assertions
constructed 200,018 random patterns with large nested bounds, lookarounds and
backreferences, 9,010 of which compiled, and searched each of those over five short
inputs with `regex_search`, `regex_match` and the token iterator, without a panic.

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
  group runs all its empty iterations;
- an unbounded repeat of a group that can match the empty string fails the empty
  iteration and backtracks into the group's other alternatives, where Perl and
  Boost accept it and leave the loop (`(?:b?|a)*` on `ba`);
- positive lookaheads are not atomic, so a failing continuation backtracks into
  the lookahead body.

All are candidates for upstream reports.

**Crate update checklist.** The bounds rest on engine internals that a
`fancy-regex` or `regex-automata` update could change. Before accepting one, rerun
the full corpus and the tests, and check in the new source:

- `compile_concat` still hands only the trailing run of sub-expressions that need
  no backtracking to an automaton, and `optimize_trailing_lookahead` still removes
  at most one trailing positive lookahead at the root (the spelling
  `(?:X)(?=)(?=)` depends on both; `bounded_spelling_runs_on_the_backtracking_machine`
  and `forced_backtracking_spends_the_budget` fail otherwise);
- `RepeatGr` and `RepeatNg` still push no branch below the minimum and one per
  iteration above it, and `RepeatEpsilon` is still selected for an unbounded
  repeat whose body can match empty (`compile_repeat`);
- the backtrack limit still counts popped branches over all start positions of one
  search, and `MAX_STACK` is still 1,000,000 with `RuntimeError::StackOverflow`
  (the message the facade maps it to, and the ceilings above);
- the analyzer's shortest-match product (`MAX_REPEAT`), the lazy DFA's fallback to
  the PikeVM (the automaton ceilings above) and `utf8_empty` in bytes mode.
