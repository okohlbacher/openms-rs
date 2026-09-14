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
lookahead bodies; and a case-insensitive backreference to a group that grows
across the haystack can take seconds on a megabyte, because the engine compares
the whole group on every attempt. Such a consumer should also expect `Unsupported` for the
constructs under [Refused constructs](#refused-constructs), for example a
backreference inside the group it refers to, `\<` in an expression that switches case
sensitivity, or a case switch inside a repeated group, where Boost's own start maps
lose matches.

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
  `[:^name:]` negates, and a reversed range is an error. A class name is looked up
  as `cpp_regex_traits::lookup_classname` does: the lookup is retried with the name
  lower-cased, and `get_default_class_id`'s table carries the one-letter aliases
  `d h l s u v w` beside the POSIX names, so `[[:ALPHA:]]` and `[[:Alpha:]]` are
  `[[:alpha:]]`, `[[:D:]]` and `[[:d:]]` are `[[:digit:]]`, `[[:h:]]` is `\h` and
  `[[:v:]]` is `\n`-`\r`; a name neither lookup finds is an error, as Boost's
  `error_ctype`. A `[` inside a bracket expression is a literal unless a `.`, `=`
  or `:` follows; at either endpoint of a range only `[.` has a meaning
  (`get_next_set_literal`), so `[A-[x]` is the range `A` to `[` plus `x` and
  `[A-[.a.]]` is the range `A` to `a`. Collating elements, equivalence classes and
  Boost's own `[[:unicode:]]` are refused (below).
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
  end of the searched range. `(?P` is not the Python prefix it looks like:
  `parse_perl_extension` consumes the `P`, reads `(?P>name)` as a recursion by
  name (refused, below) and hands everything else to the option-group parser, so
  `(?Pi)` is `(?i)`, `(?P:a|b)` is `(?:a|b)`, `(?P)` is an empty option group and
  `(?Pim-s:a)` is `(?im-s:a)`. Only the Python spellings `(?P<name>...)` and
  `(?P=name)` are errors, because the option parser rejects their `<` and `=`.
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
- **Escapes.** `\xH`, `\xHH` and `\x{H...}` up to `0x7F`, with hex digits only
  (Boost also reads white space, a sign and `0x` there, which is refused, below);
  `\a \e \f \n \r \t`;
  punctuation escapes are literals. A backreference `\1`-`\9` must name an
  existing group, which may come later in the pattern; one inside the group it
  names is refused (below). A backreference to a group that did not take part
  fails, and one under `icase` compares ASCII letters without case and every
  byte above `0x7F` exactly. `(?=)` and `(?>)` are errors; `(?!)` and an empty
  lookbehind are accepted.

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

### Engine rewrites

Before compiling an expression, the engine crates rewrite some shapes into forms
that match the same strings but choose among them differently. Boost, like Perl,
takes the first match in the order its backtracking tries them (leftmost-first),
and reports the captures of that match; the rewritten forms give other answers.
The third round's fuzzing found three rewrites in `fancy-regex` 0.19.2's optimizer
and one in `regex-syntax` 0.8.11, and the fourth review found that the last one also
applies inside atomic groups:

| Rewrite | Shape | Example: engine before the fix / Boost |
| --- | --- | --- |
| `optimize_ambiguous_concat_repeats` | `X+Y?X+` and `X+Y*X+` become `X+(?:Y{1}X+)?`, which also matches a single `X`; related shapes become `(?:X*Y{1})?X+` and `X+(?:YX*)*` | `\w{1,}b?\w{1,}` on `a`: match `0..1` / no match; `(?:a+(?:ba*)?)+$` on `abba`: `0..4` / `3..4` |
| `optimize_nested_repeats`, repeat of a group holding one repeat | `(X+)+` becomes `(X+)`, also when a backreference reads the group | `(.{1,})+\1+` on `abb`: `1..3` / `0..3` with group 1 at `1..2` |
| `optimize_nested_repeats`, without backreferences | `(X)*` whose body is one unbounded repeat becomes `(X)?`, which captures differently when that repeat is lazy | `(\w+?)*?(?<=b)` matched against all of `ab`: group 1 at `0..2` / `1..2` |
| `Hir::alternation` (`lift_common_prefix`), in what `fancy-regex` hands to `regex-automata` | `XA\|XB` becomes `X(?:A\|B)` when every branch is a concatenation, which tries `B` before the other ways of matching `X` with `A` | `[ab]+b\|[ab]+c` on `abbc`: `0..4` / `0..3`; `\w+?a{1,3}?\|\w{1,}?()` on `aba`: `0..1` / `0..3` |
| the same, in the body of an atomic group (found by the fourth review) | `fancy-regex` compiles an atomic body as if nothing followed it (`Compiler::visit(child, false)`) and hands it, or a branch or trailing run of it, to `regex-automata` whole; the group then keeps the first match of the factored form | `(?>[ab]?b\|[ab]?c)` on `bc`: `0..2` / `0..1`; `(?>(?:a\|ab)c\|(?:a\|ab)b)` on `abc`: `0..2` / `0..3`, and no full match / full match; `(?=(?>(?:a?\|ab)c\|(?:a?\|ab)b)c)` on `abc`: `0..0` / `1..1` |

The facade spells its translation so that none of them applies:

- **Alternations on the automaton.** A branch that matches no haystack character
  (`[^\x00-\x7F\x{E080}-\x{E0FF}]`) ends every alternation whose branch order an
  automaton could lose. A branch that is not a concatenation stops the prefix
  factoring, and a branch that never matches changes no answer. `fancy-regex`
  hands three kinds of sub-expression that need no backtracking to `regex-automata`
  (`Compiler::visit`, `compile_concat`): a whole expression, which is the plain
  spelling; the body of a lookaround, or its trailing run; and the body of an
  atomic group, or a branch or trailing run of it, compiled as if nothing followed.
  Elsewhere the counted spelling, which ends in assertions that need backtracking,
  is handed over only in runs of fixed size, whose ways of matching all end at the
  same position. So the plain spelling guards every alternation outside
  lookarounds, and both spellings guard every alternation inside an atomic group,
  at any depth and also inside a lookahead. A lookaround body only asks whether it
  matches, which the factoring does not change, so a lookaround's own alternations
  are not guarded. Nor is an alternation whose nearest enclosing lookaround is a
  lookbehind, even inside an atomic group: its branches must keep one width, which the extra
  branch would break, and every way of matching a body of one width ends at the same
  position, so the group's choice cannot change what follows.
- **Repeats the optimizer rewrites.** Every translation is parsed with the engine's
  parser and passed through the engine's own `optimize` (with an empty node
  appended to the root, so that only the trailing-lookahead rewrite described
  under [Work bounds](#work-bounds) is left out). If the tree changes, the
  expression is translated again with every unbounded or optional repeat (`*`,
  `+`, `?`, `{n,}`, greedy or lazy) written as `(?:X*|[^...])`, an alternation with
  the same never-matching branch, which none of the rewrites matches. If the tree
  still changes, the expression is refused; no expression of the corpus or the
  fuzzing is. Only the program the engine will run is checked: the plain spelling
  when the expression needs no backtracking, and the counted spelling always.

Both spellings cost little: the guard branch is tried once each time the
alternation fails, a shielded repeat costs one backtracking branch per pass, and
each adds one atom to the automaton count. Shielding is rare; in the full corpus
before the third round's grids, only `(?:a?)?b`, `(?:a?){0,1}b` and `(?=a*a*a*)b` needed
it, and the last now exhausts the budget at 66 `a` instead of 68 (see below). The
check adds parsing time to construction (see [Construction cost](#construction-cost)).

### Refused constructs

The facade refuses with `Error::Unsupported`, rather than approximate, every
construct below. Boost accepts them, except where a row names patterns Boost rejects
too. Refusals that exist to keep work or memory
bounded are explained under [Work bounds](#work-bounds) and
[Construction cost](#construction-cost).

| Construct | Why |
| --- | --- |
| possessive quantifiers `*+ ++ ?+ {n}+` | not translated |
| `\Q...\E`, `\K`, `\G` | not translated |
| recursion `(?R)`, `(?1)`, `(?+1)`, `(?-1)`, `(?&name)`, `(?P>name)`; conditionals `(?(...)...)`; branch reset `(?\|...)`; verbs `(*...)` | not translated |
| `\Z` | Boost's start map for `\Z` leaves out `\f`, so a leading `\Z` is never tried at a form feed (`\Z` on `"\f"` matches at 1, not 0); reproducing that bug is not worth it |
| a lazy repeat with a finite maximum at least two above its minimum (`{n,m}?` with `m >= n + 2`) of a one-byte atom (a literal, `.`, a class or bracket expression) that starts the expression: nothing precedes it but group starts and ends, flag groups that keep case sensitivity, `^ $ \b \B \< \> \A \z`, comments and whole lookarounds, the expression has no backreference, no `\|` outside a lookaround precedes it or separates the alternatives of a group around it (or of the expression), and no quantifier applies to a group around it | a Boost bug, refused rather than reproduced like `\Z`. `basic_regex_creator::probe_leading_repeat` marks such a repeat as leading. When the lazy repeat extends below its maximum, Boost records the position it got to (`unwind_char_repeat`, `unwind_short_set_repeat`, `unwind_fast_dot_repeat`, `unwind_slow_dot_repeat`), and a failed start position resumes the search behind it (`match_prefix`, `find_restart_*`), skipping start positions that match: `a{1,3}?\b` on `aaaa` matches `3..4`, and `.{3,5}?\b` on `ababab` does not match, where every start position gives `1..4` and `1..6`. A flag group that changes case sensitivity adds a `toggle_case` state and an alternation or a quantifier adds a state in front, which ends Boost's walk, so `(?i)a{1,3}?\b` and `x\|a{1,3}?\b` are translated |
| a quantifier (`{1}` and `?` included) on a group that holds, at any depth, a repeat or an alternation compiled under the other case sensitivity, such as `(?i:b.+)*C`, `(?:y\|(?i:b.+))+C` or `(?:(?i)x\|y+)*`; and `\<` anywhere in an expression that switches case sensitivity (`(?i)`, `(?-i)`, `(?i:...)` or `(?-i:...)` against the options), such as `(?i)\<A` | a Boost bug (found by the fifth review), refused like `\Z`. Before matching, Boost gives every repeat and alternation a map of the bytes that can start each of its two ways on (`basic_regex_creator::create_startmaps`), which the matcher consults before it enters an iteration, leaves a repeat or tries an alternative. The walk that builds a map follows the states that can come next and recurses (`create_startmap`) at `\<`, `\>`, and at a repeat or alternation whose map is not built yet, which it reaches when it loops back to the repeat of an enclosing group. Each recursion restarts from the case sensitivity of the state whose map is built (`bool l_icase = m_icase`), or of the options for the expression's own map, not from the one in effect where the walk stands, so a literal or set further on is looked up with the wrong case and the map drops the bytes that start it: `(?i:b.+)*C` on `bcC` matches `2..3` in Boost and `(?i:bc+)+C` does not match, `(?i)\<A` does not match `A`, and `(?i)(?:a(?-i:b.+))+c` does not match `abcC`, where every start position gives `0..3`, `0..3`, `0..1` and `0..4`. A repeat state has the case sensitivity in effect at its quantifier; Boost puts the alternation state of a group's first `\|` at the start of the group, before a scoped switch takes effect, and that of a later `\|` at the start of the alternative before it, so `(?i:a\|b)*`, `(?:(?i)a\|b)*C` and `(?:(?i:a\|b)c+)+C` are translated. After `\>` only non-word bytes remain in the map, which have no case, so `(?i)a\>` is translated too |
| `\<` or `\>` anywhere in an expression in which a quantifier applies to a group holding, at any depth, a repeat or an alternation, such as `(?:\w\w+){2}\>`, `(?:\W\W??\|a)+\>`, `(\w+?)+\>` or `\<(?:\w\w*)+` | a Boost bug (found by the fifth review), refused. At `\<` (`\>`) the start-map walk recurses and then removes every non-word (word) byte from the whole map it is filling. After the walk has looped back to an enclosing repeat, that map already holds the bytes that start another iteration, which come before the assertion, so an inner repeat can no longer leave for another iteration: `(?:\w\w+){2}\>` does not match `aaaa` in Boost, `(?:\W\W??\|a)+\>` finds `2..3` in `..a` and `(\w+?)+\>` captures group 1 at `0..4` in `aaaa`, where every start position gives `0..4`, `0..3` and `3..4`. The walk does not recurse at `\b`, `\B` or a lookaround, and no other state removes bytes, so `\b` and `\B` after the same groups, and `\<\w+\>`, `(?:\w\w){2}\>` and `(?:ab)+\>`, are translated |
| an expression in which `A_free + 2 (A_looped + L_looped + R) + 2` exceeds 100, or, when it holds a buffer end (`\z`, `` \' `` or `$` under `(?-m)`), that sum plus `2 (L_free + L_looped)`: `A` counts `$` line anchors, `\<` and `\>`, `L` alternations (one per `\|`), each free outside every repeated group and looped inside one, and `R` counts the quantifiers applied to groups | Boost throws `error_complexity` at construction when `create_startmap` recurses more than `BOOST_REGEX_MAX_RECURSION_DEPTH` (100) levels deep (found by the fifth review for `$`, `\<` and `\>` chains; by its sibling hunt for alternations and nested repeated groups). The sum bounds that depth: the walk recurses one level at an anchor and two at a repeat or alternation whose map is not built; such a state lies earlier in the expression inside a repeated group, or is an alternation all of whose ways end at a buffer end, whose map stays empty; a repeat is recursed into once per map (`set_bad_repeat`); and after a walk loops back to a repeat it stays inside that repeat's group, so a chain of recursions meets each site once, except that an anchor or an alternation whose map stays empty inside a repeated group can be met a second time after the loop back. `$` × 101 then `a`, `(?:a+)*` then `$` × 99, an alternation of 51 branches in `(?:(...)x?)+`, 40 nested `(?:...a+$)+` and `$` × 80 before `(?:\b\|\B)` × 12 and `\z` are Boost errors. The bound refuses some expressions Boost compiles (`$` × 99 or 100 then `a`, `a` then `$` × 101, `x*` then `$` × 99, 30 nested `(?:...$)+`, 50 branches in a repeated group); `$` × 98 then `a`, 49 branches and 20 nested `(?:...$)+` are translated. 7,700 random expressions of anchors, alternations, lookaheads and nested repeated groups whose sum is between 60 and 100 all compile in Boost |
| a lookbehind holding more than 1,024 alternations | Boost's `calculate_backstep` stacks every alternation on the path it walks and rejects the lookbehind when it would stack one more than `BOOST_REGEX_MAX_BLOCKS` (1,024): `(?<=` then `(?:\|)` × 1,026 then `)a` is an error in Boost. One alternation of 1,030 branches, which Boost accepts, is refused as well |
| `\x` escapes whose digits start with white space, `+`, `-` or `0x`/`0X`, such as `\x{+41}`, `\x{ 41}`, `\x{0x41}`, `\x+4`, `[\x 4]` or `\x0x` | Boost reads the digits with `std::istream >> std::hex` (`cpp_regex_traits::toi`), which skips white space and accepts a sign and a `0x` prefix: the first five are valid in Boost (`\x-0` is a NUL byte), where the facade reported a syntax error before the fifth review, and `\x0x` and `\x{0x}` are errors, which it compiled. Every other `\x` escape reads the same in both |
| octal `\0...`, `\c` control escapes, `\g` and `\k` backreferences, multi-digit backreferences `\10`, `\p`/`\P` properties, `\N`, `\R`, `\X`, `\C` | not translated |
| escape letters Boost reads as a literal or a locale class (`\y`, `\j`, `\l`, `\u`, `\L`, `\U`, `\E` ...), the same inside a bracket expression, where `[\l]`, `[\u]`, `[\L]` and `[\U]` are the lower and upper classes and `[\p]`, `[\y]`, `[\0]` are literals, and `[a-\d]` and `[\V]` | not translated |
| a `[[:name:]]` whose name no lookup finds and that is not all letters, such as Boost's word-boundary spellings `[[:<:]]` and `[[:>:]]`, `[[::]]` or `[[:al:pha:]]`, and a `[[:` that never closes with `:]` | `[[:<:]]` and `[[:>:]]` are `\<` and `\>` in Boost when the rest of the set is empty (`parse_inner_set`); the others are `error_ctype` or `error_brack` there. Refused rather than reported as an error, because Boost compiles the first two. An all-letters name no lookup finds stays an error, as Boost's `error_ctype` |
| collating elements `[[.a.]]` and equivalence classes `[[=a=]]`, the elements at either endpoint of a range (`[A-[.a.]]`, `[[.a.]-z]`, `[\n-[.-.]]`, `[!-[.].]]`) included | locale-dependent. `get_next_set_literal` opens a collating element wherever a `[` is followed by a `.`, at a range endpoint as well as at the start of a set item, and `lookup_collatename` then decides the endpoint: `[A-[.a.]]` is the range `A` to `a` and matches `[`, `\`, `]`, `^`, `_` and `` ` ``, `[A-[.ab.]]` is `error_collate` and `[A-[.tab.]]` is `error_ctype`. Found at the range end by the sixth review, where the facade read the `[` as a literal; `[` followed by anything else stays one (`[A-[x]`, `[A-[=a=]]`, `[A-[:alpha:]]` are translated) |
| Boost's own `[[:unicode:]]` class, in any case spelling and negated | `get_default_class_id` knows a `unicode` name whose mask (`mask_unicode`) no `char` can carry, so `[[:unicode:]]` matches nothing and `[[:^unicode:]]` every byte. Found by the sixth review, where the facade reported an unknown class name; refused rather than translated to an empty and a full set |
| the `x` modifier (`(?x)`) | not translated |
| group names with a character outside `[A-Za-z0-9_]` (Boost accepts `(?<n-x>...)`) | not translated |
| non-ASCII pattern bytes outside comments, and `\x` escapes above `0x7F` | a pattern byte above `0x7F` would have to split a UTF-8 character; on platforms where `char` is unsigned Boost also accepts `\x{80}`-`\x{FF}`, on signed-`char` platforms it rejects them |
| a capturing group inside a lookaround or an atomic group | Boost runs these as independent sub-matches and does not restore captures made inside them when the surrounding match backtracks past them: `(?!(a))b` or a failed branch after `(?>...()...)` leaves the group set. `fancy-regex` restores them. |
| a repeat allowing more than one repetition of a group that can match empty and captures | Boost records a final empty iteration: `(a*)+` on `"a"` gives group 1 = `1,1`, `fancy-regex` gives `0,1` |
| a repeat allowing more than one iteration (`*`, `+`, `{n,}`, and `{n}` or `{n,m}` above 1, greedy or lazy) of a group that can match the empty string, and a bounded one (`{n}`, `{n,m}` above 1) of a backreference that can | Boost ends a repeat as soon as an iteration matches the empty string, even below the minimum: it accepts the iteration and takes the exit (`repeater_count::check_null_repeat`, `perl_matcher::match_rep`). The engine differs either way. A bounded repeat iterates to its bound: `(?:\b\|a){2}b` matches `ab` in the engine and not in Boost, and `(?:a{0}\b){999999999}` ran for 16.9 s on one byte. An unbounded one (`RepeatEpsilon`) fails the empty iteration and backtracks into the group's other alternatives: `(?:b?\|a)*` on `ba` matches `0..2` in the engine and `0..1` in Boost, `(?:(?:a\|b?)*?)+b` on `abab` `0..4` and `0..2`. A backreference to a group that is closed where the backreference stands matches the same text in every iteration of its own repeat, so its unbounded repeats agree and are translated (a backreference inside its group is refused, next row) |
| a backreference inside the group it refers to, such as `(a\|\1b)x`, `(a)(b\|\2a)x` or `(?:(a\|\1b)c)+` | Boost saves a group's span when the group starts, sets only its start, and restores the saved span only when the match backtracks past that start (`perl_matcher::match_startmark`, `match_results::set_first`). A backreference inside the group therefore compares against what an abandoned attempt left: the span of an alternative that matched and failed later, so `(a\|\1b)x` matches all of `abx` and `(a\|\1a)b` all of `aab`; or, in a later iteration, the range from the new start to the previous end, which is empty when the iterations are adjacent (`(a\|b\1)+` matches all of `ab`) and never matches otherwise. `fancy-regex` treats the group as unset in the first case (no match, `1..3`) and in the second slices the haystack from the new start to the previous end, which panicked on `(?:(a\|\1b)c)+` and `acbc` |
| a repeat of a group that only asserts positions, such as `(?:$)+` | Boost never matches `(?:$)+` at all |
| a repeat of a modifier group that switches case sensitivity, such as `(?i)+a` | Boost undoes the switch when the empty repetition is abandoned, so `a` stays case-sensitive |
| a repeat allowing more than one repetition inside an atomic group or a negative lookaround, such as `(?>a+)b`, `(?!.*x)` or `(?<!(?=.*b)a)` | work bound: the engine discards the backtracking branches pushed there without counting them |
| a lookbehind wider than `MAX_LOOKBEHIND_WIDTH` (255) bytes | work bound |
| repeat bounds whose value exceeds `MAX_REPEAT` (999,999,999), and a repeat whose shortest match, multiplied out over nested repeats, exceeds it | work bound; the product also overflows `fancy-regex`'s analyzer, which panicked in builds with overflow checks on `(?:(?:a{999999999}){999999999}){999999999}` |
| nesting deeper than `MAX_GROUP_DEPTH` (48), a translation whose groups could nest 64 levels deep, and a translation longer than `MAX_TRANSLATED_BYTES` (512 KiB) | work and memory bounds, and the limit of `fancy-regex`'s parser (`MAX_RECURSION`, 64), which the translation reaches before the pattern does: a case-insensitive letter, `.` or a line anchor adds levels, and so does every counted and every shielded repeat (see [Engine rewrites](#engine-rewrites)). The translator counts, along every path and as if every repeat were counted and shielded, one level for a one-byte atom or a backreference, two for an assertion, one for a group (two for a positive lookahead), one more for a quantifier on a lookaround or a zero repeat of a capturing group, and otherwise, for a quantifier, one level above at least two when its minimum is 2 or more and one more when it is unbounded or optional; two more for the programs' own groups. `(?:` × 19, `a{2,}`, `){2,}` × 19 is translated and 20 levels are refused; none of 300,000 random nested expressions reached the engine's limit (see Evidence). A pattern longer than `MAX_PATTERN_BYTES` (64 KiB) is `InvalidValue` |

None of these constructs occurs in a pinned OpenMS expression: every family
derived from the OpenMS sources compiles in full, and none of the rules added in
the third to sixth rounds refuses one of them. The sixteen families and their
pattern counts, all with 0 refusals in both corpora: `ENZ` 30, `RNA` 17, `CLASS`
24, `MSP` 14, `CHROM` 12, `NATIVE` 8, `MRM` 5, `ANNOT` 4, `PERCOUT` 3, `TITLE` 3,
`DECOY` 2, `MZTAB` 2, `FRAG` 1, `IDX` 1, `LOOKUP` 1, `PEPXML` 1. How often each
construct was hit in the corpus is under Evidence. A translation the engine itself refuses would be
`Unsupported` too; no corpus pattern is. (`(?:(?:a{999}){999}){999}`, whose
automaton the engine refused before the second round, is searched on the
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
The branch that guards an alternation and the extra branch of a shielded repeat
(see [Engine rewrites](#engine-rewrites)) count as one atom each.
`\w{0,9999}b` has 10,000, the largest expression of the pinned sources that needs
no backtracking (`TransitionGroupPicker:PeakPickerChromatogram:(.+)`) 46. A search
runs on the automaton only when the expression has at most
`MAX_AUTOMATON_ATOMS` (4,096) atoms divided by the budget's scale factor below:
4,096 up to 64 KiB, 1,024 up to 256 KiB, 256 up to 1 MiB and 64 beyond. Every
other search runs on the backtracking machine with the counted spelling, where the
budget bounds it. So an automaton search does at most about 2.7 × 10⁸ atom steps
up to 4 MiB, and 64 per byte beyond. Measured with expressions that defeat the
lazy DFA, at each tier's limit, over random `a`/`b` haystacks (as recorded in the
second round, on a less loaded machine than the transcoded rows below):

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

An atom step is not a fixed amount of time. A haystack byte above `0x7F` is three
bytes of the transcoded haystack the engine searches, and a class that accepts such
bytes (`.`, `\W`, `[^b]`) is several automaton states, so the same expressions over
haystacks with such bytes are slower. Measured in the third round, release build, the same
machine for every row (so the `a`/`b` rows, repeated from the table above, show its
speed against the second round's):

| Expression | Atoms | Haystack | Search |
| --- | ---: | --- | ---: |
| `a.{0,2046}a.{0,2046}c` | 4,095 | 64 KiB, random `a`/`b` | 2.72 s, no match |
| `a.{0,2046}a.{0,2046}c` | 4,095 | 64 KiB, random `a`/`0xE9` | 4.51 s, no match |
| `(?:.{0,10}a){0,372}c` | 4,093 | 64 KiB, random `a`/`0xE9` | 3.59 s, no match |
| `.{0,4095}b` | 4,096 | 64 KiB of `a` | 3.53 s, no match |
| `.{0,4095}b` | 4,096 | 64 KiB of `0xE9` | 7.63 s, no match |
| `[^b]{0,4095}c` | 4,096 | 64 KiB of `0xE9` | 7.72 s, no match |
| `\W{0,4095}b` | 4,096 | 64 KiB of `0xE9` | 8.17 s, no match |
| `a.{0,510}a.{0,510}c` | 1,023 | 256 KiB, random `a`/`0xE9` | 4.82 s, no match |
| `(?:.{0,10}a){0,93}c` | 1,024 | 256 KiB, random `a`/`0xE9` | 3.97 s, no match |
| `a.{0,126}a.{0,126}c` | 255 | 1 MiB, random `a`/`0xE9` | 3.81 s, no match |
| `(?:.{0,10}a){0,23}c` | 254 | 1 MiB, random `a`/`0xE9` | 3.58 s, no match |
| `a.{0,30}a.{0,30}c` | 63 | 16 MiB, random `a`/`0xE9` | 15.7 s, no match |
| `(?:.{0,6}a){0,9}c` | 64 | 16 MiB, random `a`/`b` | 10.7 s, no match |
| `(?:.{0,6}a){0,9}c` | 64 | 16 MiB, random `a`/`0xE9` | 17.3 s, no match |

So bytes above `0x7F` make an automaton search up to about twice as slow at the same
atom count; the ceilings above scale accordingly. Boost 1.92 stops the 64 KiB
searches with `error_complexity` after 79 ms.

The bound is per search. A token iterator runs one search per match, each with its
own bound (as each has its own budget, below), so a whole iteration may take the
bound once per match: `a.{0,2045}a.{0,2046}c|a` (4,096 atoms, on the automaton) over
4 KiB of random `a`/`b` yields 2,073 tokens in 44 s, where Boost yields the same
2,073 tokens in 1.0 s; with one atom more (`a.{0,2046}a.{0,2046}c|a`) the first
search runs on the backtracking machine and stops with the budget after 15 ms.

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
refused since the second round for their answers (see [Refused constructs](#refused-constructs)).

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
has its own budget, scaled by the whole haystack, and its own automaton bound, so
an iteration's total work is the bound of one search times the number of matches.

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
lookbehind widths between two backtracking steps, plus the length of the group a
backreference names, so every search's work is bounded by its budget times that.
The translated pattern can be long: `a` × 65,535 then `b` spends one step per start
position but compares up to 64 KiB at each (1.19 s over 1 MB, above).

**Backreferences compare their group.** A case-sensitive backreference compares
bytes up to the first one that differs, as Boost's does. A case-insensitive one
does not stop there: `fancy-regex`'s `matches_literal_casei` first compares exactly
and, when that fails, checks the whole group (`is_ascii`, or UTF-8 decoding for a
transcoded haystack) before it compares without case, so every attempt costs the
group's length even when the first byte differs. A group that grows lazily across
the haystack then makes one search quadratic within its budget where Boost's stays
linear. When the backreference matches, the comparison costs the group's length in
both engines, and Boost's search is quadratic too. Measured this round, release
build, one search, Boost 1.92 on the same haystacks (its time includes reading the
input):

| Expression | Haystack | Facade | Boost |
| --- | --- | --- | --- |
| `([a-c]*?)(?i)[^a](a\1)` | 64 KiB, random `a`/`b` | 0.023 s, match `2..6` | 0.014 s, match `2..6` |
| `([a-c]*?)(?i)[^a](a\1)` | 256 KiB, random `a`/`b` | 0.29 s, match | 0.013 s, match |
| `([a-c]*?)(?i)[^a](a\1)` | 1 MiB, random `a`/`b` | 4.4 s, match | 0.041 s, match |
| `([a-c]*?)(?i)[^a](a\1)` | 4 MiB, random `a`/`b` | 70 s, match | 0.16 s, match |
| `([a-c]*?)[^a](a\1)` (case-sensitive) | 1 MiB / 4 MiB, random `a`/`b` | 0.050 s / 0.20 s, match | not measured |
| `(a*?)\1b` | 16 KiB of `a` | 0.055 s, budget | 60 s, `error_complexity` |
| `(a*?)\1b` | 32 KiB of `a` | 0.095 s, budget | 126 s, `error_complexity` |
| `(a*?)\1b` | 64 KiB of `a` | 0.16 s, budget | 246 s, `error_complexity` |
| `(a*?)\1b` | 256 KiB of `a` | 2.5 s, budget (4,000,000 steps) | not measured |
| `(a*?)\1b` | 1 MiB of `a` | 48 s, budget (16,000,000 steps) | not measured |
| `(a*?)\1b` | 4 MiB of `a` | 916 s, budget (64,000,000 steps) | not measured |

The case-sensitive `(a*?)\1b` is the per-tier ceiling of a matching comparison: the
budget grows with the haystack and so does the group, so the time to the error grows
with about the square of the haystack length from one tier to the next. A
caller-supplied expression with a backreference to a group that can grow across the
haystack, and any case-insensitive backreference to one, should expect seconds to
minutes on haystacks of megabytes.

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
| `(?=a*a*a*)b` | `a` × 300, `c` | 0.09 ms, no match | 16 ms, budget | 0.05 ms, no match | 65 `a` (67 before its repeats were shielded) |
| `(?=.*a).*b` | `a` × 4,000 | 9.8 ms, no match | 8.5 ms, budget | 0.02 ms, no match | 142 bytes |
| `(?=[A-Z]*K)[A-Z]+R` | `A` × 12,000 | 19 ms, no match | 13 ms, budget | 57 ms, no match | 1,411 bytes |

A caller-supplied filter of this shape should expect the error on inputs of a few
dozen bytes to a few kilobytes.

Measured in a release build on aarch64 macOS, the facade before the counted
spelling (854d996, as recorded in the first round), with it (measured in the second round),
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
bytes per counted repeat, 4 per positive lookahead and 12 in all; the guard branch
adds 28 bytes per alternation outside lookarounds in the plain spelling and per
alternation inside an atomic group in both; and a translation with shielded
repeats (see [Engine rewrites](#engine-rewrites)) adds 32 per unbounded or optional
repeat. `MAX_TRANSLATED_BYTES` applies to both spellings. Before compiling,
construction parses the one or two translations the engine will run and passes them
through the engine's optimizer. The engine
compiles three programs per budget (search, full match, and the non-empty retry of
the token iterator) and holds an automaton per lookaround and per run of
sub-expressions it delegates, so memory grows with the translation; an expression
searched on the automaton also builds that automaton, whose size grows with its
atoms. `MAX_TRANSLATED_BYTES` (512 KiB) caps the translation, and
`MAX_AUTOMATON_ATOMS` the automaton: an expression with more atoms is compiled only
for the backtracking machine, whose program grows with the translation, not with
the repeat bounds. Peak resident set size of constructing one expression near the
caps, release build, measured in the third round on one machine for commit cf594ef
and for cfbaf1a (this round's guards apply only inside atomic groups, which none of
these patterns has):

| Pattern | Pattern bytes | Construction, cf594ef / cfbaf1a | Peak RSS, cf594ef / cfbaf1a |
| --- | ---: | ---: | ---: |
| `^` × 12,190 | 12,190 | 87 ms / 105 ms | 65 MB / 77 MB |
| `^a\|` × 11,000 then `b` | 33,001 | 100 ms / 125 ms | 82 MB / 97 MB |
| `a` × 58,252, `icase` (the longest that fits) | 58,252 | 60 ms / 89 ms | 26 MB / 32 MB |
| `a` × 65,536 | 65,536 | 22 ms / 39 ms | 18 MB / 27 MB |
| `(?<=a)` × 10,922 | 65,532 | 10 ms / 14 ms | 18 MB / 20 MB |
| `^` × 12,190, then one search over 5 MB (scaled programs compiled) | 12,190 | 162 ms / 185 ms in all | 98 MB / 111 MB |
| `a` × 65,536, then one search over 5 MB of `a` | 65,536 | 41 ms / 57 ms in all | 27 MB / 32 MB |

The difference is the parse and optimizer pass over the translation. In the second
round, on a less loaded machine, the first five constructions took 58, 69, 42, 17 and
7 ms.

`a` × 58,253 with `icase` translates to more than 512 KiB and is refused before
the engine sees it. At c93fefd the 65,536-atom literal was compiled as an automaton
too (37 MB) and took 10 s to find its match at the start of 5 MB of `a`; since the
second round it runs on the backtracking machine (12 ms then, 18 ms on the third round's
machine).

Boost's own construction is exponential in two shapes the facade compiles in linear time
(measured with the oracle driver, fifth round). `calculate_backstep` walks every
combination of the alternatives in a lookbehind: 15, 20, 22 and 24 alternations `(?:|)` in
a row took 0.00, 0.01, 0.03 and 0.11 s, and 1,025 did not finish (1,026 are rejected at
once, see [Refused constructs](#refused-constructs)). And the start map of an
alternation all of whose ways end at a buffer end stays empty, so every later walk
recurses into it again: 40 or more `(?:\b|\B)` in a row before `\z` did not finish within
20 s. The facade answers both.

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
- **Optimizer safety net.** Construction fails with `Unsupported` ("a shape the
  engine's optimizer rewrites") if the engine's optimizer still rewrites the
  translation after its repeats were shielded (see
  [Engine rewrites](#engine-rewrites)). No pattern of the corpus or of the fuzzing
  hits it.
- **Name hashing.** Boost looks names up by a hash and would conflate two names
  with the same hash; the facade compares names.
- **Speed.** Expressions on the backtracking engine are slower than Boost (above).
  The counted spelling costs little on OpenMS's expressions; release build, same
  haystacks, facade at commit 854d996, at c93fefd and now: `(?<=[KR])(?!P)` tokens
  over a 10 MB random protein 255, 258 and 263 ms, `(?<=[KRX])` 251, 255 and
  263 ms, `(?=[DBX])` 55, 211 and 213 ms (its lookahead body runs on the
  backtracking engine since c93fefd); at 854d996 and c93fefd, `=(?<SCAN>\d+)$`
  2.4 and 2.5 ms and `^(?:Name|NAME): (.+)` with `no_mod_s` 2.6 and 2.4 ms (not
  re-measured; the second round changed their spelling only by a second trailing `(?=)`).
  The generated protease table stays the fast path for the pinned enzymes.

## Checked boundaries and evidence

**Tier 1, executed differential.** The oracle lives outside the repository in
`../oracle/boost-regex/`: `driver.cpp` (compiled with Apple clang 21.0.0,
`-std=c++17 -O2`, against the header-only Boost.Regex 1.92 in
`/opt/homebrew/include`), `gen.py` (deterministic, seed 20260913, reads the pinned
checkout), `build.sh <checkout or worktree>` (builds the driver, regenerates both
corpora, runs Boost and copies the fixture pair into `tests/data`; it takes the
compiler from `CXX`) and
`refusals.py`, plus `fuzz.py` for the random expressions described below. The
manifest records their sha256, the Boost headers' sha256 and the full corpus's
sha256. Running `build.sh` in a scratch directory holding only the five oracle
sources reproduced both corpora, both Boost outputs and the committed fixture pair
byte for byte, and every sha256 in the manifest verifies (39 pinned sources, 2
fixtures, 7 external artifacts, 10 Boost headers).

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
| Full corpus (`../oracle/boost-regex/corpus_full.txt`) | 134,075 | 11,427,044 | 9,638,702 | 6,531,674 | 0 |
| Committed fixture (`tests/data/boost_regex_corpus.txt`) | 2,416 | 228,565 | 198,124 | 162,714 | 0 |

Boost answers every case of a pattern it compiles; compared cases exclude the
patterns the facade refuses. The fixture keeps every pattern except the systematic
probe grids of the full corpus (`UNB`, `LEAD`, the `BR` grids, `OPTNEST`, `OPTCAT`,
`ALTPFX`, `ATOMALT`, `ATOMCTX`, `ATOMBR`, `SMCASE`, `SMWORD`, `SMDEPTH`, `SMBR`,
`HEX`, `COLL`, `PEXT` and `CLSNAME`), every fixed input (edge, structured and
source-derived) and every run except
extra operations marked full-only, and cuts seeded random inputs to a prefix per input
set. It is 1,123,279 bytes (287,867 corpus and 835,412 Boost output). A scratch build
with `MAX_AUTOMATON_ATOMS` set to 0, which searches every pattern with the counted
spelling on the backtracking machine, also compares all 6,531,674 cases of the full
corpus with 0 mismatches, and so does a build with overflow checks and debug
assertions, without a panic. Rerunning `gen.py` and `refusals.py` reproduces the
committed fixture pair, the full corpus and its Boost output, and the test constants
byte for byte.

Commit 9597f21, the previous round, run over this full corpus in the same harness,
refuses 38,507 patterns, gives 55,968 different answers (none of them panics) and
differs from Boost in 2,136 compile outcomes. The answers are 54,934 in `COLL` and
1,034 in `ADV` (the sixth review's probes); the compile outcomes 1,186 in `COLL`
(ranges whose end is a collating element, in both directions: Boost compiles and it
reported an invalid range, or Boost rejects and it compiled), 643 in `CLSNAME`
(class names it read as unknown), 230 in `PEXT` (`(?P` groups it read as a syntax
error) and 77 in `ADV`. Over the fixture it gives 362 different answers and 77
different compile outcomes, all in the new `ADV` probes. Earlier rounds: e6c581f
gave 12,119 different answers and 795 different compile outcomes over the fifth
round's full corpus (127,161 patterns), cfbaf1a 12,399 different answers over the
fourth round's (84,017 patterns), and cf594ef 10,206, 1,969 of them panics, over the
third round's (67,916 patterns).

Per family in the committed fixture (the families other than `SYN`, `ADV` and
`FUZZ` come from OpenMS sources; `ADV` holds the probes of the first two
independent reviews and of the work bounds, the third review's backreference and
engine-rewrite probes, the fourth review's atomic-alternation probes, the fifth
review's start-map, recursion-limit and `\x` probes, and the sixth review's
collating-element, `(?P` and class-name probes; some of them, such as `\<\w+\>`,
`\x{41}`, `[[.a.]]`, `[[=a=]]` and `(?P<n>a)`, repeat `SYN` patterns and are counted
there):

| Family | Patterns | Refused by the facade | Refused by both | Compared cases |
| --- | ---: | ---: | ---: | ---: |
| `ADV` | 599 | 179 | 49 | 23,162 |
| `ANNOT` | 4 | 0 | 0 | 1,188 |
| `CHROM` | 12 | 0 | 0 | 1,344 |
| `CLASS` | 24 | 0 | 0 | 2,768 |
| `DECOY` | 2 | 0 | 0 | 516 |
| `ENZ` | 30 | 0 | 0 | 17,269 |
| `FRAG` | 1 | 0 | 0 | 338 |
| `FUZZ` | 1,111 | 246 | 422 | 7,088 |
| `IDX` | 1 | 0 | 0 | 231 |
| `LOOKUP` | 1 | 0 | 0 | 628 |
| `MRM` | 5 | 0 | 0 | 615 |
| `MSP` | 14 | 0 | 0 | 2,310 |
| `MZTAB` | 2 | 0 | 0 | 514 |
| `NATIVE` | 8 | 0 | 0 | 2,826 |
| `PEPXML` | 1 | 0 | 0 | 171 |
| `PERCOUT` | 3 | 0 | 0 | 354 |
| `RNA` | 17 | 0 | 0 | 6,395 |
| `SYN` | 578 | 87 | 79 | 94,079 |
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
their allowed neighbours; lookbehind widths at 255 and 256; the counted
spelling of expressions that need backtracking; backreferences inside the group they
refer to (refused: the third review's examples and the shape on which the engine
panicked) beside backreferences to closed, later, unset, repeated, duplicate-named,
zero-repeated and case-insensitive groups and into lookarounds and atomic groups
(compared, over inputs with bytes above `0x7F`); the expressions the engines'
rewrites answered differently (compared); and alternations inside atomic groups
whose branches share a prefix that can match in several ways, among them the fourth
review's examples, in lookaheads and negative lookaheads, in a lookahead inside a
lookbehind, with hard parts before and after, nested and repeated, beside
backreferences (case-insensitive, over bytes above `0x7F`) and after a lazy leading
repeat, and alternations of one width in lookbehinds (compared); and the fifth review's
start-map examples (refused) beside the same groups without the repeat or with the
switch where Boost's walk keeps it, `\b`, `\B` and `$` after the same repeated groups
and `\<`, `\>` without an inner repeat (compared), `\x` escapes with white space, a sign
or `0x` (refused, or rejected by both) beside plain ones, expressions on both sides of
the start-map recursion bound (compared, refused, or rejected by both with Boost's
`error_complexity`), and translations nested at the engine's parser limit; and the
sixth review's examples: collating elements as the end of a range, negated, beside
other items and followed by a second dash (refused, or rejected by both with
`error_collate`, `error_range` or `error_ctype`) beside the literal `[` of
`[A-[x]`, `[A-[]`, `[A-[=a=]]x` and `[A-[:alpha:]]x` (compared), `(?P` groups in
every option-group, recursion and Python spelling (compared, refused for `(?P>name)`
and the `x` modifier, or rejected by both), and class names in four case spellings
with Boost's one-letter aliases and `[[:unicode:]]` (compared, refused, or rejected
by both).

The full corpus adds grids. `UNB` (19,183 patterns, the second review's probe)
repeats every body of two nullable or non-nullable parts, concatenated or
alternated, with `*`, `+`, `*?`, `+?`, `{1,}` and `{0,}`, before four suffixes,
over 19 inputs with search, full match and tokens: 15,346 refused, 219,087 cases
compared (the count includes the runs that later grids add to the patterns they share
with it, here and below). `LEAD` (18,365 patterns) puts eight atoms under ten greedy, lazy, bounded
and unbounded quantifiers before ten followers, inside twenty prefixes and wrappers
(anchors, lookarounds, alternations before and after, capturing and non-capturing
groups, case-changing and case-keeping flag groups, a quantified lookahead, a
trailing backreference), with and without `icase`: 2,506 refused, 573,294 cases
compared (the count includes the `BRLEAD` patterns that coincide with `LEAD`'s).

The third round adds grids of backreference siblings. Each puts eight group bodies
(`a`, `ab`, `a?`, `a*`, `a|b`, `ab|a`, `a|`, `a+?`) under seven quantifiers on the
backreference (none, `?`, `*`, `+`, `{2}`, `*?`, `{0,2}`) before three suffixes, in the
templates of its family, over 26 inputs with search, full match and the token
iterator (-1), unless the row says otherwise; patterns shared with an earlier grid
are counted there:

| Family | Templates | Patterns | Refused | Refused by both | Compared cases |
| --- | --- | ---: | ---: | ---: | ---: |
| `BROPEN` | a backreference inside the group it refers to (`(X\|\1Qb)`, `(\1Qb\|X)`, `(X\1Q)`, `((X)\|\1Q)`, `(X)(b\|\2Qa)`, `(?:(X\|\1Qb)c)+` and 7 more) beside closed-group neighbours (`((X)\2Q)`, `(X)\1Q`, `(?:(X)\|b)\1Q`, `(?:x\|(a)(?:X\|\1Q))`, `((X)\|\2Qb)x`) | 3,003 | 2,246 | 0 | 74,274 |
| `BRFWD` | a backreference to a later group (`\1Q(X)`, `(?:\1Qx\|(X))+`, `(\2Qa\|(b))` and 5 more) | 1,195 | 761 | 0 | 42,042 |
| `BRUNSET` | a backreference to a group unset in the current alternative (`(?:(X)\|b)\1Q`, `(X)?\1Q`, `(?:(X)\|b\1Q)+` and 7 more) | 1,510 | 327 | 0 | 100,542 |
| `BRREP` | groups under `+`, `*`, `{2}` and `*?`, then backreferenced (`(X)R\1Q`, `((X)\|b)R\2Q`, `(?:(X)b\|ac)R\1Q` and 3 more) | 3,863 | 1,269 | 0 | 202,566 |
| `BRLOOK` | backreferences in lookaheads, lookbehinds (nested lookarounds included) and atomic groups, and groups in them (`(X)(?=\1Q)`, `(X)(?>\1Q)b`, `(X)(?<=(?=\1Q).)`, `(?>(X))\1Q` and 19 more) | 3,802 | 1,603 | 335 | 146,562 |
| `BRDUP` | repeated group names with numbered backreferences (`(?<n>X)\|(?<n>b)\2Q`, `(?<n>X)(?<n>\1Qb)` and 7 more), with named lookups | 1,510 | 201 | 0 | 102,102 |
| `BRZERO` | groups repeated `{0}` or `{0,0}`, then backreferenced (6 templates) | 1,007 | 108 | 0 | 70,122 |
| `BRNULL` | eight bodies that can match empty (`a?`, `\|a`, empty, `a{0}`, `\b`, `(?=a)` ...) under six quantifiers in 9 templates (`(X)\1Q`, `(X)(?:\1b)Q`, `(?:(X)\1)Q` ...) | 1,212 | 420 | 0 | 61,776 |
| `BRICASE` | 11 atoms under 5 quantifiers in 7 templates mixing `icase`, `(?i)`, `(?-i)` and `(?i:...)` around the group and the backreference, with and without the `icase` flag, over 29 inputs of case pairs and bytes above `0x7F` (`é`/`É`, `ß`, `ı`/`İ`, `ſ`, Kelvin sign, Cyrillic), with search and tokens | 758 | 14 | 0 | 43,152 |
| `BRLEAD` | four atoms under five lazy quantifiers before four followers in 10 templates with backreferences (`R F(a)\1`, `()R F\1`, `(R)F\1` ...), with and without `icase`, over 15 inputs with search and tokens (0) | 1,480 | 0 | 0 | 45,600 |

And grids of the shapes the engines' rewrites changed (see
[Engine rewrites](#engine-rewrites)): `OPTNEST` (1,823 patterns; `(\w+?)+$` is
now an `ADV` probe) repeats a group holding one repeat, `(?:XI)O` and `(XI)O` for four
atoms, eight inner and seven outer quantifiers and four suffixes, plus `(XI)O\1`: 636
refused, 93,783 cases compared. `OPTCAT` (3,316) puts seven edge repeats around seven optional middles
(`XMY`, `(?:XM(?:x|Y))`) and repeats `(?:P(?:MR)?)`: 144 refused, 211,744 compared.
`ALTPFX` (3,692) gives the branches of an alternation a common prefix, 11 prefixes
and 8 tails in six templates: 144 refused, 234,168 compared.

The fourth round adds grids around alternations inside atomic groups, over 30
inputs (`ATOM`: short strings of `a`, `b`, `c`, `x`, `d` and a space, case pairs, and
`é` as UTF-8 and as a lone `0xE9`) with search, full match and the token iterator
(-1 and 0). Each gives the branches of an alternation a prefix that can match in
several ways (`[ab]?`, `a??`, `(?:a|ab)`, `(?:ab|a)`, `(?:|a)`, `\w?`, `.?`,
`(?:a?|ab)`, `(?i:a|ab)`, `a{0,1}`, `(?:a|b|ab)`, `(?:a(?:b|))`; eight of them in
`ATOMCTX` and `ATOMBR`) and a pair of tails (`b`/`c`, `c`/`b`, `bc`/`b`, empty/`b`,
`\b`/`b` and, in `ATOMALT`, `b`/empty, `(?=c)`/`b`, `c$`/`b`):

| Family | Templates | Patterns | Refused | Refused by both | Compared cases |
| --- | --- | ---: | ---: | ---: | ---: |
| `ATOMALT` | five body forms (`XA\|XB`, `XA\|XB\|X`, `(?:XA\|XB)d?`, `XA\|y\|XB`, `(?:XA\|XB)\|z`) in 24 contexts: `(?>S)`, followed by `c?`, `c` or `$`, after `x?` or `\b`, with `\b` or a lookahead inside before or after the body, captured, repeated `{2}`, in a lookahead, a negative lookahead and a lookahead inside a lookbehind, nested, beside a backreference, alternated, under `(?i)` outside and inside; three contexts also with `icase` | 12,951 | 168 | 480 | 1,107,270 |
| `ATOMCTX` | 21 placements (in a lookahead or negative lookahead inside a lookbehind, repeated `+`, `+?`, `*`, `{1,2}`, `??`, inside repeated groups and alternations, nested atomic groups, doubly nested lookaheads, after a lazy leading repeat or an optional group, with an empty branch), with and without `icase`, plus 12 lookbehinds holding atomic alternations of one width | 1,683 | 180 | 0 | 135,270 |
| `ATOMBR` | 18 templates with backreferences: to the group the atomic group is in, to later groups, to groups unset in the current alternative, to repeated groups, to duplicate names, case-insensitive (`(.)(?i)(?>\1XA\|\1XB)`), to groups that can match empty, from a lookahead, and after a lazy leading repeat with and without a backreference; with and without `icase` | 1,440 | 160 | 0 | 115,200 |

The refusals are the rules of earlier rounds: 268 repeats of more than one
iteration of a group that can match the empty string (repeated atomic groups whose
body can), 160 lazy leading repeats and 80 backreferences inside the group they
refer to. The 480 patterns Boost refuses put an empty lookahead `(?=)` in the body.
The backreference grids of the third round (`BROPEN`, `BRFWD`, `BRUNSET`, `BRREP`,
`BRLOOK`, `BRDUP`, `BRZERO`, `BRNULL`, `BRICASE`, `BRLEAD`, above) cover
backreferences to open, later and unset groups, into and out of lookarounds and
atomic groups, to repeated and duplicate-named groups, under `icase` over bytes
above `0x7F`, beside lazy leading repeats and to groups that can match empty; they
were rerun with this change and still compare equal. Beyond the corpus, the review's
own atomic probes (7,020 patterns, 463,320 cases; and 25,000 random shared-prefix
alternations in every group kind, 8,312,448 cases) and two scratch explorations of
this round (33,600 patterns in the `ATOMALT` shape with all contexts under both
flags, 2,659,104 cases; 3,060 patterns in the `ATOMCTX` and `ATOMBR` shapes, 241,800
cases) compare with 0 mismatches, where cfbaf1a gave 3,287, 4,029, 19,726 and 3,744
different answers.

The fifth round adds grids around Boost's start maps, its start-map recursion
limit and `\x` escapes. `SMC` has 22 inputs of `b`, `c`, `C`, `x`, `y` and `d` in case
pairs, `SMW` 21 of `a`, `A`, `.` and spaces with a lone `0xE9`, `SM` 19 short ones, and
`HEX` 21 of control bytes, `x`, `{`, `}`, signs, digits and white space:

| Family | Templates | Patterns | Refused | Refused by both | Compared cases |
| --- | --- | ---: | ---: | ---: | ---: |
| `SMCASE` | eight bodies with an inner repeat or alternation (`b.+`, `bc+`, `b(?:c\|d)`, `b(?:c\|C)`, `b\w*?`, `b[a-z]{1,2}`, `(?:y\|b.+)`, and `b` as a control) inside eight switches (`(?i)X`, `(?-i)X`, `(?i:X)`, `(?-i:X)`, `X(?i)`, `X(?-i)`, `(?:(?i)X)`, `(?:X)`), in seven group templates (plain, alternated before and after, after a literal, captured, nested with an inner quantifier, followed by `c?`) under `*`, `+`, `?`, `{2}`, `{1}` and `*?`, before `C`, `[a-z]`, `\w`, `$` or nothing, with and without `icase`, over `SMC` with search and full match | 26,878 | 8,819 | 0 | 794,596 |
| `SMWORD` | ten bodies with and without an inner repeat or alternation (`\w\w+`, `\w{2,}`, `\W\W??\|a`, `\w+?`, `a\|\W`, `\w\w`, `ab`, `\w`, `.\W?`, `(?:a\|b)b`) under seven quantifiers, in 16 placements of `\<`, `\>`, `\b`, `\B` or `$` (after, before and inside the repeated group, at the start of an alternative inside it, after a lookahead, in an atomic group and a lookahead, after `(?i)` or `(?i:a)`), over `SMW` with search, full match and tokens (-1 and 0) | 4,620 | 1,580 | 0 | 237,258 |
| `SMDEPTH` | around the recursion limit: 40 to 102 `$`, `\<` or `\>` before and after atoms, `x*`, `(?:a+)*` and `(?:a\|b)*`; alternations of 41 to 56 branches in repeated groups; 15 to 48 nested repeated groups with and without `$`; 40 to 90 `$` beside 4 to 12 alternations that end at `\z` or `(?-m)$`; lookbehinds with 10, 16, 1,026 and 1,030 alternations; over `SM` with search and tokens (-1) | 194 | 57 | 18 | 4,560 |
| `SMBR` | backreferences in the shapes above: to a group before a case-switched repeated group, from inside one, to a later or unset group, before `\<` and `\>`, after 60 `$`, and to a group that can match empty, with and without `icase`, over `SM` with search, full match and tokens (-1) | 400 | 184 | 0 | 12,312 |
| `HEX` | `\x` followed by every pair of `0 1 4 7 a A f F g x X + - space tab newline } { ,` and an optional `a` or `}`, in and outside a bracket expression, and `\x{...}` with up to three of `0 1 4 7 a A f x X + - space tab }`, also followed by `1`, over `HEX` with search | 10,975 | 942 | 8,744 | 27,363 |

The 18 `SMDEPTH` patterns both refuse are Boost's `error_complexity` (16) and the
lookbehinds of 1,026 and 1,030 alternations in a row (2). Commit e6c581f gives 10,852, 1,003 and 12 different answers in
`SMCASE`, `SMWORD` and `SMBR` and 778 different compile outcomes in `SMDEPTH` and `HEX`
(above). Beyond the corpus, the review's own grids, rerun on this change, compare with 0
mismatches and 0 different compile outcomes: `F7b` (119,680 patterns, 29,568 refused,
4,866,048 cases), `F7c` (18,368, 7,396, 268,016), `F11` (12,240, 6,468, 329,004) and
`F13` (51,402, 1,739, 908,649), where e6c581f gave 109,656, 25,338, 34,625 and 0
different answers and 584 different compile outcomes in `F13`; its `OMNI` random
grammar, seeds 1 to 4 (240,000 expressions, 98,505 refused, 19,101,825 cases), gives 0
different answers, where e6c581f gave 9, 63 and 6 on seeds 2 to 4. The review's depth
probes (91 patterns) give 0 different compile outcomes.

The sixth round adds grids around the bracket-expression and `(?` parsers. `COLLIN`
has 137 inputs (every byte, `0xE9`, `0xFF`, `0x80` and short strings with `]`, `-`
and `[.a.]`), `PEXTIN` 21 of `a`, `b`, `c`, `n` and spaces in case pairs, and `CLSIN`
133 (every ASCII byte, `0xE9` alone and as UTF-8, `0xFF`, `0x80`, `ab` and the empty
string):

| Family | Templates | Patterns | Refused | Refused by both | Compared cases |
| --- | --- | ---: | ---: | ---: | ---: |
| `COLL` | 25 set items (`[.a.]`, `[.-.]`, `[.].]`, `[...]`, `[.NUL.]`, `[.tab.]`, `[.ab.]`, truncations such as `[.`, `[.a`, `[.a.`, and the literal-`[` neighbours `[x]`, `[=a=]`, `[:alpha:]`, `[`, `[]`) as the end of a range from 12 starts (`A`, `a`, `!`, `\n`, `0`, `\x41`, `[`, `-`, `]`, `\\`, `z`, `Z`) in 11 templates (plain, negated, beside another item, before a second dash, as the start of the range, alone, twice), with and without `icase`, over `COLLIN` with search and tokens (-1) | 4,364 | 1,025 | 2,846 | 210,020 |
| `PEXT` | 37 tails after `(?P` (every option-group spelling, `>name`, `<n>a`, `=n`, `&n`, `1`, `P`, `#c`, `'n'a`, `*FAIL`, `R`, `+1`, `-1`, `?i`, truncations) in four prefixes (alone, after a named group, after a literal, inside a group) and six placements (plain, before `c`, in a non-capturing group, alternated, after `(?i)`, repeated), over `PEXTIN` with search, full match and tokens (-1) | 865 | 48 | 635 | 11,466 |
| `CLSNAME` | every name `get_default_class_id` knows and 12 it does not (`ascii`, `any`, `assigned`, `b`, `foo`, `Alnum2`, `wo rd`, `<`, `>`, empty, `^`, `al:pha`) in four case spellings, in 8 templates (plain, negated with `^` inside and outside, beside a literal or another class, as a range endpoint), with and without `icase`, over `CLSIN` with search and tokens (-1) | 1,570 | 52 | 681 | 224,316 |

`COLL`'s 2,846 and `CLSNAME`'s 681 patterns that both refuse are Boost's
`error_collate`, `error_range`, `error_brack` and `error_ctype`; `PEXT`'s 635 are the
Python spellings, the unknown option letters and the truncations. Beyond the corpus,
a hand-written grid of the same three families over its own inputs (5,958 patterns,
386,073 cases) also compares with 0 different answers and 0 different compile
outcomes, where 9597f21 gave 81,913 different answers and 2,067 different compile
outcomes.

The sibling hunt of the fifth round also looked beyond that review's examples. Two more ways
past Boost's recursion limit turned up (an alternation of 51 branches in
`(?:(...)x?)+`, and alternations that end at a buffer end, which also make Boost's
construction exponential, see [Construction cost](#construction-cost)), and so did
Boost's lookbehind alternation stack and a translation nesting past the engine's
parser limit (`(?:` × 48, `a+`, `)+` × 48, which the engine refused after the
optimizer's shielding): all four are refused now. Seeded explorations outside the
corpus: 7,700 random expressions of anchors, alternations, lookaheads and nested
repeated groups whose recursion bound is between 60 and 100 all compile in Boost; and
300,000 random nested expressions (180,000 of plain, flag and non-capturing groups
under counted, lazy and bounded repeats around anchors, sets and case-insensitive
letters, and 120,000 with positive lookaheads, repeated lookarounds and zero repeats)
never reached the engine's parser limit, with `refusals.py` agreeing on every
refusal; about 600 of the compiled ones sat at the highest bound the translator allows.
The backreference grids of the third round (`BROPEN`, `BRFWD`, `BRUNSET`, `BRREP`,
`BRLOOK`, `BRDUP`, `BRZERO`, `BRNULL`, `BRICASE`, `BRLEAD`) and the `ATOMBR` grid of the
fourth, which cover backreferences to open, later and unset groups, into and out of
lookarounds and atomic groups, to repeated and duplicate-named groups, under `icase`
over bytes above `0x7F`, beside lazy leading repeats and to groups that can match
empty, were rerun with this change; their refusals and answers are unchanged.

Over the full corpus the facade refuses 39,452 patterns: 17,091 for a repeat of more
than one iteration of a group that can match the empty string, 7,973 for a repeat of a
group holding a repeat or alternation under other case sensitivity, 2,693 for a leading
lazy repeat, 2,451 for a backreference inside the group it refers to, 1,321 for `\<` or
`\>` beside such a repeated group, 1,036 for a collating element, 761 for a `\x` escape
with stream syntax, 92 for `\<` after a case switch, 49 for `[[:unicode:]]`, 48 beyond
the start-map recursion bound, 12 for a class name that is neither a class nor all
letters, 7 for a recursion (`(?P>name)` among them), 4 for an equivalence class, 4 for
a lookbehind with more than 1,024 alternations and 2 for a translation nested past the
engine's parser limit among them.

Adversarial inputs in every family: the empty string, `\n`, `\r`, `\r\n`, `\f`,
`\v`, tab, space, NUL, `\n\r`, `\r\r\n`, `\f\n`, `\r\f`, `é`, lone `0xFF`, `0xC3`,
`0x85` and `0xA0`, U+0663, NBSP, NEL and U+2028, plus structured inputs with
separators embedded (`scan=12\r\n`, `foo\f500_12`, `[1]\r\n[2]`).

Refusals over the committed fixture, the fuzz family included, by category
(`refusals.py`, equal to the facade's; the test asserts them as
`EXPECTED_REFUSALS`). Each construct is a row of the
[refused-construct table](#refused-constructs) above, which gives its Boost-side
reason:

| Construct | Patterns |
| --- | ---: |
| a repeat inside an atomic group or a negative lookaround | 93 |
| a capturing group inside a lookaround or atomic group | 84 |
| a repeat of more than one iteration of a group that can match the empty string | 73 |
| a backreference inside the group it refers to | 29 |
| a lazy repeat with a finite maximum of a one-byte atom that starts the expression | 27 |
| a repeat of a group that can match the empty string and captures | 26 |
| a repeat of a group holding a repeat or alternation under other case sensitivity | 20 |
| `\<` or `\>` in an expression with a repeated group holding a repeat or alternation | 17 |
| a collating element | 14 |
| `\<` in an expression that switches case sensitivity | 11 |
| a repeat of a group that only asserts a position | 11 |
| a hexadecimal escape whose digits start with a sign, white space or `0x` | 10 |
| `\Z` | 7 |
| a counted repeat of a backreference that can match the empty string | 7 |
| an escape letter without a translated meaning | 7 |
| the `[:unicode:]` class | 7 |
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
| more line anchors, word-boundary assertions, alternations and repeated groups than Boost's start-map recursion limit allows | 3 |
| the `x` (extended) modifier | 3 |
| `\K` | 2 |
| a backreference with more than one digit | 2 |
| a backtracking control verb | 2 |
| a character class name that is not a POSIX class | 2 |
| a conditional expression | 2 |
| a lookbehind wider than `MAX_LOOKBEHIND_WIDTH` bytes | 2 |
| a recursive sub-expression | 2 |
| an octal escape | 2 |
| groups and repeats nested deeper than the engine's parser allows once translated | 2 |
| `\G` | 1 |
| a branch-reset group | 1 |
| a group name outside `[A-Za-z0-9_]` | 1 |
| a repeat bound above `MAX_REPEAT` | 1 |
| an equivalence class | 1 |

266 of the 512 are probes outside `FUZZ`, listed one by one in
`EXPECTED_UNSUPPORTED`; 246 are grammar-generated. The test asserts the list, these
counts and the number of compared cases. The sixth round adds two categories and
raises three: `[:unicode:]` 7, a class name that is neither a class nor all letters
2, a collating element 1 to 14, a recursion 1 to 2 (`(?P>name)`) and the `x`
modifier 1 to 3 (`(?Pim-sx:a|b)`, `(?Px)a b`). All 25 are new `ADV` probes: no
pattern of an earlier family, and no `FUZZ` pattern, changed side, and the
grammar-generated refusals stay 246. The fifth round added six categories, which
refuse 36 of its probes and 26 `FUZZ` patterns that compiled before (and compared
equal on the fixture's inputs), and one `FUZZ` pattern that the capturing-group rule
refused before now stops at an earlier rule. The backreference refusal is from the
third round: its 29 patterns are that round's 16 probes, `(a\1)` from `SYN`, 5
grammar-generated patterns and 7 patterns another rule refused before it.

**Random expressions.** `../oracle/boost-regex/fuzz.py <mode> <seed> <count>` writes
expressions from one of four grammars, each run with search, full match and the
token iterator (-1) over 34 short inputs (49 for `prefix`), and the harness that includes
`src/concept/boost_regex.rs` verbatim compared them with Boost and with
`refusals.py`, which agreed with the facade on every refusal and category. A search
that one side stops with its limit (Boost's `error_complexity`, the facade's
budget) where the other answers is counted apart from a different answer:

| Mode | Seeds and expressions | Refused | Compared cases | Different answers | Limit disagreements | cfbaf1a: different answers | cf594ef: different answers, of them panics |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `backref`: groups, lookarounds, atomic and flag groups, anchors, `\1`-`\3` as a quarter of the atoms | 1 (40,000), 2 to 5 (100,000 each) | 122,084 | 4,300,380 | 0 | 0 (3, all the facade's budget, before this round's refusals) | 0 | 19,954, 1,880 |
| `repeat`: the same with lazy and `{2,}` quantifiers and 2% backreferences | 11 to 13 (100,000 each) | 125,735 | 12,308,814 | 0 | 14 (11 Boost's, 6 the facade's) | 0 | 5,585, 1,232 |
| `easy`: no lookaround, anchor, atomic group or backreference, more alternation | 21 to 23 (100,000 each) | 135,190 | 20,766,060 | 0 | 71, all Boost's | 0 | 464, 0 |
| `prefix` (fourth round): alternations of 2 to 4 branches sharing a prefix that can match in several ways, nested up to three deep in capturing, non-capturing, atomic, lookahead, flag and named groups, in 18 contexts (optional, repeated, alternated, after a lookahead, in a lookahead inside a lookbehind, after an atomic lookbehind, beside backreferences) | 31 to 33 (100,000 each) | 152,303 | 17,425,821 | 0 | 0 | 9,163 | 37,210, 261 |

The table is the fifth round's rerun of the same seeds, which that round's refusals
changed: they refuse 402, 2,115, 11,677 and 3,358 more expressions of the four
grammars, and with them the 3 cases the facade's budget stopped in `backref` and 94 of
the cases Boost's `error_complexity` stopped in `easy`. The sixth round leaves the
table as it stands: none of the four grammars can produce a `(?P` group, a collating
element or a `[[:name:]]`, so none of its rules can fire there. Confirmed by a fresh
run of all four grammars at seed 20260915 (20,000 expressions each: 80,000
expressions, 32,412 refused, 3,638,166 compared cases, 0 different answers, 0
different compile outcomes, 6 limit disagreements in `easy`, all Boost's
`error_complexity`), whose refusals, compared cases and limit disagreements are
identical at 9597f21 and here, and by a search of the four generated corpora, which
contain no such construct. The facade's budget errors are
lookahead bodies the engine backtracks into (for example
`(?i:(?i:(ba+?|.A^)*?[ab](?=\w??.{0,}?))*?(?<=b)(?<n>[ab]))ca|(?=(?:)c|)\b` on
`aabbaabb`), described under [Work bounds](#work-bounds). The different answers at
cf594ef were the backreferences inside their groups and the engine rewrites of the
third round; those at cfbaf1a the alternations inside atomic groups of the fourth
round (every one of the 974 distinct patterns behind them contains `(?>`). A scratch
build with `MAX_AUTOMATON_ATOMS` set to 0 also gives 0 different answers on the
`prefix` seeds. The fuzzing also found
that `refusals.py` did not mirror one detail of the translator: a repeated
lookaround that the facade spells as an optional group (`(?:(?<!\w)){0,}?`) no
longer counts as a lone lookaround of the enclosing group; it does now, and no
corpus pattern was affected.

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
`backreferences_inside_their_group_are_refused` and
`engine_rewrites_keep_boost_answers` assert the third review's probes: every
backreference inside its group is refused with its category within a watchdog, and
the neighbours give Boost's answers (`((a)|\2b)x` on `abax` matches `2..4` with both
groups at `2..3`, `(?:\1b|(a))+` on `aaab` `0..4` with group 1 at `1..2`, `(.)\1`
under `icase` matches `aA` and not `éÉ`); `\w{1,}b?\w{1,}` does not match `a`,
`(.{1,})+\1+` matches all of `abb` with group 1 at `1..2`, `(\w+?)*?(?<=b)` matched
against all of `ab` captures `1..2`, and `[ab]+b|[ab]+c` finds `0..3` in `abbc`.
`atomic_alternations_keep_boost_answers` asserts the fourth review's probes:
`(?>[ab]?b|[ab]?c)` finds `0..1` in `bc` (and `2..3` after a UTF-8 `é`),
`(?>(?:a|ab)c|(?:a|ab)b)` matches all of `abc`, and under `icase` all of `Abc` but
not `aAbc`, `(?=(?>(?:a?|ab)c|(?:a?|ab)b)c)` matches at 1 in `abc`,
`(?!(?>[ab]?b|[ab]?c)c)\w` finds `1..2` in `bc`, the lookbehinds
`(?<=(?=(?>[ab]?b|[ab]?c)c)..)` and `(?<=(?>ab|a[bc]))c` find `2..2` in `bcc` and
`2..3` in `abc`, nested, repeated and backreferenced forms give Boost's spans, and the
tokens -1 and 0 of `(?>[ab]?b|[ab]?c)` over `abcc` are Boost's six.
`start_map_shapes_are_refused` asserts the fifth review's start-map probes: each is
refused with its category within a watchdog (the case-switched repeats, `\<` after a
switch, `\<` and `\>` beside repeated groups with an inner repeat, and the recursion
bound, including four expressions Boost rejects with `error_complexity` and three it
compiles), and the neighbours give Boost's answers (`(?i:bc)+C` finds `1..4` in
`aBcC`, `(?:b.+(?i))*C` finds `1..2` in `bCc`, `(?:(?i)a|b)*C` `0..2` in `bCc`, `(?i)\<A`
under `icase` `1..2` in ` A`, `(?:\W\W??|a)+\B` `0..3` in `aaaa`, `(\w+?)+$` all of
`aaaa` with group 1 at `3..4`, `(?:(?:a|...)x?)+` with 49 branches all of `abab`, and
`$` × 60 before `(?:\b|\B)` × 12 and `\z` `2..2` in `A\n`).
`escapes_and_limits_follow_boost` asserts that the `\x` escapes Boost reads with stream
syntax are refused whether Boost accepts them (`\x{+41}`, `\x-0`) or not (`\x0x`), that
plain ones give Boost's answers and syntax errors (`[\x4-\x{41}]` finds `2..3` in
`ab ab`; `\x{4 }` and `\xg` are errors), that lookbehinds with 1,026 alternations in a
row or 1,030 in one are refused, and that 20 levels of `(?:...a{2,}){2,}` are refused
while 19 compile and do not match `aaaa`, as in Boost.
`bracket_elements_and_class_names_follow_boost` and
`python_style_group_openings_follow_boost` assert the sixth review's probes: a
collating element at either endpoint of a range is refused with its category,
whether Boost compiles it (`[A-[.a.]]`, `[[.a.]-z]`, `[!-[.].]]`) or rejects it
(`[A-[.ab.]]`, `[A-[.a.]-z]`); the literal `[` gives Boost's answers (`[A-[x]`
matches `A`, `[`, `Z` and `x` and not `\` or `a`, `[A-[=a=]]x` matches all of `=]x`);
`[[:unicode:]]` is refused in every spelling and `[[:<:]]`, `[[:>:]]` with the
class-name category, while `[[:FOO:]]` stays an error; 30 class-name spellings match
exactly the bytes Boost matches out of `0 9 a z A Z _ - space tab \n \v \r 0x7f NUL
0xE9` (so `[[:D:]]` and `[[:d:]]` are `[[:digit:]]`, `[[:h:]]` is space and tab,
`[[:v:]]` is `\n`, `\v` and `\r`); and 29 `(?P` cases give Boost's spans (`(?Pi)A`
matches `a`, `(?P:a|b)c` finds `1..3` in `abc`, `(?P)` matches empty at 0, `(?Pm)^a`
finds `2..3` in `b\na`, `a(?Pi)b` matches `aB` and not `Ab`), with `(?P>name)` and
the `x` modifier refused and the Python spellings errors.

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
- `backreferences_inside_their_group_are_refused`: the third review's
  backreferences inside their groups, among them `(?:(a|\1b)c)+`, on which the engine
  panicked, are refused within a watchdog.
- Module tests `refuses_start_maps_boost_builds_wrongly`, `bounds_start_map_recursion`,
  `bounds_translated_nesting` and `refuses_hex_escapes_read_as_stream_syntax`: the case
  sensitivity the translator gives Boost's repeat and alternation states (a later `|`
  takes the one at the previous `|`), the recursion bound exactly at 100 and one site
  beyond for anchors, looped anchors, alternations in repeated groups, alternations
  before a buffer end and sequential repeated groups, and, for five deeply nested
  expressions wrapped in more and more groups, that the first refusal is the nesting
  bound and that one group less compiles every engine program.
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
This round the same probe, with `\3`, `|`, `(?:`, `+`, `{0,}?`, `{2,}`, named groups,
`\k<n>`, `[^\s\S]`, `b?` and `a+` added to its pieces, constructed 200,000 patterns
(8,660 compiled) and searched 30,000 (1,309 compiled) over short inputs and 70 KB,
without a panic and with no operation slower than 0.5 s; the random-expression runs
above found no panic either.

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
  the lookahead body;
- a case-insensitive backreference (`matches_literal_casei`) checks its whole group
  after a failed exact comparison, before it compares without case, so an attempt
  that fails at the first byte still costs the group's length
  (`([a-c]*?)(?i)[^a](a\1)` over 4 MiB: 70 s, Boost 0.16 s);
- a backreference inside its own group, in an iteration after the first, slices
  the haystack from the group's new start to its previous end and panics when the
  start is the larger (`vm.rs`, `Insn::Backref`; `(?:(a|\1b)c)+` on `acbc`);
- the optimizer's `optimize_ambiguous_concat_repeats` turns `X+Y?X+` into
  `X+(?:Y{1}X+)?`, which also matches one `X`, and `optimize_nested_repeats` turns
  `(X+)+` into `(X+)` even when a backreference reads the group, and `(X)*` into
  `(X)?`, which captures differently when `X` is a lazy repeat.

In `regex-syntax` 0.8.11, which `fancy-regex` and `regex-automata` build their
automata from: `Hir::alternation` factors a common prefix out of branches that are
all concatenations (`lift_common_prefix`), which breaks leftmost-first priority
when the prefix can match in several ways (`[ab]+b|[ab]+c` on `abbc` gives `0..4`
in the `regex` crate as well, where Perl gives `0..3`). `fancy-regex` compiles the
body of an atomic group as if nothing followed it and hands it to such an
automaton when it needs no backtracking, so the factoring also changes which match
an atomic group keeps (`(?>[ab]?b|[ab]?c)` on `bc` gives `0..2`, where Perl gives
`0..1`).

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
  the PikeVM (the automaton ceilings above) and `utf8_empty` in bytes mode;
- the parser's group nesting limit (`MAX_RECURSION`, 64), which the translator's
  nesting bound mirrors (`bounds_translated_nesting` fails if a program nests deeper
  than the bound allows);
- the optimizer: the facade calls `fancy_regex::internal::optimize`,
  `fancy_regex::internal::FLAG_UNICODE` and `Expr::parse_tree_with_flags`, which are
  outside `fancy-regex`'s stable API, and relies on `RegexBuilder` parsing with only
  the Unicode flag and on no optimizer pass matching an alternation where it matches
  a repeat (`shields_repeats_the_optimizer_would_rewrite` fails if a shielded
  spelling is rewritten again);
- `regex-syntax`'s `Hir::alternation`: whether `lift_common_prefix` still needs
  every branch to be a concatenation (`guards_alternations_on_the_automaton` and
  `guards_alternations_inside_atomic_groups` check that the guard still changes the
  answer the factoring gives); and which sub-expressions `fancy-regex` hands to
  `regex-automata` (`Compiler::visit`, `compile_concat`, `compile_lookaround_inner`).
  In 0.19.2 they are whole expressions that need no backtracking, lookaround bodies,
  atomic-group bodies and their branches and trailing runs, and fixed-size runs
  elsewhere. The facade guards the alternations of the first and of the third
  outside lookbehind widths, and relies on lookaround bodies being asked only
  whether they match and on fixed-size runs ending at one position whichever branch
  matches; any other delegated sub-expression with an alternation would need a
  guard too.
