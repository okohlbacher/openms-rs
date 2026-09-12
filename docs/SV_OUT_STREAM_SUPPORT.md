# Separated-value output stream

Source pin: `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. This increment ports
`FORMAT/SVOutStream.h` (192 lines) and `FORMAT/SVOutStream.cpp` (138 lines) to
`src/format/sv_out_stream.rs`. Seven TOPP tools depend on the header —
`MassCalculator`, `MetaProSIP`, `NucleicAcidSearchEngine`, `ProteinQuantifier`,
`RNAMassCalculator`, `SeedListGenerator` and `TextExporter` — and none of them
is ported yet; this is the writer they will all use.

The class is a `std::ostream` subclass whose `operator<<` overloads insert a
separator between items but not at the start of a line, quote strings, and
recognise `nl` and `std::endl` as line delimiters. Rust has no `operator<<`, so
the overload set becomes named methods; the separator bookkeeping — one bit,
`newline_` — is ported exactly, including the two places where the source
leaves it deliberately untouched.

## API mapping

Every public member of `SVOutStream.h`, plus the three `StringUtils.h` members
the class cannot be ported without.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `enum Newline { nl }` | `SVOutStream::newline()` | A sentinel value in C++ because `operator<<` needs an argument type; a method here. The enumerator name `nl` has no Rust counterpart. |
| `SVOutStream(const std::string& file_out, sep, replacement, quoting)` | `SVOutStream::create_with_options(path, separator, replacement, quoting)` | Creates and truncates. `Error::Io` replaces `Exception::FileNotWritable`. |
| *(the same constructor with all three defaults)* | `SVOutStream::create(path)` | Tab separator, `_` replacement, `QuotingMethod::Double`. |
| `SVOutStream(std::ostream& out, sep, replacement, quoting)` | `SVOutStream::with_options(writer, separator, replacement, quoting)` | Takes any `W: Write` by value rather than borrowing a stream. Returns `Result`: an empty separator or a line break in either string is refused, which the source does not check. |
| *(the same constructor with all three defaults)* | `SVOutStream::new(writer)` | |
| `~SVOutStream()` | `SVOutStream::finish() -> Result<W>`, plus the implicit drop | The destructor closes an owned `ofstream`. Rust drops the writer either way; `finish` is the flush whose failure can be reported. |
| `SVOutStream& operator<<(std::string str)` | `SVOutStream::write_field(&str)` | Separator, then quote or substitute. Rejects `\n` (`Error::InvalidValue` for `Exception::IllegalArgument`) and additionally `\r`. |
| `SVOutStream& operator<<(const std::string& str)` | `SVOutStream::write_field(&str)` | The same method; Rust needs no by-value/by-reference pair. The header declares only the by-value form, but the class test has a section for the reference form. |
| `SVOutStream& operator<<(const char* c_str)` | `SVOutStream::write_field(&str)` | The source implements it as `operator<<(std::string(c_str))`. |
| `SVOutStream& operator<<(const char c)` | `SVOutStream::write_char(char)` | Via `StringUtils::toStr(c)` upstream, so it is quoted like any field. A Rust `char` is a scalar value, not a byte. |
| `SVOutStream& operator<<(std::ostream& (*fp)(std::ostream&))` | `SVOutStream::end_line()` | The manipulator overload exists to catch `std::endl`; `end_line` is newline plus flush. Forwarding an arbitrary manipulator has no counterpart — see *Native differences*. |
| `SVOutStream& operator<<(enum Newline)` | `SVOutStream::newline()` | Newline without flush. |
| `template<T arithmetic> SVOutStream& operator<<(const T&)` | `SVOutStream::write_number<T: SvNumber>(T)` | Converts through the ported `StringUtils::toStr`, never quoted. |
| `template<T non-arithmetic> SVOutStream& operator<<(const T&)` | `SVOutStream::write_display<T: Display>(T)` | Separator then the value, no quoting. Rust's `Display` carries no stream precision — see *Native differences*. |
| `SVOutStream& write(const std::string& str)` | `SVOutStream::write_raw(&str)` | Verbatim, no separator, and the line state is deliberately not changed. |
| `bool modifyStrings(bool modify)` | `SVOutStream::set_modify_strings(bool) -> bool` | Returns the previous state. |
| `template<NumericT> SVOutStream& writeValueOrNan(NumericT)` | `SVOutStream::write_value_or_nan<T: SvNumber>(T)` | Finite values go to `write_number`; the rest write `nan_`/`inf_`/`-inf_` unmodified and restore the modification state. |
| `std::ostream` base class | not ported: a base-class cast is the documented way to bypass the separator bookkeeping, and `write_raw` is the controlled equivalent | |
| `protected std::ofstream* ofs_` | the `W` type parameter | No raw pointer, no manual delete. `SVOutStream<BufWriter<File>>` is the filename case. |
| `protected std::string sep_` | `SVOutStream::separator()` | Read-only; the source has no accessor. |
| `protected std::string replacement_` | `SVOutStream::replacement()` | Read-only. |
| `protected std::string nan_` | `SVOutStream::nan_text()`, default `SVOutStream::DEFAULT_NAN` | Read-only; `"nan"`, lower case, and unrelated to the `"NaN"` that `appendNumeric` produces. |
| `protected std::string inf_` | `SVOutStream::infinity_text()`, default `SVOutStream::DEFAULT_INFINITY` | Read-only; `"inf"`, with `-` prefixed for the negative case. |
| `protected OpenMS::QuotingMethod quoting_` | `SVOutStream::quoting()` | Read-only. |
| `protected bool modify_strings_` | `SVOutStream::modify_strings()` | Readable as well as settable. |
| `protected bool newline_` | `SVOutStream::at_line_start()` | Exposed because a caller mixing `write_raw` with field writers has to be able to see it. |
| `protected std::stringstream ss_` | not ported: it exists only to test a function pointer against `std::endl`, which has no Rust analogue | Its unflushed state is a latent source defect — see *Native differences*. |
| `OpenMS::QuotingMethod` (`StringUtils.h:39`) | `sv_out_stream::QuotingMethod` | `NONE`/`ESCAPE`/`DOUBLE` become `None`/`Escape`/`Double`. |
| `StringUtils::quote(s, q, method)` (`StringUtils.h:553`) | `sv_out_stream::quote(&str, QuotingMethod)` | Only `q == '"'` is offered, which is the only value this class uses. |
| `StringUtils::substitute(s, from, to)` (`StringUtils.h:484`) | `str::replace` | Not re-exported. An empty pattern returns the input unchanged upstream; here an empty separator is refused at construction instead. |
| `StringUtils::toStr(float/double)` → `NumericFormatting::appendNumeric` | `sv_out_stream::source_float_text(f64, usize)`, `sv_out_stream::source_f32_text(f32, usize)` and `SvNumber` | Full precision only. `appendNumeric` is a template, so the `float` instantiation compares `abs_val` against `T(1e-2)`/`T(1e4)` and calls `std::to_chars` at `float` width; that is `source_f32_text`, and an `f32` is never promoted to `f64` before formatting. |
| `StringUtils::appendToStrLowP` (3 fractional digits) | not ported: `SVOutStream` never calls the low-precision variant | |

Native additions with no source counterpart: `SVLimits` and its three
constants, `SVOutStream::with_limits`, `limits`, `row_fields`, `rows`, `flush`,
`DEFAULT_SEPARATOR`, `DEFAULT_REPLACEMENT`, `NumberClass`, `SvNumber`,
`F64_FIXED_DIGITS`, `F32_FIXED_DIGITS`, `SCIENTIFIC_LOWER`, `SCIENTIFIC_UPPER`,
`MAX_FIXED_DIGITS`.

## Preserved source conventions

- **The separator state machine.** A field writer emits the separator only when
  the stream is not at the start of a line, and otherwise clears the flag
  (`SVOutStream.cpp:67`). `newline_` starts true, so no line ever begins with a
  separator.
- **Both line delimiters, and the flush that distinguishes them.** `nl` writes
  `"\n"` and sets the flag; `std::endl` does the same and flushes. The header
  recommends `nl` "for improved performance", which is exactly the missing
  flush.
- **A newline in a field is an error.** The class documentation says a literal
  `"\n"` "won't be accepted"; `SVOutStream.cpp:62` throws
  `Exception::IllegalArgument` for it, and `write_field` returns
  `Error::InvalidValue`.
- **Branch order inside a field** (`SVOutStream.cpp:76`): modification off
  writes verbatim; otherwise a quoting method other than `NONE` quotes; only
  `NONE` substitutes the separator. `QuotingMethod::None` therefore never adds
  quote characters through this class, even though `StringUtils::quote` would.
- **`quote`'s substitution order.** `ESCAPE` escapes backslashes *before* quote
  characters (`StringUtils.h:557-558`). Swapping them would double-escape.
- **`write` does not reset the line state** (`SVOutStream.cpp:125`). After a raw
  comment ending in `\n` written mid-line, the next field still receives a
  leading separator at column zero. That is why the header says "use only on a
  line of its own!", and it is reproduced, with a test.
- **`writeValueOrNan` writes its own texts unmodified**, switching modification
  off around them and restoring it afterwards, so `nan`/`inf`/`-inf` are never
  quoted even under `QuotingMethod::Double`.
- **Numeric text comes from `toStr`, not from stream formatting.** The
  arithmetic overload's independence from locale and stream precision is the
  reason `SvNumber` exists rather than a blanket `Display` bound.
- **The whole `appendNumeric` spelling**: `NaN` upper case, `inf`/`-inf` lower
  case, scientific notation for a nonzero magnitude outside `[1e-2, 1e4)`,
  fixed notation with 15 (f64) or 6 (f32) digits after the point, trailing
  fractional zeros trimmed to one surviving digit, the exponent's `+` dropped
  but its two-digit zero padding kept, and a mantissa without a point given
  `.0`. `1e4` is written `1.0e04`; `5.0` never degrades to `5`.
- **Both template instantiations, at their own width.** The source's
  `appendNumeric<float>` compares `abs_val` against `T(1e-2)` and `T(1e4)` in
  float arithmetic and hands the `float` to `std::to_chars`, whose shortest
  round-trip is the shortest decimal that round-trips as a *float*. `f32` text
  therefore goes through `source_f32_text`, not through the `f64` pipeline:
  `1.23e-5f32` is `1.23e-05` (a promotion gives `1.2299999980314169e-05`),
  `12345.6f32` is `1.23456e04`, and `0.01f32` — equal to `float(1e-2)`, so not
  below it — takes the fixed branch and prints `0.01`. Test:
  `f32_text_is_formatted_at_f32_width_not_through_f64`.

## Native differences

- **No `std::ostream` base.** The source's `static_cast<std::ostream&>(*this)`
  is how it writes past its own bookkeeping, and a caller can do the same. This
  type implements no `Write`, so the only way past the separator logic is
  `write_raw`, the source's own escape hatch.
- **The manipulator overload is one method, not a forwarding path.**
  `SVOutStream.cpp:102` applies the function pointer to the member
  `std::stringstream ss_`, compares the result with `"\n"` and forwards the
  manipulator to the real stream. Only `std::endl` matters, and `end_line` is
  that case. The forwarding of other manipulators is not ported because Rust has
  none — and it carries a latent defect worth recording: `ss_` is cleared only
  inside the `== "\n"` branch, so one manipulator that writes anything else
  (`std::ends`, or any user-defined one) poisons `ss_` permanently and
  `std::endl` is never recognised again, after which every line silently starts
  with a separator. Filed as a candidate, not confirmed by execution.
- **A carriage return is rejected too.** The source tests only for `\n`, so a
  lone `\r` reaches the file and a reader treating it as a line break sees a row
  the writer counted as one.
- **The separator and replacement are validated.** An empty separator makes
  `StringUtils::substitute` a no-op and emits nothing between fields, producing
  a file that cannot be read back; a line break in either string desynchronises
  the line state. Both are `Error::InvalidValue` at construction.
- **`write_display` has no stream precision.** The constructor's
  `precision(numeric_limits<double>::digits10)` affected only the generic
  non-arithmetic overload. Rust's `Display` has no such setting, so a type
  rendering a float through `write_display` is spelled by Rust;
  `write_number` is the source-equivalent route for numbers. `write_display`
  also rejects a newline, which the source's generic overload does not check at
  all.
- **`char` is a Unicode scalar value**, so `write_char` can emit a multi-byte
  field where the source emits one byte.
- **Errors instead of exceptions**, and `finish()` instead of a destructor,
  because a Rust drop cannot report a failed flush.
- **The source parallelises nothing here**; neither does this port. No gap.

## Checked boundaries and evidence

`SVLimits` bounds one rendered field or raw chunk (16 MiB), the fields on one
line (2^20) and the lines in one file (2^30). Worst-case quoting — double the
input plus two quote characters plus the separator — is charged *before* the
rendering allocation, so exceeding a ceiling leaves the output byte-for-byte
unchanged. That worst case bounds both `quote` paths, each of which at most
doubles the field, but it does not bound the `QuotingMethod::None`
substitution, whose growth factor is `replacement.len() / separator.len()`: a
caller-chosen replacement longer than the separator can grow a field without
limit. The substituted length is therefore computed exactly — by counting
separators, which allocates nothing — and charged before `str::replace` runs,
so the ceiling is consulted before the allocation there too. Test:
`a_replacement_longer_than_the_separator_is_charged_before_the_substitution`. A field and its separator are built into one buffer and written with
a single `write_all`, so a rejected or failed field emits no separator and does
not advance the line state; the tests assert that on both the rejection and the
I/O-failure path.

No string this module did not construct is ever byte-sliced. Quoting and
separator substitution go through `str::replace` and `char`-level pushes, and
the numeric formatter only inspects text produced by `format!`. A test writes a
multi-byte separator, a multi-byte replacement and CJK and emoji fields to pin
that.

Evidence is tier 3: every expected string in `tests/sv_out_stream.rs` is a
literal transcribed from the 12 sections of `SVOutStream_test.cpp`, which is
built upstream (`executables.cmake:284`). All 12 sections are ported, including
the three marked `NOT_TESTABLE` (`SVOutStream_test.cpp:100`, `:106`, `:112`). The upstream `-1.23e45` assertion accepts three
platform spellings and the pinned `NumericFormatting` produces the third,
`-1.23e45`, which is the one reproduced here. No C++ was built or executed and
no C++ output was retained, so this is not a tier-1 differential. The numeric
formatter beyond those literals, the resource ceilings and the rejection paths
are independently derived. Source hashes and line-level anchors are in
[the provenance record](../tests/data/sv_out_stream_provenance.json).

`StringUtils` itself remains unported; `QuotingMethod`, `quote` and the numeric
conversion live here because this class cannot be ported without them, and a
`DATASTRUCTURES/StringUtils` work package should take them over.
