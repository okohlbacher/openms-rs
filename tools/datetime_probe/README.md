# Bounded DateTime reference probe

This harness compiles the pinned DateTime.cpp scientific body and its DateTime.h/HashUtils.h headers **unmodified**. It does not build the SDK. The `include` directory supplies only:

- UInt = unsigned int and export/pretty-function macros;
- StringUtils::has/prefix with the same source find/substr behavior;
- unsigned integer-to-text for error construction through std::to_string;
- ParseError/InvalidValue exception tags, without source message formatting or hierarchy.

The probe fixes the C locale and reports operation status, signed component projections, all seven source rendering formats and fallback getters. It does not call current-time functions or hash values. Hex fields preserve control characters and NUL bytes. `generate_cases.py` extracts the 19 immediate set/get or exception class-test literals and supplies bounded independent adversarial inputs. No expected value comes from native Rust. `derive_calendar_corrections.py` uses independent Python calendar.monthrange month-stepping for the 35 explicitly identified macOS timegm departures.

With `SOURCE` set to the exact 82ce5b3 checkout and `WORK` to an existing temporary directory, from the Rust repository root:

```sh
clang++ -std=c++20 -O0 -g -fsanitize=undefined -fno-sanitize-recover=undefined \
  -I tools/datetime_probe/include -I "$SOURCE/src/openms/include" \
  tools/datetime_probe/main.cpp "$SOURCE/src/openms/source/DATASTRUCTURES/DateTime.cpp" \
  -o "$WORK/probe"
python3 tools/datetime_probe/generate_cases.py --source-root "$SOURCE" \
  --probe "$WORK/probe" --work-dir "$WORK" --output-dir "$WORK/fixtures"
python3 tools/datetime_probe/derive_calendar_corrections.py --data-dir "$WORK/fixtures"
clang++ -std=c++20 tools/datetime_probe/timegm_diagnostic.cpp -o "$WORK/timegm_diagnostic"
"$WORK/timegm_diagnostic"
```

The libc diagnostic separately records timegm return/errno for sampled years -1 through i32::MAX. The captured macOS source fixture is platform-specific for the 35 early-year arithmetic rows. A different libc may produce correct calendar results there; do not replace or claim byte equality with this original fixture without recording the new environment and separating its outputs. No source/vendored/contrib file is changed by the harness.
