#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
"""Extract the exact C++ tokenizer, parser, AST and writer for text comparisons.

This is not a full SDK build. A capture-only ParseError adapter replaces the
source exception infrastructure; it records the unmodified parser's code,
position and message. No diagnostic formatting or backend behavior is probed.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    paths = ['src/openms/include/OpenMS/CHEMISTRY/ProForma.h',
             'src/openms/source/CHEMISTRY/ProForma.cpp',
             'src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h']
    header, source, strings = [(args.source / p).read_text() for p in paths]
    begin = header.index('    enum class ConversionPolicy\n')
    last = header.index('    struct OPENMS_DLLAPI CrossLinkGroup\n', begin)
    ast = header[begin:header.index('    };', last) + len('    };')]
    blocks = {'ast': ast}
    for name, marker in [('tokenizer', '// Internal ProFormaTokenizer class'),
                         ('writer', '// Internal ProFormaWriter class'),
                         ('parser', '// Internal ProFormaParserImpl class')]:
        begin = source.index('namespace detail\n', source.index(marker))
        end = source.index('} // namespace detail', begin) + len('} // namespace detail')
        blocks[name] = source[begin:end]
    begin = strings.index('    inline bool hasPrefix(const std::string& s, const std::string& prefix)')
    blocks['prefix_helper'] = strings[begin:strings.index('\n    }', begin) + len('\n    }')]
    prelude = '''#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>
#define OPENMS_DLLAPI
#define OPENMS_PRETTY_FUNCTION __func__
namespace OpenMS {
class ResidueModification;
class ProForma { public:
'''
    adapter = '''
// Probe-only exception transport: no source exception formatting/global state.
class ParseError : public std::runtime_error {
public:
  ErrorCode code;
  size_t position;
  ParseError(const char*, int, const char*, ErrorCode value, size_t pos,
             const std::string& input, const std::string& message)
    : std::runtime_error(message), code(value), position(std::min(pos, input.size())) {}
};
};
namespace StringUtils {
'''
    harness = (ROOT / 'tests/data/proforma_parser_probe.cpp').read_text()
    content = (prelude + ast + adapter + blocks['prefix_helper'] + '\n}\n'
               + blocks['tokenizer'] + '\n' + blocks['writer'] + '\n'
               + blocks['parser'] + '\n}\n' + harness)
    args.output.write_text(content)
    digest = lambda b: hashlib.sha256(b).hexdigest()
    record = {
        'sources': [{'path': p, 'sha256': digest((args.source / p).read_bytes())} for p in paths],
        'blocks': {name: digest(value.encode()) for name, value in blocks.items()},
        'translation_unit_sha256': digest(content.encode()),
        'method': 'Exact contiguous AST, tokenizer, parser, writer and StringUtils::hasPrefix blocks; standard includes/macros/enclosures plus capture-only ParseError adapter and test harness. No full SDK, exception formatting, JSON, resolution or scientific backend execution.',
    }
    args.output.with_suffix('.json').write_text(json.dumps(record, indent=2) + '\n')


if __name__ == '__main__':
    main()
