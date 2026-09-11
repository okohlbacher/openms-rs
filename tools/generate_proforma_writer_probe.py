#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
"""Assemble unchanged source AST/writer blocks for a standalone C++ text probe.

This is an extraction probe, not a full SDK build. Only unused DLL decoration
and the opaque resolved-modification pointer need declarations; no scientific
implementation is replaced. Source block hashes are emitted beside the probe.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    paths = ['src/openms/include/OpenMS/CHEMISTRY/ProForma.h',
             'src/openms/source/CHEMISTRY/ProForma.cpp']
    header, source = [(args.source / p).read_text() for p in paths]
    begin = header.index('    enum class ConversionPolicy\n')
    last = header.index('    struct OPENMS_DLLAPI CrossLinkGroup\n', begin)
    ast = header[begin:header.index('    };', last) + len('    };')]
    marker = source.index('// Internal ProFormaWriter class')
    begin = source.index('namespace detail\n', marker)
    writer = source[begin:source.index('} // namespace detail', begin) + len('} // namespace detail')]
    harness = (ROOT / 'tests/data/proforma_writer_probe.cpp').read_text()
    prelude = '''#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>
#define OPENMS_DLLAPI
namespace OpenMS {
class ResidueModification;
class ProForma { public:
'''
    args.output.write_text(prelude + ast + '\n};\n' + writer + '\n}\n' + harness)
    digest = lambda b: hashlib.sha256(b).hexdigest()
    record = {
        'sources': [{'path': p, 'sha256': digest((args.source / p).read_bytes())} for p in paths],
        'ast_block_sha256': digest(ast.encode()), 'writer_block_sha256': digest(writer.encode()),
        'translation_unit_sha256': digest(args.output.read_bytes()),
        'method': 'Exact contiguous source AST and complete writer blocks; only standalone includes, DLL macro, opaque pointer forward declaration, namespace/class enclosure and test harness added.',
    }
    args.output.with_suffix('.json').write_text(json.dumps(record, indent=2) + '\n')


if __name__ == '__main__':
    main()
