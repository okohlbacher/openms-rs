#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. and contributors
# SPDX-License-Identifier: BSD-3-Clause
"""Execute only the pinned massDeltaText_ source helper, not the OpenMS SDK.

A tiny exception stand-in allows the exact helper to compile independently.
Expected outputs come from C++ std::to_chars(fixed), never native Rust output.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import tempfile

SOURCE_PATH = 'src/openms/source/CHEMISTRY/ProForma.cpp'
SOURCE_SHA = '9d625b063aa8c52c3be022d95343eb7d0f037fbad7ce045a3c94b0b9552a891c'
REVISION = '82ce5b373c97f934ffd9b1ffd80215ca66473d0b'

def cases():
    values = {0, 1, 2, 0x000fffffffffffff, 0x0010000000000000,
              0x7fefffffffffffff}
    for value in [0.00335, 1., 10., 12345.6789, 100000000000000016384., 2.**52, 2.**53]:
        bits = struct.unpack('>Q', struct.pack('>d', value))[0]
        values.update([bits-1, bits, bits+1])
    # Deterministic integer recurrence supplies arbitrary finite binary64 values.
    state = 0x1848ad11e1053127
    for _ in range(512):
        state = (state * 6364136223846793005 + 1442695040888963407) & ((1 << 64)-1)
        bits = state & ((1 << 63)-1)
        if (bits >> 52) != 0x7ff:
            values.add(bits)
    return sorted(values | {x | (1 << 63) for x in values})

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--cxx', default='c++')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    raw = (args.source_root / SOURCE_PATH).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == SOURCE_SHA
    source = raw.decode()
    start = source.index('  std::string massDeltaText_(double mass)')
    end = source.index('\n  /**', start)
    helper = source[start:end].rstrip()
    cpp = r'''
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#define OPENMS_PRETTY_FUNCTION __func__
namespace Exception { struct InvalidValue : std::runtime_error {
InvalidValue(const char*, int, const char*, const char* m, const std::string&) : std::runtime_error(m) {}
}; }
''' + helper + r'''
int main() {
  uint64_t bits;
  while (std::cin >> std::hex >> bits) {
    double value; std::memcpy(&value, &bits, sizeof(value));
    std::cout << std::hex << bits << '\t' << massDeltaText_(value) << '\n';
  }
}
'''
    inputs = cases()
    with tempfile.TemporaryDirectory(prefix='openms-conversion-number-') as directory:
        path = Path(directory)
        (path/'probe.cpp').write_text(cpp)
        subprocess.run([args.cxx, '-std=c++17', '-O2', str(path/'probe.cpp'), '-o', str(path/'probe')], check=True)
        version = subprocess.check_output([args.cxx, '--version'], text=True).splitlines()[0]
        output = subprocess.check_output([str(path/'probe')], input=''.join(f'{x:x}\n' for x in inputs), text=True)
    rows = output.splitlines()
    assert len(rows) == len(inputs)
    tsv = '# binary64_bits\tsource_massDeltaText_fixed\n' + ''.join(f'{int(row.split(chr(9))[0],16):016x}\t{row.split(chr(9))[1]}\n' for row in rows)
    root = Path(__file__).resolve().parents[1]
    fixture = root/'tests/data/proforma_conversion_mass_text.tsv'
    if args.check:
        assert fixture.read_text() == tsv, 'source formatter fixture changed'
    else:
        fixture.write_text(tsv)
    print(json.dumps({'rows': len(rows), 'source_revision': REVISION,
        'source_sha256': SOURCE_SHA, 'source_helper_sha256': hashlib.sha256(helper.encode()).hexdigest(),
        'fixture_sha256': hashlib.sha256(tsv.encode()).hexdigest(), 'compiler':version,
        'evidence':'Executed exact extracted source massDeltaText_ with std::to_chars(fixed); exception stand-in only; no full SDK.'}, indent=2))
if __name__ == '__main__':
    main()
