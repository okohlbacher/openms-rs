#!/usr/bin/env python3
"""Project source raw byte literals through independent Python zlib/base64.

Run from any directory. This does not execute the native writer or C++ wrapper.
The raw-byte TSV is generated/attributed by the separate raw Numpress stage.
"""
import base64
from pathlib import Path
import zlib

root = Path(__file__).resolve().parent.parent
rows = ['codec\tbase64_raw_source\tbase64_zlib_projection']
for row in (root/'tests/data/numpress_source_bytes.tsv').read_text().splitlines()[1:]:
    mode, literal, raw_hex, line = row.split('\t')
    rows.append('\t'.join([mode, literal, base64.b64encode(zlib.compress(bytes.fromhex(raw_hex))).decode()]))
(root/'tests/data/numpress_coder_transport.tsv').write_text('\n'.join(rows)+'\n')
print(f'Projected three raw source literals using Python zlib {zlib.ZLIB_RUNTIME_VERSION}.')
