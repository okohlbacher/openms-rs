#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# Independent fixture bytes. No Rust output is consumed.
import base64
import pathlib
import struct
import zlib

root = pathlib.Path(__file__).resolve().parents[2]
values = [100., 200., 300.00005, 400.00010]
original = {
    'ordinary': base64.b64encode(struct.pack('<4d', *values)).decode(),
    # Unchanged MSNumpressCoder_test.cpp202/241/281 literals.
    'linear': 'QWR64UAAAADo//8/0P//f1kSgA==',
    'pic': 'ZGaMXCFQkQ==',
    'slof': 'QMVagAAAAAAZxX3ivPP8/w==',
}
terms = {'ordinary': ('MS:1000576', 'MS:1000574'),
         'linear': ('MS:1002312', 'MS:1002746'),
         'pic': ('MS:1002313', 'MS:1002747'),
         'slof': ('MS:1002314', 'MS:1002748')}
rows = ['mode\tzlib\taccession\tbase64']
for mode, text in original.items():
    for compressed in [False, True]:
        payload = base64.b64encode(zlib.compress(base64.b64decode(text))).decode() if compressed else text
        rows.append(f'{mode}\t{int(compressed)}\t{terms[mode][int(compressed)]}\t{payload}')
(root/'tests/data/mzml_normalization/payloads.tsv').write_text('\n'.join(rows)+'\n')
