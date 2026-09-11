#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# Source-expression oracle; never reads native Rust output.
import base64
import pathlib
import struct

root = pathlib.Path(__file__).resolve().parents[2]
values = [100.123456789, 200.25, -0.0]
intensities = [struct.unpack('<f', struct.pack('<f', x))[0] for x in [1.25, 2.5, -0.0]]
rows = ['mz32\tintensity32\tmass_np\tcoordinate_width\tintensity_width\tcoordinate_base64\tintensity_base64']
for mz32 in [False, True]:
    for intensity32 in [False, True]:
        for mass_np in [False, True]:
            # MzMLHandler.cpp5620-5621: source f32 preparation only without
            # any mass/time Numpress request, including ordinary fallback.
            xfmt = 'f' if mz32 and not mass_np else 'd'
            yfmt = 'f' if intensity32 and not mass_np else 'd'
            encode = lambda fmt, xs: base64.b64encode(b''.join(struct.pack('<' + fmt, x) for x in xs)).decode()
            rows.append('\t'.join(map(str, [int(mz32), int(intensity32), int(mass_np), struct.calcsize(xfmt), struct.calcsize(yfmt), encode(xfmt, values), encode(yfmt, intensities)])))
(root/'tests/data/mzml_writing/precision.tsv').write_text('\n'.join(rows)+'\n')
