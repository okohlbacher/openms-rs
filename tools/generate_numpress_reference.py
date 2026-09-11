#!/usr/bin/env python3
"""Project literal source fixtures and run an explicitly supplied C++ raw probe.

Usage: python3 tools/generate_numpress_reference.py SOURCE_ROOT PROBE_BINARY
Neither C++ source nor the implementation is modified. Requires only Python stdlib.
"""
import base64
import pathlib
import random
import struct
import subprocess
import sys

stage = pathlib.Path(__file__).resolve().parent.parent
source = pathlib.Path(sys.argv[1])
probe = pathlib.Path(sys.argv[2])
text = (source / 'src/tests/class_tests/openms/source/MSNumpressCoder_test.cpp').read_text()
literals = [('linear', 'QWR64UAAAADo//8/0P//f1kSgA=='), ('pic', 'ZGaMXCFQkQ=='), ('slof', 'QMVagAAAAAAZxX3ivPP8/w==')]
rows = ['codec\tsource_base64\traw_hex\tsource_line']
for mode, literal in literals:
    line = text[:text.index('TEST_EQUAL(out, "'+literal+'")')].count('\n') + 1
    rows.append(f'{mode}\t{literal}\t{base64.b64decode(literal, validate=True).hex()}\t{line}')
(stage / 'tests/data/numpress_source_bytes.tsv').write_text('\n'.join(rows)+'\n')

def bits(x):
    return struct.pack('>d', x).hex()

rng = random.Random(311804)
cases = []
for mode in ['linear', 'pic', 'slof', 'safe']:
    for n in [0, 1, 2, 3, 4, 5, 9, 17, 32]:
        for iteration in range(8):
            if mode == 'linear':
                data = [1000.0 + i*0.25 + rng.randint(-10,10)/8192 for i in range(n)]
                fp = [1.0, 4096.0, 1e6, 1048576.0][iteration % 4]
            elif mode == 'pic':
                data = [float(rng.randrange(0, 2_000_000_000)) + rng.choice([0.0,0.25,0.5]) for _ in range(n)]
                fp = 0.0
            elif mode == 'slof':
                data = [rng.randrange(0, 10_000_000)/128 for _ in range(n)]
                fp = [1.0,100.0,1000.0,5000.0][iteration%4]
            else:
                data = [rng.randint(-100000,100000)/256 for _ in range(n)]
                fp = 0.0
            cases.append((mode, fp, data))
# Defined source corners distinct from the ordinary positive series.
cases += [
    ('linear',1.0,[-1.0,-2.0,-3.0,-4.0]),
    ('linear',-1.0,[-100.,-200.,-300.,-400.]),
    ('linear',1.0,[4294967296.,4294967297.,4294967298.]),
    ('linear',1.0,[2147483648.,2147483648.,0.]),
    ('pic',0.0,[-0.5,-0.25,0.,0.499,0.5,2147483646.5]),
    ('slof',-1.0,[1.,1.5,2.]),
    ('safe',0.0,[0.1,1.0,1.0,0.1]),
]
commands = [f'{mode} {bits(fp)} {len(data)} '+ ' '.join(map(bits,data)) for mode,fp,data in cases]
run = subprocess.run([str(probe)], input='\n'.join(commands)+'\n', text=True, capture_output=True, check=True)
outputs = run.stdout.splitlines()
assert len(outputs) == len(cases), run.stderr
rows = ['id\tcodec\tfixed_point_bits\tinput_bits\tencoded_hex\tdecoded_bits\toptimal_linear_bits\taccuracy_0_001_bits\toptimal_slof_bits']
for i, ((mode,fp,data), output) in enumerate(zip(cases,outputs)):
    assert not output.startswith('ERROR:'), (i,mode,data,output)
    rows.append(f'{i}\t{mode}\t{bits(fp)}\t'+','.join(map(bits,data))+'\t'+output)
(stage/'tests/data/numpress_cpp_differential.tsv').write_text('\n'.join(rows)+'\n')
print(f'Projected three source literal strings and {len(cases)} executed C++ probe cases.')
