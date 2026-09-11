#!/usr/bin/env python3
"""Extract unchanged decimal arrays/name sets from the pinned XLMS class test.

No OpenMS execution or numerical recomputation. Pass the source checkout explicitly.
"""
import argparse
import hashlib
from pathlib import Path
import re

SOURCE = 'src/tests/class_tests/openms/source/TheoreticalSpectrumGeneratorXLMS_test.cpp'
SHA256 = '44b65c0954621f0d3d1f68b9efe716076172c194c0aead818214ba1a1ea31417'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_root', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    raw = (args.source_root / SOURCE).read_bytes()
    if hashlib.sha256(raw).hexdigest() != SHA256:
        raise SystemExit('Pinned source hash mismatch')
    source = raw.decode()
    arrays = re.findall(r'double result\[\] = \{([^}]+)\};', source)
    if [len(a.split(',')) for a in arrays] != [18, 17, 17]:
        raise SystemExit('Unexpected source array layout')
    masses = '# Literal source arrays; case, index, mz; source class-test tolerance 0.001 Da\n'
    masses += ''.join(f'{case}\t{i}\t{value.strip()}\n' for case, a in enumerate(arrays) for i, value in enumerate(a.split(',')))
    groups = re.findall(r'(?:  ion_names\.insert\("[^"\n]+"\);\n)+', source)
    if len(groups) != 6:
        raise SystemExit('Unexpected source name-set layout')
    names = '# Literal source allowed-name sets; case(0=linear,1=single,2=pair), alpha/beta, name\n'
    names += ''.join(f'{i // 2}\t{1 - i % 2}\t{name}\n' for i, group in enumerate(groups) for name in re.findall(r'insert\("([^"\n]+)"\)', group))
    root = Path(__file__).resolve().parents[1]
    for path, text in [('tests/data/xlms_source_masses.tsv', masses), ('tests/data/xlms_source_names.tsv', names)]:
        target = root / path
        if args.check:
            if target.read_bytes() != text.encode():
                raise SystemExit(f'Fixture differs: {path}')
        else:
            target.write_text(text)
        print(f'{path}: {len(text.splitlines()) - 1} rows, {hashlib.sha256(text.encode()).hexdigest()}')


if __name__ == '__main__':
    main()
