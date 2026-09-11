#!/usr/bin/env python3
"""Project the pinned ProForma spectrum class tests without executing chemistry.

Every assertion and input remains source text. There are no generated mass,
annotation, or peak-count expectations. Pass the exact SDK checkout explicitly.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

REVISION = '82ce5b373c97f934ffd9b1ffd80215ca66473d0b'
SOURCE = 'src/tests/class_tests/openms/source/ProFormaParser_test.cpp'
SHA256 = '53b3b52b9490df9554f2373f229136a1fe1d66935fdcf2ac494f52ab33fc2696'
FIRST, LAST = 2998, 3135
TARGET = 'tests/data/proforma_spectra_source.json'


def project(raw):
    if hashlib.sha256(raw).hexdigest() != SHA256:
        raise SystemExit('Pinned ProForma class-test source hash mismatch')
    lines = raw.decode('utf-8').splitlines(keepends=True)
    section = None
    sections = []
    for line_number in range(FIRST, LAST + 1):
        line = lines[line_number - 1]
        title = re.fullmatch(r'START_SECTION\((.*)\)\n', line)
        if title:
            if section is not None:
                raise SystemExit('Unexpected nested source test section')
            section = {
                'title': title[1], 'line_start': line_number,
                'inputs': [], 'resolution_calls': [],
                'generation_calls': [], 'assertions': [],
            }
        if section is None:
            continue
        item = re.fullmatch(
            r'  (Peptidoform(?:Ion)?) (\w+) = ProForma::(parse(?:Ion)?)\(("[^"\n]*")\);\n',
            line,
        )
        if item:
            section['inputs'].append({
                'line': line_number, 'type': item[1], 'variable': item[2],
                'parse_method': item[3], 'literal': json.loads(item[4]),
                'source': line.strip(),
            })
        if 'ProForma::resolveModifications(' in line:
            section['resolution_calls'].append({
                'line': line_number, 'source': line.strip(),
            })
        generation = re.fullmatch(
            r'  MSSpectrum (\w+) = ProForma::generateSpectrum\((\w+), (\d+), (\d+), ("[^"]*"), (true|false), (true|false)\);\n',
            line,
        )
        if generation:
            section['generation_calls'].append({
                'line': line_number, 'result': generation[1],
                'input': generation[2],
                'min_charge': int(generation[3]), 'max_charge': int(generation[4]),
                'ion_types': json.loads(generation[5]),
                'add_losses': generation[6] == 'true',
                'add_metainfo': generation[7] == 'true',
                'source': line.strip(),
            })
        assertion = re.fullmatch(r'  TEST_EQUAL\((.*), (true|false)\)\n', line)
        if assertion:
            section['assertions'].append({
                'line': line_number, 'expression': assertion[1],
                'expected': assertion[2] == 'true', 'source': line.strip(),
            })
        if line == 'END_SECTION\n':
            section['line_end'] = line_number
            section['source_text'] = ''.join(lines[section['line_start'] - 1:line_number])
            sections.append(section)
            section = None
    if section is not None or len(sections) != 10:
        raise SystemExit('Expected exactly ten complete source sections')
    assertions = sum(len(s['assertions']) for s in sections)
    generations = sum(len(s['generation_calls']) for s in sections)
    if assertions != 18 or generations != 6:
        raise SystemExit('Expected eighteen source assertions and six generation calls')
    block = ''.join(lines[FIRST - 1:LAST]).encode('utf-8')
    return {
        'source_repository': 'https://github.com/okohlbacher/OpenMS4-core',
        'source_revision': REVISION,
        'source_file': {
            'path': SOURCE, 'sha256': SHA256,
            'line_start': FIRST, 'line_end': LAST,
            'projected_block_sha256': hashlib.sha256(block).hexdigest(),
        },
        'source_license': 'BSD-3-Clause',
        'source_copyright': raw.decode('utf-8').splitlines()[0].removeprefix('// '),
        'projection_tool': 'tools/generate_proforma_spectra_reference.py',
        'evidence_kind': 'Literal source test projection; no OpenMS or native spectrum execution.',
        'assertion_strength': (
            'Predicate/issue booleans and six coarse spectrum size assertions only: '
            'nonempty, at least ten peaks, or precursor-option count not smaller. '
            'No exact spectrum masses, annotations or total counts are supplied by these tests.'
        ),
        'resolution_note': (
            'Explicit source resolveModifications calls are retained per section. '
            'Unknown-modification and chimeric negative cases test predicates/issues, '
            'not generateSpectrum exceptions.'
        ),
        'section_count': len(sections), 'assertion_count': assertions,
        'generation_call_count': generations, 'sections': sections,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_root', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = project((args.source_root / SOURCE).read_bytes())
    encoded = (json.dumps(result, indent=2, ensure_ascii=False) + '\n').encode('utf-8')
    target = Path(__file__).resolve().parents[1] / TARGET
    if args.check:
        if target.read_bytes() != encoded:
            raise SystemExit(f'Fixture differs: {TARGET}')
    else:
        target.write_bytes(encoded)
    print(f'{TARGET}: 10 sections, 18 assertions, 6 generation calls; SHA-256 {hashlib.sha256(encoded).hexdigest()}')


if __name__ == '__main__':
    main()
