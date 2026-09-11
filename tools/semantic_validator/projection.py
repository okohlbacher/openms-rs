#!/usr/bin/env python3
"""Verify pinned raw semantic-validator inputs and the literal fixture inventory.

This does not execute/import Rust or claim an executed C++ oracle. Diagnostics
are copied from the pinned SemanticValidator class test, with source line numbers.
"""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'tests/data'

def projection(manifest):
    result = []
    for row in manifest['raw_fixtures']:
        raw = (ROOT / row['native_path']).read_bytes()
        assert len(raw) == row['bytes'], row['native_path']
        assert hashlib.sha256(raw).hexdigest() == row['sha256'], row['native_path']
        item = {'path': row['native_path'], 'bytes': len(raw)}
        if row['native_path'].endswith('.xml'):
            root = ET.fromstring(raw)
            item['elements'] = sum(1 for _ in root.iter())
            item['cv_params'] = sum(1 for e in root.iter() if e.tag.split('}')[-1] == 'cvParam')
        else:
            item['term_stanzas'] = raw.count(b'[Term]')
        result.append(item)
    return {'fixtures': result, 'literal_errors': manifest['class_test_literals']['errors'],
            'literal_warnings': manifest['class_test_literals']['warnings']}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((DATA / 'semantic_validator_provenance.json').read_text())
    text = json.dumps(projection(manifest), indent=2, ensure_ascii=False) + '\n'
    output = DATA / 'semantic_validator/projection.json'
    if args.check:
        assert output.read_text() == text, 'stale source projection'
        print('Four raw fixtures and all nine source diagnostic literals verified.')
    else:
        output.write_text(text)

if __name__ == '__main__':
    main()
