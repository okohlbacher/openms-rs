#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Guard the coverage ledger against false completion and dropped SDK scope."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import core_sdk_coverage as coverage


class CoverageTests(unittest.TestCase):
    def test_names_and_source_references_never_imply_parity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'docs').mkdir()
            (root / 'src').mkdir()
            (root / 'tests/data').mkdir(parents=True)
            files = [{'path': f'src/openms/include/OpenMS/FORMAT/{name}.h',
                      'sha256': 'a' * 64, 'domain': 'FORMAT',
                      'registration': 'registered_public_header'} for name in ['Existing', 'Missing', 'Referenced']]
            (root / 'docs/core-sdk-update.json').write_text(json.dumps({
                'identity': {'current_package_revision': 'b' * 40}, 'files': files}))
            (root / 'docs/core-sdk-reviewed-apis.json').write_text(json.dumps({
                'target_revision': 'b' * 40, 'headers': {}}))
            (root / 'src/lib.rs').write_text('pub struct Existing;')
            (root / 'SOURCE_PROVENANCE.json').write_text(json.dumps({
                'source': 'src/openms/source/FORMAT/Referenced.cpp'}))
            snapshot = {'files': [{'path': 'src/topp/Example.cpp', 'sha256': 'c' * 64,
                                  'includes': ['OpenMS/FORMAT/Existing.h', 'OpenMS/FORMAT/Missing.h',
                                               'OpenMS/PRODUCT/External.h']} ]}
            with patch.object(coverage, 'ROOT', root):
                data = coverage.build(snapshot)
            self.assertEqual(data['counts']['registered_public_headers'], 3)
            self.assertEqual([h['status'] for h in data['headers']],
                             ['evidence_requires_review', 'unmapped', 'evidence_requires_review'])
            self.assertEqual(len(data['tools'][0]['open_sdk_headers']), 2)
            self.assertEqual(data['tools'][0]['outside_registered_sdk'], ['OpenMS/PRODUCT/External.h'])
            self.assertFalse(data['tools'][0]['workflow_validated'])

    def test_real_inventory_accounts_for_every_registered_header_once(self):
        snapshot = json.loads((coverage.ROOT / 'docs/topp-source-inventory.json').read_text())
        data = coverage.build(snapshot)
        inventory = json.loads((coverage.ROOT / 'docs/core-sdk-update.json').read_text())
        expected = {h['path'] for h in inventory['files'] if h['registration'] == 'registered_public_header'}
        self.assertEqual({h['header'] for h in data['headers']}, expected)
        self.assertEqual(len(data['headers']), len(expected))
        self.assertGreater(data['counts']['by_status'].get('unmapped', 0), 0)


if __name__ == '__main__':
    unittest.main()
