#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Guard the coverage ledger against false completion and dropped SDK scope."""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import check_source_citations
import core_sdk_coverage as coverage

CORE, CLI = 'b' * 40, 'd' * 40


def header(path, registration='registered_public_header'):
    return {'path': path, 'sha256': 'a' * 64, 'domain': path.split('/')[-2], 'registration': registration}


def write_fixture(root, core_files, cli_files, provenance, revisions=None):
    """A repository with a core and a cli inventory, and nothing reviewed."""
    (root / 'docs').mkdir()
    (root / 'src').mkdir()
    (root / 'tests/data').mkdir(parents=True)
    (root / 'docs/core-sdk-update.json').write_text(json.dumps({
        'identity': {'current_package_revision': CORE}, 'files': core_files}))
    registration = {}
    for item in cli_files:
        registration[item['registration']] = registration.get(item['registration'], 0) + 1
    (root / 'docs/cli-sdk-update.json').write_text(json.dumps({
        'identity': {'package': 'cli', 'current_package_revision': CLI}, 'files': cli_files,
        'summary': {'registration': registration}}))
    (root / 'docs/core-sdk-reviewed-apis.json').write_text(json.dumps({
        'target_revisions': revisions or {'core': CORE, 'cli': CLI}, 'headers': {}}))
    (root / 'SOURCE_PROVENANCE.json').write_text(json.dumps(provenance))


def tool(*includes):
    return {'files': [{'path': 'src/topp/Example.cpp', 'sha256': 'c' * 64, 'includes': list(includes)}]}


class CoverageTests(unittest.TestCase):
    def test_names_and_source_references_never_imply_parity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            files = [header(f'src/openms/include/OpenMS/FORMAT/{name}.h') for name in ['Existing', 'Missing', 'Referenced']]
            write_fixture(root, files, [], {'source': 'src/openms/source/FORMAT/Referenced.cpp'})
            (root / 'src/lib.rs').write_text('pub struct Existing;')
            snapshot = tool('OpenMS/FORMAT/Existing.h', 'OpenMS/FORMAT/Missing.h', 'OpenMS/PRODUCT/External.h')
            with patch.object(coverage, 'ROOT', root):
                data = coverage.build(snapshot)
            self.assertEqual(data['counts']['registered_public_headers'], 3)
            self.assertEqual([h['status'] for h in data['headers']],
                             ['evidence_requires_review', 'unmapped', 'evidence_requires_review'])
            self.assertEqual(len(data['tools'][0]['open_sdk_headers']), 2)
            self.assertEqual(data['tools'][0]['outside_registered_sdk'], ['OpenMS/PRODUCT/External.h'])
            self.assertFalse(data['tools'][0]['workflow_validated'])

    def test_a_cli_header_counts_against_the_tools_that_include_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cli = 'packages/cli/include/OpenMS/APPLICATIONS/'
            files = [header(cli + f'{name}.h') for name in ['Absent', 'Cited', 'Prefixed', 'Reviewed']]
            files.append(header('packages/cli/source/APPLICATIONS/Cited.cpp', 'implementation'))
            # The two ways a manifest cites the cli package: the sibling
            # checkout's path, and a package name beside a package-relative path.
            provenance = {'sources': [
                {'path': 'OpenMS4-tests/packages/cli/include/OpenMS/APPLICATIONS/Prefixed.h'},
                {'package': 'cli', 'revision': CLI, 'path': 'source/APPLICATIONS/Cited.cpp'},
                {'package': 'topp', 'path': 'include/OpenMS/APPLICATIONS/Absent.h'}]}
            write_fixture(root, [header('src/openms/include/OpenMS/APPLICATIONS/ConsoleUtils.h')], files, provenance)
            for link in ['src/cli.rs', 'tests/cli.rs', 'docs/CLI.md']:
                (root / link).write_text('')
            reviewed = json.loads((root / 'docs/core-sdk-reviewed-apis.json').read_text())
            reviewed['headers'][cli + 'Reviewed.h'] = {'status': 'partial', 'rust': ['src/cli.rs'], 'tests': ['tests/cli.rs'],
                                                       'documentation': 'docs/CLI.md', 'scope': 'Partly.'}
            (root / 'docs/core-sdk-reviewed-apis.json').write_text(json.dumps(reviewed))
            snapshot = tool('OpenMS/APPLICATIONS/Reviewed.h', 'OpenMS/APPLICATIONS/Absent.h',
                            'OpenMS/APPLICATIONS/ConsoleUtils.h', 'OpenMS/APPLICATIONS/Unregistered.h')
            with patch.object(coverage, 'ROOT', root):
                data = coverage.build(snapshot)
            self.assertEqual(data['counts']['registered_public_headers'], 5)
            self.assertEqual(data['counts']['by_package']['cli'],
                             {'registered_public_headers': 4,
                              'by_status': {'evidence_requires_review': 2, 'partial': 1, 'unmapped': 1}})
            self.assertEqual(data['target_revisions'], {'core': CORE, 'cli': CLI})
            rows = {row['header'].rsplit('/', 1)[1]: row for row in data['headers']}
            self.assertEqual({name: row['status'] for name, row in rows.items()},
                             {'Absent.h': 'unmapped', 'Cited.h': 'evidence_requires_review', 'ConsoleUtils.h': 'unmapped',
                              'Prefixed.h': 'evidence_requires_review', 'Reviewed.h': 'partial'})
            self.assertEqual(rows['Reviewed.h']['package'], 'cli')
            self.assertEqual(rows['ConsoleUtils.h']['package'], 'core')
            self.assertEqual(rows['Reviewed.h']['direct_topp_consumers'], ['Example'])
            example = data['tools'][0]
            self.assertEqual(example['outside_registered_sdk'], ['OpenMS/APPLICATIONS/Unregistered.h'])
            self.assertEqual(sorted(example['open_sdk_headers']),
                             ['OpenMS/APPLICATIONS/Absent.h', 'OpenMS/APPLICATIONS/ConsoleUtils.h', 'OpenMS/APPLICATIONS/Reviewed.h'])

    def test_two_packages_may_not_register_one_include_path(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fixture(root, [header('src/openms/include/OpenMS/APPLICATIONS/ConsoleUtils.h')],
                          [header('packages/cli/include/OpenMS/APPLICATIONS/ConsoleUtils.h')], {})
            with patch.object(coverage, 'ROOT', root), self.assertRaisesRegex(AssertionError, 'registered by core and cli'):
                coverage.build(tool('OpenMS/APPLICATIONS/ConsoleUtils.h'))

    def test_each_package_is_held_to_its_own_revision(self):
        for revisions in [{'core': CORE, 'cli': CORE}, {'core': CORE}, {'core': CORE, 'cli': CLI, 'topp': CLI}]:
            with self.subTest(revisions=revisions), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                write_fixture(root, [], [], {}, revisions)
                with patch.object(coverage, 'ROOT', root), self.assertRaises(AssertionError):
                    coverage.build(tool())

    def test_a_hand_edited_cli_registration_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fixture(root, [], [header('packages/cli/source/APPLICATIONS/TOPPBase.cpp')], {})
            with patch.object(coverage, 'ROOT', root), self.assertRaisesRegex(AssertionError, 'install rule'):
                coverage.build(tool())

    def test_cli_inventory_reads_the_pinned_objects_by_the_install_rule(self):
        with tempfile.TemporaryDirectory() as directory:
            checkout = Path(directory)

            def git(*arguments):
                return subprocess.run(['git', '-C', str(checkout), '-c', 'user.name=t', '-c', 'user.email=t@t', '-c', 'commit.gpgsign=false',
                                       *arguments], check=True, capture_output=True, text=True).stdout.strip()

            files = {'CMakeLists.txt': '\n'.join(coverage.CLI_RULES) + '\n',
                     'include/OpenMS/APPLICATIONS/Base.h': 'one\ntwo\n',
                     'include/OpenMS/APPLICATIONS/Base_defs.h': 'no final newline',
                     'source/APPLICATIONS/Base.cpp': 'a\nb\nc\n',
                     'source/APPLICATIONS/detail/Helper.cpp': 'x\n',
                     'tests/source/Base_test.cpp': '1\n2\n3\n4\n',
                     'README.md': 'r\n'}
            git('init', '-q')
            for path, text in files.items():
                (checkout / path).parent.mkdir(parents=True, exist_ok=True)
                (checkout / path).write_text(text)
            git('add', '-A')
            git('commit', '-q', '-m', 'pin')
            pin = git('rev-parse', 'HEAD')
            # The working tree moves on; the record must still describe the pin.
            (checkout / 'include/OpenMS/APPLICATIONS/Base.h').write_text('changed\n')
            (checkout / 'include/OpenMS/APPLICATIONS/Later.h').write_text('')
            record = coverage.cli_inventory(checkout, pin[:7])
            self.assertEqual(record['identity']['current_package_revision'], pin)
            self.assertEqual([(f['path'], f['registration'], f['physical_lines']) for f in record['files']],
                             [('packages/cli/include/OpenMS/APPLICATIONS/Base.h', 'registered_public_header', 2),
                              ('packages/cli/include/OpenMS/APPLICATIONS/Base_defs.h', 'registered_public_header', 1),
                              ('packages/cli/source/APPLICATIONS/Base.cpp', 'implementation', 3)])
            self.assertEqual((record['summary']['package_files'], record['summary']['class_test_cpp_files'],
                              record['summary']['class_test_cpp_lines']), (7, 1, 4))
            coverage.check_cli_inventory(record)
            (checkout / 'CMakeLists.txt').write_text('install(FILES include/OpenMS/APPLICATIONS/Base.h)\n')
            git('commit', '-q', '-am', 'narrower install')
            with self.assertRaisesRegex(AssertionError, 'no longer says'):
                coverage.cli_inventory(checkout, 'HEAD')

    def test_real_inventory_accounts_for_every_registered_header_once(self):
        snapshot = json.loads((coverage.ROOT / 'docs/topp-source-inventory.json').read_text())
        data = coverage.build(snapshot)
        expected = set()
        for record in coverage.PACKAGES.values():
            inventory = json.loads((coverage.ROOT / record).read_text())
            expected |= {h['path'] for h in inventory['files'] if h['registration'] == 'registered_public_header'}
        self.assertEqual({h['header'] for h in data['headers']}, expected)
        self.assertEqual(len(data['headers']), len(expected))
        self.assertGreater(data['counts']['by_status'].get('unmapped', 0), 0)
        # Every APPLICATIONS header a TOPP source includes is one some package
        # installs; one left outside would be a package the ledger does not read.
        self.assertEqual([i for t in data['tools'] for i in t['outside_registered_sdk'] if i.startswith('OpenMS/APPLICATIONS/')], [])

    def test_real_cli_inventory_agrees_with_the_hashes_the_manifests_recorded(self):
        """The manifests hashed the cli sources when they ported from them; the inventory must agree."""
        inventory = {f['path']: f['sha256'] for f in json.loads((coverage.ROOT / coverage.PACKAGES['cli']).read_text())['files']}
        seen = []

        def walk(value):
            if isinstance(value, dict):
                if 'sha256' in value:
                    for path in coverage.referenced_paths({key: value[key] for key in ('package', 'path') if key in value}):
                        if path in inventory:
                            seen.append((path, value['sha256']))
                for item in value.values():
                    walk(item)
            elif isinstance(value, list):
                for item in value:
                    walk(item)

        for path in [coverage.ROOT / 'SOURCE_PROVENANCE.json', *sorted((coverage.ROOT / 'tests/data').glob('*provenance.json'))]:
            walk(json.loads(path.read_text()))
        self.assertGreater(len(seen), 0)
        self.assertEqual([item for item in seen if inventory[item[0]] != item[1]], [])

    def test_real_cli_inventory_regenerates_from_the_pin(self):
        revision = json.loads((coverage.ROOT / 'docs/core-sdk-reviewed-apis.json').read_text())['target_revisions']['cli']
        checkouts = [root.parent / 'OpenMS4-tests/packages/cli' for root in check_source_citations.checkout_roots()]
        reachable = [path for path in checkouts
                     if check_source_citations.git(path, 'cat-file', '-e', f'{revision}^{{commit}}') is not None]
        if not reachable:
            self.skipTest(f'no cli checkout holding {revision[:7]} beside this repository')
        record = json.dumps(coverage.cli_inventory(reachable[0], revision), indent=2) + '\n'
        self.assertEqual(record, (coverage.ROOT / coverage.PACKAGES['cli']).read_text())


if __name__ == '__main__':
    unittest.main()
