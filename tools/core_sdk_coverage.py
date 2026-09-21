#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Build a conservative, exhaustive SDK header/TOPP dependency work ledger.

Source references and matching Rust declarations are evidence to review, never
proof of complete C++ API parity. Only explicit reviewed entries can be complete.
The cached TOPP snapshot makes --check usable without either C++ checkout.

The ledger counts the installed public headers of every package in `PACKAGES`:
core, whose registration is the `sources.cmake` union recorded by
tools/core_sdk_retarget.py, and cli, which installs its `include/` directory
wholesale and is recorded by `cli_inventory` from the pinned git objects.
"""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
INCLUDE = re.compile(r'^\s*#\s*include\s*[<"](OpenMS/[^">]+)[">]', re.MULTILINE)
# Each package's inventory, and the revision docs/core-sdk-reviewed-apis.json
# pins its reviews to under the same name in `target_revisions`.
PACKAGES = {'core': 'docs/core-sdk-update.json', 'cli': 'docs/cli-sdk-update.json'}
# Where the ledger keys a cli file. The prefix keeps a segment in front of
# `include/`, which the include-path split in build() needs, and it is the path
# the TOPP manifests already cite with the sibling checkout's name in front.
CLI_ROOT = 'packages/cli/'
CLI_CITED_AS = 'OpenMS4-tests/' + CLI_ROOT
CLI_REPOSITORY = 'https://github.com/okohlbacher/OpenMS4-cli.git'
# The two CMakeLists.txt rules the cli registration is read from. If a retarget
# finds either one changed, the rule in cli_inventory has to be re-derived.
CLI_RULES = {'install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})': 'registered_public_header',
             'file(GLOB cli_sources CONFIGURE_DEPENDS source/APPLICATIONS/*.cpp)': 'implementation'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physical_lines(data):
    """Lines as tools/core_sdk_retarget.py counts them: a last line without a newline counts."""
    return 0 if not data else data.count(b'\n') + (0 if data.endswith(b'\n') else 1)


def cli_registration(path):
    """The (kind, registration) CLI_RULES give a package-relative cli path, or None outside SDK scope."""
    if path.startswith('include/'):
        return 'cli_header', 'registered_public_header'
    if re.fullmatch(r'source/APPLICATIONS/[^/]+\.cpp', path):
        return 'cli_implementation', 'implementation'
    return None


def cli_inventory(checkout, revision):
    """Record the cli package at `revision`, read from its git objects, not its working tree.

    Registration here is not a `sources.cmake` union. CMakeLists.txt installs
    `include/` as a directory, so every file under it is an installed public
    header, and compiles `source/APPLICATIONS/*.cpp` by glob. The ninth
    installed header, OpenMSCLIConfig.h, is written by generate_export_header
    into the build tree, so it is in no revision's tree; it belongs with the
    generated configuration headers core's inventory also leaves out.
    """
    def git(*arguments):
        return subprocess.check_output(['git', '-C', str(checkout), *arguments])

    revision = git('rev-parse', '--verify', f'{revision}^{{commit}}').decode().strip()
    tree = [path for path in git('ls-tree', '-r', '-z', '--name-only', revision).decode().split('\0') if path]
    cmake = git('cat-file', 'blob', f'{revision}:CMakeLists.txt')
    for rule in CLI_RULES:
        assert rule.encode() in cmake, f'cli CMakeLists.txt at {revision} no longer says {rule}'
    assert not any(path.endswith('/OpenMSCLIConfig.h') for path in tree), 'generated header committed'
    files = []
    for path in tree:
        if cli_registration(path) is None:
            continue
        kind, registration = cli_registration(path)
        data = git('cat-file', 'blob', f'{revision}:{path}')
        files.append({'path': CLI_ROOT + path, 'kind': kind, 'domain': path.split('/')[-2],
                      'bytes': len(data), 'physical_lines': physical_lines(data),
                      'sha256': hashlib.sha256(data).hexdigest(), 'registration': registration,
                      'registration_file': CLI_ROOT + 'CMakeLists.txt'})
    tests = [git('cat-file', 'blob', f'{revision}:{path}') for path in tree
             if re.fullmatch(r'tests/source/[^/]+_test\.cpp', path)]
    by_kind = {}
    for item in files:
        entry = by_kind.setdefault(item['kind'], {'files': 0, 'physical_lines': 0, 'bytes': 0})
        entry['files'] += 1
        entry['physical_lines'] += item['physical_lines']
        entry['bytes'] += item['bytes']
    return {'schema_version': 1,
            'identity': {'package': 'cli', 'current_package_revision': revision, 'repository': CLI_REPOSITORY,
                         'checkout_kind': 'git objects at the pinned revision; no working tree is read'},
            'methodology': [
                'Regenerate from the pin with `python3 tools/core_sdk_coverage.py --cli-source <checkout> --write`; without --write the same option checks this record against the pin.',
                'Registration is read from CMakeLists.txt, not from a sources.cmake union: `install(DIRECTORY include/ ...)` installs every file under include/, unconditionally, and `file(GLOB cli_sources ... source/APPLICATIONS/*.cpp)` compiles every implementation.',
                'OpenMSCLIConfig.h is generated by generate_export_header into the build tree and installed beside the headers; like core\'s generated configuration headers it is not inventoried.',
                'Class tests, fixtures, the tool registry fixture and build files are package files outside SDK scope; they are counted in the summary, not listed.'],
            'summary': {'by_kind': dict(sorted(by_kind.items())),
                        'registration': dict(sorted(Counter(item['registration'] for item in files).items())),
                        'package_files': len(tree),
                        'class_test_cpp_files': len(tests),
                        'class_test_cpp_lines': sum(physical_lines(data) for data in tests)},
            'files': files,
            'registration_evidence': [{'path': CLI_ROOT + 'CMakeLists.txt', 'sha256': hashlib.sha256(cmake).hexdigest(),
                                       'rules': dict(CLI_RULES)}]}


def check_cli_inventory(inventory):
    """What a hand edit could break without a checkout to regenerate against."""
    assert inventory['identity']['package'] == 'cli'
    for item in inventory['files']:
        path = item['path']
        assert path.startswith(CLI_ROOT), path
        assert cli_registration(path.removeprefix(CLI_ROOT)) == (item['kind'], item['registration']), \
            f'{path} is not registered as the CMakeLists.txt rules say'
        assert re.fullmatch(r'[0-9a-f]{64}', item['sha256']), path
    assert Counter(item['registration'] for item in inventory['files']) == inventory['summary']['registration']


def topp_snapshot(source):
    source = source.resolve()
    git_root = Path(subprocess.check_output(['git', '-C', str(source), 'rev-parse', '--show-toplevel'], text=True).strip())
    revision = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    records = []
    for path in sorted(source.glob('*.cpp')):
        records.append({'path': path.relative_to(git_root).as_posix(), 'sha256': digest(path),
                        'includes': sorted(set(INCLUDE.findall(path.read_text())))})
    return {'repository_revision': revision,
            'scope': 'All local src/topp/*.cpp files; registration conditions not evaluated. Per-file hashes identify working-tree content independently of the Git revision.',
            'files': records}


def referenced_paths(value):
    if isinstance(value, str):
        if value.startswith('src/openms/') or value.startswith('src/openswathalgo/'):
            yield value
        elif value.startswith(CLI_CITED_AS):
            yield value.removeprefix('OpenMS4-tests/')
    elif isinstance(value, dict):
        # A manifest may instead name the package beside a package-relative
        # path. Only cli's is resolved: its root holds include/ and source/
        # directly, so the path means one file.
        if value.get('package') == 'cli' and isinstance(value.get('path'), str):
            yield CLI_ROOT + value['path']
        for key, item in value.items():
            yield from referenced_paths(key)
            yield from referenced_paths(item)
    elif isinstance(value, list):
        for item in value:
            yield from referenced_paths(item)


def validated_workflows():
    """Tools with executed differential evidence, from the TOPP provenance manifests.

    A manifest qualifies only when it declares tier 1 and names the upstream
    test definition it reproduces; a tool whose manifest records source review
    is not a validated workflow.
    """
    provenance = json.loads((ROOT / 'SOURCE_PROVENANCE.json').read_text())
    found = {}
    for manifest in provenance.get('topp_package_reference_manifests', []):
        data = json.loads((ROOT / manifest).read_text())
        # str(): a manifest that records the tier as a bare number states no tier
        # in the sense this rule means, and must not crash the check either.
        if 'tier 1' not in str(data.get('evidence_tier', '')) or not data.get('upstream_test_definition'):
            continue
        for path in data.get('native_implementation', []):
            if path.startswith('src/bin/') and path.endswith('.rs'):
                found[Path(path).stem] = manifest
    return found


def build(snapshot):
    reviewed = json.loads((ROOT / 'docs/core-sdk-reviewed-apis.json').read_text())
    assert set(reviewed['target_revisions']) == set(PACKAGES)
    headers = []
    for package, record in PACKAGES.items():
        inventory = json.loads((ROOT / record).read_text())
        assert reviewed['target_revisions'][package] == inventory['identity']['current_package_revision'], package
        if package == 'cli':
            check_cli_inventory(inventory)
        headers += [dict(item, package=package) for item in inventory['files']
                    if item['registration'] == 'registered_public_header']
    # Every package installs into the one OpenMS/ include namespace; core
    # already owns OpenMS/APPLICATIONS/ConsoleUtils.h beside cli's headers. A
    # tool's include names no package, so two packages registering one include
    # path would leave the ledger unable to say whose header a tool uses, and a
    # map keyed by it would silently keep only the last.
    by_include = {}
    for header in headers:
        assert '/include/OpenMS/' in header['path'], header['path']
        include = header['path'].split('/include/', 1)[1]
        assert include not in by_include, f"{include} is registered by {by_include[include]['package']} and {header['package']}"
        by_include[include] = header
    assert len({h['path'] for h in headers}) == len(headers)
    assert set(reviewed['headers']) <= {h['path'] for h in headers}
    declarations = defaultdict(set)
    for path in sorted((ROOT / 'src').rglob('*.rs')):
        for name in re.findall(r'pub\s+(?:struct|enum|type|trait)\s+(\w+)', path.read_text()):
            declarations[name].add(path.relative_to(ROOT).as_posix())
    evidence = defaultdict(set)
    manifests = [ROOT / 'SOURCE_PROVENANCE.json', *sorted((ROOT / 'tests/data').glob('*provenance.json'))]
    for path in manifests:
        for reference in referenced_paths(json.loads(path.read_text())):
            evidence[reference].add(path.relative_to(ROOT).as_posix())
    consumers = defaultdict(list)
    tools = []
    for entry in snapshot['files']:
        name = Path(entry['path']).stem
        for include in entry['includes']:
            if include in by_include:
                consumers[include].append(name)
        tools.append({'name': name, 'source_sha256': entry['sha256'],
                      'sdk_headers': [i for i in entry['includes'] if i in by_include],
                      'outside_registered_sdk': [i for i in entry['includes'] if i not in by_include]})
    rows = []
    for header in sorted(headers, key=lambda h: h['path']):
        path = header['path']
        include = path.split('/include/', 1)[1]
        implementation = path.replace('/include/OpenMS/', '/source/').removesuffix('.h') + '.cpp'
        refs = sorted(evidence[path] | evidence[implementation])
        rust = sorted(declarations[Path(path).stem])
        review = reviewed['headers'].get(path)
        status = review['status'] if review else ('evidence_requires_review' if refs or rust else 'unmapped')
        if review:
            assert status in {'complete', 'partial', 'native_equivalent'}
            for link in [*review['rust'], *review['tests'], review['documentation']]:
                assert (ROOT / link).is_file(), link
            assert review['tests'] and review['scope'], path
        rows.append({'header': path, 'package': header['package'], 'sha256': header['sha256'], 'domain': header['domain'],
                     'status': status, 'review': review, 'candidate_rust_files': rust,
                     'reference_manifests': refs, 'direct_topp_consumers': sorted(consumers[include])})
    states = {row['header'].split('/include/', 1)[1]: row['status'] for row in rows}
    assert len(states) == len(rows)
    validated = validated_workflows()
    for tool in tools:
        tool['open_sdk_headers'] = [h for h in tool['sdk_headers'] if states[h] not in {'complete', 'native_equivalent'}]
        # A direct include match never establishes behavioral tool compatibility.
        # Only an executed differential comparison against retained C++ output
        # does, which a TOPP provenance manifest records as tier 1.
        tool['workflow_validated'] = tool['name'] in validated
        if tool['workflow_validated']:
            tool['workflow_evidence'] = validated[tool['name']]
    by_package = {package: {'registered_public_headers': sum(1 for r in rows if r['package'] == package),
                            'by_status': dict(sorted(Counter(r['status'] for r in rows if r['package'] == package).items()))}
                  for package in PACKAGES}
    return {'schema_version': 1, 'target_revisions': reviewed['target_revisions'],
            'methodology': 'Every registered public SDK header of the core and cli packages is accounted for. Direct include dependencies and candidate Rust declarations guide review; neither source hashes nor names demonstrate method coverage. Transitive dependencies, conditional configurations, runtime data and tool workflows need separate validation; in particular a TOPP source that textually includes a sibling source (the FeatureLinker tools include FeatureLinkerBase.cpp) is credited only with its own OpenMS/ includes. External/product-owned headers are not added to SDK scope.',
            'counts': {'registered_public_headers': len(rows), 'by_status': dict(sorted(Counter(r['status'] for r in rows).items())),
                       'by_package': by_package,
                       'topp_sources': len(tools),
                       'validated_topp_workflows': sum(1 for t in tools if t['workflow_validated'])},
            'headers': rows, 'tools': tools}


def markdown(data):
    counts = data['counts']
    packages = counts['by_package']
    targets = ', '.join(f"{package} `{revision}`" for package, revision in data['target_revisions'].items())
    split = ' and '.join(f"{c['registered_public_headers']} {package}" for package, c in packages.items())
    lines = ['# SDK completion ledger: core and cli packages', '',
             f"Targets: {targets}. This ledger covers all **{counts['registered_public_headers']} registered public headers** of those packages ({split}) and direct includes from **{counts['topp_sources']} TOPP source files**.", '',
             'This is a work inventory, not a completion percentage. Matching declarations and source references remain unverified until each API and its behavior are reviewed. A TOPP workflow counts as validated only when an executed differential comparison against retained C++ output is recorded in its provenance manifest. Physical unregistered headers and product backends are tracked separately by the SDK source inventory.', '',
             f"Validated TOPP workflows: **{counts['validated_topp_workflows']}** of {counts['topp_sources']}, each reproducing its upstream test against retained C++ output.", '',
             '| Review state | ' + ' | '.join(packages) + ' | Headers |', '| --- |' + ' ---: |' * (len(packages) + 1)]
    lines.extend(f"| {status} | " + ' | '.join(str(c['by_status'].get(status, 0)) for c in packages.values()) + f' | {count} |'
                 for status, count in counts['by_status'].items())
    lines += ['', '## Highest fan-out open SDK dependencies', '',
              'These counts show direct consumers; they do not establish full dependency closure.', '',
              '| Header | Package | Direct TOPP consumers | State |', '| --- | --- | ---: | --- |']
    priorities = sorted((r for r in data['headers'] if r['status'] not in {'complete', 'native_equivalent'}),
                        key=lambda r: (-len(r['direct_topp_consumers']), r['header']))
    for row in priorities[:35]:
        short = row['header'].split('/include/', 1)[1]
        lines.append(f"| `{short}` | {row['package']} | {len(row['direct_topp_consumers'])} | {row['status']} |")
    lines += ['', '## Completion requirements', '',
              '1. Review every public API against the exact source revision and configuration; record a tested native implementation or an explicit standard-library equivalent.',
              '2. Finish partial containers, formats, numerical routines and domain algorithms; preserve required metadata and error behavior.',
              '3. Check transitive tool dependencies, resources and platform backends, then run ported TOPP workflows against C++ results.',
              '4. Validate supported Rust versions, feature combinations, release builds and performance on realistic inputs.', '',
              'The complete per-header and per-tool records are in [core-sdk-coverage.json](core-sdk-coverage.json). Explicit reviews are maintained in [core-sdk-reviewed-apis.json](core-sdk-reviewed-apis.json); the package inventories are [core-sdk-update.json](core-sdk-update.json) and [cli-sdk-update.json](cli-sdk-update.json); TOPP source hashes and includes are in [topp-source-inventory.json](topp-source-inventory.json).', '',
              'Regenerate with `python3 tools/core_sdk_coverage.py --write`. CI checks that the ledger agrees with the package inventories, current Rust declarations and evidence. To replace the local TOPP snapshot, pass `--topp-source /path/to/OpenMS/src/topp --write`. To check the cli inventory against its pin, pass `--cli-source /path/to/OpenMS4-cli`, and add `--write` to regenerate it.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--topp-source', type=Path)
    parser.add_argument('--cli-source', type=Path,
                        help='A cli package checkout holding the pinned revision; checks docs/cli-sdk-update.json against it, or with --write regenerates it')
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--require-complete', action='store_true',
                        help='Fail while any registered header lacks complete reviewed coverage')
    args = parser.parse_args()
    snapshot_path = ROOT / 'docs/topp-source-inventory.json'
    if args.topp_source:
        assert args.write, '--topp-source requires --write'
        snapshot_path.write_text(json.dumps(topp_snapshot(args.topp_source), indent=2) + '\n')
    if args.cli_source:
        revision = json.loads((ROOT / 'docs/core-sdk-reviewed-apis.json').read_text())['target_revisions']['cli']
        record = json.dumps(cli_inventory(args.cli_source, revision), indent=2) + '\n'
        path = ROOT / PACKAGES['cli']
        if args.write:
            path.write_text(record)
        else:
            assert path.read_text() == record, f'Stale {PACKAGES["cli"]}: regenerate with --cli-source <checkout> --write'
    result = build(json.loads(snapshot_path.read_text()))
    outputs = {'docs/core-sdk-coverage.json': json.dumps(result, indent=2) + '\n',
               'docs/CORE_SDK_COMPLETION.md': markdown(result)}
    for name, content in outputs.items():
        path = ROOT / name
        if args.write:
            path.write_text(content)
        else:
            assert path.read_text() == content, f'Stale {name}: regenerate with --write'
    print(json.dumps(result['counts'], sort_keys=True))
    if args.require_complete:
        remaining = [r for r in result['headers'] if r['status'] not in {'complete', 'native_equivalent'}]
        if remaining:
            raise SystemExit(f"SDK incomplete: {len(remaining)} headers still require review or implementation")


if __name__ == '__main__':
    main()
