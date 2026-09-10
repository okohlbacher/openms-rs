#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Build a conservative, exhaustive SDK header/TOPP dependency work ledger.

Source references and matching Rust declarations are evidence to review, never
proof of complete C++ API parity. Only explicit reviewed entries can be complete.
The cached TOPP snapshot makes --check usable without either C++ checkout.
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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from referenced_paths(key)
            yield from referenced_paths(item)
    elif isinstance(value, list):
        for item in value:
            yield from referenced_paths(item)


def build(snapshot):
    inventory = json.loads((ROOT / 'docs/core-sdk-update.json').read_text())
    headers = [item for item in inventory['files'] if item['registration'] == 'registered_public_header']
    reviewed = json.loads((ROOT / 'docs/core-sdk-reviewed-apis.json').read_text())
    assert reviewed['target_revision'] == inventory['identity']['current_package_revision']
    by_include = {h['path'].split('/include/', 1)[1]: h for h in headers}
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
        rows.append({'header': path, 'sha256': header['sha256'], 'domain': header['domain'],
                     'status': status, 'review': review, 'candidate_rust_files': rust,
                     'reference_manifests': refs, 'direct_topp_consumers': sorted(consumers[include])})
    states = {row['header'].split('/include/', 1)[1]: row['status'] for row in rows}
    for tool in tools:
        tool['open_sdk_headers'] = [h for h in tool['sdk_headers'] if states[h] not in {'complete', 'native_equivalent'}]
        # A direct include match never establishes behavioral tool compatibility.
        tool['workflow_validated'] = False
    return {'schema_version': 1, 'target_revision': reviewed['target_revision'],
            'methodology': 'Every registered public SDK header is accounted for. Direct include dependencies and candidate Rust declarations guide review; neither source hashes nor names demonstrate method coverage. Transitive dependencies, conditional configurations, runtime data and tool workflows need separate validation. External/product-owned headers are not added to Core scope.',
            'counts': {'registered_public_headers': len(rows), 'by_status': dict(sorted(Counter(r['status'] for r in rows).items())),
                       'topp_sources': len(tools), 'validated_topp_workflows': 0},
            'headers': rows, 'tools': tools}


def markdown(data):
    counts = data['counts']
    lines = ['# Core SDK completion ledger', '',
             f"Target: `{data['target_revision']}`. This ledger covers all **{counts['registered_public_headers']} registered public headers** and direct includes from **{counts['topp_sources']} TOPP source files**.", '',
             'This is a work inventory, not a completion percentage. Matching declarations and source references remain unverified until each API and its behavior are reviewed. No TOPP workflow is yet certified as port-ready. Physical unregistered headers and product backends are tracked separately by the SDK source inventory.', '',
             '| Review state | Headers |', '| --- | ---: |']
    lines.extend(f'| {status} | {count} |' for status, count in counts['by_status'].items())
    lines += ['', '## Highest fan-out open SDK dependencies', '',
              'These counts show direct consumers; they do not establish full dependency closure.', '',
              '| Header | Direct TOPP consumers | State |', '| --- | ---: | --- |']
    priorities = sorted((r for r in data['headers'] if r['status'] not in {'complete', 'native_equivalent'}),
                        key=lambda r: (-len(r['direct_topp_consumers']), r['header']))
    for row in priorities[:35]:
        short = row['header'].split('/include/', 1)[1]
        lines.append(f"| `{short}` | {len(row['direct_topp_consumers'])} | {row['status']} |")
    lines += ['', '## Completion requirements', '',
              '1. Review every public API against the exact source revision and configuration; record a tested native implementation or an explicit standard-library equivalent.',
              '2. Finish partial containers, formats, numerical routines and domain algorithms; preserve required metadata and error behavior.',
              '3. Check transitive tool dependencies, resources and platform backends, then run ported TOPP workflows against C++ results.',
              '4. Validate supported Rust versions, feature combinations, release builds and performance on realistic inputs.', '',
              'The complete per-header and per-tool records are in [core-sdk-coverage.json](core-sdk-coverage.json). Explicit reviews are maintained in [core-sdk-reviewed-apis.json](core-sdk-reviewed-apis.json); TOPP source hashes and includes are in [topp-source-inventory.json](topp-source-inventory.json).', '',
              'Regenerate with `python3 tools/core_sdk_coverage.py --write`. CI checks that the ledger agrees with the source inventory, current Rust declarations and evidence. To replace the local TOPP snapshot, pass `--topp-source /path/to/OpenMS/src/topp --write`.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--topp-source', type=Path)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--require-complete', action='store_true',
                        help='Fail while any registered header lacks complete reviewed coverage')
    args = parser.parse_args()
    snapshot_path = ROOT / 'docs/topp-source-inventory.json'
    if args.topp_source:
        assert args.write, '--topp-source requires --write'
        snapshot_path.write_text(json.dumps(topp_snapshot(args.topp_source), indent=2) + '\n')
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
            raise SystemExit(f"Core SDK incomplete: {len(remaining)} headers still require review or implementation")


if __name__ == '__main__':
    main()
