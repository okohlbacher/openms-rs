#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# Independent SHA-1, byte-offset and XML identity check for native output.
import hashlib
import re
import sys
import xml.etree.ElementTree as ET

raw = sys.stdin.buffer.read()
ns = '{http://psi.hupo.org/ms/mzml}'
root = ET.fromstring(raw)
assert root.tag == ns+'indexedmzML'
end = raw.index(b'<fileChecksum>') + len(b'<fileChecksum>')
actual = root.find(ns+'fileChecksum').text
assert actual == hashlib.sha1(raw[:end]).hexdigest(), 'SHA-1 prefix differs'
assert len(actual) == 40
index_offset = int(root.find(ns+'indexListOffset').text)
assert raw[index_offset:].startswith(b'<indexList ')
index_list = root.find(ns+'indexList')
assert int(index_list.attrib['count']) == len(index_list)
expected = []
for match in re.finditer(rb'<(spectrum|chromatogram) ', raw):
    offset = match.start()
    tag = match.group(1)
    closing = raw.index(b'>', offset)
    fragment = ET.fromstring(raw[offset:closing+1] + b'</' + tag + b'>')
    expected.append((tag.decode(), fragment.attrib['id'], offset))
actual_rows = []
seen_sections = set()
for index in index_list:
    assert index.attrib['name'] not in seen_sections
    seen_sections.add(index.attrib['name'])
    assert len(index) > 0
    tag = index.attrib['name'].encode()
    for row in index:
        offset = int(row.text)
        assert offset >= 0
        actual_rows.append((tag.decode(), row.attrib['idRef'], offset))
        assert raw[offset:].startswith(b'<' + tag + b' ')
        closing = raw.index(b'>', offset)
        fragment = ET.fromstring(raw[offset:closing+1] + b'</' + tag + b'>')
        assert fragment.attrib['id'] == row.attrib['idRef']
assert actual_rows == expected, 'index is not a complete ordered bijection'
assert len({(k, i) for k, i, _ in actual_rows}) == len(actual_rows)
assert len({o for _, _, o in actual_rows}) == len(actual_rows)
assert seen_sections == {k for k, _, _ in expected}
print(end % 64)
