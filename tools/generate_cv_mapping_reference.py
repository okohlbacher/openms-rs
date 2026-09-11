#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
"""Independent XML projection of complete source mapping fixtures; no Rust use."""
from pathlib import Path
import argparse, hashlib, xml.etree.ElementTree as ET
p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args()
base=Path(__file__).resolve().parents[1]/'tests/data/cv_mapping'
EXPECTED = {'cv_mapping_test_file.xml': '7ecf60e80a8b148377030c0b49679fdc385791acd857fec53b1f53ed1113e4b6', 'SemanticValidator_mapping.xml': '52fdb344da5250f6d99b461ecc82dc41e4e43678461c0e0be730505c1902b24a', 'ms-mapping.xml': '8407daa22c49e61fb4d6c5a965899a158bba96cb8f4783dd5fbfe3b097dc18ff', 'mzdata-mapping.xml': '2613831cc48e961d1e7e0a06a116885441501c8d11306e28cf42b53e139acfef', 'mzIdentML-mapping.xml': 'add96df84db7219120176113cf77d240b637237c4a16ee64d628461e8fb6be4a', 'TraML-mapping.xml': '141ed4eda7a183caff55e4cfececbd73e69aed46899f3e3b427d17604a6213d4'}
rows=[]
h=lambda s:s.encode().hex()
for filename in sorted(EXPECTED):
    f=base/filename
    raw=f.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == EXPECTED[filename], filename
    root=ET.fromstring(raw);ri=0
    for ref in root.iter('CvReference'):
        rows.append('\t'.join([f.name,'reference',h(ref.attrib['cvName']),h(ref.attrib['cvIdentifier'])]))
    for rule in root.iter('CvMappingRule'):
        x=rule.attrib
        rows.append('\t'.join([f.name,'rule',str(ri),h(x['id']),h(x['cvElementPath']),x['requirementLevel'],h(x['scopePath']),x['cvTermsCombinationLogic']]))
        for ti,term in enumerate(rule.iter('CvTerm')):
            x=term.attrib
            rows.append('\t'.join([f.name,'term',str(ri),str(ti),h(x['termAccession']),x.get('useTermName') or 'false',x['useTerm'],h(x['termName']),x.get('isRepeatable') or 'true',x['allowChildren'],h(x['cvIdentifierRef'])]))
        ri+=1
text='\n'.join(rows)+'\n';out=base/'projection.tsv'
if a.check:assert out.read_text()==text
else:out.write_text(text)
print(f'{len(rows)} exact complete XML record rows verified')
