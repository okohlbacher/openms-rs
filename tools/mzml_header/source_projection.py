#!/usr/bin/env python3
"""Extract source enum/CV pairs; never imports or executes Rust production."""
import argparse, hashlib, json, re
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('source', type=Path)
p.add_argument('output', type=Path)
a = p.parse_args()
rel = 'src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp'
text = (a.source / rel).read_text()
start = text.index('else if (parent_tag == "instrumentConfiguration")')
end = text.index('else if (parent_tag == "processingMethod")', start)
read = text[start:end]
rows = []
pattern = r'accession == "(MS:\d+)"[^\{]*\{([^{}]*)\}'
for m in re.finditer(pattern, read):
    body = m.group(2)
    setter = re.search(r'\.set(\w+)\((\w+)::(\w+)::(\w+)\)', body)
    if setter:
        operation, owner, enum, member = setter.groups()
        rows.append(dict(accession=m.group(1), owner=owner, enum=enum, member=member,
                         setter=operation, line=text[:start+m.start()].count('\n')+1))
wstart = text.index('void MzMLHandler::writeInstrument_(')
wend = text.index('void MzMLHandler::writePrecursor_(', wstart)
write = text[wstart:wend]
written = []
for m in re.finditer(r'(?:in|so|ma|id)\.get\w+\(\) == (\w+)::(\w+)::(\w+)\)\s*\{([^{}]*)\}', write):
    cv = re.search(r'accession=\\"(MS:\d+)\\"', m.group(4))
    if cv:
        written.append(dict(owner=m.group(1), enum=m.group(2), member=m.group(3),
                            accession=cv.group(1), line=text[:wstart+m.start()].count('\n')+1))
rmap={(x['owner'],x['enum'],x['member']):x['accession'] for x in rows}
wmap={(x['owner'],x['enum'],x['member']):x['accession'] for x in written}
assert rmap == wmap, (rmap.keys()-wmap.keys(), wmap.keys()-rmap.keys())
a.output.mkdir(parents=True,exist_ok=True)
(a.output/'instrument_enum_pairs.json').write_text(json.dumps({'source_revision':'82ce5b373c97f934ffd9b1ffd80215ca66473d0b','source':rel,'sha256':hashlib.sha256(text.encode()).hexdigest(),'oracle':'Parsed independently from literal C++ reader setters and writer branches; both maps must match. No Rust production is invoked.','reader':rows,'writer':written},indent=2)+'\n')
print(f'{len(rows)} reader and {len(written)} writer enum/accession pairs agree')
# Header enum order is a second, independent oracle for native discriminants.
lines = ['owner\tenum\tmember\taccession\tordinal']
for row in rows:
    header = (a.source / ('src/openms/include/OpenMS/METADATA/' + row['owner'] + '.h')).read_text()
    body = re.search(r'enum class ' + row['enum'] + r'\s*\{(.*?)\}', header, re.S).group(1)
    body = re.sub(r'//[^\n]*', '', body)
    members = [part.strip().split('=')[0].strip() for part in body.split(',') if part.strip()]
    lines.append('\t'.join([row['owner'], row['enum'], row['member'], row['accession'], str(members.index(row['member']))]))
(a.output / 'instrument_enum_pairs.tsv').write_text('\n'.join(lines)+'\n')
