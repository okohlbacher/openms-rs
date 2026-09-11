#!/usr/bin/env python3
"""Project the source header unchanged except declared child-count repairs.
No native Rust output supplies expected values; assertions are upstream literals.
"""
import hashlib,json,re,sys
from pathlib import Path
import xml.etree.ElementTree as ET
source=Path(sys.argv[1])/'src/tests/class_tests/openms/data/MzMLFile_1.mzML'
out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
raw=source.read_bytes(); (out/'MzMLFile_1.original.mzML').write_bytes(raw)
text=raw.decode('ascii'); prefix=text[:text.index('    <spectrumList')]
repairs=[]
for tag,child in [('sampleList','sample'),('softwareList','software'),('dataProcessingList','dataProcessing')]:
 m=re.search(r'<'+tag+r' count="(\d+)">(.*?)</'+tag+r'>',prefix,re.S)
 actual=len(re.findall(r'<'+child+r'(?=\s|>)',m.group(2)))
 if int(m.group(1))!=actual:
  start=m.start(1);end=m.end(1);prefix=prefix[:start]+str(actual)+prefix[end:]
  repairs.append({'element':tag,'old':int(m.group(1)),'new':actual})
result=(prefix+'  </run>\n</mzML>\n').encode('ascii');ET.fromstring(result)
(out/'header.mzML').write_bytes(result)
(out/'header_projection.json').write_text(json.dumps({'source_commit':'82ce5b373c97f934ffd9b1ffd80215ca66473d0b','source_sha256':hashlib.sha256(raw).hexdigest(),'projection_sha256':hashlib.sha256(result).hexdigest(),'removed':'spectrumList and all following record payloads; source run/header parameters unchanged','declared_count_repairs':repairs},indent=2)+'\n')
print(repairs)
