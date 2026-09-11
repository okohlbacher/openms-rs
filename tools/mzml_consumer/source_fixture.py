#!/usr/bin/env python3
"""Repair only stale source child counts; preserve all original scientific data.
Expected 4 spectra / 40 peaks / TIC 350 come from MzMLFile_test.cpp, not Rust.
"""
import hashlib, json, re, sys
from pathlib import Path
import xml.etree.ElementTree as ET
source = Path(sys.argv[1]) / 'src/tests/class_tests/openms/data/MzMLFile_1.mzML'
out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
raw = source.read_bytes(); text = raw.decode('ascii'); repairs = []
for tag, child in [('sampleList','sample'), ('softwareList','software'),
                   ('dataProcessingList','dataProcessing'),
                   ('binaryDataArrayList','binaryDataArray'), ('productList','product')]:
    pattern = r'(<'+tag+r'\b[^>]*\bcount=")(\d+)("[^>]*>)(.*?)(</'+tag+r'>)'
    def repair(m):
        actual = len(re.findall(r'<'+child+r'(?=\s|>)',m[4]))
        if int(m[2]) != actual:
            repairs.append({'element':tag, 'old':int(m[2]), 'new':actual})
        return m[1]+str(actual)+m[3]+m[4]+m[5]
    text = re.sub(pattern, repair, text, flags=re.S)
result = text.encode('ascii'); ET.fromstring(result)
(out/'source.mzML').write_bytes(result)
(out/'projection.json').write_text(json.dumps({
 'source_commit':'82ce5b373c97f934ffd9b1ffd80215ca66473d0b',
 'source_sha256':hashlib.sha256(raw).hexdigest(),
 'projection_sha256':hashlib.sha256(result).hexdigest(),
 'declared_count_repairs':repairs,
 'other_changes':[], 'source_test_literals':{'spectra':4,'peaks':40,'tic':350},
}, indent=2)+'\n')
print(repairs)
