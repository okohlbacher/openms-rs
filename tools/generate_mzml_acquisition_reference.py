#!/usr/bin/env python3
"""Source XML projection; never reads native algorithm output."""
from pathlib import Path
import copy
import xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[1]
NS = 'http://psi.hupo.org/ms/mzml'
ET.register_namespace('', NS)
def q(s): return '{' + NS + '}' + s
source = ET.parse(ROOT / 'tests/data/mzml_load_source_original.mzML').getroot()
out = ET.Element(q('mzML'), {'version': '1.1.0'})
files = source.find('.//' + q('sourceFileList'))
if files is not None:
    desc = ET.SubElement(out, q('fileDescription'))
    target = ET.SubElement(desc, q('sourceFileList'), {'count': files.attrib['count']})
    for file in files:
        ET.SubElement(target, q('sourceFile'), {k: file.attrib[k] for k in ['id','name','location']})
out.append(copy.deepcopy(source.find(q('referenceableParamGroupList'))))
run = ET.SubElement(out, q('run'), {'id':'acquisition_projection'})
spectra = ET.SubElement(run, q('spectrumList'), {'count':'4'})
scan_cvs = {'MS:1000016','MS:1000497','MS:1000927','MS:1002082','MS:1002083','MS:1002527','MS:1002528','MS:1002892','MS:1003057','MS:1003371','MS:1003394'}
for old in source.findall('.//' + q('spectrum')):
    item = ET.SubElement(spectra, q('spectrum'), {'id': old.attrib['id'], 'defaultArrayLength':'0'})
    for child in old:
        name = child.tag.removeprefix('{' + NS + '}')
        if name == 'referenceableParamGroupRef' or (name == 'cvParam' and child.get('accession') in ['MS:1000579','MS:1000580','MS:1000127','MS:1000128','MS:1000511','MS:1000130']):
            item.append(copy.deepcopy(child))
        elif name == 'productList':
            products = copy.deepcopy(child)
            # Literal source has count=1 but two products; its class test asserts
            # both. Only the declaration is repaired; quantity literals survive.
            products.set('count', str(len(products.findall(q('product')))))
            item.append(products)
        elif name == 'scanList':
            scans = ET.SubElement(item, q('scanList'), dict(child.attrib))
            for entry in child:
                if entry.tag != q('scan'):
                    scans.append(copy.deepcopy(entry))
                    continue
                scan = ET.SubElement(scans, q('scan'), dict(entry.attrib))
                for value in entry:
                    if value.tag in [q('userParam'),q('referenceableParamGroupRef'),q('scanWindowList')] or (value.tag == q('cvParam') and value.get('accession') in scan_cvs):
                        scan.append(copy.deepcopy(value))
ET.indent(out)
(ROOT / 'tests/data/mzml_acquisition_source_projection.mzML').write_bytes(ET.tostring(out,encoding='utf-8',xml_declaration=True) + b'\n')
