#!/usr/bin/env python3
"""Project source settings without deriving numeric expectations through Rust."""
from pathlib import Path
import copy
import xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[1]
NS = 'http://psi.hupo.org/ms/mzml'
ET.register_namespace('', NS)
def q(s): return '{' + NS + '}' + s
source = ET.parse(ROOT / 'tests/data/mzml_load_source_original.mzML').getroot()
out = ET.Element(q('mzML'), {'version': '1.1.0'})
out.append(copy.deepcopy(source.find(q('referenceableParamGroupList'))))
run = ET.SubElement(out, q('run'), {'id': 'settings_projection'})
spectra = ET.SubElement(run, q('spectrumList'), {'count': '4'})
for old in source.findall('.//' + q('spectrum')):
    item = ET.SubElement(spectra, q('spectrum'), {'id': old.attrib['id'], 'defaultArrayLength': '0'})
    for child in old:
        name = child.tag.removeprefix('{' + NS + '}')
        if name == 'referenceableParamGroupRef' or (name == 'cvParam' and child.get('accession') in ['MS:1000579','MS:1000580','MS:1000127','MS:1000128','MS:1000511','MS:1000130']):
            item.append(copy.deepcopy(child))
        elif name == 'productList':
            products = copy.deepcopy(child)
            # Pinned source declares one but contains two, also asserted by its
            # class test. Repair only the structural count for the strict reader.
            products.set('count', str(len(products.findall(q('product')))))
            item.append(products)
        elif name == 'scanList':
            scans = ET.SubElement(item, q('scanList'), {'count': child.attrib['count']})
            for scan in child.findall(q('scan')):
                new = ET.SubElement(scans, q('scan'))
                for value in scan:
                    if value.tag == q('scanWindowList') or (value.tag == q('cvParam') and value.get('accession') in ['MS:1000016','MS:1000497']): new.append(copy.deepcopy(value))
ET.indent(out)
path = ROOT / 'tests/data/mzml_settings_source_projection.mzML'
path.write_bytes(ET.tostring(out, encoding='utf-8', xml_declaration=True) + b'\n')
