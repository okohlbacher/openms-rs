"""Project source arrays/scalars without using the native reader or chemistry.
Usage: python3 tests/data/mzml_load_projection.py /path/to/source-checkout
"""
import base64, copy, hashlib, json, pathlib, re, struct, sys
import xml.etree.ElementTree as ET
source = pathlib.Path(sys.argv[1])
out = pathlib.Path(__file__).parent
rel = 'src/tests/class_tests/openms/data/MzMLFile_1.mzML'
root = ET.parse(source / rel).getroot()
for node in root.iter():
    node.tag = node.tag.split('}')[-1]
new = ET.Element('mzML', {'xmlns':'http://psi.hupo.org/ms/mzml','version':'1.1.0'})
run = ET.SubElement(new, 'run', {'id':'source-filter-projection'})
rows = ['kind\trecord\tpeak\tposition_bits\tintensity_bits']
for kind in ['spectrum', 'chromatogram']:
    originals = list(root.iter(kind))
    parent = ET.SubElement(run, kind+'List', {'count':str(len(originals))})
    for idx, old in enumerate(originals):
        record = ET.SubElement(parent, kind, {k:old.attrib[k] for k in ['id','index','defaultArrayLength']})
        if kind == 'spectrum':
            for cv in old.findall('cvParam'):
                if cv.get('accession') == 'MS:1000511': record.append(copy.deepcopy(cv))
            rt = next((cv for cv in old.iter('cvParam') if cv.get('accession') == 'MS:1000016'), None)
            if rt is not None:
                scan = ET.SubElement(ET.SubElement(record, 'scanList', {'count':'1'}), 'scan')
                scan.append(copy.deepcopy(rt))
        arrays = old.find('binaryDataArrayList')
        if arrays is None: continue
        arrays = copy.deepcopy(arrays)
        arrays.set('count', str(len(arrays)))
        decoded = {}
        for array in arrays:
            array.attrib.pop("dataProcessingRef", None)
            for user in array.findall('userParam'): array.remove(user)
            cvs = {cv.get('accession'):cv for cv in array.findall('cvParam')}
            binary = base64.b64decode(array.find('binary').text or '')
            fmt = 'd' if 'MS:1000523' in cvs else 'f'
            values = struct.unpack('<'+fmt*(len(binary)//struct.calcsize(fmt)), binary)
            for cv in ['MS:1000514','MS:1000595','MS:1000515']:
                if cv in cvs: decoded[cv] = values
        record.append(arrays)
        positions = decoded['MS:1000514' if kind=='spectrum' else 'MS:1000595']
        for i,(p,v) in enumerate(zip(positions,decoded['MS:1000515'])):
            bits = lambda x: struct.pack('>d',x).hex()
            rows.append(f'{kind}\t{idx}\t{i}\t{bits(p)}\t{bits(v)}')
ET.indent(new)
xml = ET.tostring(new,encoding='utf-8',xml_declaration=True)+b'\n'
(out/'mzml_load_source_projection.mzML').write_bytes(xml)
(out/'mzml_load_source_values.tsv').write_text('\n'.join(rows)+'\n')
# Independent current-vocabulary descendant/type projection.
terms = {}
for block in (source / 'share/OpenMS/CV/psi-ms.obo').read_text().split('[Term]')[1:]:
    values = {}; parents = []
    for line in block.splitlines():
        if line.startswith('id: '): values['id'] = line[4:]
        elif line.startswith('name: '): values['name'] = line[6:]
        elif line.startswith('is_a: '): parents.append(line[6:].split()[0])
    if 'id' in values: terms[values['id']] = (values['name'], parents, block)
def descendant(key):
    return any(p == 'MS:1000513' or descendant(p) for p in terms.get(key, ('', [], ''))[1])
canonical = ['accession\tname\ttype_mask']
for accession, (name, _, block) in terms.items():
    if descendant(accession) and accession not in ['MS:1000514', 'MS:1000515', 'MS:1000595', 'MS:1000786']:
        types = re.findall(r'binary-data-type:MS\\:([0-9]+)', block)
        mask = sum({'1000519': 1, '1000521': 2, '1000522': 4, '1000523': 8, '1001479': 16}[t] for t in types)
        canonical.append(f'{accession}\t{name}\t{mask}')
(out/'mzml_load_canonical_arrays.tsv').write_text('\n'.join(canonical)+'\n')
(out/'mzml_load_source_original.mzML').write_bytes((source/rel).read_bytes())
paths = [rel, 'src/tests/class_tests/openms/source/MzMLFile_test.cpp', 'src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp', 'src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLHandler.h', 'src/openms/source/FORMAT/OPTIONS/PeakFileOptions.cpp', 'src/openms/include/OpenMS/FORMAT/OPTIONS/PeakFileOptions.h', 'src/openms/include/OpenMS/DATASTRUCTURES/DRange.h', 'share/OpenMS/CV/psi-ms.obo', 'src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp']
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
manifest = {'source_revision':'54a232fe2cae9c590d5c997fa49d20e7769860fb','source_repository':'https://github.com/okohlbacher/OpenMS4-core','source_files':[{'path':p,'sha256':sha(source/p)} for p in paths], 'fixtures':[{'path':'tests/data/'+p,'sha256':sha(out/p)} for p in ['mzml_load_source_projection.mzML','mzml_load_source_values.tsv','mzml_load_projection.py','mzml_load_canonical_arrays.tsv','mzml_load_source_original.mzML']], 'canonical_auxiliary_terms':len(canonical)-1, 'primary_peak_rows':len(rows)-1,'source_literals':{'ms_levels':'MzMLFile_test.cpp:895-904: [1] yields RT 5.1,5.3,5.4','rt':'911-918: [5.15,5.35) yields RT 5.2,5.3','mz':'923-940: [6.5,9.5) yields counts 3,1,3,0 and m/z [7,8,9],[8],[7,8,9],[]','intensity':'945-963: [6.5,9.5) yields counts 3,1,3,0 and values [9,8,7],[8],[9,8,7],[]'},'transformation':['Keep source record order, IDs, declared primary lengths, MS-level and first scan RT CVs.','Copy all original binary array payloads and CVs verbatim in value; remove unsupported array userParams and dataProcessingRef.','Strip unrelated original header, instrument, scan, precursor, product and record metadata; regenerate minimal run wrapper.','Correct original spectrum 1 binaryDataArrayList count from 2 to its actual 4. Native strict length/count policy otherwise unchanged.','Decode original uncompressed float arrays using Python struct/base64 only; primary values stored as binary64 bits. No C++ execution and no production-derived expected values.'], 'native_boundaries':['Excluded records still decode and validate; filters do not waive syntax, finite-value, type, array or raw resource checks.','Only this explicit scientific projection is claimed readable; original rich MzMLFile_1 compatibility remains separate.','Stable equal-coordinate sorting and skip_spectra are native additions.', 'The 26 canonical non-primary names are reserved roles and serialize to their pinned CV accessions. Nonstandard input names equal to a reserved canonical name are rejected; aliases and case variants are not inferred.', 'Canonical declared binary type constraints are checked; charge-array output uses signed i32. Missing CV type restrictions remain permissive. Auxiliary units are unsupported.','Metadata-only, fill_data=false, skip_xml_checks and isolation-target precursor output are rejected.']}
(out/'mzml_load_provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(len(rows)-1)
