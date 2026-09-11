#!/usr/bin/env python3
"""Extract unchanged source payloads and decode them independently using stdlib.
Run with the pinned source root as the sole argument. No native Rust outputs.
"""
from pathlib import Path
import base64, hashlib, json, math, struct, sys, xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path(sys.argv[1])
REL='src/tests/class_tests/openms/data/MSChromatogramParquetConsumer_1_output.chrom.mzML'
original=(SOURCE/REL).read_bytes()
out=ROOT/'tests/data'; out.mkdir(exist_ok=True)
(out/'mzml_numpress_source_original.mzML').write_bytes(original)
ns={'m':'http://psi.hupo.org/ms/mzml'}
tree=ET.fromstring(original)
records=tree.findall('.//m:chromatogram',ns)
def raw(array):
 import zlib
 return zlib.decompress(base64.b64decode(array.find('m:binary',ns).text))
def linear(b):
 fp=struct.unpack('>d',b[:8])[0]
 if len(b)==8:return []
 values=[int.from_bytes(b[8:12],'little')]
 if len(b)>12:values.append(int.from_bytes(b[12:16],'little'))
 nibs=[v for byte in b[16:] for v in (byte>>4,byte&15)]
 while nibs and nibs!=[0]:
  head=nibs.pop(0); leading=head if head<=8 else head-8
  number=0 if head<=8 else ((1<<(leading*4))-1)<<(32-leading*4)
  for i in range(8-leading): number|=nibs.pop(0)<<(4*i)
  if number>=2**31:number-=2**32
  values.append(2*values[-1]-values[-2]+number)
 return [v/fp for v in values]
def slof(b):
 fp=struct.unpack('>d',b[:8])[0]
 return [math.exp(int.from_bytes(b[i:i+2],'little')/fp)-1 for i in range(8,len(b),2)]
ET.register_namespace('',ns['m'])
projection=ET.Element('{'+ns['m']+'}mzML',{'version':'1.1.0'})
run=ET.SubElement(projection,'run',{'id':'projection'})
lst=ET.SubElement(run,'chromatogramList',{'count':str(len(records))})
rows=['chromatogram\tpoint\trt_decimal\trt_f64_bits\tintensity_decimal\tintensity_f32_bits']
array_count=0
for i, record in enumerate(records):
 node=ET.SubElement(lst,'chromatogram',{k:record.attrib[k] for k in ['id','index','defaultArrayLength']})
 arrays=record.find('m:binaryDataArrayList',ns);node.append(arrays)
 a=list(arrays);assert len(a)==2;array_count+=2
 times=linear(raw(a[0]));ints=slof(raw(a[1]));assert len(times)==len(ints)==int(record.attrib['defaultArrayLength'])
 for j,(rt,iv) in enumerate(zip(times,ints)):
  rows.append(f'{i}\t{j}\t{rt!r}\t{struct.pack(">d",rt).hex()}\t{iv!r}\t{struct.pack(">f",iv).hex()}')
ET.indent(projection)
(out/'mzml_numpress_source_projection.mzML').write_bytes(ET.tostring(projection,encoding='utf-8',xml_declaration=True)+b'\n')
(out/'mzml_numpress_source_values.tsv').write_text('\n'.join(rows)+'\n')
paths=[REL,'src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLHandlerHelper.h','src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp','src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp','src/openms/include/OpenMS/FORMAT/OPTIONS/PeakFileOptions.h','src/openms/source/FORMAT/OPTIONS/PeakFileOptions.cpp','src/openms/include/OpenMS/FORMAT/MSNumpressCoder.h','src/openms/source/FORMAT/MSNumpressCoder.cpp','src/openms/source/FORMAT/MSNUMPRESS/MSNumpress.cpp','src/tests/class_tests/openms/source/MSNumpressCoder_test.cpp','share/OpenMS/SCHEMAS/mzML_1_10.xsd','share/OpenMS/CV/psi-ms.obo']
# Schema filename is read from the actual pinned tree.
if not (SOURCE/paths[-2]).exists():
 paths[-2]=str(next((SOURCE/'share/OpenMS/SCHEMAS').glob('mzML1.1.0.xsd')).relative_to(SOURCE))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'source_revision':'54a232fe2cae9c590d5c997fa49d20e7769860fb','source_files':[{'path':p,'sha256':sha(SOURCE/p)} for p in paths],
 'extraction':{'chromatograms':len(records),'binary_arrays':array_count,'points':len(rows)-1,'original':'Exact source bytes retained. Projection retains original binary array element/attribute/text content and chromatogram IDs/counts, reserializing XML with UTF-8 declaration (original declares ISO-8859-1), replacing document envelope and omitting unrelated metadata/index.', 'values':'Independent Python stdlib base64/zlib, signed-nibble linear recurrence and math.exp SLOF decode, then IEEE binary32 intensity storage. Derived values are not source class-test literals or executed C++ results. RT exact; SLOF tolerate host libm last-bit differences.', 'literal_cases':'Reuses the three unchanged MSNumpressCoder class-test base64 strings at lines 202/241/281 and existing raw/zlib projection fixtures by hash.', 'native_boundaries':['strict duplicate/conflicting compression CVs','strict decoded point count','effective f64 Numpress type repairs','integer/string unsupported except source PIC integer repair','canonical float32-only writer fallback','empty compressed payload rejected','excluded records fully validated; retained-point filtering before f32 storage']},
 'fixtures':[], 'source_semantics':{'compression_CVs':'MzMLHandlerHelper.cpp:36-77,329-363', 'precision_repairs':'MzMLHandlerHelper.cpp:154-183', 'primary_writer_fallback':'MzMLHandler.cpp:5714-5750', 'auxiliary_writer_fallback':'MzMLHandler.cpp:5798-5828', 'source_type_selection_interaction':'MzMLHandler.cpp:5621-5660; fixed native ordinary precision intentionally retained'}}
for p in [out/'mzml_numpress_source_original.mzML',out/'mzml_numpress_source_projection.mzML',out/'mzml_numpress_source_values.tsv',out/'numpress_source_bytes.tsv',out/'numpress_coder_transport.tsv',out/'mzml_1_10.xsd',Path(__file__)]:
 manifest['fixtures'].append({'path':str(p.relative_to(ROOT)),'sha256':sha(p)})
(out/'mzml_numpress_provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(manifest['extraction'])
