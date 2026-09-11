#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2002-present, OpenMS Inc.
"""Project pinned source literals/routes only; no native-generated numerical oracle."""
import argparse
import hashlib
import json
from pathlib import Path
REVISION = "82ce5b373c97f934ffd9b1ffd80215ca66473d0b"
SOURCES = [{'path': 'src/openms/include/OpenMS/FORMAT/MzMLFile.h',
  'sha256': 'afea184c76be31c629b0c5f845211e99c086d10743db49ded63604ede1b56d20'},
 {'path': 'src/openms/source/FORMAT/MzMLFile.cpp',
  'sha256': '6c64cc4549b088907c4ce85b526a56ca46918728063ad9a5eaad4a5ce355900e'},
 {'path': 'src/tests/class_tests/openms/source/MzMLFile_test.cpp',
  'sha256': 'f8b49cc35915c8e45e5ffca9ba6525f95d36eec1f24a0c7216e127510a0924fe'},
 {'path': 'src/openms/include/OpenMS/FORMAT/OPTIONS/PeakFileOptions.h',
  'sha256': '108b91216a03ad3e325620dc41c9fb2cadb02553f2e6946b0e16bf8418f7909f'},
 {'path': 'src/openms/source/FORMAT/OPTIONS/PeakFileOptions.cpp',
  'sha256': '1ebf8ea1e96d63537245efce79d24a16e81033e9281b58284ab9081524d24393'},
 {'path': 'src/openms/include/OpenMS/KERNEL/MSSpectrum.h',
  'sha256': 'ce7779212eb303dc95f704daad301a1a6574c3e690c9a16b70da26a0fbd317f3'},
 {'path': 'src/openms/source/KERNEL/MSSpectrum.cpp',
  'sha256': '0a1bae71430283397b5c3a4acb358fd4edc8e7d2df0ca93655308c6784403913'},
 {'path': 'src/tests/class_tests/openms/source/MSSpectrum_test.cpp',
  'sha256': 'a0b5d1ac569910120c039e996be8e6e05829642ba00d52b9dc96ba6ec398d16e'},
 {'path': 'src/openms/include/OpenMS/FORMAT/PeakTypeEstimator.h',
  'sha256': '039634ddb6e898e0c76428723a0a88a3069a380bf09e7421d2d2ae640fea218e'},
 {'path': 'src/tests/class_tests/openms/source/PeakTypeEstimator_test.cpp',
  'sha256': '313686c99b54d2f0cb4658836cc2610563b518774afb59355b56d9c6a26c9c10'},
 {'path': 'src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLHandler.h',
  'sha256': 'c15c097047358e3e2c3208c767d523f25e95fd18f6b1aba67a8b6383374b0287'},
 {'path': 'src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp',
  'sha256': '8546ea20bbfd478e5464001e0a3e0cb931142eed250c244f4f19113159d9198a'},
 {'path': 'src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLHandlerHelper.h',
  'sha256': 'c36f957c183551e34ada80ac5c595d4352fdc473e7820ef129189887004095b0'},
 {'path': 'src/openms/source/FORMAT/HANDLERS/MzMLHandlerHelper.cpp',
  'sha256': '339b9f224fe3b142d2af42b65869b838549810d379cd406e85609ad499ae62f9'},
 {'path': 'src/openms/include/OpenMS/FORMAT/VALIDATORS/MzMLValidator.h',
  'sha256': '7da5dad2c7c950883db37db044af658fed0e29012af4a2a7205847b64f5c708f'},
 {'path': 'src/openms/source/FORMAT/VALIDATORS/MzMLValidator.cpp',
  'sha256': '4e7ef71c2fb4c486308b4b0bbf8d79fea3f1ba310d860a7f13be13535425f776'},
 {'path': 'src/tests/class_tests/openms/source/MzMLValidator_test.cpp',
  'sha256': '64824e2fac71f7a44ea9fc9b2776a1326405550696613cf92f18277136fc41a5'},
 {'path': 'src/openms/include/OpenMS/FORMAT/VALIDATORS/SemanticValidator.h',
  'sha256': '539f5dd9331d6f343257477a871ffdd5a157374a92dc12da4d65c75c058cf858'},
 {'path': 'src/openms/source/FORMAT/VALIDATORS/SemanticValidator.cpp',
  'sha256': '99d9ee2fc48a95880b935b5ffc4a35433cad25850d4e6cebc7e319358daeac59'},
 {'path': 'src/openms/include/OpenMS/FORMAT/VALIDATORS/XMLValidator.h',
  'sha256': 'd1af4acfc0836b18a4a7969ff0c4a69dfe48845595a081a516071682de49556f'},
 {'path': 'src/openms/source/FORMAT/VALIDATORS/XMLValidator.cpp',
  'sha256': '16193b35d7de8885d1061581d9e3eac717b0bfdd0aee41c6ba37a9a5f5bde7cc'},
 {'path': 'src/openms/include/OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h',
  'sha256': 'acb4f9a230851f6dd1359cfba4e188458436e2de6122b7b75c992cb808ff9dfb'},
 {'path': 'src/openms/source/FORMAT/HANDLERS/IndexedMzMLDecoder.cpp',
  'sha256': 'fc9343bf28a6756f063f1fa91471a015d2f6d09e572b6a136372fa9d0473478c'},
 {'path': 'share/OpenMS/MAPPING/ms-mapping.xml',
  'sha256': '8407daa22c49e61fb4d6c5a965899a158bba96cb8f4783dd5fbfe3b097dc18ff'},
 {'path': 'share/OpenMS/SCHEMAS/mzML_1_10.xsd',
  'sha256': '8be24aab80f1f43a610d84745534e25abcd8734738c2699b0f67b4c242c154db'},
 {'path': 'share/OpenMS/SCHEMAS/mzML_idx_1_10.xsd',
  'sha256': 'b99b668f2d9faabbf62ab8cde1f778d126bf0213d083472590670ee64318f5b0'},
 {'path': 'src/tests/class_tests/openms/data/MzMLFile_1.mzML',
  'sha256': '076fd42e8b2281b526868c874e4e267e18c0ed67543b03d2a718d86116390bc2'},
 {'path': 'src/tests/class_tests/openms/data/MzMLFile_3_invalid.mzML',
  'sha256': '3ee7f12e2d111d4c2db39c79ce16244b45d43249ff2de902de249ae7e1e3b013'},
 {'path': 'src/tests/class_tests/openms/data/MzMLFile_4_indexed.mzML',
  'sha256': '51425c17fe6bdfa0dff926b6f0e8fc65c48752dbe9fdd9acf05616eaf0d8d4cc'},
 {'path': 'src/tests/class_tests/openms/data/PeakTypeEstimator_raw.dta',
  'sha256': '262bb2a65406892bd6ed1c0027151a863d91e7653bc487402ae58122647ea792'},
 {'path': 'src/tests/class_tests/openms/data/PeakTypeEstimator_rawTOF.dta',
  'sha256': '8419dab0821b3ea794a8fa73c1846c5b0f933564f2fc5618a65756d03b70de14'},
 {'path': 'src/tests/class_tests/openms/data/PeakTypeEstimator_peak.dta',
  'sha256': '69d552ec11db33256f09c1e403aecc45018f3d4c2553ba27112d8ebecd6d5297'},
 {'path': 'src/openms/include/OpenMS/FORMAT/HANDLERS/XMLHandler.h',
  'sha256': '43d1f74e3a061475c167bfef9a5059936889f6c9bbb0a4516200f64544443264'},
 {'path': 'src/openms/source/METADATA/MetaInfoInterface.cpp',
  'sha256': '4784d93bf3e62f4774f926520f1f9c7b0ab3598e686091c5c83ecf1569368699'},
 {'path': 'src/openms/source/METADATA/MetaInfo.cpp',
  'sha256': 'bafa2699024ef964a166234e0706bc0b29d5314abd45ff516046dd568ba42e18'},
 {'path': 'src/openms/source/DATASTRUCTURES/DataValue.cpp',
  'sha256': 'ab04af951f41dc9b767f73287557d3c5a3f7b7edcacc6611fefffe3f02a2d89b'},
 {'path': 'src/openms/source/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.cpp',
  'sha256': 'a0e5ea8e0ecc64635b3f4a453e4d4e315d4087572f974dad819c590298bc1dbd'},
 {'path': 'src/openms/source/CHEMISTRY/SpectrumAnnotator.cpp',
  'sha256': '35e2add302b5c830461539b8743cad51bee285182abf64b20b494c60faa1fe3d'},
 {'path': 'share/OpenMS/CV/psi-ms.obo',
  'sha256': '1623792d5fd37ab305bc7228ab6b51839cfde06130f9b8727bd8ce99993a27d3'},
 {'path': 'share/OpenMS/CV/unit.obo',
  'sha256': 'f87734299c881fc03e7e143d35f7ab7607ac5a1f2a8cad07603e1673ad5a15df'},
 {'path': 'src/openms/include/OpenMS/METADATA/InstrumentSettings.h',
  'sha256': '6c7e3122bfe3f6fddc39102ce91c5a87053946d642353d0f6388ab7501a836ed'},
 {'path': 'src/openms/source/METADATA/InstrumentSettings.cpp',
  'sha256': 'f412905ada7edc511bb89f974be4a694163e3a247db483999a9cf28757c64632'},
 {'path': 'src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h',
  'sha256': '1ea32c73e80d7607aa699a26a22be0ed2eca2d9d0b3249b4850b466384f2a74f'}]
SECTIONS = [
    ("src/openms/source/DATASTRUCTURES/DataValue.cpp", 452, 491, "source numeric casts and CPP-058 boundary"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 1380, 1408, "elution-time fallback before array metadata merging"),
    ("src/tests/class_tests/openms/source/MzMLFile_test.cpp", 1362, 1515, "source typed/independent-noise/PDA/pressure class-test literals"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 350, 450, "spectrum noise extraction and primary metadata merge"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 620, 755, "chromatogram roles and metadata order; CPP-053"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 1550, 1580, "absorption role"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 1677, 1729, "spectrum record CV routes"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 2273, 2420, "scan record CV routes"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 3418, 3570, "typed user parameters"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 3637, 3737, "source scalar/list serialization"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 5529, 5543, "independent-noise write order"),
    ("src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp", 5611, 5720, "primary precision, roles, units and noise codec suppression"),
    ("src/openms/source/METADATA/MetaInfo.cpp", 45, 65, "right metadata overwrites left"),
]
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    for item in SOURCES:
        data=(args.source/item["path"]).read_bytes()
        if hashlib.sha256(data).hexdigest()!=item["sha256"]:
            raise SystemExit("source hash mismatch: "+item["path"])
    sections=[]
    for path,first,last,purpose in SECTIONS:
        lines=(args.source/path).read_text().splitlines(keepends=True)
        sections.append(dict(path=path,first_line=first,last_line=last,purpose=purpose,text="".join(lines[first-1:last])))
    # Preserve complete OBO stanzas for all scalar routes, role CVs and units.
    ids = "1000285 1000504 1000505 1000527 1000528 1000618 1000619 1000796 1000797 1000798 1000502 1000011 1000015 1000512 1000803 1000616 1000800 1000880 1000617 1000515 1000821 1000820 1000786 1002743 1002744 1002745".split()
    wanted = {"MS:"+x for x in ids}
    obo = (args.source/"share/OpenMS/CV/psi-ms.obo").read_text().splitlines(keepends=True)
    terms=[]
    for i,line in enumerate(obo):
        if line.startswith("id: ") and line.strip()[4:] in wanted:
            a=i-1
            b=i+1
            while b<len(obo) and not obo[b].startswith("["):b+=1
            terms.append(dict(accession=line.strip()[4:],first_line=a+1,last_line=b,text="".join(obo[a:b])))
    if {x["accession"] for x in terms}!=wanted:raise SystemExit("missing source CV terms")
    out=dict(source_revision=REVISION,evidence="Verbatim source slices; source-reviewed evidence, no executed C++ parity claim and no Rust-derived expected output.",sources=SOURCES,sections=sections,cv_terms=terms)
    content=json.dumps(out,indent=2,ensure_ascii=False)+"\n"
    path=Path(__file__).resolve().parents[1]/"tests/data/mzml_typed_transport_source.json"
    if args.check:
        if path.read_text()!=content:raise SystemExit("source projection differs")
    else:path.write_text(content)
    print(f"verified {len(SOURCES)} source hashes; {len(sections)} source sections; {sum(x['last_line']-x['first_line']+1 for x in sections)} source lines")
if __name__=="__main__":main()
