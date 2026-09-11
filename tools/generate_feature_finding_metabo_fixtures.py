#!/usr/bin/env python3
"""Project authoritative FFM mzML inputs and featureXML expected data, without Rust.
Usage: generate_feature_finding_metabo_fixtures.py PINNED_SDK OUTPUT_DIRECTORY
"""
import base64
from pathlib import Path
import struct
import sys
import xml.etree.ElementTree as ET
import zlib


def project_input(source, destination):
    ns = {"m": "http://psi.hupo.org/ms/mzml"}
    rows = ["# scan\trt\tms_level\tmz\tintensity\tIM name or .\tIM value or ."]
    for index, scan in enumerate(ET.parse(source).findall(".//m:spectrum", ns)):
        level = scan.find("m:cvParam[@accession='MS:1000511']", ns).attrib["value"]
        time = scan.find(".//m:cvParam[@accession='MS:1000016']", ns)
        rt = float(time.attrib["value"]) * (60 if time.attrib.get("unitAccession") == "UO:0000031" else 1)
        arrays = {}
        name = "."
        for array in scan.findall(".//m:binaryDataArray", ns):
            terms = {x.attrib["accession"]: x.attrib for x in array.findall("m:cvParam", ns)}
            raw = base64.b64decode(array.find("m:binary", ns).text or "")
            if "MS:1000574" in terms:
                raw = zlib.decompress(raw)
            kind = "d" if "MS:1000523" in terms else "f"
            values = struct.unpack("<" + str(len(raw) // struct.calcsize(kind)) + kind, raw)
            if "MS:1000514" in terms:
                key = "mz"
            elif "MS:1000515" in terms:
                key = "intensity"
            else:
                key = "im"
                assert terms["MS:1000786"]["value"] == "Ion Mobility"
                name = "Ion Mobility"
            arrays[key] = values
        assert len(arrays["mz"]) == len(arrays["intensity"])
        assert name == "." or len(arrays["im"]) == len(arrays["mz"])
        if not arrays["mz"]:
            rows.append(f"{index}\t{rt!r}\t{level}\t.\t.\t{name}\t.")
        for i, (mz, intensity) in enumerate(zip(arrays["mz"], arrays["intensity"])):
            intensity = struct.unpack("<f", struct.pack("<f", intensity))[0]
            im = "." if name == "." else repr(struct.unpack("<f", struct.pack("<f", arrays["im"][i]))[0])
            rows.append(f"{index}\t{rt!r}\t{level}\t{mz!r}\t{intensity!r}\t{name}\t{im}")
    destination.write_text("\n".join(rows) + "\n")


def project_output(source, destination):
    rows = ["# F index rt mz intensity quality charge; M index name type value; H index hull x y"]
    for i, f in enumerate(ET.parse(source).findall(".//feature")):
        rt, mz = [f.find(f"position[@dim='{d}']").text for d in [0, 1]]
        rows.append("\t".join(["F", str(i), rt, mz, f.findtext("intensity"), f.findtext("overallquality"), f.findtext("charge")]))
        for meta in f.findall("UserParam"):
            rows.append("\t".join(["M", str(i), meta.attrib["name"], meta.attrib["type"], meta.attrib["value"]]))
        for h, hull in enumerate(f.findall("convexhull")):
            for point in hull.findall("pt"):
                rows.append("\t".join(["H", str(i), str(h), point.attrib["x"], point.attrib["y"]]))
    destination.write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    sdk, out = map(Path, sys.argv[1:])
    data = sdk / "src/tests/class_tests/openms/data"
    project_input(data / "FeatureFindingMetabo_input1.mzML", out / "feature_finding_metabo_source.tsv")
    project_input(data / "extracted_topp/FeatureFinderMetabo_6_input.mzML", out / "feature_finding_metabo_im_source.tsv")
    project_output(data / "FeatureFindingMetabo_output1.featureXML", out / "feature_finding_metabo_expected.tsv")
