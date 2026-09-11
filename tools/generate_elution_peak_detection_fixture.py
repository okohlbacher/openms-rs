#!/usr/bin/env python3
"""Project source ElutionPeakDetection_input1.mzML to dependency-free test scalars.

Usage: generate_elution_peak_detection_fixture.py source.mzML output.tsv
"""
import base64
from pathlib import Path
import struct
import sys
import xml.etree.ElementTree as ET
import zlib


def project(source, output):
    ns = {"m": "http://psi.hupo.org/ms/mzml"}
    rows = ["# Exact source fixture scalar projection: scan\trt\tms_level\tmz\tintensity; empty scans use .\t."]
    for index, scan in enumerate(ET.parse(source).findall(".//m:spectrum", ns)):
        level = scan.find("m:cvParam[@accession='MS:1000511']", ns).attrib["value"]
        time = scan.find(".//m:cvParam[@accession='MS:1000016']", ns)
        rt = float(time.attrib["value"]) * (60 if time.attrib.get("unitAccession") == "UO:0000031" else 1)
        arrays = {}
        for array in scan.findall(".//m:binaryDataArray", ns):
            terms = {x.attrib["accession"] for x in array.findall("m:cvParam", ns)}
            data = base64.b64decode(array.find("m:binary", ns).text or "")
            if "MS:1000574" in terms:
                data = zlib.decompress(data)
            kind = "d" if "MS:1000523" in terms else "f"
            values = struct.unpack("<" + str(len(data) // struct.calcsize(kind)) + kind, data)
            arrays["mz" if "MS:1000514" in terms else "int"] = values
        assert len(arrays["mz"]) == len(arrays["int"])
        if not arrays["mz"]:
            rows.append(f"{index}\t{rt!r}\t{level}\t.\t.")
        for mz, intensity in zip(arrays["mz"], arrays["int"]):
            intensity = struct.unpack("<f", struct.pack("<f", intensity))[0]
            rows.append(f"{index}\t{rt!r}\t{level}\t{mz!r}\t{intensity!r}")
    Path(output).write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    project(sys.argv[1], sys.argv[2])
