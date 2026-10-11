#!/usr/bin/env python3
"""Cases for MascotXMLFile::initializeLookup with a caller-supplied scan_regex
(lane F3a). Writes inputs/spectra_*.tsv (native_id, rt, precursor_mz, charge,
ms_level; one spectrum per line) and cases.tsv (name, mascotXML, spectra file
or '-', scan_regex or '<default>')."""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(HERE, "inputs")

text3 = open(os.path.join(IN, "MascotXMLFile_test_3.mascotXML")).read()
scans = sorted({int(s) for s in re.findall(r"<pep_scan_title>scan=(\d+)", text3)})
text2 = open(os.path.join(IN, "MascotXMLFile_test_2.mascotXML")).read()
mz_ints = sorted({int(s) for s in re.findall(r"<pep_scan_title>(\d+)\.", text2)})


def write(name, rows):
    with open(os.path.join(IN, name), "w") as handle:
        for row in rows:
            handle.write("\t".join(str(x) for x in row) + "\n")


# every other scan of test_3, as "scan=N"; RT N/10 s, precursor m/z N/7
write("spectra_scans.tsv",
      [(f"controllerType=0 controllerNumber=1 scan={s}", s / 10, round(s / 7, 6), 2, 2) for s in scans[::2]])
# native IDs that are the integer part of test_2's precursor m/z, every other one
write("spectra_ids.tsv", [(str(m), m * 3 + 0.25, m + 0.5, 3, 2) for m in mz_ints[::2]])
# 4000 spectra by position, for INDEX0/INDEX1 look-ups of test_3's scan numbers
write("spectra_index.tsv", [(f"index={i}", 1000 + i, 500 + i / 8, 2, 2) for i in range(4000)])

cases = [
    ("r01_default_no_spectra", "MascotXMLFile_test_2.mascotXML", "-", "<default>"),
    ("r02_rt_suffix", "MascotXMLFile_test_2.mascotXML", "-", r"_(?<RT>\d+\.\d+)$"),
    ("r03_mz_and_integer_rt", "MascotXMLFile_test_2.mascotXML", "-", r"^(?<MZ>\d+\.\d+)_(?<RT>\d+)"),
    ("r04_scan_lookup", "MascotXMLFile_test_3.mascotXML", "spectra_scans.tsv", r"scan=(?<SCAN>\d+)"),
    ("r05_default_scan_lookup", "MascotXMLFile_test_3.mascotXML", "spectra_scans.tsv", "<default>"),
    ("r06_index0", "MascotXMLFile_test_3.mascotXML", "spectra_index.tsv", r"=(?<INDEX0>\d+)"),
    ("r07_index1", "MascotXMLFile_test_3.mascotXML", "spectra_index.tsv", r"=(?<INDEX1>\d+)"),
    ("r08_native_id", "MascotXMLFile_test_2.mascotXML", "spectra_ids.tsv", r"^(?<ID>\d+)\."),
    ("r09_rt_lookup", "MascotXMLFile_test_2.mascotXML", "spectra_ids.tsv", r"_(?<RT>\d+)\."),
    ("r10_no_named_group", "MascotXMLFile_test_2.mascotXML", "-", r"scan=(\d+)"),
    ("r11_whitespace", "MascotXMLFile_test_2.mascotXML", "-", " "),
    ("r12_invalid_expression", "MascotXMLFile_test_2.mascotXML", "-", r"(?<SCAN>\d+"),
    ("r13_unmatched", "MascotXMLFile_test_3.mascotXML", "spectra_scans.tsv", r"Spectrum(?<SCAN>\d+)"),
    ("r14_scan_no_spectra", "MascotXMLFile_test_3.mascotXML", "-", r"scan=(?<SCAN>\d+)"),
]
write("cases.tsv", cases)
print(len(cases), "cases", len(scans), "scans", len(mz_ints), "mz ints")
