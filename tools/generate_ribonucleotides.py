#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Regenerate the pinned RNA record projection; no compiler or network needed.

This generator implements only the grammar present in the bundled inputs. The
library's bounded TSV/JSON readers handle arbitrary caller inputs. MODOMICS data
has separate notices in resources/rna/README.md; this software license does not
relicense the source dataset or generated record values.
"""

import argparse
import json
import re
import struct
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "resources/rna"


def element_masses():
    """Use the port's source-audited isotope literals, with the same sum order."""
    text = (ROOT / "src/chemistry/elements.rs").read_text()
    result = {}
    for symbol, body in re.findall(r'symbol: "([A-Za-z]+)".*?isotopes: &\[(.*?)\] }', text, re.S):
        isotopes = [(float(m), float(a)) for m, a in re.findall(
            r"mass: ([\deE.+-]+), abundance: ([\deE.+-]+)", body
        )]
        mono = max(isotopes, key=lambda pair: pair[1])[0]
        average = sum(mass * abundance for mass, abundance in isotopes)
        result[symbol] = (mono, average)
    assert len(result) == 84
    return result


def mass(formula, elements, average=False):
    # All bundled formulas are uncharged, natural elements and nonnegative
    # counts. Fail loudly if a future source adds a different grammar.
    parts = re.findall(r"([A-Z][a-z]*)(\d*)", formula)
    assert "".join(symbol + count for symbol, count in parts) == formula, formula
    counts = {}
    for symbol, count in parts:
        counts[symbol] = counts.get(symbol, 0) + int(count or "1")
    result = 0.0
    for symbol, count in sorted(counts.items()):
        result += elements[symbol][int(average)] * count
    return result


def json_entries(elements):
    document = json.loads((DATA / "Modomics.json").read_text())
    result = []
    assert len(document) == 424
    for key in sorted(document):
        row = document[key]
        code = row["short_name"]
        if not code:
            continue
        formula = row["formula"]
        moieties = row["reference_moiety"]
        term = 0
        if len(moieties) == 1 and len(moieties[0]) == 1:
            origin = moieties[0]
        else:
            assert len(moieties) == 4
            origin = "X"
            if code.endswith("pN"):
                term = 1
            elif code.startswith("N") and code.endswith("p"):
                term = 2
        mono = row.get("mass_monoiso")
        average = row.get("mass_avg")
        base = row.get("baseloss_formula")
        if base is None:
            if code.startswith("d"):
                base = "C5H10O4"
            elif code.endswith(("m", "m*")):
                base = "C6H12O5"
            elif code.endswith(("Ar(p)", "Gr(p)")):
                base = "C10H19O21P"
            else:
                base = "C5H10O5"
        alternatives = row["alternatives"][:2] if code.endswith(("?", "?*")) else None
        result.append(dict(
            name=row["name"], code=code, new_code=code,
            html_code=row.get("abbrev", "."), formula=formula, origin=origin,
            mono_mass=mass(formula, elements) if mono is None else float(mono),
            average_mass=0.0 if average is None else float(average),
            term_specificity=term, baseloss_formula=base, alternatives=alternatives,
        ))
    assert len(result) == 333
    return result


def custom_entries(elements):
    lines = (DATA / "Custom_RNA_modifications.tsv").read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if not line.startswith("#"))
    assert lines[header].startswith("name\tshort_name\tnew_nomenclature\t")
    result = []
    for line in lines[header + 1:]:
        parts = line.replace("\u2032", "'").split("\t")
        assert len(parts) >= 9, line
        name, original_code, new_code, original_origin, _, html, formula, mono, average = parts[:9]
        code = original_code[:-4] if original_code.endswith("QtRNA") else original_code
        origin = "G" if original_origin == "preQ0base" else original_origin if len(original_origin) == 1 else "."
        formula = "" if formula == "-" else formula
        declared = []
        for value, is_average in [(mono, False), (average, True)]:
            if value in {"", "None"}:
                declared.append(0.0)
            else:
                value = float(value)
                declared.append(mass(formula, elements, is_average) if value == 0 and formula else value)
        term, base, alternatives = 0, "C5H10O5", None
        if new_code.endswith("N"):
            if "55" in new_code or new_code == "N":
                term = 1
            elif "33" in new_code:
                term = 2
        elif original_code.startswith("d"):
            base = "C5H10O4"
        elif original_code.endswith(("m", "m*")):
            base = "C6H12O5"
        elif original_code.endswith("?"):
            assert len(parts) >= 10 and " " in parts[9]
            alternatives = [parts[9].split(" ", 1)[0], parts[9].rsplit(" ", 1)[1]]
        elif original_code in {"Ar(p)", "Gr(p)"}:
            base = "C10H19O21P"
        result.append(dict(
            name=name, code=code, new_code=new_code, html_code=html, formula=formula,
            origin=origin, mono_mass=declared[0], average_mass=declared[1],
            term_specificity=term, baseloss_formula=base, alternatives=alternatives,
        ))
    assert len(result) == 45
    return result


def rust_string(value):
    assert isinstance(value, str)
    escapes = {"\\": "\\\\", '"': '\\"', "\n": "\\n", "\r": "\\r", "\t": "\\t"}
    return '"' + "".join(escapes.get(c, "\\u{" + format(ord(c), "x") + "}" if ord(c) < 32 else c) for c in value) + '"'


def render(entries):
    text = [
        "// Generated by tools/generate_ribonucleotides.py; do not edit by hand.",
        "// Scientific data provenance and terms: resources/rna/README.md.",
        "// $Maintainer: OpenMS Rust contributors $",
        "&[",
    ]
    for entry in entries:
        fields = []
        for key in ["name", "code", "new_code", "html_code", "formula"]:
            fields.append(f"{key}: {rust_string(entry[key])}")
        origin = entry["origin"]
        assert len(origin) == 1 and origin.isascii() and origin not in "'\\"
        fields.append(f"origin: '{origin}'")
        for key in ["mono_mass", "average_mass"]:
            bits = struct.unpack(">Q", struct.pack(">d", entry[key]))[0]
            fields.append(f"{key}_bits: 0x{bits:016x}")
        fields.append(f'term_specificity: {entry["term_specificity"]}')
        fields.append(f'baseloss_formula: {rust_string(entry["baseloss_formula"])}')
        alternatives = entry["alternatives"]
        if alternatives is not None:
            assert len(alternatives) == 2
            alternatives = "Some([" + ", ".join(map(rust_string, alternatives)) + "])"
        else:
            alternatives = "None"
        fields.append("alternatives: " + alternatives)
        text.append("    EmbeddedRibonucleotide { " + ", ".join(fields) + " },")
    text.append("]")
    return "\n".join(text) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify without writing")
    args = parser.parse_args()
    elements = element_masses()
    entries = json_entries(elements) + custom_entries(elements)
    assert len(entries) == 378 and len({e["code"] for e in entries}) == 375
    codes = {entry["code"] for entry in entries}
    assert all(set(e["alternatives"] or []) <= codes for e in entries)
    output = DATA / "ribonucleotides.rs"
    expected = render(entries)
    if args.check:
        if not output.exists() or output.read_text() != expected:
            raise SystemExit("RNA projection differs; run generator without --check")
    else:
        output.write_text(expected)
    print("378 RNA records / 375 codes; projection " + ("verified" if args.check else "written"))


if __name__ == "__main__":
    main()
