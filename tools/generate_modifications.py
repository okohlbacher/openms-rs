#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Generate the native lookup data from the bundled, pinned modification XML.

The output is the OpenMS Rust Modification Table, a transformed UniMod dataset
under the Design Science License. Original XML and notices remain alongside it.
Only declared delta/neutral-loss elements are read, never Ignore or brick nodes.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "resources" / "modifications"
NS = {"u": "http://www.unimod.org/xmlns/schema/unimod_2"}


def formula(node):
    if node is None:
        return ""
    counts = {}
    for atom in node.findall("u:element", NS):
        match = re.fullmatch(r"(\d*)([A-Z][a-z]*)", atom.attrib["symbol"])
        if not match:
            raise ValueError(f"Unsupported element {atom.attrib['symbol']}")
        isotope, symbol = match.groups()
        key = f"({isotope}){symbol}" if isotope else symbol
        counts[key] = counts.get(key, 0) + int(atom.attrib["number"])
    return "".join(f"{key}{count}" for key, count in sorted(counts.items()) if count)


def generate(check=False):
    rows = []
    sources = []
    terms = {"Anywhere": "anywhere", "Any N-term": "n-term",
             "Any C-term": "c-term", "Protein N-term": "protein-n-term",
             "Protein C-term": "protein-c-term"}
    for filename in ["unimod.xml", "custom_mods.xml"]:
        path = DATA / filename
        mods = ET.parse(path).getroot().findall("u:modifications/u:mod", NS)
        sources.append({"file": filename, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "records": len(mods)})
        for mod in mods:
            delta = mod.find("u:delta", NS)
            if delta is None:
                raise ValueError("Missing modification delta")
            for spec in mod.findall("u:specificity", NS):
                losses = []
                for loss in spec.findall("u:NeutralLoss", NS):
                    composition = formula(loss)
                    if composition:
                        losses.append("@".join([composition, loss.attrib["mono_mass"],
                                                loss.attrib["avge_mass"]]))
                row = [mod.attrib["record_id"], mod.attrib["title"], mod.attrib["full_name"],
                       spec.attrib["site"], terms[spec.attrib["position"]],
                       delta.attrib["mono_mass"], delta.attrib["avge_mass"], formula(delta),
                       spec.attrib["hidden"], spec.attrib["classification"], ";".join(losses)]
                if any("\t" in field or "\n" in field or "\r" in field for field in row):
                    raise ValueError("TSV control character in source field")
                rows.append("\t".join(row))
    header = ("# OpenMS Rust Modification Table, generated 2026-09-10\n"
              "# Derived from UniMod (C) 2002-2006 Unimod, Design Science License,\n"
              "# and OpenMS custom modifications, BSD-3-Clause. See README.md.\n"
              "# id\tname\tfull_name\tsite\tterm\tmono_delta\taverage_delta\tformula\thidden\tclassification\tneutral_losses\n")
    output = (header + "\n".join(rows) + "\n").encode()
    provenance = {"source_revision": "7c029e8cdba6abab503708ecdd56f6ab55e38ce4",
                  "source_files": sources, "specificity_records": len(rows),
                  "generated_table_sha256": hashlib.sha256(output).hexdigest(),
                  "transformation": "XML records flattened by specificity; delta and neutral-loss atom counts canonicalized; original masses retained"}
    files = {"openms-rust-modifications.tsv": output,
             "provenance.json": (json.dumps(provenance, indent=2) + "\n").encode()}
    for filename, content in files.items():
        if check:
            if (DATA / filename).read_bytes() != content:
                raise ValueError(f"Generated data is out of date: {filename}")
        else:
            (DATA / filename).write_bytes(content)
    print(f"{'Verified' if check else 'Generated'} {len(rows)} modification-specificity records")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify without changing files")
    generate(check=parser.parse_args().check)
