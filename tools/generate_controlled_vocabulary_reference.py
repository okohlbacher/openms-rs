#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Independent source-loop projection; never imports or runs Rust production.

Reads byte-pinned original OBO files. Ported finite loader/index rules are
identified in CONTROLLED_VOCABULARY_SUPPORT.md. The xref_analog prefix correction
and explicit BTO decoding are the documented native corrections.
"""
import argparse
import bisect
import hashlib
import json
from pathlib import Path

FILES = [
    ("MS", "psi-ms.obo", "utf-8", "1623792d5fd37ab305bc7228ab6b51839cfde06130f9b8727bd8ce99993a27d3"),
    ("PATO", "quality.obo", "utf-8", "9c759c87833b7c966dab66bfb31bf94d7196546eb2aef84776397a0e64d245e3"),
    ("UO", "unit.obo", "utf-8", "f87734299c881fc03e7e143d35f7ab7607ac5a1f2a8cad07603e1673ad5a15df"),
    ("BTO", "brenda.obo", "cp1252", "4466c14e537ec6ae2b79d7cb22ea16ef0e79f9317d65d1f281c0003e9dbb846a"),
    ("GO", "goslim_goa.obo", "utf-8", "542711370e7406b3461902b2e0214a054237e51c8d1f2dd01f07059001cb50f8"),
]
WS = " \t\r\n"
TYPES = [
    (["string"], 0), (["integer", "int"], 1), (["decimal", "float", "double"], 2),
    (["negativeInteger"], 3), (["positiveInteger"], 4), (["nonNegativeInteger"], 5),
    (["nonPositiveInteger"], 6), (["boolean", "bool"], 7), (["date"], 8), (["anyURI"], 9),
]


def term():
    return dict(id="", name="", description="", parents=set(), children=set(),
                obsolete=False, synonyms=[], unparsed=[], xref_type=10, xref_binary=[], units=set())


def quoted(line):
    return line.split('"', 1)[-1].strip(WS).split('"', 1)[0].strip(WS)


def relationship(line, kind):
    # Source starts one byte after the keyword; fixture keywords are ASCII.
    start = line.find(kind) + len(kind) + 1
    return line[start:].split(":", 1)[0] + ":" + line.rsplit(":", 1)[-1].split("!", 1)[0].strip(WS)


def load(name, data, terms, names, headers):
    current, in_term, defined = term(), False, 0
    headers["name"] = name

    def commit():
        nonlocal defined
        if current["id"]:
            terms[current["id"]] = current
            defined += 1

    # Only LF is the source delimiter; splitlines would additionally split BTO
    # legacy controls or Unicode line separators and manufacture new records.
    for raw in data.split("\n"):
        line = raw.strip(WS)
        compact = "".join(c for c in line if c not in WS)
        if not line:
            continue
        for prefix, field in [("data-version:", "version"), ("default-namespace:", "label")]:
            if compact.startswith(prefix):
                headers[field] = line.split(":", 1)[1].strip(WS)
        if compact.startswith("remark:URL:"):
            index = line.find("http://")
            if index < 0:
                index = line.find("https://")
            if index >= 0:
                headers["url"] = line[index:].strip(WS)
        if compact.startswith("["):
            if compact.lower() == "[term]":
                commit()
                current, in_term = term(), True
            else:
                in_term = False
            continue
        if not in_term:
            continue
        if compact.startswith("id:"):
            current["id"] = line.split(":", 1)[1].strip(WS)
        elif compact.startswith("name:"):
            current["name"] = line.split(":", 1)[1].strip(WS)
        elif compact.startswith("is_a:"):
            current["parents"].add(line.split(":", 1)[1].split("!", 1)[0].strip(WS))
        elif name == "brenda" and compact.startswith(("relationship:DRV", "relationship:part_of")):
            kind = "DRV" if compact.startswith("relationship:DRV") else "part_of"
            current["parents"].add(relationship(line, kind))
        elif compact.startswith("relationship:has_units"):
            current["units"].add(relationship(line, "has_units"))
        elif compact.startswith("def:"):
            current["description"] = quoted(line)
        elif compact.startswith("synonym:"):
            current["synonyms"].append(quoted(line))
        elif compact == "is_obsolete:true":
            current["obsolete"] = True
        elif compact.startswith(("xref:value-type", "xref_analog:value-type", "relationship:has_value_type")):
            is_relation = compact.startswith("relationship:")
            expression = compact if is_relation else compact.replace("\\", "")
            for aliases, number in TYPES:
                if any(("" if is_relation else "value-type:") + "xsd:" + alias in expression for alias in aliases):
                    current["xref_type"] = number
                    break
            else:
                if is_relation and any(x in expression for x in ["MS:1002711", "MS:1002712", "MS:1002713"]):
                    current["xref_type"] = 0
        elif compact.startswith(("xref:binary-data-type", "xref_analog:binary-data-type")):
            expression = compact.replace("\\", "").split('"', 1)[0]
            current["xref_binary"].append(expression[29 if expression.startswith("xref_analog:") else 22:].strip(WS))
        else:
            current["unparsed"].append(line)
    commit()
    keys = sorted(terms)
    index = 0
    while index < len(keys):
        key = keys[index]
        record = terms[key]
        for parent in sorted(record["parents"]):
            if parent not in terms:
                terms[parent] = term()
                bisect.insort(keys, parent)
            terms[parent]["children"].add(key)
        name = record["name"]
        if name in names:
            name += record["description"]
        names.setdefault(name, key)
        # Same ++map-iterator behavior after insertions, not a frozen key list.
        index = bisect.bisect_right(keys, key)
    return defined


def encoded(value):
    if isinstance(value, set):
        value = sorted(value)
    if isinstance(value, list):
        return str(len(value)) + ":" + ",".join(x.encode().hex() for x in value)
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, int):
        return str(value)
    return value.encode().hex()


def outputs(source):
    terms, names, headers = {}, {}, dict(name="", label="", version="", url="")
    counts, byte_hashes = [], []
    for name, filename, encoding, expected in FILES:
        raw = (source / filename).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == expected, filename
        fnv = 14695981039346656037
        for byte in raw:
            fnv = ((fnv ^ byte) * 1099511628211) % 2**64
        byte_hashes.append(dict(file=filename, bytes=len(raw), sha256=expected, fnv1a64=fnv))
        count = load(name, raw.decode(encoding), terms, names, headers)
        counts.append(dict(load_name=name, definitions=count, terms=len(terms),
                           placeholders=sum(not t["id"] for t in terms.values()), names=len(names), **{k:v for k,v in headers.items() if k!="name"}))
    fields = ["id", "name", "description", "parents", "children", "obsolete", "xref_type", "synonyms", "unparsed", "xref_binary", "units"]
    lines = ["key\t" + "\t".join(fields)]
    for key in sorted(terms):
        lines.append(encoded(key) + "\t" + "\t".join(encoded(terms[key][field]) for field in fields))
    aliases = ["name\tid"] + [encoded(name)+"\t"+encoded(names[name]) for name in sorted(names)]
    summary = dict(source_revision="82ce5b373c97f934ffd9b1ffd80215ca66473d0b",
                   method="Independent Python transcription of pinned OBO/index loops, not Rust-produced expected values.",
                   corrections=["explicit Windows-1252 BTO decoding", "correct xref_analog binary prefix length"],
                   providers=counts, raw_hashes=byte_hashes)
    return {"controlled_vocabulary_projection.tsv":"\n".join(lines)+"\n",
            "controlled_vocabulary_aliases.tsv":"\n".join(aliases)+"\n",
            "controlled_vocabulary_projection.json":json.dumps(summary, indent=2)+"\n"}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=Path(__file__).resolve().parents[1]/"resources/cv")
    parser.add_argument("--check", action="store_true")
    args=parser.parse_args()
    destination=Path(__file__).resolve().parents[1]/"tests/data"
    for filename, content in outputs(args.source_dir).items():
        path=destination/filename
        if args.check:
            assert path.read_bytes() == content.encode(), filename
        else:
            path.write_bytes(content.encode())
        print(filename, len(content.encode()), hashlib.sha256(content.encode()).hexdigest())
