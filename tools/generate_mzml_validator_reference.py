#!/usr/bin/env python3
"""Extract unchanged upstream semantic-validator assertions; never invoke Rust."""
import argparse
import hashlib
import json
from pathlib import Path

SOURCE = 'src/tests/class_tests/openms/source/MzMLFile_test.cpp'
SOURCE_SHA256 = 'f8b49cc35915c8e45e5ffca9ba6525f95d36eec1f24a0c7216e127510a0924fe'
REVISION = "82ce5b373c97f934ffd9b1ffd80215ca66473d0b"


def projection(root):
    raw = (root / SOURCE).read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA256:
        raise SystemExit("source class test hash differs from the audited pin")
    lines = raw.decode("utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("START_SECTION(bool isSemanticallyValid("))
    end = next(i for i in range(start, len(lines)) if lines[i] == "END_SECTION")
    assertions = [{"line": i + 1, "source": lines[i].strip()} for i in range(start, end + 1)
                  if "TEST_EQUAL(" in lines[i]]
    calls = [a for a in assertions if "isSemanticallyValid(" in a["source"]]
    if len(calls) != 5 or len(assertions) != 15:
        raise SystemExit("source assertion surface changed")
    return {"source_revision": REVISION, "source_path": SOURCE, "source_sha256": SOURCE_SHA256,
            "lines": [start + 1, end + 1], "evidence": "Verbatim source assertions only. The two store-produced inputs are unavailable source writer outputs; they are not native writer goldens.",
            "source_text": "\n".join(lines[start:end + 1]) + "\n", "assertions": assertions,
            "cases": [
                {"input": "source store of empty experiment", "valid": True, "errors": 0, "warnings": 0, "fixture_available": False},
                {"input": "source store of loaded MzMLFile_1.mzML", "valid": True, "errors": 0, "warnings": 2, "fixture_available": False},
                {"input": "MzMLFile_1.mzML", "valid": True, "errors": 0, "warnings": 0, "fixture_available": True},
                {"input": "MzMLFile_4_indexed.mzML", "valid": True, "errors": 0, "warnings": 0, "fixture_available": True},
                {"input": "MzMLFile_3_invalid.mzML", "valid": False, "errors": 8, "warnings": 1, "fixture_available": True}]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    destination = Path(__file__).resolve().parents[1] / "tests/data/mzml_validator_source.json"
    encoded = (json.dumps(projection(args.source_root), indent=2, ensure_ascii=False) + "\n").encode()
    if args.check:
        if destination.read_bytes() != encoded:
            raise SystemExit("source projection differs")
    else:
        destination.write_bytes(encoded)


if __name__ == "__main__":
    main()
