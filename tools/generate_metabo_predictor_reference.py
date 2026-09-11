#!/usr/bin/env python3
"""Generate deterministic binary64 cases and query an explicit LIBSVM probe.

No Rust code supplies expected values. The small Python evaluator only locates
boundary neighborhoods; every recorded decision and label comes from LIBSVM.
See docs/METABO_PREDICTOR_SUPPORT.md for the reproducible compile command.
"""
import argparse
from pathlib import Path
import math
import random
import struct
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "resources/metabolite_isotope_models"


def bits(value):
    return f"{struct.unpack('<Q', struct.pack('<d', value))[0]:016x}"


def load(noise):
    lines = (DATA / f"MetaboliteIsoModelNoised{noise}.svm").read_text().splitlines()
    start = lines.index("SV")
    header = dict(line.split(" ", 1) for line in lines[:start])
    support = [(float(row[0]), [float(field.split(":")[1]) for field in row[1:]])
               for row in (line.split() for line in lines[start + 1:])]
    scale = [list(map(float, line.split())) for line in
             (DATA / f"MetaboliteIsoModelNoised{noise}.scale").read_text().splitlines()]
    return float(header["gamma"]), float(header["rho"]), support, scale


def boundary_score(model, raw):
    gamma, rho, support, scales = model
    features = [(raw[i] - scales[i][0]) / scales[i][1] for i in range(4)]
    return sum(coefficient * math.exp(-gamma * sum((a - b) ** 2 for a, b in zip(features, sv)))
               for coefficient, sv in support) - rho


def cases(noise):
    model = load(noise)
    rng = random.Random(4711 + noise)
    huge = float.fromhex('0x1.fffffffffffffp+1023')
    values = [[0., 0., 0., 0.], [-0., -0., -0., -0.], [1000., 0., 0., 0.],
              [2000., 0., 0., 0.], [-1., -1., -1., -1.], [huge] * 4,
              [-huge] * 4, [1e-300] * 4, [1., 1e200, -1e200, 0.],
              [row[0] for row in model[3]]]
    for mass in (1., 100., 300., 600., 900., 1000., 1500.):
        for ratio in (0., .05, .2, .4, .7, 1., 2.):
            values.append([mass, ratio, ratio * .6, ratio * .25])
    values += [[rng.uniform(-100., 1300.)] + [rng.uniform(-.1, 2.) for _ in range(3)] for _ in range(128)]
    for i in range(33):
        sv = model[2][i * (len(model[2]) - 1) // 32][1]
        values.append([sv[j] * model[3][j][1] + model[3][j][0] for j in range(4)])
    ordinary = values[9:]
    positive = [row for row in ordinary if boundary_score(model, row) > .01]
    negative = [row for row in ordinary if boundary_score(model, row) < -.01]
    assert positive and negative
    # Independent bisection only chooses inputs; margins stay away from the
    # platform-dependent libm last-bit region. LIBSVM supplies the final oracle.
    for i in range(12):
        left, right = positive[i % len(positive)], negative[i % len(negative)]
        low, high = 0., 1.
        for _ in range(45):
            mid = (low + high) / 2
            row = [a + mid * (b - a) for a, b in zip(left, right)]
            if boundary_score(model, row) > 0:
                low = mid
            else:
                high = mid
        for fraction in (max(0., low - 1e-6), min(1., high + 1e-6)):
            values.append([a + fraction * (b - a) for a, b in zip(left, right)])
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    inputs = [str(noise) + " " + " ".join(map(bits, row))
              for noise in (2, 5) for row in cases(noise)]
    result = subprocess.run([str(args.probe.resolve()), str(DATA)],
                            input="\n".join(inputs) + "\n", text=True,
                            capture_output=True, check=True)
    rows = result.stdout.splitlines()
    assert len(rows) == len(inputs)
    for original, row in zip(inputs, rows):
        assert row.split()[:5] == original.split()
    content = "# Executed LIBSVM v337; model raw[4] decision binary64 bits; label.\n" + result.stdout
    output = ROOT / "tests/data/metabo_predictor_libsvm.tsv"
    if args.check:
        assert output.read_text() == content, "LIBSVM reference differs"
    else:
        output.write_text(content)
    print(f"{'checked' if args.check else 'generated'} {len(rows)} executed LIBSVM cases")


if __name__ == "__main__":
    main()
