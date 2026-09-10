#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: OpenMS Rust contributors $
"""Print independent decoy RNG references as JSON; uses only Python integers.

The append-only MT19937-64 recurrence differs from the Rust engine's in-place
twist. Parameters follow MathFunctions.h's named engine and supplemental
Boost.Random 1.90. No C++, Boost binary, Rust code or third-party module is used.
The standard's 10,000th-word literal is distinguished from derived references.
"""

import json

MODULUS = 1 << 64
INDICES = [
    0, 1, 2, 3, 154, 155, 156, 157, 310, 311, 312, 313,
    467, 468, 623, 624, 935, 936, 999, 1247,
]


def words(seed, count):
    state = [seed]
    for index in range(1, 312):
        previous = state[-1]
        state.append(
            (6364136223846793005 * (previous ^ (previous >> 62)) + index)
            % MODULUS
        )
    result = []
    # Initial low 31 bits of state[0] are replaced before use; normalization
    # is irrelevant to this append-only output oracle.
    for position in range(312, 312 + count):
        combined = (
            state[position - 312] // (1 << 31) * (1 << 31)
            + state[position - 311] % (1 << 31)
        )
        value = (
            state[position - 156]
            ^ (combined // 2)
            ^ (0xB5026F5AA96619E9 if combined % 2 else 0)
        )
        state.append(value)
        value ^= (value >> 29) & 0x5555555555555555
        value ^= (value << 17) & 0x71D67FFFEDA60000
        value ^= (value << 37) & 0xFFF7EEE000000000
        value ^= value >> 43
        result.append(value)
    return result


def checksum(values):
    value = 0xCBF29CE484222325
    for word in values:
        for byte in word.to_bytes(8, "little"):
            value = ((value ^ byte) * 0x100000001B3) % MODULUS
    return value


def shuffle(text, values):
    letters = list(text)
    used = 0
    for end in range(len(letters) - 1, 0, -1):
        bucket = MODULUS // (end + 1)
        while True:
            selected = values[used] // bucket
            used += 1
            if selected <= end:
                break
        letters[end], letters[selected] = letters[selected], letters[end]
    return "".join(letters), used


def reference():
    records = []
    for seed in [0, 4711, MODULUS - 1]:
        output = words(seed, 1248)
        permutation, draws = shuffle("ABCDEFGHIJKLMNOPQRST", output)
        records.append({
            "seed": seed,
            "indices": INDICES,
            "values": [output[index] for index in INDICES],
            "fnv1a64_little_endian_1248": checksum(output),
            "shuffle_input": "ABCDEFGHIJKLMNOPQRST",
            "shuffle_output": permutation,
            "shuffle_draws": draws,
        })
    assert words(5489, 10000)[-1] == 9981545732273789042
    return {
        "method": "Derived Python integer append-only recurrence and bucket division; no C++ execution.",
        "primary_source_literal": {
            "url": "https://eel.is/c++draft/rand.predef#4",
            "seed": 5489,
            "one_based_draw": 10000,
            "value": 9981545732273789042,
        },
        "records": records,
    }


if __name__ == "__main__":
    print(json.dumps(reference(), indent=2))
