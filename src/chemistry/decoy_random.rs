// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// Copyright Jens Maurer 2000-2001
// Copyright Steven Watanabe 2010, 2011
// SPDX-License-Identifier: BSD-3-Clause AND BSL-1.0
// $Maintainer: OpenMS Rust contributors $
//
// MT19937-64 normalization and bounded-integer mapping derive from Boost.Random
// 1.90 mersenne_twister.hpp and uniform_int_distribution.hpp, distributed under
// the Boost Software License, Version 1.0 (https://www.boost.org/LICENSE_1_0.txt).
//! Private deterministic engine for OpenMS Math::RandomShuffler semantics.
//! The outer decoy generator stages this small state and its output so failed
//! draw/work checks remain atomic at the public operation boundary.

use crate::{Error, Result};

const WORDS: usize = 312;
const MIDDLE: usize = 156;
const LOWER: u64 = 0x7fff_ffff;
const UPPER: u64 = !LOWER;
const MATRIX: u64 = 0xb502_6f5a_a966_19e9;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct DecoyRandom {
    state: [u64; WORDS],
    index: usize,
}

impl DecoyRandom {
    pub(super) fn seeded(seed: u64) -> Self {
        let mut random = Self {
            state: [0; WORDS],
            index: WORDS,
        };
        random.reseed(seed);
        random
    }

    pub(super) fn reseed(&mut self, seed: u64) {
        self.state[0] = seed;
        for i in 1..WORDS {
            let previous = self.state[i - 1];
            self.state[i] = (previous ^ (previous >> 62))
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(i as u64);
        }
        self.index = WORDS;
        self.normalize();
    }

    /// Descending Fisher–Yates, as in MathFunctions.h. Every raw draw, including
    /// rejection, consumes one allowance before advancing the engine. Empty and
    /// singleton slices use no draws. On error the private slice/state may have
    /// advanced through earlier swaps; the public caller stages both.
    pub(super) fn shuffle(&mut self, data: &mut [u8], remaining_draws: &mut usize) -> Result<()> {
        for end in (1..data.len()).rev() {
            let high = u64::try_from(end)
                .map_err(|_| Error::InvalidValue("decoy shuffle range exceeds u64".into()))?;
            let selected = self.bounded(high, remaining_draws)? as usize;
            // selected <= end, so the native usize conversion cannot truncate.
            data.swap(end, selected);
        }
        Ok(())
    }

    // Boost uniform_int(0, high) uses bucket division/rejection, not modulo.
    // u128 expresses 2^64 directly; for high>=1 the bucket fits a nonzero u64.
    fn bounded(&mut self, high: u64, remaining_draws: &mut usize) -> Result<u64> {
        if high == 0 {
            return Ok(0);
        }
        let bucket = ((1_u128 << 64) / (u128::from(high) + 1)) as u64;
        loop {
            *remaining_draws = remaining_draws
                .checked_sub(1)
                .ok_or_else(|| Error::InvalidValue("decoy random draw limit exceeded".into()))?;
            let result = self.next_u64() / bucket;
            if result <= high {
                return Ok(result);
            }
        }
    }

    fn next_u64(&mut self) -> u64 {
        if self.index == WORDS {
            self.twist();
        }
        let mut value = self.state[self.index];
        self.index += 1;
        value ^= (value >> 29) & 0x5555_5555_5555_5555;
        value ^= (value << 17) & 0x71d6_7fff_eda6_0000;
        value ^= (value << 37) & 0xfff7_eee0_0000_0000;
        value ^= value >> 43;
        value
    }

    fn twist(&mut self) {
        // Intentionally in-place: the second half reads words already replaced
        // in this cycle, and the last word uses the updated first word's low bits.
        for i in 0..WORDS {
            let combined = (self.state[i] & UPPER) | (self.state[(i + 1) % WORDS] & LOWER);
            self.state[i] = self.state[(i + MIDDLE) % WORDS]
                ^ (combined >> 1)
                ^ if combined & 1 == 0 { 0 } else { MATRIX };
        }
        self.index = 0;
    }

    fn normalize(&mut self) {
        // The low 31 bits of the first seeded word are redundant. Boost canonicalizes
        // them using the inverse recurrence, then repairs an all-zero state.
        let mut value = self.state[MIDDLE - 1] ^ self.state[WORDS - 1];
        value = if value & (1_u64 << 63) == 0 {
            value << 1
        } else {
            ((value ^ MATRIX) << 1) | 1
        };
        self.state[0] = (self.state[0] & UPPER) | (value & LOWER);
        if self.state.iter().all(|&word| word == 0) {
            self.state[0] = 1_u64 << 63;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Derived reference words from an append-only integer recurrence, rather
    // than this engine's in-place cyclic twist. Reproduce them with
    // tools/generate_decoy_rng_reference.py; no C++ outputs were executed.
    // The aggregate FNV check covers all 1248 words, not only these checkpoints.
    const INDICES: [usize; 20] = [
        0, 1, 2, 3, 154, 155, 156, 157, 310, 311, 312, 313, 467, 468, 623, 624, 935, 936, 999, 1247,
    ];
    const REFERENCES: [(u64, u64, [u64; 20]); 3] = [
        (
            0,
            16540305979319870945,
            [
                2947667278772165694,
                18301848765998365067,
                729919693006235833,
                11021831128136023278,
                8347057230978218896,
                10994878445802041726,
                18081729901143002129,
                2316049702218351750,
                14165040369326021580,
                11228354904504431959,
                17661967264253682746,
                18066610840454826573,
                7563074798427102612,
                6865221602724658541,
                12220678344985132467,
                13999015252384676179,
                6179576798513002509,
                8807259459797039904,
                13588344625309223635,
                6337304555349855400,
            ],
        ),
        (
            4711,
            4366349731850474445,
            [
                14865558658767168450,
                1570543588041465562,
                11762632291112786969,
                14568021577324630227,
                6038258891483693612,
                8816989010904942752,
                11692864621701489625,
                141362168128559841,
                9079445985709285136,
                14727089087045250198,
                6705393561745281304,
                11594486355944102317,
                3644847148224172158,
                6628338981688287279,
                9719476312151176449,
                10052214408205094771,
                12851822198666506759,
                8074177430413494320,
                17341260920688064693,
                3143845266529271378,
            ],
        ),
        (
            18446744073709551615,
            13937249967666308697,
            [
                478026398904862820,
                13243134898385798468,
                709236020254955927,
                9482188692832154854,
                8350728016034303289,
                11509856835826702400,
                3572905101740335275,
                3136533006162691474,
                16547922902421472960,
                8835741269252529079,
                17926718052445221126,
                8287976176388696091,
                1829913486053131009,
                2741724640216463737,
                12758722211879373259,
                13753320293866618415,
                10699870868994113466,
                6075951796255949511,
                10740104408641833802,
                7420411914205995656,
            ],
        ),
    ];

    #[test]
    fn independent_integer_oracle_covers_four_twist_cycles_and_seed_extrema() {
        for (seed, expected_digest, checkpoints) in REFERENCES {
            let mut random = DecoyRandom::seeded(seed);
            let mut digest = 0xcbf2_9ce4_8422_2325_u64;
            let mut checkpoint = 0;
            for index in 0..1248 {
                let value = random.next_u64();
                if INDICES[checkpoint] == index {
                    assert_eq!(value, checkpoints[checkpoint], "seed{seed}, word{index}");
                    checkpoint += 1;
                }
                for byte in value.to_le_bytes() {
                    digest = (digest ^ u64::from(byte)).wrapping_mul(0x100_0000_01b3);
                }
            }
            assert_eq!(checkpoint, INDICES.len());
            assert_eq!(digest, expected_digest, "seed{seed}");
        }
    }

    #[test]
    fn primary_standard_10000th_default_word_matches() {
        // C++ working draft [rand.predef]/4:
        // https://eel.is/c++draft/rand.predef#4
        // This is a published primary-source literal, not a generated oracle.
        let mut random = DecoyRandom::seeded(5489);
        let value = (0..10_000).map(|_| random.next_u64()).last().unwrap();
        assert_eq!(value, 9_981_545_732_273_789_042);
    }

    #[test]
    fn clone_and_reseed_preserve_all_future_words() {
        for seed in [0, 4711, u64::MAX] {
            let mut random = DecoyRandom::seeded(seed);
            for _ in 0..500 {
                random.next_u64();
            }
            let mut copy = random.clone();
            assert_eq!(random, copy);
            for _ in 0..1000 {
                assert_eq!(random.next_u64(), copy.next_u64());
            }
            random.reseed(seed);
            assert_eq!(random, DecoyRandom::seeded(seed));
            random.reseed(seed.wrapping_add(1));
            assert_eq!(random, DecoyRandom::seeded(seed.wrapping_add(1)));
        }
    }

    #[test]
    fn source_fisher_yates_order_matches_independent_integer_permutations() {
        for (seed, expected) in [
            (0, "ONLTGBCHMEJQFPRIKASD"),
            (4711, "IHEDJOFCGMSTAKRPNLBQ"),
            (u64::MAX, "KEJPLCMBGQSRFDHOITNA"),
        ] {
            let mut random = DecoyRandom::seeded(seed);
            let mut data = *b"ABCDEFGHIJKLMNOPQRST";
            let mut remaining = 100;
            random.shuffle(&mut data, &mut remaining).unwrap();
            assert_eq!(data, expected.as_bytes());
            assert_eq!(remaining, 81);
            let mut raw = DecoyRandom::seeded(seed);
            for _ in 0..19 {
                raw.next_u64();
            }
            assert_eq!(random, raw);
        }
        let mut zeroes = forced(&[0, 0, 0]);
        let mut data = *b"ABCD";
        let mut remaining = 3;
        zeroes.shuffle(&mut data, &mut remaining).unwrap();
        assert_eq!(&data, b"BCDA"); // Ascending swaps would instead give DABC.
        assert_eq!(remaining, 0);
    }

    #[test]
    fn empty_singleton_and_zero_width_mapping_never_draw() {
        let mut random = DecoyRandom::seeded(4711);
        let before = random.clone();
        let mut remaining = 0;
        random.shuffle(&mut [], &mut remaining).unwrap();
        random.shuffle(&mut [b'A'], &mut remaining).unwrap();
        assert_eq!(random.bounded(0, &mut remaining).unwrap(), 0);
        assert_eq!(random, before);
        assert_eq!(remaining, 0);
    }

    #[test]
    fn bucket_boundaries_use_division_instead_of_modulo() {
        for high in [1, 2, 3, 9, 1000, u64::MAX] {
            let bucket = (1_u128 << 64) / (u128::from(high) + 1);
            for raw in [
                0,
                (bucket - 1) as u64,
                bucket as u64,
                u64::MAX - 1,
                u64::MAX,
            ] {
                let expected = u128::from(raw) / bucket;
                if expected > u128::from(high) {
                    continue;
                }
                let mut random = forced(&[raw]);
                let mut remaining = 1;
                assert_eq!(
                    random.bounded(high, &mut remaining).unwrap(),
                    expected as u64
                );
                assert_eq!(remaining, 0);
                assert_eq!(random.index, 1);
            }
        }
        // For [0,1], the top raw word maps to 1 (modulo 2 happens to agree),
        // but bucket-1 maps to 0 even though modulo 2 would be 1.
        let mut random = forced(&[(1_u64 << 63) - 1]);
        assert_eq!(random.bounded(1, &mut 1).unwrap(), 0);
    }

    #[test]
    fn rejected_word_consumes_allowance_and_advances_state_once() {
        // For inclusive range [0,2], floor(2^64/3) leaves raw u64::MAX outside
        // every complete bucket. The next zero word is accepted as 0.
        let mut random = forced(&[u64::MAX, 0]);
        let mut remaining = 2;
        assert_eq!(random.bounded(2, &mut remaining).unwrap(), 0);
        assert_eq!(remaining, 0);
        assert_eq!(random.index, 2);

        let mut random = forced(&[u64::MAX, 0]);
        let mut remaining = 1;
        let error = random.bounded(2, &mut remaining).unwrap_err();
        assert!(error.to_string().contains("draw limit"));
        assert_eq!(remaining, 0);
        assert_eq!(random.index, 1); // Failed precharge does not read second word.
        assert_eq!(random.next_u64(), 0);
    }

    #[test]
    fn exhausted_shuffle_exposes_only_successful_private_swaps_and_draws() {
        let mut random = forced(&[0, 0]);
        let mut data = *b"ABCD";
        let before = random.clone();
        assert!(random.shuffle(&mut data, &mut 0).is_err());
        assert_eq!(&data, b"ABCD");
        assert_eq!(random, before);
        let mut remaining = 1;
        assert!(random.shuffle(&mut data, &mut remaining).is_err());
        assert_eq!(&data, b"DBCA");
        assert_eq!(random.index, 1);
        assert_eq!(remaining, 0);
    }

    #[test]
    fn boost_normalization_sets_redundant_bits_and_repairs_zero_state() {
        let mut random = DecoyRandom {
            state: [0; WORDS],
            index: WORDS,
        };
        random.normalize();
        assert_eq!(random.state[0], 1_u64 << 63);
        assert!(random.state[1..].iter().all(|&word| word == 0));
        let normalized = random.clone();
        random.normalize();
        assert_eq!(random, normalized);
        for seed in [0, 4711, u64::MAX] {
            let mut random = DecoyRandom::seeded(seed);
            let before = random.clone();
            random.normalize();
            assert_eq!(random, before);
            assert_eq!(random.state[0] & UPPER, seed & UPPER);
        }
    }

    // Invert tempering solely to inject adversarial raw words into a valid
    // private read position. No injectable/public RNG abstraction is required.
    fn forced(words: &[u64]) -> DecoyRandom {
        let mut random = DecoyRandom {
            state: [0; WORDS],
            index: 0,
        };
        for (slot, &word) in random.state.iter_mut().zip(words) {
            let mut value = undo_right(word, 43, u64::MAX);
            value = undo_left(value, 37, 0xfff7_eee0_0000_0000);
            value = undo_left(value, 17, 0x71d6_7fff_eda6_0000);
            *slot = undo_right(value, 29, 0x5555_5555_5555_5555);
        }
        random
    }
    fn undo_right(value: u64, shift: u32, mask: u64) -> u64 {
        let mut original = value;
        for _ in 0..64 {
            original = value ^ ((original >> shift) & mask);
        }
        original
    }
    fn undo_left(value: u64, shift: u32, mask: u64) -> u64 {
        let mut original = value;
        for _ in 0..64 {
            original = value ^ ((original << shift) & mask);
        }
        original
    }
}
