// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// Copyright Jens Maurer 2000-2001
// Copyright Steven Watanabe 2011
// SPDX-License-Identifier: BSD-3-Clause AND BSL-1.0
// $Maintainer: OpenMS Rust contributors $
//
// The bounded-integer mapping derives from Boost.Random 1.90
// uniform_int_distribution.hpp, distributed under the Boost Software License,
// Version 1.0 (https://www.boost.org/LICENSE_1_0.txt). The MT19937-64 engine is
// the rand_mt crate (MIT OR Apache-2.0); no Boost engine code remains here.
//! Private deterministic engine for OpenMS Math::RandomShuffler semantics.
//! The outer decoy generator stages this small state and its output so failed
//! draw/work checks remain atomic at the public operation boundary.
//!
//! Raw words come from `rand_mt::Mt64`, the reference MT19937-64, in place of
//! the source's `boost::mt19937_64`. Only Boost `uniform_int`'s range mapping
//! and the source's descending Fisher–Yates loop are implemented here: `rand`'s
//! range mapping and slice shuffle choose different values and would change
//! every decoy.

use crate::{Error, Result};
use rand_mt::Mt64;

/// MT19937-64 word stream with the source shuffler's range mapping on top.
///
/// Boost's `mersenne_twister_engine::seed` also rewrites the low 31 bits of the
/// first state word and repairs an all-zero state; `Mt64` does neither, and
/// neither changes any output. The twist reads only the upper 33 bits of that
/// word before replacing it, and the repair cannot trigger after an integer
/// seed: a zero second word forces the third to be 2. Equality is unaffected
/// too, because the second seeded word alone determines the seed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DecoyRandom {
    engine: Mt64,
}

impl DecoyRandom {
    /// Engine in the state of `boost::mt19937_64(seed)`.
    pub(crate) fn seeded(seed: u64) -> Self {
        Self {
            engine: Mt64::new(seed),
        }
    }

    /// Restart the stream as `RandomShuffler::seed` does; the next word is the
    /// first word of [`DecoyRandom::seeded`] with the same seed.
    pub(crate) fn reseed(&mut self, seed: u64) {
        self.engine.reseed(seed);
    }

    /// Descending Fisher–Yates, as in MathFunctions.h. Every raw draw, including
    /// rejection, consumes one allowance before advancing the engine. Empty and
    /// singleton slices use no draws. On error the private slice/state may have
    /// advanced through earlier swaps; the public caller stages both.
    pub(super) fn shuffle(&mut self, data: &mut [u8], remaining_draws: &mut usize) -> Result<()> {
        shuffle_from(|| self.next_u64(), data, remaining_draws)
    }

    /// One raw 64-bit engine word.
    pub(crate) fn next_u64(&mut self) -> u64 {
        self.engine.next_u64()
    }
}

// The Fisher–Yates loop and range mapping read raw words through a closure so
// the private tests can script adversarial words; production passes the engine.
fn shuffle_from(
    mut next: impl FnMut() -> u64,
    data: &mut [u8],
    remaining_draws: &mut usize,
) -> Result<()> {
    for end in (1..data.len()).rev() {
        let high = u64::try_from(end)
            .map_err(|_| Error::InvalidValue("decoy shuffle range exceeds u64".into()))?;
        let selected = bounded_from(&mut next, high, remaining_draws)? as usize;
        // selected <= end, so the native usize conversion cannot truncate.
        data.swap(end, selected);
    }
    Ok(())
}

// Boost uniform_int(0, high) uses bucket division/rejection, not modulo.
// u128 expresses 2^64 directly; for high>=1 the bucket fits a nonzero u64.
fn bounded_from(
    mut next: impl FnMut() -> u64,
    high: u64,
    remaining_draws: &mut usize,
) -> Result<u64> {
    if high == 0 {
        return Ok(0);
    }
    let bucket = ((1_u128 << 64) / (u128::from(high) + 1)) as u64;
    loop {
        *remaining_draws = remaining_draws
            .checked_sub(1)
            .ok_or_else(|| Error::InvalidValue("decoy random draw limit exceeded".into()))?;
        let result = next() / bucket;
        if result <= high {
            return Ok(result);
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
        let script = [0, 0, 0];
        let mut zeroes = Words::new(&script);
        let mut data = *b"ABCD";
        let mut remaining = 3;
        shuffle_from(|| zeroes.draw(), &mut data, &mut remaining).unwrap();
        assert_eq!(&data, b"BCDA"); // Ascending swaps would instead give DABC.
        assert_eq!(remaining, 0);
        assert_eq!(zeroes.read, 3);
    }

    #[test]
    fn empty_singleton_and_zero_width_mapping_never_draw() {
        let mut random = DecoyRandom::seeded(4711);
        let before = random.clone();
        let mut remaining = 0;
        random.shuffle(&mut [], &mut remaining).unwrap();
        let mut singleton = *b"A";
        random.shuffle(&mut singleton, &mut remaining).unwrap();
        assert_eq!(
            bounded_from(|| random.next_u64(), 0, &mut remaining).unwrap(),
            0
        );
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
                let script = [raw];
                let mut words = Words::new(&script);
                let mut remaining = 1;
                assert_eq!(
                    bounded_from(|| words.draw(), high, &mut remaining).unwrap(),
                    expected as u64
                );
                assert_eq!(remaining, 0);
                assert_eq!(words.read, 1);
            }
        }
        // For [0,1], the top raw word maps to 1 (modulo 2 happens to agree),
        // but bucket-1 maps to 0 even though modulo 2 would be 1.
        let script = [(1_u64 << 63) - 1];
        let mut words = Words::new(&script);
        assert_eq!(bounded_from(|| words.draw(), 1, &mut 1).unwrap(), 0);
    }

    #[test]
    fn rejected_word_consumes_allowance_and_advances_state_once() {
        // For inclusive range [0,2], floor(2^64/3) leaves raw u64::MAX outside
        // every complete bucket. The next zero word is accepted as 0.
        let script = [u64::MAX, 0];
        let mut words = Words::new(&script);
        let mut remaining = 2;
        assert_eq!(bounded_from(|| words.draw(), 2, &mut remaining).unwrap(), 0);
        assert_eq!(remaining, 0);
        assert_eq!(words.read, 2);

        let mut words = Words::new(&script);
        let mut remaining = 1;
        let error = bounded_from(|| words.draw(), 2, &mut remaining).unwrap_err();
        assert!(error.to_string().contains("draw limit"));
        assert_eq!(remaining, 0);
        assert_eq!(words.read, 1); // Failed precharge does not read second word.
        assert_eq!(words.draw(), 0);
    }

    #[test]
    fn exhausted_shuffle_exposes_only_successful_private_swaps_and_draws() {
        let script = [0, 0];
        let mut words = Words::new(&script);
        let mut data = *b"ABCD";
        assert!(shuffle_from(|| words.draw(), &mut data, &mut 0).is_err());
        assert_eq!(&data, b"ABCD");
        assert_eq!(words.read, 0);
        let mut remaining = 1;
        assert!(shuffle_from(|| words.draw(), &mut data, &mut remaining).is_err());
        assert_eq!(&data, b"DBCA");
        assert_eq!(words.read, 1);
        assert_eq!(remaining, 0);

        // The engine-backed method charges before reading in the same way: no
        // allowance leaves the engine untouched, one allowance reads one word.
        let mut random = DecoyRandom::seeded(4711);
        let before = random.clone();
        let mut data = *b"ABCD";
        assert!(random.shuffle(&mut data, &mut 0).is_err());
        assert_eq!(random, before);
        let mut remaining = 1;
        assert!(random.shuffle(&mut data, &mut remaining).is_err());
        assert_eq!(remaining, 0);
        let mut advanced = before;
        advanced.next_u64();
        assert_eq!(random, advanced);
    }

    #[test]
    fn engine_keeps_the_size_charged_by_decoy_work_accounting() {
        // decoy_generator.rs charges size_of::<DecoyRandom>() for each staged
        // and per-product engine. Pin it to the replaced [u64; 312] + index
        // layout so a rand_mt upgrade cannot silently move decoy work limits.
        assert_eq!(size_of::<DecoyRandom>(), size_of::<([u64; 312], usize)>());
    }

    // Scripted raw words for the range mapping and Fisher–Yates loop. Reading
    // past the script panics, so a passing test also proves no extra draw.
    struct Words<'a> {
        words: &'a [u64],
        read: usize,
    }

    impl<'a> Words<'a> {
        fn new(words: &'a [u64]) -> Self {
            Self { words, read: 0 }
        }

        fn draw(&mut self) -> u64 {
            let word = self.words[self.read];
            self.read += 1;
            word
        }
    }
}
