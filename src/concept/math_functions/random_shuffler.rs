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
//! `Math::RandomShuffler`: a platform-independent seeded shuffle.
//!
//! Ported from the `RandomShuffler` class of `MATH/MathFunctions.h` at
//! `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The source wraps a
//! `boost::mt19937_64` engine and shuffles with its own descending
//! Fisher–Yates loop over Boost `uniform_int`, precisely so that a seed gives
//! the same permutation on every platform and standard library (unlike
//! `std::shuffle`). This port keeps that contract bit for bit.
//!
//! Raw words come from `rand_mt::Mt64`, the reference MT19937-64, in place of
//! the source's `boost::mt19937_64`. Only Boost `uniform_int`'s range mapping
//! and the source's descending Fisher–Yates loop are implemented here: `rand`'s
//! range mapping and slice shuffle choose different values and would change
//! every permutation.
//!
//! This type was the crate-private `DecoyRandom` of `chemistry::decoy_generator`
//! until it was promoted here; the decoy generator now uses it through the
//! crate-private draw-limited [`RandomShuffler::shuffle_limited`].

use crate::{Error, Result};
use std::convert::Infallible;

pub use rand_mt::Mt64;

/// Seeded MT19937-64 shuffler with the source's portable permutation order.
///
/// Source `Math::RandomShuffler`. Two shufflers compare equal when their engine
/// states are equal, so equal shufflers produce identical futures.
///
/// Boost's `mersenne_twister_engine::seed` also rewrites the low 31 bits of the
/// first state word and repairs an all-zero state; `Mt64` does neither, and
/// neither changes any output. The twist reads only the upper 33 bits of that
/// word before replacing it, and the repair cannot trigger after an integer
/// seed: a zero second word forces the third to be 2. Equality is unaffected
/// too, because the second seeded word alone determines the seed.
///
/// ```
/// use openms::concept::math_functions::RandomShuffler;
///
/// let mut shuffler = RandomShuffler::new(0);
/// let mut items = [1, 2, 3, 4, 5];
/// shuffler.portable_random_shuffle(&mut items);
/// let mut again = RandomShuffler::new(0);
/// let mut copy = [1, 2, 3, 4, 5];
/// again.portable_random_shuffle(&mut copy);
/// assert_eq!(items, copy);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RandomShuffler {
    engine: Mt64,
}

impl Default for RandomShuffler {
    /// Source `RandomShuffler() = default`: a default-constructed
    /// `boost::mt19937_64`, which is seeded with the standard default seed
    /// 5489 ([`RandomShuffler::DEFAULT_SEED`]).
    fn default() -> Self {
        Self::with_seed(Self::DEFAULT_SEED)
    }
}

impl RandomShuffler {
    /// Seed of a default-constructed engine (`mersenne_twister_engine::default_seed`).
    pub const DEFAULT_SEED: u64 = 5489;

    /// Source `explicit RandomShuffler(int seed)`.
    ///
    /// The `int` converts to the engine's 64-bit result type as C++ does, by
    /// sign extension, so `new(-1)` is seeded with `u64::MAX` and equals
    /// `with_seed(u64::MAX)`.
    pub fn new(seed: i32) -> Self {
        Self::with_seed(i64::from(seed) as u64)
    }

    /// Engine in the state of `boost::mt19937_64(seed)` for any 64-bit seed.
    ///
    /// Native addition: the source constructor takes an `int`, but
    /// [`RandomShuffler::seed`] accepts the full `uint64_t` range, and callers
    /// such as the decoy generator seed with values beyond `int`.
    pub fn with_seed(seed: u64) -> Self {
        Self {
            engine: Mt64::new(seed),
        }
    }

    /// Source `explicit RandomShuffler(const boost::mt19937_64& mt_rng)`: adopt
    /// an engine in whatever state it is in.
    pub fn from_engine(engine: Mt64) -> Self {
        Self { engine }
    }

    /// Read access to the engine; source public member `rng_`.
    pub fn engine(&self) -> &Mt64 {
        &self.engine
    }

    /// Mutable access to the engine; source public member `rng_`. Drawing from
    /// it advances the same stream the shuffle draws from, exactly as using
    /// `rng_` directly does in the source.
    pub fn engine_mut(&mut self) -> &mut Mt64 {
        &mut self.engine
    }

    /// Give up the shuffler and return its engine.
    pub fn into_engine(self) -> Mt64 {
        self.engine
    }

    /// Source `seed(uint64_t val)`: restart the stream. The next word equals the
    /// first word of [`RandomShuffler::with_seed`] with the same seed.
    pub fn seed(&mut self, val: u64) {
        self.engine.reseed(val);
    }

    /// One raw 64-bit engine word (source `rng_()`).
    pub fn next_u64(&mut self) -> u64 {
        self.engine.next_u64()
    }

    /// Source `portable_random_shuffle(first, last)`: shuffle `data` in place.
    ///
    /// The loop runs `i` from `len - 1` down to `1` and swaps `data[i]` with
    /// `data[j]`, where `j` is drawn from the inclusive range `[0, i]` by Boost
    /// `uniform_int`'s bucket division with rejection (not a modulo). Empty and
    /// one-element slices draw nothing. The permutation depends only on the
    /// seed and the length, never on the element type or the platform.
    ///
    /// A rejected raw word is drawn again, as in the source. Rejection needs a
    /// raw word in the final partial bucket, which has probability below
    /// `(i + 1) / 2^64` per draw, so the loop's cost is linear in the length.
    pub fn portable_random_shuffle<T>(&mut self, data: &mut [T]) {
        let engine = &mut self.engine;
        let unlimited = || -> std::result::Result<(), Infallible> { Ok(()) };
        match shuffle_from(|| engine.next_u64(), data, unlimited) {
            Ok(()) => {}
            Err(never) => match never {},
        }
    }

    /// Descending Fisher–Yates with a cumulative raw-draw allowance.
    ///
    /// Every raw draw, including a rejected one, consumes one allowance before
    /// the engine advances. On error the slice and the engine may have advanced
    /// through earlier swaps; callers that need atomicity stage both, as the
    /// decoy generator does.
    pub(crate) fn shuffle_limited<T>(
        &mut self,
        data: &mut [T],
        remaining_draws: &mut usize,
    ) -> Result<()> {
        let engine = &mut self.engine;
        shuffle_from(|| engine.next_u64(), data, || charge(remaining_draws))
    }
}

fn charge(remaining_draws: &mut usize) -> Result<()> {
    *remaining_draws = remaining_draws
        .checked_sub(1)
        .ok_or_else(|| Error::InvalidValue("random shuffle draw limit exceeded".into()))?;
    Ok(())
}

// The Fisher–Yates loop and range mapping read raw words through a closure so
// the private tests can script adversarial words; production passes the engine.
// `charge` runs before every raw draw; it is the draw allowance, or a no-op.
fn shuffle_from<T, E>(
    mut next: impl FnMut() -> u64,
    data: &mut [T],
    mut charge: impl FnMut() -> std::result::Result<(), E>,
) -> std::result::Result<(), E> {
    for end in (1..data.len()).rev() {
        // usize is at most 64 bits on every supported target, so the index and
        // the drawn value (at most `end`) convert without loss both ways.
        let selected = bounded_from(&mut next, end as u64, &mut charge)? as usize;
        data.swap(end, selected);
    }
    Ok(())
}

// Boost uniform_int(0, high) uses bucket division/rejection, not modulo.
// u128 expresses 2^64 directly; for high>=1 the bucket fits a nonzero u64.
// Boost computes brange / (high + 1), plus one when the remainder equals high;
// that is floor(2^64 / (high + 1)), which this computes directly.
fn bounded_from<E>(
    mut next: impl FnMut() -> u64,
    high: u64,
    mut charge: impl FnMut() -> std::result::Result<(), E>,
) -> std::result::Result<u64, E> {
    if high == 0 {
        return Ok(0);
    }
    let bucket = ((1_u128 << 64) / (u128::from(high) + 1)) as u64;
    loop {
        charge()?;
        let result = next() / bucket;
        if result <= high {
            return Ok(result);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn limited(remaining: &mut usize) -> impl FnMut() -> Result<()> + '_ {
        move || charge(remaining)
    }

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
            let mut random = RandomShuffler::with_seed(seed);
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
        let mut random = RandomShuffler::default();
        let value = (0..10_000).map(|_| random.next_u64()).last().unwrap();
        assert_eq!(value, 9_981_545_732_273_789_042);
    }

    #[test]
    fn clone_and_reseed_preserve_all_future_words() {
        for seed in [0, 4711, u64::MAX] {
            let mut random = RandomShuffler::with_seed(seed);
            for _ in 0..500 {
                random.next_u64();
            }
            let mut copy = random.clone();
            assert_eq!(random, copy);
            for _ in 0..1000 {
                assert_eq!(random.next_u64(), copy.next_u64());
            }
            random.seed(seed);
            assert_eq!(random, RandomShuffler::with_seed(seed));
            random.seed(seed.wrapping_add(1));
            assert_eq!(random, RandomShuffler::with_seed(seed.wrapping_add(1)));
        }
    }

    #[test]
    fn source_fisher_yates_order_matches_independent_integer_permutations() {
        for (seed, expected) in [
            (0, "ONLTGBCHMEJQFPRIKASD"),
            (4711, "IHEDJOFCGMSTAKRPNLBQ"),
            (u64::MAX, "KEJPLCMBGQSRFDHOITNA"),
        ] {
            let mut random = RandomShuffler::with_seed(seed);
            let mut data = *b"ABCDEFGHIJKLMNOPQRST";
            let mut remaining = 100;
            random.shuffle_limited(&mut data, &mut remaining).unwrap();
            assert_eq!(data, expected.as_bytes());
            assert_eq!(remaining, 81);
            let mut raw = RandomShuffler::with_seed(seed);
            for _ in 0..19 {
                raw.next_u64();
            }
            assert_eq!(random, raw);

            // The unlimited public shuffle takes the same path.
            let mut public = RandomShuffler::with_seed(seed);
            let mut data = *b"ABCDEFGHIJKLMNOPQRST";
            public.portable_random_shuffle(&mut data);
            assert_eq!(data, expected.as_bytes());
            assert_eq!(public, raw);
        }
        let script = [0, 0, 0];
        let mut zeroes = Words::new(&script);
        let mut data = *b"ABCD";
        let mut remaining = 3;
        shuffle_from(|| zeroes.draw(), &mut data, limited(&mut remaining)).unwrap();
        assert_eq!(&data, b"BCDA"); // Ascending swaps would instead give DABC.
        assert_eq!(remaining, 0);
        assert_eq!(zeroes.read, 3);
    }

    #[test]
    fn empty_singleton_and_zero_width_mapping_never_draw() {
        let mut random = RandomShuffler::with_seed(4711);
        let before = random.clone();
        let mut remaining = 0;
        random
            .shuffle_limited::<u8>(&mut [], &mut remaining)
            .unwrap();
        let mut singleton = *b"A";
        random
            .shuffle_limited(&mut singleton, &mut remaining)
            .unwrap();
        random.portable_random_shuffle::<u8>(&mut []);
        random.portable_random_shuffle(&mut singleton);
        let engine = &mut random.engine;
        assert_eq!(
            bounded_from(|| engine.next_u64(), 0, limited(&mut remaining)).unwrap(),
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
                    bounded_from(|| words.draw(), high, limited(&mut remaining)).unwrap(),
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
        let mut remaining = 1;
        assert_eq!(
            bounded_from(|| words.draw(), 1, limited(&mut remaining)).unwrap(),
            0
        );
    }

    #[test]
    fn rejected_word_consumes_allowance_and_advances_state_once() {
        // For inclusive range [0,2], floor(2^64/3) leaves raw u64::MAX outside
        // every complete bucket. The next zero word is accepted as 0.
        let script = [u64::MAX, 0];
        let mut words = Words::new(&script);
        let mut remaining = 2;
        assert_eq!(
            bounded_from(|| words.draw(), 2, limited(&mut remaining)).unwrap(),
            0
        );
        assert_eq!(remaining, 0);
        assert_eq!(words.read, 2);

        let mut words = Words::new(&script);
        let mut remaining = 1;
        let error = bounded_from(|| words.draw(), 2, limited(&mut remaining)).unwrap_err();
        assert!(error.to_string().contains("draw limit"));
        assert_eq!(remaining, 0);
        assert_eq!(words.read, 1); // Failed precharge does not read second word.
        assert_eq!(words.draw(), 0);

        // The unlimited mapping redraws after the same rejection.
        let mut words = Words::new(&script);
        let unlimited = || -> std::result::Result<(), Infallible> { Ok(()) };
        assert_eq!(bounded_from(|| words.draw(), 2, unlimited), Ok(0));
        assert_eq!(words.read, 2);
    }

    #[test]
    fn exhausted_shuffle_exposes_only_successful_private_swaps_and_draws() {
        let script = [0, 0];
        let mut words = Words::new(&script);
        let mut data = *b"ABCD";
        let mut remaining = 0;
        assert!(shuffle_from(|| words.draw(), &mut data, limited(&mut remaining)).is_err());
        assert_eq!(&data, b"ABCD");
        assert_eq!(words.read, 0);
        let mut remaining = 1;
        assert!(shuffle_from(|| words.draw(), &mut data, limited(&mut remaining)).is_err());
        assert_eq!(&data, b"DBCA");
        assert_eq!(words.read, 1);
        assert_eq!(remaining, 0);

        // The engine-backed method charges before reading in the same way: no
        // allowance leaves the engine untouched, one allowance reads one word.
        let mut random = RandomShuffler::with_seed(4711);
        let before = random.clone();
        let mut data = *b"ABCD";
        assert!(random.shuffle_limited(&mut data, &mut 0).is_err());
        assert_eq!(random, before);
        let mut remaining = 1;
        assert!(random.shuffle_limited(&mut data, &mut remaining).is_err());
        assert_eq!(remaining, 0);
        let mut advanced = before;
        advanced.next_u64();
        assert_eq!(random, advanced);
    }

    #[test]
    fn engine_keeps_the_size_charged_by_decoy_work_accounting() {
        // decoy_generator.rs charges size_of::<RandomShuffler>() for each staged
        // and per-product engine. Pin it to the replaced [u64; 312] + index
        // layout so a rand_mt upgrade cannot silently move decoy work limits.
        assert_eq!(
            size_of::<RandomShuffler>(),
            size_of::<([u64; 312], usize)>()
        );
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
