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

//! `Math::RandomShuffler`: a portable, seedable Fisher–Yates shuffle.
//!
//! Ports the `RandomShuffler` class of `MATH/MathFunctions.h` (lines 477-506 at
//! core `bc9cc12`); the header's free functions live in
//! [`crate::concept::math_functions`]. The class exists so that a shuffle gives
//! the same permutation on every platform for the same seed, which
//! `std::shuffle` does not promise. [`crate::math::ml::ransac`] draws its
//! samples through it, and `DecoyGenerator` shuffles peptide sequences with the
//! same algorithm (its private engine in `chemistry` delegates to the mapping
//! and loop defined here).
//!
//! Three pieces make the permutation source-exact:
//!
//! 1. raw words from MT19937-64 (`boost::mt19937_64` in the source,
//!    `rand_mt::Mt64` here — the reference generator, whose output the C++
//!    standard fixes);
//! 2. Boost `uniform_int(0, i)`'s bucket division with rejection, *not* modulo
//!    and not `rand`'s widening multiply;
//! 3. the source's descending loop `for i = n-1 .. 1: swap(a[i], a[d(0, i)])`.
//!
//! See `docs/RANSAC_SUPPORT.md` for the evidence.

use crate::{Error, Result};
use rand_mt::Mt64;

/// Boost's `mt19937_64::default_seed`, used by the source's default
/// constructor.
pub const DEFAULT_SEED: u64 = 5489;

/// A seedable shuffler whose permutations match the source's on every
/// platform.
///
/// Equality compares the whole engine state, so two shufflers are equal
/// exactly when they will produce the same future permutations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RandomShuffler {
    engine: Mt64,
}

impl Default for RandomShuffler {
    /// The source's `RandomShuffler()`: an engine seeded with
    /// [`DEFAULT_SEED`], as a default-constructed `boost::mt19937_64` is.
    fn default() -> Self {
        Self::new(DEFAULT_SEED)
    }
}

impl RandomShuffler {
    /// A shuffler whose engine is in the state of `boost::mt19937_64(seed)`.
    ///
    /// The source's seeding constructor takes an `int`, so a 64-bit seed
    /// passed to it is first narrowed; [`RandomShuffler::from_int_seed`]
    /// reproduces that conversion. This constructor takes the full 64-bit
    /// seed, as the source's [`seed`](RandomShuffler::seed) member does.
    pub fn new(seed: u64) -> Self {
        Self {
            engine: Mt64::new(seed),
        }
    }

    /// The source's `explicit RandomShuffler(int seed)`.
    ///
    /// `boost::mt19937_64`'s seeding constructor takes a `uint64_t`, so a
    /// negative `int` arrives sign-extended: `-1` seeds with `u64::MAX`.
    pub fn from_int_seed(seed: i32) -> Self {
        Self::new(i64::from(seed) as u64)
    }

    /// Restart the engine as `boost::mt19937_64::seed(val)` does; the next
    /// permutation is the first one [`RandomShuffler::new`] with the same seed
    /// would produce.
    pub fn seed(&mut self, val: u64) {
        self.engine.reseed(val);
    }

    /// One raw engine word.
    ///
    /// The source exposes its engine as the public member `rng_`; this is the
    /// one operation a caller can perform on it without changing its state
    /// in another way.
    pub fn next_u64(&mut self) -> u64 {
        self.engine.next_u64()
    }

    /// Shuffle `data` in place: the source's `portable_random_shuffle`.
    ///
    /// For `i` from `len - 1` down to `1`, draws `j` uniformly from `[0, i]`
    /// with Boost's `uniform_int` mapping and swaps `data[i]` with `data[j]`.
    /// Empty and one-element slices draw nothing. Each position normally costs
    /// one engine word; Boost's rejection step can cost more, with probability
    /// below `(i + 1) / 2^64` per draw.
    pub fn portable_random_shuffle<T>(&mut self, data: &mut [T]) {
        for end in (1..data.len()).rev() {
            // A slice length always fits in u64 on every supported target.
            let selected = bounded_unbudgeted(|| self.engine.next_u64(), end as u64) as usize;
            data.swap(end, selected);
        }
    }
}

/// Boost `uniform_int(0, high)` over a full-range 64-bit engine, without a
/// draw budget.
///
/// Rejection continues until a word lands in a complete bucket; the chance
/// that a word is rejected is `(2^64 mod (high + 1)) / 2^64`.
fn bounded_unbudgeted(mut next: impl FnMut() -> u64, high: u64) -> u64 {
    if high == 0 {
        return 0;
    }
    let bucket = bucket_size(high);
    loop {
        let result = next() / bucket;
        if result <= high {
            return result;
        }
    }
}

/// Boost's bucket size for `uniform_int(0, high)` over the range `[0, 2^64)`,
/// `floor(2^64 / (high + 1))`, for `high >= 1`.
///
/// Boost computes it as `brange / (range + 1)`, incremented when
/// `brange % (range + 1) == range`, with `brange = 2^64 - 1`; that is the same
/// integer. `u128` expresses `2^64` directly and the quotient fits a nonzero
/// `u64` for every `high >= 1`.
fn bucket_size(high: u64) -> u64 {
    ((1_u128 << 64) / (u128::from(high) + 1)) as u64
}

/// The source's descending Fisher–Yates loop over caller-supplied raw words,
/// with a draw budget.
///
/// Every raw draw, including a rejected one, consumes one unit of
/// `remaining_draws` before the word is read. On error `data` keeps the swaps
/// already made; callers that need atomicity stage a copy. Empty and
/// one-element slices draw nothing.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the budget is exhausted before the
/// shuffle finishes.
pub(crate) fn shuffle_from<T>(
    mut next: impl FnMut() -> u64,
    data: &mut [T],
    remaining_draws: &mut usize,
) -> Result<()> {
    for end in (1..data.len()).rev() {
        let high = u64::try_from(end)
            .map_err(|_| Error::InvalidValue("shuffle range exceeds u64".into()))?;
        let selected = bounded_from(&mut next, high, remaining_draws)? as usize;
        // selected <= end, so the usize conversion cannot truncate.
        data.swap(end, selected);
    }
    Ok(())
}

/// Boost `uniform_int(0, high)` with a draw budget; see [`shuffle_from`].
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the budget runs out.
pub(crate) fn bounded_from(
    mut next: impl FnMut() -> u64,
    high: u64,
    remaining_draws: &mut usize,
) -> Result<u64> {
    if high == 0 {
        return Ok(0);
    }
    let bucket = bucket_size(high);
    loop {
        *remaining_draws = remaining_draws
            .checked_sub(1)
            .ok_or_else(|| Error::InvalidValue("random draw limit exceeded".into()))?;
        let result = next() / bucket;
        if result <= high {
            return Ok(result);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unbudgeted_and_budgeted_paths_choose_the_same_permutation() {
        for seed in [0, 4711, u64::MAX] {
            let mut public = RandomShuffler::new(seed);
            let mut data: Vec<u32> = (0..257).collect();
            public.portable_random_shuffle(&mut data);

            let mut raw = Mt64::new(seed);
            let mut staged: Vec<u32> = (0..257).collect();
            let mut remaining = usize::MAX;
            shuffle_from(|| raw.next_u64(), &mut staged, &mut remaining).unwrap();
            assert_eq!(data, staged);
            assert_eq!(public.engine, raw);
        }
    }

    #[test]
    fn rejection_repeats_until_a_complete_bucket() {
        // For [0,2], floor(2^64/3) leaves u64::MAX outside every bucket.
        let words = [u64::MAX, u64::MAX, 7];
        let mut index = 0;
        let value = bounded_unbudgeted(
            || {
                index += 1;
                words[index - 1]
            },
            2,
        );
        assert_eq!(value, 0);
        assert_eq!(index, 3);
    }

    #[test]
    fn default_is_the_boost_default_seed() {
        assert_eq!(RandomShuffler::default(), RandomShuffler::new(5489));
        assert_eq!(
            RandomShuffler::from_int_seed(-1),
            RandomShuffler::new(u64::MAX)
        );
        assert_eq!(RandomShuffler::from_int_seed(42), RandomShuffler::new(42));
    }
}
