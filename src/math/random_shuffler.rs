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
//! core `bc9cc12`). The class exists so that a seed gives the same permutation
//! on every platform and standard library, which `std::shuffle` does not
//! promise. [`crate::math::ml::ransac`] draws its samples through it, which is
//! why it lives in `math`: that module must stay a sink of the module graph.
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
//! The crate-private `DecoyRandom` engine of the decoy generator adds a draw
//! budget on top of the same engine; it shares this module's descending loop
//! (`descending_fisher_yates`) and bucket size (`bucket_size`), so the source's
//! mapping exists once.

pub use rand_mt::Mt64;

/// A seedable shuffler whose permutations match the source's on every
/// platform.
///
/// Equality compares the whole engine state, so two shufflers are equal
/// exactly when they will produce the same future permutations.
///
/// Boost's `mersenne_twister_engine::seed` also rewrites the low 31 bits of the
/// first state word and repairs an all-zero state; `Mt64` does neither, and
/// neither changes any output: the twist reads only the upper 33 bits of that
/// word before replacing it, and the repair cannot trigger after an integer
/// seed.
///
/// ```
/// use openms::math::random_shuffler::RandomShuffler;
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
    /// The source's `RandomShuffler() = default`: a default-constructed
    /// `boost::mt19937_64`, seeded with [`RandomShuffler::DEFAULT_SEED`].
    fn default() -> Self {
        Self::with_seed(Self::DEFAULT_SEED)
    }
}

impl RandomShuffler {
    /// Boost's `mt19937_64::default_seed`, used by the default constructor.
    pub const DEFAULT_SEED: u64 = 5489;

    /// The source's `explicit RandomShuffler(int seed)`.
    ///
    /// `boost::mt19937_64`'s seeding constructor takes a `uint64_t`, so a
    /// negative `int` arrives sign-extended: `new(-1)` equals
    /// `with_seed(u64::MAX)`.
    pub fn new(seed: i32) -> Self {
        Self::with_seed(i64::from(seed) as u64)
    }

    /// A shuffler whose engine is in the state of `boost::mt19937_64(seed)`
    /// for any 64-bit seed.
    ///
    /// Native addition: the source's constructor takes an `int`, while its
    /// [`seed`](RandomShuffler::seed) member accepts the full `uint64_t` range.
    pub fn with_seed(seed: u64) -> Self {
        Self {
            engine: Mt64::new(seed),
        }
    }

    /// The source's `explicit RandomShuffler(const boost::mt19937_64&)`: adopt
    /// an engine in whatever state it is in.
    pub fn from_engine(engine: Mt64) -> Self {
        Self { engine }
    }

    /// Read access to the engine: the source's public member `rng_`.
    pub fn engine(&self) -> &Mt64 {
        &self.engine
    }

    /// Mutable access to the engine: the source's public member `rng_`.
    /// Drawing from it advances the stream the shuffle draws from, as using
    /// `rng_` directly does in the source.
    pub fn engine_mut(&mut self) -> &mut Mt64 {
        &mut self.engine
    }

    /// Give up the shuffler and return its engine.
    pub fn into_engine(self) -> Mt64 {
        self.engine
    }

    /// The source's `seed(uint64_t val)`: restart the stream; the next word is
    /// the first word of [`RandomShuffler::with_seed`] with the same seed.
    pub fn seed(&mut self, val: u64) {
        self.engine.reseed(val);
    }

    /// One raw 64-bit engine word (the source's `rng_()`).
    pub fn next_u64(&mut self) -> u64 {
        self.engine.next_u64()
    }

    /// Shuffle `data` in place: the source's `portable_random_shuffle`.
    ///
    /// For `i` from `len - 1` down to `1`, draws `j` uniformly from `[0, i]`
    /// with Boost's `uniform_int` mapping and swaps `data[i]` with `data[j]`.
    /// Empty and one-element slices draw nothing. The permutation depends only
    /// on the seed and the length, never on the element type or the platform.
    /// Each position normally costs one engine word; Boost's rejection step can
    /// cost more, with probability below `(i + 1) / 2^64` per draw, so the cost
    /// is linear in the length.
    pub fn portable_random_shuffle<T>(&mut self, data: &mut [T]) {
        // usize is at most 64 bits on every supported target, so the index and
        // the drawn value (at most `end`) convert without loss.
        let _: Result<(), std::convert::Infallible> = descending_fisher_yates(data, |end| {
            Ok(bounded(&mut self.engine, end as u64) as usize)
        });
    }
}

/// The source's descending Fisher-Yates loop,
/// `for i = n-1 .. 1: swap(a[i], a[select(i)])`, shared by every shuffle in the
/// crate that must reproduce `Math::RandomShuffler`: this module's
/// [`RandomShuffler`] and the decoy generators' budgeted engine.
///
/// `select(end)` returns the index drawn from `[0, end]`; a larger value is
/// clamped to `end`, so a faulty selector cannot index out of range. Empty and
/// one-element slices never call it. An error from `select` stops the loop with
/// the earlier swaps already applied.
pub(crate) fn descending_fisher_yates<T, E>(
    data: &mut [T],
    mut select: impl FnMut(usize) -> Result<usize, E>,
) -> Result<(), E> {
    for end in (1..data.len()).rev() {
        let selected = select(end)?;
        data.swap(end, selected.min(end));
    }
    Ok(())
}

/// Boost's `uniform_int(0, high)` bucket size over a full-range 64-bit engine:
/// `floor(2^64 / (high + 1))`, computed in `u128` so that `2^64` is
/// representable. For `high >= 1` the quotient fits a nonzero `u64`; `high` of
/// zero is the caller's special case (Boost consumes no word for it).
pub(crate) fn bucket_size(high: u64) -> u64 {
    ((1_u128 << 64) / (u128::from(high) + 1)) as u64
}

/// Boost `uniform_int(0, high)` over a full-range 64-bit engine.
///
/// Boost computes the bucket size as `brange / (range + 1)`, plus one when
/// `brange % (range + 1) == range`, with `brange = 2^64 - 1`: that is
/// `floor(2^64 / (high + 1))`, computed here directly in `u128`. A word that
/// lands in the final partial bucket is rejected and another one drawn.
fn bounded(engine: &mut Mt64, high: u64) -> u64 {
    if high == 0 {
        return 0;
    }
    let bucket = bucket_size(high);
    loop {
        let result = engine.next_u64() / bucket;
        if result <= high {
            return result;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructors_and_seeding() {
        assert_eq!(RandomShuffler::default(), RandomShuffler::with_seed(5489));
        assert_eq!(RandomShuffler::new(-1), RandomShuffler::with_seed(u64::MAX));
        assert_eq!(RandomShuffler::new(42), RandomShuffler::with_seed(42));
        let mut reseeded = RandomShuffler::new(1);
        reseeded.seed(7);
        assert_eq!(reseeded, RandomShuffler::with_seed(7));
        let engine = RandomShuffler::with_seed(9).into_engine();
        assert_eq!(
            RandomShuffler::from_engine(engine.clone()).engine(),
            &engine
        );
    }

    #[test]
    fn bucket_division_not_modulo() {
        // For [0,1] the bucket is 2^63: the word 2^63 - 1 maps to 0, although
        // modulo 2 would give 1.
        let high = 1_u64;
        let bucket = ((1_u128 << 64) / (u128::from(high) + 1)) as u64;
        assert_eq!(((1_u64 << 63) - 1) / bucket, 0);
        assert_eq!(u64::MAX / bucket, 1);
    }
}
