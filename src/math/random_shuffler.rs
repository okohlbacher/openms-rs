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
//! The type lives in [`crate::concept::math_functions::RandomShuffler`], where
//! the decoy generator's draw-budgeted shuffle is part of the same
//! implementation; this module is the name `math` code (RANSAC, the OpenSWATH
//! helpers) reaches it by, and the one place the source's mapping exists.

pub use crate::concept::math_functions::{Mt64, RandomShuffler};

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
