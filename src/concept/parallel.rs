// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Data parallelism with a determinism contract.
//!
//! OpenMS parallelises with OpenMP: 176 pragmas across the SDK, about 100 of
//! them in ANALYSIS. Their shapes are lopsided — roughly 40 `parallel for`
//! loops against some 90 `critical` sections and 11 `atomic` updates, and only
//! **two** `reduction` clauses in the whole tree.
//!
//! That distribution is the design input. Nearly every `critical` section
//! exists to serialise pushes into a shared container, which is a C++ problem
//! rather than an algorithmic one: the equivalent here is
//! [`map_collect`](crate::concept::parallel::map_collect), which returns
//! results in input order with no lock at all. The port is therefore both
//! parallel and *more* reproducible than the source, because a run's output
//! cannot depend on which thread finished first.
//!
//! # The contract
//!
//! **A parallel result must be bit-identical to the serial one.** A scientific
//! result that shifts with the thread count is not a result, and the two
//! `reduction(+:...)` clauses in the C++ are exactly where OpenMP gives that up:
//! floating-point addition is not associative, so an OpenMP sum depends on how
//! the runtime split the loop. Nothing here uses an unordered reduction.
//! [`sum_in_order`](crate::concept::parallel::sum_in_order) computes chunk
//! subtotals in parallel and then adds them in index order, so the answer is
//! fixed by the input alone.
//!
//! `tests/parallel_determinism.rs` enforces the contract by running the same
//! inputs at several thread counts and comparing bits.
//!
//! # Thread count
//!
//! TOPP tools take `-threads`, where 0 means every available core
//! (`src/cli/context.rs`). Pass that through [`Threads`](crate::concept::parallel::Threads) so a tool's setting
//! reaches the pool rather than being silently ignored.
//!
//! Without the `parallel` feature every entry point here runs serially and
//! returns the same values, which is the intended way to get a reproducible
//! single-threaded build for debugging.

/// How many worker threads a computation may use.
///
/// Mirrors the TOPP `-threads` parameter, where zero means every available
/// core. Construction is total: a negative count is clamped to one, because the
/// source treats any nonsensical value as single-threaded rather than failing a
/// tool that is otherwise ready to run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Threads(usize);

impl Threads {
    /// Every available core, as `-threads 0` requests.
    pub fn all() -> Self {
        Self(std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get))
    }

    /// A specific count, as TOPP passes it. Zero means [`Threads::all`].
    pub fn from_cli(threads: i64) -> Self {
        match threads {
            0 => Self::all(),
            n if n < 1 => Self(1),
            n => Self(usize::try_from(n).unwrap_or(1)),
        }
    }

    /// Exactly one worker: the reproducible baseline the determinism test uses.
    pub fn serial() -> Self {
        Self(1)
    }

    /// The worker count this policy permits, never zero.
    pub fn get(self) -> usize {
        self.0.max(1)
    }
}

impl Default for Threads {
    fn default() -> Self {
        Self::all()
    }
}

/// Map `f` over `items` in parallel, returning results in **input order**.
///
/// This replaces the source's `parallel for` plus `critical`-guarded push: the
/// ordering is a property of the call rather than of who finished first, so the
/// output is reproducible across runs and thread counts.
pub fn map_collect<T, U, F>(items: &[T], threads: Threads, f: F) -> Vec<U>
where
    T: Sync,
    U: Send,
    F: Fn(&T) -> U + Sync + Send,
{
    #[cfg(feature = "parallel")]
    {
        if threads.get() > 1 && items.len() > 1 {
            use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
            if let Ok(pool) = rayon::ThreadPoolBuilder::new()
                .num_threads(threads.get())
                .build()
            {
                return pool.install(|| items.par_iter().map(&f).collect());
            }
        }
    }
    let _ = threads;
    items.iter().map(f).collect()
}

/// Sum `values` with an order fixed by the input, not by the scheduler.
///
/// Chunk subtotals are computed in parallel and then added in ascending chunk
/// order, so the result depends only on `values` and `chunk`. An OpenMP
/// `reduction(+:x)` does not offer that: it sums partial results in whatever
/// order threads complete, and floating-point addition is not associative.
///
/// `chunk` is part of the contract, not a tuning knob — changing it changes the
/// summation tree and therefore the last bits. Callers that need to match a
/// previously published number must keep it fixed.
pub fn sum_in_order(values: &[f64], threads: Threads, chunk: usize) -> f64 {
    let chunk = chunk.max(1);
    #[cfg(feature = "parallel")]
    {
        if threads.get() > 1 && values.len() > chunk {
            use rayon::iter::ParallelIterator;
            use rayon::slice::ParallelSlice;
            if let Ok(pool) = rayon::ThreadPoolBuilder::new()
                .num_threads(threads.get())
                .build()
            {
                let subtotals: Vec<f64> = pool.install(|| {
                    values
                        .par_chunks(chunk)
                        .map(|part| part.iter().copied().sum::<f64>())
                        .collect()
                });
                return subtotals.iter().copied().sum();
            }
        }
    }
    let _ = threads;
    values
        .chunks(chunk)
        .map(|part| part.iter().copied().sum::<f64>())
        .sum()
}
