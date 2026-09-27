// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! General helper and condition macros: `CONCEPT/Macros.h`.
//!
//! The header declares no type and no function, only four preprocessor macros.
//! Each has a Rust counterpart here or in the standard library.
//!
//! # API mapping
//!
//! | Source (`Macros.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `STRINGIFY(a)` (`#a`) | [`std::stringify!`] |
//! | `OPENMS_THREAD_CRITICAL(name)` (`omp critical (name)`) | [`thread_critical`] |
//! | `OPENMS_PRECONDITION(condition, message)` | [`openms_precondition!`](crate::openms_precondition) |
//! | `OPENMS_POSTCONDITION(condition, message)` | [`openms_postcondition!`](crate::openms_postcondition) |
//!
//! # Critical sections
//!
//! `OPENMS_THREAD_CRITICAL(name)` expands to `#pragma omp critical (name)` in an
//! OpenMP build and to nothing otherwise. An OpenMP critical construct lets one
//! thread at a time execute its block; all critical constructs with the same
//! name exclude each other, across the whole program, while differently named
//! ones do not. [`thread_critical`] gives a closure exactly that guarantee: a
//! process-wide mutual-exclusion lock per name. In a serial program it only
//! costs an uncontended lock, where the non-OpenMP expansion costs nothing;
//! the observable behaviour is the same.
//!
//! The port's parallel algorithms do not use it: `concept::parallel` replaces
//! the source's critical sections, which mostly serialise pushes into a shared
//! container, with per-worker results merged in order, so that a parallel
//! result is bit-identical to the serial one. [`thread_critical`] is for the
//! remaining case of genuinely shared state, such as one output stream.
//!
//! # Conditions
//!
//! The source's condition macros are "enabled if debug info is enabled and
//! optimization is disabled"; otherwise "they are replaced by an empty string,
//! so they won't cost any performance". Enabled, they throw
//! `Exception::Precondition` or `Exception::Postcondition` with the message
//! `"<condition> <message>"`. The pinned Release build the port is compared
//! against compiles both away (`OPENMS_ASSERTIONS` is off, see
//! `docs/BENCHMARKS.md`), and the port never turns a Debug-only precondition
//! into a Rust refusal. The Rust counterparts are therefore
//! [`debug_assert!`]: enabled with `debug_assertions`, compiled away in an
//! optimised release build, and reporting the condition's source text and the
//! message, as the source's `# condition " " # message` does (the source also
//! stringifies the message, so a string literal keeps its quotes there). A failure panics rather than returning an error, because the
//! condition states a programming error, not a property of input data: no
//! crate function uses them on untrusted input.

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, OnceLock, PoisonError};

/// The lock of every critical-section name used so far.
type CriticalSectionLocks = BTreeMap<&'static str, Arc<Mutex<()>>>;

/// One lock per critical-section name, created on first use.
static CRITICAL_SECTIONS: OnceLock<Mutex<CriticalSectionLocks>> = OnceLock::new();

/// The lock for `name`, created on first use.
fn critical_section_lock(name: &'static str) -> Arc<Mutex<()>> {
    let registry = CRITICAL_SECTIONS.get_or_init(|| Mutex::new(BTreeMap::new()));
    // The registry only maps names to locks, so a panic while it was held
    // cannot have left it inconsistent; recover it rather than propagate.
    let mut locks = registry.lock().unwrap_or_else(PoisonError::into_inner);
    Arc::clone(locks.entry(name).or_default())
}

/// Runs `body` inside the critical section `name`: source
/// `OPENMS_THREAD_CRITICAL(name)`, that is `#pragma omp critical (name)`.
///
/// At most one thread at a time runs a `body` under a given `name`, program
/// wide; bodies under different names may run concurrently. The source's
/// unnamed `omp critical` is one more name; pass a name of your own instead.
/// Returns what `body` returns.
///
/// A panic inside `body` releases the section; the next caller still enters
/// it, since the lock guards no data of its own.
///
/// # Deadlock
///
/// Entering a section from inside a section of the same name deadlocks, as
/// nesting same-named OpenMP critical constructs does. Nesting differently
/// named sections is allowed, as in OpenMP, but two threads nesting the same
/// two names in opposite orders can deadlock, as in OpenMP.
///
/// # Examples
///
/// ```
/// use openms::concept::macros::thread_critical;
/// use std::sync::atomic::{AtomicUsize, Ordering};
///
/// static INSIDE: AtomicUsize = AtomicUsize::new(0);
/// std::thread::scope(|scope| {
///     for _ in 0..4 {
///         scope.spawn(|| {
///             thread_critical("doc_example", || {
///                 assert_eq!(INSIDE.fetch_add(1, Ordering::SeqCst), 0);
///                 INSIDE.fetch_sub(1, Ordering::SeqCst);
///             })
///         });
///     }
/// });
/// ```
pub fn thread_critical<R>(name: &'static str, body: impl FnOnce() -> R) -> R {
    let lock = critical_section_lock(name);
    let _guard = lock.lock().unwrap_or_else(PoisonError::into_inner);
    body()
}

/// Checks a precondition in builds with debug assertions: source
/// `OPENMS_PRECONDITION(condition, message)`.
///
/// The source throws `Exception::Precondition` with the condition's text and
/// the message when assertions are enabled and expands to nothing otherwise;
/// this is [`debug_assert!`] reporting the same two parts. See the
/// `concept::macros` module notes for why it panics rather than returning an
/// error.
///
/// # Examples
///
/// ```
/// let values = [1.0, 2.0];
/// openms::openms_precondition!(!values.is_empty(), "values must not be empty");
/// ```
#[macro_export]
macro_rules! openms_precondition {
    ($condition:expr, $message:expr $(,)?) => {
        debug_assert!(
            $condition,
            "Precondition failed: {} {}",
            stringify!($condition),
            $message
        )
    };
}

/// Checks a postcondition in builds with debug assertions: source
/// `OPENMS_POSTCONDITION(condition, message)`.
///
/// As [`openms_precondition!`](crate::openms_precondition), for the source's
/// `Exception::Postcondition`.
///
/// # Examples
///
/// ```
/// let sum: i32 = [1, 2, 3].iter().sum();
/// openms::openms_postcondition!(sum == 6, "the sum of 1, 2 and 3");
/// ```
#[macro_export]
macro_rules! openms_postcondition {
    ($condition:expr, $message:expr $(,)?) => {
        debug_assert!(
            $condition,
            "Postcondition failed: {} {}",
            stringify!($condition),
            $message
        )
    };
}
