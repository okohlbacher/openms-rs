// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Matching elements of two sorted containers in one forward pass:
//! `DATASTRUCTURES/MatchedIterator.h`.
//!
//! For each element of a reference container the closest element of a target
//! container is searched; if none lies within the tolerance, the reference
//! element yields nothing and the next one is tested. A typical use is walking
//! the matching peaks of an experimental and a theoretical spectrum, in m/z
//! (ppm or Da), RT, or any user-defined metric. Both containers must be sorted
//! with respect to the metric. The pass is `O(n + m)`, against `O(n log m)`
//! for a binary search per reference element, and cache-friendly.
//!
//! # API mapping
//!
//! | Source (`MatchedIterator.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `MatchedIterator<CONT_T, TRAIT, CONST_T = true>` | [`MatchedIterator<'a, T, Tr>`], over slices |
//! | `MatchedIterator<CONT_T, TRAIT, false>` (mutable target) | [`MatchedIteratorMut<'a, T, Tr>`] |
//! | the five iterator typedefs, `CONT_IT`, `CONST_CONT_IT` | `Iterator` over `(ref_idx, tgt_idx)` pairs; slices and indices |
//! | `MatchedIterator(const CONT_T& ref, const CONT_T& target, float tolerance)` | [`MatchedIterator::new`] |
//! | `MatchedIterator(ref_begin, ref_end, tgt_begin, tgt_end, tolerance)` | [`MatchedIterator::new`] over sub-slices (`&ref[a..b]`) |
//! | `MatchedIterator()` ("do not use this for anything other than assigning") | `Default` |
//! | copy constructor, `operator=` | `Clone`, assignment |
//! | `operator==`, `operator!=` | `PartialEq` |
//! | `operator*`, `operator->` | [`MatchedIterator::target_element`] (and [`MatchedIteratorMut::target_element_mut`]) |
//! | `ref()` | [`MatchedIterator::reference_element`] |
//! | `refIdx()`, `tgtIdx()` | [`MatchedIterator::ref_idx`], [`MatchedIterator::tgt_idx`] |
//! | `operator++()` | [`MatchedIterator::advance`] |
//! | `operator++(int)` | [`Iterator::next`], which yields the current pair and then advances |
//! | `static end()` | [`MatchedIterator::end`], and [`MatchedIterator::is_end`] |
//! | `struct ValueTrait`, `PpmTrait`, `DaTrait` | [`ValueTrait`], [`PpmTrait`], [`DaTrait`], all implementing [`MatchTrait`] |
//! | `TRAIT::allowedTol`, `TRAIT::getDiffAbsolute` | [`MatchTrait::allowed_tol`], [`MatchTrait::diff_absolute`] |
//! | `.getMZ()`, which `PpmTrait` and `DaTrait` require of an element | [`MzCoordinate`] |
//!
//! # Arithmetic
//!
//! The source keeps the running best distance in a `float` initialised to
//! `std::numeric_limits<float>::max()`, compares it with each new distance in
//! the distance's own type (`float` for the m/z traits, the element type for
//! [`ValueTrait`]), narrows an improvement to `float`, and accepts a match when
//! that `float` is at most the `float` tolerance from `allowedTol`. This port
//! does the same: [`MatchTrait::diff_absolute`] returns the distance as the
//! source's type widened exactly to `f64`, and the running best is an `f32`.
//! A target element at the same distance as the previous one counts as
//! "getting worse": the pass stops there, so of two equally distant targets
//! the lower-positioned one is chosen, and a duplicate target can hide a later,
//! closer one. The target cursor never moves backwards across reference
//! elements, so a target element can match several reference elements.
//!
//! # Native differences
//!
//! - The source decrements the target iterator when the very first distance
//!   of a reference element is not below `float` max (a NaN or an infinite
//!   distance); at the first target element that decrement is undefined
//!   behaviour. Here the cursor stays at the first element.
//! - Advancing an end iterator is a checked error ([`MatchedIterator::advance`]);
//!   the source's `OPENMS_PRECONDITION` is compiled away in Release builds,
//!   where the increment reads past the reference container.
//! - Dereferencing is `Option`: `None` at the end and on a default-constructed
//!   iterator, where the source dereferences a singular iterator.
//! - The source's `CONST_T = false` variant initialises its non-const
//!   iterators from `const_iterator`s, which no standard container converts,
//!   so it does not compile for `std::vector` or `MSSpectrum`.
//!   [`MatchedIteratorMut`] provides the intended mutable access.

use std::marker::PhantomData;

use crate::{Error, Result};

/// The distance metric of a [`MatchedIterator`]: source `TRAIT` template
/// argument, "two functions which return some distance metrics and accept
/// elements of the container".
pub trait MatchTrait<T: ?Sized> {
    /// The largest distance that still matches `reference`: source
    /// `TRAIT::allowedTol(float tol, const T& elem_ref)`, whose unit is the
    /// trait's (ppm, Da, seconds, ...).
    fn allowed_tol(tolerance: f32, reference: &T) -> f32;

    /// The distance between a reference and a target element: source
    /// `TRAIT::getDiffAbsolute(elem_ref, elem_tgt)`, returned exactly in
    /// `f64` from the type the source computes it in (`float` or `double`).
    fn diff_absolute(reference: &T, target: &T) -> f64;
}

/// An element with an m/z coordinate: the `.getMZ()` member function that
/// [`PpmTrait`] and [`DaTrait`] require of container elements.
pub trait MzCoordinate {
    /// The element's m/z.
    fn mz_coordinate(&self) -> f64;
}

/// Distance computed directly on the value: source `ValueTrait`.
///
/// `allowedTol` is the tolerance itself and `getDiffAbsolute` is
/// `fabs(elem_ref - elem_tgt)` in the element type, `f32` or `f64`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ValueTrait;

impl MatchTrait<f64> for ValueTrait {
    fn allowed_tol(tolerance: f32, _reference: &f64) -> f32 {
        tolerance
    }
    fn diff_absolute(reference: &f64, target: &f64) -> f64 {
        (reference - target).abs()
    }
}

impl MatchTrait<f32> for ValueTrait {
    fn allowed_tol(tolerance: f32, _reference: &f32) -> f32 {
        tolerance
    }
    fn diff_absolute(reference: &f32, target: &f32) -> f64 {
        f64::from((reference - target).abs())
    }
}

/// A ppm distance in m/z: source `PpmTrait`.
///
/// `allowedTol` is `Math::ppmToMass(tol, (float)elem_ref.getMZ())`, i.e.
/// `(tol / 1e6f) * mz` in `float`; `getDiffAbsolute` is
/// `fabs(elem_ref.getMZ() - elem_tgt.getMZ())` in `double`, returned as
/// `float`. The window is taken from the reference element's m/z.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PpmTrait;

impl<T: MzCoordinate + ?Sized> MatchTrait<T> for PpmTrait {
    fn allowed_tol(tolerance: f32, reference: &T) -> f32 {
        (tolerance / 1e6_f32) * reference.mz_coordinate() as f32
    }
    fn diff_absolute(reference: &T, target: &T) -> f64 {
        f64::from((reference.mz_coordinate() - target.mz_coordinate()).abs() as f32)
    }
}

/// A Th/Da distance in m/z: source `DaTrait`.
///
/// `allowedTol` is the tolerance itself; `getDiffAbsolute` is as for
/// [`PpmTrait`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DaTrait;

impl<T: MzCoordinate + ?Sized> MatchTrait<T> for DaTrait {
    fn allowed_tol(tolerance: f32, _reference: &T) -> f32 {
        tolerance
    }
    fn diff_absolute(reference: &T, target: &T) -> f64 {
        f64::from((reference.mz_coordinate() - target.mz_coordinate()).abs() as f32)
    }
}

/// The position of a matched pass: source `it_ref_`, `it_tgt_` and `is_end_`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Cursor {
    ref_idx: usize,
    tgt_idx: usize,
    is_end: bool,
}

impl Cursor {
    /// Source constructor body: an empty target ends at once; otherwise the
    /// first match is searched.
    fn start<T, Tr: MatchTrait<T>>(reference: &[T], target: &[T], tolerance: f32) -> Self {
        let mut cursor = Self::default();
        if target.is_empty() {
            cursor.is_end = true;
        } else {
            cursor.advance_target::<T, Tr>(reference, target, tolerance);
        }
        cursor
    }

    /// Source `advanceTarget_()`.
    fn advance_target<T, Tr: MatchTrait<T>>(
        &mut self,
        reference: &[T],
        target: &[T],
        tolerance: f32,
    ) {
        while let Some(reference_element) = reference.get(self.ref_idx) {
            let max_dist = f64::from(Tr::allowed_tol(tolerance, reference_element));
            // Forward over the target while the distance improves.
            let mut diff = f32::MAX;
            while let Some(target_element) = target.get(self.tgt_idx) {
                let d = Tr::diff_absolute(reference_element, target_element);
                if f64::from(diff) > d {
                    diff = d as f32;
                } else {
                    // Getting worse (overshot): back to the previous element.
                    self.tgt_idx = self.tgt_idx.saturating_sub(1);
                    break;
                }
                self.tgt_idx += 1;
            }
            if self.tgt_idx == target.len() {
                // Reset to the last valid element.
                self.tgt_idx = target.len().saturating_sub(1);
            }
            if f64::from(diff) <= max_dist {
                return;
            }
            self.ref_idx += 1;
        }
        self.is_end = true;
    }
}

/// Iterator over the matching pairs of a reference and a target slice: source
/// `MatchedIterator<CONT_T, TRAIT, true>`.
///
/// It yields `(ref_idx, tgt_idx)` for every reference element that has a
/// target within the tolerance, choosing the closest target (the lower
/// position on a tie). [`Self::reference_element`], [`Self::target_element`],
/// [`Self::ref_idx`] and [`Self::tgt_idx`] read the current pair, as the
/// source's accessors do.
///
/// # Examples
///
/// The upstream class test's containers:
///
/// ```
/// use openms::data_structures::matched_iterator::{MatchedIterator, ValueTrait};
///
/// let reference = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 13.0];
/// let target = [-0.01, 2.5, 3.5, 7.0, 11.0];
/// let pairs: Vec<_> = MatchedIterator::<f64, ValueTrait>::new(&reference, &target, 0.5).collect();
/// assert_eq!(pairs, [(0, 0), (2, 1), (3, 1), (4, 2), (7, 3)]);
/// ```
#[derive(Debug)]
pub struct MatchedIterator<'a, T, Tr> {
    reference: &'a [T],
    target: &'a [T],
    cursor: Cursor,
    tolerance: f32,
    metric: PhantomData<Tr>,
}

impl<'a, T, Tr: MatchTrait<T>> MatchedIterator<'a, T, Tr> {
    /// The iterator positioned at the first match: source
    /// `MatchedIterator(ref, target, tolerance)` and, over sub-slices, the
    /// iterator-range constructor.
    ///
    /// `tolerance` is the maximal distance between a matching pair, in the
    /// unit of [`MatchTrait::allowed_tol`]. An empty target, or no match at
    /// all, gives an iterator that is already at its end.
    pub fn new(reference: &'a [T], target: &'a [T], tolerance: f32) -> Self {
        Self {
            reference,
            target,
            cursor: Cursor::start::<T, Tr>(reference, target, tolerance),
            tolerance,
            metric: PhantomData,
        }
    }

    /// Advances to the next valid pair: source `operator++()`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the iterator is already at its end: the
    /// source's `Exception::InvalidIterator` precondition ("Tried to advance
    /// beyond end iterator!"). The iterator is unchanged.
    pub fn advance(&mut self) -> Result<&mut Self> {
        if self.cursor.is_end {
            return Err(Error::InvalidValue(
                "tried to advance a MatchedIterator beyond its end".into(),
            ));
        }
        self.cursor.ref_idx += 1;
        self.cursor
            .advance_target::<T, Tr>(self.reference, self.target, self.tolerance);
        Ok(self)
    }
}

impl<'a, T, Tr> MatchedIterator<'a, T, Tr> {
    /// The end iterator: source `static end()`.
    ///
    /// Every iterator that has reached its end compares equal to it.
    pub fn end() -> Self {
        Self {
            reference: &[],
            target: &[],
            cursor: Cursor {
                is_end: true,
                ..Cursor::default()
            },
            tolerance: 0.0,
            metric: PhantomData,
        }
    }

    /// Whether the iterator is at its end, i.e. equal to [`Self::end`].
    pub fn is_end(&self) -> bool {
        self.cursor.is_end
    }

    /// The current reference element: source `ref()`. `None` at the end or on
    /// a default-constructed iterator.
    pub fn reference_element(&self) -> Option<&'a T> {
        if self.cursor.is_end {
            return None;
        }
        self.reference.get(self.cursor.ref_idx)
    }

    /// The current, matching target element: source `operator*` and
    /// `operator->`. `None` at the end or on a default-constructed iterator.
    pub fn target_element(&self) -> Option<&'a T> {
        if self.cursor.is_end {
            return None;
        }
        self.target.get(self.cursor.tgt_idx)
    }

    /// Index of the current element in the reference slice: source `refIdx()`.
    ///
    /// At the end of a pass it is the reference length, as the source's
    /// `ref.end()` position.
    pub fn ref_idx(&self) -> usize {
        self.cursor.ref_idx
    }

    /// Index of the current element in the target slice: source `tgtIdx()`.
    ///
    /// At the end of a pass it still points to a valid target element (the
    /// last one reached), as in the source.
    pub fn tgt_idx(&self) -> usize {
        self.cursor.tgt_idx
    }
}

impl<T, Tr> Default for MatchedIterator<'_, T, Tr> {
    /// Source `MatchedIterator()`: an iterator over nothing that is not the
    /// end iterator, to be assigned to. It yields nothing.
    fn default() -> Self {
        Self {
            reference: &[],
            target: &[],
            cursor: Cursor::default(),
            tolerance: 0.0,
            metric: PhantomData,
        }
    }
}

impl<T, Tr> Clone for MatchedIterator<'_, T, Tr> {
    fn clone(&self) -> Self {
        Self {
            reference: self.reference,
            target: self.target,
            cursor: self.cursor,
            tolerance: self.tolerance,
            metric: PhantomData,
        }
    }
}

impl<T, Tr> PartialEq for MatchedIterator<'_, T, Tr> {
    /// Source `operator==`: two end iterators are equal and an end iterator
    /// equals no other; otherwise both positions and both underlying slices
    /// (by identity, not by value) must agree.
    fn eq(&self, other: &Self) -> bool {
        if self.cursor.is_end || other.cursor.is_end {
            return self.cursor.is_end == other.cursor.is_end;
        }
        self.cursor == other.cursor
            && std::ptr::eq(self.reference, other.reference)
            && std::ptr::eq(self.target, other.target)
    }
}

impl<T, Tr: MatchTrait<T>> Iterator for MatchedIterator<'_, T, Tr> {
    type Item = (usize, usize);

    /// The current `(ref_idx, tgt_idx)` pair, after which the iterator
    /// advances: source `operator++(int)`, the post-increment.
    fn next(&mut self) -> Option<(usize, usize)> {
        if self.cursor.is_end || self.reference.get(self.cursor.ref_idx).is_none() {
            return None;
        }
        let pair = (self.cursor.ref_idx, self.cursor.tgt_idx);
        self.cursor.ref_idx += 1;
        self.cursor
            .advance_target::<T, Tr>(self.reference, self.target, self.tolerance);
        Some(pair)
    }
}

/// A matched pass with mutable access to the target: source
/// `MatchedIterator<CONT_T, TRAIT, false>`.
///
/// The pass is the same as [`MatchedIterator`]'s;
/// [`Self::target_element_mut`] gives the current target element mutably.
/// Rust's borrow rules make this a cursor rather than an `Iterator`: at most
/// one mutable target element is borrowed at a time.
///
/// # Examples
///
/// ```
/// use openms::data_structures::matched_iterator::{MatchedIteratorMut, ValueTrait};
///
/// let reference = [1.0, 5.0];
/// let mut target = [0.9, 3.0, 5.2];
/// let mut matches = MatchedIteratorMut::<f64, ValueTrait>::new(&reference, &mut target, 0.5);
/// while let Some(element) = matches.target_element_mut() {
///     *element = -*element;
///     if matches.advance().is_err() {
///         break;
///     }
/// }
/// assert_eq!(target, [-0.9, 3.0, -5.2]);
/// ```
#[derive(Debug)]
pub struct MatchedIteratorMut<'a, T, Tr> {
    reference: &'a [T],
    target: &'a mut [T],
    cursor: Cursor,
    tolerance: f32,
    metric: PhantomData<Tr>,
}

impl<'a, T, Tr: MatchTrait<T>> MatchedIteratorMut<'a, T, Tr> {
    /// The cursor positioned at the first match; see [`MatchedIterator::new`].
    pub fn new(reference: &'a [T], target: &'a mut [T], tolerance: f32) -> Self {
        let cursor = Cursor::start::<T, Tr>(reference, target, tolerance);
        Self {
            reference,
            target,
            cursor,
            tolerance,
            metric: PhantomData,
        }
    }

    /// Advances to the next valid pair; see [`MatchedIterator::advance`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the cursor is already at its end.
    pub fn advance(&mut self) -> Result<&mut Self> {
        if self.cursor.is_end {
            return Err(Error::InvalidValue(
                "tried to advance a MatchedIterator beyond its end".into(),
            ));
        }
        self.cursor.ref_idx += 1;
        self.cursor
            .advance_target::<T, Tr>(self.reference, self.target, self.tolerance);
        Ok(self)
    }

    /// Whether the cursor is at its end.
    pub fn is_end(&self) -> bool {
        self.cursor.is_end
    }

    /// The current reference element; see [`MatchedIterator::reference_element`].
    pub fn reference_element(&self) -> Option<&T> {
        if self.cursor.is_end {
            return None;
        }
        self.reference.get(self.cursor.ref_idx)
    }

    /// The current target element; see [`MatchedIterator::target_element`].
    pub fn target_element(&self) -> Option<&T> {
        if self.cursor.is_end {
            return None;
        }
        self.target.get(self.cursor.tgt_idx)
    }

    /// The current target element, mutably: the non-const source
    /// `operator*` and `operator->`.
    pub fn target_element_mut(&mut self) -> Option<&mut T> {
        if self.cursor.is_end {
            return None;
        }
        self.target.get_mut(self.cursor.tgt_idx)
    }

    /// Index into the reference slice; see [`MatchedIterator::ref_idx`].
    pub fn ref_idx(&self) -> usize {
        self.cursor.ref_idx
    }

    /// Index into the target slice; see [`MatchedIterator::tgt_idx`].
    pub fn tgt_idx(&self) -> usize {
        self.cursor.tgt_idx
    }
}
