// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The region quadtree `FeatureOverlapFilter` queries, ported from the
//! source's bundled `extern/Quadtree` (`Quadtree.h`, `Box.h`, `Vector2.h`;
//! Pierre Vigier's MIT-licensed quadtree, vendored unchanged by OpenMS).
//!
//! The tree stores values whose axis-aligned boxes are computed on demand by a
//! caller-supplied box function. A node holds up to 16 values
//! (`Quadtree::THRESHOLD`) before it splits into four children, and nodes at
//! depth 8 (`Quadtree::MAX_DEPTH`) never split. A value stays in the deepest
//! node whose quadrant contains its box strictly; a box that straddles a
//! centre line stays in the interior node. A query visits a node's own values
//! in insertion order before its children, and the children in the order
//! north-west, north-east, south-west, south-east. That order is observable:
//! `FeatureOverlapFilter` hands candidates to a state-changing callback in it,
//! so the port keeps every step of it.
//!
//! # Precision
//!
//! Boxes are `f32`, the source's `Box<float>`, and every derived quantity
//! (right and bottom edges, centres, child boxes) is computed with the same
//! `f32` operations in the same order as the source. Intersection is strict:
//! boxes that only touch do not intersect. Containment is inclusive.
//!
//! # Native differences
//!
//! - The source keeps the box function inside the tree. Here each operation
//!   takes it as an argument, so a caller can compute boxes from data it
//!   mutates between queries, as `FeatureOverlapFilter` does.
//! - The source's `assert`s (a value box inside its node, a query box
//!   intersecting the root) are not checked. A Release build compiles them out,
//!   and the port follows that build: a value outside the root is stored where
//!   the quadrant tests put it, and a query that misses the root still tests the
//!   root's own values. No input panics; NaN compares as IEEE 754 prescribes, as
//!   it does in the source.
//! - Removing a value that is not stored returns an error; the source asserts
//!   and then writes through the end iterator.
//! - The number of stored values and of reported intersections is bounded.
//! - The source is generic over the float type; only `float` is instantiated in
//!   OpenMS, and only `f32` is ported.
//!
//! See `docs/FEATURE_OVERLAP_FILTER_SUPPORT.md` for the API mapping and the
//! evidence.

use crate::{Error, Result};
use std::ops::{Add, Div};

/// A two-dimensional `f32` vector (source `quadtree::Vector2<float>`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vector2 {
    /// Horizontal component; m/z when the tree holds features.
    pub x: f32,
    /// Vertical component; retention time when the tree holds features.
    pub y: f32,
}

impl Vector2 {
    /// A vector with the given components.
    pub const fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }
}

impl Add for Vector2 {
    type Output = Self;

    /// Component-wise sum, the source `operator+`.
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl Div<f32> for Vector2 {
    type Output = Self;

    /// Division of both components by `rhs`, the source `operator/`.
    fn div(self, rhs: f32) -> Self {
        Self::new(self.x / rhs, self.y / rhs)
    }
}

/// An axis-aligned `f32` box given by its top-left corner and its size (source
/// `quadtree::Box<float>`).
///
/// The source requires a positive width and height and checks neither; this
/// type does not check them either, so the derived edges follow IEEE 754
/// arithmetic for any value.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct QuadBox {
    /// Left edge (smallest x).
    pub left: f32,
    /// Top edge (smallest y).
    pub top: f32,
    /// Extent along x.
    pub width: f32,
    /// Extent along y.
    pub height: f32,
}

impl QuadBox {
    /// A box with the given left and top edges, width and height.
    pub const fn new(left: f32, top: f32, width: f32, height: f32) -> Self {
        Self {
            left,
            top,
            width,
            height,
        }
    }

    /// A box at `position` (its top-left corner) with the given `size`.
    pub const fn from_position_size(position: Vector2, size: Vector2) -> Self {
        Self::new(position.x, position.y, size.x, size.y)
    }

    /// The right edge, `left + width` (source `getRight`).
    pub fn right(self) -> f32 {
        self.left + self.width
    }

    /// The bottom edge, `top + height` (source `getBottom`).
    pub fn bottom(self) -> f32 {
        self.top + self.height
    }

    /// The top-left corner (source `getTopLeft`).
    pub const fn top_left(self) -> Vector2 {
        Vector2::new(self.left, self.top)
    }

    /// The centre, `(left + width / 2, top + height / 2)` (source `getCenter`).
    pub fn center(self) -> Vector2 {
        Vector2::new(self.left + self.width / 2.0, self.top + self.height / 2.0)
    }

    /// The size `(width, height)` (source `getSize`).
    pub const fn size(self) -> Vector2 {
        Vector2::new(self.width, self.height)
    }

    /// Whether `other` lies inside this box, edges included (source `contains`).
    pub fn contains(self, other: Self) -> bool {
        self.left <= other.left
            && other.right() <= self.right()
            && self.top <= other.top
            && other.bottom() <= self.bottom()
    }

    /// Whether the boxes overlap with positive area (source `intersects`).
    ///
    /// Boxes that only share an edge or a corner do not intersect, and a box
    /// of zero width or height intersects nothing, not even an identical box.
    pub fn intersects(self, other: Self) -> bool {
        !(self.left >= other.right()
            || self.right() <= other.left
            || self.top >= other.bottom()
            || self.bottom() <= other.top)
    }
}

#[derive(Clone, Debug)]
struct Node<T> {
    children: Option<[usize; 4]>,
    values: Vec<T>,
}

impl<T> Node<T> {
    const fn new() -> Self {
        Self {
            children: None,
            values: Vec::new(),
        }
    }
}

/// A region quadtree over values with `f32` boxes (source
/// `quadtree::Quadtree<T, GetBox, Equal, float>`).
///
/// Nodes live in an arena; a node freed by a merge is reused by the next
/// split. Equality of values, the source's `Equal` functor, is [`PartialEq`].
#[derive(Clone, Debug)]
pub struct Quadtree<T> {
    bounds: QuadBox,
    nodes: Vec<Node<T>>,
    free: Vec<usize>,
    len: usize,
}

impl<T> Quadtree<T> {
    /// Values a leaf holds before it splits (source `Threshold`).
    pub const THRESHOLD: usize = 16;
    /// Depth at which leaves stop splitting; the root has depth 0 (source
    /// `MaxDepth`).
    pub const MAX_DEPTH: usize = 8;
    /// Values one tree stores; [`Self::add`] refuses more. The source is
    /// unbounded.
    pub const MAX_VALUES: usize = 10_000_000;
    /// Pairs one [`Self::find_all_intersections`] call reports before it fails.
    /// The source is unbounded.
    pub const MAX_INTERSECTIONS: usize = 10_000_000;

    /// An empty tree covering `bounds` (source constructor).
    pub fn new(bounds: QuadBox) -> Self {
        Self {
            bounds,
            nodes: vec![Node::new()],
            free: Vec::new(),
            len: 0,
        }
    }

    /// The box the tree covers (source `getBox`).
    pub const fn bounds(&self) -> QuadBox {
        self.bounds
    }

    /// The number of stored values. The source has no counterpart.
    pub const fn len(&self) -> usize {
        self.len
    }

    /// Whether the tree stores no value. The source has no counterpart.
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Store `value`, whose box `get_box` computes (source `add`).
    ///
    /// A leaf below [`Self::MAX_DEPTH`] that already holds
    /// [`Self::THRESHOLD`] values splits first. In an interior node the value
    /// descends into the quadrant that contains its box strictly and stays in
    /// the node when no quadrant does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the tree already holds
    /// [`Self::MAX_VALUES`] values; the tree is unchanged.
    pub fn add<G>(&mut self, value: T, get_box: &G) -> Result<()>
    where
        G: Fn(&T) -> QuadBox,
    {
        if self.len >= Self::MAX_VALUES {
            return Err(Error::InvalidValue(format!(
                "a quadtree holds at most {} values",
                Self::MAX_VALUES
            )));
        }
        let value_box = get_box(&value);
        let mut node = 0;
        let mut depth = 0;
        let mut node_box = self.bounds;
        loop {
            match self.nodes[node].children {
                None => {
                    if depth >= Self::MAX_DEPTH || self.nodes[node].values.len() < Self::THRESHOLD {
                        self.nodes[node].values.push(value);
                        break;
                    }
                    self.split(node, node_box, get_box);
                }
                Some(children) => match quadrant(node_box, value_box) {
                    Some(i) => {
                        node = children[i];
                        depth += 1;
                        node_box = child_box(node_box, i);
                    }
                    None => {
                        self.nodes[node].values.push(value);
                        break;
                    }
                },
            }
        }
        self.len += 1;
        Ok(())
    }

    /// Every stored value whose box intersects `query_box`, in the source's
    /// query order (source `query`).
    pub fn query<G>(&self, query_box: QuadBox, get_box: &G) -> Vec<T>
    where
        T: Clone,
        G: Fn(&T) -> QuadBox,
    {
        let mut values = Vec::new();
        self.query_into(query_box, get_box, &mut values);
        values
    }

    /// As [`Self::query`], writing into `values`, which is cleared first, so a
    /// caller can reuse one buffer across queries.
    pub fn query_into<G>(&self, query_box: QuadBox, get_box: &G, values: &mut Vec<T>)
    where
        T: Clone,
        G: Fn(&T) -> QuadBox,
    {
        values.clear();
        self.query_node(0, self.bounds, query_box, get_box, values);
    }

    fn query_node<G>(
        &self,
        node: usize,
        node_box: QuadBox,
        query_box: QuadBox,
        get_box: &G,
        values: &mut Vec<T>,
    ) where
        T: Clone,
        G: Fn(&T) -> QuadBox,
    {
        let current = &self.nodes[node];
        for value in &current.values {
            if query_box.intersects(get_box(value)) {
                values.push(value.clone());
            }
        }
        if let Some(children) = current.children {
            for (i, &child) in children.iter().enumerate() {
                let child_box = child_box(node_box, i);
                if query_box.intersects(child_box) {
                    self.query_node(child, child_box, query_box, get_box, values);
                }
            }
        }
    }

    /// Remove one stored value equal to `value` (source `remove`).
    ///
    /// The value is searched where [`Self::add`] would put it now, so its box
    /// must not have changed since it was added. It is swapped with the last
    /// value of its node and popped, which reorders that node as the source
    /// does. After a removal from a leaf, each ancestor, nearest first, merges
    /// its children into itself when all of them are leaves holding, with the
    /// ancestor, at most [`Self::THRESHOLD`] values; the walk stops at the first
    /// ancestor that does not merge.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the node searched holds no equal
    /// value; the tree is unchanged. The source asserts in a Debug build and
    /// writes through the end iterator in a Release build.
    pub fn remove<G>(&mut self, value: &T, get_box: &G) -> Result<()>
    where
        T: PartialEq,
        G: Fn(&T) -> QuadBox,
    {
        let value_box = get_box(value);
        let mut path = Vec::with_capacity(Self::MAX_DEPTH + 1);
        let mut node = 0;
        let mut node_box = self.bounds;
        let from_leaf = loop {
            match self.nodes[node].children {
                None => break true,
                Some(children) => match quadrant(node_box, value_box) {
                    Some(i) => {
                        path.push(node);
                        node = children[i];
                        node_box = child_box(node_box, i);
                    }
                    None => break false,
                },
            }
        };
        let position = self.nodes[node]
            .values
            .iter()
            .position(|stored| stored == value)
            .ok_or_else(|| {
                Error::InvalidValue("the quadtree does not hold the value to remove".into())
            })?;
        self.nodes[node].values.swap_remove(position);
        self.len -= 1;
        if from_leaf {
            while let Some(ancestor) = path.pop() {
                if !self.try_merge(ancestor) {
                    break;
                }
            }
        }
        Ok(())
    }

    /// Every pair of stored values whose boxes intersect, each pair once
    /// (source `findAllIntersections`).
    ///
    /// Within a node the pair is `(later, earlier)` in the node's order. A
    /// node's values are then paired with the values of each child's subtree,
    /// child by child, as `(node value, descendant)`, before the children are
    /// searched in turn.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when more than
    /// [`Self::MAX_INTERSECTIONS`] pairs intersect.
    pub fn find_all_intersections<G>(&self, get_box: &G) -> Result<Vec<(T, T)>>
    where
        T: Clone,
        G: Fn(&T) -> QuadBox,
    {
        let mut pairs = Vec::new();
        self.intersections_in(0, get_box, &mut pairs)?;
        Ok(pairs)
    }

    fn intersections_in<G>(&self, node: usize, get_box: &G, pairs: &mut Vec<(T, T)>) -> Result<()>
    where
        T: Clone,
        G: Fn(&T) -> QuadBox,
    {
        let current = &self.nodes[node];
        for (i, later) in current.values.iter().enumerate() {
            for earlier in &current.values[..i] {
                if get_box(later).intersects(get_box(earlier)) {
                    push_pair(pairs, later, earlier)?;
                }
            }
        }
        if let Some(children) = current.children {
            for &child in &children {
                for value in &current.values {
                    self.intersections_below(child, value, get_box, pairs)?;
                }
            }
            for &child in &children {
                self.intersections_in(child, get_box, pairs)?;
            }
        }
        Ok(())
    }

    fn intersections_below<G>(
        &self,
        node: usize,
        value: &T,
        get_box: &G,
        pairs: &mut Vec<(T, T)>,
    ) -> Result<()>
    where
        T: Clone,
        G: Fn(&T) -> QuadBox,
    {
        let current = &self.nodes[node];
        for other in &current.values {
            if get_box(value).intersects(get_box(other)) {
                push_pair(pairs, value, other)?;
            }
        }
        if let Some(children) = current.children {
            for &child in &children {
                self.intersections_below(child, value, get_box, pairs)?;
            }
        }
        Ok(())
    }

    fn split<G>(&mut self, node: usize, node_box: QuadBox, get_box: &G)
    where
        G: Fn(&T) -> QuadBox,
    {
        let children = [
            self.allocate(),
            self.allocate(),
            self.allocate(),
            self.allocate(),
        ];
        let values = std::mem::take(&mut self.nodes[node].values);
        let mut kept = Vec::new();
        for value in values {
            match quadrant(node_box, get_box(&value)) {
                Some(i) => self.nodes[children[i]].values.push(value),
                None => kept.push(value),
            }
        }
        self.nodes[node].values = kept;
        self.nodes[node].children = Some(children);
    }

    fn try_merge(&mut self, node: usize) -> bool {
        let Some(children) = self.nodes[node].children else {
            return false;
        };
        let mut count = self.nodes[node].values.len();
        for &child in &children {
            if self.nodes[child].children.is_some() {
                return false;
            }
            count += self.nodes[child].values.len();
        }
        if count > Self::THRESHOLD {
            return false;
        }
        for &child in &children {
            let values = std::mem::take(&mut self.nodes[child].values);
            self.nodes[node].values.extend(values);
            self.free.push(child);
        }
        self.nodes[node].children = None;
        true
    }

    fn allocate(&mut self) -> usize {
        match self.free.pop() {
            Some(index) => index,
            None => {
                self.nodes.push(Node::new());
                self.nodes.len() - 1
            }
        }
    }
}

fn push_pair<T: Clone>(pairs: &mut Vec<(T, T)>, first: &T, second: &T) -> Result<()> {
    if pairs.len() >= Quadtree::<T>::MAX_INTERSECTIONS {
        return Err(Error::InvalidValue(format!(
            "more than {} intersecting quadtree pairs",
            Quadtree::<T>::MAX_INTERSECTIONS
        )));
    }
    pairs.push((first.clone(), second.clone()));
    Ok(())
}

/// The child box `i` of `node_box`: 0 north-west, 1 north-east, 2 south-west,
/// 3 south-east (source `computeBox`).
fn child_box(node_box: QuadBox, i: usize) -> QuadBox {
    let origin = node_box.top_left();
    let child_size = node_box.size() / 2.0;
    match i {
        0 => QuadBox::from_position_size(origin, child_size),
        1 => {
            QuadBox::from_position_size(Vector2::new(origin.x + child_size.x, origin.y), child_size)
        }
        2 => {
            QuadBox::from_position_size(Vector2::new(origin.x, origin.y + child_size.y), child_size)
        }
        _ => QuadBox::from_position_size(origin + child_size, child_size),
    }
}

/// The quadrant of `node_box` that contains `value_box` strictly, or `None`
/// (source `getQuadrant`, which returns -1).
fn quadrant(node_box: QuadBox, value_box: QuadBox) -> Option<usize> {
    let center = node_box.center();
    if value_box.right() < center.x {
        if value_box.bottom() < center.y {
            Some(0)
        } else if value_box.top >= center.y {
            Some(2)
        } else {
            None
        }
    } else if value_box.left >= center.x {
        if value_box.bottom() < center.y {
            Some(1)
        } else if value_box.top >= center.y {
            Some(3)
        } else {
            None
        }
    } else {
        None
    }
}
