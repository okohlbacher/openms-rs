// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! `boost::heap::fibonacci_heap` as the Release build's Boost 1.89 behaves:
//! which element is on top among elements the comparison does not separate.
//!
//! The heap orders items by a caller-supplied strict "less" (the max-heap of
//! `boost::heap` with `std::less`: [`top`](FibonacciHeap::top) is an item no
//! other item is greater than, except after a lazy update, see below). Items
//! are identified by `usize` handles (the node numbers
//! [`push`](FibonacciHeap::push) returns), and their values live outside the
//! heap: the caller's `less(a, b)` reads them, so a value changed in place and
//! reported with [`update_lazy`](FibonacciHeap::update_lazy) is seen at once,
//! as `boost::heap`'s `*handle` reference is.
//!
//! What is reproduced, step by step as `boost/heap/fibonacci_heap.hpp` does it
//! (`stable` off, `constant_time_size` on):
//!
//! - `push`: the node goes to the front of the root list; it becomes the top
//!   when no top exists or `less(top, node)`;
//! - `pop`: the top is erased from the root list, its children are appended to
//!   the end of the root list with their parent cleared, and the heap is
//!   *consolidated*: the roots are walked in order, merged by rank with the
//!   greater root kept (`less(n, other)` swaps), the merged-in root appended to
//!   the end of the winner's child list, and the top is the last root that is
//!   not less than the running top;
//! - `update_lazy(handle)` (value already changed): a non-root node moves to
//!   the *front* of the root list (without resetting its mark or cutting its
//!   parent's ancestors), its children are appended to the end of the root
//!   list, and it becomes the top when `less(top, node)`. A node whose value
//!   *decreased* is not removed from the top: the next [`pop`](FibonacciHeap::pop)
//!   repairs it, but a `top` before that may be stale, and `QTClusterFinder`
//!   uses that top.
//!
//! The lists are doubly linked in an arena, so every operation is
//! O(1) except consolidation, as in the C++.
//!
//! Not ported, because no OpenMS code reaches them: `increase`, `decrease`,
//! the eager `update`, `erase` and `merge`, and with them the marks of
//! cascading cuts.

const NONE: usize = usize::MAX;

#[derive(Clone, Debug)]
struct Node {
    parent: usize,
    prev: usize,
    next: usize,
    /// Head and tail of the child list.
    first_child: usize,
    last_child: usize,
    child_count: usize,
}

/// Which list a node sits in: the roots or a parent's children.
#[derive(Clone, Copy)]
enum ListId {
    Roots,
    Children(usize),
}

/// A Fibonacci heap with `boost::heap::fibonacci_heap`'s structure and tie
/// behaviour; see the [module documentation](self).
#[derive(Clone, Debug, Default)]
pub struct FibonacciHeap {
    nodes: Vec<Node>,
    first_root: usize,
    last_root: usize,
    top: usize,
    size: usize,
}

impl FibonacciHeap {
    /// An empty heap.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            first_root: NONE,
            last_root: NONE,
            top: NONE,
            size: 0,
        }
    }

    /// Whether the heap is empty.
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// The number of items.
    pub fn len(&self) -> usize {
        self.size
    }

    /// The top item (`top()`), or `None` for an empty heap.
    pub fn top(&self) -> Option<usize> {
        if self.top == NONE {
            None
        } else {
            Some(self.top)
        }
    }

    /// Remove every item (`clear()`); the handles are invalid afterwards.
    pub fn clear(&mut self) {
        self.nodes.clear();
        self.first_root = NONE;
        self.last_root = NONE;
        self.top = NONE;
        self.size = 0;
    }

    fn list_bounds(&self, list: ListId) -> (usize, usize) {
        match list {
            ListId::Roots => (self.first_root, self.last_root),
            ListId::Children(parent) => (
                self.nodes[parent].first_child,
                self.nodes[parent].last_child,
            ),
        }
    }

    fn set_list_bounds(&mut self, list: ListId, first: usize, last: usize) {
        match list {
            ListId::Roots => {
                self.first_root = first;
                self.last_root = last;
            }
            ListId::Children(parent) => {
                self.nodes[parent].first_child = first;
                self.nodes[parent].last_child = last;
            }
        }
    }

    fn list_of(&self, node: usize) -> ListId {
        let parent = self.nodes[node].parent;
        if parent == NONE {
            ListId::Roots
        } else {
            ListId::Children(parent)
        }
    }

    /// Unlink `node` from its list.
    fn unlink(&mut self, node: usize) {
        let list = self.list_of(node);
        let (first, last) = self.list_bounds(list);
        let prev = self.nodes[node].prev;
        let next = self.nodes[node].next;
        if prev != NONE {
            self.nodes[prev].next = next;
        }
        if next != NONE {
            self.nodes[next].prev = prev;
        }
        let new_first = if first == node { next } else { first };
        let new_last = if last == node { prev } else { last };
        self.set_list_bounds(list, new_first, new_last);
        if let ListId::Children(parent) = list {
            self.nodes[parent].child_count -= 1;
        }
        self.nodes[node].prev = NONE;
        self.nodes[node].next = NONE;
    }

    /// Link an unlinked `node` at the front of `list`.
    fn push_front(&mut self, list: ListId, node: usize) {
        let (first, last) = self.list_bounds(list);
        self.nodes[node].prev = NONE;
        self.nodes[node].next = first;
        if first != NONE {
            self.nodes[first].prev = node;
        }
        let new_last = if last == NONE { node } else { last };
        self.set_list_bounds(list, node, new_last);
        if let ListId::Children(parent) = list {
            self.nodes[parent].child_count += 1;
        }
    }

    /// Link an unlinked `node` at the end of `list`.
    fn push_back(&mut self, list: ListId, node: usize) {
        let (first, last) = self.list_bounds(list);
        self.nodes[node].next = NONE;
        self.nodes[node].prev = last;
        if last != NONE {
            self.nodes[last].next = node;
        }
        let new_first = if first == NONE { node } else { first };
        self.set_list_bounds(list, new_first, node);
        if let ListId::Children(parent) = list {
            self.nodes[parent].child_count += 1;
        }
    }

    /// Insert an item (`push`) and return its handle.
    pub fn push(&mut self, less: &mut impl FnMut(usize, usize) -> bool) -> usize {
        let node = self.nodes.len();
        self.nodes.push(Node {
            parent: NONE,
            prev: NONE,
            next: NONE,
            first_child: NONE,
            last_child: NONE,
            child_count: 0,
        });
        self.size += 1;
        self.push_front(ListId::Roots, node);
        if self.top == NONE || less(self.top, node) {
            self.top = node;
        }
        node
    }

    /// `add_children_to_root`: clear the parent of every child of `node` and
    /// append them to the end of the root list.
    fn add_children_to_root(&mut self, node: usize) {
        let first = self.nodes[node].first_child;
        let last = self.nodes[node].last_child;
        if first == NONE {
            return;
        }
        let mut child = first;
        while child != NONE {
            self.nodes[child].parent = NONE;
            child = self.nodes[child].next;
        }
        // Splice the whole child list after the last root.
        self.nodes[first].prev = self.last_root;
        if self.last_root != NONE {
            self.nodes[self.last_root].next = first;
        } else {
            self.first_root = first;
        }
        self.last_root = last;
        self.nodes[node].first_child = NONE;
        self.nodes[node].last_child = NONE;
        self.nodes[node].child_count = 0;
    }

    /// `pop()`: remove the top item, which the caller identifies by the
    /// handle [`top`](Self::top) returned.
    pub fn pop(&mut self, less: &mut impl FnMut(usize, usize) -> bool) {
        let top = self.top;
        if top == NONE {
            return;
        }
        // `roots.erase(top)`: the top need not be a root only after a lazy
        // update, which always makes the node a root.
        self.unlink(top);
        self.finish_erase_or_pop(top, less);
    }

    fn finish_erase_or_pop(&mut self, erased: usize, less: &mut impl FnMut(usize, usize) -> bool) {
        self.add_children_to_root(erased);
        self.size -= 1;
        if self.size != 0 {
            self.consolidate(less);
        } else {
            self.top = NONE;
        }
    }

    /// `update_lazy(handle)`: the item's value changed in place.
    pub fn update_lazy(&mut self, node: usize, less: &mut impl FnMut(usize, usize) -> bool) {
        let parent = self.nodes[node].parent;
        if parent != NONE {
            // `roots.splice(roots.begin(), parent->children, n)`.
            self.unlink(node);
            self.nodes[node].parent = NONE;
            self.push_front(ListId::Roots, node);
        }
        self.add_children_to_root(node);
        if self.top != NONE && less(self.top, node) {
            self.top = node;
        }
    }

    /// `consolidate()`.
    fn consolidate(&mut self, less: &mut impl FnMut(usize, usize) -> bool) {
        if self.first_root == NONE {
            return;
        }
        let mut aux = [NONE; 64];
        let mut it = self.first_root;
        self.top = it;
        loop {
            let mut n = it;
            it = self.nodes[n].next;
            let mut node_rank = self.nodes[n].child_count;
            if aux[node_rank] == NONE {
                aux[node_rank] = n;
            } else {
                loop {
                    let mut other = aux[node_rank];
                    if less(n, other) {
                        std::mem::swap(&mut n, &mut other);
                    }
                    // `n->children.splice(n->children.end(), <other's list>, other)`.
                    self.unlink(other);
                    self.nodes[other].parent = n;
                    self.push_back(ListId::Children(n), other);
                    aux[node_rank] = NONE;
                    node_rank = self.nodes[n].child_count;
                    if aux[node_rank] == NONE {
                        break;
                    }
                }
                aux[node_rank] = n;
            }
            if !less(n, self.top) {
                self.top = n;
            }
            if it == NONE {
                break;
            }
        }
    }
}
