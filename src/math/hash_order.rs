// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The iteration orders of the C++ hash containers that OpenMS algorithms
//! expose through their results.
//!
//! The order in which `std::unordered_map`, `std::unordered_set` and
//! `boost::unordered_map` visit their elements is unspecified, but it is
//! deterministic for a given implementation, hash function, insertion and
//! erasure history, and some OpenMS algorithms depend on it: clusters of equal
//! quality are extracted in the order their grid cells were visited,
//! floating-point distances are summed in the order a hash map lists them, and
//! the elements of a cluster are listed in a hash map's order. This module
//! reproduces two implementations, as the Linux x86_64 Release build uses
//! them:
//!
//! - [`StdHashMap`](crate::math::hash_order::StdHashMap): libstdc++ 14's `std::unordered_map<size_t, V>` (a
//!   `std::unordered_set<size_t>` is the same table with no value). The
//!   identity hash and the prime rehash policy (`_Prime_rehash_policy`, load
//!   factor 1, growth factor 2, the `__fast_bkt` table and the `__prime_list`
//!   from its 7th entry) place `key % bucket_count`; nodes sit in one singly
//!   linked list, a node entering an empty bucket at the head of the list and
//!   one entering a non-empty bucket in front of that bucket's first node;
//!   `clear` keeps the bucket count and the policy's resize threshold.
//! - [`BoostTable`](crate::math::hash_order::BoostTable): Boost 1.89's node-based `boost::unordered_map` and
//!   `boost::unordered_multimap` (`boost/unordered/detail/fca.hpp`,
//!   `implementation.hpp`): a prime bucket count (`prime_fmod_size`), the
//!   bucket of a hash `h` is `(u32(h) + u32(h >> 32)) % prime`, a new node goes
//!   to the head of its bucket's chain (an equal key's node is followed
//!   directly by the new equal-key node of a multimap), and iteration walks
//!   the *groups* of 64 buckets in reverse order of their activation, then
//!   the group that holds the end sentinel last, each group's buckets in
//!   ascending index and each bucket's chain in order.
//!
//! Both are validated against the executed C++ (`../oracle/M2-hashorder`,
//! `tests/hash_order_oracle.rs`).
//!
//! The reproductions model the *order*, not a general-purpose container: they
//! keep keys and (for `StdHashMap`) values, and lookups go through the
//! bucket structure like the original, so a lookup is as fast as a hash
//! lookup.

/// `__prime_list` of libstdc++ from its 7th entry (17), up to the largest
/// bucket count a table of at most 2^26 elements can reach.
static PRIMES: [usize; 197] = [
    17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 103, 109, 113, 127,
    137, 139, 149, 157, 167, 179, 193, 199, 211, 227, 241, 257, 277, 293, 313, 337, 359, 383, 409,
    439, 467, 503, 541, 577, 619, 661, 709, 761, 823, 887, 953, 1031, 1109, 1193, 1289, 1381, 1493,
    1613, 1741, 1879, 2029, 2179, 2357, 2549, 2753, 2971, 3209, 3469, 3739, 4027, 4349, 4703, 5087,
    5503, 5953, 6427, 6949, 7517, 8123, 8783, 9497, 10273, 11113, 12011, 12983, 14033, 15173,
    16411, 17749, 19183, 20753, 22447, 24281, 26267, 28411, 30727, 33223, 35933, 38873, 42043,
    45481, 49201, 53201, 57557, 62233, 67307, 72817, 78779, 85229, 92203, 99733, 107897, 116731,
    126271, 136607, 147793, 159871, 172933, 187091, 202409, 218971, 236897, 256279, 277261, 299951,
    324503, 351061, 379787, 410857, 444487, 480881, 520241, 562841, 608903, 658753, 712697, 771049,
    834181, 902483, 976369, 1056323, 1142821, 1236397, 1337629, 1447153, 1565659, 1693859, 1832561,
    1982627, 2144977, 2320627, 2510653, 2716249, 2938679, 3179303, 3439651, 3721303, 4026031,
    4355707, 4712381, 5098259, 5515729, 5967347, 6456007, 6984629, 7556579, 8175383, 8844859,
    9569143, 10352717, 11200489, 12117689, 13109983, 14183539, 15345007, 16601593, 17961079,
    19431899, 21023161, 22744717, 24607243, 26622317, 28802401, 31160981, 33712729, 36473443,
    39460231, 42691603, 46187573, 49969847, 54061849, 58488943, 63278561, 68460391,
];

/// `__fast_bkt` of `_Prime_rehash_policy::_M_next_bkt`.
const FAST_BKT: [usize; 14] = [2, 2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11, 13, 13];

/// No node / empty bucket.
const NONE: usize = usize::MAX;
/// The `_M_before_begin` pseudo-node of a libstdc++ table.
const BEFORE_BEGIN: usize = usize::MAX - 1;

/// A libstdc++ node.
#[derive(Clone, Debug)]
struct StdNode<V> {
    key: usize,
    value: Option<V>,
    next: usize,
}

/// The iteration order of libstdc++'s `std::unordered_map<size_t, V>` with
/// the default (identity) hash, the default prime rehash policy and unique
/// keys; see the [module documentation](self).
///
/// [`iter`](Self::iter) lists the elements in the order the C++ container's
/// iterators do after the same sequence of [`insert`](Self::insert),
/// [`get_or_insert_with`](Self::get_or_insert_with), [`remove`](Self::remove)
/// and [`clear`](Self::clear) calls.
#[derive(Clone, Debug)]
pub struct StdHashMap<V> {
    nodes: Vec<StdNode<V>>,
    free: Vec<usize>,
    /// `_M_before_begin._M_nxt`: the first node of the list.
    before_begin_next: usize,
    /// For each bucket the node *before* its first node (`BEFORE_BEGIN` for
    /// the list head), or `NONE`.
    buckets: Vec<usize>,
    bucket_count: usize,
    element_count: usize,
    /// `_Prime_rehash_policy::_M_next_resize`.
    next_resize: usize,
}

impl<V> Default for StdHashMap<V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<V> StdHashMap<V> {
    /// An empty table, as a default-constructed `std::unordered_map`: one
    /// bucket, no allocated nodes.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            free: Vec::new(),
            before_begin_next: NONE,
            buckets: vec![NONE],
            bucket_count: 1,
            element_count: 0,
            next_resize: 0,
        }
    }

    /// The number of elements (`size()`).
    pub fn len(&self) -> usize {
        self.element_count
    }

    /// Whether the table is empty (`empty()`).
    pub fn is_empty(&self) -> bool {
        self.element_count == 0
    }

    /// The bucket count (`bucket_count()`).
    pub fn bucket_count(&self) -> usize {
        self.bucket_count
    }

    fn next_of(&self, node: usize) -> usize {
        if node == BEFORE_BEGIN {
            self.before_begin_next
        } else {
            self.nodes[node].next
        }
    }

    fn set_next(&mut self, node: usize, next: usize) {
        if node == BEFORE_BEGIN {
            self.before_begin_next = next;
        } else {
            self.nodes[node].next = next;
        }
    }

    fn bucket_of(&self, key: usize) -> usize {
        key % self.bucket_count
    }

    /// `_M_find_before_node`: the node before the node with `key` in bucket
    /// `bkt`, or `NONE`.
    fn find_before(&self, bkt: usize, key: usize) -> usize {
        let mut prev = self.buckets[bkt];
        if prev == NONE {
            return NONE;
        }
        let mut p = self.next_of(prev);
        while p != NONE {
            if self.nodes[p].key == key {
                return prev;
            }
            let next = self.nodes[p].next;
            if next == NONE || self.bucket_of(self.nodes[next].key) != bkt {
                return NONE;
            }
            prev = p;
            p = next;
        }
        NONE
    }

    fn find_node(&self, key: usize) -> usize {
        let before = self.find_before(self.bucket_of(key), key);
        if before == NONE {
            NONE
        } else {
            self.next_of(before)
        }
    }

    /// The value for `key` (`find`).
    pub fn get(&self, key: usize) -> Option<&V> {
        let node = self.find_node(key);
        if node == NONE {
            None
        } else {
            self.nodes[node].value.as_ref()
        }
    }

    /// The value for `key`, mutable (`find`).
    pub fn get_mut(&mut self, key: usize) -> Option<&mut V> {
        let node = self.find_node(key);
        if node == NONE {
            None
        } else {
            self.nodes[node].value.as_mut()
        }
    }

    /// Whether `key` is present (`contains`).
    pub fn contains_key(&self, key: usize) -> bool {
        self.find_node(key) != NONE
    }

    /// `_Prime_rehash_policy::_M_next_bkt`.
    fn next_bkt(&mut self, n: usize) -> usize {
        if n < FAST_BKT.len() {
            if n == 0 {
                return 1;
            }
            self.next_resize = FAST_BKT[n];
            return FAST_BKT[n];
        }
        let index = PRIMES.partition_point(|&prime| prime < n);
        let prime = PRIMES[index.min(PRIMES.len() - 1)];
        self.next_resize = prime;
        prime
    }

    /// `_Prime_rehash_policy::_M_need_rehash` for `n_ins` insertions.
    fn need_rehash(&mut self, n_ins: usize) -> Option<usize> {
        let n_elt = self.element_count;
        if n_elt + n_ins > self.next_resize {
            // With nothing allocated yet (`_M_next_resize` 0) start at 11.
            let floor = if self.next_resize == 0 { 11 } else { 0 };
            let min_bkts = (n_elt + n_ins).max(floor);
            if min_bkts >= self.bucket_count {
                let wanted = (min_bkts + 1).max(self.bucket_count * 2);
                return Some(self.next_bkt(wanted));
            }
            self.next_resize = self.bucket_count;
            return None;
        }
        None
    }

    /// `_M_rehash(count, true_type)`: relink every node into `count` buckets.
    fn rehash(&mut self, count: usize) {
        let mut new_buckets = vec![NONE; count];
        let mut p = self.before_begin_next;
        self.before_begin_next = NONE;
        let mut bbegin_bkt = 0usize;
        while p != NONE {
            let next = self.nodes[p].next;
            let bkt = self.nodes[p].key % count;
            if new_buckets[bkt] == NONE {
                self.nodes[p].next = self.before_begin_next;
                self.before_begin_next = p;
                new_buckets[bkt] = BEFORE_BEGIN;
                if self.nodes[p].next != NONE {
                    new_buckets[bbegin_bkt] = p;
                }
                bbegin_bkt = bkt;
            } else {
                let before = new_buckets[bkt];
                let after = self.next_of(before);
                self.nodes[p].next = after;
                self.set_next(before, p);
            }
            p = next;
        }
        self.buckets = new_buckets;
        self.bucket_count = count;
    }

    /// `_M_insert_bucket_begin`.
    fn insert_bucket_begin(&mut self, bkt: usize, node: usize) {
        if self.buckets[bkt] != NONE {
            let before = self.buckets[bkt];
            let after = self.next_of(before);
            self.nodes[node].next = after;
            self.set_next(before, node);
        } else {
            self.nodes[node].next = self.before_begin_next;
            self.before_begin_next = node;
            let next = self.nodes[node].next;
            if next != NONE {
                let next_bkt = self.bucket_of(self.nodes[next].key);
                self.buckets[next_bkt] = node;
            }
            self.buckets[bkt] = BEFORE_BEGIN;
        }
    }

    fn allocate(&mut self, key: usize, value: V) -> usize {
        let node = StdNode {
            key,
            value: Some(value),
            next: NONE,
        };
        if let Some(slot) = self.free.pop() {
            self.nodes[slot] = node;
            slot
        } else {
            self.nodes.push(node);
            self.nodes.len() - 1
        }
    }

    /// `operator[](key) = value`: assign when present, otherwise insert a new
    /// node. Returns the replaced value.
    pub fn insert(&mut self, key: usize, value: V) -> Option<V> {
        let existing = self.find_node(key);
        if existing != NONE {
            return self.nodes[existing].value.replace(value);
        }
        self.insert_new(key, value);
        None
    }

    fn insert_new(&mut self, key: usize, value: V) -> usize {
        let node = self.allocate(key, value);
        if let Some(count) = self.need_rehash(1) {
            self.rehash(count);
        }
        let bkt = self.bucket_of(key);
        self.insert_bucket_begin(bkt, node);
        self.element_count += 1;
        node
    }

    /// `operator[](key)`: the value for `key`, default-inserting `make()`
    /// when absent.
    pub fn get_or_insert_with(&mut self, key: usize, make: impl FnOnce() -> V) -> &mut V {
        let mut node = self.find_node(key);
        if node == NONE {
            node = self.insert_new(key, make());
        }
        self.nodes[node]
            .value
            .as_mut()
            .expect("a live node holds a value")
    }

    /// `erase(key)`: the removed value, if any.
    pub fn remove(&mut self, key: usize) -> Option<V> {
        let bkt = self.bucket_of(key);
        let prev = self.find_before(bkt, key);
        if prev == NONE {
            return None;
        }
        let node = self.next_of(prev);
        let next = self.nodes[node].next;
        // `_M_erase(bkt, prev, node)`.
        if prev == self.buckets[bkt] {
            // `_M_remove_bucket_begin`.
            if next == NONE {
                self.buckets[bkt] = NONE;
            } else {
                let next_bkt = self.bucket_of(self.nodes[next].key);
                if next_bkt != bkt {
                    self.buckets[next_bkt] = self.buckets[bkt];
                    self.buckets[bkt] = NONE;
                }
            }
        } else if next != NONE {
            let next_bkt = self.bucket_of(self.nodes[next].key);
            if next_bkt != bkt {
                self.buckets[next_bkt] = prev;
            }
        }
        self.set_next(prev, next);
        let value = self.nodes[node].value.take();
        self.free.push(node);
        self.element_count -= 1;
        value
    }

    /// `clear()`: removes every element and keeps the bucket count.
    pub fn clear(&mut self) {
        self.nodes.clear();
        self.free.clear();
        self.before_begin_next = NONE;
        for bucket in &mut self.buckets {
            *bucket = NONE;
        }
        self.element_count = 0;
    }

    /// The elements in the container's iteration order.
    pub fn iter(&self) -> StdHashMapIter<'_, V> {
        StdHashMapIter {
            map: self,
            current: self.before_begin_next,
        }
    }

    /// The keys in iteration order.
    pub fn keys(&self) -> impl Iterator<Item = usize> + '_ {
        self.iter().map(|(key, _)| key)
    }

    /// The values in iteration order.
    pub fn values(&self) -> impl Iterator<Item = &V> + '_ {
        self.iter().map(|(_, value)| value)
    }

    /// Visit every value mutably in iteration order.
    pub fn for_each_value_mut(&mut self, mut visit: impl FnMut(usize, &mut V)) {
        let mut current = self.before_begin_next;
        while current != NONE {
            let next = self.nodes[current].next;
            let key = self.nodes[current].key;
            if let Some(value) = self.nodes[current].value.as_mut() {
                visit(key, value);
            }
            current = next;
        }
    }
}

/// Iterator of [`StdHashMap::iter`].
pub struct StdHashMapIter<'a, V> {
    map: &'a StdHashMap<V>,
    current: usize,
}

impl<'a, V> Iterator for StdHashMapIter<'a, V> {
    type Item = (usize, &'a V);

    fn next(&mut self) -> Option<Self::Item> {
        if self.current == NONE {
            return None;
        }
        let node = &self.map.nodes[self.current];
        self.current = node.next;
        node.value.as_ref().map(|value| (node.key, value))
    }
}

/// `prime_fmod_size::sizes` of Boost.Unordered (those below 2^32).
static BOOST_SIZES: [usize; 29] = [
    13, 29, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241,
    786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611,
    402653189, 805306457, 1610612741, 3221225473,
];

/// A Boost.Unordered node.
#[derive(Clone, Debug)]
struct BoostNode<K> {
    key: K,
    hash: u64,
    next: usize,
    live: bool,
}

/// The iteration order of Boost 1.89's node-based `boost::unordered_map` and
/// `boost::unordered_multimap`; see the [module documentation](self).
///
/// Nodes are numbered in insertion order and never renumbered; the caller
/// keeps its payload in a parallel vector indexed by node number and asks
/// [`order`](Self::order) for the iteration order of the live nodes. The hash
/// of a key is supplied by the caller (`boost::hash` of the key type), and
/// keys compare with `==` as the container's `key_equal` does.
#[derive(Clone, Debug)]
pub struct BoostTable<K> {
    nodes: Vec<BoostNode<K>>,
    /// First node of each bucket's chain, or `NONE`.
    buckets: Vec<usize>,
    /// One bit per bucket of a group of 64, plus the end sentinel's bit.
    bitmask: Vec<u64>,
    /// The linked groups other than the sentinel's, in order of activation;
    /// iteration visits them in reverse.
    activation: Vec<usize>,
    size: usize,
    max_load: usize,
}

impl<K: PartialEq> Default for BoostTable<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: PartialEq> BoostTable<K> {
    /// An empty table (no buckets allocated, as a default-constructed
    /// container).
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            buckets: Vec::new(),
            bitmask: Vec::new(),
            activation: Vec::new(),
            size: 0,
            max_load: 0,
        }
    }

    /// The number of live nodes.
    pub fn len(&self) -> usize {
        self.size
    }

    /// Whether the table is empty.
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// The bucket count (`bucket_count()`).
    pub fn bucket_count(&self) -> usize {
        self.buckets.len()
    }

    /// `prime_fmod_size::position`.
    fn position(&self, hash: u64) -> usize {
        let size = self.buckets.len() as u64;
        if size < (1u64 << 32) {
            let folded = (hash as u32).wrapping_add((hash >> 32) as u32);
            (u64::from(folded) % size) as usize
        } else {
            (hash % size) as usize
        }
    }

    /// `bucket_count_for`.
    fn bucket_count_for(num_buckets: usize) -> usize {
        if num_buckets == 0 {
            return 0;
        }
        for &size in &BOOST_SIZES {
            if size >= num_buckets {
                return size;
            }
        }
        BOOST_SIZES[BOOST_SIZES.len() - 1]
    }

    /// `insert_node` / `insert_node_hint` with `append_bucket_group`.
    fn link(&mut self, bucket: usize, node: usize, hint: usize) {
        if self.buckets[bucket] == NONE {
            let group = bucket / 64;
            let last = (self.buckets.len()) / 64;
            if group != last && self.bitmask[group] == 0 {
                self.activation.push(group);
            }
            self.bitmask[group] |= 1u64 << (bucket % 64);
        }
        if hint != NONE {
            self.nodes[node].next = self.nodes[hint].next;
            self.nodes[hint].next = node;
        } else {
            self.nodes[node].next = self.buckets[bucket];
            self.buckets[bucket] = node;
        }
    }

    /// `rehash_impl`: a new bucket array, nodes moved bucket by bucket in
    /// ascending bucket index and chain order, each to the head of its new
    /// chain.
    fn rehash_impl(&mut self, num_buckets: usize) {
        let old_buckets = std::mem::take(&mut self.buckets);
        self.buckets = vec![NONE; num_buckets];
        let groups = num_buckets / 64 + 1;
        self.bitmask = vec![0; groups];
        self.bitmask[num_buckets / 64] = 1u64 << (num_buckets % 64);
        self.activation.clear();
        for head in old_buckets {
            let mut p = head;
            while p != NONE {
                let next = self.nodes[p].next;
                let bucket = self.position(self.nodes[p].hash);
                self.link(bucket, p, NONE);
                p = next;
            }
        }
        self.max_load = num_buckets;
    }

    /// `table::reserve(n)` for `mlf` 1.
    fn reserve(&mut self, num_elements: usize) {
        let wanted = Self::bucket_count_for(self.size.max(num_elements));
        if wanted != self.buckets.len() {
            self.rehash_impl(wanted);
        }
    }

    fn find_in_chain(&self, bucket: usize, key: &K) -> usize {
        let mut p = self.buckets[bucket];
        while p != NONE {
            if self.nodes[p].key == *key {
                return p;
            }
            p = self.nodes[p].next;
        }
        NONE
    }

    /// The first node whose key equals `key` (`find`).
    pub fn find(&self, key: &K, hash: u64) -> Option<usize> {
        if self.buckets.is_empty() {
            return None;
        }
        let node = self.find_in_chain(self.position(hash), key);
        if node == NONE { None } else { Some(node) }
    }

    /// Insert a node with a unique key (`emplace_unique` / `operator[]`):
    /// the existing node and `false` when the key is present.
    pub fn insert_unique(&mut self, key: K, hash: u64) -> (usize, bool) {
        if let Some(node) = self.find(&key, hash) {
            return (node, false);
        }
        let node = self.nodes.len();
        self.nodes.push(BoostNode {
            key,
            hash,
            next: NONE,
            live: true,
        });
        if self.size + 1 > self.max_load {
            self.reserve(self.size + 1);
        }
        let bucket = self.position(hash);
        self.link(bucket, node, NONE);
        self.size += 1;
        (node, true)
    }

    /// Insert a node into a multimap (`emplace_equiv`): after the first node
    /// with an equal key, or at the head of its bucket.
    pub fn insert_equiv(&mut self, key: K, hash: u64) -> usize {
        let hint = self.find(&key, hash).map_or(NONE, |node| node);
        let node = self.nodes.len();
        self.nodes.push(BoostNode {
            key,
            hash,
            next: NONE,
            live: true,
        });
        if self.size + 1 > self.max_load {
            self.reserve(self.size + 1);
        }
        let bucket = self.position(hash);
        self.link(bucket, node, hint);
        self.size += 1;
        node
    }

    /// Erase a node (`erase`). Nothing happens for a dead node.
    pub fn remove(&mut self, node: usize) {
        if !self.nodes[node].live {
            return;
        }
        let bucket = self.position(self.nodes[node].hash);
        let mut p = self.buckets[bucket];
        let mut prev = NONE;
        while p != NONE && p != node {
            prev = p;
            p = self.nodes[p].next;
        }
        if p == NONE {
            return;
        }
        let next = self.nodes[node].next;
        if prev == NONE {
            self.buckets[bucket] = next;
        } else {
            self.nodes[prev].next = next;
        }
        if self.buckets[bucket] == NONE {
            let group = bucket / 64;
            self.bitmask[group] &= !(1u64 << (bucket % 64));
            if self.bitmask[group] == 0 {
                self.activation.retain(|&g| g != group);
            }
        }
        self.nodes[node].live = false;
        self.size -= 1;
    }

    /// Remove every node (`clear`, `clear_impl`): the bucket array and the
    /// load threshold stay, the groups are unlinked, and the node numbers
    /// start again at 0.
    pub fn clear(&mut self) {
        self.nodes.clear();
        for bucket in &mut self.buckets {
            *bucket = NONE;
        }
        for mask in &mut self.bitmask {
            *mask = 0;
        }
        if !self.buckets.is_empty() {
            let buckets = self.buckets.len();
            self.bitmask[buckets / 64] = 1u64 << (buckets % 64);
        }
        self.activation.clear();
        self.size = 0;
    }

    /// The key of a node.
    pub fn key(&self, node: usize) -> &K {
        &self.nodes[node].key
    }

    /// The live nodes in the container's iteration order.
    pub fn order(&self) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.size);
        if self.buckets.is_empty() {
            return out;
        }
        let buckets = self.buckets.len();
        let last = buckets / 64;
        let visit_group = |group: usize, out: &mut Vec<usize>| {
            let mut mask = self.bitmask[group];
            if group == last {
                // The end sentinel's bit is not a bucket.
                mask &= !(1u64 << (buckets % 64));
            }
            while mask != 0 {
                let bit = mask.trailing_zeros() as usize;
                mask &= mask - 1;
                let mut p = self.buckets[group * 64 + bit];
                while p != NONE {
                    out.push(p);
                    p = self.nodes[p].next;
                }
            }
        };
        for &group in self.activation.iter().rev() {
            visit_group(group, &mut out);
        }
        visit_group(last, &mut out);
        out
    }
}

/// `boost::hash<std::int64_t>`: the value itself as `size_t`.
pub fn boost_hash_i64(value: i64) -> u64 {
    value as u64
}

/// `boost::hash<double>`: the bit pattern of `value + 0` (which turns `-0.0`
/// into `+0.0`), as `size_t`.
pub fn boost_hash_f64(value: f64) -> u64 {
    (value + 0.0).to_bits()
}
