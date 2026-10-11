// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The minimum-cost-flow solver of the Wasserstein network: LEMON's
//! `NetworkSimplex` with the block search pivot rule, on a `StaticDigraph`.
//!
//! Third-party dependency of the source, not an OpenMS header: the vendored
//! `pylmcf` v0.9.8 (`d2c9c52bd67d7198ae17b389d77884357260f114`) copy of LEMON
//! 1.3's `lemon/network_simplex.h` and `lemon/static_graph.h`, which
//! `FeatureGroupingAlgorithmWNet` and `WNetMatcher` reach through `wnet`
//! v0.9.11 and `wnetalign` v0.9.8. See `docs/FEATURE_GROUPING_SUPPORT.md` and
//! `docs/THIRD_PARTY_CRATE_DECISIONS.md` for why the solver is transcribed
//! rather than taken from a crate: the *optimal flow* of a min-cost-flow
//! problem is not unique when costs tie, and which optimum the simplex stops
//! at depends on its pivot order, which the consumers observe (a consensus
//! feature is the pair of a flow's heaviest arcs). A generic solver returns a
//! different optimum on exactly the inputs where it matters.
//!
//! What the transcription keeps, step by step as the executed code does:
//!
//! - `StaticDigraph` iterates nodes and arcs from the highest id down, so the
//!   solver's internal node `i` is the graph's node `N - 1 - i`, and its
//!   arcs, taken in descending id, are stored in the *mixed* order of
//!   `arc_mixing` (`skip = max(M / N, 3)`);
//! - the initial spanning tree of artificial arcs to a root node, the
//!   heuristic initial pivots (a reverse depth-first search from the single
//!   demand node over in-arcs by descending arc id, for the single supply and
//!   demand node case, and the cheapest incoming arc per demand node
//!   otherwise), and the block search pivot rule with block size
//!   `max(floor(sqrt(M)), 10)`;
//! - leaving arc selection with the strictly-feasible-tree tie rule
//!   (`d < delta` on the first path, `d <= delta` on the second), the
//!   thread/predecessor/successor tree updates and potential updates.
//!
//! Only the equality case (supplies summing to zero) is implemented: the
//! Wasserstein network always has one source with `+total` and one sink with
//! `-total`. Arithmetic wraps where the C++ would overflow a signed integer
//! (undefined behaviour there); the network's scale factor keeps every
//! reachable cost below `2^60`, so the wrap is not reached by the callers.

use crate::error::{Error, Result};

/// Arc state: at its upper bound.
const STATE_UPPER: i8 = -1;
/// Arc state: in the spanning tree.
const STATE_TREE: i8 = 0;
/// Arc state: at its lower bound.
const STATE_LOWER: i8 = 1;
/// Tree arc direction: the arc points from the node to its parent.
const DIR_UP: i8 = 1;
/// Tree arc direction: the arc points from the parent to the node.
const DIR_DOWN: i8 = -1;

/// The largest `Value`: `std::numeric_limits<int64_t>::max()`, also the
/// "infinite" capacity, since `int64_t` has no infinity.
const MAX: i64 = i64::MAX;
const INF: i64 = MAX;

/// Index of a vector element stored as the source's `int`.
#[inline]
fn ix(i: i32) -> usize {
    debug_assert!(i >= 0);
    i as usize
}

/// A minimum-cost-flow problem on a static digraph: node `n` has `supply[n]`
/// (negative for demand), arc `a` runs from `arc_source[a]` to
/// `arc_target[a]` with capacity `upper[a]` and cost `cost[a]`. Arc and node
/// numbers are the `StaticDigraph` ids.
pub(crate) struct Problem<'a> {
    pub node_num: usize,
    pub arc_source: &'a [usize],
    pub arc_target: &'a [usize],
    pub upper: &'a [i64],
    pub cost: &'a [i64],
    pub supply: &'a [i64],
}

/// The solver's mutable state (the members of `lemon::NetworkSimplex`).
struct Simplex {
    node_num: usize,
    arc_num: usize,
    all_arc_num: usize,
    search_arc_num: usize,
    sum_supply: i64,
    /// Internal arc index of each graph arc.
    arc_id: Vec<usize>,
    /// Graph arc of each internal arc index (inverse of `arc_id`), for the
    /// real arcs; used to walk graph arcs in the source's order.
    source: Vec<i32>,
    target: Vec<i32>,
    cap: Vec<i64>,
    cost: Vec<i64>,
    supply: Vec<i64>,
    flow: Vec<i64>,
    pi: Vec<i64>,
    parent: Vec<i32>,
    pred: Vec<i32>,
    thread: Vec<i32>,
    rev_thread: Vec<i32>,
    succ_num: Vec<i32>,
    last_succ: Vec<i32>,
    pred_dir: Vec<i8>,
    state: Vec<i8>,
    dirty_revs: Vec<i32>,
    root: i32,
    in_arc: i32,
    join: i32,
    u_in: i32,
    v_in: i32,
    u_out: i32,
    v_out: i32,
    delta: i64,
}

/// Block-search pivot rule state (`BlockSearchPivotRule`).
struct BlockSearch {
    block_size: i32,
    next_arc: i32,
}

/// Solve the problem with `NetworkSimplex::run(BLOCK_SEARCH)` and return the
/// flow of every arc, by `StaticDigraph` arc id.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the supplies do not sum to zero, an
/// arc names a node outside the graph, the problem is infeasible, or it is
/// unbounded (`INFEASIBLE` and `UNBOUNDED` in the source, which the callers
/// of `run()` there ignore; the Wasserstein network is feasible and bounded
/// by construction, so neither is reachable from it).
pub(crate) fn solve(problem: &Problem<'_>) -> Result<Vec<i64>> {
    let n = problem.node_num;
    let m = problem.arc_source.len();
    if problem.arc_target.len() != m
        || problem.upper.len() != m
        || problem.cost.len() != m
        || problem.supply.len() != n
    {
        return Err(Error::InvalidValue(
            "network simplex: the arc and node arrays do not match".into(),
        ));
    }
    if n == 0 {
        return Err(Error::InvalidValue(
            "network simplex: the graph has no nodes".into(),
        ));
    }
    if i32::try_from(n + 1 + 2 * n + m).is_err() {
        return Err(Error::InvalidValue(
            "network simplex: the graph is too large for 32-bit indices".into(),
        ));
    }
    for a in 0..m {
        if problem.arc_source[a] >= n || problem.arc_target[a] >= n {
            return Err(Error::InvalidValue(
                "network simplex: an arc names a node outside the graph".into(),
            ));
        }
    }
    let mut sum: i64 = 0;
    for &s in problem.supply {
        sum = sum.wrapping_add(s);
    }
    if sum != 0 {
        return Err(Error::InvalidValue(
            "network simplex: only supplies that sum to zero are supported".into(),
        ));
    }

    let mut simplex = Simplex::new(problem);
    simplex.init();
    simplex.start()?;
    Ok((0..m).map(|a| simplex.flow[simplex.arc_id[a]]).collect())
}

impl Simplex {
    /// `reset()` followed by the `upperMap`, `costMap` and `supplyMap` calls.
    fn new(problem: &Problem<'_>) -> Self {
        let node_num = problem.node_num;
        let arc_num = problem.arc_source.len();
        let all_node_num = node_num + 1;
        let max_arc_num = arc_num + 2 * node_num;

        // StaticDigraph::first(Node) is the highest id and next decrements, so
        // the n-th visited node (internal id n) is graph node N - 1 - n.
        let node_id = |graph_node: usize| node_num - 1 - graph_node;

        let mut arc_id = vec![0usize; arc_num];
        let mut source = vec![0i32; max_arc_num];
        let mut target = vec![0i32; max_arc_num];
        let mut cap = vec![0i64; max_arc_num];
        let mut cost = vec![0i64; max_arc_num];
        let mut supply = vec![0i64; all_node_num];

        // Arcs are visited from the highest id down.
        if node_num > 1 {
            let skip = std::cmp::max(arc_num / node_num, 3);
            let mut i = 0usize;
            let mut j = 0usize;
            for a in (0..arc_num).rev() {
                arc_id[a] = i;
                source[i] = node_id(problem.arc_source[a]) as i32;
                target[i] = node_id(problem.arc_target[a]) as i32;
                cap[i] = problem.upper[a];
                cost[i] = problem.cost[a];
                i += skip;
                if i >= arc_num {
                    j += 1;
                    i = j;
                }
            }
        } else {
            for (i, a) in (0..arc_num).rev().enumerate() {
                arc_id[a] = i;
                source[i] = node_id(problem.arc_source[a]) as i32;
                target[i] = node_id(problem.arc_target[a]) as i32;
                cap[i] = problem.upper[a];
                cost[i] = problem.cost[a];
            }
        }
        for (n, &s) in problem.supply.iter().enumerate() {
            supply[node_id(n)] = s;
        }

        Self {
            node_num,
            arc_num,
            all_arc_num: 0,
            search_arc_num: 0,
            sum_supply: 0,
            arc_id,
            source,
            target,
            cap,
            cost,
            supply,
            flow: vec![0i64; max_arc_num],
            pi: vec![0i64; all_node_num],
            parent: vec![0i32; all_node_num],
            pred: vec![0i32; all_node_num],
            thread: vec![0i32; all_node_num],
            rev_thread: vec![0i32; all_node_num],
            succ_num: vec![0i32; all_node_num],
            last_succ: vec![0i32; all_node_num],
            pred_dir: vec![0i8; all_node_num],
            state: vec![0i8; max_arc_num],
            dirty_revs: Vec::new(),
            root: 0,
            in_arc: 0,
            join: 0,
            u_in: 0,
            v_in: 0,
            u_out: 0,
            v_out: 0,
            delta: 0,
        }
    }

    /// `init()` for the equality case: one artificial arc per node to or from
    /// the root, all in the initial spanning tree.
    fn init(&mut self) {
        let node_num = self.node_num;
        let arc_num = self.arc_num;
        // The caller checked that the supplies sum to zero.
        self.sum_supply = 0;
        // ART_COST = max / 2 + 1 for an exact cost type.
        let art_cost: i64 = i64::MAX / 2 + 1;

        for i in 0..arc_num {
            self.flow[i] = 0;
            self.state[i] = STATE_LOWER;
        }

        let root = node_num as i32;
        self.root = root;
        let r = ix(root);
        self.parent[r] = -1;
        self.pred[r] = -1;
        self.thread[r] = 0;
        self.rev_thread[0] = root;
        self.succ_num[r] = node_num as i32 + 1;
        self.last_succ[r] = root - 1;
        self.supply[r] = -self.sum_supply;
        self.pi[r] = 0;

        self.search_arc_num = arc_num;
        self.all_arc_num = arc_num + node_num;
        for u in 0..node_num {
            let e = arc_num + u;
            self.parent[u] = root;
            self.pred[u] = e as i32;
            self.thread[u] = u as i32 + 1;
            self.rev_thread[u + 1] = u as i32;
            self.succ_num[u] = 1;
            self.last_succ[u] = u as i32;
            self.cap[e] = INF;
            self.state[e] = STATE_TREE;
            if self.supply[u] >= 0 {
                self.pred_dir[u] = DIR_UP;
                self.pi[u] = 0;
                self.source[e] = u as i32;
                self.target[e] = root;
                self.flow[e] = self.supply[u];
                self.cost[e] = 0;
            } else {
                self.pred_dir[u] = DIR_DOWN;
                self.pi[u] = art_cost;
                self.source[e] = root;
                self.target[e] = u as i32;
                self.flow[e] = self.supply[u].wrapping_neg();
                self.cost[e] = art_cost;
            }
        }
    }

    /// The reduced cost of arc `e`, signed by its state.
    #[inline]
    fn signed_reduced_cost(&self, e: usize) -> i64 {
        i64::from(self.state[e]).wrapping_mul(
            self.cost[e]
                .wrapping_add(self.pi[ix(self.source[e])])
                .wrapping_sub(self.pi[ix(self.target[e])]),
        )
    }

    /// `findJoinNode()`.
    fn find_join_node(&mut self) {
        let mut u = self.source[ix(self.in_arc)];
        let mut v = self.target[ix(self.in_arc)];
        while u != v {
            if self.succ_num[ix(u)] < self.succ_num[ix(v)] {
                u = self.parent[ix(u)];
            } else {
                v = self.parent[ix(v)];
            }
        }
        self.join = u;
    }

    /// `findLeavingArc()`: true when the leaving arc is not the entering arc.
    fn find_leaving_arc(&mut self) -> bool {
        let in_arc = ix(self.in_arc);
        let (first, second) = if self.state[in_arc] == STATE_LOWER {
            (self.source[in_arc], self.target[in_arc])
        } else {
            (self.target[in_arc], self.source[in_arc])
        };
        self.delta = self.cap[in_arc];
        let mut result = 0;

        // The cycle from the first node to the join node.
        let mut u = first;
        while u != self.join {
            let e = ix(self.pred[ix(u)]);
            let mut d = self.flow[e];
            if self.pred_dir[ix(u)] == DIR_DOWN {
                let c = self.cap[e];
                d = if c == MAX { INF } else { c.wrapping_sub(d) };
            }
            if d < self.delta {
                self.delta = d;
                self.u_out = u;
                result = 1;
            }
            u = self.parent[ix(u)];
        }
        // The cycle from the second node to the join node.
        let mut u = second;
        while u != self.join {
            let e = ix(self.pred[ix(u)]);
            let mut d = self.flow[e];
            if self.pred_dir[ix(u)] == DIR_UP {
                let c = self.cap[e];
                d = if c == MAX { INF } else { c.wrapping_sub(d) };
            }
            if d <= self.delta {
                self.delta = d;
                self.u_out = u;
                result = 2;
            }
            u = self.parent[ix(u)];
        }

        if result == 1 {
            self.u_in = first;
            self.v_in = second;
        } else {
            self.u_in = second;
            self.v_in = first;
        }
        result != 0
    }

    /// `changeFlow(change)`.
    fn change_flow(&mut self, change: bool) {
        let in_arc = ix(self.in_arc);
        if self.delta > 0 {
            let val = i64::from(self.state[in_arc]).wrapping_mul(self.delta);
            self.flow[in_arc] = self.flow[in_arc].wrapping_add(val);
            let mut u = self.source[in_arc];
            while u != self.join {
                let p = ix(self.pred[ix(u)]);
                self.flow[p] =
                    self.flow[p].wrapping_sub(i64::from(self.pred_dir[ix(u)]).wrapping_mul(val));
                u = self.parent[ix(u)];
            }
            let mut u = self.target[in_arc];
            while u != self.join {
                let p = ix(self.pred[ix(u)]);
                self.flow[p] =
                    self.flow[p].wrapping_add(i64::from(self.pred_dir[ix(u)]).wrapping_mul(val));
                u = self.parent[ix(u)];
            }
        }
        if change {
            self.state[in_arc] = STATE_TREE;
            let leaving = ix(self.pred[ix(self.u_out)]);
            self.state[leaving] = if self.flow[leaving] == 0 {
                STATE_LOWER
            } else {
                STATE_UPPER
            };
        } else {
            self.state[in_arc] = -self.state[in_arc];
        }
    }

    /// `updateTreeStructure()`.
    fn update_tree_structure(&mut self) {
        let u_out = self.u_out;
        let u_in = self.u_in;
        let v_in = self.v_in;
        let join = self.join;
        let in_arc = ix(self.in_arc);
        let old_rev_thread = self.rev_thread[ix(u_out)];
        let old_succ_num = self.succ_num[ix(u_out)];
        let old_last_succ = self.last_succ[ix(u_out)];
        let v_out = self.parent[ix(u_out)];
        self.v_out = v_out;

        if u_in == u_out {
            // Update parent, pred, pred_dir.
            self.parent[ix(u_in)] = v_in;
            self.pred[ix(u_in)] = self.in_arc;
            self.pred_dir[ix(u_in)] = if u_in == self.source[in_arc] {
                DIR_UP
            } else {
                DIR_DOWN
            };

            // Update thread and rev_thread.
            if self.thread[ix(v_in)] != u_out {
                let mut after = self.thread[ix(old_last_succ)];
                self.thread[ix(old_rev_thread)] = after;
                self.rev_thread[ix(after)] = old_rev_thread;
                after = self.thread[ix(v_in)];
                self.thread[ix(v_in)] = u_out;
                self.rev_thread[ix(u_out)] = v_in;
                self.thread[ix(old_last_succ)] = after;
                self.rev_thread[ix(after)] = old_last_succ;
            }
        } else {
            // old_rev_thread equals v_in also means join and v_out coincide.
            let thread_continue = if old_rev_thread == v_in {
                self.thread[ix(old_last_succ)]
            } else {
                self.thread[ix(v_in)]
            };

            // Update thread and parent along the stem nodes.
            let mut stem = u_in;
            let mut par_stem = v_in;
            let mut last = self.last_succ[ix(u_in)];
            let mut after = self.thread[ix(last)];
            self.thread[ix(v_in)] = u_in;
            self.dirty_revs.clear();
            self.dirty_revs.push(v_in);
            while stem != u_out {
                // Insert the next stem node into the thread list.
                let next_stem = self.parent[ix(stem)];
                self.thread[ix(last)] = next_stem;
                self.dirty_revs.push(last);

                // Remove the subtree of stem from the thread list.
                let before = self.rev_thread[ix(stem)];
                self.thread[ix(before)] = after;
                self.rev_thread[ix(after)] = before;

                // Change the parent node and shift stem nodes.
                self.parent[ix(stem)] = par_stem;
                par_stem = stem;
                stem = next_stem;

                // Update last and after.
                last = if self.last_succ[ix(stem)] == self.last_succ[ix(par_stem)] {
                    self.rev_thread[ix(par_stem)]
                } else {
                    self.last_succ[ix(stem)]
                };
                after = self.thread[ix(last)];
            }
            self.parent[ix(u_out)] = par_stem;
            self.thread[ix(last)] = thread_continue;
            self.rev_thread[ix(thread_continue)] = last;
            self.last_succ[ix(u_out)] = last;

            // Remove the subtree of u_out from the thread list except for the
            // case when old_rev_thread equals v_in.
            if old_rev_thread != v_in {
                self.thread[ix(old_rev_thread)] = after;
                self.rev_thread[ix(after)] = old_rev_thread;
            }

            // Update rev_thread using the new thread values.
            for i in 0..self.dirty_revs.len() {
                let u = self.dirty_revs[i];
                let t = self.thread[ix(u)];
                self.rev_thread[ix(t)] = u;
            }

            // Update pred, pred_dir, last_succ and succ_num for the stem nodes
            // from u_out to u_in.
            let mut tmp_sc: i32 = 0;
            let tmp_ls = self.last_succ[ix(u_out)];
            let mut u = u_out;
            let mut p = self.parent[ix(u)];
            while u != u_in {
                self.pred[ix(u)] = self.pred[ix(p)];
                self.pred_dir[ix(u)] = -self.pred_dir[ix(p)];
                tmp_sc += self.succ_num[ix(u)] - self.succ_num[ix(p)];
                self.succ_num[ix(u)] = tmp_sc;
                self.last_succ[ix(p)] = tmp_ls;
                u = p;
                p = self.parent[ix(u)];
            }
            self.pred[ix(u_in)] = self.in_arc;
            self.pred_dir[ix(u_in)] = if u_in == self.source[in_arc] {
                DIR_UP
            } else {
                DIR_DOWN
            };
            self.succ_num[ix(u_in)] = old_succ_num;
        }

        // Update last_succ from v_in towards the root.
        let up_limit_out = if self.last_succ[ix(join)] == v_in {
            join
        } else {
            -1
        };
        let last_succ_out = self.last_succ[ix(u_out)];
        let mut u = v_in;
        while u != -1 && self.last_succ[ix(u)] == v_in {
            self.last_succ[ix(u)] = last_succ_out;
            u = self.parent[ix(u)];
        }

        // Update last_succ from v_out towards the root.
        if join != old_rev_thread && v_in != old_rev_thread {
            let mut u = v_out;
            while u != up_limit_out && self.last_succ[ix(u)] == old_last_succ {
                self.last_succ[ix(u)] = old_rev_thread;
                u = self.parent[ix(u)];
            }
        } else if last_succ_out != old_last_succ {
            let mut u = v_out;
            while u != up_limit_out && self.last_succ[ix(u)] == old_last_succ {
                self.last_succ[ix(u)] = last_succ_out;
                u = self.parent[ix(u)];
            }
        }

        // Update succ_num from v_in to join.
        let mut u = v_in;
        while u != join {
            self.succ_num[ix(u)] += old_succ_num;
            u = self.parent[ix(u)];
        }
        // Update succ_num from v_out to join.
        let mut u = v_out;
        while u != join {
            self.succ_num[ix(u)] -= old_succ_num;
            u = self.parent[ix(u)];
        }
    }

    /// `updatePotential()`: the potentials of the subtree that has been moved.
    fn update_potential(&mut self) {
        let in_arc = ix(self.in_arc);
        let u_in = ix(self.u_in);
        let sigma = self.pi[ix(self.v_in)]
            .wrapping_sub(self.pi[u_in])
            .wrapping_sub(i64::from(self.pred_dir[u_in]).wrapping_mul(self.cost[in_arc]));
        let end = self.thread[ix(self.last_succ[u_in])];
        let mut u = self.u_in;
        while u != end {
            self.pi[ix(u)] = self.pi[ix(u)].wrapping_add(sigma);
            u = self.thread[ix(u)];
        }
    }

    /// `initialPivots()`: false when a pivot is unbounded.
    fn initial_pivots(&mut self, problem_in_arcs: &[Vec<usize>]) -> bool {
        // Nodes are visited from the highest graph id down; internal ids are
        // reversed, so the supply and demand nodes are visited in descending
        // graph id, i.e. ascending internal id.
        let mut total: i64 = 0;
        let mut supply_nodes: Vec<usize> = Vec::new();
        let mut demand_nodes: Vec<usize> = Vec::new();
        for internal in 0..self.node_num {
            let curr = self.supply[internal];
            if curr > 0 {
                total = total.wrapping_add(curr);
                supply_nodes.push(internal);
            } else if curr < 0 {
                demand_nodes.push(internal);
            }
        }
        if self.sum_supply > 0 {
            total = total.wrapping_sub(self.sum_supply);
        }
        if total <= 0 {
            return true;
        }

        let mut arc_vector: Vec<usize> = Vec::new();
        if self.sum_supply >= 0 {
            if supply_nodes.len() == 1 && demand_nodes.len() == 1 {
                // A reverse graph search from the sink to the source over
                // arcs of capacity at least the total supply.
                let mut reached = vec![false; self.node_num];
                let s = supply_nodes[0];
                let t = demand_nodes[0];
                let mut stack: Vec<usize> = Vec::new();
                reached[t] = true;
                stack.push(t);
                while let Some(v) = stack.pop() {
                    if v == s {
                        break;
                    }
                    for &j in &problem_in_arcs[v] {
                        let u = ix(self.source[j]);
                        if reached[u] {
                            continue;
                        }
                        if self.cap[j] >= total {
                            arc_vector.push(j);
                            reached[u] = true;
                            stack.push(u);
                        }
                    }
                }
            } else {
                // The cheapest incoming arc of each demand node.
                for &v in &demand_nodes {
                    let mut min_cost = i64::MAX;
                    let mut min_arc: Option<usize> = None;
                    for &j in &problem_in_arcs[v] {
                        let c = self.cost[j];
                        if c < min_cost {
                            min_cost = c;
                            min_arc = Some(j);
                        }
                    }
                    if let Some(j) = min_arc {
                        arc_vector.push(j);
                    }
                }
            }
        }

        for &j in &arc_vector {
            self.in_arc = j as i32;
            if self.signed_reduced_cost(j) >= 0 {
                continue;
            }
            self.find_join_node();
            let change = self.find_leaving_arc();
            if self.delta == MAX {
                return false;
            }
            self.change_flow(change);
            if change {
                self.update_tree_structure();
                self.update_potential();
            }
        }
        true
    }

    /// `findEnteringArc()` of the block search rule.
    fn find_entering_arc(&mut self, rule: &mut BlockSearch) -> bool {
        let mut min: i64 = 0;
        let mut cnt = rule.block_size;
        let search_arc_num = self.search_arc_num as i32;
        let mut e = rule.next_arc;
        let mut found = false;
        while e != search_arc_num {
            let c = self.signed_reduced_cost(ix(e));
            if c < min {
                min = c;
                self.in_arc = e;
            }
            cnt -= 1;
            if cnt == 0 {
                if min < 0 {
                    found = true;
                    break;
                }
                cnt = rule.block_size;
            }
            e += 1;
        }
        if !found {
            e = 0;
            while e != rule.next_arc {
                let c = self.signed_reduced_cost(ix(e));
                if c < min {
                    min = c;
                    self.in_arc = e;
                }
                cnt -= 1;
                if cnt == 0 {
                    if min < 0 {
                        found = true;
                        break;
                    }
                    cnt = rule.block_size;
                }
                e += 1;
            }
        }
        if !found && min >= 0 {
            return false;
        }
        rule.next_arc = e;
        true
    }

    /// `start<BlockSearchPivotRule>()`.
    fn start(&mut self) -> Result<()> {
        let block_size = std::cmp::max((self.search_arc_num as f64).sqrt() as i32, 10);
        let mut rule = BlockSearch {
            block_size,
            next_arc: 0,
        };

        // In-arcs of every internal node by descending graph arc id, which is
        // how `StaticDigraph::InArcIt` lists them. The arcs are stored in the
        // mixed internal order, so rebuild the list from the graph order.
        let mut graph_order: Vec<usize> = vec![0; self.arc_num];
        // arc_id maps graph arc -> internal; invert it.
        for (graph_arc, &internal) in self.arc_id.iter().enumerate() {
            graph_order[internal] = graph_arc;
        }
        let mut in_arcs: Vec<Vec<usize>> = vec![Vec::new(); self.node_num];
        let mut by_graph: Vec<usize> = (0..self.arc_num).collect();
        by_graph.sort_unstable_by_key(|&internal| std::cmp::Reverse(graph_order[internal]));
        for internal in by_graph {
            in_arcs[ix(self.target[internal])].push(internal);
        }

        if !self.initial_pivots(&in_arcs) {
            return Err(Error::InvalidValue(
                "network simplex: the problem is unbounded".into(),
            ));
        }

        while self.find_entering_arc(&mut rule) {
            self.find_join_node();
            let change = self.find_leaving_arc();
            if self.delta == MAX {
                return Err(Error::InvalidValue(
                    "network simplex: the problem is unbounded".into(),
                ));
            }
            self.change_flow(change);
            if change {
                self.update_tree_structure();
                self.update_potential();
            }
        }

        for e in self.search_arc_num..self.all_arc_num {
            if self.flow[e] != 0 {
                return Err(Error::InvalidValue(
                    "network simplex: the problem is infeasible".into(),
                ));
            }
        }
        Ok(())
    }
}
