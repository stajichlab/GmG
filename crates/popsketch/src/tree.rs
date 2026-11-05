//! Neighbour-joining tree on mash distances, with branches below the noise floor collapsed to polytomies, and Newick output.
//! The tree is unrooted unless a tip is named as the root (the root is then put on that tip's branch).

#[derive(Clone, Debug)]
pub struct Tree {
    pub names: Vec<String>,
    /// nodes 0..names.len() are tips (degree 1); the rest are internal
    pub adj: Vec<Vec<(usize, f64)>>,
}

impl Tree {
    pub fn n_tips(&self) -> usize {
        self.names.len()
    }
}

fn add_edge(adj: &mut [Vec<(usize, f64)>], a: usize, b: usize, len: f64) {
    adj[a].push((b, len));
    adj[b].push((a, len));
}

/// Neighbour-joining on a symmetric n x n distance matrix (row-major). Negative branch lengths are set to zero.
pub fn neighbor_joining(d: &[f64], n: usize, names: Vec<String>) -> Tree {
    assert_eq!(names.len(), n);
    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    if n == 2 {
        add_edge(&mut adj, 0, 1, d[1].max(0.0));
        return Tree { names, adj };
    }
    if n < 2 {
        return Tree { names, adj };
    }
    let mut dist = d.to_vec();
    let mut node: Vec<usize> = (0..n).collect(); // tree node held by each slot
    let mut active: Vec<usize> = (0..n).collect();
    while active.len() > 3 {
        let m = active.len();
        let r: Vec<f64> = active
            .iter()
            .map(|&i| active.iter().map(|&k| dist[i * n + k]).sum())
            .collect();
        let (mut bi, mut bj, mut best) = (0, 1, f64::INFINITY);
        for x in 0..m {
            for y in (x + 1)..m {
                let q = (m as f64 - 2.0) * dist[active[x] * n + active[y]] - r[x] - r[y];
                if q < best {
                    (bi, bj, best) = (x, y, q);
                }
            }
        }
        let (i, j) = (active[bi], active[bj]);
        let dij = dist[i * n + j];
        let li = 0.5 * dij + (r[bi] - r[bj]) / (2.0 * (m as f64 - 2.0));
        let lj = dij - li;
        let u = adj.len();
        adj.push(Vec::new());
        add_edge(&mut adj, u, node[i], li.max(0.0));
        add_edge(&mut adj, u, node[j], lj.max(0.0));
        for &k in &active {
            if k != i && k != j {
                let v = 0.5 * (dist[i * n + k] + dist[j * n + k] - dij);
                dist[i * n + k] = v;
                dist[k * n + i] = v;
            }
        }
        node[i] = u;
        active.remove(bj);
    }
    let (a, b, c) = (active[0], active[1], active[2]);
    let (dab, dac, dbc) = (dist[a * n + b], dist[a * n + c], dist[b * n + c]);
    let u = adj.len();
    adj.push(Vec::new());
    add_edge(&mut adj, u, node[a], (0.5 * (dab + dac - dbc)).max(0.0));
    add_edge(&mut adj, u, node[b], (0.5 * (dab + dbc - dac)).max(0.0));
    add_edge(&mut adj, u, node[c], (0.5 * (dac + dbc - dab)).max(0.0));
    Tree { names, adj }
}

/// Contract every internal-internal branch shorter than `floor` (a polytomy replaces it). Returns the number collapsed.
pub fn collapse(t: &mut Tree, floor: f64) -> usize {
    let n = t.n_tips();
    let mut collapsed = 0;
    loop {
        let mut found = None;
        'scan: for u in n..t.adj.len() {
            for &(v, len) in &t.adj[u] {
                if v >= n && v > u && len < floor {
                    found = Some((u, v));
                    break 'scan;
                }
            }
        }
        let Some((u, v)) = found else { break };
        t.adj[u].retain(|&(w, _)| w != v);
        let moved: Vec<(usize, f64)> = std::mem::take(&mut t.adj[v])
            .into_iter()
            .filter(|&(w, _)| w != u)
            .collect();
        for (w, len) in moved {
            for e in t.adj[w].iter_mut() {
                if e.0 == v {
                    e.0 = u;
                }
            }
            t.adj[u].push((w, len));
        }
        collapsed += 1;
    }
    collapsed
}

#[derive(Clone, Debug)]
pub struct RNode {
    /// tip index (into `Tree::names`) for a tip
    pub tip: Option<usize>,
    /// length of the branch above this node
    pub len: f64,
    pub children: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct Rooted {
    pub nodes: Vec<RNode>,
    pub root: usize,
}

/// Root the tree: on the branch above tip `at` (halving it), or, without `at`, at the internal node of highest degree.
/// Internal nodes left without a neighbour after collapsing are never reached.
pub fn root_tree(t: &Tree, at: Option<&str>) -> Rooted {
    let n = t.n_tips();
    let mut nodes: Vec<RNode> = Vec::new();
    let push =
        |nodes: &mut Vec<RNode>, tip: Option<usize>, len: f64, parent: Option<usize>| -> usize {
            nodes.push(RNode {
                tip,
                len,
                children: Vec::new(),
            });
            let id = nodes.len() - 1;
            if let Some(p) = parent {
                nodes[p].children.push(id);
            }
            id
        };
    // (tree node, the node it was reached from, branch length, parent in the rooted tree)
    let mut work: Vec<(usize, usize, f64, usize)> = Vec::new();
    let tip_root = at
        .and_then(|name| t.names.iter().position(|x| x == name))
        .filter(|&tip| t.adj[tip].len() == 1);
    let root = match tip_root {
        Some(tip) => {
            let (p, len) = t.adj[tip][0];
            let root = push(&mut nodes, None, 0.0, None);
            push(&mut nodes, Some(tip), len / 2.0, Some(root));
            work.push((p, tip, len / 2.0, root));
            root
        }
        None => {
            let r = (n..t.adj.len())
                .filter(|&u| !t.adj[u].is_empty())
                .max_by_key(|&u| (t.adj[u].len(), std::cmp::Reverse(u)))
                .unwrap_or(0);
            let root = push(&mut nodes, (r < n).then_some(r), 0.0, None);
            for &(v, len) in t.adj[r].iter().rev() {
                work.push((v, r, len, root));
            }
            root
        }
    };
    while let Some((v, from, len, par)) = work.pop() {
        let id = push(&mut nodes, (v < n).then_some(v), len, Some(par));
        if v >= n {
            for &(w, l) in t.adj[v].iter().rev() {
                if w != from {
                    work.push((w, v, l, id));
                }
            }
        }
    }
    Rooted { nodes, root }
}

fn clean(name: &str) -> String {
    name.chars()
        .map(|c| {
            if matches!(
                c,
                '(' | ')' | ',' | ':' | ';' | '[' | ']' | '\'' | ' ' | '\t'
            ) {
                '_'
            } else {
                c
            }
        })
        .collect()
}

/// Newick text of a rooted tree. Characters that Newick reserves are replaced by underscores in tip names.
pub fn newick(r: &Rooted, names: &[String]) -> String {
    fn rec(r: &Rooted, names: &[String], i: usize, out: &mut String, top: bool) {
        let node = &r.nodes[i];
        if !node.children.is_empty() {
            out.push('(');
            for (k, &c) in node.children.iter().enumerate() {
                if k > 0 {
                    out.push(',');
                }
                rec(r, names, c, out, false);
            }
            out.push(')');
        }
        if let Some(t) = node.tip {
            out.push_str(&clean(&names[t]));
        }
        if !top {
            out.push_str(&format!(":{:.6}", node.len));
        }
    }
    let mut s = String::new();
    rec(r, names, r.root, &mut s, true);
    s.push(';');
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// path length between two tree nodes
    fn path(t: &Tree, a: usize, b: usize) -> f64 {
        let mut stack = vec![(a, usize::MAX, 0.0)];
        while let Some((v, from, d)) = stack.pop() {
            if v == b {
                return d;
            }
            for &(w, l) in &t.adj[v] {
                if w != from {
                    stack.push((w, v, d + l));
                }
            }
        }
        f64::NAN
    }

    /// ((A:1,B:2):3,(C:4,D:5):6,E:7), as pairwise path lengths
    fn additive() -> (Vec<f64>, Vec<String>) {
        let names: Vec<String> = "ABCDE".chars().map(|c| c.to_string()).collect();
        let leg = [1.0, 2.0, 4.0, 5.0, 7.0];
        let (ab, cd) = (3.0, 6.0);
        let mut d = vec![0.0; 25];
        for i in 0..5 {
            for j in 0..5 {
                if i == j {
                    continue;
                }
                let v = match (i.min(j), i.max(j)) {
                    (0, 1) | (2, 3) => leg[i] + leg[j],
                    (0, 2) | (0, 3) | (1, 2) | (1, 3) => leg[i] + leg[j] + ab + cd,
                    (a, 4) => leg[i] + leg[j] + if a < 2 { ab } else { cd },
                    _ => unreachable!(),
                };
                d[i * 5 + j] = v;
            }
        }
        (d, names)
    }

    #[test]
    fn nj_recovers_an_additive_tree_exactly() {
        let (d, names) = additive();
        let t = neighbor_joining(&d, 5, names);
        for i in 0..5 {
            for j in 0..5 {
                assert!((path(&t, i, j) - d[i * 5 + j]).abs() < 1e-9, "({i},{j})");
            }
            assert_eq!(t.adj[i].len(), 1);
        }
        // leg lengths of the tips are recovered
        assert!((t.adj[0][0].1 - 1.0).abs() < 1e-9 && (t.adj[4][0].1 - 7.0).abs() < 1e-9);
    }

    #[test]
    fn small_trees() {
        let t2 = neighbor_joining(&[0.0, 0.3, 0.3, 0.0], 2, vec!["a".into(), "b".into()]);
        assert!((path(&t2, 0, 1) - 0.3).abs() < 1e-12);
        let d3 = [0.0, 2.0, 3.0, 2.0, 0.0, 3.0, 3.0, 3.0, 0.0];
        let t3 = neighbor_joining(&d3, 3, vec!["a".into(), "b".into(), "c".into()]);
        assert!((path(&t3, 0, 1) - 2.0).abs() < 1e-9 && (path(&t3, 0, 2) - 3.0).abs() < 1e-9);
    }

    #[test]
    fn short_internal_branches_become_a_polytomy() {
        // ((A:1,B:1):0.0001,(C:1,D:1):5,E:3): the 0.0001 branch is below a floor of 0.001
        let names: Vec<String> = "ABCDE".chars().map(|c| c.to_string()).collect();
        let mut d = vec![0.0; 25];
        let leg = [1.0, 1.0, 1.0, 1.0, 3.0];
        for i in 0..5 {
            for j in 0..5 {
                if i != j {
                    let (a, b) = (i.min(j), i.max(j));
                    let extra = match (a, b) {
                        (0, 1) | (2, 3) => 0.0,
                        (0, 4) | (1, 4) => 0.0001,
                        (2, 4) | (3, 4) => 5.0,
                        _ => 5.0001,
                    };
                    d[i * 5 + j] = leg[i] + leg[j] + extra;
                }
            }
        }
        let mut t = neighbor_joining(&d, 5, names.clone());
        let internal_before = (5..t.adj.len()).filter(|&u| !t.adj[u].is_empty()).count();
        let c = collapse(&mut t, 0.001);
        let internal_after = (5..t.adj.len()).filter(|&u| !t.adj[u].is_empty()).count();
        assert_eq!(internal_before, 3);
        assert_eq!(c, 1);
        assert_eq!(internal_after, 2);
        // distances between tips are preserved except for the removed sliver
        assert!((path(&t, 0, 2) - d[2]).abs() < 0.001);
    }

    #[test]
    fn newick_is_well_formed_and_roots_on_a_tip() {
        let (d, names) = additive();
        let t = neighbor_joining(&d, 5, names.clone());
        let r = root_tree(&t, None);
        let s = newick(&r, &names);
        assert!(s.ends_with(");") && s.starts_with('('));
        assert_eq!(s.matches('(').count(), s.matches(')').count());
        for nm in &names {
            assert_eq!(s.matches(nm.as_str()).count(), 1, "{s}");
        }
        // rooted on E: the root has two children, E with half of its 7-unit branch
        let re = root_tree(&t, Some("E"));
        assert_eq!(re.nodes[re.root].children.len(), 2);
        let e = re.nodes.iter().find(|x| x.tip == Some(4)).unwrap();
        assert!((e.len - 3.5).abs() < 1e-9);
        let se = newick(&re, &names);
        assert_eq!(se.matches('(').count(), se.matches(')').count());
        // the same total tree length either way
        let total = |r: &Rooted| r.nodes.iter().map(|x| x.len).sum::<f64>();
        assert!((total(&r) - total(&re)).abs() < 1e-9);
    }

    #[test]
    fn reserved_characters_are_cleaned_from_names() {
        assert_eq!(clean("a b(c):d"), "a_b_c__d");
    }
}
