//! Elementary collapse and strong collapse of simplicial complexes.
//!
//! An **elementary collapse** removes a free face σ and its unique coface τ
//! (where dim(τ) = dim(σ) + 1). This preserves the homotopy type.
//!
//! A simplex σ is a **free face** if it is contained in exactly one simplex
//! of dimension dim(σ) + 1.
//!
//! A complex is **collapsible** if it can be reduced to a single vertex by
//! a sequence of elementary collapses. Every cone is collapsible.

use crate::complex::SimplicialComplex;
use crate::simplex::Simplex;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Result of a collapse operation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CollapseResult {
    /// The collapsed complex.
    pub complex: SimplicialComplex,
    /// Number of collapse steps performed.
    pub steps: usize,
    /// Whether the complex is collapsible (reduced to a single vertex).
    pub is_collapsible: bool,
    /// Sequence of removed (free face, coface) pairs.
    pub removed_pairs: Vec<(Vec<usize>, Vec<usize>)>,
}

/// Find a free face in the complex.
///
/// A free face is a simplex σ that has exactly one coface of dimension dim(σ)+1.
/// Returns (free_face, coface) if found.
fn find_free_face(complex: &SimplicialComplex) -> Option<(Simplex, Simplex)> {
    let max_dim = complex.max_dimension();
    if max_dim == 0 {
        return None; // Can't collapse points further (unless single vertex)
    }

    for dim in 0..max_dim {
        let simplices = complex.simplices_of_dimension(dim);
        for sigma in simplices {
            // Count cofaces of dimension dim+1
            let cofaces = complex.simplices_of_dimension(dim + 1);
            let cofaces: Vec<&Simplex> = cofaces
                .iter()
                .filter(|tau| sigma.is_face_of(tau))
                .copied()
                .collect();

            if cofaces.len() == 1 {
                return Some((sigma.clone(), cofaces.into_iter().next().unwrap().clone()));
            }
        }
    }
    None
}

/// Perform elementary collapse iteratively until no more free faces exist.
///
/// Uses a queue-based approach: find a free face, remove it and its coface,
/// repeat until no more free faces are found.
///
/// This preserves the homotopy type of the complex.
pub fn collapse(complex: &SimplicialComplex) -> CollapseResult {
    let mut current = complex.clone();
    let mut steps = 0;
    let mut removed_pairs = Vec::new();

    while let Some((sigma, tau)) = find_free_face(&current) {
        removed_pairs.push((sigma.vertices().to_vec(), tau.vertices().to_vec()));
        let sigma_key = sigma.vertices().to_vec();
        let tau_key = tau.vertices().to_vec();
        current.simplices.remove(&sigma_key);
        current.simplices.remove(&tau_key);
        current.recompute_f_vector();
        steps += 1;
    }

    let is_collapsible = current.num_simplices() == 1 && current.num_vertices() == 1;

    CollapseResult {
        complex: current,
        steps,
        is_collapsible,
        removed_pairs,
    }
}

/// Strong collapse: iteratively remove dominated vertices.
///
/// A vertex *v* is **dominated** by vertex *w* if the closed neighborhood
/// N[v] ⊆ N[w] in the 1-skeleton graph, where N[v] = {v} ∪ {u : {u,v} ∈ K}.
///
/// Strong collapse preserves the homotopy type and is often more powerful
/// than elementary collapse. Based on Barmak & Minian (2008).
pub fn strong_collapse(complex: &SimplicialComplex) -> CollapseResult {
    let mut current = complex.clone();
    let mut steps = 0;
    let mut removed_pairs = Vec::new();

    loop {
        let vertices: Vec<usize> = current
            .simplices_of_dimension(0)
            .iter()
            .map(|s| s.vertices()[0])
            .collect();

        if vertices.len() <= 1 {
            break;
        }

        // Build adjacency for 1-skeleton
        let mut adjacency: HashMap<usize, HashSet<usize>> = HashMap::new();
        for &v in &vertices {
            let mut neighbors = HashSet::new();
            neighbors.insert(v); // closed neighborhood includes self
            adjacency.insert(v, neighbors);
        }
        for edge in current.simplices_of_dimension(1) {
            let vs = edge.vertices();
            if vs.len() == 2 {
                adjacency.get_mut(&vs[0]).unwrap().insert(vs[1]);
                adjacency.get_mut(&vs[1]).unwrap().insert(vs[0]);
            }
        }

        let mut dominated: Option<usize> = None;

        'outer: for &v in &vertices {
            for &w in &vertices {
                if w == v {
                    continue;
                }
                // Check N[v] ⊆ N[w]
                if adjacency[&v].iter().all(|u| adjacency[&w].contains(u)) {
                    dominated = Some(v);
                    break 'outer;
                }
            }
        }

        match dominated {
            Some(v) => {
                // Remove all simplices containing v
                let to_remove: Vec<Vec<usize>> = current
                    .simplices()
                    .filter(|s| s.vertices().contains(&v))
                    .map(|s| s.vertices().to_vec())
                    .collect();
                removed_pairs.push((
                    vec![v],
                    to_remove.iter().last().cloned().unwrap_or_default(),
                ));
                for key in to_remove {
                    current.simplices.remove(&key);
                }
                current.recompute_f_vector();
                steps += 1;
            }
            None => break,
        }
    }

    let is_collapsible = current.num_simplices() == 1 && current.num_vertices() == 1;

    CollapseResult {
        complex: current,
        steps,
        is_collapsible,
        removed_pairs,
    }
}

/// Check if a complex is a cone (has an apex vertex connected to all others).
///
/// Cones are always collapsible.
pub fn is_cone(complex: &SimplicialComplex) -> bool {
    let vertices: Vec<usize> = complex
        .simplices_of_dimension(0)
        .iter()
        .map(|s| s.vertices()[0])
        .collect();
    let edges: HashSet<(usize, usize)> = complex
        .simplices_of_dimension(1)
        .iter()
        .filter_map(|s| {
            let v = s.vertices();
            if v.len() == 2 {
                let (a, b) = (v[0].min(v[1]), v[0].max(v[1]));
                Some((a, b))
            } else {
                None
            }
        })
        .collect();

    for &apex in &vertices {
        let others: Vec<usize> = vertices.iter().copied().filter(|&v| v != apex).collect();
        if others.is_empty() {
            return true; // single vertex
        }
        let all_connected = others.iter().all(|&v| {
            let (a, b) = (apex.min(v), apex.max(v));
            edges.contains(&(a, b))
        });
        if all_connected {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collapse_triangle_removes_all() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        let result = collapse(&k);
        assert!(result.is_collapsible);
        assert!(result.steps > 0);
    }

    #[test]
    fn collapse_hollow_triangle_not_collapsible() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1]);
        k.add(vec![1, 2]);
        k.add(vec![0, 2]);
        let result = collapse(&k);
        // No free faces in a hollow triangle (each edge has 0 cofaces of dim 2)
        assert_eq!(result.steps, 0);
        assert!(!result.is_collapsible);
    }

    #[test]
    fn collapse_tetrahedron() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2, 3]);
        let result = collapse(&k);
        assert!(result.is_collapsible);
    }

    #[test]
    fn cone_is_detected() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1]);
        k.add(vec![0, 2]);
        k.add(vec![1, 2]);
        k.add(vec![0, 1, 2]);
        // This is a cone with apex 0? No - apex must connect to all others
        // Actually a cone: vertex 0 connected to 1 and 2, and all simplices present
        // Let's make a proper cone
        let mut cone_k = SimplicialComplex::new();
        cone_k.add(vec![0, 1, 2]);
        // cone_k is a triangle with all faces - that's a cone (vertex 0 is connected to 1 and 2)
        // A cone requires an apex connected to ALL other vertices
        assert!(is_cone(&cone_k));
    }

    #[test]
    fn strong_collapse_path() {
        // Path 0-1-2-3: vertex 3 is dominated by 2 (N[3]={2,3} ⊆ N[2]={1,2,3})
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1]);
        k.add(vec![1, 2]);
        k.add(vec![2, 3]);
        let result = strong_collapse(&k);
        assert!(result.steps > 0);
    }

    #[test]
    fn collapse_preserves_vertices_count_for_point() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0]);
        let result = collapse(&k);
        // Single vertex: nothing to collapse, but it's trivially collapsible
        // Actually with our definition, a single vertex is already "a vertex"
        // is_collapsible checks num_simplices==1 && num_vertices==1
        assert_eq!(result.complex.num_simplices(), 1);
    }

    #[test]
    fn collapse_two_triangles_sharing_edge() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        k.add(vec![1, 2, 3]);
        let result = collapse(&k);
        // Should be collapsible (both triangles share edge {1,2})
        assert!(result.is_collapsible);
    }

    #[test]
    fn collapse_step_count() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        let result = collapse(&k);
        // Triangle: remove 1 face+triangle, then 1 edge+edge, then 1 edge+edge, ...
        // Steps should be positive
        assert!(result.steps >= 1);
    }

    #[test]
    fn removed_pairs_recorded() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        let result = collapse(&k);
        assert!(!result.removed_pairs.is_empty());
    }
}
