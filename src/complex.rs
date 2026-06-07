//! Simplicial complex construction and queries.
//!
//! A **simplicial complex** K is a collection of simplices satisfying:
//!
//! 1. If σ ∈ K and τ ⊂ σ, then τ ∈ K (closure under face relation).
//! 2. The intersection of any two simplices in K is either empty or a face of both.
//!
//! This module provides [`SimplicialComplex`] — a concrete representation using
//! `Vec<Simplex>` with a `HashMap` index for fast lookups.

use crate::simplex::Simplex;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// A simplicial complex: a collection of simplices closed under the face relation.
///
/// Internally stores simplices indexed by dimension for efficient queries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimplicialComplex {
    /// All simplices, indexed by their sorted vertex key.
    pub(crate) simplices: HashMap<Vec<usize>, Simplex>,
    /// Number of simplices at each dimension.
    f_vector: Vec<usize>,
}

impl SimplicialComplex {
    /// Create an empty simplicial complex.
    pub fn new() -> Self {
        Self {
            simplices: HashMap::new(),
            f_vector: vec![],
        }
    }

    /// Add a simplex and all of its faces.
    ///
    /// This ensures the complex remains closed under the face relation.
    /// Adding `{0, 1, 2}` automatically adds `{0}`, `{1}`, `{2}`, `{0,1}`,
    /// `{0,2}`, and `{1,2}`.
    pub fn add_simplex(&mut self, simplex: &Simplex) {
        for face in simplex.all_faces() {
            let key = face.vertices().to_vec();
            let dim = face.dimension();
            self.simplices.insert(key, face);
            // Ensure f_vector is large enough
            if self.f_vector.len() <= dim {
                self.f_vector.resize(dim + 1, 0);
            }
        }
        self.recompute_f_vector();
    }

    /// Add a simplex from raw vertex indices.
    pub fn add(&mut self, vertices: Vec<usize>) {
        let s = Simplex::from_vertices(vertices);
        self.add_simplex(&s);
    }

    /// Check whether the stored collection is a valid simplicial complex
    /// (closed under face relation).
    pub fn is_simplicial_complex(&self) -> bool {
        for simplex in self.simplices.values() {
            for face in simplex.boundary() {
                if !self.contains(&face) {
                    return false;
                }
            }
        }
        true
    }

    /// Does the complex contain a simplex with the given vertices?
    pub fn contains(&self, simplex: &Simplex) -> bool {
        self.simplices.contains_key(simplex.vertices())
    }

    /// Number of simplices in the complex.
    pub fn num_simplices(&self) -> usize {
        self.simplices.len()
    }

    /// Iterate over all simplices.
    pub fn simplices(&self) -> impl Iterator<Item = &Simplex> {
        self.simplices.values()
    }

    /// Get simplices of a specific dimension.
    pub fn simplices_of_dimension(&self, dim: usize) -> Vec<&Simplex> {
        self.simplices
            .values()
            .filter(|s| s.dimension() == dim)
            .collect()
    }

    /// The **f-vector**: counts of simplices by dimension.
    ///
    /// `f[k]` = number of k-simplices.
    pub fn f_vector(&self) -> &[usize] {
        &self.f_vector
    }

    /// The **Euler characteristic**: χ = Σ_k (−1)^k · f_k.
    pub fn euler_characteristic(&self) -> i64 {
        let mut chi: i64 = 0;
        for (k, &count) in self.f_vector.iter().enumerate() {
            if k % 2 == 0 {
                chi += count as i64;
            } else {
                chi -= count as i64;
            }
        }
        chi
    }

    /// The **k-skeleton**: all simplices of dimension ≤ k.
    pub fn skeleton(&self, k: usize) -> SimplicialComplex {
        let mut skel = SimplicialComplex::new();
        for simplex in self.simplices.values() {
            if simplex.dimension() <= k {
                skel.simplices
                    .insert(simplex.vertices().to_vec(), simplex.clone());
            }
        }
        skel.recompute_f_vector();
        skel
    }

    /// Number of vertices (0-simplices).
    pub fn num_vertices(&self) -> usize {
        self.simplices_of_dimension(0).len()
    }

    /// Number of edges (1-simplices).
    pub fn num_edges(&self) -> usize {
        self.simplices_of_dimension(1).len()
    }

    /// Number of triangles (2-simplices).
    pub fn num_triangles(&self) -> usize {
        self.simplices_of_dimension(2).len()
    }

    /// Maximum dimension of any simplex in the complex.
    pub fn max_dimension(&self) -> usize {
        self.f_vector.len().saturating_sub(1)
    }

    /// Connected components via 1-skeleton (graph connectivity).
    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let vertices: HashSet<usize> = self
            .simplices_of_dimension(0)
            .iter()
            .map(|s| s.vertices()[0])
            .collect();

        let mut adjacency: HashMap<usize, HashSet<usize>> = HashMap::new();
        for &v in &vertices {
            adjacency.insert(v, HashSet::new());
        }
        for edge in self.simplices_of_dimension(1) {
            let vs = edge.vertices();
            if vs.len() == 2 {
                adjacency.get_mut(&vs[0]).unwrap().insert(vs[1]);
                adjacency.get_mut(&vs[1]).unwrap().insert(vs[0]);
            }
        }

        let mut visited: HashSet<usize> = HashSet::new();
        let mut components = Vec::new();

        for &start in &vertices {
            if visited.contains(&start) {
                continue;
            }
            let mut component = Vec::new();
            let mut stack = vec![start];
            while let Some(v) = stack.pop() {
                if visited.insert(v) {
                    component.push(v);
                    if let Some(neighbors) = adjacency.get(&v) {
                        for &n in neighbors {
                            if !visited.contains(&n) {
                                stack.push(n);
                            }
                        }
                    }
                }
            }
            components.push(component);
        }
        components
    }

    pub(crate) fn recompute_f_vector(&mut self) {
        if self.simplices.is_empty() {
            self.f_vector.clear();
            return;
        }
        let max_dim = self
            .simplices
            .values()
            .map(|s| s.dimension())
            .max()
            .unwrap_or(0);
        self.f_vector = vec![0usize; max_dim + 1];
        for simplex in self.simplices.values() {
            self.f_vector[simplex.dimension()] += 1;
        }
    }
}

impl Default for SimplicialComplex {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_complex() {
        let k = SimplicialComplex::new();
        assert_eq!(k.num_simplices(), 0);
        assert!(k.is_simplicial_complex());
    }

    #[test]
    fn add_triangle_closes_under_faces() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        // 3 vertices + 3 edges + 1 triangle = 7
        assert_eq!(k.num_simplices(), 7);
        assert!(k.is_simplicial_complex());
    }

    #[test]
    fn add_edge_only() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1]);
        assert_eq!(k.num_simplices(), 3); // 2 vertices + 1 edge
    }

    #[test]
    fn f_vector_triangle() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        assert_eq!(k.f_vector(), &[3, 3, 1]);
    }

    #[test]
    fn euler_characteristic_triangle() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        // χ = V - E + F = 3 - 3 + 1 = 1
        assert_eq!(k.euler_characteristic(), 1);
    }

    #[test]
    fn euler_characteristic_two_triangles() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        k.add(vec![1, 2, 3]);
        // vertices: 0,1,2,3 = 4; edges: {0,1},{0,2},{1,2},{1,3},{2,3} = 5; faces: 2
        // χ = 4 - 5 + 2 = 1
        assert_eq!(k.euler_characteristic(), 1);
    }

    #[test]
    fn skeleton_1_of_triangle() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        let skel = k.skeleton(1);
        // 3 vertices + 3 edges = 6
        assert_eq!(skel.num_simplices(), 6);
    }

    #[test]
    fn connected_components_single() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1]);
        k.add(vec![1, 2]);
        assert_eq!(k.connected_components().len(), 1);
    }

    #[test]
    fn connected_components_disconnected() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1]);
        k.add(vec![2, 3]);
        assert_eq!(k.connected_components().len(), 2);
    }

    #[test]
    fn contains_simplex() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        assert!(k.contains(&Simplex::from_vertices(vec![0, 1])));
        assert!(!k.contains(&Simplex::from_vertices(vec![0, 3])));
    }

    #[test]
    fn max_dimension() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2, 3]);
        assert_eq!(k.max_dimension(), 3);
    }
}
