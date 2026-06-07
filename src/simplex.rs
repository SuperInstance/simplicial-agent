//! Individual simplex operations.
//!
//! A **simplex** is the most basic building block of a simplicial complex.
//! An *n*-simplex is defined by an ordered set of *n + 1* vertices:
//!
//! - 0-simplex = point (single vertex)
//! - 1-simplex = edge (two vertices)
//! - 2-simplex = triangle (three vertices)
//! - 3-simplex = tetrahedron (four vertices)
//!
//! # Conventions
//!
//! Vertices are stored as sorted `Vec<usize>`. Construction normalises the
//! vertex order and deduplicates.

use serde::{Deserialize, Serialize};
use std::collections::HashSet;

/// A simplex represented by a sorted set of vertex indices.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Simplex {
    /// Sorted, deduplicated vertex indices.
    pub(crate) vertices: Vec<usize>,
}

impl Simplex {
    /// Construct a simplex from a vector of vertex indices.
    ///
    /// Vertices are sorted and deduplicated, so `{2, 0, 1}` and `{0, 1, 2}`
    /// produce the same simplex.
    ///
    /// # Panics
    ///
    /// Panics if `vertices` is empty.
    ///
    /// # Example
    ///
    /// ```
    /// use simplicial_agent::simplex::Simplex;
    ///
    /// let s = Simplex::from_vertices(vec![2, 0, 1]);
    /// assert_eq!(s.vertices(), &[0, 1, 2]);
    /// assert_eq!(s.dimension(), 2); // triangle
    /// ```
    pub fn from_vertices(vertices: Vec<usize>) -> Self {
        assert!(
            !vertices.is_empty(),
            "a simplex must have at least one vertex"
        );
        let mut v: Vec<usize> = vertices
            .into_iter()
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        v.sort_unstable();
        Self { vertices: v }
    }

    /// Return a reference to the sorted vertex list.
    pub fn vertices(&self) -> &[usize] {
        &self.vertices
    }

    /// Dimension of the simplex: `|vertices| - 1`.
    ///
    /// A single vertex has dimension 0, an edge has dimension 1, etc.
    pub fn dimension(&self) -> usize {
        self.vertices.len().saturating_sub(1)
    }

    /// Number of vertices.
    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    /// True if the simplex has no vertices (should never happen after construction).
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Is `self` a face of `other`?
    ///
    /// σ is a face of τ iff every vertex of σ is also a vertex of τ.
    /// Every simplex is a face of itself (reflexive).
    pub fn is_face_of(&self, other: &Simplex) -> bool {
        let other_set: HashSet<_> = other.vertices.iter().copied().collect();
        self.vertices.iter().all(|v| other_set.contains(v))
    }

    /// Is `self` a **proper** face of `other`? (face and not equal)
    pub fn is_proper_face_of(&self, other: &Simplex) -> bool {
        self != other && self.is_face_of(other)
    }

    /// Boundary: all (n−1)-dimensional faces.
    ///
    /// For a simplex `{v₀, v₁, …, vₙ}`, the boundary consists of the n+1
    /// faces obtained by removing exactly one vertex.
    ///
    /// Returns an empty vec for a 0-simplex (point).
    pub fn boundary(&self) -> Vec<Simplex> {
        if self.vertices.len() <= 1 {
            return vec![];
        }
        let mut faces = Vec::with_capacity(self.vertices.len());
        for i in 0..self.vertices.len() {
            let face_verts: Vec<usize> = self
                .vertices
                .iter()
                .enumerate()
                .filter_map(|(j, &v)| if j != i { Some(v) } else { None })
                .collect();
            faces.push(Simplex::from_vertices(face_verts));
        }
        faces
    }

    /// All faces of dimension exactly `k`.
    ///
    /// Uses combinatorial selection without recursion.
    pub fn faces_of_dimension(&self, k: usize) -> Vec<Simplex> {
        if k > self.dimension() || self.vertices.len() < k + 1 {
            return vec![];
        }
        let _choose = self.vertices.len() - (k + 1);
        self.combinations(k + 1)
    }

    /// All faces of any dimension (including self).
    pub fn all_faces(&self) -> Vec<Simplex> {
        let mut result = Vec::new();
        for size in 1..=self.vertices.len() {
            result.extend(self.combinations(size));
        }
        result
    }

    /// Combinatorial (n choose k) selections of vertices, iterative.
    fn combinations(&self, k: usize) -> Vec<Simplex> {
        if k == 0 || k > self.vertices.len() {
            return vec![];
        }
        let n = self.vertices.len();
        let mut result = Vec::new();

        // Iterative combination generation using index tracking
        let mut indices: Vec<usize> = (0..k).collect();
        loop {
            let combo: Vec<usize> = indices.iter().map(|&i| self.vertices[i]).collect();
            result.push(Simplex::from_vertices(combo));

            // Find the rightmost index that can be incremented
            let mut i = k as i32 - 1;
            while i >= 0 {
                let ii = i as usize;
                if indices[ii] < n - k + ii + 1 - 1 {
                    // can increment
                    break;
                }
                i -= 1;
            }
            if i < 0 {
                break;
            }
            let ii = i as usize;
            indices[ii] += 1;
            for j in (ii + 1)..k {
                indices[j] = indices[j - 1] + 1;
            }
        }
        result
    }
}

/// Star of a simplex σ in a collection of simplices: all simplices containing σ.
pub fn star(sigma: &Simplex, simplices: &[Simplex]) -> Vec<Simplex> {
    simplices
        .iter()
        .filter(|s| sigma.is_face_of(s))
        .cloned()
        .collect()
}

/// Link of a simplex σ: faces of simplices in the star that do not intersect σ.
///
/// Lk(σ) = { τ ∈ Cl(St(σ)) : τ ∩ σ = ∅ }
pub fn link(sigma: &Simplex, simplices: &[Simplex]) -> Vec<Simplex> {
    let sigma_verts: HashSet<_> = sigma.vertices.iter().copied().collect();
    let st = star(sigma, simplices);
    // Closure of the star: all faces of all simplices in the star
    let mut closure_set: HashSet<Simplex> = HashSet::new();
    for s in &st {
        for face in s.all_faces() {
            closure_set.insert(face);
        }
    }
    closure_set
        .into_iter()
        .filter(|tau| {
            let tau_verts: HashSet<_> = tau.vertices.iter().copied().collect();
            tau_verts.intersection(&sigma_verts).count() == 0
        })
        .collect()
}

/// Closure of a set of simplices: the set plus all faces of every simplex.
pub fn closure(simplices: &[Simplex]) -> Vec<Simplex> {
    let mut set: HashSet<Simplex> = HashSet::new();
    for s in simplices {
        for face in s.all_faces() {
            set.insert(face);
        }
    }
    set.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simplex_from_vertices_sorts() {
        let s = Simplex::from_vertices(vec![3, 1, 2]);
        assert_eq!(s.vertices(), &[1, 2, 3]);
    }

    #[test]
    fn simplex_deduplicates() {
        let s = Simplex::from_vertices(vec![1, 1, 2]);
        assert_eq!(s.vertices(), &[1, 2]);
    }

    #[test]
    fn dimension_point() {
        let s = Simplex::from_vertices(vec![0]);
        assert_eq!(s.dimension(), 0);
    }

    #[test]
    fn dimension_triangle() {
        let s = Simplex::from_vertices(vec![0, 1, 2]);
        assert_eq!(s.dimension(), 2);
    }

    #[test]
    #[should_panic]
    fn empty_simplex_panics() {
        let _ = Simplex::from_vertices(vec![]);
    }

    #[test]
    fn face_relation_reflexive() {
        let s = Simplex::from_vertices(vec![0, 1]);
        assert!(s.is_face_of(&s));
    }

    #[test]
    fn edge_face_of_triangle() {
        let edge = Simplex::from_vertices(vec![0, 1]);
        let tri = Simplex::from_vertices(vec![0, 1, 2]);
        assert!(edge.is_face_of(&tri));
    }

    #[test]
    fn triangle_not_face_of_edge() {
        let edge = Simplex::from_vertices(vec![0, 1]);
        let tri = Simplex::from_vertices(vec![0, 1, 2]);
        assert!(!tri.is_face_of(&edge));
    }

    #[test]
    fn proper_face_excludes_self() {
        let s = Simplex::from_vertices(vec![0, 1]);
        assert!(!s.is_proper_face_of(&s));
    }

    #[test]
    fn boundary_of_triangle() {
        let tri = Simplex::from_vertices(vec![0, 1, 2]);
        let bd = tri.boundary();
        assert_eq!(bd.len(), 3);
        assert!(bd.contains(&Simplex::from_vertices(vec![0, 1])));
        assert!(bd.contains(&Simplex::from_vertices(vec![0, 2])));
        assert!(bd.contains(&Simplex::from_vertices(vec![1, 2])));
    }

    #[test]
    fn boundary_of_edge() {
        let edge = Simplex::from_vertices(vec![0, 1]);
        let bd = edge.boundary();
        assert_eq!(bd.len(), 2);
        assert!(bd.contains(&Simplex::from_vertices(vec![0])));
        assert!(bd.contains(&Simplex::from_vertices(vec![1])));
    }

    #[test]
    fn boundary_of_point_empty() {
        let pt = Simplex::from_vertices(vec![0]);
        assert!(pt.boundary().is_empty());
    }

    #[test]
    fn faces_of_dimension_edge_in_triangle() {
        let tri = Simplex::from_vertices(vec![0, 1, 2]);
        let edges = tri.faces_of_dimension(1);
        assert_eq!(edges.len(), 3);
    }

    #[test]
    fn all_faces_of_triangle() {
        let tri = Simplex::from_vertices(vec![0, 1, 2]);
        let faces = tri.all_faces();
        // 3 vertices + 3 edges + 1 triangle = 7
        assert_eq!(faces.len(), 7);
    }

    #[test]
    fn star_operation() {
        let edge = Simplex::from_vertices(vec![0, 1]);
        let simplices = vec![
            Simplex::from_vertices(vec![0]),
            Simplex::from_vertices(vec![1]),
            Simplex::from_vertices(vec![2]),
            Simplex::from_vertices(vec![0, 1]),
            Simplex::from_vertices(vec![1, 2]),
            Simplex::from_vertices(vec![0, 1, 2]),
        ];
        let s = star(&edge, &simplices);
        assert_eq!(s.len(), 2); // {0,1} and {0,1,2}
    }

    #[test]
    fn closure_adds_all_faces() {
        let tri = Simplex::from_vertices(vec![0, 1, 2]);
        let cl = closure(&[tri.clone()]);
        assert_eq!(cl.len(), 7);
    }

    #[test]
    fn combinations_size() {
        let s = Simplex::from_vertices(vec![0, 1, 2, 3]);
        // 4 choose 2 = 6
        assert_eq!(s.combinations(2).len(), 6);
        // 4 choose 3 = 4
        assert_eq!(s.combinations(3).len(), 4);
    }
}
