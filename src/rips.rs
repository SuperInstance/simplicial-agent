//! Vietoris-Rips complex construction from distance matrices.
//!
//! Given a finite metric space (X, d) and a threshold ε > 0, the
//! **Vietoris-Rips complex** VR(X, ε) is the simplicial complex whose
//! simplices are all finite subsets σ ⊂ X with diam(σ) ≤ ε.
//!
//! In practice we build it by:
//! 1. Adding all vertices.
//! 2. Adding edges for all pairs (i, j) with d(i, j) ≤ ε.
//! 3. Finding all cliques in the resulting graph (iterative worklist).
//!
//! # Clique Detection
//!
//! Cliques are found **iteratively** using a worklist + `HashSet`,
//! NOT recursive backtracking. Each worklist entry is a clique (sorted
//! vertex set). We extend each entry by adding vertices that are adjacent
//! to all current members.

use crate::complex::SimplicialComplex;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Configuration for building a Vietoris-Rips complex.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RipsConfig {
    /// Distance threshold ε. Pairs with distance ≤ ε become edges.
    pub epsilon: f64,
    /// Maximum simplex dimension to construct. Limits clique search.
    pub max_dimension: Option<usize>,
}

impl RipsConfig {
    /// Create config with given epsilon and no dimension limit.
    pub fn new(epsilon: f64) -> Self {
        Self {
            epsilon,
            max_dimension: None,
        }
    }

    /// Set maximum dimension.
    pub fn with_max_dimension(mut self, dim: usize) -> Self {
        self.max_dimension = Some(dim);
        self
    }
}

/// Builder for Vietoris-Rips complexes.
pub struct VietorisRipsBuilder {
    config: RipsConfig,
}

impl VietorisRipsBuilder {
    /// Create a builder with the given configuration.
    pub fn new(config: RipsConfig) -> Self {
        Self { config }
    }

    /// Build the Vietoris-Rips complex from a distance matrix.
    ///
    /// `distances[i][j]` is the distance between points i and j.
    /// The matrix must be square and symmetric with zeros on the diagonal.
    pub fn build(&self, distances: &[Vec<f64>]) -> SimplicialComplex {
        let n = distances.len();
        assert!(n > 0, "distance matrix must be non-empty");

        let mut complex = SimplicialComplex::new();

        // Step 1: Add all vertices
        for i in 0..n {
            complex.add(vec![i]);
        }

        // Step 2: Build adjacency from edges within epsilon
        let mut adjacency: HashMap<usize, HashSet<usize>> = HashMap::new();
        for i in 0..n {
            adjacency.insert(i, HashSet::new());
        }
        for (i, row) in distances.iter().enumerate().take(n) {
            for (j, &d) in row.iter().enumerate().take(n).skip(i + 1) {
                if d <= self.config.epsilon {
                    adjacency.get_mut(&i).unwrap().insert(j);
                    adjacency.get_mut(&j).unwrap().insert(i);
                }
            }
        }

        // Step 3: Iterative clique detection using worklist
        let max_clique_size = self
            .config
            .max_dimension
            .map(|d| d + 1)
            .unwrap_or(usize::MAX);

        // Worklist entries: sorted vertex sets (cliques)
        // Start with all edges as initial 2-cliques
        let mut worklist: Vec<Vec<usize>> = Vec::new();
        for i in 0..n {
            if let Some(neighbors) = adjacency.get(&i) {
                for &j in neighbors {
                    if j > i {
                        worklist.push(vec![i, j]);
                    }
                }
            }
        }

        let mut all_cliques: HashSet<Vec<usize>> = HashSet::new();
        // Track vertex cliques separately (they're always present)
        for i in 0..n {
            all_cliques.insert(vec![i]);
        }
        // Add all edge cliques
        for clique in &worklist {
            all_cliques.insert(clique.clone());
        }

        // Iteratively extend cliques
        while let Some(clique) = worklist.pop() {
            if clique.len() >= max_clique_size {
                continue;
            }

            // The last vertex in the clique
            let last = *clique.last().unwrap();

            // Find common neighbors of all clique vertices that are > last vertex
            // (to avoid duplicates, we only extend with higher-numbered vertices)
            let mut candidates: Vec<usize> = adjacency
                .get(&last)
                .map(|s| s.iter().copied().filter(|&v| v > last).collect())
                .unwrap_or_default();

            for &v in &clique {
                if v == last {
                    continue;
                }
                if let Some(neighbors) = adjacency.get(&v) {
                    candidates.retain(|&c| neighbors.contains(&c));
                }
            }

            for c in candidates {
                let mut new_clique = clique.clone();
                new_clique.push(c);
                new_clique.sort_unstable();
                if !all_cliques.contains(&new_clique) {
                    all_cliques.insert(new_clique.clone());
                    worklist.push(new_clique);
                }
            }
        }

        // Add all cliques as simplices
        for clique in &all_cliques {
            complex.add(clique.clone());
        }

        complex
    }
}

/// Convenience function: build a Vietoris-Rips complex from a distance matrix.
pub fn rips_complex(distances: &[Vec<f64>], epsilon: f64) -> SimplicialComplex {
    VietorisRipsBuilder::new(RipsConfig::new(epsilon)).build(distances)
}

/// Build a Rips complex with a maximum dimension constraint.
pub fn rips_complex_bounded(
    distances: &[Vec<f64>],
    epsilon: f64,
    max_dim: usize,
) -> SimplicialComplex {
    VietorisRipsBuilder::new(RipsConfig::new(epsilon).with_max_dimension(max_dim)).build(distances)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_distance_matrix() -> Vec<Vec<f64>> {
        // 4 points: 0-1-2-3 in a line, distance = index difference
        vec![
            vec![0.0, 1.0, 2.0, 3.0],
            vec![1.0, 0.0, 1.0, 2.0],
            vec![2.0, 1.0, 0.0, 1.0],
            vec![3.0, 2.0, 1.0, 0.0],
        ]
    }

    #[test]
    fn rips_epsilon_1_line() {
        let dist = unit_distance_matrix();
        let complex = rips_complex(&dist, 1.0);
        // Edges: {0,1}, {1,2}, {2,3}
        assert_eq!(complex.num_edges(), 3);
        // No triangles (no 3 points all within distance 1)
        assert_eq!(complex.num_triangles(), 0);
    }

    #[test]
    fn rips_epsilon_2_line() {
        let dist = unit_distance_matrix();
        let complex = rips_complex(&dist, 2.0);
        // Edges: all adjacent + {0,2}, {1,3}
        assert_eq!(complex.num_edges(), 5);
        // Triangles: {0,1,2}, {1,2,3}
        assert_eq!(complex.num_triangles(), 2);
    }

    #[test]
    fn rips_epsilon_3_line() {
        let dist = unit_distance_matrix();
        let complex = rips_complex(&dist, 3.0);
        // All pairs connected → tetrahedron
        assert_eq!(complex.num_triangles(), 4);
    }

    #[test]
    fn rips_clique_complete_graph() {
        // 3 points all at distance 0.5
        let dist = vec![
            vec![0.0, 0.5, 0.5],
            vec![0.5, 0.0, 0.5],
            vec![0.5, 0.5, 0.0],
        ];
        let complex = rips_complex(&dist, 1.0);
        assert_eq!(complex.num_triangles(), 1);
        assert_eq!(complex.num_edges(), 3);
    }

    #[test]
    fn rips_bounded_dimension() {
        let dist = vec![
            vec![0.0, 0.5, 0.5],
            vec![0.5, 0.0, 0.5],
            vec![0.5, 0.5, 0.0],
        ];
        let complex = rips_complex_bounded(&dist, 1.0, 1);
        // No triangles, only edges
        assert_eq!(complex.num_triangles(), 0);
        assert_eq!(complex.num_edges(), 3);
    }

    #[test]
    fn rips_single_point() {
        let dist = vec![vec![0.0]];
        let complex = rips_complex(&dist, 1.0);
        assert_eq!(complex.num_simplices(), 1);
        assert_eq!(complex.num_vertices(), 1);
    }

    #[test]
    fn rips_two_points_within_epsilon() {
        let dist = vec![vec![0.0, 0.5], vec![0.5, 0.0]];
        let complex = rips_complex(&dist, 1.0);
        assert_eq!(complex.num_edges(), 1);
    }

    #[test]
    fn rips_two_points_outside_epsilon() {
        let dist = vec![vec![0.0, 2.0], vec![2.0, 0.0]];
        let complex = rips_complex(&dist, 1.0);
        assert_eq!(complex.num_edges(), 0);
    }
}
