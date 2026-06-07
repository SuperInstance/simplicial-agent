//! Agent collaboration topology analysis.
//!
//! Models agents as vertices and collaboration groups as simplices in a
//! simplicial complex. Higher-dimensional simplices represent larger
//! collaboration groups.
//!
//! ## Overview
//!
//! - **0-simplices**: Individual agents.
//! - **1-simplices**: Pairwise collaborations.
//! - **2-simplices**: Three-way collaboration.
//! - **k-simplices**: (k+1)-agent collaboration.
//!
//! This lets us analyze collaboration topology using tools from algebraic
//! topology: connected components, holes (missing collaborations), coverage,
//! and collaboration dynamics via the simplicial Laplacian.

use crate::complex::SimplicialComplex;
use crate::euler::{betti_numbers, euler_characteristic};
use crate::simplex::Simplex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A single collaboration record: which agents participated and how strongly.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CollaborationRecord {
    /// Agent indices participating in this collaboration.
    pub agents: Vec<usize>,
    /// Interaction strength (0.0 to 1.0).
    pub interaction_strength: f64,
}

/// Result of collaboration topology analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CollaborationTopology {
    /// The collaboration simplicial complex.
    pub complex: SimplicialComplex,
    /// Number of connected components.
    pub connected_components: usize,
    /// Euler characteristic.
    pub euler_characteristic: i64,
    /// Betti numbers.
    pub betti_numbers: Vec<usize>,
    /// Collaboration coverage: fraction of possible pairwise collaborations realized.
    pub coverage: f64,
    /// Missing collaborations (holes): pairs of agents in the same component
    /// who are not directly connected.
    pub holes: Vec<(usize, usize)>,
    /// Simplicial Laplacian eigenvalues (approximation).
    pub laplacian_eigenvalues: Vec<f64>,
}

/// Agent collaboration complex builder.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentCollaboration {
    /// Collaboration records.
    records: Vec<CollaborationRecord>,
    /// Number of agents.
    num_agents: usize,
}

impl AgentCollaboration {
    /// Create a new collaboration model for `num_agents` agents.
    pub fn new(num_agents: usize) -> Self {
        Self {
            records: Vec::new(),
            num_agents,
        }
    }

    /// Add a collaboration record.
    pub fn add_collaboration(&mut self, agents: Vec<usize>, interaction_strength: f64) {
        assert!(
            (0.0..=1.0).contains(&interaction_strength),
            "interaction_strength must be in [0, 1]"
        );
        self.records.push(CollaborationRecord {
            agents,
            interaction_strength,
        });
    }

    /// Build the collaboration simplicial complex with a strength threshold.
    ///
    /// Only collaborations with `interaction_strength >= threshold` are included.
    pub fn build_complex(&self, threshold: f64) -> SimplicialComplex {
        let mut complex = SimplicialComplex::new();

        // Add all agents as vertices
        for i in 0..self.num_agents {
            complex.add(vec![i]);
        }

        // Add qualifying collaborations as simplices
        for record in &self.records {
            if record.interaction_strength >= threshold {
                complex.add(record.agents.clone());
            }
        }

        complex
    }

    /// Analyze the collaboration topology at the given threshold.
    pub fn analyze(&self, threshold: f64) -> CollaborationTopology {
        let complex = self.build_complex(threshold);
        let components = complex.connected_components();
        let betti = betti_numbers(&complex);
        let chi = euler_characteristic(&complex);

        // Compute coverage: realized pairs / possible pairs
        let n = self.num_agents;
        let possible_pairs = if n >= 2 { n * (n - 1) / 2 } else { 0 };
        let realized_pairs = complex.num_edges();
        let coverage = if possible_pairs > 0 {
            realized_pairs as f64 / possible_pairs as f64
        } else {
            1.0
        };

        // Find holes: agents in same component but not directly connected
        let mut holes = Vec::new();
        for component in &components {
            if component.len() < 2 {
                continue;
            }
            for i in 0..component.len() {
                for j in (i + 1)..component.len() {
                    let a = component[i];
                    let b = component[j];
                    let edge = Simplex::from_vertices(vec![a, b]);
                    if !complex.contains(&edge) {
                        holes.push((a.min(b), a.max(b)));
                    }
                }
            }
        }
        holes.sort_unstable();
        holes.dedup();

        // Compute simplicial Laplacian eigenvalues (1-skeleton / graph Laplacian)
        let laplacian_eigenvalues = compute_laplacian_eigenvalues(&complex);

        CollaborationTopology {
            complex,
            connected_components: components.len(),
            euler_characteristic: chi,
            betti_numbers: betti,
            coverage,
            holes,
            laplacian_eigenvalues,
        }
    }
}

/// Compute approximate eigenvalues of the 1-skeleton graph Laplacian.
///
/// The graph Laplacian L = D - A where D is the degree matrix and A is the
/// adjacency matrix. For small graphs, we compute eigenvalues using Jacobi
/// eigenvalue algorithm for symmetric matrices.
///
/// For larger graphs, returns a subset of extreme eigenvalues.
#[allow(clippy::needless_range_loop)]
fn compute_laplacian_eigenvalues(complex: &SimplicialComplex) -> Vec<f64> {
    let n = complex.num_vertices();
    if n == 0 {
        return vec![];
    }

    // Map vertex indices to 0..n-1
    let vertices: Vec<usize> = complex
        .simplices_of_dimension(0)
        .iter()
        .map(|s| s.vertices()[0])
        .collect();

    let mut idx_map: HashMap<usize, usize> = HashMap::new();
    for (i, &v) in vertices.iter().enumerate() {
        idx_map.insert(v, i);
    }

    // Build adjacency matrix
    let mut adj = vec![vec![0.0f64; n]; n];
    let mut degree = vec![0.0f64; n];

    for edge in complex.simplices_of_dimension(1) {
        let vs = edge.vertices();
        if vs.len() == 2
            && let (Some(&i), Some(&j)) = (idx_map.get(&vs[0]), idx_map.get(&vs[1]))
        {
            adj[i][j] = 1.0;
            adj[j][i] = 1.0;
            degree[i] += 1.0;
            degree[j] += 1.0;
        }
    }

    // Laplacian: L[i][j] = degree[i] if i==j, else -adj[i][j]
    let mut lap: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            lap[i][j] = if i == j { degree[i] } else { -adj[i][j] };
        }
    }

    // Compute eigenvalues via QR-like approach (simple power method for top eigenvalue)
    // For small matrices, just compute trace and frobenius norm as proxies
    if n <= 1 {
        return vec![0.0];
    }

    // Compute eigenvalues using iterative deflation
    let mut eigenvalues = Vec::new();
    let dim_remaining = n;
    let _remaining_lap = lap.clone();

    // Simple approach: compute characteristic polynomial eigenvalues
    // For practical purposes, return the diagonal of the Schur-like form
    // Use Jacobi-like iteration for symmetric matrices
    let mut mat = lap;
    for _ in 0..100 {
        // Find the largest off-diagonal element
        let mut max_val = 0.0f64;
        let mut max_i = 0;
        let mut max_j = 1;
        for i in 0..dim_remaining {
            for j in (i + 1)..dim_remaining {
                if mat[i][j].abs() > max_val {
                    max_val = mat[i][j].abs();
                    max_i = i;
                    max_j = j;
                }
            }
        }
        if max_val < 1e-10 {
            break;
        }

        // Jacobi rotation
        let a_ii = mat[max_i][max_i];
        let a_jj = mat[max_j][max_j];
        let a_ij = mat[max_i][max_j];
        let theta = if (a_ii - a_jj).abs() < 1e-15 {
            std::f64::consts::FRAC_PI_4
        } else {
            0.5 * (2.0 * a_ij / (a_ii - a_jj)).atan()
        };
        let c = theta.cos();
        let s = theta.sin();

        for k in 0..dim_remaining {
            if k == max_i || k == max_j {
                continue;
            }
            let new_ik = c * mat[max_i][k] + s * mat[max_j][k];
            let new_jk = -s * mat[max_i][k] + c * mat[max_j][k];
            mat[max_i][k] = new_ik;
            mat[k][max_i] = new_ik;
            mat[max_j][k] = new_jk;
            mat[k][max_j] = new_jk;
        }

        let new_ii = c * c * a_ii + 2.0 * s * c * a_ij + s * s * a_jj;
        let new_jj = s * s * a_ii - 2.0 * s * c * a_ij + c * c * a_jj;
        mat[max_i][max_i] = new_ii;
        mat[max_j][max_j] = new_jj;
        mat[max_i][max_j] = 0.0;
        mat[max_j][max_i] = 0.0;
    }

    for i in 0..dim_remaining {
        eigenvalues.push(mat[i][i]);
    }
    eigenvalues.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

    eigenvalues
}

/// Build a collaboration complex from raw records.
///
/// Convenience function: takes `Vec<(agents, strength)>` tuples.
pub fn collaboration_complex(
    num_agents: usize,
    records: &[(Vec<usize>, f64)],
    threshold: f64,
) -> SimplicialComplex {
    let mut collab = AgentCollaboration::new(num_agents);
    for (agents, strength) in records {
        collab.add_collaboration(agents.clone(), *strength);
    }
    collab.build_complex(threshold)
}

/// Analyze collaboration topology from raw records.
pub fn analyze_collaboration(
    num_agents: usize,
    records: &[(Vec<usize>, f64)],
    threshold: f64,
) -> CollaborationTopology {
    let mut collab = AgentCollaboration::new(num_agents);
    for (agents, strength) in records {
        collab.add_collaboration(agents.clone(), *strength);
    }
    collab.analyze(threshold)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_collaboration_complex() {
        let records = vec![(vec![0, 1], 0.9), (vec![1, 2], 0.8), (vec![0, 2], 0.7)];
        let complex = collaboration_complex(3, &records, 0.5);
        assert_eq!(complex.num_vertices(), 3);
        assert_eq!(complex.num_edges(), 3);
    }

    #[test]
    fn threshold_filters_collaborations() {
        let records = vec![(vec![0, 1], 0.9), (vec![1, 2], 0.3), (vec![0, 2], 0.7)];
        let complex = collaboration_complex(3, &records, 0.5);
        assert_eq!(complex.num_edges(), 2);
    }

    #[test]
    fn analyze_topology() {
        let records = vec![(vec![0, 1], 0.9), (vec![1, 2], 0.8), (vec![0, 2], 0.7)];
        let topo = analyze_collaboration(3, &records, 0.5);
        assert_eq!(topo.connected_components, 1);
        assert_eq!(topo.holes.len(), 0);
        assert!(topo.coverage > 0.99); // 3/3 edges
    }

    #[test]
    fn holes_detected() {
        // 4 agents in a line: 0-1-2-3, missing 0-2, 0-3, 1-3
        let records = vec![(vec![0, 1], 0.9), (vec![1, 2], 0.8), (vec![2, 3], 0.7)];
        let topo = analyze_collaboration(4, &records, 0.5);
        assert_eq!(topo.connected_components, 1);
        assert_eq!(topo.holes.len(), 3); // 0-2, 0-3, 1-3
    }

    #[test]
    fn coverage_calculation() {
        // 4 agents, 6 possible pairs, 3 realized
        let records = vec![(vec![0, 1], 0.9), (vec![1, 2], 0.8), (vec![2, 3], 0.7)];
        let topo = analyze_collaboration(4, &records, 0.5);
        assert!((topo.coverage - 0.5).abs() < 0.01); // 3/6 = 0.5
    }

    #[test]
    fn disconnected_components() {
        let records = vec![(vec![0, 1], 0.9), (vec![2, 3], 0.8)];
        let topo = analyze_collaboration(4, &records, 0.5);
        assert_eq!(topo.connected_components, 2);
    }

    #[test]
    fn laplacian_eigenvalues_connected() {
        let records = vec![(vec![0, 1], 0.9), (vec![1, 2], 0.8), (vec![0, 2], 0.7)];
        let topo = analyze_collaboration(3, &records, 0.5);
        // Connected graph: smallest eigenvalue should be ~0
        assert!(topo.laplacian_eigenvalues[0].abs() < 0.1);
        // Second smallest > 0 (algebraic connectivity)
        assert!(topo.laplacian_eigenvalues[1] > 0.1);
    }

    #[test]
    fn three_way_collaboration() {
        let records = vec![(vec![0, 1, 2], 0.9)];
        let complex = collaboration_complex(3, &records, 0.5);
        assert_eq!(complex.num_triangles(), 1);
    }

    #[test]
    fn mixed_dimensionality() {
        let records = vec![
            (vec![0, 1], 0.9),
            (vec![1, 2], 0.8),
            (vec![0, 1, 2], 0.7),
            (vec![2, 3], 0.6),
        ];
        let topo = analyze_collaboration(4, &records, 0.5);
        assert_eq!(topo.connected_components, 1);
    }
}
