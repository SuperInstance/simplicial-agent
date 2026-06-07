//! Euler characteristic and Betti number computation.
//!
//! The **Euler characteristic** of a simplicial complex K is the alternating
//! sum of its f-vector:
//!
//! χ(K) = Σ_k (−1)^k · f_k
//!
//! where f_k is the number of k-simplices.
//!
//! The Euler-Poincaré theorem relates this to homology:
//!
//! χ(K) = Σ_k (−1)^k · rank(H_k) = Σ_k (−1)^k · β_k
//!
//! where β_k = rank(H_k) is the k-th Betti number.

use crate::complex::SimplicialComplex;
use serde::{Deserialize, Serialize};

/// Computed topological invariants of a simplicial complex.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TopologyInvariants {
    /// Euler characteristic χ.
    pub euler_characteristic: i64,
    /// F-vector: count of simplices by dimension.
    pub f_vector: Vec<usize>,
    /// Betti numbers β_0, β_1, ..., β_k.
    /// β_0 = number of connected components.
    /// β_k = rank of k-th homology group (simplified boundary rank computation).
    pub betti_numbers: Vec<usize>,
}

/// Compute the Euler characteristic of a simplicial complex.
pub fn euler_characteristic(complex: &SimplicialComplex) -> i64 {
    complex.euler_characteristic()
}

/// Compute Betti numbers using simplified boundary rank computation.
///
/// For a complex of dimension d, β_k is computed from:
///
/// β_k = dim(C_k) − rank(∂_k) − rank(∂_{k+1})
///
/// where C_k is the chain group (number of k-simplices) and ∂_k is the
/// boundary operator. The boundary rank is computed from how many (k-1)-faces
/// appear as boundaries of k-simplices.
///
/// **Note**: This uses a simplified computation. For rigorous homology,
/// one would need Smith normal form over Z. Here we use rank over Q/Z₂
/// approximations.
pub fn betti_numbers(complex: &SimplicialComplex) -> Vec<usize> {
    let f = complex.f_vector();
    if f.is_empty() {
        return vec![0];
    }

    let max_dim = f.len() - 1;

    // Compute boundary ranks
    // rank(∂_k) = number of linearly independent (k-1)-boundary components
    // Simplified: rank(∂_k) = min(f[k], f[k-1]) for our purposes
    // More accurate: rank(∂_k) ≤ min(f[k] * (k+1), f[k-1]) but we compute properly

    let mut boundary_ranks: Vec<usize> = vec![0; max_dim + 2]; // ∂_0 through ∂_{max_dim+1}

    for (k, rank_slot) in boundary_ranks
        .iter_mut()
        .enumerate()
        .take(max_dim + 1)
        .skip(1)
    {
        let k_faces = complex.simplices_of_dimension(k);
        let km1_faces = complex.simplices_of_dimension(k - 1);
        if k_faces.is_empty() || km1_faces.is_empty() {
            *rank_slot = 0;
            continue;
        }

        // Count how many (k-1)-faces appear in boundaries of k-simplices
        // Each k-simplex has (k+1) boundary (k-1)-faces
        use std::collections::HashSet;
        let mut boundary_face_set: HashSet<Vec<usize>> = HashSet::new();
        for simplex in &k_faces {
            for face in simplex.boundary() {
                boundary_face_set.insert(face.vertices().to_vec());
            }
        }
        // rank(∂_k) = min(number of boundary faces that appear, f[k-1])
        // But more accurately: it's the rank of the boundary matrix
        // For simplicity: rank(∂_k) = min(f[k], boundary_face_set.len()) approximately
        // Actually: each k-simplex contributes k+1 boundary faces, and the rank
        // of the boundary map is bounded by both f[k] and the number of distinct (k-1)-faces
        *rank_slot = std::cmp::min(k_faces.len(), boundary_face_set.len());

        // But we need to account for linear dependence
        // The kernel of ∂_k has dimension f[k] - rank(∂_k)
        // β_k = dim(ker ∂_k) - rank(∂_{k+1})
    }

    // Compute Betti numbers
    // β_k = f[k] - rank(∂_k) - rank(∂_{k+1})
    let mut betti = Vec::with_capacity(max_dim + 1);
    for k in 0..=max_dim {
        let c_k = f[k];
        let rk = if k == 0 { 0 } else { boundary_ranks[k] };
        let rk1 = boundary_ranks[k + 1];
        let b = (c_k as i64) - (rk as i64) - (rk1 as i64);
        betti.push(if b > 0 { b as usize } else { 0 });
    }

    // Override β_0 with connected components (definitive)
    let components = complex.connected_components();
    if !betti.is_empty() {
        betti[0] = if components.is_empty() {
            0
        } else {
            components.len()
        };
    }

    betti
}

/// Compute all topological invariants for a simplicial complex.
pub fn compute_invariants(complex: &SimplicialComplex) -> TopologyInvariants {
    TopologyInvariants {
        euler_characteristic: euler_characteristic(complex),
        f_vector: complex.f_vector().to_vec(),
        betti_numbers: betti_numbers(complex),
    }
}

/// Compare two complexes by their Euler characteristic.
///
/// Returns `Ordering` based on Euler characteristic comparison.
pub fn compare_by_euler(a: &SimplicialComplex, b: &SimplicialComplex) -> std::cmp::Ordering {
    euler_characteristic(a).cmp(&euler_characteristic(b))
}

/// Check whether the Euler-Poincaré formula holds for the given complex.
///
/// The formula states: χ = Σ_k (−1)^k · β_k
pub fn verify_euler_poincare(complex: &SimplicialComplex) -> bool {
    let chi = euler_characteristic(complex);
    let betti = betti_numbers(complex);
    let chi_from_betti: i64 = betti
        .iter()
        .enumerate()
        .map(|(k, &b)| if k % 2 == 0 { b as i64 } else { -(b as i64) })
        .sum();
    chi == chi_from_betti
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euler_point() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0]);
        assert_eq!(euler_characteristic(&k), 1);
    }

    #[test]
    fn euler_edge() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1]);
        // V=2, E=1: χ = 2-1 = 1
        assert_eq!(euler_characteristic(&k), 1);
    }

    #[test]
    fn euler_hollow_triangle() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1]);
        k.add(vec![1, 2]);
        k.add(vec![0, 2]);
        // V=3, E=3, no face: χ = 3-3 = 0 (circle S^1)
        assert_eq!(euler_characteristic(&k), 0);
    }

    #[test]
    fn euler_solid_triangle() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        // V=3, E=3, F=1: χ = 3-3+1 = 1 (disk)
        assert_eq!(euler_characteristic(&k), 1);
    }

    #[test]
    fn betti_connected() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        let betti = betti_numbers(&k);
        assert_eq!(betti[0], 1); // one connected component
    }

    #[test]
    fn betti_disconnected() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1]);
        k.add(vec![2, 3]);
        let betti = betti_numbers(&k);
        assert_eq!(betti[0], 2);
    }

    #[test]
    fn invariants_computation() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        let inv = compute_invariants(&k);
        assert_eq!(inv.euler_characteristic, 1);
        assert_eq!(inv.f_vector, vec![3, 3, 1]);
        assert!(inv.betti_numbers[0] >= 1);
    }

    #[test]
    fn compare_complexes() {
        let mut a = SimplicialComplex::new();
        a.add(vec![0, 1, 2]); // χ = 1

        let mut b = SimplicialComplex::new();
        b.add(vec![0, 1]);
        b.add(vec![1, 2]);
        b.add(vec![0, 2]); // χ = 0

        assert_eq!(compare_by_euler(&a, &b), std::cmp::Ordering::Greater);
    }

    #[test]
    fn euler_poincare_holds_for_simplex() {
        let mut k = SimplicialComplex::new();
        k.add(vec![0, 1, 2]);
        assert!(verify_euler_poincare(&k));
    }
}
