# simplicial-agent

**Simplicial complexes for modeling agent collaboration topology.**

[![crates.io](https://img.shields.io/crates/v/simplicial-agent.svg)](https://crates.io/crates/simplicial-agent)
[![docs.rs](https://docs.rs/simplicial-agent/badge.svg)](https://docs.rs/simplicial-agent)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`simplicial-agent` provides concrete, zero-overhead data structures and algorithms for building, analyzing, and collapsing simplicial complexes — with a focus on modeling multi-agent collaboration patterns as topological objects.

---

## Table of Contents

- [Why Simplicial Topology for Agents?](#why-simplicial-topology-for-agents)
- [Theory](#theory)
- [Architecture](#architecture)
- [Quick Start](#quick-start)
- [Module Reference](#module-reference)
- [Code Examples](#code-examples)
- [Performance](#performance)
- [Design Decisions](#design-decisions)
- [Comparisons](#comparisons)
- [Practical Applications](#practical-applications)
- [API Reference](#api-reference)
- [References](#references)

---

## Why Simplicial Topology for Agents?

Graphs model pairwise relationships. But multi-agent collaboration is **multiway** — three agents working together is not just three pairwise collaborations. Simplicial complexes capture this:

- A **0-simplex** (point) represents a single agent.
- A **1-simplex** (edge) represents a pairwise collaboration.
- A **2-simplex** (triangle) represents a three-way collaboration — not decomposable into edges.
- A **k-simplex** represents a (k+1)-agent collaboration group.

This matters because:

1. **Holes** in the complex reveal missing collaborations — groups that should work together but don't.
2. **Euler characteristic** gives a single number summarizing collaboration structure.
3. **Betti numbers** count topological features: connected components (β₀), collaboration loops (β₁), enclosed collaboration spaces (β₂).
4. **Collapse** simplifies the topology while preserving its essential shape — identifying core collaboration patterns.
5. **Vietoris-Rips complexes** let you build topology from distance metrics (collaboration distance, communication latency, knowledge divergence).

---

## Theory

### Simplicial Complexes

A **simplicial complex** K is a finite collection of sets (called *simplices*) such that:

> K ⊂ 2^V such that σ ∈ K ∧ τ ⊂ σ ⟹ τ ∈ K

That is: if a simplex is in K, all of its subsets (faces) are also in K. This *closure under the face relation* is what makes K a complex rather than just a hypergraph.

An **n-simplex** σ = {v₀, v₁, …, vₙ} is a set of n+1 vertices. Its **dimension** is n. Its **boundary** ∂σ consists of all (n−1)-dimensional faces obtained by removing one vertex:

> ∂σ = { σ \ {vᵢ} : i = 0, 1, …, n }

The **f-vector** (f₀, f₁, …, f_d) counts simplices by dimension: f_k = |{σ ∈ K : dim(σ) = k}|.

### Vietoris-Rips Complex

Given a finite metric space (X, d) and threshold ε > 0, the **Vietoris-Rips complex** is:

> VR(X, ε) = { σ ⊂ X : diam(σ) ≤ ε }

where diam(σ) = max{d(xᵢ, xⱼ) : xᵢ, xⱼ ∈ σ} is the diameter of σ.

In practice, we build VR(X, ε) by:
1. Adding all points as 0-simplices.
2. Adding edges {xᵢ, xⱼ} for all pairs with d(xᵢ, xⱼ) ≤ ε.
3. Finding all cliques in the resulting graph.

A k-clique in the 1-skeleton corresponds to a (k−1)-simplex in the Rips complex. This follows from the definition: if all pairwise distances are ≤ ε, then the diameter of the set is ≤ ε.

### Euler Characteristic

The **Euler characteristic** of a simplicial complex K is the alternating sum:

> χ(K) = Σₖ (−1)^k · f_k

By the **Euler-Poincaré theorem**, this equals:

> χ(K) = Σₖ (−1)^k · β_k

where β_k = rank(H_k(K)) is the k-th **Betti number**. This means χ is a *topological invariant* — it depends only on the shape of K, not on how K is triangulated.

Betti numbers have concrete meanings:
- β₀ = number of connected components.
- β₁ = number of "holes" (1-dimensional loops).
- β₂ = number of enclosed voids (2-dimensional cavities).

### Elementary Collapse

An **elementary collapse** removes a free face-coface pair (σ, τ) where:
- dim(σ) + 1 = dim(τ)
- σ is contained in exactly one simplex of dimension dim(σ) + 1, namely τ.

The key theorem: elementary collapse **preserves homotopy type**. A complex K is **collapsible** if it can be reduced to a single vertex through a sequence of elementary collapses. Every **cone** (a complex with an apex vertex connected to all others) is collapsible.

**Strong collapse** (Barmak & Minian, 2008) removes dominated vertices: vertex v is dominated by w if N[v] ⊆ N[w] in the 1-skeleton. Strong collapse is often more powerful than elementary collapse.

### Simplicial Laplacian

The **graph Laplacian** of the 1-skeleton is L = D − A, where D is the degree matrix and A is the adjacency matrix. Its eigenvalues encode collaboration dynamics:

- λ₁ = 0 always (for connected graphs), with multiplicity equal to β₀.
- λ₂ (algebraic connectivity / Fiedler value) measures how well-connected the graph is.
- λ₂ → 0 means the complex is nearly disconnected — a collaboration bottleneck.

---

## Architecture

```
┌─────────────┐     ┌──────────────┐     ┌──────────────────┐
│   Agents     │────▶│  Distances   │────▶│  Rips Complex    │
│  (vertices)  │     │  (metric)    │     │  (VR(X, ε))      │
└─────────────┘     └──────────────┘     └────────┬─────────┘
                                                    │
                     ┌──────────────────────────────┘
                     │
                     ▼
            ┌──────────────────┐     ┌──────────────────┐
            │   Topology       │────▶│   Insights        │
            │   (Euler, Betti) │     │   (holes, cover)  │
            └────────┬─────────┘     └──────────────────┘
                     │
                     ▼
            ┌──────────────────┐
            │   Collapse       │
            │   (simplify)     │
            └──────────────────┘

Data Flow:
  agents → distance matrix → Vietoris-Rips complex → topology analysis → insights
  collaboration records → simplicial complex → collapse → core structure
```

---

## Quick Start

```toml
[dependencies]
simplicial-agent = "0.1"
```

```rust
use simplicial_agent::complex::SimplicialComplex;

let mut k = SimplicialComplex::new();
k.add(vec![0, 1, 2]); // Add a triangle (2-simplex with all faces)
k.add(vec![1, 2, 3]); // Add another triangle sharing edge {1,2}

println!("Euler characteristic: {}", k.euler_characteristic());
println!("F-vector: {:?}", k.f_vector());
println!("Components: {}", k.connected_components().len());
```

---

## Module Reference

| Module | Key Type | Description | Example |
|--------|----------|-------------|---------|
| `simplex` | `Simplex` | Individual simplex: vertex set, faces, boundary, star, link | `Simplex::from_vertices(vec![0,1,2])` |
| `complex` | `SimplicialComplex` | Collection of simplices closed under face relation | `k.add(vec![0,1,2])` |
| `rips` | `VietorisRipsBuilder` | Build Rips complex from distance matrix | `rips_complex(&dist, 1.0)` |
| `euler` | `TopologyInvariants` | Euler characteristic, Betti numbers, f-vector | `compute_invariants(&k)` |
| `collapse` | `CollapseResult` | Elementary and strong collapse | `collapse(&k)` |
| `collaborate` | `AgentCollaboration` | Agents as vertices, collaboration groups as simplices | `collab.analyze(0.5)` |

---

## Code Examples

### Example 1: Building a Simplicial Complex and Computing Euler Characteristic

```rust
use simplicial_agent::complex::SimplicialComplex;
use simplicial_agent::euler::{euler_characteristic, betti_numbers};

// Build a hollow tetrahedron (boundary of a 3-simplex, homeomorphic to S²)
let mut k = SimplicialComplex::new();

// Add all four triangular faces of the tetrahedron
k.add(vec![0, 1, 2]); // face without vertex 3
k.add(vec![0, 1, 3]); // face without vertex 2
k.add(vec![0, 2, 3]); // face without vertex 1
k.add(vec![1, 2, 3]); // face without vertex 0

// The complex automatically includes all sub-faces (edges and vertices)
// Total: 4 vertices + 6 edges + 4 triangles = 14 simplices
assert_eq!(k.num_simplices(), 14);
assert_eq!(k.num_vertices(), 4);
assert_eq!(k.num_edges(), 6);
assert_eq!(k.num_triangles(), 4);

// Euler characteristic: V - E + F = 4 - 6 + 4 = 2
// This is the Euler characteristic of S² (the 2-sphere)!
let chi = euler_characteristic(&k);
assert_eq!(chi, 2);

// Betti numbers: β₀ = 1 (connected), β₁ = 0 (no tunnels), β₂ = 1 (one enclosed void)
let betti = betti_numbers(&k);
assert_eq!(betti[0], 1); // one connected component

// f-vector: [4, 6, 4] (4 vertices, 6 edges, 4 triangles)
assert_eq!(k.f_vector(), &[4, 6, 4]);

println!("Hollow tetrahedron: χ = {}, β = {:?}", chi, betti);
```

### Example 2: Vietoris-Rips Complex from a Collaboration Distance Matrix

```rust
use simplicial_agent::rips::{rips_complex, rips_complex_bounded};

// Distance matrix for 5 agents
// Agents 0,1,2 are close together (distance 0.3-0.5)
// Agents 3,4 are close to each other but far from 0,1,2
// Agent 2 is the bridge between the two clusters
let distances = vec![
    vec![0.0, 0.3, 0.5, 2.0, 2.5],  // agent 0
    vec![0.3, 0.0, 0.4, 1.8, 2.3],  // agent 1
    vec![0.5, 0.4, 0.0, 1.0, 1.5],  // agent 2 (bridge)
    vec![2.0, 1.8, 1.0, 0.0, 0.3],  // agent 3
    vec![2.5, 2.3, 1.5, 0.3, 0.0],  // agent 4
];

// Build Rips complex with ε = 0.5 (tight collaboration threshold)
let tight = rips_complex(&distances, 0.5);
// Edges: {0,1}, {0,2}, {1,2}, {3,4}
// Triangle: {0,1,2}
assert_eq!(tight.num_edges(), 4);
assert_eq!(tight.num_triangles(), 1);

// With ε = 1.0 (looser threshold), agent 2 bridges to cluster {3,4}
let medium = rips_complex(&distances, 1.0);
// More edges, possibly more triangles
println!("Medium ε: {} edges, {} triangles", 
    medium.num_edges(), medium.num_triangles());

// Limit to 1-skeleton (graph only, no triangles)
let graph_only = rips_complex_bounded(&distances, 0.5, 1);
assert_eq!(graph_only.num_triangles(), 0);

// Compute topology
let chi = medium.euler_characteristic();
let components = medium.connected_components();
println!("Components: {}, χ = {}", components.len(), chi);
```

### Example 3: Collapsing a Complex and Checking Homotopy Type

```rust
use simplicial_agent::complex::SimplicialComplex;
use simplicial_agent::collapse::{collapse, strong_collapse, is_cone};

// A solid triangle (2-simplex {0,1,2} with all faces) is collapsible
let mut k = SimplicialComplex::new();
k.add(vec![0, 1, 2]);

// It's also a cone (any vertex connects to all others)
assert!(is_cone(&k));

// Elementary collapse reduces it to a single vertex
let result = collapse(&k);
assert!(result.is_collapsible);
println!("Collapsed in {} steps", result.steps);

// A hollow triangle (just the boundary, no face) is NOT collapsible
let mut hollow = SimplicialComplex::new();
hollow.add(vec![0, 1]);
hollow.add(vec![1, 2]);
hollow.add(vec![0, 2]);

let hollow_result = collapse(&hollow);
// No free faces to collapse — stays as S¹ (circle)
assert!(!hollow_result.is_collapsible);
assert_eq!(hollow_result.steps, 0);

// A tetrahedron is collapsible
let mut tet = SimplicialComplex::new();
tet.add(vec![0, 1, 2, 3]);
let tet_result = collapse(&tet);
assert!(tet_result.is_collapsible);

// Strong collapse on a path graph 0-1-2-3
let mut path = SimplicialComplex::new();
path.add(vec![0, 1]);
path.add(vec![1, 2]);
path.add(vec![2, 3]);

let strong_result = strong_collapse(&path);
// Endpoints 0 and 3 are dominated, will be removed iteratively
println!("Strong collapse: {} steps", strong_result.steps);
```

### Example 4: Agent Collaboration Topology Analysis

```rust
use simplicial_agent::collaborate::{AgentCollaboration, analyze_collaboration};

// Model 6 agents with collaboration records
// Each record: (participating agents, interaction strength)
let records = vec![
    (vec![0, 1], 0.95),    // agents 0,1 collaborate strongly
    (vec![1, 2], 0.85),    // agents 1,2 collaborate
    (vec![0, 2], 0.80),    // agents 0,2 collaborate
    (vec![0, 1, 2], 0.70), // all three in {0,1,2} collaborate together
    (vec![3, 4], 0.90),    // agents 3,4 collaborate
    (vec![4, 5], 0.60),    // agents 4,5 collaborate
    (vec![2, 3], 0.40),    // weak bridge between clusters
];

// Analyze with threshold 0.5
let topo = analyze_collaboration(6, &records, 0.5);

println!("Connected components: {}", topo.connected_components);
println!("Euler characteristic: {}", topo.euler_characteristic);
println!("Coverage: {:.1}%", topo.coverage * 100.0);

// Holes: missing collaborations between agents in the same component
if !topo.holes.is_empty() {
    println!("Collaboration holes detected:");
    for (a, b) in &topo.holes {
        println!("  agents {} and {} should collaborate but don't", a, b);
    }
}

// Laplacian eigenvalues reveal collaboration dynamics
// λ₂ (algebraic connectivity) near 0 = bottleneck
let lambda2 = topo.laplacian_eigenvalues.get(1).copied().unwrap_or(0.0);
println!("Algebraic connectivity (λ₂): {:.4}", lambda2);
if lambda2 < 0.1 {
    println!("⚠ Collaboration bottleneck detected!");
}

// Adjust threshold to see how topology changes
let strict_topo = analyze_collaboration(6, &records, 0.8);
println!("Strict threshold (0.8): {} components", strict_topo.connected_components);
```

---

## Performance

### Complexity Analysis

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Add simplex (with faces) | O(2^d) | d = simplex dimension; generates all sub-faces |
| Face closure check | O(\|K\|) | Verify every simplex has all faces present |
| F-vector | O(1) amortized | Maintained incrementally |
| Euler characteristic | O(d) | d = max dimension; just f-vector alternating sum |
| Connected components | O(V + E) | Union-find / DFS on 1-skeleton |
| Rips construction | O(n^k) worst case | k = max clique size; dominates at high ε |
| Clique detection | O(n^k · k) | Iterative worklist, not recursive |
| Elementary collapse | O(\|K\|²) | Each step scans all simplices for free faces |
| Strong collapse | O(V² · \|K\|) | Each step checks all vertex pairs |
| Betti numbers | O(\|K\|) | Simplified boundary rank computation |
| Laplacian eigenvalues | O(n³) | Jacobi iteration on n×n matrix |

### Practical Scaling

For typical agent collaboration scenarios (10–100 agents):

- **Rips construction** dominates at high ε where many simplices form. Use `max_dimension` to cap.
- **Collapse** is fast for collapsible complexes (terminates at a single vertex) but may do nothing for complexes with topological holes.
- **Collaboration analysis** is bounded by the Laplacian eigenvalue computation (Jacobi, O(n³)).

### Memory

Each simplex stores its sorted vertex list. A complex with S simplices of average dimension d uses approximately O(S · d) memory. The HashMap index adds O(S) overhead for key storage.

---

## Design Decisions

### Why `Vec<usize>` for vertices?

Simplicial complexes in agent systems map to concrete agent IDs (indices). Generic vertex types add API complexity without benefit for our use case. `usize` is the natural index type: it's what you get from enumeration, and it's what you use to index into arrays of agents.

### Why no abstract algebraic traits?

This crate targets **engineers building agent systems**, not mathematicians proving theorems. Abstract traits for chain groups, boundary operators, and homology computations would add cognitive overhead without making the core use case (collaboration topology) any easier. If you need abstract algebra, this crate gives you the concrete data you need to build your own abstractions on top.

### Why iterative clique detection?

Recursive backtracking for clique detection is elegant but:
1. Stack overflows on large inputs.
2. Harder to instrument (progress tracking, early termination).
3. No performance advantage — worst case is the same O(n^k).

Our iterative worklist approach (maintain a `HashSet` of known cliques, extend each by adjacent vertices) is:
1. Bounded by available heap memory.
2. Easy to add max-dimension pruning.
3. Naturally deduplicates via the `HashSet`.

### Why iterative collapse?

Same reasoning: an iterative queue of collapse operations is:
1. Bounded (no stack overflow).
2. Auditable (every step is recorded).
3. Composable (can interleave with other operations).

### Why Jacobi eigenvalues (not power iteration)?

The Laplacian is a real symmetric matrix. Jacobi rotation is the textbook algorithm for this case — guaranteed convergence, O(n³), and produces all eigenvalues. Power iteration only finds the dominant eigenvalue and requires deflation for the rest, which is fragile for near-degenerate spectra.

---

## Comparisons

### Simplicial vs. Cubical Complexes

| Feature | Simplicial | Cubical |
|---------|-----------|---------|
| Building blocks | Simplices (triangles, tetrahedra) | Cubes (squares, cubes) |
| Natural for | Point clouds, graphs, agent groups | Image data, grids |
| Face structure | Any subset | Remove one interval |
| Implementation | `Vec<usize>` vertices | Grid coordinates |
| Best for collaboration | ✅ Multiway relationships | ❌ Grid-like structure |

**When to use simplicial**: Your data is relational (agents, interactions), not spatial (pixels, voxels).

### Rips vs. Čech Complex

| Feature | Vietoris-Rips | Čech |
|---------|--------------|------|
| Condition | diam(σ) ≤ ε | ⋂ B(xᵢ, ε/2) ≠ ∅ |
| Computation | Pairwise distances + cliques | Geometric intersection tests |
| Approximation | Over-approximates Čech | Geometrically exact |
| Efficiency | O(n^k) clique detection | Depends on intersection oracle |

The Rips complex is computationally cheaper and only requires a distance matrix. The Čech complex requires geometric data (positions, not just distances). For agent collaboration, you usually have distances, not positions, so Rips is the natural choice.

### Collapse vs. Morse Theory

| Feature | Elementary Collapse | Discrete Morse Theory |
|---------|-------------------|----------------------|
| Operation | Remove free face + coface | Pair cells (critical/unpaired) |
| Preserves | Homotopy type | Homotopy type + homology |
| Strength | Limited (not all complexes collapse) | More general (Morse matchings) |
| Implementation | Simple iteration | Combinatorial optimization |

Collapse is simpler and sufficient for identifying core collaboration structure. Morse theory is more powerful but requires solving a matching problem.

---

## Practical Applications

### Team Topology Analysis

Build a simplicial complex from team interaction records. Higher-dimensional simplices indicate genuine multi-way collaboration. Missing simplices (holes) indicate teams that should collaborate but don't.

```rust
// 5-person team, 3 strong sub-groups
let records = vec![
    (vec![0, 1, 2], 0.9),  // backend team
    (vec![2, 3, 4], 0.8),  // frontend team  
    (vec![0, 4], 0.3),     // weak cross-team link
];
let topo = analyze_collaboration(5, &records, 0.5);
// Agent 2 is the bridge (in both teams)
// Holes: {0,3}, {1,3}, {1,4} should collaborate but don't
```

### Social Network Hole Detection

In a social network, a 1-dimensional hole (β₁ > 0) indicates a cycle of relationships with a missing "shortcut" — people who could benefit from direct connection. The collaboration Laplacian's second-smallest eigenvalue quantifies this bottleneck.

### Knowledge Coverage Gaps

Model knowledge areas as vertices and expert groups as simplices. Missing simplices indicate knowledge areas that no team covers collectively. The f-vector tells you how many coverage levels exist (individual expertise, pair expertise, team expertise).

### Collaboration Dynamics

The simplicial Laplacian eigenvalues track collaboration dynamics over time:
- Decreasing λ₂ → collaboration is strengthening.
- Increasing λ₂ → fragmentation is increasing.
- Sudden changes in the spectrum → organizational restructuring.

---

## API Reference

### `simplex` — Individual Simplex Operations

```rust
// Create a simplex
let s = Simplex::from_vertices(vec![0, 1, 2]); // 2-simplex (triangle)

// Properties
s.dimension();           // 2
s.vertices();            // &[0, 1, 2]
s.len();                 // 3

// Relations
s.is_face_of(&other);    // true if every vertex of s is in other
s.is_proper_face_of(&other);  // face but not equal

// Operations
s.boundary();            // all (n-1)-faces: [{0,1}, {0,2}, {1,2}]
s.faces_of_dimension(1); // edges: [{0,1}, {0,2}, {1,2}]
s.all_faces();           // all 2^n sub-faces

// Free functions
star(&sigma, &simplices);   // all simplices containing sigma
link(&sigma, &simplices);   // faces of star that don't intersect sigma
closure(&[sigma]);          // sigma plus all its faces
```

### `complex` — Simplicial Complex

```rust
let mut k = SimplicialComplex::new();
k.add(vec![0, 1, 2]);  // adds triangle + all faces

k.num_simplices();           // 7 (3 vertices + 3 edges + 1 triangle)
k.num_vertices();            // 3
k.num_edges();               // 3
k.num_triangles();           // 1
k.max_dimension();           // 2
k.f_vector();                // &[3, 3, 1]
k.euler_characteristic();    // 1
k.is_simplicial_complex();   // true
k.contains(&Simplex::from_vertices(vec![0, 1])); // true

k.skeleton(1);               // 1-skeleton (graph)
k.connected_components();    // vec of vertex groups
```

### `rips` — Vietoris-Rips Complex

```rust
let distances = vec![
    vec![0.0, 1.0, 2.0],
    vec![1.0, 0.0, 1.0],
    vec![2.0, 1.0, 0.0],
];

// Simple interface
let complex = rips_complex(&distances, 1.5);

// With dimension bound
let bounded = rips_complex_bounded(&distances, 1.5, 2);

// Full builder API
use simplicial_agent::rips::{VietorisRipsBuilder, RipsConfig};
let config = RipsConfig::new(1.5).with_max_dimension(2);
let complex = VietorisRipsBuilder::new(config).build(&distances);
```

### `euler` — Topological Invariants

```rust
use simplicial_agent::euler::{euler_characteristic, betti_numbers, compute_invariants};

let chi = euler_characteristic(&k);
let betti = betti_numbers(&k);
let inv = compute_invariants(&k);
// inv.euler_characteristic, inv.f_vector, inv.betti_numbers

// Compare two complexes
use simplicial_agent::euler::compare_by_euler;
let ordering = compare_by_euler(&k1, &k2);

// Verify Euler-Poincaré formula
use simplicial_agent::euler::verify_euler_poincare;
assert!(verify_euler_poincare(&k));
```

### `collapse` — Complex Simplification

```rust
use simplicial_agent::collapse::{collapse, strong_collapse, is_cone};

let result = collapse(&k);
result.is_collapsible;    // reduced to single vertex?
result.steps;             // number of collapse steps
result.complex;           // the collapsed complex
result.removed_pairs;     // (face, coface) pairs removed

let strong = strong_collapse(&k);
// same fields, but using vertex domination

is_cone(&k);  // true if any vertex connects to all others
```

### `collaborate` — Agent Topology

```rust
use simplicial_agent::collaborate::{AgentCollaboration, analyze_collaboration};

let mut collab = AgentCollaboration::new(5);
collab.add_collaboration(vec![0, 1], 0.9);
collab.add_collaboration(vec![1, 2], 0.8);
collab.add_collaboration(vec![0, 1, 2], 0.7);

let topo = collab.analyze(0.5);
topo.connected_components;      // number of groups
topo.euler_characteristic;      // topology summary
topo.betti_numbers;             // [components, holes, voids, ...]
topo.coverage;                  // fraction of pairs collaborating
topo.holes;                     // missing collaborations
topo.laplacian_eigenvalues;     // dynamics spectrum

// Convenience function
let topo = analyze_collaboration(5, &[
    (vec![0, 1], 0.9),
    (vec![1, 2], 0.8),
], 0.5);
```

---

## Crate Features

- **`serde`** (default): `Serialize` and `Deserialize` for all public types.
- No other dependencies. Zero external crates beyond `serde`.

---

## Minimum Supported Rust Version (MSRV)

Edition 2024 requires Rust 1.85+.

---

## License

MIT License. See [LICENSE](LICENSE) or <https://opensource.org/licenses/MIT>.

---

## References

1. **Hatcher, A.** (2002). *Algebraic Topology*. Cambridge University Press.  
   The standard reference for simplicial complexes, homology, and the Euler-Poincaré theorem. Available free at <https://pi.math.cornell.edu/~hatcher/AT/AT.pdf>

2. **Carlsson, G.** (2009). "Topology and Data." *Bulletin of the American Mathematical Society*, 46(2), 255–308.  
   Foundational paper on applying topological data analysis to point clouds via Vietoris-Rips and Čech complexes.

3. **Ghrist, R.** (2014). *Elementary Applied Topology*. CreateSpace.  
   Accessible treatment of applied topology including simplicial complexes, sheaves, and persistence. Chapter 1 covers complexes comprehensively.

4. **Zomorodian, A.** (2005). *Computational Topology*. In *Proceedings of the 21st Annual Symposium on Computational Geometry*.  
   Algorithms for building and simplifying simplicial complexes from data.

5. **Edelsbrunner, H. & Harer, J.** (2010). *Computational Topology: An Introduction*. American Mathematical Society.  
   Comprehensive textbook covering simplicial complexes, persistence, and Morse theory. Chapter 3 covers Rips complexes.

6. **Barmak, J.A. & Minian, E.G.** (2008). "Strong homotopy types, nerve complexes and Dunwoody's accessibility theorem." *Revista de la Unión Matemática Argentina*, 49(1), 1–19.  
   Introduces strong collapse for simplicial complexes, proving it preserves simple homotopy type.

7. **Patania, A., Petri, G. & Vaccarino, F.** (2017). "The shape of collaborations." *EPJ Data Science*, 6(1), 1–16.  
   Applies simplicial complexes to model multi-agent collaboration topology, introducing higher-order network analysis.

8. **Schaub, M.T., Benson, A.R., Horn, P., Lippner, G. & Jadbabaie, A.** (2020). "Random walks on simplicial complexes and the normalized Hodge 1-Laplacian." *SIAM Review*, 62(2), 353–391.  
   Develops the simplicial Laplacian for higher-order network analysis, including diffusion on collaboration complexes.

9. **Chung, F.** (1997). *Spectral Graph Theory*. CBMS Regional Conference Series in Mathematics, No. 92. American Mathematical Society.  
   Reference for graph Laplacian eigenvalues and their relationship to connectivity (algebraic connectivity / Fiedler value).
