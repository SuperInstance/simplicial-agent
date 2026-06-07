//! # simplicial-agent
//!
//! Simplicial complexes for modeling agent collaboration topology.
//!
//! This crate provides concrete, zero-dependency (except `serde`) data structures
//! and algorithms for building, analyzing, and collapsing simplicial complexes —
//! with a focus on modeling multi-agent collaboration patterns.
//!
//! ## Modules
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`simplex`] | Individual simplex operations |
//! | [`complex`] | Simplicial complex construction and queries |
//! | [`rips`] | Vietoris-Rips complex from distance matrices |
//! | [`euler`] | Euler characteristic and Betti numbers |
//! | [`collapse`] | Elementary and strong collapse |
//! | [`collaborate`] | Agent collaboration topology analysis |

pub mod collaborate;
pub mod collapse;
pub mod complex;
pub mod euler;
pub mod rips;
pub mod simplex;
