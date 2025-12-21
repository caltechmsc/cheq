//! Exact two-center Coulomb integrals over ns Slater-type orbitals.
//!
//! This module implements the EXACT computation of the Coulomb integral J_AB(R)
//! between two ns Slater orbital charge distributions, as required by the
//! Charge Equilibration (QEq) method of Rappé & Goddard (1991).
//!
//! This implementation delegates to the high-performance `sto-ns` crate.

use super::constants::DISTANCE_THRESHOLD_BOHR;

/// Computes the two-center Coulomb integral between two ns Slater orbital
/// charge distributions using the exact expansion.
///
/// # Arguments
///
/// * `distance_bohr` - The distance between the two atomic centers in Bohr.
/// * `n1` - Principal quantum number of the first atom.
/// * `r_cov1_bohr` - Covalent radius of the first atom in Bohr.
/// * `n2` - Principal quantum number of the second atom.
/// * `r_cov2_bohr` - Covalent radius of the second atom in Bohr.
/// * `lambda` - The global screening scaling parameter.
///
/// # Returns
///
/// The Coulomb integral value in Hartree.
#[inline]
pub fn calculate_integral(
    distance_bohr: f64,
    n1: u8,
    r_cov1_bohr: f64,
    n2: u8,
    r_cov2_bohr: f64,
    lambda: f64,
) -> f64 {
    let zeta1 = compute_slater_exponent(n1, r_cov1_bohr, lambda);
    let zeta2 = compute_slater_exponent(n2, r_cov2_bohr, lambda);

    sto_ns::sto_coulomb_integral(distance_bohr, n1, zeta1, n2, zeta2)
}

/// Computes the Slater exponent from atomic parameters.
#[inline]
fn compute_slater_exponent(n: u8, r_cov_bohr: f64, lambda: f64) -> f64 {
    if r_cov_bohr < DISTANCE_THRESHOLD_BOHR {
        return 0.0;
    }
    lambda * (2.0 * n as f64 + 1.0) / (2.0 * r_cov_bohr)
}
