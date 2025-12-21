//! Gaussian-Type Orbital (GTO) approximation for screened Coulomb interactions.
//!
//! This module implements the calculation of the effective electrostatic interaction
//! between two charge distributions approximated by Gaussian functions. This approach
//! allows for efficient analytical evaluation of the Coulomb integral, serving as a
//! faster but approximate alternative to the exact Slater-Type Orbital (STO) method.

use super::constants::DISTANCE_THRESHOLD_BOHR;
use libm::erf;
use std::f64::consts::PI;

/// Computes the screened Coulomb integral between two atoms using
/// the Gaussian-Type Orbital (GTO) approximation.
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

    let alpha1 = compute_gaussian_exponent(n1, zeta1);
    let alpha2 = compute_gaussian_exponent(n2, zeta2);

    compute_screened_potential(distance_bohr, alpha1, alpha2)
}

/// Computes the Slater exponent from atomic parameters.
#[inline]
fn compute_slater_exponent(n: u8, r_cov_bohr: f64, lambda: f64) -> f64 {
    if r_cov_bohr < DISTANCE_THRESHOLD_BOHR {
        return 0.0;
    }
    lambda * (2.0 * n as f64 + 1.0) / (2.0 * r_cov_bohr)
}

/// Computes the equivalent Gaussian exponent alpha from the Slater exponent zeta.
#[inline]
fn compute_gaussian_exponent(n: u8, zeta: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let c_n = get_fitting_coefficient(n);
    c_n * zeta * zeta / (n as f64)
}

/// Returns the empirical scaling coefficient for the GTO approximation.
///
/// These values are derived to minimize the difference between the Slater and
/// Gaussian density profiles for a given principal quantum number n.
#[inline]
fn get_fitting_coefficient(n: u8) -> f64 {
    match n {
        1 => 0.270917,
        2 => 0.098800,
        3 => 0.055600,
        4 => 0.039100,
        5 => 0.029600,
        _ => 0.27 * (n as f64).powf(-1.35),
    }
}

/// Evaluates the Coulomb integral between two Gaussian charge distributions.
#[inline]
fn compute_screened_potential(r: f64, alpha1: f64, alpha2: f64) -> f64 {
    let sum_alpha = alpha1 + alpha2;
    if sum_alpha < 1e-10 {
        return 0.0;
    }

    let alpha_eff = (alpha1 * alpha2) / sum_alpha;

    if r > DISTANCE_THRESHOLD_BOHR {
        let sqrt_alpha_eff = alpha_eff.sqrt();
        erf(sqrt_alpha_eff * r) / r
    } else {
        2.0 * (alpha_eff / PI).sqrt()
    }
}
