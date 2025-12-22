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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shielding::constants;
    use approx::assert_relative_eq;

    #[test]
    fn test_sto_integral_symmetry() {
        let r_cov_bohr_h = 0.371 / constants::BOHR_TO_ANGSTROM;
        let r_cov_bohr_c = 0.759 / constants::BOHR_TO_ANGSTROM;
        let distance_bohr = 2.0;
        let lambda = 0.5;

        let result_h_c =
            calculate_integral(distance_bohr, 1, r_cov_bohr_h, 2, r_cov_bohr_c, lambda);
        let result_c_h =
            calculate_integral(distance_bohr, 2, r_cov_bohr_c, 1, r_cov_bohr_h, lambda);

        assert_relative_eq!(result_h_c, result_c_h, epsilon = 1e-12);
    }

    #[test]
    fn test_sto_integral_large_distance() {
        let r_cov_bohr = 1.0;
        let distance_bohr = 100.0;
        let lambda = 0.5;

        let result = calculate_integral(distance_bohr, 1, r_cov_bohr, 1, r_cov_bohr, lambda);
        let expected = 1.0 / distance_bohr;

        assert_relative_eq!(result, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_sto_integral_zero_distance() {
        let r_cov_bohr = 1.0;
        let lambda = 0.5;

        let result_zero = calculate_integral(0.0, 1, r_cov_bohr, 1, r_cov_bohr, lambda);
        let result_near_zero = calculate_integral(1e-10, 1, r_cov_bohr, 1, r_cov_bohr, lambda);

        assert!(result_zero.is_finite());
        assert!(result_zero > 0.0);
        assert_relative_eq!(result_zero, result_near_zero, epsilon = 1e-6);
    }

    #[test]
    fn test_sto_vs_gto_consistency() {
        use crate::shielding::gto;

        let r_cov_bohr = 0.7;
        let distance_bohr = 1.5;
        let lambda = 0.5;
        let n = 2;

        let sto_val = calculate_integral(distance_bohr, n, r_cov_bohr, n, r_cov_bohr, lambda);
        let gto_val = gto::calculate_integral(distance_bohr, n, r_cov_bohr, n, r_cov_bohr, lambda);

        assert_relative_eq!(sto_val, gto_val, epsilon = 0.25);
    }
}
