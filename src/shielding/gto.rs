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
    let beta = (2.0 * alpha1 * alpha2 / (alpha1 + alpha2)).sqrt();

    if r > DISTANCE_THRESHOLD_BOHR {
        erf(beta * r) / r
    } else {
        2.0 * beta / PI.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shielding::constants;
    use approx::assert_relative_eq;

    #[test]
    fn test_slater_exponent_calculation() {
        let r_cov_bohr_h = 0.371 / constants::BOHR_TO_ANGSTROM;
        let result = compute_slater_exponent(1, r_cov_bohr_h, 0.5);
        let expected = (2.0 * 1.0 + 1.0) / (4.0 * r_cov_bohr_h);
        assert_relative_eq!(result, expected, epsilon = 1e-12);

        let r_cov_bohr_c = 0.759 / constants::BOHR_TO_ANGSTROM;
        let result_c = compute_slater_exponent(2, r_cov_bohr_c, 0.5);
        let expected_c = (2.0 * 2.0 + 1.0) / (4.0 * r_cov_bohr_c);
        assert_relative_eq!(result_c, expected_c, epsilon = 1e-12);

        let result_lambda_custom = compute_slater_exponent(1, r_cov_bohr_h, 0.4913);
        let expected_lambda_custom = 0.4913 * (2.0 * 1.0 + 1.0) / (2.0 * r_cov_bohr_h);
        assert_relative_eq!(
            result_lambda_custom,
            expected_lambda_custom,
            epsilon = 1e-12
        );

        assert_eq!(compute_slater_exponent(1, 0.0, 0.5), 0.0);
    }

    #[test]
    fn test_get_fitting_coefficient_values() {
        assert_eq!(get_fitting_coefficient(1), 0.270917);
        assert_eq!(get_fitting_coefficient(2), 0.098800);
        assert_eq!(get_fitting_coefficient(3), 0.055600);
        assert_eq!(get_fitting_coefficient(4), 0.039100);
        assert_eq!(get_fitting_coefficient(5), 0.029600);

        let n6 = 6u8;
        let expected_c6 = 0.27 * (n6 as f64).powf(-1.35);
        assert_eq!(get_fitting_coefficient(n6), expected_c6);
    }

    #[test]
    fn test_gaussian_exponent_calculation() {
        assert_relative_eq!(
            compute_gaussian_exponent(1, 1.0),
            get_fitting_coefficient(1) / 1.0
        );
        assert_relative_eq!(
            compute_gaussian_exponent(2, 1.0),
            get_fitting_coefficient(2) / 2.0
        );

        let zeta = 2.0;
        let alpha = compute_gaussian_exponent(1, zeta);
        let expected_alpha = get_fitting_coefficient(1) * zeta.powi(2) / 1.0;
        assert_relative_eq!(alpha, expected_alpha);

        assert_eq!(compute_gaussian_exponent(0, 1.0), 0.0);
    }

    #[test]
    fn test_integral_at_zero_distance_limit() {
        let r_cov_bohr = 0.371 / constants::BOHR_TO_ANGSTROM;
        let result_at_zero = calculate_integral(0.0, 1, r_cov_bohr, 1, r_cov_bohr, 0.5);
        let result_near_zero = calculate_integral(1e-15, 1, r_cov_bohr, 1, r_cov_bohr, 0.5);

        assert!(result_at_zero > 0.0, "Result at R=0 must be positive");
        assert!(
            !result_at_zero.is_infinite(),
            "Result at R=0 must be finite"
        );

        assert_relative_eq!(result_at_zero, result_near_zero, epsilon = 1e-9);
    }

    #[test]
    fn test_integral_at_large_distance_asymptote() {
        let r_cov_bohr_h = 0.371 / constants::BOHR_TO_ANGSTROM;
        let r_cov_bohr_c = 0.759 / constants::BOHR_TO_ANGSTROM;
        let large_distance_bohr = 1000.0;

        let result = calculate_integral(large_distance_bohr, 1, r_cov_bohr_h, 2, r_cov_bohr_c, 0.5);

        let point_charge_result = 1.0 / large_distance_bohr;

        assert_relative_eq!(result, point_charge_result, epsilon = 1e-9);
    }

    #[test]
    fn test_integral_with_known_values() {
        const LAMBDA: f64 = 0.4913;

        let r_h2_bohr = 0.74 / constants::BOHR_TO_ANGSTROM;
        let r_cov_h_bohr = 0.371 / constants::BOHR_TO_ANGSTROM;
        let j_ev_h2 = calculate_integral(r_h2_bohr, 1, r_cov_h_bohr, 1, r_cov_h_bohr, LAMBDA)
            * constants::HARTREE_TO_EV;
        assert_relative_eq!(j_ev_h2, 14.02504752, epsilon = 1e-8);

        let r_co_bohr = 1.128 / constants::BOHR_TO_ANGSTROM;
        let r_cov_c_bohr = 0.759 / constants::BOHR_TO_ANGSTROM;
        let r_cov_o_bohr = 0.669 / constants::BOHR_TO_ANGSTROM;
        let j_ev_co = calculate_integral(r_co_bohr, 2, r_cov_c_bohr, 2, r_cov_o_bohr, LAMBDA)
            * constants::HARTREE_TO_EV;
        assert_relative_eq!(j_ev_co, 5.837573337, epsilon = 1e-8);
    }

    #[test]
    fn test_integral_symmetry_property() {
        let r_cov_bohr_h = 0.371 / constants::BOHR_TO_ANGSTROM;
        let r_cov_bohr_si = 1.176 / constants::BOHR_TO_ANGSTROM;
        let distance_bohr = 3.0;
        let lambda = 0.5;

        let result_h_si =
            calculate_integral(distance_bohr, 1, r_cov_bohr_h, 3, r_cov_bohr_si, lambda);
        let result_si_h =
            calculate_integral(distance_bohr, 3, r_cov_bohr_si, 1, r_cov_bohr_h, lambda);

        assert_eq!(result_h_si, result_si_h);
    }

    #[test]
    fn test_integral_scaling_with_lambda() {
        let r_cov_bohr = 0.7 / constants::BOHR_TO_ANGSTROM;
        let distance_bohr = 1.5;

        let result_lambda_low =
            calculate_integral(distance_bohr, 1, r_cov_bohr, 1, r_cov_bohr, 0.4);
        let result_lambda_high =
            calculate_integral(distance_bohr, 1, r_cov_bohr, 1, r_cov_bohr, 0.6);

        assert!(
            result_lambda_low < result_lambda_high,
            "Smaller lambda should lead to stronger screening (smaller J_AB)"
        );
    }

    #[test]
    fn test_integral_chemical_trends() {
        let r_cov_bohr = 0.7 / constants::BOHR_TO_ANGSTROM;
        let j_at_1_angstrom = calculate_integral(
            1.0 / constants::BOHR_TO_ANGSTROM,
            1,
            r_cov_bohr,
            1,
            r_cov_bohr,
            0.5,
        );
        let j_at_3_angstroms = calculate_integral(
            3.0 / constants::BOHR_TO_ANGSTROM,
            1,
            r_cov_bohr,
            1,
            r_cov_bohr,
            0.5,
        );
        assert!(
            j_at_1_angstrom > j_at_3_angstroms,
            "Interaction should decrease with distance"
        );

        let distance_bohr = 3.0 / constants::BOHR_TO_ANGSTROM;

        let r_h_bohr = 0.371 / constants::BOHR_TO_ANGSTROM;
        let j_hh = calculate_integral(distance_bohr, 1, r_h_bohr, 1, r_h_bohr, 0.5);

        let r_na_bohr = 2.085 / constants::BOHR_TO_ANGSTROM;
        let j_nana = calculate_integral(distance_bohr, 3, r_na_bohr, 3, r_na_bohr, 0.5);

        assert!(
            j_hh > j_nana,
            "Interaction between diffuse orbitals (Na-Na) should be weaker than between compact ones (H-H) at the same distance"
        );
    }
}
