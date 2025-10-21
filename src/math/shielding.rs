use super::constants::DISTANCE_THRESHOLD_BOHR;
use libm::erf;

#[inline]
pub fn gaussian_coulomb_integral(
    distance_bohr: f64,
    n1: u8,
    r_cov1_bohr: f64,
    n2: u8,
    r_cov2_bohr: f64,
    lambda: f64,
) -> f64 {
    let zeta1 = slater_exponent(n1, r_cov1_bohr, lambda);
    let zeta2 = slater_exponent(n2, r_cov2_bohr, lambda);

    let alpha1 = equivalent_gaussian_exponent(n1, zeta1);
    let alpha2 = equivalent_gaussian_exponent(n2, zeta2);

    screened_potential(distance_bohr, alpha1, alpha2)
}

#[inline]
fn slater_exponent(n: u8, r_cov_bohr: f64, lambda: f64) -> f64 {
    if r_cov_bohr < 1e-8 {
        return 0.0;
    }
    lambda * (2.0 * n as f64 + 1.0) / (2.0 * r_cov_bohr)
}

#[inline]
fn equivalent_gaussian_exponent(n: u8, zeta: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let c_n = get_c_n_fit_coeff(n);
    c_n * zeta.powi(2) / (n as f64)
}
