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
