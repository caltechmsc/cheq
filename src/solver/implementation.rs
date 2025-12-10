//! This module implements the core `QEqSolver` for performing charge equilibration calculations.
//!
//! The `QEqSolver` encapsulates the Self-Consistent Field (SCF) iterative procedure that solves
//! the QEq equations to determine partial atomic charges. It constructs a linear system based on
//! atomic parameters and geometry, iteratively refining charges until convergence. The solver
//! integrates with the broader `cheq` architecture by using the `AtomView` trait for atom data
//! and `Parameters` for element-specific values, enabling decoupled and flexible molecular
//! simulations.

use super::options::SolverOptions;
use crate::{
    error::CheqError,
    math,
    params::{ElementData, Parameters},
    types::{AtomView, CalculationResult},
};
use faer::{Col, Mat, prelude::*};
use rayon::prelude::*;
use std::panic::{self, AssertUnwindSafe};

/// Charge dependence factor for hydrogen hardness, derived from empirical fitting.
/// This constant adjusts the hardness of hydrogen atoms based on their partial charge,
/// introducing non-linearity into the QEq equations.
const H_CHARGE_DEPENDENCE_FACTOR: f64 = 0.93475415965;

/// The main solver for charge equilibration calculations.
///
/// This struct holds references to atomic parameters and solver options, providing methods
/// to perform QEq calculations on molecular systems. It implements an iterative SCF procedure
/// to solve the non-linear charge equilibration equations.
pub struct QEqSolver<'p> {
    /// Reference to the atomic parameters used in calculations.
    parameters: &'p Parameters,
    /// Configuration options for the solver, such as convergence tolerance and iteration limits.
    options: SolverOptions,
}

impl<'p> QEqSolver<'p> {
    /// Creates a new `QEqSolver` with default options.
    ///
    /// # Arguments
    ///
    /// * `parameters` - A reference to the `Parameters` containing element data.
    ///
    /// # Returns
    ///
    /// A new `QEqSolver` instance with default `SolverOptions`.
    ///
    /// # Examples
    ///
    /// ```
    /// use cheq::{get_default_parameters, QEqSolver};
    ///
    /// let params = get_default_parameters();
    /// let solver = QEqSolver::new(&params);
    /// ```
    pub fn new(parameters: &'p Parameters) -> Self {
        Self {
            parameters,
            options: SolverOptions::default(),
        }
    }

    /// Configures the solver with custom options.
    ///
    /// This method allows setting non-default solver parameters such as tolerance and maximum
    /// iterations. It consumes the solver and returns a new instance with the updated options.
    ///
    /// # Arguments
    ///
    /// * `options` - The `SolverOptions` to apply to the solver.
    ///
    /// # Returns
    ///
    /// A new `QEqSolver` instance with the specified options.
    ///
    /// # Examples
    ///
    /// ```
    /// use cheq::{get_default_parameters, QEqSolver, SolverOptions};
    ///
    /// let params = get_default_parameters();
    /// let options = SolverOptions { tolerance: 1e-8, max_iterations: 50, lambda_scale: 1.0 };
    /// let solver = QEqSolver::new(&params).with_options(options);
    /// ```
    pub fn with_options(mut self, options: SolverOptions) -> Self {
        self.options = options;
        self
    }

    /// Solves the charge equilibration equations for a given molecular system.
    ///
    /// This method performs the SCF iterative procedure to compute partial atomic charges that
    /// equalize the chemical potential across all atoms, subject to the total charge constraint.
    /// The process involves building and solving a linear system in each iteration, with special
    /// handling for hydrogen atoms whose hardness depends on their charge.
    ///
    /// # Arguments
    ///
    /// * `atoms` - A slice of atom data implementing the `AtomView` trait.
    /// * `total_charge` - The desired total charge of the system.
    ///
    /// # Returns
    ///
    /// A `Result` containing `CalculationResult` with the computed charges and metadata on success,
    /// or a `CheqError` on failure.
    ///
    /// # Errors
    ///
    /// This function can return the following errors:
    /// * `CheqError::NoAtoms` if the atom slice is empty.
    /// * `CheqError::ParameterNotFound` if parameters for an atom's element are missing.
    /// * `CheqError::LinalgError` if the linear system cannot be solved due to numerical issues.
    /// * `CheqError::NotConverged` if the SCF procedure does not converge within the maximum iterations.
    ///
    /// # Examples
    ///
    /// ```
    /// use cheq::{get_default_parameters, QEqSolver, Atom};
    ///
    /// let params = get_default_parameters();
    /// let solver = QEqSolver::new(&params);
    /// let atoms = vec![
    ///     Atom { atomic_number: 6, position: [0.0, 0.0, 0.0] },
    ///     Atom { atomic_number: 8, position: [1.128, 0.0, 0.0] },
    /// ];
    /// let result = solver.solve(&atoms, 0.0).unwrap();
    /// assert_eq!(result.charges.len(), 2);
    /// ```
    pub fn solve<A: AtomView>(
        &self,
        atoms: &[A],
        total_charge: f64,
    ) -> Result<CalculationResult, CheqError> {
        let n_atoms = atoms.len();
        if n_atoms == 0 {
            return Err(CheqError::NoAtoms);
        }

        let element_data = self.fetch_element_data(atoms)?;
        let invariant = self.build_invariant_system(atoms, &element_data, total_charge)?;
        let has_hydrogen = !invariant.hydrogen_meta.is_empty();

        let mut work_matrix = invariant.base_matrix.clone();
        let mut charges = Col::zeros(n_atoms);
        let mut max_charge_delta = 0.0;

        for iteration in 1..=self.options.max_iterations {
            if has_hydrogen {
                for &(idx, hardness) in &invariant.hydrogen_meta {
                    work_matrix[(idx, idx)] =
                        hardness * (1.0 + charges[idx] * H_CHARGE_DEPENDENCE_FACTOR);
                }
            }

            let solve_result = panic::catch_unwind(AssertUnwindSafe(|| {
                work_matrix.partial_piv_lu().solve(&invariant.rhs)
            }));

            let solution = match solve_result {
                Ok(sol) => sol,
                Err(_) => {
                    return Err(CheqError::LinalgError(
                        "Linear system is likely singular or ill-conditioned due to input geometry.".to_string(),
                    ));
                }
            };

            let new_charges = solution.as_ref().subrows(0, n_atoms);
            max_charge_delta = new_charges
                .as_ref()
                .iter()
                .zip(charges.as_ref().iter())
                .map(|(new, old): (&f64, &f64)| (*new - *old).abs())
                .fold(0.0, f64::max);
            charges.as_mut().copy_from(&new_charges);

            if !has_hydrogen || max_charge_delta < self.options.tolerance {
                let equilibrated_potential = -solution[n_atoms];
                return Ok(CalculationResult {
                    charges: charges.as_ref().iter().cloned().collect(),
                    equilibrated_potential,
                    iterations: iteration,
                });
            }
        }

        Err(CheqError::NotConverged {
            max_iterations: self.options.max_iterations,
            delta: max_charge_delta,
        })
    }

    /// Retrieves element data for each atom from the parameters.
    ///
    /// This helper method maps each atom to its corresponding `ElementData` based on atomic number,
    /// ensuring all required parameters are available before proceeding with calculations.
    ///
    /// # Arguments
    ///
    /// * `atoms` - A slice of atom data implementing the `AtomView` trait.
    ///
    /// # Returns
    ///
    /// A vector of references to `ElementData` for each atom.
    ///
    /// # Errors
    ///
    /// Returns `CheqError::ParameterNotFound` if data for any atomic number is missing.
    fn fetch_element_data<A: AtomView>(
        &self,
        atoms: &[A],
    ) -> Result<Vec<&'p ElementData>, CheqError> {
        atoms
            .iter()
            .map(|atom| {
                let atomic_number = atom.atomic_number();
                self.parameters
                    .elements
                    .get(&atomic_number)
                    .ok_or(CheqError::ParameterNotFound(atomic_number))
            })
            .collect()
    }

    /// Builds the linear system matrix and vector for the QEq equations.
    ///
    /// This method constructs the coefficient matrix `A` and right-hand side vector `b` for the
    /// linear system `A x = b`, where `x` includes atomic charges and the Lagrange multiplier
    /// for the total charge constraint. Diagonal elements represent atomic hardness, off-diagonal
    /// elements are screened Coulomb integrals, and the last row enforces charge conservation.
    ///
    /// # Arguments
    ///
    /// * `atoms` - A slice of atom data implementing the `AtomView` trait.
    /// * `element_data` - Pre-fetched element data for each atom.
    /// * `total_charge` - The desired total charge of the system.
    /// * `current_charges` - Current estimates of atomic charges for hydrogen hardness adjustment.
    ///
    /// # Returns
    ///
    /// A tuple `(A, b)` representing the matrix and vector of the linear system.
    ///
    /// # Errors
    ///
    /// This function does not directly return errors but propagates any from underlying calculations.
    fn build_system<A: AtomView>(
        &self,
        atoms: &[A],
        element_data: &[&'p ElementData],
        total_charge: f64,
        current_charges: ColRef<f64>,
    ) -> Result<(Mat<f64>, Col<f64>), CheqError> {
        let n_atoms = atoms.len();
        let matrix_size = n_atoms + 1;

        let mut a = Mat::zeros(matrix_size, matrix_size);
        let mut b = Col::zeros(matrix_size);
        let mut a_mut = a.as_mut();
        let mut b_mut = b.as_mut();

        for i in 0..n_atoms {
            let data_i = element_data[i];

            // Adjust hardness for hydrogen based on current charge to model non-linearity
            let j_ii = if data_i.principal_quantum_number == 1 {
                data_i.hardness * (1.0 + current_charges[i] * H_CHARGE_DEPENDENCE_FACTOR)
            } else {
                data_i.hardness
            };
            a_mut[(i, i)] = j_ii;

            for j in (i + 1)..n_atoms {
                let data_j = element_data[j];

                let pos_i = atoms[i].position();
                let pos_j = atoms[j].position();
                let dist_sq: f64 = pos_i
                    .iter()
                    .zip(pos_j.iter())
                    .map(|(pi, pj)| (pi - pj).powi(2))
                    .sum();
                let distance_angstrom = dist_sq.sqrt();

                // Convert to atomic units for shielding calculations
                let distance_bohr = distance_angstrom / math::constants::BOHR_TO_ANGSTROM;
                let radius_i_bohr = data_i.radius / math::constants::BOHR_TO_ANGSTROM;
                let radius_j_bohr = data_j.radius / math::constants::BOHR_TO_ANGSTROM;

                let j_ij_hartree = math::shielding::gaussian_coulomb_integral(
                    distance_bohr,
                    data_i.principal_quantum_number,
                    radius_i_bohr,
                    data_j.principal_quantum_number,
                    radius_j_bohr,
                    self.options.lambda_scale,
                );

                // Convert back to eV for consistency with electronegativity units
                let j_ij_ev = j_ij_hartree * math::constants::HARTREE_TO_EV;
                a_mut[(i, j)] = j_ij_ev;
                a_mut[(j, i)] = j_ij_ev;
            }

            b_mut[i] = -data_i.electronegativity;
        }

        // Apply total charge constraint to the last row of the matrix
        a.col_mut(matrix_size - 1)
            .subrows_mut(0, n_atoms)
            .fill(-1.0);
        a.row_mut(matrix_size - 1).subcols_mut(0, n_atoms).fill(1.0);

        b_mut[matrix_size - 1] = total_charge;

        Ok((a, b))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::get_default_parameters;
    use crate::types::Atom;
    use approx::assert_relative_eq;
    use std::panic;

    fn get_test_solver() -> QEqSolver<'static> {
        let params = get_default_parameters();
        QEqSolver::new(params)
    }

    #[test]
    fn test_linear_system_carbon_monoxide() {
        let solver = get_test_solver();
        let atoms = vec![
            Atom {
                atomic_number: 6,
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 8,
                position: [1.128, 0.0, 0.0],
            },
        ];

        let result = solver.solve(&atoms, 0.0).unwrap();
        assert_eq!(result.iterations, 1);
        let (q_c, q_o) = (result.charges[0], result.charges[1]);
        assert!(q_o < 0.0);
        assert!(q_c > 0.0);
        assert_relative_eq!(q_c + q_o, 0.0, epsilon = 1e-9);
        assert!(
            q_o < -0.2 && q_o > -0.6,
            "Oxygen charge magnitude seems off. Got: {}",
            q_o
        );
    }

    #[test]
    fn test_nonlinear_system_water_molecule() {
        let solver = get_test_solver();
        let bond_angle_rad = 104.45f64.to_radians();
        let bond_length = 0.9575;
        let atoms = vec![
            Atom {
                atomic_number: 8,
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                position: [bond_length, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                position: [
                    bond_length * bond_angle_rad.cos(),
                    bond_length * bond_angle_rad.sin(),
                    0.0,
                ],
            },
        ];

        let result = solver.solve(&atoms, 0.0).unwrap();
        assert!(result.iterations > 1);
        let (q_o, q_h1, q_h2) = (result.charges[0], result.charges[1], result.charges[2]);

        assert!(q_o < 0.0);
        assert!(q_h1 > 0.0);
        assert_relative_eq!(q_h1, q_h2, epsilon = 1e-7);
        assert_relative_eq!(q_o + q_h1 + q_h2, 0.0, epsilon = 1e-9);

        assert!(
            q_h1 > 0.15 && q_h1 < 0.4,
            "Hydrogen charge is outside a reasonable physical range for water"
        );
    }

    #[test]
    fn test_symmetry_dihydrogen_molecule() {
        let solver = get_test_solver();
        let atoms = vec![
            Atom {
                atomic_number: 1,
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                position: [0.74, 0.0, 0.0],
            },
        ];

        let result = solver.solve(&atoms, 0.0).unwrap();
        assert_relative_eq!(result.charges[0], 0.0, epsilon = 1e-9);
        assert_relative_eq!(result.charges[1], 0.0, epsilon = 1e-9);
    }

    #[test]
    fn test_ionic_system_lithium_hydride() {
        let solver = get_test_solver();
        let atoms = vec![
            Atom {
                atomic_number: 3,
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                position: [1.595, 0.0, 0.0],
            },
        ];

        let result = solver.solve(&atoms, 0.0).unwrap();
        assert!(result.iterations > 1);
        let (q_li, q_h) = (result.charges[0], result.charges[1]);

        assert!(q_h < 0.0, "Hydrogen should be negative in LiH");
        assert!(q_li > 0.0, "Lithium should be positive in LiH");
        assert_relative_eq!(q_li + q_h, 0.0, epsilon = 1e-9);
    }

    #[test]
    fn test_total_charge_constraint_hydroxide_ion() {
        let solver = get_test_solver();
        let atoms = vec![
            Atom {
                atomic_number: 8,
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                position: [0.96, 0.0, 0.0],
            },
        ];

        let result = solver.solve(&atoms, -1.0).unwrap();
        let (q_o, q_h) = (result.charges[0], result.charges[1]);

        assert_relative_eq!(q_o + q_h, -1.0, epsilon = 1e-9);
    }

    #[test]
    fn test_chemical_trend_nacl_vs_nabr() {
        let solver = get_test_solver();

        let nacl_bond_length = 2.36;
        let nacl_atoms = vec![
            Atom {
                atomic_number: 11,
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 17,
                position: [nacl_bond_length, 0.0, 0.0],
            },
        ];

        let nabr_bond_length = 2.50;
        let nabr_atoms = vec![
            Atom {
                atomic_number: 11,
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 35,
                position: [nabr_bond_length, 0.0, 0.0],
            },
        ];

        let nacl_result = solver.solve(&nacl_atoms, 0.0).unwrap();
        let nabr_result = solver.solve(&nabr_atoms, 0.0).unwrap();

        let charge_na_in_nacl = nacl_result.charges[0];
        let charge_na_in_nabr = nabr_result.charges[0];

        assert!(
            charge_na_in_nacl > charge_na_in_nabr,
            "Chemical trend is incorrect: Charge on Na in NaCl should be greater than in NaBr."
        );

        assert!(charge_na_in_nacl > 0.0);
        assert!(charge_na_in_nabr > 0.0);
    }

    #[test]
    fn test_error_handling_no_atoms() {
        let solver = get_test_solver();
        let atoms: Vec<Atom> = vec![];
        let result = solver.solve(&atoms, 0.0);
        assert!(matches!(result, Err(CheqError::NoAtoms)));
    }

    #[test]
    fn test_error_handling_parameter_not_found() {
        let solver = get_test_solver();
        let atoms = vec![Atom {
            atomic_number: 118,
            position: [0.0, 0.0, 0.0],
        }];
        let result = solver.solve(&atoms, 0.0);
        assert!(matches!(result, Err(CheqError::ParameterNotFound(118))));
    }

    #[test]
    fn test_not_converged_error() {
        let params = get_default_parameters();
        let options = SolverOptions {
            max_iterations: 1,
            tolerance: 1e-20,
            lambda_scale: 1.0,
        };
        let solver = QEqSolver::new(params).with_options(options);
        let atoms = vec![
            Atom {
                atomic_number: 8,
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                position: [0.9575, 0.0, 0.0],
            },
            Atom {
                atomic_number: 1,
                position: [0.9575 * (-0.5), 0.9575 * (3.0_f64).sqrt() / 2.0, 0.0],
            },
        ];
        let result = solver.solve(&atoms, 0.0);
        assert!(matches!(result, Err(CheqError::NotConverged { .. })));
    }

    #[test]
    fn test_linalg_error_singular_matrix() {
        let solver = get_test_solver();
        let atoms = vec![
            Atom {
                atomic_number: 6,
                position: [0.0, 0.0, 0.0],
            },
            Atom {
                atomic_number: 6,
                position: [0.0, 0.0, 0.0],
            },
        ];
        let result = solver.solve(&atoms, 0.0);

        match result {
            Err(CheqError::LinalgError(_)) => (),
            Ok(_) => (),
            _ => panic!("Unexpected error"),
        }
    }

    #[test]
    fn test_with_options() {
        let params = get_default_parameters();
        let options = SolverOptions {
            max_iterations: 100,
            tolerance: 1e-10,
            lambda_scale: 1.0,
        };
        let solver = QEqSolver::new(params).with_options(options);
        let atoms = vec![Atom {
            atomic_number: 6,
            position: [0.0, 0.0, 0.0],
        }];
        let result = solver.solve(&atoms, 0.0).unwrap();
        assert_eq!(result.charges.len(), 1);
    }
}
