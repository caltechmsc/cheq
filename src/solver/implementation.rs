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
use std::collections::HashMap;
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
    /// let options = SolverOptions { tolerance: 1e-8, max_iterations: 50, lambda_scale: 1.0, ..SolverOptions::default() };
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
        let mut charges = Col::zeros(n_atoms);
        let hydrogen_scf = self.options.hydrogen_scf && has_hydrogen;
        let inner_iters = if hydrogen_scf {
            self.options.hydrogen_inner_iters
        } else {
            0
        };

        if !hydrogen_scf {
            let (_, equilibrated_potential) =
                self.run_single_solve(&invariant, &mut charges, hydrogen_scf)?;
            return Ok(CalculationResult {
                charges: charges.as_ref().iter().cloned().collect(),
                equilibrated_potential,
                iterations: 1,
            });
        }

        let mut max_charge_delta = 0.0;

        for iteration in 1..=self.options.max_iterations {
            if inner_iters > 0 {
                for _ in 0..inner_iters {
                    let (delta, _) =
                        self.run_single_solve(&invariant, &mut charges, hydrogen_scf)?;
                    if delta < self.options.tolerance {
                        break;
                    }
                }
            }

            let (delta, equilibrated_potential) =
                self.run_single_solve(&invariant, &mut charges, hydrogen_scf)?;
            max_charge_delta = delta;

            if !hydrogen_scf || max_charge_delta < self.options.tolerance {
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

    /// Precomputes the geometry-invariant parts of the linear system.
    ///
    /// Builds the base coefficient matrix (diagonal hardness, screened Coulomb off-diagonals,
    /// and charge conservation row/col) plus the RHS vector. Hydrogen diagonal metadata is
    /// stored so the per-iteration charge-dependent hardness update only touches a few entries.
    fn build_invariant_system<A: AtomView>(
        &self,
        atoms: &[A],
        element_data: &[&'p ElementData],
        total_charge: f64,
    ) -> Result<InvariantSystem, CheqError> {
        let n_atoms = atoms.len();
        let matrix_size = n_atoms + 1;

        let mut base_matrix = Mat::zeros(matrix_size, matrix_size);
        let mut rhs = Col::zeros(matrix_size);

        let mut hydrogen_meta = Vec::new();

        for i in 0..n_atoms {
            let data_i = element_data[i];
            base_matrix[(i, i)] = data_i.hardness;
            rhs[i] = -data_i.electronegativity;

            if data_i.principal_quantum_number == 1 {
                hydrogen_meta.push((i, data_i.hardness));
            }
        }

        let positions: Vec<[f64; 3]> = atoms.iter().map(AtomView::position).collect();

        let off_diagonals = compute_off_diagonal_rows(
            &positions,
            element_data,
            self.options.lambda_scale,
            self.options.cutoff_radius,
        );
        for (i, entries) in off_diagonals {
            for (j, val) in entries {
                base_matrix[(i, j)] = val;
                base_matrix[(j, i)] = val;
            }
        }

        base_matrix
            .col_mut(matrix_size - 1)
            .subrows_mut(0, n_atoms)
            .fill(-1.0);
        base_matrix
            .row_mut(matrix_size - 1)
            .subcols_mut(0, n_atoms)
            .fill(1.0);
        rhs[matrix_size - 1] = total_charge;

        Ok(InvariantSystem {
            base_matrix,
            rhs,
            hydrogen_meta,
        })
    }

    /// Performs a single SCF iteration to solve the linear system and update charges.
    fn run_single_solve(
        &self,
        invariant: &InvariantSystem,
        charges: &mut Col<f64>,
        hydrogen_scf: bool,
    ) -> Result<(f64, f64), CheqError> {
        let n_atoms = charges.nrows();
        let mut work_matrix = invariant.base_matrix.clone();

        if hydrogen_scf {
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
                    "Linear system is likely singular or ill-conditioned due to input geometry."
                        .to_string(),
                ));
            }
        };

        let new_charges = solution.as_ref().subrows(0, n_atoms);
        let max_charge_delta = new_charges
            .as_ref()
            .iter()
            .zip(charges.as_ref().iter())
            .map(|(new, old): (&f64, &f64)| (*new - *old).abs())
            .fold(0.0, f64::max);
        charges.as_mut().copy_from(&new_charges);

        let equilibrated_potential = -solution[n_atoms];
        Ok((max_charge_delta, equilibrated_potential))
    }
}

/// Geometry-invariant linear system components reused across SCF iterations.
struct InvariantSystem {
    base_matrix: Mat<f64>,
    rhs: Col<f64>,
    hydrogen_meta: Vec<(usize, f64)>,
}

/// Computes screened Coulomb off-diagonal entries for the upper triangle.
fn compute_off_diagonal_rows(
    positions: &[[f64; 3]],
    element_data: &[&ElementData],
    lambda_scale: f64,
    cutoff_radius: Option<f64>,
) -> Vec<(usize, Vec<(usize, f64)>)> {
    match cutoff_radius {
        None => compute_dense_off_diagonals(positions, element_data, lambda_scale),
        Some(radius) => compute_cutoff_off_diagonals(positions, element_data, lambda_scale, radius),
    }
}

/// Computes all off-diagonal entries without cutoff.
fn compute_dense_off_diagonals(
    positions: &[[f64; 3]],
    element_data: &[&ElementData],
    lambda_scale: f64,
) -> Vec<(usize, Vec<(usize, f64)>)> {
    let n_atoms = positions.len();

    (0..n_atoms)
        .into_par_iter()
        .map(|i| {
            let data_i = element_data[i];
            let pos_i = positions[i];
            let radius_i_bohr = data_i.radius / math::constants::BOHR_TO_ANGSTROM;
            let mut row_vals = Vec::with_capacity(n_atoms.saturating_sub(i + 1));

            for j in (i + 1)..n_atoms {
                let val = compute_pair_entry(
                    pos_i,
                    positions[j],
                    data_i,
                    element_data[j],
                    radius_i_bohr,
                    lambda_scale,
                );
                row_vals.push((j, val));
            }

            (i, row_vals)
        })
        .collect()
}

/// Computes off-diagonal entries with a hard cutoff radius.
fn compute_cutoff_off_diagonals(
    positions: &[[f64; 3]],
    element_data: &[&ElementData],
    lambda_scale: f64,
    cutoff_radius: f64,
) -> Vec<(usize, Vec<(usize, f64)>)> {
    let n_atoms = positions.len();
    let cell_size = cutoff_radius;
    let mut cells: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();

    for (idx, pos) in positions.iter().enumerate() {
        let key = (
            (pos[0] / cell_size).floor() as i64,
            (pos[1] / cell_size).floor() as i64,
            (pos[2] / cell_size).floor() as i64,
        );
        cells.entry(key).or_default().push(idx);
    }

    let cells = std::sync::Arc::new(cells);
    let cutoff_sq = cutoff_radius * cutoff_radius;

    (0..n_atoms)
        .into_par_iter()
        .map(|i| {
            let data_i = element_data[i];
            let pos_i = positions[i];
            let radius_i_bohr = data_i.radius / math::constants::BOHR_TO_ANGSTROM;
            let key = (
                (pos_i[0] / cell_size).floor() as i64,
                (pos_i[1] / cell_size).floor() as i64,
                (pos_i[2] / cell_size).floor() as i64,
            );

            let mut entries = Vec::new();

            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let neighbor_key = (key.0 + dx, key.1 + dy, key.2 + dz);
                        if let Some(indices) = cells.get(&neighbor_key) {
                            for &j in indices {
                                if j <= i {
                                    continue;
                                }
                                let dist_sq: f64 = pos_i
                                    .iter()
                                    .zip(positions[j].iter())
                                    .map(|(pi, pj)| {
                                        let diff = pi - pj;
                                        diff * diff
                                    })
                                    .sum();
                                if dist_sq > cutoff_sq {
                                    continue;
                                }
                                let val = compute_pair_entry(
                                    pos_i,
                                    positions[j],
                                    data_i,
                                    element_data[j],
                                    radius_i_bohr,
                                    lambda_scale,
                                );
                                entries.push((j, val));
                            }
                        }
                    }
                }
            }

            entries.sort_by_key(|(j, _)| *j);
            (i, entries)
        })
        .collect()
}

/// Computes the screened Coulomb interaction entry between two atoms.
fn compute_pair_entry(
    pos_i: [f64; 3],
    pos_j: [f64; 3],
    data_i: &ElementData,
    data_j: &ElementData,
    radius_i_bohr: f64,
    lambda_scale: f64,
) -> f64 {
    let dist_sq: f64 = pos_i
        .iter()
        .zip(pos_j.iter())
        .map(|(pi, pj)| {
            let diff = pi - pj;
            diff * diff
        })
        .sum();
    let distance_angstrom = dist_sq.sqrt();
    let distance_bohr = distance_angstrom / math::constants::BOHR_TO_ANGSTROM;
    let radius_j_bohr = data_j.radius / math::constants::BOHR_TO_ANGSTROM;

    let j_ij_hartree = math::shielding::gaussian_coulomb_integral(
        distance_bohr,
        data_i.principal_quantum_number,
        radius_i_bohr,
        data_j.principal_quantum_number,
        radius_j_bohr,
        lambda_scale,
    );

    j_ij_hartree * math::constants::HARTREE_TO_EV
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
            ..SolverOptions::default()
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
            ..SolverOptions::default()
        };
        let solver = QEqSolver::new(params).with_options(options);
        let atoms = vec![Atom {
            atomic_number: 6,
            position: [0.0, 0.0, 0.0],
        }];
        let result = solver.solve(&atoms, 0.0).unwrap();
        assert_eq!(result.charges.len(), 1);
    }

    #[test]
    fn test_cutoff_large_matches_dense() {
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

        let base_solver = get_test_solver();
        let base = base_solver.solve(&atoms, 0.0).unwrap();

        let params = get_default_parameters();
        let cutoff_solver = QEqSolver::new(params).with_options(SolverOptions {
            cutoff_radius: Some(10.0),
            ..SolverOptions::default()
        });
        let cutoff = cutoff_solver.solve(&atoms, 0.0).unwrap();

        assert_relative_eq!(cutoff.charges[0], base.charges[0], epsilon = 1e-9);
        assert_relative_eq!(cutoff.charges[1], base.charges[1], epsilon = 1e-9);
        assert_relative_eq!(cutoff.charges.iter().sum::<f64>(), 0.0, epsilon = 1e-9);
    }

    #[test]
    fn test_hydrogen_inner_iters_respects_charge() {
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

        let base = get_test_solver().solve(&atoms, -1.0).unwrap();

        let params = get_default_parameters();
        let solver = QEqSolver::new(params).with_options(SolverOptions {
            hydrogen_inner_iters: 2,
            ..SolverOptions::default()
        });
        let result = solver.solve(&atoms, -1.0).unwrap();

        assert_relative_eq!(result.charges.iter().sum::<f64>(), -1.0, epsilon = 1e-9);
        assert_relative_eq!(result.charges[0], base.charges[0], epsilon = 1e-5);
        assert_relative_eq!(result.charges[1], base.charges[1], epsilon = 1e-5);
    }

    #[test]
    fn test_hydrogen_scf_off_converges_in_one_iteration() {
        let params = get_default_parameters();
        let solver = QEqSolver::new(params).with_options(SolverOptions {
            hydrogen_scf: false,
            ..SolverOptions::default()
        });

        let bond_length = 0.9575;
        let bond_angle_rad = 104.45f64.to_radians();
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

        assert_eq!(result.iterations, 1);
        assert_relative_eq!(result.charges.iter().sum::<f64>(), 0.0, epsilon = 1e-9);
        assert_relative_eq!(result.charges[1], result.charges[2], epsilon = 1e-9);
    }
}
