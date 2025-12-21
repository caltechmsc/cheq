//! This module implements the core `QEqSolver` for performing charge equilibration calculations.
//!
//! The `QEqSolver` encapsulates the Self-Consistent Field (SCF) iterative procedure that solves
//! the QEq equations to determine partial atomic charges. It constructs a linear system based on
//! atomic parameters and geometry, iteratively refining charges until convergence. The solver
//! integrates with the broader `cheq` architecture by using the `AtomView` trait for atom data
//! and `Parameters` for element-specific values, enabling decoupled and flexible molecular
//! simulations.

use super::options::{BasisType, SolverOptions};
use crate::{
    error::CheqError,
    params::{ElementData, Parameters},
    shielding::{self, constants},
    types::{AtomView, CalculationResult},
};
use faer::{Col, Mat, prelude::*};
use rayon::prelude::*;
use std::panic::{self, AssertUnwindSafe};

/// Charge dependence factor for hydrogen hardness, derived from empirical fitting.
/// This constant adjusts the hardness of hydrogen atoms based on their partial charge,
/// introducing non-linearity into the QEq equations.
const H_CHARGE_DEPENDENCE_FACTOR: f64 = 0.93475415965;

/// A thread-safe wrapper for raw matrix access to enable parallel filling.
///
/// This struct allows multiple threads to write to disjoint parts of a matrix
/// without locking, which is safe because we ensure unique indices in the parallel iterator.
struct UnsafeMatView {
    ptr: *mut f64,
    row_stride: isize,
    col_stride: isize,
}

unsafe impl Send for UnsafeMatView {}
unsafe impl Sync for UnsafeMatView {}

impl UnsafeMatView {
    /// Writes a value to the matrix at the specified (row, col) index.
    ///
    /// # Safety
    ///
    /// The caller must ensure that:
    /// 1. The (row, col) indices are within bounds.
    /// 2. No other thread is writing to the same address simultaneously.
    unsafe fn write(&self, row: usize, col: usize, val: f64) {
        let offset = (row as isize) * self.row_stride + (col as isize) * self.col_stride;
        unsafe {
            *self.ptr.offset(offset) = val;
        }
    }
}

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
    /// use cheq::get_default_parameters;
    /// use cheq::QEqSolver;
    ///
    /// let params = get_default_parameters();
    /// let solver = QEqSolver::new(params);
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
    /// use cheq::get_default_parameters;
    /// use cheq::{QEqSolver, SolverOptions};
    ///
    /// let params = get_default_parameters();
    /// let options = SolverOptions {
    ///     max_iterations: 100,
    ///     tolerance: 1e-6,
    ///     ..Default::default()
    /// };
    ///
    /// let solver = QEqSolver::new(params).with_options(options);
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
    /// # Examples
    ///
    /// ```
    /// use cheq::get_default_parameters;
    /// use cheq::QEqSolver;
    /// use cheq::Atom;
    ///
    /// // 1. Setup parameters and solver
    /// let params = get_default_parameters();
    /// let solver = QEqSolver::new(params);
    ///
    /// // 2. Define a molecule (e.g., H2)
    /// let atoms = vec![
    ///     Atom { atomic_number: 1, position: [0.0, 0.0, 0.0] },
    ///     Atom { atomic_number: 1, position: [0.74, 0.0, 0.0] },
    /// ];
    ///
    /// // 3. Run calculation
    /// let result = solver.solve(&atoms, 0.0).unwrap();
    ///
    /// assert_eq!(result.charges.len(), 2);
    /// println!("Charges: {:?}", result.charges);
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

        let mut charges = Col::zeros(n_atoms);
        let has_hydrogen = !invariant.hydrogen_meta.is_empty();
        let hydrogen_scf = self.options.hydrogen_scf && has_hydrogen;

        let mut work_matrix = invariant.base_matrix.clone();

        if !hydrogen_scf {
            let (_, equilibrated_potential) =
                self.run_single_solve(&invariant, &mut work_matrix, &mut charges, false)?;
            return Ok(CalculationResult {
                charges: charges.as_ref().iter().cloned().collect(),
                equilibrated_potential,
                iterations: 1,
            });
        }

        let mut max_charge_delta = 0.0;

        for iteration in 1..=self.options.max_iterations {
            if iteration > 1 {
                work_matrix.copy_from(&invariant.base_matrix);
            }

            let (delta, equilibrated_potential) =
                self.run_single_solve(&invariant, &mut work_matrix, &mut charges, true)?;

            max_charge_delta = delta;

            if max_charge_delta < self.options.tolerance {
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

        let mat_view = UnsafeMatView {
            ptr: base_matrix.as_ptr_mut(),
            row_stride: base_matrix.row_stride(),
            col_stride: base_matrix.col_stride(),
        };

        let lambda = self.options.lambda_scale;
        let basis_type = self.options.basis_type;

        (0..n_atoms).into_par_iter().for_each(|i| {
            let data_i = element_data[i];
            let pos_i = positions[i];
            let radius_i_bohr = data_i.radius / constants::BOHR_TO_ANGSTROM;

            for j in (i + 1)..n_atoms {
                let pos_j = positions[j];
                let diff_sq = (pos_i[0] - pos_j[0]).powi(2)
                    + (pos_i[1] - pos_j[1]).powi(2)
                    + (pos_i[2] - pos_j[2]).powi(2);

                let dist_angstrom = diff_sq.sqrt();
                let dist_bohr = dist_angstrom / constants::BOHR_TO_ANGSTROM;

                let data_j = element_data[j];
                let radius_j_bohr = data_j.radius / constants::BOHR_TO_ANGSTROM;

                let integral_hartree = match basis_type {
                    BasisType::Gto => shielding::gto::calculate_integral(
                        dist_bohr,
                        data_i.principal_quantum_number,
                        radius_i_bohr,
                        data_j.principal_quantum_number,
                        radius_j_bohr,
                        lambda,
                    ),
                    BasisType::Sto => shielding::sto::calculate_integral(
                        dist_bohr,
                        data_i.principal_quantum_number,
                        radius_i_bohr,
                        data_j.principal_quantum_number,
                        radius_j_bohr,
                        lambda,
                    ),
                };

                let val_ev = integral_hartree * constants::HARTREE_TO_EV;

                // SAFETY: Each unordered pair (i, j) with i < j is handled only by the thread for i.
                // That thread writes (i, j) and (j, i), so no two threads write the same entries.
                unsafe {
                    mat_view.write(i, j, val_ev);
                    mat_view.write(j, i, val_ev);
                }
            }
        });

        base_matrix
            .col_mut(matrix_size - 1)
            .subrows_mut(0, n_atoms)
            .fill(1.0);
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
        work_matrix: &mut Mat<f64>,
        charges: &mut Col<f64>,
        hydrogen_scf: bool,
    ) -> Result<(f64, f64), CheqError> {
        let n_atoms = charges.nrows();

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
                    "Linear system solver panicked. Matrix might be singular.".to_string(),
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

        let equilibrated_potential = solution[n_atoms];

        Ok((max_charge_delta, equilibrated_potential))
    }
}

/// Geometry-invariant linear system components reused across SCF iterations.
struct InvariantSystem {
    base_matrix: Mat<f64>,
    rhs: Col<f64>,
    hydrogen_meta: Vec<(usize, f64)>,
}
