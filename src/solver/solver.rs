use crate::{
    error::CheqError,
    math,
    params::{ElementData, Parameters},
    solver::options::SolverOptions,
    types::{AtomView, CalculationResult},
};
use faer::{Col, ColRef, Mat, prelude::*};
use std::panic::{self, AssertUnwindSafe};

const H_CHARGE_DEPENDENCE_FACTOR: f64 = 0.93475415965;

pub struct QEqSolver<'p> {
    parameters: &'p Parameters,
    options: SolverOptions,
}

impl<'p> QEqSolver<'p> {
    pub fn new(parameters: &'p Parameters) -> Self {
        Self {
            parameters,
            options: SolverOptions::default(),
        }
    }

    pub fn with_options(mut self, options: SolverOptions) -> Self {
        self.options = options;
        self
    }

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
        let has_hydrogen = element_data
            .iter()
            .any(|data| data.principal_quantum_number == 1);

        let mut charges = Col::zeros(n_atoms);

        for iteration in 1..=self.options.max_iterations {
            let (a, b) = self.build_system(atoms, &element_data, total_charge, charges.as_ref())?;

            let solve_result =
                panic::catch_unwind(AssertUnwindSafe(|| a.partial_piv_lu().solve(&b)));

            let solution = match solve_result {
                Ok(sol) => sol,
                Err(_) => {
                    return Err(CheqError::LinalgError(
                        "Linear system is likely singular or ill-conditioned due to input geometry.".to_string(),
                    ));
                }
            };

            let new_charges = solution.as_ref().subrows(0, n_atoms);

            let diff_abs = Col::from_fn(n_atoms, |i| (new_charges[i] - charges[i]).abs());
            let max_charge_delta = diff_abs.max().unwrap();

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
            delta: Col::from_fn(n_atoms, |i| charges[i].abs()).max().unwrap(),
        })
    }

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

                let j_ij_ev = j_ij_hartree * math::constants::HARTREE_TO_EV;
                a_mut[(i, j)] = j_ij_ev;
                a_mut[(j, i)] = j_ij_ev;
            }

            b_mut[i] = -data_i.electronegativity;
        }

        a.col_mut(matrix_size - 1)
            .subrows_mut(0, n_atoms)
            .fill(-1.0);
        a.row_mut(matrix_size - 1).subcols_mut(0, n_atoms).fill(1.0);

        b_mut[matrix_size - 1] = total_charge;

        Ok((a, b))
    }
}
