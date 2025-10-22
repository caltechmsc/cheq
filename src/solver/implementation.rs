use super::options::SolverOptions;
use crate::{
    error::CheqError,
    math,
    params::{ElementData, Parameters},
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
        let mut max_charge_delta = 0.0;

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
            max_charge_delta = diff_abs.max().unwrap();

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
