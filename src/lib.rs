//! A high-performance Rust library for calculating dynamic partial atomic charges
//! using the Charge Equilibration (QEq) method.
//!
//! # Cheq: Charge Equilibration for Molecular Dynamics
//!
//! **Cheq** implements the robust and general Charge Equilibration (QEq) method
//! proposed by Rappé and Goddard (1991) to predict the self-consistent,
//! geometry-dependent partial charges of atoms in a molecular system.
//!
//! Unlike fixed-charge force fields, QEq allows atomic charges to dynamically
//! respond to changes in molecular geometry and external fields, which is
//! essential for accurate molecular dynamics (MD) simulations, particularly
//! for systems like ceramics, polymers, and biological molecules where
//! charge transfer is significant.
//!
//! The method relies solely on atomic properties (electronegativity, hardness,
//! and radius) and the molecular geometry to achieve charge neutrality.
//!
//! # Quickstart
//!
//! The primary entry point is the [`QEqSolver`], which requires a set of [`Parameters`]
//! (atomic data) and a slice of atoms with coordinates.
//!
//! ```
//! use cheq::{get_default_parameters, QEqSolver, Atom, SolverOptions};
//! use approx::assert_relative_eq;
//!
//! // 1. Set up the parameters and solver.
//! //    (Parameters include chi, J, radius, and n for each element)
//! let params = get_default_parameters();
//! let options = SolverOptions::default(); // Use default convergence criteria
//! let solver = QEqSolver::new(params).with_options(options);
//!
//! // 2. Define the molecular system (Water molecule geometry).
//! let bond_length = 0.9575;
//! let angle = 104.45f64.to_radians();
//!
//! let atoms = vec![
//!     Atom { atomic_number: 8, position: [0.0, 0.0, 0.0] },         // Oxygen
//!     Atom { atomic_number: 1, position: [bond_length, 0.0, 0.0] }, // Hydrogen 1
//!     Atom { atomic_number: 1, position: [
//!         bond_length * angle.cos(),
//!         bond_length * angle.sin(),
//!         0.0,
//!     ]},                                                           // Hydrogen 2
//! ];
//!
//! // 3. Solve for partial charges (total charge Q_tot = 0.0).
//! let result = solver.solve(&atoms, 0.0).unwrap();
//!
//! // 4. Inspect the results.
//! let q_o = result.charges[0];
//! let q_h1 = result.charges[1];
//!
//! assert!(q_o < 0.0);   // Oxygen is negative
//! assert!(q_h1 > 0.0);  // Hydrogen is positive
//! assert_relative_eq!(q_o + 2.0 * q_h1, 0.0, epsilon = 1e-9); // Charge conservation
//! println!("Charges (O, H1, H2): ({:.3}, {:.3}, {:.3})", q_o, q_h1, result.charges[2]);
//! ```

// Internal modules are private to enforce encapsulation.
mod error;
mod params;
mod shielding;
mod solver;
mod types;

// Re-export public components at the crate root for easy access.

pub use self::error::CheqError;

pub use self::types::{Atom, AtomView, CalculationResult};

pub use self::solver::{BasisType, QEqSolver, SolverOptions};

pub use self::params::{ElementData, Parameters};

// Provide access to default parameters.
use std::sync::OnceLock;
static DEFAULT_PARAMETERS: OnceLock<Parameters> = OnceLock::new();

/// Provides a thread-safe, cached reference to the default atomic parameters.
///
/// These parameters are statically compiled into the library from the `resources/qeq.data.toml`
/// file and include data for a wide range of elements, following the conventions
/// derived from the Rappé & Goddard QEq formalism.
///
/// This function ensures that the parameters are loaded only once, improving
/// performance for repeated calls.
///
/// # Returns
///
/// A static reference (`&'static Parameters`) to the default parameter set.
///
/// # Panics
///
/// Panics only if the embedded default TOML file cannot be parsed, which indicates
/// a bug in the library distribution.
pub fn get_default_parameters() -> &'static Parameters {
    DEFAULT_PARAMETERS.get_or_init(|| {
        const DEFAULT_PARAMS_TOML: &str = include_str!("../resources/qeq.data.toml");
        Parameters::load_from_str(DEFAULT_PARAMS_TOML)
            .expect("Failed to parse embedded default parameters. This is a library bug.")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_default_parameters() {
        let params1 = get_default_parameters();
        assert!(
            params1.elements.contains_key(&6),
            "Carbon (6) should be present"
        );
        assert!(
            params1.elements.contains_key(&8),
            "Oxygen (8) should be present"
        );

        let params2 = get_default_parameters();
        assert_eq!(
            params1 as *const _, params2 as *const _,
            "Subsequent calls should return a cached reference"
        );
    }
}
