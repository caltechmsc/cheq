//! This module contains the solver components for performing charge equilibration calculations.
//!
//! It includes the `QEqSolver` implementation and `SolverOptions` for configuring the solver,
//! providing the core functionality for the QEq method in the `cheq` library.

mod implementation;
mod options;

/// The `QEqSolver` struct for performing charge equilibration calculations.
///
/// This is the main solver that implements the SCF iterative procedure to compute partial
/// atomic charges based on the QEq method.
pub use implementation::QEqSolver;

/// The `SolverOptions` struct for configuring solver parameters.
///
/// This struct allows customization of convergence tolerance, maximum iterations, and
/// shielding parameters for the charge equilibration process.
pub use options::SolverOptions;
