//! This module contains the solver components for performing charge equilibration calculations.
//!
//! It includes the `QEqSolver` implementation and `SolverOptions` for configuring the solver,
//! providing the core functionality for the QEq method in the `cheq` library.

mod implementation;
mod options;

pub use implementation::QEqSolver;
pub use options::{BasisType, DampingStrategy, SolverOptions};
