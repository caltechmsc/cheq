//! This module defines configuration options for the charge equilibration solver.
//!
//! It provides the `SolverOptions` struct, which allows users to customize the convergence
//! criteria, iteration limits, and screening parameters for the QEq algorithm. These options
//! control the trade-off between computational cost and solution accuracy.

/// Configuration parameters for the charge equilibration solver.
///
/// This struct encapsulates the numerical settings that control the iterative solution process
/// in the QEq algorithm. Users can adjust these parameters to balance convergence speed with
/// accuracy for different molecular systems and computational requirements.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverOptions {
    /// The convergence tolerance for charge equilibration.
    ///
    /// The solver iterates until the root-mean-square change in atomic charges between
    /// successive iterations falls below this threshold. Smaller values yield more accurate
    /// results but require more iterations.
    pub tolerance: f64,
    /// The maximum number of iterations allowed.
    ///
    /// If convergence is not achieved within this limit, the solver will terminate and return
    /// the current best estimate. This prevents infinite loops in difficult cases.
    pub max_iterations: u32,
    /// The screening parameter scale factor.
    ///
    /// This multiplier adjusts the orbital screening strength in Coulomb integral calculations.
    /// Values typically range from 0.4 to 0.6, with smaller values corresponding to stronger
    /// screening and more compact effective potentials.
    pub lambda_scale: f64,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-6,
            max_iterations: 100,
            lambda_scale: 0.5,
        }
    }
}
