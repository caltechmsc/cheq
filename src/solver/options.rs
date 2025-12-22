//! This module defines configuration options for the charge equilibration solver.
//!
//! It provides the `SolverOptions` struct, which allows users to customize the convergence
//! criteria, iteration limits, and screening parameters for the QEq algorithm. These options
//! control the trade-off between computational cost and solution accuracy.

/// Specifies the type of basis functions used for Coulomb integrals.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BasisType {
    /// Gaussian-Type Orbitals (GTO).
    ///
    /// Uses Gaussian approximations for Slater orbitals. This allows for fast analytical
    /// integration but introduces a small approximation error.
    Gto,

    /// Slater-Type Orbitals (STO).
    ///
    /// Uses exact Slater orbitals. This is more accurate but computationally more expensive
    /// as it requires evaluating complex analytical expansions.
    #[default]
    Sto,
}

/// Strategy for damping charge updates during SCF iterations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DampingStrategy {
    /// No damping is applied (equivalent to damping = 1.0).
    ///
    /// Fastest convergence for simple, well-behaved systems, but prone to oscillation
    /// or divergence in complex cases (e.g., LiH, large systems).
    None,

    /// A fixed damping factor (0.0 < d <= 1.0).
    ///
    /// A constant mixing rate. Lower values (e.g., 0.1) are more stable but slower.
    /// Higher values (e.g., 0.8) are faster but risk instability.
    Fixed(f64),

    /// Automatically adjusts damping based on convergence behavior.
    ///
    /// Starts with an initial value. If the error increases (divergence), it reduces damping (brakes).
    /// If the error decreases significantly (convergence), it increases damping (accelerates).
    /// This is the recommended strategy for general use.
    Auto { initial: f64 },
}

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
    /// It corresponds to the `λ` parameter in the Rappe & Goddard paper.
    /// The original paper suggests a value of ~0.4913, often rounded to 0.5.
    pub lambda_scale: f64,

    /// Whether to update hydrogen hardness each iteration (nonlinear SCF term).
    ///
    /// When disabled, hydrogen uses a fixed hardness (lossy simplification). Enabled by default.
    pub hydrogen_scf: bool,

    /// The type of basis functions to use for Coulomb integrals.
    ///
    /// Defaults to `BasisType::Sto` (Slater-Type Orbitals) for maximum accuracy.
    pub basis_type: BasisType,

    /// The damping strategy for the SCF iteration.
    ///
    /// Controls how charge updates are mixed between iterations to ensure convergence.
    pub damping: DampingStrategy,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            tolerance: 1.0e-6,
            max_iterations: 2000,
            lambda_scale: 0.5,
            hydrogen_scf: true,
            basis_type: BasisType::default(),
            damping: DampingStrategy::Auto { initial: 0.4 },
        }
    }
}
