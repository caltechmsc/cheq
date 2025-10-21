//! This module defines the core types used in the cheq library for representing atoms and calculation results.
//!
//! It includes the `AtomView` trait for abstracting atom data access, the `Atom` struct for concrete atom
//! representation, and the `CalculationResult` struct for storing the outcomes of charge equilibration
//! calculations. These types form the foundation for the decoupled design that allows integration with
//! various molecular data structures.

/// A trait for viewing atom data without owning it.
///
/// This trait provides a common interface for accessing an atom's atomic number and 3D position,
/// enabling the charge equilibration solver to work with different atom representations. By decoupling
/// the solver from specific data structures, users can integrate the `cheq` library with their own
/// molecular representations without data conversion overhead.
pub trait AtomView {
    /// Returns the atomic number of the atom.
    ///
    /// The atomic number uniquely identifies the chemical element and is used to look up atomic
    /// parameters such as electronegativity and hardness.
    fn atomic_number(&self) -> u8;

    /// Returns the 3D position of the atom in Cartesian coordinates.
    ///
    /// The position is represented as an array of three `f64` values corresponding to x, y, and z
    /// coordinates. This information is crucial for calculating interatomic distances and Coulomb
    /// interactions in the charge equilibration method.
    fn position(&self) -> [f64; 3];
}

/// A concrete representation of an atom with atomic number and position.
///
/// This struct provides a simple, owned implementation of the `AtomView` trait. It can be used
/// directly for basic atom representations or as a building block for more complex atom types that
/// include additional properties like velocities or forces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Atom {
    /// The atomic number of the atom, identifying its chemical element.
    pub atomic_number: u8,
    /// The 3D position of the atom in Cartesian coordinates.
    pub position: [f64; 3],
}

impl AtomView for Atom {
    #[inline(always)]
    fn atomic_number(&self) -> u8 {
        self.atomic_number
    }

    #[inline(always)]
    fn position(&self) -> [f64; 3] {
        self.position
    }
}

/// The result of a charge equilibration calculation.
///
/// This struct encapsulates the output of a successful charge equilibration run, including the
/// computed partial atomic charges, the equilibrated chemical potential, and diagnostic information
/// about the iterative solution process.
#[derive(Debug, Clone, PartialEq)]
pub struct CalculationResult {
    /// The computed partial atomic charges for each atom in the system.
    ///
    /// Charges are stored in the same order as the input atoms. The sum of all charges equals the
    /// total system charge specified in the calculation.
    pub charges: Vec<f64>,
    /// The equilibrated chemical potential achieved at convergence.
    ///
    /// This value represents the uniform chemical potential across all atoms when the system has
    /// reached charge equilibration.
    pub equilibrated_potential: f64,
    /// The number of iterations performed to reach convergence.
    ///
    /// This provides insight into the computational effort required and can help diagnose
    /// convergence issues in difficult systems.
    pub iterations: u32,
}
