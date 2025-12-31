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
    /// The 3D position of the atom in Cartesian coordinates (Ångströms).
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

/// A point charge source external to the QEq system.
///
/// This struct represents an atom from the environment (e.g., a protein residue) that contributes
/// a fixed electrostatic potential to the QEq system but does not participate in charge
/// redistribution. Unlike QEq atoms, these charges remain constant throughout the calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointCharge {
    /// The atomic number of the external atom, used to retrieve screening parameters.
    ///
    /// This allows the solver to compute screened Coulomb integrals using the same STO/GTO
    /// formalism as internal atoms, ensuring physical consistency at short distances.
    pub atomic_number: u8,
    /// The 3D position of the external charge in Cartesian coordinates (Ångströms).
    pub position: [f64; 3],
    /// The fixed partial charge assigned to this atom (elementary charge units).
    ///
    /// This value typically comes from a classical force field (e.g., AMBER, CHARMM) and
    /// remains constant throughout the QEq calculation.
    pub charge: f64,
}

impl PointCharge {
    /// Creates a new point charge with the specified parameters.
    ///
    /// # Arguments
    ///
    /// * `atomic_number` - The atomic number of the element.
    /// * `position` - The 3D Cartesian coordinates in Ångströms.
    /// * `charge` - The fixed partial charge in elementary charge units.
    ///
    /// # Returns
    ///
    /// A new `PointCharge` instance.
    ///
    /// # Examples
    ///
    /// ```
    /// use cheq::PointCharge;
    ///
    /// let oxygen = PointCharge::new(8, [1.234, 5.678, 9.012], -0.5679);
    /// ```
    pub fn new(atomic_number: u8, position: [f64; 3], charge: f64) -> Self {
        Self {
            atomic_number,
            position,
            charge,
        }
    }
}

/// An external electrostatic potential acting on the QEq system.
///
/// This struct describes the electrostatic environment surrounding a molecular fragment
/// undergoing charge equilibration. It enables hybrid QEq/MM calculations where the QEq
/// subsystem (e.g., a ligand) responds to the electrostatic field generated by the
/// environment (e.g., a protein binding pocket) without including the environment atoms
/// in the expensive matrix diagonalization.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ExternalPotential {
    /// Collection of external point charges contributing to the electrostatic environment.
    ///
    /// Each point charge represents an atom from the environment with a fixed partial charge.
    /// The charges are typically derived from classical force fields.
    point_charges: Vec<PointCharge>,

    /// A uniform external electric field applied to the system (V/Å).
    ///
    /// The field is represented as a 3D vector `[Ex, Ey, Ez]`. The electrostatic potential
    /// The electrostatic potential contribution at position `r` is `V = -E · r` (dot product).
    uniform_field: [f64; 3],
}

impl ExternalPotential {
    /// Creates an empty external potential with no contributions.
    ///
    /// This is equivalent to computing QEq in vacuum (the default behavior).
    ///
    /// # Examples
    ///
    /// ```
    /// use cheq::ExternalPotential;
    ///
    /// let empty = ExternalPotential::new();
    /// assert!(empty.is_empty());
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates an external potential from a collection of point charges.
    ///
    /// # Arguments
    ///
    /// * `charges` - A vector of `PointCharge` representing external atoms.
    ///
    /// # Returns
    ///
    /// A new `ExternalPotential` instance containing the specified point charges.
    ///
    /// # Examples
    ///
    /// ```
    /// use cheq::{ExternalPotential, PointCharge};
    ///
    /// let charges = vec![
    ///     PointCharge::new(8, [3.0, 0.0, 0.0], -0.82),
    ///     PointCharge::new(1, [3.5, 0.7, 0.0],  0.41),
    /// ];
    ///
    /// let external = ExternalPotential::from_point_charges(charges);
    /// assert_eq!(external.point_charges().len(), 2);
    /// ```
    pub fn from_point_charges(charges: Vec<PointCharge>) -> Self {
        Self {
            point_charges: charges,
            uniform_field: [0.0, 0.0, 0.0],
        }
    }

    /// Creates an external potential from a uniform electric field.
    ///
    /// # Arguments
    ///
    /// * `field` - The electric field vector `[Ex, Ey, Ez]` in V/Å.
    ///
    /// # Returns
    ///
    /// A new `ExternalPotential` instance with the specified uniform field.
    ///
    /// # Examples
    ///
    /// ```
    /// use cheq::ExternalPotential;
    ///
    /// // Apply a field of 0.1 V/Å along the z-axis (e.g., transmembrane potential)
    /// let external = ExternalPotential::from_uniform_field([0.0, 0.0, 0.1]);
    /// ```
    pub fn from_uniform_field(field: [f64; 3]) -> Self {
        Self {
            point_charges: Vec::new(),
            uniform_field: field,
        }
    }

    /// Sets the point charges for this external potential.
    ///
    /// This method consumes and returns `self`, enabling a builder-style API.
    ///
    /// # Arguments
    ///
    /// * `charges` - A vector of `PointCharge` representing external atoms.
    ///
    /// # Examples
    ///
    /// ```
    /// use cheq::{ExternalPotential, PointCharge};
    ///
    /// let external = ExternalPotential::new()
    ///     .with_point_charges(vec![
    ///         PointCharge::new(6, [0.0, 0.0, 5.0], 0.1),
    ///     ])
    ///     .with_uniform_field([0.0, 0.0, 0.05]);
    /// ```
    pub fn with_point_charges(mut self, charges: Vec<PointCharge>) -> Self {
        self.point_charges = charges;
        self
    }

    /// Sets the uniform electric field for this external potential.
    ///
    /// This method consumes and returns `self`, enabling a builder-style API.
    ///
    /// # Arguments
    ///
    /// * `field` - The electric field vector `[Ex, Ey, Ez]` in V/Å.
    ///
    /// # Examples
    ///
    /// ```
    /// use cheq::ExternalPotential;
    ///
    /// let external = ExternalPotential::new()
    ///     .with_uniform_field([0.0, 0.0, 0.2]);
    /// ```
    pub fn with_uniform_field(mut self, field: [f64; 3]) -> Self {
        self.uniform_field = field;
        self
    }

    /// Returns a reference to the point charges.
    pub fn point_charges(&self) -> &[PointCharge] {
        &self.point_charges
    }

    /// Returns the uniform electric field vector.
    pub fn uniform_field(&self) -> [f64; 3] {
        self.uniform_field
    }

    /// Returns `true` if this external potential has no contributions.
    ///
    /// An empty potential is equivalent to vacuum conditions and will not
    /// affect the QEq calculation.
    pub fn is_empty(&self) -> bool {
        self.point_charges.is_empty()
            && self.uniform_field[0] == 0.0
            && self.uniform_field[1] == 0.0
            && self.uniform_field[2] == 0.0
    }

    /// Returns the total number of point charges.
    pub fn num_point_charges(&self) -> usize {
        self.point_charges.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_charge_new() {
        let pc = PointCharge::new(8, [1.0, 2.0, 3.0], -0.5);

        assert_eq!(pc.atomic_number, 8);
        assert_eq!(pc.position, [1.0, 2.0, 3.0]);
        assert_eq!(pc.charge, -0.5);
    }

    #[test]
    fn test_point_charge_struct_literal() {
        let pc = PointCharge {
            atomic_number: 6,
            position: [0.0, 0.0, 0.0],
            charge: 0.1,
        };

        assert_eq!(pc.atomic_number, 6);
        assert_eq!(pc.charge, 0.1);
    }

    #[test]
    fn test_point_charge_clone() {
        let pc1 = PointCharge::new(7, [1.0, 1.0, 1.0], 0.25);
        let pc2 = pc1;

        assert_eq!(pc1.charge, pc2.charge);
        assert_eq!(pc1.position, pc2.position);
    }

    #[test]
    fn test_external_potential_new() {
        let ext = ExternalPotential::new();

        assert!(ext.is_empty());
        assert_eq!(ext.num_point_charges(), 0);
        assert_eq!(ext.uniform_field(), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_external_potential_default() {
        let ext = ExternalPotential::default();

        assert!(ext.is_empty());
    }

    #[test]
    fn test_external_potential_from_point_charges() {
        let charges = vec![
            PointCharge::new(8, [1.0, 0.0, 0.0], -0.5),
            PointCharge::new(1, [2.0, 0.0, 0.0], 0.25),
        ];

        let ext = ExternalPotential::from_point_charges(charges);

        assert!(!ext.is_empty());
        assert_eq!(ext.num_point_charges(), 2);
        assert_eq!(ext.point_charges()[0].atomic_number, 8);
        assert_eq!(ext.point_charges()[1].charge, 0.25);
    }

    #[test]
    fn test_external_potential_from_uniform_field() {
        let ext = ExternalPotential::from_uniform_field([0.1, 0.2, 0.3]);

        assert!(!ext.is_empty());
        assert_eq!(ext.num_point_charges(), 0);
        assert_eq!(ext.uniform_field(), [0.1, 0.2, 0.3]);
    }

    #[test]
    fn test_external_potential_builder_pattern() {
        let charges = vec![PointCharge::new(6, [0.0, 0.0, 0.0], 0.1)];

        let ext = ExternalPotential::new()
            .with_point_charges(charges)
            .with_uniform_field([0.0, 0.0, 1.0]);

        assert!(!ext.is_empty());
        assert_eq!(ext.num_point_charges(), 1);
        assert_eq!(ext.uniform_field()[2], 1.0);
    }

    #[test]
    fn test_external_potential_is_empty_with_only_field() {
        let ext = ExternalPotential::from_uniform_field([0.0, 0.0, 0.001]);
        assert!(!ext.is_empty());
    }

    #[test]
    fn test_external_potential_is_empty_with_only_charges() {
        let ext =
            ExternalPotential::from_point_charges(vec![PointCharge::new(1, [0.0, 0.0, 0.0], 0.0)]);
        assert!(!ext.is_empty());
    }

    #[test]
    fn test_atom_view_implementation() {
        let atom = Atom {
            atomic_number: 6,
            position: [1.0, 2.0, 3.0],
        };

        assert_eq!(atom.atomic_number(), 6);
        assert_eq!(atom.position(), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_atom_copy_clone() {
        let atom1 = Atom {
            atomic_number: 8,
            position: [0.0, 0.0, 0.0],
        };
        let atom2 = atom1;

        assert_eq!(atom1.atomic_number, atom2.atomic_number);
        assert_eq!(atom1.position, atom2.position);
    }
}
