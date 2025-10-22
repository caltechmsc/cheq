//! This module defines fundamental physical and mathematical constants used throughout the cheq library.
//!
//! These constants are essential for unit conversions and numerical thresholds in charge equilibration
//! calculations, ensuring consistent handling of atomic units, energy scales, and geometric tolerances.

/// Conversion factor from Bohr radii to angstroms.
///
/// This constant represents the standard conversion factor between atomic units (Bohr) and
/// angstroms for length measurements. It is used to convert interatomic distances from the
/// internal Bohr units to angstroms for user-facing coordinates and radii.
///
/// The value is approximately 0.529 Å per Bohr radius.
pub const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;

/// Conversion factor from Hartree energy units to electron volts.
///
/// This constant provides the conversion between atomic energy units (Hartree) and electron
/// volts, which are more commonly used in chemical contexts. It is applied to convert
/// electronegativity, hardness, and other energy-related parameters between these units.
///
/// The value is approximately 27.21 eV per Hartree.
pub const HARTREE_TO_EV: f64 = 27.211_386_245_988;

/// Threshold for considering two distances equal in Bohr units.
///
/// This small numerical threshold is used in distance comparisons and geometric operations
/// to account for floating-point precision limitations. Distances smaller than this value
/// are considered negligible, preventing division by zero and ensuring numerical stability
/// in Coulomb interaction calculations.
///
/// The value corresponds to approximately 5.29 × 10⁻¹³ Å.
pub const DISTANCE_THRESHOLD_BOHR: f64 = 1e-12;
