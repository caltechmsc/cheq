//! This module provides mathematical utilities and physical constants for the cheq library.
//!
//! It contains fundamental constants for unit conversions and numerical thresholds, as well as
//! functions for computing screened electrostatic interactions between atoms. These components
//! support the charge equilibration algorithm by providing the necessary mathematical infrastructure
//! for potential energy calculations and parameter handling.

/// Physical and mathematical constants used throughout the library.
///
/// This module defines essential constants for unit conversions (Bohr to angstrom, Hartree to eV)
/// and numerical tolerances, ensuring consistent handling of different energy and length scales
/// in computational chemistry calculations.
pub mod constants;

/// Functions for computing screened Coulomb potentials.
///
/// This module implements the orbital screening model used in charge equilibration, providing
/// analytical integrals for electrostatic interactions between atoms with Gaussian orbital
/// approximations. The screening accounts for charge penetration and orbital overlap effects.
pub mod shielding;
