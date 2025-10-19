use std::path::PathBuf;
use thiserror::Error;

/// The primary error type for all fallible operations in the `cheq` library.
///
/// This enum is designed to be comprehensive, providing clear and actionable
/// information for each potential failure mode, from I/O issues to numerical
/// convergence problems. It implements `std::error::Error`, allowing it to be
/// composed with other error types in application code.
#[derive(Error, Debug)]
pub enum CheqError {
    /// Indicates that the parameters for a specific element, identified by its
    /// atomic number, could not be found in the provided `Parameters` set.
    ///
    /// This is a common error when attempting to calculate charges for a molecule
    /// containing an element not defined in the parameter file.
    #[error("Atomic parameters not found for element with atomic number: {0}")]
    ParameterNotFound(u8),

    /// Occurs when the Self-Consistent Field (SCF) iterative process fails to
    /// converge to the desired tolerance within the maximum number of allowed iterations.
    ///
    /// This can happen with unusual molecular geometries or problematic parameter sets.
    /// The final charge difference (`delta`) is provided for diagnostic purposes.
    #[error(
        "SCF failed to converge after {max_iterations} iterations. Final charge difference: {delta:.2e}"
    )]
    NotConverged {
        /// The maximum number of iterations that were performed.
        max_iterations: u32,
        /// The largest absolute change in any atomic charge during the final iteration.
        delta: f64,
    },

    /// A failure within the underlying linear algebra solver, for example,
    /// if the matrix system is singular or ill-conditioned.
    #[error("Failed to solve the linear matrix system: {0}")]
    LinalgError(String),

    /// An I/O error that occurred while attempting to read a parameter file.
    ///
    /// The path to the file and the underlying I/O error are provided for context.
    #[error("I/O error at path '{path}': {source}")]
    IoError {
        /// The path of the file that caused the I/O error.
        path: PathBuf,
        /// The underlying `std::io::Error`.
        #[source]
        source: std::io::Error,
    },

    /// An error that occurred while parsing a parameter file, typically indicating
    /// invalid TOML or a structural mismatch with the expected `Parameters` format.
    #[error("Failed to deserialize TOML parameters: {0}")]
    DeserializationError(#[from] toml::de::Error),

    /// A validation error indicating that the input slice of atoms was empty.
    /// At least one atom is required to perform a calculation.
    #[error("Input validation failed: at least one atom is required for a calculation")]
    NoAtoms,
}
