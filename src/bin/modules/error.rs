use std::path::PathBuf;

#[derive(thiserror::Error, Debug)]
pub enum CliError {
    /// Errors originating from the core cheq library calculations.
    #[error("Calculation error: {0}")]
    Calculation(#[from] cheq::CheqError),

    /// I/O errors associated with a specific file path.
    #[error("I/O error for '{}': {source}", .path.display())]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// General I/O errors not tied to a specific file.
    #[error("I/O error: {0}")]
    GenericIo(#[from] std::io::Error),

    /// Errors parsing XYZ format input.
    #[error("Failed to parse XYZ from {source_name}: {details}")]
    XyzParse {
        source_name: String,
        details: String,
    },

    /// Errors parsing custom parameter TOML files.
    #[error("Failed to parse parameters TOML: {0}")]
    ParamsParse(#[from] toml::de::Error),
}
