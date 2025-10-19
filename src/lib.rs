pub mod error;
pub mod params;
pub mod types;

use crate::params::Parameters;
use std::sync::OnceLock;

static DEFAULT_PARAMETERS: OnceLock<Parameters> = OnceLock::new();

pub fn get_default_parameters() -> &'static Parameters {
    DEFAULT_PARAMETERS.get_or_init(|| {
        const DEFAULT_PARAMS_TOML: &str = include_str!("../resources/qeq.data.toml");
        Parameters::load_from_str(DEFAULT_PARAMS_TOML)
            .expect("Failed to parse embedded default parameters. This is a library bug.")
    })
}
