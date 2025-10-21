pub mod error;
pub mod math;
pub mod params;
pub mod solver;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_default_parameters() {
        let params1 = get_default_parameters();
        assert!(
            params1.elements.get(&6).is_some(),
            "Carbon (6) should be present"
        );
        assert!(
            params1.elements.get(&8).is_some(),
            "Oxygen (8) should be present"
        );

        let params2 = get_default_parameters();
        assert_eq!(
            params1 as *const _, params2 as *const _,
            "Subsequent calls should return a cached reference"
        );
    }
}
