use super::error::CheqError;
use serde::Deserialize;
use serde::de::{self, Deserializer, MapAccess, Visitor};
use std::collections::HashMap;
use std::fmt;
use std::path::Path;

#[derive(Deserialize, Debug, Clone, PartialEq)]
pub struct ElementData {
    #[serde(rename = "chi")]
    pub electronegativity: f64,
    #[serde(rename = "j")]
    pub hardness: f64,
    pub radius: f64,
    #[serde(rename = "n")]
    pub principal_quantum_number: u8,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Parameters {
    #[serde(deserialize_with = "deserialize_element_map")]
    pub elements: HashMap<u8, ElementData>,
}

impl Parameters {
    pub fn load_from_file(path: &Path) -> Result<Self, CheqError> {
        let content = std::fs::read_to_string(path).map_err(|io_error| CheqError::IoError {
            path: path.to_path_buf(),
            source: io_error,
        })?;

        Self::load_from_str(&content)
    }

    pub fn load_from_str(toml_str: &str) -> Result<Self, CheqError> {
        toml::from_str(toml_str).map_err(CheqError::from)
    }

    pub fn new() -> Self {
        Parameters {
            elements: HashMap::new(),
        }
    }
}

impl Default for Parameters {
    fn default() -> Self {
        Self::new()
    }
}
