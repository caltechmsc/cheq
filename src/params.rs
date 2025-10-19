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

fn deserialize_element_map<'de, D>(deserializer: D) -> Result<HashMap<u8, ElementData>, D::Error>
where
    D: Deserializer<'de>,
{
    struct ElementMapVisitor;

    impl<'de> Visitor<'de> for ElementMapVisitor {
        type Value = HashMap<u8, ElementData>;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("a map from atomic number or symbol to element data")
        }

        fn visit_map<M>(self, mut map: M) -> Result<Self::Value, M::Error>
        where
            M: MapAccess<'de>,
        {
            let mut elements = HashMap::with_capacity(map.size_hint().unwrap_or(0));
            while let Some((key, value)) = map.next_entry::<String, ElementData>()? {
                let atomic_number = key.parse::<u8>().or_else(|_| {
                    element_symbol_to_atomic_number(&key)
                        .ok_or_else(|| de::Error::custom(format!("invalid element key: '{}'", key)))
                })?;
                elements.insert(atomic_number, value);
            }
            Ok(elements)
        }
    }

    deserializer.deserialize_map(ElementMapVisitor)
}
