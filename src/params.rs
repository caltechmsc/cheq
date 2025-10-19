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

fn element_symbol_to_atomic_number(symbol: &str) -> Option<u8> {
    match symbol {
        "H" => Some(1),
        "He" => Some(2),
        "Li" => Some(3),
        "Be" => Some(4),
        "B" => Some(5),
        "C" => Some(6),
        "N" => Some(7),
        "O" => Some(8),
        "F" => Some(9),
        "Ne" => Some(10),
        "Na" => Some(11),
        "Mg" => Some(12),
        "Al" => Some(13),
        "Si" => Some(14),
        "P" => Some(15),
        "S" => Some(16),
        "Cl" => Some(17),
        "Ar" => Some(18),
        "K" => Some(19),
        "Ca" => Some(20),
        "Sc" => Some(21),
        "Ti" => Some(22),
        "V" => Some(23),
        "Cr" => Some(24),
        "Mn" => Some(25),
        "Fe" => Some(26),
        "Co" => Some(27),
        "Ni" => Some(28),
        "Cu" => Some(29),
        "Zn" => Some(30),
        "Ga" => Some(31),
        "Ge" => Some(32),
        "As" => Some(33),
        "Se" => Some(34),
        "Br" => Some(35),
        "Kr" => Some(36),
        "Rb" => Some(37),
        "Sr" => Some(38),
        "Y" => Some(39),
        "Zr" => Some(40),
        "Nb" => Some(41),
        "Mo" => Some(42),
        "Tc" => Some(43),
        "Ru" => Some(44),
        "Rh" => Some(45),
        "Pd" => Some(46),
        "Ag" => Some(47),
        "Cd" => Some(48),
        "In" => Some(49),
        "Sn" => Some(50),
        "Sb" => Some(51),
        "Te" => Some(52),
        "I" => Some(53),
        "Xe" => Some(54),
        "Cs" => Some(55),
        "Ba" => Some(56),
        "La" => Some(57),
        "Ce" => Some(58),
        "Pr" => Some(59),
        "Nd" => Some(60),
        "Pm" => Some(61),
        "Sm" => Some(62),
        "Eu" => Some(63),
        "Gd" => Some(64),
        "Tb" => Some(65),
        "Dy" => Some(66),
        "Ho" => Some(67),
        "Er" => Some(68),
        "Tm" => Some(69),
        "Yb" => Some(70),
        "Lu" => Some(71),
        "Hf" => Some(72),
        "Ta" => Some(73),
        "W" => Some(74),
        "Re" => Some(75),
        "Os" => Some(76),
        "Ir" => Some(77),
        "Pt" => Some(78),
        "Au" => Some(79),
        "Hg" => Some(80),
        "Tl" => Some(81),
        "Pb" => Some(82),
        "Bi" => Some(83),
        "Po" => Some(84),
        "At" => Some(85),
        "Rn" => Some(86),
        "Fr" => Some(87),
        "Ra" => Some(88),
        "Ac" => Some(89),
        "Th" => Some(90),
        "Pa" => Some(91),
        "U" => Some(92),
        "Np" => Some(93),
        "Pu" => Some(94),
        "Am" => Some(95),
        "Cm" => Some(96),
        "Bk" => Some(97),
        "Cf" => Some(98),
        "Es" => Some(99),
        "Fm" => Some(100),
        "Md" => Some(101),
        "No" => Some(102),
        "Lr" => Some(103),
        "Rf" => Some(104),
        "Db" => Some(105),
        "Sg" => Some(106),
        "Bh" => Some(107),
        "Hs" => Some(108),
        "Mt" => Some(109),
        "Ds" => Some(110),
        "Rg" => Some(111),
        "Cn" => Some(112),
        "Nh" => Some(113),
        "Fl" => Some(114),
        "Mc" => Some(115),
        "Lv" => Some(116),
        "Ts" => Some(117),
        "Og" => Some(118),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_toml_string() -> String {
        r#"
        [elements]
        "1" = { chi = 2.20, j = 13.60, radius = 0.37, n = 1 }
        "Fe" = { chi = 1.83, j = 7.90, radius = 1.26, n = 4 }
        "O" = { chi = 3.44, j = 13.62, radius = 0.60, n = 2 }
        "#
        .to_string()
    }

    fn get_expected_parameters() -> Parameters {
        let mut elements = HashMap::new();
        elements.insert(
            1,
            ElementData {
                electronegativity: 2.20,
                hardness: 13.60,
                radius: 0.37,
                principal_quantum_number: 1,
            },
        );
        elements.insert(
            26,
            ElementData {
                electronegativity: 1.83,
                hardness: 7.90,
                radius: 1.26,
                principal_quantum_number: 4,
            },
        );
        elements.insert(
            8,
            ElementData {
                electronegativity: 3.44,
                hardness: 13.62,
                radius: 0.60,
                principal_quantum_number: 2,
            },
        );
        Parameters { elements }
    }

    #[test]
    fn test_load_from_str_valid() {
        let toml_str = create_test_toml_string();
        let params = Parameters::load_from_str(&toml_str).unwrap();
        let expected_params = get_expected_parameters();
        assert_eq!(params, expected_params);
    }

    #[test]
    fn test_load_from_str_mixed_keys() {
        let toml_str = r#"
        [elements]
        "1" = { chi = 2.20, j = 13.60, radius = 0.37, n = 1 } # Hydrogen by atomic number
        "Fe" = { chi = 1.83, j = 7.90, radius = 1.26, n = 4 }  # Iron by symbol
        "#;
        let params = Parameters::load_from_str(toml_str).unwrap();
        let mut elements = HashMap::new();
        elements.insert(
            1,
            ElementData {
                electronegativity: 2.20,
                hardness: 13.60,
                radius: 0.37,
                principal_quantum_number: 1,
            },
        );
        elements.insert(
            26,
            ElementData {
                electronegativity: 1.83,
                hardness: 7.90,
                radius: 1.26,
                principal_quantum_number: 4,
            },
        );
        assert_eq!(params.elements, elements);
    }

    #[test]
    fn test_load_from_str_invalid_toml() {
        let toml_str = "this is not valid toml";
        let result = Parameters::load_from_str(toml_str);
        assert!(matches!(result, Err(CheqError::DeserializationError(_))));
    }

    #[test]
    fn test_load_from_str_invalid_element_key() {
        let toml_str = r#"
        [elements]
        "InvalidKey" = { chi = 1.0, j = 1.0, radius = 1.0, n = 1 }
        "#;
        let result = Parameters::load_from_str(toml_str);
        assert!(result.is_err());
        let error_string = result.unwrap_err().to_string();
        assert!(error_string.contains("invalid element key: 'InvalidKey'"));
    }

    #[test]
    fn test_load_from_str_missing_field() {
        let toml_str = r#"
        [elements]
        "1" = { chi = 2.20, j = 13.60, radius = 0.37 } # Missing 'n'
        "#;
        let result = Parameters::load_from_str(toml_str);
        assert!(matches!(result, Err(CheqError::DeserializationError(_))));
    }

    #[test]
    fn test_load_from_file_valid() {
        let toml_str = create_test_toml_string();
        let mut temp_file = NamedTempFile::new().unwrap();
        write!(temp_file, "{}", toml_str).unwrap();

        let params = Parameters::load_from_file(temp_file.path()).unwrap();
        let expected_params = get_expected_parameters();
        assert_eq!(params, expected_params);
    }

    #[test]
    fn test_load_from_file_not_found() {
        let path = Path::new("non_existent_file.toml");
        let result = Parameters::load_from_file(path);
        assert!(matches!(result, Err(CheqError::IoError { .. })));
    }

    #[test]
    fn test_new_and_default() {
        let params_new = Parameters::new();
        let params_default = Parameters::default();
        assert_eq!(params_new.elements.len(), 0);
        assert_eq!(params_default.elements.len(), 0);
        assert_eq!(params_new, params_default);
    }

    #[test]
    fn test_element_symbol_to_atomic_number() {
        assert_eq!(element_symbol_to_atomic_number("H"), Some(1));
        assert_eq!(element_symbol_to_atomic_number("O"), Some(8));
        assert_eq!(element_symbol_to_atomic_number("Fe"), Some(26));
        assert_eq!(element_symbol_to_atomic_number("Og"), Some(118));
        assert_eq!(element_symbol_to_atomic_number("Xx"), None);
        assert_eq!(element_symbol_to_atomic_number("h"), None);
    }
}
