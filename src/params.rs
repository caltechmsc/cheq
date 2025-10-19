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
