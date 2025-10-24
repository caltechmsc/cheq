use super::cli::OutputFormat;
use super::error::CliError;
use cheq::Atom;
use prettytable::*;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

pub fn read_atoms(input_spec: &str) -> Result<(Vec<Atom>, String), CliError> {
    let reader: Box<dyn BufRead> = if input_spec == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else {
        let file = std::fs::File::open(input_spec).map_err(|e| CliError::Io {
            path: PathBuf::from(input_spec),
            source: e,
        })?;
        Box::new(BufReader::new(file))
    };

    let mut lines = reader.lines();

    let num_atoms_line = lines.next().ok_or_else(|| CliError::XyzParse {
        source_name: input_spec.to_string(),
        details: "Missing number of atoms line".to_string(),
    })??;
    let num_atoms: usize = num_atoms_line
        .trim()
        .parse()
        .map_err(|_| CliError::XyzParse {
            source_name: input_spec.to_string(),
            details: format!("Invalid number of atoms: {}", num_atoms_line),
        })?;

    let comment = lines.next().ok_or_else(|| CliError::XyzParse {
        source_name: input_spec.to_string(),
        details: "Missing comment line".to_string(),
    })??;

    let mut atoms = Vec::with_capacity(num_atoms);
    for (i, line) in lines.enumerate() {
        if i >= num_atoms {
            break;
        }
        let line = line.map_err(|e| CliError::XyzParse {
            source_name: input_spec.to_string(),
            details: format!("Error reading line {}: {}", i + 3, e),
        })?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(CliError::XyzParse {
                source_name: input_spec.to_string(),
                details: format!(
                    "Line {}: expected at least 4 fields, got {}",
                    i + 3,
                    parts.len()
                ),
            });
        }
        let atomic_number = parse_element(parts[0]).ok_or_else(|| CliError::XyzParse {
            source_name: input_spec.to_string(),
            details: format!("Unknown element: {}", parts[0]),
        })?;
        let x: f64 = parts[1].parse().map_err(|_| CliError::XyzParse {
            source_name: input_spec.to_string(),
            details: format!("Invalid x coordinate: {}", parts[1]),
        })?;
        let y: f64 = parts[2].parse().map_err(|_| CliError::XyzParse {
            source_name: input_spec.to_string(),
            details: format!("Invalid y coordinate: {}", parts[2]),
        })?;
        let z: f64 = parts[3].parse().map_err(|_| CliError::XyzParse {
            source_name: input_spec.to_string(),
            details: format!("Invalid z coordinate: {}", parts[3]),
        })?;
        atoms.push(Atom {
            atomic_number,
            position: [x, y, z],
        });
    }

    if atoms.len() != num_atoms {
        return Err(CliError::XyzParse {
            source_name: input_spec.to_string(),
            details: format!("Expected {} atoms, got {}", num_atoms, atoms.len()),
        });
    }

    Ok((atoms, comment))
}
