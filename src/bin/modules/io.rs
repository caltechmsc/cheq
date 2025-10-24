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

fn parse_element(s: &str) -> Option<u8> {
    if let Ok(num) = s.parse::<u8>() {
        return Some(num);
    }

    match s.to_uppercase().as_str() {
        "H" => Some(1),
        "HE" => Some(2),
        "LI" => Some(3),
        "BE" => Some(4),
        "B" => Some(5),
        "C" => Some(6),
        "N" => Some(7),
        "O" => Some(8),
        "F" => Some(9),
        "NE" => Some(10),
        "NA" => Some(11),
        "MG" => Some(12),
        "AL" => Some(13),
        "SI" => Some(14),
        "P" => Some(15),
        "S" => Some(16),
        "CL" => Some(17),
        "AR" => Some(18),
        "K" => Some(19),
        "CA" => Some(20),
        "SC" => Some(21),
        "TI" => Some(22),
        "V" => Some(23),
        "CR" => Some(24),
        "MN" => Some(25),
        "FE" => Some(26),
        "CO" => Some(27),
        "NI" => Some(28),
        "CU" => Some(29),
        "ZN" => Some(30),
        "GA" => Some(31),
        "GE" => Some(32),
        "AS" => Some(33),
        "SE" => Some(34),
        "BR" => Some(35),
        "KR" => Some(36),
        "RB" => Some(37),
        "SR" => Some(38),
        "Y" => Some(39),
        "ZR" => Some(40),
        "NB" => Some(41),
        "MO" => Some(42),
        "TC" => Some(43),
        "RU" => Some(44),
        "RH" => Some(45),
        "PD" => Some(46),
        "AG" => Some(47),
        "CD" => Some(48),
        "IN" => Some(49),
        "SN" => Some(50),
        "SB" => Some(51),
        "TE" => Some(52),
        "I" => Some(53),
        "XE" => Some(54),
        "CS" => Some(55),
        "BA" => Some(56),
        "LA" => Some(57),
        "CE" => Some(58),
        "PR" => Some(59),
        "ND" => Some(60),
        "PM" => Some(61),
        "SM" => Some(62),
        "EU" => Some(63),
        "GD" => Some(64),
        "TB" => Some(65),
        "DY" => Some(66),
        "HO" => Some(67),
        "ER" => Some(68),
        "TM" => Some(69),
        "YB" => Some(70),
        "LU" => Some(71),
        "HF" => Some(72),
        "TA" => Some(73),
        "W" => Some(74),
        "RE" => Some(75),
        "OS" => Some(76),
        "IR" => Some(77),
        "PT" => Some(78),
        "AU" => Some(79),
        "HG" => Some(80),
        "TL" => Some(81),
        "PB" => Some(82),
        "BI" => Some(83),
        "PO" => Some(84),
        "AT" => Some(85),
        "RN" => Some(86),
        "FR" => Some(87),
        "RA" => Some(88),
        "AC" => Some(89),
        "TH" => Some(90),
        "PA" => Some(91),
        "U" => Some(92),
        "NP" => Some(93),
        "PU" => Some(94),
        "AM" => Some(95),
        "CM" => Some(96),
        "BK" => Some(97),
        "CF" => Some(98),
        "ES" => Some(99),
        "FM" => Some(100),
        "MD" => Some(101),
        "NO" => Some(102),
        "LR" => Some(103),
        "RF" => Some(104),
        "DB" => Some(105),
        "SG" => Some(106),
        "BH" => Some(107),
        "HS" => Some(108),
        "MT" => Some(109),
        "DS" => Some(110),
        "RG" => Some(111),
        "CN" => Some(112),
        "FL" => Some(114),
        "LV" => Some(116),
        "TS" => Some(117),
        "OG" => Some(118),
        _ => None,
    }
}

pub fn get_writer(output_path: &Option<PathBuf>) -> Result<Box<dyn Write>, CliError> {
    match output_path {
        Some(path) => {
            let file = std::fs::File::create(path).map_err(|e| CliError::Io {
                path: path.clone(),
                source: e,
            })?;
            Ok(Box::new(BufWriter::new(file)))
        }
        None => Ok(Box::new(io::stdout())),
    }
}

pub fn write_results(
    mut writer: Box<dyn Write>,
    atoms: &[Atom],
    result: &cheq::CalculationResult,
    comment: &str,
    format: &OutputFormat,
    precision: usize,
    source_name: &str,
) -> Result<(), CliError> {
    match format {
        OutputFormat::Pretty => {
            write_pretty_table(&mut writer, atoms, result, precision, source_name)
        }
        OutputFormat::Xyz => write_xyz_charged(&mut writer, atoms, result, comment, precision),
        OutputFormat::Csv => write_csv(&mut writer, atoms, result, precision),
        OutputFormat::Json => write_json(&mut writer, atoms, result, precision),
    }
}

fn write_pretty_table(
    writer: &mut dyn Write,
    atoms: &[Atom],
    result: &cheq::CalculationResult,
    precision: usize,
    source_name: &str,
) -> Result<(), CliError> {
    let box_format = format::FormatBuilder::new()
        .column_separator('│')
        .borders('│')
        .separators(
            &[format::LinePosition::Top],
            format::LineSeparator::new('─', '┬', '╭', '╮'),
        )
        .separators(
            &[format::LinePosition::Title],
            format::LineSeparator::new('═', '╪', '╞', '╡'),
        )
        .separators(
            &[format::LinePosition::Intern],
            format::LineSeparator::new('─', '┼', '├', '┤'),
        )
        .separators(
            &[format::LinePosition::Bottom],
            format::LineSeparator::new('─', '┴', '╰', '╯'),
        )
        .padding(1, 1)
        .build();

    let no_intern_format = format::FormatBuilder::new()
        .column_separator('│')
        .borders('│')
        .separators(
            &[format::LinePosition::Top],
            format::LineSeparator::new('─', '┬', '╭', '╮'),
        )
        .separators(
            &[format::LinePosition::Bottom],
            format::LineSeparator::new('─', '┴', '╰', '╯'),
        )
        .padding(1, 1)
        .build();

    let total_charge = result.charges.iter().sum::<f64>();

    let mut title_table = Table::new();
    title_table.set_format(box_format.clone());
    title_table.add_row(row![bc->"Cheq Charge Equilibration Results"]);
    title_table.print(writer)?;
    writeln!(writer)?;

    let mut summary_table = Table::new();
    summary_table.set_format(no_intern_format);
    summary_table.add_row(row![b->"Source File:", source_name]);
    summary_table.add_row(row![b->"Total Atoms:", atoms.len()]);
    summary_table
        .add_row(row![b->"Total Charge:", format!("{:.prec$} e", total_charge, prec = precision)]);
    summary_table.add_row(row![b->"Converged Iterations:", result.iterations]);
    summary_table.add_row(row![b->"Equilibrated Potential:", format!("{:.prec$} eV", result.equilibrated_potential, prec = precision)]);
    summary_table.print(writer)?;
    writeln!(writer)?;

    let mut data_table = Table::new();
    data_table.set_format(box_format);
    data_table.set_titles(
        row![bc->"Index", bc->"Element", bc->"X (Å)", bc->"Y (Å)", bc->"Z (Å)", bc->"Charge (e)"],
    );

    for (i, (atom, &charge)) in atoms.iter().zip(result.charges.iter()).enumerate() {
        let symbol = atomic_number_to_symbol(atom.atomic_number).unwrap_or("??");
        data_table.add_row(row![
            r->i,
            l->symbol,
            r->format!("{:.prec$}", atom.position[0], prec = precision),
            r->format!("{:.prec$}", atom.position[1], prec = precision),
            r->format!("{:.prec$}", atom.position[2], prec = precision),
            r->format!("{:.prec$}", charge, prec = precision)
        ]);
    }

    data_table.print(writer)?;

    Ok(())
}

fn write_xyz_charged(
    writer: &mut dyn Write,
    atoms: &[Atom],
    result: &cheq::CalculationResult,
    comment: &str,
    precision: usize,
) -> Result<(), CliError> {
    writeln!(writer, "{}", atoms.len())?;
    writeln!(
        writer,
        "{} | QEq charges | iterations: {} | potential: {:.*}",
        comment.trim(),
        result.iterations,
        precision,
        result.equilibrated_potential
    )?;
    for (atom, &charge) in atoms.iter().zip(result.charges.iter()) {
        let symbol = atomic_number_to_symbol(atom.atomic_number).unwrap_or("??");
        writeln!(
            writer,
            "{} {:.*} {:.*} {:.*} {:.*}",
            symbol,
            precision,
            atom.position[0],
            precision,
            atom.position[1],
            precision,
            atom.position[2],
            precision,
            charge
        )?;
    }
    Ok(())
}
