use clap::{Args, Parser, ValueEnum};
use std::path::PathBuf;

const AUTHORS: &str = "Tony Kan, Ted Yu, William A. Goddard III";
const ABOUT: &str =
    "A command-line tool for calculating dynamic partial atomic charges using the QEq method.";
const COPYRIGHT: &str = "Copyright (c) 2025 California Institute of Technology, Materials and Process Simulation Center (MSC)";
const HELP_TEMPLATE: &str = "\
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}
";

#[derive(Parser)]
#[command(
    author = AUTHORS,
    version,
    about = ABOUT,
    after_help = COPYRIGHT,
    help_template = HELP_TEMPLATE,
)]
#[command(propagate_version = true)]
pub struct Cli {
    /// Input file containing molecular structure in XYZ format.
    ///
    /// Use '-' to read from standard input. The XYZ format should contain the number of atoms
    /// on the first line, a comment on the second line, followed by lines with element symbol
    /// (or atomic number) and x, y, z coordinates.
    #[arg(value_name = "INPUT")]
    pub input: String,

    #[command(flatten)]
    pub output: OutputOptions,

    #[command(flatten)]
    pub calculation: CalculationOptions,

    #[command(flatten)]
    pub solver: SolverOptions,
}

/// Options for controlling the output format and destination.
#[derive(Args)]
#[command(next_help_heading = "Output Options")]
pub struct OutputOptions {
    /// Output file path.
    ///
    /// If not specified, results are written to standard output.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Output format for the results.
    #[arg(short, long, value_enum, default_value_t = OutputFormat::Pretty)]
    pub format: OutputFormat,

    /// Number of decimal places to display for floating-point values.
    #[arg(short, long, default_value_t = 6)]
    pub precision: usize,
}

/// Options for controlling the calculation parameters.
#[derive(Args)]
#[command(next_help_heading = "Calculation Options")]
pub struct CalculationOptions {
    /// Custom parameters file in TOML format.
    ///
    /// If not specified, built-in default parameters are used.
    #[arg(short = 'P', long, value_name = "FILE")]
    pub params: Option<PathBuf>,

    /// Total charge of the molecular system.
    #[arg(short = 'q', long, default_value_t = 0.0)]
    pub total_charge: f64,
}

/// Options for controlling the solver behavior.
#[derive(Args)]
#[command(next_help_heading = "Solver Options")]
pub struct SolverOptions {
    /// Convergence tolerance for charge equilibration.
    ///
    /// The solver iterates until the RMS change in charges falls below this threshold.
    #[arg(long, default_value_t = 1e-6)]
    pub tolerance: f64,

    /// Maximum number of iterations allowed.
    #[arg(long, default_value_t = 100)]
    pub max_iterations: u32,

    /// Screening parameter scale factor.
    ///
    /// Multiplier for the orbital screening strength in Coulomb calculations.
    #[arg(long, default_value_t = 0.5)]
    pub lambda_scale: f64,
}

/// Output format for the calculation results.
#[derive(Clone, ValueEnum)]
pub enum OutputFormat {
    /// Pretty-printed table with atom indices, elements, positions, and charges.
    Pretty,
    /// XYZ format with charges appended to each atom line.
    Xyz,
    /// Comma-separated values with columns: index, element, x, y, z, charge.
    Csv,
    /// JSON object containing atoms array and metadata.
    Json,
}
