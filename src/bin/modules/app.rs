use super::cli::Cli;
use super::error::CliError;
use super::io;
use cheq::{QEqSolver, SolverOptions, get_default_parameters};
use indicatif::{ProgressBar, ProgressStyle};
use std::fs;

pub fn run(args: Cli) -> Result<(), CliError> {
    let params = if let Some(params_path) = &args.calculation.params {
        let content = fs::read_to_string(params_path).map_err(|e| CliError::Io {
            path: params_path.clone(),
            source: e,
        })?;
        toml::from_str(&content)?
    } else {
        get_default_parameters().clone()
    };

    let solver_options = SolverOptions {
        tolerance: args.solver.tolerance,
        max_iterations: args.solver.max_iterations,
        lambda_scale: args.solver.lambda_scale,
    };
    let solver = QEqSolver::new(&params).with_options(solver_options);

    let (atoms, comment) = io::read_atoms(&args.input)?;

    let source_name = if args.input == "-" {
        "stdin".to_string()
    } else {
        args.input.clone()
    };

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );
    pb.set_message("Calculating partial charges...");
    pb.enable_steady_tick(std::time::Duration::from_millis(100));

    let result = solver.solve(&atoms, args.calculation.total_charge)?;

    pb.finish_and_clear();

    let writer = io::get_writer(&args.output.output)?;
    io::write_results(
        writer,
        &atoms,
        &result,
        &comment,
        &args.output.format,
        args.output.precision,
        &source_name,
    )?;

    Ok(())
}
