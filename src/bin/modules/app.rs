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

#[cfg(test)]
mod tests {
    use super::super::cli::{CalculationOptions, OutputFormat, OutputOptions, SolverOptions};
    use super::super::error::CliError;
    use super::*;
    use std::fs;
    use std::path::Path;
    use tempfile::TempDir;

    fn write_test_xyz(path: &Path) {
        fs::write(path, "2\nHydrogen\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n").unwrap();
    }

    fn default_solver_options() -> SolverOptions {
        SolverOptions {
            tolerance: 1e-6,
            max_iterations: 50,
            lambda_scale: 0.5,
        }
    }

    #[test]
    fn test_run_generates_pretty_output() {
        let temp_dir = TempDir::new().unwrap();
        let input_path = temp_dir.path().join("input.xyz");
        write_test_xyz(&input_path);
        let output_path = temp_dir.path().join("output.txt");

        let args = Cli {
            input: input_path.to_string_lossy().into_owned(),
            output: OutputOptions {
                output: Some(output_path.clone()),
                format: OutputFormat::Pretty,
                precision: 4,
            },
            calculation: CalculationOptions {
                params: None,
                total_charge: 0.0,
            },
            solver: default_solver_options(),
        };

        run(args).unwrap();

        let output_contents = fs::read_to_string(&output_path).unwrap();
        assert!(output_contents.contains("Cheq Charge Equilibration Results"));
        assert!(output_contents.contains("Total Atoms:"));
    }

    #[test]
    fn test_run_invalid_params_file_returns_error() {
        let temp_dir = TempDir::new().unwrap();
        let input_path = temp_dir.path().join("input.xyz");
        write_test_xyz(&input_path);
        let bad_params_path = temp_dir.path().join("bad.toml");
        fs::write(&bad_params_path, "not valid toml").unwrap();

        let args = Cli {
            input: input_path.to_string_lossy().into_owned(),
            output: OutputOptions {
                output: None,
                format: OutputFormat::Pretty,
                precision: 4,
            },
            calculation: CalculationOptions {
                params: Some(bad_params_path),
                total_charge: 0.0,
            },
            solver: default_solver_options(),
        };

        let error = run(args).unwrap_err();
        assert!(matches!(error, CliError::ParamsParse(_)));
    }
}
