#[path = "modules/app.rs"]
mod app;
#[path = "modules/cli.rs"]
mod cli;
#[path = "modules/error.rs"]
mod error;
#[path = "modules/io.rs"]
mod io;

use clap::Parser;
use std::error::Error;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args = cli::Cli::parse();

    if let Err(e) = app::run(args) {
        eprintln!("Error: {}", e);

        let mut source = e.source();
        while let Some(s) = source {
            eprintln!("Caused by: {}", s);
            source = s.source();
        }

        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    }
}
